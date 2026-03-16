# ONT long-read WGS pipeline for FSHD diagnosis

**This pipeline classifies samples as FSHD1, FSHD2, or unaffected using standard whole-genome ONT sequencing on R10.4.1 chemistry, without CRISPR-Cas9 enrichment or adaptive sampling.** The input is a QC-filtered hg38-aligned BAM with modification tags (MM/ML) intact. The core constraint—WGS-only—means D4Z4-region coverage will be **~5–15× effective depth** (after multi-mapping filters) at 30× mean genome coverage, substantially below the 50–500× achieved by Cas9-targeted approaches. The pipeline compensates through flanking-anchor read selection, BLAST-based chromosome assignment, and a 100-sample non-FSHD reference panel to establish baseline distributions. All downstream coordinates reference T2T-CHM13v2.0.

---

## Step 1: targeted read extraction from hg38 BAM and realignment to T2T-CHM13v2.0

Rather than realigning the entire genome (~30× WGS = ~90 Gbp), extract only reads from FSHD-relevant regions plus unmapped reads, then realign this subset to T2T-CHM13v2.0. This reduces alignment input by **~95%**.

### 1a. Define hg38 extraction regions

```bash
cat > fshd_regions_hg38.bed << 'EOF'
chr4	189900000	190300000	4q35_D4Z4_flanking
chr10	133600000	133900000	10q26_D4Z4_flanking
chr18	2655000	2810000	SMCHD1
chr20	31350000	31400000	DNMT3B
chr1	73900000	73970000	LRIF1
chr3	195460000	195680000	3q29_D4Z4_like
chr10	41900000	42100000	10p_D4Z4_like
chr18	78010000	78200000	18q_D4Z4_like
chr21	1	500000	21p_acrocentric
chr22	1	500000	22p_acrocentric
EOF
```

The hg38 4q35 region includes a **~82 kb assembly gap** (chr4:190,030,000–190,112,000 in GRCh38) where D4Z4-containing reads often end up unmapped or poorly mapped. The 400 kb extraction window around each locus ensures we capture split/supplementary alignments whose primary mapping landed in flanking sequence but whose clipped portion spans into D4Z4. The D4Z4-like regions on chr3, chr10p, chr18, and acrocentric chromosomes are included because reads from 4q/10q D4Z4 may have been misassigned to these paralogs under hg38.

### 1b. Extract target reads + unmapped reads

```bash
# Index if not already
samtools index input_hg38.bam

# Extract reads from target regions (includes mates, supplementary, secondary)
samtools view -b -h -M -L fshd_regions_hg38.bed \
  --fetch-pairs input_hg38.bam \
  -o target_reads.bam

# Extract unmapped reads (critical: hg38 gap region causes many D4Z4 reads to be unmapped)
samtools view -b -f 4 input_hg38.bam -o unmapped_reads.bam

# Merge
samtools merge -o fshd_subset.bam target_reads.bam unmapped_reads.bam
samtools index fshd_subset.bam
```

The `-M` flag uses the BED as a multi-region iterator (faster than repeated region queries). `--fetch-pairs` ensures both reads of a pair are extracted even if only one overlaps a target region (relevant if supplementary alignments span distant locations). The unmapped read pool is essential: in hg38, reads spanning the 4q D4Z4 gap frequently fail to map; these reads often align successfully to the gapless T2T reference.

**Expected yield:** ~0.5–2 GB from a 30× WGS BAM (~60–80 GB), depending on how many reads fall in the gap regions.

### 1c. Revert to unaligned BAM preserving modification tags

```bash
samtools fastq -T MM,ML fshd_subset.bam \
  | samtools import -T MM,ML - \
  -o fshd_unaligned.bam
```

The `-T MM,ML` flag in `samtools fastq` propagates modification tags to the FASTQ header as auxiliary tag comments, and `samtools import` restores them into the unaligned BAM. This round-trip preserves per-read 5mC probabilities without re-basecalling from POD5.

**Alternative (simpler, if samtools ≥1.19):**

```bash
samtools reheader --no-PG -c 'grep -v "^@SQ"' fshd_subset.bam \
  | samtools reset -o fshd_unaligned.bam --keep-tag MM,ML
```

`samtools reset` (v1.19+) strips alignment information (CIGAR, MAPQ, coordinates) while preserving specified tags, producing a clean unaligned BAM without the FASTQ round-trip.

### 1d. Align to T2T-CHM13v2.0

**Primary aligner:** Winnowmap2 v2.03 (preferred for segmental duplications) or minimap2 ≥v2.22

#### Option A — Winnowmap2 (recommended)

```bash
meryl count k=15 output merylDB chm13v2.0.fa
meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt

winnowmap -W repetitive_k15.txt -ax map-ont -y --eqx \
  chm13v2.0.fa fshd_unaligned.bam \
  | samtools sort -@ 16 -o aligned.bam
samtools index aligned.bam
```

Winnowmap2's weighted minimizer scheme down-weights high-frequency k-mers during seeding, preventing the mapping failures minimap2 exhibits in segmental duplications. Its Minimal Confidently Alignable Substrings (MCAS) decomposition reduces allelic bias between 4q and 10q paralogs. In benchmarks on T2T chr8 segmental duplications, Winnowmap2 achieved an FNR of **1.89%** versus minimap2's **39.62%** for SVs in repeats ≥95% identity.

The `-y` flag propagates methylation MM/ML tags from the input BAM to the aligned output. The `--eqx` flag emits `=`/`X` CIGAR operators for exact match/mismatch discrimination, useful for downstream PSV analysis.

#### Option B — minimap2 ≥v2.22

```bash
minimap2 -ax map-ont -y --eqx --secondary=yes -N 50 -p 0.5 \
  chm13v2.0.fa fshd_unaligned.bam \
  | samtools sort -@ 16 -o aligned.bam
samtools index aligned.bam
```

The `--secondary=yes -N 50 -p 0.5` flags retain secondary alignments critical for downstream DuploMap processing. From v2.22, minimap2 adds supplementary minimizers in sparse repetitive regions, partially closing the gap with Winnowmap2.

### Why T2T-CHM13v2.0 over GRCh38

GRCh38 contains an **~82 kb assembly gap** downstream of the 4q D4Z4 array and incomplete subtelomeric assembly at both 4q35 and 10q26. T2T-CHM13v2.0 provides gapless, complete D4Z4 arrays: **~33 D4Z4 repeat units on chr4** (4qA haplotype) and **~32–33 units on chr10**, plus full flanking sequences. T2T-CHM13 also reveals D4Z4-like clusters on ≥10 additional chromosomes (3q, 18q, 21q, 22q, acrocentric short arms), which are absent from GRCh38 and can act as mapping sinks if not accounted for.

**Key T2T-CHM13v2.0 coordinates:**

| Region | Coordinates | Notes |
|--------|------------|-------|
| Chr4 D4Z4 array + flanking | chr4:193,282,466–193,558,261 | ~276 kb; 33 RUs in reference |
| Chr4 distal D4Z4/DUX4 target | chr4:193,540,172–193,543,634 | Methylation analysis target (DUCKS4) |
| Chr4 p13E-11/D4F104S1 region | chr4:192,667,301–192,902,247 | Upstream unique anchoring region |
| Chr10 D4Z4 array | chr10:134,615,120–134,741,148 | ~126 kb |

---

## Step 2: 4q versus 10q disambiguation

The D4Z4 arrays on 4q35 and 10q26 share **>98% sequence identity**, but three independent signals enable chromosome assignment: flanking unique sequence anchors, paralog-specific restriction-site motifs within D4Z4 units, and distal haplotype sequences.

### 2a. Extract candidate D4Z4 reads

```bash
# Extract reads mapping to either D4Z4 locus ± 50 kb flanking
samtools view -b aligned.bam \
  "chr4:193,230,000-193,610,000" \
  "chr10:134,560,000-134,800,000" \
  > d4z4_candidates.bam
```

### 2b. BLAST-based chromosome assignment (DUCKS4/Huang method)

This approach, validated by Huang et al. (2024) with **100% concordance to optical genome mapping** across 29 samples, uses NCBI BLAST to align chromosome-specific feature sequences against each candidate read.

**Feature sequences to BLAST against reads:**

- **D4F104S1** (p13E-11 probe region): 4q-specific vs 10q-specific variants
- **BlnI site motif** (CCTAGG): present in 10q D4Z4 units, absent from canonical 4q units
- **XapI site motif** (RAATTY/AATTCC): present in 4q D4Z4 units, absent from canonical 10q units
- **pLAM sequence** (260 bp): distal to last D4Z4 repeat; 4qA-specific ATTAAA polyadenylation signal vs 10qA ATCAAA

```bash
# Convert reads to FASTA
samtools fasta d4z4_candidates.bam > d4z4_reads.fa

# BLAST against feature database (4q and 10q reference features)
blastn -query d4z4_reads.fa \
  -db d4z4_features.fa \
  -outfmt "6 qseqid sseqid bitscore pident length" \
  -evalue 1e-5 \
  > blast_hits.tsv
```

**Classification logic per read:**

1. Sum bitscores for all 4q-specific feature hits → `score_4q`
2. Sum bitscores for all 10q-specific feature hits → `score_10q`
3. If `score_4q > score_10q × 1.2` → assign to **chr4**
4. If `score_10q > score_4q × 1.2` → assign to **chr10**
5. Otherwise → **ambiguous** (exclude from downstream analysis)

Count BlnI-sensitive (B+) and XapI-sensitive (X+) motifs within each read's D4Z4 repeats as confirmation: **canonical 4q reads are B−X+; canonical 10q reads are B+X−**. Mixed B−X− patterns indicate hybrid/translocated arrays requiring manual review.

### 2c. Optional enhancement: DuploMap post-processing

For reads that remain ambiguous after BLAST, DuploMap refines placement using Paralogous Sequence Variants (PSVs)—the ~1–2% of positions that differ between 4q and 10q D4Z4 copies.

```bash
# Requires minimap2 alignment with secondary alignments retained
duplomap -i d4z4_candidates.bam \
  -r chm13v2.0.fa \
  -d sd_database_chm13.bed \
  -o d4z4_duplomap.bam
```

DuploMap iteratively identifies reliable PSVs within each read, discards PSVs overlapping polymorphisms, and computes per-read probabilistic assignment scores. In benchmarks on ONT data in segmental duplications, DuploMap improved correctly mapped reads from **75.5% to 88.2%** at MAPQ ≥10. For reads >50 kb, recall reached **94.9%**.

### 2d. Building a 4q/10q PSV catalog from the 100-sample panel

Process all 100 non-FSHD WGS samples through Steps 1–2b. For reads confidently assigned to 4q or 10q (by flanking anchor + BLAST), catalog all sequence variants within D4Z4 units. Positions where all 4q-assigned reads carry one allele and all 10q-assigned reads carry the other are **high-confidence PSVs**. Expected yield: **50–150 reliable PSVs** across the 3.3 kb repeat unit, based on the ~1–2% divergence between paralogs.

Store this PSV catalog as a VCF or BED file. In diagnostic samples, use PSV allele counts per read to compute a log-likelihood ratio for 4q vs 10q assignment—analogous to how HLA*LA uses its Population Reference Graph to distinguish HLA paralogs.

---

## Step 3: D4Z4 repeat unit counting

### 3a. Per-read repeat counting via BLAST

For each chromosome-assigned read, count D4Z4 repeat units by BLASTing a single 3.3 kb D4Z4 reference unit against the read and counting non-overlapping hits.

```bash
blastn -query d4z4_single_unit.fa \
  -subject <(samtools fasta -r READ_ID d4z4_candidates.bam) \
  -outfmt "6 qseqid sseqid pident length qstart qend sstart send" \
  -evalue 1e-50 \
  | awk '{print $2}' | sort | uniq -c
```

A read is classified as **"complete"** (spanning the full array) if it contains both the D4F104S1 proximal anchor and the pLAM/haplotype distal anchor. Only complete reads provide accurate total repeat counts. Reads containing D4Z4 + pLAM but not D4F104S1 are **"partial-distal"**—they provide a minimum repeat count.

### 3b. Alternative: tandem-genotypes with LAST

tandem-genotypes (v1.9.2) natively handles macrosatellite repeats and was demonstrated on D4Z4 by its original authors.

```bash
# LAST alignment pipeline
lastdb -P8 -uNEAR chm13_db chm13v2.0.fa
last-train -P8 chm13_db reads.fq > params.par
lastal -P8 -p params.par chm13_db reads.fq | last-split > aligned.maf

# Custom BED for D4Z4 (3.3kb motif)
echo -e "chr4\t193282466\t193558261\t3300" > d4z4_loci.txt
echo -e "chr10\t134615120\t134741148\t3300" >> d4z4_loci.txt

tandem-genotypes -u 3000 d4z4_loci.txt aligned.maf > d4z4_genotypes.txt
```

tandem-genotypes reports the change in copy number relative to the reference. Since T2T-CHM13v2.0 carries ~33 units on chr4, a reported deletion of 25 units indicates a sample with **~8 units** on that allele—in the FSHD1 pathogenic range.

### 3c. WGS-specific considerations for repeat counting

At 30× mean genome coverage, expect **3–10 reads** spanning into the D4Z4 array from each flanking side on each chromosome. For **short pathogenic arrays** (1–10 repeat units = 3.3–33 kb), ultra-long ONT reads (>50 kb) can span the entire array plus flanking anchors on both sides, enabling complete repeat counting. For **normal-length arrays** (>11 units = >36 kb), most WGS reads will be partial, providing only minimum repeat counts.

**Critical threshold:** If any chr4-assigned read shows a complete array of **≤10 D4Z4 units** on a 4qA haplotype, this is a candidate FSHD1 allele. The 100-sample non-FSHD panel establishes the expected distribution: healthy controls should show **≥11 units** on all 4qA alleles.

Sequencing depth recommendation: **≥60× mean genome coverage** with an ultra-long read protocol (e.g., SQK-ULK114 kit targeting N50 >25 kb) to maximize the number of array-spanning reads. At 60×, expect ~6–20 informative reads per D4Z4 allele.

### 3d. Assembly-based validation

For samples with sufficient ultra-long coverage:

```bash
# Extract D4Z4 region reads and perform local assembly
samtools view -b aligned.bam "chr4:193,200,000-193,600,000" > chr4_d4z4_reads.bam
samtools fastq chr4_d4z4_reads.bam > chr4_d4z4.fq

flye --nano-hq chr4_d4z4.fq \
  --out-dir chr4_d4z4_assembly \
  --genome-size 400k \
  --threads 16
```

Align the assembled contig(s) back to T2T-CHM13v2.0 and count D4Z4 units via dot-plot or BLAST. Tools like **Shasta** (ONT-optimized) or **Flye** can resolve repeats shorter than the read length. D4Z4End2End demonstrated full-length array assembly up to **42 repeat units (~140 kb)** from ultra-long ONT reads.

---

## Step 4: methylation analysis at D4Z4

### 4a. Region-specific methylation extraction

```bash
# Create D4Z4 target BED
cat > d4z4_methylation_targets.bed << 'EOF'
chr4	193540172	193543634	distal_D4Z4_DUX4_chr4
chr10	134738000	134741148	distal_D4Z4_chr10
EOF

# Per-site methylation pileup
modkit pileup aligned.bam d4z4_methylation.bed \
  --cpg \
  --combine-strands \
  --combine-mods \
  --ref chm13v2.0.fa \
  --include-bed d4z4_methylation_targets.bed \
  --filter-threshold 0.5

# Per-read methylation extraction (for allele-specific analysis)
modkit extract aligned.bam d4z4_per_read.tsv \
  --include-bed d4z4_methylation_targets.bed \
  --ref chm13v2.0.fa
```

The target region **chr4:193,540,172–193,543,634** (used by DUCKS4) covers the most distal D4Z4 repeat unit including the DUX4 promoter and coding region. This is the most diagnostically informative region: it captures both the **DR1 region** (29 CpGs; focal demethylation in FSHD2) and the **DUX4-PAS region** (10 CpGs; informative for FSHD1).

### 4b. Allele-specific methylation via read-level phasing

A key advantage of ONT over bisulfite sequencing: methylation is called per-read, enabling allele-specific analysis without additional experimental steps.

```bash
# Variant call for phasing (SNPs in flanking regions)
pepper_margin_deepvariant call_variant \
  --bam aligned.bam --fasta chm13v2.0.fa \
  --output_dir variants/ --ont_r10_q20 --threads 32

# Phase and haplotag reads
whatshap phase --reference chm13v2.0.fa \
  variants/output.vcf.gz aligned.bam -o phased.vcf.gz
whatshap haplotag --reference chm13v2.0.fa \
  phased.vcf.gz aligned.bam -o haplotagged.bam

# Haplotype-separated methylation
modkit pileup haplotagged.bam hap1_meth.bed \
  --cpg --haplotype 1 --ref chm13v2.0.fa \
  --include-bed d4z4_methylation_targets.bed
modkit pileup haplotagged.bam hap2_meth.bed \
  --cpg --haplotype 2 --ref chm13v2.0.fa \
  --include-bed d4z4_methylation_targets.bed
```

For samples with two 4qA alleles (one contracted, one normal), bulk bisulfite sequencing cannot separate them. ONT per-read phasing resolves this. Huang et al. (2024) demonstrated this advantage in a mosaic sample: the contracted 2-RU allele showed **12.15% methylation** versus **76.46%** for the non-pathogenic 37-RU allele—a distinction invisible to bisulfite.

### 4c. Diagnostically relevant methylation thresholds

| Metric | Healthy | FSHD1 | FSHD2 | Source |
|--------|---------|-------|-------|--------|
| DR1 region mean CpG methylation | >40% | 25–40% (contracted allele) | **<30%** (all alleles) | Hartweck 2013 |
| DUX4-PAS CpG6 methylation | >72% | **<72%** | **<72%** | Cascella 2018 |
| FseI-equivalent single CpG | 35–55% | Intermediate | **≤25%** | Lemmers 2015 |
| Delta1 score (size-corrected) | >−21% | Variable | **≤−21%** | Lemmers 2015 |
| DR1_CpG3 | >37% | Variable | **≤37%** | Megalizzi 2022 |

**FSHD2 methylation signature:** Hypomethylation affects **all four D4Z4 arrays** (both chr4 alleles + both chr10 alleles), distinguishing it from FSHD1 where hypomethylation is predominantly on the contracted allele only.

### 4d. Methylation gradient analysis

Butterfield et al. (2023) discovered that CpG methylation across D4Z4 arrays follows an **asymmetric length-dependent gradient**: proximal (centromeric) repeats are hypomethylated, with methylation increasing toward telomeric repeats and reaching a plateau at approximately **10 repeat units**—precisely the pathogenic threshold. Within each 3.3 kb unit, a ~180-nucleotide periodicity in the DUX4 coding region reflects phased nucleosome positioning. Reproducing this gradient from WGS data confirms correct array reconstruction.

---

## Step 5: haplotype determination (4qA vs 4qB)

### 5a. pLAM/polyadenylation signal detection

**FSHD requires the permissive 4qA haplotype.** The 4qA allele contains a 260 bp pLAM region immediately distal to the last D4Z4 repeat, harboring DUX4 exon 3 with the **ATTAAA** polyadenylation signal. The 4qB allele lacks pLAM entirely, making contractions on 4qB non-pathogenic.

For each chr4-assigned read extending distal to the D4Z4 array, BLAST the pLAM reference sequence:

```bash
blastn -query plam_4qA.fa -subject read.fa -outfmt 6
```

**Classification:**

- **pLAM hit with ATTAAA** → **4qA** (permissive; FSHD-relevant)
- **pLAM hit with ATCAAA** → **10qA** (non-permissive; confirms chr10 origin)
- **No pLAM hit** + β-satellite absence → **4qB** (non-permissive)
- **pLAM hit with ATTAAA** on chr10-assigned read → rare **10q-4qA translocation** (permissive on chr10; clinically relevant)

### 5b. Sub-haplotype classification

The DUCKS4 pipeline uses an identification key (`Haplotypes_identification_key.xlsx`) incorporating three markers:

- **SSLP** (Simple Sequence Length Polymorphism) in the CLUHP-4-201 region, ~3.5 kb proximal to D4Z4
- **D4F104S1 SNPs** (p13E-11 region variants)
- **Distal A/B polymorphism** (pLAM presence/absence + β-satellite)

The dominant FSHD-associated haplotype is **4A161** (~85% of FSHD cases), which has two variants: **4A161S** (standard) and **4A161L** (1.6 kb larger terminal D4Z4 unit with alternative DUX4 3′UTR exons 3a/3b). Both are permissive.

### 5c. Using the 100-sample panel for haplotype frequency estimation

Process all 100 non-FSHD samples through Steps 1–5b to build a haplotype frequency table. Expected distribution in European populations: ~50% 4qA, ~50% 4qB on chr4; ~96% 10qA166 on chr10. This panel also identifies population-specific rare haplotypes that could complicate interpretation.

---

## Step 6: FSHD2 gene variant screening

FSHD2 requires D4Z4 hypomethylation on a 4qA background plus a pathogenic variant in a chromatin modifier gene. Screen three genes:

```bash
# PEPPER-Margin-DeepVariant for SNV/indel calling
pepper_margin_deepvariant call_variant \
  --bam aligned.bam --fasta chm13v2.0.fa \
  --output_dir fshd2_variants/ \
  --ont_r10_q20 --threads 32

# Filter to FSHD2 genes (T2T-CHM13v2.0 coordinates)
bcftools view fshd2_variants/output.vcf.gz \
  --regions chr18:2,655,860-2,806,008 \
  > smchd1_variants.vcf   # SMCHD1

bcftools view fshd2_variants/output.vcf.gz \
  --regions chr20:31,350,193-31,397,163 \
  > dnmt3b_variants.vcf   # DNMT3B
```

**Note on variant calling scope:** Since Step 1 only extracted FSHD-relevant reads from the hg38 BAM, the realigned BAM contains reads at SMCHD1/DNMT3B/LRIF1 loci (included in the extraction BED) but not genome-wide variants. PEPPER-Margin-DeepVariant will call variants only in the covered regions. If full-genome variant calling is needed for other purposes, run it on the original hg38 BAM separately.

**SMCHD1** (chr18p11.32): mutated in ~80% of FSHD2 cases. >52 pathogenic mutations described—47% splice-site, 27.5% missense, 15.5% nonsense, 10% indels. Dominant-negative missense mutations in the ATPase domain cause more severe phenotypes than haploinsufficiency mutations. **DNMT3B** and **LRIF1** account for the remaining FSHD2 cases.

Annotate variants with VEP or ClinVar. Any heterozygous likely-pathogenic variant in SMCHD1/DNMT3B/LRIF1, combined with D4Z4 hypomethylation and a 4qA allele, supports FSHD2 classification.

---

## Step 7: diagnostic classification

### Decision tree

```
1. Does the sample have a chr4 4qA allele with ≤10 D4Z4 repeat units?
   ├─ YES → Candidate FSHD1
   │   └─ Is the contracted 4qA allele hypomethylated?
   │       ├─ YES → **FSHD1** (confirm with clinical phenotype)
   │       │   └─ 1–3 units: severe/infantile onset expected
   │       │   └─ 4–7 units: classic adult onset
   │       │   └─ 8–10 units: mild/reduced penetrance (~20% non-penetrant)
   │       └─ NO → Possible non-penetrant carrier or technical artifact → review
   │
   ├─ NO (all 4qA alleles have ≥11 units) →
   │   └─ 2. Is there genome-wide D4Z4 hypomethylation?
   │       (DR1 mean <30% AND/OR Delta1 ≤ −21% AND/OR
   │        all four D4Z4 arrays show concordant hypomethylation)
   │       ├─ YES → Candidate FSHD2
   │       │   └─ 3. Is there a pathogenic SMCHD1/DNMT3B/LRIF1 variant?
   │       │       ├─ YES + 4qA present → **FSHD2**
   │       │       └─ NO → Possible novel FSHD2 gene; report as
   │       │              "D4Z4 hypomethylation, unresolved etiology"
   │       └─ NO → **Unaffected**
   │
   └─ BORDERLINE (8–10 units, found in 1–2% of controls) →
       └─ Combined assessment: methylation + clinical phenotype + SMCHD1 status
          (SMCHD1 variants can act as disease modifiers in compound FSHD1+FSHD2)
```

### Classification output fields

For each sample, report:

| Field | Value |
|-------|-------|
| Sample ID | — |
| Chr4 allele 1: haplotype | 4qA / 4qB |
| Chr4 allele 1: D4Z4 repeat units | Integer (or ≥N if partial) |
| Chr4 allele 1: DR1 mean methylation | Percentage |
| Chr4 allele 2: haplotype | 4qA / 4qB |
| Chr4 allele 2: D4Z4 repeat units | Integer (or ≥N if partial) |
| Chr4 allele 2: DR1 mean methylation | Percentage |
| Chr10 allele 1: repeat units | Integer |
| Chr10 allele 2: repeat units | Integer |
| Genome-wide D4Z4 methylation (4 arrays) | Mean percentage |
| SMCHD1/DNMT3B/LRIF1 variants | VCF entries |
| Classification | FSHD1 / FSHD2 / Unaffected / Borderline |
| Confidence | High / Medium / Low |
| Complete spanning reads (chr4) | Count per allele |

**Confidence scoring:** "High" requires ≥3 complete spanning reads per allele and concordant methylation; "Medium" requires ≥1 complete read or ≥3 partial-distal reads; "Low" indicates no spanning reads (report minimum repeat count only).

---

## Leveraging the 100 non-FSHD reference panel

The 100 non-FSHD WGS samples serve four functions across the pipeline:

**Baseline D4Z4 distributions.** Process all 100 samples through Steps 1–5 to establish reference distributions for repeat count, methylation, and haplotype frequencies. Any diagnostic sample falling outside the 100-sample normal range triggers closer review. Expected healthy distributions: 11–100 repeat units per 4q allele, DR1 methylation >40%, DUX4-PAS CpG6 >72%.

**PSV catalog construction.** Reads confidently assigned to 4q or 10q (by flanking anchors) across 100 samples yield a robust catalog of **paralog-specific variants** within D4Z4. Positions where 4q reads unanimously carry one allele and 10q reads carry another are high-confidence PSVs. This catalog directly feeds the DuploMap-style disambiguation in Step 2.

**Sniffles2 population SV panel.** Generate per-sample `.snf` files for population-level SV calling:

```bash
# Per sample
sniffles --input sample.bam --snf sample.snf --reference chm13v2.0.fa

# Joint calling across panel
sniffles --input sample1.snf sample2.snf ... sample100.snf \
  --vcf population_svs.vcf --reference chm13v2.0.fa
```

The population VCF identifies common structural variants and D4Z4-flanking polymorphisms in the healthy population. Diagnostic samples can be force-genotyped against this panel to identify rare/novel variants.

**Custom methylation classifier training.** Using the 100-sample panel's per-CpG methylation values at D4Z4, train a lightweight classifier (logistic regression or random forest via scikit-learn) to predict FSHD status. Megalizzi et al. (2022) demonstrated that just **4 CpG sites** (DUX4-PAS_CpG6, DUX4-PAS_CpG3, DR1_CpG1, DR1_CpG22) achieve 94% accuracy, 93% sensitivity, and 96% specificity for distinguishing FSHD from controls. The 100-sample panel provides the negative-class training data; published FSHD methylation values from the literature supply the positive class.

---

## ML-based tools and their applicability

No existing ML-based SV caller is specifically designed for D4Z4 macrosatellite analysis. The tools most relevant to this pipeline fall into three categories:

**SV callers with ML components.** SVision (CNN for complex SV characterization), SAVANA (random forest for somatic SVs, explicitly retrainable via `savana train`), and dysgu (gradient-boosted trees) all use ML architectures but target general structural variants, not macrosatellite repeat counting. SVision's sequence-to-image CNN could in principle be adapted to classify D4Z4 array patterns, but this would require substantial custom training data.

**Repeat-specific tools with ML.** DeepRepeat uses a CNN on raw nanopore signal for STR counting but requires legacy Albacore basecalling (incompatible with R10.4.1/Dorado). HMMSTR uses profile HMMs and supports ONT R10 data, but targets short tandem repeats (1–6 bp motifs), not 3.3 kb macrosatellites. **LongTR** (Gymrek Lab, 2024) supports ONT and outperforms other tools at long VNTRs (>500 bp), but its upper motif-size limit has not been validated for 3.3 kb units. **NanoSatellite** uses dynamic time warping on raw signal and was validated on the ABCA7 VNTR (25 bp motif, arrays up to 10 kb)—closer to D4Z4's scale but still an order of magnitude smaller.

**The most effective approach for D4Z4 remains alignment-based counting** (BLAST hit counting per read, as in DUCKS4 and Huang et al. 2024) combined with assembly-based validation. The BLAST approach is simple, validated, and directly interpretable. For ML augmentation, training a random forest on features extracted from reads at D4Z4 (repeat count, methylation level, PSV pattern, flanking sequence identity) using the 100-sample panel as negative examples could improve borderline classifications.

---

## Existing validated ONT FSHD pipelines

Three published pipelines directly inform this design:

**DUCKS4** (Löwenstern et al. 2026, Human Genomics) is the most comprehensive automated pipeline. Docker-based, it takes Dorado-basecalled BAM files, aligns to T2T-CHM13v2.0, uses NCBI BLAST for chromosome/haplotype assignment, counts D4Z4 repeats per read, calls methylation with modkit, screens FSHD2 genes, and detects mosaicism. Available at github.com/tamara-nano/ducks4. It was designed for both Cas9-targeted and adaptive-sampling data but can process WGS reads.

**D4Z4End2End** (Xiao et al. 2025, Genome Research) provides complete genetic and epigenetic D4Z4 characterization. It resolves full-length arrays up to 42 repeat units, generates consensus sequences with Racon, and uses NanoMethViz for methylation visualization. It was applied to 1000 Genomes Project and HPRC data, demonstrating feasibility with publicly available WGS datasets. Available at github.com/lucindaxiao/D4Z4End2End.

**Huang et al. (2024, J Translational Medicine)** demonstrated ONT WGS at ~29× coverage on R9.4.1 with minimap2 alignment to T2T-CHM13v2.0. They achieved **100% concordance with OGM** for pathogenic D4Z4 contractions across 16 FSHD1 cases and showed that ONT allele-specific methylation outperformed bisulfite sequencing (DUX4 upstream AUC = **0.996** vs bisulfite AUC = 0.943). This is the strongest validation of WGS-only ONT for FSHD.

---

## Practical pipeline summary

| Step | Tool(s) | Key output |
|------|---------|-----------|
| Read extraction (hg38) | samtools view/merge (region + unmapped) | FSHD-subset BAM (~2–5% of total) |
| Revert to unaligned | samtools reset --keep-tag MM,ML | Unaligned modBAM |
| Align to T2T | Winnowmap2 v2.03 (or minimap2 ≥v2.22) → T2T-CHM13v2.0 | Sorted, indexed BAM |
| D4Z4 read extraction | samtools view (region filter) | Candidate BAM |
| 4q/10q assignment | BLAST (feature sequences) + DuploMap (PSV refinement) | Per-read chromosome label |
| Repeat counting | BLAST (D4Z4 unit hits per read) + tandem-genotypes | Repeat units per allele |
| Haplotype calling | BLAST (pLAM/β-satellite) | 4qA / 4qB / 10qA per allele |
| Methylation | modkit pileup + extract | Per-CpG methylation frequencies per allele |
| FSHD2 screening | PEPPER-Margin-DeepVariant → VEP | SMCHD1/DNMT3B/LRIF1 variants |
| Classification | Decision tree + panel comparison | FSHD1 / FSHD2 / Unaffected |

---

## Known limitations of WGS-only approach

**Coverage at D4Z4 is the primary bottleneck.** Without enrichment, expect ~5–15× effective coverage at D4Z4 versus 50–500× with Cas9 targeting. This limits confidence for allele-specific methylation calls and repeat counting of large normal arrays. Mitigation: sequence to ≥60× mean genome coverage; use ultra-long read protocols (SQK-ULK114) to maximize array-spanning reads.

**Normal-length arrays (≥11 units) are difficult to fully span.** An 11-unit array is ~36 kb; a 30-unit array is ~99 kb. Standard ONT read N50 of 10–25 kb cannot span these. Only ultra-long reads (N50 >50 kb) can fully span shorter normal arrays. For most normal arrays, the pipeline reports a minimum repeat count rather than an exact count.

**Hybrid/translocated arrays (~10–20% of individuals)** carrying 10q-type units on chr4 or vice versa will show mixed BlnI/XapI patterns. DUCKS4 flags these as requiring manual review. The BLAST bitscore approach handles most cases, but complex translocations may need orthogonal confirmation.

**FSHD2 with normal repeat count but subtle methylation changes** is harder to detect at low coverage, since methylation confidence intervals widen. The 4-CpG classifier (DR1_CpG3 ≤37%, DUX4-PAS CpG6 <72%) requires ≥5 reads per allele for reliable calls.

Despite these limitations, Huang et al. (2024) demonstrated that standard ONT WGS at ~29× coverage achieves 100% concordance with optical genome mapping for detecting pathogenic D4Z4 contractions, validating this approach for clinical FSHD1 diagnosis. FSHD2 detection requires additional methylation depth and should be interpreted with appropriate confidence intervals when coverage is limited.