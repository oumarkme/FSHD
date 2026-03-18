# nf-fshd: ONT D4Z4 Repeat Unit Counting, Haplotyping, and Methylation Pipeline

A Nextflow DSL2 pipeline for FSHD diagnosis from hg38-aligned ONT BAM/CRAM files with HP tags from upstream phasing. The pipeline resolves all four D4Z4 alleles (2× chr4, 2× chr10), determines 4qA/4qB haplotype, counts D4Z4 repeat units per allele, and reports allele-specific CpG methylation at the distal DUX4 region. All analysis coordinates reference **T2T-CHM13v2.0**.

---

## Biological background

### D4Z4 and FSHD

The D4Z4 macrosatellite array on chromosome 4q35 consists of tandem 3.3 kb repeat units. Healthy individuals carry 11–100+ RU per allele. FSHD1 occurs when one allele contracts to 1–10 RU **on a permissive haplotype background** (4qA or 4qA-L), which stabilizes DUX4 mRNA via a polyadenylation signal (PAS) in the pLAM region distal to the array.

### Haplotype backgrounds

| Haplotype | Permissive? | Key feature | Population freq |
|-----------|-------------|-------------|-----------------|
| **4qA** (includes 4A161) | Yes | Contains pLAM with functional DUX4 PAS (ATTAAA) | ~50% of chr4 alleles |
| **4qA-L** (4A161L) | Yes | 1.6 kb longer distal partial D4Z4 unit; functional PAS | ~20% of European 4A161 alleles |
| **4qB** | No | Lacks pLAM entirely; DUX4 mRNA unstable | ~50% of chr4 alleles |
| **10qA** | No | 98.5% identity to 4qA but non-functional PAS (SNPs) | Most chr10 alleles |

Only contractions on 4qA or 4qA-L cause FSHD1. Contractions on 4qB are clinically silent.

### D4Z4 array sizes

| RU count | Array size | Clinical significance |
|----------|-----------|----------------------|
| 1–3 | 3.3–10 kb | Severe/infantile FSHD1 (high penetrance) |
| 4–7 | 13–23 kb | Classic FSHD1 (variable expressivity) |
| 8–10 | 26–33 kb | Borderline; reduced penetrance (~40–50%); requires methylation |
| 11–20 | 36–66 kb | Normal; may contribute to FSHD2 if modifier present |
| 21–100+ | 69–330+ kb | Normal/healthy |

---

## Technical constraints (ONT sequencing)

### Read length vs D4Z4 spanning

| Parameter | Standard ligation (LSK114, with size selection) | Ultra-long (ULK114) |
|-----------|------------------------------------------------|---------------------|
| Typical N50 | 25–30 kb | 50–100+ kb |
| Max reliable spanning | ~7–8 RU (~26 kb) | ~15–20 RU |
| Marginal spanning | 9–10 RU (~33 kb) | 20–30 RU |
| Never spanned | ≥11 RU (≥36 kb) | ≥30 RU |

A read is **"complete"** when it contains both the proximal p13E-11 anchor AND the distal pLAM anchor. With standard N50 ~25–30 kb, complete reads can only come from alleles with ≤~10 RU. Healthy alleles (≥11 RU) will never produce complete reads with standard prep.

### Coverage per allele

| Genome coverage | Reads per allele (theoretical) | Effective (subtelomeric penalty) |
|-----------------|-------------------------------|----------------------------------|
| 30× | ~7–8× | ~4–6× |
| 60× | ~15× | ~10–12× |

### HP tag reliability

- HP tags from WhatsHap are reliable in unique flanking regions adjacent to D4Z4
- Within the D4Z4 array: very few heterozygous SNPs → phasing may fail
- Most D4Z4 reads carry HP tags inherited from flanking-region phasing, but some reads (especially those entirely within the array) may lack HP tags

### Chr4 vs chr10 discrimination

- D4Z4 arrays share ~98.5–99% identity between chr4 and chr10
- Discrimination relies on unique flanking sequences and a chr10-specific BlnI restriction site within D4Z4 units
- Reads must extend into flanking regions for reliable chromosome assignment
- Haplotype detection (4qA vs 10qA) requires reads extending past the last D4Z4 unit to detect pLAM

### Flye assembly caveats

- Minimum ~10× coverage; 30–40× recommended. Per-allele coverage at D4Z4 is typically 4–12× (marginal)
- Flye may collapse identical D4Z4 copies, underestimating RU count
- Resolution: take whichever has MORE RUs — assembly or longest raw read

---

## Algorithm: complete-reads-first allele assignment

The pipeline processes chr4 and chr10 independently. After BLAST annotation, every read has: `RU_Count`, `Has_Head` (p13E-11), `Has_Tail` (pLAM), `Haplotype` (4qA/10qA/−), and `HP` tag (1/2/none).

The algorithm has three sequential phases:

```
Phase 1: SEPARATE ALLELES     → which reads belong together?
Phase 2: FIND COMPLETE READS  → which allele has spanning reads?
Phase 3: COUNT per allele     → read-level (consensus) or assembly?
```

### Phase 1: Separate alleles

Group reads into two alleles using the best available signal:

```
IF any reads on this chromosome have HP tags:
    Allele_1 = reads with HP == 1
    Allele_2 = reads with HP == 2
    unphased = reads without HP tag
    separation_method = 'hp'

ELSE (no HP tags at all):
    IF complete reads exist:
        complete_ru = median(RU_Count of complete reads)
        Allele_1 = reads with RU_Count <= complete_ru + 2
        Allele_2 = reads with RU_Count > complete_ru + 2
        separation_method = 'ru-anchor'
    ELSE:
        IF std(RU_counts) >= 1.5 AND len(reads) >= 2:
            Jenks 2-class split -> Allele_1 (shorter), Allele_2 (longer)
            separation_method = 'jenks'
        ELSE:
            Allele_1 = all reads
            Allele_2 = empty
            separation_method = 'single'
```

HP is always preferred when available — it uses flanking-region SNPs and doesn't depend on RU count differences.

### Phase 2: Assign unphased reads

Complete reads (Has_Head AND Has_Tail) span the full array and **can only come from the shorter allele** (biological constraint: healthy alleles ≥11 RU exceed ONT N50).

```
IF no complete reads:
    Skip to Phase 3 (all alleles use assembly)

IF complete reads already have HP tags:
    Already in the correct allele from Phase 1

IF complete reads lack HP tags (separation was by HP):
    Assign to the allele whose median RU is closest

All remaining unphased reads:
    IF allele has complete reads (max_RU = X) AND read RU > X + 2:
        Force-assign to the OTHER allele (ru-constraint)
    ELSE:
        minimap2 align read against HAP1 pool and HAP2 pool
        Assign to allele with more matching bases (alignment)
        Tie: assign to allele with closer median RU (ru-proximity)
```

### Phase 3: Per-allele strategy

Each allele independently gets one of two strategies:

```
For each allele:

    IF allele has complete reads -> STRATEGY: READ-LEVEL
        spanning = complete + partial-distal reads
        raw_ref = longest spanning read

        IF >= 3 spanning reads:
            Racon consensus (raw_ref as seed, polish with others)
        ELSE:
            consensus = raw_ref as-is

        reference = consensus
        RU_count = BLAST d4z4.fa against consensus (hits > 1500 bp)

    ELSE -> STRATEGY: ASSEMBLY
        Run Flye on allele reads (--nano-hq --genome-size 300k)
        max_raw_ru = max(RU_Count) among reads in this allele

        IF Flye succeeded AND Flye_RU >= max_raw_ru:
            reference = Flye contig
        ELSE:
            reference = longest raw read

        RU_count = BLAST d4z4.fa against reference (hits > 1500 bp)
```

### Example walkthroughs

**Classic FSHD1 (most common):**

| Phase | What happens |
|-------|-------------|
| Phase 1 | HP tags exist → Allele_1 (HP=1), Allele_2 (HP=2) |
| Phase 2 | Complete reads found with HP=1 → already in Allele_1 |
| Phase 3 | Allele_1: read-level (Racon consensus). Allele_2: assembly (Flye vs max raw read) |

**No HP tags, complete reads exist:**

| Phase | What happens |
|-------|-------------|
| Phase 1 | No HP → use complete reads as RU anchor. Allele_1 (≤RU+2), Allele_2 (>RU+2) |
| Phase 2 | Complete reads already in Allele_1 by construction |
| Phase 3 | Allele_1: read-level. Allele_2: assembly (or empty → single allele) |

**Both alleles short (homozygous contraction):**

| Phase | What happens |
|-------|-------------|
| Phase 1 | HP tags → Allele_1 (HP=1), Allele_2 (HP=2) |
| Phase 2 | Complete reads in BOTH alleles |
| Phase 3 | Both alleles: read-level (Racon consensus each) |

**No complete reads (chr10, or healthy chr4):**

| Phase | What happens |
|-------|-------------|
| Phase 1 | HP tags → Allele_1 (HP=1), Allele_2 (HP=2) |
| Phase 2 | No complete reads → skip |
| Phase 3 | Both alleles: assembly (Flye vs max raw read) |

---

## Pipeline steps

### Step 1: Extract FSHD-region reads (EXTRACT_FSHD_READS)

Extract reads from FSHD-associated regions in hg38 coordinates plus unmapped reads.

```
BED (hg38):
    chr4:189000000-190500000     # 4q35 D4Z4 flanking
    chr10:133000000-134000000    # 10q26 D4Z4 flanking

samtools view -b -h -M -L regions.bed -> target.bam
samtools view -b -f 4 -> unmapped.bam
samtools merge -> subset.bam
```

### Step 2: Realign to T2T-CHM13v2.0 (ALIGN_T2T)

Realign to T2T, preserving MM, ML, HP, PS tags via `samtools fastq -T` and `minimap2 -y`.

`minimap2 -ax lr:hq -y --eqx --secondary=no | samtools sort`

### Step 3: Extract D4Z4 candidate reads (EXTRACT_D4Z4)

```
BED (T2T):
    chr4:193230000-193610000
    chr10:134560000-134800000
```

Published output: `{SAMPLE_ID}.d4z4.bam` + `.bai`

### Step 4: Chromosome assignment and feature annotation (BLAST_ANNOTATE)

```
blastn d4z4.fa       -> RU counting (hits > 1500 bp)
blastn pLAM-polyA.fa -> distal anchor + 4qA/10qA haplotype
blastn P13E11.fasta  -> proximal anchor
```

**Chromosome assignment (three-tier):**

1. **BAM alignment** (initial): minimap2 alignment to T2T reference
2. **BlnI motif counting** (correction): chr4-type `CTAGGTCTAGG` vs chr10-type `CTAGGCCTAGG` — majority wins
3. **Haplotype override** (final): 4qA (ATTAAA) → force chr4; 10qA (ATCAAA) → force chr10

### Step 5: Allele assignment (ALLELE_ASSIGN)

Three-phase algorithm: Phase 1 (separate by HP/RU-anchor/Jenks), Phase 2 (assign unphased reads via RU-constraint + alignment + proximity), Phase 3 (determine read-level vs assembly strategy per allele).

### Step 6: Read-level counting (READLEVEL_COUNT)

For alleles with complete reads: extract spanning reads → Racon consensus (if ≥3 reads) → BLAST RU count.

### Step 7: Assembly-based counting (ASSEMBLY_COUNT)

For alleles without complete reads: Flye assembly → compare Flye RU with max raw-read RU → select reference with more RUs → BLAST RU count.

### Step 8: Per-allele CpG methylation (METHYLATION)

```
modkit pileup --phased --cpg --combine-strands --combine-mods
    --modified-bases 5mC
    --include-bed chr4:193540172-193543634, chr10:134738000-134741148
```

Runs independently; only needs T2T BAM from Step 2.

### Step 9: Generate report (REPORT)

Produces `repeats.tsv` and `methylation.tsv` from allele assignments, RU counts, and methylation data.

---

## Pipeline flow diagram

```
Input: hg38 BAM/CRAM (HP/PS + MM/ML tags)
    |
    v
Step 1: EXTRACT_FSHD_READS (hg38 regions)
    |
    v
Step 2: ALIGN_T2T (preserve MM/ML/HP/PS) -------.
    |                                            |
    v                                            |
Step 3: EXTRACT_D4Z4                             |
    |                                            |
    v                                            |
Step 4: BLAST_ANNOTATE                           |
(RU, anchors, BlnI chr, haplotype)               |
    |                                            |
    v                                            |
Step 5: ALLELE_ASSIGN  <-------------------------'
    |
    |  Phase 1: Separate alleles
    |      HP tags? -> group by HP
    |      No HP?   -> RU-anchor (if complete reads) or Jenks
    |
    |  Phase 2: Assign complete reads
    |      Already have HP? -> done
    |      No HP? -> assign to allele with closest median RU
    |
    |  Phase 3: Per-allele strategy
    |      Has complete reads? -> Step 6 (read-level)
    |      No complete reads?  -> Step 7 (assembly)
    |
    |               .-----------------.
    v               v                 |
Step 6:         Step 7:               |
READLEVEL       ASSEMBLY              |
(Racon          (Flye vs raw read     |
 consensus)      -> pick more RUs)    |
    |               |                 |
    '-------+-------'                 |
            |                         |
            |         .---------------'
            v         v
        Step 8: METHYLATION
        (modkit --phased from T2T BAM)
            |
            v
        Step 9: REPORT
        (repeats.tsv + methylation.tsv)
```

---

## Output

Four published files per sample:

| File | Description |
|------|-------------|
| `{SAMPLE_ID}.d4z4.bam` | D4Z4 candidate reads aligned to T2T-CHM13v2.0 |
| `{SAMPLE_ID}.d4z4.bam.bai` | BAM index |
| `{SAMPLE_ID}.repeats.tsv` | D4Z4 repeat unit counts per allele |
| `{SAMPLE_ID}.methylation.tsv` | CpG methylation at distal D4Z4/DUX4 region |

### repeats.tsv

One row per (Allele, Haplotype, Complete) group. RU_Count is the max within each group.

| Column | Description |
|--------|-------------|
| Chromosome | chr4 or chr10 |
| Allele | HAP1 / HAP2 / − (unphased) |
| Haplotype | 4qA / 10qA / − (no pLAM = 4qB or untyped) |
| RU_Count | Max RU count among reads in this group |
| Complete | complete / head-only / tail-only / incomplete |
| Count_Method | read-level / assembly / − (unphased) |
| Reference_Source | consensus / flye / raw-read / − (unphased) |
| Read_Count | Number of reads in this group |

### methylation.tsv

| Column | Description |
|--------|-------------|
| Haplotype | HP1, HP2, or Combined |
| Region | DR1 or DUX4-PAS |
| Methylation | Mean % modified across CpGs in region |
| N_CpG | Number of CpGs with data |
| Ref_Atlas (Mean (SD) %) | ONT methylation reference atlas (~45–49 healthy controls) |
| Lit_Caputo (%) | Diagnostic threshold from Caputo et al. 2022 |

**Reference values (T2T coordinates):**

| Region | T2T coordinates | Ref_Atlas | Lit_Caputo threshold |
|--------|----------------|-----------|---------------------|
| DR1 | chr4:193,541,000–193,541,600 | 54.70% (1.97%) | CpG1 ≤ 53% |
| DUX4-PAS | chr4:193,543,430–193,543,634 | 51.81% (0.23%) | CpG6 ≤ 78% |

---

## Usage

```bash
bash nf-fshd/run.sh /path/to/sample.bam
```

The `run.sh` script auto-detects the environment, builds the Apptainer container if needed, and launches Nextflow with `-resume`.

**Input:** hg38-aligned BAM/CRAM with HP/PS tags (from WhatsHap) and MM/ML tags (from Dorado).

---

## Upstream prerequisites

These steps are performed outside this pipeline:

| Step | Tool | Output |
|------|------|--------|
| Basecalling | Dorado ≥0.7.0 (SUP + `--modified-bases 5mCG_5hmCG`) | Unaligned modBAM |
| Genome-wide alignment | minimap2 | hg38-aligned BAM |
| SNP calling | Clair3 ≥1.0.10 | VCF |
| Phasing + haplotagging | WhatsHap ≥2.0 | BAM with HP/PS tags |

---

## Reference genome and query sequences

**Reference:** T2T-CHM13v2.0 (`GCF_009914755.1_T2T-CHM13v2.0_genomic.fna`)

**Key T2T-CHM13v2.0 coordinates:**

| Region | Coordinates | Notes |
|--------|------------|-------|
| Chr4 D4Z4 array | chr4:193,282,466–193,558,261 | ~276 kb; 33 RUs in reference |
| Chr4 DUX4 methylation target | chr4:193,540,172–193,543,634 | Distal D4Z4 + DUX4 promoter/coding |
| Chr4 p13E-11 anchor | chr4:192,667,301–192,902,247 | Upstream unique anchoring region |
| Chr10 D4Z4 array | chr10:134,615,120–134,741,148 | ~126 kb |

**Query sequences (`nf-fshd/query/`):**

| File | Purpose |
|------|---------|
| `d4z4.fa` | Single D4Z4 RU (~3.3 kb) for RU counting |
| `pLAM-polyA.fa` | hapA (4qA ATTAAA PAS) / hapB (10qA ATCAAA PAS) for haplotype + distal anchor; 4qB has no pLAM |
| `P13E11.fasta` | p13E-11 probe for proximal anchor |

---

## Software

### In container (`environment.sif`)

| Tool | Version | Purpose |
|------|---------|---------|
| Python | 3.12 | Scripting (pandas, pysam, jenkspy) |
| samtools | 1.23 | BAM manipulation |
| minimap2 | 2.30 | Alignment |
| BLAST+ | 2.17.0 | RU counting, anchor detection, haplotyping |
| modkit | 0.6.1 | CpG methylation |
| Flye | 2.9.6 | Local assembly |
| Racon | 1.5.0 | Consensus polishing |
| seqkit | 2.13.0 | FASTA manipulation |
| jenkspy | 0.4.1 | Jenks natural breaks (fallback allele separation) |

---

## Known limitations

1. **Read length constraint.** At N50 ~25–30 kb, complete reads are limited to alleles ≤~10 RU. Ultra-long prep extends this range but is not standard.

2. **Subtelomeric coverage.** Effective per-allele coverage at D4Z4 is 30–50% lower than genome-wide average. At 30× WGS, expect ~4–6 reads per allele.

3. **Flye tandem repeat collapse.** Flye may underestimate RU count for large arrays. The pipeline mitigates this by comparing Flye RU with max single-read RU and selecting the larger.

4. **Methylation coverage.** At 30× WGS, per-allele methylation coverage may be only 3–5 reads. ≥60× recommended for clinical-grade calls.

5. **CHM13 reference gap for 4qA.** T2T-CHM13v2.0 carries 4qB on chr4. 4qA-specific pLAM/PAS is absent from the reference. BLAST-based annotation detects pLAM in read sequence directly.

6. **Somatic mosaicism.** Mosaic samples may show multiple RU count clusters from the same allele. The pipeline does not explicitly detect mosaicism but unusual clustering patterns should be flagged for review.

7. **Cis duplications.** Two adjacent D4Z4 arrays on one allele (~1.5% of Europeans) may produce two distinct RU groups in one HP. The pipeline may misinterpret this as two alleles.

---

## References

- Caputo V et al. (2022) Methylation thresholds for FSHD diagnosis. *Clin Epigenetics*.
- Caputo V et al. (2024) Updated FSHD methylation analysis.
- Lemmers RJLF et al. (2015) D4Z4 methylation and FSHD2. *Hum Mutat*.
- Löwenstern A et al. (2026) DUCKS4 pipeline. *Hum Genomics*. github.com/tamara-nano/ducks4
- Xiao L et al. (2025) D4Z4End2End. *Genome Res*. github.com/lucindaxiao/D4Z4End2End
- Huang Y et al. (2024) ONT WGS for FSHD with 100% OGM concordance. *J Transl Med*.
