# FSHD D4Z4 Repeat Analysis Pipeline (v3)

ONT long-read WGS pipeline for FSHD D4Z4 analysis using R10.4.1 chemistry. Extracts FSHD-relevant reads from hg38 BAM, realigns to T2T-CHM13v2.0, disambiguates 4q/10q, counts D4Z4 repeat units per read, and reports methylation percentage at the distal D4Z4/DUX4 region. FSHD2 modifier gene screening (SMCHD1/DNMT3B/LRIF1) is performed separately outside this pipeline. Full design rationale in `doc/fshd_info.md`.

## Server Environment

Both machines share the same mounted filesystem. Development on WSL, production runs on GS4.

| | WSL (dev) | GS4 (production) |
|---|---|---|
| Access | local | `ssh gs4` |
| CPUs | 8 (limit: 5) | 96 (limit: 67) |
| RAM | 30 GB (limit: 21 GB) | 1.7 TB (limit: 1.2 TB) |
| GPU | — | 4× NVIDIA RTX 4090 |
| Conda env | `/home/mark/miniforge3/envs/mark/bin` | `/home/user/miniforge-pypy3/envs/mark/bin` |
| Nextflow | In conda env above | In conda env above |
| Apptainer | In conda env above | In conda env above |

CPU and RAM usage capped at **70%** on both machines.

### Storage

| Path | Location | Purpose |
|---|---|---|
| `/lab_nas_home_user/` | NAS (mounted) | Code, reference genomes, pipeline output |
| `/gene_lv/` | NAS (mounted) | Input CRAM files |
| `/raid/nf/tmp` | Local disk (fast) | Nextflow work directory (intermediate files) |

### Directory Structure

```
/lab_nas_home_user/Mark/ONT/FSHD/                  # Pipeline output root
└── {SAMPLE_ID}/                                    # e.g. WGS-ONT019-5657831-D23881/
    ├── {SAMPLE_ID}.d4z4.bam                        # D4Z4 subset BAM aligned to T2T
    ├── {SAMPLE_ID}.d4z4.bam.bai                    # BAM index
    ├── {SAMPLE_ID}.d4z4_counts.tsv                 # RU counts per read
    ├── {SAMPLE_ID}.methylation.bed                 # Per-CpG methylation (bedMethyl)
    └── {SAMPLE_ID}.report.tsv                      # Summary: RU counts + methylation %
```

Nextflow work directory: `/raid/nf/tmp` (local SSD on both WSL and GS4, cleaned after successful runs).

## Pipeline Overview

```
hg38 BAM → extract FSHD regions + unmapped → revert to unaligned
→ minimap2 to T2T-CHM13v2.0 → extract D4Z4 candidates
→ BLAST 4q/10q disambiguation → BLAST D4Z4 RU counting
→ modkit methylation pileup → report (RU counts + methylation %)
```

## Pipeline Steps

### Step 1: Targeted read extraction + T2T realignment

Extract FSHD-relevant reads from hg38-aligned BAM (~2–5% of total data), revert to unaligned, and realign to T2T-CHM13v2.0 where the D4Z4 arrays are fully resolved.

**1a. Define hg38 extraction regions:**

| Region | hg38 Coordinates | Purpose |
|---|---|---|
| 4q35 D4Z4 + flanking | chr4:189,900,000–190,215,000 | D4Z4 array + flanking (includes ~82 kb assembly gap) |
| 10q26 D4Z4 + flanking | chr10:133,500,000–133,798,000 | 10q D4Z4 homolog (poorly assembled in hg38) |
| D4Z4-like paralogs | chr3, chr10p, chr18q, acrocentric chr | Reads misassigned to these paralogs under hg38 |

Unmapped reads are also extracted — the hg38 gap at 4q35 causes many D4Z4 reads to be unmapped.

```bash
samtools view -b -h -M -L fshd_regions_hg38.bed --fetch-pairs input.bam -o target.bam
samtools view -b -f 4 input.bam -o unmapped.bam
samtools merge -o fshd_subset.bam target.bam unmapped.bam
```

**1b. Revert to unaligned BAM preserving methylation tags:**

```bash
samtools reset --keep-tag MM,ML fshd_subset.bam -o unaligned.bam
```

Requires samtools ≥1.17. Strips alignment info (CIGAR, MAPQ, coordinates) while preserving MM/ML modification tags. Handles reverse-complement orientation reversal correctly.

**1c. Align to T2T-CHM13v2.0 with minimap2:**

```bash
minimap2 -ax lr:hq -y --eqx --secondary=no \
  -t ${CORES} chm13v2.0.fa unaligned.bam \
  | samtools sort -@ 8 -o aligned.bam
samtools index aligned.bam
```

| Flag | Purpose |
|---|---|
| `-ax lr:hq` | ONT R10.4.1 Q20+ preset (not `map-ont` which is for older/noisier data) |
| `-y` | Propagates MM/ML methylation tags from input to output |
| `--eqx` | Emits `=`/`X` CIGAR operators for exact match/mismatch (useful for PSV analysis) |
| `--secondary=no` | Suppresses secondary alignments (recommended for methylation analysis) |

**Key T2T-CHM13v2.0 coordinates:**

| Region | Coordinates | Notes |
|---|---|---|
| Chr4 D4Z4 + flanking | chr4:193,282,466–193,558,261 | ~276 kb; 33 RUs in reference |
| Chr4 DUX4 methylation target | chr4:193,540,172–193,543,634 | Distal D4Z4 + DUX4 promoter/coding |
| Chr4 p13E-11 anchor | chr4:192,667,301–192,902,247 | Upstream unique anchoring region |
| Chr10 D4Z4 array | chr10:134,615,120–134,741,148 | ~126 kb |

**Reference genomes:**

| Genome | Path |
|---|---|
| hg38 | `/lab_nas_home_user/Mark/reference_genome/hg38/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set_modified.fasta` |
| T2T-CHM13v2.0 | `/lab_nas_home_user/Mark/reference_genome/t2t/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna` |

Both use chr-prefixed naming (chr1, chr2, ...). The T2T reference is the NCBI version but has standard chromosome names — compatible with all coordinates in this document.

---

### Step 2: 4q vs 10q disambiguation

The D4Z4 arrays on 4q35 and 10q26 share >98% sequence identity. Classification follows the DUCKS4 method (Löwenstern et al. 2026): BLAST chromosome-specific feature sequences against each candidate read and assign by percent identity comparison.

**2a. Extract D4Z4 candidate reads:**

```bash
samtools view -b aligned.bam \
  "chr4:193,230,000-193,610,000" \
  "chr10:134,560,000-134,800,000" \
  > d4z4_candidates.bam
samtools fasta d4z4_candidates.bam > d4z4_reads.fa
```

**2b. Build feature BLAST database and classify reads:**

Feature sequences for the BLAST database (following DUCKS4):

| Feature | Variants | Purpose |
|---|---|---|
| D4F104S1 (p13E-11 probe) | 4qA, 4qB, 10q26 | ~15 SNPs differentiate haplotypes; highest identity determines chromosome |
| D4Z4 repeat unit | chr4 (3,298 bp), chr10 (3,309 bp) | Chromosome-specific D4Z4 unit sequences from T2T |
| pLAM | 4qA, 4qB, 10q26 | Distal sequence with diagnostic PAS variants |
| BlnI motif (CCTAGG) | — | Present in 10q D4Z4 (B+), absent in 4q (B−) |
| XapI motif (AAATTCC) | — | Present in 4q D4Z4 (X+), absent in 10q (X−) |

```bash
makeblastdb -in d4z4_features.fa -dbtype nucl -out features_db
blastn -query d4z4_reads.fa -db features_db \
  -outfmt "6 qseqid sseqid pident length bitscore" \
  -evalue 1e-5 -num_threads ${CORES} \
  > blast_hits.tsv
```

**2c. Classification logic per read:**

For each read, compare BLAST percent identity against each feature variant (4q vs 10q). The variant with highest identity determines chromosome assignment. Confirm with restriction site motif pattern:
- **Canonical 4q**: B−X+ (no BlnI site, has XapI site)
- **Canonical 10q**: B+X− (has BlnI site, no XapI site)
- **Mixed B−X−**: Hybrid/translocated array — flag for manual review

Reads that cannot be confidently assigned are excluded from downstream analysis.

---

### Step 3: D4Z4 repeat counting

**3a. Per-read repeat counting via BLAST:**

BLAST a full 3.3 kb D4Z4 reference unit against each chromosome-assigned read. Count hits with alignment length >1,500 bp as valid repeat units (DUCKS4 method).

```bash
blastn -query d4z4_single_unit.fa -db reads_db \
  -outfmt "6 qseqid sseqid pident length sstart send" \
  -evalue 1e-50 -num_threads ${CORES} \
  > d4z4_counts.tsv
```

Filter: keep hits where `length > 1500`. Count remaining hits per read = D4Z4 repeat units.

**3b. Read completeness classification:**

| Category | Definition | Use |
|---|---|---|
| **Complete** | Contains both D4F104S1 (proximal) AND pLAM (distal) anchors | Accurate total repeat count |
| **Partial-distal** | Contains D4Z4 + pLAM but not D4F104S1 | Minimum repeat count |
| **Partial-proximal** | Contains D4F104S1 + D4Z4 but not pLAM | Minimum repeat count |
| **Incomplete** | D4Z4 only, no flanking anchors | Excluded from sizing |

Only **complete reads** provide accurate total repeat counts. At 30× WGS coverage, expect 3–10 reads spanning into the D4Z4 array from each side per chromosome. Short pathogenic arrays (1–10 RU = 3.3–33 kb) can be fully spanned by ultra-long ONT reads (>50 kb).

**3c. Permissive haplotype check:**

For chr4-assigned reads containing pLAM, check the polyadenylation signal at chr4:193,543,619:
- **ATTAAA** → 4qA (permissive; FSHD-relevant)
- **ATCAAA** → 10qA (confirms chr10 origin)
- **No pLAM** + beta-satellite → 4qB (non-permissive)

---

### Step 4: Methylation analysis

**4a. Per-site methylation pileup with modkit:**

```bash
modkit pileup aligned.bam d4z4_methylation.bed \
  --cpg \
  --combine-strands \
  --combine-mods \
  --ref chm13v2.0.fa \
  --include-bed d4z4_methylation_targets.bed \
  --filter-threshold C:0.75
```

Target regions in `d4z4_methylation_targets.bed`:
```
chr4    193540172    193543634    distal_D4Z4_DUX4_chr4
chr10   134738000    134741148    distal_D4Z4_chr10
```

The target region chr4:193,540,172–193,543,634 covers the most distal D4Z4 repeat unit including the DUX4 promoter and coding region. Key sub-regions:
- **DR1 region** (29 CpGs): Focal demethylation domain; most informative for FSHD2
- **DUX4-PAS region** (10 CpGs): Includes CpG6, the single most discriminating site

**4b. Reference methylation thresholds:**

| Metric | Healthy | FSHD1 (contracted allele) | FSHD2 (all arrays) | Source |
|---|---|---|---|---|
| DR1 region mean CpG methylation | >35% | Variable (hypomethylated) | <25% | Hartweck 2013, Lemmers 2015 |
| DUX4-PAS CpG6 | >78% | ≤78% | ≤78% | Cascella 2018 |
| DUX4 upstream (ONT-specific) | >46.85% | ≤46.85% | ≤46.85% | Huang et al. 2024 (AUC 0.996) |
| DUX4 gene body (ONT-specific) | >58.56% | ≤58.56% | ≤58.56% | Huang et al. 2024 |

---

## Output

Per sample, the pipeline produces:

| Field | Description |
|---|---|
| Sample ID | Sample prefix |
| Chr4 reads: RU counts | D4Z4 repeat units per read (with completeness flag) |
| Chr10 reads: RU counts | D4Z4 repeat units per read (with completeness flag) |
| Chr4 haplotype | 4qA / 4qB per read (from pLAM check) |
| D4Z4 methylation % (chr4) | Mean CpG methylation at chr4:193,540,172–193,543,634 |
| D4Z4 methylation % (chr10) | Mean CpG methylation at chr10:134,738,000–134,741,148 |
| Per-CpG methylation | bedMethyl file with per-site methylation frequencies |
| Complete spanning reads | Count per chromosome |

---

## Query Sequences (`querySeq/`)

| File | Purpose |
|---|---|
| `DBET.fasta` | Head (proximal) marker for D4Z4 array |
| `exon2.fasta` | Tail (distal) marker for D4Z4 array |
| `MIR8078.fasta` | D4Z4 repeat unit marker (83 bp; v1/v2 only, replaced by full 3.3 kb unit in v3) |
| `pLAM-polyA.fasta` | hapA (ATTAAA) / hapB (ATCAAA) for permissive check |
| `pLAM.fasta` | pLAM distal sequence |
| `P13E11.fasta` | p13E-11 probe region for 4q/10q disambiguation |
| `meth_reg.fasta` | Methylation region reference |

**New feature sequences needed for v3** (from DUCKS4 BLAST database):
- D4F104S1 variants: 4qA, 4qB, 10q26
- Full D4Z4 units: chr4 (3,298 bp), chr10 (3,309 bp) from T2T-CHM13v2.0
- pLAM variants: 4qA, 4qB, 10q26

## Test Sample

**ONT019** — confirmed FSHD patient.

| Field | Value |
|---|---|
| Sample ID | WGS-ONT019-5657831-D23881 |
| Input CRAM | `/gene_lv/NGS_gene_lv_Storage/WGS/analysis/pipeline_output/WGS-ONT/WGS-ONT019-5657831-D23881/WGS-ONT019-5657831-D23881_sup_Q13_20250310120409/OUTPUT/1.ALIGNMENT/WGS-ONT019-5657831-D23881.cram` |
| Output dir | `/lab_nas_home_user/Mark/ONT/FSHD/WGS-ONT019-5657831-D23881/` |
| Reference | hg38 (`GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set_modified`) |
| Basecaller | Dorado 0.7.4, `dna_r10.4.1_e8.2_400bps_sup@v4.3.0` |
| Modifications | 5mCG + 5hmCG (`@v1` model) |
| Alignment preset | lr:hq (minimap2 index `.lrhq.mmi`) |
| CRAM format | Reference-free (`no_ref=1`) — readable without original FASTA |
| Size | ~46 GB |

---

## Dependencies

### Conda container (`nf-fshd/Dockerfile`)

| Tool | Channel | Version | Purpose |
|---|---|---|---|
| samtools | bioconda | 1.23 | BAM extraction, reset, sort, index |
| bcftools | bioconda | 1.23 | General BAM/VCF utilities |
| minimap2 | bioconda | 2.30 | T2T alignment with `-y` methylation tag propagation |
| blast | bioconda | 2.17.0 | 4q/10q disambiguation + D4Z4 RU counting |
| ont-modkit | bioconda | 0.6.1 | ONT CpG methylation pileup |
| pysam | bioconda | 0.23.3 | BAM manipulation in Python |
| pandas | conda-forge | 3.0.1 | Data analysis |
| plotnine | conda-forge | 0.15.3 | Visualization |

## Reference Pipelines

- **DUCKS4** (Löwenstern et al. 2026, Human Genomics): Docker, T2T, BLAST + modkit. github.com/tamara-nano/ducks4
- **D4Z4End2End** (Xiao et al. 2025, Genome Research): Full array assembly + Racon consensus. github.com/lucindaxiao/D4Z4End2End
- **Huang et al. (2024, J Transl Med)**: WGS-only ONT at ~29× with 100% OGM concordance, AUC 0.996 for methylation

## Known Limitations (WGS-only)

- **Coverage**: ~5–15× effective at D4Z4 vs 50–500× with Cas9 enrichment. Mitigate with ≥60× mean genome coverage + ultra-long read protocol (SQK-ULK114)
- **Large normal arrays** (≥11 RU = ≥36 kb): Most WGS reads cannot span; pipeline reports minimum repeat count
- **Hybrid/translocated arrays** (~10–20% of individuals): Mixed BlnI/XapI patterns flagged for manual review
