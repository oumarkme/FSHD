# FSHD D4Z4 Repeat Analysis Pipeline (v3)

ONT long-read WGS pipeline for FSHD D4Z4 analysis using R10.4.1 chemistry. Extracts FSHD-relevant reads from hg38 BAM/CRAM, realigns to T2T-CHM13v2.0, counts D4Z4 repeat units per read, classifies haplotype (4qA/4qB), and reports methylation percentage at the distal D4Z4/DUX4 region. FSHD2 modifier gene screening (SMCHD1/DNMT3B/LRIF1) is performed separately outside this pipeline. Full design rationale in `doc/fshd_info.md`.

## Usage

```bash
bash nf-fshd/run.sh <input.bam|cram> [extra nextflow args...]
```

The `run.sh` script auto-detects the environment (WSL or GS4), builds the Apptainer container if needed, and launches Nextflow with `-resume`.

## Server Environment

Both machines share the same mounted filesystem. Development on WSL, production runs on GS4.

| | WSL (dev) | GS4 (production) |
|---|---|---|
| Access | local | `ssh gs4` |
| CPUs | 8 (limit: 5) | 96 (limit: 67) |
| RAM | 30 GB (limit: 21 GB) | 1.7 TB (limit: 1.2 TB) |
| GPU | — | 4x NVIDIA RTX 4090 |
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

### Output Directory Structure

```
/lab_nas_home_user/Mark/ONT/FSHD/                  # Pipeline output root
└── {SAMPLE_ID}/                                    # e.g. WGS-ONT019-5657831-D23881/
    ├── {SAMPLE_ID}.d4z4.bam                        # D4Z4 subset BAM aligned to T2T
    ├── {SAMPLE_ID}.d4z4.bam.bai                    # BAM index
    ├── {SAMPLE_ID}.repeats.tsv                     # Table 1: RU counts by haplotype
    └── {SAMPLE_ID}.methylation.tsv                 # Table 2: Methylation vs controls
```

Nextflow work directory: `/raid/nf/tmp` (local SSD, cleaned after successful runs).

---

## Pipeline Overview

```
EXTRACT_FSHD_READS ─→ REALIGN_T2T ─→ EXTRACT_D4Z4 ─→ BLAST_ANALYSIS ─→ METHYLATION ─→ REPORT
                                            │              │                  │              │
                                      d4z4.bam (pub)  d4z4_counts.tsv   methylation.bed  repeats.tsv (pub)
                                      d4z4.bam.bai                                        methylation.tsv (pub)
```

`(pub)` = published to output directory. Intermediate files stay in the Nextflow work directory.

---

## Pipeline Steps

### Step 1: EXTRACT_FSHD_READS — Targeted read extraction

Extract FSHD-relevant reads from hg38-aligned BAM/CRAM. Includes reads mapping to D4Z4 regions on chr4/chr10 and unmapped reads (the hg38 gap at 4q35 causes many D4Z4 reads to be unmapped).

**hg38 extraction regions:**

| Region | hg38 Coordinates | Purpose |
|---|---|---|
| 4q35 D4Z4 + flanking | chr4:189,900,000–190,215,000 | D4Z4 array + flanking (includes ~82 kb assembly gap) |
| 10q26 D4Z4 + flanking | chr10:133,500,000–133,798,000 | 10q D4Z4 homolog (poorly assembled in hg38) |

```bash
samtools view -b -h -M -L regions.bed input.cram -o target.bam
samtools view -b -f 4 input.cram -o unmapped.bam
samtools merge -o subset.bam target.bam unmapped.bam
```

### Step 2: REALIGN_T2T — T2T-CHM13v2.0 realignment

Convert extracted reads to FASTQ (preserving MM/ML methylation tags), realign to T2T-CHM13v2.0 where the D4Z4 arrays are fully resolved.

```bash
samtools fastq -T MM,ML subset.bam \
  | minimap2 -ax lr:hq -y --eqx --secondary=no -t ${CPUS} chm13v2.0.fa - \
  | samtools sort -o sample.t2t.bam
samtools index sample.t2t.bam
```

| Flag | Purpose |
|---|---|
| `samtools fastq -T MM,ML` | Convert BAM to FASTQ carrying methylation tags in the comment |
| `-ax lr:hq` | ONT R10.4.1 Q20+ preset |
| `-y` | Propagates MM/ML methylation tags from FASTQ comment to SAM output |
| `--eqx` | Emits `=`/`X` CIGAR operators for exact match/mismatch |
| `--secondary=no` | Suppresses secondary alignments (recommended for methylation analysis) |

### Step 3: EXTRACT_D4Z4 — D4Z4 candidate extraction

Extract reads mapping to D4Z4 array regions in T2T coordinates. This subset BAM is published for IGV review.

**T2T-CHM13v2.0 extraction regions:**

| Region | Coordinates | Notes |
|---|---|---|
| Chr4 D4Z4 + flanking | chr4:193,230,000–193,610,000 | ~380 kb; includes full D4Z4 array |
| Chr10 D4Z4 + flanking | chr10:134,560,000–134,800,000 | ~240 kb |

### Step 4: BLAST_ANALYSIS — RU counting + haplotype classification

Build a BLAST database from D4Z4 candidate reads and query with three feature sequences:

| Query file | Purpose |
|---|---|
| `d4z4.fa` | D4Z4 repeat unit (~3.3 kb) for RU counting |
| `pLAM-polyA.fa` | hapA (ATTAAA) / hapB (ATCAAA) for permissive haplotype check + distal anchor |
| `P13E11.fasta` | D4F104S1 (p13E-11) probe for proximal anchor |

**4a. D4Z4 repeat unit counting:**

BLAST the 3.3 kb D4Z4 unit against each read. Hits with alignment length >1,500 bp are counted as valid repeat units (DUCKS4 method).

**4b. Read completeness classification (DUCKS4):**

| Category | Definition | Haplotype determination |
|---|---|---|
| **Complete** | Contains both P13E11 (proximal) AND pLAM (distal) | Yes — from pLAM SNP |
| **Partial-distal** | Contains pLAM but not P13E11 | Yes — from pLAM SNP |
| **Partial-proximal** | Contains P13E11 but not pLAM | No — haplotype = `-` |
| **Incomplete** | D4Z4 only, no flanking anchors | No — haplotype = `-` |

**4c. Permissive haplotype check (chr4 only):**

For chr4-assigned reads containing pLAM, compare BLAST mismatch counts of hapA vs hapB queries:

| PAS variant | Haplotype | FSHD relevance |
|---|---|---|
| ATTAAA (hapA matches better) | 4qA | Permissive — DUX4 mRNA is polyadenylated and stable |
| ATCAAA (hapB matches better) | 4qB | Non-permissive — no stable DUX4 transcript |

Permissive check is only applied to chr4 reads. Chr10 reads matching hapB confirms 10q origin.

### Step 5: METHYLATION — Allele-specific CpG methylation

Runs after BLAST_ANALYSIS. Selects reads for methylation pileup using a tiered priority:

1. **Two haplotypes on chr4** (4qA + 4qB identified): use reads from the shorter (contracted) allele
2. **Complete 4qA reads available**: use those
3. **Any 4qA reads**: use all 4qA reads
4. **Fallback**: all chr4 reads from the D4Z4 BAM

For each tier, filters the D4Z4 BAM to selected reads and runs modkit. If the tier produces no CpG data (reads don't overlap the methylation target), the next tier is tried.

```bash
modkit pileup selected.bam methylation.bed \
  --ref chm13v2.0.fa \
  --cpg \
  --combine-strands \
  --combine-mods \
  --modified-bases 5mC \
  --include-bed meth_targets.bed
```

**Methylation target:** chr4:193,540,172–193,543,634 (distal D4Z4 + DUX4 promoter/coding region).

### Step 6: REPORT — Summary tables

Generates two output tables from the BLAST counts and methylation data.

**Table 1: `{SAMPLE_ID}.repeats.tsv`** — D4Z4 repeat analysis, grouped by chromosome, haplotype, and completeness:

| Chromosome | Haplotype | RU_Counts | Read_Count | Complete |
|---|---|---|---|---|
| chr4 | 4qA | 4 | 1 | TRUE |
| chr4 | - | 1~22 | 51 | FALSE |
| chr10 | - | 9~13 | 4 | FALSE |

- **Haplotype**: `4qA` / `4qB` for reads with pLAM (distal anchor); `-` otherwise
- **RU_Counts**: range as `min~max`, or single value if all reads agree
- **Complete**: `TRUE` only for reads with both P13E11 (proximal) and pLAM (distal) anchors

**Table 2: `{SAMPLE_ID}.methylation.tsv`** — Methylation at two sub-regions vs healthy controls:

| Region | Methylation_Pct | Control_Methylation_IQR [25%, 50%, 75%] |
|---|---|---|
| DR1 | 22.58 | [49.91, 52.29, 54.65] |
| DUX4-PAS | 30.0 | [53.13, 55.33, 59.2] |

**Methylation sub-regions (T2T-CHM13v2.0 chr4 coordinates):**

| Region | Coordinates | CpGs | Clinical significance |
|---|---|---|---|
| DR1 | chr4:193,541,000–193,541,600 | 29 | D4Z4 promoter; focal demethylation domain. <25% suggests FSHD2 |
| DUX4-PAS | chr4:193,543,430–193,543,634 | 10 | Polyadenylation signal region; includes CpG6 |

**Reference methylation thresholds:**

| Metric | Healthy | FSHD1 (contracted allele) | FSHD2 (all arrays) | Source |
|---|---|---|---|---|
| DR1 region mean CpG methylation | >35% | Variable (hypomethylated) | <25% | Hartweck 2013, Lemmers 2015 |
| DUX4-PAS CpG6 | >78% | <=78% | <=78% | Cascella 2018 |
| DUX4 upstream (ONT-specific) | >46.85% | <=46.85% | <=46.85% | Huang et al. 2024 (AUC 0.996) |
| DUX4 gene body (ONT-specific) | >58.56% | <=58.56% | <=58.56% | Huang et al. 2024 |

---

## Reference Genomes

| Genome | Path |
|---|---|
| hg38 | `/lab_nas_home_user/Mark/reference_genome/hg38/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set_modified.fasta` |
| T2T-CHM13v2.0 | `/lab_nas_home_user/Mark/reference_genome/t2t/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna` |

Both use chr-prefixed naming (chr1, chr2, ...). The T2T reference is the NCBI version with standard chromosome names.

**Key T2T-CHM13v2.0 coordinates:**

| Region | Coordinates | Notes |
|---|---|---|
| Chr4 D4Z4 array | chr4:193,282,466–193,558,261 | ~276 kb; 33 RUs in reference |
| Chr4 DUX4 methylation target | chr4:193,540,172–193,543,634 | Distal D4Z4 + DUX4 promoter/coding |
| Chr4 p13E-11 anchor | chr4:192,667,301–192,902,247 | Upstream unique anchoring region |
| Chr10 D4Z4 array | chr10:134,615,120–134,741,148 | ~126 kb |

---

## Query Sequences (`nf-fshd/query/`)

| File | Used | Purpose |
|---|---|---|
| `d4z4.fa` | Yes | D4Z4 repeat unit for RU counting (hits >1,500 bp = 1 RU) |
| `pLAM-polyA.fa` | Yes | hapA (ATTAAA) / hapB (ATCAAA) for permissive check + distal anchor |
| `P13E11.fasta` | Yes | D4F104S1 (p13E-11) probe for proximal anchor |
| `dbet.fa` | No | Previously used as proximal marker (replaced by P13E11) |
| `exon2.fasta` | No | Previously used as distal marker (replaced by pLAM) |

---

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

## Container

All Nextflow processes use a local Apptainer SIF image built from `nf-fshd/environment.def`. The `run.sh` script automatically rebuilds the SIF when the def file is newer or the SIF doesn't exist.

### Container packages (`nf-fshd/environment.def`)

| Tool | Channel | Version | Purpose |
|---|---|---|---|
| samtools | bioconda | 1.23 | BAM/CRAM extraction, FASTQ conversion, sort, index, merge |
| bcftools | bioconda | 1.23 | General BAM/VCF utilities |
| minimap2 | bioconda | 2.30 | T2T alignment with `-y` methylation tag propagation |
| blast | bioconda | 2.17.0 | D4Z4 RU counting + haplotype classification |
| ont-modkit | bioconda | 0.6.1 | ONT CpG methylation pileup (`--modified-bases 5mC`) |
| pysam | bioconda | 0.23.3 | BAM read iteration in Python |
| pandas | conda-forge | 3.0.1 | Data analysis for BLAST results and report generation |
| plotnine | conda-forge | 0.15.3 | Visualization |

---

## Reference Pipelines

- **DUCKS4** (Lowenstern et al. 2026, Human Genomics): Docker, T2T, BLAST + modkit. github.com/tamara-nano/ducks4
- **D4Z4End2End** (Xiao et al. 2025, Genome Research): Full array assembly + Racon consensus. github.com/lucindaxiao/D4Z4End2End
- **Huang et al. (2024, J Transl Med)**: WGS-only ONT at ~29x with 100% OGM concordance, AUC 0.996 for methylation

## Known Limitations (WGS-only)

- **Coverage**: ~5-15x effective at D4Z4 vs 50-500x with Cas9 enrichment. Mitigate with >=60x mean genome coverage + ultra-long read protocol (SQK-ULK114)
- **Large normal arrays** (>=11 RU = >=36 kb): Most WGS reads cannot span; pipeline reports minimum repeat count
- **Hybrid/translocated arrays** (~10-20% of individuals): Mixed BlnI/XapI patterns flagged for manual review
- **Methylation allele specificity**: When haplotype-classified reads don't overlap the methylation target region, the pipeline falls back to all chr4 D4Z4 reads, which may mix alleles
