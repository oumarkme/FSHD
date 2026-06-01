# FSHDetector: ONT D4Z4 Repeat Unit Counting, Haplotyping, and Methylation Pipeline

A Nextflow DSL2 pipeline for FSHD diagnosis from hg38-aligned ONT BAM/CRAM files with HP tags from upstream phasing. The pipeline resolves all four D4Z4 alleles (2× chr4, 2× chr10), determines 4qA-L/4qA-S/4qB haplotype, counts D4Z4 repeat units per allele via per-read BLAST annotation, validates results against biological constraints, and reports allele-specific CpG methylation at the distal DUX4 region. All analysis coordinates reference **T2T-CHM13v2.0**.

> **Just want to run it?** → **[Quick start](#quick-start)**. You need a phased, hg38-aligned ONT BAM/CRAM plus [Nextflow and Conda](#prerequisites) — no other installs.

## Contents

**Getting started**
- [Usage](#usage) — install, references, and how to run · [Prerequisites](#prerequisites) · [Reference genomes](#reference-genomes) · [Quick start](#quick-start) · [Profiles](#profiles) · [Input](#input)
- [Upstream prerequisites](#upstream-prerequisites) — how to produce the input BAM/CRAM
- [Output](#output) — files, [`repeats.tsv`](#repeatstsv), [`methylation.tsv`](#methylationtsv)

**How it works**
- [Pipeline flow diagram](#pipeline-flow-diagram) · [Pipeline steps](#pipeline-steps)
- [Algorithm](#algorithm-blast-based-per-read-annotation) — per-read BLAST annotation, chromosome & haplotype logic
- [Biological background](#biological-background) — D4Z4 & FSHD, haplotype backgrounds, array sizes
- [Technical constraints (ONT)](#technical-constraints-ont-sequencing) — read length, coverage, HP tags, chr4/chr10
- [Reference genome & query sequences](#reference-genome-and-query-sequences)
- [Software](#software) · [Known limitations](#known-limitations) · [References](#references)

---

## Usage

### Prerequisites

Only two host tools are required:

| Tool | Purpose | Install |
|------|---------|---------|
| **Nextflow** (≥23.04) | Workflow engine | https://www.nextflow.io/docs/latest/install.html |
| **Conda** (Miniforge, ships `mamba`) | Builds the analysis environment automatically | https://github.com/conda-forge/miniforge |

Nothing else needs to be installed: on the first run Nextflow builds the tool environment (samtools, minimap2, BLAST+, modkit, pysam, pandas) from `environment.yml` and caches it for later runs.

### Reference genomes

Two FASTA references are needed (place the T2T `.fai` index next to its FASTA):

| Parameter | Reference | Notes |
|-----------|-----------|-------|
| `--ref_hg38` | The hg38 build your input was aligned to | Decodes CRAM and is the coordinate basis for region extraction. **Required for CRAM input; optional for BAM.** |
| `--ref_t2t` | T2T-CHM13v2.0 (`GCF_009914755.1_T2T-CHM13v2.0_genomic.fna`) | All analysis coordinates. Index once: `samtools faidx <t2t.fna>`. |

### Quick start

```bash
git clone https://github.com/oumarkme/FSHD.git
cd FSHD

nextflow run FSHDetector/main.nf \
    --input    /path/to/sample.bam \
    --ref_hg38 /path/to/hg38.fa \
    --ref_t2t  /path/to/T2T-CHM13v2.0.fna \
    --outdir   results
```

Conda is the default engine, so the first run builds the environment (a few minutes) and reuses it afterwards. Per-sample results are written under `--outdir` (default `results/`).

A convenience wrapper is also provided — it just adds `-resume`:

```bash
bash FSHDetector/run.sh /path/to/sample.bam \
    --ref_hg38 /path/to/hg38.fa --ref_t2t /path/to/T2T-CHM13v2.0.fna
```

### Profiles

| Profile | Engine | When to use |
|---------|--------|-------------|
| *(default)* or `-profile conda` | Conda (from `environment.yml`) | Normal runs — recommended for everyone |
| `-profile test` | none (stub) | Wiring smoke-test with no data (see `nextflow.config`) |

### Input

hg38-aligned BAM/CRAM with **HP/PS tags** (from any phasing tool) and **MM/ML tags** (from basecalling). The matching `.bai`/`.crai` index must sit next to it. BAM input does not require `--ref_hg38`; CRAM input does. See [Upstream prerequisites](#upstream-prerequisites) for how to produce this file.

---

## Upstream prerequisites

These steps are performed outside this pipeline. The pipeline input is an hg38-aligned BAM/CRAM carrying both methylation tags (MM/ML) from basecalling and haplotype tags (HP/PS) from phasing. Any phasing tool that writes standard HP/PS tags is compatible.

| Step | Tool | Version | Output | Reference |
|------|------|---------|--------|-----------|
| Basecalling | Dorado (SUP + `--modified-bases 5mCG_5hmCG`) | ≥0.7.0 | Unaligned modBAM with MM/ML tags | Oxford Nanopore Technologies (2024b) |
| Genome-wide alignment | minimap2 (`-ax lr:hq`) | ≥2.28 | hg38-aligned BAM | Li (2018) |
| Variant calling | e.g., Clair3 | — | VCF (SNPs for phasing) | Zheng et al. (2022) |
| Phasing + haplotagging | e.g., LongPhase, WhatsHap | — | BAM with HP/PS tags | Lin et al. (2022); Martin et al. (2023) |

---

## Output

Seven published files per sample:

| File | Description |
|------|-------------|
| `{SAMPLE_ID}.d4z4.bam` | D4Z4 candidate reads aligned to T2T-CHM13v2.0 |
| `{SAMPLE_ID}.d4z4.bam.bai` | BAM index |
| `{SAMPLE_ID}.complete_reads.bam` | Complete reads (Has_Head AND Has_Tail) extracted from D4Z4 BAM |
| `{SAMPLE_ID}.complete_reads.bam.bai` | BAM index |
| `{SAMPLE_ID}.repeats.tsv` | D4Z4 repeat unit counts grouped by haplotype |
| `{SAMPLE_ID}.methylation.tsv` | CpG methylation at distal D4Z4/DUX4 region |
| `{SAMPLE_ID}.validation.txt` | Biological validation report with complete read details |

### repeats.tsv

One row per (Chromosome, Haplotype, Head, Tail) group. Chr4 rows appear above chr10 rows.

| Column | Description |
|--------|-------------|
| Chromosome | chr4 or chr10 |
| Haplotype | 4qA-S / 4qA-L / 4qA / 4qB / 4q / 10qA / 10qB / 10q (prefix matches chromosome) |
| RU_Count | Maximum RU count among reads in this group. Prefixed with **≥** for incomplete reads (missing Head or Tail), indicating the count is a lower bound |
| Head | TRUE if reads in this group carry a proximal marker (D4F104S1 or p13E-11) |
| Tail | TRUE if reads in this group carry a distal marker (4qA/4qB/10qA pLAM) |
| Read_Count | Number of reads in this group |

### methylation.tsv

| Column | Description |
|--------|-------------|
| Haplotype | HP1, HP2, or Combined |
| Region | DR1 or DUX4-PAS |
| Methylation | Mean % modified across CpGs in region |
| N_CpG | Number of CpGs with data |
| Ref_Atlas (Mean (SD) %) | ONT methylation reference atlas (\~45–49 healthy controls; Huang et al., 2024) |
| Ref_Atlas (IQR %) | Interquartile range [Q1, Q2, Q3] computed from all CpG REF_METH values in the reference atlas |
| Lit_Caputo (%) | Diagnostic threshold (Caputo et al., 2022; Strafella et al., 2024) |

**Reference values (T2T coordinates):**

| Region | T2T coordinates | Ref_Atlas Mean (SD) | Ref_Atlas IQR [Q1, Q2, Q3] | Lit_Caputo threshold |
|--------|----------------|---------------------|---------------------------|---------------------|
| DR1 | chr4:193,541,000–193,541,600 | 54.70% (1.97%) | [53.07%, 54.89%, 56.83%] | CpG1 ≤ 53% |
| DUX4-PAS | chr4:193,543,430–193,543,634 | 51.81% (0.23%) | [51.63%, 51.75%, 51.99%] | CpG6 ≤ 78% |

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
(RU count, anchors, BlnI+XapI chr,               |
 haplotype, validation)                          |
    |                                            |
    |                     Step 5: METHYLATION <--'
    |                     (modkit --phased from T2T BAM)
    |                            |
    '----------+-----------------'
               |
               v
           Step 6: REPORT
           (repeats.tsv + methylation.tsv)
```

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

### Step 4: Chromosome assignment, feature annotation, and validation (BLAST_ANNOTATE)

BLAST-based per-read annotation following the approach of Lowenstern et al. (2026), Huang et al. (2024), and Butterfield et al. (2023). Three BLAST queries annotate each read, followed by biological validation:

```
blastn d4z4.fa              -> RU counting (hits > 1500 bp)
blastn haplotype_markers.fa -> 6 haplotype-specific markers (proximal + distal)
blastn P13E11.fasta         -> additional proximal anchor
```

After classification, the pipeline validates results against biological constraints and flags impossible configurations in a `Flag` column and a `validation.txt` report:

| Check | Rule | Flag |
|-------|------|------|
| Haplotype count | At most 2 typed haplotypes per chromosome (excluding 4q/10q fragments) | WARN in validation.txt |
| Completeness consistency | An incomplete read longer than the shortest complete read in the same haplotype group | `longer_than_complete` |
| RU consistency | An incomplete read with more RUs than the max complete read RU in the same haplotype group | `RU_exceeds_complete` |
| Fragment RU sanity | A fragment read (4q/10q) with RU exceeding the max complete RU across all alleles on that chromosome | `RU_exceeds_all_alleles` |
| Chimeric detection | Proximal marker type ≠ distal marker type (possible translocation or chimeric read) | `chimeric` |

### Step 5: Per-allele CpG methylation (METHYLATION)

```
modkit pileup <t2t.bam> meth_out/ \
    --ref <t2t_reference.fna> \
    --cpg --combine-strands \
    --modified-bases 5mC \
    --phased \
    --include-bed meth_targets.bed
```

Target regions in `meth_targets.bed`: chr4:193,540,172–193,543,634 and chr10:134,738,000–134,741,148. Runs independently; only needs the T2T BAM from Step 2.

### Step 6: Generate report (REPORT)

Produces `repeats.tsv` and `methylation.tsv` from BLAST annotations and methylation pileup data.

---

## Algorithm: BLAST-based per-read annotation

Each read is annotated independently after BLAST alignment against three query files (D4Z4 repeat unit, 6 haplotype-specific markers, P13E11). No consensus assembly or allele separation is performed; the maximum RU count observed across reads within a (Chromosome, Haplotype, Head, Tail) group is reported.

### Chromosome and haplotype determination flow

```
Read with D4Z4 BLAST hits (RU ≥ 1, alignment > 1500 bp)
│
├─ EVIDENCE GATHERING (parallel)
│   ├─ BAM alignment → initial chr (Tier 1, weakest)
│   ├─ BlnI motif count: CTAGGTCTAGG (chr4) vs CTAGGCCTAGG (chr10)
│   ├─ XapI motif count: [AG]AATT[CT] (chr4-specific)
│   ├─ Best proximal marker: 4qA / 4qB / 10qA D4F104S1 (by bitscore)
│   ├─ Best distal marker:   4qA / 4qB / 10qA pLAM (by bitscore)
│   └─ P13E11 anchor (proximal completeness)
│
├─ CHROMOSOME ASSIGNMENT (three tiers, strongest last)
│   │
│   ├─ Tier 1: BAM alignment to T2T → chr4 or chr10 (initial)
│   │
│   ├─ Tier 2: BlnI + XapI motif voting (overrides Tier 1)
│   │   chr4_score  = BlnI_chr4 + XapI
│   │   chr10_score = BlnI_chr10
│   │   ├─ chr10_score > chr4_score → chr10
│   │   ├─ chr4_score > chr10_score → chr4
│   │   └─ tie → retain Tier 1
│   │
│   └─ Tier 3: Haplotype marker concordance (strongest)
│       ├─ proximal AND distal agree on chr → OVERRIDE to that chr
│       ├─ proximal AND distal disagree → keep Tier 1/2 (chimeric)
│       └─ single marker + motif total = 0 → use marker as tiebreaker
│
├─ HAPLOTYPE ASSIGNMENT (prefix always matches chromosome)
│   │
│   ├─ Has distal marker?
│   │   ├─ 4qA_pLAM or 10qA_pLAM → A-type
│   │   │   ├─ chr = chr4 + 4qA_pLAM → subtype by terminal D4Z4 unit:
│   │   │   │   ├─ closest D4Z4 to pLAM > 2800 bp → 4qA-L
│   │   │   │   ├─ closest D4Z4 to pLAM ≤ 2800 bp → 4qA-S
│   │   │   │   └─ no D4Z4 position data → 4qA
│   │   │   ├─ chr = chr4 + 10qA_pLAM → 4qA (translocation)
│   │   │   ├─ chr = chr10 + any A-type pLAM → 10qA
│   │   │   └─ unknown chr → Unknown
│   │   │
│   │   └─ 4qB_pLAM → B-type
│   │       ├─ chr = chr4 → 4qB
│   │       └─ chr = chr10 → 10qB
│   │
│   ├─ No distal, has proximal D4F104S1?
│   │   ├─ *A_D4F104S1 → *qA (prefix from chromosome)
│   │   └─ *B_D4F104S1 → *qB (prefix from chromosome)
│   │
│   └─ No markers at all?
│       ├─ chr4 → 4q
│       └─ chr10 → 10q
│
├─ CHIMERIC CHECK
│   └─ proximal marker type ≠ distal marker type → Chimeric = TRUE
│
├─ PERMISSIVE CHECK (sequence-dependent, not chromosome-dependent)
│   └─ Distal_Marker = 4qA_pLAM AND Chimeric = FALSE → Permissive = TRUE
│       (functional ATTAAA PAS; applies to chr4 or chr10 translocations)
│       All other cases → Permissive = FALSE
│
├─ COMPLETENESS
│   ├─ Has_Head = proximal marker hit (D4F104S1 or P13E11)
│   └─ Has_Tail = distal marker hit (4qA/4qB/10qA pLAM)
│
└─ OUTPUT: reads_annotated.tsv
    Chromosome, Haplotype, Permissive, RU_Count, Has_Head, Has_Tail,
    Proximal_Marker, Distal_Marker, Chimeric, Chr_Confidence, ...
```

### Chromosome confidence scoring

| Level | Criterion |
|-------|-----------|
| HIGH | All motifs agree (minority side score = 0), or majority ≥ 80% of total motifs |
| MEDIUM | Majority ≥ 60% but < 80% of total motifs |
| LOW | No motif evidence (BAM-only), or majority < 60% (mixed evidence) |

### 4qA-S / 4qA-L subtyping

4qA-S (standard) and 4qA-L (long) differ only in the length of the most distal (terminal) D4Z4 repeat unit adjacent to pLAM (Lemmers et al., 2018). For reads with a `4qA_pLAM` hit, the pipeline finds the D4Z4 BLAST hit closest to the pLAM position: **>2800 bp** → **4qA-L** (full-length \~3.3 kb); **≤2800 bp** → **4qA-S** (truncated \~1.8 kb). Both subtypes are permissive.

---

## Biological background

### D4Z4 and FSHD

The D4Z4 macrosatellite array on chromosome 4q35 consists of tandem 3.3 kb repeat units. Healthy individuals carry 11–100+ RU per allele. FSHD1 occurs when one allele contracts to 1–10 RU **on a permissive haplotype background** (4qA or 4qA-L), which stabilizes DUX4 mRNA via a polyadenylation signal (PAS) in the pLAM region distal to the array (Lemmers et al., 2012).

### Haplotype backgrounds

D4Z4 haplotypes are defined by three features: a proximal SSLP (simple sequence length polymorphism) allele, the D4F104S1 (p13E-11) SNP pattern (\~15 diagnostic SNPs), and the distal sequence structure. At least 17 distinct haplotypes have been identified across chr4 and chr10 (Lemmers et al., 2007, 2010; van der Vliet et al., 2015).

| Haplotype | Chr | SSLP | Distal | β-satellite | PAS | D4F104S1 | XapI | BlnI | Permissive | Freq | Reference |
|-----------|-----|------|--------|-------------|-----|----------|------|------|------------|------|-----------|
| **4A161** (4qA-S) | chr4 | 161 | pLAM (241 bp) | 6.2 kb | ATTAAA | chr4-type | + | − | **Yes** | Most common 4qA | Lemmers et al. (2002) |
| **4A161L** (4qA-L) | chr4 | 161 | pLAM (241 bp) | 6.2 kb | ATTAAA | chr4-type | + | − | **Yes** | \~20% of 4A161 | Lemmers et al. (2018) |
| **4A159** | chr4 | 159 | pLAM | present | ATTAAA | chr4-type | + | − | **Yes** | Uncommon | van der Vliet et al. (2015) |
| **4A168** | chr4 | 168 | pLAM | present | ATTAAA | chr4-type | + | − | **Yes** | Uncommon | van der Vliet et al. (2015) |
| **4A166** | chr4 | 166 | pLAM (241 bp) | 6.2 kb | ATTAAA | **chr10-like** | **−** | mixed | No | \~11% of 4qA | Lemmers et al. (2007, 2018) |
| **4B163** | chr4 | 163 | 4qB distal (181 bp) | absent | absent | chr4-type | + | − | No | \~50% of chr4 | Lemmers et al. (2002) |
| **4B168** | chr4 | 168 | 4qB distal | absent | absent | chr4-type | + | − | No | Uncommon | Lemmers et al. (2010) |
| **10A166** | chr10 | 166 | pLAM-like (241 bp) | 6.2 kb | ATCAAA | chr10-type | − | + | No | Most common 10q | van der Vliet et al. (2015) |
| **10A164** | chr10 | 164 | pLAM-like | present | ATCAAA | chr10-type | − | + | No | Common | Lemmers et al. (2010) |
| **10A176T** | chr10 | 176 | pLAM-like | present | ATTTAA | chr10-type | − | + | No | Rare | van der Vliet et al. (2015) |
| **10B161T** | chr10 | 161 | B-type distal | absent | absent | chr10-type | − | + | No | Rare | Lemmers et al. (2010) |

Column definitions:
- **SSLP**: Length (bp) of the simple sequence length polymorphism proximal to D4Z4.
- **Distal**: Sequence immediately distal to the last D4Z4 unit. A-type carries pLAM (\~241 bp, GenBank HM190178.1 / HM190187.1); B-type carries a distinct \~181 bp sequence (HM190161.1) unrelated to pLAM.
- **β-satellite**: 6.2 kb tandem repeat of 68-bp units distal to pLAM; present in A-type, absent in B-type.
- **PAS**: Polyadenylation signal in pLAM. ATTAAA is functional (stabilizes DUX4 mRNA); ATCAAA and ATTTAA are non-functional (single SNP disrupts PAS).
- **D4F104S1**: Proximal marker SNP pattern. "chr4-type" matches GU550601.1 (4A161); "chr10-type" matches GU550630.1 (10A164b); "chr10-like" = 4A166 resembles chr10 (Lemmers et al., 2007).
- **XapI**: Restriction site `RAATTY` in proximal D4Z4 unit. Present (+) on chr4-type units; absent (−) on chr10-type and 4A166.
- **BlnI**: Restriction site `CCTAGG` within D4Z4 units. Absent (−) on chr4-type units (`CTAGGTCTAGG`); present (+) on chr10-type units (`CTAGGCCTAGG`).

**4qA-S vs 4qA-L**: Both share the same SSLP (161), D4F104S1 pattern, pLAM sequence, and ATTAAA PAS. They differ only in the length of the most distal (terminal) D4Z4 repeat unit adjacent to pLAM: 4qA-S has a truncated terminal unit (\~1.8 kb, ≤2800 bp) while 4qA-L has a full-length terminal unit (\~3.3 kb, >2800 bp). Both are permissive for FSHD1 (Lemmers et al., 2018).

**4A166**: Despite having a functional ATTAAA PAS, 4A166 is non-permissive — DUX4 polyadenylation is very inefficient on this haplotype. The proximal D4F104S1 SNPs and D4Z4 XapI pattern resemble 10qA166, suggesting an evolutionary hybrid origin. The mechanism of non-permissiveness remains under investigation (Lemmers et al., 2007, 2018). The pipeline will classify 4A166 reads as chr10/10qA with a chimeric flag (proximal=10qA, distal=4qA); see Known limitations #9.

**Subtelomeric exchanges**: The 4q35 and 10q26 regions share \~98.5% identity. Interchromosomal exchanges occur in \~3–10% of the population (Lemmers et al., 2010), producing hybrid alleles (e.g., chr4 with 10qA-type pLAM or vice versa). Such exchanges can revert FSHD pathogenicity (Moriatis et al., 2024). The pipeline detects these as chimeric reads (proximal marker type ≠ distal marker type).

Only contractions on permissive haplotypes (4A161/4A161L/4A159/4A168) cause FSHD1. Contractions on all other haplotypes are clinically silent (Lemmers et al., 2012).

### D4Z4 array sizes

| RU count | Array size | Clinical significance |
|----------|-----------|----------------------|
| 1–3 | 3.3–10 kb | Severe/infantile FSHD1 (high penetrance) |
| 4–7 | 13–23 kb | Classic FSHD1 (variable expressivity) |
| 8–10 | 26–33 kb | Borderline; reduced penetrance (\~40–50%); requires methylation |
| 11–20 | 36–66 kb | Normal; may contribute to FSHD2 if modifier present |
| 21–100+ | 69–330+ kb | Normal/healthy |

---

## Technical constraints (ONT sequencing)

### Read length vs D4Z4 spanning

| Parameter | Standard ligation (LSK114, with size selection) | Ultra-long (ULK114) |
|-----------|------------------------------------------------|---------------------|
| Typical N50 | 25–30 kb | 50–100+ kb |
| Max reliable spanning | \~7–8 RU (\~26 kb) | \~15–20 RU |
| Marginal spanning | 9–10 RU (\~33 kb) | 20–30 RU |
| Never spanned | ≥11 RU (≥36 kb) | ≥30 RU |

A read is **"complete"** when it contains both the proximal p13E-11 anchor AND the distal pLAM anchor. With standard N50 \~25–30 kb, complete reads can only come from alleles with ≤\~10 RU. Healthy alleles (≥11 RU) will never produce complete reads with standard prep.

### Coverage per allele

| Genome coverage | Reads per allele (theoretical) | Effective (subtelomeric penalty) |
|-----------------|-------------------------------|----------------------------------|
| 30× | \~7–8× | \~4–6× |
| 60× | \~15× | \~10–12× |

### HP tag reliability

- HP tags are reliable in unique flanking regions adjacent to D4Z4
- Within the D4Z4 array: very few heterozygous SNPs → phasing may fail
- Most D4Z4 reads carry HP tags inherited from flanking-region phasing, but some reads (especially those entirely within the array) may lack HP tags

### Chr4 vs chr10 discrimination

- D4Z4 arrays share \~98.5–99% identity between chr4 and chr10
- Chromosome discrimination uses three tiers: BAM alignment, BlnI/XapI restriction motif voting, and haplotype marker concordance
- The BlnI motif distinguishes chr4-type (`CTAGGTCTAGG`) from chr10-type (`CTAGGCCTAGG`) D4Z4 units
- The XapI motif (`RAATTY` = `[AG]AATT[CT]`) is chr4-specific (Butterfield et al., 2023); absent in 4A166 and chr10
- Haplotype-specific D4F104S1 (proximal) and pLAM (distal) markers provide independent chromosome and haplotype evidence

---

## Reference genome and query sequences

**Reference:** T2T-CHM13v2.0 (Nurk et al., 2022; `GCF_009914755.1_T2T-CHM13v2.0_genomic.fna`)

**Key T2T-CHM13v2.0 coordinates:**

| Region | Coordinates | Notes |
|--------|------------|-------|
| Chr4 D4Z4 array | chr4:193,282,466–193,558,261 | \~276 kb; 33 RUs in reference |
| Chr4 DUX4 methylation target | chr4:193,540,172–193,543,634 | Distal D4Z4 + DUX4 promoter/coding |
| Chr4 p13E-11 anchor | chr4:192,667,301–192,902,247 | Upstream unique anchoring region |
| Chr10 D4Z4 array | chr10:134,615,120–134,741,148 | \~126 kb |

**Query sequences (`FSHDetector/query/`):**

| File | Sequences | Purpose |
|------|-----------|---------|
| `d4z4.fa` | Single D4Z4 RU (3,297 bp) | RU counting |
| `haplotype_markers.fa` | 6 haplotype-specific markers (see below) | Proximal/distal haplotype discrimination |
| `P13E11.fasta` | p13E-11 probe (220 bp) | Additional proximal anchor |

**Haplotype markers (`haplotype_markers.fa`):**

| Marker | GenBank | Length | Position | Purpose |
|--------|---------|--------|----------|---------|
| `4qA_D4F104S1` | GU550601.1 (4A161) | 489 bp | Proximal | 4qA-specific D4F104S1 |
| `4qB_D4F104S1` | GU550603.1 (4B163) | 489 bp | Proximal | 4qB-specific D4F104S1 |
| `10qA_D4F104S1` | GU550630.1 (10A164b) | 489 bp | Proximal | 10qA-specific D4F104S1 |
| `4qA_pLAM` | HM190178.1:3205-3445 | 241 bp | Distal | 4qA pLAM with ATTAAA PAS |
| `4qB_pLAM` | HM190161.1:3124-3304 | 181 bp | Distal | 4qB distal (no functional PAS) |
| `10qA_pLAM` | HM190187.1:3214-3454 | 241 bp | Distal | 10qA pLAM with ATCAAA PAS |

---

## Software

### Analysis environment

Nextflow provisions the tools below automatically from `environment.yml` (Conda) on the first run, pinned to the versions listed for reproducibility.

| Tool | Version | Purpose | Reference |
|------|---------|---------|-----------|
| Python | 3.12 | Scripting runtime | |
| samtools | 1.23 | BAM/CRAM extraction, sort, index, merge, FASTQ conversion | Danecek et al. (2021) |
| minimap2 | 2.30 | T2T realignment (Step 2 only) | Li (2018) |
| BLAST+ | 2.17.0 | RU counting, anchor detection, haplotyping, chromosome assignment | Camacho et al. (2009) |
| modkit | 0.6.1 | Phased CpG methylation pileup (5mC only) | Oxford Nanopore Technologies (2024a) |
| pysam | 0.23.3 | BAM read iteration, motif counting in Python | |
| pandas | 3.0.1 | Tabular data processing and report generation | |

### Runtime requirements (host)

| Tool | Purpose | Reference |
|------|---------|-----------|
| Nextflow | Workflow engine (DSL2) | Di Tommaso et al. (2017) |
| Conda | Builds the analysis environment from `environment.yml`; Miniforge recommended | |

---

## Known limitations

1. **Read length constraint.** At N50 \~25–30 kb, complete reads are limited to alleles ≤\~10 RU. Ultra-long prep extends this range but is not standard.

2. **Subtelomeric coverage.** Effective per-allele coverage at D4Z4 is 30–50% lower than genome-wide average. At 30× WGS, expect \~4–6 reads per allele.

3. **Per-read BLAST counting.** Per-read BLAST counting reflects only the D4Z4 content within each individual read. For alleles longer than the read length, the maximum observed RU count across all reads is a lower bound.

4. **Methylation coverage.** At 30× WGS, per-allele methylation coverage may be only 3–5 reads. ≥60× recommended for clinical-grade calls.

5. **CHM13 reference gap for 4qA.** T2T-CHM13v2.0 carries 4qB on chr4 (Nurk et al., 2022). 4qA-specific pLAM/PAS is absent from the reference. BLAST-based annotation detects pLAM in read sequence directly.

6. **Somatic mosaicism.** Mosaic samples may show multiple RU count clusters from the same allele. The pipeline does not explicitly detect mosaicism but unusual clustering patterns should be flagged for review.

7. **Cis duplications.** Two adjacent D4Z4 arrays on one allele (\~1.5% of Europeans) may produce two distinct RU groups in one HP. The pipeline may misinterpret this as two alleles.

8. **Biological validation flags.** The validation step detects inconsistencies (e.g., incomplete reads longer than complete reads, RU counts exceeding complete-read counts) but does not auto-correct them. Flagged reads should be reviewed manually. Common causes include chimeric reads, BLAST sensitivity limits, and reads from different alleles with spurious cross-haplotype BLAST hits.

9. **4A166 haplotype misclassification.** The 4A166 haplotype (\~11% of 4qA alleles) has chr10-like features: its D4F104S1 SNP pattern resembles 10qA166, and its proximal D4Z4 unit lacks the chr4-specific XapI site (Lemmers et al., 2007). As a result, 4A166 reads may be assigned to chr10 by both BAM alignment (Tier 1) and motif voting (Tier 2), and the proximal marker will best-match `10qA_D4F104S1`. However, the distal marker will match `4qA_pLAM` (ATTAAA PAS), causing the read to be flagged as **chimeric** (proximal=10qA, distal=4qA). Since 4A166 is non-permissive for FSHD (like 10qA), this misclassification has no clinical impact. A chimeric flag with proximal=10qA and distal=4qA in the output should prompt consideration of the 4A166 haplotype.

---

## References

### FSHD biology, haplotypes, and clinical

Caputo, V., Megalizzi, D., Fabrizio, C., Termine, A., Colantoni, L., Bax, C., Gimenez, J., Monforte, M., Tasca, G., Ricci, E., Caltagirone, C., Giardina, E., Cascella, R., & Strafella, C. (2022). D4Z4 methylation levels combined with a machine learning pipeline highlight single CpG sites as discriminating biomarkers for FSHD patients. *Cells*, *11*(24), Article 4114. https://doi.org/10.3390/cells11244114

Hiramuki, Y., Kure, Y., Saito, Y., Ogawa, M., Ishikawa, A., Murata, M., Hayashi, Y. K., Oya, Y., Komaki, H., Nonaka, I., Nishino, I., Inoue, M., Sasaki, M., & Mitsuhashi, S. (2018). Simultaneous measurement of the size and methylation of chromosome 4qA-D4Z4 repeats in facioscapulohumeral muscular dystrophy by long-read sequencing. *Journal of Translational Medicine*, *20*, Article 517. https://doi.org/10.1186/s12967-022-03743-7

Lemmers, R. J. L. F., de Kievit, P., Sandkuijl, L., Padberg, G. W., van Ommen, G.-J. B., Frants, R. R., & van der Maarel, S. M. (2002). Facioscapulohumeral muscular dystrophy is uniquely associated with one of the two variants of the 4q subtelomere. *Nature Genetics*, *32*(2), 235–236. https://doi.org/10.1038/ng999

Lemmers, R. J. L. F., Wohlgemuth, M., van der Gaag, K. J., van der Vliet, P. J., van Teijlingen, C. M. M., de Knijff, P., Padberg, G. W., Frants, R. R., & van der Maarel, S. M. (2007). Specific sequence variations within the 4q35 region are associated with facioscapulohumeral muscular dystrophy. *American Journal of Human Genetics*, *81*(5), 884–894. https://doi.org/10.1086/521986

Lemmers, R. J. L. F., van der Vliet, P. J., van der Gaag, K. J., Jacobs, S. J., van der Stoep, N., de Knijff, P., & van der Maarel, S. M. (2010). Worldwide population analysis of the 4q and 10q subtelomeres identifies only four discrete interchromosomal sequence transfers in human evolution. *American Journal of Human Genetics*, *86*(3), 364–377. https://doi.org/10.1016/j.ajhg.2010.01.035

Lemmers, R. J. L. F., Tawil, R., Petek, L. M., Balog, J., Block, G. J., Santen, G. W. E., Amell, A. M., van der Vliet, P. J., Almomani, R., Straasheijm, K. R., Krom, Y. D., Klooster, R., Sun, Y., den Dunnen, J. T., Helmer, Q., Donlin-Smith, C. M., Padberg, G. W., van Engelen, B. G. M., de Greef, J. C., ... van der Maarel, S. M. (2012). Digenic inheritance of an SMCHD1 mutation and an FSHD-permissive D4Z4 allele causes facioscapulohumeral muscular dystrophy type 2. *Nature Genetics*, *44*(12), 1370–1374. https://doi.org/10.1038/ng.2454

Lemmers, R. J. L. F., van der Vliet, P. J., Balog, J., Goeman, J. J., Arindrarto, W., Krom, Y. D., Straasheijm, K. R., Debipersad, R. D., Ozel, G., Sowden, J., Snider, L., Mul, K., Sacconi, S., van Engelen, B., Tapscott, S. J., Tawil, R., & van der Maarel, S. M. (2018). Deep characterization of a common D4Z4 variant identifies biallelic DUX4 expression as a modifier for disease penetrance in FSHD2. *European Journal of Human Genetics*, *26*(1), 94–106. https://doi.org/10.1038/s41431-017-0015-0

Lemmers, R. J. L. F., van der Vliet, P. J., Blatnik, A., Balog, J., Knijff, P., Tapscott, S. J., & van der Maarel, S. M. (2021). Chromosome 10q-linked FSHD identifies DUX4 as principal disease gene. *Journal of Medical Genetics*, *59*(2), 180–188. https://doi.org/10.1136/jmedgenet-2020-107474

Moriatis, J. S., Himeda, C. L., Ma, K. A. C., & Jones, P. L. (2024). Exchange of subtelomeric regions between chromosomes 4q and 10q reverts the FSHD genotype and phenotype. *Science Advances*, *10*(18), Article eadl1922. https://doi.org/10.1126/sciadv.adl1922

Strafella, C., Megalizzi, D., Trastulli, G., Proietti Piorgo, E., Colantoni, L., Tasca, G., Monforte, M., Zampatti, S., Primiano, G., Sancricca, C., Bortolani, S., Torchia, E., Ravera, B., Torri, F., Gadaleta, G., Risi, B., Caria, F., Gerardi, F., Carraro, E., ... Giardina, E. (2024). Integrating D4Z4 methylation analysis into clinical practice: Improvement of FSHD molecular diagnosis through distinct thresholds for 4qA/4qA and 4qA/4qB patients. *Clinical Epigenetics*, *16*, Article 148. https://doi.org/10.1186/s13148-024-01747-2

van der Vliet, P. J., Lemmers, R. J. L. F., Balog, J.,"; Wohlgemuth, M., Baechner, L., Tapscott, S. J., Tawil, R., & van der Maarel, S. M. (2015). A unifying genetic model for facioscapulohumeral muscular dystrophy. *Science*, *375*, 1650. (see also *Current Opinion in Neurology*, *28*, 567–575)

### FSHD long-read pipelines

Butterfield, R. J., Dunn, D. M., Duval, B., Moldt, S., & Weiss, R. B. (2023). Deciphering D4Z4 CpG methylation gradients in fascioscapulohumeral muscular dystrophy using nanopore sequencing. *Genome Research*, *33*, 1439–1454. https://doi.org/10.1101/gr.277871.123

Huang, M., Zhang, Q., Jiao, J., Shi, J., Xu, Y., Zhang, C., Zhou, R., Liu, W., Liang, Y., Chen, H., Wang, Y., Xu, Z., & Hu, P. (2024). Comprehensive genetic analysis of facioscapulohumeral muscular dystrophy by Nanopore long-read whole-genome sequencing. *Journal of Translational Medicine*, *22*, Article 451. https://doi.org/10.1186/s12967-024-05259-8

Lowenstern, T., Madritsch, S., Horner, D., Brait, N., Lafci, N. G., Schachner, A., Bujalkova, M. G., Hengstschlager, M., Dremsek, P., & Laccone, F. (2026). DUCKS4: A comprehensive workflow for Nanopore sequencing analysis of facioscapulohumeral muscular dystrophy (FSHD). *Human Genomics*, *20*, Article 21. https://doi.org/10.1186/s40246-026-00921-2

Lowenstern, T., Madritsch, S., Horner, D., Brait, N., Lafci, N. G., Schachner, A., Bujalkova, M. G., Hengstschlager, M., Dremsek, P., & Laccone, F. (2025). Rethinking genomics of facioscapulohumeral muscular dystrophy in the telomere-to-telomere era: pitfalls in the hidden landscape of D4Z4 repeats. *European Journal of Human Genetics*. https://doi.org/10.1038/s41431-025-02000-x

Xiao, L. C., Semwal, A., St John, B., Zeglinski, K., Su, S., Lancaster, J., Xue, S., Reversade, B., Ritchie, M. E., Magdinier, F., Blewitt, M. E., & Gouil, Q. (2025). D4Z4End2End: Complete genetic and epigenetic architecture of D4Z4 macrosatellites in FSHD, BAMS and reference cohorts. *medRxiv*. https://doi.org/10.1101/2025.04.24.25326320

### Reference genome

Nurk, S., Koren, S., Rhie, A., Rautiainen, M., Bzikadze, A. V., Mikheenko, A., Vollger, M. R., Altemose, N., Uralsky, L., Gershman, A., Aganezov, S., Hoyt, S. J., Diekhans, M., Logsdon, G. A., Alonge, M., Antonarakis, S. E., Borchers, M., Bouffard, G. G., Brooks, S. Y., ... Phillippy, A. M. (2022). The complete sequence of a human genome. *Science*, *376*(6588), 44–53. https://doi.org/10.1126/science.abj6987

### Bioinformatics tools

Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: Architecture and applications. *BMC Bioinformatics*, *10*, Article 421. https://doi.org/10.1186/1471-2105-10-421

Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, *10*(2), Article giab008. https://doi.org/10.1093/gigascience/giab008

Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*, *35*(4), 316–319. https://doi.org/10.1038/nbt.3820

Li, H. (2018). Minimap2: Pairwise alignment for nucleotide sequences. *Bioinformatics*, *34*(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191

Lin, J.-H., Chen, L.-C., Yu, S.-C., & Huang, Y.-T. (2022). LongPhase: An ultra-fast chromosome-scale phasing algorithm for small and large variants. *Bioinformatics*, *38*(7), 1816–1822. https://doi.org/10.1093/bioinformatics/btac058

Martin, M., Patterson, M., Garg, S., Fischer, S. O., Pisanti, N., Klau, G. W., Schonhuth, A., & Marschall, T. (2023). Read-based phasing and analysis of phased variants with WhatsHap. In *Methods in Molecular Biology* (Vol. 2590, pp. 127–138). Springer. https://doi.org/10.1007/978-1-0716-2819-5_8

Oxford Nanopore Technologies. (2024a). *Modkit* [Computer software]. https://github.com/nanoporetech/modkit

Oxford Nanopore Technologies. (2024b). *Dorado* [Computer software]. https://github.com/nanoporetech/dorado

Zheng, Z., Li, S., Su, J., Leung, A. W.-S., Lam, T.-W., & Luo, R. (2022). Symphonizing pileup and full-alignment for deep learning-based long-read variant calling. *Nature Computational Science*, *2*(12), 797–803. https://doi.org/10.1038/s43588-022-00387-x
