# FSHD D4Z4 Repeat Analysis Pipeline

ONT (Oxford Nanopore) long-read sequencing pipeline for FSHD (Facioscapulohumeral Dystrophy) D4Z4 repeat counting and methylation analysis on chromosome 4q35.

## Usage

```bash
python pipeline/fshd_v2.py -p PREFIX
```

`PREFIX` is the sample identifier (e.g., ONT019). v2 is the current version.

## Pipeline Stages

1. **BAM merge** — Find all BAMs under `/nanopore_nas/nanopore/{PREFIX}/`, merge reads from `chr4:190000000-190175845`, sort and index
2. **FASTA extraction** — Convert merged BAM to FASTA
3. **BLAST DB** — Build a nucleotide BLAST database from extracted reads
4. **Head/tail identification** — BLAST query sequences to identify reads spanning the D4Z4 array:
   - Head (proximal): DBET (v2) / P13E11 (v1)
   - Tail (distal): exon2 (v2) / pLAM (v1)
5. **D4Z4 counting** — Count D4Z4 repeats per read using MIR8078 as the repeat marker. Reads are labeled `both`/`head`/`tail`/`none` based on which flanking sequences they contain
6. **Permissive check** — BLAST pLAM-polyA (hapA/hapB) to determine 4qA (permissive) vs 4qB (non-permissive) haplotype based on a single-nucleotide difference (ATTAAA vs ATCAAA)
7. **CRAM output** — Subset BAM to reads with head or tail hits, save as CRAM
8. **Methylation** — Run `aligned_bam_to_cpg_scores` (pb-CpG-tools) on the CRAM, filter to region `chr4:190175600-190175834`, plot methylation boxplot

## v1 vs v2 Differences

| | v1 | v2 |
|---|---|---|
| Head marker | P13E11 | DBET |
| Tail marker | pLAM | exon2 |

All other logic is identical.

## Query Sequences (`querySeq/`)

| File | Purpose |
|---|---|
| `DBET.fasta` | Head (proximal) marker for D4Z4 array (v2) |
| `exon2.fasta` | Tail (distal) marker for D4Z4 array (v2) |
| `P13E11.fasta` | Head marker (v1, replaced by DBET) |
| `pLAM.fasta` | Tail marker (v1, replaced by exon2) |
| `MIR8078.fasta` | D4Z4 repeat unit marker — BLAST hit count = repeat count |
| `pLAM-polyA.fasta` | Two haplotype sequences (hapA/hapB) for permissive/non-permissive determination |
| `meth_reg.fasta` | Methylation region reference (not used in current pipeline) |

## External Dependencies (hardcoded paths)

- **samtools 1.21**: `/Everythings/misc/samtools-1.21/samtools`
- **BLAST 2.13.0+**: `/Everythings/misc/ncbi-blast-2.13.0+/bin/`
- **pb-CpG-tools v2.3.2**: `/Everythings/misc/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores`
- **Pileup model**: `/nanopore_nas/pipeline/data/pileup_calling_model.v1.tflite`

## Resource Paths

- **Input BAMs**: `/nanopore_nas/nanopore/{PREFIX}/`
- **Output**: `/nanopore_nas/nanopore_pipeline/{PREFIX}/60.FSHD/`
- **Reference genome**: `/nanopore_nas/pipeline/data/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set_modified.fasta` (GRCh38)
- **Temp directory**: `/raid/tmp/FSHD/` (cleaned up at end)

## Python Dependencies

- pandas
- plotnine
- argparse (stdlib)

## Other Directories

- `HG38_old/` — Earlier analysis notebooks (ONT019, ONT022, ONT026)
- `test/` — Experimental scripts (DBET BLAST test, T2T realignment)
- `archive/` — Unrelated old project, gitignored
- `doc/` — Background documentation and pipeline design notes

---

## v3 Pipeline Plan (T2T-CHM13v2.0 rewrite)

The current pipeline (v1/v2) operates entirely on GRCh38 with simple BLAST-based counting. The next version (v3) will realign to T2T-CHM13v2.0 and add 4q/10q disambiguation, allele-specific methylation, and FSHD2 gene screening. Full design rationale is in `doc/fshd_info.md`.

### Why rewrite

- GRCh38 has an ~82 kb assembly gap at 4q35 — many D4Z4 reads go unmapped
- Current pipeline cannot distinguish 4q from 10q D4Z4 arrays (>98% identity)
- No FSHD2 detection (requires SMCHD1/DNMT3B/LRIF1 variant screening)
- No allele-specific methylation analysis
- No per-allele repeat counting

### Step 1: Targeted read extraction + T2T realignment

1. Extract reads from hg38 BAM at FSHD-relevant regions (4q35, 10q26, SMCHD1, DNMT3B, LRIF1, D4Z4-like paralog loci) plus unmapped reads
2. Revert to unaligned BAM preserving MM/ML methylation tags (`samtools reset --keep-tag MM,ML`)
3. Align to T2T-CHM13v2.0 with Winnowmap2 (preferred) or minimap2 ≥v2.22
   - Winnowmap2 handles segmental duplications far better (FNR 1.89% vs minimap2's 39.62%)
   - `-y` flag propagates methylation tags; `--eqx` for PSV analysis

**Key T2T coordinates:**
| Region | Coordinates |
|---|---|
| Chr4 D4Z4 + flanking | chr4:193,282,466–193,558,261 |
| Chr4 DUX4 methylation target | chr4:193,540,172–193,543,634 |
| Chr4 p13E-11 anchor | chr4:192,667,301–192,902,247 |
| Chr10 D4Z4 array | chr10:134,615,120–134,741,148 |

### Step 2: 4q vs 10q disambiguation

1. Extract reads mapping to either D4Z4 locus ± 50 kb flanking
2. BLAST chromosome-specific feature sequences against each read:
   - D4F104S1 (p13E-11): 4q vs 10q variants
   - BlnI motif (CCTAGG): present in 10q, absent in 4q
   - XapI motif (RAATTY): present in 4q, absent in 10q
   - pLAM: ATTAAA (4qA) vs ATCAAA (10qA)
3. Classify by BLAST bitscore ratio (1.2× threshold)
4. Optional: DuploMap PSV-based refinement for ambiguous reads

### Step 3: D4Z4 repeat counting

1. Per-read counting via BLAST (3.3 kb D4Z4 unit as query, count non-overlapping hits)
2. Classify reads as complete (both anchors), partial-distal, or partial
3. Alternative/validation: tandem-genotypes with LAST aligner
4. Assembly-based validation with Flye for ultra-long coverage samples

**Clinical threshold:** ≤10 D4Z4 units on a 4qA allele = candidate FSHD1

### Step 4: Methylation analysis

1. Per-site methylation with modkit (replaces pb-CpG-tools) at distal D4Z4/DUX4 region
2. Allele-specific methylation via read-level phasing:
   - Variant calling with PEPPER-Margin-DeepVariant
   - Phasing with WhatsHap
   - Haplotype-separated modkit pileup
3. Key diagnostic thresholds:
   - DR1 mean CpG methylation: healthy >40%, FSHD2 <30%
   - DUX4-PAS CpG6: healthy >72%, FSHD <72%

### Step 5: Haplotype determination (4qA vs 4qB)

1. BLAST pLAM against chr4-assigned distal reads
   - ATTAAA → 4qA (permissive, FSHD-relevant)
   - ATCAAA → 10qA (confirms chr10)
   - No pLAM → 4qB (non-permissive)
2. Sub-haplotype classification using SSLP + D4F104S1 SNPs + distal A/B
3. Dominant FSHD haplotype: 4A161 (~85% of FSHD cases)

### Step 6: FSHD2 gene screening

1. Call variants at SMCHD1 (chr18), DNMT3B (chr20), LRIF1 (chr1) with PEPPER-Margin-DeepVariant
2. Annotate with VEP/ClinVar
3. SMCHD1 mutated in ~80% of FSHD2; dominant-negative ATPase mutations = more severe

### Step 7: Diagnostic classification

Decision tree output per sample:
- **FSHD1**: 4qA allele with ≤10 D4Z4 units + hypomethylation on contracted allele
- **FSHD2**: Normal repeat count + genome-wide D4Z4 hypomethylation + SMCHD1/DNMT3B/LRIF1 variant + 4qA
- **Unaffected**: All 4qA alleles ≥11 units, normal methylation

Confidence levels: High (≥3 complete spanning reads), Medium (≥1 complete or ≥3 partial), Low (minimum count only)

### New dependencies (v3)

| Tool | Purpose |
|---|---|
| Winnowmap2 v2.03 | T2T alignment (segmental dup-aware) |
| modkit | ONT methylation pileup + per-read extraction |
| PEPPER-Margin-DeepVariant | SNV/indel calling for phasing + FSHD2 screening |
| WhatsHap | Read phasing + haplotagging |
| DuploMap | PSV-based paralog disambiguation (optional) |
| tandem-genotypes + LAST | Alternative repeat counting (optional) |
| Flye/Shasta | Local assembly validation (optional) |
| Sniffles2 | Population SV panel from 100-sample reference |
| bcftools | Variant filtering |
| VEP | Variant annotation |

### 100-sample reference panel

Process 100 non-FSHD WGS samples through the full pipeline to establish:
- Baseline D4Z4 repeat count distributions (expect ≥11 units on all 4qA alleles)
- Baseline methylation ranges (DR1 >40%, DUX4-PAS CpG6 >72%)
- PSV catalog for 4q/10q disambiguation (expect 50–150 reliable PSVs per 3.3 kb unit)
- Haplotype frequency table
- Population SV panel via Sniffles2

### Reference pipelines

- **DUCKS4** (Löwenstern et al. 2026): Docker-based, T2T, BLAST + modkit. github.com/tamara-nano/ducks4
- **D4Z4End2End** (Xiao et al. 2025): Full array assembly + consensus. github.com/lucindaxiao/D4Z4End2End
- **Huang et al. (2024)**: WGS-only ONT at ~29× with 100% OGM concordance
