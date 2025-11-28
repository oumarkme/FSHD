# FSHD

A bioinformatics pipeline for analyzing **Facioscapulohumeral Muscular Dystrophy (FSHD)** using long-read sequencing data (e.g., Oxford Nanopore). This pipeline counts D4Z4 repeats and analyzes methylation levels in the 4q35 region.

## How It Works

### Overview

FSHD is caused by reduced D4Z4 repeat units on chromosome 4q35, typically dropping from 11–100 repeats to 1–10. This pipeline:

1. **Extracts** reads from the FSHD-associated region on chromosome 4
2. **Identifies** reads containing specific marker sequences (head/tail regions)
3. **Counts** D4Z4 repeat units per read using MIR8078 as a marker
4. **Determines** permissiveness (4qA vs 4qB haplotype)
5. **Analyzes** methylation levels in the DUX4 promoter region

### Workflow Steps

```
BAM files → Extract chr4 region → Create FASTA → Build BLAST DB
                                                       ↓
                                            BLAST for marker sequences
                                                       ↓
                              ┌─────────────────────────┴─────────────────────────┐
                              ↓                         ↓                         ↓
                         Head markers              Tail markers              D4Z4 counts
                         (DBET/P13E11)            (Exon2/pLAM)              (MIR8078)
                              ↓                         ↓                         ↓
                              └─────────────────────────┬─────────────────────────┘
                                                        ↓
                                              Classify reads & count D4Z4
                                                        ↓
                                              Check permissiveness (poly-A)
                                                        ↓
                                              Analyze methylation levels
                                                        ↓
                                              Generate output (TSV + plots)
```

### Key Marker Sequences

| Marker | Purpose | Location |
|--------|---------|----------|
| DBET | Head marker (v2) | Upstream of D4Z4 array |
| P13E11 | Head marker (v1) | Upstream of D4Z4 array |
| Exon 2 | Tail marker (v2) | DUX4 exon 2 |
| pLAM | Tail marker (v1) | Downstream of D4Z4 array |
| MIR8078 | D4Z4 repeat counter | Within each D4Z4 unit |
| pLAM-polyA | Permissiveness check | 4qA/4qB haplotype determination |

### Pipeline Versions

- **v1** (`fshd_v1.py`): Uses P13E11 (head) and pLAM (tail) markers
- **v2** (`fshd_v2.py`): Uses DBET (head) and Exon 2 (tail) markers

## Usage

```bash
python pipeline/fshd_v2.py -p PREFIX
```

Where `PREFIX` is your sample identifier (e.g., `WGS-ONT-A0003-sample-ID`).

### Arguments

| Argument | Description | Required |
|----------|-------------|----------|
| `-p`, `--prefix` | Sample prefix/identifier | Yes |

## Requirements

### Python Packages
- pandas
- plotnine

### System Software
- samtools (≥1.21 recommended)
- NCBI BLAST+ (`makeblastdb`, `blastn`)
- pb-CpG-tools (`aligned_bam_to_cpg_scores`)

## Output Files

| File | Description |
|------|-------------|
| `FSHD-{PREFIX}.tsv` | D4Z4 counts per read with target classification and permissiveness |
| `DUX4.cram` | Filtered alignments containing head/tail markers |
| `methylation.jpeg` | Boxplot of methylation levels in DUX4 promoter region |

### TSV Output Columns

| Column | Description |
|--------|-------------|
| Seq-id | Read identifier |
| Count | Number of D4Z4 repeats detected |
| Target | Classification: `both`, `head`, `tail`, or `none` |
| Permissive | `Permissive (4qA)`, `Non-permissive (4qB)`, or `Unknown` |

## Protocol Details

1. **Extract Reads**: From BAM files, extract reads mapped to `chr4:190000000-190175845` (GRCh38/hg38)

2. **Create FASTA**: Convert extracted reads to FASTA format

3. **Build BLAST Database**: Create a nucleotide BLAST database from the reads

4. **Identify Marker Sequences**:
   - **Tail sequences**: Exon 2 (v2) or pLAM (v1)
   - **Head sequences**: DBET (v2) or P13E11 (v1)

5. **Count D4Z4 Repeats**: Count occurrences of MIR8078 sequence in each read

6. **Determine Permissiveness**: Check for poly-A signal to distinguish 4qA (permissive) from 4qB (non-permissive) haplotype

7. **Methylation Analysis**: 
   - Extract methylation levels using `aligned_bam_to_cpg_scores`
   - Analyze region `chr4:190175600-190175834` (GRCh38/hg38, DUX4 promoter)

## Directory Structure

```
FSHD/
├── pipeline/
│   ├── fshd_v1.py      # Pipeline version 1
│   └── fshd_v2.py      # Pipeline version 2 (recommended)
├── querySeq/           # Query sequences for BLAST
│   ├── DBET.fasta      # Head marker (v2)
│   ├── P13E11.fasta    # Head marker (v1)
│   ├── exon2.fasta     # Tail marker (v2)
│   ├── pLAM.fasta      # Tail marker (v1)
│   ├── MIR8078.fasta   # D4Z4 repeat marker
│   ├── pLAM-polyA.fasta # Permissiveness check
│   └── meth_reg.fasta  # Methylation region
├── test/               # Test scripts and notebooks
└── HG38_old/          # Legacy analysis files
```

## References

- FSHD is associated with contractions of D4Z4 repeats on chromosome 4q35
- Permissive (4qA) haplotype is required for DUX4 expression
- Methylation of the D4Z4 region affects disease manifestation