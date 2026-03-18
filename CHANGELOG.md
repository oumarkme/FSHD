# Changelog

## v4 (2026-03-18)

### Major changes
- Three-phase allele assignment algorithm (HP → RU-anchor → Jenks → single)
- Alignment-based unphased read assignment (minimap2 + RU-constraint)
- Strategy-based routing: read-level (Racon consensus) vs assembly (Flye vs raw-read)
- Per-allele BLAST counting replaces fixed chr4-HAP1-only read-level counting
- Flye-vs-raw-read reference selection for assembly counting

### Output changes
- 4 published files only: d4z4.bam, d4z4.bam.bai, repeats.tsv, methylation.tsv
- repeats.tsv: grouped by (Allele, Haplotype, Complete) with RU_Count, Count_Method, Reference_Source
- methylation.tsv: includes reference atlas Mean(SD) and Caputo et al. diagnostic thresholds

### Container
- Removed: bcftools, seqtk, plotnine, whatshap (unused)
- Kept: samtools, minimap2, blast, modkit, flye, racon, seqkit, pysam, pandas, jenkspy

## v3 (2025-03-16)

- Initial tested version with HP-tag-primary allele assignment
- Fixed chr4 HAP1 read-level counting
- Assembly counting for all 4 allele combos unconditionally
- Archived in archive/nf-fshd-v3/
