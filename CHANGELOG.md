# Changelog

## v5.1 (2026-06-02)

### Added — dual-reference input (hg38 or T2T)
- Input BAM/CRAM may now be aligned to **hg38 or T2T-CHM13v2.0**. `--input_ref auto` (default) detects the build from the `@SQ` chr4 length, routing to whichever canonical length (hg38 190,214,555 / T2T 193,574,945) is closer; `--input_ref hg38|t2t` forces it.
- New `DETECT_REF` process (header-only build detection) and `SUBSET_T2T` process (subsets an already-T2T-aligned input to the D4Z4 region, replacing Steps 1–2 — no realignment). Both lanes converge on the same T2T BAM feeding the unchanged Steps 3–6.
- A T2T CRAM is decoded with `--ref_t2t`; `--ref_hg38` is now needed only to decode an **hg38 CRAM**.

### Changed
- D4Z4 region coordinates centralized in `assets/d4z4_regions.bed`, shared by `EXTRACT_D4Z4` and `SUBSET_T2T` (the hg38-flanking and methylation-target BEDs are unchanged).
- Process `stub`/`script` sections reordered (`script` first) to satisfy the strict Nextflow grammar (`nextflow lint` / language server); no behavioural change.

## v5.0 (2026-06-01)

### Renamed
- Pipeline renamed `nf-fshd` → **FSHDetector** (directory `FSHDetector/`); repository slug unchanged (`oumarkme/FSHD`).

### Algorithm (rewrite)
- BLAST-based per-read annotation replaces v4's consensus/assembly counting. Each read is annotated independently against three queries (D4Z4 repeat unit, 6 haplotype-specific markers, p13E-11).
- Three-tier chromosome assignment: T2T alignment → BlnI/XapI restriction-motif voting → haplotype-marker concordance (strongest).
- Seven-value haplotype classification (4qA-S / 4qA-L / 4qA / 4qB / 4q / 10qA / 10qB / 10q), with 4qA-S vs 4qA-L decided by terminal D4Z4 unit length.
- Sequence-based permissiveness (functional 4qA pLAM ATTAAA PAS, non-chimeric) and chimeric/translocation detection.
- Biological validation report (`validation.txt`) flagging impossible configurations (haplotype count, completeness, RU consistency, chimerism).

### Output (7 files per sample)
- `d4z4.bam(.bai)`, `complete_reads.bam(.bai)`, `repeats.tsv`, `methylation.tsv`, `validation.txt`.

### Environment & portability (public release)
- **Conda** is the only engine (`environment.yml`) — no root, manual installs, or containers required. Two profiles only: `conda` (default) and `test`.
- `run.sh` no longer hard-codes a site conda path; it is applied only when present, so external users rely on their own Nextflow/Conda.
- `--ref_hg38` is now required only for CRAM input; BAM input may omit it.
- Added a `-profile test` (`-stub-run`) wiring smoke-test, per-process `stub` blocks, `tag` directives, and `checkIfExists` input validation.

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
