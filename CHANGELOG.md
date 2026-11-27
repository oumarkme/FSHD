# Changelog

## v6.0.0 (unreleased)

A rewrite in Python and Cython that runs from one container image. The command
line and several output tables changed. Version 5.1.0 stays available as the
manuscript version (tag `v5.1.0`).

### Changed
- Nextflow, Conda profiles, and `run.sh` are replaced by the `fshdetector run`
  and `fshdetector index` commands and one Apptainer definition (`env/`).
- Reads are collected through the alignment index and screened by k-mers; the
  whole file is no longer scanned.
- Result tables: new `alleles.tsv` (two alleles per chromosome, evidence tier,
  methylation, flags), `read_alleles.tsv`, `unit_methylation.tsv`, and
  `run.json`. `reads_annotated.tsv` has a `Mixed_Strand` column.

### Added
- Supporting evidence for alleles that no read spans: reassembly from partial
  reads by unit variants (experimental) and lower bounds.
- Methylation per allele and per repeat unit, assigned from the reads of each
  allele, so it no longer depends on upstream haplotype phasing.
- Unit tests, an end-to-end test, and a read simulator
  (`benchmarks/simulate.py`).

### Fixed
- Chromosome assignment counted restriction motifs over the whole read. XapI
  sites are frequent in flanking sequence, so chr10 reads could be assigned to
  chr4. Motifs are now counted inside the repeat array only. All earlier 5.x
  versions, including 5.1.0, have this error.
- Combined methylation counted phased reads twice.
- Extraction errors were hidden, so a damaged input could give a partial read
  set without notice. The run now stops and names the region.
- Read length and motifs could come from a clipped supplementary record; the
  primary record is now used, and primary records outside the regions are
  recovered.
- minimap2 hard-clipped supplementary records during hg38 realignment.
- BLAST reported at most 500 reads per query.
- A repeat unit split into two long alignments was counted twice.
- A fold-back read was counted in both orientations, which doubled its units.
- A threshold estimated by modkit was not applied to some inputs; it is now
  passed explicitly.
- hg19 input is rejected instead of being processed as hg38.

### Changed counting and typing rules
- A poorly read unit between two counted units is counted.
- Fold-back reads are flagged and excluded from allele calls.
- A single spanning read that looks chimeric is excluded.
- Permissive calls use all reads of an allele.

## v5.2 (2026-07-29)

### Fixed — confidence-gated haplotype classification
- `BLAST_ANNOTATE` previously assigned 4qA/4qB per read from the single best-bitscore distal marker hit with no confidence gate, and forced a 4qA/4qB label from a lone proximal D4F104S1 hit even with no distal (pLAM) evidence at all. Because 4qA_pLAM and 10qA_pLAM are >95% identical to each other, this let a chr4-assigned read whose *true* best distal hit was the chr10-type marker still get forced into `4qA` — producing biologically impossible 3-allele calls on one chromosome (confirmed on one internal sample: chr4 reported `4qA`, `4qA-L`, and `4qB` simultaneously).
- Haplotype output is now exactly 3 values per chromosome: `4qA` / `4qB` / `4q` (`10qA` / `10qB` / `10q` for chr10 — `4q`/`10q` is the pre-existing untyped label, now also the fallback for low-confidence evidence, not just "no marker hit at all"). A read is only called `4qA`/`4qB` when its best distal-marker hit is >=90% identity **and** its native chromosome (marker-implied) matches the read's assigned chromosome — the one documented exception being a non-chimeric chr10 read whose best hit is `4qA_pLAM` (rare 4qA-PAS translocation onto chr10; Lemmers et al., 2021), which stays typeable/permissive. A lone proximal marker, a chimeric read, or a chromosome-inconsistent/low-identity distal hit now falls back to the untyped `4q`/`10q` label instead of guessing.
- The terminal-D4Z4-unit length subtype (previously the `4qA-S`/`4qA-L` values in the `Haplotype` field) moved to a new `Length_Subtype` column in `reads_annotated.tsv` — it is length annotation, not an independent haplotype call.
- `is_permissive()` is now gated by the same confidence check.
- On that sample this collapses chr4 from the invalid 3-typed-haplotype WARN to a clean `PASS: 2 typed haplotype(s): ['4qA', '4qB']` (4qA: RU=4, 3 concordant complete reads, permissive; 4qB: RU up to 9, fragment evidence, non-permissive).

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
