# Output files

All files are written to `--outdir`. Except for `fshd_feature_t2t.bed`, their
names start with the sample name.

| File | Content |
|------|---------|
| `alleles.tsv` | The result: two rows per chromosome |
| `read_alleles.tsv` | Reads assigned to each allele |
| `unit_methylation.tsv` | 5mC per repeat unit of each allele |
| `methylation.tsv` | DR1 and DUX4-PAS 5mC for each upstream `HP` tag and for all reads |
| `reads_annotated.tsv` | One row per read: chromosome, haplotype, unit count, array ends, markers |
| `repeats.tsv` | Largest unit count and read number for each chromosome, haplotype, and array end |
| `validation.txt` | Consistency checks of the read calls |
| `d4z4.bam` | Candidate reads |
| `complete_reads.bam` | Candidate reads that hold both array ends |
| `fshd_feature_t2t.bed` | T2T-CHM13v2.0 features of the chr4 D4Z4 locus, for a genome browser |
| `reassembly_validation.tsv` | Written with `--validate-reassembly`: reassembly tested against counts from spanning reads |
| `run.json` | Version, options, tool versions, read counts, timings, and warnings |

## `alleles.tsv`

`NA` means no data. A repeat unit (RU) is one D4Z4 unit.

| Column | Meaning |
|--------|---------|
| `sample`, `chromosome`, `allele` | Sample name, chr4 or chr10, and row 1 or 2 |
| `haplotype`, `permissive` | A or B type from the assigned reads, and whether the haplotype is permissive |
| `reported_ru`, `source` | The count to use and the evidence tier it comes from |
| `spanning_ru`, `spanning_reads` | Count from reads that hold both array ends, and their number |
| `reassembly_ru`, `weakest_link_reads` | Count from reassembly, and the reads that cover its weakest link |
| `contig_lower_bound` | Units of an assembled contig, linked from its first unit by several reads, that shares reads with the allele |
| `longest_read_ru` | Largest unit count among the assigned reads |
| `supporting_reads` | Partial reads assigned to the allele |
| `conflicting_reads` | Partial reads that match the allele but carry more units |
| `DR1_methylation`, `DR1_CpG` | Mean 5mC (%) of DR1 and the number of CpGs (chr4) |
| `DUX4-PAS_methylation`, `DUX4-PAS_CpG` | The same for the DUX4 polyadenylation signal (chr4) |
| `distal_units_methylation` | Mean 5mC (%) of the three most distal units |
| `distal_units_calls`, `distal_units_reads` | CpG calls and reads behind that value |
| `extra_alleles` | Spanning groups beyond the two reported alleles |
| `flags` | Warnings, listed below |

On T2T-CHM13v2.0, DR1 is chr4:193,541,000-193,541,600 and DUX4-PAS is
chr4:193,543,430-193,543,634 (0-based, half-open).

## Flags

| Flag | Meaning |
|------|---------|
| `experimental` | The count comes from reassembly |
| `single_read` | One read supports the allele |
| `extra_allele` | More than two spanning groups (see `extra_alleles`) |
| `split_by_barcode`, `possible_split`, `unplaced_spanning_read` | Spanning reads of one count disagree in their variants |
| `count_discordance` | A read that differs by one unit was merged into the group |
| `mixed_haplotype` | Reads of the allele carry both A-type and B-type calls |
| `assembly_conflict`, `unexplained_reads` | Partial reads contradict or do not fit the called alleles |
| `suspected_chimera_excluded` | A single spanning read looked like two joined molecules and was left out |
| `unsupported_read` | A lower bound rests on a read that no other read supports |
| `second_allele_not_seen`, `no_reads` | Why a row is empty |

## `methylation.tsv`

Each region has one row for `HP1`, `HP2`, and `Combined` reads, with the mean 5mC
(%) and the number of CpGs. The remaining columns are reference values for
comparison: an in-house ONT control atlas and the thresholds of Caputo et al.
(2022), *Cells, 11*(24), 4114. <https://doi.org/10.3390/cells11244114>
