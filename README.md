# FSHDetector

D4Z4 repeat counting, haplotyping, and methylation from Oxford Nanopore (ONT)
whole-genome reads, for research on facioscapulohumeral muscular dystrophy
(FSHD).

Version 6.0.0 · experimental, research use only, not a validated diagnostic
test · MIT license. For the version in our manuscript, see
[Versions](#versions).

FSHDetector takes one aligned BAM or CRAM file and reports two alleles on chr4
and two on chr10. Each allele has a D4Z4 repeat unit (RU) count, the evidence
behind the count, a haplotype, and CpG methylation. Everything runs from a
single container image.

## Install

You need [Apptainer](https://apptainer.org) on Linux (x86-64).

```bash
git clone https://github.com/oumarkme/FSHD.git
cd FSHD
apptainer build --fakeroot fshdetector.sif env/fshdetector.def
apptainer test fshdetector.sif
```

The test runs the unit tests and one complete analysis of a small synthetic
genome.

## Quick start

```bash
apptainer run --bind /data fshdetector.sif run \
    --input /data/sample.bam \
    --ref-t2t /data/hs1.fa \
    --outdir /data/results \
    --threads 16
```

The result is `/data/results/sample.alleles.tsv`, with two rows per chromosome.
Selected columns for the synthetic test genome are:

```text
chromosome  allele  haplotype  permissive  reported_ru  source    spanning_reads
chr4        1       4qA        yes         3            spanning  8
chr4        2       4qA        yes         6            spanning  8
chr10       1       10qA       no          4            spanning  8
chr10       2       10qA       no          8            spanning  8
```

## Usage

### Input

| Option | Requirement |
|--------|-------------|
| `--input` | Indexed BAM or CRAM, aligned to hg38 or T2T-CHM13v2.0 with `chr` contig names. Keep the `MM` and `ML` tags for methylation. Super-accuracy basecalls are recommended. |
| `--ref-t2t` | Always required. Telomere-to-telomere (T2T) CHM13v2.0 FASTA (`GCF_009914755.1`) with `chr` contig names and a `.fai` index, for example [UCSC `hs1`](https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/). |
| `--ref-hg38` | The hg38 FASTA of the alignment. Required for hg38 CRAM input. |
| `--t2t-index` | Optional minimap2 index for hg38 input, built once with the `index` command. It saves time in every run. |

Other options are `--outdir`, `--sample`, `--threads` (default 8), `--tmpdir`,
`--build`, and `--validate-reassembly`. List them with
`apptainer run fshdetector.sif run --help`.

To create the two index files:

```bash
apptainer exec --bind /data fshdetector.sif \
    python -c 'import pysam; pysam.faidx("/data/hs1.fa")'
apptainer run --bind /data fshdetector.sif index \
    --ref-t2t /data/hs1.fa --output /data/hs1.mmi
```

The genome build of the input comes from the chr4 length in its header. hg19
input is rejected. A file that cannot be decoded stops the run and names the
failing region, because a partial read set would give wrong counts.

### Reading the result

`reported_ru` comes from the strongest evidence available, named in `source`. A
weaker tier never changes a stronger one.

| `source` | Evidence | `reported_ru` |
|----------|----------|---------------|
| `spanning` | Reads that hold both ends of the array. Primary evidence. | Exact |
| `reassembly` | Partial reads tiled by the variants that differ between repeat units. Supporting, experimental. | Exact |
| `longest_read` | The longest read that the other allele cannot explain. | Lower bound, `>=N` |
| `none` | No evidence for a second allele. | `NA` |

Check the `flags` column before you use a count. [Output files, columns, and
flags](references/output.md) are described separately.

## Methodology

1. **Collect reads** (pysam 0.24.1). Read the two D4Z4 regions and the unmapped
   reads, recover full sequences of split reads, and keep reads that share a
   21-mer with a query. For hg38 input, reads from chr4:189.0-190.5 Mb and
   chr10:133-134 Mb are realigned to T2T with minimap2 2.31 (`lr:hq`).
2. **Annotate** (BLAST+ 2.17.0). Search the repeat unit, six haplotype markers,
   and the p13E-11 probe. A unit counts if its alignment exceeds 1,500 bp. A
   poorly read unit between two counted units counts too.
3. **Assign reads.** The chromosome comes from the alignment, restriction motifs
   inside the array, and the markers. A or B type needs a distal marker with at
   least 90% identity. Fold-back reads are left out.
4. **Describe units.** Learn the positions at which the repeat units of the
   sample differ, estimate the error rate at each, and reduce every unit to its
   variants.
5. **Call alleles.** Group spanning reads by count, type, and variants. Rebuild
   a missing allele from partial reads only if one array length fits, several
   reads cover every link, and random read subsets give no other length.
   Otherwise report a lower bound.
6. **Measure methylation** (modkit 0.6.4) on the reads of each allele.

### Accuracy

We simulated 8,000 samples that no threshold was tuned on
([`benchmarks/simulate.py`](benchmarks/simulate.py)): T2T-CHM13v2.0 repeat
units, 2 to 30 units per allele, 8x to 30x per allele, and a per-base error of
about 0.4% to 5.5%. The samples include mosaic alleles, chr4/chr10 unit
exchanges, chimeric reads, and misaligned reads.

| Reported value | Correct | 95% CI (Wilson) |
|----------------|---------|-----------------|
| `spanning`, 2 or more reads | 20,175 of 20,177 (99.99%) | 99.96-100 |
| `spanning`, 1 read | 3,514 of 3,599 (97.6%) | 97.1-98.1 |
| `reassembly` | 1,126 of 1,131 (99.6%) | 99.0-99.8 |
| `contig_lower_bound` holds | 11,477 of 11,509 (99.7%) | 99.6-99.8 |
| `longest_read` bound holds | 5,703 of 5,821 (98.0%) | 97.6-98.3 |

Reassembly resolved 13% of the alleles that no read spans: 21% below 1.1%
per-base error, 14% from 1.1% to 2.8%, and 2% above. Simulation does not replace
validation against an orthogonal method such as Southern blot or optical
mapping, which has not been done for this version.

On two whole-genome samples, a run took 13 to 21 s from a T2T CRAM (16 threads)
and 42 to 63 s from an hg38 CRAM (8 threads, 13 GB peak memory, index built
beforehand).

### Limitations

* An exact count depends on reads that cross the array. Longer alleles usually
  get a lower bound only.
* Identical repeat units cannot be placed, so reassembly gives no count for
  them. It rarely succeeds at low coverage.
* A count from a single spanning read is wrong more often (2.4% in simulation).
* Arrays that carry units from the other chromosome can be assigned to the wrong
  chromosome.
* No chr10 B-type marker is included, so 10qB is never called.
* Only two alleles per chromosome are reported. Further groups, as in mosaic
  samples, are listed in `extra_alleles`.
* The reference atlas values in `methylation.tsv` were not recomputed for this
  version.

## Development

```bash
pip install -e '.[test]'    # needs a C compiler
pytest                      # the end-to-end test needs BLAST+ and modkit
```

## Versions

| Version | Where | Use it for |
|---------|-------|------------|
| **v5.1.0 (manuscript version)** | tag [`v5.1.0`](https://github.com/oumarkme/FSHD/tree/v5.1.0) (commit `28f10b6`) | Running the version that our manuscript describes |
| v6.0.0 and later | `main` branch | Current work. A Python rewrite with a new command line; results differ from v5.1.0 |

To use the manuscript version, run `git checkout v5.1.0` after cloning, and set
`--input_ref hg38` or `--input_ref t2t`. Its automatic build detection can stop
with an error on files that have long headers. [`CHANGELOG.md`](CHANGELOG.md)
describes a chromosome assignment error in v5.x that version 6.0.0 fixes.

## References

* Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer,
  K., & Madden, T. L. (2009). BLAST+: Architecture and applications. *BMC
  Bioinformatics, 10*, 421. <https://doi.org/10.1186/1471-2105-10-421>
* Caputo, V., Megalizzi, D., Fabrizio, C., et al. (2022). D4Z4 methylation
  levels combined with a machine learning pipeline highlight single CpG sites
  as discriminating biomarkers for FSHD patients. *Cells, 11*(24), 4114.
  <https://doi.org/10.3390/cells11244114>
* Li, H. (2018). Minimap2: Pairwise alignment for nucleotide sequences.
  *Bioinformatics, 34*(18), 3094-3100.
  <https://doi.org/10.1093/bioinformatics/bty191>
* Nurk, S., Koren, S., Rhie, A., et al. (2022). The complete sequence of a
  human genome. *Science, 376*(6588), 44-53.
  <https://doi.org/10.1126/science.abj6987>
* Oxford Nanopore Technologies. *modkit* (version 0.6.4) [Computer software].
  <https://github.com/nanoporetech/modkit>

## Citation, license, and help

Cite as: Ou, J.-H. *FSHDetector* (version 6.0.0) [Computer software].
<https://github.com/oumarkme/FSHD>

MIT license, see [`LICENSE`](LICENSE). Report problems on the
[issue tracker](https://github.com/oumarkme/FSHD/issues). Do not attach patient
data to an issue.
