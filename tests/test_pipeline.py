"""The whole run on a small synthetic genome with known alleles."""

import array
import csv
import json
import shutil

import conftest
import numpy as np
import pysam
import pytest

from fshdetector import layout
from fshdetector import pipeline
from fshdetector import sequences

pytestmark = pytest.mark.skipif(
    not all(map(shutil.which, ("blastn", "makeblastdb", "modkit"))),
    reason="BLAST+ and modkit are not on PATH",
)

TRUTH = {"chr4": (3, 6), "chr10": (4, 8)}
FLANK = 8000
READS_PER_ALLELE = 8


def _arrays(generator) -> dict[tuple[str, int], str]:
    """Allele sequences: flank, proximal marker, units, distal marker, flank."""
    unit = next(
        iter(sequences.read_fasta(sequences.data_path("d4z4.fa")).values())
    )
    markers = sequences.read_fasta(sequences.data_path("haplotype_markers.fa"))
    result = {}
    seed = 100
    for chromosome, counts in TRUTH.items():
        prefix = "4qA" if chromosome == "chr4" else "10qA"
        for count in counts:
            units = [
                conftest.substitute(unit, 0.01, generator) for _ in range(count)
            ]
            seed += 1
            result[chromosome, count] = "".join(
                [
                    conftest.random_dna(FLANK, seed),
                    markers[f"{prefix}_D4F104S1"],
                    *units,
                    markers[f"{prefix}_pLAM"],
                    conftest.random_dna(FLANK, seed + 50),
                ]
            )
    return result


def _cpg_tags(record: pysam.AlignedSegment, generator) -> None:
    """Marks every CpG as called, about half of them methylated."""
    sequence = record.get_forward_sequence()
    cytosines = [i for i, base in enumerate(sequence) if base == "C"]
    skips = []
    probabilities = []
    previous = -1
    for rank, position in enumerate(cytosines):
        if sequence[position : position + 2] != "CG":
            continue
        skips.append(str(rank - previous - 1))
        previous = rank
        probabilities.append(250 if generator.random() < 0.5 else 4)
    record.set_tag("MM", f"C+m?,{','.join(skips)};")
    record.set_tag("ML", array.array("B", probabilities))


@pytest.fixture(name="outdir", scope="module")
def fixture_outdir(tmp_path_factory):
    work = tmp_path_factory.mktemp("run")
    generator = np.random.default_rng(5)
    arrays = _arrays(generator)
    # The reference holds the longer allele of each chromosome.
    reference = {c: arrays[c, max(n)] for c, n in TRUTH.items()}
    fasta = work / "t2t.fa"
    sequences.write_fasta(fasta, reference.items())
    pysam.faidx(str(fasta))
    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": name, "LN": len(s)} for name, s in reference.items()],
        }
    )
    records = []
    for (chromosome, count), sequence in arrays.items():
        for i in range(READS_PER_ALLELE):
            read = conftest.substitute(sequence, 0.002, generator)
            if i % 2:
                read = sequences.reverse_complement(read)
            record = pysam.AlignedSegment(header)
            record.query_name = f"{chromosome}_{count}_{i}"
            record.query_sequence = (
                sequences.reverse_complement(read) if i % 2 else read
            )
            record.flag = 16 if i % 2 else 0
            record.reference_name = chromosome
            # Align the distal ends, so that every read covers the targets.
            length = min(len(read), len(reference[chromosome]))
            soft = len(read) - length
            record.reference_start = len(reference[chromosome]) - length
            record.cigartuples = ([(4, soft)] if soft else []) + [(0, length)]
            _cpg_tags(record, generator)
            records.append(record)
    bam = work / "input.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for record in sorted(
            records, key=lambda r: (r.reference_id, r.reference_start)
        ):
            out.write(record)
    pysam.index(str(bam))

    targets = tuple(
        layout.Region(c, len(s) - FLANK - 300, len(s) - FLANK, f"{c}_distal")
        for c, s in reference.items()
    )
    chr4_end = targets[0].end
    small = layout.Layout(
        chr4_lengths={"t2t": len(reference["chr4"])},
        t2t_windows=tuple(
            layout.Region(c, 0, len(s)) for c, s in reference.items()
        ),
        hg38_windows=(),
        meth_targets=targets,
        meth_regions=(
            layout.Region("chr4", chr4_end - 300, chr4_end - 150, "DR1"),
            layout.Region("chr4", chr4_end - 150, chr4_end, "DUX4-PAS"),
        ),
    )
    options = pipeline.Options(
        input=bam,
        ref_t2t=fasta,
        outdir=work / "out",
        sample="test",
        threads=4,
        tmpdir=work,
    )
    return pipeline.run(options, small)


def _table(path):
    with open(path, encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_every_allele_is_counted_from_spanning_reads(outdir):
    rows = _table(outdir / "test.alleles.tsv")
    found = [(r["chromosome"], r["reported_ru"], r["source"]) for r in rows]
    assert found == [
        (chromosome, str(count), "spanning")
        for chromosome, counts in TRUTH.items()
        for count in counts
    ]
    assert {r["haplotype"] for r in rows} == {"4qA", "10qA"}
    assert all(int(r["spanning_reads"]) == READS_PER_ALLELE for r in rows)


def test_each_chr4_allele_has_methylation(outdir):
    rows = _table(outdir / "test.alleles.tsv")
    for row in rows[:2]:
        assert 0 < float(row["DR1_methylation"]) < 100
        assert int(row["DR1_CpG"]) > 0
    units = _table(outdir / "test.unit_methylation.tsv")
    assert len(units) == sum(sum(counts) for counts in TRUTH.values())


def test_run_record_and_outputs(outdir):
    record = json.loads((outdir / "test.run.json").read_text())
    assert record["input_build"] == "t2t"
    assert record["reads"]["candidate_reads"] == 4 * READS_PER_ALLELE
    for name in (
        "reads_annotated.tsv",
        "repeats.tsv",
        "methylation.tsv",
        "read_alleles.tsv",
        "validation.txt",
        "d4z4.bam",
        "complete_reads.bam",
    ):
        assert (outdir / f"test.{name}").stat().st_size > 0
