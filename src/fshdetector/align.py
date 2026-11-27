"""Realigns hg38 reads to T2T-CHM13v2.0 with minimap2."""

from collections.abc import Mapping
import pathlib

import pysam

from fshdetector import tools

PRESET = "lr:hq"
_CARRIED_TAGS = ("MM", "ML", "MN", "HP", "PS")


def build_index(
    reference: str | pathlib.Path, output: str | pathlib.Path, threads: int
) -> None:
    """Writes a minimap2 index of the reference."""
    tools.require("minimap2")
    tools.run(
        ["minimap2", "-x", PRESET, "-t", threads, "-d", output, reference]
    )


def _tag_text(name: str, value, kind: str) -> str:
    if name == "ML":
        return "ML:B:C," + ",".join(map(str, value))
    return f"{name}:{'i' if kind in 'cCsSiI' else kind}:{value}"


def write_fastq(
    path: pathlib.Path, records: Mapping[str, pysam.AlignedSegment]
) -> None:
    """Writes reads in sequencing orientation, with tags in the comment."""
    with open(path, "w", encoding="ascii") as out:
        for name, record in records.items():
            sequence = record.get_forward_sequence()
            qualities = record.get_forward_qualities()
            quality = (
                "".join(chr(q + 33) for q in qualities)
                if qualities is not None
                else "!" * len(sequence)
            )
            tags = [
                _tag_text(tag, *record.get_tag(tag, with_value_type=True))
                for tag in _CARRIED_TAGS
                if record.has_tag(tag)
            ]
            comment = "\t".join(tags)
            out.write(f"@{name} {comment}\n{sequence}\n+\n{quality}\n")


def realign(
    records: Mapping[str, pysam.AlignedSegment],
    target: str | pathlib.Path,
    workdir: pathlib.Path,
    threads: int,
) -> pathlib.Path:
    """Aligns reads to T2T and returns a sorted, indexed BAM.

    Args:
        records: read name to the record that holds the full sequence.
        target: T2T FASTA or its minimap2 index.
        workdir: scratch directory.
        threads: minimap2 and sort threads.
    """
    fastq = workdir / "candidates.fastq"
    write_fastq(fastq, records)
    sam = workdir / "realigned.sam"
    # -Y keeps full sequences on supplementary records.
    tools.run(
        [
            "minimap2",
            "-ax", PRESET,
            "-y", "-Y", "--eqx", "--secondary=no",
            "-t", threads,
            target,
            fastq,
        ],
        stdout_path=str(sam),
    )  # fmt: skip
    bam = workdir / "realigned.bam"
    pysam.sort("-@", str(threads), "-o", str(bam), str(sam))
    pysam.index(str(bam))
    return bam
