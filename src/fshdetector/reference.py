"""Genome build detection and reference preparation."""

from collections.abc import Mapping
import pathlib

import pysam

from fshdetector import layout as layout_module


class ReferenceError_(ValueError):
    """The input or reference does not match a supported genome build."""


def detect_build(
    lengths: Mapping[str, int], layout: layout_module.Layout
) -> str:
    """Names the genome build of an alignment from its chromosome lengths.

    Args:
        lengths: contig name to length, from the alignment header.
        layout: coordinates that list the chr4 length of each build.

    Returns:
        "t2t" or "hg38".

    Raises:
        ReferenceError_: the build is hg19, unknown, or uses other contig names.
    """
    missing = [c for c in layout_module.CHROMOSOMES if c not in lengths]
    if missing:
        raise ReferenceError_(
            f"the input header has no contig named {', '.join(missing)}; "
            "only 'chr'-prefixed hg38 and T2T-CHM13v2.0 inputs are supported"
        )
    chr4 = lengths["chr4"]
    for build, length in layout.chr4_lengths.items():
        if chr4 == length:
            if build == "hg19":
                raise ReferenceError_(
                    "the input is aligned to hg19, which is not supported; "
                    "realign it to hg38 or T2T-CHM13v2.0"
                )
            return build
    raise ReferenceError_(
        f"chr4 length {chr4} matches neither hg38 nor T2T-CHM13v2.0; "
        "set --build if you know the coordinates are compatible"
    )


def check_t2t_reference(
    path: str | pathlib.Path, layout: layout_module.Layout
) -> None:
    """Fails unless the FASTA is indexed and has T2T chr4 and chr10."""
    index = pathlib.Path(f"{path}.fai")
    if not index.is_file():
        raise ReferenceError_(
            f"{index} not found; index the reference with 'samtools faidx'"
        )
    with pysam.FastaFile(str(path)) as fasta:
        for chrom in layout_module.CHROMOSOMES:
            if chrom not in fasta.references:
                raise ReferenceError_(f"{path} has no contig named {chrom}")
        length = fasta.get_reference_length("chr4")
    if length != layout.chr4_lengths["t2t"]:
        raise ReferenceError_(
            f"{path} has chr4 length {length}, not T2T-CHM13v2.0 "
            f"({layout.chr4_lengths['t2t']})"
        )


def write_chromosome_reference(
    path: str | pathlib.Path, destination: pathlib.Path
) -> None:
    """Writes chr4 and chr10 of a reference to a small indexed FASTA.

    modkit reads the whole reference it is given, so a two-chromosome copy
    saves minutes per run.
    """
    width = 100_000
    with (
        pysam.FastaFile(str(path)) as fasta,
        open(destination, "w", encoding="ascii") as out,
    ):
        for chrom in layout_module.CHROMOSOMES:
            out.write(f">{chrom}\n")
            length = fasta.get_reference_length(chrom)
            for start in range(0, length, width):
                out.write(fasta.fetch(chrom, start, start + width) + "\n")
    pysam.faidx(str(destination))


def target_sequences(
    path: str | pathlib.Path, layout: layout_module.Layout
) -> list[str]:
    """Returns the reference sequence of each methylation target."""
    with pysam.FastaFile(str(path)) as fasta:
        return [
            fasta.fetch(region.chrom, region.start, region.end)
            for region in layout.meth_targets
        ]
