"""Genomic coordinates used by the analysis.

All intervals are 0-based and half-open, as in BED. T2T means T2T-CHM13v2.0.
"""

import dataclasses


@dataclasses.dataclass(frozen=True)
class Region:
    """A named genomic interval (0-based, half-open)."""

    chrom: str
    start: int
    end: int
    name: str = ""


@dataclasses.dataclass(frozen=True)
class Layout:
    """Coordinates and reference lengths. Tests use a smaller layout.

    Attributes:
        chr4_lengths: chr4 length of each recognized genome build.
        t2t_windows: D4Z4 loci on T2T where candidate reads are collected.
        hg38_windows: FSHD regions of hg38 whose reads are realigned to T2T.
        meth_targets: T2T intervals passed to modkit.
        meth_regions: T2T intervals that are summarized in the report.
    """

    chr4_lengths: dict[str, int]
    t2t_windows: tuple[Region, ...]
    hg38_windows: tuple[Region, ...]
    meth_targets: tuple[Region, ...]
    meth_regions: tuple[Region, ...]


CHROMOSOMES = ("chr4", "chr10")

HUMAN = Layout(
    chr4_lengths={
        "t2t": 193574945,
        "hg38": 190214555,
        "hg19": 191154276,
    },
    t2t_windows=(
        Region("chr4", 193230000, 193610000, "chr4_D4Z4"),
        Region("chr10", 134560000, 134800000, "chr10_D4Z4"),
    ),
    hg38_windows=(
        Region("chr4", 189000000, 190500000, "4q35_D4Z4_flanking"),
        Region("chr10", 133000000, 134000000, "10q26_D4Z4_flanking"),
    ),
    meth_targets=(
        Region("chr4", 193540172, 193543634, "chr4_distal"),
        Region("chr10", 134738000, 134741148, "chr10_distal"),
    ),
    meth_regions=(
        Region("chr4", 193541000, 193541600, "DR1"),
        Region("chr4", 193543430, 193543634, "DUX4-PAS"),
    ),
)
