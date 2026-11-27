"""Writes the result tables."""

from collections.abc import Iterable, Mapping, Sequence
import pathlib

from fshdetector import alleles as alleles_module
from fshdetector import annotate
from fshdetector import layout as layout_module
from fshdetector import methylation

NA = "NA"

# Healthy-control reference values (percent) shown next to each region.
# Atlas: ONT whole-genome data of about 45 to 49 controls.
ATLAS = {
    "DR1": ("54.7% (1.97%)", "[53.07%, 54.89%, 56.83%]"),
    "DUX4-PAS": ("51.81% (0.23%)", "[51.63%, 51.75%, 51.99%]"),
}
# Thresholds of Caputo et al. (2022), Cells 11:4114.
LITERATURE = {"DR1": "CpG1 <= 53%", "DUX4-PAS": "CpG6 <= 78%"}


def write_table(
    path: pathlib.Path,
    columns: Sequence[str],
    rows: Iterable[Sequence[object]],
) -> None:
    """Writes a tab-separated table; None becomes NA."""
    with open(path, "w", encoding="utf-8") as out:
        out.write("\t".join(columns) + "\n")
        for row in rows:
            out.write(
                "\t".join(NA if v is None else str(v) for v in row) + "\n"
            )


def write_reads(path: pathlib.Path, calls: Sequence[annotate.ReadCall]) -> None:
    columns = (
        "Read_ID", "Chromosome", "Chr_Confidence", "RU_Count", "Haplotype",
        "Length_Subtype", "Permissive", "Read_Length", "Has_Head", "Has_Tail",
        "Proximal_Marker", "Distal_Marker", "Chimeric", "BlnI_Chr4",
        "BlnI_Chr10", "Mixed_Strand",
    )  # fmt: skip
    write_table(
        path,
        columns,
        ([getattr(c, name.lower()) for name in columns] for c in calls),
    )


def write_repeats(
    path: pathlib.Path, calls: Sequence[annotate.ReadCall]
) -> None:
    """Largest unit count per chromosome, haplotype, and anchor pattern."""
    groups: dict[tuple, list[int]] = {}
    for c in calls:
        key = (c.chromosome, c.haplotype, c.has_head, c.has_tail)
        groups.setdefault(key, []).append(c.ru_count)
    order = {name: i for i, name in enumerate(layout_module.CHROMOSOMES)}
    rows = []
    for key in sorted(
        groups,
        key=lambda k: (order.get(k[0], 99), k[0], k[1], not k[2], not k[3]),
    ):
        chromosome, haplotype, head, tail = key
        largest = max(groups[key])
        count = str(largest) if head and tail else f">={largest}"
        rows.append(
            (chromosome, haplotype, count, head, tail, len(groups[key]))
        )
    write_table(
        path,
        ("Chromosome", "Haplotype", "RU_Count", "Head", "Tail", "Read_Count"),
        rows,
    )


def write_methylation(
    path: pathlib.Path,
    partitions: Mapping[str, Sequence[methylation.Site]],
    layout: layout_module.Layout,
) -> None:
    """DR1 and DUX4-PAS per upstream haplotype tag and for all reads."""
    rows = []
    for label, part in (
        ("HP1", "hp1"),
        ("HP2", "hp2"),
        ("Combined", "combined"),
    ):
        for region in layout.meth_regions:
            summary = methylation.summarize(partitions[part], region)
            atlas, spread = ATLAS[region.name]
            rows.append(
                (
                    label,
                    region.name,
                    summary.percent,
                    summary.n_cpg,
                    atlas,
                    spread,
                    LITERATURE[region.name],
                )
            )
    write_table(
        path,
        (
            "Haplotype", "Region", "Methylation", "N_CpG",
            "Ref_Atlas (Mean (SD) %)", "Ref_Atlas (IQR %)", "Lit_Caputo (%)",
        ),
        rows,
    )  # fmt: skip


def write_alleles(
    path: pathlib.Path,
    sample: str,
    alleles: Sequence[alleles_module.Allele],
    methylation_values: Sequence[Mapping[str, object]],
) -> None:
    """The final table: two rows per chromosome."""
    methylation_columns = (
        "DR1_methylation", "DR1_CpG", "DUX4-PAS_methylation",
        "DUX4-PAS_CpG", "distal_units_methylation", "distal_units_calls",
        "distal_units_reads",
    )  # fmt: skip
    columns = (
        "sample", "chromosome", "allele", "haplotype", "permissive",
        "reported_ru", "source", "spanning_ru", "spanning_reads",
        "reassembly_ru", "weakest_link_reads", "contig_lower_bound",
        "longest_read_ru", "supporting_reads", "conflicting_reads",
        *methylation_columns, "extra_alleles", "flags",
    )  # fmt: skip
    rows = []
    numbers: dict[str, int] = {}
    for allele, values in zip(alleles, methylation_values, strict=True):
        numbers[allele.chromosome] = numbers.get(allele.chromosome, 0) + 1
        rows.append(
            (
                sample,
                allele.chromosome,
                numbers[allele.chromosome],
                allele.haplotype or None,
                allele.permissive,
                allele.reported_ru,
                allele.source,
                allele.spanning_ru,
                allele.spanning_reads,
                allele.reassembly_ru,
                allele.weakest_link_reads,
                allele.contig_lower_bound,
                allele.longest_read_ru,
                allele.supporting_reads,
                allele.conflicting_reads,
                *(values.get(name) for name in methylation_columns),
                ",".join(allele.extra_alleles) or None,
                ",".join(allele.flags) or None,
            )
        )
    write_table(path, columns, rows)
