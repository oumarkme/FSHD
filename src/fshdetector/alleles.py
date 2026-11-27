"""Fills two allele slots per chromosome from tiered evidence.

1. ``spanning``: reads that hold both array ends. This is the primary result.
2. ``reassembly``: an exact count tiled from partial reads (experimental).
3. ``longest_read``: a lower bound from the longest read that the other allele
   cannot explain.

A lower tier never changes the count of a higher tier.
"""

from collections.abc import Mapping, Sequence
import collections
import dataclasses

from fshdetector import annotate
from fshdetector import assemble
from fshdetector import barcode
from fshdetector import layout as layout_module

SPANNING = "spanning"
REASSEMBLY = "reassembly"
LONGEST = "longest_read"
MISSING = "none"
SLOTS = 2


@dataclasses.dataclass
class Allele:
    """One reported allele slot; None means not available."""

    chromosome: str
    source: str = MISSING
    ru: int | None = None
    spanning_ru: int | None = None
    spanning_reads: int = 0
    reassembly_ru: int | None = None
    weakest_link_reads: int | None = None
    contig_lower_bound: int | None = None
    longest_read_ru: int | None = None
    supporting_reads: int = 0
    conflicting_reads: int = 0
    haplotype: str = ""
    permissive: str = "unknown"
    extra_alleles: list[str] = dataclasses.field(default_factory=list)
    flags: list[str] = dataclasses.field(default_factory=list)
    read_ids: list[str] = dataclasses.field(default_factory=list)
    template: assemble.Template | None = None

    @property
    def reported_ru(self) -> str:
        if self.ru is None:
            return "NA"
        return f">={self.ru}" if self.source == LONGEST else str(self.ru)


@dataclasses.dataclass
class _Group:
    reads: list[assemble.Fragment]
    ru: int
    flags: list[str] = dataclasses.field(default_factory=list)

    @property
    def kind(self) -> str:
        kinds = {f.kind for f in self.reads} - {""}
        return kinds.pop() if len(kinds) == 1 else ""

    def template(self) -> assemble.Template:
        rows = collections.Counter(f.units for f in self.reads)
        common = max(rows, key=lambda n: (rows[n], n))
        return assemble.Template(
            [(f, 0) for f in self.reads], common, self.kind
        )


def _clusters(
    model: barcode.Model, reads: Sequence[assemble.Fragment]
) -> list[list[assemble.Fragment]]:
    """Single-linkage groups of spanning reads with compatible barcodes."""
    parent = list(range(len(reads)))

    def root(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, first in enumerate(reads):
        for j in range(i + 1, len(reads)):
            placement = model.place(first.matrix, reads[j].matrix, [0])[0]
            if not model.informative(placement) or model.consistent(placement):
                parent[root(j)] = root(i)
    found: dict[int, list[assemble.Fragment]] = {}
    for i, read in enumerate(reads):
        found.setdefault(root(i), []).append(read)
    return list(found.values())


def group_spanning(
    model: barcode.Model,
    reads: Sequence[assemble.Fragment],
    counts: Mapping[str, int],
) -> list[_Group]:
    """Groups spanning reads of one chromosome into candidate alleles.

    Reads are grouped by unit count and A/B type, split where their barcodes
    form two supported clusters, and a single read that differs by one unit
    from a supported group with the same barcodes is merged into it.
    """
    by_count: dict[int, list[assemble.Fragment]] = {}
    for read in reads:
        by_count.setdefault(counts[read.read_id], []).append(read)

    groups = []
    for ru, members in sorted(by_count.items()):
        kinds = sorted({f.kind for f in members} - {""})
        if len(kinds) <= 1:
            subsets = [members]
        else:
            subsets = [
                [f for f in members if f.kind == kind] for kind in kinds + [""]
            ]
        for subset in filter(None, subsets):
            clusters = _clusters(model, subset)
            supported = [c for c in clusters if len(c) >= 2]
            if len(supported) >= 2:
                left_out = len(subset) - sum(map(len, supported))
                for cluster in supported:
                    group = _Group(cluster, ru, ["split_by_barcode"])
                    if left_out:
                        group.flags.append("unplaced_spanning_read")
                    groups.append(group)
            else:
                flags = ["possible_split"] if len(clusters) > 1 else []
                groups.append(_Group(list(subset), ru, flags))

    merged = []
    for group in sorted(groups, key=lambda g: len(g.reads)):
        if len(group.reads) == 1:
            single = group.reads[0]
            targets = [
                other
                for other in groups
                if other is not group
                and len(other.reads) >= 2
                and abs(other.ru - group.ru) == 1
                and assemble.kinds_agree(single.kind, other.kind)
            ]
            for target in targets:
                template = target.template()
                offsets = sorted({0, template.rows - single.units})
                if any(
                    model.consistent(p)
                    for p in model.place(
                        template.matrix, single.matrix, offsets
                    )
                ):
                    target.reads.append(single)
                    target.flags.append("count_discordance")
                    merged.append(group)
                    break
    return [g for g in groups if all(g is not m for m in merged)]


def _relations(
    model: barcode.Model,
    reads: Sequence[assemble.Fragment],
    alleles: Sequence[Allele],
) -> dict[str, list[str]]:
    """Relation of every partial read to every allele that has a template."""
    return {
        f.read_id: [
            assemble.compare_with_allele(model, f, a.template) for a in alleles
        ]
        for f in reads
    }


def _longest(
    reads: Sequence[assemble.Fragment], counts: Mapping[str, int]
) -> assemble.Fragment:
    return min(
        reads,
        key=lambda f: (
            -counts[f.read_id],
            -(f.head + f.tail),
            f.read_id,
        ),
    )


def _witnesses(
    model: barcode.Model,
    reads: Sequence[assemble.Fragment],
    everything: Sequence[assemble.Fragment],
    counts: Mapping[str, int],
) -> tuple[assemble.Fragment, bool]:
    """The longest read that other reads support, else the longest read.

    Poor and chimeric reads differ from every array, so on their own they must
    not stand for another allele.

    Returns:
        The read, and whether other reads support it.
    """
    remaining = list(reads)
    while remaining:
        read = _longest(remaining, counts)
        if not read.interrupted and assemble.corroborated(
            model, read, everything
        ):
            return read, True
        remaining = [f for f in remaining if f is not read]
    return _longest(reads, counts), False


def _describe(allele: Allele, calls: Mapping[str, annotate.ReadCall]) -> None:
    """Sets haplotype and permissiveness from the assigned reads."""
    prefix = "4q" if allele.chromosome == "chr4" else "10q"
    typed = [
        calls[r]
        for r in allele.read_ids
        if calls[r].haplotype in (prefix + "A", prefix + "B")
    ]
    if not typed:
        allele.haplotype = prefix if allele.ru is not None else ""
        return
    votes = collections.Counter(c.haplotype for c in typed)
    top = votes.most_common()
    if len(top) > 1:
        allele.flags.append("mixed_haplotype")
    if len(top) > 1 and top[0][1] == top[1][1]:
        allele.haplotype = prefix
        return
    allele.haplotype = top[0][0]
    permissive = sum(c.permissive for c in typed)
    allele.permissive = "yes" if 2 * permissive > len(typed) else "no"


def call_chromosome(
    chromosome: str,
    model: barcode.Model,
    calls: Mapping[str, annotate.ReadCall],
    spanning: Sequence[assemble.Fragment],
    partial: Sequence[assemble.Fragment],
) -> list[Allele]:
    """Fills the two slots of one chromosome.

    Args:
        chromosome: "chr4" or "chr10".
        model: scoring model of the sample.
        calls: read name to read call.
        spanning: reads of this chromosome that hold both array ends.
        partial: all other reads of this chromosome with repeat units.
    """
    counts = {r: c.ru_count for r, c in calls.items()}
    groups = group_spanning(model, spanning, counts)
    chimeras = [
        g
        for g in groups
        if len(g.reads) == 1
        and assemble.suspect_chimera(model, g.reads[0], [*spanning, *partial])
    ]
    groups = [g for g in groups if all(g is not c for c in chimeras)]
    groups.sort(key=lambda g: (-len(g.reads), g.ru))
    extra = [f"{g.ru}x{len(g.reads)}" for g in groups[SLOTS:]]
    alleles = []
    for group in groups[:SLOTS]:
        allele = Allele(
            chromosome,
            SPANNING,
            ru=group.ru,
            spanning_ru=group.ru,
            spanning_reads=len(group.reads),
            flags=sorted(set(group.flags)),
            template=group.template(),
        )
        if len(group.reads) == 1:
            allele.flags.append("single_read")
        if extra:
            allele.flags.append("extra_allele")
            allele.extra_alleles = extra
        allele.read_ids = [f.read_id for f in group.reads]
        alleles.append(allele)

    contigs: list[assemble.Template] = []
    if len(alleles) < SLOTS:
        relations = _relations(model, partial, alleles)
        supporting = [f for f in partial if "fits" in relations[f.read_id]]
        pool = [f for f in partial if "fits" not in relations[f.read_id]]
        called_arrays = [a.template for a in alleles]
        assemblies, contigs = assemble.reassemble(
            model, pool, called=called_arrays, chromosome=chromosome
        )
        # Reads that fit a called allele may also belong to the other one.
        shared = [
            f
            for f in supporting
            if any(
                assemble.compare_with_allele(model, f, c) == "fits"
                for c in contigs
            )
        ]
        if shared:
            assemblies, contigs = assemble.reassemble(
                model,
                pool + shared,
                frozenset(f.read_id for f in shared),
                called_arrays,
                chromosome,
            )
        for result in assemblies:
            if len(alleles) >= SLOTS or not result.exact:
                continue
            alleles.append(
                Allele(
                    chromosome,
                    REASSEMBLY,
                    ru=result.template.rows,
                    reassembly_ru=result.template.rows,
                    weakest_link_reads=result.weakest_link,
                    flags=["experimental"],
                    template=result.template,
                    read_ids=list(result.template.read_ids),
                )
            )

    called = list(alleles)
    relations = _relations(model, partial, called)
    by_id = {f.read_id: f for f in partial}
    lower: list[tuple[Allele, assemble.Fragment]] = []
    if len(alleles) < SLOTS:
        if called:
            open_reads = [
                f
                for f in partial
                if relations[f.read_id]
                and all(
                    r in ("differs", "overruns") for r in relations[f.read_id]
                )
            ]
        else:
            open_reads = list(partial)
        while open_reads and len(alleles) < SLOTS:
            read, supported = _witnesses(
                model, open_reads, [*spanning, *partial], counts
            )
            if not supported and (called or lower):
                break
            allele = Allele(
                chromosome,
                LONGEST,
                ru=counts[read.read_id],
                read_ids=[read.read_id],
                flags=[] if supported else ["unsupported_read"],
            )
            alleles.append(allele)
            lower.append((allele, read))
            open_reads = [
                f
                for f in open_reads
                if f is not read and assemble.reads_incompatible(model, read, f)
            ]

    # Assign partial reads that point to exactly one allele.
    anchors = {read.read_id for _, read in lower}
    explained: dict[str, list[Allele]] = {}
    for f in partial:
        homes = [
            a
            for a, relation in zip(called, relations[f.read_id], strict=True)
            if relation == "fits"
        ]
        if not homes:
            homes = [
                a
                for a, read in lower
                if read.read_id == f.read_id
                or not assemble.reads_incompatible(model, read, f)
            ]
        explained[f.read_id] = homes
        if len(homes) == 1 and f.read_id not in anchors:
            home = homes[0]
            if f.read_id not in home.read_ids:
                home.read_ids.append(f.read_id)

    unexplained = sum(
        not explained[f.read_id]
        and bool(relations[f.read_id])
        and not {"uninformative", "unclear"} & set(relations[f.read_id])
        for f in partial
    )
    for index, allele in enumerate(called):
        overrun = [
            r for r, row in relations.items() if row[index] == "overruns"
        ]
        allele.conflicting_reads = len(overrun)
        if any(not explained[r] for r in overrun):
            allele.flags.append("assembly_conflict")
    for allele in alleles:
        allele.supporting_reads = sum(r in by_id for r in allele.read_ids)
        units = [counts[r] for r in allele.read_ids]
        allele.longest_read_ru = max(units, default=None)
        bounds = [
            c.supported_rows()
            for c in contigs
            if set(c.read_ids) & set(allele.read_ids)
        ]
        allele.contig_lower_bound = max(bounds, default=None)
        if allele.source == LONGEST and allele.supporting_reads == 1:
            allele.flags.append("single_read")
        if unexplained >= 2:
            allele.flags.append("unexplained_reads")
        if chimeras:
            allele.flags.append("suspected_chimera_excluded")

    reads_present = bool(spanning) or bool(partial)
    while len(alleles) < SLOTS:
        flag = "second_allele_not_seen" if reads_present else "no_reads"
        alleles.append(Allele(chromosome, flags=[flag]))
    for allele in alleles:
        _describe(allele, calls)
    alleles.sort(key=lambda a: (a.ru is None, a.ru or 0))
    return alleles


def call_alleles(
    model: barcode.Model,
    calls: Sequence[annotate.ReadCall],
    fragments: Mapping[str, assemble.Fragment],
) -> list[Allele]:
    """Returns two alleles for chr4 followed by two for chr10."""
    by_id = {c.read_id: c for c in calls}
    result = []
    for chromosome in layout_module.CHROMOSOMES:
        own = [
            fragments[c.read_id]
            for c in calls
            if c.chromosome == chromosome and c.read_id in fragments
        ]
        spanning = [f for f in own if f.head and f.tail]
        # Flanking markers of the other chromosome make a partial read useless
        # as evidence for this one.
        partial = [
            f
            for f in own
            if not (f.head and f.tail) and f.hint in ("", chromosome)
        ]
        result += call_chromosome(chromosome, model, by_id, spanning, partial)
    return result
