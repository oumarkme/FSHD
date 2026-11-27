"""Per-read chromosome, haplotype, and repeat unit calls.

Chromosome evidence, from weakest to strongest:

1. Alignment of the primary record.
2. Restriction motif voting (BlnI-site variants and XapI) inside the repeat
   array of the read. XapI sites are common outside the array, so the rest of
   the read is not counted.
3. Haplotype markers: concordant proximal and distal markers override the
   other evidence; a single marker decides only when no motif was found.

A read is typed 4qA or 4qB (10qA or 10qB) only from a confident distal marker;
every other read stays 4q or 10q.
"""

from collections.abc import Iterable, Mapping
import dataclasses
import re

from fshdetector import barcode
from fshdetector import blast

PROXIMAL_MARKERS = frozenset({"4qA_D4F104S1", "4qB_D4F104S1", "10qA_D4F104S1"})
DISTAL_MARKERS = frozenset({"4qA_pLAM", "4qB_pLAM", "10qA_pLAM"})
MARKER_CHROMOSOME = {
    "4qA_D4F104S1": "chr4",
    "4qB_D4F104S1": "chr4",
    "10qA_D4F104S1": "chr10",
    "4qA_pLAM": "chr4",
    "4qB_pLAM": "chr4",
    "10qA_pLAM": "chr10",
}

LONG_UNIT_LENGTH = 2800  # Terminal unit longer than this is subtype L.
MIN_MARKER_IDENTITY = 85.0
MIN_P13E11_LENGTH = 200
CONFIDENT_IDENTITY = 90.0

_CHR4_BLNI = ("CTAGGTCTAGG", "CCTAGACCTAG")
_CHR10_BLNI = ("CTAGGCCTAGG", "CCTAGGCCTAG")
_XAPI = re.compile(r"[AG]AATT[CT]")


@dataclasses.dataclass
class ReadCall:
    """Calls for one read that carries at least one repeat unit."""

    read_id: str
    chromosome: str
    chr_confidence: str
    ru_count: int
    haplotype: str
    length_subtype: str
    permissive: bool
    read_length: int
    has_head: bool
    has_tail: bool
    proximal_marker: str
    distal_marker: str
    chimeric: bool
    blni_chr4: int
    blni_chr10: int
    mixed_strand: bool = False

    @property
    def complete(self) -> bool:
        """True if the read holds both ends of the array."""
        return self.has_head and self.has_tail


@dataclasses.dataclass(frozen=True)
class Motifs:
    """Restriction motif counts of one read."""

    chr4: int
    chr10: int
    xapi: int

    @property
    def total(self) -> int:
        return self.chr4 + self.chr10 + self.xapi


def count_motifs(sequence: str) -> Motifs:
    """Counts chromosome-specific motifs on both strands."""
    return Motifs(
        chr4=sum(sequence.count(motif) for motif in _CHR4_BLNI),
        chr10=sum(sequence.count(motif) for motif in _CHR10_BLNI),
        xapi=len(_XAPI.findall(sequence)),
    )


def motif_confidence(motifs: Motifs) -> str:
    """Grades how one-sided the motif evidence is."""
    if motifs.total == 0:
        return "LOW"
    chr4 = motifs.chr4 + motifs.xapi
    chr10 = motifs.chr10
    if min(chr4, chr10) == 0:
        return "HIGH"
    ratio = max(chr4, chr10) / motifs.total
    if ratio >= 0.8:
        return "HIGH"
    return "MEDIUM" if ratio >= 0.6 else "LOW"


def marker_chromosome(proximal: str, distal: str) -> str:
    """The chromosome that the markers of a read agree on, or ""."""
    found = {MARKER_CHROMOSOME[m] for m in (proximal, distal) if m}
    return found.pop() if len(found) == 1 else ""


def assign_chromosome(
    aligned: str, motifs: Motifs, proximal: str, distal: str
) -> str:
    """Combines alignment, motif, and marker evidence."""
    chromosome = aligned
    if motifs.chr10 > motifs.chr4 + motifs.xapi:
        chromosome = "chr10"
    elif motifs.chr4 + motifs.xapi > motifs.chr10:
        chromosome = "chr4"
    proximal_chromosome = MARKER_CHROMOSOME.get(proximal)
    distal_chromosome = MARKER_CHROMOSOME.get(distal)
    if proximal_chromosome and distal_chromosome:
        if proximal_chromosome == distal_chromosome:
            return proximal_chromosome
        return chromosome
    if motifs.total == 0:
        return proximal_chromosome or distal_chromosome or chromosome
    return chromosome


def is_chimeric(proximal: str, distal: str) -> bool:
    """True if the two markers come from different haplotype types."""
    if not proximal or not distal:
        return False
    return proximal.split("_")[0] != distal.split("_")[0]


def distal_confident(
    chromosome: str, distal: str, identity: float, chimeric: bool
) -> bool:
    """True if the distal marker is strong enough to type the read.

    The marker must match the assigned chromosome. The one exception is the
    4qA pLAM sequence on a chr10 read, a known translocation.
    """
    if chimeric or not distal or identity < CONFIDENT_IDENTITY:
        return False
    if MARKER_CHROMOSOME[distal] == chromosome:
        return True
    return chromosome == "chr10" and distal == "4qA_pLAM"


def classify_haplotype(chromosome: str, distal: str, confident: bool) -> str:
    """Returns 4qA, 4qB, 4q, 10qA, 10qB, 10q, or Unknown."""
    prefix = {"chr4": "4q", "chr10": "10q"}.get(chromosome)
    if prefix is None:
        return "Unknown"
    if not confident:
        return prefix
    return prefix + ("A" if distal in ("4qA_pLAM", "10qA_pLAM") else "B")


def _best_by_score(hits: Iterable[blast.Hit]) -> blast.Hit | None:
    return max(hits, key=lambda h: h.bitscore, default=None)


def _length_subtype(units: list[blast.Hit], plam_hits: list[blast.Hit]) -> str:
    """L or S from the unit next to the 4qA pLAM hit."""
    if not units or not plam_hits:
        return ""
    plam = min(plam_hits, key=lambda h: h.evalue)
    position = min(plam.sstart, plam.send)
    closest = min(units, key=lambda h: abs(min(h.sstart, h.send) - position))
    return "L" if closest.length > LONG_UNIT_LENGTH else "S"


def _group(hits: Iterable[blast.Hit]) -> dict[str, list[blast.Hit]]:
    grouped: dict[str, list[blast.Hit]] = {}
    for hit in hits:
        grouped.setdefault(hit.read_id, []).append(hit)
    return grouped


def annotate(
    read_units: Mapping[str, barcode.ReadUnits],
    unit_hits: Iterable[blast.Hit],
    marker_hits: Iterable[blast.Hit],
    p13e11_hits: Iterable[blast.Hit],
    sequences: Mapping[str, str],
    aligned_chromosomes: Mapping[str, str],
) -> list[ReadCall]:
    """Calls every read that has at least one counted repeat unit.

    Args:
        read_units: repeat units of each read; the unit count is the number of
            distinct unit positions, so a unit split into two alignments is
            counted once.
        unit_hits: hits of the repeat unit query.
        marker_hits: hits of the haplotype marker queries.
        p13e11_hits: hits of the p13E-11 query.
        sequences: read name to full sequence.
        aligned_chromosomes: contig of each mapped primary record.

    Returns:
        Calls sorted by decreasing unit count, then read name.
    """
    long_hits = _group(
        h for h in unit_hits if h.length > barcode.COUNTED_LENGTH
    )
    markers = _group(
        h for h in marker_hits if h.identity >= MIN_MARKER_IDENTITY
    )
    p13e11 = {
        h.read_id
        for h in p13e11_hits
        if h.identity >= MIN_MARKER_IDENTITY and h.length >= MIN_P13E11_LENGTH
    }
    calls = []
    for read_id, read in read_units.items():
        if not read.counted:
            continue
        hits = markers.get(read_id, [])
        proximal_hit = _best_by_score(
            h for h in hits if h.query in PROXIMAL_MARKERS
        )
        distal_hit = _best_by_score(
            h for h in hits if h.query in DISTAL_MARKERS
        )
        proximal = proximal_hit.query if proximal_hit else ""
        distal = distal_hit.query if distal_hit else ""
        sequence = sequences[read_id]
        start, end = read.span()
        motifs = count_motifs(sequence[start:end])
        chromosome = assign_chromosome(
            aligned_chromosomes.get(read_id, "Unknown"),
            motifs,
            proximal,
            distal,
        )
        chimeric = is_chimeric(proximal, distal)
        confident = distal_confident(
            chromosome,
            distal,
            distal_hit.identity if distal_hit else 0.0,
            chimeric,
        )
        haplotype = classify_haplotype(chromosome, distal, confident)
        subtype = ""
        if haplotype in ("4qA", "10qA") and distal == "4qA_pLAM":
            subtype = _length_subtype(
                long_hits[read_id],
                [h for h in hits if h.query == "4qA_pLAM"],
            )
        calls.append(
            ReadCall(
                read_id=read_id,
                chromosome=chromosome,
                chr_confidence=motif_confidence(motifs),
                ru_count=read.counted,
                haplotype=haplotype,
                length_subtype=subtype,
                permissive=confident and distal == "4qA_pLAM",
                read_length=len(sequence),
                has_head=bool(proximal) or read_id in p13e11,
                has_tail=bool(distal),
                proximal_marker=proximal,
                distal_marker=distal,
                chimeric=chimeric,
                blni_chr4=motifs.chr4,
                blni_chr10=motifs.chr10,
                mixed_strand=read.mixed_strand,
            )
        )
    return sorted(calls, key=lambda c: (-c.ru_count, c.read_id))
