"""Unit barcodes: the variants that tell repeat units apart.

Each read is cut into repeat units along the array (proximal to distal). Variant
sites are learned from the sample itself, and every unit is reduced to the
alleles it carries at those sites.

Array coordinates run from the proximal to the distal end, in the direction of
the repeat unit consensus. A read that matches the reverse strand of the
consensus is flipped.
"""

from collections.abc import Iterable, Mapping, Sequence
import dataclasses
import math
import re

import numpy as np

from fshdetector import _score
from fshdetector import blast

UNKNOWN = -1
DELETION = 4
_ALLELES = 5
_BASE_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
_BTOP = re.compile(r"\d+|[A-Z-]{2}")

MIN_FRAGMENT = 300  # Shorter hits carry too few sites to place a unit.
COUNTED_LENGTH = 1500  # A unit is counted if one hit is longer than this.
MIN_MINOR_COUNT = 4
MIN_MINOR_FRACTION = 0.01
# A real difference between units shows up on both strands; a systematic
# basecalling error mostly on one. Lower over higher minor allele frequency.
MIN_STRAND_RATIO = 0.25
HOMOPOLYMER = 3
# Alignments place the gaps of one difference in more than one way, so close
# sites rise and fall together. One site is kept per stretch, so that sites can
# be scored as independent observations.
MIN_SPACING = 10
_INDEL_WINDOW = 8
ERROR_PRIOR = 0.015
# Mismatches allowed above expectation: standard deviations plus a constant. A
# placement within the first limit is consistent; one beyond the second is
# clearly wrong; in between it is neither.
_CONSISTENT_LIMIT = (3.0, 1.0)
_PLAUSIBLE_LIMIT = (6.0, 2.0)
_ERROR_PRIOR_WEIGHT = 50.0
_ERROR_RANGE = (0.002, 0.25)
MAX_GROUP_DISAGREEMENT = 0.12


@dataclasses.dataclass
class Unit:
    """One repeat unit seen on one read.

    Attributes:
        bases: allele per consensus position (0-3 base, 4 deletion, -1 unseen).
        inserted: 1 if an insertion follows the position, 0 if not, -1 unseen.
        length: aligned length of the best hit.
        start: first array coordinate covered on the read (0-based).
        end: array coordinate after the last covered one.
    """

    bases: np.ndarray
    inserted: np.ndarray
    length: int
    start: int
    end: int


@dataclasses.dataclass
class ReadUnits:
    """The units of one read in array order; None marks a missed unit.

    Attributes:
        read_id: read name.
        reverse: the units match the reverse strand of the consensus.
        read_length: sequence length.
        units: units of the dominant orientation.
        mixed_strand: counted units occur in both orientations, as in a
            fold-back read; only the dominant orientation is kept.
    """

    read_id: str
    reverse: bool
    read_length: int
    units: list[Unit | None]
    mixed_strand: bool = False

    def span(self) -> tuple[int, int]:
        """The repeat array on the read, 0-based and half-open."""
        found = [u for u in self.units if u is not None]
        start = min(u.start for u in found)
        end = max(u.end for u in found)
        if self.reverse:
            return self.read_length - end, self.read_length - start
        return start, end

    def counted_range(self) -> tuple[int, int]:
        """Indices of the first and past the last counted unit."""
        found = [
            i
            for i, u in enumerate(self.units)
            if u is not None and u.length > COUNTED_LENGTH
        ]
        return (found[0], found[-1] + 1) if found else (0, 0)

    @property
    def counted(self) -> int:
        """Repeat units of the read.

        A unit counts if one of its hits is longer than `COUNTED_LENGTH`. A
        unit between two such units counts too: a poorly read stretch breaks
        its alignment into short pieces, or loses it, but the spacing of the
        neighbors shows that it is there.
        """
        start, end = self.counted_range()
        return end - start


def _observe(
    hit: blast.Hit, consensus: str, bases: np.ndarray, inserted: np.ndarray
) -> None:
    """Writes the alleles of one alignment where nothing is recorded yet."""
    position = hit.qstart - 1
    fresh = bases == UNKNOWN
    seen: dict[int, int] = {}
    insertions: set[int] = set()
    for token in _BTOP.findall(hit.btop):
        if token[0].isdigit():
            for p in range(position, position + int(token)):
                seen[p] = _BASE_CODE.get(consensus[p], UNKNOWN)
            position += int(token)
            continue
        query, subject = token
        if query == "-":
            insertions.add(position - 1)
            continue
        seen[position] = (
            DELETION if subject == "-" else _BASE_CODE.get(subject, UNKNOWN)
        )
        position += 1
    if not seen:
        return
    index = np.fromiter(seen.keys(), dtype=np.intp, count=len(seen))
    value = np.fromiter(seen.values(), dtype=np.int8, count=len(seen))
    keep = fresh[index]
    index = index[keep]
    bases[index] = value[keep]
    inserted[index] = 0
    added = [p for p in insertions if p >= 0 and p in seen and fresh[p]]
    inserted[added] = 1


def build_units(
    unit_hits: Iterable[blast.Hit],
    read_lengths: Mapping[str, int],
    consensus: str,
) -> dict[str, ReadUnits]:
    """Cuts every read into ordered repeat units.

    Hits whose inferred unit starts lie within half a period belong to one
    unit. The distance between neighboring starts gives the number of units in
    between, so a missed unit stays as a gap instead of closing up.

    Args:
        unit_hits: hits of the repeat unit query.
        read_lengths: read name to sequence length.
        consensus: the repeat unit query sequence.
    """
    period = len(consensus)
    grouped: dict[str, list[blast.Hit]] = {}
    for hit in unit_hits:
        if hit.length >= MIN_FRAGMENT:
            grouped.setdefault(hit.read_id, []).append(hit)

    result = {}
    for read_id, hits in grouped.items():
        length = read_lengths[read_id]
        reverse_bases = sum(h.length for h in hits if h.sstart > h.send)
        reverse = 2 * reverse_bases > sum(h.length for h in hits)
        clusters: list[tuple[int, list[tuple[blast.Hit, int, int]]]] = []
        for hit in sorted(hits, key=lambda h: -h.bitscore):
            if (hit.sstart > hit.send) != reverse:
                continue
            if reverse:
                start, end = length - hit.sstart, length - hit.send + 1
            else:
                start, end = hit.sstart - 1, hit.send
            anchor = start - (hit.qstart - 1)
            for existing, members in clusters:
                if abs(anchor - existing) < period / 2:
                    members.append((hit, start, end))
                    break
            else:
                clusters.append((anchor, [(hit, start, end)]))
        clusters.sort(key=lambda c: c[0])

        units: list[Unit | None] = []
        previous = None
        for anchor, members in clusters:
            if previous is not None:
                steps = max(1, math.floor((anchor - previous) / period + 0.5))
                units.extend([None] * (steps - 1))
            previous = anchor
            bases = np.full(period, UNKNOWN, dtype=np.int8)
            inserted = np.full(period, UNKNOWN, dtype=np.int8)
            for hit, _, _ in members:
                _observe(hit, consensus, bases, inserted)
            units.append(
                Unit(
                    bases=bases,
                    inserted=inserted,
                    length=max(m[0].length for m in members),
                    start=min(m[1] for m in members),
                    end=max(m[2] for m in members),
                )
            )
        mixed = any(
            (h.sstart > h.send) != reverse and h.length > COUNTED_LENGTH
            for h in hits
        )
        result[read_id] = ReadUnits(read_id, reverse, length, units, mixed)
    return result


@dataclasses.dataclass(frozen=True)
class Sites:
    """Variant sites of one sample.

    Attributes:
        position: consensus position (0-based) of each site.
        insertion: True where the site records an insertion after the position.
    """

    position: np.ndarray
    insertion: np.ndarray

    def __len__(self) -> int:
        return len(self.position)

    def subset(self, keep: np.ndarray) -> "Sites":
        """The sites selected by a boolean mask."""
        return Sites(self.position[keep], self.insertion[keep])


def homopolymer_mask(consensus: str) -> np.ndarray:
    """Marks homopolymer runs and one base on either side."""
    mask = np.zeros(len(consensus), dtype=bool)
    pattern = rf"(.)\1{{{HOMOPOLYMER - 1},}}"
    for run in re.finditer(pattern, consensus):
        mask[max(0, run.start() - 1) : run.end() + 1] = True
    return mask


def find_sites(reads: Iterable[ReadUnits], consensus: str) -> Sites:
    """Learns the sites where this sample's units differ.

    A site needs a minor allele with at least `MIN_MINOR_COUNT` observations
    and `MIN_MINOR_FRACTION` of all observations, with similar frequencies on
    the two strands (`MIN_STRAND_RATIO`), outside homopolymers. Sites keep
    `MIN_SPACING` between them; where candidates are closer, the one farthest
    from frequent insertions and deletions wins, then the one with the most
    minor allele observations.
    """
    period = len(consensus)
    counts = np.zeros((2, period, _ALLELES), dtype=np.int64)
    inserts = np.zeros((2, period, 2), dtype=np.int64)
    for read in reads:
        strand = int(read.reverse)
        for unit in read.units:
            if unit is None:
                continue
            index = np.flatnonzero(unit.bases >= 0)
            counts[strand, index, unit.bases[index]] += 1
            index = np.flatnonzero(unit.inserted >= 0)
            inserts[strand, index, unit.inserted[index]] += 1

    allowed = ~homopolymer_mask(consensus)
    gaps = counts[:, :, DELETION].sum(axis=0) + inserts[:, :, 1].sum(axis=0)
    indel = gaps / np.maximum(counts.sum(axis=(0, 2)), 1)
    candidates = []
    for table, is_insertion in ((counts, False), (inserts, True)):
        total = table.sum(axis=0)
        depth = total.sum(axis=1)
        minor = np.argsort(total, axis=1, kind="stable")[:, -2]
        rows = np.arange(period)
        minor_count = total[rows, minor]
        frequency = table[:, rows, minor] / np.maximum(table.sum(axis=2), 1)
        both = (frequency.min(axis=0) > 0) & (
            frequency.min(axis=0) >= MIN_STRAND_RATIO * frequency.max(axis=0)
        )
        selected = (
            allowed
            & both
            & (minor_count >= MIN_MINOR_COUNT)
            & (minor_count >= MIN_MINOR_FRACTION * np.maximum(depth, 1))
        )
        for p in np.flatnonzero(selected):
            nearby = indel[max(0, p - _INDEL_WINDOW) : p + _INDEL_WINDOW + 1]
            candidates.append(
                (
                    round(float(nearby.max()), 2),
                    -int(minor_count[p]),
                    int(p),
                    is_insertion,
                )
            )
    taken = np.zeros(period, dtype=bool)
    chosen = []
    for _, _, p, is_insertion in sorted(candidates):
        if not taken[p]:
            chosen.append((is_insertion, p))
            taken[max(0, p - MIN_SPACING + 1) : p + MIN_SPACING] = True
    chosen.sort()
    return Sites(
        np.array([p for _, p in chosen], dtype=np.intp),
        np.array([flag for flag, _ in chosen], dtype=bool),
    )


def matrix(read: ReadUnits, sites: Sites) -> np.ndarray:
    """Returns the units x sites allele matrix of a read (-1 unknown)."""
    result = np.full((len(read.units), len(sites)), UNKNOWN, dtype=np.int8)
    for row, unit in enumerate(read.units):
        if unit is None:
            continue
        result[row] = np.where(
            sites.insertion,
            unit.inserted[sites.position],
            unit.bases[sites.position],
        )
    return result


class AlleleCounts:
    """Allele counts per unit and site over placed reads."""

    def __init__(self, sites: int):
        self._counts = np.zeros((0, sites, _ALLELES), dtype=np.int32)

    def add(self, values: np.ndarray, offset: int) -> None:
        """Adds a read whose row j is unit j + offset."""
        needed = offset + values.shape[0] - self._counts.shape[0]
        if needed > 0:
            extra = np.zeros((needed, *self._counts.shape[1:]), dtype=np.int32)
            self._counts = np.concatenate([self._counts, extra])
        _score.add_counts(self._counts, values, offset)

    def consensus(self, rows: int) -> np.ndarray:
        """The most common allele of the first rows; ties are unknown."""
        return _score.consensus(self._counts, rows)


@dataclasses.dataclass(frozen=True)
class Placement:
    """Result of comparing a query with a template at one offset."""

    offset: int
    score: float
    sites: int
    mismatches: int
    consistent: bool


@dataclasses.dataclass(frozen=True)
class Scan:
    """Comparisons as arrays: one value per offset, or per pair of units.

    Attributes:
        score: log likelihood ratio.
        sites: compared sites.
        mismatches: sites with different alleles.
        informative: enough sites are compared.
        consistent: the mismatches are explained by read errors.
        plausible: the mismatches do not rule the placement out.
        units: units that share a site (offsets only).
    """

    score: np.ndarray
    sites: np.ndarray
    mismatches: np.ndarray
    informative: np.ndarray
    consistent: np.ndarray
    plausible: np.ndarray
    units: np.ndarray | None = None


class Model:
    """Scores unit comparisons with the error and allele rates of one sample.

    The score is a log likelihood ratio: both units are the same array unit,
    against two units drawn at random from this sample.

    Args:
        units: all unit rows of the sample (units x sites).
        error: per site probability that one read shows the wrong allele.
    """

    MIN_SITES = 3  # Fewer shared sites cannot support or reject a placement.

    def __init__(self, units: np.ndarray, error: np.ndarray | None = None):
        sites = units.shape[1]
        if error is None:
            error = np.full(sites, ERROR_PRIOR)
        self.error = np.ascontiguousarray(error, dtype=np.float64)
        counts = np.stack(
            [(units == a).sum(axis=0) for a in range(_ALLELES)], axis=1
        ).astype(np.float64)
        seen = counts > 0
        counts = np.where(seen, counts + 0.5, 0.0)
        totals = counts.sum(axis=1, keepdims=True)
        frequency = counts / np.maximum(totals, 1.0)
        kinds = np.maximum(seen.sum(axis=1), 2)
        e = self.error[:, None, None]
        # emission[s, truth, observed]
        emission = np.broadcast_to(
            e / (kinds[:, None, None] - 1), (sites, _ALLELES, _ALLELES)
        ).copy()
        diagonal = np.arange(_ALLELES)
        emission[:, diagonal, diagonal] = 1.0 - self.error[:, None]
        joint = np.einsum("st,stx,sty->sxy", frequency, emission, emission)
        marginal = np.einsum("st,stx->sx", frequency, emission)
        independent = marginal[:, :, None] * marginal[:, None, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.log(joint / independent)
        self.table = np.ascontiguousarray(np.nan_to_num(ratio, posinf=0.0))
        self.noise = np.ascontiguousarray(2 * self.error * (1 - self.error))

    def _judge(
        self,
        score: np.ndarray,
        sites: np.ndarray,
        mismatches: np.ndarray,
        expected: np.ndarray,
        variance: np.ndarray,
        units: np.ndarray | None = None,
    ) -> Scan:
        informative = sites >= self.MIN_SITES
        spread = np.sqrt(variance)

        def within(limit: tuple[float, float]) -> np.ndarray:
            return mismatches <= expected + limit[0] * spread + limit[1]

        consistent = informative & within(_CONSISTENT_LIMIT) & (score > 0)
        plausible = informative & within(_PLAUSIBLE_LIMIT)
        return Scan(
            score,
            sites,
            mismatches,
            informative,
            consistent,
            plausible,
            units,
        )

    def scan(
        self, template: np.ndarray, query: np.ndarray, offsets: Sequence[int]
    ) -> Scan:
        """Compares the query with the template at each offset."""
        return self._judge(
            *_score.score_offsets(
                np.ascontiguousarray(template),
                np.ascontiguousarray(query),
                self.table,
                self.noise,
                np.asarray(offsets, dtype=np.int_),
            )
        )

    def unit_scores(self, template: np.ndarray, query: np.ndarray) -> Scan:
        """Compares every template unit (rows) with every query unit."""
        return self._judge(
            *_score.score_units(
                np.ascontiguousarray(template),
                np.ascontiguousarray(query),
                self.table,
                self.noise,
            )
        )

    def place(
        self, template: np.ndarray, query: np.ndarray, offsets: Sequence[int]
    ) -> list[Placement]:
        """Like `scan`, as one object per offset."""
        scan = self.scan(template, query, offsets)
        return [
            Placement(
                int(offset),
                float(scan.score[i]),
                int(scan.sites[i]),
                int(scan.mismatches[i]),
                bool(scan.consistent[i]),
            )
            for i, offset in enumerate(offsets)
        ]

    @classmethod
    def informative(cls, placement: Placement) -> bool:
        return placement.sites >= cls.MIN_SITES

    @staticmethod
    def consistent(placement: Placement) -> bool:
        """True if the mismatches are explained by read errors."""
        return placement.consistent


def same_allele(matrices: Sequence[np.ndarray]) -> list[np.ndarray]:
    """Drops matrices that are far from the consensus of the group.

    All matrices start at the same unit. Reads of another array with the same
    unit count differ much more than read errors do.
    """
    if len(matrices) < 3:
        return list(matrices)
    rows = min(m.shape[0] for m in matrices)
    counts = AlleleCounts(matrices[0].shape[1])
    for m in matrices:
        counts.add(m[:rows], 0)
    consensus = counts.consensus(rows)
    kept = []
    for m in matrices:
        compared = (m[:rows] >= 0) & (consensus >= 0)
        differing = (m[:rows] != consensus) & compared
        if differing.sum() <= MAX_GROUP_DISAGREEMENT * compared.sum():
            kept.append(m)
    return kept


def estimate_error(groups: Iterable[Sequence[np.ndarray]]) -> np.ndarray | None:
    """Estimates per site error from reads that cover the same units.

    Args:
        groups: each group holds the matrices of reads that span one allele,
            all starting at the same unit.

    Returns:
        Error rate per site, shrunk toward `ERROR_PRIOR`, or None if no group
        has two reads.
    """
    mismatched = None
    compared = None
    for matrices in groups:
        if len(matrices) < 2:
            continue
        rows = min(m.shape[0] for m in matrices)
        stack = np.stack([m[:rows] for m in matrices])
        known = (stack >= 0).sum(axis=0)
        pairs = known * (known - 1) / 2
        same = np.zeros_like(pairs)
        for allele in range(_ALLELES):
            count = (stack == allele).sum(axis=0)
            same += count * (count - 1) / 2
        if mismatched is None:
            mismatched = np.zeros(stack.shape[2])
            compared = np.zeros(stack.shape[2])
        mismatched += (pairs - same).sum(axis=0)
        compared += pairs.sum(axis=0)
    if mismatched is None:
        return None
    # Two reads disagree about twice as often as one read is wrong.
    error = (mismatched / 2 + ERROR_PRIOR * _ERROR_PRIOR_WEIGHT) / (
        compared + _ERROR_PRIOR_WEIGHT
    )
    return np.clip(error, *_ERROR_RANGE)
