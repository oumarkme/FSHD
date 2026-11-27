"""Supporting evidence: tiles partial reads by their unit barcodes.

A read matrix is trimmed so that row 0 of a proximally anchored read is the
first unit of the array, and the last row of a distally anchored read is the
last unit.

Repeat units resemble each other, so a wrong tiling can score well. The module
therefore prefers no count to a doubtful one:

1. A read extends a contig only if it has one consistent placement over all
   contigs, shares enough units there, and agrees with the other reads that
   extend the same contig.
2. A proximal and a distal contig are joined only if one array length is
   consistent with both and every other length is contradicted.
3. A joined array is exact only if each link between neighboring units is
   covered by several reads, no read contradicts it, and random subsets of the
   reads give the same count and no other.
4. A joined array is not exact if its count equals that of an allele that
   spanning reads already show.
"""

from collections.abc import Sequence
import dataclasses

import numpy as np

from fshdetector import barcode

EXTEND_SCORE = 4.6  # ln(100): odds needed to add a read to a contig.
EXTEND_MARGIN = 2.3  # ln(10): lead over the second best placement.
# Units a read shares with a contig that it extends. A short read needs all but
# one of its units there, and every read at least two.
MIN_OVERLAP = 4
BRANCH_READS = 1  # Disagreeing reads that stop the extension of a contig.
JOIN_SCORE = 6.9  # ln(1000): floor of the join threshold.
JOIN_MARGIN = 3.0  # ln(20): floor of the join margin.
MIN_JOIN_OVERLAP = 2
_DECOYS = 300
_DECOY_QUANTILE = 0.99
_SEED = 20261126
MIN_LINK_READS = 2
MIN_BOUND_LINK_READS = 3  # Links of a contig that counts as a lower bound.
MIN_ALLELE_READS = 3
MAX_DIVERGENT_READS = 3  # Unexplained reads that share part of the array.
MIN_DIVERGING_READS = 2  # Reads that contradict a junction of a single read.
MAX_CROSSING_FRACTION = 0.1  # Reads across a junction per read against it.
MIN_CORROBORATING_READS = 2  # Support wanted for a lower bound from one read.
RESAMPLES = 20
RESAMPLE_FRACTION = 0.8
MIN_AGREEMENT = 0.3  # Fraction of read subsets that give the same count.
MAX_OTHER_COUNTS = 0  # Read subsets that may give a different count.
_SAME_UNIT_SITES = 1  # Consensus units that differ in at most this are alike.


@dataclasses.dataclass
class Fragment:
    """The part of a read that the assembly uses.

    Attributes:
        read_id: read name.
        matrix: units x sites alleles.
        head: the proximal end of the array is on the read.
        tail: the distal end of the array is on the read.
        kind: "A" or "B" from a confident distal marker, else "".
        hint: chromosome of the flanking markers, or "".
        interrupted: a unit before the last one is missing or much shorter
            than a repeat unit.
        units: repeat units of the read.
    """

    read_id: str
    matrix: np.ndarray
    head: bool
    tail: bool
    kind: str = ""
    hint: str = ""
    interrupted: bool = False

    @property
    def units(self) -> int:
        return self.matrix.shape[0]


@dataclasses.dataclass
class Template:
    """A called allele or a contig, as placed reads.

    Attributes:
        placed: (fragment, offset) pairs.
        rows: number of units.
        kind: "A", "B", or "".
    """

    placed: list[tuple[Fragment, int]]
    rows: int
    kind: str = ""
    matrix: np.ndarray = dataclasses.field(init=False)

    def __post_init__(self):
        self._counts = barcode.AlleleCounts(self.placed[0][0].matrix.shape[1])
        for fragment, offset in self.placed:
            self._counts.add(fragment.matrix, offset)
        self.matrix = self._counts.consensus(self.rows)
        if not self.kind:
            kinds = {f.kind for f, _ in self.placed} - {""}
            self.kind = kinds.pop() if len(kinds) == 1 else ""

    def add(self, fragment: Fragment, offset: int) -> None:
        """Places one more read and updates the consensus."""
        self.placed.append((fragment, offset))
        self.rows = max(self.rows, offset + fragment.units)
        self._counts.add(fragment.matrix, offset)
        self.matrix = self._counts.consensus(self.rows)
        self.kind = self.kind or fragment.kind

    def copy(self) -> "Template":
        return Template(list(self.placed), self.rows, self.kind)

    @property
    def read_ids(self) -> list[str]:
        return [f.read_id for f, _ in self.placed]

    def link_reads(self) -> list[list[str]]:
        """Reads that cover both units of each neighboring pair."""
        return [
            [
                f.read_id
                for f, offset in self.placed
                if offset <= row and offset + f.units > row + 1
            ]
            for row in range(self.rows - 1)
        ]

    def supported_rows(self) -> int:
        """Units from the first one on that are linked by several reads."""
        rows = 1
        for covering in self.link_reads():
            if len(covering) < MIN_BOUND_LINK_READS:
                break
            rows += 1
        return rows


def kinds_agree(first: str, second: str) -> bool:
    return not first or not second or first == second


def compare_with_allele(
    model: barcode.Model, fragment: Fragment, allele: Template
) -> str:
    """Relates a partial read to a complete allele.

    Returns:
        "fits" if the read can lie inside the allele, "unclear" if that cannot
        be ruled out, "overruns" if it matches the allele but has more units,
        "uninformative" if too few sites are shared, and "differs" otherwise.
    """
    if not kinds_agree(fragment.kind, allele.kind):
        return "differs"
    n, m = allele.rows, fragment.units
    if fragment.head:
        offsets = np.array([0])
    elif fragment.tail:
        offsets = np.array([n - m])
    else:
        offsets = np.arange(1 - m, n)
    scan = model.scan(allele.matrix, fragment.matrix, offsets)
    if not scan.informative.any():
        return "uninformative"
    inside = (offsets >= 0) & (offsets + m <= n)
    if (scan.consistent & inside).any():
        return "fits"
    if (scan.plausible & inside).any():
        return "unclear"
    return "overruns" if scan.consistent.any() else "differs"


def reads_incompatible(
    model: barcode.Model, first: Fragment, second: Fragment
) -> bool:
    """True if two partial reads cannot come from the same array.

    Reads that do not share an anchored end can come from separate stretches of
    one array, so only different A/B types or a clear mismatch at a shared
    anchor prove two alleles.
    """
    if not kinds_agree(first.kind, second.kind):
        return True
    for shared, offset in (
        (first.head and second.head, 0),
        (first.tail and second.tail, first.units - second.units),
    ):
        if shared:
            scan = model.scan(first.matrix, second.matrix, [offset])
            if scan.informative[0] and not scan.plausible[0]:
                return True
    return False


def _flip(fragment: Fragment) -> Fragment:
    return Fragment(
        fragment.read_id,
        np.ascontiguousarray(fragment.matrix[::-1]),
        fragment.tail,
        fragment.head,
        fragment.kind,
        fragment.hint,
        fragment.interrupted,
    )


def _seed_contigs(
    model: barcode.Model, fragments: Sequence[Fragment]
) -> list[Template]:
    """Groups the reads that start at the anchored end, longest first."""
    anchored = sorted(
        (f for f in fragments if f.head and not f.tail),
        key=lambda f: (-f.units, f.read_id),
    )
    contigs: list[Template] = []
    for fragment in anchored:
        options = []
        informative = True
        for contig in contigs:
            if not kinds_agree(fragment.kind, contig.kind):
                continue
            placement = model.place(contig.matrix, fragment.matrix, [0])[0]
            informative &= model.informative(placement)
            if model.consistent(placement):
                options.append((placement.score, contig))
        options.sort(key=lambda o: -o[0])
        if len(options) == 1 or (
            len(options) > 1 and options[0][0] - options[1][0] >= EXTEND_MARGIN
        ):
            options[0][1].add(fragment, 0)
        elif not options and informative:
            contigs.append(
                Template([(fragment, 0)], fragment.units, fragment.kind)
            )
    return contigs


def _unique_placement(
    model: barcode.Model, fragment: Fragment, contigs: Sequence[Template]
) -> tuple[float, int, int] | None:
    """The only consistent placement of a read, as (score, contig, offset)."""
    found = None
    runner_up = 0.0
    for index, contig in enumerate(contigs):
        if not kinds_agree(fragment.kind, contig.kind):
            continue
        offsets = np.arange(contig.rows)
        scan = model.scan(contig.matrix, fragment.matrix, offsets)
        hits = np.flatnonzero(scan.consistent)
        if len(hits) > 1 or (len(hits) and found is not None):
            return None
        others = scan.score[~scan.consistent]
        if len(others):
            runner_up = max(runner_up, float(others.max()))
        if len(hits):
            at = hits[0]
            needed = min(MIN_OVERLAP, max(2, fragment.units - 1))
            if scan.units[at] < needed:
                return None
            found = (float(scan.score[at]), index, int(offsets[at]))
    if found is None or found[0] < EXTEND_SCORE:
        return None
    return found if found[0] - runner_up >= EXTEND_MARGIN else None


def _anchored_contigs(
    model: barcode.Model, fragments: Sequence[Fragment]
) -> list[Template]:
    """Seeds contigs at the anchored end and extends them round by round."""
    contigs = _seed_contigs(model, fragments)
    free = sorted(
        (f for f in fragments if not f.head and not f.tail),
        key=lambda f: (-f.units, f.read_id),
    )
    stopped: set[int] = set()
    while free and contigs:
        candidates: dict[int, list[tuple[float, int, Fragment]]] = {}
        for fragment in free:
            placement = _unique_placement(model, fragment, contigs)
            if placement is None:
                continue
            score, index, offset = placement
            if (
                index in stopped
                and offset + fragment.units > contigs[index].rows
            ):
                continue
            candidates.setdefault(index, []).append((score, offset, fragment))
        added: set[str] = set()
        for index, found in sorted(candidates.items()):
            found.sort(key=lambda c: (-c[0], c[2].read_id))
            trial = contigs[index].copy()
            accepted = []
            conflicts = []
            for _, offset, fragment in found:
                placement = model.place(
                    trial.matrix, fragment.matrix, [offset]
                )[0]
                if model.consistent(placement):
                    trial.add(fragment, offset)
                    accepted.append(fragment)
                else:
                    conflicts.append((offset, fragment))
            end = contigs[index].rows
            if sum(o + f.units > end for o, f in conflicts) >= BRANCH_READS:
                # The reads disagree on what follows: keep only the inside.
                stopped.add(index)
                for _, offset, fragment in found:
                    if offset + fragment.units <= end:
                        contigs[index].add(fragment, offset)
                        added.add(fragment.read_id)
                continue
            contigs[index] = trial
            added.update(f.read_id for f in accepted)
        if not added:
            break
        free = [f for f in free if f.read_id not in added]
    contigs.sort(key=lambda c: (-len(c.placed), -c.rows))
    return contigs


def build_contigs(
    model: barcode.Model, fragments: Sequence[Fragment]
) -> tuple[list[Template], list[Template]]:
    """Builds proximal and distal contigs.

    Distal contigs are returned in array order: their last row is the last unit.
    """
    proximal = _anchored_contigs(model, fragments)
    flipped = _anchored_contigs(model, [_flip(f) for f in fragments])
    originals = {f.read_id: f for f in fragments}
    distal = []
    for contig in flipped:
        placed = [
            (originals[f.read_id], contig.rows - offset - f.units)
            for f, offset in contig.placed
        ]
        distal.append(Template(placed, contig.rows, contig.kind))
    return proximal, distal


@dataclasses.dataclass
class Join:
    """Best placement of a distal contig on a proximal contig.

    Attributes:
        units: array length implied by the placement.
        overlap: units shared by the two contigs.
        score: score of that length.
        margin: lead over the second best length (or over no overlap).
        threshold: score needed, from decoy joins.
        needed_margin: margin needed, from decoy joins.
        consistent_lengths: lengths that the barcodes do not contradict.
    """

    proximal: Template
    distal: Template
    units: int
    overlap: int
    score: float
    margin: float
    threshold: float
    needed_margin: float
    consistent_lengths: list[int]

    @property
    def accepted(self) -> bool:
        return (
            self.consistent_lengths == [self.units]
            and self.overlap >= MIN_JOIN_OVERLAP
            and self.score >= self.threshold
            and self.margin >= self.needed_margin
        )


def _best_two(scores: np.ndarray) -> tuple[int, float, float]:
    """Index and value of the best score, and the second best (at least 0)."""
    best = int(scores.argmax())
    rest = np.delete(scores, best)
    return best, float(scores[best]), max(float(rest.max(initial=0.0)), 0.0)


def score_join(
    model: barcode.Model,
    proximal: Template,
    distal: Template,
    pool_units: np.ndarray,
) -> Join | None:
    """Scores every overlap and calibrates the threshold with decoys.

    Decoys are the distal contig with its units shuffled, and random units of
    the sample in the same shape. Neither has a true placement.
    """
    lengths = np.arange(
        max(proximal.rows, distal.rows), proximal.rows + distal.rows
    )
    if not len(lengths):
        return None
    offsets = lengths - distal.rows
    scan = model.scan(proximal.matrix, distal.matrix, offsets)
    best, score, second = _best_two(scan.score)
    join = Join(
        proximal=proximal,
        distal=distal,
        units=int(lengths[best]),
        overlap=int(scan.units[best]),
        score=score,
        margin=score - second,
        threshold=JOIN_SCORE,
        needed_margin=JOIN_MARGIN,
        consistent_lengths=[int(n) for n in lengths[scan.consistent]],
    )
    if not join.accepted:
        return join
    generator = np.random.default_rng(_SEED)
    scores = np.zeros(_DECOYS)
    margins = np.zeros(_DECOYS)
    for trial in range(_DECOYS):
        if trial % 2 and len(pool_units):
            index = generator.integers(len(pool_units), size=distal.rows)
            decoy = pool_units[index]
        else:
            decoy = distal.matrix[generator.permutation(distal.rows)]
        values = model.scan(
            proximal.matrix, np.ascontiguousarray(decoy), offsets
        ).score
        _, top, runner = _best_two(values)
        scores[trial] = top
        margins[trial] = top - runner
    join.threshold = max(
        JOIN_SCORE, float(np.quantile(scores, _DECOY_QUANTILE))
    )
    join.needed_margin = max(
        JOIN_MARGIN, float(np.quantile(margins, _DECOY_QUANTILE))
    )
    return join


@dataclasses.dataclass
class Assembly:
    """A fully tiled allele and the checks that make its count exact."""

    template: Template
    weakest_link: int
    margin: float
    problems: list[str]

    @property
    def exact(self) -> bool:
        return not self.problems


def _alike(first: np.ndarray, second: np.ndarray) -> bool:
    shared = (first >= 0) & (second >= 0)
    return int((first[shared] != second[shared]).sum()) <= _SAME_UNIT_SITES


def _uncrossed_run(template: Template) -> bool:
    """True if a run of look-alike units is not crossed by any single read."""
    rows = template.rows
    start = 0
    while start < rows:
        end = start
        while end + 1 < rows and _alike(
            template.matrix[end], template.matrix[end + 1]
        ):
            end += 1
        if end > start:
            low, high = max(0, start - 1), min(rows - 1, end + 1)
            crossed = any(
                offset <= low
                and offset + f.units - 1 >= high
                and (start > 0 or f.head)
                and (end < rows - 1 or f.tail)
                for f, offset in template.placed
            )
            if not crossed:
                return True
        start = end + 1
    return False


def _has_unique_unit(model: barcode.Model, join: Join) -> bool:
    """True if one overlapped unit matches its partner far better than others.

    Neighboring units of one array can be nearly identical. Their overlap then
    scores well at several lengths, so a join is trusted only if some unit pins
    the placement.
    """
    head, tail = join.proximal.matrix, join.distal.matrix
    shift = join.units - tail.shape[0]
    units = model.unit_scores(head, tail)
    for j in range(tail.shape[0]):
        i = j + shift
        if i >= head.shape[0]:
            break
        if not units.consistent[i, j]:
            continue
        rivals = np.concatenate(
            [np.delete(units.score[:, j], i), np.delete(units.score[i], j)]
        )
        if units.score[i, j] - rivals.max(initial=0.0) >= JOIN_MARGIN:
            return True
    return False


def _fits_somewhere(
    model: barcode.Model, fragment: Fragment, template: Template
) -> bool:
    """True if the read can lie inside this array."""
    return compare_with_allele(model, fragment, template) == "fits"


def _shares_part(
    model: barcode.Model, fragment: Fragment, template: Template
) -> bool:
    """True if two or more units of the read match a stretch of the array."""
    units = model.unit_scores(template.matrix, fragment.matrix)
    good = units.consistent & (units.score >= EXTEND_SCORE / 2)
    pairs = good[:-1, :-1] & good[1:, 1:]
    return bool(pairs.any())


def suspect_chimera(
    model: barcode.Model, read: Fragment, others: Sequence[Fragment]
) -> bool:
    """True if a spanning read looks like two molecules joined together.

    Inside a real array every unit before the last one is complete. A chimeric
    read often breaks that rule at its junction. Otherwise the other reads
    show the junction: hardly any crosses it, while reads of the true arrays
    follow the chimera up to it and then differ.
    """
    if read.interrupted:
        return True
    links = read.units - 1
    if links < 1:
        return False
    crossing = np.zeros(links, dtype=int)
    diverging = np.zeros(links, dtype=int)
    for other in others:
        if other.read_id == read.read_id or other.units < 2:
            continue
        good = model.unit_scores(read.matrix, other.matrix).consistent
        crosses = np.zeros(links, dtype=bool)
        leaves = np.zeros(links, dtype=bool)
        for offset in range(2 - other.units, read.units - 1):
            first = max(0, offset)
            last = min(read.units, offset + other.units)
            rows = np.arange(first, last)
            match = good[rows, rows - offset]
            pair = match[:-1] & match[1:]
            crosses[rows[:-1][pair]] = True
            # Two matching units next to one that differs.
            for i in range(len(rows) - 2):
                window = match[i : i + 3].tolist()
                if window == [True, True, False]:
                    leaves[rows[i + 1]] = True
                elif window == [False, True, True]:
                    leaves[rows[i]] = True
        crossing += crosses
        diverging += leaves & ~crosses
    return bool(
        (
            (crossing <= MAX_CROSSING_FRACTION * diverging)
            & (diverging >= MIN_DIVERGING_READS)
        ).any()
    )


def corroborated(
    model: barcode.Model, read: Fragment, others: Sequence[Fragment]
) -> bool:
    """True if enough other reads match two or more units of the read."""
    support = 0
    for other in others:
        if other.read_id == read.read_id or other.units < 2:
            continue
        offsets = np.arange(2 - other.units, read.units - 1)
        scan = model.scan(read.matrix, other.matrix, offsets)
        support += bool((scan.consistent & (scan.units >= 2)).any())
        if support >= MIN_CORROBORATING_READS:
            return True
    return False


def check_join(
    model: barcode.Model,
    join: Join,
    fragments: Sequence[Fragment],
    ambiguous: frozenset[str] = frozenset(),
    others: Sequence[Template] = (),
    chromosome: str = "",
) -> Assembly:
    """Merges an accepted join and tests whether its count is exact.

    Args:
        model: scoring model.
        join: an accepted join.
        fragments: all partial reads that may belong to this allele.
        ambiguous: reads that also fit another allele; they cannot be the only
            support of a link.
        others: other arrays of the chromosome that can explain a read.
        chromosome: chromosome that the array is assigned to.
    """
    shift = join.units - join.distal.rows
    placed = list(join.proximal.placed)
    names = {f.read_id for f, _ in placed}
    placed += [
        (f, o + shift) for f, o in join.distal.placed if f.read_id not in names
    ]
    template = Template(
        placed,
        join.units,
        join.proximal.kind or join.distal.kind,
    )
    links = template.link_reads()
    problems = []
    weakest = min(map(len, links), default=0)
    if weakest < MIN_LINK_READS:
        problems.append("weak_link")
    if any(covering and set(covering) <= ambiguous for covering in links):
        problems.append("link_from_ambiguous_reads")
    if sum(f.units >= 2 for f, _ in placed) < MIN_ALLELE_READS:
        problems.append("few_reads")
    if _uncrossed_run(template):
        problems.append("identical_unit_run")
    if not _has_unique_unit(model, join):
        problems.append("no_unique_unit_in_overlap")
    longer = any(
        f.units > join.units
        and model.scan(
            template.matrix,
            f.matrix,
            np.arange(join.units - f.units, 1),
        ).consistent.any()
        for f in fragments
    )
    if longer:
        problems.append("longer_read_fits")
    placed_ids = set(template.read_ids)
    divergent = sum(
        f.read_id not in placed_ids
        and f.units >= 2
        and kinds_agree(f.kind, template.kind)
        and _shares_part(model, f, template)
        and not _fits_somewhere(model, f, template)
        and not any(_fits_somewhere(model, f, other) for other in others)
        for f in fragments
    )
    if any(
        other.rows == template.rows
        and model.scan(other.matrix, template.matrix, [0]).consistent[0]
        for other in others
    ):
        problems.append("same_as_another_array")
    if divergent > MAX_DIVERGENT_READS:
        problems.append("divergent_reads")
    hints = [f.hint for f, _ in template.placed if f.hint]
    foreign = sum(chromosome != "" and h != chromosome for h in hints)
    if foreign and 2 * foreign >= len(hints):
        # Units exchanged between chr4 and chr10 mislead the unit motifs.
        problems.append("markers_of_other_chromosome")
    return Assembly(template, weakest, join.margin, problems)


def _joins(
    model: barcode.Model, fragments: Sequence[Fragment]
) -> tuple[list[Join], list[Template], list[Template]]:
    """Returns accepted, unambiguous joins and both sets of contigs."""
    proximal, distal = build_contigs(model, fragments)
    pool = np.concatenate([f.matrix for f in fragments])
    joins = []
    for head in proximal[:2]:
        for tail in distal[:2]:
            if kinds_agree(head.kind, tail.kind):
                join = score_join(model, head, tail, pool)
                if join is not None and join.accepted:
                    joins.append(join)
    # A contig that joins two partners equally well is ambiguous.
    unique = [
        j
        for j in joins
        if sum(o.proximal is j.proximal for o in joins) == 1
        and sum(o.distal is j.distal for o in joins) == 1
    ]
    return unique[:2], proximal, distal


def _reproducible(
    model: barcode.Model,
    fragments: Sequence[Fragment],
    result: Assembly,
) -> bool:
    """True if random read subsets support the count and no other.

    Many subsets give no count at all, because they lose a needed read. Enough
    of them must give the same count, and at most `MAX_OTHER_COUNTS` may give
    a different one for an array that shares reads with this one.
    """
    members = set(result.template.read_ids)
    generator = np.random.default_rng(_SEED)
    size = int(RESAMPLE_FRACTION * len(fragments))
    same = 0
    different = 0
    for _ in range(RESAMPLES):
        chosen = np.sort(
            generator.choice(len(fragments), size=size, replace=False)
        )
        joins, _, _ = _joins(model, [fragments[i] for i in chosen])
        counts = {
            j.units
            for j in joins
            if 2 * len(members & set(j.proximal.read_ids + j.distal.read_ids))
            > len(set(j.proximal.read_ids + j.distal.read_ids))
        }
        same += counts == {result.template.rows}
        different += bool(counts - {result.template.rows})
        if different > MAX_OTHER_COUNTS:
            return False
    return same >= MIN_AGREEMENT * RESAMPLES


def reassemble(
    model: barcode.Model,
    fragments: Sequence[Fragment],
    ambiguous: frozenset[str] = frozenset(),
    called: Sequence[Template] = (),
    chromosome: str = "",
) -> tuple[list[Assembly], list[Template]]:
    """Rebuilds up to two alleles of one chromosome from partial reads.

    Args:
        model: scoring model of the sample.
        fragments: partial reads of the chromosome.
        ambiguous: reads that also fit an allele that is already called.
        called: arrays of the alleles that are already called.
        chromosome: chromosome that the reads are assigned to.

    Returns:
        Accepted joins with their checks, and all contigs (for lower bounds).
    """
    if not fragments:
        return [], []
    joins, proximal, distal = _joins(model, fragments)
    arrays = [
        Template(
            list(j.proximal.placed)
            + [
                (f, o + j.units - j.distal.rows)
                for f, o in j.distal.placed
                if f.read_id not in j.proximal.read_ids
            ],
            j.units,
        )
        for j in joins
    ]
    assemblies = []
    for index, join in enumerate(joins):
        others = list(called) + [a for i, a in enumerate(arrays) if i != index]
        others += [
            c
            for c in proximal[:2] + distal[:2]
            if c is not join.proximal and c is not join.distal
        ]
        result = check_join(
            model, join, fragments, ambiguous, others, chromosome
        )
        # A rebuilt count equal to a called one is more often the same allele
        # tiled again than a second allele of equal length.
        if any(array.rows == join.units for array in called):
            result.problems.append("same_count_as_called_allele")
        if result.exact and not _reproducible(model, fragments, result):
            result.problems.append("unstable_in_subsets")
        assemblies.append(result)
    return assemblies, proximal[:2] + distal[:2]
