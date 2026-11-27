"""Checks reassembly against counts known from spanning reads.

* V1 removes every spanning read and rebuilds from the partial reads.
* V2 cuts spanning reads into pieces and rebuilds from the pieces.

A wrong exact count, or a lower bound above the known count, is a failure. In
V1 a chromosome with only one known allele cannot judge any other count, which
may belong to the allele that no read spans.
"""

from collections.abc import Mapping, Sequence

import numpy as np

from fshdetector import alleles as alleles_module
from fshdetector import assemble
from fshdetector import barcode

TRIALS = 20
_SEED = 20261126
COLUMNS = (
    "test", "chromosome", "known_ru", "trial", "exact_ru", "overlap_units",
    "contig_lower_bound", "outcome",
)  # fmt: skip


def _outcome(
    exact: Sequence[int],
    bound: int,
    known: Sequence[int],
    complete: bool,
) -> str:
    """Judges a rebuild.

    Args:
        complete: every allele of the data is known, so any other result is
            wrong. If False, another count may belong to the unseen allele.
    """
    if any(ru not in known for ru in exact):
        return "wrong_exact" if complete else "unverifiable"
    if exact:
        return "correct"
    if bound > max(known):
        return "bound_exceeds_known" if complete else "unverifiable"
    return "unresolved"


def _rebuild(
    model: barcode.Model, fragments: Sequence[assemble.Fragment]
) -> tuple[list[int], list[int], int]:
    """Returns exact counts, their overlaps, and the longest contig."""
    results, contigs = assemble.reassemble(model, fragments)
    exact = [r for r in results if r.exact]
    overlaps = []
    for result in exact:
        covered = [
            (offset, offset + f.units, f.head)
            for f, offset in result.template.placed
        ]
        head_end = max((e for _, e, head in covered if head), default=0)
        tail_start = min(
            (s for s, _, head in covered if not head), default=head_end
        )
        overlaps.append(head_end - tail_start)
    return (
        [r.template.rows for r in exact],
        overlaps,
        max((c.supported_rows() for c in contigs), default=0),
    )


def cut(
    read: assemble.Fragment, generator: np.random.Generator
) -> list[assemble.Fragment]:
    """Cuts one spanning read into a head, a tail, and a middle piece."""
    units = read.units
    pieces = []
    for label, head, tail in (
        ("head", True, False),
        ("tail", False, True),
        ("mid", False, False),
    ):
        if units < 3:
            size = units - 1
        else:
            size = int(generator.integers(2, units))
        if size < 1:
            continue
        if head:
            start = 0
        elif tail:
            start = units - size
        else:
            start = int(generator.integers(0, units - size + 1))
        pieces.append(
            assemble.Fragment(
                f"{read.read_id}:{label}",
                np.ascontiguousarray(read.matrix[start : start + size]),
                head,
                tail,
                read.kind if tail else "",
            )
        )
    return pieces


def run(
    model: barcode.Model,
    alleles: Sequence[alleles_module.Allele],
    fragments: Mapping[str, assemble.Fragment],
    partial: Mapping[str, Sequence[assemble.Fragment]],
) -> list[tuple]:
    """Runs V1 and V2 on every chromosome that has a spanning allele.

    Args:
        model: scoring model of the sample.
        alleles: reported alleles.
        fragments: all reads by name.
        partial: chromosome to its reads without both ends.

    Returns:
        Rows in the order of `COLUMNS`.
    """
    rows = []
    generator = np.random.default_rng(_SEED)
    spanning = [a for a in alleles if a.source == alleles_module.SPANNING]
    for chromosome, reads in partial.items():
        counts = [a.ru for a in spanning if a.chromosome == chromosome]
        known = sorted(set(counts))
        if not known:
            continue
        exact, overlaps, bound = _rebuild(model, reads)
        rows.append(
            (
                "V1_hide_spanning",
                chromosome,
                ",".join(map(str, known)),
                1,
                ",".join(map(str, exact)) or None,
                ",".join(map(str, overlaps)) or None,
                bound,
                _outcome(exact, bound, known, len(counts) == 2),
            )
        )
    for allele in spanning:
        own = [
            fragments[r]
            for r in allele.read_ids
            if fragments[r].head and fragments[r].tail
        ]
        if len(own) < 2:
            continue
        for trial in range(1, TRIALS + 1):
            pieces = [p for read in own for p in cut(read, generator)]
            exact, overlaps, bound = _rebuild(model, pieces)
            rows.append(
                (
                    "V2_cut_spanning",
                    allele.chromosome,
                    allele.ru,
                    trial,
                    ",".join(map(str, exact)) or None,
                    ",".join(map(str, overlaps)) or None,
                    bound,
                    _outcome(exact, bound, [allele.ru], True),
                )
            )
    return rows
