"""Unit positions on reads, chromosome votes, and barcode scoring."""

import conftest
import numpy as np

from fshdetector import annotate
from fshdetector import barcode
from fshdetector import blast

PERIOD = 3300
CONSENSUS = conftest.random_dna(PERIOD, 11)
LENGTH = 30_000


def _hit(start: int, reverse: bool, length: int = PERIOD) -> blast.Hit:
    """A full-unit hit that starts at array coordinate `start` (0-based)."""
    if reverse:
        sstart = LENGTH - start
        send = sstart - length + 1
    else:
        sstart = start + 1
        send = sstart + length - 1
    return blast.Hit(
        "d4z4", "read", 99.0, length, 1, length, sstart, send, 0.0,
        float(length), str(length),
    )  # fmt: skip


def _units(hits: list[blast.Hit]) -> barcode.ReadUnits:
    return barcode.build_units(hits, {"read": LENGTH}, CONSENSUS)["read"]


def test_units_keep_their_gaps_on_both_strands():
    starts = [1000 + i * PERIOD for i in (0, 1, 3)]
    for reverse in (False, True):
        read = _units([_hit(s, reverse) for s in starts])
        assert read.reverse == reverse
        assert [u is not None for u in read.units] == [True, True, False, True]
        assert read.counted == 4
        array = (1000, 1000 + 4 * PERIOD)
        expected = (LENGTH - array[1], LENGTH - array[0]) if reverse else array
        assert read.span() == expected
        assert not read.mixed_strand


def test_short_terminal_pieces_are_not_counted():
    read = _units([_hit(1000, False), _hit(1000 + PERIOD, False, length=900)])
    assert len(read.units) == 2
    assert read.counted == 1


def test_fold_back_read_is_flagged():
    read = _units(
        [_hit(1000, False), _hit(1000 + PERIOD, False), _hit(20_000, True)]
    )
    assert read.mixed_strand
    assert read.counted == 2


def test_restriction_motifs_count_only_inside_the_array():
    # XapI sites are frequent in flanking sequence and must not vote.
    flank = "GAATTC" * 200
    array = "ACAC" * 100 + "CCTAGGCCTAGG" + "ACAC" * 724
    sequence = flank + array + flank
    unit = barcode.Unit(
        bases=np.zeros(1, np.int8),
        inserted=np.zeros(1, np.int8),
        length=len(array),
        start=len(flank),
        end=len(flank) + len(array),
    )
    read = barcode.ReadUnits("read", False, len(sequence), [unit])
    calls = annotate.annotate(
        {"read": read}, [], [], [], {"read": sequence}, {"read": "chr4"}
    )
    assert calls[0].chromosome == "chr10"


def test_allele_counts_leave_ties_unknown():
    counts = barcode.AlleleCounts(2)
    counts.add(np.array([[0, 1]], np.int8), 0)
    counts.add(np.array([[1, 1], [0, -1]], np.int8), 0)
    assert counts.consensus(2).tolist() == [[-1, 1], [0, -1]]


def test_model_separates_same_and_different_units(generator):
    array = conftest.array_matrix(6)
    model = conftest.sample_model(array)
    noisy = array[2:5].copy()
    noisy[generator.random(noisy.shape) < 0.01] ^= 1
    scan = model.scan(array, noisy, [0, 1, 2, 3])
    assert scan.consistent.tolist() == [False, False, True, False]
    assert scan.score.argmax() == 2


def test_error_estimate_follows_the_disagreement(generator):
    truth = conftest.array_matrix(20, sites=200)
    reads = []
    for _ in range(10):
        read = truth.copy()
        read[generator.random(read.shape) < 0.05] ^= 1
        reads.append(read)
    error = barcode.estimate_error([reads])
    assert 0.035 < error.mean() < 0.065
    assert barcode.estimate_error([reads[:1]]) is None
