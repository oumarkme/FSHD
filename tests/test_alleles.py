"""Reassembly and the evidence tiers of the allele calls."""

import conftest
import numpy as np

from fshdetector import alleles
from fshdetector import annotate
from fshdetector import assemble

UNITS = 12


def _read(
    array: np.ndarray, start: int, end: int, name: str
) -> assemble.Fragment:
    return assemble.Fragment(
        name,
        np.ascontiguousarray(array[start:end]),
        head=start == 0,
        tail=end == len(array),
    )


def _partial_reads(array: np.ndarray) -> list[assemble.Fragment]:
    """Overlapping reads that tile the array but never span it."""
    n = len(array)
    intervals = [(0, k) for k in (5, 6, 7, 8)]
    intervals += [(n - k, n) for k in (5, 6, 7, 8)]
    intervals += [(s, s + 6) for s in range(1, n - 6)] * 2
    return [
        _read(array, start, end, f"p{i}")
        for i, (start, end) in enumerate(intervals)
    ]


def _call(fragment: assemble.Fragment) -> annotate.ReadCall:
    return annotate.ReadCall(
        read_id=fragment.read_id,
        chromosome="chr4",
        chr_confidence="HIGH",
        ru_count=fragment.units,
        haplotype="4q",
        length_subtype="",
        permissive=False,
        read_length=1,
        has_head=fragment.head,
        has_tail=fragment.tail,
        proximal_marker="",
        distal_marker="",
        chimeric=False,
        blni_chr4=0,
        blni_chr10=0,
    )


def _chr4_alleles(reads: list[assemble.Fragment], model) -> list:
    called = alleles.call_alleles(
        model,
        [_call(f) for f in reads],
        {f.read_id: f for f in reads},
    )
    return [a for a in called if a.chromosome == "chr4"]


def test_reassembly_rebuilds_an_array_from_partial_reads():
    array = conftest.array_matrix(UNITS)
    results, _ = assemble.reassemble(
        conftest.sample_model(array), _partial_reads(array)
    )
    assert [r.template.rows for r in results if r.exact] == [UNITS]


def test_reassembly_gives_no_count_for_identical_units():
    unit = conftest.array_matrix(1)
    array = np.repeat(unit, UNITS, axis=0)
    other = conftest.array_matrix(UNITS, seed=2)
    results, _ = assemble.reassemble(
        conftest.sample_model(array, other), _partial_reads(array)
    )
    assert not [r for r in results if r.exact]


def test_reassembly_rejects_the_count_of_a_called_allele():
    array = conftest.array_matrix(UNITS)
    other = conftest.array_matrix(UNITS, seed=2)
    called = assemble.Template([(_read(other, 0, UNITS, "s"), 0)], UNITS)
    results, _ = assemble.reassemble(
        conftest.sample_model(array, other),
        _partial_reads(array),
        called=[called],
    )
    assert results
    assert all("same_count_as_called_allele" in r.problems for r in results)


def test_spanning_reads_come_before_reassembly():
    first = conftest.array_matrix(5, seed=3)
    second = conftest.array_matrix(UNITS)
    reads = [_read(first, 0, 5, f"s{i}") for i in range(3)]
    reads += _partial_reads(second)
    found = _chr4_alleles(reads, conftest.sample_model(first, second))
    assert [(a.source, a.reported_ru) for a in found] == [
        ("spanning", "5"),
        ("reassembly", str(UNITS)),
    ]
    assert found[0].spanning_reads == 3


def test_longest_read_gives_a_lower_bound():
    first = conftest.array_matrix(5, seed=3)
    second = conftest.array_matrix(UNITS)
    reads = [_read(first, 0, 5, f"s{i}") for i in range(2)]
    # Proximal reads only: no count is possible.
    reads += [_read(second, 0, k, f"p{k}") for k in (7, 8, 9)]
    found = _chr4_alleles(reads, conftest.sample_model(first, second))
    assert [(a.source, a.reported_ru) for a in found] == [
        ("spanning", "5"),
        ("longest_read", ">=9"),
    ]


def test_missing_second_allele_is_reported_as_such():
    first = conftest.array_matrix(5, seed=3)
    reads = [_read(first, 0, 5, f"s{i}") for i in range(2)]
    found = _chr4_alleles(reads, conftest.sample_model(first))
    assert [a.reported_ru for a in found] == ["5", "NA"]
    assert found[1].flags == ["second_allele_not_seen"]
    other = [
        a
        for a in alleles.call_alleles(
            conftest.sample_model(first),
            [_call(f) for f in reads],
            {f.read_id: f for f in reads},
        )
        if a.chromosome == "chr10"
    ]
    assert [a.flags for a in other] == [["no_reads"], ["no_reads"]]
