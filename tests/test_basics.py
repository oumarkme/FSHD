"""Screen, genome build detection, and methylation helpers."""

import array

import conftest
import pysam
import pytest

from fshdetector import _screen
from fshdetector import layout
from fshdetector import methylation
from fshdetector import reference
from fshdetector import sequences


def test_screen_finds_a_query_on_both_strands():
    query = conftest.random_dna(500, 1)
    screen = _screen.KmerScreen([query])
    read = conftest.random_dna(300, 2) + query[100:160]
    assert screen.matches(read)
    assert screen.matches(sequences.reverse_complement(read))
    assert not screen.matches(conftest.random_dna(5000, 3))
    assert not screen.matches("ACGT")


def test_screen_ignores_simple_repeats_of_extra_sequences():
    extra = "A" * 60 + conftest.random_dna(200, 4)
    screen = _screen.KmerScreen([], extra=[extra])
    assert not screen.matches("A" * 200)
    assert screen.matches(extra[60:])


@pytest.mark.parametrize(
    "length, build", [(193574945, "t2t"), (190214555, "hg38")]
)
def test_detect_build(length, build):
    lengths = {"chr4": length, "chr10": 1}
    assert reference.detect_build(lengths, layout.HUMAN) == build


@pytest.mark.parametrize(
    "lengths",
    [
        {"chr4": 191154276, "chr10": 1},  # hg19
        {"chr4": 5, "chr10": 1},
        {"4": 190214555, "10": 1},
    ],
)
def test_detect_build_rejects_unsupported_input(lengths):
    with pytest.raises(reference.ReferenceError_):
        reference.detect_build(lengths, layout.HUMAN)


def test_summarize_averages_covered_sites_inside_the_region():
    region = layout.Region("chr4", 10, 20, "DR1")
    sites = [
        methylation.Site("chr4", 10, 5, 20.0),
        methylation.Site("chr4", 19, 1, 60.0),
        methylation.Site("chr4", 20, 9, 100.0),  # Half-open end.
        methylation.Site("chr4", 12, 0, 100.0),  # Not covered.
        methylation.Site("chr10", 12, 3, 100.0),
    ]
    assert methylation.summarize(sites, region) == methylation.Summary(40.0, 2)
    assert methylation.summarize([], region) == methylation.Summary(None, 0)


def test_forward_interval_flips_reverse_reads():
    assert methylation.forward_interval(100, 10, 30, False) == (10, 30)
    assert methylation.forward_interval(100, 10, 30, True) == (70, 90)


def _record(sequence: str, mm: str, ml: list[int]) -> pysam.AlignedSegment:
    record = pysam.AlignedSegment()
    record.query_name = "r"
    record.query_sequence = sequence
    record.flag = 4
    record.set_tag("MM", mm)
    record.set_tag("ML", array.array("B", ml))
    return record


def test_read_cpg_calls_skips_uncalled_and_uncertain_sites():
    # CpG cytosines at 1, 5, and 8; the tag lists the first and third.
    record = _record("ACGTACGTCG", "C+m?,0,1;", [250, 5])
    assert methylation.read_cpg_calls(record, 0.8) == {1: True, 8: False}
    middle = _record("ACGTACGTCG", "C+m?,0,1;", [128, 5])
    assert methylation.read_cpg_calls(middle, 0.8) == {8: False}
    implicit = _record("ACGTACGTCG", "C+m.,0,1;", [250, 5])
    assert methylation.read_cpg_calls(implicit, 0.8)[5] is False
