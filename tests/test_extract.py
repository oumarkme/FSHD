"""Read collection through the alignment index."""

import conftest
import pysam
import pytest

from fshdetector import extract
from fshdetector import layout

WINDOW = layout.Region("chr4", 10_000, 20_000)
HEADER = {
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chr4", "LN": 500_000}, {"SN": "chr10", "LN": 500_000}],
}


class _Everything:

    def matches(self, sequence: str) -> bool:
        return True


def _record(name, start, sequence, flag=0, sa=None) -> pysam.AlignedSegment:
    record = pysam.AlignedSegment()
    record.query_name = name
    record.flag = flag
    record.query_sequence = sequence
    if flag & 4:
        record.reference_id = -1
        record.reference_start = -1
    else:
        record.reference_id = 0
        record.reference_start = start
        record.cigartuples = [(0, len(sequence))]
    if sa:
        record.set_tag("SA", sa)
    return record


def _write(path, records):
    header = pysam.AlignmentHeader.from_dict(HEADER)
    extract.write_bam(path, header, records)


def _collect(path):
    with extract.alignment(path) as source:
        return extract.collect(
            source,
            [WINDOW],
            _Everything(),
            lambda: extract.open_alignment(path),
        )


def test_collect_completes_reads_from_their_primary_record(tmp_path):
    full = conftest.random_dna(400, 1)
    path = tmp_path / "in.bam"
    _write(
        path,
        [
            _record("inside", 12_000, conftest.random_dna(300, 2)),
            # Primary far outside the window, reached through the SA tag.
            _record("split", 300_000, full),
            _record(
                "split", 15_000, full[:100], 0x800, "chr4,300001,+,400M,60,0;"
            ),
            _record("secondary", 16_000, full[:100], 0x100),
            _record("outside", 100_000, conftest.random_dna(300, 3)),
            _record("unmapped", -1, conftest.random_dna(300, 4), 4),
        ],
    )
    reads = _collect(path)
    assert sorted(reads.candidates) == ["inside", "split", "unmapped"]
    assert reads.candidates["split"] == full
    assert reads.counts["primaries_recovered_outside_windows"] == 1
    assert reads.counts["secondary_only_reads"] == 1


def test_damaged_input_fails_and_names_the_region(tmp_path):
    path = tmp_path / "in.bam"
    _write(
        path,
        [
            _record(f"r{i}", 10_000 + i, conftest.random_dna(2000, i))
            for i in range(2000)
        ],
    )
    data = bytearray(path.read_bytes())
    middle = len(data) // 2
    data[middle : middle + 2000] = bytes(2000)
    path.write_bytes(bytes(data))
    with pytest.raises(OSError, match="cannot decode chr4:10001-20000"):
        _collect(path)


def test_cram_requires_a_reference(tmp_path):
    with pytest.raises(ValueError, match="reference"):
        extract.open_alignment(tmp_path / "in.cram")
