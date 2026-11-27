"""Collects candidate reads through the alignment index.

Only the D4Z4 windows and the unmapped block are read, never the whole file.
"""

from collections.abc import Callable, Iterable, Iterator, Sequence
import concurrent.futures
import contextlib
import dataclasses
import pathlib
import threading

import pysam

from fshdetector import layout as layout_module

_NOT_PRIMARY = 0x100 | 0x800
_SAME_REQUEST = 100_000  # SA positions this close are fetched together.


@dataclasses.dataclass
class Reads:
    """Reads collected from one alignment file.

    Attributes:
        header: header of the source file.
        window_records: every alignment record that overlaps a window.
        primaries: read name to the record that holds the full sequence.
        candidates: read name to the full sequence, in sequencing orientation,
            for reads that passed the k-mer screen.
        searched_bases: bases of all screened reads, used as the BLAST
            database size.
        counts: counters for the run report.
    """

    header: pysam.AlignmentHeader
    window_records: list[pysam.AlignedSegment]
    primaries: dict[str, pysam.AlignedSegment]
    candidates: dict[str, str]
    searched_bases: int
    counts: dict[str, int]


def open_alignment(
    path: str | pathlib.Path,
    reference: str | pathlib.Path | None = None,
    threads: int = 1,
) -> pysam.AlignmentFile:
    """Opens an indexed BAM or CRAM.

    Raises:
        ValueError: a CRAM has no reference, or the index is missing.
    """
    path = str(path)
    cram = path.endswith(".cram")
    if cram and reference is None:
        raise ValueError(f"{path} is a CRAM; its reference FASTA is required")
    handle = pysam.AlignmentFile(
        path,
        "rc" if cram else "rb",
        reference_filename=str(reference) if cram else None,
        threads=max(1, threads),
    )
    if not handle.has_index():
        handle.close()
        raise ValueError(f"{path} has no index (.bai or .crai)")
    return handle


@contextlib.contextmanager
def alignment(
    path: str | pathlib.Path,
    reference: str | pathlib.Path | None = None,
    threads: int = 1,
) -> Iterator[pysam.AlignmentFile]:
    """`open_alignment` as a context that keeps the first error.

    Closing a damaged file fails too, and that message would hide the cause.
    """
    handle = open_alignment(path, reference, threads)
    try:
        yield handle
    except BaseException:
        with contextlib.suppress(OSError):
            handle.close()
        raise
    handle.close()


def _fetch(
    source: pysam.AlignmentFile, *region
) -> Iterator[pysam.AlignedSegment]:
    """Fetches a region (contig, 0-based start, end) or the unmapped reads.

    Raises:
        OSError: the records cannot be decoded. The message names the region;
            a partial read set would give wrong counts without notice.
    """
    try:
        yield from source.fetch(*region)
    except OSError as error:
        where = (
            "the unmapped reads"
            if region == ("*",)
            else f"{region[0]}:{region[1] + 1}-{region[2]}"
        )
        raise OSError(
            f"cannot decode {where} of {source.filename.decode()} ({error}); "
            "the file is damaged. Check it with `samtools view -c`."
        ) from error


def _is_primary(record: pysam.AlignedSegment) -> bool:
    return not record.flag & _NOT_PRIMARY and record.query_sequence is not None


def _recover_primaries(
    reopen: Callable[[], pysam.AlignmentFile],
    records: Iterable[pysam.AlignedSegment],
    workers: int,
) -> dict[str, pysam.AlignedSegment]:
    """Finds primary records of other loci through the SA tags.

    Each random access into a CRAM decodes a whole block and takes a noticeable
    fraction of a second, so neighboring positions share one request and the
    requests run in parallel, each thread with its own file handle.

    Args:
        reopen: opens another handle of the same alignment file.
        records: one supplementary record per read that lacks its primary.
        workers: number of parallel requests.
    """
    wanted: dict[str, set[tuple[int, str]]] = {}
    for record in records:
        if record.has_tag("SA"):
            for entry in record.get_tag("SA").rstrip(";").split(";"):
                chrom, position = entry.split(",")[:2]
                wanted.setdefault(chrom, set()).add(
                    (int(position) - 1, record.query_name)
                )
    requests: list[tuple[str, list[tuple[int, str]]]] = []
    for chrom, positions in sorted(wanted.items()):
        for item in sorted(positions):
            if (
                requests
                and requests[-1][0] == chrom
                and item[0] - requests[-1][1][-1][0] <= _SAME_REQUEST
            ):
                requests[-1][1].append(item)
            else:
                requests.append((chrom, [item]))
    if not requests:
        return {}

    local = threading.local()
    handles = []

    def fetch(request) -> list[pysam.AlignedSegment]:
        chrom, items = request
        if not hasattr(local, "handle"):
            local.handle = reopen()
            handles.append(local.handle)
        keys = set(items)
        return [
            other
            for other in _fetch(
                local.handle, chrom, items[0][0], items[-1][0] + 1
            )
            if (other.reference_start, other.query_name) in keys
            and _is_primary(other)
        ]

    try:
        with concurrent.futures.ThreadPoolExecutor(max(1, workers)) as pool:
            results = list(pool.map(fetch, requests))
    finally:
        for handle in handles:
            # Nothing was written; a close error would only hide a read error.
            with contextlib.suppress(OSError):
                handle.close()
    return {r.query_name: r for found in results for r in found}


def collect(
    source: pysam.AlignmentFile,
    windows: Sequence[layout_module.Region],
    screen,
    reopen: Callable[[], pysam.AlignmentFile],
    workers: int = 1,
) -> Reads:
    """Reads the windows and the unmapped block, then screens each read.

    Args:
        source: an open, indexed alignment file.
        windows: regions whose reads are collected.
        screen: object with ``matches(sequence) -> bool``.
        reopen: opens another handle of the same alignment file.
        workers: parallel requests for primaries outside the windows.

    Returns:
        The collected reads. Reads whose only records in the windows are
        supplementary are completed with their primary record from elsewhere.
        Reads with only secondary records belong to another locus and are
        counted but not used.
    """
    window_records = []
    for region in sorted(windows, key=lambda r: source.get_tid(r.chrom)):
        window_records.extend(
            _fetch(source, region.chrom, region.start, region.end)
        )
    primaries = {r.query_name: r for r in window_records if _is_primary(r)}
    lacking = {
        r.query_name: r
        for r in window_records
        if r.is_supplementary and r.query_name not in primaries
    }
    recovered = _recover_primaries(reopen, lacking.values(), workers)
    primaries.update(recovered)
    without_primary = {r.query_name for r in window_records} - set(primaries)

    candidates = {}
    searched_bases = 0
    for name, record in primaries.items():
        sequence = record.get_forward_sequence()
        searched_bases += len(sequence)
        if screen.matches(sequence):
            candidates[name] = sequence.upper()
    window_reads = len(primaries)

    unmapped = 0
    for record in _fetch(source, "*"):
        sequence = record.query_sequence
        if sequence is None or record.query_name in primaries:
            continue
        unmapped += 1
        searched_bases += len(sequence)
        if screen.matches(sequence):
            primaries[record.query_name] = record
            candidates[record.query_name] = (
                record.get_forward_sequence().upper()
            )

    return Reads(
        header=source.header,
        window_records=window_records,
        primaries=primaries,
        candidates=candidates,
        searched_bases=searched_bases,
        counts={
            "window_records": len(window_records),
            "window_reads": window_reads,
            "primaries_recovered_outside_windows": len(recovered),
            "secondary_only_reads": len(without_primary),
            "unmapped_reads": unmapped,
            "candidate_reads": len(candidates),
        },
    )


def write_bam(
    path: str | pathlib.Path,
    header: pysam.AlignmentHeader,
    records: Iterable[pysam.AlignedSegment],
) -> None:
    """Writes records as a coordinate-sorted, indexed BAM."""

    def order(record):
        tid = record.reference_id
        return (tid if tid >= 0 else 1 << 30, record.reference_start)

    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for record in sorted(records, key=order):
            out.write(record)
    pysam.index(str(path))
