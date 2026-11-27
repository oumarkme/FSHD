"""CpG methylation from modkit pileups and from single reads.

modkit reads the whole reference it is given. The pileup therefore runs on a
reference cut to the target regions, with record positions shifted to match;
positions are shifted back when the results are read.

The pass threshold is always given to modkit explicitly. With modkit 0.6.4 a
threshold that was only estimated was not applied to a BAM that listed every
chromosome in its header.
"""

from collections.abc import Iterable, Sequence
import dataclasses
import pathlib
import re

import pysam

from fshdetector import extract
from fshdetector import layout as layout_module
from fshdetector import tools

_THRESHOLD = re.compile(r"(?:setting threshold|pass threshold) ([0-9.]+)")
_PADDING = 100
PARTITIONS = ("hp1", "hp2", "combined")


@dataclasses.dataclass(frozen=True)
class Site:
    """Methylation of one CpG (0-based T2T start of the cytosine)."""

    chrom: str
    start: int
    n_valid: int
    percent: float


@dataclasses.dataclass(frozen=True)
class Summary:
    """Methylation of one region: unweighted mean over its CpGs.

    ``percent`` is None when the region has no covered CpG.
    """

    percent: float | None
    n_cpg: int


def summarize(sites: Iterable[Site], region: layout_module.Region) -> Summary:
    """Averages the per-CpG percentages inside a region."""
    values = [
        s.percent
        for s in sites
        if s.chrom == region.chrom
        and region.start <= s.start < region.end
        and s.n_valid > 0
    ]
    if not values:
        return Summary(None, 0)
    return Summary(round(sum(values) / len(values), 2), len(values))


def read_bedmethyl(path: pathlib.Path, offsets: dict[str, int]) -> list[Site]:
    """Reads a bedMethyl table and shifts positions back to T2T."""
    if not path.is_file():
        return []
    sites = []
    with open(path, encoding="ascii") as handle:
        for line in handle:
            f = line.split()
            sites.append(
                Site(
                    chrom=f[0],
                    start=int(f[1]) + offsets[f[0]],
                    n_valid=int(f[9]),
                    percent=float(f[10]),
                )
            )
    return sites


class Pileup:
    """Cut-down reference and records for repeated modkit runs.

    Args:
        records: alignment records on T2T; only primary records that overlap a
            target and carry modified base tags are used, as modkit does.
        reference: FASTA that holds chr4 and chr10 of T2T.
        targets: regions to pile up.
        workdir: scratch directory.
    """

    def __init__(
        self,
        records: Iterable[pysam.AlignedSegment],
        reference: str | pathlib.Path,
        targets: Sequence[layout_module.Region],
        workdir: pathlib.Path,
    ):
        self.workdir = workdir
        self.targets = tuple(targets)
        by_chrom = {t.chrom: t for t in targets}
        self.records = []
        for record in records:
            if record.is_unmapped or record.is_secondary:
                continue
            if record.is_supplementary or not record.has_tag("MM"):
                continue
            target = by_chrom.get(record.reference_name)
            if (
                target is not None
                and record.reference_start < target.end
                and record.reference_end > target.start
            ):
                self.records.append(record)
        self.offsets: dict[str, int] = {}
        self._lengths: dict[str, int] = {}
        self.reference = workdir / "targets.fa"
        self.bed = workdir / "targets.bed"
        with (
            pysam.FastaFile(str(reference)) as fasta,
            open(self.reference, "w", encoding="ascii") as out,
        ):
            for target in targets:
                chrom = target.chrom
                own = [r for r in self.records if r.reference_name == chrom]
                low = min([target.start] + [r.reference_start for r in own])
                high = max([target.end] + [r.reference_end for r in own])
                low = max(0, low - _PADDING)
                high = min(fasta.get_reference_length(chrom), high + _PADDING)
                self.offsets[chrom] = low
                self._lengths[chrom] = high - low
                out.write(f">{chrom}\n{fasta.fetch(chrom, low, high)}\n")
        pysam.faidx(str(self.reference))
        with open(self.bed, "w", encoding="ascii") as out:
            for t in targets:
                low = self.offsets[t.chrom]
                out.write(f"{t.chrom}\t{t.start - low}\t{t.end - low}\n")

    def _write_bam(
        self, path: pathlib.Path, haplotypes: dict[str, int] | None
    ) -> int:
        header = pysam.AlignmentHeader.from_dict(
            {
                "HD": {"VN": "1.6", "SO": "coordinate"},
                "SQ": [
                    {"SN": t.chrom, "LN": self._lengths[t.chrom]}
                    for t in self.targets
                ],
            }
        )
        shifted = []
        for record in self.records:
            chrom = record.reference_name
            copy = record.to_dict()
            copy["ref_pos"] = str(
                record.reference_start - self.offsets[chrom] + 1
            )
            copy["next_ref_name"] = "*"
            copy["next_ref_pos"] = "0"
            if haplotypes is not None:
                tags = [t for t in copy["tags"] if not t.startswith("HP:")]
                group = haplotypes.get(record.query_name)
                if group is not None:
                    tags.append(f"HP:i:{group}")
                copy["tags"] = tags
            shifted.append(pysam.AlignedSegment.from_dict(copy, header))
        extract.write_bam(path, header, shifted)
        return len(shifted)

    def _modkit(
        self,
        bam: pathlib.Path,
        output: pathlib.Path,
        threads: int,
        threshold: float | None,
    ) -> pathlib.Path:
        log = output.with_suffix(".log")
        command = [
            "modkit", "pileup", bam, output,
            "--ref", self.reference,
            "--cpg",
            "--combine-strands",
            "--modified-bases", "5mC",
            "--phased",
            "--include-bed", self.bed,
            "--prefix", "p",
            "--threads", threads,
            "--log-filepath", log,
        ]  # fmt: skip
        if threshold is not None:
            command += ["--filter-threshold", f"{threshold:.7f}"]
        tools.run(command)
        return log

    def estimate_threshold(self, threads: int) -> float | None:
        """Returns modkit's pass threshold for these records, or None."""
        if not self.records:
            return None
        bam = self.workdir / "estimate.bam"
        self._write_bam(bam, None)
        log = self._modkit(bam, self.workdir / "estimate", threads, None)
        found = _THRESHOLD.search(log.read_text(encoding="utf-8"))
        if found is None:
            raise tools.ToolError(f"no pass threshold found in {log}")
        return float(found.group(1))

    def run(
        self,
        name: str,
        threshold: float | None,
        threads: int,
        haplotypes: dict[str, int] | None = None,
    ) -> dict[str, list[Site]]:
        """Piles up the records, split by HP tag.

        Args:
            name: label of this run, used for file names.
            threshold: pass threshold from `estimate_threshold`.
            threads: modkit threads.
            haplotypes: read name to HP value (1 or 2) that replaces the HP
                tags of the input; None keeps the input tags.

        Returns:
            Sites of the "hp1", "hp2", and "combined" (all reads) partitions.
        """
        if not self.records or threshold is None:
            return {part: [] for part in PARTITIONS}
        bam = self.workdir / f"{name}.bam"
        self._write_bam(bam, haplotypes)
        output = self.workdir / name
        self._modkit(bam, output, threads, threshold)
        return {
            part: read_bedmethyl(output / f"p_{part}.bedmethyl", self.offsets)
            for part in PARTITIONS
        }


def has_implicit_canonical(record: pysam.AlignedSegment) -> bool:
    """True if CpGs absent from the MM tag mean "not modified"."""
    for part in str(record.get_tag("MM")).split(";"):
        code = part.split(",", 1)[0]
        if code.startswith("C+"):
            return not code.endswith("?")
    return False


def read_cpg_calls(
    record: pysam.AlignedSegment, threshold: float
) -> dict[int, bool]:
    """Calls the CpGs of one read.

    Args:
        record: primary record, which holds the full sequence.
        threshold: minimum probability of the most likely state.

    Returns:
        Position of the C in sequencing orientation to True (5mC) or False
        (unmodified or 5hmC). Positions below the threshold are omitted.
    """
    if not record.has_tag("MM"):
        return {}
    modified = record.modified_bases_forward or {}
    probabilities: dict[int, dict[str, float]] = {}
    for (base, _, code), calls in modified.items():
        if base != "C" or code not in ("m", "h"):
            continue
        for position, quality in calls:
            if quality >= 0:
                probabilities.setdefault(position, {})[code] = (
                    quality + 0.5
                ) / 256
    sequence = record.get_forward_sequence().upper()
    implicit = has_implicit_canonical(record)
    result = {}
    position = sequence.find("CG")
    while position >= 0:
        states = probabilities.get(position)
        if states is None:
            if implicit:
                result[position] = False
        else:
            m = states.get("m", 0.0)
            h = states.get("h", 0.0)
            c = max(0.0, 1.0 - m - h)
            if max(m, h, c) >= threshold:
                result[position] = m > max(h, c)
        position = sequence.find("CG", position + 1)
    return result


def interval_calls(
    calls: dict[int, bool], start: int, end: int
) -> tuple[int, int]:
    """Counts (5mC calls, all calls) inside a 0-based half-open interval."""
    inside = [value for pos, value in calls.items() if start <= pos < end]
    return sum(inside), len(inside)


def forward_interval(
    length: int, start: int, end: int, reverse: bool
) -> tuple[int, int]:
    """Converts an array-oriented interval back to read coordinates."""
    if reverse:
        return length - end, length - start
    return start, end
