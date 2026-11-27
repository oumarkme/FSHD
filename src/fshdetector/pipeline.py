"""Runs one sample from an aligned BAM or CRAM to the result tables."""

from collections.abc import Iterator, Sequence
import contextlib
import dataclasses
import datetime
import json
import logging
import pathlib
import re
import shutil
import tempfile
import time

import pysam

from fshdetector import __version__
from fshdetector import _screen
from fshdetector import align
from fshdetector import alleles as alleles_module
from fshdetector import analysis
from fshdetector import blast
from fshdetector import extract
from fshdetector import layout as layout_module
from fshdetector import methylation
from fshdetector import reference
from fshdetector import report
from fshdetector import sequences as seqlib
from fshdetector import tools
from fshdetector import validate
from fshdetector import validation

LOG = logging.getLogger("fshdetector")
_DISTAL_UNITS = 3
_PARALLEL_REQUESTS = 8
_MODEL = re.compile(r"basecall_model=(\S+)")


@dataclasses.dataclass
class Options:
    """Settings of one run."""

    input: pathlib.Path
    ref_t2t: pathlib.Path
    outdir: pathlib.Path
    ref_hg38: pathlib.Path | None = None
    t2t_index: pathlib.Path | None = None
    build: str = "auto"
    sample: str | None = None
    threads: int = 8
    tmpdir: pathlib.Path | None = None
    validate_reassembly: bool = False

    @property
    def sample_name(self) -> str:
        return self.sample or self.input.stem


class _Timer:

    def __init__(self):
        self.seconds: dict[str, float] = {}

    @contextlib.contextmanager
    def step(self, name: str) -> Iterator[None]:
        start = time.monotonic()
        yield
        self.seconds[name] = round(time.monotonic() - start, 2)
        LOG.info("%s: %.1f s", name, self.seconds[name])


def basecall_models(header: pysam.AlignmentHeader) -> list[str]:
    """Returns the basecaller models named in the read groups."""
    found = set()
    for group in header.to_dict().get("RG", []):
        match = _MODEL.search(group.get("DS", ""))
        if match:
            found.add(match.group(1))
    return sorted(found)


def _candidate_records(
    reads: extract.Reads, names: set[str]
) -> list[pysam.AlignedSegment]:
    records = [r for r in reads.window_records if r.query_name in names]
    inside = {
        (r.query_name, r.flag, r.reference_id, r.reference_start)
        for r in records
    }
    for name in names:
        primary = reads.primaries[name]
        key = (name, primary.flag, primary.reference_id,
               primary.reference_start)  # fmt: skip
        if key not in inside:
            records.append(primary)
    return records


def _unit_methylation(
    allele: alleles_module.Allele,
    sample_reads: analysis.SampleReads,
    primaries: dict[str, pysam.AlignedSegment],
    threshold: float,
) -> dict[tuple, list]:
    """5mC calls per repeat unit of one allele, keyed by unit position.

    Returns:
        (unit from proximal end or None, unit from distal end or None) to
        [methylated calls, all calls, read names].
    """
    offsets = {}
    if allele.template is not None:
        offsets = {f.read_id: o for f, o in allele.template.placed}
    length = allele.ru if allele.source != alleles_module.LONGEST else None
    profile: dict[tuple, list] = {}
    for read_id in allele.read_ids:
        fragment = sample_reads.fragments[read_id]
        if not (fragment.head or fragment.tail or read_id in offsets):
            continue
        read = sample_reads.units[read_id]
        calls = methylation.read_cpg_calls(primaries[read_id], threshold)
        for row, unit in enumerate(sample_reads.placed_units[read_id]):
            if unit is None:
                continue
            start = distal = None
            if read_id in offsets and length is not None:
                start = row + offsets[read_id] + 1
            elif fragment.head:
                start = row + 1
            if fragment.tail:
                distal = fragment.units - row
            if length is not None:
                if start is None and distal is not None:
                    start = length - distal + 1
                if start is not None:
                    distal = length - start + 1
                    if not 1 <= start <= length:
                        continue
            if start is None and distal is None:
                continue
            low, high = methylation.forward_interval(
                read.read_length, unit.start, unit.end, read.reverse
            )
            modified, total = methylation.interval_calls(calls, low, high)
            if not total:
                continue
            entry = profile.setdefault((start, distal), [0, 0, set()])
            entry[0] += modified
            entry[1] += total
            entry[2].add(read_id)
    return profile


def _methylation(
    options: Options,
    layout: layout_module.Layout,
    reads: extract.Reads,
    alleles: Sequence[alleles_module.Allele],
    sample_reads: analysis.SampleReads,
    workdir: pathlib.Path,
):
    """Region methylation by upstream tag and by allele, and unit profiles."""
    pileup = methylation.Pileup(
        reads.window_records, options.ref_t2t, layout.meth_targets, workdir
    )
    threshold = pileup.estimate_threshold(options.threads)
    upstream = pileup.run("upstream", threshold, options.threads)

    numbers = []
    seen: dict[str, int] = {}
    groups: dict[str, int] = {}
    for allele in alleles:
        seen[allele.chromosome] = seen.get(allele.chromosome, 0) + 1
        numbers.append(seen[allele.chromosome])
        for read_id in allele.read_ids:
            groups[read_id] = numbers[-1]
    by_allele = pileup.run("alleles", threshold, options.threads, groups)

    values = []
    profile_rows = []
    in_pileup = {r.query_name: r.reference_name for r in pileup.records}
    for allele, number in zip(alleles, numbers, strict=True):
        row: dict[str, object] = {}
        sites = by_allele[f"hp{number}"]
        covered = [
            r for r in allele.read_ids if in_pileup.get(r) == allele.chromosome
        ]
        for region in layout.meth_regions:
            if region.chrom != allele.chromosome or not covered:
                continue
            summary = methylation.summarize(sites, region)
            row[f"{region.name}_methylation"] = summary.percent
            row[f"{region.name}_CpG"] = summary.n_cpg
        if threshold is not None and allele.read_ids:
            profile = _unit_methylation(
                allele, sample_reads, reads.primaries, threshold
            )
            distal = [
                v
                for (_, d), v in profile.items()
                if d is not None and d <= _DISTAL_UNITS
            ]
            total = sum(v[1] for v in distal)
            if total:
                row["distal_units_methylation"] = round(
                    100 * sum(v[0] for v in distal) / total, 2
                )
                row["distal_units_calls"] = total
                row["distal_units_reads"] = len(
                    set().union(*(v[2] for v in distal))
                )
            for (start, end), (modified, calls, names) in sorted(
                profile.items(),
                key=lambda item: (item[0][0] or 10**6, -(item[0][1] or 0)),
            ):
                profile_rows.append(
                    (
                        allele.chromosome,
                        number,
                        start,
                        end,
                        round(100 * modified / calls, 2),
                        calls,
                        len(names),
                    )
                )
        values.append(row)
    return threshold, upstream, values, profile_rows, numbers


def _collect(
    path: pathlib.Path,
    reference_path: pathlib.Path | None,
    windows: Sequence[layout_module.Region],
    screen,
    threads: int,
) -> extract.Reads:
    with extract.alignment(path, reference_path, threads) as source:
        return extract.collect(
            source,
            windows,
            screen,
            lambda: extract.open_alignment(path, reference_path),
            min(_PARALLEL_REQUESTS, threads),
        )


def run(
    options: Options, layout: layout_module.Layout = layout_module.HUMAN
) -> pathlib.Path:
    """Analyzes one sample and returns the output directory."""
    started = datetime.datetime.now().astimezone()
    timer = _Timer()
    tools.require("blastn", "makeblastdb", "modkit")
    reference.check_t2t_reference(options.ref_t2t, layout)
    sample = options.sample_name
    warnings: list[str] = []
    threads = max(1, options.threads)

    with pysam.AlignmentFile(str(options.input), check_sq=False) as handle:
        header = handle.header
    build = options.build
    if build == "auto":
        lengths = dict(zip(header.references, header.lengths, strict=True))
        build = reference.detect_build(lengths, layout)
    LOG.info("sample %s, input build %s", sample, build)
    models = basecall_models(header)
    if models and not all("sup" in m for m in models):
        warnings.append(
            f"basecall model is not super-accuracy: {', '.join(models)}"
        )

    options.outdir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="fshdetector.", dir=options.tmpdir
    ) as scratch:
        work = pathlib.Path(scratch)
        results = work / "results"
        results.mkdir()

        with timer.step("build_screen"):
            queries = [
                sequence
                for name in blast.QUERIES
                for sequence in seqlib.read_fasta(
                    seqlib.data_path(name)
                ).values()
            ]
            screen = _screen.KmerScreen(
                queries,
                extra=reference.target_sequences(options.ref_t2t, layout),
            )

        with timer.step("read_input"):
            if build == "hg38":
                source_reference = options.ref_hg38
                windows = layout.hg38_windows
            else:
                source_reference = options.ref_t2t
                windows = layout.t2t_windows
            reads = _collect(
                options.input, source_reference, windows, screen, threads
            )
        counts = dict(reads.counts)

        if build == "hg38":
            tools.require("minimap2")
            with timer.step("realign_to_t2t"):
                bam = align.realign(
                    {n: reads.primaries[n] for n in reads.candidates},
                    options.t2t_index or options.ref_t2t,
                    work,
                    threads,
                )
                searched = reads.searched_bases
                reads = _collect(bam, None, layout.t2t_windows, screen, threads)
                reads.searched_bases = searched
            counts.update({f"t2t_{k}": v for k, v in reads.counts.items()})

        with timer.step("blast"):
            hits = blast.search(
                reads.candidates, reads.searched_bases, work, threads
            )

        with timer.step("annotate"):
            aligned = {
                name: record.reference_name
                for name, record in reads.primaries.items()
                if not record.is_unmapped
            }
            sample_reads = analysis.analyze(reads.candidates, aligned, hits)
        calls = sample_reads.calls
        if sample_reads.misplaced:
            warnings.append(
                f"{len(sample_reads.misplaced)} reads have an array-end marker "
                "on the wrong side; that end is ignored for allele calls"
            )

        with timer.step("alleles"):
            called = alleles_module.call_alleles(
                sample_reads.model, calls, sample_reads.fragments
            )

        with timer.step("methylation"):
            threshold, upstream, values, profile, numbers = _methylation(
                options, layout, reads, called, sample_reads, work
            )
            if threshold is None:
                warnings.append(
                    "no read with modified base tags covers the methylation "
                    "targets"
                )

        with timer.step("write_results"):
            prefix = results / sample
            report.write_reads(
                pathlib.Path(f"{prefix}.reads_annotated.tsv"), calls
            )
            report.write_repeats(pathlib.Path(f"{prefix}.repeats.tsv"), calls)
            report.write_methylation(
                pathlib.Path(f"{prefix}.methylation.tsv"), upstream, layout
            )
            report.write_alleles(
                pathlib.Path(f"{prefix}.alleles.tsv"), sample, called, values
            )
            report.write_table(
                pathlib.Path(f"{prefix}.unit_methylation.tsv"),
                (
                    "chromosome", "allele", "unit_from_proximal",
                    "unit_from_distal", "methylation", "calls", "reads",
                ),
                profile,
            )  # fmt: skip
            report.write_table(
                pathlib.Path(f"{prefix}.read_alleles.tsv"),
                ("read_id", "chromosome", "allele"),
                (
                    (read_id, allele.chromosome, number)
                    for allele, number in zip(called, numbers, strict=True)
                    for read_id in sorted(allele.read_ids)
                ),
            )
            pathlib.Path(f"{prefix}.validation.txt").write_text(
                validate.report(sample, calls), encoding="utf-8"
            )
            names = set(reads.candidates)
            records = _candidate_records(reads, names)
            extract.write_bam(f"{prefix}.d4z4.bam", reads.header, records)
            complete = {c.read_id for c in calls if c.complete}
            extract.write_bam(
                f"{prefix}.complete_reads.bam",
                reads.header,
                (r for r in records if r.query_name in complete),
            )
            shutil.copyfile(
                seqlib.data_path("fshd_feature_t2t.bed"),
                results / "fshd_feature_t2t.bed",
            )

        if options.validate_reassembly:
            with timer.step("validate_reassembly"):
                report.write_table(
                    pathlib.Path(f"{prefix}.reassembly_validation.tsv"),
                    validation.COLUMNS,
                    validation.run(
                        sample_reads.model,
                        called,
                        sample_reads.fragments,
                        sample_reads.partial_by_chromosome(),
                    ),
                )

        for warning in warnings:
            LOG.warning(warning)
        record = {
            "version": __version__,
            "sample": sample,
            "started": started.isoformat(timespec="seconds"),
            "input_build": build,
            "options": {
                key: str(value) if isinstance(value, pathlib.Path) else value
                for key, value in dataclasses.asdict(options).items()
            },
            "tools": tools.versions(),
            "basecall_models": models,
            "reads": counts,
            "calls": len(calls),
            "variant_sites": len(sample_reads.sites),
            "mean_site_error": (
                round(float(sample_reads.model.error.mean()), 5)
                if len(sample_reads.sites)
                else None
            ),
            "modkit_pass_threshold": threshold,
            "seconds": timer.seconds,
            "warnings": warnings,
        }
        pathlib.Path(f"{prefix}.run.json").write_text(
            json.dumps(record, indent=2) + "\n", encoding="utf-8"
        )
        for path in sorted(results.iterdir()):
            shutil.move(str(path), options.outdir / path.name)
    return options.outdir
