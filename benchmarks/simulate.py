"""Measures allele calls on simulated reads with known repeat counts.

Alleles are built from the repeat units and flanks of T2T-CHM13v2.0, reads are
drawn from them with ONT-like errors, and the whole analysis runs on the reads
(BLAST, unit barcodes, spanning calls, reassembly). Every reported count is
then compared with the truth.

Each sample draws its own error level, from about 0.4% (recent super-accuracy
basecalls) to about 5% (older chemistry or fast basecalls), with per-read
variation, a few very poor reads, strand-specific systematic errors, and
homopolymer errors. Some samples also carry what makes assembly go wrong:
chimeric reads, a mosaic contraction of one allele, units exchanged between
chr4 and chr10, and reads aligned to the wrong chromosome.

Outcomes per reported allele:

* exact call (spanning or reassembly): correct if it matches an unused true
  count, else wrong.
* lower bound (longest read or contig): violated if no true count reaches it.

Example:
    python benchmarks/simulate.py --ref-t2t T2T.fna --seeds 0:400 \
        --workers 40 --cache cache/ --out results.tsv
"""

import argparse
from collections.abc import Sequence
import concurrent.futures
import dataclasses
import importlib
import pathlib
import pickle
import tempfile

import numpy as np
import pysam

from fshdetector import alleles as alleles_module
from fshdetector import analysis
from fshdetector import barcode
from fshdetector import blast
from fshdetector import layout
from fshdetector import sequences as seqlib

CHROMOSOMES = layout.CHROMOSOMES
MODES = ("block", "shuffled", "uniform")
_FLANK = 40_000
# Error level of a sample relative to a base rate of about 1.1% per base.
_LEVELS = (0.4, 5.0)
_BAD_READS = 0.03  # Fraction of reads with four times the sample error.
_MOSAIC = 0.15  # Chromosomes with a third, contracted allele.
_EXCHANGE = 0.1  # Alleles holding units of the other chromosome.
_CHIMERA = 0.01  # Reads joined to the end of another read.
_BASES = np.frombuffer(b"ACGT", dtype=np.uint8)


@dataclasses.dataclass
class Library:
    """Repeat units and flanks of one chromosome, in array order."""

    proximal: str
    units: list[str]
    distal: str


def build_library(ref_t2t: str, threads: int) -> dict[str, Library]:
    """Cuts the T2T D4Z4 loci into flanks and full repeat units."""
    with pysam.FastaFile(ref_t2t) as fasta:
        windows = {
            r.chrom: fasta.fetch(r.chrom, r.start, r.end).upper()
            for r in layout.HUMAN.t2t_windows
        }
    with tempfile.TemporaryDirectory() as work:
        hits = blast.search(windows, 0, pathlib.Path(work), threads)
    consensus = next(
        iter(seqlib.read_fasta(seqlib.data_path(blast.UNIT)).values())
    )
    units = barcode.build_units(
        hits[blast.UNIT],
        {name: len(s) for name, s in windows.items()},
        consensus,
    )
    result = {}
    for chromosome, sequence in windows.items():
        read = units[chromosome]
        if read.reverse:
            sequence = seqlib.reverse_complement(sequence)
        found = [
            u
            for u in read.units
            if u is not None and u.length > barcode.COUNTED_LENGTH
        ]
        full = [
            sequence[u.start : u.end]
            for u in found
            if 3100 <= u.end - u.start <= 3500
        ]
        result[chromosome] = Library(
            sequence[max(0, found[0].start - _FLANK) : found[0].start],
            full,
            sequence[found[-1].end : found[-1].end + _FLANK],
        )
    return result


def _mutate(unit: str, count: int, rng: np.random.Generator) -> str:
    bases = list(unit)
    for position in rng.integers(0, len(bases), size=count):
        choices = [b for b in "ACGT" if b != bases[position]]
        bases[position] = choices[rng.integers(3)]
    return "".join(bases)


def make_array(
    library: Library, count: int, mode: str, rng: np.random.Generator
) -> list[str]:
    """Returns the units of one allele.

    Modes set how alike the units are: a contiguous block of reference units,
    reference units in random order, or one unit repeated. Each allele also
    gets private substitutions.
    """
    inner = library.units[1:-1]
    size = max(0, count - 2)
    if mode == "block":
        start = int(rng.integers(len(inner)))
        chosen = [inner[(start + i) % len(inner)] for i in range(size)]
    elif mode == "shuffled":
        chosen = [inner[i] for i in rng.integers(len(inner), size=size)]
    else:
        chosen = [inner[int(rng.integers(len(inner)))]] * size
    units = ([library.units[0]] + chosen + [library.units[-1]])[:count]
    private = rng.choice([0, 1, 3, 8])
    return [_mutate(u, int(rng.poisson(private)), rng) for u in units]


def add_errors(
    sequence: str,
    rng: np.random.Generator,
    systematic: np.ndarray,
    level: float,
) -> str:
    """Applies random, homopolymer, and position-specific read errors.

    Args:
        sequence: error-free read on the array strand.
        rng: random generator.
        systematic: read positions that are often misread on this strand.
        level: error rate of the sample relative to the base rate.
    """
    scale = level * rng.lognormal(0.0, 0.4)
    if rng.random() < _BAD_READS:
        scale *= 4
    bases = np.frombuffer(sequence.encode(), dtype=np.uint8).copy()
    n = len(bases)
    substitute = rng.random(n) < min(0.1, 0.004 * scale)
    often = min(0.4, 0.15 * np.sqrt(level))
    substitute[systematic[rng.random(len(systematic)) < often]] = True
    index = np.flatnonzero(substitute)
    shift = rng.integers(1, 4, size=len(index))
    codes = np.searchsorted(_BASES, bases[index])
    bases[index] = _BASES[(codes + shift) % 4]
    keep = rng.random(n) >= min(0.1, 0.004 * scale)
    run = (
        (bases[3:] == bases[2:-1])
        & (bases[2:-1] == bases[1:-2])
        & (bases[1:-2] == bases[:-3])
    )
    ends = np.flatnonzero(run) + 3
    keep[ends[rng.random(len(ends)) < min(0.4, 0.06 * scale)]] = False
    repeats = keep.astype(np.int64)
    inserted = rng.random(n) < min(0.08, 0.003 * scale)
    inserted[ends[rng.random(len(ends)) < min(0.3, 0.04 * scale)]] = True
    repeats[inserted & keep] += 1
    out = np.repeat(bases, repeats)
    extra = np.flatnonzero(np.repeat(inserted & keep, repeats))[1::2]
    out[extra] = _BASES[rng.integers(0, 4, size=len(extra))]
    return out.tobytes().decode()


@dataclasses.dataclass
class Sample:
    """A simulated sample and its truth.

    ``truth`` lists the unit counts of the two alleles of each chromosome,
    followed by the count of a mosaic allele if there is one. ``origin`` maps a
    read to its allele and to the first array unit that it covers by at least
    300 bases.
    """

    reads: dict[str, str]
    aligned: dict[str, str]
    truth: dict[str, list[int]]
    modes: dict[str, list[str]]
    spanned: dict[str, list[bool]]
    coverage: float
    level: float
    origin: dict[str, tuple[int, int]]
    chimeras: set[str]


def _design(
    libraries: dict[str, Library], chromosome: str, rng: np.random.Generator
) -> list[tuple[list[str], str, float]]:
    """Returns the arrays of a chromosome as (units, mode, read share)."""
    other = CHROMOSOMES[1 - CHROMOSOMES.index(chromosome)]
    arrays = []
    for _ in range(2):
        count = int(
            rng.integers(2, 12) if rng.random() < 0.4 else rng.integers(8, 31)
        )
        mode = MODES[int(rng.integers(len(MODES)))]
        units = make_array(libraries[chromosome], count, mode, rng)
        if rng.random() < _EXCHANGE and count >= 4:
            size = int(rng.integers(2, count - 1))
            start = int(rng.integers(1, count - size))
            donor = make_array(libraries[other], size + 2, "block", rng)
            units[start : start + size] = donor[1:-1]
        arrays.append((units, mode, 1.0))
    if rng.random() < _MOSAIC:
        units, mode, _ = arrays[int(rng.integers(2))]
        if len(units) >= 4:
            lost = int(rng.integers(1, len(units) - 2))
            start = int(rng.integers(1, len(units) - lost))
            share = float(rng.uniform(0.15, 0.5))
            arrays.append((units[:start] + units[start + lost :], mode, share))
    return arrays


def simulate(libraries: dict[str, Library], seed: int) -> Sample:
    """Draws the alleles of a sample and reads from them."""
    rng = np.random.default_rng(seed)
    coverage = float(rng.uniform(8, 30))
    median = float(rng.uniform(14_000, 30_000))
    level = float(np.exp(rng.uniform(*np.log(_LEVELS))))
    reads: dict[str, str] = {}
    aligned: dict[str, str] = {}
    truth: dict[str, list[int]] = {}
    modes: dict[str, list[str]] = {}
    spanned: dict[str, list[bool]] = {}
    origin: dict[str, tuple[int, int]] = {}
    for chromosome in CHROMOSOMES:
        library = libraries[chromosome]
        truth[chromosome] = []
        modes[chromosome] = []
        spanned[chromosome] = []
        for allele, (units, mode, share) in enumerate(
            _design(libraries, chromosome, rng)
        ):
            array = "".join(units)
            template = library.proximal + array + library.distal
            starts = np.cumsum(
                [len(library.proximal)] + [len(u) for u in units[:-1]]
            )
            # Positions that one strand tends to misread, repeated per unit.
            offsets = rng.integers(0, 3000, size=(2, 6))
            bias = [
                np.sort((starts[:, None] + offsets[s]).ravel()) for s in (0, 1)
            ]
            array_start = len(library.proximal)
            array_end = array_start + len(array)
            wanted = share * coverage * len(array)
            drawn = 0
            number = 0
            full = False
            while drawn < wanted:
                length = int(
                    np.clip(rng.lognormal(np.log(median), 0.8), 3000, 150_000)
                )
                start = int(rng.integers(-length, len(template)))
                low, high = max(0, start), min(len(template), start + length)
                overlap = min(high, array_end) - max(low, array_start)
                if overlap < 1500:
                    continue
                drawn += overlap
                full |= low <= array_start - 5000 and high >= array_end + 500
                strand = int(rng.integers(2))
                positions = bias[strand]
                positions = (
                    positions[(positions >= low) & (positions < high)] - low
                )
                read = add_errors(template[low:high], rng, positions, level)
                if strand:
                    read = seqlib.reverse_complement(read)
                name = f"{chromosome}_{allele}_{number}"
                number += 1
                reads[name] = read
                origin[name] = (
                    allele,
                    next(
                        i
                        for i, (s, u) in enumerate(
                            zip(starts, units, strict=True)
                        )
                        if min(high, s + len(u)) - max(low, s) >= 300
                    ),
                )
                wrong = rng.random() < 0.02
                aligned[name] = (
                    CHROMOSOMES[1 - CHROMOSOMES.index(chromosome)]
                    if wrong
                    else chromosome
                )
            truth[chromosome].append(len(units))
            modes[chromosome].append(mode)
            spanned[chromosome].append(full)
    names = list(reads)
    chimeras = set()
    for name in names:
        if rng.random() < _CHIMERA:
            chimeras.add(name)
            partner = reads[names[int(rng.integers(len(names)))]]
            cut = int(rng.integers(len(reads[name]) // 4, len(reads[name])))
            reads[name] = (
                reads[name][:cut]
                + partner[int(rng.integers(len(partner) // 2)) :]
            )
    return Sample(
        reads,
        aligned,
        truth,
        modes,
        spanned,
        coverage,
        level,
        origin,
        chimeras,
    )


def judge(
    truth: Sequence[int], called: Sequence[alleles_module.Allele]
) -> list[dict]:
    """Compares the two reported alleles of a chromosome with the truth."""
    unused = list(truth)
    rows = []
    # Exact calls claim true alleles first, spanning before reassembly.
    for allele in sorted(
        called, key=lambda a: a.source != alleles_module.SPANNING
    ):
        if allele.source in (
            alleles_module.SPANNING,
            alleles_module.REASSEMBLY,
        ):
            correct = allele.ru in unused
            if correct:
                unused.remove(allele.ru)
            rows.append(
                {
                    "source": allele.source,
                    "reported": allele.ru,
                    "outcome": "correct" if correct else "wrong",
                    "spanning_reads": allele.spanning_reads,
                }
            )
    bounds = sorted(
        (a for a in called if a.source == alleles_module.LONGEST),
        key=lambda a: -a.ru,
    )
    for allele in bounds:
        fitting = sorted(n for n in unused if n >= allele.ru)
        if fitting:
            unused.remove(fitting[0])
        rows.append(
            {
                "source": allele.source,
                "reported": allele.ru,
                "outcome": "valid_bound" if fitting else "violated_bound",
            }
        )
    for allele in called:
        bound = allele.contig_lower_bound
        if bound is not None:
            rows.append(
                {
                    "source": "contig_bound",
                    "reported": bound,
                    "outcome": (
                        "valid_bound"
                        if bound <= max(truth)
                        else "violated_bound"
                    ),
                }
            )
    return rows


def load(
    seed: int, libraries: dict[str, Library], cache: pathlib.Path
) -> tuple[Sample, dict[str, list[blast.Hit]]]:
    """Simulates a sample and runs BLAST, or reloads both from the cache."""
    path = cache / f"sim{seed}.pkl"
    if path.is_file():
        with open(path, "rb") as handle:
            fields, hits = pickle.load(handle)
        return Sample(**fields), hits
    sample = simulate(libraries, seed)
    with tempfile.TemporaryDirectory() as work:
        hits = blast.search(sample.reads, 0, pathlib.Path(work), 1)
    with open(path, "wb") as handle:
        pickle.dump((vars(sample), hits), handle)
    return sample, hits


def evaluate(
    seed: int, libraries: dict[str, Library], cache: pathlib.Path
) -> list[dict]:
    """Judges the calls of one simulated sample."""
    sample, hits = load(seed, libraries, cache)
    reads = analysis.analyze(sample.reads, sample.aligned, hits)
    called = alleles_module.call_alleles(
        reads.model, reads.calls, reads.fragments
    )
    error = float(np.median(reads.model.error)) if len(reads.sites) else np.nan
    rows = []
    for chromosome in CHROMOSOMES:
        own = [a for a in called if a.chromosome == chromosome]
        truth = sample.truth[chromosome]
        hidden = [
            n
            for n, full in zip(
                truth[:2], sample.spanned[chromosome][:2], strict=True
            )
            if not full
        ]
        solved = [a.ru for a in own if a.source == alleles_module.REASSEMBLY]
        shared = {
            "seed": seed,
            "chromosome": chromosome,
            "truth": ",".join(map(str, truth)),
            "modes": ",".join(sample.modes[chromosome]),
            "coverage": round(sample.coverage, 1),
            "error_level": round(sample.level, 2),
            "read_error": round(error, 4),
            "spanning_reads": 0,
        }
        rows += [shared | row for row in judge(truth, own)]
        for count in hidden:
            found = count in solved
            if found:
                solved.remove(count)
            rows.append(
                shared
                | {
                    "source": "unspanned_allele",
                    "reported": count,
                    "outcome": "reassembled" if found else "not_reassembled",
                }
            )
    return rows


def summarize(rows: list[dict]) -> str:
    counts: dict[tuple[str, str], int] = {}
    for row in rows:
        key = (row["source"], row["outcome"])
        counts[key] = counts.get(key, 0) + 1
    for row in rows:
        level = row["error_level"]
        stratum = "<1" if level < 1 else "1-2.5" if level < 2.5 else ">=2.5"
        key = (f"{row['source']}[level {stratum}]", row["outcome"])
        counts[key] = counts.get(key, 0) + 1
    return "\n".join(
        f"{source}\t{outcome}\t{n}"
        for (source, outcome), n in sorted(counts.items())
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--ref-t2t", required=True)
    parser.add_argument("--seeds", default="0:100", help="first:last")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--cache", required=True, type=pathlib.Path)
    parser.add_argument("--out", required=True, type=pathlib.Path)
    parser.add_argument(
        "--set",
        action="append",
        default=[],
        metavar="module.NAME=value",
        help="override a threshold, for example assemble.BRANCH_READS=2",
    )
    args = parser.parse_args()
    first, last = map(int, args.seeds.split(":"))
    args.cache.mkdir(parents=True, exist_ok=True)
    for setting in args.set:
        target, value = setting.split("=")
        module, name = target.split(".")
        module = importlib.import_module(f"fshdetector.{module}")
        if not hasattr(module, name):
            parser.error(f"unknown threshold: {target}")
        setattr(module, name, type(getattr(module, name))(value))
    libraries = build_library(args.ref_t2t, args.workers)
    for name, library in libraries.items():
        print(f"{name}: {len(library.units)} full reference units")
    rows = []
    with concurrent.futures.ProcessPoolExecutor(args.workers) as pool:
        jobs = [
            pool.submit(evaluate, s, libraries, args.cache)
            for s in range(first, last)
        ]
        for done, job in enumerate(concurrent.futures.as_completed(jobs), 1):
            rows += job.result()
            if done % 100 == 0:
                print(f"{done} samples done", flush=True)
    rows.sort(key=lambda r: (r["seed"], r["chromosome"], r["source"]))
    columns = list(rows[0])
    with open(args.out, "w", encoding="ascii") as out:
        out.write("\t".join(columns) + "\n")
        for row in rows:
            out.write("\t".join(str(row[c]) for c in columns) + "\n")
    print(summarize(rows))


if __name__ == "__main__":
    main()
