"""BLAST search of candidate reads against the bundled query sequences."""

from collections.abc import Mapping
import concurrent.futures
import pathlib
from typing import NamedTuple

from fshdetector import sequences as seqlib
from fshdetector import tools

UNIT = "d4z4.fa"
MARKERS = "haplotype_markers.fa"
P13E11 = "P13E11.fasta"
QUERIES = (UNIT, MARKERS, P13E11)

_FORMAT = "6 std btop"
# BLAST keeps 500 subjects per query by default, which deep samples exceed.
_ALL_SUBJECTS = 10_000_000


class Hit(NamedTuple):
    """One BLAST alignment of a query (subject coordinates are 1-based).

    ``sstart > send`` means the read matches the reverse strand of the query.
    """

    query: str
    read_id: str
    identity: float
    length: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    btop: str


def parse(path: str | pathlib.Path, names: Mapping[str, str]) -> list[Hit]:
    """Parses tabular output and maps subject identifiers back to read names."""
    hits = []
    with open(path, encoding="ascii") as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            hits.append(
                Hit(
                    query=f[0],
                    read_id=names[f[1]],
                    identity=float(f[2]),
                    length=int(f[3]),
                    qstart=int(f[6]),
                    qend=int(f[7]),
                    sstart=int(f[8]),
                    send=int(f[9]),
                    evalue=float(f[10]),
                    bitscore=float(f[11]),
                    btop=f[12],
                )
            )
    return hits


def search(
    reads: Mapping[str, str],
    database_bases: int,
    workdir: pathlib.Path,
    threads: int,
) -> dict[str, list[Hit]]:
    """Runs the three searches and returns the hits of each query file.

    Args:
        reads: read name to sequence.
        database_bases: database size used for E-values, so that screening the
            reads does not change which hits are reported.
        workdir: scratch directory.
        threads: total threads shared by the searches.
    """
    if not reads:
        return {query: [] for query in QUERIES}
    names = {f"r{i}": name for i, name in enumerate(reads)}
    fasta = workdir / "candidates.fasta"
    seqlib.write_fasta(
        fasta, ((key, reads[name]) for key, name in names.items())
    )
    database = workdir / "blastdb"
    tools.run(
        ["makeblastdb", "-in", fasta, "-dbtype", "nucl", "-out", database]
    )
    total_bases = sum(len(s) for s in reads.values())
    database_bases = max(database_bases, total_bases)

    def one(query: str, query_threads: int) -> list[Hit]:
        output = workdir / f"{query}.tsv"
        tools.run(
            [
                "blastn",
                "-db", database,
                "-query", seqlib.data_path(query),
                "-out", output,
                "-num_threads", query_threads,
                "-outfmt", _FORMAT,
                "-dbsize", database_bases,
                "-max_target_seqs", _ALL_SUBJECTS,
            ]
        )  # fmt: skip
        return parse(output, names)

    # The repeat unit search dominates, so it gets the spare threads.
    share = {query: 1 for query in QUERIES}
    share[UNIT] = max(1, threads - 2)
    workers = len(QUERIES) if threads >= len(QUERIES) else 1
    if workers == 1:
        share = {query: threads for query in QUERIES}
    with concurrent.futures.ThreadPoolExecutor(workers) as pool:
        jobs = {q: pool.submit(one, q, share[q]) for q in QUERIES}
        return {query: job.result() for query, job in jobs.items()}
