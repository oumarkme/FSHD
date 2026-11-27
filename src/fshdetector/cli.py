"""Command line interface."""

import argparse
from collections.abc import Sequence
import logging
import os
import pathlib
import sys

from fshdetector import __version__
from fshdetector import align
from fshdetector import pipeline
from fshdetector import reference
from fshdetector import tools

_DEFAULT_THREADS = 8


def _existing(value: str) -> pathlib.Path:
    path = pathlib.Path(value)
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"file not found: {value}")
    return path


def _positive(value: str) -> int:
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return number


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="fshdetector",
        description="D4Z4 repeat counting, haplotyping, and methylation "
        "from ONT long reads.",
    )
    parser.add_argument(
        "--version", action="version", version=f"fshdetector {__version__}"
    )
    commands = parser.add_subparsers(dest="command", required=True)

    run = commands.add_parser("run", help="analyze one sample")
    run.add_argument(
        "--input",
        required=True,
        type=_existing,
        help="indexed BAM or CRAM aligned to hg38 or T2T-CHM13v2.0",
    )
    run.add_argument(
        "--ref-t2t",
        required=True,
        type=_existing,
        help="T2T-CHM13v2.0 FASTA with chr names and a .fai index",
    )
    run.add_argument(
        "--ref-hg38",
        type=_existing,
        help="hg38 FASTA; needed for an hg38 CRAM",
    )
    run.add_argument(
        "--t2t-index",
        type=_existing,
        help="minimap2 index from 'fshdetector index' (hg38 input only)",
    )
    run.add_argument("--build", choices=("auto", "hg38", "t2t"), default="auto")
    run.add_argument("--outdir", type=pathlib.Path, default="fshd_results")
    run.add_argument("--sample", help="sample name; default: input file name")
    run.add_argument("--threads", type=_positive, default=_DEFAULT_THREADS)
    run.add_argument(
        "--tmpdir", type=pathlib.Path, help="scratch; default: $TMPDIR"
    )
    run.add_argument(
        "--validate-reassembly",
        action="store_true",
        help="also test reassembly against the spanning reads",
    )

    index = commands.add_parser(
        "index", help="build the minimap2 index used for hg38 input"
    )
    index.add_argument("--ref-t2t", required=True, type=_existing)
    index.add_argument("--output", required=True, type=pathlib.Path)
    index.add_argument("--threads", type=_positive, default=_DEFAULT_THREADS)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        stream=sys.stderr,
    )
    # htslib must not download reference sequences from the internet.
    os.environ.setdefault("REF_PATH", os.devnull)
    try:
        if args.command == "index":
            align.build_index(args.ref_t2t, args.output, args.threads)
            return 0
        pipeline.run(
            pipeline.Options(
                input=args.input,
                ref_t2t=args.ref_t2t,
                outdir=args.outdir,
                ref_hg38=args.ref_hg38,
                t2t_index=args.t2t_index,
                build=args.build,
                sample=args.sample,
                threads=args.threads,
                tmpdir=args.tmpdir,
                validate_reassembly=args.validate_reassembly,
            )
        )
    except (tools.ToolError, reference.ReferenceError_, ValueError,
            OSError) as error:  # fmt: skip
        logging.getLogger("fshdetector").error("%s", error)
        return 1
    return 0
