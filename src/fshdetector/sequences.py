"""Small sequence and FASTA helpers."""

from collections.abc import Iterable
import importlib.resources
import pathlib

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(sequence: str) -> str:
    """Returns the reverse complement of a DNA sequence."""
    return sequence.translate(_COMPLEMENT)[::-1]


def data_path(name: str) -> pathlib.Path:
    """Returns the path of a file shipped inside the package."""
    return pathlib.Path(str(importlib.resources.files("fshdetector"))) / (
        "data/" + name
    )


def read_fasta(path: str | pathlib.Path) -> dict[str, str]:
    """Reads a small FASTA file into a name to sequence mapping."""
    records: dict[str, list[str]] = {}
    name = None
    with open(path, encoding="ascii") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = []
            elif line:
                if name is None:
                    raise ValueError(f"{path} is not a FASTA file")
                records[name].append(line.upper())
    return {key: "".join(parts) for key, parts in records.items()}


def write_fasta(
    path: str | pathlib.Path, records: Iterable[tuple[str, str]]
) -> None:
    """Writes (name, sequence) pairs as FASTA with one line per sequence."""
    with open(path, "w", encoding="ascii") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n{sequence}\n")
