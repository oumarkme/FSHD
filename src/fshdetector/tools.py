"""Runs the external programs (BLAST+, minimap2, modkit)."""

from collections.abc import Sequence
import shutil
import subprocess


class ToolError(RuntimeError):
    """An external program is missing or exited with an error."""


def require(*names: str) -> None:
    """Raises ToolError if a program is not on PATH."""
    missing = [name for name in names if shutil.which(name) is None]
    if missing:
        raise ToolError(f"required programs not found: {', '.join(missing)}")


def run(command: Sequence[str], stdout_path: str | None = None) -> str:
    """Runs a command and fails loudly if it fails.

    Args:
        command: program and arguments, without a shell.
        stdout_path: if given, standard output is written to this file.

    Returns:
        Standard output, or an empty string when it went to a file.

    Raises:
        ToolError: the program exited with a non-zero status.
    """
    command = [str(part) for part in command]
    if stdout_path is None:
        done = subprocess.run(
            command, capture_output=True, text=True, check=False
        )
        output = done.stdout
    else:
        with open(stdout_path, "wb") as handle:
            done = subprocess.run(
                command, stdout=handle, stderr=subprocess.PIPE, check=False
            )
        output = ""
    if done.returncode != 0:
        error = done.stderr
        if isinstance(error, bytes):
            error = error.decode(errors="replace")
        raise ToolError(
            f"{command[0]} failed with status {done.returncode}:\n"
            f"{error[-3000:]}"
        )
    return output


def versions() -> dict[str, str]:
    """Returns the first version line of each external program."""
    found = {}
    for name, flag in (
        ("blastn", "-version"),
        ("minimap2", "--version"),
        ("modkit", "--version"),
    ):
        if shutil.which(name):
            found[name] = run([name, flag]).strip().splitlines()[0]
    return found
