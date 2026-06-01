#!/bin/bash
# nf-fshd: Run the FSHD D4Z4 analysis pipeline on a single sample.
# Usage: bash nf-fshd/run.sh /path/to/sample.bam
set -euo pipefail

# --- Environment (GS4 server) ---
CONDA_DIR="/home/user/miniforge-pypy3/envs/mark"
export PATH="${CONDA_DIR}/bin:$PATH"
export JAVA_HOME="${CONDA_DIR}"
export JAVA_CMD="${CONDA_DIR}/bin/java"

# --- Input validation ---
if [[ $# -lt 1 ]]; then
    echo "Usage: bash run.sh <input.bam|cram>"
    echo "  e.g. bash run.sh /raid/data/sample.bam"
    exit 1
fi

INPUT="$1"
shift

# --- Pipeline directory ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Build container if needed ---
SIF="$SCRIPT_DIR/environment.sif"
DEF="$SCRIPT_DIR/environment.def"
if [[ ! -f "$SIF" || "$DEF" -nt "$SIF" ]]; then
    echo "Building container: $SIF"
    apptainer build --force "$SIF" "$DEF"
fi

# --- Run pipeline ---
rm -f "$SCRIPT_DIR"/*.log*
nextflow -log "$SCRIPT_DIR/nf-fshd.log" run "$SCRIPT_DIR/main.nf" \
    --input "$INPUT" \
    -resume \
    "$@"
