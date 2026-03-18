#!/bin/bash
set -euo pipefail

# GS4 server environment
CONDA_DIR="/home/user/miniforge-pypy3/envs/mark"
export PATH="${CONDA_DIR}/bin:$PATH"
export JAVA_HOME="${CONDA_DIR}"
export JAVA_CMD="${CONDA_DIR}/bin/java"

# Require input BAM/CRAM file
if [[ $# -lt 1 ]]; then
    echo "Usage: bash run.sh <input.bam|cram> [extra nextflow args...]"
    echo ""
    echo "Examples:"
    echo "  bash run.sh /raid/data/sample1.bam"
    echo "  bash run.sh /raid/data/sample1.bam -resume"
    exit 1
fi

INPUT="$1"
shift

# Pipeline directory (where this script lives)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Build SIF if environment.def is newer or SIF doesn't exist
SIF="$SCRIPT_DIR/environment.sif"
DEF="$SCRIPT_DIR/environment.def"
if [[ ! -f "$SIF" || "$DEF" -nt "$SIF" ]]; then
    echo "Building container: $SIF"
    apptainer build --force "$SIF" "$DEF"
fi

# Clean old Nextflow log files
rm -f "$SCRIPT_DIR"/*log*

# Run pipeline
nextflow -log "$SCRIPT_DIR/nf-fshd.log" run "$SCRIPT_DIR/main.nf" \
    --input "$INPUT" \
    -resume \
    "$@"
