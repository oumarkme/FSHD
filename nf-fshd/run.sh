#!/bin/bash
set -euo pipefail

# Detect environment and set conda PATH + JAVA_HOME
# WSL has /home/mark, GS4 has /home/user
if [[ -d "/home/mark/miniforge3" ]]; then
    CONDA_DIR="/home/mark/miniforge3/envs/mark"
    PROFILE="wsl"
else
    CONDA_DIR="/home/user/miniforge-pypy3/envs/mark"
    PROFILE="gs4"
fi

export PATH="${CONDA_DIR}/bin:$PATH"
export JAVA_HOME="${CONDA_DIR}"
export JAVA_CMD="${CONDA_DIR}/bin/java"

# Require input file as argument
if [[ $# -lt 1 ]]; then
    echo "Usage: bash run.sh <input.bam|cram> [extra nextflow args...]"
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
    -profile "$PROFILE" \
    --input "$INPUT" \
    -resume \
    "$@"
