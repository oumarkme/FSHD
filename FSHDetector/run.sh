#!/bin/bash
# FSHDetector: Run the FSHD D4Z4 analysis pipeline on a single sample.
#
# Usage:
#   bash FSHDetector/run.sh <input.bam|cram> [extra nextflow options]
#
# Provide the references (--ref_hg38 is required for CRAM; BAM may omit it):
#   bash FSHDetector/run.sh sample.bam \
#       --ref_hg38 /path/to/hg38.fa --ref_t2t /path/to/T2T-CHM13v2.0.fna
#
# This wrapper just adds -resume. You can equally call Nextflow directly:
#   nextflow run FSHDetector/main.nf --input sample.bam \
#       --ref_hg38 hg38.fa --ref_t2t t2t.fna
set -euo pipefail

# --- Optional site environment (maintainer only) ---
# On the maintainer's server the toolchain lives in a conda env that is not on
# the default PATH. If that env is present we use it; otherwise we rely on
# whatever Nextflow/Conda are already installed (the normal case). External
# users never need this directory.
GS4_ENV="/home/user/miniforge-pypy3/envs/mark"
if [[ -d "$GS4_ENV" ]]; then
    export PATH="${GS4_ENV}/bin:$PATH"
    export JAVA_HOME="${GS4_ENV}"
    export JAVA_CMD="${GS4_ENV}/bin/java"
fi

# --- Input validation ---
if [[ $# -lt 1 ]]; then
    echo "Usage: bash run.sh <input.bam|cram> [nextflow options]"
    echo "  e.g. bash run.sh sample.bam --ref_hg38 hg38.fa --ref_t2t t2t.fna"
    exit 1
fi

INPUT="$1"
shift

# --- Dependency check ---
# Nextflow runs the pipeline; Conda (Miniforge recommended) provides the tools.
if ! command -v nextflow >/dev/null 2>&1; then
    echo "ERROR: 'nextflow' is not on your PATH." >&2
    echo "  Install: https://www.nextflow.io/docs/latest/install.html" >&2
    exit 1
fi
if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: 'conda' is not on your PATH." >&2
    echo "  FSHDetector builds its tools with Conda. Install Miniforge:" >&2
    echo "  https://github.com/conda-forge/miniforge" >&2
    exit 1
fi

# --- Pipeline directory ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Run pipeline ---
nextflow -log "$SCRIPT_DIR/FSHDetector.log" run "$SCRIPT_DIR/main.nf" \
    --input "$INPUT" \
    -resume \
    "$@"
