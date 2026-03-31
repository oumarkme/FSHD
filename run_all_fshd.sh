#!/bin/bash
set -euo pipefail

# =============================================================================
# Batch rerun: all FSHD patients through nf-fshd pipeline
# Run on GS4:  bash /lab_nas_home_user/Mark/FSHD/run_all_fshd.sh
# =============================================================================

PIPELINE="/lab_nas_home_user/Mark/FSHD/nf-fshd/run.sh"
LOGDIR="/lab_nas_home_user/Mark/FSHD/batch_logs"
mkdir -p "$LOGDIR"

CRAMS=(
  "/gene_lv/NGS_gene_lv_Storage/WGS/analysis/pipeline_output/WGS-ONT/WGS-ONT-A0002-7298388-D26478/WGS-ONT-A0002-7298388-D26478_sup_Q13_20250526115031/OUTPUT/1.ALIGNMENT/WGS-ONT-A0002-7298388-D26478.cram"
  "/gene_lv/NGS_gene_lv_Storage/WGS/analysis/pipeline_output/WGS-ONT/WGS-ONT-A0003-3616158-D26526/WGS-ONT-A0003-3616158-D26526_sup_Q13_20250526175003/OUTPUT/1.ALIGNMENT/WGS-ONT-A0003-3616158-D26526.cram"
  "/gene_lv/NGS_gene_lv_Storage/WGS/analysis/pipeline_output/WGS-ONT/WGS-ONT022-5676618-D25586/WGS-ONT022-5676618-D25586_sup_Q13_20250610180425/OUTPUT/1.ALIGNMENT/WGS-ONT022-5676618-D25586.cram"
  "/gene_lv/NGS_gene_lv_Storage/WGS/analysis/pipeline_output/WGS-ONT/WGS-ONT024-6777438-D25585/WGS-ONT024-6777438-D25585_sup_Q13_20250301181702/OUTPUT/1.ALIGNMENT/WGS-ONT024-6777438-D25585.cram"
  "/gene_lv/NGS_gene_lv_Storage/WGS/analysis/pipeline_output/WGS-ONT/WGS-ONT025-5676618-D25587/WGS-ONT025-5676618-D25587_sup_Q13_20250228104819/OUTPUT/1.ALIGNMENT/WGS-ONT025-5676618-D25587.cram"
  "/gene_lv/NGS_gene_lv_Storage/WGS/analysis/pipeline_output/WGS-ONT/WGS-ONT019-5657831-D23881/WGS-ONT019-5657831-D23881_sup_Q13_20250310120409/OUTPUT/1.ALIGNMENT/WGS-ONT019-5657831-D23881.cram"
  "/gene_lv/NGS_gene_lv_Storage/WGS/analysis/pipeline_output/WGS-ONT/WGS-ONT056-6553584-D29308/WGS-ONT056-6553584-D29308_sup_Q13_20260324125935/OUTPUT/1.ALIGNMENT/WGS-ONT056-6553584-D29308.cram"
  "/lab_nas_home_user/Mark/ONT/WGS-ONT045-6733790-D27741/WGS-ONT045-6733790-D27741.cram"
)

echo "=========================================="
echo " FSHD Batch Pipeline Run"
echo " $(date)"
echo "=========================================="
echo ""
echo "Skipped (no hg38 CRAM available):"
echo "  - WGS-ONT040-7336310-D27297"
echo "  - WGS-ONT-R0011-7311800-D27622"
echo "  - WGS-ONT-R0015-7327380-D27810"
echo ""

TOTAL=${#CRAMS[@]}
PASSED=0
FAILED=0

for i in "${!CRAMS[@]}"; do
  CRAM="${CRAMS[$i]}"
  SAMPLE=$(basename "$CRAM" .cram)
  N=$((i + 1))
  LOG="$LOGDIR/${SAMPLE}.log"

  echo "[$N/$TOTAL] $SAMPLE"
  echo "  CRAM: $CRAM"
  echo "  Log:  $LOG"

  if [[ ! -f "$CRAM" ]]; then
    echo "  SKIPPED: CRAM file not found"
    FAILED=$((FAILED + 1))
    continue
  fi

  if bash "$PIPELINE" "$CRAM" > "$LOG" 2>&1; then
    echo "  DONE ($(date))"
    PASSED=$((PASSED + 1))
  else
    echo "  FAILED (see log)"
    FAILED=$((FAILED + 1))
  fi
  echo ""
done

echo "=========================================="
echo " Batch complete: $PASSED passed, $FAILED failed out of $TOTAL"
echo " Output: /lab_nas_home_user/Mark/ONT/FSHD/"
echo " Logs:   $LOGDIR/"
echo "=========================================="
