#!/usr/bin/env bash
# Explicit Pegasus normal-queue smoke test. This does not run as part of normal docking.
set -euo pipefail

PROJECT="${PROJECT:-brd}"
QUEUE="${QUEUE:-normal}"
LOG_DIR="${LOG_DIR:-logs}"
RESULTS_ROOT="${VINA_RESULTS_ROOT:-/projectnb/cheminfo/home/${USER}/AutoDockVina-Results/smoke-tests}"
SCRATCH_ROOT="${VINA_SCRATCH_ROOT:-/tmp}"
mkdir -p "$LOG_DIR"
export RESULTS_ROOT SCRATCH_ROOT

bsub \
  -P "$PROJECT" \
  -q "$QUEUE" \
  -n 1 \
  -W 00:10 \
  -o "$LOG_DIR/normal_queue_test_%J.out" \
  -e "$LOG_DIR/normal_queue_test_%J.err" \
  '
set -euo pipefail
probe="normal_queue_probe_${LSB_JOBID}"
results_probe="$RESULTS_ROOT/.${probe}"
scratch_probe="$(mktemp -d "${SCRATCH_ROOT%/}/${probe}_XXXXXX")"
cleanup() {
  rm -rf "$scratch_probe"
  rm -rf "$results_probe"
}
trap cleanup EXIT
mkdir "$results_probe"
printf "shared-results-write=ok\\n" > "$results_probe/proof.txt"
test -s "$results_probe/proof.txt"
printf "local-scratch-write=ok\\n" > "$scratch_probe/proof.txt"
test -s "$scratch_probe/proof.txt"
echo "HOST=$(hostname)"
echo "CPUS=${LSB_DJOB_NUMPROC:-unknown}"
df -h "$RESULTS_ROOT" "$SCRATCH_ROOT"
'
