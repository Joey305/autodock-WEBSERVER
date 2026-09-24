#!/usr/bin/env bash
# Explicit Pegasus/LSF smoke test.  This does not run as part of normal docking.
set -euo pipefail

PROJECT="${PROJECT:-brd}"
QUEUE="${QUEUE:-gpu_cheminfo}"
LOG_DIR="${LOG_DIR:-logs}"
mkdir -p "$LOG_DIR"

bsub \
  -P "$PROJECT" \
  -q "$QUEUE" \
  -n 1 \
  -W 00:10 \
  -o "$LOG_DIR/gpu_cheminfo_test_%J.out" \
  -e "$LOG_DIR/gpu_cheminfo_test_%J.err" \
  'echo "HOST=$(hostname)"; echo "CPUS=${LSB_DJOB_NUMPROC:-unknown}"; sleep 5'
