#!/usr/bin/env bash
# Verbose inner runner that logs EVERYTHING and keeps tmux pane open on failure.
set -Eeuo pipefail

CONDA_ENV="${1:-docking}"
PROJECT_DIR="${2:-$PWD}"
PORT="${3:-5050}"
LOGDIR="$PROJECT_DIR/logs"
LOG="$LOGDIR/server.out"

mkdir -p "$LOGDIR"

{
  echo "========== $(date '+%Y-%m-%d %H:%M:%S') =========="
  echo "[INFO] Starting inner runner"
  echo "[INFO] CONDA_ENV=$CONDA_ENV"
  echo "[INFO] PROJECT_DIR=$PROJECT_DIR"
  echo "[INFO] PORT=$PORT"
  echo "[INFO] SHELL=$SHELL"

  # --- Conda init (try common locations) ---
  conda_init() {
    local c1="$HOME/miniconda3/etc/profile.d/conda.sh"
    local c2="$HOME/anaconda3/etc/profile.d/conda.sh"
    local c3="/opt/conda/etc/profile.d/conda.sh"
    if [[ -f "$c1" ]]; then source "$c1"; echo "[INFO] Sourced $c1"
    elif [[ -f "$c2" ]]; then source "$c2"; echo "[INFO] Sourced $c2"
    elif [[ -f "$c3" ]]; then source "$c3"; echo "[INFO] Sourced $c3"
    elif command -v conda >/dev/null 2>&1; then
      echo "[WARN] conda on PATH but conda.sh not found. Proceeding anyway."
    else
      echo "[ERROR] conda not found; install or adjust path."
      return 20
    fi
  }

  if ! conda_init; then
    echo "[FATAL] Failed to initialize conda."
    exit 20
  fi

  echo "[INFO] conda: $(command -v conda || echo 'NOT FOUND')"
  echo "[INFO] Activating env: $CONDA_ENV"
  if ! conda activate "$CONDA_ENV"; then
    echo "[FATAL] conda activate $CONDA_ENV failed."
    exit 21
  fi

  cd "$PROJECT_DIR"
  echo "[INFO] CWD=$(pwd)"
  echo "[INFO] Python: $(command -v python || true)"; python -V || true
  echo "[INFO] Gunicorn: $(command -v gunicorn || echo 'NOT FOUND')"

  export PORT="$PORT"
  export PYTHONUNBUFFERED=1

  if command -v gunicorn >/dev/null 2>&1; then
    echo "[INFO] Launching gunicorn 0.0.0.0:$PORT ..."
    # Do NOT exec; if it dies we want to keep the shell alive.
    gunicorn "app:app" -b 0.0.0.0:"$PORT" --workers 2 --threads 4 --timeout 120
  else
    echo "[WARN] gunicorn not found — running python app.py ..."
    python app.py
  fi
  rc=$?
  echo "[ERROR] Server process exited with code $rc"
  echo "[INFO] Keeping pane open. Check logs above."
  # Keep pane open so tmux session doesn’t die; infinite sleep.
  while true; do sleep 3600; done
} >>"$LOG" 2>&1
