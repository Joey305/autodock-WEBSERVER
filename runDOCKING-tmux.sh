#!/usr/bin/env bash
# Minimal + robust tmux launcher for DockingServer

set -e

SESSION_NAME="docking_flask"
CONDA_ENV="docking"
PROJECT_DIR="/mnt/c/Users/joeys/Documents/DockingServer"
PORT="5050"

# Two common conda locations (adjust if yours differs)
CONDA_SH_1="/home/joey305/miniconda3/etc/profile.d/conda.sh"
CONDA_SH_2="/home/joey305/anaconda3/etc/profile.d/conda.sh"

# Start a session with a login shell so the window exists and won't vanish
if ! tmux has-session -t "$SESSION_NAME" 2>/dev/null; then
  echo "Starting new tmux session: $SESSION_NAME"
  tmux new-session -d -s "$SESSION_NAME" "bash -l"
else
  echo "Session $SESSION_NAME already exists — restarting pane"
  tmux send-keys -t "$SESSION_NAME" C-c
fi

# Keep pane visible after process exit
tmux set-option -t "$SESSION_NAME" remain-on-exit on

# Run everything step-by-step, leaving a shell at the end
tmux send-keys -t "$SESSION_NAME" \
  "set -e" C-m \
  "if [ -f '$CONDA_SH_1' ]; then source '$CONDA_SH_1'; elif [ -f '$CONDA_SH_2' ]; then source '$CONDA_SH_2'; else echo '[FATAL] conda.sh not found'; fi" C-m \
  "conda activate '$CONDA_ENV' || { echo '[FATAL] conda activate $CONDA_ENV failed'; }" C-m \
  "cd '$PROJECT_DIR' || { echo '[FATAL] cd failed'; }" C-m \
  "export PYTHONUNBUFFERED=1 PORT='$PORT'" C-m \
  "mkdir -p logs" C-m \
  "echo \"=== \$(date '+%F %T') — launching server on :$PORT ===\" | tee -a logs/server.out" C-m \
  "if command -v gunicorn >/dev/null 2>&1; then \
      echo '[INFO] gunicorn found' | tee -a logs/server.out; \
      gunicorn 'app:app' -b 0.0.0.0:\$PORT --workers 2 --threads 4 --timeout 120 >> logs/server.out 2>&1; \
    else \
      echo '[WARN] gunicorn not found — running python app.py' | tee -a logs/server.out; \
      python app.py >> logs/server.out 2>&1; \
    fi" C-m \
  "echo '[ERROR] server process exited — keeping pane open'; bash -l" C-m

echo "Logs: $PROJECT_DIR/logs/server.out"
echo "Attach: tmux attach -t $SESSION_NAME"
echo "Detach: Ctrl-b then d"
