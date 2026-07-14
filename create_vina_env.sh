#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${ENV_NAME:-vina_env}"
YAML_FILE="${YAML_FILE:-docking.yaml}"
RECREATE="${RECREATE:-0}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --name)
      ENV_NAME="$2"
      shift 2
      ;;
    --yaml)
      YAML_FILE="$2"
      shift 2
      ;;
    --recreate)
      RECREATE=1
      shift
      ;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: ./create_vina_env.sh [--name vina_env] [--yaml docking.yaml] [--recreate]" >&2
      exit 2
      ;;
  esac
done

if ! command -v conda >/dev/null 2>&1; then
  echo "❌ conda is not available on PATH." >&2
  exit 2
fi

if [[ ! -f "$YAML_FILE" ]]; then
  echo "❌ Environment file not found: $YAML_FILE" >&2
  exit 2
fi

CONDA_BASE="$(conda info --base)"
source "$CONDA_BASE/etc/profile.d/conda.sh"

if conda env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
  if [[ "$RECREATE" == "1" ]]; then
    echo "♻️ Removing existing environment: $ENV_NAME"
    conda env remove -n "$ENV_NAME" -y
  else
    echo "ℹ️ Environment already exists: $ENV_NAME"
    echo "   Re-run with --recreate to rebuild it."
    exit 0
  fi
fi

TMP_FULL="$(mktemp "${TMPDIR:-/tmp}/vina_env_full.XXXXXX.yaml")"
TMP_NOVINA="$(mktemp "${TMPDIR:-/tmp}/vina_env_novina.XXXXXX.yaml")"
cleanup() {
  rm -f "$TMP_FULL" "$TMP_NOVINA"
}
trap cleanup EXIT

python3 - "$YAML_FILE" "$ENV_NAME" "$TMP_FULL" "$TMP_NOVINA" <<'PY'
import sys
from pathlib import Path

src = Path(sys.argv[1])
env_name = sys.argv[2]
full_dst = Path(sys.argv[3])
novina_dst = Path(sys.argv[4])

lines = src.read_text().splitlines()

def rewrite_name(input_lines):
    out = []
    replaced = False
    for line in input_lines:
        if line.startswith("name:"):
            out.append(f"name: {env_name}")
            replaced = True
        else:
            out.append(line)
    if not replaced:
        out.insert(0, f"name: {env_name}")
    return out

full_lines = rewrite_name(lines)
novina_lines = [line for line in full_lines if line.strip() != "- vina"]

full_dst.write_text("\n".join(full_lines) + "\n")
novina_dst.write_text("\n".join(novina_lines) + "\n")
PY

echo "🚀 Creating conda environment '$ENV_NAME' from $YAML_FILE"
if conda env create -f "$TMP_FULL"; then
  echo "✅ Environment created successfully."
else
  echo "⚠️ Initial environment creation failed. Retrying without vina in YAML, then installing vina separately."
  conda env remove -n "$ENV_NAME" -y >/dev/null 2>&1 || true
  conda env create -f "$TMP_NOVINA"
  conda install -n "$ENV_NAME" -c conda-forge vina -y
fi

echo "🔎 Verifying environment..."
conda run -n "$ENV_NAME" python --version
conda run -n "$ENV_NAME" obabel -V
conda run -n "$ENV_NAME" vina --version

echo
echo "✅ Ready."
echo "Activate with:"
echo "  conda activate $ENV_NAME"
