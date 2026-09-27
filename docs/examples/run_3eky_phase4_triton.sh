#!/usr/bin/env bash
# Copy/paste Triton demo: build a 3EKY approved-drug package, then optionally submit it.
set -euo pipefail

BASE="${BASE_URL:-https://autodockvina.com}"
JOB="${JOB:-3eky-phase4-triton-$(date +%s)}"
WORKDIR="${WORKDIR:-$PWD/$JOB}"
SUBMIT_JOBS="${SUBMIT_JOBS:-0}"

check_ok() {
  python -c 'import json, sys; payload=json.load(open(sys.argv[1])); ok=payload.get("ok", False); print(payload.get("message", "") if not ok else "", file=sys.stderr); sys.exit(0 if ok else 1)' "$1"
}

mkdir -p "$WORKDIR"
cd "$WORKDIR"
curl -fsS -X POST "$BASE/api/v1/headless/package" -H "Content-Type: application/json" -d "{
  \"workspace_name\":\"$JOB\",
  \"receptor\":{\"pdb_id\":\"3EKY\",\"chains\":\"A\"},
  \"bound_ligand\":{\"resname\":\"DR7\",\"chain\":\"A\",\"resi\":\"100\"},
  \"center\":{\"method\":\"same_as_bound_ligand\",\"size\":20},
  \"ligand\":{\"source\":\"curated\",\"library\":\"phase4\"},
  \"package\":{\"package_mode\":\"triton_lsf\",\"poses_conf\":64,\"poses_vina\":9}
}" > build.json
check_ok build.json

ARTIFACT=$(python -c 'import json; print(json.load(open("build.json"))["data"]["artifact"]["download_url"])')
curl -fL -o "${JOB}.zip" "$BASE$ARTIFACT"
mkdir package
unzip -q "${JOB}.zip" -d package
cd package/job

if [[ "$SUBMIT_JOBS" == "1" ]]; then
  CONFGEN_SUBMISSION=$(bsub < run_confgen_job.lsf)
  CONFGEN_JOB_ID=$(printf '%s\n' "$CONFGEN_SUBMISSION" | sed -n 's/.*<\([0-9][0-9]*\)>.*/\1/p')
  [[ -n "$CONFGEN_JOB_ID" ]] || { echo "Could not read ConfGen job ID: $CONFGEN_SUBMISSION" >&2; exit 1; }
  bsub -w "done($CONFGEN_JOB_ID)" < run_vina_job.lsf
  bjobs -u "$USER"
else
  printf '%s\n' "Prepared Triton package: $PWD"
  printf '%s\n' "Review hpc_profile.json, then submit with:"
  printf '%s\n' "  bsub < run_confgen_job.lsf"
  printf '%s\n' "  bsub < run_vina_job.lsf"
  printf '%s\n' "To submit both with an LSF dependency automatically, rerun with SUBMIT_JOBS=1."
fi
