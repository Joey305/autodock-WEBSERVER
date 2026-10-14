#!/usr/bin/env bash
# Copy/paste API demo: fetch 3EKY, upload the repository DR7 ligand, build a portable package.
set -euo pipefail

BASE="${BASE_URL:-https://autodockvina.com}"
JOB="${JOB:-3eky-dr7-$(date +%s)}"
WORKDIR="${WORKDIR:-$PWD/$JOB}"
DR7_URL="https://raw.githubusercontent.com/Joey305/autodock-WEBSERVER/main/Example/Ligand_SDF/DR7.sdf"

check_ok() {
  python -c 'import json, sys; payload=json.load(open(sys.argv[1])); ok=payload.get("ok", False); print(payload.get("message", "") if not ok else "", file=sys.stderr); sys.exit(0 if ok else 1)' "$1"
}

mkdir -p "$WORKDIR"
cd "$WORKDIR"
curl -fsSLo DR7.sdf "$DR7_URL"

curl -fsS -X POST "$BASE/api/v1/workspaces" -H "Content-Type: application/json" \
  -d "{\"workspace_name\":\"$JOB\",\"reuse\":true}" > workspace.json
check_ok workspace.json
curl -fsS -X POST "$BASE/api/v1/workspaces/$JOB/receptors/fetch" -H "Content-Type: application/json" \
  -d '{"pdb_id":"3EKY","chains":"A"}' > receptor.json
check_ok receptor.json
curl -fsS -X POST "$BASE/api/v1/workspaces/$JOB/centers/save" -H "Content-Type: application/json" \
  -d '{"method":"xyz","receptor":"Receptors/3eky.pdb","center":[21.116,29.387,11.742],"size":20}' > center.json
check_ok center.json
curl -fsS -X POST "$BASE/api/v1/workspaces/$JOB/prep/start" -H "Content-Type: application/json" \
  -d '{"remove_hets":"all","remove_chains":[],"altloc":"collapse"}' > prep.json
check_ok prep.json
curl -fsS -X POST "$BASE/api/v1/workspaces/$JOB/ligands/upload" \
  -F "mode=single" -F "file=@DR7.sdf" > ligands.json
check_ok ligands.json
curl -fsS -X POST "$BASE/api/v1/workspaces/$JOB/build" -H "Content-Type: application/json" \
  -d '{"package_mode":"portable","poses_conf":64,"poses_vina":9}' > build.json
check_ok build.json

ARTIFACT=$(python -c 'import json; print(json.load(open("build.json"))["data"]["download_url"])')
curl -fL -o "${JOB}_portable.zip" "$BASE$ARTIFACT"
printf 'Created %s/%s\n' "$WORKDIR" "${JOB}_portable.zip"
