# API Examples

Set a base URL:

```bash
BASE="${BASE_URL:-https://autodockvina.com}"
```

For local development, use:

```bash
BASE="http://127.0.0.1:5050"
```

## Health

```bash
curl -s "$BASE/api/v1/health" | python -m json.tool
```

## Create A Workspace

```bash
curl -s -X POST "$BASE/api/v1/workspaces" \
  -H "Content-Type: application/json" \
  -d '{"workspace_name":"api-example","reuse":true}' | python -m json.tool
```

## Fetch A Receptor

```bash
curl -s -X POST "$BASE/api/v1/workspaces/api-example/receptors/fetch" \
  -H "Content-Type: application/json" \
  -d '{"pdb_id":"3EKY","chains":["A"]}' | python -m json.tool
```

## Upload A Receptor

```bash
curl -s -X POST "$BASE/api/v1/workspaces/api-example/receptors/upload" \
  -F mode=single \
  -F file=@example.pdb | python -m json.tool
```

## Resolve And Save Residue Center

```bash
curl -s -X POST "$BASE/api/v1/workspaces/api-example/centers/save" \
  -H "Content-Type: application/json" \
  -d '{"method":"residue","receptor":"3eky.pdb","chain":"A","resi":"82","resname":"GLU","size":20}' \
  | python -m json.tool
```

## Resolve And Save HETATM Center

```bash
curl -s -X POST "$BASE/api/v1/workspaces/api-example/centers/save" \
  -H "Content-Type: application/json" \
  -d '{"method":"hetatm","receptor":"3eky.pdb","het":"DR7","chain":"A","resi":"100","size":20}' \
  | python -m json.tool
```

## List Bound HETATM Ligands

```bash
curl -s "$BASE/api/v1/workspaces/api-example/hetatms?receptor=3eky.pdb" \
  | python -m json.tool
```

## Save Explicit XYZ Center

```bash
curl -s -X POST "$BASE/api/v1/workspaces/api-example/centers/save" \
  -H "Content-Type: application/json" \
  -d '{"method":"xyz","receptor":"3eky.pdb","center":[11.0,21.0,31.0],"size":20}' \
  | python -m json.tool
```

## Upload Ligands

```bash
curl -s -X POST "$BASE/api/v1/workspaces/api-example/ligands/upload" \
  -F mode=single \
  -F file=@ligands.sdf | python -m json.tool
```

## Extract A Bound HETATM Ligand

```bash
curl -s -X POST "$BASE/api/v1/workspaces/api-example/ligands/extract" \
  -H "Content-Type: application/json" \
  -d '{"receptor":"3eky.pdb","resname":"DR7","chain":"A","resi":"100"}' \
  | python -m json.tool
```

## Prepare Receptors

```bash
curl -s -X POST "$BASE/api/v1/workspaces/api-example/prep/start" \
  -H "Content-Type: application/json" \
  -d '{"remove_hets":"all","remove_chains":[],"altloc":"collapse"}' | python -m json.tool
```

## Poll Status

```bash
curl -s "$BASE/api/v1/workspaces/api-example/summary" | python -m json.tool
curl -s "$BASE/api/v1/workspaces/api-example/prep/status" | python -m json.tool
```

## Build Package

```bash
curl -s -X POST "$BASE/api/v1/workspaces/api-example/build" \
  -H "Content-Type: application/json" \
  -d '{"package_mode":"portable","poses_conf":64,"poses_vina":9}' | python -m json.tool
```

## One-Call Bound-Ligand Redocking

```bash
curl -s -X POST "$BASE/api/v1/headless/package" \
  -H "Content-Type: application/json" \
  -d '{
    "workspace_name": "9g94-redock",
    "reuse": true,
    "receptor": {"pdb_id": "9G94"},
    "bound_ligand": {"resname": "A1D73", "chain": "A", "resi": "101"},
    "center": {"method": "same_as_bound_ligand", "size": 20},
    "prep": {"remove_hets": "all", "remove_chains": ["B", "C"]},
    "package": {"package_mode": "portable", "poses_conf": 64, "poses_vina": 9}
  }' | python -m json.tool
```

Guided interactive wrapper:

```bash
curl -fsSLo autodock-redock-interactive.py \
  "$BASE/api/v1/clients/headless_redock_interactive.py"

python autodock-redock-interactive.py
```

The guided wrapper fetches the receptor, shows non-water HETATM candidates with centroids, and lets you choose the ligand/center/package/download options.

Flag-driven wrapper:

```bash
curl -fsSLo autodock-redock-bound-ligand.py \
  "$BASE/api/v1/clients/headless_redock_bound_ligand.py"

python autodock-redock-bound-ligand.py \
  --base-url "$BASE" \
  --pdb-id 9G94 \
  --ligand A1D73 \
  --download-dir "$HOME/Docking"
```

Or run without saving the wrapper when using bash or zsh:

```bash
python <(curl -fsSL "$BASE/api/v1/clients/headless_redock_interactive.py")
```

Non-interactive server run:

```bash
python <(curl -fsSL "$BASE/api/v1/clients/headless_redock_interactive.py") \
  --base-url "$BASE" \
  --pdb-id 9G94 \
  --ligand A1D73 \
  --chain A \
  --resi 101 \
  --download-dir "$HOME/Docking" \
  --yes
```

## Python Skeleton

```python
import os
import requests

base = os.environ.get("BASE_URL", "http://127.0.0.1:5050")
session = requests.Session()

workspace = session.post(f"{base}/api/v1/workspaces", json={"workspace_name": "python-api", "reuse": True}).json()["data"]
job = workspace["jobname"]
session.post(f"{base}/api/v1/workspaces/{job}/receptors/fetch", json={"pdb_id": "3EKY", "chains": ["A"]}).raise_for_status()
center = session.post(f"{base}/api/v1/workspaces/{job}/centers/save", json={
    "method": "residue",
    "receptor": "3eky.pdb",
    "chain": "A",
    "resi": "82",
    "resname": "GLU",
    "size": 20,
}).json()
print(center)
```
