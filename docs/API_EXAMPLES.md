# API Examples

Set a base URL:

```bash
BASE="${BASE_URL:-http://127.0.0.1:5050}"
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
