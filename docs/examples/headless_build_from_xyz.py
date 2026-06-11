#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path

import requests


def check(response):
    payload = response.json()
    if not response.ok or not payload.get("ok", False):
        print(payload, file=sys.stderr)
        sys.exit(1)
    return payload["data"]


def main():
    parser = argparse.ArgumentParser(description="Create a portable package from an explicit XYZ center.")
    parser.add_argument("--base-url", default=os.environ.get("BASE_URL", "http://127.0.0.1:5050"))
    parser.add_argument("--workspace", default="xyz-build-example")
    parser.add_argument("--pdb-id", default="3EKY")
    parser.add_argument("--chains", default="A")
    parser.add_argument("--x", type=float, required=True)
    parser.add_argument("--y", type=float, required=True)
    parser.add_argument("--z", type=float, required=True)
    parser.add_argument("--ligand", required=True, help="Path to .sdf, .smiles, .smi, or .csv ligand input.")
    args = parser.parse_args()

    ligand_path = Path(args.ligand)
    s = requests.Session()
    job = check(s.post(f"{args.base_url}/api/v1/workspaces", json={"workspace_name": args.workspace, "reuse": True}))["jobname"]
    fetched = check(s.post(f"{args.base_url}/api/v1/workspaces/{job}/receptors/fetch", json={"pdb_id": args.pdb_id, "chains": args.chains}))
    check(s.post(f"{args.base_url}/api/v1/workspaces/{job}/centers/save", json={
        "method": "xyz",
        "receptor": fetched["rel"],
        "center": [args.x, args.y, args.z],
        "size": 20,
    }))
    with ligand_path.open("rb") as handle:
        check(s.post(
            f"{args.base_url}/api/v1/workspaces/{job}/ligands/upload",
            data={"mode": "single"},
            files={"file": (ligand_path.name, handle)},
        ))
    built = check(s.post(f"{args.base_url}/api/v1/workspaces/{job}/build", json={"package_mode": "portable", "poses_conf": 64, "poses_vina": 9}))
    print(f"workspace={job}")
    print(f"download_url={args.base_url}{built['download_url']}")


if __name__ == "__main__":
    main()
