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
    parser = argparse.ArgumentParser(description="Run the practical staged headless package workflow.")
    parser.add_argument("--base-url", default=os.environ.get("BASE_URL", "http://127.0.0.1:5050"))
    parser.add_argument("--workspace", default="headless-package-example")
    parser.add_argument("--pdb-id", default="3EKY")
    parser.add_argument("--chains", default="A")
    parser.add_argument("--center-method", choices=["residue", "hetatm", "xyz"], default="residue")
    parser.add_argument("--chain", default="A")
    parser.add_argument("--resi", default="82")
    parser.add_argument("--resname", default="GLU")
    parser.add_argument("--het", default="DR7")
    parser.add_argument("--xyz", nargs=3, type=float)
    parser.add_argument("--ligand", required=True)
    parser.add_argument("--package-mode", choices=["portable", "lsf"], default="portable")
    args = parser.parse_args()

    s = requests.Session()
    check(s.get(f"{args.base_url}/api/v1/health"))
    job = check(s.post(f"{args.base_url}/api/v1/workspaces", json={"workspace_name": args.workspace, "reuse": True}))["jobname"]
    fetched = check(s.post(f"{args.base_url}/api/v1/workspaces/{job}/receptors/fetch", json={"pdb_id": args.pdb_id, "chains": args.chains}))

    center_payload = {"method": args.center_method, "receptor": fetched["rel"], "size": 20}
    if args.center_method == "residue":
        center_payload.update({"chain": args.chain, "resi": args.resi, "resname": args.resname})
    elif args.center_method == "hetatm":
        center_payload.update({"het": args.het, "chain": args.chain, "resi": args.resi})
    else:
        if not args.xyz:
            parser.error("--xyz X Y Z is required for --center-method xyz")
        center_payload["center"] = args.xyz
    center = check(s.post(f"{args.base_url}/api/v1/workspaces/{job}/centers/save", json=center_payload))

    ligand_path = Path(args.ligand)
    with ligand_path.open("rb") as handle:
        check(s.post(
            f"{args.base_url}/api/v1/workspaces/{job}/ligands/upload",
            data={"mode": "single"},
            files={"file": (ligand_path.name, handle)},
        ))

    built = check(s.post(f"{args.base_url}/api/v1/workspaces/{job}/build", json={
        "package_mode": args.package_mode,
        "poses_conf": 64,
        "poses_vina": 9,
    }))
    print(f"workspace={job}")
    print(f"center={center['center']}")
    print(f"download_url={args.base_url}{built['download_url']}")


if __name__ == "__main__":
    main()
