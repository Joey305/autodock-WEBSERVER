#!/usr/bin/env python3
import argparse
import os
import sys

import requests


def check(response):
    payload = response.json()
    if not response.ok or not payload.get("ok", False):
        print(payload, file=sys.stderr)
        sys.exit(1)
    return payload["data"]


def main():
    parser = argparse.ArgumentParser(description="Resolve and save a docking center from a HETATM ligand instance.")
    parser.add_argument("--base-url", default=os.environ.get("BASE_URL", "http://127.0.0.1:5050"))
    parser.add_argument("--workspace", default="hetatm-center-example")
    parser.add_argument("--pdb-id", default="3EKY")
    parser.add_argument("--chains", default="A")
    parser.add_argument("--receptor", default="")
    parser.add_argument("--het", default="DR7")
    parser.add_argument("--chain", default="A")
    parser.add_argument("--resi", default="100")
    parser.add_argument("--size", type=float, default=20)
    args = parser.parse_args()

    s = requests.Session()
    job = check(s.post(f"{args.base_url}/api/v1/workspaces", json={"workspace_name": args.workspace, "reuse": True}))["jobname"]
    fetched = check(s.post(f"{args.base_url}/api/v1/workspaces/{job}/receptors/fetch", json={"pdb_id": args.pdb_id, "chains": args.chains}))
    receptor = args.receptor or fetched["rel"]
    center = check(s.post(f"{args.base_url}/api/v1/workspaces/{job}/centers/save", json={
        "method": "hetatm",
        "receptor": receptor,
        "het": args.het,
        "chain": args.chain,
        "resi": args.resi,
        "size": args.size,
    }))
    print(f"workspace={job}")
    print(f"center={center['center']} atom_count={center['matched'].get('atom_count')}")


if __name__ == "__main__":
    main()
