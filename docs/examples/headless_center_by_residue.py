#!/usr/bin/env python3
import argparse
import os
import sys

import requests


def check(response):
    try:
        payload = response.json()
    except Exception:
        response.raise_for_status()
        raise
    if not response.ok or not payload.get("ok", False):
        print(payload, file=sys.stderr)
        sys.exit(1)
    return payload["data"]


def main():
    parser = argparse.ArgumentParser(description="Resolve and save a docking center from a residue selector.")
    parser.add_argument("--base-url", default=os.environ.get("BASE_URL", "http://127.0.0.1:5050"))
    parser.add_argument("--workspace", default="residue-center-example")
    parser.add_argument("--pdb-id", default="3EKY")
    parser.add_argument("--chains", default="A")
    parser.add_argument("--receptor", default="")
    parser.add_argument("--chain", default="A")
    parser.add_argument("--resi", default="82")
    parser.add_argument("--resname", default="GLU")
    parser.add_argument("--size", type=float, default=20)
    args = parser.parse_args()

    s = requests.Session()
    job = check(s.post(f"{args.base_url}/api/v1/workspaces", json={"workspace_name": args.workspace, "reuse": True}))["jobname"]
    fetched = check(s.post(f"{args.base_url}/api/v1/workspaces/{job}/receptors/fetch", json={"pdb_id": args.pdb_id, "chains": args.chains}))
    receptor = args.receptor or fetched["rel"]
    center = check(s.post(f"{args.base_url}/api/v1/workspaces/{job}/centers/save", json={
        "method": "residue",
        "receptor": receptor,
        "chain": args.chain,
        "resi": args.resi,
        "resname": args.resname,
        "size": args.size,
    }))
    print(f"workspace={job}")
    print(f"center={center['center']} atom_count={center['matched'].get('atom_count')}")


if __name__ == "__main__":
    main()
