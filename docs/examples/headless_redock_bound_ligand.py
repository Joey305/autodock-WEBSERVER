#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from urllib.parse import urljoin

import requests


def check(response):
    try:
        payload = response.json()
    except Exception:
        print(response.text[:2000], file=sys.stderr)
        response.raise_for_status()
        raise
    if not response.ok or not payload.get("ok", False):
        print(payload, file=sys.stderr)
        sys.exit(1)
    return payload["data"]


def split_csv(value: str):
    return [item.strip() for item in (value or "").split(",") if item.strip()]


def choose_instance(instances, args):
    candidates = [item for item in instances if not item.get("is_water")]
    if args.ligand:
        candidates = [item for item in candidates if item.get("resname", "").upper() == args.ligand.upper()]
    if args.chain:
        candidates = [item for item in candidates if item.get("chain", "").upper() == args.chain.upper()]
    if args.resi:
        candidates = [item for item in candidates if str(item.get("resi", "")) == str(args.resi)]

    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        print("No matching HETATM ligand instances were found.", file=sys.stderr)
        sys.exit(1)
    if not sys.stdin.isatty():
        print("Multiple ligand instances matched; rerun with --ligand, --chain, and --resi.", file=sys.stderr)
        for item in candidates:
            print(format_instance(item), file=sys.stderr)
        sys.exit(1)

    print("Choose a bound ligand instance:")
    for idx, item in enumerate(candidates, start=1):
        print(f"{idx}. {format_instance(item)}")
    while True:
        raw = input("> ").strip()
        try:
            choice = int(raw)
        except Exception:
            continue
        if 1 <= choice <= len(candidates):
            return candidates[choice - 1]


def format_instance(item):
    center = item.get("center") or [0, 0, 0]
    return (
        f"{item.get('resname')} chain={item.get('chain')} resi={item.get('resi')} "
        f"atoms={item.get('atom_count')} center={center[0]:.3f},{center[1]:.3f},{center[2]:.3f}"
    )


def main():
    parser = argparse.ArgumentParser(description="Create a redocking package from a PDB-bound HETATM ligand.")
    parser.add_argument("--base-url", default="https://autodockvina.com")
    parser.add_argument("--workspace", default="")
    parser.add_argument("--pdb-id", default="9G94")
    parser.add_argument("--fetch-chains", default="", help="Optional RCSB chain filter. Leave empty for mmCIF-only entries.")
    parser.add_argument("--ligand", default="A1D73", help="HETATM/component ID to redock.")
    parser.add_argument("--chain", default="")
    parser.add_argument("--resi", default="")
    parser.add_argument("--size", type=float, default=20.0)
    parser.add_argument("--remove-chains", default="", help="Comma-separated chains to remove during receptor prep.")
    parser.add_argument("--package-mode", default="portable", choices=["portable", "joey_lsf", "mainak_lsf", "custom_lsf"])
    parser.add_argument("--download-dir", default="", help="Optional local folder where the built ZIP should be downloaded.")
    parser.add_argument("--no-prep", action="store_true", help="Build without running receptor PDBQT prep.")
    args = parser.parse_args()

    workspace = args.workspace or f"{args.pdb_id.lower()}-redock"
    session = requests.Session()
    check(session.get(f"{args.base_url}/api/v1/health"))
    job = check(session.post(f"{args.base_url}/api/v1/workspaces", json={"workspace_name": workspace, "reuse": True}))["jobname"]
    fetched = check(
        session.post(
            f"{args.base_url}/api/v1/workspaces/{job}/receptors/fetch",
            json={"pdb_id": args.pdb_id, "chains": split_csv(args.fetch_chains)},
        )
    )
    hetatms = check(
        session.get(
            f"{args.base_url}/api/v1/workspaces/{job}/hetatms",
            params={"receptor": fetched["rel"]},
        )
    )
    picked = choose_instance(hetatms["instances"], args)
    built = check(
        session.post(
            f"{args.base_url}/api/v1/headless/package",
            json={
                "workspace_name": job,
                "reuse": True,
                "receptor": {"rel": fetched["rel"]},
                "bound_ligand": {
                    "resname": picked["resname"],
                    "chain": picked["chain"],
                    "resi": picked["resi"],
                    "insertion_code": picked.get("insertion_code") or "",
                },
                "center": {"method": "same_as_bound_ligand", "size": args.size},
                "prep": {
                    "run": not args.no_prep,
                    "remove_hets": "all",
                    "remove_chains": split_csv(args.remove_chains),
                    "altloc": "collapse",
                },
                "package": {"package_mode": args.package_mode, "poses_conf": 64, "poses_vina": 9},
            },
        )
    )

    artifact = built["artifact"]
    download_url = urljoin(args.base_url, artifact["download_url"])
    print(f"workspace={built['jobname']}")
    print(f"ligand={format_instance(picked)}")
    print(f"center={built['center']['center']}")
    print(f"download_url={download_url}")

    if args.download_dir:
        out_dir = Path(args.download_dir).expanduser()
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / Path(artifact["zip"]).name
        with session.get(download_url, stream=True) as response:
            response.raise_for_status()
            with out_path.open("wb") as handle:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        handle.write(chunk)
        print(f"downloaded={out_path}")


if __name__ == "__main__":
    main()
