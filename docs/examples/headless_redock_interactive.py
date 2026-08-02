#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
import zipfile
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


def prompt_text(label: str, default: str = "") -> str:
    suffix = f" [{default}]" if default else ""
    value = input(f"{label}{suffix}: ").strip()
    return value or default


def prompt_float(label: str, default: float) -> float:
    while True:
        raw = prompt_text(label, str(default))
        try:
            return float(raw)
        except Exception:
            print("Please enter a number.")


def prompt_yes_no(label: str, default: bool = True) -> bool:
    hint = "Y/n" if default else "y/N"
    while True:
        raw = input(f"{label} [{hint}]: ").strip().lower()
        if not raw:
            return default
        if raw in {"y", "yes"}:
            return True
        if raw in {"n", "no"}:
            return False


def prompt_choice(label: str, choices: list[str], default: str) -> str:
    normalized = {choice.lower(): choice for choice in choices}
    while True:
        raw = prompt_text(f"{label} ({'/'.join(choices)})", default).lower()
        if raw in normalized:
            return normalized[raw]
        print(f"Choose one of: {', '.join(choices)}")


def format_instance(item):
    center = item.get("center") or [0, 0, 0]
    return (
        f"{item.get('resname')} chain={item.get('chain') or '-'} resi={item.get('resi') or '-'} "
        f"atoms={item.get('atom_count')} center={center[0]:.3f},{center[1]:.3f},{center[2]:.3f}"
    )


def filter_candidates(instances, ligand="", chain="", resi=""):
    candidates = [item for item in instances if not item.get("is_water")]
    if ligand:
        candidates = [item for item in candidates if item.get("resname", "").upper() == ligand.upper()]
    if chain:
        candidates = [item for item in candidates if item.get("chain", "").upper() == chain.upper()]
    if resi:
        candidates = [item for item in candidates if str(item.get("resi", "")) == str(resi)]
    return candidates


def choose_instance(instances, ligand="", chain="", resi=""):
    candidates = filter_candidates(instances, ligand, chain, resi)
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        print("No matching non-water HETATM instances were found.")
        candidates = filter_candidates(instances)
    if not candidates:
        print("No non-water HETATM instances were available for selection.", file=sys.stderr)
        sys.exit(1)

    print("\nAvailable non-water HETATM instances:")
    for idx, item in enumerate(candidates, start=1):
        print(f"{idx:>2}. {format_instance(item)}")

    while True:
        raw = input("Select ligand number: ").strip()
        try:
            choice = int(raw)
        except Exception:
            print("Enter a number from the list.")
            continue
        if 1 <= choice <= len(candidates):
            return candidates[choice - 1]
        print("That number is not in the list.")


def suggest_remove_chains(instances, picked) -> str:
    resname = picked.get("resname")
    picked_chain = picked.get("chain")
    chains = sorted(
        {
            item.get("chain")
            for item in instances
            if item.get("resname") == resname and item.get("chain") and item.get("chain") != picked_chain
        }
    )
    return ",".join(chains)


def parse_xyz(raw: str):
    parts = [part.strip() for part in raw.replace(",", " ").split() if part.strip()]
    if len(parts) != 3:
        raise ValueError("XYZ center needs exactly three values.")
    return [float(parts[0]), float(parts[1]), float(parts[2])]


def download_artifact(session, base_url: str, artifact: dict, download_dir: str, job: str, unpack: bool):
    out_dir = Path(download_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    download_url = urljoin(base_url, artifact["download_url"])
    out_path = out_dir / f"{job}.zip"
    with session.get(download_url, stream=True) as response:
        response.raise_for_status()
        with out_path.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
    print(f"downloaded={out_path}")
    if unpack:
        extract_dir = out_dir / job
        extract_dir.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(out_path) as zf:
            zf.extractall(extract_dir)
        print(f"unpacked={extract_dir}")


def main():
    parser = argparse.ArgumentParser(description="Interactive redocking package builder for AutoDock-Vina PrepServer.")
    parser.add_argument("--base-url", default="")
    parser.add_argument("--workspace", default="")
    parser.add_argument("--pdb-id", default="")
    parser.add_argument("--fetch-chains", default="")
    parser.add_argument("--ligand", default="")
    parser.add_argument("--chain", default="")
    parser.add_argument("--resi", default="")
    parser.add_argument("--size", type=float)
    parser.add_argument("--remove-chains", default="")
    parser.add_argument("--package-mode", default="")
    parser.add_argument("--download-dir", default="")
    parser.add_argument("--no-prep", action="store_true")
    parser.add_argument("--no-unpack", action="store_true")
    parser.add_argument("--yes", action="store_true", help="Accept defaults for prompts not supplied as flags.")
    args = parser.parse_args()

    interactive = sys.stdin.isatty() and not args.yes
    base_url = args.base_url or "https://autodockvina.com"
    pdb_id = args.pdb_id or "9G94"
    fetch_chains = args.fetch_chains
    workspace = args.workspace or f"{pdb_id.lower()}-redock"

    if interactive:
        print("AutoDock-Vina headless redocking setup\n")
        base_url = prompt_text("Base URL", base_url)
        pdb_id = prompt_text("PDB ID", pdb_id).upper()
        workspace = prompt_text("Workspace name", workspace or f"{pdb_id.lower()}-redock")
        fetch_chains = prompt_text("Fetch only chains (comma-separated, blank for all)", fetch_chains)

    session = requests.Session()
    check(session.get(f"{base_url}/api/v1/health"))
    job = check(session.post(f"{base_url}/api/v1/workspaces", json={"workspace_name": workspace, "reuse": True}))["jobname"]
    print(f"workspace={job}")

    fetched = check(
        session.post(
            f"{base_url}/api/v1/workspaces/{job}/receptors/fetch",
            json={"pdb_id": pdb_id, "chains": split_csv(fetch_chains)},
        )
    )
    print(f"receptor={fetched['rel']}")

    hetatms = check(
        session.get(
            f"{base_url}/api/v1/workspaces/{job}/hetatms",
            params={"receptor": fetched["rel"]},
        )
    )
    picked = choose_instance(hetatms["instances"], args.ligand, args.chain, args.resi)
    print(f"selected={format_instance(picked)}")

    size = args.size if args.size is not None else 20.0
    center_payload = {"method": "same_as_bound_ligand", "size": size}
    if interactive:
        center_mode = prompt_choice("Docking center", ["ligand", "xyz"], "ligand")
        if center_mode == "xyz":
            while True:
                try:
                    center_payload = {
                        "method": "xyz",
                        "center": parse_xyz(prompt_text("XYZ center as x,y,z")),
                        "size": prompt_float("Box size", size),
                    }
                    break
                except Exception as exc:
                    print(exc)
        else:
            center_payload["size"] = prompt_float("Box size", size)

    suggested_remove = suggest_remove_chains(hetatms["instances"], picked)
    remove_chains = args.remove_chains
    if interactive:
        remove_chains = prompt_text("Remove receptor chains during prep", remove_chains or suggested_remove)
    elif not remove_chains:
        remove_chains = suggested_remove

    package_mode = args.package_mode or "portable"
    if interactive:
        package_mode = prompt_choice("Package mode", ["portable", "joey_lsf", "mainak_lsf", "custom_lsf"], package_mode)

    run_prep = not args.no_prep
    if interactive:
        run_prep = prompt_yes_no("Run receptor prep now", run_prep)

    download_dir = args.download_dir or str(Path.home() / "Docking")
    unpack = not args.no_unpack
    if interactive:
        download_dir = prompt_text("Download folder", download_dir)
        unpack = prompt_yes_no("Unpack downloaded ZIP", unpack)

    print("\nBuilding package...")
    built = check(
        session.post(
            f"{base_url}/api/v1/headless/package",
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
                "center": center_payload,
                "prep": {
                    "run": run_prep,
                    "remove_hets": "all",
                    "remove_chains": split_csv(remove_chains),
                    "altloc": "collapse",
                },
                "package": {"package_mode": package_mode, "poses_conf": 64, "poses_vina": 9},
            },
        )
    )

    artifact = built["artifact"]
    download_url = urljoin(base_url, artifact["download_url"])
    print(f"center={built['center']['center']}")
    print(f"download_url={download_url}")
    download_artifact(session, base_url, artifact, download_dir, job, unpack)


if __name__ == "__main__":
    main()
