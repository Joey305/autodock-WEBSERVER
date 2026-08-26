"""PyMOL command helpers for companion ``.results.json`` docking manifests.

Run this file in PyMOL, then use ``docking_results <manifest.json>``.  The
command prints the filtered score table and optionally shows a selected pose;
it deliberately works without Qt so it remains portable across PyMOL builds.
"""
import json
from pathlib import Path


def load_results(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def rows_for(data, receptor, ligand):
    return data.get("receptors", {}).get(receptor, {}).get("ligands", {}).get(ligand, [])


def best_score(rows):
    return min((float(row["score"]) for row in rows), default=None)


def docking_results(path, receptor="", ligand="", rank=""):
    data = load_results(path)
    receptors = data.get("receptors", {})
    receptor = receptor or next(iter(receptors), "")
    ligands = receptors.get(receptor, {}).get("ligands", {})
    ligand = ligand or next(iter(ligands), "")
    rows = rows_for(data, receptor, ligand)
    print(f"Docking Results: {receptor} / {ligand}; Best: {best_score(rows):.2f} kcal/mol" if rows else "No docking rows")
    for row in rows:
        print("{rank:>3}  {score:>7.2f}  {protomer:<4} {tautomer:<4} {conformer:<4} {vina_pose}".format(**row))
    if rank:
        selected = next((row for row in rows if int(row.get("rank", 0)) == int(rank)), None)
        if selected:
            from pymol import cmd
            cmd.enable(selected["object"])
            cmd.zoom(selected["object"])


try:
    from pymol import cmd
    cmd.extend("docking_results", docking_results)
except Exception:
    pass
