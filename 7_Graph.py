#!/usr/bin/env python3
"""
plot_vina_spread.py

Create reusable histogram + box/whisker plots from Vina docking CSV outputs.

Key features:
- Normalizes ligand IDs by stripping pose suffixes like _pose31 or __pose31
- Groups by (Receptor, LigandBase) and uses ALL affinities across conformers/poses
- Produces:
  - Per-ligand hist + boxplot (within each receptor)
  - Per-receptor combined boxplot for Top-N ligands
  - Per-receptor pooled histogram for Top-N ligands
- Writes outputs into an organized folder structure.

Input CSV expected columns (case-insensitive):
  Receptor, Ligand, Pose, Binding_Affinity, OutFile, Run_Timestamp
"""

import re
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
import time
import argparse
from pathlib import Path


# -------------------------
# Ligand normalization
# -------------------------
POSE_SUFFIX_RE = re.compile(r"(?i)(?:__pose\d+|_pose\d+)$")

def normalize_ligand_id(ligand: str, outfile: str | None = None) -> str:
    """
    Convert ligand identifiers like:
      UMF504_pose31, UMF504_pose8, SOME__pose12  -> UMF504, SOME

    Falls back to OutFile-derived ligand folder name if ligand is missing.
    """
    lig = (ligand or "").strip()

    if not lig and outfile:
        try:
            # OutFile example: 9UWI/UMF504_pose31/out.pdbqt  -> ligand folder: UMF504_pose31
            lig = Path(str(outfile)).parts[-2]
        except Exception:
            lig = ""

    if not lig:
        return "UNKNOWN"

    lig = POSE_SUFFIX_RE.sub("", lig)
    return lig


# -------------------------
# Plot helpers
# -------------------------
def _annotate_best_mean_worst(ax, values: np.ndarray):
    """
    Annotate best (min), mean, worst (max) with vertical lines + text.
    Vina: more negative = better, so "best" is min(values).
    """
    if values.size == 0:
        return

    best = float(np.min(values))
    worst = float(np.max(values))
    mean = float(np.mean(values))

    # vertical reference lines
    ax.axvline(best, linestyle="--", linewidth=1.5, label=f"Best: {best:.3f}")
    ax.axvline(mean, linestyle=":", linewidth=2.0,  label=f"Mean: {mean:.3f}")
    ax.axvline(worst, linestyle="--", linewidth=1.5, label=f"Worst: {worst:.3f}")

    # place text near top
    ymax = ax.get_ylim()[1]
    ax.text(best, ymax * 0.95, f"best\n{best:.3f}", ha="center", va="top", fontsize=9)
    ax.text(mean, ymax * 0.90, f"mean\n{mean:.3f}", ha="center", va="top", fontsize=9)
    ax.text(worst, ymax * 0.95, f"worst\n{worst:.3f}", ha="center", va="top", fontsize=9)


def plot_histogram(values: np.ndarray, title: str, outpath: Path, bins: int = 30):
    plt.figure(figsize=(9, 6))
    plt.hist(values, bins=bins, edgecolor="black", linewidth=0.5, alpha=0.85)
    plt.xlabel("Docking affinity (kcal/mol)")
    plt.ylabel("Count")
    plt.title(title)

    _annotate_best_mean_worst(plt.gca(), values)

    plt.legend(frameon=True, fontsize=9)
    plt.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(outpath, dpi=300)
    plt.close()

def plot_overlay_histograms(values_by_label: dict[str, np.ndarray],
                            title: str,
                            outpath: Path,
                            bins: int = 30,
                            density: bool = True,
                            alpha: float = 0.35,
                            max_legend: int = 25):
    """
    Overlay histograms for many ligands on the same axes.

    - values_by_label: {ligand: array_of_affinities}
    - density=True makes comparisons fair when ligands have different N
    - Uses shared bins computed from global min/max across ligands
    """
    labels = list(values_by_label.keys())
    all_vals = np.concatenate([np.asarray(values_by_label[k], dtype=float) for k in labels])
    if all_vals.size == 0:
        return

    # shared bins across all ligands (critical for apples-to-apples)
    vmin = float(np.min(all_vals))
    vmax = float(np.max(all_vals))
    bin_edges = np.linspace(vmin, vmax, bins + 1)

    plt.figure(figsize=(10, 6.5))
    ax = plt.gca()

    # one histogram per ligand
    for lab in labels:
        vals = np.asarray(values_by_label[lab], dtype=float)
        if vals.size == 0:
            continue
        ax.hist(
            vals,
            bins=bin_edges,
            density=density,
            alpha=alpha,
            edgecolor="black",
            linewidth=0.3,
            label=lab
        )

    ax.set_xlabel("Docking affinity (kcal/mol)")
    ax.set_ylabel("Density" if density else "Count")
    ax.set_title(title)

    # annotate global best/mean/worst for the receptor's top ligands set
    _annotate_best_mean_worst(ax, all_vals)

    # legend handling
    if len(labels) <= max_legend:
        ax.legend(frameon=True, fontsize=9, ncol=2)
    else:
        ax.legend().remove()
        ax.text(
            0.99, 0.98,
            f"Legend hidden (N={len(labels)}). Use summary CSV for labels.",
            transform=ax.transAxes,
            ha="right", va="top",
            fontsize=9
        )

    plt.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(outpath, dpi=300)
    plt.close()



def plot_single_box(values: np.ndarray, title: str, outpath: Path):
    """
    Single-series boxplot with mean + median annotated.
    """
    plt.figure(figsize=(7.5, 6))
    ax = plt.gca()
    ax.boxplot([values], labels=[""], showfliers=False)

    median = float(np.median(values))
    mean = float(np.mean(values))
    q3 = float(np.quantile(values, 0.75))

    # annotate on plot
    y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
    y_offset = 0.02 * y_range

    ax.text(1, median + y_offset, f"median={median:.3f}", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.text(1.08, q3, f"μ={mean:.3f}", ha="left", va="center", fontsize=10)

    ax.set_ylabel("Docking affinity (kcal/mol)")
    ax.set_title(title)
    ax.grid(axis="y", linestyle="--", alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(outpath, dpi=300)
    plt.close()


def plot_multi_box(data_list, labels, title: str, outpath: Path):
    """
    Multi-ligand boxplot for a receptor (Top-N ligands), with mean annotated per ligand.
    """
    plt.figure(figsize=(max(10, 0.45 * len(labels)), 6.5))
    ax = plt.gca()

    ax.boxplot(data_list, labels=labels, showfliers=False)

    # annotate mean above each box (at Q3 height)
    y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
    y_offset = 0.02 * y_range
    for i, vals in enumerate(data_list, start=1):
        vals = np.asarray(vals, dtype=float)
        mean = float(np.mean(vals))
        q3 = float(np.quantile(vals, 0.75))
        ax.text(i, q3 + y_offset, f"μ={mean:.2f}", ha="center", va="bottom", fontsize=8)

    ax.set_ylabel("Docking affinity (kcal/mol)")
    ax.set_title(title)
    ax.grid(axis="y", linestyle="--", alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.xticks(rotation=60, ha="right", fontsize=9)
    plt.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(outpath, dpi=300)
    plt.close()


# -------------------------
# Core logic
# -------------------------
def build_summary(df: pd.DataFrame) -> pd.DataFrame:
    """
    Returns a per (Receptor, LigandBase) summary:
      n_scores, best(min), worst(max), mean, median, std
    """
    g = df.groupby(["Receptor", "LigandBase"], as_index=False)["Binding_Affinity"]
    out = g.agg(
        n_scores="count",
        best_affinity="min",
        worst_affinity="max",
        mean_affinity="mean",
        median_affinity="median",
        std_affinity="std",
    )
    # rank helper: best first (more negative)
    out = out.sort_values(["Receptor", "best_affinity", "median_affinity"], ascending=[True, True, True])
    return out
def find_csv_candidates(root: Path, pattern: str) -> list[Path]:
    root = root.expanduser().resolve()
    # Search common places: root and one level down, then recursive
    cands = list(root.glob(pattern))
    cands += list(root.glob(f"**/{pattern}"))
    # de-dupe
    uniq = {}
    for p in cands:
        try:
            uniq[str(p.resolve())] = p.resolve()
        except Exception:
            pass
    return sorted(uniq.values(), key=lambda p: p.stat().st_mtime, reverse=True)

def choose_path(prompt: str, items: list[Path]) -> Path:
    print(f"\n{prompt}")
    for i, p in enumerate(items):
        mtime = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(p.stat().st_mtime))
        print(f" [{i}] {p}   (mtime: {mtime})")
    while True:
        s = input("Enter index: ").strip()
        try:
            i = int(s)
            if 0 <= i < len(items):
                return items[i]
        except ValueError:
            pass
        print("Invalid choice. Try again.")




def main():
    ap = argparse.ArgumentParser(description="Plot histograms and boxplots for Vina docking result CSVs.")
    ap.add_argument("--csv", default=None,
                help="Input docking CSV. If omitted, you will be prompted to select one. "
                     "Use --csv auto to pick the newest match automatically.")
    ap.add_argument("--root", default=".",
                    help="Where to search for CSVs if --csv not provided (default: current dir).")
    ap.add_argument("--pattern", default="*vina_docking_scores_sorted.csv",
                    help="Glob pattern used to find candidate CSVs (default: *vina_docking_scores_sorted.csv).")
    ap.add_argument("--non-interactive", action="store_true",
                    help="If --csv is missing, do not prompt; fail instead (useful for automation).")
    ap.add_argument("--outdir", default=None, help="Output directory (default: alongside CSV, in Plots/)")
    ap.add_argument("--top-n", type=int, default=25, help="Top N ligands per receptor for combined plots (default 25)")
    ap.add_argument("--bins", type=int, default=30, help="Histogram bins (default 30)")
    ap.add_argument("--min-scores", type=int, default=1, help="Skip ligands with < min scores (default 1)")
    ap.add_argument("--skip-per-ligand", action="store_true", help="Do not write per-ligand plots")
    ap.add_argument("--skip-per-receptor", action="store_true", help="Do not write per-receptor combined plots")
    args = ap.parse_args()

    root = Path(args.root)

    if args.csv is None or str(args.csv).strip().lower() == "auto":
        cands = find_csv_candidates(root, args.pattern)
        if not cands:
            raise RuntimeError(f"No CSVs found under {root} matching pattern: {args.pattern}")

        if str(args.csv).strip().lower() == "auto":
            in_csv = cands[0]
            print(f"✅ Auto-selected newest CSV: {in_csv}")
        else:
            if args.non_interactive:
                raise RuntimeError("Missing --csv and --non-interactive set. Provide --csv explicitly.")
            in_csv = choose_path("Select a docking CSV to plot:", cands)
    else:
        in_csv = Path(args.csv).expanduser().resolve()
        if not in_csv.exists():
            raise FileNotFoundError(f"CSV not found: {in_csv}")


    outdir = Path(args.outdir).expanduser().resolve() if args.outdir else in_csv.parent / "Plots"
    outdir.mkdir(parents=True, exist_ok=True)

    # Load
    df = pd.read_csv(in_csv, low_memory=False)

    # Column normalization (case-insensitive)
    colmap = {c.lower(): c for c in df.columns}
    required = ["receptor", "ligand", "binding_affinity"]
    for r in required:
        if r not in colmap:
            raise ValueError(f"Missing required column: {r} (found columns: {list(df.columns)})")

    receptor_col = colmap["receptor"]
    ligand_col = colmap["ligand"]
    aff_col = colmap["binding_affinity"]
    outfile_col = colmap.get("outfile", None)

    # Clean / numeric
    df = df.rename(columns={receptor_col: "Receptor", ligand_col: "Ligand", aff_col: "Binding_Affinity"})
    df["Binding_Affinity"] = pd.to_numeric(df["Binding_Affinity"], errors="coerce")
    df = df.dropna(subset=["Receptor", "Binding_Affinity"])

    # LigandBase
    if outfile_col:
        df["OutFile"] = df[outfile_col].astype(str)
        df["LigandBase"] = [
            normalize_ligand_id(lig, out)
            for lig, out in zip(df["Ligand"].astype(str), df["OutFile"])
        ]
    else:
        df["LigandBase"] = [normalize_ligand_id(lig, None) for lig in df["Ligand"].astype(str)]

    # Summary table
    summary = build_summary(df)

    summary_path = outdir / f"{in_csv.stem}_ligand_summary.csv"
    summary.to_csv(summary_path, index=False)
    print(f"✅ Wrote summary: {summary_path}")

    # -------------------------
    # Per receptor plots
    # -------------------------
    if not args.skip_per_receptor:
        for receptor, sub_sum in summary.groupby("Receptor"):
            # filter ligands with enough scores
            sub_sum = sub_sum[sub_sum["n_scores"] >= args.min_scores].copy()
            if sub_sum.empty:
                continue

            # choose top N by best_affinity (most negative first)
            top = sub_sum.sort_values(["best_affinity", "median_affinity"], ascending=[True, True]).head(args.top_n)

            top_ligands = top["LigandBase"].tolist()
            sub_df = df[(df["Receptor"] == receptor) & (df["LigandBase"].isin(top_ligands))].copy()

            receptor_dir = outdir / str(receptor)
            receptor_dir.mkdir(parents=True, exist_ok=True)

            # Combined boxplot
            data_list = []
            labels = []
            for lig in top_ligands:
                vals = sub_df.loc[sub_df["LigandBase"] == lig, "Binding_Affinity"].values
                if len(vals) >= args.min_scores:
                    data_list.append(vals)
                    labels.append(lig)

            if data_list:
                out_box = receptor_dir / f"{receptor}_TOP{min(args.top_n, len(labels))}_combined_boxplot.png"
                plot_multi_box(
                    data_list,
                    labels,
                    title=f"{receptor}: Docking score spread (Top {min(args.top_n, len(labels))} ligands)",
                    outpath=out_box
                )
                print(f"📦 Boxplot: {out_box}")

                # Pooled histogram (all top ligands combined)
                pooled = np.concatenate([np.asarray(v, dtype=float) for v in data_list])
                out_hist = receptor_dir / f"{receptor}_TOP{min(args.top_n, len(labels))}_pooled_histogram.png"
                plot_histogram(
                    pooled,
                    title=f"{receptor}: Pooled docking score distribution (Top {min(args.top_n, len(labels))} ligands)",
                    outpath=out_hist,
                    bins=args.bins
                )
                print(f"📦 Histogram: {out_hist}")

                # Overlay histogram (each ligand shown separately)
                values_by_label = {labels[i]: data_list[i] for i in range(len(labels))}

                out_overlay = receptor_dir / f"{receptor}_TOP{min(args.top_n, len(labels))}_overlay_histograms.png"
                plot_overlay_histograms(
                    values_by_label,
                    title=f"{receptor}: Docking score distributions (Top {min(args.top_n, len(labels))} ligands)",
                    outpath=out_overlay,
                    bins=args.bins,
                    density=True,
                    alpha=0.35
                )
                print(f"📦 Overlay histogram: {out_overlay}")


    # -------------------------
    # Per ligand plots
    # -------------------------
    if not args.skip_per_ligand:
        for receptor, sub in df.groupby("Receptor"):
            receptor_dir = outdir / str(receptor)
            lig_dir = receptor_dir / "Ligands"
            lig_dir.mkdir(parents=True, exist_ok=True)

            for lig, lig_df in sub.groupby("LigandBase"):
                vals = lig_df["Binding_Affinity"].dropna().values.astype(float)
                if len(vals) < args.min_scores:
                    continue

                safe_lig = re.sub(r"[^A-Za-z0-9._-]+", "_", str(lig))

                out_hist = lig_dir / f"{receptor}__{safe_lig}__hist.png"
                plot_histogram(
                    vals,
                    title=f"{receptor} / {lig}: Docking score distribution (all conformers + poses)",
                    outpath=out_hist,
                    bins=args.bins
                )

                out_box = lig_dir / f"{receptor}__{safe_lig}__box.png"
                plot_single_box(
                    vals,
                    title=f"{receptor} / {lig}: Docking score spread (box/whisker)",
                    outpath=out_box
                )

        print(f"✅ Per-ligand plots written under: {outdir}/<RECEPTOR>/Ligands/")

    print("\n🎉 Done.")
    print(f"Outputs: {outdir}")
    print(f"Summary CSV: {summary_path}")


if __name__ == "__main__":
    main()
