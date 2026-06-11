#!/usr/bin/env python3
from __future__ import annotations

"""
Expanded docking analytics and report generation for Vina result CSVs.

The script preserves the original histogram/boxplot workflow while adding:
- richer ligand/state summary statistics
- consensus ranking
- hit-rate and state heatmaps
- optional PDF report generation
- a machine-readable manifest of generated outputs
"""

import argparse
from collections import defaultdict
import json
import math
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    HAVE_MATPLOTLIB = True
except Exception:  # pragma: no cover - optional plotting dependency
    plt = None
    PdfPages = None  # type: ignore[assignment]
    HAVE_MATPLOTLIB = False

from ligand_manifest import CHEMICAL_METADATA_COLUMNS, merge_ligand_metadata
from ligand_naming import parse_ligand_variant


DEFAULT_HIT_THRESHOLDS = [-8.0, -8.5, -9.0, -9.5, -10.0]
DEFAULT_STATE_HEATMAP_METRIC = "best"
DEFAULT_RANKING_METRIC = "consensus"
DEFAULT_PDF_TOP_LIGANDS = 10


def normalize_ligand_id(ligand: str, outfile: Optional[str] = None) -> str:
    lig = (ligand or "").strip()
    if not lig and outfile:
        try:
            lig = Path(str(outfile)).parts[-2]
        except Exception:
            lig = ""
    if not lig:
        return "UNKNOWN"
    return str(parse_ligand_variant(lig)["LigandBase"])


def ligand_variant_from_row(ligand: str, outfile: Optional[str] = None) -> str:
    lig = (ligand or "").strip()
    if lig:
        return lig
    if outfile:
        try:
            return Path(str(outfile)).parts[-2]
        except Exception:
            return ""
    return ""


def safe_float(x: Any, default: float = 1e9) -> float:
    try:
        return float(x)
    except Exception:
        return default


def sanitize_filename(name: str) -> str:
    cleaned = "".join(c if c.isalnum() or c in ("-", "_", ".") else "_" for c in str(name or "item"))
    cleaned = cleaned.strip("._")
    return cleaned or "item"


def threshold_label(threshold: float) -> str:
    rendered = f"{threshold:g}"
    return rendered.replace(".", "_")


def find_csv_candidates(root: Path, pattern: str) -> List[Path]:
    root = root.expanduser().resolve()
    candidates = list(root.glob(pattern)) + list(root.glob(f"**/{pattern}"))
    unique: Dict[str, Path] = {}
    for path in candidates:
        try:
            unique[str(path.resolve())] = path.resolve()
        except Exception:
            continue
    return sorted(unique.values(), key=lambda path: path.stat().st_mtime, reverse=True)


def choose_path(prompt: str, items: Sequence[Path]) -> Path:
    print(f"\n{prompt}")
    for idx, path in enumerate(items):
        mtime = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(path.stat().st_mtime))
        print(f" [{idx}] {path}   (mtime: {mtime})")
    while True:
        raw = input("Enter index: ").strip()
        try:
            idx = int(raw)
            if 0 <= idx < len(items):
                return items[idx]
        except ValueError:
            pass
        print("Invalid choice. Try again.")


def _boxplot(ax, data: Sequence[Sequence[float]], labels: Sequence[str], showfliers: bool = False):
    try:
        return ax.boxplot(data, tick_labels=labels, showfliers=showfliers)
    except TypeError:  # pragma: no cover - compatibility with older matplotlib
        return ax.boxplot(data, labels=labels, showfliers=showfliers)


def _annotate_best_mean_worst(ax, values: np.ndarray):
    if values.size == 0:
        return
    best = float(np.min(values))
    worst = float(np.max(values))
    mean = float(np.mean(values))
    ax.axvline(best, linestyle="--", linewidth=1.5, label=f"Best: {best:.3f}")
    ax.axvline(mean, linestyle=":", linewidth=2.0, label=f"Mean: {mean:.3f}")
    ax.axvline(worst, linestyle="--", linewidth=1.5, label=f"Worst: {worst:.3f}")
    ymax = ax.get_ylim()[1]
    ax.text(best, ymax * 0.95, f"best\n{best:.3f}", ha="center", va="top", fontsize=9)
    ax.text(mean, ymax * 0.90, f"mean\n{mean:.3f}", ha="center", va="top", fontsize=9)
    ax.text(worst, ymax * 0.95, f"worst\n{worst:.3f}", ha="center", va="top", fontsize=9)


def plot_histogram(values: np.ndarray, title: str, outpath: Path, bins: int = 30):
    if not HAVE_MATPLOTLIB:
        raise RuntimeError("matplotlib is required for plotting")
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


def plot_overlay_histograms(
    values_by_label: Dict[str, np.ndarray],
    title: str,
    outpath: Path,
    bins: int = 30,
    density: bool = True,
    alpha: float = 0.35,
    max_legend: int = 25,
):
    if not HAVE_MATPLOTLIB:
        raise RuntimeError("matplotlib is required for plotting")
    labels = list(values_by_label.keys())
    all_vals = np.concatenate([np.asarray(values_by_label[label], dtype=float) for label in labels])
    if all_vals.size == 0:
        return
    bin_edges = np.linspace(float(np.min(all_vals)), float(np.max(all_vals)), bins + 1)
    plt.figure(figsize=(10, 6.5))
    ax = plt.gca()
    for label in labels:
        vals = np.asarray(values_by_label[label], dtype=float)
        if vals.size == 0:
            continue
        ax.hist(vals, bins=bin_edges, density=density, alpha=alpha, edgecolor="black", linewidth=0.3, label=label)
    ax.set_xlabel("Docking affinity (kcal/mol)")
    ax.set_ylabel("Density" if density else "Count")
    ax.set_title(title)
    _annotate_best_mean_worst(ax, all_vals)
    if len(labels) <= max_legend:
        ax.legend(frameon=True, fontsize=9, ncol=2)
    else:
        legend = ax.legend()
        if legend:
            legend.remove()
        ax.text(
            0.99,
            0.98,
            f"Legend hidden (N={len(labels)}). Use summary CSV for labels.",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
        )
    plt.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(outpath, dpi=300)
    plt.close()


def plot_single_box(values: np.ndarray, title: str, outpath: Path):
    if not HAVE_MATPLOTLIB:
        raise RuntimeError("matplotlib is required for plotting")
    plt.figure(figsize=(7.5, 6))
    ax = plt.gca()
    _boxplot(ax, [values], labels=[""], showfliers=False)
    median = float(np.median(values))
    mean = float(np.mean(values))
    q3 = float(np.quantile(values, 0.75))
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


def plot_multi_box(data_list: Sequence[Sequence[float]], labels: Sequence[str], title: str, outpath: Path):
    if not HAVE_MATPLOTLIB:
        raise RuntimeError("matplotlib is required for plotting")
    plt.figure(figsize=(max(10, 0.45 * len(labels)), 6.5))
    ax = plt.gca()
    _boxplot(ax, data_list, labels=labels, showfliers=False)
    y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
    y_offset = 0.02 * y_range
    for idx, vals in enumerate(data_list, start=1):
        arr = np.asarray(vals, dtype=float)
        mean = float(np.mean(arr))
        q3 = float(np.quantile(arr, 0.75))
        ax.text(idx, q3 + y_offset, f"μ={mean:.2f}", ha="center", va="bottom", fontsize=8)
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


def plot_ranking_bar(summary_df: pd.DataFrame, receptor: str, top_n: int, outpath: Path):
    if not HAVE_MATPLOTLIB or summary_df.empty:
        return
    metrics = ["best_affinity", "median_affinity", "top5_mean_affinity", "top10_mean_affinity"]
    fig, ax = plt.subplots(figsize=(max(10, 1.2 * len(summary_df)), 6.5))
    x = np.arange(len(summary_df))
    width = 0.18
    for idx, metric in enumerate(metrics):
        ax.bar(x + (idx - 1.5) * width, summary_df[metric].astype(float).values, width=width, label=metric.replace("_affinity", ""))
    ax.set_xticks(x)
    ax.set_xticklabels(summary_df["LigandBase"], rotation=45, ha="right")
    ax.set_ylabel("Docking affinity (kcal/mol)")
    ax.set_title(f"{receptor}: ligand ranking metrics (Top {top_n})")
    ax.legend(frameon=True, fontsize=9)
    ax.grid(axis="y", linestyle="--", alpha=0.25)
    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def plot_consensus_rank(summary_df: pd.DataFrame, receptor: str, top_n: int, outpath: Path):
    if not HAVE_MATPLOTLIB or summary_df.empty:
        return
    fig, ax = plt.subplots(figsize=(max(10, 1.0 * len(summary_df)), 6))
    ax.bar(summary_df["LigandBase"], summary_df["consensus_rank_score"].astype(float).values, color="#4C78A8")
    ax.set_ylabel("Consensus rank score (lower is better)")
    ax.set_title(f"{receptor}: consensus ligand ranking (Top {top_n})")
    ax.grid(axis="y", linestyle="--", alpha=0.25)
    plt.xticks(rotation=45, ha="right")
    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def plot_heatmap(df: pd.DataFrame, title: str, outpath: Path, value_format: str = ".2f", cmap: str = "viridis"):
    if not HAVE_MATPLOTLIB or df.empty:
        return
    fig, ax = plt.subplots(figsize=(max(8, 0.8 * len(df.columns)), max(4, 0.45 * len(df.index))))
    values = df.astype(float).values
    im = ax.imshow(values, aspect="auto", cmap=cmap)
    ax.set_xticks(np.arange(len(df.columns)))
    ax.set_xticklabels(list(df.columns), rotation=45, ha="right")
    ax.set_yticks(np.arange(len(df.index)))
    ax.set_yticklabels(list(df.index))
    ax.set_title(title)
    for row_idx in range(values.shape[0]):
        for col_idx in range(values.shape[1]):
            val = values[row_idx, col_idx]
            if np.isnan(val):
                label = ""
            else:
                label = format(val, value_format)
            ax.text(col_idx, row_idx, label, ha="center", va="center", color="white" if not np.isnan(val) and val < np.nanmean(values) else "black", fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def plot_state_colored_histogram(values_by_state: Dict[str, np.ndarray], title: str, outpath: Path, bins: int = 30):
    plot_overlay_histograms(values_by_state, title, outpath, bins=bins, density=False, alpha=0.45, max_legend=20)


def _safe_float_or_nan(value: Any) -> float:
    try:
        return float(value)
    except Exception:
        return float("nan")


def _nonempty_unique_count(series: pd.Series) -> int:
    if series is None:
        return 0
    cleaned = series.fillna("").astype(str)
    cleaned = cleaned[cleaned != ""]
    return int(cleaned.nunique())


def _first_nonempty(series: pd.Series) -> str:
    for value in series:
        if pd.notna(value) and str(value).strip() != "":
            return str(value)
    return ""


def top_k_values(values: np.ndarray, k: int) -> np.ndarray:
    if values.size == 0:
        return np.asarray([], dtype=float)
    return np.sort(values)[: min(k, values.size)]


def top_k_mean(values: np.ndarray, k: int) -> float:
    selected = top_k_values(values, k)
    return float(np.mean(selected)) if selected.size else float("nan")


def top_k_median(values: np.ndarray, k: int) -> float:
    selected = top_k_values(values, k)
    return float(np.median(selected)) if selected.size else float("nan")


def affinity_metrics(values: Sequence[float], hit_thresholds: Sequence[float] = DEFAULT_HIT_THRESHOLDS) -> Dict[str, Any]:
    arr = np.asarray(list(values), dtype=float)
    arr = arr[~np.isnan(arr)]
    if arr.size == 0:
        return {}
    arr_sorted = np.sort(arr)
    quantiles = {
        "q05_affinity": float(np.quantile(arr_sorted, 0.05)),
        "q10_affinity": float(np.quantile(arr_sorted, 0.10)),
        "q25_affinity": float(np.quantile(arr_sorted, 0.25)),
        "q75_affinity": float(np.quantile(arr_sorted, 0.75)),
        "q90_affinity": float(np.quantile(arr_sorted, 0.90)),
        "q95_affinity": float(np.quantile(arr_sorted, 0.95)),
    }
    best = float(arr_sorted[0])
    worst = float(arr_sorted[-1])
    median = float(np.median(arr_sorted))
    mean = float(np.mean(arr_sorted))
    std = float(np.std(arr_sorted, ddof=1)) if arr_sorted.size > 1 else 0.0
    out: Dict[str, Any] = {
        "n_scores": int(arr_sorted.size),
        "best_affinity": best,
        "worst_affinity": worst,
        "mean_affinity": mean,
        "median_affinity": median,
        "std_affinity": std,
        **quantiles,
        "iqr_affinity": quantiles["q75_affinity"] - quantiles["q25_affinity"],
        "score_range": worst - best,
        "top3_mean_affinity": top_k_mean(arr_sorted, 3),
        "top5_mean_affinity": top_k_mean(arr_sorted, 5),
        "top10_mean_affinity": top_k_mean(arr_sorted, 10),
        "top3_median_affinity": top_k_median(arr_sorted, 3),
        "top5_median_affinity": top_k_median(arr_sorted, 5),
        "top10_median_affinity": top_k_median(arr_sorted, 10),
        "outlier_gap": median - best,
    }
    # More negative is better, so a lower robustness score indicates stronger and more consistent performance.
    out["robustness_score"] = out["top5_mean_affinity"] + out["outlier_gap"]
    for threshold in hit_thresholds:
        label = threshold_label(threshold)
        hit_count = int(np.sum(arr_sorted <= threshold))
        out[f"hit_count_below_{label}"] = hit_count
        out[f"hit_rate_below_{label}"] = float(hit_count / arr_sorted.size)
    return out


def _best_row_fields(group: pd.DataFrame) -> Dict[str, Any]:
    ordered = group.sort_values(["Binding_Affinity", "Pose"], ascending=[True, True]).reset_index(drop=True)
    best = ordered.iloc[0]
    result = {
        "best_ligand_variant": best.get("LigandVariant", ""),
        "best_state_tag": best.get("StateTag", ""),
        "best_protomer_tag": best.get("ProtomerTag", ""),
        "best_tautomer_tag": best.get("TautomerTag", ""),
        "best_conformer_tag": best.get("ConformerTag", ""),
        "best_vina_pose": best.get("Pose", ""),
        "best_outfile": best.get("OutFile", ""),
        "StateCanonicalSMILES": _first_nonempty(ordered.get("StateCanonicalSMILES", pd.Series(dtype=object))),
        "StateIsomericSMILES": _first_nonempty(ordered.get("StateIsomericSMILES", pd.Series(dtype=object))),
        "Formula": _first_nonempty(ordered.get("Formula", pd.Series(dtype=object))),
        "FormalCharge": _first_nonempty(ordered.get("FormalCharge", pd.Series(dtype=object))),
    }
    return result


def _summarize_group(group: pd.DataFrame, hit_thresholds: Sequence[float] = DEFAULT_HIT_THRESHOLDS) -> pd.Series:
    ordered = group.sort_values(["Binding_Affinity", "Pose"], ascending=[True, True]).reset_index(drop=True)
    values = ordered["Binding_Affinity"].astype(float).values
    summary = affinity_metrics(values, hit_thresholds)
    summary.update(
        {
            "n_variants": _nonempty_unique_count(ordered.get("LigandVariant", pd.Series(dtype=object))),
            "n_states": _nonempty_unique_count(ordered.get("StateTag", pd.Series(dtype=object))),
            "n_protomers": _nonempty_unique_count(ordered.get("ProtomerTag", pd.Series(dtype=object))),
            "n_tautomers": _nonempty_unique_count(ordered.get("TautomerTag", pd.Series(dtype=object))),
            "n_conformers": _nonempty_unique_count(ordered.get("ConformerTag", pd.Series(dtype=object))),
        }
    )
    summary.update(_best_row_fields(ordered))
    return pd.Series(summary)


def add_consensus_ranks(summary_df: pd.DataFrame, hit_thresholds: Sequence[float] = DEFAULT_HIT_THRESHOLDS) -> pd.DataFrame:
    if summary_df.empty:
        return summary_df
    out = summary_df.copy()
    hit_col = f"hit_rate_below_{threshold_label(-9.0)}"
    for receptor, idx in out.groupby("Receptor").groups.items():
        block = out.loc[list(idx)].copy()
        out.loc[list(idx), "rank_best_affinity"] = block["best_affinity"].rank(method="min", ascending=True)
        out.loc[list(idx), "rank_median_affinity"] = block["median_affinity"].rank(method="min", ascending=True)
        out.loc[list(idx), "rank_top5_mean_affinity"] = block["top5_mean_affinity"].rank(method="min", ascending=True)
        out.loc[list(idx), "rank_top10_mean_affinity"] = block["top10_mean_affinity"].rank(method="min", ascending=True)
        if hit_col in block.columns:
            out.loc[list(idx), "rank_hit_rate_below_-9"] = block[hit_col].rank(method="min", ascending=False)
        else:
            out.loc[list(idx), "rank_hit_rate_below_-9"] = np.nan
        out.loc[list(idx), "rank_outlier_gap"] = block["outlier_gap"].rank(method="min", ascending=True)

    rank_cols = [
        "rank_best_affinity",
        "rank_median_affinity",
        "rank_top5_mean_affinity",
        "rank_top10_mean_affinity",
        "rank_hit_rate_below_-9",
        "rank_outlier_gap",
    ]
    out["consensus_rank_score"] = out[rank_cols].mean(axis=1, skipna=True)
    for receptor, idx in out.groupby("Receptor").groups.items():
        block = out.loc[list(idx)]
        out.loc[list(idx), "consensus_rank"] = block["consensus_rank_score"].rank(method="min", ascending=True)
    return out


def build_summary(df: pd.DataFrame, hit_thresholds: Sequence[float] = DEFAULT_HIT_THRESHOLDS) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for (receptor, ligand_base), group in df.groupby(["Receptor", "LigandBase"], dropna=False):
        row = {"Receptor": receptor, "LigandBase": ligand_base}
        row.update(_summarize_group(group, hit_thresholds).to_dict())
        rows.append(row)
    summary = pd.DataFrame(rows)
    if summary.empty:
        return summary
    summary = add_consensus_ranks(summary, hit_thresholds)
    summary = summary.sort_values(["Receptor", "consensus_rank", "best_affinity", "median_affinity"], ascending=[True, True, True, True]).reset_index(drop=True)
    return summary


def summarize_by_groups(
    df: pd.DataFrame,
    group_cols: Sequence[str],
    hit_thresholds: Sequence[float] = DEFAULT_HIT_THRESHOLDS,
) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for keys, group in df.groupby(list(group_cols), dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        record = {col: key for col, key in zip(group_cols, keys)}
        record.update(_summarize_group(group, hit_thresholds).to_dict())
        rows.append(record)
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    return result.sort_values(["Receptor", "best_affinity", "median_affinity"], ascending=[True, True, True]).reset_index(drop=True)


def write_state_summary_tables(
    df: pd.DataFrame,
    outdir: Path,
    stem: str,
    hit_thresholds: Sequence[float] = DEFAULT_HIT_THRESHOLDS,
) -> Dict[str, Path]:
    summary_specs = {
        "ligand_base_summary": ["Receptor", "LigandBase"],
        "tautomer_summary": ["Receptor", "LigandBase", "TautomerTag"],
        "protomer_summary": ["Receptor", "LigandBase", "ProtomerTag"],
        "state_summary": ["Receptor", "LigandBase", "StateTag"],
        "variant_summary_state": ["Receptor", "LigandBase", "LigandVariant"],
    }
    written: Dict[str, Path] = {}
    summary_tables_dir = outdir / "summary_tables"
    summary_tables_dir.mkdir(parents=True, exist_ok=True)
    for label, group_cols in summary_specs.items():
        available = [col for col in group_cols if col in df.columns]
        if len(available) != len(group_cols):
            continue
        filtered = df.copy()
        for col in group_cols[1:]:
            filtered = filtered[filtered[col].fillna("").astype(str) != ""]
        if filtered.empty:
            continue
        summary_df = summarize_by_groups(filtered, group_cols, hit_thresholds)
        out_path = summary_tables_dir / f"{stem}_{label}.csv"
        summary_df.to_csv(out_path, index=False)
        compatibility_path = outdir / f"{stem}_{label}.csv"
        summary_df.to_csv(compatibility_path, index=False)
        written[label] = out_path
    return written


def _top_receptor_summary(summary_df: pd.DataFrame, receptor: str, top_n: int, ranking_metric: str) -> pd.DataFrame:
    sub = summary_df[summary_df["Receptor"] == receptor].copy()
    if sub.empty:
        return sub
    metric_map = {
        "consensus": "consensus_rank",
        "best": "best_affinity",
        "median": "median_affinity",
        "top5_mean": "top5_mean_affinity",
        "top10_mean": "top10_mean_affinity",
    }
    sort_col = metric_map.get(ranking_metric, "consensus_rank")
    ascending = True
    sub = sub.sort_values([sort_col, "best_affinity", "median_affinity"], ascending=[ascending, True, True]).head(top_n)
    return sub


def write_summary_csv(df: pd.DataFrame, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def write_manifest_json(path: Path, payload: Dict[str, Any]):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))


def _fig_text_page(lines: Sequence[str], title: str):
    fig = plt.figure(figsize=(8.27, 11.69))
    fig.suptitle(title, fontsize=16, fontweight="bold", y=0.98)
    fig.text(0.06, 0.94, "\n".join(lines), va="top", ha="left", fontsize=10, family="monospace")
    return fig


def _fig_table_page(df: pd.DataFrame, title: str, max_rows: int = 18):
    fig, ax = plt.subplots(figsize=(11.0, 8.5))
    ax.axis("off")
    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)
    table = ax.table(cellText=df.head(max_rows).values, colLabels=df.columns, loc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.25)
    fig.tight_layout()
    return fig


def _fig_image_page(image_path: Path, title: str):
    img = plt.imread(image_path)
    fig, ax = plt.subplots(figsize=(11.0, 8.5))
    ax.imshow(img)
    ax.axis("off")
    ax.set_title(title, fontsize=14, fontweight="bold")
    fig.tight_layout()
    return fig


def write_pdf_report(
    pdf_path: Path,
    csv_path: Path,
    summary_df: pd.DataFrame,
    variant_summary: pd.DataFrame,
    receptor_plot_index: Dict[str, Dict[str, Path]],
    receptor_ligand_plots: Dict[str, Dict[str, Dict[str, Path]]],
    manifest_payload: Dict[str, Any],
    top_n: int,
    hit_thresholds: Sequence[float],
    pdf_top_ligands: int,
):
    if not HAVE_MATPLOTLIB or PdfPages is None:
        raise RuntimeError("matplotlib with PdfPages is required for PDF generation")
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(pdf_path) as pdf:
        cover_lines = [
            f"Input CSV: {csv_path}",
            f"Generated: {manifest_payload['generated_at']}",
            f"Receptors: {summary_df['Receptor'].nunique() if not summary_df.empty else 0}",
            f"Ligand bases: {summary_df['LigandBase'].nunique() if not summary_df.empty else 0}",
            f"Ligand variants: {variant_summary['LigandVariant'].nunique() if not variant_summary.empty else 0}",
            f"Total docking scores: {manifest_payload.get('total_scores', 0)}",
            f"Top-N setting: {top_n}",
            f"Hit thresholds: {', '.join(str(th) for th in hit_thresholds)}",
        ]
        pdf.savefig(_fig_text_page(cover_lines, "Docking Analytics Report"))
        plt.close("all")

        exec_rows: List[pd.DataFrame] = []
        hit_rate_col = f"hit_rate_below_{threshold_label(-9.0)}"
        for receptor, group in summary_df.groupby("Receptor"):
            cols = [
                "Receptor",
                "consensus_rank",
                "LigandBase",
                "best_affinity",
                "median_affinity",
                "top5_mean_affinity",
                "top10_mean_affinity",
                "outlier_gap",
                "best_ligand_variant",
                "best_tautomer_tag",
                "best_conformer_tag",
                "best_vina_pose",
            ]
            if hit_rate_col in group.columns:
                cols.insert(7, hit_rate_col)
            exec_rows.append(group.sort_values(["consensus_rank", "best_affinity"]).head(min(5, len(group)))[cols])
        if exec_rows:
            executive = pd.concat(exec_rows, ignore_index=True)
            pdf.savefig(_fig_table_page(executive, "Executive Summary"))
            plt.close("all")

        preferred_plot_order = [
            "combined_boxplot",
            "ranking_bar",
            "hit_rate_heatmap",
            "overlay_histograms",
            "pooled_histogram",
            "state_heatmap",
            "consensus_rank",
            "tautomer_summary",
        ]
        for receptor, plots in receptor_plot_index.items():
            for key in preferred_plot_order:
                image_path = plots.get(key)
                if image_path and image_path.exists():
                    pdf.savefig(_fig_image_page(image_path, f"{receptor}: {key.replace('_', ' ').title()}"))
                    plt.close("all")

        for receptor, ligands in receptor_ligand_plots.items():
            count = 0
            for ligand, plot_map in ligands.items():
                for plot_key in ("hist", "box", "state_hist"):
                    image_path = plot_map.get(plot_key)
                    if image_path and image_path.exists():
                        pdf.savefig(_fig_image_page(image_path, f"{receptor} / {ligand}: {plot_key.replace('_', ' ').title()}"))
                        plt.close("all")
                count += 1
                if count >= pdf_top_ligands:
                    break

        output_lines = [
            f"Ligand summary CSV: {manifest_payload['summary_csvs'].get('ligand_summary', '')}",
            f"Variant summary CSV: {manifest_payload['summary_csvs'].get('variant_summary', '')}",
            f"State summary CSVs: {', '.join(manifest_payload['summary_csvs'].get('state_summaries', []))}",
            f"Plot directory: {manifest_payload['output_dir']}",
            f"PDF report path: {pdf_path}",
        ]
        pdf.savefig(_fig_text_page(output_lines, "Output Paths"))
        plt.close("all")


def normalize_rows(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    colmap = {c.lower(): c for c in out.columns}
    required = ["receptor", "ligand", "binding_affinity"]
    for key in required:
        if key not in colmap:
            raise ValueError(f"Missing required column: {key} (found columns: {list(out.columns)})")
    receptor_col = colmap["receptor"]
    ligand_col = colmap["ligand"]
    affinity_col = colmap["binding_affinity"]
    outfile_col = colmap.get("outfile")

    out = out.rename(columns={receptor_col: "Receptor", ligand_col: "Ligand", affinity_col: "Binding_Affinity"})
    out["Binding_Affinity"] = pd.to_numeric(out["Binding_Affinity"], errors="coerce")
    out = out.dropna(subset=["Receptor", "Binding_Affinity"])
    if outfile_col:
        out["OutFile"] = out[outfile_col].fillna("").astype(str)
    elif "OutFile" not in out.columns:
        out["OutFile"] = ""
    if "LigandVariant" not in out.columns:
        out["LigandVariant"] = [
            ligand_variant_from_row(lig, out_file)
            for lig, out_file in zip(out["Ligand"].astype(str), out["OutFile"].astype(str))
        ]
    else:
        out["LigandVariant"] = out["LigandVariant"].fillna("").astype(str)

    normalized: List[Dict[str, Any]] = []
    for _, row in out.iterrows():
        variant = row.get("LigandVariant", "") or ligand_variant_from_row(row.get("Ligand", ""), row.get("OutFile", ""))
        normalized.append(merge_ligand_metadata(str(variant), base_row=row.to_dict()))
    cooked = pd.DataFrame(normalized)
    cooked["Binding_Affinity"] = pd.to_numeric(cooked["Binding_Affinity"], errors="coerce")
    cooked = cooked.dropna(subset=["Receptor", "Binding_Affinity"])
    for col in CHEMICAL_METADATA_COLUMNS:
        if col not in cooked.columns:
            cooked[col] = ""
    for col in ("ProtomerTag", "TautomerTag", "ConformerTag", "StateTag", "LigandBase", "LigandVariant", "Ligand"):
        if col not in cooked.columns:
            cooked[col] = ""
        cooked[col] = cooked[col].fillna("").astype(str)
    if "Pose" in cooked.columns:
        cooked["Pose"] = pd.to_numeric(cooked["Pose"], errors="coerce").fillna(0).astype(int)
    else:
        cooked["Pose"] = 0
    return cooked


def generate_outputs(
    df: pd.DataFrame,
    csv_path: Path,
    outdir: Path,
    top_n: int,
    bins: int,
    min_scores: int,
    hit_thresholds: Sequence[float],
    ranking_metric: str,
    write_state_summaries: bool,
    write_ranking_plots: bool,
    write_hit_rate_plots: bool,
    write_state_heatmaps: bool,
    write_state_plots: bool,
    state_heatmap_metric: str,
    write_pdf_report_flag: bool,
    pdf_top_ligands: int,
    skip_per_ligand: bool,
    skip_per_receptor: bool,
) -> Dict[str, Any]:
    outdir.mkdir(parents=True, exist_ok=True)
    summary_tables_dir = outdir / "summary_tables"
    summary_tables_dir.mkdir(parents=True, exist_ok=True)

    summary_df = build_summary(df, hit_thresholds)
    ligand_summary_path = summary_tables_dir / f"{csv_path.stem}_ligand_summary.csv"
    write_summary_csv(summary_df, ligand_summary_path)
    write_summary_csv(summary_df, outdir / f"{csv_path.stem}_ligand_summary.csv")

    variant_summary = summarize_by_groups(df, ["Receptor", "LigandBase", "LigandVariant"], hit_thresholds)
    variant_summary_path = summary_tables_dir / f"{csv_path.stem}_variant_summary.csv"
    write_summary_csv(variant_summary, variant_summary_path)
    write_summary_csv(variant_summary, outdir / f"{csv_path.stem}_variant_summary.csv")

    state_summary_paths: Dict[str, Path] = {}
    if write_state_summaries:
        state_summary_paths = write_state_summary_tables(df, outdir, csv_path.stem, hit_thresholds)

    receptor_plot_index: Dict[str, Dict[str, Path]] = defaultdict(dict)
    receptor_ligand_plots: Dict[str, Dict[str, Dict[str, Path]]] = defaultdict(dict)
    warnings: List[str] = []

    for receptor in sorted(summary_df["Receptor"].unique()) if not summary_df.empty else []:
        receptor_dir = outdir / str(receptor)
        receptor_dir.mkdir(parents=True, exist_ok=True)
        top_summary = _top_receptor_summary(summary_df, receptor, top_n, ranking_metric)
        if top_summary.empty:
            continue
        top_ligands = top_summary["LigandBase"].tolist()
        sub_df = df[(df["Receptor"] == receptor) & (df["LigandBase"].isin(top_ligands))].copy()

        data_list: List[np.ndarray] = []
        labels: List[str] = []
        for ligand in top_ligands:
            vals = sub_df.loc[sub_df["LigandBase"] == ligand, "Binding_Affinity"].dropna().values.astype(float)
            if len(vals) >= min_scores:
                data_list.append(vals)
                labels.append(ligand)
        if data_list and HAVE_MATPLOTLIB and not skip_per_receptor:
            limit = min(top_n, len(labels))
            box_path = receptor_dir / f"{receptor}_TOP{limit}_combined_boxplot.png"
            plot_multi_box(data_list, labels, f"{receptor}: Docking score spread (Top {limit} ligands)", box_path)
            receptor_plot_index[receptor]["combined_boxplot"] = box_path

            pooled = np.concatenate([np.asarray(vals, dtype=float) for vals in data_list])
            pooled_path = receptor_dir / f"{receptor}_TOP{limit}_pooled_histogram.png"
            plot_histogram(pooled, f"{receptor}: Pooled docking score distribution (Top {limit} ligands)", pooled_path, bins=bins)
            receptor_plot_index[receptor]["pooled_histogram"] = pooled_path

            overlay_path = receptor_dir / f"{receptor}_TOP{limit}_overlay_histograms.png"
            plot_overlay_histograms({labels[idx]: data_list[idx] for idx in range(len(labels))}, f"{receptor}: Docking score distributions (Top {limit} ligands)", overlay_path, bins=bins)
            receptor_plot_index[receptor]["overlay_histograms"] = overlay_path

            if write_ranking_plots:
                rank_path = receptor_dir / f"{receptor}_TOP{limit}_ranking_bar.png"
                plot_ranking_bar(top_summary, receptor, limit, rank_path)
                receptor_plot_index[receptor]["ranking_bar"] = rank_path

                consensus_path = receptor_dir / f"{receptor}_TOP{limit}_consensus_rank.png"
                plot_consensus_rank(top_summary, receptor, limit, consensus_path)
                receptor_plot_index[receptor]["consensus_rank"] = consensus_path

            if write_hit_rate_plots:
                hit_columns = [f"hit_rate_below_{threshold_label(th)}" for th in hit_thresholds]
                heatmap_df = top_summary.set_index("LigandBase")[hit_columns].rename(columns={col: col.replace("hit_rate_below_", "<= ") for col in hit_columns})
                heatmap_path = receptor_dir / f"{receptor}_TOP{limit}_hit_rate_heatmap.png"
                plot_heatmap(heatmap_df, f"{receptor}: hit-rate heatmap (Top {limit} ligands)", heatmap_path, value_format=".2f", cmap="magma")
                receptor_plot_index[receptor]["hit_rate_heatmap"] = heatmap_path

            if write_state_heatmaps:
                state_summary = summarize_by_groups(
                    sub_df[sub_df["StateTag"].fillna("").astype(str) != ""],
                    ["Receptor", "LigandBase", "StateTag"],
                    hit_thresholds,
                )
                if not state_summary.empty:
                    metric_col = {
                        "best": "best_affinity",
                        "median": "median_affinity",
                        "top5_mean": "top5_mean_affinity",
                    }[state_heatmap_metric]
                    pivot = state_summary.pivot(index="LigandBase", columns="StateTag", values=metric_col).reindex(top_ligands)
                    if not pivot.dropna(how="all").empty:
                        state_path = receptor_dir / f"{receptor}_TOP{limit}_state_heatmap.png"
                        plot_heatmap(pivot.fillna(np.nan), f"{receptor}: {state_heatmap_metric} state heatmap (Top {limit} ligands)", state_path, value_format=".2f", cmap="coolwarm_r")
                        receptor_plot_index[receptor]["state_heatmap"] = state_path

                tautomer_summary = summarize_by_groups(
                    sub_df[sub_df["TautomerTag"].fillna("").astype(str) != ""],
                    ["Receptor", "LigandBase", "TautomerTag"],
                    hit_thresholds,
                )
                if not tautomer_summary.empty and HAVE_MATPLOTLIB:
                    fig, ax = plt.subplots(figsize=(max(10, 1.1 * len(tautomer_summary)), 6.5))
                    labels_t = [f"{row.LigandBase}:{row.TautomerTag}" for row in tautomer_summary.itertuples()]
                    ax.bar(labels_t, tautomer_summary["best_affinity"].astype(float).values, color="#72B7B2")
                    ax.set_ylabel("Best affinity (kcal/mol)")
                    ax.set_title(f"{receptor}: tautomer best-score summary")
                    ax.grid(axis="y", linestyle="--", alpha=0.25)
                    plt.xticks(rotation=60, ha="right")
                    fig.tight_layout()
                    tautomer_path = receptor_dir / f"{receptor}_TOP{limit}_tautomer_summary.png"
                    fig.savefig(tautomer_path, dpi=300)
                    plt.close(fig)
                    receptor_plot_index[receptor]["tautomer_summary"] = tautomer_path

        if not skip_per_ligand:
            lig_dir = receptor_dir / "Ligands"
            lig_dir.mkdir(parents=True, exist_ok=True)
            for ligand, ligand_df in df[df["Receptor"] == receptor].groupby("LigandBase"):
                vals = ligand_df["Binding_Affinity"].dropna().values.astype(float)
                if len(vals) < min_scores or not HAVE_MATPLOTLIB:
                    continue
                safe_lig = sanitize_filename(str(ligand))
                hist_path = lig_dir / f"{receptor}__{safe_lig}__hist.png"
                plot_histogram(vals, f"{receptor} / {ligand}: Docking score distribution (all variants + Vina poses)", hist_path, bins=bins)
                box_path = lig_dir / f"{receptor}__{safe_lig}__box.png"
                plot_single_box(vals, f"{receptor} / {ligand}: Docking score spread (box/whisker)", box_path)
                receptor_ligand_plots[receptor][ligand] = {"hist": hist_path, "box": box_path}

                if write_state_plots:
                    state_df = ligand_df[ligand_df["StateTag"].fillna("").astype(str) != ""].copy()
                    if not state_df.empty:
                        values_by_state = {
                            state: block["Binding_Affinity"].dropna().values.astype(float)
                            for state, block in state_df.groupby("StateTag")
                        }
                        if values_by_state:
                            state_hist_path = lig_dir / f"{receptor}__{safe_lig}__state_hist.png"
                            plot_state_colored_histogram(values_by_state, f"{receptor} / {ligand}: state-level score distribution", state_hist_path, bins=bins)
                            receptor_ligand_plots[receptor][ligand]["state_hist"] = state_hist_path

    manifest_payload: Dict[str, Any] = {
        "input_csv": str(csv_path),
        "output_dir": str(outdir),
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "top_n": int(top_n),
        "bins": int(bins),
        "min_scores": int(min_scores),
        "hit_thresholds": list(hit_thresholds),
        "total_scores": int(len(df)),
        "summary_csvs": {
            "ligand_summary": str(ligand_summary_path),
            "variant_summary": str(variant_summary_path),
            "state_summaries": [str(path) for path in state_summary_paths.values()],
        },
        "plot_files": sorted({str(path) for plot_map in receptor_plot_index.values() for path in plot_map.values()} | {str(path) for ligands in receptor_ligand_plots.values() for plot_map in ligands.values() for path in plot_map.values()}),
        "pdf_report": "",
        "warnings": warnings,
    }

    if write_pdf_report_flag:
        if not HAVE_MATPLOTLIB or PdfPages is None:
            warnings.append("Skipping PDF report because matplotlib PdfPages is unavailable.")
        else:
            pdf_path = outdir / f"docking_report_{csv_path.stem}.pdf"
            write_pdf_report(
                pdf_path,
                csv_path,
                summary_df,
                variant_summary,
                receptor_plot_index,
                receptor_ligand_plots,
                manifest_payload,
                top_n,
                hit_thresholds,
                pdf_top_ligands,
            )
            manifest_payload["pdf_report"] = str(pdf_path)

    manifest_path = outdir / "graph_report_manifest.json"
    write_manifest_json(manifest_path, manifest_payload)
    manifest_payload["manifest_json"] = str(manifest_path)
    return {
        "summary": summary_df,
        "variant_summary": variant_summary,
        "state_summary_paths": state_summary_paths,
        "manifest": manifest_payload,
        "manifest_path": manifest_path,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot histograms, ranking plots, and PDF reports for Vina docking result CSVs.")
    parser.add_argument("--csv", default=None, help="Input docking CSV. Use --csv auto to pick the newest match automatically.")
    parser.add_argument("--root", default=".", help="Where to search for CSVs if --csv not provided (default: current dir).")
    parser.add_argument("--pattern", default="*vina_docking_scores_sorted.csv", help="Glob pattern used to find candidate CSVs.")
    parser.add_argument("--non-interactive", action="store_true", help="If --csv is missing, do not prompt; fail instead.")
    parser.add_argument("--outdir", default=None, help="Output directory (default: alongside CSV, in Plots/)")
    parser.add_argument("--top-n", type=int, default=25, help="Top N ligands per receptor for combined plots.")
    parser.add_argument("--bins", type=int, default=30, help="Histogram bins.")
    parser.add_argument("--min-scores", type=int, default=1, help="Skip ligands with fewer than this many scores.")
    parser.add_argument("--skip-per-ligand", action="store_true", help="Do not write per-ligand plots.")
    parser.add_argument("--skip-per-receptor", action="store_true", help="Do not write per-receptor combined plots.")
    parser.add_argument("--write-state-summaries", action="store_true", help="Write ligand/protomer/tautomer/state/variant summary CSVs.")
    parser.add_argument("--write-pdf-report", action="store_true", help="Write a PDF report that wraps key plots and summary tables.")
    parser.add_argument("--pdf-top-ligands", type=int, default=DEFAULT_PDF_TOP_LIGANDS, help="Top ligands per receptor to include in the PDF appendix.")
    parser.add_argument("--write-ranking-plots", dest="write_ranking_plots", action="store_true", default=True, help="Write ranking bar and consensus rank plots.")
    parser.add_argument("--skip-ranking-plots", dest="write_ranking_plots", action="store_false", help="Skip ranking bar and consensus rank plots.")
    parser.add_argument("--write-hit-rate-plots", dest="write_hit_rate_plots", action="store_true", default=True, help="Write hit-rate heatmaps.")
    parser.add_argument("--skip-hit-rate-plots", dest="write_hit_rate_plots", action="store_false", help="Skip hit-rate heatmaps.")
    parser.add_argument("--write-state-heatmaps", dest="write_state_heatmaps", action="store_true", default=True, help="Write state heatmaps when state metadata exists.")
    parser.add_argument("--skip-state-heatmaps", dest="write_state_heatmaps", action="store_false", help="Skip state heatmaps.")
    parser.add_argument("--write-state-plots", action="store_true", help="Write per-ligand state-colored histograms.")
    parser.add_argument("--state-heatmap-metric", choices=["best", "median", "top5_mean"], default=DEFAULT_STATE_HEATMAP_METRIC, help="Metric to display in the state heatmap.")
    parser.add_argument("--hit-thresholds", nargs="+", type=float, default=DEFAULT_HIT_THRESHOLDS, help="Affinity thresholds for hit counts and hit rates.")
    parser.add_argument("--ranking-metric", choices=["consensus", "best", "median", "top5_mean", "top10_mean"], default=DEFAULT_RANKING_METRIC, help="Metric used to choose Top-N ligands per receptor.")
    parser.add_argument("--skip-state-plots", action="store_true", help="Compatibility flag; state plots are still off by default unless --write-state-plots is used.")
    parser.add_argument("--top-n-states", type=int, default=20, help="Compatibility flag reserved for future state-plot filtering.")
    return parser.parse_args()


def main():
    args = parse_args()
    root = Path(args.root)
    if args.csv is None or str(args.csv).strip().lower() == "auto":
        candidates = find_csv_candidates(root, args.pattern)
        if not candidates:
            raise RuntimeError(f"No CSVs found under {root} matching pattern: {args.pattern}")
        if str(args.csv).strip().lower() == "auto":
            in_csv = candidates[0]
            print(f"✅ Auto-selected newest CSV: {in_csv}")
        else:
            if args.non_interactive:
                raise RuntimeError("Missing --csv and --non-interactive set. Provide --csv explicitly.")
            in_csv = choose_path("Select a docking CSV to plot:", candidates)
    else:
        in_csv = Path(args.csv).expanduser().resolve()
        if not in_csv.exists():
            raise FileNotFoundError(f"CSV not found: {in_csv}")

    outdir = Path(args.outdir).expanduser().resolve() if args.outdir else in_csv.parent / "Plots"
    raw_df = pd.read_csv(in_csv, low_memory=False)
    df = normalize_rows(raw_df)

    result = generate_outputs(
        df=df,
        csv_path=in_csv,
        outdir=outdir,
        top_n=args.top_n,
        bins=args.bins,
        min_scores=args.min_scores,
        hit_thresholds=args.hit_thresholds,
        ranking_metric=args.ranking_metric,
        write_state_summaries=args.write_state_summaries,
        write_ranking_plots=args.write_ranking_plots and not args.skip_per_receptor,
        write_hit_rate_plots=args.write_hit_rate_plots and not args.skip_per_receptor,
        write_state_heatmaps=args.write_state_heatmaps and not args.skip_per_receptor,
        write_state_plots=args.write_state_plots and not args.skip_state_plots,
        state_heatmap_metric=args.state_heatmap_metric,
        write_pdf_report_flag=args.write_pdf_report,
        pdf_top_ligands=args.pdf_top_ligands,
        skip_per_ligand=args.skip_per_ligand,
        skip_per_receptor=args.skip_per_receptor,
    )

    print(f"✅ Wrote summary: {outdir / 'summary_tables' / f'{in_csv.stem}_ligand_summary.csv'}")
    print(f"✅ Wrote variant summary: {outdir / 'summary_tables' / f'{in_csv.stem}_variant_summary.csv'}")
    if args.write_state_summaries:
        for path in result["state_summary_paths"].values():
            print(f"✅ Wrote state summary: {path}")
    if result["manifest"].get("pdf_report"):
        print(f"✅ Wrote PDF report: {result['manifest']['pdf_report']}")
    print(f"✅ Wrote run manifest: {result['manifest_path']}")
    print("\n🎉 Done.")
    print(f"Outputs: {outdir}")


if __name__ == "__main__":
    main()
