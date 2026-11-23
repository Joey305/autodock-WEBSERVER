from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


MD_ROOT = Path("/Users/jxs794/Library/CloudStorage/Box-Box/PROTAC_Degrader_DB_Project_2025/Kif20a/Manuscript Pieces/KIF20A_Final_Analysis/MolecularDynamics/UMF607-Z")
RESULTS = MD_ROOT / "Analysis_Results"
OUT = RESULTS / "Figure5_Replacement"


def rmsd_panel() -> None:
    protein = pd.read_csv(RESULTS / "Protein_RMSD.csv")
    ligand = pd.read_csv(RESULTS / "UMF607-Z_Ligand_RMSD.csv")
    frame = pd.DataFrame({
        "Time (ns)": protein["Time_ps"] / 1000,
        "Protein backbone RMSD": protein["Protein_RMSD"],
        "Ligand heavy-atom RMSD": ligand["Ligand_RMSD"],
    })
    smooth = frame.set_index("Time (ns)").rolling(1001, center=True, min_periods=1).mean().reset_index()
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(7.1, 5.2), constrained_layout=True)
    styles = [("Protein backbone RMSD", "#f26b21"), ("Ligand heavy-atom RMSD", "#006747")]
    for ax, (column, color) in zip(axes, styles):
        ax.plot(frame["Time (ns)"], frame[column], color=color, alpha=0.12, linewidth=0.45)
        ax.plot(smooth["Time (ns)"], smooth[column], color=color, linewidth=1.8)
        ax.set_ylabel("RMSD (Å)")
        ax.set_title(column, loc="left", fontsize=10, fontweight="bold")
        ax.set_xlim(0, 500)
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1].set_xlabel("Simulation time (ns)")
    fig.suptitle("UMF-607–KIF20A molecular dynamics simulation", fontsize=13, fontweight="bold")
    fig.savefig(OUT / "Fig5C_UMF607_KIF20A_RMSD.png", dpi=600, bbox_inches="tight")


def contact_panel() -> None:
    contacts = pd.read_csv(RESULTS / "UNK_contact_frequency_filtered.csv")
    top = contacts.nlargest(10, "PercentFrames").sort_values("PercentFrames")
    # Portrait aspect fits the narrow lower-left position in the Figure 5 layout.
    fig, ax = plt.subplots(figsize=(4.25, 7.25), constrained_layout=True)
    bars = ax.barh(top["Residue"], top["PercentFrames"], color="#f26b21", edgecolor="#333333", linewidth=0.5)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Frames with contact (%)", fontsize=10)
    ax.set_ylabel("KIF20A residue", fontsize=10)
    ax.set_title("UMF-607–KIF20A contact frequency", fontsize=11, fontweight="bold", pad=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#d9d9d9", linewidth=0.6)
    ax.set_axisbelow(True)
    for bar, value in zip(bars, top["PercentFrames"]):
        ax.text(value + 1.2, bar.get_y() + bar.get_height() / 2, f"{value:.1f}%", va="center", fontsize=8)
    fig.savefig(OUT / "Fig5D_UMF607_contact_persistence.png", dpi=600, bbox_inches="tight")


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    rmsd_panel()
    contact_panel()
