from pymol import cmd
import csv
import os

# === CONFIGURATION ===
receptor_dir = "./Receptors"
output_csv = "vina_centers.csv"
current_receptor = None
remaining_receptors = []

# === Ensure CSV has correct headers ===
if not os.path.exists(output_csv) or os.stat(output_csv).st_size == 0:
    with open(output_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["PDB_ID", "X", "Y", "Z", "SIZE"])

def center_this():
    global current_receptor
    global remaining_receptors
    if not current_receptor:
        print("⚠️ No receptor currently loaded.")
        return
    center_from_selection(current_receptor)
    load_next_receptor()

def center_from_selection(receptor_name):
    model = cmd.get_model("pocket")
    if len(model.atom) == 0:
        print("❌ No atoms in selection 'pocket'")
        return

    x, y, z = 0, 0, 0
    for atom in model.atom:
        x += atom.coord[0]
        y += atom.coord[1]
        z += atom.coord[2]
    n = len(model.atom)
    center = (x/n, y/n, z/n)

    with open(output_csv, "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([f"{receptor_name}.pdbqt", round(center[0], 3), round(center[1], 3), round(center[2], 3), 20])

    print(f"✅ Center for {receptor_name} written: {center}")

def load_next_receptor():
    global current_receptor
    global remaining_receptors

    if not remaining_receptors:
        print("🎉 All receptors processed!")
        return

    fname = remaining_receptors.pop(0)
    filepath = os.path.join(receptor_dir, fname)
    object_name = fname.replace(".pdbqt", "")

    cmd.reinitialize()
    cmd.load(filepath, object_name)
    current_receptor = object_name
    print(f"📌 Loaded {fname}. Please define a selection named 'pocket' and then run: center_this()")

# === Prepare list of unprocessed receptors ===
receptors = [f for f in os.listdir(receptor_dir) if f.endswith(".pdbqt")]
processed = set()

if os.path.exists(output_csv):
    with open(output_csv, "r") as f:
        reader = csv.DictReader(f)
        processed = {row["PDB_ID"] for row in reader if "PDB_ID" in row}

remaining_receptors = [f for f in receptors if f not in processed]

# === Load the first receptor ===
load_next_receptor()
