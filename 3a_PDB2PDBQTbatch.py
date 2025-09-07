import os
import subprocess

# Path to your local AutoDockTools prepare script
prepare_script = os.path.join("AutoDockTools_py3", "AutoDockTools", "Utilities24", "prepare_receptor4.py")
if not os.path.exists(prepare_script):
    raise FileNotFoundError("❌ Could not find prepare_receptor4.py at expected location.")

# Prompt for input folder
print("\n📁 Available folders in current directory:")
folders = [f for f in os.listdir() if os.path.isdir(f)]
for i, folder in enumerate(folders, 1):
    print(f"{i}. {folder}")

choice = input("🔍 Enter the number of the folder containing receptor .pdb files: ")
try:
    folder_index = int(choice) - 1
    selected_folder = folders[folder_index]
except (ValueError, IndexError):
    raise ValueError("❌ Invalid folder selection.")

# Prompt for residue names to remove
residue_input = input("🧽 Enter 3-letter residue names to remove (comma-separated, e.g. SOL,ADP,HOH,NA): ")
residues_to_remove = set(res.strip().upper() for res in residue_input.split(",") if res.strip())

print(f"\n🧼 Will remove residues: {', '.join(residues_to_remove)}")

# Output folder for cleaned + converted files
output_dir = selected_folder + "_PDBQT_Converted"
os.makedirs(output_dir, exist_ok=True)
# Temporary folder for cleaned PDBs
temp_cleaned_dir = ".temp_cleaned_pdbs"
os.makedirs(temp_cleaned_dir, exist_ok=True)

# Process each file
for fname in sorted(os.listdir(selected_folder)):
    if not fname.endswith(".pdb"):
        continue

    input_path = os.path.join(selected_folder, fname)
    cleaned_pdb = os.path.join(temp_cleaned_dir, fname)

    # Clean file by removing residues
    with open(input_path, "r") as fin, open(cleaned_pdb, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                resname = line[17:20].strip().upper()
                if resname in residues_to_remove:
                    continue
            fout.write(line)

    # Convert cleaned PDB to PDBQT
    output_name = fname.replace(".pdb", "") + ".converted.pdbqt"
    output_path = os.path.join(output_dir, output_name)

    cmd = f'python "{prepare_script}" -r "{cleaned_pdb}" -o "{output_path}" -A checkhydrogens'
    print(f"\n⚙️ Converting {fname} → {output_name}")
    try:
        subprocess.run(cmd, shell=True, check=True)
        print(f"✅ Saved: {output_path}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed on {fname}: {e}")

# Cleanup temporary files
import shutil
shutil.rmtree(temp_cleaned_dir)

print(f"\n🎉 Done. All cleaned and converted receptors saved to: {output_dir}")
