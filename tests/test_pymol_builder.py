import csv
import io
import importlib
import importlib.util
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parent.parent


def load_script_module(filename: str, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, REPO_ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class PymolBuilderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.have_rdkit = importlib.util.find_spec("rdkit") is not None

    def test_infer_tmp_sdf_relative_path(self):
        module = load_script_module("5C_BuildPymolSesh.py", "pymol_builder_tmp_path")
        rel = module.infer_tmp_sdf_relative_path("DR7_p01_t02_c005")
        self.assertEqual(rel, Path("DR7") / "DR7_p01_t02" / "DR7_p01_t02_c005.sdf")

    def test_discover_score_csvs_filters_tmp_and_non_score_csvs(self):
        module = load_script_module("5C_BuildPymolSesh.py", "pymol_builder_score_csv_discovery")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            combined = root / "ALL_Docking_Results_with_provenance_TOP10000_2026-07-29_05-26-47.csv"
            combined.write_text("Receptor,Ligand,Pose,Binding_Affinity,OutFile\n", encoding="utf-8")
            per_dir = root / "Docking_Results_demo_2026-07-27_14-28-28_vina_docking_scores_sorted.csv"
            per_dir.write_text("Receptor,Ligand,Pose,Binding_Affinity,OutFile\n", encoding="utf-8")
            tmp = root / "ALL_Docking_Results_with_provenance_TOPALL_2026-07-29_05-26-01_UNSORTED.tmp.csv"
            tmp.write_text("Receptor,Ligand,Pose,Binding_Affinity,OutFile\n", encoding="utf-8")
            centers = root / "vina_centers.csv"
            centers.write_text("PDB_ID,X,Y,Z,SIZE\nrec.pdbqt,1,2,3,20\n", encoding="utf-8")

            self.assertEqual(module.discover_score_csvs(root), [combined, per_dir])

    def test_print_selection_summary_uses_rank_order_with_scores(self):
        module = load_script_module("5C_BuildPymolSesh.py", "pymol_builder_ranked_summary")
        grouped = {
            "9ba9": [
                {"LigandBase": "CONIVAPTAN", "Binding_Affinity": "-9.344"},
                {"LigandBase": "DIHYDROERGOTAMINE_MESYLATE", "Binding_Affinity": "-9.377"},
                {"LigandBase": "CONIVAPTAN", "Binding_Affinity": "-9.100"},
            ]
        }
        stream = io.StringIO()
        with mock.patch("sys.stdout", stream):
            module.print_selection_summary(grouped, "per_ligand", 15)

        output = stream.getvalue()
        first = output.index("[1] DIHYDROERGOTAMINE_MESYLATE")
        second = output.index("[2] CONIVAPTAN")
        self.assertLess(first, second)
        self.assertIn("best=-9.377  poses=1", output)
        self.assertIn("best=-9.344  poses=2", output)

    def test_manifest_lookup_returns_exact_variant_sdf(self):
        module = load_script_module("5C_BuildPymolSesh.py", "pymol_builder_manifest_lookup")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_dir = root / "Ligands_Ligands_PDBQT_8Poses_20260611_1718"
            manifest_dir.mkdir()
            tmp_root = root / "Ligands_TMP_SDF_8Poses_20260611_1718" / "DR7" / "DR7_p01_t02"
            tmp_root.mkdir(parents=True)
            sdf_path = tmp_root / "DR7_p01_t02_c005.sdf"
            sdf_path.write_text("fake sdf\n")
            with (manifest_dir / "ligand_state_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["LigandVariant", "TmpSDFFile", "SDFFile"])
                writer.writeheader()
                writer.writerow(
                    {
                        "LigandVariant": "DR7_p01_t02_c005",
                        "TmpSDFFile": "DR7/DR7_p01_t02/DR7_p01_t02_c005.sdf",
                        "SDFFile": "DR7_p01_t02_c005.sdf",
                    }
                )

            resolver = module.ReferenceResolver(root)
            resolution = resolver.resolve({"LigandVariant": "DR7_p01_t02_c005", "LigandBase": "DR7"})
            self.assertEqual(resolution["reference_sdf"], sdf_path.resolve())
            self.assertEqual(resolution["reference_lookup_method"], "manifest:TmpSDFFile")
            self.assertEqual(resolution["reference_lookup_status"], "exact_manifest_match")
            self.assertEqual(resolution["reference_lookup_warning"], "")
            self.assertEqual(resolution["manifest_row"]["LigandVariant"], "DR7_p01_t02_c005")
            self.assertTrue(resolution["is_exact_variant_match"])

    def test_reference_resolver_indexes_only_direct_ligand_manifests(self):
        module = load_script_module("5C_BuildPymolSesh.py", "pymol_builder_manifest_fast_paths")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            ligand_dir = root / "chembl_part_001_Ligands_PDBQT_64Poses_20260714_0754"
            ligand_dir.mkdir()
            with (ligand_dir / "ligand_state_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["LigandVariant", "TmpSDFFile"])
                writer.writeheader()
                writer.writerow({"LigandVariant": "Good_p01_t01_c001", "TmpSDFFile": "Good/Good_p01_t01/Good_p01_t01_c001.sdf"})

            deep_results = root / "Docking_Results_demo" / "recA" / "Bad_p01_t01_c001"
            deep_results.mkdir(parents=True)
            with (deep_results / "ligand_state_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["LigandVariant", "TmpSDFFile"])
                writer.writeheader()
                writer.writerow({"LigandVariant": "Bad_p01_t01_c001", "TmpSDFFile": "Bad.sdf"})

            resolver = module.ReferenceResolver(root)

            self.assertIn("Good_p01_t01_c001", resolver.manifest_rows)
            self.assertNotIn("Bad_p01_t01_c001", resolver.manifest_rows)

    def test_previous_corrected_outputs_are_excluded_from_reference_lookup(self):
        module = load_script_module("5C_BuildPymolSesh.py", "pymol_builder_exclude_previous_outputs")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            old_outputs = root / "11_AA_Ligands"
            old_outputs.mkdir()
            old_ref = old_outputs / "DR7_p01_t02_c005_reference.sdf"
            old_ref.write_text("fake sdf\n")
            resolver = module.ReferenceResolver(root)
            resolution = resolver.resolve({"LigandVariant": "DR7_p01_t02_c005", "LigandBase": "DR7"})
            self.assertEqual(resolution["reference_lookup_status"], "missing_exact_reference")
            self.assertFalse(resolution["is_exact_variant_match"])
            self.assertEqual(resolution["reference_sdf"], None)

    def test_base_source_fallback_is_marked_as_non_exact(self):
        module = load_script_module("5C_BuildPymolSesh.py", "pymol_builder_fallback_reference")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            ligands = root / "Ligands"
            ligands.mkdir()
            base_sdf = ligands / "DR7.sdf"
            base_sdf.write_text("fake sdf\n")
            resolver = module.ReferenceResolver(root)
            resolution = resolver.resolve(
                {
                    "LigandVariant": "DR7_p01_t02_c005",
                    "LigandBase": "DR7",
                    "SourceInput": "Ligands/DR7.sdf",
                }
            )
            self.assertEqual(resolution["reference_sdf"], base_sdf.resolve())
            self.assertEqual(resolution["reference_lookup_status"], "fallback_reference")
            self.assertTrue(resolution["is_fallback"])
            self.assertFalse(resolution["is_exact_variant_match"])

    def test_audit_reference_resolution_mode_writes_csv_without_runtime_tools(self):
        module = load_script_module("5C_BuildPymolSesh.py", "pymol_builder_reference_audit")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            csv_path = root / "scores.csv"
            receptors = root / "Receptors"
            receptors.mkdir()
            tmp_root = root / "Ligands_TMP_SDF_8Poses_demo" / "DR7" / "DR7_p01_t02"
            tmp_root.mkdir(parents=True)
            sdf_path = tmp_root / "DR7_p01_t02_c005.sdf"
            sdf_path.write_text("fake sdf\n")
            pdqt_root = root / "Ligands_Ligands_PDBQT_demo"
            pdqt_root.mkdir()
            with (pdqt_root / "ligand_state_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["LigandVariant", "TmpSDFFile"])
                writer.writeheader()
                writer.writerow(
                    {
                        "LigandVariant": "DR7_p01_t02_c005",
                        "TmpSDFFile": "DR7/DR7_p01_t02/DR7_p01_t02_c005.sdf",
                    }
                )
            with csv_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["Receptor", "Ligand", "LigandVariant", "LigandBase", "Pose", "Binding_Affinity"])
                writer.writeheader()
                writer.writerow(
                    {
                        "Receptor": "R1",
                        "Ligand": "DR7",
                        "LigandVariant": "DR7_p01_t02_c005",
                        "LigandBase": "DR7",
                        "Pose": 1,
                        "Binding_Affinity": -8.4,
                    }
                )
            old_cwd = os.getcwd()
            try:
                os.chdir(root)
                with mock.patch("sys.argv", [
                    "5C_BuildPymolSesh.py",
                    "--csv", str(csv_path),
                    "--mode", "per_ligand",
                    "--top", "5",
                    "--receptor-roots", str(receptors),
                    "--audit-reference-resolution",
                    "--non-interactive",
                ]), mock.patch.object(module, "obabel_path", side_effect=AssertionError("obabel should not be called")), mock.patch.object(module, "launch_pymol", side_effect=AssertionError("pymol should not be launched")):
                    module.main()
            finally:
                os.chdir(old_cwd)
            audit_files = sorted((root / "11_AA_PDBOUTPUT").glob("reference_resolution_audit_*.csv"))
            self.assertTrue(audit_files)
            with audit_files[-1].open() as handle:
                audit_rows = list(csv.DictReader(handle))
            self.assertEqual(audit_rows[0]["resolved_sdf_stem"], "DR7_p01_t02_c005")
            self.assertEqual(audit_rows[0]["is_exact_variant_match"], "True")

    def test_extract_selected_pose_pdbqt_uses_requested_model(self):
        module = load_script_module("5C_BuildPymolSesh.py", "pymol_builder_pose_extract")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            src = root / "out.pdbqt"
            src.write_text(
                "MODEL 1\nREMARK VINA RESULT: -7.0 0 0\nATOM A\nENDMDL\n"
                "MODEL 2\nREMARK VINA RESULT: -8.0 0 0\nATOM B\nENDMDL\n"
            )
            dest = root / "pose2.pdbqt"
            result, warning = module.extract_selected_pose_pdbqt(src, 2, dest)
            self.assertEqual(result, dest)
            text = dest.read_text()
            self.assertIn("MODEL 2", text)
            self.assertIn("ATOM B", text)
            self.assertNotIn("ATOM A", text)
            self.assertEqual(warning, "")

    def test_group_names_do_not_collide_with_receptor_object(self):
        module = load_script_module("5C_BuildPymolSesh.py", "pymol_builder_groups")
        names = module.build_pymol_names("3eky", "DR7", "DR7_p01_t02_c005", 1)
        self.assertEqual(names["receptor_obj"], "obj_3eky")
        self.assertEqual(names["receptor_group"], "grp_3eky")
        self.assertEqual(names["ligand_group"], "grp_3eky_DR7")
        self.assertNotEqual(names["receptor_obj"], names["receptor_group"])

    def _build_test_mol(self):
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.AddHs(Chem.MolFromSmiles("CC(O)N"))
        AllChem.EmbedMolecule(mol, randomSeed=0xF00D)
        return mol

    def test_strip_all_hydrogens_preserves_heavy_atoms(self):
        if not self.have_rdkit:
            self.skipTest("rdkit not installed")
        module = load_script_module("5C_BuildPymolSesh.py", "pymol_builder_strip_h")
        mol = self._build_test_mol()
        stripped = module.strip_all_hydrogens(mol)
        counts = module.atom_count_summary(stripped)
        self.assertEqual(counts["hydrogen"], 0)
        self.assertTrue(module.heavy_atoms_preserved(mol, stripped))

    def test_add_nonpolar_hydrogens_only_keeps_only_carbon_bound_h(self):
        if not self.have_rdkit:
            self.skipTest("rdkit not installed")
        module = load_script_module("5C_BuildPymolSesh.py", "pymol_builder_nonpolar_h")
        mol = self._build_test_mol()
        processed = module.add_nonpolar_hydrogens_only(mol)
        self.assertTrue(module.heavy_atoms_preserved(mol, processed))
        for atom in processed.GetAtoms():
            if atom.GetAtomicNum() != 1:
                continue
            neighbor = atom.GetNeighbors()[0]
            self.assertEqual(neighbor.GetAtomicNum(), 6)

    def test_apply_hydrogen_policy_modes(self):
        if not self.have_rdkit:
            self.skipTest("rdkit not installed")
        module = load_script_module("5C_BuildPymolSesh.py", "pymol_builder_h_modes")
        mol = self._build_test_mol()
        all_counts = module.atom_count_summary(mol)

        none_mol, none_msg = module.apply_hydrogen_policy(mol, "none")
        self.assertEqual(module.atom_count_summary(none_mol)["hydrogen"], 0)
        self.assertIn("Removed all hydrogens", none_msg)
        self.assertTrue(module.heavy_atoms_preserved(mol, none_mol))

        nonpolar_mol, nonpolar_msg = module.apply_hydrogen_policy(mol, "nonpolar")
        self.assertLess(module.atom_count_summary(nonpolar_mol)["hydrogen"], all_counts["hydrogen"])
        self.assertIn("carbon-bound hydrogens", nonpolar_msg)
        self.assertTrue(module.heavy_atoms_preserved(mol, nonpolar_mol))

        all_mol, all_msg = module.apply_hydrogen_policy(mol, "all")
        self.assertEqual(module.atom_count_summary(all_mol)["hydrogen"], all_counts["hydrogen"])
        self.assertIn("Preserved all hydrogens", all_msg)
        self.assertTrue(module.heavy_atoms_preserved(mol, all_mol))

    def test_heavy_atom_count_mismatch_is_rejected(self):
        if not self.have_rdkit:
            self.skipTest("rdkit not installed")
        module = load_script_module("5C_BuildPymolSesh.py", "pymol_builder_mapping_mismatch")
        from rdkit import Chem
        from rdkit.Chem import AllChem

        ref = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        dock = Chem.AddHs(Chem.MolFromSmiles("CCCO"))
        AllChem.EmbedMolecule(ref, randomSeed=0xBEEF)
        AllChem.EmbedMolecule(dock, randomSeed=0xCAFE)
        result = module.transfer_docked_coordinates(ref, dock, allow_partial_mcs=False)
        self.assertEqual(result["status"], "heavy_atom_count_mismatch")
        self.assertFalse(result["complete"])

    def test_exact_heavy_atom_order_mapping_is_complete(self):
        if not self.have_rdkit:
            self.skipTest("rdkit not installed")
        module = load_script_module("5C_BuildPymolSesh.py", "pymol_builder_mapping_order")
        mol = self._build_test_mol()
        result = module.transfer_docked_coordinates(mol, mol, allow_partial_mcs=False)
        self.assertEqual(result["status"], "complete")
        self.assertEqual(result["method"], "heavy_atom_order")
        self.assertTrue(result["complete"])

    def test_atom_map_audit_csv_is_written(self):
        if not self.have_rdkit:
            self.skipTest("rdkit not installed")
        module = load_script_module("5C_BuildPymolSesh.py", "pymol_builder_atom_map_audit")
        from rdkit import Chem
        from rdkit.Chem import AllChem

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            reference_sdf = root / "ref.sdf"
            selected_pdbqt = root / "pose.pdbqt"
            selected_pdbqt.write_text("MODEL 1\nENDMDL\n")
            mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
            AllChem.EmbedMolecule(mol, randomSeed=0x1234)
            result = module.transfer_docked_coordinates(mol, mol, allow_partial_mcs=False)
            out_csv = root / "atom_map.csv"
            module.write_atom_map_audit_csv(
                out_csv,
                "LigA_p01_t01_c001",
                reference_sdf,
                selected_pdbqt,
                result["method"],
                mol,
                mol,
                result["mol"],
                result["atom_map"],
            )
            self.assertTrue(out_csv.exists())
            rows = list(csv.DictReader(out_csv.open()))
            self.assertTrue(rows)
            self.assertEqual(rows[0]["ligand_variant"], "LigA_p01_t01_c001")


if __name__ == "__main__":
    unittest.main()
