import importlib.util
import tempfile
import unittest
from unittest import mock
from pathlib import Path
import importlib

import pandas as pd

from ligand_manifest import HAVE_RDKIT, build_ligand_state_metadata, write_ligand_state_manifest
from ligand_naming import parse_ligand_variant


REPO_ROOT = Path(__file__).resolve().parent.parent


def load_script_module(filename: str, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, REPO_ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class LigandNamingTests(unittest.TestCase):
    def test_parse_new_style_variant(self):
        parsed = parse_ligand_variant("obj01_p01_t02_c005")
        self.assertEqual(parsed["LigandBase"], "obj01")
        self.assertEqual(parsed["LigandVariant"], "obj01_p01_t02_c005")
        self.assertEqual(parsed["ProtomerTag"], "p01")
        self.assertEqual(parsed["TautomerTag"], "t02")
        self.assertEqual(parsed["ConformerTag"], "c005")
        self.assertEqual(parsed["StateTag"], "p01_t02")
        self.assertEqual(parsed["ConformerIndex"], 5)

    def test_parse_new_style_with_internal_delimiters(self):
        name = "rxn210____EMOL43216804____EMOL43478895_p01_t02_c005"
        parsed = parse_ligand_variant(name)
        self.assertEqual(parsed["LigandBase"], "rxn210____EMOL43216804____EMOL43478895")
        self.assertEqual(parsed["LigandVariant"], name)

    def test_parse_legacy_pose_variants(self):
        parsed = parse_ligand_variant("UMF-604_pose17")
        self.assertEqual(parsed["LigandBase"], "UMF-604")
        self.assertEqual(parsed["LegacyPoseTag"], "pose17")

        parsed = parse_ligand_variant("UMF-604__pose17")
        self.assertEqual(parsed["LigandBase"], "UMF-604")
        self.assertEqual(parsed["LegacyPoseTag"], "pose17")

    def test_parse_plain_name(self):
        parsed = parse_ligand_variant("plainLigand")
        self.assertEqual(parsed["LigandBase"], "plainLigand")
        self.assertEqual(parsed["LigandVariant"], "plainLigand")
        self.assertEqual(parsed["ConformerTag"], "")


class LigandManifestTests(unittest.TestCase):
    def test_build_manifest_metadata(self):
        if HAVE_RDKIT:
            from rdkit import Chem
            canonical = Chem.MolFromSmiles("CCO")
            state = Chem.MolFromSmiles("CCO")
        else:
            canonical = None
            state = None
        row = build_ligand_state_metadata(
            ligand_base="ethanol",
            ligand_variant="ethanol_p01_t01_c001",
            source_ligand_id="ethanol",
            source_input="input.csv",
            source_input_type="csv",
            canonical_mol=canonical,
            state_mol=state,
            pdbqt_file="ethanol_p01_t01_c001.pdbqt",
            sdf_file="ethanol_p01_t01_c001.sdf",
            tmp_sdf_file="ethanol/ethanol_p01_t01/ethanol_p01_t01_c001.sdf",
            generation_status="ok",
        )
        self.assertEqual(row["LigandBase"], "ethanol")
        self.assertEqual(row["SourceLigandID"], "ethanol")
        self.assertEqual(row["ProtomerTag"], "p01")
        self.assertEqual(row["TautomerTag"], "t01")
        self.assertEqual(row["ConformerTag"], "c001")
        self.assertEqual(row["SourceRecordName"], "")
        if HAVE_RDKIT:
            self.assertEqual(row["CanonicalSMILES"], "CCO")
            self.assertTrue(row["Formula"])
            self.assertNotEqual(row["NumAtoms"], "")


class ParseScoresTests(unittest.TestCase):
    def test_build_score_rows_preserves_variant_metadata(self):
        module = load_script_module("4_ParseScores.py", "parse_scores_module")
        rows = module.build_score_rows(
            "ReceptorA",
            "obj01_p01_t02_c005",
            Path("/tmp/out.pdbqt"),
            [(1, -7.4), (2, -6.9)],
        )
        self.assertEqual(rows[0]["Ligand"], "obj01")
        self.assertEqual(rows[0]["LigandVariant"], "obj01_p01_t02_c005")
        self.assertEqual(rows[0]["ProtomerTag"], "p01")
        self.assertEqual(rows[0]["TautomerTag"], "t02")
        self.assertEqual(rows[0]["ConformerTag"], "c005")
        self.assertEqual(rows[0]["Pose"], 1)

    def test_build_score_rows_merges_manifest_metadata(self):
        module = load_script_module("4_ParseScores.py", "parse_scores_manifest_module")
        manifest_row = {
            "LigandBase": "obj01",
            "LigandVariant": "obj01_p01_t02_c005",
            "CanonicalSMILES": "CCO",
            "Formula": "C2H6O",
            "SourceInput": "input.csv",
            "SourceInputType": "csv",
        }
        rows = module.build_score_rows(
            "ReceptorA",
            "obj01_p01_t02_c005",
            Path("/tmp/out.pdbqt"),
            [(1, -7.4)],
            manifest_row=manifest_row,
        )
        self.assertEqual(rows[0]["CanonicalSMILES"], "CCO")
        self.assertEqual(rows[0]["Formula"], "C2H6O")
        self.assertEqual(rows[0]["SourceInputType"], "csv")


class ConformerInputIdentityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.have_rdkit = importlib.util.find_spec("rdkit") is not None

    def test_filename_preserved_over_internal_sdf_title(self):
        if not self.have_rdkit:
            self.skipTest("rdkit not installed")
        module = load_script_module("1_ConformerGeneration.py", "confgen_module")
        with tempfile.TemporaryDirectory() as tmpdir:
            folder = Path(tmpdir)
            for name in ("DR7", "DR8", "DR9"):
                (folder / f"{name}.sdf").write_text(
                    "obj01\n"
                    "  RDKit          2D\n\n"
                    "  3  2  0  0  0  0            999 V2000\n"
                    "    1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                    "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                    "   -1.2990   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
                    "  1  2  1  0\n"
                    "  2  3  1  0\n"
                    "M  END\n$$$$\n"
                )
            rows = list(module.iter_mols_from_sdf_folder(folder))
            self.assertEqual([row["ligand_id"] for row in rows], ["DR7", "DR8", "DR9"])
            self.assertTrue(all(row["source_record_name"] == "obj01" for row in rows))

    def test_multirecord_sdf_uses_file_stem_and_record_index(self):
        if not self.have_rdkit:
            self.skipTest("rdkit not installed")
        module = load_script_module("1_ConformerGeneration.py", "confgen_module_multi")
        with tempfile.TemporaryDirectory() as tmpdir:
            sdf_path = Path(tmpdir) / "batch.sdf"
            sdf_path.write_text(
                "obj01\n  RDKit          2D\n\n  3  2  0  0  0  0            999 V2000\n"
                "    1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                "   -1.2990   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
                "  1  2  1  0\n  2  3  1  0\nM  END\n$$$$\n"
                "obj01\n  RDKit          2D\n\n  3  2  0  0  0  0            999 V2000\n"
                "    1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                "   -1.2990   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
                "  1  2  1  0\n  2  3  1  0\nM  END\n$$$$\n"
            )
            rows = list(module.iter_mols_from_sdf_file(sdf_path))
            self.assertEqual([row["ligand_id"] for row in rows], ["batch_rec001", "batch_rec002"])
            self.assertTrue(all(row["source_record_name"] == "obj01" for row in rows))

    def test_collision_safe_filename_assignment(self):
        if not self.have_rdkit:
            self.skipTest("rdkit not installed")
        module = load_script_module("1_ConformerGeneration.py", "confgen_module_collision")
        seen = {}
        a = module.ensure_unique_ligand_id(module.ligand_id_from_source_file(Path("DR 7.sdf")), seen, context="DR 7.sdf")
        b = module.ensure_unique_ligand_id(module.ligand_id_from_source_file(Path("DR@7.sdf")), seen, context="DR@7.sdf")
        self.assertEqual(a, "DR_7")
        self.assertEqual(b, "DR_7_002")

    def test_process_one_dir_preserves_variant_in_csv_rows(self):
        module = load_script_module("4_ParseScores.py", "parse_scores_module_process")
        with tempfile.TemporaryDirectory() as tmpdir:
            results_dir = Path(tmpdir) / "Docking_Results_Test"
            lig_dir = Path(tmpdir) / "Ligands"
            lig_dir.mkdir()
            outfile = results_dir / "ReceptorA" / "obj01_p01_t02_c005" / "out.pdbqt"
            outfile.parent.mkdir(parents=True, exist_ok=True)
            outfile.write_text(
                "REMARK VINA RESULT: -7.4 0.0 0.0\n"
                "REMARK VINA RESULT: -6.9 0.0 0.0\n"
            )
            write_ligand_state_manifest(
                lig_dir / "ligand_state_manifest.csv",
                [
                    {
                        **build_ligand_state_metadata(
                            ligand_base="obj01",
                            ligand_variant="obj01_p01_t02_c005",
                            source_ligand_id="obj01",
                            source_input="input.csv",
                            source_input_type="csv",
                            canonical_mol=None,
                            state_mol=None,
                            generation_status="ok",
                            generation_warning="",
                        ),
                        "CanonicalSMILES": "CCO",
                        "StateCanonicalSMILES": "CCO",
                        "Formula": "C2H6O",
                    }
                ],
            )
            rows = module.process_one_dir(results_dir, workers=1, heartbeat=1, fallback_crawl=True)
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["LigandVariant"], "obj01_p01_t02_c005")
            self.assertEqual(rows[0]["LigandBase"], "obj01")
            self.assertEqual(rows[0]["StateTag"], "p01_t02")
            self.assertEqual(rows[0]["CanonicalSMILES"], "CCO")


class ConcatenateScoresTests(unittest.TestCase):
    def test_enrich_row_keeps_internal_delimiters_and_legacy_pose(self):
        module = load_script_module("4C_ConcatenateScores.py", "concat_scores_module")
        row = {
            "Receptor": "ReceptorA",
            "Ligand": "rxn210____EMOL43216804____EMOL43478895",
            "Pose": "1",
            "Binding_Affinity": "-8.1",
            "OutFile": "/tmp/Docking_Results_Test/ReceptorA/rxn210____EMOL43216804____EMOL43478895_p01_t02_c005/out.pdbqt",
        }
        enriched = module.enrich_row(row)
        self.assertEqual(enriched["LigandBase"], "rxn210____EMOL43216804____EMOL43478895")
        self.assertEqual(enriched["LigandVariant"], "rxn210____EMOL43216804____EMOL43478895_p01_t02_c005")
        self.assertEqual(enriched["StateTag"], "p01_t02")

        legacy_row = {
            "Receptor": "ReceptorB",
            "Ligand": "legacyLig_pose38",
            "Pose": "2",
            "Binding_Affinity": "-7.0",
            "OutFile": "/tmp/Docking_Results_Test/ReceptorB/legacyLig_pose38/out.pdbqt",
        }
        enriched_legacy = module.enrich_row(legacy_row)
        self.assertEqual(enriched_legacy["LigandBase"], "legacyLig")
        self.assertEqual(enriched_legacy["LegacyPoseTag"], "pose38")
        self.assertIn("CanonicalSMILES", enriched_legacy)


class DockingRunnerTests(unittest.TestCase):
    def test_build_jobs_preserves_full_ligand_variant(self):
        module = load_script_module("3_Complete_batch_docking.py", "docking_runner_module")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            receptors = root / "Receptors"
            ligands = root / "Ligands"
            results = root / "Docking_Results_Test"
            receptors.mkdir()
            ligands.mkdir()
            (receptors / "ReceptorA.pdbqt").write_text("RECEPTOR\n")
            (ligands / "obj01_p01_t04_c032.pdbqt").write_text("LIGAND\n")
            write_ligand_state_manifest(
                ligands / "ligand_state_manifest.csv",
                [
                    build_ligand_state_metadata(
                        ligand_base="obj01",
                        ligand_variant="obj01_p01_t04_c032",
                        source_ligand_id="obj01",
                        source_input="input.csv",
                        source_input_type="csv",
                        canonical_mol=None,
                        state_mol=None,
                        generation_status="ok",
                        generation_warning="",
                    )
                ],
            )

            jobs = module.build_jobs(
                str(results),
                str(ligands),
                str(receptors),
                [{"PDB_ID": "ReceptorA.pdbqt", "X": "1", "Y": "2", "Z": "3"}],
                5,
                "/usr/bin/vina",
            )
            self.assertEqual(jobs[0]["ligand"], "obj01_p01_t04_c032")
            self.assertTrue(jobs[0]["output_pdbqt"].endswith("ReceptorA/obj01_p01_t04_c032/out.pdbqt"))
            self.assertEqual(jobs[0]["ligand_metadata"]["LigandVariant"], "obj01_p01_t04_c032")

    def test_missing_vina_message_is_clear(self):
        module = load_script_module("3_Complete_batch_docking.py", "docking_runner_module_missing")
        with self.assertRaises(SystemExit) as ctx:
            module.resolve_vina_executable("/definitely/missing/vina")
        self.assertIn("Cannot find AutoDock Vina executable", str(ctx.exception))

    def test_vina_env_path_is_honored(self):
        module = load_script_module("3_Complete_batch_docking.py", "docking_runner_module_env")
        with tempfile.TemporaryDirectory() as tmpdir:
            fake_vina = Path(tmpdir) / "vina"
            fake_vina.write_text("#!/bin/bash\nexit 0\n")
            fake_vina.chmod(0o755)
            with mock.patch.dict("os.environ", {"VINA_EXE": str(fake_vina)}, clear=False):
                resolved = module.resolve_vina_executable(None)
        self.assertEqual(resolved, str(fake_vina))


class GraphSummaryTests(unittest.TestCase):
    def test_build_summary_includes_richer_metrics_and_consensus_rank(self):
        module = load_script_module("7_Graphs.py", "graph_summary_module")
        df = pd.DataFrame(
            [
                {
                    "Receptor": "R1",
                    "Ligand": "LigA",
                    "LigandBase": "LigA",
                    "LigandVariant": "LigA_p01_t01_c001",
                    "Pose": 1,
                    "Binding_Affinity": -10.0,
                    "OutFile": "/tmp/R1/LigA_p01_t01_c001/out.pdbqt",
                    "ProtomerTag": "p01",
                    "TautomerTag": "t01",
                    "ConformerTag": "c001",
                    "StateTag": "p01_t01",
                },
                {
                    "Receptor": "R1",
                    "Ligand": "LigA",
                    "LigandBase": "LigA",
                    "LigandVariant": "LigA_p01_t01_c002",
                    "Pose": 2,
                    "Binding_Affinity": -9.0,
                    "OutFile": "/tmp/R1/LigA_p01_t01_c002/out.pdbqt",
                    "ProtomerTag": "p01",
                    "TautomerTag": "t01",
                    "ConformerTag": "c002",
                    "StateTag": "p01_t01",
                },
                {
                    "Receptor": "R1",
                    "Ligand": "LigA",
                    "LigandBase": "LigA",
                    "LigandVariant": "LigA_p01_t02_c001",
                    "Pose": 3,
                    "Binding_Affinity": -8.0,
                    "OutFile": "/tmp/R1/LigA_p01_t02_c001/out.pdbqt",
                    "ProtomerTag": "p01",
                    "TautomerTag": "t02",
                    "ConformerTag": "c001",
                    "StateTag": "p01_t02",
                },
                {
                    "Receptor": "R1",
                    "Ligand": "LigB",
                    "LigandBase": "LigB",
                    "LigandVariant": "LigB_p01_t01_c001",
                    "Pose": 1,
                    "Binding_Affinity": -7.5,
                    "OutFile": "/tmp/R1/LigB_p01_t01_c001/out.pdbqt",
                    "ProtomerTag": "p01",
                    "TautomerTag": "t01",
                    "ConformerTag": "c001",
                    "StateTag": "p01_t01",
                },
                {
                    "Receptor": "R1",
                    "Ligand": "LigB",
                    "LigandBase": "LigB",
                    "LigandVariant": "LigB_p01_t01_c002",
                    "Pose": 2,
                    "Binding_Affinity": -7.0,
                    "OutFile": "/tmp/R1/LigB_p01_t01_c002/out.pdbqt",
                    "ProtomerTag": "p01",
                    "TautomerTag": "t01",
                    "ConformerTag": "c002",
                    "StateTag": "p01_t01",
                },
            ]
        )

        summary = module.build_summary(df, [-8.0, -9.0])
        lig_a = summary[summary["LigandBase"] == "LigA"].iloc[0]
        lig_b = summary[summary["LigandBase"] == "LigB"].iloc[0]

        self.assertEqual(lig_a["best_affinity"], -10.0)
        self.assertEqual(lig_a["median_affinity"], -9.0)
        self.assertAlmostEqual(lig_a["top3_mean_affinity"], -9.0)
        self.assertAlmostEqual(lig_a["top5_mean_affinity"], -9.0)
        self.assertAlmostEqual(lig_a["q25_affinity"], -9.5)
        self.assertAlmostEqual(lig_a["q75_affinity"], -8.5)
        self.assertAlmostEqual(lig_a["iqr_affinity"], 1.0)
        self.assertAlmostEqual(lig_a["outlier_gap"], 1.0)
        self.assertEqual(lig_a["hit_count_below_-8"], 3)
        self.assertEqual(lig_a["hit_count_below_-9"], 2)
        self.assertAlmostEqual(lig_a["hit_rate_below_-9"], 2 / 3)
        self.assertEqual(lig_a["best_ligand_variant"], "LigA_p01_t01_c001")
        self.assertEqual(lig_a["n_states"], 2)
        self.assertEqual(lig_a["consensus_rank"], 1.0)
        self.assertGreater(lig_b["consensus_rank"], lig_a["consensus_rank"])

    def test_write_state_summary_tables(self):
        module = load_script_module("7_Graphs.py", "graph_module")
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir)
            df = pd.DataFrame(
                [
                    {
                        "Receptor": "R1",
                        "Ligand": "obj01",
                        "LigandBase": "obj01",
                        "LigandVariant": "obj01_p01_t01_c001",
                        "Pose": 1,
                        "Binding_Affinity": -8.0,
                        "OutFile": "/tmp/R1/obj01_p01_t01_c001/out.pdbqt",
                        "ProtomerTag": "p01",
                        "TautomerTag": "t01",
                        "ConformerTag": "c001",
                        "StateTag": "p01_t01",
                        "Formula": "C2H6O",
                        "FormalCharge": 0,
                        "StateCanonicalSMILES": "CCO",
                        "StateIsomericSMILES": "CCO",
                    },
                    {
                        "Receptor": "R1",
                        "Ligand": "obj01",
                        "LigandBase": "obj01",
                        "LigandVariant": "obj01_p01_t02_c005",
                        "Pose": 2,
                        "Binding_Affinity": -7.2,
                        "OutFile": "/tmp/R1/obj01_p01_t02_c005/out.pdbqt",
                        "ProtomerTag": "p01",
                        "TautomerTag": "t02",
                        "ConformerTag": "c005",
                        "StateTag": "p01_t02",
                        "Formula": "C2H6O",
                        "FormalCharge": 0,
                        "StateCanonicalSMILES": "CCO",
                        "StateIsomericSMILES": "CCO",
                    },
                ]
            )
            written = module.write_state_summary_tables(df, outdir, "scores")
            self.assertTrue(any(path.name.endswith("ligand_base_summary.csv") for path in written.values()))
            tautomer_csv = outdir / "scores_tautomer_summary.csv"
            self.assertTrue(tautomer_csv.exists())
            tautomer_df = pd.read_csv(tautomer_csv)
            self.assertIn("best_ligand_variant", tautomer_df.columns)
            self.assertIn("StateCanonicalSMILES", tautomer_df.columns)
            self.assertIn("top5_mean_affinity", tautomer_df.columns)
            self.assertIn("hit_rate_below_-9", tautomer_df.columns)

    def test_generate_outputs_writes_manifest_and_summary_tables_without_plot_dependencies(self):
        module = load_script_module("7_Graphs.py", "graph_outputs_module")
        with tempfile.TemporaryDirectory() as tmpdir:
            csv_path = Path(tmpdir) / "scores.csv"
            outdir = Path(tmpdir) / "Plots"
            df = pd.DataFrame(
                [
                    {
                        "Receptor": "R1",
                        "Ligand": "LigA",
                        "LigandBase": "LigA",
                        "LigandVariant": "LigA_p01_t01_c001",
                        "Pose": 1,
                        "Binding_Affinity": -8.4,
                        "OutFile": "/tmp/R1/LigA_p01_t01_c001/out.pdbqt",
                        "ProtomerTag": "p01",
                        "TautomerTag": "t01",
                        "ConformerTag": "c001",
                        "StateTag": "p01_t01",
                    }
                ]
            )
            df.to_csv(csv_path, index=False)

            result = module.generate_outputs(
                df=df,
                csv_path=csv_path,
                outdir=outdir,
                top_n=4,
                bins=20,
                min_scores=1,
                hit_thresholds=[-8.0, -9.0],
                ranking_metric="consensus",
                write_state_summaries=True,
                write_ranking_plots=True,
                write_hit_rate_plots=True,
                write_state_heatmaps=True,
                write_state_plots=False,
                state_heatmap_metric="best",
                write_pdf_report_flag=False,
                pdf_top_ligands=10,
                skip_per_ligand=True,
                skip_per_receptor=True,
            )

            self.assertTrue((outdir / "summary_tables" / "scores_ligand_summary.csv").exists())
            self.assertTrue((outdir / "summary_tables" / "scores_variant_summary.csv").exists())
            self.assertTrue((outdir / "graph_report_manifest.json").exists())
            self.assertIn("summary_csvs", result["manifest"])
            self.assertIn("plot_files", result["manifest"])
            self.assertEqual(result["manifest"]["pdf_report"], "")


if __name__ == "__main__":
    unittest.main()
