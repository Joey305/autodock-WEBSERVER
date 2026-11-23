import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

from ligand_manifest import build_ligand_state_metadata, write_ligand_state_manifest


ROOT = Path(__file__).resolve().parents[1]


def load_script(name: str, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, ROOT / name)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


class CompactResultsTests(unittest.TestCase):
    def test_scratch_precedence(self):
        runner = load_script("3_Complete_batch_docking.py", "compact_runner_precedence")
        self.assertEqual(runner.resolve_compact_scratch_root("/cli", "/results")[0], "/cli")
        old = __import__("os").environ.get("VINA_SCRATCH_ROOT")
        __import__("os").environ["VINA_SCRATCH_ROOT"] = "/environment"
        self.assertEqual(runner.resolve_compact_scratch_root(None, "/results")[0], "/environment")
        if old is None:
            del __import__("os").environ["VINA_SCRATCH_ROOT"]
        else:
            __import__("os").environ["VINA_SCRATCH_ROOT"] = old
        self.assertEqual(runner.resolve_compact_scratch_root(None, "/results")[0], "/results")

    def test_compact_runner_resumes_without_pose_tree_and_parser_reads_it(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            receptors = root / "Receptors"; receptors.mkdir()
            ligands = root / "Ligands"; ligands.mkdir()
            (receptors / "recA.pdbqt").write_text("RECEPTOR\n", encoding="utf-8")
            (ligands / "ligA.pdbqt").write_text("LIGAND\n", encoding="utf-8")
            (ligands / "ligB.pdbqt").write_text("LIGAND\n", encoding="utf-8")
            centers = root / "centers.csv"
            centers.write_text("PDB_ID,X,Y,Z\nrecA,1,2,3\n", encoding="utf-8")
            vina = root / "fake_vina.py"
            vina.write_text(
                "#!/usr/bin/env python3\nimport pathlib,sys\n"
                "if '--version' in sys.argv: print('Fake Vina 1.0'); raise SystemExit(0)\n"
                "cfg=pathlib.Path(sys.argv[sys.argv.index('--config')+1]).read_text()\n"
                "out=[x.split(' = ',1)[1] for x in cfg.splitlines() if x.startswith('out = ')][0]\n"
                "pathlib.Path(out).write_text('REMARK VINA RESULT: -7.5 0.0 0.0\\nREMARK VINA RESULT: -7.0 1.0 2.0\\n')\n",
                encoding="utf-8",
            )
            vina.chmod(0o755)
            results = root / "results"
            command = [sys.executable, str(ROOT / "3_Complete_batch_docking.py"), "--receptors", str(receptors),
                       "--ligands", str(ligands), "--centers_csv", str(centers), "--poses", "2",
                       "--vina-exe", str(vina), "--results-dir", str(results), "--compact-results",
                       "--workers", "1", "--max-pending", "1", "--seed", "42"]
            subprocess.run(command, cwd=root, check=True, capture_output=True, text=True)
            with (results / "compact_results.csv").open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 4)
            self.assertFalse((results / "recA").exists())
            self.assertFalse(list(results.glob("vina_job_*")))
            completion = json.loads((results / "compact_completion.json").read_text(encoding="utf-8"))
            self.assertTrue(completion["poses_retained"] is False)
            progress = json.loads((results / "vina_progress.json").read_text(encoding="utf-8"))
            self.assertEqual((progress["total"], progress["attempted"], progress["successful"], progress["failed"]), (2, 2, 2, 0))
            self.assertIn("updated_at", progress)
            subprocess.run(command, cwd=root, check=True, capture_output=True, text=True)
            with (results / "compact_results.csv").open(newline="", encoding="utf-8") as handle:
                self.assertEqual(len(list(csv.DictReader(handle))), 4)
            subprocess.run(command + ["--compact-retain-top-poses", "1"], cwd=root, check=True, capture_output=True, text=True)
            retained = json.loads((results / "compact_retained_poses.json").read_text(encoding="utf-8"))
            self.assertEqual(len(retained["retained"]), 1)
            self.assertTrue((results / retained["retained"][0]["path"]).is_file())
            with (results / "compact_results.csv").open(newline="", encoding="utf-8") as handle:
                self.assertEqual(len(list(csv.DictReader(handle))), 4)
            parser = load_script("4_ParseScores.py", "compact_score_parser")
            parsed = parser.process_compact_results(results)
            self.assertEqual(len(parsed), 4)
            self.assertEqual(parsed[0]["OutFile"], "")

    def test_retained_poses_group_parent_states_rank_ties_cleanup_and_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            receptors, ligands = root / "Receptors", root / "Ligands"
            receptors.mkdir(); ligands.mkdir()
            (receptors / "recA.pdbqt").write_text("EXACT RECEPTOR\n", encoding="utf-8")
            states = {
                "parentA_prot1_taut1_conf1": ("Parent A", "p1", "t1", "c1", -9.0),
                "parentA_prot2_taut2_conf2": ("Parent A", "p2", "t2", "c2", -10.0),
                "parentB_prot1_taut1_conf1": ("Parent B", "p1", "t1", "c1", -10.0),
                "parentC_prot1_taut1_conf1": ("Parent C", "p1", "t1", "c1", -8.0),
            }
            for state in states:
                (ligands / f"{state}.pdbqt").write_text("LIGAND\n", encoding="utf-8")
            write_ligand_state_manifest(ligands / "ligand_state_manifest.csv", [
                build_ligand_state_metadata(
                    parent, state, parent, "library.csv", "csv",
                    pdbqt_file=f"{state}.pdbqt",
                ) | {"ProtomerTag": protomer, "TautomerTag": tautomer, "ConformerTag": conformer}
                for state, (parent, protomer, tautomer, conformer, _) in states.items()
            ])
            (root / "centers.csv").write_text("PDB_ID,X,Y,Z\nrecA,1,2,3\n", encoding="utf-8")
            vina = root / "fake_vina.py"
            affinities = {state: score for state, (*_, score) in states.items()} | {"parentD_prot1_taut1_conf1": -11.0}
            vina.write_text(
                "#!/usr/bin/env python3\nimport pathlib,sys\n"
                f"scores={affinities!r}\n"
                "if '--version' in sys.argv: print('Fake Vina'); raise SystemExit(0)\n"
                "cfg=pathlib.Path(sys.argv[sys.argv.index('--config')+1]).read_text()\n"
                "lig=pathlib.Path([x.split(' = ',1)[1] for x in cfg.splitlines() if x.startswith('ligand = ')][0]).stem\n"
                "out=[x.split(' = ',1)[1] for x in cfg.splitlines() if x.startswith('out = ')][0]\n"
                "pathlib.Path(out).write_text(f'MODEL 1\\nREMARK VINA RESULT: {scores[lig]} 0.0 0.0\\nENDMDL\\nMODEL 2\\nREMARK VINA RESULT: -1.0 0.0 0.0\\nENDMDL\\n')\n",
                encoding="utf-8",
            )
            vina.chmod(0o755)
            results = root / "results"
            command = [sys.executable, str(ROOT / "3_Complete_batch_docking.py"), "--receptors", str(receptors),
                       "--ligands", str(ligands), "--centers_csv", str(root / "centers.csv"), "--poses", "2",
                       "--vina-exe", str(vina), "--results-dir", str(results), "--compact-results",
                       "--compact-retain-top-poses", "2", "--compact-retain-poses-per-parent", "2",
                       "--workers", "1", "--max-pending", "1", "--seed", "11"]
            subprocess.run(command, cwd=root, check=True, capture_output=True, text=True)
            manifest_path = results / "compact_retained_poses.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            retained = manifest["retained"]
            self.assertEqual([entry["parent_ligand"] for entry in retained], ["Parent A", "Parent A", "Parent B", "Parent B"])
            self.assertEqual(retained[0]["state_specific_ligand"], "parentA_prot2_taut2_conf2")
            self.assertEqual(retained[0]["state_specific_ligand_filename"], "parentA_prot2_taut2_conf2.pdbqt")
            self.assertEqual(retained[0]["vina_pose_rank"], 1)
            self.assertEqual(retained[0]["protomer"], "p2")
            self.assertEqual(len({entry["parent_ligand"] for entry in retained}), 2)
            self.assertEqual([entry["within_parent_rank"] for entry in retained[:2]], [1, 2])
            self.assertEqual([entry["parent_rank"] for entry in retained], [1, 1, 2, 2])
            self.assertIn("state_specific_ligand", retained[0]["tie_breakers"])
            self.assertEqual((results / retained[0]["receptor_path"]).read_text(encoding="utf-8"), "EXACT RECEPTOR\n")
            for entry in retained:
                pose = (results / entry["path"]).read_text(encoding="utf-8")
                self.assertIn(f"MODEL {entry['vina_pose_rank']}", pose)
                self.assertEqual(pose.count("MODEL "), 1)

            # A later/retry shard adds a better parent; its snapshot merges
            # with the earlier winners and removes the displaced Parent B file.
            state = "parentD_prot1_taut1_conf1"
            (ligands / f"{state}.pdbqt").write_text("LIGAND\n", encoding="utf-8")
            write_ligand_state_manifest(ligands / "ligand_state_manifest.csv", [
                build_ligand_state_metadata(parent, name, parent, "library.csv", "csv", pdbqt_file=f"{name}.pdbqt")
                for name, (parent, *_rest) in states.items()
            ] + [build_ligand_state_metadata("Parent D", state, "Parent D", "library.csv", "csv", pdbqt_file=f"{state}.pdbqt")])
            displaced = results / next(entry["path"] for entry in retained if entry["parent_ligand"] == "Parent B")
            subprocess.run(command, cwd=root, check=True, capture_output=True, text=True)
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(
                [entry["parent_ligand"] for entry in manifest["retained"]],
                ["Parent D", "Parent D", "Parent A", "Parent A"],
            )
            self.assertFalse(displaced.exists())
            self.assertTrue(all((results / entry["path"]).is_file() for entry in manifest["retained"]))

    def test_compact_runner_stages_selected_inputs_and_cleans_workspace(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            receptors = root / "Receptors"; receptors.mkdir()
            ligands = root / "Ligands"; ligands.mkdir()
            scratch = root / "staging"; scratch.mkdir()
            (receptors / "recA.pdbqt").write_text("RECEPTOR\n", encoding="utf-8")
            (ligands / "ligA.pdbqt").write_text("LIGAND\n", encoding="utf-8")
            (root / "centers.csv").write_text("PDB_ID,X,Y,Z\nrecA,1,2,3\n", encoding="utf-8")
            vina = root / "fake_vina.py"
            vina.write_text(
                "#!/usr/bin/env python3\nimport pathlib,sys\n"
                "if '--version' in sys.argv: print('Fake Vina 1.0'); raise SystemExit(0)\n"
                "cfg=pathlib.Path(sys.argv[sys.argv.index('--config')+1]).read_text()\n"
                "pathlib.Path(__file__).with_name('seen_config.txt').write_text(cfg)\n"
                "out=[x.split(' = ',1)[1] for x in cfg.splitlines() if x.startswith('out = ')][0]\n"
                "pathlib.Path(out).write_text('REMARK VINA RESULT: -7.5 0.0 0.0\\n')\n",
                encoding="utf-8",
            )
            vina.chmod(0o755)
            subprocess.run([
                sys.executable, str(ROOT / "3_Complete_batch_docking.py"),
                "--receptors", str(receptors), "--ligands", str(ligands),
                "--centers_csv", str(root / "centers.csv"), "--poses", "1",
                "--vina-exe", str(vina), "--results-dir", str(root / "results"),
                "--compact-results", "--compact-scratch-root", str(scratch),
                "--stage-inputs-to-scratch", "--workers", "1", "--max-pending", "1",
            ], cwd=root, check=True, capture_output=True, text=True)
            config = (root / "seen_config.txt").read_text(encoding="utf-8")
            self.assertIn(str(scratch / "vina_inputs_"), config)
            self.assertNotIn(str(ligands / "ligA.pdbqt"), config)
            self.assertFalse(list(scratch.glob("vina_inputs_*")))

    def test_interrupted_compact_run_resumes_without_duplicate_scores(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            receptors = root / "Receptors"; receptors.mkdir()
            ligands = root / "Ligands"; ligands.mkdir()
            (receptors / "recA.pdbqt").write_text("RECEPTOR\n", encoding="utf-8")
            (ligands / "ligA.pdbqt").write_text("LIGAND\n", encoding="utf-8")
            (ligands / "ligB.pdbqt").write_text("LIGAND\n", encoding="utf-8")
            (root / "centers.csv").write_text("PDB_ID,X,Y,Z\nrecA,1,2,3\n", encoding="utf-8")
            vina = root / "slow_fake_vina.py"
            vina.write_text(
                "#!/usr/bin/env python3\nimport pathlib,sys,time\n"
                "if '--version' in sys.argv: print('Slow Fake Vina'); raise SystemExit(0)\n"
                "cfg=pathlib.Path(sys.argv[sys.argv.index('--config')+1]).read_text()\n"
                "if 'ligB.pdbqt' in cfg: time.sleep(3)\n"
                "out=[x.split(' = ',1)[1] for x in cfg.splitlines() if x.startswith('out = ')][0]\n"
                "pathlib.Path(out).write_text('REMARK VINA RESULT: -7.0 0.0 0.0\\n')\n",
                encoding="utf-8",
            )
            vina.chmod(0o755)
            results = root / "results"
            command = [sys.executable, str(ROOT / "3_Complete_batch_docking.py"), "--receptors", str(receptors),
                       "--ligands", str(ligands), "--centers_csv", str(root / "centers.csv"), "--poses", "1",
                       "--vina-exe", str(vina), "--results-dir", str(results), "--compact-results",
                       "--workers", "1", "--max-pending", "1", "--seed", "7"]
            process = subprocess.Popen(command, cwd=root, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            results_csv = results / "compact_results.csv"
            deadline = time.time() + 10
            while time.time() < deadline and (not results_csv.exists() or results_csv.stat().st_size == 0):
                time.sleep(0.05)
            self.assertTrue(results_csv.exists(), "first score was not durably collected")
            process.terminate()
            process.wait(timeout=10)
            # Let an orphaned fake worker finish before the temporary directory is removed.
            time.sleep(3.1)
            subprocess.run(command, cwd=root, check=True, capture_output=True, text=True)
            with results_csv.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(sorted(row["Ligand"] for row in rows), ["ligA", "ligB"])

    def test_timeout_is_recorded_per_ligand_and_trace_is_written(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            receptors = root / "Receptors"; receptors.mkdir()
            ligands = root / "Ligands"; ligands.mkdir()
            (receptors / "recA.pdbqt").write_text("RECEPTOR\n", encoding="utf-8")
            (ligands / "ligA.pdbqt").write_text("LIGAND\n", encoding="utf-8")
            (root / "centers.csv").write_text("PDB_ID,X,Y,Z\nrecA,1,2,3\n", encoding="utf-8")
            vina = root / "slow_fake_vina.py"
            vina.write_text(
                "#!/usr/bin/env python3\nimport sys,time\n"
                "if '--version' in sys.argv: print('Slow Fake Vina'); raise SystemExit(0)\n"
                "time.sleep(5)\n",
                encoding="utf-8",
            )
            vina.chmod(0o755)
            results, trace = root / "results", root / "memory.tsv"
            subprocess.run([
                sys.executable, str(ROOT / "3_Complete_batch_docking.py"),
                "--receptors", str(receptors), "--ligands", str(ligands),
                "--centers_csv", str(root / "centers.csv"), "--poses", "1",
                "--vina-exe", str(vina), "--results-dir", str(results), "--compact-results",
                "--workers", "1", "--max-pending", "1", "--timeout-seconds", "1",
                "--memory-trace", str(trace),
            ], cwd=root, check=True, capture_output=True, text=True)
            with (results / "compact_failures.csv").open(encoding="utf-8") as handle:
                failures = list(csv.DictReader(handle))
            self.assertEqual(len(failures), 1)
            self.assertIn("timeout 1s exceeded", failures[0]["Error"])
            self.assertIn("failed", trace.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
