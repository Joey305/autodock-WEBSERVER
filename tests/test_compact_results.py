import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path


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


if __name__ == "__main__":
    unittest.main()
