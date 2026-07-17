import json
import tempfile
import unittest
from pathlib import Path

from hpc_profiles import (
    JOEY_LSF_PROFILE,
    build_custom_profile,
    load_packaged_profile,
    normalize_package_mode,
    render_lsf_header,
    render_setup_block,
    save_packaged_profile,
)
from lsf_templates import build_confgen_lsfs, build_vina_lsfs


class HpcProfileTests(unittest.TestCase):
    def test_mode_aliases_remain_backward_compatible(self):
        self.assertEqual(normalize_package_mode({"package_mode": "lsf"}), "joey_lsf")
        self.assertEqual(normalize_package_mode({"include_lsf": "1"}), "joey_lsf")
        self.assertEqual(normalize_package_mode({"package_mode": "custom_lsf"}), "custom_lsf")

    def test_custom_profile_is_sanitized(self):
        profile = build_custom_profile(
            {
                "lsf_email": "cluster@example.org",
                "notify_begin": "1",
                "queue": "long",
                "project": "",
                "workers": "24",
                "mem_per_core": "4096",
                "confgen_walltime": "12:00",
                "vina_walltime": "72:00",
                "conda_sh": "/apps/conda.sh",
                "conda_env": "dock",
                "vina_path": "/apps/vina",
                "python_command": "/usr/bin/python3",
                "setup_commands": "module load gcc\nmodule load cuda",
            }
        )
        self.assertEqual(profile.email, "cluster@example.org")
        self.assertFalse(profile.project)
        self.assertEqual(profile.workers, 24)
        self.assertEqual(profile.setup_commands, ("module load gcc", "module load cuda"))

    def test_render_header_omits_optional_project_and_notifications(self):
        profile = build_custom_profile({"queue": "short", "workers": "8"})
        header = render_lsf_header(profile=profile, jobname="demo", log_prefix="demo", walltime="04:00")
        self.assertNotIn("#BSUB -P", header)
        self.assertNotIn("#BSUB -u", header)
        self.assertIn("#BSUB -q short", header)

    def test_profile_round_trips_inside_package(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            save_packaged_profile(root, JOEY_LSF_PROFILE)
            loaded = load_packaged_profile(root)
        self.assertIsNotNone(loaded)
        self.assertEqual(loaded.email, JOEY_LSF_PROFILE.email)

    def test_generated_lsf_files_persist_custom_profile(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "Ligands").mkdir()
            (root / "Receptors").mkdir()
            (root / "Ligands" / "ligands.csv").write_text("smiles,id\nCCO,lig1\n", encoding="utf-8")
            (root / "vina_centers.csv").write_text("PDB_ID,X,Y,Z,SIZE\nrec.pdbqt,1,2,3,20\n", encoding="utf-8")

            profile = build_custom_profile(
                {
                    "lsf_email": "cluster@example.org",
                    "notify_end": "1",
                    "queue": "gpu",
                    "project": "",
                    "workers": "32",
                    "mem_per_core": "6000",
                    "confgen_walltime": "10:00",
                    "vina_walltime": "50:00",
                    "python_command": "/opt/python/bin/python3",
                    "setup_commands": "module load anaconda",
                }
            )

            build_confgen_lsfs(
                root,
                root,
                profile=profile,
                poses=64,
                lig_mode="1",
                lig_filetype="csv",
                csv_smiles_col="smiles",
                csv_id_col="id",
                single_sdf_rel=None,
            )
            build_vina_lsfs(root, root, profile=profile, poses=20)

            packaged = json.loads((root / "hpc_profile.json").read_text(encoding="utf-8"))
            confgen = (root / "run_confgen_job.lsf").read_text(encoding="utf-8")
            vina = (root / "run_vina_job.lsf").read_text(encoding="utf-8")

        self.assertEqual(packaged["queue"], "gpu")
        self.assertEqual(packaged["python_command"], "/opt/python/bin/python3")
        self.assertIn('PYBIN="/opt/python/bin/python3"', confgen)
        self.assertIn("module load anaconda", confgen)
        self.assertNotIn("#BSUB -P", confgen)
        self.assertIn("#BSUB -u cluster@example.org", vina)
        self.assertIn('"$PYBIN" 3_Complete_batch_docking.py', vina)


if __name__ == "__main__":
    unittest.main()
