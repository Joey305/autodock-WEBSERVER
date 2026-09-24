import json
import tempfile
import unittest
from pathlib import Path

from hpc_profiles import (
    JOEY_LSF_PROFILE,
    MAINAK_LSF_PROFILE,
    build_custom_profile,
    load_packaged_profile,
    normalize_package_mode,
    profile_for_mode,
    replace_profile,
    render_lsf_header,
    render_setup_block,
    save_packaged_profile,
    validate_hpc_profile,
)
from lsf_templates import build_compacted_sdf_html_lsf, build_confgen_lsfs, build_pymol_lsf, build_vina_lsfs


class HpcProfileTests(unittest.TestCase):
    def test_mode_aliases_remain_backward_compatible(self):
        self.assertEqual(normalize_package_mode({"package_mode": "lsf"}), "joey_lsf")
        self.assertEqual(normalize_package_mode({"include_lsf": "1"}), "joey_lsf")
        self.assertEqual(normalize_package_mode({"package_mode": "mainak_lsf"}), "mainak_lsf")
        self.assertEqual(normalize_package_mode({"package_mode": "custom_lsf"}), "custom_lsf")

    def test_mainak_mode_uses_mainak_profile(self):
        profile = profile_for_mode("mainak_lsf")
        self.assertEqual(profile, MAINAK_LSF_PROFILE)
        self.assertEqual(profile.email, "mxb2638@miami.edu")
        self.assertEqual(profile.conda_sh, "/nethome/mxb2638/miniconda3/etc/profile.d/conda.sh")

    def test_joey_pegasus_profile_has_verified_allocations(self):
        self.assertEqual(JOEY_LSF_PROFILE.profile_name, "joey_pegasus")
        self.assertEqual(JOEY_LSF_PROFILE.queue, "gpu_cheminfo")
        self.assertEqual(JOEY_LSF_PROFILE.project, "brd")
        self.assertEqual(JOEY_LSF_PROFILE.effective_vina_cpus, 16)
        self.assertEqual(JOEY_LSF_PROFILE.effective_confgen_cpus, 16)
        self.assertEqual(JOEY_LSF_PROFILE.effective_confgen_workers, 16)
        self.assertEqual(JOEY_LSF_PROFILE.vina_walltime, "240:00")
        self.assertEqual(JOEY_LSF_PROFILE.confgen_walltime, "48:00")
        validate_hpc_profile(JOEY_LSF_PROFILE)

    def test_confgen_workers_cannot_exceed_lsf_cpu_request(self):
        invalid = replace_profile(JOEY_LSF_PROFILE, confgen_cpus=16, confgen_workers=17)
        with self.assertRaisesRegex(ValueError, "cannot exceed"):
            validate_hpc_profile(invalid)

    def test_custom_profile_is_sanitized(self):
        profile = build_custom_profile(
            {
                "lsf_email": "cluster@example.org",
                "notify_begin": "1",
                "queue": "long",
                "project": "account",
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
        self.assertEqual(profile.project, "account")
        self.assertEqual(profile.workers, 24)
        self.assertEqual(profile.setup_commands, ("module load gcc", "module load cuda"))

    def test_render_header_keeps_project_and_omits_notifications_when_not_requested(self):
        profile = build_custom_profile({"queue": "short", "project": "account", "workers": "8"})
        header = render_lsf_header(profile=profile, jobname="demo", log_prefix="demo", walltime="04:00")
        self.assertIn("#BSUB -P account", header)
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
                    "project": "account",
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
            build_pymol_lsf(root, root, profile=profile)
            build_compacted_sdf_html_lsf(root, root, profile=profile)

            packaged = json.loads((root / "hpc_profile.json").read_text(encoding="utf-8"))
            confgen = (root / "run_confgen_job.lsf").read_text(encoding="utf-8")
            vina = (root / "run_vina_job.lsf").read_text(encoding="utf-8")
            pymol = (root / "run_pymol_job.lsf").read_text(encoding="utf-8")
            compacted_sdf = (root / "run_compacted_sdf_html_job.lsf").read_text(encoding="utf-8")

        self.assertEqual(packaged["queue"], "gpu")
        self.assertEqual(packaged["python_command"], "/opt/python/bin/python3")
        self.assertIn('PYBIN="/opt/python/bin/python3"', confgen)
        self.assertIn("module load anaconda", confgen)
        self.assertIn("#BSUB -P account", confgen)
        self.assertIn("#BSUB -u cluster@example.org", vina)
        self.assertIn('"$PYBIN" 3_Complete_batch_docking.py', vina)
        self.assertIn('"$PYBIN" 5C_BuildPymolSesh.py', pymol)
        self.assertIn("--non-interactive", pymol)
        self.assertIn('"$PYBIN" 5_COMPACTED_SDF_HTML.py', compacted_sdf)
        self.assertIn("--top-ligands", compacted_sdf)

    def test_joey_generation_uses_workload_specific_lsf_settings(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "Ligands").mkdir()
            (root / "Receptors").mkdir()
            (root / "Ligands" / "ligands.csv").write_text("smiles,id\nCCO,lig1\n", encoding="utf-8")
            (root / "vina_centers.csv").write_text("PDB_ID,X,Y,Z,SIZE\nrec.pdbqt,1,2,3,20\n", encoding="utf-8")
            build_confgen_lsfs(root, root, profile=JOEY_LSF_PROFILE, poses=64, lig_mode="1", lig_filetype="csv", csv_smiles_col="smiles", csv_id_col="id", single_sdf_rel=None)
            build_vina_lsfs(root, root, profile=JOEY_LSF_PROFILE, poses=20)
            confgen = (root / "run_confgen_job.lsf").read_text(encoding="utf-8")
            vina = (root / "run_vina_job.lsf").read_text(encoding="utf-8")

        self.assertIn("#BSUB -P brd", vina)
        self.assertIn("#BSUB -W 240:00", vina)
        self.assertIn("#BSUB -q gpu_cheminfo", vina)
        self.assertIn("#BSUB -n 16", vina)
        self.assertIn('#BSUB -R "span[hosts=1]"', vina)
        self.assertIn('#BSUB -R "rusage[mem=2000]"', vina)
        self.assertIn("#BSUB -W 48:00", confgen)
        self.assertIn("#BSUB -n 16", confgen)
        self.assertIn("--workers 16", confgen)


if __name__ == "__main__":
    unittest.main()
