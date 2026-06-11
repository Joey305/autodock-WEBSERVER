import importlib.util
import io
import tempfile
import unittest
import zipfile
import importlib
from pathlib import Path

from packager import normalize_ligand_tree, save_uploaded_ligand_zip


REPO_ROOT = Path(__file__).resolve().parent.parent


def load_script_module(filename: str, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, REPO_ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class _UploadStub:
    def __init__(self, payload: bytes):
        self._payload = payload

    def read(self) -> bytes:
        return self._payload


def make_zip(entries):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for name, content in entries.items():
            zf.writestr(name, content)
    return buf.getvalue()


class LigandUploadNormalizationTests(unittest.TestCase):
    def test_zip_with_top_level_ligands_folder_is_flattened(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            dest = Path(tmpdir) / "Ligands"
            payload = make_zip({
                "Ligands/DR7.sdf": "dr7",
                "Ligands/DR8.sdf": "dr8",
            })
            result = save_uploaded_ligand_zip(_UploadStub(payload), dest)
            self.assertEqual(sorted(result["accepted_files"]), ["DR7.sdf", "DR8.sdf"])
            self.assertTrue((dest / "DR7.sdf").exists())
            self.assertFalse((dest / "Ligands").exists())

    def test_zip_with_macos_artifacts_ignores_hidden_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            dest = Path(tmpdir) / "Ligands"
            payload = make_zip({
                "__MACOSX/._DR7.sdf": "junk",
                "Ligands/.DS_Store": "junk",
                "Ligands/DR7.sdf": "dr7",
            })
            result = save_uploaded_ligand_zip(_UploadStub(payload), dest)
            self.assertEqual(result["accepted_files"], ["DR7.sdf"])
            self.assertFalse((dest / "__MACOSX").exists())
            self.assertFalse((dest / ".DS_Store").exists())

    def test_zip_ignores_unsupported_files_with_warning(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            dest = Path(tmpdir) / "Ligands"
            payload = make_zip({
                "Ligands/DR7.sdf": "dr7",
                "Ligands/readme.txt": "ignore me",
            })
            result = save_uploaded_ligand_zip(_UploadStub(payload), dest)
            self.assertEqual(result["accepted_files"], ["DR7.sdf"])
            self.assertIn("Ligands/readme.txt", result["ignored_files"])
            self.assertTrue(any("readme.txt" in warning for warning in result["warnings"]))

    def test_zip_duplicate_filenames_are_renamed(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            dest = Path(tmpdir) / "Ligands"
            payload = make_zip({
                "DR7.sdf": "root",
                "Ligands/DR7.sdf": "nested",
            })
            result = save_uploaded_ligand_zip(_UploadStub(payload), dest)
            self.assertEqual(sorted(result["accepted_files"]), ["DR7.sdf", "DR7_002.sdf"])
            self.assertTrue((dest / "DR7.sdf").exists())
            self.assertTrue((dest / "DR7_002.sdf").exists())

    def test_assemble_normalization_flattens_existing_nested_workspace(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            src = Path(tmpdir) / "workspace_ligands"
            src.mkdir()
            (src / "Ligands").mkdir()
            (src / "Ligands" / "DR7.sdf").write_text("dr7")
            dest = Path(tmpdir) / "job_ligands"
            result = normalize_ligand_tree(src, dest)
            self.assertEqual(result["accepted_files"], ["DR7.sdf"])
            self.assertTrue((dest / "DR7.sdf").exists())
            self.assertFalse((dest / "Ligands").exists())

    def test_recursive_discovery_finds_nested_ligands(self):
        if importlib.util.find_spec("rdkit") is None:
            self.skipTest("rdkit not installed")
        module = load_script_module("1_ConformerGeneration.py", "confgen_discovery_module")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "Ligands"
            (root / "Ligands").mkdir(parents=True)
            (root / "Ligands" / "DR7.sdf").write_text("dr7")
            (root / "__MACOSX").mkdir()
            (root / "__MACOSX" / "._DR7.sdf").write_text("junk")
            discovered = module.discover_ligand_input_files(root, (".sdf",))
            self.assertEqual([path.name for path in discovered], ["DR7.sdf"])
            self.assertEqual(discovered[0].parent.name, "Ligands")


if __name__ == "__main__":
    unittest.main()
