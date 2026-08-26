import importlib
import io
import json
import tempfile
import unittest
from pathlib import Path
import zipfile


class PublicAccessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        tmp_root = Path(cls._tmp.name)
        db_path = tmp_root / "test.db"
        cls.workspace_root = tmp_root / "workspaces"

        import app as app_module

        app_module.Config.SQLALCHEMY_DATABASE_URI = f"sqlite:///{db_path}"
        app_module.Config.TMP_ROOT = str(cls.workspace_root)
        app_module.Config.PUBLIC_EMAIL = "public@autodock.local"
        app_module.Config.DEFAULT_PACKAGE_MODE = "portable"
        app_module.Config.ENABLE_LSF_PACKAGE = True
        cls.app_module = importlib.reload(app_module)
        cls.app_module.Config.SQLALCHEMY_DATABASE_URI = f"sqlite:///{db_path}"
        cls.app_module.Config.TMP_ROOT = str(cls.workspace_root)
        cls.app_module.Config.PUBLIC_EMAIL = "public@autodock.local"
        cls.app_module.Config.DEFAULT_PACKAGE_MODE = "portable"
        cls.app_module.Config.ENABLE_LSF_PACKAGE = True

        cls.app = cls.app_module.create_app()
        cls.client = cls.app.test_client()

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def test_home_is_public(self):
        response = self.client.get("/")
        self.assertEqual(response.status_code, 200)
        self.assertIn(b"AutoDock-Vina PrepServer", response.data)
        self.assertNotIn(b"Sign in", response.data)

    def test_results_page_offers_example_button(self):
        response = self.client.get("/results")
        self.assertEqual(response.status_code, 200)
        self.assertIn(b"Open Example", response.data)
        self.assertIn(b"/viz/example", response.data)
        self.assertIn(b"filteredEntries", response.data)
        self.assertIn(b"findEntryIndexFromFrame", response.data)
        self.assertIn(b"Upload visualization ZIP", response.data)
        self.assertIn(b"uploadVisualizationZip", response.data)

    def test_login_route_redirects_to_home(self):
        response = self.client.get("/auth/login", follow_redirects=False)
        self.assertEqual(response.status_code, 302)
        self.assertEqual(response.headers["Location"], "/")

    def test_workspace_can_be_created_without_auth(self):
        response = self.client.post("/api/workspace")
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertTrue(payload["jobname"].endswith("-public"))

    def test_frontend_api_is_public_after_workspace_creation(self):
        workspace = self.client.post("/api/workspace").get_json()
        response = self.client.get("/api/receptors/list", query_string={"jobname": workspace["jobname"]})
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["receptors"], [])

    def test_wsfile_serves_receptor_from_workspace(self):
        workspace = self.client.post("/api/workspace").get_json()
        receptor = self.workspace_root / workspace["jobname"] / "Receptors" / "3eky.pdb"
        receptor.parent.mkdir(parents=True, exist_ok=True)
        receptor.write_text("ATOM      1  N   GLY A   1      11.104  13.207  10.451  1.00 20.00           N\n")

        response = self.client.get(
            "/api/wsfile",
            query_string={"jobname": workspace["jobname"], "rel": "Receptors/3eky.pdb"},
        )
        self.assertEqual(response.status_code, 200)
        self.assertIn(b"ATOM", response.data)
        response.close()

    def test_wsfile_falls_back_to_receptors_pdb_snapshot(self):
        workspace = self.client.post("/api/workspace").get_json()
        receptor = self.workspace_root / workspace["jobname"] / "Receptors_PDB" / "3eky.pdb"
        receptor.parent.mkdir(parents=True, exist_ok=True)
        receptor.write_text("ATOM      1  N   GLY A   1      11.104  13.207  10.451  1.00 20.00           N\n")

        response = self.client.get(
            "/api/wsfile",
            query_string={"jobname": workspace["jobname"], "rel": "Receptors/3eky.pdb"},
        )
        self.assertEqual(response.status_code, 200)
        self.assertIn(b"ATOM", response.data)
        response.close()

    def test_wsfile_blocks_path_traversal(self):
        workspace = self.client.post("/api/workspace").get_json()
        response = self.client.get(
            "/api/wsfile",
            query_string={"jobname": workspace["jobname"], "rel": "../README.md"},
        )
        self.assertEqual(response.status_code, 404)

    def test_visualization_project_page_renders_workspace_manifest(self):
        workspace = self.client.post("/api/workspace").get_json()
        ws = self.workspace_root / workspace["jobname"]
        viz_dir = ws / "Docking_HTML_Viz_Project"
        viewer_dir = viz_dir / "viewers"
        viewer_dir.mkdir(parents=True)
        (viewer_dir / "entry.html").write_text("<html><body>viewer</body></html>", encoding="utf-8")
        (viz_dir / "manifest.json").write_text(
            json.dumps(
                {
                    "project_name": "Docking_HTML_Viz_Project",
                    "page_title": "Docking Viz",
                    "source_csv": "/tmp/demo.csv",
                    "entry_count": 1,
                    "entries": [
                        {
                            "receptor": "recA",
                            "ligand": "LigA",
                            "best_affinity": "-9.1",
                            "viewer_file": "viewers/entry.html",
                        }
                    ],
                    "attribution": {
                        "repo": "https://github.com/muntisa/py-VinaScope-Docking-Viewer",
                        "viewer": "https://muntisa.github.io/VinaDock-Viz/VinaDock_Viz.html",
                        "author": "Cristian R. Munteanu, PhD",
                        "affiliation": "Professor of Computer Science, University of A Coruna, RNASA-IMEDIR",
                    },
                }
            ),
            encoding="utf-8",
        )

        response = self.client.get(
            "/viz/project",
            query_string={"jobname": workspace["jobname"], "rel": "Docking_HTML_Viz_Project/manifest.json"},
        )
        self.assertEqual(response.status_code, 200)
        self.assertIn(b"Docking Viz", response.data)
        self.assertIn(b"py-VinaScope-Docking-Viewer", response.data)
        self.assertIn(b'view=full', response.data)
        self.assertIn(b'data-receptor="recA"', response.data)
        self.assertIn(b"filteredEntryModels", response.data)

        inline = self.client.get(
            "/api/wsinline",
            query_string={"jobname": workspace["jobname"], "rel": "Docking_HTML_Viz_Project/viewers/entry.html"},
        )
        self.assertEqual(inline.status_code, 200)
        self.assertIn(b"viewer", inline.data)

        viewer_file = self.client.get(
            f"/viz/file/{workspace['jobname']}/Docking_HTML_Viz_Project/viewers/entry.html"
        )
        self.assertEqual(viewer_file.status_code, 200)
        self.assertIn(b"viewer", viewer_file.data)
        self.assertIn("sandbox allow-scripts", viewer_file.headers.get("Content-Security-Policy", ""))
        self.assertEqual(viewer_file.headers.get("X-Content-Type-Options"), "nosniff")

        standalone = self.client.get(
            "/viz/project",
            query_string={"jobname": workspace["jobname"], "rel": "Docking_HTML_Viz_Project/manifest.json", "entry": 0, "view": "full"},
        )
        self.assertEqual(standalone.status_code, 200)
        self.assertIn(b"Open raw viewer", standalone.data)

    def test_visualization_zip_upload_creates_server_project(self):
        manifest = {
            "project_name": "Docking_HTML_Viz_Project_Compacted_SDF",
            "page_title": "Uploaded Docking Viz",
            "source_csv": "scores.csv",
            "entry_count": 1,
            "entries": [
                {
                    "receptor": "recA",
                    "ligand": "LigA",
                    "best_affinity": "-9.1",
                    "viewer_file": "viewers/entry.html",
                }
            ],
            "attribution": {
                "repo": "https://github.com/muntisa/py-VinaScope-Docking-Viewer",
                "viewer": "https://muntisa.github.io/VinaDock-Viz/VinaDock_Viz.html",
                "author": "Cristian R. Munteanu, PhD",
                "affiliation": "Professor of Computer Science, University of A Coruna, RNASA-IMEDIR",
            },
        }
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("Docking_HTML_Viz_Project_Compacted_SDF/manifest.json", json.dumps(manifest))
            zf.writestr("Docking_HTML_Viz_Project_Compacted_SDF/viewers/entry.html", "<html><body>uploaded viewer</body></html>")
            zf.writestr("Docking_HTML_Viz_Project_Compacted_SDF/inputs/receptor_library.js", "window.VINA_RECEPTOR_LIBRARY={};")
        buf.seek(0)

        response = self.client.post(
            "/api/results/upload",
            data={"project_zip": (buf, "Docking_HTML_Viz_Project_Compacted_SDF.zip")},
            content_type="multipart/form-data",
        )
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertTrue(payload["ok"])
        self.assertEqual(payload["entry_count"], 1)
        self.assertIn("/viz/project", payload["project_url"])

        project = self.client.get(payload["project_url"])
        self.assertEqual(project.status_code, 200)
        self.assertIn(b"Uploaded Docking Viz", project.data)
        self.assertIn(b"/viz/file/", project.data)

    def test_visualization_zip_upload_rejects_path_traversal(self):
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("../escape.html", "nope")
        buf.seek(0)
        response = self.client.post(
            "/api/results/upload",
            data={"project_zip": (buf, "bad.zip")},
            content_type="multipart/form-data",
        )
        self.assertEqual(response.status_code, 400)
        payload = response.get_json()
        self.assertEqual(payload["error"], "invalid_visualization_zip")
        self.assertIn("Unsafe ZIP path", payload["message"])

    def test_visualization_zip_upload_rejects_symlink_entry(self):
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
            info = zipfile.ZipInfo("Docking_HTML_Viz_Project/link")
            info.create_system = 3
            info.external_attr = 0o120777 << 16
            zf.writestr(info, "manifest.json")
        buf.seek(0)
        response = self.client.post(
            "/api/results/upload",
            data={"project_zip": (buf, "symlink.zip")},
            content_type="multipart/form-data",
        )
        self.assertEqual(response.status_code, 400)
        payload = response.get_json()
        self.assertEqual(payload["error"], "invalid_visualization_zip")
        self.assertIn("Symlinks are not allowed", payload["message"])

    def test_public_example_visualization_project_renders(self):
        response = self.client.get("/viz/example")
        self.assertEqual(response.status_code, 200)
        self.assertIn(b"Docking Visualization Project", response.data)
        self.assertIn(b"Open standalone", response.data)
        self.assertIn(b'view=full', response.data)

    def test_ligand_zip_upload_flattens_nested_ligands_folder(self):
        workspace = self.client.post("/api/workspace").get_json()
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("Ligands/DR7.sdf", "dr7")
            zf.writestr("__MACOSX/._DR7.sdf", "junk")
        buf.seek(0)

        response = self.client.post(
            "/api/ligands/upload",
            data={"jobname": workspace["jobname"], "mode": "zip", "file": (buf, "Ligands.zip")},
            content_type="multipart/form-data",
        )
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["accepted_files"], ["DR7.sdf"])
        lig_dir = self.workspace_root / workspace["jobname"] / "Ligands"
        self.assertTrue((lig_dir / "DR7.sdf").exists())
        self.assertFalse((lig_dir / "Ligands").exists())

    def test_ligand_folder_upload_accepts_multiple_files(self):
        workspace = self.client.post("/api/workspace").get_json()
        response = self.client.post(
            "/api/ligands/upload",
            data={
                "jobname": workspace["jobname"],
                "mode": "folder",
                "folder_name": "Ligands",
                "files": [
                    (io.BytesIO(b"dr7"), "Ligands/DR7.sdf"),
                    (io.BytesIO(b"dr8"), "Ligands/DR8.sdf"),
                    (io.BytesIO(b"junk"), "Ligands/.DS_Store"),
                ],
            },
            content_type="multipart/form-data",
        )
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(sorted(payload["accepted_files"]), ["DR7.sdf", "DR8.sdf"])
        lig_dir = self.workspace_root / workspace["jobname"] / "Ligands"
        self.assertTrue((lig_dir / "DR7.sdf").exists())
        self.assertTrue((lig_dir / "DR8.sdf").exists())
        self.assertFalse((lig_dir / ".DS_Store").exists())


if __name__ == "__main__":
    unittest.main()
