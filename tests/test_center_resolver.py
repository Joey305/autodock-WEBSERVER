import tempfile
import unittest
from pathlib import Path

from center_resolver import CenterResolutionError, resolve_center_from_file


FIXTURE = """\
HETATM    1  C1  DR7 A 100      10.000  20.000  30.000  1.00 20.00           C
HETATM    2  C2  DR7 A 100      12.000  22.000  32.000  1.00 20.00           C
ATOM      3  CA  GLU A  82       1.000   2.000   3.000  1.00 20.00           C
ATOM      4  CB  GLU A  82       3.000   4.000   5.000  1.00 20.00           C
HETATM    5  C1  ATP B 301       2.000   4.000   6.000  1.00 20.00           C
HETATM    6  C2  ATP B 301       4.000   6.000   8.000  1.00 20.00           C
HETATM    7  C1  ATP C 301      20.000  40.000  60.000  1.00 20.00           C
"""


class CenterResolverTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.path = Path(self._tmp.name) / "fixture.pdb"
        self.path.write_text(FIXTURE)

    def tearDown(self):
        self._tmp.cleanup()

    def test_explicit_xyz_validation(self):
        result = resolve_center_from_file(self.path, {"method": "xyz", "center": [1, "2.5", 3], "size": 18})
        self.assertEqual(result["center"], [1.0, 2.5, 3.0])
        self.assertEqual(result["size"], 18.0)

    def test_residue_center_by_chain_and_resi(self):
        result = resolve_center_from_file(self.path, {"method": "residue", "chain": "A", "resi": "82"})
        self.assertEqual(result["center"], [2.0, 3.0, 4.0])
        self.assertEqual(result["matched"]["atom_count"], 2)

    def test_residue_center_by_chain_resi_resname(self):
        result = resolve_center_from_file(self.path, {"method": "residue", "chain": "A", "resi": "82", "resname": "GLU"})
        self.assertEqual(result["matched"]["resname"], "GLU")
        self.assertEqual(result["center"], [2.0, 3.0, 4.0])

    def test_hetatm_center_by_resname_chain_resi(self):
        result = resolve_center_from_file(self.path, {"method": "hetatm", "het": "DR7", "chain": "A", "resi": "100"})
        self.assertEqual(result["center"], [11.0, 21.0, 31.0])
        self.assertEqual(result["matched"]["atom_count"], 2)

    def test_no_match_error(self):
        with self.assertRaises(CenterResolutionError) as ctx:
            resolve_center_from_file(self.path, {"method": "hetatm", "het": "DR7", "chain": "Z", "resi": "100"})
        self.assertEqual(ctx.exception.error, "no_atoms_matched")

    def test_ambiguous_ligand_error(self):
        with self.assertRaises(CenterResolutionError) as ctx:
            resolve_center_from_file(self.path, {"method": "hetatm", "het": "ATP"})
        self.assertEqual(ctx.exception.error, "ambiguous_selection")
        self.assertEqual(len(ctx.exception.details["candidates"]), 2)

    def test_atom_specific_match(self):
        result = resolve_center_from_file(
            self.path,
            {"method": "atom", "record": "HETATM", "resname": "DR7", "chain": "A", "resi": "100", "atom_name": "C2"},
        )
        self.assertEqual(result["center"], [12.0, 22.0, 32.0])


if __name__ == "__main__":
    unittest.main()
