import tempfile
import unittest
from pathlib import Path

import app


MMCIF_FIXTURE = """\
data_test
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . ARG A 1 59 ? -31.877 -7.496 -11.522 1.00 111.60 59 ARG A N 1
ATOM 2 C CA . ARG A 1 59 ? -31.236 -6.409 -12.316 1.00 107.83 59 ARG A CA 1
HETATM 3 C C1 . A1AKL A 1 801 ? 10.000 20.000 30.000 1.00 20.00 801 A1AKL A C1 1
HETATM 4 C C2 . A1AKL A 1 801 ? 12.000 22.000 32.000 1.00 20.00 801 A1AKL A C2 1
#
"""


class ReceptorPdbSnapshotTests(unittest.TestCase):
    def test_snapshot_from_mmcif_has_no_nul_bytes_and_aliases_long_hets(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "fixture.cif"
            dst = Path(tmp) / "fixture.pdb"
            src.write_text(MMCIF_FIXTURE, encoding="utf-8")

            alias_map = app._write_receptor_pdb_snapshot(src, dst, altloc_mode="collapse")

            text = dst.read_text(encoding="utf-8")
            self.assertNotIn("\x00", text)
            self.assertIn("ATOM", text)
            self.assertIn("HETATM", text)
            self.assertIn(" A1A ", text)
            self.assertEqual(alias_map.get("A1AKL"), "A1A")


if __name__ == "__main__":
    unittest.main()
