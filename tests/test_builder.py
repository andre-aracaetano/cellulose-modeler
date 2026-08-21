import tempfile
import unittest
from pathlib import Path

from dynamics_cellulose.builder import build_structure
from dynamics_cellulose.geometry import distance
from dynamics_cellulose.pdbio import write_pdb
from dynamics_cellulose.reference import REFERENCE
from dynamics_cellulose.validation import validate


class BuilderTests(unittest.TestCase):
    def test_experimental_cell(self):
        self.assertAlmostEqual(REFERENCE.a, 7.784, places=3)
        self.assertAlmostEqual(REFERENCE.b, 8.201, places=3)
        self.assertAlmostEqual(REFERENCE.c, 10.380, places=3)
        self.assertAlmostEqual(REFERENCE.gamma, 96.55, places=2)
        self.assertEqual(set(REFERENCE.fractional_atoms), {0, 1})

    def test_single_chain_arbitrary_dp(self):
        report = validate(build_structure(7, [1]))
        self.assertEqual(report.chains, 1)
        self.assertEqual(report.residues, 7)
        self.assertAlmostEqual(report.glycosidic_min, 1.41502, places=4)

    def test_eighteen_chain_layout(self):
        structure = build_structure(20, [2, 3, 4, 4, 3, 2])
        report = validate(structure)
        self.assertEqual(report.chains, 18)
        self.assertEqual(report.residues, 360)
        self.assertEqual(report.atoms, 18 * (21 * 20 + 3))

    def test_all_hydroxyls_and_chain_ends_are_complete(self):
        structure = build_structure(3, [1])
        residues = structure.chains[0].residues
        for residue in residues:
            names = {atom.name for atom in residue.atoms}
            self.assertTrue({"HO2", "HO3", "HO6"} <= names)
            for oxygen, hydrogen in (("O2", "HO2"), ("O3", "HO3"), ("O6", "HO6")):
                self.assertAlmostEqual(
                    distance(residue.atom(oxygen).position, residue.atom(hydrogen).position),
                    0.96,
                    places=6,
                )
        self.assertIn("HO4", {atom.name for atom in residues[0].atoms})
        self.assertNotIn("HO4", {atom.name for atom in residues[1].atoms})
        self.assertIn("HO1", {atom.name for atom in residues[-1].atoms})
        self.assertEqual(sum(len(residue.atoms) for residue in residues), 21 * 3 + 3)

    def test_valid_pdb_fields(self):
        structure = build_structure(3, [2, 2])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test.pdb"
            write_pdb(path, structure)
            atom_lines = [line for line in path.read_text().splitlines() if line.startswith("ATOM")]
            self.assertTrue(atom_lines)
            self.assertTrue(all(len(line) >= 78 for line in atom_lines))
            self.assertEqual({line[21] for line in atom_lines}, {"A", "B", "C", "D"})
            self.assertTrue(any(line.startswith("CONECT") for line in path.read_text().splitlines()))
            self.assertEqual(sum(line.startswith("TER") for line in path.read_text().splitlines()), 4)


if __name__ == "__main__":
    unittest.main()
