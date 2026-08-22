import tempfile
import unittest
from pathlib import Path
import subprocess
import sys

from cellulose_core.builder import build_structure
from cellulose_core.geometry import distance
from cellulose_core.pdbio import write_pdb
from cellulose_core.reference import REFERENCE
from cellulose_core.validation import validate


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

    def test_surface_c6_oxidation_is_reproducible_and_complete(self):
        structure = build_structure(
            20,
            [2, 3, 4, 4, 3, 2],
            oxidation_degree=0.25,
            oxidation_scope="surface", oxidation_seed=20260822,
        )
        report = validate(structure)
        oxidized = [
            (chain.number, residue.number)
            for chain in structure.chains
            for residue in chain.residues
            if residue.oxidized
        ]
        self.assertEqual(len(oxidized), 60)
        self.assertEqual(report.net_charge, -60)
        surface_chain_numbers = {1, 2, 3, 5, 6, 9, 10, 13, 14, 16, 17, 18}
        self.assertTrue(all(chain in surface_chain_numbers for chain, _ in oxidized))
        residue = next(
            residue for chain in structure.chains for residue in chain.residues if residue.oxidized
        )
        names = {atom.name for atom in residue.atoms}
        self.assertTrue({"C6", "O61", "O62"} <= names)
        self.assertFalse({"H61", "H62", "O6", "HO6", "HO62"} & names)
        self.assertAlmostEqual(
            distance(residue.atom("C6").position, residue.atom("O61").position), 1.26
        )

    def test_protonated_all_chain_oxidation(self):
        structure = build_structure(
            10,
            [2, 2],
            oxidation_degree=0.1,
            oxidation_scope="all",
            oxidation_protonated=True, oxidation_seed=7,
        )
        report = validate(structure)
        self.assertEqual(report.oxidized_sites, 4)
        self.assertEqual(report.net_charge, 0)
        oxidized = [r for c in structure.chains for r in c.residues if r.oxidized]
        self.assertTrue(all("HO62" in {atom.name for atom in r.atoms} for r in oxidized))

    def test_oxidation_degree_validation(self):
        with self.assertRaises(ValueError):
            build_structure(2, [1], oxidation_degree=0.09)

    def test_root_builder_creates_default_output_directory(self):
        root_builder = Path(__file__).resolve().parents[1] / "builder.py"
        with tempfile.TemporaryDirectory() as directory:
            result = subprocess.run(
                [sys.executable, str(root_builder), "--dp", "1"],
                cwd=directory,
                check=True,
                capture_output=True,
                text=True,
            )
            output = Path(directory) / "outputs" / "cellulose.pdb"
            self.assertTrue(output.is_file())
            self.assertIn("Atoms: 24", result.stdout)


if __name__ == "__main__":
    unittest.main()
