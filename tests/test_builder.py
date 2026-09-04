import tempfile
import unittest
from pathlib import Path
import subprocess
import sys
from collections import defaultdict

from cellulose_modeler.builder import build_structure, _paajanen_eligible_residues
from cellulose_modeler.geometry import distance
from cellulose_modeler.pdbio import write_pdb
from cellulose_modeler.reference import REFERENCE, REFERENCES
from cellulose_modeler.validation import validate


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

    def test_allomorph_reference_cells_and_chain_directions(self):
        expected = {
            "ialpha": ((10.400, 6.717, 5.962, 80.37, 118.08, 114.80), [(1, 0, 0)]),
            "ibeta": ((7.784, 8.201, 10.380, 90.0, 90.0, 96.55), [(0, 0, 1), (0, 0, 1)]),
            "ii": ((8.10, 9.03, 10.31, 90.0, 90.0, 117.10), [(0, 0, 1), (0, 0, -1)]),
            "iii-i": ((4.450, 7.850, 10.310, 90.0, 90.0, 105.10), [(0, 0, 1)]),
        }
        for key, (cell, cycles) in expected.items():
            reference = REFERENCES[key]
            self.assertEqual(reference.cycle_translations, cycles)
            for observed, target in zip(
                (reference.a, reference.b, reference.c, reference.alpha,
                 reference.beta, reference.gamma), cell
            ):
                self.assertAlmostEqual(observed, target, places=3)

    def test_all_allomorphs_build_valid_18_chain_dp20_structures(self):
        expected_links = {
            "ialpha": (1.402, 1.406),
            "ibeta": (1.415, 1.429),
            "ii": (1.393, 1.395),
            "iii-i": (1.377, 1.379),
        }
        for key, (lower, upper) in expected_links.items():
            structure = build_structure(20, [2, 3, 4, 4, 3, 2], allomorph=key)
            report = validate(structure)
            self.assertEqual(report.chains, 18)
            self.assertEqual(report.residues, 360)
            self.assertEqual(report.atoms, 7614)
            self.assertGreaterEqual(report.glycosidic_min, lower)
            self.assertLessEqual(report.glycosidic_max, upper)

    def test_cellulose_ii_contains_antiparallel_chains(self):
        structure = build_structure(20, [1, 1], allomorph="ii")
        directions = []
        for chain in structure.chains:
            directions.append(
                chain.residues[-1].atom("C1").position[2]
                - chain.residues[0].atom("C1").position[2]
            )
        self.assertGreater(directions[0], 0.0)
        self.assertLess(directions[1], 0.0)

    def test_non_ibeta_oxidation_is_rejected_until_validated(self):
        for key in ("ialpha", "ii", "iii-i"):
            with self.assertRaisesRegex(ValueError, "validated only for cellulose I-beta"):
                build_structure(20, [1], allomorph=key, oxidation_degree=0.25)

    def test_cli_rejects_unvalidated_charmm_gui_allomorphs(self):
        for allomorph in ("ii", "iii-i"):
            result = subprocess.run(
                [sys.executable, "builder.py", "--allomorph", allomorph,
                 "--charmm-gui"],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("does not have a validated CHARMM-GUI", result.stderr)

    def test_eighteen_chain_layout(self):
        structure = build_structure(20, [2, 3, 4, 4, 3, 2])
        report = validate(structure)
        self.assertEqual(report.chains, 18)
        self.assertEqual(report.residues, 360)
        self.assertEqual(report.atoms, 7614)

    def test_unmodified_18x20_has_intact_c6_on_every_residue(self):
        structure = build_structure(20, [2, 3, 4, 4, 3, 2])
        self.assertEqual(sum(len(chain.residues) for chain in structure.chains), 360)
        for chain in structure.chains:
            for residue in chain.residues:
                names = [atom.name for atom in residue.atoms]
                self.assertEqual(len(names), len(set(names)))
                self.assertTrue({"C6", "H61", "H62", "O6", "HO6"} <= set(names))
                self.assertNotIn("H63", names)
                self.assertFalse(residue.oxidized)

    def test_chain_13_even_residues_keep_primary_alcohol(self):
        chain = build_structure(20, [2, 3, 4, 4, 3, 2]).chains[12]
        self.assertEqual(chain.number, 13)
        for residue in chain.residues[1::2]:
            self.assertEqual(residue.number % 2, 0)
            names = {atom.name for atom in residue.atoms}
            self.assertTrue({"H61", "H62", "O6", "HO6"} <= names)
            self.assertNotIn("H63", names)

    def test_small_multichain_structure_has_no_parity_dependent_chemistry(self):
        structure = build_structure(4, [2, 2])
        signatures = {
            (chain.number, residue.number): {
                atom.name for atom in residue.atoms
                if atom.name in {"H61", "H62", "H63", "O6", "HO6"}
            }
            for chain in structure.chains for residue in chain.residues
        }
        self.assertEqual(set(map(frozenset, signatures.values())), {
            frozenset({"H61", "H62", "O6", "HO6"})
        })

    def test_unmodified_pdb_c6_conect_and_output_are_deterministic(self):
        structure = build_structure(20, [2, 3, 4, 4, 3, 2])
        with tempfile.TemporaryDirectory() as directory:
            first = Path(directory) / "first.pdb"
            second = Path(directory) / "second.pdb"
            write_pdb(first, structure)
            write_pdb(second, build_structure(20, [2, 3, 4, 4, 3, 2]))
            self.assertEqual(first.read_bytes(), second.read_bytes())
            lines = first.read_text().splitlines()
            atom_lines = [line for line in lines if line.startswith("ATOM")]
            self.assertEqual(len(atom_lines), 7614)
            self.assertEqual({line[17:21].strip() for line in atom_lines}, {"BGLC"})
            self.assertFalse(any(line[12:16].strip() == "H63" for line in atom_lines))

            atoms = {int(line[6:11]): line for line in atom_lines}
            adjacency = defaultdict(set)
            for line in lines:
                if line.startswith("CONECT"):
                    values = [int(line[i:i + 5]) for i in range(6, len(line), 5)]
                    adjacency[values[0]].update(values[1:])
            by_key = {
                (line[21], int(line[22:26]), line[12:16].strip()): serial
                for serial, line in atoms.items()
            }
            for residue_number in range(2, 21, 2):
                c6 = by_key[("M", residue_number, "C6")]
                o6 = by_key[("M", residue_number, "O6")]
                ho6 = by_key[("M", residue_number, "HO6")]
                self.assertIn(o6, adjacency[c6])
                self.assertIn(ho6, adjacency[o6])

    def test_charmm_gui_profile_prevents_quinovose_signature(self):
        structure = build_structure(20, [2, 3, 4, 4, 3, 2])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "charmm_gui.pdb"
            write_pdb(path, structure, charmm_gui=True)
            lines = path.read_text().splitlines()
            atoms = [line for line in lines if line.startswith("ATOM")]

            self.assertEqual(len(atoms), 6498)
            names = [line[12:16].strip() for line in atoms]
            self.assertFalse({"HO1", "HO2", "HO3", "HO4", "HO6", "H63"} & set(names))
            self.assertEqual(names.count("O6"), 360)
            self.assertEqual(names.count("H61"), 360)
            self.assertEqual(names.count("H62"), 360)

            by_key = {
                (line[21], int(line[22:26]), line[12:16].strip()): line
                for line in atoms
            }
            # Residue 1 is the original high-z end after whole-block reversal.
            self.assertGreater(float(by_key[("M", 1, "C1")][46:54]), 90.0)
            self.assertLess(float(by_key[("M", 20, "C1")][46:54]), 10.0)

            serial_key = {
                int(line[6:11]): (line[21], int(line[22:26]), line[12:16].strip())
                for line in atoms
            }
            edges = set()
            for line in lines:
                if line.startswith("CONECT"):
                    values = [int(line[i:i + 5]) for i in range(6, len(line), 5)]
                    edges.update(tuple(sorted((values[0], other))) for other in values[1:])
            key_serial = {value: key for key, value in serial_key.items()}
            for chain_id in "ABCDEFGHIJKLMNOPQR":
                for residue_number in range(1, 20):
                    expected = tuple(sorted((
                        key_serial[(chain_id, residue_number, "O4")],
                        key_serial[(chain_id, residue_number + 1, "C1")],
                    )))
                    self.assertIn(expected, edges)

    def test_charmm_gui_profile_preserves_oxidized_c6_heavy_atoms(self):
        structure = build_structure(
            20,
            [2, 3, 4, 4, 3, 2],
            oxidation_degree=0.25,
            oxidation_scope="surface",
            oxidation_seed=344,
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tocnf25_charmm_gui.pdb"
            write_pdb(path, structure, charmm_gui=True)
            atoms = [
                line for line in path.read_text().splitlines()
                if line.startswith("ATOM")
            ]
            names = [line[12:16].strip() for line in atoms]
            self.assertEqual(len(atoms), 6438)
            self.assertEqual(names.count("O6"), 300)
            self.assertEqual(names.count("O61"), 60)
            self.assertEqual(names.count("O62"), 60)
            self.assertEqual(names.count("H63"), 0)
            self.assertFalse({"HO1", "HO2", "HO3", "HO4", "HO6", "HO62"} & set(names))

    def test_charmm_gui_profile_retains_explicit_cooh_protonation(self):
        structure = build_structure(
            10,
            [2, 2],
            oxidation_degree=0.1,
            oxidation_scope="all",
            oxidation_protonated=True,
            oxidation_seed=7,
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "protonated_charmm_gui.pdb"
            write_pdb(path, structure, charmm_gui=True)
            names = [
                line[12:16].strip()
                for line in path.read_text().splitlines()
                if line.startswith("ATOM")
            ]
            self.assertEqual(names.count("HO62"), 4)
            self.assertFalse({"HO1", "HO2", "HO3", "HO4", "HO6"} & set(names))

    def test_validation_rejects_accidental_6_deoxy_residue(self):
        structure = build_structure(2, [1])
        residue = structure.chains[0].residues[1]
        residue.atoms = [
            atom for atom in residue.atoms if atom.name not in {"O6", "HO6"}
        ]
        residue.atoms.append(type(residue.atoms[0])("H63", "H", residue.atom("C6").position))
        with self.assertRaisesRegex(ValueError, "unexpected 6-deoxy H63"):
            validate(structure)

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

    def test_paajanen_18_chain_25_percent_of_eligible_sites(self):
        structure = build_structure(
            20,
            [2, 3, 4, 4, 3, 2],
            oxidation_degree=0.25,
            oxidation_scope="surface",
            oxidation_seed=344,
            oxidation_model="paajanen",
        )
        report = validate(structure)
        oxidized = [
            (chain.number, residue.number)
            for chain in structure.chains
            for residue in chain.residues
            if residue.oxidized
        ]
        self.assertEqual(structure.oxidation_eligible_sites, 120)
        self.assertEqual(report.oxidized_sites, 30)
        surface_chain_numbers = {1, 2, 3, 5, 6, 9, 10, 13, 14, 16, 17, 18}
        self.assertTrue(all(chain_number in surface_chain_numbers for chain_number, _ in oxidized))

    def test_paajanen_surface_25_percent_means_60_carboxylates(self):
        structure = build_structure(
            20,
            [2, 3, 4, 4, 3, 2],
            oxidation_degree=0.50,
            oxidation_scope="surface",
            oxidation_seed=344,
            oxidation_model="paajanen",
        )
        report = validate(structure)
        self.assertEqual(structure.oxidation_eligible_sites, 120)
        self.assertEqual(report.oxidized_sites, 60)
        self.assertEqual(report.net_charge, -60)
        self.assertEqual(report.oxidized_sites / 240, 0.25)

    def test_paajanen_selects_outward_alternating_class_per_surface_chain(self):
        structure = build_structure(20, [2, 3, 4, 4, 3, 2])
        eligible_ids = {id(residue) for residue in _paajanen_eligible_residues(structure.chains)}
        expected_parity = {
            1: 0, 2: 1, 3: 0, 5: 1, 6: 0, 9: 1,
            10: 0, 13: 1, 14: 0, 16: 1, 17: 0, 18: 1,
        }
        for chain in structure.chains:
            selected = [r for r in chain.residues if id(r) in eligible_ids]
            if chain.number in expected_parity:
                self.assertEqual(len(selected), 10)
                self.assertEqual(
                    {residue.number % 2 for residue in selected},
                    {expected_parity[chain.number]},
                )
            else:
                self.assertFalse(selected)

    def test_paajanen_requires_surface_scope(self):
        with self.assertRaisesRegex(ValueError, "requires oxidation_scope='surface'"):
            build_structure(
                20, [2, 3, 4, 4, 3, 2], oxidation_degree=0.25,
                oxidation_scope="all", oxidation_model="paajanen",
            )

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
