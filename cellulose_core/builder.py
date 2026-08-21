from __future__ import annotations

import math
from statistics import mean

from .geometry import add, cross, dot, norm, normalize, scale, sub
from .model import Atom, Chain, Residue, Structure, Vec3
from .reference import REFERENCE


def _row_base(row: int) -> tuple[int, int, int]:
    half = row // 2
    if row % 2 == 0:
        return 0, half, -half
    return 1, half, -half - 1


def _lattice_sites(layers: list[int]) -> list[tuple[int, int, int, int, int]]:
    vector_b = REFERENCE.vector_b
    b_squared = dot(vector_b, vector_b)
    bases: list[tuple[int, int, int, Vec3]] = []
    ideal_centers = []
    for row, count in enumerate(layers):
        chain_type, cell_a, cell_b = _row_base(row)
        fractional_center = (
            cell_a + 0.5 * chain_type,
            cell_b + 0.5 * chain_type,
            0.0,
        )
        base = REFERENCE.fractional_to_cartesian(fractional_center)
        along_b = dot(base, vector_b) / b_squared
        ideal_centers.append(along_b + (count - 1) / 2.0)
        bases.append((chain_type, cell_a, cell_b, base))
    target_center = mean(ideal_centers)
    sites = []
    for row, (count, base_data) in enumerate(zip(layers, bases)):
        chain_type, cell_a, cell_b, base = base_data
        along_b = dot(base, vector_b) / b_squared
        first_column = round(target_center - along_b - (count - 1) / 2.0)
        for offset in range(count):
            column = first_column + offset
            sites.append((row, column, chain_type, cell_a, cell_b + column))
    return sites


def _build_chain(
    glucose_units: int,
    chain_type: int,
    cell_a: int,
    cell_b: int,
    number: int,
    row: int,
    column: int,
) -> Chain:
    residues = [
        Residue(index, REFERENCE.residue_atoms(chain_type, index, cell_a, cell_b))
        for index in range(1, glucose_units + 1)
    ]
    next_o4 = next(
        atom
        for atom in REFERENCE.residue_atoms(chain_type, glucose_units + 1, cell_a, cell_b)
        if atom.name == "O4"
    )
    residues[-1].atoms.append(Atom("O1", "O", next_o4.position))
    return Chain(number, row, column, residues)


def _center_transverse(chains: list[Chain]) -> None:
    anchors = [chain.residues[0].atom("C1").position for chain in chains]
    shift = (-mean(p[0] for p in anchors), -mean(p[1] for p in anchors), 0.0)
    for chain in chains:
        for residue in chain.residues:
            for atom in residue.atoms:
                atom.position = add(atom.position, shift)


def _hydroxyl_position(
    oxygen: Atom, carbon: Atom, occupied: list[Vec3], bond_length: float = 0.96
) -> Vec3:
    """Return a tetrahedral O-H position with minimal steric overlap.

    Rotation around C-O is not determined by the experimental CIF. Twelve
    valid rotamers are sampled while every experimental coordinate is kept.
    """
    axis = normalize(sub(carbon.position, oxygen.position))
    coordinate_axes: tuple[Vec3, ...] = (
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
    )
    seed = min(coordinate_axes, key=lambda candidate: abs(dot(axis, candidate)))
    perpendicular = normalize(cross(axis, seed))
    second = normalize(cross(axis, perpendicular))
    angle = math.radians(109.47)
    candidates = []
    for step in range(12):
        torsion = math.radians(step * 30.0)
        radial = add(
            scale(perpendicular, math.cos(torsion)),
            scale(second, math.sin(torsion)),
        )
        direction = add(scale(axis, math.cos(angle)), scale(radial, math.sin(angle)))
        candidates.append(add(oxygen.position, scale(direction, bond_length)))

    excluded = {oxygen.position, carbon.position}
    # Only atoms close to O can constrain a 0.96 A proton. Prefiltering keeps
    # construction effectively linear even for large crystallites.
    nearby = [
        other
        for other in occupied
        if other not in excluded and norm(sub(oxygen.position, other)) < 4.0
    ]

    def clearance(position: Vec3) -> float:
        separations = [norm(sub(position, other)) for other in nearby]
        return min(separations, default=float("inf"))

    return max(candidates, key=clearance)


def _add_missing_hydroxyl_hydrogens(chains: list[Chain]) -> None:
    """Complete every hydroxyl, including both ends of each finite chain."""
    occupied = [
        atom.position
        for chain in chains
        for residue in chain.residues
        for atom in residue.atoms
    ]
    for chain in chains:
        last = len(chain.residues)
        for residue in chain.residues:
            hydroxyls = [
                ("O2", "C2", "HO2"),
                ("O3", "C3", "HO3"),
                ("O6", "C6", "HO6"),
            ]
            if residue.number == 1:
                hydroxyls.append(("O4", "C4", "HO4"))
            if residue.number == last:
                hydroxyls.append(("O1", "C1", "HO1"))
            for oxygen_name, carbon_name, hydrogen_name in hydroxyls:
                oxygen = residue.atom(oxygen_name)
                carbon = residue.atom(carbon_name)
                position = _hydroxyl_position(oxygen, carbon, occupied)
                residue.atoms.append(Atom(hydrogen_name, "H", position))
                occupied.append(position)


def build_structure(glucose_units: int, layers: list[int]) -> Structure:
    if glucose_units < 1:
        raise ValueError("glucose_units must be at least 1")
    if not layers or any(count < 1 for count in layers):
        raise ValueError("layers must contain positive chain counts")
    chains = [
        _build_chain(glucose_units, chain_type, cell_a, cell_b, index, row, column)
        for index, (row, column, chain_type, cell_a, cell_b) in enumerate(
            _lattice_sites(layers), start=1
        )
    ]
    _add_missing_hydroxyl_hydrogens(chains)
    _center_transverse(chains)
    return Structure("I-beta", list(layers), glucose_units, chains)
