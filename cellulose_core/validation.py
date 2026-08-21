from __future__ import annotations

from dataclasses import dataclass

from .geometry import distance
from .model import Structure


@dataclass(frozen=True)
class Validation:
    chains: int
    residues: int
    atoms: int
    glycosidic_min: float | None
    glycosidic_max: float | None
    chain_spacing_min: float | None


def validate(structure: Structure) -> Validation:
    expected_residues = len(structure.chains) * structure.glucose_units
    residues = sum(len(chain.residues) for chain in structure.chains)
    if residues != expected_residues:
        raise ValueError(f"residue count mismatch: expected {expected_residues}, found {residues}")
    distances = []
    for chain in structure.chains:
        for index in range(len(chain.residues) - 1):
            distances.append(
                distance(
                    chain.residues[index].atom("C1").position,
                    chain.residues[index + 1].atom("O4").position,
                )
            )
    if distances and not all(1.30 <= value <= 1.55 for value in distances):
        raise ValueError("invalid glycosidic C1-O4 distance detected")
    anchors = [chain.residues[0].atom("C1").position for chain in structure.chains]
    spacings = [
        distance(first, second)
        for index, first in enumerate(anchors)
        for second in anchors[index + 1 :]
    ]
    if spacings and min(spacings) < 5.0:
        raise ValueError("two chains occupy overlapping or invalid lattice sites")
    return Validation(
        chains=len(structure.chains),
        residues=residues,
        atoms=sum(len(residue.atoms) for chain in structure.chains for residue in chain.residues),
        glycosidic_min=min(distances) if distances else None,
        glycosidic_max=max(distances) if distances else None,
        chain_spacing_min=min(spacings) if spacings else None,
    )
