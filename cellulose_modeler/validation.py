from __future__ import annotations

from dataclasses import dataclass

from .geometry import distance
from .model import Structure


NORMAL_GLUCOSE_ELEMENTS = {
    "C1": "C", "C2": "C", "C3": "C", "C4": "C", "C5": "C", "C6": "C",
    "H1": "H", "H2": "H", "H3": "H", "H4": "H", "H5": "H",
    "H61": "H", "H62": "H",
    "O2": "O", "O3": "O", "O4": "O", "O5": "O", "O6": "O",
    "HO2": "H", "HO3": "H", "HO6": "H",
}


@dataclass(frozen=True)
class Validation:
    chains: int
    residues: int
    atoms: int
    glycosidic_min: float | None
    glycosidic_max: float | None
    chain_spacing_min: float | None
    oxidized_sites: int
    net_charge: int


def validate(structure: Structure) -> Validation:
    expected_residues = len(structure.chains) * structure.glucose_units
    residues = sum(len(chain.residues) for chain in structure.chains)
    if residues != expected_residues:
        raise ValueError(f"residue count mismatch: expected {expected_residues}, found {residues}")
    distances = []
    for chain in structure.chains:
        last_residue = len(chain.residues)
        for residue in chain.residues:
            names = [atom.name for atom in residue.atoms]
            if len(names) != len(set(names)):
                raise ValueError(
                    f"duplicate atom name in chain {chain.number}, residue {residue.number}"
                )
            atoms = {atom.name: atom for atom in residue.atoms}
            if "H63" in atoms:
                raise ValueError(
                    f"unexpected 6-deoxy H63 in chain {chain.number}, residue {residue.number}"
                )
            if not residue.oxidized:
                missing = set(NORMAL_GLUCOSE_ELEMENTS) - atoms.keys()
                if missing:
                    raise ValueError(
                        f"incomplete glucose in chain {chain.number}, residue "
                        f"{residue.number}: missing {', '.join(sorted(missing))}"
                    )
                for name, element in NORMAL_GLUCOSE_ELEMENTS.items():
                    if atoms[name].element != element:
                        raise ValueError(
                            f"wrong element for {name} in chain {chain.number}, "
                            f"residue {residue.number}"
                        )
                for first, second, lower, upper in (
                    ("C5", "C6", 1.40, 1.65),
                    ("C6", "O6", 1.30, 1.55),
                    ("C6", "H61", 0.95, 1.20),
                    ("C6", "H62", 0.95, 1.20),
                    ("O6", "HO6", 0.90, 1.05),
                ):
                    value = distance(atoms[first].position, atoms[second].position)
                    if not lower <= value <= upper:
                        raise ValueError(
                            f"implausible {first}-{second} distance in chain "
                            f"{chain.number}, residue {residue.number}: {value:.3f} A"
                        )
            if residue.number == 1 and "HO4" not in atoms:
                raise ValueError(f"missing HO4 at the start of chain {chain.number}")
            if residue.number == last_residue and not {"O1", "HO1"} <= atoms.keys():
                raise ValueError(f"incomplete reducing end of chain {chain.number}")
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
        distance((first[0], first[1], 0.0), (second[0], second[1], 0.0))
        for index, first in enumerate(anchors)
        for second in anchors[index + 1 :]
    ]
    # Cellulose III_I has a legitimate crystallographic a spacing of 4.45 A.
    # A 4.0 A guard still detects duplicate/overlapping chain lattice sites.
    if spacings and min(spacings) < 4.0:
        raise ValueError("two chains occupy overlapping or invalid lattice sites")
    return Validation(
        chains=len(structure.chains),
        residues=residues,
        atoms=sum(len(residue.atoms) for chain in structure.chains for residue in chain.residues),
        glycosidic_min=min(distances) if distances else None,
        glycosidic_max=max(distances) if distances else None,
        chain_spacing_min=min(spacings) if spacings else None,
        oxidized_sites=sum(
            residue.oxidized for chain in structure.chains for residue in chain.residues
        ),
        net_charge=-sum(
            residue.oxidized and not residue.oxidation_protonated
            for chain in structure.chains for residue in chain.residues
        ),
    )
