from __future__ import annotations

from dataclasses import dataclass, field


Vec3 = tuple[float, float, float]


@dataclass
class Atom:
    name: str
    element: str
    position: Vec3


@dataclass
class Residue:
    number: int
    atoms: list[Atom] = field(default_factory=list)
    oxidized: bool = False
    oxidation_protonated: bool = False

    def atom(self, name: str) -> Atom:
        return next(atom for atom in self.atoms if atom.name == name)


@dataclass
class Chain:
    number: int
    row: int
    column: int
    residues: list[Residue]


@dataclass
class Structure:
    allomorph: str
    layers: list[int]
    glucose_units: int
    chains: list[Chain]
    oxidation_scope: str | None = None
    oxidation_degree: float = 0.0
    oxidation_protonated: bool = False
    oxidation_seed: int | None = None
