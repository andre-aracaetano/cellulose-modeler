"""Experimental cellulose I-beta crystallographic reference.

Coordinates come from COD 4114994 (Nishiyama, Langan & Chanzy, 2002),
which the COD file declares to be in the public domain. No coordinates from
Cellulose Builder are used here.
"""

from __future__ import annotations

import math
import re
from importlib.resources import files

from .model import Atom, Vec3


LABEL_TO_PDB = {
    "C1": "C1", "C2": "C2", "C3": "C3", "C4": "C4", "C5": "C5", "C6": "C6",
    "H1": "H1", "H2": "H2", "H3": "H3", "H4": "H4", "H5": "H5",
    "H6A": "H61", "H6B": "H62",
    "O1": "O4", "O2": "O2", "O3": "O3", "O5": "O5", "O6": "O6",
}


def _number(value: str) -> float:
    """Convert a CIF number, discarding its parenthesized uncertainty."""
    return float(re.sub(r"\([^)]*\)$", "", value))


class CrystallographicReference:
    def __init__(self) -> None:
        path = files("cellulose_core").joinpath("data/COD_4114994_cellulose_Ibeta.cif")
        text = path.read_text(encoding="utf-8")
        self.a = self._tag(text, "_cell_length_a")
        self.b = self._tag(text, "_cell_length_b")
        self.c = self._tag(text, "_cell_length_c")
        self.alpha = self._tag(text, "_cell_angle_alpha")
        self.beta = self._tag(text, "_cell_angle_beta")
        self.gamma = self._tag(text, "_cell_angle_gamma")
        self.fractional_atoms = self._read_atom_sites(text)
        gamma = math.radians(self.gamma)
        self.vector_a: Vec3 = (self.a, 0.0, 0.0)
        self.vector_b: Vec3 = (self.b * math.cos(gamma), self.b * math.sin(gamma), 0.0)
        self.vector_c: Vec3 = (0.0, 0.0, self.c)

    @staticmethod
    def _tag(text: str, name: str) -> float:
        match = re.search(rf"^{re.escape(name)}\s+(\S+)", text, re.MULTILINE)
        if not match:
            raise ValueError(f"required CIF tag missing: {name}")
        return _number(match.group(1))

    @staticmethod
    def _read_atom_sites(text: str) -> dict[int, list[tuple[str, str, Vec3]]]:
        marker = "_atom_site_refinement_flags\n"
        try:
            block = text.split(marker, 1)[1].split("\nloop_", 1)[0]
        except IndexError as exc:
            raise ValueError("COD CIF atom-site loop was not found") from exc
        result: dict[int, list[tuple[str, str, Vec3]]] = {0: [], 1: []}
        for line in block.splitlines():
            fields = line.split()
            if len(fields) < 5:
                continue
            label, element = fields[0], fields[1]
            label_match = re.fullmatch(r"([CHO])([1-6])([12])([AB]?)", label)
            if not label_match:
                continue
            atom_element, atom_number, chain_digit, suffix = label_match.groups()
            chain_type = int(chain_digit) - 1
            base_label = f"{atom_element}{atom_number}{suffix}"
            if base_label not in LABEL_TO_PDB:
                continue
            result[chain_type].append(
                (
                    LABEL_TO_PDB[base_label],
                    element,
                    tuple(_number(value) for value in fields[2:5]),  # type: ignore[arg-type]
                )
            )
        if not all(result.values()):
            raise ValueError("both crystallographically independent chains are required")
        return result

    def fractional_to_cartesian(self, fractional: Vec3) -> Vec3:
        return tuple(
            fractional[0] * self.vector_a[index]
            + fractional[1] * self.vector_b[index]
            + fractional[2] * self.vector_c[index]
            for index in range(3)
        )  # type: ignore[return-value]

    def residue_atoms(
        self, chain_type: int, residue_number: int, cell_a: int, cell_b: int
    ) -> list[Atom]:
        identity = residue_number % 2 == 1
        cycle = (residue_number - 1) // 2
        target_center = 0.0 if chain_type == 0 else 0.5
        atoms = []
        for name, element, fractional in self.fractional_atoms[chain_type]:
            if identity:
                x, y, z = fractional
            else:
                x, y, z = -fractional[0], -fractional[1], fractional[2] + 0.5
            x += round(target_center - x) + cell_a
            y += round(target_center - y) + cell_b
            z += cycle
            atoms.append(Atom(name, element, self.fractional_to_cartesian((x, y, z))))
        return atoms


REFERENCE = CrystallographicReference()
