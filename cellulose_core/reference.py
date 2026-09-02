"""Independent experimental crystallographic references for cellulose.

The bundled data are taken from public COD deposits (I-alpha and I-beta) or
transcribed from the primary crystallographic publications (II and III_I).
No file from Cellulose Builder is read, imported, or required.
"""

from __future__ import annotations

import itertools
import math
import re
from dataclasses import dataclass
from importlib.resources import files

from .geometry import cross, dot, norm, normalize, scale, sub
from .model import Atom, Vec3


PDB_NAMES = {
    "C1": "C1", "C2": "C2", "C3": "C3", "C4": "C4", "C5": "C5", "C6": "C6",
    "H1": "H1", "H2": "H2", "H3": "H3", "H4": "H4", "H5": "H5",
    "H6A": "H61", "H6B": "H62", "H61": "H61", "H62": "H62",
    "O1": "O4", "O4": "O4", "O2": "O2", "O3": "O3", "O5": "O5", "O6": "O6",
}


def _number(value: str) -> float:
    return float(re.sub(r"\([^)]*\)$", "", value))


def _cell_vectors(
    a: float, b: float, c: float, alpha: float, beta: float, gamma: float
) -> tuple[Vec3, Vec3, Vec3]:
    alpha_r, beta_r, gamma_r = map(math.radians, (alpha, beta, gamma))
    vector_a: Vec3 = (a, 0.0, 0.0)
    vector_b: Vec3 = (b * math.cos(gamma_r), b * math.sin(gamma_r), 0.0)
    cx = c * math.cos(beta_r)
    cy = c * (math.cos(alpha_r) - math.cos(beta_r) * math.cos(gamma_r)) / math.sin(gamma_r)
    cz = math.sqrt(max(0.0, c * c - cx * cx - cy * cy))
    return vector_a, vector_b, (cx, cy, cz)


@dataclass(frozen=True)
class ReferenceMetadata:
    key: str
    name: str
    filename: str
    doi: str
    citation: str
    data_source: str
    chain_axis: int
    parser: str


METADATA = {
    "ibeta": ReferenceMetadata(
        "ibeta", "I-beta", "COD_4114994_cellulose_Ibeta.cif", "10.1021/ja0257319",
        "Nishiyama, Langan & Chanzy; JACS 2002, 124, 9074-9082",
        "Crystallography Open Database 4114994 (CC0)", 2, "cod-ibeta",
    ),
    "ialpha": ReferenceMetadata(
        "ialpha", "I-alpha", "COD_4114383_cellulose_Ialpha.cif", "10.1021/ja037055w",
        "Nishiyama, Sugiyama, Chanzy & Langan; JACS 2003, 125, 14300-14306",
        "Crystallography Open Database 4114383 (CC0)", 0, "cod-ialpha",
    ),
    "ii": ReferenceMetadata(
        "ii", "II", "cellulose_II_Langan2001.cif", "10.1021/bm005612q",
        "Langan, Nishiyama & Chanzy; Biomacromolecules 2001, 2, 410-416",
        "Primary-paper crystallographic coordinates transcribed into CIF", 2, "auth-chain",
    ),
    "iii-i": ReferenceMetadata(
        "iii-i", "III_I", "cellulose_IIII_Wada2004.cif", "10.1021/ma0485585",
        "Wada, Chanzy, Nishiyama & Langan; Macromolecules 2004, 37, 8548-8555",
        "Primary-paper crystallographic coordinates transcribed into CIF", 2, "auth-chain",
    ),
}

FractionalAtom = tuple[str, str, Vec3]


class CrystallographicReference:
    """A cellulose asymmetric unit expanded into covalently ordered chains."""

    def __init__(self, metadata: ReferenceMetadata) -> None:
        self.metadata = metadata
        path = files("cellulose_core").joinpath(f"data/{metadata.filename}")
        text = path.read_text(encoding="utf-8")
        self.a = self._tag(text, "_cell_length_a")
        self.b = self._tag(text, "_cell_length_b")
        self.c = self._tag(text, "_cell_length_c")
        self.alpha = self._tag(text, "_cell_angle_alpha")
        self.beta = self._tag(text, "_cell_angle_beta")
        self.gamma = self._tag(text, "_cell_angle_gamma")
        self.crystal_vectors = _cell_vectors(
            self.a, self.b, self.c, self.alpha, self.beta, self.gamma
        )
        self._frame = self._orientation_frame(metadata.chain_axis)
        self.vector_a, self.vector_b, self.vector_c = tuple(
            self._orient(vector) for vector in self.crystal_vectors
        )

        if metadata.parser == "cod-ibeta":
            raw_templates = [[atoms] for atoms in self._read_cod_ibeta(text)]
        elif metadata.parser == "cod-ialpha":
            raw_templates = [self._read_cod_ialpha(text)]
        else:
            raw_templates = [[atoms] for atoms in self._read_auth_chain(text)]

        self.templates: list[tuple[list[FractionalAtom], list[FractionalAtom]]] = []
        self.cycle_translations: list[tuple[int, int, int]] = []
        for templates in raw_templates:
            first = templates[0]
            if len(templates) == 1:
                second_unshifted = [
                    (name, element, (-position[0], -position[1], position[2] + 0.5))
                    for name, element, position in first
                ]
            else:
                second_unshifted = templates[1]
            second_shift = self._best_link_translation(first, second_unshifted)
            second = self._translate_atoms(second_unshifted, second_shift)
            cycle = self._best_link_translation(second, first)
            self.templates.append((first, second))
            self.cycle_translations.append(cycle)

        # Backward-compatible view used by the original I-beta validation.
        self.fractional_atoms = {
            index: templates[0] for index, templates in enumerate(self.templates)
        }
        self.transverse_axes = tuple(index for index in range(3) if index != metadata.chain_axis)

    @property
    def key(self) -> str:
        return self.metadata.key

    @property
    def name(self) -> str:
        return self.metadata.name

    @property
    def chain_types(self) -> int:
        return len(self.templates)

    @staticmethod
    def _tag(text: str, name: str) -> float:
        match = re.search(rf"^{re.escape(name)}\s+(\S+)", text, re.MULTILINE)
        if not match:
            raise ValueError(f"required CIF tag missing: {name}")
        return _number(match.group(1))

    def _orientation_frame(self, chain_axis: int) -> tuple[Vec3, Vec3, Vec3]:
        chain = normalize(self.crystal_vectors[chain_axis])
        transverse_index = next(index for index in range(3) if index != chain_axis)
        candidate = self.crystal_vectors[transverse_index]
        x_axis = normalize(sub(candidate, scale(chain, dot(candidate, chain))))
        y_axis = normalize(cross(chain, x_axis))
        return x_axis, y_axis, chain

    def _orient(self, vector: Vec3) -> Vec3:
        return tuple(dot(vector, axis) for axis in self._frame)  # type: ignore[return-value]

    def fractional_to_cartesian(self, fractional: Vec3) -> Vec3:
        crystal = tuple(
            sum(fractional[j] * self.crystal_vectors[j][i] for j in range(3))
            for i in range(3)
        )
        return self._orient(crystal)  # type: ignore[arg-type]

    @staticmethod
    def _atom_block(text: str) -> list[list[str]]:
        marker = "_atom_site_refinement_flags\n"
        try:
            block = text.split(marker, 1)[1].split("\nloop_", 1)[0]
        except IndexError as exc:
            raise ValueError("COD CIF atom-site loop was not found") from exc
        return [line.split() for line in block.splitlines() if line.strip()]

    def _read_cod_ibeta(self, text: str) -> list[list[FractionalAtom]]:
        result: list[list[FractionalAtom]] = [[], []]
        for fields in self._atom_block(text):
            if len(fields) < 5:
                continue
            match = re.fullmatch(r"([CHO])([1-6])([12])([AB]?)", fields[0])
            if not match:
                continue
            atom_element, atom_number, chain_digit, suffix = match.groups()
            base = f"{atom_element}{atom_number}{suffix}"
            if base in PDB_NAMES:
                result[int(chain_digit) - 1].append(
                    (PDB_NAMES[base], fields[1], tuple(_number(v) for v in fields[2:5]))  # type: ignore[arg-type]
                )
        if not all(result):
            raise ValueError("I-beta CIF must contain both independent chains")
        return result

    def _read_cod_ialpha(self, text: str) -> list[list[FractionalAtom]]:
        result: list[list[FractionalAtom]] = [[], []]
        for fields in self._atom_block(text):
            if len(fields) < 5 or fields[1] == "D":
                continue
            match = re.fullmatch(r"([CHO])([1-6])([12])([AB]?)", fields[0])
            if not match:
                continue
            atom_element, atom_number, residue_digit, suffix = match.groups()
            base = f"{atom_element}{atom_number}{suffix}"
            if base in PDB_NAMES:
                # In the P1 I-alpha deposit O11/O12 are the bridging oxygens
                # attached to the opposite glucosyl residue in the covalent
                # chain sequence, so their template assignment is exchanged.
                template = int(residue_digit) - 1
                position = tuple(_number(v) for v in fields[2:5])
                if atom_element == "O" and atom_number == "1":
                    template = 1 - template
                    if residue_digit == "2":
                        position = (position[0] - 1.0, position[1], position[2])
                result[template].append(
                    (PDB_NAMES[base], fields[1], position)  # type: ignore[arg-type]
                )
        if not all(result):
            raise ValueError("I-alpha CIF must contain both glucosyl residues")
        return result

    @staticmethod
    def _read_auth_chain(text: str) -> list[list[FractionalAtom]]:
        try:
            block = text.split("_atom_site_fract_z\n", 1)[1]
        except IndexError as exc:
            raise ValueError("custom CIF atom-site loop was not found") from exc
        groups: dict[str, list[FractionalAtom]] = {}
        for line in block.splitlines():
            fields = line.split()
            if len(fields) != 6:
                continue
            name, element, chain = fields[:3]
            if name in PDB_NAMES:
                groups.setdefault(chain, []).append(
                    (PDB_NAMES[name], element, tuple(_number(v) for v in fields[3:6]))  # type: ignore[arg-type]
                )
        if not groups:
            raise ValueError("custom CIF contains no recognized cellulose atoms")
        return list(groups.values())

    @staticmethod
    def _translate_atoms(
        atoms: list[FractionalAtom], translation: tuple[int, int, int]
    ) -> list[FractionalAtom]:
        return [
            (name, element, tuple(position[i] + translation[i] for i in range(3)))  # type: ignore[arg-type]
            for name, element, position in atoms
        ]

    @staticmethod
    def _position(atoms: list[FractionalAtom], name: str) -> Vec3:
        return next(position for atom_name, _, position in atoms if atom_name == name)

    def _best_link_translation(
        self, donor: list[FractionalAtom], acceptor: list[FractionalAtom]
    ) -> tuple[int, int, int]:
        c1 = self.fractional_to_cartesian(self._position(donor, "C1"))
        o4 = self._position(acceptor, "O4")
        candidates = itertools.product(range(-2, 3), repeat=3)
        return min(
            candidates,
            key=lambda shift: norm(sub(c1, self.fractional_to_cartesian(
                tuple(o4[i] + shift[i] for i in range(3))  # type: ignore[arg-type]
            ))),
        )

    def site_translation(self, row: int, column: int) -> tuple[int, int, int]:
        result = [0, 0, 0]
        result[self.transverse_axes[0]] = row
        result[self.transverse_axes[1]] = column
        return tuple(result)  # type: ignore[return-value]

    def residue_atoms(
        self,
        chain_type: int,
        residue_number: int,
        site_translation: tuple[int, int, int],
        glucose_units: int,
    ) -> list[Atom]:
        template_index = (residue_number - 1) % 2
        cycle_index = (residue_number - 1) // 2
        cycle = self.cycle_translations[chain_type]
        cycles = (glucose_units + 1) // 2
        start = tuple(-cycles * min(component, 0) for component in cycle)
        shift = tuple(
            site_translation[i] + start[i] + cycle_index * cycle[i] for i in range(3)
        )
        return [
            Atom(
                name,
                element,
                self.fractional_to_cartesian(
                    tuple(fractional[i] + shift[i] for i in range(3))  # type: ignore[arg-type]
                ),
            )
            for name, element, fractional in self.templates[chain_type][template_index]
        ]


REFERENCES = {key: CrystallographicReference(metadata) for key, metadata in METADATA.items()}
REFERENCE = REFERENCES["ibeta"]


def get_reference(allomorph: str) -> CrystallographicReference:
    try:
        return REFERENCES[allomorph]
    except KeyError as exc:
        choices = ", ".join(REFERENCES)
        raise ValueError(f"unknown cellulose allomorph {allomorph!r}; choose {choices}") from exc
