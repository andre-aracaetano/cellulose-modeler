from pathlib import Path
import argparse
import math

LAYERS = [1]
GLUCOSE_UNITS = 16
OUTPUT = "cellulose.pdb"

ATOM_ORDER = [
    "C1", "H1",
    "C5", "H5", "O5",
    "C2", "H2", "O2", "HO2",
    "C3", "H3", "O3", "HO3",
    "C4", "H4", "O4",
    "C6", "H61", "H62", "O6", "HO6"
]

FIRST_GLUCOSE = {
    "C1": (-0.355, 0.694, 1.227),
    "H1": (-0.068, 1.202, 2.176),
    "O1": (-1.581, 0.046, 1.461),
    "HO1": (-1.797, -0.496, 0.699),
    "C5": (0.784, 2.415, 0.033),
    "H5": (1.056, 2.874, 1.013),
    "O5": (-0.438, 1.681, 0.192),
    "C2": (0.710, -0.317, 0.847),
    "H2": (0.419, -0.796, -0.119),
    "O2": (0.812, -1.362, 1.808),
    "HO2": (0.756, -2.208, 1.357),
    "C3": (2.038, 0.387, 0.661),
    "H3": (2.335, 0.861, 1.626),
    "O3": (3.058, -0.537, 0.318),
    "HO3": (3.870, -0.285, 0.766),
    "C4": (1.889, 1.459, -0.390),
    "H4": (1.608, 1.007, -1.369),
    "O4": (3.055, 2.270, -0.554),
    "C6": (0.537, 3.520, -0.970),
    "H61": (1.460, 4.123, -1.111),
    "H62": (0.244, 3.093, -1.954),
    "O6": (-0.491, 4.400, -0.529),
    "HO6": (-1.337, 3.949, -0.593)
}

LINK_AB = {
    "C1": (-0.640326, 0.206353, -1.236193),
    "H1": (0.038968, -0.049458, -2.080527),
    "C5": (-1.758719, 1.826376, -2.591480),
    "H5": (-1.075145, 1.550141, -3.430021),
    "O5": (-1.088961, 1.563181, -1.350946),
    "C2": (-1.869572, -0.669743, -1.262797),
    "H2": (-2.535055, -0.380520, -0.414612),
    "O2": (-1.513685, -2.029284, -1.056477),
    "HO2": (-2.050870, -2.388266, -0.345311),
    "C3": (-2.620231, -0.477833, -2.562488),
    "H3": (-1.950590, -0.744948, -3.413824),
    "O3": (-3.749721, -1.332214, -2.636446),
    "HO3": (-3.774421, -1.742767, -3.504249),
    "C4": (-3.012710, 0.973418, -2.679453),
    "H4": (-3.694435, 1.259286, -1.843876),
    "O4": (-3.620784, 1.316081, -3.926533),
    "C6": (-2.038963, 3.310335, -2.659636),
    "H61": (-2.549751, 3.562266, -3.613847),
    "H62": (-2.698479, 3.620933, -1.819909),
    "O6": (-0.837217, 4.068987, -2.599025),
    "HO6": (-0.486636, 4.027977, -1.706032)
}

LINK_BA = {
    "C1": (-0.639688, 0.252726, -1.226760),
    "H1": (0.055902, 0.084780, -2.080424),
    "C5": (-1.832114, 1.899476, -2.481341),
    "H5": (-1.134174, 1.711077, -3.331667),
    "O5": (-1.152766, 1.590998, -1.256528),
    "C2": (-1.821542, -0.682682, -1.318484),
    "H2": (-2.501731, -0.484602, -0.455187),
    "O2": (-1.396510, -2.032942, -1.200135),
    "HO2": (-1.903021, -2.459267, -0.503777),
    "C3": (-2.577863, -0.442405, -2.605632),
    "H3": (-1.891070, -0.611686, -3.468774),
    "O3": (-3.657972, -1.349445, -2.754534),
    "HO3": (-3.552123, -1.829254, -3.578515),
    "C4": (-3.042476, 0.992317, -2.626761),
    "H4": (-3.731853, 1.187683, -1.771225),
    "O4": (-3.680994, 1.381606, -3.844407),
    "C6": (-2.187265, 3.369139, -2.457690),
    "H61": (-2.710114, 3.652499, -3.397437),
    "H62": (-2.861853, 3.594421, -1.603320),
    "O6": (-1.025418, 4.183545, -2.347656),
    "HO6": (-0.684362, 4.121249, -1.452950)
}

TERMINAL_HO4 = (-0.311682, -0.343048, -0.840799)


def add(a, b):
    return (
        a[0] + b[0],
        a[1] + b[1],
        a[2] + b[2]
    )


def sub(a, b):
    return (
        a[0] - b[0],
        a[1] - b[1],
        a[2] - b[2]
    )


def mul(v, scalar):
    return (
        v[0] * scalar,
        v[1] * scalar,
        v[2] * scalar
    )


def dot(a, b):
    return (
        a[0] * b[0] +
        a[1] * b[1] +
        a[2] * b[2]
    )


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    )


def norm(v):
    return math.sqrt(dot(v, v))


def normalize(v):
    length = norm(v)

    if length == 0:
        raise ValueError("Zero-length vector")

    return mul(v, 1.0 / length)


def distance(a, b):
    return norm(sub(a, b))


def make_frame(anchor, atom1, atom2):
    e1 = normalize(
        sub(atom1, anchor)
    )

    v2 = sub(
        atom2,
        anchor
    )

    v2 = sub(
        v2,
        mul(
            e1,
            dot(v2, e1)
        )
    )

    e2 = normalize(v2)

    e3 = normalize(
        cross(e1, e2)
    )

    return e1, e2, e3


def local_to_global(local, anchor, e1, e2, e3):
    return add(
        anchor,
        add(
            mul(e1, local[0]),
            add(
                mul(e2, local[1]),
                mul(e3, local[2])
            )
        )
    )


def build_next_residue(previous, link_template):
    anchor = previous["O4"]

    e1, e2, e3 = make_frame(
        anchor,
        previous["C4"],
        previous["C3"]
    )

    residue = {}

    for atom_name, local_position in link_template.items():
        residue[atom_name] = local_to_global(
            local_position,
            anchor,
            e1,
            e2,
            e3
        )

    return residue


def add_terminal_ho4(residue):
    anchor = residue["O4"]

    e1, e2, e3 = make_frame(
        anchor,
        residue["C4"],
        residue["C3"]
    )

    residue["HO4"] = local_to_global(
        TERMINAL_HO4,
        anchor,
        e1,
        e2,
        e3
    )


def build_chain(glucose_units):
    if glucose_units < 1:
        raise ValueError(
            "glucose_units must be at least 1"
        )

    first = {
        name: tuple(xyz)
        for name, xyz in FIRST_GLUCOSE.items()
    }

    residues = [first]

    for residue_index in range(1, glucose_units):
        previous = residues[-1]

        if residue_index % 2 == 1:
            link_template = LINK_AB
        else:
            link_template = LINK_BA

        next_residue = build_next_residue(
            previous,
            link_template
        )

        residues.append(
            next_residue
        )

    add_terminal_ho4(
        residues[-1]
    )

    return residues


def atom_order_for_residue(index, total):
    order = ATOM_ORDER.copy()

    if index == 0:
        position = order.index("H1") + 1

        order[position:position] = [
            "O1",
            "HO1"
        ]

    if index == total - 1:
        position = order.index("O4") + 1

        order.insert(
            position,
            "HO4"
        )

    return order


def element_from_atom_name(name):
    if name.startswith("H"):
        return "H"

    if name.startswith("O"):
        return "O"

    if name.startswith("C"):
        return "C"

    raise ValueError(
        f"Unknown atom element for {name}"
    )


def pdb_atom_name(name):
    if len(name) < 4:
        return f" {name:<3}"

    return name[:4]


def format_atom_line(
    serial,
    atom_name,
    resid,
    xyz
):
    element = element_from_atom_name(
        atom_name
    )

    return (
        f"ATOM  "
        f"{serial:5d} "
        f"{pdb_atom_name(atom_name)} "
        f"BGLC "
        f"{resid:4d}    "
        f"{xyz[0]:8.3f}"
        f"{xyz[1]:8.3f}"
        f"{xyz[2]:8.3f}"
        f"{1.00:6.2f}"
        f"{0.00:6.2f}"
        f"      "
        f"CARB"
        f"{element:>2s}"
    )


def write_pdb(path, residues):
    lines = [
        "REMARK  Generated by dynamic_cellulose",
        "REMARK  Single beta-1,4-glucan chain"
    ]

    serial = 1
    total = len(residues)

    for index, residue in enumerate(
        residues
    ):
        resid = index + 1

        order = atom_order_for_residue(
            index,
            total
        )

        for atom_name in order:
            xyz = residue[atom_name]

            lines.append(
                format_atom_line(
                    serial,
                    atom_name,
                    resid,
                    xyz
                )
            )

            serial += 1

    lines.append(
        f"TER   {serial:5d}"
        f"      BGLC "
        f"{total:4d}"
    )

    lines.append("END")

    Path(path).write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8"
    )


def validate_chain(residues):
    print()
    print("Validation")

    glycosidic_distances = []

    for i in range(
        len(residues) - 1
    ):
        d = distance(
            residues[i]["O4"],
            residues[i + 1]["C1"]
        )

        glycosidic_distances.append(
            d
        )

        print(
            f"BGLC {i + 1:2d} O4 -> "
            f"BGLC {i + 2:2d} C1 = "
            f"{d:.3f} A"
        )

    if glycosidic_distances:
        print()
        print(
            "Glycosidic O4-C1 range: "
            f"{min(glycosidic_distances):.3f} - "
            f"{max(glycosidic_distances):.3f} A"
        )

    first = residues[0]
    last = residues[-1]

    print(
        "Initial C1-O1: "
        f"{distance(first['C1'], first['O1']):.3f} A"
    )

    print(
        "Initial O1-HO1: "
        f"{distance(first['O1'], first['HO1']):.3f} A"
    )

    print(
        "Terminal O4-HO4: "
        f"{distance(last['O4'], last['HO4']):.3f} A"
    )


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--glucose",
        type=int,
        default=GLUCOSE_UNITS
    )

    parser.add_argument(
        "--output",
        default=OUTPUT
    )

    return parser.parse_args()


def main():
    args = parse_args()

    residues = build_chain(
        args.glucose
    )

    write_pdb(
        args.output,
        residues
    )

    print(
        f"Glucose units: {len(residues)}"
    )

    print(
        f"Atoms: "
        f"{sum(len(residue) for residue in residues)}"
    )

    print(
        f"Output: {Path(args.output).resolve()}"
    )

    validate_chain(
        residues
    )


if __name__ == "__main__":
    main()