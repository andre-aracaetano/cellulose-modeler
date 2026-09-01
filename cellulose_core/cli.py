from __future__ import annotations

import argparse
from pathlib import Path

from .builder import build_structure
from .pdbio import write_pdb
from .validation import validate


PRESETS = {
    "single": [1],
    "18": [2, 3, 4, 4, 3, 2],
    "36": [3, 4, 5, 6, 6, 5, 4, 3],
}


def parse_layers(text: str) -> list[int]:
    try:
        layers = [int(value.strip()) for value in text.split(",")]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("use comma-separated integers, e.g. 2,3,4,4,3,2") from exc
    if not layers or any(value < 1 for value in layers):
        raise argparse.ArgumentTypeError("all layer sizes must be positive")
    return layers


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        prog="dynamics-cellulose",
        description="Generate cellulose I-beta chains and configurable crystallites as PDB.",
    )
    layout = result.add_mutually_exclusive_group()
    layout.add_argument(
        "--layers",
        type=parse_layers,
        metavar="N,N,...",
        help="chains in each crystal row from top to bottom, e.g. 2,3,4,4,3,2",
    )
    layout.add_argument(
        "--preset",
        choices=sorted(PRESETS),
        help="predefined transverse geometry (default: single)",
    )
    result.add_argument(
        "--glucose",
        "--dp",
        type=int,
        default=20,
        dest="glucose",
        metavar="N",
        help="glucose units in every cellulose chain (default: 20)",
    )
    result.add_argument(
        "--output",
        "-o",
        type=Path,
        default=Path("outputs/cellulose.pdb"),
        help="output PDB path (default: outputs/cellulose.pdb)",
    )
    result.add_argument(
        "--no-conect",
        action="store_true",
        help="omit PDB CONECT records (not recommended for carbohydrate readers)",
    )
    result.add_argument(
        "--charmm-gui",
        action="store_true",
        help=(
            "prepare CHARMM-GUI input: omit hydroxyl H coordinates and order "
            "residues for O4(i)-C1(i+1) recognition"
        ),
    )
    result.add_argument(
        "--oxidation",
        type=float,
        metavar="FRACTION",
        help="random C6 oxidation degree from 0.1 to 1.0 (disabled by default)",
    )
    result.add_argument(
        "--oxidation-scope",
        choices=("surface", "all"),
        default="surface",
        help="oxidize boundary chains only or all chains (default: surface)",
    )
    result.add_argument(
        "--oxidation-model",
        choices=("random-surface", "paajanen"),
        default="random-surface",
        help=(
            "eligible-site model: all C6 in the requested scope, or alternating "
            "C6 sites on surface chains following Paajanen et al. (default: random-surface)"
        ),
    )
    result.add_argument(
        "--oxidation-state",
        choices=("deprotonated", "protonated"),
        default="deprotonated",
        help="generate COO- or COOH groups (default: deprotonated)",
    )
    result.add_argument(
        "--seed",
        type=int,
        help="random seed for reproducible oxidation-site selection",
    )
    result.add_argument(
        "--interactive", "-i", action="store_true", help="ask for DP, layers, and output path"
    )
    return result


def main() -> None:
    args = parser().parse_args()
    if args.interactive:
        args.glucose = int(input("Glucose units per chain (DP): ").strip())
        layout_text = input(
            "Layers top-to-bottom (e.g. 2,3,4,4,3,2; blank for one chain): "
        ).strip()
        args.layers = parse_layers(layout_text) if layout_text else [1]
        output_text = input(f"Output PDB [{args.output}]: ").strip()
        if output_text:
            args.output = Path(output_text)
    if args.glucose < 1:
        raise SystemExit("error: --glucose/--dp must be at least 1")
    if args.oxidation is not None and not 0.1 <= args.oxidation <= 1.0:
        raise SystemExit("error: --oxidation must be between 0.1 and 1.0")
    if args.oxidation_model == "paajanen" and args.oxidation_scope != "surface":
        raise SystemExit("error: --oxidation-model paajanen requires --oxidation-scope surface")
    layers = args.layers if args.layers is not None else PRESETS[args.preset or "single"]
    structure = build_structure(
        args.glucose,
        layers,
        oxidation_degree=args.oxidation or 0.0,
        oxidation_scope=args.oxidation_scope,
        oxidation_protonated=args.oxidation_state == "protonated",
        oxidation_seed=args.seed,
        oxidation_model=args.oxidation_model,
    )
    report = validate(structure)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_pdb(
        args.output,
        structure,
        conect=not args.no_conect,
        charmm_gui=args.charmm_gui,
    )
    print(f"Allomorph: cellulose {structure.allomorph}")
    print(f"Layers: {','.join(map(str, structure.layers))}")
    print(f"Chains: {report.chains}")
    print(f"Glucose units per chain: {structure.glucose_units}")
    print(f"Residues: {report.residues}")
    print(f"Atoms: {report.atoms}")
    if report.oxidized_sites:
        print(f"Oxidation model: {structure.oxidation_model}")
        print(f"Eligible C6 sites: {structure.oxidation_eligible_sites}")
        print(f"Oxidized C6 sites: {report.oxidized_sites}")
        print(f"Nominal net charge: {report.net_charge:+d} e")
    if report.glycosidic_min is not None:
        print(
            f"Glycosidic C1-O4 range: {report.glycosidic_min:.3f} - "
            f"{report.glycosidic_max:.3f} A"
        )
    if report.chain_spacing_min is not None:
        print(f"Minimum chain-axis spacing: {report.chain_spacing_min:.3f} A")
    if args.charmm_gui:
        print("Output profile: CHARMM-GUI (hydroxyl H coordinates omitted)")
    print(f"Output: {args.output.resolve()}")
