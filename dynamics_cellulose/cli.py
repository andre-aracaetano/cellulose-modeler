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
    layout.add_argument("--layers", type=parse_layers, help="chains per row, top-to-bottom")
    layout.add_argument("--preset", choices=sorted(PRESETS))
    result.add_argument("--glucose", "--dp", type=int, default=20, dest="glucose")
    result.add_argument("--output", "-o", type=Path, default=Path("cellulose.pdb"))
    result.add_argument("--no-conect", action="store_true", help="omit explicit PDB CONECT records")
    result.add_argument("--interactive", "-i", action="store_true", help="ask for specifications")
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
    layers = args.layers if args.layers is not None else PRESETS[args.preset or "single"]
    structure = build_structure(args.glucose, layers)
    report = validate(structure)
    write_pdb(args.output, structure, conect=not args.no_conect)
    print(f"Allomorph: cellulose {structure.allomorph}")
    print(f"Layers: {','.join(map(str, structure.layers))}")
    print(f"Chains: {report.chains}")
    print(f"Glucose units per chain: {structure.glucose_units}")
    print(f"Residues: {report.residues}")
    print(f"Atoms: {report.atoms}")
    if report.glycosidic_min is not None:
        print(
            f"Glycosidic C1-O4 range: {report.glycosidic_min:.3f} - "
            f"{report.glycosidic_max:.3f} A"
        )
    if report.chain_spacing_min is not None:
        print(f"Minimum chain-axis spacing: {report.chain_spacing_min:.3f} A")
    print(f"Output: {args.output.resolve()}")
