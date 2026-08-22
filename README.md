# Dynamics Cellulose Maker

`dynamics_cellulose` is a dependency-free Python generator for finite cellulose
Iβ chains and crystallites with arbitrary degrees of polymerization and
transverse row profiles. The user-facing launcher is `builder.py`; reusable
implementation modules live in `cellulose_core/`, and generated structures go
to `outputs/` by default.

Version 0.3 adds optional, reproducible C6 oxidation to protonated carboxylic
acid (`COOH`) or deprotonated carboxylate (`COO-`) groups. Oxidation can be
restricted to surface chains or applied across the complete crystallite.

## Scientific basis and provenance

The program constructs cellulose Iβ directly from the experimental fractional
coordinates deposited as [Crystallography Open Database entry
4114994](https://www.crystallography.net/cod/4114994.html). The bundled CIF
declares its contents to be in the public domain. Generation never reads
coordinates produced by Cellulose Builder or CHARMM-GUI; legacy PDBs retained
under `examples/` serve only as output-comparison fixtures.

Experimental source:

> Yoshiharu Nishiyama, Paul Langan, and Henri Chanzy. “Crystal Structure and
> Hydrogen-Bonding System in Cellulose Iβ from Synchrotron X-ray and Neutron
> Fiber Diffraction.” *Journal of the American Chemical Society* 124 (2002),
> 9074–9082. <https://doi.org/10.1021/ja0257319>.

The implementation independently parses the deposited CIF, applies its
`P 1 1 21` symmetry, converts fractional coordinates through the experimental
unit cell, repeats the crystallographic motif to the requested DP, places the
two independent chain types on discrete Iβ sites, and writes the selected row
profile as PDB.

See [the detailed data provenance](cellulose_core/data/PROVENANCE.md) and
[`CITATION.cff`](CITATION.cff).

## Quick start

Download the project, open a terminal in its root directory, and run:

```bash
python builder.py --dp 20 --preset 18
```

The resulting file is written to `outputs/cellulose.pdb`. The directory is
created automatically when needed. On systems where Python 3 is exposed as
`python3` (commonly Linux and macOS), use `python3 builder.py` instead.

An arbitrary transverse arrangement is a comma-separated list of chain counts
from top to bottom:

```bash
python builder.py \
  --dp 16 \
  --layers 3,5,7,7,5,3 \
  --output outputs/custom_cellulose.pdb
```

For one chain with any positive DP:

```bash
python builder.py --dp 37 --preset single -o outputs/chain_DP37.pdb
```

Question-and-answer mode:

```bash
python builder.py --interactive
```

Randomly oxidize 25% of the C6 sites on surface chains to carboxylate groups:

```bash
python builder.py --dp 20 --preset 18 --oxidation 0.25 \
  --oxidation-scope surface --oxidation-state deprotonated --seed 20260822 \
  -o outputs/cnf_18x20_ox25_surface_deprotonated.pdb
```

Oxidation is disabled unless `--oxidation` is supplied. Surface chains are the
chains in the top and bottom rows plus both boundary chains of every intermediate
row. The requested fraction is applied to all C6 sites eligible under the chosen
scope, rounded to the nearest whole site. Selection is random; use `--seed` to
make it reproducible. Deprotonated `COO-` is the default oxidation state, while
`--oxidation-state protonated` generates neutral `COOH`.

Available presets are `single`, `18` (`2,3,4,4,3,2`) and `36`
(`3,4,5,6,6,5,4,3`). Presets are conveniences only; `--layers` accepts any
positive row profile that fits within classic PDB limits.

## Command-line options

Run `python builder.py --help` at any time to display the built-in command
reference. All currently available options are listed below.

| Option | Meaning | Default |
| --- | --- | --- |
| `--dp N`, `--glucose N` | Degree of polymerization: the number of glucose units placed in **each** cellulose chain. `N` must be a positive integer. | `20` |
| `--layers N,N,...` | Custom transverse crystal geometry. Each integer is the number of parallel chains in one row, ordered from top to bottom. For example, `2,3,4,4,3,2` creates 6 rows and 18 chains. All values must be positive. | Not set |
| `--preset single` | Build one isolated cellulose chain. | Used when neither `--preset` nor `--layers` is supplied |
| `--preset 18` | Build the 18-chain geometry `2,3,4,4,3,2`. | Not set |
| `--preset 36` | Build the 36-chain geometry `3,4,5,6,6,5,4,3`. | Not set |
| `--output PATH`, `-o PATH` | Select the generated PDB filename and directory. Missing parent directories are created automatically. | `outputs/cellulose.pdb` |
| `--no-conect` | Omit all PDB `CONECT` records. This produces a smaller file but is not recommended for CHARMM-GUI or other carbohydrate readers that rely on explicit bonds. | Connectivity included |
| `--oxidation FRACTION` | Enable random C6 oxidation at a degree from `0.1` to `1.0`. | Disabled |
| `--oxidation-scope surface\|all` | Restrict eligible C6 sites to boundary chains or include every chain. | `surface` |
| `--oxidation-state deprotonated\|protonated` | Generate `COO-` or neutral `COOH`. | `deprotonated` |
| `--seed N` | Fix random selection for reproducible structures. | Random |
| `--interactive`, `-i` | Ask interactively for DP, layers, and output path instead of requiring those values on the command line. A blank layers answer selects one chain. | Disabled |
| `--help`, `-h` | Display the command reference and exit. | — |

`--layers` and `--preset` describe the same transverse arrangement and are
therefore mutually exclusive. The total number of residues is
`DP × total chains`; for `--dp 20 --layers 2,3,4,4,3,2`, this is
`20 × 18 = 360` glucose residues.

## CHARMM-GUI compatibility

PDB files generated by this version have been tested as direct input for
[CHARMM-GUI Glycan Reader & Modeler](https://www.charmm-gui.org/?doc=input/glycan).
Upload `outputs/cellulose.pdb` there as a PDB carbohydrate structure. The
reader can identify the glucose rings and β(1→4) connectivity, after which the
CHARMM-GUI workflow can generate the processed PDB and the CHARMM topology/PSF
files needed to prepare a molecular-dynamics system.

The Glycan Reader & Modeler interface can also edit the recognized cellulose
chains before building the system. This includes changing individual sugar
units, extending a chain, adding branches, selecting linkage positions and
anomeric configurations, and applying the chemical modifications made
available by CHARMM-GUI. These operations are performed in CHARMM-GUI after
upload; `dynamics_cellulose` currently supplies the experimentally based
starting crystal geometry.

Keep explicit connectivity enabled for this workflow (the default). The
CHARMM-GUI documentation explains that carbohydrate recognition uses the bond
graph supplied by PDB `CONECT` records. Compatibility described here applies
to the present cellulose Iβ output and the current Glycan Reader & Modeler
workflow; future project versions will expand and automate preparation for
additional simulation applications.

## Output guarantees

- experimental cellulose Iβ coordinates rather than accumulated geometric
  chain growth;
- fixed-column PDB records readable by VMD;
- unique chain IDs and segment names;
- explicit intramolecular and β-1,4 glycosidic `CONECT` records;
- chemically complete hydroxyls and terminal O4-H/O1-H groups for every finite chain;
- `TER` records separating independent cellulose chains;
- validation of residue count, C1–O4 distances and lattice spacing;
- provenance and DOI embedded as PDB `REMARK` records.

## Hydrogens and force fields

COD 4114994 contains carbon-bound hydrogens but does not locate every hydroxyl
hydrogen. The generator preserves all deposited heavy-atom and C-H coordinates
and models the missing O-H bonds at 0.96 A with tetrahedral H-O-C geometry. It
samples hydroxyl rotamers to avoid steric overlap. Internal residues receive
HO2, HO3, and HO6; each finite chain also receives HO4 at its non-reducing end
and O1-HO1 at its reducing end. Thus a DP *n* chain has the complete neutral
formula C6*n* H(10*n*+2) O(5*n*+1).

The output uses CHARMM carbohydrate atom names and the `BGLC` residue name,
plus explicit `CONECT` and chain `TER` records. It is an all-atom molecular
structure, but a PDB never contains force-field atom types, charges, angles, or
dihedrals; generating a PSF still requires a carbohydrate topology and the
β(1→4) linkage/terminal patches in the selected preparation program. Hydroxyl
orientations should be relaxed during the normal minimization/equilibration.

## C6 oxidation model

Oxidation converts the primary C6 alcohol into a planar carboxyl group. The
original `H61`, `H62`, `O6`, and `HO6` atoms are removed and replaced with
`O61` and `O62`, matching the atom naming observed in CHARMM-GUI carbohydrate
output. Protonated groups additionally contain `HO62`. Deprotonated sites each
have a nominal charge of `-1 e`; the PDB records composition and connectivity,
while partial atomic charges and force-field parameters remain the
responsibility of the subsequent topology/PSF preparation step.

For `surface` scope, eligible chains are those in the top and bottom rows plus
the two boundary chains in every intermediate row. For the 18-chain profile,
this selects 12 surface chains and leaves 6 interior chains unchanged. The
oxidation degree is calculated over all glucose residues in the eligible
chains. The nearest whole number of sites is sampled without replacement.

The generated coordinates are an initial chemical geometry based on the
CHARMM-GUI reference supplied during development: approximately 1.26 A C6-O
bonds and trigonal-planar carboxyl geometry. As with hydroxyl orientations,
oxidized structures should undergo normal force-field assignment,
minimization, and equilibration before production dynamics.

## Current scope and limits

Version 0.3 implements finite cellulose Iβ and optional C6 oxidation. Other allomorphs require their own
experimental datasets and symmetry implementations; they must not be
approximated by changing Iβ lattice dimensions.

Classic PDB supports at most 62 distinct one-character chain IDs, 99,999 atoms,
and coordinates below 10,000 Å in the fields used here. A future mmCIF writer
can remove those serialization limits.

## Tests

```bash
python -m unittest discover -s tests -v
```

The checks include odd-DP chains, the 18-chain profile, fixed-column PDB
formatting, explicit connectivity, and independent VMD reads.
