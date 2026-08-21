# dynamics_cellulose

`dynamics_cellulose` is a dependency-free Python generator for finite cellulose
Iβ chains and crystallites with arbitrary degrees of polymerization and
transverse row profiles.

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

See [the detailed data provenance](dynamics_cellulose/data/PROVENANCE.md) and
[`CITATION.cff`](CITATION.cff).

## Quick start

Run from this directory without installing anything:

```bash
PYTHONPATH=. python3 -m dynamics_cellulose --dp 20 --preset 18 -o cellulose_18.pdb
```

An arbitrary transverse arrangement is a comma-separated list of chain counts
from top to bottom:

```bash
PYTHONPATH=. python3 -m dynamics_cellulose \
  --dp 16 \
  --layers 3,5,7,7,5,3 \
  --output custom_cellulose.pdb
```

For one chain with any positive DP:

```bash
PYTHONPATH=. python3 -m dynamics_cellulose --dp 37 --preset single -o chain_DP37.pdb
```

Question-and-answer mode:

```bash
PYTHONPATH=. python3 -m dynamics_cellulose --interactive
```

Available presets are `single`, `18` (`2,3,4,4,3,2`) and `36`
(`3,4,5,6,6,5,4,3`). Presets are conveniences only; `--layers` accepts any
positive row profile that fits within classic PDB limits.

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

## Current scope and limits

Version 0.2 implements finite cellulose Iβ. Other allomorphs require their own
experimental datasets and symmetry implementations; they must not be
approximated by changing Iβ lattice dimensions.

Classic PDB supports at most 62 distinct one-character chain IDs, 99,999 atoms,
and coordinates below 10,000 Å in the fields used here. A future mmCIF writer
can remove those serialization limits.

## Tests

```bash
PYTHONPATH=. python3 -m unittest discover -s tests -v
```

The checks include odd-DP chains, the 18-chain profile, fixed-column PDB
formatting, explicit connectivity, and independent VMD reads.
