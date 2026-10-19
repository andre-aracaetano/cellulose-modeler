# Experimental C6 sulfate-half-ester model

## Chemical interpretation

Sulfuric-acid hydrolysis of cellulose can esterify surface-accessible hydroxyls
while removing amorphous material. The product is a sulfate half-ester, not a
direct carbon-sulfur substituent:

```text
cellulose-C6-OH + H2SO4 -> cellulose-C6-O-SO3H + H2O
cellulose-C6-O-SO3H -> cellulose-C6-O-SO3- + H+
```

The protonated model is neutral. The deprotonated model has nominal charge
`-1 e` per sulfate group. In aqueous CNC dispersions the anionic form is the
usual practical model above the strongly acidic regime.

## Evidence and scope

Sulfuric-acid-produced CNCs are experimentally reported to carry covalently
bound, negatively charged sulfate groups on their surfaces:

- Lin and Dufresne, *Nanoscale* 6 (2014), 5384,
  <https://doi.org/10.1039/C3NR06761K>.
- Gu et al., *Carbohydrate Polymers* 92 (2013), 1809-1816,
  <https://doi.org/10.1016/j.carbpol.2012.10.078>.

Liquid-state NMR has directly identified sulfate half-esters at C6 in CNC, but
the same study notes that low-abundance C2/C3 substitution cannot be excluded:

- King et al., *Biomacromolecules* 19 (2018), 2708-2720,
  <https://doi.org/10.1021/acs.biomac.8b00295>.

A review of CNC end-group and surface chemistry likewise states that a
definitive preferred localization of the sulfate half-esters has not been
established for all materials:

- Heise et al., *Angewandte Chemie International Edition* 60 (2021),
  66-87, <https://doi.org/10.1002/anie.202002433>.

Cellulose Modeler therefore implements a clearly labeled **C6 structural
hypothesis**, rather than claiming complete experimental regioselectivity.

## Atomistic transformation

For every selected residue, the model:

1. retains `C6`, `H61`, `H62`, and the experimental `O6` coordinate;
2. omits the alcohol hydrogen `HO6`;
3. places `S6` 1.60 Å from `O6` using the least-crowded sampled C-O-S rotamer;
4. places `O61`, `O62`, and `O63` tetrahedrally around sulfur;
5. uses 1.47 Å terminal S-O distances in the deprotonated model;
6. uses a 1.60 Å S-OH bond and a 0.96 Å O-H bond for `O63-HO63` in the
   protonated model.

These are initial covalent geometries. The crystallographic cellulose atoms
are not altered, and the new sulfate coordinates are not an optimized quantum
or force-field minimum.

## Site selection

`--sulfation-scope surface` selects top and bottom rows plus the two boundary
chains of each intermediate row. For the standard 18-chain profile, this is 12
surface chains. `--sulfation FRACTION` samples that eligible population without
replacement; `--seed` makes the realization reproducible.

`--sulfation-scope all` includes internal chains. It is provided for controlled
computational comparisons, not as the default representation of heterogeneous
sulfuric-acid hydrolysis.

`--sulfation-model paajanen` applies the existing TEMPO eligibility algorithm
without modification: each surface chain contributes its outward-facing
alternating C6 class, all candidates are pooled, and the requested number is
sampled globally. It does not divide the requested count equally between
chains. In an 18-chain, DP20 fibril there are 120 such candidates, so selecting
60 sulfates requires `--sulfation 0.50` and corresponds to 25% of the 240 C6
sites on surface chains.

## Example

```bash
python3 builder.py --allomorph ibeta --dp 20 --preset 18 \
  --sulfation 0.50 --sulfation-scope surface --sulfation-model paajanen \
  --sulfation-state deprotonated --seed 344 --charmm-gui \
  -o outputs/cnc_ibeta_18x20_sulfate25_seed344_charmm_gui.pdb
```

For this system, 240 C6 sites belong to surface chains and 60 are selected.
The nominal total charge contributed by sulfate is `-60 e` before counterions.

## CHARMM-GUI and force-field status

The output contains explicit `O6-S6`, `S6-O61`, `S6-O62`, and `S6-O63`
`CONECT` records. The `--charmm-gui` profile retains all sulfate heavy atoms
and retains `HO63` when the protonated state is requested.

This output is intentionally labeled experimental. Successful visual parsing
or residue recognition in CHARMM-GUI does not guarantee that a suitable
carbohydrate sulfate patch, partial charges, atom types, angles, dihedrals, and
nonbonded parameters were assigned. Those details and the generated PSF must be
inspected before minimization or molecular dynamics.
