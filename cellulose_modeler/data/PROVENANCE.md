# Crystallographic data provenance

Cellulose Modeler contains independent crystallographic references for four
cellulose allomorphs. The generator never reads or imports coordinates produced
by Cellulose Builder or CHARMM-GUI.

## Cellulose Iα

`COD_4114383_cellulose_Ialpha.cif` was downloaded directly from the
Crystallography Open Database (COD), where it is distributed under CC0:

- COD: <https://www.crystallography.net/cod/4114383.cif>
- Yoshiharu Nishiyama, Junji Sugiyama, Henri Chanzy, and Paul Langan.
  “Crystal Structure and Hydrogen Bonding System in Cellulose Iα from
  Synchrotron X-ray and Neutron Fiber Diffraction.” *Journal of the American
  Chemical Society* 125 (2003), 14300–14306.
  <https://doi.org/10.1021/ja037055w>

The deposited model is a one-chain triclinic `P1` unit cell. Its crystallographic
`a` vector is the chain axis. Cellulose Modeler preserves the two experimental glucosyl
residues forming the cellobiose repeat and rotates the final Cartesian frame so
the chain axis is written along PDB `z`.

## Cellulose Iβ

`COD_4114994_cellulose_Ibeta.cif` was downloaded directly from COD under CC0:

- COD: <https://www.crystallography.net/cod/4114994.cif>
- Yoshiharu Nishiyama, Paul Langan, and Henri Chanzy. “Crystal Structure and
  Hydrogen-Bonding System in Cellulose Iβ from Synchrotron X-ray and Neutron
  Fiber Diffraction.” *Journal of the American Chemical Society* 124 (2002),
  9074–9082. <https://doi.org/10.1021/ja0257319>

The model contains two crystallographically independent parallel chains in a
monoclinic `P2_1` unit cell. The deposited symmetry operation generates the
second glucosyl residue of each cellobiose repeat.

## Cellulose II

`cellulose_II_Langan2001.cif` records the fractional coordinates and unit-cell
parameters reported in the primary high-resolution synchrotron study:

- Paul Langan, Yoshiharu Nishiyama, and Henri Chanzy. “X-ray Structure of
  Mercerized Cellulose II at 1 Å Resolution.” *Biomacromolecules* 2 (2001),
  410–416. <https://doi.org/10.1021/bm005612q>

The coordinates were transcribed into a minimal, human-readable CIF maintained
by this project. The original model has two independent antiparallel chains in
`P2_1`. Cellulose Modeler explicitly expands one chain in `+c` and the other in `-c`,
including the crystallographic lateral translation required to preserve the
C1–O4 linkage. The X-ray model does not experimentally locate hydroxyl
hydrogens; Cellulose Modeler models those positions and declares them as modeled.

The earlier neutron refinement of cellulose II is also scientifically relevant
to its hydrogen-bonding alternatives, but it is not the coordinate source used
by the present implementation:

- Paul Langan, Yoshiharu Nishiyama, and Henri Chanzy. “A Revised Structure and
  Hydrogen-Bonding System in Cellulose II from a Neutron Fiber Diffraction
  Analysis.” *Journal of the American Chemical Society* 121 (1999), 9940–9946.
  <https://doi.org/10.1021/ja9916254>

## Cellulose III_I

`cellulose_IIII_Wada2004.cif` records the fractional coordinates and unit-cell
parameters from the primary synchrotron X-ray and neutron study:

- Masahisa Wada, Henri Chanzy, Yoshiharu Nishiyama, and Paul Langan.
  “Cellulose III_I Crystal Structure and Hydrogen Bonding by Synchrotron X-ray
  and Neutron Fiber Diffraction.” *Macromolecules* 37 (2004), 8548–8555.
  <https://doi.org/10.1021/ma0485585>

The coordinates were transcribed into a minimal project-owned CIF. The model is
a one-chain monoclinic `P2_1` structure with parallel chains. Cellulose Modeler applies
the experimental screw symmetry to generate the cellobiose repeat.

## Coordinate treatment common to all references

The implementation:

1. reads experimental fractional heavy-atom and carbon-bound-hydrogen positions;
2. applies the appropriate `P1` or `P2_1` expansion;
3. uses a general triclinic unit-cell matrix;
4. orders residues by actual covalent C1–O4 connectivity;
5. preserves parallel or antiparallel chain direction;
6. rotates only the global Cartesian frame to place the chain axis along `z`;
7. creates finite reducing and non-reducing ends; and
8. models missing hydroxyl hydrogens at 0.96 Å with a steric rotamer search.

Experimental heavy atoms are not minimized or altered by the generator.
Force-field software remains responsible for atom types, partial charges,
bonded parameters, minimization, and equilibration.

## Scientific limitation of C6 oxidation

The current oxidation implementations, particularly the Paajanen alternating
surface-site model, are validated only for cellulose Iβ. Cellulose Modeler rejects
oxidation requests for Iα, II, and III_I rather than silently transferring an
Iβ-specific surface hypothesis to different crystal packings.
