# Crystallographic data provenance

`COD_4114994_cellulose_Ibeta.cif` was downloaded from the Crystallography Open
Database at <https://www.crystallography.net/cod/4114994.cif>.

The CIF itself states that all data on the COD site have been placed in the
public domain by the contributors. Its experimental source is:

Yoshiharu Nishiyama, Paul Langan, and Henri Chanzy. “Crystal Structure and
Hydrogen-Bonding System in Cellulose Iβ from Synchrotron X-ray and Neutron
Fiber Diffraction.” *Journal of the American Chemical Society* 124 (2002),
9074–9082. <https://doi.org/10.1021/ja0257319>.

COD identifier: 4114994.

`dynamics_cellulose` reads the fractional coordinates, applies the deposited
space-group operation, converts them through the experimental unit-cell matrix,
and constructs finite chains and crystallites. The distributed data directory
contains no coordinate output copied from Cellulose Builder or CHARMM-GUI.

The deposited model contains carbon-bound hydrogens but not all hydroxyl
hydrogens. The generator preserves the experimental heavy-atom and C-H
coordinates and clearly labels its added hydroxyl H coordinates as modeled.
They use a 0.96 A O-H bond, tetrahedral H-O-C geometry, and a steric rotamer
search. A force-field backend remains responsible for atom types, charges, and
bonded parameters before molecular dynamics.
