# Experimental cellulose allomorphs and CHARMM-GUI validation

## Implemented references

| CLI key | Crystal system | Space group | Independent chains | Chain packing | Primary source |
|---|---|---|---:|---|---|
| `ialpha` | triclinic | P1 | 1 | parallel | Nishiyama et al. 2003, <https://doi.org/10.1021/ja037055w> |
| `ibeta` | monoclinic | P2_1 | 2 | parallel | Nishiyama et al. 2002, <https://doi.org/10.1021/ja0257319> |
| `ii` | monoclinic | P2_1 | 2 | antiparallel | Langan et al. 2001, <https://doi.org/10.1021/bm005612q> |
| `iii-i` | monoclinic | P2_1 | 1 | parallel | Wada et al. 2004, <https://doi.org/10.1021/ma0485585> |

## Experimental unit cells

| Allomorph | a (Å) | b (Å) | c (Å) | α (°) | β (°) | γ (°) | Crystallographic chain axis |
|---|---:|---:|---:|---:|---:|---:|---|
| Iα | 10.400 | 6.717 | 5.962 | 80.37 | 118.08 | 114.80 | a |
| Iβ | 7.784 | 8.201 | 10.380 | 90.00 | 90.00 | 96.55 | c |
| II | 8.100 | 9.030 | 10.310 | 90.00 | 90.00 | 117.10 | c |
| III_I | 4.450 | 7.850 | 10.310 | 90.00 | 90.00 | 105.10 | c |

Dynamics rotates the global Cartesian frame so the chain axis is written along
PDB `z`; it does not deform or minimize the experimental heavy-atom geometry.

## Local 18-chain, DP20 validation

All three new allomorphs generate 18 chains, 360 residues, 7,614 atoms in the
complete standalone structure. Each contains exactly 342 C1–O4 glycosidic
links. The hydroxyl-H-free profile contains 6,498 atoms, but it is not offered
for cellulose II or III_I because no validated Glycan Reader workflow is
available for those allomorphs.

| Allomorph | C1–O4 range (Å) | Minimum transverse chain-axis separation (Å) | Minimum interchain heavy-atom contact (Å) |
|---|---:|---:|---:|
| Iα | 1.402–1.406 | 5.260 | 2.770 |
| II | 1.393–1.395 | 4.366 | 2.526 |
| III_I | 1.378 | 4.450 | 2.617 |

No generated PDB contains an interchain covalent bond. Cellulose II alternates
chains whose covalent residue order advances along `+z` and `-z`; Iα and III_I
chains are parallel.

## Files for external CHARMM-GUI validation

```text
outputs/cellulose_ialpha_18x20_charmm_gui.pdb
```

Upload each supported candidate through Glycan Reader using **CHARMM** input interpretation. The
expected initial recognition result is:

- 18 glycan segments;
- 20 β-D-glucose residues per segment;
- 360 total glucose residues;
- no glucuronate, quinovose, or 6-deoxy residue;
- β(1→4) connectivity throughout every chain.

After CHARMM-GUI generates PDB/PSF files, validation must additionally confirm:

1. 360 `BGLC` residues and zero modified residues;
2. total solute charge zero before ions;
3. 342 glycosidic O4–C1 bonds and no intersegment covalent bond;
4. intact C6–H61/H62–O6 chemistry in every residue;
5. heavy-atom RMSD of zero at PDB precision relative to the input;
6. preservation of parallel directions for Iα and III_I; and
7. preservation of the crystallographic chain directions.

Recognition as the same β(1→4)-glucose sequence is necessary but not sufficient:
the allomorph is encoded by three-dimensional packing, not by a distinct CHARMM
residue name. Therefore, coordinate preservation and chain direction must be
checked after topology generation.

## Cellulose-II interoperability result

The standalone cellulose-II model passed the local crystallographic and bond
graph checks: 18 antiparallel chains, DP20, 360 glucose residues, 342
intrachain O4–C1 links, no interchain covalent bonds, and a 1.393–1.395 Å
glycosidic-distance range. Nevertheless, CHARMM-GUI Glycan Reader 3.7
reproducibly detached the terminal glucose of the same physical chain as
`HETA/BGL`, leaving 359 `BGLC` residues in the generated PDB/PSF when the
heterogen was not selected.

Tests changing the PDB chain identifier, SEGID, record order, coordinate
rounding below 0.001 Å, and even isolating the affected chain did not remove
the split. Treating the detached glucose with a ligand force field or merely
renaming it to `BGLC` would not recreate the missing β(1→4) patch. For that
reason cellulose II is supported only as a standalone PDB in this release;
`--allomorph ii --charmm-gui` fails explicitly rather than emitting an input
that has not passed topology validation.

## Oxidation restriction

C6 oxidation is deliberately rejected for `ialpha`, `ii`, and `iii-i` in this
release. The Paajanen eligibility model remains specific to the validated Iβ
surface hypothesis. Support for oxidation of another allomorph requires a
separate surface and C6-orientation analysis.
