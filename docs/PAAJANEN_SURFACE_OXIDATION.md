# Paajanen-inspired surface oxidation

## Purpose and scientific scope

The `paajanen` oxidation model generates an atomistic hypothesis for
TEMPO-oxidized cellulose Iβ surfaces. It is intended for studies in which C6
carboxylates should be restricted to solvent-facing, alternating sites on the
crystallite boundary. It does not simulate the TEMPO reaction, predict reaction
kinetics, or reconstruct a unique experimental charge distribution.

This distinction matters: the implementation is a final-structure builder for
classical simulation, not a reactive model.

## Why cellulose Iβ and C6?

The crystal coordinates originate from the experimental cellulose Iβ structure
of Nishiyama, Langan, and Chanzy. Iβ is the dominant native cellulose allomorph
in higher plants and is therefore a defensible structural starting point for
plant-derived cellulose nanofibrils. Real samples can additionally contain
defects, disordered regions, twisting, variable widths, and chain ends that the
ideal crystal model does not reproduce.

TEMPO-mediated oxidation is regioselective toward the primary C6 alcohol. The
idealized chemical change represented by the builder is:

```text
C6-CH2OH  ->  C6-CHO  ->  C6-COOH/COO-
```

The carbon C6 is retained. In the default deprotonated state, each oxidized site
has nominal formal charge -1. Oxidation changes hydrogen bonding, hydration,
counterion association, and electrostatic interactions, but it does not itself
break a β(1→4) glycosidic bond in this non-reactive model.

Experimental TEMPO-oxidation studies support preferential modification of
accessible crystalline surfaces while much of the cellulose-I core remains
ordered. Dense surface oxidation has also yielded evidence for alternating
glucose/glucuronic-acid sequences, consistent with the two-glucosyl cellobiose
repeat.

## Why use the Paajanen eligibility model?

Paajanen *et al.* modeled TEMPO-oxidized cellulose nanofibrils by allowing
carboxylate substitution at every second hydroxymethyl C6 group of each surface
chain. They studied several functionalization levels and multiple spatial charge
distributions at intermediate levels. Their results showed that equal numbers
of carboxylates do not necessarily produce equal interactions: local charge
density, counterion placement, hydration, and the relative alignment of charged
regions all matter.

This makes the model useful for Cellulose Modeler because it provides:

1. an explicit definition of the eligible-site denominator;
2. a structural connection to the alternating cellobiose repeat;
3. exclusion of internal chains from a surface-oxidation hypothesis;
4. reproducible random sampling through `--seed`; and
5. a natural way to generate independent charge-pattern replicas.

It should not be described as universally correct for every TEMPO-oxidized
sample. It is the more physically defensible option implemented here when the
stated hypothesis is **alternating, surface-limited oxidation of an ideal Iβ
crystallite**. Defects, amorphous regions, reaction conditions, fibril morphology,
and experimental history can lead to distributions not captured by this model.

## Why a global odd/even rule is incorrect

Successive glucosyl residues are related by the cellobiose repeat, so their C6
groups form two alternating orientational classes. However, residue-number
parity is only a PDB convention. A class that faces outward on one fibril face
can face inward on the opposite face.

Applying the same odd or even class to every chain can therefore bury some
carboxylates in the crystallite. An internal formal charge would artificially
favor water and counterion penetration and could distort crystal cohesion and
adsorbate accessibility.

Cellulose Modeler consequently chooses the class separately for every surface chain.
For chain `i`, it defines an outward transverse direction from the fibril center
to the chain and evaluates each residue's normalized C5→C6 vector:

```text
p_j = unit(C5->C6)_j dot unit(outward)_i
```

The mean projection is calculated for the two alternating classes. The class
with the larger mean outward projection becomes eligible. This projection rule
is a Cellulose Modeler adaptation motivated by the surface chemistry; it is **not** an
atom-selection equation reported in the Paajanen paper.

For the regular 18-chain `2,3,4,4,3,2` profile at DP20, the eligible classes are:

| Surface chain | Eligible residue class |
|---:|:---|
| 1 | even |
| 2 | odd |
| 3 | even |
| 5 | odd |
| 6 | even |
| 9 | odd |
| 10 | even |
| 13 | odd |
| 14 | even |
| 16 | odd |
| 17 | even |
| 18 | odd |

The center-to-chain approximation is suitable for this convex regular profile.
Strongly twisted, concave, or irregular morphologies may require local surface
normals or a more detailed neighbor-based surface definition.

## Why SASA is not the reaction criterion

Solvent-accessible surface area (SASA) measures geometrical accessibility to a
probe. It does not directly measure TEMPO/oxoammonium reaction free energy or
the orientation required for chemical conversion. SASA also depends on probe
radius, atomic radii, the chosen static conformation, thresholds, and chain-end
exposure. Oxidation and hydration can themselves change C6 conformations and
therefore alter accessibility.

SASA can be valuable as a diagnostic to compare faces or detect buried choices,
but the current implementation does not claim that a water-probe SASA is a
chemical reactivity model.

## Degree definitions for the 18-chain DP20 fibril

The `2,3,4,4,3,2` fibril contains:

- 18 chains × 20 residues = 360 total anhydroglucose units (AGUs);
- 12 surface chains × 20 residues = 240 AGUs on surface-classified chains;
- one alternating class per surface chain = 12 × 10 = 120 eligible C6 sites.

In `paajanen` mode, `--oxidation` is the fraction of the **120 eligible sites**,
not the fraction of all 360 AGUs.

Thus, keeping the three denominators explicit:

| Request | Carboxylates | Eligible fraction | Surface-chain C6 fraction | All-AGU fraction |
|---|---:|---:|---:|---:|
| `--oxidation 0.25` | 30 | 25% | 12.5% | 8.33% |
| `--oxidation 0.50` | 60 | 50% | 25% | 16.67% |
| `--oxidation 0.75` | 90 | 75% | 37.5% | 25% |
| `--oxidation 1.00` | 120 | 100% | 50% | 33.33% |

For the selected surface-25% system:

```text
60 / 240 = 25% of C6 sites belonging to surface chains
60 / 120 = 50% of Paajanen-eligible C6 sites
60 / 360 = 16.67% of all AGUs
```

For a distinct global-25% system:

```text
90 / 360 = 25% of all AGUs
90 / 120 = 75% of Paajanen-eligible C6 sites
90 / 240 = 37.5% of C6 sites belonging to surface chains
```

Using 162.14 g/mol as the approximate unmodified AGU molar mass gives:

```text
1000 / 162.14 = 6.17 mmol AGU/g
0.25 * 6.17 = 1.54 mmol COO-/g
```

This conversion is approximate because oxidation and counterion form change the
sample mass convention.

## Reproducible commands

Twenty-five percent of eligible sites, equivalent to 8.33% global oxidation:

```bash
python builder.py --dp 20 --preset 18 \
  --oxidation-model paajanen --oxidation-scope surface \
  --oxidation 0.25 --oxidation-state deprotonated \
  --seed 344 --charmm-gui \
  -o outputs/tocnf_18x20_paajanen_eligible25_seed344.pdb
```

Twenty-five percent of surface-chain C6 sites, equivalent to 50% of eligible
sites and used for the validated surface-25% structure:

```bash
python builder.py --dp 20 --preset 18 \
  --oxidation-model paajanen --oxidation-scope surface \
  --oxidation 0.50 --oxidation-state deprotonated \
  --seed 344 --charmm-gui \
  -o outputs/tocnf_18x20_paajanen_surface25_seed344_charmm_gui.pdb
```

Twenty-five percent global oxidation, equivalent to 75% of eligible sites:

```bash
python builder.py --dp 20 --preset 18 \
  --oxidation-model paajanen --oxidation-scope surface \
  --oxidation 0.75 --oxidation-state deprotonated \
  --seed 344 --charmm-gui \
  -o outputs/tocnf_18x20_paajanen_global25_seed344.pdb
```

The PDB records the oxidation model, requested fraction, eligible-site count,
oxidized-site count, protonation state, and random seed in `REMARK` records.

## Statistical recommendation

Partially oxidized surfaces should be studied with multiple seeds. Keeping 90
carboxylates fixed while changing the seed isolates sensitivity to the spatial
charge distribution. A single seed is one microscopic realization, not a
statistical representation of the experimental material.

## Main references

- Nishiyama Y, Langan P, Chanzy H. *Crystal Structure and Hydrogen-Bonding
  System in Cellulose Iβ from Synchrotron X-ray and Neutron Fiber Diffraction*.
  J Am Chem Soc. 2002;124:9074–9082.
  <https://doi.org/10.1021/ja0257319>
- Saito T, Isogai A. *TEMPO-Mediated Oxidation of Native Cellulose*.
  Biomacromolecules. 2004;5:1983–1989.
  <https://doi.org/10.1021/bm0497769>
- Hirota M, Furihata K, Saito T, Kawada T, Isogai A.
  *Glucose/Glucuronic Acid Alternating Copolysaccharides Prepared from
  TEMPO-Oxidized Native Celluloses by Surface Peeling*.
  Angew Chem Int Ed. 2010;49:7670–7672.
  <https://doi.org/10.1002/anie.201003848>
- Isogai A, Saito T, Fukuzumi H. *TEMPO-Oxidized Cellulose Nanofibers*.
  Nanoscale. 2011;3:71–85. <https://doi.org/10.1039/C0NR00583E>
- Paajanen AT, Sonavane Y, Ignasiak D, Ketoja JA, Maloney T, Paavilainen S.
  *Atomistic Molecular Dynamics Simulations on the Interaction of
  TEMPO-Oxidized Cellulose Nanofibrils in Water*.
  Cellulose. 2016;23:3449–3462.
  <https://doi.org/10.1007/s10570-016-1076-x>
