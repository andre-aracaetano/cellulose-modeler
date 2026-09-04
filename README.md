# Cellulose Modeler

`cellulose_modeler` is a dependency-free Python generator for finite cellulose
Iα, Iβ, II, and III_I chains and crystallites with arbitrary degrees of
polymerization and transverse row profiles. The user-facing launcher is `builder.py`; reusable
implementation modules live in `cellulose_modeler/`, and generated structures go
to `outputs/` by default.

Version 0.7 provides the renamed Cellulose Modeler package and command-line
interface. Version 0.6 added experimental Iα, II, and III_I references and a general
allomorph-aware crystal engine. Version 0.4 added the validated `--charmm-gui`
output profile. It prevents an
observed carbohydrate-recognition ambiguity without changing the cellulose
heavy-atom structure. Version 0.3 introduced optional, reproducible C6
oxidation to protonated carboxylic acid (`COOH`) or deprotonated carboxylate
(`COO-`) groups.

## Scientific basis and provenance

The program constructs every allomorph from an independently bundled
experimental crystallographic reference. Iα and Iβ use the public-domain COD
deposits 4114383 and 4114994. Cellulose II and III_I use minimal CIFs transcribed
from their primary high-resolution crystallographic publications. Generation
never reads coordinates produced by Cellulose Builder or CHARMM-GUI.

Experimental sources:

- **Iα:** Nishiyama, Sugiyama, Chanzy, and Langan, *JACS* 125 (2003),
  14300–14306. <https://doi.org/10.1021/ja037055w>; COD 4114383.
- **Iβ:** Nishiyama, Langan, and Chanzy, *JACS* 124 (2002), 9074–9082.
  <https://doi.org/10.1021/ja0257319>; COD 4114994.
- **II:** Langan, Nishiyama, and Chanzy, *Biomacromolecules* 2 (2001),
  410–416. <https://doi.org/10.1021/bm005612q>.
- **III_I:** Wada, Chanzy, Nishiyama, and Langan, *Macromolecules* 37 (2004),
  8548–8555. <https://doi.org/10.1021/ma0485585>.

The implementation parses the selected CIF, applies `P1` or `P2_1` symmetry,
uses a general triclinic cell matrix, repeats the crystallographic motif to the
requested DP, and writes the selected row profile as PDB. It preserves the
parallel packing of Iα, Iβ, and III_I and the antiparallel origin/center chains
of cellulose II.

See [the detailed data provenance](cellulose_modeler/data/PROVENANCE.md) and
[`CITATION.cff`](CITATION.cff). The geometry and external validation protocol
are recorded in [`docs/ALLOMORPH_VALIDATION.md`](docs/ALLOMORPH_VALIDATION.md).

## Quick start

Download the project, open a terminal in its root directory, and run:

```bash
python builder.py --dp 20 --preset 18
```

Select another experimental allomorph with `--allomorph`:

```bash
python builder.py --allomorph ialpha --dp 20 --preset 18 --charmm-gui \
  -o outputs/cellulose_ialpha_18x20_charmm_gui.pdb
python builder.py --allomorph ii --dp 20 --preset 18 \
  -o outputs/cellulose_ii_18x20.pdb
python builder.py --allomorph iii-i --dp 20 --preset 18 \
  -o outputs/cellulose_iii_i_18x20.pdb
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

Use the alternating surface-site model of Paajanen *et al.* and oxidize 25%
of all C6 sites on the 12 surface chains. Because only one alternating class
is eligible, this requires selecting 50% of the Paajanen-eligible sites:

```bash
python builder.py --dp 20 --preset 18 --oxidation 0.50 \
  --oxidation-scope surface --oxidation-model paajanen \
  --oxidation-state deprotonated --seed 344 --charmm-gui \
  -o outputs/tocnf_18x20_paajanen_surface25_seed344_charmm_gui.pdb
```

For a PDB that will be uploaded to CHARMM-GUI, add `--charmm-gui`:

```bash
python builder.py --dp 20 --preset 18 --charmm-gui \
  -o outputs/cnf_18x20_charmm_gui.pdb

python builder.py --dp 20 --preset 18 --oxidation 0.25 \
  --oxidation-scope surface --oxidation-state deprotonated --seed 344 \
  --charmm-gui -o outputs/tocnf25_18x20_charmm_gui.pdb
```

Oxidation is disabled unless `--oxidation` is supplied. Surface chains are the
chains in the top and bottom rows plus both boundary chains of every intermediate
row. With the default `random-surface` model, the requested fraction is applied
to every C6 in the chosen scope. With `paajanen`, it is applied only to the
alternating eligible sites described below. Counts are rounded to the nearest
whole site. Selection is random; use `--seed` to make it reproducible.
Deprotonated `COO-` is the default oxidation state, while `--oxidation-state
protonated` generates neutral `COOH`.

Available presets are `single`, `18` (`2,3,4,4,3,2`) and `36`
(`3,4,5,6,6,5,4,3`). Presets are conveniences only; `--layers` accepts any
positive row profile that fits within classic PDB limits.

## Command-line options

Run `python builder.py --help` at any time to display the built-in command
reference. All currently available options are listed below.

| Option | Meaning | Default |
| --- | --- | --- |
| `--allomorph ialpha\|ibeta\|ii\|iii-i` | Select the experimental cellulose crystal form. | `ibeta` |
| `--dp N`, `--glucose N` | Degree of polymerization: the number of glucose units placed in **each** cellulose chain. `N` must be a positive integer. | `20` |
| `--layers N,N,...` | Custom transverse crystal geometry. Each integer is the number of parallel chains in one row, ordered from top to bottom. For example, `2,3,4,4,3,2` creates 6 rows and 18 chains. All values must be positive. | Not set |
| `--preset single` | Build one isolated cellulose chain. | Used when neither `--preset` nor `--layers` is supplied |
| `--preset 18` | Build the 18-chain geometry `2,3,4,4,3,2`. | Not set |
| `--preset 36` | Build the 36-chain geometry `3,4,5,6,6,5,4,3`. | Not set |
| `--output PATH`, `-o PATH` | Select the generated PDB filename and directory. Missing parent directories are created automatically. | `outputs/cellulose.pdb` |
| `--no-conect` | Omit all PDB `CONECT` records. This produces a smaller file but is not recommended for CHARMM-GUI or other carbohydrate readers that rely on explicit bonds. | Connectivity included |
| `--charmm-gui` | Preserve heavy atoms and C-H atoms, omit hydroxyl-H coordinates for topology-based reconstruction, and order residues for `O4(i)-C1(i+1)` recognition. | Disabled |
| `--oxidation FRACTION` | Enable random C6 oxidation at a degree from `0.1` to `1.0`. | Disabled |
| `--oxidation-scope surface\|all` | Restrict eligible C6 sites to boundary chains or include every chain. | `surface` |
| `--oxidation-model random-surface\|paajanen` | Use every C6 in the selected scope, or alternating C6 sites on surface chains. | `random-surface` |
| `--oxidation-state deprotonated\|protonated` | Generate `COO-` or neutral `COOH`. | `deprotonated` |
| `--seed N` | Fix random selection for reproducible structures. | Random |
| `--interactive`, `-i` | Ask interactively for DP, layers, and output path instead of requiring those values on the command line. A blank layers answer selects one chain. | Disabled |
| `--help`, `-h` | Display the command reference and exit. | — |

`--layers` and `--preset` describe the same transverse arrangement and are
therefore mutually exclusive. The total number of residues is
`DP × total chains`; for `--dp 20 --layers 2,3,4,4,3,2`, this is
`20 × 18 = 360` glucose residues.

## CHARMM-GUI compatibility

The PDB generator is fully functional as a standalone structure builder. Its
outputs, including non-oxidized cellulose and C6-oxidized nanofibrils, have
been successfully recognized by both
[CHARMM-GUI PDB Reader & Manipulator](https://www.charmm-gui.org/?doc=input/pdbreader)
and [Glycan Reader & Modeler](https://www.charmm-gui.org/?doc=input/glycan).
The readers identify the glucose rings, β(1→4) connectivity, chain boundaries,
and supported C6 modifications. Structures processed this way have produced
CHARMM-formatted PDB and PSF files and have been used successfully in short
molecular-dynamics tests.

This validated statement currently applies to Iβ (including the supported
Paajanen oxidation workflow) and to non-oxidized Iα. Cellulose II and III_I
remain available as standalone crystallographic PDBs, but `--charmm-gui` is
rejected for both because no validated Glycan Reader workflow is available.
For cellulose II specifically, Glycan Reader 3.7 repeatedly detached the
terminal glucose of one chain in the tested 18-chain, DP20 antiparallel fibril
as `HETA/BGL`. Selecting that entry
would produce an isolated ligand rather than the required β(1→4)-linked chain.
This is an interoperability limitation, not evidence that the cellulose-II or
III_I coordinates or bond graphs are chemically invalid.

The validated workflow for molecular dynamics is:

1. Generate the desired cellulose structure with `--charmm-gui`.
2. Keep explicit `CONECT` records enabled, as they describe the carbohydrate
   bond graph used during recognition.
3. Upload the generated PDB and select **CHARMM** as its input format in
   CHARMM-GUI (rather than generic PDB interpretation).
4. Inspect the recognized chains, linkages, terminal groups, and oxidation
   states.
5. Ask CHARMM-GUI to generate the processed coordinates, PSF, and the topology,
   parameter, and input files appropriate for the user's selected simulation
   engine.
6. Perform the normal minimization and equilibration before production
   dynamics.

> **Recognition note:** During validation, the generic reader occasionally
> interpreted intact C6 alcohols as `BQUI` (beta-D-quinovose) even though the
> source PDB contained `O6`, `HO6`, `H61`, and `H62` and no `H63`. This is an
> ambiguity at the interface between coordinate input and automatic glycan
> recognition, not evidence that CHARMM-GUI changes cellulose at random. The
> validated `--charmm-gui` profile avoids it by leaving hydroxyl-H placement to
> the force-field topology and using the residue direction expected by the
> reader. Always inspect the proposed sequence before building a PSF. Manual
> `BQUI` to glucose replacement is discouraged because it can silently change
> a neighboring linkage from `1→4` to `1→1`.

The Glycan Reader & Modeler interface can also edit the recognized cellulose
chains before building the system. This includes changing individual sugar
units, extending a chain, adding branches, selecting linkage positions and
anomeric configurations, and applying the chemical modifications made
available by CHARMM-GUI. These operations are performed in CHARMM-GUI after
upload; Cellulose Modeler currently supplies the experimentally based
starting crystal geometry.

Keep explicit connectivity enabled for this workflow (the default). The
CHARMM-GUI documentation explains that carbohydrate recognition uses the bond
graph supplied by PDB `CONECT` records. A standalone PDB does not contain the
force-field atom types, partial charges, angles, dihedrals, or parameters held
by a PSF and its associated force-field files. This separation does not make
the current generator incomplete: it deliberately produces the molecular
structure, while CHARMM-GUI performs simulation-specific topology preparation.

## Roadmap: automated CHARMM-GUI preparation

The next planned integration is an optional link to the official CHARMM-GUI
API. The intended workflow will submit the PDB generated locally, request PDB
Reader or Glycan Reader processing, monitor the remote job, and download the
resulting CHARMM-formatted coordinates and topology package. This would remove
the current manual upload step without changing the standalone or
dependency-free PDB generation workflow.

Authentication, job-status, and result-download endpoints are publicly
documented by CHARMM-GUI. However, the public documentation does not currently
specify the submission endpoint and request fields for PDB Reader or Glycan
Reader. The project has requested this information from the CHARMM-GUI team.
API integration will remain a roadmap item until an official, supportable
submission workflow is available; no private endpoint or fragile browser
automation will be treated as a production interface.

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

The standard output uses CHARMM carbohydrate atom names and the `BGLC` residue name,
plus explicit `CONECT` and chain `TER` records. It is an all-atom molecular
structure, but a PDB never contains force-field atom types, charges, angles, or
dihedrals; generating a PSF still requires a carbohydrate topology and the
β(1→4) linkage/terminal patches in the selected preparation program. Hydroxyl
orientations should be relaxed during the normal minimization/equilibration.
With `--charmm-gui`, only hydroxyl hydrogens are omitted from the upload PDB;
CHARMM-GUI adds them from the selected topology.

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
this selects 12 surface chains and leaves 6 interior chains unchanged. With
the default `random-surface` model, oxidation degree is calculated over
all glucose residues in the eligible chains. The nearest whole number of sites
is sampled without replacement. The `paajanen` denominator is described below.

The generated coordinates are an initial chemical geometry based on the
CHARMM-GUI reference supplied during development: approximately 1.26 A C6-O
bonds and trigonal-planar carboxyl geometry. As with hydroxyl orientations,
oxidized structures should undergo normal force-field assignment,
minimization, and equilibration before production dynamics.

### Paajanen alternating surface-site model

The optional `--oxidation-model paajanen` mode follows the eligibility model
reported by A. Paajanen, Y. Sonavane, D. Ignasiak, J. A. Ketoja, T. Maloney,
and S. Paavilainen, “Atomistic molecular dynamics simulations on the
interaction of TEMPO-oxidized cellulose nanofibrils in water,” *Cellulose* 23
(2016), 3449–3462, <https://doi.org/10.1007/s10570-016-1076-x>.

That study randomly substituted carboxylates among every second C6
hydroxymethyl group of each surface chain. This is consistent with the
alternating glucosyl environments of the cellulose Iβ 2_1 helical repeat and
with experimental reports of alternating glucose/glucuronate structures on
oxidized cellulose-I microfibril surfaces. Because the solvent-facing class is
not the same residue-number parity on opposite sides of a fibril,
Cellulose Modeler chooses it separately for each surface chain. It compares
the mean projection of the two alternating classes' C5→C6 vectors onto that
chain's outward transverse direction and retains the class pointing farther
out of the crystallite. No SASA or reaction-energy calculation is implied.
This outward-facing rule is a documented geometric adaptation made by
Cellulose Modeler; it is not presented as an atom-selection algorithm
published by Paajanen *et al.*

In this mode, `--oxidation FRACTION` uses the alternating sites—not every C6
on a surface chain—as its denominator. For the 18-chain, DP20 profile, 12
surface chains each contribute 10 alternating sites, giving 120 eligible C6
sites. Thus `--oxidation 0.25` creates 30 carboxylates (25% of 120, or 8.33%
of all 360 anhydroglucose units). Different `--seed` values create independent
spatial distributions at the same functionalization level, as required for
studying the configuration sensitivity emphasized by Paajanen *et al.* The
mode is defined only for `--oxidation-scope surface`.

For the 18-chain, DP20 profile, 25% oxidation of the 240 C6 sites belonging to
surface chains means 60 carboxylates. This is 50% of the 120 alternating
Paajanen-eligible sites, so the correct request is `--oxidation 0.50`.

By contrast, 25% oxidation relative to all 360 anhydroglucose units
means 90 carboxylates. Because 90 is 75% of the 120 Paajanen-eligible sites,
generate that model with `--oxidation 0.75`, not `--oxidation 0.25`.

See [`docs/PAAJANEN_SURFACE_OXIDATION.md`](docs/PAAJANEN_SURFACE_OXIDATION.md)
for the scientific rationale, degree definitions, validation and limitations.

## Current scope and limits

Version 0.6 implements experimental Iα, Iβ, II, and III_I structures with
locally tested crystallographic geometry and PDB connectivity. The CHARMM-GUI
profile is supported only for Iα and Iβ; II and III_I remain standalone PDB
outputs and are not advertised as Glycan Reader compatible.

C6 oxidation is currently restricted to Iβ. In particular, the Paajanen model
is not transferred automatically to other allomorphs because their C6
orientations, exposed faces, and, for cellulose II, chain directions differ.

Classic PDB supports at most 62 distinct one-character chain IDs, 99,999 atoms,
and coordinates below 10,000 Å in the fields used here. A future mmCIF writer
can remove those serialization limits.

## Tests

```bash
python -m unittest discover -s tests -v
```

The checks include odd-DP chains, the 18-chain profile, fixed-column PDB
formatting, explicit connectivity, intact C6 chemistry, the historical CARM
regression, CHARMM-directed residue ordering, hydroxyl-H omission, and 25%
C6-oxidized output. They also verify the Paajanen eligible-site count,
chain-specific outward alternating classes and rejection of incompatible
`all` scope.

The full investigation and numerical validation are recorded in
[`docs/CHARMM_GUI_RECOGNITION_REPORT.md`](docs/CHARMM_GUI_RECOGNITION_REPORT.md).
