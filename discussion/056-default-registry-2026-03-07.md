# 056 — Default atom-type registry

Status: In Progress
Date: 2026-03-07
Relates: [217](217-rhea-participant-failures-2026-08-30.md),
[240](240-validation-constraint-ops-2026-10-18.md),
[189](189-property-test-organization-2026-08-06.md)

The [HW extraction](#hw-extraction-spike-2026-10-31) records hydrogen selections,
noble-gas evidence, halogen selections, chalcogen extraction and the
[pnictogen work table](#pnictogen-work-table). The noble-gas audit covers earlier and
newly added entries; the supplied primary-paper PDFs have been reviewed. The selected
hydrogen, noble-gas, halogen, chalcogen and pnictogen passes are complete. Unselected
candidates remain in the extraction tables. The [chalcogen work table](#chalcogen-work-table)
includes the completed audit of existing entries.
The selected pnictogen additions are implemented: 95 states, 134 fixtures and the
approved Counts bounds/table changes. The 79 earlier pnictogen entries have also
received the citation/comment audit. SbCl4−/BiCl4− use the accepted ionic/donor
description; the other extracted cases are explicitly unselected.
The March planning and implementation notes below are historical; later updates to
the tracking table are dated.

This document tracks incremental, reviewable expansion of the default atom type registry.
Registry entries are treated as code: explicit scope, per-entry review, and matching conformance coverage.

## Project Status (2026-11-02)

**Last curation check (2026-11-02):** 1041 resolution conformance cases, 403 valence
unit tests, focused Clippy and formatting passed. The seven follow-up fixtures agree
across both models/policies and satisfy the valence invariants. The preceding
90-state batch also passed workspace tests, 215 graph properties, 3150 Python tests
(2 skipped), workspace Clippy and formatting. The Counts accounting correction in
[240](240-validation-constraint-ops-2026-10-18.md#counts-electron-accounting) is implemented:
107 snapshots now give the selected donor/multicenter assignments under both policies.

Each EDN input runs through AtomTyping and Counts with both Strict and MostSaturated policies. The selected
noble-gas HW additions and six xenon-anion replacements are implemented. Perxenate
uses Xe#v8 as an accepted assumption. All Rn compound entries, neutral Xe#m8
and the XeO4 adduct entries/fixtures are removed; MeCN→XeO3 uses only the Xe=O description.
Registry references are complete for this slice: HW, selected primary sources and PW assignment markers.
All 34 retained noble-gas states and the original three F/N companion states have conformance coverage.
Hydrogen H05–H09 and the B/C/F/Cr/Rh companions are implemented, with seven new fixtures;
all nine current H states have coverage, including H#c+#m0 in HF2−/H3+ and
H#c-#m2 in [HCr2(CO)10]−. The retained halogen additions contribute 55 states,
60 new fixtures and two revised fixtures (ClF3 and IF5). All 55 added states have
conformance coverage, including the six oxyfluoride cases in five valence graphs
and H2F3−, branched H3F4− and H4F5− as F−·nHF with hydrogen bonds. The latter
replace the three aggregate-specific multicenter F states and their fixtures.
The hydrogen-bond group also includes neutral H2O, HF and NH3 dimers using existing types.
The donor–acceptor pass retains 19 states with thirteen fixtures after restoring
F#c-#n#d3 and [F(XeOF4)3]− on 2026-11-01. All 19 states are covered. Terminal F−/Cl−/Br−/I−
donors are general ligand states, tested with Rh(III) complexes using the existing Rh state.
Terminal At− donors, [At(py)2]+ and HIO3·DMSO are not selected for the default.
The retained Cl+/Br+/I+ mechanism states have three new atomic-ion fixtures.
Further extracted candidates are recorded below.
The [citation and default-scope audit](#whole-registry-audit-2026-10-31) reviewed
all 112 halogen states, including 37 pre-existing entries. Cl#c+#n3#m0 remains omitted.
The initial omission of F#c-#n#d3 was reversed when extending the general
triply bridging halide convention to the chalcogen bromide/iodide aggregates.
The chalcogen pass adds 40 chalcogen states and four N/Br−/I− companions.
The additions cover the selected halide and nitride patterns, aromatic Se4²⁺/Te4²⁺ and
multicenter O3/O3−. Zundel H5O2+ has both 3c2e and 3c4e descriptions, each with
explicit, implicit and mixed terminal hydrogens; five O states are added.
H7O3+, H9O4+, H13O6+ and H3O2− have four hydrogen-bond compatibility fixtures
using existing atom types; the H13O6+ fixture retains the existing Zundel core.
All 44 additions have AtomTyping coverage in 109 new fixtures and revised SF4.
The final citation/coverage audit adds 22 of those fixtures using existing types,
including separate aromatic and Kekulé furan inputs. All 100 retained chalcogen
states now have citations and conformance coverage under both AtomTyping policies.
O2+, O2− and S2− have paired atom-localized and bond-charge/spin fixtures.
S3− has localized and multicenter fixtures, adding one terminal S state.
SeBr6²− and TeCl6²−/TeBr6²− add two states with three trans 3c4e units.
SeX4²− (X = Cl, Br) and TeX4²− (X = Cl, Br, I) add two states with two trans
3c4e units, explicitly supported by K84/K91/P85; SeBr4²− and TeCl4²− have fixtures.
SO, SeO and TeO have triplet bond-spin fixtures using existing divalent states.
SF5 adds S#u#v#m4: one axial single bond, two basal 3c4e units and one unpaired electron.
S4N4 and Se4N4 have localized transannular-bond descriptions. S4N4 also retains
one HW resonance contributor with two 3c4e π units. The optional 12-electron
whole-cage descriptions and their three additional atom types are excluded from the default.
SN+, S2N+ and S4N2 have fixtures using existing states; S4N2 has one N–S–N 3c4e π unit.
S2N2 has aromatic and singlet-diradical multicenter fixtures, also using existing states.
Seven aromatic nitride ions and radical E3N2+ (E = S, Se) have nine fixtures;
S3N3− adds N#c-#n#v2#a2. Se4N2 is excluded following its reassignment to Se3N2Cl2.
S3N−, S4N− and SN2²− have three fixtures using existing states; the two chains
use six- and eight-electron π systems, respectively, and SN2²− uses N−=S=N−.
Localized S8²⁺/Se8²⁺ and aromatic S4²⁺/S6²⁺ reuse existing states.
S5+/S6+/S7+ have multicenter fixtures using existing S#n#v2#m2, with the accepted PW allocation.
Te6²⁺/Te6⁴⁺ have fixtures for localized resonance contributors using existing Te states (Be94/Bu82).
Singlet phenyloxenium tests existing O#c+#n2#v; H13 replaces the OH+ registry example.
TeF3+ and the restored explicit-H H3Te+ fixture cover the existing Te#c+#n#v3 state.
Their distinct experimental evidence is recorded in Ch19; no atom type is added.
Restoring F#c-#n#d3 adds one further state and restores its fixture.
One duplicate hydronium row is removed. The chalcogens now have 100 unique states;
the current halogen total is 110 after the HF-aggregate pruning.
The general aromaticity default is now permissive Hückel; named Daylight/MDL presets
retain their restrictions. The living ValenceTable adds aromatic contributions for
Sb and Te to support Te⁺ and Te in Te4²⁺.
The approved Se=14 and Te=16 electron-capacity limits and their scope are recorded
in the [expansion rationale](#electron-capacity-expansion-2026-11-01).
PoF6 is omitted: the supplied Weinstock report establishes volatility without a formula.

The pnictogen pass adds 17 N, 23 P, 16 As, 22 Sb and 17 Bi states, all covered by
134 new fixtures. Of these, 132 give identical determined graphs across all four
resolver/policy combinations; the two NH inputs with unspecified spin remain
underdetermined under Strict AtomTyping and select the existing singlet under
MostSaturated. All 534 determined outputs satisfy the valence invariants.
Eight older snapshots gain the triplet N candidate under Strict AtomTyping;
their MostSaturated outputs are unchanged. The 79 earlier pnictogen states are
retained with checked citations and corrected examples. Three stale Python
expectations are corrected for the earlier Counts donor accounting and conservative
H-count retention during Strict export. No resolver code changes.

### Conformance suite breakdown

| Category | Files | Description |
|---|---|---|
| `atoms/` | 108 | Free atoms (ground/excited states), H through Hs |
| `atomic_ions/` | 137 | Atomic ions and Hg₂²⁺ |
| `hydrides/` | 29 | Hydrides with explicit or implicit H |
| `inorganic_small/` | 103 | Small inorganic molecules and salts |
| `inorganic_ions/` | 173 | Inorganic ions |
| `donor_acceptor/` | 65 | Donor–acceptor complexes with dative or multicenter bonds |
| `hydrogen_bonds/` | 10 | Resolution compatibility with explicit hydrogen-bond relations |
| `hydrocarbons/` | 53 | C/H compounds |
| `functional_groups/` | 74 | Organic heteroatom groups (halides, amines, etc.) |
| `carbonyl_derivatives/` | 62 | Non-aromatic carbonyl compounds (aldehydes, ketones, acids, esters, amides, acyl derivatives) |
| `heterocycles/` | 53 | Non-aromatic O/N/S heterocycles, crown ethers |
| `aromatic/` | 82 | Aromatic carbocycles, heterocycles and ions; aromatic and Kekulé inputs |
| `hypervalent/` | 49 | Expanded-octet and multicenter species |
| `transition_metals/` | 7 | Transition-metal complexes |
| `stereo_tetrahedral/` | 22 | Tetrahedral stereochemistry |
| `stereo_cis_trans/` | 14 | Cis/trans stereochemistry |
| **Total** | **1041** | |

## March 2026 record

### Key code files

- **Registry**: `umol-models-graph/config/default-registry.toml` — per-element atom type spec lists
- **Valence table**: `umol-models-graph/config/default-valence-table.toml` — element `outer_electrons` + `allowed_valences`
- **Test harness**: `umol-models-graph/tests/resolution/mod.rs` — rstest-based conformance runner with insta snapshots
- **Test data**: `umol-models-graph/tests/resolution/data/<category>/<compound>.toml`
- **Snapshots**: `umol-models-graph/tests/resolution/snapshots/`
- **Resolver**: `umol-models-graph/src/graph_ir/resolver.rs`
- **Valence logic**: `umol-models-graph/src/graph_ir/valence.rs`, `config_data.rs`
- **Atom type spec/query**: `umol-models-graph/src/graph_ir/atom_type.rs`
- **Tracking doc**: `discussion/56-default-registry-2026-03-07.md` (this file)

### Conformance execution notes

- Run the conformance suite with feature flag enabled:
  - `cargo test --features conformance --test resolution`
- Keep the compile-time refresh marker comment in `umol-models-graph/tests/resolution/mod.rs` up to date when adding/removing conformance files (current marker: `v6`).
- Validation criterion for accepted conformance snapshots:
  - `atom_typing.success: true`
  - `counts.success: true`

## Tracking Table

Pnictogen review and chalcogen/N/F companion rows updated 2026-11-02; Br/I/Xe companion rows updated
2026-11-01; hydrogen, other noble-gas and other halogen rows updated 2026-10-31.
Other rows retain the March record.

Status values:

- `todo`: not started
- `in_progress`: currently being edited
- `review`: implemented, awaiting chemistry review
- `done`: accepted

| Element | Ground state | Discrete ions | Typical valences | Conformance added | Status | Notes |
|---|---|---|---|---|---|---|
| H | yes | yes (+1, −1) | yes (selected) | yes (all nine states; see [hydrogen fixtures](#conformance-coverage)) | done | H05–H09 implemented; selected hydrogen scope complete |
| He | yes | no | n/a | yes (`atoms/he.edn`) | done | HW review complete; zero-lone-pair description retained |
| Li | yes | yes (+1) | no | yes (`atoms/li.toml`, `atomic_ions/li+1.toml`) | done |  |
| Be | yes | yes (+2) | no | yes (`atoms/be.toml`, `atomic_ions/be+2.toml`) | done |  |
| B | yes | no | partial | yes (`atoms/b.edn`, `hydrides/bh3.edn`, `inorganic_small/b2h6.edn`) | done | Hydrogen companion addition: B#v2#m |
| C | yes | no | partial | yes (`atoms/c.edn`, `donor_acceptor/hcr2-co10-.edn`) | done | Hydrogen companion addition: C#c-#v3#d |
| N | yes | yes (-3) | partial | yes (`atoms/n.edn`, `atomic_ions/n-3.edn`, `donor_acceptor/`, `inorganic_small/`, `aromatic/`) | in_progress | 52 current states; 17 additions and both NO− fixtures implemented; existing-entry citation audit complete |
| O | yes | yes (-2) | yes (selected) | yes (all 29 states) | done | Thirteen states added; duplicate hydronium removed; citation and default-scope audit complete; aromatic and Kekulé furan covered |
| F | yes | yes (-1) | partial | yes (selected HW states) | done | Three higher HF-aggregate states removed; H-bond fixtures use existing types. Triply bridging F− restored 2026-11-01; terminal donor covered by [RhF6]3− |
| Ne | yes | no | n/a | yes (`atoms/ne.edn`) | done | HW review complete; zero-lone-pair description retained |
| Na | yes | yes (+1) | no | yes (`atoms/na.toml`, `atomic_ions/na+1.toml`) | done |  |
| Mg | yes | yes (+2) | no | yes (`atoms/mg.toml`, `atomic_ions/mg+2.toml`) | done |  |
| Al | yes | yes (+3) | no | yes (`atoms/al.toml`, `atomic_ions/al+3.toml`) | done |  |
| Si | yes | no | no | yes (`atoms/si.toml`) | done |  |
| P | yes | yes (-3) | partial | yes (`atoms/p.edn`, `atomic_ions/p-3.edn`, `hydrides/`, `hypervalent/`) | in_progress | 44 current states; 23 additions implemented; existing-entry citation audit complete |
| S | yes | yes (-2) | yes (selected) | yes (all 22 states) | done | Five states added; citation and default-scope audit complete; S4N4 retains localized and HW three-center descriptions |
| Cl | yes | yes (+1, -1) | partial | yes (selected HW states) | done | 16 states added; terminal donor covered by [RhCl(NH3)5]2+ |
| Ar | yes | no | n/a | yes (`atoms/ar.edn`) | done | HW review complete; zero-lone-pair description retained |
| K | yes | yes (+1) | no | yes (`atoms/k.toml`, `atomic_ions/k+1.toml`) | done |  |
| Ca | yes | yes (+2) | no | yes (`atoms/ca.toml`, `atomic_ions/ca+2.toml`) | done |  |
| Sc | yes | yes (+3) | no | yes (`atoms/sc.toml`, `atomic_ions/sc+3.toml`) | done |  |
| Ti | yes | yes (+3,+4) | no | yes (`atoms/ti.toml`, `atomic_ions/ti+3.toml`, `atomic_ions/ti+4.toml`) | done |  |
| V | yes | yes (+2,+3) | no | yes (`atoms/v.toml`, `atomic_ions/v+2.toml`, `atomic_ions/v+3.toml`) | done |  |
| Cr | yes | yes (+2,+3) | partial | yes (`atoms/cr.edn`, `atomic_ions/cr+2.edn`, `atomic_ions/cr+3.edn`, `donor_acceptor/hcr2-co10-.edn`) | done | Hydrogen companion addition: Cr#n3#t5#m0 |
| Mn | yes | yes (+2,+3) | no | yes (`atoms/mn.toml`, `atomic_ions/mn+2.toml`, `atomic_ions/mn+3.toml`) | done |  |
| Fe | yes | yes (+2,+3) | no | yes (`atoms/fe.toml`, `atomic_ions/fe+2.toml`, `atomic_ions/fe+3.toml`) | done |  |
| Co | yes | yes (+2,+3) | no | yes (`atoms/co.toml`, `atomic_ions/co+2.toml`, `atomic_ions/co+3.toml`) | done |  |
| Ni | yes | yes (+2) | no | yes (`atoms/ni.toml`, `atomic_ions/ni+2.toml`) | done |  |
| Cu | yes | yes (+1,+2) | no | yes (`atoms/cu.toml`, `atomic_ions/cu+1.toml`, `atomic_ions/cu+2.toml`) | done |  |
| Zn | yes | yes (+2) | no | yes (`atoms/zn.toml`, `atomic_ions/zn+2.toml`) | done |  |
| Ga | yes | yes (+3) | no | yes (`atoms/ga.toml`, `atomic_ions/ga+3.toml`) | done |  |
| Ge | yes | no | no | yes (`atoms/ge.toml`) | done |  |
| As | yes | yes (-3) | partial | yes (`atoms/as.edn`, `atomic_ions/as-3.edn`, `hydrides/`, `inorganic_ions/`) | in_progress | 27 current states; 15 additions implemented; existing arsole types retained with paired aromatic/Kekulé fixtures; citation audit complete |
| Se | yes | yes (-2) | yes (selected) | yes (all 22 states) | done | Six states added; citation and default-scope audit complete, including localized Se4N4 |
| Br | yes | yes (+1, -1) | partial | yes (selected HW states) | done | 18 halogen-pass additions plus Br#c-#n#d3 for binary chalcogen bromides |
| Kr | yes | no | yes (selected) | yes (`atoms/kr.edn`, `inorganic_small/krf2.edn`, `inorganic_ions/krf+.edn`, `donor_acceptor/hcn-krf+.edn`) | done | HW review complete; all four retained states covered |
| Rb | yes | yes (+1) | no | yes (`atoms/rb.toml`, `atomic_ions/rb+1.toml`) | done |  |
| Sr | yes | yes (+2) | no | yes (`atoms/sr.toml`, `atomic_ions/sr+2.toml`) | done |  |
| Y | yes | yes (+3) | no | yes (`atoms/y.toml`, `atomic_ions/y+3.toml`) | done |  |
| Zr | yes | yes (+4) | no | yes (`atoms/zr.toml`, `atomic_ions/zr+4.toml`) | done |  |
| Nb | yes | yes (+3) | no | yes (`atoms/nb.toml`, `atomic_ions/nb+3.toml`) | done |  |
| Mo | yes | yes (+2,+3) | no | yes (`atoms/mo.toml`, `atomic_ions/mo+2.toml`, `atomic_ions/mo+3.toml`) | done |  |
| Tc | yes | no | no | yes (`atoms/tc.toml`) | done | No ion added in current scope |
| Ru | yes | yes (+2,+3) | no | yes (`atoms/ru.toml`, `atomic_ions/ru+2.toml`, `atomic_ions/ru+3.toml`) | done |  |
| Rh | yes | yes (+2,+3) | partial | yes (`atoms/rh.edn`, `atomic_ions/rh+2.edn`, `atomic_ions/rh+3.edn`, Rh–H donor fixture below) | done | Hydrogen companion addition: Rh#c+3#n3#t6 |
| Pd | yes | yes (+2) | no | yes (`atoms/pd.toml`, `atomic_ions/pd+2.toml`) | done |  |
| Ag | yes | yes (+1,+2) | partial | yes (`atoms/ag.edn`, `atomic_ions/ag+1.edn`, `atomic_ions/ag+2.edn`, `donor_acceptor/ag2-i2-4-2+.edn`) | done | I2-ligand companion: Ag#c+#n5#t3 |
| Cd | yes | yes (+2) | no | yes (`atoms/cd.toml`, `atomic_ions/cd+2.toml`) | done |  |
| In | yes | yes (+3) | no | yes (`atoms/in.toml`, `atomic_ions/in+3.toml`) | done |  |
| Sn | yes | no | no | yes (`atoms/sn.toml`) | done |  |
| Sb | yes | yes (-3) | partial | yes (`atoms/sb.edn`, `atomic_ions/sb-3.edn`, `hydrides/`, `inorganic_ions/`) | in_progress | 29 current states; 20 additions and capacity 12 → 14 implemented; existing-entry citation audit complete |
| Te | yes | yes (-2, -1, +1) | yes (selected) | yes (all 20 states) | done | Eleven states added, including aromatic Te4²⁺; citation and default-scope audit complete |
| I | yes | yes (+1, -1) | partial | yes (selected HW states) | done | 22 halogen-pass additions plus I#c-#n#d3 for binary TeI4 aggregates |
| Xe | yes | no | yes (selected) | yes (`atoms/xe.edn`, `inorganic_small/xe*.edn`, `inorganic_ions/xe*.edn`, `donor_acceptor/`) | done | All 26 retained states covered after restoring [F(XeOF4)3]− |
| Cs | yes | yes (+1) | no | yes (`atoms/cs.toml`, `atomic_ions/cs+1.toml`) | done |  |
| Ba | yes | yes (+2) | no | yes (`atoms/ba.toml`, `atomic_ions/ba+2.toml`) | done |  |
| Tl | yes | yes (+1,+3) | no | yes (`atoms/tl.toml`, `atomic_ions/tl+1.toml`, `atomic_ions/tl+3.toml`) | done |  |
| Pb | yes | yes (+2,+4) | no | yes (`atoms/pb.toml`, `atomic_ions/pb+2.toml`, `atomic_ions/pb+4.toml`) | done |  |
| Bi | yes | yes (+3, −3) | partial | yes (`atoms/bi.edn`, `atomic_ions/`, `hydrides/`, `inorganic_ions/`, `aromatic/`) | in_progress | 17 current states; 15 additions and capacity 12 → 14 implemented; existing-entry citation audit complete |
| Po | yes | yes (+2, +4, −2) | yes (selected) | yes (all seven states) | done | Five states added; citation and default-scope audit complete; PoF6 omitted for lack of composition evidence |
| At | yes | yes (+1, -1) | partial | yes (selected HW states) | done | Nine states added; terminal At− donor and [At(py)2]+ not selected |
| Rn | yes | no | none selected | yes (`atoms/rn.edn`) | done | HW review complete; unsupported compound states removed; Rn#n4 retained |
| La | yes | yes (+3) | no | yes (`atoms/la.toml`, `atomic_ions/la+3.toml`) | done |  |
| Ce | yes | yes (+3,+4) | no | yes (`atoms/ce.toml`, `atomic_ions/ce+3.toml`, `atomic_ions/ce+4.toml`) | done |  |
| Pr | yes | yes (+3,+4) | no | yes (`atoms/pr.toml`, `atomic_ions/pr+3.toml`, `atomic_ions/pr+4.toml`) | done | Counts requires `Pr` valence-table support for +4 |
| Nd | yes | yes (+3) | no | yes (`atoms/nd.toml`, `atomic_ions/nd+3.toml`) | done |  |
| Pm | yes | yes (+3) | no | yes (`atoms/pm.toml`, `atomic_ions/pm+3.toml`) | done |  |
| Sm | yes | yes (+2,+3) | no | yes (`atoms/sm.toml`, `atomic_ions/sm+2.toml`, `atomic_ions/sm+3.toml`) | done |  |
| Eu | yes | yes (+2,+3) | no | yes (`atoms/eu.toml`, `atomic_ions/eu+2.toml`, `atomic_ions/eu+3.toml`) | done |  |
| Gd | yes | yes (+3) | no | yes (`atoms/gd.toml`, `atomic_ions/gd+3.toml`) | done |  |
| Tb | yes | yes (+3,+4) | no | yes (`atoms/tb.toml`, `atomic_ions/tb+3.toml`, `atomic_ions/tb+4.toml`) | done |  |
| Dy | yes | yes (+3) | no | yes (`atoms/dy.toml`, `atomic_ions/dy+3.toml`) | done |  |
| Ho | yes | yes (+3) | no | yes (`atoms/ho.toml`, `atomic_ions/ho+3.toml`) | done |  |
| Er | yes | yes (+3) | no | yes (`atoms/er.toml`, `atomic_ions/er+3.toml`) | done |  |
| Tm | yes | yes (+3) | no | yes (`atoms/tm.toml`, `atomic_ions/tm+3.toml`) | done |  |
| Yb | yes | yes (+2,+3) | no | yes (`atoms/yb.toml`, `atomic_ions/yb+2.toml`, `atomic_ions/yb+3.toml`) | done |  |
| Lu | yes | yes (+3) | no | yes (`atoms/lu.toml`, `atomic_ions/lu+3.toml`) | done |  |
| Hf | yes | yes (+4) | no | yes (`atoms/hf.toml`, `atomic_ions/hf+4.toml`) | done |  |
| Ta | yes | no | no | yes (`atoms/ta.toml`) | done |  |
| W | yes | no | no | yes (`atoms/w.toml`) | done |  |
| Re | yes | no | no | yes (`atoms/re.toml`) | done |  |
| Os | yes | no | no | yes (`atoms/os.toml`) | done |  |
| Ir | yes | no | no | yes (`atoms/ir.toml`) | done |  |
| Pt | yes | yes (+2) | no | yes (`atoms/pt.toml`, `atomic_ions/pt+2.toml`) | done |  |
| Au | yes | yes (+1,+3) | no | yes (`atoms/au.toml`, `atomic_ions/au+1.toml`, `atomic_ions/au+3.toml`) | done | `Au3+` conformance pinned to low-spin |
| Hg | yes | yes (+2, [Hg2]2+) | no | yes (`atoms/hg.toml`, `atomic_ions/hg+2.toml`, `atomic_ions/hg+2_dimer.toml`) | done | Dimer encoded as two `Hg+1` atoms with single bond |
| Fr | yes | yes (+1) | no | yes (`atoms/fr.toml`, `atomic_ions/fr+1.toml`) | done |  |
| Ra | yes | yes (+2) | no | yes (`atoms/ra.toml`, `atomic_ions/ra+2.toml`) | done |  |
| Ac | yes | yes (+3) | no | yes (`atoms/ac.toml`, `atomic_ions/ac+3.toml`) | done |  |
| Th | yes | yes (+3,+4) | no | yes (`atoms/th.toml`, `atomic_ions/th+3.toml`, `atomic_ions/th+4.toml`) | done |  |
| Pa | yes | yes (+4,+5) | no | yes (`atoms/pa.toml`, `atomic_ions/pa+4.toml`, `atomic_ions/pa+5.toml`) | done | Valence table outer electrons adjusted to 5 |
| U | yes | yes (+3,+4,+5,+6) | no | yes (`atoms/u.toml`, `atomic_ions/u+3.toml`, `atomic_ions/u+4.toml`, `atomic_ions/u+5.toml`, `atomic_ions/u+6.toml`) | done | Valence table outer electrons adjusted to 6 |
| Np | yes | yes (+3,+4,+5,+6) | no | yes (`atoms/np.toml`, `atomic_ions/np+3.toml`, `atomic_ions/np+4.toml`, `atomic_ions/np+5.toml`, `atomic_ions/np+6.toml`) | done | Valence table outer electrons adjusted to 6 |
| Pu | yes | yes (+3,+4,+5,+6,+7) | no | yes (`atoms/pu.toml`, `atomic_ions/pu+3.toml`, `atomic_ions/pu+4.toml`, `atomic_ions/pu+5.toml`, `atomic_ions/pu+6.toml`, `atomic_ions/pu+7.toml`) | done | Valence table outer electrons adjusted to 7 |
| Am | yes | yes (+2,+3,+4) | no | yes (`atoms/am.toml`, `atomic_ions/am+2.toml`, `atomic_ions/am+3.toml`, `atomic_ions/am+4.toml`) | done |  |
| Cm | yes | yes (+3,+4) | no | yes (`atoms/cm.toml`, `atomic_ions/cm+3.toml`, `atomic_ions/cm+4.toml`) | done |  |
| Bk | yes | yes (+3,+4) | no | yes (`atoms/bk.toml`, `atomic_ions/bk+3.toml`, `atomic_ions/bk+4.toml`) | done |  |
| Cf | yes | yes (+2,+3,+4) | no | yes (`atoms/cf.toml`, `atomic_ions/cf+2.toml`, `atomic_ions/cf+3.toml`, `atomic_ions/cf+4.toml`) | done |  |
| Es | yes | yes (+2,+3) | no | yes (`atoms/es.toml`, `atomic_ions/es+2.toml`, `atomic_ions/es+3.toml`) | done |  |
| Fm | yes | yes (+2,+3) | no | yes (`atoms/fm.toml`, `atomic_ions/fm+2.toml`, `atomic_ions/fm+3.toml`) | done |  |
| Md | yes | yes (+2,+3) | no | yes (`atoms/md.toml`, `atomic_ions/md+2.toml`, `atomic_ions/md+3.toml`) | done |  |
| No | yes | yes (+2) | no | yes (`atoms/no.toml`, `atomic_ions/no+2.toml`) | done |  |
| Lr | yes | yes (+3) | no | yes (`atoms/lr.toml`, `atomic_ions/lr+3.toml`) | done |  |
| Rf | yes | no | no | yes (`atoms/rf.toml`) | done |  |
| Db | yes | no | no | yes (`atoms/db.toml`) | done |  |
| Sg | yes | no | no | yes (`atoms/sg.toml`) | done |  |
| Bh | yes | no | no | yes (`atoms/bh.toml`) | done |  |
| Hs | yes | no | no | yes (`atoms/hs.toml`) | done |  |

Extend this table as new elements are considered.

### Organic compounds

| Class | Scope | Conformance dir | Status | Notes |
|-------|-------|-----------------|--------|-------|
| Hydrocarbons | Alkanes, alkenes, alkynes, dienes/cumulenes, cycloalkanes, unsaturated alicycles, fused/spiro/bridged bicycles, polycycles | `hydrocarbons/` | done | 51 files, `implicit_h=true`, `?{CHn}` |
| Alcohols / ethers | R–OH, R–O–R′ | — | todo | |
| Carbonyl derivatives | Aldehydes, ketones, carboxylic acids, esters, amides | `carbonyl_derivatives/` | review | 62 files, `implicit_h=true`, non-aromatic |
| Organic halides | R–F, R–Cl, R–Br, R–I | — | todo | |
| Other | Amines, thiols, etc. | — | todo | |
| Functional groups | Halides, azides, cyanates, thiocyanates, alcohols, ethers, amines, nitroso, nitro, sulfoxides, hydrazines, diazo/azo | `functional_groups/` | review | 73 files, `implicit_h=true` |
| Heterocycles | Non-aromatic O/N/S heterocycles, crown ethers, bridged/fused | `heterocycles/` | review | 52 files, `implicit_h=true` |

Extend this table as new organic classes are added.

## Conformance Policy

Per element, add explicit resolution inputs that cover:

- Free atom query (with explicit `implicit_h` policy)
- At least one bonded case when chemically meaningful
- Relevant charge states (for ion stage)
- At least one "should fail" case where useful to constrain over-broad matching

Suggested naming pattern:

- Atomic ground state: `tests/resolution/data/atoms/<element>.toml`
- Atomic valence/excited variants: `tests/resolution/data/atoms/<element>_<motif>.toml`
- Ionic ground state: `tests/resolution/data/atomic_ions/<element><charge>.toml`
- Ionic valence/excited variants: `tests/resolution/data/atomic_ions/<element><charge>_<motif>.toml`

Conventions:

- `<element>`: lowercase element symbol/name token used in file naming.
- `<charge>`: signed integer with explicit sign, e.g. `+1`, `-2` (no underscore between element and charge).
- `<motif>`: lowercase, underscore-separated descriptor; may encode occupation and optional multiplet.

## Categories

The [conformance suite breakdown](#conformance-suite-breakdown) lists the current
categories. Donor–acceptor adducts use `donor_acceptor/`, including charged adducts.
The `hydrogen_bonds/` fixtures check resolution of graphs carrying explicit Hbd
relations, including their preservation alongside localized and multicenter bonds.

The former `implicit_h=true/false` setting was the precursor to atom defaults.
EDN inputs use MoleculeDefaults with per-input `:config-overrides`; hydrogen
handling follows those defaults and explicit atom predicates.

## Current Conformance Set

- `atoms/`: `h`, `he`, `li`, `be`, `b`, `c`, `n`, `o`, `f`, `ne`, `na`, `mg`, `al`, `si`, `p`, `s`, `cl`, `ar`, `k`, `ca`, `sc`, `ti`, `v`, `cr`, `mn`, `fe`, `co`, `ni`, `cu`, `zn`, `ga`, `ge`, `as`, `se`, `br`, `kr`, `rb`, `sr`, `y`, `zr`, `nb`, `mo`, `tc`, `ru`, `rh`, `pd`, `ag`, `cd`, `in`, `sn`, `sb`, `te`, `i`, `xe`, `cs`, `ba`, `tl`, `pb`, `bi`, `po`, `at`, `rn`, `la`, `ce`, `pr`, `nd`, `pm`, `sm`, `eu`, `gd`, `tb`, `dy`, `ho`, `er`, `tm`, `yb`, `lu`, `hf`, `ta`, `w`, `re`, `os`, `ir`, `pt`, `au`, `hg`, `fr`, `ra`, `ac`, `th`, `pa`, `u`, `np`, `pu`, `am`, `cm`, `bk`, `cf`, `es`, `fm`, `md`, `no`, `lr`, `rf`, `db`, `sg`, `bh`, `hs`
- `atomic_ions/`: `h+1`, `li+1`, `be+2`, `mg+2`, `al+3`, `n-3`, `o-2`, `f-1`, `na+1`, `p-3`, `s-2`, `cl-1`, `k+1`, `ca+2`, `sc+3`, `ti+3`, `ti+4`, `v+2`, `v+3`, `cr+2`, `cr+3`, `mn+2`, `mn+3`, `fe+2`, `fe+3`, `co+2`, `co+3`, `ni+2`, `cu+1`, `cu+2`, `zn+2`, `ga+3`, `as-3`, `se-2`, `br-1`, `rb+1`, `sr+2`, `y+3`, `zr+4`, `nb+3`, `mo+2`, `mo+3`, `ru+2`, `ru+3`, `rh+2`, `rh+3`, `pd+2`, `ag+1`, `ag+2`, `cd+2`, `in+3`, `sb-3`, `te-2`, `i-1`, `cs+1`, `ba+2`, `tl+1`, `tl+3`, `pb+2`, `pb+4`, `bi+3`, `po-2`, `at-1`, `la+3`, `ce+3`, `ce+4`, `pr+3`, `pr+4`, `nd+3`, `pm+3`, `sm+2`, `sm+3`, `eu+2`, `eu+3`, `gd+3`, `tb+3`, `tb+4`, `dy+3`, `ho+3`, `er+3`, `tm+3`, `yb+2`, `yb+3`, `lu+3`, `hf+4`, `pt+2`, `au+1`, `au+3`, `hg+2`, `hg+2_dimer`, `fr+1`, `ra+2`, `ac+3`, `th+3`, `th+4`, `pa+4`, `pa+5`, `u+3`, `u+4`, `u+5`, `u+6`, `np+3`, `np+4`, `np+5`, `np+6`, `pu+3`, `pu+4`, `pu+5`, `pu+6`, `pu+7`, `am+2`, `am+3`, `am+4`, `cm+3`, `cm+4`, `bk+3`, `bk+4`, `cf+2`, `cf+3`, `cf+4`, `es+2`, `es+3`, `fm+2`, `fm+3`, `md+2`, `md+3`, `no+2`, `lr+3`
- `hydrides/`: `h2`, `bh3`, `ch4`, `nh3`, `h2o`, `hf`, `sih4`, `ph3`, `h2s`, `hcl`, `ash3`, `h2se`, `hbr`, `sbh3`, `h2te`, `hi`
- `inorganic_small/`: neutral small inorganics (explicit H only): e.g. `bf2`, `bf3`, `bn`, `c2`, `ch3`, `ch2_singlet`, `ch2_triplet`, `sih3`, `sih2_*`, `nh2`, `ph2`, `ash2`, `sbh2`, `oh`, `sh`, `seh`, `teh`, `br2`, `cl2`, `cn-`, `f2`, `i2`, `n2`, `o2`, etc. Added nitrogen oxides (NO, NO2, N2O, N2O3, N2O4), carbon compounds (CO2, CS2, COS, HCN), oxyacids/hydroxides (HNO2, HNO3, HClO, H2CO3, H3BO3, HOBr, HOI), N-H compounds (N2H4, NH2OH, N2H2), peroxides (H2O2, H2S2), halides (NF3, NCl3, PF3, PCl3, SF2, SiF4, SiCl4, BCl3, ClF, BrF), and other (O3, HNCO, HOCN, HSCN, HNCS).
- `inorganic_ions/`: molecular ions (explicit H only): e.g. `bcl4-`, `cf3-`, `ccl3+`, `nh2-`, `nh4+`, `ph2-`, `ph4+`, `ash2-`, `ash4+`, `sbh2-`, `sbh4+`, `oh-`, `oh3+`, `sh-`, `sh3+`, `seh-`, `seh3+`, `teh-`, `teh3+`, `f-`, `cl-`, `br-`, `i-`, `sif3-`, `sicl3+`, etc.
- `hydrocarbons/`: `ethane`, `propane`, `butane`, `isobutane`, `pentane`, `isopentane`, `neopentane`, `hexane`, `ethene`, `propene`, `but-1-ene`, `but-2-ene`, `2-methylpropene`, `2-methylbut-2-ene`, `ethyne`, `propyne`, `but-1-yne`, `but-2-yne`, `allene`, `buta-1,3-diene`, `penta-1,4-diene`, `butatriene`, `cyclopropane`, `cyclobutane`, `cyclopentane`, `cyclohexane`, `cycloheptane`, `cyclopropene`, `cyclobutene`, `cyclopentene`, `cyclopentadiene`, `cyclohexene`, `cyclohexa-1,3-diene`, `cyclohexa-1,4-diene`, `bicyclo-1.1.0-butane`, `bicyclo-2.1.0-pentane`, `bicyclo-2.2.0-hexane`, `bicyclo-3.3.0-octane`, `bicyclo-4.3.0-nonane`, `decalin`, `spiropentane`, `spirohexane`, `spiroheptane`, `spirononane`, `bicyclo-1.1.1-pentane`, `norbornane`, `bicyclo-2.2.2-octane`, `tetrahedrane`, `prismane`, `cubane`, `adamantane`
- `functional_groups/`: methyl/ethyl × {fluoride, chloride, bromide, iodide, azide, isocyanide, cyanate, thiocyanate, isothiocyanate, nitrite, nitrate, thiol, selenol} (26); {propyl, isopropyl, butyl, isobutyl, tert-butyl, cyclopropyl, cyclohexyl}-chloride (7); dichloromethane, chloroform, carbon-tetrachloride, 1,1-dichloroethane, 1,2-dichloroethane, 1,2-dichlorocyclohexane, 1-chloro-2-hydroxyethane (7); methanol, ethanol, propan-1-ol, propan-2-ol, butan-1-ol, butan-2-ol, 2,2-dimethylpropan-2-ol, hydroxycyclopropane, hydroxycyclohexane (9); dimethylether, diethylether, 1,2-dimethoxyethane (3); methylamine, ethylamine, dimethylamine, trimethylamine, tetramethylammonium, trimethylamine-oxide (6); nitrosomethane, nitrosoethane, nitromethane, nitroethane, nitrosocyclohexane, nitrocyclohexane (6); dimethyl-sulfoxide, methylhydrazine, 1,1-dimethylhydrazine, 1,2-dimethylhydrazine, diazomethane, dimethyldiazene (6). Total: 73 files, `implicit_h=true`.
- `carbonyl_derivatives/`: non-aromatic carbonyl compounds (all `implicit_h=true`): aldehydes, ketones, dicarbonyls/hydroxy-carbonyls, monocarboxylic and dicarboxylic acids, acyl halides, anhydrides, carboxylates, esters/lactones/thioesters, amides/lactams/imides, and carbonyl-derived C=N species (oxime/hydrazone), plus pyruvic acid. Total: 62 files.
- `hypervalent/`: expanded-octet species (all `implicit_h=false`): SO2, SO3, H2SO3, H2SO4, SF4, SF6, H3PO3, H3PO4, P4O10, PF5, PCl5, POCl3, HClO2, HClO3, HClO4, ClF3, IF5, SO4(2-), PO4(3-), ClO4(-). Total: 20 files.
- `heterocycles/`: non-aromatic heterocycles (all `implicit_h=true`). O-saturated (oxirane, dioxirane, oxetane, 1,2-dioxetane, 1,3-dioxetane, tetrahydrofuran, tetrahydropyran); N-saturated (aziridine, azetidine, pyrrolidine, piperidine); S-saturated (thiirane, thietane, tetrahydrothiophene, thiane); O/S-unsaturated (2H-oxete, 2H-thiete, 2H-pyran, 4H-pyran, 2H-thiopyran, 4H-thiopyran); N-unsaturated (2H-azirine, 2,3-dihydroazete, azete, 3-pyrroline, 2-pyrroline); functionalized (sulfolane, 2-oxazolidone, hydantoin); multi-O/S (1,3-dioxolane, 1,2,3-trioxolane, 1,4-dioxane, 1,3-dithiane, 1,4-dithiane, 1,3,5-trithiane); multi-N (1,3-diazetidine, pyrazolidine, imidazolidine, 2-pyrazoline, 2-imidazoline); mixed (1,2-oxazetidine, 1,2-oxathiolane, 1,3-oxathiolane, morpholine, thiomorpholine); bridged/fused (pyrrolizidine, quinuclidine, 1-azaadamantane); crown ethers (12-crown-4, 15-crown-5, 18-crown-6). Total: 52 files.

## Decision Log

### 2026-03-07

- Agreed staged rollout: ground states -> discrete ions -> typical non-metal valences.
- Agreed registry expansion must be paired with conformance additions.
- Agreed registry specs are reviewed with code-level rigor (small, explicit, traceable changes).
- Agreed conformance naming/categories: `atoms/`, `atomic_ions/`, `hydrides/` and file naming rules above.
- Added `implicit_h` control to conformance inputs; counts strategy now honors this via `enable_implicit_hydrogens`.
- Expanded coverage to period-3 atoms (`Na`..`Ar`) and matching ions (`Na+`, `Al3+`, `P3-`, `S2-`, `Cl-`), plus `Be2+`.
- Expanded coverage to period-4 and period-5 s/p blocks, 3d block, and 4d block with ions: `Y3+`, `Zr4+`, `Nb3+`, `Mo2+`, `Mo3+`, `Ru2+`, `Ru3+`, `Rh2+`, `Rh3+`, `Pd2+`, `Ag+`, `Ag2+`, `Cd2+`.
- Added missing 3d atom conformance files (`atoms/sc.toml` through `atoms/zn.toml`) and refreshed snapshots.
- Added 6s/6p conformance coverage: atoms `Cs`, `Ba`, `Tl`, `Pb`, `Bi`, `Po`, `At`, `Rn`; ions `Cs+`, `Ba2+`, `Tl+`, `Tl3+`, `Pb2+`, `Pb4+`, `Bi3+`, `Po2-`, `At-`.
- Added 4f conformance coverage (`La`..`Lu`) with ions: all `3+`, plus `Ce4+`, `Pr4+`, `Tb4+`, and `Sm2+`, `Eu2+`, `Yb2+`.
- Added 5d conformance coverage (`Hf`..`Hg`) with ions: `Hf4+`, `Pt2+`, `Au+`, `Au3+`, `Hg2+`, and `[Hg2]2+` as `atomic_ions/hg+2_dimer.toml`.
- Pinned spin-sensitive transition-metal ion queries in conformance inputs (high-spin/low-spin where requested) and aligned `Pr` counts support (`default-valence-table`: `outer_electrons=5`, `allowed_valences=[3,4]`) so all conformance snapshots resolve with success.
- Added conformance coverage for `Fr`..`Lr` atoms and all corresponding registry ion charge states.
- Added neutral-atom conformance coverage for `Rf`, `Db`, `Sg`, `Bh`, `Hs`.
- Updated actinide counts input defaults by setting valence-table `outer_electrons`: `Pa=5`, `U=6`, `Np=6`, `Pu=7`.
- Added simple non-metal hydrides as atom-like conformance inputs (`?{EHn}` / `?{EH}`) with `implicit_h=true` and no explicit bonds.
- Added matching hydride atom-type specs to the default registry for atom-typing validation (`BH3`, `CH4`, `NH3`, `OH2`, `FH`, `SiH4`, `PH3`, `SH2`, `ClH`, `AsH3`, `SeH2`, `BrH`, `SbH3`, `TeH2`, `IH`).
- Counts-based typing now uses explicit bond-order sum for `v` (implicit hydrogens no longer inflate valence display), resolving outputs like `{O/2H2}` instead of `{O/2H2v2}`.
- Atoms/ions/hydrides rollout pass completed; table statuses are all `done`.
- Added covalent-state rollout phase in `covalent/` with `implicit_h=false` for all cases (all hydrogens explicit as atoms/bonds), following Hill+charge file naming (`CH3-.toml`, `NH4+.toml`) and `_isomer` suffix for spin-isomer cases (`CH2_singlet`, `CH2_triplet`).
- Covered required non-metal motif families with mixed substituents (`H`, `F`, `Cl`) across B/C/Si, pnictogens, chalcogens, and halogens, including radicals, anions, and cations.
- Added only required registry states for covalent motifs (`C/1v2`, `Si+v3`, `Si-/1v3`, `Si/1^2v2`, `P/1^1v2`, `P-/2v2`, `As/1^1v2`, `As+v4`, `As-/2v2`, `Sb/1^1v2`, `Sb+v4`, `Sb-/2v2`, `Se+/1v3`, `Te/2^1v1`, `Te+/1v3`, `Te-/3v1`).
- Conformance parity verified: `atom_typing` and `counts` both succeed on all new `covalent/` snapshots; no valence-table adjustments were required in this phase.
- Mixed implicit/explicit hydrogen representations are deferred to a future phase.
- Added `hydrocarbons/` conformance category (51 files) with `implicit_h=true` using `?{CHn}` atom queries and C-C bonds only. Covers alkanes (8), alkenes (6), alkynes (4), dienes/cumulenes (4), cycloalkanes (5), unsaturated alicycles (7), fused bicycles (6), spiro bicycles (4), bridged bicycles (3), and polycycles (4: tetrahedrane, prismane, cubane, adamantane). No new registry entries required. All 51 tests pass for both `atom_typing` and `counts` strategies.
- Restructured conformance categories: renamed `ions/` → `atomic_ions/` (single-atom ions). Split `covalent/` → `inorganic_small/` (neutral small inorganics, including sextet species) and `inorganic_ions/` (molecular ions). Policy: `hydrides/` use implicit H only (one-atom graphs); `inorganic_small/`, `inorganic_ions/`, and hydrocarbons use explicit H where present (hydrocarbons keep implicit H in queries). Compromise between fine-grained compound classes and a single mixed category.
- Added `functional_groups/` conformance category (64 files) with `implicit_h=true`. Covers organic halides (methyl/ethyl × 4 halogens + 7 alkyl chlorides + 7 poly/mixed chloro compounds = 39), azides, isocyanides, cyanates, thiocyanates, isothiocyanates, nitrites, nitrates, thiols, selenols (methyl/ethyl × 9 = 18), alcohols (9), ethers (3), amines (6), nitroso (3), nitro (3). Inputs reviewed: formal charges on N+/O-/C-/N- where required (nitro, nitrate, azide, isocyanide, TMAO, tetramethylammonium). All 64 counts-strategy results pass. 26 atom-typing results show `ValenceAmbiguous` (registry disambiguation work remaining).
- Fixed `compute_implicit_hydrogens` to use RDKit-style effective-element lookup via `element.shift(-charge)` for allowed-valence comparison. Previously, charged atoms with specific allowed valences (e.g. N+ needing C's `[4]` instead of N's `[3]`) would fail to match.
- Extended Cl `allowed_valences` to `[1, 3, 5, 7]` in `default-valence-table.toml` to support hypervalent chlorine compounds (HClO2, HClO3, HClO4, ClF3).
- Added 36 new `inorganic_small/` conformance files: nitrogen oxides (NO, NO2, N2O, N2O3, N2O4), carbon compounds (CO2, CS2, COS, HCN), oxyacids/hydroxides (HNO2, HNO3, HClO, H2CO3, H3BO3, HOBr, HOI), N-H compounds (N2H4, NH2OH, N2H2), peroxides (H2O2, H2S2), halides (NF3, NCl3, PF3, PCl3, SF2, SiF4, SiCl4, BCl3, ClF, BrF), other (O3, HNCO, HOCN, HSCN, HNCS). All `implicit_h=false` with explicit H atoms and formal charges.
- Added 17 new `hypervalent/` conformance files: sulfur (SO2, SO3, H2SO3, H2SO4, SF4, SF6), phosphorus (H3PO3, H3PO4, P4O10, PF5, PCl5, POCl3), chlorine (HClO2, HClO3, HClO4, ClF3), iodine (IF5). All `implicit_h=false`.
- Resolution results: all 53 new files resolve successfully with `counts` strategy. `atom_typing` fails on 3 files (NO2, IF5, HClO4) due to missing registry entries for radical N (`{N/1^1v3}`), hypervalent I (`{I/1v5}`), and heptavalent Cl (`{Clv7}`). N2O4 input was corrected to symmetric O=N+(O-)–N+(O-)=O structure; both strategies now succeed.
- Identified 5 new atom type specs needed for full `atom_typing` coverage of new compounds: `{Cl/2v3}` (Cl valence 3), `{Cl/1v5}` (Cl valence 5), `{Clv7}` (Cl valence 7), `{I/1v5}` (I valence 5), `{N/1^1v3}` (radical N in NO2).
- Added `heterocycles/` conformance category (53 files) with `implicit_h=true`. Covers O-saturated (7), N-saturated (5), S-saturated (4), O/S-unsaturated (6), N-unsaturated (5), functionalized (4: sulfolane, succinimide, 2-oxazolidone, hydantoin), multi-O/S (6), multi-N (5), mixed heteroatom (5), bridged/fused (3: pyrrolizidine, quinuclidine, 1-azaadamantane), crown ethers (3: 12-crown-4, 15-crown-5, 18-crown-6). All 53 counts-strategy results pass. 5 atom-typing results show `ValenceAmbiguous` on carbonyl O (same registry disambiguation issue as functional_groups).
- Added new `carbonyl_derivatives/` conformance category (62 files, non-aromatic) covering aldehydes, ketones, acids, esters, amides, acyl halides, anhydrides, carboxylates, lactones/thioesters, and C=N derivatives tied to carbonyl chemistry.
- Added 6 `functional_groups/` compounds requested for S/N chemistry expansion: dimethyl-sulfoxide, methylhydrazine, 1,1-dimethylhydrazine, 1,2-dimethylhydrazine, diazomethane, dimethyldiazene.
- Moved `beta-lactam` and `succinimide` from `heterocycles/` into `carbonyl_derivatives/` to keep carbonyl-centered compounds in one category.

## HW extraction spike (2026-10-31)

For each description, distinguish evidence for the species from evidence for the
bonding used to construct its encoding. Where HW is insufficient, follow its
references to primary papers (through reviews where necessary), or search independently.
Record the supporting passage and flag unresolved assignments. Electron conservation
and conformance checks do not establish source support for a bonding description.

Registry comments use HW volume/page references as the baseline. A source table at
the top of the registry defines short keys; primary references supplement HW only
where they add information used in an assignment. PW means “Present work —
calculations and adopted assignments,” marked as `assignment: PW` on the relevant
examples. Methods, results and assumptions are recorded here.

Default multicenter entries are limited to common patterns or selected molecules
for which no simpler description, multicenter or otherwise, is available. Delocalized
molecular orbitals or a possible electron allocation alone do not justify additional
atom types. Hydrogen-bond, donor–acceptor and aromatic descriptions retain their own categories.

**Multicenter counting limitation (2026-11-02):** the graph retains individual
relations and electron contributions, but registry matching uses the aggregate
per-atom contribution. Thus #m6 does not require three 3c4e units with the atom
in the central role. Describing the number and kind of units and each atom's role
could constrain counting more usefully. This is a future model-design question;
the Pn24 multicenter entries are accepted for now under the current model.

**Multicenter scope review (2026-11-02):** reviewed all 60 multicenter states in
hydrogen, halogens and noble gases, together with their B/Cr/N companions and
current fixture usage.

| Group | Disposition | Reason |
| --- | --- | --- |
| H: 3 states | All retained | H#m covers borane bridges; H#c+#m0 covers shared protons/H3+; H#c-#m2 covers bridging hydrides. These are recurring small bonding patterns; the B and Cr companion states support the selected bridges |
| Halogens: 37 states | 34 retained; 3 F states removed | F#c-#n2#m4, F#c-#n#m6 and F#c-#m8 were used only by H2F3−, H3F4− and H4F5−. HW 498, Fig. 12.8 and Iv17 support simpler F−·nHF hydrogen-bond descriptions. The remaining entries cover trihalides, small fluoride/oxyfluoride patterns and the common Br/I bis(pyridine) reagent family; its N companion remains |
| Noble gases: 20 states | All retained | Compact Kr/Xe fluoride and oxyfluoride families use the selected multicenter conventions; donor adducts retain those cores. No additional states arise from increasing an aggregate's size |

The three HF aggregates now use existing F#c-#n4, F#n3#v and H#v in
hydrogen_bonds/ fixtures. HF2− retains its common shared-proton description and
F#c-#n3#m2. XeO3F2 has explicit 3c4e discussion in W75, p. 1477; XeO2F4 remains
an inferred application of the selected even-fluoride convention, rather than a
separate aggregate pattern.

Neutral H2O, HF and NH3 dimers have hydrogen_bonds/{h2o,hf,nh3}-dimer.edn
fixtures. Each has one explicit Hbd contact and ordinary intramolecular bonds;
their lone-pair counts are left for resolution. No further atom types are needed.

### Hydrogen

Source: [HW 103, Band 1, chapter VIII](<../materials/valence/Holleman-Wiberg - Anorganische Chemie - 103. Auflage - Band 1.pdf#page=320>),
printed pp. 282–324, PDF pp. 320–362. The chapter's bonding cross-references were
also followed: HW I, pp. 174–176, 179, 1236, and
[HW II, pp. 2068–2073](<../materials/valence/Holleman-Wiberg - Anorganische Chemie - 103. Auflage - Band 2.pdf#page=572>).
Page references below are printed pages. Atom descriptions use AtomTypeRegistry::raise_defaults.
H01–H09 are implemented: the four earlier entries and five selected additions.
The other extracted examples are accounted for below.

#### Entries and selection status

| ID | Atom description | Example / source evidence | HW | Registry status |
| --- | --- | --- | --- | --- |
| H01 | H#v | H2 and ordinary covalent E–H bonds; HW explicitly describes two-center/two-electron bonding | I 289, 311 | Present |
| H02 | H#u | Atomic hydrogen | I 295–298 | Present |
| H03 | H#c+ | Proton | I 290 | Present |
| H04 | H#c-#n | Ionic hydrides: LiH, NaH, CaH2, etc. | I 290, 306–308 | Present; conformance gap filled |
| H05 | H#h | H2 represented by one explicit H carrying one implicit H; distinct from H01 under the loading defaults | I 287–290 | Selected; implemented |
| H06 | H#m | Neutral-H construction of the B–H–B bridges in B2H6; HW explicitly assigns three-center/two-electron bonds | I 179, 311, 1236 | Selected; implemented |
| H07 | H#c-#d | Terminal hydride donor, e.g. [RhH(NH3)5]2+; HW explicitly describes coordination through the H− lone pair | II 2068, 2071 | Selected; implemented |
| H08 | H#c+#m0 | Shared-proton construction of HF2− with two F#c-#n3#m2 partners, [2,0,2]; also H3+ with two H06 partners, [1,0,1] | I 174–176, 290–291, 311, 320 | Selected; implemented; charge placement is a constructed convention |
| H09 | H#c-#m2 | [HCr2(CO)10]−: closed three-center/two-electron Cr–H–Cr bridge (D79, p. 19, Fig. 2), using the selected H− convention | II 2069–2071 | Selected; implemented |

The atom strings are constructed from the source descriptions. H06 uses neutral B/H
bookkeeping: the two B–H–B units can carry contributions [1,1,0] and [0,1,1],
with the middle atom H. Each B then has two ordinary B–H bonds and total #m1.
H09 uses contributions [0,2,0] for Cr–H–Cr. D79 includes the partial Cr–Cr
bonding in this multicenter unit. The companion entries B#v2#m, Rh#c+3#n3#t6,
F#c-#n3#m2, Cr#n3#t5#m0 and C#c-#v3#d are implemented for the fixtures below;
the CO donors retain the C−/O+ triple-bond convention used for free CO.

**D79:** Darensbourg and Incorvia, *Inorg. Chem.* **18** (1979), 18–22;
[DOI](https://doi.org/10.1021/ic50191a004),
[PDF](../materials/valence/Darensbourg-1979-bridging-carbonyl-hydrides.pdf#page=2).

#### Registry and ValenceTable comparison

The AtomTyping registry covers H01–H09. After the Counts accounting correction
(2026-11-02), both models agree in all seven fixtures below, including the donated
pair in H07 and the multicenter contributions in H06/H09. The existing H target [1]
already covers H05; ValenceTable is unchanged.

#### Other extracted examples

| Example | H-state coverage / disposition | HW |
| --- | --- | --- |
| H(OH2)2+; H(NH3)2+; associated HF/H2O/NH3 | H08 covers the shared proton under the selected charge convention. Full encodings require companion-atom curation, rather than another H entry | I 174–176, 311 |
| W(H2)(CO)3(PCy3)2 and other η2-H2 complexes | H06 supplies the neutral-H contribution to HW's three-center/two-electron σ unit. The remaining full-complex work concerns metal/ligand assignments | I 291, 303; II 2071–2073 |
| H4Co4(C5H5)4; [H2Ru8(CO)21]2−; Mg2NiH4; HNb6I11; HRu6(CO)18 | H coordination numbers 3–6 do not themselves require further atom states. H09 covers a hydride contribution to a multicenter unit; individual metal/cluster assignments remain outside this hydrogen slice | I 303; II 2069–2070 |
| HeH+, CH5+; H2+, H4+, Hn+ (n = 5–11, 13, 15) | Retained as extraction results. No further default H entry is proposed from these occurrence examples alone | I 290–291 |

Metallic H and interstitial metal hydrides are described as extended phases; no discrete
atom state was extracted from their stoichiometries (HW I 289, 308–309).
H/D/T and ortho/para H2 introduce no additional electronic valence entries (I 298–302).

#### Conformance coverage

Existing cases cover atomic H, H+ and explicit H2. Seven fixtures now cover the additions
and the earlier H− gap; their lone-pair/spin fields are left for resolution.

| Added fixture | Coverage |
| --- | --- |
| [Implicit H2](../umol-graph/tests/resolution/data/hydrides/h2-implicit.edn) | H05 |
| [H−](../umol-graph/tests/resolution/data/atomic_ions/h-.edn) | H04 |
| [B2H6](../umol-graph/tests/resolution/data/inorganic_small/b2h6.edn) | H06 and B#v2#m; explicit terminal H, two three-center bridges |
| [[RhH(NH3)5]2+](../umol-graph/tests/resolution/data/donor_acceptor/pentaamminehydridorhodium-2+.edn) | H07 and Rh#c+3#n3#t6; explicit N–H bonds reuse N#v3#d |
| [HF2−](../umol-graph/tests/resolution/data/inorganic_ions/hf2-.edn) | H08 and F#c-#n3#m2; three-center/four-electron [2,0,2] |
| [H3+](../umol-graph/tests/resolution/data/inorganic_ions/h3+.edn) | H06 and H08; three-center/two-electron [1,0,1] |
| [[HCr2(CO)10]−](../umol-graph/tests/resolution/data/donor_acceptor/hcr2-co10-.edn) | H09, Cr#n3#t5#m0 and C#c-#v3#d; three-center/two-electron [0,2,0] |

Both models resolve all seven under Strict and MostSaturated. In [HCr2(CO)10]−,
the bridging H− and every CO carbon now receive #n0, as required by the supplied
multicenter and dative contributions. The correction is recorded in
[240](240-validation-constraint-ops-2026-10-18.md#counts-electron-accounting).

### Noble Gases

Source: [Holleman/Wiberg 103, Band 1, chapter XI](<../materials/valence/Holleman-Wiberg - Anorganische Chemie - 103. Auflage - Band 1.pdf#page=497>),
printed pp. 459–474, PDF pp. 497–512. Page references below are printed pages.
Atom strings use AtomTypeRegistry::raise_defaults. This audit covers all **34 current
noble-gas entries**, plus the F/N companions. Earlier and newly added entries use the
same evidence criterion. These are evidence findings; registry selections remain with
the user.

Current tally: **28 supported**, **6 inferred** and **0 open**.
F#n3#m and N#v3#d are supported. F#c-#n2#d2 is also selected for the
Xe2F11+ bridge as a constructed donor–acceptor description.

#### Current entries: evidence review

**Supported** means the indicated source description supports the constructed atom
state. **Inferred** means selected under the construction convention below, without
a verified individual bonding assignment. **Open** identifies incomplete bonding
evidence. References are listed below.

For neutral/cationic Xe, the selected construction convention assigns #v2 to each
Xe=O, #m<2m> to even fluoride sets, and #v to a localized Xe–F unit as in XeF+,
combined as needed. For the retained monocations, one Xe–F contributes #v and the
remaining even F set contributes #m. Source-specific localized bonds are retained.
For the selected anions, all F contribute to #m, with #v2 for each Xe=O.

| Atom descriptions | Evidence status and example | Source |
| --- | --- | --- |
| He; Ne; Ar; Kr#n4; Xe#n4; Rn#n4 | Supported: free atoms, retaining the user's He/Ne/Ar zero-lone-pair convention | HW 461 |
| Kr#n3#m2; Xe#n3#m2 | Supported: KrF2/XeF2, three-center/four-electron description; Kr mapping by analogy | HW 465, 469; C64, pp. 1445–1448 |
| Xe#n2#m4; Xe#n#m6 | Supported by XeF4/XeF6 | HW 467–468; C64, pp. 1445–1448 |
| Kr#c+#n3#v; Xe#c+#n3#v | Supported: KrF+, XeF+, XeCl+ | HW 466, 469 |
| Kr#c+#n3#v#t; Xe#c+#n3#v#t | Supported: nitrile donor adducts; HW explicitly gives HCN→KrF+ | HW 466, 469, 473; B07, p. 3608 |
| Xe#c+#n2#v#m2 | Supported: M73 explicitly describes axial three-center bonding, an equatorial electron-pair bond and two Xe lone pairs in XeF3+. K14 independently supports C6F5XeF2+ | HW 467, 473; M73, pp. 1714, 1716; K14, pp. 11656–11657 |
| Xe#c+#n#v#m4 | Supported: G10 explicitly describes XeF5+ with two equatorial three-center/four-electron bonds, one covalent axial Xe–F bond and one Xe lone pair | HW 468; G10, p. 265 |
| Xe#c+#n#v#t#m4 | Supported as a constructed donor–acceptor description: H10 describes Xe2F11+ as two XeF5+ units sharing a fluoride ion; each Xe accepts one pair from that F−, retaining the G10 core description | HW 468; G10, p. 265; H10, pp. 3504, 3513 |
| Xe#c-#n2#m5 | Supported: XeF5−, all five equivalent F in the multicenter unit, with two Xe lone pairs | HW 467; C91, pp. 3351, 3360 |
| Xe#c-#n#m7 | Inferred, selected: XeF7− with one Xe lone pair and all seven F in #m7. P64/P66 establish the salts; E96/G10 establish structures. The Wiberg comparison below supports the selection | HW 468; P64; P66; E96, pp. 1123–1124; G10, pp. 265–266 |
| Xe#c-2#n#m8 | Inferred, selected: XeF8²− with one Xe lone pair and all eight F in the multicenter description, by analogy with XeF2/XeF+ and the D4d geometry. A96 establishes the species, geometry and inactive lone pair; the calculation check below supports the lone pair | HW 468; A96, pp. 400–401 |
| Xe#n#v6; Xe#n#v6#t | Supported: XeO3 and its 1:1 nitrile adduct; B72 explicitly allows Xe=O, and G16 uses Xe=O bonds, one Xe lone pair and N-lone-pair donation | HW 469–470; B72, p. 1127; G16, pp. 13780–13782 |
| Xe#v8 | Supported by XeO4: B72 explicitly permits four Xe=O bonds; G70 favors the double-bond description from bond lengths. Perxenate [XeO6]4− is also selected, using the two-double/four-single encoding as an accepted assumption | HW 470–471; B72, p. 1127; G70, pp. 814–815 |
| Xe#n2#v2#m2 | Supported as a formal description of XeOF2: B07 includes Xe=O valence structures and explicitly describes F–Xe–F as three-center/four-electron | HW 471; B07, pp. 3601, 3607 |
| Xe#n#v2#m4 | Supported as a formal description of XeOF4: HW supplies Xe=O; W75 explicitly uses Rundle three-center/four-electron fluoride bonding in this series | HW 471; W75, p. 1477; H11, p. 666 |
| Xe#n#v2#t#m4 | Supported as a constructed donor–acceptor description: each XeOF4 in [F(XeOF4)3]− accepts one pair from F−, retaining its fluoride/oxo core | HW 471, 499; W75, p. 1477 |
| Xe#n#v4#m2 | Supported as a formal description of XeO2F2: P06 supports the two Xe=O contributors; W75 supplies the three-center/four-electron F–Xe–F description | HW 471; W75, p. 1477; P06, pp. 1529, 1532; SI p. S2 |
| Xe#v6#m2 | Inferred: XeO3F2 with three Xe=O and an F–Xe–F multicenter unit. CH71 establishes geometry and W75 supports the fluoride assignment; the Xe=O assignments use the construction convention | HW 471; CH71, pp. 1506–1507; W75, p. 1477 |
| Xe#v4#m4 | Inferred: XeO2F4 with two Xe=O and four multicenter F contributions. H82 establishes detection and volatility; the bonding assignment uses the construction convention | HW 465, 471; H82, pp. 687–688 |
| Xe#c+#n#v5 | Supported as a formal description of XeO2F+: P06 gives two Xe=O bonds, one Xe–F bond and one Xe lone pair; its bonding analysis includes this Lewis structure | HW 471; P06, pp. 1519, 1531 |
| Xe#c+#n#v3#m2 | Supported: M93 explicitly assigns axial three-center/four-electron F–Xe–F and equatorial two-center/two-electron Xe–F, and discusses the Xe=O multiple bond in XeOF3+ | HW 471; M93, pp. 388, 391 |
| Xe#c-#n2#v2#m3 | Inferred, selected: XeOF3− with Xe=O, two Xe lone pairs and all three F in #m3. B10 distributes fluoride bonding over all three F; the Wiberg comparison below supports this construction | HW 471; B10, p. 10940; SI p. S9 |
| Xe#c-#n#v2#m5 | Supported as a formal description of XeOF5−: one Xe=O, one Xe lone pair and all five equatorial F in a six-center/ten-electron unit | HW 471; C95, pp. 1873–1874; E95, pp. 1586–1587; B00, p. 62 |
| Xe#c-#n#v4#m3 | Inferred, selected: XeO2F3− with two Xe=O, one Xe lone pair and all three F in #m3. C88 establishes the cis assignment; the Wiberg comparison below supports the selection | HW 471; C88, pp. 3767–3768 |
| F#n3#m; N#v3#d | Supported: terminal F in the multicenter fluorides; nitrile donor N | HW 465–466, 469; C64; B07; G16 |
| F#c-#n2#d2 | Inferred, selected: bridging F− donates one pair to each of two XeF5+ units in Xe2F11+ | HW 468; H10, pp. 3504, 3513; PW |

#### Primary sources read

PDFs are retained in [materials/valence](../materials/valence/); text extracts and
page images are temporary scratch intermediates.

- **C64:** Coulson, “The Nature of the Bonding in Xenon Fluorides,” *J. Chem. Soc.* (1964), 1442–1454. [DOI](https://doi.org/10.1039/JR9640001442).
- **P64:** Peacock, Selig and Sheft, “Complex Fluoroxenates(VI),” *Proc. Chem. Soc.* (1964), 285. [PDF](../materials/valence/Peacock-1964-fluoroxenates.pdf).
- **P66:** Peacock, Selig and Sheft, “Complex Compounds of Xenon Hexafluoride with the Alkali Fluorides,” *J. Inorg. Nucl. Chem.* **28** (1966), 2561–2567. [DOI](https://doi.org/10.1016/0022-1902(66)80380-9).
- **G70:** Gundersen, Hedberg and Huston, “Molecular Structure of Xenon Tetroxide, XeO4,” *J. Chem. Phys.* **52** (1970), 812–815. [DOI](https://doi.org/10.1063/1.1673060).
- **CH71:** Claassen and Huston, “Matrix-Isolation Raman and Infrared Spectra, Molecular Structure, and Force Constants for XeO3F2,” *J. Chem. Phys.* **55** (1971), 1505–1507. [DOI](https://doi.org/10.1063/1.1676271).
- **B72:** Bartlett et al., “The Crystal Structure of Xenon(II) Fluoride Fluorosulfate, FXeOSO2F,” *Inorg. Chem.* **11** (1972), 1124–1127. [DOI](https://doi.org/10.1021/ic50111a042).
- **M73:** McKee, Zalkin and Bartlett, “Crystal Structure of [XeF3+][Sb2F11−],” *Inorg. Chem.* **12** (1973), 1713–1717. [DOI](https://doi.org/10.1021/ic50126a001); [PDF](../materials/valence/McKee-1973-XeF3-cation.pdf).
- **W75:** Willett, LaBonville and Ferraro, “Normal Coordinate Treatment of XeO2F2,” *J. Chem. Phys.* **63** (1975), 1474–1478. [DOI](https://doi.org/10.1063/1.431510).
- **H82:** Huston, “Chemical and Physical Properties of Some Xenon Compounds,” *Inorg. Chem.* **21** (1982), 685–688. [DOI](https://doi.org/10.1021/ic00132a043).
- **S85:** Stein, “New Evidence that Radon is a Metalloid Element: Ion-exchange Reactions of Cationic Radon,” *J. Chem. Soc., Chem. Commun.* (1985), 1631–1632. [DOI](https://doi.org/10.1039/C39850001631).
- **C88:** Christe and Wilson, “New Syntheses and Properties of XeO2F2, Cs+XeO2F3−, and NO2+[XeO2F3·nXeO2F2]−,” *Inorg. Chem.* **27** (1988), 3763–3768. [DOI](https://doi.org/10.1021/ic00294a019).
- **C91:** Christe et al., “The Pentafluoroxenate(IV) Anion, XeF5−: The First Example of a Pentagonal Planar AX5 Species,” *JACS* **113** (1991), 3351–3361. [DOI](https://doi.org/10.1021/ja00009a021).
- **M93:** Mercier et al., “The Oxotrifluoroxenon(VI) Cation: X-ray Crystal Structure of XeOF3+SbF6− and a Solution 17O and 129Xe Nuclear Magnetic Resonance Study of the 17,18O-Enriched XeOF3+ Cation,” *Inorg. Chem.* **32** (1993), 386–393. [DOI](https://doi.org/10.1021/ic00056a009).
- **C95:** Christe et al., “On the Structure of the [XeOF5]− Anion and of Heptacoordinated Complex Fluorides Containing One or Two Highly Repulsive Ligands or Sterically Active Free Valence Electron Pairs,” *Inorg. Chem.* **34** (1995), 1868–1874. [DOI](https://doi.org/10.1021/ic00111a039).
- **E95:** Ellern and Seppelt, “XeOF5−, an Anion with Pentagonal-Pyramidal Structure,” *Angew. Chem. Int. Ed. Engl.* **34** (1995), 1586–1587. [DOI](https://doi.org/10.1002/anie.199515861). Reports five approximately coplanar F and apical O; assigns partial Xe–O double-bond character and a trans lone pair.
- **A96:** Adam, Ellern and Seppelt, “Structural Principles of the Coordination Number Eight: WF8²−, ReF8²−, and XeF8²−,” *Chem. Eur. J.* **2** (1996), 398–402. [DOI](https://doi.org/10.1002/chem.19960020408); [PDF](../materials/valence/Adam-1996-XeF8-structure.pdf).
- **E96:** Ellern, Mahjoub and Seppelt, “Structures of XeF7− and Xe2F13−,” *Angew. Chem. Int. Ed. Engl.* **35** (1996), 1123–1125. [DOI](https://doi.org/10.1002/anie.199611231).
- **B00:** Baran, “Mean Amplitudes of Vibration of the Pentagonal Pyramidal XeOF5− and IOF5²− Anions,” *J. Fluorine Chem.* **101** (2000), 61–63. [DOI](https://doi.org/10.1016/S0022-1139(99)00194-3).
- **P06:** Pointner, Suontamo and Schrobilgen, “Syntheses and X-ray Crystal Structures… and Computational Studies of the XeO2F+ and FO2XeFXeO2F+ Cations and Related Species,” *Inorg. Chem.* **45** (2006), 1517–1534. [DOI](https://doi.org/10.1021/ic0510260); [SI](https://doi.org/10.1021/ic0510260.s002), p. S2 and Table S8.
- **B07:** Brock et al., “XeOF2, F2OXeN≡CCH3, and XeOF2·nHF: Rare Examples of Xe(IV) Oxide Fluorides,” *JACS* **129** (2007), 3598–3611. [DOI](https://doi.org/10.1021/ja0673480).
- **B10:** Brock et al., “XeOF3−, an Example of an AX3YE2 Valence Shell Electron Pair Repulsion Arrangement…,” *JACS* **132** (2010), 10935–10943. [DOI](https://doi.org/10.1021/ja103730c); [SI](https://doi.org/10.1021/ja103730c.s001), p. S9, Scheme S1.
- **G10:** Grant et al., “Heats of Formation of XeF3+, XeF3−, XeF5+, XeF7+, XeF7−, and XeF8 from High Level Electronic Structure Calculations,” *Inorg. Chem.* **49** (2010), 261–270. [DOI](https://doi.org/10.1021/ic901956g); [SI](../materials/valence/Grant-2010-xenon-fluoride-ions-SI.pdf), Table SM-3.
- **H10:** Hughes, Mercier and Schrobilgen, “Syntheses, Raman Spectra, and X-ray Crystal Structures of [XeF5][μ-F(OsO3F2)2] and [M][OsO3F3] (M = XeF5+, Xe2F11+),” *Inorg. Chem.* **49** (2010), 3501–3515. [DOI](https://doi.org/10.1021/ic100062x).
- **H11:** Hughes et al., “A Raman Spectroscopic Study of the XeOF4/XeF2 System and the X-ray Crystal Structure of α-XeOF4·XeF2,” *J. Fluorine Chem.* **132** (2011), 660–668. [DOI](https://doi.org/10.1016/j.jfluchem.2011.05.010).
- **P13:** Papakondylis, “Xenon oxides, sulfides, and oxysulfides,” *Comput. Theor. Chem.* **1015** (2013), 15–20. [DOI](https://doi.org/10.1016/j.comptc.2013.03.030).
- **K14:** Koppe et al., “Xenon(IV)–Carbon Bond of [C6F5XeF2]+; Structural Characterization and Bonding…,” *Inorg. Chem.* **53** (2014), 11640–11661. [DOI](https://doi.org/10.1021/ic501831j).
- **G16:** Goettel et al., “Syntheses and Structures of Xenon Trioxide Alkylnitrile Adducts,” *Angew. Chem. Int. Ed.* **55** (2016), 13780–13783. [DOI](https://doi.org/10.1002/anie.201607583).

Full text not retrieved (2026-10-31):

- **P88:** Popov, A. I.; Fadeeva, N. E.; Kiselev, Yu. M., “Oxidimetric determination of active fluorine in fluorides used as fluorinating and oxidizing agents,” *Zhurnal Analiticheskoi Khimii* **43**(3) (1988), 465–471. User-supplied reference from the XeF8²− compound searches; no online full text located, bonding content unassessed.

#### XeF8²− calculation check (2026-10-31)

PySCF 2.14.0, RHF/def2-TZVPP with the 28-electron Xe def2 ECP, at the
user-supplied MP2 geometry: D4d, Xe–F = 2.04889054 Å, upper-square polar angle
57.32284694°, lower square reflected in z and rotated 45°. XeF2 (linear) and
XeF+ references use that same distance without optimization. Wiberg indices use
the spin-summed density in PySCF's orthonormal NAO basis.

| Species | Xe–F Wiberg index per contact | Sum at Xe |
| --- | --- | --- |
| XeF8²− | 0.38480 | 3.07842 |
| XeF2 | 0.48208 | 0.96417 |
| XeF+ | 0.93872 | 0.93872 |

XeF8²−: 298 basis functions, 100 explicit electrons; converged RHF energy
−1123.50405686658 Eh. HOMO/LUMO energies agree with the supplied TURBOMOLE
output within 1 × 10⁻⁶ Eh. PySCF's Xe NAO charge is +3.70911, versus +3.66218
in TURBOMOLE; the population analyses are not numerically identical.
TURBOMOLE gives the Xe 5s NAO an occupation of 1.99841 electrons. The user selected
Xe#c-2#n#m8 by analogy with XeF2/XeF+ and the D4d geometry; the indices are supporting
comparisons, not a rule assigning integer #v/#m values.

#### Remaining-anion calculation checks (2026-10-31)

The same RHF/def2-TZVPP/ECP28 and PySCF NAO Wiberg calculation gives:

| Anion | Geometry used | Xe–F Wiberg indices | Xe–O Wiberg indices |
| --- | --- | --- | --- |
| XeF7− | G10, SI Table SM-3, C3v CCSD(T)/aVTZ | cap 0.36048; near triangle 0.43215 × 3; far triangle 0.47983 × 3 | — |
| XeOF3− | B10, Table 3, C2v MP2/aug-cc-pVTZ(-PP) | Fa 0.39454 × 2; Fb trans to O 0.20643 | 1.11069 |
| XeO2F3− | New converged MP2/def2-TZVPP geometry | axial 0.39019, 0.39016; equatorial 0.24292 | O trans to F 1.06450; other O 1.13495 |

E96 establishes a C3v crystal structure of CsXeF7. G10's gas-phase C2v and C3v
structures differ by only 0.02 kcal/mol; the calculation here samples C3v only.
B10 establishes XeOF3− by Raman spectroscopy and calculation, without a crystal
structure. C88 assigns cis XeO2F3− from vibrational spectra; no atomic coordinates
were found in the reviewed sources.

Published geometry parameters reconstructed for these checks (Å, degrees):
XeF7− cap r = 2.0441; near triangle r = 2.0046, polar angle 75.4; far triangle
r = 1.9667, polar angle 131.3, azimuth offset 60. XeOF3− r(Xe–O) = 1.821,
r(Xe–Fa) = 2.060, r(Xe–Fb) = 2.129, O–Xe–Fa = 94.8, O–Xe–Fb = 180.
The latter reconstruction gives Fa–Xe–Fa = 170.4°, versus the independently
rounded 170.3° in B10. All three RHF calculations converged; density normalization,
NAO orthogonality and closed-shell idempotency checks passed.
B10 SI Table S7 labels its quantities NBO bond orders, without identifying them
as Wiberg indices; they are not used as the numerical reference here.

XeO2F3− uses an unconstrained MP2/def2-TZVPP optimization with Xe ECP28 and all
explicit electrons correlated, starting from C88's cis topology with estimated
distances. It converged in 11 evaluations: E = −778.15641465 Eh, maximum gradient
component 0.00004959 Eh/bohr; no frequency calculation was performed. Coordinates
used (Å, rigidly aligned; O trans to F first, equatorial F last):

```text
Xe   0.00000000   0.00000000   0.00000000
O    1.76385519   0.00000000   0.00000000
O   -0.28197632   1.70408983   0.00000000
F   -0.11278732  -0.11252632   1.99680396
F   -0.11285876  -0.11151882  -1.99689469
F   -2.07556951  -0.26762219   0.00005425
```

Compared with the preliminary calculation at evaluation 7, the largest Wiberg-index
change is 0.00276. The selected registry replacements put all F contributions in
#m and retain #v2 per Xe=O.

#### Selections applied

Perxenate [XeO6]4− retains Xe#v8 as an accepted assumption. Its existing
[conformance input](../umol-graph/tests/resolution/data/inorganic_ions/xeo64-.edn)
uses two Xe=O and four Xe–O− bonds, with neutral Xe. This settles the final
noble-gas example; the atom state already has independent support from XeO4.

F#c-#n2#d2 is added for the Xe2F11+ bridge, donating one pair to each
Xe#c+#n#v#t#m4 center. Its lone-pair count is left for resolution in the
[conformance input](../umol-graph/tests/resolution/data/donor_acceptor/xe2f11+.edn).

| Implemented replacement | Replaced entry | Example |
| --- | --- | --- |
| Xe#c-#n2#m5 | Xe#c-#n2#v#m4 | XeF5− |
| Xe#c-#n#v2#m5 | Xe#c-#n#v3#m4 | XeOF5− |
| Xe#c-2#n#m8 | Xe#c-2#n#v2#m6 | XeF8²− |
| Xe#c-#n#m7 | Xe#c-#n#v#m6 | XeF7− |
| Xe#c-#n2#v2#m3 | Xe#c-#n2#v3#m2 | XeOF3− |
| Xe#c-#n#v4#m3 | Xe#c-#n#v5#m2 | XeO2F3− |

These mappings retain the formal charge on Xe and use the existing F#n3#m
description. Each n-fluoride multicenter unit contains 2n electrons. The states
conserve electrons under the registry defaults; multicenter descriptions remain
outside the Counts comparison.

Rn#c+2#n3 is excluded by the user: S85, pp. 1631–1632, proposes Rn²⁺ but does
not distinguish it from RnF+. The default entry and ionic RnF2 conformance fixture
are removed. Rn#n3#m2, Rn#n2#m4 and Rn#n#m6 are also removed for insufficient
support for their compound descriptions (HW 465, 469; S85); only Rn#n4 remains.

Xe#m8 is removed: its oxide/perxenate examples lack support for an all-multicenter
description (HW 470–471; P13). Neutral XeF8 is unestablished in HW and was not an
example for this entry. Obsolete oxide examples are removed from retained fluoride
entries; B07 supersedes HW's polymeric XeOF2 interpretation.

Xe#t#m8 and Xe#v8#t and their MeCN→XeO4 fixtures are removed by the user's
decision: no verified XeO4 adduct was found; G16, p. 13780, eq. 2, instead gives
MeCN→XeO3 + ½O2. N#v3#d is retained for supported N-donor complexes.

MeCN→XeO3 retains only Xe#n#v6#t, supported by G16. The conflicting
Xe#n#t#m6 entry is removed.

#### Further extracted examples and qualifications

| Examples | Source detail | HW pages |
| --- | --- | --- |
| Xe2F3+, Kr2F3+, Xe2F11+, Xe2F13−; oligomeric XeF6 | Bridging/association; Kr2F3+ is asymmetric | 466, 468–469 |
| XeCl2, XeCl4, XeBr2, XeClF | XeCl2 in a cold condensate; XeCl2/4 and XeBr2 as decay products; XeClF mentioned | 465, 469 |
| HeH+, EgH+, HArF, HKrCl, HKrCN, HXeCl/Br/I/CN/NC/C2H, XeH2 | Gas-phase and matrix chemistry; HW describes HArF as H–Ar+ / F− | 472 |
| Ne₂⁺, Ar₂⁺, Kr₂⁺, Xe₂⁺ | Mass spectrometry; an Xe₂⁺ salt is also reported | 472 |
| Kr(OTeF5)2; Xe(OTeF5)2/4/6; XeF(OAc), Xe(OAc)2; XeF[N(SO2F)2], Xe[N(SO2F)2]2 | Further O/N ligand examples | 472–473 |
| C6F5Xe+ and MeCN adduct; C6F5XeF/Cl/CN, C6F5XeOCO(C6F5), Xe(C6F5)2, C6F5XeF2+ | Xe–C examples | 473 |

Not selected: Xe#n3#d for AuXe4²+/Au–Xe/Hg–Xe complexes (HW 473);
Xe#c+#n3#m for Xe₂⁺ (trial encoding, HW 472); Xe#c+2#n3#t2 for
Xe(FH)2²+ (proposed interpretation in HW 466).

He#n, Ne#n4 and Ar#n4 are excluded by the user: retain the zero-lone-pair
descriptions and current valence limits. The trial descriptions He#c+#v,
Ne#c+#n3#v, Ar#c+#n3#v, Ne#n3#m2, Ar#n3#m2 and Ne#c+#n3#m are not selected.

HW marks ArF2, KrF4, neutral XeF8 and XeF7+ as unestablished; ArF+ and XeOF6
are questioned. XeO, XeO2 and Xe(OH)2 occur as possible intermediates.

#### Implementation checks

F#n3#m, F#c-#n2#d2, N#v3#d and the selected noble-gas descriptions are in the default registry.
All noble-gas entries carry HW volume/page references. After the selections above,
the noble-gas increment retains 26 conformance inputs, including
[XeO3](../umol-graph/tests/resolution/data/inorganic_small/xeo3.edn),
[XeOF4](../umol-graph/tests/resolution/data/inorganic_small/xeof4.edn),
[XeF5−](../umol-graph/tests/resolution/data/inorganic_ions/xef5-.edn),
[XeOF5−](../umol-graph/tests/resolution/data/inorganic_ions/xeof5-.edn),
[XeF8²−](../umol-graph/tests/resolution/data/inorganic_ions/xef8-2.edn) and
[XeF7−](../umol-graph/tests/resolution/data/inorganic_ions/xef7-.edn).
The existing [XeOF3−](../umol-graph/tests/resolution/data/inorganic_ions/xeof3-.edn)
and [XeO2F3−](../umol-graph/tests/resolution/data/inorganic_ions/xeo2f3-.edn) inputs
also use the selected replacements. Nine further fixtures cover KrF2, KrF+,
HCN→KrF+, XeF4, XeF6, XeF+, XeF3+, XeF5+ and Xe2F11+. All **33 retained
noble-gas states** and the three F/N companions are now exercised. Both AtomTyping
policies match the specified states, including two lone pairs on the bridging F−
and one on each Xe in Xe2F11+.
All **711 conformance cases and 366 valence unit tests pass**,
preserving the 685 cases preceding the HW increment.
The perxenate encoding resolves under both policies with AtomTyping and Counts.
No noble-gas selections remain open.

After the Counts accounting correction (2026-11-02), its multicenter and dative
outputs also match the selected AtomTyping assignments under both policies.
ValenceTable is unchanged.

### Halogens

Source: [HW 103, Band 1, chapter XII](<../materials/valence/Holleman-Wiberg - Anorganische Chemie - 103. Auflage - Band 1.pdf#page=513>),
printed pp. 475–552, PDF pp. 513–590. HW references below are printed Band I pages.
Atom descriptions use AtomTypeRegistry::raise_defaults. Charges and contributions
are constructed assignments; X expands only over the elements named in its row.

The registry currently has **110 halogen descriptions**: F 12, Cl 27, Br 28, I 32, At 11.
Existing monovalent, atomic and halide entries supply the baseline (HW 477–506).
The existing Cl/Br/I states X#n2#v3, X#n#v5 and X#v7 have independent oxoacid/oxide
examples (HW 512–546); revising their fluoride examples need not remove these states.
The existing X#c+#n2#v2 states cover XF2+; S20 directly supports localized ClF2+.
The bare Cl+/Br+/I+ entries are retained for electrophilic halogenation mechanism
representations; no bare F+ entry is selected.

#### Whole-registry audit (2026-10-31)

All 112 entries were reviewed, including the 37 pre-existing entries. Binary molecules
and simple ions remain in scope, including short-lived species. The pruning criterion
is an extra atom state needed only for an obscure complex or complicated ion; a large
example alone does not justify removing a state with simpler uses. Small oxyfluorides
and their mononuclear ions remain included.

AtomTyping scans the element/charge bucket; an unspecified charge uses all entries
for that element. Bucket sizes after the 2026-11-02 multicenter scope review are:

| Element | Charge 0 | Charge +1 | Charge −1 | Total |
| --- | ---: | ---: | ---: | ---: |
| F | 4 | 3 | 5 | 12 |
| Cl | 13 | 5 | 9 | 27 |
| Br | 14 | 5 | 9 | 28 |
| I | 15 | 7 | 10 | 32 |
| At | 6 | 1 | 4 | 11 |

**Removed from the default for now:**

| State | Reason for omission | Bucket reduction |
| --- | --- | --- |
| Cl#c+#n3#m0 | Specialized bis(pyridine)chloronium coordination. R20 establishes the 3c4e interpretation; P23, pp. 2326–2327 supplies solid-state structures obtained at low temperature. Supported, but less broadly useful than the Br/I reagent families | Cl+: 6 → 5 |

All other states are retained under this criterion. In particular, Cl#c-#n#d3 covers
the binary SeCl4/TeCl4 aggregates (HW 711), I#c+#n2#v2#t2 covers I2Cl6, and the
terminal/doubly bridging halide states are general ligand roles. Br/I acceptor states
also cover pentahalides; I2 donation and Br/I bis(pyridine) reagent families remain
selected (HW 489, 493–495; P23, p. 2325).

**Bridge-scope revision, 2026-11-01:** F#c-#n#d3 is restored alongside the new
Br−/I− triply donating states. The earlier exclusion judged the particular fluoride
complexes too elaborate, but that does not consistently exclude a general ligand
role while retaining the corresponding binary halide aggregates. HW 499 explicitly
gives [F(XeOF4)3]− and the binary [F(SbF3)3]−; the earlier SbF5 transcription was wrong.
The subsequent multicenter scope review removes three HF-aggregate states, leaving
110 states; neutral buckets are unchanged. These counts are not timing results.

Citation corrections cover the previously unreferenced baseline entries, the shared
Li96 key for the old trihalogen-cation DOI, and C73/S20 for the explicit localized
bonding evidence. At#n3#v now cites AtBr/AtI (HW 506), IOF3/IOF4− cite HW 545,
and ClF6− includes HW 498. The bare Cl+/Br+/I+ states are labelled PW mechanism
fragments. I#c-#n#v6 and I#c-#v8 also cover IO3− and IO4− with charge on I under
the adopted formal-charge convention; their usefulness is not confined to the rarer
examples previously listed. The two corresponding fixtures and snapshots are
removed; the ValenceTable is unchanged.

#### Consolidated proposals

**Accepted** means selected by the user or explicitly described in reviewed literature,
awaiting implementation; **Proposed** means one surviving assignment awaiting selection;
**Open** identifies an unresolved assignment; **Present** means implemented;
**Not selected** records an explicit omission from the default registry.
Explicit literature assignments are accepted automatically within the selected scope. Ask the user only where
there is ambiguity, controversy or no literature support. List inaccessible papers
before considering an inferred assignment. Earlier row IDs are retained for reference;
Open rows carry no asserted atom description where the bonding is unsettled.

| Rows | Atom descriptions | Status | Examples; evidence or open point |
| --- | --- | --- | --- |
| L01 | F#c+#n2#v2; F#c+#n2#v#h; F#c+#n2#h2 | Present | H2F+ with explicit, mixed or implicit H; HW 498–499 |
| L02 | At#c-#n4; At#c+#n3 | Present | Astatides; At+ precipitation as At2S and cathodic deposition; tracer chemistry, HW 490, 492, 499, 502 |
| L03 | At#n2#v3; At#n#v5; At#v7; At#c-#n#v6; At#c-#v8 | Present | HAtO2/HAtO3/HAtO4; separate At− entries for AtO3−/AtO4− with charge on At and At=O bonds; HW 512–513 |
| L11 | At#n3#h | Present | HAt in the HX redox series; HW 502 supplies no preparation or structure |
| L04, oxo examples | Cl#c+#n#v4; I#c+#n#v4 | Present | ClO2+, IO2+ units; HW 540, 543, 545 |
| L04, valence six | X#c+#v6 (X = Cl, Br, I) | Present | XF6+: Le04, pp. 6917–6919 describes two-center, highly polar single X–F bonds. ClO2F2+: C73, p. 1362 supports two Cl=O and two localized Cl–F bonds. I(OH)6+ independently supports I#c+#v6 (HW 528) |
| L06, oxo example | I#c-#n#v6 | Present | H2IO4−, with two I=O and two I–OH bonds; HW 527 |
| L08 | X#n2#u#v2 (X = Cl, Br, I) | Present | XO radicals; formal X=O assignment with unpaired electron on X; HW 530–531 discusses delocalization |
| L09 | X#n#u#v4; X#u#v6 (X = Br, I) | Present | XO2/XO3 radicals; HW 531. Corresponding Cl states are already present |
| L10 | N#c-2#n3#v | Present | User-added N entry covers [ClO3N]2− with three Cl=O and one Cl–N; HW 546 |
| L05/M01/M02 | X#n3#m; X#c-#n3#m2 (X = Cl, Br, I) | Present | X3−, ClF2−/BrF2−/IF2−, BrCl2−/ICl2−/IBr2−; HW 393, 493–494, 510–511 and L97 support the trihalide construction. Contributions [1,2,1] and central charge are our convention. Both F states already exist |
| M02, At | At#c-#n3#m2 | Present; by analogy | AtBr2−/AtI2− occur in HW 510; the analogous trihalide assignment is selected. No direct At bonding analysis reviewed |
| L04, tetrafluorides | X#c+#n#v2#m2 (X = Cl, Br, I) | Present | XF4+: two localized equatorial X–F bonds and an axial 3c4e F–X–F unit. Constructed assignment from C01, p. 6347, Table 22 and discussion of axial 3c4e contributions; replaces the all-localized fluoride proposal |
| L06/M03, tetrahalides | X#c-#n2#m4 (X = Cl, Br, I) | Present | XF4−: two 3c4e F–X–F units, each with contributions [1,2,1]; C68, pp. 125, 132–133, and CN73, pp. 61–62 support the assignment. Extension to square-planar BrCl4−/ICl4−/IBr4− is by analogy (HW 510–511). Terminal F#n3#m, Cl#n3#m and Br#n3#m are present |
| L06/M04, hexafluorides | X#c-#n#m6 (X = Cl, Br, I) | Present; by analogy | XF6−: three 3c4e F–X–F units, each [1,2,1], extending the selected XeF6 assignment. CW89, pp. 3276–3277; D08, pp. 5488, 5493; M22, pp. 10–11 support one central lone pair. D08 invokes 3c4e contributions in iodine fluorides but does not specify this all-multicenter partition. HW 498 includes ClF6−; HW 510–511 covers BrF6−/IF6− |
| L07/M05 | I#c-#v8 | Present | IF8−: eight localized I–F bonds, no lone pair; selected no-lone-pair convention. Square-antiprismatic D4d structure: HW 510–511; D08, p. 5487, Fig. 1 and pp. 5489–5490 |
| M06, trifluorides | X#n2#v#m2 (X = Cl, Br, I) | Present; Br/I by analogy | XF3: one localized equatorial X–F bond, one axial 3c4e F–X–F unit [1,2,1], two lone pairs. Sh73, pp. 5873–5875 explicitly supports this bond partition for ClF3. HW 508–509 gives the common T-shaped structures; D08, pp. 5487, 5493 supports the axial/equatorial distinction in IF3 and discusses 3c4e contributions |
| M07, pentafluorides | X#n#v#m4 (X = Cl, Br, I) | Present; Br by analogy | XF5: one lone pair, one localized apical X–F bond, two basal 3c4e F–X–F units, each [1,2,1]. C68, pp. 132–133 supports this partition for ClF5; C93A, pp. 2705–2706 explicitly invokes basal multicenter bonding and a more covalent apical bond in IF5. M22, p. 10 supports the apical/basal distinction in BrF5; HW 508–509 |
| M07, IF7 | I#v7 | Present; retained | IF7: seven localized I–F bonds; no central lone pair to reserve s character. Selected convention applies the lone-pair criterion in Le04, pp. 6917–6918. HW 508–509; D08, p. 5487, Fig. 1 |
| M08, Cl/Br | X#n#v3#m2 (X = Cl, Br) | Present | ClOF3/BrOF3: X=O, one localized equatorial X–F bond, one axial 3c4e F–X–F unit [1,2,1], one lone pair. Explicit bonding descriptions in OC82, p. 275 and C78, pp. 1537–1538; E02 confirms the structures; HW 544–545 |
| M08, I | I#n#v3#m2 | Present; by analogy | IOF3: same partition as ClOF3/BrOF3. ET74, pp. 175, 178–179 and CA71, p. 2299 establish the molecular structure and axial/equatorial distinction, but do not explicitly assign 3c4e bonding; E02, p. 1992 concurs on the structure |
| M09, Cl | Cl#v7 | Present | ClO2F3: two Cl=O and three localized Cl–F bonds, no lone pair. Selected no-lone-pair convention supersedes C73D's axial 3c4e partition; HW 545–546 |
| M09, I | I#v7 | Present | IO2F3: localized I–F bonds, no lone pair. S77 establishes the oxygen-bridged crystal dimer: each I has one terminal I=O, three I–F and two ordinary I–O bridge bonds. The monomer above 100 °C has two I=O bonds and the same atom state (HW 546). Selected no-lone-pair convention |
| M10 | I#v7 | Present | IOF5: I=O and five localized I–F bonds, no lone pair; selected no-lone-pair convention. C93O, pp. 9656–9658 gives comparable axial/equatorial I–F lengths and force constants; HW 546 |
| L06/M12, IOF4− | I#c-#n#v2#m4 | Present | IOF4−: I=O, one lone pair and two basal 3c4e units [1,2,1]. C93A, pp. 2705–2706 explicitly invokes multicenter bonding for the four basal I–F bonds; E02, pp. 1995–1996 compares the structures |
| L07/M13, IOF6− | I#c-#v8 | Present | IOF6−: I=O and six localized I–F bonds, no lone pair; selected no-lone-pair convention. C93A establishes the structure; its equatorial 6c10e partition (pp. 2704–2706) is not adopted |
| L06/M11 | X#c-#n#v4#m2 (X = Cl, Br) | Present; Br by analogy | XO2F2−: two X=O bonds, one lone pair and an axial F–X–F 3c4e unit [1,2,1]. C72D, pp. 38–39 explicitly supports this partition for Cl. GS76, p. 2393 supports the corresponding Br structure; the same partition is adopted. HW 545 |
| L06/M12, Cl/Br | X#c-#n#v2#m4 (X = Cl, Br) | Present; by convention | XOF4−: X=O, one lone pair and two basal 3c4e units [1,2,1], extending the selected IOF4− partition. C78T, pp. 537–538 compares the three anions; E02, pp. 1995–1997 supports their square-pyramidal structures and central lone pair. HW 545 |
| L07/M13, IO2F4− | I#c-#v8 | Present; existing state | cis-/trans-IO2F4−: two I=O and four localized I–F bonds under the selected no-lone-pair convention. C81, p. 2109 supports I=O and distinguishes both isomers; HW 546. One fixture covers their shared valence graph |
| D01/D10, Br/I | X#c+#n3#m0 (X = Br, I); N#v2#a#m2 | Present | [X(py)2]+: R20, pp. 7983–7986 explicitly describes N–X–N 3c4e bonding; selected contributions [2,0,2]. Supersedes the two-dative-bond proposal. HW 491–492 lists Br/I examples |
| D01/D10, Cl | Cl#c+#n3#m0 | Not selected | [Cl(py)2]+: supported by R20, pp. 7983–7986 and P23, pp. 2326–2327; omitted from the default for now |
| D01/D10, At | Unassigned | Not selected | [At(py)2]+ occurs in HW 492; no extension by analogy for the default registry |
| D02 | I#n3#v#t | Present | Me3N→I2 charge-transfer complex; HW 494–495. N#v3#d already exists |
| D03, F/Cl/Br/I | X#c-#n3#d (X = F, Cl, Br, I) | Present | Monodentate X− donation; HW 499. Fixtures: [RhF6]3− (HW II 2012) and [RhX(NH3)5]2+ for Cl/Br/I (HW II 2015–2016); existing Rh#c+3#n3#t6 |
| D03, At | At#c-#n3#d | Not selected | HW 499 includes astatido ligands; omitted from the default scope |
| D04/D06, halide bridges | X#c-#n2#d2 (X = Cl, Br, I) | Present | Two X− donations; HW 499. Extends the selected F#c-#n2#d2 convention in Xe2F11+. Covered by I2Cl6, Br5− and I5− |
| D05, Cl | Cl#c-#n#d3; Te#c+#n#v3#t3 | Present | Te4Cl16 fixture. HW 499 gives triply bridging donors; HW 711 explicitly describes TeCl3+/Cl− cubane units |
| D05, F | F#c-#n#d3 | Present; restored 2026-11-01 | General triply bridging fluoride; [F(XeOF4)3]− and [F(SbF3)3]−, HW 499. Xe fixture restored |
| D05, Xe companion | Xe#n#v2#t#m4 | Present | XeOF4 acceptor, covered by restored [F(XeOF4)3]−. HW 471, 499; W75, p. 1477 |
| D06, I2 ligand | I#n2#v#d; Ag#c+#n5#t3 | Present | HW 489 explicitly describes I2 donation. M16, pp. 9263–9264, Fig. 1 gives discrete [Ag2(I2)4]2+: two bridging and two terminal I2 ligands, three I donors per Ag. The fixture records the additional I···I halogen contact as noncovalent |
| D08, I2Cl6 dimer | I#c+#n2#v2#t2; Cl#c-#n2#d2 | Present | Two [ICl2]+ cores joined by two doubly donating Cl− bridges, applying the settled halide-donor convention. One I acceptor state added; the Cl state exists. HW 508–509 and BW54, pp. 421–423 establish the planar dimer. The discrete encoding is our adopted assignment |
| D07, fluoride aggregates | Unassigned | Open | Additional intermolecular contacts in BrF3/IF3; the monomer states and fixtures already exist. HW 508–509 |
| D08/D11, HIO3 adduct | I#n#v5#t | Not selected | HIO3·Me2SO: selected O→I encoding retains both I=O bonds and DMSO S=O. HW 527 reports base addition; N02 depicts O attachment but does not establish this bond partition. No I acceptor entry or fixture; O#n#v2#d is now independently included for aqua donors |
| D09 | Unassigned | Open | ClO2·py: establish the attachment before assigning Cl#t; HW 538 reports the stable adduct |
| Existing atomic cations | X#c+#n3 (X = Cl, Br, I) | Present; retained | Needed for electrophilic halogenation mechanism representations. Three atomic-ion fixtures added. No bare F+ entry selected |
| Further polyhalogen cations | Unassigned | Open | Br2+/I2+, Cl4+, I4²+, Br5+/I5+, I15³+: full encodings and default inclusion; HW 490–493 gives bond order 1.5 for X2+ and multicenter links in larger cations |
| Pentahalide anions | X#c-#n2#d2; X#n3#v#t (X = Br, I) | Present | Br5− and I5− as X− donating to two X2 molecules, following HW 493–494 and the selected halide-donor convention. One Br2 acceptor state added; the other states exist |
| Higher polyhalide anions | Unassigned | Open | Polyiodides beyond I5−: individual encodings and default inclusion; HW 493–494 |
| Higher HF aggregates | Existing F#c-#n4; F#n3#v; H#v | Present; hydrogen-bond descriptions | H2F3−, branched H3F4− and H4F5− as F−·nHF, with ordinary H–F bonds and Hbd contacts to F−. HW 498, Fig. 12.8; Iv17 reports asymmetric hydrogen bonds. Three hydrogen_bonds/ fixtures replace the multicenter fixtures; the three aggregate-specific F states are removed under the scope rule |
| Further oxides | Unassigned | Open | ClO4, Cl2O2+, polymeric IO+/IO2+/I3O6+ and iodine oxides: full encodings and default inclusion; HW 531, 539, 542–544 |
| Weak associates | Unassigned | Open | Br2·toluene; F2/Cl2/Br2/ClF/BrCl with H2O, H2S or NH3: whether to include a coordination bond; HW 494–495 supplies occurrence/association evidence |

The localized X#c-#n3#v2 proposal for trihalides is removed in light of the
3-center/4-electron evidence. Other unverified localized/multicenter partitions
are replaced by the single Open row for their compound family. Halide bridges use
one X− donation convention; D06 retains only its independent neutral-I2 example.

The inferred all-dative aluminium fixtures and the Hg(OH)At fixture are withdrawn,
including their Al/O/Hg companion entries. Default inclusion remains a case-by-case
decision; the extraction may retain compounds omitted from the default.

IF7, IF8−, IOF6−, IOF5, IO2F3 and ClO2F3 use localized bonds by the selected
no-lone-pair convention: without a central lone pair to reserve s character,
no multicenter unit is assigned. These are pragmatic graph assignments;
the earlier mixed and 6c10e proposals for these compounds are superseded.

F3− (15 K matrix, decomposition at 40 K) and X2− (X = F, Cl, Br, I;
irradiated glasses/matrices) require no additional states beyond the trihalide rows;
no dedicated default additions are proposed from these examples (HW 493–494).
HW marks ClF7, BrF7, F3+ and I(OAc)7 as unestablished, and I7+, I(OSO2F)5
and XO3+ as possible or questioned. No default entries are proposed for them
(HW 491–492, 508–510, 545).

Primary sources reviewed; retained PDFs in materials/valence, extracts in scratch:

- **Iv17:** Ivlev et al., *Z. Anorg. Allg. Chem.* **643** (2017), 1436–1443;
  [DOI](https://doi.org/10.1002/zaac.201700228), [author manuscript](../materials/valence/Ivlev-2017-sodium-hydrogen-fluorides.pdf).
  Crystal structures of NaH2F3, NaH3F4 and NaH4F5, including neutron diffraction
  for NaH2F3; all F–H···F bridges are asymmetric. The fixtures use the reported
  F−·nHF description with Hbd contacts.
- **C72D:** Christe and Curtis, *Inorg. Chem.* **11** (1972), 35–39;
  [DOI](https://doi.org/10.1021/ic50107a009), [PDF](../materials/valence/Christe-1972-ClO2F2-anion.pdf#page=4).
  Pages 38–39 explicitly describe Cl=O bonds and an axial 3c4e F–Cl–F unit in ClO2F2−.
- **GS76:** Gillespie and Spekkens, *J. Chem. Soc., Dalton Trans.* (1976), 2391–2396;
  [DOI](https://doi.org/10.1039/DT9760002391), [PDF](../materials/valence/Gillespie-1976-bromine-oxyfluoride-anions.pdf#page=3).
  Page 2393 assigns BrO2F2− two axial F atoms and equatorial O atoms and lone pair from its Raman spectrum.
- **C78T:** Christe et al., *Inorg. Chem.* **17** (1978), 533–538;
  [DOI](https://doi.org/10.1021/ic50181a003), [PDF](../materials/valence/Christe-1978-XOF4-anions.pdf#page=5).
  Pages 537–538 compare the square-planar XF4 groups and X–O/X–F force constants in ClOF4−, BrOF4− and IOF4−.
- **C81:** Christe, Wilson and Schack, *Inorg. Chem.* **20** (1981), 2104–2114;
  [DOI](https://doi.org/10.1021/ic50221a034), [PDF](../materials/valence/Christe-1981-tetrafluoroperiodates.pdf#page=6).
  Page 2109 supports I=O bonds in IO2F4− and identifies both isomers by vibrational and 19F NMR spectra.
- **M16:** Malinowski, Himmel and Krossing, *Angew. Chem. Int. Ed.* **55** (2016), 9262–9266;
  [DOI](https://doi.org/10.1002/anie.201603913),
  [PDF](../materials/valence/Malinowski-2016-diiodine-coordination.pdf#page=2).
  Figure 1 establishes terminal/bridging I2 ligands in discrete Ag complexes;
  p. 9264 explicitly describes I2→Ag donation.
- **BW54:** Boswijk and Wiebenga, *Acta Cryst.* **7** (1954), 417–423;
  [DOI](https://doi.org/10.1107/S0365110X54001260),
  [PDF](../materials/valence/Boswijk-1954-I2Cl6-structure.pdf#page=6).
  Planar I2Cl6 dimer with two chlorine bridges; the proposed ionic resonance
  description does not settle the discrete bond partition.
- **R20:** Reiersolmoen et al., *Chem. Sci.* **11** (2020), 7979–7990;
  [DOI](https://doi.org/10.1039/D0SC02076A),
  [PDF](../materials/valence/Reiersolmoen-2020-three-center-halogen-bonds.pdf).
  Pages 7983–7986 explicitly describe 3c4e bonding in the Cl/Br/I bis(pyridine) cations.
- **P23:** Pröhm et al., *Chem. Sci.* **14** (2023), 2325–2329;
  [DOI](https://doi.org/10.1039/D2SC06757A),
  [PDF](../materials/valence/Prohm-2023-pyridine-chloronium.pdf#page=2).
  Pages 2326–2327 establish [Cl(py)2]+ in the solid state and retain its 3c4e description.
- **N02:** Nicolaou, Montagnon and Baran, *Angew. Chem. Int. Ed.* **41** (2002), 1386–1389;
  [DOI](https://doi.org/10.1002/1521-3773%2820020415%2941%3A8%3C1386%3A%3AAID-ANIE1386%3E3.0.CO%3B2-V),
  [PDF](../materials/valence/Nicolaou-2002-HIO3-DMSO.pdf).
  Scheme 1B depicts HIO3·DMSO; the crystallographic evidence on p. 1388 concerns I2O5·DMSO.
- **S77:** Smart, *J. Chem. Soc., Chem. Commun.* (1977), 519–520;
  [DOI](https://doi.org/10.1039/C39770000519),
  [PDF](../materials/valence/Smart-1977-IO2F3-dimer.pdf).
  Establishes the oxygen-bridged IO2F3 crystal dimer, superseding GK76's trimer
  assignment for the solid; no explicit F-bond partition.
- **C93A:** Christe et al., *J. Am. Chem. Soc.* **115** (1993), 2696–2706;
  [DOI](https://doi.org/10.1021/ja00060a019),
  [PDF](../materials/valence/Christe-1993-IOF6-anion-bonding.pdf#page=9).
  Pages 2704–2706 propose the equatorial 6c10e model for IOF6−; this partition is
  not adopted. They also describe the selected basal multicenter bonding in IOF4−
  and IF5.
- **T12:** Thanthiriwatte et al., *Inorg. Chem.* **51** (2012), 10966–10982;
  [DOI](https://doi.org/10.1021/ic301438b),
  [PDF](../materials/valence/Thanthiriwatte-2012-halogen-oxyfluorides.pdf).
  Structural and energetic comparisons include IO2F3 and IOF5; no explicit
  multicenter partition for either. Predicted species are distinguished from
  experimentally characterized compounds.
- **C73D:** Christe and Curtis, *Inorg. Chem.* **12** (1973), 2245–2251;
  [DOI](https://doi.org/10.1021/ic50128a007),
  [PDF](../materials/valence/Christe-1973-ClO2F3-bonding.pdf#page=7).
  Pages 2250–2251 explicitly describe localized equatorial bonds and one axial
  3c4e unit in ClO2F3. The selected no-lone-pair convention supersedes this partition.
- **GK76:** Gillespie and Krasznai, *Inorg. Chem.* **15** (1976), 1251–1256;
  [DOI](https://doi.org/10.1021/ic50160a002),
  [PDF](../materials/valence/Gillespie-1976-IO2F3-structure.pdf#page=4).
  Pages 1254–1255 propose an oxygen-bridged cyclic trimer in the solid, melt and
  solution; S77 subsequently establishes a dimer in the crystal. No explicit F-bond
  partition.
- **C93O:** Christe, Curtis and Dixon, *J. Am. Chem. Soc.* **115** (1993), 9655–9658;
  [DOI](https://doi.org/10.1021/ja00074a034),
  [PDF](../materials/valence/Christe-1993-IOF5-structure.pdf#page=2).
  Comparable axial/equatorial I–F lengths and force constants in IOF5; no explicit
  multicenter partition. The IOF6− comparison on p. 9657 distinguishes its mainly
  covalent axial I=O/I–F bonds from its more ionic equatorial bonds.
- **OC82:** Oberhammer and Christe, *Inorg. Chem.* **21** (1982), 273–275;
  [DOI](https://doi.org/10.1021/ic00131a050),
  [PDF](../materials/valence/Oberhammer-1982-ClOF3-structure.pdf#page=3).
  Page 275 explicitly describes Cl=O, an equatorial Cl–F single bond, one lone pair
  and axial 3c4e bonding in ClOF3.
- **C78:** Christe, Curtis and Bougon, *Inorg. Chem.* **17** (1978), 1533–1539;
  [DOI](https://doi.org/10.1021/ic50184a029),
  [PDF](../materials/valence/Christe-1978-BrOF3-bonding.pdf#page=6).
  Pages 1537–1538 support Br=O and explicitly explain the weaker axial Br–F bonds
  using 3c4e bonding; the equatorial Br–F bond is stronger.
- **ET74:** Edwards and Taylor, *J. Fluorine Chem.* **4** (1974), 173–179;
  [DOI](https://doi.org/10.1016/S0022-1139%2800%2982511-7),
  [PDF](../materials/valence/Edwards-1974-IOF3-structure.pdf#page=3).
  Pages 175, 178–179 establish IOF3 molecules with weaker intermolecular contacts;
  the oxygen position is reassigned. No explicit 3c4e description.
- **CA71:** Carter and Aubke, *Inorg. Chem.* **10** (1971), 2296–2301;
  [DOI](https://doi.org/10.1021/ic50104a041),
  [PDF](../materials/valence/Carter-1971-iodine-oxyfluorides.pdf#page=4).
  Page 2299 distinguishes equatorial and axial I–F stretches in IOF3 and identifies
  one iodine lone pair; no explicit 3c4e assignment.
- **E02:** Ellern et al., *Z. Anorg. Allg. Chem.* **628** (2002), 1991–1999;
  [DOI](https://doi.org/10.1002/1521-3749%28200209%29628%3A9%2F10%3C1991%3A%3AAID-ZAAC1991%3E3.0.CO%3B2-1),
  [PDF](../materials/valence/Ellern-2002-halogen-oxyfluoride-structures.pdf#page=2).
  Pages 1992–1995 confirm ClOF3/BrOF3 structures and their relation to IOF3;
  pp. 1995–1996 establish square-pyramidal BrOF4− with one bromine lone pair.
- **C68:** Christe and Sawodny, *Z. Anorg. Allg. Chem.* **357** (1968), 125–133;
  [DOI](https://doi.org/10.1002/zaac.19683570303),
  [PDF](../materials/valence/Christe-1968-ClF4-anion.pdf#page=8).
  Pages 125 and 132–133 support two 3c4e units in ClF4−.
  Pages 132–133 also apply the 3c4e model to the basal bonds of ClF5 and
  distinguish its more covalent apical bond.
- **CN73:** Christe and Naumann, *Inorg. Chem.* **12** (1973), 59–62;
  [DOI](https://doi.org/10.1021/ic50119a016),
  [PDF](../materials/valence/Christe-1973-IF4-anion.pdf#page=4).
  Table III (p. 61) compares ClF4−, BrF4− and IF4−; the discussion on p. 62
  interprets all three using semiionic 3c4e bonds and identifies two central lone pairs.
- **CW89:** Christe and Wilson, *Inorg. Chem.* **28** (1989), 3275–3277;
  [DOI](https://doi.org/10.1021/ic00316a001),
  [PDF](../materials/valence/Christe-1989-BrF6-IF6-lone-pairs.pdf#page=2).
  Pages 3276–3277 establish octahedral BrF6− and distorted IF6− in solution,
  both fluxional on the NMR timescale, with inactive/active lone pairs respectively.
  The paper does not specify a 3c4e bond partition.
- **C73:** Christe, Wilson and Curtis, *Inorg. Chem.* **12** (1973), 1358–1362;
  [DOI](https://doi.org/10.1021/ic50124a029),
  [PDF](../materials/valence/Christe-1973-difluoroperchloryl-cation.pdf#page=5).
  Page 1362 depicts the two-Cl=O resonance form and identifies predominantly
  covalent Cl–F bonds in ClO2F2+.
- **C01:** Christe, Zhang, Sheehy and Bau, *J. Am. Chem. Soc.* **123** (2001), 6338–6348;
  [DOI](https://doi.org/10.1021/ja003347a),
  [PDF](../materials/valence/Christe-2001-tetrafluorohalogen-cations.pdf#page=10).
  Page 6347 attributes weaker axial bonds to strong 3c4e contributions in
  ClF4+/BrF4+/IF4+; the discrete registry partition above is the selected assignment.
- **C93:** Christe, Curtis and Dixon, *J. Am. Chem. Soc.* **115** (1993), 1520–1526;
  [DOI](https://doi.org/10.1021/ja00057a044),
  [PDF](../materials/valence/Christe-1993-IF7-bonding.pdf#page=6).
  Pages 1525–1526 propose a six-center/ten-electron equatorial system and two
  more covalent axial bonds in IF7. This partition is not adopted for IF7.
- **D08:** Dixon et al., *Inorg. Chem.* **47** (2008), 5485–5494;
  [DOI](https://doi.org/10.1021/ic800021h),
  [PDF](../materials/valence/Dixon-2008-iodine-fluorides.pdf#page=9).
  Pages 5486–5488 identify experimentally known ClF6−/BrF6− and calculate Oh
  minima, versus C3v for IF6−. Page 5493 attributes iodine-fluoride bond-length
  differences to 3c4e contributions and discusses the central lone pairs;
  three separate 3c4e units in XF6− are our adopted assignment.
- **L97:** Landrum, Goldberg and Hoffmann, *J. Chem. Soc., Dalton Trans.* (1997),
  3605–3613; [DOI](https://doi.org/10.1039/A703736H),
  [PDF](../materials/valence/Landrum-1997-trihalide-bonding.pdf).
  Table 1 and the mixed-trihalide analysis support electron-rich three-center bonding.
- **Le04:** Lehmann et al., *Inorg. Chem.* **43** (2004), 6905–6921;
  [DOI](https://doi.org/10.1021/ic040015o),
  [PDF](../materials/valence/Lehmann-2004-hexafluorohalogen-cations.pdf#page=13).
  The bonding discussion on pp. 6917–6918 and conclusion on p. 6919 support six
  localized X–F bonds in XF6+ (X = Cl, Br, I), hence X#c+#v6.
- **M22:** Möbs et al., *Chem. Eur. J.* **28** (2022), e202202466;
  [DOI](https://doi.org/10.1002/chem.202202466),
  [PDF](../materials/valence/Mobs-2022-BrF6-IF6-bonding.pdf#page=12).
  Pages 10–11 and Figs. 13/15 identify a central lone pair in both BrF6− and IF6−;
  the former is 4s-like, the latter predominantly 5s with a small 5p contribution.
- **S20:** Scheibe et al., *Eur. J. Inorg. Chem.* (2020), 4483–4496;
  [DOI](https://doi.org/10.1002/ejic.202000845),
  [PDF](../materials/valence/Scheibe-2020-difluorochloronium.pdf#page=12).
  Page 4493 distinguishes localized and 3-center/4-electron Cl–F bonding and
  supports localized ClF2+; the remaining fluoride assignments are listed above.
- **Sh73:** Shaw, Carroll and Thomas, *J. Am. Chem. Soc.* **95** (1973), 5870–5875;
  [DOI](https://doi.org/10.1021/ja00799a006),
  [PDF](../materials/valence/Shaw-1973-ClF3-photoelectron.pdf#page=4).
  Pages 5873–5875 explicitly describe equatorial 2c2e and axial 3c4e bonding in
  ClF3, supported by F 1s binding energies, charge analysis and the MO discussion.

#### Registry and ValenceTable comparison

All 55 retained added states are non-redundant under the registry raise defaults and have
conformance coverage: 60 new fixtures and revised ClF3/IF5 fixtures. Both AtomTyping
policies produce the selected assignments in all 62 fixtures, including the IO2F3
monomer and crystal dimer.

No ValenceTable extension is needed for the localized assignments, including
I#c-#v8 and the separate At− oxoanions. Counts MostSaturated agrees with AtomTyping
on all 32 localized fixtures and the three HF-aggregate hydrogen-bond fixtures.
Both models and policies agree on the latter without ties, preserving all Hbd
contacts. Counts Strict leaves H2F+ underdetermined in the
two fixtures with unspecified H counts; MostSaturated selects the intended H2F+.

The 27 multicenter fixtures now agree under both models and policies after the
[Counts correction](240-validation-constraint-ops-2026-10-18.md#counts-electron-accounting)
on 2026-11-02.
For BrO with both H and spin unspecified, Counts Strict remains underdetermined and
MostSaturated adds H; the BrO fixture supplies H=0.

The donor–acceptor pass retains 19 non-redundant states after the 2026-11-01 fluoride
restoration. Thirteen fixtures cover all 19:
Br/I bis(pyridine), Me3N→I2, Te4Cl16, [Ag2(I2)4]2+, I2Cl6,
Br5−, I5−, [RhF6]3−, [RhX(NH3)5]2+ for Cl/Br/I and restored [F(XeOF4)3]−.
Both AtomTyping policies produce the selected assignments. The terminal At− donor entry is omitted.
No ValenceTable change is made for this pass.
The three retained bare Cl+/Br+/I+ states have separate atomic-ion fixtures;
both AtomTyping and Counts resolve each as X#c+#n3 under both policies.

### Chalcogens

Source: [HW 103, Band I, chapter XIII](<../materials/valence/Holleman-Wiberg - Anorganische Chemie - 103. Auflage - Band 1.pdf#page=591>),
printed pp. 553–738. References below use printed pages; descriptions use
AtomTypeRegistry::raise_defaults. Extraction covers the chapter; the selected
implementation and audit of existing entries are complete (2026-11-02).

Forty chalcogen states and four N/Br/I companions are added. One duplicate
O#c+#n#h3 is removed, leaving 100 unique chalcogen states: O 29, S 22, Se 22, Te 20,
Po 7. F#c-#n#d3 is also restored under the bridge-scope revision above.

#### Chalcogen work table

**Present** means implemented; **Accepted** means selected with implementation work
remaining; **Open** means the assignment is unresolved; **Not selected** means omitted
from the default; **Deferred** means held outside the default pending reconsideration.
Citation or fixture work on an accepted entry does not reopen its
inclusion. Row IDs are retained; unassigned families have no proposed atom description.

Charged diatomics admit both atom-localized and bond-based charge/spin descriptions.
Resolution preserves the supplied placement and fills missing atom fields. The
bond-based fixtures use the neutral E2 double-bond reference with bond #c±#s2;
the integer order is a graph encoding, not the fractional MO bond order.

| ID | Atom descriptions | Status | Examples; evidence or remaining work |
| --- | --- | --- | --- |
| Ch01 | O#n2#v2; O#n2#v#h; O#n2#h2; O#c-#n3#v; O#c-#n3#h; O#c-2#n4; O#c+#n#v3; O#c+#n#v2#h; O#c+#n#v#h2; O#c+#n#h3 | Present | Ordinary O, hydroxide, peroxide and hydronium states; HW 565–567, 572–573, 590–600. Duplicate hydronium removed, OX2H+ comment corrected; O2 fixture supplies triplet bond spin |
| Ch02 | O#n#v2#d; O#n#v#h#d; O#n#h2#d; O#c-#n2#v#d; O#c-#n2#h#d | Present; 5 added | Terminal aqua and hydroxido donors. HW I 572; II 2015 gives [Rh(H2O)6]3+ and [Rh(OH)(H2O)5]2+; all five states covered |
| Ch02, bridging neutral O | O#v2#d2; O#v#h#d2; O#h2#d2 | Deferred | Neutral O bridging two acceptors: [water](https://doi.org/10.1107/S2056989017005564), [methanol](https://doi.org/10.1016/S1387-7003(03)00143-6), [ether](https://doi.org/10.1016/j.poly.2010.04.036). Two donated pairs, no remaining lone pair under the adopted donor convention. Available abstracts/structural excerpts reviewed; PDFs not obtained. Default inclusion deferred for scope and per-O matching cost |
| Ch03 | E#n2#v2; E#n2#v#h; E#n2#h2; E#c-#n3#v; E#c-#n3#h; E#c-2#n4 (E = S, Se, Te) | Present; Te#c-#n3#h added | Chains, hydrides and chalcogenide ions; HW 615–628, 699–704. HTe− has explicit and implicit H fixtures |
| Ch04 | E#n#v4; E#v6 (E = S, Se, Te) | Present; 2 Te states added | Ordinary oxo convention. TeO3²−, TeF6 and Te(OH)6: HW 703, 710, 723–724. Solid H2TeO3 is not a structural example: HW leaves its structure unknown |
| Ch05 | E#n#v2#m2 (E = S, Se, Te) | Present; 3 added | EF4: two equatorial single bonds and one axial 3c4e unit [1,2,1]. HW 633, 707, 710–711; C01, p. 6347. Three fixtures |
| Ch06 | E#c-#n#v#m4 (E = S, Se, Te) | Present; 3 added | EF5−: one lone pair, one localized axial E–F bond and two basal 3c4e units [1,2,1]. HW 633, 707, 711; C72P, p. 1682. Three fixtures |
| Ch07 | E#c+#n#v3 (E = S, Se, Te); S#c+#v5 | Present; S#c+#v5 added | Pyramidal EX3+ (HW 633, 707, 711); tetrahedral SOF3+ with S=O and three S–F bonds (HW 664) |
| Ch08 | Se#c+#n#v3#t3; Br#c-#n#d3; I#c-#n#d3 | Present; 3 added | Binary Se4Cl16/Te4Br16/epsilon-Te4I16 fixtures. HW 707, 711 explicitly gives [EX3+X−]4 cubanes; extends the existing Te4Cl16 donor convention |
| Ch09 | Po#n2#v2; Po#n2#v#h; Po#n2#h2; Po#c+2#n2; Po#c+4#n | Present; 5 added | H2Po (HW 704), Po2+/Po4+ solution and salt chemistry (HW 725); five fixtures |
| Ch10 | O#c+#n#u#v2; O#n2#v#m; O#n#v2#m2; existing O#n2#v2, O#n2#u#v, O#c-#n3#v | Present; 3 added | HW 565–567, Fig. 13.6. Four O2± fixtures cover atom-localized and bond-charge/spin encodings. O3/O3− use two O–O singles and a multicenter π system with contributions [1,2,1]; O3− has system charge −1 and doublet spin (PW). Two new fixtures; both models agree under both policies |
| Ch11, Zundel | O#n#v2#m2; O#n#v#h#m2; O#n#h2#m2; O#n2#v2#m0; O#n2#v#h#m0; O#n2#h2#m0; existing H#c+#m0 | Present; 5 added | Both 3c2e (M17) and 3c4e (HW 174–176, 311, 593–594) descriptions accepted. O–H–O contributions [2,0,0] retain one/two lone pairs on the O atoms; [2,0,2] retains one on each (PW). Six fixtures cover explicit, implicit and mixed terminal H |
| Ch11, hydration | Existing H2O, H3O+, OH− and H5O2+ states | Present; no new states | H7O3+ and H9O4+: hydronium plus two/three waters; H13O6+: existing Zundel core plus four waters; H3O2−: localized hydroxide/water description. HW 593–594, Fig. 13.17. Four hydrogen_bonds/ fixtures use Hbd for the hydration contacts and preserve those relations under resolution |
| Ch12 | S#n2#u#v; S#c-#n3#v; S#n2#v#m; S#n#v2#m2; E#n2#v2 and E#c+#n#v3 (E = S, Se, Te); E#n#v2#a2 and E#c+#n#v2#a (E = S, Se, Te) | Present; selected slice complete | S2− has atom-localized and bond-charge/spin fixtures (HW 622; PW). S3− has localized S−–S–S• and multicenter π [1,2,1] with system charge −1/doublet spin (HW 622, ozonide analogy; PW); terminal S#n2#v#m added. Localized S8²⁺/Se8²⁺ use six divalent E and two trivalent E+ atoms (M12). Te6²⁺ has six ring singles and one transannular single in one contributor (Be94, Fig. 4); Te6⁴⁺ has two Te3 triangles and two inter-ring singles (Bu82). HW assigns 6π aromaticity to E4²⁺ (621, 700–701) and 10π to S6²⁺ (621); three missing Se/Te states added. S5+/S6+/S7+ use ring singles and 9/11/13-electron multicenter systems with system charge +1 and doublet spin (KP04; allocation: PW). All have fixtures; both models agree |
| Ch13 | Te#c-#v7; Te#c-2#v8 | Present; 2 added | TeF7−/TeOF5− and TeF8²−/TeOF6²−. HW 710, 725; C91H/C93H/C93T. Selected no-central-lone-pair / ordinary-oxo convention (PW); four fixtures |
| Ch14, hexahalides | E#c-2#n#m6 (E = Se, Te) | Present; 2 added | SeBr6²−, TeCl6²− and TeBr6²−: three trans 3c4e units [1,2,1], one central lone pair. HW 707, 711–712; A07, p. 18; N14, §4.2.1; contributions/charge placement: PW. Three fixtures; both models agree |
| Ch14, tetrahalides | E#c-2#n2#m4 (E = Se, Te) | Present; 2 added | SeX4²− (X = Cl, Br), TeX4²− (X = Cl, Br, I): two trans 3c4e units [1,2,1], two central lone pairs. HW 708, 712; explicit 3c4e assignments: K84, p. 1637; K91, pp. 23–24; P85, pp. 254–255. Contributions/charge placement: PW. SeBr4²− and TeCl4²− fixtures; both models agree |
| Ch15 | S#u#v#m4; existing E#n2#v2 (E = O, S, Se, Te) | Present; 1 added | SF5: axial single bond, two basal 3c4e units [1,2,1], one unpaired S electron and no lone pair (HW 634; C93R, pp. 1133–1134; contributions: PW). SO/SeO/TeO use a double bond with #s3 (HW 646, 714–715; NIST). Four fixtures; both models agree. Existing EH radical citations completed in Ch20 |
| Ch16, E4N4 | E#c+#n#v3; N#c-#n2#v2 (E = S, Se); existing HW-contributor states | Present; existing states | Localized: eight E–N singles and two E–E singles (M12). Earlier S4N4 HW contributor with two N–S–N 3c4e units retained |
| Ch16, E4N4 whole cage | E#c+#n#v2#m; N#c-#n#v2#m2 (E = S, Se) | Not selected | Optional 12-electron, eight-center description (AP24; Se analogy and contributions: PW). Three entries and two fixtures removed 2026-11-02: simpler descriptions already cover these molecules |
| Ch16, small nitrides | S#c+#n#v3; N#n#v3; N#c+#v4; S#n2#v2; S#n#v2#m2; N#n#v2#m | Present; existing states | SN+: S+≡N; S2N+: S=N+=S (HW 683). S4N2: six-membered ring, six singles and one N–S–N 3c4e π unit [1,2,1] (HW 675, 682; contributions: PW). Three fixtures; both models agree |
| Ch16, S2N2 | S#n#v2#a2; N#n#v2#a; S#n#v2#m2; N#n#v2#m | Present; existing states | Two alternative encodings of four S–N sigma bonds and six π electrons: aromatic (K18), or nonaromatic multicenter [2,1,2,1]#u2#s1 carrying HW's two singlet-coupled radical electrons (HW 681; allocation: PW). One lone pair per atom; radical sites unspecified. Both have fixtures; both models agree |
| Ch16, aromatic nitride ions | E#n#v2#a2; E#c+#n#v2#a (E = S, Se); N#n#v2#a; N#c-#n#v2#a2 | Present; N− state added | HW 683–686, 726–727: S2N3+, S3N2²⁺ and Se3N2²⁺ (6π); S4N3+, S4N4²⁺ and S3N3− (10π); S5N5+ (14π). Seven fixtures, one lone pair per atom; formal-charge placement: PW |
| Ch16, radical nitride cations | E#n#v2#m2 (E = S, Se); N#n#v2#m | Present; existing states | E3N2+: five-membered E–N–E–E–N ring, five sigma singles, seven-electron π system [2,1,2,2,1]#c+#s2. HW 683–684, 726–727; contributions/system charge: PW. Two fixtures; both models agree |
| Ch16, Se4N2 | Unassigned | Not selected; reassigned | HW 725–726 repeats the original formula. Siivari et al. (1992), pp. 1518–1519, identify the material as Se3N2Cl2 by elemental analysis, reaction stoichiometry and independent synthesis with matching IR; full paper reviewed. No Se4N2 fixture or registry example |
| Ch16, acyclic SxN− | S#n2#v#m; S#n#v2#m2; N#n#v2#m | Present; existing states | S3N−: S–S–N–S, three singles and six π electrons [1,2,1,1]#c− (B82). S4N−: S–S–N–S–S, four singles and eight π electrons [1,2,1,2,1]#c− (C80, pp. 5780–5781). HW 685; graph allocation: PW. Two fixtures; both models agree |
| Ch16, SN2²− | S#n#v4; N#c-#n2#v2 | Present; existing states | N−=S=N−: HW 686; V18, Fig. 1a and Table 2. One fixture; both models agree |
| Ch17 | Unassigned | Not selected in this pass | Larger polyatomic clusters and specialized adducts: HW 620–622, 699–703, 707–713, 725–727. Se2N2/SeSN2 occur as adducts, not established free-molecule examples. No states added solely for elaborate ions/complexes |
| Ch18 | Unassigned | Not selected; unsupported | PoF6: HW 706 lists it, but W60 establishes only a volatile Po fluoride without composition determination |
| Ch19, oxenium | O#c+#n2#v | Present; cited and covered | Singlet phenyloxenium (H13, Scheme 2 and p. 9080) replaces the OH+ registry example. One fixture uses the existing state; both models agree |
| Ch19, Te examples | Te#c+#n#v3 | Present; H3Te+ fixture restored | TeF3+ (HW 711) and gas-phase H3Te+ (KJ86, pp. 117–118) have separate fixtures using the same existing state. HW 705 reports H2Te decomposition in HF/SbF5; the full papers distinguish this failed salt preparation from gas-phase protonation |
| Ch20, remaining existing entries | Existing atomic, aromatic and remaining hydride/radical states | Present; audit complete | Atomic configurations: HW 328; OH, H3S+, H3Se+: HW 584, 625, 705; SH/SeH/TeH: NIST; heterocycles: I98/Q15. Twenty-two fixtures complete coverage of all 100 retained chalcogen states under both AtomTyping policies. Furan has aromatic and Kekulé fixtures. No further registry states or scope removals |

**H3Te+/H3Se+ evidence correction (2026-11-02).** Both full papers reviewed.
[Minkwitz, Kornath and Sawodny (1992)](../materials/valence/Minkwitz-1992-selenonium-hexafluoroantimonate.pdf),
doi:[10.1002/anie.199206431](https://doi.org/10.1002/anie.199206431), pp. 643–644:
isolated H3Se+SbF6− and its deuterated analogue, characterized by IR/Raman spectra;
stable for about a week at 195 K, decomposing within minutes at 213 K.
H2Te decomposes in HF/SbF5 at 190 K; no telluronium salt was isolated.
KJ86: [Karpas and Jasien (1986)](../materials/valence/Karpas-1986-hydrogen-telluride-proton-affinity.pdf),
doi:[10.1016/0168-1176(86)87046-X](https://doi.org/10.1016/0168-1176(86)87046-X), pp. 117–118:
gas-phase proton-transfer bracketing and equilibria establish PA(H2Te) = 175.8 ± 0.5 kcal/mol,
above H2Se (171.3); calculated H3Te+ has H–Te–H = 93.6° and Te–H = 1.654 Å.
This supports restoring the explicit-H H3Te+ fixture alongside TeF3+, using the existing
Te#c+#n#v3 state; both resolvers and both policies agree without ties.
H3Se+ retains its four explicit/mixed/implicit-H fixtures.

**Deferred O donor check (2026-11-02).** A trial registry with these three entries
changed no candidate sets, resolved molecules or reports across 906 fixtures under
Strict and MostSaturated. All 556 O atoms have complete donated-pair counts of 0 or 1.
The neutral-O bucket grows from 17 to 20 entries. A small release-build comparison
(20 alternating measured pairs after warm-up; parsing, input cloning and startup
excluded) gave median MostSaturated resolution increases of 4.5% for water, 5.0% for
[Rh(H2O)6]3+ and 0.8–1.4% for methanol, dimethyl ether and 18-crown-6. These are
fixture-specific measurements, not a general workload estimate. The default is unchanged.

Potential optimization: reject incompatible candidate families earlier using complete
structural readings, such as donated-pair counts that exclude all #d2 entries.
Any approach must preserve partial-input semantics, admitted candidates and tie-break
outcomes; the mechanism remains an open design question.

#### Primary sources

C01: [Christe et al. (2001)](../materials/valence/Christe-2001-tetrafluorohalogen-cations.pdf),
doi:[10.1021/ja003347a](https://doi.org/10.1021/ja003347a), p. 6347, Table 22 and discussion.
C72P: [Christe et al. (1972)](../materials/valence/Christe-1972-chalcogen-pentafluoride-anions.pdf),
doi:[10.1021/ic50113a046](https://doi.org/10.1021/ic50113a046), p. 1682, right-column bonding discussion.
W60: [Weinstock and Chernick (1960)](../materials/valence/Weinstock-1960-volatile-polonium-fluoride.pdf),
doi:[10.1021/ja01500a080](https://doi.org/10.1021/ja01500a080), pp. 4116–4117.
C91H: [Christe et al. (1991)](../materials/valence/Christe-1991-high-coordinate-anions.pdf),
doi:[10.1039/C39910000837](https://doi.org/10.1039/C39910000837), pp. 837–840; p. 839 explicitly excludes a central lone pair for TeF8²−.
C93H: [Christe et al. (1993)](../materials/valence/Christe-1993-TeF7.pdf),
doi:[10.1021/ja00074a011](https://doi.org/10.1021/ja00074a011), pp. 9461–9467.
C93T: [Christe et al. (1993)](../materials/valence/Christe-1993-TeOF6-TeOF5.pdf),
doi:[10.1021/ic00071a020](https://doi.org/10.1021/ic00071a020), pp. 4089–4093.
These six papers are available locally and have been read.

M12: [Moilanen et al. (2012)](../materials/valence/Moilanen-2012-transannular-chalcogen-interactions.pdf),
doi:[10.1021/ct300627p](https://doi.org/10.1021/ct300627p), bonding analysis and
[SI, List S2](../materials/valence/Moilanen-2012-transannular-chalcogen-interactions-SI.pdf).
The E4N4 NBO Lewis reference has E+/N− atoms, eight E–N singles, two E–E singles,
one lone pair per E and two per N. The S8²⁺/Se8²⁺ references have one transannular
E–E single; its endpoints are E#c+#n#v3, and the other six atoms are E#n2#v2.
The E8²⁺ NBO description is discussed in the main paper's Bonding analysis section.

KP04: [Krossing and Passmore (2004)](../materials/valence/Krossing-2004-sulfur-cations.pdf),
doi:[10.1021/ic0207303](https://doi.org/10.1021/ic0207303), p. 1004 and
[SI, pp. 2–3](../materials/valence/Krossing-2004-sulfur-cations-SI.pdf), reviewed.
S5+/S6+/S7+ are rings with distributed charge/spin; an explicit σ/π partition is
not supplied. The accepted PW convention uses ring singles, one lone pair per S
and a 9/11/13-electron multicenter system, respectively. Every S uses existing
S#n#v2#m2; system charge +1 and doublet spin give one unpaired electron on the
multicenter relation. Three fixtures implement this allocation.

Bu82: [Burns et al. (1982)](../materials/valence/Burns-1982-polyatomic-ion-bonding.pdf),
doi:[10.1021/ic00132a065](https://doi.org/10.1021/ic00132a065), pp. 804–805, reviewed.
For Te6⁴⁺ the authors explicitly give two bonds between the Te3 rings, represented
by three resonance forms. The fixture uses one contributor: four Te#c+#n#v3
and two Te#n2#v2 atoms. The paper's planar aromatic Te6²⁺ is a prediction,
preceding the distorted-prismatic structure in HW 700–701; it is not used.

Be95: [Beck (1995)](../materials/valence/Beck-1995-Te6-dication.pdf),
doi:[10.1002/cber.19951280104](https://doi.org/10.1002/cber.19951280104), pp. 24–25,
reviewed. Te6²⁺ has six ring Te–Te singles (270–274 pm) in a boat conformation
and two weak transannular links (about 328 pm). The paper describes two positive
charges distributed over four Te atoms; the salt's paramagnetism belongs to W(V).

Be94: [Beck (1994)](../materials/valence/Beck-1994-tellurium-polycations.pdf),
doi:[10.1002/anie.199401631](https://doi.org/10.1002/anie.199401631), p. 166,
Fig. 4, reviewed. Its two explicit resonance contributors each have one transannular
single, Te(2)–Te(6) or Te(3)–Te(5), with positive charges at that bond's endpoints.
The fixture uses the former: two Te#c+#n#v3 and four Te#n2#v2 atoms, all existing states.

AP24: [Al-Yassiri and Puchta (2024)](../materials/valence/Al-Yassiri-2024-S4N4-bonding.pdf),
doi:[10.1016/j.comptc.2024.114828](https://doi.org/10.1016/j.comptc.2024.114828), §3.3.
The 12-electron, eight-center model was reviewed. The E:1/N:2 allocation and Se
analogue are PW assignments, excluded from the default because simpler descriptions
are available. The earlier S4N4 fixture remains one HW resonance contributor, with two
N–S–N 3c4e π units (HW 675, 677, Figs. 13.48/13.50; contributions: PW).
[Gopinathan and Whitehead (1975)](../materials/valence/Gopinathan-1975-S4N4-localized-orbitals.pdf),
doi:[10.1139/v75-185](https://doi.org/10.1139/v75-185), pp. 1344–1347, also supports
the localized S+/N− reference.

M17: [Molčanov et al. (2017), H5O2+](../materials/valence/Molcanov-2017-Zundel-bonding.pdf),
doi:[10.1039/C7CE00501F](https://doi.org/10.1039/C7CE00501F), pp. 3898–3900.
The paper assigns 3c2e character from charge-density analysis. The fixture's
[2,0,0] contributions and asymmetric lone-pair placement are one PW allocation,
not a reported orbital decomposition. Both it and the HW-based [2,0,2] encoding
account for 16 valence electrons; neither is preferred by AtomTyping.

A07: [Aragoni et al. (2007)](../materials/valence/Aragoni-2007-three-center-bonding.pdf),
doi:[10.1155/2007/17416](https://doi.org/10.1155/2007/17416), p. 18, Conclusions (iv),
explicitly includes SeBr6²− in the orthogonal-3c4e description.
N14: [Närhi et al. (2014)](../materials/valence/Narhi-2014-hexahalotellurates.pdf),
doi:[10.1016/j.saa.2013.09.063](https://doi.org/10.1016/j.saa.2013.09.063), §4.2.1
(accepted manuscript p. 11), assigns 3c4e bonding through Te 5p orbitals and a
retained 5s lone pair for TeCl6²−/TeBr6²−.
[Nakamura et al. (1963)](../materials/valence/Nakamura-1963-hexahaloselenates.pdf),
doi:[10.1021/ic50005a018](https://doi.org/10.1021/ic50005a018), pp. 63–64,
supports p-only bonding and the retained s² pair for Se/Te hexahalides, without
an explicit three-center partition. Its inferred partial charges are not the
fixture formal charges.

K84: [Krebs, Schäffer and Pohl (1984)](../materials/valence/Krebs-1984-bromoselenates.pdf#page=5),
doi:[10.1515/znb-1984-1202](https://doi.org/10.1515/znb-1984-1202), p. 1637.
K91: [Krebs et al. (1991)](../materials/valence/Krebs-1991-chloroselenates.pdf#page=7),
doi:[10.1002/zaac.19915920103](https://doi.org/10.1002/zaac.19915920103), pp. 23–24.
Both explicitly describe linear X–Se–X 3c4e units in SeBr4²− and SeCl4²−, respectively.
P85: [Pohl, Saak and Krebs (1985)](../materials/valence/Pohl-1985-tetrahalotellurates.pdf#page=4),
doi:[10.1515/znb-1985-0218](https://doi.org/10.1515/znb-1985-0218), pp. 254–255,
explicitly assigns 3c4e bonding in TeX4²− (X = Cl, Br, I), with bond orders near 0.5.
The [1,2,1] contributions and central −2 formal charge follow the adopted PW convention.

NIST: Huber/Herzberg diatomic constants in the NIST Chemistry WebBook, SRD 69,
doi:[10.18434/T4D303](https://doi.org/10.18434/T4D303):
[SeO](https://webbook.nist.gov/cgi/cbook.cgi?ID=C12640890&Mask=1000) and
[TeO](https://webbook.nist.gov/cgi/cbook.cgi?ID=C13451177&Mask=1000).
Both have X³Σ− ground states. The fixtures use the established O2 graph convention:
two lone pairs per atom, double bond, and triplet spin on the bond (PW).
HW 646 explicitly gives SO double-bond character and a triplet ground state.

The same NIST compilation supplies the doublet ground states of
[SH](https://webbook.nist.gov/cgi/inchi?ID=C13940211&Mask=1000),
[SeH](https://webbook.nist.gov/cgi/cbook.cgi?ID=C13940222&Mask=186C) and
[TeH](https://webbook.nist.gov/cgi/inchi?ID=C13940368&Mask=1000), supporting the existing
single-bond radical descriptions. Atomic ns²np⁴ configurations use HW I, p. 328,
Table 9.1; OH and pyramidal H3S+/H3Se+ use HW I, pp. 584, 625 and 705, respectively.

I98: [Moss (1998), IUPAC fused-ring nomenclature](https://iupac.qmul.ac.uk/fusedring/app2.html),
doi:[10.1351/pac199870010143](https://doi.org/10.1351/pac199870010143), Appendix 2,
entries 71, 78 and 83: furan, thiophene and selenophene structures.
Furan retains both ordinary O#n2#v2 and aromatic O#n#v2#a2 descriptions.
Q15: [Quiñonero (2015)](../materials/valence/Quinonero-2015-pyrylium.pdf#page=2),
doi:[10.3390/molecules200711632](https://doi.org/10.3390/molecules200711632),
p. 11633, Scheme 1: pyrylium and thiopyrylium, supporting the existing E#c+#n#v2#a states.

C93R: [Christe et al. (1993)](../materials/valence/Christe-1993-fluorosulfate-radicals.pdf#page=5),
doi:[10.1021/ja00056a043](https://doi.org/10.1021/ja00056a043), pp. 1133–1134.
The EPR/LDF discussion places the SF5 unpaired electron mainly in an axial sulfur
sp orbital. Its comparison with trans-SF4O− supports two basal 3c4e units and an
axial covalent S–F bond. S#u#v#m4 has no lone pair; [1,2,1] contributions are PW.

K18: [Karadakov et al. (2018)](../materials/valence/Karadakov-2018-S2N2-aromaticity.pdf),
doi:[10.1002/chem.201804292](https://doi.org/10.1002/chem.201804292), introduction and conclusions.
The six-electron π space is consistent with aromaticity and some singlet-diradical character.
[Tuononen et al. (2005)](../materials/valence/Tuononen-2005-chalcogen-nitride-bonding.pdf),
doi:[10.1021/jp052502a](https://doi.org/10.1021/jp052502a), likewise finds mainly aromatic
S2N2 with 6% diradical character; its “2π” description counts the bonding π electrons.

[Siivari, Chivers and Laitinen (1992), “Se3N2Cl2, a Novel Selenium–Nitrogen Chloride:
Reinvestigation of ‘Se4N2’”](../materials/valence/Siivari-1992-Se4N2-reassignment.pdf),
doi:[10.1002/anie.199215181](https://doi.org/10.1002/anie.199215181),
Angew. Chem. Int. Ed. 31, 1518–1519. Full paper reviewed. Elemental analysis and
reaction stoichiometry support Se3N2Cl2; an independent synthesis gives the same
IR spectrum. The paper leaves ionic versus covalent forms of Se3N2Cl2 open.

B82: [Bojes et al. (1982)](https://doi.org/10.1021/ja00382a017), pp. 4839–4840;
author-posted full text reviewed. Explicit six-electron π system in S3N−;
reference contributions are one electron per terminal S and per N, and two per internal S.
C80: [Chivers et al. (1980)](../materials/valence/Chivers-1980-S4N-anion.pdf),
doi:[10.1021/ja00538a013](https://doi.org/10.1021/ja00538a013), pp. 5773, 5780–5781, Table V.
Full paper reviewed: eight π electrons in S4N−, two lone pairs per terminal S and
one per internal S. The S(2)···S(4) contact has essentially cancelling bonding and
antibonding contributions (p. 5781); it adds no edge to the fixture.
The fixtures retain one N lone pair and put −1 on the π system (PW).

V18: [Valjus et al. (2018)](../materials/valence/Valjus-2018-triatomic-sulfur-nitrogen-anions.pdf),
doi:[10.1139/cjc-2018-0047](https://doi.org/10.1139/cjc-2018-0047), Fig. 1a and Table 2.
N−=S=N− is an explicit Lewis contributor; the generalized S–N Wiberg index is 1.88.
The cited [Borrmann et al. (2004) crystal-structure paper](https://doi.org/10.1002/ejic.200300954)
is not yet available locally; the fixture's bonding evidence comes from V18.

H13: [Hanway et al. (2013), phenyloxenium](../materials/valence/Hanway-2013-phenyloxenium.pdf),
doi:[10.1021/ja403370k](https://doi.org/10.1021/ja403370k), pp. 9079–9080, reviewed.
Scheme 2 gives Ph–O+; the spectroscopy and calculations support the singlet
assignment. The fixture resolves O#c+#n2#v, retaining the phenyl aromatic system.

#### Electron-capacity expansion (2026-11-01)

**Approved:** retain Element::valence_capacity = 14 for Se and 16 for Te. This
follows the calibrated-bound approach in [096](096-valence-resolution-plan-2026-05-21.md):
the former Se=12 relied on the assumption that no covered compound needed 14.
The calibration examples here are binary compounds within the selected default
scope. Complex compounds alone do not justify automatic expansion; future changes
to element bounds require consultation separately from accepting registry entries.

For these descriptions, with H=0 and no aromatic contribution, the Counts ceiling
is applied to `valence_electrons − charge + v + 2t`:

| Element | Capacity | Calibration description | Electron count | Evidence |
| --- | --- | --- | --- | --- |
| Se | 12 → 14 | Se#c+#n#v3#t3 in Se4Cl16 | (6−1) + 3 + 2×3 = 14 | HW 707, 711 |
| Te | 14 → 16 | Te#c-2#v8 in TeF8²− | (6+2) + 8 = 16 | HW 710; C91H; selected encoding: PW |

Both descriptions satisfy electron conservation. These are the smallest capacities
that accommodate them, not established exhaustive chemical maxima. TeOF6²− uses the
same Te state and count (C93T).

The bounds apply across all charges/environments in Counts resolution, classification
and implicit-H inference used by export. AtomTyping matching and registry size are
unchanged. ValenceInvariants validation checks conservation without this ceiling;
its enumeration uses capacity/2 as a lone-pair bound. The 820-case conformance suite
changed only the three Counts fixture outcomes below; it is not an exhaustive audit
of newly admitted inputs.

#### Conformance and Counts comparison

One hundred and nine new fixtures and revised SF4 cover all 44 new states under both
AtomTyping policies, with identical assignments. The final audit adds 22 fixtures
using existing states; all 100 retained chalcogen states are covered under both
AtomTyping policies. All four configurations agree without ties in these 22 cases.
Aromaticity's minimum-electron
tie-break is reported for S4N4²⁺ and S5N5+ under AtomTyping, and for all six new
aromatic nitride cations under Counts; S3N3− needs no tie-break. The restored
[F(XeOF4)3]− fixture additionally covers restored F#c-#n#d3 and the existing Xe acceptor.
The explicit-H H3Te+ fixture is restored alongside TeF3+ after reviewing KJ86;
both use the existing Te#c+#n#v3 state (Ch19).

| Coverage | Fixture paths under tests/resolution/data |
| --- | --- |
| EF4 and TeF6 | hypervalent/{sf4,sef4,tef4,tef6}.edn |
| SF5 radical | hypervalent/sf5.edn |
| EF5−, SOF3+, tellurite, implicit HTe− | inorganic_ions/{sf5-,sef5-,tef5-,sof3+,teo3-2,teh-implicit-}.edn |
| Higher-coordinate Te anions; TeF3+ and H3Te+ | inorganic_ions/{tef7-,tef8-2,teof5-,teof6-2,tef3+,teh3+}.edn |
| Se/Te hexahalide dianions | inorganic_ions/{sebr6-2,tecl6-2,tebr6-2}.edn |
| Se/Te tetrahalide dianions | inorganic_ions/{sebr4-2,tecl4-2}.edn |
| Aqua / hydroxido donors, explicit and implicit H | donor_acceptor/rh-h2o6*.edn; donor_acceptor/rh-oh*.edn |
| Binary cubanes | donor_acceptor/{se4cl16,te4br16,te4i16}.edn |
| Po hydride and ions | hydrides/h2po*.edn; atomic_ions/{po+2,po+4}.edn |
| Charged diatomics, atom-localized and bond-charge/spin descriptions | inorganic_ions/{o2+-atom,o2+-bond,o2--atom,o2--bond,s2--atom,s2--bond}.edn |
| Triplet chalcogen monoxides | inorganic_small/{so,seo,teo}.edn |
| Ozone and ozonide, multicenter π descriptions | inorganic_small/o3-multicenter.edn; inorganic_ions/o3-.edn |
| Trisulfide radical anion, localized and multicenter descriptions | inorganic_ions/{s3-localized-,s3-multicenter-}.edn |
| Zundel ion, 3c2e and 3c4e descriptions with explicit, implicit and mixed terminal H | inorganic_ions/h5o2-3c2e*.edn; inorganic_ions/h5o2-3c4e*.edn |
| Hydrogen-bond compatibility using existing atom types | hydrogen_bonds/{h7o3+,h9o4+,h13o6+,h3o2-}.edn |
| S4N4 HW contributor | inorganic_small/s4n4.edn |
| Small nitrides, existing states | inorganic_ions/{sn+,s2n+}.edn; inorganic_small/s4n2.edn |
| S2N2, aromatic and singlet-diradical multicenter descriptions | aromatic/s2n2.edn; inorganic_small/s2n2-diradical.edn |
| Aromatic nitride ions | aromatic/{s2n3+,s3n2-2+,se3n2-2+,s4n3+,s4n4-2+,s5n5+,s3n3-}.edn |
| Radical nitride cations | inorganic_ions/{s3n2+,se3n2+}.edn |
| Acyclic nitride anions | inorganic_ions/{s3n-,s4n-,sn2-2}.edn |
| S4N4 / Se4N4, localized descriptions | inorganic_small/{s4n4,se4n4}-localized.edn |
| Aromatic chalcogen dications | aromatic/{s4,s6,se4,te4}-2+.edn |
| Cyclic sulfur radical cations | inorganic_ions/{s5+,s6+,s7+}.edn |
| Localized S8²⁺/Se8²⁺ | inorganic_ions/{s8,se8}-localized-2+.edn |
| Te6²⁺/Te6⁴⁺ resonance contributors | inorganic_ions/te6-localized-{2+,4+}.edn |
| Singlet phenyloxenium | functional_groups/phenyloxenium.edn |
| Existing EH radicals, implicit H | inorganic_small/{oh,sh,seh}-implicit.edn |
| Existing hydronium, sulfonium and selenonium states, implicit/mixed H | inorganic_ions/{oh3,sh3,seh3}-{implicit,mixed-h1,mixed-h2}+.edn |
| Existing hydrochalcogenide anions, implicit H | inorganic_ions/{oh,sh,seh}-implicit-.edn |
| Existing H2S/H2Se/H2Te states, mixed H | hydrides/{h2s,h2se,h2te}-mixed.edn |
| Existing selenium valences | inorganic_small/seo2.edn; hypervalent/sef6.edn |
| Existing heteroaromatic states and Kekulé notation | aromatic/{pyrylium,furan,furan-kekule}.edn |

The capacity expansion removes Counts rejections of Se4Cl16, TeF8²− and TeOF6²−.
After the [accounting correction](240-validation-constraint-ops-2026-10-18.md#counts-electron-accounting)
on 2026-11-02, Counts agrees with AtomTyping in all 111 new/revised chalcogen cases.
The forty previously incorrect donor/multicenter assignments are corrected, including
the Zundel core in H13O6+. The restored [F(XeOF4)3]− fixture also agrees. All four
hydration fixtures preserve their Hbd relations and require no registry or table changes.
N−#a2 uses the existing isoelectronic O aromatic contribution in Counts.
ValenceTable covalence targets remain unchanged. The Te dianion fixtures have H
fixed to zero and use the approved element-capacity adjustment.

#### Aromaticity default and table entries (2026-11-01)

ChemistryModel::default now uses AromaticityModel::permissive: Hückel perception
with ElementScope::Any. Python delegates to the Rust default, and both conformance
valence models use this same default. The named Daylight and MDL presets retain
their element scopes and rules.

The living ValenceTable adds aromatic_valences = [1] for Sb and [2] for Te.
Counts uses the isoelectronic Sb row for Te⁺, so both additions are required to
infer Te4²⁺ from aromatic membership alone. They also apply to other inputs selecting
those element rows. The Te4²⁺ fixture resolves to contributions [1,2,1,2] and one
lone pair per atom under both models and both tie-break policies. Frozen SMILES/MDL
tables, saturation targets and electron-capacity bounds are unchanged by this step.
The E4²⁺ fixtures place the two positive charges on opposite atoms. With those
charges fixed, the two one-electron contributors are nonadjacent and have no
Kekulé matching. Kekulization owns independent success and no_matching fixture groups;
its corpus error test asserts NoMatching explicitly. The resolution corpus is unchanged.
Other test-organization and generator corrections are tracked in
[189](189-property-test-organization-2026-08-06.md#graph-test-follow-ups-2026-11-01).

### Pnictogens

Source: [HW 103, Band I, chapter XIV](<../materials/valence/Holleman-Wiberg - Anorganische Chemie - 103. Auflage - Band 1.pdf#page=777>),
printed pp. 739–985; chapter references pp. 986–993. First extraction table:
2026-11-02. Implemented after discussion: N 52, P 44, As 28, Sb 31, Bi 19 states
(95 additions). The selected entries, parameter changes and fixtures are implemented;
the other extracted cases are explicitly unselected. No selection remains open in
this chapter pass.

#### Pnictogen work table

**Present** identifies types that preceded this pass.
**Implemented** identifies selected additions and fixture work now in the repository.
**Proposed** identifies additions for discussion; **Open** identifies an unsettled
assignment, evidence or default-scope question. **No addition proposed** retains
the extraction result without expanding the default. State lists use the concrete,
zero-valence registry defaults; each value of E denotes a separate element entry.
HW references below are printed pages in Band I. PW identifies an adopted graph
allocation or extension of an existing convention, not an explicit HW assignment.

| ID | Atom descriptions / gap | Disposition | Examples, evidence and fixture candidates |
| --- | --- | --- | --- |
| Pn01 | Existing ordinary N, N+, N−, N²− and N³− states | Present | N2, NH3, N2H4, N2H2, NH4+, NH2−, NH²−, N³−, nitrogen halides and oxo compounds; HW 747–758, 768–781, 788–845. NH2/NO/NF2 and localized NO2 radicals also use existing states (799, 806, 811). Citations checked; existing fixtures retained |
| Pn02 | N#n#u2#v; N#n#u2#h | Implemented: 2 | Triplet NH, XN and NOH; HW 779 explicitly distinguishes singlet/triplet NH, 800 gives triplet XN, 821 triplet NOH. Existing N#n2#v / N#n2#h retain the singlet. Fixtures: explicit/implicit NH with and without supplied spin |
| Pn03 | Existing N/O and accepted Pn02 states | Implemented: two fixtures; no further type | NO−: HW 806 gives a double-bond reference; 822 explicitly gives triplet spin. Selected fixtures: N−=O with bond #u2#s3 uses existing types; N(2 unpaired electrons)–O− uses accepted N#n#u2#v and existing O#c-#n3#v. Both resolve under AtomTyping and Counts, Strict and MostSaturated. Charge/spin allocations are PW. N#c-#n#u2#v2 is not selected |
| Pn04 | N#v2#h#d; N#v#h2#d; N#h3#d | Implemented: 3 | Complete the existing N#v3#d family for ammine/amine donors; HW 757–758. Fixtures: the same ammine complex with explicit, mixed and implicit N–H bonds |
| Pn05 | N#c-#n#v2#d; N#c-#n#v#h#d; N#c-#n#h2#d | Implemented: 3 | Terminal amido donors, retaining one lone pair. HW 757–758, Fig. 14.6b; formal N− donor convention follows the selected anionic-ligand convention. Fixtures use the Be(NH2)2 motif, with all H representations |
| Pn06 | N#c-#v2#d2; N#c-#v#h#d2; N#c-#h2#d2 | Implemented: 3 | Bridging amido donors: two donated pairs, one to each acceptor. HW 757–758, Fig. 14.6c, [R2AlNR2]2. Fixtures use the Al amido dimer motif with R = H on N; Al acceptor attributes are supplied explicitly |
| Pn07 | N#c-2#n#v#d2; N#c-2#n#h#d2; N#c-2#v#d3; N#c-2#h#d3; N#c-3#n#d3; N#c-3#d4 | Implemented: 6 | Imido: two donated pairs and one remaining lone pair, or three donated pairs and no lone pair; explicit/implicit H variants. Nitrido: three donated pairs and one lone pair, or four donated pairs and none. HW 757–758, Fig. 14.6d–k explicitly describes these donor counts, including Cl3VNR, [Cl3VN]− and [{RuCl4(H2O)}2N]3−. Formal N²−/N³− donor convention follows the selected anionic-ligand convention. Fixtures cover the V imido/nitrido, Al imido and Ru nitrido motifs; metal attributes are supplied explicitly. All fit the N electron capacity of 8; no neutral-N bucket additions |
| Pn08 | Existing N#n#v3, N#c+#v4, N#c-#n2#v2 and aromatic N/N− states | Present; N5− fixtures implemented | Localized azide, HN3, H2N3+, N5+, N2F+, NOF3 and N2O5; HW 746–747, 773–774, 800, 813, 827. S17 establishes pentazolate rings in CsN5 by diffraction/Raman evidence and explicitly assigns aromaticity, extending HW 748. Paired aromatic/Kekulé N5− fixtures use existing types; localized charge placement is PW |
| Pn09 | Optional additional N multicenter states | No addition proposed | HW explicitly describes two 3c4e π units in azide (747), but existing localized states suffice. N2−/N2²− (747) belong with the charged-diatomic graph checks; their existence does not justify a new family of multicenter N types |
| Pn10 | Bi#n#v3; Bi#n#v2#h; Bi#n#v#h2; Bi#n#h3 | Implemented: 4; P/As/Sb counterparts present | BiH3, molecular BiX3 and RnBiH3−n; HW 949–954, 979–980. Fixtures: BiH3 with explicit/mixed/implicit H and BiCl3. Existing E#n#v3 also covers localized P4/As4, E4O6 and simple E–E skeletons (849, 944, 964, 971) |
| Pn11 | E#c+#v3#h; E#c+#v2#h2; E#c+#v#h3; E#c+#h4 (E = As, Sb); Bi#c+#v4 | Implemented: 9 | AsH4+/SbH4+ salts: HW 951; complete the H variants of existing As/Sb#c+#v4. BiR4+ is established in 980–981; use an organobismuth example, not hypothetical BiH4+ |
| Pn12 | E#c-#n2#v#h; E#c-#n2#h2 (E = As, Sb); Bi#c-#n2#v2 | Implemented: 5 | AsH2−/SbH2−: HW 951. Bi2²− has a Bi=Bi bond (948), and R2Bi− salts are described in 979. P hydride-anion variants already exist; no BiH2− example is assumed |
| Pn13 | E#c-2#n3#v (E = P, As, Sb, Bi); P#c-2#n3#h | Implemented: 5 | Singly bonded termini in E2⁴−/chain pentelides: HW 861, 948; PH²− in Li2PH: 871. Fixtures: P2⁴−, PH²−, and As2⁴−/Sb2⁴−/Bi2⁴− fragments from CaAs/Ca11Sb10/Ca11Bi10 (948) |
| Pn14 | Bi#c-3#n4 | Implemented: 1 | Ionic bismuthides M3Bi; HW 948. Completes the existing N/P/As/Sb E³− family |
| Pn15 | P#n#u2#v; P#n#u2#h | Implemented: 2 | Triplet PH is explicit in HW 866, Table 14.9 footnote a; triplet PR intermediates in 937. The existing P#n#u#v2 family covers PR2 radicals (870, 936) |
| Pn16 | P#v3#d; P#v2#h#d; P#v#h2#d; P#h3#d | Implemented: 4 | PH3, PR3, PF3 and PCl3 donors; HW 871, 888, 936. Fixtures cover PH3 H variants and PMe3 in Cr(CO)5 adducts, with Cr acceptor attributes supplied explicitly |
| Pn17 | E#v3#d; E#v2#h#d; E#v#h2#d; E#h3#d (E = As, Sb) | Implemented: 8 | Arsane/stibane donors and their hydride/organic variants; HW 951, 979. BiPh3 donation is also reported, but described as rare (979): no Bi donor addition proposed at this stage |
| Pn18 | As#v5; Sb#v5; Bi#v5; existing P#v5 | Implemented: 3 | Monomeric EF5, E(V) organics and localized oxo descriptions; HW 953, 958–959, 966, 980–981. Continues the existing P#v5/no-central-lone-pair convention (PW). Use the molecular phase; solid BiF5 is polymeric. No alternative #m encoding proposed |
| Pn19 | P#v4#h; P#v3#h2; P#v2#h3 | Implemented: 3 | HPF4, H2PF3, H3PF2 (886, 891), plus P–H oxo compounds H3PO2/H3PO3 (903, 905–908). Completes the explicitly supported H variants of P#v5; does not add PH5 |
| Pn20 | E#c-#v6 (E = P, As, Sb, Bi); P#c-#v5#h; P#c-#v4#h2 | Implemented: 6 | PF6−/PCl6−, AsF6−/AsCl6−, SbF6−/SbCl6− and BiF6−; HW 891–892, 958–959. HPF5− and trans-H2PF4− contain P-bound H and form by fluoride addition to HPF4 and H2PF3 (891); explicit-H graphs use the accepted P#c-#v6. Localized, no-lone-pair continuation (PW), also applicable to Sb(OH)6− (972). Counts capacity is exactly 12; no table-target extension needed with H supplied |
| Pn21 | Sb#c-2#v7; Bi#c-2#v7 | Implemented: 2 and capacity 12 → 14 | SbF7²−/BiF7²− are established pentagonal bipyramids (958–959). The adopted localized convention requires 14 electrons. Sb/Bi valence_capacity is raised from 12 to 14; rationale below |
| Pn22 | P#n#v3#t; P#v5#t; Sb#v5#t | Implemented: 3 | Me3N→PCl3 (888), Me3N→PF5 / Me3P→PF5 (891), O-bound SO2→SbF5 (958): explicit Lewis-acid/base adduct descriptions in HW. One accepted pair; Counts capacities suffice. Fixtures: Me3N→PCl3, Me3N→PF5 and SO2→SbF5, using existing N/O donor types |
| Pn23 | E#c+#v4#t2 (E = P, Sb, Bi); existing bridging halide donor types | Implemented: 3; assignment PW | P2Cl10 (892), cyclic (SbF5)4 and polymeric BiF5 (953–954): four terminal halides and two bridge contacts per E. Extends the selected halide-anion donor convention; no extra multicenter types. Counts electron count is 12. Fixtures cover P2Cl10, (SbF5)4 and a BiF5 chain fragment cut at boundary F atoms. Ordinary phosphonium inputs retain the same candidates and resolution results: their complete #t0 reading excludes #t2 before tie-breaking |
| Pn24, P/As tetrahalides | P#c-#n#v2#m2; As#c-#n#v2#m2 | Implemented (accepted for now): 2 | PF4−: C94P, p. 2857 explicitly assigns axial 3c4e bonding; PCl4−/AsCl4−: W02, p. 473 supports the same partition by NBO analysis. Two equatorial single bonds, one lone pair, one axial 3c4e unit; contributions [1,2,1] and central charge follow the selected convention (PW). Existing terminal F/Cl types suffice. Counts electron count: 8. Fixtures: PF4−, PCl4− and AsCl4−; HW 887–889, 955–957 |
| Pn24, asymmetric PBr4− | Accepted P#n#v3#t (Pn22); existing Br#c-#n3#d | Implemented: fixture; no further type | [Et4N][PBr4]: C94P, p. 2858 explicitly describes PBr3 with a loosely attached Br−, unlike the more symmetric [nPr4N] salt. Use Br−→PBr3 under the selected halide-donor convention; no extra PBr4-specific type |
| Pn24, SbCl5²− | Sb#c-2#n#v#m4 | Implemented (accepted for now): 1 | WJ74, p. 1680 explicitly supports basal 3c4e bonding in square-pyramidal SbCl5²−. One ordinary apical bond, two basal 3c4e units and one lone pair; contributions [1,2,1] and central charge follow the selected convention (PW). Counts electron count: 8. Fixture: SbCl5²−; HW 957. RC72 establishes SbF5²− geometry and a lone pair, but not this partition. No SbF5²− multicenter fixture is selected from that evidence alone |
| Pn24, SbCl6³−/SbBr6³− | Sb#c-3#n#m6 | Implemented (accepted for now): 1 | SJ73, pp. 212–213 explicitly applies Rundle's polyhalide bonding model to both hexahalides: Sb 5p bonding and a retained 5s lone pair. Encode three trans 3c4e units; contributions [1,2,1] and central charge follow the selected convention (PW). Existing terminal Cl/Br types suffice. Counts electron count: 8. Fixtures: SbCl6³− and SbBr6³−; the chloride structure is measured in SJ73; HW 957 |
| Pn24, BiCl6³− | Bi#n#v3#t3 | Implemented: 1 | Lo17, manuscript p. 4 explicitly permits BiCl3 with three coordinated Cl−. Use three ordinary Bi–Cl bonds and three Cl− donations, retaining one Bi lone pair; the existing Cl− donor type suffices. This donor description meets the multicenter inclusion rule without another Bi multicenter state. Counts electron count: 14, covered by the accepted Bi expansion in Pn21. Fixture: BiCl6³−; HW 957 |
| Pn24, SbCl4−/BiCl4− | Sb#c+3#n#t4; Bi#c+3#n#t4 | Implemented: 2 donor states | W02, p. 473 favors E³⁺ surrounded by four Cl− ions over the P/As 3c4e description. Accepted four terminal Cl−→E³⁺ bonds apply the existing halide-donor convention (PW); existing Cl#c-#n3#d suffices. Central electron count 2 + 4×2 = 10, net molecular charge −1. Both fixtures resolve identically under both models/policies and satisfy the invariants; no bound/table changes or new multicenter states |
| Pn24, other halogenoanions | No further states selected | No addition proposed | HW 955–957 and WK92 describe further EX4−/EX5²−/EX6³− structures and aggregates without establishing additional selected bond partitions. Retain the literature-supported multicenter assignments above; do not extrapolate new Sb/Bi multicenter types across the family. Larger aggregate ions do not justify additional states |
| Pn25 | E#n#v2#a (E = P, Sb, Bi); N/As counterparts present | Implemented: 3 | Phosphinine, stibinine, bismine; HW 938, 985. Paired aromatic/Kekulé fixtures. Counts companion: Bi aromatic contribution 1 added for inputs supplying only aromatic membership |
| Pn26 | E#c-#n#v2#a2 (E = P, As, Sb, Bi); N counterpart present | Implemented: 4 | P5− (883) and aromatic EC4R4− rings (938, 985–986). Localized charge on E is the graph allocation (PW). Paired aromatic/Kekulé fixtures. Counts companion: Po aromatic contribution 2 added for Bi−, which reads the isoelectronic chalcogen row |
| Pn27 | P#c+#n#v2#a0 | Implemented: 1 | Phosphirenylium, explicitly called aromatic in HW 937–938, Fig. 14.91a. Two C π electrons, P+ retains its lone pair. Paired aromatic/Kekulé fixtures. Counts companion: aromatic contribution 0 added to the Si row used for P+ |
| Pn28 | P#v3#a2; P#v2#h#a2; existing As#v3#a2 / As#v2#h#a2 | Implemented: 2; retain existing As | EC76, pp. 4365, 4368–4370 supports aromatic descriptions of phospholes and arsoles, including pyramidal geometries. M02, pp. 1208, 1211–1213 finds weak, substituent-dependent phosphole aromaticity; pyramidality alone does not classify the series. Qualifies HW 937–938, 985. Both aromatic and localized notations selected, with P–H explicit/implicit and P–Me fixtures and corresponding As coverage. The 6π graph allocation uses E contribution 2 (PW). Counts already admits P/As contribution 2. No extension to Sb/Bi proposed |
| Pn29 | Existing P#n#v3, P#v5, P#c+#v4 and N neutral/anionic states | Present | PN, phosphazanes, iminophosphoranes and cyclophosphazenes; HW 927–934. Localized and charge-separated contributors use existing states. No cage/ring-wide multicenter type proposed merely to reproduce delocalization |
| Pn30 | As#c+#n#v2; Sb#c+#n#v2; Bi#c+#n2 | Implemented: 3 | Accepted AsF2+/SbF2+ cation states: two ordinary bonds and one lone pair (955); EF3-stabilized aggregates require their own topology. Bi+ occurs in [Bi][Bi9][HfCl6]3 (947); the accepted two-pair allocation is PW. Three fixtures resolve identically under both models/policies; no bounds/table changes |
| Pn31 | All 79 earlier N/P/As/Sb/Bi states | Audited: retained | HW references filled throughout; AsR2 radicals: HW 983. G18 supplies explicit NBO support for Sb#n#u#v2. H81 supports singlet PH2+/HPF+/PF2+. F18 supports the singlet nitrenium entries; the unrelated NO+ example is removed. Z22 distinguishes pyrrolyl π and σ states: the existing N#u#v2#a2 example is qualified as the σ-radical B(2A1) reading, not the ground state. Corrected PH2 and R3PH+ formula comments; no type changes from this audit |
| Pn32 | Additional cluster-, aggregate- or coordination-specific states | No addition proposed | Larger P/As/Sb/Bi cages, HP4+, polynuclear halide ions, higher-coordinate nitrido complexes and Bi hydrolysis clusters: HW 855–863, 947–948, 957, 976. Many ordinary E–E skeletons already reuse the types above; do not add a multicenter state for every cluster size |
| Pn33 | Gas-phase atomic/polyatomic ions; unestablished hydrides | No addition proposed | Mass-spectrometer Nn+/Pn+ series (746, 860) do not by themselves qualify for the default. PH5, AsH5, SbH5, BiH5 are not established examples (865, 951); no extra H variants inferred from their formulas |

Ninety-five additions are implemented under the concrete, zero-valence defaults.
All have AtomTyping fixture coverage and satisfy the scalar electron balance.
The 134 fixtures include paired aromatic/Kekulé inputs and explicit, mixed and
implicit H variants. In fixtures of pnictogen donors, transition-metal acceptor
attributes are supplied; this pass adds no transition-metal states.

Counts comparison: the two Pn21 states and Pn24's Bi#n#v3#t3 required the
implemented Sb/Bi capacity of 14. The living ValenceTable now admits Bi aromatic
contribution 1, Po contribution 2 (for Bi−), and Si contribution 0 (for P+),
covering Pn25–Pn27 when only aromatic membership is supplied. Saturation targets
and frozen tables are unchanged; targets are not hard maxima.

**Implemented capacity expansion (Pn21, 2026-11-02).** Sb/Bi valence_capacity
is raised from 12 to 14 to admit localized SbF7²−/BiF7²−: 5 − (−2) + 7 = 14 electrons.
As in the Se/Te expansion, established binary compounds justify the general bound;
it is not restricted to these two registry entries. Counts already uses the I/At
saturation targets for Sb²−/Bi²−, so no ValenceTable target change is needed.

For Pn02, the earlier triplet omission avoided ambiguity before the tie-break was
developed. MostSaturated now prefers the existing singlet when spin is unspecified;
explicit triplet fixtures supply #u2#s3 and resolve under both models/policies.
Counts likewise needs the supplied unpaired count for this reading. No resolver
changes are part of this pass.

The focused existing benchmark comparison (30 samples, 1 s warm-up, 2 s measurement)
measured pyridine admission at 5.65 µs (+3.5%). Full pyridine resolution was
24.43 µs, within the benchmark's noise threshold. Octane and the charged case
showed no detected full-resolution regression. These three inputs do not establish
performance across the registry.

#### Pnictogen primary sources

- **C94P:** Christe et al., “Tetrafluorophosphite, PF4−, Anion,” *J. Am. Chem. Soc.* **116** (1994), 2850–2858. [DOI](https://doi.org/10.1021/ja00086a018); [PDF](../materials/valence/Christe-1994-PF4-anion.pdf).
- **W02:** Walker et al., “The Thermochemistry of Group 15 Tetrachloride Anions,” *J. Am. Soc. Mass Spectrom.* **13** (2002), 469–476. [DOI](https://doi.org/10.1016/S1044-0305%2802%2900369-0); [PDF](../materials/valence/Walker-2002-tetrachloride-anions.pdf).
- **Lo17:** Loor et al., “Synthesis of Sb2Se3 and Bi2Se3 Nanoparticles in Ionic Liquids at Low Temperatures and Solid State Structure of [C4C1Im]3[BiCl6],” *Z. Anorg. Allg. Chem.* **643** (2017), 60–68. [DOI](https://doi.org/10.1002/zaac.201600325); [accepted manuscript](../materials/valence/Loor-2017-BiCl6.pdf).
- **WJ74:** Wismer and Jacobson, “Crystal Structure of Potassium Pentachloroantimonate(III),” *Inorg. Chem.* **13** (1974), 1678–1680. [DOI](https://doi.org/10.1021/ic50137a026); [PDF](../materials/valence/Wismer-1974-K2SbCl5.pdf).
- **RC72:** Ryan and Cromer, “Structures of the MF5E Type. The Crystal Structure of Ammonium Pentafluoroantimonate(III),” *Inorg. Chem.* **11** (1972), 2322–2324. [DOI](https://doi.org/10.1021/ic50116a006); [PDF](../materials/valence/Ryan-1972-ammonium-SbF5.pdf).
- **WK92:** Wheeler and Kumar, “Stereochemically Active or Inactive Lone Pair Electrons in Some Six-Coordinate, Group 15 Halides,” *J. Am. Chem. Soc.* **114** (1992), 4776–4784. [DOI](https://doi.org/10.1021/ja00038a049); [PDF](../materials/valence/Wheeler-1992-group15-lone-pairs.pdf).
- **SJ73:** Schroeder and Jacobson, “Crystal Structure of Hexaamminecobalt Hexachloroantimonate(III),” *Inorg. Chem.* **12** (1973), 210–213. [DOI](https://doi.org/10.1021/ic50119a046); [PDF](../materials/valence/Schroeder-1973-SbCl6.pdf).
- **S17:** Steele et al., “High-Pressure Synthesis of a Pentazolate Salt,” *Chem. Mater.* **29** (2017), 735–741. [DOI](https://doi.org/10.1021/acs.chemmater.6b04538); [preprint](../materials/valence/Steele-2017-pentazolate-preprint.pdf). Preprint pp. 7–10: aromatic N5− rings, diffraction and Raman evidence for CsN5; no new N type required.
- **EC76:** Epiotis and Cherry, “On the Aromaticity of Phospholes and Arsoles,” *J. Am. Chem. Soc.* **98** (1976), 4365–4370. [DOI](https://doi.org/10.1021/ja00431a003); [PDF](../materials/valence/Epiotis-1976-phospholes-arsoles.pdf).
- **M02:** Mattmann et al., “De-aromatizing Phosphole,” *J. Org. Chem.* **67** (2002), 1208–1213. [DOI](https://doi.org/10.1021/jo0108156); [PDF](../materials/valence/Mattmann-2002-phosphole.pdf).

- **F18:** Falvey, “Discrete Existence of Singlet Nitrenium Ions Revisited: Computational Studies of Non-Aryl Nitrenium Ions and Their Rearrangements,” *ACS Omega* **3** (2018), 10418–10432. [DOI](https://doi.org/10.1021/acsomega.8b01038); [full-text XML](../materials/valence/Falvey-2018-nitrenium.xml). Tables 1–2 distinguish the singlet entries from triplet ground states.
- **G18:** Ganesamoorthy et al., “From Stable Sb- and Bi-Centered Radicals to a Compound with a Ga=Sb Double Bond,” *Nat. Commun.* **9** (2018), 87. [DOI](https://doi.org/10.1038/s41467-017-02581-2); [PDF](../materials/valence/Ganesamoorthy-2018-Sb-radicals.pdf). NBO analysis on p. 4 supports two bonds, one lone pair and one unpaired electron at Sb.
- **H81:** Harrison, “Electronic Structure of the Phosphenium Ions PH2+, HPF+, and PF2+,” *J. Am. Chem. Soc.* **103** (1981), 7406–7413. [DOI](https://doi.org/10.1021/ja00415a002); [PDF](../materials/valence/Harrison-1981-phosphenium.pdf). Singlet states with two bonds and a lone pair, pp. 7410–7412.
- **Z22:** Zhang, Yuan and Wang, “Probing the Electronic Structure and Spectroscopy of Pyrrolyl and Imidazolyl Radicals Using High-Resolution Photoelectron Imaging of Cryogenically Cooled Anions,” *Phys. Chem. Chem. Phys.* **24** (2022), 6505–6514. [DOI](https://doi.org/10.1039/D2CP00189F); [accepted manuscript](../materials/valence/Zhang-2022-pyrrolyl.pdf). Section 3.1.3 assigns the observed B(2A1) pyrrolyl state to σ-electron detachment; the registry's atom-local allocation is PW.

Unavailable review: Fisher and Norman, “The Structures of the Group 15 Element(III)
Halides and Halogenoanions,” *Adv. Inorg. Chem.* **41** (1994), 233–271;
[DOI](https://doi.org/10.1016/S0898-8838%2808%2960173-7), cited by HW I, p. 990.
The user supplied the abstract but also cannot access the full chapter. The abstract
establishes structural families, not a bond partition; it adds no accepted assignment.
