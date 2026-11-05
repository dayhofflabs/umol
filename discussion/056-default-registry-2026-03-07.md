# 056 — Default atom-type registry

Status: In Progress
Date: 2026-03-07
Relates: [217](217-rhea-participant-failures-2026-08-30.md),
[240](240-validation-constraint-ops-2026-10-18.md),
[189](189-property-test-organization-2026-08-06.md)

This document tracks expansion of the default atom-type registry through source review,
explicit selection and conformance coverage. The [HW extraction](#hw-extraction-spike-2026-10-31)
contains the inclusion rules, evidence and detailed chapter tables.

## Project Status (2026-11-04)

**Selected scope implemented: group chapters through Ni/Pd/Pt, chapter XXXII,
and chapters XXXV–XXXVI (lanthanides and actinides).** Every added state has source
citations and conformance coverage. Unselected and deferred candidates remain in
the linked tables.
Companion states belong to elements outside the chapter being reviewed.

| Chapter / group | Registry | Conformance fixtures | Model changes and scope notes |
| --- | --- | --- | --- |
| [Hydrogen](#hydrogen) | 5 H additions and 5 B/C/F/Cr/Rh companions; all 9 retained H states covered | 7 new | Terminal hydride donors and selected 3c2e/3c4e descriptions; ValenceTable unchanged |
| [Noble gases](#noble-gases) | All 34 retained states and 3 F/N companions covered | Selected assignments covered, including xenon-anion replacements | Perxenate uses the accepted Xe#v8 assumption. Rn compounds, neutral Xe#m8 and XeO4 adducts excluded |
| [Halogens](#halogens) | 110 retained states after citation and scope audit | 60 valence, 13 donor–acceptor and 3 atomic-ion fixtures; 2 revised | HF aggregates use hydrogen bonds. General halide donors retained; At donor extensions and specialized adducts excluded |
| [Chalcogens](#chalcogen-work-table) | 40 additions and 4 N/Br/I companions; 100 unique chalcogen states after later companion work and deduplication | 109 new; SF4 revised | Capacities Se 14, Te 16. Zundel alternatives retained; larger water clusters use hydrogen bonds. Optional whole-cage nitride types excluded |
| [Pnictogens](#pnictogen-work-table) | 95 additions: N 17, P 23, As 16, Sb 22, Bi 17; 79 earlier states retained | 134 new | Sb/Bi capacity 14; aromatic Counts entries for Bi, Po and Si. SbCl4−/BiCl4− use the selected ionic/donor description |
| [Tetrels](#tetrel-work-table) | C 12, Si 28 and Ge/Sn/Pb 93 additions; Al and S companions | C 13, Si 30, Ge/Sn/Pb 98 new | Sn/Pb capacity 14; aromatic Counts entries for Al, Ge and Ga |
| [Triels](#triel-work-table) | B 27, Al 37, Ga 39, In 32, Tl 22 additions; 4 C companions and the pyridine-donor N state | B 32, heavy triels 137 new; earlier borolide correction covered separately | Al capacity 12; ValenceTable unchanged |
| [Alkaline earths](#alkaline-earth-work-table) | 76 additions: Be 20, Mg 22, Ca 12, Sr 11, Ba 11 | 77 new; 3 Be-amido inputs revised | Capacities Be/Mg 12, Ca/Sr 16, Ba 18; Ra retains 8. KMgH3 reuses existing ions; bridging water deferred |
| [Alkali metals](#alkali-metal-work-table) | 34 additions: Li 9, Na 7, K 7, Rb 5, Cs 6; 13 earlier states retained | 35 new | Capacities Li/Na/Rb 16, K 20, Cs 24; Fr retains 8. Specialized carbon-contact/Cp assemblies excluded; bridging water deferred |
| [Cu, Ag, Au](#copper-group-work-table) | 77 additions: Cu 29, Ag 18, Au 30; 11 earlier states retained, including free-ion Au³⁺ triplet | 78 new | Capacities Cu/Ag 21; Au retains 18. Ambiguous Cs2AgF6 adds no states; η²-H2 donors deferred |
| [Zn, Cd, Hg](#zinc-group-work-table) | 53 additions: Zn 19, Cd 17, Hg 17; 1 μ₃-alkoxide O companion; 7 earlier states retained | 52 new; 4 ZnMe2 inputs revised | Hg capacity 22; Zn/Cd limits unchanged. Selected hydride bridges included; Hg oxidation states above II excluded |
| [Sc, Y, La, Ac](#scandium-group-work-table) | 79 additions: Sc 28, Y 25, La 26; 8 earlier states retained | 82 new | Sc/Y/La capacity 20 after XXXII; Ac retains 18; ValenceTable unchanged. Uncertain hydrides, H2 adducts and formula-only cases remain unselected |
| [Ti, Zr, Hf](#titanium-group-work-table) | 118 additions: Ti 47, Zr 44, Hf 27; 8 earlier states retained | 135 new | Capacities Zr/Hf 20; Ti retains 18. Both ZrH2 occupations included; broad borohydride MC states and specialized assemblies excluded |
| [V, Nb, Ta](#vanadium-group-work-table) | 160 additions: V 57, Nb 48, Ta 55; 8 earlier states retained | 165 new; 1 nitride input revised | Capacity 18 and ValenceTable unchanged; no new MC types. Four uncertain Nb/Ta(II) phosphine occupations omitted |
| [Cr, Mo, W](#chromium-group-work-table) | 260 additions: Cr 72, Mo 94, W 94; 2 hydroxido O companions; 11 earlier states retained | 262 new | Capacity 18 and ValenceTable unchanged. Two Mo/W hydride-bridge MC states included; 13 W phosphine-hydride occupations remain open |
| [Mn, Tc, Re](#manganese-group-work-table) | 243 additions: Mn 62, Tc 60, Re 121; 9 earlier states retained | 242 new | Mn capacity 21; Tc/Re/Bh retain 18; ValenceTable unchanged. Four Re states support separate 3c2e hydride bridges; uncertain occupations excluded |
| [Fe, Ru, Os](#iron-group-work-table) | 210 additions: Fe 81, Ru 57, Os 68; 3 S and 1 O companions; 17 earlier states retained | 223 new; 2 Fe halide spin inputs revised; 1 source comment corrected | Capacities Fe 21, Os 20; Ru/Hs retain 18; ValenceTable unchanged. One compact hydride-bridge MC state; uncertain occupations and specialized Fe(1-norbornyl)₄ states excluded |
| [Co, Rh, Ir](#cobalt-group-work-table) | 195 additions: Co 63, Rh 68, Ir 64; 13 earlier states retained, 12 sourced and covered | 199 new; 3 supplied-spin inputs revised | Capacities Co 23, Rh/Ir 19; ValenceTable unchanged; no new MC or companion types. Existing Co#c+#n3#u#v remains source-unresolved; other unselected cases below |
| [Ni, Pd, Pt](#nickel-group-work-table) | 173 additions: Ni 59, Pd 44, Pt 69 and 1 O companion; all 7 earlier states retained and covered | 190 new; 2 bare-ion snapshots revised | Capacities Ni/Pd 20; Pt retains 18; ValenceTable unchanged. Eight local hydride-bridge MC states included; PdF₆− omitted |
| [Transition-metal classes, XXXII](#transition-metal-classes-chapter-xxxii) | 21 additions: Cr 4, Mo 2, W 4; Sc/Ti/Fe/Y/Zr/Nb/Ru/La/Hf/Ta/Os 1 each | 21 new; all added states covered | Sc/Y/La capacity 20 for bis-COT anions. ValenceTable unchanged; no new MC or companion types. Ti/Zr/Hf NMR-supported pairing accepted |
| [Lanthanides, XXXV](#lanthanide-work-table) | 232 additions, including 3 MC and 3 La companions; 34 earlier states retained and sourced | 246 new: 218 in `lanthanides/`, 25 in `organometallics/`, 3 in `multicenter_bonds/`; all added states covered | Capacities Dy 27, Ho/Tm 30, Er/Yb 31; Nd stays 26. Nd charge maximum 4, Dy unpaired maximum 6. ValenceTable unchanged; Nd crown pair excluded |
| [Actinides, XXXVI](#actinide-work-table) | 243 additions, including 6 hydride-bridge MC states; 51 earlier states retained and sourced | 243 new: 216 in `actinides/`, 21 in `organometallics/`, 6 in `multicenter_bonds/`; all added states covered | Capacities Np/Cf 27, Pu 28; No unpaired maximum 1. ValenceTable unchanged; UF₉³⁻/UF₁₀⁴⁻ unselected |

### Shared changes and verification

| Item | Current status |
| --- | --- |
| Resolution coverage | Each EDN input runs through AtomTyping and Counts with Strict and MostSaturated policies. Coverage does not require every extracted example to become a fixture |
| Counts accounting | [240](240-validation-constraint-ops-2026-10-18.md#counts-electron-accounting) implemented; 107 snapshots now give the selected donor/multicenter assignments under both policies |
| Aromaticity | Default is permissive Hückel; named Daylight/MDL presets retain their restrictions. Sb/Te aromatic Counts entries support Te⁺/Te; later chapter additions are listed above |
| Registry lookup | Element/valence grouping implemented; API and measurements [below](#registry-lookup-2026-11-02) |
| Convey-model test | Uses MostSaturated; both models omit H counts for BH3 and H2O |
| Deferred shared work | Multiple neutral-water/alcohol/ether donations (Ch02) and η²-H2 donor descriptions remain deferred |
| Verification after XXXVI | 3,602 available resolution cases, 3,594 graph unit tests, 223 graph property tests and 316 chemistry tests pass; formatting and Clippy pass. The 108 Beyer/Walter inputs recorded below were absent from this checkout |
| Sources | Retained sources are archived under `materials/valence/`; chapter source tables record access limitations. `scratch/` contains intermediates only |

The Decision Log below is the historical March record; later tracking-table updates are
dated.

### Conformance suite breakdown

Counts include the 108 Beyer/Walter chapter-4 inputs recorded below; these were absent
from the checkout used for the current verification.

| Category | Files | Description |
|---|---|---|
| `atoms/` | 108 | Free atoms (ground/excited states), H through Hs |
| `atomic_ions/` | 217 | Atomic ions and Hg₂²⁺ |
| `organometallics/` | 332 | Metal–carbon bonds, carbonyls, haptic ligands and related organometallic descriptions |
| `organoelement/` | 16 | Organic B, Si–Pb, As–Bi and Se–Po compounds |
| `organosulfur/` | 4 | Carbon–sulfur compounds; broader organic selection is being reorganized separately |
| `donor_acceptor/` | 46 | Molecular Lewis donor–acceptor adducts |
| `multicenter_bonds/` | 178 | Explicit multicenter bonding, including small hydride bridges |
| `metal_metal_bonds/` | 64 | Ordinary or dative metal–metal bonds |
| `hydrogen/` | 2 | Hydrogen molecules |
| `noble_gases/` | 6 | Noble-gas compounds outside the bonding classes above |
| `halogens/` | 54 | Halogen compounds |
| `chalcogens/` | 84 | O, S, Se, Te and Po compounds |
| `n_group/` | 100 | N, P, As, Sb and Bi compounds |
| `c_group/` | 101 | C, Si, Ge, Sn and Pb compounds |
| `b_group/` | 93 | B, Al, Ga, In and Tl compounds |
| `alkali_earths/` | 42 | Be, Mg, Ca, Sr and Ba compounds |
| `alkali_metals/` | 27 | Li, Na, K, Rb and Cs compounds |
| `cu_group/` | 68 | Cu, Ag and Au compounds |
| `zn_group/` | 40 | Zn, Cd and Hg compounds |
| `sc_group/` | 65 | Sc, Y and La compounds |
| `ti_group/` | 95 | Ti, Zr and Hf compounds |
| `v_group/` | 132 | V, Nb and Ta compounds |
| `cr_group/` | 206 | Cr, Mo and W compounds |
| `mn_group/` | 182 | Mn, Tc and Re compounds |
| `fe_group/` | 181 | Fe, Ru and Os compounds |
| `co_group/` | 179 | Co, Rh and Ir compounds |
| `ni_group/` | 156 | Ni, Pd and Pt compounds |
| `lanthanides/` | 218 | Chapter-XXXV compounds outside the bonding classes above |
| `actinides/` | 216 | Chapter-XXXVI compounds outside the bonding classes above |
| `hydrogen_bonds/` | 10 | Resolution compatibility with explicit hydrogen-bond relations |
| `hydrocarbons/` | 82 | C/H compounds |
| `functional_groups/` | 102 | Organic heteroatom groups (halides, amines, etc.) |
| `carbonyl_derivatives/` | 86 | Non-aromatic carbonyl compounds (aldehydes, ketones, acids, esters, amides, acyl derivatives) |
| `heterocycles/` | 53 | Non-aromatic O/N/S heterocycles, crown ethers |
| `aromatic/` | 105 | Aromatic carbocycles, heterocycles and ions; aromatic and Kekulé inputs |
| `stereo_tetrahedral/` | 22 | Tetrahedral stereochemistry |
| `stereo_cis_trans/` | 14 | Cis/trans stereochemistry |
| `metabolites/` | 9 | Metabolites and related natural compounds |
| `natural_products/` | 5 | Terpenes, steroids and alkaloids |
| `reactive_intermediates/` | 10 | Organic radicals, ions and related intermediates |
| **Total** | **3710** | |

### Registry lookup (2026-11-02)

Implemented: one pattern vector per element, grouped by localized valence, with
direct valence offsets. `lookup(element, valence: u8)` and
`lookup_by_element(element)` both return slices. Insertion order and duplicates
are preserved within each valence group; charge is matched within that group.
Admission, conformance and implicit-H inference use `evidenced_complete` valence.

The TOML numeric keys now denote `#v`. The reindexing preserved all 763 then-current
rows across 108 elements and their citations. Registry entries require literal
valence; TOML and `registry!` retain their zero-valence defaults. Python exposes the same lookup
names with detached atom-form results. Content hashes use each stored row once.

Verification: workspace tests, including 1083 resolution cases with unchanged
snapshots and all 223 graph property-suite cases; 3157 Python tests (2 skipped);
workspace Clippy and formatting. Reaction roundtrip, mapping and exact H omission
are tested separately; mapping and H-omission tests use explicit small registries.

Same-process release comparison against the preserved element/charge lookup:
50 samples, 500 ms warm-up and 2 s measurement per case; registry loading and
input cloning excluded. Candidate lists and resolved graphs agreed on all six
inputs. Median full-resolution times:

| Input | Charge lookup, µs | Valence offsets, µs | Time reduction |
|---|---:|---:|---:|
| Octane | 16.84 | 14.99 | 11.0% |
| C64 chain | 129.65 | 114.36 | 11.8% |
| Benzene | 19.11 | 18.16 | 5.0% |
| Pyridine | 24.24 | 22.91 | 5.5% |
| CH3+ | 1.345 | 1.351 | −0.5% |
| NH4+ | 1.418 | 1.413 | 0.4% |

Admission alone improved by 17.5–24.5% on the four neutral inputs; the two charged
single-atom cases were effectively unchanged.

## Code and test layout

### Key code files

- **Registry**: `umol-graph/config/default-registry.toml` — per-element atom types grouped by valence
- **Valence table**: `umol-graph/config/default-valence-table.toml` — Counts element data
- **Test harness**: `umol-graph/tests/resolution/resolution_suite.rs` — rstest-based conformance runner with insta snapshots
- **Test data**: `umol-graph/tests/resolution/data/<category>/<compound>.edn`
- **Snapshots**: `umol-graph/tests/resolution/snapshots/`
- **Resolver**: `umol-graph/src/ops/resolve.rs`, `umol-graph/src/ops/resolve/valence.rs`
- **Valence logic**: `umol-graph/src/ops/valence.rs`; `atom_typing.rs`, `counts.rs`, `registry.rs` and `table.rs` under `umol-graph/src/ops/valence/`

### Conformance execution notes

- Run the conformance suite with its feature flag enabled:
  - `cargo test -p umol-graph --features conformance --test resolution`
- Each snapshot records AtomTyping and Counts under Strict and MostSaturated. All 2,703
  inputs succeed under MostSaturated. 296 fail under at least one Strict configuration,
  mostly organic inputs whose hydrogen counts are left to resolution: 271 of the 278 inputs in
  `hydrocarbons/`, `functional_groups/`, `carbonyl_derivatives/`, `heterocycles/` and the two
  stereo categories.

## Tracking Table

Carbon/silicon review, Al companion, pnictogen review and chalcogen/N/F companion rows updated 2026-11-02; Br/I/Xe companion rows updated
2026-11-01; hydrogen, other noble-gas and other halogen rows updated 2026-10-31.
Group-1/2 rows updated 2026-11-03 for the approved batches; Cu/Ag/Au rows record
the completed selected scope on the same date, as do Zn/Cd/Hg and the O companion.
Sc/Y/La/Ac rows record the completed selected scope on 2026-11-03.
Ti/Zr/Hf rows record the completed selected scope on 2026-11-03.
V/Nb/Ta/Db rows record the completed selected scope on 2026-11-03.
Cr/Mo/W/Sg rows record the completed selected scope on 2026-11-04.
Mn/Tc/Re/Bh rows record the completed selected scope on 2026-11-04.
Elements with chapter-XXXII additions include that completed scope on 2026-11-04.
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
| Li | yes | yes (+1) | yes (selected) | yes (9 additions covered) | done | 12 current states; selected HW scope, citation audit and capacity 16 implemented |
| Be | yes | yes (+2) | yes (selected) | yes (20 additions covered) | done | 22 current states; approved additions, citation audit and capacity 12 implemented; remaining hydride/borohydride descriptions checked for reuse, with no further default additions proposed |
| B | yes | no | yes (selected) | yes (27 additions covered; Tr02–Tr11, Tr35) | done | 50 current states; selected boron pass and existing-entry citation audit complete; borolide charge corrected and unsupported B–H row removed |
| C | yes | yes (−4) | yes (selected) | yes (12 tetrel additions and 4 carborane companions covered) | done | 52 current states; selected tetrel/carborane scope implemented; earlier entries retained with checked citations |
| N | yes | yes (-3) | partial | yes (`atoms/n.edn`, `atomic_ions/n-3.edn`, `donor_acceptor/`, `n_group/`, `aromatic/`) | in_progress | 53 current states; 17 pnictogen additions, both NO− fixtures and Tr37's general pyridine donor implemented; existing-entry citation audit complete |
| O | yes | yes (-2) | yes (selected) | yes (all 32 states) | done | Chalcogen selection and citation audit complete; Zg24's μ₃-alkoxide and Crg37's two μ₂-hydroxido donor entries are implemented and covered |
| F | yes | yes (-1) | partial | yes (selected HW states) | done | Three higher HF-aggregate states removed; H-bond fixtures use existing types. Triply bridging F− restored 2026-11-01; terminal donor covered by [RhF6]3− |
| Ne | yes | no | n/a | yes (`atoms/ne.edn`) | done | HW review complete; zero-lone-pair description retained |
| Na | yes | yes (+1, −1) | yes (selected) | yes (7 additions covered) | done | 9 current states; selected HW scope, citation audit and capacity 16 implemented |
| Mg | yes | yes (+2) | yes (selected) | yes (22 additions covered) | done | 26 current states; additions, citation audit and capacity 12 implemented; KMgH3's ionic description reuses existing states |
| Al | yes | yes (+3) | yes (selected) | yes (37 additions covered) | done | 41 current states; selected HW scope and existing-entry citation audit complete; capacity 12 for five-/six-coordinate donor descriptions |
| Si | yes | yes (−4) | yes (selected) | yes (28 additions covered by 30 fixtures) | done | 39 current states; selected HW scope complete; all 11 earlier entries retained with checked citations; Al zero-π Counts companion added |
| P | yes | yes (-3) | partial | yes (`atoms/p.edn`, `atomic_ions/p-3.edn`, `n_group/`, `multicenter_bonds/`) | in_progress | 44 current states; 23 additions implemented; existing-entry citation audit complete |
| S | yes | yes (-2) | yes (selected) | yes (all 23 states) | done | Six states added, including the Tt27 thioether donor companion; citation and default-scope audit complete; S4N4 retains localized and HW three-center descriptions |
| Cl | yes | yes (+1, -1) | partial | yes (selected HW states) | done | 16 states added; terminal donor covered by [RhCl(NH3)5]2+ |
| Ar | yes | no | n/a | yes (`atoms/ar.edn`) | done | HW review complete; zero-lone-pair description retained |
| K | yes | yes (+1, −1) | yes (selected) | yes (7 additions covered) | done | 9 current states; selected HW scope, citation audit and capacity 20 implemented |
| Ca | yes | yes (+2) | yes (selected) | yes (12 additions covered) | done | 14 current states; additions, citation audit, seven-donor aqua fixture and capacity 16 implemented; water bridges deferred with Ch02 |
| Sc | yes | yes (+3) | yes (selected) | yes (29 additions covered) | done | 31 current states; both earlier states retained; selected HW scope and citation audit complete; capacity 20 for bis-COT donors |
| Ti | yes | yes (+2,+3,+4) | yes (selected) | yes (48 additions covered) | done | 51 current states; all 3 earlier states retained; selected HW scope and citation audit complete; capacity 18 unchanged |
| V | yes | yes (+2,+3,+4,+5) | yes (selected) | yes (57 additions covered) | done | 61 current states; all 4 earlier states retained; selected HW scope and citation audit complete; capacity 18 unchanged |
| Cr | yes | yes (+2,+3,+4,+5,+6) | yes (selected) | yes (76 additions covered) | done | 81 current states; all 5 earlier states retained; selected HW scope and citation audit complete; capacity 18 unchanged |
| Mn | yes | yes (+2 through +7) | yes (selected) | yes (62 additions covered) | done | 68 current states; all 6 earlier states retained; citation audit complete; capacity 21 implemented for high-coordinate Mn(II) |
| Fe | yes | yes (selected) | yes (selected) | yes (82 additions covered) | done | 94 current states; selected group and XXXII scopes implemented; capacity 21 |
| Co | yes | yes (+2,+3) | no | yes (`atoms/co.toml`, `atomic_ions/co+2.toml`, `atomic_ions/co+3.toml`) | done |  |
| Ni | yes | yes (+2) | no | yes (`atoms/ni.toml`, `atomic_ions/ni+2.toml`) | done |  |
| Cu | yes | yes (+1,+2) | yes (selected) | yes (29 additions covered) | done | 32 current states; selected HW scope, citation audit and capacity 21 implemented |
| Zn | yes | yes (+2) | yes (selected) | yes (19 additions covered) | done | 21 current states; selected HW scope and citation audit complete; capacity 18 unchanged |
| Ga | yes | yes (+1,+3) | yes (selected) | yes (39 additions covered) | done | 41 current states; selected HW scope and existing-entry citation audit complete |
| Ge | yes | yes (−4) | partial | yes (30 additions covered) | done | 31 current states; selected HW scope and existing-entry audit complete |
| As | yes | yes (-3) | partial | yes (`atoms/as.edn`, `atomic_ions/as-3.edn`, `n_group/`) | in_progress | 27 current states; 15 additions implemented; existing arsole types retained with paired aromatic/Kekulé fixtures; citation audit complete |
| Se | yes | yes (-2) | yes (selected) | yes (all 22 states) | done | Six states added; citation and default-scope audit complete, including localized Se4N4 |
| Br | yes | yes (+1, -1) | partial | yes (selected HW states) | done | 18 halogen-pass additions plus Br#c-#n#d3 for binary chalcogen bromides |
| Kr | yes | no | yes (selected) | yes (`atoms/kr.edn`, `multicenter_bonds/krf2.edn`, `noble_gases/krf+.edn`, `donor_acceptor/hcn-krf+.edn`) | done | HW review complete; all four retained states covered |
| Rb | yes | yes (+1, −1) | yes (selected) | yes (5 additions covered) | done | 7 current states; selected HW scope, citation audit and capacity 16 implemented |
| Sr | yes | yes (+2) | yes (selected) | yes (11 additions covered) | done | 13 current states; additions, citation audit, seven-/eight-donor aqua fixtures and capacity 16 implemented; water bridges deferred with Ch02 |
| Y | yes | yes (+3) | yes (selected) | yes (26 additions covered) | done | 28 current states; both earlier states retained; selected HW scope and citation audit complete; capacity 20 for bis-COT donors |
| Zr | yes | yes (+3,+4) | yes (selected) | yes (45 additions covered) | done | 47 current states; both earlier states retained; both ZrH2 occupations, Schwartz's reagent and selected donor states implemented; capacity 20 implemented |
| Nb | yes | yes (+3,+4,+5) | yes (selected) | yes (49 additions covered) | done | 51 current states; both earlier states retained; selected HW scope and citation audit complete; uncertain II phosphine states omitted; capacity 18 unchanged |
| Mo | yes | yes (+2,+3,+4,+5,+6) | yes (selected) | yes (96 additions covered) | done | 100 current states; all 4 earlier states retained; selected HW scope and citation audit complete; capacity 18 unchanged |
| Tc | yes | yes (+2 through +7) | yes (selected) | yes (60 additions covered) | done | 61 current states; earlier atom retained; selected scope and citation audit complete; capacity 18 unchanged |
| Ru | yes | yes (selected) | yes (selected) | yes (58 additions covered) | done | 61 current states; selected group and XXXII scopes implemented; capacity 18 |
| Rh | yes | yes (+2,+3) | partial | yes (`atoms/rh.edn`, `atomic_ions/rh+2.edn`, `atomic_ions/rh+3.edn`, Rh–H donor fixture below) | done | Hydrogen companion addition: Rh#c+3#n3#t6 |
| Pd | yes | yes (+2) | no | yes (`atoms/pd.toml`, `atomic_ions/pd+2.toml`) | done |  |
| Ag | yes | yes (+1,+2) | yes (selected) | yes (18 additions covered) | done | 22 current states; selected HW scope, citation audit and capacity 21 implemented; ambiguous Cs2AgF6 unselected |
| Cd | yes | yes (+2) | yes (selected) | yes (17 additions covered) | done | 19 current states; selected HW scope and citation audit complete; capacity 18 unchanged |
| In | yes | yes (+1,+3) | yes (selected) | yes (32 additions covered) | done | 34 current states; selected HW scope and existing-entry citation audit complete |
| Sn | yes | yes (+2, −4) | partial | yes (31 additions covered) | done | 32 current states; selected HW scope and existing-entry audit complete; capacity 14 for six-donor Sn(II) |
| Sb | yes | yes (-3) | partial | yes (`atoms/sb.edn`, `atomic_ions/sb-3.edn`, `n_group/`) | in_progress | 29 current states; 20 additions and capacity 12 → 14 implemented; existing-entry citation audit complete |
| Te | yes | yes (-2, -1, +1) | yes (selected) | yes (all 20 states) | done | Eleven states added, including aromatic Te4²⁺; citation and default-scope audit complete |
| I | yes | yes (+1, -1) | partial | yes (selected HW states) | done | 22 halogen-pass additions plus I#c-#n#d3 for binary TeI4 aggregates |
| Xe | yes | no | yes (selected) | yes (`atoms/xe.edn`, `noble_gases/xe*.edn`, `multicenter_bonds/xe*.edn`, `donor_acceptor/`) | done | All 26 retained states covered after restoring [F(XeOF4)3]− |
| Cs | yes | yes (+1, −1) | yes (selected) | yes (6 additions covered) | done | 8 current states; selected HW scope, citation audit and capacity 24 implemented |
| Ba | yes | yes (+2) | yes (selected) | yes (11 additions covered) | done | 13 current states; additions, citation audit, eight-/nine-donor aqua fixtures and capacity 18 implemented; water bridges deferred with Ch02 |
| Tl | yes | yes (+1,+3) | yes (selected) | yes (22 additions covered) | done | 25 current states; selected HW scope and existing-entry citation audit complete |
| Pb | yes | yes (+2, +4, −4) | partial | yes (32 additions covered) | done | 35 current states; selected HW scope and existing-entry audit complete, including both ionic states; capacity 14 for six-donor Pb(II) |
| Bi | yes | yes (+3, −3) | partial | yes (`atoms/bi.edn`, `atomic_ions/`, `n_group/`, `aromatic/`) | in_progress | 17 current states; 15 additions and capacity 12 → 14 implemented; existing-entry citation audit complete |
| Po | yes | yes (+2, +4, −2) | yes (selected) | yes (all seven states) | done | Five states added; citation and default-scope audit complete; PoF6 omitted for lack of composition evidence |
| At | yes | yes (+1, -1) | partial | yes (selected HW states) | done | Nine states added; terminal At− donor and [At(py)2]+ not selected |
| Rn | yes | no | none selected | yes (`atoms/rn.edn`) | done | HW review complete; unsupported compound states removed; Rn#n4 retained |
| La | yes | yes (+3) | yes (selected) | yes (30 additions covered) | done | 32 current states; both earlier states retained; includes three XXXV companions; capacity 20 for bis-COT donors |
| Ce | yes | yes (+3,+4) | yes (selected) | yes (25 additions covered) | done | 28 current states; 3 earlier states retained and sourced; Ce atom term retained and sourced to ASD |
| Pr | yes | yes (+3,+4) | yes (selected) | yes (15 additions covered) | done | 18 current states; 3 earlier states retained and sourced |
| Nd | yes | yes (+2,+3,+4) | yes (selected) | yes (17 additions covered) | done | 19 current states; 2 earlier states retained and sourced; charge maximum 4; crown pair excluded |
| Pm | yes | yes (+3) | yes (selected) | yes (8 additions covered) | done | 10 current states; 2 earlier states retained and sourced |
| Sm | yes | yes (+2,+3) | yes (selected) | yes (20 additions covered) | done | 23 current states; 3 earlier states retained and sourced; one hydride-bridge MC state |
| Eu | yes | yes (+2,+3) | yes (selected) | yes (15 additions covered) | done | 18 current states; 3 earlier states retained and sourced |
| Gd | yes | yes (+3) | yes (selected) | yes (10 additions covered) | done | 12 current states; 2 earlier states retained and sourced |
| Tb | yes | yes (+3,+4) | yes (selected) | yes (15 additions covered) | done | 18 current states; 3 earlier states retained and sourced |
| Dy | yes | yes (+2,+3,+4) | yes (selected) | yes (21 additions covered) | done | 23 current states; 2 earlier states retained and sourced; capacity 27; unpaired maximum 6 |
| Ho | yes | yes (+3) | yes (selected) | yes (11 additions covered) | done | 13 current states; 2 earlier states retained and sourced; capacity 30 |
| Er | yes | yes (+3) | yes (selected) | yes (18 additions covered) | done | 20 current states; 2 earlier states retained and sourced; capacity 31; one hydride-bridge MC state |
| Tm | yes | yes (+2,+3) | yes (selected) | yes (16 additions covered) | done | 18 current states; 2 earlier states retained and sourced; capacity 30 |
| Yb | yes | yes (+2,+3) | yes (selected) | yes (22 additions covered) | done | 25 current states; 3 earlier states retained and sourced; capacity 31 |
| Lu | yes | yes (+3) | yes (selected) | yes (16 additions covered) | done | 18 current states; 2 earlier states retained and sourced; one hydride-bridge MC state; f14 core retained |
| Hf | yes | yes (+3,+4) | yes (selected) | yes (28 additions covered) | done | 30 current states; both earlier states retained; selected HW scope and citation audit complete; capacity 20 implemented; tentative HfH3 unselected |
| Ta | yes | yes (+3,+4,+5) | yes (selected) | yes (56 additions covered) | done | 57 current states; earlier state retained; selected HW scope and citation audit complete; uncertain II phosphine states omitted; capacity 18 unchanged |
| W | yes | yes (+2,+3,+4,+5,+6) | yes (selected) | yes (98 additions covered) | done | 99 current states; earlier atom retained; selected HW scope and citation audit complete; capacity 18 unchanged; 13 unselected phosphine-hydride occupations remain open |
| Re | yes | yes (+2 through +7) | yes (selected) | yes (121 additions covered) | done | 122 current states; earlier atom retained; selected scope and citation audit complete, including four hydride-bridge states; capacity 18 unchanged |
| Os | yes | yes (selected) | yes (selected) | yes (69 additions covered) | done | 70 current states; selected group and XXXII scopes implemented; capacity 20 |
| Ir | yes | no | no | yes (`atoms/ir.toml`) | done |  |
| Pt | yes | yes (+2) | no | yes (`atoms/pt.toml`, `atomic_ions/pt+2.toml`) | done |  |
| Au | yes | yes (−1,+1,+3) | yes (selected) | yes (30 additions and existing Au³⁺ triplet covered by 31 new fixtures) | done | 34 current states; all four earlier states retained; selected HW scope and citation audit complete; capacity 18 unchanged |
| Hg | yes | yes (+2, [Hg2]2+) | yes (selected) | yes (17 additions covered) | done | 20 current states; all three earlier states retained; selected HW scope, citation audit and capacity 22 implemented; Hg(IV) unselected |
| Fr | yes | yes (+1) | no additions selected | yes (`atoms/fr.edn`, `atomic_ions/fr+1.edn`) | done | Two existing states retained with checked HW citations; capacity unchanged |
| Ra | yes | yes (+2) | ionic descriptions retained | existing atom/ion fixtures | done | 2 current states retained with checked HW citations; no additional state selected |
| Ac | yes | yes (+3) | yes (selected) | yes (3 additions covered) | done | 5 current states; 2 earlier states retained and sourced |
| Th | yes | yes (+3,+4) | yes (selected) | yes (27 additions covered) | done | 30 current states; 3 earlier states retained and sourced; three hydride-bridge MC states |
| Pa | yes | yes (+3,+4,+5) | yes (selected) | yes (20 additions covered) | done | 23 current states; 3 earlier states retained and sourced |
| U | yes | yes (+3 through +6) | yes (selected) | yes (54 additions covered) | done | 59 current states; 5 earlier states retained and sourced; three hydride-bridge MC states |
| Np | yes | yes (+3 through +6) | yes (selected) | yes (39 additions covered) | done | 44 current states; 5 earlier states retained and sourced; capacity 27 |
| Pu | yes | yes (+3 through +7) | yes (selected) | yes (33 additions covered) | done | 39 current states; 6 earlier states retained and sourced; capacity 28 |
| Am | yes | yes (+2 through +6) | yes (selected) | yes (26 additions covered) | done | 30 current states; 4 earlier states retained and sourced |
| Cm | yes | yes (+3,+4) | yes (selected) | yes (10 additions covered) | done | 13 current states; 3 earlier states retained and sourced |
| Bk | yes | yes (+3,+4) | yes (selected) | yes (12 additions covered) | done | 15 current states; 3 earlier states retained and sourced |
| Cf | yes | yes (+2,+3,+4) | yes (selected) | yes (8 additions covered) | done | 12 current states; 4 earlier states retained and sourced; capacity 27 |
| Es | yes | yes (+2,+3) | yes (selected) | yes (4 additions covered) | done | 7 current states; 3 earlier states retained and sourced |
| Fm | yes | yes (+2,+3) | yes (selected) | yes (2 additions covered) | done | 5 current states; 3 earlier states retained and sourced |
| Md | yes | yes (+2,+3) | yes (selected) | yes (2 additions covered) | done | 5 current states; 3 earlier states retained and sourced |
| No | yes | yes (+2,+3) | yes (selected) | yes (3 additions covered) | done | 5 current states; 2 earlier states retained and sourced; unpaired maximum 1 |
| Lr | yes | yes (+3) | no additions selected | existing atom/ion fixtures | done | 2 current states; 2 earlier states retained and sourced; f14 core retained |
| Rf | yes | no | no | yes (`atoms/rf.edn`) | done | Existing calculated atom occupation retained; citation checked against HW II, pp. 1538, 2358; compound chemistry deferred to XXXVII |
| Db | yes | no | no | yes (`atoms/db.edn`) | done | Existing calculated atom occupation retained; citation checked against HW II, p. 1538; compound chemistry deferred to XXXVII |
| Sg | yes | no | no | yes (`atoms/sg.edn`) | done | Existing calculated atomic occupation retained with HW citation; no additions |
| Bh | yes | no | no | yes (`atoms/bh.edn`) | done | Existing calculated atomic occupation retained with checked HW citation; no additions |
| Hs | yes | no | no | yes (`atoms/hs.toml`) | done |  |

Extend this table as new elements are considered.

### Organic compounds

Chapters refer to Beyer/Walter, 25th edition. Polycyclic, charged and stereo-defined
compounds from any chapter may also fall under the classes below.

| Class | Scope | Conformance dir | Chapters | Status | Notes |
|-------|-------|-----------------|----------|--------|-------|
| Hydrocarbons | Alkanes, alkenes, alkynes, dienes/cumulenes, cycloalkanes, unsaturated alicycles, fused/spiro/bridged bicycles, polycycles, non-aromatic annulenes | `hydrocarbons/` | 4, 5, 25 | in_progress | 82 files, including the `k33` and `k5` graph-topology inputs |
| Functional groups | Halides, azides, cyanates, thiocyanates, alcohols, thiols, ethers, thioethers, amines, nitroso, nitro, sulfoxides, hydrazines, diazo/azo; N, P, S and Si groups to be extended | `functional_groups/` | 6–10, 16 | in_progress | 102 files |
| Carbonyl derivatives | Aldehydes, ketones, carboxylic acids, esters, amides, other acid derivatives, carbonic acid derivatives | `carbonyl_derivatives/` | 12–16 | in_progress | 86 files, non-aromatic |
| Heterocycles | Non-aromatic heterocycles, crown ethers, bridged/fused | `heterocycles/` | 8, 9, 28–36 | review | 53 files |
| Aromatic | Aromatic carbocycles, heterocycles and ions | `aromatic/` | 18–25, 28–35 | in_progress | 105 files, including non-carbon rings |
| Stereochemistry | Tetrahedral and cis/trans inputs | `stereo_tetrahedral/`, `stereo_cis_trans/` | any | review | 36 files |
| Organometallics | Main-group and transition-metal organyls | `organometallics/` (noted) | 11 | todo | Created after the HW pass finishes |
| Organoelement compounds | Organic P, As, Si and B compounds | `organoelement/` (noted) | 10 | todo | Created after the HW pass finishes |
| Reactive intermediates | Carbocations, carbanions, radicals, carbenes, nitrenes, ylides | `reactive_intermediates/` | any | in_progress | 10 files; expanded in a separate reaction-mechanism step |
| Metabolites | Carbohydrates, amino acids, peptides, nucleic acid components, lipids, cofactors, pigments, vitamins, hydroxy acids | `metabolites/` | 13, 15, 17, 28, 35, 37, 38 | in_progress | 9 files; vitamins occur across chapters |
| Natural products | Terpenes, steroids, alkaloids | `natural_products/` | 26, 27, 36 | in_progress | 5 files |

Extend this table as new organic classes are added.

## Conformance Policy

Per element, add explicit resolution inputs that cover:

- Free atom query (with explicit hydrogen handling in `:config-overrides`)
- At least one bonded case when chemically meaningful
- Relevant charge states (for ion stage)
- At least one "should fail" case where useful to constrain over-broad matching

Suggested naming pattern:

- Atomic ground state: `tests/resolution/data/atoms/<element>.edn`
- Atomic valence/excited variants: `tests/resolution/data/atoms/<element>_<motif>.edn`
- Ionic ground state: `tests/resolution/data/atomic_ions/<element><charge>.edn`
- Ionic valence/excited variants: `tests/resolution/data/atomic_ions/<element><charge>_<motif>.edn`

Conventions:

- `<element>`: lowercase element symbol/name token used in file naming.
- `<charge>`: signed integer with explicit sign, e.g. `+1`, `-2` (no underscore between element and charge).
- `<motif>`: lowercase, underscore-separated descriptor; may encode occupation and optional multiplet.

Organic additions follow Beyer/Walter, 25th edition, chapter by chapter from chapter 4.
Coverage aims at all compound classes and functional groups, those of N, P, S and Si in
particular, organometallic compounds and a good selection of compounds found in nature.
Not every compound encountered is added, and a fixture need not add an atom state. These
are judgment calls that favor inclusion: the suite is feature-gated, so its size is not a
constraint. Primary literature is consulted where the book is insufficient. Proposed
fixtures, atom types and ValenceTable entries are reviewed by chapter, then added in batch.

Each compound has one input carrying all its relevant features, including aromaticity and
stereo, with hydrogens implicit. Separate Kekulé inputs are not paired with natural
compounds; Kekulé encodings of aromatic compounds and heterocycles may be added later.

## Categories

The [conformance suite breakdown](#conformance-suite-breakdown) lists current
directories. Bonding classes take precedence where they help compare related
examples: `organometallics/`, `organoelement/`, `organophosphorus/`, `organosulfur/`,
`donor_acceptor/`, `multicenter_bonds/`, `hydrogen_bonds/` and `metal_metal_bonds/`.
`donor_acceptor/` contains molecular Lewis adducts; ordinary coordination complexes
belong with their element group. Explicit Hbd relations belong in `hydrogen_bonds/`.
Remaining inorganic cases use the chapter groups, with `b_group/`, `c_group/` and
`n_group/` for groups 13–15 and `sc_group/` through `zn_group/` for transition metals.
Categories are practical browsing groups; no exclusive chemical taxonomy is implied.

The 2026-11-04 reorganization moved 2,380 inputs and their snapshots without changing
their contents or the harness. Explicit/implicit H and ordinary/donor alternatives
stay together where they exercise the same class. The existing organic directories
are being reorganized separately; `organophosphorus/` has no inputs in this pass.

The former `implicit_h=true/false` setting was the precursor to atom defaults.
EDN inputs use MoleculeDefaults with per-input `:config-overrides`; hydrogen
handling follows those defaults and explicit atom predicates.

## Current Conformance Set

Inorganic compound directories and current counts are listed in the [breakdown above](#conformance-suite-breakdown).

- `atoms/`: `h`, `he`, `li`, `be`, `b`, `c`, `n`, `o`, `f`, `ne`, `na`, `mg`, `al`, `si`, `p`, `s`, `cl`, `ar`, `k`, `ca`, `sc`, `ti`, `v`, `cr`, `mn`, `fe`, `co`, `ni`, `cu`, `zn`, `ga`, `ge`, `as`, `se`, `br`, `kr`, `rb`, `sr`, `y`, `zr`, `nb`, `mo`, `tc`, `ru`, `rh`, `pd`, `ag`, `cd`, `in`, `sn`, `sb`, `te`, `i`, `xe`, `cs`, `ba`, `tl`, `pb`, `bi`, `po`, `at`, `rn`, `la`, `ce`, `pr`, `nd`, `pm`, `sm`, `eu`, `gd`, `tb`, `dy`, `ho`, `er`, `tm`, `yb`, `lu`, `hf`, `ta`, `w`, `re`, `os`, `ir`, `pt`, `au`, `hg`, `fr`, `ra`, `ac`, `th`, `pa`, `u`, `np`, `pu`, `am`, `cm`, `bk`, `cf`, `es`, `fm`, `md`, `no`, `lr`, `rf`, `db`, `sg`, `bh`, `hs`
- `atomic_ions/`: `h+1`, `li+1`, `be+2`, `mg+2`, `al+3`, `n-3`, `o-2`, `f-1`, `na+1`, `p-3`, `s-2`, `cl-1`, `k+1`, `ca+2`, `sc+3`, `ti+3`, `ti+4`, `v+2`, `v+3`, `cr+2`, `cr+3`, `mn+2`, `mn+3`, `fe+2`, `fe+3`, `co+2`, `co+3`, `ni+2`, `cu+1`, `cu+2`, `zn+2`, `ga+3`, `as-3`, `se-2`, `br-1`, `rb+1`, `sr+2`, `y+3`, `zr+4`, `nb+3`, `mo+2`, `mo+3`, `ru+2`, `ru+3`, `rh+2`, `rh+3`, `pd+2`, `ag+1`, `ag+2`, `cd+2`, `in+3`, `sb-3`, `te-2`, `i-1`, `cs+1`, `ba+2`, `tl+1`, `tl+3`, `pb+2`, `pb+4`, `bi+3`, `po-2`, `at-1`, `la+3`, `ce+3`, `ce+4`, `pr+3`, `pr+4`, `nd+3`, `pm+3`, `sm+2`, `sm+3`, `eu+2`, `eu+3`, `gd+3`, `tb+3`, `tb+4`, `dy+3`, `ho+3`, `er+3`, `tm+3`, `yb+2`, `yb+3`, `lu+3`, `hf+4`, `pt+2`, `au+1`, `au+3`, `hg+2`, `hg+2_dimer`, `fr+1`, `ra+2`, `ac+3`, `th+3`, `th+4`, `pa+4`, `pa+5`, `u+3`, `u+4`, `u+5`, `u+6`, `np+3`, `np+4`, `np+5`, `np+6`, `pu+3`, `pu+4`, `pu+5`, `pu+6`, `pu+7`, `am+2`, `am+3`, `am+4`, `cm+3`, `cm+4`, `bk+3`, `bk+4`, `cf+2`, `cf+3`, `cf+4`, `es+2`, `es+3`, `fm+2`, `fm+3`, `md+2`, `md+3`, `no+2`, `lr+3`
- `hydrocarbons/`: `ethane`, `propane`, `butane`, `isobutane`, `pentane`, `isopentane`, `neopentane`, `hexane`, `ethene`, `propene`, `but-1-ene`, `but-2-ene`, `2-methylpropene`, `2-methylbut-2-ene`, `ethyne`, `propyne`, `but-1-yne`, `but-2-yne`, `allene`, `buta-1,3-diene`, `penta-1,4-diene`, `butatriene`, `cyclopropane`, `cyclobutane`, `cyclopentane`, `cyclohexane`, `cycloheptane`, `cyclopropene`, `cyclobutene`, `cyclopentene`, `cyclopentadiene`, `cyclohexene`, `cyclohexa-1,3-diene`, `cyclohexa-1,4-diene`, `bicyclo-1.1.0-butane`, `bicyclo-2.1.0-pentane`, `bicyclo-2.2.0-hexane`, `bicyclo-3.3.0-octane`, `bicyclo-4.3.0-nonane`, `decalin`, `spiropentane`, `spirohexane`, `spiroheptane`, `spirononane`, `bicyclo-1.1.1-pentane`, `norbornane`, `bicyclo-2.2.2-octane`, `tetrahedrane`, `prismane`, `cubane`, `adamantane`, `k33`, `k5`; [Beyer/Walter additions](#chapter-4-saturated-hydrocarbons). Total: 82 files.
- `functional_groups/`: methyl/ethyl × {fluoride, chloride, bromide, iodide, azide, isocyanide, cyanate, thiocyanate, isothiocyanate, nitrite, nitrate, thiol, selenol} (26); {propyl, isopropyl, butyl, isobutyl, tert-butyl, cyclopropyl, cyclohexyl}-chloride (7); dichloromethane, chloroform, carbon-tetrachloride, 1,1-dichloroethane, 1,2-dichloroethane, 1,2-dichlorocyclohexane, 1-chloro-2-hydroxyethane (7); methanol, ethanol, propan-1-ol, propan-2-ol, butan-1-ol, butan-2-ol, 2,2-dimethylpropan-2-ol, hydroxycyclopropane, hydroxycyclohexane (9); dimethylether, diethylether, 1,2-dimethoxyethane (3); methylamine, ethylamine, dimethylamine, trimethylamine, tetramethylammonium, trimethylamine-oxide (6); nitrosomethane, nitrosoethane, nitromethane, nitroethane, nitrosocyclohexane, nitrocyclohexane (6); dimethyl-sulfoxide, methylhydrazine, 1,1-dimethylhydrazine, 1,2-dimethylhydrazine, diazomethane, dimethyldiazene (6); ethyleneglycol, glycerol, methanediol (3); phenyloxenium (1); [Beyer/Walter additions](#chapter-4-saturated-hydrocarbons) (28). Total: 102 files.
- `carbonyl_derivatives/`: non-aromatic carbonyl compounds: aldehydes, ketones, dicarbonyls/hydroxy-carbonyls, monocarboxylic and dicarboxylic acids, acyl halides, anhydrides, carboxylates, esters/lactones/thioesters, amides/lactams/imides, and carbonyl-derived C=N species (oxime/hydrazone), plus pyruvic acid; [Beyer/Walter additions](#chapter-4-saturated-hydrocarbons) (24). Total: 86 files.
- `heterocycles/`: non-aromatic heterocycles. O-saturated (oxirane, dioxirane, oxetane, 1,2-dioxetane, 1,3-dioxetane, tetrahydrofuran, tetrahydropyran); N-saturated (aziridine, azetidine, pyrrolidine, piperidine); S-saturated (thiirane, thietane, tetrahydrothiophene, thiane); O/S-unsaturated (2H-oxete, 2H-thiete, 2H-pyran, 4H-pyran, 2H-thiopyran, 4H-thiopyran); N-unsaturated (2H-azirine, 2,3-dihydroazete, azete, 3-pyrroline, 2-pyrroline); functionalized (sulfolane, 2-oxazolidone, hydantoin); multi-O/S (1,3-dioxolane, 1,2,3-trioxolane, 1,4-dioxane, 1,3-dithiane, 1,4-dithiane, 1,3,5-trithiane); multi-N (1,3-diazetidine, pyrazolidine, imidazolidine, 2-pyrazoline, 2-imidazoline); mixed (1,2-oxazetidine, 1,2-oxathiolane, 1,3-oxathiolane, morpholine, thiomorpholine); bridged/fused (pyrrolizidine, quinuclidine, 1-azaadamantane); crown ethers (12-crown-4, 15-crown-5, 18-crown-6); cages and polyaza rings (urotropin, tetrahydrooctaazanaphthalene). Total: 53 files.

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

## Beyer/Walter extraction

Source: [Schirmeister, Schmuck and Wich, Beyer/Walter Organische Chemie, 25. Auflage](<../materials/valence/Schirmeister T. et al . - Beyer Walter Organische_Chemie - 25. Auflage.pdf>).
Page references are printed pages; PDF page = printed page + 22. Working data for each
chapter is under `scratch/bw-chNN/`.

### After the HW pass

- Make parsing and resolution settings available in conformance inputs. `:config-overrides`
  sets raising defaults only; the harness fixes the valence tie-break, running every input
  under Strict and MostSaturated. Until then, inputs keep hydrogens implicit, so most
  organic inputs do not resolve under Strict.
- Create `organometallics/` and `organoelement/`, and add the held compounds.

### Chapter 4: saturated hydrocarbons

Implemented; pages 99–127. All 108 inputs resolve under MostSaturated in AtomTyping and
Counts, so the chapter requires no new atom types or ValenceTable entries. Hydrogens are
implicit; the three hydrogen-free inputs also resolve under Strict, and the other 105 do
not. All 2,811 resolution cases pass with existing snapshots unchanged.

| Category | Count | Compounds (printed page) |
|---|---:|---|
| `hydrocarbons/` | 29 | heptane (100), octane (104), decane (100), dodecane (103), heptadecane (102), hectane (100), 2-2-dimethylhexane (104), 2-2-3-3-tetramethylbutane (104), di-tert-butylethyne (124), cyclodecane (110), cycloeicosane (111), sterane (118), norbornene (123), dicyclopentadiene (123), tetrahydrodicyclopentadiene (123), bicyclo-3.3.1-non-1-ene (120), bicyclo-4.2.1-non-1-ene (120), bicyclo-4.2.1-non-6-ene (120), twistane (120), diamantane (120), dodecahedrane (122), pagodane (125), dewar-benzene (121), benzvalene (121), 1-1-1-propellane (121), tetra-tert-butyltetrahedrane (124), 1-2-dimethylenecyclopropane (126), 3-3-dimethylene-1-1-bi-cyclobutylidene (126), trimethylene-trispirododecane (126) |
| `functional_groups/` | 28 | propyl-iodide (105), 1-bromo-3-chloropropane (110), 1-3-dibromopropane (111), 1-4-dibromobutane (111), allyl-chloride (110), methallyl-dichloride (123), cyclobutyl-bromide (113), cyclobutanol (114), cyclopropylmethanol (116), cyclobutylamine (116), cyclopropylmethylamine (116), cyclopropanone-hydrate (114), cyclohexyl-hydroperoxide (117), 1-1-2-2-tetrafluorocyclobutane (112), octafluorocyclobutane (112), tetrachlorotetrafluorocyclobutane (112), 1-3-diiodobicyclo-1.1.1-pentane (126), bicyclo-1.1.1-pentane-1-3-dithiol (126), 1-chloro-7-7-dimethylnorbornane (125), 1-nitronorbornane (126), 1-3-5-7-tetranitrocubane (127), octanitrocubane (127), pentaerythritol (122), cyclodecyl-tosylate (110), amantadine (125), memantine (125), rimantadine (125), isodrin (125) |
| `carbonyl_derivatives/` | 24 | cyclobutanone (114), cyclopentanone (114), cycloheptanone (115), cyclodecanone (110), cyclopent-2-enone (115), 2-chlorocyclohexanone (114), 2-bromocyclopentadienone (124), tetra-tert-butylcyclopentadienone (124), tetramethylcyclobutane-1-3-dione (112), 3-tert-butylfuran-2-5-dione (124), cyclobutanecarboxylic-acid (113), cyclobutane-1-1-dicarboxylic-acid (112), ethyl-cyclopentanecarboxylate (114), diethyl-malonate (111), diethyl-cyclopropane-1-1-dicarboxylate (112), diethyl-cyclobutane-1-1-dicarboxylate (113), adipic-acid (114), pimelic-acid (115), suberic-acid (115), spiro-3.3-heptane-2-6-dicarboxylic-acid (118), bicyclo-1.1.1-pentane-1-3-dicarboxylic-acid (123), 1-3-diacetylbicyclo-1.1.1-pentane (123), dimethyl-cubane-1-4-dicarboxylate (124), cubane-1-3-5-7-tetracarboxylic-acid (127) |
| `aromatic/` | 3 | ketamine (117), alpha-truxillic-acid (113), beta-truxinic-acid (113) |
| `metabolites/` | 9 | the eight inositol diastereomers, with both chiro enantiomers: cis, epi, allo, myo, muco, neo, scyllo, 1D-chiro, 1L-chiro (109–110) |
| `natural_products/` | 5 | curcumenol (107), chrysanthemic-acid (107), r-muscone (111), civetone (111), cholesterol (119) |
| `reactive_intermediates/` | 10 | ethyl-radical (106), propyl-radical (106), butyl-radical (106), cyclobutyl-cation (116), cyclopropylmethyl-cation (116), cyclobutyldiazonium (116), cyclopropylmethyldiazonium (116), 7-7-dimethylnorborn-1-yl-cation (125), 7-7-dimethylnorborn-1-yl-radical (125), 7-7-dimethylnorborn-1-yl-anion (125) |
| **Total** | **108** | |

Checks:

- PubChem named each uncharged, closed-shell structure as intended, with the same formula
  that umol resolves. It treats radicals and cations as their parent compounds, so those
  hydrogen counts were checked directly.
- `canonical_eq` found no duplicate among 1,432 existing inputs after hydrindane was dropped
  as `hydrocarbons/bicyclo-4.3.0-nonane`.
- `:stereo-atoms` input matches the `#T` input of `stereo_tetrahedral/l-alanine`.
- Inputs were lowered from SMILES with the fixture defaults, so `#h0` on bracket atoms
  survives. Resolved outputs were unchanged by that correction.
- The nine inositols follow [IUPAC, Nomenclature of Cyclitols, Eur. J. Biochem. 57 (1975), 1–7](<../materials/valence/IUPAC-1975-nomenclature-of-cyclitols.pdf>),
  doi:10.1111/j.1432-1033.1975.tb02270.x, formulas (1)–(9) and Recommendation I-10. Inputs
  rebuilt from those fractions and numbering directions equal the fixtures under
  `canonical_eq`, and the nine are pairwise distinct with only the chiro pair chiral.
- 1D- and 1L-chiro-inositol share the fraction 1,2,4/3,5,6 and differ in numbering direction.
  The book's Abb. 4.21 uses the older enantiomer-specific prefixes: (+) 1,2,5/3,4,6 and
  (−) 1,2,4/3,5,6. Its (+) equals 1D and its (−) equals 1L.
- NIST WebBook 2D files agree for myo, scyllo, allo, epi, muco and chiro. Its 3D myo file is
  all-cis, contradicting that record's name and descriptor.

Limits:

- The DSL has only tetrahedral and cis/trans stereo. Axial and planar chirality are not
  planned for now, so spiro[3.3]heptane-2,6-dicarboxylic acid and bicyclo[3.3.1]non-1-ene
  carry no stereo.
- The truxillic and truxinic acid structures carry no stereo because the book's scheme gives
  none. Ketamine is the racemate.
- [1.1.1]Propellane uses an ordinary central bond; the book only reports a proposed 3c2e bond.
- The 7,7-dimethylnorborn-1-yl cation is included although the book presents it as not
  forming.
- Bissecododecahedradiene and propellane polymers are omitted: the first needs its structure
  confirmed, the second is not a single compound.
- Held for `organometallics/`: methylsodium, diethylzinc, methyllithium, tert-butyllithium,
  cyclobutylmagnesium bromide, 4-bornanelithium and carboplatin.
- Book misprints: "Padogan" for pagodane, and C7H14 for heptane.

### Chapter 5: unsaturated aliphatic hydrocarbons

Proposed, awaiting review; the suite is unchanged. Pages 129–174. All 237 proposed inputs
resolve under MostSaturated in AtomTyping and Counts, so the chapter requires no new atom
types or ValenceTable entries. Inputs are in `scratch/bw-ch05/proposed/`.

| Category | Count | Compounds (printed page) |
|---|---:|---|
| `hydrocarbons/` | 29 | 2-methylbut-1-ene (129), 3-methylbut-1-ene (129), pent-1-ene (129), 2-methylpent-1-ene (131), dec-1-ene (131), hex-1-ene (131), oct-1-ene (131), 2-3-dimethylbut-2-ene (135), isododecane (141), 5-ethylidenenorbornene (143), isooctane (147), cyclooctane (149), 2-3-dimethylbuta-1-3-diene (150), isoprene (150), penta-1-3-diene (154), 3-methylhexa-1-5-diene (156), hepta-1-5-diene (156), hexa-1-5-diene (156), dimethylfulvene (158), fulvalene (158), 6-methylfulvene (159), bicyclo-3.2.0-hepta-2-6-diene (159), cycloheptatriene (159), norcaradiene (159), cyclooctatetraene (160), cyclododecane (162), fulvene (167), isohexane (167), 18-annulene (171) |
| `functional_groups/` | 71 | 1-2-dibromoethane (130), isopropyl-hydrogen-sulfate (132), 1-2-dichloropropane (133), 1-chloropropan-2-ol (133), acrylonitrile (133), allyl-alcohol (133), ethanolamine (133), ethylenediamine (133), propylene-glycol (133), vinyl-chloride (133), vinylidene-dichloride (133), 1-2-dibromopropane (135), 1-bromopropane (135), 2-3-dibromo-2-3-dimethylbutane (135), 2-iodopropane (135), allyl-bromide (135), diacetyl-peroxide (136), ethyl-hydrogen-sulfate (136), 1-chloro-2-nitrosoethane (137), 1-chloro-2-nitrosoethane-dimer (137), peracetic-acid (137), triethylamine (139), triethylammonium (139), aibn (140), cumene-hydroperoxide (140), thiophenol (141), dibenzoyl-peroxide (142), tetrafluoroethene (144), 1-1-2-trichloroethane (145), isobutyl-vinyl-ether (145), vinyl-alcohol (146), 1-2-dibromocyclopropane (147), cyclopropylamine (147), cyclopropyltrimethylammonium-hydroxide (147), isobutyl-alcohol (147), cyclobutyltrimethylammonium-hydroxide (148), cyclopentanol (148), performic-acid (148), cyclopentyl-bromide (149), 3-4-dibromobut-1-ene (152), 3-4-dichlorobut-1-ene (152), chloroprene (152), 1-chlorobut-2-ene (153), 3-chlorobut-1-ene (153), 3-methylbutane-1-3-diol (157), acetonitrile (157), butane-1-3-diol (157), bicyclo-2.2.2-oct-5-ene-2-2-3-3-tetracarbonitrile (159), bicyclo-3.2.0-hept-2-en-6-ol (159), bicyclo-3.2.0-hept-2-en-6-yl-methanesulfonate (159), tetracyanoethylene (159), tmeda (160), 7-8-dibromobicyclo-4.2.0-octa-2-4-diene (161), 1-4-dibromobut-2-yne (164), 2-3-dichloroprop-1-ene (164), tetrabromomethane (166), 6-6-dibromobicyclo-3.1.0-hex-2-ene (167), 1-1-2-2-tetrachloroethane (170), pentan-1-ol (170), tetrachloroethene (170), trichloroethene (170), 1-1-diiodoethane (171), vinyl-iodide (171), but-3-yn-2-ol (172), ethyl-vinyl-ether (172), hex-3-yne-2-5-diol (172), adiponitrile (173), but-2-yne-1-4-diol (173), butane-1-4-diol (173), hexamethylenediamine (173), propargyl-alcohol (173) |
| `carbonyl_derivatives/` | 32 | acetone-cyanohydrin (133), acrolein (133), acrylic-acid (133), butyraldehyde (133), diacetone-alcohol (133), isobutyl-methyl-ketone (133), mesityl-oxide (133), methyl-methacrylate (133), vinyl-acetate (133), acetaldol (157), n-methyl-2-pyrrolidone (157), cyclopentadienecarboxylic-acid (158), divinyl-ketone (158), benzyl-methyl-carbonate (159), bicyclo-3.2.0-hept-2-en-6-one (159), ketene (159), monomethyl-carbonate (159), chloro-ethyl-ketene (160), chloro-isopropyl-ketene (160), chloro-methyl-ketene (160), cyclododecanone (162), cyclododecanone-oxime (162), laurolactam (162), carbon-suboxide (163), dihydrophthalic-acid (167), dihydrophthalic-anhydride (167), bicyclo-2.2.0-hex-5-ene-2-3-dicarboxylic-acid (168), bicyclo-2.2.0-hex-5-ene-2-3-dicarboxylic-anhydride (168), vinylacetic-acid (168), methyl-9-10-dibromostearate (170), vinyl-stearate (172), methyl-vinyl-ketone (173) |
| `heterocycles/` | 7 | epichlorohydrin (133), propylene-oxide (133), 1-2-4-trioxolane (139), tetrathiafulvalene (158), benzvalene-derived-azo-compound (167), benzvalene-ptad-adduct (167), ptad (167) |
| `aromatic/` | 26 | cumene (133), ethylbenzene (133), phenol (133), styrene (133), m-chlorobenzoic-acid (137), m-chloroperbenzoic-acid (137), hydroquinone (141), pyrogallol (141), methyl-benzoate (149), benzaldehyde (158), benzophenone (158), hexaphenylfulvalene (158), phenylfulvene (158), phenylacetaldehyde (161), terephthalaldehyde (161), sa-1-3-di-naphthalen-1-yl-1-3-diphenylpropa-1-2-diene (164), diphenylbutadiyne (165), 1-1-6-6-tetraphenylhexa-2-4-diyne-1-6-diol (166), 1-6-diiodo-1-1-6-6-tetraphenylhexa-2-4-diyne (166), diaryl-penta-diynol (166), tetraaryl-9-cumulene (166), tetraaryl-octatetrayne-diol (166), tetraphenylbutatriene (166), tetraphenylhexapentaene (166), 7-phenylbicyclo-4.2.0-octa-1-3-diene (167), phthalic-acid (167) |
| `stereo_cis_trans/` | 16 | e-but-2-ene (129), e-pent-2-ene (129), z-pent-2-ene (129), e-cycloheptene (149), z-cycloheptene (149), e-1-4-dibromobut-2-ene (152), e-1-4-dichlorobut-2-ene (152), z-1-4-dibromobut-2-ene (152), e-z-hexa-2-4-diene (153), z-z-cycloocta-1-5-diene (160), e-e-e-cyclododeca-1-5-9-triene (162), e-e-z-cyclododeca-1-5-9-triene (162), 1-4-bis-4-nitrobiphenyl-2-2-diyl-butatriene (165), 1-4-dibromo-1-4-diphenylbutatriene (165), 2-3-dichloro-1-1-4-4-tetraphenylbut-2-ene (166), e-1-2-dichloroethene (170) |
| `stereo_tetrahedral/` | 14 | cis-9-oxabicyclo-6.1.0-nonane (137), trans-9-oxabicyclo-6.1.0-nonane (137), cis-cyclopentane-1-2-diol (149), r-r-trans-cyclopentane-1-2-diol (149), s-s-trans-cyclopentane-1-2-diol (149), s-1-phenylethylamine (150), cis-3-4-dimethylcyclobutene (154), trans-3-4-dimethylcyclobutene (154), endo-6-chloro-6-methylbicyclo-3.2.0-hept-2-en-7-one (160), exo-6-chloro-6-methylbicyclo-3.2.0-hept-2-en-7-one (160), bicyclo-4.2.0-octa-2-4-7-triene (162), 7-phenylbicyclo-4.2.0-octa-1-4-diene (167), r-rasagiline (174), r-selegiline (174) |
| `metabolites/` | 12 | 1-aminocyclopropane-1-carboxylic-acid (146), dihydrosterculic-acid (148), lactobacillic-acid (148), sterculic-acid (148), calcitriol (162), cholecalciferol (162), previtamin-d3 (163), methyl-oleate (169), oleic-acid (169), stearolic-acid (170), stearic-acid (172), dopamine (174) |
| `natural_products/` | 5 | cocaine (159), pseudopelletierine (160), 7-dehydrocholesterol (162), nystatin-a1 (163), ethinylestradiol (174) |
| `reactive_intermediates/` | 25 | propyl-cation (134), isopropyl-cation (135), tetramethylethylene-bromonium (135), 1-bromopropan-2-yl-radical (136), acetoxy-radical (136), carbonyl-oxide (139), phenylthiyl-radical (141), 3-amino-1-cyanoethyl-anion (142), ethyl-cation (146), ethyloxonium (146), cyclopropenylidene (148), cyclopropenylium (148), allyl-cation (151), 4-bromobut-2-enyl-cation (152), methylallyl-cation (153), pentadienyl-radical (155), allyl-radical (156), 3-hydroxy-1-1-dimethylpropyl-cation (157), hydroxymethyl-cation (157), cyclooctatetraene-dianion (160), 4-bromo-1-4-diphenylbutatrienyl-cation (165), isobenzene (166), dibromocarbene (167), acetylide (168), vinyl-cation (170) |
| **Total** | **237** | |

Of 423 extracted entries, 81 duplicate existing inputs under `canonical_eq` and 67 are simple
inorganic reagents left to the HW pass. 30 are held for `organometallics/` and
`organoelement/`. Tropone and tropolone are omitted as Kekulé forms of `aromatic/` inputs;
the amide anion, hydrogen sulfate, methyl radical and methylene as existing or inorganic
species; atropine until chapter 36, because its C3 configuration is undefined in PubChem.

Checks:

- Every compound name found in the page text is on its listed page. PubChem names were
  compared with the intended structures; it treats radicals and most cations as parent
  compounds and misreads the cyclopropenylium SMILES, so those were checked directly.
- PubChem has no record for 14 inputs, so these rest on the book's drawings alone.

Limits and items to review:

- Resolution rejects `#C` on a cumulene's central bond, so the book's E/Z cumulene pairs are
  one stereo-free input each. Axial and planar chirality are likewise unencoded.
- The [9]cumulene family uses 2,6-di-tert-butylphenyl as drawn on p. 166. The literature
  compounds are believed to carry 3,5-di-tert-butylphenyl; this is not checked against a
  primary source.
- The chloroketene adducts follow the drawn 7-one on p. 160, opposite to the book's own
  neighbouring scheme and to the known 6-one regiochemistry.
- isooctane, isododecane and isohexane are assumed single isomers; the book gives no
  structure.
- Book misprints: Tetrathiofulvalen, Tetraethylmethylendiamin for TMEDA, Mesitylchlorid for
  mesityl oxide, Latobacillinsäure, and cis-/trans-1,2-Cyclobuten.

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

For non-binary molecules, a new default atom type serving only a small set or an
individual molecule has a high inclusion threshold. Literature support establishes
the description but does not by itself justify default inclusion. Check reuse of
existing types and require a strong use case before proposing another entry.
Binary compounds remain broadly in scope. This applies to all bonding categories,
in addition to the multicenter restriction below.

Within this scope, source-supported ordinary localized-valence states and ions
with ligand donations are accepted under the standing instruction. Multicenter
descriptions, complicated ternary compounds and capacity changes require explicit
user review.

Do not add speculative occupation states prospectively. A justified extension of
an established assignment by analogy is acceptable, with restraint in adding new
rows on analogy alone. Suggested electronic structures that use existing or
approved types may be covered by conformance fixtures; identify adopted analogies
as such rather than presenting them as established source assignments.

Default multicenter entries are limited to common patterns or selected molecules
for which no simpler description, multicenter or otherwise, is available. Delocalized
molecular orbitals or a possible electron allocation alone do not justify additional
atom types. Hydrogen-bond, donor–acceptor and aromatic descriptions retain their own categories.
Use the smallest adequate multicenter units, preserving localized and dative bonds
where they provide an adequate description. Whole-cage descriptions remain allowed
where no adequate smaller model is available. Reuse common per-atom contributions
across cages rather than adding scaffold-specific types.

For descriptions requiring many multicenter bonds, require a common compound
class before adding default states: the aggregate #m value does not constrain
the number, kind or arrangement of those bonds.

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
| [Implicit H2](../umol-graph/tests/resolution/data/hydrogen/h2-implicit.edn) | H05 |
| [H−](../umol-graph/tests/resolution/data/atomic_ions/h-.edn) | H04 |
| [B2H6](../umol-graph/tests/resolution/data/multicenter_bonds/b2h6.edn) | H06 and B#v2#m; explicit terminal H, two three-center bridges |
| [[RhH(NH3)5]2+](../umol-graph/tests/resolution/data/co_group/pentaamminehydridorhodium-2+.edn) | H07 and Rh#c+3#n3#t6; explicit N–H bonds reuse N#v3#d |
| [HF2−](../umol-graph/tests/resolution/data/multicenter_bonds/hf2-.edn) | H08 and F#c-#n3#m2; three-center/four-electron [2,0,2] |
| [H3+](../umol-graph/tests/resolution/data/multicenter_bonds/h3+.edn) | H06 and H08; three-center/two-electron [1,0,1] |
| [[HCr2(CO)10]−](../umol-graph/tests/resolution/data/multicenter_bonds/hcr2-co10-.edn) | H09, Cr#n3#t5#m0 and C#c-#v3#d; three-center/two-electron [0,2,0] |

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
[conformance input](../umol-graph/tests/resolution/data/noble_gases/xeo64-.edn)
uses two Xe=O and four Xe–O− bonds, with neutral Xe. This settles the final
noble-gas example; the atom state already has independent support from XeO4.

F#c-#n2#d2 is added for the Xe2F11+ bridge, donating one pair to each
Xe#c+#n#v#t#m4 center. Its lone-pair count is left for resolution in the
[conformance input](../umol-graph/tests/resolution/data/multicenter_bonds/xe2f11+.edn).

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
[XeO3](../umol-graph/tests/resolution/data/noble_gases/xeo3.edn),
[XeOF4](../umol-graph/tests/resolution/data/multicenter_bonds/xeof4.edn),
[XeF5−](../umol-graph/tests/resolution/data/multicenter_bonds/xef5-.edn),
[XeOF5−](../umol-graph/tests/resolution/data/multicenter_bonds/xeof5-.edn),
[XeF8²−](../umol-graph/tests/resolution/data/multicenter_bonds/xef8-2.edn) and
[XeF7−](../umol-graph/tests/resolution/data/multicenter_bonds/xef7-.edn).
The existing [XeOF3−](../umol-graph/tests/resolution/data/multicenter_bonds/xeof3-.edn)
and [XeO2F3−](../umol-graph/tests/resolution/data/multicenter_bonds/xeo2f3-.edn) inputs
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
| Ch15 | S#u#v#m4; existing E#n2#v2 (E = O, S, Se, Te) | Present; 1 added | SF5: axial single bond, two basal 3c4e units [1,2,1], one unpaired S electron and no lone pair (HW 634; C93R, pp. 1133–1134; contributions: PW). SO/SeO/TeO use a double bond with #s3 (HW 646, 714–715; HH79). Four fixtures; both models agree. Existing EH radical citations completed in Ch20 |
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
| Ch20, remaining existing entries | Existing atomic, aromatic and remaining hydride/radical states | Present; audit complete | Atomic configurations: HW 328; OH, H3S+, H3Se+: HW 584, 625, 705; SH/SeH/TeH: HH79; heterocycles: I98/Q15. Twenty-two fixtures complete coverage of all 100 retained chalcogen states under both AtomTyping policies. Furan has aromatic and Kekulé fixtures. No further registry states or scope removals |
| Ch21, thioether donor | S#n#v2#d | Implemented: 1 | SnCl4(SMe2)2 (R05, pp. m1852–m1853); two S–C bonds, one donated pair and one remaining lone pair. Covered by the Tt27 adduct fixture; no table changes |

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

HH79: Huber and Herzberg, *Molecular Spectra and Molecular Structure IV: Constants
of Diatomic Molecules* (1979),
doi:[10.1007/978-1-4757-0961-2](https://doi.org/10.1007/978-1-4757-0961-2).
The entries were checked through their NIST Chemistry WebBook reproduction:
[SeO](https://webbook.nist.gov/cgi/cbook.cgi?ID=C12640890&Mask=1000) and
[TeO](https://webbook.nist.gov/cgi/cbook.cgi?ID=C13451177&Mask=1000).
Both have X³Σ− ground states. The fixtures use the established O2 graph convention:
two lone pairs per atom, double bond, and triplet spin on the bond (PW).
HW 646 explicitly gives SO double-bond character and a triplet ground state.

The same compilation supplies the doublet ground states of
[SH](https://webbook.nist.gov/cgi/inchi?ID=C13940211&Mask=1000),
[SeH](https://webbook.nist.gov/cgi/cbook.cgi?ID=C13940222&Mask=186C) and
[TeH](https://webbook.nist.gov/cgi/inchi?ID=C13940368&Mask=1000), supporting the existing
single-bond radical descriptions. Atomic ns²np⁴ configurations use HW I, p. 328,
Table 9.1; OH and pyramidal H3S+/H3Se+ use HW I, pp. 584, 625 and 705, respectively.

I98: [Moss (1998), IUPAC fused-ring nomenclature](../materials/valence/Moss-1998-fused-ring-appendix-2.html),
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

B82: [Bojes et al. (1982)](../materials/valence/Bojes-1982-S3N-anion.pdf),
doi:[10.1021/ja00382a017](https://doi.org/10.1021/ja00382a017), pp. 4839–4840;
full paper reviewed. Explicit six-electron π system in S3N−;
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
| EF4 and TeF6 | chalcogens/tef6.edn; multicenter_bonds/{sef4,sf4,tef4}.edn |
| SF5 radical | multicenter_bonds/sf5.edn |
| EF5−, SOF3+, tellurite, implicit HTe− | chalcogens/{sof3+,teh-implicit-,teo3-2}.edn; multicenter_bonds/{sef5-,sf5-,tef5-}.edn |
| Higher-coordinate Te anions; TeF3+ and H3Te+ | chalcogens/{tef7-,tef8-2,teof5-,teof6-2,tef3+,teh3+}.edn |
| Se/Te hexahalide dianions | multicenter_bonds/{sebr6-2,tecl6-2,tebr6-2}.edn |
| Se/Te tetrahalide dianions | multicenter_bonds/{sebr4-2,tecl4-2}.edn |
| Aqua / hydroxido donors, explicit and implicit H | co_group/rh-h2o6*.edn; co_group/rh-oh*.edn |
| Binary cubanes | chalcogens/{se4cl16,te4br16,te4i16}.edn |
| Po hydride and ions | chalcogens/h2po*.edn; atomic_ions/{po+2,po+4}.edn |
| Charged diatomics, atom-localized and bond-charge/spin descriptions | chalcogens/{o2+-atom,o2+-bond,o2--atom,o2--bond,s2--atom,s2--bond}.edn |
| Triplet chalcogen monoxides | chalcogens/{so,seo,teo}.edn |
| Ozone and ozonide, multicenter π descriptions | multicenter_bonds/o3-multicenter.edn; multicenter_bonds/o3-.edn |
| Trisulfide radical anion, localized and multicenter descriptions | chalcogens/s3-localized-.edn; multicenter_bonds/s3-multicenter-.edn |
| Zundel ion, 3c2e and 3c4e descriptions with explicit, implicit and mixed terminal H | multicenter_bonds/h5o2-3c2e*.edn; multicenter_bonds/h5o2-3c4e*.edn |
| Hydrogen-bond compatibility using existing atom types | hydrogen_bonds/{h7o3+,h9o4+,h13o6+,h3o2-}.edn |
| S4N4 HW contributor | multicenter_bonds/s4n4.edn |
| Small nitrides, existing states | chalcogens/{sn+,s2n+}.edn; multicenter_bonds/s4n2.edn |
| S2N2, aromatic and singlet-diradical multicenter descriptions | aromatic/s2n2.edn; multicenter_bonds/s2n2-diradical.edn |
| Aromatic nitride ions | aromatic/{s2n3+,s3n2-2+,se3n2-2+,s4n3+,s4n4-2+,s5n5+,s3n3-}.edn |
| Radical nitride cations | multicenter_bonds/{s3n2+,se3n2+}.edn |
| Acyclic nitride anions | chalcogens/sn2-2.edn; multicenter_bonds/{s3n-,s4n-}.edn |
| S4N4 / Se4N4, localized descriptions | chalcogens/{s4n4,se4n4}-localized.edn |
| Aromatic chalcogen dications | aromatic/{s4,s6,se4,te4}-2+.edn |
| Cyclic sulfur radical cations | multicenter_bonds/{s5+,s6+,s7+}.edn |
| Localized S8²⁺/Se8²⁺ | chalcogens/{s8,se8}-localized-2+.edn |
| Te6²⁺/Te6⁴⁺ resonance contributors | chalcogens/te6-localized-{2+,4+}.edn |
| Singlet phenyloxenium | functional_groups/phenyloxenium.edn |
| Existing EH radicals, implicit H | chalcogens/{oh,sh,seh}-implicit.edn |
| Existing hydronium, sulfonium and selenonium states, implicit/mixed H | chalcogens/{oh3,sh3,seh3}-{implicit,mixed-h1,mixed-h2}+.edn |
| Existing hydrochalcogenide anions, implicit H | chalcogens/{oh,sh,seh}-implicit-.edn |
| Existing H2S/H2Se/H2Te states, mixed H | chalcogens/{h2s,h2se,h2te}-mixed.edn |
| Existing selenium valences | chalcogens/seo2.edn; chalcogens/sef6.edn |
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

- **F18:** Falvey, “Discrete Existence of Singlet Nitrenium Ions Revisited: Computational Studies of Non-Aryl Nitrenium Ions and Their Rearrangements,” *ACS Omega* **3** (2018), 10418–10432. [DOI](https://doi.org/10.1021/acsomega.8b01038); [full-text XML](../materials/valence/Falvey-2018-nitrenium/Falvey-2018-nitrenium.xml), with figures and SI alongside it. Tables 1–2 distinguish the singlet entries from triplet ground states.
- **G18:** Ganesamoorthy et al., “From Stable Sb- and Bi-Centered Radicals to a Compound with a Ga=Sb Double Bond,” *Nat. Commun.* **9** (2018), 87. [DOI](https://doi.org/10.1038/s41467-017-02581-2); [PDF](../materials/valence/Ganesamoorthy-2018-Sb-radicals.pdf). NBO analysis on p. 4 supports two bonds, one lone pair and one unpaired electron at Sb.
- **H81:** Harrison, “Electronic Structure of the Phosphenium Ions PH2+, HPF+, and PF2+,” *J. Am. Chem. Soc.* **103** (1981), 7406–7413. [DOI](https://doi.org/10.1021/ja00415a002); [PDF](../materials/valence/Harrison-1981-phosphenium.pdf). Singlet states with two bonds and a lone pair, pp. 7410–7412.
- **Z22:** Zhang, Yuan and Wang, “Probing the Electronic Structure and Spectroscopy of Pyrrolyl and Imidazolyl Radicals Using High-Resolution Photoelectron Imaging of Cryogenically Cooled Anions,” *Phys. Chem. Chem. Phys.* **24** (2022), 6505–6514. [DOI](https://doi.org/10.1039/D2CP00189F); [accepted manuscript](../materials/valence/Zhang-2022-pyrrolyl.pdf). Section 3.1.3 assigns the observed B(2A1) pyrrolyl state to σ-electron detachment; the registry's atom-local allocation is PW.

Unavailable review: Fisher and Norman, “The Structures of the Group 15 Element(III)
Halides and Halogenoanions,” *Adv. Inorg. Chem.* **41** (1994), 233–271;
[DOI](https://doi.org/10.1016/S0898-8838%2808%2960173-7), cited by HW I, p. 990.
The user supplied the abstract but also cannot access the full chapter. The abstract
establishes structural families, not a bond partition; it adds no accepted assignment.

### Tetrels

Source: [HW 103, Band I, chapter XV](<../materials/valence/Holleman-Wiberg - Anorganische Chemie - 103. Auflage - Band 1.pdf#page=1032>),
printed pp. 994–1206; chapter references pp. 1206–1215. Initial extraction:
2026-11-02. Baseline: C 36, Si 11, Ge 1, Sn 1, Pb 3 registry states.
Tt02–Tt05 are implemented: 12 C states and one Al companion, with 13 fixtures.
The 28 silicon additions are implemented with 30 fixtures and Al's zero-π Counts
entry. The 93 Ge/Sn/Pb additions are implemented with 98 fixtures, covering Zintl
states, tetrelylene acceptors, tetrelyl donors, neutral donor adducts,
aromatic/localized rings, ordinary neutral/charged states, halide bridges and
dibridged hydrides. The adduct fixtures also cover the Ch21 thioether donor companion.
Tt28–Tt31 complete the selected additions: Si 1, Ge 1, Sn 5 and Pb 5, with the
approved Sn/Pb capacity increase to 14. No selected tetrel assignment remains open.

#### Tetrel work table

**Implemented** identifies selected additions now in the repository.
**Proposed** identifies candidate additions for discussion; **Present** identifies
existing types. **Open** identifies missing evidence or an unsettled assignment;
**No addition proposed** records an extracted family outside the proposed expansion.
States use the concrete, zero-valence registry defaults. Each listed E denotes a
separate entry; explicit, mixed and implicit H descriptions require separate states.
HW pages refer to Band I unless stated otherwise. **PW** marks our graph allocation,
distinct from the literature's structural or bonding assignment. Examples guide
coverage of selected states; they are not a requirement to fixture every compound.

| ID | Atom descriptions / gap | Disposition | Examples and evidence |
| --- | --- | --- | --- |
| Tt01 | Existing 36 C states | Audited: retained | HW citations filled throughout. CO removed from the C+ example; C− donor comment broadened. G12 experimentally identifies CH's X²Π ground state and excited a⁴Σ− state. B04 supplies the chemical counting scheme for exocyclic-carbonyl C#v4#a0 (umol allocation: PW); existing 2-pyridone/tropone fixtures retained |
| Tt02 | C#c-4#n4; C#c-2#n2#v2 | Implemented: 2 | Methanides Al4C3/Be2C and allenides Li4C3/Mg2C3; HW 1021–1022. Terminal C²−=C=C²− allocation. Fixtures: atomic C⁴− and C3⁴−. Acetylide C#c-#n#v3 already exists |
| Tt03 | C#v2#d; C#v#h#d; C#h2#d | Implemented: 3 | Neutral carbene donors, including NHCs; HW 1044–1045, II 2165–2167. Three OsCl(NO)(PPh3)2(CH2) fixtures cover explicit/mixed/implicit carbene H, using Fig. 32.64g's donor description. Os occupation supplied; NO+ allocation is PW. No additional aromatic donor type selected |
| Tt04 | C#c-#v2#h#d; C#c-#v#h2#d; C#c-#h3#d | Implemented: 3 | Organyl donors, HW 1060–1062; selected anionic-donor convention (PW). Four ZnMe2 fixtures cover all H representations, including existing C#c-#v3#d; Zn(II) occupation is inferred with Zg03 |
| Tt05 | C#v3#m; C#v2#h#m; C#v#h2#m; C#h3#m; companion Al#v2#m | Implemented: 5 | Al2Me6: HW 1024 explicitly assigns 3c2e Al–C–Al bridges. Neutral contributions [1,1,0] and [0,1,1] follow the existing diborane allocation (PW). Four fixtures cover every bridge-C H representation and the Al companion |
| Tt06 | Existing Si states and Ge/Sn/Pb atomic states | Audited: retained | SiH4/SiX4: HW 1082–1097; partial-H derivatives: 1086; SiH3 radical and singlet/triplet SiH2: 1090; SiX2: 1093; SiX3−: 1087, 1143; SiX3+: 1144. Atomic Si's triplet ground state: 1066. Ge/Sn/Pb s²p² ground configurations: 328, Table 9.1. Corrected Si#v3#h's example to SiX3H |
| Tt07 | Si#n#v#h; Si#n#h2; Si#u2#v#h; Si#u2#h2; Si#u#v2#h; Si#u#v#h2; Si#u#h3 | Implemented: 7 | Complete singlet/triplet SiH2 and SiH3-radical H representations; HW 1076–1078, 1090; WA03 observes SiH3. Explicit-H Si#n#v2, Si#u2#v2 and Si#u#v3 already exist |
| Tt08 | Si#c-#n#v2#h; Si#c-#n#v#h2; Si#c-#n#h3 | Implemented: 3 | SiH3−, HW 1087; WA03 Table 1. Explicit-H Si#c-#n#v3 already exists |
| Tt09 | E#c-2#n2#v2; E#c-3#n3#v; E#c-4#n4 (E = Si, Ge, Sn, Pb) | Implemented: 12 | Zintl chains, E2⁶− dumbbells and isolated E⁴− ions; HW 1069–1072, 1165. Ordinary (8−N) descriptions. Ge/Sn/Pb: nine fixtures cover isolated ions, dumbbells and finite chain counting fragments. E⁻ in E4⁴− tetrahedra uses existing Si#c-#n#v3 or Tt20 |
| Tt10 | Si#c-#v3#m2; Si#c-2#m6 | Implemented: 2 | SiF5−: three equatorial ordinary bonds and one axial 3c4e unit; SiF6²−: three 3c4e units. HW 1082 explicitly assigns these partitions; structures 1092–1095. F–Si–F contributions [1,2,1] are PW; F#n3#m exists. No additional fully localized alternatives proposed |
| Tt11 | Si#v4#t; Si#v4#t2 | Implemented: 2 | SiF4·NH3 and SiF4·2NH3, HW 1094; SiCl4 donor adducts, 1096. These are neutral donor adducts, distinct from the bare fluoroanions in Tt10 |
| Tt12 | Si#c+#v3#t; Si#c+#v3#t2 | Implemented: 2 | Donor-stabilized silylium ions, including Me3Si(py)+ and Ph3Si(dipy)+; HW 1144, Fig. 15.89 |
| Tt13 | E#n#v2#t (E = Si, Ge, Sn, Pb) | Implemented: 4 | Tetrelylene acceptors retain one lone pair: HW 1081, 1151–1152, 1173, 1200. Fixtures: IPr→SiCl2 (G09), IMes→GeI2 (A93; neutral donor encoding: PW), SnCl2(NMe3) (HW 1173), cyclic disilylated plumbylene–PEt3 (A12, compound 2). No higher acceptor counts added |
| Tt14 | E#c-2#n#v2#a2 (E = Si, Ge, Sn) | Implemented: 3 | Silole/germole/stannole dianions, HW 1157, 1204–1205; atom-local charge allocation is PW. Aromatic Si/Ge E5⁶− and E6¹⁰− descriptions (HW 1071; vS96) are excluded from the default scope; E#c-#n#v2#a is not retained. Localized dianions use Tt09 |
| Tt15 | Si#v3#a; Ge#v3#a | Implemented: 2 | Substituted silabenzene/silanaphthalene and germabenzene/germanaphthalene; HW 1157, 1204–1205. Aromatic and localized notations; localized types exist for both Si and Ge (Tt18) |
| Tt16 | Si#c+#v3#a0; Ge#c+#v3#a0 | Implemented: 2 | Aromatic 2π E3R3+ rings, HW 1154–1155, 1203. One E+ contributes zero π electrons; two neutral E atoms use Tt15. Atom-local charge allocation is PW. Localized Ge3R3+ uses the implemented Ge#v4 (Tt18) and Ge#c+#v3 (Tt23) |
| Tt17 | Si#c-#v3#a2 | Implemented: 1 | 6π Si4R4²− ring, HW 1154; weakly aromatic silole monoanion, 1157. Localized Si#c-#n#v3 already exists; π/charge allocation is PW |
| Tt18 | E#v4; E#v3#h; E#v2#h2; E#v#h3; E#h4 (E = Ge, Sn, Pb) | Implemented: 15 | EH4, higher hydrides, EX4, ER4 and ordinary multiple-bond descriptions; HW 1167–1169, 1170–1176, 1191–1193, 1202–1206. WA03 establishes PbH4 experimentally. Fifteen EH4 fixtures cover all H representations. Counts has missing-H limitations for Sn/Pb, below |
| Tt19 | E#n#v2; E#n#v#h; E#n#h2 (E = Ge, Sn, Pb) | Implemented: 9 | Singlet EH2, molecular EX2, EO and ER2; HW 1076–1077, 1168–1169, 1199; WA03 establishes EH2, including singlet PbH2 (Table 2). Nine EH2 fixtures cover all H representations |
| Tt20 | E#c-#n#v3 (E = Ge, Sn, Pb); E#c-#n#v2#h; E#c-#n#v#h2; E#c-#n#h3 (last three: E = Ge, Pb) | Implemented: 9 | R3E−, EX3− and E4⁴− tetrahedra: HW 1194, 1173–1174, 1165. GeH3−: HW 1168; PbH3−: WA03 Table 1 and p. 6585. Eight GeH3−/PbH3− fixtures cover all H representations; SnCl3− covers Sn. No SnH3− inference from WA03 |
| Tt21 | E#u#v3; E#u#v2#h; E#u#v#h2; E#u#h3 (E = Ge, Sn, Pb) | Implemented: 12 | ER3 radicals, HW 1193, 1195; EH3 radicals, WA03. Twelve EH3 fixtures cover all H representations. Counts has missing-H limitations for Sn/Pb, below |
| Tt22 | E#n#u#v; E#n#u#h (E = Si, Ge, Sn, Pb) | Implemented: 8 | EH doublets, HW 1077; WA03 discusses the observed monohydrides and gives PbH's doublet state in Table 2. Six Ge/Sn/Pb fixtures join the two Si fixtures, covering explicit and implicit H |
| Tt23 | E#c+#v3 (E = Ge, Sn, Pb) | Implemented: 3 | R3E+ cores, HW 1195; uncoordinated Ge/Sn examples and donor-stabilized species. Me3Sn+/Me3Pb+ fixtures use the solution cation cores (NBS74). Ge is covered by localized Ge3R3+; Si counterpart already exists |
| Tt24 | Sn#c+#n#v; Pb#c+#n#v | Implemented: 2 | SnOH+/PbOH+, HW 1181, 1185; RPb+, 1201. Divalent cations retain one lone pair; two hydrolysis-species fixtures omit solvent |
| Tt25 | Sn#c+2#n | Implemented: 1 | Sn(II) ions and salts, HW 1181–1182; atomic-ion fixture added. Existing Pb#c+2#n and its four-mode fixture retained; HW 1185 citation added |
| Tt26 | E#c-#v3#d (E = Si, Ge, Sn, Pb) | Implemented: 4 | Terminal tetrelyl donors, HW 1143, 1194; selected anionic-donor convention (PW). Ge/Sn/Pb fixtures: Pd(PEt3)2(GePh3)2 (B67), [Ni(CO)3(SnPh3)]− (HW 1194), trans-Pt(PEt3)2(PbPh3)2 (D68). Metal occupations supplied. No implicit-H donor variants |
| Tt27 | E#v4#t; E#v4#t2 (E = Ge, Sn, Pb) | Implemented: 6 | Tetrel(IV) donor adducts; HW 1170, 1175, 1195. One donor: GeCl4(NMe3) (BW72), Ph3SnCl(OPPh3) and Ph3PbBr(OPPh3) (E92). Two donors: GeF4(PMe3)2 (K21), SnCl4(SMe2)2 (R05), PbPh2Cl2(bipy) (Ca01). Six fixtures cover these states and companion S#n#v2#d (Ch21) |
| Tt28 | Sn#c+2#v2#t4; Pb#c+2#v2#t4 | Implemented: 2 | Solid SnF4/PbF4 sheets: two terminal F and four bridging F per E, HW 1170, Fig. 15.109e. Ordinary terminal bonds plus four F−→E donations (PW), using existing F#c-#n2#d2; two sheet-fragment fixtures |
| Tt29 | Sn#c+#v3#t2; Pb#c+#v3#t2 | Implemented: 2 | R3EX chains, including Me3SnF/Ph3PbCl; HW 1192–1193, Fig. 15.118a. Two bridging X− donors per R3E+ center (PW), using existing halide donor states; two chain-fragment fixtures |
| Tt30 | EX5−/EX6²− (E = Ge, Sn, Pb); Sn2X5− | Present; donor graphs accepted | EX4 plus one/two X− donors uses existing E#v4#t/#t2; HW 1175–1176 explicitly describes halide addition, H05 treats GeCl5−/SnCl5− as Lewis acid–base adducts. Sn2F5−/Sn2Cl5− uses existing Sn#n#v2#t and bridging X− donors (HW 1174). Integer partitions are PW; seven trial graphs agree across both models/policies. KB99 also gives a three-center description of SnX6²−; no extra multicenter row is selected where the donor graph suffices |
| Tt30a | Sn#c+2#n#t4 | Implemented: 1 | SnCl4²−: Sn(II) ns² core with four chloride donors (PW); HW 1174, Fig. 15.110h. One fixture agrees across both models/policies |
| Tt30b | Sn#c+2#n#t6; Pb#c+2#n#t6 | Implemented: 2; capacity 14 | SnBr6⁴− and PbCl6⁴−/PbBr6⁴−: E(II) ns² core with six halide donors (PW); HW 1174; S19 confirms the retained ns² pair in Cs4SnBr6/Cs4PbBr6. Both bromide fixtures agree across both models/policies and pass the valence invariants with the approved Sn/Pb capacity increase to 14 |
| Tt31 | E#n#v#m (E = Si, Ge, Sn, Pb); Pb#n#h#m | Implemented: 5 | E2H2: one E–E bond, two 3c2e H bridges and one lone pair per E; HW 1090 explicitly gives this Si2H2 description, WA03 establishes the hydrides experimentally, LKF05 p. 6294 corroborates the E–E σ bond and E–H donation. Pb2H4 uses two terminal Pb–H bonds and two 3c2e bridges, with no direct Pb–Pb bond (T90 p. 2136); the same #v state covers explicit H, and #h covers implicit H. Neutral contributions [1,1,0]/[0,1,1] follow diborane; no new H state. Seven fixtures agree across both models/policies |
| Tt32 | Higher acceptor counts and seven-/eight-coordinate complexes | No addition proposed | HW 1081, 1152, 1167. Extracted examples do not yet justify additional default states beyond the common adduct families above |
| Tt33 | Existing Pb#c+4 citation | Audited: retained | HW 1187 explicitly uses Pb⁴⁺ in the schematic acid dissolution of PbO2 and discusses Pb(IV) acetate/sulfate. The existing state is retained for this ionic salt description, with its citation filled; no new Sn⁴⁺ state proposed. The C#a0 citation is resolved in Tt01 |
| Tt34 | Metal–carbene/carbyne/carbido and haptic ligands | Present / no further addition proposed | Ordinary M=C/M≡C descriptions reuse C#v4; carbene donation is Tt03. HW II 2165–2167, 2175 ff. Haptic multidonor accounting remains separate work; no new carbon families proposed merely to recast every metal–ligand interaction |
| Tt35 | C2S2/C2OS and other binary carbon oxides/sulfides | Present; triplet graphs checked | HW 1030–1042, including linear triplet S=C=C=S and O=C=C=S (1042). Assigning #s3 to the central double bond (PW) resolves identically under both models/policies and satisfies the valence invariants, using existing C#v4 and O/S#n2#v2. No new type required |
| Tt36 | Larger Zintl/Wade clusters, metal-cluster carbon and fulleride families | No addition proposed | HW 1014–1023, 1071, 1165–1166, 1197. No cluster-size-specific multicenter entries. Ordinary tetrahedra are covered above; the aromatic Zintl rings are excluded in Tt14 |
| Tt37 | Further excited-state and gas-ionized families | No addition proposed | HW 1077 includes calculated excited hydride states. Do not complete every spin/ionization ladder from those listings |

The retained extraction lists 134 concrete candidates: C 12, Si 28, Ge 30, Sn 31,
Pb 32 and one Al companion. All 134 states are implemented and covered by 141
fixtures. Tt27 also adds the Ch21 sulfur companion, outside this count.
No addition duplicates an earlier entry.

Tt28–Tt31 implementation (2026-11-02): fourteen fixtures cover all twelve additions.
All four resolver/policy combinations agree; all 56 determined outputs satisfy
the valence invariants and need no tie-break. All 1168 earlier snapshots are
unchanged. Sheet/chain fixtures are finite fragments cut at boundary halides;
all tetrel coordination spheres are complete. The seven reviewed Tt30 graphs
use existing types and remain trial-only; no additional fixture is required for them.

Tt27 implementation (2026-11-02): six primary-source fixtures cover all six
acceptor additions and the Ch21 thioether companion. Each graph resolves identically
across both resolvers/policies; all 24 outputs need no tie-break and satisfy the
valence invariants. The trial left admission, outputs and reports unchanged for all
1162 earlier fixtures; their implemented snapshots are also unchanged.
No table or capacity changes.

Tt26 implementation (2026-11-02): all three Ge/Sn/Pb donor additions have fixture
coverage. Each fixture agrees across both resolvers/policies, requires no tie-break
and satisfies the valence invariants. The trial left admission, outputs and reports
unchanged for all 1159 earlier fixtures; their snapshots are also unchanged.
No table or capacity changes.

Tt14–Tt16 implementation (2026-11-02): four aromatic types plus Ge#v4
(Tt18) and Ge#c+#v3 (Tt23), with Counts Ge `[1]` and shifted Ga `[0]` aromatic
valences. Paired aromatic/localized graphs for tetraphenylgermole and -stannole
dianions, Tbt-germabenzene and Ge3(SitBu3)3+ cover all six types. All 32 outputs
agree across models/policies and satisfy the valence invariants. The prior trial
left admission, outputs and reports unchanged for all 1095 earlier fixtures;
implementation also leaves their snapshots unchanged. HW 1203–1205 explicitly
supplies the aromatic descriptions; charge placement follows the selected Si convention (PW).

Carbon verification (2026-11-02): all 13 added states have fixture coverage. Each of
the 13 fixtures agrees across both models/policies; all 52 outputs satisfy the valence invariants.
The prior trial found unchanged admission/results/reports for all 1041 earlier
fixtures, of which 393 contain C; implementation also leaves their snapshots unchanged.
Neutral C grows from 23 to 30 entries, C− from 7 to 10, with one new entry each for
C²− and C⁴−. The seven neutral additions require donation or multicenter membership;
ordinary graphs supply zero/nonmember readings and exclude them. Runtime cost was
not measured.

Silicon implementation (2026-11-02): all 27 added states have conformance coverage.
All 29 new fixtures agree across both models/policies; all 116 outputs satisfy the
valence invariants. Tetramethylsilole dianion, Si3R3+ and Si4R4²− have aromatic and
localized fixtures.
IPr→SiCl2 covers the silylene acceptor. The prior trial left admission, outputs and
reports unchanged for all 1054 earlier fixtures (nine contain Si); implementation
also leaves their snapshots unchanged. Si grows from 11 to 38 states: neutral
9→22, +1 1→4, −1 1→7, plus three −2, one −3 and one −4 states. Runtime cost was
not measured.

Counts needs `#u2` for triplet SiH2; `#s3` alone does not control its
lone-pair/unpaired-electron split.
For aromatic Si4R4²−, unspecified spin/contributions let Counts and the default
MinElectronCount aromatic policy select a 2π assignment with four unpaired electrons.
Supplying `#u0` on the ring atoms selects the intended 6π assignment under both
models/policies. The triplet and Si4R4²− fixtures supply those accepted inputs.

#### Counts comparison

Tt30b raises Sn/Pb capacity from 12 to 14 for the accepted binary hexahalide-anion
description: six accepted pairs plus the retained ns² pair require 14 electrons.
Counts rejected these graphs at 12 despite consistent electron accounting; the
approved bound admits them generally. The maximum localized valence stays 4,
and target covalences stay [2,4]. The other selected states need no capacity change.
Counts already includes Ge contribution 1 and Al/Ga contribution 0 for the shifted
Si+/Ge+ lookups.
The aromatic dianions use the existing shifted S/Se/Te contribution-2 rows.
The permissive aromaticity model already admits these elements; no change to a
named restrictive model is needed.

With implicit H unspecified, the first Sn/Pb target at or above the supplied ordinary
valence limits the H budget. Their current targets [2,4] exclude Tt18's v2h2/vh3/h4
and Tt21's u v2h/u vh2/u h3 descriptions in that reading. Supplied H counts or explicit
H bonds avoid this restriction. Changing the targets would also affect divalent
defaults and is not proposed without considering that behavior. Element charge/H
enumeration bounds are a separate limitation of invariant enumeration, not a
Counts rejection of supplied values.

#### Primary sources beyond HW

- **G09:** Ghadwal et al., “Lewis Base Stabilized Dichlorosilylene,” *Angew. Chem. Int. Ed.* **48** (2009), 5683–5686. [DOI](https://doi.org/10.1002/anie.200901766); [SI](../materials/valence/Ghadwal-2009-dichlorosilylene-SI.pdf). Supplied SI pp. 2–3 and 7 give the synthesis and crystallographic data for IPr–SiCl2 (compound 1); HW 1151–1152 supplies the donor/acceptor description.
- **A93:** Arduengo, Dias, Calabrese and Davidson, “A Carbene-GeI2 Adduct: Model of the Non-Least-Motion Pathway for Dimerization of Singlet Carbenes,” *Inorg. Chem.* **32** (1993), 1541–1542. [DOI](https://doi.org/10.1021/ic00061a004); [PDF](../materials/valence/Arduengo-1993-GeI2-carbene.pdf). Compound 2 is the characterized IMes–GeI2 adduct; Scheme I and the pyramidal Ge geometry support carbene donation. The fixture uses the neutral donor–acceptor encoding (PW).
- **A12:** Arp et al., “Dispersion Energy Enforced Dimerization of a Cyclic Disilylated Plumbylene,” *J. Am. Chem. Soc.* **134** (2012), 6409–6415. [DOI](https://doi.org/10.1021/ja300654t); [archived full-text XML](../materials/valence/Arp-2012-plumbylene/Arp-2012-plumbylene.xml), with figures and SI alongside it. Compound 2 and Fig. 1 supply the plumbylene–PEt3 fixture; the crystallographic discussion and NBO analysis explicitly identify P→Pb donation. Full text and structure figure reviewed.
- **WA03:** Wang and Andrews, “Infrared Spectra of Group 14 Hydrides in Solid Hydrogen: Experimental Observation of PbH4, Pb2H2, and Pb2H4,” *J. Am. Chem. Soc.* **125** (2003), 6581–6587. [DOI](https://doi.org/10.1021/ja029862l); [PDF](../materials/valence/Wang-Andrews-2003-group14-hydrides.pdf). Matrix-IR evidence for hydrides, including PbH4, PbH3− and bridged Pb2H2/Pb2H4; Table 2 gives the calculated electronic states used in the spectral assignments.
- **NBS74:** McNesby, *Measures for Air Quality (1972–1973): Annual Report—FY 1973*, NBS Technical Note 828 (1974). [PDF](../materials/valence/McNesby-1974-air-quality-TN828.pdf#page=131). Brinckman and Iverson's project report, p. 123, identifies aqueous Me3Sn+/Me3Pb+ chemistry by proton NMR and laser-Raman spectroscopy; supports the two cation-core fixtures.
- **B67:** Brooks, *Germyl and Hydrido Palladium Complexes*, doctoral thesis, Durham University (1967). [Repository](https://etheses.durham.ac.uk/id/eprint/8826/); [PDF](../materials/valence/Brooks-1967-germyl-palladium-thesis.pdf). Pages 106–108 establish isolated Pd(PEt3)2(GePh3)2; the fixture uses the selected anionic-donor allocation (PW).
- **D68:** Deganello, Carturan and Belluco, “Platinum–Lead Bonded Complexes,” *J. Chem. Soc. A* (1968), 2873–2878. [DOI](https://doi.org/10.1039/J19680002873); [PDF](../materials/valence/Deganello-1968-Pt-Pb-complexes.pdf). Pages 2873–2874 and 2877 establish trans-Pt(PEt3)2(PbPh3)2 by synthesis, analysis and IR. The fixture's PbPh3−→Pt²⁺ allocation is PW.
- **BW72:** Bilton and Webster, “Crystal Structure of Germanium Tetrachloride–Trimethylamine,” *J. Chem. Soc., Dalton Trans.* (1972), 722–724. [DOI](https://doi.org/10.1039/DT9720000722); [PDF](../materials/valence/Bilton-1972-GeCl4-NMe3.pdf). Discrete GeCl4(NMe3), with axial NMe3 in a trigonal bipyramid.
- **K21:** King, Levason and Reid, “Neutral and cationic germanium(IV) fluoride complexes with phosphine coordination,” *Dalton Trans.* **50** (2021), 17751–17765. [DOI](https://doi.org/10.1039/D1DT03339E); [PDF](../materials/valence/King-2021-GeF4-phosphine-adducts.pdf). Pages 17752–17753 establish octahedral trans-GeF4(PMe3)2 by synthesis and X-ray diffraction.
- **E92:** Eppley et al., “Structure and equilibria in triorganolead halide adduct formation,” *J. Organomet. Chem.* **431** (1992), 133–142. [DOI](https://doi.org/10.1016/0022-328X(92)80112-B); [PDF](../materials/valence/Eppley-1992-organotin-organolead-adducts.pdf). Pages 140–141 establish five-coordinate Ph3SnCl(OPPh3) and Ph3PbBr(OPPh3), with axial O-bound phosphine oxide.
- **R05:** Ruth, Lerner and Bolte, “cis-trans-Tetrachlorobis(dimethyl sulfide-κS)tin,” *Acta Cryst. E* **61** (2005), m1852–m1853. [DOI](https://doi.org/10.1107/S1600536805026346); [PDF](../materials/valence/Ruth-2005-SnCl4-SMe2-adduct.pdf). Octahedral SnCl4(SMe2)2 with two S-bound donors; also supports Ch21.
- **Ca01:** Casas et al., “The crystal structure of dichlorodiphenylbipyridine lead(IV), [PbPh2Cl2(bipy)],” *Main Group Metal Chemistry* **24** (2001), 297–298. [DOI](https://doi.org/10.1515/MGMC.2001.24.5.297); [PDF](../materials/valence/Casas-2001-PbPh2Cl2-bipyridine.pdf). Discrete distorted-octahedral complex with chelating 2,2′-bipyridine; Fig. 1 supplies the graph.
- **LKF05:** Lein, Krapp and Frenking, “Why Do the Heavy-Atom Analogues of Acetylene E2H2 (E = Si–Pb) Exhibit Unusual Structures?”, *J. Am. Chem. Soc.* **127** (2005), 6290–6299. [DOI](https://doi.org/10.1021/ja042295c); [PDF](../materials/valence/Lein-2005-group14-dihydrides.pdf). Page 6294 and Fig. 3b give one E–E σ bond and two E–H donor interactions for the doubly bridged structure. Corroborates the HW Si2H2 graph; the neutral contribution allocation is PW.
- **T90:** Trinquier, “Double Bonds and Bridged Structures in the Heavier Analogues of Ethylene,” *J. Am. Chem. Soc.* **112** (1990), 2130–2137. [DOI](https://doi.org/10.1021/ja00162a014); [PDF](../materials/valence/Trinquier-1990-group14-tetrahydrides.pdf). Page 2136 explicitly assigns diborane-like 3c2e bridges and omits direct E–E bonding in the bridged heavier E2H4 structures; WA03 supplies experimental Pb2H4 support.
- **H05:** Hao et al., “Effect of Substituents on the Strength of A–Cl− (A = Si, Ge, and Sn) Bonds in Hypervalent Systems: ACl5−, ACl4F−, and A(CH3)3Cl2−,” *J. Phys. Chem. A* **109** (2005), 2026–2034. [DOI](https://doi.org/10.1021/jp040743x); [PDF](../materials/valence/Hao-2005-tetrel-pentachloroanions.pdf). Experimental dissociation energies and calculations use AX4/X− Lewis acid–base fragments; supports that graph partition without an additional multicenter state.
- **KB99:** Kravchenko and Buslaev, “Electronic effects of substitution and reconstruction of the coordination polyhedron in adducts of Main Group IV element halides,” *Russ. Chem. Rev.* **68** (1999), 709–726. [DOI](https://doi.org/10.1070/RC1999v068n09ABEH000497); [PDF](../materials/valence/Kravchenko-1999-tetrel-halide-bonding.pdf). Page 711 discusses three equivalent three-center bonds in SnX6²−; recorded as a literature alternative, with no extra default state proposed.
- **S19:** Shi, Han, Chen and Du, “Impact of metal ns² lone pair on luminescence quantum efficiency in low-dimensional halide perovskites,” *Phys. Rev. Materials* **3** (2019), 034604. [DOI](https://doi.org/10.1103/PhysRevMaterials.3.034604); [accepted manuscript](../materials/valence/Shi-2019-tetrel-halide-lone-pairs.pdf). Section III.A identifies retained Sn(II)/Pb(II) ns² lone pairs in Cs4SnBr6/Cs4PbBr6. Six halide-donor bonds are our graph allocation (PW).
- **vS96:** von Schnering et al., “Hückel Arenes with Ten π Electrons: Cyclic Zintl Anions Si6¹⁰− and Ge6¹⁰−, Isosteric to P6⁴− and As6⁴−,” *Angew. Chem. Int. Ed. Engl.* **35** (1996), 984–986. [DOI](https://doi.org/10.1002/anie.199609841); [PDF](../materials/valence/von-Schnering-1996-Si-Ge-aromatic-rings.pdf). Reviewed; these aromatic Zintl-ring assignments do not justify additional default-registry states and are excluded in Tt14.
- **B04:** Balaban, Oniciu and Katritzky, “Aromaticity as a Cornerstone of Heterocyclic Chemistry,” *Chem. Rev.* **104** (2004), 2777–2812. [DOI](https://doi.org/10.1021/cr0306790); [PDF](../materials/valence/Balaban-2004-heterocyclic-aromaticity.pdf). Page 2783 discusses pyridone tautomerism; pp. 2788–2791 define zero-electron Z-type ring positions and include carbons bearing exocyclic C=O, with tropone in Scheme 27. The C#v4#a0 allocation retaining the localized C=O bond is PW. This replaces the inappropriate RDKit citation without changing the existing entry or fixtures.
- **G12:** Goebbert, “Photoelectron imaging of CH−,” *Chem. Phys. Lett.* **551** (2012), 19–25. [DOI](https://doi.org/10.1016/j.cplett.2012.08.065); [PDF](../materials/valence/Goebbert-2012-CH-photoelectron.pdf). Page 19 gives the X²Π configuration with a paired nonbonding 3σ orbital and one unpaired π electron, and the a⁴Σ− configuration with three unpaired electrons. Page 21 and Fig. 3 identify both states experimentally. These support the existing single-bond CH descriptions C#n#u#v and C#u3#v, including their implicit-H variants.

### Triels

Source: [HW 103, Band I, chapter XVI](<../materials/valence/Holleman-Wiberg - Anorganische Chemie - 103. Auflage - Band 1.pdf#page=1254>),
printed pp. 1216–1421; chapter references pp. 1422–1428. Initial extraction:
2026-11-03. Baseline: B 24, Al 4, Ga 2, In 2, Tl 3 registry states.
The selected additions and Tr09 correction/removal are implemented. Current totals:
B 50, Al 41, Ga 41, In 34 and Tl 25. The four carborane companions bring C to 52;
the aromatic pyridine donor brings N to 53.

For boron, default multicenter types describe reusable per-atom contributions across
local bonds and whole-cage systems. Whole-cage descriptions are accepted, especially
where no adequate smaller model suffices. Cage size, closo/nido/arachno classification
and choice of localization do not each require a separate atom type. HW's local
B–H–B/B–B–B descriptions and its Wade counting both supply representative examples.

#### Triel work table

**Implemented** marks accepted additions; **Not selected** records exclusions from
this pass. States use the concrete, zero-valence registry defaults; each listed E
expands to a separate element. H variants are
listed explicitly. HW references below are printed Band I pages. **PW** identifies
our electron/charge allocation, rather than an assignment stated by HW.
Examples guide coverage of selected states; not every extracted molecule needs a fixture.
Citation keys refer to the [primary sources](#triel-primary-sources) below.

| ID | Atom descriptions / gap | Disposition | Examples and evidence |
| --- | --- | --- | --- |
| Tr01 | Original 35 B/Al/Ga/In/Tl states; 34 after Tr09 | Citation audit complete | The two original borolide entries failed electron conservation; Tr09 corrects one and removes the other. All 34 retained states pass. B ordinary BH3/BX3/BR3 and BH4−/BX4− families: HW 1226, 1248, 1280–1285, 1315–1316. BH/BF and BH2/BO descriptions: 1286, 1287, 1318, 1321. Existing aromatic B: boroxine 1291; borazine 1297, 1303; boratabenzene/borepin 1320–1321. Al/Ga/In/Tl atoms: 1216; Al³⁺: 1335, 1345; Al#v3: 1339, 1345, 1366; Al#v2#m: 1366; heavy E³⁺ and Tl⁺: 1387, 1399. Ce03 supports singlet BH2− and thus B#c-#n#v2 plus its two H variants. S06 independently supports the localized boryl-anion description with a B lone pair; neither source requires replacing these states |
| Tr02 | B#n#h | Implemented: 1 | Implicit-H counterpart of existing B#n#v for BH; HW 1321. BF's retained s² pair is explicit at 1286 |
| Tr03 | B#v#h#m; B#h2#m; B#v#m2; B#h#m2 | Implemented: 4 | Complete the H representations of B#v2#m and add the one-localized-bond/two-MC-electron pattern. B2H6, B4H10 and B5H9; HW 1236–1238, Figs. 16.17–16.18 explicitly separate BH/BB 2c2e from BHB/BBB 3c2e bonds. Use [1,1,0] permutations for the 3c2e units (PW); B5H9 also uses Tr04. B#v#m2/B#h#m2 also supplies the neutral BH vertices in Tr06's cage descriptions |
| Tr04 | B#v3#m0; B#v2#h#m0; B#v#h2#m0; B#h3#m0 | Implemented: 4 | Zero-contribution participants in local 3c2e units: B5H9 (HW 1237), B3H8− (1254), and B2H7− = H3B–H–BH3 (1243). For B2H7−, use existing H#c-#m2 between neutral BH3 units, contributions [0,2,0] (PW). These states are distinct from non-MC B#v3 and its H variants |
| Tr05 | B#c-#v2#m2; B#c-#v#h#m2; B#c-#h2#m2 | Implemented: 3 | B3H8−, HW 1254, Eq. 16.12: one ordinary B–B bond, six terminal B–H bonds and two B–H–B bridges. Put −1 on the B without an ordinary B–B bond; it contributes one electron to each bridge, each bridge H contributes one, and the other B contributes zero (PW). The latter two B sites use Tr04 |
| Tr06 | B#c-#v#m3; B#c-#h#m3; C#v#m3; C#h#m3; C#v2#m2; C#v#h#m2 | Implemented: 6 | HW 1234–1235, 1265–1277: reusable BH/CH skeletal contributions. For B6H6²−/B12H12²−, place −1 on two B vertices, contributing three electrons each; remaining neutral BH vertices contribute two (Tr03). C2B10H12 uses two neutral CH contributions of three and ten BH contributions of two. Retaining the ordinary C–C bond in ortho-C2B10H12 leaves two MC electrons per C (M18, Table 1, p. 4669). Integer charge/contribution allocations and whole-cage grouping: PW. Explicit/implicit H variants tested; no dependence on cage size |
| Tr07 | B#v3#t; B#v2#h#t; B#v#h2#t; B#h3#t | Implemented: 4 | BH3·NH3/BH3·NMe3 and BF3·OEt2, with mixed H representations; HW 1242–1243, 1281, 1284 explicitly describes donor adducts. Also covers B2Cl4·2D (1285), BR3·D (1315) and the corresponding Lewis-acid functionality |
| Tr08 | B#c+#v2#t2; B#c+#v#h#t2; B#c+#h2#t2 | Implemented: 3 | [BH2(NH3)2]⁺ in diammoniate of diborane, HW 1243; [Ph2B(bipy)]⁺, 1318. Two ordinary bonds and two donated pairs; not a bare two-coordinate borenium proposal |
| Tr09 | Correct B#c-#v3#a2 to B#c-2#v3#a2; remove B#c-#v2#h#a2; add localized companion B#c-2#n#v3 | Implemented: 1 correction + 1 removal + 1 addition | Borolide dianion, HW 1320–1321, Fig. 16.88. Three σ contributions plus two aromatic electrons require five B-assigned electrons, hence B²−. The corrected entry is covered by 1-(dimethylamino)borolide in aromatic/nme2-borolide-2.edn. The localized counterpart retains a B lone pair and two C=C bonds. The B–H entry is removed because reviewed sources do not establish the dianion; B11's neutral NHC adduct instead uses Tr07 |
| Tr10 | B#v2#a#t | Implemented: 1 | Donor-stabilized borabenzene, PMe3→BC5H5; HW 1320. Aromatic graph with B contributing one π electron and accepting one pair (PW); the localized graph uses Tr07's B#v3#t |
| Tr11 | B#c-#u#v3 | Implemented: 1 | Triarylborane radical anions BPh3•−/BMes3•−, HW 1315. Atom-localized charge/spin description; no additional delocalized relation or B–H variants proposed |
| Tr12 | Additional boryl/borylene donor and cluster-adduct states | Not selected | HW 1317–1319, 1322 describes metal complexes and low-valent derivatives. Existing ordinary B states and Tr07–Tr08 cover simpler selected descriptions. S06, Scheme 1 and pp. 113–114, supports retaining the existing B− lone-pair state through the ionic boryl−/Li+ description; the contact dimer does not require extra donor/acceptor states for this example |
| Tr13 | Al: #v2#h, #v#h2, #h3; Ga/In/Tl: #v3, #v2#h, #v#h2, #h3 | Implemented: 15 | Ordinary trivalent centers and monomeric EH3, with explicit/mixed/implicit H. AlH3: HW 1338–1339; GaH3/InH3/TlH3: 1388–1389 (matrix species distinguished from polymers). EX3/ER3: 1393–1394, 1405–1407. Existing Al#v3 retained |
| Tr14 | E#n#v; E#n#h (E = Al, Ga, In, Tl) | Implemented: 8 | Singlet monohydrides EH and monohalides EX; HW 1336, 1339, 1370, 1388–1389, 1395–1396. Also covers ordinary monoorganyl centers; 1407–1408 explicitly retains the lone pair |
| Tr15 | E#u#v2; E#u#v#h; E#u#h2 (E = Al, Ga, In, Tl) | Implemented: 12 | EH2 radicals, explicit/mixed/implicit H. HW 1372 gives AlHMe and alanyl radicals; 1388–1389 gives heavy EH2. Ce03, Tables II–IV and pp. 84–85, explicitly treats neutral AlH2 and GaH2 as ²A1 radicals, resolving the Al#u#h2 source gap |
| Tr16 | E#c-#v4; E#c-#v3#h; E#c-#v2#h2; E#c-#v#h3; E#c-#h4 (E = Al, Ga, In, Tl) | Implemented: 20 | Tetrahydridoaluminate/gallate/indate/thallate: HW 1342–1343, 1392. EX4−/ER4− and substituted alanates: 1343, 1346, 1368, 1394. E#c-#v4 also covers ordinary tetrahedral E− sites |
| Tr17 | Ga#c+#n; In#c+#n | Implemented: 2 | Chemical E(I) counterions in mixed-valent halides such as E⁺EX4−; HW 1387, 1393–1395. Tl#c+#n already exists |
| Tr18 | E#v3#t; E#v2#h#t; E#v#h2#t; E#h3#t (E = Al, Ga, In); Tl#v3#t | Implemented: 13 | EH3·D for Al/Ga/In: HW 1340, 1390; EX3·D and ER3·D: 1346, 1368, 1394, 1407. No TlH3·D extension: HW 1390 discusses that as a prospect, not an obtained adduct |
| Tr19 | E#v3#t2; E#v2#h#t2; E#v#h2#t2; E#h3#t2 (E = Al, Ga, In) | Implemented: 12; Al capacity 12 | AlH3(NMe3)2, GaH3(NMe3)2 and InH3(PCy3)2; HW 1340, 1390. The Al/Ga halide analogues and InX3·2D are also explicit at 1346, 1394. No automatic Tl analogue |
| Tr20 | In#v3#t3 | Implemented: 1 | InX3·3D, D = ether/amine, octahedral; HW 1394. Three ordinary In–X bonds and three neutral donor bonds |
| Tr21 | E#c+3#t6 (E = Al, Ga, In, Tl) | Implemented: 4; Al capacity 12 | [E(H2O)6]³⁺, HW 1345, 1399. The same acceptor states support six halide donors in EF6³− and InCl6³−/TlCl6³−/TlBr6³− (1337, 1344, 1394), and six hydride donors in AlH6³− (1343). Ionic donor allocation: PW. No extra MC or localized #v6 family |
| Tr22 | In#c+3#t5; Tl#c+3#t5 | Implemented: 2 | InCl5²− and TlCl5²−/TlBr5²−, HW 1394. Five halide donors around E³⁺ (PW); reuse existing terminal halide-donor types |
| Tr23 | E#c+#v2#t2 (E = Al, Ga, In); E#c+#v#h#t2 and E#c+#h2#t2 (E = Al, Ga) | Implemented: 7 | E2X6 dimers: HW 1337, 1345–1346, 1394. Two EX2⁺ units joined by two X− donors, each donating to both E centers (PW; existing halide-bridge convention). H variants: Al/Ga amido/halide bridges, 1340, 1390–1391. Al-amido fixtures already supply the Al occupation explicitly; the added states also permit inference |
| Tr24 | Al#c+2#v#t3; Al#c+2#h#t3 | Implemented: 2 | Threefold imido bridges in [RAlNR′]4 and [HAlNR′]4; HW 1368–1369, Fig. 16.112l. AlR²⁺/NR′²− donor allocation (PW), already used by the Al-imido fixtures with supplied Al occupation. No new N states needed |
| Tr25 | Al#v#h#m; Al#h2#m; Ga#v2#m; Ga#v#h#m; Ga#h2#m | Implemented: 5 | Dibridged Al2H6/Ga2H6: HW 1339, 1389, Fig. 16.126a; neutral [1,1,0] bridge allocations as for diborane (PW). Al#v2#m already exists for Al2Me6. HW qualifies the In2H6 structure as tentative; no In analogue selected on that basis |
| Tr26 | Al#m3 | Implemented: 1 | Polymeric α-AlH3, HW 1339: six Al–H–Al bridges at each Al; the recurring local 3c2e pattern, with three total Al contributions per center (PW). The fixture resolves the central Al in a local coordination fragment with supplied cut-boundary Al occupations. No whole-polymer relation or dependence on aggregate size |
| Tr27 | E#c-2#n#v3 (E = Al, Ga, In, Tl) | Implemented: 4 | Tetrahedral E4⁸− Zintl units, HW 1384: Al/Ga in Sr8Al7/Sr8Ga7/Ba8Ga7, In/Tl in Na2In/Na2Tl. Ordinary E–E edges and one lone pair per E²− (PW). E− diamond-like networks reuse Tr16's #v4; no additional cage-aromatic states |
| Tr28 | E#c+#n#t3 (E = Al, Ga, In, Tl) | Implemented: 4 | Cp*Al and gaseous CpGa/CpIn/CpTl monomers; HW 1371, 1408–1409. Six-electron Cp− donation to E⁺ retaining s². The metal's accepted-pair count resolves; donor-side multidonor accounting remains separate work. No MC reinterpretation or polymer-contact additions |
| Tr29 | E#v#d (E = Ga, In, Tl); Al#c+#t3#d | Implemented: 4 | Trielylene Lewis bases: RE→B(C6F5)3, HW 1408; Cp*Al→B(C6F5)3, 1371. Donation consumes the retained lone pair; the Al entry also retains Cp donation. No EH donor series inferred |
| Tr30 | Tl#c+#v2; In#c+#v2#t4 | Implemented: 2 | Dialkylthallium cations and [Me2In(H2O)4]⁺, HW 1407. [Me2Ga(H2O)2]⁺ uses Tr23. The weakly hydrated, linear Me2Tl⁺ description does not require a ladder of Tl acceptor counts |
| Tr31 | Ga#c-#v3#a; Ga#v3#a0 | Implemented: 2 | Gallatabenzene (HW 1411–1412) and Ga3R3²− rings (1416–1417). For the latter, two Ga− sites contribute one π electron each and the neutral Ga contributes zero (PW). Localized contributors use Tr13/Tr16; no additional MC state |
| Tr32 | Larger borane/heteroborane adducts, metallaboranes and triel clusters | Not selected | HW 1269–1278, 1322–1326, 1373–1376, 1413–1421. Retain the extraction without generating types from every cluster electron count. B4Cl4's explicitly stated BBB 3c2e bonding (1280) can reuse the neutral local patterns above; Al/Ga/In tetrahedranes and larger metalloids do not by themselves justify another default series |
| Tr33 | Higher coordination and specialized donor combinations | Not selected | Crown-bound Al/Ga/In/Tl ions, larger borohydride adducts and nine-coordinate solid TlF3; HW 1342, 1372, 1394–1395, 1407. Solid TlF3 can retain the ionic Tl³⁺/F− description using existing states; no nine-donor row or per-complex acceptor ladder is proposed |
| Tr34 | E#c-#n#v2; E#c-#n#v#h; E#c-#n#h2 (E = Al, Ga) | Implemented: 6 | Calculated singlet AlH2− and GaH2−, Ce03, Tables II–VII and pp. 84–86. Two ordinary E–H bonds and one E lone pair; explicit/mixed/implicit H. These are computationally supported binary anions, not claimed isolated salts. HW 1416–1417 also describes [RGa–GaR]²− with a single Ga–Ga bond and retained lone pairs, using Ga#c-#n#v2. No extrapolation to In/Tl |
| Tr35 | B#c-#u2#v2; B#c-#u2#v#h; B#c-#u2#h2 | Implemented: 3 | Calculated triplet BH2− (³B1), Ce03, Tables I, IV–V and pp. 81–84. Retain the existing singlet descriptions as well: the final ZPE-corrected estimate favors singlet by only about 0.01 eV, with ordering sensitive to the calculation. Two ordinary B–H bonds and two unpaired electrons; no MC relation needed |
| Tr36 | E#n#m (E = Al, Ga, In, Tl) | Implemented: 4 | H-bridged E2H2, HW 1339, 1389, Fig. 16.126c. Two local 3c2e bridges with [1,1,0]/[0,1,1] contributions and one retained lone pair per E (PW); existing H#m suffices. The photogenerated HE=EH isomers use ordinary Tr13 states |
| Tr37 | N#v2#a#d | Implemented: 1 companion | General pyridine-type aromatic N donor; Tu25, Fig. 2, supplies InCl3(py)3. Donation consumes the lone pair while retaining one aromatic electron. Paired aromatic/Kekulé fixtures; the latter uses existing N#v3#d |

#### Counts comparison and coverage

The reviewed set contains 163 descriptions: 162 additions and one replacement.
All are implemented and pass the valence invariants. The original 35 triel states
included two borolide entries with orbital/electron counts 8/7; Tr09 corrects one
and removes the other.

The boron batch has 32 fixtures covering all 31 B/C additions, plus the earlier
aromatic borolide fixture. Its 128 outputs agree across both models/policies and
satisfy the valence invariants. Strict atom typing of the earlier BN fixture now
exposes B− lone-pair/two-unpaired-electron alternatives; MostSaturated retains its
previous molecule and records the tie-break. Counts is unchanged.

The heavy-triel batch has 137 fixtures covering all 130 triel additions and the
pyridine-donor N companion. All 548 outputs agree across both models/policies
and satisfy the valence invariants; all 1215
preceding snapshots are unchanged. The donor fixtures use characterized ligands:
InCl3(py)3, TlPh3·NMe3 and ArE→B(C6F5)3 (E = Ga, In, Tl), with
Ar = 2,6-Dipp2-C6H3 and Dipp = 2,6-iPr2-C6H3. InCl3(py)3 has paired aromatic/Kekulé
fixtures. The α-AlH3 fixture resolves the central Al in a local coordination
fragment with supplied cut-boundary occupations. Cp fixtures test metal acceptance;
donor-side multidonor accounting remains separate work.

| Check | Finding / disposition |
| --- | --- |
| Actual-element capacity | Al capacity increased from 8 to 12. The four Tr19 states for AlH3·2D occupy 10 electrons, and Al#c+3#t6 for [Al(H2O)6]³⁺ occupies 12. These selected five-/six-coordinate adducts and ions have consistent electron accounting; the former bound excluded them. The approved bound applies generally, with no molecule-specific exception. Maximum localized valence, charge/H enumeration bounds and saturation targets are unchanged |
| Aromatic ValenceTable entries | No new entry needed. Corrected B²− reads N's aromatic contribution 2; B− reads C's 1; neutral borabenzene with one acceptor reads N's 1. Ga− reads Ge's 1 and neutral Ga reads the existing Ga 0 |
| Free-H completion, ordinary Tl | The first-target rule for Tl [1,3] cannot infer Tl#v#h2, Tl#h3, Tl#u#v#h or Tl#u#h2 with free H at those localized valences. Supplied H counts permit these descriptions; the lower target is retained |
| Free-H completion, InH4−/TlH4− | Shifted Sn/Pb [2,4] lookup cannot infer #v2#h2, #v#h3 or #h4 for either anion with free H at those valences. Supplied H counts permit them |
| Free-H completion, EH3·2D | Shifted Cl/Br/I [1,3,…] lookup cannot infer #v#h2#t2 or #h3#t2 for Al/Ga/In with free H. Supplied H counts permit them |
| Enumeration bounds | Maximum-valence, charge and implicit-H enumeration bounds do not reject supplied literal fields; no blanket increases made |
| Conformance | Every added state is covered, including H variants and aromatic/localized ring descriptions. Coverage checks bookkeeping; it does not distinguish cage geometries or establish bonding assignments |

#### Triel primary sources

All nine papers below are available in full and retained in materials/valence.
The supplied ME63, Ha03, Wr02 and Wr05 PDFs were checked on 2026-11-03; their
experimental and bonding sections confirm the selected donor fixtures and
assignments. No registry or fixture changes were needed.

| Key | Source | Contribution beyond HW |
| --- | --- | --- |
| Ce03 | Černušák et al., “Group IIIa Hydrides XH2 and XH2− (X = B, Al, Ga): Electron Affinities and Singlet-Triplet Splittings Revisited,” *Collect. Czech. Chem. Commun.* **68** (2003), 75–88. [DOI](https://doi.org/10.1135/cccc20030075), [PDF](../materials/valence/Cernusak-2003-group13-dihydrides.pdf) | Doublet neutral dihydrides, singlet anions and the nearly degenerate BH2− triplet. Supports Tr01, Tr15, Tr34–Tr35; distinguishes calculated states from experimental characterization |
| S06 | Segawa, Yamashita and Nozaki, “Boryllithium: Isolation, Characterization, and Reactivity as a Boryl Anion,” *Science* **314** (2006), 113–115. [DOI](https://doi.org/10.1126/science.1131914), [PDF](../materials/valence/Segawa-2006-boryllithium.pdf) | Crystallography, NMR and reactivity support the B− lone-pair description in Tr01/Tr12. The supplied PDF includes the correction to the gallium reference and SI parameter count |
| M18 | Melichar, Hnyk and Fanfrlík, “A systematic examination of classical and multi-center bonding in heteroborane clusters,” *Phys. Chem. Chem. Phys.* **20** (2018), 4666–4675. [DOI](https://doi.org/10.1039/C7CP07422K), [PDF](../materials/valence/Melichar-2018-heteroborane-bonding.pdf) | Local 3c2e and 4c2e IBOs, ordinary bonds where found (including ortho-carborane C–C), and differences from NBO localization. Integer atom contributions and grouping several bonds into a whole-cage relation are PW choices |
| B11 | Braunschweig, Chiu, Kupfer and Radacki, “NHC-Stabilized 1-Hydro-1H-borole and Its Nondegenerate Sigmatropic Isomers,” *Inorg. Chem.* **50** (2011), 4247–4249. [DOI](https://doi.org/10.1021/ic200533e), [PDF](../materials/valence/Braunschweig-2011-hydroborole.pdf) | The neutral 2-BH adduct has tetracoordinate B and alternating C–C bonds (p. 4248). An NHC→B graph uses Tr07's B#v2#h#t or B#v3#t with explicit H. The paper does not establish a free B–H borolide dianion; its Scheme 1 B–H borole examples are metal-bound |
| Tu25 | Tuna et al., *J. Turkish Chem. Soc. A* **12** (2025), 15–22. [DOI](https://doi.org/10.18596/jotcsa.1590991), [PDF](../materials/valence/Tuna-2025-indium-pyridine.pdf) | Fig. 2 and Table 2 reproduce the InCl3(py)3 structure of Jeffs et al. (1984), [DOI](https://doi.org/10.1107/S010827018400785X), supplying a specific Tr20 example |
| ME63 | Maher and Evans, *J. Chem. Soc.* (1963), 5534–5539. [DOI](https://doi.org/10.1039/JR9630005534), [PDF](../materials/valence/Maher-1963-thallium-amines.pdf) | Pages 5535 and 5537 establish isolated TlPh3·NMe3; elemental analysis and NMR confirm the 1:1 adduct used in Tr18 |
| Ha03 | Hardman et al., *J. Am. Chem. Soc.* **125** (2003), 2667–2679. [DOI](https://doi.org/10.1021/ja028579j), [PDF](../materials/valence/Hardman-2003-gallium-borane-adducts.pdf) | Compound 10, ArGa→B(C6F5)3 with Ar = 2,6-Dipp2-C6H3; pp. 2676–2677 explicitly describe Ga lone-pair donation |
| Wr02 | Wright et al., *J. Am. Chem. Soc.* **124** (2002), 8538–8539. [DOI](https://doi.org/10.1021/ja026285s), [PDF](../materials/valence/Wright-2002-indium-borane-adducts.pdf) | Page 8538 and Fig. 2 identify ArInB(C6F5)3, compound 2, with the same Ar group; p. 8539 explicitly discusses In donor action |
| Wr05 | Wright et al., *J. Am. Chem. Soc.* **127** (2005), 4794–4799. [DOI](https://doi.org/10.1021/ja0432259), [PDF](../materials/valence/Wright-2005-thallium-borane-adducts.pdf) | Compound 3, ArTl→B(C6F5)3 with the same Ar group; Fig. 3 gives its structure, and p. 4799 explicitly assigns Tl lone-pair donation to B |

Tr09's unsupported B–H borolide entry is removed. HW's characterized borolides have
a non-H B substituent, and B11 does not establish the B–H dianion.
Further low-valent donor and specialized cluster types are not selected merely
because papers about them exist.

### Alkaline earths

Source: [HW 103, Band I, chapter XVII](<../materials/valence/Holleman-Wiberg - Anorganische Chemie - 103. Auflage - Band 1.pdf#page=1467>),
printed pp. 1429–1480; chapter references pp. 1481–1482. Extraction: 2026-11-03.
Baseline: Be 2, Mg 4, Ca 2, Sr 2, Ba 2 and Ra 2 states. The approved 76 additions
bring the totals to Be 22, Mg 26, Ca 14, Sr 13, Ba 13 and Ra 2. Selected scope complete.

#### Alkaline-earth work table

**Implemented** identifies approved additions. Unselected and deferred descriptions
are marked explicitly. States use the concrete, zero-valence registry defaults.
HW references are printed Band I pages; **PW** marks our charge/contribution
allocation. Examples identify coverage for the states, not a requirement to
turn every extracted compound into a fixture.

| ID | Atom descriptions / gap | Disposition | Examples and evidence |
| --- | --- | --- | --- |
| Ae01 | Existing 14 Be/Mg/Ca/Sr/Ba/Ra states | Retained; citations updated | Neutral s² atoms: HW 1429. E²⁺ ions: 1432–1434, 1443–1444, 1458. Existing Mg#v2: gaseous MgO/MgX2 and monomeric MgR2, 1443, 1452, 1460. Existing Mg#u#v: MgH doublet, MH23, Table 1. No new isolated E⁺ ion series |
| Ae02 | Be#v2; Ca#v2; Sr#v2; Ba#v2 | Implemented: 4 | Molecular EX2 and EO: HW 1432–1433, 1458–1460. BeR2 and CaR2 monomers: 1439–1440, 1476. Ordinary two-bond descriptions; Mg#v2 already exists |
| Ae03 | Be#u#v | Implemented: 1 | BeCl, HW 1431–1432. One ordinary bond and one unpaired electron; molecular subhalide, not a bare ionization state |
| Ae04 | Be#c-#v3; Be#c-2#v4; Mg#c-#v3; Mg#c-2#v4 | Implemented: 4 | Trigonal [BetBu3]− and tetrahedral [BeMe4]²−, HW 1440; monomeric bulky MgR3− and MgMe4²−, 1453. Ordinary bonds with formal charge on the metal (PW). Hydride formulas alone are not used to infer all H variants; see Ae21 |
| Ae05 | Be#v2#t; Be#v2#t2 | Implemented: 2 | BetBu2(OMe2), BePh2(OMe2)2, BeMe2(PMe3)2 and BeMe2(bipy), HW 1440; BeX2·2D, 1432, 1436. Two ordinary Be bonds plus one/two donor pairs |
| Ae06 | Mg#v2#t; Mg#v2#t2; Mg#v2#t3; Mg#v2#t4 | Implemented: 4 | (Me3Si)2CHMgCl(OEt2), EtMgBr(OEt2)2/MgMe2(tmeda), MeMgBr(THF)3 and Mg(GeMe3)2(DME)2; HW 1451–1452, Fig. 17.13. MgBr2(THF)4 independently supports #t4, 1446. #t3/#t4 require capacity 10/12 |
| Ae07 | Ca#v2#t2; Ca#v2#t4 | Implemented: 2 | Ca[CH(SiMe3)2]2(dioxane)2, Ca(SnMe3)2(THF)4 and MesCaI(THF)4; HW 1476. Ordinary Ca–R/Ca–I bonds plus solvent donors (PW). #t4 requires capacity 12 |
| Ae08 | E#c+#v#t2; E#c+#v#t3 (E = Be, Mg) | Implemented: 4 | Halide-bridged Be2Cl4 and [Me3SiMgBr]2: HW 1433, 1452; donor-solvated [tBuBeHal(OEt2)]2 and [EtMgBr(NEt3)]2: 1440, 1452. One ordinary terminal bond, two bridging X− donations, and zero/one neutral donor (PW). Be#c+#v#t3 also covers [Be(OH)(H2O)3]⁺, 1435. Be amido oligomers, 1435, reuse these states and Ae10 |
| Ae09 | Mg#c+#v | Implemented: 1 | RMg⁺ in the Grignard dissociation scheme, HW 1452–1453, Eq. 17.1. Supports the explicitly donor-omitted reaction representation; not a claim that the solution cation is uncoordinated |
| Ae10 | Be#c+2#t2; Be#c+2#t4 | Implemented: 2 | Two amido donors in monomeric Be(NR2)2, HW 1435, Fig. 17.5c (ionic allocation: PW); four donors in [Be(H2O)4]²⁺/[Be(NH3)4]²⁺, 1432, 1436. #t4 also covers the four halide-donor description of polymeric BeX2 and η5/η1 BeCp2, 1432–1433, 1439. The three Be-amido fixtures now infer Be's occupation and cover #t2 |
| Ae11 | Mg#c+2#t4; Mg#c+2#t5; Mg#c+2#t6 | Implemented: 3 | Four Cl− donors in MgCl4²− (PW), HW 1446; [Mg(OAsMe3)5]²⁺, 1443; [Mg(H2O)6]²⁺/[Mg(NH3)6]²⁺, 1443–1446. MgCp2, 1452, uses two six-electron Cp− donations and the same #t6. #t5/#t6 require capacity 10/12 |
| Ae12 | Ca#c+2#t6; Sr#c+2#t6; Ba#c+2#t6 | Implemented: 3 | Two η5-Cp*− ligands in E(C5Me5)2, HW 1476, Fig. 17.19a: six donated pairs in the existing haptic convention (PW). Ca–EDTA, 1475, provides six ordinary donor sites. Capacity 12 required; no new metal MC type |
| Ae13 | Ca#c+2#t8 | Implemented: 1 | [Ca(H2O)8]²⁺ in CaO2·8H2O, explicitly eight-coordinate at HW 1458. Eight water donors require capacity 16 |
| Ae14 | Be#m2; Mg#m2 | Implemented: 2 | Polymeric BeH2 and (BeMe2)n: HW 1433–1434, 1439; (MgR2)n: 1452, 1475. Recurring local E–H–E/E–C–E 3c2e units, with two total metal contributions per site (PW). S03, p. 3622, identifies solid BeH2 as a three-dimensional network of connected BeH4 tetrahedra, not linear chains. Use existing H#m and C#h3#m or their explicit-H counterparts; no aggregate-size-specific type |
| Ae15 | Be#v#m | Implemented: 1 | (BeEt2)2 and (tBuBeMe)3, HW 1439–1440, Fig. 17.9i,k. One ordinary Be–R bond and one total contribution to the two local alkyl bridges (PW) |
| Ae16 | Be#v#m#t; Be#h#m#t; Mg#v#m#t | Implemented: 3 | [BeH2(NMe3)]2 and [BeMe2(NMe3)]2, HW 1435, 1440; [MgEt2(OEt2)]2 and [HMgCl(THF)]2, 1444, 1452. One ordinary terminal bond, two local 3c2e bridges and one donor per metal (PW). Be's terminal hydride has explicit/implicit H variants; HMgCl has bridging H, not an implicit terminal H |
| Ae17 | Be#v2#m0; Mg#v2#m0 | Implemented: 2 | [R2E(μ-H)]2²−, HW 1440, Fig. 17.10r, and 1453, Fig. 17.14m. Two neutral ER2 units joined by two H− bridges; each bridge uses [0,2,0] and existing H#c-#m2 (PW). Zero-contribution MC membership is distinct from ordinary E#v2 |
| Ae18 | Mg#c-#v2#m | Implemented: 1 | [Mg2R6]²−, HW 1453, Fig. 17.14m, explicitly compared with Al2R6. Two ordinary Mg–R bonds and one total Mg contribution to two alkyl 3c2e bridges; −1 on each Mg (PW). Existing neutral bridging-carbon types suffice |
| Ae19 | E#v#h; E#h2 (E = Be, Mg, Ca, Sr, Ba) | Implemented: 10 | Molecular EH2: S03, pp. 3622, 3626–3627; WA04M, pp. 11514–11516; WA04H, pp. 11506–11509. Linear BeH2/MgH2 and bent or nearly linear Ca/Sr/Ba monomers all use two ordinary M–H bonds. These are the mixed/implicit-H counterparts of E#v2; distinct from bulk ionic hydrides |
| Ae20 | E#u#h (E = Be, Mg, Ca, Sr, Ba); E#u#v (E = Ca, Sr, Ba) | Implemented: 8 | MH23, Introduction and Table 1, identifies X²Σ⁺ ground states for all five monohydrides. Full text reviewed; the ordinary M–H graph places the remaining unpaired electron on M (PW). Be's explicit-bond radical is in Ae03; Mg's already exists. No Ra extension inferred |
| Ae21, Be | LiBeH3, Li2BeH4 | Existing states; no further default additions proposed | Wu81's explicit-H BeH3−/BeH4²− graphs reuse Be#c-#v3 and Be#c-2#v4. Hu08, pp. 234105-3–4, 234105-7, describes BeH4 chains with ionic Li+. Two ordinary terminal Be–H bonds and H− bridges [0,2,0] reuse Be#v2#m0 and H#c-#m2 (allocation: PW); the interior-site check resolves under both models/policies. Ba90 re-indexes old powder data and explicitly leaves its perovskite proposals hypothetical; no extra type follows from that paper. Li-bound monomer encodings and additional implicit-H variants remain unselected |
| Ae21, Mg | KMgH3 | Existing states; checked | HW 1445; Fo07, section IV, p. 214118-2, assigns ionic bonding with hydride anions; section V, p. 214118-5, finds effective charges close to nominal ionic values. The K+ + Mg²⁺ + 3H− description reuses existing states and resolves under both models/policies. No isolated MgH3− or lattice MC state is needed |
| Ae22 | Be#c+2#t6 | Implemented: 1; capacity 12 | BeCp*2: Co03, pp. 4455–4456, Fig. 3 and Table 2 establish the η5/η5 crystal structure; pp. 4457–4458 analyze ionic and donor–acceptor interactions between Be²⁺ and Cp− fragments. Two three-pair Cp− donations use the existing haptic convention (PW), not a measured Be shell population. The η5/η1 compounds reuse Ae10's #t4; HW 1439–1440 |
| Ae23, gas | Be(BH4)2 | Existing states; checked | HW 1435: HB(μ-H)3Be(μ-H)3BH. Six local 3c2e bridges reuse Be#m2, B#v#m2/B#h#m2 and H#m (allocation: PW). Explicit/implicit terminal-H graphs resolve under both models/policies |
| Ae23, adduct | Be(BH4)2·D | Existing states; tentative HW structure checked | HW 1435, Fig. 17.6e proposes a Be–B–B ring with three H bridges, one terminal Be–H, two terminal B–H per B and one donor at Be. It reuses Be#v#m#t/Be#h#m#t and B#v2#m/H variants (allocation: PW). NMe3 example checked with explicit/implicit terminal H; this establishes representability of HW's tentative graph |
| Ae23, solid/Cp | Solid Be(BH4)2; CpBeBH4 | No additional default types proposed | HW 1435, 1439–1440; vS08, Fig. 1. The tested neutral local-3c2e polymer description needs B#m3 at its chain B sites; Be#m2 already exists. The Cp−/neutral-B, double-bridge description needs Be#c+#m#t3. Both pass Counts and the invariants; these limited non-binary uses alone do not justify the extra default entries |
| Ae24, seven donors | Ca#c+2#t7; Sr#c+2#t7 | Implemented: 2 | BPS26, sections 3.2–3.3, assigns seven-coordinate aqua ions in concentrated CaCl2/SrCl2 solutions. Seven water donations (PW) use existing O types for the monomeric description and occupy 14 electrons, within the approved Ca/Sr capacity 16 |
| Ae24, eight donors | Sr#c+2#t8; Ba#c+2#t8 | Implemented: 2 | PS95 redetermines [Sr(H2O)8](OH)2 and [Ba(H2O)8](OH)2: eight terminal waters, with hydrogen bonding to the outer shell. Ordinary aqua ions supply a general use independent of the bridging-water lattices. Existing O donor types suffice; capacity 16 required |
| Ae24, nine donors | Ba#c+2#t9 | Implemented: 1 | BPS26, section 3.4, fits nine-coordinate hydrated Ba²⁺ in saturated BaBr2 solution. Nine water donations (PW), using existing O types, require the approved Ba capacity 18 |
| Ae24, bridges | Ca/Sr#c+2#t9; O#v2#d2 and its H variants | Deferred with Ch02 | HW 1474–1475: nine-coordinate hydrated Ca/Sr/Ba chains and eight-coordinate SrCl2·2H2O/BaCl2·2H2O layers with bridging water. The latter have four chloride and four water neighbours per metal. No new O entry or extra metal state is proposed for these lattices in this pass |
| Ae25 | MgEt2(18-crown-6), [RMg(15-crown-5)]⁺, inverse crowns, larger metallocene assemblies | Not selected | HW 1452–1453, 1476–1477. Retain the examples without proposing Mg#v2#t6, Mg#c+#v#t5 or additional acceptor/MC ladders for the specialized assemblies |
| Ae26 | Ionic salts, extended hydrides, suboxides/subnitrides and Ra compounds | No additional state proposed | HW 1458–1474. Existing E²⁺ states cover the ionic descriptions. Do not turn lattice coordination or metallic electrons into extra MC states. Ra remains at its existing atom/Ra²⁺ descriptions in this pass |
| Ae27 | Be#h#m; E#v#m and E#h#m (E = Mg, Ca, Sr, Ba); Sr#m2; Ba#m2 | Implemented: 11 | WA04M, pp. 11516–11518, establishes HMg(μ-H)2MgH; WA04H, pp. 11507–11509, establishes HCa(μ-H)2CaH and HSr(μ-H)3Sr/HBa(μ-H)3Ba, and summarizes the analogous Be dimer. Use local 3c2e bridges with H#m (PW): each metal contributes 1 across the two bridges of E2H4 (E = Be, Mg, Ca); the terminal-H and bare sites of the Sr/Ba dimers contribute 1 and 2 across three bridges. Be#v#m already appears in Ae15. The cyclic Mg3H6/Mg4H8 structures, WA04M pp. 11514, 11517–11518, reuse Mg#v#m/H variants |
| Ae28 | Mg#u#m | Implemented: 1 | HMg(μ-H)2Mg radical: WA04M, p. 11519 and Tables 3–4. Two local 3c2e bridges; terminal-H Mg uses Mg#v#m/H, and the other Mg contributes 1 with one unpaired electron (PW). The paper's linear Mg–Mg–H and H–Mg–Mg–H structures use ordinary Mg–Mg bonds and reuse Mg#u#v, Mg#v2 and Ae19's Mg#v#h; no further type |

#### Alkaline-earth Counts comparison

All 76 additions are distinct from the previous registry and covered by inputs
that leave the target atom's occupation to be inferred. The 80 affected fixtures
(77 new, three revised) produce identical, determined outputs under AtomTyping
and Counts with both policies; all 320 outputs satisfy the valence invariants.

| Limit | Implemented change / verification |
| --- | --- |
| Actual-element capacity | Approved bounds are Be/Mg 12, Ca/Sr 16 and Ba 18; Ra stays at 8. Selected donor/haptic graphs exceed the previous bound 8. The final Sr/Ba increases cover ordinary seven-/eight-/nine-donor aqua ions without adding oxygen types. These counts follow the selected graph convention; they do not assert a physical shell population |
| ValenceTable entries | Unchanged: all six have target_covalences = [2]. No aromatic additions. The selected states resolve with supplied H counts without target or aromatic-table changes |
| Enumeration bounds | Unchanged: maximum localized valence 2 and charge interval [0,2] restrict enumeration, not supplied #v3/#v4 or negative charge |
| Coverage | All Ae02–Ae20, Ae22, Ae24 and Ae27–Ae28 additions are covered. Existing H/C/N/O/halide ligand types supply the companions. Polymer/network fixtures infer an interior site's occupation with supplied boundary occupations. Cp* fixtures use the existing haptic convention; general multidonor accounting remains separately deferred |

The Be follow-up checked seven further inputs using the default registry: BeH3−,
Li2BeH4 with separate Li+ counterions, a LiBeH3 chain fragment, and explicit/implicit
terminal-H forms of gaseous Be(BH4)2 and HW's tentative NMe3 adduct. All 28 outputs
agree across both models/policies and satisfy the invariants. The two additional
solid/Cp graph probes confirm the specific registry gaps recorded in Ae23; Counts
resolves both under both policies. These checks add no registry rows or permanent
fixtures; the approved states already have conformance coverage.

Ae24 adds five donor_acceptor fixtures. All 20 outputs agree across both models
and policies, satisfy the invariants and need no tie break. The 1424 earlier
snapshots are unchanged. The KMgH3 ionic probe also agrees across all four
configurations and satisfies the invariants; it needs no additional registry row.

#### Alkaline-earth primary-source checks

The sources below are reviewed and archived in materials/valence.

| Key | Source | Needed for |
| --- | --- | --- |
| S03 | Shayesteh et al., “Infrared emission spectra of BeH2 and BeD2,” *J. Chem. Phys.* **118** (2003), 3622–3627. [DOI](https://doi.org/10.1063/1.1539850), [PDF](../materials/valence/Shayesteh-2003-BeH2-BeD2.pdf) | Reviewed: gaseous linear BeH2, pp. 3622, 3626–3627; solid BeH2 network, p. 3622. Ae14, Ae19 |
| WA04M | Wang and Andrews, “Infrared Spectra of Magnesium Hydride Molecules, Complexes, and Solid Magnesium Dihydride,” *J. Phys. Chem. A* **108** (2004), 11511–11520. [DOI](https://doi.org/10.1021/jp046410h), [PDF](../materials/valence/Wang-Andrews-2004-Mg-hydrides.pdf) | Reviewed: MgH2 and Mg2H4, pp. 11514–11517; cyclic oligomers, pp. 11517–11518; Mg–Mg bonds and the Mg2H3 radical, pp. 11518–11520. Ae19–Ae20, Ae27–Ae28 |
| WA04H | Wang and Andrews, “Metal Dihydride (MH2) and Dimer (M2H4) Structures in Solid Argon, Neon, and Hydrogen (M = Ca, Sr, and Ba): Infrared Spectra and Theoretical Calculations,” *J. Phys. Chem. A* **108** (2004), 11500–11510. [DOI](https://doi.org/10.1021/jp046046m), [PDF](../materials/valence/Wang-Andrews-2004-Ca-Sr-Ba-hydrides.pdf) | Reviewed: monomers, pp. 11506–11509; isotope evidence for dibridged Ca2H4 and tribridged Sr2H4/Ba2H4, pp. 11507–11509. Ae19, Ae27 |
| MH23 | Pang et al., “Theoretical Investigation of Spectroscopic Properties of the Alkaline-Earth-Metal Monohydrides toward Laser Cooling and Magneto-Optical Trapping,” *ACS Omega* **8** (2023), 19391–19401. [DOI](https://doi.org/10.1021/acsomega.3c00352), [PDF](../materials/valence/Pang-2023-alkaline-earth-monohydrides.pdf) | Reviewed: X²Σ⁺ monohydride ground states and experimental references, pp. 19392–19393, Table 1. Ae20 |
| Co03 | Conejo et al., “Synthesis, Solid-State Structure, and Bonding Analysis of the Beryllocenes [Be(C5Me4H)2], [Be(C5Me5)2], and [Be(C5Me5)(C5Me4H)],” *Chem. Eur. J.* **9** (2003), 4452–4461. [DOI](https://doi.org/10.1002/chem.200304876), [PDF](../materials/valence/Conejo-2003-beryllocenes.pdf) | Reviewed: η5/η5 BeCp*2 crystal structure, pp. 4455–4456; ionic/donor–acceptor fragment analysis, pp. 4457–4458. Ae22; the η5/η1 compounds reuse Ae10 |
| Wu81 | Würthwein, Krogh-Jespersen and Schleyer, “Structures of Complex Beryllium Hydrides and Fluorides, LiBeX3 and Li2BeX4,” *Inorg. Chem.* **20** (1981), 3663–3667. [DOI](https://doi.org/10.1021/ic50225a016), [PDF](../materials/valence/Wuerthwein-1981-complex-beryllium-hydrides-fluorides.pdf) | Reviewed: calculated BeX3−/BeX4²− anions, pp. 3664–3665; LiBeX3 and Li2BeX4 monomers, pp. 3664–3667. Ae21 |
| Ba90 | Bastide, “Crystal structure of lithium beryllium hydrides LiBeH3, Li2BeH4,” *Solid State Commun.* **74** (1990), 355–358. [DOI](https://doi.org/10.1016/0038-1098(90)90502-3), [PDF](../materials/valence/Bastide-1990-lithium-beryllium-hydrides.pdf) | Reviewed: revised indexing and lattice parameters, pp. 356–357; hypothetical status stated in the conclusion, p. 357. No bond-electron assignment or new registry support |
| Hu08 | Hu et al., “Crystal structure prediction of LiBeH3 using ab initio total-energy calculations and evolutionary simulations,” *J. Chem. Phys.* **129** (2008), 234105. [DOI](https://doi.org/10.1063/1.3021079), [PDF](../materials/valence/Hu-2008-LiBeH3.pdf) | Reviewed: corner-sharing BeH4 chains, pp. 234105-3–4; ionic Li+/chain interaction, pp. 234105-6–7, 234105-11. The local H− bridge allocation is PW |
| vS08 | van Setten, de Wijs and Brocks, “First-principles calculations of the crystal structure, electronic structure, and thermodynamic stability of Be(BH4)2,” *Phys. Rev. B* **77** (2008), 165115. [DOI](https://doi.org/10.1103/PhysRevB.77.165115), [reviewed preprint](../materials/valence/van-Setten-2008-BeBH4-preprint.pdf) | Preprint p. 2, Fig. 1 distinguishes chain and pendant BH4 groups in the solid, consistent with HW 1435. Local 3c2e allocation is PW |
| PS95 | Persson, Sandström, Yokoyama and Chaudhry, “Structure of the Solvated Strontium and Barium Ions in Aqueous, Dimethyl Sulfoxide and Pyridine Solution, and Crystal Structure of Strontium and Barium Hydroxide Octahydrate,” *Z. Naturforsch. A* **50** (1995), 21–37. [DOI](https://doi.org/10.1515/zna-1995-0105), [PDF](../materials/valence/Persson-1995-Sr-Ba-solvation.pdf) | Eight-coordinate Sr/Ba hydroxide octahydrate crystal structures; solution measurements in Table 9. Supports Ae24's #t8 states independently of water bridges |
| BPS26 | Bajnóczi, Sigfridsson Clauss and Persson, “Structure of Concentrated Aqueous Solutions of Alkaline Earth Metal Chlorides and Bromides,” *J. Solution Chem.* **55** (2026), 851–865. [DOI](https://doi.org/10.1007/s10953-025-01559-9), [PDF](../materials/valence/Bajnoczi-2026-alkaline-earth-solutions.pdf) | Sections 3.2–3.4: seven-coordinate Ca/Sr and nine-coordinate Ba aqua-ion descriptions from solution scattering/EXAFS; also distinguishes concentrated Ca water bridges. Supports Ae24's #t7/#t9 states |
| Fo07 | Fornari, Subedi and Singh, “Structure and dynamics of perovskite hydrides AMgH3 (A = Na, K, Rb) in relation to the corresponding fluorides: A first-principles study,” *Phys. Rev. B* **76** (2007), 214118. [DOI](https://doi.org/10.1103/PhysRevB.76.214118), [PDF](../materials/valence/Fornari-2007-perovskite-hydrides.pdf) | Reviewed: ionic hydride description, section IV, p. 214118-2; effective charges near nominal ionic values, section V and Table IV, p. 214118-5. Existing K+, Mg²⁺ and H− states suffice for Ae21 |

### Alkali metals

Source: [HW 103, Band I, chapter XVIII](<../materials/valence/Holleman-Wiberg - Anorganische Chemie - 103. Auflage - Band 1.pdf#page=1521>),
printed pp. 1483–1531; chapter references pp. 1531–1533. Extraction: 2026-11-03.
Baseline: Li 3, Na/K/Rb/Cs/Fr 2 each, 13 states in total.

#### Alkali-metal work table

All 34 approved additions are implemented: Li 9, Na 7, K 7, Rb 5 and Cs 6;
none for Fr. This includes the three higher-crown states in Ak15.
States use the concrete, zero-valence registry defaults. HW references are
printed Band I pages; **PW** marks our graph allocation. No new MC entry is added.
Default examples are aqua ions, ordinary donor solvates, crowns/cryptands and
standard organolithium compounds. Specialized carbon-contact arrangements and
Cp stacks are retained as extracted cases but do not justify further types.

| ID | Atom descriptions / gap | Disposition | Examples and evidence |
| --- | --- | --- | --- |
| Ak01 | Existing M#u and M#c+ (Li–Fr); Li#v | Retained; citations updated | Atomic s¹ occupations: HW 1483, 1496; M⁺: 1485, 1504–1505. Li2, LiX and monomeric LiR: 1484, 1487, 1494. All 13 existing entries audited |
| Ak02 | M#v (M = Na, K, Rb, Cs) | Implemented: 4 | M2 vapor molecules with a filled σ orbital, HW 1500; gaseous MX monomers, 1508. Ordinary one-bond descriptions |
| Ak03 | M#h (M = Li, Na, K, Rb, Cs) | Implemented: 5 | Molecular MH, distinct from bulk ionic hydrides. ADD09, section 2 and Table 2, treats the X¹Σ⁺ ground states of all five molecules and compares with spectroscopy. Implicit-H counterparts of M#v; ordinary M–H allocation: PW |
| Ak04 | M#c-#n (M = Na, K, Rb, Cs) | Implemented: 4 | Isolated alkalide anions in characterized salts, HW 1505–1506. Closed s² description. HW 1485 does not establish condensed-phase Li−; no Li− or Fr− extension proposed |
| Ak05 | Li#c+#t; Li#c+#t2; Li#c+#t3; Li#c+#t4 | Implemented: 4 | One solvent donor per Li in (LiMe)4·4THF; two in (LiPh·TMEDA)2; three in LiR·3THF; four in separated [LiD4]⁺ and [Li(H2O)4]⁺. HW 1485, 1492–1494. Li⁺/R− ionic contacts plus neutral donors (PW); carbon retains its anionic lone pair |
| Ak06 | Li#v#t3 | Implemented: 1 | Ordinary Li–C bond plus three N/O donors in monomeric LiR·3THF or LiPh·PMDTA, HW 1493–1494, Fig. 18.5k. Supports the bonded monomer notation alongside Ak05's ion-pair notation |
| Ak07 | Li#c+#t5 | Implemented: 1; minimum capacity 10 | LiNCS(12-crown-4), HW 1485: four O donors and N-bound NCS−; existing nitrile N and terminal S− types suffice |
| Ak08 | M#c+#t6 (M = Li, Na, K) | Implemented: 3; minimum capacity 12 | Li(DME)3⁺: MP12, p. 429 and SI Table S4, with 75 six-coordinate DME-solvate structures. Na(H2O)6⁺: MP12, pp. 431, 436; Na(H2O)(benzo-15-crown-5)⁺ and K(18-crown-6)⁺: HW 1525. K–valinomycin, 1526, and K(THF)6⁺, MP12 p. 430, reuse the K state. No Cp example is needed |
| Ak09 | Na#c+#t2; Na#c+#t3; Rb#c+#t3 | Not selected | The collected examples are [Na(TMEDA)]⁺ in the LiPh4³− salt (HW 1494), (NaPh·PMDTA)2 and [RbCH(SiMe3)2·PMDTA]2 (1529). Their split between N donation and ionic carbon contacts alone does not justify these extra default counts. |
| Ak10 | Na#c+#t5; Na#c+#t7 | Implemented: 2; minimum capacity 10/14 | Na(15-crown-5)⁺ and Na(H2O)(18-crown-6)⁺, HW 1525. Five ether donors, or six ether donors plus one terminal water |
| Ak11 | Rb#c+#t7; Cs#c+#t7 | Implemented: 2; minimum capacity 14 | Rb/Cs complexes of 21-crown-7, explicitly described at HW 1525. Seven ordinary ether-donor pairs |
| Ak12 | M#c+#t8 (M = Na, K, Rb, Cs) | Implemented: 4; minimum capacity 16 | Na/K crypt-222 salts, HW 1506; Rb crypt-222, 1525. MMW73, pp. 392–394, explicitly gives six O and two N donors for Rb/Cs crypt-222, with outer-sphere SCN− and water. MP12, p. 436, also selects eight-water descriptions for hydrated Rb⁺/Cs⁺. Existing N/O donor types suffice |
| Ak13 | (LiR)n, heavier MR aggregates, amides, alkoxides and molecular/lattice MX or MH aggregates | Existing ionic states; no aggregate-specific MC additions | HW 1490–1495, 1507–1508, 1529–1531. HW 1490 explicitly describes ionic Li⁺/R− subunits and noncovalent aggregation; 1529 gives the heavier organyls' ionic description. Ion contacts preserve the aggregate topology (PW). (LiMe)4 was checked with the current registry; donor solvates additionally use the new acceptor states |
| Ak14 | Li#c-#v2; Li#c-3#v4; Na#c-#v2 | Not selected | Localized alternatives for LiTsi2−, LiPh4³− and NaTsi2−, HW 1493–1494, 1529. Ionic contact graphs use existing M⁺/C− types; the Li anion-donor descriptions can reuse Ak05. These specialized examples justify no extra metal-charge states or fixtures |
| Ak15 | Li#c+#t8; K#c+#t10; Cs#c+#t12 | Implemented: 3; minimum capacities 16/20/24 | Ordinary crown donors: Li(12-crown-4)2⁺, HW 1485; K(benzo-15-crown-5)2⁺ and K(dibenzo-30-crown-10)⁺, 1525; Cs(18-crown-6)2⁺, 1506. Eight, ten and twelve O donors, respectively. These remain in the ordinary crown family; their capacities are implemented |
| Ak16 | Bridging neutral water donors, including O#v2#d2 and H variants | Deferred with Ch02 | Hydrated alkali-halide/sulfate/carbonate chains and dimers, HW 1523–1524. Do not infer terminal aqua coordination from a hydrate's stoichiometric water count. No additional O types proposed here |
| Ak17 | Metal-local states for electrides, ligand-reduced cryptates, metal clusters/suboxides, higher bare cations and further Fr states | No additions selected | HW 1500, 1504–1506, 1516–1517, 1527, 1530–1531. Keep interstitial/ligand electrons distinct from metal-local occupations. The current Fr atom/cation descriptions are retained |
| Ak18 | Cs#c+#t6; Li#c+#t7; Cp/arene assemblies as additional default examples | Not selected | LiCp2−, LiCp(TMEDA), LiCp–crown, NaCp2−, KCp chains, Cs2Cp3− and CsTsi–benzene; HW 1485, 1493–1494, 1529–1530. No further row or fixture is justified solely by these carbon-ligand arrangements. Li/Na/K #t6 have independent solvent/crown support in Ak08 |
| Ak19 | K#c+#t7 | Implemented: 1; minimum capacity 14 | K(H2O)7⁺: MP12, pp. 431, 436, adopts sevenfold hydration from measured distances and their comparison with structural/thermodynamic data. Seven ordinary water donations (PW); existing O#h2#n#d suffices. This is the paper's selected integer description, with the fitted hydration number fixed during final refinement |

#### Alkali-metal Counts comparison

All six ValenceTable entries retain target_covalences = [1]. No aromatic-table
entry is needed. The approved capacities admit the selected donor descriptions.

| Limit | Finding / implemented change |
| --- | --- |
| Actual-element capacity | M⁺#tN occupies 2N electrons in this counting model. Sixteen added states exceeded the former bound of 8. Capacities are now Li/Na/Rb 16, K 20 and Cs 24; Fr remains 8. The bounds support the approved aqua, solvate, crown and cryptand donor graphs, including Ak15's eight-/ten-/twelve-donor crowns; they apply generally without molecule-specific exceptions |
| Saturation targets | Retain [1]. With supplied H counts, the ordinary/donor states need no target change. Li#v#t3 occupies exactly 8 electrons |
| Enumeration bounds | Unchanged. Supplied #h is accepted for K/Rb/Cs despite their zero free-H enumeration bound; the molecular-hydride fixtures confirm this |
| Conformance | All 34 additions have inferred coverage in 35 new fixtures. Their 140 outputs agree across both models/policies, satisfy the valence invariants and require no tie breaks. All 1429 earlier snapshots are unchanged |

The fixtures cover explicit/implicit LiH–CsH, the four alkalides, the selected
solvates and aqua ions, crowns and crypt-222 complexes. (LiMe)4(THF)4 and
(LiPh·TMEDA)2 use ionic carbon contacts; LiPh(PMDTA) has separate ion-pair and
ordinary Li–C fixtures. Each metal's lone-pair and unpaired-electron counts are inferred.

#### Alkali-metal primary sources

- **ADD09:** Aymar, Deiglmayr and Dulieu, “Systematic trends in electronic properties
  of alkali hydrides,” *Can. J. Phys.* **87** (2009), 543–556.
  [DOI](https://doi.org/10.1139/P09-005),
  [reviewed and archived preprint](../materials/valence/Aymar-2009-alkali-hydrides-preprint.pdf).
  Preprint pp. 3–4 and Table 2 (p. 21): molecular LiH–CsH ground states and
  experimental comparison. Supports Ak03; no excited-state or molecular-cation
  series proposed from the additional calculations.
- **MMW73:** Moras, Metz and Weiss, “Etude structurale des cryptates. II. Structure
  cristalline et moléculaire des ‘cryptates de rubidium et de césium’,”
  *Acta Cryst. B* **29** (1973), 388–395.
  [DOI](https://doi.org/10.1107/S0567740873002645),
  [reviewed and archived PDF](../materials/valence/Moras-1973-Rb-Cs-cryptates.pdf).
  Page 392 explicitly identifies all eight ligand heteroatoms as bound to the
  cation, with inward-directed lone pairs; p. 394 describes eightfold coordination
  and absence of additional metal–anion/water contacts. Supports Rb/Cs#c+#t8
  with six ether-O and two tertiary-amine-N donations (Ak12).
- **MP12:** Mähler and Persson, “A Study of the Hydration of the Alkali Metal Ions
  in Aqueous Solution,” *Inorg. Chem.* **51** (2012), 425–438.
  [DOI](https://doi.org/10.1021/ic2018693),
  [reviewed and archived PDF](../materials/valence/Mahler-2012-alkali-hydration.pdf),
  [supporting information](../materials/valence/Mahler-2012-alkali-hydration-si.pdf).
  Page 429 surveys six-coordinate Li ether solvates; SI Table S4 identifies the
  75 structures as DME solvates, supplying the Li(DME)3⁺ fixture. Pages 431 and 436 support
  Na/K/Rb/Cs aqua descriptions with 6/7/8/8 donors. These are inferred integer
  hydration numbers supported by the authors' analysis, not direct counts of
  individual solution molecules. Ak08, Ak12 and Ak19.

### Copper group

Source: [HW 103, Band II, chapter XXII](<../materials/valence/Holleman-Wiberg - Anorganische Chemie - 103. Auflage - Band 2.pdf#page=190>),
printed pp. 1686–1743; chapter references pp. 1744–1746. Extraction: 2026-11-03.
Chapters XIX–XXI supply transition-metal, coordination and solid-state background;
the next element chapter after the alkalis is Cu/Ag/Au. Rg is treated separately by HW.
Baseline: Cu 3, Ag 4 and Au 4 states. All 77 selected additions below are implemented;
current totals are Cu 32, Ag 22 and Au 34. All 11 previous states are retained.
The existing-entry citation audit, capacity changes and conformance coverage are complete.

#### Copper-group work table

M = Cu, Ag, Au; X = halide; D = neutral donor. States use concrete, zero-valence
defaults. HW page references are to Band II. **PW** marks our graph allocation:
ordinary M–X bonds leave X neutral, whereas X−→M donation retains the ligand charge.
Where both appear, the states support those two supplied graph notations.
The accepted policy supports neutral and ionic electron-allocation conventions
for the selected TM families, including hydride and organyl ligands. Neutral
donor ligands remain pair donors; additional mixed partitions need their own justification.
Examples identify the supported families; they are not a requirement to make a fixture
for every compound. **Capacity** is the minimum required by the current Counts formula.

| ID | Atom descriptions / gap | Disposition | Examples and evidence |
| --- | --- | --- | --- |
| Cg01 | Existing M#n5#u; M#c+#n5; Cu/Ag#c+2#n4#u; Ag#c+#n5#t3; Au#c+3#n4 | Retained: 10; citations updated | Atomic d¹⁰s¹: HW 1538, Table 19.1. Cu(I/II): 1693–1694; Ag(I/II): 1714, 1719–1720; Au(I): 1728. Free-ion Au³⁺ d⁸ singlet levels: Z23, Table 2. Ag's #t3 also covers three-donor coordination beyond the existing [Ag2(I2)4]²⁺ example (M16): AgI(PR3)2 in the iodide-donor notation, HW 1714. The remaining existing Au entry is Cg32 |
| Cg02 | M#n5#v; M#n5#h | Implemented: 6 | Molecular MH: HW 1695, 1714, 1733. M2 molecules: 1731. AuF monomer: 1733. Ordinary single bonds; explicit and implicit H |
| Cg03 | M#c+#n5#t; M#c+#n5#t2; M#c+#n5#t3; M#c+#n5#t4 | Implemented: 11; Ag #t3 exists | One to four donor pairs. Cu(CO)n⁺ (n = 1, 3, 4), [Cu(NH3)2]⁺, [Cu(CN)3]²− and [Cu(py)4]⁺: HW 1694–1697. Ag(CO)n⁺ (n = 1–4), [Ag(NH3)2]⁺, AgI(PR3)2, [Ag(PPh3)4]⁺ and Ag(H2O)4⁺: 1713–1718. LAu⁺, Au(CO)n⁺ (n = 1–3) and Au(PR3)n⁺ (n = 2–4): 1734, 1742. Maximum capacity 18 |
| Cg04 | Cu#n5#v#t3; Ag#n5#v#t2; Ag#n5#v#t3; Au#n5#v#t; Au#n5#v#t2; Au#n5#v#t3 | Implemented: 6 | One ordinary terminal M–X/M–R bond plus donors: CuCl(PR3)3 and CuCl(CO)(H2O)2, HW 1696; AgI(PR3)2 and AgX(PR3)3, 1714, 1717; AuCl(PR3)n (n = 1–3) and RAuPR3, 1734, 1742–1743. Maximum capacity 18 |
| Cg05 | M#c-#n5#v2 | Implemented: 3 | Ordinary-bond notation for M(CN)2− and MR2−: HW 1693–1697, 1707–1708, 1715–1717, 1722, 1734, 1742–1743. CuCl2−/CuBr2− also supported by WWS01; all three linear MH2− anions reuse these explicit-H states (AW03H, pp. 11753–11756; Au also HW 1733). Donor notation reuses Cg03 |
| Cg06 | Cu#c-2#n5#v3; Cu#c-3#n5#v4 | Implemented: 2 | Ordinary-bond [Cu(CN)3]²− and [Cu(CN)4]³−, HW 1697, Fig. 22.5e,g; four-coordinate CuCl4³−, 1696. Donor notation reuses Cg03. Capacity 16/18 |
| Cg07 | Cu#n4#u#v2; Cu#n4#u#v2#t2; Cu#n4#u#v2#t4 | Implemented: 3; #t4 needs capacity 21 | Molecular CuCl2/CuBr2 doublets: WWS01, pp. 7388, 7391–7392. CuCl2D2/CuBr2(NH3)2 and CuCl2(H2O)4: HW 1700–1701. Two ordinary bonds with zero, two or four additional donor pairs |
| Cg08 | Cu#c+2#n4#u#t4; Cu#c+2#n4#u#t5; Cu#c+2#n4#u#t6 | Implemented: 3; capacities 17/19/21 | Cu(NH3)4²⁺, CuCl5³−, Cu(H2O)6²⁺ and Cu(NH3)6²⁺: HW 1694, 1700, 1705. Halide-donor allocation: PW. The Cu acetate dimer, 1705, provides another five-donor example without assigning an ordinary Cu–Cu bond |
| Cg09 | Cu#c-2#n4#u#v4; Cu#c-3#n4#u#v5; Cu#c-4#n4#u#v6 | Implemented: 3; capacities 17/19/21 | Ordinary-bond CuCl4²−, CuCl5³− and CuF6⁴−: HW 1700–1701. These are the corresponding halide graph alternatives to Cg08 |
| Cg10 | Cu#c+3#n4#t4; Cu#c-#n4#v4 | Implemented: 2 | Diamagnetic Cu(OH)4−, HW 1707: four OH− donors or four ordinary Cu–O bonds (PW). Other square-planar Cu(III) complexes reuse these counts |
| Cg11 | Cu#c+3#n3#u2#t6; Cu#c-3#n3#u2#v6 | Implemented: 2; capacity 20 | CuF6³−, explicitly high-spin d⁸ at HW 1707; AW69, pp. 1895–1897, assigns the ³A2g ground state. Six F− donors or six ordinary Cu–F bonds (PW). Supply #u2 when checking Counts; its unconstrained occupation completion pairs the eight remaining electrons |
| Cg12 | Cu#c+4#n3#u#t6; Cu#c-2#n3#u#v6 | Implemented: 2; capacity 19 | CuF6²− in Cs2CuF6, explicitly low-spin paramagnetic d⁷ at HW 1707. Donor and ordinary-bond notations (PW) |
| Cg13 | Ag#c+2#n4#u#t4; Ag#c+2#n4#u#t6 | Implemented: 2; capacities 17/21 | Ag(py)4²⁺ and Ag(bipy)2²⁺: HW 1720–1721; Ag(H2O)6²⁺ in acidic solution: 1713, 1720. Six-coordinate AgF6⁴− provides binary-ion coverage, 1719 |
| Cg14 | Ag#c-2#n4#u#v4; Ag#c-4#n4#u#v6 | Implemented: 2; capacities 17/21 | Ordinary-bond AgF4²− and AgF6⁴−, isolated units in BaAgF4/Ba2AgF6: HW 1719. Donor notation reuses Cg13 |
| Cg15 | Ag#c+3#n4#t4; Ag#c-#n4#v4 | Implemented: 2 | Diamagnetic square-planar AgF4−: HW 1721. Donor and ordinary-bond notations (PW); Ag(CF3)4−, 1722, independently supports the ordinary-bond state |
| Cg16 | Ag#c+3#n3#u2#t6; Ag#c-3#n3#u2#v6 | Implemented: 2; capacity 20 | Octahedral AgF6³− in Cs2KAgF6: HW 1714, 1721. AW69, pp. 1895–1897, Table II, assigns a d⁸ ³A2g ground state, supporting #n3#u2; it also cites a 2.6 BM moment. Donor and ordinary-bond notations (PW). Supply #u2 for Counts as in Cg11 |
| Cg17 | Cs2AgF6: Ag(IV) versus mixed Ag(III/V) | Not selected; assignment unresolved | HW 1721 calls this Ag(IV). S78, pp. 247–250, reports weak magnetism and proposes paired Ag(III)/Ag(V), supported by diamagnetic Cs2Ga0.5Ag0.5F6. The authors explicitly leave that interpretation unconfirmed; attempted CsAgF6 synthesis failed. This source does not settle a new Ag(IV), paired six-coordinate Ag(III), or Ag(V) state |
| Cg18 | Au#c-#n6 | Implemented: 1 | Auride Au− in CsAu, NMe4Au and BaAu2: HW 1739–1740. Closed d¹⁰s² count; capacity 12 |
| Cg19 | Au#n4#v3; Au#n4#v3#t | Implemented: 2 | Molecular AuCl3/AuBr3 and AuMe3: HW 1736, 1743; AuCl3(H2O) and AuMe3(PPh3): 1736–1737, 1743. Ordinary terminal bonds, with one additional donor in the adducts |
| Cg20 | Au#c+#n4#v2#t2 | Implemented: 1 | Au2Cl6/Au2Br6 and the primary AuF3 chain: HW 1736–1737. Two ordinary terminal Au–X bonds and two bridging X− donations per Au (PW). Also Me2Au(H2O)2⁺ and [Me2AuBr]2, 1743–1744 |
| Cg21 | Au#c+3#n4#t4; Au#c-#n4#v4 | Implemented: 2 | Au(NH3)4³⁺ and donor notation for AuX4−, HW 1737–1738; ordinary AuX4−/Au(CN)4− and AuMe4−, 1737–1738, 1743. Explicit AuH4− also reuses the latter, 1733 |
| Cg22 | Au#n3#v5; Au#c-#n3#v6; Au#c+5#n3#t6; Au#c+#n3#v4#t2 | Implemented: 4; maximum capacity 18 | AuF5 monomer; ordinary-bond and F−-donor AuF6−; Au2F10 with four ordinary terminal bonds and two bridging F− donations per Au (PW). HW 1738–1739, including the diamagnetic Au(V) description and dimer structure |
| Cg23 | Au#c-#n5#v#h; Au#c-#n5#h2; Au#c-#n4#v3#h; Au#c-#n4#v2#h2; Au#c-#n4#v#h3; Au#c-#n4#h4 | Implemented: 6 | All mixed/implicit-H counterparts of linear AuH2− and square-planar AuH4−, HW 1733, Fig. 22.23d,e. Explicit-H entries are Cg05/Cg21. AW03H, pp. 11754–11757, supports the singlet structures with isotope-resolved IR and calculations |
| Cg24 | Cu#n5#m | Implemented: 1 | (CuCH2SiMe3)4 and related CuR aggregates: HW 1708 explicitly describes local Cu–C–Cu 3c2e bridges. One total Cu contribution across its two bridges (PW); existing C#v3#m and its H variants suffice. No aggregate-size-dependent state |
| Cg25 | Au#n4#u#v2 | Implemented: 1 | Molecular AuCl2/AuBr2, observed as short-lived gas-phase binaries: HW 1735, in the d⁹ Au(II) section. Separate from solid AuCl2, which HW identifies as mixed Au(I/III). Also explicit-H AuH2: AW03H, pp. 11754, 11757–11758, gives the bent ²B2 state and matrix IR evidence |
| Cg26 | Cu/Ag#c-#n5#v#h; Cu/Ag#c-#n5#h2; Au#n4#u#v#h; Au#n4#u#h2 | Implemented: 6 | Mixed/implicit-H variants of linear singlet CuH2−/AgH2− and bent doublet AuH2. AW03H, pp. 11753–11755, 11757–11758: isotope-resolved matrix IR plus calculated structures/states. Explicit-H states are accepted in Cg05/Cg25. No additions for the calculated CuH4−/AgH4− or HW's calculated Au2H6 |
| Cg27 | Cu#c+#n5#t#m0 | Implemented: 1; allocation PW | Stryker's reagent [(PPh3)CuH]6: HW 1694–1696. Fr24, pp. 2, 4–5, supports Cu(I) and edge-bridging hydrides, revising the older face-capping picture. Use six local Cu–H−–Cu bridges with contributions [0,2,0] and one P donation per Cu; existing H#c-#m2 and P#v3#d suffice |
| Cg28 | Au#c+3#n4#t5 | Implemented: 1; capacity 18 | Five-coordinate intermediate in halide substitution at AuX4−, HW 1737–1738. Five anionic donor pairs (PW). Its inclusion rests on the substitution motif; rare isolated chelate complexes are not needed as the justification |
| Cg29 | Au#c+2#n4#u#t4; Au#n4#v3#t2; Au#c+#n4#v2#t4 | Not selected | Au(II) dithiolate/phthalocyanine complexes and five-/six-coordinate Au(III) chelates: HW 1729, 1735, 1737. These limited non-binary examples do not presently justify extra default states; the six-coordinate diarsine description would additionally require capacity 20 |
| Cg30 | Au#n5#m; an analogous Ag organyl MC state | Not selected | HW explicitly assigns 3c2e bridges in (AuMes)5, 1743; (AgMes)6 is listed at 1722 without a corresponding bond assignment. No extra row solely for these individual aggregates. Ordinary AuR/AgR donor adducts and MR2− are already accepted |
| Cg31 | Extended solids, neutral matrix carbonyls, larger metal/chalcogenide clusters and additional bare ions | No further default additions selected | HW 1695–1699, 1702–1706, 1713–1722, 1731–1732, 1739–1742. Ionic salts reuse existing ions; accepted donor counts cover many halide/ligand bridges. No all-metal MC ladder, secondary metallophilic contacts, or scaffold-specific state series. Stryker's reagent is the separately justified Cg27 case |
| Cg32 | Existing Au#c+3#n3#u2 | Retained; covered | Free Au³⁺ (spectroscopic Au IV) has a 5d⁸ ³F₄ ground state: Z23, Table 2. The proposed removal incorrectly applied the paired configuration of Au(III) complexes to the free ion. The existing singlet fixture remains; atomic_ions/au+3-triplet.edn supplies #u2 and infers #n3 |
| Cg33 | η²-H2 adducts (H2)MH and (H2)AuH3 | Deferred; no new states selected | HW 1733 and AW03H, pp. 11753, 11756, distinguish the H2 adducts from ordinary hydrides. Their H2-donor accounting remains separate work; no MC replacement introduced in this pass |

#### Copper-group Counts comparison and coverage

Cu and Ag now have capacity 21; Au retains 18. None has a direct ValenceTable entry.
For these nonaromatic states the checked arithmetic is
`E = 11 − charge + h + v + 2t` and `2n + u = 11 − charge − h − v − 2d − m`.
All 77 additions satisfy the latter and are distinct from the 11 retained entries.

| Constraint | Implemented selection |
| --- | --- |
| Cu capacity | Increased 18 → 21. This admits nine states previously excluded: Cg07 #t4; Cg08 #t5/#t6; Cg09 #v5/#v6; both Cg11 states; both Cg12 states. The maxima come from ordinary six-coordinate Cu(II) examples, including its aqua ion |
| Ag capacity | Increased 18 → 21. This admits four states previously excluded: Cg13 #t6 and Cg14 #v6 need 21; both Cg16 Ag(III) states need 20 |
| Au capacity | All 30 accepted states fit the current 18. The unselected six-coordinate Au(III) chelate is not used to enlarge the limit |
| ValenceTable | No target or aromatic entry change selected. Metal implicit H is zero except for the deliberate MH/MH2−/AuH2/AuH4− variants; supplied #h is supported despite the zero free-H enumeration bound |
| Other bounds | The registry grows its valence index for supplied #v values, so max_valence is not an admission limit. Literal charges, including the negative formal charges above, are accepted independently of the free-charge enumeration bounds. No bound change selected for these inputs |
| Spin evidence | Cg11/Cg16 require supplied #u2 for Counts to infer #n3. No geometry-dependent spin inference or resolver change is proposed |
| Coverage | 78 fixtures: Cu 29, Ag 18, Au 31. All 77 added states infer the metal lone-pair count; the four Cu/Ag d⁸ triplet states supply #u2. The extra fixture covers the existing Au³⁺ triplet. Explicit/mixed/implicit H variants and the selected ordinary/donor notations are covered; the two Cu MC states use (CuCH2SiMe3)4 and Stryker's reagent |

All four resolver/policy combinations agree on the new outputs. Independent checks
matched 88 metal sites to their expected registry states (352 checks) and checked
electron conservation for every resolved molecule. The 1,464 earlier snapshots are
unchanged; the full resolution suite now has 1,542 passing cases.
The graph library (3,594), graph properties (223), chemistry library (306), other
graph integration tests and doc tests pass. Clippy with warnings denied, formatting
and diff whitespace checks pass.

#### Copper-group primary sources

| Key | Source and access | Role |
| --- | --- | --- |
| WWS01 | Wang, Wang, Brown, Schwerdtfeger, Schröder and Schwarz, *J. Chem. Phys.* **114** (2001), 7388–7395. [DOI](https://doi.org/10.1063/1.1362289), [reviewed PDF](../materials/valence/Wang-2001-CuCl2-CuBr2.pdf) | Molecular CuCl2/CuBr2 and their anions; pp. 7388, 7391–7392 distinguish the d⁹ neutral and closed-shell anion descriptions. Cg05/Cg07 |
| Fr24 | Fritsch et al., *Chem. Eur. J.* (2024), e202400357. [DOI](https://doi.org/10.1002/chem.202400357), [reviewed PDF](../materials/valence/Fritsch-2024-Cu-hydride-Xray.pdf) | Stryker reagent's edge-bridging H and Cu(I), pp. 2, 4–5. Its reference 13 is Bennett et al., *Inorg. Chem.* **53** (2014), 2963–2967, [DOI](https://doi.org/10.1021/ic402736t), the earlier structural revision. Fr24 itself was downloaded and archived; Cg27's graph allocation is PW |
| AW69 | Allen and Warren, *Inorg. Chem.* **8** (1969), 1895–1901. [DOI](https://doi.org/10.1021/ic50079a018), [reviewed PDF](../materials/valence/Allen-Warren-1969-hexafluoro-anions.pdf) | CuF6³−/AgF6³− d⁸ ³A2g ground states: pp. 1895–1897, Table II. Cg11/Cg16 |
| S78 | Sorbe, Grannec, Portier and Hagenmuller, “Sur deux nouveaux composés fluorés de l'argent: Cs2AgF6 et Cs2Ga0.5Ag0.5F6,” *J. Fluorine Chem.* **11** (1978), 243–250. [DOI](https://doi.org/10.1016/S0022-1139(00)82444-6), [reviewed PDF](../materials/valence/Sorbe-1978-silver-fluorides.pdf) | Weak magnetism, proposed paired Ag(III)/Ag(V), and the Ga/Ag comparison: pp. 247–250. The mixed-valence interpretation remains a hypothesis, not a settled replacement for HW's Ag(IV). Cg17 |
| AW03H | Andrews and Wang, *J. Am. Chem. Soc.* **125** (2003), 11751–11760. [DOI](https://doi.org/10.1021/ja036307q), [reviewed PDF](../materials/valence/Andrews-Wang-2003-coinage-metal-hydrides.pdf) | Linear singlet MH2− and square-planar singlet AuH4−: pp. 11753–11757, Tables 2–3. Bent doublet AuH2: pp. 11757–11758. Isotopic IR supports the assignments; CuH4−/AgH4− are calculated but not observed. Cg05/Cg21/Cg23/Cg25/Cg26/Cg33 |
| Z23 | Zainab et al., “Critically Evaluated Atomic Data for Au IV Spectrum” (2023). [DOI](https://doi.org/10.3847/1538-4365/acd117), [reviewed PDF](../materials/valence/Zainab-2023-Au-IV.pdf) | Free Au³⁺ d⁸ triplet ground state and excited singlet levels, Table 2. Existing states Cg01/Cg32 |

### Zinc group

Source: [HW 103, Band II, chapter XXIII](<../materials/valence/Holleman-Wiberg - Anorganische Chemie - 103. Auflage - Band 2.pdf#page=251>),
printed pp. 1747–1781; chapter references pp. 1782–1783. Extraction: 2026-11-03.
Baseline: Zn 2, Cd 2 and Hg 3 states. The accepted additions are implemented:
53 metal states (Zn 19, Cd 17, Hg 17) and one O companion. Cn belongs to HW's
later transactinoid chapter.

#### Zinc-group work table

M = Zn or Cd; X = halide; R = organyl; D = neutral donor. States use concrete,
zero-valence defaults. Page references are to HW II. The current model counts two
valence electrons for Zn/Cd and twelve for Hg; these entries retain that convention.
**PW** marks our graph allocation: ordinary M–X/M–R bonds or anionic pair donors,
following the selected copper-group policy. Neutral donors remain dative.
The source establishes the compound/connectivity; the paired counts follow this allocation.

| ID | Atom descriptions / gap | Disposition | Examples and evidence |
| --- | --- | --- | --- |
| Zg01 | Existing Zn#n, Zn#c+2; Cd#n, Cd#c+2; Hg#n6, Hg#c+2#n5, Hg#c+#n5#v | All 7 retained; citations updated | Atomic configurations: HW 1747, 1751, 1769. Zn²⁺/Cd²⁺: 1753–1754; Hg²⁺/Hg2²⁺: 1767–1768, 1770–1772. The existing Hg(I) state also covers RHg⁺, 1781 |
| Zg02 | Zn#v2; Cd#v2; Hg#n5#v2 | Implemented: 3 | Ordinary bonds in MR2 and RMX, HW 1763–1764; HgX2, Hg(CN)2, HgR2 and RHgX, 1772, 1775, 1780–1781. Also X–Hg–Hg–X, 1771. Explicit-H MH2 reuses these states (SYB05, pp. 4710–4711), as do linear HZnZnH/HCdCdH (WA04Z, p. 11011) |
| Zg03 | M#c+2#t2; M#c+2#t3; M#c+2#t4; M#c+2#t5; M#c+2#t6 | Implemented: 10 | Two to six donor pairs. MR2/MR2D/MR2D2: HW 1763 (R− allocation: PW); MX3−/MX4²−: 1756; Zn(NH3)4²⁺, M(H2O)6²⁺ and M(NH3)6²⁺: 1753–1754, 1758. Five donors: Zn(acac)2(H2O), 1754, and CdCl5³−, 1756. These counts also cover the cited bridging-ligand and bioinorganic coordination without protein-specific types |
| Zg04 | Hg#c+2#n5#t2; Hg#c+2#n5#t3; Hg#c+2#n5#t4; Hg#c+2#n5#t5; Hg#c+2#n5#t6 | Implemented: 5; capacity increase approved and implemented | Hg(NH3)2²⁺/Hg(NH3)4²⁺, HW 1773; HgX3−/HgX4²−, 1768, 1772–1775; discrete HgCl5³−, 1773; Hg(H2O)6²⁺, 1768. Halide-donor allocation: PW. Six water ligands are explicitly stated, rather than inferred from hydrate stoichiometry |
| Zg05 | M#v2#t; M#v2#t2; Hg#n5#v2#t2 | Implemented: 5 | Two ordinary bonds plus neutral donors. R2MD/R2MD2, HW 1763; ZnCl2(NH3)2 and ZnCl2(H2O)2, 1753, 1755–1756; HgX2D2 (D = R2S/R3P), 1775, Fig. 23.10f. Hg2Cl6²− also uses Hg#n5#v2#t2 with two ordinary terminal bonds and two bridging Cl− donations per Hg (PW), 1772, Fig. 23.8a |
| Zg06 | M#c-#v3; M#c-2#v4; Hg#c-#n5#v3; Hg#c-2#n5#v4 | Implemented: 6 | Ordinary-bond descriptions of MX3−/MX4²−, HW 1756; MR3−/MR4²−, 1763; HgX3−/HgX4²− and Hg(CN)4²−, 1768, 1772–1775. Charges are formal graph allocations (PW); donor alternatives reuse Zg03/Zg04 |
| Zg07 | Cd#c-3#v5; Hg#c-3#n5#v5 | Implemented: 2 | Ordinary-bond descriptions of discrete CdCl5³− and HgCl5³−, HW 1756, 1773. Donor alternatives reuse Zg03/Zg04 |
| Zg08 | Zn#c-#v3#t | Implemented: 1 | Zn(OH)3(H2O)−, HW 1758: three ordinary Zn–OH bonds and one water donor. The all-donor graph reuses Zn#c+2#t4 (PW) |
| Zg09 | Zn#c+#v; Cd#c+#v | Implemented: 2 | Zn2²⁺/Cd2²⁺ with an ordinary M–M bond, HW 1762; Cd2[AlCl4]2 is isolated. Hg2²⁺ already has the corresponding #n5 state. Ordinary X–M–M–X uses Zg02 |
| Zg10 | Hg#c+#n5#v#t; Hg#c+#n5#v#t2 | Implemented: 2 | RHg(H2O)⁺/RHgD⁺ and RHg(bipy)⁺, HW 1781: one ordinary Hg–C bond and one/two neutral donor pairs. The first state also covers H2O→Hg–Hg←OH2²⁺, 1771, and the halide-donor description of Hg2X2 (PW) |
| Zg11 | Zn#c+#v#t2; Zn#c+#v#t3; Cd#c+#v#t3 | Implemented: 3 | [MeZnNPh2]2 and [RZnX]4/[RZnOR]4, HW 1764, Fig. 23.6d,f: ordinary Zn–R plus two/three anionic bridge donors (PW). ZnCp2/CdCp2 have one σ-Cp and one η⁵-Cp, 1763, 1781; η⁵-Cp− donates three pairs (PW), reusing #v#t3. Zn2Cp*2 with its ordinary Zn–Zn bond also reuses the Zn state, 1763, Fig. 23.5c. Existing N/halide companions suffice; the alkoxide companion is Zg24. Cp donor-side multidonor accounting remains separate, as in Tr28/Ae12 |
| Zg12 | Hg#c+#n5#v#t3 | Implemented: 1 | Hg2X4D2 and terminal sites of Hg3X6D2, HW 1775, Fig. 23.10g,h: one ordinary terminal Hg–X bond, two bridging X− donations and one neutral donor (PW). The central Hg in Hg3X6D2 reuses Hg#c+2#n5#t4 |
| Zg13 | Hg3²⁺/Hg4²⁺ | Reuse; no further type proposed | HW 1778 gives (nearly) linear Hg chains. Ordinary single bonds, terminal Hg#c+#n5#v and internal Hg#n5#v2 conserve the charge/electrons (PW). No cluster-size-dependent MC states |
| Zg14 | ZnN2⁴−; HgO2²−, HgS2²−, HgS4⁶− | Reuse; no further type proposed | HW 1761 draws N=Zn=N; an ordinary-double-bond graph uses Zn#c-2#v4 and two existing N#c-#n2#v2 states (PW). HW 1776–1777 gives isolated linear O–Hg–O/S–Hg–S and tetrahedral HgS4 units. Ordinary single bonds use existing terminal O−/S− states and Zg02/Zg06 |
| Zg15 | HgO/HgS chains; [HgNH2]⁺/[HgOH]⁺ chains; [Hg2N]⁺ network | Reuse; no further type proposed | HW 1773, 1776–1777 explicitly gives ordinary skeletal Hg–O/Hg–N single bonds. Hg#n5#v2 and existing O/N/S states suffice. Secondary interchain contacts need not be counted as skeletal bonds |
| Zg16 | Zn(BH4)2 | No additional type proposed | HW 1755. The separated-ion description uses Zn²⁺ and BH4−. A specific coordinated borohydride graph would need its own structural evidence; the formula alone does not select a Zn MC state |
| Zg17 | M#v#h; M#h2; Hg#n5#v#h; Hg#n5#h2 | Implemented: 6 | ZnH2/CdH2/HgH2, HW 1755. SYB05, pp. 4710–4711, resolves their rotational spectra, establishes linear closed-shell ¹Σg⁺ molecules and explicitly assigns strong covalent bonds. Mixed/implicit-H variants of Zg02; M#v#h also covers linear HZnZnH/HCdCdH, WA04Z p. 11011 |
| Zg18 | Zn#c+#v#t#m0 | Implemented: 1; bridge allocation PW | [RZnH·py]3, R = Et or Ph: HW 1764, Fig. 23.6e; KBK80, pp. 3–4, Fig. 1. Cryoscopy establishes trimers and IR supports bridging H; the six-membered Zn3H3 ring is the authors' structural proposal, not a diffraction result. One ordinary Zn–R bond, one pyridine donation and two local Zn–H−–Zn bridges per Zn, each [0,2,0], reuse H#c-#m2. The H− pair cannot supply two independent ordinary/dative bonds; no whole-ring MC state |
| Zg19 | LiZnH3, Li2ZnH4, Li3ZnH5 | Evidence reviewed; no added states proposed | HW 1755 lists solid products; II 2069, Table 32.1 footnote d lists ZnH3−/ZnH4²−/ZnH5³−. W77, pp. 36, 41–42, explicitly leaves their structures unresolved and reports unsuccessful attempts to repeat the Li3ZnH5 preparation. No isolated mononuclear unit or bridge pattern is established; no extra hydride/MC states proposed from stoichiometry alone |
| Zg20 | Hg#n4#v4; Hg#c+4#n4#t4 | Not selected: Hg oxidation states above II | HW 1770, Table 23.4 footnote b lists a prediction. WRK07, pp. 8371–8374, assigns HgF4 from matrix IR and calculations as square-planar, low-spin d⁸. Ordinary Hg–F bonds or four F− donors (PW) give these two descriptions. Ro08, pp. 4594, 4599, 4603–4604, found no convincing HgF4 in Ar at 10–14 K; it explicitly regards this as broadly consistent with WRK07's principally Ne evidence at 4 K. Both descriptions are excluded from the default selection |
| Zg21 | Zn#c+2#t8; Cd#c+2#t7; Hg#c+2#n5#t8 | Not proposed for the default registry | Zn(NO3)4²−, Cd(quinoline)2(NO3)2(H2O), HW 1754; eight-coordinate Hg(NO2)4²−/Hg(NO3)4²−, 1768, 1777. Extracted, but these particular higher-coordinate complexes do not presently justify additional default rows. Hg #t8 would require capacity 26 |
| Zg22 | Additional mixed ordinary/donor partitions | No additional types proposed | Zn(acac)2(H2O), HgCl2(terpy) and the higher-coordinate chelates in HW 1754, 1768 are representable with the donor counts above. Their existence does not require adding every mixed partition of the same coordination environment |
| Zg23 | Extended halide/chalcogenide/nitride solids, amalgams and mixed-metal Hg clusters | No additional types proposed | HW 1755–1761, 1770, 1778–1779. Ionic descriptions and selected ordinary/donor states cover the simple salts and motifs. No formal-ion ladder from metallic partial charges, lattice coordination-number ladder, or separate types for the Rh/Os–Hg clusters |
| Zg24 | O#c-#v#d3 | Implemented: 1 companion | μ₃-Alkoxide in [RZnOR]4, HW 1764, Fig. 23.6f: one ordinary O–R bond and three O−→Zn donations (PW). Supports the alkoxide-bridge family, independent of the organyl substituents. This is distinct from the deferred neutral-water/ether #v2#d2 states; no implicit-OH counterpart proposed here |
| Zg25 | M#u#v; M#u#h; Hg#n5#u#v; Hg#n5#u#h | Implemented: 6 | ZnH/CdH: WA04Z, pp. 11010–11011, Tables 5–6, identifies the ²Σ radicals spectroscopically. HgH: SYB05, pp. 4710–4711; Mo01, Table III, explicitly gives the ²Σ⁺ ground state and compares calculated and experimental constants. Ordinary single bond, explicit/implicit H; no additional metal-cluster states |
| Zg26 | Solid ZnH2/CdH2; small (MH2)n aggregates | Evidence reviewed; no further type proposed | WA04Z, p. 11012, supports hydrogen-bridged solid polymers but does not establish a specific network topology. The observed molecular dimers/trimers are assigned as van der Waals aggregates; no evidence was found for the calculated D2h dibridged dimers. Retain existing ionic salt descriptions and Zg02/Zg17 monomer states; no inferred polymer motif or extra MC ladder |

#### Zinc-group Counts comparison and coverage

For these nonaromatic states, the implemented arithmetic is
`E = ve − charge + h + v + 2t` and
`2n + u = ve − charge − h − v − 2d − m`, with `ve = 2` for Zn/Cd and `12` for Hg.
All 53 metal additions are distinct from the existing states and satisfy the
electron balance. The O companion uses `ve = 6` and requires capacity 8, within its
current limit. Every addition has inferred conformance coverage.

| Constraint | Implemented scope |
| --- | --- |
| Zn/Cd capacity | Retained at 18; the added states require at most 12 |
| Hg capacity | Approved and implemented: 18 → 22 for the explicitly six-coordinate aqua ion, counting ten d electrons and six donated pairs. Hg#c+2#n5#t5 and Hg#c-3#n5#v5 require 20; Hg#c+2#n5#t6 requires 22. All other selected states fit 18 |
| ValenceTable | No direct Zn/Cd/Hg entries or target/aromatic additions. Explicit-H and nonhydride fixtures default implicit H to zero; selected implicit-H variants supply #h/#h2. Zg18 needs capacity 4. All additions fit the current capacities |
| Other bounds | Unchanged. Supplied valences/charges are not restricted by the free-enumeration bounds; the registry valence index grows for the actual entries |
| Coverage | All 54 additions are inferred by independent resolution fixtures: 52 new inputs and four ZnMe2 inputs revised to infer Zn's occupation. Ordinary/donor and explicit/implicit-H variants are covered where selected. The seven baseline atom/ion cases are retained |

All four resolver/policy combinations agree for each new or revised input. The
expected atom states and valence invariants were checked independently for all
four outputs. The 1,542 earlier snapshots are unchanged; the full resolution suite
now has 1,594 passing cases. General Cp donor-side multidonor accounting remains
separate work, as recorded in Zg11.
The graph library (3,594), graph properties (223), chemistry library (306), other
graph integration tests and doc tests pass, as do Clippy with warnings denied,
formatting and diff whitespace checks.

#### Zinc-group primary-source checks

The five supplied PDFs and W77 are reviewed and archived. Mo01 was obtained as an
author manuscript to verify HgH's ground-state assignment. No download remains
necessary for the selected scope. Preserved sources are in `materials/valence/`;
`scratch/hw-zinc-review/` holds intermediates only.

| Key | Source | Evidence used |
| --- | --- | --- |
| SYB05 | Shayesteh, Yu and Bernath, “Gaseous HgH2, CdH2, and ZnH2,” *Chem. Eur. J.* **11** (2005), 4709–4712. [DOI](https://doi.org/10.1002/chem.200500332), [PDF](../materials/valence/Shayesteh-2005-group12-dihydrides.pdf) | Rotational spectra establish linear, closed-shell dihydrides; explicit covalent-bond discussion, pp. 4710–4711. Zg02/Zg17 |
| WA04Z | Wang and Andrews, “Infrared Spectra of Zn and Cd Hydride Molecules and Solids,” *J. Phys. Chem. A* **108** (2004), 11006–11013. [DOI](https://doi.org/10.1021/jp046414m), [PDF](../materials/valence/Wang-Andrews-2004-Zn-Cd-hydrides.pdf) | ZnH/CdH doublets and linear HMMH: pp. 11010–11011, Tables 5–6; molecular aggregates versus bridged solids: p. 11012. Zg02/Zg17/Zg25/Zg26 |
| KBK80 | de Koning, Boersma and van der Kerk, “Synthesis, characterization and properties of some organozinc hydride complexes,” *J. Organomet. Chem.* **195** (1980), 1–12. [DOI](https://doi.org/10.1016/S0022-328X(00)93381-3), [PDF](../materials/valence/de-Koning-1980-organozinc-hydrides.pdf) | Cryoscopy, NMR and IR support the trimer/bridge proposal; pp. 3–4, Fig. 1. The [0,2,0] allocation is PW. Zg18 |
| WRK07 | Wang, Andrews, Riedel and Kaupp, “Mercury Is a Transition Metal: The First Experimental Evidence for HgF4,” *Angew. Chem. Int. Ed.* **46** (2007), 8371–8375. [DOI](https://doi.org/10.1002/anie.200703710), [PDF](../materials/valence/Wang-2007-HgF4.pdf) | IR bands at 703 cm⁻¹ in Ne and 682 cm⁻¹ in Ar; square-planar low-spin d⁸ assignment, pp. 8371–8374. Zg20 |
| Ro08 | Rooms, Wilson, Harvey, Bridgeman and Young, “Mercury–fluorine interactions: a matrix isolation investigation of Hg⋯F2, HgF2 and HgF4 in argon matrices,” *Phys. Chem. Chem. Phys.* **10** (2008), 4594–4605. [DOI](https://doi.org/10.1039/B805608K), [PDF](../materials/valence/Rooms-2008-mercury-fluorides.pdf) | Negative Ar result and explicit comparison with WRK07's lower-temperature/Ne experiment; pp. 4594, 4599, 4603–4604. Zg20 |
| W77 | Watkins, *Preparations, Solution Composition, and Reactions of Complex Metal Hydrides and Ate Complexes of Zinc, Aluminum, and Copper*, PhD thesis, Georgia Institute of Technology (1977). [Repository PDF](https://repository.gatech.edu/bitstreams/e81805c5-5e1f-4d9b-88ef-891aad8840c5/download), [archived PDF](../materials/valence/Watkins-1977-complex-metal-hydrides-thesis.pdf) | Reviewed Part I, pp. 35–42: distinct complex-hydride products but unresolved structures; Li3ZnH5 preparation not reproduced. Zg19 |
| Mo01 | Mosyagin, Titov, Eliav and Kaldor, “GRECP/RCC calculation of the spectroscopic constants for the HgH molecule and its cation” (2001), [author manuscript](https://arxiv.org/abs/physics/0101047v4), [PDF](../materials/valence/Mosyagin-2001-HgH.pdf) | Table III explicitly identifies ground-state HgH as (σ²σ¹) ²Σ⁺ and compares calculations with experimental spectroscopy. Zg25 |

### Scandium group

Source: [HW 103, Band II, chapter XXIV](<../materials/valence/Holleman-Wiberg - Anorganische Chemie - 103. Auflage - Band 2.pdf#page=288>),
printed pp. 1784–1791; references p. 1792. Extraction: 2026-11-03.
Baseline: two states each for Sc, Y, La and Ac, all retained. The 79 accepted
additions (Sc 28, Y 25, La 26) are implemented and covered by 82 new fixtures.
The four supplied primary papers are reviewed and archived. The selected scope is
complete; lanthanoids and actinoids have later HW chapters.

#### Scandium-group work table

M = Sc, Y or La unless specified; R = organyl. States use concrete, zero-valence
defaults. Page references are to HW II. Ordinary M–R/M–N bonds and anionic pair
donors follow the selected transition-metal graph convention (PW); neutral donors
remain dative. Their counts describe the supplied graph, not measured populations.

| ID | Atom descriptions / gap | Disposition | Examples and evidence |
| --- | --- | --- | --- |
| Sg01 | M#n#u and M#c+3, including Ac | All 8 retained; citations checked | HW 1788–1790: free atoms s²d¹, trivalent ions d⁰. Y's atomic comment is corrected to s²d¹ |
| Sg02 | Sc#c+3#t6; Y#c+3#t8; La#c+3#t9 | Implemented: 3 | Sc(H2O)6³⁺, Y(H2O)8³⁺ and La(H2O)9³⁺, 1788, 1790. These are explicit aqua coordination spheres, not hydrate stoichiometries. The Y charge printed in the equation on p. 1790 conflicts with the explicit Y(III)/Y(H2O)8³⁺ description on p. 1788; use +3 |
| Sg03 | Y#c+3#t6; La#c+3#t6; Y#c+3#t7; La#c+3#t8 | Implemented: 4 | Six-coordinate M(NCS)6³−, seven-coordinate Y(acac)3(H2O), eight-coordinate La(bipyO2)4³⁺, 1788. These donor counts also serve other ligands with the same local description. Sc's six-pair state is Sg02 |
| Sg04 | M#v3 | Implemented: 3 | MR3 and M{N(SiMe3)2}3, 1791: three ordinary M–C or M–N bonds. Alkyl examples include CH2tBu, CH2SiMe3 and CH(SiMe3)2; HW explicitly describes the pyramidal tris-amides. Does not assert that all phenyl or solid halide examples are monomers |
| Sg05 | M#c+3#t3 | Implemented: 3 | Donor alternatives for Sg04, using three R− or amido donors (PW). Existing C/N companion states suffice; no new ligand state |
| Sg06 | M#v3#t2; M#c+3#t5 | Implemented: 6 | MR3(THF)2, 1788, 1791: three ordinary M–R bonds plus two neutral donors, or three R− donors plus two neutral donors (PW). HW gives trigonal-bipyramidal coordination |
| Sg07 | Sc#c-3#v6 | Implemented: 1 | Ordinary-bond alternative for ScF6³−/ScCl6³−, 1787, 1790. HW explicitly gives octahedral ScF6³−. The halide-donor alternative reuses Sc#c+3#t6; existing F/Cl states suffice |
| Sg08a | La#c-3#v6 | Implemented: 1 | LaCl6³−: M02G, pp. m186–m187, establishes discrete, slightly distorted octahedral anions by diffraction. Six ordinary La–Cl bonds (PW); the chloride-donor graph reuses La#c+3#t6 from Sg03 |
| Sg08b | YF4−, LaF4−, LaF6³−; LiLaPh4 | Connectivity unresolved; no further types selected from formulae | HW 1790–1791. The hydride additions independently supply Y/La#c-#v4 and Y/La#c+3#t4; Sg08a/Sg03 supply La's six-bond/six-donor counts. These do not establish the fluoride or phenyl structures. No formula-only fixture or additional acceptor count |
| Sg09a | M#n#v; M#n#h | Implemented: 6 | ScH/YH/LaH: WCA02, pp. 9216–9217, Table 4, assigns ¹Σ⁺ ground states using gas-phase spectroscopy and calculations. Ordinary single bond and explicit/implicit H; one remaining electron pair |
| Sg09b | M#u#v2; M#u#v#h; M#u#h2, M = Y or La | Implemented: 6 | Observed YH2/LaH2: WCA02, pp. 9218–9220, Table 5, assigns bent ²A1 molecules. Bare ScH2 was calculated but not observed; no Sc counterparts selected from that calculation alone |
| Sg09c | M#c+#v2; M#c+#v#h; M#c+#h2 | Implemented: 9 | ScH2⁺/YH2⁺/LaH2⁺: WCA02, pp. 9218–9221, Table 6. Matrix IR, isotopic substitution and charge-selective experiments support the molecular cations; calculated bent ¹A1 ground states |
| Sg09d | M#v2#h; M#v#h2; M#h3 | Implemented: 9; explicit-H M#v3 reused from Sg04 | ScH3/YH3/LaH3: WCA02, p. 9221, Table 7. Observed molecular trihydrides with singlet assignments; Sc's calculated planar/pyramidal distinction does not change these counts |
| Sg09e | M#c-#v4; M#c-#v3#h; M#c-#v2#h2; M#c-#v#h3; M#c-#h4 | Implemented: 15 | ScH4−/YH4−/LaH4−: WCA02, pp. 9222–9223, Table 8. Isotope-resolved spectra and electron-capture experiments support the anions; tetrahedral ¹A1 structures. All explicit/mixed/implicit-H variants |
| Sg09f | Sc#c+#u#v; Sc#c+#u#h | Implemented: 2 | ScH⁺: WCA02, pp. 9216–9217, Table 4. Matrix assignment compared with gas-phase work and the calculated ²Δ ground state. No extrapolated YH⁺/LaH⁺ states |
| Sg09g | M#c+3#t2; M#c+3#t4 | Implemented: 6 | Hydride-donor alternatives for the closed-shell MH2⁺ and MH4− ions, using H−→M³⁺ (PW). The molecular/electronic-state evidence is WCA02; the allocation is our existing donor convention. MH3's donor graph reuses Sg05 |
| Sg09h | (H2)MH2, (H2)2MH2 and related H2 adducts | Deferred; no new MC states | WCA02, pp. 9221–9224, distinguishes side-bound H2 from ordinary hydrides. H2-donor accounting remains separate, as in Cg33 |
| Sg10a | Sc#c+2#v#t6; Sc#c+3#t7 | Implemented: 2 | Cp*2ScX, X = Cl, Me, CH2Ph or Ph: TB84, pp. 2–4, establishes monomers and d⁰/14-electron Sc(III). Two η⁵-Cp− ligands donate three pairs each; X is ordinary or anionic donor (PW). Q24's crystallographic Cp-substituted Sc(III) iodide and phenyl structures also support these counts |
| Sg10b | Sc#c+2#v#t7; Sc#c+3#t8 | Implemented: 2 | Cp*2ScMe(THF), TB84 pp. 3–4: one additional neutral donor. Sc#c+3#t8 also covers the parent [Cp2Sc(μ-Cl)]2 cited on pp. 2–3, with two bridging Cl− donors per Sc (PW); Y/La eight-pair states already occur in Sg02/Sg03. No additional bridge-specific type |
| Sg10c | Sc#c+2#u#t6 | Implemented: 1 | Isolated (C5H2tBu3)2Sc: Q24, pp. 3281, 3285–3287, establishes Sc(II), one unpaired electron and two Cp rings by diffraction, magnetism, EPR and calculations. Existing η⁵-Cp− donation convention (PW); 13 electrons at Sc |
| Sg10d | Sc#c+#u#v#t6; Sc#c+2#u#t7 | Supported; not selected for the default registry | Q24, pp. 3283–3285, establishes [Cp″2ScI]− with d¹ Sc(II). Ordinary-I or I−-donor graphs require these two extra states (PW). Recorded separately from the implemented six-pair scandocene state; no default expansion solely for this more specialized anion |
| Sg10e | [Cp*2ScH]x; Cp*2ScH(THF); other substituted Y/La Cp derivatives | No further types selected in this pass | TB84, pp. 3–4, 10–11, explicitly leaves the base-free hydride's molecularity uncertain and calls its THF adduct only very likely monomeric. No hydride-specific #h or MC addition inferred from these formulations. Other monomeric Y/La derivatives need their own structural evidence before further ordinary/donor combinations |
| Sg11 | ScCp3/LaCp3 chains; zero-valent M(η6-arene)2 (M = Sc, Y) | No additional default types selected in this pass | HW 1791 gives mixed-hapticity Cp bridges and bulky bis-arene examples. No per-chain acceptor ladder or neutral-metal spin assignment from oxidation state alone. Ordinary molecular metallocene cases remain in Sg10 |
| Sg12 | Hydrated ScOH²⁺, Sc(OH)2⁺ and hydroxo-bridged Sc oligomers | No additional types selected from formulae | HW 1790 lists compositions without complete hydration/bridge graphs. No additional mixed ordinary/donor patterns are inferred from these formulae |
| Sg13 | Sc/Y/La borohydrides | Existing ionic description; no new MC states selected | HW 1789 cross-refers to I 1248 and II 2067–2069. Those passages list M(BH4)3, whereas p. 1789 prints M(BH3)3. M³⁺ plus three BH4− uses existing states; the cross-references do not establish a particular Sc/Y/La bridge graph |
| Sg14 | Metallic MH2, YH3/LaH3, MI2/MO/MS and low-valent cluster halides | No extra states selected from bulk stoichiometry | HW 1787, 1789–1790 distinguishes conduction electrons from formal M²⁺ and describes extended cluster structures. No isolated M²⁺ ladder, arbitrary whole-lattice MC states or fictitious molecular MH2 graphs |
| Sg15 | Oxides, hydroxides, nitrides, carbides, borides, oxosalts and nitridoscandates | Existing ionic states; no new lattice-specific rows selected | HW 1790–1791. Li3ScN2 contains a network of corner-linked ScN4 tetrahedra, not an isolated ScN2³− molecule; other lattice coordination numbers likewise do not establish discrete molecular units |
| Sg16 | Nine-coordinate Sc nitrate; ten-coordinate Y/La and twelve-coordinate La complexes | No additional default rows selected | HW 1788: Sc(NO3)5²−, Y(NO3)5²−, La(EDTA)(H2O)4 and hydrated lanthanum sulfate. These particular complexes do not yet justify extending the default counts; the common aqua ions are Sg02 |
| Sg17 | Ac compounds | Existing atom/Ac³⁺ states retained; citations checked | HW 1789–1791 lists halides and other binary salts but supplies no discrete molecular or aqua coordination sphere for further Ac types. No extrapolation of the Sc/Y/La donor series |

#### Scandium-group Counts comparison and coverage

The group-3 review used three valence electrons and capacity 18 for all four elements.
The later [chapter-XXXII selection](#transition-metal-classes-chapter-xxxii) raises
Sc/Y/La capacity to 20 for bis-COT anions; Ac retains 18.
For the 79 additions, `E = 3 − charge + h + v + 2t` and
`2n + u = 3 − charge − h − v − 2d`; none has an MC contribution.
All are distinct from the existing entries. That batch required capacities from 3 to
18 and did not change the limits. The hydride singlet/doublet occupations and
Sc(II) doublet are taken from the reviewed sources, rather than from free-ion spin.
There are no direct Sc/Y/La/Ac ValenceTable entries. Explicit-H and nonhydride
graphs default implicit H to zero; deliberate implicit-H variants supply #h.
No target or aromatic table addition was needed.

All 79 added states have inferred occupation coverage in 82 independent fixtures:
27 donor–acceptor complexes, 50 hydride inputs, two molecular ions and three
ordinary trialkyls. The three extra fixtures cover explicit-H MH3 alongside the
trialkyls that share M#v3. All eight earlier atom/ion fixtures are retained.
Cp fixtures use the existing haptic convention; general donor-side multidonor
accounting remains separate, as in Tr28/Ae12/Zg11.

All four resolver/policy outputs agree for each new input. Expected metal states
and the valence invariants were checked independently for all four outputs,
including the retained atom/ion cases. All 1,676 resolution cases pass; the 1,594
earlier snapshots are unchanged. Graph library tests (3,594), properties (223),
the remaining graph integration tests and doc tests pass, as do Clippy with
warnings denied, formatting and diff whitespace checks.

#### Scandium-group primary-source checks

All four supplied PDFs are reviewed and preserved in `materials/valence/`.
`scratch/hw-scandium-review/` contains intermediates only. No additional paper is
needed for the implemented selections above; unresolved structures and unselected extensions
remain explicit in the work table.

| Key | Source | Evidence used |
| --- | --- | --- |
| WCA02 | Wang, Chertihin and Andrews, “Matrix Infrared Spectra and DFT Calculations of the Reactive MHx (x = 1, 2, and 3), (H2)MH2, MH2⁺, and MH4− (M = Sc, Y, and La) Species,” *J. Phys. Chem. A* **106** (2002), 9213–9225. [DOI](https://doi.org/10.1021/jp026166z), [PDF](../materials/valence/Wang-Chertihin-Andrews-2002-group3-hydrides.pdf) | Tables 4–8 and pp. 9216–9224 distinguish molecular hydrides, ground-state occupations and H2 adducts. Bare ScH2 is explicitly not observed. Sg09 |
| TB84 | Thompson and Bercaw, “Some aspects of the chemistry of alkyl and hydride derivatives of permethylscandocene,” *Pure Appl. Chem.* **56** (1984), 1–11. [DOI](https://doi.org/10.1351/pac198456010001), [PDF](../materials/valence/Thompson-Bercaw-1984-permethylscandocene.pdf) | Monomeric alkyl/chloro/aryl derivatives and THF adducts, pp. 2–4; no ground-state agostic methyl interaction, pp. 4–5, 11. The conclusion explicitly qualifies the hydride molecularity. Sg10a/b/e |
| M02G | Matsumoto et al., “Tris(1-ethyl-3-methylimidazolium) hexachlorolanthanate,” *Acta Cryst. C* **58** (2002), m186–m187. [DOI](https://doi.org/10.1107/S0108270101020601), [PDF](../materials/valence/Matsumoto-2002-LaCl6.pdf) | Single-crystal structure establishes discrete octahedral LaCl6³−. Sg08a |
| Q24 | Queen et al., “Synthesis of Crystallographically Characterizable Bis(cyclopentadienyl) Sc(II) Complexes: (C5H2tBu3)2Sc and {[C5H3(SiMe3)2]2ScI}1−,” *J. Am. Chem. Soc.* **146** (2024), 3279–3292. [DOI](https://doi.org/10.1021/jacs.3c11922), [PDF](../materials/valence/Queen-2024-Sc-II-scandocene.pdf) | Diffraction, magnetism, EPR and calculations support d¹ Sc(II) in the isolated scandocene and iodido anion, pp. 3281–3287. The putative Cp*2Sc solution product is not structurally established and its calculated hyperfine coupling disagrees with experiment. Sg10a/c/d |

The earlier MH4− report, Wang and Andrews, *J. Am. Chem. Soc.* **124** (2002),
7610–7613, [DOI](https://doi.org/10.1021/ja020112l), was not obtained. WCA02 supplies
the needed structures, electronic states and experimental assignments; it is not
an outstanding download requirement for this batch.

### Titanium group

Source: [HW 103, Band II, chapter XXV](<../materials/valence/Holleman-Wiberg - Anorganische Chemie - 103. Auflage - Band 2.pdf#page=297>),
printed pp. 1793–1817; references p. 1818. Extraction: 2026-11-03.
Baseline: Ti 3, Zr 2, Hf 2 and Rf 1 states. All eight earlier states are retained
with checked HW citations. Rf chemistry is
referred to chapter XXXVII on p. 1793.

#### Titanium-group work table

M = Ti, Zr or Hf unless restricted; R = organyl; D = neutral pair donor.
States use concrete, zero-valence defaults. Page references are to HW II.
Ordinary M–X/M–R bonds, anionic donor alternatives and Cp− donation follow PW's
selected transition-metal convention. All **118 selected additions are implemented**
(Ti 47, Zr 44, Hf 27), including both ZrH2 occupations (Tg18g), the two
Schwartz's-reagent states (Tg17) and two ZrCp4 states admitted by the approved
capacity 20. Counts exclude states already listed above in the table.
Unselected rows are not included.

| ID | Atom descriptions / gap | Disposition | Examples and evidence |
| --- | --- | --- | --- |
| Tg01 | Ti/Zr/Hf/Rf#n#u2; Ti#c+3#u; Ti/Zr/Hf#c+4 | Retained: all 8; citation audit implemented | HW 1538, Table 19.1, gives the free-atom s²d² triplets, with Rf explicitly calculated (also 2358). HW 1797, 1811 gives the ionic states. Free-atom/ion occupations are separate from complex occupations |
| Tg02 | M#v4; M#c+4#t4 | Implemented: 6 | Tetrahedral molecular MX4, 1799, 1813, and MR4, 1807, 1817. Four ordinary bonds or four anionic donors (PW). Gas-phase molecular halides are distinguished from polymeric solids |
| Tg03 | Ti#v4#t; Ti#v4#t2; Ti#c+4#t5; Ti#c+4#t6 | Implemented: 4 | TiMe4·OEt2, TiMe4·2OEt2, 1806; TiCl4·2D, 1799. Ordinary alkyl/halide bonds plus neutral donation, or all-donor alternatives. TiOCl2(NMe3)2, 1797, also uses Ti#v4#t2 with localized Ti=O |
| Tg04 | M#c-2#v6; M#c+4#t6 | Implemented: 5 new; Ti donor state reused | Octahedral MX6²−, 1799, 1811, 1814. Ordinary halide bonds or six X− donors; the Zr/Hf hexamethyl anions on 1817 reuse these counts |
| Tg05 | Ti#c-#v5; Zr#c-#v5; Zr#c+4#t5 | Implemented: 3; Ti donor state reused | TiMe5−, 1806–1807, Fig. 25.4a; discrete ZrCl5−, 1811. Ordinary or anionic-donor graphs; no extrapolated Hf five-coordinate entry |
| Tg06 | Zr#c-3#v7; Zr#c-4#v8; Zr/Hf#c+4#t7; Zr/Hf#c+4#t8 | Implemented: 6 | ZrF7³− and ZrF8⁴−, 1811, 1814. Seven/eight-donor Zr/Hf chelates MCl(acac)3, M(acac)4 and M(ox)4⁴−, 1814–1815; Zr4(OH)8(H2O)16⁸⁺, 1815–1816, reuses Zr's eight-donor state. No hydration count inferred from salt stoichiometry |
| Tg07 | Ti#c+4#t7; Ti#c+4#t8 | Implemented: 2 | Ti(O2)F5³− and TiCl(S2CNMe2)3, 1797; Ti(NO3)4, 1799, Fig. 25.2c. Seven/eight pair donors, with peroxide O–O retained; no MC description needed |
| Tg08 | M#v4#t3; M#v4#t4; Zr/Hf#v4#t2 | Implemented: 8 | Mixed ordinary/donor chelates: TiCl(S2CNMe2)3 and Ti(NO3)4, 1797–1799; Zr/Hf MCl2(acac)2, MCl(acac)3 and M(acac)4, 1814–1815. One ordinary anionic-donor bond and one neutral donation per chelate (PW). These are local counts shared across the ligand families |
| Tg09 | Ti#u#v3; Ti#c+3#u#t3 | Implemented: 2 | Tris-organyls Ti{CH(SiMe3)2}3, 1807, and tris-amides, 1797. Ordinary bonds or anionic donors leave one Ti electron |
| Tg10 | Ti#c+3#u#t5; Ti#c+3#u#t6; Ti#u#v3#t2; Ti#u#v3#t3; Ti#c+#u#v2#t4; Ti#c+2#u#v#t5; Ti#c-3#u#v6 | Implemented: 7 | TiBr3(NMe3)2, 1797; TiCl3(THF)3/amine adducts, Ti(H2O)6³⁺, TiCl2(H2O)4⁺, TiCl(H2O)5²⁺ and TiX6³−, 1803–1804. HW explicitly assigns one unpaired electron to the Ti(III) complexes. All-donor six-coordinate alternatives share Ti#c+3#u#t6 |
| Tg11 | Ti#c+2#u2; Zr#c+3#u; Hf#c+3#u | Implemented: 3 | Ionic descriptions of the lower halides, 1804–1805, 1812–1814. CS79, p. 20, establishes free Ti²⁺ as 3d², ³F2; this does not assign the spin of every Ti(II) complex. The trivalent Zr/Hf ions have one valence electron |
| Tg12 | M#c+2#v2#t6 | Implemented: 3; M#c+4#t8 reused | Cp2MX2 and Cp2MR2, 1807, 1817: two η⁵-Cp− donors contribute six pairs; X/R is ordinary or anionic donor (PW). TiCp4 and HfCp4 have two η⁵ and two η¹ rings and reuse these counts |
| Tg13 | Ti#c+2#u#v#t6; Ti#c+3#u#t7; Ti#c+3#u#t8 | Implemented: 3 | Monomeric Cp*2TiX and dimeric [Cp2TiX]2, X = Cl, Br, I, 1807. Ordinary terminal X or X− donation; the dimer uses two bridging X− donors at each Ti. No extra MC type for halide bridges |
| Tg14 | Zr#c+3#v#t6; Zr#c+3#h#t6 | Implemented: 2; Zr#c+4#t7 reused | Cp2ZrR⁺ and Cp2ZrH⁺ in the polymerization scheme, 1817–1818: ordinary R/H, implicit-H variant, or anionic donor. Solvent/olefin adducts with one further donated pair reuse Zr#c+4#t8 in the all-donor description |
| Tg15 | Ti#c-2#n3#t6 | Implemented: 1 | Ti(CO)6²−, 1797, 1805 n. 5. HW explicitly identifies six-coordinate Ti(−II) as closed shell; six CO pair donors and three remaining pairs give 18 electrons |
| Tg16 | Ti#u#m3; Zr#m4; Hf#m4 | Unselected: no default additions | Ti(BH4)3, 1803, and Zr/Hf(BH4)4, 1812: HW specifies the B–H–M bridges; the neutral local-3c2e contribution allocation is PW. Aggregate #m3/#m4 does not distinguish the individual bridge patterns. These examples do not meet the common-class inclusion threshold for adding such broad MC states; extraction retained |
| Tg17 | Zr#c+3#v#t6#m0; Zr#c+4#t7#m0 | Implemented: 2 | Schwartz's reagent [Cp2ZrCl(μ-H)]2, HW 1817; J19, p. 1508, Fig. 1, establishes the dimer and two hydride bridges. Widespread synthetic use justifies inclusion. Existing H#c-#m2 supplies each local 3c2e bridge (allocation: PW); Cl is ordinary or a Cl− donor |
| Tg18a | Ti#u3#v; Ti#u3#h; Ti#c+#u3#t | Implemented: 3 | TiH: CA94, pp. 8325–8326, identifies the molecule and assigns a quartet ground state. Ordinary/implicit H or H− donation (PW); three unpaired Ti electrons |
| Tg18b | Hf#n#u#v; Hf#n#u#h; Hf#c+#n#u#t | Implemented: 3 | HfH: CA95's matrix band is tentative, but RB94, pp. 74–75, 78–79, independently observes HfH/HfD and assigns the lower state as X²Δ3/2. One remaining pair and one unpaired electron in the ordinary-H or H−-donor allocation (PW) |
| Tg18c | Ti/Zr#u#v3; #u#v2#h; #u#v#h2; #u#h3; #c+3#u#t3 | Implemented: 8 new; Ti#u#v3 and Ti#c+3#u#t3 reused | TiH3 and ZrH3: CA94, pp. 8325–8326, and CA95, pp. 15007–15008, identify the trihydrides by mixed-isotope spectra. Three ordinary bonds leave one metal electron; the all-hydride-donor alternative is PW |
| Tg18d | M#v4; M#v3#h; M#v2#h2; M#v#h3; M#h4; M#c+4#t4 | Implemented: 12 new H variants; Tg02 states reused | TiH4, ZrH4 and HfH4: CA94, pp. 8326–8327, and CA95, pp. 15008–15010, establish tetrahedral molecular tetrahydrides. Explicit/mixed/implicit H and the all-H−-donor alternative (PW). CA95 supplies the ZrH4/HfH4 evidence missing from HW 1812 |
| Tg18e | Ti#u2#v2; Ti#u2#v#h; Ti#u2#h2; Ti#c+2#u2#t2; Hf#n#v2; Hf#n#v#h; Hf#n#h2; Hf#c+2#n#t2 | Implemented: 8 | MCS96, pp. 872–873, Table 1: TiH2's lowest states are both triplets; the singlet is 32.6 kcal/mol higher at TZP CCSD. DB91, pp. 3981–3982, Tables II/V/VI: HfH2 is predominantly closed-shell ¹A1; the lowest triplet-derived level remains 0.40 eV higher after spin–orbit coupling. Ordinary/implicit H or H− donation (PW). MCS96, pp. 877–879, also suggests TiH2·H2 for the matrix bands; the bare-TiH2 occupation rests on its calculation |
| Tg18f | ZrH and HfH3 | Unselected: tentative identification | CA95, pp. 15008–15009, explicitly qualifies these assignments. The follow-up found no independent firm identification to replace those tentative bands; no entries proposed |
| Tg18g | Zr#n#v2; Zr#n#v#h; Zr#n#h2; Zr#c+2#n#t2; Zr#u2#v2; Zr#u2#v#h; Zr#u2#h2; Zr#c+2#u2#t2 | Implemented: 8; both occupations accepted | WDB90, p. 152, Table 2: MRSDCI places ¹A1 lowest, with ³B1/³A1 at +0.086/+0.145 eV; CASSCF reverses the ordering. The authors explicitly leave the ground-state identification uncertain. Both paired and unpaired descriptions are included, with ordinary/implicit H or H− donation (PW), without selecting a definitive ground state |
| Tg19a | Ti#c+2#u#h#t6 | Implemented: 1; explicit-H and all-donor states reuse Tg13 | Cp*2TiH, HW 1807. L95, pp. 125–127, 139–140: EPR and XRD with terminal H located/refined. Ordinary H uses Ti#c+2#u#v#t6; H− donation uses Ti#c+3#u#t7 (PW) |
| Tg19b | Zr/Hf#c+2#v#h#t6; Zr/Hf#c+2#h2#t6 | Implemented: 4; explicit-H and all-donor states reuse Tg12/Tg06 | Monomeric Cp*2MH2, C00, chapter 1, pp. 3–6. Mixed/implicit H variants of M#c+2#v2#t6; two H− donors use M#c+4#t8 (PW). Hf monomer structure established, but H positions were not located. These dihydrides suffice without assuming monomeric solid Schwartz reagent |
| Tg20a | Ti#u2#v2#t4; Ti#c+2#u2#t6 | Implemented: 2 | TiCl2(tmeda)2, HW 1805; Z20 establishes two Cl/four N coordination and S = 1 by magnetometry, HFEPR and MCD. Ordinary chlorides or six pair donors (PW) |
| Tg20b | Ti#c+2#u2#t6 | Reuse Tg20a; no extra state | Cp*2Ti, HW 1808. Bercaw's reported magnetic moment supports two unpaired electrons; the values are reproduced explicitly in [the 2025 titanocene study](https://doi.org/10.1039/D5DT00480B). The [1974 PDF](https://doi.org/10.1021/ja00823a012) is unavailable; the state has independent primary support from Z20 |
| Tg20c | Ti#c+2#n#t8 | Implemented: 1 | Cp2Ti(CO)2: HW 1808, Fig. 25.6h supplies the graph. M61, p. 1288, supports a monomer with terminal CO by IR and gives a closed-shell electron-count argument; M79, p. 3425, explicitly describes a diamagnetic monomer. Two Cp− and two CO donors leave one metal pair (PW). M61 reports no magnetic measurement |
| Tg21 | Zr#c+3#v#t9; Zr#c+4#t10 | Implemented: 2; capacity 20 approved and implemented | ZrCp4, 1817: three η⁵-Cp− rings plus one η¹ ring. Ordinary η¹ Zr–C or Cp− pair donation gives 20 electrons under PW. This bounds expansion admits the selected donor convention without changing ligand roles |
| Tg22a | Ti(CN)6³−; older Ti(CN)7⁴− and lower cyanide formulae | Hexacyanide reuses Tg10; other formulae unselected | E97 establishes octahedral C-bound Ti(CN)6³− and d¹ by XRD/EPR/magnetometry: Ti#c-3#u#v6 or Ti#c+3#u#t6 (PW). Its review distinguishes inconsistent older six-/seven-/eight-cyanide assignments. No new type inferred from Ti(CN)7⁴−, K4Ti(CN)4 or K5[M(CN)5] stoichiometry alone |
| Tg22b | Ti/Zr/Hf bipyridyl complexes | Reuse Ti#c+3#u#t6 and Zr/Hf#c+4#t6 | B13's primary abstract explicitly assigns TiIII d¹ and ZrIV/HfIV with ligand-centered reduction. No formal Ti(0)/Ti(−I) occupation ladder. Ligand-radical graph examples are unselected in this pass; no companion MC types proposed |
| Tg22c | Zr/Hf#c-2#n3#t6 | Implemented: 2; Ti state reused from Tg15 | Zr/Hf(CO)6²−: EC90, pp. 6022, 6024, Tables II/III, gives octahedral structures, a low-spin d⁶ rationale and closely similar ¹³C NMR spectra across Ti/Zr/Hf. These support the paired assignment: six CO donors and three metal pairs (PW), 18 electrons. No separate magnetic measurement is reported |
| Tg22d | Ti#n#u2#t6 | Implemented: 1 | Neutral Ti(CO)6, HW 1805 n. 5. D20, p. 10493, calculates a D3d triplet; one pair and two unpaired metal electrons with six CO donors (PW). This is a calculated occupation, not an experimental spin determination. The paper's Ti(CO)7 and Zr/Hf(CO)8 are outside this extracted batch |
| Tg22e | Ti/Zr/Hf#n2#t6 | Implemented: 3 | Bis(η⁶-arene) complexes, HW 1808, 1811, 1817. C87 explicitly reports diamagnetic Ti/Zr/Hf bis(1,3,5-tri-tert-butylbenzene) complexes in its primary abstract. Two neutral six-electron donors leave two metal pairs (PW); no new MC description. PDF unavailable for archiving |
| Tg22f | Zr/Hf#n2#t7 | Implemented: 2 | Zr/Hf(arene)2(PMe3), HW 1811. CG79, p. 127, reports isolated compounds and resolved ¹H NMR spectra with ordinary chemical shifts and P–H couplings. Paired occupation is inferred from these data; no magnetic measurement is reported. Two η⁶-arene donors plus PMe3 leave two metal pairs (PW), 18 electrons; the bent geometry is the authors' proposal |
| Tg22g | Remaining lower-valent Zr/Hf organyls | Unselected: no further default types | HW 1817 lists bulky substituted-Cp/alkyl derivatives. The located [Lappert et al. (1981)](https://doi.org/10.1039/DT9810000805) abstract concerns anionic Ti/Zr reduction products, not a general neutral Zr/Hf triorganyl family. It does not justify extrapolated neutral occupations |
| Tg22h | Zr#c+3#u#t9 | Implemented: 1 | Cp3Zr: L95, pp. 99–101, establishes a monomer with three η⁵ rings, Zr-coupled EPR and μeff = 1.64 BM. This contradicts HW 1817's statement that ZrCp3 was not obtainable. Three Cp− donors and one unpaired electron give 19 electrons (PW), within the approved capacity 20 |
| Tg23 | Alkoxide, halide and hydroxo aggregates with ordinary anionic donor bridges | Reuse implemented local donor counts; no extra scaffold types | Ti(OR)4 and [TiCl4(OPCl3)]2, 1799; hydrated Zr/Hf fluorides, 1814; Zr hydroxo tetramer and basic salts, 1815–1816. Existing O−/halide donor types provide the bridges. Ti(IV) aqueous formulae on 1801 do not establish an unhydrolyzed Ti(H2O)6⁴⁺ ion |
| Tg24 | Centered Zr6 clusters, extended low-valent solids, N2-bridged assemblies and unusual Cp rearrangement products | Extracted; no additional default types proposed | HW 1805–1808, 1813–1817. Bulk stoichiometry does not select a molecular graph; specialized assemblies do not justify an extra occupation/MC ladder in this pass. Ti(AlH4)4 on 1798 also lacks a complete bridge graph |

#### Titanium-group limits and coverage

All 118 implemented additions satisfy the nonbonding
budget and capacities: Ti 18, Zr/Hf 20. Tg21 requires 20; Tg22h requires 19;
the other accepted states require at most 18. ZrCp4 prompted the increase to 20.
D20 also observes both Zr(CO)8 and Hf(CO)8, each counting 20 with eight CO
pair donors under our convention; the approved bound is therefore 20 for both
Zr and Hf. Capacity and default-registry inclusion are separate decisions; the
octacarbonyl states remain outside this batch. No target-covalence or aromatic
ValenceTable addition is identified for graphs with zero or explicit implicit H.
All added states have inferred fixture coverage and pass full-graph valence invariants.

For the paired/unpaired ZrH2 alternatives with fixed H and donor counts,
AtomTyping's MostSaturated favors #n; Strict retains both if occupation is unspecified.
Supplied #u2 selects the unpaired alternative. Counts pairs the two remaining
electrons when both #n and #u are absent; a triplet input must supply occupation
evidence such as #u2. This follows the current resolver behavior, not an energy
ordering of the molecular states.

The 135 fixtures cover all 118 additions: eight ZrH2 inputs and 127 inputs for
the remaining 110 states. They cover ordinary/donor alternatives, explicit/mixed/
implicit hydrides, chelates, Cp/arene donors, carbonyls and both Schwartz's-reagent
descriptions. Occupation is inferred from supplied charge, H count and spin evidence;
ZrH2 uses #u0 or #u2 to distinguish the requested occupation. Both resolvers and
both policies agree for every input, with independently checked expected metal
states and valence invariants. Existing ligand types suffice.

Verification: all 3,594 graph unit tests, 223 property tests and 1,811 resolution
cases pass, as do the other graph integration tests and doctests. The 1,684 earlier
snapshots are unchanged. Existing multidonor accounting limitations remain outside
this chapter's scope.

#### Titanium-group primary-source checks

| Key | Source | Access and evidence |
| --- | --- | --- |
| CS79 | Corliss and Sugar, “Energy Levels of Titanium, Ti I through Ti XXII,” *J. Phys. Chem. Ref. Data* **8** (1979), 1–62. [DOI](https://doi.org/10.1063/1.555591), [PDF](../materials/valence/Corliss-Sugar-1979-Ti-energy-levels.pdf) | Reviewed and archived. Ti²⁺ ground configuration and term, p. 20; Tg11 |
| J19 | Jones et al., “Characterization of Reactive Organometallic Species via MicroED,” *ACS Cent. Sci.* **5** (2019), 1507–1513. [DOI](https://doi.org/10.1021/acscentsci.9b00403), [PDF](../materials/valence/Jones-2019-Schwartz-reagent-MicroED.pdf) | Reviewed and archived from the authors' Caltech copy. Dimeric Schwartz reagent and bridging-hydride locations, p. 1508, Fig. 1; Tg17 |
| CA94 | Chertihin and Andrews, “Reactions of Laser Ablated Ti Atoms with Hydrogen during Condensation in Excess Argon. Infrared Spectra of the TiH, TiH2, TiH3, and TiH4 Molecules,” *J. Am. Chem. Soc.* **116** (1994), 8322–8327. [DOI](https://doi.org/10.1021/ja00097a045), [PDF](../materials/valence/Chertihin-Andrews-1994-Ti-hydrides.pdf) | Supplied PDF reviewed and archived. TiH quartet, TiH3 doublet and tetrahedral TiH4; the TiH2 triplet assignment is qualified in the discussion. Tg18a/c/d/e |
| CA95 | Chertihin and Andrews, “Reactions of laser-ablated Zr and Hf atoms with hydrogen. Matrix infrared spectra of the MH, MH2, MH3, and MH4 molecules,” *J. Phys. Chem.* **99** (1995), 15004–15010. [DOI](https://doi.org/10.1021/j100041a014), [PDF](../materials/valence/Chertihin-Andrews-1995-Zr-Hf-hydrides.pdf) | Supplied PDF reviewed and archived. ZrH3 and both tetrahydrides supported; dihydride spin orderings differ between calculations. ZrH, HfH and HfH3 matrix assignments are tentative. Tg18b–f |
| RB94 | Ram and Bernath, “Fourier transform emission spectroscopy of HfH and HfD,” *J. Chem. Phys.* **101** (1994), 74–79. [DOI](https://doi.org/10.1063/1.468130), [PDF](../materials/valence/Ram-Bernath-1994-HfH-HfD.pdf) | Supplied PDF reviewed and archived. Independent HfH/HfD observation and X²Δ3/2 assignment, pp. 74–75, 78–79; Tg18b |
| MCS96 | Ma, Collins and Schaefer, “Periodic Trends for Transition Metal Dihydrides MH2, Dihydride Dihydrogen Complexes MH2·H2, and Tetrahydrides MH4 (M = Ti, V, and Cr),” *J. Am. Chem. Soc.* **118** (1996), 870–879. [DOI](https://doi.org/10.1021/ja951376t), [PDF](../materials/valence/Ma-Collins-Schaefer-1996-TiH2.pdf) | Supplied PDF reviewed and archived. TiH2 triplet occupation, pp. 872–873, Table 1; the ³A1/³B1 distinction needs no separate atom type. The proposed matrix-band reassignment to TiH2·H2, pp. 877–879, is not used to add an adduct state. Tg18e |
| WDB90 | Wang, Das and Balasubramanian, “Potential energy surfaces for insertion of Zr into H2,” *Mol. Phys.* **69** (1990), 147–159. [DOI](https://doi.org/10.1080/00268979000100101), [PDF](../materials/valence/Wang-Das-Balasubramanian-1990-ZrH2.pdf) | Supplied PDF reviewed and archived. Singlet/triplet minima and explicitly unresolved ordering, p. 152, Table 2; Tg18g |
| DB91 | Das and Balasubramanian, “Potential energy surfaces for insertion of hafnium atoms into hydrogen,” *J. Phys. Chem.* **95** (1991), 3979–3984. [DOI](https://doi.org/10.1021/j100163a018), [PDF](../materials/valence/Das-Balasubramanian-1991-HfH2.pdf) | Supplied PDF reviewed and archived. HfH2 singlet ground state at CASSCF/SOCI and with spin–orbit coupling, pp. 3981–3982, Tables II/V/VI; Tg18e |
| Z20 | Zolnhofer et al., “Electronic Structure and Magnetic Properties of a Titanium(II) Coordination Complex,” *Inorg. Chem.* **59** (2020), 6187–6201. [DOI](https://doi.org/10.1021/acs.inorgchem.0c00311), [PDF](../materials/valence/Zolnhofer-2020-TiII-tmeda.pdf) | Reviewed and archived from the authors' FSU copy. TiCl2(tmeda)2 structure and experimental S = 1; Tg20a |
| L95 | Lukens, *Trivalent Metallocene Chemistry of Some Uranium, Titanium, and Zirconium Complexes*, Berkeley thesis (1995), LBL-37646. [Repository](https://digital.library.unt.edu/ark:/67531/metadc621122/), [PDF](../materials/valence/Lukens-1995-metallocene-thesis.pdf) | Reviewed and archived. Original Cp*2TiH EPR/XRD and Cp3Zr measurements; Tg19a/Tg22h. Supplies the needed hydride data without the unavailable [1998 article](https://doi.org/10.1021/om980601n) |
| C00 | Chirik, *Ancillary Ligand Effects on Fundamental Transformations in Metallocene Catalyzed Olefin Polymerization*, Caltech thesis (2000), chapter 1. [DOI](https://doi.org/10.7907/2qx2-9753), [PDF](../materials/valence/Chirik-2000-metallocene-thesis.pdf) | Reviewed and archived. Terminal Zr/Hf metallocene dihydrides; Tg19b. Contains the study published as [Chirik, Day and Bercaw (1999)](https://doi.org/10.1021/om9810293); that journal PDF is no longer required |
| E97 | Entley et al., “The Hexacyanotitanate Ion,” *J. Am. Chem. Soc.* **119** (1997), 6251–6258. [DOI](https://doi.org/10.1021/ja962773m), [PDF](../materials/valence/Entley-1997-Ti-hexacyanide.pdf) | Reviewed and archived from Girolami's group site. Ti(CN)6³− XRD/EPR/magnetometry and older cyanide-assignment review; Tg22a |
| D20 | Deng et al., “Filling a Gap: The Coordinatively Saturated Group 4 Carbonyl Complexes TM(CO)8 (TM = Zr, Hf) and Ti(CO)7,” *Chem. Eur. J.* **26** (2020), 10487–10500. [DOI](https://doi.org/10.1002/chem.201905552), [PDF](../materials/valence/Deng-2020-group4-carbonyls.pdf) | Reviewed and archived. Calculated neutral Ti(CO)6 triplet, p. 10493; Tg22d. Observed Zr/Hf octacarbonyls support capacity 20 under our donor counting. No MO-based reassignment of donated-pair counts is adopted |
| B13 | Bowman et al., *Inorg. Chem.* **52** (2013), 2242–2256. [DOI](https://doi.org/10.1021/ic302799s) | Primary abstract reviewed; full PDF unavailable. Explicit metal/ligand redox allocation suffices for metal-state reuse in Tg22b; no ligand graph inferred |
| C87 | Cloke et al., “Synthesis of bis(η-1,3,5-tri-t-butylbenzene) sandwich complexes of titanium, zirconium, and hafnium,” *J. Chem. Soc., Chem. Commun.* (1987), 1667–1668. [DOI](https://doi.org/10.1039/C39870001667) | Primary abstract explicitly reports diamagnetism for all three complexes; Tg22e. PDF download failed; archive copy remains desirable |
| M61 | Murray, “A Metal Carbonyl Compound of Titanium,” *J. Am. Chem. Soc.* **83** (1961), 1287–1289. [DOI](https://doi.org/10.1021/ja01467a007), [PDF](../materials/valence/Murray-1961-titanocene-dicarbonyl.pdf) | Supplied PDF reviewed and archived. Cp2Ti(CO)2 synthesis, terminal-CO IR bands and closed-shell electron-count argument, p. 1288; Tg20c. No magnetic measurement |
| M79 | McPherson et al., “Electronic properties of bis(η5-cyclopentadienyl)titanium 2,2′-bipyridyl,” *J. Am. Chem. Soc.* **101** (1979), 3425–3430. [DOI](https://doi.org/10.1021/ja00507a001) | Indexed p. 3425 explicitly describes Cp2Ti(CO)2 as a diamagnetic monomer; Tg20c. This is a literature statement in the introduction, not this paper's new measurement. Full PDF unavailable; not an assignment blocker |
| EC90 | Ellis and Chi, *J. Am. Chem. Soc.* **112** (1990), 6022–6025. [DOI](https://doi.org/10.1021/ja00172a017), [PDF](../materials/valence/Ellis-Chi-1990-group4-hexacarbonylates.pdf) | Supplied PDF reviewed and archived. Octahedral Ti/Zr/Hf(CO)6²− structures, low-spin d⁶ rationale and comparative ¹³C NMR, pp. 6022, 6024; Tg22c. No separate magnetic measurement |
| CG79 | Cloke and Green, “Synthesis of zerovalent bis-(η-arene)(trimethylphosphine)-hafnium and -zirconium compounds using metal vapours,” *J. Chem. Soc., Chem. Commun.* (1979), 127–128. [DOI](https://doi.org/10.1039/C39790000127), [PDF](../materials/valence/Cloke-Green-1979-Zr-Hf-bisarenes-PMe3.pdf) | Supplied PDF reviewed and archived. Zr/Hf(arene)2(PMe3) preparation, resolved ¹H NMR and proposed bent geometry, p. 127; Tg22f. The paired occupation is our inference from the reported spectroscopy |

Intermediates are in `scratch/hw-titanium-review/`; retained source PDFs are in
`materials/valence/`. The available full papers and the explicitly identified
abstract/page excerpts have been reviewed; 16 PDFs are archived. Source research,
selection, registry implementation and conformance coverage are complete for the
selected scope above.

### Vanadium group

Source: [HW 103, Band II, chapter XXVI](<../materials/valence/Holleman-Wiberg - Anorganische Chemie - 103. Auflage - Band 2.pdf#page=323>),
printed pp. 1819–1841; references pp. 1842–1843. Extraction: 2026-11-03.
Baseline: V 4, Nb 2, Ta 1 and Db 1 states. All eight are retained with checked
registry citations. Db chemistry is referred to chapter XXXVII on p. 1819.

#### Vanadium-group work table

M = V, Nb or Ta unless restricted; X = halide; D = neutral pair donor.
States use concrete, zero-valence defaults. Page references are to HW II.
Ordinary M–X/M–R/M–H bonds, anionic donor alternatives and Cp− donation follow
the selected PW convention; terminal M=O remains localized. All 160 selected
additions below (V 57, Nb 48, Ta 55) are **implemented and covered**. Counts exclude
states already listed in earlier rows. Deferred and omitted rows add no states.

| ID | Atom descriptions / gap | Disposition | Examples and evidence |
| --- | --- | --- | --- |
| Vg01 | V/Ta/Db#n#u3; Nb#u5; V#c+2#u3; V/Nb#c+3#u2; V#c+2#u#v2 | Retained all 8; citations checked | Free atoms: 1538, Table 19.1, including calculated Db. V²⁺/V³⁺: 1828–1830; Nb³⁺: 1833, 1837. Vanadyl VO²⁺: 1826; its comment distinguishes it from dioxovanadium VO2⁺ |
| Vg02 | M#c+4#u; M#c+5; Ta#c+3#u2 | Implemented: 7 | Ionic descriptions of IV/V oxides and halides, 1821–1822, 1825–1826, 1833–1837; Ta(III), 1837. Free-ion occupations remain separate from complex occupations |
| Vg03 | M#v5; M#c+5#t5 | Implemented: 6 | Molecular VF5 and Nb/TaX5, 1822, 1834–1836; Nb/TaMe5, 1841. Ordinary terminal bonds or five anionic donors. V#v5 also covers localized VOX3, orthovanadate/protonated vanadates, V2O7⁴− and cyclic V4O12⁴−, 1822–1824; no VCl5/VBr5 or VMe5 example inferred |
| Vg04 | M#c-#v6; M#c+5#t6 | Implemented: 6 | VF6− and Nb/TaX6−, 1822, 1835; MMe6−, 1831, 1841. Six ordinary bonds or six anionic donors |
| Vg05 | Nb/Ta#c-2#v7; #c-3#v8; #c+5#t7; #c+5#t8 | Implemented: 8 | MF7²− and MF8³−, 1835. Seven/eight ordinary bonds or halide donors. M(O2)4³−, 1833, reuses the eight-coordinate states with the peroxide O–O bonds retained |
| Vg06 | Nb/Ta#v5#t; #v5#t2; V#v5#t2 | Implemented: 5; donor alternatives reuse earlier rows | Nb/TaX5·D, 1835, and MMe5(dmpe), 1841. VOCl3·2D, 1822, retains V=O; its halide-donor description is in Vg09 |
| Vg07 | Nb/Ta#c+#v4#t2 | Implemented: 2; all-halide-donor forms reuse Vg04 | M2X10 dimers and (MF5)4 rings, 1834–1836, Fig. 26.7. Four ordinary terminal halides and two bridging X− donors at each metal; existing X−#d2 states suffice. No MC addition |
| Vg08 | V#c+#v4; V#c+#v4#t4 | Implemented: 2 | Dioxovanadium VO2⁺ and [VO2(H2O)4]⁺, 1821, 1824: two localized V=O bonds, with four water donors in the hydrate |
| Vg09 | M#c+3#v2#t3; V#c+3#v2#t5; Nb#c+3#v2#t4; Nb#c+3#v2#t5; Nb#c+3#v2#t6 | Implemented: 7 | MOX3 and VOCl3·2D, 1822, 1835–1836; NbOF4−/NbOF5²−/NbOF6³−, 1836. Localized M=O plus halide/neutral donors. Ordinary Nb–F alternatives reuse Vg04–Vg05 |
| Vg10 | V#c-#v6#t; V#c-3#v8; V#c+5#t8 | Implemented: 3 | [VO(O2)2(NH3)]− and V(O2)4³−, 1825. Peroxide O–O bonds retained; ordinary V–O or peroxide O− donors. VO2(O2)2³− reuses V#c-3#v8; the oxoperoxo donor form retains V=O and reuses Vg09 |
| Vg11 | V/Nb#u#v4; #c+4#u#t4 | Implemented: 4 | Molecular VCl4/VBr4 and VR4, 1825, 1831; tetrahedral Nb(NEt2)4, 1833. Ordinary bonds or anionic donors, leaving one metal electron. No free TaX4 molecule inferred from its solid formula |
| Vg12 | M#c-2#u#v6; M#c+4#u#t6 | Implemented: 6 | V(IV)X6²−, 1825–1826, and Nb/Ta(IV)X6²−, 1836. Six ordinary bonds or six halide donors, d¹ |
| Vg13 | V#u#v4#t2; M#u#v4#t4; M#c+4#u#t8 | Implemented: 7 | VCl4·2D and MCl4(diars)2, 1822, 1825–1826, 1836. Mixed ordinary/donor and all-donor forms; shared four/eight-donor counts. V#u#v4#t2 also covers localized VO(acac)2 |
| Vg14 | Nb#c-3#u#v7; Nb#c+4#u#t7; Nb#c-4#u#v8 | Implemented: 3; eight-donor form reused | NbF7³− and Nb(CN)8⁴−, 1833, 1836, 1838. One remaining electron; ordinary or anionic-donor descriptions. KG75, p. 2490, confirms the IV cyanide's paramagnetism; the paired III cyanide is in Vg30 |
| Vg15 | V#c+2#u#v2#t4; V#c+2#u#v2#t5; V#c-#u#v5; V#c+2#u#v2#t3; V#c-3#u#v7 | Implemented: 5; bare vanadyl retained | [VO(H2O)5]²⁺, VO(acac)2, VO(OH)3− and VO(CN)5³−, 1822, 1826–1827. Localized V=O; the remaining ligands are ordinary or donors. VO(CN)5³− is distinguished from the explicitly unknown V(CN)6²− |
| Vg16 | V#u2#v3; V#c+3#u2#t3; V#c-#u2#v4; V#c+3#u2#t4 | Implemented: 4 | V[N(SiMe3)2]3 and VCl4−/VBr4−, 1822, 1829. Trigonal-planar/tetrahedral V(III) examples; ordinary or donor descriptions |
| Vg17 | V#u2#v3#t2; #u2#v3#t3; #c+#u2#v2#t4; #c+3#u2#t5; #c+3#u2#t6; #c-3#u2#v6; #c+2#u2#v#t5 | Implemented: 7 | VX3·2D, VX3·3D, VX2D4⁺, [V(H2O)6]³⁺, VX6³− and [V(OH)(H2O)5]²⁺, 1828–1830. HW gives the V(III) magnetic moments supporting two unpaired electrons |
| Vg18 | V#c-4#u2#v7; V#c+3#u2#t7 | Implemented: 2 | V(CN)7⁴−, 1829: explicitly two unpaired electrons and seven-coordinate V; ordinary C-bound cyanides or CN− donors |
| Vg19 | V#u3#v2#t4; V#c+2#u3#t6; V#c-4#u3#v6 | Implemented: 3 | V(II)X2D4, [V(H2O)6]²⁺ and V(CN)6⁴−, 1829–1830. Three unpaired electrons are explicit in HW. VPh6⁴− reuses the six-ordinary-bond state |
| Vg19b | V#u3#v2#t4; V#c+2#u3#t6 | Reuse Vg19; no additional state | VCl2(dmpe)2 and VMe2(dmpe)2: G85, p. 1340, Table 2, gives 3.7 ± 0.1 BM for both and explicitly lists three unpaired V electrons. The paper's collective “low-spin” description does not mean #n#u for d³; its ESR analysis is S = 3/2, pp. 1341, 1343 |
| Vg20 | Nb/Ta#u2#v3#t3; #c+3#u2#t6 | Implemented: 4 | NbCl3(py)3, NbBr3(PMe2Ph)3 and TaX3(PMe3)3, 1837: HW explicitly assigns two unpaired electrons. Ordinary halides plus neutral donors, or all donors |
| Vg21 | M#c-#n3#t6; M#c-3#n4#t5; M#c+#n2#t7 | Implemented: 9 | M(CO)6−, M(CO)5³− and CpM(CO)4, 1822, 1830, 1833, 1841, 2111. HW's closed-shell/low-spin carbonyl descriptions; CO donates one pair and Cp− three pairs (PW) |
| Vg22 | V#n2#u#t6 | Implemented: 1 | V(CO)6, 1822, 1830, 2109–2111: low-spin d⁵, one unpaired electron. V(arene)2 reuses the state; HW 2207, Table 32.25, gives 1.68 BM |
| Vg23 | M#c+2#u#v2#t6; V#c+3#v2#t6; Nb/Ta#c+2#v3#t6; #c+5#t9 | Implemented: 8; other all-donor states reused | Cp2MX2, [Cp2VCl2]⁺ and Cp2Nb/TaX3, 1831, 1841. Two η⁵-Cp− donors plus ordinary X or X− donation. Nb/TaCp4 has two η⁵ and two η¹ rings and reuses the IV counts; Cp2MeTa=CH2 reuses the V counts |
| Vg24 | V#c+2#u3#t6 | Reuse Vg19 | Vanadocene, HW 1830. J12 establishes S = 3/2 by EPR/HFEPR; use this evidence rather than the inconsistent 2.84 BM printed in HW 2190, Table 32.22 |
| Vg25 | V#u3#v2; V#u3#v#h; V#u3#h2; V#c+2#u3#t2 | Implemented: 4 | Molecular VH2: MCS96, pp. 872–874, Table 2, calculates a quartet ground state; VLW95, pp. 4367, 4370–4372, establishes S = 3/2 by ESR. Ordinary, mixed/implicit H or H− donors. The geometry discussion does not change the local occupation |
| Vg26 | V#c+2#v3#t3 | Implemented: 1 | [Cl3VN]−, HW I, pp. 757–758, Fig. 14.6g: three ordinary V–Cl bonds plus the existing N³− three-pair donor. The revised `v_group/vcl3-nitride-.edn` infers the metal occupation. The fully ordinary V≡N graph reuses Vg04 |
| Vg27 | Ta#v5#t4; #v4#h#t4; #v3#h2#t4; #v2#h3#t4; #v#h4#t4; #h5#t4; Ta#u#v3#h#t4; #u#v2#h2#t4; #u#v#h3#t4; #u#h4#t4 | Implemented: 10; explicit tetrahydride and H−-donor states reused | TaH5(dmpe)2 and TaH4(dmpe)2, 1834, 2068–2069. Five/four terminal hydrides and four phosphine donors; d⁰/d¹. Explicit, mixed and implicit H; all-H− alternatives reuse Ta#c+5#t9 and Ta#c+4#u#t8 |
| Vg28a | Nb/Ta#u#v3#h; #u#v2#h2; #u#v#h3; #u#h4; Ta#u#v4; Ta#c+4#u#t4 | Implemented: 10; Nb explicit/donor forms reuse Vg11 | Molecular NbH4 and TaH4: WA11, pp. 14179–14182, Tables 5–6, assigns doublet D2d structures. Ordinary, mixed/implicit H or four H− donors; no change to the already accepted coordinated tetrahydrides |
| Vg28b | Nb/Ta#u3#v2; #u3#v#h; #u3#h2; #c+2#u3#t2 | Implemented: 8 | Molecular NbH2/TaH2: WA11 establishes the molecules. VLW95, p. 4367, reports calculated NbH2 ⁴B1 and three competing TaH2 quartet states (⁴B1/⁴A2/⁴B2); their ordering does not affect #u3. This is theoretical occupation evidence, not an ESR determination for these dihydrides. Ordinary, mixed/implicit H or H− donors |
| Vg28c | VH2(H2), NbH4(H2)4 and TaH4(H2)4 | Deferred; no new states | WA11, pp. 14178–14182, explicitly distinguishes the intact H2 ligands from terminal hydrides. Keep with Cg33/Sg09h: H2-donor accounting is separate work; no MC replacement or VH4 assignment inferred from these formulae |
| Vg29a | Nb#u#v3#h#t4; Nb#u#v2#h2#t4 | Implemented: 2; explicit, Ta and all-anionic-donor forms reused | NbCl2H2(dmpe)2, TaCl2H2(dmpe)2 and TaCl2H2(PMe3)4: L84, pp. 1719–1721, establishes terminal hydrides and one unpaired electron by ESR and moments of 1.72–1.74 BM. Nb mixed/implicit H adds two states; Vg13/Vg27 cover the rest |
| Vg29b | Nb/Ta#n#u#v2#t4; #c+2#n#u#t6 | Omitted: occupation insufficiently established | Nb/TaCl2L4 (L = PMe3 or one dmpe P donor), HW 1837. L83, p. 4474, footnote 12, gives 1.45/1.49 BM for TaCl2(PMe3)4; DW77, p. 1134, gives ESR g = 1.943 for TaCl2(dmpe)2. Neither paper explicitly assigns the occupation. L84/P25/BLH81 do not supply the Nb(II) occupation. User decision 2026-11-03: omit all four uncertain states |
| Vg29c | Ta#n2#v#t6; Ta#n2#h#t6 | Implemented: 2; Ta#c+#n2#t7 reused | TaX(CO)2(dmpe)2, X = Cl, Br, H, Me or CN: DW77, pp. 1134–1136. One ordinary X/H bond and six neutral donor pairs; the paired d⁴ assignment follows the resolved ¹H/³¹P NMR and the selected closed-shell carbonyl convention. The hydride has an implicit-H form; all-anionic-donor forms reuse Vg21 |
| Vg29d | Ta#n#v3#t4; #n#v2#h#t4; #n#v#h2#t4; #c+3#n#t7 | Implemented: 4 | TaClH2(PMe3)4: L83, p. 4475, explicitly identifies the diamagnetic Ta(III) monomer and terminal hydrides. Paired d²; ordinary, mixed/implicit H or all-anionic-donor forms. Distinct from the six-coordinate unpaired III complexes in Vg20 |
| Vg29e | Nb#n2#v#t6 | Implemented: 1; Nb#c+#n2#t7 reused | NbCl(CO)2(dmpe)2: BLH81, pp. 794–795, 798, gives the synthesis, IR, mass spectrum and expected phosphine ¹H NMR signals. The paired d⁴ assignment follows the selected closed-shell carbonyl convention, as for Ta in Vg29c; it is not a separately measured spin assignment. One ordinary Nb–Cl bond and six neutral donor pairs, or Cl− donation. HW 1837's other MX·6D carbonyls reuse these counts |
| Vg30 | Nb#c-5#n#v8; Nb#c+3#n#t8 | Implemented: 2 | Nb(CN)8⁵−, HW 1838; KG75, p. 2490, explicitly states that the III salts are diamagnetic. Eight ordinary C-bound cyanides or eight CN− donors, with the two remaining metal electrons paired |
| Vg31 | Nb/Ta#v5#t3; Ta#c+3#v2#t6 | Implemented: 3; Nb#c+3#v2#t6 reused | Nb2X9³− (X = Cl, Br, I) and Ta2Cl9³−: HW 1837; S98, pp. 6802–6804, 6806, explicitly assigns a metal–metal double bond with two unpaired bonding electrons. Encode M=M as `2#u2#s3` (PW), with three terminal ordinary halides and three bridging X− donors per metal, or all six halides as donors. No atom-local unpaired electron or MC type is added |
| Vg32 | V#c+2#u2#v#t6; V#c+3#n#t8 | Implemented: 2; V#c+3#u2#t7 reused | Cp2VCl: FS77, pp. 46–47, explicitly assigns two unpaired electrons from 2.80(5) BM. Ordinary terminal Cl or Cl− donation; the paper extends the occupation to Cp2VR. Its pp. 51–52 separately identify diamagnetic [Cp2VL2]⁺ complexes (L = CO/phosphine), giving the paired eight-accepted-pair state |
| Vg33 | V#n2#v#t6; V#n2#h#t6 | Implemented: 2; V#c+#n2#t7 reused | KH74, pp. 201–202, establishes diamagnetic VH(PF3)6, including the hydride resonance and molecular ion. HW 1828's VH(PF3)5 formula is incorrect. Ordinary/implicit H plus six phosphane donors, or H− donation |
| Vg34 | Vanadates/thiovanadates, extended binary solids, polyoxo-/metal clusters | Reuse where local states suffice; no additional default types selected | VO4³−/VS4³−/VSe4³− and small linked tetrahedral vanadates reuse localized #v5, 1823–1825. No aggregate MC state from bulk oxides, nitrides, carbides, borides, M3/M6 halide clusters or the large polyoxometalates, 1827–1830, 1837–1841. M2X9³− is covered separately above |
| Vg35 | VH2·2BH3; V2(form)4; Np3V–N2–VNp3; rearranged NbCp2 dimer; further bis-arene/bipyridyl redox states | Unselected in this pass | Extraction retained from 1828–1831, 1833, 1841 and cross-references. No broad borohydride MC state, scaffold-specific metal state or formal metal-redox ladder inferred from these formulae. The common Cp and carbonyl classes are selected above |

#### Vanadium-group limits and coverage

All 160 implemented additions satisfy the nonbonding budget and the existing
capacity 18; none is excluded by that bound. No new target-covalence or aromatic
ValenceTable entry is identified for these graphs with zero or supplied implicit H.
Counts uses the entry shifted by `2t − charge`, while capacity belongs to the
actual element. With occupation absent it pairs the remaining electrons; the
unpaired V(II/III), Nb/Ta(III) and quartet molecular MH2 inputs therefore need supplied
#u3/#u2 for agreement with AtomTyping. The 166 fixture inputs cover every addition,
including mixed/implicit hydrides and both selected ligand descriptions. All four
resolver/policy outputs agree and pass the valence-invariants validator; 688 checks
verify the expected metal states. The Nb2Cl9³−/Ta2Cl9³− fixtures retain triplet M=M
bond occupation with either ordinary terminal halides or all-halide donation.

The suite adds 165 independent inputs and revises [Cl3VN]− to infer its metal
occupation. All 3,594 graph unit tests, 223 property tests and 1,976 resolution
cases pass; all 1,811 earlier snapshots are unchanged. Existing H−, halide,
peroxide, water, amine, phosphine, carbonyl and Cp donor types suffice. The existing
multidonor accounting issue remains outside this work. No resolver, capacity or
ValenceTable changes were needed.

#### Vanadium-group primary-source checks

| Key | Source | Status and use |
| --- | --- | --- |
| J12 | Jackson et al., “Vanadocene de Novo: Spectroscopic and Computational Analysis of Bis(η⁵-cyclopentadienyl)vanadium(II),” *Organometallics* **31** (2012), 8265–8274. [DOI](https://doi.org/10.1021/om300892y), [PDF](../materials/valence/Jackson-2012-vanadocene.pdf) | Author copy obtained from the National High Magnetic Field Laboratory, reviewed and archived; pp. 8265–8267 establish the quartet. Vg24 |
| MCS96 | Ma, Collins and Schaefer, *J. Am. Chem. Soc.* **118** (1996), 870–879. [DOI](https://doi.org/10.1021/ja951376t), [PDF](../materials/valence/Ma-Collins-Schaefer-1996-TiH2.pdf) | Existing archived paper; VH2 section reviewed, pp. 872–874, Table 2. Vg25 |
| WA11 | Wang and Andrews, “Tetrahydrometalate Species VH2(H2), NbH4, and TaH4: Matrix Infrared Spectra and Quantum Chemical Calculations,” *J. Phys. Chem. A* **115** (2011), 14175–14183. [DOI](https://doi.org/10.1021/jp2076148), [PDF](../materials/valence/Wang-Andrews-2011-group5-hydrides.pdf) | Supplied PDF reviewed and archived. NbH4/TaH4 doublets, Tables 5–6; VH2(H2) quartet and intact H2 adducts, pp. 14179–14182. Bare NbH2/TaH2 occupations are not labelled. Vg28 |
| L84 | Luetkens, Elcesser, Huffman and Sattelberger, “Paramagnetic hydride complexes of niobium(IV) and tantalum(IV),” *Inorg. Chem.* **23** (1984), 1718–1726. [DOI](https://doi.org/10.1021/ic00180a017), [PDF](../materials/valence/Luetkens-1984-Nb-Ta-hydrides.pdf) | Supplied PDF reviewed and archived. Terminal M(IV) hydrides and d¹ occupation, pp. 1719–1721; no explicit spin assignment for the M(II) precursors. Vg29 |
| KG75 | Kiernan and Griffith, “Studies on transition-metal cyano-complexes. Part I. Octacyanoniobates(III), -niobates(IV), -molybdates(V), and -tungstates(V),” *J. Chem. Soc., Dalton Trans.* (1975), 2489–2494. [DOI](https://doi.org/10.1039/DT9750002489), [PDF](../materials/valence/Kiernan-Griffith-1975-octacyanoniobates.pdf) | Supplied PDF reviewed and archived. Page 2490 explicitly distinguishes diamagnetic III from paramagnetic IV salts. Vg14/Vg30 |
| S98 | Stranger, McGrady and Lovell, “Metal–Metal Bonding in d¹d¹ and d²d² Bioctahedral Dimer Systems,” *Inorg. Chem.* **37** (1998), 6795–6806. [DOI](https://doi.org/10.1021/ic971511y), [PDF](../materials/valence/Stranger-McGrady-Lovell-1998-metal-metal-bonding.pdf) | Supplied PDF reviewed and archived. Nb2X9³−/Ta2Cl9³− double bond: paired σ electrons and two parallel-spin electrons in the δπ bonding orbitals, pp. 6802–6804, 6806. The weakly coupled V2Cl9³− quintet is distinct. Vg31 |
| FS77 | Fieselmann and Stucky, “The crystal structure of dicyclopentadienylvanadium monochloride and its implications for the structures of other d² dicyclopentadienyl compounds,” *J. Organomet. Chem.* **137** (1977), 43–54. [DOI](https://doi.org/10.1016/S0022-328X(00)80841-4), [PDF](../materials/valence/Fieselmann-Stucky-1977-Cp2VCl.pdf) | Supplied PDF reviewed and archived. Two unpaired electrons in Cp2VCl, pp. 46–47; paired [Cp2VL2]⁺ complexes, pp. 51–52. Vg32 |
| KH74 | Kruck and Hempel, “Hexakis(trifluorophosphane)vanadium Hydride,” *Angew. Chem. Int. Ed. Engl.* **13** (1974), 201–202. [DOI](https://doi.org/10.1002/anie.197402011), [PDF](../materials/valence/Kruck-Hempel-1974-VH-PF3-6.pdf) | Supplied PDF reviewed and archived. VH(PF3)6, diamagnetism and seven-coordinate V; corrects HW's five-PF3 formula. Vg33 |
| P25 | Pampaloni et al., “Oxidative Chlorination Reactions of [Nb(CO)6]− in the Presence of N- and P-Donors,” *Eur. J. Inorg. Chem.* (2025), e202500026. [DOI](https://doi.org/10.1002/ejic.202500026), [PDF](../materials/valence/Pampaloni-2025-Nb-carbonyl-chlorination.pdf) | Author manuscript obtained from Bologna's repository and archived. Confirms NbCl2(dmpe)2 preparation; no explicit spin assignment. Vg29b omitted |
| L83 | Luetkens, Huffman and Sattelberger, “Development of low-valent tantalum chemistry,” *J. Am. Chem. Soc.* **105** (1983), 4474–4475. [DOI](https://doi.org/10.1021/ja00351a058), [PDF](../materials/valence/Luetkens-1983-low-valent-Ta.pdf) | Supplied PDF reviewed and archived. TaCl2(PMe3)4 moments, p. 4474, footnote 12, without an explicit spin assignment; diamagnetic terminal TaClH2(PMe3)4, p. 4475. Vg29b/Vg29d |
| DW77 | Datta and Wreford, “Seven-coordinate complexes of tantalum(I),” *Inorg. Chem.* **16** (1977), 1134–1137. [DOI](https://doi.org/10.1021/ic50171a031), [PDF](../materials/valence/Datta-Wreford-1977-TaI-complexes.pdf) | Supplied PDF reviewed and archived. Ta(II) ESR, p. 1134; Ta(I) carbonyl structures, hydride identification and resolved NMR, pp. 1134–1136. The paired I occupation is inferred from the spectroscopy; no explicit II spin assignment. Vg29b/Vg29c |
| VLW95 | Van Zee, Li and Weltner, “V, Nb, and Ta hydride molecules in deuterium and rare-gas matrices: Infrared and electron spin resonance spectra,” *J. Chem. Phys.* **102** (1995), 4367–4374. [DOI](https://doi.org/10.1063/1.469485), [PDF](../materials/valence/Van-Zee-Li-Weltner-1995-group5-hydrides.pdf) | Supplied PDF reviewed and archived. Calculated NbH2/TaH2 quartet states summarized on p. 4367; direct VH2 quartet and NbH4/TaH4 doublet ESR evidence. Vg25/Vg28 |
| G85 | Girolami, Wilkinson, Galas, Thornton-Pett and Hursthouse, “Synthesis and properties of the divalent 1,2-bis(dimethylphosphino)ethane (dmpe) complexes MCl2(dmpe)2 and MMe2(dmpe)2 (M = Ti, V, Cr, Mn, or Fe),” *J. Chem. Soc., Dalton Trans.* (1985), 1339–1348. [DOI](https://doi.org/10.1039/DT9850001339), [PDF](../materials/valence/Girolami-1985-divalent-dmpe-complexes.pdf) | Supplied PDF reviewed and archived. Three unpaired V electrons and 3.7 ± 0.1 BM for VCl2(dmpe)2/VMe2(dmpe)2, p. 1340, Table 2; quartet ESR treatment, pp. 1341, 1343. Vg19b reuses Vg19 |
| BLH81 | Burt, Leigh and Hughes, “Carbonyl complexes of niobium and tantalum and evidence of a niobium dinitrogen complex,” *J. Chem. Soc., Dalton Trans.* (1981), 793–799. [DOI](https://doi.org/10.1039/DT9810000793), [PDF](../materials/valence/Burt-Leigh-Hughes-1981-Nb-Ta-carbonyls.pdf) | Supplied PDF reviewed and archived. NbCl(CO)2(dmpe)2 and its Br/N3 substitutions, pp. 794–795, 798; NMR supports the selected closed-shell carbonyl convention. The reported ESR concerns Nb(IV), not the Nb(II) phosphine occupation. Vg29e |

Intermediates and the deduplication/count audit are in `scratch/hw-vanadium-review/`.
Retained papers are in `materials/valence/`; neither the work table nor its evidence
depends on keeping the scratch directory.

### Chromium group

HW II, chapter XXVII, pp. 1844–1898 (chemistry through p. 1895), with the cited
chapter-XXXII cross-references. Updated 2026-11-04. The accepted set of **260 new
metal states: Cr 72, Mo 94, W 94**, plus **two O bridge companions**, is implemented
with inferred conformance coverage.
All 11 earlier Cr/Mo/W/Sg states are retained with checked citations. Thirteen W
phosphine-hydride occupations remain open by decision and are outside the selected
counts, as are the other open rows.

#### Chromium-group work table

E = Cr, Mo, W; M = Mo, W. A prefix applies to every following suffix in its cell.
Counts are new states after reuse across rows. The 260-state selection was
accepted and implemented on 2026-11-04; fixture verification is recorded below.
Ionic ligand allocation follows the convention used in the preceding
transition-metal chapters. Oxo bonds remain ordinary E=O unless stated otherwise.

| ID | States / extracted family | Status; new states | Examples and evidence |
| --- | --- | --- | --- |
| Crg01 | Cr#u6; Cr#c+2#u4; Cr#c+2#n#u2; Cr#c+3#u3; Cr#n3#t5#m0; Mo#u6; Mo#c+2#u4; Mo#c+2#n#u2; Mo#c+3#u3; W#n#u4; Sg#n#u4 | Retained; 11 existing | Atomic occupations: HW II, p. 1538; Cr(II) high-/low-spin and Cr(III), pp. 1861–1864; Mo(II/III), pp. 1879, 1891. Existing Cr bridge: pp. 2069, 2137 and D79. Citations/comments updated; no existing occupation is removed |
| Crg02 | E#c+4#u2; #c+5#u; #c+6; W#c+2#u4; #c+3#u3 | Implemented; 11 | Ionic states for the chapter's II–VI chemistry; d⁴/d³/d²/d¹/d⁰ free-ion occupations. These do not prescribe the occupations of coordinated complexes. HW 1849, 1871, 1891 |
| Crg03 | E#v6; M#c+6#t6 | Implemented; 5 | Chromate/dichromate, molybdate/tungstate and localized E(VI) oxo graphs; MF6, WCl6/WBr6, MMe6, alkylidene/alkylidyne and terminal nitrido derivatives. HW 1849–1853, 1872–1873, 1880, 1893–1894. MF6 has ordinary or six F−-donor descriptions; CrF6 is not an example |
| Crg04 | M#c-#v7; #c-2#v8; #c+6#t7; #c+6#t8 | Implemented; 8 | MF7− and MF8²−, 1872; MMe7−, 1894. Ordinary or anionic-donor ligands. M(O2)4²−, 1882, reuses the eight-coordinate states while retaining each O–O bond |
| Crg05 | E#c+4#v2#t4; #c+2#v4#t2 | Implemented; 6 | EO X4 and EO2 X2: ordinary oxo bonds plus X− donors. Cr, 1853, 1855; Mo/W, 1873. Fully ordinary graphs reuse E#v6 |
| Crg06 | Cr#c-#v7; #c-2#v8; #c+2#v4#t4 | Implemented; 3 | CrO(O2)2(OH)− and halogenochromates, 1853–1854; CrO2F4²−, 1849. Ordinary ligand bonds, or F− donation with the two Cr=O bonds retained |
| Crg07 | Cr#v6#t; #v6#t2; #c+4#v2#t5; #c+4#v2#t6 | Implemented; 4 | CrO(O2)2·D and CrO(O2)2(bipy), 1853–1854: one/two neutral donor pairs, with ordinary Cr–peroxo bonds or O− donation. The latter allocates one negative charge to each peroxide O, retaining O–O |
| Crg08 | M#c+2#v4#t4; #v6#t2; #v6#t3; #c+4#v2#t6; Mo#c+2#v4#t3; #c-3#v9; W#v6#t | Implemented; 11 | MoO2(H2O)4²+, MoO2Cl2(H2O)2, MO2F4²−, MoO2F3−, MoO3F3³−, WO3X−/WO3F2²−, MO3(dien), and MO(O2)2(H2O)2. HW 1870, 1873, 1881–1882. Other ordinary/all-halide-donor forms reuse preceding rows. No state is inferred from the inconsistent WO3F3²− formula printed on p. 1881 |
| Crg09 | E#u#v5; #c+5#u#t5; #c-#u#v6; #c+5#u#t6 | Implemented; 12 | CrF5/CrF6−; molecular MCl5 and MX6−; MMe5. HW 1854–1855, 1872, 1894: d¹. Ordinary or anionic-donor ligand descriptions |
| Crg10 | M#c-3#u#v8; #c+5#u#t8; #c+#u#v4#t2 | Implemented; 6 | MF8³− and M2Cl10, 1872. Each chloride dimer has four terminal ordinary bonds and two bridging Cl− donors; the all-chloride-donor graph reuses M#c+5#u#t6. HW explicitly excludes an M–M bond in these dimers. The eight-ligand states also cover M(CN)8³−, 1876 |
| Crg11 | E#c+3#u#v2#t3; #c+3#u#v2#t5; #c-2#u#v7; Cr#c+3#u#v2#t4; Mo#c-#u#v6#t | Implemented; 11 | EO X3, EO X5²−; CrOX4−; MoO(H2O)5³+ and MoOCl4(H2O)−. Cr, 1854; Mo/W, 1873, 1879. One unpaired electron, ordinary E=O, with ordinary halides or X− donors |
| Crg12 | Cr#c-3#u#v8; #c+5#u#t8 | Implemented; 2 | Cr(O2)4³−, 1855–1856: isolated paramagnetic Cr(V) tetraperoxochromate. Ordinary Cr–O or peroxide-O− donation; O–O bonds retained |
| Crg13 | E#u2#v4; #c+4#u2#t4; #c-2#u2#v6; #c+4#u2#t6 | Implemented; 12 | Tetrahedral CrX4 and EAr4; octahedral EX6²−. HW 1854–1856, 1866, 1872, 1894 and its d² ligand-field treatment. SDV25 explicitly calculates triplet tetraaryl ground states. Seven-coordinate IV fluorides remain separate below |
| Crg14 | M#c-4#n#v8; #c+4#n#t8 | Implemented; 4 | M(CN)8⁴−: explicitly diamagnetic d², 1876; KG75. Ordinary C-bound cyanides or CN− donors. The d¹ M(CN)8³− states reuse Crg10 |
| Crg15 | Cr#u2#v4#t3; #c+4#u2#t7 | Implemented; 2 | Cr(O2)2(NH3)3, 1856. CND01, pp. 599–600, determines S = 1 independently from susceptibility and heat capacity. Ordinary Cr–peroxo bonds or peroxide-O− donation; three NH3 donors. No occupation is transferred to Cr(O2)2(CN)3³− without evidence |
| Crg16 | Cr#c-#u3#v4; #c+3#u3#t4; #u3#v3#t2; #c+3#u3#t5 | Implemented; 4 | Four- and five-coordinate Cr(III): CrCl4− and CrCl3(NMe3)2; HW 1860–1861. Ordinary anionic ligands or their donor forms. The trigonal trisamide is separated below |
| Crg17 | Cr#c+3#u3#t6; #c+2#u3#v#t5; #c+#u3#v2#t4; #u3#v3#t3; #c-#u3#v4#t2; #c-2#u3#v5#t; #c-3#u3#v6; M#u3#v3#t3; #c-3#u3#v6; #c+3#u3#t6 | Implemented; 13 | Cr(NH3)6−k Xk (k = 0–6), aquochloro ions, CrCl3(THF)3, Cr(acac)3 and Reinecke's salt; 1860–1861. MCl3(py)3/MCl6³−, 1891; Mo(H2O)6³+, 1879. Three unpaired electrons. Cr hydroxo dimers reuse the metal state but their OH− donor graph needs Crg37 |
| Crg18 | Cr#u4#v2#t2; #c+2#u4#t4; #c+#u4#v#t4; #c+2#u4#t5; #c-3#u4#v5; #c+2#u4#t6; #c-4#u4#v6 | Implemented; 7 | High-spin Cr(II): Cr(acac)2/CrCl2(MeCN)2, CrBrL+ (L = four N donors), Cr(CN)5³−, aqua/ammine ions and discrete CrI6⁴− in Tl4CrI6. HW 1862–1864 explicitly distinguishes the four-unpaired class |
| Crg19 | Cr/Mo#n#u2#v2#t4; #c+2#n#u2#t6; Cr#c-4#n#u2#v6 | Implemented; 5 | CrCl2(diars)2/CrCl2(dmpe)2 and Cr(CN)6⁴−, 1862, 1864; G85 confirms two unpaired electrons for the dmpe complexes. MoCl2(PR3)4 belongs to HW's two-unpaired class, 1891. CrCp2 reuses the Cr six-pair donor state, 2190 |
| Crg20 | Mo#c-4#n#u#v7; #c+3#n#u#t7 | Implemented; 2 | Mo(CN)7⁴− is a doublet, not an example of the generic three-unpaired Mo(III) statement on HW 1891. R71, p. 51, gives S = ½, susceptibility and EPR; MDP19 explicitly uses S = ½. Ordinary cyanides or CN− donors |
| Crg21 | E#n3#t6; #c-2#n4#t5; #c-#n3#v#t5; #c+#n2#u#t6 | Implemented; 12 | E(CO)6, E(CO)5²−, E2(CO)10²− (ordinary E–E single bond), E(C6H6)2 and its monocation. HW 1861, 1871, 1895, 2111, 2207. Neutral carbonyl/bis-arene class paired; bis-arene monocations have one unpaired electron. No additional metal MC state for the haptic arene ligands |
| Crg22 | Cr#v6#t3; E#c+2#v4#t5; reuse E#v6#t2, E#c+2#v4#t4 and E#c-2#v8 | Implemented; 4 | Quadruple-bonded E2(OAc)4, its axial donor adducts, Cr2Me8⁴−, M2Cl8⁴− and Mo2(H2O)8⁴+. HW 1864–1866, 1880, 1891–1892. Ordinary E–E order 4; mixed ordinary/donor acetate or all-O donation. No MC description |
| Crg23 | Mo#v6#t; M#c+3#v3#t3; #c+3#v3#t4; #c+3#v3#t6; other ordinary forms reused | Implemented; 7 | Triple-bonded M2(OR)6/M2(NR2)6 and their two-donor adducts; W2Cl9³−; Mo2(OH)2(H2O)8⁴+ (OH− companions in Crg37). HW 1880, 1891–1893. Ordinary M–M order 3, with ordinary terminal ligands or anionic donors. W#v6#t and E#v6#t3 already occur above |
| Crg24 | Cr#u4#v2; #u4#v#h; #u4#h2; #c+2#u4#t2 | Implemented; 4 | Molecular CrH2: MCS96, pp. 873–874, Table 3, identifies the quintet ground state. Its calculated low-spin states are about 59 kcal/mol higher, not competing ground-state candidates. Ordinary, mixed/implicit H or two H− donors |
| Crg25 | Cr#n#v4#t4; #n#v3#h#t4; #n#v2#h2#t4; #n#v#h3#t4; #n#h4#t4; #c+4#n#t8 | Implemented; 6 | CrH4(dmpe)2 is explicitly diamagnetic, HW 1861–1862. Four terminal hydrides and four phosphine donors; ordinary, mixed/implicit H or all-H− donation. Distinct from the triplet tetrahalides |
| Crg26 | M#v5#h#t3; #v4#h2#t3; #v3#h3#t3; #v2#h4#t3; #v#h5#t3; #h6#t3; #c+6#t9; explicit M#v6#t3 reused | Implemented; 14 | MH6(PR3)3, HW 1871: d⁰ terminal hexahydrides. HC90, pp. 12–14, supplies WH6(PMe2Ph)3 and MoH6(PCy3)3 as fixture examples. Explicit, mixed/implicit H or H− donors. Nine accepted pairs fit capacity 18 |
| Crg27 | Mo#n3#t5#m0; W#n3#t5#m0 | Implemented; 2 | HMo2(CO)10−/HW2(CO)10−, HW 2069, 2137. Adopted the existing Cr convention: H− contribution [0,2,0] to the 3c2e bridge, with no separate M–M bond (PW). Cr is covered by D79/PW |
| Crg28a | W#c+2#n#u2#t6; Mo counterpart reused | Implemented; 1 | MoCp2/WCp2: C83, pp. 3614, 3617–3619, establishes paramagnetism and uses a triplet-derived ground manifold. Its spin–orbit doublet is not an S = ½ assignment. Cp− ligands donate six pairs in total; no metal MC state |
| Crg28b | M#c+2#n#v2#t6; #c+2#n#v#h#t6; #c+2#n#h2#t6; all-H− form M#c+4#n#t8 reused | Implemented; 6 | Cp2MH2: HW 2199–2200 gives the isolated 18-electron dihydrides; C83 uses them as photochemical precursors and reports resolved precursor NMR. Paired d² occupation in the standard closed-shell description; explicit, mixed/implicit H or H− donors |
| Crg28c | M#c+3#v2#h#t6; #c+3#v#h2#t6; #c+3#h3#t6; explicit M#c+3#v3#t6 and all-H− M#c+6#t9 reused | Implemented; 6 | Cp2MH3+: HW 2201 identifies protonation products of Cp2MH2. Ordinary terminal hydrides leave d⁰; explicit, mixed/implicit H or H− donors |
| Crg28d | Mo#n2#v2#t5; #n2#v#h#t5; #n2#h2#t5; donor Mo#c+2#n2#t7 supported in Crg30a | Implemented; 3 | MoH2(PMe3)5: Ly84, pp. 695–696, establishes two terminal hydrides and resolved ¹H/³¹P NMR down to 190 K. Paired d⁴ assignment accepted from that evidence (PW); no magnetic measurement or explicit electronic multiplicity is reported. The Crg30a donor state is also implemented |
| Crg28e | E#c-#n3#h#t5; explicit E#c-#n3#v#t5 and H−-donor E#n3#t6 reused | Implemented; 3 | HM(CO)5−, HW 2136–2137: implicit-H counterparts of the supported terminal carbonyl-hydride description. No extra donor state |
| Crg28f | Mo#n#v4#t4; #n#v3#h#t4; #n#v2#h2#t4; #n#v#h3#t4; #n#h4#t4; donor Mo#c+4#n#t8 reused | Implemented; 5 | MoH4(PR3)4, HW 1871. P71 explicitly reports diamagnetism from magnetic measurements for PMePh2, PEtPh2 and dppe complexes. HC90, p. 10, supports terminal hydrides by structure and relaxation measurements. Ordinary/mixed/implicit H or H− donors |
| Crg28g | W#n2#v2#t5; #n2#v#h#t5; #n2#h2#t5; W#n#v4#t4; #n#v3#h#t4; #n#v2#h2#t4; #n#v#h3#t4; #n#h4#t4; W#c+#n#v3#t4; #c+#n#v2#h#t4; #c+#n#v#h2#t4; #c+#n#h3#t4; #c+4#n#t7 | Open by decision; 13, excluded from supported count | WH2(PMe3)5, WH4(PMe3)4 and WH3(PMe3)4+: Ch81, pp. 1205, 1207–1209, gives resolved NMR coupling but no magnetic measurement or explicit electronic multiplicity. Paired assignments were left open. W#c+2#n2#t7 and W#c+4#n#t8 have independent cyanide support; their use for these hydrides remains open |
| Crg28h | W#c+#v5#t4; #c+#v4#h#t4; #c+#v3#h2#t4; #c+#v2#h3#t4; #c+#v#h4#t4; #c+#h5#t4; donor W#c+6#t9 reused | Implemented; 6 | WH5(PMePh2)4+, HW 1872. HC90, p. 11, explicitly identifies the pentahydride and supports its classical description with T1 = 179 ms at 240 K; ordinary hydrides leave d⁰. This is distinct from Ch81's PMe3 trihydride product formed with H2 loss |
| Crg28i | H2E(CO)5 | Unselected | HW 2137 marks these carbonyl dihydrides as not isolated. No additional default states selected |
| Crg29a | E#u5#v; #u5#h; #c+#u5#t | Implemented; 9 | CrH, MoH, WH: WA03Cr, p. 574; WA05Mo, Table 2; WA02W, Table 2. Observed monohydrides, ⁶Σ⁺ ground-state assignments; one ordinary H or H− donation |
| Crg29b | M#u4#v2; #u4#v#h; #u4#h2; #c+2#u4#t2 | Implemented; 8 | MoH2/WH2: WA05Mo and WA02W, Table 3, identify bent ⁵B2 ground states against the isotope spectra. CrH2 is already in Crg24 |
| Crg29c | W#u3#v3; #u3#v2#h; #u3#v#h2; #u3#h3; #c+3#u3#t3 | Implemented; 5 | WH3: WA02W, Table 4 and pp. 6725–6726, assigns the observed trihydride to pyramidal ⁴A1. Not the higher-energy H2 adduct |
| Crg29d | M#u2#v3#h; #u2#v2#h2; #u2#v#h3; #u2#h4; explicit M#u2#v4 and H−-donor M#c+4#u2#t4 reused | Implemented; 8 | MoH4/WH4: WA05Mo and WA02W, Table 5, give the tetrahedral triplet assignments used to interpret their observed spectra. The Cr paper instead observes H2 adducts of CrH2 |
| Crg29e | M#v5#h; #v4#h2; #v3#h3; #v2#h4; #v#h5; #h6; explicit M#v6 and H−-donor M#c+6#t6 reused | Implemented; 12 | MoH6/WH6: WA05Mo, Table 6; WA02W, Table 7. Observed hexahydrides with singlet calculations. MoH6 geometry is less decisively assigned than WH6, but both considered Mo geometries have the same occupation and ordinary-bond graph |
| Crg29f | CrH3/CrH4/CrH6; MoH3; WH5; hydride–H2 adducts | Unselected, or deferred for intact H2 | WA03Cr explicitly observes CrH/CrH2 and their H2 adducts, not the calculated higher ordinary hydrides. WA05Mo does not identify MoH3 experimentally; WA02W calculates WH5 but does not identify it in the spectra. No additions from calculation tables alone; intact-H2 cases remain with Crg33 |
| Crg30a | Mo#c-5#n2#v7; #c+2#n2#t7 | Implemented; 2 | Mo(CN)7⁵−: D77, p. 1076, explicitly concludes diamagnetism and a ¹A1′ ground state from susceptibility. Ordinary C-bound cyanides or CN− donors; distinct from doublet Mo(CN)7⁴− |
| Crg30b | W#c-5#n2#v7; #c+2#n2#t7; W#c-4#n#v7#h; #c+3#n#v#t7; #c+3#n#h#t7; explicit W#c-4#n#v8 and all-donor W#c+4#n#t8 reused | Implemented; 5 | W(CN)7⁵− and HW(CN)7⁴−, HW 1876. SG81, pp. 1886–1887, explicitly reports diamagnetism for both; NMR coupling and H/D vibrational shifts establish the W–H unit. Ordinary CN/H bonds or CN−/H− donors; the hydride's ionic allocation leaves paired d² |
| Crg31a | EF7³− | Open occupation/coordination; not counted | HW 1855, 1872. The formulas alone do not establish a discrete seven-coordinate ion or its d² occupation. No specific primary structure/spin source located; no transfer from octahedral EX6²− |
| Crg31b | Cr(O2)2(CN)3³− | Open occupation; not counted | HW 1856; S65, p. 71, establishes seven-coordinate, pentagonal-bipyramidal Cr(IV), two peroxide groups and three C-bound cyanides. Neither S65 nor B77, pp. 1295–1297, reports magnetic data or an electronic occupation. The longer *Ark. Kemi* reports cited by B77 remain unlocated. No transfer of the ammonia-peroxide triplet |
| Crg32 | Mo#u2#v4#t3; #c+3#u2#v#t6 | Implemented single-bond graph; 2 | Mo2Cl9³−: HW 1893 gives a Mo–Mo single bond. MSL97, p. 6270, describes the experimentally relevant 2.6–2.8 Å region by one σ bond and four weakly coupled δπ electrons, consistent with paramagnetism. This gives two unpaired electrons per Mo, with ordinary terminal Cl/bridging Cl− donors or all-Cl− donation. The paper does not assert this partition for every salt; molecular antiferromagnetic coupling does not replace the local #u2 by #n |
| Crg33 | η²-H2 complexes, including W(H2)(CO)3(PR3)2 | Deferred with the earlier H2-donor accounting work | HW 1872, 2071–2073. The classical two-hydride form and the intact-H2 form are distinct; no aggregate MC replacement |
| Crg34 | Extended binary oxides/halides/sulfides, carbides, nitrides, borides; highly charged nitridometallates; larger polyoxo-, sulfido- and halide clusters | Extracted; no additional default states selected | HW 1860, 1863, 1874–1890. Retain the distinction between bulk formulas and finite molecular graphs. The low-spin d³ CrN3⁶− statement on p. 1863 is recorded here, not used to replace the ordinary Cr(III) quartet. No broad cluster MC states |
| Crg35 | CrF6; MoCl6; WMe4/WMe8²−; Cr(O2)4²−; specialized ligand rearrangements and higher assemblies | Unselected | HW itself questions the first four assignments (1855, 1872, 1894); Cr(O2)4²− is only a postulated transient, 1854–1855. No extra states from exotic scaffolds or hypothetical species |
| Crg36 | Cr#u3#v3; #c+3#u3#t3 | Implemented; 2 | Cr(NiPr2)3: A68, p. 495, gives μeff = 3.80 BM, temperature-independent over 298–123 K. BHN71, p. 412, identifies the zero-field-split d³ ESR, g∥ = 2 and g⊥ = 4 at 130 K. These support the quartet; the room-temperature g = 1.975 signal is an oxygen-complex impurity |
| Crg37 | O#c-#n#v#d2; O#c-#n#h#d2 | Implemented; 2 O companions | μ₂-hydroxido donors with explicit/implicit H, HW 1860, 1880, 1893; adopted donor allocation (PW). One pair donated to each metal, one pair retained. Both Mo2(OH)2(H2O)8⁴+ fixtures infer the O states and pass the electron-count validator in all four resolver/policy combinations; O capacity 8 unchanged |

#### Chromium-group limits and coverage

Expansion under the registry's concrete zero defaults and deduplication against
the 11 earlier states gives 260 implemented new metal states. All satisfy the
nonbonding budget and capacity 18. The two O companions fit capacity 8.
The ValenceTable is unchanged: Counts uses the table shifted by `2t − charge` and
the actual element's capacity. Supplied #u
preserves unpaired occupations when Counts' default pairing would differ.
Crg27 adds the only two metal MC states; no aromatic metal state was added.

All 260 new metal states and both O companions have inferred coverage in 262
independent fixtures: 170 `donor_acceptor`, 41 `hydrides`, 32 `inorganic_ions`,
11 `atomic_ions` and 8 `inorganic_small`. Their inputs omit target-atom lone pairs,
localized valence and accepted-pair counts. AtomTyping and Counts, each with Strict and
MostSaturated, infer the expected states and agree on the complete output graphs.
The independent audit checks 288 target atoms across those four configurations
(1,152 checks), including electron conservation. Crg37's explicit-/implicit-H
fixtures also cover both Mo atoms; neither supplies O lone pairs or donation counts.

Verification: 3,594 graph unit tests, 223 graph property tests and 2,238 resolution
conformance cases pass. All 1,976 pre-chapter snapshots are unchanged. The
existing-entry citation/comment audit is complete; the 13 open W occupations,
Crg31 source checks and intact-H2 accounting remain outside this selection.

#### Chromium-group primary sources

| Key | Source | Access and use |
| --- | --- | --- |
| MCS96 | Ma, Collins and Schaefer, *J. Am. Chem. Soc.* **118** (1996), 870–879. [DOI](https://doi.org/10.1021/ja951376t), [PDF](../materials/valence/Ma-Collins-Schaefer-1996-TiH2.pdf) | Archived paper; CrH2 ground-state discussion and Table 3 reviewed. Crg24 |
| G85 | Girolami et al., *J. Chem. Soc., Dalton Trans.* (1985), 1339–1348. [DOI](https://doi.org/10.1039/DT9850001339), [PDF](../materials/valence/Girolami-1985-divalent-dmpe-complexes.pdf) | Archived paper; CrCl2(dmpe)2/CrMe2(dmpe)2 have two unpaired electrons, p. 1340. Crg19 |
| KG75 | Kiernan and Griffith, *J. Chem. Soc., Dalton Trans.* (1975), 2489–2494. [DOI](https://doi.org/10.1039/DT9750002489), [PDF](../materials/valence/Kiernan-Griffith-1975-octacyanoniobates.pdf) | Archived paper; octacyanide oxidation states/occupations. Crg14 |
| R71 | Rossman, Caltech PhD thesis (1971), part I, especially p. 51. [Repository](https://thesis.caltech.edu/9999/), [PDF](../materials/valence/Rossman-1971-heptacyanomolybdate-thesis.pdf) | Downloaded from Caltech and archived; explicit Mo(CN)7⁴− S = ½, with susceptibility/EPR. Crg20 |
| MDP19 | Magott, Dunbar and Pinkowicz, “Correlating magnetic anisotropy with [Mo(CN)7]4− geometry of MnII–MoIII magnetic frameworks,” *Dalton Trans.* (2019). [DOI](https://doi.org/10.1039/C9DT02164G), [PDF](../materials/valence/Magott-Dunbar-Pinkowicz-2019-heptacyanomolybdate.pdf) | Downloaded from NSF PAR and archived; explicitly assigns S = ½ to Mo(CN)7⁴−. The framework itself is not a proposed fixture |
| SDV25 | Sauza-de la Vega et al., “Designing molecular qubits: computational insights into first-row and group 6 transition metal complexes,” *Chem. Sci.* **16** (2025), 12896–12905. [DOI](https://doi.org/10.1039/D5SC02544C), [archived full text](../materials/valence/Sauza-de-la-Vega-2025-group6-qubits.html) | Full text reviewed via PMC and archived as HTML; Table 2 explicitly assigns calculated triplet Cr/Mo/W tetraaryl ground states. PDF retrieval failed; the HTML contains the evidence |
| D79 | Darensbourg and Incorvia, *Inorg. Chem.* **18** (1979), 18–22. [DOI](https://doi.org/10.1021/ic50191a004), [PDF](../materials/valence/Darensbourg-1979-bridging-carbonyl-hydrides.pdf) | Existing archived source for the accepted chromium hydride-bridge convention; HW explicitly extends the structural class to Mo/W. Crg27 |
| S98 | Stranger, McGrady and Lovell, *Inorg. Chem.* **37** (1998), 6795–6806. [DOI](https://doi.org/10.1021/ic971511y), [PDF](../materials/valence/Stranger-McGrady-Lovell-1998-metal-metal-bonding.pdf) | Existing archived paper; d³d³ comparison and reference to MSL97. Crg32 |
| CND01 | Cage, Nguyen and Dalal, “Three-dimensional antiferromagnetism in triamminodiperoxychromate, Cr(NH3)3(O2)2,” *Solid State Commun.* **119** (2001), 597–601. [DOI](https://doi.org/10.1016/S0038-1098(01)00296-4), [PDF](../materials/valence/Cage-Nguyen-Dalal-2001-CrIV-peroxo-magnetism.pdf) | Supplied PDF reviewed and archived. Susceptibility and magnetic entropy independently establish S = 1, pp. 599–600. Crg15 |
| C83 | Cox et al., “Electronic structure of molybdenocene and tungstenocene: detection of paramagnetism by magnetic circular dichroism in argon matrixes,” *Inorg. Chem.* **22** (1983), 3614–3620. [DOI](https://doi.org/10.1021/ic00166a025), [PDF](../materials/valence/Cox-1983-molybdenocene-tungstenocene.pdf) | Reviewed and archived. MCD establishes paramagnetism; the analysis uses a triplet manifold with spin–orbit splitting, pp. 3617–3619. The spin–orbit doublet does not imply one unpaired electron. Crg28a–b |
| Ly84 | Lyons et al., “Synthesis, X-ray crystal structure, and reactions of dihydridopentakis(trimethylphosphine)molybdenum(II),” *J. Chem. Soc., Dalton Trans.* (1984), 695–700. [DOI](https://doi.org/10.1039/DT9840000695), [PDF](../materials/valence/Lyons-1984-molybdenum-phosphine-dihydride.pdf) | Reviewed and archived. Two terminal hydrides and pentagonal-bipyramidal geometry by diffraction; resolved ¹H/³¹P NMR, pp. 695–696. Paired occupation accepted as PW; the paper gives no magnetic measurement or explicit electronic multiplicity. Crg28d |
| D77 | Drew, Mitchell and Pygall, “Preparation and properties of the heptacyanomolybdate(II) anion and the crystal and molecular structures of its sodium and potassium salts,” *J. Chem. Soc., Dalton Trans.* (1977), 1071–1077. [DOI](https://doi.org/10.1039/DT9770001071), [PDF](../materials/valence/Drew-Mitchell-Pygall-1977-heptacyanomolybdate-II.pdf) | Supplied PDF reviewed and archived. Explicit diamagnetism and ¹A1′ assignment for Mo(CN)7⁵−, p. 1076. Crg30a |
| WA03Cr | Wang and Andrews, “Chromium Hydrides and Dihydrogen Complexes in Solid Neon, Argon, and Hydrogen,” *J. Phys. Chem. A* **107** (2003), 570–578. [DOI](https://doi.org/10.1021/jp026930h), [PDF](../materials/valence/Wang-Andrews-2003-chromium-hydrides.pdf) | Supplied PDF reviewed and archived. CrH ⁶Σ⁺ and CrH2 ⁵B2, Table 2 and pp. 574–575. Observed higher-H species are H2 adducts; no observed CrH3 or ordinary CrH4 assignment. Crg24, Crg29a/f |
| WA05Mo | Wang and Andrews, “Matrix Infrared Spectra and Density Functional Theory Calculations of Molybdenum Hydrides,” *J. Phys. Chem. A* **109** (2005), 9021–9027. [DOI](https://doi.org/10.1021/jp053591u), [PDF](../materials/valence/Wang-Andrews-2005-molybdenum-hydrides.pdf) | Supplied PDF reviewed and archived. Observed MoH/MoH2/MoH4/MoH6 interpreted with sextet/quintet/triplet/singlet calculations, Tables 2, 3, 5, 6. MoH3 calculations alone do not establish an observed species. Crg29 |
| WA02W | Wang and Andrews, “Neon Matrix Infrared Spectra and DFT Calculations of Tungsten Hydrides WHx (x = 1–4, 6),” *J. Phys. Chem. A* **106** (2002), 6720–6729. [DOI](https://doi.org/10.1021/jp025920d), [PDF](../materials/valence/Wang-Andrews-2002-tungsten-hydrides.pdf) | Supplied PDF reviewed and archived. Observed WH/WH2/WH3/WH4/WH6 interpreted with sextet/quintet/quartet/triplet/singlet calculations, Tables 2–5, 7. WH5 is calculated but not identified experimentally. Crg29 |
| BHN71 | Bradley, Hursthouse and Newing, “Crystal structure and electron spin resonance spectrum of tris-di-isopropylaminatochromium(III),” *J. Chem. Soc. D* (1971), 411–412. [DOI](https://doi.org/10.1039/C2971000411B), [PDF](../materials/valence/Bradley-Hursthouse-Newing-1971-chromium-trisamide.pdf) | Supplied PDF reviewed and archived. Trigonal structure and zero-field-split d³ ESR, g∥ = 2, g⊥ = 4 at 130 K, p. 412; the room-temperature signal is an oxygen impurity. Crg36's quartet is the interpretation of these ESR data |
| B77 | Baran, “Die IR-Spektren von [Cr(O2)2(NH3)3] und K3[Cr(O2)2(CN)3],” *Monatsh. Chem.* **108** (1977), 1295–1299. [DOI](https://doi.org/10.1007/BF01046442), [PDF](../materials/valence/Baran-1977-chromium-peroxides.pdf) | Retrieved from the author's university repository, reviewed and archived. Structure, C-bound cyanide and Cr(IV) supported; no spin assignment. Crg31b |
| SG81 | Soares and Griffith, “Studies on transition-metal cyano-complexes. Part 3. Heptacyano-complexes of tungsten, molybdenum, and vanadium, and a hydrido-heptacyano-complex of tungsten,” *J. Chem. Soc., Dalton Trans.* (1981), 1886–1891. [DOI](https://doi.org/10.1039/DT9810001886), [PDF](../materials/valence/Soares-Griffith-1981-heptacyanides.pdf) | Reviewed and archived. Explicit diamagnetism for both W salts, p. 1886; W–H NMR coupling and isotope-sensitive vibrations, p. 1887. Crg30b |
| MSL97 | McGrady, Stranger and Lovell, “Broken-symmetry and approximate spin-projected potential energy curves for bimetallic systems,” *J. Phys. Chem. A* **101** (1997), 6265–6272. [DOI](https://doi.org/10.1021/jp9625568), [PDF](../materials/valence/McGrady-Stranger-Lovell-1997-d3-dimers.pdf) | Reviewed and archived. Single σ bond with four weakly coupled δπ electrons in the longer-bonded Mo2Cl9³− salts, p. 6270; distinct from the computed short-distance triple-bond regime. Crg32 |
| A68 | Alyea, Basi, Bradley and Chisholm, “Low co-ordination numbers for chromium in chromium(III) compounds,” *Chem. Commun.* (1968), 495. [DOI](https://doi.org/10.1039/C19680000495), [PDF](../materials/valence/Alyea-1968-chromium-amides.pdf) | Reviewed and archived. Cr(NiPr2)3 μeff = 3.80 BM, constant over 298–123 K. Crg36 |
| P71 | Pennella, “Tetrahydrido-complexes of molybdenum,” *Chem. Commun.* (1971), 158. [DOI](https://doi.org/10.1039/C2971000158A), [PDF](../materials/valence/Pennella-1971-molybdenum-tetrahydrides.pdf) | Reviewed and archived. Explicit diamagnetism from magnetic measurements; four hydrides from NMR integration and H/D vibrational shifts. Crg28f |
| Ch81 | Chiu et al., “Reactions of hexamethyltungsten(VI) in the presence of trimethylphosphine,” *J. Chem. Soc., Dalton Trans.* (1981), 1204–1211. [DOI](https://doi.org/10.1039/DT9810001204), [PDF](../materials/valence/Chiu-1981-tungsten-phosphine-hydrides.pdf) | Reviewed and archived. WH2(PMe3)5, WH4(PMe3)4 and WH3(PMe3)4+; resolved NMR but no explicit electronic multiplicity. Protonation gives the trihydride plus H2, pp. 1208–1209. Occupations left open in Crg28g |
| HC90 | Hlatky and Crabtree, “Molybdenum and tungsten phosphine polyhydrides,” *Inorg. Synth.* **27** (1990), 8–14. [DOI](https://doi.org/10.1002/9780470132586.ch2), [PDF](../materials/valence/Hlatky-Crabtree-1990-Mo-W-phosphine-polyhydrides.pdf) | Retrieved from a university-hosted copy; complete chapter reviewed and archived. Mo tetrahydride structure/relaxation, p. 10; classical WH5(PMePh2)4+ and its relaxation evidence, p. 11; WH6(PMe2Ph)3 and MoH6(PCy3)3, pp. 12–14. Crg26, Crg28f/h |
| S65 | Stomberg, “Crystal structure of potassium diperoxotricyanochromate(IV),” *Nature* **205** (1965), 71. [DOI](https://doi.org/10.1038/205071a0), [PDF](../materials/valence/Stomberg-1965-potassium-diperoxotricyanochromate.pdf) | Supplied PDF reviewed and archived. Seven-coordinate, pentagonal-bipyramidal structure; O–O 1.440/1.455 Å and Cr–C 2.081–2.102 Å. No magnetic data or electronic occupation reported. Crg31b |

All supplied PDFs are reviewed and archived. B77 and HC90 were additionally
obtained from university repositories. Crg31 records the outstanding
structure/occupation evidence. Crg28g remains open by decision.
Intermediates and the count audit are in `scratch/hw-chromium-review/`; retained
source material is in `materials/valence/`.

### Manganese group

HW II, chapter XXVIII, pp. 1899–1933 (chemistry through p. 1931), with the cited
chapter-XXXII cross-references. Selected scope implemented 2026-11-04:
**243 additions: Mn 62, Tc 60, Re 121**, covered by 242 fixtures. Mn capacity
18 → 21 is implemented; ValenceTable is unchanged. All nine earlier Mn/Tc/Re/Bh
states are retained with checked citations.

#### Manganese-group work table

E = Mn, Tc, Re; M = Tc, Re. A prefix applies to every following suffix in its cell.
Only new states are listed after reuse across rows; counts include all displayed
explicit/mixed/implicit-H forms. Entries use the concrete, zero-valence registry
defaults. “Implemented” records the selection approved and added on 2026-11-04. Ordinary oxo bonds and anionic-ligand donation follow the preceding
chapters' graph conventions.

| ID | States / extracted family | Status; new states | Examples and evidence |
| --- | --- | --- | --- |
| Mng01 | E#n#u5; Bh#n#u5; Mn#c+2#u5; Mn#c+2#n2#u; Mn#c+3#u4; Mn#c+3#n#u2; Mn#u3#v4 | Retained; 9 existing | Atomic occupations: HW II, p. 1538; Mn(II/III) occupations, pp. 1906–1907, 1910; Mn(IV), 1907, 1914–1915. The neutral #v4 comment now identifies MnF4/MnR4 |
| Mng02 | E#c+4#u3; #c+5#u2; #c+6#u; #c+7 | Implemented; 12 | IV–VII ionic occupations, d³/d²/d¹/d⁰. HW 1901, 1907, 1911–1913, 1917, 1921–1922. Free-ion occupations do not determine those of coordinated complexes. |
| Mng03 | M#c+2#u5; #c+3#u4 | Implemented; 4 | II/III free-ion occupations, d⁵/d⁴, for the chapter's lower-valent chemistry. HW 1917, 1923. The corresponding Mn states already exist. |
| Mng04 | Mn#u5#v2; #c+2#u5#t2; #c-#u5#v3; #c+2#u5#t3; #c-2#u5#v4; #c+2#u5#t4 | Implemented; 6 | Monomeric MnR₂, [MnMes₃]−, [MnMe₄]²−; MnX₄²− and tetrahedral high-spin Mn(CN)₄²−. HW 1903–1904, 1906–1907, 1913–1914. Ordinary bonds or anionic-ligand donors. |
| Mng05 | Mn#u5#v2#t2; #u5#v2#t4; #c-#u5#v3#t3; #c+2#u5#t6; #c-4#u5#v6 | Implemented; 5 | MnX₂L₂, MnCl₂·4H₂O, MnBr₂(dmpe)₂, [Mn(OH)₃(H₂O)₃]−, Mn(H₂O)₆²+/Mn(NH₃)₆²+, MnCl₆⁴−. HW 1903, 1905–1907; G85, pp. 1340–1341, explicitly establishes high-spin MnBr₂(dmpe)₂/MnI₂(dmpe)₂. |
| Mng06 | Mn#c+2#u5#t7; #c-2#u5#v4#t3; #c+2#u5#t8; #c-6#u5#v8 | Implemented; 4 | [Mn(edta)(H₂O)]²− and [Mn(NO₃)₄]²−: seven/eight donors, or ordinary anionic O bonds with the other donors retained. HW 1902, 1906–1907. Counts electron counts 19/21 use the implemented Mn capacity 21. |
| Mng07 | Mn#c+2#n2#u#t6; #c-4#n2#u#v6; #n2#u#v2#t4 | Implemented; 3 | Low-spin Mn(CN)₆⁴−, Mn(CNR)₆²+, MnMe₂(dmpe)₂ and MnCp*₂; HW 1907, 1914; G85, p. 1343, explicitly distinguishes the low-spin methyl from the high-spin halides. The haptic Cp graph reuses the six-accepted-pair state. |
| Mng08 | Mn#c+3#u4#t5; #c-2#u4#v5; #c+3#u4#t6; #u4#v3#t3; #c-3#u4#v6; #c-#u4#v4; #c+3#u4#t4 | Implemented; 7 | MnCl₅²−, Mn(H₂O)₆³+, Mn(acac)₃, Mn(OH)₆³− and explicitly high-spin MnMe₄−/MnMe₅²−. HW 1907–1910, 1914. Anionic donors or ordinary ligand bonds; acac retains its neutral O donor. |
| Mng09 | Seven-coordinate Mn(NO₃)₃(bipy) | No addition; uncertain occupation | HW 1902 gives the geometry; 1910 discusses predominantly high-spin Mn(III) but does not assign this particular complex. No prospective coordination/occupation state is selected for this specific mixed-ligand example |
| Mng10 | Mn#c-3#n#u2#v6; #c+3#n#u2#t6 | Implemented; 2 | Low-spin Mn(CN)₆³−, HW 1907–1908, 1910. Paired plus two unpaired d electrons, distinct from the high-spin III family. |
| Mng11 | Mn#c+4#u3#t4; #c-2#u3#v6; #c+4#u3#t6; #u3#v4#t2; #c+#u3#v3#t3 | Implemented; 5 | MnF₄/MnR₄; MnF₆²−, MnCl₆²−, Mn(CN)₆²−, MnMe₆²−; MnMe₄(PMe₃)₂ and CpMnR₃. HW 1907–1908, 1914–1915. Existing Mn#v4#u3 supplies the neutral four-bond state. |
| Mng12 | Mn#u2#v5; #c-3#u2#v8; #u#v6; #c-2#u#v8; #c+2#u#v4#t2 | Implemented; 5 | MnO₄³−/MnO₄²−, explicitly paramagnetic d²/d¹, and MnO₂Cl₂. HW 1904, 1907, 1911–1912. Ordinary oxo graphs, including all-double-bond charge placement; the last state retains Mn=O and makes the chlorides donors. MnOCl₃ reuses #u2#v5; see Mng30a. |
| Mng13 | E#v7; #c-#v8; #c+#v6#t | Implemented; 9 | EO₄−, E₂O₇ and EO₃X (X = F/Cl), HW 1912–1913, 1921, 1924–1925. Ordinary E=O; the last state makes X− a donor. C₆F₅Mn(NtBu)₃ and MeMO₃ reuse #v7 (1915, 1930). |
| Mng13a | Tc#c+#v6 | Implemented; 1 | TcO₃+: F82, pp. 5304, 5306, characterizes the d⁰ trioxocation by ⁹⁹Tc/¹⁷O NMR in HF and reports TcO₃AsF₆ stable to pumping at room temperature. Three ordinary Tc=O bonds. |
| Mng14 | E#c+#n3#t6; #n3#v#t5; #n3#h#t5; #c-#n4#t5; #c-3#n5#t4; #c-5#n3#v6 | Implemented; 18 | E₂(CO)₁₀, EH(CO)₅, E(CO)₅−, E(CO)₄³− and E(CN)₆⁵−; HW 1902–1903, 1914, 1923–1924, 2111–2113, 2135–2137. GKB78, pp. 1414–1415, confirms the Re(CN)₆⁵− formulation. The six-pair cation is also the CN−/H−-donor allocation. ReCl(PMe₃)₅ reuses #n3#v#t5 (1923). |
| Mng14b | Re#c-#n3#v2#t4; #c-#n3#v#h#t4; #c-#n3#h2#t4 | Implemented; 3 | ReH₂(CO)₄−, HW 2136–2137: ordinary terminal H, all explicit/mixed/implicit forms. HRe₂(CO)₉− reuses the explicit and mixed-H states with an ordinary Re–Re bond; H−-donor descriptions reuse Mng14. |
| Mng15 | M#u#v6; #c-#u#v7; #c-2#u#v8; #c+6#u#t6; #c+6#u#t7; #c+6#u#t8 | Implemented; 12 | MF₆, MF₇− and MF₈²−, d¹; HW 1921. ReCl₆, ReMe₆ and ReMe₈²− reuse these states (1921, 1930). TcCl₆ is not used as evidence. |
| Mng16 | M#c-2#u3#v6; #c+4#u3#t6; #u3#v4#t2 | Implemented; 6 | Octahedral MX₆²− and MX₄(PR₃)₂, d³; HW 1922. The eight-coordinate diars case remains outside this selection. |
| Mng17 | M#c-#u2#v6; #c+5#u2#t6 | Implemented; 4 | Octahedral MF₆−, d²; HW 1922 and its octahedral ligand-field treatment. No spin-paired alternative is inferred; the polymeric/dimeric pentahalides are treated in Mng34/Mng34a. |
| Mng18 | Re#c+#v6; #c+7#t6; #c+7#t7; #c+7#t8 | Implemented; 4 | ReF₆+, ReF₇ and ReF₈−, HW 1921. Ordinary F bonds or F− donors; neutral #v7 and negative #v8 are reused. No MC alternative. |
| Mng19 | Re#c+5#v2#t5; #c+5#v2#t6; #c+3#v4#t3; #c+#v6#t2; #c+3#v4#t4; #c+#v6#t3 | Implemented; 6 | ReOF₅/ReOF₆−, molecular ReO₂Cl₃ and its chloride-bridged dimer, ReO₃Cl₃²−; HW 1921. ReO₃(pmdta)+ and CpReO₃ reuse the last state (1925, 1930). Oxo bonds stay ordinary; bridging Cl− uses the existing donor convention. |
| Mng20 | Re#c+4#u#v2#t4; #c+4#u#v2#t5 | Implemented; 2 | ReOCl₄/ReOBr₄ and ReOX₅−, d¹; HW 1921. Ordinary halide descriptions reuse #u#v6 and #c-#u#v7. |
| Mng21 | Re#v7#t; #v7#t2; #c+2#v5#t3; #c+4#v3#t5 | Implemented; 4 | Re₂O₇(D)₂, MeReO₃(D)/(D)₂, methylrhenium mono-/diperoxo adducts and O[ReO(O₂)₂(H₂O)]₂; HW 1924–1925, 1930–1931. S15, Fig. 2 (5), supplies the concrete 4,4′-dimethyl-2,2′-bipyridine adduct for #v7#t2. Ordinary peroxo bonds or peroxide-O− donors; O–O bonds retained. No cluster-specific MC state. |
| Mng22 | M#c-2#v9; #c-2#v8#h; #c-2#v7#h2; #c-2#v6#h3; #c-2#v5#h4; #c-2#v4#h5; #c-2#v3#h6; #c-2#v2#h7; #c-2#v#h8; #c-2#h9; #c+7#t9 | Implemented; 22 | MH₉²−, explicitly diamagnetic with nine terminal hydrides; HW 1919–1920. All explicit/mixed/implicit-H forms and H− donors. #c-2#v9 also covers O₃M≡N²− (1927); ReO₃Cl₃²− reuses it. Maximum electron count 18. |
| Mng23 | Re#c-#v8#t; #c-#v7#h#t; #c-#v6#h2#t; #c-#v5#h3#t; #c-#v4#h4#t; #c-#v3#h5#t; #c-#v2#h6#t; #c-#v#h7#t; #c-#h8#t; #v6#h#t2; #v5#h2#t2; #v4#h3#t2; #v3#h4#t2; #v2#h5#t2; #v#h6#t2; #h7#t2 | Implemented; 16 | ReH₈(PR₃)− and classical ReH₇(PR₃)₂ with alkyl phosphines; HW 1919–1920. All H-count splits; H−-donor forms reuse Re#c+7#t9. The aryl-phosphine ReH₅(H₂)(PR₃)₂ description is excluded from this row. |
| Mng24 | Re#n#v5#t3; #n#v4#h#t3; #n#v3#h2#t3; #n#v2#h3#t3; #n#v#h4#t3; #n#h5#t3; #c+5#n#t8; #n2#v3#t4; #n2#v2#h#t4; #n2#v#h2#t4; #n2#h3#t4; #c+3#n2#t7 | Implemented; 12 | Classical ReH₅(PR₃)₃ and ReH₃(PR₃)₄; HW 1919 and the closed-shell 18-electron statement on 1920. The latter page misprints ReH₄(PR₃)₃. RJ23, pp. 3, 9–10, supports the classical hydrides: trihydride H positions are inferred from the Tc analogue, pentahydride H positions located and refined. ReH(PR₃)₅ reuses Mng14. |
| Mng25 | Re#c-4#n2#v7 | Implemented; 1 | Diamagnetic, low-spin d⁴ Re(CN)₇⁴−, HW 1923 and explicitly GKB78, p. 1412; the CN−-donor state Re#c+3#n2#t7 is reused. The Tc homologue is supported separately in Mng32. |
| Mng26 | Tc#c+3#v4#t4; Tc#c-#v8#t; Tc#c+3#v4#t5; Re#c+3#v4#t5 | Implemented; 4 | M₂X₈²−, including an axial donor per M; explicit quadruple M–M bond. HW 1922–1923 and 1928, including p. 1922 footnote 6. Re₂Me₈²− reuses the unadducted states (1929). |
| Mng27 | Reuses preceding rows | Implemented; 0 | Binary Re₃X₉ trimers and Re₃X₉L₃; HW 1922, 1928. Two ordinary Re=Re bonds per metal and two bridging X− donors; no aggregate MC encoding. The gas-phase trimer and adduct structures are explicitly distinguished from the extended solid. |
| Mng28 | Re#c+4#v3#t6; #n#v5#t2; #c+2#n#v3#t4 | Implemented; 3 | Re₂Cl₉− and Re₂Cl₄(PR₃)₄ with explicit Re≡Re bonds; HW 1922–1923, 1928, Table 28.3. For the Re(II) dimer, #n is the residual paired two-electron budget after its triple bond. Bridging Cl− uses existing types. |
| Mng29a | Mn#u6#v; #u6#h; #c+#u6#t | Implemented; 3 | MnH, ⁷Σ⁺: WA03MnRe, pp. 4086–4087, matrix IR and comparison with gas-phase spectroscopy. Explicit/implicit H and H− donor. |
| Mng29b | Mn#u5#v#h; #u5#h2 | Implemented; 2 | MnH₂, sextet: WA03MnRe, p. 4087, isotope-resolved IR agrees with CCSD(T)/DFT frequencies. Geometry differences do not change the occupation. Explicit-H and H−-donor states reuse Mng04; these are the mixed/implicit-H forms. |
| Mng29c | Mn#c-#n#u4#v2; #c-#n#u4#v#h; #c-#n#u4#h2; #c+#n#u4#t2 | Implemented; 4 | MnH₂−, quintet: WA03MnRe, pp. 4085, 4087–4088, Table 3 and isotope-resolved IR assignment, also tied to prior gas-phase photodetachment. All H splits and H− donors. |
| Mng29d | Re#n#u4#v; #n#u4#h; #c+#n#u4#t | Implemented; 3 | ReH, ⁵Σ⁺: WA03MnRe, p. 4088, both DFT methods favor this over the septet and its frequency agrees with matrix IR. Explicit/implicit H and H− donor. |
| Mng29e | Re#u5#v2; #u5#v#h; #u5#h2; #c+2#u5#t2 | Implemented; 4 | ReH₂, sextet ⁶A₁: WA03MnRe, p. 4088. Both functionals' sextet frequencies fit the isotope-resolved spectrum; the paper explicitly rejects the nearby calculated quartet's frequencies. All H splits and H− donors; no quartet proposal. |
| Mng29f | Re#n#u#v4; #n#u#v3#h; #n#u#v2#h2; #n#u#v#h3; #n#u#h4; #c+4#n#u#t4 | Implemented; 6 | ReH₄, doublet ²A₁: WA03MnRe, pp. 4086, 4088–4089, Table 6 and reaction 5; isotope-resolved stretching and bending bands support the assignment. All H splits and H− donors. ReH₃ is calculated but not identified and supplies no proposal. |
| Mng30 | CpMnO₂ | No addition; uncertain occupation | HW 1915. No reviewed primary occupation assignment; no prospective state is selected. MnOCl₃ and the terminal nitrido complexes are resolved below |
| Mng30a | Mn#c+3#u2#v2#t3 | Implemented; 1 | MnOCl₃: JH75, p. 1269, explicitly assigns the triplet ³A₂ ground state; pp. 1270–1271 compare the calculated transitions with its optical spectrum. Ordinary Mn=O and three Cl− donors. The all-ordinary graph reuses Mn#u2#v5 in Mng12; HW 1907. |
| Mng30b | Mn#n#v5#t2; #c+2#n#v3#t4 | Implemented; 2 | Mn(salen)N and saldpen analogues: C98, p. 3108, explicitly assigns the paired ¹A′ ground state. Ordinary Mn≡N; two imine donors and either ordinary Mn–O bonds or two phenolate-O− donors. HW 1911's terminal nitrido-Mn(V) class. |
| Mng31 | Re#c+3#n#v2#t5; #c+#n#v4#t4 | Implemented; 2 | Paired Re(V): ReOCl₃(PPh₃)₂ and ReNCl₂(PPh₃)₂, FL71, p. 2691, Table 2; trans-ReO₂(en)₂+, GW66, p. 1227. Ordinary Re=O/Re≡N; ordinary or donating Cl−. The ordinary #n#v5#t2 and nitrido donor #c+2#n#v3#t4 states reuse Mng28. |
| Mng31a | Re#c-2#n#v7; #c-#n#v6#t | Implemented; 2 | ReOCl₅²− and ReOCl₄(PPh₃)−: FL71, pp. 2689–2691, establishes paired d²; earlier large ReOCl₅²− moments arose from ReCl₆²− contamination. Chloride-donor graphs reuse Mng31. No assignment to bare ReOCl₄− is inferred from the phosphine adduct. |
| Mng31b | ReNCl₄²−; remaining Tc/Re(V) oxocyanides | Fixture candidates; 0 new states | HW 1924, 1927–1928: ReNCl₄²−, ReO(OH)(CN)₄²− and Re₂O₃(CN)₈⁴−. TJD80/T79 establish TcO(CN)₅²− and TcO(OMe)(CN)₄²− without an explicit occupation assignment. GKB78 confirms the Re oxo-bridged dimer's formulation without assigning its occupation. Paired fixture descriptions may use analogy with the accepted oxo/nitrido assignments, reusing Mng28/Mng31/Mng31a/Mng31c. No prospective occupation states. The dioxocyanides and ReOX₄− are resolved in Mng31e–g |
| Mng31c | Tc#c-#n#v6; #c+3#n#v2#t4; #c-2#n#v7; #c+3#n#v2#t5 | Implemented; 4 | TcOCl₄− and TcOCl₅²−, paired d²: T79, pp. 48, 74, 87–88, reports the former's susceptibility and summarizes diamagnetism of the latter. J13, pp. 1098–1101, establishes their five-/six-coordinate structures and analyzes the electronic levels. Ordinary Tc=O; ordinary or donating Cl−. |
| Mng31d | Tc#c+2#n#v3#t4; #c+3#u#v3#t4; #c+4#u#v2#t4 | Implemented; 3 | TcNX₄²− (paired d²), TcNX₄− (d¹), X = Cl/Br, and TcOCl₄ (d¹): B98, pp. 2303–2307, explicitly assigns the occupations using spectroelectrochemistry and calculations. The reduced nitrido ions and neutral TcOCl₄ are characterized in solution; TcNX₄− salts are isolated. Ordinary graphs reuse Mng31c/Mng15. No neutral TcOBr₄ proposal. |
| Mng31e | Re#c-3#n#v8 | Implemented; 1 | trans-ReO₂(CN)₄³−: W84, pp. 83, 148, 152, gives the measured compound and explicitly assigns the diamagnetic ¹A₁g[(b₂g)²] ground state in its spectroscopic analysis. Two ordinary Re=O and four Re–CN bonds; the CN−-donor graph reuses Re#c+#n#v4#t4 from Mng31. HW 1924. |
| Mng31f | Tc#c-3#n#v8; #c+#n#v4#t4 | Implemented; 2 | trans-TcO₂(CN)₄³−: F82, pp. 5303, 5305, explicitly identifies diamagnetic d² Tc(V) and reports its ⁹⁹Tc NMR spectrum. Two ordinary Tc=O bonds with either four ordinary Tc–CN bonds or four CN− donors; TJD80 supplies the structure. Both allocations count 18 electrons. |
| Mng31g | Re#c-#n#v6; #c+3#n#v2#t4 | Implemented; 2 | ReOX₄−, X = Cl/Br/I: CL66, pp. 9–11, characterizes the salts; p. 14 explicitly assigns paired d² to the square-pyramidal family, based on the measured susceptibility of the bromide. Ordinary Re=O with ordinary halide bonds or X− donors; both count 14 electrons. The solvent adducts reuse Mng31/Mng31a. |
| Mng32 | Tc#c-4#n2#v7; #c+3#n2#t7 | Implemented; 2 | Tc(CN)₇⁴−: TJD80, pp. 1994–1995, establishes the isolated salt, pentagonal-bipyramidal structure by IR/Raman and diamagnetic d⁴ singlet ground state by susceptibility. Resolves HW 1923's tentative existence. Ordinary CN bonds or CN− donors. |
| Mng33 | Cp₂TcH | Fixture candidate; 0 new states | HW 1928–1930; FS67, p. 93, establishes a monomeric terminal hydride by IR, NMR and mass spectroscopy, with a proposed bent structure analogous to Cp₂ReH. No electronic occupation is explicitly assigned. A paired d⁴ fixture by analogy with LH76's Cp₂ReH assignment can use H− donation and accepted Tc#c+3#n2#t7 (Mng32). The additional localized-H states Tc#c+2#n2#v#t6 and Tc#c+2#n2#h#t6 are not selected. FN59's anomalous TcI₂(diars)₂ moment remains unused |
| Mng33a | Tc#n2#v3#t4; #n2#v2#h#t4; #n2#v#h2#t4; #n2#h3#t4 | Implemented; 4 | TcH₃(PPh₃)₄: JEA21, pp. 3096–3097, explicitly assigns paired d⁴, supported by resolved NMR; X-ray diffraction locates all three classical hydrides. All explicit/mixed/implicit-H forms; the H−-donor state reuses Mng32. The mixed PPh₃/PMe₃ nonclassical hydride is not included. |
| Mng33b | Re#c+2#n2#u#t6 | Implemented; 1 | ReCp*₂: HW 1929; M24, accepted manuscript pp. 3–4, 9–10, explicitly assigns the low-spin d⁵ ²E₂ ground state, consistent with earlier magnetic and spectroscopic measurements. Two six-electron Cp*− donors; two pairs and one unpaired electron on Re. Counts 17 electrons; no MC state. |
| Mng33c | Tc#n2#u#v2#t4; #c+2#n2#u#t6 | Implemented; 2 | TcX₂(dppe)₂, X = Cl/Br: FH70, pp. 454, 459–460, explicitly assigns octahedral low-spin d⁵ (²T₂), supported by magnetic moments and optical spectra. Two ordinary Tc–X bonds and four P donors, or six anionic/neutral donors; 17 electrons. FN59 leaves TcI₂(diars)₂'s anomalous moment unresolved; that compound is not used as occupation evidence. |
| Mng33d | Re#n2#u#v2#t4 | Implemented; 1 | ReX₂(diars)₂, X = Cl/Br/I: FN59, p. 1039, explicitly assigns octahedral low-spin d⁵ with one unpaired electron. Two ordinary Re–X bonds and four As donors; the six-donor form reuses Mng33b. Both allocations count 17 electrons. |
| Mng33e | M#c+#n#u2#v2#t4; #c+3#n#u2#t6 | Implemented; 4 | [TcX₂(dppe)₂]+ (X = Cl/Br) and [ReX₂(diars)₂]+ (X = Cl/Br/I): FH70, pp. 454, 458–459, assigns low-spin d⁴ (³T₁); FN59, p. 1039, explicitly relates the Re moments to two unpaired electrons with strong spin–orbit coupling. One pair plus two unpaired electrons. Two ordinary M–X bonds and four P/As donors, or six donors; 16 electrons. Distinct from the paired seven-coordinate d⁴ states. |
| Mng33f | Re#n#u2#v3#t3 | Implemented; 1 | mer-ReX₃(PR₃)₃, X = Cl/Br: NMR68, p. 420, gives the ³T₁g(t₂g⁴) parent term and its spin–orbit-coupled J = 0 ground level; CLM69, pp. 1674, 1677–1678, explicitly applies the t₂g⁴/second-order-paramagnetism description to ReCl₃(PR₃)₃. GL71, pp. 2230–2231, confirms the chloride/bromide magnetic behaviour. The parent occupation gives one pair and two unpaired electrons; J = 0 does not supply a #n2 alternative. Three ordinary Re–X bonds and three P donors; the six-donor form reuses Mng33e. Both count 16 electrons. |
| Mng33g | Re#c+2#n2#v#t6; #c+2#n2#h#t6 | Implemented; 2 | Cp₂ReH: GPW58, pp. 3916–3918, 3921, establishes the terminal hydride by NMR/IR and diamagnetism over 77–510 K. LH76, p. 1731, diagram 4, explicitly fills two nonbonding orbitals, supporting paired d⁴. Two six-electron Cp− donors and explicit/implicit H; the H−-donor state reuses Mng24. All allocations count 18 electrons. |
| Mng34 | Polymeric MF₅ and dimeric ReBr₅ | Fixture candidates; 0 new states | HW 1921–1922 establishes the connectivity and excludes Re–Re bonding in the pentahalide dimers. HP60, pp. 1101–1102, reports ReF₅ μeff = 1.41 μB at 295 K and antiferromagnetic interactions without assigning local spin. S = 1 fixture descriptions may use analogy with Re₂Cl₁₀ (Mng34a): all-halide donation reuses M#c+5#u2#t6 (Mng17); Re also has the mixed ordinary/donor state in Mng34a. No new Tc#c+#u2#v4#t2 or MC state is selected |
| Mng34a | Re#c+#u2#v4#t2 | Implemented; 1 | Re₂Cl₁₀: HW 1921–1922 and V24, preprint pp. 3–4, 7, 10, explicitly establish S = 1 per d² Re and no Re–Re bond. Four ordinary terminal Re–Cl bonds and two bridging Cl− donors per Re; the all-Cl−-donor form reuses Mng17. No MC state. |
| Mng35a | Re#c+4#n#v#t4#m0; Re#c+4#v3#t4#m0 | Implemented; 2 | Re₂H₈(PR₃)₄: four bridging and four terminal hydrides, two P donors per Re. HW 1919–1920, Fig. 28.7; FRW89, p. 3208, explicitly allows ordinary single- and triple-Re–Re descriptions. The monodentate-phosphine family includes PMe₃, PEt₃ and PMe₂Ph. Separate Re–H–Re 3c2e units and terminal H− donors; charge/contribution allocation: PW. No whole-dimer MC unit |
| Mng35b | Re#c+2#n#v3#t3#m0; Re#c+2#n#v3#t4#m0 | Implemented; 2 | Re₂H₄(PMe₂Ph)₄[P(OCH₂)₃CEt]₂: three bridging hydrides, one terminal hydride, three P donors per Re and an ordinary Re≡Re bond. GHC82, p. 2319, locates all hydrides and explicitly assigns the triple bond. This is a specific mixed-phosphine/phosphite example of HW 1920, Fig. 28.7c's core, not evidence for every [ReH₂(PR₃)₃]₂ formulation. Separate 3c2e bridges and terminal H− donation; charge/contribution allocation: PW |
| Mng36 | Eight-coordinate Tc/Re(IV) diars complexes; Mn₁₂/Mn₈₄ aggregates, higher Tc/Re clusters, solid nitridometallates and Mg₃ReH₇ | Outside initial default selection | HW 1904, 1909–1910, 1918–1920, 1923–1927. No per-scaffold MC entries or isolated finite molecules inferred from extended structures. An eight-donor d³ state would count 19 electrons; it is not used to propose a Tc/Re capacity increase |
| Mng37 | Intact-H₂ adducts; mixed-valent dimers with formal half-integer M–M orders; reduced phthalocyanines/noninnocent dithiolates | Deferred representation/occupation; not counted | HW 1902, 1920, 1922, 1928. Keep ligand-based redox and H₂ distinct from ordinary hydride/metal-ion assignments |
| Mng38 | TcF₇, TcCl₆, M(CN)₈³−/M(CN)₉²−, ReR₅/ReR₇; [MnBr(pmdta)]+ as a five-coordinate example | No addition from these citations | HW 1921, 1924, 1930 marks the first species unknown or uncertain. On p. 1902, the printed Mn formula has four donor atoms but is labelled trigonal-bipyramidal; the page image confirms the inconsistency |

#### Manganese-group limits and coverage

For these accepted ordinary/dative states, the Counts electron count is
`7 − charge + v + h + 2t`; the nonbonding budget is `7 − charge − v − h = 2n + u`.
Every added state passes that arithmetic, is distinct under the registry defaults,
and is covered by inferred resolution in the 242 new fixtures.

The implemented **Mn capacity 18 → 21** covers Mng06: 19 for [Mn(edta)(H₂O)]²− and
21 for [Mn(NO₃)₄]²− in either listed allocation. The remaining 239 states fit 18.
Tc/Re/Bh retain 18; in particular, MH₉²− counts 18, despite its nine localized H
bonds. The four Re MC states in Mng35a/b also fit 18.

There is no default ValenceTable entry for Mn/Tc/Re/Bh. Counts consults the shifted
element's table and uses the actual element's capacity; a supplied literal H count
bypasses saturation-target enumeration. The fixtures supply H counts and cover all
selected explicit/mixed/implicit-H states without table changes. Existing free-ion
states are preserved independently of complex occupations.
Mng31b/Mng33/Mng34 are optional fixture examples using accepted types and explicitly
identified analogies; they do not expand the 243-state selection. Registry membership
and a passing fixture do not establish the actual occupation of those compounds.

Mng35a/b use the existing H#c-#m2 with one Re–H–Re 3c2e relation per bridge,
contributions [0,2,0]. Terminal hydrides use existing H#c-#d; P donors and the
ordinary Re–Re bond remain separate. Residual Re pairs follow this allocation,
not an independently assigned ionic occupation. Only these four Re#…#m0 states
are accepted; localized/implicit terminal-H variants are not added to this selection.
Three conformance fixtures cover the two octahydride graphs and the actual mixed-ligand
tetrahydride graph. AtomTyping and Counts agree under Strict and MostSaturated,
and the electron-count validator passes. All fit Re capacity 18.

Verification: all 243 additions are inferred in the new snapshots, with 1,018
resolved atom-state/invariant checks. MostSaturated AtomTyping and both Counts
policies agree throughout. Strict AtomTyping also agrees except for the unpinned
TcOCl₄− and ReOBr₄− inputs, which retain paired/triplet alternatives; MostSaturated
selects the documented paired state.
All 2,238 earlier inputs and snapshots are unchanged. The 2,480-case conformance
suite, 3,594 graph unit tests, 223 graph property tests and 307 chemistry tests pass;
formatting and Clippy for both changed crates pass.

#### Manganese-group primary-source checks

| Key | Source | Status and use |
| --- | --- | --- |
| S15 | Stolte et al., “Preliminary toxicity and ecotoxicity assessment of methyltrioxorhenium and its derivatives,” *Green Chem.* **17** (2015), 1136–1144. [DOI](https://doi.org/10.1039/C4GC01919A), [publisher HTML](https://pubs.rsc.org/en/content/articlehtml/2015/gc/c4gc01919a), [archived publisher text](../materials/valence/Stolte-2015-MTO-adducts-publisher-text.txt) | Publisher HTML reviewed; Fig. 2 (5) identifies the 4,4′-dimethyl-2,2′-bipyridine adduct used for Mng21's two-donor fixture. Journal and accepted-manuscript PDFs return 403; publisher text archived. No additional atom type |
| G85 | Girolami et al., *J. Chem. Soc., Dalton Trans.* (1985), 1339–1348. [DOI](https://doi.org/10.1039/DT9850001339), [PDF](../materials/valence/Girolami-1985-divalent-dmpe-complexes.pdf) | Existing archived PDF reviewed for Mn: high-spin halides, pp. 1340–1341; low-spin MnMe₂(dmpe)₂, p. 1343. Mng05/Mng07 |
| WA03MnRe | Wang and Andrews, “Matrix Infrared Spectra and Density Functional Theory Calculations of Manganese and Rhenium Hydrides,” *J. Phys. Chem. A* **107** (2003), 4081–4091. [DOI](https://doi.org/10.1021/jp034392i), [PDF](../materials/valence/Wang-Andrews-2003-Mn-Re-hydrides.pdf) | Reviewed and archived. MnH septet, MnH₂ sextet, MnH₂− quintet, pp. 4085–4088; ReH quintet, ReH₂ sextet and ReH₄ doublet, pp. 4086–4089. Quartet ReH₂ frequencies explicitly rejected against experiment; ReH₃ not identified. Mng29a–f |
| GW66 | Grove and Wilkinson, “Oxo-complexes of rhenium(V),” *J. Chem. Soc. A* (1966), 1224–1230. [DOI](https://doi.org/10.1039/J19660001224), [PDF](../materials/valence/Grove-Wilkinson-1966-rhenium-V-oxo-complexes.pdf) | Reviewed and archived. Explicit diamagnetism for ReOCl₃(PPh₃)(Me₂SO), p. 1224, and trans-ReO₂(en)₂+, p. 1227. Its inconclusive ReOCl₅²− magnetism is resolved by FL71. Mng31 |
| TJD80 | Trop, Jones and Davison, “Technetium cyanide chemistry: synthesis and characterization of technetium(III) and -(V) cyanide complexes,” *Inorg. Chem.* **19** (1980), 1993–1997. [DOI](https://doi.org/10.1021/ic50209a032), [PDF](../materials/valence/Trop-Jones-Davison-1980-technetium-cyanides.pdf) | Reviewed and archived. Tc(CN)₇⁴− isolated and diamagnetic, with explicit singlet assignment, pp. 1994–1995. Tc(V) oxocyanide connectivity established, pp. 1995–1997; occupation not explicitly assigned. Mng31b/Mng32 |
| RJ23 | Roca Jungfer and Abram, “[ReH₃(PPh₃)₄] – A Key Compound in the Rhenium Hydride Chemistry,” *Chem. Eur. J.* **29** (2023), e202203317. [DOI](https://doi.org/10.1002/chem.202203317), [PDF](../materials/valence/Roca-Jungfer-Abram-2023-rhenium-phosphine-hydrides.pdf) | Reviewed and archived. Classical ReH₃(PPh₃)₄ supported by IR/DFT and comparison with Tc, p. 3; its H positions were not located crystallographically. ReH₅(PPh₃)₃ polymorphs, pp. 9–10, include located/refined H positions. Pairing in Mng24 is supported by HW's closed-shell assignment, not an electronic-spin measurement in RJ23. Ref. 31 identifies JEA21 for Tc |
| FL71 | Fergusson and Love, “The nature of the [ReOCl₅]²− anion,” *Aust. J. Chem.* **24** (1971), 2689–2693. [DOI](https://doi.org/10.1071/CH9712689), [PDF](../materials/valence/Fergusson-Love-1971-oxopentachlororhenate.pdf) | Reviewed and archived. Pure Cs/Rb salts are spin-paired; ReCl₆²− contamination explains earlier large moments, pp. 2689–2690. Table 2, p. 2691, also supports the paired oxo/nitrido phosphine complexes. Mng31/Mng31a |
| J13 | Johnstone et al., “X-ray Crystallographic and First-Principles Theoretical Studies of K₂[TcOCl₅] and UV/Vis Investigation of the [TcOCl₅]²− and [TcOCl₄]− Ions,” *Eur. J. Inorg. Chem.* (2013), 1097–1104. [DOI](https://doi.org/10.1002/ejic.201201346), [PDF](../materials/valence/Johnstone-2013-technetium-oxochlorides.pdf) | Reviewed and archived. Five-/six-coordinate structures and electronic levels, pp. 1098–1101; no new magnetic measurement. T79 supplies the explicit paired assignment. Supplement download failed; not required for Mng31c |
| JEA21 | Jungfer, Elsholz and Abram, “Technetium Hydrides Revisited: Syntheses, Structures, and Reactions of [TcH₃(PPh₃)₄] and [TcH(CO)₃(PPh₃)₂],” *Organometallics* **40** (2021), 3095–3112. [DOI](https://doi.org/10.1021/acs.organomet.1c00274), [PDF](../materials/valence/Jungfer-Elsholz-Abram-2021-technetium-hydrides.pdf) | Reviewed and archived. Explicit paired d⁴ interpretation supported by NMR, p. 3096; located classical H ligands, p. 3097. Mixed-phosphine product has a separate nonclassical description, p. 3098. Mng33a; monohydride reuses Mng14 |
| T79 | Trop, “Synthesis, characterization, and reactivity of technetium and rhenium complexes in intermediate oxidation states,” MIT PhD thesis (1979). [Repository](https://hdl.handle.net/1721.1/37642), [PDF](../materials/valence/Trop-1979-technetium-rhenium-thesis.pdf) | Retrieved, reviewed and archived. TcOCl₄− susceptibility, p. 48, and explicit singlet assignment, pp. 87–88; earlier TcOCl₅²− diamagnetism summarized on p. 74. Mng31c |
| B98 | Baldas et al., “Spectroelectrochemical and computational studies of tetrachloro and tetrabromo oxo- and nitrido-technetium(V) and their TcVI counterparts,” *J. Chem. Soc., Dalton Trans.* (1998), 2303–2314. [DOI](https://doi.org/10.1039/a802559b), [PDF](../materials/valence/Baldas-1998-technetium-oxo-nitrido-halides.pdf) | Reviewed and archived. Explicit paired Tc(V) and d¹ Tc(VI), pp. 2303–2307. Confirms Mng31d; no additional state |
| V24 | Vorobyova et al., “Anisotropy of exchange interactions in honeycomb ladder compound ReCl₅,” *Mater. Chem. Phys.* **332** (2025), 130215. [DOI](https://doi.org/10.1016/j.matchemphys.2024.130215), [2024 preprint](https://arxiv.org/abs/2411.02965), [archived preprint](../materials/valence/Vorobyova-2024-rhenium-pentachloride-preprint.pdf) | Preprint reviewed and archived. Explicit S = 1, L = 1, J = 2 assignment from magnetic entropy, pp. 4, 10; no Re–Re bond, p. 7. Mng34a uses the spin occupation, not the J multiplicity |
| JH75 | Jasinski and Holt, “Electronic structure and optical properties of manganese(V) oxytrichloride, MnOCl₃,” *Inorg. Chem.* **14** (1975), 1267–1273. [DOI](https://doi.org/10.1021/ic50148a013), [PDF](../materials/valence/Jasinski-Holt-1975-manganese-oxytrichloride.pdf) | Reviewed and archived. Explicit triplet ground state, p. 1269; optical-spectrum assignment, pp. 1270–1271. Mng30a adds the chloride-donor state; Mng12 already supplies the ordinary graph |
| GKB78 | Griffith, Kiernan and Brégeault, “Studies on transition-metal cyano-complexes. Part 2. Unsubstituted cyanorhenates, [Re(CN)x]n−, and cyanorhenates with thio-, seleno-, and nitrosyl ligands,” *J. Chem. Soc., Dalton Trans.* (1978), 1411–1417. [DOI](https://doi.org/10.1039/DT9780001411), [PDF](../materials/valence/Griffith-Kiernan-Bregeault-1978-rhenium-cyanides.pdf) | Reviewed and archived. Confirms low-spin d⁴ Re(CN)₇⁴−, p. 1412, and Re(CN)₆⁵−, pp. 1414–1415. Oxocyanide formulations discussed on p. 1415 without an explicit occupation assignment. Reassigned higher-valent hexacyanide and mixed-product preparations supply no additional state |
| C98 | Chang, Connick, Low, Day and Gray, “Electronic Structures of Nitridomanganese(V) Complexes,” *Inorg. Chem.* **37** (1998), 3107–3110. [DOI](https://doi.org/10.1021/ic970598k), [PDF](../materials/valence/Chang-1998-manganese-nitrides.pdf), [SI](../materials/valence/Chang-1998-manganese-nitrides-SI.pdf) | Main text and SI reviewed and archived. Explicit paired ¹A′ ground state, p. 3108; square-pyramidal terminal-nitrido structures, p. 3109. Mng30b |
| W84 | Winkler, “Spectroscopy and Photochemistry of Metal-Oxo Complexes,” Caltech PhD thesis (1984), chapter 4. [DOI](https://doi.org/10.7907/DE5E-2Y69), [PDF](../materials/valence/Winkler-1984-metal-oxo-thesis.pdf) | Retrieved and archived. ReO₂(CN)₄³− preparation, p. 83; explicit diamagnetic ¹A₁g[(b₂g)²] assignment and spectral interpretation, pp. 148, 152. Mng31e. Provides the primary evidence despite the blocked Winkler–Gray 1985 article download |
| F82 | Franklin et al., “Chemical applications of technetium-99 NMR spectroscopy: preparation of novel technetium(VII) species and their characterization by multinuclear NMR spectroscopy,” *J. Am. Chem. Soc.* **104** (1982), 5303–5306. [DOI](https://doi.org/10.1021/ja00384a007), [PDF](../materials/valence/Franklin-1982-technetium-NMR.pdf) | Reviewed and archived. Explicit diamagnetic d² TcO₂(CN)₄³− assignment and ⁹⁹Tc NMR, pp. 5303, 5305: Mng31f. TcO₃+ and its AsF₆− salt, pp. 5304, 5306: Mng13a. Also confirms TcH₉²− in Mng22. Its tentative Tc₂O₅F₄/TcO₂F₃ structures supply no additional state from this paper alone |
| M24 | May et al., “Electronic Structure and Photophysics of Low Spin d⁵ Metallocenes,” *Inorg. Chem.* **63** (2024), 1858–1866. [DOI](https://doi.org/10.1021/acs.inorgchem.3c03451), [PDF](../materials/valence/May-2024-low-spin-metallocenes.pdf) | Accepted manuscript retrieved from NSF PAR, reviewed and archived. Explicit ²E₂ assignments for MnCp*₂ and ReCp*₂, manuscript pp. 3–4, 9–10, 19; calculations agree with prior magnetic/spectroscopic assignments. Mng33b; Mn reuses Mng07 |
| FN59 | Fergusson and Nyholm, “New Oxidation States of Technetium,” *Nature* **183** (1959), 1039–1040. [DOI](https://doi.org/10.1038/1831039a0), [PDF](../materials/valence/Fergusson-Nyholm-1959-technetium-oxidation-states.pdf) | Reviewed and archived. ReX₂(diars)₂ explicitly low-spin d⁵; [ReX₂(diars)₂]+ has two unpaired d electrons with strong spin–orbit coupling, p. 1039. Mng33d/e. The anomalous TcI₂(diars)₂ moment is explicitly unresolved and supplies no occupation |
| FH70 | Fergusson and Hickford, “Studies in the chemistry of technetium. Complexes of technetium-(II), -(III), and -(IV) with 1,2-bis(diphenylphosphino)ethane,” *Aust. J. Chem.* **23** (1970), 453–461. [DOI](https://doi.org/10.1071/CH9700453), [PDF](../materials/valence/Fergusson-Hickford-1970-technetium-diphosphine.pdf) | Reviewed and archived. Measured moments, p. 454; explicit low-spin d⁴/³T₁ for [TcX₂(dppe)₂]+, pp. 458–459, and low-spin d⁵/²T₂ for TcX₂(dppe)₂, pp. 459–460 (X = Cl/Br). Mng33c/e. TcBr₄(PPh₃)₂ quartet assignment confirms Mng16, pp. 459–460 |
| RS65 | Randall and Shaw, “Assignment of the proton magnetic resonance spectra of some tris-(diethylphenylphosphine) complexes of Re(III), Re(IV), Ir(III) and Rh(III),” *Molecular Physics* **10** (1965), 41–48. [DOI](https://doi.org/10.1080/00268976600100061), [PDF](../materials/valence/Randall-Shaw-1965-rhenium-phosphine-NMR.pdf) | Reviewed and archived. Shifted sharp NMR resonances and meridional structure for ReCl₃(PMe₂Ph)₃/ReCl₃(PEt₂Ph)₃, pp. 44–47. Its contact/pseudocontact interpretation is revised by NMR68/CLM69; occupation support for Mng33f comes from those later papers |
| CR62 | Chatt and Rowe, “Complex compounds of tertiary phosphines and a tertiary arsine with rhenium(V), rhenium(III), and rhenium(II),” *J. Chem. Soc.* (1962), 4019–4033. [DOI](https://doi.org/10.1039/JR9620004019), [PDF](../materials/valence/Chatt-Rowe-1962-rhenium-phosphines.pdf) | Reviewed and archived. ReCl₃(PEt₂Ph)₃ preparation, pp. 4031–4032; the blanket Re(III) diamagnetism statement on p. 4026 is corrected by C68/GL71 and supplies no occupation. Paired Re(V) oxophosphines, pp. 4026–4028, confirm the existing Mng31 proposal |
| C68 | Chatt, Leigh, Mingos and Paske, “Complexes of osmium, ruthenium, rhenium, and iridium halides with some tertiary monophosphines and monoarsines,” *J. Chem. Soc. A* (1968), 2636–2641. [DOI](https://doi.org/10.1039/J19680002636), [PDF](../materials/valence/Chatt-1968-phosphine-arsine-halides.pdf) | Reviewed and archived. Explicitly corrects earlier diamagnetism of mer-ReCl₃(PR₃)₃, p. 2638; Table 4, p. 2640, reports moments and compositions. The retained bromide diamagnetism claim is corrected by GL71. Mng33f |
| GL71 | Gunz and Leigh, “Magnetic susceptibilities of some rhenium(III) and osmium(IV) halide complexes; preparation of some new d⁴ complexes,” *J. Chem. Soc. A* (1971), 2229–2233. [DOI](https://doi.org/10.1039/J19710002229), [PDF](../materials/valence/Gunz-Leigh-1971-rhenium-osmium-magnetism.pdf) | Reviewed and archived. Tables 1–2, pp. 2230–2231, establish solid/solution and largely temperature-independent paramagnetism, including mer-ReBr₃(PEt₂Ph)₃; p. 2231 rejects its earlier diamagnetism report. NMR68/CLM69 supply the electronic interpretation used in Mng33f. Table 1's separate diamagnetic fac-ReCl₃(PMe₂Ph)₃ observation supplies no additional occupation |
| CLM69 | Chatt, Leigh and Mingos, “Configurations of some complexes of rhenium, ruthenium, osmium, rhodium, iridium, and platinum halides with mono(tertiary phosphines) and mono(tertiary arsines),” *J. Chem. Soc. A* (1969), 1674–1680. [DOI](https://doi.org/10.1039/J19690001674), [PDF](../materials/valence/Chatt-1969-phosphine-arsine-configurations.pdf) | Reviewed and archived. Explicit t₂g⁴ ReCl₃(PR₃)₃ assignment, p. 1674; meridional phosphine environments, pp. 1677–1678. The J = 0/J = 1 second-order Zeeman interpretation on p. 1678 revises RS65's contact/pseudocontact account. Mng33f |
| NMR68 | Chatt, Leigh, Mingos, Randall and Shaw, “Nuclear magnetic resonance in complexes exhibiting second-order paramagnetism,” *Chem. Commun. (London)* (1968), 419–420. [DOI](https://doi.org/10.1039/C19680000419), [PDF](../materials/valence/Chatt-1968-second-order-paramagnetism.pdf) | Reviewed and archived. Re(III)/Os(IV) t₂g⁴ complexes, p. 419; p. 420 explicitly derives the J = 0 ground level from the ³T₁g parent by spin–orbit coupling, applying the same reasoning to Re. Supports the #n#u2 parent occupation in Mng33f |
| CL66 | Cotton and Lippard, “Chemical and Structural Studies of Rhenium(V) Oxyhalide Complexes. II. M[ReX₄O] and M[ReX₄OL] Complexes from KReO₄,” *Inorg. Chem.* **5** (1966), 9–16. [DOI](https://doi.org/10.1021/ic50035a003), [PDF](../materials/valence/Cotton-Lippard-1966-rhenium-oxyhalides.pdf) | Reviewed and archived. ReOX₄− salts (X = Cl/Br/I), pp. 9–11; p. 14 reports the bromide's susceptibility and explicitly assigns paired d² to the square-pyramidal family. Mng31g |
| GPW58 | Green, Pratt and Wilkinson, “Biscyclopentadienylrhenium hydride,” *J. Chem. Soc.* (1958), 3916–3922. [DOI](https://doi.org/10.1039/JR9580003916), [PDF](../materials/valence/Green-Pratt-Wilkinson-1958-rhenocene-hydride.pdf) | Reviewed and archived. Direct Re–H bond established by NMR/IR, pp. 3916–3918; p. 3921 gives diamagnetism over 77–510 K. LH76 supplies the explicit paired occupation used in Mng33g |
| FS67 | Fischer and Schmidt, “Dicyclopentadienyltechnetium Hydride,” *Angew. Chem. Int. Ed. Engl.* **6** (1967), 93. [DOI](https://doi.org/10.1002/anie.196700931), [PDF](../materials/valence/Fischer-Schmidt-1967-technetocene-hydride.pdf) | Reviewed and archived. Molecular weight, mass spectrum and 10:1 Cp/H NMR integration establish monomeric Cp₂TcH; IR establishes Tc–H. Bent geometry proposed by analogy with Re. No magnetic measurement or explicit electronic occupation; Mng33 permits a fixture by analogy without new states |
| HP60 | Hargreaves and Peacock, “New Fluorides and Oxyfluorides of Rhenium,” *J. Chem. Soc.* (1960), 1099–1103. [DOI](https://doi.org/10.1039/JR9600001099), [PDF](../materials/valence/Hargreaves-Peacock-1960-rhenium-fluorides.pdf) | Reviewed and archived. ReF₅ preparation and magnetic susceptibility, pp. 1100–1102; antiferromagnetic interactions are inferred, without a local-spin assignment. Mng34 permits fixtures by analogy without new states. Its d¹ ReOF₄ discussion concerns a different compound and is not evidence for the pentafluoride occupation |
| LH76 | Lauher and Hoffmann, “Structure and Chemistry of Bis(cyclopentadienyl)-MLₙ Complexes,” *J. Am. Chem. Soc.* **98** (1976), 1729–1742. [DOI](https://doi.org/10.1021/ja00423a017), [author copy](https://roaldhoffmann.com/sites/default/files/fromd6/109s_0.pdf), [PDF](../materials/valence/Lauher-Hoffmann-1976-bent-metallocenes.pdf) | Retrieved, reviewed for the Cp₂MH assignment and archived. Diagram 4 and accompanying text, p. 1731, explicitly show two filled nonbonding orbitals for d⁴ Cp₂ReH; p. 1738 identifies the 18-electron configuration. Mng33g. Tc is not explicitly treated |
| FRW89 | Fanwick, Root and Walton, “Reactions of the Dirhenium(II) Complexes Re₂Cl₄(PR₃)₄ with Lithium Aluminum Hydride. A Convenient Synthetic Route to the Dirhenium Octahydride Complexes Re₂H₈(PR₃)₄,” *Inorg. Chem.* **28** (1989), 3203–3209. [DOI](https://doi.org/10.1021/ic00315a024), [PDF](../materials/valence/Fanwick-Root-Walton-1989-dirhenium-octahydrides.pdf) | Reviewed and archived. Monodentate-phosphine family and four-bridge structure; p. 3208 explicitly permits single- and triple-Re–Re interpretations. Mng35a. The distinct dppm two-bridge structure is not included in this proposal |
| GHC82 | Green, Huffman and Caulton, “Could One Have Five? Synthesis and Structure of Re₂H₄(PMe₂Ph)₄[P(OCH₂)₃CEt]₂ and Its Conjugate Acid,” *J. Am. Chem. Soc.* **104** (1982), 2319–2320. [DOI](https://doi.org/10.1021/ja00372a040), [PDF](../materials/valence/Green-Huffman-Caulton-1982-dirhenium-hydrides.pdf) | Reviewed and archived. Three bridging and one terminal hydride located/refined; explicit Re≡Re assignment, p. 2319. Mng35b uses the actual neutral mixed-ligand compound; the protonated cation is not included |
| DAH79 | Dedieu, Albright and Hoffmann, “Some Hydrido-Bridged Transition-Metal Dimers and Their Unsupported Analogues. Speculations on Pentuple Bonding and Pentuple Bridging,” *J. Am. Chem. Soc.* **101** (1979), 3141–3151. [DOI](https://doi.org/10.1021/ja00506a001), [PDF](../materials/valence/Dedieu-Albright-Hoffmann-1979-hydrido-bridged-dimers.pdf) | Retrieved, reviewed for Re₂H₈L₄ and archived. Figure 4, pp. 3144–3146, treats four bridging hydrides as H₄⁴− in a delocalized orbital analysis. This supports the hydride allocation; separate 3c2e graph units remain PW, not a unique partition established by the paper |

Intermediates, page images and the deduplication/count audit are in
`scratch/hw-manganese-review/`. Retained papers belong in `materials/valence/`;
the table and its citations do not depend on retaining scratch.

### Iron group

HW II, chapter XXIX, pp. 1934–1988 (chemistry through p. 1986), with chapter-XXXII
cross-references. Extraction and primary-source review started 2026-11-04.
**Implemented: 210 additions — Fe 81, Ru 57, Os 68, S 3 and O 1.**
These comprise 209 ordinary-bond/donor states and one hydride-bridge MC state;
Fe capacity 21 and Os capacity 20 are implemented. The two Fe(1-norbornyl)₄ states
are excluded because they are supported here only by this specialized example.
All additions and the 17 retained states have conformance coverage; citations are updated.
Additional research is concluded for this selection. End-on FeOO remains unselected;
uncertain spin orderings that reuse accepted types do not add prospective states.

#### Iron-group work table

E = Fe, Ru, Os; M = Ru, Os. Element prefixes apply to every suffix in the same
cell; a slash-separated prefix applies to each named element. New-state counts
exclude existing entries and reuse across rows, and include every displayed H-count
split. States use the concrete, zero-valence defaults. Oxo/nitrido bonds are ordinary
multiple bonds; anionic-ligand donation follows the preceding chapters' convention.

| ID | States / extracted family | Status; new states | Examples and evidence |
| --- | --- | --- | --- |
| Fg01 | E/Hs#n2#u4; Fe#c+2#n#u4; Fe#c+2#n2#u2; Fe#c+2#n3; Fe#c+3#u5; Fe#c+3#n#u3; Fe#c+3#n2#u; Ru#c+2#n3; Ru#c+3#n2#u; Fe#c+#u5#v2; Fe#u5#v3; Fe#c-#u5#v4; Fe#c-2#u5#v5; Fe#c-3#u5#v6 | Retained; 17 existing; citations updated | Atomic occupations: HW 1538; Fe spin occupations: 1962–1965; Ru low-spin II/III: 1984. Formula-charge comments corrected to FeX₂+, FeX₃, FeX₅²− and FeX₆³−. FeCl₃, FeCl₄−, Fe(N₃)₅²− and FeF₆³− support #v3–#v6 (1949–1950, 1953). B00Fe, pp. 38–39 and Table 4-3 (p. 83), explicitly calculates sextet bare Fe(OH)₂+ and FeCl₂+, supporting existing #c+#u5#v2. B00Fe supplies the bare hydroxide fixture citation; HW 1957 describes the hydrated aqueous species. |
| Fg02 | Ru/Os#c+2#n#u4; #c+3#u5 | Implemented; 4 | Free-ion d⁶/d⁵ ground occupations, distinct from the low-spin complexes below. ASD lists ⁵D₄ and ⁶S₅/₂ for both elements. Existing Fe high-spin ions already supply these occupations. |
| Fg03 | Fe#c-2#n#u4#v4; #c+2#n#u4#t4; #c+2#n#u4#t6; #n#u4#v2#t4; #c-2#n#u4#v4#t2; #c-4#n#u4#v6 | Implemented; 6 | High-spin FeCl₄²−, Fe(H₂O)₆²+/Fe(NH₃)₆²+, FeCl₂(H₂O)₄, FeCl₄(H₂O)₂²− and Fe(OH)₆⁴−. HW 1947, 1949–1950, 1957, 1964. Ordinary anionic-ligand bonds or ligand donation. |
| Fg04 | Fe#c+3#u5#t3; #c+3#u5#t4; #c+3#u5#t5; #c+3#u5#t6; #u5#v3#t2; #c+2#u5#v#t5; #c+#u5#v2#t4; #u5#v3#t3; #c-#u5#v4#t2; #c-2#u5#v5#t; #c+#u5#v2#t2 | Implemented; 11 | High-spin Fe[N(SiMe₃)₂]₃, FeCl₄−, FeCl(acac)₂, Fe(H₂O)₆³+, hydrated chlorides/hydroxides, FeF₅(H₂O)²− and Fe₂Cl₆. HW 1947, 1949–1950, 1957–1958, 1962–1964. Fe₂Cl₆ uses two ordinary terminal Cl bonds and two bridging Cl− donors per Fe; Fe₂Cl₉³− reuses the three-terminal/three-donor state. Existing unadducted #v3/#v4/#v6 states are reused where applicable. |
| Fg05 | Fe/Ru/Os#c+2#n3#t6; #c-4#n3#v6; #n3#v2#t4; #c+3#n2#u#t6; #c-3#n2#u#v6 | Implemented; 15 | Low-spin II/III hexacyanides, neutral N-donor complexes, and paired II dihalide/phosphine or diorganyl complexes. HW 1949–1952, 1962–1965, 1974, 1976–1977, 1980, 1984–1986. E(NH₃)₆²+ is evidence for Ru/Os, not low-spin Fe. Ferro-/rutheno-/osmocene reuse the six-accepted-pair II state with two six-electron Cp− donors (2189–2192). |
| Fg06 | Fe#n2#u2#v2#t2; #c+2#n2#u2#t4; #n2#u2#v2#t4; #c+2#n2#u2#t6 | Implemented; 4 | Intermediate-spin Fe(II): FeTPP and Fe(phen)₂(ox). HW explicitly assigns two unpaired electrons, 1947, 1962–1963, 1965. Two ordinary anionic N/O bonds with remaining donors, or all donors. |
| Fg07 | Fe#n#u3#v3#t2; #c+3#n#u3#t5 | Implemented; 2 | Intermediate-spin FeCl(S₂CNR₂)₂, explicitly three unpaired electrons; HW 1962–1963. Two ordinary anionic S bonds and one ordinary Cl bond, with two neutral S donors, or all five donors. |
| Fg08 | Fe#n#u4#v2#t3; #c+2#n#u4#t5; #c+#n2#u#v2#t4 | Implemented; 3 | Five-coordinate high-spin deoxyheme and six-coordinate low-spin ferric cytochrome; HW 1965–1967. Ferrous cytochrome reuses Fg05. Two ordinary porphyrin N bonds with remaining donors, or all donors; no new haptic/MC state. |
| Fg09 | Fe#n3#u#v#t4; #c+#n3#u#t5; #c+2#n#u2#v2#t4; #c+4#n#u2#t6 | Implemented; 4 | FeI(CO)₂(diars), explicitly low-spin Fe(I), and [FeCl₂(diars)₂]²+, explicitly triplet Fe(IV). HW 1962, 1969. Ordinary halide bonds or halide donation. |
| Fg10 | Fe#n2#v4; #c+4#n2#t4 | Not selected; 0 | Tetrahedral Fe(1-norbornyl)₄ is explicitly diamagnetic/paired d⁴ (HW 1970). Excluded from the default registry: only this specialized organometallic example supports these states in the extraction. |
| Fg11 | Fe#u2#v6; #c-2#u2#v8 | Implemented; 2 | FeO₄²−: HW 1969 explicitly gives two unpaired electrons. Two Fe=O plus two Fe–O− bonds, or four Fe=O with charge on Fe. |
| Fg12 | Fe#c+3#u5#t7; #c-#u5#v4#t3; #c+3#u5#t8; #c-#u5#v4#t4 | Implemented; 4 | Fe(edta)(H₂O)− and Fe(NO₃)₄−: HW 1947 and high-spin discussion 1963. K71, p. 554, confirms four nearly symmetric bidentate nitrates and compares the nonbonding occupation to isoelectronic Mn(II). Seven/eight donors, or four ordinary anionic O bonds with the other donors retained. Counts gives 19/21 electrons; Fe capacity 18 → 21 implemented. |
| Fg13a | S#c-2#n2#d2; #c-2#n#d3 | Implemented; 2 | μ₂-/μ₃-sulfide donors in the common Fe₂S₂/Fe₄S₄ families; HW 1961, 1967. Fe(II/III) local high-spin four-donor states reuse Fg03/Fg04. S²− donation is the adopted ionic graph convention; no cluster-specific MC state. |
| Fg13b | O#c-2#n2#d2 | Implemented; 1 | μ₂-oxide donor alternative for [Cl₃Fe–O–FeCl₃]²− and the hydrated Fe(III) oxide bridge, HW 1950, 1957, 1964. Existing O#n2#v2 supplies the ordinary bridge. Local high-spin Fe(III) is retained despite antiferromagnetic coupling. |
| Fg13c | S#c-#n2#v#d | Implemented; 1 | General terminal thiolate/dithiocarbamate donor, needed by the selected Fe–S cluster and all-donor FeCl(S₂CNR₂)₂ graphs. HW 1961–1962; anionic-donor allocation: PW. |
| Fg14 | Fe/Ru/Os#n4#t5; #c-2#n5#t4 | Implemented; 6 | E(CO)₅ and E(CO)₄²−; HW 1947, 1970, 1984, 2111, 2135–2136. Explicit paired carbonyl occupations; Fe(CO)₄²− is Collman's reagent. |
| Fg15 | Fe/Ru/Os#n3#v#h#t4; #n3#h2#t4 | Implemented; 6 | Classical EH₂(PR₃)₄ and EH₂(CO)₄, all explicit/mixed/implicit-H forms; HW 1948, 1974, 2136–2137. Explicit octahedral 18-electron descriptions. H− donors reuse E#c+2#n3#t6; the fully explicit state also reuses Fg05. |
| Fg16 | Fe/Ru/Os#c-4#n3#v5#h; #c-4#n3#v4#h2; #c-4#n3#v3#h3; #c-4#n3#v2#h4; #c-4#n3#v#h5; #c-4#n3#h6 | Implemented; 18 | EH₆⁴− in Mg₂FeH₆, Mg₂RuH₆ and Mg₂OsH₆: isolated octahedral 18-electron/low-spin units; HW 1948, 1973–1974. All seven H-count splits; explicit state and H−-donor state reuse Fg05. |
| Fg17 | Fe/Ru/Os#c-#n4#v#t4; #c-#n4#h#t4 | Implemented; 6 | EH(CO)₄−, HW 2136–2137. Explicit/implicit terminal H; H−-donor allocation reuses E#n4#t5. |
| Fg18 | Fe#n2#v4#t3 | Implemented; 1 | Fe₂(CO)₉ and the two bridge-bearing sites of Fe₃(CO)₁₂; HW 2111–2114, Fig. 32.28. Ordinary Fe–Fe bonds and C(=O)(–Fe)₂ bridges; terminal CO remains a donor. The remaining Fe₃ site and Ru₃/Os₃(CO)₁₂ reuse #n3#v2#t4. No new MC state. |
| Fg19 | Os#n2#v4#t3; #n2#v3#h#t3; #n2#v2#h2#t3; #n2#v#h3#t3; #n2#h4#t3; #c+4#n2#t7; #n#v6#t2; #n#v5#h#t2; #n#v4#h2#t2; #n#v3#h3#t2; #n#v2#h4#t2; #n#v#h5#t2; #n#h6#t2; #c+6#n#t8 | Implemented; 14 | Classical OsH₄(PR₃)₃ and OsH₆(PR₃)₂, all H-count splits and H− donors; HW 1974, 2071. Their 18-electron occupations are paired d⁴/d². The nonclassical Ru homologues do not support corresponding terminal-hydride states. |
| Fg20 | Ru/Os#n3#v2#t3; #c+2#n3#t5 | Implemented; 4 | Five-coordinate MCl₂(PPh₃)₃, HW 1976–1977, low-spin II discussion 1984. RuHCl(PPh₃)₃ reuses the explicit state; its implicit H is listed separately. |
| Fg20a | Ru#n3#v#h#t3 | Implemented; 1 | RuHCl(PPh₃)₃, HW 1976–1977: one implicit hydride and one ordinary chloride; five-donor ionic alternative reuses Fg20. |
| Fg21 | Ru/Os#n2#u#v3#t3 | Implemented; 2 | MCl₃(PR₃)₃, HW 1976, 1984. Explicit low-spin d⁵; six-donor alternative reuses Fg05. |
| Fg21a | Ru#c+2#n2#u#v#t5; #c+#n2#u#v2#t4; #c-#n2#u#v4#t2; #c-2#n2#u#v5#t | Implemented; 4 | [RuClₙ(H₂O)₆₋ₙ]³⁻ⁿ, n = 1–5, HW 1976; n = 0, 3, 6 reuse earlier rows. The ordinary Cl bonds preserve each isolated/aquated charge state; all-donor graphs reuse Fg05. |
| Fg22 | Ru/Os#c-2#n#u2#v6; #c+4#n#u2#t6 | Implemented; 4 | MX₆²−: HW 1983–1984 explicitly gives low-spin d⁴ with two unpaired electrons. The spin–orbit-coupled J = 0 limit is not a paired-d⁴ alternative. Ordinary halides or donors. |
| Fg22a | Os#n#u2#v4#t2 | Implemented; 1 | OsCl₄(PR₃)₂/OsCl₄(AsR₃)₂: NMR68, pp. 419–420, explicitly derives the J = 0 level from a ³T₁g parent. Four ordinary Cl bonds and two P/As donors; the all-donor state reuses Fg22. |
| Fg23 | Ru/Os#c-#u3#v6; #c+5#u3#t6; #u2#v6; #c+6#u2#t6 | Implemented; 8 | MF₆− (V, three unpaired) and MF₆ (VI, two unpaired); HW 1974–1976, 1982–1983. Octahedral species with explicit occupations, not an inference about gas-phase MF₅. |
| Fg24 | Ru#c-2#u2#v8 | Implemented; 1 | RuO₄²−, tetrahedral d² with two unpaired electrons; HW 1978–1979, 1982. All-double-bond allocation; two ordinary Ru=O plus two Ru–O− bonds reuse #u2#v6. |
| Fg25 | Ru/Os#c-2#n#v8; #c+2#n#v4#t4; #c-#n#v7; #c+3#n#v3#t4; #c+3#n#v3#t5 | Implemented; 10 | Paired VI oxo/nitrido families MO₂X₄²−, MNX₄− and MNX₅²−; HW 1978–1979, 1982–1983 explicitly gives M=O/M≡N and the paired occupation. Ordinary or donating halides. The MO π-MO discussion adds no MC alternative. OsO₂(en)₂²+ reuses the dioxo/donor state. |
| Fg26 | Ru/Os#u#v7; #c-#u#v8; #v8 | Implemented; 6 | MO₄− (VII, one unpaired) and MO₄ (VIII, paired d⁰); HW 1978–1979, 1982. Ordinary oxo charge allocations. The ordinary OsOF₅ graph also uses #u#v7; unsupported OsF₇ is not an example. |
| Fg27 | Os#c+5#u#v2#t5; #c+2#v6#t2; #c+4#v4#t4; #c-#v9; #v8#t | Implemented; 5 | Donor-F alternatives for OsOF₅ (VII), OsO₃F₂ and OsO₂F₄ (VIII); OsO₃N− with Os≡N and the OsO₄·D adduct. HW 1974–1975, 1978, 1982; SS06 confirms monomeric OsOF₅. Ordinary oxyfluoride graphs reuse Fg26. |
| Fg27a | OsF₇; monomeric OsO₂F₃ | Not selected; 0 | Do not select Os#c+7#u#t7 for OsF₇ or Os#c+3#u#v4#t3 for monomeric OsO₂F₃. SS06, pp. 7931–7934, finds only OsF₆ from attempted OsF₇ synthesis and resolves solid OsO₂F₃ as mixed-valence [OsO₃F₂·OsOF₄]₂/∞. The possible gas-phase monomer remains tentative. No replacement oligomer-specific state is proposed. |
| Fg27b | OsOF₆ | Not selected; 0 | Do not select Os#c+6#v2#t6 for OsOF₆. C93Os, pp. 11279–11284, corrects the reported product to cis-OsO₂F₄ using composition, diffraction and spectroscopy; already covered by Fg26/Fg27. |
| Fg28 | Os#c-2#v10; #v8#t2 | Implemented; 2 | cis-OsO₄F₂²− and cis-OsO₄(OH)₂²−; HW 1978, 1982 explicitly gives d⁰ and the six-coordinate structures. Four ordinary Os=O plus two ordinary anionic-ligand bonds, or two anionic donors. Both count 20 electrons; Os capacity 18 → 20 implemented. |
| Fg29 | Reuses preceding rows | Implemented; 0 | Os₂X₈²−: explicit diamagnetic Os≡Os (HW 1976, 1984–1985, Fig. 29.18d). Ordinary Os#c-#n#v7 and halide-donor Os#c+3#n#v3#t4 reuse Fg25. No MC bond or axial-donor extension. |
| Fg30 | Ru#n#v6; #c+3#n#v3#t3; #n3#v2#t2; #c+#n3#v#t3 | Implemented; 4 | Ru₂R₆ and Ru₂Cl₂(PMe₃)₄: HW 1984–1985, Table 29.5, gives paired cluster configurations and Ru≡Ru/Ru–Ru respectively. Ordinary metal–metal bonds; ordinary R/Cl bonds or donors. R = CH₂CMe₃/CH₂SiMe₃. |
| Fg31 | Ru#n2#v4#t2; #c+2#n2#v2#t4 | Implemented; 2 | Common Grubbs-type RuCl₂(=CHR)(NHC)(PR₃) precatalysts. S19Ru, pp. 2–3, explicitly uses Ru=C and a singlet ground state. Ordinary Ru=C and Cl bonds with two neutral donors, or Cl− donation. No scaffold-specific type or MC encoding. |
| Fg32 | Fe#c-2#n3#v4#t2 | Implemented; 1 | trans-[Fe(CN)₄(CO)₂]²−: R01, p. 6934, establishes the formula, structure and low-spin Fe(II) assignment; C02 confirms both cis/trans structures. Four ordinary CN bonds and two CO donors; all-donor graph reuses Fg05. HW 1953's CN₂ formula is erroneous. K03's Fe(0) [Fe(CN)₂(CO)₃]²− is a different compound, not a correction to the Fe(II) example. |
| Fg32a | Fe#c-3#n3#v5#t; #c-2#n2#u#v5#t | Implemented; 2 | Fe(II/III) pentacyanide complexes Fe(CN)₅D, HW 1953, with low-spin occupations from 1962–1964. Five ordinary CN bonds and one donor; all-donor forms reuse Fg05. Nitroprusside fits the Fe(II) state with HW's explicit NO+ allocation; the distinct brown-ring nitrosyl is covered in Fg39b. |
| Fg33 | Ru#n#u4#v2; #n2#v4 | Implemented; 2 | Molecular RuO, RuO₂ and RuO₃, with ordinary Ru=O bonds. M16Ru, pp. 74–78, 80–84, 94–95, supports quintet RuO and singlet RuO₂/RuO₃ using measured structures/spectra and correlated calculations. The RuO₂ triplet is explicitly higher. Solid rutile RuO₂ does not supply the molecular assignment. |
| Fg34 | Fe#n#u4#v2 | Implemented; 1 | Molecular FeO, ⁵Δ ground state: W96, pp. 5298–5299, Fig. 5; C96, pp. 5266–5267, matrix and gas-phase vibrational evidence. Ordinary Fe=O; the source's Fe-derived nonbonding orbitals supply one pair and four unpaired electrons. |
| Fg35 | Fe#n3#v2#t3#m0 | Implemented; 1 | HFe₂(CO)₈−: CB78, pp. 2314–2317, confirms one symmetric hydride bridge, two μ-CO and three terminal CO per Fe, and explicitly discusses three-center M–H–M bonding. Accepted H− contributions [0,2,0] use the adopted convention (PW); μ-CO retains the ordinary representation of Fg18. This graph counts no additional Fe–Fe bond; that partition is not determined by the crystallographic distance. HW 2137–2138, Fig. 32.44g. One MC unit, not a whole-cluster description. |
| Fg36 | Bare RuH/OsH, RuH₂/OsH₂, MH₃/MH₄ and H₂ adducts | Not selected / deferred; 0 | WA09 leaves the earlier isolated RuH assignment provisional (p. 560); bare OsH and RuH₂/OsH₂ were not observed. No terminal MH₃/MH₄ assignment is selected; the observed MH₃ compositions are (H₂)MH. H₂ complexes retain the separate donor-model deferral; no MC substitute is proposed. The observed Fe species are below. |
| Fg36a | Fe#n2#u3#v; #n2#u3#h; #c+#n2#u3#t | Implemented; 3 | FeH, observed quartet ⁴Δ ground state: WA09, pp. 551, 553–554, Table 4. Ordinary explicit/implicit H and H− donation. |
| Fg36b | Fe#n#u4#v#h; #n#u4#h2; #c+2#n#u4#t2 | Implemented; 3 | FeH₂, observed quintet with gas-phase ⁵Δ assignment: WA09, pp. 551, 553–554, Table 4. Ordinary H-count splits and H− donation; the fully explicit state reuses Fg34. |
| Fg36c | Fe#c-#n2#u3#v2; #c-#n2#u3#v#h; #c-#n2#u3#h2; #c+#n2#u3#t2 | Implemented; 4 | FeH₂−: WA09, Table 4, explicitly calculates ⁴Δg for this previously photodetachment-observed anion (p. 559, ref. 44); not observed in WA09's matrix experiment. Ordinary H-count splits and H− donation. |
| Fg37 | FeMes₂ | Not selected; 0 | K94 establishes [Fe₂Mes₄] with two terminal and two bridging mesityls. Its magnetic analysis requires coupled centers and supports direct Fe–Fe interaction; it does not establish monomeric Fe#n2#u2#v2. No additional aryl-bridge or cluster-specific state is selected. Ru/Os tetraorganyls are settled below. |
| Fg37a | Ru#c+4#n2#t4; Os#n2#v4; #c+4#n2#t4 | Implemented; 3 | MR₄ (R = cyclohexyl, o-tolyl; also RuMes₄): HW 1985; PRE12, pp. 5711–5712, explicitly identifies the characterized tetrahedral compounds as diamagnetic d⁴ and describes two nonbonding pairs. Four ordinary M–C bonds or anionic C donors. Ru#n2#v4 reuses Fg33; no MC state. |
| Fg38 | End-on FeOO | Not selected; no new state | C96, pp. 5269–5270, establishes the observed triplet superoxo isomer but does not fix the local Fe/O spin allocation. No Fe-local occupation is inferred from the total triplet. Higher iron oxides are settled in Fg38i–Fg38k; no calculated-only tetrahedral/diperoxo FeO₄ states or weak O₂ adducts are added. |
| Fg38a | Fe#n#u2#v4; #c-#n2#u#v4 | Implemented; 2 | OFeO and OFeO− with two ordinary Fe=O bonds. HT12, Fig. 1 and pp. 3092–3095, assigns neutral ³Δg and anionic ²Δg ground states using RASPT2/RCCSD(T) and the photoelectron spectrum. The Fe-derived nonbonding occupations supply #n#u2 and #n2#u. Higher quartet/quintet/sextet alternatives are not selected; the separate FeOO/Fe(O₂) isomers are not covered by this assignment. |
| Fg38b | Os#c+4#n#v2#t4 | Implemented; 1 | OsOX₄ (X = Cl, F): HR57, pp. 25, 29, explicitly pairs the two remaining d electrons of OsOCl₄; H88, pp. 64–65, independently assigns dxy² to monomeric OsOF₄. Ordinary Os=O plus four X− donors; the all-ordinary Os#n#v6 state reuses Fg39c. |
| Fg38c | Reuses preceding rows | Implemented; 0 | Cyclic Fe(O₂): C96, p. 5269, assigns the observed isomer as a quintet peroxide from isotopic vibrational data. An O–O single bond with two ordinary Fe–O bonds uses Fg34; two O− donors use Fg36b and existing O#c-#n2#v#d. No new state or MC bond. |
| Fg38d | Os#n#u2#v4 | Implemented; 1 | Molecular OsO₂ with two ordinary Os=O bonds. Z00, pp. 3461–3462, identifies the matrix-isolated molecule and assigns a ³B₁ ground state with both BP86 and B3LYP; the singlet lies 12.9/14.5 kcal/mol higher. The two nonbonding unpaired electrons and one pair give this state. |
| Fg38e | Ru/Os#c-#n2#u#v4 | Implemented; 2 | RuO₂− and OsO₂− with two ordinary M=O bonds and charge on M. Z00, Table 5 and pp. 3463–3464, assigns the observed anions as ²A₁ and ²B₁ respectively; both functionals agree. Closed-shell O atoms leave two pairs and one unpaired electron on M. |
| Fg38f | Reuses preceding rows | Implemented; 0 | Molecular OsO₃: Z00, Table 6 and p. 3462, identifies the matrix-isolated species and assigns a planar ¹A₁′ ground state. Three ordinary Os=O bonds reuse Fg39c; no new state. |
| Fg38g | Ru#c+2#n#v4#t2 | Implemented; 1 | Observed singlet (O₂)MO₂ peroxide isomers: Z00, Table 7 and pp. 3462–3463, explicitly gives M(VI), an O–O bond and two terminal oxo groups. Ordinary bonds reuse M#n#v6. Making the two peroxo O− atoms donors adds only the Ru state; Os reuses Fg39c and both use existing O#c-#n2#v#d. No MC bond. |
| Fg38h | Ru#u2#v6#t2 | Implemented; 1 | RuO₃(OH)₂²−: A06, pp. 1921–1922, explicitly selects the triplet for the characterized five-coordinate ion; inclusion of spin–orbit coupling leaves its predominantly triplet character. Three ordinary Ru=O bonds plus two OH− donors add #u2#v6#t2; the all-ordinary graph reuses Fg24. Both count 18 electrons. OsO₂(OH)₄²− remains the paired six-coordinate case in Fg25. |
| Fg38i | Fe#n#v6; #c+2#n#v4#t2 | Implemented; 2 | Observed singlet (O₂)FeO₂: C96, pp. 5270–5271, gives the structure and singlet assignment; W96, pp. 5299–5300, independently supports two terminal Fe=O bonds and an O–O unit. Two Fe=O plus two ordinary Fe–O bonds, or two peroxide O− donors. No tetrahedral FeO₄ state is inferred. |
| Fg38j | Fe#c+2#u4#v2#t2 | Implemented; 1 | Observed quintet (O₂)FeO: C96, p. 5270, gives the C₂v isomer and multiplicity. One Fe=O plus a cyclic peroxide; the ordinary graph reuses Fg39f. Two peroxide O− donors add the second state. No MC bond. |
| Fg38k | Reuses preceding rows | Implemented; 0 | Three-terminal-oxo FeO₃: W96, p. 5300, supports the constitution. X08, pp. 15884–15886, finds singlet and triplet structures separated by only 0.020 eV and explicitly leaves their ordering unresolved. Both supplied-spin graphs reuse independently accepted Fg38i/Fg11 states; no prospective occupation is added. |
| Fg39 | Ru/Os μ-O/μ-N and carboxylate dimers | Research concluded | Explicit descriptions and state reuse are recorded in Fg39h–Fg39k. No additional MC state is required. |
| Fg39a | Fe#c+2#u4#v2; #c+2#u4#v2#t5 | Implemented; 2 | Ferryl FeO²+ and [FeO(H₂O)₅]²+: HW 1969 identifies the hydrated intermediate; L08, pp. 89–93, explicitly gives S = 2 for both bare and hydrated FeO²+. Ordinary Fe=O, with five water donors in the hydrated graph. No extrapolation to Fe(V) or different ligand fields. |
| Fg39b | Fe#c+#n2#u3#t6 | Implemented; 1 | Brown-ring [Fe(H₂O)₅(NO)]²+: C23, Fig. 1a–c and pp. 16864–16865, explicitly compares Fe(I)/NO+, high-spin Fe(II)/NO and high-spin Fe(III)/triplet NO−. All three are representable. N→Fe with N≡O+ uses existing N#v3#d and O#c+#n#v3; the NO and NO− graphs carry #c-#u and #c-2#u2#s3 on N–O (PW). Only quartet Fe(I) adds a type; 19 electrons fit accepted Fe capacity 21. Local spins are supplied; exchange coupling is not inferred by resolution. |
| Fg39c | Os#n#v6; #c+2#n#v4#t2 | Implemented; 2 | Common Os(VI) dioxoglycolate ester in HW 1978, Fig. 29.15: two Os=O plus two Os–O bonds, or two alkoxide donors. H19 main text and SI p. S9 explicitly select low spin; Fig. SI5 compares singlet/triplet energies. Bis-pyridine/TMEDA adducts reuse #n#v6#t2 and #c+2#n#v4#t4. No MC description or further calculated-intermediate states are selected. |
| Fg39d | Fe#u3#v5 | Implemented; 1 | FeO₄³−: one Fe=O and three Fe–O− bonds. S16 confirms K₃FeO₄; N18 SI Table SI11 explicitly uses a quartet (S = 3/2, ⟨S²⟩ = 3.810), consistent with W56's 3.6–3.7 μB. Closed-shell O atoms leave three unpaired electrons on Fe in this graph. Count: 13 electrons. |
| Fg39e | Fe#c-3#u3#v8 | Implemented; 1 | Alternative all-double-bond FeO₄³− graph for Fg39d, placing charge on Fe. Same quartet occupation; count 19 fits the accepted Fe capacity 21. |
| Fg39f | Fe#u4#v4 | Implemented; 1 | FeO₄⁴− in Na₄FeO₄: four Fe–O− bonds. S16 confirms the salt; N18 identifies tetrahedral units and SI Table SI11 explicitly uses a quintet (S = 2, ⟨S²⟩ = 6.045). Closed-shell O atoms leave four unpaired electrons on Fe. Count: 12 electrons. This does not extend to HW's nonmolecular Sr/Ba ferrates(IV). |
| Fg39g | Fe#c-4#u4#v8 | Implemented; 1 | Alternative all-double-bond FeO₄⁴− graph for Fg39f, placing charge on Fe. Same quintet occupation; count 20 fits the accepted Fe capacity 21. |
| Fg39h | Ru/Os#c+4#n2#t6 | Implemented; 2 | [M₂OCl₁₀]⁴−: HW 1980, 1983 footnote and 1984 explicitly gives π-induced pairing; D20Ru, pp. 889–890/Fig. 5, supplies a closed-shell Ru–O–Ru MO description. Adopt the all-donor graph with two paired metal nonbonding pairs and bridging O²−#n2#d2 (PW); each metal accepts five chloride pairs and one oxide pair. This does not replace the distinct #n#u2 occupations of MX₆²−. No MC state or new O state. |
| Fg39i | Ru/Os#c-3#n2#v7; #c+2#n2#v2#t5 | Implemented; 4 | [M₂NCl₁₀]⁵− and [M₂NCl₈(H₂O)₂]³−: HW 1983–1984 explicitly gives diamagnetic M=N=M; DZ20, p. 5, independently identifies Ru=N double bonds in the hydrated dimer. Ordinary M=N+=M uses existing N#c+#v4, with the paired nonbonding remainder on M (PW). Five ordinary Cl bonds give the first state; five terminal donors give the second. No new N or MC state. |
| Fg39j | Reuses preceding rows | Implemented; 0 | Ru₂(OAc)₄(THF)₂ and mixed-valent Ru₂(O₂CR)₄X₂−: HW 1984–1985 gives BO 2/2.5 and two/three unpaired electrons; G26 confirms the latter occupation for R = n-Pr, X = Br. Adopt M–M strings 2#u2#s3 and 3#c-#u3#s4 respectively, with four carboxylate donors and an axial donor per metal (PW). These reuse Fg39i/Fg25. Bond charge avoids assigning inequivalent II/III metal sites; no MC bond. |
| Fg39k | Reuses preceding rows | Implemented; 0 | Os₂(OAc)₄Cl₂: GP13, archived PDF pp. 5, 15–16, gives σ²π⁴δ²δ*¹π*¹ and two unpaired electrons, supported by magnetic studies and calculations. Four carboxylate donors and one Cl− donor per Os, with Os–Os string 3#u2#s3 (PW), reuse Fg25. Reduced Os₂(OAc)₄Cl₂− also needs no new atom type if charge/spin are bond-carried; no unique reduced-ion multiplicity is selected and no additional fixture is required. |
| Fg40 | Larger carbonyl/hydride cages, heterometallic assemblies, polymeric Mg₂RuH₄/Mg₃RuH₃, borohydrides, η²-H₂ complexes | Unselected / deferred | HW 1948, 1973–1974, 1984–1986, 2113–2138. No additional per-cage MC state, and no isolated molecular ion inferred from an extended solid's formula. The established η²-H₂ donor limitation remains separate work. Small ordinary-bond carbonyl clusters are covered above. |

#### Limits and coverage

Capacities are **Fe 21 and Os 20**; Ru/Hs retain 18. Fe covers the seven/eight-donor
high-spin complexes (19/21 electrons), the ferrate(V/IV) charge alternatives (19/20)
and the Fe(I) brown-ring description (19). Os covers the six-coordinate Os(VIII)
oxo/fluoride states in Fg28 (20). ValenceTable targets are unchanged; supplied H
counts remain necessary where Counts has no saturation target.

**223 new conformance fixtures** cover all 210 additions and complete coverage of
all 17 retained states. They comprise 137 donor–acceptor, 36 inorganic-ion,
25 hydride, 17 small-molecule and 8 atomic-ion inputs. Every new fixture resolves
under AtomTyping and Counts with both policies; exact expected atom types,
agreement of all four outputs and valence invariants were checked independently.
The three deoxyheme/cytochrome coordination models explicitly identify their capped
ligands; the source-supported local occupations are supplied to resolution.

Two existing Fe halide inputs now supply their known high-spin occupations;
the bare Fe(OH)₂+ source comment is corrected. The spin-unspecified Ru²⁺/Ru³⁺
inputs now record Strict ambiguity and the existing low-spin MostSaturated result.
Four existing snapshots change. The full **2,703-case** suite, 3,594 graph unit tests,
223 graph property tests and 307 chemistry tests pass; formatting and Clippy pass.

#### Iron-group sources

| Key | Source | Access and contribution |
| --- | --- | --- |
| ASD | NIST Atomic Spectra Database, ground-state tables, retrieved 2026-11-04: [Ru](../materials/valence/NIST-ASD-2026-ru-ground-states.html), [Os](../materials/valence/NIST-ASD-2026-os-ground-states.html); [database](https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html) | Reviewed and archived. Free Ru²⁺/Os²⁺: d⁶ ⁵D₄; Ru³⁺/Os³⁺: d⁵ ⁶S₅/₂. Fg02 does not substitute complex occupations for free ions. |
| K71 | King, Logan, Morris and Wallwork, “Eight-co-ordinate Iron(III) Complex: Crystal Structure of Tetraphenylarsonium Tetranitratoferrate(III),” *Chem. Commun.* (1971), 554. [DOI](https://doi.org/10.1039/C29710000554), [PDF](../materials/valence/King-1971-tetranitratoferrate.pdf) | User supplied; reviewed and archived. Discrete Fe(NO₃)₄−, four nearly symmetric bidentate nitrates; spherical nonbonding occupation compared with Mn(NO₃)₄²−. Fg12. |
| NMR68 | Chatt et al., “Nuclear magnetic resonance in complexes exhibiting second-order paramagnetism,” *Chem. Commun.* (1968), 419–420. [DOI](https://doi.org/10.1039/C19680000419), [PDF](../materials/valence/Chatt-1968-second-order-paramagnetism.pdf) | Already archived; OsCl₄L₂ passage rechecked. Explicit ³T₁g parent and spin–orbit-coupled J = 0 level; Fg22a. |
| S19Ru | Śliwa, Mitoraj, Sagan and Handzlik, “Formation of active species from ruthenium alkylidene catalysts—an insight from computational perspective,” *J. Mol. Model.* **25** (2019), 331. [DOI](https://doi.org/10.1007/s00894-019-4202-5), [PDF](../materials/valence/Sliwa-2019-ruthenium-alkylidene-catalysts.pdf) | Retrieved, reviewed and archived. Structures, pp. 1–2; explicit singlet assignment and comparison with triplet/quintet energies, p. 3. Fg31. |
| M16Ru | Miradji, “Quantum Modelling of Ruthenium Chemistry in the field of Nuclear Power Plant Safety,” Lille 1 PhD thesis (2016). [Institutional copy](https://www.irsn.fr/sites/default/files/documents/larecherche/formation_recherche/theses/theses-soutenues/psn-res/2016-These-Miradji.pdf), [PDF](../materials/valence/Miradji-2016-ruthenium-thesis.pdf) | Relevant oxide sections reviewed and archived. Quintet RuO and singlet RuO₂/RuO₃; ordinary double-bond assignments, pp. 74–78, 80–84, 94–95. Fg33. Calculated-only oxyhydroxide families are not selected. |
| W96 | Wu, Desai and Wang, “Observation and Photoelectron Spectroscopic Study of Novel Mono- and Diiron Oxide Molecules: FeOᵧ− (y = 1–4) and Fe₂Oᵧ− (y = 1–5),” *J. Am. Chem. Soc.* **118** (1996), 5296–5301. [DOI](https://doi.org/10.1021/ja954153b), [author copy](https://bpb-us-w2.wpmucdn.com/sites.brown.edu/dist/0/196/files/2020/09/39.pdf), [PDF](../materials/valence/Wu-1996-iron-oxides.pdf) | Retrieved, relevant monoiron passages reviewed and archived. FeO occupation diagram, p. 5299; FeO₃ structural evidence. The paper explicitly leaves detailed higher-oxide electronic assignments to further calculations. Fg34/Fg38. |
| C96 | Chertihin et al., “Reactions of Laser-Ablated Iron Atoms with Oxygen Molecules in Condensing Argon. Infrared Spectra and Density Functional Calculations of Iron Oxide Product Molecules,” *J. Phys. Chem.* **100** (1996), 5261–5273. [DOI](https://doi.org/10.1021/jp953198w), [PDF](../materials/valence/Chertihin-1996-iron-oxides.pdf) | User supplied; relevant monoiron sections reviewed and archived. Observed isomers, isotopic assignments and explicit spin discussion, pp. 5266–5270; cyclic quintet peroxide, p. 5269. Does not establish observed tetrahedral FeO₄. Fg34/Fg38/Fg38c. |
| AA01 | Aullón and Alvarez, “Theoretical Study of Bonding, Structure, and Vibrational Spectra of the [Fe₂(CO)₈]²− Anion and Its Derivatives,” *Organometallics* **20** (2001), 818–826. [DOI](https://doi.org/10.1021/om000548e), [PDF](../materials/valence/Aullon-Alvarez-2001-iron-carbonylate.pdf) | User supplied; reviewed and archived. Parent carbonylate isomers and methylene derivatives; it does not analyze HFe₂(CO)₈−. The unbridged parent uses Fe#c-#n4#v#t4 from Fg17 with an ordinary Fe–Fe bond. |
| SS06 | Shorafa and Seppelt, “Osmium(VII) Fluorine Compounds,” *Inorg. Chem.* **45** (2006), 7929–7934. [DOI](https://doi.org/10.1021/ic0608290), [PDF](../materials/valence/Shorafa-Seppelt-2006-osmium-fluorides.pdf) | Reviewed and archived. OsOF₅ confirmed; OsF₇ synthesis not reproduced; solid OsO₂F₃ is mixed-valence, with a gas-phase monomer left tentative. Fg26–Fg27a. |
| CB92/C93Os | Christe and Bougon, “Osmium tetrafluoride dioxide, OsF₄O₂: a new osmium(VIII) oxide fluoride,” *Chem. Commun.* (1992), 1056a, [DOI](https://doi.org/10.1039/C3992001056A); Christe et al., “Osmium tetrafluoride dioxide, cis-OsO₂F₄,” *J. Am. Chem. Soc.* **115** (1993), 11279–11284, [DOI](https://doi.org/10.1021/ja00077a029), [PDF](../materials/valence/Christe-1993-osmium-tetrafluoride-dioxide.pdf) | C93Os relevant sections reviewed and full PDF archived; CB92 abstract reviewed. Corrects reported OsOF₆ to cis-OsO₂F₄ using composition, electron diffraction and spectroscopy; no OsOF₆-specific state selected. Fg27b. |
| R01 | Rauchfuss et al., “The Influence of Cyanide on the Carbonylation of Iron(II): Synthesis of Fe-SR-CN-CO Centers Related to the Hydrogenase Active Sites,” *J. Am. Chem. Soc.* **123** (2001), 6933–6934. [DOI](https://doi.org/10.1021/ja015948n), [PDF](../materials/valence/Rauchfuss-2001-cyanocarbonyl.pdf) | Reviewed and archived. Formula, structure and explicit low-spin Fe(II) assignment of trans-[Fe(CN)₄(CO)₂]²−, p. 6934. Fg32. |
| C02 | Contakes et al., “Preparative and Structural Studies on the Carbonyl Cyanides of Iron, Manganese, and Ruthenium: Fundamentals Relevant to the Hydrogenases,” *Inorg. Chem.* **41** (2002), 1670–1678. [DOI](https://doi.org/10.1021/ic011030l), [PDF](../materials/valence/Contakes-2002-cyanocarbonyls.pdf) | Retrieved, relevant iron sections reviewed and archived. Both cis/trans-[Fe(CN)₄(CO)₂]²− and [Fe(CN)₅(CO)]³−. Fg32/Fg32a. |
| K03 | Kayal and Rauchfuss, “Protonation Studies of the New Iron Carbonyl Cyanide trans-[Fe(CO)₃(CN)₂]²−: Implications with Respect to Hydrogenases,” *Inorg. Chem.* **42** (2003), 5046–5048. [DOI](https://doi.org/10.1021/ic034455b), [PDF](../materials/valence/Kayal-2003-cyanocarbonyl.pdf) | Retrieved, reviewed and archived. Establishes the distinct Fe(0) tricarbonyl dicyanide; distinguishes it from Fg32's Fe(II) compound. |
| HT12 | Hendrickx and Tran, “On the Electronic and Geometric Structures of FeO₂−/0 and the Assignment of the Anion Photoelectron Spectrum,” *J. Chem. Theory Comput.* **8** (2012), 3089–3096. [DOI](https://doi.org/10.1021/ct300516y), [PDF](../materials/valence/Hendrickx-Tran-2012-iron-dioxide.pdf) | Reviewed and archived. Neutral triplet and anionic doublet ground states; Fig. 1 and pp. 3092–3095 give the occupations used in Fg38a. |
| CB78 | Chin and Bau, “Crystal Structure of [(Ph₃P)₂N]+[HFe₂(CO)₈]−,” *Inorg. Chem.* **17** (1978), 2314–2317. [DOI](https://doi.org/10.1021/ic50186a061), [PDF](../materials/valence/Chin-Bau-1978-diiron-hydride.pdf) | Reviewed and archived. Symmetric hydride bridge and explicit three-center M–H–M discussion. The H− electron allocation and exclusion of an additional Fe–Fe bond in Fg35 are the accepted graph convention. |
| L08 | Louwerse, *Computational Modeling of Oxidation Catalysis: Studies concerning Fenton’s reaction*, Vrije Universiteit thesis (2008), chapter 5, pp. 85–104. [Thesis copy](https://www.scm.com/wp-content/uploads/Thesis_Manuel_Louwerse.pdf), [PDF](../materials/valence/Louwerse-2008-oxidation-catalysis-thesis.pdf); reprints Louwerse and Baerends (2007), [DOI](https://doi.org/10.1039/B613182D) | Relevant chapter reviewed and archived. Explicit S = 2 for bare and hydrated ferryl, pp. 89–93. Fg39a. |
| C23 | Coates, Banerjee and Odelius, “Simulations of the Aqueous ‘Brown-Ring’ Complex Reveal Fluctuations in Electronic Character,” *Inorg. Chem.* **62** (2023), 16854–16866. [DOI](https://doi.org/10.1021/acs.inorgchem.3c02320), [PDF](../materials/valence/Coates-2023-brown-ring.pdf) | Retrieved, relevant assignments/conclusions reviewed and archived. Fig. 1a–c supplies three explicit metal/NO charge and spin descriptions. Fg39b retains all three; the bond-carried NO charge/spin allocation is PW. |
| H19 | Hussein et al., “Mechanism of Os-Catalyzed Oxidative Cyclization of 1,5-Dienes,” *J. Org. Chem.* **84** (2019), 15173–15183. [DOI](https://doi.org/10.1021/acs.joc.9b02174), [PDF](../materials/valence/Hussein-2019-osmium-cyclization.pdf), [SI](../materials/valence/Hussein-2019-osmium-cyclization-SI.pdf), [SI source](https://doi.org/10.1021/acs.joc.9b02174.s001) | Relevant main-text and SI sections reviewed; both PDFs archived. Explicit low-spin Os(VI) glycolate and spin comparison (main p. 15175; SI p. S9) support Fg39c. |
| HR57 | Hepworth and Robinson, “Simple and complex oxyhalides of ruthenium and osmium,” *J. Inorg. Nucl. Chem.* **4** (1957), 24–29. [DOI](https://doi.org/10.1016/0022-1902(57)80104-3), [PDF](../materials/valence/Hepworth-Robinson-1957-oxyhalides.pdf) | Reviewed and archived. OsOCl₄ is diamagnetic (p. 25); p. 29 explicitly attributes this to pairing of its two d electrons. Fg38b. Does not establish an OsOF₄ occupation. |
| W56 | Wahl, Klemm and Wehrmeyer, “Über einige Oxokomplexe von Übergangselementen,” *Z. anorg. allg. Chem.* **285** (1956), 322–336. [DOI](https://doi.org/10.1002/zaac.19562850325), [PDF](../materials/valence/Wahl-1956-oxocomplexes.pdf) | Relevant ferrate(V) section reviewed and archived. The measured 3.6–3.7 μB supports three unpaired electrons, but sample identity was tentative, pp. 329–330. S16 confirms K₃FeO₄ and N18 independently uses the quartet occupation; Fg39d/Fg39e. |
| WA09 | Wang and Andrews, “Infrared Spectra and Theoretical Calculations for Fe, Ru, and Os Metal Hydrides and Dihydrogen Complexes,” *J. Phys. Chem. A* **113** (2009), 551–563. [DOI](https://doi.org/10.1021/jp806845h), [PDF](../materials/valence/Wang-Andrews-2009-iron-group-hydrides.pdf) | Reviewed and archived. Observed FeH/FeH₂ and previously observed FeH₂−, with explicit ground-state assignments (Table 4). Separates H₂ complexes from calculated-only bare hydrides; Fg36–Fg36c. |
| K94 | Klose et al., “Magnetic Properties Diagnostic for the Existence of Iron(II)–Iron(II) Bonds in Dinuclear Complexes Which Derive from Stepwise Insertion Reactions on Unsupported Iron–Aryl Bonds,” *J. Am. Chem. Soc.* **116** (1994), 9123–9135. [DOI](https://doi.org/10.1021/ja00099a030), [PDF](../materials/valence/Klose-1994-diiron-mesityl.pdf) | Relevant structure and magnetic sections reviewed and archived. Fe₂Mes₄, not a monomer; coupled-center analysis and evidence for direct Fe–Fe interaction, pp. 9131–9133. No new local occupation selected; Fg37. |
| PRE12 | Poblador-Bahamonde, Raynaud and Eisenstein, “Structures of d⁴ MH₃X: a Computational Study of the Influence of the Metal and the Ligands,” *Inorg. Chem.* **51** (2012), 5705–5715. [DOI](https://doi.org/10.1021/ic3001448), [institutional copy](https://archive-ouverte.unige.ch/unige:26053), [PDF](../materials/valence/Poblador-Bahamonde-2012-d4-tetracoordinate.pdf) | Relevant sections reviewed and archived. Explicitly identifies characterized Ru/Os tetraorganyls as diamagnetic d⁴ (p. 5711) and describes two nonbonding pairs; Fg37a. The original 1987 synthesis paper, DOI 10.1039/DT9870000557, returned 403; PRE12 supplies the explicit occupation used here. |
| S16 | Machalová Šišková et al., “High-valent iron (FeVI, FeV, and FeIV) species in water: characterization and oxidative transformation of estrogenic hormones,” *Phys. Chem. Chem. Phys.* **18** (2016), 18802–18810. [DOI](https://doi.org/10.1039/C6CP02216B), [accepted manuscript](../materials/valence/Siskova-2016-high-valent-iron.pdf) | Characterization sections reviewed and manuscript archived. Manuscript pp. 2–3 report 94% Fe(IV) in Na₄FeO₄ and 93% Fe(V) in K₃FeO₄, with >90% phase purity for each; residual Fe is Fe(III). Supplies identity, with spin assignments from N18. The unarchived SI is not needed for these main-text results. Fg39d–Fg39g. |
| N18 | Novák et al., “Transformations of ferrates(IV,V,VI) in liquids: Mössbauer spectroscopy of frozen solutions,” *Phys. Chem. Chem. Phys.* **20** (2018), 30247–30256. [DOI](https://doi.org/10.1039/C8CP05952G), [PDF](../materials/valence/Novak-2018-ferrate-transformations.pdf), [SI](../materials/valence/Novak-2018-ferrate-transformations-SI.pdf) | Relevant characterization and calculation sections reviewed; main/SI archived. Tetrahedral Na₄FeO₄ and K₃FeO₄; SI Table SI11 (PDF p. 14) explicitly uses quintet FeO₄⁴− and quartet FeO₄³− in the calculations. These are calculated multiplicities, distinct from the Mössbauer line-shape labels. Fg39d–Fg39g. |
| B00Fe | Bärsch, *Reactivity and Thermochemistry of Gaseous Iron Compounds of Different Valency*, TU Berlin thesis (2000). [Repository copy](https://webdoc.sub.gwdg.de/ebook/diss/2003/tu-berlin/diss/2000/baersch_susanne.pdf), [PDF](../materials/valence/Baersch-2000-gaseous-iron-thesis.pdf) | Relevant gas-phase calculations reviewed and thesis archived. Explicit sextet Fe(OH)₂+ (pp. 38–39) and FeCl₂+ (Table 4-3, p. 83) support the existing state in Fg01. |
| Z00 | Zhou et al., “Infrared Spectra and Density Functional Calculations of MO₂, MO₃, (O₂)MO₂, MO₄, MO₂− (M = Re, Ru, Os) and ReO₃−, ReO₄− in Solid Neon and Argon,” *J. Phys. Chem. A* **104** (2000), 3457–3465. [DOI](https://doi.org/10.1021/jp993026p), [PDF](../materials/valence/Zhou-2000-ruthenium-osmium-oxides.pdf) | Relevant Ru/Os sections reviewed and PDF archived. Observed OsO₂, RuO₂−/OsO₂−, OsO₃ and (O₂)MO₂, with explicit calculated occupations in Tables 5–7 and pp. 3461–3464; Fg38d–Fg38g. |
| H88 | Hope, Levason and Ogden, “Spectroscopic studies on matrix isolated osmium pentafluoride oxide, osmium tetrafluoride oxide, and osmium difluoride trioxide,” *J. Chem. Soc., Dalton Trans.* (1988), 61–65. [DOI](https://doi.org/10.1039/DT9880000061), [PDF](../materials/valence/Hope-1988-osmium-oxyfluorides.pdf) | Reviewed and archived. Monomeric OsOF₄ vibrational assignment, pp. 62–63; explicit dxy² ground occupation used to interpret its electronic spectrum, pp. 64–65. Fg38b. |
| A06 | Avdei et al., “Ab initio study of the comparative stability of ruthenium(VI) and osmium(VI) oxo hydroxo complexes with different coordination numbers,” *Russ. J. Inorg. Chem.* **51** (2006), 1920–1924. [DOI](https://doi.org/10.1134/S0036023606120138), [PDF](../materials/valence/Avdei-2006-ruthenium-osmium-oxohydroxo.pdf) | Reviewed and archived. Triplet RuO₃(OH)₂²− and singlet OsO₂(OH)₄²−, pp. 1921–1922; spin–orbit treatment retains their predominant spin character. Fg38h/Fg25; no states selected for the uncharacterized coordination alternatives. |
| X08 | Xue et al., “Experimental and Theoretical Study of the Reactions between Small Neutral Iron Oxide Clusters and Carbon Monoxide,” *J. Am. Chem. Soc.* **130** (2008), 15879–15888. [DOI](https://doi.org/10.1021/ja8023093), [author copy](https://bernsteinlab.colostate.edu/213.pdf), [PDF](../materials/valence/Xue-2008-neutral-iron-oxides.pdf) | Relevant monoiron sections reviewed and PDF archived. The small singlet/triplet FeO₃ separation and explicitly unresolved ordering, pp. 15884–15886, support retaining both supplied-spin descriptions through existing selections. Fg38k. |
| D20Ru | Dissanayaka Mudiyanselage et al., “Li₄Ru₂OCl₁₀·10H₂O: crystal structure, magnetic properties and bonding interactions in ruthenium-oxo complexes,” *Acta Cryst. B* **76** (2020), 884–891. [DOI](https://doi.org/10.1107/S2052520620010914), [repository copy](https://www.osti.gov/servlets/purl/1851622), [PDF](../materials/valence/Mudiyanselage-2020-ruthenium-oxo.pdf) | Reviewed and archived. Ru(IV) established by XPS; the body, Fig. 5 and conclusion explicitly describe paired occupied MOs, pp. 889–890. Fg39h uses this description, distinct from the isolated MX₆²− occupation. |
| DZ20 | Dzhabieva et al., “Biomimetic Water Oxidation Catalyzed by a Binuclear Ruthenium (IV) Nitrido-Chloride Complex with Lithium Counter-Cations,” *Biomimetics* **5** (2020), 3. [DOI](https://doi.org/10.3390/biomimetics5010003), [PDF](../materials/valence/Dzhabieva-2020-ruthenium-nitrido.pdf) | Relevant structural section reviewed and PDF archived. Discrete Ru₂NCl₈(H₂O)₂³−, a linear centrosymmetric bridge and explicit Ru=N double-bond assignment, pp. 4–5. Paired occupation comes from HW 1983–1984. Fg39i does not include the later tetranuclear catalyst. |
| G26 | Grenda et al., “Electron and spin density distributions and magnetic anisotropy in mixed-valent diruthenium(V) tetracarboxylates assessed via crystallographic and theoretical analyses,” *Dalton Trans.* **55** (2026), 7173–7183. [DOI](https://doi.org/10.1039/D6DT00342G), [repository copy](https://hal.science/hal-05613872), [PDF](../materials/valence/Grenda-2026-diruthenium.pdf) | Reviewed and archived from HAL. Polarized-neutron/X-ray measurements and calculations support equivalent Ru sites, the quartet and σ²π⁴δ²(δ*π*)³, pp. 7173–7174, 7178–7179. The atom/bond allocation in Fg39j is PW; the calculated effective bond order is not substituted for the nominal graph order. |
| GP13 | Gracia and Patmore, “Structural, spectroscopic and theoretical studies of diosmium(III,III) tetracarboxylates,” *Dalton Trans.* **42** (2013), 13118–13125. [DOI](https://doi.org/10.1039/C3DT51449H), [repository copy](https://eprints.hud.ac.uk/id/eprint/18802/), [PDF](../materials/valence/Gracia-Patmore-2013-osmium-carboxylates.pdf) | Relevant magnetic and electronic-structure sections reviewed; manuscript archived. PDF pp. 5 and 15–16 give the triplet σ²π⁴δ²δ*¹π*¹ description and distinguish it from the singlet alternative. Fg39k; no ligand-scaffold-specific state. |

Intermediates and the candidate/count audit are in `scratch/hw-iron-review/`.
Retained sources are archived in `materials/valence/`.

### Cobalt group

HW II, chapter XXX, pp. 1989–2022, with chapter-XXXII cross-references.
**Accepted selection implemented: 195 additions — Co 63, Rh 68, Ir 64.**
Capacities are Co 23 and Rh/Ir 19; no new MC or companion types, and no ValenceTable
change. All 13 earlier states remain: 12 have sources and conformance coverage;
Co38 remains source-unresolved. The other open and excluded cases below add no states.

#### Cobalt-group work table

Element prefixes apply to every suffix in the same cell; slash-separated prefixes
apply to each named element. Counts exclude existing states and reuse across rows,
and include every displayed H-count split. States use concrete, zero-valence defaults.
Anionic-ligand donation follows the adopted convention of the preceding chapters.
The capacity changes and all 195 additions are approved.

| ID | States / extracted family | Status; new states | Examples and evidence |
| --- | --- | --- | --- |
| Co01 | Co/Rh/Ir#n3#u3; Co#c+2#n2#u3; Co/Rh#c+2#n3#u; Co/Rh#c+3#n3; Co#c+3#n#u4; Co#c-2#n3#u#v4; Co#c-3#n3#u#v5; Rh#c+3#n3#t6 | Retained; 12 existing; citations updated | Atomic occupations: HW 1538 and ASD. Co(II/III): HW 1994, 2000–2002; Rh(II/III): 2017. Corrected Co #v4/#v5 example charges: Co(CN)₄²−/Co(CN)₅³−. Rh's #t6 hydride-donor state is supported by HW 2068, 2071 and by the six-donor III complexes here. The thirteenth existing state is unresolved in Co38. |
| Co02a | Co/Rh#c+#n3#u2 | Implemented; 2 | Free M+ ions: d⁸ ³F₄, ASD. These do not replace paired d⁸ complex occupations. |
| Co02b | Rh/Ir#c+2#n2#u3; #c+3#n#u4 | Implemented; 4 | Free M²+/M³+: d⁷ ⁴F₉/₂ and d⁶ ⁵D₄, ASD. Existing Co states already cover both; retain the existing low-spin Rh states. |
| Co02c | Ir#c+#n2#u4 | Implemented; 1 | Free Ir+ is 5d⁷6s, ⁵F₅ (ASD), not the d⁸ triplet of free Co+/Rh+. |
| Co03 | Co#c-2#n2#u3#v4; #c+2#n2#u3#t4; #n2#u3#v2#t2 | Implemented; 3 | High-spin CoCl₄²−, Co(H₂O)₄²+, CoCl₂D₂ (D = pyridine/phosphine). HW 1994, 2001–2002; K12 also identifies S = 3/2 CoX₂(PPh₃)₂. Ordinary anionic-ligand bonds or donation. |
| Co04 | Co#c+2#n2#u3#t6; #n2#u3#v2#t4; #c-4#n2#u3#v6 | Implemented; 3 | High-spin Co(H₂O)₆²+/Co(NH₃)₆²+, Co(acac)₂(H₂O)₂ and Co(OH)₆⁴−; HW 1994, 1996–1997, 2001–2002. Each graph counts 19 electrons, within Co capacity 23. |
| Co05 | Co#c+2#n3#u#t4; #c+2#n3#u#t5; #n3#u#v2#t2 | Implemented; 3 | Low-spin Co(CN)₄²−/Co(CN)₅³−, Co(CNR)₅²+, Co(salen) and Co(pc); HW 1994, 2001–2002. The ordinary cyanide states already exist. Salen/pc use two ordinary anionic N/O bonds plus two donors, or all donors. |
| Co06 | Co#c+2#n3#u#t6; #c-3#n3#u#v5#t | Implemented; 2 | Low-spin Co(terpy)₂²+ and Co(CN)₅(H₂O)³−, HW 1997, 2002. Cobaltocene also uses the first state with two six-electron Cp− donors (2189–2195). These graphs count 19 electrons. |
| Co07 | Co#c+3#n3#t6; #c+2#n3#v#t5; #c+#n3#v2#t4; #n3#v3#t3; #c-3#n3#v6 | Implemented; 5 | Paired Co(III): hexaaqua/ammine/cyanide, CoX(NH₃)₅²+, CoCl₂(NH₃)₄+, Co(acac)₃ and cobalamins; HW 1994, 1996–1997, 2000–2004. Cobalamin's corrin/axial anions can be ordinary bonds or donors. The Co(III) centers of the peroxide/superoxide examples reuse these states; ligand radical allocation does not require a metal MC state. |
| Co08 | Co#c-3#n#u4#v6; #c+3#n#u4#t6 | Implemented; 2 | High-spin CoF₆³−, HW 1994, 2000: four unpaired electrons. Ordinary F bonds or F− donation. CoF₃(H₂O)₃ is settled separately in Co30. |
| Co09 | Co#c-4#n4#v5; #n4#v#t4 | Implemented; 2 | Co₂(CN)₁₀⁶− and Co₂(CN)₈⁸− with ordinary Co–Co bonds; HW 1994–1995, 2002 explicitly gives pairing by dimerization. Co(CN)₅⁴− also uses Co#c-4#n4#v5. Co₂(CN)₁₀⁶− with ordinary CN bonds reuses Co07. |
| Co10 | Co#n4#h#t4; #c+#n4#t5; Rh/Ir#n4#v#t4; #n4#h#t4; #c+#n4#t5 | Implemented; 8 | MH(PR₃)₄ and HM(CO)₄: paired 18-electron monohydrides; HW 1992, 2009–2010, 2068–2069, 2136–2137. Explicit/implicit H and H− donation. M₂(CO)₈ with four terminal CO per M reuses the ordinary #v#t4 state for its M–M bond (2112). |
| Co11a | Co#c-4#n4#v4#h; #c-4#n4#v3#h2; #c-4#n4#v2#h3; #c-4#n4#v#h4; #c-4#n4#h5; Rh/Ir#c-4#n4#v5; #c-4#n4#v4#h; #c-4#n4#v3#h2; #c-4#n4#v2#h3; #c-4#n4#v#h4; #c-4#n4#h5 | Implemented; 17 | MH₅⁴− in Mg₂CoH₅ and Ca₂/Sr₂MH₅ (M = Rh, Ir): square-pyramidal, paired d⁸; HW 1992, 2009–2010, 2069. All six H-count splits; H− donors reuse Co10. |
| Co11b | Rh#c-3#n4#v4; #c-3#n4#v3#h; #c-3#n4#v2#h2; #c-3#n4#v#h3; #c-3#n4#h4; #c+#n4#t4 | Implemented; 6 | Square-planar RhH₄³− in Li₃RhH₄, HW 2009–2010, 2069; Rh(I) diamagnetism, 2017. All five H-count splits and H− donors. |
| Co11c | Rh#c-3#n3#v6; #c-3#n3#v5#h; #c-3#n3#v4#h2; #c-3#n3#v3#h3; #c-3#n3#v2#h4; #c-3#n3#v#h5; #c-3#n3#h6; Ir#c-3#n3#v6; #c-3#n3#v5#h; #c-3#n3#v4#h2; #c-3#n3#v3#h3; #c-3#n3#v2#h4; #c-3#n3#v#h5; #c-3#n3#h6; #c+3#n3#t6 | Implemented; 15 | Octahedral MH₆³− in Li₃IrH₆ and Na₃MH₆; HW 2009–2010, 2069. Paired III occupation, all seven H-count splits and H− donors. Rh's six-accepted-pair state already exists. |
| Co12 | Co#n3#u2#v#t3; #c+#n3#u2#t4 | Implemented; 2 | CoX(PR₃)₃, X = Cl/Br/I: tetrahedral S = 1; K12, pp. 4954–4956, 4962–4963, HFEPR and electronic-structure analysis. Ordinary halide or donation. Distinct from paired square-planar Rh/Ir(I). |
| Co13a | Co#c+#n4#t4; #n4#v#t3 | Implemented; 2 | Square-planar cob(I)alamin: K12, p. 4954, explicitly gives diamagnetic Co(I); HW 2003–2004 gives the corrin family. Four donors, or one ordinary anionic corrin N bond and three donors. |
| Co13b | Co#c+#n3#u#v#t3; #c+#n3#u#v#t4 | Implemented; 2 | Four-/five-coordinate cob(II)alamin, L09: EPR and electronic-structure analysis explicitly give low-spin d⁷, S = 1/2. One ordinary anionic corrin N bond and three/four donors. All-donor Co²+ states reuse Co05. |
| Co14 | Co#c+2#n2#u3#t7; #c-2#n2#u3#v4#t3; #c+2#n2#u3#t8; #c-2#n2#u3#v4#t4 | Implemented; 4 | Co(NO₃)₄²−: C16 establishes seven/eight O donors and S = 3/2. All O donors, or four ordinary anionic O bonds plus the other donors. Counts 21/23 electrons; covered by Co capacity 23. Eight coordination is also in HW 1991. |
| Co15 | Co#c-2#n2#u#v6; #c+4#n2#u#t6 | Implemented; 2 | CoF₆²−: low-spin Co(IV), one unpaired electron; HW 2004. Ordinary F bonds or F− donation. |
| Co16 | Rh/Ir#c+2#n3#v#t5; #c+#n3#v2#t4; #n3#v3#t3; #c-#n3#v4#t2; #c-2#n3#v5#t | Implemented; 10 | Paired M(III), [MClₖ(H₂O)₆₋ₖ]³⁻ᵏ (k = 0–6), corresponding ammines and hexacyanides; HW 2011–2014, 2017. The full substitution series is explicit. Ordinary anionic-ligand bonds or all donors. |
| Co17a | Rh/Ir#c-2#n2#u#v6; #c+4#n2#u#t6 | Implemented; 4 | MX₆²−, X = F/Cl/Br: low-spin d⁵, one unpaired electron; HW 2012, 2017. Ordinary bonds or donors. |
| Co17b | Ir#c+#n2#u#v3#t3; #n2#u#v4#t2; #c-#n2#u#v5#t | Implemented; 3 | IrCl₃(H₂O)₃+, IrCl₄(H₂O)₂ and IrCl₅(H₂O)−, HW 2012, 2017; same low-spin IV occupation. Their all-donor form reuses Co17a. |
| Co18 | Rh/Ir#c-#n#u2#v6; #c+5#n#u2#t6; #u3#v6; #c+6#u3#t6 | Implemented; 8 | MF₆−: octahedral d⁴ with two unpaired electrons (HW 2012). Molecular MF₆: octahedral d³ (2010–2011; t₂g³ occupation under HW's octahedral splitting, 1594). Ordinary F bonds or F− donors; no hypothetical IrF₇. |
| Co19a | Rh#n4#v#t3; Ir#n4#v#t3; #c+#n4#t4 | Implemented; 3 | Paired square-planar M(I): Wilkinson's RhCl(PPh₃)₃ and Vaska's IrCl(CO)(PPh₃)₂, HW 2012–2013, 2017. Ordinary Cl or chloride donation. The five-coordinate CO adduct of Vaska's complex reuses Co10. |
| Co19b | Rh#c+#n4#t3; #c-#n4#v2#t2 | Implemented; 2 | T-shaped Rh(PPh₃)₃+ (2013), square-planar Rh(CN)₄³− (2014), and the Monsanto catalyst RhI₂(CO)₂− (2021). Explicit paired Rh(I) description; all-donor alternatives reuse Co19a. |
| Co20 | Rh#n3#u#v2#t2; #c+2#n3#u#t4 | Implemented; 2 | RhX₂(PR₃)₂, square-planar low-spin Rh(II), one unpaired electron; HW 2012, 2017. No unsupported Ir analogue. |
| Co21a | Rh/Ir#n3#v2#h#t3; #n3#v#h2#t3 | Implemented; 4 | MH₂Cl(PR₃)₃ and IrH₂Cl(CO)(PPh₃)₂: explicit classical III dihydrides formed by oxidative addition, HW 2013, 2017. All H splits with one ordinary Cl; all-donor M³+ forms reuse Co16. |
| Co21b | Rh#n3#v3#t2; #n3#v2#h#t2; #n3#v#h2#t2; #c+3#n3#t5; #c-#n3#v4#t | Implemented; 5 | Five-coordinate RhH₂Cl(PPh₃)₂ (2013), and the acyl Rh(III) intermediate in the Monsanto cycle (2021). Classical H/Cl or acyl/I bonds plus neutral donors; five-donor alternative. No new type for the six-coordinate oxidative-addition intermediate, which reuses Co16. |
| Co22 | Rh#n4#h#t3 | Implemented; 1 | RhH(CO)(PPh₃)₂, HW 2013: paired Rh(I) monohydride, explicit/implicit H. H− donors reuse Co19a. |
| Co23 | Co#c-3#n3#v5#h | Implemented; 1 | MH(CN)₅³−, HW 1994 and 2014. Implicit-H counterpart; explicit H and H− donors reuse the paired III states above. |
| Co24 | Ir#n2#u#v3#h#t2; #n2#u#v2#h2#t2 | Implemented; 2 | IrH₂Cl₂(PiPr₃)₂, HW 2068, explicitly a rare paramagnetic 17-electron hydride. Ordinary H/Cl with two P donors; all H splits. The IV all-donor state already occurs in Co17a. |
| Co25a | Rh#n#u#v6 | Implemented; 1 | RhO₄²−: HW 2014 gives one unpaired electron. Two Rh=O plus two Rh–O− bonds retain that occupation; the proposed tetrahedral geometry is qualified in HW. |
| Co25b | Rh#c-2#n#u#v8 | Implemented; 1 | The same RhO₄²− with four ordinary Rh=O bonds and charge on Rh. Counts 19 electrons; Rh capacity 18 → 19 is approved. This is a charge-allocation alternative, not a different spin assignment. |
| Co26 | Co/Rh/Ir#c-#n5#t4; #c-3#n6#t3 | Implemented; 6 | M(CO)₄− and the superreduced M(CO)₃³−, HW 2131–2135. The mononuclear carbonylate table and closed-shell description support these paired 18-electron states. No cluster-specific MC state. |
| Co27 | Co/Rh#n2#v5#t2 | Implemented; 2 | Co₂(CO)₈ and Co₄/Rh₄(CO)₁₂, HW 2112–2114, Fig. 32.28: ordinary M–M bonds and C(=O)(–M)₂ bridges. Terminal CO remains a donor. The #v3 state is reused; basal bridge-bearing tetrahedral-cluster sites add #v5#t2. All-terminal Ir₄(CO)₁₂ reuses Ir#n3#v3#t3. |
| Co28 | Reuses preceding rows | Implemented; 0 | Rh₂(H₂O)₁₀⁴+ and common Rh₂(O₂CR)₄L₂: explicit diamagnetic Rh–Rh single bond, HW 2018. Five donors per Rh, or two ordinary carboxylate O bonds and three donors. Mixed-valent Rh₂⁵+ remains in Co35 until its bond-carried charge/spin graph is worked through. |
| Co29 | Co#n3#u2#v; #n3#u2#h; #c+#n3#u2#t | Implemented; 3 | Molecular CoH: G05, pp. 2–3 and 51–65, establishes X³Φ and gives one bonding pair plus three nonbonding pairs and two unpaired electrons. Explicit/implicit ordinary H and H− donation. |
| Co30 | Co#n#u4#v3#t3 | Implemented; 1 | High-spin CoF₃(H₂O)₃: N84, p. 2722, explicitly lists it with CoF₆³−, citing Cotton/Wilkinson (1980), p. 775; its own NMR experiment concerns a different compound. This supports HW 2000 over the contradictory low-spin wording on 1994. Three ordinary F bonds and three water donors; all-donor form reuses Co08. |
| Co31a | Co/Rh/Ir#n3#u#v2; #n3#u#v#h; #n3#u#h2; #c+2#n3#u#t2 | Implemented; 12 | Matrix-isolated MH₂: VH92, pp. 8123, 8126–8129, establishes bent S = 1/2 species with predominantly metal-local spin; WA02 independently supports RhH₂ (²A₁). All H splits and H− donation. CoH₂ is the observed doublet, not a claim to have resolved its conflict with calculated quartet ground states; no prospective quartet H variants are selected. |
| Co31b | Rh#n3#u2#v; #n3#u2#h; #c+#n3#u2#t | Implemented; 3 | Molecular RhH: WA02, pp. 3706–3709, Table 2, gives the observed X³Δ molecule with three pairs and two unpaired electrons after the ordinary Rh–H bond. Explicit/implicit H and H− donation. Its distinct H₂ adducts do not justify further hydride states. |
| Co31c | Ir#n3#u2#v; #n3#u2#h; #c+#n3#u2#t | Implemented; 3 | Molecular IrH: L14, pp. 39–41, finds a predominantly two-open-shell X³Φ₄ configuration from the hyperfine constants. Use three nonbonding pairs and two unpaired electrons, with explicit/implicit ordinary H or H− donation (PW). The reported four-open-shell admixture does not add a separate occupation. |
| Co31d | RhH₃ and Rh hydride–H₂ adducts | No new states selected | WA02 observes (H₂)RhH₃, not bare RhH₃; the other intact-H₂ adducts also do not establish additional classical hydride states. H₂-donor accounting remains outside this pass. |
| Co32a | Co#n2#u3#v2 | Implemented; 1 | Molecular CoO with an ordinary Co=O bond: C97, pp. 8795–8796, 8798, retains the experimentally established X⁴Δ state. Two nonbonding pairs and three unpaired electrons; the erroneous calculated sextet ordering is not adopted. |
| Co32b | Co/Rh/Ir#n2#u#v4 | Implemented; 3 | Molecular O=M=O: VO92, pp. 7247, 7249–7251, ESR establishes doublets with predominantly metal-local spin and little oxygen spin. Two ordinary M=O bonds. Later bent RhO₂ calculations in G11 change no atom count; these are distinct from peroxide/superoxide isomers. |
| Co32c | Ir#u#v8 | Implemented; 1 | IrO₄ with four ordinary Ir=O bonds: G09Ir, pp. 7880–7882, identifies the matrix-isolated D₂d ²A₁ molecule and explicitly gives d¹ Ir(VIII). Counts 17 electrons; no capacity change. |
| Co32d | Rh#c+2#n#u#v4#t2; Ir#n#u#v6; #c+2#n#u#v4#t2 | Implemented; 3 | Observed (η²-O₂)MO₂ peroxide isomers: G09Ir, p. 7881, and G11, pp. 672, 675, explicitly assign ²A₂ and a closed-shell peroxide plus two oxo groups. Two ordinary M=O bonds and two M–O bonds, or two peroxide O− donors. Rh's all-ordinary state reuses Co25a; existing O#c-#n2#v#d suffices. |
| Co32e | Rh#c+4#n2#u#t4 | Implemented; 1 | Observed D₂d bis-peroxide Rh(η²-O₂)₂: G11, pp. 672–674, assigns ²B₁, distinct from its D₂h quartet bis-superoxide isomer. Two O–O single bonds and four ordinary Rh–O bonds reuse Co32b; four peroxide O− donors add the displayed state. No MC bond. |
| Co32f | Rh#c-#n3#v4 | Implemented; 1 | Observed RhO₂− with two ordinary Rh=O bonds and charge on Rh: G11, Fig. 7b and Table 6, assigns ¹A₁. Three nonbonding pairs; no prospective triplet alternative. |
| Co32g | Cyclic Co(O₂) | Covered by implemented types; 0 | C97, pp. 8796, 8798–8799, establishes a peroxide-like ring but leaves doublet/quartet ordering unresolved. An O–O single bond and two ordinary Co–O bonds use Co31a/Co32a for the respective supplied occupations. No additional prospective occupation or donor state is selected. |
| Co32h | CoOO; (η¹-O₂)MO₂ (M = Rh, Ir); Rh(η²-O₂) and the D₂h bis-superoxide Rh(η²-O₂)₂ | Open local occupations; no new states selected | C97/G09Ir/G11 establish these constitutions and discuss total multiplicities, but do not fix the metal/ligand spin allocation. CA99 supplies charges and total multiplicities, not the missing local spins; its earlier bis-dioxygen interpretation does not replace G11's distinction between isomers. G09Ir SI Fig. S9 shows spin density only for IrO₄, not the end-on isomer. |
| Co32i | Molecular RhO₃ and IrO₃ | Open; occupation sources still to locate | HW 2014 identifies the volatile trioxides. The reviewed papers do not settle the three-terminal-oxo occupations; oxide-peroxide MO₄ isomers are not substitutes. No hypothetical RhO₄ state is selected: G11 reports no experimental identification. |
| Co32j | Rh#n2#u3#v2 | Implemented; 1 | Molecular RhO: CA99, pp. 4846–4847, identifies the matrix-isolated monoxide and gives the experimentally supported X⁴Σ− ground state. An ordinary Rh=O graph assigns two nonbonding pairs and three unpaired electrons to Rh (PW). |
| Co33a | CoO₄³− | Open; occupation source still to locate | HW 2004–2005 identifies the Co(V) anion. Co(I)/Co(IV) occupations are settled separately in Co33f/Co33g; neither establishes this d⁴ occupation. |
| Co33b | RhO₄³− | Open; occupation source still to locate | HW 2014 identifies the orange Rh(V) solution and its disproportionation but does not specify the occupation. The one-unpaired-electron measurement there concerns RhO₄²−, not this ion. |
| Co33c | CoN₂⁵− | Triplet occupation supported; bond allocation open | B19, §VII, explicitly gives Co(I) d⁸, three pairs and two unpaired electrons, in LiSr₂[CoN₂]/LiCa₂[CoN₂]. The paper does not assign integer Co–N bond orders or donor contributions; no atom state has yet been selected from that occupation alone. |
| Co33d | CoN₃⁶− in reported Ca₃CoN₃ | Not selected; identity challenged | HW 2000 lists the planar anion. B07, p. 165, reassigns the earlier reported Ca₃CoN₃ powder pattern to Ca₅[CoN₂]₂ with CaO impurities. No prospective Co(III) state is selected from this disputed identification; the independently established CoN₂⁵− family remains in Co33c. |
| Co33e | CoO₃⁴−, CoO₄⁶− and Co₄O₉¹⁰− | Covered by implemented types; 0 | M98, pp. 3198–3201, establishes d⁷ quartet Co(II) in the isolated oxoanions and supports the same local occupation in the linked Co₄O₉¹⁰− units. One Co=O per Co with remaining ordinary Co–O bonds in CoO₃⁴−/Co₄O₉¹⁰−, or four Co–O− bonds in CoO₄⁶−, uses Co#c-2#n2#u3#v4 from Co03 (PW graph allocation). All-double-bond CoO₃⁴− reuses Co04 within the approved Co capacity. |
| Co33f | Co#c-3#n3#u2#v4 | Implemented; 1 | CoO₂³−: S01 establishes the near-linear anion in Rb₃CoO₂; S03Co, §3.6.4, pp. 97–98, explicitly fits S = 1 Co(I), g = 2.19, μeff = 3.10 μB. Two ordinary Co=O bonds as in HW 1999 give three pairs and two unpaired electrons, with charge on Co (PW). Counts 16 electrons; existing O states suffice. |
| Co33g | Co#u5#v4 | Implemented; 1 | CoO₄⁴− in Na₄CoO₄: K15, manuscript pp. 1–3/Figs. 3–4, explicitly gives tetrahedral d⁵ S = 5/2. Four ordinary Co–O− bonds preserve that occupation (PW); existing O#c-#n3#v suffices. Counts 13 electrons; no capacity change. |
| Co34 | CoH₂(PPh₃)₂ | Open; source for this exact formulation still to locate | HW 1992 lists this phosphine dihydride. The sourced CoH₃(PR₃)₃ and CoH₂(PMe₃)₄+ compounds do not establish its occupation or hydride/H₂ partition. |
| Co34a | Co#c+#n3#v#h#t4; #c+#n3#h2#t4 | Implemented; 2 | cis-CoH₂(PMe₃)₄+: P90, pp. 1043, 1045–1046, assigns two terminal hydrides from X-ray, IR/isotope and NMR evidence. The paired six-coordinate Co(III) occupation follows HW 1992/2000–2002 and reuses Co07 for explicit H and H− donation; only mixed/implicit-H variants are new. This assignment is specific to PMe₃, not every CoH₂(PR₃)₄+ formula. |
| Co34b | Ir#c+#n3#v#h#t4; #c+#n3#h2#t4 | Implemented; 2 | IrH₂(CO)₂(PPh₃)₂+: HC88, Table II, p. 4128 and p. 4130, explicitly establishes a classical octahedral d⁶ dihydride. Paired Ir(III) follows HW 2017; explicit H and H− donation reuse Co16/Co11c. Only mixed/implicit-H variants are new. |
| Co34c | Ir#n2#v5#t2; #n2#v4#h#t2; #n2#v3#h2#t2; #n2#v2#h3#t2; #n2#v#h4#t2; #n2#h5#t2; #c+5#n2#t7 | Implemented; 7 | IrH₅(PR₃)₂, HW 2010: HC88, Tables I–II, establishes the classical PCy₃ complex; JA16, §§2–3.1 and Fig. 1A/D, uses closed-shell calculations for the neutral PPh₃/PMe₃ complexes. Its bare-IrH₅ model also recovers five localized 2c2e bonds (Fig. 2A). The ordinary graph leaves two nonbonding pairs on Ir (PW); all H-count splits and H− donation count 18 electrons. |
| Co34d | Co#n3#v2#h#t3; #n3#v#h2#t3; #n3#h3#t3 | Implemented; 3 | CoH₃(PR₃)₃, PR₃ = PPh₃, PEtPh₂, PEt₂Ph: S68, p. 130/Table V, explicitly assigns octahedral Co(III) with a trans H–Co–H pair from IR/H–D substitution; SV19, Table 3, retains the classical hydride treatment. Reuses Co07's paired III occupation, explicit-H and H−-donor states; three mixed/implicit-H variants are new. |
| Co34e | Ir#n3#h3#t3 | Implemented; 1 | fac-IrH₃(PPh₂Me)₃: B93, pp. 3359, 3361–3362, directly locates three terminal H atoms by neutron diffraction, Ir–H = 1.627(4) Å. L90, p. 858/Table VII, independently classifies IrH₃(PMe₂Ph)₃ as classical (T₁ = 465 ms). Reuses paired octahedral Ir(III), Co16/Co21a; only the all-implicit-H state is new. |
| Co34f | IrH₂(η²-H₂)(PMe₂Ph)₃+ | Assignment supported; H₂-donor accounting deferred; 0 | L90, pp. 858–859/Table VII, assigns two terminal hydrides and one intact H₂ ligand from T₁ = 20 ms and reactivity. Retain the H₂ donor description under the shared deferral; this does not support a four-terminal-hydride state. |
| Co35 | Mixed-valent Rh₂⁵+ carboxylates | Open; graph check | HW 2018 gives the Rh₂⁵+ bond occupation; work through its bond-carried charge/spin allocation before deciding whether another type is needed. |
| Co35a | Rh/Ir#c+2#n3#u#t6 | Implemented; 2 | Monomeric RhCp₂ and IrCp₂: KW67, pp. 186–187, establishes S = 1/2. Formal M(II) with two six-electron Cp− donors uses the same graph allocation as Co06 (PW), counting 19 electrons; Rh/Ir capacity 18 → 19 is approved. KW67 discusses ligand character of the singly occupied orbital; #u on M is the formal allocation, not a claim that the spin density is confined to M. |
| Co36 | Co#c+#n3#u2#t6 | Implemented; 1 | Six-coordinate Co(I) polypyridyl complexes Co(tpy)₂+, Co(tbpy)₃+ and Co(phen)₃+: E15, pp. 12005–12007, Table 4 and p. 12016, establishes neutral ligands and d⁸ S = 1 by crystallography, magnetic measurements and electronic-structure analysis. Six N donors give 20 electrons; covered by Co capacity 23. No ligand-radical or MC state is needed. |
| Co36a | Co#n4#u#t4 | Implemented; 1 | Co(PMe₃)₄: K71Co, p. 343, establishes the monomer and explicitly assigns tetrahedral d⁹ from μeff = 1.70 μB and EPR. Four P donors, four nonbonding pairs and one unpaired electron; 17 electrons. Distinct from the paired HCo(PMe₃)₄ in Co10. |
| Co36b | Rh#n4#u#t4 | Implemented; 1 | Rh(dppe)₂: K24, pp. 5–6 and Figs. 3–4, treats the observed paramagnetic Rh(0) species as d⁹ with one singly occupied orbital. Four P donors give the displayed 17-electron state (PW). The +I/−I complexes reuse Co19a/Co26; partial spin delocalization does not require a new MC type. |
| Co37 | Co(1-norbornyl)₄/CoR₄+, high-nuclearity carbonyl/hydride clusters, specialized Rh/Ir oxo/nitrido assemblies | Not selected for this pass | HW 2005, 2015–2019, 2115–2117, 2137–2138. These do not justify additional scaffold-specific default states. Small ordinary-bond carbonyl clusters are covered in Co10/Co27. Hypothetical IrF₇ is also excluded (2011–2012). |
| Co38 | Existing Co#c+#n3#u#v | Open; 1 existing state | Corrected the old CoX− comment to CoX+. A supporting bare low-spin CoX+ example has not yet been established. Do not silently retain an unsupported citation, substitute a ligated example, or remove the state before completing that check. |

#### Implementation and coverage

Co capacity 23 covers the eight-donor tetranitrate graph (Co14); the seven-donor
form counts 21, Co04/Co06 count 19, and the six-donor Co(I) triplet counts 20.
Rh/Ir capacity 19 covers the formal M(II) metallocene donors; Rh also needs 19 for
the all-double-bond RhO₄²− graph. These approved increases change graph-count limits,
not free-ion configurations. ValenceTable remains unchanged: fixed H counts use
the existing Counts route for these elements.

| Coverage | Result |
| --- | --- |
| New fixtures | 199: donor_acceptor 112, hydrides 51, inorganic_ions 18, atomic_ions 9, inorganic_small 9 |
| Added states | All 195; each input resolves to the expected state under AtomTyping/Counts × Strict/MostSaturated and passes ValenceInvariantsValidator |
| Retained states | All 12 supported states covered; Co38 retained with its source gap explicit |
| Existing inputs | Rh²⁺/Rh³⁺ low-spin inputs and cobalt porphine now supply their intended spin; their existing snapshots are unchanged. Separate free-ion inputs retain ground-state occupations |
| Corrin fixtures | Peripheral substituents replaced by H, following L09's computational model; Co(I/II) occupations and donor allocations follow Co13 |
| Regression checks | 2,902 available conformance cases, 3,594 graph unit tests, 223 graph property tests and 310 chemistry tests pass; formatting and Clippy pass |

Open source/graph cases remain Co32h–i, Co33a–c, Co34, Co35 and Co38.
The H₂-donor descriptions remain deferred; no unsupported assignments were added.

#### Cobalt-group sources and access

| Key | Source | Contribution / access |
| --- | --- | --- |
| ASD | NIST Atomic Spectra Database, ground configurations/terms. Archived query outputs: [Co](../materials/valence/NIST-ASD-2026-co-ground-states.html), [Rh](../materials/valence/NIST-ASD-2026-rh-ground-states.html), [Ir](../materials/valence/NIST-ASD-2026-ir-ground-states.html); [live database](https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html) | Atomic and free-ion occupations in Co01/Co02, checked independently of complex spin states. |
| K12 | Krzystek et al., “High Spin Co(I): High-Frequency and -Field EPR Spectroscopy of CoX(PPh₃)₃ (X = Cl, Br),” *Inorg. Chem.* **51** (2012), 4954–4964. [DOI](https://doi.org/10.1021/ic202185x), [author copy](https://nationalmaglab.org/library/publications/NHMFL_Publication-6133.pdf), [PDF](../materials/valence/Krzystek-2012-cobaltI.pdf) | Retrieved and reviewed. Explicit triplet phosphine/halide complexes; paired corrin and five-coordinate Co(I) comparison. Co12/Co13a. |
| C16 | Chen et al., “Slow Magnetic Relaxations in Cobalt(II) Tetranitrate Complexes,” *Inorg. Chem.* **55** (2016), 12603–12617. [DOI](https://doi.org/10.1021/acs.inorgchem.6b01544), [repository copy](https://www.osti.gov/servlets/purl/1341566), [PDF](../materials/valence/Chen-2016-cobalt-tetranitrate.pdf) | Retrieved and reviewed. Manuscript pp. 11–12: seven/eight coordination; pp. 14–18: S = 3/2 and spectroscopic magnetic analysis. Co14. |
| G05 | Gordon, *Spectroscopy of selected metal-containing diatomic molecules*, Waterloo PhD thesis (2005). [Repository copy](https://dspacemainprd01.lib.uwaterloo.ca/server/api/core/bitstreams/17bd1ee2-64d8-4615-8f65-7f7ae4acc04a/content), [PDF](../materials/valence/Gordon-2005-metal-hydrides-thesis.pdf) | Relevant CoH sections reviewed; pp. 2–3 give the orbital occupation, pp. 51–65 report the spectroscopy. Co29. |
| L09 | Liptak et al., “Spectroscopic and Computational Characterization of the Base-off Forms of Cob(II)alamin,” *J. Phys. Chem. B* **113** (2009), 5245–5254. [DOI](https://doi.org/10.1021/jp810136d), [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC3127260/), [archived HTML](../materials/valence/Liptak-2009-cobII.html) | Full text retrieved and relevant EPR/occupation sections reviewed. Explicit four-/five-coordinate low-spin Co(II), Co13b. The PDF request returned a challenge page; the complete HTML was available. |
| N84 | Navon and Kläui, “⁵⁹Co NMR of a cobalt(III) spin-crossover compound,” *Inorg. Chem.* **23** (1984), 2722–2725. [DOI](https://doi.org/10.1021/ic00185a036), [PDF](../materials/valence/Navon-Klaui-1984-cobalt-spin.pdf) | Reviewed and archived. Introductory high-spin assignment for CoF₃(H₂O)₃, p. 2722, citing Cotton/Wilkinson; not an experimental study of that hydrate. Co30. |
| VH92 | Van Zee et al., “Electron-spin resonance of the Co, Rh, and Ir dihydride molecules,” *J. Chem. Phys.* **97** (1992), 8123–8129. [DOI](https://doi.org/10.1063/1.463433), [PDF](../materials/valence/VanZee-1992-dihydrides.pdf) | Reviewed and archived. ESR-observed metal-centered doublets for all three MH₂ molecules; CoH₂ conflict with calculated quartet ground states is retained explicitly. Co31a. |
| WA02 | Wang and Andrews, “Infrared Spectra of Rhodium Hydrides in Solid Argon, Neon, and Deuterium with Supporting Density Functional Calculations,” *J. Phys. Chem. A* **106** (2002), 3706–3713. [DOI](https://doi.org/10.1021/jp013624f), [PDF](../materials/valence/Wang-Andrews-2002-rhodium-hydrides.pdf) | Reviewed and archived. RhH triplet, RhH₂ doublet, and separation of bare hydrides from intact-H₂ adducts, pp. 3706–3709, 3711–3713. Co31a/Co31b/Co31d. |
| A12 | Adam et al., “Laser spectroscopy of iridium monohydride and iridium monodeuteride,” *Chem. Phys. Lett.* **535** (2012), 21–25. [DOI](https://doi.org/10.1016/j.cplett.2012.03.038), [PDF](../materials/valence/Adam-2012-iridium-hydride.pdf) | Reviewed and archived. Ω = 4 ground level; the discussion on pp. 23–24 leaves the two-/four-open-shell mixture unresolved. L14 is the experimental follow-up, Co31c. |
| VO92 | Van Zee et al., “Cobalt, Rhodium, and Iridium Dioxide Molecules and Walsh-Type Rules,” *J. Phys. Chem.* **96** (1992), 7247–7251. [DOI](https://doi.org/10.1021/j100197a022), [PDF](../materials/valence/VanZee-1992-dioxides.pdf) | Reviewed and archived. ESR assigns doublet O=M=O molecules with predominantly metal-local spin; does not assign peroxide/superoxide isomers. Co32b. |
| C97 | Chertihin et al., “Reactions of Laser-Ablated Cobalt Atoms with O₂: Infrared Spectra of Cobalt Oxides in Solid Argon,” *J. Phys. Chem. A* **101** (1997), 8793–8802. [DOI](https://doi.org/10.1021/jp972022m), [PDF](../materials/valence/Chertihin-1997-cobalt-oxides.pdf) | Reviewed and archived. CoO quartet; separate OCoO, cyclic peroxide and end-on superoxide descriptions, pp. 8795–8799. Co32a/Co32g/Co32h. |
| G09Ir | Gong et al., “Formation and Characterization of the Iridium Tetroxide Molecule with Iridium in the Oxidation State +VIII,” *Angew. Chem. Int. Ed.* **48** (2009), 7879–7883. [DOI](https://doi.org/10.1002/anie.200902733), [PDF](../materials/valence/Gong-2009-iridium-tetroxide.pdf), [SI](../materials/valence/Gong-2009-iridium-tetroxide-SI.pdf) | Main paper and SI reviewed and archived. Explicit d¹ IrO₄ and doublet peroxide (η²-O₂)IrO₂, p. 7881, Co32c/Co32d. SI Table S3 confirms the peroxide doublet below the quartet; Fig. S9 supplies spin density only for IrO₄. No local occupation is added for (η¹-O₂)IrO₂. |
| G11 | Gong et al., “Is rhodium tetroxide in the formal oxidation state VIII stable? A quantum chemical and matrix isolation investigation of rhodium oxides,” *Theor. Chem. Acc.* **129** (2011), 667–676. [DOI](https://doi.org/10.1007/s00214-011-0919-7), [PDF](../materials/valence/Gong-2011-rhodium-oxides.pdf) | Reviewed and archived. Doublet oxo-peroxide and bis-peroxide, singlet RhO₂−; separate superoxide isomers and no experimentally identified RhO₄. Fig. 7, pp. 672–675; Co32d–Co32i. |
| L14 | Linton et al., “Hyperfine structure and Stark effect in the [18.0] 4–X³Φ₄ transition of iridium hydride, IrH and [18.2] 4–X³Φ₄ transition of iridium deuteride, IrD,” *J. Mol. Spectrosc.* **304** (2014), 37–42. [DOI](https://doi.org/10.1016/j.jms.2014.08.003), [PDF](../materials/valence/Linton-2014-iridium-hydride.pdf) | Reviewed and archived. Hyperfine analysis and conclusions, pp. 39–41, identify a predominantly two-open-shell configuration with significant four-open-shell admixture. Co31c uses the dominant occupation; the admixture is not a second registry state. |
| B19 | Ballé et al., “Ferromagnetic ordering of linearly coordinated Co ions in LiSr₂[CoN₂],” *Phys. Rev. B* **99** (2019), 094422. [DOI](https://doi.org/10.1103/PhysRevB.99.094422), [author manuscript](https://arxiv.org/abs/1811.00867), [PDF](../materials/valence/Balle-2019-nitridocobaltate-author-manuscript.pdf) | Retrieved and relevant structural/magnetic sections reviewed. §VII, manuscript pp. 12–14, gives the d⁸ triplet ground configuration. Co33c. |
| M98 | Möller, “Synthesis, Structures, Magnetic Properties, and Absorption Spectra of the Alkaline Oxocobaltates(II): Li₆CoO₄, Na₄CoO₃, and Na₁₀Co₄O₉,” *Chem. Mater.* **10** (1998), 3196–3201. [DOI](https://doi.org/10.1021/cm9802907), [PDF](../materials/valence/Moller-1998-oxocobaltates.pdf) | Reviewed and archived. Structural, absorption and magnetic analysis, pp. 3198–3201, supports high-spin d⁷ Co(II) in tetrahedral/trigonal-planar oxide environments; Co33e reuses proposed types. Does not settle the Co(I/IV/V) anions in Co33a. |
| K71Co | Klein, “Tetrakis(trimethylphosphane)cobalt(0): Preparation and Reactions,” *Angew. Chem. Int. Ed. Engl.* **10** (1971), 343. [DOI](https://doi.org/10.1002/anie.197103431), [PDF](../materials/valence/Klein-1971-cobalt0-phosphine.pdf) | Reviewed and archived. Molecular-weight, magnetic and EPR evidence explicitly supports monomeric tetrahedral d⁹ Co(PMe₃)₄; Co36a. |
| CA99 | Citra and Andrews, “Reactions of Laser-Ablated Rhodium Atoms with O₂. Infrared Spectra and DFT Calculations for RhO, ORhO, (O₂)RhO₂, Rh₂O₂, Rh(O₂) and (O₂)Rh(O₂) in Solid Argon,” *J. Phys. Chem. A* **103** (1999), 4845–4854. [DOI](https://doi.org/10.1021/jp990864a), [PDF](../materials/valence/Citra-Andrews-1999-rhodium-oxides.pdf) | Reviewed and archived. RhO quartet, pp. 4846–4847, supports Co32j; peroxide discussion, pp. 4850–4851, corroborates Co32d. Charges and total multiplicities do not settle Co32h's local spins. G11 supplies the later distinction between bis-peroxide/bis-superoxide isomers. |
| B07 | Bendyna, Höhn and Kniep, “Crystal structure of pentacalcium bis(dinitridocobaltate(I)), Ca₅[CoN₂]₂, and a note on ‘Ca₃CoN₃’,” *Z. Kristallogr. NCS* **222** (2007), 165–166. [DOI](https://doi.org/10.1524/ncrs.2007.0068), [PDF](../materials/valence/Bendyna-2007-dinitridocobaltate.pdf) | Reviewed and archived. Reassigns the 1993 Ca₃CoN₃ powder pattern to Ca₅[CoN₂]₂/CaO, p. 165; Co33d. |
| K24 | Kick et al., “Understanding Ligand Effects on Bielectronic Transitions: Chemo- and Electroreduction of Rhodium Bis(Diphosphine) Complexes to Low Oxidation States,” *Angew. Chem. Int. Ed.* **63** (2024), e202408356. [DOI](https://doi.org/10.1002/anie.202408356), [PDF](../materials/valence/Kick-2024-rhodium-diphosphines.pdf) | Relevant sections reviewed and PDF archived, pp. 5–6/Figs. 3–4: d⁹ Rh(0), a singly occupied orbital, and prior EPR observation; Co36b. |
| E15 | England et al., “Molecular and Electronic Structures of Homoleptic Six-Coordinate Cobalt(I) Complexes of 2,2′:6′,2″-Terpyridine, 2,2′-Bipyridine, and 1,10-Phenanthroline. An Experimental and Computational Study,” *Inorg. Chem.* **54** (2015), 12002–12018. [DOI](https://doi.org/10.1021/acs.inorgchem.5b02415), [PDF](../materials/valence/England-2015-cobaltI-polypyridyl.pdf) | Relevant structural, magnetic and electronic-structure sections reviewed; PDF archived. Neutral polypyridyl ligands and Co(I) d⁸ triplet, pp. 12005–12007, Table 4 and p. 12016; Co36. Co(II/III) results are consistent with the occupations already represented in Co04/Co06/Co07. |
| P90 | Pérés et al., “Reaction of CH₂O with [Co(PMe₃)₄]X (X = PF₆ and BPh₄). Structural Characterization of [Co(CO)(PMe₃)₄]PF₆ and the Nonstoichiometric Complex [CoH₂₋₂ₓ(CO)ₓ(PMe₃)₄]BPh₄,” *Organometallics* **9** (1990), 1041–1047. [DOI](https://doi.org/10.1021/om00118a021), [PDF](../materials/valence/Peres-1990-cobalt-dihydride-carbonyl.pdf) | Reviewed and archived. X-ray and IR/isotope/NMR evidence supports cis-octahedral Co(III) with two terminal hydrides, pp. 1043, 1045–1046. The paired occupation uses HW’s existing six-coordinate Co(III) description; Co34a. |
| HC88 | Hamilton and Crabtree, “An NMR Method for Distinguishing Classical from Nonclassical Structures in Transition-Metal Polyhydrides,” *J. Am. Chem. Soc.* **110** (1988), 4126–4133. [DOI](https://doi.org/10.1021/ja00221a005), [PDF](../materials/valence/Hamilton-Crabtree-1988-polyhydride-NMR.pdf) | Relevant Ir hydride sections reviewed; PDF archived. Tables I–II and p. 4130 establish classical IrH₅(PCy₃)₂ and IrH₂(CO)₂(PPh₃)₂+. The dihydride reuses paired octahedral Ir(III); JA16 supplies the computational paired description for the pentahydride family. Co34b/Co34c. |
| JA16 | Jimenez-Izal and Alexandrova, “σ-Aromaticity in polyhydride complexes of Ru, Ir, Os, and Pt,” *Phys. Chem. Chem. Phys.* **18** (2016), 11644–11652. [DOI](https://doi.org/10.1039/C5CP04330A), [advance-publication PDF](../materials/valence/Jimenez-Izal-2016-polyhydrides.pdf) | Relevant methods, hydride analysis and figures reviewed; PDF archived. §§2–3.1/Fig. 1A,D give closed-shell neutral IrH₅(PR₃)₂ calculations; the explicitly singlet bare-IrH₅ model also supports five localized 2c2e bonds (Fig. 2A). Co34c uses the paired ordinary-bond graph (PW); no aromatic/MC type is required. Main-text evidence suffices; SI not needed for this assignment. |
| S01 | Sofin and Jansen, “Rb₃CoO₂, a Novel Oxocobaltate(I) Synthesized via the Azide/Nitrate Route,” *Z. Anorg. Allg. Chem.* **627** (2001), 2115–2117. [DOI](https://doi.org/10.1002/1521-3749%28200109%29627%3A9%3C2115%3A%3AAID-ZAAC2115%3E3.0.CO%3B2-L), [PDF](../materials/valence/Sofin-2001-rubidium-oxocobaltate.pdf) | Supplied scan/OCR reviewed against the pages and archived. Near-linear CoO₂³−, Table 3, p. 2117; no magnetic measurement in this paper. S03Co supplies the occupation for Co33f. |
| S03Co | Sofin, *Neue Alkalioxometallate über die Azid/Nitrat-Route*, Stuttgart PhD thesis (2003). [Repository PDF](https://elib.uni-stuttgart.de/server/api/core/bitstreams/50ed890c-50bd-4ae2-b745-57662b32ef90/content), [PDF](../materials/valence/Sofin-2003-alkali-oxometallates-thesis.pdf) | Retrieved and archived; §3.6.4, pp. 97–98, reviewed. Impurity-corrected susceptibility of Rb₃CoO₂ is fitted as high-spin Co(I), S = 1, with g = 2.19 and μeff = 3.10 μB. The same d⁸ triplet assignment is explicit for the CoO₂³− unit in Na₅[CoO₂]CO₃, §IV.6.4, p. 117. Co33f. |
| K15 | Kanungo et al., “Na₄IrO₄: Square-Planar Coordination of a Transition Metal in d⁵ Configuration due to Weak On-Site Coulomb Interactions,” *Angew. Chem. Int. Ed.* (2015). [DOI](https://doi.org/10.1002/anie.201411959), [author manuscript](https://arxiv.org/abs/1502.02275), [PDF](../materials/valence/Kanungo-2015-iridate-cobaltate-author-manuscript.pdf) | Retrieved and archived. The Na₄CoO₄ comparison explicitly assigns high-spin d⁵, S = 5/2; manuscript pp. 1–3/Figs. 3–4. Co33g uses that occupation with ordinary Co–O− bonds (PW). |
| S68 | Sacco and Rossi, “Nitrogen fixation. I. Hydrido- and nitrogen-complexes of cobalt,” *Inorg. Chim. Acta* **2** (1968), 127–132. [DOI](https://doi.org/10.1016/S0020-1693(00)87009-7), [PDF](../materials/valence/Sacco-1968-cobalt-hydrides.pdf) | Reviewed and archived. CoH₃(PR₃)₃ characterization and H/D infrared data, p. 130/Table V, support an octahedral classical trihydride. The explicit diamagnetism statement on p. 132 concerns CoH(N₂)(PPh₃)₃, not the trihydride. Co34d reuses the accepted paired Co(III) occupation. |
| SV19 | Schnieders et al., “Metal Hydride Vibrations: The Trans Effect of the Hydride,” *Inorg. Chem.* **58** (2019), 12467–12479. [DOI](https://doi.org/10.1021/acs.inorgchem.9b02302), [author manuscript](https://hdl.handle.net/1807/105558), [PDF](../materials/valence/Schnieders-2019-hydride-vibrations-author-manuscript.pdf) | Retrieved and archived; Table 3 and accompanying discussion reviewed. CoH₃(PRAr₂)₃/CoH₃(PAr₃)₃ and IrH₃(PAr₃)₃ are treated with classical trans-hydride force constants. Corroborates Co34d/Co34e; no separate spin measurement claimed. |
| B93 | Bau et al., “Neutron diffraction study of fac-[IrH₃(PPh₂Me)₃]·MeOH,” *J. Chem. Soc., Dalton Trans.* (1993), 3359–3362. [DOI](https://doi.org/10.1039/DT9930003359), [PDF](../materials/valence/Bau-1993-iridium-trihydride.pdf) | Reviewed and archived. Neutron diffraction directly locates three normal terminal hydrides, with distorted octahedral coordination and Ir–H = 1.627(4) Å. Co34e. |
| L90 | Lundquist et al., “Reactivity of the Molecular Hydrogen Complex [IrH₄(PMe₂Ph)₃]BF₄ toward Olefins,” *J. Am. Chem. Soc.* **112** (1990), 855–863. [DOI](https://doi.org/10.1021/ja00158a054), [PDF](../materials/valence/Lundquist-1990-iridium-dihydrogen.pdf) | Reviewed and archived. pp. 858–859/Table VII distinguish classical IrH₃(PMe₂Ph)₃ (T₁ = 465 ms) from IrH₂(η²-H₂)(PMe₂Ph)₃+ (20 ms), supported by H₂ loss/reactivity. Co34e/Co34f. |
| KW67 | Keller and Wawersik, “Spektroskopische Untersuchungen an Komplexverbindungen VI. EPR-Spektren von (C₅H₅)₂Rh und (C₅H₅)₂Ir,” *J. Organomet. Chem.* **8** (1967), 185–188. [DOI](https://doi.org/10.1016/S0022-328X(00)84718-X), [PDF](../materials/valence/Keller-1967-metallocene-EPR.pdf) | Reviewed and archived. pp. 186–187 assign molecular doublets to monomeric RhCp₂/IrCp₂; the final orbital discussion suggests substantial ligand character. Co35a uses formal M(II) donor counting (PW), independently of the spatial spin distribution. |

Supplied PDFs are archived in `materials/valence/`.
Extraction text, download logs and the deduplicated candidate/electron-count audit
are in `scratch/hw-cobalt-review/`.

### Nickel group

HW II, chapter XXXI, pp. 2023–2065, with hydride, carbonyl and metallocene
cross-references in chapter XXXII. **173 additions implemented: Ni 59, Pd 44, Pt 69,
and one O companion; all seven existing metal states retained and covered.** Eight
states cover three local hydride-bridge patterns. Source review, registry changes
and conformance coverage are complete. Ni/Pd capacity is 20; Pt retains 18.
PdF₆− is omitted because its occupation is unresolved.

#### Nickel-group work table

Element prefixes apply to each suffix in their cell. Counts exclude existing states
and reuse across rows, and include all displayed explicit/mixed/implicit H forms.
States use concrete, zero-valence defaults. **Implemented** entries are in the
registry and covered by resolution fixtures.
**PW** marks the adopted graph allocation,
including anionic-ligand donation and the ordinary bridging-carbonyl convention.
D denotes a neutral two-electron donor. No MC state is selected in Ni02–Ni36.

| ID | States / extracted family | Status; new states | Examples and evidence |
| --- | --- | --- | --- |
| Ni01 | Ni/Pt#n4#u2; Pd#n5; Ni#c+2#n3#u2; Ni/Pd/Pt#c+2#n4 | Retained; 7 covered | Atomic ground occupations and free Ni²+ triplet: ASD and HW 1538. Paired d⁸ complex occupations: HW 2033–2035, 2055–2056. The existing paired Pd²+/Pt²+ states remain alongside the free-ion triplets in Ni02; none is replaced. |
| Ni02 | Ni: #c+#n4#u, #c+3#n2#u3; Pd: #c+#n4#u, #c+2#n3#u2, #c+3#n2#u3; Pt: #c+#n4#u, #c+2#n3#u2, #c+3#n2#u3 | Implemented; 8 | Free M+ and M³+: d⁹ ²D₅/₂ and d⁷ ⁴F₉/₂; free Pd²+/Pt²+: d⁸ ³F₄. ASD; distinct from paired complex occupations. |
| Ni03 | Pt: #c-2#n6 | Implemented; 1 | Cs₂Pt: Pt²−, filled 5d¹⁰6s² shell; HW 2042–2043. Ba₂Pt contains the same ion plus electride electrons. |
| Ni04 | Ni: #n5#t2, #n5#t3, #n5#t4; Pd: #n5#t2, #n5#t3, #n5#t4; Pt: #n5#t2, #n5#t3, #n5#t4 | Implemented; 9 | Paired M(0): Ni(PR₃)₂/₃/₄, M(PR₃)₂/₃/₄ (Pd/Pt), Ni(CO)₄ and corresponding PF₃ complexes; HW 2026, 2032, 2038, 2056, 2109–2111. Olefin/N₂ donors reuse the same accepted-pair totals. |
| Ni05 | Pd: #n5#t | Implemented; 1 | Monophosphine Pd(0) in the aryl-halide coupling cycle, HW 2064. One donor; no scaffold-specific state. |
| Ni06 | Ni: #c-2#n4#v4, #c+2#n4#t4, #n4#v2#t2 | Implemented; 3 | Paired four-coordinate Ni(II): Ni(CN)₄²−, NiR₂D₂, NiX₂D₂, Ni(dmg)₂; HW 2029–2035. Ordinary anionic-ligand bonds or donation. The halide-bridged [R₂Ni(μ-X)₂NiR₂]²− graph reuses #n4#v2#t2 (2036). |
| Ni07 | Ni: #c-3#n4#v5, #c+2#n4#t5, #n4#v2#t3, #c+#n4#v#t4 | Implemented; 4 | Paired five-coordinate Ni(II): Ni(CN)₅³−, Ni(CN)₂(PR₃)₃, NiBr(As₄)+ and NiBr(NP₃)+; HW 2026, 2029, 2033. No paired six-coordinate family inferred. |
| Ni08 | Ni: #c-2#n3#u2#v4, #c+2#n3#u2#t4, #n3#u2#v2#t2 | Implemented; 3 | Triplet tetrahedral Ni(II): NiCl₄²− and NiX₂D₂; HW 2033–2035 explicitly distinguishes the paramagnetic and diamagnetic forms. |
| Ni09 | Ni: #c+2#n3#u2#t5, #c+#n3#u2#v#t4 | Implemented; 2 | Triplet five-coordinate Ni(II): Ni(OAsMe₃)₅²+ and NiBr{N(CH₂CH₂NMe₂)₃}+; HW 2033. |
| Ni10 | Ni: #c+2#n3#u2#t6, #n3#u2#v2#t4, #c-4#n3#u2#v6 | Implemented; 3 | Triplet octahedral Ni(II): aqua/ammine, NiCl₂(H₂O)₄, Ni(acac)₂(H₂O)₂, Ni(OH)₆⁴−; HW 2028, 2030–2033. Each graph counts 20 electrons. Nickelocene reuses the first state with two six-electron Cp− donors (2036–2037, 2194–2195). |
| Ni11 | Ni: #n4#u#v#t3, #c+#n4#u#t4, #c-3#n4#u#v4 | Implemented; 3 | Doublet Ni(I): NiX(PR₃)₃ and Ni(CN)₄³−; HW 2026–2027, 2029, 2032. The localized O=Ni=O³− formula in HW 2031 reuses the last state; no extra oxide state inferred. |
| Ni12 | Ni: #n4#u#v#t4, #c+#n4#u#t5 | Implemented; 2 | Five-coordinate NiI(NP₃), NP₃ = N(CH₂CH₂PR₂)₃; HW 2026, 2032 identifies one unpaired electron. Both graphs count 19 electrons. |
| Ni13 | Ni: #n3#u#v3#t2, #c+3#n3#u#t5, #c+3#n3#u#t6, #c+#n3#u#v2#t4, #c-3#n3#u#v6 | Implemented; 5 | Low-spin Ni(III): NiX₃(PR₃)₂, NiCl₂(diars)₂+, NiF₆³− and Ni(phen)₃³+; HW 2027, 2036. Explicit doublet occupation. Six-coordinate graphs count 19. The earlier high-spin solid NiF₃ interpretation is reviewed in Ni37. |
| Ni14 | Ni: #c-2#n3#v6, #c+4#n3#t6 | Implemented; 2 | Paired Ni(IV), NiF₆²−; HW 2036. NiCp*₂²+ reuses the donor state. Solid NiF₄ is not used to assert a molecular four-coordinate occupation. |
| Ni15 | Ni: #c-4#n5#v4, #c-4#n5#v3#h, #c-4#n5#v2#h2, #c-4#n5#v#h3, #c-4#n5#h4 | Implemented; 5 | NiH₄⁴− in Mg₂NiH₄: tetrahedral, paired d¹⁰; HW 2027, 2069. All explicit/mixed/implicit H counts. Ni(CN)₄⁴− independently supports the all-explicit valence; H− donation reuses Ni04. |
| Ni16 | Ni: #n4#v#h#t2, #n4#v#h#t3 | Implemented; 2 | Implicit H in square-planar NiHX(PR₃)₂ and five-coordinate NiH(CN)(PR₃)₃; HW 2027, 2071. Explicit-H and all-donor states reuse Ni06/Ni07. |
| Ni17 | Ni: #c+#n4#v#t3 | Implemented; 1 | Ni₂(CN)₆⁴−: paired Ni–Ni single bond, HW 2029. Cyanide donors give this state; ordinary Ni–CN bonds reuse Ni06. |
| Ni18 | Ni: #c-#n5#v#t2 | Implemented; 1 | NiMe(η²-C₂H₄)₂−, HW 2037: ordinary Ni–Me plus two olefin donations. Methyl-anion donation reuses Ni04. Other Ni olefin, allyl and Cp examples reuse the ordinary/donor families above. |
| Ni19 | Pd: #c+2#n4#t4, #c+#n4#v#t3, #n4#v2#t2, #c-#n4#v3#t, #c-2#n4#v4; Pt: #c+2#n4#t4, #c+#n4#v#t3, #n4#v2#t2, #c-#n4#v3#t, #c-2#n4#v4 | Implemented; 10 | Paired square-planar M(II): MClₖD₄₋ₖ²⁻ᵏ, cyanides, aqua/ammines, cisplatin, MR₂D₂, Zeise salt, acetate and polysulfide complexes; HW 2047–2056, 2060–2064. Full chloro/aqua/ammine substitution series; ordinary anionic-ligand bonds or donors. |
| Ni20 | Pd: #c+2#n4#t5, #c+#n4#v#t4; Pt: #c+2#n4#t5, #c+#n4#v#t4 | Implemented; 4 | Paired five-coordinate M(II): MIL′+ (four As donors), HW 2056. CpMX(PR₃) and CpM(allyl) use the same accepted-pair totals (2061). |
| Ni21 | Pd: #c+2#n4#t6, #n4#v2#t4 | Implemented; 2 | Diamagnetic octahedral PdCl₂(diars)₂, HW 2056. Ordinary chloride or Cl− donation; both count 20 electrons. No six-donor Pt(II) analogue inferred. |
| Ni22 | Pd: #c+4#n3#t6, #n3#v4#t2, #c-2#n3#v6 | Implemented; 3 | Paired Pd(IV): PdX₆²−, Pd(CN)₆²−, PdCl₄(py)₂ and octahedral organopalladium(IV) compounds; HW 2047–2049, 2055, 2061. No unsupported Pd(IV) aqua/ammine substitution series. |
| Ni23 | Pt: #c+4#n3#t6, #c+3#n3#v#t5, #c+2#n3#v2#t4, #c+#n3#v3#t3, #n3#v4#t2, #c-#n3#v5#t, #c-2#n3#v6 | Implemented; 7 | Paired octahedral Pt(IV): full chloro/aqua/ammine series, Pt(CN)₆²−, Pt(OH)₆²−, PtMe₆²−, PtMe₄D₂, Me₃PtD₃+ and [Me₃PtX]₄; HW 2047–2049, 2053–2055, 2061–2062. Ordinary anionic ligands or donation. |
| Ni24 | Pd: #c+#n4#v#t3; Pt: #c+#n4#v#t3 | Implemented; 0 | Paired M(I) dimers with ordinary M–M single bonds: M₂(CNMe)₆²+ and M₂Cl₂(dppm)₂; HW 2056–2057. Ordinary terminal Cl plus M–M reuses Ni19; Cl− donation gives these states. |
| Ni25 | Pt: #c+3#n3#v#t4, #c+3#n3#v#t5, #c-#n3#v5, #c-#n3#v5#t | Implemented; 2 | Pt₂(SO₄)₄L₂²−, with L absent or H₂O: explicit Pt–Pt single bond and diamagnetic dimer; HW 2056–2057, Fig. 31.18a. Four ordinary Pt–O sulfate bonds or four O− donors per Pt; optional axial water. The wider Pt(III) paddlewheel family provides reuse. |
| Ni26 | Pd: #c-#n3#v5, #c+#n3#v3#t2; Pt: #c-#n3#v5, #c+#n3#v3#t2 | Implemented; 3 | M₂Cl₄(μ-CO)₂²−, HW 2057: ordinary M–M and bridging M–C bonds with C=O, following the adopted μ-CO graph convention. Terminal chloride ordinary or donated. No MC type; one Pt state already occurs in Ni25. |
| Ni27 | Pd: #c-2#n5#v2, #c-2#n5#v#h, #c-2#n5#h2 | Implemented; 3 | Linear PdH₂²− in Li₂/Na₂/Cs₂PdH₂, HW 2044, 2069; paired d¹⁰. K₃PdH₃ contains this ion plus separate H−. All three H-count splits; donor form reuses Ni04. |
| Ni28 | Pt: #c-2#n4#v4, #c-2#n4#v3#h, #c-2#n4#v2#h2, #c-2#n4#v#h3, #c-2#n4#h4 | Implemented; 4 | Square-planar PtH₄²−, HW 2044, 2069. K₃PtH₅ contains this ion plus separate H−, not PtH₅³−. All five H-count splits; explicit/donor forms reuse Ni19. |
| Ni29 | Pt: #c-2#n3#v6, #c-2#n3#v5#h, #c-2#n3#v4#h2, #c-2#n3#v3#h3, #c-2#n3#v2#h4, #c-2#n3#v#h5, #c-2#n3#h6 | Implemented; 6 | Octahedral PtH₆²−, HW 2044, 2069. All seven H-count splits; explicit/donor forms reuse Ni23. |
| Ni30 | Pt: #n4#v2#t2, #n4#v#h#t2, #n4#h2#t2 | Implemented; 2 | cis/trans-PtH₂(PR₃)₂ and PtHX(PR₃)₂: paired planar Pt(II), HW 2044, 2048. All H-count splits; terminal hydride donors reuse Ni19. |
| Ni31 | Pd: #n4#v#h#t2 | Implemented; 1 | PdHCl(PEt₃)₂ and the PdHX(PR₃)₂ catalytic intermediate, HW 2048, 2063. Implicit H; explicit and all-donor descriptions reuse Ni19. No unsupported PdH₂(PR₃)₂ analogue. |
| Ni32 | Pt: #n3#v3#h#t2, #n3#v2#h2#t2 | Implemented; 2 | Classical PtHCl₃(PR₃)₂ and PtH₂Cl₂(PR₃)₂, HW 2048. Mixed/implicit H; explicit and donor alternatives reuse Ni23. |
| Ni33 | Pd: #n4#v2#t, #c+2#n4#t3 | Implemented; 2 | Three-coordinate 14-electron Pd(II): Wacker intermediate, Fig. 31.24q, and ArPdX(PtBu₃) in the coupling cycle; HW 2063–2064. Two ordinary anionic ligands plus one donor, or all donors. |
| Ni34 | Pt: #c+#n3#v3#t2, #c+4#n3#t5 | Implemented; 1 | Me₃PtD₂+ formulas in HW 2061; paired Pt(IV) occupation from HW 2055. P05 explicitly identifies five-coordinate Me₃Pt(diphosphine)+ intermediates in reductive elimination. Ordinary methyl bonds or methyl-anion donors. |
| Ni35 | Pt: #n4#u#v, #c+#n4#u#t, #n3#u2#v2, #c+2#n3#u2#t2, #n3#u#v3, #c+3#n3#u#t3, #n2#u2#v4, #c+4#n2#u2#t4, #n3#v4, #c+4#n3#t4, #n2#u#v5, #c+5#n2#u#t5, #n2#v6, #c+6#n2#t6 | Implemented; 14 | Molecular PtF₁–₆: S21, Table 2 and pp. 2–5, gives doublet/triplet/doublet/(triplet and singlet)/doublet/singlet respectively. Both PtF₄ states are spectroscopically supported. Ordinary Pt–F or F− donation; no MC states. |
| Ni36 | Ni: #n4#u#v, #c+#n4#u#t, #n3#u2#v2, #c+2#n3#u2#t2, #n2#u3#v3, #c+3#n2#u3#t3, #n2#u2#v4, #c+4#n2#u2#t4 | Implemented; 8 | Molecular NiF₁–₄: L21F, Scheme 1/Table 1, gives doublet, triplet, quartet and triplet ground states. Ordinary Ni–F or F− donation, with formal metal occupation (PW); the paper reports appreciable ligand spin density, not an atom-local spin measurement. |
| Ni37 | R-NiF₃: reuse Ni10/Ni14 | Resolved; 0 | S99, pp. 561–563, neutron diffraction at 2 K supports predominantly NiII/NiIV, including low-spin d⁶ Ni(IV). No extra Ni#c+3#n2#u3#t6 from Z95’s earlier high-spin interpretation. The H/P polymorph occupations remain unselected; the R result is not generalized to them. |
| Ni38a | Pd: #n4#u#v, #c+#n4#u#t, #n3#u2#v2, #c+2#n3#u2#t2, #n2#u2#v4, #c+4#n2#u2#t4 | Implemented; 6 | Molecular PdF, PdF₂ and PdF₄: W16, Table 1, pp. 1110–1111, 1121–1122. Doublet, triplet and scalar-parent triplet; the SOC A₁ label of PdF₄ is not a singlet assignment. Ordinary Pd–F or F− donation (PW). PdF₃ identification is ambiguous; no PdF₅/₆ selected. |
| Ni38b | Pd: #c+2#n3#u2#t6 | Implemented; 1 | Octahedral triplet Pd(II) in PdF₂ and mixed-valent PdII[PdIVF₆]: HW 2045 explicitly gives two unpaired electrons for the latter. Bridging F− donors use existing #d2/#d3 states; paired Pd(IV) reuses Ni22. Counts 20 electrons. |
| Ni38c | NiF₃+; solid NiF₄; tentative molecular PdF₃/PdF₅/PdF₆ | Not selected for additional types | Z95 does not settle solid NiF₄ coordination/spin or NiF₃+ coordination. W16 leaves the PdF₃/PdF₆ band assignment ambiguous and does not establish PdF₅. No speculative occupation series. Molecular NiF₄ is independently covered in Ni36. |
| Ni39a | Ni: #n4#u#h; Pd: #n4#u#h; Pt: #n4#u#h | Implemented; 3 | Implicit H in doublet NiH/PdH/PtH; L97Ni p. 2056, A01Pd Table 3, A01Pt pp. 1561–1562. Explicit-H and H−-donor states reuse Ni35/Ni36/Ni38a. |
| Ni39b | Ni: #n4#v2, #n4#v#h, #n4#h2, #c+2#n4#t2 | Implemented; 4 | Bent singlet NiH₂ observed in L97Ni, pp. 2057–2059; the paper leaves computed singlet/triplet ground-state ordering unresolved. Select the observed bent singlet, not extra triplet H-count variants. Three H-count splits and H− donors. |
| Ni39c | Pd: #c-#n4#u#v2, #c-#n4#u#v#h, #c-#n4#u#h2, #c+#n4#u#t2 | Implemented; 4 | Linear doublet PdH₂−: A01Pd Table 8 and pp. 3060–3062; ELF explicitly identifies two ordinary Pd–H bonds. Three H-count splits and H− donors. The bare neutral PdH₂/PdH₃ molecules were not observed. |
| Ni39d | Pt: #n4#v2, #n4#v#h, #n4#h2, #c+2#n4#t2 | Implemented; 4 | Bent singlet PtH₂: A01Pt Table IV, pp. 1562–1565. Three H-count splits and H− donors. |
| Ni39e | Pt: #n3#u#v3, #n3#u#v2#h, #n3#u#v#h2, #n3#u#h3, #c+3#n3#u#t3 | Implemented; 3 | Observed Cₛ PtH₃ doublet: A01Pt Fig. 6/Table V, pp. 1562–1565. Three ordinary Pt–H bonds, with all H-count splits. Explicit-H and donor forms reuse Ni35. The calculated (H₂)PtH isomer and higher intact-H₂ adducts remain deferred. |
| Ni39f | Pd: #c+#n4#u#m0 | Implemented; 1 | Rhombic singlet (PdH)₂: A01Pd pp. 3061–3062/Table 9 explicitly finds two Pd–H–Pd 3c2e basins, no Pd–Pd basin. PW graph: two H− bridges with contributions [0,2,0], two Pd(I) d⁹ centers whose unpaired electrons couple to a singlet; no Pd–Pd edge. Existing H#c-#m2. |
| Ni39g | Intact-H₂ complexes, including Pd(H₂)₁–₃, (H₂)PdH and PtH₂(H₂) | Deferred; 0 | A01Pd/A01Pt distinguish these from dissociated hydrides. Retain the shared η²-H₂ donor deferral; no ordinary-H-count or MC substitutes. Interstitial PdHₓ is not a molecular MHₙ graph. |
| Ni40a | Ni: #n3#v4, #c+4#n3#t4 | Implemented; 2 | Linear singlet O=Ni=O: C97Ni Table 2, pp. 3113–3114. Ordinary double bonds or two-pair O²− donors (PW); existing O#c-2#n2#d2. NiO triplet reuses Ni36; PtO triplet and OPtO singlet reuse Ni35 (B99 Tables 2–3). |
| Ni40b | Ni: #n4#v2, #c+2#n4#t2; Pd: #n4#v2, #c+2#n4#t2; Pt: #n4#v2, #c+2#n4#t2 | Implemented; 2 | Cyclic singlet M(O₂): C97Ni Table 2/pp. 3114–3115 and B99 Tables 1–2/pp. 5459–5460. PW closed-shell peroxide graphs: O–O and two M–O single bonds, or two O− donors. Ni/Pt reuse Ni39b/Ni39d; Pd adds two states. Source supplies the cyclic structure and singlet; the integer bond/charge allocation is our convention. |
| Ni40c | Remaining O₂ adducts/oxide clusters, ONiO−, trigonal PtO₃ and calculated OPdO/NiO₃ | Not selected for additional types | C97Ni/B99 distinguish the isomers. Do not assign all molecular spin to the metal in superoxo/O₂ complexes. ONiO− and trigonal PtO₃ identifications are tentative (C97Ni p. 3115; B99 p. 5461); OPdO/NiO₃ are unobserved calculated structures. No further default types selected from these cases. |
| Ni41a | Pt: #c-#n2#u#v6, #c+5#n2#u#t6 | Implemented; 2 | PtF₆−: BL62 pp. 5260–5261 explicitly assigns nonbonding dπ⁵. Ordinary Pt–F or F− donation. The O₂+ counterion uses existing types. |
| Ni41b | PdF₆−; candidate Pd#c-#n2#u#v6 / Pd#c+5#n2#u#t6 | Omitted; 0 | F76 confirms the octahedral anion and reports susceptibility close to O₂PtF₆, but does not explicitly assign the d⁵ occupation. The displayed doublet is an analogy with PtF₆−, not a directly established local assignment; omitted from this selection. |
| Ni42 | Pt: #c-2#n4#v4#m0, #c-2#n4#v3#h#m0, #c-2#n4#v2#h2#m0, #c-2#n4#v#h3#m0, #c-2#n4#h4#m0, #c+2#n4#t4#m0 | Implemented; 6 | [Pt₂H₉]⁵−: B95 pp. 899–900 establishes four terminal hydrides per Pt, one linear Pt–H–Pt bridge and temperature-independent diamagnetism. PW bridge [0,2,0] on H−; no Pt–Pt bond (3.704 Å). Five terminal-H splits plus all-terminal-H− donors; existing H#c-#m2. |
| Ni43 | NiN₂⁴− in Sr₂NiN₂ | Not selected for additional types | K00, pp. 547–550, gives isolated linear units, metallic/Pauli behavior and a band description. This does not settle a unique atom-local Ni occupation; no new nitride state inferred from formal Ni(II). |
| Ni44a | Ni: #n3#u#v3, #c+3#n3#u#t3; Pd: #n3#u#v3, #c+3#n3#u#t3 | Implemented; 4 | OMF₂ triplets: L21O Tables 1/3, SI pp. 7–10, 33 and Fig. S2.17 explicitly support oxyl character. PW integer localization: O#n2#u#v with an ordinary M–O single bond, or O•− donation. One metal unpaired electron plus one on O. Pt reuses Ni35; quartet OMF/OPtF₃ and triplet FOPdF also reuse fluoride states. No MC. |
| Ni44b | O: #c-#n2#u#d | Implemented; 1 | Terminal O•− donor companion for the Ni44a all-donor graphs and OMF/OPtF₃. L21O explicitly assigns oxyl radical character; the integer donor allocation is PW. Ordinary M–O• graphs already use O#n2#u#v. |
| Ni45a | Ni: #n5#t3#m0 | Implemented; 1 | HNi₂(CO)₆−: L79 pp. C42–C43 gives the Ni₂(CO)₆ framework and infers bridging H from steric energy maps; P82 pp. 93–94, 98 treats the closed-shell cluster; I00 p. 1989 explicitly identifies a three-center Ni–H–Ni bond. PW graph: H− bridge [0,2,0], paired Ni(0) centers, three CO donors each and no separate Ni–Ni edge. Existing H#c-#m2; no additional capacity. |
| Ni45b | Ni₂(CO)₆²− | Not selected for additional types | L79 p. C43 gives indirect solution evidence, not an isolated or structurally characterized dianion. P82 pp. 93–94 explicitly uses an assumed geometry and distinguishes HNi₂(CO)₆− as the characterized dinuclear species. No new state from the uncertain dianion or calculated dihydride. |
| Ni46 | Pt₂ units in Ca₃Pt₂/Ba₃Pt₂ | Not selected for additional types | K04 pp. 1391–1393 supports Pt–Pt covalency within a metallic band structure. Its Bader charges (−1.2/−1.5 per Pt) are partial charges, not formal atomic charges; neither a discrete Pt₂³− ion nor a local pair/unpaired occupation is established. No additional Pt type or MC substitute selected. |
| Ni47 | PdCl₅³−, PdBr₅³−, PdBr₆⁴− | No additional default types selected | HW 2047 lists these formulas without settling coordination versus longer contacts. Paired five-/six-donor types are independently supported in Ni20/Ni21; these examples are not needed for their inclusion or coverage. |
| Ni48 | Higher carbonyl/thiolate clusters, triple-decker Ni₂Cp₃+, Pt₁₂O₈(SO₄)₁₂⁴−, Pt₄ acetate and ligand-specific hydride/borohydride assemblies | Not selected for new default types | HW 2029, 2036–2038, 2052, 2057, 2070–2071, 2131. No per-scaffold MC types. The selected M–M dimers and ordinary/donor families remain in Ni17/Ni24–Ni26. Weak stacking contacts and extended metallic structures do not add resolution bonds. |
| Ni49 | [PtI(dmpe)₂(κ-I₂)]+ | Not selected for new default types | M01, p. 1642, establishes six-coordinate Pt(II) and explicitly describes Pt d(z²) → I₂ σ* donation. Four phosphine donors and one I− donor give 18 electrons in that graph; the outgoing pair is already part of Pt’s eight electrons. This specialized complex supplies capacity evidence, not additional registry states. |

#### Nickel-group limits and coverage

| Item | Finding / next action |
| --- | --- |
| Capacity | Ni/Pd capacity raised from 18 to 20. Octahedral Ni(II) in Ni10 and Pd(II) in Ni21/Ni38b require 20; Ni12/Ni13 also require 19. Pt retains 18. The six-coordinate M01 complex has five incoming donor pairs and one Pt→I₂ donation, so its graph also counts 18. This changes the graph-count limit, not free-ion occupations. |
| ValenceTable | Unchanged. Selected graphs have fixed H counts and no metal aromatic contribution; Counts resolves them without Ni/Pd/Pt table entries. |
| State audit | 173 distinct additions: 172 metal states and one O companion. Each satisfies E − charge − valence − implicit H − 2 × donated pairs − MC contribution = 2 × lone pairs + unpaired electrons (E = 10 for Ni/Pd/Pt, 6 for O). Capacity audit uses E − charge + valence + implicit H + 2 × accepted pairs. |
| MC review | All eight MC states implemented and covered. Six Pt₂H₉⁵− states: five explicit/mixed/implicit terminal-H variants of Pt#c-2#n4#v4#m0, plus Pt#c+2#n4#t4#m0. Ni#n5#t3#m0 for HNi₂(CO)₆− and Pd#c+#n4#u#m0 for two separate 3c2e bridges in (PdH)₂ complete the set. Each bridge uses existing H#c-#m2 with [0,2,0]. No whole-cluster MC bond. The integer charge/contribution allocations are PW graph choices. |
| Other graph selections | Ni40b uses the closed-shell peroxide allocation for cyclic singlet M(O₂). Ni44a/b uses the atom-local oxyl allocation with accepted O#c-#n2#u#d for its all-donor alternative. Fractional spin populations in L21O are evidence for oxyl character, not literal integer graph attributes. |
| Conformance | 190 new fixtures: 12 atomic ions, 39 hydrides, 22 small inorganic molecules, 20 inorganic ions and 97 donor–acceptor graphs. All 173 added states are inferred with AtomTyping and Counts under Strict and MostSaturated; all 760 resolved graphs pass valence invariants. Both Ni(II) and PtF₄ occupations have distinguishing inputs. All seven retained metal states are covered. |
| Existing cases | Inputs preserved. Bare Pd²⁺/Pt²⁺ snapshots now leave the paired/triplet choice unresolved under Strict AtomTyping; MostSaturated retains the paired result. Separate supplied-occupation fixtures cover both states. Other existing snapshots are unchanged. |
| Verification | 3,092 resolution cases, 223 graph property tests, 3,594 graph unit tests and 312 chemistry tests pass; formatting and Clippy pass. Registry audit confirms exactly the accepted 173 additions, no removed states, correct valence buckets and charge ordering. |
| Outstanding scope | Selected nickel-group scope complete; no implementation or source work remains. PdF₆− is omitted; no further occupation search is pending. Other tentative or locally unresolved cases are explicitly unselected. η²-H₂ and specialized aggregates remain deferred/unselected. Ds belongs to the later superheavy-element chapter. |

#### Nickel-group sources and access

HW references above are printed Band II pages. Retrieved sources are archived under
`materials/valence/`; extraction text and candidate audits are temporary files under
`scratch/hw-nickel-review/`. All sixteen requested papers and the oxyfluoride SI
are present, reviewed and archived; the duplicate Kowach download is one source.
No required source check remains outstanding.

| Key | Source | Contribution / access |
| --- | --- | --- |
| ASD | NIST Atomic Spectra Database, ground configurations/terms. Archived query outputs: [Ni](../materials/valence/NIST-ASD-2026-ni-ground-states.html), [Pd](../materials/valence/NIST-ASD-2026-pd-ground-states.html), [Pt](../materials/valence/NIST-ASD-2026-pt-ground-states.html); [live database](https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html) | Reviewed free-atom/ion occupations, Ni01/Ni02. Kept separate from complex spin states. |
| S21 | Senges et al., “Matrix Isolation Spectroscopic and Relativistic Quantum Chemical Study of Molecular Platinum Fluorides PtFₙ (n = 1–6) Reveals Magnetic Bistability of PtF₄,” *Chem. Eur. J.* **27** (2021), 13642–13650. [DOI](https://doi.org/10.1002/chem.202102055), [PDF](../materials/valence/Senges-2021-platinum-fluorides.pdf) | Main paper reviewed and archived. Table 2 and manuscript pp. 2–5 support Ni35, including both observed PtF₄ occupations and closed-shell PtF₆. SI download failed; it is not needed for these main-text assignments. |
| L21F | Li et al., “Searching for Monomeric Nickel Tetrafluoride: Unravelling Infrared Matrix Isolation Spectra of Higher Nickel Fluorides,” *Angew. Chem. Int. Ed.* **60** (2021), 6391–6394. [DOI](https://doi.org/10.1002/anie.202015501), [repository copy](https://publikationen.bibliothek.kit.edu/1000130872/107155258), [PDF](../materials/valence/Li-2021-nickel-fluorides.pdf) | Reviewed and archived. Scheme 1/Table 1 explicitly give NiF₁–₄ ground terms; Ni36. The ligand-spin discussion qualifies the formal metal allocation. Its introduction identifies the later bulk R-NiF₃ reassignment, Ni37. |
| Z95 | Žemva et al., “Thermodynamically Unstable Fluorides of Nickel: NiF₄ and NiF₃ Syntheses and Some Properties,” *J. Am. Chem. Soc.* **117** (1995), 10025–10034. [DOI](https://doi.org/10.1021/ja00145a013), [PDF](../materials/valence/Zemva-1995-nickel-fluorides.pdf) | Retrieved and relevant structural/magnetic sections reviewed. Distinguishes the three solid NiF₃ forms; high-spin interpretation at pp. 10030–10031, with mixed-valence contributions not excluded. Solid NiF₄ and NiF₃+ remain incompletely characterized here. |
| L21T | Li, *Matrix Isolation and IR Spectroscopic Study on Late-Transition Fluorides, Metal Oxo Complexes and Actinide Halides*, FU Berlin dissertation (2021). [Repository PDF](https://d-nb.info/1232726710/34), [PDF](../materials/valence/Li-2021-molecular-fluorides-thesis.pdf) | Overview reviewed and archived. The public thesis removes the published articles and their SI; L21O and its SI have now been reviewed separately. |
| P05 | Procelewska et al., “Mechanistic information on the reductive elimination from cationic trimethylplatinum(IV) complexes to form carbon-carbon bonds,” *Inorg. Chem.* **44** (2005), 7732–7742. [DOI](https://doi.org/10.1021/ic050478%2B), [abstract](https://pubmed.ncbi.nlm.nih.gov/16241122/), [archived PubMed record](../materials/valence/Procelewska-2005-five-coordinate-abstract.xml) | Abstract reviewed and archived: explicitly distinguishes the five-coordinate cationic intermediate from a six-coordinate solvent adduct. Supports coordination in Ni34; its paired occupation comes from HW, not an inferred magnetic measurement in this paper. |
| S99 | Shen et al., *C. R. Acad. Sci. Paris, Série IIc* **2** (1999), 557–563, R-NiF₃ structure and mixed nickel fluorides. [DOI](https://doi.org/10.1016/S1387-1609%2800%2988566-7), [PDF](../materials/valence/Shen-1999-nickel-fluoride.pdf) | Reviewed and archived. Neutron diffraction supports the R-NiF₃ mixed-valent description, Ni37; no new high-spin Ni(III) type. |
| L97Ni | Li et al., magneto-infrared study of matrix-isolated NiH/NiH₂, *J. Chem. Phys.* **106** (1997), 2055. [DOI](https://doi.org/10.1063/1.473342), [PDF](../materials/valence/Li-1997-nickel-hydrides.pdf) | Reviewed and archived. Doublet NiH and observed bent singlet NiH₂; calculated singlet/triplet ordering remains unresolved in this paper. Ni39a/b. |
| A01Pd | Andrews, Wang, Alikhani and Manceron, “Observed and Calculated Infrared Spectra of Pd(H₂)₁,₂,₃ Complexes and Palladium Hydrides in Solid Argon and Neon,” *J. Phys. Chem. A* **105** (2001), 3052. [DOI](https://doi.org/10.1021/jp003721t), [PDF](../materials/valence/Andrews-2001-palladium-hydrides.pdf) | Reviewed and archived. PdH/PdH₂−, intact-H₂ complexes, and the two explicitly assigned 3c2e bridges of (PdH)₂ (pp. 3061–3062/Table 9). Ni39a/c/f/g. |
| A01Pt | Andrews, Wang and Manceron, “Infrared spectra and density functional calculations of platinum hydrides,” *J. Chem. Phys.* **114** (2001), 1559–1566. [DOI](https://doi.org/10.1063/1.1333020), [PDF](../materials/valence/Andrews-2001-platinum-hydrides.pdf) | Reviewed and archived. PtH doublet, PtH₂ singlet and observed PtH₃ doublet; distinguishes the latter from the calculated (H₂)PtH isomer. Ni39a/d/e/g. |
| C97Ni | Citra et al., nickel/O₂ matrix-isolation study, *J. Phys. Chem. A* **101** (1997), 3109–3118. [DOI](https://doi.org/10.1021/jp970019x), [PDF](../materials/valence/Citra-1997-nickel-oxides.pdf) | Reviewed and archived. NiO/ONiO and cyclic peroxide structures/states; separates the other O₂/oxide isomers and tentative ONiO−. Ni40. |
| B99 | Bare et al., platinum/palladium–O₂ matrix-isolation study, *J. Phys. Chem. A* **103** (1999), 5456–5462. [DOI](https://doi.org/10.1021/jp990823e), [PDF](../materials/valence/Bare-1999-platinum-palladium-oxides.pdf) | Reviewed and archived. PtO/OPtO and cyclic Pd/Pt(O₂) states; PtO₃ identification remains tentative and OPdO was not observed. Ni40. |
| W16 | Wilson et al., “A Matrix Isolation and Computational Study of Molecular Palladium Fluorides: Does PdF₆ Exist?” *Inorg. Chem.* **55** (2016), 1108–1123. [DOI](https://doi.org/10.1021/acs.inorgchem.5b02273), [PDF](../materials/valence/Wilson-2016-palladium-fluorides.pdf) | Reviewed and archived. PdF/PdF₂/PdF₄ assignments, scalar versus SOC labels, and limits of the higher-fluoride identifications. Ni38a/c. |
| BL62 | Bartlett and Lohmann, “Fluorides of the noble metals. Part II. Dioxygenyl hexafluoroplatinate(V),” *J. Chem. Soc.* (1962), 5253. [DOI](https://doi.org/10.1039/JR9620005253), [PDF](../materials/valence/Bartlett-Lohmann-1962-platinum-fluoride.pdf) | Reviewed and archived. Explicit Pt(V) nonbonding dπ⁵ assignment, pp. 5260–5261. Ni41a. |
| F76 | Falconer et al., “Dioxygenyl hexafluoropalladate(V),” *J. Inorg. Nucl. Chem.*, supplement (1976), 59–60. [DOI](https://doi.org/10.1016/0022-1902%2876%2980595-7), [PDF](../materials/valence/Falconer-1976-palladium-fluoride.pdf) | Reviewed and archived. Octahedral PdF₆− and bulk susceptibility; no explicit d⁵ local occupation. Ni41b omitted because the local occupation remains unresolved. |
| B95 | Bronger and à Brassard, “Li₅Pt₂H₉, a Complex Hydride Containing Isolated [Pt₂H₉]⁵− Ions,” *Angew. Chem. Int. Ed. Engl.* **34** (1995), 898–900. [DOI](https://doi.org/10.1002/anie.199508981), [PDF](../materials/valence/Bronger-1995-lithium-platinum-hydride.pdf) | Reviewed and archived. Terminal/bridging hydrides and diamagnetism of Li₅Pt₂H₉. Ni42; integer bridge contributions are PW. |
| K00 | Kowach et al., “An Unusual Metallic Nitride: Sr₂NiN₂,” *J. Solid State Chem.* **154** (2000), 542–550. [DOI](https://doi.org/10.1006/jssc.2000.8877), [PDF](../materials/valence/Kowach-2000-strontium-nickel-nitride.pdf) | Reviewed and archived. Isolated NiN₂ units and metallic band description; no unique local occupation selected. Ni43. |
| L21O | Li et al., “Molecular oxofluorides OMFₙ of nickel, palladium and platinum: oxyl radicals with moderate ligand field inversion,” *Inorg. Chem. Front.* **8** (2021), 1215–1228. [DOI](https://doi.org/10.1039/D0QI01151G), [PDF](../materials/valence/Li-2021-oxofluorides.pdf), [SI](../materials/valence/Li-2021-oxofluorides-SI.pdf) | Main paper and SI reviewed and archived. Explicit oxyl character, OMF/OMF₂/OPtF₃ multiplicities, and FOPdF d⁸ assignment. Ni44a/b. |
| I00 | Ignatyev et al., “Binuclear Homoleptic Nickel Carbonyls,” *J. Am. Chem. Soc.* **122** (2000), 1989–1994. [DOI](https://doi.org/10.1021/ja9914083), [PDF](../materials/valence/Ignatyev-2000-binuclear-nickel-carbonyls.pdf) | Retrieved and archived; p. 1989 identifies the Ni–H–Ni bridge in HNi₂(CO)₆− and its primary sources. Its calculated neutral Ni₂(CO)₅–₇ structures are not added as characterized compounds. |
| L79 | Longoni, Manassero and Sansoni, “Synthesis and crystal structure of [(Ph₃P)₂N][Ni₂(CO)₆H],” *J. Organomet. Chem.* **174** (1979), C41–C44. [DOI](https://doi.org/10.1016/S0022-328X(00)91513-4), [PDF](../materials/valence/Longoni-1979-dinickel-carbonyl-hydride.pdf) | Reviewed and archived. Ni₂(CO)₆ framework, bridging-H inference from steric energy maps, and indirect evidence for the unisolated dianion. Ni45a/b. |
| P82 | Pacchioni, Fantucci and Valenti, “The electronic structure of carbonyl metal clusters: I. An INDO investigation of mono-, di- and trinuclear nickel complexes,” *J. Organomet. Chem.* **224** (1982), 89–105. [DOI](https://doi.org/10.1016/S0022-328X(00)82570-X), [PDF](../materials/valence/Pacchioni-1982-nickel-carbonyl-clusters.pdf) | Reviewed and archived. Closed-shell HNi₂(CO)₆− calculation and protonation analysis; explicitly separates characterized and assumed cluster structures. Ni45a/b. H− contribution allocation is PW, not the paper’s formal Ni(−I)/H+ bookkeeping. |
| K04 | Karpov, Wedig and Jansen, “Pt-Dumbbells in Ba₃Pt₂: Interplay of Geometric and Relativistic Effects on Pt–Pt Bonding,” *Z. Naturforsch. B* **59** (2004), 1387–1394. [DOI](https://doi.org/10.1515/znb-2004-11-1204), [PDF](../materials/valence/Karpov-2004-platinum-dumbbells.pdf) | Reviewed and archived. Pt–Pt covalency, metallic bands and partial-charge analysis; no definite local occupation for a discrete Pt₂³− graph. Ni46 remains unselected. |
| M01 | Makiura et al., “An unusual six-co-ordinate platinum(II) complex containing a neutral I₂ ligand,” *Chem. Commun.* (2001), 1642–1643. [DOI](https://doi.org/10.1039/B103648N), [PDF](../materials/valence/Makiura-2001-platinum-diiodine.pdf) | Reviewed and archived. X-ray structure and XPS establish six-coordinate Pt(II); p. 1642 explicitly describes Pt→I₂ donation. Supports retaining capacity 18 for this graph, Ni49; no additional default state selected. |

### Transition-metal classes (chapter XXXII)

Selected scope implemented 2026-11-04. HW II, pp. 2067–2217 covers inorganic complexes,
carbonyls and related ligands, and organometallic families; pp. 2217–2223 contains
the bibliography. Much of the inorganic/carbonyl evidence already supports the
completed group tables. The cross-group selection and source checks are complete.

**21 additions implemented** (Cr 4, Mo 2, W 4;
Sc/Ti/Fe/Y/Zr/Nb/Ru/La/Hf/Ta/Os 1 each), covered by 21 independent fixtures.
Sc/Y/La capacity is raised from 18 to 20. ValenceTable is unchanged; no new MC or
C/N/O companion state is added. The NMR-supported Ti/Zr/Hf paired occupations
are accepted assignments.

#### Chapter-XXXII work table

Element prefixes apply to each suffix. States use concrete, zero-valence defaults.
The selected graphs retain the existing Cp− six-electron donor and terminal-CO
donor conventions. Bridging CO uses ordinary M–C bonds and C=O. These integer
graph allocations are PW; HW supplies the structures and bonding descriptions.

| ID | Family / states | Disposition | Evidence and scope |
| --- | --- | --- | --- |
| Tm01 | Hydrides, halides, cyanides, oxides, nitrides and N₂/O₂ complexes | Reuse / evidence | HW 2067–2108 reinforces the group selections. Ordinary ligand bonds, anionic donation and the selected hydride bridges cover the reviewed common motifs. Extended solids and higher nuclearity cages do not supply molecular types merely from stoichiometry. |
| Tm02 | Carbonyls; CS/CSe/CTe, CNR and PF₃ substitutions | Mainly existing states | HW 2108–2150. The reviewed substitutions reuse existing metal acceptor totals and C/N/P/S/Se/Te ligand states. The Ti addition and unsupported neutral Co(CNR)₅ are separated in Tm13a/b. |
| Tm03 | Cr: #c-4#n5#t4 | Implemented; 1 | Cr(NO)₄ with four NO⁺ donors. HW 2150–2151 explicitly gives the NO⁺ description and diamagnetism. Existing N#v3#d and O#c+#n#v3 suffice. |
| Tm04 | Other selected nitrosyls | Existing states | HW 2151–2157: Co(CO)₃NO, Fe(CO)₂(NO)₂, Mn(CO)(NO)₃, Mn(CO)₄NO, V(CO)₅NO, and the Fe/Co nitrosyl-halide dimers reuse existing types. The NO− description of bent [Co(NH₃)₅NO]²⁺ also does. HW 2154 gives exceptions to assigning NO charge from geometry alone; no universal linear/bent rule is adopted. |
| Tm05 | Cr/W: #n2#v2#t5 | Implemented; 2 | Localized M=C in (CO)₅M=C(R)(OR), HW 2165–2167, Figs. 32.63c/32.64. Mo's corresponding state already exists. The alternative neutral-carbene donation uses existing M#n3#t6 and C#v2#d. |
| Tm06 | W: #n#v4#t4 | Implemented; 1 | X–W≡CR with four CO donors, HW 2168, Fig. 32.65b; ordinary halide and carbyne bonds. Cr/Mo counterparts already exist. The separate Schrock W(OR)₃≡CR example uses existing W#v6. |
| Tm07 | Cr/Mo/W: #c+#n2#v#t6 | Implemented; 3 | [CpM(CO)₃]₂: one ordinary M–M bond, three terminal CO donors per metal. HW 2202–2203, Fig. 32.100f explicitly describes the single bond. |
| Tm08 | Cr/Mo/W: #c+#n#v3#t5 | Implemented; 3 | [CpM(CO)₂]₂: ordinary M≡M and two terminal CO donors per metal. HW 2202 explicitly assigns the triple bond; Fig. 32.100g. |
| Tm09 | Fe/Ru: #c+#n2#v3#t4; Os: #c+#n3#v#t5 | Implemented; 3 | [CpM(CO)₂]₂, HW 2202–2203, Fig. 32.100i. Fe/Ru: two μ-CO bridges, one terminal CO per metal and M–M. Os: no bridging CO, as explicitly noted in the figure; M–M plus two terminal CO per metal. |
| Tm10a | Olefins, dienes, allyls, Cp and arenes | Existing coverage / graph choices | HW 2175–2216 adds examples for existing acceptor totals: Ni(olefin)₃, Ni(COD)₂, Fe(CO)₃(diene), M(CO)₃(arene) (Cr/Mo/W), and the common metallocenes. Hapticity is represented by dative relations; it does not itself require MC states. Mixed hapticity must be read from the particular structure. |
| Tm10b | Zr/Hf: #c+2#n#t8 | Implemented; 2 | Cp₂M(CO)₂, HW 2201–2203, Table 32.24/Fig. 32.100a: two Cp− and two terminal CO donors, with HW's 18-electron count. TB76, pp. 298–300, establishes monomers by MS/IR/NMR; sharp Cp ¹H signals at τ 4.89/4.97 support the paired allocation alongside the existing Ti counterpart (Tg20c). Pairing is inferred; the paper reports no magnetic measurement or explicit electronic multiplicity. |
| Tm11 | Fe: #n3#u2#t6 | Not selected | Fe(η⁶-C₆Me₆)₂: HW 2180 explicitly gives 20 electrons and two unpaired electrons; see also 2207. Fits Fe capacity 21, but excluded from the default additions because the present support is this narrow example. Fe(C₆H₆)₂ has η⁶/η⁴ coordination and reuses Fe#n4#t5. |
| Tm12 | Sc/Y/La: #c+3#t10 | Implemented; 3 | [M(η⁸-C₈H₈)₂]−, HW 2212, 2214: two COT²− ten-electron donors to M³⁺. Capacity raised from 18 to 20, without MC or new ligand types. HW 2214 explicitly lists these three true sandwich structures; Ac's formula on 2212 alone is not used to extend the assignment. |
| Tm13a | Ti: #n2#t7 | Implemented; 1 | Ti(CO)₅(dmpe), Ti(CO)₃(dmpe)₂ and Ti(CO)₄(trmpe), HW 2109; CE88, pp. 164–165, 168–170. Seven neutral pair donors and two remaining pairs give 18 electrons. Resolved ¹H couplings and multinuclear NMR support the paired allocation, which is inferred rather than an explicit spin assignment; no magnetic measurement is reported. The state covers all three isolated complexes, including the more stable tri-/tetracarbonyls. |
| Tm13b | Neutral Co(CNR)₅ | Omitted; occupation unsupported | HW 2150 lists the formula without specifying R or an occupation. No primary evidence establishing that species and its occupation was obtained; no state is added prospectively. |
| Tm13c | Nb/Ta: #n2#u#t6 | Implemented; 2 | HW 2207; CDG83, pp. 1152–1155, 1157: ESR explicitly supports the ²A₁, e₂⁴a₁¹ occupation for Nb(arene)₂ and Ta(C₆H₆)₂; Nb photoelectron spectra corroborate it. Nb examples include benzene, toluene and mesitylene. Two neutral six-electron arene donors give 17 electrons under the existing graph convention. |
| Tm14 | Bare [CpM(NO)₂]⁺; Fe(NO)₄/Co(NO)₃; isolated carbene formulas | Not selected | LW84's opening page says bare [CpCr(NO)₂]⁺ could not be isolated; its MeCN adduct reuses Cr#n3#t6. No bare Mo/W state is inferred by analogy. HW 2150 questions the homoleptic Fe/Co nitrosyl identifications. Cr(CO)₄(carbene) and Ni(CO)₃(carbene)₂ on 2165 require primary verification before deriving extra acceptor/occupation states. |
| Tm15 | η²-H₂, agostic C–H and other σ-bond donors | Deferred model work | HW 2071–2074, 2169–2174 explicitly describes donation from an existing bond. Counting that pair as an independent atomic lone-pair donation would double-count it. No atom-type or MC workaround is added. |
| Tm16 | Higher clusters, multideckers and specialized redox endpoints | Not selected | HW 2081–2084, 2113–2121, 2197–2203, 2212–2216. No per-scaffold MC types or automatic state ladders from the listed formulas. The small dimers in Tm07–Tm09 retain explicit ordinary metal–metal bonds. |

The capacity expansion admits two ten-electron COT²− donors to Sc/Y/La³⁺ in Tm12:
20 electrons, with no change to their three valence electrons or to Ac's capacity.
The other states count 18 electrons, except Nb/Ta bisarenes at 17.

All 21 added states have inferred occupation coverage in 21 donor–acceptor fixtures;
both metal sites are checked in each dimer. AtomTyping and Counts agree under
Strict and MostSaturated for every new input, with passing electron-count validation
(120 atom-state checks across the four configurations). All 3,092 earlier snapshots
remain unchanged. The full 3,113-case conformance suite, 223 graph property tests,
3,594 graph unit tests and 312 chemistry tests pass, as do formatting and Clippy.
No selected work remains outstanding for this chapter.

#### Chapter-XXXII sources and access

HW citations are printed Band II pages. Intermediates are in
`scratch/hw-chapter32-review/`; retained source material is under `materials/valence/`.

| Key | Source | Access / contribution |
| --- | --- | --- |
| CE88 | Chi, Frerichs, Stein, Blackburn and Ellis, “Highly Reduced Organometallics. 20. The First General Synthetic Route to Zerovalent Titanium Carbonyls. Synthesis and Properties of Ti(CO)₅(Me₂PCH₂CH₂PMe₂),” *J. Am. Chem. Soc.* **110** (1988), 163–171. [DOI](https://doi.org/10.1021/ja00209a026), [PDF](../materials/valence/Chi-1988-titanium-carbonyls.pdf) | Supplied full paper reviewed and archived. Seven-coordinate Ti(CO)₅(dmpe) is established by ¹³C/³¹P coupling, pp. 168–169; its detailed geometry remains tentative. Ti(CO)₃(dmpe)₂ and Ti(CO)₄(trmpe) independently share the seven-donor count. NMR supports the adopted paired occupation in Tm13a; no explicit electronic multiplicity or magnetic measurement. |
| CDG83 | Cloke, Dix, Green, Perutz and Seddon, “Electron spin resonance and photoelectron studies of neutral bis(η-arene)metal compounds (metal = Ti, V, Nb, Ta, Mo, and W),” *Organometallics* **2** (1983), 1150–1159. [DOI](https://doi.org/10.1021/om50003a014), [PDF](../materials/valence/Cloke-1983-bisarenes.pdf) | Supplied full paper reviewed and archived. ESR measurements and simulations, pp. 1152–1154, and the explicit ²A₁, e₂⁴a₁¹ assignment on p. 1157 settle Tm13c. Photoelectron spectra independently support Nb's occupation, pp. 1155–1156; Ta's photoelectron spectrum was not obtained because of thermal decomposition. |
| LW84 | Legzdins and Wassink, “Organometallic nitrosyl chemistry. 21. Redox properties of bis[(η⁵-cyclopentadienyl)dinitrosylchromium] and related complexes,” *Organometallics* **3** (1984), 1811–1817. [DOI](https://doi.org/10.1021/om00090a006), [archived first page](../materials/valence/Legzdins-Wassink-1984-nitrosyl-redox-first-page.pdf) | Only p. 1811 retrieved, from a repository PDF appended to the preceding article. Its abstract distinguishes the unisolated bare cation from the MeCN adduct; sufficient for the restricted Tm14 finding, not a full-paper review. |
| TB76 | Thomas and Brown, “Zirconocene dicarbonyl and hafnocene dicarbonyl, first carbonyl derivatives of zirconocene and hafnocene,” *J. Organomet. Chem.* **111** (1976), 297–301. [Repository PDF](https://deepblue.lib.umich.edu/bitstream/handle/2027.42/21771/0000165.pdf?sequence=1), [archived PDF](../materials/valence/Thomas-Brown-1976-metallocene-dicarbonyls.pdf) | Supplied full paper reviewed and archived. MS/IR/NMR and isolation data establish the two monomeric dicarbonyls, pp. 298–300. Sharp Cp ¹H NMR signals support Tm10b's adopted paired occupation; no explicit electronic multiplicity or magnetic measurement. |

### Lanthanides

Reviewed 2026-11-04: HW II, chapter XXXV, pp. 2288–2310; bibliography 2310–2311.
**Implemented: 232 additions**, including 3 La companions and 3 hydride-bridge
MC states; all 34 earlier Ce–Lu bare-atom/ion states are retained. The capacity
changes and Nd/Dy enumeration bounds are implemented. The two Nd crown states are
excluded; Nd capacity remains 26. The 246 new fixtures are now split across
`lanthanides/` (218),
`organometallics/` (25) and `multicenter_bonds/` (3). Actinide fixtures use
`actinides/` and the applicable bonding categories.

#### Occupations used in the table

Each `M` pattern below expands to the listed elements with the corresponding
occupation suffix inserted. For example, Dy(III) `M#c+3#t9` means
`Dy#c+3#n2#u5#t9`. These are concrete states under registry zero-valence defaults;
empty occupation means `#n0#u0`, not an unspecified occupation.

| Element | II | III | IV |
| --- | --- | --- | --- |
| Ce | — | #u | empty |
| Pr | — | #u2 | #u |
| Nd | #u4 | #u3 | #u2 |
| Pm | — | #u4 | — |
| Sm | #u6 | #u5 | — |
| Eu | #u7 | #u6 | — |
| Gd | — | #u7 | — |
| Tb | — | #n#u6 | #u7 |
| Dy | #n3#u4 | #n2#u5 | #n#u6 |
| Ho | — | #n3#u4 | — |
| Er | — | #n4#u3 | — |
| Tm | #n6#u | #n5#u2 | — |
| Yb | #n7 | #n6#u | — |
| Lu | — | empty | — |

HW 2297–2299 supplies the III occupations and their weak ligand-field dependence;
2298–2301 supplies the related II/IV assignments. ASD verifies Nd(II/IV), Dy(II/IV)
and Tm(II). Eu(III)'s J = 0 does not remove its six unpaired f electrons. Lu retains
the existing three-valence-electron convention with 4f¹⁴ in the core; La(III) is
also empty. No alternative spin ladder is generated from magnetic moments alone.

All 34 earlier states remain supported. Ce's existing `Ce#n#u2#s` is retained:
NIST gives ¹G₄, whereas HW 2228 prints ³H₄; the registry cites ASD for the term.
The other existing comments now cite HW 2228, 2297–2301.

#### Lanthanide work table

Counts are new states after deduplication against the pre-chapter registry and
earlier rows. Ordinary anionic-ligand bonds versus anionic donation are the adopted graph
allocations (PW), as in preceding chapters. Cp− contributes three accepted pairs;
COT²− contributes five. These haptic relations are dative, not MC.

| ID | Family / states | Disposition; new | Evidence and scope |
| --- | --- | --- | --- |
| Ln01 | Nd#c+2#u4; Dy#c+2#n3#u4; Tm#c+2#n6#u; Nd#c+4#u2; Dy#c+4#n#u6 | Implemented; 5 | Chemically occurring II/IV ions, HW 2300–2304; ASD occupations. Other II/IV ions already exist. No additional gas-phase-only ionization ladder. |
| Ln02 | Ce–Lu(III): M#v3; M#c+3#t3 | Implemented; 28 | Trihalide formula units and the tris(amido)/triorganyl family, HW 2301, 2303–2304, 2309. The formula-unit graph does not assert a three-coordinate crystal structure. |
| Ln03 | Ce–Lu(III): M#c-3#v6; M#c+3#t6 | Implemented; 28 | Hexahalides, HW 2303. Hexamethyl anions of Pr/Nd/Sm/Er/Tm/Yb/Lu reuse these states (2309). |
| Ln04 | Ce–Lu(III): M#c+3#t8; M#c+3#t9 | Implemented; 28 | [Ln(en)₄]³⁺ and aqua ions, HW 2306–2307. P08Ln, pp. 3061–3062, explicitly supports an eight/nine-water equilibrium for Ho–Lu; Pm uses HW's series assignment. Fractional occupancies are not fractional graph bonds. |
| Ln05 | Tb–Lu(III): M#v3#t3 | Implemented; 7 | Monomeric Ln(dpm)₃, HW 2307. Each chelate has one ordinary Ln–O bond and one O donor; all-donor descriptions reuse Ln03. The lighter solid dimers are not assigned this monomer graph. |
| Ln06 | Ce–Lu(III): M#v3#t4; M#c+3#t7 | Implemented; 28 | Ln(dpm)₃L, L = H₂O/ROH/R₃PO, HW 2307: the seven-coordinate shift-reagent/adduct family. The series extension follows HW's generic family, not separate structural determinations for every element. |
| Ln07 | La/Ce/Sm/Gd/Tb/Dy/Er/Yb/Lu(III): M#c+#v2#t5; M#c+3#t7 | Implemented; 10 | [LnX₂(THF)₅]⁺, B01Ln p. 3177 explicitly lists the nine elements. Nine mixed descriptions plus La's donor state; the other donor states reuse Ln06. |
| Ln08 | Ce/Eu/Ho/Er(III): M#c-2#v5#t5; M#c+3#t10 | Implemented; 8 | [Ln(NO₃)₅]²⁻, ten O contacts, HW 2307. These are the listed pentanitrate ions, not a coordination ladder for every lanthanide. |
| Ln09 | Ce(III): M#c-3#v6#t6; M#c+3#t12 | Implemented; 2 | [Ce(NO₃)₆]³⁻, twelve O contacts, HW 2301, 2307; 25 counted electrons. |
| Ln10 | Pr/Sm/Gd/Tb/Dy(III): M#c-#v4#t5; Eu/Yb(III): M#c-#v4#t4 | Implemented; 7 | [Ln(EDTA)(H₂O)₃]⁻ and [Ln(EDTA)(H₂O)₂]⁻, respectively, HW 2307. Four ordinary carboxylate bonds, two N donors and three/two waters. All-donor forms reuse Ln04. |
| Ln11 | Nd#u3#v3#t9; Nd#c+3#u3#t12 | Not selected; 0 | Nd(18-crown-6)(NO₃)₃, HW 2307: six crown O and six nitrate O contacts; 27 counted electrons. Crown coordination is a common class, but this particular twelve-coordinate extension currently rests on one listed Nd example. |
| Ln12 | Er/Tm/Yb/Lu(III): M#v3#t2; M#c+3#t5 | Implemented; 8 | Ln(CH₂SiMe₃)₃(THF)₂, HW 2309. Ordinary organyl bonds or three organyl-anion donors. |
| Ln13 | Sm/Er/Yb/Lu(III): M#c-#v4; M#c+3#t4 | Implemented; 8 | [LnR₄]⁻, R = CH₂SiMe₃, tBu or 2,6-Me₂C₆H₃, HW 2309. |
| Ln14 | Nd/Sm/Eu/Dy/Tm/Yb(II): M#v2; M#c+2#t2 | Implemented; 12 | Binary dihalide formula-unit descriptions, HW 2303–2304; III electride-like CeI₂/PrI₂/GdI₂ explicitly excluded. PW graph allocation, not a claim that the crystalline dihalides contain isolated two-coordinate molecules. |
| Ln15 | II: M#v2#t(k−2); M#c+2#tk. k = 7: Nd/Sm/Dy/Tm/Yb; k = 8: Sm/Dy; k = 6: Yb | Implemented; 16 | Solvated diiodides, HW 2307; B01Ln pp. 3176–3177. Nd/Sm/Dy THF₅ and Tm mixed DME/THF give k = 7; Dy/Sm(DME)₃ give k = 8; Yb(THF)₄ gives k = 6. Yb(η²-DME)₂(η¹-DME) gives k = 7. |
| Ln16 | II: Sm/Yb M#c+2#t6/#t7/#t8; Eu M#c+2#t6/#t7 | Implemented; 4 | Cp*₂Ln and donor adducts: HW 2308, E86Ln, T80Ln. Four new states after Ln15. E86Ln establishes Eu, not HW's printed Er; no Er(II) state. T80Ln's uncertain Eu ether of crystallization does not justify Eu#c+2#u7#t8. |
| Ln17 | Ce/Pr/Tb(IV): M#v4; M#c+4#t4 | Implemented; 6 | Binary tetrafluoride formula-unit graphs, HW 2303–2304; f⁰/f¹/f⁷ occupations. As in Ln14, the graphs do not replace the extended UF₄-type crystal structures. |
| Ln18 | Ce/Pr(IV): M#c-2#v6; M#c+4#t6 | Implemented; 4 | CeCl₆²⁻ and PrF₆²⁻, HW 2302–2304. Use CeCl₆²⁻ for the discrete octahedral example; HW describes eight-coordinate Ce in (NH₄)₂CeF₆. |
| Ln19 | Ce/Pr/Nd/Tb/Dy(IV): M#c-3#v7; M#c+4#t7 | Implemented; 10 | Pentagonal-bipyramidal [LnF₇]³⁻, HW 2304. Nd/Dy(IV) are supported by these fluorocomplexes, not inferred neutral NdF₄/DyF₄ molecules. |
| Ln20 | Ce(IV): M#c-4#v8; M#c+4#t8 | Implemented; 2 | CeF₈⁴⁻, square antiprismatic, HW 2304. |
| Ln21 | Ce(IV): M#v4#t4; M#c-2#v6#t6; M#c+4#t12 | Implemented; 3 | Ce(acac)₄ and [Ce(NO₃)₆]²⁻, HW 2307. The acac all-donor description reuses Ln20. Twelve nitrate O donors count 24 electrons. |
| Ln22 | La/Ce/Nd(III): M#c+#v2#t6 | Implemented; 3 | Cp*LaI₂(THF)₃, Cp*CeI₂(THF)₃ and CpNdCl₂(THF)₃, HW 2310. One Cp− donation plus three THF donors; two ordinary halide bonds. All-donor states reuse the existing and Ln04 entries. |
| Ln23 | Er/Lu(III): M#c+2#v#t7 | Implemented; 2 | Cp₂Ln(tBu)(THF), E82Ln pp. 2009–2010; the ordinary organyl bond plus two Cp− and THF donations. The donor alternative reuses Ln04. |
| Ln24 | Sm#c+3#u5#t6#m0; Er#c+3#n4#u3#t7#m0; Lu#c+3#t7#m0 | Implemented; 3 | [Cp*₂Sm(μ-H)]₂, E83Ln; [Cp′₂Ln(THF)(μ-H)]₂, Ln = Er/Lu, E82Ln. Two separate three-center bridges, each [M,H⁻,M] = [0,2,0] (PW), existing H#c-#m2; no M–M bond. Er/Sm H positions were not located directly; the authors' assignments also use spectroscopy/chemical evidence. |
| Ln25 | Ce(COT)₂⁻; SmCp*₃; Cp₂LnCl dimers and their donor monomers | Reuse | HW 2308–2310. Ce bis-COT uses Ln08's Ce#c+3#u#t10; monomeric SmCp*₃ uses Ln04. Two chloride donors per metal in the chloride dimer use the III #t8 state and existing bridging Cl− types. |
| Ln26 | Extended hydrides, subhalides, oxides, chalcogenides, pnictides, carbides and borides | Ions / no further types | HW 2302–2306. No molecular hydride, whole-lattice MC or cluster-specific state is inferred from bulk stoichiometry. Mixed-valent/conduction-electron phases remain distinct from molecular Ln(II) formulas. |
| Ln27 | Neutral bisarenes; Sm–N₂ adduct; mixed Ln/Pt complex; Cp₄Ln⁻; Ce(IV) Cp₃OR; detailed hydrolysis products | Not selected for this pass | HW 2306, 2308–2310. Occupations, haptic allocations or complete donor counts are not established sufficiently here to justify extra default states. Formula alone does not supply them. No papers are outstanding for the selected rows; these cases would require a separately chosen extension. |
| Ln28 | Agostic and σ-H₂ intermediates | Deferred model work | HW 2309–2310. Existing bond-donation limitation; no MC substitute. |

B01Ln corrects HW's uniform seven-coordinate description of the DME solvates:
DyI₂(DME)₃ has three bidentate DME ligands (eight contacts), while TmI₂(DME)₃
has two bidentate and one monodentate DME (seven). Only those observed donor counts
are used. P08Ln supports integer eight-/nine-water graphs; its fractional crystal
occupancies remain source information rather than new atom states.

#### Limits and conformance coverage

Counts includes retained f electrons in its capacity check. For these III donor
states, the count is `f electrons + 2 × accepted pairs`. Both the ordinary and
all-donor encodings of a family give the same count. The implemented capacity
increases admit the selected families below. These are model limits, not new
claims about physical atomic capacities.

| Element | Capacity required | Limiting selected family | Disposition |
| --- | --- | --- | --- |
| Dy | 26 → 27 | III nonaaqua / EDTA: 9 + 18 | Implemented |
| Ho | 26 → 30 | III pentanitrate: 10 + 20 | Implemented |
| Er | 26 → 31 | III pentanitrate: 11 + 20 | Implemented |
| Tm | 26 → 30 | III nonaaqua: 12 + 18 | Implemented |
| Yb | 26 → 31 | III nonaaqua: 13 + 18 | Implemented |
| Nd | 26 → 27 | III crown trinitrate: 3 + 24 | Not selected; capacity remains 26 |

Ce–Tb and Lu otherwise fit 26; the three La additions fit 20. Lu's lower electron
count follows its existing f¹⁴-core convention. The three selected MC metal states
fit their capacities. Nd's fallback charge maximum (3 → 4) and Dy's fallback
unpaired maximum (5 → 6) are implemented to align invariant enumeration with
Ln01/Ln19; explicit values bypass those fallback bounds. Valence-electron counts
and ValenceTable are unchanged.

All 232 additions have representative coverage, including both metal sites of
the three hydride dimers. These inputs supply the source unpaired count while
leaving paired occupation to resolution; AtomTyping and Counts agree under both
Strict and MostSaturated. Fourteen additional trihalide inputs leave both occupation
counts unknown: AtomTyping infers the selected fⁿ state, while Counts uses its
minimum unpaired count by parity. All 246 fixtures pass, with 996 independent
metal-state and molecule electron-invariant checks across the four configurations.
Existing bare-atom/ion fixtures stay in place. No new ligand types are needed;
the known haptic-donor-side limitation remains deferred. All 3,113 earlier snapshots
remain unchanged. The full conformance, unit and property suites, formatting and
Clippy pass. No selected work remains outstanding for this chapter.

#### Lanthanide sources and access

All six primary papers are available, reviewed and archived. HW references above
are printed Band II pages. Extraction and the expanded candidate/count audit are
intermediates in `scratch/hw-lanthanide-review/`; the tables here are the durable record.

| Key | Source | Assignment-specific evidence |
| --- | --- | --- |
| E86Ln | Evans, Hughes and Hanusa, *Organometallics* **5** (1986), 1285–1291. [DOI](https://doi.org/10.1021/om00138a001), [PDF](../materials/valence/Evans-1986-Sm-Eu-metallocenes.pdf) | SmCp*₂ and EuCp*₂ structures; explicit f⁶/f⁷ occupations, p. 1286. Corrects the Er label in HW 2308. |
| T80Ln | Tilley et al., *Inorg. Chem.* **19** (1980), 2999–3003. [DOI](https://doi.org/10.1021/ic50212a031), [archived author manuscript](../materials/valence/Tilley-1980-Eu-Yb-metallocenes.pdf), [OSTI](https://www.osti.gov/servlets/purl/1030602) | Eu/Yb metallocene ether adducts; Eu f⁷ and diamagnetic Yb f¹⁴, manuscript pp. 10–11. The additional Eu ether may be lattice solvent; it is not counted as a confirmed second donor. |
| B01Ln | Bochkarev et al., *Angew. Chem. Int. Ed.* **40** (2001), 3176–3178. [DOI](https://doi.org/10.1002/1521-3773%2820010903%2940%3A17%3C3176%3A%3AAID-ANIE3176%3E3.0.CO%3B2-Y), [PDF](../materials/valence/Bochkarev-2001-Nd-Tm-diiodides.pdf) | NdI₂(THF)₅ and TmI₂(DME)₂(THF) structures; explicit six-/seven-/eight-coordinate comparisons and III dihalide-THF cations, pp. 3176–3177. II occupations use HW/ASD; this paper is not a spin determination. |
| P08Ln | Persson et al., “Hydration of Lanthanoid(III) Ions in Aqueous Solution and Crystalline Hydrates Studied by EXAFS Spectroscopy and Crystallography,” *Chem. Eur. J.* **14** (2008), 3056–3066. [DOI](https://doi.org/10.1002/chem.200701281), [PDF](../materials/valence/Persson-2008-lanthanide-hydration.pdf) | Eight-/nine-water equilibrium for Ho–Lu, pp. 3061–3062; no abrupt coordination-number break at Gd. |
| E82Ln | Evans et al., “Organolanthanide Hydride Chemistry. 1,” *J. Am. Chem. Soc.* **104** (1982), 2008–2014. [DOI](https://doi.org/10.1021/ja00371a035), [PDF](../materials/valence/Evans-1982-lanthanide-hydrides.pdf) | Er/Lu/Y hydride dimers; directly located bridge H in the Y structure, inferred Er bridges, and Lu spectroscopic/chemical characterization. The precursor organyls also support Ln23. Integer bridge allocation is PW. |
| E83Ln | Evans, Bloom, Hunter and Atwood, “Organolanthanide Hydride Chemistry. 3,” *J. Am. Chem. Soc.* **105** (1983), 1401–1403. [DOI](https://doi.org/10.1021/ja00343a071), [PDF](../materials/valence/Evans-1983-samarium-hydride.pdf) | Sm hydride dimer, Sm(III) evidence, and explicit exclusion of a Sm–Sm bond, p. 1402. H sites were not located; bridge assignment combines the dimer structure, spectroscopy and comparison. Integer allocation is PW. |
| ASD | NIST atomic data / Atomic Spectra Database. [Ce ground state](https://physics.nist.gov/PhysRefData/Handbook/Tables/ceriumtable1.htm), [ion query](https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html); archived [Ce](../materials/valence/NIST-Ce-atomic-data.html), [Nd](../materials/valence/NIST-Nd-ground-states.html), [Dy](../materials/valence/NIST-Dy-ground-states.html), [Tm](../materials/valence/NIST-Tm-ground-states.html) | Ce's existing singlet term; Nd²⁺ f⁴/⁵I₄, Nd⁴⁺ f²/³H₄, Dy²⁺ f¹⁰/⁵I₈, Dy⁴⁺ f⁸/⁷F₆, Tm²⁺ f¹³/²F₇/₂. Retrieved 2026-11-04. |

### Actinides

Reviewed 2026-11-04: HW II, chapter XXXVI, pp. 2312–2346; bibliography
2347–2348. **Implemented: 243 new states and 243 fixtures.**
The unselected rows below are excluded from that count. Element-organized
fixtures are in `actinides/`; organometallic and explicit MC inputs belong in
their respective bonding categories. Added-state coverage, rather than a fixture
for every extracted compound, remains the requirement.

#### Occupations used in the table

Insert the appropriate suffix into each `M` pattern, as in the lanthanide table.
For example, Np(IV) `M#c+4#t12` expands to `Np#c+4#u3#t12`.
Empty means zero pairs and zero unpaired electrons under the registry defaults.

| Element | II | III | IV | V | VI |
| --- | --- | --- | --- | --- | --- |
| Ac | — | empty | — | — | — |
| Th | — | #u | empty | — | — |
| Pa | — | #u2 | #u | empty | — |
| U | — | #u3 | #u2 | #u | empty |
| Np | — | #u4 | #u3 | #u2 | #u |
| Pu | — | #u5 | #u4 | #u3 | #u2 |
| Am | #u7 | #u6 | #u5 | #u4 | #u3 |
| Cm | — | #u7 | #u6 | — | — |
| Bk | — | #n#u6 | #u7 | — | — |
| Cf | #n3#u4 | #n2#u5 | #n#u6 | — | — |
| Es | #n4#u3 | #n3#u4 | — | — | — |
| Fm | #n5#u2 | #n4#u3 | — | — | — |
| Md | #n6#u | #n5#u2 | — | — | — |
| No | #n7 | #n6#u | — | — | — |
| Lr | — | empty | — | — | — |

HW 2228, 2320–2331 supplies the configurations and chemical oxidation states;
ASD checks the missing ions and existing configuration comments. These are the
retained-electron occupations used for graph counting, not total-J assignments.
In particular, C19An finds a nonmagnetic crystal-field singlet within predominantly
³H₄ U(IV); it does not justify replacing f² `#u2` with a paired `#n` state.
Lr retains the existing f¹⁴-core convention. Th(III)'s single electron gives `#u`
whether a particular complex favors 5f or 6d character.

The 51 earlier Ac–Lr rows are retained. Comment corrections:
U atom s²d¹f³; Np atom s²d¹f⁴; Am atom s²f⁷; Pu⁵⁺ labeled Pu rather than Np;
free Th³⁺ 5f¹ rather than 6d¹; Lr atom 7s²7p¹ with f¹⁴ in the core. HW/ASD
citations replace the unspecific or missing references. Existing Pu⁷⁺ `#u`
is retained as the formal f¹ chemical occupation; An25's Np(VII) compounds use f⁰.
No new bare Np⁷⁺ row is inferred from this chemical assignment.

#### Actinide work table

Counts are new states after deduplication against the existing registry and earlier
rows. Ordinary anionic-ligand bonds and anionic donation use the adopted allocation
(PW). Cp− contributes three accepted pairs, COT²− five, and allyl− two; these are
dative descriptions. Actinyl donor alternatives retain the ordinary M=O core.
Binary formula-unit graphs follow the earlier chapter convention and do not assert
the coordination number of an isolated molecule or its crystal lattice.

| ID | Family / states | Disposition; new | Evidence and scope |
| --- | --- | --- | --- |
| An01 | Pa#c+3#u2; Am#c+5#u4; Am#c+6#u3; No#c+3#n6#u | Implemented; 4 | Chemically occurring ions, HW 2325, 2328–2330, 2345; ASD occupations. No gas-phase ionization ladder. |
| An02 | Am/Cf/Es/Fm/Md/No(II): M#v2; M#c+2#t2 | Implemented; 12 | Dihalide formula units, explicitly identified as a family in HW 2336; table 36.12 gives Am/Cf examples. No Th(II) assignment from metallic ThI₂. |
| An03 | Ac and U–Es(III): M#v3; M#c+3#t3 | Implemented; 18 | Trihalide formula units, HW 2324, 2336–2337; tris(amido)/triorganyl descriptions reuse these where present. No Th/Pa(III) assignment from conducting iodides. |
| An04 | Th–Cf(IV): M#v4; M#c+4#t4 | Implemented; 18 | Binary tetrafluorides, HW 2336–2337. ThO₂, U(NPh₂)₄ and ordinary tetraorganyl descriptions reuse these counts. |
| An05 | Pa/U/Np(V): M#v5; M#c+5#t5 | Implemented; 6 | Binary pentahalide formula units, HW 2336–2338. Their solid-state bridges are not represented by an invented MC state. |
| An06 | U/Np/Pu(VI): M#v6; M#c+6#t6 | Implemented; 6 | Molecular UF₆, UCl₆, NpF₆ and PuF₆, HW 2336–2337; six terminal ligands. |
| An07 | U–Bk(III): M#c-3#v6; M#c+3#t6 | Implemented; 12 | Octahedral hexahalides, HW 2324, 2326, 2339. |
| An08 | Th/Pa/U/Np/Pu/Am/Cm/Bk(IV): M#c-2#v6; M#c+4#t6 | Implemented; 16 | Hexachloro/bromo/iodo complexes of Th–Pu/Bk, HW 2339; AmF₆²⁻ and CmF₆²⁻ in fluoride-rich solution, 2327. Other solid AnF₆²⁻ salts can have shared F and higher coordination. |
| An09 | U(IV): M#c-4#v8; M#c+4#t8 | Implemented; 2 | Discrete UF₈⁴⁻ in (NH₄)₄UF₈, C19An; f² occupation. HW 2339's generic fluoride formulas are not expanded into every element/coordination combination. |
| An10 | Pa/U/Np/Pu(V): M#c-#v6; M#c+5#t6. Pa/U/Np(V): M#c-3#v8; M#c+5#t8 | Implemented; 14 | Hexafluoride and cubic Na₃AnF₈ examples, HW 2329, 2339. Counterion-dependent heptafluoride structures remain outside these molecular assignments. |
| An11 | U(VI): M#c-2#v8; M#c+6#t8 | Implemented; 2 | UF₈²⁻, HW 2339; HS00An pp. 70–71, Fig. 3 and Table 1 establish eight terminal U–F contacts in single-crystal (NO)₂UF₈, with square-antiprismatic geometry. |
| An12 | Ac(III): M#c+3#t9. U–Cf(III): M#c+3#t8; M#c+3#t9 | Implemented; 15 | Aqua ions, HW 2326. F21An supports Ac nonaqua and reviews Am/Cm hydration. U–Cf follows HW's eight/nine-water family; no separate claim of an individual structure determination for every member. |
| An13 | Th/Pa/U/Np/Pu/Bk(IV): M#c+4#t8; M#c+4#t9 | Implemented; 11 | Aqua ions, HW 2326–2327; U#c+4#u2#t8 already counted in An09. No corresponding free aqua complexes inferred for fluoride-stabilized Am/Cm/Cf(IV). |
| An14 | Th–Pu(IV): M#c-4#v8; M#v4#t4; M#c+4#t8 | Implemented; 9 | [An(NCS)₈]⁴⁻ and An(acac)₄, HW 2324, 2328, 2344. Th/U/Np/Pu dithiocarbamates reuse these chelate counts. |
| An15 | Th/U/Np/Pu(IV): M#c-2#v6#t6; M#c+4#t12 | Implemented; 8 | [An(NO₃)₆]²⁻: six bidentate nitrates, HW 2324, 2344. Uses the implemented Np/Pu capacity extensions. |
| An16 | Th/U/Pu(IV): M#c-6#v10; M#c+4#t10 | Implemented; 6 | [An(CO₃)₅]⁶⁻, ten O contacts, HW 2344. Two ordinary anionic O bonds per carbonate, or their donor alternatives. |
| An17 | Th(IV): M#v4#t7; M#c+4#t11; M#v4#t6 | Implemented; 3 | Th(NO₃)₄(H₂O)₃ and Th(NO₃)₄(OPPh₃)₂, HW 2324, 2328. Only coordinated water counts; the second all-donor state reuses An16. |
| An18 | U(IV): M#v4#t; M#c+4#t5; M#v4#t2 | Implemented; 3 | [U₂(NEt₂)₈] has five N contacts per U; UCl₄(OPR₃)₂ has six, HW 2324, 2328. Amido bridges use existing N donor types; no MC needed. |
| An19 | Pa(V): M#c+2#v3; M#c+#v4 | Implemented; 2 | PaO(OH)²⁺ and PaO(OH)₂⁺ cores, HW 2328. PaO(OH)₄⁻ reuses An10's ordinary state. Hydration counts are not supplied by these formulas. |
| An20 | U/Np/Pu/Am(V): M#c+#v4. U/Np/Pu/Am(VI): M#c+2#v4 | Implemented; 8 | Actinyl cores O=M=O, HW 2328–2329. HW's O≡U≡O²⁺ alternative uses U#v6 from An06 and existing O⁺ triple-bond types; it adds no state. |
| An21 | U/Np/Pu/Am(VI): M#c+2#v4#t5 | Implemented; 4 | Actinyl pentaaqua ions, HW 2329–2330. |
| An22 | U/Np/Pu/Am(V): M#c-#v6; M#c+#v4#t2. VI: M#c-2#v8; M#c+2#v4#t4 | Implemented; 12 | AnO₂(OH)₂⁻ and AnO₂(OH)₄²⁻, HW 2328–2330. Hydroxido bonds ordinary or donated; oxo bonds remain ordinary. |
| An23 | U(VI): M#c-3#v9 | Implemented; 1 | UO₂(NCS)₅³⁻, HW 2324; donor alternative reuses An21. UO₂Cl₄²⁻/Br₄²⁻ reuse An22. |
| An24 | U(VI): M#v6#t4; M#c-#v7#t3. U/Np/Pu/Am(V): M#c-5#v10; M#c+#v4#t6. VI: M#c-4#v10; M#c+2#v4#t6 | Implemented; 18 | UO₂(NO₃)₂(H₂O)₂ and UO₂(NO₃)₃⁻, HW 2343. U/NpO₂(O₂)₃⁴⁻, 2342–2343. Actinyl tricarbonate V/VI family, 2344, expanded to the U–Am actinyl series. Ordinary carbonate/peroxide bonds and donor alternatives reuse the same states. |
| An25 | Np#c-3#v10; Np#c-#v8#t2 | Implemented; 2 | [NpO₄(OH)₂]³⁻, B24An pp. 16170–16172, 16182: four equatorial oxo and two axial hydroxido ligands, supported in solution and the cobalt-hexammine salt. Replaces HW 2330–2331's dioxo formulation. Ordinary Np=O with ordinary or donated OH⁻ is PW; no MC. Protonation intermediates and Pu analogues add no prospective rows. |
| An26 | Th/Pa(III): M#c+3#t9. U/Np/Pu/Am(III): M#c+3#t10 | Implemented; 6 | Tris-Cp family Ac/Th–Cf, Cp₃An(THF) for U/Np/Pu, and [An(COT)₂]⁻ for U/Np/Pu/Am, HW 2345–2346. Other #t9 states reuse An12. D17An establishes monomeric substituted tris-Cp Np; the parent Cp₃Np crystal has additional bridging contacts. |
| An27 | Pa(IV): M#c+4#t12. Pa/Np(IV): M#c+4#t10 | Implemented; 3 | Cp₄An (Th/Pa/U/Np) and An(COT)₂ (Th–Pu), HW 2345–2346. Other counts reuse An15/An16. D17An confirms four η⁵ rings in Cp₄Np; no MC encoding. |
| An28 | Th/U/Np(IV): M#c+3#v#t9. Th/U(IV): M#c+2#v2#t6. U(III): M#c+2#v#t9 | Implemented; 6 | Cp₃AnX, Cp₂AnX₂ and Cp₃UCl⁻, HW 2345–2346. Use established halide/organyl representatives; all-donor alternatives already counted. Hydride topology is reviewed separately in An33. |
| An29 | Th(IV): M#v4#t2; M#c-3#v7; M#c+4#t7 | Implemented; 3 | Th/U R₄(dmpe) and ThMe₇³⁻, HW 2345. U's corresponding phosphine state is in An18; U hexorganyl anions reuse An08. |
| An30 | Th(III): M#c+3#t8 | Implemented; 1 | [Th(η³-allyl)₄]⁻, HW 2345. Neutral Th/U tetraallyls reuse IV #t8; four two-pair allyl donations. |
| An31 | U(V): M#c+3#v2#t9 | Implemented; 1 | (C₅H₄Me)₃U=NAr family, HW 2346: ordinary U=N plus three Cp donations. No whole-assembly MC state. |
| An32 | Th#v3#h; U#u2#v3#h | Implemented; 2 | Terminal HM[N(SiMe₃)₂]₃, HW 2345; S81An explicitly establishes monomeric hydrides and the f² U occupation. Explicit-H and all-donor states reuse An04. |
| An33a | Th/U(IV): M#c+3#v#t6#m0; M#c+3#h#t6#m0; M#c+4#t7#m0 | Implemented; 6 | [Cp*₂M(μ-H)H]₂, F81An pp. 6660–6662: two separate M–H–M bridges and one terminal hydride per metal. Th bridge H located by neutron diffraction; U dimer assigned by the authors from spectroscopy. Each bridge uses existing H⁻ and contributions [0,2,0] (PW). Terminal H is ordinary, implicit or donated. Th chlorohydride dimer reuses these states. |
| An33b | Th#c+2#v#h#t6 | Implemented; 1 | Monomeric Cp*₂Th(OtBu)H, F81An p. 6663: implicit terminal H with ordinary Th–O; explicit-H and all-donor forms reuse An28/An27. |
| An34a | U#c-#v7; U#c+6#t7 | Implemented; 2 | UF₇⁻, HW 2339; B76An pp. 2533–2534 identifies the ion in HF solutions by Raman spectroscopy and fluoride-addition equilibria. Seven ordinary U–F bonds or seven F⁻ donations use PW. Exact geometry and counterion-dependent distortions remain distinct from this connectivity assignment; no MC state or solvent-adduct type. |
| An34b | UF₉³⁻ | Not selected | HW 2339; I81An establishes K₃UF₉/Rb₃UF₉ compositions, U(VI) and orthorhombic unit cells, but supplies no atomic positions or discrete nine-coordinate assignment. No new state or fixture from the stoichiometry alone. |
| An34c | UF₁₀⁴⁻ | Not selected by user | D84An and M90An p. 47 interpret the spectroscopy as a monomeric, distorted UF₁₀⁴⁻ ion; single-crystal attempts failed, and the powder data give a unit cell only. Leave this formulation out for now. U#c-4#v10 is independently covered by An24; no U#c+6#t10 addition. |
| An35 | Cp₃U⁺; Cp₃U(CO); further III/IV organyls and oxohalides | Reuse where occupation and graph are established | HW 2339–2346. These need not all become fixtures. A new representation is not added solely to mirror each formula in the chapter. |
| An36 | Borohydrides, extended binary phases and hydrolysis aggregates; inverse U arenes, U(C₇H₇)₂⁻, Cp₄An⁻, Cp₂U(MeCN)₅²⁺, matrix U(CO)₆ | Not selected for new default states | HW 2336–2346. No whole-lattice/cage MC states or speculative occupations. D17An also corrects HW's U(II) description of the inverse arene to U(III)/arene²⁻. Fm–Lr bonded trihalide/aqua ladders are not added from the generic series statements. |

The 243 additions comprise Ac 3, Th 27, Pa 20, U 54, Np 39, Pu 33, Am 26,
Cm 10, Bk 12, Cf 8, Es 4, Fm 2, Md 2 and No 3. Lr adds none.
The six MC additions in An33a describe separate 3c2e hydride bridges; no
metal–metal bond or whole-assembly MC encoding is used. F81An's U(III) hydride
aggregate has unresolved size/topology, and monomeric Cp*₂ThH₂ appears only as
an exchange-mechanism possibility; neither adds states. Existing companion types
suffice for every fixture.

#### Limits and coverage

For these states, Counts checks retained electrons plus twice the ordinary bond
orders, implicit H and accepted pairs. The implemented extensions are:

| Element | Capacity | Limiting family | Status |
| --- | --- | --- | --- |
| Np | 26 → 27 | IV hexanitrate and tetrakis-Cp: 3 + 24 | Implemented |
| Pu | 26 → 28 | IV hexanitrate: 4 + 24 | Implemented |
| Cf | 26 → 27 | III nonaqua: 9 + 18 | Implemented |

All other selected states fit the previous capacities. No's unpaired-electron
enumeration maximum is 0 → 1 for No(III), An01; its capacity remains 26.
ValenceTable and valence-electron counts are unchanged.

The 243 independent fixtures cover every added state: 216 in `actinides/`,
21 in `organometallics/` and six in `multicenter_bonds/`. Source occupations are
supplied where needed to distinguish them from Counts' minimum-unpaired fallback.
All four model/policy combinations give the selected states and pass the valence
invariants check: 1,008 metal-site checks across 252 sites. All 3,359 earlier
resolution snapshots are unchanged; all 3,602 available resolution cases pass.
Graph unit (3,594), graph property (223) and chemistry (316) tests pass, as do
formatting and Clippy. The selected chapter scope is complete.

#### Actinide sources and access

Available primary files are archived below. Intermediates, including the expanded
state list, are in `scratch/hw-actinide-review/`. All requested papers are reviewed
and archived. No source-access blocker remains for this chapter; An34b/c record
the unselected higher-fluoride cases after full-paper review.

| Key | Source | Evidence / access |
| --- | --- | --- |
| F21An | Ferrier et al., “Advancing understanding of actinide(III) (Ac, Am, Cm) aqueous complexation chemistry,” *Chem. Sci.* (2021). [DOI](https://doi.org/10.1039/D1SC00233C), [archived full text](../materials/valence/Ferrier-2021-actinide-hydration.xml) | Full primary text reviewed. Ac nine-water assignment and Am/Cm comparisons; does not support an extra Ac eleven-water type from the older measurement. |
| D17An | Dutkiewicz et al., “Reduction chemistry of neptunium cyclopentadienide complexes: from structure to understanding,” *Chem. Sci.* **8** (2017), 2553–2561. [DOI](https://doi.org/10.1039/C7SC00034K), [author manuscript](../materials/valence/Dutkiewicz-2017-neptunium-cyclopentadienides-author-manuscript.pdf), [final full text](../materials/valence/Dutkiewicz-2017-neptunium-cyclopentadienides.xml) | Full text reviewed. Np(IV) Cp₄ and Np(III) Cp/donor structures; distinguishes monomeric substituted tris-Cp from parent Cp₃Np's additional contacts. |
| C19An | Chemey et al., “Electronic, Magnetic, and Theoretical Characterization of (NH₄)₄UF₈, a Simple Molecular Uranium(IV) Fluoride,” *Inorg. Chem.* **58** (2019), 637–647. [DOI](https://doi.org/10.1021/acs.inorgchem.8b02800), [PDF](../materials/valence/Chemey-2019-octafluorouranate-IV.pdf) | Full paper reviewed. Discrete UF₈⁴⁻; f²/³H₄ and crystal-field splitting, rather than a new paired occupation. |
| S81An | Simpson, Turner and Andersen, “Preparation and Hydrogen-Deuterium Exchange of Alkyl and Hydride Bis(trimethylsilyl)amido Derivatives of the Actinide Elements,” *Inorg. Chem.* **20** (1981), 2991–2995. [DOI](https://doi.org/10.1021/ic50223a047), [author manuscript](../materials/valence/Simpson-1981-actinide-monohydrides-author-manuscript.pdf) | Full manuscript reviewed; pp. 3–4 establish monomeric terminal Th/U hydrides, diamagnetic Th and paramagnetic f² U. |
| HS00An | Hwang and Seppelt, “The structures of ReF₈⁻ and UF₈²⁻,” *J. Fluorine Chem.* **102** (2000), 69–72. [DOI](https://doi.org/10.1016/S0022-1139(99)00248-1), [PDF](../materials/valence/Hwang-Seppelt-2000-octafluoroanions.pdf) | Full paper reviewed. Single-crystal (NO)₂UF₈ has a discrete square-antiprismatic anion, with eight U–F distances of 207.0–215.1 pm, pp. 70–71. |
| F81An | Fagan et al., “Synthesis and properties of bis(pentamethylcyclopentadienyl) actinide hydrocarbyls and hydrides,” *J. Am. Chem. Soc.* **103** (1981), 6650–6667. [DOI](https://doi.org/10.1021/ja00412a021), [PDF](../materials/valence/Fagan-1981-actinide-hydrocarbyls-hydrides.pdf) | Full paper reviewed. Hydride dimers and terminal alkoxohydride, pp. 6660–6663; explicit Th(IV) f⁰ / U(IV) f², p. 6665. Th neutron structure and U spectroscopic assignment are distinguished in An33a; integer bridge allocation is PW. |
| B24An | Benthin et al., “Probing the protonation and reduction of heptavalent neptunium with computational guidance,” *Dalton Trans.* **53** (2024), 16170–16185. [DOI](https://doi.org/10.1039/D4DT01706D), [PDF](../materials/valence/Benthin-2024-neptunium-VII.pdf) | Full paper reviewed. Reproduces the tetraoxo cobalt-hexammine salt and supplies complementary solution/solid Raman evidence; confirms An25's four-oxo/two-hydroxido formulation. |
| D84An | Družina, Milićev and Slivnik, “Synthesis and characterization of a new uranium(VI) compound with a co-ordination number larger than eight, ammonium decafluorouranate(VI),” *J. Chem. Soc., Chem. Commun.* (1984), 363–364. [DOI](https://doi.org/10.1039/C39840000363), [PDF](../materials/valence/Druzina-1984-ammonium-decafluorouranate.pdf) | Full paper reviewed. Composition, U(VI) oxidation and diamagnetism established; Raman supports coordination above eight. Distorted UF₁₀⁴⁻ is the authors' interpretation, without a structure determination. |
| B76An | Bougon, Charpin, Desmoulin and Malm, “Structural study and properties of the alkali metal, nitrosyl, and ammonium hepta- and octafluorouranates(VI),” *Inorg. Chem.* **15** (1976), 2532–2540. [DOI](https://doi.org/10.1021/ic50164a043), [PDF](../materials/valence/Bougon-1976-hepta-octafluorouranates-VI.pdf) | Full paper reviewed. Solution Raman and fluoride-addition equilibria identify UF₇⁻ and UF₈²⁻, pp. 2533–2534. NMR shows exchange, not a unique geometry; pp. 2537–2539 discuss temperature/counterion-dependent distortions. Ordinary/donor graph allocation is PW. |
| I81An | Iwasaki, Ishikawa, Ohwada and Fujino, “Syntheses and X-ray observations of new uranium fluoride complexes K₃UF₉ and Rb₃UF₉,” *Inorg. Chim. Acta* **54** (1981), L193–L194. [DOI](https://doi.org/10.1016/S0020-1693(00)95431-8), [PDF](../materials/valence/Iwasaki-1981-nonafluorouranates-VI.pdf) | Full paper reviewed. Composition, U(VI) absorption and indexed powder patterns; no atomic-coordinate structure or molecular connectivity assignment. |
| M90An | Milićev, Družina and Volavšek, “Synthesis, thermal decomposition and powder diffraction data of ammonium decafluorouranate(VI), (NH₄)₄UF₁₀,” *Polyhedron* **9** (1990), 47–51. [DOI](https://doi.org/10.1016/S0277-5387(00)84246-1), [PDF](../materials/valence/Milicev-1990-ammonium-decafluorouranate.pdf) | Full paper reviewed. Repeats the spectroscopic monomeric UF₁₀⁴⁻ interpretation and reports unsuccessful single-crystal attempts, p. 47. Powder indexing and decomposition data add no atomic-coordinate determination. |
| ASD | NIST Atomic Spectra Database. [Ion query](https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html); archived [Th](../materials/valence/NIST-Th-ground-states.html), [Pa](../materials/valence/NIST-Pa-ground-states.html), [Np](../materials/valence/NIST-Np-ground-states.html), [Am](../materials/valence/NIST-Am-ground-states.html), [No](../materials/valence/NIST-No-ground-states.html), [Lr](../materials/valence/NIST-Lr-ground-states.html) | Retrieved 2026-11-04. Free Th³⁺ 5f¹, Pa³⁺ 5f², Am⁵⁺/Am⁶⁺ f⁴/f³, No³⁺ f¹³, and Lr's 7p ground-state electron. Higher-charge free-ion core-hole configurations are not chemical oxidation-state assignments. |
