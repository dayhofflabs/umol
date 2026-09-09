# 217 — Rhea participant rejection census

Status: In Progress
Date: 2026-08-30
Relates: [153](153-format-parsing-outstanding-tasks-2026-07-18.md),
[168](168-api-hygiene-2026-07-27.md),
[174](174-aromatic-hydrogen-resolution-2026-07-31.md),
[216](216-canonicalization-performance-2026-08-30.md),
[data types guide](../docs/development/data-types.md)

## Purpose

Doc 216 used the ChEBI participants of Rhea release 141 to construct a molecule corpus for
canonicalization measurements. This document classifies every participant that did not reach the
stored resolved cohort, first for release 141 and then for release 142. Abstract participants,
unsupported chemistry models, and explicitly unknown stereochemistry are legitimate exclusions from
a corpus that requires concrete molecules; the census separates those from defects.

The census established one implementation defect, the endpoint-blind reading of V2000 wedges, and
found the resolver's diagnostics insufficient to classify its outcomes without returning to the
source records. This document owns the wedge representation contract and its implementation. The
remaining findings are recorded as evidence with the open work listed at the end; no admission
policy, chemistry-model extension, or raising design is settled here.

## Release-141 census

The census used the same inputs and configuration as doc 216 S0a:

| Input | Release-141 evidence |
| --- | --- |
| Rhea release | 141, dated 2026-06-10 |
| Molfile archive | `rhea-mol.tar.gz`, 5,172,329 bytes, xxh3-128 `ef5cb65c0becec63f39da049e3639566` |
| participant list | `chebiId_name.tsv`, 589,145 bytes, xxh3-128 `2e530735f3d6e0b185a3db1e71929bb7` |
| CTfile configuration | `CtfileIoConfig::basic()` |
| chemistry model | `ValenceModel::mdl()`, Daylight aromaticity, default stereo model |
| resolver configuration | `ResolveConfig::default()` |

The manifest contains one outcome for every listed participant. Messages were normalized by stage,
the first rejected CTfile line was inspected for every parse rejection, and every reported stereo
atom was correlated with its raw atom parity and incident wedge records, retaining the raw first
bond atom because it is the pointed end of a V2000 wedge. The census records the first observed
failure under one configuration, not every defect in a record; fixing an earlier failure can expose
a later one.

The manifest reports 12,950 listed participants, 14,048 Molfiles, 12,942 joined participants, and
10,107 resolved records:

| Stage | Rejected | Complete message classification |
| --- | ---: | --- |
| source read | 8 | missing listed Molfile |
| CTfile parse | 1,464 | 1,448 extended atom symbols; 16 `M  ZZC` editor properties |
| TableIR raise | 56 | 55 wedge conflicts; 1 tetrahedral ligand-count failure |
| resolution underdetermined | 378 | explicit undetermined atom or bond stereo |
| resolution contradictory | 937 | 894 tetrahedral stereo failures; 43 valence failures |
| **Total** | **2,843** | all listed participants not stored in the resolved cohort |

There were no `resolve_execution`, `resolve_non_ground`, unexpected, or aromatization failures. All
10,107 resolved records also aromatized successfully.

### Source and parse exclusions

The eight source-read failures are participant-list entries without matching archive members:
CHEBI 4194, 24431, 59720, 84139, 131859, 137328, 140159, and 229467. They are generic classes such
as “a D-hexose”. The archive separately contains 1,106 Molfiles that are not in the participant
list; those are outside this census.

The first rejected atom line accounts for 1,448 parse outcomes:

| First atom symbol | Count | Source meaning |
| --- | ---: | --- |
| `R` | 985 | R group |
| `R#` | 416 | unlabeled R group |
| `R1` | 39 | numbered R group |
| `X` | 4 | halogen wildcard |
| `A` | 2 | heavy-atom wildcard |
| `hv` | 1 | photon pseudoatom |
| `e` | 1 | electron pseudoatom |

These are expected exclusions under `CtfileIoConfig::basic()`. The extended parser has
representations for R groups, query wildcards, and pseudoatoms under the corresponding flags. CHEBI
685, an R-group phospholipid, and CHEBI 10545, an electron, are representative records. The
remaining 16 parse outcomes stop on `M  ZZC`, an ACD/ChemSketch label property that the parser
accepts only under `EDITOR_EXTENSIONS`; CHEBI 11302 and CHEBI 11909 are representative.

### TableIR raise findings

Fifty-five records report `inconsistent wedge bonds at atom N`; one record reports a tetrahedral
atom with two ligands. Comparing the reported atom with the raw first and second atoms of every
incident wedge separates the cases:

| Reported-site relation to raw wedge | Count | Interpretation |
| --- | ---: | --- |
| only the wide endpoint | 20 | wedge applied to an endpoint it does not describe |
| both a narrow and a wide endpoint | 31 | an unrelated incoming wedge pollutes an actual focus |
| only the narrow endpoint | 4 | multiple outgoing wedges disagree under the current geometric reading |
| explicit atom parity on a two-neighbor atom | 1 | source tetrahedral marker has no valid frame |

The first two rows are the wedge endpoint defect. CHEBI 15393, (1R,2S,4R)-borneol, has two wedges
whose raw first atoms are distinct stereo sites and whose common second atom was treated as another
focus. CHEBI 58540 has one incoming and one outgoing wedge at the reported atom; only the outgoing
wedge belongs to that atom. The four narrow-only conflicts are CHEBI 34596, 35306, 219836, and
231826: two raw outgoing wedges imply opposite configurations under the current coordinate
calculation. CHEBI 80541 supplies explicit parity on an atom with two neighbors.

### Resolution findings

All 378 underdetermined outcomes contain a source marker whose configuration is explicitly unknown:

| Source markers in one record | Records |
| --- | ---: |
| atom parity 3 only | 363 |
| double-bond stereo code 3 only | 12 |
| both | 3 |

The raised `#T` or `#C` assertion contains an undetermined coset and the stereo resolver correctly
publishes no concrete result. CHEBI 7 is representative of atom parity 3 and CHEBI 23316 of
crossed double-bond stereo. Every returned `ResolveReport` has an empty `unresolved` map and 29
also have no tie-break entries; the report does not state that partial stereo caused the outcome.

Of the 894 tetrahedral stereo contradictions, 886 name the raw wide endpoint of a wedge: the
source describes stereo at the wedge's first atom, but raise also created a `#T` assertion at its
second atom. CHEBI 40, (+)-pinoresinol, is representative. Five name the raw narrow endpoint of an
outgoing wedge on a carbon with three neighbors and localized-bond valence four, and three come
from explicit atom parity on the same three-neighbor, valence-four shape. Those eight sites are
trigonal under the localized source structure; CHEBI 32446 and CHEBI 34463 are representative.

All 43 `no matching valence state` outcomes contain a transition metal outside the frozen MDL
counts table:

| Transition metal present | Records |
| --- | ---: |
| Fe | 20 |
| Co | 18 |
| Mn | 2 |
| V | 1 |
| Cu | 1 |
| Cr | 1 |

The set includes isolated ions such as CHEBI 29034, iron and manganese oxides, and coordination
complexes such as CHEBI 4991 and CHEBI 16304. The contradiction message does not name the atom,
element, charge, or valence that failed even though the counts module has `CountsMismatch`
vocabulary. No resolver outcome reports an aromaticity contradiction.

### Root-cause disposition

| Disposition | Records | Included causes |
| --- | ---: | --- |
| expected boundary, source, or selected-model exclusion | 1,893 | 8 missing records, 1,464 basic-parse exclusions, 378 unknown stereo outcomes, 43 transition-metal valence exclusions |
| confirmed umol wedge-endpoint defect | 937 | 51 raise conflicts and 886 stereo contradictions |
| source/perception compatibility requiring an independent reading | 13 | 4 conflicting outgoing-wedge depictions, 1 parity site with two neighbors, 8 tetrahedral markers on trigonal carbon |
| **Total** | **2,843** | complete release-141 participant rejection census |

## Wedge endpoint defect

V2000 wedge direction is endpoint-relative: the CTfile specification places the pointed end of a
stereo bond at the first atom of the bond line. At the census baseline TableIR documented this rule
on the wedge enum, but `Bond::new` and `ExtendedBond::new` normalize the endpoints through
`AtomPair::new`, so the retained `first()` is the smaller atom id rather than the source first
atom. More decisively, `wedge_bond_neighbors` treats every incident `Up` or `Down` wedge as
applying to the atom currently being examined; it never checks which endpoint bears the wedge. The
same gap applies to CXSMILES `w:`, `wU:`, and `wD:` entries, which supply an explicit atom endpoint
that the parser checks for incidence and then discards.

### Representation contract

The wedge is a depiction-level feature of the external boundary; the graph IR has no notion of a
wedge endpoint and gains none. TableIR keeps the endpoint pair-relative, as `BondDonation` already
is, so a wedge without an endpoint is unrepresentable and no separate incidence check exists.

- Type and role: `Bond` and `ExtendedBond` remain open external-boundary records. The `wedge`
  field becomes `Option<BondWedge>`, where `BondWedge` is a struct with public fields
  `orientation: BondOrientation` and `taper: BondTaper`.
- `BondOrientation` is the current five-variant enum under its new name: `Up` (MOL code 1),
  `Down` (code 6), `Either` (code 4, CXSMILES `w:`), `EitherUp` (`wU:`), `EitherDown` (`wD:`).
- `BondTaper::{Widening, Narrowing}` describes the wedge's width from the stored `AtomPair`
  first endpoint toward its second. `Widening` puts the pointed end at `first()`; `Narrowing`
  at `second()`. `BondTaper::flip` exchanges the variants.
- Accessor: `Bond::narrow_endpoint` and `ExtendedBond::narrow_endpoint` return
  `Option<u32>`, the pointed endpoint when a wedge is present.
- Construction: the existing constructors create no wedge; producers assign the field.
- Conversions: CTfile parsing derives the taper from source endpoint order before
  normalization. CXSMILES derives it from its explicit atom after the existing incidence
  check. `From<Bond>` and `TryFrom<ExtendedBond>` copy the value.
- Transformation: `update_atoms` receives the new ids of the old first and second endpoints;
  when their order reverses it flips the taper together with the donation and leaves the
  orientation unchanged.
- First consumer: raise reads a definite (`Up`/`Down`) wedge only at its narrow endpoint. The
  wide endpoint receives no tetrahedral assertion by incidence. `Either`, `EitherUp`, and
  `EitherDown` remain unread by raise, as at the baseline.
- Failure and absence: `None` means no wedge. The taper adds no error or panic boundary; the
  existing outgoing-wedge conflict error remains.
- Laws: two flips are identity; non-reversing `update_atoms` preserves the wedge; reversing
  `update_atoms` flips only the taper; basic/extended conversion preserves the wedge exactly.
- Field order: `Bond` and `ExtendedBond` group their fields as endpoints and bond kind, electronic
  state, stereo bond markers (`stereo`, `direction`), the stereo atom depiction (`wedge`), then
  source notation fields, instead of interleaving `wedge` between `stereo` and `direction`.
- Rust/Python boundary: no Python binding exposes these types. No graph-IR change.

Open before implementation: the taper variant names, whether `flip` is also offered on
`BondWedge`, whether a `wide_endpoint` accessor accompanies `narrow_endpoint`, and the exact
field order of the two bond records and of the atom records.

### Verification

Unit coverage: taper flip; `update_atoms` with and without endpoint reversal, asserting the whole
bond; `narrow_endpoint` for absent, widening, and narrowing wedges; CTfile bond lines with the
pointed end at the lower- and at the higher-numbered atom in both parsers; CXSMILES wiggly entries
naming either endpoint; conversion preservation; and raise cases for a wide endpoint with three or
four neighbors, two stereo sites sharing a wide endpoint, a site with one incoming and one outgoing
wedge, and a wedge whose pointed end is the higher-numbered atom.

Corpus evidence: the release-142 basic pipeline was rerun under
[wedge-endpoint](../scratch/rhea-census/rhea-142/wedge-endpoint/run-manifest.json) with
the same inputs, parser configurations, MDL counts-valence model, resolution configuration, and
aromaticity model as the census below.

| First outcome | Basic before | Basic after | Basic + editor before | Basic + editor after |
| --- | ---: | ---: | ---: | ---: |
| Resolved, concrete, and aromatized | 10,115 | 11,051 | 10,124 | 11,064 |
| Raise failure | 56 | 5 | 56 | 5 |
| Resolution contradictory | 937 | 52 | 944 | 55 |
| Resolution underdetermined | 378 | 378 | 378 | 378 |
| Parse failure | 1,465 | 1,465 | 1,449 | 1,449 |
| Missing source | 8 | 8 | 8 | 8 |

No previously accepted participant is lost. All 51 wedge-attributed raise conflicts and 885 of the
886 wide-endpoint contradictions are accepted; CHEBI:146007 moves from a spurious assertion at atom
4 to a contradiction at atom 6 (zero-based). The editor-enabled run additionally accepts its four
wide-endpoint contradictions. The remaining basic outcomes are the 43 valence exclusions, nine
tetrahedral contradictions (five narrow-endpoint sites, three parity sites on trigonal carbon,
CHEBI:146007), the four narrow-only wedge conflicts, and the two-ligand parity record.
Extended-parser acceptance is unchanged.
[Transitions](../scratch/rhea-census/rhea-142/wedge-endpoint/transitions.json) and the
per-configuration changed-outcome files retain the per-record evidence.

The 936 newly accepted records were probed directly
([raise-sites.jsonl](../scratch/rhea-census/rhea-142/wedge-endpoint/raise-sites.jsonl)).
In every record the raised tetrahedral sites are exactly the narrow endpoints of its Up/Down
wedges together with its parity-1/2 atoms, every raised coset is literal, resolution is
Determined, and the stereo atom count equals the site count. At the baseline these records had
carried assertions at 1,230 wide endpoints with three or four neighbors; 918 records contained at
least one, and the other 18 are mixed incoming/outgoing records whose only wide endpoint was
itself a site.

| Baseline category (basic) | Records | Outcome now |
| --- | ---: | --- |
| tetrahedral contradiction at a wide endpoint | 886 | 885 accepted; CHEBI:146007 fails at its narrow endpoint instead |
| wedge conflict at a wide endpoint only | 20 | accepted |
| wedge conflict at a narrow and a wide endpoint | 31 | accepted |
| wedge conflict at a narrow endpoint only | 4 | unchanged, raise conflict |
| tetrahedral contradiction at a narrow endpoint | 5 | unchanged |
| tetrahedral contradiction at a parity site | 4 | unchanged (3 contradictions, 1 raise ligand-count error) |
| transition-metal valence | 43 | unchanged |
| explicitly unknown stereo | 378 | unchanged |
| parse and missing source | 1,473 | unchanged |

CHEBI:40, CHEBI:15393, and CHEBI:58540 are catalogued as ChEBI parsing examples in the MOL parsing
suite. The umol-io unit, conformance, and property suites and clippy pass at the stage gate.

### CTfile stereo reading contract

The reader follows the published CTfile specification; the two additions below are stated as
such.

- Tetrahedral stereo comes from wedges only. An `Up` or `Down` wedge asserts `#T` with a literal
  coset at its narrow endpoint; an `Either` wedge (code 4), and the CXSMILES `w:`, `wU:`, `wD:`
  marks, assert `#T` with an undetermined coset at their narrow endpoint, following Appendix A ("a
  center with an Either bond has a parity value of 3"). Both readings require three or four
  ligands at the endpoint, as today. Atom parity is not read ("Ignored when read"). No tetrahedral
  configuration is derived from coordinates, in 2D or 3D; a record without wedges carries no `#T`.
- Double-bond stereo with code 0 is derived from the atom-block x, y, z coordinates for every
  double bond whose two ends each have two distinguishable ligands, the predicate
  `cis_trans_capable` already applies to SMILES. Code 3 asserts `#C` with an undetermined coset,
  as today. Coordinates are consulted only when the double bond's neighborhood carries no SMILES
  directional marks. Ring membership plays no part in the reader.
- The coordinate reading is supplementary and asserts nothing the drawing does not show. A
  code-0 double bond gets a `#C` only when the coordinates settle it: coincident bond atoms, a
  substituent on the bond axis, or both substituents of one atom on one side of it, as a drawn
  cage places them, produce no constraint, not an undetermined one and not a rejection. The
  reading is symmetric in the substituents: each atom's pair must lie on opposite sides, and one
  substituent of each atom then decides cis or trans, so the result does not depend on atom or
  bond numbering. A wedge whose endpoint ligands have all-zero or collinear projected positions
  raises `RaiseError::DegenerateWedgeGeometry { atom }`. The atom-block parsers keep all-zero
  coordinates; positions are absent only under `IGNORE_POSITIONS`.
- The parsers reject stereo codes the specification does not define for the bond order: codes 1,
  4, 6 on anything but a single bond, code 3 on anything but a double bond, and any nonzero code on
  other orders. Unconditional, in every preset.
- geometric-core gains `same_side_of_axis(axis_start, axis_end, first, second) -> Option<bool>`:
  whether two points lie on the same side of the line through the axis, using their
  components perpendicular to it, `None` when either component is below the tolerance constant
  `AXIS_SIDE_TOLERANCE`. One formula serves 2D and 3D.
- The stereo model, not the reader, decides which ring double bonds are stereo bonds.
  `StereoModel` gains `stereo_bond_minimum_ring_size: u32`, default 8, the smallest
  ring in which a trans double bond is isolable. Stereo perception drops a `#C` assertion on a
  double bond whose smallest ring, from `RingSet::bond_smallest_ring_size` under the default ring
  model, is below the threshold: no stereo bond is created and the assertion is cleared to
  undetermined in the edit plan, the same edit `reset_stereo_constraints` uses. It is not a
  contradiction. The parameter is bound in Python with the rest of `StereoModel`.
- Addition to the specification: an `Either` wedge whose narrow endpoint is an atom of exactly
  one double bond is read as an unknown configuration of that double bond, `#C` with an
  undetermined coset, the reading code 3 gives. The specification defines code 4 for tetrahedral
  centers only; the census has 44 records drawn with the wavy-bond convention for unknown E/Z and
  none of them carries code 3 on the double bond. No `#T` is raised at such an endpoint.
- Addition to the specification: the coordinate reading produces no `#C` for a double bond whose
  end carries only a further double bond, the cumulated ends of allenes and heterocumulenes,
  where the substituents of the two ends are not coplanar. The SMILES directional reading and
  `cis_trans_capable` are unchanged.
- CXSMILES `c:` and `t:` are read only once their ligand frame is settled from ChemAxon's
  published definition; the vendored CDK reader ignores them and the RDKit reader applies a
  "CX ordering" the repository does not document. Until then they stay unread.
- TableIR supplies the neighbor table the readings consult. `Molecule::atom_neighbors()` and
  `ExtendedMolecule::atom_neighbors()` build an `AtomNeighbors` in one pass over the bond block:
  per atom, its `Neighbor { atom, bond }` entries in bond order, a duplicate bond listed twice, a
  self-bond once, an endpoint outside the record skipped; `neighbors(atom)` and `degree(atom)`
  (distinct atoms) read it. The raise builds the table once per record and passes it to every
  reading; nothing rescans the bond list per atom or per bond. The table is built on request
  because the records are open to mutation.

Effects measured on the repository catalogs: two catalog records (`openbabel/culgi_10.mol`,
`indigo/bold-bonds.mol`) carry a wedge code on a double bond and move to the invalid category; the
parser tests that expect wedge codes on double or triple bonds become rejection rows. The census
effects are under Stereo reading verification below; 3D records without wedges and the 23
reversed-wedge sites carry no tetrahedral stereo, both future lint items.

### Stereo reading verification

The release-142 basic and basic+editor pipelines were rerun under
[stereo-reading](../scratch/rhea-census/rhea-142/stereo-reading/) with the inputs and models of
the census; "before" is the wedge-endpoint run above.

| First outcome | Basic before | Basic after | Basic + editor before | Basic + editor after |
| --- | ---: | ---: | ---: | ---: |
| Resolved, concrete, and aromatized | 11,051 | 11,367 | 11,064 | 11,380 |
| Raise failure | 5 | 4 | 5 | 4 |
| Resolution contradictory | 52 | 49 | 55 | 52 |
| Resolution underdetermined | 378 | 66 | 378 | 66 |
| Parse failure | 1,465 | 1,465 | 1,449 | 1,449 |
| Missing source | 8 | 8 | 8 | 8 |

Outcome transitions in the basic configuration, identical in the editor configuration
([comparison](../scratch/rhea-census/rhea-142/stereo-reading/comparison.txt),
[per-record diagnostics](../scratch/rhea-census/rhea-142/stereo-reading/new-failures.txt)):

| Transition | Records | Cause |
| --- | ---: | --- |
| underdetermined to accepted | 357 | atom parity 3 no longer read; the six other parity-3 records carry an `Either` wedge and stay Underdetermined |
| contradictory to accepted | 3 | parity on a trigonal carbon no longer read (CHEBI 34463, 79808, 79866) |
| raise failure to accepted | 1 | parity on a two-neighbor carbon no longer read (CHEBI:80541) |
| accepted to underdetermined | 45 | 41 `Either` wedges at a double-bond atom, now `#C` undetermined; 4 `Either` wedges at a tetrahedral atom, now `#T` undetermined (CHEBI 17901, 71494, 141215, 144485) |

The 66 Underdetermined records are the 15 with bond stereo code 3, 43 further records with an
`Either` wedge at a double-bond atom, and 8 with an `Either` wedge at a tetrahedral atom (the four
above and CHEBI 25164, 57513, 57850, 57878, whose parity-3 flags sat on the same sites). No
record newly fails as contradictory or at the raise: coordinate-derived `#C` produced no
inconsistency anywhere, and the drawings that settle nothing, hyponitrous acid drawn as a line
(CHEBI:14428) and three ring double bonds with a substituent on the axis (CHEBI 63464, 63466,
234435), produce no constraint.
The 21 records that lose parity-only stereo atoms are those identified under MOL atom parity
reading. The remaining stereo rejections are unchanged: the four narrow-endpoint wedge conflicts
(CHEBI 34596, 35306, 219836, 231826), and the six tetrahedral contradictions at allene termini,
biaryl axes, and double-bonded ring carbons (CHEBI 32446, 63663, 145950, 146007, 231497, 231926);
the editor configuration adds its three taxane wedges.

The all-record probe ([raise-sites-all.jsonl](../scratch/rhea-census/rhea-142/stereo-reading/raise-sites-all.jsonl))
raised 31,504 `#C` assertions in the 11,367 accepted records, all literal; 3,261 records have at
least one stereo bond after resolution, 9,595 stereo bonds in all. RDKit
2025.03.2 read the same files as independent evidence
([arbitration](../scratch/rhea-census/rhea-142/stereo-reading/rdkit-cis-trans-arbitration.txt)),
with both readings transported to the lowest-indexed substituent at each end:

| Raised `#C` assertions | Count | Reading |
| --- | ---: | --- |
| stereo bond in umol and RDKit, same configuration | 8,928 | agreement at every shared bond; no disagreement |
| stereo bond in umol only, chain bond | 628 | two ligands of one end are constitutionally equivalent at 626 (610 with one substituent at the other end); one imine (CHEBI:16825) and one azo bond (CHEBI:17903) RDKit leaves unassigned |
| stereo bond in umol only, aromatic macrocycle bond | 39 | meso bonds of ten porphyrins and chlorins, smallest ring 16, Kekulé-drawn; stereo resolution precedes aromaticity perception |
| no stereo bond, ring below 8 | 21,909 | 19,353 Kekulé bonds of aromatic five- to seven-membered rings, 2,552 double bonds in five- to seven-membered rings, 4 heme meso bonds in a six-membered ring through Fe (CHEBI:61715) |
| stereo bond in RDKit only | 0 | |

The 9,595 stereo bonds are the sum of the first three rows; the 21,909 skipped assertions are
cleared by the stereo model's ring threshold. The zero-disagreement result rests on RDKit's
coordinate reading being independent of umol's, which the parity finding above already used.
The `stereo_bond_sites` field of the probe lists the stereo bonds per record; the
`cis_trans_arbitration.py` and `stereo_reading_diagnostics.py` scripts of the research runner
reproduce the comparison and the diagnostics file.

### Implementation plan

All work is in `umol-io`; no other crate names the wedge types. S0 is one breaking stage because
the field type change invalidates every producer and the one consumer; each subitem carries its
tests and the tree is green at the end of the stage.

S0 — TableIR wedge representation

- S0a `table_ir::bond`: `BondTaper` with `flip`. Tests: `test_bond_taper_flip`. Additive. Completed.
- S0b `table_ir::bond`: rename the orientation enum to `BondOrientation`; add the `BondWedge`
  struct with `flip`; change the `wedge` field of `Bond` and `ExtendedBond` to `Option<BondWedge>`;
  reorder both records and their `From`/`TryFrom` impls to the settled field order; `update_atoms`
  flips the wedge with `BondWedge::flip` beside the donation; add `narrow_endpoint` and
  `wide_endpoint` to both records. Tests: `test_bond_wedge_flip`, `update_atoms` tables asserting
  whole bonds with identity tables split out, `narrow_endpoint` and `wide_endpoint` tables, and
  conversion rows carrying a wedge. Breaking [dep: S0a]. Completed.
- S0c `ctfile::parser::convert`, `ctfile::parser::bond`: the stereo/direction converter returns
  `Option<BondOrientation>`; both bond line parsers derive the taper from source endpoint order
  and assign a `BondWedge`. Tests: the converter table and parser rows with the pointed end at the
  higher-numbered atom. Migration [dep: S0b]. Completed.
- S0d `smiles::parser::cx`: `CxEntry::WigglyBonds` carries `BondOrientation`; both apply sites
  derive the taper from the explicit atom after the incidence check. Tests: wiggly rows naming
  either endpoint. Migration [dep: S0b]. Completed.
- S0e `table_ir::raise::utils`: `wedge_bond_neighbors` keeps only wedges whose narrow endpoint
  is the examined atom and reads the orientation; the stale exclusion comment in
  `raise_tetrahedral_stereo` is removed. Tests: reduced MOL reproducers in the existing
  tetrahedral tables for a wide endpoint with three or four neighbors, two sites sharing a wide
  endpoint, one incoming and one outgoing wedge at a site, and a wedge pointed at the
  higher-numbered atom. Migration [dep: S0b]. Completed.

Gate: `cargo test -p umol-io --features conformance,proptest` and
`cargo clippy -p umol-io --all-targets --features conformance,proptest -- -D warnings`.

S1 — Parsing examples

- S1a `tests/mol_parsing`: add CHEBI_40, CHEBI_15393, and CHEBI_58540 to `data_raw/chebi/`,
  classify them with `classify_mol_files` (all three are `molecule`), place the copies under
  `data/molecule/chebi/`, record their snapshots, and remove `tests/raise/`. Additive. Completed.

Gate: `cargo test -p umol-io --features conformance --test mol_parsing`.

S2 — Corpus verification

- S2a Rebuild the census runner against the workspace, rerun the basic and basic+editor
  configurations on the release-142 inputs, and record the outcome table under Verification
  above [dep: S0]. Completed.

Critical path: S0a → S0b → S0c, S0d, S0e → S2a. S1 is independent of S0.

Reopened for the CTfile stereo reading contract. Each stage ends green; breaking subitems carry
their test migration.

S3 — Parser stereo codes by bond order

- S3a `ctfile::parser::bond`: both bond line parsers reject a stereo code not defined for the
  parsed bond order at column 9. Tests: the rows `len21_double_wedge_up`, `len21_double_wedge_down`,
  `len21_double_wedge_either`, `len12_double_wedge_up`, `len21_triple_ignored_stereo`,
  `len21_triple_empty_fields` (if its code is nonzero), `len_12_double_wedge_up`, `len_21`,
  `len_21_reaction_center`, and the `len21`/`len12` rows of both block tables become error rows
  with the column; new rows cover code 3 on a single bond. Breaking [dep: none]. Completed.
- S3b `tests/mol_parsing`: reclassify `openbabel/culgi_10.mol` and `indigo/bold-bonds.mol` into
  `data/invalid/`, record their snapshots. [dep: S3a]. Completed.

Gate: `cargo test -p umol-io --features conformance,proptest`.

S4 — Parity removal

- S4a `table_ir::raise`, `table_ir::raise::utils`: `raise_tetrahedral_stereo` ignores
  `chirality` when the frame is `LastNeighborAway`; `last_neighbor_away_ordering` is removed;
  the `ChiralityFrame::LastNeighborAway` doc comment states the field is retained as parsed and
  not read. Tests: delete `mol_parity_clockwise` and `CHIRAL_PARITY_MOL`; point
  `test_parse_mol_to_ir_stereo::tetrahedral` at `CFCLBRI_SINGLE_WEDGE_MOL` site 1. Breaking
  [dep: none]. Completed.

S5 — Either wedges

- S5a `table_ir::raise::utils`, `table_ir::raise`: `Either`, `EitherUp`, and `EitherDown` at
  their narrow endpoint yield `#T` with an undetermined coset after the ligand-count check; a
  definite and an either wedge at one endpoint is a wedge conflict. Tests: MOL rows with code 4,
  CXSMILES rows with `w:`, `wU:`, `wD:`, a two-ligand code-4 error row, a mixed-wedge conflict
  row. The predicate is `has_either_wedge` in the raise utilities. Additive [dep: S4a].
  Completed.

S6 — Geometry primitive

- S6a `umol-geometric-core`: `same_side_of_axis` and `AXIS_SIDE_TOLERANCE`. Tests: same side,
  opposite side, in 2D and 3D, degenerate axis, degenerate substituent, tolerance boundary.
  Additive [dep: none]. Completed.

S7 — Stereo-model ring threshold

- S7a `umol-graph::ops::model`: `StereoModel::stereo_bond_minimum_ring_size`, default 8; Python
  binding of the field. Tests: default, construction, Python parity. Additive [dep: none].
  Completed.
- S7b `umol-graph::ops::stereo`, `ops::resolve::stereo`: perception consults the ring set and
  skips `#C` on ring double bonds below the threshold; the resolver clears those assertions.
  Tests: cyclohexene with `#C` yields no stereo bond and a cleared constraint; cyclooctene yields a
  stereo bond; threshold 0 keeps cyclohexene; a resolution-suite fixture for each. Additive
  [dep: S7a]. Completed. The perception reports the skipped assertions in
  `StereoDerivation::skipped_stereo_bonds`; the fixture is
  `stereo_cis_trans/cyclohexene-asserted.edn`.

Gate: `cargo test -p umol-graph --features conformance` and the Python suite.

S8 — Cis/trans from coordinates

- S8a `table_ir::raise::utils`: a positions-based counterpart of `cis_trans_side` building the
  same `StereoBondAtom`, chosen when the bond's neighborhood has no directional marks and the
  molecule has positions; degenerate geometry raises the new errors; the wedge path raises
  `DegenerateWedgeGeometry` on all-zero or collinear positions instead of computing. Tests:
  fumarate and maleate literals in 2D and in 3D, a ring double bond asserted as the spec says,
  code 3 unchanged, all-zero coordinates rejected, a collinear substituent rejected, SMILES
  directional rows unchanged. Breaking for the raise's MOL output [dep: S6a, S7b]. Completed;
  the parsers' former collapse of all-zero coordinates to absent positions is removed. A drawing
  that does not settle the configuration yields no constraint, with rows for zero coordinates, a
  substituent on the axis, the same drawing with a reversed bond line, and a folded projection.
  The ring gate of S7 now also
  applies to the resolver's partial-constraint check, through
  `StereoPerception::skipped_stereo_bonds`, so an undetermined assertion on a small-ring double
  bond is cleared instead of making the whole resolution Underdetermined; a cyclohexene row with
  an undetermined coset covers it in perception and resolver tests.

Gate: `cargo test -p umol-io --features conformance,proptest` and `cargo test -p umol-graph`.

S9 — Corpus verification

- S9a Rerun the census basic and basic+editor configurations and the all-record probe; record
  the outcome table, the stereo bond counts, and any new contradictions from coordinate-derived
  `#C` under Verification [dep: S8a]. Completed: Stereo reading verification above; no new
  contradiction; RDKit agrees at every shared stereo bond.

S10 — Covalently bonded transition-metal states

- S10a `umol-graph/config/default-registry.toml`: rows for the formal-charge encodings in the
  census, nonbonding electrons from valence electrons minus charge minus covalence; one spin
  state per encoding for now, Fe(III) high-spin and Co(II) low-spin. Each row's lone pairs,
  unpaired electrons, and multiplicity are to be checked against Holleman-Wiberg before commit.

  ```toml
  [Fe]   # Fe(III), d5, high-spin only for now: nonbonding 5 = #u5 (sextet)
  -3 = ["Fe #c-3 #v6 #u5"]
  -2 = ["Fe #c-2 #v5 #u5"]
  -1 = ["Fe #c- #v4 #u5"]
   0 = ["Fe #v3 #u5"]          # added to the existing key
   1 = ["Fe #c+ #v2 #u5"]
  [Co]   # Co(II), d7, low-spin only for now: nonbonding 7 = #n3 #u (doublet)
  -3 = ["Co #c-3 #v5 #n3 #u"]
  -2 = ["Co #c-2 #v4 #n3 #u"]
   1 = ["Co #c+ #v1 #n3 #u"]
  [Mn]   # Mn(IV), d3: nonbonding 3 = #u3 (quartet)
   0 = ["Mn #v4 #u3"]                                                     # added to the existing key
  [V]    # V(IV) vanadyl, d1: nonbonding 1 = #u (doublet)
   2 = ["V #c+2 #v2 #u"]                                                  # added to the existing key
  ```

  Tests: registry parse and resolution-suite fixtures under the atom-typing model. Additive
  [dep: none]. Completed: seven fixtures under `resolution/data/transition_metals/`, real
  compounds with explicit ligand hydrogen counts so that only the metal state is open:
  hexafluoridoferrate(3-), tetrachloridoferrate(1-), dihydroxidoiron(1+), pentacyanocobaltate(3-),
  cobalt(II) porphine, manganese dioxide, vanadyl. Each resolves to the intended metal state in
  all four cells, atom-typing and counts, most-saturated and strict, with no tie-breaks. The
  Fe #c-2 #v5, Fe #v3, and Co #c+ #v1 rows have no fixture. Cobalt porphine is Kekulé-written
  with two N+ ring atoms, the way ChEBI draws metalloporphyrins; the aromatic-written input of
  the same molecule is contradictory (open item 12). The frozen MDL counts table the census
  selects is unchanged.
- The counts model is untouched: the MDL table is frozen by its header, and the counts model's
  isoelectronic charge shift has no reading for a formal-charge metal encoding.

S11 — Line-level parser errors are cuts

- S11a `ctfile::parser::{atom, bond, counts, legacy_atom_list, properties, utils}`: a field or
  line error inside a fixed-width record is `ErrMode::Cut` carrying the column from the point of
  failure; the block parsers stop converting backtracks into cuts. `Backtrack` remains only where
  an alternative exists: block dispatch in `parser.rs`, the `alt` sites in `rgroup` and `sgroup`,
  and the `opt` tails of multi-line properties. Today 73 line and field sites backtrack, which no
  caller can use, because a fixed-width line admits no alternative parse. Tests: the existing
  error tables keep their columns; new rows show that a line error is reported at its column
  without a later alternative being tried. Breaking [dep: none]; independent of S3 to S10.

Deferrable: CXSMILES `c:`/`t:` reading, gated on the ChemAxon definition.

Critical path: S6a → S7a → S7b → S8a → S9a; S3, S4, S5 are independent of it and of each other
except S5 after S4.

## Bond stereo code converter

An earlier revision of this document reported that bond stereo code 2 reached the `unreachable!`
arm of `convert_bond_stereo_direction_code`. At the current parser both bond line parsers validate
the field against `0 | 1 | 3 | 4 | 6` and reject any other value at column 9 before conversion,
so that arm is not reachable from parsed bytes. No change is required. The archive contains no
bond stereo code outside 0, 1, 3, 4, and 6 in any of its raw V2000 records.

## Release-142 census

Rhea 142, released 2026-09-02, is the current baseline. This rerun uses individual MOL files, not
SMILES. The release downloads live under `materials/databases/elixir/rhea-142`; the standalone
research runner and every analysis output live under
[scratch/rhea-census](../scratch/rhea-census/README.md). Both directories are ignored by git, and
the runner changes no production behavior. The upstream
[CTfile README](https://ftp.expasy.org/databases/rhea/ctfiles/README.txt) identifies chebi.sdf.gz
as the ChEBI snapshot synchronized to Rhea. SHA-256 fingerprints and download URLs are in
[sources.json](../materials/databases/elixir/rhea-142/sources.json).

| Download | Bytes |
| --- | ---: |
| rhea-mol.tar.gz | 5,173,355 |
| chebiId_name.tsv | 589,500 |
| chebi.sdf.gz | 69,543,987 |

There are 12,959 distinct listed participants. The archive holds 14,058 ChEBI MOL files and 264
polymer MOL files. Of the listed participants, 12,951 join to ChEBI MOL files and the same eight
ids as above are missing. Every joined ChEBI MOL block matches its Rhea MOL after excluding only
the title line and line endings, and the synchronized SDF supplies the same 12,951 participants.

The runner uses source revision cd902742b20eec4cf4aac65c0d0e5afa8e93fee9 and applies
`ValenceModel::mdl()`, MostSaturated selection, Daylight aromaticity, the default stereo model,
and `ResolveConfig::default()`. Each participant has one outcome per configuration in
[census.jsonl](../scratch/rhea-census/rhea-142/census.jsonl).

### Parser acceptance

| Parser configuration | Parsed MOL files | Parse failures | Missing files |
| --- | ---: | ---: | ---: |
| Basic | 11,486 | 1,465 | 8 |
| Extended | 12,932 | 19 | 8 |
| Basic + editor | 11,502 | 1,449 | 8 |
| Extended + editor | 12,948 | 3 | 8 |
| Extended + editor + pseudoatoms | 12,949 | 2 | 8 |

Every basic-parse success also parses with the extended parser. The 1,446 additional records
accepted by the default extended parser contain these parsed feature combinations:

| Additional input features | Records | Examples |
| --- | ---: | --- |
| R groups only | 1,437 | CHEBI:685 |
| Heavy-atom wildcard only | 2 | CHEBI:13193, CHEBI:17499 |
| Halogen wildcard only | 3 | CHEBI:16042, CHEBI:70856, CHEBI:85638 |
| R groups and halogen wildcards | 4 | CHEBI:17792, CHEBI:18060, CHEBI:137405, CHEBI:137406 |

The 19 extended parse failures are the 16 editor properties plus CHEBI:10545 (electron),
CHEBI:30212 (photon), and CHEBI:57503 (an R-group anthocyanidin glycoside). Enabling editor
properties accepts the same 16 files in either parser; nine of them resolve and aromatize, and
seven expose tetrahedral contradictions, four naming only a wide wedge endpoint and three naming
a narrow endpoint. Enabling pseudoatoms admits CHEBI:10545 as `Pseudoatom("e")` with charge -1.
CHEBI:30212 and CHEBI:57503 then fail at column 69, where a trailing zero remains beyond the
69-character atom fields. Offsets here and in the ledger are zero-based.

The extended inputs were inventoried at atom and attachment level in
[semantic-inventory.json](../scratch/rhea-census/rhea-142/semantic-inventory.json) and
[semantic-features.jsonl](../scratch/rhea-census/rhea-142/semantic-features.jsonl):

| Observed feature | Scope in the 1,446 records |
| --- | --- |
| R-group sites | 2,550 sites in 1,441 records |
| R-group attachment degree | 2,499 degree-one sites; 51 degree-two sites in 51 records |
| R-group labels | 1,993 unlabeled sites; 557 labeled sites, labels 1 through 9 |
| Repeated numbered label in one record | CHEBI:194022 has two sites labeled 1 |
| Heavy-atom wildcard | Two sites: isolated A and A bonded to two explicit H atoms |
| Halogen wildcard | Seven sites: one isolated anion and six degree-one sites |
| Explicit extended atom fields beyond the symbols | Charge -1 on CHEBI:16042's halogen wildcard |
| R-group label properties | `M  RGP` in 272 records |
| Explicit R-group logic or attachment-order properties | No `M  LOG`, `M  APO`, or `M  AAL` |
| Other extended payload | No extended bond features or molecule-level S-group/R-group records |

CHEBI:28965 and CHEBI:76307 have degree-two R sites. Parsed `RGroup` values contain 2,039
occurrence lists equal to `[GreaterThan(0)]` and 511 empty lists: `RGroup::new` supplies the first
form and applying `M  RGP` replaces the atom symbol with a minimal `RGroup` and an empty list.
These are producer-path differences, not two source assertions. Forty-eight of the extended
records have stereo-coded bonds incident to an extended site, five contain explicit unknown stereo
markers, and two contain transition metals (Hg in CHEBI:83725, Co in CHEBI:140785).

### Basic pipeline outcomes

| First outcome | Basic | Basic + editor |
| --- | ---: | ---: |
| Missing source | 8 | 8 |
| Parse failure | 1,465 | 1,449 |
| Raise failure | 56 | 56 |
| Resolution underdetermined | 378 | 378 |
| Resolution contradictory | 937 | 944 |
| Resolved, concrete, and aromatized | 10,115 | 10,124 |
| **Total participants** | **12,959** | **12,959** |

There were no parser panics, resolution execution failures, determined-but-non-concrete results,
or aromatization failures. The basic first-symbol breakdown is R 985, R# 417, R1 39, X 4, A 2,
hv 1, and e 1, plus 16 records stopping on `M  ZZC`. Relative to release 141 the basic run has
eight more accepted participants and one more atom-line rejection; the other stage totals are
unchanged.

The stereo and valence classifications were recomputed from the release-142 raw atom and bond
records and reproduce the release-141 breakdown exactly: 20 wide-only and 31 mixed wedge raise
conflicts, four narrow-only conflicts, one two-ligand parity failure, 886 wide-only tetrahedral
contradictions, eight other tetrahedral contradictions, 363/12/3 unknown-stereo records, and the
same six-metal valence table.
[Classification counts](../scratch/rhea-census/rhea-142/classification-counts.json),
[category members](../scratch/rhea-census/rhea-142/category-members.json), and the
[per-participant evidence](../scratch/rhea-census/rhea-142/classification.jsonl) preserve
exact ids, source fingerprints, first-failure lines, and wedge endpoint relationships.

### Resolver behavior on partial input

Both `CountsValence::admit` and `AtomTypingValence::admit` first scan every atom and return
`Underdetermined` with an empty completion map if any element is non-literal, including finite
and complement sets. A standalone
[partial-IR probe](../scratch/rhea-census/census/src/bin/partial-probe.rs) constructs IR
directly; its [results](../scratch/rhea-census/rhea-142/partial-probe.jsonl) establish:

- An isolated neutral carbon with open hydrogen count completes under the MDL counts model.
- Adding a disconnected unspecified atom, bonding the carbon to an unspecified atom, or adding a
  disconnected finite element set makes both valence admission engines return an empty
  underdetermined result, and composite resolution leaves the molecule unchanged.
- The default atom-typing model retains plural carbon candidates on the isolated control under
  its Strict policy; the non-literal element suppresses even that candidate report.
- `AtomView::total_hydrogens` treats every non-literal neighbor element as potentially H,
  including a finite set that excludes H.

The determinacy gates of the resolver entry points:

| Resolver / entry point | Specific dependency | Current consequence |
| --- | --- | --- |
| [ValenceResolver::admit](../umol-graph/src/ops/resolve/valence.rs) | Valence, donated-pair, and accepted-pair incidence assertions are checked with DerivedComplete. | Any Underdetermined check stops admission for the whole molecule. |
| [CountsValence::admit](../umol-graph/src/ops/valence/counts.rs) | Every atom must have a literal element. | One non-literal element returns an empty underdetermined completion map for the whole molecule. |
| CountsValence::admitted_completions / candidate_states | Local admission skips undetermined charge and non-literal incident valence. Enumeration assumes literal element and charge. CountsInput reads accepted pairs with unwrap_or(0). | An atom may be skipped even after global admission succeeds. |
| [AtomTypingValence::admit](../umol-graph/src/ops/valence/atom_typing.rs) | Every atom must have a literal element. | Same global stop as counts. Local registry lookup permits non-literal charge. |
| [AromaticityResolver::select](../umol-graph/src/ops/resolve/aromaticity.rs) | A non-literal positive aromatic contribution outside the candidate carrier stops selection. | Whole-selection Underdetermined. MAX_ASSIGNMENTS is a separate search-limit exit. |
| AromaticityResolver::plan / [AromaticityPerception::derive](../umol-graph/src/ops/aromaticity.rs) | Standalone derivation checks for non-literal positive aromatic-valence assertions and undetermined stored electron counts. | Whole-derivation Underdetermined; standalone resolution produces no edits. |
| [HmoAromaticity::build_calculator](../umol-graph/src/ops/aromaticity/hmo.rs) | Participating atoms require literal elements and supplied electron contributions. | UndeterminedAtom propagates to an underdetermined aromaticity outcome. |
| [StereoResolver::plan](../umol-graph/src/ops/resolve/stereo.rs) | Any tetrahedral or cis/trans assertion that is non-ground and non-vacuous. | Stops all stereo work before frame derivation, returning an empty underdetermined edit plan. |
| [StereoPerception::derive_stereo_atom / derive_stereo_bond](../umol-graph/src/ops/stereo.rs) | Literal site elements; when a virtual ligand is needed, literal implicit-H counts and possibly lone-pair counts. | Missing local values return None, which perception classifies as a stereo failure. |
| [BondsResolver::plan](../umol-graph/src/ops/resolve/bonds.rs) | No determinacy gate. | Fills charge and unpaired-electron defaults while other fields can remain partial. |
| [MulticenterBondsResolver::plan](../umol-graph/src/ops/resolve/multicenter.rs) | Multicenter-valence incidence assertions checked with DerivedComplete. | Any Underdetermined check stops defaults for all multicenter bonds. |

[Resolver::resolve](../umol-graph/src/ops/resolve.rs) adds three publication barriers: unresolved
atom candidates stop processing before constitution edits are applied; later underdetermined
phases stop the operation; and the final `is_concrete` check controls assignment back to the
caller. All such exits leave the caller's molecule unchanged, including when an internal working
copy has completed valence; `test_resolver_resolve_later_underdetermined` covers this behavior.
[ingest::interpret_molecule](../umol-graph/src/ingest.rs) turns `Underdetermined` into an ingestion
error.

### MOL atom parity reading

The CTfile specification marks the V2000 atom parity field "Ignored when read". The raise reads it
nevertheless, with precedence over wedges, into a tetrahedral constraint in the LastNeighborAway
frame. A scratch build in which the raise ignores parity under that frame was compared with the
current tree on every accepted record; the evidence is under
[parity-unread](../scratch/rhea-census/rhea-142/parity-unread/).

- At 15,404 sites carrying both a parity flag and an outgoing wedge, the parity reading and the
  wedge reading give different cosets at 15,033. Of the 371 agreements, 348 are sites with an
  explicit hydrogen neighbor that is not the highest-indexed neighbor; the MDL definition numbers a
  hydrogen ligand last, the raise orders by index only, and the two errors cancel there
  ([cross-tab](../scratch/rhea-census/rhea-142/parity-unread/parity-vs-wedge-by-site.txt),
  [hydrogen rule](../scratch/rhea-census/rhea-142/parity-unread/parity-vs-wedge-crosstab.txt)).
- RDKit, reading only wedges, reproduces umol's expected cosets on umol's own four-neighbor wedge
  fixtures and then agrees with umol's wedge reading at all 750 four-neighbor parity sites
  examined, while disagreeing with the parity reading at 634 of them
  ([arbitration](../scratch/rhea-census/rhea-142/parity-unread/rdkit-arbitration.txt)).
- From the definition: with ligand 4 behind and ligands 1, 2, 3 counterclockwise (parity 2), the
  configuration is the one SMILES writes as `@` over the ascending ligand order, which is umol
  coset 0. The raise maps parity 1 to coset 0, and the `ChiralityFrame` doc comment states the
  same inverted correspondence. The `CHIRAL_PARITY_MOL` fixture encodes it.

The current tree therefore derives an inverted configuration from every parity flag that is not
accompanied by an out-of-order explicit hydrogen. Because the resolver sees a self-consistent
assertion either way, no census outcome exposed this. It does not explain the CHEBI:40
cross-format discrepancy: with parity unread both routes still produce four stereo atoms and
remain canonically unequal.

Decision: the raise stops reading atom parity, following the specification (S4). Measured effect of
the scratch build in the basic configuration, with `Either` wedges still unread: 363 unknown-stereo
records, the three parity-on-trigonal-carbon contradictions, and the two-neighbor parity record
become accepted; no record newly fails; the 15 records with bond stereo code 3 remain
Underdetermined. The two unit tests that encoded parity reading,
`test_raise_tetrahedral_stereo::mol_parity_clockwise` and `test_parse_mol_to_ir_stereo::tetrahedral`
with the `CHIRAL_PARITY_MOL` fixture, were removed; no conformance suite, umol-graph test, or Python
test depended on it. Twenty-one accepted records lose 23 stereo atoms: at each of those sites the
only wedge is drawn from a terminal substituent toward the center, so the narrow end carries no site
and the parity flag was the only readable marker. These reversed wedges are a first candidate for
boundary linting. The final measurement, with `Either` wedges read, is under Stereo reading
verification.

### Paired MOL and SMILES input for CHEBI:40

With the wedge defect removed, and equally with parity unread, both CHEBI:40 inputs resolve
concretely with four stereo atoms, but their complete canonical IR values are not equal. The MOL
route retains Kekulé double bonds alongside aromatic systems; the SMILES route has single localized
bonds in those systems. `AromaticityPerceiver::add_systems` adds aromatic constraints without
changing localized bond orders. Removing atom/bond constraints and setting aromatic-system localized
bonds to order one in a diagnostic copy does not eliminate the discrepancy. The `wedge-pair-probe`
bin of the research runner reproduces the comparison; the cause is open.

### Remaining non-accepted records

Basic configuration after the stereo reading plan, 1,592 of 12,959 participants:

| Group | Records | Direct cause |
| --- | ---: | --- |
| Missing source | 8 | participant-list entry without a MOL file |
| Extended atom symbols | 1,449 | R 985, R# 417, R1 39, X 4, A 2, hv 1, e 1; rejected by the basic atom-line parser |
| Editor property | 16 | `M  ZZC` outside `EDITOR_EXTENSIONS` |
| Explicitly unknown stereo | 66 | bond stereo code 3 (15), `Either` wedge at a double-bond atom (43), `Either` wedge at a tetrahedral atom (8); raised as undetermined cosets, so the stereo resolver returns Underdetermined |
| Transition-metal valence | 43 | no counts row for Fe, Co, Mn, V, Cu, or Cr under the frozen MDL model; the atom-typing registry now holds the covalent states (S10) |
| Axial stereo drawn with wedges | 5 | allene termini (CHEBI 32446, 34596, 35306) and biaryl axes (CHEBI 145950, 146007): a wedge whose narrow end is a trigonal carbon |
| Wedge at a double-bonded ring carbon | 3 | CHEBI 63663, 231497, 231926: one Up wedge at a carbon with a ring double bond, three neighbors, valence four |
| Two outgoing wedges forming an inconsistent projection | 2 | CHEBI 219836, 231826: each wedge alone reads a definite coset; the two differ |
| **Total** | **1,592** | |

The editor-enabled configuration adds three more wedges at a double-bonded ring carbon
(CHEBI 11302, 15208, 50436, all taxanes) and clears four wide-endpoint contradictions.

Stereo causes in detail:

- Allene termini and biaryl axes: the wedge's narrow endpoint has three neighbors and localized
  valence four. In CHEBI 34596 and 35306 the terminus carries one Up and one Down wedge that
  disagree under the tetrahedral reading, so raise reports a conflict; in CHEBI 32446 the two
  readings coincide, raise succeeds, and the stereo resolver finds no fourth ligand. The biaryl
  wedge sits on the axis bond with an aromatic carbon as narrow end; for CHEBI:146007 this is the
  site that previously received a spurious assertion from the wide end. umol has no axial stereo
  class, so these records cannot be admitted faithfully without one.
- Double-bonded ring carbons: the sole wedge is on a substituent of an alkene carbon in a
  bridged or medium ring. The tetrahedral reading has no fourth ligand.
- Inconsistent projections: in a projected tetrahedron with two plain in-plane bonds, both
  out-of-plane bonds must project into the exterior angle of the plain pair. In CHEBI 219836 (site
  C6, neighbors at -102°, 2°, 49°, 148°, wedges to 2° and 148°) and CHEBI 231826 (site C0,
  neighbors at -157°, -90°, 20°, 150°, wedges to -157° and 20°) one wedged bond projects into the
  interior angle, so no tetrahedron satisfies both wedges. Each wedge alone is readable; umol
  reports the disagreement rather than selecting one.
- Unknown stereo: the 15 code-3 records and the `Either` wedges are read as the source states
  them, undetermined cosets. They are accepted only once the resolver admits undetermined
  stereo input (open item 5).

Transition-metal encodings among the 43 valence records (45 metal atoms):

| Encoding | Metal atoms |
| --- | ---: |
| Co, charge -2, four single bonds | 15 |
| Fe, charge -3, six single bonds | 10 |
| Fe, charge -1, four single bonds | 6 |
| Fe, charge -2, five single bonds | 2 |
| Fe, charge +1, two single bonds | 2 |
| Co, charge -3, five single bonds | 1 |
| isolated ions Fe3+, Mn2+, Cu2+, Co2+, Cr3+ | 5 |
| V, charge +2, one double bond | 1 |
| Co, charge +1, one single bond | 1 |
| Fe, charge 0, one single and one double bond | 1 |
| Mn, charge 0, two double bonds | 1 |

The census selects the frozen MDL counts table, which has no row for these encodings. Under the
atom-typing model the S10 registry rows resolve the `resolution/data/transition_metals/` fixtures;
the census configuration was not changed, so the 43 records stay in this table.

## Open work

1. Resolver diagnostics: identify the underdetermined phase and cause in `ResolveReport`, and
   retain atom-local context for counts-valence mismatch, preserving the distinction between
   semantic `Solution` outcomes and operational `Result` errors.
2. The CTfile stereo reading contract is implemented (S3 to S9); S11, line-level parser errors
   as cuts, and the CXSMILES `c:`/`t:` reading remain. Axial stereo (allene, biaryl) is a
   separate stereo class to enable. Wedges at double-bonded ring carbons and inconsistent
   two-wedge projections have no reading in the specification and stay rejected.
3. Boundary linting: report uninterpretable or reversed marks with a proposed fix that the caller
   accepts or rejects before raise, instead of changing input. Candidates from the census: the 23
   reversed wedges, wedges at trigonal carbons, the two inconsistent projections, hyponitrous acid
   drawn as a line (CHEBI:14428), and ring double bonds drawn with a substituent on the axis
   (CHEBI 63464, 63466, 234435). Doc 036 holds the earlier SMILES diagnostics taxonomy.
4. Raising extended MOL features: R groups and wildcard atoms into the existing partial IR forms
   (`ElementForm::Undetermined`, element sets, open attachments). The inventory above bounds the
   observed input; no mapping or public conversion surface is selected. Doc 153 T3 and T5 own the
   related boundary-storage and raise items.
5. Partial-input resolution: the gates listed above stop all work on the first non-literal
   element and never publish a partially refined molecule. Publication and reporting semantics
   are undecided. The 66 unknown-stereo records wait on this.
6. Transition-metal valence: the 45 failing metal atoms all encode the oxidation state as formal
   charge plus covalence. The atom-typing registry holds one spin state per encoding (S10:
   Fe(III) high-spin, Co(II) low-spin, Mn(IV), vanadyl); further spin states and the
   Holleman-Wiberg check of the electron counts are open. The census configuration still selects
   the frozen MDL counts table.
7. The CHEBI:40 cross-format canonical discrepancy.
8. Editor property acceptance (`M  ZZC`) and overlong atom lines (CHEBI:30212, CHEBI:57503);
   presets are not changed for them.
9. Pseudoatom symbols: the specification defines periodic-table symbols, L, A, Q, *, LP, and R#.
   The extended preset parses all of them (LP under `ELECTRONS`, R# under `RGROUPS`) plus the
   vendor wildcards X and M under `WILDCARDS`, which the vendored specification does not define.
   `PSEUDOATOMS` is in no preset but `EXTENDED_MAX`, so `e` and `hv` are rejected by the basic,
   extended, and lenient presets already. No change.
10. Parser selection: no strategy decides when the basic or the extended parser is used.
    Candidates are always-extended with conversion to the basic record when no extended feature
    is present, basic with fallback to extended, or content sniffing. Doc 153 T3 owns the merge of
    the two result types. Design discussion needed.
11. Stereo bonds the stereo model keeps although no configuration is isolable: 628 chain double
    bonds whose one end carries two constitutionally equivalent ligands, and 39 Kekulé-drawn meso
    bonds of porphyrins and chlorins in a 16-membered aromatic macrocycle. `cis_trans_capable`
    distinguishes ligands by identity, and stereo resolution precedes aromaticity perception. The
    SMILES path behaves the same; whether the stereo model should exclude either class is a
    stereo-model question, not a reader question.
12. Aromatic-written metalloporphyrins: cobalt(II) porphine with aromatic ring bonds, two N+
    and two neutral ring nitrogens all single-bonded to Co, is contradictory in every resolution
    cell (aromaticity inconsistency at a ring carbon), while the Kekulé-written input of the same
    molecule resolves. Porphine itself resolves in both forms. The census heme and cobalamin
    records are Kekulé-written MOL, so their aromatization after the metal rows is untested.
