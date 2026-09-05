# 217 — Rhea participant rejection census

Status: Proposed
Date: 2026-08-30
Relates: [153](153-format-parsing-outstanding-tasks-2026-07-18.md),
[168](168-api-hygiene-2026-07-27.md),
[216](216-canonicalization-performance-2026-08-30.md)

## Purpose

Doc 216 used the ChEBI participants of Rhea release 141 to construct a varied molecule corpus for
canonicalization measurements. This document classifies every participant that did not reach the
stored resolved cohort. The purpose is diagnostic evidence about that fixed census, not maximal
ingestion coverage: abstract participants, unsupported chemistry models, and explicitly unknown
stereochemistry are legitimate reasons for exclusion from a corpus that requires concrete
molecules.

The analysis nevertheless found two implementation defects and incomplete resolver diagnostics.
Those findings leave focused future work, so the document is `Proposed` rather than
`Informational`. It does not change the corpus definition or the canonicalization work in doc 216.

## Frozen input and method

The census used the same inputs and configuration as doc 216 S0a:

| Input | Release-141 evidence |
| --- | --- |
| Rhea release | 141, dated 2026-06-10 |
| Molfile archive | `rhea-mol.tar.gz`, 5,172,329 bytes, xxh3-128 `ef5cb65c0becec63f39da049e3639566` |
| participant list | `chebiId_name.tsv`, 589,145 bytes, xxh3-128 `2e530735f3d6e0b185a3db1e71929bb7` |
| CTfile configuration | `CtfileIoConfig::basic()` |
| chemistry model | `ValenceModel::mdl()`, Daylight aromaticity, default stereo model |
| resolver configuration | `ResolveConfig::default()` |

The manifest contains one outcome for every listed participant. The analysis normalized messages
by stage, inspected the exact first rejected CTfile line for every parse rejection, and correlated
every reported stereo atom with its raw atom parity and incident wedge records. The CTfile's first
bond atom was retained separately during that analysis because it is the focus, or narrow end, of a
V2000 wedge. Representative source records were then inspected within each complete classification.

The stored manifest reports 12,950 listed participants, 14,048 Molfiles, 12,942 joined participants,
and 10,107 resolved records. The arithmetic closes exactly:

| Stage | Rejected | Complete message classification |
| --- | ---: | --- |
| source read | 8 | missing listed Molfile |
| CTfile parse | 1,464 | 1,448 extended atom symbols; 16 `M  ZZC` editor properties |
| TableIR raise | 56 | 55 wedge conflicts; 1 tetrahedral ligand-count failure |
| resolution underdetermined | 378 | explicit undetermined atom or bond stereo |
| resolution contradictory | 937 | 894 tetrahedral stereo failures; 43 valence failures |
| **Total** | **2,843** | all listed participants not stored in the resolved cohort |

There were no `resolve_execution`, `resolve_non_ground`, unexpected, or aromatization failures.
All 10,107 resolved records also aromatized successfully under the selected transformation.

## Source and parse exclusions

### Missing listed records

The eight source-read failures are participant-list entries without matching archive members:
CHEBI 4194, 24431, 59720, 84139, 131859, 137328, 140159, and 229467. They are generic classes such
as “a D-hexose” and “a hydroperoxyeicosatetraenoate,” not evidence of an I/O failure while reading an
existing file. The archive separately contains 1,106 Molfiles that are not in the participant list;
those are outside this participant census.

### Basic-parser boundary

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

These are expected exclusions under `CtfileIoConfig::basic()`, whose optimized TableIR result holds
ordinary atoms. The extended parser has representations for R groups, query wildcards, and
pseudoatoms under the corresponding flags, but those values would not be concrete molecule inputs
for this canonicalization corpus. CHEBI 685, an R-group phospholipid, and CHEBI 10545, an electron,
are representative records.

The remaining 16 parse outcomes stop on `M  ZZC`, an ACD/ChemSketch label property. The repository
already parses this editor extension when `EDITOR_EXTENSIONS` is enabled; the selected basic preset
intentionally does not. CHEBI 11302 and CHEBI 11909 are representative. Enabling that one property
could allow these records to proceed to later stages, but maximizing the admission count is not a
goal of this census.

## TableIR raise findings

Fifty-five records report `inconsistent wedge bonds at atom N`; one record reports a tetrahedral
atom with two rather than three or four ligands. Comparing the reported atom with the raw first and
second atoms of every incident wedge separates the cases:

| Reported-site relation to raw wedge | Count | Interpretation |
| --- | ---: | --- |
| only the wide endpoint | 20 | wedge applied to an endpoint it does not describe |
| both a narrow and a wide endpoint | 31 | an unrelated incoming wedge pollutes an actual focus |
| only the narrow endpoint | 4 | multiple outgoing wedges disagree under the current geometric reading |
| explicit atom parity on a two-neighbor atom | 1 | source tetrahedral marker has no valid frame |

The first two rows are manifestations of the TableIR endpoint defect described below. CHEBI 15393
has two wedges whose raw first atoms are distinct stereo sites and whose common second atom is
incorrectly treated as another focus. CHEBI 58540 has one incoming and one outgoing wedge at the
reported atom; only the outgoing wedge belongs to that atom.

The four narrow-only conflicts are CHEBI 34596, 35306, 219836, and 231826. In each, two raw outgoing
wedges imply opposite configurations under the current coordinate calculation. The census alone
does not distinguish inconsistent source depictions from a remaining geometric-perception defect;
they require an independent CTfile stereo oracle. CHEBI 80541 supplies explicit parity on an atom
with only two neighbors and is directly non-realizable as tetrahedral stereo.

## Resolution findings

### Explicitly undetermined stereo

All 378 underdetermined outcomes contain a source marker for stereo whose configuration is
explicitly unknown:

| Source markers in one record | Records |
| --- | ---: |
| atom parity 3 only | 363 |
| double-bond stereo code 3 only | 12 |
| both | 3 |

The raised `#T` or `#C` assertion therefore contains an undetermined coset and the stereo resolver
correctly cannot publish a concrete result. CHEBI 7 is representative of atom parity 3 and CHEBI
23316 is representative of crossed/unknown double-bond stereo.

Every returned `ResolveReport` has an empty `unresolved` atom-completion map; 29 also have no
tie-break entries. The report does not state that partial stereo caused the outcome. The complete
classification required returning to the source records, so the report is insufficient as a
pipeline diagnostic even though the outcome is semantically correct.

### Tetrahedral stereo contradictions

Of the 894 tetrahedral stereo contradictions:

- 886 name the raw wide endpoint of a wedge. The source describes stereo at the wedge's first atom,
  but raise also creates a `#T` assertion at its second atom. CHEBI 40 is representative.
- five name the raw narrow endpoint of an outgoing wedge on a carbon with three neighbors and
  localized-bond valence four;
- three come from explicit atom parity on the same three-neighbor, valence-four carbon shape.

The latter eight sites are trigonal, not tetrahedral, under the localized source structure. CHEBI
32446 is representative of the outgoing-wedge group; CHEBI 34463 is representative of the atom
parity group. They are source-marker/model incompatibilities, not aromaticity failures.

### Valence contradictions

All 43 `no matching valence state` outcomes contain a transition metal that is outside the frozen
MDL counts table:

| Transition metal present | Records |
| --- | ---: |
| Fe | 20 |
| Co | 18 |
| Mn | 2 |
| V | 1 |
| Cu | 1 |
| Cr | 1 |

The set includes isolated ions such as CHEBI 29034, iron and manganese oxides, and coordination
complexes such as CHEBI 4991 and CHEBI 16304. This is selected-model scope, not malformed TableIR
and not an aromaticity or stereo failure. The generic contradiction message does not name the atom,
element, charge, or valence that failed even though that information is useful and the counts
module already has `CountsMismatch` vocabulary.

No resolver outcome reports an aromaticity contradiction. Aromaticity therefore contributed no
separate rejection category in this release/configuration pair.

## Confirmed wedge-endpoint defect

V2000 wedge direction is endpoint-relative: the first atom in the bond record is the narrow end and
the stereo focus. TableIR documents this rule on `BondWedge`, but `Bond::new` and
`ExtendedBond::new` immediately normalize the endpoints through `AtomPair::new`. The retained
`first()` is then the smaller atom id, not necessarily the source first atom. More decisively,
`wedge_bond_neighbors` treats every incident `Up` or `Down` wedge as applying to the atom currently
being examined; it never checks which endpoint bears the wedge.

This explains 937 rejections across two reported stages: 51 raise conflicts and 886 later stereo
contradictions. The same representation gap was previously noted for CXSMILES `w:` entries, which
supply an explicit atom endpoint but currently retain only the bond-wide `BondWedge` value.

The future data-type contract is settled only to the following extent:

- TableIR bonds are open external-boundary records, not graph-IR localized bonds.
- Their topological endpoints may remain an unordered normalized pair, but a present wedge must
  retain its focus endpoint as part of the source representation.
- CTfile parsing, CXSMILES application, basic/extended conversion, atom-id updates, and TableIR raise
  must preserve that endpoint without inferring it from atom-id order.
- Raise reads a wedge only at its retained focus. The other endpoint does not acquire a `#T`
  assertion merely by incidence.
- Independently supplied malformed endpoint references fail at the first conversion or consumer
  that requires them; they do not become indexing panics.

The public representation, field name, constructor surface, and compatibility strategy remain
unsettled. No new public seam should be selected merely to make the fix convenient. The complete
`Bond` and `ExtendedBond` construction and conversion surfaces must be reconciled before an
implementation plan is written.

## Panic audit

The release-141 run did not panic. The builder does not wrap individual records in
`catch_unwind`; it processed all 12,950 listed participants and wrote the complete manifest. The
archive also contains no bond stereo code outside `0`, `1`, `3`, `4`, and `6` in any of its 14,048
raw V2000 records.

There is nevertheless a reachable malformed-input panic in the same parser path. Both basic and
extended bond parsers accept an integer stereo field and pass it directly to
`convert_bond_stereo_direction_code`, whose default arm is `unreachable!`. A standalone probe with
bond stereo code 2 reached:

```text
internal error: entered unreachable code: invalid stereo/direction code: 2
```

That source-controlled value must produce `ParseError`, not panic. The raise and resolver path also
contains `expect` and `unreachable!` sites, but this audit found them downstream of explicit
frame, non-empty-candidate, store-membership, or prior-policy invariants. No other malformed-record
route to a panic was demonstrated. This is a focused audit of the exercised ingestion path, not a
claim that every parser configuration or resolver model is panic-free.

## Root-cause disposition

The complete rejection population partitions as follows:

| Disposition | Records | Included causes |
| --- | ---: | --- |
| expected boundary, source, or selected-model exclusion | 1,893 | 8 missing records, 1,464 basic-parse exclusions, 378 unknown stereo outcomes, 43 transition-metal valence exclusions |
| confirmed umol wedge-endpoint defect | 937 | 51 raise conflicts and 886 stereo contradictions |
| source/perception compatibility requiring an oracle | 13 | 4 conflicting outgoing-wedge depictions, 1 parity site with two neighbors, 8 tetrahedral markers on trigonal carbon |
| **Total** | **2,843** | complete release-141 participant rejection census |

The 13-record oracle group is kept separate because the structural reason for rejection is known,
but assigning fault to source data or umol perception is not justified without an independent
reading of the same CTfile stereo descriptors.

## Proposed follow-up scope

The findings support four focused work units, but their public design and sequence are not yet an
implementation plan:

1. Settle and implement endpoint-bearing wedge representation across basic and extended TableIR,
   CTfile and CXSMILES ingestion, id updates, conversions, and raise. Verification must include
   swapped numeric endpoint order, incoming wedges that must not create stereo, adjacent true
   stereo sites, and redundant consistent and inconsistent outgoing wedges.
2. Replace the invalid bond-stereo `unreachable!` path with an exact `ParseError` and table-driven
   parser cases for every unsupported code in range. Malformed external bytes must not panic.
3. Make resolver diagnostics identify the underdetermined phase and cause, and retain atom-local
   context for counts-valence mismatch. The design must preserve the distinction between semantic
   `Solution` outcomes and operational `Result` errors.
4. Compare the 13 residual stereo-marker records with an independent CTfile stereo implementation,
   then rerun this exact release/configuration census after the selected fixes. Counterfactual
   admission counts are evidence, not a requirement to accept abstract participants or chemistry
   outside the chosen models.

Supporting query/R-group ingestion, extending the MDL valence model to transition metals, or
turning this census into a completeness campaign is explicitly outside this document.
