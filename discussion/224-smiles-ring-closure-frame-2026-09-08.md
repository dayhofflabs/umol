# 224 — SMILES ring-closure stereo frame

Status: In Progress
Date: 2026-09-08
Relates: [104](104-stereochemistry-implementation-plan-2026-05-31.md),
[153](153-format-parsing-outstanding-tasks-2026-07-18.md),
[217](217-rhea-participant-failures-2026-08-30.md)

## Purpose

The SMILES tetrahedral reading can be inverted at stereocenters that carry ring-closing digits.
Opening-order bond storage is retained for its source-order meaning and usefulness when
reconstructing SMILES. The proposed finalizer publishes complete stereo frames that raising can map
directly to graph IR. The work separates behavior-preserving builder streamlining, the stereo-frame correction, and limited
permutation tests. General permutation/randomization testing is out of scope.
The 2026-09-09 impact review below preserves the evidence behind this decision.

## Finding

OpenSMILES §3.8.2 defines the frame of a tetrahedral center as the order of its bonds in the
string: the bond to the preceding atom first; an implicit hydrogen next, or first when the center
opens the SMILES; then, left to right after the bracket, each ring-closure digit at its position,
each branch, and the continuing chain. `@` lists the neighbors after the first anticlockwise when
viewed from the first toward the center, `@@` clockwise. "One exception to the atom order is when
these atoms are bonded to the chiral center via a ring bond. In these cases, it is the order of
the bonds to these atoms that should be considered."

umol takes the order from the TableIR bond list (`first_neighbor_toward_ordering` through
`AtomNeighbors`). The SMILES builder reserves a ring-closure bond's slot at its opening digit and
fills it at the closing (doc 104, 2026-06-07). At the opening atom the list order is the written
order. At the closing atom the ring bond precedes the bond to the preceding atom, one
transposition away from the rule, so the coset is inverted there.

Minimal pair, same isomer by the rule: `C[C@H]1CCCCO1` has the frame (CH3, H, O, CH2) with `@`;
`O1CCCC[C@@H]1C` has (CH2, H, O, CH3) with `@@`, and swapping CH3 and CH2 turns `@@` into `@`.
umol reads the second in the frame (O, H, CH2, CH3), yielding the mirror image. Asked directly,
umol reports the pair as canonically unequal and reports `C[C@H]1CCCCO1` equal to its mirror
image `O1CCCC[C@H]1C`. Controls whose digits open rings at the center, such as
`C[C@]1([H])CCCCO1` against `[H][C@@]1(C)CCCCO1`, read correctly.

The census surfaced it through CHEBI:40, (+)-pinoresinol, whose Rhea SMILES closes rings at
two of its four centers; the MOL wedge reading of the same record agrees with an independent 3D
reading at all four centers (doc 217, "Paired MOL and SMILES input for CHEBI:40").
Three resolution fixtures were lowered from affected SMILES and encode the wrong isomers:
`stereo_tetrahedral/cis-1-2-dichlorocyclohexane.edn` (from `Cl[C@H]1CCCC[C@H]1Cl`) holds the
chiral trans isomer, `trans-1-2-dichlorocyclohexane.edn` the meso cis isomer, and
`alpha-d-glucopyranose.edn` the C4 epimer, alpha-D-galactopyranose. A stereocenter at a
ring-closing digit can be read incorrectly, including in the usual writing of pyranoses.

## No bond order carries the frame

For `X A1(B) ... Y Z1 D` with the ring bond between the opening atom A and the closing atom Z,
the rule needs the ring bond before A–B at A and after Y–Z at Z, while A–B is written before
Y–Z. The three inequalities form a cycle, so no single position in one bond list satisfies both
ends whenever the opening atom has any bond written after its digit. Open order serves the
opening atom, close order the closing atom. The frame therefore has to be recorded by the parser,
the only place that sees it. This does not erase the boundary meaning of bond encounter order: the
builder retains opening order, while private sparse records recover the closing-end order for the
complete frame written at finalization.

## Agreed design — 2026-09-09

The design and plan are approved for the S0 experiments. Production parser changes remain pending
review of the revised prototype and subsequent approval. The earlier public ring-position and bond-table-boundary proposals are superseded by
complete explicit stereo frames; their experimental results are retained below as historical evidence.

### Parser state and finalization

Retain opening-order bond storage, the existing pending-ring table, and the optional CX completion-rank
map. Extra ring-encounter bookkeeping is limited to closing digits at supported tetrahedral-marked
atoms. Marked traversal roots are recorded separately as specified below. The sparse prototype records the atom, reserved ring-bond slot, and current bond-table
length; equal boundaries retain digit order. These records are private parser state, not TableIR data.
There is no per-atom counter vector or accumulator updated for every ordinary bond.

Finalization reconstructs the complete ordered ligand frame for each supported tetrahedral descriptor,
including centers without ring closures. It inserts and identifies virtual ligands according to the
input convention defined below, associates the configuration with that frame, and discards the
private ordering records. An opening-only or acyclic stereo center needs an output frame but no extra ring bookkeeping.
Ordinary atoms and unmarked ring atoms need neither. This distinction matters for performance:
private ring records scale with affected closing digits; published frames scale with stereo sites.

The output frame uses the same participant-frame meaning as graph IR: a site, ordered actual/virtual
ligands, and a configuration interpreted in that order. It does not carry instructions for interpreting
SMILES, a local-position patch, or a bond-vector boundary. It need not select graph IR's representative
frame. Explicit frames remove the need for a source-convention flag at consumption.

Finalization is part of the parse benchmark. Reconstruct only the needed frames, skip that work when
there are no supported descriptors, and measure the collection strategy before choosing it. Neither
repeated whole-molecule scans nor a newly allocated all-atom neighbor table is assumed cheap. No generic
accumulator framework or chemistry-model dependency is introduced into the low-level parser.

### TableIR and raising

TableIR carries complete operative frames separately from raw source descriptors. Basic/extended
conversion preserves the frames, participant order, ligand kinds, and associated configurations.
Empty construction supplies no operative frames. The public surface and consumption rules are
specified below; implementation stages sequence this design rather than decide it.

Raising maps an operative frame directly into `MoleculeEntries::stereo_atoms`, preserving its ligand
order and associated configuration. `Molecule::try_from_entries` remains the authoritative graph-IR
integrity gate. There is no SMILES ring-order reconstruction, virtual-ligand insertion, or transport
into sorted-neighbor order in the new explicit-frame path.

Remove ChiralityFrame and both parsers' assignments to it. MOL atom parity remains raw metadata and
never becomes operative merely because it is present. The existing MOL wedge/coordinate path and its
unknown/conflicting-stereo handling remain unchanged. ConfigurationScope, directional-bond stereo,
CX indexing, and other source metadata are not removed with ChiralityFrame.

The current `table_ir::raise::raise_tetrahedral_stereo` publishes a `#T` atom constraint rather than
an explicit stereo entity. The new path changes that representation as well as fixing the source
frame. The old 516-case comparisons therefore do not predict all new raw output differences or prove
equivalent resolution behavior. Check explicit-entity handling in `ops::resolve::stereo` and compare
resolved semantics; do not silently retain a duplicate `#T` assertion to imitate old outputs.

### Public frame contract

Follow the existing TableIR convention: boundary-owned records, plain indices, and literal values.
Define these public types in `table_ir::stereo`:

```rust
pub struct StereoAtom {
    pub atom: u32,
    pub ligands: Vec<StereoLigand>,
    pub winding: Winding,
}

pub enum StereoLigand {
    Atom(u32),
    ImplicitHydrogen,
    LonePair,
}

pub enum Winding {
    Clockwise,
    CounterClockwise,
}
```

Add `pub stereo_atoms: Vec<StereoAtom>` to both Molecule and ExtendedMolecule. The atom and actual
ligand indices refer to the owning atom vector. Virtual ligands belong to StereoAtom.atom; their
anchor is not stored redundantly. The first ligand points toward the viewer; the remaining three
ligands wind clockwise or counterclockwise as specified. This is a tetrahedral record, not a generic
configuration-expression carrier. The ordered frame has the same meaning as the graph-IR frame,
but no graph-IR ids, forms, constraints, or entry tuples are stored in TableIR.

The record has public fields and no additional constructor or validator. StereoAtom derives Clone,
Debug, PartialEq, and Eq; the enums also derive Copy. Export the three types through the existing
`table_ir::stereo` module and TableIR's existing stereo re-export. This is an open boundary record,
not a promise that independently assembled frames are integrity-valid. A Vec retains incomplete
or oversized input frames for checked raising; it does not imply arbitrary tetrahedral arity.

`empty()` initializes an empty stereo_atoms vector. The existing Molecule → ExtendedMolecule
conversion moves the records, and the checked reverse conversion clones them without reordering,
validating, or repairing them. Existing source-wrapper construction stays unchanged. Complete
struct literals must supply the new field. Other source stereo kinds are not added by this work.

Keep Atom.chirality and ExtendedAtom.chirality as raw source metadata. StereoAtom.winding in its ordered
ligand frame is the authoritative operative value; raising never reinterprets, compares against,
or derives it from the raw field. Editing raw metadata alone does not change an explicit frame.
SMILES finalization preserves the token and emits the associated row. MOL parsing preserves parity
but emits no row. Remove ChiralityFrame, both molecule fields, all assignments/conversions, and its
existing TableIR export. ConfigurationScope is unchanged. No Python signature or new Python
constructor is added.

### Finalization rules

Every `@`, `@@`, `@TH1`, or `@TH2` atom receives a row, including syntactically accepted inputs with
an incomplete frame. Emit records in atom-index order. Set winding to CounterClockwise for
`@`/`@TH1` and Clockwise for `@@`/`@TH2`. Winding is relative to the emitted source frame; do not
sort that frame or transport the configuration into another frame during raising.

Actual ligands are StereoLigand::Atom values containing the neighboring atom's u32 index, in written
bond order recovered from opening-order storage plus the private closing records. Preserve source
incidences rather than silently deduplicating invalid topology. Explicit `[H]` neighbors are ordinary
actual ligands. Virtual ligands are anchored at the stereo site and follow these rules:

- A bracket hydrogen count contributes that many ImplicitHydrogen entries. The parser already stores
  an omitted bracket count as zero. Repeated entries share the record's implicit site anchor; do not invent
  distinct virtual occurrences to make an invalid descriptor pass integrity checking.
- With exactly three actual incident ligands and bracket hydrogen count zero, emit one LonePair
  entry. This is the declared fourth participant of the supported three-coordinate tetrahedral
  convention, not an inference that the atom has a chemically valid lone pair. Do not write
  Atom.lone_pairs, infer electron counts, or call a valence model. Resolution checks the supplied
  entity against the chemically derived frame and can reject it.
- Otherwise add no unspecified padding. A frame of the wrong arity remains in the open TableIR
  record and is rejected at checked raising. Supplying a descriptor never causes a missing atom or
  hydrogen to be invented.

Insert virtual entries before the first actual ligand when the marked atom starts a source traversal,
and immediately after its incoming chain ligand otherwise. The hydrogen placement follows the
[OpenSMILES tetrahedral convention](https://github.com/opensmiles/OpenSMILES/blob/master/opensmiles.asciidoc).
Using the same slot for the lone pair retains umol's existing virtual-slot convention; it is not a
claim that OpenSMILES specifies lone-pair placement completely. Record marked traversal-root ids
privately when the parser has no preceding atom, including after a dot. Do not infer this from global
atom index zero or from graph connectivity after ring closure. This adds information only for marked
roots, not a per-atom vector. Discard root and ring-order records after finalization.

The old constraint path deliberately left virtual identity to resolution (doc 104 B8). Direct typed
frames make that identity explicit without asserting a hydrogen/lone-pair count on the atom. This
is a boundary-representation decision, and possible changes in invalid-input diagnostics or
constraint-versus-entity failure policy are part of the prototype comparison, not hidden compatibility
requirements. An unresolved chemical model result remains downstream; it is not a parser decision.

For example, `N[C@H](F)Cl` emits `[Atom(0), ImplicitHydrogen, Atom(2), Atom(3)]` at atom 1
with CounterClockwise winding; `[C@H](F)(Cl)Br` emits
`[ImplicitHydrogen, Atom(1), Atom(2), Atom(3)]` at atom 0 with the same winding.
`C[S@](=O)CC` uses `[Atom(0), LonePair, Atom(2), Atom(3)]` at atom 1. Virtual-ligand anchors are
supplied from the record's atom field during raising. Frame construction does not itself certify
any of these centers as chemically stereogenic.

Non-tetrahedral raw descriptors, including AL/SP/TB/OH, remain parseable where already supported
but produce no operative row in this task, preserving the current scope. Existing lexical range
errors such as `@TH3` remain parser errors. Do not conflate unsupported raw descriptors with supported
but incomplete tetrahedral descriptors: `[C@]` produces a tetrahedral row with an empty ligand list,
so its descriptor is retained and checked raising fails rather than succeeding with stereo erased.

### Raising, absence, and errors

Raising converts each TableIR record into a `MoleculeEntries::stereo_atoms` entry without changing
its frame. Convert atom indices to graph-IR AtomId; convert actual ligands to graph-IR StereoLigand
values of kind Atom and virtual ligands to their respective kinds anchored at the record's atom.
Map CounterClockwise to tetrahedral coset 0 and Clockwise to coset 1 using StereoAtomForm::new.
These are representation conversions, not reconstruction or frame reordering.

Call `Molecule::try_from_entries` to publish the result. Raising must not index a supplied
site/ligand before that publication gate, sort or deduplicate the frame, repair incomplete records,
or synthesize duplicate `#T` constraints. No new RaiseError variant is required: existing
RaiseError::MoleculeEntries preserves integrity errors for references, duplicate sites/ligands,
incidence, and arity. An invalid kind or out-of-range coset cannot be supplied through these TableIR
types. Examples with a single defect include:

| Input | Required checked-raise result |
| --- | --- |
| `[C@]` | MoleculeEntries wrapping StereoLigandArity: tetrahedral, expected 4, actual 0. Parsing itself remains successful. |
| `[C@H2](F)Cl` | MoleculeEntries wrapping DuplicateStereoLigand for the repeated implicit hydrogen at the site; no distinct virtual identities are manufactured. |
| Independently supplied row naming an absent atom | MoleculeEntries wrapping InvalidReference. |
| Two rows with the same valid site | MoleculeEntries wrapping StereoAtomSitesDuplicate. |

For any atom named as an explicit stereo site, that row is authoritative and the atom-level wedge
interpretation is skipped. Duplicate or malformed explicit rows still fail at the gate; do not
fall back to wedges to repair them. Without a row, retain the existing wedge/coordinate path,
including unknown and conflicting stereo handling. Raw metadata neither activates nor suppresses
that path. This also defines independently assembled TableIR containing both explicit rows and
wedges: the explicit atom row takes precedence. Directional-bond stereo is a separate mechanism
and continues to be interpreted normally.

An explicit frame declares participants/configuration, not chemical validity. Existing model scope,
valence, stereo perception, and stereo-entity failure policies apply during resolution. Do not change
those policies to reproduce old constraint-path diagnostics. The prototype must expose differences;
a broader resolver change would require its own decision.

### Preservation contract and public-symbol inventory

Basic/extended conversion preserves records exactly. Atom remapping changes the atom field and
actual ligand indices; virtual ligands follow their owning site without a separate anchor field.
Bond-vector permutation alone does not change an atom frame. Reordering the tetrahedral ligand
frame must also transport winding: odd permutations flip it and even permutations preserve it.
Use the existing graph-IR frame-transport API when testing the raised form. No new production
permutation operation or independent frame algebra is introduced for tests.

| Public surface | Change |
| --- | --- |
| StereoAtom, StereoLigand, Winding | Add the TableIR-owned record and enums above in table_ir::stereo and its existing re-export. |
| Molecule.stereo_atoms, ExtendedMolecule.stereo_atoms | Add Vec<StereoAtom>. |
| Molecule::empty, ExtendedMolecule::empty | Initialize it empty; signatures unchanged. |
| Existing basic/extended conversions | Preserve complete rows; signatures and existing conversion errors unchanged. |
| ChiralityFrame and both chirality_frame fields | Remove, including the wildcard TableIR export and all callers. |
| Atom.chirality, ExtendedAtom.chirality | Retain raw metadata; no operative interpretation in raising. |
| TryIntoIr for the TableIR molecule | Convert the boundary records to graph-IR entries and retain checked publication; signature unchanged. |
| RaiseError | Reuse MoleculeEntries; keep errors still used by the unchanged wedge/directional paths. |
| SMILES/MOL parse APIs and Python bindings | Signatures unchanged; SMILES adds rows, MOL leaves them empty. |
| Graph-IR types and constructors | Reuse as-is; no new graph-IR symbols or invariant changes. |

These are design decisions made before implementation. S0 measures and exercises them; it does not
select their public shape or error contract.

A complete explicit frame also addresses the targeted ChiralityFrame redesign in doc 153 T4. Broader
format metadata work and T9 arrangement-number verification remain outside this task. In particular,
this plan does not expand support to other stereo kinds.

## Scope and verification properties

The three work units remain distinct: stereo correction, limited permutation coverage, and
behavior-preserving opening-order builder cleanup. No closing-order migration, SMILES writer,
general randomization suite, or broader resolver redesign is included.

### Focused stereo and permutation coverage

Coverage includes both tetrahedral atom stereo and directional-bond cis/trans stereo. The latter
continues through the existing directional-bond interpretation and `#C` constraint path; adding its
permutation tests does not migrate it to the new explicit atom-frame mechanism.

Retain exact parser expectations and independent expected frames. Cover the known equivalent/mirror
pair, explicit and virtual hydrogen, source-defined lone-pair cases, opening-only centers, closing
centers, mixed digits, several closures at one center, branch/digit order, `%nn`, digit reuse,
dot components and reaction sides, and direction/donation/CX controls. Check virtual-ligand placement
at component starts explicitly under the rule specified above; do not transplant the old
`atom_idx > 0` heuristic into finalization. Classify changes at later traversal roots separately
from closing-ring corrections in the impact comparison.

The limited transformation properties are:

- **Equivalent atom-stereo source presentations:** use a small named set of ring-label, branch-order, and
  opening/closing presentations with independently specified atom correspondence and stereo meaning.
  Adjust descriptor parity for the known ligand permutation. Compare raised explicit structures
  using remapping-aware frame equivalence, then check canonical agreement after ingestion. The
  mirror control must remain distinct; fingerprints alone are not a stereo oracle.
- **Atom-stereo bond storage order:** on small explicit-frame TableIR fixtures, apply selected ring-bond moves,
  swaps, and a reversal while transporting any bond references. Leave atom ids and explicit ligand
  frames unchanged. Compare raised values under the induced entity mapping and canonical resolved
  results. Cover both opening and closing centers. The earlier requirement to preserve non-ring
  written order is no longer needed for the atom stereo frame; restrict fixtures or transport other
  order-sensitive metadata so the test changes only storage order.
- **Equivalent bond-stereo source presentations:** use a bounded set of independently specified
  E/Z pairs. Vary traversal direction, branch order and the choice of marked substituent, adjusting
  `/` and `\` for the resulting directed bonds. Include simultaneous reversal of both directional
  signs, marked ring bonds expressed at opening versus closing, and a small conjugated example
  whose directional bond participates in two double-bond readings. Retain opposite-E/Z controls
  (for an isolated alkene, reversing just one side's sign), plus unmarked and one-sided-marker
  cases that must not acquire a determined configuration. Check basic and extended parsing.
- **Bond-stereo storage and endpoint order:** on the same small TableIR cases, move the double bond
  and its directional neighbors within the bond vector, transporting bond references. Separately
  reverse selected stored bond endpoints and invert the direction of each reversed directional
  bond. Preserve the represented geometry; compare the raised cis/trans constraints in their
  corresponding frames and the resolved stereo bonds under the induced entity mapping. Do not
  compare raw coset numbers across changed endpoint/ligand frames. Include exact dangling-direction,
  contradictory-direction, and mismatched ring-direction controls with expected remapped ids or
  source positions where applicable. These controls exercise the existing parser/raise failures.
- **Conversion preservation:** basic/extended conversion preserves complete frames and descriptor
  association exactly. Where a tested atom relabeling is needed to express a source correspondence,
  remap each stereo record's site as well as actual ligands, so its virtual anchors follow. For
  stereo-bond cases, preserve directional
  annotations and their endpoint meaning through the basic/extended conversion, then check the same
  cis/trans interpretation after raising and ingestion.

These are bounded fixture transformations, not new arbitrary graph/SMILES generators. Keep the
operational domain and comparison relation explicit; use existing `framed_eq_under` and canonical
comparison APIs where applicable. General permutation suites for matching, reactions, persistence,
depiction, or other downstream tasks remain out of scope.

### Output and performance evidence

Compare original production behavior, the earlier corrected-opening reference, and the revised
explicit-frame prototype separately. Record parser acceptance, bond tables, published frames,
raw raised/serialized output, resolution acceptance and diagnostics, and canonical resolved results.
Separate intended stereo corrections from constraint-to-entity representation changes. Preserve the
original canonicalization bound and explicit algorithm/context when reusing the 516-case sample;
report skips separately. Do not use prior digest equality as a required law for the new raw form.

Benchmark basic and extended parsing through finalization and result destruction, raising, and
representative complete ingestion separately. Include ring-free non-stereo input, ordinary rings,
acyclic stereo, opening-only stereo, affected closures, mixed digits, a long prefix, repeated stereo
sites, and CX controls. Record allocations/state where practical, without assigning a slowdown to an
allocator or counter from timing alone. The old short measurements are not a performance guarantee,
representative molecular distribution, or an agreed regression threshold.

## Earlier sparse-boundary prototype — 2026-09-09

This experiment exposed private bond-table boundaries on TableIR and reconstructed frames during
raising. That public shape is superseded by the agreed full-frame finalization design. Its results
measure the earlier design, not the revised finalization or direct stereo-entity publication.

The experiment used an isolated source copy and a separate Cargo target. No production source,
fixture, or test was changed. The prototype did not implement complete public-input validation or
migrate the existing test literals; it was not a production implementation.

- All 516 comparison rows exactly match the earlier corrected-opening prototype's recorded
  digests: bond tables, raising, resolution, bounded canonicalization, fingerprints, and error
  details. This reuses the same domain and canonicalization limit described in the impact review;
  it is evidence of equivalent correction, not an independent oracle for all stereo semantics.
- Twelve exact-record cases pass through both basic and extended parsing: ring-free and unmarked
  controls, opening-only stereo, one closing digit, mixed digits in both orders, two closures at
  one boundary, a branch before closure, closure after returning to an earlier marked atom,
  `%99`, digit reuse across components, and an explicit tetrahedral marker.
- A separate release timing harness compared the original source with the sparse boundary design.
  Each sample repeats parsing and result destruction; each case uses seven 50 ms samples. The
  binaries ran original/sparse/sparse/original, without concurrent builds. The table gives median
  nanoseconds from the second pass; percentage ranges cover both passes. These are short local
  microbenchmarks, not Criterion confidence estimates or claims about representative molecules.

| Input | Basic original → sparse (ns) | Basic change, both passes | Extended original → sparse (ns) | Extended change, both passes |
| --- | ---: | ---: | ---: | ---: |
| Chain, 10 C | 413 → 418 | +0.9–1.2% | 615 → 605 | −1.6 to −0.8% |
| Six-membered ring | 382 → 391 | +1.6–2.3% | 479 → 492 | +1.2–2.7% |
| Fused rings | 464 → 474 | −2.4 to +2.1% | 680 → 694 | +1.4–2.0% |
| Adamantane | 514 → 504 | −2.0 to −0.1% | 723 → 741 | +0.7–2.5% |
| Stereo opening only | 374 → 371 | −0.9 to −0.3% | 522 → 525 | +0.6–2.0% |
| Stereo closing | 374 → 404 | +5.9–8.0% | 519 → 557 | +5.6–7.2% |
| Mixed opening/closing | 458 → 490 | +5.0–6.9% | 671 → 704 | +4.5–4.9% |
| 1,000 C prefix, then stereo closing | 14,515 → 14,675 | +1.1–1.3% | 32,746 → 33,085 | +1.0–1.3% |
| 10 stereo-closing components | 1,621 → 1,812 | +11.8–12.1% | 2,827 → 2,910 | +2.7–2.9% |
| 100 stereo-closing components | 12,346 → 13,251 | +5.4–7.3% | 25,890 → 25,958 | −0.5 to +0.3% |

The earlier 13–24% broad slowdown is absent from the unrecorded cases in this experiment; their
changes stay within approximately 3%. Affected inputs still pay for sparse recording and allocation:
about 6–8% in the one-closure case, and up to 12% in the deliberately repeated basic-parser case.
These timings do not isolate allocator cost or establish zero overhead. The long-prefix case avoids
the scan-based design's extra dependence on the already-parsed prefix. No assumptions about real
molecular size or stereo density are needed for the constant-time recording design.

The earlier 13–24% measurements below used a different harness; absolute times from the two
experiments should not be compared as paired measurements.

## Impact review — 2026-09-09

This review compared the former closing-order proposal with a corrected opening-order alternative.
The decision above retains opening order; the measurements remain evidence, not pending migration
work. Closing order is a small parser change with broad representation-output consequences. The
stereo-frame correction works with either order. The experiment found substantial bond-id and
serialized-output movement, concentrated expectation migrations, and no additional canonical
changes from choosing closing order within the measured domain. The observed failures are specific expectation changes, not demonstrated non-canonical behavior.
The limited measured timing difference does not determine the storage-order decision.

### Method and limits

Reviewed commit `be941902c7ce4ca3e1a2fbe7994bad25dd2f5334`, with a clean working tree. The
prototype used an isolated source archive. No production implementation or fixture was changed in
the working tree. Three variants distinguished the effects:

- **Original:** current parser and current frame reader.
- **Opening:** ring-end positions and the corrected frame reader, retaining opening-order bonds.
- **Closing:** the same positions and frame reader, appending bonds at closure and making CX index
  translation identity.

The prototype is an experiment, not completion of the former S0–S2 stages. In particular it retains the optional bond
buffer and legacy bookkeeping to isolate insertion timing, and does not settle validation of
independently assembled ring records. Existing output expectations were retained. Test-only
adaptations supplied the new empty field in complete struct literals and excluded the new
`ring_closures` field from old parser equality checks. A copied geometry harness needed its data
path canonicalized to resolve filesystem aliases.

Rust runs used `cargo test --workspace --features
umol-io/conformance,umol-io/proptest,umol-graph/conformance --no-fail-fast`, with the repository
Python 3.13 environment active. Python comparisons loaded freshly built wheels from separate
extracted directories using `PYTHONPATH`; the installed environment was not replaced.
Feature-gated depiction/CoordGen suites, graph property suites beyond the selected features,
fuzzing, and external application/corpus persistence were not measured.

### Existing tests and snapshots

| Comparison | Result |
| --- | --- |
| Original Rust baseline | 29,270 passed, 0 failed, 9 ignored |
| Corrected opening order | No existing-behavior failures after the harness adaptations |
| Corrected closing order | 29,052 passed, 218 failed, 9 ignored |
| Opening-order Python | 1,480 passed, 2 skipped |
| Closing-order Python | 1,471 passed, 9 failed, 2 skipped |
| Existing snapshots | No changed snapshots in either corrected variant |

The opening run initially had 181 copied-path failures and one new-metadata comparison failure;
the two targeted rechecks passed all 181 geometry cases and all 16 wildcard cases. Those are
harness adjustments, not parser regressions. No closing-order expectation was accepted or
rewritten to make the comparison green.

The 218 closing-order Rust failures are concentrated in four source files:

| Source | Cases | What moves |
| --- | ---: | --- |
| `umol-io/src/smiles/parser/tests/basic.rs` | 100 | Exact bond lists in ring, bond, component, bracket, wildcard, and lenient cases |
| `umol-io/src/smiles/parser/tests/extended.rs` | 99 | The corresponding extended-parser expectations |
| `umol-io/src/table_ir/raise.rs` | 5 | Four cis/trans tests select the former double-bond index; one diagnostic reports a moved bond id |
| `umol-graph/src/ingest.rs` | 14 | Exact resolved molecule/reaction output, including aromatic and expanded-element-scope cases |

The four cis/trans failures return `Ok(None)` because the supplied index now names another bond;
they do not demonstrate changed double-bond stereochemistry. The diagnostic changes from
`DanglingBondDirection { bond: 6 }` to `{ bond: 8 }`.

The nine Python failures are in `test_molecule.py` and `test_reaction.py`: two selenium-ring bond
list checks, three molecule aromaticity-policy cases, one reaction resolve-config case, and three
reaction aromaticity-policy cases. No binding signature changes are required, but existing Python
results and expectations do change. The earlier “no binding changes” consequence is insufficient
as an impact summary.

The unchanged snapshots are weak evidence for numbering stability:

- SMILES conformance snapshots contain success/error information and sum formula, atom count, and
  bond count. They do not contain the bond table or raised stereochemistry.
- Resolution snapshots read committed EDN through `TestInput::from_edn_str`. They do not reparse the
  SMILES that originally generated a fixture, so an unregenerated wrong isomer remains invisible.

Regenerating the complete 35-example `verify_stereo` catalog into separate temporary directories
changed **11 input EDN files** between corrected opening and corrected closing order: the three
already identified as incorrect, cis/trans-decalin, L-ascorbic acid, R-1-phenylethanol, R-methyloxirane,
cyclohexene, and E/Z-cyclooctene. This measures potential regeneration churn; it does not require
regenerating all eleven to fix the three known incorrect fixtures.

### Downstream output comparison

The deterministic sample took every eighteenth entry of the 9,151 lexicographically sorted
nonempty OpenSMILES fixture inputs and added seven explicit controls: 516 inputs total.
It is a repository sample, not an estimate of production prevalence. Outputs were compared using
deterministic digests; IR output was serialized through `MoleculeDsl` with concrete defaults,
rather than hashing internal graph/debug storage. Repeated opening runs produced identical rows.

| Output, corrected opening versus corrected closing | Changed / compared |
| --- | ---: |
| Atom counts | 0 / 516 |
| Bond-table order | 437 / 516 |
| Bond contents after sorting by endpoints | 0 / 516 |
| Raise success/failure | 0 / 516 |
| Raised EDN among successful raises | 436 / 513 |
| Resolve success/failure | 0 / 516 |
| Resolved EDN among successful resolutions | 404 / 479 |
| Canonical resolved EDN, at most 40 atoms | 0 / 459 |
| ECFP radius 2 and Morgan radius 2 feature sets | 0 / 479 |
| Resolution error text/details | 2 / 37 |

Twenty successfully resolved inputs exceeded the canonicalization size bound. The explicit
canonicalization context used Nauty with `para_stereo: false`. The two changed resolution errors
still reject the same inputs for aromatic-assertion discharge; their reported bond ids move from
1 to 0 and from 7 to 6.

Original versus corrected opening order changes **37 raised outputs and 32 canonical resolved
outputs**, with no bond-table changes. This separates the intended frame correction from optional
renumbering. The fingerprint results also stay unchanged across that correction: these particular
fingerprints do not establish stereo preservation. The minimal equivalent/mirror pair changes
exactly as described in the finding.

The propagation follows the current code directly. TableIR raising iterates the bond vector and
appends localized, dative, and noncovalent entries to their respective entity tables. Their relative
orders therefore affect ids. `ingest::interpret_molecule` resolves that result without selecting a
canonical frame. Consumers of raw bond ids, exact equality, EDN, reaction deltas/correspondences,
error references, and traversal output must consequently be audited. Canonical output is the
appropriate semantic comparison, but its agreement here does not establish stability of unmeasured
layout, matching enumeration, persisted evidence, or external callers.

### Parsing cost

Standalone release benchmark binaries were retained for all three variants. After compilation and
tests finished, the same nine existing cases were run sequentially in original/opening/closing and
reverse order, with 20 samples, 0.1 s warmup, and 0.5 s measurement per case. Rounded mean
nanoseconds from the final pass:

| Existing benchmark | Original | Corrected opening | Corrected closing |
| --- | ---: | ---: | ---: |
| Basic chain, 10 carbons | 440 | 522 | 524 |
| Basic six-membered ring | 399 | 459 | 465 |
| Basic fused-ring example 1 | 464 | 575 | 582 |
| Basic adamantane | 506 | 616 | 620 |
| Extended chain, 10 carbons | 617 | 725 | 737 |
| Extended six-membered ring | 475 | 556 | 562 |
| Extended fused-ring example 1 | 681 | 808 | 802 |
| Extended adamantane | 741 | 855 | 842 |
| Extended CX multicenter ferrocene | 788 | 889 | 888 |

Opening versus closing timing differs by less than 2% in this experiment; no material speed
advantage is established. Both straightforward prototypes cost approximately 13–24% more than the
original on these cases. The shared per-atom counting and ring-record allocation deserve attention,
including their cost on ring-free input. Removal of the optional bond buffer and bookkeeping was
not benchmarked here and is no longer proposed. These are measurements of the old full-record prototypes, not the agreed finalization design.
The CX example is multicenter annotation, not a bond-index-remapping stress benchmark.

### Review findings and their disposition

- Correct recording/raising does not require changing the bond order. Opening order can retain
  current ids; closing order removes reservation and CX translation. Opening order is now selected
  for its boundary meaning; the performance difference does not drive that choice.
- The former S1 green gate was invalid: changing order before updating the frame reader changes
  existing opening-center stereo readings. That ordering-first plan is withdrawn. The stereo fix
  can be implemented with current storage order independently of builder streamlining. Benchmarks
  belong at the first parser change.
- RingClosure and RingClosureEnd are public open records. The design must specify what raising does
  with invalid bond/atom references, mismatched endpoints, duplicate end positions, out-of-range
  positions, or absent required frame information. Parser-established consistency does not establish
  a contract for arbitrary TableIR construction. No additional wrapper type is implied.
- With multiple ring ends, reconstruct the frame by removing all ring neighbors and inserting the
  ends in increasing written position. Unspecified successive moves can disturb previously placed
  positions. The assertion that *every* closing center is inverted is too broad: parity depends on
  the number/order of closing digits. `N1CCC2CCC[C@H]12` gives the same canonical result before and
  after correction; moving the preceding bond past two closing bonds is an even permutation.

The all-ring design used in this impact review and the later public sparse-boundary design are both
superseded by the full-frame finalization design above. The implementation plan separates the
approved baseline work from the revised prototype and production implementation, which have not started.

Use separate Cargo target directories when comparing source variants: shared targets reused a stale
benchmark from the other source tree until the package was explicitly rebuilt. The reported original
binary was rebuilt and its distinct digest checked before the final timing pass.

## S0a baseline results — 2026-09-09

The local experimental harness pins original revision
`be941902c7ce4ca3e1a2fbe7994bad25dd2f5334`, the earlier corrected-opening patch,
its own dependency lockfile, the 516-input corpus,
and 39 focused inputs with literal proposed frames and directional-stereo controls. Sources and
Cargo targets are isolated per variant. The corrected-opening baseline is the earlier all-ring-record,
per-atom-count prototype, **not** the sparse complete-frame implementation planned in S0b.

Both variants pass all five reference-test instances. The tests check independent frame construction,
checked-constructor errors, expected-frame chemical equivalence and mirror distinction, and existing
bounded directional-bond equivalence with opposite E/Z and conjugated EE/EZ controls. They do not
assert that a baseline parser already produces the proposed TableIR records. Missing-site and duplicate-site
errors are tested through independently assembled graph-IR entries; malformed source frames have
literal arity/duplicate-ligand expectations for the new path.

The new corpus captures reproduce both historical captures byte for byte. Comparing original with
corrected-opening again gives **37 changed raised outputs, 35 changed resolved outputs, and 32 changed
canonical resolved outputs**. Ordered and sorted bond tables, acceptance/error results, and paired
fingerprint digests are unchanged. Canonicalization still uses Nauty with para-stereo disabled and
skips inputs above 40 atoms; this is a measurement bound, not a parser restriction. The corpus is a
pinned diagnostic sample, not evidence that all downstream operations are invariant.

Focused captures additionally read resolved stereo in each independently specified ligand frame.
The original disagrees on the closing pair, reversed mixed-digit order, large-label and reused-label
closing cases, and the later dot-component root. Corrected-opening fixes those five closing cases;
it still disagrees on `C.[C@H](F)(Cl)Br`, because it retains the old global-index virtual-slot heuristic.
That remaining root correction is already covered by the approved finalization design. The two-closing
case agrees in both baselines; a closing digit does not necessarily imply an odd frame permutation.

The captures retain exact syntax, raise and ingestion failures. Both baselines currently reject
`[C@]` and `[C@H2](F)Cl` with TetrahedralLigandCount; the new reference entries instead require
StereoLigandArity and DuplicateStereoLigand, respectively. Directional dangling/conflict and mismatched
ring-direction diagnostics, unsupported raw stereo and TH range errors are retained for S0b comparison.
Explicit-reference raised EDN omits the old path's source `#T` constraints. Full resolution already
discharges satisfied assertions through Resolver::plan_discharge, even though the stereo sub-pass
does not reset them by default. Resolved raw-output differences can instead reflect different
ligand frames with transported cosets; compare canonical results or read in a common frame.

The harness and raw captures are temporary experimental artifacts; this section records the findings
without depending on local artifact links. Measurements used Rust 1.96.0, release builds on
macOS 15.7.3/aarch64, dependency default features, and no Python extension. Ingestion used
ChemistryModel defaults with ValenceModel::smiles() and ResolveConfig::default().
Timings cover basic and extended
parsing through existing finalization, raising from a reused table, and full ingestion. Long-prefix
and stereo-dense inputs omit ingestion timing; the CX control measures extended parsing only. Each
operation has seven 50 ms samples in each of two passes, with variant order reversed for the second
pass and no concurrent builds. Timing includes result destruction and per-invocation clock checks.

Selected pass medians, in microseconds (ranges span the two passes):

| Case / operation | Original | Corrected-opening | Change across passes |
| --- | ---: | ---: | ---: |
| Ring-free chain10 / parse | 0.434–0.443 | 0.534–0.538 | +21.6–23.2% |
| Closing center / parse | 0.380–0.387 | 0.472–0.480 | +22.1–26.5% |
| Mixed digits / parse | 0.472–0.475 | 0.597–0.601 | +25.5–27.4% |
| Long prefix / parse | 14.374–14.391 | 15.417–15.508 | +7.3–7.8% |
| Ring-free chain10 / ingest | 20.125–20.288 | 20.496–20.638 | +1.0–2.5% |
| Closing center / ingest | 17.424–17.514 | 17.424–17.556 | −0.5–+0.8% |

Across the selected basic-parse cases, corrected-opening costs about 7–27% more. Across the seven
measured ingestion cases, changes range from −1.1% to +2.5%. This preserves the distinction between
parser overhead and its downstream share of total work. The ring-free overhead is consistent with
the old prototype's shared bookkeeping, but these measurements do not isolate allocations from
counters or predict S0b's finalization cost. They do not reopen the semantic choice of opening-order
bond storage. The two-pass ranges describe observed variation, not statistical confidence bounds.

Production source and production test expectations remain unchanged. S0a is complete; the subsequent
S0b results are recorded below. Production stages have not started; the S0 approval gate remains
open for review of the revised prototype results.

## S0b full-frame prototype results — 2026-09-09

The isolated experiment implements the settled TableIR-owned StereoAtom/StereoLigand/Winding
shape, sparse private closing/root records, complete finalization, checked direct raising, and
removal of ChiralityFrame. Production code and existing production test expectations are unchanged.
Experimental code, test additions, captures, patch and measurements remain local experimental
artifacts. This section records the reviewable findings without local artifact links.

### Implementation and correctness

The builder creates an empty output record when it reads a supported descriptor. Ordinary atoms do
not get counters or records. Only affected closing digits add private records, and only marked
traversal roots add root ids. Finalization walks the bond table once, finds relevant output rows by
binary search, and delays affected incidences until their recorded boundaries. A temporary index
sorted by reserved bond slot supports those skips while the original closing records preserve digit
order at equal boundaries. Virtual ligands are then inserted in the specified root/incoming slot.
There is no new all-atom neighbor table in parsing. The existing raising neighbor table is unchanged.

The measured bond pass and closure-index construction cost O(B log(S + 1) + B log(R + 1) +
R log(R + 1)) for B bonds, S marked sites, and R affected closures. Virtual insertion additionally
visits the S output rows and looks up their marked-root ids. Emitted ligand storage and private
records scale with S and R. Raising uses a linear explicit-site membership search per atom in this
prototype, which is another O(A S) cost to revisit if dense stereo cases matter. No new public
constructor, validator, source-format flag, or graph-IR type was added.

Fourteen release test instances pass. The 39 S0a focused inputs are joined by three independently
specified cases covering two marked roots, a ring with both ends marked, and two openings at one
marked center. Tests compare exact basic/extended rows with literal frames, preserve rows through
conversion and reaction parsing, compare checked raising with independently assembled entries,
and exercise missing sites/ligands, non-neighbor ligands, duplicate sites, malformed arity and
repeated virtual hydrogen. Explicit rows override wedge interpretation even when malformed; raw
chirality changes do not affect them. Disabling the tetrahedral model still produces the existing
stereo-entity failure. Every successful resolved focused reference agrees in its specified frame.

The branch-return fixture exposed an implementation error during development: Bond normalizes
endpoint order, so the closing site cannot be assumed to be the second stored endpoint. The final
prototype uses the recorded site to identify its opposite neighbor and to delay the correct
incidence. The fixture passes without changing its expectation.

MOL controls retain raw V2000 parity and emit no explicit rows. Basic and extended raising agree;
parity-only, definite-wedge, unknown-wedge and conflicting-wedge raised/resolved results match both
baselines exactly. A V3000 control receives the same existing InvalidCountsLine rejection in all
variants; it does not establish V3000 support. Directional-bond equivalence and opposite E/Z controls
continue to pass through the unchanged mechanism. General permutation suites remain outside scope;
production bounded permutation coverage is still S3.

### Downstream differences and diagnostics

On the same pinned 516-case corpus, comparing the final prototype with either baseline gives:

| Observation | Versus original | Versus corrected-opening |
| --- | ---: | ---: |
| Changed ordered or sorted bond tables | 0 | 0 |
| Changed raised outputs | 101 | 101 |
| Changed resolved outputs | 96 | 96 |
| Changed canonical resolved outputs | 32 | 0 |
| Changed parse/raise/ingestion acceptance | 0 | 0 |
| Changed paired fingerprint digests | 0 | 0 |
| Changed ingestion diagnostics | 1 | 1 |

The wider raised-output changes are the direct entity representation replacing source constraints.
The resolved differences reflect retained source ligand frames and their equivalent cosets instead
of the old resolver's chosen frames. Satisfied `#T` constraints were already discharged by full
resolution; their retention was incorrectly stated in the S0a write-up and is corrected above.
Canonicalization still uses Nauty, disables para-stereo and skips inputs above 40 atoms. The counts
support limited downstream impact on this sample; they are not a general invariance guarantee.

The sole changed corpus diagnostic is `[H]C(C)=[C@]=C([H])I`: the two-ligand tetrahedral descriptor
now fails checked raising with StereoLigandArity (expected 4, actual 2), replacing
TetrahedralLigandCount (count 2). Both paths reject it. The separate later-component-root fixture
now has the correct winding; both older baselines disagreed there. No successful focused frame
mismatch remains.

Additional controls confirm the intentional failure-boundary changes:

- `[C@H2](F)Cl` fails with DuplicateStereoLigand rather than the old TetrahedralLigandCount.
- `[C@H](F)(Cl)(Br)I` fails checked raising with arity 5; it cannot reach model resolution with
  an oversized explicit frame.
- `[C@](F)(Cl)Br` and `[B@](F)(Cl)Br` retain the declared lone-pair ligand through raising, then
  fail with StereoAtomFailure under the default model/policy. The parser does not infer a lone-pair
  count or erase the descriptor to make either input acceptable.
- Error/Keep/Remove now act through `stereo_atom_failure` for these explicit entities.
  Changing only `tetrahedral_stereo_failure` no longer changes their outcome. Keep retains the
  entity; Remove removes it. These are the existing policies, not new resolver behavior.

### Complete-path cost

The timing harness, lockfile, toolchain/platform, model and canonicalization bound remain as in S0a.
All three variants were measured anew after compilation, in original/corrected/prototype order and
then the reverse order. Each pass has seven 50 ms samples per operation. The allocation probe is a
separate binary and was not used for timing. Table entries are pass-median ranges in microseconds:

| Case / operation | Original | Corrected-opening | Full-frame prototype |
| --- | ---: | ---: | ---: |
| Ring-free chain10 / parse | 0.445–0.460 | 0.518–0.536 | 0.434–0.448 |
| Unmarked ring6 / parse | 0.396–0.399 | 0.472–0.481 | 0.388–0.403 |
| Acyclic stereo / parse | 0.274 | 0.295–0.300 | 0.343–0.345 |
| Opening center / parse | 0.394–0.399 | 0.459–0.468 | 0.457–0.478 |
| Closing center / parse | 0.382–0.383 | 0.468–0.476 | 0.504–0.701 |
| Mixed digits / parse | 0.476–0.479 | 0.581–0.611 | 0.606–0.612 |
| Long prefix / parse | 14.161–14.331 | 15.441–15.455 | 17.209–17.531 |
| 10 stereo components / parse | 1.664–1.682 | 2.133–2.160 | 2.753–2.804 |
| 100 stereo components / parse | 12.487–12.507 | 14.915–15.096 | 28.727–28.802 |
| Acyclic stereo / ingest | 11.120–11.217 | 10.751–11.042 | 11.196–11.530 |
| Opening center / ingest | 17.388–17.390 | 17.329–17.492 | 17.339–18.057 |
| Closing center / ingest | 17.335–17.480 | 17.836–18.255 | 17.427–17.656 |
| Mixed digits / ingest | 21.637–22.016 | 21.993–22.042 | 22.274–22.290 |

Unmarked chain/ring parsing stays within about 3% of the original in these passes, rather than paying
the old prototype's bookkeeping allocation. Complete frames nevertheless add real cost at marked
sites. Small stereo inputs generally add about 15–29% to basic parsing, while the closing case ranges
from +32% to +84% and has substantial pass-to-pass variation. The synthetic 100-component case is
about 2.3 times the original. Neither that input nor the long-prefix case establishes typical scale.
Extended parsing of the dense case costs about 60–62% more; the small marked cases cost about 10–24%
more. Raising the small stereo inputs costs about 10–24% more as explicit entries pass integrity
checking. Full ingestion across the seven measured inputs ranges from −2.7% to +3.8%; ingestion of
the long-prefix and dense inputs was not measured. The CX extended-parse control differs by about
+2%. These are observations with run variation, not confidence intervals or a prediction of production
workloads; the unstable closing-case timing warrants attention in subsequent performance work.

Allocation/reallocation call counts for basic parsing provide a separate structural check:

| Case | Original | Corrected-opening | Full-frame prototype |
| --- | ---: | ---: | ---: |
| Ring-free chain10 | 4 / 2 | 5 / 4 | 4 / 2 |
| Unmarked ring6 | 5 / 1 | 7 / 2 | 5 / 1 |
| Acyclic stereo | 5 / 0 | 6 / 0 | 7 / 0 |
| Closing center | 5 / 1 | 7 / 2 | 9 / 1 |
| 100 stereo components | 5 / 8 | 7 / 21 | 108 / 18 |

The prototype adds no allocation calls on the unmarked chain/ring controls. Their total requested
bytes grow by a fixed 96 bytes in this builder lifecycle because the molecule record is larger.
Marked-site allocations include the output row vector and each frame's ligand vector; affected
closures additionally allocate private records and the finalization index. The dense case therefore
exposes genuine per-frame allocation and search costs, not an all-atom counter vector. The probe
counts allocation requests, not peak or retained memory, and does not isolate CPU time by allocation.

S0b is complete. The prototype supports the agreed semantics and isolates the expected diagnostic and
policy changes; performance is not uniformly free. Review the marked/dense-case costs before approving
production work. At S0b closeout S1–S4 had not started; subsequent authorization and progress are
recorded below. The separate opening-order builder cleanup remains S4.

## S1a completion — 2026-09-09

The user approved proceeding with S1a after the S0b results. Production TableIR now defines
StereoAtom, StereoLigand and Winding in its existing stereo module, with the agreed public fields,
u32 indices, variants and derives. The existing TableIR wildcard export exposes them; no additional
export or constructor was introduced. Each new type has a one-line rustdoc summary; detailed
comments and examples were removed at the user's request.

All 3,422 umol-io library tests passed for S1a. Library Clippy with warnings denied,
workspace formatting and diff whitespace checks pass. Storage, conversions and checked raising
remain S1b. Parser emission and flag retirement remain S2; neither current parsing nor
raising behavior changes in S1a. The full S1 conformance/property gate remains due at the end of S1.

## S1b completion — 2026-09-09

Both TableIR molecule types now store stereo_atoms and preserve the records through conversions.
Raising preserves supplied ligand order and winding, anchors virtual ligands at their site, and
publishes through the existing checked graph-IR construction. Explicit records take precedence
over raw chirality and wedges; malformed records return integrity errors without falling back.

Conversion, frame, malformed-record and MOL controls pass, as does the full
`cargo test -p umol-io --features conformance,proptest` gate. All-target Clippy with those features
and warnings denied, formatting, and diff whitespace checks pass. Existing parser expectations
required no changes. Parser frame emission and ChiralityFrame removal remain S2.

## S2a completion — 2026-09-09

Both builders now retain marked traversal-root ids and private ring-closure tuples
(atom, reserved bond slot, bond-table length). Only supported tetrahedral closing sites add
closure records; equal boundaries preserve digit order. The parser supplies root status from
its preceding-atom state, including after dots. Ordinary atoms, unmarked closures and opening-only
sites add no closure records. No per-atom vector or ordinary-bond accumulator was added.
Records are cleared at component finalization; emitting full frames remains S2b.

The initial 28 new builder cases and the full IO conformance/property gate passed, with no existing
expectation changes. All-target Clippy with those features and warnings denied, formatting,
and diff whitespace checks pass. Opening slots, diagnostics, direction/donation handling,
CX completion ranks and spans use their existing paths.

The existing smiles_parsing benchmarks were measured against S1 commit
`bc5ab5247052fdb1f9a9dfddfa777f3d37492383`, using Rust 1.96.0, the default features and release
bench profile. No builds or tests ran concurrently with measurements. The baseline and first
comparison used 20 samples, 0.2 s warmup and 0.5 s measurement; a repeat used 40 samples,
0.5 s warmup and 1 s measurement. Criterion's repeat point estimates relative to baseline were:

| Control | Basic | Extended |
| --- | ---: | ---: |
| chain/c_10 | +0.5% | +1.8% |
| rings/c6 | −2.2% | +0.4% |
| ring_stereo/dir_up_open | +0.7% | +1.3% |
| brackets/brkt_C_50 | +6.6% | +0.6% |

The basic bracket-chain slowdown persisted in both comparisons (+6.8% and +6.6%); the initial
basic directed-ring slowdown (+5.0%) did not. These unmarked controls allocate no new records,
so the bracket-chain result is not a record-allocation cost.

Follow-up isolation identified an inlining change: the root-recording branch makes basic
MoleculeEditor::on_atom an out-of-line call at every bracket atom. S1 inlines it; removing only
root recording restores inlining, while removing only closure recording does not. Extended
on_atom was already out of line. Two reverse-order rounds (50 samples, 0.5 s warmup, 2 s measurement)
measured S1 at 1.34–1.35 µs, S2a at 1.41–1.44 µs, root recording removed at 1.35–1.36 µs, and
closure recording removed at 1.39–1.41 µs. A further isolated experiment retained all S2a logic
and forced basic on_atom to inline: 1.32–1.33 µs versus 1.41 µs for S2a between those runs.
This supports changed inlining as the dominant cost, not sparse-record allocation. Forced inlining
remains experimental; production retains the existing inline annotation. Full-finalization
measurements remain due in S2b.

The approved follow-up moves the supported-root check into both bracket-atom parsing branches,
with a small on_stereo_root method retaining private storage. on_atom has its original signature
and body; the boolean parameter is gone. Generated code confirms basic on_atom and the root push
inline without forcing the annotation. Alternating the updated and saved S1 executables with the
same 50-sample settings measured 1.340–1.357 µs for the update and 1.339–1.340 µs for S1, recovering
the earlier slowdown. The remaining 14 builder cases cover closure recording and component cleanup;
full IO conformance/property tests and all-target Clippy pass again without expectation changes.

## Implementation plan

S0a and S0b are complete, with their results above. The design contract is specified above.
S0 establishes measurements and exercises the revised prototype;
S1 and S2a are complete. S2b–S4 have not started. No mutating git operation
or commit is implied by the reviewable subitems below.

Stages end green, including required caller/expectation migrations. Each code subitem includes its
own focused tests; a public breaking change and the migration restoring green stay in the same stage.
The temporary coexistence of old and new paths in S1 is only to sequence the migration, not a permanent
compatibility surface. Do not add source-format checks as a replacement for the retiring flag.

### S0 — Baselines and revised prototype

- **S0a — Reference cases and measurements (completed 2026-09-09)** — additive (green), experimental harness. Retain the
  original/corrected-opening baselines and encode the explicit-frame success/error cases specified
  above, together with the bounded directional-bond equivalence and E/Z controls. Benchmark production
  parsing and ingestion before changes, including finalization cost. Pin source revisions, corpus,
  feature/config choices, and canonicalization limits. Use an isolated source copy and separate Cargo
  target; no production test expectations are rewritten. **[dep: none]**
- **S0b — Revised prototype (completed 2026-09-09)** — additive (green in the production tree), isolated experiment.
  Implement the specified private ring/root bookkeeping, public entry shape, finalization rules, and
  direct checked raising. Exercise MOL controls without activating parity. Measure the full path and
  classify differences from both baselines, including malformed frames and existing entity-failure
  policies. Record the measurements and observed behavior in this document for production approval.
  **[dep: S0a]**

**Gate:** the prototype implements the design above, passes its independently specified cases, and
has reviewed full-finalization timings and acceptance/diagnostic differences. A contradiction in the
specified design returns to design review; it is not permission to choose another contract inside an
implementation subitem. S1 requires production approval, not merely a compiling prototype.

### S1 — Add the explicit-frame boundary and checked consumer

- **S1a — TableIR stereo records (completed 2026-09-09)** — additive (green). Add StereoAtom, StereoLigand, and Winding
  in `table_ir::stereo` with the specified fields, variants, derives, and exports. Document the
  tetrahedral viewing convention, implicit virtual anchors, and open-record contract. Add exact
  representation cases without adding constructors or graph-IR expression machinery.
  **[dep: S0b approval]**
- **S1b — Storage, conversions, and direct raise (completed 2026-09-09)** — breaking (red→green within the subitem).
  Add Vec<StereoAtom> to Molecule and ExtendedMolecule; migrate empty/conversion paths and complete
  literals. In `table_ir::raise`, convert indices, ligands, and winding to graph-IR entries, preserving
  frame order, then use the existing checked publication boundary. Apply the explicit site precedence
  rule without indexing unchecked references. Add exact conversion/frame-preservation and malformed
  frame error cases, plus MOL wedge/parity controls. Keep the old parser-produced descriptor path
  working until S2; the new field must not be writable yet silently ignored by raising. Do not
  duplicate a frame as a `#T` constraint. **[dep: S1a]**

**Gate:** TableIR conversion/raising tests and `cargo test -p umol-io --features conformance,proptest`;
existing parser behavior stays unchanged, while independently supplied explicit frames are supported.
Graph-IR construction and its invariants remain unchanged.

### S2 — Finalize SMILES frames and retire the interpretation flag

- **S2a — Sparse encounter state (completed 2026-09-09)** — additive (green). Update the basic and extended builders'
  ring-closing paths to retain only the specified private records, and record marked traversal roots
  when there is no preceding atom. Preserve opening slots, pending ring diagnostics, direction/donation
  handling, CX completion ranks, and source spans. Carry focused record-order/root tests and benchmark ordinary/ring-free controls at this change. **[dep: S1b]**
- **S2b — Complete SMILES finalization** — breaking (red→green within the subitem). Finalize full
  frames and associated configurations for supported descriptors in basic, extended, and reaction-side
  paths; discard private records. Migrate parser construction utilities and exact expected values
  together. Add independently specified ligand-frame cases, the equivalent/mirror regression, and
  all specified failure/unsupported-input cases. Migrate downstream exact output expectations that
  necessarily change at this switch, with a recorded explanation per class. **[dep: S2a]**
- **S2c — Remove ChiralityFrame** — breaking (red→green within the subitem). Remove the enum,
  molecule fields, conversions, parser assignments, and old source-descriptor interpretation branch;
  migrate all imports, literals, and tests. Remove helpers used solely for the retired SMILES frame
  reconstruction. Retain helpers still used by wedges or directional stereo, raw MOL parity, and
  ConfigurationScope. Verify that raw parity with no operative frame remains unread and that no
  replacement source-format dispatch was introduced. **[dep: S2b]**
- **S2d — Ingestion and known fixture corrections** — additive (green) with required expectation
  migrations. In `umol-graph::ingest` and its resolution fixtures, exercise direct-entity resolution
  and the known equivalent/mirror cases. Regenerate the three identified incorrect isomer fixtures
  only after independently checking their intended identities. Review `verify_stereo` output changes
  as evidence, not authorization for wholesale fixture regeneration. Include existing stereo failure
  and mismatch policy controls; report newly exposed differences without widening resolver scope.
  **[dep: S2b, S2c]**

**Gate:** all affected tests pass, including `cargo test -p umol-io --features conformance,proptest`
and `cargo test -p umol-graph --features conformance`. If a fixture migration is required for S2b
to be green, perform that migration in S2b rather than deferring a failing test to S2d. Re-run the
parse/finalization, raise, and ingestion measurements. Record explicit-frame output movement apart
from changes in stereo meaning; do not advertise the old 37/32 counts as the new expected totals.

### S3 — Limited permutation coverage and downstream verification

- **S3a — Boundary equivalence laws** — additive (green). Add the bounded source-presentation and
  bond-storage transformations specified above beside the owning `umol-io` parser/raise tests.
  Reuse existing comparison APIs and literal expected frames; do not widen production visibility or
  add general generators to support tests. Include opening/closing centers, virtual ligands,
  multiple ring digits, conversion preservation, and mirror-negative controls. Add the separate
  stereo-bond source, storage, and endpoint transformations above, including opposite-E/Z and
  under-specified/error controls; test directional interpretation rather than treating it as an atom
  frame. **[dep: S2d]**
- **S3b — Ingestion and Python coverage** — additive (green). Exercise the same small semantic cases
  through molecule and reaction ingestion in `umol-graph`, and the supported Python molecule/reaction
  surface. Verify explicit atom frames survive resolution and directional-bond cases retain the
  expected E/Z distinction through the existing `#C` resolution path. Verify errors reach the existing
  boundary and MOL behavior stays unchanged. Use fingerprints only as supplementary observations. Build the extension
  from this source before Python checks. **[dep: S3a]**
- **S3c — Impact review and documentation** — additive (green). Run the final gates below, repeat the
  bounded output comparison, and retain timings including finalization. Update the relevant public
  rustdoc and this record with actual outcomes and limits. Add a narrow dated qualification to doc 104
  about the old bond-list-frame claim; reconcile doc 153 T4's targeted scope and doc 217's link to the
  completed correction. Keep unrelated tasks and unmeasured workloads explicit. **[dep: S3b]**

**Gate:** the core stereo correction and bounded coverage are complete and green, with every observed
semantic/acceptance difference explained. A performance or semantic concern returns to review; broad
snapshot replacement is not the acceptance criterion.

### S4 — Separate builder streamlining (deferrable)

- **S4a — Finalization lifecycle** — breaking private rewire (red→green within the subitem), with
  behavior-preserving results. Review the duplicated builders' `molecules` vector and sole
  `on_component_end` call from `finish`; direct consuming finalization may replace the one-element
  vector and subsequent `take_ring_bonds` call. Migrate parser wrappers and test utilities together,
  preserving empty input, dot components, reaction sides, spans/errors, frames, and CX indices.
  Measure against the S3 implementation so cleanup cost is separated from stereo correction.
  **[dep: S3c]**
- **S4b — Cleanup verification and closeout** — additive (green). Check exact output identity against
  the corrected pre-cleanup implementation, run applicable gates, and record the timing comparison.
  Do not remove `bond_table`, `ring_table`, `closed_bonds`, or the CX map merely because another
  computation could reproduce them. No shared-builder framework or other parser refactor is implied.
  **[dep: S4a]**

**Gate:** exact behavior preservation and green tests. The core fix does not depend on S4. If cleanup
is deferred, record that scope disposition explicitly before marking the overall discussion completed.

### Verification commands and critical path

Use narrow checks while iterating, then these final gates:

```sh
cargo fmt --all
cargo test -p umol-io --features conformance,proptest
cargo test -p umol-graph --features conformance,proptest
source umol-py/.venv/bin/activate && python --version
source umol-py/.venv/bin/activate && cargo test --workspace
source umol-py/.venv/bin/activate && cargo clippy --workspace --all-targets -- -D warnings
source umol-py/.venv/bin/activate && maturin develop --manifest-path umol-py/Cargo.toml
source umol-py/.venv/bin/activate && pytest -q umol-py/tests
```

Confirm the activated interpreter is Python 3.13 before any workspace/PyO3 compilation. Run the
selected `smiles_parsing` benchmark cases and the revised S0 raise/ingestion harness after builds have
finished. Feature-gated IO and graph property/conformance targets are explicit above; the default
workspace test alone does not cover them. Existing property suites are regression gates, not an
expansion into a new randomized-testing project. Additional graph-IR or other feature-specific gates
are required only if approved implementation changes actually touch those components.

Critical path: **S0 baseline/prototype approval → S1 checked boundary → S2 parser/flag migration →
S3 bounded coverage and downstream review**. S4 is independent cleanup scheduled afterward to isolate
its effects. S0, S1 and S2a are complete. Full frame finalization and subsequent
subitems remain pending.
