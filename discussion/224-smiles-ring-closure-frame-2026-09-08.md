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

The architectural exploration is complete. The cursor architecture with pending opening-order bond
slots is selected for both basic and extended parsing, preserving source-order semantics for upcoming
SMILES roundtripping. Retain the selected pending finalizer and vector reuse. The measured extended
cost and retained-capacity limitations are accepted for this integration; closing order is a diagnostic
counterexample, not an implementation alternative. S2b6 fixtures and baseline verification are
complete; S2b7 cursor/pending integration is next, before the remaining flag removal and limited
permutation work. No further tuning round is queued, and production integration has not started.

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

## S2b completion — 2026-09-09

Basic, extended and reaction-side SMILES now emit complete StereoAtom records for supported
tetrahedral descriptors. Recording stays in the bracket parsing branch; on_atom retains its original
body and optimized inlining. Finalization recovers actual ligand order from opening-order bonds and
private closing records, inserts typed virtual ligands, preserves winding, and discards private
state. Unmarked molecules skip finalization. Self-loops retain both source incidences; parallel
incidences and repeated virtual ligands are preserved for checked raising rather than repaired.

Thirty-one new parser cases cover exact frames, roots, ring-digit ordering, virtual ligands,
unsupported descriptors, malformed topology, both representations, conversions and reaction sides.
The equivalent/mirror ingestion regression passes. Existing lexical range-error tests remain green.
The wildcard literal now includes its incomplete explicit frame; SMILES-only `#T` helper cases were
replaced by frame coverage. Aggregate raising expects an entity without a duplicate constraint, and
repeated-hydrogen/arity failures expect MoleculeEntries integrity errors. Three graph ingestion
expectations changed from TetrahedralLigandCount to wrapped StereoLigandArity. No snapshot or isomer
fixture regeneration was required. MOL wedge/parity and directional-stereo paths are unchanged.

Both gates pass: `cargo test -p umol-io --features conformance,proptest` (3,469 library tests,
10,032 SMILES, 2,253 MOL, 407 SDF, six layout and six property tests) and
`cargo test -p umol-graph --features conformance` (984 library tests and all integration suites,
including 683 resolution cases). All-target Clippy with the same features and warnings denied,
formatting and diff whitespace checks pass.

Against S2a commit `57ddf339433ae2055fa5c51d642507b42a414472`, the same 516-input census changes
101 raised, 96 resolved and 32 canonical outputs. Bond tables, acceptance, fingerprints and size-skip
counts are unchanged. The sole diagnostic change is the existing two-ligand tetrahedral allene input
now reaching the frame-arity gate. Successful focused results agree with their independently
specified canonical references. The existing entity-versus-constraint policy controls reproduce
S0b's documented differences; resolver policy was not changed. Canonicalization retains the existing
40-atom bound and context.

Release measurements used Rust 1.96.0 and the same eleven inputs and seven 50 ms samples per operation
as S0b. Two before/after comparisons ran without concurrent builds or tests. Small marked-input basic
parsing increased 16–48%, extended parsing 12–19%, and raising 9–19%; small marked-input ingestion
ranged from −9.4% to +5.3%. Timing variation also affected unchanged controls (up to +13% ingestion),
so these are observed ranges rather than isolated causal estimates. The 100-center input increased
107–121% in basic parsing, 53–60% in extended parsing, and 21–23% in raising, consistent with the dense
stereo concern recorded in S0b. Full frame storage and the bond pass are included in these timings.

The separate alternating bracket-chain check used 40 Criterion samples, 0.5 s warmup and 1 s
measurement. Basic parsing measured 1.350–1.365 µs versus S2a's 1.335–1.351 µs; extended parsing
measured 2.269–2.280 µs versus 2.315–2.341 µs. The earlier lost-inlining regression has not returned.
ChiralityFrame and the remaining legacy interpretation branch are retained for removal in S2c;
no new public types or signatures were introduced in S2b.

## S2b1 completion — 2026-09-09

Allocation capture and Rust source attribution are verified in the existing Debian ARM64 VM.
The probe parses `N[C@@H](C)C(=O)O` 1,000 times through Smiles::parse, dropping each result.
It uses the current working-tree parser, Rust 1.90.0, and an optimized release build with full
debug information. This is a tool-capability check, not the S2b3 performance baseline.

The macOS installation provides cargo-heaptrack 0.1.0, but not the separate heaptrack executable
that it launches. For Linux, Heaptrack 1.4.0 and its missing dependencies were unpacked into the
experimental area without changing installed Debian packages; the same Cargo wrapper version was
built there. Experimental source, scripts, lockfile, and reports remain in `scratch/`.

The complete Linux workflow succeeds: `cargo heaptrack`, raw-data interpretation with
`heaptrack_interpret`, and report generation with `heaptrack_print`. This wrapper defaults to raw
capture, so interpretation is a separate required step while the matching executable and debug
information are available. The report resolves inline/caller stacks to parser source, including
MoleculeEditor::on_stereo_atom at `builder.rs:172`, finish_stereo at `builder.rs:64`, and bond-table
collection at `builder.rs:403`. It records 8,010 allocation-function calls for the whole process;
the on_stereo_atom allocation stack accounts for 1,000 calls. Process totals include runtime
startup, and instrumented runtime/RSS are not clean parser timing or per-parse memory measurements.

Native macOS sampling and CPU-counter capture are also verified. Instruments Allocations no longer
crashes at startup after the Xcode update, but still fails to attach with both 2-second and 10-second
recording windows. Its cause remains undiagnosed; no additional attachment investigation is needed
for this subitem because the Linux allocation path supplies the required attribution. S2b2 establishes
the workload set before drawing optimization conclusions.

## S2b2 completion — 2026-09-09

The profiling set contains 244 inputs: 52 controlled cases and 192 real-source inputs. Generation,
validation, frozen inputs/output digests, and the operation protocol are experimental artifacts in
`scratch/`. Selection is deterministic and independent of parser acceptance:

- 64 rows from `materials/formats/opensmiles/examples/ZINC.FL.smi`, selected by lowest SHA-256 of
  the complete source line. The source contains 37,661 rows and no stereo markers; selected inputs
  have 7–21 atoms.
- 64 unique atom-stereo participants and 64 unique directional-bond participants from
  `materials/atom_mapping/rhea_reactions.csv`, selected by lowest SHA-256 of participant SMILES.
  Directional markers take precedence for grouping; first occurrence supplies reaction/row provenance.
  The source pools contain 5,811 and 4,067 unique participants respectively. Selected inputs have
  6–130 and 6–150 atoms, with up to 29 and 30 marked atom sites respectively.

These are coverage samples, not estimates of production frequencies. Report each source group
separately with equal input weights and individual results; do not mix synthetic controls into
real-source throughput. Source-file hashes, row identifiers, working-tree source digest, Cargo
lockfile digest, and compiler/profile configuration are pinned in the manifest. The input TSV has
SHA-256 `400120aba6257ab8cf24122dddce1bae068feb1a4134b300f31cb8edc071b93e`.

Controlled families cover 16/64/256-atom chains and bracket density, 64-atom chain/comb branching,
serial/nested rings with 1/2/4/8 closures, and equivalent ring labels 1/9/%10/%99. Stereo-density
families hold 20-site or 100-site chain topology, atom/bond counts, brackets, and branches fixed,
varying only which sites carry @; selected front/spread variants vary placement. Six small controls
cover internal/root stereo atoms, ring-opening/closing sites, and cis/trans directional bonds.
The independent structural counts pass. Ring topology and bracket spelling controls are not claims
of equivalent resolved chemistry.

All 244 inputs parse and raise successfully through both paths on Linux ARM64 with Rust 1.90.0,
default features, and optimized release code with full debug information. Raised graph IR agrees
for 243 inputs. The retained Rhea input `[1*][C@@H]1O[C@@H]2COP(=O)([O-])O[C@H]2[C@H]1[2*]`
gains NotAromatic constraints on its two isotope-labelled wildcard atoms through the extended path.
It is recorded separately with both outputs, not silently excluded or repaired here. Output digests
are movement checks, not independent semantic or canonicalization oracles.

S2e revisits the basic/extended discrepancy for isotope-labelled wildcard atoms after the profiling
and flag migration work. Neither output is designated correct by the current comparison.

Four operations are specified: basic parse, extended parse, basic parse plus raise, and extended
parse plus checked conversion to basic TableIR plus raise. Each iteration includes destruction of
its result and intermediate tables. Input loading/generation, validation, hashing, resolution, and
canonicalization are outside measurement. Mac Rust 1.96.0 and Linux Rust 1.90.0 remain separate
baselines; toolchains stay fixed within comparisons. No timings or optimization conclusions are
claimed by this subitem. S2b3 implements measurement against this frozen workload set.

## S2b3 completion — 2026-09-09

The frozen 244-input set now has clean native timings, scoped allocation measurements, eight native
CPU profiles, four Linux Heaptrack captures, and two native CPU-counter captures. Experimental
runners, binaries, metadata, and full per-input results remain in `scratch/`. Production code is
unchanged. Mac validation exactly matches the S2b2 Linux outputs, including the separate wildcard
case. The input digest remains the one recorded in S2b2.

Timings use the M1 Pro, macOS 15.7.3, Rust 1.96.0, default features, and optimized release code with
debug information. Seven rotated rounds use calibrated approximately 20 ms batches; clock reads
surround batches, and output/intermediate destruction is included. Builds and other profiling runs
do not overlap timing. The table gives median batch-average microseconds per input, with equal input
weights within each source group. The wildcard discrepancy is measured individually and excluded
from the 63-input atom-stereo group comparison.

| Source group | Inputs | Basic parse | Extended parse | Basic parse + raise | Extended parse + conversion + raise |
| --- | ---: | ---: | ---: | ---: | ---: |
| ZINC | 64 | 0.698 | 1.051 | 5.235 | 6.767 |
| Rhea atom stereo | 63 | 1.642 | 2.271 | 12.735 | 15.960 |
| Rhea directional bonds | 64 | 1.667 | 2.329 | 13.749 | 17.207 |

Extended parsing is 38–51% slower across these groups. This compares representation paths, not
implementation revisions. The median full sample range across 988 case/operation rows is 4.6%;
30 rows exceed a 20% range. Group min/max ranges are about 4–6% of their medians. These are observed
ranges, not confidence intervals, and small individual differences require focused remeasurement.

The fixed-topology controls separate molecule growth from increasing marked-site count:

| Available sites | Marked sites | Basic parse, µs | Basic parse + raise, µs | Parse allocations / reallocations | Parse peak live bytes |
| --- | ---: | ---: | ---: | --- | ---: |
| 20 | 0 | 1.289 | 9.397 | 5 / 4 | 23,960 |
| 20 | 1 | 1.471 | 10.316 | 7 / 4 | 24,264 |
| 20 | 10 | 1.944 | 15.864 | 16 / 6 | 26,232 |
| 20 | 20 | 2.343 | 21.916 | 26 / 7 | 28,504 |
| 100 | 0 | 4.674 | 38.132 | 5 / 6 | 112,280 |
| 100 | 10 | 6.219 | 44.926 | 16 / 8 | 114,552 |
| 100 | 100 | 9.886 | 102.153 | 106 / 11 | 133,976 |

Ten marked sites on the 20-site chain add 0.655 µs to basic parsing and 6.467 µs to the combined
operation. Fully marking the 100-site chain adds 5.212 µs to parsing. Atom/bond counts stay fixed;
adding markers still changes input byte length and can affect capacity estimates. These differences
include the cost of publishing additional stereo information, not just avoidable bookkeeping.

All 976 scoped allocation records agree exactly between Mac Rust 1.96.0 and Linux Rust 1.90.0.
Each records one warmed operation through destruction, checks zero remaining measured live bytes,
and balances allocations against deallocations. Reallocations are counted separately. Gross requested
bytes include the new size of each successful allocation/reallocation; peak is logical requested
live Rust heap size, not RSS or the allocator's internal/transient memory.

| Source group | Basic parse allocations / reallocations, mean | Basic gross / peak, KiB | Extended gross / peak, KiB |
| --- | --- | --- | --- |
| ZINC | 5.84 / 2.48 | 6.87 / 6.09 | 16.95 / 15.09 |
| Rhea atom stereo | 13.22 / 3.92 | 17.15 / 15.37 | 42.90 / 38.83 |
| Rhea directional bonds | 10.84 / 4.19 | 18.35 / 16.20 | 46.18 / 41.20 |

Extended parsing makes the same allocation/reallocation counts for these inputs but requests about
2.5 times the bytes. Basic parse plus raise averages 72.16, 143.63, and 158.64 allocation/reallocation
calls respectively. Heaptrack resolves parser allocation stacks to buffer construction, stereo-frame
ligand growth, and final table collection. The three parse-only captures reproduce the scoped count
totals plus the same 7,622 whole-process setup calls. Combined-path stacks also show allocations in
Molecule::try_from_entries and hash-table construction. Heaptrack process peaks include manifest
loading; they are not used as per-operation peaks. Rust GlobalAlloc counts exclude direct native
library allocations. No Linux elapsed-time comparison is claimed.

Samply profiles contain 3,960–4,000 samples each in the 0.5–4.5 s window of five-second runs.
Allocator-library leaf shares are 42.0% for basic ZINC parsing, 29.9% for Rhea atom-stereo parsing,
and 29.6% for Rhea directional-bond parsing. MoleculeEditor::finish has 16.6% of leaf samples in the
atom-stereo group and 28.7% in the fully marked 100-site control. This includes inlined finalization
work; it is not an isolated measurement of finish_stereo. Extended atom-stereo parsing has 7.2%
of leaf samples in memmove. Both combined paths expose hashing and graph-construction costs beyond
the parser. Full function/library distributions are retained for S2b4 attribution.

Instruments CPU Counters captured the 20-site controls with zero and ten markers. Cycle-weighted
instruction-bandwidth categories in the 1–4.5 s window are delivery 8.2%/11.8%, discarded 2.1%/2.9%,
processing 11.7%/8.6%, and useful 77.9%/76.6%. Execution is almost entirely on performance cores.
These are one capture per input, not elapsed-time percentages or per-parse cycle counts; they guide
follow-up investigation rather than establish the cause of the timing difference.

S2b4 now explains which costs are necessary and which work can be removed. No optimization or
acceptance of the previously recorded implementation regressions is implied by this baseline.

## S2b4 completion — 2026-09-09

Code inspection, layout/capacity measurements, lookup-count replay, and 20 bounded Callgrind captures
explain the main S2b3 costs without production changes. Native timings remain those from S2b3;
instruction counts below describe the Linux ARM64/Rust 1.90.0 binary, not M1 cycles or estimated
microsecond savings. Experimental sources and reports remain in `scratch/`.

**Scanning and table construction.** Both inner parsers pass input byte length as atom capacity
and byte length minus one as pending-bond capacity. Atom rows occupy 104 bytes in basic TableIR and
288 bytes in extended TableIR; bond rows occupy 40 and 96 bytes. Their Option<Bond> counterparts
have the same sizes, so reserved opening slots add no per-row padding in these builds. Retaining
opening-order semantics is not the source of that basic/extended size difference.

The 64-atom chain reserves 64 atom slots; the fully bracketed 64-atom chain reserves 320. Its basic
parse peak rises from 12,632 to 49,496 requested live bytes with the same four allocations and four
reallocations. In the real-source groups, aggregate reserved/used atom-slot ratios are 1.86 for ZINC,
2.76 for Rhea atom stereo, and 2.43 for Rhea directional bonds. Atom capacity survives in the output.
These are capacity costs, not extra parsed atoms. Large extended rows amplify reservation, row
initialization, movement, and destruction; allocation count alone misses this effect.

Bracket recognition is a separate cost. The chain/bracket control increases from 17,890 to 29,664
instructions per parse; the bracket parser accounts for 10,761 exclusive instructions in the latter.
Native time rises from 1.169 to 1.713 µs. This supports treating token recognition separately from
capacity management rather than attributing the entire difference to allocation.

**Final table construction.** The parser validates that all ring slots are closed, then drains the
optional bond table through flatten/collect into a new bond vector. For the 64-atom chain this creates
the final bond allocation and four growth reallocations despite the pending table already having
capacity for every bond. A temporary molecules vector adds another allocation solely to return one
molecule, which the wrapper immediately pops. These are ownership/lifecycle choices, distinct from
the required atom and bond outputs. The already-planned S4 lifecycle work touches this area; no
representation or iterator rewrite has been made as part of attribution.

**Ring handling.** The ring table grows to the largest label and is scanned for unmatched rings at
the end. Relabelling the same 64-atom ring from 1 to %99 increases peak requested memory by 4,416
bytes and instructions from 18,297 to 19,656; native time changes from 1.224 to 1.297 µs. Of the extra
bytes, 576 follow from the longer input's atom/pending-bond reservations and 3,840 from ring-table
capacity. Eight serial or nested closures cost 20,790/21,079 instructions and 1.400/1.502 µs. Ring
matching and diagnostics are required, but storage proportional to label range is an implementation
choice. These cases do not make ring handling the dominant explanation for stereo regressions:
the density controls contain no rings or stereo-closure records.

**Stereo finalization.** With no marked sites the extra bond pass is skipped. With any sites,
finish_stereo scans every bond and searches the sorted parser-private frame list for both endpoints.
Closure-record lookup/sorting adds work when such records exist. The ring-free controls isolate the
site searches and per-frame output construction:

| Available / marked sites | Endpoint-search comparisons | Total instructions/parse | Exclusive builder-finish instructions/parse |
| --- | ---: | ---: | ---: |
| 20 / 0 | 0 | 18,630 | 414 |
| 20 / 1 | 82 | 21,001 | 2,085 |
| 20 / 10 | 410 | 28,923 | 5,838 |
| 20 / 20 | 492 | 36,407 | 7,466 |
| 100 / 0 | 0 | 75,669 | 1,374 |
| 100 / 1 | 402 | 84,275 | 9,285 |
| 100 / 10 | 2,010 | 105,892 | 24,558 |
| 100 / 100 | 3,216 | 177,723 | 43,422 |

On the 202-atom control, adding one marked site adds 8,606 instructions, of which 7,911 are in
builder finalization itself. This is predominantly the newly activated whole-bond pass, not a ring
record or per-atom counter allocation. With 100 marked sites, allocator exclusive instructions also
grow from 3,264 to 42,157. Every site's ligand Vec allocates separately; the current tetrahedral
frames fit its initial capacity of four, so implicit-H insertion does not add ligand reallocations
in these controls. Complete frames and their ordering are required; the whole-bond search pass and
one heap buffer per frame are implementation/storage choices. Changing the latter requires review
of the settled TableIR storage contract, not an incidental parser patch.

**Raising and graph publication.** The combined operation has additional costs beyond parsing:

- TableIR raising always builds AtomNeighbors as a Vec of per-atom Vecs. On the connected
  multi-atom controls, this adds an outer allocation and one neighbor buffer per atom before output construction.
- raise_tetrahedral_stereo linearly searches the explicit-frame list for every atom to enforce
  explicit-frame precedence. The 20-site/10-marker and 100-site/100-marker controls perform 375 and
  15,250 comparisons respectively. This is O(A*S) work separate from the parser's O(B log S) searches.
  Public TableIR frame order is not guaranteed sorted, so replacing this with an unchecked binary
  search would be incorrect. Legacy interpretation bookkeeping also remains until S2c.
- Raising materializes graph-IR ligand Vecs; relation construction then copies them into flat
  participant storage and builds incidence indexes. The integrity pass materializes ligand_frame
  again and uses small HashSets to check ligand and participant uniqueness, plus sets for bond/site
  uniqueness. Hashing and allocation in the profiles have concrete consumers here. Atom constraints
  themselves use inline SmallVec storage, not per-atom hash maps.
- The consuming ExtendedMolecule-to-Molecule conversion nevertheless iterates by reference and
  clones atom/bond rows, stereo frames, and other carried fields before dropping the extended input.
  Ownership transfer is therefore not currently used to avoid these copies.

The required contract is faithful conversion and representation integrity, including references,
frame arity, uniqueness, and incidence. The current temporary collections, repeated materialization,
and hash-set mechanisms are not themselves semantic requirements. S2b5 must distinguish parser-local
work from changes owned by TableIR conversion, graph-IR publication, or graph-core relation storage;
none of this authorizes skipping checks or widening the wildcard correction.

**Tooling qualification.** The initial function-filtered Callgrind capture was rejected after
inconsistent attribution. Valgrind's explicit start/stop client requests now delimit the batch. Inclusive call attribution
still produces impossible percentages with this installed toolchain, even with inline-info reading
disabled, so only exclusive counts are used and their sum is checked against the reported total.
Zero/100/200 iterations of the 20-site/10-marker case record 23/2,892,324/5,781,924 instructions;
per-call totals agree within 0.05% between the nonempty runs. Builder-finish totals include its inlined
work but exclude allocator callees. They are not interchangeable with native CPU-time percentages.

S2b5 can now rank opportunities and estimate realistic savings from these measured costs. This
subitem selects no optimization implementation and changes no public contract.

## S2b5 completion — 2026-09-09

Historical analysis: the ranking below is no longer the active work queue. See the holistic
accumulator design before the implementation plan.

The ranked review is complete; experiment selection remains open. No production code changed in
S2b1–S2b5, and the S2b regressions remain unresolved. The ranking below prioritizes explaining and
reducing that regression, then broader savings without changing the settled representation.

**Savings scale.** Basic parsing averages 0.698–1.667 µs in the three real-source groups; their
combined parse/raise times are 5.235–13.749 µs. Several microseconds cannot be removed from their
average parse alone. Larger individual inputs have more room: the 202-atom, 100-marker control takes
9.886 µs to parse and 102.153 µs including raising. These selected groups do not establish a
production workload distribution.

Multiplying native sample shares by clean timings gives a rough scale for the work under review,
not predicted savings or a strict bound: sampling and timing used separate runs, and changing one
operation can affect others. Builder-finish leaf work corresponds to about 0.27 µs per Rhea
atom-stereo input, 0.30 µs in the 20-site/10-marker control, and 2.84 µs in the 100-marker control.
All allocator-library leaf work corresponds to approximately 0.29–0.49 µs per real-source parse and
1.96 µs in the 100-marker control. Only part of either cost is removable. Candidate budgets overlap
and must not be added; Callgrind instruction shares are not substituted for time shares.

| Rank | Opportunity and affected paths | Evidence and plausible scale | Scope and risk |
| --- | --- | --- | --- |
| 1 | Reduce stereo finalization searches in both parsers | One marker activates a full bond pass; on the 202-atom control, finalization accounts for 7,911 of the additional 8,606 instructions. Dense controls also repeat endpoint searches. The measured finish work is a few tenths of a microsecond on smaller cases and several microseconds on the dense control. | Private parser algorithm, medium complexity. Preserve exact ligand order, opening slots, delayed closing ligands, roots, and virtual ligands. Avoid adding state or allocations to unmarked atoms. |
| 2 | Remove avoidable finalization allocations and growth | The 64-atom chain grows its final bond Vec four times, then allocates a vector to return one molecule. Benefits apply even without stereo. Expect fractions of a microsecond on the real-source means, within the shared allocator cost above; larger-input movement needs measurement. | Small ownership/lifecycle changes, relatively low semantic risk. Overlaps S4; selecting an earlier experiment does not silently reschedule production S4. Preserve CX records and error behavior. |
| 3 | Use ownership in extended-to-basic conversion | The consuming conversion clones rows and carried fields. Rhea atom-stereo combined time exceeds the basic path by 3.225 µs, of which 0.629 µs is the parse-only difference. The remaining roughly 2.6 µs is a budget to investigate, not an isolated conversion timing or expected saving. | Bounded TableIR conversion experiment, low algorithmic complexity. Move owned fields where possible; different row representations still require conversion. Preserve failure behavior and the separately tracked wildcard discrepancy. |
| 4 | Reduce temporary collections in raising and graph publication | Neighbor buffers allocate per atom; explicit-frame precedence costs O(A*S); frame materialization and small uniqueness HashSets add work. Rhea atom-stereo combined time exceeds parsing by roughly 11.1 µs. Allocator leaf work alone scales to about 3.2 µs of the combined operation. This is the strongest broader route to multi-microsecond savings, but no individual mechanism has that entire budget. | Separate experiments in their owning modules, medium to high scope. Public frames need not be sorted. Preserve wedge/directional interpretation and all checked-publication invariants; no unchecked shortcut. |
| 5 | Reduce byte-length-based overreservation | Real-source atom capacity is 1.86–2.76 times used slots; bracketed basic chain peak is 49,496 versus 12,632 bytes. Extended rows amplify the memory cost. A memory reduction is well motivated; elapsed savings are unquantified and extra growth or a counting pass could make time worse. | Parser capacity policy, medium performance risk. Measure retained capacity, peak bytes, and reallocations together. No extra full scan or new general allocator machinery is justified yet. |
| 6 | Reduce per-frame ligand allocations | The dense control has one ligand buffer per marked site and high allocator cost. It shares the 1.96 µs dense-parser allocator budget with other allocations. The smaller real-source means offer much less absolute room. | Higher design cost: changing public Vec storage reopens the settled TableIR contract. Review only after private algorithm and ownership experiments establish what remains. |

For an end-to-end optimization project, rank 4 has more absolute headroom than ranks 1–3. It is
ranked later here because the immediate issue is the parser regression and its scope is broader.
The 20-site/10-marker case remains a required control; a dense-case improvement cannot excuse its
regression or hide it in a corpus mean.

Bracket scanning and ring-label storage are lower priorities for this objective. Full bracketing
adds 0.544 µs to the 64-atom basic control, including both recognition and representation costs;
changing ring label 1 to %99 adds 0.073 µs. Neither currently supports a broad scanner rewrite or
changing opening-order storage. These measurements also do not justify treating all extra time for
marked stereo as removable: the corrected parser produces additional required frame data.

**Proposed first experiments, for review.** Start with ranks 1 and 2 as independent candidates.
For rank 1, first compare the existing endpoint lookup with a sparse lookup strategy over marked
sites while retaining the same replay and output representation. This isolates lookup cost from a
larger rewrite; it may leave the one-marker full-scan cost largely unchanged. Any subsequent design
that gathers ligands during parsing must account explicitly for ring slots and prove that ordinary
atoms incur no new retained state. For rank 2, measure final-bond reservation separately from direct
single-molecule return so their benefits are distinguishable. Do not assume a consuming iterator
will reuse an allocation without checking generated behavior and allocation counts.

Rank 3 concerns the explicit conversion operation only. ExtendedMolecule should be processed in
its extended representation; requiring conversion to basic Molecule defeats its purpose. Removing
that architectural requirement is distinct from making an explicitly requested conversion efficient.
The measured extended parse/conversion/raise lane characterizes the conversion path, not the intended
general path for extended data. Its timing does not establish the cost of direct extended processing.

For the conversion experiment, first time conversion of prebuilt extended tables with setup excluded,
then compare ownership transfer and the measured combined path. Rank 4 needs similarly separate
measurements of raising and checked graph publication before choosing between neighbor storage,
precedence lookup, or tiny-set validation. The broader extended-processing defect is recorded here
for separate design work; this review does not authorize an implementation or add it to the stereo
correction. No common collection framework or public storage change is proposed by this review.

**Comparison protocol.** Keep candidates isolated under `scratch/`, with pinned inputs, source,
compiler, features, and allocator matching a freshly measured control. Run clean interleaved A/B
batches; report absolute microseconds and individual sparse/dense controls alongside group means.
Measure allocation count, growth, requested bytes, and retained capacity separately from timing.
The S2b3 median sample range of 4.6% is insufficient to establish small gains: use longer repeated
focused measurements where a candidate's effect is comparable to that variation. Reject conclusions
that depend on an instrumented timing or disappear when control/candidate order changes.

Before interpreting speedups, compare exact successful outputs and acceptance/error behavior,
including roots, ring-opening/closing centers, virtual ligands, directional bonds, reaction sides,
and CX indexing. Preserve the known wildcard difference until S2e. Promotion requires the owning
subitem's correctness gates and an explicit review of any regression; aggregate speed alone is
insufficient. Keep S2c flag removal separate from optimization comparisons, rebasing and remeasuring
both sides if its implementation precedes an experiment. Experiment selection and any resulting
production-plan change remain for review before S2c.

## Independent rank-1 and rank-2 experiments — 2026-09-09

The first two ranks were authorized for isolated experiments. Four candidates were compared with a
fresh control copied from the current S2b implementation; none includes another candidate's changes:

- **Rank 1, map:** replace frame binary searches with a standard HashMap from marked atom index to
  frame index. Keep bond replay, closing records, and ligand storage unchanged.
- **Rank 1, existing chirality:** skip endpoint searches when the existing atom chirality is absent.
  Keep binary searches for potentially marked endpoints and the same replay. This adds no atom field,
  accumulator, or allocation; unsupported descriptors still pass through the existing frame lookup.
- **Rank 2, reservation:** allocate the final bond Vec with pending-bond count as capacity before
  extending it. Keep the pending table and one-element result vector.
- **Rank 2, direct return:** return Option<Molecule> or Option<ExtendedMolecule> instead of allocating
  the one-element result vector. Keep bond collection and separate CX ring-record extraction.
  This isolates result allocation; it is not the full consuming-builder design contemplated in S4.

**Validation and measurement.** Each candidate preserves all 244 pinned validation records, including
basic/extended tables, raised outputs, and the known wildcard discrepancy. Each passes 3,469 unit
tests, 10,032 SMILES conformance cases, and six property tests with seed 224. No expected output was
rewritten. The copied direct-return builder tests only explicitly discard the changed return value.
All 976 scoped allocation records for the control match S2b3; the existing-chirality candidate's
records match the control exactly.

Native clean timing uses Rust 1.96.0, release optimization with debug information, the same dependency
lock across candidates, and the S2b2 inputs. Fourteen controls and the three real-source groups are
measured for both parsers; five workloads also exercise each combined path, making 44 operation/input
groups. Seven rounds rotate workload and variant order, reversing variant order on alternate rounds.
Each measurement has three approximately 80 ms batches following calibration, with destruction
included and input loading excluded. Builds, tests, and allocation instrumentation are separate.
Sources, patches, binary hashes, full samples, and allocation reports remain experimental in `scratch/`.

The tables report control medians of round medians and median paired candidate-minus-control changes,
in microseconds. Negative changes are faster. Paired changes need not equal differences between
independently calculated medians. These are measurements, not confidence intervals.

| Basic parse workload | Control, µs | Map Δ | Existing chirality Δ | Reservation Δ | Direct return Δ |
| --- | ---: | ---: | ---: | ---: | ---: |
| ZINC mean | 0.714 | +0.004 | +0.010 | −0.098 | −0.055 |
| Rhea atom-stereo mean | 1.648 | +0.646 | −0.061 | −0.144 | −0.048 |
| Rhea directional-bond mean | 1.665 | +0.497 | −0.005 | −0.172 | −0.040 |
| 20 sites / 10 markers | 2.010 | +0.778 | −0.112 | −0.224 | −0.086 |
| 100 sites / 0 markers | 4.657 | −0.212 | −0.230 | +0.163 | −0.021 |
| 100 sites / 10 markers | 6.238 | +2.191 | −0.787 | −0.231 | −0.159 |
| 100 sites / 100 markers | 10.108 | +5.049 | −0.336 | −0.305 | −0.079 |

| Extended parse workload | Control, µs | Map Δ | Existing chirality Δ | Reservation Δ | Direct return Δ |
| --- | ---: | ---: | ---: | ---: | ---: |
| ZINC mean | 1.053 | −0.014 | −0.003 | −0.161 | −0.079 |
| Rhea atom-stereo mean | 2.300 | +0.634 | −0.076 | −0.180 | −0.079 |
| Rhea directional-bond mean | 2.370 | +0.509 | −0.065 | −0.201 | −0.075 |
| 20 sites / 10 markers | 2.718 | +0.767 | −0.076 | −0.285 | −0.095 |
| 100 sites / 0 markers | 7.852 | +0.075 | +0.017 | −0.461 | −0.129 |
| 100 sites / 10 markers | 9.379 | +2.444 | −0.552 | −0.486 | −0.130 |
| 100 sites / 100 markers | 13.689 | +4.877 | −0.389 | −0.368 | +0.108 |

**Rank 1 outcome.** Reject the standard-map candidate: construction, allocation, and hashed lookups
cost substantially more than the binary searches they replace. It adds one allocation per marked
parse, including 280 requested bytes at ten markers and 2,184 at 100 markers. This rejects this
lookup mechanism, not the broader opportunity to improve finalization.

Reusing chirality is more promising. The basic 100-site/10-marker saving is 0.787 µs (12.6%), with
all seven paired changes between −0.881 and −0.686 µs. The shorter ten-marker control saves 0.112 µs
(5.6%); a separate nine-round ABBA/BAAB comparison confirms a median saving of 0.101 µs, with every
round faster. The full bond scan and per-frame buffers remain, so this is a limited lookup improvement.
Small changes also occur on unmarked inputs, where the changed replay is not executed: the basic
64-atom chain median increases 0.025 µs. The focused ZINC comparison changes by +0.008 µs, with
paired values from −0.011 to +0.033 µs. These effects are not attributed to the new predicate;
generated-code effects and measurement variation have not been separated.

**Rank 2 outcome.** Reservation removes final-bond growth: the basic 64-atom chain goes from four
allocations/four reallocations to four/zero, gross requested bytes from 15,032 to 12,592, and peak
live bytes from 12,632 to 12,592. At 100 markers, total reallocations fall from eleven to five;
remaining growth belongs to other buffers. Across the real-source means it saves 0.098–0.172 µs in
basic parsing and 0.161–0.201 µs in extended parsing. Direct return removes one allocation and 896
requested bytes for basic output or 1,856 for extended output; its group-mean savings are smaller.

Neither allocation reduction guarantees a timing improvement for every input. Reservation's basic
100-site/unmarked slowdown was +0.163 µs in the main run. A separate nine-round ABBA/BAAB check
measures +0.062 µs, with all paired changes positive (approximately +0.001 to +0.160 µs). Direct
return's extended 100-marker slowdown repeats at +0.097 µs median, but with a wider range crossing
zero (−0.145 to +0.684 µs; seven of nine pairs slower). These countercases remain unattributed;
the allocation improvements do not justify declaring either candidate uniformly faster.

**Combined paths and disposition.** Reservation's basic parse-plus-raise group changes are −0.111,
−0.140, and −0.095 µs for ZINC, Rhea atom stereo, and Rhea directional bonds. The extended combined
lane still includes explicit conversion to basic TableIR. Several combined measurements have large
outliers: one dense basic control produced paired differences of roughly −40 µs across candidates.
Do not interpret small dense end-to-end differences as established gains. All raw samples are retained;
no outlier is silently removed from these medians.

The existing-chirality check and bond reservation merit further review; the map does not. Direct
return has a clear allocation benefit but smaller and mixed timing effects. No candidates were
combined, no production change was promoted, and these comparisons do not remeasure against S2a or
resolve the original regression. Before promotion, review the unmarked and dense countercases and
whether the measured benefits justify each change. A larger rank-1 design would need to address the
remaining bond pass, not merely substitute another lookup container. Rank 3 and later opportunities,
S2c, and production S4 remain unimplemented by this experiment.

## Holistic accumulator design — draft, 2026-09-09

This is the active design discussion. The S2b5 ranking and its rank-1/rank-2 experiments are historical
measurements, not the work queue. The task is to choose the parser's complete accumulation and
publication model before another implementation or optimization experiment. The remaining plan below
needs reconciliation after this design is settled; its separate S4 cleanup boundary does not constrain
this review.

**Ownership.** Use one private assembly algorithm with statically specialized basic and extended row
construction. The input loop recognizes tokens and advances the byte cursor; the assembler owns
structural state and all updates caused by a token. Basic atoms/bonds and extended atoms/bonds are
constructed directly in their respective representations. No conversion through basic Molecule,
runtime format dispatch, stored token stream, event bus, or public builder interface is needed.
Sharing the assembly algorithm must not force ordinary atoms through the larger bracket/stereo path.

The scanner remains responsible for lexical spelling, lexical range errors, input termination, and
source positions. Parentheses, groups, pending bonds, dot transitions, and ring matching belong to
one structural state owner. Moving code into methods alone is not an optimization; the point is to
stop scattering coordinated writes across token branches and subsequent reconstruction passes.

**What TableIR requires.** The inventory separates output data from the temporary references needed
to assemble it:

| Output | When it becomes available | Accumulation and publication |
| --- | --- | --- |
| Atom rows: element/symbol, aromaticity, isotope, charge, implicit H, raw descriptor, class, span | Atom token | Append directly to the target atom Vec. Preserve the basic/extended distinctions, including the separately tracked wildcard issue. No generic graph object or duplicate atom table. |
| Bond rows: endpoints, order, direction/donation, span; defaults for other fields | Ordinary atom attachment; ring closure for a previously reserved row | Reserve opening-order slots. Complete ordinary rows immediately and ring rows when both specifications are known. Normalize endpoint-relative direction/donation exactly once using the existing bond construction rules. |
| Explicit stereo frames | Source order is known incrementally; some endpoints remain unknown until ring closure | Retain ordered bond-slot indices only for marked sites; finalization converts those indices into full StereoAtom records. |
| Positions, labels, values, radicals, lone-pair annotations, wedges, cis/trans annotations, coordinate/hydrogen bonds, multicenter bonds | CX suffix | Apply the existing checked CX updates after body completion. These do not need accumulators in the atom loop. |
| Extended query fields, atom properties, ligand order, link nodes, SGroups, fragment groups, enhanced/bicyclo stereo, configuration scope | Extended CX suffix | Retain the extended annotation phase and its own grouping state. Its LigandOrder field is distinct from the explicit tetrahedral frame. |
| Molecule/reaction comments and properties, CTfile-only data, other absent fields | Not supplied by the current SMILES path | Keep existing empty/absent values; create no speculative collectors. Source format is known by the boundary, with existing empty-input behavior preserved. |
| Reaction sides and atom mapping | Body boundaries and atom-class tokens | Each side owns its local atom/bond namespace. The reaction wrapper owns the final mapping and CX distribution across sides. |

Directional bonds remain bond data. This design does not add eager E/Z perception to parsing or
transfer the raise's graph-integrity responsibilities into the assembler. MOL parsing and ignored
parity remain separate. The target retains the agreed removal of ChiralityFrame, not another source
interpretation flag.

**Live traversal state.** Replace the current atom index with a small private cursor containing its
atom index, aromatic status for implicit-bond selection, and an optional pending-frame index.
Branch frames save and restore that cursor together with their existing grammar/error state. Dot
and group transitions follow the current grammar; a dot resets the current cursor, not the ring table
or completed output. Reaction-side termination ends the whole body and its local namespace.

There is no cursor vector indexed by every atom. Cursors live only at the current position and in
the branch stack. A marked atom's frame index is assigned when its bracket is parsed and follows
that cursor through branch returns. This eliminates both binary and hash lookups for stereo sites.
The additional cursor fields do enlarge live branch entries; that cost must be measured, not described
as free.

**Accumulators.** The proposed storage is:

| State | Contents and lifetime |
| --- | --- |
| Atom Vec | Final target rows, moved into the result. |
| Bond storage | Selected: Vec<Option<Bond>> or Vec<Option<ExtendedBond>>, with unresolved physical opening slots. Sparse rings remain a measured alternative, not the active design. |
| Current cursor, branch stack, pending bond specification, group state | Live structural parsing state. No molecule-wide adjacency or per-atom counter array. |
| Open-ring table and open-ring count | Existing label-indexed opening records: opening atom, reserved bond slot, opening specification and diagnostic positions. Count changes only at ring opening/closure. No new stereo fields in these records. |
| Pending stereo Vec | One record per supported marked site: atom, winding, is_root, and ordered bond indices. No separate roots vector, closure replay log, site map, or atom-degree counter. |
| Optional CX ring records | Existing closing-rank/opening-slot pairs, only when CX tracking is needed; passed privately to annotation processing and then discarded. |

Use SmallVec<[usize; 4]> for each pending frame's bond indices as the initial proposed storage: four
is the supported tetrahedral degree, not an assumed molecule size. Larger lists spill rather than
truncate, preserving malformed input for the existing raise error. This adds no per-frame heap buffer
for ordinary valid tetrahedral input before publication; the public Vec<StereoLigand> still needs its
own allocation. Larger pending records and overlap with final output can increase peak memory, so
inline storage is a design candidate to measure, not a claimed allocation-free parser.

**One incidence operation for ordinary bonds and rings.** A pending stereo record contains bond
slots, not partly filled public ligands or a ring-specific placeholder variant:

- On ordinary attachment, append the new bond slot to the previous atom's frame if marked and to
  the new atom's frame if marked. For the new atom this is its first actual incidence.
- On a ring opening digit, reserve a bond slot and append its index to the current atom's frame if
  marked. The unresolved endpoint remains in the bond accumulator, which already has to represent it.
- On a ring closing digit, complete that same bond slot and append its index to the current closing
  atom's frame if marked. The opening frame needs no update: it already names the stable slot.

Thus both ends receive the same bond index at their respective source positions. A self-loop inserts
the same index twice in one frame; parallel bonds remain separate entries that can yield repeated
neighbors. Neither is silently repaired. Bond-slot indices stay valid across Vec reallocations and
ring-label reuse; no pointers into frame or bond buffers are retained.

For the existing case `C[C@H]1CCCCO1`, the center records bond slots `[0, 1, 2]`: incoming bond, ring
opening, then continuation. Closing the ring fills slot 1 without moving it. Finalization reads those
three bonds to obtain actual ligands `[0, 6, 2]`, then emits `[Atom(0), ImplicitHydrogen, Atom(6),
Atom(2)]`. It never visits the other bonds to find the center's neighbors. A closing center appends
the reserved slot when its closing digit is read, giving the other required source ordering naturally.

**Callback boundary.** These are direct parser-to-assembler operations, not allocated event objects:

| Input operation | Complete assembler responsibility |
| --- | --- |
| Organic atom or wildcard | Append the target row; attach to the current cursor using the pending or implicit bond; record that bond only for a marked parent; advance the cursor and update branch/group occupancy. Keep this path small. |
| Bracket atom | The same attachment operation, plus classify its supported descriptor once, create its pending frame if needed, and record its incoming incidence. The lexer does not separately call on_stereo_atom. |
| Bond token | Check placement and retain the pending specification and diagnostic position. This does not yet create a bond row. |
| Ring token | Consume the pending specification; open or close the label, resolve conflicting specifications, update the bond slot, the current marked frame, and optional CX record. |
| Open/close parenthesis | Maintain branch/group grammar, save or restore the complete cursor, and report existing placement/empty/unbalanced errors. Parentheses are not all ordinary branches. |
| Dot | Validate the separator and reset the current cursor while preserving body-wide rings and output. |
| Body finish | Validate outstanding syntax state in the existing error order, consume accumulators, and return one body plus any private CX records. |

Atom construction and attachment are one externally visible assembler operation even if their small
row-writing helpers remain separate internally. Do not add a second callback chain for stereo after
ordinary bond emission. Basic and extended emitters use the same structural algorithm with concrete
row constructors. Private compile-time specialization is appropriate; a general parser framework or
public row-builder trait is not proposed.

**CX completion order without a hot-path counter.** Let B be pending-bond Vec length and O the number
of currently open rings. Every slot is either a completed bond or belongs to one open ring, so the
number of completed bonds is B − O. When a ring closes, after decrementing O its zero-based completion
rank is B − O − 1. Record that rank with the opening slot if CX tracking is enabled. Ordinary bond
append needs no closed_bonds increment. Multiple closures at one atom and label reuse preserve the
same identity. Keep the label table and existing CX translation algorithm; they solve different
problems from stereo ordering.

**Reaction mapping.** The reaction wrapper can own the final BTreeMap before parsing its sides.
When a bracket atom supplies a class, its emission operation appends its local index directly to the
reactant or product list for that class. A private optional mapping destination supplies the side;
standalone molecules and agents retain the row's class without updating the reaction map. This needs
no per-atom class census, no temporary class list, and no final scan of both sides. Duplicate classes
and index order are preserved. The currently supported CX updates do not change atom classes. The
mapping remains local to the parse operation and is discarded if any later body or suffix fails.

**Consuming finalization.** Finish the body once, not at each dot:

1. Preserve the existing error priority: trailing bond, unclosed parenthesis/group, then unresolved
   ring with the current diagnostic selection. No partial result or incomplete frame is published.
2. Consume the pending bond Vec in opening order with a one-to-one mapping to final bonds. No open
   rings establishes that every slot is Some. Use this producer invariant rather than flattening
   away a missing slot. An owning exact-length iterator avoids the current drain/flatten growth;
   allocation reuse is compiler behavior to verify, not a Rust API guarantee. This remains at worst
   a linear B pass, with no stereo searches inside it. No fabricated endpoints or unsafe buffer
   reinterpretation are needed.
3. Visit only pending stereo records. For each named bond slot, read the endpoint opposite the site.
   Emit virtual hydrogens at position zero for a root or after the incoming ligand otherwise; apply
   the existing three-actual-ligand lone-pair rule when there is no implicit H. Read H metadata only
   for these sites. Construct the full public ligand Vec in order, retaining repeated/oversized
   frames for checked raising. Discard the bond-index lists and is_root bits.
4. Move atom rows and final stereo/bond vectors into one concrete Molecule or ExtendedMolecule.
   Move optional CX ring records alongside it as private parse output. Remove the molecules vector,
   on_component_end lifecycle, and separate post-finish take_ring_bonds extraction.
5. The outer wrapper parses, remaps, and applies CX annotations using final counts and the existing
   per-representation update functions; the reaction wrapper first distributes entries across its
   three sides. Drop the temporary CX maps before returning public TableIR.

This separates the final bond representation conversion from stereo construction. Stereo work is
O(L + S + H), where L is the number of recorded incidences at marked sites, S is marked-site count,
and H is emitted virtual-hydrogen count. A large molecule with one marked center does not cause a
stereo traversal of all B bonds. Ordinary emission still needs a predictable check for a marked
parent; this design does not promise zero additional instructions on unmarked input. No new general
adjacency structure is needed.

**Boundary contract.** The public symbol delta is empty: Molecule, ExtendedMolecule, StereoAtom,
StereoLigand, Winding, parser entrypoints, conversions, and Python interfaces retain their settled
surfaces. TableIR remains an open carrier; successful parsing resolves its private ring references
but does not establish graph-IR stereo integrity or chemistry validity. The checked raise remains
the first consumer enforcing frame arity, uniqueness, and incidence. Frame bond indices, cursor
indices, pending slots, and parsing flags never enter public TableIR. Internal direct indexing relies
on append-only producer-owned namespaces; input rejected by the current syntax checks must retain
its ParseError path rather than become an indexing panic. Representable malformed stereo still
reaches the existing checked-raise failure boundary. Exact frame order, winding, opening-order bonds, directional and
CX interpretation, and current failure boundaries are preservation requirements. ChiralityFrame
removal is the already-agreed separate public migration, not an extra flag in this design.

**What still needs review and measurement.** This is a proposed whole-assembler design, not an
implementation plan or a speedup claim. Capacity hints remain hints rather than per-atom counters:
current byte-length reservation can include bracket syntax, CX text, and later reaction sides and
must not be treated as an atom estimate. Keep its cost visible when comparing the complete design;
do not introduce an extra counting pass or guessed molecule-size threshold as an unnoticed detail.
Measure cursor/branch-stack size, pending-frame storage, output allocation/ownership, and generated
code for the ordinary emitter together. A small method boundary is not proof of inlining.

The prototype comparison below exercises this ownership/callback/storage design against corrected
S2b, using existing independent frame, ring, directional, CX, reaction and malformed-input coverage.
Lexical scanning and target representation remain fixed for attribution. The comparison reports
parse and combined costs, sparse and dense centers, branch depth, memory, and output/error agreement. The S2b5 ranking
is not resumed if one implementation choice needs revision.

### Greenfield comparison and current hypothesis

Opening order remains the chosen semantics. The earlier less-than-2% opening/closing timing
comparison retained the optional bond buffer and legacy bookkeeping; it did not measure the full
structural simplification available with closing order. Conversely, the S2a/S2b regression compares
two opening-order implementations, and the later stereo-density controls contain no rings. Neither
establishes that opening order itself costs twice as much. The user's willingness to accept roughly
10% overhead is a tradeoff preference, not a measured result or a prediction for this design.

Two independent design agents were given the output and ordering requirements without the current
parser, this discussion, or the profiles. One preferred dense ordinary bonds plus sparse ring
records; the other proposed publishing a completed prefix while buffering behind unresolved rings.
These are independent structural proposals, not correctness or performance validation.

The strongest hypothesis is shared across the alternatives: append final atom rows directly, carry
the optional stereo-frame index with the current/branch cursor, and record incidences only at marked
sites. Final stereo construction visits those recorded incidences, not every bond. Publish one body
by consuming its accumulators. Basic and extended rows remain separate direct targets. Lexical
scanning, source semantics, directional-bond handling, CX updates, and failure boundaries stay fixed.

Opening-order storage has at least three safe designs:

| Candidate | Assembly and finalization | Principal tradeoff |
| --- | --- | --- |
| Pending vector | Reserve physical opening slots in Vec<Option<Bond>>; fill rings and consume into final rows. | Simplest stable addressing. All bonds participate in final unwrapping; allocation reuse needs verification. |
| Ordinary rows plus sparse rings | Append complete ordinary rows to Vec<Bond>; reserve logical positions in a separate ring-row vector. Merge once at finish. | Ring-free bodies transfer their bond vector directly. Ring-containing bodies need a destination allocation and linear merge, with overlapping storage. |
| Completed prefix plus pending suffix | Append directly while no hole blocks publication; buffer behind the earliest unresolved ring and flush ready prefixes. | Can avoid buffering portions of a body, but adds routing/flush state. An early long-lived ring buffers almost the whole body. |

For the sparse-ring candidate, the next logical bond index is ordinary-row count plus ring-row count.
Opening a ring appends its logical index and an unresolved row to the ring vector; the active label
record names that entry. Closing fills it without changing either count. Ring entries remain in
opening order even when labels are reused or closures occur in reverse order. At finish, walk final
indices: take a ring row when its saved index is next, otherwise take the next ordinary row. This is
O(B + R), needs no sorting or repeated insertion, and preserves exactly the same final bond indices.
Here R counts ring bonds, not just simultaneously open labels. No incomplete public bond is fabricated.

Use the same pending stereo bond-index lists in both initial candidates. Their entries name logical
final indices, so physical slot reservation is unnecessary. Merge/unwrap bonds before materializing
frames. The independent alternative of direct ligand drafts with exact opening-slot patches is also
valid in principle, but introduces a second kind of deferred reference and extra ring metadata.
Bond-index lists reuse information already needed for opening order and keep ring completion
independent of stereo. Neither representation requires a map or counter for ordinary atoms.

For either initial candidate, B is the total logical bond count and O is the number of open labels;
B - O is the completed-bond count. Existing sparse CX closing-rank/opening-index records therefore
remain sufficient. The sparse-ring merge changes storage, not CX numbering or stereo semantics.

The proposed comparison is bounded to two complete assembler prototypes: the cursor design with
the pending vector, and the same design with ordinary rows plus sparse rings. The first establishes
what the new collection algorithm achieves with straightforward storage; the second tests whether
pending storage can be restricted to rings profitably. The prefix/suffix design is not an initial
prototype: its extra state lacks a demonstrated benefit over a single merge. The two isolated
prototypes were authorized for construction and comparison, with results below. This does not
authorize production promotion or resume the S2b5 ranking.

Compare each with corrected S2b for equal-output performance, and report S2a separately to show the
remaining total cost of the correction. Neither baseline isolates the marginal price of opening
versus closing order. Measure absolute microseconds and ratios, allocations and peak live bytes;
keep ring-free sparse/dense stereo, long-lived and overlapping rings, branch depth, CX and reaction
cases visible alongside corpus results. Required public ligand allocations remain part of the cost.
A near-2x total regression remains a problem to explain, not evidence that the semantic preference
must be expensive. Roughly 10% was an acceptance preference, not a prediction; the results below
show where the candidates approach it and where substantial overhead remains.

### Prototype comparison — 2026-09-09

Both proposed assemblers were constructed in isolated copies. Production source is unchanged.
Both retain opening-order bonds and emit the same corrected TableIR and raised outputs as current
S2b in the measured domain. The shared cursor architecture improves parsing, but the storage results
do not identify one universal winner: pending generally leads on the basic corpus samples and has a
lower transient peak for rings; sparse leads on several extended-parser cases.

**Implemented scope.** Both use one statically specialized basic/extended assembler owning the
current cursor, branch/group state, pending bond specification, ring matching, and output vectors.
Bracket emission creates sparse marked-site records and collects reaction mapping directly. Organic
atom emission stays separate. Stereo finalization consumes ordered bond-index lists, with no full-bond
stereo search, site map, separate roots vector, or closing-record replay. CX completion ranks derive
from total bond count minus open-ring count. Owning finalization returns one molecule and private CX
records. The existing source flag and basic/extended error-offset distinctions are preserved.

The candidates differ only in the builder's bond storage: physical optional slots versus complete
ordinary rows plus sparse ring rows and a linear final merge. Sparse ring records have separate
logical bond and physical ring-row indices. No public API, lexical grammar, target representation,
MOL behavior, directional-bond interpretation, or extended-to-basic conversion was changed.

**Correctness and tooling.** Each candidate passes 3,469 library tests, 10,032 SMILES conformance
cases, six existing property tests with seed 224, and nine new exact architecture cases: 13,516
cases total. The nine cases exercise reverse closures, reused labels after branch restoration,
marked roots joined across a dot, nested marked branches, CX bond remapping, reaction namespaces,
and terminal/root implicit H in incomplete frames. Both representations are checked directly.
Fourteen private builder cases were adapted to the new state; existing expected-value helpers were
changed to construct rows independently of the removed builder callbacks. Existing external test
inputs and expectations were retained. All-target Clippy with conformance/proptest and warnings
denied passes for both candidates. No snapshots were updated.

All 250 expanded workload validation records match current S2b exactly, including complete basic
and extended tables and raised results. The three original regression workloads also match, for
253 records total. General randomized/permutation coverage was not added; this experiment does not
close the remaining planned permutation work.

Timing uses the same release profile, Rust 1.96.0, lockfile, shared workspace dependencies, and
harness across isolated variants. Clean timing and allocation-counting binaries are separate.
Compilation, tests, lint, validation, and allocation counting finished before timing. The main run
covers 66 operation/workload pairs in seven rotated rounds, with three calibrated approximately
80 ms batches per pair/variant per round. Parse timing includes destruction and excludes input
loading/process startup. The original three regression inputs were measured afterward in a separate
seven-round comparison with the same binaries. Tables give median round-median microseconds;
paired changes are computed within rounds rather than from rounded table entries.

The additional S2a reference restores all umol-io source from
`57ddf339433ae2055fa5c51d642507b42a414472`, rebuilt with the same current manifest, dependencies,
lockfile, compiler and harness. Only parsing is timed against that reference. S2a does not emit the
same corrected full frames, so it measures remaining total correction cost, not equal-output speed
or the marginal price of opening order. All variants retain opening order.

**Basic parsing (µs).**

| Input | S2a reference | Current S2b | Cursor / pending | Cursor / sparse |
| --- | ---: | ---: | ---: | ---: |
| ZINC plain sample | 0.700 | 0.708 | 0.482 | 0.556 |
| Rhea atom-stereo sample | 1.104 | 1.650 | 1.099 | 1.169 |
| Rhea bond-stereo sample | 1.238 | 1.651 | 1.116 | 1.166 |
| 42 atoms / 10 marked sites | 1.251 | 1.967 | 1.446 | 1.425 |
| 202 atoms / 10 marked sites | 4.425 | 6.057 | 4.749 | 4.532 |
| 202 atoms / 100 marked sites | 4.296 | 10.001 | 7.413 | 7.264 |
| Original 10 stereo components | 1.772 | 2.782 | 1.699 | 2.020 |
| Original 100 stereo components | 13.043 | 28.844 | 14.378 | 15.827 |
| Original 1000-atom prefix | 14.172 | 17.680 | 13.059 | 14.134 |

**Extended parsing (µs).**

| Input | S2a reference | Current S2b | Cursor / pending | Cursor / sparse |
| --- | ---: | ---: | ---: | ---: |
| ZINC plain sample | 1.040 | 1.033 | 0.896 | 0.914 |
| Rhea atom-stereo sample | 1.737 | 2.262 | 1.922 | 1.823 |
| Rhea bond-stereo sample | 1.913 | 2.333 | 2.040 | 1.969 |
| 42 atoms / 10 marked sites | 1.963 | 2.660 | 2.384 | 2.178 |
| 202 atoms / 10 marked sites | 7.634 | 9.305 | 9.045 | 8.151 |
| 202 atoms / 100 marked sites | 7.836 | 13.381 | 12.020 | 10.673 |
| Original 10 stereo components | 3.020 | 4.065 | 3.327 | 3.302 |
| Original 100 stereo components | 26.641 | 42.358 | 32.824 | 30.452 |
| Original 1000-atom prefix | 32.708 | 35.823 | 36.241 | 30.096 |

**Interpretation and countercases.** On the original basic-parser cases, pending reduces ten stereo
components by 1.075 µs (38.7%) and 100 components by 14.442 µs (50.1%), with every paired round faster
than S2b. Relative to S2a, those cases are respectively 4.1% faster and 10.0% slower. The original
basic long-prefix case saves 4.688 µs versus S2b and is 7.4% faster than S2a. The prior doubling is
therefore not necessary to retain opening order on those inputs.

That does not establish a universal 10% correction cost. On the 202-atom/100-marker acyclic control,
basic pending remains about 73% above S2a and sparse about 69%; all variants retain opening order.
For the original extended ten-component case both candidates are about 9% above S2a. At 100
components pending remains 22.8% above S2a and sparse 14.2%. Required full-frame construction and
per-site ligand allocations remain; this comparison does not assign every remaining microsecond
to unavoidable output work.

Pending also has new extended-parser regressions against current S2b. Its 64-bracket-atom case is
slower in all seven rounds, paired median +0.246 µs (+8.3%); the unmarked 202-atom control is slower
in all seven, +0.801 µs (+10.5%), paired range +0.612–1.134 µs. The extended 256-atom chain is slower
in every round as well. These are counterexamples to promoting pending unchanged for both targets.
Sparse avoids these regressions in this run. Both candidates improve the sampled reaction and CX
cases; pending leads on those cases in both representations. The extended difference has not been
attributed to specific generated instructions; identical allocations on ring-free input do not
establish identical execution cost. No third implementation or forced-inlining experiment was added.

Basic parse-plus-raise improves by roughly 3–5% on the three corpus samples, substantially less than
standalone parsing. The extended combined lane shows small mixed changes and does not establish a
consistent end-to-end improvement. That lane retains the existing explicit extended-to-basic
conversion before raising; it is not a proposed architecture for consuming ExtendedMolecule.

**Allocation, peak, and retained memory.** Counts are scoped Rust allocator requests; byte values
exclude allocator overhead and are not RSS measurements. A separate post-timing probe measures
live requested bytes immediately after parsing, before dropping the result. Peak and retained
memory answer different questions:

| Basic-parser input | Metric | Current S2b | Cursor / pending | Cursor / sparse |
| --- | --- | ---: | ---: | ---: |
| 64-atom chain | Allocations / reallocations | 4 / 4 | 2 / 0 | 2 / 0 |
| 64-atom ring | Allocations / reallocations | 5 / 4 | 3 / 0 | 5 / 0 |
| 64-atom ring | Peak / retained bytes | 13,080 / 9,424 | 9,624 / 9,464 | 12,408 / 9,424 |
| 202 atoms / 100 markers | Peak / retained bytes | 133,976 / 100,944 | 126,008 / 125,816 | 126,008 / 125,816 |
| Original 100 stereo components | Peak / retained bytes | 268,200 / 204,152 | 226,344 / 226,184 | 257,320 / 194,264 |

The safe consuming optional-bond collection reuses its allocation on this compiler/build: pending
has no second bond allocation. Sparse allocates both its ring-row vector and a final merged bond
vector when rings occur. Its linear merge is consequently a real cost rather than an automatic
improvement over optional slots. No unsafe reinterpretation is used; allocation reuse is not a Rust
API guarantee.

Reuse also changes retained capacity. On the dense acyclic control, current S2b returns bond capacity
256 for 201 bonds; both candidates return capacity 801, inherited from the unchanged byte-length
reservation. Their final stereo vector has capacity 224 for 100 frames, versus current S2b's 128.
Basic retained requests increase 24.6% despite the lower peak; extended retained requests increase
from 262,848 to 318,240 bytes (21.1%). On the original 100-ring-component case, pending returns bond
capacity 1498 for 700 bonds, while sparse's merge returns capacity 700 and current S2b returns 1024.
Thus pending's lower transient peak does not imply a smaller returned value. The unchanged atom
reservation is oversized in these inputs too, but is shared by the variants.

**Disposition.** The cursor-based collection hypothesis is supported: these complete prototypes
remove substantial measured cost while preserving opening-order semantics and corrected outputs.
Pending slots are now selected for both basic and extended targets, retaining the cursor architecture.
The extended countercases and retained-capacity measurements remain evidence for finalization work;
they do not reopen the storage choice or resume the old ranked substitutions. S2c and the remaining
production implementation/permutation work have not advanced.

### Finalization design — draft

**Settled foundation.** The parser owns final atom rows, opening-ordered optional bond rows, sparse
pending stereo records, traversal/ring state, and optional CX remapping records. Finalization consumes
one body. Pending slots are selected for both targets; no sparse merge or storage selector is added.
The public symbol delta remains empty. The following publication/capacity choices were proposed
for the experiment recorded below; they have not been promoted to production.

**Validate, then discard completed traversal state.** Preserve the existing error priority: pending
bond, unclosed branch/group, unclosed ring. The existing open-ring count makes the successful ring
check constant-time; only when it is nonzero is the label table scanned to select the existing latest
opening diagnostic. A zero count follows from the opening/closure updates, not a separate new census.
After validation, release the branch stack and ring table before constructing final frames. Current
cursor, pending syntax, and the mapping destination borrow have no further role. Keep the atom rows,
bond rows, pending stereo records, source/legacy metadata, and CX records until their consumers finish.
No partial result is published on a parse error.

**Consume the bond rows.** Use the owning, one-to-one optional-row conversion in opening order.
Successful ring validation establishes that every slot is filled. Preserve the producer assertion;
do not filter missing rows or reinterpret the allocation with unsafe code. Allocation reuse is an
observed optimization, not the semantic contract. This remains a B-row conversion, separate from
stereo construction; inspect its generated checks and payload movement before claiming it is free
or prescribing a lower-level rewrite. No second bond allocation or unconditional shrinking is
proposed as the default publication path.

**Materialize only marked frames.** Let S be the number of pending frames. Allocate the destination
StereoAtom vector with requested capacity S and append records explicitly, avoiding reuse of the
larger PendingStereo allocation as an oversized public frame vector. This intentionally trades one
additional outer-vector allocation for predictable requested output storage; its time and peak cost
must be measured. With no marked frames, create no stereo allocation or traversal.

For each pending record, read its atom's H count and actual incidence count K. The virtual ligand
count is the supplied H count, or one lone pair under the existing K=3/H=0 rule. Request K plus that
virtual count for the ligand vector, and emit in three segments:

1. For a nonroot, consume and emit the first actual incidence (the incoming bond).
2. Emit the virtual hydrogens, or the applicable lone pair.
3. Consume and emit the remaining actual incidences.

A root starts at segment 2. Each actual incidence names a completed bond; emit its opposite endpoint.
The nonroot's first incidence exists by atom attachment, without an adjacency search. This removes
the prototype's per-incidence insertion-position branch and terminal-H special case. Root zero-degree
frames, nonroots with only their incoming bond, multiple virtual H, repeated incidences, self-loops,
and oversized frames preserve their supplied meanings. No fixed-four truncation or arity repair is
introduced. Unsupported raw descriptors still produce no supported frame. Each completed public
record carries only its atom, full ligand vector, and winding; temporary bond indices and is_root
are consumed. Checked raising retains responsibility for graph-IR frame integrity.

**Publish and annotate.** Move atom rows, completed bond rows, and full stereo frames into the concrete
basic or extended result. Forward only the private CX records needed by the outer annotation phase.
Reaction-side mapping is already accumulated in its wrapper-owned map; finalization does not rescan
classes. Existing CX remapping/application follows body completion, and its temporary records are
then released. Preserve current source metadata, empty-input behavior, and the separately planned
legacy-flag removal. No parser record enters public TableIR.

**Capacity policy.** The initial proposal distinguishes buffers by ownership rather than applying
blanket compaction:

| Buffer | Proposed publication policy |
| --- | --- |
| Final atom rows | Move the existing vector. |
| Pending bond rows | Consume through the safe conversion; retain allocation reuse when available. |
| Pending stereo records | Consume into a separately allocated final vector requested for S records, then release temporary storage. |
| Per-frame ligands | Allocate once with requested capacity K plus virtual count. |
| Branch/ring matching state | Release after syntax validation. |
| CX remapping records | Retain through annotation only. |

This removes the incidental oversized stereo-vector result, but does not claim to solve oversized
atom/bond reservations. The existing byte-length hints, including later reaction sides and CX text,
remain an explicit accumulation-capacity question. Do not silently substitute guessed input-size
thresholds, a counting pass, unconditional shrink_to_fit, or a public compactness option in the
finalizer. If compact atom/bond publication is required, compare its additional allocation/movement
against reuse as an explicit policy decision. Retained capacity is measured resource behavior, not a
new public representation invariant.

**Review boundary.** First review this consuming lifecycle and the proposed distinct stereo-output
allocation. A subsequent bounded prototype should preserve exact outputs/errors and measure total
parse time, finalization work, allocations, peak requested bytes, and retained requested bytes. Keep
the original stereo components, ordinary inputs, and extended countercases visible. The sparse
prototype is a reference, not another candidate to retune. No production change or new implementation
stage is authorized by this design draft.

### Finalization experiment — 2026-09-09

**Scope and correctness.** Compared the selected pending prototype with a variant changing only
Assembler::finish: count-gated ring validation, early release of branch/ring storage, three-segment
ligand emission, and a fresh outer stereo vector requested for S frames. Bond conversion and the
parse loop were unchanged. Basic and extended results remain direct targets. This comparison is
against the pending cursor prototype, not current production S2b or the older S2a reference.

The full proposal passed 3,469 library tests, nine exact architecture cases, 10,032 conformance
cases, and six property tests (seed 224): 13,516 tests with unchanged expectations. All-target
Clippy, including the conformance/property features, passed with warnings denied. All 253 validation
records matched the pending prototype, including complete basic/extended tables and raised outputs.
Malformed and unusual frames retain the existing behavior; no new repair or arity restriction was
introduced.

**Timing.** Saved clean release binaries on the same macOS machine, Rust 1.96.0; seven rotated rounds,
three calibrated approximately 80 ms batches per variant per workload per round, including result
destruction. The main run covered 72 operation/workload pairs: basic/extended parsing, selected
parse-and-raise groups, reaction/CX lanes, and the original stereo component/long-prefix cases.
There were no concurrent builds or profiling runs during timing.

The fresh-vector proposal increased basic parse time broadly, including unmarked inputs: chain64
by 0.083 µs (9.4%), the sampled Rhea atom-stereo group by 0.095 µs (8.6%), original 10 stereo
components by 0.196 µs (11.8%), and original 100 components by 1.352 µs (9.5%). Each was slower in
all seven rounds. Extended results were generally smaller and mixed; combined parse-and-raise
measurements established no consistent improvement.

A second diagnostic retained the new validation, early release, and ligand sequence, but restored
the owning map/collect conversion for the outer stereo vector. Its 1,844 parser tests, nine exact
architecture cases, library Clippy, and all 253 validation records passed. The full conformance and
property suites were not rerun for this diagnostic. It was compared with both other variants over
20 parse workload/target pairs, using the same seven-round protocol.

The following values are median round times from that three-way diagnostic, in µs. Parenthesized
changes are medians of paired round differences against pending, rather than subtraction of the
rounded displayed medians.

| Target / workload | Pending | Full proposal, fresh vector | New sequence, vector reuse |
| --- | ---: | ---: | ---: |
| Basic / chain64, no stereo | 0.874 | 0.955 (+0.080) | 0.874 (−0.000) |
| Basic / ring label 99 | 0.963 | 1.031 (+0.067) | 0.945 (−0.019) |
| Basic / 42 atoms, 10 markers | 1.423 | 1.573 (+0.144) | 1.448 (+0.021) |
| Basic / 202 atoms, no markers | 4.252 | 4.289 (+0.042) | 3.993 (−0.259) |
| Basic / 202 atoms, 100 markers | 7.354 | 8.024 (+0.649) | 7.633 (+0.272) |
| Basic / sampled Rhea atom stereo | 1.097 | 1.189 (+0.091) | 1.109 (+0.007) |
| Basic / sampled Rhea bond stereo | 1.098 | 1.168 (+0.074) | 1.099 (+0.004) |
| Basic / original 10 stereo components | 1.650 | 1.840 (+0.191) | 1.690 (+0.044) |
| Basic / original 100 stereo components | 13.981 | 15.349 (+1.369) | 14.319 (+0.322) |
| Extended / 42 atoms, 10 markers | 2.378 | 2.476 (+0.098) | 2.421 (+0.043) |
| Extended / 202 atoms, 100 markers | 12.079 | 12.283 (+0.315) | 12.281 (+0.137) |
| Extended / sampled Rhea atom stereo | 1.901 | 1.947 (+0.047) | 1.948 (+0.045) |
| Extended / original 100 stereo components | 32.078 | 32.459 (+0.336) | 32.381 (+0.185) |

Reuse removes most of the basic regression, but the remaining marked synthetic/original cases are
still slower in all seven rounds: 1.5% for 42 atoms/10 markers, 3.7% for 202 atoms/100 markers, and
2.7%/2.3% for the original 10/100 components. The ordinary Rhea basic differences are small and
mixed across rounds. The unmarked 202-atom improvement is repeatable here, but is not a general
speedup across workloads. Unmarked chain64's fresh-vector regression cannot be explained by an
extra stereo allocation: there are no frames. The outer-vector expression affects performance
beyond allocation, plausibly through generated code; this experiment does not isolate that mechanism.
Nor does it separately attribute the residual costs to validation, early release, or ligand emission.

**Memory.** These are requested Rust allocation bytes, not RSS or allocator overhead. All tracked
allocations were released after result destruction. Atom and bond capacities were unchanged.

| Basic workload | Outer stereo capacity, pending → fresh | Retained bytes, pending → fresh | Peak bytes, pending → fresh |
| --- | ---: | ---: | ---: |
| 42 atoms, 10 markers | 28 → 10 | 23,064 → 22,488 | 23,256 → 23,384 |
| 202 atoms, 100 markers | 224 → 100 | 125,816 → 121,848 | 126,008 → 129,016 |
| Original 100 stereo components | 224 → 100 | 226,184 → 222,216 | 226,344 → 229,384 |

The fresh destination adds one allocation on marked inputs. Saving 3,968 retained bytes in the
100-frame cases costs approximately 3 KB of peak requested memory, since the destination and pending
records coexist. The reuse diagnostic restores pending's retained capacities and allocation counts;
early traversal release reduces the original 100-component peak from 226,344 to 226,184 bytes.
Most reserved space remains in the unchanged atom/bond vectors, so exact outer stereo capacity does
not resolve the larger capacity question.

**Disposition.** The proposed sequence preserves tested semantics, but the complete fresh-vector
version is not supported as the default: it trades time, another allocation, and higher peak memory
for a modest retained-size reduction. Reuse is the better measured publication choice, although
these results do not establish the remaining sequence changes as an optimization. Keep the selected
pending architecture and its existing finalizer as the measured reference while reviewing whether
the simpler ligand sequence warrants its residual cost. No production changes were made. The
atom/bond capacity policy remains open; sparse bond assembly is not reopened by these results.

### Bounded closing-order counterexample and remeasurement — 2026-09-09

**Boundary.** One closing-order implementation, one static diagnostic pass over the extended
ring-free paths, and one coordinated timing run. No tuning variants, inlining experiment, new
capacity heuristic, or second search over architectures. Pending opening order remains the selected
production direction until the counterexample is reviewed. These experiments make no production
changes and do not complete the remaining implementation/permutation work.

**Counterexample.** Starting with the selected cursor/pending prototype, ordinary bonds append as
complete rows at atom attachment and ring bonds append as complete rows at closure. The assembler
has no optional bond rows, opening bond reservations, final bond conversion, open-ring count, or CX
ring-remapping vector. Ring matching remains necessary. Only a marked ring-opening atom reserves a
position in its private ligand list; its open-ring record retains the stereo-record/position pair,
and closure supplies the completed bond index. Existing frame finalization and outer-vector reuse
are unchanged. CX closure-order indices address physical bond rows directly. Basic and extended
parsing remain separate direct targets. Runtime source changes are confined to the builder.

**Correctness and visible ordering changes.** All 10,032 unchanged conformance cases and six unchanged
property tests pass (seed 224). Nine architecture cases pass with independently enumerated
closing-order bond rows and unchanged expected stereo frames. Library Clippy passes. The initial
library run has 3,263 passing tests and 206 failures:

- 194 compare complete tables. Reconstructing both assertion sides and sorting only complete bond
  rows makes every comparison equal, preserving all other fields.
- Twelve address positions directly: two private bond-index vectors, four first-bond order
  expectations, one wildcard bond vector, four cis/trans target bond indices, and one error's bond
  index. Explicit scratch expectation changes make all twelve pass. The second run has 3,275 passes
  and precisely the 194 full-table permutation failures above. The suite was not rewritten to hide
  these changes, and the unmodified library suite is not claimed to pass.

All 266 normalized validation records agree between opening and closing: complete basic/extended
tables, full frames, raised values after row normalization, reaction sides/mapping, and sampled CX
bond donations. These comprise the existing 250 workloads, three original regression inputs, five
additional ring/frame cases, and eight inputs exercising positional assertions. Normalization sorts
complete bond rows with all their fields, preserving parallel-row multiplicity. It is a representation
comparison; it does not by itself verify native-order raising. CX checks cover donation placement,
not every nested annotation that can contain a bond reference.

A separate native-order downstream check covers 254 inputs in each target (508 records): the existing
253 inputs plus the reverse-closure CX donation case. All raise successfully; 334 raw raised strings
change with bond order. For the 171 inputs with at most 40 atoms, all 342 resolution records agree,
and all 316 canonical outputs from determined resolutions agree. The remaining 26 resolutions are
underdetermined in both variants. The 83 larger inputs (166 target records) are explicitly skipped
before resolution/canonicalization, keeping the downstream check bounded; they still undergo raising
and the separate table/frame comparison. Canonicalization uses Nauty with para_stereo disabled.
Extended downstream checks use the existing explicit conversion to basic; this is a validation route,
not a proposed mandatory model conversion.

**Extended ring-free diagnostic.** A single static pass over the previously measured ARM64 release
binaries identifies concrete work in pending's optional-row representation. Its extended successful
conversion loop has 29 instructions per bond: read the discriminant/payload, stage the 96-byte row
through the stack, check, and write it back. The basic loop has seven instructions and does not copy
its 40-byte payload. These are static loop counts, not measured retired instructions or time fractions.
Allocation reuse does not remove the extended copy. Ordinary extended bond insertion also differs:
pending publishes fields individually, while the sparse reference's final transfer uses paired vector
loads/stores. Closing storage removes optional-row publication/conversion without another diagnostic
variant.

Current S2b already copies optional bonds through drain/flatten/collect into a new vector, so pending
has not added a whole pass relative to that baseline. The regression is a net balance of copying,
allocation, and cursor/control-flow changes. Pending also outlines the organic extended-atom
constructor where control initializes fields directly; sparse shares that call, and fully bracketed
inputs bypass it. It therefore cannot explain either all ring-free regressions or pending/sparse
differences. The diagnostic identifies removable work, not a complete attribution. No live profiling
or further tuning was needed for these observations.

**Measurement protocol.** One coordinated seven-round comparison of saved clean S2a, current S2b,
pending, and closing release binaries, using Rust 1.96.0 and matching timing-harness sources and probe
lockfiles. The three existing binary hashes match the earlier comparison. There are 56 operation/input
group pairs: 46 basic/extended parse pairs include all four variants; six sampled parse-and-raise and
four reaction/CX pairs compare the three corrected variants. Each round rotates task/variant order
and records three calibrated approximately 80 ms batches. Parsing includes result destruction.
Builds, tests, downstream checks, and allocator probes completed before timing; no concurrent build
or profiling experiment ran. S2a remains an older-semantics cost reference, not an equal-output oracle.

**Results (µs).** Values are medians of round medians. The last column is the median paired change
from pending to closing; it need not equal subtraction of the rounded marginal medians.

| Input | S2a reference | Current S2b | Pending opening | Closing | Paired closing − pending |
| --- | ---: | ---: | ---: | ---: | ---: |
| Basic / sampled Rhea atom stereo | 1.120 | 1.679 | 1.129 | 1.097 | -0.029 |
| Basic / original 10 stereo components | 1.772 | 2.771 | 1.692 | 1.663 | -0.048 |
| Basic / original 100 stereo components | 12.945 | 28.815 | 14.291 | 14.240 | -0.051 |
| Basic / 202 atoms, 100 markers | 4.376 | 10.044 | 7.517 | 7.375 | -0.107 |
| Basic / original 1000-atom prefix | 14.194 | 17.738 | 13.076 | 12.528 | -0.570 |
| Extended / sampled Rhea atom stereo | 1.772 | 2.314 | 1.962 | 1.716 | -0.242 |
| Extended / sampled Rhea bond stereo | 1.966 | 2.374 | 2.111 | 1.849 | -0.265 |
| Extended / original 10 stereo components | 3.039 | 4.075 | 3.354 | 2.888 | -0.470 |
| Extended / original 100 stereo components | 26.790 | 42.690 | 32.867 | 28.460 | -4.509 |
| Extended / 64 bracket atoms | 2.983 | 2.961 | 3.215 | 2.893 | -0.311 |
| Extended / 256-atom chain | 8.059 | 8.048 | 8.777 | 6.646 | -1.760 |
| Extended / 202 atoms, no markers | 7.940 | 7.985 | 8.682 | 8.012 | -0.620 |
| Extended / 202 atoms, 100 markers | 8.012 | 13.697 | 12.305 | 11.010 | -1.337 |
| Extended / original 1000-atom prefix | 32.809 | 35.895 | 42.254 | 28.652 | -13.718 |

The coordinated run confirms the large cursor gain against S2b. Pending saves 1.062 µs (38.2%) on
the original basic ten-component case and 14.449 µs (50.2%) on the hundred-component case, faster in
all seven rounds. Against S2a, those are respectively 4.4% faster and 10.7% slower. Closing changes
the hundred-component result by only −0.051 µs (−0.4%), faster in five of seven rounds, with paired
range −0.349 to +0.166 µs. The basic long-prefix case saves a further 0.570 µs (4.4%) with closing;
most other basic differences are much smaller, commonly 1–4%. The basic dense acyclic control still
costs about 7.4–7.5 µs in either cursor design against S2a's 4.4 µs: optional bond storage does not
explain that remaining correction cost.

Extended closing is a substantive counterexample. The sampled Rhea atom/bond-stereo groups save
0.242/0.265 µs (12.4%/12.7%), faster in all seven rounds. The original ten/hundred-component cases
save 0.470/4.509 µs (both 13.8%), again in all seven. The hundred-component paired range is −5.207
to −4.158 µs. These gains are not limited to densely marked input: closing removes the fully bracketed
regression, and substantially improves ordinary chains. Pending's regressions against current S2b
recur in all seven rounds for 64 bracket atoms (+0.237 µs, 8.0%), the 256-atom chain (+0.793 µs,
9.9%), and the unmarked 202-atom control (+0.716 µs, 9.0%). Closing is near control on that last
case, rather than establishing that all cursor overhead has disappeared.

The extended long-prefix case deserves its full range: closing saves a paired median 13.718 µs
(32.5%), but individual paired savings range from 6.968 to 14.268 µs. Pending's marginal median is
42.254 µs here versus 36.241 µs in the earlier architecture run. This is not evidence of a source
change: the pending timing binary hash is unchanged. It is a reason to report this case as variable,
not turn its largest point estimate into a general percentage. Closing is faster in all seven rounds;
no further run was started to chase the variation.

The combined sampled parse-and-raise lanes show much smaller differences. Basic closing changes are
mixed and below 0.05 µs in paired medians; no consistent gain is established there. Extended closing
saves 0.115–0.259 µs, about 1.5–1.7%, with all rounds faster in two groups and six of seven in the
third. Reaction parsing saves only 0.019/0.031 µs (basic/extended); the sampled CX case saves
0.059/0.111 µs. Standalone parser gains should not be presented as whole-pipeline percentages.

**Memory.** Post-parse requested capacities and retained bytes match pending on the measured inputs;
closing does not compact the over-reserved atom, bond, or stereo vectors. Ring-free ordinary/stereo
allocation and reallocation counts also match. The open-ring record grows to hold an optional
stereo-record/position pair: original ten/hundred stereo-component peaks rise by 64 bytes, and the
label-99 control rises by 1,600 bytes. On the basic hundred-component case, retained bytes stay
226,184 and peak changes from 226,344 to 226,408. The sampled CX case drops from nine allocations
to seven when bond remapping disappears, reducing peak by 64 bytes in both targets. Counts are
requested Rust allocation sizes, not RSS or allocator overhead; no tracked bytes remain after drop.

**Disposition.** This completes the bounded experiment set. The cursor result holds, and basic
opening-order performance is competitive with this actual closing-order implementation. Extended
closing has meaningful additional headroom, including multiple-microsecond differences on larger
inputs, and pending's ring-free countercases are reproduced. Thus opening versus closing cannot be
called performance-neutral, nor can all remaining cost be assigned to required stereo work. The
static evidence and direct-row counterexample support the optional-row machinery as a real target;
they do not isolate every saved microsecond or prove that any further opening-order implementation
would necessarily pay the same price.

**Settled selection.** Retain the cursor/pending design and opening order for both targets. Its
source-order semantics support upcoming SMILES roundtripping; closing order is excluded on that
basis despite its measured extended-parser advantage. Keep the selected pending finalizer, including
safe owning bond conversion and outer stereo-vector reuse. Do not bundle the experimental fresh
stereo vector or three-segment ligand rewrite into integration. The measured extended overhead and
capacity limitations are acknowledged rather than treated as universally resolved. Additional
finalization or extended-path improvements may be considered separately without reopening this
selection or blocking integration. The bounded architectural exploration is closed; the following
subitems cover production integration and the remaining correction/coverage work.

### S2b6 integration fixtures and baseline — 2026-09-09

Nine fixed cases now live in the owning parser test files: seven stereo assembly cases, one reaction
assembly case, and one CX assembly case. They assert opening-ordered bond endpoints and complete
frames for both direct targets, plus side-local classes/mapping and CX bond order/donation where
applicable. Expected values are local literals. The cases cover reverse closures, restored label
reuse, connected marked roots, nested marked branches, terminal/repeated/root virtual H, reaction
namespaces, and reverse-closure CX indices. No general permutation generators or parser changes
were introduced.

Before editing tests, the production IO source matched the coordinated-run snapshot. Rust 1.96.0,
the workload manifests, saved S2b/pending timing binaries, and probe lockfiles also matched. The 250
main and three original full output records remain identical between corrected S2b and pending.
The coordinated timing and allocation/peak/retained-memory results above are therefore reused as
the S2b7 integration baseline; no new performance run was needed. The extended bracket64, chain256,
and unmarked 202-atom countercases remain part of that baseline.

| Baseline artifact | SHA-256 |
| --- | --- |
| Current S2b timing binary | `2e34bacae3e8b666aa0c50f80669755d47e6dbc2afbbba2c5909e630dd043019` |
| Selected pending timing binary | `541d97a27bb85d48222f36782fc33c0794a8e510c609a1da1f407eb39f7aeab2` |
| Main workload manifest | `83e8c3c3f9bf42f4c9c98b86ec126ec08837d964f0588ed60a0d64633c8136b6` |
| Original regression manifest | `106941601b7c9bf273dadc01897a321c4460d0138a798283cad77e54377bcdde` |

Verification passed: all nine focused cases; the full IO conformance/property gate with seed 224
(3,478 library, 6 layout, 2,253 MOL, 407 SDF, 10,032 SMILES conformance, and 6 property tests);
library/test Clippy with both features and warnings denied; formatting and diff checks. The final
focused rerun also passes after matching the existing byte-input/map-literal conventions. Runtime
source is unchanged. S2b6 is complete; S2b7 remains unstarted.

### Maintained conformance and benchmark inputs — 2026-09-09

The approved follow-up adds 191 OpenSMILES conformance fixtures and their snapshots from the fixed
ZINC/Rhea sample. One of the 192 distinct inputs already existed and is reused. Fixture headers
preserve source paths, rows, and identifiers; no existing snapshots were changed. All new snapshot
categories and atom/bond counts were checked against the frozen validation records. The nine S2b6
unit cases remain the exact frame/index tests; this addition does not implement S3's permutation laws.

The maintained SMILES benchmark now includes 18 synthetic controls for stereo scaling,
long-prefix behavior, unmarked chains/brackets, stereo density, serial/nested rings, ring-label size,
and branch depth. Both targets parse directly. The real-input cohorts are defined as local literals
in the benchmark file. They retain the original timing split: 64 ZINC, 63 Rhea atom-stereo, and 64 Rhea bond-stereo inputs, plus the single known wildcard
boundary-difference input as a separate diagnostic cohort. The earlier shorthand of three 64-input
groups described the input pools, not the timed success cohorts. The discrepancy itself remains S2e.

Real-input parse-and-raise cohorts and the reaction/CX controls are also retained. Setup is outside
timing and destruction inside; corpus throughput counts molecules while iteration time covers the
whole cohort. The extended combined lane explicitly includes conversion to basic. Allocation/peak/
retained-memory instrumentation remains separate from Criterion. Input and measurement conventions
are recorded here; no separate benchmark directory or README is introduced.

All 191 new conformance cases pass. All 56 new benchmark lanes pass Criterion's test mode, and the
benchmark passes Clippy with warnings denied. All 18 synthetic controls and all corpus selections
match the frozen exploration inputs exactly. Formatting and diff checks pass. No new timing study
was run, no runtime parser source changed, and S2b7 remains next.

## Implementation plan

S0, S1, S2a, S2b, and S2b1–S2b5 are complete. The subsequent bounded experiments and architecture
selection are also complete. S2b6 is complete; S2b7–S2b8 integrate and verify the selected
cursor/pending implementation. S2c–S2e and S3 retain the remaining flag migration, known corrections,
and limited permutation scope. These
remaining subitems have not started. S2b7 subsumes the former S4 cleanup; earlier S4 references
describe the original sequencing. No mutating git operation or commit is implied.

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
  handling, CX completion ranks, and source spans. Carry focused record-order/root tests and benchmark
  ordinary/ring-free controls at this change. **[dep: S1b]**
- **S2b — Complete SMILES finalization (completed 2026-09-09)** — breaking (red→green within the subitem). Finalize full
  frames and associated configurations for supported descriptors in basic, extended, and reaction-side
  paths; discard private records. Migrate parser construction utilities and exact expected values
  together. Add independently specified ligand-frame cases, the equivalent/mirror regression, and
  all specified failure/unsupported-input cases. Migrate downstream exact output expectations that
  necessarily change at this switch, with a recorded explanation per class. **[dep: S2a]**
- **S2b1 — Verify allocation profiling (completed 2026-09-09)** — additive (green), experimental harness in `scratch/`.
  Run the installed `cargo heaptrack` against a parser workload and verify allocation capture and
  Rust source attribution. Investigate Instruments attachment only if it supplies needed additional
  evidence. Record usable tools and limitations. **[dep: S2b]**
- **S2b2 — Establish profiling workloads (completed 2026-09-09)** — additive (green), experimental harness in `scratch/`.
  Combine real SMILES inputs with controlled cases varying atom count, branching, rings, bracket
  atoms, and stereo density independently where possible. Separate basic/extended parsing and
  parse-only/parse-plus-raise workloads. Pin inputs and build configuration, validate outputs outside
  measurement, and define whether destruction is included. **[dep: S2b1]**
- **S2b3 — Collect the attributed baseline (completed 2026-09-09)** — additive (green), experimental harness in `scratch/`.
  Measure clean elapsed time separately from CPU sampling, allocation instrumentation, and targeted
  CPU-counter recording. Report absolute microseconds and relative differences, allocation counts,
  cumulative bytes, peak live memory, and scaling across the selected workloads. Keep baseline and
  candidate comparisons repeatable with a small runner and measurement protocol reusable by other
  modules. **[dep: S2b2]**
- **S2b4 — Explain parser costs (completed 2026-09-09)** — additive (green), analysis of `umol-io::smiles::parser` and its
  builders. Attribute costs to scanning, table construction, ring handling, stereo finalization, and
  output construction. Use Callgrind where sampling cannot distinguish competing explanations;
  separate necessary representation costs from avoidable work. **[dep: S2b3]**
- **S2b5 — Review optimization opportunities (completed 2026-09-09)** — additive (green), analysis in this document.
  Rank candidates by measured evidence, affected workloads, plausible savings in microseconds,
  implementation complexity, and correctness risks. The deliverable is the baseline and ranked
  analysis, with no production parser changes. Select bounded optimization experiments only after
  review; implementation is not authorized by these exploration subitems. **[dep: S2b4]**
- **S2b6 — Integration fixtures and baseline (completed 2026-09-09)** — additive (green), `umol-io::smiles::parser` tests
  and parser benchmarks. Bring the independently specified cursor cases into the owning test files:
  reverse closures, label reuse after branch restoration, marked roots across a dot, nested marked
  branches, terminal/root virtual H, reaction-side mapping, and CX bond annotations. Define expected
  rows and frames locally following each file's conventions. Keep these fixed regression cases
  distinct from S3's permutation laws. Pin current S2b and the selected pending prototype's exact
  output and timing/allocation baselines before the rewire. Reuse the completed coordinated run after
  checking source/input/build identity; remeasure only if those inputs have drifted. Keep the extended
  ring-free countercases visible. No exploratory data or code goes in `materials/`.
  Run the focused tests and IO conformance/property gate. **[dep: S2b5, settled pending selection]**
- **S2b7 — Cursor/pending parser integration** — breaking private rewire (red→green within the
  subitem), `umol-io::smiles::parser` and `parser::builder`. Integrate the statically specialized
  basic/extended assembler, atom cursor, branch restoration, and direct final atom rows. Keep organic
  atom emission separate from marked bracket handling; retain private incidence records only for
  supported marked sites. Reserve optional bond rows at opening and complete them at closure. Remove
  the superseded full-bond stereo scan, site lookup, separate root collection, and closing-record
  replay. Use the selected consuming finalizer with safe owning bond conversion and stereo-vector
  reuse; preserve current capacity hints. Migrate molecule/reaction wrappers and private builder tests
  together, including direct mapping accumulation and private CX records. Preserve exact rows, frames,
  source spans, empty/dot behavior, diagnostics, and both direct target representations. Keep
  ChiralityFrame until S2c. No public symbol changes, closing-order path, storage selector, new
  compaction policy, or experimental finalization rewrite. Include S2b6's regressions and the IO
  conformance/property gate; fit the production files' conventions rather than copying harness
  scaffolding. **[dep: S2b6]**
- **S2b8 — Integration verification and performance handoff** — additive (green), parser/raise
  validation, benchmarks, and this record. Compare the integrated parser against the selected pending
  prototype and corrected S2b: full basic/extended tables and frames, errors, reaction mapping, CX
  annotations, and raised outputs. Require exact identity for this opening-order rewire, with focused
  regressions for any discrepancy. Run the affected IO and graph suites. After builds finish, repeat
  the established bounded timing and allocation/peak/retained-memory comparison, keeping S2a solely
  as an older-semantics cost reference. Verify integration preserves the selected prototype's gains
  and explicitly retain its known extended/capacity limitations. Investigate an integration mismatch
  within this rewire; do not start another architecture or tuning search. Record the handoff before
  removing the flag. **[dep: S2b7]**
- **S2c — Remove ChiralityFrame** — breaking (red→green within the subitem). Remove the enum,
  molecule fields, conversions, parser assignments, and old source-descriptor interpretation branch;
  migrate all imports, literals, and tests. Remove helpers used solely for the retired SMILES frame
  reconstruction. Retain helpers still used by wedges or directional stereo, raw MOL parity, and
  ConfigurationScope. Verify that raw parity with no operative frame remains unread and that no
  replacement source-format dispatch was introduced. **[dep: S2b8]**
- **S2d — Ingestion and known fixture corrections** — additive (green) with required expectation
  migrations. In `umol-graph::ingest` and its resolution fixtures, exercise direct-entity resolution
  and the known equivalent/mirror cases. Regenerate the three identified incorrect isomer fixtures
  only after independently checking their intended identities. Review `verify_stereo` output changes
  as evidence, not authorization for wholesale fixture regeneration. Include existing stereo failure
  and mismatch policy controls; report newly exposed differences without widening resolver scope.
  **[dep: S2b, S2c]**
- **S2e — Wildcard constraint discrepancy** — green at completion. In `umol-io`, trace where the
  extended path introduces NotAromatic for the isotope-labelled wildcard case recorded in S2b2.
  Establish the intended behavior from the boundary contracts before making a correction; return
  unresolved semantics to review. Correct the discrepancy and add focused basic/extended conversion
  and raising regression coverage, recording any output movement separately from the stereo fix.
  **[dep: S2d]**

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
  frame. **[dep: S2d, S2e]**
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

### Optional optimization scope

The former S4 consuming-finalization cleanup is subsumed by S2b7; it is no longer a separate later
rewire. Further finalization and extended-parser improvements, including allocation-capacity policy,
are optional follow-up topics with no selected implementation design. They are not prerequisites for
this plan, and no new optimization subitems or exploration rounds are scheduled. Preserve opening
order and the public full-frame representation in any separately reviewed follow-up.

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

Remaining critical path: **S2b7 cursor/pending integration →
S2b8 verification → S2c flag removal → S2d ingestion/fixtures → S2e wildcard discrepancy →
S3a limited atom/bond permutation laws → S3b downstream/Python coverage → S3c impact review**.
All stages end green. The wildcard discrepancy remains later required work, not deferred work.
Optional optimization is outside this critical path; the core deliverable does not wait for it.
