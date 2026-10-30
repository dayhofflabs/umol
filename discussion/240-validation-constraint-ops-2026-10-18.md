# 240 — Validation and constraint operations

Status: In Progress
Date: 2026-10-18
Relates: [166](166-molecule-ops-2026-07-27.md),
[118](118-validator-architecture-2026-06-20.md),
[195](195-molecule-constraint-matching-2026-08-12.md),
[236](236-operation-ownership-copy-costs-2026-10-09.md),
[237](237-release-0-10-scope-2026-10-15.md),
[242](242-stereo-transformers-2026-10-29.md)

## Scope

S13, S14, and S16 have moved here from doc 166 and are outside the 0.10 release.
They cover invariant validation and stereo ligand construction, constraint
relocation, and explicit constraint discharge. Their subitem identifiers are
retained so existing references remain usable. S13c's stereo-perception ligand
construction is implemented; the remaining validation and constraint work is
outstanding.
S15a, the counts-resolver valence-capacity correction, remains in doc 166 and
in scope for 0.10.

The relocation and discharge contracts below carry forward the agreed design.
S13a records the rejected grouping and the replacement proposal separately.
This work adds no eager validation, fallible attribute mutation, or automatic
constraint relocation/discharge. New molecule-level constraint validation,
including unpaired-electron coupling, remains in doc 195.

## Count and bond-order validation

The former EntityFieldInvariantsValidator proposal was rejected because its
scope and name grouped unrelated checks. Integer-width checks are excluded:
this work will not enforce i8/u8 limits or replace them with charge-range checks.
Nonnegativity is an intended invariant; the existing electron-conservation
check does not currently establish it. These predicates must remain distinct.

The replacement proposal is to add nonnegativity checks for implicit-hydrogen
and lone-pair counts to ValenceInvariantsValidator, and check individual
localized and dative bond orders separately. A balanced conservation equation
can conceal compensating negative counts; a sum of incident orders can likewise
conceal an invalid individual order. A negative literal should not be hidden
by another undetermined field or an earlier undetermined atom.

The proposed separate BondOrderInvariantsValidator and its diagnostic API have
not been settled. The current DSL specification restricts localized orders to
1–4; that domain needs explicit consideration separately from nonnegativity.
There is no corresponding upper bound of four for dative orders. Extending
ValenceInvariantsValidator would also affect the composite Validator that uses
it; the former promise of unchanged composite behavior cannot be carried over
to that proposal without review.

## Stereo-ligand validator API

Expose the following stateless validator and outcome types from
umol_graph::ops::validate. It has no configuration or constructor method and
is not added to the composite Validator or Python by S13b. The uninhabited
error type retains the common validator result contract.

```rust
pub struct StereoLigandInvariantsValidator;

pub enum StereoLigandInvariantsContradiction {
    MissingActualLigands {
        entity: Entity,
        bearer: AtomId,
        ligands: Vec<AtomId>,
    },
    ImplicitHydrogenCountMismatch {
        entity: Entity,
        bearer: AtomId,
        implicit_hydrogens: i64,
        virtual_hydrogens: usize,
    },
    LonePairCountMismatch {
        entity: Entity,
        bearer: AtomId,
        lone_pairs: i64,
        virtual_lone_pairs: usize,
    },
}

pub enum StereoLigandInvariantsError {}

impl StereoLigandInvariantsValidator {
    pub fn validate(
        &self,
        molecule: &Molecule,
    ) -> Result<
        Solution<(), StereoLigandInvariantsContradiction>,
        StereoLigandInvariantsError,
    >;
}
```

## Stereo-ligand invariant validation

Add a standalone `StereoLigandInvariantsValidator` to check the
model-independent requirement that a stereo entity lists every ligand of its
site. For each stereo-atom frame:

- its actual-atom ligands are exactly the site's localized-bond neighbors;
- its number of virtual implicit-hydrogen ligands equals the site's literal
  `implicit_hydrogens` count;
- its number of virtual lone-pair ligands equals the site's literal
  `lone_pairs` count.

Apply the same rules independently to each endpoint block of a stereo bond,
excluding the opposite site endpoint from that block's actual-atom ligands.
An explicit hydrogen is an ordinary actual-atom ligand: it is present by atom
id in the frame, is not included in the implicit-hydrogen count, and is not a
virtual ligand.

A non-literal inherent count makes the corresponding equality
underdetermined. A literal count mismatch or a missing or extraneous
actual-atom ligand produces a typed
`StereoLigandInvariantsContradiction` naming the stereo entity, bearer atom,
and mismatch. Representation integrity already limits each virtual-ligand kind
to zero or one occurrence within a stereo-atom frame or stereo-bond endpoint.
Consequently, an inherent implicit-hydrogen or lone-pair count greater than one
contradicts an existing stereo entity because that entity cannot represent all
of the atom's ligands. This validator imposes no corresponding restriction on
an atom that is not the bearer of a stereo entity.

The check is tier 2, not `Molecule` integrity. A ligand list that disagrees with
the site's neighbors or counts remains structurally interpretable and does not
make ordinary graph-IR operations unsafe; it is a model-independent semantic
contradiction.
The existing integrity rule allowing at most one equal virtual ligand within a
stereo-atom frame or stereo-bond endpoint remains unchanged. Producers of
stereo entities, including stereo perception, must construct ligand lists that
agree with those neighbors and counts. They collect every actual ligand and
add a virtual ligand for each implicit-hydrogen or lone-pair count of exactly
one. A count greater than one precludes producing the stereo entity; it is not
truncated to one virtual ligand.

`HydrogenUnfolder` enforces the implicit-hydrogen requirement directly as an
operation precondition and returns `UnfoldHydrogensError` on disagreement. It
does not require a prior validator call and does not inspect lone-pair support,
which its transformation does not change.

S13b implements this standalone validator; hydrogen transformations do not
depend on a prior validation call.

## Constraint relocation

Decision revised 2026-10-18: move lifting and inlining out of the inherent
Molecule API in graph IR and into separate Transformer implementations under
umol-graph's ops::transform. The engine names are ConstraintLifter and
ConstraintInliner. Transform supplies lift_constraints / lift_constraints_into
and inline_constraints / inline_constraints_into; Python follows the same
receiver-method design as the other transformations. Neither operation runs
automatically during resolution, construction, or serialization.

Lifting must relocate every meaningful inline assertion or return an error.
Leaving an unliftable assertion inline is not successful lifting. Existing
molecule-level entries remain, and successful lifting leaves every inline store
empty. The established exception is kindless stereo: its inline constraints
are vacuous and are discarded without interpreting or normalizing their
payloads. A known kind with an undetermined coset retains meaningful constraints
and follows the ordinary lifting path.

Entity attributes permit values that molecule-level stereo wrappers cannot
represent: a known kind can disagree with the ligand-frame arity, or an inline
constraint can contain an incompatible permutation or ligand position. Lifting
reports the existing representation-integrity failure when such a payload
cannot be relocated. It does not repair the frame, reinterpret the payload,
or require a determined coset. An empty inline store creates no wrapper and
requires no additional stereo check. Attribute mutation remains infallible.

Inlining must move every bare top-level entity assertion into its respective
store or return an error. Partial inlining is not a successful result: an
assertion that cannot be moved must not be left at molecule level while the
others are inlined. Entries with the same key combine by meet, including
entries moved earlier in the same operation; bottom rejects the whole operation.
Logical subtrees, relational constraints, and molecule predicates remain at
molecule level in their stored order; they are outside the relocation's scope.
Inlining is not a traversal that extracts leaves from logical expressions.

A molecule-level stereo assertion carries a kind; its inline counterpart gets
that meaning from the entity's configuration. The current implementation drops
the carried kind. For a meaningful assertion, the replacement must require
agreement before relocation. A different or undetermined entity kind cannot
absorb that assertion without changing its meaning and must cause an error.
Do not change the entity's kind or silently discard the assertion to make
inlining succeed.

Both engines preserve topology, entity ids, inherent fields, configurations,
and ligand frames. Consuming failure drops the input; in-place failure leaves
the input unchanged. Preparation produces a checked candidate before replacing
the caller's borrowed molecule. Discharge is the separate operation for removing
proven redundancy.

Preparation decision, 2026-10-18: use the existing public graph-IR editor
publication boundary from umol-graph. Perform the complete relocation in a
MoleculeEditor and call finish once to obtain the checked candidate. Its full
integrity check can reject an unliftable payload before the candidate is
published to the caller. Do not call probe followed by finish and repeat the
check, or replay the relocation after successfully preparing the result.

The consuming path moves the input into the editor. The _into path prepares a
working copy and replaces the caller's molecule only after finish succeeds.
Cloning shares backing storage; modifying shared stores can copy their
contents. Both paths pay for a full integrity check, including unchanged
topology and overlays. These costs belong in the bounded comparison for S14.

No broader public tier-1 checking API or duplicate integrity implementation is
needed for this work. MoleculeEditor::finish returns MoleculeIntegrityError;
lifting translates the corresponding stereo-admission failures into
LiftConstraintsError. S14 replaces the former local kindless-panic correction
with this migration.

### Relocation errors — agreed 2026-10-18

Expose these operation errors from umol-graph's ops::transform:

```rust
pub enum LiftConstraintsError {
    StereoKindSiteMismatch {
        entity: Entity,
        kind: StereoKind,
    },
    StereoLigandArity {
        entity: Entity,
        kind: StereoKind,
        expected: usize,
        actual: usize,
    },
    StereoPermutationDegree {
        entity: Entity,
        expected: usize,
        actual: usize,
    },
    StereoLigandPositionOutOfRange {
        entity: Entity,
        position: usize,
        degree: usize,
    },
}

pub enum InlineConstraintsError {
    Contradiction {
        entity: Entity,
    },
    StereoKindMismatch {
        entity: Entity,
        asserted: StereoKind,
        actual: Option<StereoKind>,
    },
}
```

LiftConstraintsError retains the names and diagnostic fields of the existing
stereo-admission failures. InlineConstraintsError::Contradiction reports a
same-key meet reaching bottom. StereoKindMismatch reports a meaningful
molecule-level assertion whose carried kind differs from the target entity's
kind; actual is None when that entity has an undetermined kind.

The consuming and in-place method pairs use the same operation error. Python
maps every LiftConstraintsError to InvalidStructureError,
InlineConstraintsError::Contradiction to ContradictionError, and
InlineConstraintsError::StereoKindMismatch to ModelConversionError, retaining
the Rust diagnostic message.

## Explicit constraint discharge

Add an explicitly invoked Transformer under `umol-graph/src/ops/transform` that
removes constraints entailed by independently derived molecular information.
This retains the useful redundancy removal of the automatic closing resolver
pass removed in [236](236-operation-ownership-copy-costs-2026-10-09.md), while
keeping discharge an operation the caller requests.

Start in umol-graph, reusing the model-independent constraint-evaluation
machinery. Move a reusable implementation into graph IR only if a later need
justifies it; no graph-IR operation is part of the initial design.

The intended scope includes entity-local assertions and molecule-level logical
constraints. Remove an assertion only when independently derived information
entails it: compatibility is insufficient, and evidence incorporating the
assertion cannot establish its redundancy. Preserve undecidable assertions and
all topology, entity identities, inherent fields, and relation data.

Discharge can act on a partially specified molecule and requires no chemistry
model or prior resolution. The broader definition of Transformer is outside
this work item.

### Agreed discharge contract — 2026-10-18

The parameterless engine is ConstraintDischarger, separate from ConstraintLifter
and ConstraintInliner. It implements Transformer through transform,
transform_into, and transform_iter. Transform exposes discharge_constraints and
discharge_constraints_into on Molecule; both take no operation arguments.
Python exposes the same named methods directly on Molecule, with the ownership
contract of the other transformations. Its contract is:

- Use the open derived reading. An existing overlay can establish a positive
  assertion, but absence alone does not justify deleting a negative assertion.
  Intrinsic fields and actual topology remain usable evidence.
- Exclude all stored assertions from the evidence, both inline and at molecule
  level. This prevents self-justification, mutual justification, and dependence
  on the order in which assertions are considered.
- Require entailment rather than compatibility. A derived set {2, 3} proves
  an assertion allowing {1, 2, 3}, but not one allowing {3, 4}. Use the existing
  lattice matching contract; derived information need not be literal.
- Remove a complete stored entry only when that evidence proves it true.
  Evaluate logical expressions to decide whether the whole entry is redundant;
  retain the original expression when it is not proved. Do not rewrite its
  internal tree as part of this operation.
- Retain false, undecidable, and unsupported assertions. A false assertion is
  available for an explicit validator to report; discharge does not repair it.
  Inconsistent evidence proves nothing, including under negation.
- Apply the established vacuity rule to kindless inline stereo payloads:
  discard them without inspecting their contents. Meaningful stereo assertions
  remain unless independent information proves them; do not use another stored
  assertion as that proof.

ConstraintDischarger is infallible: Transformer::Error is std::convert::Infallible.
An inability to establish redundancy retains the assertion. The Transform
methods therefore return Result<Self, Infallible> and Result<(), Infallible>;
Python returns the transformed molecule or None for the _into method, with no
operation-specific exception.

Use the fixed ring selection prescribed by the constraint contract. Connected
components use the existing Bfs algorithm internally. There is no discharge
configuration or caller-selected algorithm. Ring and component results are
computed lazily and reused within the invocation.

The current ConstraintInvariantsValidator cannot be reused wholesale: its
stereo comparison reads stored assertions, and relational atom predicates use
the closed derived reading. Reuse individual independent readings and evaluation
code where their contracts agree. The discharge evaluator must preserve the
chosen open reading through nested predicates and reuse ring/component results
within one invocation rather than perceiving them separately for each entry.
This does not require a new general-purpose reasoning framework or a change to
the validator's existing behavior.

The current derived readings for constraints on stereo entities are
undetermined. Non-vacuous assertions in those families therefore remain unless
independent information proves them; this operation does not introduce a new
stereo perception or symmetry calculation.

The required laws are preservation of meaning, topology, entity ids, inherent
fields, and relation data; idempotence; and independence from constraint
traversal order. Exact cases must distinguish proved redundancy from mere
compatibility, false predicates, missing information, contradictory readings,
logical negation, and duplicate assertions. S16 implements and verifies these
contracts.

## Transformation interface and module organization

Extend the Transform trait implemented for Molecule in umol-graph. Its existing
receiver interface and Python ownership contract were completed in doc 166 S12.
Transformer remains the interface for the operation engines. Add no free-standing
molecule functions or Python transformer class hierarchy.

| Engine | Consuming Rust method | In-place method | Error |
| --- | --- | --- | --- |
| ConstraintLifter | lift_constraints | lift_constraints_into | LiftConstraintsError |
| ConstraintInliner | inline_constraints | inline_constraints_into | InlineConstraintsError |
| ConstraintDischarger | discharge_constraints | discharge_constraints_into | std::convert::Infallible |

All methods take no operation arguments. Consuming Rust methods return
Result<Self, E>; _into methods take &mut self and return Result<(), E>.
Consuming rejection drops the input, while in-place rejection preserves it.
Each invocation applies and publishes independently; a later failed operation
does not undo an earlier successful call.

Python exposes the same names on Molecule. Base methods preserve the receiver
and return an independently mutable molecule by cloning the Rust receiver before
consuming application. The clone shares backing storage until mutation requires
copying. The _into methods mutate the receiver and return None, with the
established view/iterator invalidation behavior. The relocation error mappings
are specified above; discharge has no operation-specific exception.

Use ops::transform::constraints for all three engines, with a tests submodule
and one child module per transformer type. Retain type re-exports from
ops::transform. This follows the existing subject-based hydrogens and aromaticity
modules. Create the constraints module with whichever engine is implemented first.

ConstraintLifter and ConstraintInliner prepare checked candidates and can reject
while doing so. Keep their working-copy publication contract; do not replace it
with Molecule::mutate or introduce a second planning pass. ConstraintDischarger
filters stores directly after making removal decisions; do not introduce an
Edits batch just to use Molecule::mutate.

## Implementation plan

The original S13, S14, and S16 identifiers are retained. References qualified
166:S12a and 166:S12c denote the completed Rust receiver interface and Python
ownership/lifecycle work in [166](166-molecule-ops-2026-07-27.md).

### S13 — Invariant validation and stereo ligand construction

S13b's public additions are StereoLigandInvariantsValidator,
StereoLigandInvariantsContradiction, and StereoLigandInvariantsError with the
signatures and payloads above. S13a's replacement API remains a proposal; the
rejected EntityFieldInvariantsValidator inventory is not an implementation target.

- **S13a — Count and bond-order validation** — design review before implementation.
  Settle the replacement described above, including its public diagnostics and
  the composite-validator consequences, before specifying implementation modules.
  Do not restore integer-width checks or the rejected EntityField grouping.
  Preserve the planned coverage of partial inputs, concrete contradictions
  despite other unknowns, traversal order, unaffected fields, and unchanged
  input storage. Exact cases and properties must follow the accepted predicates.

- **S13b — StereoLigandInvariantsValidator** — additive (green).
  Modules: new `umol-graph/src/ops/validate/stereo_ligand.rs` and
  `umol-graph/src/ops/validate.rs`. [dep: none; existing Molecule views]
  Implement the validator and its two outcome types using the neighbor and
  virtual-ligand count checks above. Check stereo atoms and each stereo-bond
  endpoint block, excluding the opposite site endpoint. Count explicit
  hydrogens as actual atom ligands. Compare virtual hydrogen/lone-pair counts
  with the bearer's inherent counts; non-literal counts contribute
  underdetermination without concealing a concrete contradiction elsewhere.
  Use the existing neighborhood and ligand accessors. Rely on molecule
  integrity for references, incidence, and duplicate-ligand admission rather
  than adding another whole-molecule integrity pass.

  Establish exact cases for matching ligand lists, omitted actual neighbors,
  explicit hydrogen ligands, zero/one virtual ligands, and count mismatches,
  including counts greater than one. Cover both stereo entity kinds, both
  endpoint blocks, their allowed stored orders, no stereo entities, and
  non-literal counts. Exercise every contradiction variant with exact payloads
  and unchanged-input assertions. Construct fixtures through public molecule
  APIs without invalid references or bypassing integrity. Add rustdoc and run
  the focused unit tests; this check is independent of S13a.

- **S13c — Stereo ligand construction** — additive (green).
  Module: `umol-graph/src/ops/stereo.rs`, with affected resolver and conformance
  tests. [dep: none; existing stereo perception and the settled ligand contract]
  Make derive_stereo_atom and bond_side_ligands construct the actual and virtual
  ligand list. Read the inherent hydrogen and lone-pair counts;
  literal zero contributes no virtual ligand, one contributes the corresponding
  ligand, and any other or non-literal count prevents a determined construction.
  Collect every actual neighbor, excluding the opposite site for a stereo-bond
  endpoint; reject a count inconsistent with the selected kind's frame size.
  Include both virtual kinds when each occurs once. Preserve the existing
  actual-neighbor order and append hydrogen before lone pair, matching the
  ligand-kind order. Existing None and StereoInconsistency outcomes suffice;
  add no validation pass, public symbol, or new policy.

  Cover omitted extra neighbors, hydrogen/lone-pair counts zero, one, greater
  than one and non-literal, both virtual kinds together, and both stereo-bond
  endpoints. Retain the successful ordinary stereo cases and check the existing
  resolver error/keep/remove policy behavior on unsuccessful construction.
  Run focused perception, resolution, and stereo-conformance tests.

  **Completed — 2026-10-29.** StereoPerception now collects every actual
  neighbor and adds virtual hydrogen and lone pair independently for literal
  counts of one. Both counts must be literal zero or one, including when actual
  neighbors alone fill the frame. The ligand list must contain four
  positions for Th or two per Ct endpoint; the opposite endpoint is excluded.
  Ct endpoint frames sharing an actual ligand are unrepresentable and rejected.
  Existing neighbor order is preserved, with hydrogen before lone pair.

  Exact cases cover explicit H atoms, omitted extra neighbors, zero/one,
  negative, greater-than-one and non-literal counts, both virtual kinds together,
  both Ct endpoints, and preservation of literal and undetermined cosets.
  Resolver tests check Error/Keep/Remove for both assertions and stored entities,
  consuming/in-place agreement, and unchanged input on Error and Keep.
  Stereo-conformance tests cover the new successful frames and failures.
  The stereo-discharge property now supplies explicit zero lone-pair counts
  at its bearers while retaining its discharge, preservation and idempotence
  assertions and unrelated open inherent fields.

  Validation passed: 2,779 graph unit tests, 186 graph property tests, 683
  resolution-conformance tests, and all-target/all-feature Clippy for graph
  and graph IR. Public signatures and outcome types are unchanged; no standalone
  validator was added. S13b and S13d remain outstanding. The S1 benchmark
  measurements are recorded directly in
  [242](242-stereo-transformers-2026-10-29.md#s1--shared-stereo-foundations).

- **S13d — Validator and producer properties and verification** — additive (green).
  Modules: `umol-graph/tests/property.rs`, new
  `umol-graph/tests/property/validate.rs`, and the S13 public rustdoc/exports.
  [dep: S13a, S13b, S13c]
  Generate integrity-valid molecules with independently chosen field values
  and ligand lists that agree or disagree with the neighbors and counts. Check
  agreement with the declared numeric domains and independently constructed
  neighbor and virtual-ligand counts, preservation of the input, and independence
  from fields outside each validator's scope. Check outcome classification under atom/entity
  relabeling and allowed ligand-frame permutations. Check transported
  diagnostic ids on cases with one violation; do not require an unchanged
  first-reported contradiction when several violations are reordered.
  Check that frames produced by StereoPerception satisfy
  StereoLigandInvariantsValidator, with independently generated neighbor and
  virtual-ligand counts.

  Run the affected validator unit-test groups and the new proptest suite
  explicitly, then review public symbols and rustdoc against S13b's inventory
  and the accepted S13a design. Check the agreed composite-validator behavior;
  preserve existing validation categories and hydrogen operation preconditions.
  Review traversal and allocation for repeated work; use a bounded measurement
  only if a concrete cost question needs evidence for acceptance.

  Run the workspace/Python gate specified in [166, S10d](166-molecule-ops-2026-07-27.md)
  once after this related group, applying python-build before commands that
  include umol-py. Share the final gate with other stages completed together,
  retaining each subitem's focused checks and API review. Record the results
  and mark the completed subitems in this document.

S13a, S13b, and S13c are independent; all precede S13d. S13a needs its design
review before implementation. No mutation, transformer, or Python
implementation is a prerequisite.

### S14 — Constraint relocation transformers

The former S14a, an infallible local correction to Molecule::lift_constraints,
is replaced by the relocation design above. Do not implement that superseded
subitem. The public additions are ConstraintLifter, ConstraintInliner,
LiftConstraintsError, InlineConstraintsError, and the two Transform method
pairs specified above. Both engines are parameterless. Preparation uses the
existing checked editor boundary; no new integrity-checking API is required.

- **S14a — Relocation fixtures and baseline** — additive (green).
  Modules: `umol-graph/benches/transform.rs` and the existing graph-IR molecule
  tests. [dep: none; run before removing the inherent methods]
  Establish exact successful fixtures for lifting and inlining and capture a
  bounded baseline through the current methods. Include ordinary literal and
  undetermined forms, sparse assertions on larger topology, same-key meets,
  and a small non-literal/stereo supplement. Check the expected relocation
  outside timing; do not treat the existing kindless panic or kind-dropping
  behavior as a reference result.

  Distinguish uniquely owned inputs from inputs whose backing storage remains
  shared with a retained clone. Parsing and setup stay outside timing. Record
  the fixture dimensions, source revision, commands, compiler, and timings for
  the S14d comparison. There is no existing in-place method to serve as an
  equivalent _into baseline.

- **S14b — ConstraintLifter and LiftConstraintsError** — additive (green).
  Modules: new `umol-graph/src/ops/transform/constraints.rs` and
  `umol-graph/src/ops/transform.rs`. [dep: S14a]
  Implement the engine, its settled error variants, and all three Transformer
  methods. Move consuming input into the editor; prepare a working copy for
  _into. Drain all eight inline stores, preserve existing molecule entries,
  discard kindless stereo payloads uninspected, and call finish once. Translate
  the specified stereo-admission errors; do not add a catch-all public error.
  A successful prepared candidate is the result, not a plan to replay.

  Clarify Transformer rustdoc's preparation/publication distinction: preparing
  a candidate may mutate owned working storage and reject at finish; the
  caller's borrowed molecule is replaced only after preparation succeeds.
  This changes neither the trait signatures nor the broader transformation
  classification.

  Add exact cases for every entity kind, empty stores, existing top-level
  entries, known stereo kinds with undetermined cosets, and malformed kindless
  payloads. Exercise every error variant with a meaningful payload and verify
  that an empty stereo store requires no wrapper admission. Put a movable
  assertion before and after an unliftable one; both orders must reject with
  unchanged _into input. Check consuming/in-place agreement and iterator
  success, rejection, and source independence. Add rustdoc and run focused tests.

- **S14c — ConstraintInliner and InlineConstraintsError** — additive (green).
  Modules: `umol-graph/src/ops/transform/constraints.rs` and
  `umol-graph/src/ops/transform.rs`. [dep: S14a, S14b preparation contract]
  Implement the engine, its two error variants, and all three Transformer
  methods. Move bare top-level entity assertions into their stores, meeting
  same-key entries and checking meaningful stereo-kind agreement. Preserve
  every non-relocated entry and its stored order. Use the same checked-candidate
  publication contract; do not extract leaves from logical trees.

  Cover all eight entity kinds, existing inline entries, repeated top-level
  keys, compatible meets and bottom, matching known kinds with undetermined
  cosets, different kinds, and kindless targets. Check both positions of a
  rejecting assertion among movable entries, exact errors, unchanged _into
  input, and preservation of logical, relational, and molecule predicates.
  Check iterator behavior and source independence. Add rustdoc and run the
  focused tests alongside S14b.

- **S14d — Relocation cost comparison** — additive (green).
  Modules: `umol-graph/benches/transform.rs` and this record. [dep: S14b, S14c]
  Compare the consuming engines with S14a on equivalent successful inputs,
  before removing the old methods or adding bindings. Invoke the new engines
  explicitly so the old inherent methods cannot mask them. Measure _into
  separately, including unique/shared backing and the working-copy cost.
  Sparse-constraint cases expose the cost of checking unchanged topology.
  If the result leaves acceptance unclear, separate that full finish check
  from relocation work in a bounded diagnostic run; inspect repeated checks,
  cloning, normalization, and meets before further timing.

  Report actual added work and any remaining cost rather than promising an
  expected reduction. Do not add an unchecked publication path to recover it.
  Correct avoidable work and record the acceptance decision before the API
  cutover. Retain these results for the final gate; rerun only when subsequent
  changes alter the measured path or leave a concrete concern unresolved.

- **S14e — Relocation API cutover and property migration** — breaking (red→green).
  Modules: `umol-graph-ir/src/ir/molecule.rs`, its molecule tests,
  `umol-graph-ir/tests/property/edit.rs`, `umol-graph/src/ops/transform.rs`,
  `umol-graph/src/ops/resolve.rs`, `umol-graph/tests/property.rs`, new
  `umol-graph/tests/property/constraint_relocation.rs`, and the S14 benchmark.
  [dep: S14d, 166:S12a]
  Remove the two inherent graph-IR methods and migrate their callers while
  adding the two Transform method pairs. Inherent methods take precedence
  over extension-trait methods; do not leave the old implementation masking
  the new receiver interface. Migrate the benchmark and resolver rustdoc too.
  Move operation tests and their necessary generators to umol-graph, retaining
  graph-IR storage and serialization tests. Add no reverse crate dependency.

  Preserve the existing drain, deposit, and paired-idempotence laws, extending
  their generated coverage to stereo entities. Separate liftable/compatible
  domains from deliberate rejection domains; do not silently discard errors
  in a success property. Check individual idempotence, consuming/_into result
  and error agreement, unchanged input on rejection, and preservation of
  topology, ids, inherent fields, frames, and untouched constraint entries.
  Check successful relocation roundtrips by normalized per-entity assertions
  and unchanged non-entity entries, not raw ordering or restoration of vacuous
  payloads. Derive expected meets from the generated assertions, not by calling
  the opposite relocation engine as the sole oracle.

  Run both engines' unit tests, graph and graph-IR property targets explicitly,
  and graph rustdoc. Keep the breaking removal and all consumer migrations in
  this subitem so the tree is green at its end.

- **S14f — Python constraint relocation methods** — additive (green).
  Modules: `umol-py/src/molecule.rs`, `umol-py/src/error.rs`, and the existing
  molecule and molecule-view Python tests. [dep: S14e, 166:S12c]
  Bind both method pairs through Transform, with the recorded exception
  classes and Rust diagnostic messages. Reuse doc 166 S12's copying and accessor
  lifecycle. Cover relocation, same-key meets, each exception mapping,
  source/result independence, atomic rejection, and view/iterator behavior
  on success, rejection, and no-op execution. Rebuild with python-build before
  focused Python tests; retain Python 3.9-compatible syntax.

- **S14g — Relocation documentation and verification** — additive (green).
  Modules: operation rustdoc, the Python API guide, and this record.
  [dep: S14e, S14f]
  Reconcile the public-symbol diff, operation errors, documentation, migrated
  properties, all callers, and S14d's accepted cost results. Run the gate
  specified in [166, S10d](166-molecule-ops-2026-07-27.md) for the completed
  related group, sharing it with S16e if implemented together.
  Record verification results and mark the completed subitems.

The engine and cost path is S14a → S14b → S14c → S14d. 166:S12a joins it at
S14e, and 166:S12c at S14f. The engines do not require hydrogen, validator, or
counts changes; only their receiver-method integration depends on the completed
doc 166 S12. All seven subitems are required.

### S16 — Explicit constraint discharge

The public additions are the parameterless ConstraintDischarger and Transform's
discharge_constraints / discharge_constraints_into pair, with the corresponding
Python methods. Its error is Infallible. There is no discharge configuration,
public evaluation context, new chemistry model, or graph-IR operation.

- **S16a — ConstraintDischarger and independent evaluation** — additive (green).
  Modules: `umol-graph/src/ops/transform/constraints.rs`,
  `umol-graph/src/ops/transform.rs`, `umol-graph/benches/transform.rs`, and
  existing constraint-evaluation modules only where compatible logic is
  actually shared. [dep: none; existing views and constraint contracts]
  Start with exact behavior fixtures, then implement the private evaluation
  context and its engine consumer together. Exclude all stored assertions
  from evidence, use open derived readings throughout nested predicates, and
  distinguish proved true, proved false, and lack of proof. Bottom evidence
  cannot become proved false and then true under Not. Use entailment for
  non-literal readings; overlap is insufficient. Preserve the established
  logical and relational predicate semantics, including empty cases.

  Reuse compatible reading/evaluation code without changing existing validator
  results. Do not delegate wholesale to ConstraintInvariantsValidator or its
  closed relational atom evaluation. Use the prescribed ring selection and
  Bfs internally; one invocation owns lazy ring/component results reused by
  top-level and nested predicates.

  Implement all three Transformer methods. Prepare removal decisions before
  publication, then filter each affected store once using existing storage
  operations. Avoid repeated whole-list removal searches or normalization of
  retained entries. Discard kindless inline stereo payloads uninspected; keep
  meaningful stereo assertions unless independent evidence proves them.
  Preserve every retained payload and its order. No fallible candidate repair
  or working-copy rollback is needed for this infallible removal operation.

  Exact cases cover inline and molecule-level entries across all entity kinds;
  literal and set entailment versus mere overlap; true, false, unknown, and
  bottom readings; existing versus absent overlays; vacuous and meaningful
  stereo; whole logical trees including negation; relational atom predicates;
  ring and connectivity predicates; and duplicate assertions across stores.
  Verify unsupported predicates remain, nested predicates retain the open
  reading, and unproved trees are unchanged. Include empty/no-op cases and
  mixed retained/removed entries. Test consuming/_into/iterator results and
  source independence. Run focused engine tests and any validator tests touched
  by shared implementation changes. Establish the benchmark correctness
  fixtures for S16e alongside these cases. Add the engine's contract rustdoc.

- **S16b — Discharge properties** — additive (green).
  Modules: `umol-graph/tests/property.rs` and new
  `umol-graph/tests/property/constraint_discharge.rs`. [dep: S16a]
  Generate partial molecules with assertions independent of their intrinsic
  and topological data. Check idempotence, consuming/_into agreement, source
  independence, representation integrity, and exact preservation of topology,
  ids, inherent fields, configurations, and frames. Surviving entries must be
  an unchanged subsequence of each original store. Reordering or duplicating
  assertions must not change which assertions independent evidence proves;
  compare corresponding survivors, preserving retained multiplicity.

  Check preservation of meaning with a bounded independent oracle: generate
  finite numeric/Boolean domains and small topology with directly enumerated
  facts, then compare admissible assignments before and after discharge.
  Include partial overlap, contradictions, missing evidence, logical negation,
  and assertions stored in different scopes. Do not use the production
  discharge evaluator or the closed validator as the semantic oracle. Keep
  the oracle's finite evidence scope explicit; exact S16a cases cover the
  other predicate families. All aggregate generation uses public construction
  and mutation without widening production visibility. Run the property target
  explicitly and document the laws with the public engine.

- **S16c — Discharge methods on Transform** — additive (green).
  Module: `umol-graph/src/ops/transform.rs`. [dep: S16a, 166:S12a]
  Add the two argument-free receiver methods with Result<Self, Infallible>
  and Result<(), Infallible>, delegating to ConstraintDischarger. Document
  explicit invocation and partial-molecule support. Check public method results
  and independent publication in an ordinary sequence; resolution and other
  transformations must not invoke discharge automatically. Run focused graph
  tests and rustdoc. No composite transformer or additional method family is
  part of this integration.

- **S16d — Python discharge methods** — additive (green).
  Modules: `umol-py/src/molecule.rs` and the existing molecule and molecule-view
  Python tests. [dep: S16c, 166:S12c]
  Bind discharge_constraints and discharge_constraints_into through Transform.
  Add no config, Python engine class, or operation-specific exception. Cover
  proved redundancy, retained unknown/false assertions, whole logical entries,
  source/result independence, and the established _into view/iterator
  invalidation on changed and no-op results. Rebuild the extension with
  python-build and run the focused Python tests with Python 3.9-compatible
  test syntax.

- **S16e — Discharge cost and combined verification** — additive (green).
  Modules: `umol-graph/benches/transform.rs`, operation rustdoc, the Python API
  guide, and this record. [dep: S16b, S16c, S16d]
  Use S16a's correctness fixtures to measure a bounded set of common literal
  and undetermined assertions, mixed retained/removed lists, and repeated
  ring/component predicates, including nested forms. Vary assertion count
  separately from topology size to expose
  repeated list scans or repeated perception. Compare consuming/_into behavior
  with unique/shared backing. Parsing and setup stay outside timing; a small
  set/expression supplement must not dominate the workload.

  There is no previous explicit discharge implementation to benchmark against.
  Use absolute cost, scaling, and inspection of computation/allocation counts
  to decide whether repeated work needs correction. Inspect that ring/component
  data is computed at most once per invocation and not at all when unused.
  Extend measurements only for an actionable cost question; do not require
  broad resolver/transformer benchmark runs for this separate operation.

  Review the full diff, public symbols, documented laws, property domains,
  Python ownership, and unchanged validator behavior. Run the gate specified
  in [166, S10d](166-molecule-ops-2026-07-27.md) once for the completed related
  group, consolidating other completed stages' gates when implemented together.
  Record results and the cost decision; close this document only when all its
  required work passes.

S16a → S16b supplies semantic verification; S16a and 166:S12a → S16c, then 166:S12c
joins at S16d. These paths meet at S16e. S16 does not require relocation,
hydrogen resolution, or the new validators to supply evidence. All five
subitems are required.

### Sequencing, measurements, and test gates

- Resolve S13a's replacement API before its implementation and before finalizing
  its part of S13d. S13b and S13c retain their independent contracts.
- Preserve the S14a baseline before replacing the inherent methods. Add both
  engines and settle the measured cost in S14d first; remove the old methods,
  migrate callers and properties, and install the receiver methods together
  in S14e.
- Candidate preparation and caller publication remain distinct in S14: one
  finish check, no probe plus finish, and no replay of the checked result.
  S16's evaluator and first consumer stay together, without an unused public
  framework or intermediate validator behavior change.
- Required cost work is S14a/S14d for relocation's full check and copying, and
  S16e for discharge traversal and perception reuse. Forwarding methods need no
  independent timing. S13 needs measurements only for a concrete unresolved cost.
- Every implementation subitem carries focused exact tests. S13d, S14e, and
  S16b carry the generated laws. Organize tests by operation and method, splitting
  long modules by subject, without stage-named modules or new external fixture
  collections. Existing conformance examples serve regression coverage;
  ordinary tests cover distinct behavior.
- Run affected unit/Python tests and feature-gated property targets while
  iterating. Run the full gate specified in doc 166 S10d once after each related
  implementation group, not after every subitem. It includes formatting,
  all-feature workspace tests and doctests, property/conformance suites,
  strict lint/rustdoc, MSRV checking, and rebuilt Python tests. A shared gate
  does not replace each subitem's scope/API review or erase failed focused tests.

## Completion criteria

- S13a's replacement is agreed, implemented, and verified against its explicit
  predicates without integer-width checks.
- The standalone stereo-ligand validator checks ligand lists against actual
  neighbors and inherent hydrogen/lone-pair counts without imposing a chemistry
  model. Stereo perception constructs ligand lists satisfying those checks.
- ConstraintLifter and ConstraintInliner replace the inherent graph-IR methods.
  Lifting discards vacuous kindless stereo payloads and either relocates all
  meaningful assertions or errors. Inlining moves every bare top-level entity
  assertion or errors, preserving stereo-kind meaning. Both operations leave
  an in-place input unchanged on rejection.
- Explicit constraint discharge is infallible, removes only independently
  proven redundancies and vacuous payloads, and retains false, undecidable,
  unsupported, or inconsistently derived assertions. It preserves molecular
  structure and meaning on partial as well as concrete inputs, is idempotent,
  and is independent of constraint traversal order.
- Transform and Python expose the agreed method pairs with their recorded
  ownership, error, and publication behavior. Required tests and bounded cost
  comparisons pass before closeout.
