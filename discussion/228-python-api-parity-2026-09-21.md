# 228 — Python API semantic parity review

Status: In Progress
Date: 2026-09-21
Relates: [213](213-editor-overlay-storage-2026-08-27.md),
[192](192-python-api-type-roles-2026-08-09.md),
[181](181-python-boundary-ownership-2026-08-03.md),
[Python API guide](../docs/development/python-api.md),
[data-type guide](../docs/development/data-types.md)

## Scope and urgency

Review the public Python API for operations and semantics introduced in bindings
without an owning Rust API. Python is the primary application interface; demand
for convenience methods must not create a second implementation of domain
operations. Python may omit Rust APIs and adapt their syntax for ergonomics. It
must preserve their semantic, reference, failure, and lifecycle contracts.

The immediate trigger is Python-only collection extension for Edits and Deltas.
S0 removed that wrapper-owned extension surface; any replacement must delegate to a
settled operation on the Rust container. Doc 213 is paused for this urgent
review. It retains incremental construction of one Edits sequence and leaves
independent-batch composition unresolved; resolving that composition is not a
prerequisite for the review or removal. If resumed, its assembly and handle
semantics must move from Deltas lowering into Edits methods, with lowering
restructured to use them. Do not implement rebasing in umol-py or preserve the
current methods as compatibility paths.

This document records an urgent review and process correction, including the
first bounded passes below and a bounded P1–P6 implementation plan. It is not a
completed API audit. The broader ownership/type-role review remains in 192 and
copy-cost review in 181.

## Strict review criteria — settled 2026-09-21

Every Python/Rust difference must be justified, without exceptions. Cover types,
methods and names as well as semantics, lifecycle, ownership, aliasing, cloning,
snapshots, and error behavior. The implementer and reviewer supply the evidence;
the user is not required to disprove an asserted PyO3 or ABI necessity.

- Identify the public Rust counterpart, precise difference, concrete Python
  requirement, considered alternatives, and supporting evidence.
- Demonstrate claimed PyO3/ABI restrictions for the relevant version and code
  path. A general Rust/Python ownership mismatch does not establish that copying
  is necessary.
- Trace copying through the boundary and its consumers, including conversions,
  snapshot creation, nested values, and later copy-on-write. Distinguish sharing
  from owned-data copies; evaluate ownership transfer, borrowing, and live access
  where relevant. A consuming Rust operation made reusable through copying is a
  lifecycle deviation, not merely an implementation detail.
- Existing implementations, tests, normative text, and previous acceptance are
  evidence of history, not justification. No deviation is grandfathered or made
  acceptable by repeating it elsewhere.
- Record each inspected difference as justified, unjustified, or unresolved,
  with evidence. Uninspected surfaces remain explicitly unreviewed. Do not call
  a sampled family or the whole API compliant because selected cases passed.
- Use an API inventory to direct bounded adversarial reviews and spot checks.
  Start with reference semantics, consuming operations, snapshots, clones, and
  publication boundaries; follow a confirmed pattern across its shared binding
  machinery. Summaries should state the broken contract, consequence, and owning
  fix in a few lines, leaving detailed evidence accessible separately.

These criteria are now explicit in AGENTS.md and the Python API guide. S0 repairs
the bounded issues recorded in its outcome below; the remaining audit and fixes
are unfinished. Doc 213 remains paused.

### Consumption and mutation — clarified 2026-09-21

Preserve consuming Rust operations in Python. Possible input reuse is not a
justification for unconditional cloning; callers explicitly copy when they need
reuse. Existing editor/journal consumed states demonstrate that Python consumption
is feasible here. A departure still requires specific evidence. The first-pass
review gave hypothetical reuse too much weight by leaving this principle open.
The alias and invalidation contract is settled below; failure ordering remains a
detail to specify before implementation.

Mutable access must update backing storage; immutable observations must reject
mutation, including through nested containers. An explicitly requested independent
copy may own mutable storage. A property returning a silently disconnected mutable
container violates the rule; describing it as detached supplies no justification.
These criteria are also recorded in the Python API guide.

### Consumption and invalidation errors — settled 2026-09-21

Consumption proceeds with outstanding read-only accessors. It leaves all aliases
of the owner consumed and invalidates its views, nested accessors, and iterators.
Subsequent access raises an explicit lifecycle error; no automatic snapshot or
copy preserves the invalidated access. Explicit independent copies made before
consumption remain usable.

| Exception | Meaning | Base |
| --- | --- | --- |
| ConsumedError | Access to an owned object whose contents have been consumed | RuntimeError |
| InvalidatedViewError | Access through a view or iterator whose backing storage is unavailable because its owner was consumed | RuntimeError |

Use these names across the bindings. Messages identify the object or accessor and
its consumed owner, for example: "Edits has been consumed" and "Edit view is
invalid because its owning Edits has been consumed". S0 implements both exception
exports and uses ConsumedError for existing editor/journal consumption. Owner-backed
entry invalidation remains for later stages. Further lifecycle machinery is not
required for this design decision.

### Reusable access pattern — settled 2026-09-21

Owner-backed access is a general binding pattern, applicable where Rust exposes
borrowed data, not a special mechanism limited to Edit/Delta. An accessor retains
the Python owner and a location, briefly borrowing and checking the owner on each
operation. Nested access follows the same contract. Read-only access rejects
mutation; mutable access updates storage; explicit copies have independent
ownership; consumption invalidates dependent accessors.

Adoption is per type. Owned-only values and accessor-only views remain valid
designs. Combining both roles in one Python class needs a specific justification
and an explicit public mutation contract. An internal owned-or-access representation
is possible but does not imply automatic copy-on-write. The illustrative storage
enum is not an approved universal abstraction. The Python API guide owns this
reusable pattern; this decision does not expand the current implementation scope
beyond P1–P6 or settle the concrete Edit/Delta wrapper layout.

### Shared entry variant classes — settled 2026-09-21

Owned Edit/Delta entries and their read-only accessors share the existing Python
variant classes and field names. Avoid a second parallel hierarchy of accessor
variants; variant inspection and pattern matching use one interface. The owned
versus owner-backed storage distinction stays internal.

A read-only readonly property reports mutation permissions. Batch entry accessors
and their nested access are read-only. Explicit copy() produces an independent
owned entry of the same variant; mutation does not trigger implicit copying.
Readonly describes access permissions, not ownership in general. The bounded
experiment below verifies the class/storage approach for one nested variant;
production replacement of the generated wrappers remains outstanding.

### Entry-access feasibility experiment — 2026-09-21

**Result: feasible for the tested shape.** A standalone embedded-Python executable
used PyO3 0.29.0 (the repository's locked version), abi3-py39, Python 3.13.15 from
umol-py/.venv, and the actual graph-IR Edits, Edit::AddAtoms, and AtomForm types.
It did not change or load the production umol-py extension. Experimental source
is in scratch/228-entry-access; this section preserves the result when scratch
is deleted. No timing or workspace test campaign was run.

The prototype used a PyO3 Edit base class and one AddAtoms subclass, registered
as Edit.AddAtoms. Both construction and batch access produce that same subclass.
The subclass holds either an owned Vec<AtomForm> or a Python Edits owner plus
entry index. The batch holds Option<graph_ir::Edits>. A field access briefly
borrows the owner and matches the stored Rust variant; it never retains a Rust
reference across Python calls. Explicit PyClassInitializer construction works
with the pinned PyO3 version.

The atoms property returns a sequence accessor; indexing that sequence returns
an atom accessor. Both retain the entry wrapper and read its backing Rust data.
The atom accessor permits mutation of an owned entry and rejects it for a batch
entry. It is a minimal element-field probe, not a complete AtomForm binding.
The iterator holds the batch, position, initial length, and exhaustion state.
Its next call creates just the requested AddAtoms wrapper. Consumption takes the
actual Rust batch and consumes its iterator; this is a move probe, not molecular
application or transaction testing.

| Executed case | Observed result |
| --- | --- |
| Construct owned AddAtoms and index a batch entry | Same Python class; both are Edit instances; readonly distinguishes permissions |
| Match Edit.AddAtoms(atoms=values) and Edit.AddAtoms(values) | Both keyword and positional class patterns work |
| Create iterator over two entries | No entry wrapper created and no payload copy |
| Request first entry | Exactly one entry wrapper created; no payload copy |
| Compare exposed diagnostic payload addresses | Accessed entry uses the batch's original atom-vector storage |
| Read nested sequence and atom fields | Original values; no payload copy |
| Assign nested atom element or replace sequence item | TypeError for read-only access; sequence has no append method |
| Call entry.copy() | One explicit vector/payload clone; same variant class, readonly false, different storage address |
| Mutate copied entry's nested atom | Copy changes; original batch entry does not |
| Append 128 entries after creating iterator | Existing accessor remains valid; iterator yields only its original remaining entry |
| Consume batch while aliases/accessors remain | Owner aliases raise ConsumedError; entry, nested sequence, nested atom, copy request, and active iterator raise InvalidatedViewError |
| Read independent copy after owner consumption | Copy remains usable |
| Advance an already exhausted iterator after consumption | Remains exhausted |
| Delete batch variable and collect Python garbage while retaining nested atom | Nested accessor retains owner and stays usable |

Entry-wrapper and explicit-payload-copy counters, plus live payload-address
comparisons, support these assertions. The counters are probe instrumentation,
not a general allocator measurement. Source inspection confirms that the only
payload clone in these access paths is the explicit copy method; Python wrapper
allocation and scalar conversion still occur on demand.

Reproduction command, while the disposable source exists:

```sh
source umol-py/.venv/bin/activate
cargo run --offline --manifest-path scratch/228-entry-access/Cargo.toml
```

The final run compiled without Rust warnings and all assertions passed. This
establishes that payload copying is not required to preserve the selected variant
interface. It requires replacing generated stored-field getters with owner-aware
getters; it is not a change to an enum annotation alone. A complete Delta nesting
chain, equality/repr, constructor compatibility, structural field replacement,
batch copying, and application failure ordering were not implemented or verified
by this bounded experiment. It does not prescribe a generic wrapper framework or
claim that the full migration is complete.

## Confirmed trigger

- Rust Edits documents ordered execution, initial-host Id handles, and per-kind
  creation ordinals New(n). Its public surface deliberately excludes batch
  concatenation. See [Edits](../umol-graph-ir/src/ir/edit.rs).
- Python Edits.extend accepts another Edits or raw entries, snapshots them, and
  pushes them unchanged. It does not distinguish independent creation namespaces
  from entries already written in the destination namespace. See
  [binding](../umol-py/src/edit.rs).
- For example, let A add carbon, and B add nitrogen then remove its New(0).
  Unchanged concatenation executes add carbon, add nitrogen, remove New(0): the
  removal targets carbon. Preserving B's creation reference requires New(1).
  This was reproduced through molecule application in the first pass below.
- The Python extension test checks list equality, copying, and self-extension;
  it never applies the combined batch to establish which atom the handle names.
  See [test_edits_extend](../umol-py/tests/test_edit.py).
- Python Deltas.extend also implements collection extension in the wrapper by
  repeatedly calling Rust push; Rust Deltas has no corresponding bulk operation.
  Deltas in one shared reaction id frame need no creation-handle rebasing, so the
  Edits corruption example does not establish a Deltas reference defect. The
  shared issue is ownership of the operation. See
  [binding](../umol-py/src/delta.rs) and [Rust Deltas](../umol-graph-ir/src/ir/delta.rs).
- Python append delegates directly to Rust push. A Python spelling of that
  existing operation is distinct from independent-batch composition; do not
  classify it as a missing Rust operation merely because the names differ.
- At the start of this review, the Python API guide explicitly endorsed extension
  without rebasing identifiers or New handles. The defective behavior had been
  incorporated into guidance as well as code and tests. The guide was corrected
  on 2026-09-21; S0 removed the wrappers. Checking wrappers against
  the previous guidance alone would have preserved the defect.

These observations establish the contract and verification gap. They do not
establish who approved the change or why; no approval-history claim is made.

## Review groups — settled 2026-09-21

Cover the existing Python API, including returned types, generated methods, and
feature-dependent surfaces. Missing Rust bindings are not findings or requests
for additions. Assign each Python type and method one primary group; cross-link
its dependencies. Additional surfaces discovered during inventory must receive
a group rather than be excluded.

The grouping follows domain operations, even where their implementation spans
crates. Generic witnesses belong to graph-core and molecule-specific witnesses
to graph-ir; review their conversions across that boundary. CoordGen and nauty
are currently exposed through higher-level operations/selectors, so trace those
paths without expanding this into a full native-backend audit. Shared binding
machinery may share evidence, with its members and exceptions identified.

| Domain | Review group | Coverage status |
| --- | --- | --- |
| chem | Chemistry values and Python element namespace | Unreviewed |
| graph-core | Exposed algorithms/selectors and generic correspondence, compaction, remapping | Unreviewed |
| perm | Exposed permutation operations and conversions | Unreviewed |
| coordgen | Python-visible layout selection and backend boundary | Unreviewed |
| nauty | Python-visible algorithm selection and backend boundary | Unreviewed |
| graph-ir | Leaf forms, expressions, literal access, predicates | Unreviewed |
| graph-ir | Atom, including form/update/views and entity constraints | Unreviewed |
| graph-ir | Bond, including form/update/views and entity constraints | Unreviewed |
| graph-ir | Dative bond, including form/update/views and entity constraints | Unreviewed |
| graph-ir | Aromatic system, including form/update/views and entity constraints | Unreviewed |
| graph-ir | Multicenter bond, including form/update/views and entity constraints | Unreviewed |
| graph-ir | Noncovalent bond, including form/update/views and entity constraints | Unreviewed |
| graph-ir | Stereo atom, including form/update/views and entity constraints | Unreviewed |
| graph-ir | Stereo bond, including form/update/views and entity constraints | Unreviewed |
| graph-ir | Top-level constraints | Unreviewed |
| graph-ir | Molecule construction/access, structural queries, combine/split/extract | Unreviewed |
| graph-ir | Molecule editing and transactions | First pass and challenge complete; ownership/error questions open |
| graph-ir | Edit, Edits, handles, and ConstraintEdit | First pass and challenge complete; findings and open questions below |
| graph-ir | Lattice operations and normalization | Unreviewed |
| graph-ir | Reframe and frame transport | Unreviewed |
| graph-ir / graph | Canonicalization and its configuration | Unreviewed |
| graph-ir | Delta, Deltas, field changes, and lowering boundary | First pass and challenge complete; findings and open questions below |
| graph-ir | Molecule correspondence | Unreviewed |
| graph-ir | Molecule compaction | Unreviewed |
| graph-ir | Molecule remapping | Unreviewed |
| graph-ir | Substructure matching and its configuration/results | Unreviewed |
| graph-ir | Reaction and ReactionSpan construction/access/conversion | Unreviewed |
| graph-ir | Reaction application/composition and returned iterators | Unreviewed |
| graph-ir | DSL parse/render, defaults, and molecule/reaction metadata | Unreviewed |
| graph | Models, registries, tables, policies, and configuration | Unreviewed |
| graph | Resolution and its outcomes/reports | Unreviewed |
| graph | Projection | Unreviewed |
| graph | Transformations | Unreviewed |
| graph | Fingerprints and their values/configuration | Unreviewed |
| io | Parse/raise/ingest/interpret | Unreviewed |
| io | Format/export/convey | Unreviewed |
| io | Depiction and its configuration/results | Unreviewed |
| Python boundary | Exceptions, argument adapters, generated methods, exports, feature availability | Inventory in progress; semantic review pending |

Each entity group includes nested ownership, mutation, iteration, equality, and
hashing. Operation configuration and returned iterators/results are reviewed with
their owner, even when methods are implemented on Molecule or Reaction in Python.
An inventory entry with no current binding is marked as such, not treated as
missing implementation. No group is complete merely because a representative
method or a related group passed review.

### Initial assignments

Three first-pass reviewers covered Edits, Deltas, and molecule editing/transactions.
The coordinator inventoried exports and cross-group ownership. A fourth reviewer
independently challenged their findings. This first wave is complete, with its
adjudication below; remaining groups are unreviewed. No production fixes were made.

### First-pass source and verification basis

Review code: commit 1c2d0482dec5f459fc5d75e51ecd04f3e84e6b68, with no
production-source edits in the working tree. All reviewers use one isolated copy
of the tracked source plus current AGENTS.md, Python API guide, review skill, and
discussion decisions. The uncommitted policy changes govern the review; proposed
213 APIs are distinguished from currently implemented Rust behavior.

Source manifest SHA-256:
62fa6e09b9ededb0400e7bc5b1f8541cdf885906bfd1f45fad983cf1e49eb5ed.
The snapshot's Python API guide SHA-256 is
c8777f861481fbd98a2324c49600539a2fefb3b89fdfe3a31f9f9d42880f0029;
its python-api-review skill SHA-256 is
d9a1ce04712d895df8c06443daba322f8e72a9e2042ad6dc79cc493988fa2856.

The coordinator built one wheel from that copy with Python 3.13.15 and maturin
build --offline, using the snapshot umol-py/Cargo.toml and its pyproject features
(including depiction). The cp39-abi3 macOS arm64 wheel SHA-256 is
6d1dc4cf3d2e5cbe944360f45d42d6b2eee73a938e91902f5358b6dc5da2cf86.
Probes import the wheel from an isolated extraction directory; the primary
checkout's installed extension was not replaced. No workspace test suite ran.
Concrete probe cases and outcomes are retained with the adjudicated findings;
temporary snapshot/build/report files are disposable.

The initial static export inventory locates 222 direct class registrations and
additional returned iterator classes; field-change and stereo declarations also
use generated bindings. These are inventory leads, not semantic coverage claims.

## Review method and coverage

Use the separate [python-api-review skill](../.agents/skills/python-api-review/SKILL.md).
It reviews Python/Rust operation contracts with independent challenge, rather
than applying the Rust review-cycle's general hygiene areas or defaulting an
unresolved justification to acceptance. Review one bounded operation family per
pass and retain the coverage and dispositions here. The first bounded reviews and
their challenge are complete; the full API remains unreviewed.

Inventory every exported Python constructor, method, property, operator, and
lifecycle transition, including surfaces generated by binding macros. For each,
record the owning public Rust symbol, any Python adaptation and its concrete
usability reason, and whether the binding preserves reference namespaces,
mutation/aliasing, information, failure behavior, and publication guarantees.
An intentionally omitted Rust operation is not a parity violation.

Start with operations involving ids, handles, correspondences, collection
composition, editors, transactions, reactions, and nested mutable values. Then
cover the remaining public surface. Inspect producers and consumers as well as
wrapper signatures. Separate confirmed defects from unreviewed candidates and
legitimate syntax adaptation; a wrapper loop alone does not prove a defect.

For each finding, record the violated semantic contract, a concrete triggering
input and observable consequence, the owning Rust layer, affected Python callers,
and the proposed removal or delegation. Exercise reference-bearing operations
through their consumers: equality of the wrapper's assembled list is insufficient.
Do not add Rust methods merely to justify an unnecessary Python convenience.

## First pass — Edits, Deltas, and molecule editing

Reviewers: review_edits, review_deltas, and review_editing. The separate
challenge_first_pass reviewer checked their claims, strongest defenses, and
source traces, reran their probes, and added cases for nested removal lists,
Delta Remove properties, issued handle kinds, and error structure. All assertions
passed. The coordinator's adjudication follows; it does not approve fixes or
certify every inventoried method independently.

| Finding | Verdict after challenge | Deviation |
| --- | --- | --- |
| P1: Edits.extend | Confirmed for independent creation namespaces; not a resolver defect | Unjustified |
| P2: Deltas.extend | Confirmed operation-ownership issue; no wrong-id result established | Unjustified under the settled removal decision |
| P3: Mutable disconnected Edit/Delta properties | Confirmed, including nested lists and Add/Remove variants | Unjustified |
| P4: Eager iteration | Confirmed; lazy read-only access and errors on invalidation are settled | Unjustified; implementation pending |
| P5: Deltas.normalized_eq preclones | Confirmed equivalent borrowed path without preliminary copies | Unjustified |
| P6: Application copies Edits | Confirmed lifecycle/copy difference; disposition settled after challenge by the consumption clarification | Unjustified; preserve consumption |

Issued handle kind erasure, immutable entry roles,
repeated value conversions, normalization naming, and structured errors remain
unresolved below. Failure to prove a defect does not justify a deviation.

### Coverage

| Family | Inspected surface | Limits |
| --- | --- | --- |
| Edits | New and all eight typed-handle adapters; ConstraintEdit; all 33 Edit variants and both conversion directions; construction, parse/render, comparison, append/extend, indexing/iteration; all creation, removal, update, and molecule-constraint methods | Nested leaf/form/constraint implementations and general DSL/default behavior remain with their groups. Runtime probes cover specific defects, not every variant. |
| Deltas | All eight field-change families; all eight entity-delta families; ConstraintDelta and nine Delta variants; constructors, fields, inverse, equality/repr and conversion arms; Deltas collection, normalization/comparison and returned iterator; reference boundary in lowering | Full matching/lowering algorithm, reaction aggregate behavior, and nested leaf/constraint semantics remain unreviewed. |
| Molecule editing | Three Molecule entry points, all 22 MoleculeEditor methods, and both Transaction rollback methods; current Rust owners and input-copy paths | General molecule construction/access, witness internals, and shared exception policy are separate groups. Proposed 213 APIs are not treated as implemented. |

Direct builder delegation, typed-id conversion selected by parameter, Python
negative indexing, and finite dictionary mapping for ConstraintEdit have concrete
counterparts and purposes. Python-only collection composition and the additional
copy paths below require separate dispositions. The inventory does not approve
all nested types merely because their containing operation was inspected.

### P1 — Edits.extend joins incompatible creation namespaces

Confirmed finding: semantic corruption; unjustified deviation. Python
[edit.rs](../umol-py/src/edit.rs), lines 1344–1360 and 1403–1409, copies entries
without translation. Rust [Edits](../umol-graph-ir/src/ir/edit.rs), lines 440–453,
assigns New ordinals within one sequence and deliberately excludes concatenation.

Executed example: A adds carbon; B adds nitrogen and removes its New(0). B alone
leaves an empty molecule. A.extend(B) applied to an empty molecule leaves nitrogen,
so B's removal targeted A's carbon. Preserving B's reference would leave carbon.

The strongest defense is that raw entries may intentionally use the destination
namespace. That justifies push/append, not treating another independently built
Edits as the same namespace. Tests that compare the flattened list do not check
this consequence. Remove the wrapper extension and its adapters, as already
settled; independent-batch composition remains unresolved in 213. Migrate
extension tests and callers without inventing Python-owned rebasing.

### P2 — Deltas.extend introduces wrapper-owned composition

Confirmed finding: operation ownership; unjustified under the settled removal
decision. Python [delta.rs](../umol-py/src/delta.rs), lines 1758–1788 and
1824–1828, snapshots the operand and implements extension by repeated Rust push.
Rust [Deltas](../umol-graph-ir/src/ir/delta.rs), lines 3164–3207, owns individual
push and collection construction, but no bulk composition operation.

The strongest defense is valid: same-frame ordered concatenation needs no New
rebasing and can be expressed by repeated push. There is no demonstrated wrong-id
result here. The remaining issue is the exposed composition/self-extension
contract being defined in the binding, contrary to the settled boundary.
Remove the wrapper method/adapters; a future Rust bulk operation needs its own
justified contract. Existing extension tests establish behavior, not its necessity.

### P3 — Immutable Edit/Delta properties return mutable disconnected lists

Confirmed finding: misleading mutation; unjustified deviation. The Python
API guide requires property writes to affect the owner or fail. Frozen outer
entries do not make their list-valued properties immutable.

Affected Edit sequences include additions, participant/ligand vectors, outer
removal lists, and nested vectors in removal tuples: [edit.rs](../umol-py/src/edit.rs),
lines 380–427 and 440–537. Delta Add/Remove getters expose variable participants
for dative, aromatic, multicenter, stereo-atom, and stereo-bond families:
[delta.rs](../umol-py/src/delta.rs), lines 702–715, 838–849, 979–990, 1261–1274,
and 1418–1431. Fixed endpoint tuples do not have this particular problem.

Executed Edit probes: clearing AddAtoms.atoms succeeds, but application still
adds carbon; clearing RemoveTopology.atoms succeeds, but application still
removes atom zero. Executed Delta probes append to each of the five Add-family
participant lists; rereading the property and passing the delta through Rust-backed
Deltas both retain the original participants. Independent challenge also reproduced
this for all five Remove families and for a nested donor list in RemoveDativeBonds.

The strongest defense is that these are detached values and the outer object is
frozen. That does not justify accepting apparent writes through a property;
read-only form payloads only protect those forms, not the enclosing lists.
Proposed correction: immutable sequence output recursively for immutable entry
properties, preserving ordinary list acceptance at input. Do not introduce live
mutable entry views. Tests explicitly permitting silent detached writes need to
change with the corrected contract. Other families' properties remain unreviewed.

### P4 — Edits and Deltas iteration materializes every entry before the first next

Confirmed finding: eager conversion defeats iterator laziness. The subsequent
design discussion settles lazy iteration and prefers thin read-only access to
backing Rust entries over copied observations wherever possible.
Rust Edits::iter and Deltas::iter return borrowed slice iterators. Python
[edit.rs](../umol-py/src/edit.rs), lines 1317–1341, and
[delta.rs](../umol-py/src/delta.rs), lines 1730–1755, first convert the entire
container into Python entries with nested payload copies. Merely creating an
iterator pays for the unused suffix. This is a source-traced allocation path,
not a measured workload speedup.

The first pass considered lifetime independence and stable snapshots as defenses.
Neither establishes a requirement for eager materialization. Current Python
containers own their Rust entries and expose append-only mutation; retained input
aliases cannot replace stored entries. Executed probes confirm that current
iterators exclude later appends.

The intended direction is an owner-retaining iterator with a position and initial
length, producing a thin read-only Python accessor only when next requests an
entry. Nested access must also be read-only; freezing an outer wrapper while
returning mutable disconnected lists is insufficient. Creating the iterator must
not convert entries, and yielding one must not copy the remaining entries. Lazy
per-entry copying alone would fix eagerness but would not satisfy the preference
for thin access where backing storage can be accessed safely. Review indexing
and iteration together for their returned-entry contract.

Retaining storage changes its lifetime; conversion on access changes conversion
failure timing. Consumption with outstanding accessors is settled above: it
invalidates them, and later access raises InvalidatedViewError. These lifetime
differences are not reasons to restore eager snapshots. No claim of universally
lower retained memory is made, and no unsafe long-lived Rust reference or public
accessor hierarchy is prescribed. No fix is implemented.

This finding covers Edits.__iter__/EditIter and Deltas.__iter__/DeltaIter. Other
iterator families remain unreviewed. The general review criterion is thin, cheap,
lazy access faithful to Rust; Python's importance to applications does not justify
separate semantics. The Python API guide records that criterion.

### P5 — Deltas.normalized_eq clones before a borrowed comparison

Confirmed finding: unnecessary owned copies; unjustified deviation.
The owned conversion supplied at [delta.rs](../umol-py/src/delta.rs), lines
1844–1849, clones both GraphIrDeltas operands. The shared macro in
[lattice.rs](../umol-py/src/lattice.rs), lines 23–30, then passes references to
Rust Normalize::normalized_eq, which checks structural equality first and clones
for normalization only when needed ([traits.rs](../umol-graph-ir/src/ir/traits.rs),
lines 164–172). The wrapper defeats that early path even for self-comparison.

The strongest defense is reusing one conversion closure for both consuming
normalization and borrowed comparison. Deltas already exposes a borrowed Rust
reference, so that convenience does not require these preclones. Proposed
correction: delegate comparison through borrowed references; review other macro
instantiations individually before generalizing. Evidence is the allocation path,
not an executed allocation counter or timing claim.

### P6 — Application copies Edits instead of consuming it

Six Python entry points clone Edits before consuming Rust execution:
Molecule.apply/tracked_apply and editor apply/tracked_apply/transact/tracked_transact.
See [molecule.rs](../umol-py/src/molecule.rs), lines 286 and 298, and
[transaction.rs](../umol-py/src/transaction.rs), lines 88, 102, 116, and 131.
The cost follows batch contents: Vec entries, participant/removal vectors,
ConstraintEdit handle vectors, and owned expression/constraint trees. It is not
a receiver-sized clone or only an Arc increment. All six executed probes preserve
a reusable, mutable batch.

The first pass left this unresolved because reuse has actual test consumers.
The subsequent consumption clarification settles the disposition: those consumers
can explicitly copy; their existence does not justify copying every application.
Preserve consumption. No PyO3 limitation requiring copying has been established.
Aliases raise ConsumedError after consumption; outstanding views and iterators
raise InvalidatedViewError. Failure ordering and explicit-copy migration still
need concrete design; no fix is implemented.

### Unresolved deviations and handoffs

- **New kind erasure.** Returned Rust handles are typed per entity kind; Python
  returns the same New ordinal type for all eight. Parameter position supplies
  kind on input, but an issued wrong-kind token loses Rust's protection. The
  executed challenge passes an issued bond New(0) as an atom removal handle;
  application removes carbon and leaves nitrogen. This proves lost kind
  protection, not lost batch identity: neither Rust nor Python handles carry the
  latter. The convenience and protection tradeoff remains unresolved; no new
  handle type hierarchy is proposed by this pass.
- **Immutable entries and reusable value operations.** Independently owned Rust
  Edit/Delta values have mutable public fields; Python selects immutable standalone
  entries. That role itself requires justification. Conditional on that role,
  isolating a caller's writable input form has a concrete purpose. It does not
  justify every later copy: Edit comparison reconstructs owned Rust operands,
  and Delta/field-change inverse reconstructs and then copies from an owned
  result. The immutable role and conversion details remain unresolved; possible
  reuse alone cannot justify replacing a consuming operation with implicit copies.
- **Normalization naming.** Rust already has borrowed normalized as well as
  consuming normalize. The former provides the retained-source/owned-result
  behavior selected by Python Deltas.normalize; Python has not invented that
  operation. Its spelling remains unresolved in the lattice group. This does
  not justify the separate preliminary copies in normalized_eq.
- **Error detail.** [error.rs](../umol-py/src/error.rs), lines 88–106, preserves
  exception categories/messages but discards structured Rust reasons and fields.
  An executed HandleOutOfRange case retains only a message, without kind/index/count
  attributes or a cause. Loss of nested RollbackFailed detail is source-traced,
  not a reproduced rollback failure. A compact exception taxonomy is a possible
  justification; no recovery bug or need for one class per Rust variant is claimed.
  The boundary infrastructure group owns this question once for all consumers.
- **Reaction's live deltas.** [reaction.rs](../umol-py/src/reaction.rs), lines
  390–403 and 564–571, exposes mutable deltas and later reconstructs/validates a
  Rust Reaction. This is a source-traced candidate for the reaction aggregate
  group, not an adjudicated finding or completed reaction review.

### Examined behavior matching current Rust

Editor snapshot/tracked_snapshot delegate the existing Rust operations and add
no second Python snapshot. Molecule.edit delegates current Rust borrowed editing;
its shared storage and correspondence construction are Rust behavior. Their
future removal/change belongs to 213, not a newly alleged wrapper invention.

Executed alias probes confirm editor apply/tracked_apply/build/tracked_build
consume the wrapper on success; apply failures also consume it, while the tested
transact failures restore it. Rust also admits RollbackFailed; these probes do not
establish restoration for that branch. Failed snapshot retains a repairable draft;
failed build consumes it. Rollback/tracked_rollback restore the matching original and consume
the journal, including through aliases. These existing Option-backed wrappers
demonstrate that explicit consumption is possible in Python; they do not settle
the lifetime of Edits arguments or every other wrapper.

### Durable reproduction cases

The following cases were executed against the isolated wheel described above:

```python
from umol import AtomForm, Edit, Edits, Molecule, New

a, b = Edits(), Edits()
a.add_atom(AtomForm.parse("C"))
h = b.add_atom(AtomForm.parse("N"))
b.remove_topology([h], [])
assert h == New(0)
assert Molecule().apply(b) == Molecule()
a.extend(b)
assert Molecule().apply(a) == Molecule.parse('{:atoms ["N"]}')

entry = Edit.AddAtoms(atoms=[AtomForm.parse("C")])
entry.atoms.clear()
assert Molecule().apply(Edits([entry])) == Molecule.parse('{:atoms ["C"]}')
removal = Edit.RemoveTopology(atoms=[0], bonds=[])
removal.atoms.clear()
assert Molecule.parse('{:atoms ["C"]}').apply(Edits([removal])) == Molecule()

sequence = Edits([entry])
iterator = iter(sequence)
sequence.append(Edit.AddAtoms(atoms=[AtomForm.parse("N")]))
assert list(iterator) == [entry]
```

For the five Delta probes, Add inputs were: dative donors [0], acceptor 1;
aromatic/multicenter atoms [0,1] with electron counts [1,1]; stereo atom/bond
site 0 with one actual ligand at atom 1 and Th0/Ct0 attributes. These are open
Delta carriers, not published molecule-integrity examples. Appending atom 2
(or its actual StereoLigand) to each returned list leaves both the getter and
Deltas([wrapped_delta])[0]._0 unchanged. Deltas iteration likewise excludes an
entry appended after iterator creation. Challenge repeated the five cases with
Remove variants and clear, with the same disconnected-write outcome. Clearing
RemoveDativeBonds.removes[0][1] likewise left donors [0,1] both in the entry and
after passage through Edits.

For the kind-erasure probe, one batch added carbon, nitrogen, and a bond between
them, then used the returned bond New(0) as the atom argument to remove_topology.
Application left nitrogen. For the error probe, removing atom zero from an empty
molecule raised TransactionError with args containing only
"atom handle 0 is out of range for 0 entries"; kind/index/count attributes were
absent and __cause__ was None.

The editing probe applied one add-N batch to a carbon molecule through all six
entry points, reused it, then appended add-O and verified C/N/O. Alias and error
probes used an out-of-range atom handle 7 after an add, and a duplicate localized
bond to trigger final integrity failure. All 20 recorded lifecycle observations
matched the outcomes above. No whole-workspace test or allocation benchmark was
run. These cases/results are durable; temporary probe scripts and logs are not
required to retain the findings.

## Proposed process improvements

- The strict parity rule is recorded in AGENTS.md and the binding guidance.
  Enforce it at review; further documentation alone does not establish compliance.
- Require each added or changed Python operation to identify its public Rust
  counterpart and explain any deviation during API review. New domain semantics
  must be designed on the Rust owner before binding work.
- Extend the existing data-type contract review to check the binding against the
  Rust contract and its consumers. Include nested references, ownership, errors,
  and operation boundaries, not just matching type shapes.
- Require behavioral parity cases for semantic adaptations. For batch composition,
  use independent creation namespaces and apply the result; test both successful
  results and relevant failure behavior. Avoid tests that merely bless wrapper
  mechanics or copied implementation choices.
- Review changes to the Python guide against the owning Rust semantics. A local
  guide amendment cannot authorize a domain operation absent from Rust or override
  the user's restriction on Python-only API growth.

Use the existing contract and review process; do not introduce another runtime
layer or automatically expose every Rust method in Python. Settle the process
changes and findings before scheduling the broader implementation.

## Implementation plan — P1–P6

Scope: remove wrapper-owned extension; replace disconnected entry properties and
eager entry iteration with thin access; remove Deltas comparison preclones; and
transfer Edits into the six existing application entry points without cloning.
This implements the current Rust operations, not doc 213's future editor and
transaction API. S0 is complete; S1–S6 remain outstanding.

### Contract and scope controls

- Preserve existing variant names, constructor field names, and class-pattern
  inspection. Owned entries and accessors use the same variant classes. Readonly
  reports permissions; copy() produces independent ownership. Nested access must
  obey the same contract. Do not preserve generated storage layout at the expense
  of those semantics.
- Indexing and iteration access original storage. Iterators exclude later appends,
  allocate wrappers only on demand, and remain exhausted once exhaustion has been
  reported. Consumption invalidates active accessors without waiting for them.
- For application, complete argument conversion, acquire required borrows, and
  check both receiver and batch availability before taking either input. Once
  Rust execution begins, the batch remains consumed even on failure. Receiver
  restoration/destruction follows the current Rust method, including its actual
  rollback error contract.
- Retain single-sequence constructors, append, and domain construction methods.
  Never replace independent Edits composition with an unchanged append loop.
  No new composition API, handle-kind hierarchy, general exception redesign,
  inverse/normalization lifecycle redesign, or unrelated binding audit is included.
- Nested forms, constraints, and field-change wrappers are dependencies where
  reached through Edit/Delta access. Change only what is required for faithful
  access, copying, and affected consumers; this does not mark their broader review
  groups complete. Remove obsolete copying adapters when their last caller moves.
- The prototype proves one variant, not every boundary. Before the affected
  subitem, reconcile the exact public symbols and contracts for owned field
  replacement and retained nested accessors, and for entry arguments passed to
  constructors/append. If existing decisions do not determine a behavior, settle
  that specific point with the user; do not silently introduce copies, setters,
  invalidation rules, or additional public types. This plan supplies sequencing,
  not permission to invent missing semantics.

Every subitem includes its affected Rust/Python callers, documentation, and
focused tests. Breaking changes may be temporarily red while editing, but each
subitem should finish green; every stage must finish green. No commits are
authorized by this plan.

### S0 — Independent corrections and error vocabulary

- [x] **S0a — Lifecycle exceptions.** Modules: error.rs, lib.rs, Python package
  exports, transaction.rs. Add ConsumedError and InvalidatedViewError as RuntimeError
  subclasses; migrate existing editor/journal consumed-state errors to the common
  vocabulary. Additive exports plus breaking error refinement, green with callers
  migrated. Verify exact classes, alias behavior, and unchanged Rust failure
  outcomes. [dep: none]
- [x] **S0b — Borrowed Deltas comparison (P5).** Modules: delta.rs, lattice.rs only
  as needed. Delegate normalized_eq through borrowed GraphIrDeltas; leave normalize
  unchanged. Green implementation change. Verify equal, unequal, equal-normal-form,
  and contradiction cases against Rust; inspect the path for preliminary clones.
  Do not generalize other macro users without evidence. [dep: none]
- [x] **S0c — Extension removal (P1/P2).** Modules: edit.rs, delta.rs and their
  Rust/Python tests. Remove both extend methods, EditsExtend, DeltasExtend, and
  ResolvedDeltasExtend. Delete extension/self-extension tests; change the incidental
  single-entry extension in the delta property test to append. Breaking, green
  after migration. Recheck repository callers and preserve construction/append
  coverage, including application of a single sequence with New references.
  [dep: none]

#### S0 outcome — 2026-09-21

Both lifecycle exceptions are public RuntimeError subclasses. Existing
MoleculeEditor and Transaction consumed-state failures now raise ConsumedError,
with the owning type named in the message; their Rust lifecycle is unchanged.
InvalidatedViewError is exported for later accessor implementation and has no
production accessor producer yet.

Deltas.normalized_eq now directly borrows both stored Rust Deltas. Its former
shared-macro use was replaced with two local methods so normalization retains
its existing behavior without imposing owned conversion on comparison. The shared
macro and its other users are unchanged. Source review confirms no preliminary
clone on the comparison path; no timing claim is made.

Removed Edits.extend, Deltas.extend, all three extension adapter types, the unused
Deltas mutable-conversion helper, and extension-specific tests. The incidental
single-entry caller uses append. Repository search found no remaining binding
callers. The maintained single-sequence case adds carbon and nitrogen, removes
the issued nitrogen handle, and verifies that application leaves carbon.

Verification: Python 3.13.15, rebuilt editable extension with graph/depiction
features; 9 transaction, 264 delta, and 18 edit Rust unit tests passed. All 262
Python cases in test_import.py, test_edit.py, test_delta.py, and test_transaction.py
passed. They cover exact exception exports/inheritance and alias failures,
comparison equality/inequality, normal-form fusion, contradictory inputs, and
receiver preservation. Package formatting and staged/unstaged diff checks passed;
the complete S0 diff was reviewed against its scope. No workspace suite, MSRV
check, or timing benchmark ran. Edits application still copies its input until S5;
P3/P4 accessor and iterator changes remain outstanding.

### S1 — Batch ownership and concrete access support

- [ ] **S1a — Edits storage and explicit copy.** Module: edit.rs, Edits. Introduce
  the consumed-state-capable storage needed by accessors, with one checked route
  to its Rust contents. Add copy() preserving entries and creation counters.
  Keep application transfer for S5. Additive API/internal rewire, green. Verify
  copied batches are independent and subsequent creation handles are correct.
  [dep: S0a, S0c]
- [ ] **S1b — Deltas explicit copy.** Module: delta.rs, Deltas. Add copy() with
  independent owned contents; do not add an artificial consuming operation to
  Deltas for symmetry. Additive, green. Verify order, duplicates, and independent
  subsequent mutation. [dep: S0c]
- [ ] **S1c — Scoped owner-backed access.** Modules: entity.rs and concrete
  Edit/Delta storage support. Implement the minimal owner/location and short-borrow
  access needed by the proven variant pattern, including nested permission and
  invalidation propagation. Additive internal support, green; no public generic
  framework. Verify original-storage access, owner retention, and append growth
  without holding references across Python calls. Bring the relevant prototype
  cases into maintained tests; scratch is not a test dependency.
  [dep: S0a, S1a, S1b]

### S2 — Nested payload access

- [ ] **S2a — Entity-form access.** Modules: entity.rs and the eight entity-form
  binding families. Enable owner-backed reads and explicit owned copies for forms
  reached through entries; preserve the existing Python form classes. Migrate
  boundary consumers that currently require an unconditional owned conversion
  to scoped access where their Rust operation borrows. Breaking internal/API-role
  rewire, green with callers. Verify nested reads do not copy whole forms and
  supported writes either affect backing storage or raise TypeError.
  [dep: S1c]
- [ ] **S2b — Nested constraints and sequences.** Modules: constraint bindings,
  sequence-valued payload getters, and their concrete access support. Replace
  mutable disconnected outputs on these entry paths with owner-backed access;
  carry permissions and invalidation through nested containers. Breaking return
  behavior, green with callers. Test retained grandchildren, indexing, lazy
  iteration, invalidation, and explicit copy independence. Do not implement an
  unrelated collection API. [dep: S2a]
- [ ] **S2c — Field changes and ConstraintEdit.** Modules: delta.rs field_change
  generation and edit.rs ConstraintEdit. Provide the nested access/copy support
  required by outer entries, preserving variant and field inspection. Breaking
  wrapper rewire, green with callers. Test both scalar and structured old/new
  payloads plus reference-bearing constraint fields. Existing inverse semantics
  are not redesigned by this item. [dep: S2a, S2b]

### S3 — Edit entries and lazy Edits access

- [ ] **S3a — Edit variant storage (P3).** Module: edit.rs, Edit and its 33
  variants/conversions. Replace generated stored-field wrappers with the shared
  owned/accessor variant interface; implement readonly and copy(). Preserve
  keyword/positional matching and constructor names. Avoid whole-entry rebuilding
  for reads/comparisons where the owning Rust values can be borrowed. Breaking,
  green with constructor/conversion callers migrated. Exercise every variant's
  fields, including nested removal tuples and stereo factors, and preserve exact
  Rust conversion values. [dep: S2c]
- [ ] **S3b — Edits indexing and iteration (P4).** Module: edit.rs, Edits and
  EditIter. Return owner-backed entry accessors; replace eager Vec materialization
  with initial-length lazy traversal. Breaking result ownership, green with tests
  migrated. Verify no entry conversion at iterator creation, one demanded wrapper
  per next, negative indexing, later-append exclusion, read-only nested access,
  and stable exhaustion. Add counters only in test instrumentation; no public
  diagnostic methods. [dep: S3a]

### S4 — Delta entries and lazy Deltas access

Each family subitem below replaces its generated owned-field representation with
the shared variant interface, migrates its conversions/callers, and tests every
Add/Remove/Modify shape it actually exposes. These are breaking rewires, green
when complete. Preserve exact participant factors, field-change values, and
existing inverse results; test readonly, independent copy, and variant matching.

| Subitem | Module/type family | Dependencies |
| --- | --- | --- |
| S4a | delta.rs: AtomDelta | [dep: S2c] |
| S4b | delta.rs: BondDelta | [dep: S2c] |
| S4c | delta.rs: DativeBondDelta | [dep: S2c] |
| S4d | delta.rs: AromaticSystemDelta | [dep: S2c] |
| S4e | delta.rs: MulticenterBondDelta | [dep: S2c] |
| S4f | delta.rs: NoncovalentBondDelta | [dep: S2c] |
| S4g | delta.rs: StereoAtomDelta | [dep: S2c] |
| S4h | delta.rs: StereoBondDelta | [dep: S2c] |
| S4i | delta.rs: ConstraintDelta | [dep: S2b] |

- [ ] **S4j — Outer Delta.** Module: delta.rs, Delta. Connect all nine family
  payloads through thin nested access while preserving outer variant classes and
  fields, including _0. Breaking rewire, green. Verify owner retention through
  outer/family/form/constraint chains and no recursive payload materialization
  merely to inspect an outer variant. [dep: S4a, S4b, S4c, S4d, S4e, S4f, S4g,
  S4h, S4i]
- [ ] **S4k — Deltas indexing and iteration (P4).** Module: delta.rs, Deltas and
  DeltaIter. Use the same lazy read-only contract as Edits; remove eager iterator
  snapshots. Breaking result ownership, green with affected reaction/test callers
  migrated. Verify indexing, initial-length traversal, nested read-only access,
  and no payload clones on access. Preserve S0b's borrowed comparison path.
  [dep: S0b, S4j, S3b]

### S5 — Consume Edits at application boundaries (P6)

- [ ] **S5a — Molecule application.** Module: molecule.rs, apply/tracked_apply.
  Take Edits after argument/availability checks and pass the Rust batch by value;
  remove unconditional cloning. Breaking argument lifecycle, green with callers
  migrated. Update deliberate reuse to explicit pre-application copy(). Verify
  successful and failed application consume the batch, invalidated nested views
  raise the exact error, and invalid input state causes no premature consumption.
  Keep the current Rust molecule receiver semantics. [dep: S3b]
- [ ] **S5b — Editor application and transactions.** Module: transaction.rs,
  apply/tracked_apply/transact/tracked_transact. Check receiver and batch first,
  then transfer Edits without copying. Breaking argument lifecycle, green with
  all callers migrated. Verify both consumed-input orders, successful application,
  ordinary edit failure, final publication where applicable, journal rollback,
  and receiver preservation/destruction against the Rust contract. Do not claim
  guaranteed restoration for RollbackFailed. Include aliases, active iterators,
  retained nested accessors, and surviving explicit copies. [dep: S0a, S5a]

### S6 — Bounded integration and closeout

- [ ] **S6a — P1–P6 integration gate.** Modules: affected bindings, exports,
  tests, and documentation. Check the full changed public surface against the
  approved names/roles; remove unused snapshot/conversion adapters. Review paths
  through reaction and workflow consumers for accidental copies or changed
  semantics without expanding their audit scope. Verify combined entry access,
  copying, application, and invalidation through public Python APIs. Record exact
  coverage and remaining limitations; keep other findings unresolved and other
  groups unreviewed. Update this record/index and resume the broader review only
  after this bounded work is complete. [dep: S0b, S0c, S4k, S5b]

### Verification and critical path

Use Python 3.13 from umol-py/.venv for every PyO3 build and Python test. During
implementation run affected Rust unit filters and focused Python modules/cases;
rebuild the extension before testing changed bindings. Do not run the full
workspace suite after each subitem or stage. A stage gate covers its affected
bindings and migrated consumers, including Python behavior, not compilation alone.

At S6a run formatting/diff checks, the complete umol-py Rust test target, rebuilt
Python tests, strict Clippy and rustdoc for the affected package with the relevant
graph/depiction configuration. Exercise applicable feature-gated binding tests
explicitly. Broaden into other crates only if implementation actually changes
them. Any applicable pinned MSRV check belongs to this final gate, not iteration.

Copy/laziness verification starts with the first real access implementation:
source traces and focused instrumentation must distinguish wrappers/scalars from
payload clones. No timing campaign is required to demonstrate removal of an eager
copy or to compare against an implementation with the wrong semantics.

Critical path: S0a/S0c → S1 → S2 → S3 → S5, with S2 → S4 and S0b joining at S6.
S0b is independent; the Delta family conversions share file ownership and should
not be edited concurrently without isolation. No stage above is optional for
P1–P6 completion. Remaining review groups and independent-batch composition are
separate deferred work, not extra stages in this plan.
