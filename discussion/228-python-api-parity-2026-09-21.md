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
settled operation on the Rust container. Doc 213 resumes after the bounded
lazy-iteration correction below, without waiting for the full audit or an ownership
migration. It retains incremental construction of one Edits sequence and leaves
independent-batch composition unresolved; resolving that composition is not a
prerequisite for the review or removal. If resumed, its assembly and handle
semantics must move from Deltas lowering into Edits methods, with lowering
restructured to use them. Do not implement rebasing in umol-py or preserve the
current methods as compatibility paths.

This document records an urgent review and process correction, including the
first bounded passes below, retained S0 corrections, and the withdrawal of the
S1–S6 ownership/access plan. The sole next correction is lazy iteration with
unchanged yielded-result ownership, followed by a return to doc 213. This is not a
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
are unfinished. Only the bounded lazy-iteration correction precedes returning to
doc 213; the remaining audit and ownership questions are not prerequisites.

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

### Child access after structural replacement — settled 2026-09-21

Structural replacement invalidates all existing child accessors, including
container accessors and their descendants. Subsequent access raises
InvalidatedViewError; newly obtained accessors expose the replacement. Ordinary
mutation within an existing value remains visible through valid accessors.
Explicit independent copies survive. Do not distinguish container accessors from
item accessors, preserve selected children, or implicitly copy their contents.
This settles the retained-child replacement contract for S2; implementation is
still outstanding.

### Entry input ownership — reopened 2026-09-21

The initially approved transfer policy below was reopened because Edit/Delta
payloads reach most entity forms: its lifecycle machinery would be widespread.
The bounded comparison below evaluates that cost against construction copies.
The policy has not been implemented. Further ownership/access implementation is
unscheduled, and S1's preparatory changes were reverted as recorded below. The
former S1–S6 plan is withdrawn; the lazy-iteration proposal does not depend on it.

Entry constructors consume owned payloads into the entry. Batch constructors and
append likewise consume owned entries. Callers use explicit copy() when they need
to retain an input or supply a read-only accessor. A read-only accessor cannot be
consumed in place; construction does not silently copy it.

Convert and check every input, including ownership, availability, and required
borrows, before transferring any of them. Repeated aliases of the same owned
input cannot supply two transfers and must fail before either transfer occurs.
Successful transfer leaves input aliases consumed and invalidates their child
accessors; independent copies remain usable. This is ownership transfer, not
additional molecular validation.

#### Bounded construction comparison — 2026-09-21

**Decision supported:** the tested transfer path needs simple consumed-state
checks and ordinary PyO3 borrows, not a new lifecycle framework. This initially
supported recommending transfer with explicit copies for reuse. That recommendation
understated the delivery and review cost of migrating the connected type surface
consistently; the experiment does not settle that cost. The saving on this small
fixture is modest, and the evidence does not establish an application-wide
speedup. Neither ownership policy is selected. No additional experiment is proposed.

One standalone PyO3 probe compares two implementations of the same path:
AtomForm → Edit::AddAtoms → Edits → Molecule::apply. Inputs are the existing
test_atom_form_parse fixtures C#R(6) and O#n2, giving two forms including a ring
constraint. Both implementations use the actual graph-IR types, owner-backed
read-only observations, and consuming batch application. The copy policy retains
input forms and entries; the transfer policy takes their contents. Initial form
parsing and Python argument conversion are outside the measurement.

| Boundary | Copy policy | Transfer policy |
| --- | --- | --- |
| Construct AddAtoms | Two form clones; 1 Rust allocation requesting 384 bytes | No payload clones; 2 allocations requesting 416 bytes |
| Insert entry into Edits | One entry clone (copies both forms); 2 allocations requesting 1,408 bytes | No payload clones; 1 allocation requesting 1,024 bytes |
| Apply batch | No wrapper payload clones; 26 allocations requesting 2,452 bytes | Same |
| Explicit copy of first stored form | One form clone; no allocation | Same |

Thus transfer avoids four form copies and the 384-byte vector clone on insertion,
but uses a 32-byte temporary vector to retain input borrows during construction.
Both paths make three allocation requests across construction/insertion;
transfer requests 352 fewer bytes. The copy path borrows and copies each input
in turn, since it does not need to preserve inputs against partial consumption.
This fixture's form cloning needs no heap allocation of
its own; allocation count alone would miss the form copies. A thread-local gated
Rust allocator counts allocation/reallocation requests and requested bytes, not
Python allocations, peak memory, or retained memory. Clone counts instrument
explicit wrapper clone calls; the entry clone includes its two contained forms.
No elapsed-time measurements were needed for this structural comparison.

On this macOS arm64 build, AtomForm and Option<AtomForm> both occupy 192 bytes;
Edit and Option<Edit> both occupy 256 bytes. These are Rust storage sizes, not
complete Python object sizes or portable layout guarantees. The copy version
needs consumed state only in Batch; the transfer version also needs it in Atom
and Entry. Each access checks availability, and each transfer takes the Option.
The multi-input constructor borrows all inputs mutably and checks every Option
before taking any. Retaining the borrows makes duplicate aliases fail before
mutation without a separate alias registry or rollback journal.

Checks passed for identical applied molecules, live observations under the copy
policy, invalidation under transfer, invalidation after batch application,
surviving explicit copies, duplicate-input rejection, and an unavailable second
input leaving the first intact. Explicit reuse under transfer required one form
copy (no allocation for this fixture). The probe uses RuntimeError for lifecycle
failures; production exception names were not reimplemented. Both versions share
the same scoped-reader structure. No universal access framework was introduced.

The breadth remains real: production Edit payloads reference all eight entity
forms, all eight field-change families, and ConstraintEdit; their leaf forms and
constraint payloads extend the consuming-input chain. This experiment does not
measure the full migration, nested structural replacement, or all failure shapes.
Those facts limit the conclusion; they do not open another research task.

Reproduction while scratch exists: activate umol-py/.venv and run
`cargo run --offline --manifest-path scratch/228-entry-access/Cargo.toml --bin ownership`.
The probe is src/bin/ownership.rs plus ownership.py in that directory, using
Python 3.13.15 and PyO3 0.29.0 with abi3-py39. All assertions passed. The production
extension and implementation were unchanged. This section preserves the fixture,
method, results, limitations, and recommendation after scratch is deleted.

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

S7a/S7b update: Edits and Deltas iteration now convert one copied entry per next().
Thin entry access remains unscheduled.
The following preserves the original finding and its broader design discussion.

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
accessor hierarchy is prescribed. S7 addresses eagerness only, preserving the
existing copied-result contract.

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

## Completed and withdrawn implementation work

S0 remains complete. S1 was implemented and reverted. The former S1–S6 plan is
withdrawn, not paused for automatic resumption: it coupled local copying/access
corrections to a broad ownership migration without an acceptable delivery and
review boundary. Its entity-form, nested-access, entry-storage, consuming-application,
and integration stages are not executable work items. No replacement migration
schedule is approved. The records below preserve completed work and its reversal.

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

### S1 historical outcome — reverted 2026-09-21

The following records the implementation and its verification before reversal;
it does not describe the current code.

Edits stores optional Rust contents and checks availability on reads and writes.
Edits.copy preserves entries and all eight creation counters; Deltas.copy
preserves order and duplicates. These are the only new public methods. Deltas
has no artificial consumed state. Existing application callers now propagate
unavailable-batch errors; editor application checks the batch before taking the
editor. Application still clones the batch until S5.

Private EditAccess and DeltaAccess retain a concrete batch owner and entry index.
Their scoped readers borrow original storage immutably, including nested payloads,
and release the borrow before returning. EditAccess translates owner consumption
to InvalidatedViewError before invoking the reader. Maintained Rust tests verify
original payload addresses after append growth, owner retention, borrow release,
copy independence, and invalidation. Consumption is exercised internally until
S5 provides its production path; no test-only public method was added.

Indexing uses this scoped support but still materializes Python entries. Python
nested permission propagation and retained child wrappers belong to S2, followed
by thin entry access and lazy iteration in S3/S4. No change to entity.rs was needed
for the concrete batch support; its form integration remains in S2.

Verification: Python 3.13.15; rebuilt editable extension with graph/depiction.
All 450 Rust unit cases across edit, delta, transaction, and molecule modules and
385 Python cases across import, edit, delta, transaction, and molecule modules
passed. Copy tests cover subsequent New ordinals for every entity kind and an
applied batch containing bulk creation and removal. Package formatting and diff
checks passed; the full S1 diff was reviewed. No workspace suite or MSRV gate ran.

### S1 reversal — 2026-09-21

Removed the complete S1 production/test changes: optional Edits storage and its
availability checks, private EditAccess/DeltaAccess support, both batch copy()
methods, application-caller adaptations, and their tests. The six affected source
and test files were restored exactly to their pre-S1 contents. S0 remains intact,
including its lifecycle exception names, removed extension methods, and borrowed
Deltas.normalized_eq comparison. Existing editor/transaction consumption is
unchanged. No new consuming operation or preparatory ownership machinery remains
from S1. The subsequent design discussion and experiment are retained; neither
constitutes approval to resume the migration.

Reversal verification: all six source/test files match the pre-S1 revision
exactly; 442 focused Rust cases and 375 Python cases passed across the edit,
delta, transaction, molecule, and Python import suites. The extension was rebuilt
under Python 3.13.15. Package formatting and diff checks passed. No workspace
suite or MSRV check ran.

## Bounded proposal — lazy Edits/Deltas iteration

This is the sole next correction agreed on 2026-09-21. It removes eager work while
preserving the current ownership of yielded results. It does not implement thin
entry access, eliminate constructor/application copies, or complete P3/P4/P6.
The withdrawn S1–S6 stages are not prerequisites.

### Contract and implementation boundary

- Calling iter() performs no entry conversion or payload copy.
- Each next() converts only the requested entry and returns the same independent
  Python variant/value representation as today. Indexing is unchanged.
- The iterator retains its batch and captures the initial length. Later appends
  are excluded. Separate iterators advance independently; exhaustion is permanent.
- The existing containers are append-only. The iterator briefly borrows the batch
  when converting an entry; no Rust reference survives between Python calls.
- No consuming operation, optional owner storage, new exception, public accessor
  hierarchy, or nested-payload rewrite is introduced. Existing application and
  constructor semantics remain unchanged.

The change is confined to edit.rs (Edits iteration/EditIter) and delta.rs
(Deltas iteration/DeltaIter), their affected callers, and focused tests. Replace
prebuilt vectors of Python entries with a retained batch, position, and initial
length. Reuse the existing entry conversions at next(). Conversion failure occurs
when that entry is requested rather than at iterator creation; advance only after
successful conversion, so failure does not silently skip an entry.

The iterator keeps the batch alive until the iterator is dropped, including any
later-appended entries it will not yield. This replaces eager copied-entry storage
with owner retention; it is not a guarantee of lower retained memory in every
usage pattern.

### Acceptance and stopping point

Verify exact yielded values and variant classes, zero conversions at iterator
creation, one conversion per successful next(), initial-length exclusion,
permanent exhaustion, independent iterators, owner retention, and independence of
yielded results. Use focused Rust tests and rebuilt Python edit/delta tests under
the repository Python 3.13 environment, then package formatting and diff review.
Source inspection or private test instrumentation establishes laziness; no public
diagnostics, timing campaign, workspace suite, or migration framework is needed.

After this correction, return to doc 213. Constructor ownership, thin nested
access, consuming batch application, and the remaining audit findings stay
unresolved and unscheduled. Their unresolved status does not endorse existing
copies or authorize another investigation automatically.

## Implementation plan — bounded lazy iteration

S7 is a new, independent stage; numbering does not reactivate withdrawn S1–S6.
S7a and S7b are complete; S7c remains. Each subitem ends green and
includes its tests and method documentation. There are no preparatory stages or
new public types/methods. The existing Edits.__iter__, EditIter.__next__,
Deltas.__iter__, and DeltaIter.__next__ are the affected Python operations.

### S7 — Lazy iteration with existing result ownership

- [x] **S7a — Edits iteration.** Module: edit.rs, Edits/EditIter and edit_iter;
  tests: edit.rs unit cases and test_edit.py. Replace the eager vector of Python
  entries with a retained Edits owner, current position, and captured initial
  length. Remove edit_iter's eager collection; call the existing Edit conversion
  only from next(). Advance after successful conversion and return exhaustion
  permanently at the captured bound. Migrate affected internal callers and remove
  unused imports. Existing signatures at the Python call site and result
  ownership are preserved; evaluation/error timing becomes lazy. Green with
  callers and tests. Verify empty and populated batches, exact values/classes,
  independent iterator positions, append exclusion before and after exhaustion,
  owner retention after deleting the batch variable, and independent yielded
  results. Inspect the constructor and next() paths to establish zero entry
  conversions at construction and exactly one per successful next(); do not add
  public diagnostics. Run focused edit Rust/Python tests against the rebuilt
  extension. [dep: retained S0; none of withdrawn S1–S6]
- [x] **S7b — Deltas iteration.** Module: delta.rs, Deltas/DeltaIter and delta_iter;
  tests: delta.rs unit cases and test_delta.py. Apply the same concrete iterator
  shape and conversion timing to Deltas, preserving existing Delta variant
  construction and the borrowed normalized_eq path. Remove the eager collection
  and migrate direct iterator callers. Preserve Python call signatures and result
  ownership; green with callers and tests. Cover the same observable cases as
  S7a using actual Delta family payloads, and verify returned nested values retain
  their current independence and permissions. Run focused delta Rust/Python tests
  against the rebuilt extension. Reuse the agreed contract, not a new generic
  iterator framework. [dep: S7a for the reviewed implementation pattern]
- [ ] **S7c — Combined verification and handoff.** Modules: both iterator
  implementations, affected tests, and discussion/status records. Review the full
  change against the scope: no altered indexing, construction, append,
  application, entry storage, or consumption. Run the edit/delta Rust unit filters
  and rebuilt Python edit/delta modules together, package formatting, and
  git diff --check. Record that lazy conversion is implemented while yielded
  payload copying remains. Mark S7 complete and update 213 and the status index
  to resume the mutation-API design. Other Python findings remain open, so do not
  mark the full audit completed. Green closeout; no new API. [dep: S7a, S7b]

Critical path: S7a → S7b → S7c. Both iterator changes are required; no further
stage is implied. Use Python 3.13 from umol-py/.venv for all PyO3 builds/tests and
rebuild with maturin develop before Python checks. The focused combined gate is
sufficient for this two-module change; no workspace suite, benchmark campaign,
or MSRV campaign is scheduled. Constructor ownership, nested access, consuming
application, and broader review remain unscheduled.

### S7a outcome — 2026-09-21

EditIter retains Edits, a position, and the initial length. Edits.__iter__ performs
no entry conversion; next() calls the existing conversion for one entry and
advances only on success. Source review establishes laziness without additional
instrumentation. Yielded values keep their existing independent storage and
permissions. Later appends are excluded and exhausted iterators stay exhausted.
The iterator keeps its owner alive until dropped, including later-appended data.
No optional batch storage, consumed-state checks, or entry-access wrappers were
introduced. Deltas iteration remains for S7b.

All 19 focused edit Rust cases and 30 Python edit cases passed under Python
3.13.15 after rebuilding the extension. New cases cover empty/populated traversal,
variant classes, independent results and iterator positions, owner retention,
append exclusion, exhaustion, and borrow-failure retry without skipping an entry.
Package formatting and diff checks passed. The complete S7a source/test diff was
reviewed against the restored pre-S1 files, separately from the pending S1
reversal. No workspace suite, benchmark, or MSRV check ran.

### S7b outcome — 2026-09-21

DeltaIter now uses the same retained-owner, position, and initial-length structure
as EditIter. Iterator creation performs no entry conversion; each successful
next() converts exactly one Delta through the existing conversion and then
advances. The eager delta_iter collection is removed. Constructor, append,
indexing, normalized_eq, normalization, and copied-result semantics are unchanged.
No consuming-state machinery or shared iterator framework was introduced.

All 265 focused delta Rust cases and 113 Python delta cases passed under Python
3.13.15 after rebuilding the extension. Tests cover empty/populated batches,
Atom/Constraint variant results, independent nested objects and their existing
read-only permissions, iterator independence, retained owners, append exclusion,
permanent exhaustion, and borrow-failure retry. Source review confirms laziness;
the complete diff, package formatting, and diff checks passed. S7c remains for
combined verification and the return to 213. No broader test or timing campaign ran.
