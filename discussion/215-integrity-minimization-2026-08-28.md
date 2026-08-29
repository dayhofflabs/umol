# 215 — Integrity closure and minimization

Status: In Progress
Date: 2026-08-28
Relates: [211](211-relation-frames-and-api-2026-08-26.md),
[214](214-aggregate-frame-semantics-2026-08-28.md),
[data-type guide](../docs/development/data-types.md),
[integrity guide](../docs/development/integrity.md),
[nomenclature guide](../docs/development/nomenclature.md)

## Purpose

The admission-test audit after doc [214](214-aggregate-frame-semantics-2026-08-28.md) S0d found
that the integrity inventory is close to minimal but not yet minimal, and that the implementation
does not yet satisfy the closed-container execution model stated by the development guides. This
document owns the narrow correction: remove or loosen checks that do not protect a coherent
representation, close the remaining public mutation escape, and delete defensive rechecks once the
aggregate boundary establishes the contract.

This work precedes doc 214 S0e and the later frame-transport stages. Those stages should rely on the
final integrity domain and should not preserve error variants, parser behavior, or defensive checks
that this audit has already identified as unnecessary.

`ReactionSpan` construction preserves an explicitly supplied `Modified` tag even when its two
values are semantically equivalent. Collapsing that redundant tag belongs to explicit
normalization or to a producer, such as `superimpose`, whose contract is to derive a standardized
span rather than faithfully retain independently supplied entries.

## Admission result

The inventory in the integrity guide has 29 error rows. Against the current executable public
surface:

- 24 rows pass the admission test without qualification;
- `ReactionIntegrityError::Lhs` passes only while public molecule mutation can compromise the lhs;
- `DuplicateParticipant` passes except for its dative donor/acceptor cross-factor case;
- `ParticipantFrameMismatch` rejects a coherent removal whose explicit local frame can be aligned
  uniquely with its source frame;
- `EntityCountMismatch` has no independently constructible failure state; and
- `DativeBondsParallel` rejects more than the dative identity and lookup contract require.

`ReactionIntegrityError::Lhs` becomes redundant when that mutation surface is corrected and every
published `Molecule` is actually closed.

The audit does not reopen the rest of the inventory. References, positional collection lengths,
the other fixed entity identities, stereo frame and domain checks, structured reaction incidence,
and reaction-span projection checks all prevent a concrete panic, non-unique promised lookup,
incoherent stored value, or incorrect transport.

## Contract sheet

**Type and role:** `Molecule`, `Reaction`, and `ReactionSpan` are closed aggregate graph-IR values.
Their integrity checks establish only the minimum representation contract needed for every ordinary
operation to interpret them coherently.

**Open carrier or operation-issued value:** `MoleculeEntries`, independently assembled `Deltas`,
and `ReactionSpanEntries` are open carriers. `MoleculeEditor` may hold transient invalid state.
Published aggregates are closed and operation-issued correspondences and transactions retain their
existing provenance contracts.

**Intrinsic representation invariants:** stored references resolve; required positional payloads
match their participant frames; entity lookup keys promised to be singular are unique; stereo
frames and values lie in their supported domains; reaction removals use their source entity's
structured incidence and carry an explicit local participant frame that aligns uniquely with the
source frame; and both reaction-span projections form molecules. Atom/bond attribute parallelism
remains structurally guaranteed by construction rather than dynamically checked.

For dative bonds, the corrected intrinsic contract is factor-aware:

- donors are pairwise distinct because they occupy one participant frame;
- the acceptor may also occur as a donor because it belongs to the other distinguished factor; and
- two dative entities may share an acceptor and one or more donors, but their complete
  `(acceptor, donor multiset)` identity keys must differ.

**Contextual properties and supplied context:** correspondence coverage, reaction
materializability, host-dependent DPO conditions, and model inputs remain with the first operation
that requires them.

**Semantic predicates and validators:** chemistry, physical invariants, conformance, constraint
satisfaction, and stereogenicity remain outside integrity.

**Public constructors:** checked and asserted constructors continue to establish the same accepted
domain. Checked mutation of integrity-sensitive fields establishes the same contract transactionally
and leaves the source unchanged on error. The authoritative integrity implementation remains
crate-private for use by constructors, editors, and transactional candidate publication; a closed
published aggregate has no meaningful public `check_integrity` failure domain.

**Conversions and preserved information:** construction and faithful conversion do not normalize,
repair, or discard a representable distinction unless that transformation has a separate explicit
contract. `ReactionSpan::{try_from_entries, from_entries}` therefore preserve an explicit
equivalent `Modified` entry. `Normalize for ReactionSpan` collapses it to `Unchanged`; reframing and
canonicalization inherit that collapse because normalization is their first pipeline step.

**Explicit transformations:** normalization, reframing, remapping, and canonicalization remain
named operations. Trusted transformations preserve integrity by construction and property tests,
not by rechecking their closed inputs and outputs.

**First public consumer requiring each contextual property:** unchanged from the data-type guide.
This work removes only checks of intrinsic properties already established by the aggregate's own
publication boundary.

**Failure, absence, and panic behavior:** independently assembled input receives the owning typed
integrity error; asserted producers panic on the same rejected domain. Closed aggregate operations
do not retain unreachable integrity variants or panic conditions. Checked remapping still returns
`None` for an unsuitable independently supplied correspondence, not for an impossible invalid
closed source.

**Algebraic, preservation, or roundtrip properties:** every public publisher produces an
integrity-valid aggregate; checked mutation is atomic; trusted remapping and canonicalization
preserve integrity; and removing runtime rechecks does not change successful results.

**Rust/Python boundary:** Python live mutation must use the same checked Rust operations and map an
integrity rejection to `ValueError` without publishing the candidate or partially changing the
molecule.

## Integrity-rule corrections

### Remove `EntityCountMismatch`

`MoleculeEntries` supplies atom and bond rows, not an independently supplied `Graph`.
`Molecule::try_from_entries` derives the graph node count from `atoms.len()` and its edge set from
the bond rows. Public molecule mutation can replace forms but cannot resize the graph and attribute
tables independently. The mismatch is therefore producible only by a bug in an internal publisher.

Remove `EntityCountMismatch` and its two branches from `Molecule::check_integrity`. Trusted internal
publishers must preserve graph/table parallelism by construction and cover that preservation in
their own tests; a runtime integrity error is not the implementation-bug boundary.

### Narrow dative cross-entity uniqueness

The current `DativeBondsParallel` check records every `(acceptor, donor)` incidence and rejects a
second dative that shares any one pair. That is stronger than dative identity. Relation coincidence,
family lookup, DSL structural naming, pushout, and correspondence induction all use the complete
acceptor plus donor multiset.

Replace the incidence-wide prohibition with uniqueness of the complete
`(acceptor, donor multiset)` key. Two entries such as donors `{a, b}` to acceptor `x` and donors
`{a, c}` to the same acceptor are coherent distinct entities. Two entries with the same acceptor
and the same donor multiset remain invalid even if their stored donor sequences differ. The error
name should describe duplicate complete identity rather than imply that any shared incidence is
parallel; `DativeBondsIdentical` is the direct peer of `MulticenterBondsIdentical`.

### Make dative participant uniqueness factor-aware

The current dative branch passes the flattened donor-plus-acceptor sequence to
`check_unique_participants`. That correctly rejects repeated donors but also rejects a donor equal
to the acceptor. The latter lies across distinguished birelation factors: the roles, identity key,
incidence, and donor-frame action remain unambiguous.

Check uniqueness within the donor factor only. Continue to reject a repeated donor. Accept an atom
that appears once as donor and once as acceptor; any chemistry judgment about that state remains
lazy.

### Accept compatible reaction removal frames

An overlay `Remove` repeats its source site and structured participant incidence, but its participant
sequence is also an explicit local coordinate frame for its recorded attributes. Requiring that
sequence to equal the source sequence rejects a representation that already carries all information
needed to interpret it. Construction must preserve the supplied frame rather than silently
normalizing the removal into the source frame.

After repeated participants are prohibited, equal structured incidence determines one local
removal-to-source action for every family:

- dative bonds permute donors while the acceptor remains the distinguished other factor;
- aromatic systems and multicenter bonds permute their participant sets;
- noncovalent bonds use the identity or endpoint swap;
- stereo atoms permute their ligand frame; and
- stereo bonds permute within endpoint blocks and may swap the two complete blocks.

For stereo bonds, moving individual ligands between endpoint blocks without swapping the complete
blocks changes structured incidence and remains `IncidenceMismatch`. It is not another frame of the
same entity.

Removal matching first verifies structured incidence, derives the local action, transports the
recorded attributes into the source frame, and only then applies the caller's value relation. The
existing molecule-editor removal comparisons already follow this shape. A non-equivalent transported
old value remains an ordinary delta-continuity failure; it is not malformed representation.

Aggregate reaction reframing still exposes one owning frame and one representative action per
entity id. The lhs frame owns an existing entity and the unique `Add` frame owns a created entity.
Field and constraint deltas are stated in that owning frame. A `Remove` may use another compatible
local frame, so reaction transport composes its derived local-to-owner alignment with the supplied
owner action before transporting the removal sequence and payload. Plain `reframe` derives and
consumes that composition without allocating a per-delta action vector; `reframe_with_action`
returns only the per-entity owning actions.

Consequently, an isolated kind-specific delta is a transport-only consumer of a supplied local
action, but `Deltas` cannot blindly implement transport under `OverlayFrameActions` without the
owning frames held by `Reaction`. Reaction-level transport performs the contextual dispatch.

Remove `ParticipantFrameMismatch` and its exact-order check. `IncidenceMismatch` remains the error
for a different site, participant multiset, distinguished factor, or stereo-bond endpoint-block
assignment.

## Close the remaining molecule mutation escape

The guides say that every public mutation of `Molecule` preserves integrity, but four public Rust
methods violate that statement:

- `aromatic_system_mut` exposes an unrestricted `&mut AromaticSystemForm`;
- `modify_aromatic_systems` replaces every form without validation;
- `multicenter_bond_mut` exposes an unrestricted `&mut MulticenterBondForm`; and
- `modify_multicenter_bonds` replaces every form without validation.

Each form contains a literal electron-count vector whose length must equal its participant count.
The mutation surface can therefore create a `Molecule` rejected by its own integrity check after
publication. The Python aromatic and multicenter electron setters and whole-form `__setitem__`
paths use the same unchecked methods. Other field and constraint setters also depend on the broad
mutable view even when their particular mutation cannot compromise integrity.

Restrict the raw full-form mutation paths to graph IR and add
`try_modify_aromatic_system`, `try_modify_aromatic_systems`,
`try_modify_multicenter_bond`, and `try_modify_multicenter_bonds`, following the existing stereo
mutation contract. Each callback operates on a private candidate, runs the authoritative gate, and
commits only on success while returning the exact `MoleculeIntegrityError` otherwise. Rust internal
callers may retain crate-private raw access. Python live views route every whole-form and field
mutation through the checked operations; failure preserves the old molecule and reports
`ValueError`.

## Remove defensive rechecks

`Reaction` and `ReactionSpan` have private fields and no public mutation after publication. Twelve
post-publication integrity-check expressions are therefore redundant now:

- seven reaction expressions: the DSL check immediately after asserted `Reaction::new`, the check
  in `application_deltas`, and the five expressions reached through reaction canonicalization,
  hashing, and equality; and
- five reaction-span expressions: source and result checks in `try_remap`, plus the three checks in
  canonicalization and canonical-key construction.

Remove those checks and rely on checked/asserted construction plus preservation tests. In
particular, `ReactionSpan::try_remap` retains only the independently supplied correspondence
conditions; remapping the closed source preserves span integrity by construction.

Closing aromatic and multicenter mutation makes seven more checks redundant:

- one molecule expression in `try_remap` and five molecule canonicalization/key expressions; and
- the lhs molecule check inside `Reaction::check_integrity`.

This yields 19 removable defensive check expressions in the final closed model: 12 immediately and
seven after molecule closure. Constructor gates, transactional candidate checks, editor
publication, and reaction-span side projection remain.

The removal includes the public error and documentation consequences. Delete integrity error
variants, conversions, `None` conditions, and panic claims whose only producer was a defensive
closed-source check. Do not retain unreachable diagnostics for API symmetry. Restrict
`Molecule::check_integrity`, `Reaction::check_integrity`, and `ReactionSpan::check_integrity` to
graph IR after their external callers migrate. Constructors retain those authoritative internal
implementations even though a published value can no longer fail them.

## Reaction lhs after molecule closure

Today `ReactionIntegrityError::Lhs` is necessary because a caller can compromise a molecule through
the mutation escape and then pass it to `Reaction::try_new`. Once every public molecule publisher
and mutation preserves integrity, `Reaction::try_new` receives a closed value and the wrapper has no
failure domain.

Remove the lhs recheck and `ReactionIntegrityError::Lhs` after, not before, molecule closure. Update
reaction application, DSL and Python error translation, canonicalization errors, rustdoc, and exact
tests so no branch remains for an impossible invalid lhs.

## Reaction-span construction and derived standardization

`ReactionSpan::try_from_entries` currently calls `normalize_reaction_span_entries` before reference
and side validation. That function converts an `EntitySpan::Modified { lhs, rhs }` to `Unchanged`
when the two forms are equivalent. The conversion is not an integrity check: both variants have a
coherent stored interpretation, and the admission test identifies no panic, ambiguous lookup, or
incorrect transport caused by preserving the supplied tag.

`Modified` makes no assertion beyond carrying the two side values. When those values are equal or
`normalized_eq`, the entry is semantically a no-op and its normal form is `Unchanged`, but the raw
tag remains a representable distinction. Remove `normalize_reaction_span_entries` from checked and
asserted construction. `ReactionSpanEntries` is an open carrier, and its constructors must preserve
an explicitly supplied equivalent `Modified` entry after checking integrity.

`Normalize for ReactionSpan` performs the collapse. `Reframe` and `Canonicalize` also collapse the
entry because each includes normalization as its first pipeline prefix. This is lazy
standardization: independently supplied data is retained until a named operation requests its
semantic normal form.

Operation-generated values may already use that normal form when the operation's result contract
permits it. `ReactionSpan::superimpose` derives a valid span from two molecules and a
correspondence; it is not a faithful conversion of `ReactionSpanEntries`. For paired values that
are `normalized_eq`, it emits `Unchanged` carrying the lhs value. This preserves the exact lhs
projection while the rhs projection retains its documented semantic equivalence under the induced
correspondence. The operation need not normalize the selected lhs payload itself merely to obtain
that classification.

## Retained checks

One reviewed candidate remains deliberately unchanged:

- `AromaticSystemsOverlap` is required by singular `AtomView::aromatic_system` and
  `aromatic_system_id` and by aromatic valence, all of which select one owning system. Supporting
  overlap requires a plural-membership and derived-value redesign, not a cheap integrity
  relaxation.

No other integrity row is reopened by this work.

## Documentation and verification

Update the data-type guide and integrity guide with the corrected dative identity, factor-aware
participant rule, structurally guaranteed graph/table parallelism, actual mutation closure, and the
removal of defensive checks and exact removal-frame rejection. Remove stale claims that dense
remapping can receive an invalid closed source. Keep the nomenclature guide aligned where its
relation identity, integrity, or frame-transport descriptions encode the old rules. Doc 214
separately owns the remaining removal of repeated-virtual-ligand orbit terminology.

Each changed behavior receives its own focused test:

- accept two dative bonds that share an acceptor and donor but have different complete donor sets;
- reject two datives with the same acceptor and donor multiset in different stored frames;
- accept a donor equal to the acceptor while still rejecting a repeated donor;
- accept reordered removals for all six overlay families and transport their payloads into the
  owning frame before comparison;
- reject a stereo-bond removal that changes endpoint-block incidence as `IncidenceMismatch`;
- exercise successful aromatic and multicenter checked mutation, exact length failure, and rollback;
- exercise the corresponding Python success, `ValueError`, and rollback behavior;
- verify that checked and asserted span construction preserve an equivalent `Modified` entry;
- verify that the current reaction-span normalization and canonicalization path collapses it to
  `Unchanged` without relying on constructor normalization;
- verify that `superimpose` emits `Unchanged` for `normalized_eq` paired values while preserving
  exact lhs projection and the documented rhs equivalence; and
- verify the corresponding exact and semantic projection relations after span/reaction
  conversions.

Preservation properties cover graph/table parallelism and integrity across every trusted molecule,
reaction, and reaction-span publisher, remapping, and canonicalization path. They replace defensive
runtime checks rather than merely accompanying them. Existing algebraic properties in doc 214
remain the authority for normalization, reframing, and canonicalization behavior. Reaction
properties additionally establish that a compatible reordered removal has the same materialization,
application, and reversal behavior as its source-framed restatement. Doc 214 S0r-S0t own aggregate
reaction reframing and the composed local-to-owner-to-target action properties after the common
frame-action vocabulary exists.

## Scope boundary and completion

This document does not add chemistry validation, relax aromatic ownership or other fixed entity
identities, redesign canonicalization search, or implement the frame-action work in doc 214. It
changes only the minimum eager contract and its enforcement surface, including acceptance and
transport of compatible removal-local frames.

The work is complete when every remaining integrity check passes the admission test, every public
publisher and mutation preserves the contract, no closed aggregate operation defensively rechecks
that contract, the reaction-span construction behavior is explicit, the Rust and Python error
surfaces contain no unreachable integrity cases, and the normative guides and focused tests agree.

## Staged implementation plan

This plan is inserted between doc 214 S0d and S0e. It inherits doc 214's fifteen-failure transport
ledger and must not add another failure. The slow canonicalization integration and complete
feature-gated property targets are checkpoint gates, not per-subitem gates: each subitem runs its
dependency-local tests and relevant compile, rustdoc, format, and lint checks, while S3 runs the
slow gates once for the complete work item. Every behavioral change below carries its focused test
in the same subitem.

### S0 — Add the checked mutation vocabulary

- **S0a — aromatic-system transactional mutation** (`umol-graph-ir/src/ir/molecule.rs`, molecule
  tests): add `try_modify_aromatic_system` and `try_modify_aromatic_systems` on the existing private
  candidate-and-commit kernel without changing the current raw mutation surface yet. Give each new
  method its own success, exact invalid-reference or electron-count-length failure, and rollback
  cases. **Additive (green).** [dep: none] **Done.**
- **S0b — multicenter-bond transactional mutation** (`umol-graph-ir/src/ir/molecule.rs`, molecule
  tests): add `try_modify_multicenter_bond` and `try_modify_multicenter_bonds` with the same candidate,
  exact-error, and rollback contract. Give each new method its own focused cases rather than sharing
  the aromatic tests as indirect evidence. **Additive (green).** [dep: none] **Done.**

  `Molecule` now exposes the four checked operations over the existing private candidate-and-commit
  kernel. Each singular operation reports an exact unavailable-id or electron-count-length error;
  each family-wide operation reports the exact length error; every failure leaves the original
  molecule structurally unchanged. The focused checked-mutation run passed all ten new cases and
  the five existing stereo/constraint cases. Graph-IR doctests and strict graph-IR library clippy
  passed, and formatting is clean. The slow integration and property checkpoints were not run.

S0 ends with the replacement Rust vocabulary available while every existing caller still compiles.

### S1 — Correct the accepted integrity domain and close mutation

- **S1a — molecule integrity minimization** (`umol-graph-ir/src/ir/molecule/integrity.rs`, molecule
  integrity tests and property scenarios, direct error consumers): remove `EntityCountMismatch` and
  its unreachable graph/table branches; replace incidence-wide `DativeBondsParallel` with
  `DativeBondsIdentical` over the complete `(acceptor, donor multiset)` key; and check duplicate
  participants within the donor factor rather than across donor and acceptor factors. Migrate every
  exhaustive error match. Exact cases accept shared incidences with distinct complete identities
  and donor-equals-acceptor, reject reordered duplicate identities and repeated donors, and retain
  graph/table preservation properties for trusted publishers. **Breaking (red→green).** [dep: none]
- **S1b — aromatic-system mutation closure** (`umol-graph-ir/src/ir/molecule.rs`,
  `umol-py/src/{aromatic,constraint/aromatic}.rs`, Rust and Python tests): restrict
  `aromatic_system_mut` and `modify_aromatic_systems` to graph IR, migrate every Python live-view
  field setter, constraint mutation, and collection `__setitem__` through the S0a checked operation,
  and map rejection to `ValueError`. Cover successful field and whole-form replacement, invalid
  electron-count length, exact exception, and unchanged owner state after failure. Internal graph-IR
  transforms may retain the crate-private raw path. **Breaking (red→green).** [dep: S0a]
- **S1c — multicenter-bond mutation closure** (`umol-graph-ir/src/ir/molecule.rs`,
  `umol-py/src/{multicenter,constraint/multicenter}.rs`, Rust and Python tests): apply the S1b
  closure and migration to multicenter forms, including every live field, constraint, and
  `__setitem__` path. Cover the same exact success, `ValueError`, and rollback contract independently.
  **Breaking (red→green).** [dep: S0b]
- **S1d — compatible reaction removal frames** (`umol-graph-ir/src/ir/reaction.rs`,
  `ir/reaction/integrity.rs`, reaction property scenarios, DSL diagnostics): remove
  `ParticipantFrameMismatch` and the exact-order validation pass while retaining entity-specific
  structured-incidence checks and diagnostic precedence. Use the existing frame-aware editor
  equality/application kernels to interpret compatible local removal frames. Exact cases cover
  accepted nonidentity frames for all six overlay families, nonuniform transported payloads,
  created-entity owners, rejection of a stereo-bond cross-block move as `IncidenceMismatch`, and
  agreement of materialization, application, and reversal with an owner-framed restatement.
  **Breaking (red→green).** [dep: none]
- **S1e — reaction-span raw/normal separation** (`umol-graph-ir/src/ir/{reaction_span,delta,
  canonicalize}.rs`, reaction-span and canonicalization tests): remove
  `normalize_reaction_span_entries` from checked/asserted construction. Make the current private
  reaction-span normalization path normalize both `Modified` sides and collapse them when they are
  equivalent, so canonicalization does not depend on constructor normalization; leave the later
  public `Normalize`/`Reframe` trait migration to doc 214 S0s. Preserve `EntitySpan::superimpose` as
  the deriving kernel that emits `Unchanged(lhs)` for equivalent sides. Test raw constructor
  preservation, normalization/canonicalization collapse, exact lhs projection, semantic rhs
  projection, and conversion roundtrips. **Breaking (red→green).** [dep: none]

S1 ends with every accepted aggregate still coherent, every public molecule mutation transactional,
compatible removal-local frames admitted, and span construction representation-preserving. The
inherited transport ledger is unchanged.

### S2 — Remove closed-source defenses and unreachable errors

- **S2a — molecule closed-source cleanup** (`umol-graph-ir/src/ir/molecule/{remapping}.rs`,
  `ir/canonicalize.rs`, graph-IR/graph/Python canonicalization callers and tests): remove the source
  integrity branch from dense molecule remapping and the five molecule canonicalization/key checks.
  Remove the now-unreachable integrity arm from `MoleculeCanonicalizeError`, its conversions and
  Python mapping, and stale `None`, error, and panic documentation. Replace malformed-published-value
  tests with constructor-boundary and trusted-publisher preservation cases. **Breaking
  (red→green).** [dep: S1a, S1b, S1c]
- **S2b — reaction closed-source cleanup** (`umol-graph-ir/src/ir/reaction.rs`,
  `ir/reaction/integrity.rs`, `ir/canonicalize.rs`, reaction DSL and Python canonicalization callers and
  tests): remove `ReactionIntegrityError::Lhs` and its check after molecule closure; delete the
  redundant post-`Reaction::new`, application, canonicalization, hashing, and equality checks and
  every error branch they alone made reachable. Remove the integrity arm from
  `ReactionCanonicalizeError` and migrate Rust/Python error translation and exact tests while
  retaining operation-specific materialization, DPO, and contradiction failures. **Breaking
  (red→green).** [dep: S1a, S1b, S1c, S1d]
- **S2c — reaction-span closed-source cleanup** (`umol-graph-ir/src/ir/reaction_span.rs`,
  `ir/canonicalize.rs`, Rust/Python canonicalization callers and tests): remove source and result
  integrity checks from `ReactionSpan::try_remap` and the three canonicalization/key checks. Remove
  the unreachable integrity arm from `ReactionSpanCanonicalizeError`, update Python translation and
  rustdoc, and cover unsuitable correspondences separately from preservation of a valid closed
  source and result. **Breaking (red→green).** [dep: S1e]
- **S2d — internal integrity gates and publisher preservation** (`umol-graph-ir` aggregate
  constructors, editors, transformations, integration/property tests): restrict
  `Molecule::check_integrity`, `Reaction::check_integrity`, and `ReactionSpan::check_integrity` to
  graph IR after their external consumers are gone. Keep them as the single implementations behind
  checked/asserted construction and transactional publication. Migrate direct public-validator tests
  to exact checked-boundary cases and preservation properties covering editors, remapping,
  canonicalization, reaction/span conversion, composition, reversal, and application publishers.
  **Breaking (red→green).** [dep: S2a, S2b, S2c]

S2 ends with the closed-container execution model enforced: malformed open carriers fail at their
publication boundary, while operations on published aggregates contain no defensive integrity moat
or unreachable public diagnostic.

### S3 — Documentation, cross-surface verification, and closeout

- **S3a — normative and API documentation** (`docs/development/{data-types,integrity,
  nomenclature}.md`, graph-IR rustdoc, Python docs): update the integrity inventory and failure
  justifications, dative identity and factor rules, checked mutation names and rollback behavior,
  compatible removal-local frames, crate-private integrity gates, closed-source remapping and
  canonicalization contracts, and raw-versus-normal reaction-span construction. Remove every stale
  public error, panic, and defensive-validation claim. **Additive (green).** [dep: S2d]
- **S3b — checkpoint verification and lifecycle closeout** (workspace, doc 215, status index): run
  formatting, graph-IR unit and doctests, the affected property modules, strict relevant clippy,
  workspace all-target compilation, and the Python extension and Python tests under the repository
  Python 3.13 environment. Run the slow canonicalization integration, complete feature-gated
  property target, and workspace test gate once here; require exact agreement with doc 214's
  inherited ledger and no new failure. Record the evidence, mark every subitem **Done**, then change
  doc 215 and its status row to `Completed` only when the full scope is implemented. **Additive
  (green).** [dep: S3a]

Critical path: S0a/S0b → S1b/S1c → S2a/S2b → S2d → S3a → S3b. S1a, S1d, and S1e may proceed after
the current doc 214 S0d baseline and converge before their S2 consumers. No stage is deferrable:
doc 214 S0e and later frame transport rely on the complete corrected integrity and error surface.
