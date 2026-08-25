# 209 — Level and hash semantics

Status: Proposed
Date: 2026-08-25
Relates: [186](186-molecule-canonicalization-2026-08-05.md),
[208](208-canonicalization-scaling-2026-08-24.md),
[data-type guide](../docs/development/data-types.md),
[nomenclature guide](../docs/development/nomenclature.md)

## Purpose

Canonicalization optimization exposed unresolved relationships among structural hashing,
normalization, equivalence, canonical keys, and description levels. This document separates public
semantic changes from optimizations that must preserve the current semantics. The public boundary
is settled below; allocation and key representation remain non-blocking implementation questions.

The immediate scope is the shared `Canonicalize` implementation for `Molecule`, `Reaction`, and
`ReactionSpan`. This includes private effective-level inspection for their stored forms; an
independent redesign of any of the three representations remains out of scope.

`Molecule: Hash` is present on the development trunk but has not appeared in a release. The hash
scope must therefore be settled and implemented before the next release. Neither its current
implementation nor its numeric output carries a released compatibility obligation.

## Settled public boundary

Canonicalization is complete-only in the public API. The public `DescriptionLevel`,
`Molecule::description_level()`, `canonicalize_by`, `canonical_hash_by`, and `canonical_eq_by`
introduced during doc 208 are to be removed from Rust and Python. No public `normalize_by`,
`equiv_by`, or `equiv_under_by` operations are added.

The public `Canonicalize` contract retains:

```text
canonicalize
canonicalize_with_correspondence
canonical_hash
canonical_eq
```

The canonicalizer retains a private `CanonicalizeLevel` with `Topology`, `Constitution`,
`Structure`, and `Full`. Private inspection selects the lowest level required by each stored
aggregate:

- a molecule requires the greatest level of its populated entity families and constraints;
- a reaction requires the greater level of its lhs and every delta;
- a reaction span requires the greatest level of its populated entity-span families and constraint
  spans.

Atoms and localized bonds have base level `Topology`; non-stereo overlays have base level
`Constitution`; stereo entities have base level `Structure`. An explicit constraint delta or span,
a `ModifyConstraint` delta, or any carried entity form with a non-empty inline constraint store
requires `Full`. An entity span inspects every carried side, including both sides of `Modified`.
These inspectors are canonicalization machinery, not properties or operations of the graph-IR
representations. There is no public way to request, inspect, or force a level.

Each aggregate's `canonicalize` and `canonicalize_with_correspondence` use the private level selected
for their operand. `canonical_eq` uses the greater private level required by its two operands. These
reductions are exact only when they preserve the complete operation's value and failure behavior.

This also confines normalization and equivalence to their complete existing meanings. Complete
canonicalization normalizes the complete returned aggregate. `equiv` and `equiv_under` continue to
compare complete graph-IR values; they do not acquire level-dependent variants.
`canonical_eq` retains its existing complete comparison, contradiction, and integrity behavior;
this document reopens its allocation strategy, not its relation.

For `Molecule`, `canonicalize_with_correspondence` must satisfy:

```text
canonicalize_with_correspondence(x).0 == canonicalize(x)
x.remap(canonicalize_with_correspondence(x).1) == canonicalize(x)
```

It need not reproduce the correspondence that an inaccessible forced-`Full` path would select
among symmetry-equivalent alternatives. That comparison is useful only as internal evidence.

## Stability boundary

Canonicalization is deterministic within a fixed umol version and context, but canonical forms,
selected correspondences, and canonical hash values are not stable between versions. They are not
persistent identifiers.

The private level selection, typed comparison schema, and search implementation may therefore
change between versions without violating the public contract. A future level-selecting API, if
ever justified, would be a new operation with newly stated semantics; the complete-only API makes
no advance promise about such projections.

## Canonical hash

`canonical_hash` remains public and hashes the complete canonical aggregate. The accepted baseline
is explicit complete canonicalization followed by the aggregate's structural `Hash`. It may
materialize the canonical aggregate.

Avoiding that materialization later is an optimization, not a new operation. A virtual canonical
view or streaming hash is acceptable only if it preserves the result and failure behavior of
canonicalizing and then hashing under the same release, context, and hasher. The hash remains a
non-persistent, collision-prone value with no cross-release numeric stability guarantee.

There is no public level-constrained structural or canonical hash. Replacing the derived structural
`Hash` implementation is unnecessary unless later measurements justify it for the complete
operation. Allocation, key representation, and complete-operation performance belong to doc 208
and do not block this API correction.

## Staged implementation plan

Every stage ends with a green workspace. This plan removes the accidental level surface and installs
the private dispatch needed by doc 208; it does not select or implement a hash-specific
optimization.

### S0 — Establish private aggregate dispatch

#### S0a — Restore the private canonicalization level and leaf inspection

**Module:** `umol-graph-ir/src/ir/canonicalize.rs` and its unit tests.

Add private `CanonicalizeLevel::{Topology, Constitution, Structure, Full}` and private
inspection for deltas and entity spans. Give each entity family its base level, inspect both carried
forms of a modified span, and raise any inline or explicit constraint change to `Full`. Keep the
public `DescriptionLevel` surface temporarily so this subitem is additive and green.

**Tests and evidence:** Use module-local `rstest` tables covering every delta family and operation,
every entity-span position, both sides of `Modified`, all inline constraint stores, and explicit
constraint deltas and spans. Assert the private containment order without adding a public query.

**Change class:** additive private infrastructure (green).

**Dependencies:** [dep: none]

#### S0b — Derive the effective level of each aggregate

**Module:** `umol-graph-ir/src/ir/canonicalize.rs` and its unit tests.

Add private aggregate inspectors. A molecule takes the maximum required by its entity collections
and molecule constraints. A reaction takes the maximum of its lhs and ordered deltas. A reaction
span takes the maximum across all entity-span collections and its constraint spans. Binary
selection takes the maximum required by the two operands. Inspection is total and does not
validate, normalize, apply deltas, or materialize a reaction span.

**Tests and evidence:** Cover the four levels for each aggregate, including a reaction whose lhs
and delta require different levels and a modified entity span whose two sides require different
levels. Assert that dense id remapping and delta inversion preserve the selected level.

**Change class:** additive private infrastructure (green).

**Dependencies:** [dep: S0a]

#### S0c — Route complete aggregate operations through the private level

**Module:** `umol-graph-ir/src/ir/canonicalize.rs`, its unit tests, and aggregate canonicalization
properties.

Route `canonicalize` and `canonicalize_with_correspondence` for `Molecule`, `Reaction`, and
`ReactionSpan` through their selected private levels. Route each `canonical_eq` through the greater
level required by its operands. `canonical_hash` continues to canonicalize and then structurally
hash, so it inherits unary dispatch without a separate path.

**Tests and evidence:** Compare selected and forced-full internal results for representative
molecules, reactions, and reaction spans at every private level. Assert exact canonical aggregates,
unchanged contradiction and integrity behavior, aggregate-specific correspondence transport laws,
complete canonical-hash invariance, and asymmetric `canonical_eq` cases where only one operand
requires the higher level. Do not freeze the forced-full correspondence itself. Assert
`canonical_hash(x) == hash(canonicalize(x))` directly for all three aggregates.

**Change class:** semantics-faithful private dispatch (green).

**Dependencies:** [dep: S0b]

### S1 — Remove the public level surface

#### S1a — Retire the Python consumers

**Module:** `umol-py/src/canonicalize.rs`, molecule bindings, package exports, and Python tests.

Remove Python `DescriptionLevel`, `Molecule.description_level`, and every `canonicalize_by` and
`canonical_eq_by` method from molecule, reaction, and reaction-span bindings. Remove their package
exports and signature inventory entries. The Rust provider still exists during this subitem, so the
Python package can return to green before the Rust surface is removed.

**Tests and evidence:** Rebuild the extension with the repository Python 3.13 environment and run
the focused import, molecule, reaction, and reaction-span tests. Assert the complete operations and
their existing exceptions rather than replacing removed projection tests.

**Change class:** breaking Python API removal (green after its caller and test migration).

**Dependencies:** [dep: S0c]

#### S1b — Retire the Rust level API

**Module:** `umol-graph-ir/src/ir/canonicalize.rs`, `molecule.rs`, and the graph-IR root exports.

Remove `canonicalize_by`, `canonical_hash_by`, and `canonical_eq_by` from `Canonicalize` and all
aggregate implementations. Remove public `DescriptionLevel`, `Molecule::description_level`, and the
root re-export. Convert every retained internal level parameter to private `CanonicalizeLevel` and
remove reaction-projection helpers that existed only for the public reduced operations.

**Tests and evidence:** Compile the graph-IR library and confirm that `Canonicalize` exposes only
`canonicalize`, `canonicalize_with_correspondence`, `canonical_hash`, and `canonical_eq`. The next
subitem migrates test and benchmark callers, so this breaking subitem may be red within S1.

**Change class:** breaking Rust API removal (red until S1c).

**Dependencies:** [dep: S1a]

#### S1c — Migrate Rust tests, properties, and benchmarks

**Module:** graph-IR canonicalization unit and property tests and
`umol-graph-ir/benches/canonicalize.rs`.

Remove representation-level and public projection tests. Retain complete canonicalization, hash,
equality, renumbering, contradiction, integrity, and correspondence laws. Keep forced-level
comparisons only inside the canonicalization module as private implementation evidence. Replace
benchmark calls to removed methods with the complete public operations needed for the target to
compile; doc 208 owns the expanded performance matrix.

**Tests and evidence:** Run graph-IR unit tests, the feature-gated canonicalization properties, and
the canonicalization benchmark build. Search Rust source outside the canonicalization module for
the removed enum and methods.

**Change class:** caller and verification migration (restores Rust green).

**Dependencies:** [dep: S1b]

#### S1d — Align the living development guides

**Module:** `docs/development/data-types.md`, `docs/development/nomenclature.md`, and
`docs/development/property-tests.md`.

Remove the description-level entry, suffix inventory member, retired-name rule, public trait
methods, and level-specific property claims. Describe aggregate canonicalization, canonical hash,
and canonical equality only through their complete operations. Retain structural domain and
incidence-level terminology where it remains independently meaningful.

**Tests and evidence:** Search the living guides for `DescriptionLevel`, `description_level`, and
the removed `_by` names; every remaining occurrence must describe another current API rather than
the retired surface. Run `git diff --check`.

**Change class:** documentation migration (green).

**Dependencies:** [dep: S1c]

### S2 — Verify and close the semantic correction

#### S2a — Run the cross-language verification gate

**Module:** the graph-IR canonicalization surface, `umol-py`, and all workspace callers.

Run formatting, graph-IR unit and feature-gated property tests, the canonicalization benchmark
build, graph-IR linting, the Python 3.13 build and tests, and workspace test and lint gates. Audit
the public Rust and Python inventories for the complete-only surface. `CanonicalizeLevel` may occur
only in canonicalization-private implementation and module-local tests.

**Tests and evidence:** Every gate passes; complete canonicalization and transport properties remain
green; no public documentation, export, binding, example, or benchmark requires a level selector.

**Change class:** verification only (green).

**Dependencies:** [dep: S1d]

#### S2b — Reconcile the discussion records

**Module:** docs 208 and 209 and `discussion/000-status.md`.

Record the implemented private-dispatch and API-removal outcome. Resume doc 208 from the
complete-only API, with canonical-hash measurement and any allocation work remaining there. Mark
doc 209 `Completed` only after S2a and keep cross-version stability explicitly unsupported.

**Tests and evidence:** Discussion links and statuses agree, the doc-208 next action names current
private and public surfaces, and `git diff --check` passes.

**Change class:** closeout documentation (green).

**Dependencies:** [dep: S2a]

### Dependency summary

The critical path is
`S0a -> S0b -> S0c -> S1a -> S1b -> S1c -> S1d -> S2a -> S2b`.
No stage is deferrable within doc 209. Canonical-hash, key-allocation, orbit-pruning, and prefix-
pruning performance work remains in doc 208 and is not a completion condition here.
