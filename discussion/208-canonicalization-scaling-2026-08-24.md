# 208 — Canonicalization scaling

Status: In Progress
Date: 2026-08-24
Relates: [186](186-molecule-canonicalization-2026-08-05.md),
[205](205-mapping-test-corpus-2026-08-20.md),
[207](207-reaction-network-spike-2026-08-24.md),
[209](209-normalization-canonical-semantics-2026-08-25.md),
[211](211-relation-frames-and-api-2026-08-26.md),
[214](214-aggregate-frame-semantics-2026-08-28.md)

## Purpose

The reaction-network experiment has exposed aggregate molecule canonicalization, rather than
matching or network bookkeeping, as its dominant cost. This document scopes an immediate
investigation of that behavior before the classification-network campaign in doc 207. The goal is
to determine why the current exact canonicalizer scales poorly on small concrete hydrogen-rich
molecules and to select a principled correction.

This is the third downstream expansion from the molecular-data experiment: storage led to the atom-
mapping corpus, which led to reaction-network generation, which has now supplied a realistic
canonicalization workload. That provenance is a reason to keep the work sharply bounded.
Canonicalization is already a central graph-IR operation and should be corrected there if the
problem is general; the reaction-network crate must not acquire a corpus-specific identity scheme.

## Semantic correction from docs 209 and 214

Doc 209 implemented private effective-level dispatch and proposed removing description levels and
level-selecting canonicalization methods from the public Rust and Python APIs, but it was superseded
before its aggregate semantics and removal stages began. Completed doc
[214](214-aggregate-frame-semantics-2026-08-28.md) implemented that remaining scope: aggregate
normalization and reframing are complete, and public canonicalization is complete-only. The
completed S1 and S2 work below records the interim surface by which the branch reached this point;
it is not the current public API. Nested description levels remain private exact-search machinery.
S3b is now the next subitem and must preserve complete canonicalization, equality, correspondence,
and hash semantics.

## Current evidence

The extended carbon--hydrogen closures in doc 207 use concrete atoms and ordinary bonds, with no
stereo entities, relational constraints, aromatic systems, multicenter bonds, or noncovalent
entities. The profiled generation phases are:

| Case | Atoms in seed | Products canonicalized | Total | Product canonicalization |
| --- | ---: | ---: | ---: | ---: |
| Ethane | 8 | 855 | 0.368 s | 0.349 s |
| Propane | 11 | 17,929 | 143.173 s | 142.791 s |

Matching and reaction construction together account for less than 0.2% of the propane run.
Canonicalization accounts for 99.73%, averaging about 8 ms per produced derivation, versus about
0.4 ms in the ethane closure. The increase per call is far larger than the change in molecule size.
This is sufficient to stop treating canonicalization as incidental overhead.

The current implementation provides a concrete starting explanation. Complete `Molecule`
canonicalization uses graph-IR's typed individualization/refinement search. For the `Full` level it
unconditionally disables automorphism-orbit pruning and does not use prefix pruning. When a refined
entity cell remains non-singleton, it recursively individualizes every candidate in that cell.
Although the backend's canonical labels do not select the final graph-IR numbering, the current
branch-order function still invokes the automorphism backend at each unresolved search node to order
those candidates.

Consequently, a symmetric molecule can explore many equivalent graph-IR branches while repeatedly
calling the backend for branch-order hints. This observation locates the first problem in the
library-ordered complete-key search; it is not evidence that nauty's canonical labeling is itself
the slow operation. Doc 186 disabled full-search orbit pruning because stereo and other covariant
frames can make a projected id automorphism insufficient to prove equality of complete leaf keys.
The reaction-network values do not carry those features, so the reason for the general fallback may
not apply to this workload.

The generator currently canonicalizes every successful derivation before it can determine whether
the product is a new flask. Repeated or symmetry-related rule applications therefore multiply the
cost, but they do not create it: the per-call scaling must be understood before application-orbit
reduction or parallel execution can be evaluated honestly.

The existing `umol-graph-ir/benches/canonicalize.rs` Criterion target already supplies a seven-case
aggregate benchmark. It measures incidence construction, complete remapping, and topology,
constitution, structure, para-stereo structure, and full canonicalization. A 2026-08-24 optimized
quick run produced the following operation ranges for its two feature-free ordinary-bond cases:

| Case | Topology | Constitution | Structure | Para structure | Full |
| --- | ---: | ---: | ---: | ---: | ---: |
| Naphthalene | 65.4-66.0 us | 64.6-66.2 us | 160.0-163.3 us | 162.9-165.0 us | 196.3-197.0 us |
| Disconnected rings | 90.7-92.2 us | 91.7-92.0 us | 7.19-7.27 ms | 7.28-7.35 ms | 9.24-9.26 ms |

For the disconnected rings, incidence construction at every level is about `0.70 us` and complete
remapping is about `3.16 us`. The roughly 78-fold structure/topology and 100-fold full/topology
ratios therefore arise in the search rather than in carrier construction or final transport. This
case already demonstrates the blanket-pruning defect independently of the reaction-network loop.
The existing benchmark does not retain reaction-network products or report search-tree counters;
those are the two additions needed for the investigation rather than a new benchmark framework.

## Settled design boundary

### Private canonicalization level

Doc 209 restores private `CanonicalizeLevel::{Topology, Constitution, Structure, Full}` inside the
canonicalizer and removes every public level selector. Private inspection selects the lowest level
required by a molecule, reaction, or reaction span from the entity families and constraints carried
by that representation. Complete unary operations use that level, while complete binary equality
uses the greater level required by either operand.

This is an exact search reduction, not a public projection. `IncidenceLevel` remains a distinct
public carrier selector because its `Full` value contains every structural entity but no
constraints.

### Stereo-preserving automorphism orbits

The structure-level counterexample does not justify disabling orbit pruning for every structure
search. Let `G` be the automorphism group of the structural incidence graph and let

```text
H = { g in G | g applied to every stereo site, ligand frame, and configuration preserves stereo }
```

Orbit pruning at structure level is sound over the orbits of `H`. A non-trivial restriction occurs
only when a residual structural automorphism acts non-trivially on a stereo frame. Equivalent
substituents at a prochiral site are the immediate example, but global symmetry may also exchange
stereo entities or equivalent ligand environments. A stereo center whose ligands and related
stereo entities are already distinguished by structural refinement does not require exhaustive
branching merely because a stereo entity is present.

The first implementation may retain the backend's full occurrence-node generators, apply each
generator to the stereo atoms and bonds, discard generators that do not preserve the stereo
payload, and build pruning orbits from the retained generators. The accepted generators form a
subgroup of `H`, so this may under-prune when a product of individually rejected generators lies in
`H`, but it cannot discard a required branch. A full stabilizer calculation is a later optimization
only if the conservative subgroup leaves material search work. Full-level pruning for genuinely
constrained molecules additionally requires the corresponding action on normalized constraints;
the description-level reduction avoids imposing that unresolved case on unconstrained molecules.

## Investigation boundary

The first investigation must isolate complete canonicalization from the reaction-network loop. It
needs representative slow products, not only the ethane and propane seeds, because electronic-state
changes and disconnected intermediates may alter the residual symmetry classes. For each benchmark
case it should record at least:

- incidence construction, normalization, refinement, backend, leaf-key, and final-remapping time;
- refinement calls, backend calls, leaves visited, and branches removed by each sound pruning rule;
- residual cell sizes after structure refinement; and
- timings and canonical results at topology, constitution, structure, and full levels.

The existing private canonical-search statistics are a suitable starting point for research
instrumentation. They are not a reason to expose a permanent diagnostics API. Slow cases should be
retained as canonicalization benchmarks and correctness fixtures once their behavior is
understood.

Every candidate change must preserve exact canonical forms and the correspondence transport law
under dense entity renumbering. A source-to-canonical correspondence may change to a
symmetry-equivalent representative; the public contract does not freeze one arbitrary permutation
inside a symmetry class. Performance alone cannot justify substituting a different equality
relation. Comparison against the current implementation is useful for fixtures that finish, while
independent permutation laws remain necessary so the current slow search is not treated as an
infallible oracle.

## Candidate directions

The investigation should distinguish the following possibilities rather than combining them into
one optimization patch:

1. **Private-level reduction.** Use the canonicalizer's private aggregate-level selection to avoid
   searching empty higher-level key sections. Verify exact canonical forms and correspondence
   action against the unreduced search on bounded dense renumberings and retained slow benchmark
   cases.
2. **Stereo-preserving automorphism action.** Retain the occurrence-node action, restrict pruning
   to generators that preserve stereo atoms and bonds, and measure whether the resulting sound
   subgroup recovers the useful orbits. Compute the full stabilizer only if generator filtering is
   measurably insufficient.
3. **Typed prefix pruning.** The search already has a prefix-pruning seam but the production
   canonicalizers supply a predicate that never prunes. A branch may be rejected only by a typed
   key prefix shared by every completion of its current ordered partition. Prefix construction stops
   at the first component whose source row, referenced target id, stereo frame, or constraint
   position is unresolved. The branch is worse only when this guaranteed prefix is
   lexicographically greater than the incumbent prefix of the same length. This conservative rule
   may miss pruning opportunities but cannot reject an improving completion.
4. **Complete-operation materialization.** After private-level reduction, compare complete
   canonicalization, canonicalization with correspondence, canonical hash, and canonical equality.
   A hash-specific path is justified only if materializing and structurally hashing the canonical
   aggregate remains material. It must reproduce `hash(canonicalize(x))`; hashing a different
   canonical key would change the operation rather than optimize it.

Topology-level canonicalization is an especially useful diagnostic for the current network domain:
its molecules contain no description above topology, so full canonicalization has exactly the same
non-empty typed key. Private-level reduction handles this as canonicalizer dispatch rather than a
reaction-network-specific identity path. The unreduced search remains useful during verification
but is not a distinct semantic operation.

## Relationship to the corpus campaign

The classification corpus contains 159 manifested networks through six non-hydrogen atoms and 25
additional seven-atom GraphML networks without matching rows in the current manifests. Those are
independent jobs and can be run as six to eight single-threaded processes once canonicalization has
an acceptable execution path. Internal reaction-network multithreading is therefore not an
immediate requirement.

The campaign output remains the doc-207 artifact pair: durable QRS GraphML plus faithful endpoint
mapping records stored through the doc-201 substrate or another explicit database. The
canonicalization investigation must not redesign that artifact, QRS selection, endpoint symmetry
classification, or the atom-mapping objective.

## Entity lookup index cost

Doc [211](211-relation-frames-and-api-2026-08-26.md) relocates entity lookup from graph-core's
`find_by_participants` to the entity-family types, keyed by whatever integrity establishes as each
family's uniqueness key. That settles the semantics but leaves the index shape open, and the shapes
differ per family:

| family | uniqueness key | index available today |
| --- | --- | --- |
| aromatic system | any member atom | incidence is already exactly right — systems are atom-disjoint |
| multicenter bond | atom set | incidence gives candidates; an atom may be in several |
| noncovalent bond | unordered pair | incidence on one endpoint plus a check |
| dative bond | acceptor and donor set | incidence on the acceptor plus a donor-set check |
| stereo atom, stereo bond | the site alone | none — incidence indexes ligands as well as sites |

Stereo is the case worth measuring. Its union incidence index is over-inclusive for lookup: a ligand
atom belongs to every stereo entity it participates in, and adjacent stereocentres routinely make
each other ligands, so `incident(atom)` returns several entries and a filter to the site position
exists only to undo that. A direct site-to-id map would remove the filter at the cost of a second
structure to build and keep valid across construction, remapping, and compaction.

The DSL namespace already answered the same question for itself with a `HashMap` keyed by what it
calls the canonical participant key, so there is a precedent to measure against rather than a design
to invent.

This is a cost question, not a semantic one; doc 211 settles the semantics and does not depend on
the outcome.

## Staged implementation plan

Each stage ends with a green build. Search counters are temporary research instrumentation, and
comparison switches remain private test infrastructure; this work does not add a public profiling
API. Before the S6 verification gate, all statistics types, result fields, and counter updates must
be removed from ordinary builds. Any retained instrumentation must compile only under `cfg(test)`.
The implementation order separates feature reduction, orbit pruning, and prefix pruning so their
correctness and performance effects remain attributable.

### S0 — Retain benchmark cases and establish the search baseline

#### S0a — Retain representative scaling cases **Done**

**Module:** `umol-graph-ir/benches/canonicalize.rs`.

Retain a small set of slow, structurally distinct products from the ethane and propane reaction
networks alongside the existing naphthalene and disconnected-ring cases. The set must include a
feature-free connected molecule, a feature-free disconnected molecule, and a symmetry-heavy
electronic-state variant. Record the molecule size, projected description level, and provenance
needed to reproduce each value without depending on the experimental reaction-network crate.

This is additive benchmark infrastructure with no public API change.

**Tests and evidence:** Parse or construct every retained value through an existing graph-IR
boundary, assert that it is valid, and run the existing Criterion groups at all four description
levels. The benchmark must report the retained cases individually rather than hiding them in an
aggregate.

The benchmark target retains the exact native molecule DSL for these three values. The benchmark
ids describe the scaling behavior; provenance remains attached data:

| Benchmark id | Reaction-network provenance | Atoms | Bonds | Components | Populated level |
| --- | --- | ---: | ---: | ---: | --- |
| `feature_free_connected` | Extended C/H propane seed, flask 0 | 11 | 10 | 1 | Topology |
| `feature_free_disconnected` | Extended C/H propane product, flask 72 | 11 | 9 | 2 | Topology |
| `symmetry_heavy_radicals` | Extended C/H ethane product, flask 99 | 8 | 1 | 7 | Topology |

Every value is parsed through `MoleculeDsl`, checked for graph-IR integrity, and asserted to contain
only atoms and localized bonds: no overlays, stereo entities, inline constraints, or molecule-level
constraints. The connected seed anchors the ordinary connected shape, while the two retained
products cover disconnected and symmetry-heavy electronic-state shapes without a dependency on the
experimental reaction-network crate.

An optimized Criterion quick run on 2026-08-24 produced these individual operation ranges:

| Benchmark id | Topology | Constitution | Structure | Para structure | Full |
| --- | ---: | ---: | ---: | ---: | ---: |
| `feature_free_connected` | 95.5-98.7 us | 99.1-101.1 us | 2.890-2.892 ms | 2.834-2.876 ms | 3.728-3.836 ms |
| `feature_free_disconnected` | 91.3-93.9 us | 92.8 us | 3.989-3.997 ms | 3.899-3.952 ms | 5.092-5.136 ms |
| `symmetry_heavy_radicals` | 47.0 us | 47.4-47.7 us | 10.77-10.92 ms | 10.48-10.55 ms | 15.57-15.84 ms |

The results confirm the intended separation: topology and constitution remain below `0.11 ms`,
while empty higher-level searches cost up to `15.84 ms`. These timings are baseline evidence, not
test thresholds.

**Depends on:** none.

#### S0b — Extend private search accounting **Done**

**Module:** `umol-graph-ir/src/ir/canonicalize.rs` and
`umol-graph-ir/src/ir/canonicalize/tests.rs`.

Extend the existing private `CanonicalSearchStats` with backend-call counts and enough residual-cell
information to distinguish refinement, backend ordering, leaf comparison, orbit pruning, and prefix
pruning. Keep the accounting confined to canonical-search internals and module-local tests. It may
remain on the shared internal search path during the investigation, but it is temporary and subject
to the S6a release-code removal gate. Criterion continues to measure timings through the public
operation; do not add diagnostics to `Canonicalize` or expose a public search result wrapper for the
benchmark harness.

This is additive private instrumentation with no public API change.

**Tests and evidence:** Add example-based `rstest` cases for a singleton partition, a symmetric
partition, and an orbit-pruned partition. Assert exact counter relationships only where they are
algorithmic invariants; timings remain benchmark evidence rather than test assertions.

`CanonicalSearchStats` now records initial residual entity-cell sizes, refinement calls,
branch-order calls, backend calls, visited leaves, leaf comparisons, and both pruning counts. The
private branch-order callback reports whether it invoked the backend, allowing the collector to
count ordering calls without double-counting automorphism output reused for orbit pruning. The
module-local table establishes these exact base cases:

| Case | Residual cells | Refinements | Branch orders | Backend | Leaves | Comparisons | Prefix-pruned | Orbit-pruned |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Singleton | `[]` | 1 | 0 | 0 | 1 | 0 | 0 | 0 |
| Symmetric, backend-ordered | `[2]` | 3 | 1 | 1 | 2 | 1 | 0 | 0 |
| Symmetric, orbit-pruned | `[2]` | 2 | 1 | 1 | 1 | 0 | 0 | 1 |

The accounting remains temporary release-path code during the investigation and is still subject
to the mandatory S6a removal gate.

**Depends on:** S0a.

#### S0c — Lock the semantic baseline **Done**

**Module:** `umol-graph-ir/src/ir/canonicalize/tests.rs` and the canonicalization property-test
module.

Record the exact canonical aggregate for each retained benchmark case and the source-to-result
transport law for its returned correspondence. Do not freeze an arbitrary correspondence
representative when multiple symmetry-equivalent representatives satisfy that law. Record baseline
timings and search counters in this document before changing dispatch or pruning.

This is additive verification work with no public API change.

**Tests and evidence:** Use example tests for retained expected aggregates and property tests for
dense entity renumberings. The properties must compare the transported aggregate and canonical key,
not the raw permutation chosen inside a symmetry class.

The retained native DSL values are already the exact canonical aggregates at topology,
constitution, structure, and full levels. The module-local table parses each value independently,
asserts exact aggregate and selected-key equality at all four levels, and checks that the full
source-to-result correspondence transports the source to that aggregate. It asserts semantic
transport through `equiv_under`; it does not record or compare a particular correspondence image.

The property suite now generates independently shuffled dense bijections in all eight entity
namespaces instead of exercising only the reverse numbering. Over integrity-valid generated
molecules, it compares the complete canonical aggregate, hashes, and level-specific canonical
equality after transport, and checks the correspondence law for successful complete
canonicalization.

Together with the S0a timings, the retained-case search counters establish this pre-change
baseline. Topology and constitution have identical accounting for these feature-free values, as do
structure and full:

| Case | Level | Residual cells | Refinements | Branch orders | Backend | Leaves | Comparisons | Prefix-pruned | Orbit-pruned |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Connected | Topology / constitution | `[6, 2, 2, 6, 2, 2]` | 6 | 5 | 5 | 1 | 0 | 0 | 10 |
| Connected | Structure / full | `[6, 2, 2, 6, 2, 2]` | 271 | 127 | 127 | 144 | 143 | 0 | 0 |
| Disconnected | Topology / constitution | `[4, 4, 2, 4, 4]` | 6 | 5 | 5 | 1 | 0 | 0 | 10 |
| Disconnected | Structure / full | `[4, 4, 2, 4, 4]` | 329 | 137 | 137 | 192 | 191 | 0 | 0 |
| Radicals | Topology / constitution | `[6, 2]` | 7 | 6 | 6 | 1 | 0 | 0 | 16 |
| Radicals | Structure / full | `[6, 2]` | 2,677 | 1,237 | 1,237 | 1,440 | 1,439 | 0 | 0 |

No retained counter collector or printing path was added. The values above were read from the S0b
private search result during the focused baseline run; ordinary tests retain only semantic
assertions, and the S6a release-code removal gate remains unchanged.

**Depends on:** S0b.

### S1 — Introduce the representation-owned description level

#### S1a — Add `DescriptionLevel` and molecule description inspection **Done**

**Module:** `umol-graph-ir/src/ir/molecule.rs`, the graph-IR crate root, and their unit tests.

Add the public `DescriptionLevel` enum with the settled ordering and add
`pub fn Molecule::description_level(&self) -> DescriptionLevel`. Determine the result from
collection lengths: topology alone, any non-stereo overlay, any stereo atom or bond, and any inline
or molecule-level constraint. The operation neither validates nor interprets values.

This is an additive public API change.

**Public surface:**

- add `DescriptionLevel::{Topology, Constitution, Structure, Full}`;
- derive `Clone`, `Copy`, `Debug`, `PartialEq`, `Eq`, `PartialOrd`, and `Ord`;
- add `Molecule::description_level()`;
- introduce no constructor, conversion failure, or new error type.

**Tests and evidence:** Use an `rstest` table covering the four levels, including both inline and
molecule-level constraints. Test the documented order and `min`/`max` behavior directly. Add a
property test stating that dense id renumbering does not change the description level.

`DescriptionLevel` now lives beside `Molecule` and is re-exported from `umol_graph_ir::ir` with the
settled derived containment order. `Molecule::description_level()` checks all non-stereo overlay and
stereo collections, every entity family's inline constraint store, and the molecule-level
constraint store without validating or interpreting their values.

The module-local `rstest` table covers ordinary topology, each non-stereo overlay family, stereo
atoms and bonds, every inline entity-constraint family, and a molecule-level constraint. A separate
table pins the adjacent ordering and `min`/`max` laws. The property suite promotes the existing
independent dense permutations of all eight entity namespaces to a shared molecule strategy and
verifies that description inspection is invariant under those renumberings; the existing
canonicalization properties continue to pass against the same strategy.

**Depends on:** S0c.

### S2 — Replace the canonicalization-specific selector

#### S2a — Migrate the Rust canonicalization API **Done**

**Module:** `umol-graph-ir/src/ir/canonicalize.rs`, aggregate implementations, crate exports,
benches, tests, and all Rust callers in the workspace.

Replace `CanonicalizeLevel` with `DescriptionLevel` in `Canonicalize`, `Molecule`, `Reaction`, and
`ReactionSpan` canonicalization APIs. Remove `CanonicalizeLevel`; do not retain an alias or a second
selector with identical variants. `IncidenceLevel` remains unchanged. Migrate the graph-IR crate and
non-Python Rust callers here; the binding crate completes the workspace migration in S2b.

This is a breaking public rename. The graph-IR crate must be green at the end of the subitem, and
the whole stage becomes green when S2b has migrated the binding crate. There is no completed stage
with duplicate public names.

**Public surface:**

- remove `CanonicalizeLevel`;
- use `DescriptionLevel` in `canonicalize_by`, `canonical_hash_by`, and `canonical_eq_by`;
- preserve the existing return values and failure behavior of every canonicalization operation.

**Tests and evidence:** Migrate existing example and property tests without weakening their laws.
Compile-check all Rust call sites and run the graph-IR test and benchmark targets.

`CanonicalizeLevel` has been removed from the Rust API without an alias. `DescriptionLevel`, owned
by the molecule representation, now selects the level used by `Canonicalize` for `Molecule`,
`Reaction`, and `ReactionSpan`; `IncidenceLevel` is unchanged. Every non-Python Rust caller, unit
test, property test, and canonicalization benchmark has been migrated. The development data-type
and nomenclature guides now use the representation-owned term as well. The Python wrapper remains
the deliberate S2b boundary.

The migration passes `cargo check --workspace --exclude umol-py`, the complete graph-IR unit suite
(`6096` passed, `3` ignored), all `18` focused canonicalization properties, the canonicalization
benchmark build, and graph-IR Clippy with all targets and the property-test feature.

**Depends on:** S1a.

#### S2b — Migrate the Python selector and feature query **Done**

**Module:** `umol-py/src/canonicalize.rs`, molecule bindings, package exports, type stubs or Python
surface declarations, and Python tests.

Expose the Rust `DescriptionLevel` as Python `DescriptionLevel`, remove Python
`CanonicalizeLevel`, and add the molecule description-level query. Boundary conversions remain
`from_rust` and `to_rust`; do not add a second Python spelling or retain a deprecated duplicate
during this experimental migration.

This is a breaking Python rename plus one additive query.

**Public surface:**

- remove Python `CanonicalizeLevel`;
- add Python `DescriptionLevel` with the four Rust variants;
- add `Molecule.description_level()` returning `DescriptionLevel`;
- preserve all existing canonicalization return and exception behavior.

**Tests and evidence:** Update import, enum-conversion, molecule-query, and canonicalization tests.
Build against the repository Python 3.13 environment and run the focused Python suite.

The Python binding now exports `DescriptionLevel` as the only description selector and exposes
`Molecule.description_level()`. The enum conversion remains a direct total correspondence with the
Rust variants, and every molecule, reaction, and reaction-span canonicalization method now accepts
the renamed value without changing its result or exception surface. Package exports and public
operation signatures have been migrated; `CanonicalizeLevel` is absent from `umol-py`.

The migration passes the `umol-py` Rust suite (`1632` passed, `2` ignored), a fresh editable build
against Python 3.13.15, the complete Python suite (`1324` passed, `2` skipped), and `umol-py` Clippy
over all targets with warnings denied.

**Depends on:** S2a.

### S3 — Lower empty description levels exactly

#### S3a — Route aggregate operations through their effective level **Done**

**Module:** `umol-graph-ir/src/ir/canonicalize.rs` and aggregate canonicalization implementations.

Doc 209 implemented private canonicalization-level inspection and aggregate dispatch. Completed doc
[214](214-aggregate-frame-semantics-2026-08-28.md) supplies the molecule, reaction, and
reaction-span frame semantics and complete-only public API. The three aggregate roots now inspect
all semantic content that can raise their private `DescriptionLevel`: molecule entities and
constraints, the reaction left-hand side and every delta variant, and both sides of every reaction-
span entry. Unary canonicalization, canonicalization with correspondence, and canonical hashing use
that effective level. Binary canonical equality uses the greater effective level of its operands.

Exact table cases cover all four levels for `Molecule`, `Reaction`, and `ReactionSpan`. They verify
that effective dispatch agrees with forced `Full` canonicalization, correspondence transport
reconstructs the canonical result, canonical hashing hashes that result, and level selection is
renumbering-invariant where the aggregate admits direct remapping. A separate table uses
asymmetric topology/constitution, constitution/structure, and structure/full operands for all
three aggregate roots. Each pair is unequal in both operand orders at the greater level but would
compare equal at the lower level, so the cases distinguish greater-level dispatch from an
incorrect lower-level or left-biased implementation. Reaction-span coverage also includes a
`Modified` entry whose right side alone raises the effective level.

**Dependency satisfied by:** completed doc 214.

**Done.** Focused graph-IR tests pass for all aggregate dispatch, asymmetric equality, entity-span,
and delta-level cases.

#### S3b — Verify the reduction under renumbering and retained workloads **Done**

**Module:** canonicalization property tests and `umol-graph-ir/benches/canonicalize.rs`.

Apply the existing dense-renumbering generators to feature-free and partially featured molecules.
With doc 214's public level removal complete, benchmark the retained cases through
`canonicalize`, `canonicalize_with_correspondence`, `canonical_hash`, and `canonical_eq` for an
equal remapping and a structurally unequal input. The existing pre-change explicit-level results
provide the forced-full and lower-level timing baselines without retaining a public forcing API.

This is additive verification and benchmark work with no public API change.

**Tests and evidence:** Property-test equality of canonical aggregates, hashes, and transported
results under renumbering. For binary equality, generate operands with different populated feature
levels. Record search counters and complete-operation Criterion results in this document; do not
assert wall-clock thresholds in tests. Use the result to decide whether a separate canonical-hash
materialization optimization is justified.

The focused molecule properties now distinguish feature-free values from constitution- and
structure-bearing values rather than relying on the broad aggregate generator to sample those
domains. Independently shuffled dense correspondences preserve the complete canonical aggregate
and hash, and each operation-issued correspondence transports its own source to that aggregate.
The adjacent topology/constitution, constitution/structure, and structure/full equality cases
apply an independent dense renumbering to the higher-level operand and remain unequal in both
operand orders.

An optimized Criterion quick run on 2026-08-29 measured the complete public operations on the three
retained scaling cases. The equal comparison uses a reverse dense renumbering; the unequal
comparison adds one disconnected oxygen atom so that the inputs differ structurally.

| Case | Canonicalize | With correspondence | Hash | Equal comparison | Unequal comparison |
| --- | ---: | ---: | ---: | ---: | ---: |
| `feature_free_connected` | 93.553-93.659 us | 92.947-93.440 us | 94.266-94.381 us | 182.20-184.26 us | 185.89-186.15 us |
| `feature_free_disconnected` | 89.246-89.467 us | 89.678-91.409 us | 90.582-90.757 us | 175.49-176.36 us | 177.09-177.48 us |
| `symmetry_heavy_radicals` | 45.994-46.327 us | 46.152-46.356 us | 46.828-46.906 us | 86.400-86.407 us | 88.773-88.973 us |

The public feature-free path selects topology search and reproduces its private search accounting:

| Case | Residual cells | Refinements | Branch orders | Backend | Leaves | Comparisons | Prefix-pruned | Orbit-pruned |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Connected | `[6, 2, 2, 6, 2, 2]` | 6 | 5 | 5 | 1 | 0 | 0 | 10 |
| Disconnected | `[4, 4, 2, 4, 4]` | 6 | 5 | 5 | 1 | 0 | 0 | 10 |
| Radicals | `[6, 2]` | 7 | 6 | 6 | 1 | 0 | 0 | 16 |

Canonical hashing is within `0.6-1.5 us` of canonicalization on these cases. Its remaining cost is
the shared canonical search rather than a distinct hashing or materialization bottleneck, so S3b
does not introduce a separate canonical-hash optimization. The feature-free canonicalization range
is now `45.994-93.659 us`, compared with the pre-dispatch full-search range of `3.728-15.84 ms`.

**Dependency satisfied by:** completed doc 214.

### S4 — Restore sound orbit pruning at structure level

#### S4a — Establish source-generator sufficiency **Done**

**Module:** the graph-IR automorphism adapter and canonicalization search internals.

Use the existing source-node generator projection. Molecule integrity makes every complete ligand
value within a stereo frame distinct. Mapping each ligand's atom id while retaining its ligand kind
therefore determines a unique permutation between the mapped source frame and the target frame.
The occurrence-node action carries no additional stereo information and is not retained.

No graph-core or public API change is required.

**Tests and evidence:** The adapter tests retain exact source-orbit, source-canonical-label, and
source-generator projection checks. S4b's atom- and bond-stereo cases demonstrate that the projected
source action uniquely determines the required frame action, including a whole-stereo-site exchange.

**Depends on:** S0b and S3b.

#### S4b — Filter generators by stereo preservation **Done**

**Module:** graph-IR canonicalization stereo-action helpers and their unit tests.

Apply each source generator to every stereo entity, site, and ligand-bearing atom. Preserve ligand
kinds, derive the unique action from the mapped and target frames, and retain the generator only
when the transported normalized configuration equals the target configuration. The feature-free
fast path returns the generator vector without visiting it.

This is additive private algorithm work with no public API change.

**Tests and evidence:** The exact table covers a feature-free non-trivial generator, an ordinary
distinct-ligand stereocenter, an explicit-hydrogen prochiral transposition, a global symmetry
exchanging two stereo sites, and a stereo-bond endpoint-block exchange. The prochiral transposition
is rejected because it changes the tetrahedral configuration; both non-trivial global exchanges are
retained. Structure search filters the automorphism output already needed for backend branch
ordering, without an additional backend call. Orbit pruning remains disabled until S4c rebuilds
source orbits from the retained subgroup.

**Depends on:** S4a.

#### S4c — Build safe subgroup orbits and enable pruning

**Module:** graph-IR canonicalization search.

Build source-entity orbits from the retained generator subgroup and enable orbit pruning for
`Structure` searches. Feature-free searches use the full structural orbit information. Constrained
`Full` searches remain exhaustive in this scope. If the retained-generator subgroup under-prunes a
representative stereo case, record that fact; do not silently introduce a full stabilizer
algorithm into this subitem.

This changes private search behavior while preserving exact results.

**Tests and evidence:** Compare pruning enabled and disabled on bounded exhaustive cases, forward
and reverse candidate orders, retained stereo examples, para-stereo examples, and dense
renumberings. Assert exact canonical aggregates and transport laws, then record leaves, backend
calls, orbit-pruned branches, and timings for the retained benchmark cases.

**Depends on:** S4b.

### S5 — Enable conservative typed prefix pruning

#### S5a — Represent guaranteed complete-key prefixes

**Module:** graph-IR typed canonical-key and search internals.

Add a private representation for the ordered prefix guaranteed to be shared by every completion of
a search partition. Define comparison against an incumbent prefix of the same length. Prefix
construction must stop at the first unresolved source row or referenced id; it must not fill unknown
positions with guessed minima.

This is additive private algorithm infrastructure with no public API change.

**Tests and evidence:** Use small hand-enumerated partitions to assert that the reported prefix is a
prefix of every descendant leaf key. Cover an empty prefix, a partially fixed topology prefix, and a
branch whose guaranteed prefix is already worse than the incumbent.

**Depends on:** S0b, S3b, and S4c. Prefix pruning does not require orbit pruning, but this ordering
keeps the measured effects attributable.

#### S5b — Extend prefixes across typed feature sections

**Module:** graph-IR aggregate key builders and canonicalization search.

Extend guaranteed-prefix extraction through constitution and structure sections when their source
rows and references are fixed. Treat stereo frames and constraints conservatively: stop before the
first unresolved covariant reference or normalized constraint position. The implementation may
therefore produce a shorter prefix at higher levels without changing its meaning.

This is additive private algorithm work with no public API change.

**Tests and evidence:** Enumerate all completions of bounded topology, overlay, stereo, and
constraint partitions and assert that the extracted value prefixes every completed key. Include
disconnected and symmetry-heavy cases.

**Depends on:** S5a.

#### S5c — Connect the production pruning seam

**Module:** graph-IR canonicalization search and aggregate canonicalizers.

Replace the production no-op prefix predicate with the guaranteed-prefix comparison at every
canonicalization level. Keep an internal disabled path for tests and benchmark attribution, not as
a public algorithm selector.

This changes private search behavior while preserving exact results.

**Tests and evidence:** Compare prefix pruning enabled and disabled on retained benchmark cases, bounded
exhaustive partitions, both candidate orders, and dense renumberings. Assert exact canonical
aggregates and correspondence transport, then record prefix-pruned branches, leaves, backend calls,
and timings. The stage is not complete merely because the counter becomes non-zero.

**Depends on:** S4c and S5b.

### S6 — Integrate, measure, and dispose of follow-on work

#### S6a — Run the graph-IR correctness and performance gate

**Module:** the complete graph-IR canonicalization surface.

Run formatting, focused unit tests, feature-gated canonicalization property tests, graph-IR linting,
and the canonicalization benchmark. Reconcile the public API inventory against the implementation
and confirm that no diagnostic type or duplicate selector escaped into the public surface. Before
running the performance gate, remove search accounting from the ordinary build: delete it or place
the statistics type, instrumented result fields and entry points, and every counter update behind
`cfg(test)`. The release benchmark must exercise the uninstrumented production search.

This is verification only.

**Tests and evidence:** All graph-IR canonicalization tests and properties pass; every retained
benchmark case preserves the exact aggregate and transport law; benchmark results and search
counters are recorded in this document. A source audit and release build confirm that
`CanonicalSearchStats` and its updates are absent from the release-compiled path. This is a
completion gate, not deferred cleanup.

**Depends on:** S3b, S4c, and S5c.

#### S6b — Run the Python and workspace regression gate

**Module:** `umol-py` and all workspace users of canonicalization.

Activate `umol-py/.venv`, confirm Python 3.13, rebuild the extension, run the Python tests, and run
the workspace test and lint gates. Search source, tests, examples, and documentation for the removed
public `DescriptionLevel` and `_by` surface. `DescriptionLevel` may remain only in
canonicalization-private code and module-local tests.

This re-verifies the complete-only selector contract across language boundaries after the search
changes.

**Tests and evidence:** The Python and workspace gates pass with complete-only canonicalization and
no public selector. Any fixture or generated artifact that embeds the removed surface is either
migrated or explicitly disposed of.

**Depends on:** S6a.

#### S6c — Re-run the reaction-network workloads

**Module:** `experimental/reaction-network` benchmark or reporting path; no production dependency
from graph-IR to the experimental crate.

Re-run the ethane and propane extended-rule closures with the same seeds, rules, and limits used for
the current evidence. Confirm identical flask, adjacency, and transformation counts before
interpreting timing differences. Reversibility need not be repeated if the network artifact is
unchanged.

This is external evidence, not a graph-IR API change.

**Tests and evidence:** Record end-to-end generation time, product-canonicalization time, and calls
per produced derivation, and pair them with the final module-local counter runs on retained graph-IR
benchmark cases. A larger case is justified only if the propane run no longer identifies the
remaining bottleneck.

**Depends on:** S6a.

#### S6d — Reconcile documentation and lifecycle status

**Module:** this document, `discussion/000-status.md`, and the nomenclature and development guides.

Document the implemented private `CanonicalizeLevel` dispatch, complete-only API, pruning boundary,
benchmark evidence, and any deliberately unresolved behavior. Mark this document Completed only
after the public inventory and all verification gates agree with the implementation, including the
S6a removal of release-path search accounting.

This is documentation and lifecycle work.

**Tests and evidence:** Documentation examples name only current APIs, `git diff --check` passes,
and the status row describes completed rather than proposed scope.

**Depends on:** S6a, S6b, and S6c.

### Dependency summary

The remaining critical path is
`S3b -> S4a -> S4b -> S4c -> S5a -> S5b -> S5c -> S6a -> S6b/S6c -> S6d`.
S5 follows S4 deliberately so benchmark deltas remain attributable even though their private
infrastructure is largely independent. No remaining implementation stage is optional within this
scope.

## Potential issues after the current scope

The following concerns remain separate work. They should be reconsidered using the final S6
measurements rather than being folded into the pruning implementation:

1. **Direct backend numbering.** Determine whether a typed adapter can use backend canonical labels
   as the final graph-IR numbering while preserving typed ordering, covariant frames, and exact
   aggregate semantics, instead of using the backend only for branch ordering and automorphisms.
2. **Incremental canonicalization.** Determine whether a product obtained by a local reaction edit
   can reuse refinement, partition, or canonical-label information from its source without making
   molecule identity depend on derivation history.
3. **Repeated-call reduction.** Determine whether canonical-result caching, pre-canonical duplicate
   detection, or symmetry-adapted rule application materially reduces calls after the single-call
   search is corrected. This is distinct from making one canonicalization faster.
