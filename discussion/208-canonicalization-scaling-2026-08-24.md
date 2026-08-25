# 208 — Canonicalization scaling

Status: Proposed
Date: 2026-08-24
Relates: [186](186-molecule-canonicalization-2026-08-05.md),
[205](205-mapping-test-corpus-2026-08-20.md),
[207](207-reaction-network-spike-2026-08-24.md)

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

### Molecular feature level

`FeatureLevel` is a representation-owned nested level, not canonicalization-specific operational
configuration. Define it beside `Molecule` and re-export it from the graph-IR root:

```rust
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum FeatureLevel {
    Topology,
    Constitution,
    Structure,
    Full,
}
```

The derived order is the documented containment order
`Topology < Constitution < Structure < Full`. `Molecule::feature_level()` returns the lowest level
containing every populated part of the molecule: non-stereo overlay counts raise the result to
`Constitution`, stereo-atom or stereo-bond counts raise it to `Structure`, and any inline or
molecule-level constraint raises it to `Full`. The operation reads the lengths of collections that
are always present; it does not inspect feature values, validate chemistry, or canonicalize.

Replace the operation-owned `CanonicalizeLevel` with `FeatureLevel` throughout aggregate
canonicalization. The effective search level for a molecule is
`requested.min(molecule.feature_level())`. This is an exact reduction: every typed key section above
the molecule's feature level is empty. It changes neither the selected key nor the returned complete
molecule and correspondence. `IncidenceLevel` remains a distinct carrier selector because its
`Full` value contains every structural entity but no constraints.

Unary operations use the effective level above. Binary canonical equality must retain features
present on either side, so its effective level is
`requested.min(lhs.feature_level().max(rhs.feature_level()))`. This preserves the requested
distinction when only one operand contains a higher-level feature. Reaction and reaction-span APIs
adopt `FeatureLevel` as their selector in this scope, but do not acquire aggregate feature-level
inference or automatic lowering.

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
the feature-level reduction avoids imposing that unresolved case on unconstrained molecules.

## Investigation boundary

The first investigation must isolate complete canonicalization from the reaction-network loop. It
needs representative slow products, not only the ethane and propane seeds, because electronic-state
changes and disconnected intermediates may alter the residual symmetry classes. For each witness it
should record at least:

- incidence construction, normalization, refinement, backend, leaf-key, and final-remapping time;
- refinement calls, backend calls, leaves visited, and branches removed by each sound pruning rule;
- residual cell sizes after structure refinement; and
- timings and canonical results at topology, constitution, structure, and full levels.

The existing private canonical-search statistics are a suitable starting point for research
instrumentation. They are not a reason to expose a permanent diagnostics API. Slow witnesses should
be retained as canonicalization benchmarks and correctness fixtures once their behavior is
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

1. **Feature-level reduction.** Use the molecule's representation-owned `FeatureLevel` to avoid
   searching empty higher-level key sections. Verify exact canonical forms and correspondence
   action against the unreduced search on bounded dense renumberings and retained slow witnesses.
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

Topology-level canonicalization is an especially useful diagnostic for the current network domain:
its molecules contain no feature above topology, so full canonicalization has exactly the same
non-empty typed key. The `FeatureLevel` reduction makes this fact part of general graph-IR behavior
rather than a reaction-network-specific identity path. The unreduced search remains useful during
verification but is not a distinct semantic operation.

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

## Staged implementation plan

Each stage ends with a green build. Search counters and comparison switches remain private test and
benchmark infrastructure; this work does not add a public profiling API. The implementation order
separates feature reduction, orbit pruning, and prefix pruning so their correctness and performance
effects remain attributable.

### S0 — Retain witnesses and establish the search baseline

#### S0a — Retain representative canonicalization witnesses

**Module:** `umol-graph-ir/benches/canonicalize.rs` and its nearest benchmark fixture module.

Retain a small set of slow, structurally distinct products from the ethane and propane reaction
networks alongside the existing naphthalene and disconnected-ring cases. The set must include a
feature-free connected molecule, a feature-free disconnected molecule, and a symmetry-heavy
electronic-state variant. Record the molecule size, populated feature level, and provenance needed
to reproduce each value without depending on the experimental reaction-network crate.

This is additive benchmark infrastructure with no public API change.

**Tests and evidence:** Parse or construct every retained value through an existing graph-IR
boundary, assert that it is valid, and run the existing Criterion groups at all four feature levels.
The benchmark must report the retained cases individually rather than hiding them in an aggregate.

**Depends on:** none.

#### S0b — Extend private search accounting

**Module:** `umol-graph-ir/src/ir/canonicalize.rs` and
`umol-graph-ir/src/ir/canonicalize/tests.rs`.

Extend the existing private `CanonicalSearchStats` with backend-call counts and enough residual-cell
information to distinguish refinement, backend ordering, leaf comparison, orbit pruning, and prefix
pruning. Thread the counters only through internal search entry points used by module-local tests.
Criterion continues to measure timings through the public operation; do not add diagnostics to
`Canonicalize` or expose a public search result wrapper for the benchmark harness.

This is additive private instrumentation with no public API change.

**Tests and evidence:** Add example-based `rstest` cases for a singleton partition, a symmetric
partition, and an orbit-pruned partition. Assert exact counter relationships only where they are
algorithmic invariants; timings remain benchmark evidence rather than test assertions.

**Depends on:** S0a.

#### S0c — Lock the semantic baseline

**Module:** `umol-graph-ir/src/ir/canonicalize/tests.rs` and the canonicalization property-test
module.

Record the exact canonical aggregate for each retained witness and the source-to-result transport
law for its returned correspondence. Do not freeze an arbitrary correspondence representative when
multiple symmetry-equivalent representatives satisfy that law. Record baseline timings and search
counters in this document before changing dispatch or pruning.

This is additive verification work with no public API change.

**Tests and evidence:** Use example tests for retained expected aggregates and property tests for
dense entity renumberings. The properties must compare the transported aggregate and canonical key,
not the raw permutation chosen inside a symmetry class.

**Depends on:** S0b.

### S1 — Introduce the representation-owned feature level

#### S1a — Add `FeatureLevel` and molecule feature inspection

**Module:** `umol-graph-ir/src/ir/molecule.rs`, the graph-IR crate root, and their unit tests.

Add the public `FeatureLevel` enum with the settled ordering and add
`pub fn Molecule::feature_level(&self) -> FeatureLevel`. Determine the result from collection
lengths: topology alone, any non-stereo overlay, any stereo atom or bond, and any inline or
molecule-level constraint. The operation neither validates nor interprets feature values.

This is an additive public API change.

**Public surface:**

- add `FeatureLevel::{Topology, Constitution, Structure, Full}`;
- derive `Clone`, `Copy`, `Debug`, `PartialEq`, `Eq`, `PartialOrd`, and `Ord`;
- add `Molecule::feature_level()`;
- introduce no constructor, conversion failure, or new error type.

**Tests and evidence:** Use an `rstest` table covering the four levels, including both inline and
molecule-level constraints. Test the documented order and `min`/`max` behavior directly. Add a
property test stating that dense id renumbering does not change the feature level.

**Depends on:** S0c.

### S2 — Replace the canonicalization-specific selector

#### S2a — Migrate the Rust canonicalization API

**Module:** `umol-graph-ir/src/ir/canonicalize.rs`, aggregate implementations, crate exports,
benches, tests, and all Rust callers in the workspace.

Replace `CanonicalizeLevel` with `FeatureLevel` in `Canonicalize`, `Molecule`, `Reaction`, and
`ReactionSpan` canonicalization APIs. Remove `CanonicalizeLevel`; do not retain an alias or a second
selector with identical variants. `IncidenceLevel` remains unchanged. Migrate the graph-IR crate and
non-Python Rust callers here; the binding crate completes the workspace migration in S2b.

This is a breaking public rename. The graph-IR crate must be green at the end of the subitem, and
the whole stage becomes green when S2b has migrated the binding crate. There is no completed stage
with duplicate public names.

**Public surface:**

- remove `CanonicalizeLevel`;
- use `FeatureLevel` in `canonicalize_by`, `canonical_hash_by`, and `canonical_eq_by`;
- preserve the existing return values and failure behavior of every canonicalization operation.

**Tests and evidence:** Migrate existing example and property tests without weakening their laws.
Compile-check all Rust call sites and run the graph-IR test and benchmark targets.

**Depends on:** S1a.

#### S2b — Migrate the Python selector and feature query

**Module:** `umol-py/src/canonicalize.rs`, molecule bindings, package exports, type stubs or Python
surface declarations, and Python tests.

Expose the Rust `FeatureLevel` as Python `FeatureLevel`, remove Python `CanonicalizeLevel`, and add
the molecule feature-level query. Boundary conversions remain `from_rust` and `to_rust`; do not add
a second Python spelling or retain a deprecated duplicate during this experimental migration.

This is a breaking Python rename plus one additive query.

**Public surface:**

- remove Python `CanonicalizeLevel`;
- add Python `FeatureLevel` with the four Rust variants;
- add `Molecule.feature_level()` returning `FeatureLevel`;
- preserve all existing canonicalization return and exception behavior.

**Tests and evidence:** Update import, enum-conversion, molecule-query, and canonicalization tests.
Build against the repository Python 3.13 environment and run the focused Python suite.

**Depends on:** S2a.

### S3 — Lower empty feature levels exactly

#### S3a — Route molecule operations through their effective level

**Module:** `umol-graph-ir/src/ir/canonicalize.rs` and molecule canonicalization implementations.

Add private effective-level selection and use it for all molecule canonicalization entry points.
Unary canonicalization and hashing use `requested.min(molecule.feature_level())`; binary canonical
equality uses `requested.min(lhs.feature_level().max(rhs.feature_level()))`. Unqualified operations
continue to request `Full` before this reduction. Reaction and reaction-span operations use the new
selector type but retain their current search-level behavior.

This changes internal dispatch but not the public semantic result.

**Tests and evidence:** Add example tests for every requested/available level pair and asymmetric
binary cases where only one operand carries an overlay, stereo entity, or constraint. Assert that
lowered and explicitly unlowered internal searches return the same canonical aggregate and satisfy
the same transport law.

**Depends on:** S2a.

#### S3b — Verify the reduction under renumbering and retained workloads

**Module:** canonicalization property tests and `umol-graph-ir/benches/canonicalize.rs`.

Apply the existing dense-renumbering generators to feature-free and partially featured molecules.
Benchmark the retained feature-free witnesses through the ordinary `Full` API and through explicit
lower-level searches so the expected collapse in work is visible.

This is additive verification and benchmark work with no public API change.

**Tests and evidence:** Property-test equality of canonical aggregates, hashes, and transported
results under renumbering. For binary equality, generate operands with different populated feature
levels. Record search counters and Criterion results in this document; do not assert wall-clock
thresholds in tests.

**Depends on:** S3a.

### S4 — Restore sound orbit pruning at structure level

#### S4a — Preserve the backend occurrence-node action

**Module:** the graph-IR automorphism adapter and canonicalization search internals.

Retain full backend generator permutations through the graph-IR adapter instead of projecting them
immediately to source entity ids. Keep the existing source projection for canonical labels and
public correspondences. No graph-core API change is required if its existing full-node generators
are sufficient.

This is an additive internal representation change with no public API change.

**Tests and evidence:** Add adapter tests containing direct bonds, subdivided edge occurrences,
stereo-atom ligand occurrences, and stereo-bond endpoint occurrences. Assert that projection agrees
with the current source action while the retained action covers every occurrence node.

**Depends on:** S0b, S3b, and S4c. The algorithm does not require orbit pruning, but this ordering
keeps the measured effects attributable.

#### S4b — Filter generators by stereo preservation

**Module:** graph-IR canonicalization stereo-action helpers and their unit tests.

Implement the settled predicate that applies a full generator to each stereo site, ordered ligand
frame, and configuration. Retain only generators that preserve the complete stereo payload. The
feature-free fast path accepts the backend action without constructing stereo checks.

This is additive private algorithm work with no public API change.

**Tests and evidence:** Use example tests for an ordinary distinct-ligand stereocenter, a prochiral
counterexample, a global symmetry exchanging stereo sites, and a stereo bond. Include at least one
case in which a structural generator must be rejected and one in which a non-trivial generator is
retained.

**Depends on:** S4a.

#### S4c — Build safe subgroup orbits and enable pruning

**Module:** graph-IR canonicalization search.

Build source-entity orbits from the retained generator subgroup and enable orbit pruning for
`Structure` searches. Feature-free searches use the full structural orbit information. Constrained
`Full` searches remain exhaustive in this scope. If the retained-generator subgroup under-prunes a
representative stereo witness, record that fact; do not silently introduce a full stabilizer
algorithm into this subitem.

This changes private search behavior while preserving exact results.

**Tests and evidence:** Compare pruning enabled and disabled on bounded exhaustive cases, forward
and reverse candidate orders, retained stereo examples, para-stereo examples, and dense
renumberings. Assert exact canonical aggregates and transport laws, then record leaves, backend
calls, orbit-pruned branches, and timings for the retained witnesses.

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

**Depends on:** S0b and S3b.

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

**Tests and evidence:** Compare prefix pruning enabled and disabled on retained witnesses, bounded
exhaustive partitions, both candidate orders, and dense renumberings. Assert exact canonical
aggregates and correspondence transport, then record prefix-pruned branches, leaves, backend calls,
and timings. The stage is not complete merely because the counter becomes non-zero.

**Depends on:** S4c and S5b.

### S6 — Integrate, measure, and dispose of follow-on work

#### S6a — Run the graph-IR correctness and performance gate

**Module:** the complete graph-IR canonicalization surface.

Run formatting, focused unit tests, feature-gated canonicalization property tests, graph-IR linting,
and the canonicalization benchmark. Reconcile the public API inventory against the implementation
and confirm that no diagnostic type or duplicate selector escaped into the public surface.

This is verification only.

**Tests and evidence:** All graph-IR canonicalization tests and properties pass; every retained
witness preserves the exact aggregate and transport law; benchmark results and search counters are
recorded in this document.

**Depends on:** S3b, S4c, and S5c.

#### S6b — Run the Python and workspace migration gate

**Module:** `umol-py` and all workspace users of canonicalization.

Activate `umol-py/.venv`, confirm Python 3.13, rebuild the extension, run the Python tests, and run
the workspace test and lint gates. Search source, tests, examples, and documentation for the removed
`CanonicalizeLevel` name.

This verifies the breaking selector migration across language boundaries.

**Tests and evidence:** The Python and workspace gates pass with `FeatureLevel` as the only public
selector name. Any fixture or generated artifact that embeds the old spelling is either migrated or
explicitly disposed of.

**Depends on:** S2b and S6a.

#### S6c — Re-run the reaction-network witnesses

**Module:** `experimental/reaction-network` benchmark or reporting path; no production dependency
from graph-IR to the experimental crate.

Re-run the ethane and propane extended-rule closures with the same seeds, rules, and limits used for
the current evidence. Confirm identical flask, adjacency, and transformation counts before
interpreting timing differences. Reversibility need not be repeated if the network artifact is
unchanged.

This is external evidence, not a graph-IR API change.

**Tests and evidence:** Record end-to-end generation time, product-canonicalization time, and calls
per produced derivation, and pair them with the final module-local counter runs on retained graph-IR
witnesses. A larger case is justified only if the propane run no longer identifies the remaining
bottleneck.

**Depends on:** S6a.

#### S6d — Reconcile documentation and lifecycle status

**Module:** this document, `discussion/000-status.md`, and the nomenclature and development guides.

Document the implemented `FeatureLevel` contract, selector removal, effective-level rules, pruning
boundary, benchmark evidence, and any deliberately unresolved behavior. Mark this document
Completed only after the public inventory and all verification gates agree with the implementation.

This is documentation and lifecycle work.

**Tests and evidence:** Documentation examples name only current APIs, `git diff --check` passes,
and the status row describes completed rather than proposed scope.

**Depends on:** S6a, S6b, and S6c.

### Dependency summary

The critical path is
`S0a -> S0b -> S0c -> S1a -> S2a -> S3a -> S3b -> S4a -> S4b -> S4c -> S5a ->`
`S5b -> S5c -> S6a -> S6b/S6c -> S6d`.
`S2b` may proceed after `S2a` while the Rust search work continues, but it must finish before S6b.
S5 follows S4 deliberately so benchmark deltas remain attributable even though their private
infrastructure is largely independent. No implementation stage is optional within this scope.

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
