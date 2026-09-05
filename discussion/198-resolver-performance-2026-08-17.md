# 198 — Resolver performance profile

Status: In Progress
Date: 2026-08-17
Relates: [123](123-ast-allocation-survey-2026-06-21.md),
[148](148-validated-transactions-operations-2026-07-13.md),
[160](160-relevant-cycle-algs-2026-07-25.md),
[196](196-aromatic-assignment-selection-2026-08-13.md),
[199](199-open-container-integrity-2026-08-18.md)

## Motivation and outcome

The ChEBI corpus run established that the resolver handles real-world input broadly: of
201,337 nonempty SMILES records, 189,305 resolved to determined molecules, 11,478 remained
underdetermined, 528 were contradictory, 20 failed SMILES parsing, and 6 failed conversion
from the SMILES boundary representation to the graph IR. The initial baseline run took about
44 seconds in a release build on the development machine. Parsing alone took well below one
second, so parser throughput did not explain the end-to-end cost.

This document records the measurements that localize that cost and collects candidate
optimization directions. Its intended outcome is a benchmark-backed optimization effort that
preserves resolver semantics, stage-specific edit errors, algorithm selection, and the checked
`Molecule` publication boundary.

The current benchmark and prioritized follow-up appear first. The historical section at the bottom
preserves the baseline evidence, selected changes, and their gates.

## Boundaries

- No chemistry-model, resolution-policy, or `Solution` outcome changes are proposed. A public
  execution-error variant may be removed only when the revised lifecycle makes that failure
  unreachable.
- No silent algorithm default or runtime dependency on a comparison toolkit is proposed.
- The RDKit measurements provide scale and a separate parse-only reference; behavioral or
  performance parity is not a requirement.
- The multicore measurements use independent corpus shards in separate processes. Neither
  library is internally parallel on the measured path.
- ChEBI remains a local external corpus. A checked-in regression suite should use small,
  redistributable cases promoted from corpus findings.
- Parser optimization remains out of scope because it is not a leading cost. Valence optimization
  was outside the initial pass; the current profile brings candidate generation into follow-up scope.

## Current benchmark snapshot

The current measurements cover the optimization build through the checked immutable-application
failure correction on 2026-08-18. They use the workspace lockfile, Cargo release profile, ChEBI
release 254, and the Apple M1 Pro host described in the historical baseline below. The comparison
environment uses RDKit 2026.3.4. Inputs were loaded before the timed loops.

Three umol passes retained the exact corpus digest
`550d420e000f3d9df17c339e22a08063` and counts
`189305,11478,528,20,6,0`. Their median full-corpus outer time is 20.231 seconds, compared with
21.713 seconds for the recorded one-worker RDKit full-ingest run. The local difference is 1.482
seconds, or 6.8%. The paths perform different work, so this places umol below the recorded RDKit
median for these specific local workloads without establishing a portable or stable lead.

## Current umol stage decomposition

The median columns from the three current release passes decompose as follows. The remainder is the
difference between the median outer time and the displayed stage medians.

| Stage | Time (s) | Share of outer time |
| ----- | -------: | ------------------: |
| SMILES parse | 0.324 | 1.6% |
| Boundary representation to graph IR | 2.058 | 10.2% |
| Resolver construction/setup | 0.014 | less than 0.1% |
| `Resolver::resolve` | 17.565 | 86.8% |
| Loop, drops, and unmeasured remainder | about 0.269 | about 1.3% |
| Total | 20.231 | 100.0% |

The resolver remains the dominant stage. The front end, including graph-IR construction, now takes
about 2.38 seconds and is not the immediate optimization target.

## Current sampling profile

The Samply run of the linear consuming merge took 20.963 seconds and collected 18,166 samples under
`Resolver::resolve`. Its main inclusive regions and direct resolver children were:

No fresh profile was collected for the checked immutable-application correction because the corpus
workload does not call `Molecule::apply`; the resolver profile below therefore remains the current
evidence for the exercised path.

| Region | Share of resolver samples | Interpretation |
| ------ | ------------------------: | -------------- |
| Aromaticity selection | 33.7% | Largest remaining resolver phase |
| Ring enumeration | 19.6% | Nested within aromaticity selection |
| Checked publication | 15.6% | `build`, `try_build`, and integrity checking |
| Valence admission | 14.7% | Mostly counts-valence candidate generation |
| Checked edit application | 10.5% | Direct resolver child, separate from publication |
| Discharge planning | 6.7% | Closing assertion evaluation and edit planning |
| Perception from reused rings | 4.6% | Nested within aromaticity selection |
| Localized-bond planning | 3.0% | Direct late-phase planning cost |
| Stereo planning | 2.8% | Direct late-phase planning cost |

The rows overlap where noted and must not be summed. The largest named exclusive Rust leaves were
`RelevantCycleAnalysis::new` at 4.1%, `HashMap::insert` at 3.9%,
`BuildHasher::hash_one` at 3.2%, and `AtomConstraintsForm::find` at 2.6%; allocator leaves were also
prominent. The corresponding inclusive stacks place `RelevantCycleAnalysis::new` at 13.8%,
`CountsValence::candidate_states` at 10.1%, and `Molecule::check_integrity` at 14.0%. The former
`merge_overlapping_systems` hotspot no longer appears among the leading named leaves. Allocation,
freeing, hash-table growth, vector collection, and copy-on-write materialization remain distributed
across the remaining operations instead of forming one independent phase.

Ring enumeration remains expensive, but it no longer explains most aromaticity-selection time.
Checked publication and checked application are both large; changing either requires phase-level
attribution and preservation of the publication boundary. Valence candidate generation is now
comparable to each lifecycle region and is no longer negligible. The remaining aromaticity cost is
distributed across candidate formation, ring enumeration, hashing, allocation, and selection
instead of being dominated by overlap merging.

## Prioritized follow-up

1. Attribute checked publication and edit application to individual resolver boundaries. Separate
   integrity walks, editor creation, copy-on-write materialization, and edit execution before
   deciding whether another phase boundary can be removed or a data path should change.
2. Profile valence candidate generation in isolation. Determine whether the 10.1% candidate-state
   stack is dominated by repeated constraint lookup, `NumForm` cloning and normalization, candidate
   collection, or genuinely necessary state enumeration before selecting an optimization.
3. Revisit relevant-cycle construction. Its 13.8% inclusive and 4.1%
   exclusive shares still justify investigation, but the profile exposes no second redundant whole
   pass; any change remains subject to doc 160's conformance and comparative gates.
4. Examine discharge only after the larger regions. Its 6.7% share may contain repeated derived-view
   construction and constraint traversal, but it is not yet large enough to justify weakening the
   explicit closing validation step.
5. Promote small molecules that reproduce the valence, publication/application, and relevant-cycle
   costs into checked-in Criterion cases. Retain the pathological aromatic candidate cases as a
   targeted regression. ChEBI remains the acceptance benchmark; the small cases provide iteration
   speed and regression localization.

This is a profile-derived investigation order, not an S0/S1 implementation plan. Each item should be
remeasured after the preceding change because the relative shares will move. The questions that
remain are expressed by these investigation targets; the earlier open questions have either been
answered by the completed work or incorporated here.

## Verification and benchmark gates

Each optimization should preserve:

- the determined, underdetermined, contradictory, parse-error, and conversion-error corpus counts;
- the exact resolved molecules for determined inputs, not only their classification;
- existing resolver fixtures, semantic property tests, and ring-algorithm conformance tests;
- configured algorithm selection and error behavior;
- source-preserving failure behavior, transaction rollback guarantees, and representation-integrity
  guarantees.

Performance evidence should use a release build, the same lockfile and corpus snapshot, input loaded
outside the timed region, at least three passes, and both wall timing and a fresh profile. The first
four-worker scaling result is already adequate; optimization decisions should use single-process
timings so scheduler effects do not obscure library costs. No change should be accepted solely
because an inclusive sampling percentage falls.

## Historical benchmark and optimization record

The remainder of this document preserves the baseline measurements, initial findings, and sequence
of implemented changes. These values describe earlier builds and are not the current performance
summary. They remain useful for attribution and for reconstructing the cumulative improvement.

### Baseline benchmark snapshot

The measurements were taken at commit `9cfa0228a3c6` with the workspace lockfile, the Cargo
release profile, and ChEBI release 254. The host was an Apple M1 Pro with eight performance and
two efficiency cores, running macOS 15.7.3. The comparison environment used RDKit 2026.3.4.
The parse-only and stage timers used input already loaded in memory. The sampling profile included
input loading, which is identified separately below.

Timings below are local observations, not portable performance guarantees. Small differences
should be remeasured on the same host and input before being interpreted.

### Baseline process-level scaling

The full-ingest workloads were split into equal corpus shards and executed in independent
processes. `busy cores` is aggregate process CPU time divided by wall time; occupancy divides
that value by the worker count. Each curve measures its library's ordinary full processing path,
not equal work across the two libraries; the purpose of this table is scaling within each curve.

| Workers | RDKit wall (s) | RDKit speedup | RDKit occupancy | umol wall (s) | umol speedup | umol occupancy |
| ------- | -------------- | ------------- | --------------- | ------------- | ------------ | -------------- |
| 1       | 21.713         | 1.00x         | —               | 45.612        | 1.00x        | —              |
| 2       | 11.083         | 1.96x         | 97.1%           | 22.298        | 2.05x        | 98.4%          |
| 4       | 5.702          | 3.81x         | 97.0%           | 11.797        | 3.87x        | 97.9%          |
| 8       | 3.508          | 6.19x         | 87.5%           | 6.506         | 7.01x        | 93.6%          |

The robust result is near-linear scaling from one through four processes. Eight processes still
help substantially, but the comparative efficiency is not stable enough to explain as a library
property. An earlier three-repeat curve gave RDKit 7.18x and umol 6.74x at eight workers; the
diagnostic run above reverses that order. Eight workers occupy all eight performance cores and
leave little scheduling, frequency, cache, thermal, or straggler headroom. Because workers share
no library state, these results provide no evidence of a resolver lock or other shared internal
serialization.

### Baseline parse-only comparison

Both parse-only loops used the same in-memory strings and the median of three passes. umol called
the Rust `Smiles::parse_bytes` path. RDKit called its Python SMILES API with `sanitize=False` and
`removeHs=False`.

| Parser | Median (s) | Records/s | Accepted | Rejected |
| ------ | ---------- | --------- | -------- | -------- |
| umol   | 0.265      | 760,813   | 201,317  | 20       |
| RDKit  | 3.171      | 63,501    | 201,337  | 0        |

On this workload the direct umol parser is about 12.0x faster. This is not a conformance or
equal-work claim: the accepted languages and result representations differ, and the RDKit call
crosses the Python boundary. The scale of the difference does establish that SMILES tokenization
and syntax construction are not the current umol bottleneck. Adding conversion from the SMILES
boundary representation to graph IR brings the umol front end to about 2.79 seconds, close to but
still below the unsanitized RDKit parse time.

### Baseline umol stage decomposition

Three release-mode corpus passes gave a median outer time of 44.18 seconds. A dedicated outer
parse loop replaces the inflated per-record `Instant` measurement for the parsing row.

| Stage | Time (s) | Share of outer time |
| ----- | -------- | ------------------- |
| SMILES parse | 0.265 | 0.6% |
| Boundary representation to graph IR | 2.52 | 5.7% |
| Resolver construction/setup | 0.017 | less than 0.1% |
| `Resolver::resolve` | 40.86 | 92.5% |
| Loop, drops, and unmeasured remainder | about 0.52 | about 1.2% |

Direct timers around the first resolver phases further divided the median resolver time:

| Resolver region | Time (s) | Share of resolver time |
| --------------- | -------- | ---------------------- |
| Valence admission | 2.59 | 6.3% |
| Aromaticity selection | 18.81 | 46.0% |
| Later phases and resolver framework | 19.46 | 47.6% |

Valence resolution is therefore not predominant. Aromaticity is the largest individually timed
phase, while approximately half the resolver time remains in later phases and common lifecycle
work.

### Baseline sampling profile

A symbolized sampling run collected 25,799 samples during one corpus pass. The profiled outer
time was 47.22 seconds; input loading accounted for about 8.5% of all samples and
`Resolver::resolve` for 84.84%. The following values normalize inclusive stack samples to the
resolver root.

| Inclusive stack | Share within `Resolver::resolve` |
| --------------- | ------------------------------- |
| `AromaticityResolver::select` | 45.0% |
| `RingSet::enumerate` | 33.0% |
| Relevant-cycle enumeration/visitation | 30.5% |
| Vismara `RelevantCycleAnalysis::new` | 23.8% |
| `Molecule::check_integrity` | 16.7% |
| `MoleculeEditor::{try_build,build}` | 15.5% |
| `MoleculeEditor::transact` | 10.8% |
| `Transaction::append` | 8.8% |
| Valence admission | 6.1% |
| Discharge | 2.8% |
| Stereo planning | 1.4% |
| Bond planning | 1.3% |

These are inclusive and overlapping percentages. In particular, Vismara analysis is nested under
ring enumeration, and integrity checking is nested under molecule publication; the rows must not
be summed.

Large exclusive leaves included memory moves, allocation and freeing, hash construction and
insertion, Vismara analysis, breadth-first search, and aromatic-system merging. The profile
therefore identifies an allocation/hash/memory-traffic cluster in addition to the graph algorithms.
It does not by itself prove that any particular clone or allocation is unnecessary.

### Initial findings

#### Candidate rings are recomputed for every complete aromatic assignment

`AromaticityResolver::select` computes a `RingSet` once to form candidate components and claim
candidates. At each complete assignment it then calls `AromaticityPerception::find_systems`.
`find_systems` unconditionally calls `candidate_rings`, which repeats ring enumeration before it
dispatches to the rule's existing `find_from_rings` implementation.

The molecule topology and perception configuration do not change while these assignments are
tested; only the electron-contribution closure changes. Reusing the already computed `RingSet`
should therefore remove repeated Vismara work without changing which systems are evaluated. This
is the clearest first optimization candidate because source inspection explains the largest
profile stack directly.

The API shape is not settled. A crate-internal operation that accepts an existing `RingSet` is
preferable to widening the public surface solely for this call site.

#### Checked intermediate publication repeats whole-molecule integrity work

`MoleculeEditor::try_build` constructs a `Molecule` and runs the full representation-integrity
check. `build` delegates to that checked path. `Resolver::resolve` publishes and then re-edits
several intermediate molecules across placement, constitution, stereo, ordinary bonds,
multicenter bonds, and discharge. Determined inputs therefore pay for several complete integrity
walks.

One integrity pass over all 201,311 successfully raised corpus molecules took about 0.77 seconds.
Repeated publication amplifies that cost; the sampling profile attributes 16.7% of resolver stack
samples to integrity checking. A candidate design should retain checked public and final
publication while avoiding redundant full checks for internal phase states whose provenance
establishes the same invariants. The lifecycle and rollback contract must be reconciled with doc
148 before an unchecked or phase-specific internal path is introduced.

#### Integrity checks allocate for tiny fixed participant sets

The common `check_unique_participants` helper creates a fresh `HashSet` for every entity. It is
called for every ordinary bond even though an ordinary bond always has exactly two participants.
A direct inequality check suffices for that case. Other small fixed-arity families may admit
pairwise comparison or small local storage, leaving a `HashSet` only for genuinely variable-arity
relations.

This will not eliminate the repeated integrity passes, but it targets the allocation and hashing
leaves within them without weakening the invariant.

#### Phase publication, copy-on-write materialization, and journals form a second cost cluster

`Molecule::edit` clones the graph and constraints and clones the `Arc` handles for entity storage.
The resolver repeatedly cycles through `build`, `edit`, `transact`, and `Transaction::append` while
retaining immutable snapshots for planning. Subsequent mutation may materialize shared storage,
and the combined rollback journal grows by extending its undo vector. The profile places these
operations beside substantial allocation, move, and free activity.

The profile establishes the cluster, not the avoidability of each operation. Follow-up measurement
should distinguish necessary rollback state from incidental copying and reallocation. Candidate
changes include fewer phase publications, a consuming editor transition after a snapshot is no
longer needed, and journal capacity management.

#### Vismara optimization should follow removal of redundant invocations

Vismara relevant-cycle analysis is individually expensive, but the resolver currently invokes it
more often than the semantics require. Reusing the `RingSet` may change both its absolute cost and
the corpus distribution of remaining hot cases. Algorithm-level changes belong only after that
reuse is implemented and the workload is reprofiled. Any later alternative remains subject to doc
160's correctness and comparative-evaluation requirements.

### Implementation record

#### Candidate-ring reuse (2026-08-17)

`AromaticityPerceiver::find_systems_from_rings` now performs perception against an existing
crate-internal `RingSet`. The public `find_systems` operation retains its ordinary call shape and
constructs that ring set before delegating. `AromaticityResolver::select` passes the `RingSet` it
already uses for component formation and claim candidates through every complete assignment, so
assignment evaluation no longer repeats relevant-cycle enumeration.

The definition-level property reference deliberately continues to call public `find_systems`.
The optimized resolver path therefore remains cross-checked against a path that constructs rings
independently. The focused aromaticity suites and the resolver search-independence property pass.

The same seven-case `--quick` Criterion run before the change, after ring reuse, and after the
fixed-pair integrity change below gave these median observations:

| Case | Baseline (us) | Ring reuse (us) | Ring reuse + fixed pairs (us) | Combined change |
| ---- | ------------: | --------------: | ----------------------------: | --------------: |
| methane | 4.141 | 4.100 | 4.154 | +0.3% |
| octane | 23.351 | 22.881 | 20.806 | -10.9% |
| benzene | 38.448 | 31.377 | 28.540 | -25.8% |
| pyridine | 39.397 | 31.222 | 29.428 | -25.3% |
| naphthalene | 74.332 | 59.081 | 57.555 | -22.6% |
| quinoline | 73.090 | 60.109 | 60.468 | -17.3% |
| purine | 128.480 | 81.410 | 82.190 | -36.0% |

These short runs establish the scale and localization of the change, not a portable performance
claim. A subsequent 100-sample run of the combined implementation measured 28.67 us for benzene,
29.36 us for pyridine, 56.35 us for naphthalene, 56.66 us for quinoline, and 73.63 us for purine.

The full release-mode ChEBI runner was then rebuilt through the workspace lockfile and applied to
all 201,337 nonempty release-254 SMILES records:

| Pass | Outer (s) | Parse (s) | Raise (s) | Setup (s) | Resolve (s) |
| ---: | --------: | --------: | --------: | --------: | ----------: |
| 1 | 33.864 | 0.359 | 2.137 | 0.022 | 30.910 |
| 2 | 32.754 | 0.353 | 2.144 | 0.019 | 29.838 |
| 3 | 30.330 | 0.342 | 2.120 | 0.016 | 27.531 |
| Median | 32.754 | 0.353 | 2.137 | 0.019 | 29.838 |

The median outer time fell from 44.18 to 32.75 seconds (25.9%), and the median resolver time fell
from 40.86 to 29.84 seconds (27.0%). Every pass retained the baseline classification counts:
189,305 determined, 11,478 underdetermined, 528 contradictory, 20 parse failures, 6 graph-IR
conversion failures, and no resolver execution errors. The baseline runner did not retain resolved
molecule outputs, so exact molecule-by-molecule equivalence remains outstanding, as does a fresh
sampling profile. The local external-corpus harnesses are retained under the ignored
`materials/benchmarks/` work area for the remainder of this optimization effort.

#### Fixed-pair integrity checks (2026-08-17)

Ordinary and noncovalent bonds now detect a repeated endpoint by direct equality instead of
allocating and populating a fresh `HashSet`. Their exact `DuplicateParticipant` diagnostics are
unchanged. Existing integrity cases cover both self-loop families. The incremental quick-run
effect was largest for ordinary localized structures (about 9% for octane, 5-9% for the single-ring
cases) and indistinguishable from noise for the larger fused cases. Variable-arity and stereo
participant checks still use the general helper; changing those requires size-distribution evidence.

The repeated-publication candidate has not been implemented. `Resolver` and `MoleculeEditor` live
in different crates, so a local unchecked-build shortcut would widen the graph-IR API and weaken
the publication contract. Work on that item requires the phase-lifecycle design described below.

#### Aromaticity activation gate (corrected 2026-08-18)

`AromaticityResolver::select` returns the valence carrier unchanged before constructing candidate
rings unless the molecule contains an atom `Aromatic(_)` assertion, a localized- or dative-bond
`Lit(true)` aromatic assertion, or an existing aromatic-system entity. `#a!` excludes aromaticity
and `#a*` is vacuous; neither requires selection. The same applies to false and undetermined bond
forms. The initial presence-only implementation was overbroad because the raised SMILES
representation writes an atom aromaticity constraint even for non-aromatic atoms.

A corpus diagnostic over the post-valence carrier found that 189,698 ChEBI molecules reach this
stage, of which 111,303 require aromatic processing and 78,395 take the corrected early exit. The
exact-output corpus gate retained `550d420e000f3d9df17c339e22a08063` and counts
`189305,11478,528,20,6,0`.

Three release-mode passes measured:

| Pass | Outer (s) | Parse (s) | Raise (s) | Setup (s) | Resolve (s) |
| ---: | --------: | --------: | --------: | --------: | ----------: |
| 1 | 21.429 | 0.328 | 2.051 | 0.016 | 18.776 |
| 2 | 21.456 | 0.326 | 2.051 | 0.014 | 18.801 |
| 3 | 21.368 | 0.323 | 2.044 | 0.014 | 18.726 |
| Median | 21.429 | 0.326 | 2.051 | 0.014 | 18.776 |

Relative to the consuming-application build, median resolver time fell by 11.5% and median outer
time by 10.3%. Relative to the original snapshot, the cumulative median reductions are 54.0% for
resolver time and 51.5% for outer time.

A fresh sampling profile after this change leaves aromaticity selection at 38.7% of resolver
samples and ring enumeration at 18.1%. Checked publication is 14.1%, valence admission 13.3%, edit
application 9.9%, and discharge planning 6.3%. `merge_overlapping_systems` is the largest named
exclusive Rust leaf at 4.5%, followed by `RelevantCycleAnalysis::new` at 3.7%. These inclusive
categories overlap and must not be summed.

#### Linear consuming aromatic-system merge (2026-08-18)

The previous merger compared every accepted aromatic candidate with every later candidate using
`HashSet::is_disjoint`, then rebuilt each connected component through a `HashMap` of candidate
indices. Corpus instrumentation found 204,883 merger calls. Only 35 calls had more than 50
candidates, but those calls accounted for 27.17 million of the 28.68 million pair comparisons. The
largest call had 4,656 candidates and performed 10.84 million comparisons to produce one 30-atom
system. The quadratic tail therefore came from a very small number of chemically difficult inputs.

The replacement records the first candidate containing each atom in a dense owner vector and
unions a later candidate with that owner. A second pass consumes each candidate `HashSet` directly
into its union-find component. Its work is proportional to the candidate memberships apart from
near-constant union-find operations; it performs no candidate-pair scan and no grouping hash map.
Individual and fused-ring candidate generation is unchanged.

A temporary private selector compared the old merger, this linear consuming merger, and a direct
union implementation within the same release build. Rotated full-corpus passes measured:

| Merger | Median outer (s) | Median resolve (s) |
| ------ | ---------------: | -----------------: |
| Pairwise quadratic | 21.656 | 18.976 |
| Linear consuming | 20.267 | 17.604 |
| Direct union | 20.302 | 17.609 |

All variants retained digest `550d420e000f3d9df17c339e22a08063` and counts
`189305,11478,528,20,6,0`. On a 17-record high-candidate subset, ten repeats gave median resolver
times of 1.461 seconds for the quadratic merger, 0.146 seconds for the linear merger, and 0.140
seconds for direct union. Direct union's small advantage on that subset did not carry into the full
corpus and required restructuring perception to union accepted candidates during generation. The
linear consuming implementation was therefore retained and the temporary selector and the other
implementations were removed.

Three clean passes of the final implementation measured:

| Pass | Outer (s) | Parse (s) | Raise (s) | Setup (s) | Resolve (s) |
| ---: | --------: | --------: | --------: | --------: | ----------: |
| 1 | 20.092 | 0.332 | 2.083 | 0.016 | 17.401 |
| 2 | 20.305 | 0.333 | 2.081 | 0.014 | 17.608 |
| 3 | 20.181 | 0.328 | 2.066 | 0.014 | 17.507 |
| Median | 20.181 | 0.332 | 2.081 | 0.014 | 17.507 |

Relative to the corrected activation-gate build, median resolver time fell by 6.8% and median outer
time by 5.8%. Relative to the original snapshot, the cumulative median reductions are 57.2% for
resolver time and 54.3% for outer time. The post-linear profile summarized at the beginning of this
document no longer contains the merger as a leading named leaf. The six direct component-merging
cases, all 975 `umol-graph` library tests, the 256-case resolver search-independence property, all
675 resolution-conformance cases, and strict all-target clippy pass.

#### Exact-output corpus gate (2026-08-17)

The local ChEBI harness now computes an order-sensitive XXH3-128 digest over every nonempty corpus
record. Each frame includes the record index, outcome class, payload length, and payload. Determined
outcomes carry the full positional-EDN rendering of the resolved molecule; failures and unresolved
outcomes carry their exact class. The gate therefore detects molecule changes that aggregate outcome
counts cannot.

Two fresh pre-change runs produced the same baseline:

```text
digest=550d420e000f3d9df17c339e22a08063
counts=189305,11478,528,20,6,0
```

The transaction-chain, phase-publication, Vismara, and consuming-application changes below and the
activation-gate and linear-merger changes above all reproduce that exact digest and those counts.
This closes the exact molecule-by-molecule comparison left outstanding by the earlier count-only
corpus runs.

#### Resolver transaction chain (2026-08-17)

The resolver formerly flattened each successful phase transaction into one growing rollback journal
with `Transaction::append`. `append` consumed the later journal, so it did not deep-clone undo
payloads, but every inline `Undo` value still moved into the combined allocation and capacity growth
could move the accumulated prefix again. Sampling attributed 12.0% of the preceding resolver profile
to this inclusive stack beside substantial `memmove`, allocation, and free activity.

The resolver now retains the six operation-issued transactions separately and rolls them back in
reverse transaction order. This produces the same global undo order without copying undo entries or
changing the public `Transaction` contract. The focused late-underdetermination rollback test and the
exact corpus gate pass.

Fresh three-pass measurements immediately before and after the change were:

| Build | Median outer (s) | Median resolve (s) |
| ----- | ---------------: | -----------------: |
| Flattened journal | 29.871 | 27.154 |
| Transaction chain | 28.457 | 25.755 |
| Change | -4.7% | -5.2% |

The fresh profile contains no `Transaction::append` samples. The combined publication, transaction,
and journal stack fell from 40.0% to 30.1% of resolver samples. Dropping the retained transaction
chain accounts for 3.0%, so successful-journal storage remains a measurable cost.

`Undo` is 624 bytes on the measured target, while `Edit` is 256 bytes and `Transaction` is a 24-byte
`Vec` header. The wide enum explains why unnecessary bulk movement was expensive and may affect
journal construction and destruction more generally. Before changing its public representation,
follow-up work should measure corpus entry counts by variant and allocated bytes. Selectively boxing
rare reconstruction-heavy variants is preferable to boxing every common field undo if the measured
distribution supports it.

#### Shared late-phase publication (2026-08-17)

Stereo planning, localized-bond defaults, and multicenter-bond defaults read independent domains of
the same post-constitution molecule. Their transactions write disjoint entity fields or overlays.
The resolver now plans all three from that one immutable snapshot, then applies their checked
transactions in the original order. This removes the checked `build`/`edit` boundaries after stereo
and localized bonds while preserving error precedence and rollback classification. No unchecked or
additional public graph-IR operation was introduced.

The exact corpus gate and all 969 `umol-graph` library tests pass. Relative to the transaction-chain
build, the median outer time fell from 28.457 to 26.992 seconds (5.1%) and median resolver time fell
from 25.755 to 24.276 seconds (5.7%). The fresh profile reduced checked publication from 15.7% to
11.1% and the combined lifecycle stack from 30.1% to 26.3% of resolver samples.

#### Vismara ordered-distance pass (2026-08-17)

For each root in a biconnected component, `ShortestPathDag::vismara` formerly ran one BFS over the
component and a second BFS over the degree-ordered preceding subgraph. It also reconstructed the
component-membership and degree scratch arrays for every root or component.

The second BFS is equivalent to a dynamic pass over the first BFS visitation order: a permitted
vertex retains its full-component distance exactly when an already included neighbour occurs at the
preceding distance. The implementation now uses that pass and reuses component membership and degree
buffers. Cycle-family formation, selection, visitation, and output order are unchanged.

The 19 focused Vismara unit tests and 86 cycle property, literature, and captured-corpus tests pass;
the latter include the externally captured normalized relevant-cycle corpus. The exact ChEBI digest
also passes. Relative to the shared-phase build, median outer time fell from 26.992 to 26.621 seconds
(1.4%) and median resolver time fell from 24.276 to 23.925 seconds (1.4%). The new profile leaves ring
enumeration at 21.9% and `RelevantCycleAnalysis::new` at 15.2% of resolver samples, so ring internals
remain important but no longer have an equally obvious redundant whole pass.

Across the transaction-chain, shared-phase, and Vismara changes, median resolver time fell from
27.154 to 23.925 seconds (11.9%) and median outer time from 29.871 to 26.621 seconds (10.9%). Relative
to the original snapshot in this document, the cumulative median reductions are 41.4% for resolver
time and 39.7% for outer time.

#### Consuming edit application (2026-08-17)

Resolver phases formerly applied edits with `MoleculeEditor::transact`, retained every successful
undo journal, and replayed those journals on every non-determined exit. That rollback was not needed
for caller-visible atomicity: resolution operates on an editor cloned from the input molecule and
does not replace the input until a determined result has passed the final checked publication gate.

`MoleculeEditor::apply(self, Edits) -> Result<Self, TransactionError>` now provides checked,
non-journaled application. It consumes the editor, returns the modified editor on success, and drops
the inaccessible partial draft on failure. `transact` retains its mutable-borrowed, atomic,
rollback-journaled contract for callers that need to preserve and continue using an editor.
`transact_unchecked` has been removed. The immutable `Molecule::apply` convenience now uses the
consuming path before building its result, so it no longer creates and discards an undo journal.

`Resolver::resolve` now keeps the supplied molecule unchanged on underdetermination, contradiction,
or execution error by discarding its private draft. The rollback-only `ResolveRollbackCause` type
and `ResolveError::RollbackFailed` variant have consequently been removed. Placement,
constitution, and discharge retain their checked molecule publication boundaries; this change does
not weaken representation-integrity checking.

The exact ChEBI corpus gate retained
`550d420e000f3d9df17c339e22a08063` and counts
`189305,11478,528,20,6,0`. Three release-mode passes measured:

| Pass | Outer (s) | Parse (s) | Raise (s) | Setup (s) | Resolve (s) |
| ---: | --------: | --------: | --------: | --------: | ----------: |
| 1 | 23.950 | 0.318 | 2.057 | 0.018 | 21.285 |
| 2 | 23.811 | 0.317 | 2.042 | 0.014 | 21.165 |
| 3 | 23.881 | 0.317 | 2.048 | 0.014 | 21.226 |
| Median | 23.881 | 0.317 | 2.048 | 0.014 | 21.226 |

Relative to the preceding Vismara build, median resolver time fell by 11.3% and median outer time
by 10.3%. Relative to the original snapshot, the cumulative median reductions are 48.1% for
resolver time and 45.9% for outer time. Focused graph-IR application tests, the generated
`apply`/`transact` result-equivalence property, the graph-IR library suite (6,064 passed, 3 ignored),
all 966 graph library tests, and strict clippy pass.

#### Checked immutable-application failure (2026-08-18)

`Molecule::apply` formerly used asserted `MoleculeEditor::build` after a successful edit batch.
Caller-built edits could therefore pass transaction checks but panic when final integrity checking
rejected the result. The operation now uses `try_build` and returns a `MoleculeApplyError` that
distinguishes transaction failure from molecule-integrity publication failure. This does not alter
the successful publication work or `Resolver::resolve`, which uses `MoleculeEditor::apply` directly.

The exact ChEBI corpus gate retained
`550d420e000f3d9df17c339e22a08063` and counts
`189305,11478,528,20,6,0`. Three release-mode passes measured:

| Pass | Outer (s) | Parse (s) | Raise (s) | Setup (s) | Resolve (s) |
| ---: | --------: | --------: | --------: | --------: | ----------: |
| 1 | 20.314 | 0.331 | 2.068 | 0.016 | 17.635 |
| 2 | 20.231 | 0.324 | 2.058 | 0.014 | 17.565 |
| 3 | 19.954 | 0.317 | 2.036 | 0.014 | 17.325 |
| Median | 20.231 | 0.324 | 2.058 | 0.014 | 17.565 |

Relative to the preceding linear-merger measurements, median outer time changed from 20.181 to
20.231 seconds (0.2%) and median resolver time from 17.507 to 17.565 seconds (0.3%). These small
differences are ordinary run variation, and the benchmark does not execute the changed
`Molecule::apply` path.

### Initial optimization order

1. Reuse the precomputed candidate `RingSet` throughout aromatic assignment evaluation.
   Implemented 2026-08-17; the full corpus classification, exact-output, timing, and fresh-profile
   gates pass.
2. Skip aromaticity selection when the input has no positive aromatic atom or bond assertion and no
   stored aromatic system. Corrected 2026-08-18 so `#a!` and `#a*` take the fast path; the exact
   corpus gate and release timing gate pass.
3. Reduce redundant checked resolver-phase publication without weakening the published `Molecule`
   contract or rollback behavior. Implemented for the independent stereo, localized-bond, and
   multicenter-bond phases; earlier phase boundaries remain dependency-bearing.
4. Make participant uniqueness checks allocation-free for fixed and small arities. Implemented for
   the two fixed-pair families; other small arities remain evidence-dependent.
5. Measure and reduce avoidable phase cloning, copy-on-write materialization, and rollback-journal
   work. Resolver journals are eliminated through consuming checked application; remaining work is
   limited to phase snapshots, copy-on-write materialization, and repeated integrity publication.
6. Reprofile; optimize or replace parts of Vismara analysis only if it remains a leading exclusive
   cost. The redundant ordered-distance BFS is removed; further work requires a more specific
   benchmark-backed target.

This was the working order used for the first optimization pass. The current profile and priorities
at the beginning of this document replace it for subsequent work.
