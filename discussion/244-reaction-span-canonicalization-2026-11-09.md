# 244 — Reaction-span canonicalization aligned with molecules

Status: In Progress
Date: 2026-11-09
Relates: [186](186-molecule-canonicalization-2026-08-05.md),
[208](208-canonicalization-scaling-2026-08-24.md),
[237](237-release-0-10-scope-2026-10-15.md)

## Purpose

`ReactionSpan` has the same eight entity kinds and typed incidences as `Molecule`, with each entity
value lifted to an `EntitySpan`. Molecule canonicalization selects a private description level and
runs the search defined for that level. Reaction-span canonicalization selects the level by the same
rule and then runs one exhaustive search at every level. `Reaction` canonicalizes through its span.

This document scopes the work that gives reaction spans the molecule search of each description
level, with every technique of that search that applies to a span. It also tries orbit pruning at
Full for both aggregates, under the guide's condition for constraints, which doc 208 left
unresolved. Upstream PR 29 raised both points. It is reviewed below and is not
merged; the work is implemented separately.

The public `Canonicalize` surface does not change, and `canonical_eq` relates the same pairs.

## Review of PR 29

PR 29 (`perf/canonical-orbit-pruning`) branches from the 0.9.0 merge; `main` is 367 commits ahead.
It makes four changes:

1. reaction-span search prunes orbits below Full, with a span copy of the stereo generator filter;
2. reaction-span search prunes at Full when every constraint is an atom constraint other than #T;
3. molecule search prunes at Full under the same condition;
4. benchmarks, tests, and a discussion document.

Findings, with the review points of doc 237 named:

- **Base.** Canonicalization on `main` has since moved to `tracked_canonicalize_at` and
  `canonical_key_at` on each aggregate, with per-level functions that take the owned molecule. Two
  of the functions the PR patches no longer exist. Its tests build a #T constraint with a coset,
  which no longer exists. Its document number collides with 231.
- **Search shape.** The PR's span search keeps the constitution-level initial partition at Structure
  and Full, node-id branch order at every level, the complete adapter at Topology, and a remap and
  reframe of the whole span at every leaf. It does not read `para_stereo`. These are the points
  where it differs from the molecule search.
- **Selected value and remapping (237).** The PR keeps node-id order so that pruning selects the
  leaf the exhaustive search selects: same entity order, canonical value, and remapping. The
  argument holds when every retained generator fixes the cells of the current partition and
  preserves the leaf key; the backend request colors by cell, and key preservation is what the
  generator filters check. The guide requires that pruning and branch order leave the
  representative unchanged, which the minimum key gives. It does not fix the remapping within a
  symmetry class.
- **Generator filters (237).** The span stereo filter is a copy of the molecule filter that checks
  both sides under one frame action. The Full-level filter compares atom constraints only, so it
  runs behind an eligibility test over the constraints of all eight entity kinds. That test exists
  twice, once per aggregate. A filter that compares the constraints of every entity kind needs no
  eligibility test.
- **Coverage beyond identity spans (237).** Every span in the PR's tests has `Unchanged` entities
  only; one case adds a molecule-level constraint span. Its benchmarks use identity spans. No case
  has an added, removed, or modified entity, or a stereo entity whose sides differ, so the two-sided
  parts of the filters are not exercised.
- **Regressions on stereo cases (237).** The PR reports +6% to +13% on three stereo cases. Its span
  search starts from the constitution-level partition. Atoms that only stereo distinguishes stay in
  one cell and become branch points, and each branch point makes a backend request that prunes
  nothing. The molecule search refines by stereo descriptors before it branches.
- **#T and bond constraints (237).** The PR excludes #T because its value was a coset relative to a
  ligand frame, and treats any bond constraint, including #C, as ineligible. On `main`, #T and #C
  hold three frame-independent states, so neither exclusion corresponds to stored data.

## Settled basis

The [data-type guide](../docs/development/data-types.md) fixes the following for aggregate
canonicalization of `Molecule`, `ReactionSpan`, and `Reaction`:

- The typed order defines the canonical frame. Pruning, branch order, and the automorphism algorithm
  are operational and must not change the representative.
- Stored stereo entities take part through one stereo-sensitive refinement from the
  constitution-level partition. With para-stereo refinement enabled, the refinement is repeated
  until the ordered partition is stable.
- Constraints do not distinguish structural orbits. Complete canonicalization minimizes the pair of
  structure key and constraint key. Orbits may discard a tied frame only when the orbit action
  carries the complete participant-frame and constraint action; otherwise tied frames are searched
  exactly.
- Reaction spans reuse the molecule comparison schema and lift each entity row through the span
  tags. Canonicalization emits the lhs-anchored form.
- Canonical representatives may change between 0.x releases.

## Current state

### Molecule

| Level | Search graph | Initial partition | Orbit pruning | Branch order | Leaf key |
| --- | --- | --- | --- | --- | --- |
| Topology | Atoms, with bonds outside the largest bond-form class as colored subdivision vertices | Ranked entity colors | All generators | Backend canonical | Atom rows in search order, then sorted bond rows |
| Constitution | Adapter over the Constitution incidence graph | Constitution descriptors | All generators | Backend canonical | Adds sorted rows of the non-stereo domain |
| Structure | Adapter over the Full incidence graph | Structure partition: stereo descriptors over the constitution classes, repeated to a fixpoint when `para_stereo` is set | Generators that preserve the site, ligand frame, and configuration of every stereo entity | Backend canonical | Adds sorted stereo rows in the minimum frame |
| Full | As Structure, on the normalized molecule | Structure partition | None | Node id | Adds constraint blocks of the remapped and reframed molecule |

The table describes `tracked_canonicalize_at`. `canonical_key_at`, which `canonical_eq` and
`canonical_hash` use, repeats the setup of each level. At Structure it runs without the generator
filter and without pruning, in node-id order.

### Reaction span

- Level selection, entity colors, and incidence colors already lift the molecule's through the span
  tags.
- The search graph is the adapter over the incidence graph of the level, including at Topology.
- The initial partition is the ranked colors at Topology and the constitution descriptors at every
  other level. There is no stereo-sensitive refinement, and `para_stereo` is not read.
- There is no orbit pruning. Branch order is node id.
- The span tags have the positions `Unchanged = 0`, `Added = 1`, `Removed = 2`, `Modified = 3`, and
  the tag is the first component of every span color and row. Search order therefore places added
  entities between unchanged and removed ones. A separate step moves them behind the lhs-present
  entities of each kind.
- Each leaf derives an lhs-anchored remapping from the search order, remaps the whole span, reframes
  it by its representative action, and reads the key from the result in stored order. Non-atom
  entities therefore take their canonical order from the search partition, where molecules order
  them by row.
- `tracked_canonicalize_at` normalizes at every level; `canonical_key_at` normalizes at Full only.
- Both operations share one search.

### Measurements

Commit 5049b9796, one quick Criterion run, single-threaded, Nauty, Apple M1 Pro. `canonical_eq` on
an equal pair; each span is the molecule's reaction that adds one nitrogen atom.

| Case | Molecule | Reaction span |
| --- | ---: | ---: |
| `feature_free_connected` | 75.5 µs | 4.47 ms |
| `feature_free_disconnected` | 77.6 µs | 5.97 ms |
| `symmetry_heavy_radicals` | 60.4 µs | 22.7 ms |

On the small stereo cases, one molecule `canonicalize` and one span `canonical_eq`, which runs two
key searches, take 71.7 µs and 80.3 µs for `meso_dichlorobutane`, 134.6 µs and 108.3 µs for
`para_stereo_trichloropentane`, and 427 µs and 1.05 ms for `para_stereo_cascade`. The operations
differ, so this is not a like-for-like comparison. It shows that the molecule search has no per-call
advantage on these cases. The effect of the change on them is to be measured and is not predicted
here.

PR 29 reports 38 to 315 times faster canonicalization of the identity spans of the three
feature-free cases. Those numbers come from its older base and were not reproduced.

## Target

### Search by description level

At each description level the reaction-span search is the molecule search of that level: the same
search graph, initial partition, generators for orbit pruning, and leaf key as in the molecule table
above. Entity values are replaced by entity spans. Normalization happens at Full only.

| Step | Molecule | Reaction span |
| --- | --- | --- |
| Entity and incidence colors | Normalized inherent fields; electron contributions on aromatic and multicenter incidences | The same values per side under the span tag, as today |
| Topology search input | Molecular graph | Union graph |
| Stereo descriptor | Minimum over the frame permutations of the kind of site class, ligand classes, and configuration | Minimum over the same permutations of site class, ligand classes, lhs configuration, and rhs configuration, under the span tag |
| Stereo generator filter | One frame action carries the configuration to that of the image | The same on each side, by one frame action |
| Leaf rows | Entity rows of the comparison schema | The same rows lifted through the span tag |
| Order within an entity kind | Atoms in search order; other kinds by row | The same; the span tag order places lhs-present entities first |
| Constraint blocks at Full | Read from the remapped and reframed molecule | Read from the remapped and reframed span |

Notes on the lift:

- **Union graph.** Two bond entries on the same atom pair cannot share a color. An equal color means
  an equal span tag, which places both entries on one side, and each side is a molecule. Direct
  edges from the largest bond color class therefore keep the search graph simple.
- **Stereo kind.** Span integrity requires two kinded sides of a modified stereo entity to agree on
  the kind, so one set of frame permutations acts on both configurations. An undetermined side is
  unchanged by every permutation. When both sides are undetermined, the ligand classes are sorted,
  as in the molecule descriptor.
- **Generators and span tags.** Span tags are part of the colors, so every generator maps an entity
  to one with the same tag.

### Span tag order

The span tags are revised to `Unchanged = 0`, `Removed = 1`, `Modified = 2`, `Added = 3`. Only the
position of `Added` is forced: the lhs-anchored form places added entities after the lhs-present
entities of each kind. The three lhs-present tags keep the relative order they have today.
Constraint spans use `Unchanged`, `Removed`, and `Added`.

The tag is the first component of every span color, stereo descriptor, and row. The cells of the
initial partition therefore hold the lhs-present entities of each kind before the added ones.
Refinement and individualization split cells in place, and rows are sorted with the tag first.
Every leaf order is then lhs-anchored, and the span derives its remapping from the entity order as
the molecule does. The separate lhs-anchoring step is removed.

### Full level

At Full the leaf key is the structure key followed by the constraint blocks. A generator that
preserves the structure key preserves the complete key when it also carries every constraint.

- No inline constraint of an atom, bond, dative bond, aromatic system, multicenter bond, or
  noncovalent bond depends on a participant frame. #T and #C hold the states undetermined, negative,
  and positive. For these entity kinds the condition is an equal constraint set on the entity and
  its image.
- Ligand-symmetry, fluxionality, and topicity constraints of stereo entities name ligand positions.
  `StereoAtomForm` and `StereoBondForm` transport them together with the configuration. Comparing
  the whole transported form with the form of the image covers them with the frame action the
  stereo filter already derives.
- An aggregate with molecule-level constraints is searched exactly, as today.
- For a span, each condition holds on each side.

No Full-level search prunes today. This work tries pruning over the generators that pass, for
`ReactionSpan` and `Molecule` alike. Branch order follows the existing rule: backend canonical order
where the search requests backend automorphisms, node-id order otherwise.

### Key search

Reaction-span `tracked_canonicalize_at` and `canonical_key_at` share one search per level, as they
do today. `canonical_eq` and `canonical_hash` of spans and reactions therefore prune as
`canonicalize` does.

### Molecule side

`Molecule` keeps its two search setups. Two changes are tried:

- the Structure branch of `Molecule::canonical_key_at` prunes with the stereo generator filter, as
  `canonicalize_structure` does;
- Full-level pruning, in `canonicalize_full` and in the Full branch of `Molecule::canonical_key_at`.

### Code boundary

The molecule functions move into a submodule `canonicalize::molecule` and the reaction-span
functions into a parallel submodule `canonicalize::reaction_span`. A function moves when it takes a
`Molecule` or a `ReactionSpan`. The impl blocks of each aggregate move with its functions, and the
`Reaction` impl blocks sit with the span's. Functions in the two submodules carry the same names
without a prefix. Level selection, the comparison schema, and the search, partition, and adapter
code stay in `canonicalize` and are used as they are. The tests are split the same way.

`CompactTopologyCarrier` is renamed `TopologySearchInput`. Each submodule builds it with a function
`topology_search_input`.

No generic search over both aggregates, shared private trait, or merged function is introduced.
Parallel code between the molecule and span sides, and between the key search and the
canonical-form search, is accepted. Merging is left to an adversarial code review after this work is
running. Apart from the move and the rename, existing molecule functions change only for the two
trials above.

## Consequences

- Canonical representatives, tracked remappings, and canonical hashes of reaction spans and
  reactions change. The sources are the order of non-atom entities, the stereo-sensitive initial
  partition, and the reading of `para_stereo`.
- Each search remains an exact minimum of the typed key over a renumbering-invariant set of leaves,
  so `canonical_eq` relates the same pairs. Tests must establish this.
- Canonical forms of molecules do not change. The tracked remapping of a symmetric molecule with
  constraints may change to another remapping that yields the same canonical form.
- The exact cases in `umol-graph-ir/tests/canonicalization/reaction_span.rs` record span
  representatives and may change.
- The span tag positions in the data-type guide and the exact tag-position test change with the tag
  order. Every span comparison key changes with them.
- Molecule `canonical_eq` and `canonical_hash` return the same results. Their cost at Structure and
  Full changes if the two trials are kept.

## Tests and benchmarks

- Pruned and exhaustive span search select the same key at every level and under both branch
  orders.
- Span cases cover added, removed, and modified entities of each kind; stereo entities with a
  changed configuration, with one undetermined side, and present on one side only; and inline
  constraints that differ between sides.
- The span stereo descriptor is invariant under frame permutations acting on both sides.
- At Full, cases cover inline constraints on each entity kind, frame-relative constraints on stereo
  entities, and molecule-level constraints, for both aggregates.
- Canonical spans are lhs-anchored at every level.
- The molecule key search selects the same key pruned and exhaustive at Structure and Full.
- The existing reaction-span and reaction properties hold under both `para_stereo` settings.
- Benchmarks cover span `canonicalize`, `canonical_hash`, and `canonical_eq` over the corpus, with
  non-identity spans and inline-constraint cases, compared with `main` case by case, including the
  small stereo cases.
- Benchmarks cover molecule `canonical_eq` and `canonical_hash` on the stereo and inline-constraint
  cases, compared with `main`.

### Reference comparison

The current exhaustive span search is the reference. It runs from a separate checkout of `main`
under `scratch/`. No second version of the search is kept in the tree.

- A temporary test in the reference checkout runs `reaction_span_canonical_pair_strategy` with a
  fixed seed. It writes each pair as DSL text with the reference `canonical_eq` answer for both
  `para_stereo` settings.
- Hand-written pairs are added for symmetric spans, stereo that differs between sides, and inline
  constraints.
- A standalone program under `scratch/`, built against the working tree, reads the file, calls the
  public `canonical_eq` on each pair, and compares the answers.
- The report gives the number of pairs by answer and prints every mismatching pair as text.
- Timing uses the same Criterion cases in both checkouts.

Only equality answers are compared, because representatives and keys change. The comparison covers
the pairs in the file and nothing beyond them.

## Search study

The Structure search of the plan corrects no error: the search it replaces is exact, and the
equality answers are unchanged. It measured slower than that search on five of the six stereo cases
of the corpus, and it costs what the molecule search costs (S3d). This study compares search
variants on reaction spans and molecules. Its outcome decides the search of each level, and the
staged plan is revised from it.

### Variants

A search is described by three independent choices. Each value of a choice has a name, and a
combination is named by its partition and its enumeration, as in constitution-exhaustive.

**Initial partition**

- **Constitution partition.** The constitution descriptors, refined. A stereo entity enters with its
  kind; its configuration does not enter.
- **Structure partition.** Constitution classes from a separate refinement over the Constitution
  incidence graph, then stereo descriptors over these classes, refined again. The para-stereo rounds
  are a context option of this partition.

**Leaf enumeration**

- **Exhaustive.** Every leaf is visited, in node-id order, without a backend call.
- **Pruned.** One branch per orbit of the generators that preserve every stereo entity, with a
  backend call at every branching node and backend canonical order.

**Leaf key**

- **Remapped leaf.** The aggregate is remapped and reframed at each leaf, and the key is read from
  the result.
- **Row leaf.** The key rows are built from the leaf order.

### Combinations in the code

| Code | Partition | Enumeration | Leaf |
| --- | --- | --- | --- |
| Reaction span at Structure, before S3d | Constitution | Exhaustive | Remapped |
| Reaction span at Structure, S3d | Structure | Pruned | Remapped |
| Molecule `canonicalize` at Structure | Structure | Pruned | Row |
| Molecule `canonical_eq` and `canonical_hash` at Structure | Structure | Exhaustive | Row |
| Molecule at Full | Structure | Exhaustive | Row for the structure key, remapped for the constraint blocks |
| Reaction span at Full | Constitution | Exhaustive | Remapped |

The reaction-span operations share one search, so its rows hold for `canonicalize`, `canonical_eq`,
and `canonical_hash`.

### Candidates

None of the following is implemented.

- **Constitution-pruned.** Pruning with the stereo generator check from the constitution partition,
  without stereo descriptors.
- **Staged partition.** The constitution partition, with stereo descriptors computed only when
  stereo entities remain in shared cells.
- **Switched enumeration.** Exhaustive when a bound on the leaf count from the refined partition is
  small, pruned otherwise.
- **Root-pruned enumeration.** One backend call at the root. Deeper nodes prune with the root
  generators that fix the individualized nodes.
- **Backend labelling.** The order is taken from the backend's canonical labelling instead of the
  minimum of the typed key. The backend then defines the representative. Its exactness is not
  checked.
- **Deferred leaf key.** The key of a leaf is built when the leaf is first compared with another
  one. A `canonicalize` whose search ends with a single leaf then builds no key.

### Cases

The study uses the nine cases of the benchmark corpus that run at Structure: its six stereo cases
and three cases added for the study. Each added case has stereo atoms and a symmetric stereo-blind
graph, and every automorphism of that graph preserves the stereo atoms.

| Case | Compound | Atoms | Stereo atoms | Automorphisms of the stereo-blind graph |
| --- | --- | ---: | ---: | ---: |
| `symmetric_stereo_triisopropanolamine` | Triisopropanolamine as a heavy-atom skeleton, three equally configured stereocentres | 13 | 3 | 6 |
| `symmetric_stereo_butanol` | Butan-2-ol with explicit hydrogens | 15 | 1 | 72 |
| `symmetric_stereo_dimethylbutanol` | 3,3-Dimethylbutan-2-ol with explicit hydrogens | 21 | 1 | 7,776 |

A span is the reaction of the molecule that adds one nitrogen atom.

### Method

- **Code paths.** Criterion medians with a 1-second warm-up and 3 seconds of measurement, Nauty,
  single-threaded, Apple M1 Pro, `para_stereo` off. The constitution-exhaustive span search runs
  from a worktree at commit b8102f5a3, the structure-pruned span search from the S3d working tree,
  and the molecule searches from the working tree.
- **Stages.** Temporary ignored tests in that worktree call the stage functions of each aggregate in
  the order of the library, time each stage, and read the statistics of the search. A time is the
  smallest mean of five batches of about 0.1 s. The four combinations of partition and enumeration
  are composed from the existing functions. Constitution-pruned for both aggregates and
  constitution-exhaustive for the molecule exist only in these tests. The tests are not committed.
- **Check.** The composed totals lie 0% to 6% below the timings of the corresponding code paths, and
  11% below on one case.

### Results

Visited leaves, with backend calls in parentheses. The counts are the same for a molecule and its
span.

| Case | Constitution-exhaustive | Constitution-pruned | Structure-exhaustive | Structure-pruned |
| --- | ---: | ---: | ---: | ---: |
| `tetrahedral_stereo` | 1 | 1 (0) | 1 | 1 (0) |
| `cis_trans_stereo_bond` | 1 | 1 (0) | 1 | 1 (0) |
| `mixed_atom_and_bond_stereo` | 1 | 1 (0) | 1 | 1 (0) |
| `meso_dichlorobutane` | 2 | 2 (1) | 1 | 1 (0) |
| `para_stereo_trichloropentane` | 2 | 2 (1) | 2 | 2 (1) |
| `para_stereo_cascade` | 8 | 8 (5) | 2 | 2 (1) |
| `symmetric_stereo_triisopropanolamine` | 6 | 1 (2) | 6 | 1 (2) |
| `symmetric_stereo_butanol` | 72 | 1 (5) | 72 | 1 (5) |
| `symmetric_stereo_dimethylbutanol` | 7,776 | 1 (8) | 7,776 | 1 (8) |

Molecule `canonicalize`, composed from the stage times, row leaf:

| Case | Constitution-exhaustive | Constitution-pruned | Structure-exhaustive | Structure-pruned |
| --- | ---: | ---: | ---: | ---: |
| `tetrahedral_stereo` | 21.3 µs | 21.3 µs | 40.0 µs | 40.0 µs |
| `cis_trans_stereo_bond` | 20.4 µs | 20.4 µs | 32.2 µs | 32.2 µs |
| `mixed_atom_and_bond_stereo` | 40.2 µs | 40.0 µs | 70.9 µs | 70.9 µs |
| `meso_dichlorobutane` | 52.3 µs | 55.5 µs | 67.0 µs | 67.0 µs |
| `para_stereo_trichloropentane` | 75.6 µs | 79.7 µs | 124.2 µs | 128.5 µs |
| `para_stereo_cascade` | 778.9 µs | 845.8 µs | 395.1 µs | 406.7 µs |
| `symmetric_stereo_triisopropanolamine` | 214.9 µs | 75.3 µs | 270.0 µs | 130.1 µs |
| `symmetric_stereo_butanol` | 1.27 ms | 86.2 µs | 1.33 ms | 116.8 µs |
| `symmetric_stereo_dimethylbutanol` | 168 ms | 188.2 µs | 169 ms | 226.6 µs |

Reaction-span `canonicalize`, composed from the stage times, remapped leaf:

| Case | Constitution-exhaustive | Constitution-pruned | Structure-exhaustive | Structure-pruned |
| --- | ---: | ---: | ---: | ---: |
| `tetrahedral_stereo` | 23.5 µs | 23.6 µs | 43.4 µs | 43.4 µs |
| `cis_trans_stereo_bond` | 27.0 µs | 27.0 µs | 40.2 µs | 40.3 µs |
| `mixed_atom_and_bond_stereo` | 45.9 µs | 46.0 µs | 78.7 µs | 78.8 µs |
| `meso_dichlorobutane` | 45.4 µs | 49.2 µs | 68.0 µs | 68.0 µs |
| `para_stereo_trichloropentane` | 61.0 µs | 65.7 µs | 113.1 µs | 118.1 µs |
| `para_stereo_cascade` | 506.9 µs | 572.3 µs | 350.1 µs | 362.5 µs |
| `symmetric_stereo_triisopropanolamine` | 169.5 µs | 75.9 µs | 227.7 µs | 135.1 µs |
| `symmetric_stereo_butanol` | 1.40 ms | 97.9 µs | 1.45 ms | 131.2 µs |
| `symmetric_stereo_dimethylbutanol` | 190 ms | 208.7 µs | 192 ms | 250.7 µs |

Stage times of the molecule. The classes pass and the stereo descriptors are parts of the structure
partition. The remapped leaf is the span's.

| Case | Constitution partition | Structure partition | Classes pass | Stereo descriptors | Row leaf | Remapped leaf | Backend call |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `tetrahedral_stereo` | 4.1 µs | 22.9 µs | 6.5 µs | 14.3 µs | 8.7 µs | 7.4 µs | 1.2 µs |
| `cis_trans_stereo_bond` | 4.9 µs | 16.7 µs | 8.1 µs | 6.1 µs | 6.0 µs | 8.2 µs | 1.2 µs |
| `mixed_atom_and_bond_stereo` | 9.6 µs | 40.3 µs | 15.5 µs | 20.2 µs | 14.2 µs | 13.1 µs | 1.9 µs |
| `meso_dichlorobutane` | 5.6 µs | 39.7 µs | 7.6 µs | 28.0 µs | 15.4 µs | 9.4 µs | 1.5 µs |
| `para_stereo_trichloropentane` | 8.7 µs | 57.3 µs | 10.8 µs | 42.0 µs | 22.3 µs | 11.8 µs | 3.3 µs |
| `para_stereo_cascade` | 28.9 µs | 175.7 µs | 33.4 µs | 122.9 µs | 69.8 µs | 33.7 µs | 10.0 µs |
| `symmetric_stereo_triisopropanolamine` | 11.0 µs | 65.8 µs | 18.6 µs | 42.3 µs | 24.3 µs | 15.1 µs | 5.4 µs |
| `symmetric_stereo_butanol` | 11.9 µs | 42.5 µs | 23.4 µs | 15.0 µs | 13.0 µs | 14.2 µs | 7.7 µs |
| `symmetric_stereo_dimethylbutanol` | 15.5 µs | 52.7 µs | 32.2 µs | 15.2 µs | 15.9 µs | 18.3 µs | 21.2 µs |

Code paths, measured with Criterion:

| Case | Span `canonicalize`, constitution-exhaustive | Span `canonicalize`, structure-pruned | Molecule `canonicalize`, structure-pruned | Molecule `canonical_hash`, structure-exhaustive |
| --- | ---: | ---: | ---: | ---: |
| `tetrahedral_stereo` | 24.2 µs | 45.2 µs | 40.6 µs | 39.0 µs |
| `cis_trans_stereo_bond` | 27.4 µs | 41.0 µs | 33.1 µs | 31.7 µs |
| `mixed_atom_and_bond_stereo` | 46.5 µs | 79.8 µs | 72.2 µs | 71.5 µs |
| `meso_dichlorobutane` | 46.7 µs | 69.3 µs | 68.0 µs | 65.9 µs |
| `para_stereo_trichloropentane` | 62.3 µs | 124.9 µs | 130.6 µs | 124.4 µs |
| `para_stereo_cascade` | 523.8 µs | 382.3 µs | 415.4 µs | 401.1 µs |
| `symmetric_stereo_triisopropanolamine` | 171.9 µs | 143.5 µs | 132.1 µs | 270.1 µs |
| `symmetric_stereo_butanol` | 1.40 ms | 146.8 µs | 117.9 µs | 1.29 ms |
| `symmetric_stereo_dimethylbutanol` | 194 ms | 262.0 µs | 233.2 µs | 166 ms |

### Findings

- **Enumeration.** Exhaustive enumeration visits one leaf per automorphism of the stereo-blind
  graph on the three symmetric cases. It takes 170 to 270 µs with 6 automorphisms, 1.3 to 1.5 ms
  with 72, and 170 to 190 ms with 7,776. Pruned enumeration visits one leaf on these cases and takes
  75 to 250 µs. Where the refined partition is discrete, pruning makes no backend call and costs
  nothing; this holds on three cases. Where a backend call finds nothing to prune, it costs 3 to
  5 µs on `meso_dichlorobutane` and `para_stereo_trichloropentane`.
- **Partition.** The structure partition costs 12 to 147 µs more than the constitution partition.
  The stereo descriptors take about 14 µs for a tetrahedral stereo atom and 6 µs for a cis-trans
  stereo bond, and up to twice that for a modified stereo entity of a span. The separate classes
  pass takes 6 to 33 µs. The structure partition lowers the number of visited leaves on two cases,
  `meso_dichlorobutane` and `para_stereo_cascade`, and lowers the total time on the second only.
- **Structure-pruned against constitution-exhaustive.** This is the change of S3d, measured on the
  span code paths. It is 14 to 63 µs slower on the five small stereo cases of the corpus, a factor
  of 1.5 to 2.0. It is faster by a factor of 1.4 on `para_stereo_cascade` and by factors of 1.2, 10,
  and 740 on the symmetric cases with 6, 72, and 7,776 automorphisms.
- **Constitution-pruned.** It is as fast as constitution-exhaustive on three cases, 3 to 5 µs slower
  on two, and 9% to 13% slower on `para_stereo_cascade`. On the three symmetric cases it is faster
  by factors of 2 to 3, 14, and 900. It is faster than structure-pruned on eight cases, by 12 to
  59 µs, and slower on `para_stereo_cascade`: 846 against 407 µs for the molecule and 572 against
  362 µs for the span.
- **Molecule key search.** `canonical_eq` and `canonical_hash` of molecules run
  structure-exhaustive. `canonical_hash` takes 270 µs, 1.29 ms, and 166 ms on the three symmetric
  cases, against 132, 118, and 233 µs for `canonicalize`.
- **Leaf.** On the four cases with two or more stereo atoms, the row leaf of the molecule costs 1.6
  to 2.1 times the remapped leaf of the span. On the other five cases it costs 0.7 to 1.2 times as
  much.
- **Para-stereo.** With `para_stereo` on, the structure partition runs a second round on three
  cases and costs 102, 318, and 111 µs for the molecule, against 57, 176, and 66 µs. The second
  round removes one further leaf, on `para_stereo_cascade`.
- **Representative.** The two partitions select different minimum keys on four of the nine cases, so
  a change of partition changes canonical representatives. The data-type guide prescribes the
  stereo-sensitive refinement of the structure partition.

### Trial

A trial in the worktree at commit b8102f5a3 makes two changes. It is not committed.

- **Constitution-pruned at Structure** for both aggregates, in `canonicalize`, `canonical_eq`, and
  `canonical_hash`.
- **Deferred leaf key** for reaction spans, at every level. The key of a leaf is built when the leaf
  is first compared with another one, so a `canonicalize` whose search ends with a single leaf
  builds no key. `canonical_eq` and `canonical_hash` return the key and always build it.

The deferred key does not carry over to the molecule as it stands. The row leaf of the molecule
produces the order of the non-atom entities together with the key, so `canonicalize` needs the leaf
even when it is the only one. Deferring the molecule key would require taking that order from the
search, as the span does, and reading the key from the remapped molecule in stored order.

In the worktree the library tests, the two integration tests, and the canonicalization properties
pass without a change to an expected value. The reference comparison reports no mismatch on the
10,028 pairs. The canonical forms and tracked remappings of all 20,056 spans of the pair file are
identical to those of the committed search.

Reaction spans, Criterion medians. The committed search is constitution-exhaustive at Structure and
exhaustive at Full, and prunes below Structure.

| Case | Level | `canonicalize` committed | `canonicalize` trial | `canonical_eq` committed | `canonical_eq` trial |
| --- | --- | ---: | ---: | ---: | ---: |
| `ordinary_naphthalene` | Topology or Constitution | 37.8 µs | 25.9 µs | 62.8 µs | 61.0 µs |
| `disconnected_rings` | Topology or Constitution | 57.3 µs | 44.3 µs | 101.3 µs | 97.6 µs |
| `feature_free_connected` | Topology or Constitution | 60.4 µs | 46.8 µs | 113.9 µs | 111.4 µs |
| `feature_free_disconnected` | Topology or Constitution | 61.2 µs | 46.4 µs | 114.1 µs | 110.5 µs |
| `symmetry_heavy_radicals` | Topology or Constitution | 40.1 µs | 32.4 µs | 71.9 µs | 70.1 µs |
| `tetrahedral_stereo` | Structure | 24.2 µs | 15.4 µs | 36.4 µs | 35.3 µs |
| `cis_trans_stereo_bond` | Structure | 27.4 µs | 18.2 µs | 42.4 µs | 41.1 µs |
| `mixed_atom_and_bond_stereo` | Structure | 46.5 µs | 31.1 µs | 73.3 µs | 71.4 µs |
| `meso_dichlorobutane` | Structure | 46.7 µs | 49.5 µs | 78.8 µs | 81.7 µs |
| `para_stereo_trichloropentane` | Structure | 62.3 µs | 65.8 µs | 105.5 µs | 112.0 µs |
| `para_stereo_cascade` | Structure | 523.8 µs | 563.4 µs | 1.01 ms | 1.12 ms |
| `symmetric_stereo_triisopropanolamine` | Structure | 171.9 µs | 59.0 µs | 326.9 µs | 127.7 µs |
| `symmetric_stereo_butanol` | Structure | 1.40 ms | 82.1 µs | 2.82 ms | 177.8 µs |
| `symmetric_stereo_dimethylbutanol` | Structure | 194 ms | 185.3 µs | 382 ms | 398.8 µs |
| `overlay_heavy` | Full | 56.7 µs | 35.8 µs | 93.8 µs | 89.2 µs |
| `frame_relative_stereo_constraint` | Full | 26.4 µs | 16.2 µs | 43.5 µs | 41.7 µs |

Molecules, Criterion medians. The committed search is structure-pruned for `canonicalize` and
structure-exhaustive for `canonical_hash`.

| Case | `canonicalize` committed | `canonicalize` trial | `canonical_hash` committed | `canonical_hash` trial |
| --- | ---: | ---: | ---: | ---: |
| `tetrahedral_stereo` | 40.6 µs | 21.2 µs | 39.0 µs | 17.6 µs |
| `cis_trans_stereo_bond` | 33.1 µs | 20.2 µs | 31.7 µs | 17.5 µs |
| `mixed_atom_and_bond_stereo` | 72.2 µs | 40.0 µs | 71.5 µs | 35.8 µs |
| `meso_dichlorobutane` | 68.0 µs | 55.5 µs | 65.9 µs | 51.7 µs |
| `para_stereo_trichloropentane` | 130.6 µs | 80.8 µs | 124.4 µs | 74.5 µs |
| `para_stereo_cascade` | 415.4 µs | 830.5 µs | 401.1 µs | 813.6 µs |
| `symmetric_stereo_triisopropanolamine` | 132.1 µs | 75.7 µs | 270.1 µs | 70.0 µs |
| `symmetric_stereo_butanol` | 117.9 µs | 85.7 µs | 1.29 ms | 82.5 µs |
| `symmetric_stereo_dimethylbutanol` | 233.2 µs | 190.3 µs | 166 ms | 184.5 µs |

- **Reaction-span `canonicalize`.** It is 19% to 39% faster on the ten cases whose committed search
  already ends with one leaf, at every level, from the deferred key. It is 6% slower on
  `meso_dichlorobutane` and `para_stereo_trichloropentane` and 8% slower on `para_stereo_cascade`,
  from backend calls that prune nothing. On the three symmetric stereo cases it takes 59, 82, and
  185 µs, against 172 µs, 1.40 ms, and 194 ms.
- **Reaction-span `canonical_eq`.** It is within 5% of the committed search on those ten cases, 4%
  to 11% slower on the three cases above, and takes 128, 178, and 399 µs on the symmetric stereo
  cases, against 327 µs, 2.82 ms, and 382 ms.
- **Molecule.** `canonicalize` is 18% to 48% faster on eight cases and `canonical_hash` 22% to 55%
  faster on five; on the three symmetric cases `canonical_hash` takes 70, 83, and 185 µs, against
  270 µs, 1.29 ms, and 166 ms. On `para_stereo_cascade` both operations take twice as long, about
  0.82 ms against 0.41 ms.
- **Composed totals.** The trial timings agree with the totals composed from the stage times to
  within 6%.

### Trial of the deferred leaf for molecules

A second trial in the same worktree gives the molecule the leaf of the span at all four levels. The
entity order of a leaf is the search order, the key is read from the remapped and reframed molecule
with every row in stored order, and the key is deferred. At Topology the search orders atoms only,
so the bonds follow in the order of their endpoint pairs. Partition and enumeration are those of the
first trial: constitution-pruned at Structure, unchanged at the other levels. The trial is not
committed.

The stored-order key and the two functions around it take about 190 lines. The row-leaf functions
of the molecule and the minimum stereo frame, about 430 lines, are no longer called.

In the worktree the two integration tests and the canonicalization properties pass, and the
reference comparison of spans reports no mismatch. `canonical_eq` and the equality of
`canonical_hash` give the same answer as the committed search on 20,056 pairs of molecules, the lhs
and rhs molecules of the span pairs, under both `para_stereo` settings. One library test changes its
expected value: `test_canonicalize_constitution` asserts an exact canonical molecule, whose bonds
now stand in search order instead of the order of their endpoint pairs. Its atoms keep their
numbers.

Molecules, Criterion medians. The row-leaf columns are the first trial, which ran the Structure
cases only. The committed timings of the Topology and Full cases are those of the S0b baseline; the
molecule search has not changed since.

| Case | Level | `canonicalize` committed | Row leaf | Deferred leaf | `canonical_hash` committed | Row leaf | Deferred leaf |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `ordinary_naphthalene` | Topology | 21.3 µs |  | 17.4 µs | 21.0 µs |  | 22.6 µs |
| `disconnected_rings` | Topology | 34.4 µs |  | 29.2 µs | 33.3 µs |  | 35.5 µs |
| `feature_free_connected` | Topology | 37.8 µs |  | 32.1 µs | 38.6 µs |  | 40.5 µs |
| `feature_free_disconnected` | Topology | 40.4 µs |  | 35.2 µs | 41.4 µs |  | 43.3 µs |
| `symmetry_heavy_radicals` | Topology | 31.1 µs |  | 28.6 µs | 30.8 µs |  | 33.1 µs |
| `tetrahedral_stereo` | Structure | 40.6 µs | 21.2 µs | 12.3 µs | 39.0 µs | 17.6 µs | 15.5 µs |
| `cis_trans_stereo_bond` | Structure | 33.1 µs | 20.2 µs | 14.2 µs | 31.7 µs | 17.5 µs | 18.2 µs |
| `mixed_atom_and_bond_stereo` | Structure | 72.2 µs | 40.0 µs | 25.4 µs | 71.5 µs | 35.8 µs | 32.9 µs |
| `meso_dichlorobutane` | Structure | 68.0 µs | 55.5 µs | 40.6 µs | 65.9 µs | 51.7 µs | 37.3 µs |
| `para_stereo_trichloropentane` | Structure | 130.6 µs | 80.8 µs | 54.4 µs | 124.4 µs | 74.5 µs | 50.4 µs |
| `para_stereo_cascade` | Structure | 415.4 µs | 830.5 µs | 496.4 µs | 401.1 µs | 813.6 µs | 478.0 µs |
| `symmetric_stereo_triisopropanolamine` | Structure | 132.1 µs | 75.7 µs | 50.6 µs | 270.1 µs | 70.0 µs | 59.2 µs |
| `symmetric_stereo_butanol` | Structure | 117.9 µs | 85.7 µs | 71.7 µs | 1.29 ms | 82.5 µs | 81.1 µs |
| `symmetric_stereo_dimethylbutanol` | Structure | 233.2 µs | 190.3 µs | 170.9 µs | 166 ms | 184.5 µs | 184.1 µs |
| `overlay_heavy` | Full | 91.8 µs |  | 67.2 µs | 88.1 µs |  | 75.0 µs |
| `frame_relative_stereo_constraint` | Full | 46.7 µs |  | 34.1 µs | 44.5 µs |  | 39.1 µs |

- **`canonicalize`.** It is 8% to 19% faster than the committed search on the five Topology cases,
  27% faster on the two Full cases, and 1.4 to 3.3 times as fast on eight Structure cases. On
  `para_stereo_cascade` it takes 496 µs against 415 µs; the row leaf took 831 µs there.
- **`canonical_hash`.** It is 5% to 8% slower than the committed search on the five Topology cases,
  where the key is always built and the remapped leaf costs more than the topology row leaf.
  `canonical_eq` is 7% to 13% slower there. At Structure it is as fast as the row leaf or up to 41%
  faster, apart from 4% on `cis_trans_stereo_bond`, and it is 12% to 15% faster than the committed
  search on the two Full cases.
- **Order of non-atom entities.** With the search order, a canonical molecule at Constitution and
  above lists bonds, overlays, and stereo entities in the order of the search partition, as a
  canonical span does, while a molecule at Topology keeps its bonds in the order of their endpoint
  pairs. The order then depends on the level, which the content of the molecule selects.

### Variants of the second trial

Four variants followed in the same worktree, each on top of the one before.

- **Order by participants.** The first variant orders the non-atom entities of every level by
  their participants: bonds and noncovalent bonds by endpoint pair, dative bonds by donors and
  acceptor, aromatic systems and multicenter bonds by their atoms, stereo atoms by site, and stereo
  bonds by site bond. Entities of one kind on the same participants keep their search order. This
  is the order molecules have today and a function of the atom numbering alone. With it every
  library test passes without a changed expected value, and the integration tests, the properties,
  and both equality comparisons hold. Its cost is below the measurement noise on fifteen of the
  sixteen cases, between -1.0 and +0.7 µs, and 9 to 11 µs, about 2%, on `para_stereo_cascade`,
  which compares eight leaves.
- **Order by participants for spans.** The second variant orders, for reaction spans, the non-atom
  entities of each kind with the lhs-present ones first and each group by its participants. The
  canonical form stays lhs-anchored. The properties pass and the reference comparison reports no
  mismatch. The canonical form of 3,666 of the 20,056 spans of the pair file changes, in the order
  of non-atom entities and the site indices that follow from it; the atoms keep their numbers. Three
  exact cases change their expected value: the stereo-bond row of
  `test_reaction_span_tracked_canonicalize_at` and two cases of
  `tests/canonicalization/reaction_span.rs`. Against the span trial with search order,
  `canonicalize` and `canonical_hash` take 0.3 to 3.3 µs more on fifteen cases, 1% to 8%, and 13 to
  45 µs more on `para_stereo_cascade`, 4% to 8%. The two timings come from separate runs.
- **Topology key read directly.** The third variant reads the Topology key of a molecule from the
  entity order without remapping. Topology has no participant frame, so this key equals the key of
  the remapped molecule; a debug assertion checked the equality on every Topology canonicalization
  of the library tests and the properties. With it `canonical_hash` and `canonical_eq` take 0.94 to
  1.02 times the time of the committed search on the five Topology cases, and `canonicalize` stays
  4% to 17% faster than the committed search.
- **Span search at Topology.** The fourth variant gives the span the two things the molecule has at
  Topology: the topology search input, with atoms and the bonds outside the largest bond color class
  as vertices, and the key read directly from the entity order. The same checks hold, with the debug
  assertion on the direct key, and no canonical form of the pair file changes relative to the
  variant with the order by participants. On four of the five Topology cases `canonicalize` is 16%
  to 24% faster than that variant and 34% to 45% faster than the committed search;
  `symmetry_heavy_radicals` is 3% slower than the variant and 17% faster than the committed search.
  `canonical_hash` and `canonical_eq` are 9% to 35% faster than both on all five. The span then
  takes 1.1 to 1.2 times as long as the molecule on these cases, with one atom more.

Reaction spans at Topology, Criterion medians, for the committed search, the third variant, and the
fourth:

| Case | `canonicalize` committed | Third | Fourth | `canonical_hash` committed | Fourth | `canonical_eq` committed | Fourth |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `ordinary_naphthalene` | 37.8 µs | 26.7 µs | 20.9 µs | 34.1 µs | 23.7 µs | 62.8 µs | 40.5 µs |
| `disconnected_rings` | 57.3 µs | 44.8 µs | 33.9 µs | 53.1 µs | 36.8 µs | 101.3 µs | 66.3 µs |
| `feature_free_connected` | 60.4 µs | 46.6 µs | 37.2 µs | 57.6 µs | 42.2 µs | 113.9 µs | 80.3 µs |
| `feature_free_disconnected` | 61.2 µs | 47.9 µs | 40.1 µs | 57.4 µs | 44.4 µs | 114.1 µs | 85.1 µs |
| `symmetry_heavy_radicals` | 40.1 µs | 32.1 µs | 33.2 µs | 37.6 µs | 33.5 µs | 71.9 µs | 63.9 µs |

### Trial of constraints as fields

A third trial in the same worktree, on top of the fourth variant, treats inline constraints and
positional values as the search treats fields. It is not committed.

- **Inline constraints.** No inline constraint selects Full. The level follows from the entities
  alone: atoms and bonds give Topology; a dative bond, aromatic system, multicenter bond, or
  noncovalent bond gives Constitution; a stereo entity gives Structure. The constraints of an entity
  enter its color and its key row at the level where the entity appears. A stereo entity is colored
  by its frame-independent constraints, and its row carries all of them.
- **Positional values.** Electron counts and dative roles are incidence colors. They now also split
  the search partition, which is built from the colors.
- **Frame-relative values.** The stereo generator check transports the whole stereo form, with
  configuration and constraints, by the frame action it derives, and compares it with the form of
  the image.
- **Full.** Full is reached with molecule-level constraints only. Both aggregates search it from the
  constitution partition with pruning. A generator is kept when the aggregate remapped by it equals
  the aggregate under framed equality.
- **Normalization.** A reaction span is normalized before the search at Full only.

The canonicalization properties pass. The reference comparison reports no mismatch on the 10,028
pairs, and `canonical_eq` and the equality of `canonical_hash` of molecules give the same answer as
the committed search on the 20,056 pairs. Of the 40,112 canonical span forms of the pair file,
counted under both `para_stereo` settings, 16 change.

Temporary tests renumber constrained aggregates and compare the canonical form, the hash, and
`canonical_eq` of each renumbering with the original, for the molecule and for its span. Every
renumbering gives the form and the hash of the original. Different placements of a constraint are
unequal, and a placement moved to an equivalent entity is equal.

| Aggregate | Placements | Renumberings of each |
| --- | ---: | ---: |
| Six-ring with valence constraints on atoms and aromatic constraints on bonds | 11 | 16 |
| Benzene and biphenyl with and without `#e6` on their aromatic systems; five-ring with electron counts 2, 1, 1, 1, 1 | 6 | 1 |
| 2,3-Dichlorobutane with stereogenicity and topicity constraints on one or both centres | 8 | 6 |
| Six-ring with one molecule-level charge-sum constraint | 6 | 15 |

Twenty-seven library tests change, and each asserts a rule the trial replaces. Nineteen expect a
constraint to select Full. Five expect a key or a color below Full to leave constraints out. One
asserts an exact canonical molecule, in which the constrained atom moves from first to second. One
compares the Constitution search with the minimum over every relabeling. One asserts a visited-leaf
count.

Criterion medians of the three scaling cases with `#v1` on every hydrogen, before and after:

| Case | Molecule `canonicalize` | Molecule `canonical_eq` | Span `canonicalize` | Span `canonical_eq` |
| --- | ---: | ---: | ---: | ---: |
| `feature_free_connected` | 2.19 ms → 35.0 µs | 4.39 ms → 81.1 µs | 2.67 ms → 37.4 µs | 5.29 ms → 89.1 µs |
| `feature_free_disconnected` | 2.80 ms → 37.8 µs | 5.66 ms → 86.2 µs | 3.41 ms → 39.0 µs | 6.87 ms → 93.7 µs |
| `symmetry_heavy_radicals` | 11.59 ms → 31.0 µs | 23.24 ms → 64.6 µs | 14.10 ms → 31.6 µs | 28.10 ms → 69.2 µs |

A constrained molecule takes 3% to 6% longer than the same molecule without constraints. Molecule
`canonicalize` of the temporary cases, as the smallest mean of five batches:

| Case | Before | After |
| --- | ---: | ---: |
| Benzene with `#e6` | 154 µs | 24 µs |
| Biphenyl with `#e6` on both aromatic systems | 195 µs | 55 µs |
| Biphenyl with `#e6` on one aromatic system | 197 µs | 48 µs |
| Five-ring with electron counts 2, 1, 1, 1, 1 | 67 µs | 18 µs |
| 2,3-Dichlorobutane with a constraint on a centre, seven placements | 88 µs | 22 to 51 µs |
| Six-ring, charge-sum constraint over the whole molecule | 122 µs | 29 µs |
| Six-ring, charge-sum constraint over atom 0, or atoms 0 and 3 | 123 µs | 108 µs |
| Six-ring, charge-sum constraint over atoms 0, 1, or 0, 2, or 0, 1, 2 | 123 µs | 171 to 186 µs |

- **Cases without constraints.** The 21 Topology benchmarks lie between 0.98 and 1.04 times the time
  before. Of 128 benchmarks over the corpus, 126 lie between 0.95 and 1.06; three that first
  measured 1.05 to 1.06 returned 1.00 in a second run.
- **Cases that ran at Full.** `overlay_heavy` takes 32.1 µs against 71.3 µs and
  `frame_relative_stereo_constraint` 14.4 µs against 35.0 µs.
- **Pruned Full.** It is four times as fast as the exhaustive search when the molecule-level
  constraint keeps the symmetry of the ring, and up to 1.5 times slower when the constraint breaks
  it. The generator check remaps the whole aggregate, which costs about as much as a leaf, and the
  leaves are then visited as well.
- **Normalization.** With the span normalized at Full only, span `canonicalize` takes 0.88 to 1.02
  times the time before on 35 benchmarks, 1 to 3 µs less on most. All canonical forms and tracked
  remappings of the pair file are identical before and after.

### Trial of colors

A fourth trial, in the worktree `scratch/244-color-trial` on the commit of S6b, replaces the encoded
keys with colors read from the forms. It is not committed. It changes `canonicalize` of molecules
only; `canonical_eq`, `canonical_hash`, reaction spans, and Full keep the committed search. It has
three steps, each on top of the one before.

- **Atoms and bonds.** The molecule is normalized once, in place, before the search. The color of an
  atom or a bond is a borrowed view of its whole form, and the colors are ranked into dense
  integers. The key of a leaf is the integer of each atom in the order of the leaf, then the
  endpoint images and the integer of each bond. The result is remapped and not reframed.
- **Partition.** The incidence descriptors of `constitution_partition_descriptors` and of
  `structure_partition_descriptors` are the exact incidence colors, so dative roles and electron
  counts split the search partition.
- **Overlays and stereo entities.** Constitution and Structure rank the colors of all eight kinds
  and of the incidences, build the partition from them, and read the key rows from the entity order
  of the leaf: the participants in their new numbers, in the representative order of the entity,
  and the integers. A stereo row holds the form transported by the representative action of its
  renumbered ligands and normalized. The stereo generator check transports the whole form. No inline
  constraint selects Full.

In debug builds the trial runs the committed path beside its own for every molecule without inline
constraints and compares the canonical form and the remapping. In the equality contract of the
property suite, `canonicalize` runs the trial path and `canonical_eq` the committed one.

| Step | Library tests | Compared with the committed path | Properties |
| --- | --- | --- | --- |
| Atoms and bonds | 12,138 pass | Equal in every Topology run | 555 pass |
| Partition | 12,136 pass, 2 change | Not compared | 555 pass |
| Overlays and stereo entities | 12,134 pass, 4 change | Equal in all 238 runs | 555 pass |

The four changed tests assert rules the trial replaces.
`test_canonicalize_constitution_entity_kind_minimum` compares the search with the minimum over every
relabeling. The modified aromatic system row of `test_canonical_candidate` expects 6 and 3 visited
leaves, where the search now visits 2 and 1. `test_canonicalize_constitution_excluded_data` and the
Full row of `test_molecule_canonicalize_description_level` expect an inline constraint to stay out
of the key below Full and to select Full.

Criterion medians of molecule `canonicalize`, in the worktree, before and after each step:

| Case | Before | Atoms and bonds |
| --- | ---: | ---: |
| `ordinary_naphthalene` | 17.2 µs | 12.7 µs |
| `disconnected_rings` | 29.3 µs | 24.0 µs |
| `feature_free_connected` | 33.1 µs | 26.3 µs |
| `feature_free_disconnected` | 36.3 µs | 29.6 µs |
| `symmetry_heavy_radicals` | 30.2 µs | 26.5 µs |
| Carbon path, 8 atoms | 11.5 µs | 8.1 µs |
| Carbon path, 16 atoms | 19.8 µs | 13.2 µs |
| Carbon path, 32 atoms | 43.0 µs | 28.3 µs |
| Carbon path, 64 atoms | 101.7 µs | 74.3 µs |
| Carbon path, 77 atoms | 133.4 µs | 100.1 µs |
| Carbon path, 128 atoms | 286.5 µs | 234.7 µs |

| Case | Before | Overlays and stereo entities |
| --- | ---: | ---: |
| `tetrahedral_stereo` | 13.0 µs | 8.3 µs |
| `cis_trans_stereo_bond` | 14.7 µs | 8.9 µs |
| `mixed_atom_and_bond_stereo` | 26.0 µs | 15.1 µs |
| `meso_dichlorobutane` | 41.9 µs | 20.1 µs |
| `para_stereo_trichloropentane` | 56.6 µs | 27.7 µs |
| `para_stereo_cascade` | 528.4 µs | 257.6 µs |
| `symmetric_stereo_triisopropanolamine` | 51.1 µs | 37.3 µs |
| `symmetric_stereo_butanol` | 73.1 µs | 58.8 µs |
| `symmetric_stereo_dimethylbutanol` | 177.5 µs | 161.2 µs |

- **Controls.** `canonical_eq` and `canonical_hash` of the five Topology cases stay within 3.3% in
  the first step. `canonical_eq` of the eleven cases above Topology and span `canonicalize` of two
  stay within 5.5% in the third.
- **Partition.** The second step moves 24 benchmarks by −4.9% to +4.7%, in no consistent direction.
  No benchmark case has dative roles or electron counts that distinguish otherwise equal atoms.
- **Full.** `overlay_heavy` and `frame_relative_stereo_constraint` carry a molecule-level constraint
  and run at Full, which the trial leaves alone. They take 67.6 µs against 68.2 µs and 34.8 µs
  against 34.6 µs. No case of the corpus runs at Constitution, so that level is checked by the
  tests and not timed.

Three readings of the code accompany the trial.

- **Spans and reactions.** `canonical_candidate` of the span runs every level with the steps of the
  molecule searches. Each span color is the molecule color lifted through the span by
  `normalized_entity_span_key`: the tag orders first, and a `Modified` entity whose two sides have
  equal keys counts as `Unchanged`. `ReactionSpan::tracked_canonicalize_at` already normalizes the
  span once before the search at every level. A reaction is converted to its span, the span is
  canonicalized, and the result is converted back; only the level selection is its own. With
  `PartialOrd` and `Ord` removed from `EntitySpan`, the workspace without `umol-py` compiles, and
  `umol-py` only constructs entity spans.
- **Molecule-level constraints.** `Constraints::normalize` sorts and deduplicates by the derived
  order. The constraints of a remapped and reframed aggregate are therefore already the sorted
  sequence that the Full key builds from their encodings.
- **Encoding.** `CanonicalKeyValue`, the encoding functions, and `CanonicalComparisonKey` are used
  in `canonicalize.rs`, its two submodules, and their tests only.

### Decisions

- **Structure search.** Molecules and reaction spans use constitution-pruned at Structure, for
  `canonicalize`, `canonical_eq`, and `canonical_hash`. `para_stereo_cascade` is slower with it than
  with structure-pruned, by a factor of two with the row leaf and of 1.2 with the leaf decided
  below. It is an edge case and does not hold up the gain on the other cases.
- **No selectable variant.** The search is fixed. No setting selects a partition or an enumeration,
  and no setting changes the canonical form.
- **Representatives.** A change of canonical representatives is acceptable: canonicalization is
  consistent within a version only.
- **Guide.** The data-type guide is revised where a faster search requires it. Its prescription of
  the stereo-sensitive refinement is replaced with the constitution partition at Structure. Its
  rules that constraints do not distinguish structural orbits, and that complete canonicalization
  minimizes the structure key before the constraint key, are replaced by the decisions on colors
  and inline constraints below. Its canonical comparison schema is removed with the encoding.
- **Para-stereo setting.** `para_stereo` is removed from `CanonicalizeContext`. The setting is read
  only by the structure partition, which no search uses once Structure and Full start from the
  constitution partition. The search itself separates groups that differ only in configuration, so
  molecules with para-stereo centres are canonicalized exactly without it.
  `StereoModel::para_stereo` of `umol-graph` governs stereo recognition and validation and is not
  affected.
- **Leaf.** Both aggregates use one leaf. The search orders the atoms. The key of a leaf is read
  from the entity order of the leaf, without remapping the aggregate below Full, and it is built
  when the leaf is first compared with another one. A `canonicalize` whose search ends with a
  single leaf builds no key. The molecule gives up its row leaf.
- **Order of non-atom entities.** In a canonical molecule and in a canonical span, the non-atom
  entities of each kind are ordered by their participants, at every level: bonds and noncovalent
  bonds by endpoint pair, dative bonds by donors and acceptor, aromatic systems and multicenter
  bonds by their atoms, stereo atoms by site, and stereo bonds by site bond. Entities of one kind on
  the same participants keep their search order. A span places the lhs-present entities of each kind
  before the added ones and orders each group this way. The order is a function of the atom
  numbering and does not depend on the level. Canonical molecules keep the order they have today;
  the canonical form of spans changes.
- **Topology.** At Topology both aggregates search the topology search input.
- **Canonical form.** The canonical form of an aggregate is the leaf of minimum key among the
  leaves of its search. It is the same for every renumbering of the entities, every participant
  frame, and every equivalent encoding of a value. It is not required to be the minimum over all
  renumberings, and the canonical form of a constrained aggregate is not required to extend that of
  the same aggregate without constraints.
- **Colors.** `ReframeKind` classifies a field or a constraint by its behavior under a
  participant-frame action as `Independent`, `Positional`, or `Relative`. The color of an entity is
  the frame-independent part of its normal form, as a borrowed view of the form. For atoms, bonds,
  dative bonds, and noncovalent bonds it is the whole form. For aromatic systems and multicenter
  bonds it leaves out the electron counts. For stereo entities it is the stereo kind and the
  stereogenicity constraint. A positional value enters the colors of the incidences of its
  participants. A frame-relative value enters no color. The colors of each kind are ranked once into
  dense integers, and the search uses the integers. The search partition is built from them at
  every level, so electron counts and dative roles split it.
- **Color trait.** `AsColor` has the associated type `Color<'a>`, bounded by `Ord` and `Hash`, and
  the method `as_color`. It is implemented on the eight forms and not on the views, and it is
  visible within the crate only. The colors are `AtomColor`, `BondColor`, `DativeBondColor`,
  `AromaticSystemColor`, `MulticenterBondColor`, `NoncovalentBondColor`, `StereoAtomColor`, and
  `StereoBondColor`, each beside its form, with derived equality, order, and hash. The trait states
  two properties. Invariance: the color of a normal form is equal before and after `reframe_by` for
  every compatible action. Completeness: a form without frame-relative parts is determined by its
  color and its positional parts, so the colors of an entity and of its incidences determine the
  entity below Structure. Pruning at Topology and Constitution relies on completeness. The
  frame-relative parts of a stereo form are outside its color; the key row and the generator check
  cover them.
- **Incidence colors.** The color of an incidence is its `Incidence`, as the incidence graph stores
  it, together with the positional part of the owning form for that participant. Only aromatic
  systems and multicenter bonds have such a part: the electron count of the atom, read by
  `ElectronCountsForm::at` with the `AtomPosition` of the atom, and `None` when the counts are
  undetermined. `Incidence` keeps its cases and takes its order from a derive, which equals the
  order written out today.
- **Span colors.** `AsColor` is implemented for `EntitySpan<T>` for every `T` that implements it.
  Each side gives its color, and a `Modified` entity whose two colors are equal gives `Unchanged`.
  The variants of `EntitySpan` are reordered to `Unchanged`, `Removed`, `Modified`, `Added`, so
  that the derived order places added entities last, as the lhs-anchored form requires. A comment
  at the enum states this reason. Nothing else reads the derived order of `EntitySpan`. The
  positional part of an incidence is lifted through the span in the same way.
- **Inline constraints.** They contribute as inherent fields do, to the color of their entity and to
  its key row. No inline constraint selects Full. The level is Topology with atoms and bonds only,
  Constitution with a dative bond, aromatic system, multicenter bond, or noncovalent bond, and
  Structure with a stereo entity.
- **Frame-relative values.** The configuration and the ligand-symmetry, fluxionality, and topicity
  constraints of a stereo entity enter its key row with the form transported into the frame of the
  leaf: the representative action of the renumbered ligands is applied to the form, which is then
  normalized. Pruning uses a generator only when one frame action carries the whole form of every
  stereo entity, on each side of a span, to the form of its image.
- **Full.** Molecule-level constraints select Full. Both aggregates search it from the constitution
  partition with pruning, and a generator is used only when it carries the complete aggregate onto
  itself. The check remaps the aggregate and compares it under framed equality. Where a
  molecule-level constraint breaks the symmetry of the aggregate, this costs more than the
  exhaustive search. Such constraints are not expected, and the cost is accepted. The key of a
  Full leaf is the key of Structure followed by the molecule-level constraints of the aggregate
  remapped and reframed in the entity order of the leaf, in their derived order.
- **Normalization.** An aggregate is normalized once, before its colors are read, at every level:
  in place in `canonicalize`, and on a copy in `canonical_eq` and `canonical_hash`. No value is
  normalized when a color or a key row is read, except a stereo form after its transport into the
  frame of a leaf. No wrapper type carries the normal form.
- **Kinds and queries.** `ReframeKind` sits beside `FrameTransport`. Each of the seven forms other
  than `AtomForm` gains its field-kind enum, and the eight field-kind enums derive `VariantArray`.
  The field-kind enums and the eight constraint-key enums have `reframe_kind`, and a property checks
  it against `reframe_by`.
- **Key rows.** A key row is the participants of the entity in their new numbers, in the
  representative order of the entity, followed by the integers of its colors. An atom gives its
  integer. A bond or a noncovalent bond gives its endpoint pair and its integer. A dative bond
  gives its sorted donors, its acceptor, and its integer. An aromatic system or a multicenter bond
  gives its sorted atoms, the integer of the incidence of each atom in that order, and its integer.
  A stereo entity gives its site, its ligands in the frame of the leaf, and its transported form.
  No form is read for a row except a stereo form. Keys are compared only within one run, and their
  stability across code changes is not a design consideration.
- **Equality and hash.** `canonical_eq` and `canonical_hash` compare and hash the key of the
  minimum leaf and never build the canonical form. The integers are ranks within one aggregate, so
  the key carries the distinct colors in rank order beside the rows.
- **Encoding.** `CanonicalKeyValue`, its encoding functions, the explicit positions, and the
  per-field key functions planned before the trial of colors are removed or not written. The order
  of every value is the derived order of its type.

The search of each level, for both aggregates:

| Level | Selected by | Searched graph | Initial partition | Enumeration | Generator check | Leaf key |
| --- | --- | --- | --- | --- | --- | --- |
| Topology | Atoms and bonds only | Topology search input | Entity colors | Pruned | None | Read from the entity order |
| Constitution | A dative bond, aromatic system, multicenter bond, or noncovalent bond | Adapter over the Constitution incidence graph | Colors | Pruned | None | Read from the entity order |
| Structure | A stereo entity | Adapter over the Full incidence graph | Colors | Pruned | The whole form of every stereo entity | Read from the entity order, with transported stereo forms |
| Full | A molecule-level constraint | Adapter over the Full incidence graph | Colors | Pruned | The whole aggregate | As Structure, with the molecule-level constraints of the remapped aggregate |

Every leaf key is deferred.

## Staged implementation plan

S0 to S3 were planned before the search study. S4 to S12 implement the decisions of the study, which
replace the Target section where the two differ. Their changes exist as an uncommitted trial in the
worktree `scratch/244-canonicalization-molecule-structure`. Each subitem moves the part of the trial
it names into the tree and adds its tests. The subitems from S6c on were replanned after the trial
of colors and implement its decisions; their trial is in the worktree `scratch/244-color-trial`.

Work proceeds one subitem at a time. After each subitem the narrow tests run and work stops for
review. After each subitem that changes search behavior, the reference comparison runs first.
Expected values of exact cases are shown, old beside new, before any is edited. Each subitem reports
its changed files and a suggested commit message; no git command is run.

After S0, molecule code is in `umol-graph-ir/src/ir/canonicalize/molecule.rs` and span and reaction
code in `umol-graph-ir/src/ir/canonicalize/reaction_span.rs`, each with its own `tests.rs`. A span
function takes the name of its molecule counterpart. From S1 on, a bare function or test name means
the one in `reaction_span` unless the subitem says otherwise. Tests follow the test-writing skill
and sit in source order.

S0 relocated about 1,450 lines of molecule code and about 650 lines of span code without changing
them, and split the 6,594-line tests file. S4 to S8 change the three canonicalization source files
by about 1,400 added and 630 removed lines, the size of the trial. They also remove about 950 lines
of functions that lose their last caller, with their tests. This estimate predates the trial of
colors and does not cover the subitems from S6c on.

### S0 — Baseline, reference, and module split

No behavior changes in this stage.

#### S0a — Benchmark cases **Done**

**Module:** `umol-graph-ir/benches/canonicalize.rs`. Additive (green). [dep: none]

Three benchmark functions add 90 cases. A modified span changes the charge of atom 0 and, in the
stereo cases, flips the configuration of the first stereo atom, or of the first stereo bond when
there is no stereo atom. Clippy on the benchmark target is clean, and every new case runs.

- `canonicalize/reaction_span/canonicalize` and `canonicalize/reaction_span/canonical_hash`, over
  the spans the existing `canonical_eq_equal` group uses: each corpus molecule plus one added
  nitrogen atom.
- `canonicalize/reaction_span_modified/canonicalize` and `.../canonical_eq_equal`: each corpus
  molecule with one modified atom and, for the stereo cases, one changed stereo configuration.
- `canonicalize/inline_constraint/molecule` and `.../reaction_span`, each with `canonicalize` and
  `canonical_eq_equal`: the three scaling cases with `#v1` on every hydrogen, which run at Full.
- `canonicalize/complete/canonical_eq_equal` and `canonicalize/complete/canonical_hash`: molecules
  of the complete corpus, which covers the stereo cases.

#### S0b — Reference checkout **Done**

**Location:** `scratch/244-canonicalization-reaction-span`, not committed. [dep: S0a]

Baseline at commit 5049b9796: Criterion medians with a 1-second warm-up and 3 seconds of
measurement, one run, single-threaded, Nauty, Apple M1 Pro, load average about 3. All 155 medians
are in `scratch/244-compare/baseline/reference.tsv`, which is not committed. One median, the
modified span of `cis_trans_stereo_bond`, is from a second run in the reference checkout. Each span
adds one nitrogen atom to the molecule; a modified span is described in S0a.

| Case | Molecule `canonicalize` | Molecule `canonical_eq` | Span `canonicalize` | Span `canonical_eq` | Modified span `canonicalize` |
| --- | ---: | ---: | ---: | ---: | ---: |
| `ordinary_naphthalene` | 21.3 µs | 36.0 µs | 72.6 µs | 134.4 µs | 29.5 µs |
| `disconnected_rings` | 34.4 µs | 60.0 µs | 4.43 ms | 8.76 ms | 406.1 µs |
| `overlay_heavy` | 91.8 µs | 166.7 µs | 57.2 µs | 93.0 µs | 56.6 µs |
| `tetrahedral_stereo` | 41.9 µs | 72.4 µs | 24.8 µs | 36.3 µs | 25.4 µs |
| `cis_trans_stereo_bond` | 34.7 µs | 58.1 µs | 28.1 µs | 42.7 µs | 29.2 µs |
| `mixed_atom_and_bond_stereo` | 72.3 µs | 134.7 µs | 48.0 µs | 73.8 µs | 47.7 µs |
| `frame_relative_stereo_constraint` | 46.7 µs | 86.4 µs | 26.0 µs | 43.6 µs | 27.7 µs |
| `meso_dichlorobutane` | 68.7 µs | 127.4 µs | 46.7 µs | 79.8 µs | 33.6 µs |
| `para_stereo_trichloropentane` | 129.5 µs | 243.0 µs | 62.1 µs | 105.1 µs | 44.0 µs |
| `para_stereo_cascade` | 415.9 µs | 782.8 µs | 520.2 µs | 1.01 ms | 295.3 µs |
| `feature_free_connected` | 37.8 µs | 72.6 µs | 2.19 ms | 4.44 ms | 375.2 µs |
| `feature_free_disconnected` | 40.4 µs | 75.6 µs | 2.83 ms | 5.66 ms | 708.3 µs |
| `symmetry_heavy_radicals` | 31.1 µs | 58.1 µs | 11.12 ms | 22.11 ms | 1.81 ms |

With `#v1` on every hydrogen, which puts the case at Full:

| Case | Molecule `canonicalize` | Molecule `canonical_eq` | Span `canonicalize` | Span `canonical_eq` |
| --- | ---: | ---: | ---: | ---: |
| `feature_free_connected` | 2.16 ms | 4.30 ms | 2.62 ms | 5.45 ms |
| `feature_free_disconnected` | 2.78 ms | 5.67 ms | 3.52 ms | 6.96 ms |
| `symmetry_heavy_radicals` | 11.32 ms | 22.52 ms | 14.22 ms | 28.49 ms |

The present span `canonicalize` is slower than the molecule's on six cases: `ordinary_naphthalene`,
`disconnected_rings`, `para_stereo_cascade`, and the three feature-free cases. It is faster on the
other seven, although each span has one atom more.

The user-level cargo configuration sends every checkout to one target directory. The reference
checkout gets its own target directory, `/Users/dr/.cargo-target-244`, set in its
`.cargo/config.toml`. The
S0a benchmark file is copied into the reference checkout; it uses public API only. The benchmark
groups of S0a and the existing reaction and reaction-span groups run in the reference checkout, and
the numbers are recorded here as the baseline.

#### S0c — Reference pairs **Done**

**Location:** the reference checkout and `scratch/244-compare`, not committed. [dep: S0b]

A temporary test beside `reaction_span_canonical_pair_strategy` in the reference checkout draws
10,000 pairs with a fixed seed. It writes one pair per line as DSL text, with the reference
`canonical_eq` answer for both `para_stereo` settings. A span whose text does not parse back to an
equal span is counted and reported. Hand-written pairs are appended as text. `scratch/244-compare`
is a standalone program with its own `[workspace]` table and path dependencies on the working tree.
It reads the file and compares the public `canonical_eq` answers.

**Check:** run against the unchanged working tree, the program reports no mismatch.

`scratch/244-compare/pairs.tsv` holds 10,000 generated pairs and 28 handwritten ones. Identical
pairs are skipped when drawing, because `canonical_eq` answers them without a search: 26,541 draws
gave 16,541 identical pairs. Of the generated pairs 4,276 are equal and 5,724 unequal; of the
handwritten ones 18 and 10. Every span parses back to an equal span, a second run reproduces the
file, and the two `para_stereo` settings agree on every pair in the reference. Against the unchanged
working tree the program reports no mismatch, and it does report a deliberately flipped answer.

#### S0d — Molecule submodule **Done**

Move every function that takes a `Molecule`, together with `stereo_refinement_descriptor`, the
`impl Molecule` block, and `impl Canonicalize for Molecule`, into `canonicalize/molecule.rs`.
`molecule_canonicalize_level` and `molecule_counts`, which already carry the aggregate in their
names, stay in `canonicalize.rs`. Move each test that uses a moved function or a private `Molecule`
method into `canonicalize/molecule/tests.rs`. Fixtures and test support used by more than one tests
file stay in `canonicalize/tests.rs` and are made visible to the others. No function body and no
test body changes. Breaking while the move is in progress (red→green). [dep: S0c]

**Check:** the test names and their count are the same before and after, and every test passes. The
reference comparison reports no mismatch.

Twenty functions and three impl blocks, 1,540 lines, moved to `canonicalize/molecule.rs` as whole
items. The module imports what it uses explicitly. 54 tests, 4 fixtures, and 10 support functions
moved to `canonicalize/molecule/tests.rs`; 53 tests and 5 fixtures stay. Nine support functions and
the test-only `AutomorphismAdapter::automorphisms` are used from both tests files; they stay in
`canonicalize/tests.rs` as `pub(super)`. Every line of the old source file outside its import block
is present unchanged in one of the two source files, and the tests files differ from the old one
only in those ten signatures. The crate lists the same 11,967 library tests, 345 of them for
canonicalization, and all pass, as do the integration tests and the canonicalization properties.
Clippy over all targets and features is clean. The reference comparison reports no mismatch.

#### S0e — Reaction-span submodule **Done**

Move every function that takes a `ReactionSpan`, the `impl ReactionSpan` and `impl Reaction` blocks,
and their `Canonicalize` impls into `canonicalize/reaction_span.rs`.
`reaction_span_canonicalize_level`, `reaction_canonicalize_level`, and `reaction_span_counts` stay
in `canonicalize.rs`. The moved functions drop the prefix: `entity_color_key`,
`canonical_candidate`, `comparison_key`, `lhs_present`. `reaction_span_entity_keys` becomes
`initial_color_keys`, the name of its molecule counterpart. Move the span and reaction tests that
use a moved function or a private method into `canonicalize/reaction_span/tests.rs`;
`test_reaction_span_entity_keys_incidence` becomes `test_initial_color_keys_incidence`. Nothing else
in a function body or a test body changes. Breaking while the move is in progress (red→green).
[dep: S0d]

**Check:** as S0d.

Six functions and four impl blocks, 815 lines, moved to `canonicalize/reaction_span.rs`. Beyond the
five renames, four inline paths such as `super::dative::DativeBondForm` gained one `super::`, and
rustfmt rewrapped two lines whose length changed. Five tests moved to
`canonicalize/reaction_span/tests.rs`: the incidence-key test, the two description-level tests, and
the two tests of the span's private methods. Span and reaction tests that go through the public
trait only stay in `canonicalize/tests.rs`. The checks of S0d hold: the same 11,967 library tests
with one renamed, all passing; the integration tests and canonicalization properties pass; clippy
is clean; the reference comparison reports no mismatch.

#### S0f — Topology search input **Done**

Rename `CompactTopologyCarrier` to `TopologySearchInput`. Its constructor becomes the function
`topology_search_input` in `canonicalize::molecule`, and `test_compact_topology_carrier_new` becomes
`test_topology_search_input`. Breaking while the rename is in progress (red→green). [dep: S0d]

**Check:** as S0d.

The type, its constructor, its doc comments, one `unreachable!` message, and the test with its local
variable are renamed. Ignoring whitespace, the constructor body changes only in that message and in
the type name of the returned value. The checks of S0d hold, with one test renamed.

### S1 — Span tag order

#### S1a — Tag positions **Done**

`SpanTagPosition`: `REMOVED = 1`, `MODIFIED = 2`, `ADDED = 3`. The span tag sentence of the
data-type guide changes with it. Breaking for span keys and representatives (red→green).
[dep: S0e]

**Tests:** the rows of `test_span_tag_position`. The exact cases of
`tests/canonicalization/reaction_span.rs`, the span cases of the unit tests, and the Python
canonicalization tests of `test_reaction_span.py` and `test_reaction.py` are rerun; any changed
expected value is shown before it is edited.

The three positions and the rows of `test_span_tag_position` changed, and the guide states the new
positions. No other expected value changed: the library tests, the exact cases, the canonicalization
properties, and the six Python canonicalization tests pass as they are. The reference comparison
reports no mismatch.

#### S1b — Remapping from the search order **Done**

Add `remapping_from_order`, the counterpart of the molecule function of that name.
`ReactionSpan::tracked_canonicalize_at` and the leaf of `canonical_candidate` use it. Delete
`lhs_anchored_remapping_from_order` and `lhs_present`. No canonical output changes relative to S1a
(green). [dep: S1a]

**Tests:** `test_remapping_from_order`, with the rows of the molecule test of that name.
`test_reaction_span_tracked_canonicalize_at`: spans with added entities of each kind stored before
lhs-present ones, with the exact lhs-anchored result at each level.

The two tests, 14 cases, were written first and pass both with the reordering step and without it.
The canonical form and the tracked remapping of all 20,056 spans of the pair file, under both
`para_stereo` settings, are identical before and after the change. The stage gate passes, including
the two Python test files.

### S2 — Orbit pruning at Topology and Constitution

#### S2a — Search options as a parameter **Done**

`canonical_candidate` takes `options: CanonicalSearchOptions` and returns the
`CanonicalSearchResult`, as the molecule's per-level functions take their options. The two span
methods pass the options of today: no pruning, node-id order. No behavior change (green).
[dep: S1b]

Each of the two span methods builds the options itself. The canonical forms and tracked remappings
of all spans of the pair file are identical before and after, and the tests, clippy, and the
reference comparison pass.

#### S2b — Pruning on **Done**

The two span methods pass pruning with backend canonical order at Topology and Constitution. The
search graph, the partition, and the leaf stay as they are. [dep: S2a]

**Tests:** `test_canonical_candidate`: for each row, the candidate with pruning equals the candidate
without it in key, and the visited-leaf counts of both are exact. Rows: three
equivalent unchanged atoms; three equivalent added atoms; three equivalent removed atoms; two
unchanged atoms and one modified atom; equivalent atoms under an unchanged and under a modified
dative bond and aromatic system at Constitution.

`test_canonical_candidate`, eight cases, passes its options itself and was written before the
switch. Seven rows visit one leaf with pruning, against 6 or 2 without it. In the modified aromatic
system row the electron counts single out one of three atoms, and the search visits 6 leaves without
pruning and 3 with it: positional electron counts are automorphism colors and do not split the
constitution partition, so the three atoms share a cell and two of them share an orbit.

The canonical forms and tracked remappings of all spans of the pair file are identical to S2a. The
reference comparison reports no mismatch, and the stage gate passes, including the two Python test
files.

Medians of the five cases that run at Topology or Constitution, against the S0b baseline:

| Case | Span `canonicalize` | Span `canonical_eq` | Modified span `canonicalize` |
| --- | ---: | ---: | ---: |
| `ordinary_naphthalene` | 72.6 → 37.8 µs | 134.4 → 62.8 µs | 29.5 → 29.6 µs |
| `disconnected_rings` | 4.43 ms → 57.3 µs | 8.76 ms → 101.3 µs | 406.1 → 49.4 µs |
| `feature_free_connected` | 2.19 ms → 60.4 µs | 4.44 ms → 113.9 µs | 375.2 → 52.3 µs |
| `feature_free_disconnected` | 2.83 ms → 61.2 µs | 5.66 ms → 114.1 µs | 708.3 → 52.8 µs |
| `symmetry_heavy_radicals` | 11.12 ms → 40.1 µs | 22.11 ms → 71.9 µs | 1.81 ms → 33.0 µs |

Span `canonical_hash` follows `canonicalize`. The eight cases with stereo entities are within 4% of
the baseline in every span group. The span `canonicalize` of the five cases is now 1.3 to 1.8 times
the molecule's.

### S3 — Stereo-sensitive search at Structure

#### S3a — Span stereo descriptor **Done**

Add `stereo_refinement_descriptor`, the counterpart of the molecule function of that name. It takes
the site class, the ligand classes, and the entity span of the stereo form, and minimizes over the
frame permutations of the kind with both configurations transported. Additive (green). [dep: S0e]

**Tests:** `test_stereo_refinement_descriptor_frame_invariant` for every stereo kind and for
unchanged, added, removed, and modified entities, including one undetermined side;
`test_stereo_refinement_descriptor_error` for a frame of the wrong degree;
`test_stereo_refinement_descriptor` with exact values.

The function is generic over the stereo form and takes a closure that reads the configuration, so
stereo atoms and stereo bonds share it. Each side has the row shape of the molecule descriptor, and
the two sides sit under the normalized span tag. Until S3d calls the function it carries
`#[cfg_attr(not(test), expect(dead_code))]`, which clippy rejects once the function is used.

`test_stereo_refinement_descriptor`, ten cases, was added to the two planned tests. Frame invariance
also holds for a descriptor that minimizes each side separately or leaves out the rhs, so it does
not fix the value. The row with two ligands of one class and exchanged configurations on the two
sides requires one permutation to act on both sides. The frame-invariance test has 36 cases, six
span shapes for each of the six kinds, and the error test has six. Clippy, the library tests, and
the two integration tests pass.

#### S3b — Span structure partition **Done**

Add `constitution_entity_classes`, `structure_partition_descriptors`,
`para_stereo_partition_descriptors`, and `structure_partition`, the counterparts of the four
molecule functions. Additive (green). [dep: S3a]

**Tests:** `test_structure_partition`: round counts with and without `para_stereo` for the span of
the molecule para-stereo fixture and for a span whose rhs configuration alone splits a cell; one
round for a span without stereo entities.

The four functions follow the molecule functions line by line, with the span accessors and the span
stereo descriptor in place of the molecule ones. They carry the unused-function attribute of S3a
until S3d. The fixture `para_stereo_canonicalization_molecule` moved unchanged from the molecule
tests to the shared test support, so that both test modules use it.

`test_structure_partition`, four cases, asserts the round count and the cells of the stereo atoms,
as a set of cells, without and with `para_stereo`:

| Span | Without `para_stereo` | With `para_stereo` |
| --- | --- | --- |
| Identity span of the para-stereo fixture | 1 round; the two inner stereo atoms share a cell, the eight outer ones form four pairs | 2 rounds; ten single cells |
| Two equal fragments, both stereo atoms modified from coset 0 to coset 1 | 1 round; one cell of two | 2 rounds; one cell of two |
| The same, with the rhs of the second stereo atom undetermined | 1 round; two single cells | 1 round; two single cells |
| No stereo entities | 1 round | 1 round |

The fixture rows give the round counts of the molecule test. Clippy, the library tests, and the two
integration tests pass.

#### S3c — Span stereo generator check **Done**

Add `generator_preserves_stereo` and `retain_stereo_preserving_generators`, the counterparts of the
two molecule functions. One frame action must carry the configuration of each present side to that
of the image. Additive (green). [dep: S0e]

**Tests:** `test_retain_stereo_preserving_generators`: the rows of the molecule table as unchanged
spans; a generator that preserves the lhs configuration and not the rhs one; a generator between a
modified and an unchanged stereo entity; a stereo-bond endpoint exchange.

The two functions follow the molecule functions. The configuration check takes the lhs and the rhs
of the entity and of its image in turn, with the one frame action derived from the ligands. A side
present on the entity and absent on the image, or the reverse, fails the check. The functions carry
the unused-function attribute of S3a until S3d.

`test_retain_stereo_preserving_generators` has nine cases: the five molecule rows as unchanged
spans, and four span rows.

| Span row | Generator | Retained |
| --- | --- | --- |
| Stereocentre with two hydrogens, lhs coset 0, rhs undetermined | Hydrogen exchange | No |
| The same, lhs undetermined, rhs coset 0 | Hydrogen exchange | No |
| Two equal fragments, one stereo atom unchanged, the other modified from coset 0 to coset 1 | Fragment exchange | No |
| Stereo bond modified from coset 0 to coset 1 | Endpoint exchange | Yes |

The first two rows together require both sides to be checked. Clippy, the library tests, and the
two integration tests pass.

#### S3d — Structure search **Superseded by S4b**

The Structure arm of `canonical_candidate` starts from `structure_partition` and runs
`canonical_search_with_generator_filter` with the span generator check. The span methods pass
pruning with backend canonical order at Structure. Breaking for representatives of spans with stereo
entities (red→green). [dep: S2b, S3b, S3c]

**Tests:** Structure rows of `test_canonical_candidate`: three stereocentres with equal
configurations; one with a different lhs configuration; one with a different rhs configuration
only. Stereo rows of `test_reaction_span_tracked_canonicalize_at` under both `para_stereo`
settings.

The Structure arm is wired as the molecule's `canonicalize_structure`, and the unused-function
attributes of S3a to S3c are removed. No expected value of an existing test changed: the library
tests, the exact cases, the canonicalization properties, and the two Python test files pass as they
are. The reference comparison reports no mismatch.

Of the 20,056 spans of the pair file, two changed their canonical form, the two spans of one equal
pair, under both `para_stereo` settings. Two removed bonds exchange their positions and the site of
the stereo bond follows; the atoms keep their numbers. All other canonical forms and tracked
remappings are identical to S2b.

`test_canonical_candidate` has four Structure rows, with the visited leaves without and with
pruning:

| Row | Without pruning | With pruning |
| --- | ---: | ---: |
| Three equal fragments with equal stereocentres | 6 | 1 |
| The same, lhs configuration of one differs | 2 | 1 |
| The same, rhs configuration of one differs | 2 | 1 |
| Stereocentre with two hydrogens, lhs undetermined, rhs coset 0 | 2 | 2 |

In the last row the hydrogen exchange is an automorphism of the colored graph that does not preserve
the rhs configuration, so the generator check removes it and both leaves are visited.
`test_reaction_span_tracked_canonicalize_at` runs all eight rows under both `para_stereo` settings.

Medians of the six cases that run at Structure, against the S0b baseline:

| Case | Span `canonicalize` | Molecule `canonicalize` | Span `canonical_eq` | Modified span `canonicalize` |
| --- | ---: | ---: | ---: | ---: |
| `tetrahedral_stereo` | 24.8 → 46.8 µs | 41.9 µs | 36.3 → 80.8 µs | 25.4 → 62.5 µs |
| `cis_trans_stereo_bond` | 28.1 → 42.5 µs | 34.7 µs | 42.7 → 72.3 µs | 29.2 → 46.5 µs |
| `mixed_atom_and_bond_stereo` | 48.0 → 83.9 µs | 72.3 µs | 73.8 → 151.8 µs | 47.7 → 98.6 µs |
| `meso_dichlorobutane` | 46.7 → 72.7 µs | 68.7 µs | 79.8 → 128.7 µs | 33.6 → 87.0 µs |
| `para_stereo_trichloropentane` | 62.1 → 124.5 µs | 129.5 µs | 105.1 → 233.9 µs | 44.0 → 115.8 µs |
| `para_stereo_cascade` | 520.2 → 391.0 µs | 415.9 µs | 1.01 ms → 755.3 µs | 295.3 → 327.3 µs |

Span `canonicalize` is 1.5 to 2.0 times slower than the baseline on five cases and 25% faster on
`para_stereo_cascade`. The modified spans are 1.1 to 2.6 times slower. Span `canonicalize` is now
0.94 to 1.22 times the molecule's on the six cases. The cases at Topology, Constitution, and Full
are within 4% of S2b.

This subitem is not committed. The [search study](#search-study) replaces it: S4b wires the
Structure search from the constitution partition and removes the functions of S3a and S3b.

### S4 — Structure search

#### S4a — Symmetric stereo benchmark cases **Done**

**Module:** `umol-graph-ir/benches/canonicalize.rs`. Additive (green). [dep: S0a]

The corpus gains the three cases of the study: `symmetric_stereo_triisopropanolamine`,
`symmetric_stereo_butanol`, and `symmetric_stereo_dimethylbutanol`.

Clippy on the benchmark target is clean, and the three cases run.

#### S4b — Span Structure search **Done**

The Structure arm of `canonical_candidate` starts from the constitution descriptors and runs
`canonical_search_with_generator_filter` with `retain_stereo_preserving_generators`. The span
methods pass pruning with backend canonical order at Structure. This replaces the uncommitted wiring
of S3d. Delete the span's `stereo_refinement_descriptor`, `constitution_entity_classes`,
`structure_partition_descriptors`, `para_stereo_partition_descriptors`, and `structure_partition`,
with their tests. Canonical forms are unchanged (green). [dep: S2b, S3c]

**Tests:** Structure rows of `test_canonical_candidate`, with exact visited-leaf counts without and
with pruning: three equal fragments with equal stereocentres; the same with a different lhs
configuration on one; the same with a different rhs configuration on one; a stereocentre with two
hydrogens, lhs undetermined and rhs coset 0. Stereo rows of
`test_reaction_span_tracked_canonicalize_at`.

The Structure arm is wired as stated, and both span methods prune at every level below Full. The
five functions of S3a and S3b are removed with their four tests and the imports only they used, and
the unused-function attributes of S3c are removed.

The canonical forms and tracked remappings of all spans of the pair file are identical to S2b, and
the reference comparison reports no mismatch. The library tests, the two integration tests, and the
canonicalization properties pass, and clippy is clean.

`test_canonical_candidate` has four Structure rows, with the visited leaves without and with
pruning:

| Row | Without pruning | With pruning |
| --- | ---: | ---: |
| Three equal fragments with equal stereocentres | 6 | 1 |
| The same, lhs configuration of one differs | 6 | 3 |
| The same, rhs configuration of one differs | 6 | 3 |
| Stereocentre with two hydrogens, lhs undetermined, rhs coset 0 | 2 | 2 |

In the second and third rows the three stereo atoms share a color, which holds the kind and not the
configuration, so the three fragments share their cells. The exchange of the two equal fragments is
the one generator the check keeps, and the pruned search visits three leaves. The stereo rows of
`test_reaction_span_tracked_canonicalize_at` pass unchanged.

Medians of the nine cases that run at Structure, against the committed search:

| Case | Span `canonicalize` | Span `canonical_eq` |
| --- | ---: | ---: |
| `tetrahedral_stereo` | 24.2 → 23.4 µs | 36.4 → 35.7 µs |
| `cis_trans_stereo_bond` | 27.4 → 26.7 µs | 42.4 → 41.4 µs |
| `mixed_atom_and_bond_stereo` | 46.5 → 45.3 µs | 73.3 → 72.6 µs |
| `meso_dichlorobutane` | 46.7 → 48.9 µs | 78.8 → 83.9 µs |
| `para_stereo_trichloropentane` | 62.3 → 64.9 µs | 105.5 → 112.9 µs |
| `para_stereo_cascade` | 523.8 → 561.5 µs | 1.01 → 1.12 ms |
| `symmetric_stereo_triisopropanolamine` | 171.9 → 74.8 µs | 326.9 → 133.1 µs |
| `symmetric_stereo_butanol` | 1.40 ms → 96.7 µs | 2.82 ms → 179.1 µs |
| `symmetric_stereo_dimethylbutanol` | 194 ms → 206.8 µs | 382 ms → 400.8 µs |

The three cases whose search ends with one leaf are within 3% of the committed search.
`meso_dichlorobutane` and `para_stereo_trichloropentane` are 4% to 7% slower and
`para_stereo_cascade` 7% to 11% slower, from backend calls that prune nothing, as the study
recorded. The three symmetric cases take 75 to 207 µs for `canonicalize`.

#### S4b.1 — Contradictory stereo configurations **Done**

The color of a stereo entity reads the kind of the normalized configuration, for stereo atoms and
stereo bonds of both aggregates. A contradictory configuration is then reported where the colors are
built. Until this subitem the colors did not normalize the configuration: `canonical_hash` and
`canonical_eq` of a reaction span panicked on a contradictory one at the first leaf, and the
molecule searches reported it from the structure partition, which S4c takes out of the Structure
search. Colors of satisfiable inputs are unchanged (green). [dep: S4b]

**Tests:** `test_initial_colors_error` of the molecule, as a table with the charge row, a stereo
atom, and a stereo bond. `test_initial_color_keys_error` of the span: an unchanged stereo atom, a
modified one with a contradictory rhs, and a stereo bond.
`test_reaction_span_canonical_hash_error`.

The library tests, the two integration tests, and the canonicalization properties pass. The
reference comparison reports no mismatch. The canonical forms and tracked remappings of the pair
file and the molecule answers are identical to S4b.

#### S4c — Molecule Structure search **Done**

`canonicalize_structure` and the Structure branch of `Molecule::canonical_key_at` start from the
constitution descriptors and run `canonical_search_with_generator_filter` with
`retain_stereo_preserving_generators`, with pruning and backend canonical order. The data-type
guide's prescription of the stereo-sensitive refinement is replaced. Breaking for representatives
of molecules with stereo entities (red→green). [dep: S0d, S4b.1]

**Tests:** `test_molecule_canonical_key_at`: for the molecule stereo fixtures and the symmetric
stereo case with six automorphisms, the Structure key equals the exhaustive minimum over the
constitution partition, with exact visited-leaf counts. The existing Structure tests are rerun; any
changed expected value is shown first.

Both Structure searches of the molecule start from the constitution descriptors and prune with the
stereo generator check. The guide states this rule and limits the stereo-sensitive refinement to the
complete search of a molecule with constraints, which keeps it until S7b.

No expected value of an existing test changed. `canonical_eq` and the equality of `canonical_hash`
of molecules give the same answer as the committed search on the 20,056 pairs, and the reference
comparison of spans reports no mismatch. The library tests, the two integration tests, and the
canonicalization properties pass, and clippy is clean.

`test_molecule_canonical_key_at` has four rows. In each the Structure key equals the key of the
exhaustive search from the constitution partition, and the pruned search selects the same key. The
visited leaves without and with pruning, beside those of the structure partition:

| Fixture | Constitution partition | Structure partition |
| --- | ---: | ---: |
| Stereo atom | 1, 1 | 1, 1 |
| Stereo bond | 1, 1 | 1, 1 |
| Carbon with four equal carbon ligands and a tetrahedral configuration | 24, 24 | 24, 24 |
| Para-stereo fixture | 8, 8 | 2, 2 without `para_stereo`; 1, 1 with it |

On the third fixture every generator the backend returns exchanges two ligands and fails the stereo
check, so pruning removes nothing under either partition.

Medians of the nine cases that run at Structure, against the committed search:

| Case | `canonicalize` | `canonical_hash` | `canonical_eq` |
| --- | ---: | ---: | ---: |
| `tetrahedral_stereo` | 40.6 → 21.5 µs | 39.0 → 18.5 µs | 73.3 → 34.6 µs |
| `cis_trans_stereo_bond` | 33.1 → 20.9 µs | 31.7 → 18.0 µs | 57.8 → 33.2 µs |
| `mixed_atom_and_bond_stereo` | 72.2 → 40.6 µs | 71.5 → 37.6 µs | 134.1 → 69.5 µs |
| `meso_dichlorobutane` | 68.0 → 57.3 µs | 65.9 → 52.8 µs | 127.1 → 102.7 µs |
| `para_stereo_trichloropentane` | 130.6 → 81.8 µs | 124.4 → 77.7 µs | 241.8 → 153.2 µs |
| `para_stereo_cascade` | 415.4 → 862.5 µs | 401.1 → 841.6 µs | 780.2 µs → 1.76 ms |
| `symmetric_stereo_triisopropanolamine` | 132.1 → 77.4 µs | 270.1 → 72.5 µs | 529.2 → 141.1 µs |
| `symmetric_stereo_butanol` | 117.9 → 87.6 µs | 1.29 ms → 84.9 µs | 2.56 ms → 166.0 µs |
| `symmetric_stereo_dimethylbutanol` | 233.2 → 195.8 µs | 166 ms → 193.2 µs | 329 ms → 390.9 µs |

Eight cases are faster, by factors of 1.2 to 1.9 for `canonicalize`. `para_stereo_cascade` takes 2.1
to 2.2 times as long, the cost the study recorded for this search with the row leaf; S5d replaces
the leaf.

### S5 — One leaf for both aggregates

#### S5a — Deferred key **Done**

**Module:** `canonicalize.rs`. Additive (green). [dep: none]

Add `DeferredKey`. It holds the entity order of a leaf and builds the comparison key when the leaf
is first compared. `into_key` builds the key if no comparison did, and `into_built` returns it only
if one did.

**Tests:** `test_deferred_key_into_key`, `test_deferred_key_into_built`, and
`test_deferred_key_cmp`, with a counting build function: a key is built once, on the first
comparison or on `into_key`, and not at all for a leaf that is never compared.

`DeferredKey` stands before `TopologySearchInput` and compares by its built key. It carries the
unused-code attribute until S5b uses it. The three tests have seven cases; the library tests pass
and clippy is clean.

#### S5b — Span entity order and deferred leaf **Done**

Add `entity_order`: the atoms in search order, then the entities of each further kind with the
lhs-present ones first and each group by its participants, ties in search order. The leaf of
`canonical_candidate` takes its order from it and holds a `DeferredKey`. `canonical_candidate` takes
`defer_key`; `tracked_canonicalize_at` passes it and builds no key for a single leaf. Breaking for
span representatives (red→green). [dep: S4b, S5a]

**Tests:** `test_entity_order`: for each kind, added entities stored before lhs-present ones, and
two entities on the same participants. `test_canonical_candidate` runs each row with and without
`defer_key`. Three exact cases change and are shown first: the stereo-bond row of
`test_reaction_span_tracked_canonicalize_at`, and `case_3_stereo` and `case_1_lifecycle` of
`tests/canonicalization/reaction_span.rs`. The two Python test files are rerun.

`entity_order` and the deferred leaf are in place, and `DeferredKey` no longer carries the
unused-code attribute. `canonical_candidate` returns the key when it is requested, or when a
deferred search compared the selected leaf.

`test_entity_order` has nine rows under a reversed search order: bonds by endpoint pair; for each of
the seven non-atom kinds an added entity stored before an lhs-present one, with participants that
would place it first; and a removed and an added bond on one atom pair. Integrity admits no two
entities of a kind on the same participants on one side, so the last row is the case in which equal
participants occur. `test_canonical_candidate` runs every row with the key requested and deferred: a
deferred search returns the key only when it visited more than one leaf, and the same entity order.

The three exact cases changed as named. In each the atoms keep their numbers:

| Case | Bonds before | Bonds after |
| --- | --- | --- |
| Stereo-bond row of `test_reaction_span_tracked_canonicalize_at` | `0-2, 1-3, 4-6, 5-7, 0=1, 4=5`; stereo bonds on bonds 4 and 5 | `0=1, 0-2, 1-3, 4=5, 4-6, 5-7`; stereo bonds on bonds 0 and 3 |
| `case_3_stereo` | `[0 3] [1 3] [0 4] [1 4] [0 5] [2 5] [0 6] [2 6] [1 2]`; stereo bond on bond 8 | `[0 3] [0 4] [0 5] [0 6] [1 2] [1 3] [1 4] [2 5] [2 6]`; stereo bond on bond 4 |
| `case_1_lifecycle` | removed `[1 2]`, modified `[0 2]`, added `[0 3]` | modified `[0 2]`, removed `[1 2]`, added `[0 3]` |

The reference comparison reports no mismatch, and the molecule answers are identical. Of the 40,112
canonical span forms of the pair file, 7,332 change; all are identical to those of the trial. The
library tests, the two integration tests, the canonicalization properties, and the two Python test
files pass, and clippy is clean.

Against S4b, span `canonicalize` is 10% to 30% faster on 13 of the 16 corpus cases, whose search
ends with one leaf and builds no key, and unchanged on `meso_dichlorobutane`,
`para_stereo_trichloropentane`, and `para_stereo_cascade`. Span `canonical_eq` is 1% to 5% slower
on all 16 cases, the cost of the order by participants.

#### S5c — Molecule leaf functions **Done**

Add `entity_order`, `comparison_key`, and `leaf_key` in `canonicalize::molecule`. `entity_order`
gives the order of every entity of a level for an atom order, with the non-atom entities by their
participants. `comparison_key` reads the key of a molecule with every row in stored order.
`leaf_key` remaps and reframes the molecule in an entity order and reads its key. Additive (green).
[dep: S0d]

**Tests:** `test_entity_order`; `test_comparison_key`, with the exact keys of the present key tests
at each level; `test_leaf_key`, which equals the key of the present candidate functions on their
rows.

The three functions are added with the unused-code attribute, which S5d removes. `entity_order` is
that of the trial. `comparison_key` builds every row in the layout of the canonical comparison
schema of the data-type guide, which is the layout of the present candidate functions. The key of
the trial nested the same components differently; its layout is not taken over, so the keys and
hashes of molecules do not change in S5d.

`test_entity_order` has eight rows under a reversed search order: bonds; for each of the six
further kinds two entities stored against the order of their participants; and a row at
Constitution, which leaves a stereo atom out. `test_comparison_key` has seven rows with the exact
rows of bonds, of each of the four non-stereo overlay kinds, and of a stereo atom in its stored
frame, and one row in which the stereo block is absent below Structure. `test_leaf_key` has six
rows at the four levels: for two search orders the entity order and the key equal those of the
present candidate function of the level. The library tests pass and clippy is clean.

#### S5d — Molecule searches on the leaf **Done**

`canonicalize_topology`, `canonicalize_constitution`, `canonicalize_structure`, `canonicalize_full`,
and `Molecule::canonical_key_at` take the order of a leaf from `entity_order` and hold a
`DeferredKey` over `leaf_key`. Delete `topology_bond_key_rows`, `topology_candidate`,
`constitution_candidate`, `structure_candidate`, `complete_candidate`,
`canonical_kinded_stereo_frame`, and `CanonicalStereoFrame`, with their tests; rows that assert a
key move to `test_comparison_key`. Canonical molecules are unchanged (red→green).
[dep: S4c, S5a, S5c]

**Tests:** the tests of the four per-level functions and of `canonical_key_at` are rerun without a
changed expected value.

The four per-level functions and `canonical_key_at` take the order of a leaf from `entity_order` and
hold a `DeferredKey` over `leaf_key`. The five row-leaf functions, `canonical_kinded_stereo_frame`,
and `CanonicalStereoFrame` are removed, about 470 lines.

No expected value changed. In the molecule tests the three key wrappers over the removed functions
give way to one support function over `entity_order` and `leaf_key`, and 32 calls go through it;
the exact keys and canonical molecules these tests assert are as before. The present key tests stay
in place and are not moved into `test_comparison_key`: they assert the key of a search order, which
`comparison_key` does not take. `test_leaf_key` compares the key of a leaf with the comparison key
of the molecule remapped and reframed in the entity order.
`test_canonical_kinded_stereo_frame_error` is removed with its function.

The library tests, the two integration tests, the canonicalization properties, and the Python suite
pass, and clippy is clean. The reference comparison reports no mismatch, the molecule answers are
identical, and the canonical forms and tracked remappings of the spans are identical to S5b.

Medians of the molecule, against S4c for the nine Structure cases and against the S0b baseline for
the other seven, whose search had not changed since:

| Case | `canonicalize` | `canonical_hash` | `canonical_eq` |
| --- | ---: | ---: | ---: |
| `ordinary_naphthalene` | 21.3 → 17.2 µs | 21.0 → 24.0 µs | 36.0 → 43.4 µs |
| `disconnected_rings` | 34.4 → 29.1 µs | 33.3 → 37.6 µs | 60.0 → 69.2 µs |
| `feature_free_connected` | 37.8 → 33.4 µs | 38.6 → 43.3 µs | 72.6 → 81.6 µs |
| `feature_free_disconnected` | 40.4 → 36.7 µs | 41.4 → 45.6 µs | 75.6 → 84.9 µs |
| `symmetry_heavy_radicals` | 31.1 → 29.6 µs | 30.8 → 34.6 µs | 58.1 → 65.0 µs |
| `overlay_heavy` | 91.8 → 68.5 µs | 88.1 → 81.3 µs | 166.7 → 155.7 µs |
| `frame_relative_stereo_constraint` | 46.7 → 34.9 µs | 44.5 → 40.1 µs | 86.4 → 77.4 µs |
| `tetrahedral_stereo` | 21.5 → 12.9 µs | 18.5 → 16.5 µs | 34.6 → 30.7 µs |
| `cis_trans_stereo_bond` | 20.9 → 14.9 µs | 18.0 → 19.5 µs | 33.2 → 35.5 µs |
| `mixed_atom_and_bond_stereo` | 40.6 → 26.5 µs | 37.6 → 34.7 µs | 69.5 → 64.0 µs |
| `meso_dichlorobutane` | 57.3 → 42.7 µs | 52.8 → 39.3 µs | 102.7 → 75.5 µs |
| `para_stereo_trichloropentane` | 81.8 → 57.7 µs | 77.7 → 53.9 µs | 153.2 → 102.7 µs |
| `para_stereo_cascade` | 862.5 → 527.2 µs | 841.6 → 531.8 µs | 1.76 → 1.06 ms |
| `symmetric_stereo_triisopropanolamine` | 77.4 → 52.5 µs | 72.5 → 62.5 µs | 141.1 → 117.5 µs |
| `symmetric_stereo_butanol` | 87.6 → 73.6 µs | 84.9 → 86.1 µs | 166.0 → 163.0 µs |
| `symmetric_stereo_dimethylbutanol` | 195.8 → 176.5 µs | 193.2 → 194.1 µs | 390.9 → 381.8 µs |

`canonicalize` is 5% to 40% faster on all sixteen cases, and `para_stereo_cascade` returns to 1.3
times its time before S4c. On the five Topology cases `canonical_hash` is 10% to 14% slower and
`canonical_eq` 12% to 21% slower: these always build the key, and the remapped leaf costs more
there than the topology row leaf did. The study recorded this cost; S5e reads the Topology key
directly. `cis_trans_stereo_bond` is 7% to 8% slower in the two key operations, and the other ten
cases are unchanged or up to 40% faster in them.

#### S5e — Molecule Topology key **Done**

Add `topology_key`, which reads the Topology key from the entity order without remapping. The
Topology searches use it, with a debug assertion that it equals `leaf_key`. Results are unchanged
(green). [dep: S5d]

**Tests:** `test_topology_key`: equal to `leaf_key` on the Topology rows, including a bond outside
the largest bond color class.

`topology_key` builds the atom and bond rows of the schema from the entity order, and both Topology
searches use it under the debug assertion. `test_topology_key` has four rows, each under two search
orders. The library tests, the two integration tests, and the canonicalization properties pass with
the assertion active, and clippy is clean. The reference comparison reports no mismatch and the
molecule answers are identical.

Medians of the five Topology cases for the S0b baseline, S5d, and this subitem:

| Case | `canonicalize` | `canonical_hash` | `canonical_eq` |
| --- | ---: | ---: | ---: |
| `ordinary_naphthalene` | 21.3, 17.2, 17.3 µs | 21.0, 24.0, 21.4 µs | 36.0, 43.4, 35.0 µs |
| `disconnected_rings` | 34.4, 29.1, 29.1 µs | 33.3, 37.6, 34.4 µs | 60.0, 69.2, 60.0 µs |
| `feature_free_connected` | 37.8, 33.4, 33.2 µs | 38.6, 43.3, 40.1 µs | 72.6, 81.6, 74.1 µs |
| `feature_free_disconnected` | 40.4, 36.7, 36.2 µs | 41.4, 45.6, 42.9 µs | 75.6, 84.9, 78.7 µs |
| `symmetry_heavy_radicals` | 31.1, 29.6, 29.8 µs | 30.8, 34.6, 32.1 µs | 58.1, 65.0, 59.8 µs |

`canonical_hash` returns to 1.02 to 1.04 times the baseline and `canonical_eq` to 0.97 to 1.04
times; `canonicalize` keeps the gain of S5d.

#### S5f — Span Topology search **Done**

Add `topology_search_input` over the union graph and `topology_key` in
`canonicalize::reaction_span`. The Topology arm of `canonical_candidate` searches the topology
search input and reads the key from the entity order, with the debug assertion of S5e. Canonical
forms are unchanged relative to S5b (green). [dep: S0f, S5b]

**Tests:** `test_topology_search_input`: exact graph, colors, node sources, and partition, including
a removed and an added bond on one atom pair. `test_topology_key`: equal to the key of the remapped
span on the Topology rows.

`topology_search_input` and `topology_key` are in `canonicalize::reaction_span`, and the Topology
arm of `canonical_candidate` uses them under the debug assertion. The other arms rank their colors
and build their adapter themselves. `topology_initial_partition` had no other user and is removed.

`test_topology_search_input` has four rows: a direct edge, an added atom, a modified bond outside
the largest bond color class, and a removed and an added bond on one atom pair, of which the removed
one is the direct edge. `test_topology_key` has five rows, each under two search orders.
`test_canonical_candidate` gains a Topology row with the removed and the added bond.

The library tests, the two integration tests, and the canonicalization properties pass with the
assertion active, and clippy is clean. The reference comparison reports no mismatch, the 40,112
canonical span forms are identical to those of S5b, and the molecule answers are identical.

Medians of the five Topology cases for S5b and this subitem:

| Case | `canonicalize` | `canonical_eq` |
| --- | ---: | ---: |
| `ordinary_naphthalene` | 27.4, 21.3 µs | 62.5, 40.5 µs |
| `disconnected_rings` | 46.9, 34.0 µs | 100.9, 66.4 µs |
| `feature_free_connected` | 46.3, 37.5 µs | 114.1, 81.0 µs |
| `feature_free_disconnected` | 46.6, 40.2 µs | 115.3, 86.7 µs |
| `symmetry_heavy_radicals` | 32.2, 32.7 µs | 72.4, 64.4 µs |

Span `canonicalize` is 14% to 28% faster on four cases and unchanged on `symmetry_heavy_radicals`.
Span `canonical_eq` is 11% to 35% faster on all five. `overlay_heavy` and `tetrahedral_stereo`,
which do not take the Topology arm, are within 3% of S5b.

### S6 — Constraints and positional values in the colors

#### S6a — Field kinds **Done**

Add `BondFieldKind`, `DativeBondFieldKind`, `AromaticSystemFieldKind`, `MulticenterBondFieldKind`,
`NoncovalentBondFieldKind`, `StereoAtomFieldKind`, and `StereoBondFieldKind`, each beside its form,
with the variants of its `*FieldChange` enum. The eight field-kind enums derive `VariantArray`. The
key-and-kind entry of the nomenclature guide names them. The enums carry no behavior and no test.
Additive (green). [dep: none]

The seven enums stand beside their forms and are exported from `ir` as `AtomFieldKind` is. The eight
field-kind enums derive `VariantArray`. The key-and-kind entry of the nomenclature guide names them,
and the `*Kind` count of its suffix table rises by seven. Clippy is clean.

#### S6b — Reframe kind **Done**

**Module:** `ir/traits.rs`, the modules of the eight forms, and the modules of `ir/constraint`.
Additive (green). [dep: S6a]

Add `ReframeKind` with `Independent`, `Positional`, and `Relative` beside `FrameTransport`, and
`reframe_kind` on each of the eight field-kind enums and the eight constraint-key enums. The
electron counts of aromatic systems and multicenter bonds are `Positional`. The configuration of a
stereo entity and its ligand-symmetry, fluxionality, and topicity keys are `Relative`. Every other
field kind and constraint key is `Independent`. The nomenclature guide gains the term.

`reframe_kind` replaces `uses_participant_frame`. A constraint form is asked through its key. Each
of the six overlay deltas has `reframe_kind`: a `ModifyField` has the kind of its field, through
`kind` on the eight `*FieldChange` enums, a `ModifyConstraint` that of its key, and an `Add` or a
`Remove` is `Positional`, or `Relative` for a stereo entity.

**Tests:** one `test_<enum>_reframe_kind` table per enum with a row per variant. A property over
`reframe_by`: a field or a constraint whose kind is `Independent` is unchanged by every frame action
of its entity, and a `Positional` field is permuted with the participants.

`ReframeKind` stands after `FrameTransport` in `ir/traits.rs` and is exported from `ir`. Each of the
sixteen enums has `reframe_kind`, written as a match over every variant. No
`uses_participant_frame` remains: the constraint containers and the molecule-level constraints ask
the constraint key, and the callers in `ir/delta.rs` and `ir/reaction.rs` ask the delta.

Sixteen `test_<enum>_reframe_kind` tables hold 48 rows. The tables of the six constraint keys with a
participant frame replace the six `uses_participant_frame` tables row for row. Six
`test_<kind>_delta_reframe_kind` tables hold 30 rows, and eight `test_<kind>_field_change_kind`
tables hold 19.

Ten properties in `tests/property/frame.rs` check `reframe_kind` against `reframe_by`: the fields of
dative bonds, aromatic systems, multicenter bonds, and noncovalent bonds, and the constraints of the
six forms with a participant frame. The stereo forms have no `Independent` or `Positional` field.
With `Electrons` or the topicity key set to `Independent`, the corresponding properties fail.

The library tests, the property suite, and the integration tests pass, and clippy is clean. The
nomenclature guide has the entry *Reframe kind*.

#### S6c — Canonical equality pairs **Superseded by S6d**

`molecule_canonical_pair_strategy` was widened so that its derived molecule replaces one field or
one constraint of one entity of any of the eight kinds. The change was measured and is not
committed. With one value left out of the key, `test_molecule_canonical_eq_contract` failed on every
run for an atom field and for a bond field. For an aromatic-system field and for a noncovalent
constraint value it failed on one of three runs of 256 cases, and for the constraint value it
passed a run of 2,000 cases. The claim the widened pairs were to assert is the completeness property
of `AsColor`, which S6d to S6f assert form by form. The uncommitted change to the strategy is
discarded.

#### S6d — Color trait and the colors of four forms **Done**

**Module:** `ir/traits.rs`, `ir/atom.rs`, `ir/bond.rs`, `ir/dative.rs`, `ir/noncovalent.rs`.

Add `AsColor` beside `ReframeKind`, visible within the crate, with the associated type `Color<'a>`
bounded by `Ord` and `Hash` and the method `as_color`. Its documentation states invariance and
completeness. Add `AtomColor`, `BondColor`, `DativeBondColor`, and `NoncovalentBondColor`, each a
borrowed view of the whole form with derived equality, order, and hash, and implement `AsColor` on
the four forms. Additive (green). [dep: S6b]

**Tests:** `test_atom_form_as_color`, `test_bond_form_as_color`, `test_dative_bond_form_as_color`,
and `test_noncovalent_bond_form_as_color`. Each has a row of two equal forms with equal colors and,
for every field kind and for one constraint, a row of two forms that differ in it with unequal
colors. The trait is internal, so these tables assert its properties; the properties of canonical
forms cover them on generated aggregates.

`AsColor` follows `ReframeKind` in `ir/traits.rs`. Each color follows the last implementation of
its form: `AtomColor` and `BondColor` after `From<&str>`, `DativeBondColor` after `From<&str>`, and
`NoncovalentBondColor` after `FrameTransport`. The four tables have 8, 5, 3, and 3 rows, and the
12,157 library tests pass. The trait and the four colors have no caller outside the tests yet, so
the library build reports five dead-code warnings, which end with S6j for atoms and bonds and with
S6k for dative and noncovalent bonds. No attribute silences them.

#### S6e — Electron count at a position; colors of aromatic systems and multicenter bonds **Done**

**Module:** `ir/electrons.rs`, `ir/aromatic.rs`, `ir/multicenter.rs`.

Add `ElectronCountsForm::at`, visible within the crate. It takes an `AtomPosition` within the atom
list of the entity and returns the electron count of that atom, or `None` when the counts are
undetermined. Add `AromaticSystemColor` and `MulticenterBondColor`, which borrow the charge, the
unpaired electrons, and the constraints, and implement `AsColor` on the two forms. Additive (green).
[dep: S6d]

**Tests:** `test_electron_counts_form_at` with undetermined counts and with each position of a
literal list. `test_aromatic_system_form_as_color` and `test_multicenter_bond_form_as_color`: forms
that differ in the charge, the unpaired electrons, or a constraint have unequal colors; forms that
differ in the electron counts only have equal colors; a form and its image under a frame action
have equal colors.

`at` follows `permute` in the inherent implementation of `ElectronCountsForm`. It indexes the
literal list, so a position past the end of the list panics. The check of the list length against
the atom list stays with the caller, where `entity_color_key` has it today and where S6k keeps it.
Each color follows `FrameTransport` of its form and holds three references in the field order of
the form. `test_electron_counts_form_at` has 4 rows. Each color table has 6 rows: equal forms, the
charge, the unpaired electrons, a constraint, the electron counts, and the image under
`reframe_by`. The 12,173 library tests pass. The library build reports three more dead-code
warnings, for `at` and the two colors, which end with S6k.

#### S6f — Colors of stereo atoms and stereo bonds **Done**

**Module:** `ir/stereo.rs`.

Add `StereoAtomColor` and `StereoBondColor`, which hold the stereo kind and borrow the
stereogenicity constraint, and implement `AsColor` on the two forms. Additive (green). [dep: S6d]

**Tests:** `test_stereo_atom_form_as_color` and `test_stereo_bond_form_as_color`: forms of different
kinds, or with different stereogenicity constraints, have unequal colors; forms that differ in the
configuration, or in a ligand-symmetry, fluxionality, or topicity constraint, have equal colors; a
form and its image under a frame action have equal colors.

Each color holds the kind as `Option<StereoKind>`, read with `StereoConfigurationForm::kind`, and a
reference to the `StereogenicityForm` that `stereogenicity` of the constraints returns, which is
undetermined when the form has no such constraint. Each color follows the inherent implementation
of its form that holds `update` and `difference_to`, outside the `stereo_element!` macro. The tests
are in `ir/stereo/tests/atom_form.rs` and `ir/stereo/tests/bond_form.rs`. The stereo atom table has
9 rows: equal forms, two kinds, a kind against an undetermined configuration, the stereogenicity,
the coset, the three frame-relative constraints, and the image under `reframe_by`. The stereo bond
table has the same rows without the row of two kinds, since a stereo bond has one kind. The 12,190
library tests pass. The library build reports two more dead-code warnings, for the two colors,
which end with S6k.

#### S6g — Orders from the declaration **Done**

**Module:** `ir/delta.rs`, `ir/incidence.rs`.

Reorder the variants of `EntitySpan` to `Unchanged`, `Removed`, `Modified`, `Added`, with a comment
at the enum that states why `Added` is last. Implement `AsColor` for `EntitySpan<T>` for every `T`
that implements it: each side gives its color, and a `Modified` entity with equal colors gives
`Unchanged`. `Incidence` derives `PartialOrd` and `Ord` in place of its written order and of
`position`. Additive (green): no other code reads the order of `EntitySpan`, and the derived order
of `Incidence` equals the written one. [dep: S6d]

**Tests:** `test_entity_span_cmp` with the four tags in order. `test_entity_span_as_color` with a
row per variant and a `Modified` row with equal colors. `test_incidence_cmp` is rerun unchanged.

The reason for the position of `Added` is a "Semantic properties" paragraph in the documentation
of `EntitySpan`, since the derived order is public. The color of `EntitySpan<T>` is
`EntitySpan<T::Color<'a>>`, and the implementation follows `superimpose`. `Incidence` loses
`position` and its two written implementations, and its documentation states the declaration
order. `test_entity_span_cmp` has 3 rows, one per adjacent pair of tags, each with the larger
payload on the smaller tag. `test_entity_span_as_color` has 5 rows over spans of
`AromaticSystemForm`. It compares two colors with `cmp`, so each row fixes the place of a variant
in the order, and the `Modified` row with equal colors compares equal to `Unchanged`.
`test_incidence_cmp` passes unchanged. No expected value changed: the 12,198 library tests, the
integration tests of the crate, and the 555 properties pass. The library build reports the same
ten dead-code warnings as after S6f.

With this subitem `EntitySpan` and `ConstraintSpan` moved, with their implementations and 11
tests, from `ir/delta.rs` to the new module `ir/span.rs`. The text is unchanged, the public paths
are unchanged, and every import inside the crate names `span`.

#### S6h — Color ranking and the deferred key

**Module:** `ir/canonicalize.rs`.

Add `rank_colors`. It sorts the colors of one kind, assigns dense integers from a given first
integer, and returns the integer after the last one used. `DeferredKey` becomes generic over the
type of its key. Additive (green). [dep: none]

**Tests:** `test_rank_colors` with equal colors, distinct colors, an empty list, and a first integer
above zero. The three `test_deferred_key_*` tests are rerun unchanged.

#### S6h.1 — Benchmark cases for Constitution and Full

**Module:** `umol-graph-ir/benches/canonicalize.rs`. Additive (green). [dep: S4a]

The corpus gains four cases, so that Constitution and the generator check at Full are timed.

- `overlay_unconstrained`: the molecule of `overlay_heavy` without its stereo entities and without
  its constraints. It runs at Constitution.
- `aromatic_electron_counts`: a five-ring of carbon atoms whose aromatic system has the electron
  counts 2, 1, 1, 1, 1. It runs at Constitution, and the count distinguishes otherwise equal atoms.
- `symmetric_charge_sum`: a six-ring with a charge-sum constraint over the whole molecule. It runs
  at Full, and the constraint keeps the symmetry of the ring.
- `asymmetric_charge_sum`: the same ring with the constraint over two adjacent atoms. It runs at
  Full, and the constraint breaks the symmetry.

Their medians are recorded here before S6i, with the level at which every case of the corpus runs.
Today `ordinary_naphthalene`, `disconnected_rings`, `feature_free_connected`,
`feature_free_disconnected`, and `symmetry_heavy_radicals` run at Topology. `overlay_heavy` and
`frame_relative_stereo_constraint` carry a molecule-level constraint and run at Full. The other
nine cases run at Structure, and none runs at Constitution. The inline-constraint groups run at
Full until S8b and at the level of their entities after it.

#### S6i — Search partition from the colors

**Module:** `ir/canonicalize.rs`, `ir/canonicalize/molecule.rs`.

The incidence descriptors of `constitution_partition_descriptors` and of the molecule's
`structure_partition_descriptors` become the exact incidence colors, so dative roles and electron
counts split the search partition of both aggregates. The comment that gives the minimum over every
relabeling as the reason for coarser descriptors goes. Breaking for representatives of aggregates
in which these values distinguish otherwise equal atoms (red→green). [dep: none]

**Tests:** shown first: the modified aromatic system row of `test_canonical_candidate` visits 2
leaves without pruning and 1 with it, against 6 and 3, and
`test_canonicalize_constitution_entity_kind_minimum`, which compares the search with the minimum
over every relabeling, is removed. A row with a five-ring whose aromatic system has the electron
counts 2, 1, 1, 1, 1. The two Constitution cases are timed against S6h.1.

#### S6j — Molecule `canonicalize` at Topology

**Module:** `ir/canonicalize/molecule.rs`.

Add `colors`, which ranks the colors of the atoms and bonds of a normalized molecule with
`rank_colors`, and `color_key`, which reads the key of a leaf from its entity order: the integer of
each atom, then the endpoint images and the integer of each bond. The key is still the key of a
leaf, so its type keeps the name `CanonicalComparisonKey`: the molecule submodule defines its own,
and the encoded type of that name in `canonicalize.rs` is referred to by its path until S10a
removes it. `canonicalize_topology` normalizes the molecule once, searches the topology search
input built from the integers with deferred keys, and returns the remapped molecule. While the
encoded Topology path remains for `canonical_key_at`, a `debug_assert!` compares the two results
for molecules without inline constraints. Canonical forms are unchanged (green). [dep: S6d, S6h]

**Tests:** `test_colors` and `test_color_key` with exact integers for the fixtures of
`test_topology_key`. The Topology benchmarks of molecule `canonicalize` are timed against S5.

#### S6k — Molecule `canonicalize` at Constitution and Structure

**Module:** `ir/canonicalize/molecule.rs`.

`colors` ranks the colors of all eight kinds and of the incidences. It keeps the checks that need
the entity: the length of a literal electron list, the stereo kind, and the compatibility of a
stereo form with its ligand frame. `color_key` reads the rows of the six further kinds as decided
under Key rows, with each stereo form transported by `stereo_atom_representative_action` or
`stereo_bond_representative_action` and normalized. `canonicalize_constitution` and
`canonicalize_structure` normalize once, build the adapter and the search partition from the
integers, and search with deferred keys. `generator_preserves_stereo` transports the whole
stereo form. The `debug_assert!` of S6j covers both levels. Canonical forms are unchanged (green).
[dep: S6e, S6f, S6i, S6j]

**Tests:** rows of `test_colors` and `test_color_key` for each of the six kinds, with exact
integers and exact transported stereo forms. Rows of `test_retain_stereo_preserving_generators`: a
topicity constraint on one of two equivalent stereocentres, and the same constraint on both. The
stereo benchmarks of molecule `canonicalize` are timed against S5, and the two Constitution cases
against S6i.

### S7 — Full

#### S7a — Aggregate generator check

Add `remapping_from_generator` and `retain_preserving_generators` in each submodule. The first
builds the remapping that sends every entity to its image under a generator. The second keeps a
generator when the aggregate remapped by it equals the aggregate under framed equality. Additive
(green). [dep: none]

**Tests:** `test_remapping_from_generator` with exact remappings.
`test_retain_preserving_generators`: a six-ring with a charge-sum constraint over the whole
molecule, over one atom, and over two opposite atoms, with the ring rotations and reflections as
generators; a span with a constraint span on one side.

#### S7b — Molecule `canonicalize` at Full

**Module:** `ir/canonicalize/molecule.rs`.

`canonicalize_full` ranks the colors with `colors`, builds the search partition from the integers,
and runs the pruned search in backend canonical order with `retain_preserving_generators` as the
generator filter. The key of a Full leaf is the key of Structure followed by the
molecule-level constraints of the molecule remapped and reframed in the entity order of the leaf.
`canonical_key_at` keeps the encoded exhaustive search at Full until S8a. Breaking for
representatives of molecules that reach Full (red→green). [dep: S6k, S7a]

**Tests:** Full rows of `test_color_key`. Rows of `test_molecule_canonical_eq`: the six-ring
placements of S7a against renumbered copies and against each other. `overlay_heavy`,
`frame_relative_stereo_constraint`, the two charge-sum cases, and the inline-constraint benchmarks
of molecule `canonicalize` are timed against S6.

### S8 — Molecule key and level

#### S8a — Molecule key on colors

**Module:** `ir/canonicalize/molecule.rs`.

The molecule's `CanonicalComparisonKey` gains the distinct colors of each kind and of the
incidences, in rank order, beside its rows. `Molecule::canonical_key_at` takes a normalized
molecule, runs the search of its level as `canonicalize` does, and returns the key of the minimum
leaf, which borrows from the molecule.
`canonical_eq` and `canonical_hash` normalize a copy of each molecule, keep two
contradictions equal, and compare or hash the keys. The `debug_assert!` of S6j goes. Delete the
molecule's `comparison_key`, `topology_key`, `leaf_key`, `initial_color_keys`, `entity_color_key`,
`constraint_blocks`, `stereo_refinement_descriptor`, `constitution_entity_classes`,
`structure_partition_descriptors`, `para_stereo_partition_descriptors`, and `structure_partition`,
with their tests. Breaking for the keys and the hash values of molecules (red→green). [dep: S7b]

**Tests:** rows of `test_molecule_canonical_key_at` at each level with exact keys and exact
visited-leaf counts; the Full rows equal the exhaustive minimum. The exact keys of the deleted key
tests are shown first, beside the new keys of the same fixtures. `test_molecule_canonical_eq` and
its siblings are rerun unchanged. The benchmarks of molecule `canonical_eq` and `canonical_hash`
are timed against S7.

#### S8b — Molecule level from the entities

**Module:** `ir/canonicalize.rs`.

`molecule_canonicalize_level` reads no inline constraint: a molecule-level constraint gives Full, a
stereo entity Structure, a dative bond, aromatic system, multicenter bond, or noncovalent bond
Constitution, and atoms and bonds alone Topology. Breaking for representatives and keys of
molecules with inline constraints (red→green). [dep: S8a]

**Tests:** shown first: the molecule rows that expect a constraint to select Full, in
`test_molecule_canonicalize_description_level`, `test_canonicalize_canonical_eq_description_level`,
and the two excluded-data tests. Rows of `test_molecule_canonical_eq`: a six-ring with a valence
constraint on one atom, or an aromatic constraint on one bond, against a renumbered copy and
against another placement; biphenyl with `#e6` on one aromatic system against a renumbered copy and
against the other system; 2,3-dichlorobutane with a stereogenicity or a topicity constraint on one
centre against a renumbered copy, against the same constraint on the other centre, and against
another ligand pair. The two Python test files are rerun. The inline-constraint benchmarks are
timed against S8a.

### S9 — Reaction spans and reactions

#### S9a — Span colors and key rows

**Module:** `ir/canonicalize/reaction_span.rs`.

Add `colors`, `color_key`, and the span submodule's own `CanonicalComparisonKey`, the counterparts
of the molecule's.
The color of a span entity is the color of its `EntitySpan` of forms, and the positional part of an
incidence is lifted through the span in the same way. A stereo row holds the `EntitySpan` of the
transported forms. Additive (green). [dep: S6g, S6k]

**Tests:** `test_colors` and `test_color_key` with exact integers for the fixtures of
`test_canonical_candidate`, among them a `Modified` aromatic system whose sides differ in the
electron counts only, and an added atom.

#### S9b — Span search on colors

**Module:** `ir/canonicalize/reaction_span.rs`.

`canonical_candidate` ranks the colors with `colors`, builds the topology search input or the
adapter and the search partition from the integers, and searches with deferred keys at every
level. At Full it runs the pruned search in backend canonical order with
`retain_preserving_generators`, and the key adds the constraint spans of the span remapped and
reframed in the entity order of the leaf. `generator_preserves_stereo` transports the whole stereo
form on each side. `ReactionSpan::canonical_key_at` takes a normalized span. `canonical_eq` and
`canonical_hash` of spans and of reactions normalize a copy and compare or hash the keys.
Delete the span's `comparison_key`, `topology_key`, `initial_color_keys`, and `entity_color_key`,
with their tests. Breaking for the keys and the hash values of spans and reactions, and for
representatives of spans that reach Full (red→green). [dep: S7a, S9a]

**Tests:** rows of `test_canonical_candidate` at each level with exact keys and exact
visited-leaf counts without and with pruning; at Full a constraint span that keeps the symmetry and
one that breaks it. The exact keys of the deleted key tests are shown first. Rows of
`test_retain_stereo_preserving_generators` with a constraint on one side of a span only. The span
and reaction benchmarks are timed against S8.

#### S9c — Span and reaction level from the entities

**Module:** `ir/canonicalize.rs`.

`reaction_span_canonicalize_level`, `delta_canonicalize_level`, `entity_span_canonicalize_level`,
and `reaction_canonicalize_level` read no inline constraint, and `level_with_constraints` is
removed. Breaking for representatives and keys of spans and reactions with inline constraints
(red→green). [dep: S8b, S9b]

**Tests:** shown first: the rows of `test_delta_canonicalize_level`,
`test_entity_span_canonicalize_level`, `test_reaction_span_canonicalize_description_level`, and
`test_reaction_canonicalize_description_level` that expect a constraint to select Full. Rows of
`test_reaction_span_canonical_eq_renumbering` with the constrained fixtures of S8b lifted to spans,
and with a constraint on one side of a span only.

### S10 — Encoding

#### S10a — Remove the encoding

**Module:** `ir/canonicalize.rs`.

Delete `CanonicalKeyValue`, the encoded `CanonicalComparisonKey`, the positioned keys and their
positions, the span key and its tag positions, the encoding functions, `atom_inherent_fields`,
`bond_inherent_fields`, `InitialColorKey`, `rank_initial_colors`,
`constitution_partition_descriptors`, and the span key helpers, with their tests. The key type of
each aggregate keeps the name `CanonicalComparisonKey`. No behavior changes (green). [dep: S9c]

**Tests:** the library tests, the property suite, and the integration tests pass.

#### S10b — Guides

**Module:** `docs/development/data-types.md`, `docs/development/nomenclature.md`,
`ir/incidence.rs`.

The data-type guide loses its canonical comparison schema and its rule on explicit positions. Its
passage on constraints and structural orbits is replaced by the decisions on colors, inline
constraints, and the search of each level. The nomenclature guide's entry *Canonical and
canonicalize* describes a color as the frame-independent part of a form and the integer it is
ranked into. Its entry *Normalize* no longer names `Normalized<T>`, which does not exist in the
code. The documentation of `Incidence` and of the two `incidence_graph` functions no longer refers
to the schema. [dep: S10a]

### S11 — Para-stereo setting

#### S11a — Remove the setting

Remove `para_stereo` from `CanonicalizeContext`. `umol-graph` stops passing
`StereoModel::para_stereo` to it, and the `umol-py` binding, the tests, the properties, and the
benchmarks drop the field. The data-type guide and the Python API guide follow. Breaking for the
public API (red→green). [dep: S8a]

**Tests:** `test_canonicalize_config_context` of `umol-graph` and `test_canonicalize_context` of
`umol-py` lose their `para_stereo` rows. The Python suite is rerun.

### S12 — Closeout

#### S12a — Final comparison and measurements

The reference comparison runs on the full pair file. The benchmark groups of S0a run in the
reference checkout and in the tree back to back. The table is recorded here. [dep: S11a]

#### S12b — Gates

`cargo +nightly fmt --all`; `cargo clippy --workspace --all-targets -- -D warnings`;
`cargo test --workspace` with all features; the graph-IR property suite; the Python suite under the
python-build skill. [dep: S12a]

#### S12c — Remove the reference setup

Remove `/Users/dr/.cargo-target-244`, `/Users/dr/.cargo-target-244-1`, `scratch/244-compare`,
`scratch/244-compare-1`, and the pair file. The three worktrees are removed by the user.
[dep: S12b]

### Stage gate

After every stage: the canonicalization and reaction-span unit tests, the `canonicalization` and
`reaction_span` integration tests, the canonicalization properties, graph-IR clippy, and the
reference comparison. Stages that change representatives also run the two Python test files.

### Critical path

S4 → S5 → S6 → S7 → S8 → S9 → S10 → S11 → S12. No stage is deferrable: S6 to S9 carry the colors
and the keys of both aggregates, S10 removes the encoding that nothing uses after S9, and S11
removes a setting that nothing reads after S8a.
