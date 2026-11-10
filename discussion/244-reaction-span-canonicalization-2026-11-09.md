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

## Staged implementation plan

Work proceeds one subitem at a time. After each subitem the narrow tests run and work stops for
review. After each subitem that changes search behavior, the reference comparison runs first.
Molecule functions are moved in S0 and edited only in S6. Expected values of exact cases are shown,
old beside new, before any is edited. Each subitem reports its changed files and a suggested commit
message; no git command is run.

After S0, molecule code is in `umol-graph-ir/src/ir/canonicalize/molecule.rs` and span and reaction
code in `umol-graph-ir/src/ir/canonicalize/reaction_span.rs`, each with its own `tests.rs`. A span
function takes the name of its molecule counterpart. From S1 on, a bare function or test name means
the one in `reaction_span` unless the subitem says otherwise. Tests follow the test-writing skill
and sit in source order.

The plan adds roughly 900 lines of span code and removes about 290 lines of span code. The estimate
comes from the sizes of the molecule counterparts. S0 also relocates about 1,450 lines of molecule
code and about 650 lines of span code without changing them, and splits the 6,594-line tests file.

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

#### S3c — Span stereo generator check

Add `generator_preserves_stereo` and `retain_stereo_preserving_generators`, the counterparts of the
two molecule functions. One frame action must carry the configuration of each present side to that
of the image. Additive (green). [dep: S0e]

**Tests:** `test_retain_stereo_preserving_generators`: the rows of the molecule table as unchanged
spans; a generator that preserves the lhs configuration and not the rhs one; a generator between a
modified and an unchanged stereo entity; a stereo-bond endpoint exchange.

#### S3d — Structure search

The Structure arm of `canonical_candidate` starts from `structure_partition` and runs
`canonical_search_with_generator_filter` with the span generator check. The span methods pass
pruning with backend canonical order at Structure. Breaking for representatives of spans with stereo
entities (red→green). [dep: S2b, S3b, S3c]

**Tests:** Structure rows of `test_canonical_candidate`: three stereocentres with equal
configurations; one with a different lhs configuration; one with a different rhs configuration
only. Stereo rows of `test_reaction_span_tracked_canonicalize_at` under both `para_stereo`
settings.

#### S3e — Properties under both para-stereo settings

**Module:** `umol-graph-ir/tests/property/reaction/span/canonicalize.rs` and
`.../reaction/canonicalize.rs`. Each property draws `para_stereo` and builds its context from it.
The laws are unchanged. Additive (green). [dep: S3d]

### S4 — Full level

#### S4a — Full search from the structure partition

The Full arm of `canonical_candidate` starts from `structure_partition`. It stays exhaustive in
node-id order. Breaking for representatives of spans with constraints and stereo entities
(red→green). [dep: S3b]

**Tests:** Full rows of `test_canonical_candidate` and `test_reaction_span_tracked_canonicalize_at`:
inline constraints that differ between sides, a constraint span, and stereo entities with
constraints.

#### S4b — Normalization at Full only

`ReactionSpan::tracked_canonicalize_at` normalizes the span at Full only. No canonical output
changes (green). [dep: S4a]

**Tests:** rows of `test_reaction_span_tracked_canonicalize_at` with non-normal values below Full,
including a modified entity with equal sides, and a contradictory value at each level.

### S5 — Leaf rows and topology search input

Each level replaces the remap and reframe of the whole span at every leaf by rows built from the
atom order. The rows keep the shapes `comparison_key` builds today. Non-atom entities take their
order from their rows. Each wiring subitem is breaking for representatives (red→green).

#### S5a — Topology functions

Add `topology_search_input` over the union graph, and `topology_bond_key_rows` and
`topology_candidate`. Additive (green). [dep: S0f, S1b]

**Tests:** `test_topology_search_input`: exact graph, colors, node sources, and partition, including
a removed and an added bond on one atom pair. `test_topology_candidate`: exact key and entity order.

#### S5b — Topology search

The Topology arm uses the topology search input and the topology candidate. [dep: S2b, S5a]

#### S5c — Constitution functions

Add `constitution_candidate`. Additive (green). [dep: S5a]

**Tests:** `test_constitution_candidate`: exact key and entity order for each
non-stereo overlay kind, unchanged, added, removed, and modified, with participant-indexed electron
counts on both sides.

#### S5d — Constitution search

The Constitution arm uses the constitution candidate. [dep: S5b, S5c]

#### S5e — Structure functions

Add `structure_candidate`. The stereo rows use `canonical_kinded_stereo_frame` on each side, which
selects one ligand frame for both. Additive (green). [dep: S5c]

**Tests:** `test_structure_candidate`: exact key and entity order for stereo atoms and stereo bonds,
unchanged, added, removed, modified, and with one undetermined side.

#### S5f — Structure search

The Structure arm uses the structure candidate. [dep: S3d, S5d, S5e]

#### S5g — Full functions and search

Add `constraint_blocks`, which is the constraint part of `comparison_key` unchanged, and
`complete_candidate`, the counterpart of the molecule function of that name. The Full arm uses the
complete candidate. Delete `comparison_key`. [dep: S4b, S5f]

**Tests:** `test_constraint_blocks` and `test_complete_candidate`: exact keys for inline constraints
on each entity kind and for constraint spans.

### S6 — Trials

Each trial is kept or reverted on its measurement against the S0b baseline.

#### S6a — Molecule key search at Structure

The Structure branch of `Molecule::canonical_key_at` runs `canonical_search_with_generator_filter`
with `retain_stereo_preserving_generators`, with pruning and backend canonical order. Results are
unchanged (green). [dep: S0d]

**Tests:** `test_molecule_canonical_key_at` in the molecule tests: for the molecule stereo fixtures,
the key equals the exhaustive minimum.

#### S6b — Full-level pruning for reaction spans

Add `generator_preserves_constraints` and `retain_constraint_preserving_generators` in
`canonicalize::reaction_span`. A span without constraint spans searches at Full with the stereo and
constraint checks and with pruning. [dep: S3c, S5g]

**Tests:** `test_retain_constraint_preserving_generators`: an inline constraint on each entity kind,
on every equivalent entity and on one of them; a constraint on one side only; frame-relative
constraints on stereo entities. Full rows of `test_canonical_candidate` with exact leaf counts, and
one row with a constraint span that stays exhaustive.

#### S6c — Full-level pruning for molecules

Add `generator_preserves_constraints` and `retain_constraint_preserving_generators` in
`canonicalize::molecule`. A molecule without molecule-level constraints searches at Full with the
stereo and constraint checks and with pruning, in `canonicalize_full` and in the Full branch of
`Molecule::canonical_key_at`. Canonical forms are unchanged (green). [dep: S6a]

**Tests:** `test_retain_constraint_preserving_generators` in the molecule tests, with the molecule
rows of S6b. Full rows
comparing the pruned search with the exhaustive minimum.

### S7 — Closeout

#### S7a — Final comparison and measurements

The reference comparison runs on the full pair file. The benchmark groups of S0a run in both
checkouts back to back. The table is recorded here. [dep: S5g, S6]

#### S7b — Gates

`cargo +nightly fmt --all`; `cargo clippy --workspace --all-targets -- -D warnings`;
`cargo test --workspace` with all features; the graph-IR property suite; the Python suite under the
python-build skill. [dep: S7a]

#### S7c — Remove the reference setup

Remove `/Users/dr/.cargo-target-244`, `scratch/244-compare`, and the pair file. The reference
worktree is removed by the user. [dep: S7b]

### Stage gate

After every stage: the canonicalization and reaction-span unit tests, the `canonicalization` and
`reaction_span` integration tests, the canonicalization properties, graph-IR clippy, and the
reference comparison. Stages that change span representatives also run the two Python test files.

### Critical path

S0 → S1 → S2 → S3 → S4 → S5 → S7. S2 to S4 give the pruned, stereo-sensitive search with
about 400 lines of new code. S5 removes the per-leaf remap and reframe and is measured on its own. S6 is
deferrable: the span work is complete without it, and S6a can run at any point after S0.
