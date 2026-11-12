# 245 — Cost of canonical equality

Status: In Progress
Date: 2026-11-12
Relates: [244](244-reaction-span-canonicalization-2026-11-09.md),
[data types guide](../docs/development/data-types.md)

## Purpose

`canonical_eq` of `Molecule` and of `ReactionSpan` canonicalizes both values and compares the
canonical forms. Every pair that is not structurally equal costs about two canonicalizations,
whether it is equal or not.

This document scopes a `canonical_eq` that decides most pairs with less work:

- unequal pairs are rejected by checks ordered by cost;
- the remaining pairs are decided by a remapping between the two molecules, which is checked: the
  one that the refined cells give, or one found by a search of the second molecule;
- canonical forms are compared only for a pair at two different description levels.

The aim is to narrow the set of pairs that reach the comparison of canonical forms. When N
molecules are compared with each other, on the order of N of the N² pairs are equal or nearly so,
and the others should end at a cheap check. It is not an aim that every comparison costs one
canonicalization or less.

The document covers `Molecule` and `ReactionSpan`, and `Reaction` through its reaction span. The
text speaks of molecules. Each statement holds for reaction spans with the span counterpart of a
named function, and differences are stated where they exist. For values whose canonicalization
succeeds, `canonical_eq` relates the same pairs as before. `canonical_hash` is not changed.

## Current state

- `canonical_eq` returns true for structurally equal molecules. Otherwise it clones and
  canonicalizes both and compares the results with `==`. Two contradictions compare equal. A
  contradiction and a canonical form do not.
- `normalized_eq` and `framed_eq` have the same shape: the structural shortcut, then both sides
  cloned and normalized or reframed, with two contradictions equal.
- `Molecule::framed_eq_under` clones the receiver, remaps it, and calls `framed_eq`. A reaction
  span has no such method; its generator check at Full clones the span, calls `try_remap`, and
  calls `framed_eq`.
- Derived `==` compares the fields in order, and the length of a container before its contents.
- The counts of the eight entity kinds are read by the private `molecule_counts` and
  `reaction_span_counts` in `ir/canonicalize.rs`. Normalization, reframing, and remapping do not
  change them.
- The color of a span entity is the color of its entity span, and the color of a span incidence
  is its role and the span of the electron counts.
- `Reaction` canonicalizes through its reaction span.
- `initial_colors` assigns dense integers to the colors of each entity kind and of the
  incidences, in the order of the colors. `OrderedPartition::refine` splits a cell by the
  neighbor counts of its nodes and orders the parts by those counts. The integers and the order
  of the cells do not depend on the entity ids.

### Measurements

Shares of one molecule canonicalization, from a sampling profile of the benchmark binary. The
remainder of a row is not attributed. Reaction spans were not profiled; the trial measures
them directly.

| Case | Total | Search | of which backend | Remap | Normalize | Reframe | Colors, adapter, incidence |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `tetrahedral_stereo` | 7.9 µs | 19% | 0% | 29% | 15% | 13% | 19% |
| `overlay_heavy` | 18.2 µs | 24% | 0% | 30% | 12% | 13% | 17% |
| `aromatic_electron_counts` | 11.7 µs | 45% | 8% | 29% | 9% | 4% | 9% |
| `ordinary_naphthalene` | 12.2 µs | 51% | 35% | 19% | 15% | 0% | 12% |
| `feature_free_connected` | 25.3 µs | 77% | 63% | 9% | 7% | 0% | 5% |
| `symmetric_stereo_dimethylbutanol` | 149.1 µs | 90% | 60% | 2% | 2% | 2% | 3% |
| `asymmetric_charge_sum` | 135.3 µs | 96% | 10% | 1% | 0% | 0% | 1% |

Over the twenty cases of the benchmark corpus, remapping takes 1.4% to 32% and incidence
construction 0.2% to 3.3%.

In the first two cases the search takes a fifth to a quarter, and remapping, normalization, and
reframing take more than half. In the last three the search takes most of the time. In
`asymmetric_charge_sum`, which runs at Full, the generator check takes 37% and the leaf keys 29%.

## Design

### Order of the checks

1. Different entity counts give false.
2. Structural equality gives true.
3. Both molecules are normalized and `initial_colors` is built for each. A different multiset of
   colors in an entity kind or in the incidences gives false. When the multisets agree, an
   integer stands for the same color in both molecules.
4. The initial partitions of both molecules are refined. Cells that differ in number or in size,
   position by position, give false.
5. When every cell holds one node, the cells give a remapping between the nodes of the two
   incidence graphs, and the check of that remapping decides the pair.
6. Otherwise one remapping is tried along the path of first nodes. A remapping that passes the check
   gives true.
7. Otherwise the second molecule is searched. Each leaf of the search gives a remapping, and a
   remapping that passes the check gives true. When no leaf passes, the answer is false.

A remapping relates two molecules when `framed_eq_under` holds under it. The quantities of
steps 1, 3, and 4 are the same for two molecules that a remapping relates, and a remapping that
passes the check relates them. In step 5 the remapping is the only one that can relate the two
molecules, so its failure shows that none exists. In step 7 the search visits every leaf of the
second molecule up to the automorphisms of that molecule, so a remapping that relates the two
molecules appears at a visited leaf.

Steps 3 to 7 run when both molecules select the same description level. With equal entity counts
the levels differ only when one molecule has molecule-level constraints and the other has none,
or one span has constraint spans and the other has none. Such a pair is compared through
canonical forms. A value that fails normalization or the colors is treated by the rule for
contradictions below.

The search and `initial_colors` are used as they are. The search receives a leaf function that
checks the remapping of a leaf, in place of the one that builds the leaf key. Its value is the same
for two leaves that an automorphism of the second molecule relates, which is what the search
requires of a key.

### Remapping

When every cell of the refined partitions holds one node, the cells give the remapping, and the
check decides the pair either way. This is the case for molecules without symmetry.

When cells with several nodes remain, the first molecule is followed along one path: the first
such cell is taken, its first node is individualized, the partition is refined, and this repeats
until every cell holds one node. The second molecule is followed along the same path for as long
as its cells agree with those of the first. When they agree to the end, the two orders give a
remapping, and a remapping that passes the check gives true.

When the cells diverge or that remapping fails, the second molecule is searched with the search and
the generator filter of its description level. Each leaf is paired by position with the end of
the path of the first molecule, and the remapping is checked.

For a symmetric molecule with stereo entities or molecule-level constraints, the first remapping
can fail for an equal pair: the cells are those of the colored incidence graph, which has more
symmetry than the molecule. The search then reaches a remapping that passes.

### Check of a remapping

| Level | Check |
| --- | --- |
| Topology, Constitution | Each node maps to a node with the same integer, and the incidences between two nodes map to incidences with the same integers between the image nodes. |
| Structure | The above, and for each stereo entity: the site maps to the site of the image, and the form transported to the ligand order of the image equals the form of the image. |
| Full | `framed_eq_under` with the remapping. |

- **Topology and Constitution.** The check reads the two incidence graphs and the integers of
  step 3 and makes no copy. The color of an entity is its form, the color of an incidence is its
  role and electron count, and the search at these levels runs without a generator filter. The
  colored incidence graph therefore holds what the comparison needs.
- **Structure.** `generator_preserves_stereo` makes this comparison for one molecule and one
  generator, and its span counterpart for one span. Each takes a second value as the target and
  is renamed `remapping_preserves_stereo`. The generator filter passes the same value twice.
- **Full.** No comparison of molecule-level constraints or constraint spans under a remapping
  exists without a copy, and none is written here.

The check decides in both directions, so an error in it gives a wrong answer. `framed_eq_under`
defines the answer: the check and `framed_eq_under` agree on every remapping.

### Entity counts

`pub fn entity_counts(&self) -> [usize; EntityKind::COUNT]` on `Molecule` and on `ReactionSpan`
returns the eight entity counts, indexed by entity kind. It replaces the private
`molecule_counts` and `reaction_span_counts`. It has no per-kind form and no Python binding.

`==` stays derived. In the other three equalities the count comparison runs first and saves:

| Equality | Saved when the counts differ |
| --- | --- |
| `normalized_eq` | Two clones and normalizations. |
| `framed_eq` | Two clones and reframings. |
| `canonical_eq` | Two canonicalizations. |

### Functions

Two names are public. `entity_counts` is described above. `ReactionSpan::framed_eq_under` is
added with the contract of `Molecule::framed_eq_under`; the generator check of spans at Full
calls it in place of its inline clone, `try_remap`, and `framed_eq`. The other functions are
private, and no type is added.

| Name | Where | Role |
| --- | --- | --- |
| `canonical_eq_at(&self, other, level, context) -> bool` | method on `Molecule` and on `ReactionSpan` | Steps 3 to 7 at one description level, beside `tracked_canonicalize_at`. `canonical_eq` calls it after the count comparison. |
| `same_colors` | both submodules | Step 3: whether the color integers of two normalized values stand for the same colors, each equally often. |
| `same_entity_color` | both submodules | Whether one entity of each value has the same color. |
| `incidence_color` | both submodules | The color of one incidence: its role and electron count. Moved out of `initial_colors`, which calls it. |
| `color_representatives` | `ir/canonicalize.rs` | Checks that two lists of color integers hold each integer equally often, and gives one position per integer in each. |
| `remapping_exists` | `ir/canonicalize.rs` | Steps 4 to 7 for two adapters. The caller supplies the check and the generator filter of the level. |
| `same_cells` | `ir/canonicalize.rs` | Step 4: whether two partitions agree in the number, size, and color of their cells. |
| `remapping_preserves_incidences` | `ir/canonicalize.rs` | The check at Topology and Constitution, beside `remapping_preserves_stereo`. |

Existing canonicalization code changes in two places: the stereo function takes a target and is
renamed, and the incidence color moves out of `initial_colors` into `incidence_color`. The search
is called through `canonical_search_impl` and is not changed.

### Reaction

`Reaction::canonical_eq` converts both reactions to reaction spans and compares the spans. When
a conversion fails, the rule for contradictions below applies. In the trial this agrees with the
equality of canonical reactions on the unit and property tests of reaction canonicalization.

### Contradictions

At present two contradictory values compare equal under `normalized_eq`, `framed_eq`, and
`canonical_eq`, whatever their entity counts. A false then needs one of the two values known to
be free of contradiction, which a count comparison ahead of normalization does not have.

The rule changes for `Molecule`, `ReactionSpan`, and `Reaction`. Structurally equal values are
equal under each equality. Otherwise:

| Equality | Both values pass the stage | One fails | Both fail |
| --- | --- | --- | --- |
| `normalized_eq` | Normal forms are equal. | false | false |
| `framed_eq` | Reframed forms are equal. | false | `normalized_eq` |
| `canonical_eq` | Canonical forms are equal. | false | `framed_eq` |

Each equality is reflexive, symmetric, and transitive, and each implies the next. A contradictory
value is equal to the values that share its normal form or its reframed form, where it has one,
and otherwise only to a structurally equal value. A contradictory molecule and a renumbered copy
of it are unequal.

The step back to the finer equality is needed because a value can pass one stage and fail a later
one. An aromatic system or multicenter bond whose literal electron counts differ in length from
its atom list, and a stereo entity whose form does not fit its ligand count, are found while the
colors are built. Without the step back, two such values with equal normal forms would be equal
under `normalized_eq` and unequal under `canonical_eq`.

Forms keep the present rule: two contradictory forms are equal under `normalized_eq`. The three
aggregates override the defaults of `Normalize::normalized_eq` and `Reframe::framed_eq`.

With this rule the count comparison is the first step for every input in the three equalities.

What states the present rule for the aggregates and changes with it:

- the data types guide and the nomenclature guide, in the passages on the totalization of
  equality for two intrinsic contradictions;
- `test_molecule_canonical_eq_contradiction` and `test_reaction_canonical_eq_contradiction` in
  `ir/canonicalize/tests.rs`, whose pairs of contradictory values become unequal;
- the property tests of `canonical_eq` on contradictions in `tests/property/molecule`,
  `tests/property/reaction`, and `tests/property/reaction/span`.

### Stored summaries

A summary of the composition that is stored per molecule is not part of this design. Inside one
call, step 3 rejects the same pairs for the same work. For comparisons over a stored collection,
`canonical_hash` decides more pairs and exists. In the trial a reject at step 3 costs 0.03 to
0.70 of one canonicalization, which is the cost a stored summary would have to improve on.

## Benchmarks

The equality benchmarks of the molecule and reaction-span groups hold equal pairs, apart from
the molecule scaling cases. Unequal pairs are added for both, one kind for each step that can
reject or fail:

- different entity counts;
- equal entity counts and different colors;
- equal colors and different constitution;
- equal constitution and different stereo configuration.

## Trial

The flow above was written in a scratch worktree for molecules and reaction spans and measured on
the twenty cases of the benchmark corpus. A span case is the corpus molecule with a charge change
and, where it has stereo entities, a configuration change. Each case gives an equal pair, by
reversed entity ids, and one unequal pair for each kind listed under Benchmarks. Three pairs of
carbon skeletons with different constitution were added for molecules.

Cost in canonicalizations of one value. The present `canonical_eq` costs about 2.0 for every
pair.

| Pair | Decided at | Molecules | Spans |
| --- | --- | ---: | ---: |
| Different entity counts | step 1 | about 10 ns | about 10 ns |
| Different colors | step 3 | 0.03 to 0.60 | 0.09 to 0.70 |
| Different constitution, cells differ | step 4 | 0.66 to 0.84 | not built |
| No symmetry, below Full, equal or different in stereo | step 5 | 1.13 to 1.27 | 1.07 to 1.29 |
| No symmetry, Full, equal or different in stereo | step 5 | 1.44 to 1.98 | 1.39 to 1.98 |
| Symmetric, equal, first remapping passes | step 6 | 0.29 to 1.03 | 0.35 to 1.05 |
| Symmetric, equal, first remapping fails | step 7 | 1.26 to 1.37 | 1.28 to 1.59 |
| Symmetric, unequal past the cells | step 7 | 1.17 to 2.07 | 1.28 to 1.38 |

Number of pairs decided at each step:

| Pair | Step 5 | Step 6 | Step 7 |
| --- | ---: | ---: | ---: |
| Molecules, equal | 6 | 10 | 4 |
| Molecules, different in stereo | 4 | | 6 |
| Molecules, different constitution past the cells | | | 1 |
| Spans, equal | 9 | 8 | 3 |
| Spans, different in stereo | 6 | | 4 |

Findings:

- No pair costs more than the present `canonical_eq`. The cost is 0.15 to 0.98 of the present
  cost for molecules and 0.17 to 0.97 for spans.
- No pair of the trial reaches the comparison of canonical forms.
- The result agrees with the present `canonical_eq` on all 143 pairs. The unit and property tests
  of molecule and span canonicalization pass against the trial code, apart from
  `test_molecule_canonical_eq_contradiction`, whose expected value changes with the rule for
  contradictions.
- A trial that stopped after step 6 and compared canonical forms for the remaining pairs cost
  2.2 to 3.5 for those pairs, above the present cost. Step 7 replaces that comparison.

Costs that the design leaves in place:

- At Full the check copies, remaps, and reframes, which gives the 1.4 to 2.0 of the table.
- A symmetric unequal pair pays the first remapping and a complete search. The two highest costs
  are 2.04 and 2.07, against 2.56 and 2.79 at present.
- The search visits its remaining leaves after one remapping has passed.

## Staged implementation plan

The flow exists as an uncommitted trial in the worktree `scratch/245-canonical-eq`, under trial
names. Each subitem moves the part it names into the tree under the settled names and adds its
tests. The trial changes three source files by about 830 added and 60 removed lines and the
benchmark file by about 270 added lines.

Work proceeds one subitem at a time. After each subitem the tests of `umol-graph-ir` run and work
stops for review. Expected values of existing tests are shown, old beside new, before any is
edited. Each subitem reports its changed files and a suggested commit message; no git command is
run. Tests follow the test-writing skill and sit in source order. A function that exists for both
aggregates is changed in both submodules in the same subitem.

### S0 — Entity counts, `framed_eq_under` for spans, and benchmark pairs

No result of an equality changes in this stage.

#### S0a — `entity_counts` **Done**

**Module:** `ir/molecule.rs`, `ir/reaction_span.rs`, `ir/canonicalize.rs`, both submodules,
`benches/canonicalize.rs`.

Add `pub fn entity_counts(&self) -> [usize; EntityKind::COUNT]` to `Molecule` and to
`ReactionSpan`. Remove `molecule_counts` and `reaction_span_counts` from `ir/canonicalize.rs`;
`remapping_from_order` and `remapping_from_generator` of both submodules call the method. The
`molecule_counts` of the benchmark file goes the same way. Additive, then a rewire that changes
no behavior (green).

**Tests:** `test_molecule_entity_counts` and `test_reaction_span_entity_counts`, with a value
that has entities of every kind and one that has atoms only.

The tests have one row per entity kind, each a compound with entities of that kind only beside
atoms and bonds, so that every position of the array is read: helium, ethane, ammonia borane,
benzene, diborane, the water dimer, bromochlorofluoroiodomethane, and but-2-ene. The span test
adds a row with an added atom and one with a removed atom, which count in the union frame. Five
test helpers that called the removed functions call the method. The tests of `umol-graph-ir`
pass with all features.

#### S0b — `ReactionSpan::framed_eq_under` **Done**

**Module:** `ir/reaction_span.rs`, `ir/canonicalize/reaction_span.rs`.

Add `pub fn framed_eq_under(&self, other: &Self, remapping: &MoleculeRemapping) -> bool` with the
contract of `Molecule::framed_eq_under`. The span's `retain_preserving_generators` calls it in
place of its inline clone, `try_remap`, and `framed_eq`. Additive, then a rewire that changes no
behavior (green).

**Tests:** rows that mirror the `test_molecule_framed_eq_under_*` unit tests on spans, and the
identity, inverse, and composition properties that `tests/property/molecule/comparison.rs` holds
for molecules.

Three unit tests. `test_reaction_span_framed_eq_under` compares the span of the hydrolysis of
fluoromethane with itself under identity, with a renumbered copy under the remapping of the
renumbering and under identity, and with a renumbered copy that adds sulfur for oxygen.
`test_reaction_span_framed_eq_under_frame` compares two spans of furan whose aromatic system
lists its atoms in different orders. `test_reaction_span_framed_eq_under_count_error` has a row
per entity kind with a remapping one entry too long.

Four properties in `tests/property/reaction/span/remap.rs`, over generated spans with two
generated remappings: a span is equal to its remapped copy under the remapping; under identity
the result is that of `framed_eq`; inverting the remapping reverses the comparison; and two
successive remappings compose. The tests of `umol-graph-ir` pass with all features.

#### S0c — Benchmark pairs and baseline **Done**

**Module:** `benches/canonicalize.rs`.

Add the pairs of the trial for molecules and for spans: an equal pair, and unequal pairs with
different entity counts, different colors, and a different stereo configuration, for each corpus
case; and the three pairs of carbon skeletons with different constitution for molecules. Each
pair is asserted against `canonical_eq` when the benchmark starts. The groups are run once and the
times recorded here as the baseline. Additive (green). [dep: S0a]

The run uses `--warm-up-time 0.5 --measurement-time 1.5`.

The pairs are in `bench_canonical_eq_unequal`, in seven groups: `canonical_eq_unequal_counts`,
`_colors`, `_constitution`, and `_stereo` under `canonicalize/complete`, and
`canonical_eq_unequal_counts`, `_colors`, and `_stereo` under
`canonicalize/reaction_span_modified`. The equal pairs are the existing `canonical_eq_equal`
groups. The stage gate passes: the tests of S0b on unchanged sources, and clippy with all
targets.

Baseline, in µs, from one run of every group whose name contains `canonical_eq`, 215 benchmarks:

```text
cargo bench -p umol-graph-ir --bench canonicalize -- \
    --warm-up-time 0.5 --measurement-time 1.5 --save-baseline pre-245 'canonical_eq'
```

Criterion keeps the run as the baseline `pre-245` in the cargo target directory.

Molecules, `canonicalize/complete`:

| Case | Equal | Counts | Colors | Stereo |
| --- | ---: | ---: | ---: | ---: |
| `ordinary_naphthalene` | 25.8 | 25.6 | 19.2 |  |
| `disconnected_rings` | 47.5 | 48.2 | 41.6 |  |
| `overlay_heavy` | 38.3 | 38.1 | 37.4 | 37.7 |
| `overlay_unconstrained` | 27.3 | 27.6 | 26.9 |  |
| `aromatic_electron_counts` | 24.3 | 24.8 | 20.8 |  |
| `symmetric_charge_sum` | 43.8 | 44.7 | 35.0 |  |
| `asymmetric_charge_sum` | 275.5 | 275.4 | 162.4 |  |
| `tetrahedral_stereo` | 16.7 | 16.8 | 16.6 | 16.6 |
| `cis_trans_stereo_bond` | 17.6 | 17.9 | 17.5 |  |
| `mixed_atom_and_bond_stereo` | 29.7 | 29.4 | 29.3 | 29.7 |
| `frame_relative_stereo_constraint` | 17.5 | 17.7 | 17.8 | 17.5 |
| `meso_dichlorobutane` | 39.5 | 39.7 | 31.2 | 36.1 |
| `para_stereo_trichloropentane` | 54.1 | 55.1 | 42.4 | 56.6 |
| `para_stereo_cascade` | 515.8 | 525.1 | 513.9 | 518.4 |
| `symmetric_stereo_triisopropanolamine` | 71.4 | 74.0 | 63.6 | 90.1 |
| `symmetric_stereo_butanol` | 115.3 | 117.5 | 115.3 | 115.3 |
| `symmetric_stereo_dimethylbutanol` | 313.0 | 320.7 | 322.7 | 319.2 |
| `feature_free_connected` | 52.9 | 53.3 | 53.2 |  |
| `feature_free_disconnected` | 62.5 | 59.5 | 53.9 |  |
| `symmetry_heavy_radicals` | 52.3 | 53.6 | 46.2 |  |

Molecules with different constitution, `canonicalize/complete/canonical_eq_unequal_constitution`:

| Case | Time |
| --- | ---: |
| `hexane_2_methylpentane` | 13.5 |
| `cyclohexane_two_cyclopropanes` | 26.2 |
| `decalin_bicyclo_5_3_0_decane` | 21.8 |

Reaction spans, `canonicalize/reaction_span_modified`:

| Case | Equal | Counts | Colors | Stereo |
| --- | ---: | ---: | ---: | ---: |
| `ordinary_naphthalene` | 18.3 | 23.8 | 17.7 |  |
| `disconnected_rings` | 40.4 | 46.9 | 35.6 |  |
| `overlay_heavy` | 44.2 | 43.5 | 44.3 | 44.2 |
| `overlay_unconstrained` | 29.8 | 30.0 | 29.5 |  |
| `aromatic_electron_counts` | 25.6 | 25.6 | 21.8 |  |
| `symmetric_charge_sum` | 28.9 | 38.0 | 22.7 |  |
| `asymmetric_charge_sum` | 55.2 | 176.5 | 35.9 |  |
| `tetrahedral_stereo` | 19.2 | 18.6 | 19.1 | 19.4 |
| `cis_trans_stereo_bond` | 20.3 | 20.1 | 20.1 |  |
| `mixed_atom_and_bond_stereo` | 33.3 | 33.1 | 33.0 | 33.9 |
| `frame_relative_stereo_constraint` | 21.3 | 20.4 | 21.2 | 21.3 |
| `meso_dichlorobutane` | 25.8 | 33.8 | 25.6 | 25.7 |
| `para_stereo_trichloropentane` | 34.7 | 46.0 | 33.6 | 35.0 |
| `para_stereo_cascade` | 305.6 | 419.8 | 304.9 | 304.9 |
| `symmetric_stereo_triisopropanolamine` | 118.7 | 98.2 | 95.0 | 118.1 |
| `symmetric_stereo_butanol` | 126.8 | 124.1 | 128.4 | 125.5 |
| `symmetric_stereo_dimethylbutanol` | 310.4 | 319.4 | 315.7 | 310.5 |
| `feature_free_connected` | 43.6 | 50.4 | 43.9 |  |
| `feature_free_disconnected` | 49.6 | 56.5 | 45.6 |  |
| `symmetry_heavy_radicals` | 43.0 | 49.4 | 37.1 |  |

Reactions, `canonicalize/reaction`, and their spans, `canonicalize/reaction_span`:

| Case | Reaction, equal | Reaction, unequal | Span, equal |
| --- | ---: | ---: | ---: |
| `ordinary_naphthalene` | 46.2 | 47.3 | 30.2 |
| `disconnected_rings` | 70.3 | 71.3 | 53.6 |
| `overlay_heavy` | 80.6 | 78.5 | 43.5 |
| `overlay_unconstrained` | 55.8 | 53.2 | 30.2 |
| `aromatic_electron_counts` | 41.4 | 41.3 | 26.7 |
| `symmetric_charge_sum` | 64.9 | 63.1 | 47.4 |
| `asymmetric_charge_sum` | 325.8 | 315.6 | 297.4 |
| `tetrahedral_stereo` | 33.0 | 32.7 | 18.6 |
| `cis_trans_stereo_bond` | 35.6 | 35.4 | 20.0 |
| `mixed_atom_and_bond_stereo` | 56.2 | 56.7 | 33.3 |
| `frame_relative_stereo_constraint` | 37.4 | 37.2 | 20.2 |
| `meso_dichlorobutane` | 60.8 | 60.2 | 43.2 |
| `para_stereo_trichloropentane` | 79.4 | 79.3 | 58.5 |
| `para_stereo_cascade` | 603.6 | 575.4 | 543.2 |
| `symmetric_stereo_triisopropanolamine` | 101.2 | 100.6 | 77.0 |
| `symmetric_stereo_butanol` | 143.2 | 141.2 | 121.3 |
| `symmetric_stereo_dimethylbutanol` | 347.6 | 344.8 | 323.8 |
| `feature_free_connected` | 74.1 | 72.9 | 58.6 |
| `feature_free_disconnected` | 79.9 | 78.7 | 67.8 |
| `symmetry_heavy_radicals` | 67.9 | 67.4 | 58.4 |

Scaling cases and inline constraints:

| Case | Scaling, equal | Scaling, unequal | Inline, molecule | Inline, span |
| --- | ---: | ---: | ---: | ---: |
| `feature_free_connected` | 53.2 | 53.0 | 54.8 | 59.3 |
| `feature_free_disconnected` | 59.8 | 59.0 | 61.1 | 65.2 |
| `symmetry_heavy_radicals` | 52.4 | 54.7 | 52.8 | 57.3 |

### S1 — Rule for contradictions and count comparison

This stage changes the result of the three equalities for contradictory aggregates and nothing
else. Each subitem carries the tests that change with it (red→green).

#### S1a — `normalized_eq` and `framed_eq` of the aggregates **Done**

**Module:** `ir/molecule.rs`, `ir/reaction_span.rs`, `ir/reaction.rs`, `ir/traits.rs`.

`Molecule`, `ReactionSpan`, and `Reaction` override `Normalize::normalized_eq` and
`Reframe::framed_eq` by the table of the Contradictions section. `Molecule` and `ReactionSpan`
compare `entity_counts` first. A reaction has no entity counts and starts with `==`. The rustdoc
of the two trait methods states the default, which forms keep, and the rule of the three
aggregates. The rustdoc of both `framed_eq_under` methods loses the sentence on totalization.
Breaking for contradictory aggregates (red→green). [dep: S0a, S0b]

**Tests:** for each aggregate and each of the two equalities, rows for: two values that pass
with equal and with unequal forms; one that fails; two that fail and are structurally equal;
two that fail with equal normal forms; two that fail otherwise; and, for molecules and spans,
different entity counts. A property per aggregate, over strategies that include contradictory
values, requires that `normalized_eq` implies `framed_eq`. Existing tests whose expected value
changes are found by the test run and shown old beside new.

No unit test stated the present rule for these two equalities. Three properties did:
`test_molecule_canonical_eq_contradiction`, `test_reaction_canonical_eq_contradiction`, and
`test_reaction_span_canonical_eq_contradiction`. Each takes three structurally different values
that fail normalization.

| Assertion | Old | New |
| --- | --- | --- |
| `normalized_eq` of each pair of the three | true | false |
| `framed_eq` of each pair of the three | true | false |
| `normalized_eq` and `framed_eq` of a value and its clone | not stated | true |

The `canonical_eq` assertions of the three properties are unchanged until S1b.

New unit tables: `test_molecule_normalized_eq_contradiction`,
`test_molecule_framed_eq_contradiction`, `test_reaction_span_normalized_eq`,
`test_reaction_span_framed_eq`, `test_reaction_normalized_eq`, and `test_reaction_framed_eq`, 36
rows in all. The rows `both_fail_equal_normal_forms` and `both_fail_different_frames` use a
cyclopropenyl cation whose aromatic system has two electron counts for three atoms, which
normalizes and fails to reframe.

No separate property for the implication was added. On contradictory values the three properties
above state the exact result of both equalities, which contains it, and on satisfiable values the
standardization properties already state it. The tests of `umol-graph-ir` pass with all features.

#### S1b — `canonical_eq` of the aggregates **Done**

**Module:** `ir/canonicalize.rs`, both submodules.

`canonical_eq` of `Molecule` and of `ReactionSpan` compares `entity_counts` first. In all three
aggregates two failed canonicalizations give `self.framed_eq(other)`. Every other pair is still
compared through canonical forms. The rustdoc of `Canonicalize::canonical_eq` states the rule.
Breaking for contradictory aggregates (red→green). [dep: S1a]

**Tests:** expected values that change:

| Test | Old | New |
| --- | --- | --- |
| `test_molecule_canonical_eq_contradiction`, one-atom against two-atom molecule | true | false |
| `test_reaction_canonical_eq_contradiction`, two reactions with different deltas | true | false |
| `test_molecule_canonical_eq_contract`, `test_reaction_span_canonical_eq_contract`, `test_reaction_canonical_eq_contract`: expected value when both fail | `true` | `left.framed_eq(&right)` |
| The same three properties: the hashes when both fail | require `expected` | require nothing |

The three `test_*_canonical_eq_contradiction` properties are read at the subitem and their
expected values shown old beside new. New rows cover the `canonical_eq` line of the table and
different entity counts. The property of S1a gains that `framed_eq` implies `canonical_eq`.

The two unit tests are now tables of the rule, and `test_reaction_span_canonical_eq_contradiction`
is new, 16 rows in all. The molecule test no longer asserts that no remapping relates its two
contradictory molecules, which illustrated the old rule. The reaction row
`both_fail_equal_reframed_forms` holds two reactions whose deltas do not materialize and whose
left sides list one aromatic system in two atom orders; they are equal.

In the three contradiction properties the `canonical_eq` of each pair of the three values
changes from true to false, and a value is asserted equal to its clone. The three contract
properties passed unchanged before the edit, since their strategies do not produce two failing
values that differ; their expected value and their hash clause are changed as the table says.
As in S1a no separate property for the implication was added. The tests of `umol-graph-ir` pass
with all features.

#### S1c — Guides **Done**

**Module:** `docs/development/data-types.md`, `docs/development/nomenclature.md`,
`docs/development/property-tests.md`.

Rewrite the passages on the totalization of equality for two intrinsic contradictions to the
rule of the Contradictions section. The property-tests guide states the old rule in its section
on the standardization pipeline. `python-api.md` and `operations.md` are read for statements of
the old rule. No code. [dep: S1b]

The data types guide states the rule after its sentence on `framed_eq_under` and in its paragraph
on the totality of equality. The nomenclature guide states it under the equality ladder, where
`framed_eq_under` is now named for `Molecule` and `ReactionSpan`. The property-tests guide states
it in the section on the standardization pipeline. `python-api.md` and `operations.md` say
nothing about the equality of contradictory values and are unchanged. The sentence of the data
types guide on edit removal, which compares stored entries with `framed_eq`, is unchanged, since
entries keep the default rule.

Two sentences of the rustdoc of the `Canonicalize` trait that spoke of failure totalization were
reworded with the guides. No behavior changed in this subitem.

Stage gate of S1: the tests of S1b on unchanged code, clippy with all targets, and the Python
tests, 3,420 passed.

### S2 — Two changes to existing canonicalization code

No result changes in this stage.

#### S2a — `incidence_color` **Done**

**Module:** both submodules.

The color of one incidence, its role and electron count, moves out of the closure in
`initial_colors` into `incidence_color`, which `initial_colors` calls. No behavior changes
(green).

**Tests:** `test_incidence_color` in both submodules, with a bond endpoint, an aromatic atom with
an electron count, and for spans an atom whose electron count differs between the sides. The
tests of `initial_colors` are unchanged.

`incidence_color(value, incidence_graph, edge)` returns the role and the electron count; the
body is the former closure, moved without change. The molecule test has seven rows: a bond
endpoint of ethane, the donor and the acceptor of ammonia borane, the oxygen and a carbon of
furan, an aromatic atom of a system without literal counts, and a boron of a three-center bond.
The span test has six rows on furan: a bond endpoint, and the oxygen of an aromatic system that
is unchanged, modified, added, and removed, and a carbon whose count is equal on both sides of a
modified system. The tests of `umol-graph-ir` pass with all features.

#### S2b — `remapping_preserves_stereo` **Done**

**Module:** both submodules.

`generator_preserves_stereo` takes a second value as the target and is renamed
`remapping_preserves_stereo`. `retain_stereo_preserving_generators` passes the same value twice.
No behavior of the search changes (green).

**Tests:** `test_remapping_preserves_stereo` in both submodules: a value against itself under a
generator that keeps and one that breaks a configuration; a value against a renumbered copy under
the remapping of the renumbering; and against a copy with one configuration inverted. The tests
of `retain_stereo_preserving_generators` are unchanged.

`remapping_preserves_stereo(source, target, incidence_graph, images)` reads the image entity
from `target`; the comparison itself is unchanged. The loop variables are renamed, since `source`
and `target` now name the two values. The molecule test has seven rows: chlorotrifluoromethane
under a rotation of its three fluorine atoms and under a transposition of two;
bromochlorofluoroiodomethane under a remapping that moves the site, against a renumbered copy,
and against a renumbered copy with the other configuration; and but-2-ene against itself and
against the other configuration. The span test has the same seven rows on unchanged spans and
two rows on a stereo atom whose configuration is modified, in the same and in the reversed
direction. The tests of `umol-graph-ir` pass with all features, and clippy with all targets is
clean, which is the stage gate of S2.

### S3 — Remapping search for molecules

For molecules without contradiction, `canonical_eq` relates the same pairs after this stage as
before it.

#### S3a — Shared functions

**Module:** `ir/canonicalize.rs`.

Add `color_representatives`, `same_cells`, `remapping_preserves_incidences`, and
`remapping_exists`. `remapping_exists` covers steps 4 to 7 and calls the search through
`canonical_search_impl`. Additive. The functions have no caller outside tests until S3b, so the
dead-code warnings of this subitem end there. [dep: S2b]

**Tests:** one test per function. `test_color_representatives` takes pairs of integer lists.
`test_same_cells` and `test_remapping_preserves_incidences` take two molecules at Constitution;
the second also compares its result with `framed_eq_under` under the same remapping.
`test_remapping_exists` has a row for each way a pair is decided: cells differ, single-node cells
with a remapping that passes and one that fails, the first path, a search that finds a remapping,
and a search that finds none.

#### S3b — `Molecule::canonical_eq_at`

**Module:** `ir/canonicalize/molecule.rs`.

Add `same_entity_color`, `same_colors`, and the method `canonical_eq_at`. `canonical_eq` calls it
when both molecules select the same description level and compares canonical forms only when they
do not. Changes no result (green). [dep: S1b, S2a, S3a]

**Tests:** `test_same_entity_color` and `test_same_colors`. `test_molecule_canonical_eq` gains
rows that reach each step, with the compounds of the trial: different colors; hexane and
2-methylpentane; bromochlorofluoroiodomethane and its mirror image; the naphthalene skeleton
renumbered; meso-2,3-dichlorobutane renumbered and with one center inverted; cyclohexane and two
cyclopropanes; and a ring with a charge-sum constraint.
`test_molecule_canonical_eq_remapping_witness` gains the rows small enough for its enumeration of
remappings. The properties of S1b hold unchanged.

### S4 — Remapping search for reaction spans and reactions

#### S4a — `ReactionSpan::canonical_eq_at`

**Module:** `ir/canonicalize/reaction_span.rs`.

The span counterpart of S3b: `same_entity_color`, `same_colors`, `canonical_eq_at`, and the call
from `canonical_eq`. At Full the check is `framed_eq_under` of S0b. Changes no result (green).
[dep: S0b, S1b, S2a, S2b, S3a]

**Tests:** the span counterparts of the tests of S3b, on the span cases of the trial. A test
enumerates the remappings of small spans with `framed_eq_under`, as
`test_molecule_canonical_eq_remapping_witness` does for molecules.

#### S4b — `Reaction::canonical_eq` through reaction spans

**Module:** `ir/canonicalize/reaction_span.rs`.

`Reaction::canonical_eq` converts both reactions to reaction spans and calls `canonical_eq` of
the spans. Two failed conversions give `self.framed_eq(other)`; one gives false. Two spans that
convert and both fail are decided by the span. Changes no result (green). [dep: S4a]

**Tests:** `test_reaction_canonical_eq` and the properties of S1b hold unchanged. If the
property that `framed_eq` implies `canonical_eq` fails for reactions, work stops and the case is
raised.

### S5 — Final comparison and closeout

#### S5a — Final comparison

Run the groups of S0c with the same times and `--baseline pre-245`, and record them beside the
baseline of S0c. A pair that costs more than its baseline is raised before closeout. [dep: S4b]

#### S5b — Closeout

**Module:** `docs/development/data-types.md`, this document, `discussion/000-status.md`.

The guide describes `canonical_eq` by the pairs it relates and no longer by the comparison of
canonical forms. Run the workspace tests with all features, the Python tests, and clippy on the
workspace. Set the status of this document and its index row. [dep: S5a]

### Stage gate

After every stage: `cargo test -p umol-graph-ir --all-features` and
`cargo clippy -p umol-graph-ir --all-targets -- -D warnings`. S1 also runs the Python tests, since
it changes results that Python reaches.

For a subitem that changes no behavior, one build with the tests of `umol-graph-ir` is the whole
check. No benchmark runs outside S0c and S5a.

### Critical path

S0 → S1 → S2 → S3 → S4 → S5. S1 settles which pairs are equal before S3 and S4 change how the
answer is found, so the tests of S1 guard those two stages. No stage is deferrable.
