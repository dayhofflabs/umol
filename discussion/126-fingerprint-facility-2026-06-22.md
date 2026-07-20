# 126 — Fingerprint facility: design

**Competed as of 2026-06-23**: WL, ECFP, fold→`BitFp`, Morgan (RDKit-exact),
substructure direction A (unhashed) and direction B (RDKit `PatternFingerprint`
replica for acyclic + single rings), the counts layer (`CountedFeatureSet` /
`SignedFeatureSet` / `featurize_counted`), and reaction FPs (`Difference` /
`DisjointUnion`). Deferred with their blockers tracked elsewhere: slice-5 B
ring-junction templates (doc 128 — derived-predicate matching), slice 7 DRFP and
slice 8 BRIDGIT (both need SMILES rendering / canonicalization).

## Goals

A versatile fingerprinting facility for molecules and reactions serving three
purposes:

1. **Deduplication / similarity** — uniqueness checks and similarity (Tanimoto,
   Dice, …) over molecules and reactions.
2. **Feature representation** — fixed-width count/bit vectors as ML features.
3. **Prescreening** — bit-fingerprint prefilters for substructure library screens
   (the role RDKit's pattern fingerprint plays today).

Not boiling the ocean: the *framework* must accommodate a broad fingerprint
landscape, but the concrete deliverables now are:

- **a. ECFP** (Rogers & Hahn 2010).
- **b. Morgan** as RDKit currently produces it — pinned to a named RDKit revision,
  replicated and frozen (RDKit's output has not been stable across versions).
- **c. BRIDGIT** — circular fingerprints over substructure-matched reactive sites
  (Hadadi et al. 2019).
- **d. DRFP** — differential reaction fingerprint (Probst 2022).
- **e. Stable WL** — built on the existing `umol-graph-core` refinement.
- **f. Substructure (unhashed)** — explicit subgraph keys for exact prescreening
  (RDKit pattern-fingerprint role).

## What already exists (build on, don't rebuild)

- **`umol-graph-core::algorithms::refine`** — 1-WL color refinement with a
  pluggable, frozen `RefinementHash` (`RefinementXxh3Scheme`, named variants
  `albatross`/`bullfinch`). It already exposes `coloring_at(round)`,
  `features()`, `counts()`, `graph_hash()` — i.e. per-round colorings and an
  order-independent digest. Its own doc states "invariants/folding/dedup live in
  the caller", which is exactly the fingerprint layer. **This is the engine for
  circular (ECFP/Morgan/WL) fingerprints.**
- **`umol-graph-core::algorithms`** also has `auto` (automorphism / exact
  canonical form), `cycles`, `coloring`, `subiso`, `mcs`. Canonicalization feeds
  DRFP shingle keys and ECFP duplicate removal.
- **`umol-ast::ast::substructure`** — `substructure_matches(host, strategy,
  subiso)`. The substrate for BRIDGIT (reactive-site environments) and the
  substructure fingerprint.
- **`umol-ast::ast::reaction`** — a reaction type already exists (combinator
  substrate for reaction fingerprints).
- **Atom/bond invariant accessors** on the views: `element`, `isotope_mass`,
  `charge`, `valence`, `degree`, `total_hydrogens`, `is_in_ring`, ring perception
  (`rings_cache`), bond order/aromatic. The Daylight ECFP seed invariants are all
  computable.
- The legacy `morgan` impl is already gone; only the legacy-gated
  `umol-graph/benches/morgan.rs` remains, to be replaced.

## Landscape review — omissions worth accommodating

From the deliverable list plus a materials and literature/web review, the families
the framework should be able to express (most not built now, but the API must not
preclude them):

- **Path / subgraph topological** (RDKit `RDKitFP`, Daylight) — hashed enumerated
  paths/subgraphs; distinct from circular and from the unhashed substructure FP.
- **Atom-pair (AP)** and **topological-torsion (TT)** (RDKit) — standard for
  similarity and ML; atom-typed pairs by topological distance / 4-atom paths.
- **Functional-class circular (FCFP)** — ECFP with functional-class atom
  invariants (donor/acceptor/aromatic/halogen/…) instead of Daylight invariants.
- **Molprint2D/3D** — Bender's circular fingerprints (sibling of ECFP/BRIDGIT).
- **MACCS / structural keys** — fixed dictionary of SMARTS patterns; a *keyed*
  substructure fingerprint (generalizes deliverable f to a named dictionary).
- **MinHash / LSH variants** — MHFP6 and **MAP4** (Probst & Reymond; same author
  as DRFP). MinHash of circular/atom-pair shingles; the right tool for purpose 3
  (similarity) at very large scale. Notable given the reaction-network scale.
- **Pharmacophore 2-/3-point** — feature-pair/triplet fingerprints.
- **Reaction**: **difference fingerprint** (RDKit; signed product − reactant), and
  **Condensed Graph of Reaction (CGR)** as a reaction→graph transform that lets
  *any* molecular fingerprint apply to a reaction.

These cluster into a few feature *shapes* (circular environments, paths/pairs,
keyed patterns, MinHash-of-shingles) over a few *encodings* (sparse multiset, bit,
count) with optional *reaction combinators* — which is what the architecture below
is organized around.

## Proposed architecture

Two output flows over one set of primitives, not a single fixed pipeline:

```
                              atom + bond invariant (circular seed; enumerated, no Custom)
                                        │
MoleculeAst / ReactionAst ──► Featurizer (named algorithm enum)
                                        │
        ┌───────────────────────────────┴───────────────────────────────┐
   length-agnostic                                              width-bearing
   → FeatureSet  ──► [ReactionCombinator] ──► Encoder(width) ──► BitFp/CountFp
   (raw identifiers)                          (fold primitive)
                                                                → folds inline (width is
                                                                  a Featurizer param), emits
                                                                  BitFp directly (e.g. ECFP)

ops over any output: tanimoto / dice / cosine (similarity); is_subset (screening)
```

### Components

1. **Invariants (circular seed), enums with explicit variants only — no `Custom`.**
   On-the-fly invariant flexibility trades reproducibility for unclear gain; if a
   need appears we add a named variant.
   - `AtomInvariant`: `Ecfp` (Daylight: heavy-degree, element Z, charge,
     attached-H count, ring membership, isotope mass — Rogers & Hahn Table 1),
     `Fcfp` (functional class), `Morgan` (RDKit's set).
   - `BondInvariant`: `Order`, `OrderAromatic`, `OrderAromaticRing`.
   Invariants are *inputs* to the circular featurizers (the radius-0 identifier),
   not featurizers themselves.

2. **`Featurizer` — the algorithm, an enum (explicit dispatch, per the
   algorithm-transparency convention; no trait-object zoo).** Substructure-based
   methods are simply another family of variant.
   - `Ecfp { radius, atom_invariant, bond_invariant, width: Option<usize>, hash }`
     — width-bearing (`Some` → folds inline; `None` → sparse `FeatureSet`).
   - `Morgan { radius, …, rdkit_revision }` — pinned RDKit revision.
   - `Wl { rounds, scheme }` — thin wrapper over `refine` (length-agnostic).
   - `Drfp { radius, … }` — reaction-native (length-agnostic shingles).
   - `Bridgit { radius, reactive_site_query, … }`.
   - `Substructure { enumeration }` — unhashed structural keys.
   - (future: `AtomPair`, `TopologicalTorsion`, `Path`, `Maccs`, `Mhfp`, `Map4`.)

   **Output is heterogeneous:** a length-agnostic featurizer yields a `FeatureSet`;
   a width-bearing one (ECFP with width) yields a folded `BitFp`/`CountFp`
   directly. `featurize(...) -> FeaturizerOutput { Sparse(FeatureSet) | Bits(BitFp) }`.

   **Combining featurizers** is common (e.g. ECFP ⊕ atom-pair as ML features). It
   lives *above* a single `Featurizer`: a `&[Featurizer]` whose outputs are merged.
   Two merge modes — feature-level tagged union (each featurizer's identifiers
   namespaced by its identity to avoid cross-family collision, then encoded as one
   vector) for similarity/dedup, and encoded-vector concatenation for ML features.

3. **`FeatureSet` — the length-agnostic sparse intermediate.** Concrete
   representation: a **sorted `Vec<Id>`** of the distinct identifiers (the binary
   set), plus a parallel sorted `Vec<(Id, u32)>` of counts when needed (RDKit
   `SparseIntVect` equivalent). `Id` is the hash width (`u32`/`u64`/`u128`).
   Tanimoto / Dice / cosine / union / intersection / subset are linear merges over
   the sorted vectors — one allocation, cache-friendly.

   *Not roaring.* Roaring compresses *large, clustered* integer sets — its
   `2^16`-chunk array/bitmap/run containers assume locality. Fingerprint identifiers
   are frozen hashes spread uniformly over the whole `2^width` space in sets of only
   tens–hundreds, so there is no clustering: every container degenerates to one
   element and you pay per-element structure overhead with zero compression (and
   `u64` via `RoaringTreemap = BTreeMap<u32, RoaringBitmap>` allocates ~one bitmap
   per element). A sorted `Vec` wins on both memory and speed here. Roaring is kept
   only as a deferred option for a *large, sparse, folded* fingerprint (a dense,
   small universe) — which none of the a–f deliverables need.

   The `identifier` is a frozen hash for hashed fingerprints, or a canonical
   structural key for the unhashed substructure fingerprint. Per-feature provenance
   (covered atoms, radius) is available from the refinement loop (see Primitives) for
   bit-info and ECFP dedup.

4. **`Encoder` (length-dependent) → `BitFp<N>` / `CountFp<N>`.** Applies the shared
   **fold primitive** (`identifier mod N`; binary OR for bits, sum for counts).
   Folding collisions are inherent and lossy; an optional bit-info map records which
   identifiers landed on each index. The *same* fold primitive is used inline by
   width-bearing featurizers — so the linear (`Featurizer → FeatureSet → Encoder`)
   and iterative (fold-as-you-refine) flows share one implementation, not two.

5. **Operations:** `tanimoto`/`dice`/`cosine` over `BitFp`/`CountFp`/`FeatureSet`;
   `is_subset(query, target)` for prescreening (query bits ⊆ target bits).

6. **Reactions:** `ReactionCombinator` over per-component `FeatureSet`s:
   `Difference` (RDKit signed), `SymmetricDifference` (DRFP), `RoleTagged`. DRFP is
   reaction-native (shingles from both sides, symmetric difference, no atom map);
   `Difference`/`RoleTagged` compose any molecular featurizer and use the atom map.
   See "Reaction model" below for the semantics assumed now.

### Composable primitives (not one-off impls per algorithm)

The featurizers are compositions of a small, well-defined primitive set:

1. **Provenance-aware circular refinement** — per `(atom, round)` yields an
   identifier *and* the covered atom set. Feeds ECFP / Morgan / FCFP / BRIDGIT.
2. **Fold** — `identifier → index` with OR/sum accumulation. Feeds the Encoder and
   the width-bearing featurizers.
3. **Canonical subgraph key** — via `auto` (canonical form). Feeds DRFP shingles,
   the unhashed substructure keys, and any unlabeled-graph features.
4. **Substructure enumeration / matching** — feeds the substructure fingerprint and
   BRIDGIT reactive sites.
5. **Sparse set (sorted `Vec<Id>`) + similarity/subset ops** — linear-merge
   Tanimoto / Dice / cosine / subset over sorted identifier vectors; shared by all
   outputs.

`refine` already provides per-round colorings and a frozen hash; primitive 1 needs
it to *also* expose the covered atom set (the BFS frontier accumulated per round).
Preference is to extend `refine` with that (a composable addition) rather than
write a fingerprint-private refinement; if the featurization loop's needs exceed
what fits cleanly there, a thin shared circular-refinement primitive lives in the
fingerprint module — still shared, never per-algorithm.

### Additional primitives for the broader landscape

A rough pass over the landscape-review methods (grounded against the molintern
Python reference — `featurizer.py`, `cgr_graph.py`) surfaces five more shared
primitives plus two generalizations; the framework needs these to reach the wider
list without per-method code:

6. **Atom/bond typing (invariants)** — typed labels from topology (Daylight,
   functional-class, AP/TT atom-type, pharmacophore feature class). Already the
   `AtomInvariant`/`BondInvariant` config; named as a primitive because AP/TT and
   pharmacophore reuse it *outside* the circular seed.
7. **Pairwise topological distances** — bounded all-pairs BFS. Feeds atom-pair,
   MAP4, and pharmacophore feature pairs/triples.
8. **Bounded subgraph / path enumeration** — enumerate connected subgraphs / linear
   paths up to size *L*. Distinct from primitive 4 (matching a *given* pattern via
   subiso): this *generates* all bounded subgraphs. Feeds the RDKit path/subgraph
   FP, topological-torsion (length-4 paths), the unhashed substructure FP (f), and
   the radius-bounded shingles of DRFP/MHFP (an atom environment is a radius-bounded
   subgraph).
9. **MinHash / LSH** — a MinHash signature of a feature set (+ optional LSH
   bucketing). Feeds MHFP6 and MAP4; an *alternative encoder* to fold for
   large-scale similarity.
10. **CGR construction** — reaction → condensed graph with dynamic (Δ bond-order)
    edge labels, keyed by the atom map (stereo-aware; cf. `cgr_graph.py`). Lets any
    *molecular* featurizer apply to a reaction.

Two generalizations this implies:

- **Encoder is `Fold(N)` *or* `MinHash(k)`** — the width-dependent step has two
  members, not one.
- **Reactions have two handling modes**: (a) combine per-component feature sets
  (`Difference`/`SymmetricDifference`/`RoleTagged`); (b) **transform then
  featurize** via CGR (primitive 10). DRFP is mode (a) over shingles; CGR-Morgan is
  mode (b).

One dependency (not a fingerprint primitive): **SMARTS query patterns** — MACCS's
166-key dictionary, pharmacophore feature definitions, and richer reactive-site
queries need atom/bond query expressivity beyond `MoleculeAst`-as-pattern. Those
keyed/feature fingerprints wait on that, as richer reaction fingerprints wait on
the `reaction.rs` redesign.

### Coverage (method → primitives)

| method | primitives |
|---|---|
| ECFP / Morgan / FCFP / Molprint | 6, 1, 2, 5 |
| WL | `refine`, 2, 5 |
| DRFP | 8 (radius env), 3 (canonical SMILES), sym-diff, 2, 5 |
| BRIDGIT | 4, 1, 5 |
| Substructure (f, unhashed) | 8, 3, 5 |
| Atom-pair | 6, 7, 2, 5 |
| Topological-torsion | 6, 8 (len-4 paths), 2, 5 |
| RDKit path / subgraph | 8, 3/hash, 2, 5 |
| MACCS / keyed | 4 + SMARTS dict, 5 |
| MHFP6 | 8, 9, 5 |
| MAP4 | 6, 7, 8, 9, 5 |
| pharmacophore 2/3-pt | 6 (SMARTS feature), 7, tuple-enum, 2, 5 |
| CGR reaction FP | 10 + any molecular featurizer |
| reaction difference FP | any molecular featurizer + `Difference` |

Beyond the deliverables, the only genuinely new machinery is primitives 6–10, the
`MinHash` encoder, the CGR transform, and a SMARTS query layer — each additive,
none forcing a change to the core (circular refinement, fold, canonical key,
roaring set).

### How each deliverable maps

| deliverable | featurizer | primitives | notes |
|---|---|---|---|
| a. ECFP | `Ecfp` (width-bearing or sparse) | 1, 2, 5 | dedup per Rogers & Hahn (see below) |
| b. Morgan | `Morgan` | 1, 2, 5 | pin a revision, verify vs RDKit, freeze |
| c. BRIDGIT | `Bridgit` | 1, 4, 5 | reactive-site environments |
| d. DRFP | `Drfp` + `SymmetricDifference` | 3, 5 | reaction shingles; no atom map |
| e. WL | `Wl` | `refine` directly, 2, 5 | already frozen schemes |
| f. Substructure | `Substructure` (unhashed) | 3/4, 5 | exact `is_subset` screening |

### ECFP duplicate removal — iterative refinement carries the provenance

Rogers & Hahn remove structurally-redundant features (two features whose covered
atom set is identical, tie-broken by lower iteration then identifier). Agreed that
iterative refinement is the simpler path — and it gives the provenance for free:
during the circular refinement, an atom's covered atom set at round *r* is exactly
the accumulated BFS frontier, which the loop already walks. So track
`(identifier, covered_set)` per `(atom, round)` and apply the atom-set dedup
*inline*, rather than as a separate provenance-bookkeeping pass. This is primitive 1
above; it is also what produces RDKit-style bit-info. (This is the concrete reason
primitive 1 needs `refine` to expose the frontier.)

## Reaction model (interim, with a redesign dependency)

`umol-ast::ast::reaction` is underdeveloped: `ReactionRuleAst` (lhs/rhs +
`Vec<(AtomId, AtomId)>` atom map) is framed as a template only, but since
`MoleculeAst` is homoiconic it can hold both a concrete reaction and a rule. A
proper treatment — rename to `ReactionAst`; a partial, possibly non-injective
assignment (`Vec<(Option<id>, Option<id>)>` to cover SMIRKS *and* RDKit reaction
SMARTS, which allows non-injective templates); stereo overlays and overlay handling
generally — is a separate, larger design task and warrants its own doc.

For fingerprints now we assume **SMIRKS semantics**: a partial atom map, injective
in both directions. That suffices to define reaction handling:
- **DRFP** needs no atom map (symmetric difference of both sides' shingles).
- **Difference / RoleTagged** use the partial injective map.
- **BRIDGIT** keys off reactive sites (substructure matches), within the same map
  semantics.
Richer reaction fingerprints (non-injective templates, stereo-aware reaction
features) wait on the `reaction.rs` redesign.

## Substructure fingerprint — labeled vs. unlabeled

Deliverable f definitely needs **labeled `MoleculeAst` substructure features**:
enumerated labeled subgraphs (paths / rings / atom environments) or matches of a
pattern set, reduced to canonical structural keys, with `is_subset` for exact
prescreening (soundness requires every query feature to be derivable in the target,
which constrains the enumeration).

Open: whether we also need **unlabeled-graph features** (topology only, atom/bond
types stripped). Not required by deliverables a–f; the "how" is the same canonical
subgraph-key primitive (primitive 3) with labels removed (uniform coloring). Defer
unless a concrete use (e.g. topology-only screening, or reaction-template graph
matching) is identified. The exact enumeration set is pending a review of the
literature (`materials/subgraphs`: RDKit substructure search, pattern-only
heuristics) and existing impls (`materials/codes`, esp. RDKit's pattern FP).

## Layering / placement

- **Pure-graph** parts (WL refinement, circular recoloring, canonical keys) stay in
  `umol-graph-core` (`refine` already there).
- **Chemical** parts (atom/bond invariants, reactions, substructure-driven and
  reaction fingerprints, encodings, similarity) live in a new top-level
  **`umol-graph::fingerprint`** module (`umol-graph/src/fingerprint.rs` +
  `fingerprint/` submodules; no `mod.rs`).
- Replace `umol-graph/benches/morgan.rs` with `umol-graph/benches/fingerprint.rs`
  covering at least ECFP/Morgan over the conformance corpus.

## Reproducibility / freezing

- Hashing goes through frozen, named schemes (the `RefinementXxh3Scheme` pattern),
  never ambient `Hash`/`RandomState`, so results are bit-stable across machines and
  runs.
- **Morgan**: pin an RDKit revision, build offline fixtures (the subiso-conformance
  pattern: stored expected outputs, no FFI at test time), replicate, then freeze.
  If machine/salt dependence is found in the reference, that is itself the finding
  — freeze our own scheme and document the divergence.
- WL is already frozen; ECFP/DRFP/BRIDGIT each get a named scheme + fixtures.

## Resolved (this round)

- **Naming**: `Featurizer` (not "Extractor"); `Featurizer` is the algorithm enum;
  no `Custom` invariants — enumerated only.
- **Output is heterogeneous**: length-agnostic featurizers → `FeatureSet`;
  width-bearing (ECFP with width) → folded `BitFp`/`CountFp` inline. Encoders are
  length-dependent; one shared fold primitive serves both flows.
- **Non-folded representation**: a **sorted `Vec<Id>`** of identifiers + a parallel
  sorted `Vec<(Id, u32)>` for counts. *Not roaring* — hash identifiers are uniform
  over `2^width`, defeating roaring's clustering assumption (one degenerate container
  per element); a sorted `Vec` is smaller and faster. Roaring deferred to a possible
  large-sparse-*folded* FP only.
- **Width — fold vs. hash**: the unfolded `FeatureSet` has **no fixed width**;
  length is the feature count. A fold width `N'` is introduced *only* at the
  `Encoder` stage when a fixed-length dense vector is required (ML features,
  purpose 2); dedup and similarity (purposes 1, 3) operate on the unfolded set and
  never fold. The only width that matters there is the **hash bit-length**, and only
  for collisions — see open question 1.
- **ECFP dedup**: iterative refinement with inline provenance (covered set = BFS
  frontier); no separate provenance pass.
- **Reaction (interim)**: SMIRKS semantics (partial, doubly-injective map); richer
  templates/stereo wait on the `reaction.rs` redesign.
- **Composability**: featurizers compose the shared primitive set (5 core + 6–10
  for the broader landscape); no per-algorithm one-offs.
- **Provenance source**: adjusting `refine` to expose the per-round
  covered-atom-set frontier, or adding a separate richer-output function, are both
  acceptable (decided at implementation).

## Open questions

1. **Identifier (hash) width.** The unfolded `Vec` has no *fold*-width parameter;
   the live choice is the *hash* bit-length, which governs collisions (birthday
   bound ≈ `2^(width/2)` distinct items: `u32` ~65k, `u64` ~4.3e9, `u128` ~1.8e19).
   Two regimes: **per-feature ids** inside one fingerprint tolerate collisions
   (slight similarity noise) → `u64` ample, `u32` fine for Morgan/RDKit parity;
   **`graph_hash` for dedup** must near-uniquely identify a molecule across the whole
   reaction network → at 1e8–1e9+ nodes `u64` carries single-digit-% collision risk,
   so **`u128` for dedup hashes**. Confirm per algorithm.
2. **Substructure-FP enumeration set** + whether **unlabeled-graph features** are
   needed (and for what). Pending the literature/impl review noted above.
3. **Invariant exactness.** The precise ECFP Daylight tuple, RDKit Morgan set, FCFP
   classes; whether chirality/stereo participates at the seed stage.
4. **`reaction.rs` redesign** (rename → `ReactionAst`, `Option/Option` possibly
   non-injective assignment, stereo + overlay semantics) — its own design doc;
   prerequisite for reaction fingerprints beyond the SMIRKS-semantics interim.

## Implementation plan — slice 1: `FeatureSet` substrate + WL (built)

The thinnest end-to-end vertical: the unfolded binary substrate plus the one
featurizer that wraps `refine` directly. Validates `FeatureSet` + similarity ops
with a real featurizer before ECFP adds the invariant tuple, width/fold, and
`refine` provenance.

### Excluded from this slice

- **No folding / `BitFp` / `CountFp` / `Encoder`.** WL is length-agnostic → it emits
  a `FeatureSet` only. Folding enters with the first width-bearing featurizer (ECFP)
  and the ML-feature encoder.
- **No `FeaturizerOutput` enum.** One output shape (sparse), so `featurize` returns
  `FeatureSet` directly; the `Sparse | Bits` enum lands when `BitFp` does.
- **No counts / cosine.** `FeatureSet` is binary (presence only). A count
  fingerprint (multiset) is a separate representation, added when a count method
  (cosine, count-Tanimoto, ML vectors) actually needs it.
- **No reaction path, no benches.** Benches (`fingerprint.rs` replacing `morgan.rs`)
  come with ECFP/Morgan.

### Module layout (`umol-graph`, no `mod.rs`)

Modeled on `ValenceResolver` (enum-of-structs, no trait — featurizers are far more
alike than transformers, so the trait buys nothing):

- `src/fingerprint.rs` — re-exports only (`FeatureSet`, `Featurizer`,
  `FingerprintError`, `WlFeaturizer`).
- `src/fingerprint/feature_set.rs` — `FeatureSet<Id>` + similarity methods.
- `src/fingerprint/featurizer.rs` — `Featurizer` enum (`Wl(WlFeaturizer)`) with the
  `featurize` *method* (ground gate + dispatch); `FingerprintError`.
- `src/fingerprint/wl.rs` — `WlFeaturizer` struct with its own `featurize` method.
- `src/lib.rs` — `pub mod fingerprint;`.

No new crate dependencies.

### `FeatureSet<Id>`

- Generic over `Id: Copy + Ord` (zero-cost; WL → `u64`, Morgan `u32` / dedup `u128`
  reuse it). A single field `ids: Vec<Id>`, sorted and duplicate-free.
- **Sorted `Vec`, not a map**: identifiers are uniform hashes, so the ops are merges
  and a contiguous sorted array is cache-optimal (a `BTreeMap`/`HashMap` adds
  tree/hash indirection over the same scan). And nothing builds a map — `refine`'s
  `features()` is already sorted+unique, so the WL path wraps it zero-copy via
  `from_sorted_unique`; `from_features` sorts+dedups arbitrary input in place.
- Methods: `tanimoto` / `dice` by linear merge of the two sorted runs (optimal for
  comparable sizes); `is_subset(query, target)` binary-searches the query into the
  shrinking tail of the target — O(|q|·log|t|), which beats a merge when the query
  is much smaller (the prescreening regime). At-scale prescreening across many
  targets wants an inverted index *above* this type — later.

### WL featurizer

- **Precondition: ground molecule.** `Featurizer::featurize` checks `mol.is_ground()`
  once (shared by every featurizer) and returns `Err(FingerprintError::NotGround)`
  otherwise — no coercion. Groundness is the molecule-level contract; "resolved" is
  only one path to it, hence `NotGround`, not tied to resolution. Given that,
  `WlFeaturizer::featurize` is infallible and its seeds read concrete literals.
- Operates on the **raw atom graph** (`mol.raw_graph()`: atoms = nodes, localized
  bonds = edges). Overlay relations do not enter the topology; bond aromaticity is
  only a seed field. Overlay-aware (incidence-graph) WL is a later config.
- Seeds (bit-packed `u64`, distinct tuples stay distinct before the scheme rehashes):
  atom = (`atomic_number`, `charge`<<8, `implicit_hydrogens`<<24); bond =
  (`order`, `is_in_aromatic_system`<<16). Concrete values via `.as_lit()`.
- Calls `raw_graph().refine(node_label, edge_label,
  WeisfeilerLehman { rounds, scheme })` under the frozen
  `RefinementXxh3Scheme<RefinementWidth64>` (`albatross`), and wraps
  `Refinement::features()` (sorted+unique over all rounds) into `FeatureSet<u64>`.
- **Not** seeded from `ConstitutionColoring`: its `color` digests via std
  `DefaultHasher` (SipHash), deterministic within a run but **not stable across Rust
  versions** — disqualifying for a frozen fingerprint. The seed packs raw field
  *values* (stable) and lets the frozen scheme do all hashing.

### Tests

- `feature_set.rs`: exact-scalar `tanimoto`/`dice` and `is_subset` over small literal
  sets; `from_features` (sort+dedup) and `from_sorted_unique` asserting exact ids.
- `wl.rs`: a pinned exact-id fixture for ethane under `albatross` (reproducibility
  anchor); an `(a, b, equal)` table covering order-independence (relabeled propane →
  equal) and discrimination (ethane≠propane, propane≠isopropyl cation). `#[fixture]`
  for the shared `WlFeaturizer`.
- `featurizer.rs`: dispatch forwards to the inner struct; non-ground → `NotGround`.

## Implementation plan — slice 2: ECFP (Rogers & Hahn 2010)

Faithful to the R&H *algorithm*, not RDKit. The paper specifies **no** hash function
or seed — only "a reasonable hash, uniform into 2³² (or 2⁶⁴)" (p.746). So we use
`xxh3_64` with a dedicated ECFP seed; that is a free, frozen choice, not the WL
`refine`/albatross machinery. Graph algorithms live in `umol-graph-core`; chemistry
and dedup in the domain crate.

### `umol-graph-core` additions

- `algorithms/bfs.rs`: `Graph::bfs_distances(source, max_depth) -> Vec<u32>`
  (`u32::MAX` = beyond `max_depth` or unreachable). General-purpose.
- `algorithms/circular.rs`: the R&H circular refinement,
  `Graph::circular_refine(node_label, edge_label, radius, seed) -> Vec<Vec<u64>>`
  (`[round][node]` identifier). Round 0 = `xxh3_64(seed, node_label(n))`; round *r* =
  `xxh3_64(seed, [r, prev_id, sorted (edge_label, neighbor_prev_id)…])` — the salted
  array of p.744. Generic over the label closures and seed (so FCFP/Morgan reuse it
  with a different seed/invariant); separate from WL `refine`.

### `umol-graph::fingerprint` additions

- `EcfpFeaturizer { radius, seed }` → `Featurizer::Ecfp(..)` → binary `FeatureSet<u64>`.
- **Daylight seed invariant** (paper-exact): `heavy_atom_degree`, `heavy_atom_valence`,
  `atomic_number`, `isotope_mass`, `charge`, `total_hydrogens`, `is_in_ring`, bit-packed
  into the node label. Bond label = (order, aromatic).
- Per (atom, round 0..=radius): identifier from `circular_refine`; covered **bond set**
  (rounds ≥1) from `bfs_distances` — edges with an endpoint at distance ≤ *r*−1.
- **Dedup** (p.745): round-0 identifiers go straight into the set (the set is the
  identifier-dedup). Rounds ≥1 get structural dedup by bond set — keep the feature
  with the smallest `(round, identifier)` per bond set (R&H rules 1 & 2). The result
  is round-0 ids ∪ survivors.
- ECFP_{2·radius} naming; no stereo at seed; no fold/counts (slice 3).

### Tests

- **R&H conformance, hash-independent**: butyramide `CCCC(=O)N` feature *counts* at
  radius 0/1/2/3 = **5 / 11 / 14 / 14** (paper Fig 8: ECFP_0/2/4/6) — validates invariants +
  iteration + dedup without depending on the hash.
- Pinned exact-id anchor under the ECFP seed.
- `bfs_distances` unit tests; order-independence + discrimination, as for WL.

## Slice roadmap (deliverables a–f)

| Slice | Deliverable | Status |
|---|---|---|
| 1 | e. Stable WL | **built** |
| 2 | a. ECFP | **built** |
| 3 | (shared) fold → `BitFp` | **built** |
| 4 | b. Morgan (frozen RDKit replica) | designed below |
| 5 | f. Substructure | A (unhashed) **built**; B (RDKit replica) **built** for acyclic + single rings (6 fixtures bit-exact); junction templates blocked on ring-bond matching |
| 6 | reaction FPs from ECFP/Morgan (Difference / DisjointUnion) | **built** (counts pulled in) |
| 7 | d. DRFP | designed; **blocked on canonical-SMILES rendering** (deferred) |
| 8 | c. BRIDGIT | designed; deferred (needs SMILES rendering + reactive-site queries) |

Order rationale: the shared fold (3) gives the already-built circular FPs a
fixed-width output; Morgan (4) reuses the ECFP machinery; Substructure (5)
introduces bounded subgraph enumeration + canonical keys (reused by DRFP);
**reactions enter at 6** with the ECFP/Morgan-based reaction FPs, which establish
the reaction-combinator framework; DRFP (7) reuses that framework with
reaction-native shingles + symmetric difference; BRIDGIT (8) builds on substructure
matching + reactions. Each slice is detailed (signatures, decisions) at
implementation time, as slices 1–2 were.

### Slice 3 — fold → `BitFp` (shared) — built

The width-bearing output: fold a `FeatureSet` to a fixed-width bit fingerprint by
`id % width`, OR-accumulating. Applies to WL/ECFP/Morgan; serves purpose 2
(fixed-width ML features) and fixed-width prescreening.

- Built: `BitFp` (runtime width, backed by `bitvec::BitVec<u64, Lsb0>`) in
  `fingerprint/bit_fp.rs`; `FeatureSet::<u32/u64/u128>::fold(width) -> BitFp` (a
  method, not an `Encoder` enum); `width`/`get`/`count_ones`/`tanimoto`/`dice`/`is_subset`
  over `BitFp`. Similarity ops run on `as_raw_slice()` (word-level popcount); `repeat`
  zeros the buffer and `set` only touches live bits, so padding stays zero and the
  raw-word ops are exact. Equal width asserted.
- Storage: chose `bitvec` (already a direct dep via umol-shared) over a hand-rolled
  `Vec<u64>` (don't-reinvent rule) and over `fixedbitset` (`Vec`-backed → not
  const-constructible, so it can't replace the isotope `const BitArray` → would mean
  a second bit crate, no gain). No crate offers non-allocating dense word-popcount,
  so the similarity ops are hand-written on the raw slice regardless.
- Decisions resolved: **runtime width** (not const-generic); **no `Encoder` enum
  yet** — a direct `fold` method, because `MinHash`'s output is a different type and
  would not fit a single-return encoder; introduce `Encoder`/`minhash` when MinHash
  lands. Deferred: bit-info provenance map; **count fingerprints** (`CountFp` +
  counted feature source) — first needed by slice 6's `Difference`.
- Deps: none new.

### Slice 4 — Morgan (b): frozen RDKit replica

Circular refinement + bond-set dedup (the ECFP machinery) but with RDKit's
*invariant* and *hash*, pinned to a named RDKit revision and frozen — the point is
to nail down RDKit's moving target so a result is reproducible. RDKit's invariant
differs from R&H (we verified: total degree not heavy degree, no heavy valence,
delta-mass), so this is genuinely a distinct featurizer, not ECFP with a flag.

- Builds: the RDKit connectivity invariant (atomic number, total degree, total
  #H, charge, delta-mass, ring); a Morgan hash path for `circular_refine` (a new
  `CircularRefinementAlgorithm` variant or a hash abstraction); offline RDKit
  fixtures (stored expected outputs, no FFI at test time — the subiso-conformance
  pattern).
- Decision (**key**): **bit-exact RDKit replica** (replicate boost `hash_combine`
  + fixtures — heavier, but reproduces existing RDKit Morgan/ECFP6 fingerprints and
  freezes the unstable reference) **vs RDKit-invariant + our `xxh3`** (lighter,
  frozen, but does not reproduce RDKit's bits). The doc-stated intent ("as RDKit
  produces it, replicated") is the bit-exact path; confirm given the cost.
- Deps: `circular_refine` (built) + a hash path; RDKit fixtures.

### Slice 5 — Substructure (f): two directions on a shared enumeration primitive

The literature/impl review (open question 2) found that RDKit's *sound*
substructure-screening FP is `PatternFingerprint` — a **fixed library of ~13 small
templates** (single atom; 1 bond; 2–3-atom paths; branched 4–5-atom; rings 3–6;
ring-junction motifs) matched by subgraph isomorphism, keyed by **element +
bond-type only** (the embedding-monotone invariants; degree/ring-count are
deliberately excluded), dropping any feature that touches a query/uncertain atom so
`query ⊆ target` has no false negatives. The path-based "RDKit FP"
(`Fingerprints.cpp`) is a *similarity* FP and is **not** subset-sound (it bakes
within-subgraph degree and atom counts into the key).

Two directions are built, sharing one core enumeration primitive:

**Core primitive — `umol-graph-core`, bounded enumeration.** Independently
justified (topological torsion = length-3-bond paths, BRIDGIT, future RDKit path
FP), so it is built here regardless of the FPs. New algorithm module; output = edge
sets (`Vec<Vec<EdgeId>>`) — a feature is an edge set, since a ring and its spanning
path differ only in edges. Selector enums proxy to named algorithms
(algorithm-transparency convention):
- `enumerate_paths(max_length, PathEnumerationAlgorithm::Dfs)` — simple paths ≤
  `max_length` bonds; bounded DFS, deduped by endpoint orientation.
- `enumerate_connected_subgraphs(max_size, SubgraphEnumerationAlgorithm::Esu)` —
  connected (branched) subgraphs ≤ `max_size` bonds; **Wernicke ESU** (*Efficient
  detection of network motifs*, 2006 — exclusive-neighborhood expansion,
  duplicate-free).
- Plus **edge-induced subgraph extraction** `Graph::edge_induced_subgraph(&[EdgeId])
  -> Embedding` (today only node-induced exists, which cannot isolate a path from a
  chorded ring).

**Direction B — RDKit `PatternFingerprint` replica (`umol-graph::fingerprint`).**
Bit-exact, pinned to RDKit 2026.03.3 (same revision as Morgan), offline fixtures —
for head-to-head benchmark comparison.
- RDKit's fixed template set expressed as `MoleculeAst` patterns (no SMARTS query
  layer — the constraints already exist): element `*` = `ElementAst::Undetermined`
  (the match wildcard); `[R]`/ring-size = `AtomConstraint::RingMembership(RingScope::
  {All, Size(s)})`; `@` (ring bond) = `BondConstraint::RingMembership`; `~` =
  `BondAst` with undetermined order.
- Matched via existing `substructure_matches`; key = RDKit gboost hash over per-atom
  element + per-bond type (aromatic collapsed) in template order; query-touching
  matches dropped; fold `% width`.
- Validate via fixtures (Morgan discipline): ring perception must align with
  RDKit's; match multiplicity (`uniquify=false` occurrence/count bits) must match.
  A leaner topological matcher is deferred unless efficiency demands it.

**B status (this round).** `PatternFingerprinter` → `BitFp(2048)` in
`fingerprint/pattern.rs`; keying via a domain-private copy of the frozen RDKit boost
hash (its consolidation with core's identical copy is an **open placement decision**
— the hash is RDKit-specific, not graph-generic, so it must not enter the core
public API); output built by collecting raw 32-bit feature+count ids and reusing
`FeatureSet::fold(2048)`.
Bit-exact vs RDKit 2026.03.3 for **6 fixtures**: ethanol, benzene, pyridine, furan,
chlorobenzene, cyclohexane — validating the keying, the chained count bits, match
multiplicity, query-index `mv` order, ring-perception alignment (incl. the pyridine
fix), and cycle embedding. The small-ring templates are `*`-atom cycles (no `[R]`):
a cycle already forces ring membership, so matches are identical and there is no
ring-fact dependency. **Blocked:** the four junction templates (RDKit indices 9–12,
fused rings) use `@` ring **bonds** in non-cyclic patterns, which cannot be reduced
to cycle topology. Naphthalene confirms the gap precisely — umol's 117 bits are a
strict subset of RDKit's 126 (9 missing, 0 wrong), i.e. exactly the junction
features. Ring-bond / `[R]` matching needs ring-membership facts evaluated during
substructure matching, but `derive_constraints` (umol-ast) does not synthesize them
for atoms/bonds — a gap vs the spec, which classifies `#R` as a derived predicate
evaluated against the target graph. Closing it (synthesize ring facts in
`derive_constraints`) is a load-bearing umol-ast change that unblocks the junction
templates and ring/`@`/`[R]` queries generally.

**Direction A — experimental unhashed screen (`umol-graph::fingerprint`).**
Generative, collision-free.
- Enumerate paths + connected subgraphs (the core primitive) up to size L.
- Each subgraph → **exact canonical key**: edge-induced sub-`Graph` → nauty
  `canonical_labeling(node_color)` (color = embedding-monotone chemical invariants)
  → serialize `node_count │ colors (canonical order) │ edges (u,v,bond_label)
  sorted` into a `Vec<u8>`. Collision-free by construction (it is the canonical form
  serialized, not a hash).
- Keys collected into `FeatureSet`; **`Id` bound relaxed `Copy+Ord → Clone+Ord`**
  (the key is a variable-length `Vec<u8>`; safe for existing `FeatureSet<u64>`).
- Screen with `is_subset`; soundness from embedding-monotone labels only (element,
  bond type — not degree/ring membership), so every query feature is derivable in
  any superstructure target.

Decisions resolved: B via `MoleculeAst` (no query layer); paths = DFS, subgraphs =
ESU; exact `Vec<u8>` key with `Clone + Ord`. No finite automaton exists for the
general substructure check (subgraph iso is NP-complete; Courcelle tree automata
give linear time only for fixed patterns on bounded-treewidth graphs, with
impractical automaton size). Deferred: canonical-SMILES key (vs byte key) until
DRFP needs it; a leaner topological matcher for B; unlabeled-graph features.
- Deps: new core enumeration + edge-induced extraction + canonicalization.

### Slice 6 — Reaction FPs from ECFP/Morgan (difference / role-tagged)

Apply a *molecular* featurizer (ECFP/Morgan) to each reaction component, then
combine across roles — the "reaction difference fingerprint" (RDKit) and a
role-tagged variant. Mode (a) of the reaction model with a molecular featurizer
(vs DRFP's reaction-native shingles).

- Builds: reaction consumption over `ReactionRuleAst` (lhs/rhs, doc-127 interim);
  per-component featurization (run the chosen molecular featurizer on each reactant
  and product); a `ReactionCombinator` enum — `Difference` (signed: Σ product
  features − Σ reactant features) and `RoleTagged`. The molecular featurizer is a
  parameter.
- `RoleTagged` = tag each feature by its side (reactant / product) and union into
  one set; equivalently, when folded it is the **concatenation**
  `[reactant_fp ‖ product_fp]` (the role tag is just which block a feature lands
  in). Unlike `Difference` it does **not** cancel: a feature unchanged by the
  reaction appears in both blocks, so the unchanged scaffold is retained (larger but
  information-preserving). "Role" = the side (all reactants unioned, all products
  unioned), not the individual component. (The name `RoleTagged` is interim — `Role`
  is not used elsewhere in the codebase yet.) -> DisjointUnion
- Note: `Difference` is computed from per-component fingerprints summed by role —
  **no atom map needed** (correcting the earlier "uses the atom map" wording;
  atom-mapped / CGR handling is mode (b), separate and later). `RoleTagged` needs
  only the role too.
- Decisions: **`Difference` is signed/count-valued**, so it needs the count
  representation (the deferred `CountedFeatureSet`) — this slice either pulls counts
  in or ships `RoleTagged` (binary) first and defers `Difference`.
- Deps: ECFP (built) / Morgan (slice 4) + the reaction type + counts (for
  `Difference`). Establishes the reaction-combinator framework slice 7 reuses.

### Slice 7 — DRFP (d) - **Deferred**

Reaction-native (Probst 2022): radius-bounded subgraph *shingles* from both sides
as canonical-SMILES keys, combined by **symmetric difference** — no atom map.

- Builds: shingle extraction (reuse slice-5 enumeration → canonical SMILES); the
  `SymmetricDifference` reaction combinator (added to slice-6's `ReactionCombinator`)
  → `FeatureSet` (then optionally folded).
- Decisions: **canonical-SMILES dependency** (umol-io writer + a canonical atom
  order — confirm it exists / what it needs); reaction representation stays the
  doc-127 interim.
- Deps: slice-5 enumeration + canonical SMILES + slice-6 reaction framework.

### Slice 8 — BRIDGIT (c) **Deferred**

Circular fingerprints over substructure-matched reactive sites (Hadadi 2019):
fingerprint the environment around the reaction center.

- Builds: reactive-site queries via existing substructure matching, centering
  circular refinement on matched sites; reaction-center environment aggregation.
- Decisions: reactive-site query expressivity (may need a SMARTS-level query layer
  beyond `MoleculeAst`-as-pattern); centering/aggregation scheme.
- Deps: substructure matching (built) + `circular_refine` (built) + reactions;
  possibly the SMARTS query layer.

### Cross-cutting (threaded through the slices, not standalone)

- **Count fingerprints** (`CountedFeatureSet` + `CountFp`, and `cosine`) — for
  purpose-2 count vectors and count-Tanimoto. The first concrete need is slice 6's
  `Difference` combinator (signed counts); it ripples back into the featurizers
  (they must track multiplicity), so its scope is decided there.
- **SMARTS query layer** — MACCS keys, pharmacophore features, richer reactive
  sites. Deferred; gates the keyed/feature fingerprints and richer BRIDGIT.

### Decisions wanted before their slice

1. **Morgan (slice 4)**: resolved — bit-exact RDKit replica, built.
2. **Substructure (slice 5)**: resolved — two directions (B = RDKit `PatternFingerprint`
   replica via `MoleculeAst` patterns; A = experimental unhashed generative screen) on
   a shared core path+subgraph enumeration primitive (DFS + ESU). Exact `Vec<u8>` key,
   `FeatureSet` `Id: Clone + Ord`.
3. **Reaction FPs (slice 6)**: ship `RoleTagged` (binary) first, or pull in the
   count representation now so `Difference` lands with it.

(Slice-3 `BitFp` width is resolved: runtime.)

## Scope guardrails

- Deliverables a–f now; the omissions (AP/TT, path, MACCS, MHFP/MAP4,
  pharmacophore, CGR, difference-FP) are accommodated by the architecture but not
  implemented in this pass.
- Enum-dispatched algorithms with explicit names and citations; frozen hashing;
  offline reference fixtures for anything claiming external compatibility (Morgan).

## References

- Rogers & Hahn 2010, *Extended-connectivity fingerprints* —
  `materials/fingerprints/Rogers and Hahn 2010 - Extended-connectivity fingerprints.pdf`.
- Hadadi et al. 2019 (BRIDGIT) —
  `materials/fingerprints/hadadi-et-al-2019-…pdf`, `pnas.1818877116.*`.
- Probst et al. 2022 (DRFP) — `materials/fingerprints/Probst-DRFP-2022.pdf`.
- RDKit fingerprint generators — `materials/codes/rdkit` (Code/GraphMol/Fingerprints).
- molintern (Python reference; DRFP/RDKit wrapping, reaction handling, CGR) —
  `/Users/dr/Dayhoff/molintern/molintern` (`featurizer.py`, `cgr_graph.py`).
- Orsi & Reymond 2024, fingerprint comparison; scikit-fingerprints (Gugler/… 2024)
  as an API-shape reference; MAP4 (Capecchi, Probst & Reymond 2020); MHFP (Probst &
  Reymond 2018).
