# 123 — AST allocation/clone survey and prioritized worklist (2026-06-21)

Status: In Progress
Date: 2026-06-21
Relates: [237](237-release-0-10-scope-2026-10-15.md)

## Current worklist

- Adopt the measured removal of recursive matcher allocations: S0a and S0b below.
- Adopt the measured incidence-vector reservation in S0c, then run the shared gate.

Scope review, 2026-10-19: the predicate fast path and uncommon stereo-form
optimizations are not selected for 0.10. The host Vec<Cow> buffers and
MoleculeEmbedding are gone; returned-correspondence allocation is not part of
this work. Mutation and transactions are closed here. The reference-only
constraint-reading proposal was superseded by the semantic readings in 236.
The historical rankings and mutation-audit proposals below are not additional
implementation tasks.

## Method

Static pass over `umol-ast/src/ast`: every `.clone()` and every owned-returning accessor,
classified **read** vs **build** per the doc 122 principle (read paths must not clone;
`meet`/`join`/`canonicalize`/builder/`transact`/`edit`/`derive_constraints`/`from_*` own
legitimately). Cross-referenced with the callgrind hot spots (doc 121).

The initial `cg.out` numbers (`ValueAst::matches` 6.86%,
`AtomConstraints::matches` 4.85%, `ValueAst::clone` 2.24%) predate the empty-pattern
short-circuit, the `AromaticValenceAst`/`MulticenterValenceAst` cheap matches, and the
`Cow` host-borrow, so they overstate the predicate cost. The fresh profiles
later in this document supersede them.

## Already addressed

- Empty-pattern short-circuit in `AtomConstraints`/`BondConstraints::matches`.
- Allocation-free `matches`: `ValueAst`, `ElementAst`, `IsotopeMassAst`, `AromaticValenceAst`,
  `MulticenterValenceAst` (and `SpinStateAst`/`AtomAst`/`BondAst` via the derive's field-wise
  `matches`).
- `host_match_targets` borrows via `Cow` in the no-derive path (no per-atom host clone for
  element/bond patterns).
- Constraint accessors return `Option<&T>` rather than cloning stored values.

## Prioritized findings (read paths)

**P1 — Constraint accessors returned owned values. Completed.**
`constraint/atom.rs:505–605` — 14 accessors (`valence`, `total_valence`, `degree`,
`total_degree`, `ring_degree`, `ring_valence`, `total_hydrogens`, `donated_pairs`,
`accepted_pairs`, `ring_count`, `aromatic_valence`, `multicenter_valence`,
`tetrahedral_stereo`); `constraint/bond.rs:149,157` (`ring_count`, `cis_trans_stereo`). Each
cloned the stored value and was called by collection `matches` for constrained
patterns. Implemented using doc 122 Option A:
`-> Option<&ValueAst>` / `Option<&…Ast>`; `matches` compares by reference; `meet` clones at
the point it builds.

**P2 — Stereo value-type `matches` still meet-derived.** `TetrahedralStereoAst`,
`CisTransStereoAst` (the `stereo.rs:297` Lattice macro), `StereoConfigurationAst`
(`stereo.rs:185`). Reached via the `tetrahedral_stereo()`/`cis_trans_stereo()` accessors for
*stereo* patterns; the meet-derived default allocates. **Fix:** cheap, coset-aware `matches`
(deferred from the value-enum pass — needs care with coset frames). Cold path (stereo
patterns only), so below P1.

**P3 — `MoleculeEmbedding::from_match` per-match allocation.** `embedding.rs:81,96` build
`host_atoms`/`host_bonds` `Vec`s per accepted match; `restrict`/`restore` (`:187–222`)
allocate `HashSet`s. Scales with match count (e.g. `branched` ≈ 17k matches over the 9k
corpus). Build path, but on the per-match hot loop. **Fix:** assess buffer reuse / lazy
embedding; weight depends on match density (the fresh profile decides).

**P4 — `incidence_graph` construction.** `incidence.rs:301,340` (+ the Levi-graph build):
allocated per `substructure_matches` call on the `Incidence` strategy. **Fix:** assess; only
the `Incidence` strategy pays it.

## Not the anti-pattern (legitimate owned)

- **View compute-accessors** (`view/atom.rs:151–457`, 17 owned `-> ValueAst`): these
  *compute* derived values from topology (valence = Σ bond orders, etc.) — they construct,
  so owned is correct. Invoked only via `derive_constraints` (a build path, gated off for
  unconstrained patterns). Not a read-path clone.
- **Build/mutation clones** (deprioritized): `molecule/transact.rs` (83), `builder.rs` (18),
  `edit.rs` (12), `rewrite.rs` (8); `meet`/`join`/`canonicalize` in `value.rs`/`atom.rs`;
  `derive_constraints`; `from_*` constructors; `*/tests.rs`. All construct or mutate — owning
  is correct.

## Fresh profile (post-fix, HEAD 329486988, `prof 300 10`)

Callgrind, branched/Vf2Rdkit. Total **119.6M Ir**, down from the pre-fix 198M (−40%);
match cost ≈ 118M → ~40M (−66%). Validated: `AtomConstraints::matches` 8.4% → 2.8% (now
just the `is_empty` short-circuit per candidate, no iteration), `ValueAst::clone` 3.3% →
1.3%, `AromaticValenceAst::meet` gone from the match path.

Two structural facts from this run:

- **Parse dominates (~⅔).** `xxhash` is a fixed 3,582,598 Ir regardless of reps — one-time
  SMILES raise/resolve (~267k Ir/molecule); ditto `write_str`. As match shrank, parse's
  share grew. Use `prof 300 50` (or a collect-toggle) to isolate match. Parse/raise is a
  separate, amortizable target.
- **Match cost is now call-volume, not allocation.** `ValueAst::matches` (10.8M) dominates
  the match path — the ~5 direct atom + ~4 bond value-field comparisons per candidate, each
  trivial (`Undetermined → true`) but called field-count × candidates. `vf2rdkit_search` is
  ~2.7M. So the residual ~2–3× vs RDKit is field-wise generality (6+ fields/candidate vs
  element + bond order), not the old meet/clone waste.

## Re-ranked worklist

- **R1 — skip trivially-true pattern fields in field-wise `matches`.** New top match cost.
  For element-only patterns, charge/h/lone_pairs/spin are `Undetermined`; don't call
  `matches` on them. Derive could emit `pat.f.is_undetermined() || pat.f.matches(t.f)`, or
  the bigger version, a specialized ground/element fast path. Highest now.
- **Completed baseline — constraint accessor references (doc 122 Option A).**
  The change is not visible here because element/bond patterns short-circuit
  past the accessors; it benefits constrained queries.
- **R2 — per-call `host_match_targets` `Vec<Cow>` allocation** (one per `substructure_matches`
  call). Borrow the host directly in the no-derive closure to drop the Vec.
- **P2 (stereo value-type matches), P3 (embedding), P4 (incidence)** — unchanged priority;
  cold or strategy-specific.
- **Parse/raise** — now the largest absolute cost; its own investigation (resolver +
  aromaticity perception during SMILES raise), amortized in screening.

## Parse-clean profile (`prof 300 50`, HEAD 329486988)

Total 276.6M Ir; parse ≈ 80M (29%, down from ~67% at reps 10), match ≈ 196M (71%).
`ValueAst::clone` is **out of the top entirely** — the predicate is fully allocation-free;
fixes completely validated. No single dominant cost; match splits three ways:

- **Field-wise predicate ~53M** — `ValueAst::matches` 33.6M (12%, top) + `AtomAst`/
  `ElementAst`/`IsotopeMass` matches. Cheap per call, high volume (6+ fields × candidates).
- **Search + graph ~39M** — `vf2rdkit_search`/`_next` (18.6M) + `Graph::neighbors` (~16M).
  Partly inherent (any matcher walks adjacency); RDKit pays this too.
- **Allocations ~52M (parse+match)** — `_int_free` #2 (6.96%). With `clone` gone, the
  match-side allocation is **per-call buffers**: VF2 working vectors (graph-core), the host
  `Vec<Cow>`, the embedding. Reusable, not inherent.

## Final ranking (diminishing returns — the 4× is banked)

- **L1 — predicate fast path** (~27% of match). Ground/element matching skips trivially-true
  fields; must be matching-specific (target = satisfiable host), *not* the generic derive
  (`matches(Undetermined, ⊥)` is `false`, so `is_undetermined(pat) || …` breaks the law).
- **L2 — per-call buffer reuse / prepared-matcher** (~19% allocation bucket). VF2 allocates
  working state per call (the `_int_free` #2); a prepared-matcher / buffer-reuse API in
  graph-core helps *every* subiso algorithm. doc 104 flagged these benches as the driver for
  exactly this decision. The host `Vec<Cow>` (one per call) folds in here.
- **L3 — `Graph::neighbors` / search** (~20% of match). Mostly the algorithm; least
  recoverable (RDKit pays it).
- **Completed — constraint accessor references (doc 122 Option A).** This
  affects constrained patterns; the profiling corpus short-circuits past it.
- **Parse/raise (~80M one-time)** — largest absolute, amortized in screening; separate.

## In-place mutation review (build/mutate paths) — intent

The read/build split (doc 122) is too coarse for the build side. The honest **trichotomy**:

- **read** → references (doc 122), no clone;
- **transform-in-place** → `&mut self`, no allocation, when the caller owns the value and
  does not need the old one (`A + B → C` becomes `A.op(B)` mutating `A → A'`);
- **build-new** → allocate, *only* when old and new must coexist.

Build-path clones are **not** monolithically necessary — many are clone-then-modify that a
`&mut` would avoid. The in-place infra partly exists already: `Lattice::narrow_from(&mut
self, other)` is the in-place `meet`; `canonicalize(self)` is by-move (no input clone);
`canonical(&self)` is `Cow` (clones only on change). The work is auditing whether call sites
*use* them.

Prime target: **the resolver** (next critical review) — does its constraint narrowing
`narrow_from` in place, or `meet`-then-assign? Then `molecule/transact.rs` (83 clones —
distinguish rollback snapshots, which are necessary for transactionality, from gratuitous
clone-then-modify), `edit.rs`, `builder.rs`.

Distinct from doc 122 (read-path): different correctness surface (mutation of an otherwise
value-semantic `MoleculeAst`, kept to well-delimited contexts). Tracked here; pursued after
doc 122. **Interning is explicitly out of scope** — not a destination now; the APIs are
designed without regard for it.

## Accepted matcher experiment — 2026-10-19

The experiment compared main e3b2348225aadf8e73b2c8507fca883d9534f30f with
two internal changes in umol-graph-core: depth-based VF2 restoration and
allocation-free Vf2Rdkit candidate traversal. Both are accepted for adoption.
They preserve the matching algorithms, candidate order, predicate calls,
anchoring, and visitor termination behavior. Public signatures, ownership,
and failure behavior are unchanged.

VF2 replaces its per-branch restoration vectors with another pass over the
selected query and target nodes' adjacency, clearing only terminal markers
introduced at the current depth. This trades the saved lists for neighbor
reads; it does not scan the whole graph. Vf2Rdkit reads a mapped parent's
adjacency directly, or the target node range for a component root, instead of
collecting a candidate vector at every recursive step.

Measurements used Rust 1.99.1 on macOS 15.7.3 arm64, release builds, the same
dependency lockfile and benchmark sources, and separate binaries for timing
and allocation counting. Three paired passes alternated execution order.
Each Criterion case used 30 samples, 100 ms warm-up, and 400 ms measurement.
The seven controls were 3-atom paths in 8- and 64-atom carbon chains, an
absent-element query in the 64-atom chain, a path in a branched 4-atom host,
a six-member ring in itself and in two fused six-member rings, and one
3-atom multicenter overlay with contributions [2, 0, 1]. Ordinary fields used
literal or undetermined forms. These controls do not define workload frequencies.

Median reductions in complete molecule-matching time, including returned
correspondences and their destruction:

| Case | GraphAndOverlays / VF2 | GraphAndOverlays / Vf2Rdkit | Incidence / VF2 | Incidence / Vf2Rdkit |
| --- | ---: | ---: | ---: | ---: |
| 8-atom chain | 8.7% | 7.5% | 15.8% | 8.6% |
| 64-atom chain | 3.1% | 2.8% | 6.2% | 3.9% |
| Six-member ring | 10.7% | 14.7% | 23.2% | 16.9% |
| Fused rings | 10.9% | 13.9% | 21.8% | 18.7% |

The 64-atom chain's core matcher allocations fell from 504 to 6 for VF2 and
from 193 to 4 for Vf2Rdkit; per-call working storage remains. The only slower
full-matching case was VF2 with GraphAndOverlays on the multicenter control:
0.3%, about 0.01 microseconds, with paired changes of +0.04% to +0.40%.
The no-match VF2 Incidence case was near unchanged: paired changes ranged
from -1.84% to +1.39%. The experiment passed 111 focused unit tests, six
matching properties, exact comparisons with saved baseline correspondences,
formatting, and graph-core Clippy with warnings denied.

The separately measured incidence construction for the 64-atom chain and
its pattern took about 2.17 microseconds. This was roughly 1% of the successful
Incidence matching call and 73–79% of the corresponding no-match call.
That supports keeping construction separate from the accepted matcher changes;
it does not select a cache, prepared matcher, or other construction API.

## Accepted incidence-construction experiment — 2026-10-19

The subsequent experiment isolated one change to Molecule::incidence_graph:
count the selected participant occurrences and reserve the endpoint and incidence
vectors before filling them. Fixed-degree counts use arithmetic; variable-degree
overlays require an extra pass over rows to read list lengths. The existing
construction loops and Graph::new remain unchanged. This change is accepted for
adoption. ReactionSpan::incidence_graph is outside the selected scope.

Both binaries used main e3b2348225aadf8e73b2c8507fca883d9534f30f plus the
accepted matcher changes above. The compiler, platform, separate timing/allocation
binaries, and three alternating paired passes matched the preceding experiment.
The 51 cases covered complete matching with both algorithms and strategies,
separate pattern/host construction, and empty, mixed-overlay, and 32-multicenter-row
construction controls at all three incidence levels. Timing used 30 samples,
100 ms warm-up, and 400 ms measurement. The controls retain the preceding
experiment's literal/undetermined forms and do not define workload frequencies.

| Measured operation | Median paired time reduction |
| --- | ---: |
| 64-atom host construction | 16.8% |
| 8-atom host construction | 21.4% |
| Mixed-overlay Full construction | 23.4% |
| 32 multicenter rows, Constitution construction | 16.1% |
| Complete no-match call, Incidence / VF2 | 12.0% |
| Complete no-match call, Incidence / Vf2Rdkit | 11.1% |

The 64-atom host's eight allocations remain, while its nine reallocations become
zero. Cumulative requested bytes fall from 6,920 to 5,790; this is allocation
traffic, not retained or peak memory. Most successful enumeration timings changed
by less than 1%. A focused repeat of chain8, ring6, and fused10 used three paired
passes with one-second measurements. The initially slower fused10/VF2 Incidence
median changed from +0.88% to +0.03%; timing excursions also occurred in unchanged
GraphAndOverlays controls. The evidence supports construction and rejection gains,
without establishing a general successful-enumeration speedup.

Empty construction consistently costs about 1–5 ns more, or 3–10%, without saving
allocations. The added counting work is retained with that measured trade-off.
Exact comparisons preserved 51 complete incidence graphs and the seven fixtures'
complete correspondences across both algorithms and strategies. The experiment
passed 85 focused unit tests, two matching properties, formatting, and graph-IR
Clippy with warnings denied.

## Implementation plan

### S0 — Matching allocations

Adopt the tested implementations in the working branch. The experiments are
complete in separate worktrees; these subitems remain open until adoption
and the shared gate are complete. The changes stay inside existing implementations;
there are no new public symbols, DSL forms, Python methods, or guide contracts.

#### S0a — VF2 terminal restoration

- [ ] **Additive (green). [dep: none]**
- Module: umol-graph-core, algorithms/subgraph_isomorphism.rs, Vf2State::search.
- Remove restore_q and restore_t. After recursive return, revisit the affected
  adjacency and clear entries whose terminal marker equals the exiting depth,
  then decrement depth. Preserve earlier markers, including the seeded anchor,
  and restore state before propagating ControlFlow::Break.
- Carry over the experiment's property assertion in
  tests/property/subgraph_isomorphism.rs: anchored results equal the restriction
  of unanchored Ullmann results to the supplied anchor. Retain the existing
  algorithm cross-validation, planted matches, visitor/enumerator agreement,
  disconnected-query, branching, cycle, and early-termination coverage.

#### S0b — Vf2Rdkit candidate traversal

- [ ] **Additive (green). [dep: none]**
- Module: umol-graph-core, algorithms/subgraph_isomorphism.rs,
  vf2_rdkit::vf2rdkit_search.
- Replace the collected candidate vector with traversal of the borrowed
  adjacency slice or target node range. Preserve candidate order and
  repetitions, degree and predicate checks, edge consistency, and state
  restoration on normal and early return. Both anchored and unanchored entry
  points continue to use this function.
- Run the existing matching unit cases and cross-algorithm properties,
  including the anchored restriction assertion from S0a when both are present.
  The two implementations are independent; apply them in the listed order.

#### S0c — Incidence vector reservation

- [ ] **Additive (green). [dep: none]**
- Module: umol-graph-ir, ir/incidence.rs, Molecule::incidence_graph.
- Compute the exact selected incidence count before construction: two per localized
  bond; under Constitution or Full, donors plus acceptor per dative bond, participant
  lengths for aromatic and multicenter entities, and two per noncovalent bond;
  under Full, site plus ligand count per stereo atom or stereo bond.
- Initialize both edges and incidences with that capacity. Count variable-degree
  rows through their existing exact-size iterators' lengths, without traversing
  their participants. Retain the existing filling loops and Graph::new call.
- Preserve entity-block order, edge and role alignment, selected-level behavior,
  and parallel participant occurrences, including shared virtual-ligand anchors.
  Keep the change confined to the molecule producer tested in the experiment.
- Run the existing exact incidence cases for all levels and roles, including
  parallel occurrences, and the substructure unit and cross-strategy/cross-algorithm
  property suites. No new fixture files or property generators are needed.
- This implementation is independent of S0a/S0b; apply it after them so the combined
  result matches the measured baseline plus candidate.

#### S0 shared gate — after S0a, S0b, and S0c

Run the focused graph-core matching, graph-IR incidence, and graph-IR matching
tests and their explicitly enabled property suites:

```sh
cargo test -p umol-graph-core --features proptest subgraph_isomorphism
cargo test -p umol-graph-ir --features proptest --lib ir::incidence
cargo test -p umol-graph-ir --features proptest substructure
```

Then run formatting, workspace tests, and workspace Clippy with warnings denied
once for the combined changes, using the repository's Python 3.13 environment
for workspace commands that include umol-py. Review the complete diff against
the accepted matcher and incidence changes and the anchored-matching property
assertion.

The recorded experiments settle the performance decisions for the tested code.
Repeat affected paired measurements only if adoption changes the implementation
or exposes a new performance concern.

Execution order: S0a, S0b, S0c, then the shared gate. The three implementations
are independent; the final gate depends on all three. No further optimization
experiment or stage remains selected in this document.
