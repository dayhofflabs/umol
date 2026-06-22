# 123 — AST allocation/clone survey and prioritized worklist (2026-06-21)

## Method

Static pass over `umol-ast/src/ast`: every `.clone()` and every owned-returning accessor,
classified **read** vs **build** per the doc 122 principle (read paths must not clone;
`meet`/`join`/`canonicalize`/builder/`transact`/`edit`/`derive_constraints`/`from_*` own
legitimately). Cross-referenced with the callgrind hot spots (doc 121).

**Runtime weights are pending a refresh.** The `cg.out` numbers (`ValueAst::matches` 6.86%,
`AtomConstraints::matches` 4.85%, `ValueAst::clone` 2.24%) predate the empty-pattern
short-circuit, the `AromaticValenceAst`/`MulticenterValenceAst` cheap matches, and the
`Cow` host-borrow — so they overstate the predicate cost. Re-run `prof 300 10` under
callgrind and re-rank. The static priority below is primary until then.

## Already addressed

- Empty-pattern short-circuit in `AtomConstraints`/`BondConstraints::matches`.
- Allocation-free `matches`: `ValueAst`, `ElementAst`, `IsotopeMassAst`, `AromaticValenceAst`,
  `MulticenterValenceAst` (and `SpinStateAst`/`AtomAst`/`BondAst` via the derive's field-wise
  `matches`).
- `host_match_targets` borrows via `Cow` in the no-derive path (no per-atom host clone for
  element/bond patterns).

## Prioritized findings (read paths)

**P1 — Constraint accessors return owned (`v.clone()`), the read-stored anti-pattern.**
`constraint/atom.rs:505–605` — 14 accessors (`valence`, `total_valence`, `degree`,
`total_degree`, `ring_degree`, `ring_valence`, `total_hydrogens`, `donated_pairs`,
`accepted_pairs`, `ring_count`, `aromatic_valence`, `multicenter_valence`,
`tetrahedral_stereo`); `constraint/bond.rs:149,157` (`ring_count`, `cis_trans_stereo`). Each
clones the stored value. Called by the collection `matches` per candidate for *constrained*
patterns (the short-circuit only covers empty patterns). **Fix: doc 122 Option A** —
`-> Option<&ValueAst>` / `Option<&…Ast>`; `matches` compares by reference; `meet` clones at
the point it builds. This is the agreed design; the top read-path item.

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
- **P1 (was top) — constraint accessor references (doc 122 Option A).** Confirmed *not*
  visible here because element/bond patterns short-circuit past the accessors; it bites
  *constrained* SMARTS queries. Still the right structural fix, just not what this corpus
  exercises.
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
- **P1 — constraint accessor references (doc 122 Option A)** — for *constrained* patterns
  only (this corpus short-circuits past it).
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
