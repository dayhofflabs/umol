# Molecule DSL parsing — profiling and iterative optimization

Date: 2026-04-24
Status: Completed

## Motivation

Parsing benchmarks (`umol-ast/benches/parsing.rs`) show the DSL layer at
18–63 MiB/s, 2–5× slower than raw EDN parsing (~110 MiB/s). Size inspection
showed `AtomConstraints` at 440 B — 65 % of every `AtomAst` — and indole
(9 atoms) writes ~4 KB of constraint-slot zero-fill per parse even though no
atoms declare constraints. That made `AtomConstraints` a plausible target
for footprint reduction.

## 0. Initial guesses

The six candidate optimizations listed at the start of this investigation,
ranked by expected impact. These were derived from type-size inspection and
data-flow reading; none had profile data behind them.

1. **Make `AtomConstraints` empty-by-default / lazy.** Box it
   (`Option<Box<AtomConstraints>>`) or switch to an inline small-vec /
   `SmallMap<AtomConstraintKind, AtomConstraint>` that takes 0 bytes when
   empty. Gains fall through every `AtomAst::default()`, `from_atom_dsl`,
   and atom clone in the pipeline. For inputs without per-atom constraints
   (most molecules), this eliminates ~65 % of the `AtomAst` bulk and 11
   branches in every default construction.

2. **Drop the DSL intermediate for atoms/bonds during streaming.** The
   stream path currently does `AtomDsl` → `AtomAst` (via `FromAst`). For
   integer/keyword-indexed `:atoms` / `:bonds`, the DSL value is
   structurally trivial (just the atom string + id); parsing directly into
   `AtomAst` without an intermediate `AtomDsl` skips one 680-B construction
   per atom. Tree path can't share this, but tree already pays the `Edn`
   allocation so stream is the priority.

3. **`BondAst` at 152 bytes** has a similar story. Check whether
   `BondConstraints` / any `Vec` field dominates and apply the same
   empty-by-default treatment. 11 bonds in indole × 152 B = ~1.7 KB per
   parse.

4. **Metadata bimap on hot paths.** `bimap 0.6.3` double-hashes on every
   insert and lookup. For ref resolution, a single `IndexMap<String, u32>`
   (keyword → index) plus a reverse `Vec<String>` built only on demand is
   strictly cheaper and keeps insertion linear. Indole has 9 atom keyword
   inserts + 10 bond inserts + 10 bond endpoint lookups + 9
   aromatic-atoms-list lookups — double hashing all of that is visible at
   20 MiB/s.

5. **String allocations in `ElementAst::Bind`, `ValueAst::Expr(Var)`, bond
   types.** Every `AtomDsl` → `AtomAst` lift that carries a bind-name or
   var-name does an owned `String`. Test molecules don't use these, but
   pattern queries will, and `&str` → `Cow<'a, str>` on AST would avoid
   allocation when the input is a leaf string (e.g., ref ids during
   resolution).

6. **`Arc<Vec<AtomAst>>` construction order.** The streaming path builds a
   `Vec<AtomAst>` then wraps in `Arc`. If the parser can pre-reserve
   capacity (from the first `:atoms` vector header size), you avoid
   realloc as atoms accumulate. Check that `read_vec` uses the header
   count for `Vec::with_capacity`.

## 1. Baseline measurements (main @ f3f112c4)

### Type sizes

```
AtomAst                = 680 bytes
  AtomConstraints        = 440 bytes   ([Option<AtomConstraint>; 11])
  AtomConstraint         = 40 bytes
  Option<AtomConstraint> = 40 bytes
  ValueAst               = 32 bytes
  Expr                   = 32 bytes
  ElementAst             = 48 bytes
  IsotopeAst             = 32 bytes
BondAst                = 152 bytes
MoleculeAst            = 80 bytes
```

### Parsing benchmarks — umol-edn (raw `read_string`)

| Input            | Size   | Time    | Throughput |
| ---------------- | ------ | ------- | ---------- |
| `molecule_small` | 41 B   | 362 ns  | ~108 MiB/s |
| `molecule_large` | ~600 B | 5.01 µs | ~114 MiB/s |

### Parsing benchmarks — umol-ast `molecule_dsl`

| Input              | Size  | Tree     | Stream   | Tree thrpt | Stream thrpt |
| ------------------ | ----- | -------- | -------- | ---------- | ------------ |
| `small`            | 39 B  | 1.62 µs  | 1.26 µs  | 22 MiB/s   | 28 MiB/s     |
| `benzene`          | 145 B | 5.39 µs  | 4.05 µs  | 27 MiB/s   | 36 MiB/s     |
| `indole`           | 314 B | 17.73 µs | 15.31 µs | 18 MiB/s   | 20 MiB/s     |
| `with_constraints` | 253 B | 6.19 µs  | 4.01 µs  | 41 MiB/s   | 63 MiB/s     |

Observations:
- Streaming wins 1.2–1.5× vs tree across molecule inputs — saves the
  intermediate `Edn` walk but still pays DSL lift.
- Indole is slowest per byte (keyword ids + bimap lookups).
- Raw-EDN ceiling is ~110 MiB/s; DSL conversion cost is 2–5× raw parse.

## 2. Current assumptions

Revised model after the Section 0 #1 experiment (results in Section 8).

**Parser throughput is not dominated by AtomAst construction cost.**

Indole's per-atom budget is ~1.70 µs (15.31 µs / 9 atoms). The
`AtomConstraints` change captured the full expected ~40 ns/atom memset
savings — ~140 ns total across 9 atoms, observed as a ~1.8 % throughput
delta. The other ~98 % of the per-atom budget is spent elsewhere and has
not been measured.

**Hypothesis-driven size reduction has diminishing returns for parser
throughput.** Each subsequent size win will yield comparable (1–3 %) gains
because the arithmetic is bounded by memory-bandwidth physics at these
working-set sizes. Other phases (hash ops, allocations, parser state) are
the real candidates for the next iteration.

**Footprint wins still matter — for downstream workloads.** Cloning
`AtomAst`, matching large molecules, and reaction networks with millions of
atoms are memory-bandwidth-bound and scale with `sizeof(AtomAst)`. The 2×
`AtomAst` shrink from #1 pays off there even though parsing didn't show it.

### Status of the Section 0 guesses

| # | Target                               | Status     | Expected parser-throughput impact |
| - | ------------------------------------ | ---------- | --------------------------------- |
| 1 | `AtomConstraints` footprint          | Done       | ~2 % (measured, Section 8)        |
| 2 | Drop `AtomDsl` intermediate (stream) | Open       | Unknown; worth measuring          |
| 3 | `BondAst` footprint                  | Open       | Small (same class as #1)          |
| 4 | `Metadata` bimap replacement         | Open       | Potentially significant           |
| 5 | `String` → `Cow<'a, str>`            | Open       | Zero for current benches          |
| 6 | Pre-reserve Vec capacity             | Unverified | Small if not already done         |

### Working hypothesis for the next measurement pass

Profile indole stream parse with `samply` (see discussion with user about
profiling workflow). Expected dominant frames, in rough likelihood order:
`bimap` hash operations, `String` allocation, winnow atom-string parsing,
EDN byte scanning. The profile — not more guesses — determines the next
target.

### Policy

No further hypothesis-driven refactors in this area without profile data
behind them. Each candidate change should have a phase-isolation benchmark
that measures improvement at the phase it targets, independent of overall
throughput movement — so we can tell whether the change worked at its
target phase even if end-to-end throughput is bounded by something else.

## 3. Decision

Replace `AtomConstraints::slots: [Option<AtomConstraint>; 11]` with
`SmallVec<[AtomConstraint; 2]>`, kept kind-sorted with at most one entry per
`AtomConstraintKind`.

### Why Option A (SmallVec, inline N=2)

- Resolution sets at least `Valence`; sometimes also `AromaticValence`.
  Other kinds (`Degree`, `Connectivity`, `RingSize`, …) are rare.
- `N=2` covers the post-resolution common case with zero heap; parse-time
  atoms have 0 constraints, also inline.
- Canonical kind order maintained via binary-search insert, so `PartialEq`
  and `Hash` stay positional-on-sorted-contents — behaviour-equivalent to
  the array shape for every existing test.

### Rejected alternatives

- `Vec<AtomConstraint>`: equally simple but forces a heap allocation for the
  single-constraint post-resolution case (every atom).
- `u16 bitset + Vec`: faster `contains`, but same heap-allocation problem
  and more internal state to maintain.
- `Box<[Option<AtomConstraint>; 11]>` lazy: 0 heap when empty, but 440 B
  allocation whenever any constraint is set — wasteful for the 1–2-entry
  case.

## 4. Expected outcomes

| | Before | Target |
| ----------------- | ------ | ------ |
| `AtomConstraints` | 440 B  | ~88 B  |
| `AtomAst`         | 680 B  | ~328 B |
| Default zero-fill | 440 B  | 8 B    |

Parsing throughput hypothesis (proved wrong in Section 8):
- `indole/stream`: 20 MiB/s → 30–40 MiB/s
- `benzene/stream`: 36 MiB/s → 40–50 MiB/s
- `small/stream`: ~28 MiB/s (allocation is not the bottleneck at this size)

## 5. Scope

- `umol-ast/Cargo.toml`: add `smallvec = "1.15.1"` (already a workspace dep via
  umol-graph).
- `umol-ast/src/ast/constraint/atom.rs`: replace struct + rewrite methods.
  Invariant preserved internally: `entries` sorted ascending by `kind() as u8`.
- All `AtomConstraints` call sites unchanged (opaque API).

## 6. Out of scope for this change

- `BondAst` (152 B) — `BondConstraints` is already `Vec`-backed.
- `Metadata` bimap replacement.
- Direct `AtomDsl`-skipping lift in the streaming parser.
- `ValueAst` / `Expr` boxing experiments.

Decide on these after re-measuring.

## 7. Plan

1. Add `smallvec` dep to `umol-ast/Cargo.toml`.
2. Rewrite `AtomConstraints` in `umol-ast/src/ast/constraint/atom.rs`.
3. Run `cargo test -p umol-ast` — all existing tests must pass unchanged.
4. Re-measure type sizes (temporary test, then delete).
5. Re-run `cargo bench -p umol-ast --bench parsing` full sweep; append after
   numbers to this doc.
6. Mark status Completed or open follow-up discussion for the next target.

## 8. Results

Implemented in this pass. All 1280 umol-ast tests pass unchanged.

### Type sizes (after)

| Type              | Before | After | Δ    |
| ----------------- | ------ | ----- | ---- |
| `AtomConstraints` | 440 B  | 96 B  | −78% |
| `AtomAst`         | 680 B  | 336 B | −51% |
| `AtomConstraint`  | 40 B   | 40 B  | —    |
| `BondAst`         | 152 B  | 152 B | —    |

Default zero-fill per atom: 440 B → 8 B (just the SmallVec len tag).

### Parsing benchmarks (after — `cargo bench -p umol-ast --bench parsing`)

| Input                     | Before    | After     | Δ       |
| ------------------------- | --------- | --------- | ------- |
| `small/tree`              | 22 MiB/s  | 23 MiB/s  | +3–10%  |
| `small/stream`            | 28 MiB/s  | 29 MiB/s  | +3%     |
| `benzene/tree`            | 27 MiB/s  | 28 MiB/s  | <1%     |
| `benzene/stream`          | 36 MiB/s  | 37 MiB/s  | +2–3%   |
| `indole/tree`             | 18 MiB/s  | 18 MiB/s  | +1%     |
| `indole/stream`           | 20 MiB/s  | 21 MiB/s  | +1.8%   |
| `with_constraints/tree`   | 41 MiB/s  | 41 MiB/s  | +1.8%   |
| `with_constraints/stream` | 63 MiB/s  | 64 MiB/s  | +2.8%   |

### Interpretation

Footprint target landed fully; throughput moved 1–3 %, not 30–40 MiB/s as
hypothesized in Section 4. The captured ~140 ns/parse matches the physical
ceiling for the 432-B-per-atom memset savings we targeted — the prediction
mis-estimated the fraction of per-atom work represented by that memset.

The updated model of the system is recorded in Section 2.

## 9. Iteration 2 — eliminate `atom_only_metadata` cloning

### 9.1 Profile data (after iteration 1, `samply` on indole stream)

Top frames under `MoleculeDsl::from_edn_str`:

| Share | Frame                                                                |
| ----- | -------------------------------------------------------------------- |
| 73 %  | `MoleculeInput::into_ast`                                            |
| 24 %  | `read_molecule_input` (byte parsing)                                 |
| 3 %   | `MoleculeDsl::drop`                                                  |

Inside the 73 % `into_ast` slice:

| Share | Frame                                                                |
| ----- | -------------------------------------------------------------------- |
| 19 %  | `Map::into_iter().collect::<Result<Vec<AtomIdx>, ParseError>>()` (aromatic-atoms resolution) |
| 18 %  | `IndexMap<AtomIdx, String>::clone()`                                 |
| 16 %  | `Metadata::drop()`                                                   |
| 3.5 % | `MoleculeAst::new()`                                                 |

Inside the 24 % byte-parsing slice: 11 % atoms, 9 % bonds, 3 % aromatic systems.

### 9.2 Root cause

`MoleculeInput::into_ast` was calling `atom_only_metadata(&atom_ids)` inside
every ref-resolution site — 29 times per indole parse (20 bond endpoints,
9 aromatic atoms). Each call cloned the full `atom_ids:
IndexMap<AtomIdx, String>` into a throwaway `Metadata`, then dropped it.
Ref resolution further did a linear `.iter().find(...)` over that cloned
map for the `Id(name)` case — while the reverse map
`atom_id_to_idx: IndexMap<String, AtomIdx>` was **already built** at
`molecule.rs:765–771` and sitting in scope, unused by the resolver.

### 9.3 Fix

- Added `$name::resolve(self, count, id_to_idx: &IndexMap<String, $idx>) -> Result<$idx, ParseError>` to the `define_ref!` macro
  (`constraint.rs:141`). O(1) `id_to_idx.get(&name)` instead of linear scan.
  Existing `into_ast(self, count, &Metadata)` left intact for constraint
  resolution (which still uses `ResolveContext` / `Metadata`).
- Updated 8 call sites in `MoleculeInput::into_ast` (bond endpoints, dative
  donor/acceptor, aromatic atoms, multicenter atoms, noncovalent endpoints)
  to call `resolve(atom_count, &atom_id_to_idx)` — no `Metadata` constructed.
- Removed the now-dead `atom_only_metadata` helper.

### 9.4 Results

All 1280 umol-ast tests pass unchanged. Benchmark deltas (vs iteration 1):

| Input                     | Iter 1    | Iter 2    | Δ       |
| ------------------------- | --------- | --------- | ------- |
| `small/tree`              | 23 MiB/s  | 23 MiB/s  | +2.5%   |
| `small/stream`            | 29 MiB/s  | 31 MiB/s  | +4.6%   |
| `benzene/tree`            | 28 MiB/s  | 30 MiB/s  | +8.7%   |
| `benzene/stream`          | 37 MiB/s  | 42 MiB/s  | +10.4%  |
| **`indole/tree`**         | 18 MiB/s  | 30 MiB/s  | **+42%** |
| **`indole/stream`**       | 21 MiB/s  | 41 MiB/s  | **+97%** |
| `with_constraints/tree`   | 41 MiB/s  | 44 MiB/s  | +6.5%   |
| `with_constraints/stream` | 64 MiB/s  | 71 MiB/s  | +9.3%   |

Indole stream: 15.17 µs → 7.71 µs. Hypothesis (2× throughput from
eliminating the 53 % share) verified essentially exactly.

Scales with ref count as expected: `small` has no keyword refs (minimal
win), `benzene` has integer-indexed bonds/aromatic (small win from endpoint
validation only), `indole` has all-keyword refs (maximal win).

### 9.5 Next target

Re-profile `indole/stream` after this change; the new hot frames set the
next target. Candidates still open from Section 0: #2 (drop `AtomDsl`
intermediate in streaming), #4 (bimap replacement for aliases, likely
minor impact now since alias resolution isn't in the entity-loop hot path),
#6 (verify `read_vec` pre-reserves). Expect the share of byte-parsing
(`read_vec` sites) to rise proportionally now that `into_ast` is cheaper.

## 10. Iteration 3 — `parse_atom` fast path for bare element symbols

### 10.1 Profile data (after iteration 2, `samply` on indole stream)

Total time: 7.71 µs. Two near-equal halves:

| Share | Frame                                                                   |
| ----- | ----------------------------------------------------------------------- |
| 47 %  | `read_molecule_input` (byte parsing)                                    |
| 47 %  | `MoleculeInput::into_ast` (resolution)                                  |
| 3 %   | `MoleculeDsl::drop`                                                     |

Inside `read_molecule_input`:

| Share | Frame                                                   |
| ----- | ------------------------------------------------------- |
| 21 %  | `read_vec AtomEntryInput`                               |
|       | — 11 % `parse_atom` (winnow run on every atom-string)   |
|       | — 3 % `AtomDsl::new`                                    |
| 20 %  | `read_vec BondEntryInput`                               |
| 6 %   | `read_vec AromaticSystemEntryInput`                     |

Inside `into_ast`, frames were all 4–8 % — no single dominant hotspot.

### 10.2 Root cause

Every atom-string, including the overwhelming majority that are a bare
element symbol ("C", "N", "O", "Cl", …), ran through the full winnow
parser at `dsl/atom.rs:105`. The winnow state machine (whitespace,
alternation, predicate repeat) is overkill for a 1- or 2-byte input.

### 10.3 Fix

Added `parse_bare_element` pre-check at the top of `parse_atom`:

```rust
pub fn parse_atom(input: &str) -> Result<AtomDsl, ParseError> {
    if let Some(dsl) = parse_bare_element(input) {
        return Ok(dsl);
    }
    atom.parse(input).map_err(|e| e.into_inner())
}
```

`parse_bare_element` accepts the input only if it is exactly 1 ASCII
upper byte or exactly an ASCII upper + ASCII lower byte pair and the
pair is a valid `Element` symbol per `Element::from_symbol_bytes`.
Everything else — whitespace, wildcards, predicates, sets, binds —
falls through unchanged to the winnow parser.

### 10.4 Results

All 1280 umol-ast tests pass.

| Input                     | Iter 2    | Iter 3    | Δ       |
| ------------------------- | --------- | --------- | ------- |
| `small/tree`              | 23 MiB/s  | 27 MiB/s  | +15%    |
| `small/stream`            | 31 MiB/s  | 37 MiB/s  | +18%    |
| `benzene/tree`            | 30 MiB/s  | 35 MiB/s  | +14%    |
| `benzene/stream`          | 42 MiB/s  | 51 MiB/s  | +19%    |
| `indole/tree`             | 30 MiB/s  | 34 MiB/s  | +11%    |
| **`indole/stream`**       | 41 MiB/s  | 47 MiB/s  | **+13%** |
| `with_constraints/tree`   | 44 MiB/s  | 46 MiB/s  | +4%     |
| `with_constraints/stream` | 71 MiB/s  | 76 MiB/s  | +7%     |

Indole stream: 7.71 µs → 6.71 µs. One microsecond saved across 9 atoms ≈
111 ns/atom of winnow overhead eliminated. Captures the expected
~11 % share plus a little of the 3 % `AtomDsl::new` sibling frame.

### 10.5 Cumulative

Since the baseline at the top of Section 1:

| Input                     | Baseline  | Iter 3    | Δ        |
| ------------------------- | --------- | --------- | -------- |
| `small/stream`            | 28 MiB/s  | 37 MiB/s  | +32%     |
| `benzene/stream`          | 36 MiB/s  | 51 MiB/s  | +42%     |
| **`indole/stream`**       | 20 MiB/s  | 47 MiB/s  | **+135%** |
| `with_constraints/stream` | 63 MiB/s  | 76 MiB/s  | +21%     |

### 10.6 Opaque Metadata refactor (landed alongside this iteration)

Metadata fields made private; added accessor methods (`atom_id`,
`bond_id`, `atom_alias_for`, `iter_atom_aliases`, etc.) and a crate-
private `MetadataBuilder` used by the parser and test fixtures.
Alias-path types changed from `AtomDsl` (336 B stack value) to
`Box<AtomDsl>` (8 B pointer) — localized to alias storage
(`Metadata.atom_aliases`, `MetadataBuilder.atom_aliases`,
`MoleculeInput.atom_aliases`, local `alias_table` in `into_ast`).
`MetadataBuilder::add_atom_alias(&mut self, name, Box<AtomDsl>) ->
Result<(), ParseError>` replaces the earlier `try_add_atom_alias` that
imitated `BiMap::insert_no_overwrite`'s `Result<(), (L, R)>` shape.
`ResolveContext` inherits opacity through `metadata: &Metadata` —
counts kept as `pub` (primitives).

No measurable effect on these benches (indole has no aliases); the
refactor unlocks future storage changes without changing call sites.

### 10.7 Next target

Re-profile indole stream after this change. Candidates still open:
#2 (drop `AtomDsl` intermediate during streaming) — likely smaller
payoff now that `parse_bare_element` skips the intermediate for the
common case. #6 (verify `read_vec` pre-reserves from header count).
Expect `read_molecule_input` share to drop below 47 % and `into_ast`
to rise relatively.

## 11. Iteration 4 — Metadata id-storage experiment (negative result)

### 11.1 Motivation

Section 0 #3 / a revised form of target (b) from earlier: the
`Metadata.atom_ids` IndexMap and siblings pay a hash-table cost on
every lookup and insert. Hypothesis: for common inputs where many
maps are empty, short-circuiting with `is_empty()` would save the
hash-empty cost; and for inputs with few ids, replacing IndexMap with
`Vec<(Idx, String)>` would avoid hashing entirely.

### 11.2 New bench infrastructure

Added `benches/rendering.rs` (separate target from `benches/parsing.rs`)
and a shared `benches/fixtures.rs` carrying the existing molecule
fixtures plus three new large-molecule fixtures parameterized by id
density:

- `large_no_ids`: 100 integer-indexed carbons, 99 bonds.
- `large_all_ids`: 100 carbons, each with `[:a{i} "C"]` keyword id.
- `large_partial_ids`: 100 atoms, one id every ten (10 ids total).

Render bench times `dsl.to_edn()` — the id-lookup-heavy path used on
AST → EDN round-trips.

### 11.3 Baselines (current IndexMap, after iteration 3)

Parse (stream):

| Input              | Time    | Thrpt     |
| ------------------ | ------- | --------- |
| small              | 0.94 µs | 37 MiB/s  |
| benzene            | 2.84 µs | 51 MiB/s  |
| indole             | 6.70 µs | 46 MiB/s  |
| with_constraints   | 3.30 µs | 76 MiB/s  |
| large_no_ids       | 21 µs   | 64 MiB/s  |
| large_all_ids      | 35 µs   | 56 MiB/s  |
| large_partial_ids  | 22 µs   | 62 MiB/s  |

Render (`to_edn`):

| Input              | Time     | Thrpt      |
| ------------------ | -------- | ---------- |
| small              | 0.33 µs  | 108 MiB/s  |
| benzene            | 1.33 µs  | 109 MiB/s  |
| indole             | 3.29 µs  | 95 MiB/s   |
| with_constraints   | 2.25 µs  | 113 MiB/s  |
| large_no_ids       | 13.70 µs | 110 MiB/s  |
| large_all_ids      | 24.40 µs | 89 MiB/s   |
| large_partial_ids  | 17.20 µs | 92 MiB/s   |

### 11.4 Experiment: `is_empty` short-circuit in the 6 id accessors

```rust
pub fn atom_id(&self, idx: AtomIdx) -> Option<&str> {
    if self.atom_ids.is_empty() {
        return None;
    }
    self.atom_ids.get(&idx).map(String::as_str)
}
```

(applied to all 6 accessors)

Render deltas:

| Input              | Baseline   | After      | Δ       |
| ------------------ | ---------- | ---------- | ------- |
| small              | 108 MiB/s  | 108 MiB/s  | +0.6 %  |
| benzene            | 109 MiB/s  | 111 MiB/s  | +1.9 %  |
| indole             | 95 MiB/s   | 97 MiB/s   | +1.3 %  |
| with_constraints   | 113 MiB/s  | 113 MiB/s  | −0.3 %  |
| large_no_ids       | 110 MiB/s  | 113 MiB/s  | +2.1 %  |
| large_all_ids      | 89 MiB/s   | 89 MiB/s   | −0.3 %  |
| **large_partial_ids** | **92 MiB/s** | **87 MiB/s** | **−5.4 %** |

Parse deltas within noise (±1 %).

### 11.5 Interpretation

The change was **reverted**. Two findings invalidate the premise:

1. **IndexMap's empty-map fast path is already fast.** The
   hashbrown-backed `get` bails on `capacity == 0` before hashing;
   measured savings on `large_no_ids` (all maps empty) were
   ~1.76 ns/call — far below the ~25 ns/call I predicted.

2. **The added branch costs on non-empty lookups.** On
   `large_partial_ids`, `atom_ids` is populated (10 entries) but other
   id maps (`bond_ids`, etc.) are empty. The `atom_id` accessor runs
   on a non-empty map for 298 lookups; each pays the extra
   `is_empty` check before hashing. Net regression of ~5 %.

The "no-ids optimization" hypothesis does not hold as a distinct
improvement separable from storage choice — a good hash-map already
has this built in.

### 11.6 Re-evaluation of `Vec<(Idx, String)>` storage

Using the now-measured per-call cost of ~36 ns on a populated IndexMap
(derived from `large_no_ids` → `large_all_ids` delta / 298 calls):

| Input (render)        | IndexMap now | Vec prediction | Δ         |
| --------------------- | ------------ | -------------- | --------- |
| `large_no_ids` (K=0)  | 110 MiB/s    | ~110 MiB/s     | neutral   |
| `large_partial_ids` (K=10) | 92 MiB/s | ~120 MiB/s     | +30 %     |
| `large_all_ids` (K=100) | 89 MiB/s   | ~30 MiB/s      | **−66 %** |

Vec wins decisively for K ≤ 15 but creates a cliff past K ≈ 50
(linear scan on every atom_id call during render). Since the real
workload distribution is unknown and id-heavy patterns
(pharmacophores, anchor-tagged protein fragments) are plausible, a
storage choice that degrades badly past a threshold is not safe.

### 11.7 Decision

**Keep IndexMap.** Drop target (b) entirely.

- The IndexMap is already near the ceiling on empty-map access.
- Vec optimizes a ~10 % win on mid-density inputs at the cost of a
  potential ~3× regression on id-heavy inputs.
- Current DSL throughput (46 MiB/s indole stream, 89–113 MiB/s render
  across the suite) is within 20 % of the raw-EDN ceiling for the
  worst case, at the ceiling for the best.

The new benchmark infrastructure (`benches/rendering.rs` and the
`large_*` fixtures) stays — useful for regression tracking on any
future `Metadata`-shape changes.

### 11.8 Status of remaining Section 0 items

| # | Target                              | Status |
| - | ----------------------------------- | ------ |
| 1 | AtomConstraints footprint           | Done (iter 1) |
| 2 | Drop AtomDsl intermediate stream    | Open; expected smaller impact post-iter-3 |
| 3 | BondAst footprint                   | Open; low priority |
| 4 | Metadata bimap alias replacement    | Done structurally (opaque API + Box<AtomDsl>); no perf work |
| 5 | String → Cow on AST                 | Open; relevant for pattern queries, not benches |
| 6 | Pre-reserve Vec capacity            | Unverified |
| b | Metadata id-map storage change      | Dropped (iter 4) |

## 12. Iteration 5 — `repr(transparent)` DSL newtypes, zero-copy render

### 12.1 Audit finding

Six render sites in `molecule.rs::render_*` clone entity AST data only to wrap
in a throwaway DSL newtype for a `&self`-only method call (`Display` /
`ToEdn`):

- `render_atom_entry` line 620: `AtomDsl(atom.clone())` — 336 B per atom.
- `render_bonds` line 644: `BondDsl(view.data.clone()).to_edn()` — 152 B per bond.
- `render_dative`, `render_aromatic`, `render_multicenter`,
  `render_noncovalent`: similar clones per entity.

For `large_all_ids` render (100 atoms + 99 bonds): ~50 KB of memcpy per
render, all of it thrown away.

### 12.2 Fix

Add `#[repr(transparent)]` to the six entity DSL newtypes
(`AtomDsl`, `BondDsl`, `DativeBondDsl`, `AromaticSystemDsl`,
`MulticenterBondDsl`, `NoncovalentBondDsl`) and a `from_ref` method:

```rust
#[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct AtomDsl(pub AtomAst);

impl AtomDsl {
    pub fn from_ref(ast: &AtomAst) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(ast as *const AtomAst as *const Self) }
    }
}
```

Render sites become `XDsl::from_ref(&view.data).to_edn()` instead of
`XDsl(view.data.clone()).to_edn()` — a zero-cost reference cast.

### 12.3 Results

All 1296 umol-ast tests pass. Parse benches move within ±2 % (noise).
Render deltas:

| Input                    | Baseline   | After      | Δ         |
| ------------------------ | ---------- | ---------- | --------- |
| small                    | 108 MiB/s  | 136 MiB/s  | +25 %     |
| benzene                  | 109 MiB/s  | 132 MiB/s  | +19 %     |
| indole                   | 95 MiB/s   | 110 MiB/s  | +13 %     |
| with_constraints         | 113 MiB/s  | 118 MiB/s  | +4 %      |
| **large_no_ids**         | 110 MiB/s  | 144 MiB/s  | **+28 %** |
| large_all_ids            | 89 MiB/s   | 104 MiB/s  | +17 %     |
| **large_partial_ids**    | 92 MiB/s   | 118 MiB/s  | **+35 %** |

### 12.4 Interpretation

`large_no_ids` and `large_partial_ids` saw the biggest wins as predicted —
they have the largest per-entity-count × per-entity-size clone budget. All
fixtures are now above 100 MiB/s render, most above 110 MiB/s — at or above
the raw-EDN parse ceiling (~110 MiB/s). `with_constraints` got the smallest
gain (+4 %) because its render time is dominated by constraint-tree
serialization, not per-entity dispatch.

The change is localized: 6 struct attributes + 6 tiny `impl` blocks + 6
render-site rewrites. No change to public API; existing `XDsl(ast)` tuple
construction still compiles.

### 12.5 Next

Remaining Section 0 candidates still open: #2 (drop AtomDsl intermediate
in streaming parse), #3 (BondAst footprint), #5 (Cow on AST — relevant
only for pattern queries), #6 (pre-reserve Vec capacity). Parse path
still has room: `id.clone()` sites in the dedup-plus-insert pattern
(Section 11 audit finding C) are worth considering if parse becomes the
next target.

## 13. Iteration 6 — audit point B: `IntoAst` per-entity clones

### 13.1 Motivation

Audit finding B (Section 11 follow-up): `IntoAst for MoleculeDsl` at
`molecule.rs:562` cloned each entity while walking `atoms_mut()`,
`bonds_mut()`, etc., even though the slot was owned exclusively and
about to be overwritten. Pattern was:

```rust
for atom in ast.atoms_mut() {
    *atom = AtomDsl(atom.clone()).into_ast(&cfg.atom)?;
}
```

### 13.2 New bench target

Added `benches/conversion.rs` covering the `FromAst` and `IntoAst`
paths on `MoleculeDsl` — neither is exercised by the parse/render
benches (those bypass the trait via `MoleculeInput::into_ast`).

### 13.3 Baselines (current code, before fix)

`FromAst` (AST → DSL):

| Input               | Time    | Thrpt       |
| ------------------- | ------- | ----------- |
| small               | 0.40 µs | 87 MiB/s    |
| benzene             | 0.96 µs | 151 MiB/s   |
| indole              | 1.29 µs | 244 MiB/s   |
| with_constraints    | 0.68 µs | 372 MiB/s   |
| large_no_ids        | 9.9 µs  | 153 MiB/s   |
| large_all_ids       | 9.9 µs  | 220 MiB/s   |
| large_partial_ids   | 9.9 µs  | 160 MiB/s   |

`IntoAst` (DSL → AST):

| Input               | Time     | Thrpt       |
| ------------------- | -------- | ----------- |
| small               | 0.61 µs  | 58 MiB/s    |
| benzene             | 1.66 µs  | 87 MiB/s    |
| indole              | 2.56 µs  | 122 MiB/s   |
| with_constraints    | 0.97 µs  | 262 MiB/s   |
| large_no_ids        | 24.2 µs  | 62 MiB/s    |
| large_all_ids       | 25.3 µs  | 86 MiB/s    |
| large_partial_ids   | 24.6 µs  | 64 MiB/s    |

`IntoAst` is ~2.5× slower than `FromAst` on the large inputs — a clear
signal for the per-entity clone hypothesis.

### 13.4 Fix

Replace `atom.clone()` with `std::mem::take(atom)` in the four entity
loops of `IntoAst for MoleculeDsl`. The slot is owned exclusively and
about to be overwritten, so moving out and leaving `Default::default()`
in place is sound. `FromAst for MoleculeDsl` does not need the same
treatment — its entity loop uses `AtomDsl::from_ast(&atom, cfg)` which
takes `&atom` and has no inner clone.

### 13.5 Results

All 1296 umol-ast tests pass.

| Input                    | Baseline  | After     | Δ     |
| ------------------------ | --------- | --------- | ----- |
| small                    | 58 MiB/s  | 57 MiB/s  | noise |
| benzene                  | 87 MiB/s  | 93 MiB/s  | +7 %  |
| indole                   | 122 MiB/s | 130 MiB/s | +5.7% |
| with_constraints         | 262 MiB/s | 272 MiB/s | +4.4% |
| large_no_ids             | 62 MiB/s  | 64 MiB/s  | +3.5% |
| large_all_ids            | 86 MiB/s  | 89 MiB/s  | +3.8% |
| large_partial_ids        | 64 MiB/s  | 68 MiB/s  | +5.6% |

### 13.6 Interpretation

The win is smaller than the ~20 % I predicted from byte-count napkin
math. Two factors probably explain the gap: (a) the derived `Clone` on
`AtomAst` is already efficient (post-iter-1 the slotmap is a 96 B
`SmallVec`, not the old 440 B fixed array), and (b) the compiler likely
RVOs part of the write-then-overwrite into a single store.

The 900 ns saved on `large_no_ids` works out to ~4.5 ns per entity —
below a cache-line write. The change is still worth keeping: idiomatic,
touches correctness-neutral code, and the bench is now a regression net
for future edits to this path.

### 13.7 Point B closed

`IntoAst` now uses `mem::take`. `FromAst` uses the already-optimal
`&atom`-taking lower. No further action on point B.
