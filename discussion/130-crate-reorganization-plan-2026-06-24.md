# 130 — Crate reorganization implementation plan

Implements the decisions in doc 129. Three independent workstreams:

- **A** — frozen-hash relocation + generic `circular_refine` (§iii).
- **B** — split `umol-shared` into `umol-chem` + `umol-utils` (§iv).
- **C** — new `umol-geometric-graph` bridge crate (§ii).

Decision **i** (no umbrella crate) is a no-op. The `umol-io` coordinate→stereo
perception stays put (moving it would cycle, per 129 §ii).

## Independence and ordering

The three workstreams touch disjoint code and can land in any order; none unblocks
another. Suggested order by containment: **A → B → C** (A is localized and exactly
verifiable against the RDKit fixtures; B is a mechanical mass-rename; C is the only one
with net-new code). Each lands as one atomic change; a transient red tree between them
is acceptable.

## A. Frozen-hash relocation (§iii)

Verified prerequisites: `graph-core` has no in-crate non-test caller of the refinement;
`refine` is already generic; `ECFP_SEED` already lives in `umol-graph`; `xxhash_rust` is
confined to `graph-core/src/algorithms/refine.rs`.

### graph-core changes (`src/algorithms/refine.rs`, `src/lib.rs`, `Cargo.toml`)

1. Add the recipe trait:
   ```rust
   pub trait CircularRefinementHash {
       fn seed_hash(&self, components: &[u32]) -> u64;
       fn combine(&self, round: u32, current: u64, neighbors: &[(u32, u64)]) -> u64;
   }
   ```
2. Make circular refinement generic, mirroring `RefinementAlgorithm<H>`:
   `enum CircularRefinementAlgorithm<H> { Ec { radius: u32, scheme: H } }`;
   `circular_refine<H: CircularRefinementHash>` / `circular_refine_ec<H>` call only
   `scheme.seed_hash(..)` and `scheme.combine(..)`. (`derive(Copy, PartialEq, Eq)` stays
   valid — the recipe types in B are `Copy`/`Eq`.)
3. Delete from `graph-core`: `EcScheme` + its impl, `gboost_combine`, `gboost_hash`,
   `RefinementXxh3Scheme`, `RefinementWidth` (+ `RefinementWidth64`/`128`),
   `RefinementAggregation`, `ALBATROSS_SEED`, `BULLFINCH_SEED`. Move the xxh3-specific
   tests (`refine.rs` ~590/603/657 and the `default_scheme()` test helper ~518) out with
   them. Keep the `CountingScheme` test impl — it keeps the generic `refine`/
   `circular_refine` covered in core via the trait seam.
4. `Cargo.toml`: drop `xxhash-rust`.
5. `lib.rs` re-export becomes:
   ```rust
   pub use algorithms::refine::{
       CircularRefinementAlgorithm, CircularRefinementHash, Refinement, RefinementAlgorithm,
       RefinementHash, RefinementRounds,
   };
   ```
   (removed: `EcScheme`, `RefinementAggregation`, `RefinementWidth`, `RefinementWidth64`,
   `RefinementWidth128`, `RefinementXxh3Scheme`; added: `CircularRefinementHash`.)

### umol-graph changes

1. New module `umol-graph/src/hash.rs` (`pub mod hash;` in `lib.rs`) — the single
   transparent home for frozen hashing. It holds:
   - `gboost_combine` / `gboost_hash` (one copy);
   - the seeds `ALBATROSS_SEED`, `BULLFINCH_SEED`, and `ECFP_SEED` (relocated from
     `fingerprint/ecfp.rs`);
   - `RefinementXxh3Scheme`, `RefinementWidth` (+ `64`/`128`), `RefinementAggregation`
     (moved verbatim; `xxhash-rust` is already a `umol-graph` dependency);
   - recipe types implementing `CircularRefinementHash`:
     `struct Morgan;` (gboost bodies) and `struct RogersHahn { seed: u64 }` (xxh3 bodies),
     lifted from the old `EcScheme` arms.
2. Update callers:
   | File | Change |
   | --- | --- |
   | `fingerprint/wl.rs` | import `RefinementXxh3Scheme`/`RefinementWidth64` from `crate::hash`; `Refinement`/`RefinementAlgorithm`/`RefinementRounds` still from `umol_graph_core` |
   | `fingerprint/ecfp.rs` | `EcScheme::RogersHahn { seed }` → `crate::hash::RogersHahn { seed }`; `ECFP_SEED` now in `crate::hash` |
   | `fingerprint/morgan.rs` | `EcScheme::Morgan` → `crate::hash::Morgan` |
   | `fingerprint/pattern.rs` | delete its private `gboost_combine`; use `crate::hash::gboost_combine` |

### Tests (A)
xxh3-scheme tests + the bit-exact `test_gboost_hash` move into `umol-graph/src/hash.rs`.
Morgan/ECFP/WL/pattern fixtures in `umol-graph` must remain bit-exact. `graph-core`
retains algorithm coverage via `CountingScheme`.

## B. Split `umol-shared` → `umol-chem` + `umol-utils` (§iv)

Largest churn, fully mechanical. Per the 129 split: `umol-chem` gets `element` (incl.
`MAX_ATOMIC_NUMBER`, the `e!` macro), `isotope`, `isotope_data`, `occupation`, `spin`,
`configuration`, `units`, and the error **enums** (`ElementError`, `IsotopeError`,
`OccupationError`, `SpinStateError`, `DataError`); `umol-utils` gets the `UmolError`
**trait** and `solution`/`Solution`.

### New crates
- `umol-utils`: `src/error.rs` (the `UmolError` trait only), `src/solution.rs`; `lib.rs`
  re-exports. Deps: only what those two modules use today.
- `umol-chem`: the listed modules verbatim; `src/error.rs` holds the five enums, which
  `impl umol_utils::UmolError` ⇒ **`umol-chem` depends on `umol-utils`** (the only edge).
  Deps: carve `umol-shared`'s current deps (bitvec, phf, regex, serde, strum, thiserror)
  by module use. Move `umol-shared/benches/elements.rs` → `umol-chem/benches`.
- Delete `umol-shared`.

### Import rewrite (mechanical; every site enumerated in the 129 recon)
All occurrences are on separate `use` lines — no brace-group mixes the two halves — so
this is per-path substitution:
- `umol_shared::{element,isotope,occupation,spin,configuration,units}` → `umol_chem::…`
- `umol_shared::e!` / `umol_shared::e` → `umol_chem::e`
- `umol_shared::solution::` → `umol_utils::solution::`
- `umol_shared::error::UmolError` → `umol_utils::error::UmolError`
  (no code imports the error enums by path, so `error::` always means the trait here.)

### Dependent `Cargo.toml` edits
| Crate | `umol-shared` → |
| --- | --- |
| `umol-geometric` | `umol-chem` |
| `umol-params` | `umol-chem` |
| `umol-graph` | `umol-chem` + `umol-utils` |
| `umol-io` | `umol-chem` + `umol-utils` |
| `umol-ast` | `umol-chem` + `umol-utils` |

Files with both halves (split the two `use` lines across the new crates):
`graph/ops/{aromaticity,invariant,valence/atom_typing,valence/counts}.rs`,
`io/table_ir/raise.rs`, `io/ctfile/error.rs`, `io/smiles/parser.rs` (+ a few test mods).

## C. Bridge crate `umol-geometric-graph` (§ii)

### Create the crate
`umol-geometric-graph`, deps: `umol-graph`, `umol-geometric`, `umol-ast`, `umol-chem`.
(Depends on `umol-graph` so it can resolve/process the perceived molecule downstream;
this is also what makes absorbing the io stereo perception a cycle — hence that stays in
`io`.)

### Move bond perception (mechanical)
Move `umol-geometric/src/bond_perception.rs` → `umol-geometric-graph/src/`; drop
`pub mod bond_perception;` from `umol-geometric/src/lib.rs`. Rewrite its imports
`use crate::{algorithms::optimization::…, molecule::Molecule}` →
`use umol_geometric::{algorithms::optimization::…, molecule::Molecule}`. No visibility
changes needed — `algorithms::optimization` and `molecule::Molecule` are already `pub`.
Its `#[cfg(test)]` module moves with it (it is the only current caller).

### Net-new code: geom → `MoleculeAst`
Add the conversion that makes perception consumable downstream: from a
`umol_geometric::molecule::Molecule` (elements + coordinates) run `perceive_bonds`, then
build a `MoleculeAst` (elements as atoms, the `Vec<(usize,usize,u8)>` as bonds). This is
the one piece that is written, not moved. Resolution/ops on the result are the caller's
job (via the `umol-graph` dependency).

### Out of scope
The graph→3D embedding (doc 071 distance-geometry port) is a separate follow-on; this
plan only relocates perception and adds the geom→AST direction.

## Workspace `Cargo.toml`
`members`: remove `umol-shared`; add `umol-chem`, `umol-utils`, `umol-geometric-graph`.

## Verification
- A: `cargo test -p umol-graph-core -p umol-graph` (esp. the RDKit bit-exact Morgan/
  ECFP/pattern fixtures and WL tests); `cargo build` the `umol-graph` benches.
- B: `cargo build --workspace && cargo test --workspace` (compilation is the proof for
  the path rewrite).
- C: `cargo test -p umol-geometric-graph` — relocated bond-perception tests pass in the
  new home; a test for the geom→`MoleculeAst` conversion.
- Whole workspace green + `cargo clippy --workspace --all-targets`.

## Risks / non-mechanical points
1. **C net-new conversion** — geom→`MoleculeAst` is the only authored code; everything
   else in A–C is a move or a path rewrite.
2. **A generic enum** — `CircularRefinementAlgorithm<H>` gaining a type parameter; the
   `Copy`/`Eq` derives must still hold (they do for `Morgan`/`RogersHahn`).
3. **B error-module split** — the `UmolError` trait and the error enums separate across
   the two crates, creating the single `umol-chem → umol-utils` edge.
