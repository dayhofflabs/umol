# 106 — Extract `umol-io` and `umol-geometric-core` from `umol-graph`

Status: Done (steps 1–3) · 2026-06-07 · open item: relocate io benches/bins to `umol-io`

## Problem

MOL/V2000 encodes tetrahedral stereochemistry canonically with **wedges** (an up/down
bond, narrow end at the center) read against the 2D depiction. Resolving a wedge to a
configuration is **geometric**: a wedge fixes the out-of-plane direction of one bond, but
the handedness is the in-plane angular order of the *other* substituents, which lives only
in the coordinates. Up/down plus the neighbor-list order is two-to-one onto configuration —
two drawings with the same atom order and the same wedge but mirror-image in-plane layouts
are opposite enantiomers (see doc 104, the wedge RESOLVED note).

Two constraints collide:

- Wedges are common and canonical for MOL; they cannot be left unsupported.
- `umol-graph` must not perform geometric perception. Coordinates are passed through, not
  interpreted; depending on `umol-geometric` is out of the question; scattered geometry
  helpers in `umol-graph` are a rot vector.

Doc 104's Phase B papered over this: it labeled the wedge path "a symbolic descriptor"
while its pseudocode read `mol.positions` and computed a 2D winding determinant
(`wedge_winding`). That is a coordinate pass mislabeled — it contradicts the same doc's
"positions pass through / no geometric perception" rule, and the implementation copied it.
This doc corrects that by **relocating** the geometry, not deleting the capability.

## Observation

MOL/SDF is not a graph format — it carries 2D/3D coordinates and wedge depictions. The
3D-in-`umol-graph` tension was tolerable only while coordinates were pure pass-through.
Wedge resolution breaks that. The format-I/O layer (parsers + the `TableIR` boundary type +
the TableIR→AST raise, including whatever geometry a format needs) does not belong in the
graph/chemistry crate.

## Decision

Extract two crates.

- **`umol-geometric-core`** — geometric primitives only; leaf crate (dep: `umol-shared`),
  distinct from the heavy `umol-geometric` (symmetry, tensors). Modules:
  - `orientation` — `signed_volume(a, b, c, d) -> f64`: the signed volume of four points
    (3×3 determinant of the edge vectors); its **sign** is the orientation / handedness.
    Callers take the sign; the primitive stays a plain scalar.
  - `plane` — `complementary_direction`: the open in-plane direction (the renamed
    `leftover_inplane` — the negated mean of unit vectors to the in-plane neighbors).
  - `Point3D`.
- **`umol-io`** — SMILES/SMARTS/MOL/SDF parsing, the `TableIR` boundary type, and `raise`
  (TableIR → *unresolved* `MoleculeAst`, including wedge → `#T` via `umol-geometric-core`).
  Deps: `umol-ast`, `umol-perm`, `umol-geometric-core`, `umol-shared`.

`umol-graph` keeps `ops` (resolver / valence / aromaticity / model) and the graph
algorithms (Morgan, substructure, …), plus the resolve-on-parse convenience wrappers.

## Dependency direction: `umol-graph` → `umol-io`

io is the input boundary; graph consumes its output — the same role `umol-ast` already
plays for graph. The reverse (`umol-io` → `umol-graph`) is wrong:

- io needs nothing from graph. SMILES aromatic flags are lexical, MOL reads explicit bond
  orders; aromaticity *perception* and valence resolution are `ops`, applied at resolve
  time, not parse time.
- `raise` references `TableIR`, so it must live where `TableIR` lives (io). If `TableIR`
  moved to io while `raise` stayed in graph, raise would force graph → io, and io → graph
  would then cycle.

So `raise` moves to io with `TableIR`. It is mechanical: `raise.rs` imports only
`umol_ast`, `umol_perm`, `position`, and `table_ir` — no `crate::ops`.

## The one seam: resolve-on-parse

The only place io-output meets graph-`ops` today is the resolving convenience
(`parse_smiles_bytes_with` → `ops::Resolver`). Split it:

- **io** keeps the pure parses: `parse_*_to_table_ir`, `parse_*_to_ast` (syntax → TableIR →
  unresolved AST).
- **graph** keeps the resolving wrappers (`parse_smiles` = `io::parse_smiles_to_ast` +
  `Resolver`), since it owns `ops`.

The result is a DAG:

```
umol-graph → umol-io → { umol-ast, umol-perm, umol-geometric-core, umol-shared }
umol-graph → { umol-ast, umol-io, umol-shared }
umol-geometric-core → umol-shared
```

## Wedges after the split

`raise` in io may read `positions` (io is allowed geometry, via `umol-geometric-core`) and
produce `#T` for the wedge case; `umol-graph` never sees a coordinate. Wedges are supported,
the geometry is contained in a named crate, there is no `umol-geometric` dependency, and no
helpers are loose in graph. The `(None, _)` wedge arm and `tetrahedral_wedge_ordering`
survive in io; `wedge_winding` **splits** — its signed-volume determinant becomes
`orientation::signed_volume` in `umol-geometric-core`, and its StereoLigand-aware assembly
(lift the wedged point ±z, place the virtual ligand via `plane::complementary_direction`,
map the orientation sign → source index) stays io-side and calls the primitive. Moved and
refactored, not deleted.

## Error model (doc 065 conformance)

The split conforms the error types to doc 065's three tiers rather than carry the current
mixed model into `umol-io`. State today: `UmolError` lives in `umol-shared`; `ops` already
conforms (`ResolverError` / `ResolverContradiction` are tier-2 dispatch enums `#[from]`-
wrapping the `Valence` / `Aromaticity` / `Bonds` / `MulticenterBonds` sub-concerns). The one
non-conformant holdout is `SmilesError` / `CtfileError`, which `#[from]`-wrap both the io
`ParseError` (tier-1) and the ops `ResolverError` (tier-2) in a single enum — mixing concerns
across crates, which 065 forbids (no `#[from]` across tiers; cross-module = box). They are
**deleted**.

Target:

- **`umol-io`** — `smiles::ParseError` and `ctfile::ParseError` stay as tier-1 flat
  sub-concern enums and gain `impl UmolError`. Pure parse fns return `ParseError` directly. A
  tier-2 `IoError` dispatch is **deferred** — nothing currently spans formats.
- **`umol-graph`** — `ResolverError` and `ResolverContradiction` gain `impl UmolError`; a unit
  `ResolveUnderdetermined` error is added in `ops::resolver` (`impl UmolError`). The
  resolve-combining fns cross crates (io parse + ops resolve), so per tier-3 they return
  `Box<dyn UmolError>`: `Determined → Ok`, `Underdetermined → Err(Box::new(ResolveUnderdetermined))`,
  `Contradictory(c) → Err(Box::new(c))`, io `ParseError` boxed via `?`.
- The parse-and-raise fns that never resolve (`parse_*_to_ast`) are retyped from the deleted
  format errors to `ParseError` (their resolve variants were unreachable).

## Implementation plan

1. **`umol-geometric-core`** — *done 2026-06-07.* New leaf crate; `Point3D` + `all_zero`
   moved from `umol-graph/src/position.rs`; `orientation::signed_volume` and
   `plane::complementary_direction` added; `raise.rs::wedge_winding` refactored to call them
   (its `StereoLigand` assembly + `tetrahedral_wedge_ordering` stay raise-side); `position`
   imports repointed; `position.rs` deleted; `pub mod position` removed.

2. **Conform errors + decouple the resolve layer** — *done 2026-06-07*
   in-place in `umol-graph` (single crate, verifiable). Leaves `io`/`table_ir` ops-free.
   1. `impl UmolError` for `smiles::ParseError`, `ctfile::ParseError`, `ResolverError`,
      `ResolverContradiction`; add unit `ResolveUnderdetermined` in `ops::resolver`
      (`impl UmolError`).
   2. New `umol-graph/src/parse.rs`: move the eight resolve-combining fns (`parse_smiles`,
      `parse_mol`, each with `_bytes` / `_with` / `_bytes_with`) here, returning
      `Result<MoleculeAst, Box<dyn UmolError>>` with the outcome mapping above.
   3. Delete `SmilesError` / `CtfileError`; retype `parse_*_to_ast` to `ParseError`; the pure
      parse fns (`*_to_table_ir*`, `extended`, `sdf`, `reaction`) keep `ParseError`.
   4. Fix re-exports (`io/smiles.rs`, `io/ctfile.rs`, `lib.rs`) and callers.
   5. Verify: production `grep crate::ops src/io src/table_ir` → none; build + tests green.

3. **Create `umol-io`; move the ops-free modules.** — *done 2026-06-07.*
   1. New crate `umol-io` (deps: `umol-ast`, `umol-perm`, `umol-geometric-core`,
      `umol-shared`, plus io's external crates).
   2. Move `io/`, `table_ir/` (incl. `raise.rs`), `span.rs`, `diagnostics.rs` →
      `umol-io/src/`; wire `lib.rs`. Internal `crate::{io,table_ir,span,diagnostics}` paths
      stay valid; external refs (`umol_ast`, `umol_perm`, `umol_geometric_core`) unchanged.
   3. Pure-raise tests move with `raise.rs`; the three resolve-pipeline raise tests
      (`*_counts_resolve`, `*_resolver_*`) move to `umol-graph`.
   4. `umol-graph` depends on `umol-io`; `parse.rs` (and any graph code that used
      `crate::{table_ir,io}`) reference `umol_io::*`. (`ops` references none — verified.)
   5. Move the parse-only suites (`smiles_parsing`, `mol_parsing`, `sdf_parsing`,
      `smiles_property`) to `umol-io/tests/`; `resolution` stays in `umol-graph/tests/`.
   6. Build + test both crates; conformance + clippy green.

   *Results.* `umol-io/src/lib.rs` keeps the `io` module (paths are `umol_io::io::smiles::…`),
   so every internal `crate::io::…` path was untouched — only `umol-graph`'s `parse.rs`/bins/
   benches repointed `crate::io`→`umol_io::io`. The three resolve-pipeline tests landed in
   `umol-graph/src/parse.rs`'s test module (`counts_model` fixture + `METHANE_MOL`/
   `BENZENE_AROMATIC_MOL` duplicated there). `umol-io` dev-deps: `insta`, `serde`, `proptest`,
   `rstest`, `float-cmp`, `pretty_assertions`. `umol-graph` shed the io-only deps (`bstr`,
   `fast-float2`, `nom`*, `indexmap`, `smallvec`, `strum`, `bitflags`, `map-macro`, `num`,
   `itertools`, `bimap`, `index_vec`, `umol-perm`, `umol-geometric-core`) and the dead
   `proptest` feature/dev-dep. Tests: umol-io lib 3171 + mol 2250 + sdf 407 + smiles 10017 +
   property 3; umol-graph lib 169 + resolution 617 — no tests lost (3171+169 = the prior 3340).
   *`nom`/`umol-edn`/`clap` remain in `umol-graph` only for the io benches (`mol_parsing`,
   `smiles_parsing`) and the io bins (`test_*`, `classify_*`), which still live there — see
   open item.

   **Open item — io benches/bins still in `umol-graph`.** `benches/{mol_parsing,smiles_parsing}`
   and `src/bin/{test_mol_file,test_smiles,classify_mol_files,classify_sdf_files,
   classify_smiles_strings}` use only `umol_io` now; they keep `nom`/`clap` alive in
   `umol-graph` purely as io tooling (a lopsided dep). Relocating them to `umol-io` would let
   `umol-graph` drop `nom` (and `clap`, modulo the `molecule_dsl_parsing` bench's `umol-edn`).
   Deferred — not in this step's scope.

## Verification (2026-06-07) — clean, no cycles

Done read-only against the tree. The split is viable; the only crossings are reverse
(io→ops), and they sort into buckets that move to graph with the resolve layer.

1. **`ops` ↔ `table_ir` / `span` / `position`: none.** `ops/` references no `table_ir`,
   `span`, or `position` — it operates only on `MoleculeAst`. So moving io out creates **no**
   forward `umol-graph → umol-io` edge *from ops*; that edge exists only for the relocated
   resolve wrappers/errors. (Originally "enumerate TableIR uses in ops" — there are zero.)
2. **`umol-ast → umol-perm`: yes** (`umol-ast/Cargo.toml`). io→perm and graph→perm are
   cycle-free.
3. **`cx` (CXSMILES) → `ops`: none.** Parse-side only.
4. **Test suites split by what they exercise:** `smiles_parsing`, `mol_parsing`,
   `sdf_parsing`, `smiles_property` are parse-only (no `Resolver`/`resolve`; they call
   `parse_*_to_table_ir`) → `umol-io`. `resolution` resolves (8 resolver refs) → `umol-graph`.

Reverse couplings (io→ops) to relocate with the resolve layer (all expected):
- the resolve-on-parse wrappers in `smiles/parser.rs` + `ctfile/parser.rs` → graph.
- `SmilesError` / `CtfileError` (their `Resolve*` variants reference `ops::resolver`) →
  graph; the pure `ParseError` in each file stays in io. The two-enum split already exists.
- `raise.rs`'s `crate::ops` refs are **`#[cfg(test)]` only**; production raise is ops-free.
  Its three resolve-pipeline tests (`*_counts_resolve`, `*_resolver_*`) move to graph; the
  pure-raise tests stay in io.

Module ownership of today's `umol-graph`:

| module | → crate |
| --- | --- |
| `io`, `table_ir`, `span`, `diagnostics` | `umol-io` |
| `position` | `umol-geometric-core` |
| `ops` | stays `umol-graph` |

After the move `umol-graph` is reduced to `ops` + the relocated parse-resolve wrappers +
`Smiles`/`Ctfile` errors. This smallness is a consequence of the bottom-up implementation
strategy, not a mis-scoping — there is plenty of graph/chemistry functionality still to be
added here. The name stays `umol-graph`.

## Relationship to doc 104

Supersedes Phase B's "wedge + 2D depiction = symbolic descriptor" framing (104 §B, the
prose near "symbolic descriptor it consumes here" and the `mol.positions` /
`wedge_winding(pos, …)` pseudocode): the wedge path is geometric and lives in `umol-io` /
`umol-geometric-core`, not `umol-graph`. The `#T`/`#C` raise design is otherwise unchanged;
only its home crate and the wedge-geometry honesty change.
