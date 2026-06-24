# 129 — Crate organization and feature flags

Pre-release review of workspace layout. Four questions: (i) feature-gating the two
domains, (ii) a graph↔geometric conversion crate, (iii) a home for frozen hashes,
(iv) splitting `umol-shared`. This doc records the verified current state, a
boundary principle, ecosystem comparisons, and the open decisions for each.

## Current reality (verified)

14 members. Two crates have **no in-workspace consumers** — `umol-graph` and
`umol-geometric`. They are the two independent domain roots and share no edge.

Reverse-dependency map (X ← Y means Y depends on X):

```
umol-shared         <- geometric, graph, io, ast, params
umol-edn            <- graph, ast            umol-edn-macros   <- edn
umol-ast            <- graph, io             umol-ast-macros   <- ast
umol-graph-core     <- graph, ast            umol-perm         <- graph, io, ast
umol-io             <- graph                 umol-geometric-core <- io
umol-params         <- geometric, graph
umol-msym           <- geometric             umol-msym-sys     <- msym
```

Facts that reshape the questions:

- **`umol-msym-sys` is reachable only through `umol-geometric`.** A consumer that
  depends on `umol-graph` already compiles zero msym / C code. The graph/geometric
  isolation the feature-flag idea targets **already exists at the crate level**.
- **`umol-msym-sys` needs a C compiler (`cc`) + the `libmsym` git submodule, not
  clang/LLVM/bindgen.** The FFI in `umol-msym-sys/src/lib.rs` is hand-written;
  `build.rs` compiles 17 `libmsym/src/*.c` files via `cc`. Real cost (a C build +
  submodule), but no LLVM toolchain.
- **No `umol-graph` ↔ `umol-geometric` conversion code exists.** Grep both ways is
  empty; neither Cargo.toml names the other. The only cross-model bridge today is
  coordinate→stereo perception in `umol-io/src/table_ir/raise/utils.rs` (via
  `umol-geometric-core` 3D primitives). Coordinate→bond-order perception is a
  separate, isolated facility inside `umol-geometric` (`bond_perception.rs`).
- **`gboost` (RDKit `boost::hash_combine`) is duplicated** in
  `umol-graph-core/src/algorithms/refine.rs:483` and
  `umol-graph/src/fingerprint/pattern.rs:38` (byte-identical, intentional). It sits
  in `graph-core` only because the domain fingerprint schemes `EcScheme::{Morgan,
  RogersHahn}` live there (`refine.rs:441-480`) and `Morgan` calls it.
- **`xxh3` frozen seeds** also live in `graph-core/refine.rs` (`ALBATROSS_SEED`,
  `BULLFINCH_SEED`, the Rogers–Hahn ECFP seed). Separately, `umol-graph` uses
  `const_xxh3` to content-hash the valence registry/table — that is internal cache
  invalidation, a different category from reproducibility-frozen hashing.

## A boundary principle

A crate boundary should earn its keep through at least one of:

1. **Compile-cost isolation** — fences off a heavy or optional dependency (C/FFI
   builds, BLAS, large proc-macro trees). This is the only driver that a module
   cannot satisfy.
2. **Independent publishable API** — a unit you want versioned and depended on on
   its own.
3. **Cycle-breaking** — lets two crates that must not depend on each other
   interoperate through a third.

Grouping that is *only* about semantic tidiness is better served by **modules**, not
crates. Splitting for aesthetics alone produces the "microtome" the review wants to
avoid; the cost is real (more `Cargo.toml`s, version churn, feature plumbing).

## Ecosystem patterns

- **Umbrella facade + features.** `tokio` (one crate, `rt`/`net`/`fs`/`full`),
  `bevy` (umbrella `bevy` re-exports `bevy_*` subcrates; subsystems are features),
  `polars` (umbrella re-exports `polars-core`/`-io`/`-lazy` behind a large feature
  matrix). The subcrates are still published; the umbrella exists for one dependency
  line + one version, and gates heavy optional pieces behind `dep:`-features.
- **`*-sys` split for FFI.** `openssl`/`openssl-sys`, `git2`/`libgit2-sys`.
  `umol-msym`/`umol-msym-sys` already follows this correctly.
- **Cross-cutting grab-bag.** rust-analyzer's `stdx`, Servo's `util`, the common
  `*-common` / `*-support`. A small generic-infra crate is normal, not a failure —
  it stays acceptable by being leaf-level, dependency-light, and free of domain
  types.
- **Shared vocabulary types.** Conventionally `*-core` or `*-types` (e.g.
  `gix-hash`, `gix-object`; the rustc data-structure crates). `*-data` in the
  ecosystem usually means data *tables/assets*, which is why naming shared *types*
  `*-data` reads wrong.

## i. Feature-gating the domains

The toolchain-isolation goal is already met by the per-crate split: depend on
`umol-graph` and msym-sys never compiles. Feature flags `[graph, geometric]` only
become meaningful if a **single umbrella `umol` crate** is introduced that
re-exports both domains; then `umol = { features = ["geometric"] }` would gate the
msym-sys build via an optional dependency:

```toml
[features]
graph     = ["dep:umol-graph"]
geometric = ["dep:umol-geometric"]   # the only path that pulls umol-msym-sys
```

Decision is whether to ship an umbrella at all:

| Option | Consequence |
| --- | --- |
| **No umbrella** (publish `umol-graph`, `umol-geometric` separately) | Users pick the crate they need; isolation is automatic; no feature plumbing. Two dependency lines, two version numbers for a user who wants both. |
| **Umbrella `umol` + `graph`/`geometric` features** | One dependency, one version; discoverable; `geometric`-off cleanly omits the C build (Cargo `dep:`-feature unification handles it). Adds a crate and a feature surface to maintain. |

Note: an umbrella does **not** reduce the graph user's compile cost below the
no-umbrella case — it only matches it. Its value is ergonomic (single dep/version),
which is why polars/bevy ship one despite per-crate publishing.

## ii. Graph↔geometric conversion crate

Not speculative. Bond perception (`umol-geometric/src/bond_perception.rs`) **is** the
geom→graph conversion: it derives connectivity and bond orders from coordinates. It
had no target type until `MoleculeAst` existed; now it does. It currently sits inside
`umol-geometric`, which means any downstream graph processing of perception results
would couple `umol-geometric` to `umol-graph`. Moving it into a bridge crate that
depends on both keeps each domain free of the other (cycle-breaking, driver 3) and is
what makes the perception result consumable downstream.

The reverse direction (graph→3D *embedding*) is not in-tree, but an existing distance-
geometry implementation plus a port plan exist (doc 071). So the crate is small now
but has present code to relocate and a concrete planned addition — it earns its
boundary today.

Scope is **forced, not chosen**: the bridge holds the model-level conversions (bond
perception relocated from `umol-geometric`; doc 071 embedding later) and **cannot**
absorb `umol-io`'s coordinate→stereo perception. That perception runs inside `io`'s
raise step (`table_ir/raise/utils.rs`), and `umol-graph` depends on `umol-io` while the
bridge depends on `umol-graph` — so moving it up into the bridge closes a cycle:

```
io → geometric-graph → graph → io
```

It must stay at or below `io` in the dependency order (it already depends only on
`umol-geometric-core`, which sits below `io`). So perception is legitimately split:
coordinate→stereo stays in `io`; model→model bond perception lives in the bridge.

Name (decided): **`umol-geometric-graph`**, spelled out (no `geom` abbreviation;
`umol-geometric` / `umol-geometric-core` unchanged). `umol-conv` was rejected (vague,
abbreviated); `umol-bridge` (generic — only worth it for many model pairs, not
expected near-term); `umol-embed` named only the graph→3D half.

## iii. Home for frozen hashes

The xxh3 hashes are **not** one category, and most do not move. Verified state:

- `graph-core` holds **no** ECFP constant. `ECFP_SEED` already lives in `umol-graph`
  (`fingerprint/ecfp.rs:15`); `EcScheme::RogersHahn { seed }` just calls
  `xxh3_64_with_seed(_, seed)` with a domain-supplied seed. The only frozen *formula*
  in `graph-core` is **gboost** (the `EcScheme::Morgan` arm).
- **No in-core, non-test caller of color/circular refinement.** `subiso.rs`'s
  "refine" is Ullmann candidate-matrix refinement, unrelated. Every caller of
  `refine` / `circular_refine` / `EcScheme` / `RefinementXxh3Scheme` is a `umol-graph`
  fingerprint, so relocating `EcScheme` rewires no in-core code.
- `refine` is already generic (`refine<H: RefinementHash>`, `refine.rs:69`), but
  `circular_refine` takes a concrete `EcScheme` enum (`refine.rs:346`). That asymmetry
  is the only reason gboost is stuck in core.

gboost is categorically different from the xxh3 uses: it is a hand-rolled formula
whose sole purpose is bit-exact reproduction of an external tool (RDKit) — a domain
contract with no place in a graph-math crate. The xxh3 uses are invocations of a
standard hash library; the only "frozen" part is the **seed constants**, which are
nothing-up-my-sleeve values, not external contracts.

Principled end state (mirrors the existing `refine<H>` design):

- `graph-core` keeps the **mechanism**: generic `refine<H>`, a generic
  `circular_refine<R>` parameterized over an EC-recipe trait (the change that
  symmetrizes it with `refine`), the `RefinementHash` trait, and
  `RefinementXxh3Scheme` as a seed-agnostic xxh3 building block.
- `umol-graph` owns the **frozen recipes**: `EcScheme::{Morgan, RogersHahn}` impls,
  the gboost formula (one copy — kills the `pattern.rs` duplicate), and all named
  seeds (`ECFP_SEED` already here; WL seeds per the decision below), in one domain
  module, e.g. `umol-graph/src/fingerprint/hash.rs`.

**Decided.**

1. `circular_refine` becomes generic: `circular_refine<H: CircularRefinementHash>`
   (mirrors `refine<H: RefinementHash>`). `CircularRefinementAlgorithm::Ec` keeps only
   `radius`; the hash recipe is the type parameter.
2. **`graph-core` keeps only the trait seam + generic algorithms; everything else
   moves to `umol-graph::hash`.** Verified the generic algorithms reference the schemes
   solely through trait methods, so the whole xxh3 apparatus moves, not just the seeds:
   - Stays in `graph-core`: the `RefinementHash` and `CircularRefinementHash` traits;
     the generic `refine`/`refine_wl` and `circular_refine`/`circular_refine_ec`; and
     the algorithm-selection types `Refinement<C>`, `RefinementAlgorithm<H>`,
     `RefinementRounds`, `CircularRefinementAlgorithm`.
   - Moves to `umol-graph::hash`: `RefinementXxh3Scheme`, `RefinementWidth` (+
     `RefinementWidth64/128`), `RefinementAggregation`, the seeds
     (`ALBATROSS`/`BULLFINCH`, `ECFP_SEED` relocated from `fingerprint/ecfp.rs`), the
     gboost formula (one copy — kills the `pattern.rs` duplicate), and the `Morgan` /
     `RogersHahn` `CircularRefinementHash` impls.
   - Consequence: `graph-core` loses its `xxhash-rust` dependency entirely
     (`RefinementWidth64/128` are its only xxh3 users); the frozen-hash dependency
     surface concentrates in `umol-graph`. Core keeps algorithm test coverage via the
     existing `CountingScheme` trait-seam test; the xxh3-specific tests move with the
     scheme.
3. **No `umol-hash` crate.** Only `umol-graph` consumes any of this; `umol-geometric`
   has none. If a second consumer ever needs the seeds, move them then.
4. **Valence registry/table `const_xxh3` content hashes** are cache-keying, not a
   reproducibility contract; they stay put and do **not** join `umol-graph::hash`.

## iv. Splitting `umol-shared`

`umol-shared` today mixes two unrelated concerns:

| Chemistry vocabulary (types) | Generic infrastructure |
| --- | --- |
| `element`, `isotope`, `occupation`, `spin`, `configuration` (+ private `isotope_data`) | `error` (`UmolError` trait + domain error enums), `solution` (`Solution<T,C>` three-valued algebra), `units` (length/angle/time in atomic units) |

The discomfort is the mixing, and that the chemistry items are real shared *types*,
not a data store — so `umol-data` mis-names them. The principled split lifts the
**vocabulary** out into its own clearly-named crate (API-surface driver: e.g.
`umol-params` needs only `Element`, not `Solution`/`units`), leaving a small,
honestly-generic infra crate. The grab-bag does not disappear, but it becomes
*purely generic*, which is the form the ecosystem accepts.

**Decided: split into two crates.** A *vocabulary* crate and a generic *infra*
crate (`umol-utils`).

Placement:

| Goes to vocabulary crate | Goes to `umol-utils` |
| --- | --- |
| `element`, `isotope`, `occupation`, `spin`, `configuration` | `error`'s generic `UmolError` trait |
| `units` (length/angle/time) — sits with chemistry; e.g. half-life needs time units | `solution` (`Solution<T,C>`) |
| type-specific error enums (`ElementError`, `IsotopeError`, `OccupationError`, `SpinStateError`, `DataError`) — co-located with the types they describe | |
| private `isotope_data` table | |

Dependency direction: vocabulary → `umol-utils` (the chem error enums implement
`UmolError`). `umol-utils` is the leaf. `umol-utils` ends up small (an error trait + a
three-valued algebra) but both are genuinely cross-cutting generic infra — the
acceptable, purely-generic grab-bag.

`*-data` naming: reserve "data" for the actual tables (`isotope_data`, and
`umol-params`' radii / PPP parameters), never for the shared *types*.

Vocabulary-crate name — **`umol-chem`** (decided). Notes:

- `umol-core` is **out**: this workspace already uses `-core` to mean *math/algorithmic
  primitives* (`umol-graph-core`, `umol-geometric-core`); reusing it for chemistry
  vocabulary would collide with that established meaning.
- `umol-chem` — the hesitation ("everything here is chemistry") does not bite: the
  sibling crates are named by *model/role* (`graph`, `geometric`, `ast`, `io`,
  `params`), none competes for "chem", so `umol-chem` reads cleanly as the foundational
  chemistry vocabulary. Mild overlap: `umol-params` is also chemistry, but it is
  distinguished as *parameters* vs *types*.
- `umol-types` — fits the Rust convention (`*-types`, `sp-core`) but is neutral/SWE-
  flavored and drops the signal that these are chemical primitives.

## Decisions

All resolved.

- **ii scope** — the bridge **cannot** absorb `umol-io`'s coordinate→stereo perception:
  `io` uses it in raise, `graph` depends on `io`, and the bridge depends on `graph`, so
  relocating it closes the cycle `io → geometric-graph → graph → io`. It stays in `io`
  (depending only on `umol-geometric-core`). The bridge holds model→model conversions
  only.

Settled:

- **Names** — vocabulary crate `umol-chem`; infra crate `umol-utils`; bridge crate
  `umol-geometric-graph` (spelled out; `umol-geometric` / `umol-geometric-core`
  unchanged — no `geom` abbreviation).
- **i** — no umbrella crate for now (may revisit after colleague input; the per-crate
  split already isolates the msym C build).
- **ii structure** — create the bridge crate now (bond perception moves in; doc 071
  embedding to follow).
- **iii** — generic `circular_refine<H: CircularRefinementHash>`; all frozen
  seeds/formulas → `umol-graph::hash`; `graph-core` left frozen-constant-free and
  loses its `xxhash-rust` dependency.
- **iv** — `umol-shared` splits into `umol-chem` (vocabulary incl. `units` + the
  type-specific error enums) and `umol-utils` (`UmolError` trait, `Solution<T,C>`);
  vocabulary → `umol-utils`.
