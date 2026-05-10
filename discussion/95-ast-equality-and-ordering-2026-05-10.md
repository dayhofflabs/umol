# AST equality and ordering — current model and open questions

**Status**: Active (current decision recorded; deeper redesign deferred)
**Date**: 2026-05-10
**Trigger**: Need for deterministic alias rendering in `Metadata`'s `BiMap`. Considered switching to `BiBTreeMap` (which requires `Ord` on the value type, `Box<AtomDsl>`).

## Current equality model on AST types

`AtomAst` and friends derive `PartialEq, Eq, Hash` via standard `#[derive]`. Equality is therefore purely **structural**: field-by-field, variant-by-variant comparison. The fields that carry expression payloads (`charge: ValueAst`, `isotope_mass: IsotopeAst`, etc.) inherit the same structural derive. `simplify()` is a separate explicit operation; equality never invokes it.

### Consequences

| Pair | Structurally equal? | Semantically equal? |
|---|---|---|
| `ValueAst::Lit(5)` vs `ValueAst::Expr(Expr::Lit(5))` | no | yes |
| `Expr::BinOp(Lit(2), Add, Lit(3))` vs `Lit(5)` | no | yes |
| `LitSet(vec![3,3])` vs `Lit(3)` | no | yes |
| `"C #h2"` parsed twice from the same source | yes | yes |
| `"C #h2"` vs `"C #h(2)"` parsed | no | yes |

The structural model surfaces in: `MoleculeAst::PartialEq`, hashing in `HashMap` keys, the alias bijectivity check (`BiMap::insert_no_overwrite`), proptest round-trip assertions, etc.

## Decision (2026-05-10)

**Keep structural equality for now.** Cascade `#[derive(Ord, PartialOrd)]` through the AST chain so `BiBTreeMap` becomes usable as the alias storage. Ordering inherits the same structural model: lexicographic by variant declaration order, then field order. The order is opaque (no chemistry meaning) but stable, which is all the alias case needs.

**Why now**: Switching `Metadata::atom_aliases` to `BiBTreeMap` removes the need for a render-time sort pass and gives free deterministic iteration. The Ord cascade is mechanical work; cost at runtime is negligible for the small-n alias case.

**Why not change equality semantics**: Moving to simplification-aware equality is a bigger redesign that affects every `==` site on AST types — not justified by the narrow alias-rendering need. Pulled out into the open questions below for separate consideration.

## Open questions (deferred)

### 1. Should AST equality be simplification-aware?

Three candidate models:

| Model | `Lit(5) == Expr(Lit(5))`? | Cost | Notes |
|---|---|---|---|
| **Structural** (current) | no | free | "Same atom written two ways isn't equal" |
| **Simplified-structural**: simplify both, compare structurally | yes | one simplify call per compare | Patterns with `Var("x")` still distinct from `Var("y")`; depends on canonicalization completeness; hash needs simplified form |
| **Semantic**: same admissible ground value set | yes | constraint solving | Decidable but expensive; conflates pattern with ground |

The current alias bijectivity check accepts `:a "C #h2"` and `:b "C #h(2)"` as distinct atoms because their parsed AtomAsts differ structurally — almost certainly not desired from a chemistry standpoint. Switching to simplified-structural would fix this case without going to full semantic equality.

### 2. Where else does structural equality leak through?

Inventory needed before committing to an equality redesign:
- `MoleculeAst::PartialEq` (used by proptest round-trip)
- `HashMap` / `HashSet` keys built from AST values
- Constraint dedup paths (currently structural via discriminant matching)
- Any registry / cache keyed by AST shape
- The `simplify_values` / `simplify` story — would equality call simplify, or expect callers to canonicalize first?

### 3. What's the "right" canonical form?

If equality is simplification-aware, `simplify()` must be the canonical form on both sides. This requires:
- `simplify()` being idempotent (it currently is for `ValueAst`)
- `simplify()` covering every variant where multiple structurally-distinct forms collapse to one (e.g., should `LitSet(vec![3,3])` simplify to `Lit(3)`? Currently no.)
- Handling unsimplifiable cases (`Var("x")`, `Bind`, `Ref`) without crashing.

### 4. Ordering semantics

Even with derived `Ord`:
- Lexicographic order on enums depends on declaration order. Reordering variants in a refactor silently changes ordering. Worth a comment near each ordered enum.
- Floating-point fields (none in AST today, but if they appear) break `Ord` derivation. Use `OrderedFloat` or similar at the boundary.

## Action items resolved by this doc

- [x] `Metadata::atom_aliases` switched from `BiMap` (= `BiHashMap`) to `BiBTreeMap`.
- [x] `#[derive(Ord, PartialOrd)]` cascaded through `AtomAst` chain (atom, value, expr, spin, constraints, …).
- [x] `render_atom_aliases` reverted to the simple iterate-and-emit form (BTreeMap iteration is sorted).

## Action items deferred

- [ ] Audit equality / hashing sites for any place where structural equality is unintentional.
- [ ] Decide on simplification-aware equality (whole-codebase change; needs separate discussion).
- [ ] Define and document the canonical form for each AST type (`simplify()` contract).

## See also

- `feedback_no_semantic_change_for_tests.md` — surfaced this discussion (proptest exposed BiMap iteration non-determinism; the wrong fix was patching the renderer to sort silently).
