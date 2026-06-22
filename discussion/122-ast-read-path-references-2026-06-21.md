# 122 — AST read-path zero-copy: constraint accessors and reference discipline (2026-06-21)

## Problem

Profiling the substructure matcher (doc 104 E6(4), doc 121) showed the per-candidate
predicate cloning values it only reads. The root is the constraint accessor API
(`constraint/atom.rs:481`):

```rust
pub fn valence(&self) -> ValueAst {
    match self.get(AtomConstraintKind::Valence) {
        Some(AtomConstraint::Valence(v)) => v.clone(),   // owned copy
        _ => ValueAst::Undetermined,
    }
}
```

`AtomConstraints::matches` calls ~12 such accessors on **both** self and target per
candidate, each returning an **owned** `ValueAst` (`v.clone()`). A read-only comparison
clones up to ~24 values per candidate. The accessors return owned because they were built
for `meet` — which *consumes* owned lattice elements to build a new one — and `matches`
reuses them.

This is a design flaw in load-bearing types (`AtomConstraints`/`BondConstraints`, and by
extension `AtomAst`/`MoleculeAst`): the read path is not zero-copy.

Already addressed (do not re-litigate): empty-pattern short-circuit in the collection
`matches`; cheap `matches` on `ValueAst`/`ElementAst`/`IsotopeMassAst`/`AromaticValenceAst`/
`MulticenterValenceAst`; `host_match_targets` borrows via `Cow`. Those remove most of the
clone for *unconstrained* patterns. This doc is the structural fix for *constrained*
patterns and the general principle.

## Principle

Separate **read** from **build**:

- Read paths (`matches`, queries, views) return **references**; they never own.
- Build paths (`meet`, `join`, `canonicalize`, the builder/`with_*`) own — cloning is
  legitimate there because they construct new values.

`matches` is read-only, so it must not allocate. `meet` allocating is fine and expected.

## Option A — reference-returning accessors, keep `Vec<AtomConstraint>` storage

`valence(&self) -> Option<&ValueAst>` (and likewise for every field); `None` ≡
`Undetermined`. `matches` compares by reference (an absent side compares against a
`const ValueAst::UNDETERMINED` or via `Option`). `meet` clones at the point it builds the
result. Storage and the `AtomConstraint` enum / `iter()` / `add()` / DSL / serde are
unchanged.

- **Pros:** small, localized; no change to the constraint model, DSL, or serde; preserves
  the sparse `SmallVec` storage (most atoms carry 0–2 constraints, so dense would waste
  memory); removes the read-path clone.
- **Cons:** keeps the linear scan per accessor (cheap at ≤ ~13 constraints, but O(n) ×
  fields if `matches` stays accessor-driven — see "matches shape" below); needs a const
  `Undetermined` or `Option`-threaded comparison.

## Option B — struct-of-fields storage

Replace `AtomConstraints(SmallVec<[AtomConstraint; N]>)` with a struct of typed fields
(`valence: ValueAst, total_valence: ValueAst, …, ring_memberships: Vec<(RingScope, ValueAst)>`,
or `Option` per field). O(1) field access by reference, no scan, no enum wrapping.

- **Pros:** O(1) zero-copy field access; no scan; the densest, fastest read path.
- **Cons:** large blast radius — `AtomConstraint` enum, `iter()`/`add()`, DSL parse/render
  (`[:v …]` etc.), EDN serde, canonicalize, and every constraint consumer change. Dense
  storage: ~13 `ValueAst` per atom even when unconstrained (most atoms), a memory
  regression vs. the sparse `SmallVec` for typical molecules. Reaction-network workloads
  (100k–1B nodes) make per-atom memory matter.

## `matches` shape (orthogonal to A/B)

Independent of storage, `matches` should be **pattern-driven**: iterate the *pattern's*
constraints (few, often none) and look each up in the target by reference, rather than
scanning all ~12 fields on both sides. With Option A this is `self.iter()` over pattern
constraints + `target.get(kind)` by reference; with Option B it's per-field. Either way it
makes the cost proportional to the pattern's constraints, not the fixed field count, and
composes with the empty-pattern short-circuit already in place.

## Broader scope

The same anti-pattern likely recurs wherever a view/accessor returns owned data for a
read. Worth auditing under the same principle: `AtomView`/`BondView` accessors, the
`MoleculeAst` query surface, and `derive_constraints` (which builds an owned
`AtomConstraints` per host atom — a build path, but invoked from a read path; the
pattern-driven gate already skips it when unneeded). The goal is that a substructure match
over an immutable `MoleculeAst` is provably zero-copy on the read side.

## Open questions (for the design discussion)

1. Accessor return type: `Option<&ValueAst>` vs `&ValueAst` defaulting to a const
   `Undetermined`?
2. Storage: Option A (sparse `SmallVec`, ref accessors) or Option B (struct-of-fields)?
   The sparse-vs-dense memory tradeoff under reaction-network scale is the crux.
3. Confirm the read/build split: `meet`/`join`/`canonicalize` keep owning; only the read
   accessors + `matches` go to references.
4. Scope: just the constraint collections, or a broader reference-discipline pass over the
   `AtomAst`/`MoleculeAst` view/accessor surface?
5. Does `ValueAst` itself need anything, or is it fine (the issue is the API around it, not
   the enum — `Lit` clones cheaply; only `LitSet`/`Term`/`Predicate` allocate, and the read
   path will stop cloning them)?

No implementation until these are settled.
