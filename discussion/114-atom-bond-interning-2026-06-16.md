# 114 · Atom/bond interning (deferred)

Status: Deferred (forward-looking; not part of the 113 AST redesign)
Date: 2026-06-16
Relates: 113 (lazy canonicalization, `Canonical<T>`)

## Idea

A molecule holds many structurally-identical atoms/bonds (every H, every aromatic C),
and reaction networks generate huge numbers of near-identical molecules. Interning stores
each distinct atom/bond value **once**, shared by handle (`Rc<AtomAst>`), keyed on
canonical equality (`HashMap<Canonical<AtomAst>, Rc<AtomAst>>`) — so dedup/equality become
handle comparison and memory is shared across molecules. It matches EDN's value-identity
semantics. Interned values are immutable; "editing" is copy-on-write + re-intern.

## Evolution from the current storage model

`MoleculeAst` already shares structurally via `Arc<Vec<AtomAst>>` — whole-vector COW.
Interning is the same idea at finer granularity: `Vec<Rc<AtomAst>>` + a canonical-keyed
pool, giving **per-atom cross-molecule** sharing instead of whole-vector sharing (and
finer COW: an edit clones one atom, not the whole vector). It is additive — reads are
unaffected (views build from `&AtomAst`; `Rc` derefs; no pool needed for reads, only at
intern time), and the AST types and the lazy equality model do not change — so it can be
introduced cleanly later.

## Sticking point: `atoms_mut()`

Single-element `atom_mut(id)` survives: its internals become `Rc::make_mut` + re-intern on
the view guard's `Drop`, with call sites unchanged. But the bulk
`atoms_mut()`/`bonds_mut()` returning `impl Iterator<Item = &mut AtomAst>` cannot compose
with shared values — there is no per-item commit point to re-intern, and you cannot hand
out `&mut` into a shared `Rc`. Those must become closure/replace-based before interning
lands. Forward-looking discipline until then: do not expand the raw `&mut`-iterator APIs;
route edits through `*_mut(id)` or whole-value replacement.
