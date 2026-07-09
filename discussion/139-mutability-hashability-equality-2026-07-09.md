# 139 · Mutability, hashability, and equality across Rust and Python

Status: Informational (rationale + open items)
Date: 2026-07-09
Relates: 114 (atom/bond interning — the deferred home for the open items here), 137
(Python bindings — where the mutable surface was built), 113 (AST canonical equality
and lattice — the two-tier equality this rests on)

## Question

The Python bindings (137) moved to a mutable-container model: settable atom fields, a
mutable molecule of atoms, live constraint views that write through in place. Because
Python couples value-equality to unhashability, the mutable containers became unhashable.
Two questions this raises:

1. Did the Python side over-commit to a "mutate everything" metaphor, since Python objects
   are not immutable by default?
2. Is the Rust side's immutable / hashable / equality balance principled?

## The invariant everything turns on

`a == b` must imply `hash(a) == hash(b)`, and a hash must not change over an object's
lifetime. A type that defines *value* equality over mutable state therefore cannot be
hashable. Python enforces this at the language level: defining `__eq__` on a class sets
`__hash__` to `None` (unhashable) unless `__hash__` is also defined. `@dataclass(frozen=…)`
and PyO3's `#[pyclass(eq)]` both mirror that rule.

## Rust: value semantics with two-tier equality; hashing is orthogonal to mutability

The premise that "Rust imposes immutability" is not quite right. The value types
(`ValueAst`, `ElementAst`, `AtomConstraintAst`, `AtomConstraintsAst`) are **mutable** —
they have `with_*` builders and the container has `set`/`remove`/`update`. What Rust
imposes is **value semantics with a two-tier equality**, stated in the `ValueAst` doc
comment (`umol-ast/src/ast/value.rs`):

> Equality is lazy: derived `Eq`/`Hash`/`Ord` are structural ("same tree"); semantic
> equality is `Canonicalize::canonical_eq`, which compares canonical forms.

So:

- `==` / `Hash` / `Ord` are **derived and structural** — "same tree." Being the same
  derive, `Eq` and `Hash` are mutually consistent by construction; there is no
  equality/hash mismatch.
- Semantic equality (`LitSet{4}` denoting the same value as `Lit(4)`) is a **separate,
  explicit** `canonical_eq` (design 113). It is opt-in, for when meaning rather than
  representation is wanted.

The reason Rust keeps these types **hashable and mutable at once** is that the borrow
checker removes the footgun: you cannot obtain `&mut` to a value while it is a live
`HashSet`/`HashMap` key (the collection owns it and lends only `&`). In Rust, therefore,
**hashability and mutability are orthogonal** — no coupling is needed, and
`AstAtomConstraintsAst` legitimately derives `Hash` while also exposing `set`/`remove`.

The one exception is `MoleculeAst`: hand-written `PartialEq` that excludes `rings_cache`,
`Eq`, and no `Hash`. Adding a `Hash` that mirrors the `PartialEq` field set is the
fold-back recorded in 137.

**Assessment.** The Rust balance is principled and internally consistent. The outstanding
items are completeness, not balance: `MoleculeAst`'s missing `Hash`, and whether every
type that should carry `canonical_eq` uniformly does.

## Python: the coupling, and what it costs

Python cannot make the borrow-checker guarantee, so it couples value-equality to
unhashability. The bindings made the containers mutable with value-`__eq__`, so
`AtomAst` / `MoleculeAst` / `AtomConstraintsAst` are now **unhashable, with no frozen
counterpart**. The value *leaves* (`ValueAst`, …) kept their immutable, hashable form.

The standard-library convention is a mutable type paired with a frozen counterpart —
`list`/`tuple`, `set`/`frozenset`, `bytearray`/`bytes` (`dict` has no stdlib twin; a
`frozendict`, PEP 416, was rejected). Measured against that convention, the current Python
surface has the **frozen half only at the leaves**: the containers exist in the mutable
form alone.

## Assessment: the direction is right; the frozen/canonical molecule value is deferred

Two separate things sit inside "did we go too far":

- **The editing surface being mutable is correct.** In-place setters and a mutable
  container are what Python users expect (RDKit `RWMol`, `atom.SetFormalCharge(…)`). A
  builder-only immutable surface would be the un-Pythonic choice.
- **The absent frozen/canonical molecule *value* is the real open item.** The target
  domain is reaction networks (100k–1B+ nodes) over molecular graphs. If molecules are
  graph nodes, they want to be hashable, canonical, and dedup-able — an immutable molecule
  value. A mutable-only surface does not offer that form.

So the honest reading: the mutability *direction* is sound; shipping *no* hashable/canonical
molecule value is what would be incomplete — if the surface stopped here.

## Not locked in — the resolution belongs to interning (114)

Three reasons the mutable-container turn does not cement a "mutate everything" metaphor:

1. **Nothing hashable/immutable is precluded.** The Rust model already supports it
   (structural `Hash` + `canonical_eq` + `Canonicalize`); the binding has simply not
   surfaced a frozen molecule value. Adding one is additive.
2. **The consumer does not exist yet.** The molecule DSL and network layer are deferred, so
   a frozen/canonical form is future work, not a missed requirement.
3. **This is the interning question.** Interning (114) canonicalizes and deduplicates
   atoms/molecules, keyed on canonical equality (`HashMap<Canonical<AtomAst>, Rc<AtomAst>>`),
   so identity becomes value-identity and nodes hash by handle; 114 already states that
   interned values are immutable and "editing" is copy-on-write + re-intern. Deferring the
   handle-equality findings (137, findings 5/6) to interning parks this along the same axis.

The principled resolution, when the network layer lands, is **not** "make everything
immutable." It is to keep the mutable editing surface and **add** a hashable, canonical
molecule value, by one of:

- a `freeze()` / snapshot that returns a frozen hashable value (the `list`→`tuple`,
  `bytearray`→`bytes` twin);
- interned canonical molecules with identity hashing (114 directly);
- a canonical key (DSL string or canonical hash) that is the network node, leaving the
  object mutable.

That choice is the interning discussion's to make. The guardrail to set now: the
mutable-only container surface is the **incomplete current state, not the final metaphor**,
so it is not later mistaken for a settled decision.

## Adjacent open decision: Python `==` — structural or canonical

Python `==` currently mirrors Rust's **structural** equality, not `canonical_eq` — so
`ValueAst.Lit(4) == ValueAst.LitSet({4})` is `False` even though the two are semantically
equal. That is consistent with Rust, but a Python user may expect `==` to mean semantic
equality. Whether Python's `==` stays structural (with `canonical_eq` exposed separately)
or becomes the canonical one is a separate decision, entangled with the same "what is a
molecule value" question, and best settled alongside interning.
