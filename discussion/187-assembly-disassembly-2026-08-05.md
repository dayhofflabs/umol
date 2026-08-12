# 187 — Multimolecular reaction convenience

Status: Proposed
Date: 2026-08-05
Relates: [005](005-mutability-2025-03-01.md),
[151](151-python-molecule-workflows-2026-07-13.md)

## Motivation

A reaction applies to one `Molecule`, while a chemical reaction may consume several disconnected
reactant molecules and produce a different number of disconnected product molecules. The primitive
graph-IR operations already express the required pipeline:

1. combine the reactants by disjoint union;
2. apply the reaction to the combined host; and
3. split each resulting right-hand side into connected components.

This work adds one product-oriented convenience operation for that pipeline in Rust and Python. It
does not replace `Reaction::apply`: callers that need the complete `ReactionDerivation`, including
the host-to-product correspondence, continue to use the lower-level operation.

The convenience operation also exposes a problem in the current ownership model. Rust returns an
opaque iterator borrowing both the reaction and host. Python cannot retain those Rust borrows and
therefore implements a second application iterator over owned reaction and host snapshots. A
multi-reactant operation cannot return the Rust iterator at all because it would borrow the temporary
combined host. Reaction application, its Python binding, and the combined operation must therefore
be designed as one ownership unit rather than accumulating another adapter around the borrowed
iterator.

## Settled API

### Owned reaction application

`Reaction::apply` returns a concrete owned one-shot iterator:

```rust
pub fn apply(
    &self,
    host: &Molecule,
    match_config: SubstructureMatchConfig,
) -> Result<ReactionApplicationIter, ApplyPreconditionError>;
```

`ReactionApplicationIter` owns the reaction and host snapshots, the eagerly enumerated
correspondences, and any normalized application state needed to construct derivations lazily. Its
fields are private and it has no independent public constructor. Its item type remains:

```rust
Result<ReactionDerivation, ApplyError>
```

This preserves the existing eager-correspondence/lazy-derivation compromise. Constructing the
iterator checks reaction-wide preconditions before returning. Match-local rejection remains an
internal skip. A non-rejection application failure is yielded once and terminates the iterator; it
is not silently discarded.

The iterator has snapshot semantics. Mutating a source `Molecule` or Python `Reaction` after
application starts cannot change that application. `Molecule` clones share their large tables
through `Arc`, so taking the Rust snapshots is predominantly reference-count incrementing rather
than a deep copy. The implementation should nevertheless move a freshly combined host directly into
the iterator rather than clone it merely to satisfy an internal borrowed helper.

Python wraps this graph-IR iterator instead of maintaining a separate correspondence loop and
calling `apply_at` itself. Python retains its current one-shot iterator behavior.

### Product-oriented convenience trait

Rust adds the capability trait:

```rust
pub trait React {
    fn react(
        &self,
        reaction: &Reaction,
        match_config: SubstructureMatchConfig,
    ) -> Result<ReactionProductsIter, ApplyPreconditionError>;
}
```

`React` is implemented for `Molecule` and `[Molecule]`. A `Vec<Molecule>` uses the slice
implementation through dereferencing. There is no blanket implementation for arbitrary
`IntoIterator`: `react` is a non-consuming operation over reusable reactants, while a blanket
iterator implementation would introduce consuming receiver and coherence complications. A caller
with another iterable may collect it before applying the operation.

`ReactionProductsIter` is an owned lazy adapter over `ReactionApplicationIter` with item type:

```rust
Result<Vec<Molecule>, ApplyError>
```

For one molecule, `react` applies the reaction and splits every successful derivation's right-hand
side. For a slice, it first combines the reactants in slice order, applies the reaction to that
combined host, and splits every successful right-hand side. The two implementations therefore have
one result contract. The lower-level `Reaction::apply` remains available when discarding the
derivation and its correspondence would lose information the caller needs.

The Rust match config is mandatory. `umol-graph-ir` must not acquire a hidden matching-algorithm
default for this convenience operation.

### Python surface

Python exposes the same owned lazy application and product-oriented convenience semantics. The
application config is keyword-only and optional:

```python
reactants.react(reaction, config=None)
```

`None` selects the existing documented `ReactionApplicationConfig` default; the binding lowers the
result to the explicit Rust `SubstructureMatchConfig`. Both a single `Molecule` and an iterable of
`Molecule` values must be accepted without eagerly materializing all derivations or products. The
exact Python spelling for the iterable entry point remains to be selected because Python has no
trait implementation on the built-in `list` type.

## Semantic contract

- Combining is disjoint union only. It performs no gluing, matching, resolution, or validation
  beyond the contracts of the existing operations.
- Input order determines the dense id ranges in the combined host. It is preserved as an operational
  input; the convenience layer does not sort or canonicalize reactants.
- Product multiplicity is structural. Joining reactants may reduce the number of connected
  components, while cleavage may increase it. The result is consequently one `Vec<Molecule>` per
  successful reaction match rather than one product per input reactant.
- Splitting uses the existing conservative `Molecule::split` semantics: any supported relation keeps
  its participants in one component.
- Match enumeration order and per-product component order are inherited from `Reaction::apply` and
  `Molecule::split`. The convenience layer neither deduplicates nor reorders them.
- An empty reactant slice forms the existing empty combined molecule. Whether a reaction applies to
  it is determined by ordinary reaction preconditions and matching; the convenience operation does
  not add a special rejection.
- `ReactionApplicationIter` and `ReactionProductsIter` are operation-issued values. Callers do not
  construct them independently or supply alternate hosts or correspondences after construction.
- Reaction-wide failures remain the outer `ApplyPreconditionError`. Failures arising while realizing
  a selected match remain `ApplyError` iterator items. Ordinary match rejection is not an error item.
- The result for a molecule slice is definitionally equivalent to the manual
  `Molecule::combine_all` → `Reaction::apply` → `ReactionDerivation::rhs` →
  `Molecule::split` pipeline, including match order, product-component order, and errors.

## Scope boundaries

This work does not add `Reactants` or `Products` wrapper types, retain a second Python application
driver, infer atom mappings, canonicalize products, deduplicate symmetric matches, or change the
matching algorithms. It does not remove the explicit combine/split APIs or the derivation-oriented
`Reaction::apply` operation.

## Open question

- Select the Python spelling for the iterable-reactant entry point while keeping the single-molecule
  method and iterable operation semantically parallel.
