# 187 — Multimolecular reaction convenience

Status: In Progress
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
molecule.react(reaction, config=None)
Molecule.react_all(reactants, reaction, config=None)
```

`None` selects the existing documented `ReactionApplicationConfig` default; the binding lowers the
result to the explicit Rust `SubstructureMatchConfig`. `Molecule.react` handles one reusable
molecule; the static `Molecule.react_all` accepts any iterable of molecules, including an empty
iterable. The latter parallels `Molecule.combine_all` rather than introducing a module-level
function or a wrapper for Python's built-in `list`. Neither operation eagerly materializes all
derivations or products.

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

## Staged implementation plan

### S0 — Owned reaction application (`umol-graph-ir`)

- **S0a — Owned application iterator.** Add the public `ReactionApplicationIter` with private
  fields and no public constructor. Its private construction path takes owned `Reaction` and
  `Molecule` snapshots plus the mandatory `SubstructureMatchConfig`, checks reaction-wide
  preconditions, eagerly enumerates correspondences, and retains normalized application state.
  Implement the one-shot iterator contract over `Result<ReactionDerivation, ApplyError>`: skip
  match-local rejection, yield a fatal application error once, then terminate. Add exact tests for
  successful enumeration, rejection, terminal failure, exhaustion, and snapshot independence.
  This is additive. [dep: none] **Done.** The owned operation-issued iterator stores the reaction,
  host, normalized deltas, correspondence iteration state, and terminal-failure state. The existing
  opaque `Reaction::apply` return delegates to it so its private constructor is exercised in normal
  library builds; S0b retains the public migration to the concrete named return type. Exact cases
  cover every stated iterator boundary.
- **S0b — `Reaction::apply` ownership migration.** Change `Reaction::apply` to return the concrete
  owned `ReactionApplicationIter`, cloning the borrowed reaction and host only to establish its
  documented snapshot semantics. Remove the borrowed closure iterator and migrate every Rust
  caller, test, rustdoc example, and benchmark without changing item order or error behavior. Add a
  generated equivalence property between the public iterator and explicit matching followed by
  `apply_at`, including match rejection and the terminal-error rule. This is breaking red-to-green.
  [dep: S0a] **Done.** `Reaction::apply` now returns the named owned iterator without a lifetime in
  its signature. Every workspace consumer compiles unchanged against the concrete return type. A
  256-case property compares the complete stream with explicit correspondence enumeration and
  `apply_at`; the separate malformed-reaction property retains exact fatal-error and terminal
  exhaustion coverage.

### S1 — Product-oriented Rust operation (`umol-graph-ir`)

- **S1a — Product iterator.** Add the public `ReactionProductsIter` as an owned lazy adapter over
  `ReactionApplicationIter`, with private fields and no public constructor. Each successful
  derivation yields the right-hand side's conservative connected components as `Vec<Molecule>`;
  application errors pass through unchanged and split correspondences are intentionally discarded.
  Cover zero, one, and multiple components, error forwarding, order preservation, and one-shot
  exhaustion. This is additive. [dep: S0b]
- **S1b — `React` capability.** Add and export `React::react` with the mandatory
  `SubstructureMatchConfig`, implemented for `Molecule` and `[Molecule]`. The molecule
  implementation delegates to owned reaction application. The slice implementation uses
  `Molecule::combine_all` in slice order and moves the fresh combined host directly into the
  application iterator without an intermediate clone; an empty slice uses the existing empty
  molecule. `Vec<Molecule>` participates through slice dereferencing, with no blanket
  `IntoIterator` implementation. Add exact cases and generated properties proving equality with
  the manual combine → apply → split pipeline for results, order, multiplicity, and both error
  channels. This is additive. [dep: S1a]

### S2 — Python ownership and convenience surface (`umol-py`)

- **S2a — Application iterator delegation.** Replace the Python correspondence/application driver
  with a thin non-constructible wrapper around the Rust `ReactionApplicationIter`; make
  `Reaction.apply` call the Rust operation directly after lowering the optional keyword-only
  `ReactionApplicationConfig`. Preserve the existing Python one-shot iteration, exception mapping,
  defaults, signatures, result ownership, order, and snapshot semantics. Remove the duplicated
  `apply_at` loop and its owned reaction/host/correspondence storage. This is breaking internally
  but green at subitem completion. [dep: S0b]
- **S2b — Product iterator binding.** Add a non-constructible Python wrapper around
  `ReactionProductsIter`, yielding owned `list[Molecule]` values lazily and translating the same
  precondition and item errors as reaction application. Test iteration, exhaustion, errors,
  ownership, component order, and multiplicity. This is additive. [dep: S1a, S2a]
- **S2c — `Molecule.react` and `Molecule.react_all`.** Expose
  `molecule.react(reaction, *, config=None)` and
  `Molecule.react_all(reactants, reaction, *, config=None)`. The static operation accepts any Python
  iterable, preserves its order, and accepts an empty iterable, following `Molecule.combine_all`'s
  conversion behavior. Both lower `None` through `ReactionApplicationConfig::default` and delegate
  to the Rust `React` implementations. Cover signatures, generators and lists, invalid iterable
  members, single/multiple/empty reactants, multiple matches, errors, and agreement with the manual
  Python combine → apply → split pipeline. This is additive. [dep: S1b, S2b]

### S3 — Documentation and closeout

- **S3a — Public contracts and examples.** Audit exports and rustdoc for `ReactionApplicationIter`,
  `ReactionProductsIter`, and `React`. State snapshot ownership, eager-correspondence/lazy-emission,
  operation-issued construction, error placement, ordering, multiplicity, and manual-pipeline
  equivalence. Update current Rust and Python examples plus the permanent data-type, nomenclature,
  and property-test guides where these contracts are user-facing; do not change the author-managed
  whitepaper. This is additive. [dep: S1b, S2c]
- **S3b — Repository-wide verification and status.** Run formatting, strict workspace clippy,
  workspace tests, the affected graph-IR property targets at the agreed larger case count, Python
  3.13 integration tests, and affected fuzz builds. Confirm the borrowed application iterator and
  duplicated Python application driver are gone, then mark this document completed and update
  `000-status.md`. This is additive. [dep: S3a]

The critical path is S0a → S0b → S1a → S1b → S2a/S2b → S2c → S3. No stage is deferrable: the
owned application iterator is required by both the Rust multi-reactant operation and the single
Python application path, and closeout includes removal of the duplicate Python driver.
