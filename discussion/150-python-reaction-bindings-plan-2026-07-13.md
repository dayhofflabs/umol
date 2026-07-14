# 150 · Python bindings for `Deltas` and `ReactionAst` (plan)

Status: **ACTIVE — S4 complete; S5a is next**
Date: 2026-07-13
Relates: 131–135 (reaction semantics and implementation), 137 (Python binding
strategy), 139 (Python mutability/equality), 140 (entity-binding template)

## Scope

Bind the resolved reaction data model to Python: the complete `Delta` algebra,
its `Deltas` collection, and `ReactionAst`. The core deliverable is a reaction
that Python can construct structurally, parse/render, inspect and mutate through
its owned `lhs` and `deltas`, normalize, reverse, build from two sides, compose,
and apply to host molecules. `ReactionDerivation` is a fully owned Rust value,
so the Python surface wraps the same value model directly.

This follows doc 140's rule: bind a complete semantic slice, not only the easy
variants. In particular, `Delta::Constraint` makes the currently unbound
molecule-level constraint tree a prerequisite. Omitting it would make parsed
reactions readable only until they contained a molecule constraint.

Out of scope for the core deliverable: `ReactionSpanAst` and its per-entity
`EntitySpan<T>` views, reaction metadata/alias preservation (`ReactionDsl`),
network/interning types, and the general transaction/edit vocabulary.

## Rust facts this rests on

- `Deltas(Vec<Delta>)` is an owned, ordered input collection. It exposes
  `new`/`len`/`is_empty`/`iter`/`iter_mut`/`push`; `FromIterator`; and consuming
  `Canonicalize`. Canonicalization folds operations per entity, rejects
  discontinuous or structurally inconsistent chains, then sorts the result into
  a unique normal form. The stored input order is therefore meaningful before
  canonicalization but not afterward.
- `Delta` has nine families: atom, bond, dative, aromatic, multicenter,
  noncovalent, stereo atom, stereo bond, and molecule constraint. Every delta is
  closed under consuming `inverse()`.
- The six non-stereo entity deltas share four operation shapes: `Add`, `Remove`,
  `ModifyField`, `ModifyConstraint`. Structural participants occur only on
  `Add`/`Remove`; ids and chemistry values are separate.
- Stereo atom/bond add three relative operations (`Apply`, `Swap`, `Mirror`) and
  carry `StereoKind` context on relative operations and constraint changes.
- `ModifyField` delegates to eight public old/new field-change enums. All their
  leaf values are already bound by doc 140.
- `ConstraintDelta::{Add, Remove}` carries the recursive top-level `Constraint`
  tree, not one of the already-bound inline entity-constraint enums. Binding it
  pulls in `RelationalConstraint`, `MoleculeConstraint`, `SubPatternAnchor`, and
  the recursive `Constraint::{And, Or, Not}` representation.
- `ReactionAst { lhs, deltas }` is an owned value. Its directly usable operations
  are `new`, `from_sides`, `canonicalize`, `reverse`, `compose`, and its EDN
  `FromStr`/`Display` surface. `to_reaction_span` needs a span wrapper.
- `apply_at` takes a `MoleculeCorrespondence` and returns an owned
  `ReactionDerivation`; `apply` consumes a `SubgraphIsomorphismAlgorithm` and
  iterates owned derivations. Each derivation owns `lhs`, `rhs`, and its
  correspondence.

### Rust ownership model (settled)

`ReactionDerivation { lhs: MoleculeAst, rhs: MoleculeAst, comap }` is
non-generic and fully owned. `apply_at` clones the host once when it constructs
a successful result; `reverse` and `chain` return ordinary owned values. The
Python binding wraps this type directly.

## Settled Python surface

### Structural delta vocabulary

Use the native PyO3 complex-enum mechanism established in docs 137/140. Python
gets one structural mirror per public Rust enum:

- field changes: `AtomFieldChange`, `BondFieldChange`,
  `DativeBondFieldChange`, `AromaticSystemFieldChange`,
  `MulticenterBondFieldChange`, `NoncovalentBondFieldChange`,
  `StereoAtomFieldChange`, `StereoBondFieldChange`;
- entity deltas: `AtomDelta`, `BondDelta`, `DativeBondDelta`,
  `AromaticSystemDelta`, `MulticenterBondDelta`, `NoncovalentBondDelta`,
  `StereoAtomDelta`, `StereoBondDelta`, plus `ConstraintDelta`;
- the sum type `Delta`.

Variants and payload order mirror Rust verbatim and support construction,
attribute reads, `isinstance`, and pattern matching. Ids are bare Python `int`s,
matching the completed entity bindings (the obsolete Python id-newtype decision
in early doc 137 does not return). Structural participant payloads preserve their
Rust representation exactly (delta equality is structural, so sorting or
set-conversion here would be lossy):

- fixed arrays: `tuple[int, int]`;
- `Vec<AtomId>` participants: `list[int]`;
- stereo ligand frames: ordered `list[StereoLigand]`;
- dative roles remain separate `donors` and `acceptor` fields.

The mirrors have value equality, are unhashable, and freeze the variant and its
direct fields. Each exposes `inverse() -> Self`. Nested hold-the-value entity AST
payloads remain mutable Python objects, as the same Rust values are mutable behind
`&mut`; they are nevertheless snapshotted on construction, so mutating the source
`AtomAst`/`BondAst` later does not mutate an already-built delta.

### `Deltas`

`Deltas` is a hold-the-value, mutable, value-equal and unhashable pyclass:

- `Deltas(entries=())`, accepting any iterable of `Delta`;
- `len(deltas)`, `bool(deltas)`, integer `deltas[i]` with negative-index
  normalization, and iteration, all returning detached `Delta` snapshots;
- `append(delta)` and `extend(iterable)` as the Python spelling of Rust `push`
  and repeated `push`; no item assignment/removal, because the Rust API exposes
  neither and arbitrary mid-chain mutation obscures delta continuity;
- `canonicalize() -> Deltas`, non-mutating because the Rust operation consumes
  by value; contradiction raises `ContradictionError`.

There is deliberately no `Deltas.inverse()`: reversing a whole reaction also has
to remap ids as entities move between lhs and rhs frames, and that operation
already lives on `ReactionAst::reverse`.

### Molecule-level constraints (prerequisite surface)

Bind `SubPatternAnchor`, `RelationalConstraint`, `MoleculeConstraint`, and the
recursive `Constraint` as structural mirrors, then a mutable `Constraints`
container and live `ConstraintsView`. Wire `MoleculeAst.constraints` and a
`constraints=` keyword on `MoleculeAst.from_parts`. This makes the full payload
of `ConstraintDelta` constructible and makes a parsed reaction's constrained LHS
observable and editable.

This is a prerequisite, not a general redesign of constraints. The container
uses the existing sequence semantics (`len`/integer indexing/iteration,
`append`/`clear`; live write-through on the molecule). It does not pretend the
recursive tree is the keyed mapping used by inline entity constraints.

### `ReactionAst`

The Python wrapper is an owned **component facade**, not a detached-snapshot
getter over a held Rust `ReactionAst`:

- construction snapshots the input `MoleculeAst` and `Deltas` into private
  `Py<MoleculeAst>` and `Py<Deltas>` components;
- `.lhs` and `.deltas` return those same owned components, so
  `reaction.lhs.atoms[0].charge = ...` and `reaction.deltas.append(...)` both
  write through to the reaction;
- Rust operations rebuild an `AstReactionAst` from component snapshots and wrap
  returned values back into fresh owned components.

This preserves the mutable-container direction of doc 139 without the misleading
behavior in which a nested mutation silently edits a temporary clone. It also
preserves Rust's consuming constructor semantics: mutating the objects passed to
`ReactionAst(lhs, deltas)` afterward does not alias the reaction. Conversion code
must snapshot all right-hand-side components before taking a mutable Python
borrow, following doc 140's self-alias regression rule.

Surface:

- `ReactionAst(lhs=None, deltas=None)` (`None` means an empty component, avoiding
  mutable Python defaults) and structural `==` (unhashable);
- live `.lhs` and `.deltas` getters plus whole-component setters that snapshot;
- `ReactionAst.parse(text)`, `str(reaction)`, and constructor-style `repr`;
- `canonicalize() -> ReactionAst` and `reverse() -> ReactionAst`, mapping
  contradictions to `ContradictionError`;
- `from_sides(lhs, rhs, atom_pairs) -> ReactionAst`, where `atom_pairs` is an
  iterable of `(lhs_atom, rhs_atom)` integer pairs and side sizes are inferred;
  Rust still receives a `Correspondence<NodeId>` internally;
- `compose(other, scope=CompositionScope.RcAnchored) -> list[ReactionAst]`, with
  a simple `CompositionScope` mirror.

`from_sides` takes pairs rather than exposing a one-use atom-correspondence
wrapper. A reusable `MoleculeCorrespondence` wrapper lands with application,
where all eight correspondence families become observable.

### Errors and verification

Add `ContradictionError` now; do not fold normalization failure into
`ValueError` or `ParseError`. The application stage adds a typed `ApplyError`
mapping with the dangling atom id preserved. Every subitem below
carries Rust conversion/unit tests and Python construction, pattern-match, and
behavior tests. A stage is green only when `cargo test -p umol-py`, the pytest
suite, clippy, and formatting pass.

## Staged implementation plan

### S0 — Owned `ReactionDerivation` (Rust)

- **S0a — DONE — owned derivation value**
  (`umol-ast/src/ast/reaction_derivation.rs`): `ReactionDerivation` owns `lhs`,
  `rhs`, and `comap`; its constructor, accessors, `reverse`, and `chain` use
  ordinary owned signatures. It derives `PartialEq`/`Eq`. Unit tests cover
  owned independence, reverse, chain, and abstraction back to `ReactionAst`.
  **Breaking.** `[dep: —]`
- **S0b — DONE — producer and caller migration** (`umol-ast/src/ast/reaction.rs`,
  `compose.rs`, tests, and workspace callers): `apply_at` clones the matched host
  into each successful derivation; `apply_at` and `apply` return owned
  `ReactionDerivation` values while preserving match enumeration and DPO
  behavior. Focused 9, complete `umol-ast` 4,589, workspace check, and clippy all
  pass. **Breaking migration (red→green).** `[dep: S0a]`
- **S0c — DONE — API documentation cleanup** (affected Rust API and design docs):
  documentation describes `ReactionDerivation` as an owned value and matches the
  implemented signatures. Documentation-only. **Additive (green).** `[dep: S0b]`

### S1 — Complete the constraint payload foundation

The Python constraint bindings live under the singular `umol-py/src/constraint/`
module, grouped by domain: `atom`, `bond`, `dative`, `aromatic`, `multicenter`,
`noncovalent`, `stereo`, `molecule`, and `ring`. Entity modules retain the owned
entity values and molecule-backed entity views. Aromatic- and multicenter-valence
payloads are atom constraints and therefore live in `constraint/atom.rs`; ring
scope and membership payloads live in `constraint/ring.rs`.

- **S1a — DONE — constraint leaves and anchors**
  (`umol-py/src/constraint/molecule.rs`): the following exact type closure is
  implemented in dependency order:

  1. **DONE — `SubPatternAnchor` mirror struct** — eight `(target, pattern)` id-pair
     collections: `atoms`, `bonds`, `dative_bonds`, `aromatic_systems`,
     `multicenter_bonds`, `noncovalent_bonds`, `stereo_atoms`, and
     `stereo_bonds`. Constructor/getters use bare `int` pairs. Its
     `from_ast`/`to_ast` bridge is required by `MoleculeConstraint::SubPattern`.
  2. **DONE — `RelationalConstraint` structural mirror enum** — all 31 Rust variants,
     grouped by referenced entity:
     - dative (8): `DativeBondDonors`, `DativeBondDonor`,
       `DativeBondContainsAllDonors`, `DativeBondAllDonors`,
       `DativeBondAnyDonor`, `DativeBondAcceptor`,
       `DativeBondAcceptorSatisfies`, `DativeBondParallels`;
     - aromatic (5): `AromaticSystemAtoms`, `AromaticSystemContains`,
       `AromaticSystemContainsAll`, `AromaticSystemAllAtoms`,
       `AromaticSystemAnyAtom`;
     - multicenter (5): `MulticenterBondAtoms`, `MulticenterBondContains`,
       `MulticenterBondContainsAll`, `MulticenterBondAllAtoms`,
       `MulticenterBondAnyAtom`;
     - noncovalent (3): `NoncovalentBondEnds`,
       `NoncovalentBondContains`, `NoncovalentBondEndsSatisfy`;
     - stereo atom (5): `StereoAtomSite`, `StereoAtomContains`,
       `StereoAtomLigands`, `StereoAtomAllLigands`,
       `StereoAtomAnyLigand`;
     - stereo bond (5): `StereoBondSite`, `StereoBondContains`,
       `StereoBondLigands`, `StereoBondAllLigands`,
       `StereoBondAnyLigand`.
  3. **DONE — `MoleculeConstraint` structural mirror enum** — all five Rust variants:
     `ChargeSum { atoms, sum }`, `SpinSum { atoms, spin }`,
     `BondOrderSum { bonds, sum }`, `Connected { atoms }`, and
     `SubPattern { anchor, pattern }`. `None` on an atom/bond subset continues
     to mean the whole molecule; `Some([])` remains a distinct empty subset.

  **Reused bound children:** `AtomConstraintAst` for the 12 predicate-bearing
  relational variants (13 payload values because `EndsSatisfy` carries two),
  `ValueAst`, `SpinStateAst`, and `MoleculeAst`. The eight Rust id newtypes map
  to bare Python `int`; S1a adds no Python id classes. Shared conversion support
  is `into_py_variant`/`variant_repr`. `MoleculeAst::from_inner` is available at
  runtime so `MoleculeConstraint::SubPattern` can wrap a nested pattern. The
  `constraint::molecule` compiles, while native-module registration and
  `python/umol/__init__.py` exports remain deferred to S1d.

  **Implemented verification:** two anchor tests (constructor/getters and one
  populated round trip covering all eight families), one round-trip case for
  every relational variant (31), and six molecule-constraint cases covering all
  five variants plus the `None` versus `Some([])` distinction. Nested atom
  predicates, the ordered noncovalent predicate pair, spin/value payloads, and a
  nested subpattern all round-trip. The 39 focused tests, `cargo check`, clippy,
  rustfmt, and `git diff --check` pass. Python import-level construction,
  field-access, and pattern-match tests land with registration in S1d.
  **Additive (green).** `[dep: —]`
- **S1b — DONE — recursive `Constraint` mirror** (`constraint/molecule.rs`): add one
  native complex-enum mirror with the following exact 13-variant closure:

  1. **DONE — Ordinary entity leaves (6)** — the entity id is a bare Python `int`; the
     second field reuses the already-bound structural constraint enum:
     - `Atom(int, AtomConstraintAst)` ↔
       `Constraint::Atom(AtomId, ast::AtomConstraintAst)`;
     - `Bond(int, BondConstraintAst)` ↔
       `Constraint::Bond(BondId, ast::BondConstraintAst)`;
     - `DativeBond(int, DativeBondConstraintAst)` ↔
       `Constraint::DativeBond(DativeBondId, ast::DativeBondConstraintAst)`;
     - `AromaticSystem(int, AromaticSystemConstraintAst)` ↔
       `Constraint::AromaticSystem(AromaticSystemId,
       ast::AromaticSystemConstraintAst)`;
     - `MulticenterBond(int, MulticenterBondConstraintAst)` ↔
       `Constraint::MulticenterBond(MulticenterBondId,
       ast::MulticenterBondConstraintAst)`;
     - `NoncovalentBond(int, NoncovalentBondConstraintAst)` ↔
       `Constraint::NoncovalentBond(NoncovalentBondId,
       ast::NoncovalentBondConstraintAst)`.
  2. **DONE — Stereo entity leaves (2)** — preserve the extra geometry discriminator;
     ids remain bare `int`, `StereoKind` reuses the bound value enum, and the
     final field reuses the corresponding structural constraint enum:
     - `StereoAtom(int, StereoKind, StereoAtomConstraintAst)` ↔
       `Constraint::StereoAtom(StereoAtomId, ast::StereoKind,
       ast::StereoAtomConstraintAst)`;
     - `StereoBond(int, StereoKind, StereoBondConstraintAst)` ↔
       `Constraint::StereoBond(StereoBondId, ast::StereoKind,
       ast::StereoBondConstraintAst)`.
  3. **DONE — Aggregate leaves (2)** — reuse the S1a mirrors directly:
     - `Relational(RelationalConstraint)` ↔
       `Constraint::Relational(ast::RelationalConstraint)`;
     - `Molecule(MoleculeConstraint)` ↔
       `Constraint::Molecule(ast::MoleculeConstraint)`.
  4. **DONE — Recursive combinators (3)** — recursive children are native variant
     instances, not base-enum objects:
     - `And(list[Constraint])` ↔ `Constraint::And(Vec<Constraint>)`;
     - `Or(list[Constraint])` ↔ `Constraint::Or(Vec<Constraint>)`;
     - `Not(Constraint)` ↔ `Constraint::Not(Box<Constraint>)`.

  The Rust representation is therefore `Py<...>` for every nested complex-enum
  child, `Vec<Py<Constraint>>` for `And`/`Or`, and `Py<Constraint>` for `Not`.
  `from_ast` must call `into_py_variant` at every such edge—including entity,
  relational, molecule, and recursive children—so Python receives the concrete
  variant subtype with `_0`/`_1`/`_2`, `__match_args__`, and `match` support;
  `Py::new` would create an unusable base-enum instance. `to_ast` recursively
  borrows each child and rebuilds owned Rust values. Equality delegates to the
  complete Rust tree; `variant_repr` uses arity 2 for ordinary leaves, arity 3
  for stereo leaves, and arity 1 for aggregate/combinator variants. The mirror
  is value-equal and unhashable.

  **Implemented verification:** one round-trip row covers each of the exact 13
  variants, with distinct stereo kinds, representative child constraints, and
  explicit empty `And([])` and `Or([])` cases. A deep tree combines entity,
  relational, and molecule leaves beneath `And`/`Or`/`Not` and exercises the
  concrete Python variant fields, exact recursive repr, equality, and class
  pattern matching from Rust-side PyO3. The focused S1a/S1b suite has 53 cases;
  import-level Python coverage remains in S1d when the type is registered.
  **Additive (green).** `[dep: S1a]`
- **S1c — DONE — molecule constraint container** (`constraint/molecule.rs`): follow the
  entity constraint containers directly: the same value-container/live-view
  split, shared snapshot iterator, `Update`/resolved-update pair, whole-container
  `Arg`, and resolve-before-write discipline. The storage is an ordered sequence,
  so the principal semantic difference is that `update` appends every resolved
  entry in order and preserves duplicates instead of deduplicating by key.

  1. **DONE — Shared conversion and iteration support** — add the direct sequence
     analogue of the entity containers' iterator builder: one shared
     negative-index normalizer and one `ConstraintIter` pyclass holding an owned
     snapshot of concrete Python `Constraint` variants. Both container forms use
     this same code. Test positive, negative, and out-of-range normalization plus
     snapshot iteration. **Additive (green).** `[dep: S1b]`
  2. **DONE — Owned `Constraints` value** — add a mutable, value-equal, unhashable
     `#[pyclass]` around `ast::Constraints`, in the same hold-the-value shape as
     `AtomConstraintsAst`/`BondConstraintsAst` and their peers. Expose the
     sequence equivalents of their operations: constructor from entries, exact
     recursive repr, `append`, `clear`, `__len__`, integer `__getitem__`, and
     `__iter__`; return detached concrete variants. Provide crate-private
     `inner`, `inner_mut`, and `from_inner` bridges. Test empty/populated
     construction, structural equality and repr, append/clear, negative indexing,
     bounds errors, and detached iteration. **Additive (green).** `[dep: S1c.1]`
  3. **DONE — Live `ConstraintsView`** — add the direct analogue of the entity
     `*ConstraintsView` pyclasses, backed here by `Py<MoleculeAst>` and using
     `read`/`with_mut` over `MoleculeAst::{constraints,constraints_mut}`. Mirror
     the owned container's operations and reuse the same index and iterator
     helpers. Test write-through, observation of later molecule changes,
     negative indexing parity, and detached iteration. **Additive (green).**
     `[dep: S1c.1]`
  4. **DONE — Entity-style RHS-first inputs and update** — add
     `ConstraintsUpdate::{Container, View, Entries}`,
     `ResolvedConstraintsUpdate`, and `ConstraintsArg::{Container, View}` with
     the same responsibilities as the entity-container counterparts. Add
     `update` to both targets. Resolve the complete RHS before the target write
     borrow; `ResolvedConstraintsUpdate::apply` differs only in applying every
     entry with `push`, retaining order and duplicates. Test all three RHS forms
     and both self-alias cases: owned-from-self and view-from-the-same-view each
     append exactly one snapshot without a PyO3 double-borrow panic. **Additive
     (green).** `[dep: S1c.2, S1c.3]`

  **Implemented verification:** S1c.1 has eight index/iterator cases. S1c.2 has
  15 owned-container cases covering construction, equality, repr, mutation,
  indexing, detached iteration, and conversion. S1c.3 has nine live-view cases
  covering count repr, duplicate-preserving append, clear, observation of later
  molecule mutation, indexing, and detached iteration. S1c.4 adds six cases for
  both whole-container inputs, all three update RHS forms on both targets, and
  owned/view self-aliasing; updates preserve order and duplicates. The focused
  S1 suite has 91 cases and the complete `umol-py` library suite has 632; clippy,
  rustfmt, and `git diff --check` pass. Native registration and Python
  import-level coverage remain in S1d. `[dep: S1b]`
- **S1d — DONE — `MoleculeAst` wiring** (`umol-py/src/molecule.rs`, `lib.rs`,
  `python/umol/__init__.py`): expose `mol.constraints`, add the keyword-only
  `constraints=()` input to `from_parts`, and register/export the S1 classes.
  The property returns a live `ConstraintsView` and accepts either `Constraints`
  or another live view on assignment, using the same snapshot-before-write
  discipline as the entity containers. `ConstraintIter` remains an internal
  implementation type. Existing calls remain valid.

  **Implemented verification:** three Rust wiring tests cover `from_parts`, live
  store access, whole-container replacement, and self-view assignment. Eight
  Python tests cover public imports, anchor fields, nested class-pattern matching,
  duplicate-preserving sequence behavior, constructor input, live mutation, and
  assignment from owned, external-view, and self-view inputs. The complete
  `umol-py` suites pass with 635 Rust tests and 387 Python tests; workspace clippy,
  rustfmt, and `git diff --check` pass. **Additive (green).** `[dep: S1c]`

### S2 — Old/new field-change vocabulary

- **S2a — DONE — ordinary entity field changes** (`umol-py/src/delta.rs`, `lib.rs`,
  `python/umol/__init__.py`): add six separately named native complex-enum
  mirrors. Every variant is a struct variant with read-only, named `old` and
  `new` fields of the same bound payload type, so construction supports keywords,
  attribute reads use `.old`/`.new`, and positional or named class patterns work.
  These types do not themselves implement lattice semantics and therefore keep
  the Rust names without an `Ast` suffix. Use ordinary `#[pyclass]`, not
  `#[pyclass(frozen)]`; PyO3's generated complex-enum fields are already
  getter-only. All six mirrors are value-equal, unhashable, have an exact named
  repr, and expose non-mutating `inverse() -> Self` returning the concrete
  variant subtype through `into_py_variant`.

  1. **DONE — S2a.1 — shared field-change machinery and `AtomFieldChange`** — create
     `delta.rs`, wire the private module into `lib.rs`, and add a local macro for
     the common named-variant representation, AST conversion, equality, repr,
     and inverse surface while still emitting a separately named pyclass. Prove
     the complete payload vocabulary with all six atom variants:
     `Element { old: ElementAst, new: ElementAst }`,
     `IsotopeMass { old: IsotopeMassAst, new: IsotopeMassAst }`,
     `Charge { old: ValueAst, new: ValueAst }`,
     `ImplicitHydrogens { old: ValueAst, new: ValueAst }`,
     `LonePairs { old: ValueAst, new: ValueAst }`, and
     `Spin { old: SpinStateAst, new: SpinStateAst }`. Register/export the class.
     Rust tests cover every variant's conversion, named fields/repr, equality,
     `inverse`, and double inverse; Python tests cover keyword construction,
     attribute reads, scalar and structured class patterns, and concrete-subtype
     preservation across `inverse`. **Additive (green).** `[dep: —]`

     **Implemented verification:** six Rust round-trip rows cover the full atom
     variant closure; two equality rows, six named-field/repr rows, and six
     inverse/double-inverse rows exercise the generated common surface. Four
     Python tests cover keyword construction, read-only `old`/`new`, positional
     and named class patterns, exact repr, value equality, and concrete subtype
     preservation. The complete `umol-py` suites pass with 655 Rust tests and
     391 Python tests; workspace clippy and rustfmt pass.
  2. **DONE — S2a.2 — bond field changes** — invoke the established machinery for
     `BondFieldChange::{Order, Charge, Spin}` over `ValueAst`, `ValueAst`, and
     `SpinStateAst`, and for the single
     `DativeBondFieldChange::Order { old: ValueAst, new: ValueAst }` variant.
     Register/export both classes. Add exhaustive Rust conversion/inverse rows
     and Python construction, named-field, repr, and match coverage for both
     enums. **Additive (green).** `[dep: S2a.1]`

     **Implemented verification:** three Rust round-trip rows and three
     inverse/double-inverse rows cover every `BondFieldChange` variant; one row
     of each kind covers `DativeBondFieldChange::Order`. Four Python tests cover
     construction, named fields, exact repr, positional and named class patterns,
     value equality, and concrete inverse subtypes for both enums. The focused
     delta suites pass with 28 Rust tests and eight Python tests; the complete
     `umol-py` suites pass with 663 Rust tests and 395 Python tests. Workspace
     clippy and rustfmt pass.
  3. **DONE — S2a.3 — delocalized field changes** — add the symmetric three-variant
     `AromaticSystemFieldChange` and `MulticenterBondFieldChange` mirrors:
     `Electrons { old: ElectronCountsAst, new: ElectronCountsAst }`,
     `Charge { old: ValueAst, new: ValueAst }`, and
     `Spin { old: SpinStateAst, new: SpinStateAst }`. Register/export both
     classes. Test every Rust round trip/inverse and exercise an `Electrons`
     variant from Python so the multi-value electron-count payload and named
     class-pattern surface are covered. **Additive (green).** `[dep: S2a.1]`

     **Implemented verification:** three Rust round-trip rows and three
     inverse/double-inverse rows cover every variant of each enum, including
     undetermined and concrete electron vectors. Four Python tests cover both
     registered classes, named fields, exact repr, positional and named class
     patterns, electron-vector payloads, and concrete inverse subtypes. The
     focused delta suites pass with 40 Rust tests and 12 Python tests; the
     complete `umol-py` suites pass with 675 Rust tests and 399 Python tests.
     Workspace clippy and rustfmt pass.
  4. **DONE — S2a.4 — noncovalent field change and closure verification** — add and
     register/export `NoncovalentBondFieldChange::Kind {
     old: NoncovalentBondKindAst, new: NoncovalentBondKindAst }`. Test its Rust
     conversion/inverse and Python named-field/match surface, then run a closure
     matrix over all 17 S2a variants to verify exact repr, value equality,
     concrete subtype returns, and `inverse().inverse() == original`. The
     complete Rust and Python suites, clippy, rustfmt, and `git diff --check`
     form the stage gate. **Additive (green).**
     `[dep: S2a.1, S2a.2, S2a.3]`

     **Implemented verification:** one Rust round-trip row and one
     inverse/double-inverse row cover `NoncovalentBondFieldChange::Kind`; its
     Python test covers named fields and class matching. A 17-row Python closure
     matrix covers every S2a variant's exact repr, value inequality after one
     inversion, concrete inverse subtype, and equality after double inversion.
     The focused delta suites pass with 42 Rust tests and 30 Python tests; the
     complete `umol-py` suites pass with 677 Rust tests and 417 Python tests.
     Workspace clippy, rustfmt, and `git diff --check` pass.

  **Critical path:** `S2a.1 → {S2a.2, S2a.3} → S2a.4`. No S2a subitem is
  deferrable: S3a requires the atom/bond closure and S3b requires the dative,
  aromatic, multicenter, and noncovalent closure.
- **S2b — DONE — stereo field changes** (`umol-py/src/delta.rs`, `lib.rs`,
  `python/umol/__init__.py`): add the two separately named native complex-enum
  mirrors. Each has the single struct variant
  `Configuration { old: StereoConfigurationAst, new: StereoConfigurationAst }`
  and follows the completed S2a contract exactly: keyword construction,
  read-only `.old`/`.new`, positional and named class patterns, value equality,
  unhashability, exact named repr, and non-mutating `inverse() -> Self` returning
  the concrete variant subtype. Reuse the S2a field-change macro rather than
  introducing stereo-specific representation machinery. These field-change
  types do not implement lattice semantics, so they retain their Rust names
  without an `Ast` suffix and remain ordinary, non-frozen pyclasses.

  1. **DONE — S2b.1 — stereo configuration adapter and atom field change** — implement
     the local S2a macro's payload adapter for `StereoConfigurationAst`, using
     its existing Python-aware `from_ast`/`to_ast` conversion, then invoke the
     macro for `StereoAtomFieldChange::Configuration` and register/export the
     class. Rust tests cover round-trip conversion, named fields/repr, equality,
     inverse, and double inverse. Include both top-level `Undetermined` and
     entity-appropriate `Kinded(Tetrahedral, ...)` payloads, with both
     `StereoCosetAst::Undetermined` and a concrete coset, so the two kinds of
     uncertainty remain distinct. Python tests cover keyword construction,
     read-only fields, exact repr, positional and named class patterns, value
     equality, and concrete-subtype preservation across inverse. **Additive
     (green).** `[dep: S2a.1]`

     **Implemented verification:** two Rust round-trip rows, two equality rows,
     two named-field/repr rows, and two inverse/double-inverse rows distinguish
     unknown geometry from a tetrahedral configuration with an open or concrete
     coset. Four Python tests cover keyword construction, read-only fields,
     unhashability, exact repr, positional and named class patterns, value
     equality, and concrete inverse subtype. The focused delta suites pass with
     50 Rust tests and 34 Python tests; the complete `umol-py` suites pass with
     685 Rust tests and 421 Python tests. Workspace clippy, rustfmt, and
     `git diff --check` pass.

  2. **DONE — S2b.2 — bond field change and stereo closure verification** — invoke the
     established machinery for `StereoBondFieldChange::Configuration` and
     register/export the class. Give it the same Rust and Python surface coverage
     as S2b.1, using top-level `Undetermined` and entity-appropriate
     `Kinded(CisTrans, ...)` configurations with open and concrete cosets. Add a
     closure matrix spanning both S2b classes and all represented uncertainty
     states; verify exact repr, value equality, concrete inverse subtype, and
     `inverse().inverse() == original`. The complete Rust and Python suites,
     workspace clippy, rustfmt, and `git diff --check` form the stage gate.
     **Additive (green).** `[dep: S2b.1]`

     **Implemented verification:** two Rust round-trip rows, two equality rows,
     two named-field/repr rows, and two inverse/double-inverse rows cover the
     cis-trans bond form with unknown geometry, an open coset, and a concrete
     coset. Four direct Python tests cover the bond class's construction,
     read-only fields, unhashability, exact repr, positional and named class
     patterns, value equality, and concrete inverse subtype. A four-row Python
     closure matrix spans both stereo field-change classes and both uncertainty
     transitions. The focused delta suites pass with 58 Rust tests and 42 Python
     tests; the complete `umol-py` suites pass with 693 Rust tests and 429 Python
     tests. Workspace clippy, rustfmt, and `git diff --check` pass.

  **Critical path:** `S2a.1 → S2b.1 → S2b.2`. No S2b subitem is
  deferrable: S3c requires both stereo field-change mirrors.

### S3 — Per-family resolved deltas

- **S3a — DONE — atom and bond deltas** (`umol-py/src/delta.rs`, `atom.rs`,
  `bond.rs`, `lib.rs`, `python/umol/__init__.py`): add two separately named
  native complex-enum bindings with the exact Rust variants and named payloads:

  - `AtomDelta::{Add { id: int, ast: AtomAst }, Remove { id: int, ast:
    AtomAst }, ModifyField { id: int, change: AtomFieldChange },
    ModifyConstraint { id: int, old: AtomConstraintAst | None, new:
    AtomConstraintAst | None }}`;
  - `BondDelta::{Add { id: int, atoms: tuple[int, int], ast: BondAst }, Remove {
    id: int, atoms: tuple[int, int], ast: BondAst }, ModifyField { id: int,
    change: BondFieldChange }, ModifyConstraint { id: int, old:
    BondConstraintAst | None, new: BondConstraintAst | None }}`.

  Both bindings follow the S2 contract: keyword construction, read-only named
  fields, positional and named class patterns, structural value equality,
  unhashability, exact named repr, and non-mutating `inverse() -> Self` returning
  the concrete inverse variant subtype. IDs remain bare Python integers. Bond
  participants remain an ordered two-tuple exactly as stored; they occur only on
  `Add`/`Remove`. `AtomAst`/`BondAst` inputs are snapshotted when the delta is
  constructed, but the stored `.ast` child remains a live mutable object within
  that delta. These types do not implement lattice semantics and retain their
  Rust names without an `Ast` suffix.

  1. **DONE — S3a.1 — conversion cleanup and `AtomDelta`** — remove the
     `FieldValueMirror` trait and give every field-change class its own explicit
     `from_ast`/`to_ast` implementation. The field-change macro is limited to the
     actually uniform Python enum declaration, equality, repr, and inverse
     surface. Define `AtomDelta` and its conversion directly, with no generic
     entity-delta macro or participant-conversion trait. Promote
     `AtomAst::from_inner` from its test-only gate. Its concrete Python-field
     storage clones a constructor `AtomAst` into a fresh `Py<AtomAst>` and field
     reads return that stored child. Register/export the class.

     Rust tests cover every variant's conversion and inverse/double-inverse,
     equality, and named fields/repr; constraint rows include `None` on either
     side and two present values. Python tests cover keyword construction,
     read-only fields, positional and named patterns, exact repr, unhashability,
     concrete inverse subtypes (`Add` ↔ `Remove`, both modify variants retain
     their subtype), and both halves of the ownership contract: later mutation
     of the constructor source does not affect the delta, while mutation through
     `delta.ast` does. **Additive (green).** `[dep: S2a.1]`

  2. **DONE — S3a.2 — `BondDelta` and fixed-pair participants** — promote
     `BondAst::from_inner` from its test-only gate and define `BondDelta` and its
     `from_ast`/`to_ast` conversion directly. `Add`/`Remove` map Rust
     `[AtomId; 2]` to Python `tuple[int, int]` inline, without sorting or
     set-conversion; the modify variants have no participant field. Register and
     export the class. Give it the same Rust conversion/equality/repr/inverse and
     Python construction/field/match/ownership coverage as `AtomDelta`, including
     both optional-constraint directions. Assert that participant order survives
     conversion and inverse and that Python observes a tuple. **Additive
     (green).** `[dep: S2a.2, S3a.1]`

  3. **DONE — S3a.3 — atom/bond closure verification** — add an eight-row Python matrix
     spanning every variant of both classes. Verify exact repr, structural value
     equality, concrete inverse subtype, and `inverse().inverse() == original`;
     retain targeted rows for optional constraints, payload snapshot isolation,
     live stored payloads, and ordered bond participants. Run the complete Rust
     and Python suites, workspace clippy, rustfmt, and `git diff --check` as the
     stage gate. **Additive (green).** `[dep: S3a.1, S3a.2]`

     **Implemented verification:** 36 Rust rows cover conversion, equality,
     exact repr, inverse, and double inverse across all variants and optional
     constraint directions. Eight direct Python tests cover construction,
     read-only fields, pattern matching, source snapshot isolation, live stored
     entity values, ordered tuple participants, and inverse behavior. The
     eight-row closure matrix covers every atom/bond delta variant. The focused
     delta suites pass with 94 Rust tests and 58 Python tests; the complete
     `umol-py` suites pass with 729 Rust tests and 445 Python tests. Workspace
     clippy, rustfmt, and `git diff --check` pass.

  **Critical path:** `S2a.1 → S3a.1 → S3a.2 → S3a.3`, with `S2a.2` also feeding
  S3a.2. No S3a subitem is deferrable: S4a requires both resolved delta bindings.
- **S3b — non-stereo overlay deltas** (`umol-py/src/delta.rs`, `dative.rs`,
  `aromatic.rs`, `multicenter.rs`, `noncovalent.rs`, `lib.rs`,
  `python/umol/__init__.py`): add four separately named native complex-enum
  bindings with the exact Rust variants and named payloads:

  - `DativeBondDelta::{Add { id: int, donors: list[int], acceptor: int, ast:
    DativeBondAst }, Remove { id: int, donors: list[int], acceptor: int, ast:
    DativeBondAst }, ModifyField { id: int, change: DativeBondFieldChange },
    ModifyConstraint { id: int, old: DativeBondConstraintAst | None, new:
    DativeBondConstraintAst | None }}`;
  - `AromaticSystemDelta::{Add { id: int, atoms: list[int], ast:
    AromaticSystemAst }, Remove { id: int, atoms: list[int], ast:
    AromaticSystemAst }, ModifyField { id: int, change:
    AromaticSystemFieldChange }, ModifyConstraint { id: int, old:
    AromaticSystemConstraintAst | None, new: AromaticSystemConstraintAst |
    None }}`;
  - `MulticenterBondDelta::{Add { id: int, atoms: list[int], ast:
    MulticenterBondAst }, Remove { id: int, atoms: list[int], ast:
    MulticenterBondAst }, ModifyField { id: int, change:
    MulticenterBondFieldChange }, ModifyConstraint { id: int, old:
    MulticenterBondConstraintAst | None, new: MulticenterBondConstraintAst |
    None }}`;
  - `NoncovalentBondDelta::{Add { id: int, atoms: tuple[int, int], ast:
    NoncovalentBondAst }, Remove { id: int, atoms: tuple[int, int], ast:
    NoncovalentBondAst }, ModifyField { id: int, change:
    NoncovalentBondFieldChange }, ModifyConstraint { id: int, old:
    NoncovalentBondConstraintAst | None, new: NoncovalentBondConstraintAst |
    None }}`.

  All four bindings follow the completed S3a contract: keyword construction,
  read-only named fields, positional and named class patterns, structural value
  equality, unhashability, exact named repr, and non-mutating `inverse() ->
  Self` returning the concrete inverse variant subtype. IDs and the dative
  acceptor remain bare Python integers. Rust `Vec<AtomId>` participants become
  ordered Python lists without sorting or deduplication; the noncovalent
  `[AtomId; 2]` remains an ordered Python two-tuple. Participant payloads occur
  only on `Add`/`Remove`. Entity-AST constructor inputs are snapshotted into a
  fresh stored child, while the stored `.ast` remains live and mutable through
  the delta. Implement each enum and its `from_rust`/`to_rust` conversion
  explicitly; reuse the existing repr and PyO3 variant helpers, but introduce no
  generic entity-delta macro or conversion trait. These types do not implement
  lattice semantics and retain their Rust names without an `Ast` suffix.

  1. **DONE — S3b.1 — `DativeBondDelta` and directed participants** — add the private
     snapshotting `DativeBondDeltaAstValue` field wrapper and define all four
     variants directly. Convert `DativeBondId`/`AtomId` to bare integers,
     preserve donor order and multiplicity in `list[int]`, and keep `acceptor`
     separate. Wire `DativeBondFieldChange` and optional
     `DativeBondConstraintAst` payloads through their explicit conversions.
     Register/export the class. Rust tests cover every variant's conversion,
     equality, exact repr, inverse, and double inverse, including `None` on each
     constraint side. Python tests cover keyword construction, read-only fields,
     positional and named matching, source snapshot isolation, live stored AST
     mutation, donor order/multiplicity, list representation, and concrete
     inverse subtypes. **Additive (green).** `[dep: S2a.2, S3a.3]`

     **Implemented verification:** six Rust round-trip rows and six
     inverse/double-inverse rows cover every variant and all optional-constraint
     directions; two equality rows preserve donor order, and four repr rows
     cover the complete variant surface. Four Python tests cover keyword
     construction, read-only fields, positional and named matching, source
     snapshot isolation, live stored AST mutation, donor order and duplicates,
     list representation, unhashability, and concrete inverse subtypes. The
     focused delta suites pass with 112 Rust tests and 62 Python tests; the
     complete `umol-py` suites pass with 747 Rust tests and 449 Python tests.
     `umol-py` clippy, rustfmt, and `git diff --check` pass.

  2. **DONE — S3b.2 — `AromaticSystemDelta` and ordered member atoms** — add the private
     snapshotting `AromaticSystemDeltaAstValue` field wrapper and define the
     four variants directly. Preserve the variable-length member-atom vector as
     an ordered Python list and wire `AromaticSystemFieldChange` plus optional
     `AromaticSystemConstraintAst` payloads. Register/export the class. Give it
     the complete Rust conversion/equality/repr/inverse matrix and the same
     Python construction, matching, ownership, optional-constraint, participant
     order/multiplicity, and inverse coverage as S3b.1. **Additive (green).**
     `[dep: S2a.3, S3a.3]`

     **Implemented verification:** six Rust round-trip rows and six
     inverse/double-inverse rows cover every variant and all optional-constraint
     directions; two equality rows preserve member-atom order, and four repr
     rows cover the complete variant surface. Four Python tests cover keyword
     construction, read-only fields, positional and named matching, source
     snapshot isolation, live stored AST mutation, member order and duplicates,
     list representation, unhashability, and concrete inverse subtypes. The
     focused delta suites pass with 130 Rust tests and 66 Python tests; the
     complete `umol-py` suites pass with 765 Rust tests and 453 Python tests.
     `umol-py` clippy, rustfmt, and `git diff --check` pass.

  3. **DONE — S3b.3 — `MulticenterBondDelta` and ordered member atoms** — add the
     private snapshotting `MulticenterBondDeltaAstValue` field wrapper and
     define the four variants directly. Preserve the variable-length
     member-atom vector as an ordered Python list and wire
     `MulticenterBondFieldChange` plus optional
     `MulticenterBondConstraintAst` payloads. Register/export the class. Give it
     the complete Rust conversion/equality/repr/inverse matrix and the same
     Python construction, matching, ownership, optional-constraint, participant
     order/multiplicity, and inverse coverage as S3b.2. **Additive (green).**
     `[dep: S2a.3, S3a.3]`

     **Implemented verification:** six Rust round-trip rows and six
     inverse/double-inverse rows cover every variant and all optional-constraint
     directions; two equality rows preserve member-atom order, and four repr
     rows cover the complete variant surface. Four Python tests cover keyword
     construction, read-only fields, positional and named matching, source
     snapshot isolation, live stored AST mutation, member order and duplicates,
     list representation, unhashability, and concrete inverse subtypes. The
     focused delta suites pass with 148 Rust tests and 70 Python tests; the
     complete `umol-py` suites pass with 783 Rust tests and 457 Python tests.
     `umol-py` clippy, rustfmt, and `git diff --check` pass.

  4. **DONE — S3b.4 — `NoncovalentBondDelta` and fixed-pair participants** — add the
     private snapshotting `NoncovalentBondDeltaAstValue` field wrapper and
     define the four variants directly. Map `[AtomId; 2]` to `tuple[int, int]`
     inline without sorting, and wire `NoncovalentBondFieldChange` plus optional
     `NoncovalentBondConstraintAst` payloads. Register/export the class. Give it
     the complete Rust conversion/equality/repr/inverse matrix and Python
     construction, matching, ownership, optional-constraint, ordered-tuple, and
     inverse coverage. **Additive (green).** `[dep: S2a.4, S3a.3]`

     **Implemented verification:** six Rust round-trip rows and six
     inverse/double-inverse rows cover every variant and all optional-constraint
     directions; two equality rows preserve endpoint order, and four repr rows
     cover the complete variant surface. Four Python tests cover keyword
     construction, read-only fields, positional and named matching, source
     snapshot isolation, live stored AST mutation, endpoint order, tuple
     representation, unhashability, and concrete inverse subtypes. The focused
     delta suites pass with 166 Rust tests and 74 Python tests; the complete
     `umol-py` suites pass with 801 Rust tests and 461 Python tests. `umol-py`
     clippy, rustfmt, and `git diff --check` pass.

  5. **DONE — S3b.5 — non-stereo overlay closure verification** — add a sixteen-row
     Python matrix spanning every variant of all four classes. Verify exact
     repr, structural value equality, concrete inverse subtype, inequality
     after one non-identity inversion, and `inverse().inverse() == original`;
     retain targeted rows for optional constraints, payload snapshot isolation,
     live stored payloads, directed dative participants, ordered variable-length
     participant lists, and the ordered noncovalent pair. Run the complete Rust
     and Python suites, workspace clippy, rustfmt, and `git diff --check` as the
     stage gate. **Additive (green).**
     `[dep: S3b.1, S3b.2, S3b.3, S3b.4]`

     **Implemented verification:** the sixteen-row Python closure matrix covers
     `Add`, `Remove`, `ModifyField`, and `ModifyConstraint` for each of
     `DativeBondDelta`, `AromaticSystemDelta`, `MulticenterBondDelta`, and
     `NoncovalentBondDelta`. Every row verifies exact repr, structural value
     equality through double inversion, the concrete inverse subtype, and
     inequality after one non-identity inversion. The targeted ownership,
     optional-constraint, directed-participant, ordered-list, duplicate-member,
     and ordered-tuple tests remain in place. The focused delta suite passes
     with 90 Python tests; the complete `umol-py` suites pass with 801 Rust
     tests and 477 Python tests. Workspace clippy, rustfmt, and
     `git diff --check` pass.

  **Critical path:** `{S2a.2 → S3b.1, S2a.3 → {S3b.2, S3b.3}, S2a.4 →
  S3b.4} → S3b.5`, with S3a.3 supplying the established entity-delta surface to
  S3b.1–S3b.4. No S3b subitem is deferrable: S4a requires all four bindings.
- **S3c — stereo deltas** (`umol-py/src/delta.rs`, `stereo.rs`, `lib.rs`,
  `python/umol/__init__.py`): add two separately named native complex-enum
  bindings with the exact Rust variants and named payloads:

  - `StereoAtomDelta::{Add { id: int, site: int, ligands:
    list[StereoLigand], ast: StereoAtomAst }, Remove { id: int, site: int,
    ligands: list[StereoLigand], ast: StereoAtomAst }, ModifyField { id: int,
    change: StereoAtomFieldChange }, ModifyConstraint { id: int, kind:
    StereoKind | None, old: StereoAtomConstraintAst | None, new:
    StereoAtomConstraintAst | None }, Apply { id: int, kind: StereoKind,
    permutation: Permutation }, Swap { id: int, kind: StereoKind }, Mirror {
    id: int, kind: StereoKind }}`;
  - `StereoBondDelta::{Add { id: int, site: int, ligands:
    list[StereoLigand], ast: StereoBondAst }, Remove { id: int, site: int,
    ligands: list[StereoLigand], ast: StereoBondAst }, ModifyField { id: int,
    change: StereoBondFieldChange }, ModifyConstraint { id: int, kind:
    StereoKind | None, old: StereoBondConstraintAst | None, new:
    StereoBondConstraintAst | None }, Apply { id: int, kind: StereoKind,
    permutation: Permutation }, Swap { id: int, kind: StereoKind }, Mirror {
    id: int, kind: StereoKind }}`.

  Both bindings follow the completed overlay contract: keyword construction,
  read-only named fields, positional and named class patterns, structural value
  equality, unhashability, exact named repr, and non-mutating `inverse() ->
  Self` returning the concrete inverse variant subtype. IDs and sites remain
  bare Python integers. Ligand frames become ordered Python lists without
  sorting or deduplication and occur only on `Add`/`Remove`. Entity-AST
  constructor inputs are snapshotted into a fresh stored child, while the
  stored `.ast` remains live and mutable through the delta. `ModifyConstraint`
  preserves `kind` independently while exchanging `old` and `new`. `Apply`
  retains its kind and replaces its permutation with the inverse; `Swap` and
  `Mirror` are involutions. Implement each enum and its `from_rust`/`to_rust`
  conversion explicitly, without a shared stereo-delta macro or conversion
  trait. These types do not implement lattice semantics and retain their Rust
  names without an `Ast` suffix.

  1. **DONE — S3c.1 — `StereoAtomDelta` and atom-centered ligand frames** — promote
     the `stereo_value!`-generated `StereoAtomAst` owned-AST constructor to
     production use, add the private snapshotting `StereoAtomDeltaAstValue`
     field wrapper, and define all seven variants directly. Convert
     `StereoAtomId`, the atom site, and ligand-frame members to the existing
     Python scalar/value types while preserving frame order and multiplicity.
     Wire `StereoAtomFieldChange`,
     optional-kind `StereoAtomConstraintAst` changes, and `Permutation` through
     their existing explicit conversions. Register and export the class.

     Rust tests cover every variant's round trip, equality, exact repr, inverse,
     and double inverse; constraint rows cover `kind=None` and `Some`, `None` on
     either constraint side, and two present constraints. Use a non-involutive
     permutation for `Apply` so the inverse-image assertion is meaningful, and
     verify `Swap` and `Mirror` are self-inverse. Python tests cover keyword
     construction, read-only fields, positional and named matching, source
     snapshot isolation, live stored AST mutation, ordered duplicate-preserving
     ligand lists, optional kind/constraints, permutation degree and image, and
     concrete inverse subtypes. **Additive (green).** `[dep: S2b.2]`

     **Implemented verification:** nine Rust round-trip rows and nine
     inverse/double-inverse rows cover all seven variants, both optional-kind
     states, every optional-constraint direction, non-involutive permutation
     inversion, and the `Swap`/`Mirror` involutions. Three equality rows cover
     ordered ligand frames and distinct permutations; seven repr rows cover the
     complete variant surface. Seven Python tests cover keyword construction,
     read-only fields, positional and named matching, source snapshot isolation,
     live stored AST mutation, ordered duplicate-preserving ligand lists,
     optional kind independent of optional constraints, permutation degree and
     image, unhashability, and concrete inverse subtypes. The focused delta
     suites pass with 194 Rust tests and 97 Python tests; the complete `umol-py`
     suites pass with 829 Rust tests and 484 Python tests. `umol-py` clippy,
     rustfmt, and `git diff --check` pass.

  2. **DONE — S3c.2 — `StereoBondDelta` and bond-centered ligand frames** — promote
     the generated `StereoBondAst` owned-AST constructor to production use, add
     the private snapshotting `StereoBondDeltaAstValue` field wrapper, and
     define all seven variants directly. Map the bond site to a bare integer, preserve the
     ligand frame as an ordered duplicate-preserving list, and wire
     `StereoBondFieldChange`, optional-kind `StereoBondConstraintAst` changes,
     and `Permutation` through separate `StereoBondDelta` conversions. Register
     and export the class. Give it the same exhaustive Rust conversion,
     equality, repr, and inverse matrix and the same Python construction,
     matching, ownership, optional-kind/constraint, ligand-frame, permutation,
     and inverse coverage as S3c.1, using bond-appropriate stereo kinds and
     constraints. **Additive (green).** `[dep: S2b.2, S3c.1]`

     **Implemented verification:** nine Rust round-trip rows and nine
     inverse/double-inverse rows cover all seven variants, both optional-kind
     states, every optional-constraint direction, non-involutive permutation
     inversion, and both involutions. Three equality rows cover ordered ligand
     frames and distinct permutations; seven repr rows cover the complete
     variant surface. Seven Python tests cover construction, read-only fields,
     matching, snapshot/live-child ownership, ordered duplicate-preserving
     ligand lists, optional kind and constraints, permutation degree and image,
     unhashability, and concrete inverse subtypes. The focused delta suites pass
     with 222 Rust tests and 104 Python tests; the complete `umol-py` suites pass
     with 857 Rust tests and 491 Python tests. `umol-py` clippy, rustfmt, and
     `git diff --check` pass.

  3. **DONE — S3c.3 — stereo-delta closure verification** — add a fourteen-row Python
     matrix spanning all seven variants of both classes. Verify exact repr,
     structural value equality, concrete inverse subtype, and double inversion
     on every row. Require inequality after one inversion for `Add`, `Remove`,
     non-identity `ModifyField`/`ModifyConstraint`, and non-involutive `Apply`;
     require equality after one inversion for the self-inverse `Swap` and
     `Mirror` rows. Retain targeted rows for source snapshot isolation, live
     stored payloads, optional kind independent of optional constraints,
     ordered duplicate-preserving ligand frames, and permutation image/degree.
     Run the complete Rust and Python suites, workspace clippy, rustfmt, and
     `git diff --check` as the stage gate. **Additive (green).**
     `[dep: S3c.1, S3c.2]`

     **Implemented verification:** the fourteen-row Python closure matrix covers
     every variant of `StereoAtomDelta` and `StereoBondDelta`. Every row verifies
     exact repr, structural equality through double inversion, and the concrete
     inverse subtype. The ten transforming rows verify inequality after one
     inversion; the four `Swap`/`Mirror` rows verify equality after one inversion.
     The focused delta suite passes with 118 Python tests; the complete
     `umol-py` suites pass with 857 Rust tests and 505 Python tests. Workspace
     clippy, rustfmt, and `git diff --check` pass.

  **Critical path:** `S2b.2 → S3c.1 → S3c.2 → S3c.3`. No S3c subitem is
  deferrable: S4a requires both complete seven-variant bindings.
- **S3d — DONE — molecule constraint delta** (`umol-py/src/delta.rs`, `lib.rs`,
  `python/umol/__init__.py`): add `ConstraintDelta::{Add, Remove}` over the
  recursive S1 `Constraint` binding. Expose both variants as named-field
  classes, `Add { constraint: Constraint }` and
  `Remove { constraint: Constraint }`, with keyword construction, a read-only
  `.constraint`, positional and named class patterns, structural value equality,
  unhashability, exact named repr, and non-mutating `inverse() -> Self` returning
  the concrete opposite variant subtype. Snapshot the constructor's recursive
  constraint tree into a fresh stored `Constraint`; subsequent access to
  `.constraint` returns that stored child, so changes made through its nested
  bound payloads remain visible to conversion and inversion. Implement
  `from_rust` and `to_rust` directly for `ConstraintDelta`, using a private field
  wrapper only for snapshotting and exposing the stored child. The type has no
  lattice semantics and therefore retains the Rust name without an `Ast` suffix.

  1. **DONE — S3d.1 — `ConstraintDelta` binding and recursive payload ownership** — add
     the private snapshotting field wrapper and the two-variant public enum;
     implement exact repr, equality, inverse, and explicit Rust conversion;
     register and export the class. Rust tests cover both conversion directions,
     equality, exact repr, inverse, and double inverse for `Add` and `Remove`,
     using both an entity leaf and a recursive Boolean constraint tree. Python
     tests cover keyword construction, read-only fields, positional and named
     matching, source snapshot isolation, access to the stored child, live
     nested-payload changes, unhashability, and concrete inverse subtypes.
     **Additive (green).** `[dep: S1d]`

     **Implemented verification:** two Rust round-trip rows, three equality rows,
     two exact-repr rows, and two inverse/double-inverse rows cover both variants,
     an entity leaf, and a recursive Boolean tree. Four Python tests cover public
     import, keyword construction, read-only fields, positional and named
     matching, constructor snapshot isolation, stable access to the stored child,
     live mutation through a nested subpattern molecule, unhashability, and
     concrete inverse subtypes. The focused suites pass with nine Rust cases and
     four Python tests; the complete `umol-py` suites pass with 866 Rust tests and
     509 Python tests. `umol-py` clippy, rustfmt, and `git diff --check` pass.

  2. **DONE — S3d.2 — constraint-delta closure verification** — add a two-row Python
     matrix spanning `Add` and `Remove`, with distinct leaf and recursive
     constraints. Every row verifies exact repr, structural equality, the
     concrete inverse subtype, inequality after one inversion, and equality
     after double inversion. Retain the targeted ownership tests from S3d.1 and
     run the complete Rust and Python suites, workspace clippy, rustfmt, and
     `git diff --check` as the stage gate. **Additive (green).** `[dep: S3d.1]`

     **Implemented verification:** the two-row Python closure matrix covers
     `Add` with an entity leaf and `Remove` with a recursive Boolean constraint.
     Both rows verify exact repr, structural inequality after one inversion, the
     concrete opposite variant subtype, and equality after double inversion.
     The focused constraint-delta suite passes with six Python tests, including
     the four targeted S3d.1 ownership and matching tests. The complete
     `umol-py` suites pass with 866 Rust tests and 511 Python tests. Workspace
     clippy, rustfmt, and `git diff --check` pass.

  **Critical path:** `S1d → S3d.1 → S3d.2`. Neither S3d subitem is
  deferrable: S4a requires the complete `ConstraintDelta` binding.

### S4 — `Delta` and `Deltas`

- **S4a — DONE — top-level `Delta` sum** (`umol-py/src/delta.rs`, `lib.rs`,
  `python/umol/__init__.py`): add the tuple-shaped
  `Delta::{Atom, Bond, DativeBond, AromaticSystem, MulticenterBond,
  NoncovalentBond, StereoAtom, StereoBond, Constraint}` binding over the nine
  completed family-delta classes. Each variant takes exactly one positional
  child in Rust declaration order, exposes the corresponding read-only `_0`
  field, and supports positional class patterns. The supplied family-delta
  object becomes the stored child rather than being copied again, so nested
  entity AST changes remain live through the top-level sum; the family-delta
  constructors remain the snapshot boundary for their AST inputs, and S4b's
  `Deltas` container will snapshot whole `Delta` entries. Give `Delta`
  structural value equality, unhashability, exact positional repr, and a
  non-mutating `inverse() -> Self`. Implement the nine-arm `from_rust` and
  `to_rust` dispatch directly, without a conversion trait or dispatch macro.
  `Delta` has no lattice semantics and retains the Rust name without an `Ast`
  suffix.

  1. **DONE — S4a.1 — `Delta` binding and nine-family conversion dispatch** — define
     the complete enum, explicit conversions, equality, repr, and inverse;
     register and export it. Rust tests cover round trips and inverse/double
     inverse for all nine variants, plus representative equality and exact-repr
     rows. Targeted Python tests cover public import, construction, the read-only
     direct field, positional matching, child identity, live nested entity-AST
     mutation, structural equality, unhashability, and a detached inverse whose
     outer family variant is retained while its child is inverted. **Additive
     (green).** `[dep: S3a.3, S3b.5, S3c.3, S3d.2]`

     **Implemented verification:** nine Rust round-trip rows and nine
     inverse/double-inverse rows cover every outer family; three equality rows
     distinguish equal values, outer variants, and nested children; three exact-
     repr rows cover ordinary entity, stereo, and recursive-constraint payloads.
     Two Python tests cover public import, positional construction and matching,
     read-only direct fields, retained child identity, live nested AST mutation,
     structural equality, unhashability, and detached nested inversion while the
     outer family remains unchanged. The focused suites pass with 24 Rust cases
     and two Python tests; the complete `umol-py` suites pass with 890 Rust tests
     and 513 Python tests. `umol-py` clippy, rustfmt, and `git diff --check` pass.

  2. **DONE — S4a.2 — top-level delta closure verification** — add a nine-row Python
     matrix with one representative `Add` from each family. Every row verifies
     exact repr, structural equality, the retained concrete outer variant,
     the concrete nested `Remove` inverse subtype, inequality after one
     inversion, and equality after double inversion. Add one exhaustive class-
     pattern traversal whose nine arms descend through `Delta` into the expected
     family-delta type and payload, proving that every dispatch arm is observable
     from Python. Retain the S4a.1 ownership tests and run the complete Rust and
     Python suites, workspace clippy, rustfmt, and `git diff --check` as the stage
     gate. **Additive (green).** `[dep: S4a.1]`

     **Implemented verification:** the nine-row Python closure matrix covers one
     `Add` from every outer family. Every row verifies exact repr, the concrete
     outer variant before and after inversion, the concrete nested `Remove`
     subtype, inequality after one inversion, and structural equality after
     double inversion. A separate exhaustive traversal class-matches through all
     nine `Delta` variants into the corresponding family `Add` and verifies its
     structural payload. The focused top-level delta suite passes with 12 Python
     tests, including the two targeted S4a.1 ownership tests. The complete
     `umol-py` suites pass with 890 Rust tests and 523 Python tests. Workspace
     clippy, rustfmt, and `git diff --check` pass.

  **Critical path:** `S3a.3, S3b.5, S3c.3, S3d.2 → S4a.1 → S4a.2`. Neither
  S4a subitem is deferrable: S4b requires the complete nine-family sum binding.
- **S4b — DONE — Python `Deltas` container** (`umol-py/src/delta.rs`, `error.rs`,
  `lib.rs`, `python/umol/__init__.py`): bind the Rust `Deltas` value with the same
  ownership discipline as `Constraints`. `Deltas(entries=())` snapshots an
  iterable of `Delta` values into an owned Rust container, preserving insertion
  order and duplicates. `append` and `extend` snapshot their complete inputs
  before mutating the target; `deltas.extend(deltas)` therefore appends exactly
  one copy of the original sequence without a PyO3 double-borrow. Integer
  indexing, negative indexing, and iteration return detached `Delta` snapshots,
  and `__len__` supplies ordinary Python truthiness. Give the container
  structural value equality, unhashability, and exact constructor-style repr.
  Unlike `Constraints`, it has no live view, `clear`, item assignment/removal, or
  general `update`: the Rust API exposes append-only input construction before
  normalization. `canonicalize() -> Deltas` clones the stored Rust value before
  invoking consuming Rust canonicalization, returns a fresh container, leaves
  the source unchanged, and maps `Contradiction` to the dedicated public
  `ContradictionError`.

  1. **DONE — S4b.1 — `ContradictionError` boundary** (`umol-py/src/error.rs`, `lib.rs`,
     `python/umol/__init__.py`) — add the public exception and a focused mapping
     function from `umol_ast::ast::Contradiction`, preserving the Rust message
     `"reached a contradiction"`. Register/export it independently of the
     container. Rust and Python tests verify the exact exception class, message,
     and public import. **Additive (green).** `[dep: none]`

     **Implemented verification:** one Rust test verifies that the shared mapper
     produces the exact `ContradictionError` type and the message
     `"reached a contradiction"`. One Python test verifies the public import,
     identity with `umol.ContradictionError`, exception inheritance, and exact
     message. The complete `umol-py` suites pass with 891 Rust tests and 524
     Python tests. `umol-py` clippy, rustfmt, and `git diff --check` pass.

  2. **DONE — S4b.2 — owned sequence, conversion, and detached reads**
     (`umol-py/src/delta.rs`, `lib.rs`, `python/umol/__init__.py`) — add the
     Rust-backed `Deltas` pyclass, private snapshot `DeltaIter`, negative-index
     resolver, explicit Rust conversion helpers, default-empty iterable
     constructor, exact repr, equality, `__len__`, `__getitem__`, and `__iter__`;
     register/export the public container while keeping the iterator internal.
     Rust tests cover empty and populated conversion, repr/equality, positive and
     negative indexing, range errors, and detached iteration. Python tests cover
     default and iterable construction, order and duplicate preservation,
     truthiness, unhashability, exact repr, range errors, and prove that mutating
     a value returned by indexing or iteration does not change the container.
     **Additive (green).** `[dep: S4a.2]`

     **Implemented verification:** seven Rust index-resolution rows cover valid
     positive/negative positions and exact range errors; 15 container rows cover
     empty and populated construction/conversion, equality, exact repr, length,
     positive/negative indexing, range errors, and detached iteration. Three
     Python tests cover the default-empty and iterable constructors, insertion
     order and duplicates, truthiness, equality, unhashability, exact repr,
     constructor snapshot isolation, positive/negative indexing, exact range
     errors, and detached indexed/iterated values. The complete `umol-py` suites
     pass with 913 Rust tests and 527 Python tests. `umol-py` clippy, rustfmt, and
     `git diff --check` pass.

  3. **DONE — S4b.3 — append and RHS-first extend** (`umol-py/src/delta.rs`) — add
     `append(delta)` and an entity-container-style resolved extend input that
     accepts another `Deltas` or any iterable of `Delta`. Convert the complete
     RHS to owned Rust deltas before borrowing the target for mutation; append in
     source order without sorting or deduplication. Rust and Python tests cover
     single append, container and iterable extension, source snapshot isolation,
     order and duplicate preservation, and self-extension producing exactly two
     copies of the original sequence without aliasing or a double-borrow panic.
     **Additive (green).** `[dep: S4b.2]`

     **Implemented verification:** four Rust tests cover append, extension from
     another container, extension from entries, and self-extension. Three Python
     tests verify return values, insertion order, duplicates, snapshot isolation
     from both container and iterable sources, and exact doubling under
     `deltas.extend(deltas)` without aliasing or a double-borrow panic. The
     focused `Deltas` suites pass with 19 Rust cases and six Python tests; the
     complete `umol-py` suites pass with 917 Rust tests and 530 Python tests.
     `umol-py` clippy, rustfmt, and `git diff --check` pass.

  4. **DONE — S4b.4 — non-mutating canonicalization and container closure**
     (`umol-py/src/delta.rs`) — expose `canonicalize() -> Deltas` by cloning the
     held Rust collection, invoking `Canonicalize`, and mapping failure through
     S4b.1. Test field-change fusion, add/remove cancellation, canonical ordering
     across families, idempotence, and source non-mutation; test a discontinuous
     old/new chain raises the exact `ContradictionError` without changing the
     source. Add a Python closure matrix covering construction, detached reads,
     append, extend, and canonicalize over representative ordinary, stereo, and
     constraint deltas. Run the complete Rust and Python suites, workspace
     clippy, rustfmt, and `git diff --check` as the stage gate. **Additive
     (green).** `[dep: S4b.1, S4b.3]`

     **Implemented verification:** three Rust normalization rows cover field-change
     fusion, add/remove cancellation, and canonical ordering across entity
     families; one Rust error test verifies the exact `ContradictionError` class
     and message for a discontinuous field-change chain. Every successful row
     verifies idempotence and source non-mutation, and the error row also verifies
     source non-mutation. Five Python normalization rows cover the same ordinary
     delta semantics plus stereo involution cancellation and constraint
     multiplicity, and verify a fresh result, exact normalized contents,
     idempotence, and source non-mutation. The retained construction, detached
     read, append, and extend tests make the focused closure suite 12 Python
     cases. The complete `umol-py` suites pass with 921 Rust tests and 536 Python
     tests. Workspace clippy, rustfmt, and `git diff --check` pass.

  **Critical path:** `S4a.2 → S4b.2 → S4b.3 → S4b.4`, with S4b.1 joining at
  S4b.4. No S4b subitem is deferrable: S5a requires the complete owned container
  and its canonicalization error contract.

### S5 — `ReactionAst` owned component facade

- **S5a — facade and conversions** (`umol-py/src/reaction.rs`): add
  `ReactionAst { lhs: Py<MoleculeAst>, deltas: Py<Deltas> }`, snapshotting
  constructor, live component getters, snapshotting setters, structural equality,
  `to_ast`/`from_ast`, and repr scaffolding. Test constructor non-aliasing, live
  nested writes through both components, whole-component self-assignment, and
  Rust round trips. **Additive (green).** `[dep: S1d, S4b]`
- **S5b — reaction DSL shortcut** (`reaction.rs`): add `parse`, `str`, and
  constructor-style `repr` through Rust `FromStr`/`Display`; reuse `ParseError`.
  Test representative add/remove/modify/stereo/molecule-constraint reactions and
  canonical render/parse stability. **Additive (green).** `[dep: S5a]`
- **S5c — canonicalize and reverse** (`reaction.rs`): return fresh component
  facades and map `Contradiction`; verify the source is unchanged, normal forms
  are idempotent, reverse twice round-trips, and live components on results remain
  writable. **Additive (green).** `[dep: S5a]`

### S6 — Side construction and composition

- **S6a — `from_sides`** (`reaction.rs`): convert integer atom pairs to
  `Correspondence<NodeId>` using the two molecule sizes, reject duplicate or
  out-of-range pairs with a Python argument error before entering the AST, and
  wrap `ReactionAst::from_sides`. Test partial maps, additions/removals, overlay
  induction, and invalid pairs. **Additive (green).** `[dep: S5a]`
- **S6b — composition scope and compose** (`reaction.rs`): mirror the simple
  `CompositionScope` enum and expose `compose`, defaulting to `RcAnchored` at the
  Python signature while passing the explicit Rust argument. Test empty/no-match,
  one admissible composite, full-vs-anchored result counts, and source
  non-mutation. **Additive (green).** `[dep: S5c]`
- **S6c — reaction-data registration and contract** (`lib.rs`,
  `python/umol/__init__.py`, `tests/test_reaction.py`): export the reaction data
  vocabulary and exercise construct → normalize → render/parse → reverse →
  compose from Python. **Additive (green).** `[dep: S6a, S6b]`

### S7 — Reaction application (required)

- **S7a — correspondence and algorithm mirrors**
  (`umol-py/src/correspondence.rs`): bind read-only per-family correspondence
  views (mated pairs and left/right exposed ids), `MoleculeCorrespondence`, and
  `SubgraphIsomorphismAlgorithm` including `ArcMatch(path_length)`. Add the direct
  `umol-graph-core` dependency required at the binding boundary. **Additive
  (green).** `[dep: S0b]`
- **S7b — Python `ReactionDerivation`** (`umol-py/src/reaction.rs`): wrap the
  fully owned Rust value from S0 and expose lhs/rhs/comap/atom-map, `reverse`,
  `chain`, and `to_reaction`. The derivation is an immutable result; molecule-side
  getters return explicit owned snapshots so Python cannot mutate a side and
  silently invalidate its correspondence. Test independence from the input host
  and every operation. **Additive (green).** `[dep: S0b, S5a, S7a]`
- **S7c — `apply_at` and `apply`** (`reaction.rs`, `error.rs`): expose exact-match
  application over `MoleculeCorrespondence` and all-match application over the
  algorithm mirror, returning owned derivations; map every `ApplyError` variant to
  a typed Python exception with structured fields where present. Cover dangling,
  inconsistent, structural-conflict, transaction-failure, zero-match, and
  multi-match cases. **Additive (green).** `[dep: S7b]`
- **S7d — full public registration and end-to-end application contract**
  (`lib.rs`, `python/umol/__init__.py`, `tests/test_reaction.py`): export the
  derivation, correspondence, algorithm, and application-error types; exercise
  parse/construct → apply → inspect atom map → reverse/chain → recover reaction
  from Python. **Additive (green).** `[dep: S6c, S7c]`

At S7 the required deliverable is complete. Application is not deferrable.

## Critical path and deferrals

The reaction-data path is:

`S1a → S1b → S1c → S1d → S3d → S4a → S4b → S5a → S5c → S6b → S6c`

The ownership/application path is:

`S0a → S0b → S7a → S7b → S7c → S7d`

The paths join at S7b through its dependency on S5a and at the final S7d
end-to-end contract through S6c. `S2a/S2b → S3a/S3b/S3c` proceeds alongside the
constraint path and joins at S4a. `S5b` and `S6a` proceed alongside
`S5c → S6b` once S5a lands.

Only `ReactionSpanAst`/`EntitySpan<T>` and reaction metadata/alias preservation
are deferrable. `ReactionDerivation` remains part of the required owned Python
surface.

S0 is the sole breaking stage and must restore a green workspace before S1.
S1–S7 are additive, and every stage must end green.
