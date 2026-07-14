# 150 · Python bindings for `Deltas` and `ReactionAst` (plan)

Status: **ACTIVE — S0 complete; S1a is next**
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
- **S1b — recursive `Constraint` mirror** (`constraint/molecule.rs`): add one
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
  3. **Aggregate leaves (2)** — reuse the S1a mirrors directly:
     - `Relational(RelationalConstraint)` ↔
       `Constraint::Relational(ast::RelationalConstraint)`;
     - `Molecule(MoleculeConstraint)` ↔
       `Constraint::Molecule(ast::MoleculeConstraint)`.
  4. **Recursive combinators (3)** — recursive children are native variant
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

  **Implemented verification:** one round-trip case for each of the six ordinary
  and two stereo entity leaves, including distinct stereo kinds and their bound
  structural constraint payloads. The focused S1a/S1b suite has 47 cases.

  **Tests:** one round-trip row for every variant (13), with distinct stereo
  kinds and representative child constraints; explicit empty `And([])` and
  `Or([])` cases; and one deep tree combining entity, relational, and molecule
  leaves beneath `And`/`Or`/`Not`. Exercise native variant construction, field
  access, repr, equality, and recursive pattern matching from Rust-side PyO3;
  import-level Python coverage remains in S1d when the type is registered.
  **Additive (green).** `[dep: S1a]`
- **S1c — molecule constraint container** (`constraint/molecule.rs`): add
  `Constraints`, `ConstraintsView`, and their iterators with sequence operations
  and safe RHS-first update helpers. Test owned mutation, live molecule
  write-through, negative indexing, and self-aliasing. **Additive (green).**
  `[dep: S1b]`
- **S1d — `MoleculeAst` wiring** (`umol-py/src/molecule.rs`, `lib.rs`,
  `python/umol/__init__.py`): expose `mol.constraints`, add the keyword-only
  `constraints=()` input to `from_parts`, and register/export the S1 classes.
  Existing calls remain valid. **Additive (green).** `[dep: S1c]`

### S2 — Old/new field-change vocabulary

- **S2a — ordinary entity field changes** (`umol-py/src/delta.rs`): mirror the
  atom, bond, dative, aromatic, multicenter, and noncovalent field-change enums,
  preferably through a local macro that still emits separately named pyclasses.
  Add conversion and `inverse` tests for every variant and Python match tests for
  at least one scalar and one multi-field family. **Additive (green).** `[dep: —]`
- **S2b — stereo field changes** (`delta.rs`): mirror stereo-atom and stereo-bond
  `Configuration { old, new }`, reusing `StereoConfigurationAst`; cover kinded
  and undetermined configurations. **Additive (green).** `[dep: —]`

### S3 — Per-family resolved deltas

- **S3a — atom and bond deltas** (`delta.rs`): add `AtomDelta` and `BondDelta`
  with all four operation variants, participants on bond add/remove, conversion,
  pattern-match, and inverse tests. **Additive (green).** `[dep: S2a]`
- **S3b — non-stereo overlay deltas** (`delta.rs`): add dative, aromatic,
  multicenter, and noncovalent mirrors. Cover each participant shape and each
  family-specific constraint payload. **Additive (green).** `[dep: S2a]`
- **S3c — stereo deltas** (`delta.rs`): add stereo-atom and stereo-bond mirrors,
  including `Apply`/`Swap`/`Mirror`, ordered ligand frames, optional kind on
  constraint changes, permutation inversion, and all inverse cases. A macro may
  generate the symmetric halves. **Additive (green).** `[dep: S2b]`
- **S3d — molecule constraint delta** (`delta.rs`): add
  `ConstraintDelta::{Add, Remove}` over the S1 recursive mirror with conversion,
  match, and inverse tests. **Additive (green).** `[dep: S1b]`

### S4 — `Delta` and `Deltas`

- **S4a — top-level `Delta` sum** (`delta.rs`): add all nine variants,
  conversion dispatch, `inverse`, structural value equality, registration/export,
  and an exhaustive Python match test that traverses one value from every
  family. **Additive (green).** `[dep: S3a, S3b, S3c, S3d]`
- **S4b — Python `Deltas` container** (`umol-py/src/delta.rs`): implement the
  owned sequence surface, `append`/`extend`, and consuming-style `canonicalize`;
  add `ContradictionError` mapping in `error.rs`. Test canonical
  fusion/cancellation, discontinuous-chain failure, negative indexing,
  iteration, and input snapshot semantics. **Additive (green).** `[dep: S4a]`

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
