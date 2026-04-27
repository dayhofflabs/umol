# umol-graph engines / config / data restructure

Date: 2026-04-25

## Context

The trigger was wiring the AST into `umol-graph` after the `ast/` and `dsl/` modules were factored out into `umol-ast`. The legacy `crate::ast::*` and `crate::dsl::*` modules in `umol-graph` overlap with `umol-ast` (same type names, stale shapes) and need to go; the `api` module wraps the legacy types and was disabled in `lib.rs` (not deleted) for the refactor.

## Three groups of objects

**Data**: `MoleculeDsl` (`umol_ast::dsl`) → `MoleculeAst` (`umol_ast::ast`) → `Molecule` (`umol_graph`). `Molecule` = AST + positions + caches + (metadata?). Metadata-on-Molecule is unsettled (round-trip preservation vs. purely semantic object).

**Engines** (in `ops/`): constructed `Engine::new(&config).op(&data) -> Result<Solution<Data>, Error>`. Engines are:
- `Resolver` — composite of `ValenceResolver`, `AromaticityResolver`, future `ChiralityResolver`.
- `Validator` — composite of tier-2 invariant checkers (`ElectronInvariantValidator`, `SpinCouplingValidator`) plus `ConstraintValidator`.
- `Matcher` — for substructure + transformations (kekulization, aromatization). Deferred; file kept on disk, disabled in `ops.rs`.

**Config**: top-level `ChemistryModel` (renamed from `Chemistry`; `ModelChemistry` rejected because the `Model` prefix breaks down for sub-configs). Wraps `ValenceModel` and `AromaticityModel` directly. No `ResolveConfig` wrapper. One shared chemistry model for all ops; split if it gets messy.

## Hard rules

- Engine and config are separate types. No "config is the engine" pattern (rejected `ValenceTheory` shape where the enum carried both the strategy choice and its config).
- No `Molecule::new`, no `Molecule::resolve`. `Molecule` is opaque output of an engine; its only constructor is `pub(crate) Molecule::from_inner`.
- No CSP-flavored "Theory" naming on chemistry types. `Solution` / `Determined` / `Underdetermined` / `Contradictory` stay — they are CSP terms and there is no clean replacement.
- Narrowing is one-shot. No fixpoint iteration in scope; revisit if it ever becomes necessary.
- AST does not validate. Tier-1 (structural) is enforced at `MoleculeAst::new`; tier-2 and constraint↔data agreement are validator concerns.
- Validators are always-on (no `ValidateConfig` on `ChemistryModel`). Turning them off makes downstream resolvers accept garbage.

## Engine signatures

```rust
ValenceResolver::resolve(&self, ast: &mut MoleculeAst)
    -> Result<Solution<(), ValenceContradiction>, ValenceError>;

AromaticityResolver::resolve(&self, ast: &mut MoleculeAst)
    -> Result<Solution<(), AromaticityContradiction>, AromaticityError>;

ElectronInvariantValidator::validate(&self, ast: &MoleculeAst)
    -> Result<Solution<(), ElectronInvariantContradiction>, ElectronInvariantError>;

SpinCouplingValidator::validate(&self, ast: &MoleculeAst)
    -> Result<Solution<(), SpinCouplingContradiction>, SpinCouplingError>;

ConstraintValidator::validate(&self, ast: &MoleculeAst)
    -> Result<Solution<(), ConstraintContradiction>, ConstraintError>;

// Standalone-atom validation (atom-typing registry use, no topology):
Validator::validate_atom(&self, ast: &AtomAst) -> Result<Solution<(), C>, E>;
```

Resolvers mutate in place (`T = ()`, AST is `&mut`). Validators borrow (`&MoleculeAst`). No chaining via an AST trait — would inflate the AST surface. Composition is statement-level:

```rust
valence_resolver.resolve(&mut ast)?;
aromaticity_resolver.resolve(&mut ast)?;
electron_validator.validate(&ast)?;
```

`Solution<T, C>`:

```rust
pub enum Solution<T, C> {
    Determined(T),
    Underdetermined(T),
    Contradictory(C),
}
```

For validators, `Determined` and `Underdetermined` are both successful outcomes. Non-ground inputs yield `Underdetermined` trivially. Ground-checking is a separate concern, not the validator's responsibility.

`Contradictory(C)` carries diagnostic info typed per engine. `Err(EngineError)` is reserved for setup or parameter-table gaps that aren't chemistry contradictions (e.g., `HmoMissingParameters` for an element not in the Van-Catledge table).

Per-engine `Contradiction` and `Error` types:
- `ValenceContradiction` / `ValenceError`
- `AromaticityContradiction` / `AromaticityError`
- `ElectronInvariantContradiction` / `ElectronInvariantError`
- `SpinCouplingContradiction` / `SpinCouplingError`
- `ConstraintContradiction` / `ConstraintError`

Composites union them via `From` impls. Stringly diagnostics rejected.

## ops/ layout

```
ops.rs                module decl; matcher disabled
ops/
  config.rs           ChemistryModel + ValenceModel + AromaticityModel
  solution.rs         Solution<T, C>
  resolver.rs         composite Resolver + ResolverContradiction + ResolverError
  validator.rs        composite Validator + ElectronInvariantValidator
                      + SpinCouplingValidator + ConstraintValidator
                      + their Contradiction/Error types
  valence.rs          ValenceResolver enum + ValenceContradiction + ValenceError
  valence/
    atom_typing.rs
    counts.rs
    registry.rs
    table.rs
  aromaticity.rs      AromaticityResolver enum + AromaticityContradiction + AromaticityError
  aromaticity/
    hueckel_rule.rs
    hmo.rs
    clar.rs
```

```rust
struct Resolver { valence: ValenceResolver, aromaticity: AromaticityResolver }

struct Validator {
    electron_invariant: ElectronInvariantValidator,
    spin_coupling: SpinCouplingValidator,
    constraint: ConstraintValidator,
}

enum ValenceResolver { AtomTyping(AtomTypingValenceResolver), Counts(CountsValenceResolver) }
enum AromaticityResolver { HueckelRule(HueckelRuleResolver), Hmo(HmoResolver), Clar(ClarResolver) }
```

Validators are bare structs (no strategies). May revisit as enums later if they grow strategies.

**Removed in the restructure**: `ops/evaluate.rs` (depended on disabled `api`), `ops/chemistry.rs` (folded into `config.rs`), `ops/propagate.rs` (split into per-validator files), `ops/resolve.rs` (renamed `resolver.rs`), `ops/validate.rs` (renamed `validator.rs`). Existing `ops/valence.rs` (registry+tables only) decomposes into the `valence/` subdir.

**Disabled** in `ops.rs`: `ops/matcher.rs` (file kept on disk).

## Validator three-way split

`ConstraintValidator` is a separate beast from the tier-2 physics validators. Splitting constraint validation across the 6 entity types is a recipe for missing cases — centralize in one validator.

| Validator | Reads | Atom mode | Molecule mode |
|---|---|---|---|
| `ElectronInvariantValidator` | atom-intrinsic + topology incidence | constraints only (no topology) | constraints + topology |
| `SpinCouplingValidator` | `(unpaired, multiplicity)` parity | per-atom | per-atom over molecule |
| `ConstraintValidator` | constraints vs. data | n/a | every constraint vs. its topological/data counterpart |

`ConstraintValidator` also handles molecule-scope assertions like `:connected` — checked when present, never globally enforced. This is preferred over global topology validation, which would force every molecule to be connected.

(Naming: `ConstraintValidator` was preferred over `ConsistencyValidator` mid-discussion, then `ConsistencyValidator` was floated again near the end since the validator's scope grows past pure constraint-checking once length-mismatch detection is included. Unsettled.)

## Data-model change: per-atom electron contributions

The asymmetry made visible by the validator design: localized and dative bonds have two sources (constraints + topology) that must agree; aromatic systems and multicenter bonds had only one (the system-level `electrons: ValueAst` field). Resolution:

```rust
pub struct AromaticSystemAst {
    pub charge: ValueAst,                 // system-level
    pub spin: SpinStateAst,               // system-level
    pub electrons: Vec<ValueAst>,         // dense, parallel to relation atoms; Undetermined permitted
    pub constraints: AromaticSystemConstraints,  // gains ElectronCount(ValueAst)
}
```

Same shape on `MulticenterBondAst`.

Rationale (Cp−): each carbon is `C#c0#h1#a1` with no localized charge; the system-level `:charge -1` lives on the system. Same logic for tropylium C7H7+. This avoids the SMILES charge-localization problem while still allowing precise per-atom counting (the carbons are equivalent, just like benzene's). Spin is included for symmetry — applies to triplet arenes and similar excited-state structures.

System-level `electrons` was *not* moved out of the AST entirely — it stays as a per-atom vec (the `electrons: Vec<ValueAst>` above), and the *total* assertion moves to the constraint variant `AromaticSystemConstraint::ElectronCount(ValueAst)` (inlinable like `#e<n>`).

The `var relation indices are sorted and immutable per relation`, so a parallel `Vec<ValueAst>` keyed by position is sound — entry `i` belongs to atom `i` of the relation's `atoms` vec.

### DSL surface

```edn
{:atoms [0 1 2 3 4 5] :electrons [1 1 1 1 1 1] :type "#e6"}
```

Two flat vecs, position-aligned. Per-element follows the existing ValueAst EDN convention:

| ValueAst | EDN form |
|---|---|
| `Lit(n)` | bare `Edn::Int` — `1` |
| `Undetermined` | keyword — `:undetermined` |
| `LitSet(xs)` | bare vector — `[1 2 3]` |
| `Expr(_)` | quoted string — `"?n + 1"` |

`#e<n>` entity-string syntax stays — now writes the inline `ElectronCount` constraint, not the (now-vec) field.

`AromaticSystemDefaults::electrons` and `MulticenterBondDefaults::electrons` removed from `dsl/config.rs`. Defaulting a constraint by config is the wrong shape — defaults belong on field types where the field has a single canonical "no value" form, not on per-atom vecs whose Undetermined entries are the natural default.

### Validation semantics

- Length mismatch (`electrons.len() != atoms.len()`): not enforced at construction. Detected by validator, surfaced as `Contradictory` with a length-mismatch diagnostic.
- `Undetermined` entries: AST permits any subset. Per-atom electron equation returns `Underdetermined` for that atom. Sum-vs-`ElectronCount` returns `Underdetermined` if any entry is Undetermined.
- All-Undetermined: omit `:electrons` key entirely; parser auto-fills `vec![Undetermined; atoms.len()]`.

## AtomView back-reference (deferred)

Tier-2 validators need per-atom topology aggregates: σ valence sum, dative donor/acceptor counts, π contribution to aromatic systems, multicenter contribution. The first two require walking `graph` + the relation sets (back-ref to `&MoleculeAst` or equivalent). The last two are now reads on the per-atom `electrons` vec of the relation, no back-ref needed.

Ergonomics question: methods on `MoleculeAst` taking `AtomIdx` (`ast.valence(idx)`) vs. methods on a fattened `AtomView` carrying `&MoleculeAst` (`ast.atom(idx).valence()`). Decision: AtomView with back-ref. An attempt to land this in `umol-ast` was reverted; the placement and naming questions are open.

## Status at crash

**umol-ast (committed `fd543f68 Add atom contributions to aromatic systems and multcenter bonds`)**:
- `AromaticSystemAst.electrons: ValueAst` → `Vec<ValueAst>`
- `MulticenterBondAst.electrons: ValueAst` → `Vec<ValueAst>`
- `AromaticSystemConstraint::ElectronCount(ValueAst)` and multicenter analogue
- DSL `:electrons [...]` key on system / multicenter map entries (parse + render, tree + streaming)
- `dsl/config.rs`: `electrons` removed from defaults
- 1557 lib + 19 integration tests pass

**umol-graph (committed `9f81498a Start cleaning up umol-graph`)**:
- `api`, `ast`, `dsl` modules disabled in `lib.rs` (not deleted)
- `ops/` restructure designed; `ops/solution.rs` rewritten with `Solution<T, C>`; `ops/evaluate.rs` and `ops/matcher.rs` disabled
- Remainder of the layout (`ops/config.rs`, `ops/resolver.rs`, `ops/validator.rs`, `ops/valence.rs`, `ops/aromaticity.rs` plus subdirs) not yet implemented

## Status at close (2026-04-27)

The bottom-up implementation plan finished in nine phases. Final state:

- **Phase 0 — quarantine.** `lib.rs` and `ops.rs` reduced to the live core;
  bins / tests / benches gated behind a `legacy` Cargo feature.
- **Phase 1 — configs.** `ops/config.rs` (`ChemistryModel`, `ValenceModel`,
  `AromaticityModel`, `ElementScope`, `RingLimits`, `ConfigError`).
  `ops/valence/registry.rs`, `ops/valence/table.rs`. `AtomTypeRegistry` now
  stores `Vec<AtomAst>` directly via the new umol-ast `IntoAst` path.
- **AtomView extension.** `umol-ast::AtomView` gained an `&'a MoleculeAst`
  back-ref and the five chemistry method pairs per scheme A1
  (`bond_order_sum` / `valence_constraint`, `donated_pairs` /
  `donated_pairs_constraint`, `accepted_pairs` /
  `accepted_pairs_constraint`, `aromatic_contribution` /
  `aromatic_valence_constraint`, `multicenter_contribution` /
  `multicenter_valence_constraint`), plus `is_in_aromatic_system`. Other
  views deliberately untouched.
- **Phases 2 + 3 — validators.** `ops/validator.rs` with
  `ElectronInvariantValidator` (full impl, both `validate` and
  `validate_atom`), `SpinCouplingValidator` (stub), `ConstraintValidator`
  (stub), `EntityStructureValidator` (length checks), composite
  `Validator`. Method signatures take `impl AsRef<MoleculeAst>`. The
  electron-invariant equation keeps both `orbital_count` and
  `electron_count` as the two independent total-electron-per-atom
  accountings; renamed from the original `orbital`/`source` after
  surfacing the `thiserror` `source` field collision and clarifying that
  each side is meaningful in isolation. `impl AsRef<MoleculeAst> for
  MoleculeAst` added in umol-ast for transparent `&MoleculeAst` /
  `&Molecule` interop later.
- **Phase 4 — aromaticity resolver.** `ops/aromaticity.rs` plus the three
  algorithm modules ported off `crate::ast::*` to `umol_ast::ast::*`.
  Algorithm output is now `Vec<(Vec<AtomIdx>, AromaticSystemAst)>` directly
  — the legacy `AromaticSystem` perception-output wrapper is gone.
  Per-atom π contributions land in `AromaticSystemAst.electrons`. Bond-side
  `BondConstraint::Aromatic` is added by the dispatcher per induced bond.
  `RingFamily::InducedBenzenoid` proved unnecessary once the legacy
  coronene fixture (which had spurious cross-slice triangles) was rebuilt
  with real coronene topology — Vismara's relevant cycle basis returns the
  7 hexagonal faces directly. Each algorithm carries one error enum;
  `AromaticityResolver::resolve` classifies variants into Solution /
  Underdetermined / Err.
- **Phase 5 — valence resolver.** `ops/valence/atom_typing.rs`,
  `ops/valence/counts.rs`, private `shared.rs` for narrowing helpers.
  `NormalValenceTable` folded into `ValenceTable` via
  `ValenceEntry.normal_valence: Option<u8>`; the standalone
  `default-normal-valence-table.toml` deleted, the corresponding entries
  merged into `default-valence-table.toml`. Each algorithm has one error
  enum; composite `ValenceResolver` lifts variants via `From`.
- **Phase 6 — composite resolver.** `ops/resolver.rs` with `Resolver`,
  `ResolverContradiction`, `ResolverError`. One-shot per pass — no
  `ResolverCell`, no fixpoint, no `Progress`. `ops/chemistry.rs`,
  `ops/error.rs`, `ops/propagate.rs`, `ops/resolve.rs`, `ops/validate.rs`
  deleted.
- **Phase 7 — io re-enabled.** `umol-graph/src/table_ir/lift.rs` —
  `IntoAst<MoleculeAst> for &TableMolecule` (and per-atom and per-bond
  analogues) replaces the legacy `MoleculeAst::from_table_molecule`. SMILES
  and CTfile parser entry points return `Result<MoleculeAst, ParseError>`;
  resolver outcomes (Determined / Underdetermined / Contradictory) are
  classified at the parser boundary into the parser's error type.
- **Phase 8 — bins / tests / benches.** Ungated for everything except the
  resolution conformance suite and three benches (`morgan`,
  `molecule_dsl_parsing`, `substructure`) that depend on deleted modules.
  `cargo test -p umol-graph --features conformance`: 3253 lib + 12,677
  conformance assertions pass.
- **Phase 9 — cleanup.** Deleted `umol-graph/src/ast/`, `dsl/`, `api/`,
  module headers, `ops/evaluate.rs`, `ops/matcher.rs`. Cleaned up the
  comments in `lib.rs` and `ops.rs`.

Settled side decisions:

- **AtomView naming**: scheme A1 (`<topology>` / `<topology>_constraint`).
- **Validator naming**: `ConstraintValidator` (over `ConsistencyValidator`).
- **EntityStructureValidator** kept as a fourth, separate sub-validator —
  not folded into `ConstraintValidator`. Owns length-mismatch checks and
  any other structural invariants that aren't constraints.
- **ResolverCell removed.** Topology-invariance during narrowing means a
  cell type is unnecessary; `&mut MoleculeAst` is the only resolver
  argument. The `RingSet` cache moved onto `MoleculeAst` itself per the
  user's framing that AST is mutable during narrowing.
- **`Solution` retained.** CSP-flavored vocabulary (`Determined`,
  `Underdetermined`, `Contradictory`) stayed because there's no clean
  domain replacement and the three-valued shape is load-bearing.
- **No internal `Fault` union.** Each algorithm declares one error enum;
  the dispatcher classifies variants. No "Outcome" wrapper, no separate
  "RejectionReason" type per algorithm. Setup-vs-chemistry split happens
  at the dispatcher boundary, never inside an algorithm.

Items deferred to other docs:

- **Chemist-facing `Molecule` wrapper** (with caches and metadata
  question), **`Pattern`** with matcher caches, **`ReactionRule`**, the
  matcher port, transformations as ops, tier-3 validators, builder API,
  resolution conformance suite port. All filed under doc 86's
  "Outstanding" section.

Status: **Completed.**
