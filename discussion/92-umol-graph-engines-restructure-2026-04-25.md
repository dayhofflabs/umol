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

## Outstanding

- Implement the `ops/` layout above. Bottom-up order: `solution.rs` (done) → `config.rs` → `validator.rs` (validators don't mutate, easier first) → `valence.rs` + `valence/*` → `aromaticity.rs` + `aromaticity/*` → composites in `resolver.rs`.
- Resolve `ConstraintValidator` vs. `ConsistencyValidator` naming.
- AtomView back-ref — placement, naming, and minimum method set.
- Metadata on `Molecule`: round-trip preservation vs. purely semantic object.
- Re-enable `api`, then port to the new `ops/` shapes.
- `Matcher` re-enable (deferred to substructure / transformations work).
