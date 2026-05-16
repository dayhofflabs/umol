# Molecule AST API

Working notes on the public API around `MoleculeAst` and its derived views. Focus: where derived data is cached, how it surfaces to consumers, and which procedural operations should also have constraint counterparts in step 9 of doc 80.

## Terminology

**Bond** in this codebase means *localized two-atom bond* of positive integer order. This includes single, double, triple, quadruple bonds — i.e., not only the σ component, but also the π and (for quadruple bonds) δ components of a multiply-bonded pair. Localized-multicenter aggregates and delocalized aromatic systems are not bonds.

The codebase does not name them `LocalizedBond` for space economy, but the meaning is exactly that. In particular:

- **Not "sigma bond"** — a double bond is one `Bond`, not one σ + one π split into two entries. Counts like `bond_order_sum` and `connectivity` are localized-bond aggregates over this single entity.
- **Not "covalent bond"** — `DativeBond`, atoms participating in `AromaticSystem`, `MulticenterBond`, and even some `NoncovalentBond` types are *also* covalent in the chemistry sense. The covalent-vs-noncovalent dimension is orthogonal to the bond-vs-other-relation dimension.

The four **overlays** (`DativeBond`, `AromaticSystem`, `MulticenterBond`, `NoncovalentBond`) are explicitly not "bonds" in this vocabulary. They are separate relation types with their own participant lists and feature data, layered on top of the localized-bond graph. The umbrella term *overlay* is a codebase-internal usage — accurate (each is a typed n-ary relation over atoms, sitting above the basic atom + bond skeleton) and applicable uniformly to all four without contorting against any one of them.

**Structural assumptions on overlay × localized-bond overlap (provisional):**

- `AromaticSystem` is overlaid on both atoms *and* localized bonds — its `bonds()` method returns the localized bonds connecting its participants. Removing a participating atom *or* a participating bond invalidates the aromatic system's perceptual basis.
- `MulticenterBond` and `NoncovalentBond` should *not* share atom-pairs with localized bonds — these are alternative bonding representations, not enrichments of an existing bond. Not enforced at construction today; revisit if a use case surfaces that requires coexistence.
- `DativeBond` *may* coexist independently with localized bonds between the same atoms (e.g., a covalent backbone bond + a separate dative pair); mutation of one does not affect the other.

These assumptions drive the cascade rules for `MoleculeBuilder::remove` (§"Cascade semantics on `MoleculeBuilder::remove`").

Anywhere existing names use "sigma" or "covalent" as a substitute for "bond", they are misnamed. Flagged for the naming pass; do not propagate the misuse.

## Design context

Two anchors shape the design:

1. **General-chemistry scope.** umol targets organometallics, clusters, mixed-valence systems, and inorganic compounds alongside organics — see `project_general_chemistry_scope.md`. Cache decisions, the multicenter/dative/noncovalent relation types, and the invariant tiers (below) all reflect this scope rather than drug-like assumptions. RDKit-style "sanitization" rejects valid compounds outside that scope; we deliberately don't.

2. **Long-lived molecules.** Molecules enter once (parse/deser/transform), get queried and transformed many times. This justifies per-molecule cached views and influences the resolver's design — it can spend O(n) work per molecule that pays off across many queries.

## Concurrency

### `MoleculeAst` (current)

`MoleculeAst` carries a single-slot ring cache (`RingCache(Option<Box<RingCacheEntry>>)`, `umol-ast/src/ast/molecule.rs:54`) that is populated via `&mut self` in `mol.rings()` (canonical: Vismara relevant cycles, max_size 22). The cache is a plain `Option`, not `OnceLock` — mutation requires exclusive access; there is no interior mutability. Thread-safety story:

- `MoleculeAst: Send + Sync` per Rust's normal rules (no `Cell`/`RefCell`/raw pointers).
- Sharing across threads is via clone (cheap: most fields are `Arc`-wrapped) plus per-thread cache population, or via external synchronization (`Mutex` / `RwLock`) if a single shared instance is needed.
- `&MoleculeAst` in any thread is safe; cache writes need exclusive access by construction.
- The cache is excluded from `PartialEq` / `Hash` via custom impls, so `Molecule == Molecule` and hashing are independent of cache state.

The cache location on `MoleculeAst` is the current state; this may change if a future `Molecule` wrapper takes over caching responsibility (see §"AST-vs-API layering").

### `Molecule` / `Pattern` (when they return)

For the future chemist-facing wrappers (per §"AST-vs-API layering"):

- `Molecule: Send + Sync` via `Arc<MoleculeInner>`.
- Wrapper-side caches use `OnceLock` (init-once, read-many) so cloning is an `Arc` bump and clones share cache state.
- View contents (`RingSet`, `MatchTarget`, `MorganTarget`) must themselves be `Sync` — pure owned data, no interior mutability.
- `Pattern` follows the same shape. Coordinate annotations (when their location is settled) also must be `Sync`.

## Data structure hierarchy

### Usage assumption

Molecules are long-lived. They enter the system from external sources (SMILES/MOL parsing, DSL deserialization), get resolved once, and then are queried and transformed many times. Transformations (kekulization, tautomer enumeration, reactions) produce new `MoleculeAst` values that re-enter the same lifecycle; per doc 80 line 371 the cache does not transfer across rewrites.

### What's stable vs. changing during resolution

| Data | Status |
|---|---|
| Atoms set, bonds set, dative/noncovalent/multicenter relations | **fixed** before resolution starts |
| `AtomAst` / `BondAst` attribute fields | **narrowing** (`Undetermined` → `Lit`) |
| `aromatic_systems` relation | **appending** (perception discovers and adds tuples) |
| `constraints` vec | **draining** (entries discharged as solved) |

Topology is invariant across resolution. Topology-derived views (`RingSet`, biconnected components, distance matrix, atom degree, induced subgraphs) are valid throughout resolution and beyond. Attribute-derived views (Morgan-style fingerprints, packed `MatchTarget`) are valid only once attributes are concrete.

### Two long-lived consumer types

The chemist-facing/algebraic split is one axis; ground/non-ground is an orthogonal axis. Patterns (SMARTS-style substructure queries) are not transient resolution artifacts — they enter the system from external sources (SMARTS strings, DSL pattern syntax), persist for many matches against many targets, and have their own derived caches. Both `Molecule` and `Pattern` are long-lived `MoleculeAst` consumers.

| | `Molecule` | `Pattern` |
|---|---|---|
| Ground invariant | required | not required |
| Source | SMILES/MOL/DSL → resolve | SMARTS/DSL pattern syntax |
| Primary operation | "what is true about this?" | "does this occur in X?" |
| Useful cached views | `RingSet`, `BiconnectedComponents`, `DistanceMatrix`, `MatchTarget`, `MorganTarget`, fingerprint indexes | match-side scaffolding: per-atom constraint index, sub-pattern dependency graph, recursion order, packed pattern adjacency for VF2 |
| Lifecycle ends with | new molecule via transformation | new pattern via composition (rare) |

Cache contents are not the same. `Pattern` does not need `MorganTarget`; `Molecule` does not need a sub-pattern dependency graph. Sharing a single `Inner` shape would mean carrying dead `OnceLock` slots in both directions.

### Resolved hierarchy

```
MoleculeAst                       algebraic, ground or partial, single-slot ring cache
   │                              (cache location may move to wrapper if Molecule returns)
   ├─ Molecule                    ground invariant; chemistry-side caches
   ├─ Pattern                     no ground requirement; matcher-side caches
   └─ ResolverCell (private)      transient; transfers caches into Molecule on success
```

Three public types, but each represents a persistent concept the codebase cannot avoid:

- `Molecule` is "what does a chemist hold?"
- `Pattern` is "what does a SMARTS query become?"
- `MoleculeAst` is "what does the parser produce, the serializer consume, the transformation rewrite?"

Construction paths into the chemist-facing type:

- `Solver::resolve(ast, cfg) -> Result<Molecule, _>` — partial input, runs propagation, transfers populated topology caches into the result.
- `Molecule::new(ast) -> Result<Molecule, GroundError>` — already-ground input, only checks the invariant. For EDN literals and transformation results that the producer guarantees ground. Caches start empty; populate lazily on first access.

Both converge in a private `Molecule::from_inner(MoleculeInner { ast, ring_set, ... })`.

`Pattern` construction is simpler:

- `Pattern::new(MoleculeAst) -> Result<Pattern, PatternError>` — validates pattern well-formedness (sub-pattern anchors in range, constraints reference valid relations); no propagation. Caches start empty; populate lazily on first match.

## Operation boundaries: what each IO / transformation operates on

Current state (no chemist-facing wrapper exists; all ops operate on `MoleculeAst`).

### DSL IO

EDN parser produces `MoleculeDsl` as a boundary type — surface form with metadata (entity ids, atom aliases). `MoleculeDsl` raises to `MoleculeAst` via `IntoAst<MoleculeAst>` / `FromAst` in a separate step. Reverse: `MoleculeAst` lowers to `MoleculeDsl` via `FromAst<MoleculeAst>`, then renders to EDN. The two-step shape exists so the algebraic `MoleculeAst` carries no surface-form metadata.

### SMILES / MOL IO

Read: `SMILES/MOL → TableIR → MoleculeAst`. The lift step is `impl TryIntoAst<MoleculeAst> for &TableMolecule` (`umol-graph/src/table_ir/lift.rs:34`) — uses the same `IntoAst` trait family as the DSL → AST step. Per-atom and per-bond analogues exist in the same module.

After lift, `MoleculeAst` is passed to `Resolver::resolve(&MoleculeAst, &ChemistryModel) -> Solution` for narrowing. There is no `Molecule` wrapper to return into; resolved state lives on the (mutated) `MoleculeAst` itself. Validators (`ElectronInvariantValidator`, etc.) take `&MoleculeAst` directly.

Write: `MoleculeAst → TableIR → text`. No reverse `TryIntoTableIR` shortcut today; writers that need the table form go through whatever lowering exists per format.

### Transformation ops (kekulize, aromatize)

Operate on `&mut MoleculeAst` in place (or take `&MoleculeAst` and produce a new `MoleculeAst`). Implemented in `umol-graph/src/ops/transformer/`:

- **`Aromatizer`** (`ops/transformer/aromatizer.rs`) — perception + assignment.
- **`Kekulizer`** (`ops/transformer/kekulizer.rs`) — picks a Kekulé assignment via maximum matching on the aromatic subgraph.

Tautomer enumeration and canonical SMILES are not implemented at this point. The original `Molecule → Molecule` signatures from earlier in this doc are projections for the future chemist-facing wrapper, not current API.

### Reactions

`ReactionRuleAst` exists in `umol-ast` (`ast/reaction.rs`). Mechanical apply lives in `umol-ast/src/ast/molecule/rewrite.rs`: takes a `&ReactionRuleAst` and a target `&MoleculeAst` plus an `Assignment`, produces a rewritten `MoleculeAst`. No chemist-facing `ReactionRule` wrapper. Re-resolution after rewrite is up to the caller.

Semantic well-formedness of a full reaction rule (LHS + RHS + correspondence) — whether applying the rule to *any* matching target yields a chemically valid product — is **not** validatable statically in the general case. It depends on target-specific valence, charges, stereo. A narrow class (mass-balanced electron-pushing rules) is provably ground-preserving; the general case is not.

### Coordinates

Representation of coordinates in `MoleculeAst` is **undecided**. Likely belongs at a higher level than the AST (a separate annotation carried on a future `Molecule` wrapper, or a sidecar payload). The open question is how MOL / CXSMILES → DSL conversion threads the coordinate references — needed because DSL would have to either carry the coordinates or reference them externally, and that decision affects the boundary between AST and any wrapper.

Out of scope for the current AST API work. Status today: `umol-geometric` owns geometric primitives but no integration path with `MoleculeAst` is wired.

## API tiering: where MoleculeAst surfaces

### Two-tier surface

| Tier | Types | Audience | Entry points |
|---|---|---|---|
| 1 (chemist) | `Molecule`, `Pattern`, `ReactionRule`, `ReactionResult` | most users | `parse_smiles`, `parse_smarts`, `parse_smirks`, `Molecule::from_edn_str`, `Pattern::from_edn_str`, methods on those types |
| 2 (algebraic) | `MoleculeAst`, `AtomAst`, `BondAst`, relation types, `Solver`, `ResolverCell`, `ReactionRuleAst` | algorithm devs, custom transformations, debugging, programmatic builders, custom resolver config | `MoleculeAst::from_edn_str`, `parse_smiles_to_ast`, `Molecule::as_ast`, builder API, direct AST construction |

Tier 2 is public and documented, but tutorials and examples lead with tier 1.

### `Molecule` view API mirrors `MoleculeAst`

All view types from `MoleculeAst` (`AtomView`, `BondView`, `RingView`, etc.) are re-exposed on `Molecule` so users do not need to call `.as_ast()` for routine read access. The `Molecule` API is `MoleculeAst`'s view surface plus methods that consume cached precomputed attributes (`mol.morgan_fingerprint(r)`, `mol.find_matches(&pattern)`, `mol.rings()`).

### Where AST surfaces even at tier 1

1. **Equality and hashing.** Cache slots do not participate in identity. `impl PartialEq for Molecule` is `self.ast == other.ast`, `impl Hash for Molecule` hashes the AST. The AST is the canonical identity; the wrapper is sugar + caching. HashMap-of-molecules works without users thinking about AST, but the underlying definition is AST equality.

2. **EDN roundtrip is AST-level by design.** Per `project_molecule_dsl_roundtrip.md`, `MoleculeAst ↔ EDN` fidelity is non-negotiable. Implementation lives on `MoleculeAst::to_edn_str` / `MoleculeAst::from_edn_str`; `Molecule::to_edn_str` delegates via `self.as_ast().to_edn_str()`.

3. **Programmatic construction.** A builder API most naturally produces an AST, then resolves to a `Molecule`:
   ```rust
   let mol = MoleculeBuilder::new()
       .add_atom(Element::C)
       .add_atom(Element::C)
       .add_bond(0, 1, BondAst::aromatic())
       .build()?           // -> MoleculeAst
       .resolve(&cfg)?;    // -> Molecule
   ```
   AST is observable mid-build for users who want to inspect or mutate before resolving.

### Parser entry-point convention

Each parser exposes:

```rust
// Tier 1: parse + resolve with default config
pub fn parse_smiles(s: &str) -> Result<Molecule, SmilesError>;
pub fn parse_smarts(s: &str) -> Result<Pattern, SmartsError>;
pub fn parse_smirks(s: &str) -> Result<ReactionRule, SmirksError>;

// Tier 2: parse only, return AST (custom resolver, debugging, AST transforms)
pub fn parse_smiles_to_ast(s: &str) -> Result<MoleculeAst, SmilesError>;
pub fn parse_smarts_to_ast(s: &str) -> Result<MoleculeAst, SmartsError>;
pub fn parse_smirks_to_ast(s: &str) -> Result<ReactionRuleAst, SmirksError>;

// Configurable: explicit resolver config
pub fn parse_smiles_with(s: &str, cfg: &ResolverConfig) -> Result<Molecule, SmilesError>;
```

### Implication for naming

The two-tier framing supports:
- `MoleculeAst` — algebraic, second-tier; `Ast` suffix makes the role explicit.
- `Molecule` — chemist-facing, first-tier; the name a chemist expects.
- `Pattern` — query-facing, first-tier.

`MoleculeAst` is not hidden, just not what users type day-to-day.

## MoleculeAst operation taxonomy

This section consolidates the public operation surface on `MoleculeAst` and its views, organizing by intent. Each row carries a state tag:

- **Impl** — exists in code, verified.
- **Designed** — shape settled in a doc, no impl yet.
- **Open** — no decision; surfaces a real need from doc / ops survey but the API shape is undecided.

Names tagged `[naming TBD]` are working-name placeholders deferred to a separate naming pass after the surface stabilizes. References point to authoritative specs (other discussion docs) or call sites (file:line where the need recurs).

### Design principles for the AST surface

Three cross-cutting rules the taxonomy assumes throughout.

**1. Molecule-scope indices, not graph ids.** Public methods on `MoleculeAst` and its views speak `AtomId` / `BondId` / `DativeBondId` / `AromaticSystemId` / `MulticenterBondId` / `NoncovalentBondId`. Graph-level identifiers (`NodeId`, `EdgeId`, `RelationId`) belong to `umol-graph-core` and stay behind the `raw_graph() -> &Graph` escape hatch. Every method that wraps a graph primitive performs the index-space conversion at the boundary; no pure pass-through.

**2. View field naming.** Inside entity views (`AtomView`, `BondView`, `DativeBondView`, etc.) the current fields `data: &EntityAst` and `ast: &MoleculeAst` (the back-pointer just added) are misnamed. Settled convention:

- `data: &EntityAst` → `ast: &EntityAst` — public field. The *entity's* AST is the natural meaning of "ast" inside the view. Used as the escape hatch for whole-struct operations (cloning, comparison, low-level inspection). Not duplicated as a method; the field IS the access path.
- `ast: &MoleculeAst` → `molecule: &MoleculeAst` — private field. Internal navigation back-pointer; not part of the public surface.

Mixed field-vs-method surface on `AtomView` (and the other entity views):

- **Fields**: `pub id: AtomId`, `pub ast: &EntityAst` — structural primitives of the view. Direct field access for both.
- **Methods**: per-field accessors (`element()`, `charge()`, ...), derived readers (`valence()`, `degree()`, ...), predicates (`is_in_ring()`, ...), cross-entity navigation (`neighbors()`, `rings()`, `aromatic_system()`, ...). Methods for everything computed or naming-an-operation.

The line: field for "what this view *is* (id + entity-AST handle)"; method for "what this view *does* (reads, derivations, predicates, navigation)".

Mut views (`AtomViewMut`, etc.) follow the same shape with `pub ast: &mut EntityAst`. Auto-deref handles both read and write through the field; no `ast()` / `ast_mut()` method pair needed.

**3. Topology-only caches OK on the AST; attribute-derived caches go on the wrapper.** Topology is invariant across resolution (per §"What's stable vs. changing during resolution"); a cache that depends only on `(atoms, bonds)` (current ring cache, future `BiconnectedComponents`, `DistanceMatrix`) stays valid through every in-place mutation and is safe on `MoleculeAst`. Caches that depend on attribute concretion (Morgan fingerprint, packed `MatchTarget`) need a ground-only home — the future `Molecule` wrapper. The current `MoleculeAst::rings` slot fits the topology-only rule; no immediate need to relocate.

**4. Cross-molecule index misuse is not validated.** Index types (`AtomId`, `BondId`, `DativeBondId`, `AromaticSystemId`, `MulticenterBondId`, `NoncovalentBondId`, `RingId`) are bare newtypes over `u32`; they carry no provenance tag identifying which molecule they were obtained from. Methods that consume an index obtained from a different `MoleculeAst`, view, or `RingSet` silently produce wrong answers (out-of-range panic at best, semantically incorrect chemistry at worst) rather than rejecting the call at the boundary. Examples:

- Passing one molecule's `AtomId` to another molecule's `atom()` / `atom_view`.
- Passing a `&RingSet` enumerated against molecule A into a parametric ring query on molecule B's atom view (a real hazard for the parametric API in §"Ring access").
- Indexing into one molecule's `electrons` array with a position derived from another molecule's aromatic-system view.

The pattern matches RDKit and the rest of the cheminformatics ecosystem. Lifetime-branded indices (`'mol` parameter, `GhostCell`-style typestate) would catch the misuse at compile time but at substantial ergonomic cost; the consequence of misuse is incorrect chemistry, not memory unsafety, and we accept the hazard. Callers are responsible for keeping index provenance consistent.

**5. `ValueAst` arithmetic carries `Option<i64>` semantics; `Undetermined` ≡ "no constraint".** Two interlocking simplifications that unify derived-quantity computation and constraint inspection:

*Arithmetic*. `ValueAst` supports `Add` / `Sub` / `Mul` / `Div` (binary and scalar). Operations on two `Lit`s yield a `Lit`; anything involving a non-`Lit` (`Undetermined`, `LitSet`, `Expr`) collapses to `Undetermined`. The Option-isomorphism is:

| `Option<i64>` | `ValueAst` |
|---|---|
| `Some(n)` | `Lit(n)` |
| `None` | `Undetermined` |
| `lhs.zip(rhs).map(\|(a,b)\| a+b)` | `lhs + rhs` (collapses non-Lit) |

This isn't symbolic arithmetic — `Expr(charge + 4) + Lit(2)` collapses to `Undetermined`, not to `Expr(charge + 6)`. Reflective of the resolver's "narrow to `Lit` before computing" assumption; full symbolic resolution is a separate future pass.

*Constraint absence ≡ `Undetermined`*. Per-kind constraint accessors on `*Constraints` stores return `ValueAst` (not `Option<ValueAst>`). "No constraint of this kind asserted" returns `ValueAst::Undetermined`, identical to "constraint asserted as `Undetermined`". Justified by the existing render rule (Undetermined constraints elide on output, so storage-presence is already not load-bearing) and the wildcard semantics of `Undetermined` in `*Ast::matches` (memory: `feedback_undetermined_is_wildcard`).

Together these unify derived and constraint-asserted reads:

```rust
view.valence()                    -> ValueAst  // Lit or Undetermined
view.constraints().valence()      -> ValueAst  // Lit, LitSet, Expr, or Undetermined
view.constraints().valence().matches(&view.valence())  // direct comparison; always works
```

Marker constraints (those without a value payload — `BondConstraint::Aromatic`, `InRing`, etc.) stay `bool`-shaped: `bond.constraints().aromatic() -> bool`, `atom.in_ring() -> bool`. The ValueAst-collapse rule applies only to value-bearing constraint kinds.

*Edge cases for `ValueAst` arithmetic*:

- **Divide-by-zero**: panic. Programmer error; no chemistry-level scenario justifies a silent collapse.
- **Overflow**: panic. Chemistry counts don't approach `i64` bounds in practice; reaching the bound indicates a bug, not a degenerate input.
- **Negative results**: allowed. `Lit(1) - Lit(4) = Lit(-3)`. Whether a negative result is meaningful for a given quantity (e.g., a negative bond count is nonsense) is a semantic concern at the call site, not an arithmetic one.

*Parallel ASTs*:

- **`ImplicitHydrogensAst`** gets the same arithmetic and collapse rule. `Normal` (the "default-implicit-H" meta-state) collapses to `Undetermined` under arithmetic, identical to bare `Undetermined`.
- **`IsotopeAst`** does *not* get arithmetic. Isotope masses behave more like enum tags than numerics — adding mass numbers across atoms is meaningless outside nuclear chemistry, which is out of scope. Revisit if nuclear chemistry becomes a real consumer.

**6. Mutation goes through an `Edit` vocabulary; checked transactions are the load-bearing primitive.** Caller-facing mutations compose as `Vec<Edit>`. Checked mutation goes through `MoleculeBuilder::transact(Vec<Edit>) -> Result<Transaction, TransactionError>`, which records realized `Undo` entries and rolls back by reverse replay on failure. `transact_unchecked(Vec<Edit>) -> ()` is the separate no-journal path for generated edits known to be correct.

This buys five capabilities that direct-mutation can't:

- **Atomicity**: a batch either applies entirely or rolls back (no partial-mutation states for callers to clean up).
- **Transient invariant violation**: intermediate states inside a transaction can be temporarily tier-1-invalid (DPO needs this — between L\\K removal and R\\K addition the molecule is structurally incomplete).
- **Reified cascade**: overlay-removal cascades induced by topology removal are recorded in realized `Undo`; rollback undoes them.
- **Serializable / replayable / undoable mutation**: `Edit` remains pure caller-facing mutation data; checked application produces the realized `Undo` data needed for exact rollback.
- **Validation gate**: tier-1 invariants checked on commit, not per-Edit.

Detail design lives under §"Mutation operations" → §"Edit vocabulary and transactions"; doc 43 has the prior-art patch design that this builds on.

### Naming conventions

Captured for reference after the long naming pass. Most rules emerge from decisions documented in the surrounding sections; consolidated here for cross-referencing.

#### Method-name conventions

- **`is_<predicate>()`** for boolean checks: `is_in_ring()`, `is_in_aromatic_system()`, `is_ground()`, `is_empty()`. Plain `<predicate>()` (no `is_` prefix) is reserved for non-boolean readers.
- **`<entity_attribute>()`** for scalar/value readers on a view: `valence()`, `ring_count()`, `charge()`. Singular method name matches the attribute label; the return type carries the arity (`ValueAst`, `usize`, etc.).
- **`<plural_entity>()`** for collection iteration: `atoms()`, `bonds()`, `rings()`, `aromatic_systems()`. Iterator returns `Iterator<EntityView<'_>>`.
- **`<modifier>_<entity>(args)`** for filtered/qualified collections: `overlapping_rings(set)`, `connecting_bond(a, b)`, `induced_bonds(atoms)`, `dative_bonds()` (on `atom_view`). The modifier qualifies which entities are returned.
- **Both modifier-entity (`connecting_bond`) and entity-relation (`bond_between`) patterns exist in English** but we standardize on modifier-entity for the molecule-AST surface to reduce the patterns callers must learn.

#### Arity and return-type conventions

- **Views by default; `_id` / `_ids` suffix for index-returning companions.** `mol.connecting_bond(a, b) -> Option<BondView>`; companion `mol.connecting_bond_id(a, b) -> Option<BondId>`. Collections: `relation_view.atoms() -> Iterator<AtomView>`; companion `relation_view.atom_ids() -> &[AtomId]`. The default is the chaining-friendly form; index access is the explicit opt-in.
- **`Option<T>` for at-most-one results.** Used where the data model enforces the constraint (per-atom `aromatic_system()`, `connecting_bond(a, b)`). Iterator-returning methods naturally cover 0–* arities; `Option` for 0–1 makes the constraint legible from the signature.
- **`Iterator<T>` for collections** (0–* including 0–1 when the data model doesn't enforce uniqueness).
- **`Vec<T>` only when allocation is integral to the operation** (e.g., `transact(Vec<Edit>) -> Transaction`). Otherwise prefer `impl Iterator<Item = T>` for laziness.

#### Suffix conventions

- **`_with(args)`**: default-having method with explicit parameter override. `mol.rings_with(family, max_size, filter)` (vs `mol.rings()` which uses canonical defaults). Pattern matches Rust stdlib (`Vec::with_capacity`, `HashMap::with_hasher`).
- **`_from(&SourceCollection)`**: scope-restricted to a caller-supplied collection. `atom_view.rings_from(&RingSet)`, `atom_view.is_in_ring_from(&RingSet)`. Reads as "drawn from this set". (Mixed use of `_in` is acceptable but `_from` is preferred for clarity.)
- **`_all(kind)`**: multi-valued retrieval on stores with the `is_unique()`-aware add path. `*Constraints::get_all(kind) -> Iterator<&Constraint>`.
- **`for_<entity>(idx)` / `by_<discriminant>(kind)`**: iterator-filter methods on the `Constraints` collection. `mol.constraints().for_atom(idx).by_kind(kind)`. Composable.
- **`_id` / `_ids`**: index-returning companion of a view-returning method, per P1 above.

#### Collection-type conventions

Two distinct shapes for "give me this set of entities", chosen by structural fit:

- **`*Views<'a>` wrapper struct** for canonical molecule storage. Borrowed window over `MoleculeAst`'s entity vec; provides `.iter()` (yields `*View`), `.count()` (O(1)), `.ids()` (yields `*Id`), `.get(idx) -> Option<*View>`. Examples: `mol.atoms() -> AtomViews<'_>`, `mol.bonds() -> BondViews<'_>`, `mol.dative_bonds() -> DativeBondViews<'_>`, plus the three other overlay collections.

- **Bare `impl Iterator<Item = *View<'_>>`** for derived / filtered streams. One-shot; no O(1) count or random access. Examples: `atom_view.neighbors()`, `atom_view.dative_bonds()` (atom → relations is a filter over the full dative-bond set), `atom_view.rings()`, `aromatic_system_view.atoms()`.

The split reflects cost: `*Views` over canonical storage gives O(1) count and random access for free; derived streams don't have that property. Wrapping every derived iterator in a `Views`-like struct would add a type per accessor without buying anything.

**Exceptions where the collection type is named differently** — principled, tracking structural differences:

- **`RingSet`** owns its data (it's the perception artifact); `*Views` types borrow from molecule storage. Naming distinction follows ownership distinction.
- **`Constraints`** (the molecule-list flat-vec) holds mixed-kind entries; `*Views` types are kind-uniform with idx-keyed random access. The flat-vec needs different access patterns (iterator-filter, not idx lookup).
- **Neighbors** has no wrapper type — bare iterator only. Could grow a `Neighbors<'a>` wrapper if O(1) count becomes a hot-path concern; not justified by current usage.

#### DSL symbol conventions

- **Lowercase** for atom-intrinsic per-atom properties: `#c` (charge), `#i` (isotope), `#h` (implicit H), `#l` (lone pairs), `#u` (unpaired), `#s` (multiplicity), `#v` (localized valence), `#a` (aromatic valence), `#m` (multicenter valence), `#d` (donated pairs), `#t` (accepted pairs).
- **Uppercase** for derived totals / SMARTS-compatible aggregates: `#D` (degree), `#X` (total degree / SMARTS `X`), `#H` (total hydrogens), `#R` (ring count), `#V` (total valence per principle 5).
- **Lowercase, SMARTS-faithful** for narrow-scope quantities: `#r` (ring size, SMARTS `r<n>`), `#x` (ring degree, SMARTS `x<n>`).
- **Lowercase, new (no SMARTS analog)**: `#y` (ring valence, new for systematic naming alongside `#x`).
- **Symbol sugar**: `+` for "at least 1" (`#R+`, `#a+`, `#m+`), `!` for "not" or zero (`#a!`, `#R!`). See per-predicate parser definitions in `umol-ast/src/dsl/predicates.rs`.

### Read operations

#### Pure-graph (adjacency) reads

Operations on the underlying graph (atoms = nodes, bonds = edges) that ignore entity feature data. Per principle 1, all operations use molecule-scope indices; the raw `Graph` is reachable only via `raw_graph()`.

| Operation | State | Notes |
|---|---|---|
| `mol.graph() -> GraphView<'_>` | **Impl** | molecule-scope graph operations; constructed at `views.rs:682` |
| `mol.raw_graph() -> &Graph` | **Impl** | escape hatch for callers that want graph-id-space algorithms |
| `GraphView::degree(atom: AtomId) -> usize` | **Impl** | `views.rs:691` |
| `GraphView::connected_components(alg) -> Vec<Vec<AtomId>>` | **Impl** | `views.rs:695` |
| `GraphView::biconnected_components(alg) -> Vec<Vec<AtomId>>` | **Impl** | `views.rs:703` |
| `GraphView::shortest_cycle_through_bond(bond, alg) -> Option<usize>` | **Impl** | `views.rs:714` |
| `GraphView::connected_components_in(&[AtomId]) -> Vec<Vec<AtomId>>` | **Open** | restricted to a node subset (π-subgraph case in `ops/aromaticity/hmo.rs:85-107`) |

Atom-local adjacency reads live on `AtomView`, not on `MoleculeAst` directly:


| Operation | State | Notes |
|---|---|---|
| `atom_view.neighbors() -> impl Iterator<Item = NeighborView<'_>>` | **Impl (reshape pending)** | `views.rs:103`; new shape: `NeighborView { atom_idx, bond_idx, molecule }` — all fields private; `nbr.atom() -> AtomView` and `nbr.bond() -> BondView` accessor methods (`#[inline]`); raw indices accessible via `nbr.atom().idx` / `nbr.bond().idx`. View construction stays lazy so iteration cost is one back-pointer copy per item, lookups deferred to whichever side the caller actually consumes |

Localized-bond lookup by endpoint pair and induced-bond enumeration are categorized as inter-entity reads (atoms-in, bond-out); see §"Inter-entity derived reads" below.

#### Internal relation reads

Per-relation accessors on `DativeBondView`, `AromaticSystemView`, `MulticenterBondView`, `NoncovalentBondView`. Read the participant atom list, derived bond list (atom-pairs in the relation that also have a localized bond), and any per-participant parallel data (`electrons` on aromatic / multicenter).

| Operation | State | Notes |
|---|---|---|
| `relation_view.atoms() -> impl Iterator<Item = AtomView<'_>>` | **Impl (return-type pending)** | participant list as views (P1: views by default); currently returns `&[AtomId]`. `_ids` companion: `relation_view.atom_ids() -> &[AtomId]` retains direct slice access |
| `relation_view.bonds() -> impl Iterator<Item = BondView<'_>>` | **Impl (return-type pending)** | localized bonds induced by participants as views; `_ids` companion: `relation_view.bond_ids() -> impl Iterator<Item = BondId>` |
| `aromatic_view.electrons[pos]`, `multicenter_view.electrons[pos]` | **Impl** | per-participant parallel array |
| `view.constraints() -> &EntityConstraints` | **Open** | per principle 2, expose as direct method on the view so callers don't reach through to the inner AST; the returned store has per-kind named accessors (see below) |

Degree-style and order-summing aggregates on relation views (`AromaticSystemView::degree()`, `::heavy_atom_degree()`, `::valence()`) are categorized as inter-entity derived reads; see that section below.

#### Entity feature reads

Per-entity feature data on `AtomAst`, `BondAst`, `DativeBondAst`, `AromaticSystemAst`, `MulticenterBondAst`, `NoncovalentBondAst`. Accessed via the `*View` types returned from `mol.atoms() / bonds() / dative_bonds() / aromatic_systems() / multicenter_bonds() / noncovalent_bonds()`.

| Operation | State | Notes |
|---|---|---|
| `mol.atoms() -> AtomViews`, `mol.bonds()`, ... | **Impl** | view collections; `*Views::iter()`, `get(idx)`, `count()` |
| `view.ast` (field, public) | **Settled** | direct field access; escape hatch for whole-entity-AST operations (cloning, structural comparison, low-level inspection). Replaces the proposed `view.ast()` method — field is the access path, no redundant accessor. For mut views, `view_mut.ast: &mut EntityAst` handles both read and write via auto-deref, eliminating the `ast()` / `ast_mut()` method pair that the methods-route would require |
| Per-field accessors on each view | **Open** | `#[inline]` getters returning `&T` for AST fields, by-value for `Copy` primitives; pattern below |

##### Per-field accessor pattern

For each entity view, add one `#[inline]` method per field on the inner AST returning `&T` (or `T` if `Copy`-cheap). Callers stop reaching through `.data` / `.ast` for routine reads; the field access stays available only via `view.ast()` for whole-struct operations.

The accessor methods sit uniformly alongside the derived-quantity readers from §"Inter-entity derived reads" — at the call site, `atom_view.charge()` (stored), `atom_view.total_valence()` (derived from incident bond orders + implicit H), and `atom_view.in_ring()` (derived against the canonical ring set) all look identical. Callers don't have to mentally distinguish "direct field" from "computed quantity"; the implementation distinction is hidden, which also leaves room to migrate a field to a derived quantity (or vice versa) without breaking call sites.

| View | Accessors |
|---|---|
| `AtomView` | `element()`, `isotope_mass()`, `charge()`, `implicit_hydrogens()`, `lone_pairs()`, `spin()`, `constraints()` |
| `BondView` | `order()`, `charge()`, `spin()`, `constraints()` |
| `DativeBondView` | `acceptor_slot()`, `order()`, `constraints()` |
| `AromaticSystemView` | `electrons()`, `charge()`, `spin()`, `constraints()` |
| `MulticenterBondView` | `electrons()`, `charge()`, `spin()`, `constraints()` |
| `NoncovalentBondView` | `kind()`, `constraints()` |

Return-type convention: `&T` uniformly for AST-typed fields (`&ValueAst`, `&ElementAst`, `&IsotopeAst`, `&SpinStateAst`, `&EntityConstraints`). `Copy` primitives (if any future field) return by value. `#[inline]` everywhere — semantically a field read, optimizes to zero-cost.

The `constraints()` row already appeared under §"Internal relation reads" and §"Constraint reads"; both refer to the same accessor.

Mutation-side parallel (`atom_view_mut.charge_mut() -> &mut ValueAst`, etc.) deferred — current mutation through `view_mut.data.<field>` (or `view_mut.ast.<field>` post-rename) is workable; revisit if an ergonomic complaint surfaces.

#### Cross-entity navigation

Reverse direction: from an atom (or bond), enumerate the relations that contain it.

**Working convention on shape (provisional).** The iterator-returning method is the primary form; count queries use `.count()` directly. For boolean ("is in any relation of this type") tests, the canonical iterator idiom is `.next().is_some()`, which is correct but ergonomically awkward; named boolean shortcuts (`atom_view.is_in_aromatic_system()`, and symmetric methods for dative / multicenter / noncovalent) are kept as sugar over the iterator emptiness check.

This is provisional. If the relation count grows or the symmetric boolean methods multiply, revisit before release — possibly drop the named booleans and standardize on the iterator idiom, or introduce a small extension trait giving `is_empty()` on iterators.

All iteration items are `*View`, not `*Id`, consistent with the `NeighborView` reshape convention — back-pointers are already on the views; yielding indices would force callers to re-look up the views every iteration step.

| Operation | State | Notes |
|---|---|---|
| `atom_view.aromatic_system() -> Option<AromaticSystemView<'_>>` | **Open** | at-most-one (per doc 52 perception design); singular Option is the right shape — singular Iterator over 0–1 elements would be a footgun (one-character-typo collision with the plural-iterator forms for other overlays) |
| `atom_view.is_in_aromatic_system() -> bool` | **Impl** | binary shortcut (`atom_view.aromatic_system().is_some()`); kept because hot |
| `atom_view.dative_bonds() -> impl Iterator<Item = DativeBondView<'_>>` | **Open** | possibly many |
| `atom_view.multicenter_bonds() -> impl Iterator<Item = MulticenterBondView<'_>>` | **Open** | possibly many |
| `atom_view.noncovalent_bonds() -> impl Iterator<Item = NoncovalentBondView<'_>>` | **Open** | possibly many |
| `bond_view.aromatic_system() -> Option<AromaticSystemView<'_>>` | **Open** | bond is also in at-most-one aromatic system by perception design; symmetric Option |
| `bond_view.is_in_aromatic_system() -> bool` | **Open** | symmetric boolean shortcut |
| `dative_bond_view.aromatic_system()`, etc. (relation → other relations) | **Open** | stub out for now, mirror the atom-side shape if a need arises |

The forward direction (`relation.atoms()`) is fully impl. Reverse-direction recurs in valence and aromaticity perception. The aromatic-system singular Option asymmetry tracks the data-model constraint (perception assigns each atom / bond to at most one aromatic system); if that ever relaxes, the API breaks loudly rather than silently (Option-to-Iterator type change).

#### Ring access (perception-derived collection)

Rings are not stored on the AST — they are computed from graph topology by an enumerator. The result (`RingSet`) is a perception artifact whose shape depends on caller-chosen parameters: ring family, maximum ring size, optional atom filter. Two distinct uses:

- **Canonical-rings semantics**: a single fixed choice of (family, max_size) that all view-side ring methods and constraint evaluators (`InRing`, `RingSize`, `RingCount`) refer to. Hard-coded so that two molecules with identical structure never give different answers to "is this atom in a ring" depending on hidden configuration.
- **Non-canonical procedural use**: aromaticity perception algorithms (HMO, Clar) and other future passes (chirality perception in multicyclic systems) that want a specific family / size cap / atom filter. Available through `mol.rings_with(family, max_size, filter)` returning an owned `RingSet`; uncached, never surface in the DSL or constraint vocabulary.

##### Canonical commitment

The canonical ring set is **Vismara relevant cycles** (the unique, ordering-independent set of all cycles that participate in at least one minimum cycle basis), enumerated up to a fixed maximum size of **22 atoms**.

- *Why Vismara*: well-defined and unique (atom-ordering / algorithm-independent), unlike SSSR. Practical on highly connected systems (linear in molecule size after BCC decomposition), unlike enumerating all simple cycles (C60 has ~10^25 simple cycles).
- *Why 22*: porphyrin (an 18-π-electron macrocycle with formal ring size 16) plus headroom for common metalloporphyrin frameworks. Crown ethers, macrocyclic peptides past this size go through the non-canonical procedural path.
- *Why hard-coded*: the canonical commitment is part of the SMARTS-frontend contract. Configurable canonical params would create the RDKit-class problem where two implementations both claim "in a ring" but disagree silently.

The existing variant `RingFamily::Induced` is misleadingly named — Vismara relevant cycles are not necessarily induced (chordless) cycles in the graph-theoretic sense. **Rename to `RingFamily::Relevant`** in the naming pass.

Non-canonical configuration sits on `ChemistryModel` (for per-resolver-pass choices); the DSL has no parameter for it. A schema-like DSL extension for non-canonical ring constraints could be considered later if a real use case emerges, but is deliberately out of scope now.

##### Caching strategy

**Cache the canonical answer only.** The single slot on `MoleculeAst` holds the canonical `RingSet` (Vismara relevant cycles, max_size 22, no atom filter) — populated lazily on first access. Non-canonical or filtered enumerations bypass the cache entirely and return owned `RingSet` to the caller, who manages amortization themselves.

This resolves the caching-vs-filtering asymmetry: the filter is a closure (`impl Fn(AtomId) -> bool`) and isn't hashable, so it can't participate in a cache key anyway. Rather than pretend the parametric cache covers the no-filter subset of a more general API, the API names the actual contract — cache covers canonical; everything else is caller-managed.

Callers that need a non-canonical ring set across many queries (a perception algorithm doing thousands of filtered ring lookups) hold their own `RingSet`:

```rust
let pi_rings = mol.rings_with(RingFamily::Relevant, 22, |a| pi_atoms.contains(&a));
// reuse pi_rings across local queries
```

Same effect as a parametric cache; lifetime is explicit.

##### API surface

| Operation | State | Notes |
|---|---|---|
| `RingEnumerator::enumerate(&MoleculeAst, family, max_size, filter) -> RingSet` | **Impl** | underlying procedural primitive |
| `mol.rings() -> &RingSet` | **Designed** | canonical (Vismara relevant cycles, max_size 22, no filter); cached single-slot on `MoleculeAst`; `&mut self` for lazy init |
| `mol.rings_with(family: RingFamily, max_size: usize, atom_filter: impl Fn(AtomId) -> bool) -> RingSet` | **Designed** (rename from `enumerate_rings`) | uncached, owned; caller manages amortization |
| `RingSet::iter()`, `get`, `count`, `ids`, `family`, `scope`, `max_ring_size` | **Impl** | collection accessors |
| `RingSet::contains_atom`, `contains_bond`, `atom_smallest_ring_size`, `bond_smallest_ring_size` | **Impl** | per-atom / per-bond derived |
| `RingSet::ring_count_at(atom) -> usize` | **Open** (designed) | symmetric with `atom_smallest_ring_size`; `RingCount(ValueAst)` constraint counterpart designed not impl |
| `RingSet::overlapping_rings(&[AtomId]) -> impl Iterator<Item = RingView>` | **Open** | rings sharing any atom with the subset |
| `RingSet::overlapping_atoms(&[AtomId]) -> impl Iterator<Item = AtomView<'_>>` | **Open** | atoms in the subset that participate in at least one ring (P1: views). `_ids` companion available if needed |
| `RingSet::overlapping_bonds(&[BondId]) -> impl Iterator<Item = BondView<'_>>` | **Open** | bonds in the subset that participate in at least one ring |
| `RingSet::relation`, `are_spiro/fused/bridged`, `*_neighbors`, `fused_component(s)`, `shared_atoms`, `shared_bonds` | **Impl** | pairwise ring-vs-ring; **no constraint variants** by design |

Canonical view-side sugar (entity-centered queries against the canonical ring set):

| Operation | State | Notes |
|---|---|---|
| `atom_view.is_in_ring() -> bool` | **Open** | boolean predicate; equivalent to `atom_view.rings().next().is_some()` |
| `atom_view.ring_count() -> ValueAst` | **Open** | count of canonical rings containing this atom; paired with `RingCount(ValueAst)` constraint |
| `atom_view.ring_size() -> impl Iterator<Item = usize>` | **Open** | multi-valued: sizes of all canonical rings containing this atom; paired with `RingSize(ValueAst)` constraint (interpretation B) |
| `atom_view.rings() -> impl Iterator<Item = RingView<'_>>` | **Open** | iterator over canonical rings containing this atom; yields `RingView`, not `RingId` |
| `bond_view.is_in_ring()`, `ring_count`, `ring_size`, `rings()` | **Open** | symmetric; bond-side analogs |
| `aromatic_system_view.overlapping_rings()`, `.overlapping_atoms()`, `.overlapping_bonds()` | **Open** | canonical-keyed sugar over the `RingSet` overlap queries |

`smallest_ring_size()` is dropped — callers compose `atom_view.ring_size().min()`. The multi-valued `ring_size()` exposes all containing ring sizes; the smallest is one of many possible reductions.

##### Ring size constraint semantics (interpretation B)

`RingSize(ValueAst)` and DSL `#r<n>` follow **interpretation B**: the constraint matches iff the atom (or bond) is in *some* canonical ring whose size matches the `ValueAst` pattern. Examples for an indene 5/6-ring junction atom:

- `RingSize(Lit(5))` ↔ `#r5`: matches ✓ (junction atom is in a 5-ring).
- `RingSize(Lit(6))` ↔ `#r6`: matches ✓ (junction atom is also in a 6-ring).
- `RingSize(LitSet([3,4]))`: doesn't match (not in any 3- or 4-ring).
- `RingSize(Expr("r >= 6"))`: matches ✓ (in a 6-ring).

Equivalent reader-side formulation:

```rust
ring_size_pattern.matches_any(atom.rings().map(|r| r.len()))
```

**Divergence note**: Daylight's SMARTS reference phrases `r<n>` as "in *smallest* SSSR ring of size n", which reads as a stricter interpretation A ("the smallest containing ring has size n"). In practice most implementations (RDKit, OpenBabel) implement interpretation B silently; we adopt B as the explicit, documented choice. Under B, an indene junction atom matches both `#r5` and `#r6`, which is what users expect from looking at the structure.

`atom_smallest_ring_size()` is kept as a separate utility — useful for chemical environment classification ("is this a 5- or 6-membered ring atom?") — but it's **not** the constraint counterpart to `RingSize`. The constraint counterpart is the rings iterator + the pattern's `matches_any` against ring sizes.

Translation to Vismara relevant cycles: the smallest containing ring size is graph-theoretically invariant (same answer under SSSR or Vismara), but "atom is in some ring of size n" can differ in fused/symmetric systems where Vismara includes additional rings beyond a given SSSR choice (cubane, etc.). For those systems our answers may differ from a specific SMARTS engine's SSSR-based answer — by design, since Vismara is the canonical-rings commitment.

##### Ring size as a multi-valued constraint kind

Interpretation B makes `RingSize` the first **multi-valued** constraint kind on atoms (and bonds): an indene 5/6-ring junction atom legitimately satisfies *both* `RingSize(Lit(5))` and `RingSize(Lit(6))` simultaneously. Asserting both as separate constraints on the same atom is meaningful (conjunction by default) — the atom must be in some 5-ring AND some 6-ring.

This is the only multi-valued kind in the current `AtomConstraint` / `BondConstraint` surface; every other kind (Valence, Degree, Charge, RingCount, etc.) is single-valued per atom.

Storage adjustment (no Vec-shape restructure required):

```rust
impl AtomConstraint {
    fn is_unique(&self) -> bool {
        match self {
            AtomConstraint::RingSize(_) => false,
            _ => true,
        }
    }
}

impl AtomConstraints {
    pub fn add(&mut self, c: AtomConstraint) -> Option<AtomConstraint> {
        if c.is_unique() {
            // existing last-wins-per-kind path
        } else {
            // append within the kind cluster (entries stay grouped by kind)
        }
    }
}
```

Same applies to `BondConstraint::RingSize`. The constraint stores keep their sorted-by-kind invariant, generalized to sorted-and-clustered-by-kind. `AromaticSystemConstraints` already uses this pattern; the change unifies the implementations.

Store-side accessors:

- `view.constraints().ring_count() -> ValueAst` (single-valued)
- `view.constraints().ring_sizes() -> impl Iterator<Item = &ValueAst>` (plural; multi-valued)

The plural-vs-singular naming signals the data-model arity. Other single-valued kinds keep singular names returning `ValueAst` with Undetermined ≡ no constraint (principle 5).

Surface-form composition (DSL and EDN):

| Form | Meaning |
|---|---|
| DSL `[C#r5]` | atom in some 5-ring (one entry) |
| DSL `[C#r{5,6}]` | atom in some 5- OR 6-ring (`RingSize(LitSet([5, 6]))`, disjunction, one entry) |
| DSL `[C#r5#r6]` | atom in 5-ring AND 6-ring (two entries, conjunction) |
| EDN `{:ring-size 5}` | `RingSize(Lit(5))` (one entry) |
| EDN `{:ring-size [5 6]}` | `RingSize(LitSet([5, 6]))` (LitSet, disjunction, one entry — existing serialization, not conjunction) |
| EDN at molecule level: `[{:atom [0 {:ring-size 5}]} {:atom [0 {:ring-size 6}]}]` inside `:constraints` | two entries on atom 0, conjunction (the `:constraints` vec is a flat conjunction at top level) |
| EDN atom-inline via atom DSL string `"#r5#r6"` | conjunction (parser produces two entries inside the atom store) |

No new per-kind EDN sugar for conjunction. The vector-under-kind form `{:ring-size [5 6]}` is already taken by `LitSet` (disjunction) and stays that way. Conjunction emerges from:

- Repeated DSL predicates inside the atom string (`#r5#r6`),
- Multiple entries in the molecule-level `:constraints` vec (which is already a conjunction at top level).

Edit vocabulary implication: `SetAtomConstraint` works cleanly for unique kinds (replace by kind). For multi-valued kinds, we need explicit `Add` / `Remove` semantics — the user specifies which `RingSize(value)` is being added or removed, not just the kind. Detail design lives with the Edit-vocab impl.

Parametric view-side sugar (entity-centered queries against a caller-supplied `RingSet`). Suffix convention: `_from(&RingSet)` — reads as "drawn from this ring set":

| Operation | State | Notes |
|---|---|---|
| `atom_view.is_in_ring_from(&RingSet) -> bool` | **Open** | non-canonical boolean predicate |
| `atom_view.ring_count_from(&RingSet) -> ValueAst` | **Open** | count against caller-supplied set |
| `atom_view.ring_size_from(&RingSet) -> impl Iterator<Item = usize>` | **Open** | multi-valued: containing-ring sizes from the supplied set |
| `atom_view.rings_from(&RingSet) -> impl Iterator<Item = RingView<'_>>` | **Open** | iterator over rings from the supplied set containing this atom |
| Bond-side parametric analogs (`bond_view.is_in_ring_from(...)` etc.) | **Open** | symmetric |
| `aromatic_system_view.overlapping_rings_from(&RingSet)`, `.overlapping_atoms_from(&RingSet)`, `.overlapping_bonds_from(&RingSet)` | **Open** | non-canonical analog of the canonical-keyed overlap queries |

The view-side parametric methods are sugar — functionally they call the `RingSet`-side methods with the view's idx. The redundancy is intentional: it centers the object of interest at the call site (`atom_view.in_ring_set(&rs)` reads naturally when the caller already has the view in hand).

##### Settled: `RingView` context

`RingSet` carries a `&MoleculeAst` back-pointer; `RingView<'a>` carries a `&'a RingSet<'a>` back-pointer (and reaches the molecule transitively via `ring_view.set.molecule`). Smallest change that resolves the freestanding-RingView problem:

- Per-ring queries that need ring-vs-ring context (`is_spiro`, `is_fused`, `shared_atoms` between this ring and another) are methods on `RingView`, dispatching through `&RingSet`.
- Per-ring queries that need atom/bond data (degree of a ring atom, neighbors outside the ring) reach the molecule through `set.molecule`.
- `RingView` stays as the iteration item type for `ring_set.iter()` and the random-access return type for `ring_set.get(idx)` — keeps the per-ring API consistent with the rest of the View pattern.

#### Inter-entity derived reads

Properties computed across multiple entities. Two organizing principles:

**Constraint ↔ derived-quantity symmetry.** Every atom/bond-local constraint variant has a corresponding derived-quantity reader on the view, and the two should be name-paired (the reader name is the lower-cased / unprefixed constraint name). This enables direct read-vs-assert comparisons during validation and matcher dispatch. Per §"Procedural vs declarative split", multi-entity topological relations (ring-vs-ring spiro/fused) stay procedural with no constraint variant.

**Degree / valence variants.** Six aggregates over incident bonds, distinguished along two axes:

- *Order weighting*: count bonds equally (`degree` family) vs. sum bond orders (`valence` family).
- *Hydrogen inclusion*: exclude both explicit and implicit H (`heavy_*`), include explicit H only (the bare name, equals what's stored in the graph), or include both explicit and implicit H (`total_*`).

| Reader | Counts | Includes explicit H | Includes implicit H | Includes aromatic | Includes multicenter | SMARTS / DSL |
|---|---|---|---|---|---|---|
| `degree` | bonds, each = 1 | yes | no | no | no | `D` / `#D` |
| `total_degree` | bonds, each = 1; co-participants per multicenter | yes | yes | no (no new neighbors) | yes (`multicenter_degree`) | `X` / `#X` (extends SMARTS for multicenter — equal to SMARTS `X` where multicenter absent) |
| `heavy_atom_degree` | bonds, each = 1 | no | no | no | no | — |
| `valence` | localized bond orders | yes | no | no | no | — / `#v` |
| `total_valence` | electron-sharing contributions | yes | yes (each implicit H = 1) | yes (`aromatic_valence`) | yes (`multicenter_valence`) | — / `#V` |
| `heavy_atom_valence` | localized bond orders | no | no | no | no | — |

The aromatic-vs-multicenter asymmetry inside the *total* family is principled: aromatic systems overlay localized bonds (atoms in an aromatic system are already neighbors via localized bonds, so no new connections appear), while multicenter bonds, per the no-overlap rule (§"Terminology" → structural assumptions), connect atoms that *aren't* localized-bond neighbors. So multicenter contributes new "neighbor count" entries; aromatic doesn't.

`connectivity` is dropped as a reader name; the quantity is the same as `total_degree` and the systematic name fits the family. **DSL surface symbol `#X` is unchanged** — it still maps to the SMARTS `X<n>` semantics; the rename is purely internal (constraint variant `Connectivity → TotalDegree`, reader method `connectivity() → total_degree()`). DSL users see no change.

##### `total_valence` definition (full electron-sharing sum)

`total_valence` = `valence` + `implicit_hydrogens` + `aromatic_valence` + `multicenter_valence`. DSL symbol `#V`. This is the comprehensive electron-sharing contribution at the atom — every term is electrons this atom donates to a shared interaction. Excludes dative (`donated_pairs`, `accepted_pairs`) and non-covalent participations on purpose: those aren't electron-sharing in the same sense.

This corresponds to one column of the per-atom electron-accounting invariant from doc 52 §"Three valence types":

```
charge + lone_pairs + unpaired + total_valence + 2·donated + 2·accepted = outer_electrons
```

Diverges from SMARTS `v<n>` for aromatic atoms where the per-atom aromatic contribution is a donated lone pair (pyrrole N, furan O, thiophene S):

| Atom | Localized | Implicit H | Aromatic | Multicenter | `total_valence` | Textbook valence |
|---|---|---|---|---|---|---|
| C in benzene | 2 | 1 | 1 | 0 | **4** | 4 |
| N in pyridine | 2 | 0 | 1 | 0 | **3** | 3 |
| N in pyridinium / N-methylpyridinium | 3 | 0 | 1 | 0 | **4** | 4 |
| **N in pyrrole** | 2 | 1 | **2** | 0 | **5** | 3 |
| O in water | 0 | 2 | 0 | 0 | **2** | 2 |
| O in oxonium | 0 | 3 | 0 | 0 | **3** | 3 |
| O in furan | 2 | 0 | **2** | 0 | **4** | 2 |

The pyrrole/furan-style cases (lone-pair donors) give `total_valence` exceeding the textbook bond-count valence by exactly the donated pair (2 vs 0 for the aromatic term). This is not a bug — it reflects electron accounting precisely. Textbook valence buries the donated lone pair; our `total_valence` exposes it. Pattern matching wanting the bond-count form composes `valence + implicit_hydrogens + 1[if in aromatic system]`; `total_valence` is the electron-count form.

##### `total_degree` parallel definition

```
total_degree = degree + implicit_hydrogens + multicenter_degree
```

Where `multicenter_degree(a)` = Σ over multicenter bonds containing `a` of (atom count − 1) — the count of multicenter co-participants which, per the no-overlap rule, aren't already localized-bond neighbors.

Equals SMARTS `X<n>` for molecules without multicenter bonds. Extends SMARTS in cases SMARTS can't address (multicenter participation), introducing no incompatibility for SMARTS-compatible molecules.

Aromatic systems contribute to `total_valence` (electron donation) but *not* to `total_degree` (no new bonded neighbors — aromatic atoms are connected via localized bonds). This asymmetry inside the "total" family is principled, not arbitrary.

Heavy-atom variants have no SMARTS counterpart by design; constraint variants for them are not added preemptively. Per the "don't need to be limited by SMARTS" rule (memory: avoid software-development conventions when they don't fit) — constraint variants can be added later if a use case surfaces.

Per principle 5, value-bearing derived reads return `ValueAst` (collapsed to `Lit` or `Undetermined`). Per-kind constraint accessors on the entity stores return `ValueAst` (where "no constraint" is indistinguishable from `Undetermined`). Marker derived reads return `bool`.

Atom-side derived readers:

| Derived reader | State | Constraint counterpart | Returns | Notes |
|---|---|---|---|---|
| `atom_view.valence()` | **Impl** (rename from `bond_order_sum`) | `view.constraints().valence()` ↔ `Valence` | `ValueAst` | sum of incident `Bond.order`, no implicit H |
| `atom_view.total_valence()` | **Open** | `TotalValence` (DSL `#V`) | `ValueAst` | full electron-sharing sum: `valence + implicit_H + aromatic_valence + multicenter_valence`; diverges from SMARTS `v<n>` for aromatic lone-pair-donors (pyrrole N=5, furan O=4 vs textbook 3 and 2 — see §"Degree / valence variants") |
| `atom_view.heavy_atom_valence()` | **Designed** | — | `ValueAst` | bond-order sum, no implicit/explicit H |
| `atom_view.degree()` | **Impl** | `Degree` (SMARTS `D`) | `ValueAst` | incident bond count, no implicit H |
| `atom_view.total_degree()` | **Impl** (rename from `connectivity`) | `TotalDegree` (renamed from `Connectivity`) | `ValueAst` | `degree` + implicit H; DSL symbol `#X` unchanged (still maps to SMARTS `X<n>` semantics) |
| `atom_view.heavy_atom_degree()` | **Designed** | — | `ValueAst` | bond count, no implicit/explicit H |
| `atom_view.total_hydrogens()` | **Impl** | `TotalHydrogens` (SMARTS `H`) | `ValueAst` | implicit + explicit H |
| `atom_view.donated_pairs()`, `accepted_pairs()` | **Impl** | `DonatedPairs`, `AcceptedPairs` | `ValueAst` | dative-bond accounting |
| `atom_view.multicenter_degree()` | **Open** | — (reader-only for now) | `ValueAst` | sum, across all multicenter bonds this atom is in, of (co-participant count). Per the no-overlap rule, these aren't localized-bond neighbors. Building block for `total_degree`; constraint variant deferred until a consumer arrives |
| `atom_view.aromatic_valence()` | **Impl** (rename from `aromatic_contribution`) | `AromaticValence` | `AromaticValenceAst` | per-atom delocalized-electron count under the atom's aromatic system |
| `atom_view.multicenter_valence()` | **Impl** (rename from `multicenter_contribution`) | `MulticenterValence` | `MulticenterValenceAst` | symmetric for multicenter bonds |
| `atom_view.in_ring()` | **Open** | `InRing` (implicit via `RingCount >= 1`) | `bool` | canonical ring set (Vismara, max_size 22) |
| `atom_view.ring_count()` | **Open** | `RingCount(ValueAst)` | `ValueAst` | count of canonical containing rings |
| `atom_view.rings()` | **Open** | `RingSize(ValueAst)` matches against this | `Iterator<RingView>` | interpretation B for `RingSize`: constraint matches iff some ring in this iterator has size matching the pattern (see §"Ring access" → "Ring size constraint semantics") |
| `atom_view.smallest_ring_size()` | **Open** | — (utility, not paired) | `ValueAst` | smallest containing ring size; chemistry-classification helper |
| `atom_view.ring_degree()` | **Open** | `RingDegree(ValueAst)` (rename from `RingConnectivity`; SMARTS `x<n>`) | `ValueAst` | count of incident ring bonds, each = 1 |
| `atom_view.ring_valence()` | **Open** | `RingValence(ValueAst)` (new; DSL `#y`) | `ValueAst` | sum of bond orders of incident ring bonds; no SMARTS analog |

Bond-side derived readers:

| Derived reader | State | Constraint counterpart | Returns | Notes |
|---|---|---|---|---|
| `bond_view.endpoints() -> [AtomId; 2]` | **Open** | — | array | sugar over `Graph::edge_endpoints` |
| `bond_view.atoms() -> impl Iterator<Item = AtomView<'_>>` | **Open** | — | iter | yield endpoint atoms as views |
| `bond_view.in_ring()`, `smallest_ring_size()`, `ring_count()` | **Open** | `InRing` / `RingSize` / `RingCount` | `bool` / `ValueAst` / `ValueAst` | canonical-rings sugar |
| `bond_view.is_in_aromatic_system()` | **Open** | — | `bool` | symmetric with `atom_view.is_in_aromatic_system()` |

Relation-view derived readers:

| Derived reader | State | Returns | Notes |
|---|---|---|---|
| `aromatic_system_view.electron_count()` | **Open** | `ValueAst` | sum over `electrons[]`; `Undetermined` if any participant electron is non-Lit; pairs with `AromaticElectronCount` constraint |
| `aromatic_system_view.atom_count()`, `bond_count()` | **Open** | `usize` | participant counts |
| `aromatic_system_view.overlapping_atoms(&[AtomId])`, `overlapping_bonds(&[BondId])` | **Open** | iter of views | intersection with caller-supplied set; P1 (views); `_ids` companions available |
| `aromatic_system_view.overlapping_rings()` | **Open** | `Iterator<RingView>` | canonical-rings sugar (§"Ring access") |
| `multicenter_bond_view.electron_count()` | **Open** | `ValueAst` | parallel; pairs with `MulticenterElectronCount` |
| `multicenter_bond_view.atom_count()`, `bond_count()` | **Open** | `usize` | participant counts |
| `multicenter_bond_view.overlapping_atoms(&[AtomId])` | **Open** | `Iterator<AtomView>` | overlapping bonds and rings *provisionally* not added — multicenter bonds typically don't share bonds with other relations; revisit if a use case surfaces |
| `dative_bond_view.atom_count()` | **Open** | `usize` | not always 2 — multi-donor dative bonds are real |
| `dative_bond_view.overlapping_atoms(&[AtomId])`, `overlapping_bonds(&[BondId])`, `overlapping_rings()` | **Open** | iter of views | full overlap surface (think borazine: dative bonds participating in aromatic ring systems) |
| `noncovalent_bond_view.*` | not added | — | wait for a consumer |

Pairwise ring-vs-ring and molecule-level entries:

| Reader | State | Constraint counterpart | Returns | Notes |
|---|---|---|---|---|
| `RingSet::relation`, `are_spiro/fused/bridged`, `*_neighbors`, `fused_component(s)`, `shared_atoms`, `shared_bonds` | **Impl** | **none** (procedural-only by design) | — | aromaticity-perception internal; no SMARTS analog — see §"Procedural vs declarative split" |
| Total charge / spin | **Designed (constraint propagator stubs)** | `TotalCharge`, `TotalSpin` | `ValueAst` / `SpinStateAst` | propagators stubbed in `ConstraintValidator` per §"Outstanding" |
| `bonds().connecting(a: AtomId, b: AtomId) -> Option<BondView<'_>>` | **Impl** | — | `Option<BondView>` | atoms-in, bond-out; lives on `BondViews`, not `MoleculeAst`. `_id` companion: `bonds().connecting_id(a, b) -> Option<BondId>` |
| `bonds().induced(atoms: &[AtomId]) -> Vec<BondView<'_>>` | **Impl** | — | `Vec<BondView>` | both endpoints in subset; lives on `BondViews`. `_ids` companion: `bonds().induced_ids(atoms) -> Vec<BondId>` |

##### Comparison idiom

For value-bearing constraints (the common case):

```rust
view.constraints().valence().matches(&view.valence())
```

Always well-typed; `Undetermined` on either side is wildcard-matched. No `Option` unwrapping, no `_constraint` suffix on the view, no asymmetric numeric vs. AST types.

For marker constraints (`Aromatic`, `InRing`):

```rust
view.constraints().aromatic() == view.aromatic()  // both bool
view.in_ring() == view.constraints().in_ring()    // both bool
```

#### Constraint reads

Per doc 87 the constraint architecture has two scopes: per-entity inline (`*Constraints` stores attached to each entity AST) and molecule-list (a flat `Vec<Constraint>` for combinators, relational, and molecule-scope predicates). The two are read very differently:

##### Per-entity inline (`AtomConstraints` etc.)

Each entity-kind constraint variant is single-valued per entity (the store enforces this via `add()`'s last-wins-per-kind). The natural shape is per-kind named accessors:

| Operation | State | Notes |
|---|---|---|
| `AtomConstraints::valence() -> ValueAst` | **Open** | single-valued; `Undetermined` if no constraint or if asserted as `Undetermined` (principle 5) |
| `AtomConstraints::degree()`, `total_degree()`, `heavy_atom_degree()`, `total_hydrogens()`, `donated_pairs()`, `accepted_pairs()`, `ring_count()`, `ring_degree()`, `ring_valence()` | **Open** | single-valued; all return `ValueAst` |
| `AtomConstraints::ring_sizes() -> impl Iterator<Item = &ValueAst>` | **Open** | **multi-valued** (per §"Ring size as a multi-valued constraint kind"); plural name signals arity |
| `AtomConstraints::aromatic_valence()` | **Open** | returns `AromaticValenceAst` (its own Undetermined-bearing AST) |
| `AtomConstraints::multicenter_valence()` | **Open** | returns `MulticenterValenceAst` |
| `BondConstraints::aromatic() -> bool` | **Open** | marker; true iff asserted |
| `BondConstraints::ring_count() -> ValueAst` | **Open** | single-valued |
| `BondConstraints::ring_sizes() -> impl Iterator<Item = &ValueAst>` | **Open** | multi-valued; same rationale as atom side |
| `DativeBondConstraints::*`, `AromaticSystemConstraints::*`, `MulticenterBondConstraints::*` | **Open** | symmetric, one accessor per kind |
| `*Constraints::contains(kind)`, `get(kind)`, `add(c)`, `remove(kind)`, `retain`, `clear`, `iter` | **Impl** | generic store API; stays for code that walks all entries |
| `*Constraints::get_all(kind) -> impl Iterator<Item = &Constraint>`, `remove_all(kind) -> Vec<Constraint>` | **Open** | multi-valued generic accessors; new for the `is_unique`-aware add path |

Per-kind singular accessors are sugar over `get(kind)`, with Undetermined-unification baked in. Plural accessors walk the kind cluster directly. `is_unique` predicate on each constraint variant (`AtomConstraint`, `BondConstraint`) drives `add()` behavior — same shape as `AromaticSystemConstraints`.

##### Molecule-list (`mol.constraints()`)

A `Vec<Constraint>` where the same constraint kind can appear multiple times with different idx targets, plus combinators (`And` / `Or` / `Not`), `Relational`, and `Molecule`-scope predicates. Per-kind named accessors don't fit (multi-valued, indexed by entity reference). The natural shape is iterator-filtering:

| Operation | State | Notes |
|---|---|---|
| `mol.constraints() -> &Constraints` | **Impl** | molecule-list read |
| `Constraints::iter()`, `len()`, `as_slice()` | **Impl** | flat-vec walks |
| `Constraints::for_atom(AtomId)`, `for_bond(BondId)`, etc. | **Open** | filter by entity reference; recurring pattern in `ops/validator/invariant.rs`; per naming conventions §"Suffix conventions" |
| `Constraints::by_kind(ConstraintKind)` | **Open** | filter by variant discriminant; same recurring pattern |

Naming convention (working): `for_*` for "constraints referencing this entity", `by_*` for "constraints matching this discriminant". Both return iterators (consistent with §"Cross-entity navigation" working convention — iterators primary, no dedicated `count_*` / `has_*` methods). The two filters compose: `constraints.for_atom(idx).by_kind(AtomConstraintKind::Valence)`.

#### Molecule-scope state predicates

Boolean queries about the molecule as a whole. Each is a one-liner over an existing iterator / counter, but recurring enough in ops code (and in the future SMILES / MOL writers) to be worth named methods.

| Operation | State | Notes |
|---|---|---|
| `mol.is_ground() -> bool` | **Impl** | all entity attributes are concrete (`molecule.rs:525-528`); gates the resolver-output assumption |
| `mol.is_empty() -> bool` | **Open** | zero atoms; sugar over `mol.atoms().count() == 0` |
| `mol.has_constraints() -> bool` | **Open** | true if any per-entity inline `*Constraints` is non-empty OR `mol.constraints()` is non-empty |
| `mol.has_overlays() -> bool` | **Open** | umbrella: true if any of the four overlays is non-empty. Useful as a "is this molecule topology-only?" check — false means the molecule is a pure atom + localized-bond skeleton, no aromatic systems / multicenter / dative / noncovalent enrichment |
| `mol.has_dative_bonds() -> bool` | **Open** | sugar over `mol.dative_bonds().count() > 0` |
| `mol.has_aromatic_systems() -> bool` | **Open** | sugar over `mol.aromatic_systems().count() > 0` |
| `mol.has_multicenter_bonds() -> bool` | **Open** | sugar over `mol.multicenter_bonds().count() > 0` |
| `mol.has_noncovalent_bonds() -> bool` | **Open** | sugar over `mol.noncovalent_bonds().count() > 0` |

#### Metadata: lives on `MoleculeDsl`, not on `MoleculeAst`

Entity ids (`atom0`, `bond1`, ...) and atom aliases (`al0` → `AtomAst`) are surface-form metadata carried by `MoleculeDsl`, not by `MoleculeAst`. The split is deliberate (per doc 94 / IO ergonomics): `MoleculeAst` carries algebraic content only; `MoleculeDsl` carries the boundary representation with metadata. Reads of metadata happen on the DSL side via `MoleculeDsl::metadata() -> &Metadata` and `Metadata::iter_atom_aliases() / iter_*_ids()`. There is no metadata reader on `MoleculeAst` by design.

### Mutation operations

Per principle 6, the load-bearing checked mutation primitive is `MoleculeBuilder::transact(Vec<Edit>) -> Result<Transaction, TransactionError>`. The tables in the subsections below describe the operations *semantically*; in the implementation, user-facing mutation sugar routes through the `Edit` vocabulary where atomicity or rollback matters. The detailed revised transaction design lives in §"Mutation API revision: undo-journal transactions".

#### Edit vocabulary and transactions

This subsection is the semantic summary. The older self-inverting-`Edit` / `Action` / snapshot design was superseded during the Phase-8 revision.

Current shape:

- `Edit` is caller-facing and composable.
- Bulk topology variants are primitive: `AddAtoms`, `AddBonds`, `RemoveTopology`.
- Single atom/bond additions and removals are sugar over the bulk variants.
- Overlay relation edits remain single-item for now.
- `Undo` is the realized physical rollback vocabulary.
- `Transaction` wraps `Vec<Undo>` and rollback consumes it.
- Checked `transact` records `Undo` and reverse-replays on failure.
- `transact_unchecked` is a separate no-journal, panic-on-invalid path for known-correct generated edits.

##### Symbolic refs (intra-batch dependency)

A later edit in a batch may depend on an entity created by an earlier edit. The Ref types lift the absolute-index requirement only inside `Edit`:

```rust
enum AtomRef            { Id(AtomId),            New(usize) }
enum BondRef            { Id(BondId),            New(usize) }
enum DativeBondRef      { Id(DativeBondId),      New(usize) }
enum AromaticSystemRef  { Id(AromaticSystemId),  New(usize) }
enum MulticenterBondRef { Id(MulticenterBondId), New(usize) }
enum NoncovalentBondRef { Id(NoncovalentBondId), New(usize) }
```

Resolution at apply time: `Id(_)` is used directly; `New(N)` resolves against an applicator-maintained created-entity table. Type mismatch (e.g., `AtomRef::New(N)` where created entity N is a bond, not an atom) fails checked `transact` and panics in `transact_unchecked`.

Refs appear only inside `Edit`. `Undo`, view methods, remapping methods, and restore payloads use absolute `*Id` types.

##### Transaction API

```rust
impl MoleculeBuilder {
    pub fn transact(&mut self, edits: Vec<Edit>)
        -> Result<Transaction, TransactionError>;

    pub fn transact_unchecked(&mut self, edits: Vec<Edit>);
}

impl Transaction {
    pub fn rollback(self, builder: &mut MoleculeBuilder) -> Result<(), TransactionError>;
    pub fn undos(&self) -> &[Undo];
}
```

Semantics:

- **Atomic**: success → all edits applied; failure → state restored to pre-transaction.
- **Validation timing**: tier-1 invariants checked on commit; intermediate states inside a transaction may transiently violate them.
- **Cascade integration**: `RemoveTopology` drops invalid overlays; the realized `Undo` records the dropped overlays for rollback.
- **Rollback invariant**: reverse replay of each `Undo` restores the exact id coordinate system that existed before the corresponding forward edit.

##### Settled design choices

- **Field-enum `Set*` variants.** Single outer Edit variant per entity (`SetAtomField`, `SetBondField`, ...), wrapping a per-entity `*FieldChange` enum that carries the field discriminant + type-correct `old`/`new`:

  ```rust
  enum AtomFieldChange {
      Element            { old: ElementAst,            new: ElementAst },
      IsotopeMass        { old: IsotopeAst,            new: IsotopeAst },
      Charge             { old: ValueAst,              new: ValueAst },
      ImplicitHydrogens  { old: ImplicitHydrogensAst,  new: ImplicitHydrogensAst },
      LonePairs          { old: ValueAst,              new: ValueAst },
      Spin               { old: SpinStateAst,          new: SpinStateAst },
  }

  enum Edit {
      SetAtomField { idx: AtomRef, change: AtomFieldChange },
      // ...
  }
  ```

  Outer `Edit` stays manageable; per-field type discrimination preserved inside `*FieldChange`. Symmetric `BondFieldChange`, `DativeBondFieldChange`, `AromaticSystemFieldChange`, `MulticenterBondFieldChange`, `NoncovalentBondFieldChange` for the other entity kinds.

- **Rollback mechanism: realized `Undo` journal.** Checked `transact` records the payloads, remaps, cascades, and constraint updates needed to restore dense storage. It does not snapshot the whole molecule.

- **Pre-commit validation hook.** `tx.validate() -> Result<(), Vec<Error>>` is callable at any point during a transaction; commit calls it implicitly before applying. Lets callers do "try this, check, decide whether to commit" workflows without forcing commit-only validation.

##### Cost model and batching guidance

For checked `transact`, per-transaction overhead breaks down as `N × apply + N × undo-capture + 1 × commit-validation` for an N-Edit transaction. There is no whole-molecule snapshot. Undo-capture cost is proportional to the actual changed entities and any constraints affected by remapping.

This means **per-Edit overhead amortizes to almost nothing as long as callers batch into transactions**. Hot-path callers (resolver, propagator loops) should aggregate mutations into one transaction per natural unit of work (one fixpoint iteration, one propagator pass) rather than wrapping each mutation in its own transaction. The natural batch boundary is whatever the caller's algorithm uses anyway.

Pattern to avoid: `for narrowing in narrowings { mol.transact(vec![narrowing])? }` — pays N validations and builds N transactions.
Pattern to use: `mol.transact(narrowings.collect())?` — pays 1 validation and returns 1 transaction.

##### Unchecked escape hatch

`MoleculeBuilder::transact_unchecked(Vec<Edit>) -> ()` is a documented escape hatch that skips undo capture and rollback. Intended for callers that produce known-correct batches by construction:

- Reaction rules whose RHS is verified at rule-definition time.
- Built-in transformations (`kekulize`, `aromatize`) that compute the result and prove its correctness internally.
- Test fixtures and replay paths from logged Edit sequences.

The unchecked path is a separate implementation, not `transact(edits).unwrap()`. It constructs no `Undo`; invalid refs or failed preconditions panic. The contract: caller asserts the edit stream is correct and the resulting molecule is tier-1-valid.

Not the load-bearing primitive — checked `transact` is the default. `transact_unchecked` exists for the small slice of callers where correctness can be argued at the rule / transform level and undo allocation is avoidable.

#### Entity feature attribute mutators

In-place mutation of entity feature data without changing topology.

| Operation | State | Notes |
|---|---|---|
| `mol.atom_mut(idx) -> AtomViewMut` | **Impl** | direct in-place edit on the AST |
| `mol.bond_mut(idx) -> BondViewMut` | **Impl** | |
| `mol.dative_bond_mut(idx)`, `aromatic_system_mut`, `multicenter_bond_mut`, `noncovalent_bond_mut` | **Impl** | symmetric |
| `MoleculeBuilder::atom_mut(idx)` (Phase-2 of DPO) | **Designed** (doc 90 "Missing on MoleculeBuilder") | for K-entity attribute mutation during a structural edit; not yet on builder |
| `MoleculeBuilder::bond_mut(idx)`, etc. | **Designed** | same |

#### Adjacency mutators

Topology edits via `MoleculeBuilder` with copy-on-write of the Arc-shared CSR. Changes return an `IdRemapping` so callers can reindex external data.

| Operation | State | Notes |
|---|---|---|
| `MoleculeBuilder::add_atom() -> AtomId` | **Impl** | |
| `MoleculeBuilder::add_bond(a, b, BondAst) -> BondId` | **Impl** | |
| `MoleculeBuilder::remove(&[AtomId], &[BondId]) -> IdRemapping` | **Impl** | low-level dense removal primitive; transaction rollback captures the richer removed payloads and `ConstraintUpdate` before calling it |
| `MoleculeBuilder::remove_atom(AtomId)`, `remove_bond(BondId)` | **Impl** | sugar over `remove` |

##### Cascade semantics on `MoleculeBuilder::remove`

**R3 — cascade-and-report.** `remove` deletes the requested atoms and bonds, then drops any overlay that becomes invalid as a result, returning the dropped overlays alongside the index remapping. The caller is responsible for restoring chemistry validity (re-perception, re-adding compensating overlays, charge balancing) — the AST does not attempt to preserve overlays automatically because:

- A "shrunk" aromatic system (atom removed, system kept with fewer atoms) carries no useful semantic — the original perception assumed all atoms; the partial set isn't an aromatic system at the new shape.
- Silent re-perception would conflate caller-asserted aromaticity with algorithm-detected aromaticity, losing the user's intent.
- "Refuse" would force every reaction / DPO call to do an upfront overlay-cleanup pass, which is friction without benefit when the caller knows what they're doing.

The chosen middle ground (R3 from the cascade-options table): execute the cascade but record what was dropped, so checked transactions can roll it back and callers inspecting `Undo` can see what changed.

What gets dropped:

| Mutation | Overlays dropped |
|---|---|
| Atom removed | every overlay whose participant list contains that atom (all four overlay types) |
| Bond removed | every aromatic system whose `bonds()` includes that bond — *only* aromatic systems, since multicenter/noncovalent overlays don't share localized bonds (provisional structural rule, §"Terminology"), and dative bonds coexist independently of localized bonds |

Return type: cascade is reported through the realized `Undo` for `RemoveTopology`. The undo entry stores the dropped overlay payloads, the forward `IdRemapping`, the `UndoRemapping`, and the `ConstraintUpdate` needed to restore the previous dense coordinate system.

Index remapping (old→new atom/bond indices after the removal) remains `IdRemapping`, primarily for downstream consumers that hold ids.

##### Pre-mutation predicate

To support callers that prefer to clean up overlays explicitly before calling `remove` (rather than relying on the cascade-and-report path), the atom view exposes a binary check:

| Operation | State | Notes |
|---|---|---|
| `atom_view.is_in_overlays() -> bool` | **Open** | umbrella; atom-scope mirror of `mol.has_overlays()`. True if the atom participates in any of the four overlays; backed by the existing per-overlay-kind iterators in §"Cross-entity navigation" |

#### Whole-relation mutators

Add or remove an entire dative-bond / aromatic-system / multicenter-bond / noncovalent-bond entry on the builder.

| Operation | State | Notes |
|---|---|---|
| `add_dative_bond(donor, acceptor, DativeBondAst) -> DativeBondId` | **Impl** | doc 90 |
| `add_aromatic_system(atoms, AromaticSystemAst) -> AromaticSystemId` | **Impl** | |
| `add_multicenter_bond(atoms, MulticenterBondAst) -> MulticenterBondId` | **Impl** | |
| `add_noncovalent_bond(ends, NoncovalentBondAst) -> NoncovalentBondId` | **Impl** | |
| `remove_dative_bonds(&[DativeBondId])` | **Designed** (doc 90) | wired into builder per doc 92 but DPO Phase-3 integration deferred |
| `remove_aromatic_systems(&[AromaticSystemId])` | **Designed** | same status |
| `remove_multicenter_bonds(&[MulticenterBondId])` | **Designed** | same status |
| `remove_noncovalent_bonds(&[NoncovalentBondId])` | **Designed** | same status |

#### Internal relation mutators

**Resolved by Edit + transactions.** Internal mutation of an overlay (add/remove an atom *within* an existing aromatic system, multicenter bond, etc., keeping `electrons` arrays in sync) is expressed as a transaction containing `Edit::RemoveAromaticSystem { ... }` followed by `Edit::AddAromaticSystem { ... }` with the desired new shape. The transaction is atomic, so the molecule never observes the in-between state. No dedicated `aromatic_system_add_atom` / `aromatic_system_remove_atom` API needed.

This matches doc 90's DPO assumption (relations are immutable participant tuples; modify by remove-and-readd) and pays the same conceptual cost as elsewhere — every modification is an Edit batch, not an in-place tweak. The cascade question that previously sat under this heading (when an atom is removed from the molecule, what happens to overlays it participates in) is settled separately as R3 under §"Cascade semantics on `MoleculeBuilder::remove`".

The convenience wrapper for this pattern (a single method that takes an overlay-idx and a desired-new-shape, expanding to the Remove + Add transaction internally) can be added later if usage shows it earns its keep; not part of the core API.

#### Constraint mutators

| Operation | State | Notes |
|---|---|---|
| Entity-inline: `view.data.constraints.add/remove/retain/clear` | **Impl** | |
| Molecule-list: `mol.constraints_mut().push/retain/clear` | **Impl** | |
| `mol.lift_constraints()` | **Impl** | inline → molecule-list (per-entity store drained) |
| `mol.inline_constraints()` | **Impl** | molecule-list → inline (top-level narrow leaves only; combinators / relational / molecule-scope preserved) |
| `MoleculeBuilder::constraints_mut()` | **Designed** | builder constraint surface thinner than AST; doc 90 §"Constraint access" |

### Cross-cutting

#### Subgraph and projection

| Operation | State | Notes |
|---|---|---|
| `induced_subgraph(&[AtomId]) -> MoleculeEmbedding<'_>` | **Impl** | `MoleculeEmbedding` borrows `&MoleculeAst` and carries per-entity-type local→parent index maps plus a parent→local inverse for atoms. `embedding.extract()` materializes the standalone sub-`MoleculeAst` on demand; `embedding.edits()` derives the `Vec<Edit>` (a single `RemoveTopology` over the complement). Same type will be produced by subgraph isomorphism matching when that lands. Parallel `umol_graph_core::Embedding<'_>` does the topology-only version |
| `bonds().induced(&[AtomId]) -> Vec<BondView<'_>>` | **Impl** | on `BondViews`; lighter alternative when only the bond list is needed; used by `AromaticSystemView` and `ops/transformer/kekulizer.rs`. `_ids` companion `bonds().induced_ids(atoms) -> Vec<BondId>` for raw-index access |
| Dense atom-ordering projection (atom → 0..n_subset map) | **Impl** | `MoleculeEmbedding::local_atom(parent_atom) -> Option<AtomId>` (O(1)); HMO's hand-rolled `atom_to_idx` HashMap removed in favor of this |

#### Index-space conversions

| Operation | State | Notes |
|---|---|---|
| `IdRemapping` (returned by `MoleculeBuilder::remove`) | **Impl** | old → new for all six relation kinds |
| Subgraph remapping (returned by `induced_subgraph`) | **Impl** | `MoleculeEmbedding<'a>` borrows `&MoleculeAst` and carries local→parent maps for all six entity kinds plus a parent→local inverse for atoms; same type intended for future subgraph isomorphism results |
| Constraint remapping under structural edit | **Impl** | piggy-backs on `IdRemapping` |

#### Pattern-vs-ground scope per operation

Some operations are meaningful only on ground molecules (electron-accounting reads, tier-2 invariant checks); others are pattern-meaningful (substructure search, constraint propagation). The taxonomy above is shared between scopes — the asymmetry is in *invariants*, not API shape:

- All read operations (entity, adjacency, relation, derived, constraint) work on both `MoleculeAst` ground and partial.
- All mutation operations work on both.
- Tier-2 invariant verification (`ElectronInvariantValidator`, `ConstraintValidator`) is ground-only.
- Cache-backed reads (`mol.rings(...)`) work on both, but pattern usage rarely benefits since pattern molecules are small and short-lived.

See doc 86 §"Asymmetries worth naming" for the framing of which operation classes flow only into `Molecule` (resolved-only) vs. both.

#### Iteration order contracts

The implicit-but-load-bearing guarantees today:

| Iterator | Order contract |
|---|---|
| `mol.atoms().iter()` | by `AtomId` ascending |
| `mol.bonds().iter()` | by `BondId` ascending |
| `mol.neighbors(atom)` | by `NodeId` ascending (CSR neighbors are sorted; relied on by `connecting_bond` binary-search and `induced_bonds` set-intersection) |
| `mol.dative_bonds().iter()` etc. (relations) | by relation index ascending |
| `relation_view.atoms()` | sorted ascending by `AtomId` (set semantics) |
| `metadata.iter_atom_aliases()` | sorted by alias name (`BiBTreeMap` iteration); see doc 95 |
| `mol.constraints().iter()` | insertion order (flat `Vec`); not sorted |
| `RingSet::iter()` | enumeration order (depends on ring-finding algorithm; not user-specified) |

These should be documented on the iterator-returning methods themselves. Anywhere an op constructs a `BTreeMap` or `sort_by_key` on AST output is a candidate for upgrading the iterator's order contract.

### What stays out (intentionally not on the AST)

Operations the AST does not surface; they live on higher tiers (`Molecule`) or in separate ops modules.

| Operation | Lives in | Notes |
|---|---|---|
| Aromaticity perception | `umol-graph/src/ops/aromaticity*` | discovery is procedural (doc 80 line 194); the AST receives the resulting aromatic-system tuples via `add_aromatic_system` |
| Ring enumeration | `RingEnumerator` (`umol-graph/src/ast/rings.rs`) | procedural; canonical output cached on `MoleculeAst` via `mol.rings()` (single slot); non-canonical / filtered via `mol.rings_with(...)` (uncached, owned) |
| Substructure matching (VF2) | matcher (deleted at end of doc 92, awaits port; doc 80 step 9) | operates on `Pattern` against `MatchTarget`; not a `MoleculeAst` method |
| Transformations: `kekulize`, `aromatize`, `tautomers`, `apply_reaction` | `umol-graph/src/ops/transformer*` | take `&Molecule`, return `Molecule` or `Vec<Molecule>`; see §"Transformation ops" |
| Validation: `ElectronInvariantValidator`, `ConstraintValidator`, `EntityStructureValidator`, `SpinCouplingValidator` | `umol-graph/src/ops/validator*` | take `&MoleculeAst` (`AsRef`); not methods on the AST |
| Resolution: `Solver`, `ChemistryModel` | `umol-graph/src/ops/resolver*` | top-level; `Solver::resolve(ast, model) -> Result<Molecule, _>` |
| Canonical SMILES, Morgan fingerprint | `Molecule` methods | tier-1 chemist surface, not AST |

### Summary of gaps surfaced by this taxonomy

The Open rows above represent real recurring needs from the ops survey or doc 90 that have not been designed. Concretely, before the next round of MoleculeAst API work:

1. **Cross-entity reverse navigation** (atom → containing relations) — recurring scan in valence and aromaticity ops; modest design.
2. **Adjacency-convenience** — `bonds().connecting(a, b)` and `bonds().induced(atoms)` settled and impl on `BondViews`; `GraphView::connected_components(alg)` impl; `connected_components_in(atoms)` pending.
3. **`induced_subgraph` shape** — settled as `MoleculeEmbedding<'_>` (borrow + index maps); impl on `MoleculeAst`.
4. **Internal relation mutators** — settled. Express as a `Remove*` + `Add*` Edit pair inside a transaction; no dedicated in-place API. See §"Internal relation mutators".
5. **Cascade behavior of topology removal** — settled as R3 (cascade-and-record). `RemoveTopology` drops invalid overlays; realized `Undo` records the dropped overlay payloads, remaps, and constraint updates needed for rollback. `atom_view.is_in_overlays()` umbrella predicate to land at the same time.
6. **`ring_count_at`** + `RingCount(ValueAst)` constraint — designed, not impl.
7. **Iteration-order contracts** — implicit; should be documented on each iterator-returning method.
8. **Constraint find-by-kind / find-by-idx queries** — recurring `for c in mol.constraints().iter() { match c { ... } }` pattern.
9. **`RingView` context** — settled. `RingSet` gets `&MoleculeAst` back-pointer; `RingView` gets `&RingSet` back-pointer; molecule reachable transitively via `view.set.molecule`. See §"Ring access" → "Settled: RingView context".
10. **Canonical ring view-side sugar** — `atom_view.in_ring()`, `smallest_ring_size()`, `ring_count()`, `rings() -> Iterator<RingView>`, plus bond-side and aromatic-system-overlap analogs, plus their parametric `*_in(&RingSet)` versions. Designed in §"Ring access"; not implemented.
11. **`RingFamily::Induced` → `RingFamily::Relevant`** rename — the existing enum variant for Vismara relevant cycles is misnamed. Mechanical, but on the naming-pass list.
12. **Per-field accessors on entity views** — `#[inline]` getters returning `&T` (or value where naturally `Copy`) for every field on each entity AST. ~25 methods; mechanical. See §"Entity feature reads" → "Per-field accessor pattern".
13. **Per-kind accessors on entity inline `*Constraints` stores** — `AtomConstraints::valence() -> ValueAst`, etc., plus markers as `bool`. Per principle 5: `Undetermined` ≡ "no constraint asserted"; single-valued per-kind accessor returns `ValueAst` with no `Option` wrapping. Multi-valued kinds (currently just `RingSize` on atom and bond) use plural-named accessors returning iterators (`ring_sizes() -> Iterator<&ValueAst>`); `is_unique` predicate on the constraint variant drives `add()` behavior. See §"Constraint reads" and §"Ring size as a multi-valued constraint kind".
14. **`ValueAst` arithmetic surface** — `Add` / `Sub` / `Mul` / `Div` (binary + scalar) with Lit/Undetermined collapse semantics per principle 5; panic on div-by-zero and overflow; negative results allowed. Same scheme extends to `ImplicitHydrogensAst` (with `Normal` collapsing to `Undetermined`). `IsotopeAst` excluded — isotope masses are enum-like, not numeric. Mechanical; impl-only.
15. **Naming alignment renames on derived readers and constraints** — workspace-wide `*Idx → *Id` rename across `umol-ast`, `umol-graph`, ops/, tests: `AtomIdx → AtomId`, `BondIdx → BondId`, `DativeBondIdx → DativeBondId`, `AromaticSystemIdx → AromaticSystemId`, `MulticenterBondIdx → MulticenterBondId`, `NoncovalentBondIdx → NoncovalentBondId`, `RingIdx → RingId`; `IdxRemapping → IdRemapping`; `AtomRef::Idx → AtomRef::Id` and symmetric for the other Ref enums. Aligns molecule-scope identifiers with the `umol-graph-core` `NodeId`/`EdgeId`/`RelationId` convention; method names `.ids()`, `_id`, `_ids` become coherent with the type names. Plus reader/constraint renames: `bond_order_sum → valence`, `aromatic_contribution → aromatic_valence`, `multicenter_contribution → multicenter_valence`, `connectivity → total_degree` (reader) and `Connectivity → TotalDegree` (constraint variant; **DSL symbol `#X` unchanged**, still maps to SMARTS `X<n>` semantics — rename is internal only, DSL users see no change), `RingConnectivity → RingDegree` (constraint variant; SMARTS `x<n>` reads unweighted per Daylight spec). Pairs each derived reader with its constraint counterpart; on the naming-pass list with the field-rename items. Plus: `#R!` DSL sugar for `RingCount(Lit(0))` (atom and bond sides; symmetric with `#R+` for `RingCount(Expr("r >= 1"))`; canonical output for `Lit(0)`); `#R0` accepted on input. Parser + renderer changes localized to `umol-ast/src/dsl/predicates.rs` (`ring_count` and `fmt_ring_count`). New `RingValence` constraint variant (order-weighted ring bonds, no SMARTS analog; DSL symbol `#y`) with view-side reader `atom_view.ring_valence()`. New `TotalValence` constraint variant + `atom_view.total_valence()` reader (full electron-sharing sum `valence + implicit_H + aromatic_valence + multicenter_valence`; DSL symbol `#V`; replaces earlier "`valence + implicit H` only" definition; diverges from SMARTS `v<n>` for aromatic lone-pair-donors — see §"Degree / valence variants"). New `atom_view.multicenter_degree()` reader (reader-only for now, no constraint, no DSL; sum across multicenter bonds of co-participant count) and parallel `total_degree` definition extension: `total_degree = degree + implicit_H + multicenter_degree` — equals SMARTS `X<n>` for molecules without multicenter bonds, extends for those with.
16. **Relation-view derived readers** — `aromatic_system_view.{electron_count, atom_count, bond_count, overlapping_atoms, overlapping_bonds, overlapping_rings}`, `multicenter_bond_view.{electron_count, atom_count, bond_count, overlapping_atoms}` (overlapping_bonds and overlapping_rings provisionally not added — multicenter bonds don't typically share with other relations), `dative_bond_view.{atom_count, overlapping_atoms, overlapping_bonds, overlapping_rings}` (full overlap surface — borazine case). Noncovalent-bond view derived readers not added preemptively. See §"Inter-entity derived reads".
17. **`NeighborView` reshape** — all-private-fields shape: `{ atom_idx, bond_idx, molecule }` with `nbr.atom() -> AtomView` and `nbr.bond() -> BondView` accessor methods. Indices accessible via `nbr.atom().idx`. View construction lazy; cached `data: &BondAst` field dropped (callers go through `.bond()` accessor).
18. **Molecule-scope state predicates** — `is_empty`, `has_constraints`, `has_overlays` (topology-only check), plus per-overlay `has_dative_bonds` / `has_aromatic_systems` / `has_multicenter_bonds` / `has_noncovalent_bonds`. `is_ground` already impl. See §"Molecule-scope state predicates".
19. **Terminology**: "overlays" adopted as the umbrella term for the four typed n-ary relations (`DativeBond`, `AromaticSystem`, `MulticenterBond`, `NoncovalentBond`). Codebase-internal usage; documented in §"Terminology".
20. **Edit vocabulary** (principle 6) — `Edit` enum with bulk topology primitives (`AddAtoms`, `AddBonds`, `RemoveTopology`) plus single-item sugar; overlay edits remain single-item unless a concrete bulk need appears. `AtomRef`/`BondRef`/four overlay Refs with `Id(_)` and `New(usize)` variants remain the intra-batch dependency mechanism. See §"Mutation API revision: undo-journal transactions".
21. **Transaction primitive** — `MoleculeBuilder::transact(Vec<Edit>) -> Result<Transaction, TransactionError>` is the checked, atomic, undo-journaled path. It returns a `Transaction` wrapper over realized `Undo` entries; rollback consumes the wrapper and reverse-replays the journal. `transact_unchecked(Vec<Edit>) -> ()` is a separate implementation for known-correct generated edits: no `Undo` allocation, no rollback capability, panic on invalid input.
22. **Field-typed `Set*` Edit variants** — settled as field-change enums (`SetAtomField { change: AtomFieldChange }`, etc.) rather than per-field top-level variants. The current 8c field-change enums are kept and extended as needed.
23. **Rollback mechanism** — settled as inverse replay over realized `Undo`, not molecule snapshot restore. Dense storage remains canonical; removals record `IdRemapping`, `UndoRemapping`, removed payloads, and `ConstraintUpdate` so reverse replay restores the exact previous id coordinate system.

### Mutation API revision: undo-journal transactions

Phase 8a-c have already landed the first edit vocabulary pass: `Edit`, `Action`, symbolic refs, and field-change enums. The remaining Phase 8 work revises the transaction core before later mutation phases build on it. The revision replaces snapshot rollback and `Action` with a realized undo journal that can restore dense storage without tombstones or globally stable ids.

#### Core split

- `Edit` is the caller-facing requested mutation vocabulary.
- `Undo` is the realized physical journal vocabulary.
- `Transaction` is an opaque wrapper around `Vec<Undo>`; rollback consumes it so individual undo entries are not applied out of order.
- `IdRemapping` remains the forward public remap for downstream consumers that hold ids.
- `UndoRemapping` is the inverse/restoration view used by rollback, with `atom()`, `bond()`, and relation methods as needed.

The central rollback invariant:

> When rolling back in reverse order, undoing action `i + 1` must exactly restore the id coordinate system that existed immediately after action `i`.

This is what makes dense storage compatible with inverse replay. Public `Edit` inverses are not enough: an append-only `AddAtom` cannot undo removal of a non-tail atom without changing ids. `Undo` is allowed to use private dense restore operations that are not exposed as normal edits.

#### `Edit` vocabulary

Bulk topology edits are the primitive shape because DPO naturally removes and adds sets of entities:

- `AddAtoms(Vec<AtomAst>)`
- `AddBonds(Vec<AddBond>)`, where `AddBond` carries `a: AtomRef`, `b: AtomRef`, and `ast: BondAst`
- `RemoveTopology { atoms: Vec<AtomRef>, bonds: Vec<BondRef> }`

Single-atom and single-bond add/remove calls are sugar over these bulk variants. Separate `RemoveAtom` / `RemoveBond` top-level variants are not needed for the core vocabulary.

Overlay edits stay single-item for now:

- `AddDativeBond`, `RemoveDativeBond`
- `AddAromaticSystem`, `RemoveAromaticSystem`
- `AddMulticenterBond`, `RemoveMulticenterBond`
- `AddNoncovalentBond`, `RemoveNoncovalentBond`

Bulk overlay edits can be added later if a concrete caller needs them. `RemoveTopology` still cascades overlays whose participants are removed, and that cascade is recorded in the realized `Undo`.

Field and constraint edits use the already-landed field-change vocabulary:

- `SetAtomField { idx: AtomRef, change: AtomFieldChange }`
- `SetBondField { idx: BondRef, change: BondFieldChange }`
- parallel overlay field changes
- `Set*Constraint`, `Add*Constraint`, `Remove*Constraint`
- molecule-level constraint push/pop or equivalent list edits

`Ref::New(n)` resolution is an edit-application concern, not an undo-journal indexing concern. The transaction applicator should maintain a created-entity table for symbolic refs so future internal merging of undo entries does not require `Undo[n]` to correspond to `Edit[n]`.

#### `Undo` vocabulary

`Undo` records how to return from the post-edit state to the pre-edit state. It is not constrained to look like `Edit`.

Required families:

- remove added topology: ids and payloads for atoms/bonds appended by `AddAtoms` / `AddBonds`
- restore removed topology: removed atoms, removed bonds, cascaded overlay removals, `IdRemapping`, `UndoRemapping`, and `ConstraintUpdate`
- remove added overlay relation: realized id, participants, payload
- restore removed overlay relation: original id, participants, payload, any relation-id remap, and `ConstraintUpdate`
- field rollback: inverse `*FieldChange`
- constraint rollback: `ConstraintUpdate`

`IdRemapping` is exposed per relevant `Undo` entry. Aggregate remapping and undo-entry merging are deferred optimizations; bulk edit primitives make per-entry remaps acceptable for the first version.

#### Constraint restoration

`Constraints::remap` is lossy because constraints that reference removed entities are dropped. Rollback therefore needs a patch, not a snapshot:

- `ConstraintUpdate` records constraints dropped or rewritten by a structural edit.
- It should preserve original positions when that is naturally available and cheap.
- Exact position preservation is a nice-to-have, not a semantic requirement: constraints are rendered deterministically, and non-unique entries of the same kind are conjunctive.

The consistency rule is still patch-based restoration rather than whole-constraint-list snapshotting.

#### Transaction APIs

Checked path:

```rust
pub fn transact(&mut self, edits: Vec<Edit>) -> Result<Transaction, TransactionError>
```

`transact` applies edits in order, records realized `Undo`, and on failure rolls back by reverse-replaying the collected journal before returning the apply error. If rollback itself fails, surface a combined/internal rollback error rather than silently leaving partial state.

Rollback:

```rust
impl Transaction {
    pub fn rollback(self, builder: &mut MoleculeBuilder) -> Result<(), TransactionError>;
    pub fn undos(&self) -> &[Undo];
}
```

`undos()` is read-only inspection; mutation through individual undo entries is not part of the public API.

Unchecked path:

```rust
pub fn transact_unchecked(&mut self, edits: Vec<Edit>)
```

`transact_unchecked` is for generated edit streams known to be correct, such as DPO and built-in transformations after their own prechecks. It has a separate non-journaled implementation, constructs no `Undo`, and panics on invalid refs or failed preconditions. It must not be implemented as `transact(edits).unwrap()`.

#### Private restore machinery

Rollback needs private builder operations that can rebuild dense storage into a previous coordinate system:

- remove appended atoms/bonds/relations by realized id, validating expected payloads where useful
- expand atoms/bonds through `UndoRemapping`
- reinsert removed atoms and bonds at their original dense ids
- restore cascaded overlay relations at their original relation ids
- inverse-remap surviving constraints, then reapply `ConstraintUpdate`

These operations remain private to `MoleculeBuilder`; exposing "insert at dense id" would let callers create invalid states.

## Implementation phases

Ordered phasing for the doc-86 work, with dependencies and parallelization noted.

### Phase 1 — Workspace renames (mechanical)

rust-analyzer "rename symbol" sweeps plus a handful of manual edits.

- `*Idx → *Id`: `AtomIdx`, `BondIdx`, `DativeBondIdx`, `AromaticSystemIdx`, `MulticenterBondIdx`, `NoncovalentBondIdx`, `RingIdx`
- `IdxRemapping → IdRemapping`
- `AtomRef::Idx → AtomRef::Id` and four overlay-ref enum variants
- `RingFamily::Induced → Relevant`
- View field renames: `data → ast`, `ast → molecule`
- Reader renames: `bond_order_sum → valence`, `connectivity → total_degree`, `aromatic_contribution → aromatic_valence`, `multicenter_contribution → multicenter_valence`
- Constraint variant renames: `Connectivity → TotalDegree`, `RingConnectivity → RingDegree`

**Completion**: `cargo test --workspace --tests` green. **Dependencies**: none. **Risk**: low. **Done**

### Phase 2 — Per-field accessors + per-kind constraint accessors

- `#[inline]` per-field methods on all 6 entity read-views (~25 methods); same on mut-views as read-only borrows (no `_mut` companions per the field-write decision)
- Per-kind named accessors on `AtomConstraints`, `BondConstraints`, `DativeBondConstraints`, etc.
- `is_unique()` predicate on `AtomConstraint` and `BondConstraint` (RingSize → false, all others true)
- Multi-valued `add()` path mirroring `AromaticSystemConstraints`
- Plural `ring_sizes()` accessor (multi-valued) + `get_all(kind)` / `remove_all(kind)` generics

**Completion**: new accessors callable; tests cover per-kind accessors and multi-valued add semantics. **Dependencies**: phase 1. **Risk**: low. **Done**

### Phase 3 — `ValueAst` arithmetic + Undetermined unification

- `Add`/`Sub`/`Mul`/`Div` impls on `ValueAst` (Lit/Undetermined collapse; panic on div-by-zero and overflow)
- Same on `ImplicitHydrogensAst` (`Normal → Undetermined` under arithmetic)
- Per-kind constraint accessors return `ValueAst` (Undetermined ≡ no constraint)
- Update existing call sites in `ops/*` to use the new operator surface

**Completion**: arithmetic tests cover collapse; existing valence computations migrated. **Dependencies**: phases 1, 2. **Risk**: low. **Done**

### Phase 4 — Cross-entity navigation

- `atom_view.aromatic_system() -> Option<AromaticSystemView>` (singular)
- `atom_view.dative_bonds()`, `multicenter_bonds()`, `noncovalent_bonds()` returning view iterators
- `atom_view.rings()`, `bond_view.rings()` returning `Iterator<RingView>`
- `atom_view.is_in_ring()`, `is_in_aromatic_system()`, `is_in_overlays()`
- `bond_view.is_in_aromatic_system()`, `aromatic_system()`

**Completion**: all navigation methods callable; ops migration shows reverse navigation works. **Dependencies**: phase 1. **Risk**: low–medium. **Done**

### Phase 5 — Derived readers (atom + bond + relation views)

- Atom-side: `valence`, `total_valence`, `degree`, `total_degree`, `heavy_atom_degree`, `heavy_atom_valence`, `ring_count`, `ring_size`, `ring_degree`, `ring_valence`, `multicenter_degree`, `aromatic_valence`, `multicenter_valence`, `donated_pairs`, `accepted_pairs`, `total_hydrogens`
- Bond-side: `is_in_ring`, `ring_count`, `ring_size`, `is_in_aromatic_system`, `endpoints()`, `atoms()` (view-yielding)
- Relation views: `electron_count()`, `atom_count()`, `bond_count()`, `overlapping_*`
- Add `*_ids()` accessors for constituent atoms/bonds.

**Completion**: each reader has a test; existing ops behavior unchanged. **Dependencies**: phases 1, 2, 3, 4. **Risk**: medium. **Done**

### Phase 6 — New constraint variants + DSL surface

- `TotalValence(ValueAst)` + DSL `#V`
- `RingValence(ValueAst)` + DSL `#y`
- `RingDegree(ValueAst)` (renamed from `RingConnectivity`; DSL `#x` unchanged)
- `#R!` DSL sugar for `RingCount(Lit(0))` (parser + renderer)
- Multi-valued DSL `#r5#r6` (parser exempt from duplicate-detection for non-unique kinds)
- Multi-valued EDN preserved: vector under kind already means LitSet; conjunction via multiple molecule-list entries or repeated DSL predicates

**Completion**: proptest roundtrips for each new symbol. **Dependencies**: phases 1, 2, 5. **Risk**: medium. **Done**

### Phase 7 — Ring access redesign

- `mol.rings() -> &RingSet` (canonical, cached, Vismara/22)
- `mol.rings_with(family, max_size, filter) -> RingSet` (uncached, owned; renamed from `enumerate_rings`)
- View-side sugar: canonical `is_in_ring`, `rings`, `ring_count`, `ring_size` already added in phases 4+5
- Parametric `*_from(&RingSet)` variants

**Completion**: ring access tests pass; cache behavior verified single-slot. **Dependencies**: phases 1, 4. **Risk**: low. **Done**

### Phase 8 — Edit vocabulary + transaction primitive **Done**

Largest single phase; 8a-c are implemented from the first pass, but 8d+ are revised by §"Mutation API revision: undo-journal transactions".

- **8a** `Edit` enum core variants (Add/Remove for atoms and bonds; SetAtomField for common fields) **Done**
- **8b** `Action` enum + `*Ref` enums (AtomRef, BondRef, four overlay refs) **Done**; `Action` is superseded by `Undo` in 8e.
- **8c** `*FieldChange` enums for each entity kind **Done**
- **8d** Revise `Edit` around bulk topology primitives: `AddAtoms`, `AddBonds`, `RemoveTopology`; keep single topology helpers as sugar; keep overlay edits single-item for now. **Done**
- **8e** Replace `Action` with `Undo`; add `Transaction`, `UndoRemapping`, `ConstraintUpdate`, and removed-entity payload structs. **Done**
- **8f** Add private dense restore machinery on `MoleculeBuilder` for rollback: restore-at-old-id topology, restore overlay relations, inverse-remap constraints, and apply `ConstraintUpdate`. **Done**
- **8g** Implement checked `MoleculeBuilder::transact(Vec<Edit>) -> Result<Transaction, TransactionError>` with realized undo journaling and reverse-replay rollback on failure. **Done**
- **8h** Implement non-journaled `MoleculeBuilder::transact_unchecked(Vec<Edit>) -> ()` as a separate direct-apply path; panic on invalid generated edits; do not construct `Undo`. **Done**
- **8i** Integrate topology cascade: `RemoveTopology` removes incident bonds and overlays as needed; realized `Undo` records cascades rather than reporting `Action::Cascaded`. **Done**
- **8j** Expose per-undo `IdRemapping` for downstream consumers holding ids; defer aggregate remapping and undo-entry merging. **Done**
- **8k** Existing convenience methods (`add_atom`, `add_bond`, `remove`, relation mutators) remain low-level builder primitives or are reimplemented as sugar over the appropriate checked/unchecked edit paths where atomicity is required. **Done**
- **8l** Test rollback invariant across dense non-tail removals, DPO-style bulk edits, cascaded overlay removal, constraint updates, `Ref::New` resolution, and unchecked panic behavior. **Done**; proptests behind the `proptest` feature gate cover randomized transaction cases.

**Completion**: full Edit/Undo vocab covered; checked transaction tests prove atomicity without snapshots; unchecked path avoids undo allocation; DPO-style bulk add/remove works. **Dependencies**: phases 1–7 for data structures. **Risk**: high (dense restore correctness, constraint patching, cascade behavior).

#### Phase 8 implementation plan

Implementation order is load-bearing. Land the data model and pure remapping helpers before touching transaction application; then add the checked path; then add the no-journal unchecked path.

##### Step 1 — Preserve first-pass coverage before rewrite **Done**

- Keep the existing 8a-c tests as migration scaffolding while changing types.
- Add failing tests that describe the new contract before replacing the old snapshot transaction:
  - removing a non-tail atom in a checked transaction and rolling back restores exact atom ids, bond ids, relation ids, and equality with the pre-transaction molecule
  - `RemoveTopology` of one atom cascades incident bonds and overlays into the `Undo`
  - `transact_unchecked` applies a valid generated edit stream without returning a `Transaction`
- Do not add compatibility shims for the old `Action` API; Phase 8 is not a shipped public surface.

##### Step 2 — Revise `Edit` in `umol-ast/src/ast/edit.rs` **Done**

- Replace topology single variants with bulk primitives:
  - `AddAtoms { atoms: Vec<AtomAst> }`
  - `AddBonds { bonds: Vec<AddBond> }`
  - `RemoveTopology { atoms: Vec<AtomRef>, bonds: Vec<BondRef> }`
- Add small payload structs:
  - `AddBond { a: AtomRef, b: AtomRef, ast: BondAst }`
  - later, if needed, analogous payload structs for relation additions; do not add bulk overlay variants now.
- Keep overlay relation variants single-item:
  - `AddDativeBond` / `RemoveDativeBond`
  - `AddAromaticSystem` / `RemoveAromaticSystem`
  - `AddMulticenterBond` / `RemoveMulticenterBond`
  - `AddNoncovalentBond` / `RemoveNoncovalentBond`
- Keep `Set*Field`, `Set*Constraint`, `Add*Constraint`, `Remove*Constraint`, and molecule-constraint edits.
- Remove `Edit::inverse`; rollback belongs to `Undo`.
- Add constructors for single-item ergonomics, e.g. `Edit::add_atom(ast)`, `Edit::add_bond(a, b, ast)`, `Edit::remove_atom(idx)`, `Edit::remove_bond(idx)`, all producing bulk topology variants.
- Preserve `AtomRef`, `BondRef`, and overlay refs; `New(usize)` indexes the created-entity table, not the undo vector.

##### Step 3 — Add transaction and undo data types **Done**

Add the new public transaction surface:

```rust
pub struct Transaction {
    undo: Vec<Undo>,
}

impl Transaction {
    pub fn rollback(self, builder: &mut MoleculeBuilder) -> Result<(), TransactionError>;
    pub fn undos(&self) -> &[Undo];
}
```

Replace `Action` with `Undo`:

```rust
pub enum Undo {
    RemoveAddedTopology {
        atoms: Vec<AddedAtom>,
        bonds: Vec<AddedBond>,
    },
    RestoreTopology {
        atoms: Vec<RemovedAtom>,
        bonds: Vec<RemovedBond>,
        overlays: RemovedOverlays,
        remapping: IdRemapping,
        undo_remapping: UndoRemapping,
        constraint_update: ConstraintUpdate,
    },
    RemoveAddedDativeBond(AddedDativeBond),
    RestoreRemovedDativeBond(RemovedDativeBond),
    // same single-item pattern for aromatic, multicenter, noncovalent
    SetAtomField { id: AtomId, change: AtomFieldChange },
    SetBondField { id: BondId, change: BondFieldChange },
    // same field pattern for overlays
    ApplyConstraintUpdate(ConstraintUpdate),
}
```

Payload structs carry realized ids plus enough data to validate rollback:

- `AddedAtom { id, ast }`
- `AddedBond { id, endpoints, ast }`
- `RemovedAtom { id, ast }`
- `RemovedBond { id, endpoints, ast }`
- `Removed*Relation { id, participants, ast }`
- `RemovedOverlays` groups cascaded dative, aromatic, multicenter, and noncovalent removals.

Add `Undo::id_remapping() -> Option<&IdRemapping>` so downstream consumers can inspect per-undo forward remaps without needing to match every variant manually.

##### Step 4 — Extend remapping support **Done**

Keep `IdRemapping` as forward old-to-new mapping for downstream consumers.

Add `UndoRemapping` as the rollback-side inverse view:

- `UndoRemapping::from(&IdRemapping)` or `IdRemapping::undo_remapping()`
- `atom(post: AtomId) -> AtomId`
- `bond(post: BondId) -> BondId`
- relation methods as needed

`UndoRemapping` maps surviving current ids back into the pre-removal coordinate system. Removed ids are restored from the `Removed*` payloads, not looked up through `UndoRemapping`.

Implementation detail:

- For atom/bond ids, invert the dense "remove sorted ids and shift left" operation.
- For relation ids, use the same dense inverse over the relation-specific removed-id lists already carried by `IdRemapping`.
- Keep these helpers pure and unit-tested before using them inside rollback.

##### Step 5 — Add `ConstraintUpdate` **Done**

Add patch-based restoration for molecule-level constraints:

```rust
pub struct ConstraintUpdate {
    dropped: Vec<DroppedConstraint>,
    rewritten: Vec<RewrittenConstraint>,
}
```

The exact internal shape can be tuned during implementation, but it must support:

- constraints dropped because they referenced removed atoms, bonds, or relations
- constraints rewritten by `IdRemapping`
- inverse rewrite during rollback
- original positions where cheaply available

Preferred implementation path:

- add a `Constraints::remap_with_update(&mut self, remap: &IdRemapping) -> ConstraintUpdate` helper
- internally walk the original list with positions
- for each constraint:
  - if remap succeeds and the constraint changes, store a rewrite record
  - if remap fails, store a dropped record
- rollback inverse-remaps surviving constraints, then reinserts dropped constraints by recorded position when possible

Do not snapshot the whole `Constraints` list as the normal path.

##### Step 5a — Directional update/remap APIs **Done**

Step 6 exposed two small refinements before the restore helpers become load-bearing transaction code:

- Add rollback-facing methods on `ConstraintUpdate`, e.g. `rollback_into(&mut Constraints)`, so the direction is explicit. `ConstraintUpdate` records the forward effects of a remap; rollback applies the inverse.
- Add explicit relation-remap context for single relation restore undos. Either:
  - make `RestoreRemovedDativeBond` / `RestoreRemovedAromaticSystem` / etc. carry the relation-only `UndoRemapping`, or
  - add a small relation-only remap payload alongside the removed relation payload.
- Update the single-relation `Undo` variants to use that explicit context instead of letting builder restore helpers synthesize ad hoc remappings from the removed id.
- Keep public `Edit` unchanged; this is a realized-undo/private-restore refinement.

This step should happen before Step 6 is wired into `transact`, but it does not require expanding the caller-facing mutation vocabulary.

##### Step 6 — Add private dense restore operations **Done**

Add private builder methods in `umol-ast/src/ast/molecule/builder.rs` or a private helper module:

- `restore_topology(removed_atoms, removed_bonds, removed_overlays, undo_remapping, constraint_update)`
- `remove_added_topology(added_atoms, added_bonds)`
- `restore_dative_bond(removed)`
- `remove_added_dative_bond(added)`
- same single-relation helpers for aromatic, multicenter, noncovalent

Restoring topology rebuilds dense storage rather than inserting into vectors in place:

- create a new atom vector of pre-removal length
- place removed atoms at their original ids
- map surviving current atoms back through `UndoRemapping`
- rebuild graph endpoints in old bond-id order from surviving bonds plus removed bonds
- rebuild bond data in old bond-id order
- rebuild overlay relation sets in old relation-id order
- restore molecule constraints through `ConstraintUpdate`

These helpers must be private. They are not mutation API.

##### Step 7 — Implement checked apply **Done**

In `umol-ast/src/ast/molecule/transact.rs`, replace snapshot rollback with journaled apply:

```rust
pub fn transact(&mut self, edits: Vec<Edit>) -> Result<Transaction, TransactionError> {
    let mut journal = Vec::new();
    let mut created = CreatedEntities::new();
    for edit in edits {
        match self.apply_edit_journaled(edit, &mut created) {
            Ok(undo) => {
                created.record(&undo);
                journal.push(undo);
            }
            Err(apply_error) => {
                if let Err(rollback_error) = rollback_journal(self, journal) {
                    return Err(TransactionError::RollbackFailed {
                        apply: Box::new(apply_error),
                        rollback: Box::new(rollback_error),
                    });
                }
                return Err(apply_error);
            }
        }
    }
    Ok(Transaction { undo: journal })
}
```

Actual error handling should preserve the original apply error unless rollback fails; if rollback fails, return a transaction error variant that carries both contexts or at least makes the rollback failure explicit.

Journaled apply rules:

- `AddAtoms`: apply as appends, record `Undo::RemoveAddedTopology`.
- `AddBonds`: resolve refs through `CreatedEntities`, apply as appends, extend `RemoveAddedTopology` for the edit.
- `RemoveTopology`: resolve refs, capture all removed payloads and constraints, call low-level dense `remove`, build `IdRemapping` and `UndoRemapping`, record `RestoreTopology`.
- field edits: validate old value, set new value, record inverse field change.
- constraint edits: validate precondition, mutate, record `ConstraintUpdate`.
- overlay add/remove: record single-item added/removed relation undo.

`CreatedEntities` records realized ids for every created entity in edit order. It is separate from `Vec<Undo>` so future undo merging does not break `Ref::New`.

##### Step 8 — Implement rollback **Done**

`Transaction::rollback(self, builder)` reverse-iterates `undo`.

Rules:

- rollback is allowed only through `Transaction`, not through public `Undo::apply`
- each undo validates obvious preconditions before mutating when validation is cheap
- rollback of a failed transaction and user-requested rollback use the same private function
- after each undo, the id coordinate system must equal the system before that forward edit

Core tests should assert equality after every rollback, not only final counts.

##### Step 9 — Implement unchecked apply as a separate path **Done**

Add `transact_unchecked(Vec<Edit>) -> ()` with direct, non-journaled application.

Rules:

- no `Undo` allocation
- no removed-payload capture unless needed to perform the mutation itself
- panic on invalid refs, out-of-range ids, or old-state mismatches
- do not call `transact(...).unwrap()`
- reuse small pure helpers for ref resolution if they do not force undo capture

Unchecked `RemoveTopology` can call the low-level dense removal directly after resolving ids; it does not compute `ConstraintUpdate` except for whatever constraint remapping the builder already performs.

##### Step 10 — Reconcile low-level builder APIs **Done**

Keep low-level builder primitives where they are useful internally:

- `add_atom`, `add_bond`, overlay `add_*`
- `remove(&[AtomId], &[BondId]) -> IdRemapping`
- overlay `remove_*`
- `atom_mut`, `bond_mut`, relation mut views

Direct builder methods are preferred for structural removal. Do not add `entity_mut().remove()`; mutable views edit entity fields/constraints, while removal belongs to the owning builder because dense deletion can remap many unrelated ids.

Document that these are non-transactional builder operations. Higher-level atomic mutation goes through `transact`; generated DPO paths may use `transact_unchecked`.

Do not force all convenience methods through checked `transact` if doing so would allocate `Undo` for simple builder construction.

##### Step 11 — Test gates **Done**

Minimum focused tests for Phase 8 completion:

- `test_molecule_builder_transact_add_atoms`
- `test_molecule_builder_transact_add_bonds_new_refs`
- `test_molecule_builder_transact_remove_topology`
- `test_transaction_rollback`
- `test_transaction_rollback_non_tail_atom`
- `test_transaction_rollback_non_tail_bond`
- `test_transaction_rollback_cascaded_overlays`
- `test_transaction_rollback_constraint_update`
- `test_molecule_builder_transact_error_rolls_back`
- `test_molecule_builder_transact_unchecked`
- `test_molecule_builder_transact_unchecked_error`
- `test_undo_remapping_atom`
- `test_undo_remapping_bond`
- relation-id equivalents for overlay removals

Property-style tests can come after the focused cases:

- generate small molecules, remove arbitrary atom/bond subsets through checked `transact`, then rollback and assert exact equality
- compare checked `transact` and `transact_unchecked` final states for known-valid edit batches

##### Step 12 — Migration cleanup **Done**

- Remove `Action` once all tests and call sites are migrated to `Undo`.
- Remove old `Edit::inverse` tests.
- Update doc examples and discussion references from `Action::Cascaded` to `Undo::RestoreTopology`.
- Run `cargo test -p umol-ast --tests`; then `cargo test --workspace --tests` after downstream call sites are migrated.

### Phase 9 — Mutation API completion **Done**

- `MoleculeBuilder::atom_mut(id) -> AtomBuilderViewMut` etc. for all entity kinds (DPO Phase-2) **Done**
- `MoleculeBuilder::remove_*` overlay-mutators wired in (DPO Phase-3) **Done**
- Internal relation mutation via Edit batch (no new API; documented usage pattern) **Done**
- `atom.is_in_overlays()` umbrella predicate **Done**

**Completion**: reaction-rule application end-to-end on a simple example — gated by DPO transformer in `umol-graph/ops/transformer/`, which is not yet present (only `aromatizer.rs` and `kekulizer.rs` exist there). The `umol-ast` surface required for it is in place. **Dependencies**: phase 8. **Risk**: medium.

### Phase 10 — Subgraph + projection **Done**

- `induced_subgraph(atoms) -> MoleculeEmbedding<'_>` **Done**; borrows `&MoleculeAst`, carries six local→parent maps + parent→local inverse for atoms, `extract()` materializes the sub-AST on demand, `edits()` derives the `RemoveTopology`. Same type intended for future subgraph isomorphism results.
- `bonds().induced(atoms) -> Vec<BondView>` / `bonds().induced_ids(atoms) -> Vec<BondId>` **Done** (on `BondViews`, not `MoleculeAst`)
- `GraphView::connected_components(alg)` **Done**; `GraphView::connected_components_in(atoms, alg)` **Done**
- `bonds().connecting(a, b) -> Option<BondView>` / `bonds().connecting_id(a, b) -> Option<BondId>` **Done** (on `BondViews`, not `MoleculeAst`)

**Completion**: subgraph extraction works for kekulizer / HMO callers. **Dependencies**: phase 8 (uses Edit). **Risk**: medium.

### Phase 11 — State predicates + small tail items **Done**

- `mol.is_empty()`, `has_constraints()` (molecule-scope), `has_overlays()`, `has_dative_bonds()`, `has_aromatic_systems()`, `has_multicenter_bonds()`, `has_noncovalent_bonds()` **Done**
- `NeighborView` reshape: private `{ atom_id, bond_id, molecule }` + `.atom() -> AtomView` / `.bond() -> BondView` methods; `atom_id()` / `bond_id()` accessors **Done**
- Drop `ast: &BondAst` field from NeighborView; callers route through `.bond().ast` **Done** (aromatizer, rewrite, atom-view methods, tests migrated)

**Completion**: predicates callable; NeighborView calls in ops/ migrated. **Dependencies**: phases 1, 4. **Risk**: low.

### Phase 12 — Validation tier-2 stubs filled in

- `TotalCharge`, `TotalSpin`, `AromaticElectronCount` evaluators in `ConstraintValidator`
- `SpinCouplingValidator` real implementation (currently stub)
- `ConstraintValidator` cross-checks

**Completion**: per-test molecules pass tier-2 validation. **Dependencies**: phases 1–6. **Risk**: low–medium (chemistry-dependent).

### Phase 13 — Ops migration + proptest update

- Migrate `ops/aromaticity/*`, `ops/valence/*`, `ops/transformer/*`, `ops/validator/*` to the new API surface
- Update proptest generators and assertions for renamed types and new constraint variants
- Verify all existing benchmarks still pass

**Completion**: `cargo test --workspace --tests` green; benchmarks unchanged or improved. **Dependencies**: all earlier phases. **Risk**: medium.

#### Resolution vs. validation

Every op that handles constraint-vs-state interaction splits into two distinct passes. Conflating them was the original sin behind the `satisfies` indirection (Phase 14, Step 5) and the legacy permissive `atom_is_aromatic`.

- **Resolution** translates *declared constraints* into *derived state* using a model. It runs once when lowering from the declarative AST into the resolved hierarchy. Inputs: the `*Constraints` containers carrying user-supplied assertions. Outputs: concrete state on the entity (e.g. aromatic-system membership, atom-type assignment, ring set). The model is the input table or algorithm (e.g. aromaticity model, valence registry, ring-perception strategy); resolution may consult constraints as hints (pinned π count, declared aromaticity) but is not driven by them — declarations are *constraints to satisfy*, not *state to propagate*. Once resolved, the derived state is the source of truth; constraints remain on the entity but are no longer consulted to answer "what is the atom's aromaticity / valence / ring membership."
- **Validation** compares *declared constraints* against *resolved state* for consistency. It runs after resolution. For each constraint of kind K, derive the effective state value at K from the entity and call `constraint.matches(&effective)` — pure lattice algebra (Phase 14, Step 5). Pattern `Undetermined` (vacuous) trivially matches. Mismatches surface as op-specific errors (`AromaticityValidationError`, `ValenceValidationError`, ...). Validation never mutates state.

Concretely per op:

- `ops/aromaticity`: `AromaticityPerception::resolve` (constraints → memberships); `AromaticityPerception::validate` (memberships ↔ constraints check). `view.is_in_aromatic_system()` is the post-resolution accessor; `view.constraints().aromatic_valence()` returns the declared constraint untouched. The two never combine via disjunction.
- `ops/valence`: `counts.rs` and `atom_typing.rs` are resolution. `ops/validator/valence.rs` (when it lands) is validation; it checks the resolved candidate against the declared `Valence` / `TotalValence` / etc. constraints.
- `ops/transformer` (kekulize, aromatize): mutates resolved state. After mutation, run validation to confirm declared constraints still hold.
- `ops/validator`: pure validation; never resolves.

Implementation rules:

- Resolution APIs return resolved state; they do not return validation verdicts. If resolution detects an unsatisfiable declared constraint (no resolution exists), that's a resolution error (e.g., `CountsError::NoValidValenceState`), not a validation error.
- Validation APIs return `Result<(), ValidationError>`; they never mutate. The pattern is always `for each constraint c in entity: c.matches(&derived_state_at(c.kind()))`.
- Neither pass calls `satisfies` (it does not exist). Neither pass mixes "is the constraint asserted?" with "is the state consistent with it?" — those are answered by `entity.constraints().get(...)` and `constraint.matches(&derived)` respectively.

### Phase 14 — Lattice ops + valence resolver cleanup

Adds a `Lattice` trait covering the standard partial-order/lattice operations on AST refinement types, replaces the misnamed `narrow_atom`/`lift_constraints` flow in the valence resolvers with derived `Lattice` operations, and deletes `umol-graph/src/ops/valence/shared.rs` entirely.

**Lattice trait** (lives in `umol-ast/src/ast/traits.rs` alongside `AsLit`):

- `Lattice` trait with `is_top`, `is_bottom`, `matches(&target)`, `meet(&other) -> Option<Self>`, `join(&other) -> Self`, `narrow_from(&mut self, &other) -> bool`, `widen_with(&mut self, &other) -> bool`
- Implement on `ValueAst`, `IsotopeAst`, `ImplicitHydrogensAst`, `SpinStateAst`, `ElementAst`, `AromaticValenceAst`, `MulticenterValenceAst`
- `Undetermined` is top; ground variants are bottom; `Normal`/`Natural` sentinels sit between `Undetermined` and `Lit`
- `meet(Lit, Lit) = Some(Lit)` iff equal; `meet(LitSet, Lit)` iff Lit ∈ Set; `meet(LitSet, LitSet)` returns normalized intersection
- `LitSet` canonical form: first-occurrence preserving, in-place dedup on the `Vec<i64>`, no auxiliary structure
- `Expr` only meets/joins with itself (syntactic equality) or `Undetermined`; no symbolic evaluation
- `SpinStateAst::meet` is pure field-wise (no physics validation; matches the existing `is_ground` policy)
- `#[derive(Lattice)]` proc-macro for struct types (`AtomAst`, `BondAst`, ...) generates field-by-field composition
- `le()` deliberately omitted to avoid conflict with derived `PartialOrd` (which is structural); the underlying tension with derived `PartialOrd`/`Ord` on refinement types is parked for a follow-on cleanup

**`AsLit` trait** (`umol-ast/src/ast/traits.rs`):

- `pub trait AsLit { type Lit; fn as_lit(&self) -> Option<Self::Lit>; /* derived methods with defaults */ }`
- Derived defaults: `as_lit_ok_or`, `as_lit_ok_or_else`, `as_lit_or`, `as_lit_or_else`, `as_lit_expect`
- Implementations: `ValueAst → i64`, `ElementAst → Element`, `IsotopeAst → u32`, `ImplicitHydrogensAst → i64`, `SpinStateAst → SpinState` (physics-validated; strictly narrower than `is_ground`), `AromaticValenceAst → i64` (extracts inner of `Aromatic(_)`), `MulticenterValenceAst → i64` (same for `Multicenter(_)`)
- Replaces the per-type inherent `as_lit_or*` family on the four existing types; eliminates ~50 lines of boilerplate

**View additions** (`umol-ast/src/ast/views/atom.rs`):

- ~~`AtomView::satisfies(&AtomConstraint)`~~ — **abandoned** (see Step 5). The "constraint holds against atom state" check is purely lattice algebra: `q.matches(&effective)` where `effective` is derived per-kind from the atom view. Each validator inlines its own derivation + match.

- `matches` lives on the `Lattice` trait (Step 5). Required method, no default impl; each AST type provides its own direct impl.

**`NormalImplicitHydrogensTable`** (new config type):

- Split `normal_valence` and `aromatic_normal_valence` out of `ValenceTable` into a new struct
- Methods: `normal_valence_for(element, charge) -> Option<u8>`, `aromatic_normal_valence_for(element, charge) -> Option<u8>`, `infer(element, charge, valence, is_aromatic) -> Option<u8>`
- `ValenceModel::AtomTyping { registry, normal_implicit_hydrogens }`
- `ValenceModel::Counts { table, normal_implicit_hydrogens, allow_implicit_hydrogens }`
- New TOML: `default-normal-implicit-hydrogens.toml`; `default-valence-table.toml` loses the two fields
- `ValenceTable::compute_implicit_hydrogens` stays where it is (uses `allowed_valences`)

**Delete `shared.rs`** — per-item disposition:

- `charge_or_zero` → inline `atom.charge.as_lit_or(0) as i8`
- `aromatic_pi_pinned` → `AromaticValenceAst::as_lit()` or `ValueAst::matches`
- `atom_is_aromatic` → resolved by perception. Atoms in aromatic systems are members via the `aromatic_systems` relation; `view.is_in_aromatic_system()` is the ground-truth predicate post-perception. The legacy "system OR declared" disjunction is replaced by perception's two-pass structure (resolve reads declarations as inputs; validate cross-checks against memberships).
- `ground_spin_state` → `SpinStateAst::as_lit()`
- `spin_state_undetermined` → existing `SpinStateAst::is_undetermined()`
- `value_matches`, `spin_matches` → existing `matches` on each type
- `base_atom_compatible` → three chained `matches` inline
- `pattern_constraints_compatible` → `.iter().all(...)` inline
- `atom_constraint_holds` → deleted entirely. Constraint-vs-state checking is owned by validation passes (e.g., `AromaticityPerception::validate`), which inline the relevant `q.matches(&effective)` calls directly.
- `narrow_value`, `narrow_atom` → `Lattice::narrow_from` (via derive on struct types)
- `lift_constraints`, `narrow_atom_constraint`, `narrowable` → `AtomConstraints::narrow_with(&other: &AtomConstraints) -> bool` (bulk form)
- `try_build_candidate`, `resolve_unpaired_lone_pairs` → private fns in `counts.rs` (counts-only consumers)
- `atom_dative_counts` → inline at single callsite
- `infer_normal_*_implicit_hydrogens` → `NormalImplicitHydrogensTable::infer`
- `bond_order_lit`, `atom_view` → delete (dead code)
- `AtomCandidate` → dissolved; `candidates_for` returns `Vec<AtomAst>` with constraints embedded in `.constraints`; callsite reduces to `atom_mut.narrow_from(&candidate)`
- `atom_typing.rs` locals (`pattern_implicit_h_compatible`, `collect_pattern_constraints`) — kept in place for now; review in a follow-on pass

**Error payloads**: `CountsError`/`AtomTypingError` carry `ValueAst` (charge/valence) directly; all ASTs derive `Debug`.

**Tests**: existing `atom_constraint_holds` tests (shared.rs:440-538) deleted with the function (no `satisfies` to port to). Aromaticity perception's `validate` will grow its own tests when implemented (separate phase, not part of Phase 14). New tests for `Lattice` impls — exhaustive variant cross-products for `meet`/`join` on `ValueAst`, sentinel-bearing types, and tagged sum types.

**Completion**: `umol-graph/src/ops/valence/shared.rs` deleted; `Lattice` trait covers all value-type and entity ASTs; `cargo test --workspace --tests` green. **Dependencies**: phase 3 (`ValueAst` arithmetic), phase 13 (ops migration of `shared.rs` callers to the new API surface). **Risk**: medium-high — touches many AST types and both valence resolvers; lattice cross-product tables for `ValueAst` are the load-bearing correctness check.

#### Phase 14 implementation plan

Sequencing chosen so each step lands independently with `cargo test --workspace --tests` green and no temporary shims.

##### Step 1 — `AsLit` trait + extensions **Done**

- Add `AsLit` trait to `umol-ast/src/ast/traits.rs` with associated `type Lit`, required `fn as_lit(&self) -> Option<Self::Lit>`, and default impls for `as_lit_or` / `as_lit_or_else` / `as_lit_ok_or` / `as_lit_ok_or_else` / `as_lit_expect`
- Convert the four existing inherent implementations (`ValueAst`, `ElementAst`, `IsotopeAst`, `ImplicitHydrogensAst`) to `AsLit` impls; remove the boilerplate inherent methods
- Add `AsLit` impls for `SpinStateAst` (`type Lit = SpinState`, validates physics parity via `SpinState::try_new`), `AromaticValenceAst` (`type Lit = i64`, extracts inner of `Aromatic(_)`), `MulticenterValenceAst` (`type Lit = i64`)
- Update consumer modules with `use umol_ast::ast::AsLit;` so the trait methods are in scope at call sites
- Per-type unit tests covering the variant cross-product (including `SpinStateAst` parity-invalid case)

##### Step 2 — `Lattice` trait + manual impls (value-type ASTs) **Done**

- Add the `Lattice` trait to `umol-ast/src/ast/traits.rs` (alongside `AsLit`)
- Trait surface: `is_undetermined`, `is_ground`, `meet`, `join` required; `narrow_from`, `widen_with` default impls derived from `meet`/`join` + `PartialEq`
- Remove inherent `is_undetermined` / `is_ground` on each AST type — they become trait methods. Callers add `use Lattice;` to keep the call sites working
- `matches` stays inherent-only — not in the trait
- Manual impls for `ValueAst`, `IsotopeAst`, `ImplicitHydrogensAst`, `SpinStateAst`, `ElementAst`, `AromaticValenceAst`, `MulticenterValenceAst`
- Exhaustive `meet`/`join` tests over variant cross-products — load-bearing correctness check
- Inline `LitSet` first-occurrence dedup helper (private to `traits.rs`)

##### Step 3 — `Lattice` derive macro for struct types **Done**

- Hand-rolled `Lattice` impls on the six entity types (`AtomAst`, `BondAst`, `DativeBondAst`, `MulticenterBondAst`, `NoncovalentBondAst`, `AromaticSystemAst`) — proc-macro path deferred; no `umol-ast-macros` crate exists yet
- Added `Lattice` impls on the six constraint containers (`AtomConstraints`, `BondConstraints`, `DativeBondConstraints`, `MulticenterBondConstraints`, `NoncovalentBondConstraints`, `AromaticSystemConstraints`) — container `meet` is per-kind merge with inner value `meet` and vacuous-entry pruning; `join` is per-kind intersection
- Removed inherent `is_ground` / `is_undetermined` on the six entity types; trait now owns them. Added `use Lattice;` to consumer modules (resolver, atom_typing, counts, molecule.rs internal)
- Non-Lattice field policies confirmed: `DativeBondAst::acceptor_slot` (u8) and `NoncovalentBondAst::kind` (NoncovalentBondKindAst, equality-only) → `meet` requires equality (else `None`); `join` widens to `Self::default()` on mismatch
- `Vec<ValueAst>` length-mismatch (`AromaticSystemAst::electrons`, `MulticenterBondAst::electrons`) → `meet` `None`; `join` `Self::default()`
- Container `is_undetermined`: every entry's value is `is_undetermined` (vacuous on empty); `meet`/`join` outputs prune undetermined-valued entries (canonical: no vacuous entries)
- Tests: representative cross-product tests on `BondAst` (meet, join, narrow_from positive/negative cases); 8063 total tests pass

##### Step 4 — `AtomConstraints::narrow_with(&other) -> bool` **Done**

- Subsumed by Step 3's `Lattice` impl on `AtomConstraints`: `narrow_with` is the in-place version of `meet`, which is `Lattice::narrow_from` (default impl)
- Added 14 tests per constraint kind (empty/empty, add new kind, narrows undetermined to lit, lit/lit match preserved, lit/lit mismatch → `None`, multi-kind combines, aromatic-valence narrows, aromatic-valence not vs aromatic → `None`, RingSize union, RingSize dedup, vacuous-entry pruning, `narrow_from` extends/no-change/contradiction-leaves-unchanged, `join` keeps shared kinds, `join` widens value)
- Found and fixed a pre-existing bug in `AtomConstraints::get_all` and `remove_all`: they used `binary_search` (`find`) which can return any matching index in a multi-entry cluster (e.g., multiple `RingSize` entries). Switched to `partition_point` to find the leftmost cluster start.

##### Step 5 — Move `matches` into `Lattice` trait **Done**

(Supersedes the earlier `AtomView::satisfies` work, which was walked back: `satisfies` was added, recognized as conflating multiple semantics — legacy permissive fallback, system-OR-declared disjunction, vacuous-`Undetermined` handling — and deleted. The membership-vs-constraint interaction lives only inside validation passes; validators inline `q.matches(&effective)` per-kind directly. Atom typing (Step 8) uses lattice `meet`/`narrow_from` directly. No `satisfies` to call anywhere.)

What was actually done in this step:

- Added `fn matches(&self, target: &Self) -> bool` as a required method on the `Lattice` trait, with no default impl. Semantically the partial-order check: `pattern.matches(target) ⇔ pattern.meet(target) == Some(target)` (modulo canonicalization of set-valued representations — canonicalization will be added later if a caller needs it; existing impls keep their current loose set-matching semantics for now).
- Moved every inherent `pub fn matches` into its `impl Lattice for ...` block (no behavioral change for the leaf/value types). Impls touched: `ValueAst`, `ElementAst`, `IsotopeAst`, `ImplicitHydrogensAst`, `SpinStateAst`, `AromaticValenceAst`, `MulticenterValenceAst`, `AtomAst`, `BondAst`, `DativeBondAst`, `MulticenterBondAst`, `NoncovalentBondAst`, `AromaticSystemAst`.
- `NoncovalentBondKindAst::matches` stays inherent — that type doesn't implement `Lattice` (its parent treats `kind` as a structural anchor inline).
- Entity-level `matches` now field-wise-checks `constraints` too (previously omitted because there was no `Constraints::matches`). `DativeBondAst::matches` also now enforces `acceptor_slot` equality — consistent with `meet`, fixes an asymmetry where the prior inherent `matches` ignored the anchor.
- Added `Lattice::matches` impls on the constraint containers (which already had `Lattice`): `AtomConstraints`, `BondConstraints`, `AromaticSystemConstraints`, `MulticenterBondConstraints`, `DativeBondConstraints`, `NoncovalentBondConstraints` (trivially `true`; uninhabited inner type). Semantics: field-wise per-kind via existing accessors (absent constraint = vacuous `Undetermined`, matches anything); multi-valued kinds (`RingSize` on atom and bond) check that every pattern entry is matchable by some target entry.
- New tests: `test_aromatic_valence_ast_matches`, `test_multicenter_valence_ast_matches`, `test_atom_constraints_matches`, `test_bond_constraints_matches`. Extended `test_atom_ast_matches`, `test_bond_ast_matches`, `test_dative_bond_ast_matches`, `test_multicenter_bond_ast_matches`, `test_aromatic_system_ast_matches`, `test_noncovalent_bond_ast_matches` with cases covering previously-untested dimensions: spin, `acceptor_slot`, constraints presence/absence/value-mismatch, set-pattern kind matching.
- Validation: 2604 umol-ast lib tests pass (up from 2550 — 54 new); 8126 workspace tests pass.

##### Step 6 — Config split: `NormalImplicitHydrogensTable`

- Create `umol-graph/src/ops/normal_implicit_hydrogens.rs`
- Strip `normal_valence`, `aromatic_normal_valence` from `ValenceEntry`/`ValenceTable`
- Update `ValenceModel` variants: `AtomTyping { registry, normal_implicit_hydrogens }`, `Counts { table, normal_implicit_hydrogens, allow_implicit_hydrogens }`
- Split TOML: write `default-normal-implicit-hydrogens.toml`; strip the two fields from `default-valence-table.toml`
- Re-resolve the static `LazyLock` defaults
- Touch resolver constructors to accept the new table

This temporarily breaks `shared.rs`'s `infer_normal_*` calls — those go away in Steps 7–8.

##### Step 7 — Migrate `counts.rs`

- Move `try_build_candidate`, `resolve_unpaired_lone_pairs` from `shared.rs` into `counts.rs` as private fns; rewrite their internal uses of `infer_normal_aromatic_implicit_hydrogens` to call `NormalImplicitHydrogensTable::infer`
- The aromatic-branch decision (legacy `atom_is_aromatic`) becomes `view.is_in_aromatic_system()` alone. Inline `AromaticValence(Aromatic(_))` declarations are *constraints*, not state — they don't drive counts' branch decision. Perception runs first and resolves declarations to memberships; counts reads memberships.
- Inline `aromatic_pi_pinned` via `AromaticValenceAst::as_lit()` (still needed for the pinned-π case where counts walks constraints to find user-supplied pin values for candidate filtering)
- Dissolve `AtomCandidate`: `candidates_for` returns `Vec<AtomAst>` with synthetic `Valence(Lit(v))`/`AromaticValence(Aromatic(Lit(a)))` constraints pushed into each candidate's `.constraints`
- Callsite collapses to `atom_mut.narrow_from(&candidate)`
- `CountsError::NoValidValenceState` carries `ValueAst` for charge/valence

##### Step 8 — Migrate `atom_typing.rs`

Atom typing's job is classification + merge. The clean flow:

- For each candidate pattern, **compatibility check** via the lattice: `pattern.meet(&atom_view.ast).is_some()` (or equivalent `pattern.matches(...)` form once that lifts into the trait). Constraints are part of the entity Lattice; absent constraint on atom = vacuous Undetermined, so meet succeeds when the pattern's constraints are compatible (atom doesn't have to pre-declare them).
- **Merge**: `atom_mut.narrow_from(&pattern)` — narrows base fields and merges constraints in one lattice operation (per Q1).
- No `satisfies` call. No `pattern_constraints_compatible`. No `base_atom_compatible`.
- Replace `atom_dative_counts(view)` with two `as_lit_*` calls inline at the single use site
- Replace `infer_normal_implicit_hydrogens(...)` with `normal_implicit_hydrogens.infer(...)`
- Dissolve `AtomCandidate`
- `AtomTypingError::NoMatchingPattern` carries `ValueAst` for charge
- `pattern_implicit_h_compatible`, `collect_pattern_constraints` evaporate — pattern constraints flow through `narrow_from`; implicit-H pattern check folds into the lattice meet on `ImplicitHydrogensAst`

##### Step 9 — Delete `shared.rs`

- Remove `mod shared;` from `umol-graph/src/ops/valence.rs`
- Delete `umol-graph/src/ops/valence/shared.rs`
- `cargo test --workspace --tests`; `cargo build --workspace --all-targets`
- Verify: `grep -rn "shared::" umol-graph/src/ops/valence/` returns empty

##### Validation gates per step

- Each step lands green: `cargo test --workspace --tests`, no new clippy warnings
- Step 2 also: exhaustive variant cross-product tests for `meet`/`join` on `ValueAst` (foundation for everything else)
- Step 5 supersedes the earlier `satisfies` work (which was deleted). Regression coverage for `atom_constraint_holds` evaporates with `shared.rs`; equivalent checks are reimplemented per validator (e.g., `AromaticityPerception::validate` when that lands). New tests landed for `matches` on every Lattice impl (leaf types already had them; constraint containers and the constraint-dimension cases on entity types are new).
- Step 9 also: `grep` confirms zero `shared::` references remain

##### Risks

- **Lattice derive macro** (Step 3): if `umol-ast-macros` doesn't exist, this step doubles in scope (new proc-macro crate setup). Hand-rolled impls on the ~6 entity types is the fallback (~150 lines).
- **TOML split** (Step 6): config file format change is breaking. Contained — only the two default TOML files are involved; no external consumers yet.
- **`ImplicitHydrogensAst::Normal` sentinel in `meet`** (Step 2): semantic is `Normal ⊐ Lit(n)`, so `meet(Normal, Lit(n)) = Lit(n)` — assumes the H-demand table is consulted *elsewhere*, not at meet time. Worth a doc comment on the impl.

### Parallelization

After phase 1, two tracks:

- **Read-side**: phases 2 → 3 → 4 → 5 → 6 → 7 → 11
- **Write-side**: phase 8 (once read-side has stabilized enough) → 9 → 10

Phase 13 ops migration can begin incrementally as phase 2 lands, continuing through phases 3–7.

### Rough effort estimate

| Phase | Estimate |
|---|---|
| 1 | 1–2 days |
| 2 | 2–3 days |
| 3 | 1–2 days |
| 4 | 2–3 days |
| 5 | 3–5 days |
| 6 | 2–3 days |
| 7 | 1 day |
| 8 | 2–3 weeks (revised undo-journal design; subdivided) |
| 9 | 3–5 days |
| 10 | 2–4 days |
| 11 | 1–2 days |
| 12 | varies; chemistry-dependent |
| 13 | 1–2 weeks |
| 14 | 1–2 weeks |

Net: roughly 7–12 weeks of focused work, parallelizable to ~5–7 weeks with two tracks.

### Per-phase validation gate

Each phase lands with:

- New tests covering its additions
- `cargo test --workspace --tests` green
- No new clippy warnings
- Ops modules either migrated (in their respective phases) or passing via temporary shims

### `Undo` replaces `Action` for realized rollback

The first Phase-8 pass used `Action` as a lightweight record of what happened (`*Added`, `Done`, `Cascaded`). That is enough to report realized ids, but not enough to roll back dense storage without a snapshot. Removing a non-tail atom or bond compacts ids; undoing that operation cannot be expressed as a normal append-only `Edit::AddX`.

**Resolution:** `Action` is replaced by `Undo`, a physical rollback journal entry. `Undo` records removed payloads, `IdRemapping`, `UndoRemapping`, cascaded overlay removals, and `ConstraintUpdate` as needed. `Transaction` wraps `Vec<Undo>` and rollback consumes the wrapper, reverse-replaying entries in order.

`Edit` remains the caller-facing mutation vocabulary. It does not need to be self-inverting, and `Edit::inverse` is no longer part of the transaction design. The checked path builds `Undo`; the unchecked path applies trusted edits directly and builds no journal.

## AST-vs-API layering: parked considerations

The ring-view discussion in §"Ring access" surfaced a broader question: what would a chemist-facing API tier look like, given that the previous attempt was unwound? This section preserves the design considerations for when the question becomes actionable. Nothing here is a current plan — the actionable space right now is the AST API itself (the operation taxonomy above).

### State after the doc 92 restructure

The chemist-facing wrappers (`Molecule`, `Pattern`) were deleted because the previous layering was unsettled — operations crossed back and forth between AST and API tiers without a principled separating line. The crates were split (`umol-ast` for the algebraic tier, `umol-graph` for everything else) to enforce a hard physical boundary while the design re-stabilizes. Current boundary:

- **`umol-ast`** holds `MoleculeAst`, `AtomAst` and the entity ASTs, the constraint vocabulary, the DSL/EDN serialization, `Metadata`. No graph algorithms, no perception, no ring enumeration.
- **`umol-graph`** holds `RingSet` / `RingEnumerator`, all `ops/*` (validators, resolvers, transformers, aromaticity perception, valence narrowing).
- **No chemist-facing wrapper exists.** Parsers return `MoleculeAst`; ops take `&MoleculeAst` and `&mut MoleculeAst`.

This is workable for the algebraic surface and is where progress happens at the moment.

### What the wrappers would address (when they return)

Current state: `MoleculeAst` holds a single-slot ring cache (`umol-ast/src/ast/molecule.rs:54`); other ground-only caches (`DistanceMatrix`, `BiconnectedComponents`, `MatchTarget`, `MorganTarget`, `Coordinate` annotations) have no home. Three concerns the wrappers would address:

1. **Additional ground-only caches.** Adding more cache slots directly to `MoleculeAst` means every transient pattern, mid-resolve partial structure, and serialization payload pays for slots it never uses. The wrapper would carry the new slots; the existing ring slot on the AST may stay where it is or move into the wrapper at the same time, an open call.
2. **The chemist's first-tier API** (the `parse_smiles → mol.morgan_fingerprint(2)` shape from §"API tiering") has no type to live on. Every call currently goes through raw `MoleculeAst`.
3. **The Pattern role** (matcher target with well-formedness checks and matcher-side scaffolding caches) has no home. Patterns are currently raw `MoleculeAst` with `SubPattern` constraints.

None of these is blocking now. (1) is fine until a second cache slot needs adding. (2) and (3) gate work that hasn't started (matcher port, Morgan, transformations).

### Why the previous attempt failed

Pre-restructure (per doc 92), `Molecule` / `Pattern` lived in `umol_graph::api` and held in-tree storage that overlapped with what the AST also stored. Operations migrated unprincipledly between tiers — a method might appear on `Molecule`, then later move to `MoleculeAst`, then later split back. The boundary was defined by "this needs caching" or "this is what users type" rather than by a structural rule about what kind of data the type owns. The result: bidirectional drift, duplicated state, no single source of truth.

### Candidate separating line (for when the question becomes actionable)

Sketch, not a plan:

- **`MoleculeAst`** owns *facts authored by the user or producer* — the algebraic content. Atoms, bonds, relations, their feature data, constraints, metadata. Equality and hashing are over this content. Parsers produce, serializers consume, transformations rewrite.
- **`Molecule`** owns *resolved interpretations* — facts plus the perception artifacts that interpret them under a specific model. `Arc<MoleculeInner>` with the `MoleculeAst` plus cache slots (added as the first consumer arrives, not speculatively per §"Cache slots on Molecule"). `is_ground()` required at construction.
- **`Pattern`** owns *query interpretations* — facts plus matcher-side scaffolding. `is_ground()` not required.

Pitfalls to avoid based on the previous attempt:

- **Cache-ownership rule: topology-only on the AST, attribute-derived on the wrapper.** Topology is invariant across resolution (per §"What's stable vs. changing during resolution"); a cache that depends only on `(atoms, bonds)` is safe to live on `MoleculeAst` because every in-place mutation (`atom_mut`, `bond_mut`, narrowing, lift/inline) leaves it valid, and structural edits go through `MoleculeBuilder` which produces a fresh AST anyway. The existing single-slot ring cache fits this rule. Caches that depend on attribute concretion (Morgan fingerprint, packed `MatchTarget`, anything that bakes in element / charge / order) become invalid the moment a single attribute narrows, and belong on the wrapper where ground-ness is enforced. This rule is the operational version of principle 3 in the taxonomy.
- **Caches opaque to AST equality.** `Molecule == Molecule` is `self.ast == other.ast`. The wrapper is sugar + cache; the AST is identity. (Already true today — `MoleculeAst::PartialEq` excludes the ring cache via custom impl.)
- **Placement rule = what data does the method touch.** If an operation only reads the algebraic content, it goes on `MoleculeAst` (or its views). Wrapper methods consume cached views (`mol.morgan_fingerprint(r)`, `mol.find_matches(&pattern)`); topology-only caches stay reachable via the AST (`mol.rings(...)` and any future topology-derived cache).

### Open questions to address when the time comes

Not for now; preserved so they aren't re-discovered later.

1. **Where the ring cache lives.** Default per the topology-vs-attribute rule above: stays on `MoleculeAst` (it's topology-only). `Molecule` carries no ring slot and forwards to `MoleculeAst::rings`. Revisit only if a concrete consumer surfaces that wants ground-only ring semantics distinct from raw-topology rings.
2. **Solver return type.** `Resolver::resolve(&mut MoleculeAst, &ChemistryModel)` currently mutates the AST in place and returns `Solution`. If `Molecule` returns, the resolver might switch to `(&MoleculeAst, ...) -> Molecule` and transfer the ring set it computed during aromaticity perception. This is the cache-transfer mechanic in §"Cache-transfer mechanics".
3. **Pattern wrapper timing.** `Pattern` carries no caches today (matcher is deleted). Likely deferred until matcher port lands.
4. **Transformation re-resolution.** Transformations take `&Molecule` and return `Molecule` per §"Transformation ops". Cache slots in the result start empty. Reaffirm vs. transferring caches across rewrites (doc 80 line 371: no transfer).
5. **`RingView` back-pointer if `RingSet` moves to `Molecule`.** Either `&Molecule` or `&MoleculeAst` would work; settled at that point.

## Invariants: essential vs model-dependent

Three tiers with different enforcement points.

### Tier 1: structural validity (enforced at `MoleculeAst::new`)

The minimum that lets the AST mean anything. Construction-time errors.

- **Index validity.** Every bond endpoint references an existing atom; every relation tuple references existing atoms. No dangling refs.
- **Relation arity.** Bond = 2 atoms; multicenter ≥ 3; aromatic system ≥ 1.
- **Element / isotope existence.** Type-checked via the `Element` / `Isotope` enums; you cannot construct a fictional element.
- **Bond order ≥ 0.** No upper cap (Cr–Cr quintuple bonds are real).

Holds for ground and partial ASTs alike. A pattern with `Undetermined` element still has valid topology.

### Tier 2: physics invariants (enforced as constraints, verified by `Solver::resolve`)

Universal — they follow from electron and angular-momentum conservation. Hold for every chemical system regardless of model.

- **Per-atom electron count.** Z − charge = bonded_electrons + 2·lone_pairs + unpaired.
- **Total charge consistency.** Sum of atom charges = asserted molecular total (if asserted via `Derived { TotalCharge }`).
- **Total spin coupling.** Per-atom unpaired electrons couple consistently to asserted total spin.
- **Unpaired electrons/spin multiplicity consistency.** = u % 2 + 1 <= M <= u + 1
- **Aromatic-system electron count.** Discovered count consistent with per-atom `aromatic_valence` values that produced it.

Expressed as `MoleculeConstraint::Derived` entries. Verified during propagation. A `Molecule` whose tier-2 constraints fail is a resolution error.

These checks may have been deleted along with `graph_ir::Molecule` / `MoleculeBuilder`. Restoration is mechanical; flag for the migration plan.

### Tier 3: model-dependent rules (NOT invariants — opt-in only)

Conventions that hold for "typical" organic chemistry but exclude valid compounds outside that scope. RDKit's sanitization fails exactly here: it rejects organometallics, hypervalent main-group compounds, multicenter-bonded systems, and anything outside the Daylight aromaticity model.

Explicitly **not** enforced:

- **Octet rule.** SF6, BF3, B2H6 violate it.
- **"Normal" valence tables.** Cr in (η6-C6H6)2Cr has non-standard bonded-electron count. Every transition-metal complex breaks valence-table assumptions.
- **Specific aromaticity model.** Daylight, MDL, Hückel, Clar all give different answers — the choice is configuration.
- **Connected-component constraint.** Salts and ion pairs are validly multi-component.
- **Charge / oxidation-state bounds.** Drug-like convention only.

These live in opt-in `validator` modules or as configurable Solver strategies. They never gate `Molecule::new`.

### The dividing line

> **Physics-required invariants are essential; chemistry conventions are not.**

Direct consequence of the general-chemistry scope. RDKit conflates tier 2 and tier 3 and rejects everything outside drug-like assumptions; we separate them deliberately.

### Where each tier lives

| Type | Tier 1 | Tier 2 | Tier 3 | `is_ground()` |
|---|---|---|---|---|
| `MoleculeAst` | enforced at `new` | not checked | not checked | not required |
| `Molecule` | inherited | enforced at `Solver::resolve` | not checked | required |
| `Pattern` | inherited | **not enforced** (patterns can be intentionally non-physical) | not checked | not required |

`Pattern` skipping tier 2 is deliberate: a SMARTS query "carbon with charge ≠ 0" is a valid pattern even though no specific charge is given. The matcher evaluates against a ground target where tier 2 already holds.

### Edge cases on the tier 2 / tier 3 boundary

- **Spin reachability.** Whether a given total spin is reachable from a set of unpaired electrons depends on the coupling scheme assumed (LS vs jj vs intermediate). The *consistency* of the count is tier 2; enforcing a *specific* coupling scheme is tier 3.
- **Aromatic electron-count rule.** "4n+2" is Hückel-specific. "Aromatic system has *some* integer electron count consistent with its atoms" is tier 2; "has *4n+2* electrons" is tier 3.

Default both to tier 3 unless a use case forces stricter checks.

### Ground-ness of variable-bearing AST forms

Decision (Phase 14 / 2026-05-15): pattern variables (`Ref`, `Bind`, free `Var` inside an `Expr`) are **never ground**, regardless of whether their candidate set or expression body has collapsed to a single value. The `is_ground` / `as_lit` surface treats binding constructs as carrying an open variable role even when the value space has been narrowed.

Concretely, applied consistently across all AST types:

| Type | Variable/binding form | Ground? |
|---|---|---|
| `ValueAst` | `Expr(e)` where `e` contains a free `Var` | no — `evaluate_checked(&Bindings::new())` returns `None` |
| `IsotopeAst` | `Expr(e)` with free `Var` | no — same slow path |
| `ImplicitHydrogensAst` | `Expr(e)` with free `Var` | no — same slow path |
| `ElementAst` | `Ref(_)` | no |
| `ElementAst` | `Bind { id, set }` (any set size, including singleton) | no |

Singleton non-binding sets *do* collapse: `ValueAst::LitSet([5])` and `ElementAst::Set([C])` are ground. The distinction is binding role (carries an `id` that downstream constraints may reference) vs. value enumeration (semantically a unary alternative).

This is the conservative position. The alternative — singleton `Bind` ⇒ ground, mirroring `LitSet → Lit` — remains defensible and can be revisited if a use case demands it.

## Deferred

- **Migration plan from current state to the proposed hierarchy.** To be sketched between doc 80 points 9 and 10. The current code has no `Molecule` type; `MoleculeAst` is the only public surface; parsers return `MoleculeAst`. The migration must order the `Molecule`/`Pattern` introduction, parser API conversion, and tier-2 invariant restoration.
- **Pattern cache contents.** Only ring view is needed for the immediate use cases; matcher-side scaffolding (per-atom constraint index, sub-pattern dependency graph, packed pattern adjacency) is addressed when step 9 lands.
- **`ReactionRule` / `ReactionRuleAst` parallel.** Mirror of the `Molecule` / `MoleculeAst` split for reactions. Addressed when doc 80 step 10 lands.

## Open questions

(To be filled in as the rest of this discussion progresses.)

## Implementation status (2026-04-27)

The pre-restructure tier-1 prototype (in-tree `umol_graph::api::Molecule`,
`MoleculePattern`, the bespoke `MoleculeAst`/`MoleculePattern` storage with
`OnceLock<RingSet>`, and `ResolverCell`) was deleted along with the
`umol_graph::ast`/`dsl`/`api` trees during the doc 92 restructure. What's
actually live now:

- [x] **`MoleculeAst`** in `umol-ast` — algebraic, `Arc`-wrapped per-relation
  storage; ground or partial; **carries a single-slot
  last-request-wins ring cache** (`MoleculeAst::rings(family, max_size)`,
  init/replace via `&mut self`). Topology-invariance during narrowing keeps
  the cache valid across in-place mutations; structural edits go through the
  builder, which produces a fresh `MoleculeAst` with an empty cache. Cache
  field excluded from `PartialEq` / `Eq` via a newtype wrapper.
- [x] **Tier-1 structural invariants** enforced at `MoleculeAst::new`.
- [x] **Per-atom electron-count tier-2 invariant.**
  `ElectronInvariantValidator` (`umol-graph/src/ops/validator.rs`) runs over
  every atom; uses local atom constraints first and falls back to topology
  via the AtomView chemistry-method pairs (scheme A1 — `bond_order_sum` /
  `valence_constraint`, `donated_pairs` / `donated_pairs_constraint`, etc.).
  Standalone-atom mode (`validate_atom`) reads constraints only.
- [x] **Entity-structure validator.** `EntityStructureValidator` checks
  `electrons.len() == atoms.len()` on `AromaticSystemAst` and
  `MulticenterBondAst`.
- [x] **Composite `Validator`** with sub-validators
  `ElectronInvariantValidator`, `SpinCouplingValidator` (stub),
  `ConstraintValidator` (stub), `EntityStructureValidator`. Each declares
  its own `Contradiction` and `Error` types; composite unions via `From`.
  Methods take `impl AsRef<MoleculeAst>` so `&MoleculeAst` and (future)
  `&Molecule` work interchangeably.
- [x] **Composite `Resolver` + `ChemistryModel`.** `ChemistryModel { valence,
  aromaticity }` — no `ResolveConfig` wrapper; one shared model. `Resolver`
  composes `ValenceResolver` (AtomTyping / Counts) and `AromaticityResolver`
  (HueckelRule / Hmo / Clar). Each algorithm carries one error enum;
  dispatcher classifies variants into `Solution::Contradictory` (chemistry),
  `Solution::Underdetermined` (e.g. HMO `UndeterminedAtom`), or `Err`
  (setup-level — currently only HMO `MissingParameters`).
- [x] **TableIR → MoleculeAst lift.** `IntoAst<MoleculeAst> for &Molecule`
  (and the per-atom and per-bond analogues) in
  `umol-graph/src/table_ir/lift.rs`. Same conversion vocabulary as
  umol-ast's DSL → AST lifts. `LiftError` empty for now; reserves a fail
  surface for future strict checks.
- [x] **Parser entry-points** return `MoleculeAst` (not `api::Molecule`):
  `parse_smiles`, `parse_smiles_with(... &ChemistryModel)`,
  `parse_smiles_to_ast`, plus the ctfile equivalents. Configuration takes
  `&ChemistryModel`.

## Outstanding

The chemist-facing tier-1 wrappers and several tier-2 invariants still need
to land:

- [ ] **`Molecule` chemist-facing wrapper.** Holds an `Arc<MoleculeInner>`
  with the AST plus chemistry-side cache slots. Doc 92's open question on
  metadata-on-`Molecule` (round-trip preservation vs. purely semantic
  object) is unsettled; revisit when the type is reintroduced.
- [ ] **Cache slots on `Molecule`.** `DistanceMatrix`,
  `BiconnectedComponents`, `MatchTarget`, `MorganTarget`. Add as their first
  consumer arrives, not speculatively. The `RingSet` cache already lives on
  `MoleculeAst` itself (doc 92 settled this); other Molecule-side caches are
  for ground-only views that don't make sense on patterns.
- [ ] **Coordinate annotations on `Molecule`.** Optional per-atom
  `Coordinate` payload propagated through MOL / CXSMILES roundtrip; stored
  but never recomputed.
- [ ] **`Pattern` chemist-facing wrapper** with matcher-side scaffolding
  (per-atom constraint index, sub-pattern dependency graph, recursion
  order, packed pattern adjacency) and well-formedness checks at
  construction. Lands alongside doc 80 step 9 (matcher port).
- [ ] **`ReactionRule` / `ReactionRuleAst` split.** Mirror of the
  `Molecule` / `MoleculeAst` split for reactions. Doc 80 step 10.
- [ ] **Tier-1 parser entry-points returning chemist-facing types.**
  `parse_smiles → Molecule`, `parse_smarts → Pattern`, `parse_smirks →
  ReactionRule`. Today the parsers return `MoleculeAst`; chemist-facing
  variants land once `Molecule` and `Pattern` exist.
- [ ] **Remaining tier-2 propagators.** `TotalCharge`, `TotalSpin`, and
  `AromaticElectronCount` exist as `MoleculeConstraint` variants but the
  `ConstraintValidator` is a stub. Wire evaluators in once their per-engine
  work lands.
- [ ] **Per-entity spin-coupling propagator.** `SpinCouplingValidator` is a
  stub returning `Determined`; check is `multiplicity = unpaired − 2k + 1`
  for some `k ∈ 0..=unpaired/2` on any entity carrying a `SpinStateAst`
  (atom, aromatic system, multicenter bond). Parity rule lives in
  `umol_shared::spin::SpinState`; lift into the validator alongside
  `ElectronInvariant`.
- [ ] **`ConstraintValidator` cross-checks.** Constraint-vs-topology
  agreement across entity types, plus molecule-scope constraints
  (`:connected`, etc.). Currently a stub returning `Determined`.
- [ ] **Matcher port** + the `evaluate.rs` it carried — both deleted at the
  end of doc 92. Doc 80 step 9.
- [ ] **Transformation ops** — `kekulize`, `aromatize`, `tautomers`,
  `to_canonical_smiles`, `apply_reaction` — with the signatures and result
  types from §"Transformation ops" and §"Result types".
- [ ] **Tier-3 model-dependent validators** (octet, normal-valence tables,
  drug-like charge bounds, connectedness). Opt-in `validator` modules;
  never gate construction.
- [ ] **Builder API** producing `MoleculeAst` then resolving to `Molecule`.
- [ ] **Resolution conformance suite port.** `tests/resolution/*` is gated
  behind `cargo` feature `legacy` and currently broken (uses deleted
  `umol_graph::ast`/`dsl`/old `Chemistry`/`ValenceTheory`). Either rewrite
  on the new `ChemistryModel` + `Resolver` API or delete and grow new
  conformance coverage. Same call applies to the gated benches:
  `morgan.rs`, `molecule_dsl_parsing.rs`, `substructure.rs`.
