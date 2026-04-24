# Molecule AST API

Working notes on the public API around `MoleculeAst` and its derived views, started while finishing the `Ring` → `RingView` refactor on `feature/relational-ast`. Focus: where derived data is cached, how it surfaces to consumers, and which procedural operations should also have constraint counterparts in step 9 of doc 80.

## Design context

Two anchors shape the design:

1. **General-chemistry scope.** umol targets organometallics, clusters, mixed-valence systems, and inorganic compounds alongside organics — see `project_general_chemistry_scope.md`. Cache decisions, the multicenter/dative/noncovalent relation types, and the invariant tiers (below) all reflect this scope rather than drug-like assumptions. RDKit-style "sanitization" rejects valid compounds outside that scope; we deliberately don't.

2. **Long-lived molecules.** Molecules enter once (parse/deser/transform), get queried and transformed many times. This justifies per-molecule cached views and influences the resolver's design — it can spend O(n) work per molecule that pays off across many queries.

## Concurrency

`Molecule` must be safe to hand to per-thread workers (doc 80 line 365). Requirements:

- `Molecule: Send + Sync`. The `Arc<MoleculeInner>` shape supports this.
- All cached views use `OnceLock`, not `Cell` / `RefCell`. Init-once, read-many.
- View contents (`RingSet`, `MatchTarget`, `MorganTarget`) must themselves be `Sync`. Today they hold owned data with no interior mutability.
- Cloning `Molecule` is an `Arc` bump; all clones share the same cache.

`Pattern` follows the same pattern. Coordinate annotations (below) also must be `Sync` — pure data, no mutation.

## Current state of the ring API

`umol-graph/src/ast/rings.rs` was reorganized to match the `views.rs` namespace pattern:

- `Ring` is now private storage (atoms vec + bonds vec).
- `RingView<'a>` is the public read shape: `pub idx: RingIdx`, plus `atoms()`/`bonds()` returning slices and `len()`/`is_empty()`/`shared_atoms()`/`shared_bonds()`.
- `RingSet` fields are private; accessors `family()`, `scope()`, `max_ring_size()`, `count()`, `ids()`, `iter() -> impl Iterator<Item = RingView<'_>>`, `get(idx) -> Option<RingView<'_>>`.
- Per-atom / per-bond accessors keep `atom_smallest_ring_size`, `bond_smallest_ring_size`, `contains_atom`, `contains_bond`.
- Pairwise relation surface: `relation`, `are_spiro`, `are_fused`, `are_bridged`, `spiro_neighbors`, `fused_neighbors`, `bridged_neighbors`, `fused_component(s)`, `shared_atoms`, `shared_bonds`, `graph()`.

`RingIndex` was renamed to `RingIdx` for consistency with the rest of the index types (`AtomIdx`, `BondIdx`, etc.). All consumers updated; 3608 lib tests pass.

`MoleculeAst` itself does not expose `rings()` — the only way to get a `RingSet` today is to construct a `RingEnumerator` and call `enumerate(&ast)`. That returns an owned `RingSet` and the caller is responsible for caching it. `Solver` currently builds one inside `resolve` and does not retain it across stratification.

## Where the `RingSet` cache should live

Doc 80 line 369 commits to: **`MoleculeAst` itself holds no cache.** Reasons restated there: interior mutability does not compose cleanly with equality, hashing, or constraint narrowing; ungrounded `MoleculeAst` values (patterns, mid-resolve partial structures) do not benefit from the same views.

The cache belongs on `GroundMolecule`, in the shape doc 80 already sketched at lines 340–361:

```rust
struct GroundMoleculeInner {
    ast: MoleculeAst,
    ring_set: OnceLock<RingSet>,
    distance_matrix: OnceLock<DistanceMatrix>,
    biconnected: OnceLock<BiconnectedComponents>,
    match_target: OnceLock<MatchTarget>,
    morgan_target: OnceLock<MorganTarget>,
}
```

Step 3 of the migration plan (line 406) deferred this — `GroundMolecule` is currently a plain newtype `(MoleculeAst)` with no `Arc`/`OnceLock`. Step 9 is the natural moment to lift it: the matcher recursion needs ring access in its inner loop, and a per-target cache is what keeps `Derived { InRing | RingSize | RingCount }` from re-enumerating per query.

## Two use cases driving the cache decision

### Aromaticity perception

Today: `Solver::resolve` allocates a fresh `RingSet` inside the aromaticity stratum and discards it. The valence → aromaticity → re-valence loop pays for at most one enumeration, so the cache is not load-bearing for *resolution*. It becomes load-bearing the moment a downstream pass (matcher, constraint verifier, fingerprinter) wants the same ring data without re-enumerating.

After step 9: aromaticity perception calls `ground.rings()` and the cache survives across the stratification *and* across subsequent matcher queries on the same molecule.

### SMARTS-style ring attributes

SMARTS distinguishes:

| SMARTS atom predicate | Meaning | Backed by |
|---|---|---|
| `R` | atom is in any ring | `RingSet::contains_atom` |
| `R<n>` | atom is in exactly *n* SSSR rings | `atom_to_rings[a].len()` (accessor missing today) |
| `r<n>` | smallest containing ring has size *n* | `RingSet::atom_smallest_ring_size` |

These are **topology-derived atom attributes** in the discussion 79 taxonomy: derived from graph structure, not stored on `AtomAst`. They are evaluated against a target during matching, not pinned on a query atom.

For a query batch over 10k molecules, hitting `ring_count_at` once per atom for every SMARTS check means the build cost must be amortized. `OnceLock` per `GroundMolecule` does that without a global cache.

One missing accessor: `ring_count_at(atom: AtomIdx) -> usize`. Trivial — `atom_to_rings.get(&atom).map_or(0, Vec::len)`.

## Procedural vs declarative split

Doc 80 line 194 commits to: *"Discovery is procedural; verification and narrowing are constraint-based."* Apply that rule to each `RingSet` method:

| Operation | Stays procedural | Becomes `DerivedPred` |
|---|---|---|
| `RingEnumerator::enumerate` | yes (graph search) | — |
| `iter` / `get` / `count` / `ids` | yes (cache walk) | — |
| `contains_atom`, `contains_bond` | — | `InRing` (already in D3) |
| `atom_smallest_ring_size` | — | `RingSize(ValueAst)` (already in D3) |
| `ring_count_at` (new accessor) | — | new `RingCount(ValueAst)` |
| `relation`, `are_spiro/fused/bridged` | yes (perception-internal) | — |
| `spiro_neighbors`, `fused_neighbors`, `bridged_neighbors` | yes (perception-internal) | — |
| `fused_component(s)` | yes (perception-internal) | — |
| `shared_atoms`, `shared_bonds` | yes (perception-internal) | — |

The asymmetry is principled: SMARTS attributes are atom- or bond-local properties derived from ring membership. Pairwise ring-vs-ring queries are aromaticity-internal — no SMARTS analog, no constraint slot. Adding constraint variants for `Spiro(a, b)` / `Fused(a, b)` would invent vocabulary nobody writes.

## Ring topology in patterns: no special constraint needed

A "back-closure" in SMILES (`C1CCCCC1`'s second `1`) is a notational convenience for the string form. After parsing it is an ordinary bond entry whose endpoints both already exist. `MoleculeAst` stores bonds in `Arc<Vec<BondAst>>` indexed off `graph: Graph`; nothing distinguishes a tree edge from a cycle-closing edge.

To assert "atom 0 lies in some 6-ring" as a sub-pattern:

```rust
SubPattern {
    anchor: 0,
    pattern: MoleculeAst {
        atoms: [C, C, C, C, C, C],
        bonds: [
            (0, 1), (1, 2), (2, 3), (3, 4), (4, 5),
            (5, 0),  // cycle-closing bond — same shape as any other
        ],
        ..
    }
}
```

The `(5, 0)` entry is what SMILES would call a back-closure. From the matcher's standpoint it is an ordinary bond constraint: VF2 must find six target atoms with all six bonds present, including `(t5, t0)`. The "ring" is encoded entirely in the pattern's bond list — no `RingPattern` constraint variant is needed.

So three ways to express ring membership in queries, picking the cheapest evaluator each time:

1. `R` / `r<n>` / `R<n>` → `DerivedPred::InRing | RingSize | RingCount` against the cached `RingSet`. O(1) bit/vec lookup.
2. "Some 6-ring touches this atom" → 6-cycle `SubPattern`. VF2 search.
3. "Atom is in the same ring as atom X" → `SubPattern` that anchors both atoms in a cycle.

(1) is for SMARTS atom predicates; (2)–(3) are for genuine topology assertions that go beyond SSSR membership.

## Implications for step 9 (SubPattern + matcher recursion)

Concrete additions on top of D3 in doc 80:

1. **`MatchTarget` carries ring access.** Either it wraps `GroundMolecule` (so `target.rings()` resolves through the cache) or it holds an `Arc<RingSet>` cloned at construction. The matcher dispatches `DerivedPred::InRing | NotInRing | RingSize | RingCount` against this view in a single hot-path call.

2. **Add `RingCount(ValueAst)` to `DerivedPred`.** Symmetric with `RingSize`; closes the SMARTS `R<n>` gap.

3. **Lift `GroundMolecule` to the `Arc<Inner>` + `OnceLock` shape from doc 80 lines 340–361.** Step 3 of the migration deferred this; step 9 needs it.

4. **Solver loop materializes views via `target.rings()`.** Per doc 80 lines 232–237, view materialization is step 1 of the propagate phase. Currently `Solver::resolve` builds rings inline; the migration unit is "stop allocating, call `ground.rings()`."

What does **not** change: `RingEnumerator` stays procedural, the pairwise relation surface stays consumer-local to aromaticity, and there is no `RingPattern` constraint variant.

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
MoleculeAst                       algebraic, immutable, ground or partial, no caches
   │
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

### Asymmetries worth naming

1. **Resolution flows into `Molecule` but not `Pattern`.** Patterns are intentionally underdetermined; running a resolver on them would either fail or invent constraints the user never wrote.
2. **Transformations produce `MoleculeAst`, not `Molecule`.** The producer wraps to `Molecule::new` if the result is ground (kekulization), or returns `MoleculeAst` if the result is itself a family or partial structure.
3. **A ground `MoleculeAst` can be wrapped as `Pattern`** — useful for "exact-match this molecule" queries. Cheap rewrap; cache slots start empty because pattern caches differ from molecule caches.
4. **`ResolverCell`'s topology cache transfer only targets `Molecule`.** The resolver is not called on patterns.

### Cache-transfer mechanics

`ResolverCell` owns a `OnceLock<RingSet>` that the aromaticity stratum populates. On finalize, the cell moves the `OnceLock` value into `MoleculeInner`'s slot. No extra copy; `RingSet` is moved by value.

Topology caches not populated during resolution (e.g., a registry-only valence pass that never touched rings) leave the corresponding `OnceLock` in the resulting `Molecule` unset. The first chemistry query pays for it. Pay-on-use is preserved.

Solver per-molecule state (candidate sets, propagation queue, fixpoint counters) lives only in `ResolverCell` and is dropped at finalize. None of it leaks into `Molecule`.

### Naming

`GroundMolecule` is precise but reads like compiler vocabulary to chemists. With the recently-deleted `graph_ir::Molecule` gone, the name `Molecule` is free. Recommended pairing:

- `MoleculeAst` — algebraic, the input/output of parsers, transformations, and serializers.
- `Molecule` — chemist-facing, ground invariant, cached views.

Symmetric with the rest of the AST-vs-resolved distinction in the codebase. Final decision deferred to the implementation step.

## Operation boundaries: what each IO / transformation operates on

### DSL IO

EDN parser produces `MoleculeAst`. Caller wraps via `Molecule::new(ast)` (ground required, errors if not) or `Pattern::new(ast)` (no ground requirement, well-formedness checks only). `is_ground()` is paid at the wrap.

### SMILES / MOL IO

Read: `SMILES/MOL → TableIR → MoleculeAst → Solver::resolve → Molecule`.

TableIR → MoleculeAst is the clean handoff. The lowering does not know whether the result is ground; the resolver decides. Skipping the intermediate and going `TableIR → ResolverCell` directly would couple file-format parser with solver strategy.

If `Solver::resolve` receives an already-ground input (canonical SMILES with all atoms fully specified), it detects this up front and returns `Molecule::from_inner` with empty caches — same return type, zero cost.

Write: `Molecule → MoleculeAst (via Molecule::as_ast) → TableIR → text`. No separate "table view" on `Molecule` needed; the AST is the canonical serialization shape, and after resolution all fields are concrete. Writers that need derived data (e.g., CTAB ring block) query the `Molecule` cache.

### Transformation ops (kekulize, aromatize, tautomerize)

These are rewrites that preserve ground-ness by construction. The resolver is not needed on the output — the operation is responsible for producing a valid ground structure.

| Op | Signature | Resolver pass? |
|---|---|---|
| `kekulize` | `&Molecule → Molecule` (pick one canonical) | no |
| `kekulize_all` | `&Molecule → Vec<Molecule>` | no |
| `aromatize` | `&Molecule → Molecule` | no — perception is part of the op |
| `tautomers` | `&Molecule → Vec<Molecule>` | no — proton shifts preserve valence by construction |

Start with `Vec<Molecule>` for families. E-graph for tautomers is deferred (doc 80 line 270).

Cache implication: per doc 80 line 371, a new `Molecule` from a rewrite starts with empty `OnceLock` caches. The cost of proving a transferred cache is still correct outweighs the rebuild.

### Reactions

Pattern well-formedness is validatable at `Pattern::new` — sub-pattern anchors in range, constraints reference existing relations, no cyclic sub-pattern references. Standard AST-level structural checks.

Semantic well-formedness of a full `ReactionRule` (LHS + RHS + correspondence) — whether applying the rule to *any* matching target yields a chemically valid product — is **not** validatable statically in the general case. It depends on target-specific valence, charges, stereo. A narrow class (mass-balanced electron-pushing rules) is provably ground-preserving; the general case is not.

Flow:

```
rule: &ReactionRule + target: &Molecule
    │
    ├─► LHS Pattern matches target → Assignment
    ├─► Rule applied mechanically → MoleculeAst (possibly non-ground on RHS-introduced atoms)
    └─► Solver::resolve on the result → Molecule
```

Re-resolution cost is amortized over VF2 match + rewrite. Regiochemistry: one rule + one target can produce multiple products via different match assignments.

### Coordinates as annotations

3D and 2D coordinates are owned by `umol-geometric`, not `umol-graph`. The `Molecule` type carries coordinates only as **annotations**, not as first-class chemistry data:

- Source: MOL files and CXSMILES carry coordinates; the parser must propagate them.
- Storage on `Molecule`: an optional per-atom `Coordinate` payload, populated only if the input had coordinates.
- Operations: pass-through only. `Molecule` does not interpret, normalize, or recompute coordinates. Round-trip (MOL → Molecule → MOL) preserves them faithfully.
- No conformers. Multiple-conformer support is a `umol-graph → umol-geometric` conversion via distance geometry (doc 71), not a `Molecule` concern.

This keeps `umol-graph` graph-only at its core while letting MOL/CXSMILES roundtrip through.

### Result types

No generic `OpResult`. Metadata shapes are heterogeneous:

| Op class | Returns | Rationale |
|---|---|---|
| `kekulize`, `aromatize` | `Molecule` | Single product, no meaningful metadata |
| `to_canonical_smiles` | `String` | Query, not a transformation |
| `kekulize_all`, `tautomers` | `Vec<Molecule>` | Single-type family, no correspondence |
| `apply_reaction` | `Vec<ReactionResult>` | Per-match metadata (assignment, atom mapping, rule ref) must travel with the product |
| `match_substructure` | `Vec<Assignment>` | Matches are not new molecules |

`ReactionResult` carries products + atom mapping + rule reference + assignment. A generic `OpResult<M>` parameterized by metadata is premature; wait for a second op with the same shape before abstracting.

## Summary table

| Stage | Works on | Notes |
|---|---|---|
| EDN / DSL parse | produces `MoleculeAst` | caller wraps to `Molecule` or `Pattern` |
| SMILES / MOL read | `TableIR → MoleculeAst → Solver → Molecule` | clean boundary at AST |
| SMILES / MOL write | `Molecule → MoleculeAst → TableIR → text` | no extra table view |
| Kekulize / aromatize | `&Molecule → Molecule` or `Vec<Molecule>` | no resolver pass |
| Tautomerize | `&Molecule → Vec<Molecule>` | E-graph later |
| React | `&ReactionRule + &Molecule → Vec<ReactionResult>` | re-resolve product internally |
| Canonicalize | `&Molecule → String` | query, not transformation |

## API tiering: where MoleculeAst surfaces

Walking the four 3-month use cases through a candidate first-tier API confirms that `MoleculeAst` does not appear in the user-visible code path of any of them:

| Use case | First-tier code |
|---|---|
| SMILES → Morgan | `let mol = parse_smiles(s)?; let fp = mol.morgan_fingerprint(2);` |
| DSL → Morgan | `let mol = Molecule::from_edn_str(text)?; let fp = mol.morgan_fingerprint(2);` |
| SMILES + SMARTS → annotated DSL | `parse_smiles → Molecule`, `parse_smarts → Pattern`, `mol.find_matches(&pattern)`, `mol.to_edn_with_match(m)` |
| SMILES + SMIRKS → product SMILES | `parse_smirks → ReactionRule`, `apply_reaction(&rule, &mol) → Vec<ReactionResult>`, each `result.products: Vec<Molecule>` |

`MoleculeAst` is therefore second-tier: visible in the public API but not the type a normal user reaches for first.

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

`parse_smiles` uses `Solver::default()`. Power users compose the raw parser with a custom `Solver`.

### Implication for naming

The two-tier framing supports:
- `MoleculeAst` — algebraic, second-tier; `Ast` suffix makes the role explicit.
- `Molecule` — chemist-facing, first-tier; the name a chemist expects.
- `Pattern` — query-facing, first-tier.

`MoleculeAst` is not hidden, just not what users type day-to-day.

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

## Deferred

- **Migration plan from current state to the proposed hierarchy.** To be sketched between doc 80 points 9 and 10. The current code has no `Molecule` type; `MoleculeAst` is the only public surface; parsers return `MoleculeAst`. The migration must order the `Molecule`/`Pattern` introduction, parser API conversion, and tier-2 invariant restoration.
- **Pattern cache contents.** Only ring view is needed for the immediate use cases; matcher-side scaffolding (per-atom constraint index, sub-pattern dependency graph, packed pattern adjacency) is addressed when step 9 lands.
- **`ReactionRule` / `ReactionRuleAst` parallel.** Mirror of the `Molecule` / `MoleculeAst` split for reactions. Addressed when doc 80 step 10 lands.

## Open questions

(To be filled in as the rest of this discussion progresses.)

## Implementation status (2026-04-17)

- [x] **`MoleculeAst`** (`ast/molecule.rs`) — algebraic, immutable, `Arc`-wrapped per-relation storage; ground or partial; no caches.
- [x] **`Molecule`** (`api/molecule.rs`) — `Arc<MoleculeInner>`; ground invariant enforced at `Molecule::new`; `OnceLock<RingSet>` cache; view API mirroring `MoleculeAst` (atoms, bonds, dative/noncovalent/aromatic/multicenter relations, neighbors, graph, `bond_order_sum`, `dative_bond_order_sums`, `is_in_aromatic_system`); EDN roundtrip via `to_edn_str` / `from_edn_str`.
- [x] **`ResolverCell` cache transfer.** Topology cache (`OnceLock<RingSet>`) transfers from `ResolverCell` into the resulting `Molecule` on finalize.
- [x] **Tier-1 structural invariants** enforced at `MoleculeAst::new` (index validity, relation arity, element/isotope existence, bond order ≥ 0).
- [x] **Per-atom electron-count tier-2 invariant.** `ElectronInvariant` propagator (`unify/propagate.rs`) runs on every atom in `Validator::validate` and the matcher post-filter. Theory-independent; equation per the orbital-side / source-side balance.

## Outstanding

- [ ] **`Pattern` type with matcher caches.** Today `MoleculePattern` (`api/pattern.rs`) is a thin `Arc`-wrapped AST. Doc 86 calls for a long-lived `Pattern` carrying matcher-side scaffolding (per-atom constraint index, sub-pattern dependency graph, recursion order, packed pattern adjacency) and `Pattern::new` well-formedness checks. Lands with doc 80 step 9.
- [ ] **Tier-1 parser entry-points.** `parse_smiles → Molecule`, `parse_smarts → Pattern`, `parse_smirks → ReactionRule`, plus tier-2 `*_to_ast` and configurable `parse_*_with` variants. Current parsers return `MoleculeAst`.
- [ ] **Remaining tier-2 invariants.** `TotalCharge`, `TotalSpin`, and `AromaticElectronCount` exist as `MoleculeConstraint` variants but have no evaluator. Once doc 87's propagator evaluators land, wire them into `Validator::validate`.
- [ ] **Per-entity spin-coupling invariant.** `SpinCouplingInvariant` stub exists in `ops/propagate.rs`; check is `multiplicity = unpaired − 2k + 1` for some `k ∈ 0..=unpaired/2` on any entity carrying a `SpinStateAst` (atom, aromatic system, multicenter bond). Parser's tier-2 leak was removed (`dsl/predicates.rs::apply_spin_pair` no longer validates; `SpinStateAst::validate` and the matching `ParseError` variants are gone). Wire the propagator into `Validator::validate` + matcher post-filter alongside `ElectronInvariant`.
- [ ] **Additional `Molecule` cache slots.** `DistanceMatrix`, `BiconnectedComponents`, `MatchTarget`, `MorganTarget` — add as their first consumer arrives, not speculatively.
- [ ] **`Pattern` cache slots.** Per-atom constraint index, sub-pattern dependency graph, packed pattern adjacency. Land with step 9.
- [ ] **`ReactionRule` / `ReactionRuleAst`.** Mirror of the `Molecule` / `MoleculeAst` split for reactions. Doc 80 step 10.
- [ ] **Coordinate annotations on `Molecule`.** Optional per-atom `Coordinate` payload propagated through MOL / CXSMILES roundtrip; `Molecule` stores but never recomputes.
- [ ] **Transformations as ops** — `kekulize`, `aromatize`, `tautomers`, `to_canonical_smiles`, `apply_reaction` — with the signatures and result types from §"Transformation ops" and §"Result types". `kekulize_all` and `tautomers` return `Vec<Molecule>`; `apply_reaction` returns `Vec<ReactionResult>`.
- [ ] **Tier-3 model-dependent validators** (octet, normal-valence tables, drug-like charge bounds, connectedness). Opt-in `validator` modules; never gate `Molecule::new`.
- [ ] **Builder API** producing `MoleculeAst` then resolving to `Molecule` (tier-1 surface example in §"Where AST surfaces even at tier 1").
