# Unified constraint-based molecule AST

Date: 2026-04-10

## Context

The kekulization algorithm needed test inputs that were painful to construct in Rust. That pain motivated a better DSL. The DSL design converged on homoiconicity — ground terms and patterns share notation — which required a format with first-class tagged literals and nesting, which led to EDN. While finishing the resolution framework (valence, atom typing, aromaticity; stereo not yet started), it became clear that each resolution pass is just adding constraints to a partial representation and discharging them. The mechanism the DSL uses to express patterns — constraints over a multi-relational schema — is what the resolver needs.

Pattern matching and resolution are the same engine. This doc records the decisions that flow from recognizing that identity.

Discussion thread: docs 16 (hypergraphs), 42 (relational-molecular-structure-representation), 60 (molecule-builder-dsl), 79 (pattern-language-design).

## Central commitment

**Unified constraint solving as the destination.** Pattern matching, rule application, valence resolution, atom typing, aromaticity perception, and any future reasoning over molecular structure all run through one engine that consumes `MoleculeConstraint` values and produces narrowed `MoleculeAst` values. Two engines — one for "patterns" and one for "resolution" — is the status quo; unifying now is cheaper than unifying later. Fail-fast stance: commit now, benchmark early, back out only if performance is catastrophic.

**`MoleculeAst` is the only molecule type.** No separate `Molecule`, `MoleculePattern`, or `MoleculeBuilder` storage. Homoiconicity expressed in code: whatever you parse, query, resolve, or emit is a `MoleculeAst`. A ground molecule is a `MoleculeAst` whose attribute fields are all concrete and whose constraints vec contains only authored assertions (and for the fully-resolved case, none at all).

These two commitments only work together. Separate resolution from matching and you need two types for storage (pattern vs. ground). Unify storage but split the engine and you have two pieces of code doing the same thing. The refactor is all-or-nothing.

## Decisions summary

| ID | Decision | Resolution |
|---|---|---|
| D1 | Unified constraint solving? | Yes. Commit now; fail fast if performance fails. |
| D2 | `MoleculeAst` as THE type? | Yes. `GroundMolecule` newtype for ground-only APIs. |
| D3 | `MoleculeConstraint` starter set? | `SubPattern`, `Derived`, `Matcher`, `And`/`Or`/`Not`. |
| D4 | Atom identity in constraints? | `usize` (wrapped `AtomIdx`). Labels stay in parsers. |
| D5 | Sub-patterns: where? | In `constraints` vec. Bindings scope lexically inward. |
| D6 | Injective matching default? | Yes. Non-injective via `MatcherFlag`. |
| D7 | Reaction rules as separate top-level? | Yes. `ReactionRuleAst { lhs, rhs, guards }`. |
| D8 | Resolution emits constraints? | Yes. Each pass maps to a set of `MoleculeConstraint` emissions. |
| D9 | Solver backend architecture? | Dispatch by variant to specialized backends. Not a generic CSP engine. |

## Integration with existing code

Several pieces of this design collide with existing symbols. Flagging them here; the nomenclature pass is deferred per the general preference for non-abbreviated names, and can be done in a cleanup PR once the constraint machinery is in place.

- **`AtomIndex` / `BondIndex`** already exist in `graph_ir::molecule` as `pub type AtomIndex = NodeIndex<u32>;` / `pub type BondIndex = EdgeIndex<u32>;`. They are petgraph handles tied to a specific `StableGraph` instance and cannot be directly repurposed as indices into `MoleculeAst::atoms` without a translation layer. The constraint-side index is a distinct concept and is currently a bare `usize` indexing into `MoleculeAst::atoms` (and the parallel relation vecs). Long-term resolution: keep the petgraph aliases local to the view layer (rename to something like `GraphNodeIdx` / `GraphEdgeIdx` once specialized views land) and use `AtomIndex` / `BondIndex` as the AST-level newtypes wrapping the current `usize`.

- **`ValueAst`** already exists in `umol_shared` with `Bindings = HashMap<String, i64>`, `Expr`, `ArithOp`, `RelOp`, and `capture` / `evaluate` / `evaluate_bool`. It is the scalar constraint language used by atom and bond specs in the DSL. **Reuse, do not reinvent.** The `Derived` constraints that compare integer attributes (valence sum, total charge, ring size) accept a `ValueAst` rather than a fresh `(CompareOp, i64)` pair. `ValueAst::capture` already implements the unification-on-match behavior that D5 needs for sub-pattern bindings, so the sub-pattern scoping story and the scalar-guard story share one substrate. Spin is the exception — `TotalSpin` carries `SpinStateAst`, not `ValueAst`, because the multiplicity and unpaired-electron fields are coupled and must be expressible as a single ground value or capture pair.

- **`TableIR`** is the real input boundary for SMILES and MOL. Parsers produce `TableIR::Molecule`, which is currently lowered via `MoleculeBuilder::from_table_molecule`. After the migration, the lowering target is `MoleculeAst`. The DSL path is a second source that also lands in `MoleculeAst`. D8's parse step below reflects both paths.

- **`resolve_valence_with`** (`graph_ir::resolution`) is the current reference for valence resolution. It dispatches `ValenceMatcher::candidates_for(...)` per atom for one of two strategies: `AtomTyping { registry }` (registry-based spec matching) or `Counts { table, allow_implicit_hydrogens }` (RDKit-style element-parametrized constraints). Both strategies stay; the migration rewrites only the dispatch layer to emit constraints instead of calling `set_atom_candidates` directly.

- **`resolve_aromaticity_with`** is the current reference for aromaticity perception. Three algorithms are implemented and selected via `AromaticityStrategy`: Hückel rule, Clar sextet, Hückel MO theory. Each picks a `RingFamily` (`Simple` or `InducedBenzenoid`) and runs `AromaticityModel::aromatic_systems(...)`. These stay as the *discover* step in the discover-then-verify pattern described in D8. The verification step is what becomes a constraint pass.

## D3 — Constraint enum starter set

```rust
pub enum MoleculeConstraint {
    SubPattern { anchor: usize, pattern: Box<MoleculeAst> },
    Derived { predicate: DerivedPred, refs: RelationRefs },
    Matcher(MatcherFlag),
    And(Vec<MoleculeConstraint>),
    Or(Vec<MoleculeConstraint>),
    Not(Box<MoleculeConstraint>),
}

pub struct RelationRefs {
    pub atoms: Vec<usize>,
    pub bonds: Vec<usize>,
    pub dative_bonds: Vec<usize>,
    pub aromatic_systems: Vec<usize>,
    pub multicenter_bonds: Vec<usize>,
    pub noncovalent_bonds: Vec<usize>,
}

pub enum DerivedPred {
    // Scalar predicates — accept ValueAst for match/bind/guard semantics.
    // ValueAst::Lit(n) means "=n"; ValueAst::Expr(rel) means "satisfies rel"
    // with variables bound to the computed scalar. Unification with outer
    // bindings is handled by ValueAst::capture.
    TotalCharge(ValueAst),          // sum over atom charges
    TotalSpin(SpinStateAst),        // from per-atom spin coupling
    ValenceSum(ValueAst),           // σ-bond-order sum at an atom
    AromaticElectronCount(ValueAst),// Hückel count on a ring system
    RingSize(ValueAst),             // smallest containing ring size

    // Structural — boolean, backed by cached views on the target.
    InRing,
    NotInRing,
    InRelation(RelationSym),
    NotInRelation(RelationSym),
}

pub enum MatcherFlag {
    Injective,      // default, not emitted explicitly
    NonInjective,   // allow collapse (graph homomorphism mode)
    Induced,        // no extra edges between matched atoms
}
```

`Derived` carries `RelationRefs` rather than a single `atoms: Vec<usize>` so a predicate can scope over any sort (e.g., `AromaticElectronCount` over an aromatic-system tuple, `RingSize` over a ring closure encoded as a multicenter set). Molecule-wide aggregates like `TotalCharge` and `TotalSpin` use `RelationRefs::default()` (all sorts empty) by convention — the predicate is its own scope hint.

`TotalSpin` carries a `SpinStateAst` rather than `ValueAst` because the spin state has two coupled fields (unpaired electrons and multiplicity) and the round-trip surface for `:spin` on `MoleculeAst` already speaks `SpinStateAst`. `SpinStateAst` admits `Wildcard`, `Lit(SpinState)`, and `Pair { unpaired, multiplicity }` forms — `Pair` is what allows `?u` and `?m` capture variables.

`And`, `Or`, and `Not` are pure structural combinators on `MoleculeConstraint`. They are recursive (like `SubPattern`'s `Box<MoleculeAst>`) and have no semantic ordering. The relational matcher will discharge them by recursing — there is no separate Boolean-CSP backend.

`RelationSym` is an enum identifying one of the six core relations (`Atoms`, `Bonds`, `DativeBonds`, `AromaticSystems`, `MulticenterBonds`, `NoncovalentBonds`).

**Deferred (no slots yet):** `Path` (Kleene), external-solver constraints, probabilistic constraints, registry-backed `AtomTypeCandidate` (will live in solver state, not the AST). Add when a concrete use case arrives; do not pre-build.

## D5 — Variable bindings and sub-patterns

Sub-patterns inherit the parent's lexical scope. Matcher threads a binding environment through the recursion.

```rust
struct MatchEnv {
    bindings: HashMap<VarId, Value>,
}

fn match_pattern(
    query: &MoleculeAst,
    target: &MatchTarget,
    env: &mut MatchEnv,
) -> Option<Assignment>
```

**Scoping rules:**

- **Outer bindings are visible inside sub-patterns.** If `?h` is bound in the parent, the sub-pattern's value-DSL can reference `?h` on any attribute slot.
- **Sub-pattern bindings do not escape upward.** Bindings introduced inside a sub-pattern are discarded when the sub-pattern returns. Recursive matcher calls pass a borrowed-mut env that is snapshot-and-restored at the sub-pattern boundary.
- **Re-binding unifies.** If the sub-pattern's value-DSL binds `?h` where `?h` is already bound in the parent, the two must agree or the sub-pattern fails to match. No shadowing, no name-collision error — unification.
- **Anchor identity is implicit.** `SubPattern { anchor, pattern }` means: atom 0 of the sub-pattern is pinned to `anchor`. The matcher starts the recursive call with that assignment fixed; the sub-pattern's atom 0 does not re-match.

This matches how SMARTS recursive patterns work conceptually (local check on an atom's neighborhood, inheriting outer context) and how Datalog handles shared variables across rules.

**Example** — `[C;$(C(=O)O)]` (a carbon that is the carbonyl C of a carboxylic acid):

```edn
{:atoms [#atom "C"]
 :constraints [{:sub-pattern
                {:anchor 0
                 :pattern {:atoms [#atom "C" #atom "O" #atom "O"]
                           :bonds [[0 1 :double] [0 2 :single]]}}}]}
```

Sub-pattern atom 0 is the same as parent atom 0 (the anchor); atoms 1 and 2 are new.

**What is NOT supported in the first version:**
- Sub-patterns that reference parent atoms other than the anchor (multi-anchored sub-patterns).
- Sub-pattern results flowing back out (e.g., "bind the matched O of the sub-pattern and reference it in the outer pattern").

Both would be additive later. SMARTS does not have them either.

## D8 — Constraint flow through resolution

Each current resolution pass becomes a constraint producer. Sketch of what each pass emits and how the solver discharges it.

### Parse / lower

Two input paths:
- **SMILES / MOL** → `TableIR::Molecule` (existing parsers, unchanged) → `MoleculeAst`. The `TableIR` → `MoleculeAst` lowering replaces the current `MoleculeBuilder::from_table_molecule`.
- **DSL** (EDN-based) → `MoleculeAst` directly via `FromEdn`.

Output: `MoleculeAst` with attribute fields set where known, `None` elsewhere. Constraints vec populated with any assertions from the input (`TotalCharge` / `TotalSpin` from top-level `:charge` / `:spin`, explicit `SubPattern`s, user guards).

### Valence resolution

Input: `MoleculeAst` where some atoms have underdetermined `(unpaired, lone_pairs)`.

Current reference: `resolve_valence_with`, which already handles both `AtomTyping { registry }` and `Counts { table, allow_implicit_hydrogens }` strategies. The existing `ValenceMatcher::candidates_for` is reused as a backend.

Emits, per atom:
- `Derived { ValenceSum(ValueAst::Lit(v)), [a] }` where `v` is the σ-bond-order sum to non-H neighbors (computable from the `bonds` vec).
- A finite-domain narrowing constraint over `(unpaired, lone_pairs)` from the registry or counts table (stored in the solver state, not the AST).

Solver discharges: filters per-atom candidate sets against the electron invariant and valence sum. If one candidate survives, atom fields concrete. If multiple survive, Hund's rule selects or the candidates persist as `Family<Sum>`.

### Atom typing

Input: `MoleculeAst` with some atoms narrowed to candidate sets.

Emits, per atom: `Derived { AtomTypeCandidate(registry_ref), [a] }`.

Solver discharges: registry lookup reduces the candidate set to registry-valid specs. Narrowing cascades with valence resolution — they are actually one pass in the solver because the constraints share state.

### Aromaticity perception

Input: `MoleculeAst` with resolved atoms.

Current reference: `resolve_aromaticity_with`, which already implements three algorithms — Hückel rule, Clar sextet, Hückel MO theory — selected via `AromaticityStrategy`. Each picks a `RingFamily` (`Simple` or `InducedBenzenoid`) and runs `AromaticityModel::aromatic_systems`. All three stay.

This is a discovery step, not just a check. The pass runs graph analysis (`RingEnumerator` over the selected ring family) to PROPOSE aromatic ring sets, then emits:
- New `aromatic_systems` tuples (one per candidate ring system, populated in the core relations vec).
- `Derived { AromaticElectronCount(ValueAst::Lit(n)), ring_atoms }` for each proposed system, with `n` the Hückel / Clar / HMO count produced by the chosen algorithm.

Per-atom aromatic valence is written to `AtomAst` directly rather than emitted as a constraint — aromaticity perception only emits constraints for things that span tuples.

Solver discharges: verifies the electron count invariant on each proposed system (using the same algorithm that proposed it, in verification mode); discards failing proposals by removing the tuple from `aromatic_systems`; propagates accepted aromatic valence to per-atom `AtomAst` fields.

The graph analysis stage is NOT in the constraint engine — it is a pre-pass that uses petgraph and then emits constraints for verification. **Discovery is procedural; verification and narrowing are constraint-based.** This is the general pattern: any pass that needs to SEARCH the graph runs as a procedural step and emits constraints as output.

### Stereochemistry (future)

Will emit `Derived { Chirality(a, label), [a] }` constraints computed from CIP priority rules over the graph and attributes. Solver validates; does not invent new stereo from thin air. Not implemented; not in the three-month window.

### Final state

A "resolved" `MoleculeAst` has:
- All attribute fields concrete (or wrapped in `Family<Sum>` for genuine ambiguity).
- Constraints vec empty, or containing only user-authored assertions that passed validation (e.g., a `:charge 0` that the user wrote and the solver verified).
- Invariants checked.

There is no separate "build" step. The resolved state IS a `MoleculeAst`.

## D9 — Solver backend architecture

Not a generic CSP engine. Dispatch by constraint variant to a specialized backend.

| Constraint kind | Backend | Stage |
|---|---|---|
| Per-tuple attribute predicates | Relational matcher (VF2-style) | match |
| Core-relation tuple presence | Relational matcher | match |
| `SubPattern` | Recursive matcher call | match |
| `Derived { ValenceSum, ... }` | Finite-domain propagator | propagate |
| `Derived { AtomTypeCandidate }` | Registry filter | propagate |
| `Derived { Aromatic* }` | Hückel verifier | verify |
| `Derived { TotalCharge / TotalSpin }` | Arithmetic over tuples | verify |
| `And / Or / Not` | Recursive matcher dispatch | match |
| `Derived { InRing, RingSize }` | Cached view lookup | verify |
| `Derived { Distance }` | Cached distance matrix | verify |
| `Matcher(...)` | Configures backtracker | match-setup |
| Future: optimization constraint | Lagrangian / ILP | optimize |
| Future: probabilistic | Belief propagation | distribute |

### Solver loop

```
1. Materialize views on the target molecule:
   - ring membership bitset
   - distance matrix (on demand)
   - biconnected components
   - atom degree array
   - cached per-relation hash indexes

2. Propagate: iterate over constraints, narrowing attribute fields
   and candidate sets. Fixpoint when no change.

3. Decide: for remaining ambiguity, apply a preference rule
   (Hund, low-spin default, canonical Kekulé) or wrap as Family<Sum>.

4. Validate: check that all discharged constraints are actually
   satisfied. Unsatisfied = error.
```

Passes 2–4 loop if decisions unblock further propagation. Terminate when no change or a configured iteration cap is hit.

### Why not generic

A generic CSP frontend would force bond perception's Lagrangian dual and pattern matching's backtracking search through the same interface. They have nothing in common at the algorithm level. Shared frontend costs more than it saves. Better: constraint variants are a dispatch table; each backend owns its variants and runs in its own pass. The frontend is the `MoleculeConstraint` enum; the backend is a `match` expression inside the solver.

### Views are the performance trick

Derived views (ring bitset, distance matrix, biconnected components, atom degree array) are computed once per target molecule and cached. Every `Derived { InRing, [a] }` constraint becomes a bitset lookup. Every `Derived { Distance(...), [a, b] }` becomes a matrix access. This is what makes the constraint engine competitive with RDKit's hand-rolled algorithms.

### Incrementality

Deferred. Each molecule is resolved once and queried many times; incremental constraint addition is rare in practice. If it becomes necessary, view caches invalidate and affected constraints re-run. Not in the three-month window.

## Rewrites: beyond match-and-resolve

Not every molecule-level operation is a constraint-satisfaction problem. Constraint solving narrows attributes and candidate sets on a fixed topology; **rewrites** change the topology or the representation itself and produce new molecules (or families of molecules) from old ones.

Operations that are rewrites, not resolutions:

- **Kekulization / aromatization.** Kekulization replaces an `aromatic_systems` tuple with a specific cycle cover of localized double bonds (attribute-level change plus relation-tuple removal). Aromatization is the inverse. Both preserve atom identity, so DPO machinery is overkill — an in-place transform over the `bonds` vec and `aromatic_systems` vec suffices. The constraint engine can verify the output (e.g., the resulting Kekulé structure satisfies each atom's valence invariant) but does not drive the rewrite itself.
- **Tautomer enumeration.** Produces a family of molecules related by proton and π-bond migration. E-graph saturation (`egg`-style) is a genuinely interesting fit: tautomers share most structure and an e-graph avoids re-deriving the same local transforms. Worth prototyping after the migration lands; not in the three-month window.
- **Reactions.** Full DPO graph rewriting with LHS → RHS rules. `ReactionRuleAst` is the top-level type (D7). Rule application produces new `MoleculeAst` values. The matcher from the constraint engine is reused for the LHS-finding step; the rewrite step is separate procedural code.

The implication for the doc: **do not try to force rewriting into the constraint variant enum.** The solver engine (D9) handles match + narrow + verify. Rewrites call the matcher for their LHS step and then run specialized rewrite code. Kekulization and aromatization are the rewrites that must work in the three-month window; tautomers and reactions are scheduled after.

## Performance: matching RDKit

User-stated concern: if this is principled but 100× slower than RDKit, it is not useful.

### Honest analysis

`MoleculeAst` with `Option` fields and an almost-always-empty constraints vec pays a small per-access tax vs. RDKit's dense `ROMol`:
- `Option<T>` adds a discriminant byte (cache-line impact at molecule scale: negligible).
- Attribute access is `unwrap()` (branch, trivially predicted after resolution).
- Pattern matching on `ElementExpr::Concrete(e)` vs. direct `e` (one extra branch).

Estimated overhead on tight attribute-access loops: 1.1–1.5×. **Not 100×**, but not free.

### Where RDKit's speed comes from

- Specialized per-algorithm data structures (`ROMol` has packed neighbor lists, precomputed rings, stashed atom type indices).
- 25 years of micro-optimization on the hot paths (substructure match, fingerprints, SSSR).
- No constraint-engine overhead — everything is direct C++ method calls.

### The strategy: specialized views

Hot loops do not run on `MoleculeAst` directly. They run on specialized views built once per target:

```
MoleculeAst (canonical)
    ↓ build once, amortize across many queries
MatchTarget / MorganTarget / RingTarget / ...
    ↓ hot loop runs here
results
```

The views are dense, packed, and purpose-built. Building one takes O(n) for n atoms and is amortized across many queries against the same target. This is exactly how RDKit organizes `ROMol` internally (the public API hides it behind atom/bond iterators).

### Expected outcomes

- **Morgan fingerprints**: specialized view (packed atom hashes + neighbor lists). Within 1.5× of RDKit; potentially faster with careful Rust.
- **Substructure search**: specialized `MatchTarget` with packed adjacency and per-atom element/degree. VF2 runs on the packed view. Should match RDKit.
- **SSSR**: petgraph cycle basis built once per target. Comparable to RDKit.
- **SMILES parse**: current parser runs at < 1 μs / molecule on the existing conformance corpus. Already good; string-processing-bound, Rust nom holds up well. Not a concern.
- **SMILES emit**: not yet measured; not expected to be a bottleneck.
- **Resolution (valence, aromaticity)**: one-time per molecule. 2× slowdown here is invisible because it happens once at parse.

**Where we expect to be faster**: anything involving the constraint engine (multi-relation matching, rule application, derived predicates, mixed numeric + structural queries). RDKit does not have this layer; we do not pay to emulate it.

**Where we expect to be slower**: tight numeric loops over atom attributes that do not benefit from views. Estimate: 1.2–1.5×.

**Where we accept the loss**: anything where RDKit's bespoke optimization matters and we do not have time to match it. SSSR at extreme scale, perhaps.

### The real risk

Not that the approach is slow in theory. That we do not build the views, use `MoleculeAst` directly in hot loops, and absorb the 1.5× tax everywhere. **Discipline**: every algorithm called more than O(n) times per molecule uses a specialized view, not `MoleculeAst` directly. Reviewable as a code-review rule.

### Fail-fast checkpoint

Before committing to the full migration, benchmark Morgan fingerprints on a 10k-molecule set. The existing SMILES parsing conformance suite already contains ~10k molecules and is the natural starting corpus; larger sets can be pulled from ChEMBL if needed.
- Implementation A: directly over `MoleculeAst` (`Option` access, no view).
- Implementation B: with a `MorganTarget` view built from `MoleculeAst`.
- Baseline: RDKit via pyo3 or equivalent.

**If B is within 2× of RDKit**: proceed. **If not**: introduce a packed immutable `MoleculeView` type used on all hot paths, essentially an internal-only RDKit-style representation derived from `MoleculeAst`. This is the escape hatch; designing it is part of the fail-fast plan.

## Caching derived views on immutable molecules

Molecules, once grounded, are immutable. Derived views (ring set, distance matrix, biconnected components, packed `MatchTarget`, `MorganTarget`, per-relation hash indexes) are pure functions of the ground AST and can be cached inside `GroundMolecule` with interior mutability.

```rust
pub struct GroundMolecule {
    inner: Arc<GroundMoleculeInner>,
}

struct GroundMoleculeInner {
    ast: MoleculeAst,
    ring_set: OnceLock<RingSet>,
    distance_matrix: OnceLock<DistanceMatrix>,
    biconnected: OnceLock<BiconnectedComponents>,
    match_target: OnceLock<MatchTarget>,
    morgan_target: OnceLock<MorganTarget>,
}

impl GroundMolecule {
    pub fn rings(&self) -> &RingSet {
        self.inner.ring_set
            .get_or_init(|| RingSet::compute(&self.inner.ast))
    }
    // ... one accessor per cached view
}
```

Key properties:
- **Init-once, thread-safe.** `std::sync::OnceLock` handles concurrent access without a mutex on the hot path after the first init.
- **Cloning is cheap.** `Arc` bump; all cached views are shared across clones. A pattern-matching pass that hands out a molecule to per-thread workers still sees the same cache.
- **Per-instance, not global.** No external cache keyed by pointer identity or canonical hash. Lifetimes follow the `GroundMolecule` instance.
- **Pay-on-use.** Unused views are never computed. Morgan fingerprinting does not pay for ring enumeration; ring-based queries do not pay for Morgan hashes.

APIs that benefit from caching take `&GroundMolecule`. APIs that only need the structure (EDN serialization, pretty-printing) take `&MoleculeAst` directly and bypass the cache layer. `MoleculeAst` itself holds no cache — interior mutability on a homoiconic AST does not compose cleanly with equality, hashing, or constraint narrowing, and ungrounded `MoleculeAst` values (patterns, mid-resolve partial structures) do not benefit from the same views anyway.

For rewrites (kekulization, aromatization, reactions): the output is a new `MoleculeAst`, which is wrapped in a new `GroundMolecule` if the caller needs caching again. The cache does not transfer across rewrites. Invalidation on mutation is not needed, because there is no mutation — rewrites produce new values.

## Ground-term ergonomics

Concern: constant checking for ground-ness if everything is `MoleculeAst`.

Solution: newtype wrapper with once-checked invariant.

```rust
pub struct GroundMolecule(MoleculeAst);

impl GroundMolecule {
    pub fn new(ast: MoleculeAst) -> Result<Self, GroundError> {
        if ast.is_ground() { Ok(Self(ast)) } else { Err(GroundError::NotGround) }
    }
    pub fn as_ast(&self) -> &MoleculeAst { &self.0 }
}
```

APIs that require ground input take `&GroundMolecule`. The check happens once at construction. Inside hot loops, you have a `&GroundMolecule` and call unchecked accessors on the wrapped AST.

`is_ground()` is O(n): all `AtomAst` and `BondAst` fields concrete, constraints vec empty or contains only satisfied assertions.

Homoiconicity is preserved at the AST level; type safety is preserved at the API level. The check cost is paid once per molecule, not once per access.

## Migration plan

Ordered steps. Each step leaves the tree green.

1. **Add `constraints: Vec<MoleculeConstraint>` and the enum with the D3 starter set** to `MoleculeAst`. Constraints vec always empty for now. Parser, serializer, resolution unchanged. ~2 days. *(Done.)*
2. **Fold `charge` and `spin` fields into constraint emissions** at parse time and remove the struct fields. ~3 days. *(Done 2026-04-11.)*
   - `:charge` / `:spin` keys parse into `MoleculeConstraint::Derived { TotalCharge | TotalSpin, RelationRefs::default() }` and emit back the same way through `MoleculeAst`.
   - The `MoleculeBuilder` carrier for explicit charge/spin (`set_charge` / `set_spin` plus the `MolecularChargeMismatch` / `MolecularSpinIncompatible` validation in `build()`) was removed in the same step. It was never load-bearing — no fixture in the conformance corpus or unit tests asserted a molecule-level `:charge` / `:spin`. Without an asserted value, `build()` always derives charge from atom/bond/aromatic/multicenter sums and picks the unique compatible total spin (or errors with `MolecularSpinIncomplete` if multiple are compatible). The validator handler the original step called for is therefore moot at this stage; reintroduce it when a downstream pass actually wants to assert something about the resolved molecule.
   - Ancillary in this step: introduce `SpinStateAst` (Wildcard / Lit / Pair) so atom and bond `spin` fields carry the same coupled (unpaired, multiplicity) shape that `TotalSpin` does; widen `Derived` from `atoms: Vec<usize>` to `refs: RelationRefs`; add `And` / `Or` / `Not` combinators (no consumer yet — they exist so future predicates can compose without another schema change).
   - Shared-crate cleanup pulled in alongside the ancillary work (2026-04-11): the top-level `umol` crate was merged into `umol-shared` (the only member was the `UmolError` trait with its `as_any()` downcast hook, which now lives in `umol_shared::error`); the `units` module was factored into `units::{length, angle, time}` with a unified physical-quantity API (each type uses a `T::new(value, Unit::X)` unit-enum constructor with named-constructor sugar like `Length::bohr` / `Angle::degrees` / `Time::seconds` — all `const fn` so the isotope static tables can be populated directly); `HalfLife` was collapsed into `Time` and moved from its own file into `isotope.rs`; `umol-shared` switched to module-qualified re-exports (Option B — callers write `umol_shared::element::Element`, not flat `umol_shared::Element`) and the `e!` / `iso!` / `occ!` / `spin!` macros became hygienic via `$crate::module::Type`. No API-level behavior change; purely a reorganization to stop the shared crate from accumulating a flat top-level surface as it grows.
3. **Implement `GroundMolecule` newtype** and migrate ground-requiring APIs. ~3 days. *(Done 2026-04-12.)*
   - Plain newtype `GroundMolecule(MoleculeAst)` with `new` / `as_ast` / `into_ast`. No `Arc`, no `OnceLock` view caches — deferred to step 5 (benchmark checkpoint decides whether cached views are needed, and if so whether they live on the molecule or in a separate `MoleculeViews` wrapper).
   - `is_ground()` predicates added bottom-up: `ValueAst`, `ElementAst`, `IsotopeAst`, `HydrogenAst`, `AromaticValenceAst`, `AtomAst`, `BondAst`, `DerivedPred`, `MoleculeConstraint::is_ground_assertion()`, `MoleculeAst`.
   - `GroundError` unit struct in `ast::error`.
   - `Molecule::to_ground()` added on the graph_ir side (infallible — `Molecule` is ground by construction).
   - `AromaticAst` renamed to `AromaticValenceAst` to disambiguate from `AromaticSystem`.
4. **Relational matcher over core relations** (no constraint vec support yet). Enough for basic substructure search. ~2 weeks. *(Done 2026-04-12.)*
   - `matches` methods added to `ElementAst`, `IsotopeAst`, `HydrogenAst`, `AromaticValenceAst` (pattern-against-ground matching) in `umol-shared/src/atom_ast.rs`.
   - `AtomAst::matches_ground` and `BondAst::matches_ground` for field-by-field matching with None-semantics (None on query = unconstrained; Some on query vs None on target = no match).
   - New module `ast::matcher` with `MatchTarget`, `MatchQuery`, `Assignment`, `find_matches`.
   - VF2 via petgraph `subgraph_isomorphisms_iter` over localized bonds only (directed graph with bidirectional edges; graph construction encapsulated so composite topology graph is an internal refactor).
   - Post-filters: dative bonds (directed), noncovalent bonds (directed), aromatic systems (mapped atoms subset of single target system), multicenter bonds (mapped atoms subset of single target bond — handles overlapping membership, e.g. B2H6).
   - Not implemented: constraints vec, variable bindings, `MatcherFlag`, `DerivedPred` evaluation, `SubPattern` recursion (step 9).
5. **Benchmark checkpoint**: Morgan fingerprints, `MoleculeAst` vs. view vs. RDKit. If within 2×, proceed. If not, design the packed `MoleculeView` escape hatch before continuing. ~3 days. *(Done 2026-04-12.)*
   - ECFP algorithm (Rogers & Hahn 2010): 7 initial Daylight invariants (heavy degree, heavy valence, atomic number, atomic mass, charge, H count, ring flag), iterative hashing with sorted (bond_order, neighbor_id) pairs, duplicate structure removal via bond-set tracking, dead atom pruning.
   - New module `ast::morgan` with `MorganFingerprint`, `MorganTarget`, `morgan_direct`, `morgan_view`.
   - Criterion benchmark over 9,120 SMILES from conformance corpus (ECFP4, radius 2):

     | Implementation | Total | Per molecule |
     |---|---|---|
     | `morgan_direct` (MoleculeAst) | 126 ms | 13.8 μs |
     | `morgan_view` (MorganTarget, pre-built) | 97 ms | 10.6 μs |
     | `morgan_view` (with build) | 128 ms | 14.1 μs |
     | RDKit (C++, via Python) | 132 ms | 14.7 μs |

   - **Result: proceed.** Direct MoleculeAst access is at parity with RDKit (1.05×); pre-built view is 1.4× faster. No `MoleculeView` escape hatch needed. Unoptimized Rust — constant-factor improvements remain (allocation reuse, inline bitsets, pre-sorted adjacency).
   - SMILES parsing comparison on the same corpus (not directly related to the AST migration but recorded here for reference):

     | Implementation | Total | Per molecule |
     |---|---|---|
     | umol `parse_smiles` (nom) | 6.4 ms | 0.7 μs |
     | RDKit `MolFromSmiles` (no sanitize) | 85 ms | 9.3 μs |
     | RDKit `MolFromSmiles` (sanitized) | 446 ms | 49.7 μs |

   - umol SMILES parsing is 13× faster than RDKit's parse-only path.
6. **Migrate valence resolution and atom typing** to emit constraints; solver runs finite-domain propagation. Delete the old bespoke resolver. ~2 weeks. *(Done 2026-04-14.)*
   - Prerequisite refactor (doc 83 § "Absent vs undetermined"): rename `Wildcard` → `Undetermined` across all AST types (done 2026-04-13); remove `Option<>` wrappers from `AtomAst` and `BondAst` fields — all fields now carry their `Undetermined` default directly (done 2026-04-14); `is_ground()` returns `false` for undetermined fields (done 2026-04-14); `needs_narrowing` deleted, replaced by `!is_ground()`.
   - `solver.rs` module exists with `Solution<T>`, `Progress`, `ValenceStrategy` (AtomTyping / Counts), `Solver` with `resolve`, `validate`, `filter`. 15 tests passing.
   - `narrow_atom` now also resolves `isotope_mass` and returns `bool` to prevent infinite loops.
   - Topology helpers (`bond_order_sum`, `dative_bond_order_sums`, `is_in_aromatic_system`) moved to methods on `MoleculeAst`; `charge_or_zero` on `AtomAst`.
   - `find_matches_with` in `ast::matcher` composes `find_matches` + `solver.filter` as a post-filter.
   - Index types: `AtomIdx`, `BondIdx`, etc. as newtypes over `usize` with `index_vec::IndexVec` for ergonomic indexing. `atoms![]` macro for test construction.
   - Old resolver deleted (done 2026-04-14): `graph_ir::valence` module removed; `ValenceStrategy` enum consolidated in `solver.rs` and re-exported from `graph_ir::config`; `resolve_valence_with` now calls private `valence_candidates` which dispatches on `ValenceStrategy` variants using `AtomPattern`-based candidate generation directly (no AST bridge). 621 conformance tests pass, 3796 lib tests pass.
7. **Migrate aromaticity perception** to discover-then-verify. Graph analysis stays; perception output is constraint emissions. ~1 week. *(Discovery half done 2026-04-15.)*
   - `Graph::induced_subgraph` added to umol-graph-core for filtered ring enumeration over atom subsets.
   - `RingEnumerator::enumerate_ast` added: filters atoms by aromatic valence (optional), builds induced subgraph, runs BCC + cycle enumeration per component, maps back to original indices.
   - `find_from_ast` on all three perception models (HueckelRule, HMO, Clar) and `AromaticityModel::aromatic_systems_ast` dispatch.
   - `MoleculeAst::set_aromatic_systems` replaces the `VarRelationSet` wholesale with discovered systems.
   - `AromaticityConfig` (strategy + ring enumeration) added to `Solver`. Resolve loop: valence→aromaticity→re-valence stratification. `Solver::resolve` returns `Result<Solution<()>, AromaticityError>`.
   - Temporary index bridge: `Ring`/`RingSet`/`AromaticSystem` still use petgraph `AtomIndex`; AST methods convert via `AtomIndex::new(node_id.index())`. Goes away when GraphIR is removed (step 8).
   - Not done: constraint emission (`Derived { AromaticElectronCount(n), ring_atoms }`) and verification mode (solver discarding failing proposals). These require `DerivedPred::AromaticElectronCount` and constraint evaluation infrastructure (step 9 scope). Currently all proposed systems are accepted unconditionally.
   - 621 conformance tests pass, 3796 lib tests pass.
8. **Delete `graph_ir::Molecule`.** `MoleculeAst` is the only molecule type. ~3 days.
9. **Add `SubPattern` constraint and matcher recursion.** ~1 week.
10. **Add DPO rule application**. `ReactionRuleAst` as top-level type, rule-apply as a transform over `MoleculeAst`. ~1–2 weeks.

Total: ~8–10 weeks. Fits inside the three-month window with buffer.

## What this doc does NOT commit to

- **Generic CSP frontend.** D9 explicitly rejects this.
- **MSO₁ or MSO₂ evaluators.** Deferred. Derived structural predicates cover the immediate needs.
- **Path / Kleene queries.** Deferred. `Path` not a stub variant.
- **Probabilistic matching** (belief propagation, `Family<Product>`). Deferred.
- **Open-world schema registry.** Closed-world Rust types only.
- **Schema-parametric AST** (doc 79 §"Relation schemas as organizing principle"). Not required for the three-month window; revisit if reaction networks or crystals demand a second schema.
- **Backwards compatibility with the current `Molecule` API.** No concessions; rewrite callers.
- **Incremental constraint solving.** Full re-solve on change; revisit if profiling demands it.

## Follow-ups

- **Doc 81 (or similar)**: detailed sketch of how each current resolution pass in `graph_ir/resolution.rs` maps to constraint emissions. Write when the matcher lands and resolution migration begins — not before, since the matcher implementation will inform the constraint shape.
- **Benchmark harness** for the fail-fast checkpoint (step 5). Needs a 10k-molecule set (ChEMBL sample or equivalent) and a pyo3 bridge to RDKit for the baseline.
- **Update memory**: this doc supersedes several earlier framings; the key one is that `MoleculeAst` subsumes `Molecule`, `MoleculeBuilder`, and `MoleculePattern` at the storage level. The feedback memory about "transient resolution flags belong on `MoleculeBuilder`" needs to be revised in light of `MoleculeBuilder` disappearing as a type.
