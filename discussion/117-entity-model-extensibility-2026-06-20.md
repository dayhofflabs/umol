# Entity model extensibility — dative/haptic split and new entity types

Status: **Active / analysis.** No code authorized. Scoping doc, not an implementation plan.
Date: 2026-06-20.
Trigger: the `DativeBond` entity conflates two-atom dative bonds and haptic bonds, which are
topologically heterogeneous. This raises the broader question of which binding situations deserve
their own entity and what adding one costs.

## Governing principle

An entity (or relation) earns its place when its participants are a **localized, enumerable subset of
atoms** and graph-or-relation notation is the right abstraction. It does not when the structure is:

- **QM-delocalized** — every atom contributes to every relation (molecular orbitals over the whole
  framework, band structure). Encodable as an all-atoms relation, but the abstraction carries no
  locality and buys nothing.
- **Periodic / infinite** — infinite lattices, periodic coordination polymers, extended solids (no
  finite participant set). Finite oligomers and non-periodic polymers are a separate, unsettled
  question — see §4.
- **Ensemble / superposition** — resonance structures, tautomer mixtures, crystallographic disorder.
  These are a *collection of graphs*, a different axis from a single molecular graph.

The goal is to make more chemistry representable where locality holds — not to maximize entity count.

## 1. Dative/haptic split — why / why not

Current `DativeBond` stores both a two-atom donor→acceptor bond (B←N in borazine) and a haptic bond
(η⁵ Cp→Fe) in one type:
`FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, DativeBondAst>` — one fixed acceptor +
a variable donor set.

### Keep combined
- **Electron budget**: `#t` / `#d` (donated / accepted pairs) describes both uniformly; the
  pair-transfer model is the shared chemistry that motivated the merge.

### Split
- **Topology**: a two-atom dative bond is a single (directed) edge — it can sit in a ring exactly like
  a localized bond. A haptic bond is a hyperedge (one metal + n contiguous ligand atoms), where ring
  membership is ill-defined. Concrete symptom today: hyperedge entities expose no `ring_*` derived
  predicates, so the two-atom subcase is dragged to the lowest topological common denominator and
  cannot expose ring predicates even though it is a genuine edge.
- **Geometry**: haptic ligands are only approximately, or not at all, rotationally symmetric; the
  geometric/stereo treatment differs from a point-to-point dative bond.
- **Redox-active ligands**: non-innocent ligands donate *unpaired electrons*, not whole pairs, which
  the pair-based dative model cannot express.

### Storage mapping if split
Both targets are already-existing relation shapes:
- two-atom dative → a fixed 2-arity directed relation (`FixedRelationSet`-class, ordered/directed), an
  ordinary edge with a donor→acceptor orientation.
- haptic → the current `FixedVarBirelationSet` (1 metal + variable hapto set), a hyperedge alongside
  `multicenter_bonds`.

The split is clean at the storage layer; the cost is the per-entity surface (§3), not new plumbing.

## 2. Other entity types — scaffold

Shallow enumeration, no deep dive. "Existing" = already an entity; "candidate" = fits the locality
principle and lacks a home; "out" = excluded by §0.

| situation | natural relation shape | status |
|---|---|---|
| localized bond | graph edge | existing (`bonds`) |
| two-atom dative | directed edge | existing within `dative_bonds`; **split target** |
| haptic / π-coordination (ηⁿ) | hyperedge: 1 metal + n ligand atoms | existing within `dative_bonds`; **split target** |
| multicenter (3c-2e: B–H–B, carboranes) | variable hyperedge | existing (`multicenter_bonds`) |
| aromatic / delocalized π ring system | variable hyperedge | existing (`aromatic_systems`) |
| noncovalent (H-bond, halogen bond, …) | fixed 2-arity | existing (`noncovalent_bonds`) |
| agostic interaction (M···H–C) | 3-atom relation, or overlay on a C–H bond | candidate |
| charge-transfer / EDA / π-complex | relation over two fragments | candidate (arity/definition fuzzy) |
| metal–metal bonds (incl. multiple, δ) | edge with extended order/character | candidate — may extend `bonds` rather than be new |
| bridging-ligand classification (μ₂, μ₃) | derived from ordinary bonds | candidate — likely derived predicate, not an entity |
| redox-active / non-innocent ligand state | per-ligand electron/spin annotation | candidate — likely attribute/constraint, not a topological entity |
| delocalized MOs / band structure | all-atoms relation | out (no locality) |
| coordination polymer / extended solid | periodic | out (infinite) |
| resonance / tautomer / disorder ensembles | set of graphs | out (different axis) |

Note the recurring fork: several candidates are better served as **derived predicates** or **constraints**
on existing entities (μ-classification, M–M order, redox state) than as new entity types. Each candidate
needs the entity-vs-attribute question answered before it warrants §3 effort.

## 3. Effort of adding a new entity type

What a new entity touches today ("the whole nine yards"):

1. **AST type** (`ast/<entity>.rs`): struct/enum + fields, `Canonicalize` + `Lattice` derives,
   `Hash`/`Ord`, constructors, `is_ground`/`matches`.
2. **Storage** (`molecule.rs`): a relation-set field (choosing among the existing primitives —
   `FixedRelationSet<N>`, `VarRelationSet`, `FixedVarBirelationSet`), constructor params + wiring,
   `Clone`, accessors; plus an id type (`id.rs`) and incidence wiring (`incidence.rs`).
3. **Views** (`ast/view/<entity>.rs`): `<Entity>Views` namespace + `<Entity>View` + derived predicates;
   registration on the molecule/graph view.
4. **Constraints** (`ast/constraint/<entity>.rs`): constraint enum, `key()`, `Canonicalize`, `Ord`,
   container; any molecule-level relational constraints referencing it.
5. **Inter-entity relations**: derived relational predicates and relational constraints linking the
   entity to *other* entity types — `incident`, `is_in_aromatic_system`, coincident/incident between
   bonds and stereo bonds, etc. A new entity must define its relation to every existing entity it can
   interact with.
6. **Surface syntax** (`dsl/<entity>.rs`): `FromStr`/`Display` and EDN `FromEdn`/`ToEdn` on `<Entity>Dsl`,
   predicate tags; defaults in `dsl/config.rs`; molecule-DSL wiring.
7. **Algorithm touchpoints**: anything iterating all entities or computing canonical forms — coloring,
   matching/embedding, remap, reaction, symmetry.
8. **Tests + conformance corpus**.

### What is mechanical vs essential
- **Mechanical** (reusable / boilerplate): the storage primitives already exist; relation-set plumbing,
  id/incidence wiring, view-namespace registration, and much DSL/EDN boilerplate are repetitive and
  could be factored (macro or trait-driven entity scaffold) without changing semantics.
- **Essential** (cannot be generated): the derived predicates, canonicalization/meet rules, surface
  syntax design, and inter-entity relations are where the actual modeling lives — these are per-entity
  by nature.

The inter-entity relations (5) are the dominant long-tail cost and the least scaffoldable part: the
number of entity-pair relations grows quadratically in entity count, the relations are heterogeneous (a
single pair can carry several distinct relations — e.g. a bond and a stereo bond can be *coincident* or
merely *incident*), and some pairs have no clean definition at all (aromatic system ↔ multicenter bond).
Each new entity adds a row and a column of bespoke relational modeling, most of which cannot be generated.

So the realistic lever is reducing the mechanical surface (so the split and future candidates cost
mainly the essential design work), not eliminating entity types or auto-generating their semantics.
Whether to build that scaffold before the dative/haptic split, or do the split by hand first and extract
the pattern from it, is the open decision.

## Open decisions

- Confirm the dative/haptic split and the two target storage shapes (§1).
- For each §2 candidate: entity vs derived-predicate vs constraint.
- Whether to invest in an entity-scaffolding mechanism (§3) before or after the first split.
