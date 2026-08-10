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

### Current validation boundary

The `#d`/`#t` projection is defined cleanly only for the binary donor→acceptor case. Validation may
cross-check binary dative bond order against the donor and acceptor projections. The corresponding
multi-donor validation path remains a stub/TODO: it must not establish zero, `Underdetermined`, or an
apportionment rule as lasting semantics while the current entity still conflates coordination/haptic
bonding with binary dative bonding.

Do not implement the split merely to complete that validator. First reduce or organize the mechanical
cost of extending the entity-kind set—constraints, DSL, views, deltas, edits, undos, updates, remapping,
matching, and bindings—so that the chemical distinction does not require another bespoke sweep of the
entire stack. The validation stub is therefore an explicit dependency on the broader extensibility
decision, not evidence that multi-donor donation equals zero.

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

## 4. Case study — adding a constraint variant (noncovalent `#I` Intramolecular)

Empirical companion to §3, one axis over: §3 measured the cost of a new *entity*; this measures the
cost of a new *constraint variant* on an existing entity. Concrete trigger (2026-07-12): the noncovalent
bond's constraint enum was left **uninhabited** (`enum NoncovalentBondConstraintAst {}`) — a deferral
("bigger fish to fry"), not a structural claim that noncovalent bonds admit no constraints. It surfaced
while binding the Python noncovalent slice (doc 140 B5): pyo3 rejects a zero-variant `#[pyclass] enum`
(`"#[pyclass] can't be used on enums without any variants"`), so the empty enum cannot be mirrored and
the container cannot carry a mapping API. Rather than build a degenerate Python stub around the emptiness,
the right move is to inhabit the constraint that should have been there — which also **unblocks B5** (the
surface becomes an ordinary 1-key Boolean constraint) and documents the feature-addition cost.

### The variant

- `NoncovalentBondConstraintAst::Intramolecular(BooleanAst)`; key `NoncovalentBondConstraintKey::Intramolecular`.
- String DSL tag `#I`: `#I` / `#I+` → true (intramolecular), `#I!` → false (intermolecular),
  `#I*` → undetermined. EDN: single-key map `{:intramolecular <bool>}`.
- Structural twin: dative/bond `Aromatic(BooleanAst)` (`#a`) — same Boolean value, same `+`/`!`/`*`
  combinator, same single-`#`-tag shape — reduced to a **1-variant** container, whose shape is the
  aromatic-system `ElectronCount` container (single-value accessor `intramolecular()`, `matches` via the
  accessor). So: *bond-`Aromatic` value & DSL* + *aromatic-`ElectronCount` container*.

### Scope map (verified against the AST, 8-reader scoping pass — to reconfirm at implementation)

| Subsystem | Files | Change |
|---|---|---|
| Constraint element + Key | `ast/constraint/noncovalent.rs` | Add the `Intramolecular(BooleanAst)` variant + `Intramolecular` key; fill the ~15 empty `match self {}`/`match key {}` bodies (`key`/`compact`/`remap`/`Canonicalize`/`Lattice`×6) with the `Aromatic`-arm mirror; add ctor `intramolecular(b)` + `as_undetermined` (noncovalent omits it today). |
| Constraint container | `ast/constraint/noncovalent.rs` | Real bodies for `find`/`set`/`compare_and_set`/`extend`/`update` (mirror aromatic); container `Canonicalize`/`Lattice` (currently hardcoded trivial); `FromIterator` (currently **ignores** its iterator); missing `IntoIterator`/`From` impls; new `intramolecular() -> BooleanAst` accessor. `get`/`contains`/`remove`/`retain`/`clear`/`take`/`iter`/`len` already delegate — go live unchanged. |
| Value type | `ast/noncovalent.rs` | `with_constraints` is currently a **no-op** (`_constraints`); there is **no** `with_constraint` singular — fix/add both; check `update`/`into_ground`. `From<&str>`/`RelationData`/`new`/`from_kind`/`with_kind` unaffected. |
| Molecule-level `Constraint` | `ast/constraint/molecule.rs` | The top-level `Constraint::NoncovalentBond(id, payload)` variant **already exists** (line 45) — no add. Fill its `canonicalize` arm; `is_vacuous` doc; `inline_constraints` arm. |
| Delta | `ast/delta.rs` | `EntityPatch for NoncovalentBondDelta`: real `apply_constraint` (no-op → `compare_and_set`) + `diff_constraints`. `EntityFold::constraint_key` already `constraint.key()` — no edit. |
| Edits + transact + undo | `ast/edit.rs`, `ast/transact.rs` | **New** `Edit::ModifyNoncovalentBondConstraint` variant; `apply_edit` + `apply_edit_checked` (undo) dispatch arms; new `apply_modify_noncovalent_bond_constraint` helper. |
| Reactions | `ast/reaction.rs` | Deltas→Edits lowering arm (`Delta::NoncovalentBond(ModifyConstraint)`); delta→EDN partial-render fold. `apply` has-conflict check already generic. |
| Reaction span | `ast/reaction_span.rs` | **Zero changes** — `materialize`/classification folds `ModifyField \| ModifyConstraint` and replays through the generic `apply_noncovalent_change`; goes live automatically once delta stops returning empty. |
| String DSL | `dsl/noncovalent.rs` | **New** `NoncovalentBondPredicate` enum + `noncovalent_bond_predicate` parser + `apply_predicates` + `constraint_tag`/`fmt_constraint` (mirror bond `#a`); wire into `noncovalent_bond`, `fmt_noncovalent_bond_form`, `PartialNoncovalentBondDsl::fmt`. |
| EDN serde | `dsl/noncovalent.rs`, `dsl/constraint.rs` | Inhabit the DSL boundary enum `NoncovalentBondConstraintDsl` (also `{}` today) with `Intramolecular(BooleanDsl)` + `from_ir`/`into_ir` (don't exist) + `FromEdn`/`ToEdn`; fill `read_noncovalent_bond_constraint_dsl` + `ConstraintDsl::{to_edn,from_ir,into_ir}` noncovalent arms. |
| Errors | `dsl/error.rs` | **New** `ParseError::{UnknownNoncovalentBondPredicate, DuplicateNoncovalentBondPredicate}` (mirror the dative pair). |
| Property tests | `tests/property/strategies.rs`, `tests/property/lattice.rs` | New `noncovalent_bond_constraint(s)_strategy`; wire into `noncovalent_bond_form_strategy`; new lattice-law proptests for the constraint + value (a today-vacuous coverage slot that becomes meaningful). |
| Config | `dsl/config.rs` | `NoncovalentBondDefaults` — likely no change (verify). |

**Genuinely new** (everything else fills an empty `match {}` or copies a peer): the `Intramolecular`
AST variant + key, the `NoncovalentBondConstraintDsl::Intramolecular` variant, the two `ParseError`
variants, the one `Edit` variant, and the string-DSL predicate machinery. No new *design* — every piece
has an exact `Aromatic(#a)` / `ElectronCount` precedent.

### Staged impl plan

- **A — inhabit the AST layer (atomic red→green; ~11 files).** Adding the variant breaks *every*
  exhaustive `match` on the three currently-uninhabited enums (`NoncovalentBondConstraintAst`, `…Key`,
  `NoncovalentBondConstraintDsl`) at once, so they must all be filled in one edit to restore compilation.
  Sequence within: constraint element+container → value type → molecule `Constraint` arms → delta
  apply/diff → new `Edit` variant + transact dispatch/undo/helper → reaction lowering + delta→EDN render
  → DSL boundary type + `from_ir`/`into_ir`/`FromEdn`/`ToEdn` + `ConstraintDsl` arms + reader +
  string-DSL **render** (`fmt`). `reaction_span.rs` unchanged. Green when the workspace compiles and the
  container + lattice unit tests pass. Each edit mirrors a peer — no essential design.
- **B — string-DSL parse `#I` (additive/green).** The predicate parser/enum/applier + the two
  `ParseError` variants + parse/roundtrip tests. Purely additive (doesn't break compilation; makes
  `"#I"`/`"#I+"`/`"#I!"`/`"#I*"` parse). `[dep: A]`
- **C — property tests (additive/green).** Noncovalent constraint strategies + lattice-law proptests.
  `[dep: A]`
- **D — resume B5 (Python noncovalent slice).** Now a standard slice: the *constraint* half is the
  bond-`Aromatic` Boolean shape in a 1-key container (`intramolecular` getter/setter), the *view* half is
  bond-shaped (2-atom `(first, second)` pair, `connecting(a, b)`, `atom_ids` 2-tuple) plus the new
  `NoncovalentBondKindAst`/`NoncovalentBondKind` leaf. The zero-variant-pyclass blocker is gone.
  `[dep: A (+B for a Python parse test)]`

Critical path: **A → (B ∥ C) → D**.

### Observation — the cost of a constraint variant vs an entity (§3 companion)

- **Wide but shallow.** ~11 files, but every touch either fills an empty `match {}` with a body copied
  from a peer or adds a variant mirroring `Aromatic(#a)`. New *design* surface ≈ zero. Contrast §3's
  entity-addition, which is O(entity²) bespoke relational modeling. **Lesson: put new chemistry on an
  existing entity's constraint axis when locality allows — it is dramatically cheaper than a new entity.**
- **Atomicity is the only real friction.** The scaffold chose *uninhabited enum + empty `match self {}`*
  as the placeholder. That compiles cleanly today, but the day a variant lands it breaks all exhaustive
  matches simultaneously → one large red→green edit that cannot be incrementally green. The alternative
  placeholders (a first real variant from the start, or an `Option`-free always-inhabited enum) trade a
  slightly-less-clean today for a much smaller edit later.
- **The generic layers held.** `reaction_span.rs` needing zero changes, and `Constraint::NoncovalentBond`
  / `EntityFold::constraint_key` / most container delegators already being real, show the materialization
  and molecule-scope layers are correctly generic over the constraint payload. The friction is
  concentrated in the per-entity leaf + its DSL — exactly where a scaffolding macro (§3) would pay off:
  it would turn this ~11-file spread into a one-line variant addition.
- **Extensibility verdict.** Adding a value-only constraint variant to an existing entity is a bounded,
  mechanical, mirror-a-peer task whose cost is proportional to the number of `match` sites, not to
  conceptual difficulty. This is the encouraging half of the extensibility story (the discouraging half
  is §3's quadratic entity-pair relations).

### Progress — Stage A done (2026-07-12)

Stage A (inhabit the AST layer) is complete: workspace compiles, 4542 umol-ast tests pass. The edit
touched 11 files, all mirroring a peer. Refinements against the scope map found at implementation:

- **Delta became *less* special-cased, not more.** The scope map planned "real `apply_constraint` +
  `diff_constraints`" on the hand-written `EntityPatch for NoncovalentBondDelta`. But the hand-written
  block existed *only* because the constraint was uninhabited (a comment said the `diff_field_ops!`
  macro's constraint loop "would be unreachable"). Since the macro's field handling (`Kind => kind`) is
  byte-identical to the hand-written `apply_field`/`diff_field`, the whole block collapsed to the shared
  `diff_field_ops!(…, { Kind => kind })` macro + a one-line `apply_constraint` — exactly the peer shape.
  So the variant *removed* a special case rather than adding one: the cleaner "fill the absence
  consistently" outcome, not "handle the absence specially."
- **Two match sites were outside the scope map.** `ast/molecule.rs` (top-level, not
  `ast/constraint/molecule.rs`) `inline_constraints` had an empty `match inner {}` to fill; its `lift`
  counterpart (`lift_constraints`) was already generic. The map's "molecule-level `Constraint`" row named
  only `ast/constraint/molecule.rs`. Lesson: the inline↔molecule-scope constraint plumbing is split
  across *two* `molecule.rs` files.
- **Undo needed no new variant.** `transact.rs` (at `ast/molecule/transact.rs`, not `ast/transact.rs`)
  encodes constraint-edit undo as `Undo::ApplyEdit(Box::new(<inverse Edit>))`, so the only new symbol was
  the forward `Edit::ModifyNoncovalentBondConstraint` + its two dispatch arms + one apply helper.
- **`is_vacuous` doc was stale.** The `Constraint::NoncovalentBond` arm already delegated to
  `is_undetermined`; only a doc comment claiming it "always non-vacuous" needed correcting.

The remaining `reaction_span.rs`-zero-changes and `EntityFold::constraint_key`-unchanged predictions
held. Net confirmation of the verdict: the friction was the atomic red→green fan-out (~15 `match` sites),
not any design.

- **B done (2026-07-12).** String-DSL `#I` parse: `NoncovalentBondPredicate` enum +
  `noncovalent_bond_predicate` parser + `apply_predicates` + `constraint_tag`, wired into
  `noncovalent_bond`; the two `ParseError` variants; `PartialNoncovalentBondDsl::fmt` now renders
  undetermined explicitly (`#I*`) so a constraint→undetermined change round-trips. Purely additive —
  4551 umol-ast tests pass, workspace green. Mirrored the dative `#a` predicate machinery exactly; the
  only asymmetry is noncovalent's single predicate vs dative's two.
- **C done (2026-07-12).** `noncovalent_bond_constraint(s)_strategy` + wired into
  `noncovalent_bond_form_strategy`; `test_noncovalent_bond_form_lattice_laws` +
  `test_noncovalent_bond_constraints_lattice_laws`. The high-value side effect: the *existing*
  molecule/reaction/edit/substructure round-trip proptests now generate noncovalent `#I` constraints and
  drive them through the whole Stage A pipeline (delta ↔ edit, transact + undo, reaction lowering, EDN
  render, string parse) — 111 property tests pass, so the Stage A threading is round-trip- and
  lattice-law-clean under fuzzing, not just at the hand-written cases. Also filled a pre-existing gap:
  there was no entity-level `NoncovalentBondAst` lattice test at all (the bond was kind-only until now).
- **D pending.** Resume B5 (Python noncovalent slice). `[dep: A (+B for a Python parse test)]`.

## Open decisions

- Confirm the dative/haptic split and the two target storage shapes (§1).
- For each §2 candidate: entity vs derived-predicate vs constraint.
- Whether to invest in an entity-scaffolding mechanism (§3) before or after the first split.
- §4: resolved — proceeding. Stage A done (2026-07-12); B (string parse) / C (property tests) next, then
  D (resume B5).
