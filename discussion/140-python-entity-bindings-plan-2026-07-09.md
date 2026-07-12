# 140 · Python bindings for the remaining entity ASTs (plan)

Status: Active — **B1 · Bond slice DONE incl. the view half** (value + WET constraint surface +
`BondView`/`BondViews` + molecule-backed constraint views; Python `from_atoms` replaced by
`from_atoms_and_bonds`; Rust 201 / Python 237 green, clippy clean). **B2 · dative FULLY PLANNED**
(staged, settled 2026-07-11 — see *B2 · Dative — staged impl plan* below): Part A = the Rust
`*ViewMut` prerequisite (doc 137 pt6, in full) unblocking every overlay view half; Part B = the dative
Python slice; constructor becomes atoms-positional `from_parts`. Ready to implement Part A (S0). See
*Bond slice — build state + detailed plan* at the end for the WET template.
Date: 2026-07-09
Relates: 137 (atom slice — the template being mirrored), 139 (mutability/hashing/equality
balance), 114 (interning — where stereo/handle-identity deferrals live)

## Scope and order

Bind the seven remaining entity ASTs to Python, mirroring the atom binding (137), in this
order: **bonds → dative → aromatic → multicenter → noncovalent → stereo-atom →
stereo-bond**. The order is increasing complexity: bond defines the template; dative,
aromatic, multicenter, noncovalent reuse it with narrowing/widening; the two stereo
overlays are a separate, larger sub-project.

Excludes the already-deferred items (handle `__eq__` / interning, molecule DSL `parse`/
`str`, `append`/`extend`, slicing, `MoleculeAst` `Hash`).

## Rust facts this rests on

Storage in `MoleculeAst` (all `Arc`-wrapped, copy-on-write via `Arc::make_mut`); each family
has the uniform accessor trio `xs() -> XViews` / `x(id) -> XView` (panics if absent; `xs().
get(id)` for `Option`) / `x_mut(id)`, and `MoleculeAst: Index<XId>` → `&XAst`:

| Family | storage | id | participants |
|---|---|---|---|
| atoms | `Vec<AtomAst>` | `AtomId(u32)` ↔ NodeId | — (is a node) |
| bonds | `Vec<BondAst>` | `BondId(u32)` ↔ EdgeId | 2 atoms, from the graph edge (undirected); edge-indexed by atom pair |
| dative | `FixedVarBirelationSet<NodeId,Ordered,1, NodeId,Unordered, DativeBondAst>` | `DativeBondId` | 1 acceptor + N donors (N>1 = haptic); donors unordered |
| aromatic | `VarRelationSet<NodeId,Unordered, AromaticSystemAst>` | `AromaticSystemId` | variable unordered atom set |
| multicenter | `VarRelationSet<NodeId,Unordered, MulticenterBondAst>` | `MulticenterBondId` | variable unordered atom set |
| noncovalent | `FixedRelationSet<NodeId,Unordered, NoncovalentBondAst, 2>` | `NoncovalentBondId` | fixed unordered atom pair |
| stereo-atom | `FixedVarBirelationSet<NodeId,Ordered,1, StereoLigand,Ordered, StereoAtomAst>` | `StereoAtomId` | site atom + ordered ligand frame; ≤1 per site |
| stereo-bond | `FixedVarBirelationSet<EdgeId,Ordered,1, StereoLigand,Ordered, StereoBondAst>` | `StereoBondId` | site bond + ordered ligand frame; ≤1 per site |

Critical: **no entity value struct carries its atom ids** — endpoints/participants live in the
graph (bonds) or the relation set (everything else). The `*Ast` value is only chemistry
attributes + a `constraints` field. (Implementation note, not an API difference: today
`x_mut(id)` returns a `*ViewMut` guard for atoms/bonds but a bare `&mut XAst` for the six
relations, so the binding setters dereference one level less for relations. 137 fold-back #6
makes `*ViewMut` uniform — an interning prerequisite per 114 — after which every setter routes
the same way; the bindings should assume the uniform form.)

## A. The shared template (generalized from the atom binding)

The atom binding (137) is six component groups. Each entity reuses them:

1. **Value type `XAst`** — standalone, mutable, wrapping `AstXAst`. Settable field
   properties; `asdict`; `constraints` getter → live `XConstraintsView`; `#[pyclass(eq)]`
   (value-equal, unhashable, per 139); `__repr__`. (`parse`/`__str__` only where a DSL
   exists — deferred with the molecule DSL for now.)
2. **Handle `XView`** — on `MoleculeAst`: id + participants + settable value fields routing
   through `x_mut(id)`; `constraints` → live `XConstraintsView`; `asdict`; `id`; `__repr__`.
3. **Collection `XViews`** — `mol.<family>` → `__len__` / `__getitem__` (negative-index
   normalized, per finding 137-p3-1) / `__iter__`; id-indexed.
4. **Constraints** — value `XConstraintsAst` + live `XConstraintsView`, the mapping surface
   (`set`/`pop`/`update`/`get`/`keys`/`values`/`items`/`__getitem__`/`__delitem__`/
   `__contains__` + per-key accessors + `asdict` + subscript proxies where keyed), mutable +
   unhashable.
5. **Key `XConstraintKey`** (hashable) + constraint value enum `XConstraintAst`.
6. **Coercions** — `*Arg` unions where a field takes literal-or-mirror (e.g. `order`:
   `int | ValueAst`), plus `XConstraintsArg` (value|view snapshot) and `ConstraintsUpdate`.

Structural variations from the atom template:

- **Participants are read-only.** Atoms had no participants and fully-settable fields. Every
  other entity's view exposes its endpoints/participants (`atom_ids` / `acceptor`+`donors` /
  `site`+`ligands`) as **reads only** — topology editing is out of scope here (it is a
  builder concern). Settable = the chemistry value fields + constraints only.
- **Two addressing modes.** Bonds are dense-id-indexed (like atoms) **and** edge-indexed by
  atom pair (`of(first, second)`). The five relation families are dense-id-indexed **and**
  participant/incidence-addressed (`of(...)`, `incident(atom)`, `induced(atoms)`;
  stereo adds `coincident(site)` under the ≤1-per-site invariant).
- **Whole-value replace via `__setitem__`.** `mol.bonds[i] = BondAst(...)` (and the relation
  analogues) replace the value while keeping topology, since the value carries no endpoints.
  Implemented through `x_mut(id)`.
- **Molecule construction is a prerequisite for every view half.** `MoleculeAst.from_atoms`
  takes atoms only; to exercise any `XView` a molecule must carry bonds/relations. How
  Python builds topology — mirror `from_atoms_and_bonds`, add `add_bond`/`add_*` methods, or
  a builder — is the top cross-cutting **design call** (below). The standalone value type
  `XAst` needs no molecule and can land first for every entity.
- **Constraint containers vary in richness:** bond (3 keys), dative (2), aromatic /
  multicenter (1: `ElectronCount`), noncovalent (none — uninhabited today), stereo (4 keys +
  overlay specifics).

### Participant lookup and how it is expressed

`mol.<family>` is an id-indexed sequence of views (like `mol.atoms`). Participant lookup is
**methods on the collection**, not the indexing scheme (participants are a secondary index,
and unordered/variable factors make poor `__getitem__` keys):

- `of(…) -> XView | None` — the exact-participants lookup ("the entity *of* these
  participants"; `find_by_participants` / `find_edge`; Rust rename per 137 item 4).
- `incident(atom) -> list[XView]` — relations touching an atom.
- `induced(atoms) -> list[XView]` — relations wholly within an atom set.
- stereo only: `coincident(site) -> XView | None` (the ≤1-per-site lookup).

No subscript sugar (`mol.bonds[a, b]`) — only `[id]` plus these methods. `incident`/`induced`
return **lists** (results are small).

`of`'s call convention tracks the factor shape: positional for a fixed symmetric pair
(`bonds.of(0, 1)`), one iterable for a variable set (`aromatic_systems.of({…})`), keyword roles
for a birelation (`dative_bonds.of(acceptor=…, donors=…)`, `stereo_atoms.of(site=…,
ligands=[…])`).

**Lookup inputs are order-free but multiplicity-sensitive.** `Graph::find_edge` is
order-independent (undirected adjacency is symmetric) and `find_by_participants` sorts each
factor, so argument order never matters — but the match is a **multiset** (length included).
Pass any iterable; use a `set` only where participants are distinct (atom sets, donors). The
`ligands` arg must be an ordered sequence (list/tuple), not a `set` — a stereo frame can hold
*repeated* virtual ligands (two implicit H / lone pairs are identical `StereoLigand` values),
which a set would dedupe into a false miss.

**The relation/birelation split is not a Python concept.** It surfaces only as (a) `of`'s
call shape — one participant group for relations (`of(atoms)`), two role groups for
birelations (`of(first, second)` for bonds/noncovalent, `of(acceptor=…, donors=…)` /
`of(site=…, ligands=…)`) — and (b) the view's participant accessors below.

### Participant return types (reads on the view)

The type mirrors the factor's arity/ordering, per role:

| Entity | participant read → Python type |
|---|---|
| bond, noncovalent | `atom_ids -> tuple[int, int]` (**sorted**; the pair is undirected) |
| aromatic, multicenter | `atom_ids -> frozenset[int]` (variable, unordered) |
| dative | `acceptor -> int`, `donors -> frozenset[int]`; role composite as a top-level tuple `(donors, acceptor)` |
| stereo-atom / -bond | `site -> int`, `ligands -> tuple[StereoLigand]` (ordered frame); role composite as a top-level tuple `(site, ligands)` |

Rule: fixed unordered pair → **sorted tuple**; variable unordered set → **frozenset**; ordered
sequence → **tuple**; arity-1 role → **bare int**; a birelation's by-role whole → **top-level
tuple** of the roles. (This refines the earlier "ordering marker picks the type" note: the
fixed pairs read as a sorted tuple, not a frozenset, for ergonomics.)

## B. Per-entity specifics

### B1 · Bond — defines the template

- Value `BondAst { order: ValueAst, charge: ValueAst, spin: SpinStateAst, constraints }` —
  all leaves already bound; reuse `ValueAst`/`SpinStateAst`.
- Constraints `BondConstraintAst { Aromatic(BooleanAst), CisTransStereo(CisTransStereoAst),
  RingMembership(RingMembershipAst) }`; key `{ Aromatic, CisTransStereo,
  RingMembership(RingScope) }`. **New leaves: `BooleanAst`, `CisTransStereoAst`.** `RingMembership`
  is a keyed multimap over `RingScope` (reuse the atom `ring_count` / `ring_size_count[…]`
  proxy pattern verbatim).
- View `BondView`: read-only `atom_ids -> (int, int)`; settable `order`/`charge`/`spin`;
  `constraints`. Collection `mol.bonds`: id-indexed sequence + `connecting(a, b)` +
  `__setitem__` (value replace).
- Design calls: expose `spin` on a bond (unusual); `CisTransStereo` surface (reuse the atom
  `TetrahedralStereo.Cw/Ccw` enum idea as a cis/trans enum, or expose the mirror). `order`/
  `charge` accept `ValueArg` (already solved for atoms).

### B2 · Dative

- Value `DativeBondAst { order: ValueAst, constraints }` — no charge, no spin.
- Constraints `{ Aromatic(BooleanAst), RingMembership(RingMembershipAst) }`; key
  `{ Aromatic, RingMembership(RingScope) }`. No new leaves (`BooleanAst` from B1).
- View `DativeBondView`: read-only `acceptor -> int` + `donors -> list[int]` (unordered set;
  N>1 = haptic), `atom_ids` = donors-then-acceptor; settable `order`; `constraints`.
  Collection `mol.dative_bonds`: id-indexed + `connecting(donors, acceptor)` / `incident`.
- Design calls: the acceptor/donor split surface (a scalar `acceptor` + a `donors` list),
  and that donors are an unordered multiset.

#### B2 · Dative — staged impl plan (settled 2026-07-11)

Two design calls resolved: constructor is **atoms-positional** `from_parts(atoms, *, bonds=[],
dative=[])`; dative constraint surface is now a clean mirror of bond's (the ring-accessor optionality
discrepancy was fixed on the Rust side — `aromatic() -> BooleanAst` non-optional, `ring_count()`/
`ring_size_count() -> Option<&ValueAst>`). Dative-entry tuple is `(donors: list[int], acceptor: int,
DativeBondAst)` (donors-first).

**Blocker found + folded in:** the overlay `*_mut(id)` accessors return a bare `&mut XAst` (no view
guard) — non-uniform with `atom_mut`/`bond_mut` and the 114 interning-guard hole (doc 137 pt6). The
overlay Python view halves would edit through that bare ref. So **Part A (Rust) implements 137 pt6 in
full** as the prerequisite, then **Part B** binds dative. Feasibility of `*ViewMut`: confirmed — the
relation sets keep `participants*(&self)` and `data_mut(&mut self)` reachable in sequence, so the
accessor reads participants → copies to owned → the `&self` borrow ends → `data_mut` (as `BondViewMut`
already owns its `[AtomId;2]`).

**Part A — uniform `*ViewMut` + retire bulk `*s_mut()` (Rust-internal; implements doc 137 pt6). No umol-py.**

- **S0a — DONE** — six `*ViewMut` structs added (view/*.rs), exported via view.rs + ast.rs.
- **S0b — DONE** — the six `MoleculeAst::X_mut(id)` return their `*ViewMut` (`Arc::make_mut` → read
  participants to owned → `data_mut`); the **29 `.ast` sites** migrated (molecule.rs 13, aromaticity.rs 2,
  validate/stereo.rs 7, symmetry.rs 1, molecule/tests.rs 6). The 15 `MoleculeEditor`/`transact.rs` +
  editor-test hits already used `.ast`.
- **S0a/b consolidation — DONE** — a full **view-record convention** pass (surfaced reviewing the new
  ViewMuts) unified every `*ViewMut`/`*EditorView`/`*EditorViewMut` across all entities: field order
  **`id, parts, ast`** (`molecule` last on read `*View`); **no `pub(crate)` fields anywhere**; **pub-field
  bundles, struct-literal construction, no `new`, no accessors** for `*ViewMut` and the fixed/stereo
  editor views. The read `*View` is untouched (gold standard: `pub id/ast`, private raw participants +
  accessor, private `molecule`). The three variable-arity editor views (dative/aromatic/multicenter)
  keep a **private borrowed `&[NodeId]`** (no clone on the edit path) behind `atom_ids()` + a
  `pub(crate) fn new` — the only `pub(crate)` left, on a method not a field. `BondViewMut.atoms()`
  misnaming removed. Workspace 10,851 tests green, clippy + fmt clean.
- **S0c — DONE** — eight `MoleculeAst::map_*(&mut self, f: impl FnMut(XAst) -> XAst)` added (atoms/bonds
  + 6 relations; body `for slot in <iter>_mut() { *slot = f(mem::take(slot)); }` — no `&mut XAst`
  escapes; the container owns re-interning post-`f`, the 114 hook). The **22 bulk sites** rewritten
  (dsl/molecule.rs 16 raise/lower loops → `map_*(|x| …)` one-liners; resolve 2 + tests 4 → `map_*(|mut x|
  { …; x })`); the eight `Xs_mut()` deleted; four stale `test_*_mut` fns renamed `test_*_map_*`.
- **S0d — DONE** — workspace build 0/0, clippy + fmt clean, 10,851 tests green. **Part A (doc 137 pt6)
  complete.** Next: **Part B (dative Python slice), S1a.**
- Interning (114): the `*ViewMut` is the future `Drop` re-intern-guard slot — **no `Drop` added now**;
  `map_*` is the container's post-`f` re-intern hook. Both are prerequisite shape only.

**Part B — Python dative slice (`umol-py/src/dative.rs`), mirrors `bond.rs` minus CisTransStereo + charge/spin. After S0.**

- **S1a** *(additive/green)* — `DativeBondAst` value pyclass (`#[new] new(order: ValueArg, *,
  constraints=None)`; `parse`/`__str__`/`__repr__`; order+constraints getters/setters; `asdict`;
  `inner`/`inner_mut`/`from_inner`) + the WET constraint surface (verbatim bond mirror, 2 keys):
  `DativeBondConstraintKey`/`DativeBondConstraintAst` (`Aromatic`, `RingMembership`),
  `DativeBondConstraintsAst` (whole mapping API), `DativeBondConstraintsView`,
  `DativeBondConstraintsBacking { Molecule{owner, id: DativeBondId}, DativeBond(Py<DativeBondAst>) }`,
  `DativeBondRingSizeCounts`/`DativeBondRingSizeBacking`, the 3 iterators,
  `DativeBondConstraintsArg`/`Update`, `dative_bond_constraints_asdict`. Reuse `BooleanAst`/
  `RingMembershipAst`/`RingScope`/`ValueAst` leaves. The Molecule-backed view's write-through is
  `mol.inner_mut().dative_bond_mut(id).ast.constraints` — the `.ast` is S0's ViewMut. Register value +
  constraint pyclasses in `lib.rs` + `__init__.py`; Rust unit tests. `[dep: S0a]`
- **S1b** *(breaking→green)* — rename Python `MoleculeAst.from_atoms_and_bonds` → `from_parts(atoms,
  *, bonds=[], dative=[])`; wire `dative` → `MoleculeParts.dative` (entries `(donors, acceptor,
  Py<DativeBondAst>)`). Migrate **27 Python test sites** (test_atom/test_constraint/test_molecule/
  test_bond) + the umol-py Rust test call sites. `[dep: S1a — dative entries carry Py<DativeBondAst>]`
- **S1c** *(additive/green)* — `DativeBondView` (`id`, `acceptor -> int`, `donors -> list[int]`
  read-only, `atom_ids -> tuple` donors-then-acceptor, settable `order`/`constraints`, `asdict`) +
  `DativeBondViews` (`mol.dative_bonds`: `__len__`/`__getitem__`/`__setitem__` value-replace/`__iter__`
  + `connecting(donors, acceptor)` + `incident(atom)`) + `DativeBondViewIter` +
  `resolve_dative_bond_index` (DativeBondId is RelationId-backed but contiguous for fresh molecules →
  mirror `resolve_bond_index`). Register view pyclasses; Rust unit + `tests/test_dative.py`; maturin
  rebuild + pytest. `[dep: S1a, S1b]`

Critical path: S0a → S0b → S1a → S1b → S1c (S0c parallel to S0b; S0d gates Part B). Part A unblocks
**every** remaining overlay view half (B3–B6 + stereo), not just dative.

### B3 · Aromatic system

- Value `AromaticSystemAst { electrons: ElectronCountsAst, charge: ValueAst, spin:
  SpinStateAst, constraints }`.
- Constraints `{ ElectronCount(ValueAst) }`; key `{ ElectronCount }` (single-valued).
- **New leaf: `ElectronCountsAst { Undetermined | Lit(list[int]) }`** — a per-participant
  positional vector over the (unordered) atom set.
- View: read-only `atom_ids` (variable set) + derived `bond_ids`; settable `electrons`/
  `charge`/`spin`; `constraints`. Collection `mol.aromatic_systems`: id-indexed +
  `incident` / `induced`.
- Design calls: `electrons` is a **positional** `list[int]` aligned to an **unordered**
  participant set — how each count pairs with its atom in Python (a list aligned to
  `atom_ids`, or a `{atom: count}` mapping). This is the central aromatic/multicenter call.

### B4 · Multicenter bond

- **Structurally byte-for-byte identical to aromatic** (`electrons`/`charge`/`spin`/
  `constraints`; `ElectronCount(ValueAst)` constraint). One template serves both; only the
  type names differ. Same `electrons`-pairing design call as B3.

### B5 · Noncovalent bond

- Value `NoncovalentBondAst { kind: NoncovalentBondKindAst, constraints }`.
- **New leaves: `NoncovalentBondKindAst { Undetermined | Lit(NoncovalentBondKind) }` and
  `NoncovalentBondKind { HydrogenBond, HalogenBond, ChalcogenBond, Ionic, VanDerWaals }`**
  (a fieldless value enum → hashable, per finding 137-p3-2).
- Constraints: the enum and key are **uninhabited today** (scaffolding only). Design call:
  omit the constraints surface entirely (no container/view/key), since no variant can exist,
  vs. bind an always-empty stub.
- View: read-only atom pair (fixed 2, unordered); settable `kind`. Collection
  `mol.noncovalent_bonds`: id-indexed + `connecting(a, b)`.
- Design calls: `kind` wildcard (`Undetermined`) → an optional/nullable surface + coercion;
  the uninhabited-constraints decision above.

### B6 / B7 · Stereo atom / stereo bond (the overlay — a larger sub-project)

- Value `StereoAtomAst` / `StereoBondAst { configuration: StereoConfigurationAst, constraints:
  Stereo{Atom,Bond}ConstraintsAst }`. Site and ligands are **not** in the value — they are
  the birelation's factor-1 (site) and factor-2 (ordered `StereoLigand` frame).
- **New leaves:** `StereoConfigurationAst { Undetermined | Kinded(StereoKind, StereoCosetAst)
  }`; `StereoKind { Tetrahedral, CisTrans, Axial, SquarePlanar, TrigonalBipyramidal,
  Octahedral }`; `StereoLigand { atom_id, kind: StereoLigandKind }`; `StereoLigandKind { Atom,
  ImplicitHydrogen, LonePair }`; `StereoLigandPosition(u32)`; the stereo constraint values
  `LigandSymmetryAst` / `FluxionalityAst` / `TopicityAst` / `StereogenicityAst` and their
  helper leaves (`OrientedLigandPermutation`, `LigandPermutation`, `StereoLigandPair`,
  `TopicityRelationAst`); the derived `Topicity` / `Stereogenicity` enums; `StereoConfiguration`
  (ground). `StereoCosetAst` / `StereoTerm` / `Permutation` are already bound.
- Overlay structure: keyed on site (atom `NodeId` / bond `EdgeId`) + ordered `StereoLigand`
  frame; ≤1 per site; the coset is stated **relative to the stored ligand frame**.
- View `StereoAtomView` / `StereoBondView`: read-only `site_id` (atom / bond) + ordered
  `ligands` (incl. virtual) + `kind` + `coset`; settable `configuration`; `constraints`.
  Collection: id-indexed + `coincident(site)` + `connecting(site, ligands)`.
- Design calls (many — see the stereo block in Part C). This slice is large enough to warrant
  its own doc/staging after the covalent entities land; it also coincides with the deferred
  finding 137-p3-4 (`Permutation.image` shape) "top-level stereo" revisit.

## C. Leaf types needing design calls

| Leaf | Used by | Rust shape | The call |
|---|---|---|---|
| `BooleanAst` | bond, dative (`Aromatic`) | `Undetermined \| Lit(bool)` | **Resolved (2026-07-09):** full `BooleanAst` mirror (`Undetermined`/`Lit(bool)`) is the representation; `True`/`False` coerces on assignment (a `BooleanArg`: `bool → Lit`, mirror passthrough). |
| `CisTransStereoAst` | bond (`CisTransStereo`) | `Undetermined \| NotStereo \| Stereo(StereoCosetAst)` | **Resolved:** full mirror (like `TetrahedralStereoAst`) **plus** a simple `CisTransStereo` enum `.Z`/`.E` for assignment, mirroring `TetrahedralStereo.Ccw/Cw`. `Z → Lit(0)`, `E → Lit(1)` (dsl/stereo.rs:336ff). |
| `ElectronCountsAst` | aromatic, multicenter | `Undetermined \| Lit(list[int])` | **Resolved:** mirror whose `Lit` variant is a Python `list[int]` (positional, aligned to participant order); no `{atom: count}` dict. |
| `NoncovalentBondKindAst` + `NoncovalentBondKind` | noncovalent | mirror `Undetermined \| Lit(NoncovalentBondKind)`; fieldless enum `{HydrogenBond, HalogenBond, ChalcogenBond, Ionic, VanDerWaals}` | **Resolved (2026-07-09):** names follow Rust **exactly** (`NoncovalentBondKind.HydrogenBond`, no shortening); a simple hashable value enum (137-p3-2 pattern). `kind` is the full `NoncovalentBondKindAst` mirror (`Undetermined` as a variant, not `None` — the Rust structure, matching `atom.element: ElementAst`), assigned via the bare enum. |
| noncovalent `constraints` | noncovalent | uninhabited enum/key | **Resolved:** empty stub — an always-empty container/view (`len 0`, empty iteration, `asdict → {}`, no inhabited keys). |
| `StereoConfigurationAst` / `StereoKind` | stereo | `Undetermined \| Kinded(kind, coset)`; 6-kind enum | The coset carries no kind — configuration (kind+coset) must be the surfaced unit; any coset op threads `StereoKind`. |
| `StereoLigand` / `StereoLigandKind` | stereo | `{atom_id, kind}`; `{Atom, ImplicitHydrogen, LonePair}` | Surface `StereoLigand{atom_id, kind}` vs. bare atom ids; how virtual ligands (bearing-atom id) appear. |
| `StereoLigandPosition` (vs `ParticipantPosition`) | stereo | `u32` frame index | Which position newtype (if any) surfaces to Python. |
| Stereo constraint values (`LigandSymmetryAst`/`FluxionalityAst`/`TopicityAst`/`StereogenicityAst` + helpers) | stereo | per-site constraint bag, 4 keys | A whole sub-surface — bind now vs. defer within the stereo slice. |
| Stereo overlay addressing | stereo | id vs `(site, ligand-multiset)` content lookup (order-independent) | How a relation is addressed from Python; whether ligand order is supplied for reframe (`permutation_for`/`coset_for`/`transform_frame`); note this is handle-identity-adjacent (interning, 114). |

The covalent leaf calls are resolved (2026-07-09); **all stereo rows remain open**, deferred
to the stereo slice (its own doc), alongside the overlay structure and the constraint
sub-surface.

**Naming: internal type/enum/member names follow the Rust side exactly** (`NoncovalentBondKind.
HydrogenBond`, never a shortened `HBond`). The only Python-only additions are the assignment-
convenience enums with no Rust equivalent — `TetrahedralStereo.Ccw/Cw` and `CisTransStereo.Z/E`
(named for the chemistry / DSL keywords) — plus the `E` element shorthand. A terse `E.H`-style
shorthand for enum members (e.g. the noncovalent kinds) is a possible future ergonomic —
**deferred** ("can wait").

**Cross-cutting call #1 (gates every view half): molecule topology construction — resolved
via the Rust edit-API slice (086, 2026-07-09).** The Rust side gets a uniform
`MoleculeParts` constructor (`#[derive(Default)]`, one field per family — no atoms/bonds
privilege), `MoleculeAst::join(other) -> Remapping` (disjoint concatenation; result may be
disconnected), and `MoleculeAst::split() -> Vec<MoleculeAst>` (partition by
connected-via-any-relation; conservative — remove a bond first to split finer);
`from_atoms_and_bonds`, positional `from_parts`, `has_overlays`, and `is_in_overlays` are
retired. Python mirrors: `MoleculeParts`-style construction, `mol.join(other)`,
`mol.split()`, and none of the overlay-privileging surface. The incremental `MoleculeBuilder`
(uniform per-family `add_*`) stays.

**`MoleculeParts` + the overlay-surface removal landed in Rust 2026-07-11** (086 L2199-2211 — the
two pieces that unblock every remaining view half). Consequence for this plan: each entity's view
half is now buildable — construct a molecule carrying that entity via `MoleculeAst::from_parts` and
read/write it through the view. The Python molecule constructor currently exposes only atoms+bonds
(`MoleculeAst.from_atoms_and_bonds`, routed through the Rust `from_parts`); it grows one keyword-arg
per family (and is renamed to the `from_parts`/`MoleculeParts` mirror, de-privileging atoms/bonds)
as each entity's Python binding lands — do that rename at the **dative** slice (B2), the first
overlay entity, not speculatively now. `join`/`split` Python mirrors remain future work and are not
prerequisites for the per-entity view halves.

## Build shape (once the calls are settled)

Not a staged plan yet — the design calls above must resolve first (per the design-before-plan
convention). The shape it will take, following the given order and the staged-impl-plan rule
(additive first, breaking/rewiring last):

- Standalone value types (`XAst`) + their new leaves + constraint containers land first for
  each entity — additive, no molecule dependency, testable in isolation (as `AtomAst` was).
- The molecule-construction surface (cross-cutting call #1) lands next; it unblocks all view
  halves.
- View halves (`XView`/`XViews`) follow per entity, in the given order — bond first (defines
  the pattern), then dative, then aromatic + multicenter (one shared template), then
  noncovalent (value-only, no constraints).
- The two stereo overlays are a separate slice/doc after the covalent entities, carrying the
  bulk of the design calls (frame-relative cosets, ligand frames, the constraint sub-surface),
  and align with the deferred "top-level stereo" work.

## Bond slice — build state + detailed plan (2026-07-11)

Pick-up point after a manual compaction. Building **B1 · Bond** as the additive value-type slice
(value + new leaves + constraint surface; the `BondView`/`BondViews` molecule-handle half is deferred
on the topology-construction prerequisite).

### Factoring decision (resolved)

The constraint sub-surface is the one real call. Chosen: **generalize nothing structural now — WET
per-entity mirror**, exactly as umol-ast eschews DRY. Reason: cleanly generalizing the *backing* +
`ring_size_count` proxy needs a `ConstraintFamily` trait, a `dyn`-backed proxy, **and** retrofitting the
working atom binding — the "quite large / future pass" case the user named. Iterators are the trivial
first target when that future generalization pass happens. (User: "if in doubt, follow umol-ast WET.")

### Done (compiling green; `cargo build -p umol-py --features graph`)

- `umol-py/src/boolean.rs` — `BooleanAst { Undetermined() | Lit(bool) }` + `BooleanArg` (`bool → Lit`,
  mirror passthrough). Leaf of the bond/dative `Aromatic` constraint.
- `umol-py/src/stereo.rs` — `CisTransStereoAst { Undetermined() | NotStereo() | Stereo(Py<StereoCosetAst>) }`
  (line-for-line copy of `TetrahedralStereoAst`, `stereo.rs:199`) + `CisTransStereo { Z, E }` (`Z → Ct0`,
  `E → Ct1`; copy of `TetrahedralStereo`, `stereo.rs:250`).
- `lib.rs` — `mod boolean`; `use` + `add_class` for `BooleanAst` / `CisTransStereoAst` / `CisTransStereo`.
- Transient dead-code warnings on the leaves' `from_ast`/`to_ast`/`BooleanArg` — they feed the surface;
  clear when it lands. Not a defect (the slice is atomic: leaves + value + constraints interdepend).

### Rust API being bound (names verified 2026-07-11)

- `BondAst` (`umol-ast/src/ast/bond.rs`): `{ order: ValueAst, charge: ValueAst, spin: SpinStateAst,
  constraints: BondConstraintsAst }`; `new(order: ValueAst)`, `from_order(u8)`,
  `with_order/with_charge/with_spin(impl Into<…>)`, `From<&str>`/`FromStr` (DSL) + `Display`. → mirror the
  `AtomAst` value type (3 value fields, `parse`/`__str__`/`__repr__`).
- `BondConstraintAst` (`umol-ast/src/ast/constraint/bond.rs`): `{ Aromatic(BooleanAst),
  CisTransStereo(CisTransStereoAst), RingMembership(RingMembershipAst) }`; ctors `aromatic(impl Into<
  BooleanAst>)` (l.24), `cis_trans_stereo(impl Into<CisTransStereoAst>)` (l.28), `ring_membership(RingScope,
  impl Into<ValueAst>)` (l.32).
- `BondConstraintKey`: `{ Aromatic, CisTransStereo, RingMembership(RingScope) }`.
- `BondConstraintsAst` — container: `new`, `set(c)` (l.219), `remove(key)` (l.263), `get(key)` (l.214),
  `contains(key)` (l.210), `len` (l.202), `iter` (l.301), `extend`/`update` (mirror atom). Per-key reads:
  **`aromatic() -> BooleanAst` (by value, NON-optional — l.163; differs from the atom `Option<&…>` pattern,
  so the Python `aromatic` getter is non-optional, returning `Undetermined` when unset)**;
  `cis_trans_stereo() -> Option<&CisTransStereoAst>` (l.170); `ring_count() -> Option<&ValueAst>` (l.190);
  `ring_size_count(u8) -> Option<&ValueAst>` (l.194).

### Reuse (do NOT re-bind)

`ValueAst`/`ValueArg` (value.rs), `SpinStateAst` (atom.rs), `RingScope`/`RingMembershipAst`
(constraint.rs), `BooleanAst`/`BooleanArg` (boolean.rs), `CisTransStereoAst`/`CisTransStereo` (stereo.rs).

### Templates to copy (rename Atom→Bond, keep only the 3 keys)

- Value: `atom.rs` `AtomAst` block (`#[pyclass(eq)]` wrapper; `#[new]` kw-only fields; `parse`/`__str__`/
  `__repr__`; getters/setters for order/charge/spin; `constraints` getter→view / setter→`BondConstraintsArg`;
  `asdict`; `apply_fields`; `inner`/`inner_mut`/`#[cfg(test)] from_inner`). `order`/`charge` take `ValueArg`.
- Constraint surface: `constraint.rs` — `AtomConstraintAst`→`BondConstraintAst` (3 variants),
  `AtomConstraintKey`→`BondConstraintKey` (3), `AtomConstraintsAst`→`BondConstraintsAst` (whole mapping API:
  `new`/`__repr__`/`set`/`pop`/`update`/`__len__`/`__iter__`/`keys`/`values`/`items`/`get`/`__getitem__`/
  `__delitem__`/`__contains__` + per-key `aromatic`/`cis_trans_stereo`/`ring_count` getters+setters +
  `ring_size_count` proxy getter + `asdict`), `AtomConstraintsView`→`BondConstraintsView`,
  the 3 iterators, `RingSizeCounts`/`RingSizeBacking`→bond copies, `atom_constraints_asdict`→bond,
  `ConstraintsArg`/`ConstraintsUpdate`→`BondConstraintsArg`/`BondConstraintsUpdate`.

### Steps

1. `CisTransStereoArg` in `stereo.rs` — copy `TetrahedralStereoArg` (constraint.rs:176): `False →
   NotStereo`, `CisTransStereo → coset`, `CisTransStereoAst` passthrough. (For the `cis_trans_stereo` setter.)
2. `umol-py/src/bond.rs` (new) — `BondAst` value + the whole bond constraint surface (keep WET-separate from
   the atom surface). **`BondConstraintsBacking` = `Bond(Py<BondAst>)` only for now** — add the
   `Molecule{owner, id: BondId}` arm with the deferred `BondView` half (else dead-code). Same for the bond
   `RingSizeBacking`: `Bond(...)` + `Value(Py<BondConstraintsAst>)` now, `Molecule{id: BondId}` later.
3. Register every bond pyclass in `lib.rs`.
4. Tests (in `bond.rs`, mirror the constraint.rs tests): `BondConstraintAst` roundtrip; `BondConstraintsAst`
   len/contains/keys/values/items/get; the `ring_size_count` proxy; `BondAst` field round-trip + `parse`.
   Confirm the leaf dead-code warnings are gone.
5. **umol-py tests SIGABRT on libpython** — use the venv: `source umol-py/.venv/bin/activate` (see memory).

### View half — LANDED (2026-07-11), not deferred after all

The AST fully supports bonds by stable `BondId` today (`molecule.bonds()`, `bond_mut(id) -> BondViewMut
{ ast: &mut BondAst }`, `BondView { ast, atom_ids() }`), so nothing was actually blocked. The only real
gate was that Python could not construct a molecule *with* bonds — resolved by replacing the `from_atoms`
backstop with **`MoleculeAst.from_atoms_and_bonds(atoms, bonds=[(i, j, BondAst), …])`** (the exact Rust
counterpart). Landed:

- `molecule.rs`: `from_atoms_and_bonds` (replaces `from_atoms`; all 23 Python call sites migrated),
  `bonds` getter → `BondViews`, `__repr__` now `atoms=N, bonds=M`.
- `bond.rs`: `BondView` (id, read-only `atom_ids -> (int,int)`, order/charge/spin get+set, constraints
  get+set, asdict, repr) + `BondViews` (`mol.bonds`: len/getitem/setitem-value-replace/iter +
  `connecting(a,b)`) + `resolve_bond_index` + `BondViewIter`. `Molecule{owner, id: BondId}` arm added to
  both `BondConstraintsBacking` and `BondRingSizeBacking`.
- Registered `BondView`/`BondViews` in `lib.rs` + `__init__.py`.
- Tests: +11 Rust view tests in `bond.rs`, +18 Python tests in `test_bond.py`, `test_molecule.py`
  updated (repr + `from_atoms_and_bonds` construction).

Full parity with the atom slice reached. `MoleculeParts`/`join`/`split` (cross-cutting call #1's broader
resolution) remain future work but are **not** a prerequisite for the per-entity view halves — the peers
(dative/aromatic/…) follow this same `from_atoms_and_bonds`-backed pattern.

### Still deferred

- Structural generalization of backing/iterators/proxy — the future DRY pass (WET per-entity for now).

### Landed (2026-07-11) — bond slice complete

- `umol-py/src/bond.rs` (new): `BondAst` value type (`order`/`charge`/`spin`/`constraints`, `parse`/
  `__str__`/`__repr__`/`asdict`) + full WET constraint surface — `BondConstraintAst`/`BondConstraintKey`
  (3 keys), `BondConstraintsAst` (whole mapping API), `BondConstraintsView`, 3 iterators,
  `BondRingSizeCounts`/`BondRingSizeBacking`, `bond_constraints_asdict`, `BondConstraintsArg`/`Update`.
  `aromatic` getter is **non-optional** (returns `Undetermined` when unset), matching the Rust accessor.
  Backing = `Bond(Py<BondAst>)` only; molecule arm deferred with `BondView`.
- `stereo.rs`: added `CisTransStereoArg` (`False → NotStereo`, `CisTransStereo → coset`, passthrough).
- Naming call (user): symmetric rename `RingSizeCounts`→`AtomRingSizeCounts` (+`Backing`/`Iter`); bond
  copy is `BondRingSizeCounts`.
- `lib.rs` + `python/umol/__init__.py`: registered/exported all bond pyclasses **and** the previously
  unexported leaves `BooleanAst`/`CisTransStereo`/`CisTransStereoAst`.
- Tests: 34 Rust unit tests in `bond.rs`; Python `tests/test_bond.py` (mirrors test_atom + test_constraint,
  standalone-backed only) + `CisTransStereo*` cases in `tests/test_stereo.py`. Rust suite 190 green;
  Python suite 219 green; clippy clean.
- Note for the peers: the bond `RingScope` proxy is per-entity (WET). When the future DRY pass lands, the
  iterators are the trivial first target; the backing enum + `ring_size_count` proxy are the hard part.
