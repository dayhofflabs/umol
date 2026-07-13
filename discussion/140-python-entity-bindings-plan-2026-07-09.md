# 140 · Python bindings for the remaining entity ASTs (plan)

Status: Active — **B1 · Bond slice DONE incl. the view half** (value + WET constraint surface +
`BondView`/`BondViews` + molecule-backed constraint views; Rust/Python green, clippy clean).
**B2 · dative COMPLETE**; **B3 · aromatic COMPLETE** (+ new `ElectronCountsAst` leaf, `electrons.rs`);
**B4 · multicenter COMPLETE** (`umol-py/src/multicenter.rs`: value + WET constraint surface +
`MulticenterBondView`/`MulticenterBondViews` + molecule-backed constraint views; `mol.multicenter_bonds`
accessor; `from_parts` `multicenter=[]` kwarg; reuses the `ElectronCountsAst` leaf, no new foundation;
9 multicenter pyclasses registered; 344 Rust unit + 332 pytest green, clippy/fmt clean). **B5 ·
noncovalent COMPLETE** (2026-07-12): the uninhabited-constraint + pyo3-zero-variant-enum blocker was
resolved by inhabiting the constraint upstream (`Intramolecular(BooleanAst)`, `#I` — doc 117 §4 stages
A–C), then binding the full slice (`umol-py/src/noncovalent.rs`: kind leaf, 1-key Boolean constraint
container + live view, `NoncovalentBondAst` value pyclass, `NoncovalentBondView`/`NoncovalentBondViews`,
`mol.noncovalent_bonds` accessor, `from_parts(noncovalent=[])` kwarg). 417 Rust unit + 363 pytest green,
clippy/fmt clean. **A concurrent self-alias `RefCell` double-borrow panic** (surfaced by an adversarial
review of the noncovalent slice) was fixed in noncovalent AND swept across all peer slices (atom/bond/
dative/aromatic/multicenter) — resolve-before-borrow; regression tests per entity. The atom-constraint
container/view machinery was also moved `constraint.rs` → `atom.rs` (renamed `Constraints*` →
`AtomConstraints*`); `constraint.rs` is now the shared constraint value/scope leaves.
See the *B5 · Noncovalent bond* section below, doc 117 §4, and *B3 · Aromatic* / *B4 · Multicenter —
staged impl plan* for the WET template. **Remaining: B6/B7 · stereo** (the overlay sub-project).
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
  complete.** Next: **Part B (dative Python slice), S1a.** *(S1a now DONE — see below; next is S1b.)*
- Interning (114): the `*ViewMut` is the future `Drop` re-intern-guard slot — **no `Drop` added now**;
  `map_*` is the container's post-`f` re-intern hook. Both are prerequisite shape only.

**Part B — Python dative slice (`umol-py/src/dative.rs`), mirrors `bond.rs` minus CisTransStereo + charge/spin. After S0.**

- **S1a — DONE** *(additive/green)* — `DativeBondAst` value pyclass (`#[new] new(order: ValueArg, *,
  constraints=None)`; `parse`/`__str__`/`__repr__`; order+constraints getters/setters; `asdict`;
  `inner`/`inner_mut`/`from_inner`) + the WET constraint surface (verbatim bond mirror, 2 keys):
  `DativeBondConstraintKey`/`DativeBondConstraintAst` (`Aromatic`, `RingMembership`),
  `DativeBondConstraintsAst` (whole mapping API), `DativeBondConstraintsView`,
  `DativeBondConstraintsBacking`, `DativeBondRingSizeCounts`/`DativeBondRingSizeBacking`, the 3
  iterators, `DativeBondConstraintsArg`/`Update`, `dative_bond_constraints_asdict`. Reuse `BooleanAst`/
  `RingMembershipAst`/`RingScope`/`ValueAst` leaves. Registered the 6 value + constraint pyclasses in
  `lib.rs` + `__init__.py`; 29 Rust unit tests; `import umol` + full 237-test pytest suite green.
  **Staging correction:** the `Molecule` backing arm of both `DativeBondConstraintsBacking` and
  `DativeBondRingSizeBacking` moved to **S1c** — its only constructor is `DativeBondView::constraints`
  (S1c), so at S1a it would be an unconstructed variant (dead-code warning in the non-test lib build;
  bonds avoided this only because `BondView` shared the stage). S1a backings are `DativeBond`-only
  (`DativeBondRingSizeBacking` also keeps `Value`); S1c additively re-adds `Molecule` with the
  molecule-backed write-through `mol.inner_mut().dative_bond_mut(id).ast.constraints` (S0's ViewMut)
  and the two molecule-backed unit tests. `[dep: S0a]`
- **S1b — DONE** *(breaking→green)* — renamed Python `MoleculeAst.from_atoms_and_bonds` →
  `from_parts(atoms, *, bonds=[], dative=[])` (`bonds` now keyword-only); wired `dative` →
  `MoleculeParts.dative` (Python entry `(donors: list[int], acceptor: int, DativeBondAst)` →
  `(Vec<AtomId>, AtomId, DativeBondAst)`). Migrated the 29 Python test-file occurrences across
  test_atom/test_constraint/test_molecule/test_bond (global rename + `bonds=` at the 3
  positional-bonds sites); no Rust test called the pymethod (test modules use the AST-side `from_parts`
  directly). Added a Rust unit test `test_molecule_ast_from_parts` covering the new `dative` wiring,
  asserted via `inner().dative_bonds()` (Python has no `mol.dative_bonds` accessor until S1c, so the
  dative path has no Python observable yet — its pytest coverage lands with `tests/test_dative.py` in
  S1c). Rust 231 unit + 237 pytest green, clippy/fmt clean. `[dep: S1a]`
- **S1c — DONE** *(additive/green)* — `DativeBondView` (`id`, `acceptor -> int`, `donors -> list[int]`
  read-only, `atom_ids -> tuple` donors-then-acceptor, settable `order`/`constraints`, `asdict` =
  {order, constraints}) + `DativeBondViews` (`mol.dative_bonds`: `__len__`/`__getitem__`/`__setitem__`
  value-replace/`__iter__` + `connecting(donors, acceptor)` + `incident(atom)`) + `DativeBondViewIter`
  + `resolve_dative_bond_index` (mirrors `resolve_bond_index`). Re-added the `Molecule` backing arm to
  `DativeBondConstraintsBacking` + `DativeBondRingSizeBacking` (deferred from S1a) — `DativeBondView`
  constructs it — plus the two molecule-backed constraint-view unit tests. Added `mol.dative_bonds`
  accessor (`molecule.rs`), registered `DativeBondView`/`DativeBondViews` (`lib.rs` + `__init__.py`),
  11 new Rust unit tests + `tests/test_dative.py` (29 pytest). 242 Rust unit + 266 pytest green,
  clippy/fmt clean. **B2 · dative slice COMPLETE.** `[dep: S1a, S1b]`

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
  participant set — how each count pairs with its atom in Python. **Resolved (2026-07-11):
  positional list leaf** — `ElectronCountsAst { Undetermined() \| Lit(list[int]) }` is a 1:1
  mirror of the AST enum (the `ValueAst` leaf pattern), and `electrons` is a `list[int]`
  positionally aligned to `atom_ids` (index i ↔ i-th atom in the participant set). No
  `{atom: count}` mapping. Shared verbatim by B4 (multicenter).

#### B3 · Aromatic — staged impl plan (settled 2026-07-11)

Mirrors the completed dative slice, adapted to aromatic specifics. **All-additive** — the breaking
`from_atoms_and_bonds` → `from_parts` rename (dative's S1b) already landed, so adding an `aromatic`
kwarg is additive; the tree stays green after **every** subitem (no red anywhere). Deltas vs dative:
(1) a **new foundation leaf** `ElectronCountsAst`/`ElectronCountsArg` in its own `electrons.rs`, shared
verbatim by B4; (2) the constraint surface is **simpler** — one key `ElectronCount`, so **no** ring
proxy / `RingSizeCounts` / `ring_count` / `ring_size_count`; (3) participants are a **single unordered
atom set** → the view exposes `atom_ids` only (no acceptor/donor), collection lookups are
`connecting(atoms)` + `incident(atom)`; (4) value fields are `electrons`/`charge`/`spin`/`constraints`
(bond-shaped, `electrons` replaces `order`); (5) `from_parts` gains `aromatic=[]` — **new** (S1b wired
only `dative`). Naming: the value `electrons` (per-atom `ElectronCountsAst` vector) is distinct from the
constraint `electron_count` (total-π `ValueAst`) — both exist, no clash.

**S0 — new leaf `ElectronCountsAst` (`umol-py/src/electrons.rs`).** Foundation; own module because B4
reuses it (no lopsided dep).
- **S0a — DONE** *(additive)* — `ElectronCountsAst` pyclass enum `{ Undetermined(), Lit(Vec<i64>) }`
  (`boolean.rs` leaf pattern): `to_ast`, `__eq__`/`__hash__`/`__repr__`, `as_lit() -> Option<list[int]>`.
  New module `umol-py/src/electrons.rs`; registered in `lib.rs` + `__init__.py`; 5 Rust unit tests
  (roundtrip + as_lit); `import umol` + 266 pytest green. `from_ast` gated `#[cfg(test)]` (test-only
  consumer at S0a, `from_inner` precedent — avoids a dead-code warning); **S1c un-gates it** when the
  value pyclass calls it. `[dep: none]`
- **S0b — DONE** *(additive)* — `ElectronCountsArg` (`FromPyObject` enum `{ Lit(Vec<i64>),
  Ast(Py<ElectronCountsAst>) }`) + `to_ast(py)`, mirroring `ValueArg`; a bare `list[int]` coerces to
  `Lit`. Gated `#[cfg(test)]` (first consumer is the S1c value setter — avoids a dead-code warning)
  with a `to_ast` mapping test (both variants); **S1c un-gates it** with `from_ast`. `[dep: S0a]`

**S1 — value + WET constraint surface (`umol-py/src/aromatic.rs`).**
- **S1a — DONE** *(additive)* — `AromaticSystemConstraintKey { ElectronCount() }` +
  `AromaticSystemConstraintAst { ElectronCount(Py<ValueAst>) }` (new module `umol-py/src/aromatic.rs`):
  `key`/`__eq__`/`__hash__`/`__repr__` + `from_ast`/`to_ast`. One key only — no `RingScope`/
  `RingMembership`; the unit key makes `Key::from_ast` infallible (no `py`/`into_py_variant`) and the
  constraint `key()` getter returns the key directly (not `PyResult`). `AromaticSystemConstraintAst::
  from_ast` + the `into_py_variant` import are `#[cfg(test)]`-gated (no non-test consumer until S1b —
  avoids a dead-code warning); **S1b un-gates them**. Registered both pyclasses (`lib.rs` +
  `__init__.py`); 4 Rust unit test cases; `import umol` + 266 pytest green, clippy/fmt clean.
  Adversarial verification workflow (3 review dims → verify) confirmed 1 finding (test-ordering: key
  test before constraint-AST tests, per the skill's group-by-definition-order rule) — fixed; 3 other
  findings refuted as intended. `[dep: ValueAst leaf]`
- **S1b — DONE** *(additive)* — `AromaticSystemConstraintsAst` container: the uniform mapping API
  (`new`/`__repr__`/`set`/`pop`/`update`/`__len__`/`__iter__`/`keys`/`values`/`items`/`get`/
  `__getitem__`/`__delitem__`/`__contains__`) + `electron_count` getter/setter + `asdict`
  (`{"electron_count"}`) + `AromaticSystemConstraintsUpdate { Container, Entries }` + the 3 iterators +
  `aromatic_system_constraints_asdict`. **No** ring proxy / `RingSizeCounts` / `RingSizeBacking`.
  Adjustments vs the plan: `AromaticSystemConstraintsUpdate::View` and `AromaticSystemConstraintsArg`
  **deferred to S1c** (both reference the S1c view type); the container's `inner_mut` **omitted** (its
  only dative user was the value-backed ring proxy, absent here) — only `inner()` + cfg-test
  `from_inner`; un-gated `AromaticSystemConstraintAst::from_ast` + `into_py_variant` (iterators consume
  them). Registered the container; 14 helper-free container test cases; 32 Rust unit + 266 pytest
  green, clippy/fmt clean. Adversarial verification workflow confirmed 3 test-quality findings
  (summary-stat-only assertions in `new`/`update_entries`; missing `get` non-None-default branch) —
  all fixed; 2 findings refuted (a delegating-`keys()` coverage nit, and getitem/delitem ordering that
  faithfully mirrors dative). `[dep: S1a]`
- **S1c — DONE** *(additive)* — `AromaticSystemAst` value pyclass: `new(electrons: ElectronCountsArg,
  *, charge=None, spin=None, constraints=None)`, `parse`/`__str__`/`__repr__`, getters/setters
  `electrons`(→`ElectronCountsArg`)/`charge`(→`ValueArg`)/`spin`(→`SpinStateAst`)/`constraints`,
  `asdict` (`{electrons, charge, spin, constraints}`), `inner`/`inner_mut`/`from_inner`; plus
  `AromaticSystemConstraintsView` (live handle, full mapping API + `electron_count` getter/setter, no
  ring proxy) + `AromaticSystemConstraintsBacking { AromaticSystem(Py<AromaticSystemAst>) }` —
  **Molecule arm deferred to S3** (its only constructor is `AromaticSystemView`). Added (deferred from
  S1b): `AromaticSystemConstraintsArg { Container, View }` + the `AromaticSystemConstraintsUpdate::View`
  variant/`apply` arm. **Un-gated `ElectronCountsAst::from_ast` + `ElectronCountsArg` (S0)** — the value
  consumes them. Registered value + view (`lib.rs` + `__init__.py`); 13 value/view test cases; 281 Rust
  unit + 266 pytest green, clippy/fmt clean. Adversarial verification workflow confirmed 3 low
  test-quality findings (misnamed `new_constraints` test; `spin` asserted via `inner()` not the getter
  + uncovered `new` spin kwarg; value-test ordering `new`-before-`parse`) — all fixed; 1 refuted.
  `[dep: S0, S1a, S1b]`

**S2 — `from_parts` wiring (`umol-py/src/molecule.rs`).**
- **S2a — DONE** *(additive)* — added `aromatic=Vec::new()` kwarg to `from_parts` (`(atoms, *,
  bonds=[], dative=[], aromatic=[])`), wired to `MoleculeParts.aromatic` (Python entry `(list[int],
  AromaticSystemAst)` → `(Vec<AtomId>, AromaticSystemAst)`). Additive kwarg — **no** Python test-site
  migration. Extended `test_molecule_ast_from_parts` to cover the aromatic wiring (asserts via
  `inner().aromatic_systems()` count + `atom_ids`); Python smoke confirms the kwarg is accepted (no
  `mol.aromatic_systems` accessor until S3). 281 Rust unit + 266 pytest green, clippy/fmt clean.
  Trivial mechanical mirror of dative's S1b wiring — no verification workflow. `[dep: S1c]`

**S3 — views (`umol-py/src/aromatic.rs` + `molecule.rs`).**
- **S3a — DONE** *(additive)* — `AromaticSystemView`: `id`, read-only `atom_ids -> tuple` (member set,
  no acceptor/donor split), settable `electrons`/`charge`/`spin`/`constraints`, `asdict`
  (`{electrons, charge, spin, constraints}`); write-through via `aromatic_system_mut(id)`. **Absorbed
  S3c:** the `Molecule` backing arm's only constructor is `AromaticSystemView::constraints`, so re-adding
  it (+ `read`/`with_mut` Molecule arms) and the molecule-backed constraint-view test landed here, not
  in a separate stage. Not yet registered (S3d); `benzene` test helper + 8 view tests + molecule-backed
  test; 55 aromatic Rust tests green, clippy/fmt clean. Verification workflow confirmed 1 low finding
  (molecule-backed test drove `set_electron_count` but was named `..._set_...`; changed it to drive
  `set`, matching the dative template + its name) — fixed. `[dep: S1c]` **(S3c absorbed here.)**
- **S3b — DONE** *(additive)* — `resolve_aromatic_system_index` (mirrors `resolve_bond_index`) +
  `AromaticSystemViews` (`mol.aromatic_systems`: `__len__`/`__repr__`/`__getitem__`/`__setitem__`
  value-replace/`__iter__` + `connecting(atoms)` single-set-lookup + `incident(atom)`) +
  `AromaticSystemViewIter`. `AromaticSystemViews::new` + the accessor + registration deferred to S3d, so
  the 7 collection tests construct the collection via struct-literal (no dead code). 62 aromatic Rust
  tests green, clippy/fmt clean. Verification workflow confirmed 1 low finding (the `_repr` test was
  ordered last but `__repr__` is the 2nd-declared method — moved it after `_len_and_getitem`) — fixed;
  1 refuted. `[dep: S3a]`
- **~~S3c~~ — absorbed into S3a** (the Molecule backing arm must land with its constructor,
  `AromaticSystemView`). `[dep: S1c, S3a]`
- **S3d — DONE** *(additive)* — `AromaticSystemViews::new` + `mol.aromatic_systems` accessor
  (`molecule.rs`) + registered `AromaticSystemView`/`AromaticSystemViews` (`lib.rs` + `__init__.py`) +
  `tests/test_aromatic.py` (33 pytest); maturin rebuilt. 296 Rust unit + 299 pytest green, clippy/fmt
  clean. Verification workflow confirmed 2 low test-coverage findings (missing `__delitem__`-absent-key
  KeyError test; `test_..._new` dropped the empty-default-constraints assertion) — both fixed; 5
  refuted. `[dep: S3a, S3b]`

**B3 · aromatic slice COMPLETE** — value + WET constraint surface + `AromaticSystemView`/`Views` +
`mol.aromatic_systems` + the new `ElectronCountsAst` leaf (`electrons.rs`, shared by B4); 9 aromatic
pyclasses registered; 296 Rust unit + 299 pytest green.

Critical path: **S0 → S1 → S2 → S3** (linear). No deferrable stages; no red stages (fully additive).
B4 (multicenter) reuses S0's `electrons.rs` verbatim and is otherwise byte-for-byte this plan.

### B4 · Multicenter bond

- **Structurally byte-for-byte identical to aromatic** (`electrons`/`charge`/`spin`/
  `constraints`; `ElectronCount(ValueAst)` constraint). One template serves both; only the
  type names differ. Same `electrons`-pairing design call as B3.

#### B4 · Multicenter — staged impl plan (settled 2026-07-12)

**Design verified byte-for-byte identical to the completed aromatic slice (B3)** against the AST:
`MulticenterBondAst { electrons: ElectronCountsAst, charge: ValueAst, spin: SpinStateAst, constraints }`
(ctors `new(electrons)`/`from_electrons`/`with_charge`/`with_spin`/`with_constraint`, `FromStr`+`Display`);
one constraint key `MulticenterBondConstraintAst::ElectronCount(ValueAst)` / `MulticenterBondConstraintKey::
ElectronCount` with an `electron_count()` container accessor; `MulticenterBondView` is a **single unordered
atom set** (`atom_ids`/`electrons`/`charge`/`spin`/`constraints`, no acceptor/donor); collection
`MulticenterBondViews` has `count`/`ids`/`get`/`contains`/`connecting(atoms)`/`incident(atom)`;
`MoleculeParts.multicenter = Vec<(Vec<AtomId>, MulticenterBondAst)>`; accessors `multicenter_bonds()` /
`multicenter_bond_mut()` / `MulticenterBondId`. So B4 is the **aromatic slice with a rename**:
`Aromatic`→`Multicenter`, `AromaticSystem`→`MulticenterBond`, `mol.aromatic_systems`→`mol.multicenter_bonds`,
`aromatic=[]`→`multicenter=[]`, module `aromatic.rs`→`multicenter.rs`. Concrete per-entity module (WET —
no dedup with aromatic, mirroring the separate AST types). Fully additive; no red stages.

**Three improvements over B3's staging** (learned while building B3):
1. **No S0.** The `ElectronCountsAst`/`ElectronCountsArg` leaf already exists (`electrons.rs`, built +
   un-gated in B3) — B4 reuses it directly, no gating.
2. **Merge constraint key+enum+container into one stage (S1a)** so `MulticenterBondConstraintAst::from_ast`
   has its container-iterator consumer immediately — **no `#[cfg(test)]` gating** (B3 split these across
   S1a/S1b and had to gate).
3. **The `Molecule` backing arm lands in S3a with the view** (its only constructor is
   `MulticenterBondView::constraints`) — not a separate stage (B3's S3c was absorbed into S3a anyway).

Naming note (as in B3): the value `electrons` (per-atom `ElectronCountsAst` vector) is distinct from the
constraint `electron_count` (total `ValueAst`) — both exist, no clash.

**S1 — value + WET constraint surface (`umol-py/src/multicenter.rs`).**
- **S1a — DONE** *(additive)* — `MulticenterBondConstraintKey { ElectronCount() }` +
  `MulticenterBondConstraintAst { ElectronCount(Py<ValueAst>) }` (unit key ⇒ infallible `Key::from_ast`,
  `key()` getter returns the key directly) + `MulticenterBondConstraintsAst` container (uniform mapping API:
  `new`/`__repr__`/`set`/`pop`/`update`/`__len__`/`__iter__`/`keys`/`values`/`items`/`get`/`__getitem__`/
  `__delitem__`/`__contains__` + `electron_count` getter/setter + `asdict` `{"electron_count"}`) +
  `MulticenterBondConstraintsUpdate { Container, Entries }` + the 3 iterators +
  `multicenter_bond_constraints_asdict`. **No** ring proxy; **no** container `inner_mut` (only `inner()` +
  cfg-test `from_inner`). `from_ast` **not** gated (container consumes it — the B4 improvement). Verbatim
  rename of aromatic's constraint surface (its test-quality fixes carried over). Registered the 3 constraint
  pyclasses (`lib.rs` + `__init__.py`); 17 constraint tests (key/enum/container); 21 Rust unit + 299 pytest
  green, clippy/fmt clean. Verification: the review workflow failed as **infra** (both runs — agents did the
  review then returned `null` on structured-output emission); substituted a manual rename-slip scan (zero
  `Aromatic` leaks) + method-inventory/test-parity diff vs aromatic (identical modulo the 5 deferred view
  tests). `[dep: existing ValueAst]`
- **S1b — DONE** *(additive)* — `MulticenterBondAst` value pyclass (`new(electrons: ElectronCountsArg, *,
  charge=None, spin=None, constraints=None)`, `parse`/`__str__`/`__repr__`, getters/setters
  `electrons`(→`ElectronCountsArg`)/`charge`(→`ValueArg`)/`spin`(→`SpinStateAst`)/`constraints`, `asdict`
  `{electrons, charge, spin, constraints}`, `inner`/`inner_mut`/`from_inner`) + `MulticenterBondConstraintsView`
  (live handle, full mapping API + `electron_count` getter/setter, no ring proxy) +
  `MulticenterBondConstraintsBacking { MulticenterBond(Py<MulticenterBondAst>) }` — **Molecule arm deferred
  to S3a** — + `MulticenterBondConstraintsArg { Container, View }` + the `MulticenterBondConstraintsUpdate::View`
  variant/`apply` arm. Reuses the existing `ElectronCountsAst`/`ElectronCountsArg` (no gating). Verbatim
  rename of aromatic's S1c. Registered value + view (`lib.rs` + `__init__.py`); 9 value + 4 view tests;
  36 Rust unit + 299 pytest green, clippy/fmt clean. Verification: the review workflow's structured-output
  emission reliably dies on the rename-diff prompt shape (failed on S1a twice), so substituted the manual
  rename-slip scan (zero leaks) + value/view method-inventory + test-parity diff vs aromatic (identical
  modulo the deferred S3 view/collection + the one molecule-backed view test → S3a). `[dep: S1a]`

**S2 — `from_parts` wiring (`umol-py/src/molecule.rs`).**
- **S2a — DONE** *(additive)* — added `multicenter=Vec::new()` kwarg to `from_parts` (`(atoms, *, bonds=[],
  dative=[], aromatic=[], multicenter=[])`), wired to `MoleculeParts.multicenter` (Python entry `(list[int],
  MulticenterBondAst)` → `(Vec<AtomId>, MulticenterBondAst)`). Additive — **no** test-site migration.
  Extended `test_molecule_ast_from_parts` (asserts via `inner().multicenter_bonds()` count + `atom_ids`);
  Python smoke confirms the kwarg is accepted (no `mol.multicenter_bonds` accessor until S3c). 13 molecule
  Rust tests + 299 pytest green, clippy/fmt clean. Trivial mechanical mirror — no workflow. `[dep: S1b]`

**S3 — views (`umol-py/src/multicenter.rs` + `molecule.rs`).**
- **S3a — DONE** *(additive)* — `MulticenterBondView` (`id`, read-only `atom_ids -> tuple`, settable
  `electrons`/`charge`/`spin`/`constraints`, `asdict`; write-through via `multicenter_bond_mut(id)`) +
  **re-added the `Molecule` backing arm** to `MulticenterBondConstraintsBacking` (+ `read`/`with_mut`
  Molecule arms) — `MulticenterBondView::constraints` constructs it — + the molecule-backed constraint-view
  unit test (drives `set`, matching the aromatic S3a fix). `three_center_bond` fixture (3 borons + one
  3-center bond). Not registered yet (S3c). 44 multicenter Rust tests green, clippy/fmt clean.
  Verification: manual parity — 0 `Aromatic` leaks; `MulticenterBondView` pymethod inventory **identical**
  to `AromaticSystemView` (all 12); 7=7 view tests + the molecule-backed test present. `[dep: S1b]`
- **S3b — DONE** *(additive)* — `resolve_multicenter_bond_index` (mirrors `resolve_bond_index`) +
  `MulticenterBondViews` (`mol.multicenter_bonds`: `__len__`/`__repr__`/`__getitem__`/`__setitem__`
  value-replace/`__iter__` + `connecting(atoms)` + `incident(atom)`) + `MulticenterBondViewIter`. `new` +
  accessor + registration deferred to S3c → the 7 collection tests construct via struct-literal (no dead
  code); `_repr` test placed 2nd (per the aromatic S3b fix). 51 multicenter Rust tests green, clippy/fmt
  clean. Verification: manual parity — 0 leaks; `MulticenterBondViews` pymethod inventory **identical** to
  `AromaticSystemViews` (all 7); 7=7 collection tests. `[dep: S3a]`
- **S3c — DONE** *(additive)* — `MulticenterBondViews::new` + `mol.multicenter_bonds` accessor
  (`molecule.rs`) + registered `MulticenterBondView`/`MulticenterBondViews` (`lib.rs` + `__init__.py`) +
  `tests/test_multicenter.py` (31 pytest, renamed from `test_aromatic.py`); maturin rebuilt. 344 Rust
  unit + 332 pytest green, clippy/fmt clean. Verification: manual parity — 0 leaks; 31=31 Python tests
  vs `test_aromatic.py`. `[dep: S3a, S3b]`

**B4 · multicenter slice COMPLETE** — value + WET constraint surface + `MulticenterBondView`/`Views` +
`mol.multicenter_bonds`, reusing the S0 `ElectronCountsAst` leaf (no new foundation); 9 multicenter
pyclasses registered; 344 Rust unit + 332 pytest green. Verbatim rename of the aromatic slice with the
staging improvements (no S0, merged constraint surface = no gating, backing-with-view). The review
workflow's structured-output emission reliably dies on the rename-diff prompt shape (S1a twice), so the
multicenter stages were verified by clean build/clippy/fmt + full test suites + a manual rename-slip +
method/test-parity diff vs aromatic (semantic correctness inherited from the adversarially-verified
aromatic slice).

Critical path: **S1 → S2 → S3** (linear). No S0, no gating, no deferrable/red stages. After B4, B5–B7
(noncovalent, stereo) remain; the `electrons.rs` leaf is now shared by aromatic + multicenter.

### B5 · Noncovalent bond

- **PAUSED (2026-07-12) pending an AST change — see doc 117 §4.** Scoping B5 surfaced that the noncovalent
  constraint enum is uninhabited, and pyo3 rejects a zero-variant `#[pyclass] enum`, so the constraint
  surface can't be mirrored. Resolution (settled): don't work around the emptiness — inhabit the
  constraint upstream with `NoncovalentBondConstraintAst::Intramolecular(BooleanAst)` (`#I`), threading it
  through the whole molecule-AST machinery (doc 117 §4 has the scope + staged plan A–C). B5 then resumes as
  stage D of that plan: a **standard 1-key Boolean-constraint slice** (constraint half = bond-`Aromatic`
  shape with an `intramolecular` getter/setter; view half = 2-atom bond-shaped) + the `NoncovalentBondKind`
  leaf. The two design calls below are thereby **resolved**: constraints are the real 1-key surface (not a
  stub/omit); `kind` is the ElectronAst-style leaf (settled 2026-07-12).
- Value `NoncovalentBondAst { kind: NoncovalentBondKindAst, constraints }`.
- **New leaves: `NoncovalentBondKindAst { Undetermined | Lit(NoncovalentBondKind) }` and
  `NoncovalentBondKind { HydrogenBond, HalogenBond, ChalcogenBond, Ionic, VanDerWaals }`**
  (a fieldless value enum → hashable, per finding 137-p3-2). Python: `NoncovalentBondKind` pyclass enum +
  `NoncovalentBondKindAst { Undetermined() \| Lit(NoncovalentBondKind) }` + `NoncovalentBondKindArg`,
  mirroring `ElementAst`/`Element` (settled 2026-07-12).
- Constraints: uninhabited today; **being inhabited** with `Intramolecular(BooleanAst)` (doc 117 §4), after
  which the Python surface is a normal 1-key Boolean constraint (bind the full mapping API — no stub/omit).
- View: read-only atom pair (fixed 2, unordered); settable `kind`. Collection
  `mol.noncovalent_bonds`: id-indexed + `connecting(a, b)`.
- Design calls: `kind` wildcard (`Undetermined`) → an optional/nullable surface + coercion;
  the uninhabited-constraints decision above.

#### B5 — staged impl plan (resumed 2026-07-12, after doc 117 §4 A–C)

Composite of three peers: constraint half ≈ **aromatic** (single-key container, single-value
`intramolecular()` accessor; Boolean value reuses the already-bound `BooleanAst`/`BooleanArg`); view/
collection half ≈ **bond** (fixed unordered 2-atom pair — `atom_ids` 2-tuple, `connecting(a, b)`,
`incident(atom)`); kind leaf ≈ **`ElementAst`/`Element`**. The value pyclass is the thinnest of any slice
(`kind` + `constraints`, no order/charge/spin). New file `umol-py/src/noncovalent.rs`; the whole slice is
**additive** (new file + additive `molecule.rs` kwarg/accessor) — no red→green. Dead-code discipline
carries over from B2–B4: an unused `pub(crate)` constructor gets `#[cfg(test)]` until its non-test
consumer lands; a backing-enum `Molecule` arm lands with the view that constructs it.

- **S0 — kind leaf.**
  - **S0a** `noncovalent.rs`: `NoncovalentBondKind` (fieldless value pyenum, 5 variants, hashable) +
    `NoncovalentBondKindAst { Undetermined() | Lit(NoncovalentBondKind) }` (`as_lit`/`__eq__`/`__hash__`/
    `__repr__` + `from_ast`/`to_ast`) + `NoncovalentBondKindArg` (bare-enum | `Ast` passthrough); register
    both pyclasses + export; leaf tests. **Additive.** `[dep: —]`
- **S1 — constraint half** (independent of S0). Merged into one subitem per the B3/B4 lesson (avoids
  `#[cfg(test)]`-gating `from_ast`; the container's `from_ast` has no non-test consumer until S3's view,
  but key+element+container ship together so the mapping API exercises it).
  - **S1a** `noncovalent.rs`: `NoncovalentBondConstraintKey` (`Intramolecular`) + `NoncovalentBondConstraintAst`
    (`Intramolecular(BooleanAst)`, `key`/`__eq__`/`__hash__`/`__repr__`) + `NoncovalentBondConstraintsAst`
    container (uniform mapping API `__len__`/`__getitem__`/`__setitem__`/`__delitem__`/`__contains__`/
    `__iter__` + `intramolecular` getter/setter + `set`/`asdict`) + `NoncovalentBondConstraintsArg`/
    `…Update` + the three iters (`…ConstraintIter`/`…KeyIter`/`…ItemsIter`); register + export; container
    tests. **Additive.** `[dep: —]` (Boolean leaf already bound)
  - **S1b** `noncovalent.rs`: `NoncovalentBondConstraintsView` + `NoncovalentBondConstraintsBacking` with
    the **own-value arm only** (`Noncovalent(Py<NoncovalentBondConstraintsAst>)`) + `intramolecular`
    getter/setter + mapping API mirroring the container; register + export; view tests. **Additive.**
    `[dep: S1a]`
- **S2 — value pyclass.**
  - **S2a** `noncovalent.rs`: `NoncovalentBondAst` — `new(kind, *, constraints=[])`, `parse`/`__str__`/
    `__repr__` (DSL; `#I` now parses, per doc 117 §4 B), `kind`/`constraints` getters+setters (via
    `NoncovalentBondKindArg` / `NoncovalentBondConstraintsArg`), `asdict`, `inner`/`from_inner`; **ungate**
    the S0a/S1a constructors it consumes; register + export; value tests incl. a `#I` parse roundtrip.
    **Additive.** `[dep: S0a, S1a, S1b]`
- **S3 — view + collection + molecule wiring.**
  - **S3a** `noncovalent.rs`: `NoncovalentBondView` (read-only `id`/`atom_ids` 2-tuple, `kind`/`constraints`
    get+set, `asdict`, `__repr__`) + `resolve_noncovalent_bond_index` + **add** the
    `NoncovalentBondConstraintsBacking::Molecule` arm here (with the view that constructs it); tests.
    **Additive.** `[dep: S2a, S1b]`
  - **S3b** `noncovalent.rs`: `NoncovalentBondViews` collection (`__len__`/`__getitem__`/`__setitem__`/
    `__iter__` + `connecting(a, b)` + `incident(atom)`) + `NoncovalentBondViewIter`; register + export;
    collection tests. **Additive.** `[dep: S3a]`
  - **S3c** `molecule.rs`, `lib.rs`, `__init__.py`, `tests/test_noncovalent.py`: `mol.noncovalent_bonds`
    accessor + `from_parts(noncovalent=[([a, b], NoncovalentBondAst), …])` kwarg; extend
    `test_molecule_ast_from_parts`; final registration/exports check. **Additive.** `[dep: S3b]`

Critical path: **(S0a ∥ S1a→S1b) → S2a → S3a → S3b → S3c.** No deferrable stages — the slice is small and
every subitem is on the path to the `mol.noncovalent_bonds` surface.

#### B5 progress (2026-07-12)

- **S0a done** — kind leaf (`NoncovalentBondKind` + `NoncovalentBondKindAst` + `…Arg`). Registered/exported.
- **S1a done** — constraint half (key + `Intramolecular(BooleanAst)` element + container + `…Update{Container,Entries}` + 3 iters). `…Update::apply` returns `PyResult` so the view arm lands additively.
- **S1b folded into S2 (staging correction).** The plan had S1b as a standalone container-backed view `[dep: S1a]` preceding S2a. That is a **dependency error**: the peer's constraints view is backed on the *entity* (`Backing::AromaticSystem(Py<AromaticSystemAst>)`), and the entity's `constraints` getter *returns* that view — the two are mutually referential, so the view cannot precede the value pyclass. A container-backed view (as S1b drafted it) would also be redundant, since the S1a container is already a mutable pyclass. Correction: the entity-backed `NoncovalentBondConstraintsView` + backing (own-value arm; Molecule arm still deferred to S3a) + the `View` variant on `…Update`/`…Arg` all ship **with** the value pyclass in S2. So real dependency is **S2 `[dep: S0a, S1a]`**, view included; the plan's `S2 [dep: …S1b]` inverted it.
- **S2 done** — `NoncovalentBondAst` value pyclass (`new`/`parse`/`__str__`/`__repr__`/`kind`/`constraints`/`set_constraints`/`asdict`/`inner`/`inner_mut`/`from_inner`) + entity-backed constraints view (full mapping API + write-through) + `…Arg`; ungated the S0a kind `from_ast`/`Arg`. 42 Rust unit tests, clippy/fmt clean, Python write-through verified (`b.constraints.intramolecular = True` re-reads through the bond; `Hbd#I` round-trips).
- **Lesson for the remaining slices (stereo B6/B7):** an entity-backed constraints view and its owning value pyclass are one atomic unit — never stage the view before the value pyclass. Fix the general B-order note if reused.
- **Adversarial review of S2 (4-lens workflow) found 2 real panics the green build missed** (parity + fidelity lenses were clean — the mirror is faithful). Both are self-aliasing `RefCell` double-borrow panics: `bond.constraints = bond.constraints` (the `set_constraints` setter held `&mut self` across `value.to_ast()`, whose `View` arm re-borrows the same bond) and `bond.constraints.update(bond.constraints)` (the view `update` held `with_mut` across `apply()`, whose `View` arm re-reads the bond). **Fixed in noncovalent** by resolving every Python read to owned data *before* the write borrow: `NoncovalentBondConstraintsUpdate::resolve` → `ResolvedNoncovalentBondConstraintsUpdate::apply`, and `set_constraints`/container `update` take `slf: Py<Self>` and snapshot before `borrow_mut`. Regression tests added for all four self-alias paths (view/container `update`-self, `set_constraints` self + from-view). 46 Rust tests, clippy/fmt clean, Python repro confirms no panic.
- **SYSTEMATIC peer sweep done (2026-07-12).** The same self-alias panic existed in every peer slice; fixed all of them with the resolve-before-borrow transform: `<E>ConstraintsUpdate::apply` → `resolve` + `Resolved<E>ConstraintsUpdate::apply`; container `update` and value `set_constraints` take `slf: Py<Self>` and read before `borrow_mut`; the constraints-view `update` resolves before `with_mut`. Files: `atom.rs` + **`constraint.rs`** (atom's constraint container/view live in the shared `constraint.rs`, not `atom.rs` — a grep-by-`atom.rs` miss the end-to-end Python repro caught), `bond.rs`, `dative.rs`, `aromatic.rs`, `multicenter.rs`. Regression tests added per entity (container update-self, value set-constraints-self/from-view, view update-self). **Already-safe, left alone:** the molecule-view `set_constraints(&self)` (single `owner.borrow_mut()…= value.to_ast()?` — Rust evaluates the RHS read before the LHS place borrow); ring-size sub-containers (setters take `ValueArg`, cannot alias). **Verified:** 405 Rust unit tests, 332 pytest, clippy/fmt clean, and a Python repro proving `x.constraints = x.constraints` / `x.constraints.update(x.constraints)` / `cs.update(cs)` are no-ops (not panics) for all 6 entities.
- **Lesson (applies to B6/B7 stereo):** any pyo3 method holding a write borrow (via `&mut self`, or `with_mut`) across a call that re-reads a possibly-aliased entity (`value.to_ast()`/`other.apply()` with a `View` arm) self-alias-panics; the fix is always resolve-all-reads-to-owned-data first. Write the stereo constraint surface in the resolved-before-borrow form from the start.
- **S3a/b/c done (2026-07-12) — B5 COMPLETE.** S3a: molecule-embedded `NoncovalentBondView` (`id`, 2-tuple `atom_ids`, `kind` get/set, `constraints` get/set, `asdict`) + the `Backing::Molecule` arm on the constraints view. S3b: `NoncovalentBondViews` collection (`__len__`/`__getitem__`/`__setitem__`/`__iter__` + `connecting(a, b)` fixed-pair + `incident(atom)`) + `resolve_noncovalent_bond_index` + `NoncovalentBondViewIter`. S3c: `mol.noncovalent_bonds` accessor + `NoncovalentBondViews::new` + `from_parts(noncovalent=[([a, b], NoncovalentBondAst), …])` kwarg + `tests/test_noncovalent.py`. The view's `set_constraints` uses the safe RHS-first single-assignment form; a molecule-backed self-alias regression is included. **Verified: 417 Rust unit tests, 363 pytest, clippy/fmt clean.**
- **Slice tally:** new file `noncovalent.rs` binds the full surface — kind leaf (`NoncovalentBondKind`/`NoncovalentBondKindAst`), 1-key Boolean constraint (`Intramolecular`) container + live view, value pyclass (`NoncovalentBondAst`), molecule view + collection, molecule wiring. Pyclasses registered/exported; 31 pytest + ~58 Rust unit tests.

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
- **Write the constraint surface in the resolve-before-borrow form from the start** (self-alias panic
  lesson from the B5 sweep — see the B5 progress note). The stereo constraint bag has four keys
  (`LigandSymmetryAst`/`FluxionalityAst`/`TopicityAst`/`StereogenicityAst`), so the same three methods
  recur: the value/overlay `set_constraints`, the container `update`, and the constraints-view `update`.
  In every one, resolve each `…ConstraintsArg`/`…ConstraintsUpdate` (its `Container`/`View` arms read a
  possibly-aliased entity via `read`) to owned data **before** taking any write borrow — i.e.
  `set_constraints`/container `update` take `slf: Py<Self>` and snapshot before `borrow_mut`; the view
  `update` calls `other.resolve(py)?` before `with_mut`. Do NOT write the naive `&mut self` +
  `value.to_ast(py)?` / `with_mut(|cs| other.apply(...))` form — it panics on `x.constraints =
  x.constraints`. Add the four self-alias regression tests per entity. The molecule-view
  `set_constraints(&self)` single-assignment form (`owner.borrow_mut()…= value.to_ast()?`) is safe as-is
  (RHS evaluated before the LHS place borrow); mirror it, don't "fix" it.

### B6/B7 — resolved design (2026-07-12)

The open calls are now decided (supersedes the "all stereo rows remain open" note in §C):

1. **Constraint sub-surface — bind all four values now** (no deferral): `LigandSymmetryAst`, `FluxionalityAst`, `TopicityAst`, `StereogenicityAst`, on the helper leaves `LigandPermutation`, `OrientedLigandPermutation`, `StereoLigandPair`, `TopicityRelationAst`, plus the derived `Topicity`/`Stereogenicity` enums. (Deferring all four would recreate the empty-pyclass-enum blocker; binding all avoids it.)
2. **No reframe API** — read/write the stored configuration coset as-is; `permutation_for`/`coset_for`/`transform_frame` are NOT surfaced (interning/114-adjacent, later if ever).
3. **Overlay addressing** — id-indexed like the covalent slices, plus content lookups `coincident(site)` (the ≤1 relation at a site) and `connecting(site, ligands)` (ordered ligand frame supplied). The AST views already expose `coincident_id`/`connecting_id`/`ids`/`get`/`count`/`contains`.
4. **Ligands surface as `StereoLigand{atom_id, kind}`** (mirror). A virtual ligand (`ImplicitHydrogen`/`LonePair`) reads out as `{atom_id = bearing atom, kind = …}` — the pair disambiguates. `atom_id` is only non-trivial for stereo *bonds* (two-atom site), but the same `StereoLigand` shape is used for both entities — keep the parallel, don't special-case.
5. **Configuration setter shorthands** — assignment accepts the already-bound `TetrahedralStereo.Ccw/Cw` → `Kinded(Tetrahedral, Lit)` and `CisTransStereo.Z/E` → `Kinded(CisTrans, Lit)`, plus a full `StereoConfigurationAst.Kinded(kind, coset)`. Axial/SquarePlanar/TrigonalBipyramidal/Octahedral have no shorthand (none exists); they use the full form.
6. **B6 + B7 together**, macro-driven. The two entities (site = atom `NodeId` vs bond `EdgeId`) are fully parallel — the AST already generates the constraint types via `relation_ast!`; the Python constraint surface + value pyclass + view + collection are generated per entity by a Rust macro. Leaves are shared singletons (bound once).

### B6/B7 — staged impl plan

New file(s): the stereo leaves extend `umol-py/src/stereo.rs` (which already holds `StereoCosetAst`/`StereoTerm`/`Permutation`/`TetrahedralStereoAst`/`CisTransStereoAst` + the `TetrahedralStereo`/`CisTransStereo` shorthands); the per-entity constraint/value/view surface is macro-generated (in `stereo.rs` or a sibling). Every stage is **additive** — no red→green. `slf: Py<Self>`/resolve-before-borrow throughout the constraint surface (per the note above).

- **S0 — stereo leaves (foundation).**
  - **S0a — DONE** *(additive)* — `StereoKind` (6-variant) + `StereoLigandKind` (3) + `Topicity` (3) + `Stereogenicity` (3), fieldless hashable pyenums (`NoncovalentBondKind` pattern), `from_ast`/`to_ast` `#[cfg(test)]`-gated. Registered + exported + roundtrip tests. `[dep: —]`
  - **S0b — DONE** *(additive)* — `StereoLigand` pyclass `{atom_id: u32, kind}` (`#[pyclass(eq, hash, frozen, from_py_object)]`, getters, eval-able `__repr__`, `new`; `from_ast`/`to_ast` gated). `[dep: S0a]`
  - **S0c — DONE** *(additive)* — `StereoConfigurationAst` (`Undetermined() | Kinded(StereoKind, Py<StereoCosetAst>)`, `kind`/`coset` getters, eq/hash/repr via `to_ast`) + `StereoConfigurationArg` (Th/Ct shorthand → `Kinded` + `StereoConfigurationAst` passthrough). Ungated `StereoKind::to_ast` (live). `StereoConfiguration` ground **deferred** (no read/write consumer). `[dep: S0a]`
- **S1 — constraint value leaves.** `[dep: S0]`
  - **S1a — DONE** *(additive)* — helper leaves `Orientation` (Proper/Improper pyenum mirroring `umol_perm::Orientation`, needed by `OrientedLigandPermutation` — not called out in the original list), `LigandPermutation` (wraps the bound `Permutation`, structural eq/hash), `OrientedLigandPermutation` `{permutation, orientation}`, `StereoLigandPair` `{first, second}` u32 (`new` delegates to `AstStereoLigandPair::new` for identical normalization; `from_ast` live), `TopicityRelationAst` (`Undetermined|Lit|LitSet|NotSet` pyclass enum, `ValueAst` pattern, `as_lit`). Added `PartialOrd, Ord` to the pyo3 `Topicity` derive for `BTreeSet<Topicity>`; ungated `Topicity::to_ast` (live consumer `TopicityRelationAst::to_ast`). 75 Rust stereo tests + Python check green, clippy/fmt clean. `[dep: S0b]`
  - **S1b** the four constraint values `LigandSymmetryAst`, `FluxionalityAst`, `TopicityAst`, `StereogenicityAst` (+ Args) — `Topicity`/`Stereogenicity` values wrap the derived enums; `LigandSymmetry`/`Fluxionality` ride the permutation helpers. `[dep: S1a]`
- **S2 — constraint container + view (macro, per entity).** `[dep: S1]`
  - **S2a** Rust macro `stereo_constraint_surface!{ StereoAtom | StereoBond }` → `<E>ConstraintKey` (4 keys; ligand-pair sub-key where the value is pair-scoped — confirm at build) + `<E>ConstraintAst` (4-variant) + `<E>ConstraintsAst` container (mapping API + per-key accessors) + `<E>ConstraintsArg`/`…Update`/`Resolved…` (resolve-before-borrow) + iters; instantiate for both entities. `[dep: S1b]`
  - **S2b** macro → `<E>ConstraintsView` + backing with the **own-value arm only** (Molecule arm deferred to S4a). `[dep: S2a]`
- **S3 — value pyclass (macro, per entity).** `[dep: S0c, S2]`
  - **S3a** macro `stereo_value!{ … }` → `StereoAtomAst`/`StereoBondAst` `{configuration, constraints}`: `new`, `configuration` get/set (via `StereoConfigurationArg`), `constraints` get (→ view) / set, `asdict`, `inner`/`inner_mut`/`from_inner`. (`parse`/`__str__` only if the AST entity has `FromStr`/`Display` — the overlay likely has no standalone string form; confirm, else omit.) `[dep: S0c, S2a, S2b]`
- **S4 — view + collection (macro, per entity).** `[dep: S3]`
  - **S4a** macro → `StereoAtomView`/`StereoBondView`: read-only `id`, `site_id` (atom/bond u32), `ligands` (list of `StereoLigand`, from the AST `ligands()`→`StereoLigandView`), `kind`, `coset` (from configuration), settable `configuration`, `constraints` get/set (RHS-first single-assignment `set_constraints`), `asdict`; **add** the `Backing::Molecule` arm to the constraints view (with the view that constructs it). `[dep: S3a, S2b]`
  - **S4b** macro → `StereoAtomViews`/`StereoBondViews`: `__len__`/`__getitem__`/`__setitem__`/`__iter__` + `coincident(site)` + `connecting(site, ligands)` + `resolve_stereo_{atom,bond}_index` + iter. `[dep: S4a]`
- **S5 — molecule wiring.** `[dep: S4]`
  - **S5a** `mol.stereo_atoms`/`mol.stereo_bonds` accessors + `<E>Views::new` + `from_parts(stereo_atoms=[(site, [ligands], StereoAtomAst)], stereo_bonds=[…])` kwargs; extend `MoleculeAst.__repr__` (add `stereo_atom`/`stereo_bond` to the non-zero-extras loop); register/export; `tests/test_stereo.py`. `[dep: S4b]`

Critical path: **S0 → S1 → S2 → S3 → S4 → S5** (linear; the macro lands stereo-atom + stereo-bond together at each of S2–S5). No deferrable stages given decision 1 (all four constraint values) and decision 6 (both entities). Build-time confirmations (non-blocking): the constraint key structure (pair sub-keys), whether the value pyclass has a string DSL, and the `StereoLigandView`→`StereoLigand` mapping.

## C. Leaf types needing design calls

| Leaf | Used by | Rust shape | The call |
|---|---|---|---|
| `BooleanAst` | bond, dative (`Aromatic`) | `Undetermined \| Lit(bool)` | **Resolved (2026-07-09):** full `BooleanAst` mirror (`Undetermined`/`Lit(bool)`) is the representation; `True`/`False` coerces on assignment (a `BooleanArg`: `bool → Lit`, mirror passthrough). |
| `CisTransStereoAst` | bond (`CisTransStereo`) | `Undetermined \| NotStereo \| Stereo(StereoCosetAst)` | **Resolved:** full mirror (like `TetrahedralStereoAst`) **plus** a simple `CisTransStereo` enum `.Z`/`.E` for assignment, mirroring `TetrahedralStereo.Ccw/Cw`. `Z → Lit(0)`, `E → Lit(1)` (dsl/stereo.rs:336ff). |
| `ElectronCountsAst` | aromatic, multicenter | `Undetermined \| Lit(list[int])` | **Resolved:** mirror whose `Lit` variant is a Python `list[int]` (positional, aligned to participant order); no `{atom: count}` dict. |
| `NoncovalentBondKindAst` + `NoncovalentBondKind` | noncovalent | mirror `Undetermined \| Lit(NoncovalentBondKind)`; fieldless enum `{HydrogenBond, HalogenBond, ChalcogenBond, Ionic, VanDerWaals}` | **Resolved (2026-07-09):** names follow Rust **exactly** (`NoncovalentBondKind.HydrogenBond`, no shortening); a simple hashable value enum (137-p3-2 pattern). `kind` is the full `NoncovalentBondKindAst` mirror (`Undetermined` as a variant, not `None` — the Rust structure, matching `atom.element: ElementAst`), assigned via the bare enum. |
| noncovalent `constraints` | noncovalent | `Intramolecular(BooleanAst)`, 1 key | **Superseded (2026-07-12):** the earlier "empty stub" is obsolete — the constraint was inhabited upstream (doc 117 §4). Bind the **real** 1-key Boolean surface: `NoncovalentBondConstraintsAst` with an `intramolecular` getter/setter (aromatic-`ElectronCount` container shape, Boolean value via the bound `BooleanAst`). |
| `StereoConfigurationAst` / `StereoKind` | stereo | `Undetermined \| Kinded(kind, coset)`; 6-kind enum | The coset carries no kind — configuration (kind+coset) must be the surfaced unit; any coset op threads `StereoKind`. |
| `StereoLigand` / `StereoLigandKind` | stereo | `{atom_id, kind}`; `{Atom, ImplicitHydrogen, LonePair}` | Surface `StereoLigand{atom_id, kind}` vs. bare atom ids; how virtual ligands (bearing-atom id) appear. |
| `StereoLigandPosition` (vs `ParticipantPosition`) | stereo | `u32` frame index | Which position newtype (if any) surfaces to Python. |
| Stereo constraint values (`LigandSymmetryAst`/`FluxionalityAst`/`TopicityAst`/`StereogenicityAst` + helpers) | stereo | per-site constraint bag, 4 keys | A whole sub-surface — bind now vs. defer within the stereo slice. |
| Stereo overlay addressing | stereo | id vs `(site, ligand-multiset)` content lookup (order-independent) | How a relation is addressed from Python; whether ligand order is supplied for reframe (`permutation_for`/`coset_for`/`transform_frame`); note this is handle-identity-adjacent (interning, 114). |

The covalent leaf calls are resolved (2026-07-09); **the stereo rows are now resolved
(2026-07-12) inline** — see *B6/B7 — resolved design* and *B6/B7 — staged impl plan* above (constraint
sub-surface, overlay addressing, ligand/config representation, and B6+B7 scope all decided).

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
