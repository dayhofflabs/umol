# 222 — Depicting a supplied layout

Status: Proposed
Date: 2026-09-08
Revised: 2026-09-10
Relates: [221](221-depiction-api-2026-09-03.md),
[225](225-depiction-problems-2026-09-09.md),
[Python API guide](../docs/development/python-api.md),
[data-type guide](../docs/development/data-types.md)

## Purpose

Doc 221 sealed the depiction API around generated layouts: `Depict::depict_with` lays a molecule
out with the configured algorithm and lowers that layout in one step, and the lowering operation
that accepts an independently supplied `MoleculeLayout` became crate-private. It kept
`MoleculeLayout` public as a coordinate representation that can be generated and edited, but with
no public way to draw the edited result.

A structure editor needs exactly that missing step. The user arranges a drawing — moves an atom,
rotates a fragment, keeps a stored arrangement across a redraw — and the redraw must use the
arranged coordinates rather than regenerate them. Regeneration discards the arrangement. This
document reopens the one item 221 closed and binds the layout representation in Python.

The first revision of this document settled a public free function `depict_molecule`, a frozen
Python `MoleculeLayout` with a functional `with_position`, and left reactions CoordGen-only. Review
of the implementing pull request asked for a different shape on five points and for the editor's
requirements to be stated against it before another implementation round. This revision is that
design. The implementation on the pull request's branch is the reviewed shape and will be replaced
once the surface below is settled; its frame-check and round-trip tests carry over.

## Reversed portion of doc 221

Doc 221 settled that `MoleculeLayout`'s "public existence does not require a public operation for
lowering an arbitrary supplied layout into a depiction", removed the layout-frame mismatch from
`MoleculeDepictionError` "because the public depiction operation no longer accepts an independently
supplied layout", and stated that the layout API "is not an alternate public constructor for
`Depiction`". Those three statements are reversed here. Everything else in 221 stands: `Depict`,
`DepictConfig`, the opaque `Depiction`, `render_svg`, the `depiction` feature boundary, and the
reaction path are unchanged in kind, though the trait and the reaction path gain operations below.
The public free function `layout_molecule` predates 221 and is retired here (§ 3).

## Proposed surface

Line references are to `main` at `98934189b`.

### 1. `depict_layout` on the `Depict` trait

Today `pub trait Depict { type Error; fn depict(&self); fn depict_with(&self, &DepictConfig); }`
(`umol-io/src/depict.rs:45`–`66`). `Molecule`'s lowering is the crate-private
`fn depict(&Molecule, &MoleculeLayout)` (`umol-io/src/depict/molecule.rs:53`), which reads positions
with an `expect` justified by frame agreement (`:398`). `Reaction::depict_with`
(`umol-io/src/depict/reaction.rs:84`–`115`) materializes the span, calls `layout_molecule` per side
(`:94`, `:97`), lowers each side with `molecule::depict`, then `compose_sides` (`:29`) places the
sides and the arrow.

Proposed:

```rust
pub trait Depict {
    type Layout;
    type Error;

    fn layout(&self) -> Result<Self::Layout, Self::Error> {
        self.layout_with(&DepictConfig::default())
    }
    fn layout_with(&self, config: &DepictConfig) -> Result<Self::Layout, Self::Error>;

    fn depict(&self) -> Result<Depiction, Self::Error> {
        self.depict_with(&DepictConfig::default())
    }
    fn depict_with(&self, config: &DepictConfig) -> Result<Depiction, Self::Error> {
        self.depict_layout(&self.layout_with(config)?)
    }
    fn depict_layout(&self, layout: &Self::Layout) -> Result<Depiction, Self::Error>;
}
```

`Molecule::Layout = MoleculeLayout`, `Reaction::Layout = ReactionLayout` (§ 2). `depict_with` has a
default body, so generated-layout and supplied-layout depiction converge on `depict_layout` by
construction: there is one lowering per type and the supplied path is the only path. Neither
implementor overrides `depict_with`. No public free function is added; the first revision's
`depict_molecule` is withdrawn. `DepictConfig` (`umol-io/src/depict.rs:25`) stays configuration for
layout generation only.

`Molecule::depict_layout` runs the verification of § 4 and then the existing lowering, whose `expect`
at `molecule.rs:398` becomes justified by that verification for every caller.
`Reaction::depict_layout` lowers each side with the molecule operation against the reaction
layout's side layouts and composes them where the layout places them (§ 2).

Adding an associated type and a required method is a breaking change for external `Depict`
implementors. `git grep 'impl Depict for'` on `98934189b` finds exactly two —
`umol-io/src/depict/molecule.rs:116` and `umol-io/src/depict/reaction.rs:84` — both in this crate,
and the trait exists to be implemented by umol's two depictable types, not by downstream code.

Feature boundary: `impl Depict for Molecule` and `MoleculeDepictionError::Layout` are gated on
`coordgen` (`molecule.rs:115`, `:129`–`131`), and `depiction = ["coordgen"]` (`umol-io/Cargo.toml:20`–`21`).
The supplied-layout path needs no CoordGen. This design keeps the trait impls under `depiction` as
today — splitting the feature so `depict_layout` compiles without a backend is possible but not
asked for by any consumer, and is left as an open question.

### 2. `ReactionLayout`

```rust
pub struct ReactionLayout {
    lhs: MoleculeLayout,
    rhs: MoleculeLayout,
    arrow_start: Point2D,
    arrow_end: Point2D,
}
```

One reaction coordinate system: every side position and both arrow endpoints are in it. The side
frames are the materialized sides' atom frames (`to_reaction_span()`, `reaction.rs:89`) — the
frames `depict_layout` checks against. A reaction built by `Reaction::from_sides` materializes to
sides in the input molecules' frames (observed from Python on 0.8.0: a 7 → 6 atom reaction renders
`reaction-lhs/atom/0..6` and `reaction-rhs/atom/0..5`); this document relies on that observation
for `from_sides` reactions only and does not claim it for every reaction.

**Construction.**

- `ReactionLayout::try_new(lhs: MoleculeLayout, rhs: MoleculeLayout, arrow_start: Point2D,
  arrow_end: Point2D) -> Result<Self, ReactionLayoutError>` — explicit; rejects a non-finite arrow
  endpoint (`ReactionLayoutError::NonFiniteArrow`) and a coincident start and end
  (`ReactionLayoutError::DegenerateArrow`); no frame check, because like `MoleculeLayout`
  (`layout.rs:62`, "carries no chemical attributes") the reaction layout is an open carrier not bound
  to a reaction. The side layouts are already finite by their own construction (`layout.rs:73`).
- `ReactionLayout::arrange(lhs: MoleculeLayout, rhs: MoleculeLayout) -> Result<Self,
  ReactionLayoutError>` — the arrangement step as its own operation: umol's relative placement of
  two *supplied* side layouts and its arrow. This is what `side_offset` (`reaction.rs:145`–`162`)
  and the constants `ARROW_HALF_LENGTH`/`SIDE_ARROW_GAP` (`:22`–`23`) do inside `compose_sides`
  today, lifted out of depiction into layout generation as the review asks. The returned layout
  has each side translated by its offset and the arrow at `(-0.75, 0)`–`(0.75, 0)`, as at `:57`–`60`.
  `Reaction::layout_with` (§ 3) is "generate both sides, then `arrange`"; an editor that already
  holds arranged sides calls `arrange` itself (§ Editor requirements, c).

**Editing.** `lhs_mut(&mut self) -> &mut MoleculeLayout`, `rhs_mut(&mut self)`, and
`set_arrow(&mut self, start: Point2D, end: Point2D) -> Result<(), ReactionLayoutError>` with the
same finiteness and degeneracy rules as `try_new`. Side positions are edited through the side
layout's own `set_position` (`layout.rs:117`), so there is one position-editing operation, not
three. Read access: `lhs(&self)`, `rhs(&self)`, `arrow_start()`, `arrow_end()`.

**Depiction semantics.** `Reaction::depict_layout` does not reposition. Each side is lowered in its
own frame and its items are translated by nothing — the positions already are reaction
coordinates; the arrow is drawn from `arrow_start` to `arrow_end`. The correspondence-pair number
offsets, which `mapping_index_offset` (`reaction.rs:232`) derives from an atom's neighbour
directions, remain derived from the supplied geometry at depiction time and are not stored in the
layout: they are rendering policy, not arrangement.

### 3. One layout-generation surface; `layout_molecule` retired

Today `pub fn layout_molecule(&Molecule, MoleculeLayoutAlgorithm) -> Result<MoleculeLayout,
LayoutError>` (`layout.rs:33`) is the only generator, and `Reaction` has none — its sides are laid
out inside `depict_with`. The review names a public free generator beside trait methods as an
inconsistency to correct rather than extend.

Proposed: generation lives on the same trait as depiction, as `layout` / `layout_with` in § 1.
`Molecule::layout_with` is `layout_molecule`'s body; `Reaction::layout_with` materializes the span,
lays out both sides with the config's algorithm, and `arrange`s them. `layout_molecule` becomes
`pub(crate)` (the CoordGen adapter in `umol-io/src/layout/coordgen.rs` is its only caller besides
depiction) or is folded into the impl; `LayoutError` (`layout.rs:45`) stays as the backend error and
maps into `MoleculeDepictionError::Layout` as it already does (`reaction.rs:95`, `:98`) — so
`Self::Error` is one type per implementor for generation and depiction alike.

Alternative considered: a separate `Layout` trait (`trait Layout { type Layout; type Error; fn
layout_with(...) }`) with `Depict: Layout`. Rejected for now: it splits `Self::Layout` and
`Self::Error` across two traits that must agree, and the point of § 1 is that the generated and
supplied paths are provably the same operation — which is easiest to read when `depict_with`'s
default body calls two methods of the same trait. A consumer that wants layout without depiction
still has it: `layout_with` does not require calling `depict_layout`.

Python: `Molecule.layout(*, algorithm=MoleculeLayoutAlgorithm.CoordGen())` as in the first revision,
and `Reaction.layout(*, algorithm=...)` returning a `ReactionLayout`; `RuntimeError` if the
backend fails, `ContradictionError` if the reaction cannot be materialized (the existing mapping of
`ReactionDepictionError::Materialization`, `reaction.rs:120`–`123`).

### 4. Supplied-layout verification, no repair

Today the only check a supplied layout can get is `check_frame`, atom count alone
(`layout.rs:145`–`156`), and doc 225 § 1 requires that it "must not silently acquire stereochemical
interpretation". The cis/trans predicate on generated coordinates lives in the CoordGen wrapper:
`validate_cis_trans_geometry` (`umol-coordgen-sys/src/lib.rs:360`) with `relative_side` / `half_plane`
(`:389`, `:416`) at a relative tolerance of `1e-6` (`:23`), rejecting `DegenerateCisTransGeometry`
(`:123`) and a side mismatch. `umol-geometric-core` exports `same_side_of_axis`, `signed_volume`,
`AXIS_SIDE_TOLERANCE` and `complementary_direction` (`umol-geometric-core/src/lib.rs:7`–`9`). The
selection of which stereo bonds are definite cis/trans, and in which ligand frame, already exists in
the layout adapter: `cis_trans_bond` (`umol-io/src/layout/coordgen.rs:60`) turns a `StereoBondView`
with a literal coset into a site bond, one actual ligand per end, and a `SameSide` / `OppositeSide`
relation (`:40`–`44`).

Proposed: a named verification operation at the depiction boundary,

```rust
impl Molecule {
    pub fn verify_layout(&self, layout: &MoleculeLayout) -> Result<(), MoleculeDepictionError>;
}
```

that `depict_layout` runs first and that an editor may call alone. It performs, in order:

1. **Frame** — `check_frame`, unchanged; `MoleculeDepictionError::LayoutFrame(MoleculeLayoutError)`
   as in the first revision.
2. **Definite cis/trans agreement** — for every stereo bond with a literal cis/trans coset, select
   the site and one actual ligand per end exactly as `cis_trans_bond` does (that selection moves to
   a shared helper; it is molecular interpretation and stays in `umol-io`), then classify the two
   ligands against the site axis with `same_side_of_axis`. Three outcomes, each its own variant
   carrying the **bond id**: agreement (continue); `CisTransMismatch { bond, expected, found }`
   where `expected` is `CisTransConfiguration` (`umol-graph-ir/src/ir/stereo.rs:1427`) — the
   coordinates draw the opposite relation; `CisTransDegenerate { bond, ligand }` — a ligand is
   within tolerance of the axis, so neither relation is drawn. Collinear covers both a ligand on the
   axis line and a zero-length site bond.
3. **Finite derived geometry** — the quantities the lowering and the SVG writer compute from
   positions must themselves be finite: bond unit vectors (a zero-length bond has none —
   `BondGeometry { bond }`), label anchor bounds and the mask rectangles derived from them, wedge
   polygons (base half-width and tip from the bond direction), and the pair-number offsets of
   `mapping_index_offset` (unit vectors of every incident bond). Finite input coordinates do not
   guarantee these: two atoms at `±1e308` give an infinite difference; coincident atoms give a
   `NaN` direction. The predicates — finite difference, non-zero length, normalizable direction —
   go in `umol-geometric-core`; `umol-io` names which entity failed
   (`NonFiniteGeometry { entity }`, with `entity` a bond, stereo atom or atom id).
4. **Tetrahedral wedges** — `tetrahedral_wedges` retained as today, `TetrahedralGeometry {
   stereo_atom }` unchanged (`molecule.rs:133`–`135`).

Nothing is repaired: no coordinate is moved, no stereo descriptor is re-perceived from geometry, no
degenerate ligand is nudged. A rejection names the entity and leaves both molecule and layout as
supplied.

Tolerance: the wrapper's `1e-6` relative to axis length times ligand length and geometric-core's
`AXIS_SIDE_TOLERANCE` (perpendicular distance against axis length) are, as 225 says, "not
interchangeable numerical contracts". This design uses whichever geometric-core predicate 225
settles on and does not choose between the two; until 225 settles it, `verify_layout` calls
`same_side_of_axis` as exported and the difference is confined to one call site.

`Reaction::depict_layout` runs `verify_layout` per side against the materialized side molecules, then
checks the arrow (finite, non-degenerate — already guaranteed by `ReactionLayout`'s constructors,
re-checked cheaply), mapping side failures through `LhsDepiction` / `RhsDepiction`
(`reaction.rs:124`–`129`) unchanged.

Test cases, as the review lists: for one definite cis/trans bond, coordinates that match, that are
opposite, and that put a ligand collinear with the site axis; a zero-length site bond; finite
coordinates whose derived quantities overflow (`1e308` and `-1e308` on one bond) and coincident
atoms; a frame mismatch; a supplied layout equal to the generated one depicting byte-identically to
`depict()`; a tetrahedral centre whose supplied coordinates admit no wedge. Each error variant is
asserted by pattern, with its entity id.

### 5. Python: mutable `MoleculeLayout`, `ReactionLayout`, error mapping

The Python API guide already gives the contract: "a Rust operation that deliberately mutates through
`&mut self` normally mutates the same Python object" (`docs/development/python-api.md:116`). `Deltas`
and `Edits` follow it — `#[pyclass(eq)]` without `frozen`, `append(&mut self, ...)` mutating the
wrapped Rust value (`umol-py/src/delta.rs:1792`, `:1820`; `umol-py/src/edit.rs:1365`, `:1399`).
`MoleculeLayout` follows the same pattern; the first revision's frozen class and functional
`with_position` are withdrawn.

| Surface | Behavior |
| --- | --- |
| `MoleculeLayout(positions)` | `positions` is a sequence of `(x, y)` floats in ascending atom-id order; `ValueError` on a NaN or infinite coordinate. Unchanged. |
| `.positions` | A **tuple** of `(x, y)` tuples in ascending atom-id order — immutable, so an assignment into it fails loudly rather than mutating a disconnected list. |
| `.position(atom_id)` | `(x, y)`; `IndexError` outside the frame. Unchanged. |
| `.set_position(atom_id, position)` | Mutates this layout in place; `None`. `ValueError` for an out-of-frame id or a non-finite position, and the layout is unchanged on failure (the Rust `set_position` returns before writing, `layout.rs:117`–`131`). The review's implementation, verbatim. |
| `len(layout)`, `==` | Atom count; structural equality. `__hash__` is `None` — a mutable value is unhashable, as `Deltas` and `Edits` are. |
| `Molecule.layout(*, algorithm=...)` | Generated layout; `RuntimeError` if the backend fails. |
| `Molecule.verify_layout(layout)` | `None` or raises; the mapping below. |
| `Molecule.depict_layout(layout)` | `verify_layout` then the lowering; the mapping below. |

`set_position` accepts any finite in-frame point and asks nothing of the molecule: the layout is an
open carrier and does not know its molecule. The stereo and geometry verdict comes only from
`verify_layout` / `depict_layout`. An editor may therefore make several moves and verify once.

`ReactionLayout` is bound the same way — `#[pyclass(eq)]`, unhashable — with
`ReactionLayout(lhs, rhs, arrow_start, arrow_end)` and `ReactionLayout.arrange(lhs, rhs)` (a static
method), `.arrow_start`, `.arrow_end` as tuples, `.set_arrow(start, end)`, and `.lhs` / `.rhs`.
The side accessors are the one design choice: the guide's property rules "apply recursively" and
"returning a live view solely to avoid a copy is not sufficient reason to introduce aliasing"
(`python-api.md:156`–`166`). Proposed: **`.lhs` / `.rhs` return copies**, and side positions are
edited through `set_lhs_position(atom_id, position)` / `set_rhs_position(atom_id, position)` on the
reaction layout itself, so a write always goes to the owner and a read never aliases it. A live
`MoleculeLayoutView` would need the `*View` role and an owner lifetime for a type whose whole state is
a `Vec<Point2D>`; the copy is cheaper to specify and to hold correct. Positions cross the boundary as
tuples; no `Point2D` class.

**Error mapping.** On the supplied path every rejection is a caller-argument error:

| Rust | Python | Why |
| --- | --- | --- |
| `LayoutFrame`, `CisTransMismatch`, `CisTransDegenerate`, `NonFiniteGeometry`, and `TetrahedralGeometry` **when raised from `verify_layout` / `depict_layout`** | `ValueError`, message naming the entity id | The caller supplied coordinates that do not fit the molecule. |
| `Layout` (backend failure) and `TetrahedralGeometry` **when raised from `depict` / `depict_with`** | `RuntimeError` | Operational failure of the backend or of umol's own generated geometry; the caller supplied nothing. |
| `ReactionDepictionError::Materialization` | `ContradictionError` | Unchanged. |
| `ReactionLayoutError` (arrow) | `ValueError` | Argument error at construction or `set_arrow`. |

This changes `TetrahedralGeometry`'s mapping on the supplied path from the first revision's
`RuntimeError`: under a supplied layout a centre that admits no wedge is the coordinates' fault, not
the backend's. The mapping is therefore by operation, not by variant alone; the binding for
`depict_layout` maps the whole `MoleculeDepictionError` to `ValueError` except `Layout`, which it
cannot produce.

## Editor requirements and conflicts

The consumer is moltracer's canvas: a browser page over one record that stores a molecule, per-atom
coordinates, and a left-hand-side snapshot with an atom correspondence, and draws with umol. Each
item states the requirement, its evidence, whether it conflicts with the surface above, and the
stance taken.

### (a) Storage and in-place mutation

moltracer stores coordinates per atom on its record and will build a `MoleculeLayout(positions)` per
draw; nothing in it threads one layout object through an edit session. In-place mutation costs it
nothing and the positions constructor is what it needs. **No conflict.** Requirement: the
`MoleculeLayout(positions)` constructor and `.positions` stay.

### (b) Incremental editing leaves an atom unpositioned

Record coordinates follow edits through umol's `MoleculeCorrespondence`; a removed atom's position
drops, survivors keep theirs. An **added** atom has no position. With a full-frame layout and no
repair, the stored arrangement cannot be drawn until something places the new atom — and § 4 says
umol will not. **This is the one true conflict** between "no implicit repair" and an editor's
workflow, and it is resolved on the editor's side: moltracer places a new atom a bond length from
its anchor in a free direction (geometry, no chemistry) before building the layout. The alternative
— constrained generation in umol, fixing the placed atoms and laying out the rest — would be a
`layout_with` variant taking partial positions; CoordGen has fixed-atom support natively
(UNVERIFIED for umol's wrapper). Recorded as an open question, **not requested**.

### (c) A reaction from two arranged sides

moltracer holds the LHS snapshot, the current molecule, their atom correspondence, and will hold the
LHS arrangement from mark time and the RHS arrangement live — two `MoleculeLayout`s in the sides'
own frames. It wants umol's relative placement and arrow; it does not want to compute them. **No
conflict provided `ReactionLayout::arrange` exists** (§ 2). A `ReactionLayout` constructible only by
full generation would force the editor to re-derive the arrow gap and side offsets itself, which is
depiction policy it should not own.

### (d) Refusal, not repair, on the page

A drag that puts a definite cis/trans geometry on the wrong side, or collapses it, is refused at
draw time; the page cannot repair and must tell the chemist **which bond** to move back. Requirement:
error variants carrying entity ids, surfaced in the Python exception message (§ 4, § 5). **No
conflict** — it is what § 4 specifies. Note, not asked: umol's SVG marks tetrahedral centres
(`molecule/stereo-atom/N`) but has no marker for alkene geometry, so the page cannot pre-warn which
bonds are stereo-definite without reading the molecule's stereo data. `verify_layout` as a public
operation lets it check a drop before drawing.

### (e) Multi-component molecules

A merged record is one multi-component `Molecule`; `Molecule.layout()` on `CC(=O)O.OCC` places the
fragments apart (observed on 0.8.0: acetic acid's bonds at `y ≈ 0`, ethanol's at `y ≈ 3`, one
`viewBox`). Requirement: none beyond that. A supplied layout with **overlapping** fragments is legal
geometry — overlap is neither a stereo nor a finiteness failure — and is **not** refused; a chemist
who drags one fragment onto another sees exactly that. **No conflict.**

### (f) What the editor does not need

Live per-mouse-move validation (one redraw per drop suffices, and `verify_layout` is cheap); a
`Point2D` class in Python (tuples stay); a functional `with_position` in Rust or Python (the first
revision's open question, withdrawn); constrained generation (b) now.

## Rejected alternatives

- **A layout field on `DepictConfig`.** `Reaction::depict_with` lays out two sides from one config;
  a single supplied layout has no meaning there, and `DepictConfig` would carry a value that is not
  a configuration.
- **A module-level `umol.layout_molecule`.** The Python surface has no free functions; every
  operation is a method on the value it acts on. The generator belongs on `Molecule` and `Reaction`.
- **A public `depict_molecule` / `depict_reaction` free function** (the first revision). Withdrawn
  per review: the operations belong on `Depict`, and a free function beside a trait method is the
  inconsistency `layout_molecule` already shows.
- **A frozen Python `MoleculeLayout` with `with_position`** (the first revision). Withdrawn: the
  mutation contract exists (`python-api.md:116`) and `Deltas` / `Edits` already follow it.
- **A `Point2D` class.** Positions cross the boundary as `(x, y)` tuples; nothing in the bindings
  needs a richer point type yet.
- **A separate `Layout` trait.** § 3.
- **Live `.lhs` / `.rhs` views on `ReactionLayout`.** § 5.

## Consequences

- Rust exact tests: `depict_layout` on a generated layout equals `depict` byte for byte, for a
  molecule and for a reaction; the verification cases of § 4; `ReactionLayout::arrange` reproduces
  the offsets and arrow `compose_sides` produces today (a fixture reaction depicts identically before
  and after the lift); `try_new` and `set_arrow` rejections.
- Python tests: construction, lookup, in-place update and its failures leaving the object unchanged;
  `.positions` is a `tuple`; `hash(layout)` raises `TypeError`; the generator for both types; the
  round trip that generates, moves one atom, depicts, and observes only the incident bond endpoints
  move; the frame mismatch and each § 4 rejection as `ValueError` whose message contains the entity
  id; `ReactionLayout` side copies do not alias.
- The two tests from the first implementation that survive unchanged in intent are the frame check
  and the round trip.
- Release notes are left to the release owner.

## Open questions

- Constrained layout generation (fixed atoms, the rest generated) — the editor's (b); not requested.
- `ReactionLayout.lhs` / `.rhs` as copies with `set_*_position` on the owner (proposed) versus a
  declared live view.
- The degeneracy tolerance for the cis/trans predicate — deferred to 225 § 1.
- Whether `verify_layout` is public or only `depict_layout`'s first step. Proposed public, for (d).
- Whether `depict_layout` should compile without the `coordgen` feature (§ 1).
