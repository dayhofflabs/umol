# 222 — Depicting a supplied layout

Status: In Progress
Date: 2026-09-08
Relates: [221](221-depiction-api-2026-09-03.md),
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

## Reversed portion of doc 221

Doc 221 settled that `MoleculeLayout`'s "public existence does not require a public operation for
lowering an arbitrary supplied layout into a depiction", removed the layout-frame mismatch from
`MoleculeDepictionError` "because the public depiction operation no longer accepts an independently
supplied layout", and stated that the layout API "is not an alternate public constructor for
`Depiction`". Those three statements are reversed here. Everything else in 221 stands: `Depict`,
`DepictConfig`, the opaque `Depiction`, `render_svg`, the `depiction` feature boundary, and the
reaction path are unchanged.

## Settled surface

### Rust, `umol-io`

- `umol_io::depict::depict_molecule(molecule: &Molecule, layout: &MoleculeLayout) ->
  Result<Depiction, MoleculeDepictionError>` is a public free function re-exported beside
  `MoleculeDepictionError`, parallel to `layout::layout_molecule`. It checks frame agreement with
  `MoleculeLayout::check_frame` and then performs the existing lowering.
  `Depict::depict_with` for `Molecule` calls it; a generated layout satisfies the check trivially.
- `MoleculeDepictionError::LayoutFrame(MoleculeLayoutError)` is re-added. It is reachable only
  from `depict_molecule` and carries `MoleculeLayoutError::FrameSizeMismatch`.

The crate-private lowering previously read positions with an `expect` whose justification was
frame agreement established by the single generating caller. A public entry point taking a
caller-supplied layout must establish that agreement itself; without the check, a Python call could
abort the interpreter. The check is a contextual-validity precondition of the consuming operation,
in the sense of the data-type guide: `MoleculeLayout` is an open carrier not bound to a molecule,
so the first public operation that combines the two establishes their agreement.

`MoleculeLayout`, `MoleculeLayoutError`, `layout_molecule`, `LayoutError`, `DepictConfig`,
`Depict`, `Depiction`, and reaction depiction are untouched. Reactions remain CoordGen-only.

### Python, `umol-py`

`MoleculeLayout` is bound as a frozen newtype over the Rust layout with `from_rust`/`to_rust`:

| Surface | Behavior |
| --- | --- |
| `MoleculeLayout(positions)` | `positions` is a list of `(x, y)` floats in ascending atom-id order; `ValueError` on a NaN or infinite coordinate. |
| `.positions` | List of `(x, y)` tuples in ascending atom-id order. |
| `.position(atom_id)` | `(x, y)`; `IndexError` outside the frame. |
| `.with_position(atom_id, position)` | A new layout with one position replaced; the receiver is unchanged. `ValueError` for an out-of-frame id or non-finite position. |
| `len(layout)`, `==` | Atom count; structural equality. |
| `Molecule.layout(*, algorithm=MoleculeLayoutAlgorithm.CoordGen())` | The Python form of `layout_molecule`. `RuntimeError` if the backend fails. |
| `Molecule.depict_with_layout(layout)` | The Python form of `depict_molecule`. `ValueError` for `LayoutFrame`, `RuntimeError` for `TetrahedralGeometry`. |

The Python layout is immutable. `set_position` is deliberately not bound: a mutating Python API is
still to be designed, and the existing Python surface has no mutable value types. `with_position`
gives editors a functional update in the meantime.

`umol-py` gains `umol-geometric-core` as an optional dependency under the `depiction` feature to
construct `Point2D` from Python tuples. The dependency was already compiled transitively through
`umol-io`.

### Error mapping

`LayoutFrame` maps to `ValueError` because it is a caller-supplied argument that does not fit the
receiver — the same class the bindings use for `MoleculeLayout` construction, `Correspondence`
construction, and integrity failures. `Layout` and `TetrahedralGeometry` remain `RuntimeError`,
as before: they are operational failures of the backend or the geometry, not argument errors.

`position` raises `IndexError` rather than `ValueError` because it is a positional lookup into a
dense frame and follows Python sequence semantics, matching the other id lookups in the bindings.

## Rejected alternatives

- **A layout field on `DepictConfig`.** `Reaction::depict_with` lays out two sides from one config;
  a single supplied layout has no meaning there, and `DepictConfig` would carry a value that is not
  a configuration.
- **A module-level `umol.layout_molecule`.** The Python surface has no free functions; every
  operation is a method on the value it acts on. The generator therefore belongs on `Molecule`.
- **Binding `set_position` directly.** Would introduce the first mutable value type in `umol-py`
  ahead of a settled mutating API.
- **A `Point2D` class.** Positions cross the boundary as `(x, y)` tuples; nothing in the bindings
  needs a richer point type yet.

## Consequences

- Rust exact tests for `depict_molecule` cover frame mismatch and agreement with `Depict::depict`
  on a generated layout.
- Python tests cover construction, lookup, functional update, and their failures; the generator;
  and a round trip that generates a layout, moves one atom, depicts, and observes the incident bond
  endpoints move in the SVG while every other endpoint is unchanged. The frame-mismatch test
  asserts a `ValueError` is raised rather than the interpreter aborting.
- Release notes are left to the release owner; the repository has no unreleased section at the
  time of this change.

## Open questions

- The shape of a mutating layout API in Python, if one is wanted.
- Whether the Rust `MoleculeLayout` should gain a `with_position` companion so the two languages
  offer the same functional update.
