# 221 — Depiction API

Status: In Progress
Date: 2026-09-03
Relates: [220](220-readable-depiction-2026-09-02.md)

## Purpose

Doc 220 established readable molecule and reaction drawings, but its experimental operation path
was promoted to public API without a separate API review. The rendering behavior is feature-complete;
the public interface is not ready for release.

This work must define one coherent Rust API, preserve that operation's meaning in Python, and make
the resulting SVG usable outside notebook display. It does not reopen layout quality, drawing
semantics, styling, color, or renderer selection.

## Current surface

With the `coordgen` feature, `umol_io::depict::Depict` is an extension trait implemented for
`Molecule` and `Reaction`. Its only method,
`depict_with(MoleculeLayoutAlgorithm) -> Result<Depiction, _>`, performs layout and format-neutral
depiction.

The lower-level Rust path is a collection of free functions:

- `layout::layout_molecule` generates a `MoleculeLayout`;
- `depict::molecule::depict` combines a molecule and layout;
- `depict::reaction::depict_from_sides` combines two molecules, two layouts, and an atom
  correspondence;
- `depict::reaction::depict_from_sides_with` also generates both layouts; and
- `svg::render` converts a `Depiction` to `String`.

`Depiction` is a local, operation-issued boundary type with private aggregate construction and
publicly inspectable scene items and source references.

Python uses the same `depict_with` name for a different operation: it performs Rust depiction and
SVG rendering, then returns a Python-only `Svg`. That type exposes its text only through
`_repr_svg_`, the Jupyter rich-display protocol.

## Findings

### Operation and feature boundaries

`Depict` is already the extension trait required to provide method syntax for graph-IR
`Molecule` and `Reaction` values owned by another crate. The problem is not the use of an
extension trait. The trait is coupled to `MoleculeLayoutAlgorithm` and therefore disappears with
the CoordGen implementation feature.

`depict_with` also implies an ordinary `depict` operation that does not exist. The name does not
state whether the supplied value is an algorithm, layout, or broader configuration. The high-level
operation needs an operation-specific configuration boundary. Its initial configurable field is
the layout algorithm, whose only current variant is CoordGen.

### Boundary ownership

The separation of layout, format-neutral depiction, and SVG rendering remains correct. Requiring
the primary path to compose unrelated free functions does not follow from that separation.
`Depiction` is local to `umol-io` and is the centralized boundary between lowering and rendering;
it can own the operation that produces SVG.

The explicit-layout and explicit-side functions expose implementation stages rather than an
independently settled public capability. In particular, `depict_from_sides` accepts independently
supplied molecules, layouts, and a correspondence and therefore needs a contextual
`DepictFromSidesError`. A `Reaction` already contains the lhs and resolved deltas needed to
materialize both sides and recover their correspondence. Generated layouts likewise establish
their molecule frames. The high-level reaction path therefore does not need the independent-input
operation or its error type.

### Rust and Python semantics

The same public method name must not return a format-neutral scene in Rust and a rendered format in
Python. Either Python exposes the same depiction operation and result role, or its combined
convenience operation receives an explicitly SVG-specific name.

A rendered SVG value must expose its text through an ordinary API. `_repr_svg_` is an additional
display protocol, not the data-access interface.

## Settled direction

- Keep layout, depiction, and rendering as distinct operations and representations.
- Use `Depict` as the extension trait implemented by `Molecule` and `Reaction`.
- Add `depict()`, which uses `DepictConfig::default()`, and make `depict_with` accept
  `&DepictConfig` rather than a layout-algorithm selector.
- Make `Depiction` an opaque, operation-issued result. Its scene items, bounds, geometry,
  typography, and source-reference carriers are implementation details, not public data types or
  inspection methods.
- Add the inherent `Depiction::render_svg() -> String`. The normal Rust path is
  `value.depict()?.render_svg()` rather than another free-function call.
- Remove the public molecule-lowering, independent-side composition, and free SVG-rendering
  functions. Remove `DepictFromSidesError` entirely; it has no remaining public operation to
  describe.
- Retain `MoleculeLayout` as a separate public coordinate representation. Its public existence does
  not require a public operation for lowering an arbitrary supplied layout into a depiction.
- Give equal Rust and Python operation names equal semantics and result roles.
- Bind `DepictConfig` and the opaque `Depiction` in Python. Python `depict()` and `depict_with()`
  return `Depiction`; `render_svg()` returns ordinary SVG text, and `_repr_svg_()` delegates to the
  same rendering operation for notebook display. Remove the Python-only `Svg` result.
- Gate the complete high-level depiction surface on CoordGen while it is the only available layout
  backend. With the feature disabled, `Depict`, `DepictConfig`, `Depiction`, and their implementations
  are absent rather than present with no usable backend. The independent layout representation may
  remain available outside that gate.
- Do not add color schemes, renderer selection, or new drawing behavior in this work.
- Correct the API before its first release; no compatibility layer for the unreleased surface is
  required.

## Public contract

### `Depict`

**Role:** extension operation producing a format-neutral `Depiction` from a molecule or reaction.

**Operations:** `depict()` uses `DepictConfig::default()`. `depict_with(&DepictConfig)` uses the
supplied operation configuration. Both are implemented for `Molecule` and `Reaction`.

**Failure boundary:** layout and lowering failures must remain distinguishable at the operation that
cannot produce a depiction. Reaction materialization may additionally fail because `Reaction`
construction does not require its deltas to materialize a two-sided span. A reaction depiction
requires no separately supplied rhs or correspondence.

**Feature boundary:** the trait and its implementations are available only when the high-level
depiction capability has its CoordGen backend.

### `DepictConfig`

**Role:** operation-specific configuration selected by `depict_with` and defaulted by `depict`.

**Boundary:** `layout_algorithm` is public and defaults to
`MoleculeLayoutAlgorithm::CoordGen`. A one-variant algorithm enum is retained so the operational
choice remains explicit and can grow without changing the configuration shape. The initial API
does not invent choices for color or rendering.

### `Depiction`

**Role:** operation-issued, format-neutral scene retaining drawing order, geometry, and source
references.

**Boundary:** aggregate construction and the complete scene IR remain private. Public consumers do
not receive `DepictionItem`, item records, `WedgeKind`, `DepictionReference`, `Bounds`, `items()`, or
`bounds()`.

**Rendering:** `render_svg()` returns the complete SVG document fragment as `String`. Rendering
remains a distinct operation, but its public entry point belongs to the boundary value it consumes.

### Low-level layout and side composition

`MoleculeLayout` and its checked coordinate operations remain a separate public capability. The
current `depict::molecule::depict`, `depict::reaction::depict_from_sides`,
`depict::reaction::depict_from_sides_with`, and `svg::render` functions are removed from the public
surface.

The internal reaction pipeline materializes a `ReactionSpan`, projects its lhs and rhs, recovers the
atom correspondence, lays out both sides, lowers them, and composes the result. Span materialization
establishes correspondence-frame agreement, and layout generation establishes layout-frame
agreement. Final side composition is therefore infallible after the preceding fallible operations;
`DepictFromSidesError` and its frame-mismatch variants are not retained.

Public molecule and reaction depiction errors remain because they are the associated errors of the
public `Depict` implementations. Their variants describe failures reachable through those
operations rather than independently supplied private lowering inputs.

### Python

Python binds the same `DepictConfig` and operation-issued `Depiction` roles. `Molecule.depict()`,
`Molecule.depict_with(config)`, `Reaction.depict()`, and `Reaction.depict_with(config)` return the
Python `Depiction`. `Depiction.render_svg()` returns `str`, which can be written with ordinary Python
file APIs. `Depiction._repr_svg_()` returns the same SVG for notebook rich display.

Published Python artifacts enable the depiction capability. Builds without the Python `depiction`
feature do not expose these types or methods.

## Target public surface

With the Rust `coordgen` feature, `umol_io::depict` exposes only `Depict`, `DepictConfig`,
`Depiction`, `MoleculeDepictionError`, and `ReactionDepictionError`. `Depiction` has no public
constructor, conversion, scene accessor, or structural-equality contract; its public operation is
`render_svg()`. `DepictConfig` has public `layout_algorithm: MoleculeLayoutAlgorithm` and implements
`Default` with CoordGen selected.

`MoleculeDepictionError` retains failures for CoordGen layout and unsupported definite
tetrahedral display geometry. Its current layout-frame mismatch is removed because the public
depiction operation no longer accepts an independently supplied layout. `ReactionDepictionError`
retains materialization, lhs-depiction, and rhs-depiction failures. The independent-side functions
and `DepictFromSidesError` are removed rather than re-exported or replaced.

The separate `umol_io::layout` surface remains public, including `MoleculeLayout`,
`MoleculeLayoutError`, explicit `layout_molecule`, `MoleculeLayoutAlgorithm`, and `LayoutError`.
This API produces and edits coordinates; it is not an alternate public constructor for
`Depiction`.

With the Python `depiction` feature, `umol` exports `MoleculeLayoutAlgorithm`, `DepictConfig`, and
`Depiction`. The current `Svg` Python class is removed. Depiction methods are present on `Molecule`
and `Reaction` only with that feature.

## Implementation plan

### S0 — Add the replacement vocabulary

#### S0a — Rust depiction configuration and default operation **Done**

**Module:** `umol-io::depict`

Add `DepictConfig` with public `layout_algorithm`, retaining the one-variant
`MoleculeLayoutAlgorithm`, and default it to CoordGen. Add `Depict::depict()` as the default-config
operation while retaining the current `depict_with` signature temporarily. Test that `depict()`
agrees exactly with the current explicit CoordGen path for representative molecule and reaction
inputs and preserves their exact error results.

**Change:** additive (green).  
**Dependencies:** none.

#### S0b — SVG rendering on `Depiction` **Done**

**Modules:** `umol-io::depict`, `umol-io::svg`

Add `Depiction::render_svg() -> String`, initially delegating to the existing renderer. Move the
renderer contract and exact SVG assertions to the method while the free function remains as a
temporary implementation entry point.

**Change:** additive (green).  
**Dependencies:** none.

#### S0c — Python replacement boundary types **Done**

**Modules:** `umol-py::depict`, `umol-py::lib`, `umol` package exports

Add Python `DepictConfig` with a public one-variant `layout_algorithm` constructor argument and
getter, defaulting to CoordGen, and an opaque, non-user-constructible `Depiction` owning the Rust
value. Expose `render_svg() -> str` and `_repr_svg_()`, with exact tests that both operations return
the same complete SVG. Register and export the new classes alongside the old `Svg` temporarily.

**Change:** additive (green).  
**Dependencies:** [dep: S0a, S0b].

**Stage gate:** the focused CoordGen-enabled `umol-io` suite passes. With `umol-py/.venv` activated
and confirmed as Python 3.13, the depiction-enabled `umol-py` Rust tests pass, followed by
`maturin develop --features depiction` and the focused Python depiction and import tests.

**S0 evidence:** `cargo test -p umol-io --features coordgen` passed 3,405 unit and 15 integration
tests. Under Python 3.13.15, the depiction-enabled `umol-py` Rust suite passed 1,646 tests with two
ignored; a fresh `maturin develop --features depiction` completed, and the focused Python depiction
and import suite passed 76 tests.

### S1 — Cut over the primary Rust and Python operations

#### S1a — Config-based Rust `Depict` implementations **Done**

**Modules:** `umol-io::depict`, `umol-io::depict::molecule`,
`umol-io::depict::reaction`, `umol-io::bin::depict_dsl`

Change `Depict::depict_with` to accept `&DepictConfig`; have `depict()` delegate through
`DepictConfig::default()`. Update the `Molecule` and `Reaction` implementations so the config owns
the high-level layout choice, and update the example binary to use the method chain ending in
`render_svg()`. Preserve the current molecule layout, tetrahedral geometry, reaction
materialization, and side-specific failure meanings.

**Change:** breaking (red until S1b migrates the Python caller).  
**Dependencies:** [dep: S0a, S0b].

#### S1b — Rust-equivalent Python methods **Done**

**Modules:** `umol-py::depict`, `umol-py::lib`, `umol` package exports and Python depiction tests

Add no-argument `depict()` and change `depict_with` to accept `DepictConfig` on both `Molecule` and
`Reaction`; both return Python `Depiction`. Remove the Python-only `Svg` result, including its native
registration, package export, and tests. Retain `MoleculeLayoutAlgorithm` as the public type of
`DepictConfig.layout_algorithm`, and retain the present Python exception classes and mapping for the
corresponding Rust operation errors.

**Change:** breaking migration (restores green).  
**Dependencies:** [dep: S0c, S1a].

**Stage gate:** the CoordGen-enabled Rust API tests and the depiction-enabled Rust/Python binding
suites pass. Python tests cover molecule and reaction `depict`, explicit default configuration,
ordinary SVG text access, file-write-compatible text, notebook display, and exact materialization
and lowering failures.

**S1 evidence:** `cargo test -p umol-io --features coordgen` passed 3,405 unit and 15 integration
tests. Under Python 3.13.15, `cargo test -p umol-py --lib --features depiction` passed 1,647 tests
with two ignored; a fresh depiction-enabled `maturin develop` completed, and the focused Python
depiction and import suite passed 78 tests. Clippy passed with warnings denied for all targets in
both affected crates. A workspace source search found no remaining Python `Svg` symbol or caller
passing a layout algorithm directly to `Depict::depict_with`.

### S2 — Seal the scene representation

#### S2a — Reaction-owned side composition and errors **Done**

**Modules:** `umol-io::depict::reaction`, `umol-io::depict::molecule`

Replace `depict_from_sides` and `depict_from_sides_with` with an internal reaction-composition path
whose inputs come only from `Reaction::to_reaction_span` and generated side layouts. Remove
`DepictFromSidesError` and its frame-mismatch cases. Retain exact reaction tests for materialization,
lhs lowering, rhs lowering, correspondence indices, arrows, stereo, and bond orders, but construct
the subject through `Reaction::depict` rather than independently supplied sides.

Remove `MoleculeDepictionError::LayoutFrame`; the independent public `MoleculeLayoutError` remains
unchanged for the public layout API. Reconcile every remaining public depiction-error variant with
a reachable `Depict` failure before completing the subitem.

**Change:** breaking (red until S2c migrates the remaining external callers).  
**Dependencies:** [dep: S1a].

#### S2b — Opaque `Depiction` and private lowering/rendering **Done**

**Modules:** `umol-io::depict`, `umol-io::depict::molecule`,
`umol-io::depict::reaction`, `umol-io::svg`, `umol-io::lib`

Remove the public scene accessors and external visibility of `DepictionItem`, all item records,
`AtomLabel`, `WedgeKind`, `DepictionReference`, and `Bounds`. Make molecule lowering, reaction
composition, and SVG rendering internal modules and functions, re-exporting only the target public
depiction symbols. Gate that complete surface on `coordgen`; leave the independent layout module
and coordinate types under their existing feature behavior. Remove `Depiction`'s current derived
`Clone`, `Debug`, and `PartialEq` implementations; no copying, debug-inspection, or scene-equality
contract has been settled for the opaque result.

Keep exact item-order, bounds, typography, masking, reference-encoding, arrow, aromatic, and stereo
tests inside the owning crate. Do not add test-only public accessors.

**Change:** breaking (red until S2c migrates the remaining external tests and benchmarks).  
**Dependencies:** [dep: S2a].

#### S2c — Unit-test and benchmark migration **Done**

**Modules:** `umol-io` depiction property suite and SVG benchmark

Move the tetrahedral scene-law tests into the owning crate so they can inspect the private IR
without widening production visibility. Exercise the two tetrahedral cosets and all 24 ligand
frames deterministically, and cover the finite rotation and reflection classes with targeted unit
tests. Rewrite the SVG benchmark to hold opaque depictions produced outside the timed loop and
benchmark `Depiction::render_svg()`; replace its independently supplied reaction sides with an
actual `Reaction`. Keep the separate public layout benchmark unchanged.

**Change:** breaking caller migration (green).  
**Dependencies:** [dep: S2b].

**Stage gate:** `umol-io` tests pass with `coordgen` and `proptest`; the SVG and layout benchmarks
compile and run; clippy passes for all `umol-io` targets with both features. A workspace search
finds no external use of retired item types, accessors, side-composition functions, free SVG
rendering, `DepictFromSidesError`, or Python `Svg`.

**S2 evidence:** `cargo test -p umol-io --features 'coordgen proptest'` passed 3,397 unit tests,
including the moved tetrahedral depiction tests, plus 15 layout and six SMILES property
tests. Quick Criterion runs completed for both the opaque SVG-rendering benchmark and the unchanged
layout benchmark. Clippy passed with warnings denied for all `umol-io` targets under both features;
the feature-disabled crate and the depiction-enabled `umol-py` Rust suite also passed. A workspace
source search found no external use of the retired scene types, inspection methods, side-composition
functions, free renderer, `DepictFromSidesError`, or Python `Svg`.

### S3 — Release-facing verification and closeout

#### S3a — Public documentation and release surface

**Modules:** Rust rustdoc, Python docstrings and exports, `RELEASE_NOTES.md`

Describe only the final method-based Rust and Python APIs, including feature availability,
configuration defaulting, reaction materialization failure, SVG text access, and notebook display.
Replace release-note descriptions of the unreleased interface rather than documenting a migration
from it. Audit every public symbol and failure variant against the target list above; do not retain
compatibility aliases or hidden public construction seams.

**Change:** additive documentation and surface audit (green).  
**Dependencies:** [dep: S1b, S2c].

#### S3b — Final verification and document closeout

**Modules:** workspace verification and this discussion record

Run `cargo +nightly fmt --all -- --check`; the focused CoordGen/proptest tests, clippy, and
benchmarks; and the depiction-enabled `umol-py` Rust and Python suites using the repository Python
3.13 environment and a fresh `maturin develop`. Then run the Python-activated workspace test and
clippy gates. Inspect generated public documentation and the final feature-disabled build. Record
the resulting commands and evidence, mark every plan subitem complete, and close doc 221 only when
the implementation and status index agree.

**Change:** verification and closeout (green).  
**Dependencies:** [dep: S3a].

The critical path is S0a/S0b -> S0c -> S1a -> S1b -> S2a -> S2b -> S2c -> S3a -> S3b. No stage is
deferrable for the 0.7.1 depiction release; S0 is additive preparation, while S1 and S2 perform the
two independent public-surface cutovers and restore a green workspace at each stage boundary.

## Release consequence

The 0.7.1 release notes and publication preparation are provisional until this contract is
implemented. The `umol-perm` cleanup remains internal and is not part of this scope.
