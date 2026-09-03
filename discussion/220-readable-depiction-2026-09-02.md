# 220 — Readable molecular depiction

Status: In Progress
Date: 2026-09-02
Relates: [201](201-molecular-data-first-steps-2026-08-19.md)

## Purpose

The first depiction path established the separation between graph IR, two-dimensional layout,
format-neutral depiction, SVG rendering, and Python rich display. Its output is inspectable but is
not yet a generally readable structure diagram. Bonds run through visible atom labels, aromaticity
is represented by unrelated circles at atom and bond positions, and every stereo entity is reduced
to the same circle marker without showing its configuration.

This work should make ordinary molecule and reaction SVGs readable without turning `umol-io` into a
general drawing application. The intended visual direction is a restrained black-and-white,
ChemDraw-like bond-line diagram. CoordGen remains the layout backend; this work concerns the
projection from graph IR into drawing geometry and the fixed SVG realization of that geometry.

## Boundaries

The required outcome is:

- ordinary carbon skeletons are legible and visible atom labels do not have bonds drawn through
  them;
- localized single, double, and triple bonds have conventional geometry;
- an explicit `AromaticSystem` remains distinguishable as a system without using centered ring
  circles or pretending that a selected Kekule assignment is the stored representation;
- supported definite atom and bond stereo configurations are visible rather than marked only as
  stereo sites; and
- SVG output remains deterministic, scalable, monochrome, and linked to graph-IR entities through
  its existing structured references.

This scope does not include colored element palettes, raster output, interactive editing,
highlighting, user-defined themes, correspondence-aware reaction alignment, or exhaustive drawing
conventions for every graph-IR overlay and constraint. It also does not change CoordGen selection
into an implicit default.

## Current evidence

The current implementation has three relevant structural limitations.

1. `AtomItem` carries an unstructured label at an atom center while `BondItem` carries center-to-
   center endpoints. The SVG renderer has neither endpoint visibility nor label extents, so it
   cannot shorten a bond around a visible label.
2. Aromaticity is lowered to one dotted circular `MarkerItem` for every aromatic atom and bond.
   Those marks neither trace a system nor distinguish neighboring aromatic systems.
3. Stereo atoms and stereo bonds are lowered to identical circular markers. The lowering reads the
   stereo site but not the kind, ligand frame, or coset. The CoordGen adapter likewise projects only
   atomic number, localized-bond endpoints, and bond order, even though the vendored backend can
   accept tetrahedral and E/Z information and can choose stereo display bonds.

The checked-in reference implementations separate chemical interpretation from final drawing but
span a wide complexity range:

- RDKit measures atom-label glyph rectangles and intersects every bond line with them. This handles
  long and structured labels accurately but requires a font-metrics subsystem.
- CDK constructs atom-symbol outlines and convex hulls before generating bonds. Its delocalized
  display uses either ring donuts or an offset dashed bond stroke; its own documentation warns that
  donuts can mislead fused systems and recommends a Kekule display where that representation is
  appropriate.
- Indigo maintains label bounds, detailed bond-end geometry, and separate stereo bond geometry. It
  is a useful quality reference but not a reasonable implementation-size target.
- Open Babel shortens a bond by a fixed fraction when an endpoint has a visible label, suppresses
  ordinary carbon labels, and renders selected wedge, hashed-wedge, and unknown-stereo bond styles.
  The fixed shortening is less exact but demonstrates a small implementation that remains readable.

SVG also supplies a narrower solution that does not require measuring text. A mask can remove the
parts of the bond layer beneath duplicates of the atom-label glyphs, expanded slightly for
clearance. The masked bond is absent rather than overpainted with a presumed background color, and
browser font substitution affects the mask and visible label identically.

These implementations do not supply aromatic semantics for umol. In particular, a conventional
ring circle or a renderer-selected Kekule form would erase the identity and possibly the extent of
the explicit aromatic systems represented by graph IR.

## Preliminary direction

### Fixed style before configuration

Do not add an empty or single-purpose depiction configuration merely to reserve an API shape. The
current work has one intended style, and fixed nominal-bond-length ratios are sufficient for it.
`Depict::depict_with` should continue to receive the explicit layout algorithm directly, while
`svg::render` should continue to render one stable house style.

The monochrome color should remain `currentColor` rather than hard-coded per-element colors. This
keeps the drawing black on an ordinary page while allowing an embedding document to select one
foreground color. It does not constitute an atom color scheme, and no palette API is proposed.

A config becomes justified only when this work identifies at least two independently useful
operational choices that callers actually need. Renderer geometry such as font size, bond width,
multiple-bond separation, wedge width, and viewport margin belongs together if that point is
reached; chemical acceptance or aromatic interpretation does not belong in that config.

### Atom labels and bond clearance

The first target should use ordinary bond-line visibility: omit an undecorated carbon symbol and
its implicit hydrogens when the carbon is represented by the skeleton, while retaining labels for
heteroatoms and carbons decorated by an isotope, nonzero charge, or nonzero unpaired-electron
count. An isolated neutral carbon is the exceptional visible skeleton atom and is labeled `C`, not
`CH4`. Isotopes, hydrogen counts on otherwise visible labels, and charges remain visible when their
underlying fields are literal. The isotope precedes the element symbol; charge and literal
unpaired-electron counts follow the atom, with the latter rendered as radical dots. Lone pairs and
multiplicity remain outside this cut.

For the current SVG output, bonds should remain center-to-center in `Depiction` and the SVG renderer
should mask the bond layer beneath visible atom labels. The mask should use duplicates of the same
text with enough expansion to provide clearance and cover counters and gaps within a label. It must
not paint a background-colored box. This delegates glyph realization to the browser and handles all
bond shapes uniformly without adding font measurement or label bounds to the public scene model.

This is deliberately an SVG presentation fix, not a general solution to label collision in a
format-neutral depiction. Other renderers may eventually justify explicit label geometry, but that
problem is outside the present scope.

Atom-local radical dots belong beside the visible atom label or, for an omitted skeleton carbon,
beside the corresponding vertex. Charge and unpaired electrons borne by an `AromaticSystem` must
not be assigned visually to one member atom. Their initial annotation anchor should be the
unweighted geometric centroid of the member-atom positions. If that point overlaps drawing
geometry or lies outside a concave system contour, the annotation should move to the nearest clear
interior position, with a deterministic position just outside the upper-right of the contour as a
fallback.

### Aromatic systems

Remove aromatic circle markers. The leading experiment is one inward, offset dashed contour around
the outer boundary of each explicit aromatic system. Internal fusion bonds receive no dash. This
uses one system-level contour rather than independently choosing a side for every bond, avoiding
the ambiguous offset side on a bond shared by two fused rings. The contour retains
`AromaticSystemId` references and does not claim a renderer-selected alternating bond assignment.

For a crossing-free planar layout, determine the boundary geometrically from the aromatic induced
subgraph: sort incident system edges by polar angle, walk the embedding faces using directed edges,
identify the unbounded face by signed area, and offset that face boundary inward. Because an
aromatic system is biconnected, its unbounded face has a cycle boundary in a valid planar embedding.
The convention is unambiguous even when several aromatic systems occur because graph IR permits an
atom to belong to at most one of them.

This is an experiment, not yet a universal convention. If a crossed, degenerate, self-intersecting,
or cage-like projection does not yield a suitable system boundary, the aromatic decoration is
omitted rather than replaced with a misleading ring-local mark. C60 demonstrates the limit: the
planar outer face may describe only one projected ring and therefore fail to communicate the
system's full extent. Ring-local circles in a three-dimensional cage projection do not preserve the
graph-IR system semantics either. Atom- or bond-local aromatic constraints without a corresponding
`AromaticSystem` are assertions rather than resolved system membership and remain outside the
first readable projection.

### Stereo

Remove the generic stereo circles. The graph-IR representation is sufficient to depict definite
tetrahedral and cis/trans configurations; no additional chemical interpretation or CIP assignment
is required. This conclusion applies to materialized stereo entities with literal ground cosets,
not to an unresolved configuration expression.

A tetrahedral stereo atom supplies its site, an ordered four-ligand frame, its tetrahedral kind,
and its coset. Implicit hydrogen and lone-pair positions are explicit virtual ligands borne by the
site rather than values the renderer must infer. Molecule integrity checks the frame arity,
uniqueness, and ligand incidence. `StereoAtomView::coset_for` can reframe the configuration into a
renderer-selected ligand order. The existing TableIR raise already computes the coset from one
front or back wedge, the two-dimensional ligand positions, a complementary virtual-ligand
direction, and signed volume. Depiction can perform the inverse: choose one explicit incident bond
and select solid or hashed direction so that the resulting winding reproduces the stored coset.
The choice of which eligible bond carries the wedge is an aesthetic and collision-avoidance choice,
not missing stereochemical semantics.

A cis/trans stereo bond supplies the double-bond site and an ordered four-ligand frame partitioned
as two ligands at one endpoint followed by two at the other. Integrity permits either orientation
of the site bond but enforces the two-by-two endpoint partition. After reframing through
`StereoBondView::coset_for`, coset zero is the same-side class named Z and coset one is the
opposite-side class named E in that ligand frame. Thus the renderer can select one explicit ligand
on each endpoint and pass their required same-side or opposite-side relation to the layout backend.

Unlike a tetrahedral wedge, definite cis/trans stereo is carried by the substituent coordinates.
The present CoordGen boundary cannot preserve it because it accepts only atomic numbers, bond
endpoints, and localized bond orders. The vendored backend already has atom- and bond-stereo input
and can constrain E/Z layout, so the leading direction is a minimal extension of that projection
and native boundary. The returned coordinates must then be checked against the requested
configuration; drawing an ordinary double bond over unconstrained coordinates is not an acceptable
fallback. Tetrahedral wedges may be selected from the returned geometry in graph-IR lowering, or a
later boundary refinement may return CoordGen's chosen display bond, but the SVG renderer should
receive an already selected wedge and merely realize it.

Undetermined E/Z will not use crossed double bonds. `StereoCoset::Undetermined`, a set, or an open
term is omitted and must never be lowered to a definite wedge or definite E/Z geometry. Axial,
square-planar, trigonal-bipyramidal, and octahedral configurations are likewise omitted rather than
being misrepresented as tetrahedral wedges. The supported first-cut stereo scope is therefore only
materialized tetrahedral and cis/trans entities with a literal ground coset.

## Public contract

No public depiction configuration is added. `MoleculeLayout` remains an open, finite coordinate
carrier in one dense atom frame. `Depiction` remains an operation-issued scene with no public
aggregate constructor. Graph-IR lowering selects chemical display semantics and source references;
the SVG renderer realizes those items without reading graph IR. Python continues to expose only the
explicit layout-algorithm choice and an already rendered `Svg` value.

### Native cis/trans input

**Type and role:** `umol_coordgen_sys::CisTransBond` is the minimal bond-stereo input accepted by
the native coordinate operation. It contains `bond`, the index of the double bond in the supplied
`Bond` slice; `first_ligand`, an actual substituent of `Bond::atom_0`; `second_ligand`, an actual
substituent of `Bond::atom_1`; and `relation: SideRelation`.
`SideRelation::{SameSide, OppositeSide}` names only the requested relative placement of that
selected pair. It does not encode CIP priority, E/Z nomenclature, or a graph-IR ligand frame.

**Open carrier or operation-issued value:** `CisTransBond` is an open FFI-input carrier with public
fields, like `Bond`. Its indices have no meaning without the atom and bond slices supplied to
`generate_coordinates`.

**Intrinsic representation invariants:** Rust construction establishes only valid field types and
one of the two closed relation values. `CisTransBond` is `#[repr(C)]`; `SideRelation` is
`#[repr(u8)]` with `SameSide = 0` and `OppositeSide = 1`. These representations do not prove any
index or incidence property.

**Contextual properties and supplied context:** at `generate_coordinates`, the site-bond index and
both atom indices must be in frame, the site must have order two, each ligand must be distinct from
the site endpoints and connected to its named endpoint by a supplied bond, and at most one
`CisTransBond` may name a site bond. These checks occur before entering C++. The C entry point also
rejects a numeric relation other than the two values declared by the header.

**Semantic predicates and validators:** private safe-boundary validation establishes the preceding
properties. After CoordGen returns finite points, the signed half-planes of the two selected
ligands about the site-bond axis must realize `SameSide` or `OppositeSide` as requested. A fixed,
scale-relative tolerance classifies a collinear ligand or zero-length site bond as inconsistent
geometry rather than choosing a side arbitrarily.

**Public constructors:** no constructor is added. Callers use the public fields. The public
`generate_coordinates` signature gains `cis_trans_bonds: &[CisTransBond]`; this is an intentional
breaking change to the experimental native boundary, not a parallel compatibility entry point.
The C header mirrors the record, closed relation values, and added input slice without introducing
another semantic representation.

**Conversions and preserved information:** the graph-IR adapter considers only a materialized
`StereoBondView` of kind `CisTrans` with `StereoCoset::Lit`. It selects one actual ligand at each
endpoint, asks `coset_for` for that concrete four-position frame, and converts coset zero/Z to
`SameSide` and coset one/E to `OppositeSide`. Reversing the site orientation and swapping the two
selected ligands preserves the request. Undetermined, set-valued, term-valued, and non-cis/trans
configurations are normal absence for this projection and are not errors.

**Failure, absence, and panic behavior:** invalid indices, order, incidence, or duplicate site
records return a specific `CoordgenError`. Finite output that is degenerate or has the wrong
relative placement returns `CoordgenError` as well. None of these cases panic or silently discard a
definite supported request.

**Algebraic, preservation, or roundtrip properties:** atom order and output frame remain preserved;
empty stereo input retains the present coordinate-generation behavior; endpoint reversal with the
corresponding ligand reversal is equivalent; and every successful stereo request satisfies the
same signed-half-plane predicate used by the postcondition check.

**Rust/Python boundary:** these native input records are not Python API. Python reaches them only
through graph-IR molecule depiction.

### Depiction geometry and fallibility

**Type and role:** ordinary multiple-line bonds remain `BondItem`. `WedgeItem` carries `tip`,
`base`, `kind: WedgeKind::{Solid, Hashed}`, and structured references; the tip is the stereocenter
in an issued depiction. `DashedContourItem` carries ordered `points`, `closed`, and structured
references. `DepictionItem` gains `Wedge` and `DashedContour` variants. The obsolete `MarkerItem`,
`MarkerKind`, and `DepictionItem::Marker` are removed after all marker producers have been replaced.

**Open carrier or operation-issued value:** the individual item structs are open presentation
carriers with public fields, consistently with the existing items. `Depiction` remains
operation-issued and keeps its private aggregate constructor.

**Intrinsic representation invariants:** the open item structs do not claim contextual validity.
Graph-IR lowering issues wedges with distinct finite endpoints and contours with finite points in
drawing order. Bounds include both wedge endpoints and every contour point. A reaction translation
moves every such point and changes molecule references to the corresponding reaction-side frame
without changing item kind, order, or other references.

**Public constructors and conversions:** no item constructor or public `Depiction` constructor is
added. Molecule lowering converts one selected localized single-bond item into a wedge item and
adds the stereo-atom reference; it does not draw an ordinary bond beneath the same wedge. Aromatic
lowering emits at most one closed dashed contour per explicit system and carries its
`AromaticSystemId`. The SVG renderer maps the supplied wedge and contour geometry to fixed-style
SVG and remains infallible.

**Contextual consumer and failure boundary:** molecule lowering is the first consumer combining a
`Molecule` with a `MoleculeLayout`. A new public `MoleculeDepictionError` has
`Layout(LayoutError)`, `LayoutFrame(MoleculeLayoutError)`, and
`TetrahedralGeometry { stereo_atom: StereoAtomId }` variants; the first remains feature-gated with
CoordGen. `depict` and `Depict for Molecule` use this one boundary. `DepictFromSidesError` replaces
its separate layout/frame variants with `LhsDepiction(MoleculeDepictionError)` and
`RhsDepiction(MoleculeDepictionError)`, retaining its two correspondence-frame variants.
`ReactionDepictionError` retains `Materialization` and replaces its layout variants with the same
lhs/rhs depiction variants. Python maps materialization contradictions as it does now and maps all
layout or depiction failures to `RuntimeError`.

**Literal extraction and omission:** depiction needs only an operation-specific literal subset,
not a completely ground molecule. Nonliteral atom label fields are omitted individually; radical
dots read the literal unpaired-electron `count` directly and do not require or display a literal
multiplicity. Only materialized tetrahedral and cis/trans entities with a single literal coset are
depicted.
Undetermined, set-valued, term-valued, or unsupported stereo and atom- or bond-local aromatic
assertions are intentionally omitted. A supported definite stereo entity is never silently
omitted: native E/Z mismatch is a `CoordgenError`, and unusable tetrahedral geometry is a
`MoleculeDepictionError`. Failure to construct a trustworthy aromatic outer contour is ordinary
absence because omission is the settled non-misrepresentation policy for that experimental mark.

**Semantic properties:** tetrahedral wedge selection, when interpreted by the same winding rule as
TableIR raise, reproduces the stored coset in the selected ligand frame. Cis/trans coordinates
satisfy the requested relative-side predicate. Aromatic contours reference exactly one explicit
system and never introduce a localized Kekule assignment. Rendering preserves item references and
is deterministic for an equal depiction.

**Rust/Python boundary:** the new item and error types remain Rust API. No Python scene-model
constructors are added; `Molecule.depict_with` and `Reaction.depict_with` retain their signatures.

## Settled first-cut decisions

- renderer constants remain fixed and monochrome through `currentColor`;
- ordinary skeleton carbons are omitted, but an isolated neutral carbon is labeled `C`;
- isotope precedes the element, while charge and radical dots follow the atom;
- bonds are cleared behind labels with an SVG mask, never a background-colored patch;
- unsuitable aromatic contours, local aromatic assertions, unresolved stereo, and unsupported
  stereo kinds are omitted;
- definite cis/trans inconsistency and unusable definite tetrahedral geometry are errors; and
- mask clearance, contour offset, dash pattern, wedge dimensions, and deterministic wedge-bond
  ranking are fixed implementation constants tuned against the evidence fixtures, not open public
  configuration.

## Staged implementation plan

### S0 — Evidence corpus and additive scene vocabulary

- [x] **S0a — Benchmark fixtures** `[dep: none]`: extend `umol-io/benches/layout.rs` with paired
  literal Z/E cases and `umol-io/benches/svg.rs` with labeled atoms, tetrahedral stereo, a fused
  aromatic system, and a reaction-sized mixture. Record the pre-change Criterion baseline before
  changing either algorithm. The cases are measurement inputs, not claims that the current marker
  output is correct.

  The pre-change Criterion baseline on 2026-09-02 used CoordGen 3.0.2. Intervals are Criterion's
  reported 95% confidence intervals:

  | Layout case | Time |
  | --- | ---: |
  | acyclic/asymmetric_tree_8 | 62.924–63.180 µs |
  | cyclic/cyclooctane | 10.606–10.794 µs |
  | aromatic/benzene | 8.050–8.170 µs |
  | disconnected/mixed_components_8 | 9.856–9.899 µs |
  | underdetermined/wildcard_path_8 | 73.562–75.446 µs |
  | cis_trans/z_but_2_ene | 4.272–4.299 µs |
  | cis_trans/e_but_2_ene | 4.237–4.245 µs |
  | mapping_hard_tail/high_symmetry_complete_7 | 1.300–1.327 ms |
  | mapping_hard_tail/repeated_components_3x2 | 7.850–8.028 µs |

  | SVG case | Time |
  | --- | ---: |
  | chain/8 | 3.632–3.645 µs |
  | chain/128 | 67.342–68.316 µs |
  | representative/labeled_atoms | 0.856–0.859 µs |
  | representative/tetrahedral_stereo | 1.769–1.775 µs |
  | representative/fused_aromatic | 9.998–10.049 µs |
  | representative/mapped_reaction | 6.265–6.350 µs |

  Commands: `cargo bench -p umol-io --bench layout --features coordgen -- --noplot` and
  `cargo bench -p umol-io --bench svg -- --noplot`.
- [x] **S0b — Scene vocabulary** `[dep: none]`: add `WedgeKind`, `WedgeItem`,
  `DashedContourItem`, and the two additive `DepictionItem` variants. Extend bounds, reference
  access, reaction translation, and SVG realization in the same subitem so every new variant has
  complete consumers. Retain the old marker variant temporarily while its producers still exist.
- [x] **S0c — Contract gate** `[dep: S0a, S0b]`: inspect every changed public symbol against the two
  contract sheets above; verify that `Depiction::from_items` remains private and no config,
  unchecked aggregate constructor, Python scene type, or compatibility layout entry point has
  appeared.

S0 is green when the existing depiction behavior is preserved, new item rendering and translation
have exact Rust tests, touched tests use `rstest`, and the two benchmark baselines are recorded.

S0 completed on 2026-09-02. The public-symbol audit found only the settled open item carriers and
enum variants: `WedgeKind`, `WedgeItem`, `DashedContourItem`, `DepictionItem::Wedge`, and
`DepictionItem::DashedContour`. `Depiction::from_items` remains private; no rendering config,
aggregate constructor, Python scene type, or compatibility layout entry point was added. Exact
bounds, references, reaction translation, solid/hashed wedge SVG, and open/closed contour SVG tests
pass. Verification used `cargo test -p umol-io -q`,
`cargo test -p umol-io --features coordgen -q`, and
`cargo clippy -p umol-io --all-targets --features coordgen -- -D warnings`.

### S1 — Cis/trans-aware CoordGen layout

- [x] **S1a — Native input boundary** `[dep: S0]`: add `SideRelation` and `CisTransBond`, change the
  Rust and C coordinate-generation signatures, validate index/order/incidence/uniqueness before
  FFI, validate the C relation discriminator, and update all empty-stereo callers. Add exact ABI
  tests for every invalid-input category.
- [x] **S1b — Backend application and postcondition** `[dep: S1a]`: attach the selected substituents
  and relation to CoordGen's bond stereo record, establish its absolute stereo before coordinate
  generation, and verify the returned half-planes in the safe Rust wrapper. Test same-side,
  opposite-side, endpoint reversal, determinism, degeneracy classification, and mismatch
  detection without using CoordGen itself as the assertion oracle.
- [x] **S1c — Graph-IR projection** `[dep: S1b]`: in `umol-io/src/layout/coordgen.rs`, select actual
  endpoint ligands, reframe with `StereoBondView::coset_for`, and emit native requests only for
  `CisTrans` plus `StereoCoset::Lit`. Add layout fixtures proving Z and E geometry, reversed stored
  site orientation, implicit-H frames, and omission of undetermined, set, term, and unsupported
  stereo forms.
- [x] **S1d — Layout evidence** `[dep: S1c]`: rerun the layout benchmark cases and retain the comparison
  with the S0 baseline. Treat success rate and stereo postcondition failures as evidence alongside
  time per molecule.

S1 is green when `umol-coordgen-sys` native tests and `umol-io` layout tests pass with the CoordGen
feature and every successful literal cis/trans fixture satisfies an independently computed
relative-side predicate.

The S1 layout benchmark completed all 9/9 groups, including both definite cis/trans cases, with no
layout or stereo-postcondition errors. The intended stereo projection and postcondition cost is
about 0.4 us per four-atom definite cis/trans molecule; the remaining movement is ordinary
whole-benchmark run variation rather than a shared added cost:

| Case | S1 time per molecule, 95% interval | Criterion change from S0, 95% interval |
| --- | ---: | ---: |
| acyclic/asymmetric_tree_8 | 61.817–62.027 us | -2.04% to -1.49% |
| cyclic/cyclooctane | 10.496–10.782 us | -1.32% to +0.99% |
| aromatic/benzene | 8.012–8.044 us | -2.51% to +0.50% |
| disconnected/mixed_components_8 | 9.941–10.103 us | +1.24% to +3.02% |
| underdetermined/wildcard_path_8 | 71.151–71.578 us | -6.65% to -4.04% |
| cis_trans/z_but_2_ene | 4.605–4.878 us | +7.74% to +10.25% |
| cis_trans/e_but_2_ene | 4.619–4.637 us | +8.77% to +9.18% |
| mapping_hard_tail/high_symmetry_complete_7 | 1.293–1.301 ms | -4.15% to -0.46% |
| mapping_hard_tail/repeated_components_3x2 | 7.648–7.670 us | -4.51% to -2.58% |

Command: `cargo bench -p umol-io --bench layout --features coordgen -- --noplot`.

S1 completed on 2026-09-02. The public-symbol audit found only the settled
`SideRelation`, `CisTransBond`, extended `CoordgenError`, and breaking
`generate_coordinates` input boundary. The record remains an open `#[repr(C)]` carrier; contextual
index, order, incidence, uniqueness, and output-geometry invariants are checked by the one safe
operation. No constructor, compatibility entry point, generic graph-IR projection helper, or
Python scene API was added. Verification used
`cargo test -p umol-coordgen-sys --features native -q`,
`cargo test -p umol-io --features coordgen -q`, and Clippy with `-D warnings` for both crates and
all feature-relevant targets.

### S2 — Atom labels and bond clearance

- [x] **S2a — Label projection** `[dep: S0]`: replace the current all-elements label rule with skeleton
  carbon suppression and the isolated-carbon exception. Preserve literal isotopes and implicit H
  on visible labels, append literal charge and radical dots, and read the unpaired-electron count
  independently of the undisplayed multiplicity. Omit each nonliteral field without requiring
  complete atom groundness. Cover heteroatoms, decorated carbons, isolated carbon, multi-digit
  values, and omitted/open fields with exact item tests.
- [x] **S2b — SVG label mask** `[dep: S2a]`: add one deterministic SVG mask derived from the same atom
  glyphs used visibly and apply it to molecular bond, wedge, and contour strokes. Use mask expansion
  to cover glyph interiors and provide clearance without assuming a page background. Verify XML
  structure, escaping, references, reactions, nonwhite embedding semantics, and exact empty output.
- [x] **S2c — Rendering evidence** `[dep: S2b]`: tune only the fixed mask expansion against the S0
  label fixtures, then rerun the SVG benchmark to quantify the extra mask and duplicate-text cost.

S2 is green when ordinary carbon skeletons remain connected and legible, every visible label clears
underlying bonds without a background-colored element, and rendering remains deterministic.

The fixed mask expansion is 0.30 nominal bond lengths. Rendering the label duplicate with a 0.60
stroke width fills glyph counters and leaves visible clearance on a nonwhite background without
adding a visible background-colored element. The S2 SVG benchmark completed all 6/6 groups:

| Case | S2 time, 95% interval | Criterion change from S0, 95% interval |
| --- | ---: | ---: |
| chain/8 | 2.538–2.545 us | -30.28% to -29.91% |
| chain/128 | 50.036–51.035 us | -25.97% to -24.00% |
| representative/labeled_atoms | 1.472–1.477 us | +72.21% to +73.13% |
| representative/tetrahedral_stereo | 2.841–2.865 us | +59.42% to +60.64% |
| representative/fused_aromatic | 8.390–8.431 us | -16.21% to -15.66% |
| representative/mapped_reaction | 7.604–7.646 us | +18.64% to +21.02% |

The carbon-only cases are faster because S2a removes their atom text nodes. Cases retaining labels
pay for one mask plus one duplicate text glyph per visible atom; the smallest labeled fixture adds
about 0.62 us. Command: `cargo bench -p umol-io --bench svg -- --noplot`.

S2 completed on 2026-09-02. The public-symbol audit found no new type, constructor, configuration,
or error boundary. `Depiction` remains operation-issued, label geometry remains renderer-private,
and the existing `depict` and `render` operations implement the settled projection and SVG-output
changes directly. Element is the only required atom literal; each isotope, hydrogen, charge, and
unpaired-electron field has independent omission semantics, and multiplicity is not inspected.
Exact label, XML, reference, reaction, mask, nonwhite-background, and empty-output tests pass.
Verification used `cargo test -p umol-io --features coordgen -q`,
`cargo clippy -p umol-io --all-targets --features coordgen -- -D warnings`, and
`cargo +nightly fmt --all -- --check`.

### S3 — Definite stereo depiction

- [x] **S3a — Tetrahedral wedge selection and error boundary** `[dep: S0]`: introduce
  `MoleculeDepictionError` and migrate direct molecule depiction, reaction-side propagation, the
  `Depict` implementations, and Python error mapping within this breaking subitem. Collect literal
  tetrahedral stereo atoms, choose distinct explicit display bonds deterministically across
  adjacent stereocenters, and use the supplied two-dimensional coordinates to select `Solid` or
  `Hashed` so the inverse winding reproduces the stored coset after `coset_for`. Replace the
  selected ordinary `BondItem` rather than overlaying it. Return `MoleculeDepictionError` when the
  geometry cannot establish a valid wedge.
- [x] **S3b — Stereo omission and references** `[dep: S1, S3a]`: stop producing generic stereo markers.
  Verify that literal cis/trans stereo is visible in the S1 coordinates, literal tetrahedral stereo
  carries both bond and stereo-atom references, and unresolved or unsupported stereo emits neither
  a definite mark nor an error. Exercise reaction translation and lhs/rhs error propagation.
- [x] **S3c — Stereo laws and evidence** `[dep: S3b]`: add a feature-gated
  `depiction_property` test target for ligand reframing and wedge-to-coset recovery, plus
  correctness fixtures for virtual H/lone-pair frames, rings, adjacent stereocenters, and both
  wedge kinds. Rerun the layout and SVG benchmark subsets.

S3 is green when every supported definite stereo entity is either represented consistently or
reported as an error, while every unresolved or unsupported stereo form is omitted deliberately.

S3 completed on 2026-09-03. Literal tetrahedral configurations now replace one deterministically
selected explicit single bond with a solid or hashed wedge carrying both bond and stereo-atom
references. Selection reframes the stored configuration into an actual-ligands-first display
order, validates the emitted winding against the TableIR convention, and assigns distinct bonds
across adjacent stereocenters. Unusable definite geometry returns
`MoleculeDepictionError::TetrahedralGeometry`; layout and depiction failures propagate through the
corresponding molecule, reaction-side, reaction, and Python `RuntimeError` boundaries. Definite
cis/trans remains visible through the S1 coordinate constraint, while unresolved and unsupported
stereo emits no generic marker.

The feature-gated `depiction_property` target exercises stored-frame invariance over all
degree-four ligand permutations and independently recovers both tetrahedral cosets from wedges
over translated, scaled, rotated, and reflected nondegenerate layouts. Exact fixtures cover solid
and hashed wedges, virtual hydrogen, lone pair, a ring site, adjacent stereocenters, reaction-side
translation, and lhs/rhs errors. The public-symbol audit found only the settled
`MoleculeDepictionError` and the settled breaking variants of `DepictFromSidesError` and
`ReactionDepictionError`. No public constructor, depiction configuration, scene re-export, or
Python scene API was added; `Depiction` remains operation-issued.

The S3 evidence rerun completed all 9/9 layout groups and 6/6 SVG groups. Layout has no S3 algorithm
change and remained within ordinary run variation. The real wedge adds about 0.64 us to the
tetrahedral SVG fixture relative to S2; its absolute render time remains about 3.49 us.

| Layout case | S3 time per molecule, 95% interval |
| --- | ---: |
| acyclic/asymmetric_tree_8 | 61.529–61.859 us |
| cyclic/cyclooctane | 10.401–10.479 us |
| aromatic/benzene | 7.948–8.033 us |
| disconnected/mixed_components_8 | 9.690–9.737 us |
| underdetermined/wildcard_path_8 | 70.698–71.123 us |
| cis_trans/z_but_2_ene | 4.566–4.601 us |
| cis_trans/e_but_2_ene | 4.569–4.624 us |
| mapping_hard_tail/high_symmetry_complete_7 | 1.279–1.286 ms |
| mapping_hard_tail/repeated_components_3x2 | 7.683–7.780 us |

| SVG case | S3 time, 95% interval | Criterion change from preceding run, 95% interval |
| --- | ---: | ---: |
| chain/8 | 2.449–2.458 us | -3.70% to -3.21% |
| chain/128 | 49.219–50.168 us | -3.91% to -0.85% |
| representative/labeled_atoms | 1.428–1.446 us | -3.48% to -2.63% |
| representative/tetrahedral_stereo | 3.481–3.499 us | +22.73% to +23.63% |
| representative/fused_aromatic | 8.036–8.142 us | -4.48% to -3.72% |
| representative/mapped_reaction | 7.390–7.423 us | -2.97% to -2.41% |

Commands: `cargo bench -p umol-io --bench layout --features coordgen -- --noplot` and
`cargo bench -p umol-io --bench svg -- --noplot`. Verification used the `umol-io` default and
CoordGen test suites, the feature-gated `depiction_property` target, Clippy over all `umol-io`
targets with CoordGen and proptest, `umol-py` Rust tests and Clippy with depiction in the repository
Python 3.13 environment, the rebuilt focused Python depiction tests, and
`cargo +nightly fmt --all -- --check`.

### S4 — Aromatic-system contours and annotations

- **S4a — Outer-face extraction** `[dep: S0]`: build the explicit system's induced edge geometry,
  reject segment crossings and degenerate rotations, sort incident half-edges by angle, walk all
  faces, select the unbounded face by signed area, and offset it inward. Treat a system with no
  degree-two member in its induced graph as cage-like for this first cut, so C60 is omitted rather
  than decorated as one apparent ring. Emit one closed `DashedContourItem` or no item; never fall
  back to atom/bond dots, ring circles, or a Kekule assignment.
- **S4b — System annotations** `[dep: S4a, S2]`: render literal aromatic-system charge and unpaired
  electrons as text anchored first at the unweighted member centroid, moved through a fixed
  deterministic set of clear interior candidates, then to the upper-right exterior fallback.
  Atom-local radicals remain attached to atom labels or skeleton vertices. Nonliteral system fields
  are omitted individually.
- **S4c — Aromatic fixtures and marker removal** `[dep: S4a, S4b]`: cover a single ring, fused rings,
  multiple disjoint systems, concave and crossed layouts, and the C60 conformance molecule as a
  stress case. Remove aromatic marker production, omit local aromatic assertions, then remove
  `MarkerItem`, `MarkerKind`, and `DepictionItem::Marker` after all Rust, reaction, SVG, and Python
  consumers have migrated. Rerun the aromatic SVG benchmark subset.

S4 is green when trustworthy systems receive exactly one referenced contour, unsuitable systems
receive none, and no generic marker API or marker-shaped aromatic/stereo output remains.

### S5 — Integrated verification and closeout

- **S5a — End-to-end fixtures** `[dep: S2, S3, S4]`: exercise representative molecules and a mapped
  reaction through graph IR, CoordGen, format-neutral depiction, SVG, and Python rich display.
  Assert semantic SVG structure and exact small stable fragments; inspect generated SVGs from
  temporary output for readability without checking generated artifacts into `materials/`.
- **S5b — Performance and quality gate** `[dep: S5a]`: compare final Criterion results with S0,
  format the workspace, run native CoordGen and feature-enabled `umol-io` tests and clippy, then run
  the depiction-enabled Python tests through the repository Python 3.13 environment. Run the full
  workspace gate after the narrow checks pass.
- **S5c — API and document closeout** `[dep: S5b]`: repeat the public-symbol contract audit, update
  rustdoc to current behavior, record the implemented fixture and benchmark evidence here, change
  this document and `discussion/000-status.md` to `Completed`, and leave deferred renderer
  generalization or additional stereo conventions as separately scoped work rather than implied
  behavior.

The critical path is S0 -> S1 -> S3 -> S5. S2 and S4 depend on the S0 scene contract but not on the
native stereo implementation; they may be developed independently before S5. Nothing after S0 is
purely optional for the readable-depiction outcome. Color schemes, public rendering configuration,
font metrics, non-SVG collision geometry, unknown-stereo marks, local aromatic assertions, lone-pair
and multiplicity labels, and additional stereo kinds are explicitly deferrable follow-up work.
