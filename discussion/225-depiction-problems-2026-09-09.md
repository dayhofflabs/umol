# 225 — Depiction problems

Status: Proposed
Date: 2026-09-09
Relates: [220](220-readable-depiction-2026-09-02.md),
[221](221-depiction-api-2026-09-03.md),
[224](224-smiles-ring-closure-frame-2026-09-08.md),
[153](153-format-parsing-outstanding-tasks-2026-07-18.md)

## Scope

Collect depiction problems identified while visually verifying stereo fixtures after doc 224.
These are separate from parser work and fixture identity corrections. This document records
current evidence, required direction, and unresolved design; it is not an implementation plan.

## 1. Coordinate verification and geometry ownership

The CoordGen wrapper currently rejects returned coordinates that do not satisfy the requested
cis/trans relation. In umol-coordgen-sys/src/lib.rs, generate_coordinates calls
validate_cis_trans_geometry after native execution and the finite-point check. The wrapper's
relative_side and half_plane implement the geometric predicate and choose a relative tolerance
of 1e-6. Both degenerate geometry and a definite mismatch become CoordgenError variants.

Rejection of incorrect definite stereo was approved during doc 220. That does not establish
the wrapper as the appropriate owner of depiction policy or geometric algorithms. Returning
the native coordinates and deciding whether they are suitable for depiction are separate concerns.

### Required boundary

All umol geometry operations belong in umol-geometric-core, not umol-io and certainly not
umol-coordgen-sys. This includes the predicates, numerical tolerances, and geometric constructions
used in verification and rendering. Vendored CoordGen retains its own native algorithms; this is
about where umol implements geometry around the backend.

The io layer may select graph-IR entities and ligand frames, invoke geometric operations, compare
their results with the requested stereochemistry, and apply depiction policy. It must not implement
the geometric calculations itself. A backend-independent verification operation associated with
layout is a possible orchestration point. MoleculeLayout::check_frame currently checks only atom
count, and must not silently acquire stereochemical interpretation.

umol-geometric-core already provides same_side_of_axis, signed_volume, and complementary_direction.
The wrapper duplicates side classification, but its tolerance is relative to axis length times
ligand length, whereas same_side_of_axis compares perpendicular distance with axis length.
These are not interchangeable numerical contracts. Settle the intended degeneracy criterion before
consolidating them; simply moving the wrapper formula would leave two competing definitions.

### Wrapper inspection

Inspection covers the umol-owned Rust wrapper, C++ shim, header, and build configuration;
it is not an audit of all vendored CoordGen internals.

| Behavior | Assessment |
| --- | --- |
| Returned cis/trans geometry classification, tolerance, and rejection | Geometry and depiction policy in the wrong layer; separate from native execution. |
| Bond endpoint, stereo-site, and ligand bounds; raw pointer and enum checks | Native boundary validation, not drawing-quality policy. |
| Stereo site must be double; ligands must be distinct from site endpoints and incident on their designated endpoint; duplicate stereo sites rejected | Additional semantic input restrictions in the safe wrapper. They are not merely memory-safety checks. Reconcile each with the native stereo-input contract when narrowing the wrapper. |
| Non-finite returned points rejected | Numerical output validation, not a stereo-quality judgment; explicitly settle whether finite output remains a wrapper guarantee. |
| Mutex, allocation/exception translation, ownership cleanup, input-order coordinate extraction | Native execution and resource handling. |
| C++ stereo assignment before coordinate generation | Translation of supplied native stereo input; no additional output geometry check or coordinate repair found in the shim. |

No second drawing-quality rejection or coordinate-repair pass was found in the umol-owned wrapper.
The IO adapter's generic atom/bond projection and bond-length normalization are outside the wrapper;
they must not be confused with native behavior. The inspection does not establish that every
existing wrapper restriction was individually sanctioned.

The same ownership problem extends into umol-io/src/depict/molecule.rs: segment distance,
orientation/intersection predicates, polygon containment, area, and boundary offsets are implemented
there. Tetrahedral wedge interpretation already calls geometric-core for volume and complementary
direction. Audit and consolidate the remaining geometry, including SVG stroke construction and
layout transformations, under the same boundary; leave chemical selection and rendering policy
with their owning operations.

### Cyclooctene evidence

E-cyclooctene resolves successfully, but CoordGen's nine-atom MACROCYCLE threshold excludes its
double bond from native stereo handling. The wrapper rejects the resulting OppositeSide mismatch.
Changing the threshold to eight makes E render, but Z then fails the SameSide check. With that check
bypassed, Z has a nearly collinear ring bond, about 0.005 degrees from the double-bond line.
The macro also selects ring-layout algorithms; lowering it is not an isolated stereo fix.
Both experiments were reverted. Doc 153 retains the fixture verification record.

Remaining decisions: the verification surface and diagnostics; which operation requires it;
the treatment of near-degenerate coordinates; and how to handle a backend that cannot supply
suitable geometry without silently misrepresenting stereo.

## 2. Wedges and hashes on substituents

Stereo marks on ring bonds or the main chain make the fixtures substantially harder to read.
The required drawing convention is to place wedges and hashes on substituent bonds, keeping the
ring and main chain plain. Geometrically decodable stereo alone is not adequate readability.

Current selection is in umol-io/src/depict/molecule.rs. tetrahedral_candidates considers real
ligands attached by localized single bonds and retains geometrically valid wedge choices.
It ranks them by whether the other endpoint is a tetrahedral stereo site, then by bond index.
assign_distinct_wedge allocates distinct bonds across centers. There is no ring or main-chain
preference, and the index tie-break can select an inconvenient scaffold bond.

Settle substituent selection and what constitutes the displayed main chain, including fused rings
and adjacent centers competing for a bond. Also settle the case with no suitable explicit
substituent, including whether an implicit hydrogen should be shown for stereo depiction.
Preserve the stereo frame when selecting a different marked bond; do not merely move the glyph.
Use the reviewed ring and chain fixtures to assess readability as well as encoded configuration.

## 3. Orientation of multiatom labels

When a bond connects from the right, orient the label so the attachment atom is next to the bond:
HO–, H₃C–, and ⁻O– rather than the corresponding right-facing forms. Preserve subscripts,
superscripts, and chemical grouping; this is not character-string reversal.

Currently molecule::atom_label receives no layout information. It appends hydrogen to the element
symbol, places the hydrogen count in right_subscript, and charge/radicals in right_superscript.
svg::render_atom emits left superscript, base, right subscript, right superscript in fixed order,
with the whole text centered. This cannot represent H₃C correctly by merely swapping the base text.

The label representation and placement need to express attachment direction and positioned groups.
Update label bounds and bond masking with the chosen arrangement so bonds end at the attachment
atom, not the hydrogen or charge. Settle ambiguous directions and multiple connections, and verify
left/right examples with hydrogen counts, charges, and isotopes. This concerns visible labels;
it does not require exposing every skeletal carbon as a methyl label.
