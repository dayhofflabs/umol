# umol 0.9.0

Changes since 0.8.0. This release adds molecular and reaction SMILES output,
depiction from supplied layouts, and SVG rendering controls. It corrects stereo
parsing and separates isotope completion from valence resolution. It includes
breaking Rust and Python API changes.

## SMILES output

- Rust adds export::export_smiles and export::export_reaction_smiles in
  umol-graph. Their _with variants accept explicit IO, chemistry-model, and
  resolve configurations, in the same order as the ingestion functions.
- Python adds Molecule.to_smiles and Reaction.to_reaction_smiles, returning
  strings. Both accept the keyword options io_config, chemistry_model, and
  resolve_config used by their corresponding input methods.
- Default input and output use OpenSMILES syntax, the SMILES valence model, and
  Natural isotope policy. Output supports ordinary atoms, isotopes, charges,
  bracket-H radical notation, neutral closed-shell aromatic systems,
  tetrahedral stereo, and local double-bond stereo.
- Graph reaction output assigns one-based atom-map labels in correspondence
  order. Original numeric labels are not retained; unmatched atoms are unlabelled.
- The Rust export::Convey trait converts a graph model to an output boundary
  after projecting a private copy. Smiles and ReactionSmiles add from_table_ir
  construction and render/render_with methods. Boundary construction retains
  the table as supplied; rendering checks whether it can be expressed as SMILES.

## Stereo parsing and resolution

- Correct tetrahedral SMILES interpretation at ring-closing stereocenters,
  including the ordering of explicit atoms, implicit hydrogens, and lone pairs.
  Basic and extended parsers now retain explicit stereo frames instead of
  reconstructing them from bond-table order.
- Normalize slash/backslash evidence into explicit double-bond frames. Conflicting
  markers and markers outside the supported local double-bond domain fail during
  parsing.
- MOL wedges retain their narrow endpoint when atom indices are reordered.
  Either wedges and CXSMILES wiggly-bond annotations preserve unknown stereo.
  MOL atom parity remains parsed metadata and is no longer used to infer
  tetrahedral configuration.
- MOL double-bond configurations are derived from suitable supplied coordinates.
  Bond stereo codes inconsistent with the bond order now fail parsing.
- StereoModel.stereo_bond_minimum_ring_size defaults to 8. Resolution clears
  cis/trans assertions on ring double bonds below that threshold. A value of 0
  permits realization at every ring size.

## Isotopes and graph projection

- IsotopeResolver handles isotope completion independently of valence resolution.
  ResolveConfig gains isotope: Strict leaves unspecified composition unresolved;
  Natural fills it with natural composition. Explicit masses, sets, and
  variables are preserved.
- General resolution defaults to Strict. SMILES ingestion and export convenience
  functions select Natural explicitly. Valence resolution no longer supplies
  isotope defaults.
- Resolver::project and ProjectFlags expose stereo, aromaticity, valence, and
  isotope projection within graph IR. Only a determined result replaces the
  supplied molecule. Projection preserves implicit hydrogen counts; Convey
  omits a boundary H count only when inference under the selected model and
  policy reproduces it.

## Supplied layouts and SVG rendering

- The Rust Depict trait adds layout/layout_with, verify_layout, and depict_layout
  for Molecule and Reaction. Generated and supplied layouts use the same
  depiction path. Verification checks the atom frame, definite cis/trans
  geometry, and tetrahedral wedge selection without moving coordinates.
- Python exposes MoleculeLayout and ReactionLayout, plus layout, verify_layout,
  and depict_layout on molecules and reactions. Layouts are mutable; use
  set_position, set_lhs_position, set_rhs_position, and set_arrow for edits.
  MoleculeLayout.positions is an immutable snapshot; reaction side getters
  return copies.
- SvgConfig and Depiction.render_svg_with expose the atom-label mask identifier
  and reaction mapping-index text size in Rust and Python. Use distinct mask_id
  values when embedding multiple SVGs in one document.

## Graph traversal

- Add callback-based depth-first and breadth-first traversal, event collectors,
  and a connected-component visitor. Visitors support immediate termination
  through ControlFlow. Connected components use the shared breadth-first path.

## Migration from 0.8.0

- Add isotope to explicit Rust ResolveConfig values. Select IsotopePolicy::Natural
  (Python IsotopePolicy.Natural) when unspecified isotopes should resolve to
  natural composition. Passing a default resolve configuration explicitly to
  SMILES ingestion now selects Strict, unlike omitting that configuration.
- AtomTypeRegistry entries must have undetermined isotopes, literal elements,
  and literal charges in -128..=127. Remove explicit isotope defaults such as
  #i= from custom registry patterns. Rust adds try_from_atoms and try_add;
  asserted from_atoms/add panic on invalid entries. Python from_atoms raises
  ValueError.
- Add stereo_bond_minimum_ring_size to explicit StereoModel construction in
  Rust and Python; use 8 to match StereoModel's default. StereoDerivation.atoms
  and .bonds become .stereo_atoms and .stereo_bonds; it also reports
  skipped_stereo_bonds.
- Replace Rust layout::layout_molecule with Depict::layout or
  Depict::layout_with. Custom Depict implementations now supply Layout,
  layout_with, verify_layout, and depict_layout.
- BondWedge is now a record with orientation: BondOrientation and taper:
  BondTaper. Preserve the pointed endpoint when constructing or renumbering
  wedges.
- Replace TraversalAlgorithm with NeighborhoodAlgorithm. Graph neighborhood
  limits and returned distances, circular-refinement radii and rounds, and
  Rust ECFP/Morgan radii use usize instead of u32. Update custom
  CircularRefinementHash::combine implementations accordingly.
- Update error handling for the revised parsing and conversion boundaries.
  Direction-marker failures formerly reported during graph conversion now
  arise during parsing; Python reports them as ParseError.

umol remains an alpha library. Canonical representatives may change between
0.x releases and are not persistent identifiers.
