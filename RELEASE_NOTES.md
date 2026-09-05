# umol 0.8.0

Changes since 0.7.0. This release adds SVG depiction, composable mutation tracking,
and a unified parser implementation. It includes breaking Rust and Python API changes.

## SVG depiction

- With the Rust `umol-io/depiction` feature, `Molecule` and `Reaction` implement the
  `Depict` extension trait. `depict()` uses the default `DepictConfig`;
  `depict_with()` accepts an explicit configuration.
- Both return an opaque `Depiction`. Its `render_svg()` method returns a complete
  SVG document suitable for writing directly to a file.
- Python exposes the same methods. `Depiction.render_svg()` returns SVG text,
  and `_repr_svg_()` supplies Jupyter rich display. Published Python distributions
  enable depiction.
- The initial black-and-white renderer supports isotope labels, implicit hydrogen
  subscripts, charge and radical superscripts, ordinary bond orders, tetrahedral
  wedges and hashes, cis/trans geometry, aromatic-system contours, reaction
  arrows, and reaction correspondence indices.
- The new `umol-coordgen-sys` crate vendors CoordGen 3.0.2 for two-dimensional
  layout. Its `native` feature requires a C++11 compiler; Rust builds without
  coordinate generation do not compile the native source.

## Tracked operations and reaction application

- Optional provenance is returned by `tracked_` companions such as `tracked_apply`,
  `tracked_remove`, `tracked_extract`, and `tracked_split`. Bare and tracked
  operations produce the same ordinary result and resulting mutable state.
  Checked companions retain the outer `try_` prefix, as in `try_tracked_build`.
- Editors accumulate initial-to-current correspondences across additions,
  removals, and edits. Tracked snapshot and finalization expose this provenance;
  tracked transactions and rollback expose forward and reverse correspondences.
  Tracking stores ids and counts without cloning molecular payloads.
- Correspondences compose with intermediate-count checks. Compactions now retain
  source counts, and remappings enforce dense bijections. Borrowed conversions
  widen either specialized carrier to correspondence for mixed-operation chains.
- Reaction application returns products, tracked products, realized reactions,
  or realized reaction spans directly. `ReactionDerivation` is removed.
- Python exposes the corresponding high-level carriers and method pairs.

### Migration from 0.7.0

- `combine` and `combine_all` return only the combined molecule; `combine_from`
  returns unit. Their per-entity append order determines the input mappings.
- `split` returns component molecules. `tracked_split` returns each component
  with a **source-to-component** correspondence, reversing the old direction.
- Graph/editor removal, relation compaction, pushout, pullback, and pushout
  complement use bare/tracked pairs. Use the tracked form where the previous
  return included mapping data. Categorical mappings keep their original direction.
- Replace uses of `ReactionDerivation.rhs` with the product returned by `apply`.
  Use `tracked_apply` for `(product, correspondence)`, `apply_to_reaction` for
  realized reactions, or `apply_to_reaction_span` for realized spans. Rust also
  exposes corresponding supplied-match methods; these return
  `Result<Option<T>, ApplyError>`, with `Ok(None)` for an inapplicable match.
- `canonicalize_with_correspondence` becomes `tracked_canonicalize` and
  returns `MoleculeRemapping`. Molecule and reaction-span renumbering and framed
  comparison accept that bijective carrier. Non-bijective reference transport
  uses correspondence through `map`/`try_map`; `IdRemapping` and
  correspondence-to-remapping narrowing are removed.
- `reframe_with_action` becomes `tracked_reframe`; it still returns the reframed
  value and its participant-frame action, not an entity-id remapping.
- Construct `Remapping<Id>` with a permutation of its dense domain;
  construction now returns an error for repeated or out-of-range images.
  Compaction construction takes the source count and removed ids and is checked.
  `empty()` denotes a zero-sized domain; identity requires a declared size.
- Correspondence `compose` and `compose_all` are now checked operations in Rust.
  Python reports incompatible intermediate counts as exceptions.
- Graph-core modules are named `remap` and `compact`; graph-IR exposes compaction
  separately from remapping. Pushout result carriers contain the mappings, with
  the resulting object returned separately.

## Parsing

- CTfile and CXSMILES parsing now use Winnow, alongside EDN parsing; nom is
  removed. The workspace uses Winnow 1.x.
- CTfile's public interface consists of complete MOL/SDF parsing entry points.
  Low-level parser combinators are private. Update direct users of those helpers
  to the complete-input APIs.
- CTfile `ParseError` now reports format-level errors without nom error types.
  Variants and fields have changed, including separate atom/bond index errors
  and zero-based physical byte columns for counts lines and SDF headers.
- Unsupported basic-parser flags return `UnsupportedBasicParseFlags` instead
  of triggering an assertion.

umol remains an alpha library. Canonical representatives may change between
0.x releases and are not persistent identifiers.
