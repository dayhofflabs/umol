# umol 0.11.0

Changes since 0.10.0. This release revises stereo constraints and graph-derived
readings, adds stereo transformations, expands valence data, and reworks molecular
and reaction-span canonicalization and equality. It includes breaking Rust,
Python, DSL, and custom atom-type registry changes.

## Stereo representation and readings

- Separate site-kind constraints from configurations. Atom #T and bond #C now
  constrain tetrahedral or cis/trans stereo kind only, with undetermined,
  negative and positive forms. Configurations belong to explicit stereo atoms
  and bonds, whose ordered ligand frames give them meaning.
- Add graph-derived chirality, stereogenicity, ligand-symmetry and topicity
  readings for stereo entities. Stored arrangements remain representable even
  when nonstereogenic or prostereogenic; an entity's presence does not assert
  stereogenicity. Topicity distinguishes constitutionally different ligands
  from stereotopic relationships.
- Distinguish stereogenicity in the current configuration from a site's role
  across possible configurations of the molecule. StereogenicSiteKind identifies
  sites whose stereoisomerism follows from constitution alone (Constitutional),
  sites that become stereogenic for some configurations of the other sites
  (Configurational), and Prostereogenic or Nonstereogenic sites. One query returns
  the classifications for all stereo entities, sorted by site, omitting entities
  whose kind is undetermined.
- Add an explicit Virtual ligand kind for a frame position whose meaning is
  unspecified. Resolution identifies it as implicit hydrogen or a lone pair
  only when the identification is unique, preserving its frame position.
- Remove fluxionality constraints from stereo atoms and bonds, including their
  DSL and Python forms. Fluxional behavior requires a description above a
  single molecular graph.

## Stereo operations and graph symmetry

- Add StereoAdder, StereoRemover and StereoEnumerator, with corresponding
  molecule methods in Rust and Python. Addition creates missing tetrahedral
  and cis/trans entities and selects Constitutional and Configurational sites
  by default. Removal clears selected entities and their corresponding site
  assertions. Enumeration assigns configurations to existing entities in a
  supplied site order and retains symmetry-equivalent assignments.
- Expose reusable GraphSymmetry auxiliary data and molecule-bound query views.
  The data contains no molecule copy; queries refer to the supplied unchanged
  molecule. Rust and Python can reuse it for orbit and stereo readings.

## Canonicalization and equality

- Improve performance of molecule canonicalization and canonical equality.
- Use description levels in reaction-span canonicalization and equality.
- Correct normalized, framed and canonical equality for intrinsically
  contradictory aggregates.

## Valence data and registries

- Expand the built-in atom-type registry across alkali and alkaline-earth
  elements, transition metals, lanthanides and actinides, and add reactive
  organic intermediates and associated conformance cases.
- Index custom registries by element and localized valence rather than charge.
  Registry construction checks literal elements, charges and localized valences,
  agreement with TOML section keys, and undetermined isotope state. Charge
  remains part of each atom pattern.

## Migration from 0.10.0

- Replace configuration-valued #T/#C assertions with kind-only constraints.
  Put configuration values on explicit stereo entities together with their
  ligand frames. Update callers for the revised stereo-form constructors and
  remove fluxionality constraints.
- Replace Stereogenicity.Symmetric with Nonstereogenic and Prochiral with
  Prostereogenic. Use StereogenicSiteKind when selecting sites across possible
  configurations; stereogenicity in the current configuration is a different
  property.
- Compute GraphSymmetry with its configuration and use the molecule-bound view,
  or retain the auxiliary data for queries on the unchanged molecule. Recompute
  it after changing that molecule.
- Regroup TOML registry entries by localized valence. In Rust and Python, use
  lookup(element, valence) or lookup_by_element(element) for registry queries.
  Rust valence errors are exposed through ops::valence.
- Replace removed canonical-key calls with the canonicalization, canonical
  equality or canonical hashing operation appropriate to the intended reading.
  Remove the obsolete CanonicalizeContext.para_stereo field.

umol remains an alpha library. Canonical representatives may change between
0.x releases and are not persistent identifiers.
