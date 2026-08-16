# umol 0.6.0

umol 0.6.0 is the first public alpha release of the Rust workspace and the
`umol-py` Python distribution. The Python distribution is installed with
`pip install umol-py` and imported as `umol`.

## Included functionality

- A molecular graph IR with typed relational overlays for aromatic systems,
  stereochemistry, coordination, multicenter bonding, and other structure
  outside the localized-bond graph.
- Lattice-valued attributes shared by concrete molecules, structural
  patterns, and reaction rules.
- Parsing and rendering of umol notation, plus SMILES and MOL/SDF input
  boundaries.
- Chemistry-model-driven resolution and validation, including Python's
  non-mutating `Molecule.resolve()` and its three-valued solution; substructure
  matching, canonical comparison, fingerprints, and graph algorithms.
- Reaction construction, application, splitting, combining, and composition.
- Rust libraries and Python bindings for the high-level molecule and reaction
  workflows described in the whitepaper.

## Known limitations

- Molecule-scope constraints on patterns are rejected with an explicit error;
  they are not yet evaluated during matching.
- The default stereo model supports the tetrahedral and cis/trans behavior
  demonstrated in the whitepaper. Higher stereo kinds such as allene,
  square-planar, and octahedral are staged off by default.

## Compatibility and chemistry-data updates

The 0.6 line follows Cargo's pre-1.0 compatibility convention: compatible
changes ship as 0.6.x, while breaking public API changes require 0.7.0.

Conformance-suite growth, atom-typing registry rows, and additions to living
default valence tables are compatible 0.6.x updates. Such additions may change
resolution outcomes by adding candidates; candidate sets are the monotone
observable. Frozen preset tables are never changed in place: a revised preset
receives a new name and at least a minor version bump.
