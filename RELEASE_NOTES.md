# umol 0.6.0

umol 0.6.0 is the first public alpha release of the Rust workspace and the
`umol-py` Python distribution. The Python distribution is installed with
`pip install umol-py` and imported as `umol`.

## Features

- Molecular graph IR with typed relational overlays for aromatic systems,
  stereochemistry, coordination, multicenter bonding, and other structure
  outside the localized-bond graph.
- Lattice-valued attributes shared by concrete molecules, structural
  patterns, and reaction rules.
- Parsing and rendering of umol notation, plus SMILES and MOL/SDF ingestion.
- Chemistry-model-driven resolution and validation, substructure
  matching, canonical comparison, fingerprints, and graph algorithms.
- Reaction IR based on molecular graphs, including application and composition.
- Rust libraries and Python bindings for the high-level molecule and reaction
  workflows described in the whitepaper.

## Compatibility and chemistry-data updates

The 0.6 line follows Cargo's pre-1.0 compatibility convention: 0.6.x releases
maintain compatibility, while breaking public API changes are included 0.7.0
and higher releases.

Conformance-suite growth, atom-typing registry rows, and additions to living
default valence tables are compatible 0.6.x updates. Such additions may change
resolution outcomes by adding candidates; candidate sets are the monotone
observable. Frozen preset tables are never changed in place: a revised preset
receives a new name and at least a minor version bump.
