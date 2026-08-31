# umol

umol is an experimental molecular structure and manipulation library. The umol
graph model extends the ordinary molecular graph with typed relational overlays for
features that do not fit into a graph of atoms and localized bonds, including
aromatic systems, stereochemistry, coordination, and multicenter bonds. Its
attributes form a lattice, allowing one vocabulary to represent concrete
molecules, structural patterns, and reactions.

umol is implemented as a Rust library with Python bindings. The key Python
types are `Molecule` and `Reaction`: a molecule pattern is also a `Molecule`,
and a reaction rule is also a `Reaction`. Structures can be read from SMILES,
MOL and SDF files, and umol's own notation.

## Getting started

Install the Python distribution `umol-py`; it provides the import package
`umol`:

```console
pip install umol-py
```

Read a molecule from SMILES:

```python
from umol import Molecule

mol = Molecule.from_smiles("CCO")  # ethanol: 3 atoms, 2 bonds
```

Version 0.7.0 is an alpha release, and its interfaces may still change. Canonical
representatives are not persistent identifiers and may change between 0.x releases. The
[whitepaper](https://github.com/dayhofflabs/umol/blob/main/docs/umol-whitepaper.pdf)
introduces the model and gives a fuller Python and Rust primer. See the
[0.7.0 release notes](https://github.com/dayhofflabs/umol/blob/main/RELEASE_NOTES.md)
for the included functionality and known limitations.

## License

Except where a crate states otherwise, umol is available under either the
[Apache License 2.0](https://github.com/dayhofflabs/umol/blob/main/LICENSE-APACHE)
or the [MIT license](https://github.com/dayhofflabs/umol/blob/main/LICENSE-MIT),
at your option. The sources and licenses of vendored third-party code are
recorded in the
[third-party provenance document](https://github.com/dayhofflabs/umol/blob/main/docs/third-party-code.md).
