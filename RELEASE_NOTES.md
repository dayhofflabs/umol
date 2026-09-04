# umol 0.7.1

umol 0.7.1 adds compatible SVG depiction APIs for graph-IR molecules and reactions.

## SVG depiction

- With the Rust `umol-io/depiction` feature, `Molecule` and `Reaction` implement the `Depict`
  extension trait. `depict()` uses the default `DepictConfig`; `depict_with()` accepts an explicit
  configuration. CoordGen is the initial and default layout algorithm.
- Both operations return an opaque `Depiction`. `Depiction::render_svg()` returns a complete SVG
  document that can be written directly to a file.
- Python `Molecule` and `Reaction` provide the equivalent `depict()` and `depict_with()` methods.
  Python `Depiction.render_svg()` returns ordinary SVG text, while `_repr_svg_()` provides Jupyter
  rich display. Published Python distributions enable this capability.
- The initial black-and-white renderer depicts ordinary bond orders, atom labels and annotations,
  tetrahedral wedges, cis/trans geometry, explicit aromatic-system contours, reaction arrows, and
  reaction correspondence indices.
