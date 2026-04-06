# Geometric transformation types

umol lacks a general geometric transformation type that can be applied to molecular coordinates.

- SE(3) rigid body transformations (rotation + translation)
- O(3) symmetry operations (rotation/reflection, origin-centered) — currently on `SymmetryOp::matrix` with `transform_point`
- Non-rigid transformations for vibrational modes

nalgebra provides `Isometry3` for SE(3). Question is whether to use it directly or wrap it in a umol-geometric type that connects to `Molecule`.

- Should add explicit method to transform molecule into principal axes system.
