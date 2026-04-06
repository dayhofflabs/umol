//! Ergonomic Rust wrapper for libmsym (molecular symmetry).

pub(crate) mod context;
mod detect;
mod error;
mod point_group;
mod types;

pub use detect::{detect_symmetry, symmetrize, symmetrize_to, SymmetryResult};
pub use error::Error;
pub use point_group::{Irrep, PointGroup};
pub use types::{
    EquivalenceSet, Geometry, SchoenfliesLabel, SymmetryCenter, SymmetryOp, SymmetryOpKind,
    SymmetryOpOrientation, Thresholds,
};
