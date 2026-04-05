//! Ergonomic Rust wrapper for libmsym (molecular symmetry).

pub(crate) mod context;
mod detect;
mod error;
mod point_group;
mod types;

pub use detect::{detect_symmetry, symmetrize, SymmetryResult};
pub use error::Error;
pub use point_group::PointGroup;
pub use types::{
    BasisFunction, BasisKind, CharacterTable, EquivalenceSet, Geometry, Irrep, PointGroupKind,
    Salc, SubrepresentationSpace, SymmetryCenter, SymmetryOp, SymmetryOpKind,
    SymmetryOpOrientation, Thresholds,
};
