//! Ergonomic Rust wrapper for libmsym (molecular symmetry).

mod basis;
pub(crate) mod context;
mod detect;
mod error;
pub(crate) mod linear;
mod point_group;
mod types;

pub use basis::{BasisFunction, BasisKind, CartesianAxis, IrrepBasis, Salc, SalcBasis};
pub use detect::{
    compute_salcs, detect_symmetry, generate_symmetry_images, symmetrize, SymmetryResult,
};
pub use error::Error;
pub use point_group::{CharacterTableDisplay, Irrep, PointGroup};
pub use types::{
    EquivalenceSet, Geometry, SchoenfliesLabel, SymmetryCenter, SymmetryOp, SymmetryOpKind,
    SymmetryOpOrientation, Thresholds,
};
