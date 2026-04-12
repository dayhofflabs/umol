//! Ergonomic Rust wrapper for libmsym (molecular symmetry).

mod basis;
pub(crate) mod context;
mod detect;
mod error;
pub(crate) mod linear;
mod matrix_rep;
mod point_group;
mod subgroup;
mod types;

pub use basis::{BasisFunction, BasisKind, CartesianAxis, IrrepBasis, Salc, SalcBasis};
pub use detect::{
    compute_salcs, detect_symmetry, generate_symmetry_images, lower_symmetry, symmetrize,
    SymmetryDescentResult, SymmetryResult,
};
pub use error::Error;
pub use matrix_rep::MatrixRep;
pub use point_group::{CharacterTableDisplay, Irrep, PointGroup, ReductionError, SymmetryOp};
pub use subgroup::{correlation_table, CorrelationTable, SubgroupInfo};
pub use types::{
    EquivalenceSet, Geometry, SchoenfliesLabel, SymmetryCenter, SymmetryOpKind,
    SymmetryOpOrientation, Thresholds,
};
