//! Ergonomic Rust wrapper for libmsym (molecular symmetry).

mod basis;
pub(crate) mod context;
mod detect;
mod error;
pub(crate) mod irrep;
pub(crate) mod linear;
mod matrix_rep;
mod point_group;
mod subgroup;
mod thresholds;
mod types;

pub use basis::{BasisFunction, BasisKind, CartesianAxis, IrrepBasis, Salc, SalcBasis};
pub use detect::{
    compute_salcs, detect_symmetry, generate_symmetry_images, lower_symmetry, symmetrize,
    SymmetryDescentResult, SymmetryResult,
};
pub use error::{MsymError, ParseError, ReductionError};
pub use irrep::Irrep;
pub use matrix_rep::MatrixRep;
pub use point_group::{
    CharacterTableDisplay, PointGroup, SymmetryOp, SymmetryOpKind, SymmetryOpOrientation,
};
pub use subgroup::{correlation_table, CorrelationTable, Subgroup};
pub use thresholds::Thresholds;
pub use types::{EquivalenceSet, MolecularShape, SchoenfliesSymbol, SymmetryCenter};
