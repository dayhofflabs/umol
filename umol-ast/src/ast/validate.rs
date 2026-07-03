//! Tier-1 (data-integrity) and tier-2 (invariant) validators for molecule and reaction ASTs.
//!  
//! - [`entity`] (data-integrity): structural shape checks on molecule ASTs.
//! - [`constraint`] (data-integrity): cross-entity / molecule-scope constraint evaluation.
//! - [`dpo`] (invariant): DPO reaction invariant (dangling-freedom).

pub mod constraint;
pub mod dpo;
pub mod entity;

pub use constraint::{ConstraintContradiction, ConstraintError, ConstraintValidator};
pub use dpo::{DpoContradiction, DpoError, DpoValidator};
pub use entity::{EntityStructureContradiction, EntityStructureError, EntityStructureValidator};
