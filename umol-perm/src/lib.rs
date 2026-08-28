//! Pure permutation and coset algebra for umol stereochemistry.
//!
//! The dense coset index this crate computes reproduces the OpenSMILES
//! arrangement number (`@TH`/`@AL`/`@SP`/`@TB`/`@OH`).

mod class;
mod coset;
mod error;
mod group;
mod oriented;
mod permutation;

pub use class::{ClassKey, Coset};
pub use coset::CosetSpace;
pub use error::{ParseClassKeyError, PermutationError};
pub use group::PermutationGroup;
pub use oriented::{Orientation, OrientedPermutation, OrientedPermutationGroup};
pub use permutation::{Permutation, MAX_DEGREE};
