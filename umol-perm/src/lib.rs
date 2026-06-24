//! Pure permutation and coset algebra for umol stereochemistry.
//!
//! The dense coset index this crate computes reproduces the OpenSMILES
//! arrangement number (`@TH`/`@AL`/`@SP`/`@TB`/`@OH`).

mod class;
mod coset;
mod group;
mod oriented;
mod permutation;

pub use class::{space, ClassKey, Coset};
pub use coset::CosetSpace;
pub use group::PermutationGroup;
pub use oriented::{Orientation, OrientedPermutation, OrientedPermutationGroup};
pub use permutation::Permutation;
