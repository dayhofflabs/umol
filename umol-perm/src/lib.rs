//! Pure permutation and coset algebra for umol stereochemistry.
//!
//! The dense coset index this crate computes reproduces the OpenSMILES
//! arrangement number (`@TH`/`@AL`/`@SP`/`@TB`/`@OH`).

mod coset;
mod group;
mod permutation;

pub use coset::CosetSpace;
pub use group::PermutationGroup;
pub use permutation::Permutation;
