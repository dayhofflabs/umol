//! Permutation and coset algebra for umol.
//!
//! [`DynPermutation`] supports arbitrary-degree actions, while the dense coset index over bounded
//! [`Permutation`] values reproduces the OpenSMILES arrangement number
//! (`@TH`/`@AL`/`@SP`/`@TB`/`@OH`).

mod class;
mod coset;
mod dynamic;
mod error;
mod group;
mod oriented;
mod permutation;

pub use class::ClassKey;
pub use coset::{Coset, CosetSpace};
pub use dynamic::DynPermutation;
pub use error::{ParseClassKeyError, PermutationError};
pub use group::PermutationGroup;
pub use oriented::{Orientation, OrientedPermutation, OrientedPermutationGroup};
pub use permutation::{Permutation, MAX_DEGREE};
