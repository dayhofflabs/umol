//! Engines that operate on a `Molecule`: resolvers (valence, aromaticity)
//! and validators (electron-count, spin-coupling, constraints, entity
//! structure). Each engine has the shape `Engine::new(&model).op(molecule)`
//! returning either `Solution<T, Contradiction>` for chemistry outcomes or
//! `Err(Error)` for setup-level failures.

pub mod aromaticity;
pub mod canonicalize;
pub mod invariant;
pub mod model;
pub mod resolve;
pub mod stereo;
pub mod transform;
pub mod valence;
pub mod validate;
