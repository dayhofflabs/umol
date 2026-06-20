//! Engines that operate on a `MoleculeAst`: resolvers (valence, aromaticity)
//! and validators (electron-count, spin-coupling, constraints, entity
//! structure). Each engine has the shape `Engine::new(&model).op(ast)`
//! returning either `Solution<T, Contradiction>` for chemistry outcomes or
//! `Err(Error)` for setup-level failures.

pub mod aromaticity;
pub mod invariant;
pub mod model;
pub mod resolve;
pub mod transform;
pub mod valence;
pub mod validate;
