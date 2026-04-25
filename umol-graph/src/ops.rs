//! Engines that operate on a `MoleculeAst`: resolvers (valence, aromaticity)
//! and tier-2 validators (electron-count, spin-coupling). Each engine has
//! the shape `Engine::new(&model).op(ast)` returning either
//! `Solution<T, Contradiction>` for chemistry outcomes or
//! `Err(Error)` for setup-level failures.

pub mod aromaticity;
pub mod chemistry;
pub mod error;
pub mod propagate;
pub mod resolve;
pub mod solution;
pub mod validate;
pub mod valence;

// Disabled during the engine/config cleanup:
// - `evaluate` evaluated constraints against the disabled `api::Molecule`
//   wrapper; constraint evaluation will return as a separate pass against
//   `MoleculeAst` directly.
// - `matcher` is out of scope this iteration and depends on the disabled
//   `api` types.
// pub mod evaluate;
// pub mod matcher;
