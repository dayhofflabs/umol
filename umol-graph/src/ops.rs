//! Engines that operate on a `MoleculeAst`: resolvers (valence, aromaticity)
//! and tier-2 validators (electron-count, spin-coupling, constraints, entity
//! structure). Each engine has the shape `Engine::new(&model).op(ast)`
//! returning either `Solution<T, Contradiction>` for chemistry outcomes or
//! `Err(Error)` for setup-level failures.

pub mod aromaticity;
pub mod config;
pub mod resolver;
pub mod solution;
pub mod valence;
pub mod validator;

// `evaluate` and `matcher` stay disabled until the matcher work picks up.
// pub mod evaluate;
// pub mod matcher;
