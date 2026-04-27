//! Engines that operate on a `MoleculeAst`: resolvers (valence, aromaticity)
//! and tier-2 validators (electron-count, spin-coupling, constraints, entity
//! structure). Each engine has the shape `Engine::new(&model).op(ast)`
//! returning either `Solution<T, Contradiction>` for chemistry outcomes or
//! `Err(Error)` for setup-level failures.

pub mod aromaticity;
pub mod config;
pub mod solution;
pub mod validator;
pub mod valence;

// Disabled during the engine/config restructure (doc 92). Files kept on disk
// while their content is migrated, phase by phase:
// - `chemistry`, `error` -> `config.rs`, per-engine error types (phase 1+)
// - `propagate` -> `validator.rs`, `valence/*`, `valence.rs` (phase 3, 5)
// - `aromaticity`, `aromaticity/*` -> `aromaticity.rs` + subdir (phase 4)
// - `resolve`, `validate` -> `resolver.rs`, `validator.rs` (phase 3, 6)
// - `evaluate`, `matcher` were already disabled; remain so until the matcher
//   work picks up.
// pub mod aromaticity;
// pub mod chemistry;
// pub mod error;
// pub mod propagate;
// pub mod resolve;
// pub mod validate;
// pub mod evaluate;
// pub mod matcher;
