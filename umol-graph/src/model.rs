//! Chemistry-facing types: [`Molecule`] and [`Pattern`].
//!
//! Both wrap a [`MoleculeAst`] with usage-specific invariants and a cache layer
//! shared via `Arc`. `Molecule` enforces the ground invariant at construction
//! and lazily computes derived views (rings, etc.). `Pattern` accepts any
//! `MoleculeAst` and serves as the substructure-query type for the matcher.
//!
//! [`MoleculeAst`]: crate::ast::molecule::MoleculeAst

pub mod molecule;
pub mod pattern;

pub use molecule::Molecule;
pub use pattern::Pattern;
