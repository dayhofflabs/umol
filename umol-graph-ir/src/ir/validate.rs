//! Transitional reaction checks retained by graph-IR operations.

pub mod dpo;
pub mod reaction;

pub use dpo::{check_reaction_dpo, DpoContradiction};
pub use reaction::ReactionIntegrityError;
