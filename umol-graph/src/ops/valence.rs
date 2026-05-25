//! Valence-resolution algorithm primitives. The dispatch wrapper
//! `ValenceResolver` lives in `ops/resolver/valence.rs`.

pub mod atom_typing;
pub mod counts;
pub mod counts_new;
pub mod invariants;
pub mod normal_valence;
pub mod registry;
pub mod table;

pub use atom_typing::{AtomTypingError, AtomTypingValenceResolver};
pub use counts::{CountsError, CountsValenceResolver};
pub use counts_new::{CountsNewResolver, ValenceScheme};
pub use invariants::{Mismatch, ValenceInvariants};
pub use normal_valence::{NormalValenceEntry, NormalValenceTable};
pub use registry::AtomTypeRegistry;
pub use table::{ValenceEntry, ValenceTable};
