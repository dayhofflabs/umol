//! Valence-resolution algorithm primitives. The dispatch wrapper
//! `ValenceResolver` lives in `ops/resolver/valence.rs`.

pub mod atom_typing;
pub mod compare;
pub mod counts;
pub mod registry;
pub mod table;

pub use atom_typing::{AtomTypingError, AtomTypingValence};
pub use counts::{CountsError, CountsValence};
pub use registry::AtomTypeRegistry;
pub use table::{ValenceEntry, ValenceTable};
