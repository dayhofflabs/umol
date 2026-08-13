//! Valence algorithm primitives. The dispatch wrappers `ValenceResolver` and
//! `ValenceConformanceValidator` live in `ops/resolve/valence.rs` and
//! `ops/validate/valence.rs`.

pub mod atom_typing;
pub mod compare;
pub mod completion;
pub mod counts;
pub mod registry;
pub mod table;

pub use atom_typing::{AtomTypingError, AtomTypingMismatch, AtomTypingValence};
pub use completion::AtomCompletions;
pub use counts::{CountsError, CountsMismatch, CountsValence};
pub use registry::AtomTypeRegistry;
pub use table::{ValenceEntry, ValenceTable};
