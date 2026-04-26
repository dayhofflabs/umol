//! Valence resolver and its supporting data.
//!
//! Phase 1: only `registry` and `table` are wired. The actual resolver
//! variants (`AtomTypingValenceResolver`, `CountsValenceResolver`) and the
//! dispatching `ValenceResolver` enum land in phase 5.

pub mod registry;
pub mod table;

pub use registry::AtomTypeRegistry;
pub use table::{NormalValenceTable, ValenceEntry, ValenceTable};
