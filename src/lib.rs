// Main library exports and documentation

pub mod core;
pub mod element;
pub mod models;
pub mod validation;

pub use element::Element;
pub use core::error::{Error, Result};
pub use validation::ValidationSet;

// Re-export commonly used graph model types
pub use models::graph::{
    Atom,
    Bond,
    Molecule,
    Builder,
    Fragment,
    Query,
    Template,
};
