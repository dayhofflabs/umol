// Main library exports and documentation

pub mod element;
pub mod error;
pub mod graph;
pub mod io;
pub mod traits;

pub use element::Element;
pub use error::Error;
pub use traits::{AtomLink, AtomSite};
