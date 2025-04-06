// Main library exports and documentation

pub mod core;
pub mod element;
// pub mod models;

pub use element::Element;
pub use core::error::{Error, Result};
