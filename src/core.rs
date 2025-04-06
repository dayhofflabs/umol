//! Core functionality for umol
//! 
//! This module provides the fundamental abstractions for molecular modeling:
//! - Entity and Model traits for representing molecular systems
//! - Property system for calculating molecular properties
//! - Plugin system for extending functionality
//! - Testing utilities for verifying implementations

pub mod algebra;
pub mod conversion;
pub mod entity;
pub mod error;
pub mod examples;
pub mod instance;
pub mod io;
pub mod model;
pub mod operation;
pub mod plugin;
pub mod property;
pub mod serde;
pub mod testing;
#[cfg(test)]
pub mod tests;

pub use error::{Error, Result};
pub use model::Capability;
pub use property::Property;
pub use entity::Entity;
pub use model::Model;
pub use instance::Instance;
pub use conversion::{ConvertTo, ConvertToWithMetadata, ConversionMetadata};

// Re-export plugin-related types
pub use plugin::{
    Plugin, Registry, PluginRequirements,
    ModelProvider, PropertyProvider, ConversionProvider,
    FormatProvider, OntologyProvider,
}; 