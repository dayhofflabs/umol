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

pub use algebra::{Aggregate, Ensemble, Process};
pub use conversion::{ConversionMetadata, ConvertTo, ConvertToWithMetadata};
pub use entity::{Entity, Relation};
pub use error::{Error, Result};
pub use instance::Instance;
pub use io::{FileSystem, FormatDetector, ReadSeek};
pub use model::{Capability, Model};
pub use operation::{ConversionOperation, Operation};
pub use plugin::{
    ConversionCompute, ConversionDefinition, ModelProvider, Plugin, PluginRequirements,
    PropertyCompute, PropertyDefinition, Registry,
};
pub use property::Property;
pub use serde::{FormatVersion, ModelSerializer, SerializableModel};
pub use testing::{ModelTest, PropertyTest};
