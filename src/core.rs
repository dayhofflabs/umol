//! Core module for umol.
//! 
//! This module provides the fundamental abstractions for molecular modeling:
//! - Entity and Relation traits for representing molecular entities
//! - Model trait for representing molecular systems
//! - Property traits for calculating molecular properties
//! - Conversion traits for converting between models
//! - Error types for error handling

pub mod entity;
pub mod model;
pub mod algebra;
pub mod instance;
pub mod property;
pub mod conversion;
pub mod error;
pub mod testing;
pub mod serde;
pub mod io;

// Re-export commonly used types
pub use entity::{Entity, Relation};
pub use model::{Model, Capability};
pub use algebra::{Ensemble, Aggregate};
pub use instance::{Instance, Operation};
pub use property::{Property, MolecularProperty, EnergyProperty, StructuralProperty};
pub use conversion::{ConvertTo, ConvertToWithMetadata, ConversionMetadata};
pub use error::{Error, Result};
pub use testing::{ModelTest, PropertyTest};
pub use serde::{SerializableModel, FormatVersion};
pub use io::ReadSeek; 