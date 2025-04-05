//! Core abstractions for molecular modeling.
//! 
//! This module provides the fundamental building blocks of umol:
//! 
//! # Core Concepts
//! 
//! - **Entities**: Abstract representations of molecular objects
//! - **Models**: Concrete implementations of molecular systems
//! - **Properties**: Calculations that can be performed on models
//! - **Instances**: Combinations of entities and their model representations
//! 
//! # Plugin System
//! 
//! The core module also provides a plugin architecture for extending umol:
//! 
//! - Plugin registration and management
//! - Capability system for feature detection
//! - Lazy loading of components
//! - Model conversion framework
//! 
//! # Error Handling
//! 
//! Comprehensive error types and handling for:
//! 
//! - Plugin operations
//! - Model operations
//! - Property calculations
//! - IO operations

mod entity;
mod model;
mod algebra;
mod instance;
mod property;
mod conversion;
mod error;
mod io;

#[cfg(test)]
mod testing;

// Re-export commonly used types
pub use entity::{Entity, Relation};
pub use model::{Model, Capability};
pub use algebra::{Ensemble, Aggregate};
pub use instance::{Instance, Operation};
pub use property::Property;
pub use conversion::{ConvertTo, ConvertToWithMetadata, ConversionMetadata};
pub use error::{Error, Result}; 