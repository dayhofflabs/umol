//! Core functionality for umol
//!
//! This module provides the fundamental abstractions for molecular modeling:
//! - Entity and Model traits for representing molecular systems
//! - Property system for calculating molecular properties
//! - Plugin system for extending functionality
//! - Testing utilities for verifying implementations

mod algebra;
mod capability;
mod conversion;
mod entity;
pub mod error;
mod instance;
pub mod io;
mod model;
mod operation;
mod property;

pub use algebra::{Aggregate, Ensemble, Process};
pub use capability::Capability;
pub use conversion::{ConversionMetadata, ConvertTo, ConvertToWithMetadata};
pub use entity::{Entity, Relation};
pub use error::{Error, Result};
pub use instance::Instance;
pub use io::{FileSystem, FormatDetector, ReadSeek};
pub use model::Model;
pub use operation::{ConversionOperation, Operation};
pub use property::Property;
