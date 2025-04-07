//! Core functionality for umol
//!
//! This module provides the fundamental abstractions for chemical modeling:
//! - Entity, Model, and Instance traits for representing chemical systems
//! - Capability and Property types for calculating chemical properties
//! - Basic IO, logging, and error handling utilities

mod algebra;
mod capability;
mod conversion;
mod entity;
pub mod error;
mod instance;
pub mod io;
pub mod logging;
mod model;
mod operation;
mod property;

pub use algebra::{Aggregate, Ensemble, Process};
pub use capability::Capability;
pub use conversion::{ConversionMetadata, ConvertTo, ConvertToWithMetadata};
pub use entity::{AsEntity, Entity, Relation};
pub use error::{Error, Result};
pub use instance::Instance;
pub use io::{FileSystem, FormatDetector, ReadSeek};
pub use model::{AsModel, Model};
pub use operation::{ConversionOperation, Operation};
pub use property::{Property, PropertySpec};
