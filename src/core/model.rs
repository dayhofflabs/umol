//! Core model traits and types.
//! 
//! This module defines the fundamental abstractions for molecular models:
//! - Model trait for representing molecular systems
//! - Capability system for describing model features
//! - Basic model operations and validations

use std::collections::HashSet;
use std::fmt;
use crate::core::Result;

/// A capability that a model can provide
#[derive(Debug, Clone, Hash, Eq, PartialEq)]
pub struct Capability {
    /// The namespace of the capability
    pub namespace: String,
    /// The name of the capability
    pub name: String,
    /// The version of the capability
    pub version: u32,
}

impl Capability {
    /// Create a new capability
    pub fn new(namespace: &str, name: &str, version: u32) -> Self {
        Self {
            namespace: namespace.to_string(),
            name: name.to_string(),
            version,
        }
    }
}

impl fmt::Display for Capability {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}", self.name, self.version)
    }
}

/// A trait for molecular models
pub trait Model {
    /// The type of data stored in this model
    type Data;
    
    /// Get a reference to the model's data
    fn data(&self) -> &Self::Data;
    
    /// Get the capabilities provided by this model
    fn capabilities(&self) -> HashSet<Capability>;
    
    /// Check if this model provides a specific capability
    fn has_capability(&self, capability: &Capability) -> bool {
        self.capabilities().contains(capability)
    }
    
    /// Validate the model
    fn validate(&self) -> Result<()> {
        Ok(())
    }
}

/// Trait for converting between models
pub trait ConvertTo<M: Model> {
    /// Convert this model to another model type
    fn convert_to(&self) -> Result<M>;
}

/// Trait for converting between models with parameters
pub trait ConvertToWithMetadata<M: Model> {
    /// The type of parameters needed for conversion
    type Params;
    
    /// Convert this model to another model type with parameters
    fn convert_to_with_metadata(&self, params: &Self::Params) -> Result<(M, ConversionMetadata)>;
}

/// Metadata about a conversion
#[derive(Debug, Default)]
pub struct ConversionMetadata {
    /// The accuracy of the conversion
    pub accuracy: f64,
    /// The time taken for the conversion
    pub time_taken: std::time::Duration,
} 