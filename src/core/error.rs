//! Error types and handling.
//! 
//! This module defines the error types used throughout umol:
//! - Core error types and results
//! - Plugin-related errors
//! - Model and property errors
//! - IO and conversion errors
//! - Error conversion traits

use thiserror::Error;
use std::error::Error as StdError;
use std::fmt;
use semver::Version;
use crate::core::Capability;

/// Core error types for the molecular modeling framework
#[derive(Debug, Error)]
pub enum Error {
    #[error("Entity error: {0}")]
    Entity(EntityError),

    #[error("Model error: {0}")]
    Model(ModelError),

    #[error("Conversion error: {0}")]
    Conversion(ConversionError),

    #[error("Operation error: {0}")]
    Operation(OperationError),

    #[error("Property error: {0}")]
    Property(PropertyError),

    #[error("Validation error: {0}")]
    ValidationError(String),

    #[error("Multiple errors occurred: {0:?}")]
    Multiple(Vec<Error>),

    #[error(transparent)]
    Other(#[from] Box<dyn std::error::Error + Send + Sync>),

    // Plugin-related errors
    #[error("Invalid capability format: {0}")]
    InvalidCapability(String),
    
    #[error("Required plugin not found: {0}")]
    MissingPlugin(String),
    
    #[error("Plugin {plugin} version mismatch: required {required}, found {found}")]
    PluginVersionMismatch {
        plugin: String,
        required: Version,
        found: Version,
    },
    
    #[error("Required capability not found: {0}")]
    MissingCapability(Capability),
    
    // Component-related errors
    #[error("Property not found: {0}")]
    PropertyNotFound(String),
    
    #[error("Model not found: {0}")]
    ModelNotFound(String),
    
    #[error("No conversion found from {0} to {1}")]
    ConversionNotFound(String, String),
    
    #[error("Format not found: {0}")]
    FormatNotFound(String),
    
    // Initialization errors
    #[error("Component initialization failed: {0}")]
    ComponentInitError(#[from] Box<dyn std::error::Error + Send + Sync>),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::InvalidCapability(cap) => 
                write!(f, "Invalid capability format: {}", cap),
            Error::MissingPlugin(name) => 
                write!(f, "Required plugin not found: {}", name),
            Error::PluginVersionMismatch { plugin, required, found } => 
                write!(f, "Plugin {} version mismatch: required {}, found {}", 
                    plugin, required, found),
            Error::MissingCapability(cap) => 
                write!(f, "Required capability not found: {}", cap.to_string()),
            Error::PropertyNotFound(name) => 
                write!(f, "Property not found: {}", name),
            Error::ModelNotFound(name) => 
                write!(f, "Model not found: {}", name),
            Error::ConversionNotFound(from, to) => 
                write!(f, "No conversion found from {} to {}", from, to),
            Error::FormatNotFound(name) => 
                write!(f, "Format not found: {}", name),
            Error::ComponentInitError(e) => 
                write!(f, "Component initialization failed: {}", e),
            Error::Other(e) => 
                write!(f, "Error: {}", e),
        }
    }
}

impl StdError for Error {
    fn source(&self) -> Option<&(dyn StdError + 'static)> {
        match self {
            Error::ComponentInitError(e) => Some(e.as_ref()),
            Error::Other(e) => Some(e.as_ref()),
            _ => None,
        }
    }
}

/// Errors related to chemical entities
#[derive(Error, Debug)]
pub enum EntityError {
    #[error("Invalid entity: {0}")]
    Invalid(String),
    
    #[error("Entity not found: {0}")]
    NotFound(String),
    
    #[error("Invalid relationship between entities: {0}")]
    InvalidRelation(String),
    
    #[error("Entity operation failed: {0}")]
    OperationFailed(String),
}

/// Errors related to models
#[derive(Error, Debug)]
pub enum ModelError {
    #[error("Missing required capability: {0:?}")]
    MissingCapability(crate::core::Capability),
    
    #[error("Invalid model state: {0}")]
    InvalidState(String),
    
    #[error("Model operation failed: {0}")]
    OperationFailed(String),
}

/// Errors related to model conversions
#[derive(Error, Debug)]
pub enum ConversionError {
    #[error("Incompatible models: from {from} to {to}")]
    IncompatibleModels { from: String, to: String },
    
    #[error("Information loss during conversion: {0}")]
    InformationLoss(String),
    
    #[error("Conversion failed: {0}")]
    Failed(String),
}

/// Errors related to operations on instances
#[derive(Error, Debug)]
pub enum OperationError {
    #[error("Invalid operation parameters: {0}")]
    InvalidParameters(String),
    
    #[error("Operation not supported: {0}")]
    NotSupported(String),
    
    #[error("Operation failed: {0}")]
    Failed(String),
}

/// Errors related to property calculations
#[derive(Error, Debug)]
pub enum PropertyError {
    #[error("Missing required capability for property calculation: {0:?}")]
    MissingCapability(crate::core::Capability),
    
    #[error("Invalid property parameters: {0}")]
    InvalidParameters(String),
    
    #[error("Property calculation failed: {0}")]
    CalculationFailed(String),
    
    #[error("Property not available: {0}")]
    NotAvailable(String),
}

/// Result type for core operations
pub type Result<T> = std::result::Result<T, Error>; 