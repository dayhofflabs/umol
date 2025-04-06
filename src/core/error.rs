//! Error types and handling.

use crate::core::serde::FormatVersion;
use crate::core::Capability;
use crate::Element;
use thiserror::Error;
use std::error::Error as StdError;
use std::collections::HashSet;
use crate::core::Model;

/// Errors related to serialization and deserialization
#[derive(Error, Debug)]
pub enum SerializationError {
    #[error("Invalid format: {0}")]
    InvalidFormat(String),

    #[error("Version mismatch: expected {expected}, found {found}")]
    VersionMismatch {
        expected: FormatVersion,
        found: FormatVersion,
    },

    #[error("Missing required field: {0}")]
    MissingField(String),

    #[error("Invalid field value: {0}")]
    InvalidFieldValue(String),

    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),

    #[error("Serialization error: {0}")]
    SerdeError(#[from] serde_json::Error),
}

/// Errors related to validation
#[derive(Error, Debug)]
pub enum ValidationError {
    #[error("Invalid structure: {0}")]
    InvalidStructure(String),

    #[error("Missing required component: {0}")]
    MissingComponent(String),

    #[error("Invalid component: {0}")]
    InvalidComponent(String),

    #[error("Multiple validation errors: {0:?}")]
    Multiple(Vec<ValidationError>),
}

/// Core error types for the molecular modeling framework
#[derive(Debug, Error)]
pub enum Error {
    /// An error occurred during model operations
    #[error(transparent)]
    Model(#[from] ModelError),

    /// An error occurred during property calculations
    #[error(transparent)]
    Property(#[from] PropertyError),

    /// An error occurred during format operations
    #[error(transparent)]
    Format(#[from] FormatError),

    /// An error occurred during conversion operations
    #[error(transparent)]
    Conversion(#[from] ConversionError),

    /// An error occurred during operation execution
    #[error(transparent)]
    Operation(#[from] OperationError),

    /// An error occurred during plugin operations
    #[error(transparent)]
    Plugin(#[from] PluginError),

    /// A dependency is missing
    #[error("Missing dependency: {0}")]
    MissingDependency(String),

    /// A property was not found
    #[error("Property not found: {0}")]
    PropertyNotFound(String),

    /// A model was not found
    #[error("Model not found: {0}")]
    ModelNotFound(String),

    /// A conversion was not found
    #[error("Conversion not found from {0} to {1}")]
    ConversionNotFound(String, String),

    /// Invalid charge value
    #[error("Invalid charge {charge} for element {element}, must be between {min} and {max}")]
    InvalidCharge {
        element: Element,
        charge: i8,
        min: i8,
        max: i8,
    },

    /// Invalid number of unpaired electrons
    #[error("Invalid number of unpaired electrons {unpaired} for element {element}, maximum is {max}")]
    InvalidUnpairedElectrons {
        element: Element,
        unpaired: u8,
        max: u8,
    },

    /// An error occurred during entity operations
    #[error(transparent)]
    Entity(#[from] EntityError),

    /// An error occurred during validation
    #[error(transparent)]
    Validation(#[from] ValidationError),

    /// An error occurred during serialization
    #[error(transparent)]
    Serialization(#[from] SerializationError),

    /// An error occurred during element operations
    #[error(transparent)]
    Element(#[from] ElementError),

    /// Multiple errors occurred
    #[error("Multiple errors occurred: {0:?}")]
    Multiple(Vec<Error>),

    /// Other errors
    #[error(transparent)]
    Other(#[from] Box<dyn StdError + Send + Sync>),
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
    MissingCapability(Capability),

    #[error("Invalid model state: {0}")]
    InvalidState(String),

    #[error("Model operation failed: {0}")]
    OperationFailed(String),

    #[error("Model not found: {0}")]
    NotFound(String),
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

    #[error("No conversion found from {0} to {1}")]
    NotFound(String, String),
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
    MissingCapability(Capability),

    #[error("Invalid property parameters: {0}")]
    InvalidParameters(String),

    #[error("Property calculation failed: {0}")]
    CalculationFailed(String),

    #[error("Property not available: {0}")]
    NotAvailable(String),

    #[error("Property not found: {0}")]
    NotFound(String),
}

/// Errors related to plugins and capabilities
#[derive(Error, Debug)]
pub enum PluginError {
    #[error("Invalid capability format: {0}")]
    InvalidCapability(String),

    #[error("Required plugin not found: {0}")]
    MissingPlugin(String),

    #[error("Plugin {plugin} version mismatch: required {required}, found {found}")]
    VersionMismatch {
        plugin: String,
        required: String,
        found: String,
    },

    #[error("Required capability not found: {0}")]
    MissingCapability(Capability),

    #[error("Component initialization failed: {0}")]
    ComponentInit(String),

    #[error("Required dependency not found: {0}")]
    MissingDependency(String),
}

/// Errors related to chemical elements
#[derive(Error, Debug)]
pub enum ElementError {
    #[error("Invalid charge {charge} for element {element}, must be between {min} and {max}")]
    InvalidCharge {
        element: Element,
        charge: i8,
        min: i8,
        max: i8,
    },

    #[error(
        "Invalid number of unpaired electrons {unpaired} for element {element}, maximum is {max}"
    )]
    InvalidUnpairedElectrons {
        element: Element,
        unpaired: u8,
        max: u8,
    },
}

/// Errors related to format operations
#[derive(Error, Debug)]
pub enum FormatError {
    #[error("Format not found: {0}")]
    NotFound(String),

    #[error("Invalid format: {0}")]
    Invalid(String),

    #[error("Format operation failed: {0}")]
    Failed(String),
}

/// Result type for core operations
pub type Result<T> = std::result::Result<T, Error>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_display() {
        let error = Error::Model(ModelError::NotFound("test".to_string()));
        assert_eq!(format!("{}", error), "Model not found: test");

        let error = Error::Property(PropertyError::CalculationFailed("test error".into()));
        assert_eq!(format!("{}", error), "Property calculation failed: test error");

        let error = Error::Format(FormatError::NotFound("test".to_string()));
        assert_eq!(format!("{}", error), "Format not found: test");

        let error = Error::Conversion(ConversionError::NotFound("source".to_string(), "target".to_string()));
        assert_eq!(format!("{}", error), "No conversion found from source to target");

        let error = Error::Operation(OperationError::Failed("test error".to_string()));
        assert_eq!(format!("{}", error), "Operation failed: test error");

        let error = Error::Plugin(PluginError::MissingPlugin("test".to_string()));
        assert_eq!(format!("{}", error), "Required plugin not found: test");

        let error = Error::Element(ElementError::InvalidCharge {
            element: Element::H,
            charge: 2,
            min: -1,
            max: 1,
        });
        assert_eq!(
            format!("{}", error),
            "Invalid charge 2 for element H, must be between -1 and 1"
        );

        let errors = vec![
            Error::Model(ModelError::NotFound("error 1".to_string())),
            Error::Property(PropertyError::CalculationFailed("error 2".into())),
        ];
        let error = Error::Multiple(errors);
        assert!(format!("{}", error).contains("Multiple errors occurred"));
    }
}

// Update the helper functions to use the appropriate error types
pub fn verify_capabilities(model: &impl Model, required: &[Capability]) -> Result<()> {
    for cap in required {
        if !model.has_capability(cap) {
            return Err(ModelError::MissingCapability(cap.clone()).into());
        }
    }
    Ok(())
}

pub fn test_model_capabilities(model: &impl Model) -> Result<()> {
    let mut caps = HashSet::new();
    caps.insert(Capability::new("core", "has_atoms", 1));

    for cap in &caps {
        if !model.has_capability(cap) {
            return Err(ModelError::MissingCapability(cap.clone()).into());
        }
    }
    Ok(())
}
