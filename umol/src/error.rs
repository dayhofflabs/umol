//! Error types and handling.

use crate::{Capability, Model};
use std::collections::HashSet;
use std::error::Error as StdError;
use thiserror::Error;

/// umol error types
#[derive(Debug, Error)]
pub enum Error {
    /// Model operations errors
    #[error(transparent)]
    Model(#[from] ModelError),

    /// Property calculations errors
    #[error(transparent)]
    Property(#[from] PropertyError),

    /// Format operations errors
    #[error(transparent)]
    Format(#[from] FormatError),

    /// Conversion operations errors
    #[error(transparent)]
    Conversion(#[from] ConversionError),

    /// Operation execution errors
    #[error(transparent)]
    Operation(#[from] OperationError),

    /// Entity operations errors
    #[error(transparent)]
    Entity(#[from] EntityError),

    /// Validation errors
    #[error(transparent)]
    Validation(#[from] ValidationError),

    /// Serialization errors
    #[error(transparent)]
    Serialization(#[from] SerializationError),

    /// Data processing errors
    #[error(transparent)]
    Data(#[from] DataError),

    /// Multiple errors occurred
    #[error("Multiple errors occurred: {0:?}")]
    Multiple(Vec<Error>),

    /// Other errors
    #[error(transparent)]
    Other(#[from] Box<dyn StdError + Send + Sync>),
}

/// Errors related to serialization and deserialization
#[derive(Error, Debug)]
pub enum SerializationError {
    #[error("Parse error: {0}")]
    ParseError(String),

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

/// Errors related to chemical entities
#[derive(Error, Debug)]
pub enum EntityError {
    #[error("Invalid entity: {0}")]
    Invalid(String),

    #[error("Entity not found: {0}")]
    NotFound(String),

    #[error("Invalid relationship between entities: {0}")]
    InvalidRelation(String),
}

/// Errors related to models
#[derive(Error, Debug)]
pub enum ModelError {
    #[error("Missing required capability: {0:?}")]
    MissingCapability(Capability),

    #[error("Invalid model state: {0}")]
    InvalidState(String),

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

/// Errors related to operations on stuff
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

/// Errors related to data validation and operations
#[derive(Error, Debug)]
pub enum DataError {
    #[error("Invalid charge {charge} for element {symbol}, must be between {min} and {max}")]
    InvalidCharge {
        symbol: String,
        charge: i8,
        min: i8,
        max: i8,
    },
    #[error(
        "Invalid number of unpaired electrons {unpaired} for element {symbol}, maximum is {max}"
    )]
    InvalidUnpairedElectrons {
        symbol: String,
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

/// umol result type
pub type Result<T> = std::result::Result<T, Error>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_display() {
        let error = ModelError::NotFound("test".to_string());
        assert_eq!(format!("{}", error), "Model not found: test");

        let error = PropertyError::CalculationFailed("test error".into());
        assert_eq!(
            format!("{}", error),
            "Property calculation failed: test error"
        );

        let error = FormatError::NotFound("test".to_string());
        assert_eq!(format!("{}", error), "Format not found: test");

        let error = ConversionError::NotFound("source".to_string(), "target".to_string());
        assert_eq!(
            format!("{}", error),
            "No conversion found from source to target"
        );

        let error = OperationError::Failed("test error".to_string());
        assert_eq!(format!("{}", error), "Operation failed: test error");

        let error = DataError::InvalidCharge {
            symbol: "H".to_string(),
            charge: 2,
            min: -1,
            max: 1,
        };
        assert_eq!(
            format!("{}", error),
            "Invalid charge 2 for element H, must be between -1 and 1"
        );

        let errors: Vec<Error> = vec![
            ModelError::NotFound("error 1".to_string()).into(),
            PropertyError::CalculationFailed("error 2".into()).into(),
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
    caps.insert(Capability::local("has_atoms", 1));

    for cap in &caps {
        if !model.has_capability(cap) {
            return Err(ModelError::MissingCapability(cap.clone()).into());
        }
    }
    Ok(())
}
