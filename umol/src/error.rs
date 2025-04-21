//! Error types and handling.

use crate::{Capability, Model};
use std::collections::HashSet;
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
    #[error("Invalid element: {0}")]
    InvalidElement(String),

    #[error("Invalid bond order: {0}")]
    InvalidBondOrder(String),

    #[error("Invalid bond donation: {0}")]
    InvalidBondDonation(String),

    #[error("Invalid covalent bond: {0}")]
    InvalidCovalentBond(String),

    #[error("Invalid charge {charge} for element {symbol}, must be between {min_charge} and {max_charge}")]
    InvalidCharge {
        symbol: String,
        charge: i8,
        min_charge: i8,
        max_charge: i8,
    },
    #[error(
        "Invalid number of unpaired electrons {unpaired_electrons} for element {symbol}, maximum is {max_unpaired_electrons}"
    )]
    InvalidUnpairedElectrons {
        symbol: String,
        unpaired_electrons: u8,
        max_unpaired_electrons: u8,
    },

    #[error("Invalid number of implicit hydrogens {implicit_hydrogens} for element {symbol}, maximum is {max_implicit_hydrogens}")]
    InvalidImplicitHydrogens {
        symbol: String,
        implicit_hydrogens: u8,
        max_implicit_hydrogens: u8,
    },

    #[error("Invalid valence {valence} for element {symbol}, maximum is {max_valence}")]
    InvalidValence {
        symbol: String,
        valence: u8,
        max_valence: u8,
    },

    #[error("Invalid occupation: {0}")]
    InvalidOccupation(String),

    #[error("Invalid spin state: {}{}{}", 
        unpaired_electrons.map(|electrons| format!("{} unpaired electrons", electrons)).unwrap_or_default(),
        multiplet_name.clone().unwrap_or_default(),
        multiplicity.map(|multiplicity| format!("{} multiplicity", multiplicity)).unwrap_or_default(),
    )]
    InvalidSpinState {
        unpaired_electrons: Option<u8>,
        multiplet_name: Option<String>,
        multiplicity: Option<u8>,
    },

    #[error("Invalid valence state: {0}")]
    InvalidValenceState(String),

    #[error("Duplicate original atom index found: {0}")]
    DuplicateAtomIndex(usize),

    #[error("Bond references original atom index {0} which was not provided")]
    MissingAtomIndex(usize),

    #[error("Invalid valence atom: {0}")]
    InvalidValenceAtom(String),

    #[error("Invalid valence bond: {0}")]
    InvalidValenceBond(String),
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
            min_charge: -1,
            max_charge: 1,
        };
        assert_eq!(
            format!("{}", error),
            "Invalid charge 2 for element H, must be between -1 and 1"
        );
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
