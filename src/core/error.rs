//! Error types and handling.

use crate::core::serde::FormatVersion;
use crate::core::Capability;
use crate::Element;
use thiserror::Error;

/// Core error types for the molecular modeling framework
#[derive(Debug, Error)]
pub enum Error {
    #[error(transparent)]
    Entity(#[from] EntityError),

    #[error(transparent)]
    Model(#[from] ModelError),

    #[error(transparent)]
    Conversion(#[from] ConversionError),

    #[error(transparent)]
    Operation(#[from] OperationError),

    #[error(transparent)]
    Property(#[from] PropertyError),

    #[error(transparent)]
    Plugin(#[from] PluginError),

    #[error(transparent)]
    Element(#[from] ElementError),

    #[error("Validation error: {0}")]
    Validation(String),

    #[error("Multiple errors occurred: {0:?}")]
    Multiple(Vec<Error>),

    #[error(transparent)]
    Other(#[from] Box<dyn std::error::Error + Send + Sync>),

    #[error(transparent)]
    Serialization(#[from] SerializationError),
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

/// Errors related to serialization
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

/// Result type for core operations
pub type Result<T> = std::result::Result<T, Error>;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::Capability;

    #[test]
    fn test_entity_error_display() {
        let error = EntityError::NotFound("benzene".into());
        assert_eq!(format!("{}", error), "Entity not found: benzene");

        let error = EntityError::Invalid("invalid structure".into());
        assert_eq!(format!("{}", error), "Invalid entity: invalid structure");
    }

    #[test]
    fn test_model_error_display() {
        let error = ModelError::MissingCapability(Capability::new(
            "energy",
            "1.0.0",
            "Energy calculation capability",
        ));
        assert_eq!(format!("{}", error), "Missing required capability: energy");

        let error = ModelError::InvalidState("invalid coordinates".into());
        assert_eq!(
            format!("{}", error),
            "Invalid model state: invalid coordinates"
        );
    }

    #[test]
    fn test_conversion_error_display() {
        let error = ConversionError::IncompatibleModels {
            from: "MMFF94".into(),
            to: "UFF".into(),
        };
        assert_eq!(
            format!("{}", error),
            "Incompatible models: from MMFF94 to UFF"
        );

        let error = ConversionError::InformationLoss("stereochemistry".into());
        assert_eq!(
            format!("{}", error),
            "Information loss during conversion: stereochemistry"
        );
    }

    #[test]
    fn test_error_conversion() {
        // Test conversion from EntityError to Error
        let entity_error = EntityError::NotFound("benzene".into());
        let error: Error = entity_error.into();
        assert_eq!(format!("{}", error), "Entity not found: benzene");

        // Test conversion from ModelError to Error
        let model_error = ModelError::MissingCapability(Capability::new(
            "energy",
            "1.0.0",
            "Energy calculation capability",
        ));
        let error: Error = model_error.into();
        assert_eq!(format!("{}", error), "Missing required capability: energy");

        // Test conversion from ConversionError to Error
        let conversion_error = ConversionError::IncompatibleModels {
            from: "MMFF94".into(),
            to: "UFF".into(),
        };
        let error: Error = conversion_error.into();
        assert_eq!(
            format!("{}", error),
            "Incompatible models: from MMFF94 to UFF"
        );
    }

    #[test]
    fn test_validation_error() {
        let error = Error::Validation("Invalid structure".into());
        assert_eq!(format!("{}", error), "Validation error: Invalid structure");
    }

    #[test]
    fn test_multiple_errors() {
        let errors = vec![
            Error::Entity(EntityError::NotFound("benzene".into())),
            Error::Model(ModelError::MissingCapability(Capability::new(
                "energy",
                "1.0.0",
                "Energy calculation capability",
            ))),
        ];
        let error = Error::Multiple(errors);
        assert!(format!("{}", error).contains("Entity not found: benzene"));
        assert!(format!("{}", error).contains("Missing required capability: energy"));
    }

    #[test]
    fn test_operation_error() {
        let error = OperationError::InvalidParameters("Invalid coordinates".into());
        assert_eq!(
            format!("{}", error),
            "Invalid operation parameters: Invalid coordinates"
        );

        let error = OperationError::NotSupported("Optimization".into());
        assert_eq!(
            format!("{}", error),
            "Operation not supported: Optimization"
        );
    }

    #[test]
    fn test_property_error() {
        let error = PropertyError::MissingCapability(Capability::new(
            "energy",
            "1.0.0",
            "Energy calculation capability",
        ));
        assert_eq!(
            format!("{}", error),
            "Missing required capability for property calculation: energy"
        );

        let error = PropertyError::CalculationFailed("Convergence failed".into());
        assert_eq!(
            format!("{}", error),
            "Property calculation failed: Convergence failed"
        );
    }

    #[test]
    fn test_plugin_error() {
        let error = PluginError::MissingPlugin("forcefield".into());
        assert_eq!(
            format!("{}", error),
            "Required plugin not found: forcefield"
        );

        let error = PluginError::VersionMismatch {
            plugin: "forcefield".into(),
            required: "1.0.0".into(),
            found: "0.9.0".into(),
        };
        assert_eq!(
            format!("{}", error),
            "Plugin forcefield version mismatch: required 1.0.0, found 0.9.0"
        );
    }

    #[test]
    fn test_element_error() {
        let error = ElementError::InvalidCharge {
            element: Element::C,
            charge: 5,
            min: -4,
            max: 4,
        };
        assert_eq!(
            format!("{}", error),
            "Invalid charge 5 for element C, must be between -4 and 4"
        );

        let error = ElementError::InvalidUnpairedElectrons {
            element: Element::O,
            unpaired: 3,
            max: 2,
        };
        assert_eq!(
            format!("{}", error),
            "Invalid number of unpaired electrons 3 for element O, maximum is 2"
        );
    }
}
