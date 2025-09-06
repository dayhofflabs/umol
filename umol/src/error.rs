//! Error types and handling.

use thiserror::Error;

use crate::Capability;

/// umol error types
#[derive(Debug, Error)]
pub enum Error {
    /// Model operations errors
    #[error(transparent)]
    Model(#[from] ModelError),

    /// Property calculations errors
    #[error(transparent)]
    Property(#[from] PropertyError),

    /// Parsing errors
    #[error(transparent)]
    Parse(#[from] ParseError),

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

    #[error("Invalid isotope: {0}")]
    InvalidIsotope(String),

    #[error("Invalid occupation: {0}")]
    InvalidOccupation(String),

    #[error("Invalid spin state: {0}")]
    InvalidSpinState(String),

    #[error("Invalid atom: {0}")]
    InvalidAtom(String),

    #[error("Invalid atom spec: {0}")]
    InvalidAtomSpec(String),

    #[error("Invalid atom charge: {0}")]
    InvalidAtomCharge(String),

    #[error("Invalid atom lone pair specification: {0}")]
    InvalidAtomLonePairs(String),

    #[error("Invalid atom donated pair specification: {0}")]
    InvalidAtomDonatedPairs(String),

    #[error("Invalid atom accepted pair specification: {0}")]
    InvalidAtomAcceptedPairs(String),

    #[error("Invalid atom unpaired electron specification: {0}")]
    InvalidAtomUnpairedElectrons(String),

    #[error("Invalid atom spin multiplicity: {0}")]
    InvalidAtomMultiplicity(String),

    #[error("Invalid atom implicit hydrogen specification: {0}")]
    InvalidAtomImplicitHydrogens(String),

    #[error("Invalid atom valence: {0}")]
    InvalidAtomValence(String),

    #[error("No matching atom spec found: {0}")]
    NoAtomSpec(String),

    #[error("Multiple matching atom specs found: {0}")]
    MultipleAtomSpecs(String),

    #[error("Invalid bond: {0}")]
    InvalidBond(String),

    #[error("Invalid bond order: {0}")]
    InvalidBondOrder(String),

    #[error("Invalid bond donation: {0}")]
    InvalidBondDonation(String),

    #[error("Invalid bond spec: {0}")]
    InvalidBondSpec(String),

    #[error("No matching bond spec found: {0}")]
    NoBondSpec(String),

    #[error("Multiple matching bond specs found: {0}")]
    MultipleBondSpecs(String),

    #[error("Missing atom index: {0}")]
    MissingAtomIndex(usize),

    #[error("Duplicate atom index: {0}")]
    DuplicateAtomIndex(usize),

    #[error("Missing bond index: {0}")]
    MissingBondIndex(usize),

    #[error("Duplicate bond index: ({0}, {1})")]
    DuplicateBondIndex(usize, usize),

    #[error("Loop bond: {0}")]
    LoopBond(usize),

    #[error("Missing property {0} for atom {1}")]
    MissingAtomProperty(String, usize),

    #[error("Missing property {0} for bond {1}")]
    MissingBondProperty(String, usize),

    #[error("Invalid conformer definition: {0}")]
    InvalidConformer(String),

    #[error("Invalid fragment: {0}")]
    InvalidFragment(String),

    #[error("Invalid feature: {0}")]
    InvalidFeature(String),

    #[error("Invalid MOL file format: {0}")]
    InvalidMolFormat(String),

    #[error("Invalid SDF file format: {0}")]
    InvalidSdfFormat(String),
}

/// Parsing errors
#[derive(Error, Debug)]
pub enum ParseError {
    #[error("Unknown format: {0}")]
    UnknownFormat(String),

    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),

    #[error("Parse failed at line {0} ('{1}'): {2}")]
    Failed(usize, String, String),

    #[error("Invalid: {0}")]
    Invalid(String),

    #[error("Incomplete: {0}")]
    Incomplete(String),
}

/// umol result type
pub type Result<T> = std::result::Result<T, Error>;

#[cfg(test)]
mod tests {
    use rstest::*;

    use super::*;

    #[rstest]
    #[case(ModelError::NotFound("test".to_string()).into(), "Model not found: test")]
    #[case(PropertyError::CalculationFailed("test error".to_string()).into(), "Property calculation failed: test error")]
    #[case(ParseError::UnknownFormat("test".to_string()).into(), "Unknown format: test")]
    #[case(ConversionError::NotFound("source".to_string(), "target".to_string()).into(), "No conversion found from source to target")]
    #[case(OperationError::Failed("test error".to_string()).into(), "Operation failed: test error")]
    #[case(DataError::InvalidAtomCharge(format!("{}", 2)).into(), "Invalid atom charge: 2")]
    #[case(DataError::InvalidAtomUnpairedElectrons(format!("{}", 2)).into(), "Invalid atom unpaired electron specification: 2")]
    fn test_error_display(#[case] error: Error, #[case] expected: &str) {
        assert_eq!(format!("{}", error), expected);
    }
}
