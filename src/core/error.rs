use thiserror::Error;

/// Core error types for the molecular modeling framework
#[derive(Error, Debug)]
pub enum Error {
    #[error("Entity error: {0}")]
    Entity(#[from] EntityError),

    #[error("Model error: {0}")]
    Model(#[from] ModelError),

    #[error("Conversion error: {0}")]
    Conversion(#[from] ConversionError),

    #[error("Operation error: {0}")]
    Operation(#[from] OperationError),

    #[error("Property error: {0}")]
    Property(#[from] PropertyError),

    #[error("Validation error: {0}")]
    ValidationError(String),

    #[error("Multiple errors occurred: {0:?}")]
    Multiple(Vec<Error>),

    #[error(transparent)]
    Other(#[from] Box<dyn std::error::Error + Send + Sync>),
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