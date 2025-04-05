// Core conversion trait

use crate::error::Error;
use crate::core::model::Model;

/// A trait for converting between different molecular models
pub trait ConvertTo<M: Model> {
    /// Convert this model to another model
    fn convert_to(&self) -> Result<M, Error>;
}

/// A trait for converting between different molecular models with metadata
pub trait ConvertToWithMetadata<M: Model> {
    /// The type of metadata associated with the conversion
    type Metadata;

    /// Convert this model to another model with metadata
    fn convert_to_with_metadata(&self) -> Result<(M, Self::Metadata), Error>;
}

/// Metadata about a model conversion
pub struct ConversionMetadata {
    // Add fields as needed
}

/// Uncertainty information for a model conversion
pub struct ConversionUncertainty {
    // Add fields as needed
}

/// Explicit conversion between models
pub trait Conversion<Source: Model, Target: Model> {
    type Error;
    
    fn convert(&self, source: &Source) -> Result<Target, Self::Error>;
    fn uncertainty(&self) -> ConversionUncertainty;
    fn metadata(&self) -> ConversionMetadata;
}
