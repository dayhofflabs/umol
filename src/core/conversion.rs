// Core conversion traits

use crate::core::{Error, Model, Result};
use std::collections::HashMap;

/// A trait for basic conversion between molecular models
pub trait ConvertTo<M: Model> {
    /// Convert this model to another model
    fn convert_to(&self) -> Result<M>;
}

/// Metadata about a model conversion
#[derive(Debug, Clone, Default)]
pub struct ConversionMetadata {
    /// Uncertainty information for the conversion
    pub uncertainty: Option<f64>,
    /// Any additional metadata as string key-value pairs
    pub attributes: HashMap<String, String>,
}

/// A trait for conversions that require parameters and provide metadata
pub trait ConvertToWithMetadata<M: Model> {
    /// Parameters required for the conversion
    type Params: Default;
    
    /// Convert with specific parameters, returning both result and metadata
    fn convert_to_with_metadata(
        &self,
        params: &Self::Params
    ) -> Result<(M, ConversionMetadata)>;
}

// Convenience implementation - if something implements ConvertToWithMetadata,
// it can also do basic conversion
impl<T, M> ConvertTo<M> for T 
where 
    T: ConvertToWithMetadata<M>,
    M: Model,
{
    fn convert_to(&self) -> Result<M> {
        // Use default parameters and discard metadata
        let (model, _) = self.convert_to_with_metadata(&Default::default())?;
        Ok(model)
    }
} 