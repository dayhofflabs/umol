//! Conversion traits and utilities.
//!
//! This module provides traits and utilities for converting between models:
//! - Basic conversion between models
//! - Parameterized conversions with metadata

use crate::core::{Model, Result};
use std::collections::HashMap;

/// A trait for converting between models
pub trait ConvertTo<M2: Model> {
    /// Convert this model to another model
    fn convert_to(&self) -> Result<M2>;
}

/// Metadata describing model conversion
#[derive(Debug, Clone, Default)]
pub struct ConversionMetadata {
    /// Metadata as string key-value pairs
    pub attributes: HashMap<String, String>,
}

/// A trait for conversions that require parameters and provide metadata
pub trait ConvertToWithMetadata<M: Model> {
    /// Parameters required for the conversion
    type Params: Default;

    /// Convert with specific parameters, returning both result and metadata
    fn convert_to_with_metadata(&self, params: &Self::Params) -> Result<(M, ConversionMetadata)>;
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
