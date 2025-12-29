//! Conversion traits and utilities.
//!
//! This module provides traits and utilities for converting between models:
//! - ConvertTo: Basic conversion between models
//! - ConvertToWithMetadata: Parameterized conversions with metadata

use std::collections::HashMap;

use crate::{Model, Result};

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

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use serde::{Deserialize, Serialize};

    use super::*;
    use crate::Capability;

    // Test models for conversion
    #[derive(Debug, Clone, Serialize, Deserialize)]
    struct SourceModel {
        data: SourceData,
    }

    #[derive(Debug, Clone, Serialize, Deserialize)]
    struct SourceData {
        value: i32,
    }

    impl Model for SourceModel {
        type Data = SourceData;

        fn data(&self) -> &Self::Data {
            &self.data
        }

        fn capabilities(&self) -> HashSet<Capability> {
            let mut caps = HashSet::new();
            caps.insert(Capability::local("source", 1));
            caps
        }
    }

    #[derive(Debug, Clone, Serialize, Deserialize)]
    struct TargetModel {
        data: TargetData,
    }

    #[derive(Debug, Clone, Serialize, Deserialize)]
    struct TargetData {
        value: i32,
        processed: bool,
    }

    impl Model for TargetModel {
        type Data = TargetData;

        fn data(&self) -> &Self::Data {
            &self.data
        }

        fn capabilities(&self) -> HashSet<Capability> {
            let mut caps = HashSet::new();
            caps.insert(Capability::local("target", 1));
            caps
        }
    }

    // Parameters for conversion
    #[derive(Debug, Default)]
    struct ConversionParams {
        process: bool,
    }

    // Implementation of ConvertToWithMetadata
    impl ConvertToWithMetadata<TargetModel> for SourceModel {
        type Params = ConversionParams;

        fn convert_to_with_metadata(
            &self,
            params: &Self::Params,
        ) -> Result<(TargetModel, ConversionMetadata)> {
            let mut metadata = ConversionMetadata::default();
            metadata.attributes.insert(
                "conversion_type".to_string(),
                "source_to_target".to_string(),
            );

            let target = TargetModel {
                data: TargetData {
                    value: self.data.value,
                    processed: params.process,
                },
            };

            Ok((target, metadata))
        }
    }

    #[test]
    fn test_basic_conversion() {
        let source = SourceModel {
            data: SourceData { value: 12 },
        };

        let target = source.convert_to().unwrap();
        assert_eq!(target.data.value, 12);
        assert!(!target.data.processed); // Should use default params
    }

    #[test]
    fn test_conversion_with_metadata() {
        let source = SourceModel {
            data: SourceData { value: 12 },
        };

        let params = ConversionParams { process: true };
        let (target, metadata) = source.convert_to_with_metadata(&params).unwrap();

        assert_eq!(target.data.value, 12);
        assert!(target.data.processed);
        assert_eq!(
            metadata.attributes.get("conversion_type"),
            Some(&"source_to_target".to_string())
        );
    }

    #[test]
    fn test_conversion_chain() {
        // Test that we can chain conversions through the convenience implementation
        let source = SourceModel {
            data: SourceData { value: 12 },
        };

        let target1 = source.convert_to().unwrap();
        let target2 = source.convert_to().unwrap();

        assert_eq!(target1.data.value, target2.data.value);
        assert_eq!(target1.data.processed, target2.data.processed);
    }

    #[test]
    fn test_metadata_attributes() {
        let source = SourceModel {
            data: SourceData { value: 12 },
        };

        let params = ConversionParams { process: true };
        let (_, metadata) = source.convert_to_with_metadata(&params).unwrap();

        assert_eq!(metadata.attributes.len(), 1);
        assert!(metadata.attributes.contains_key("conversion_type"));
    }
}
