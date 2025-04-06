//! Serialization and deserialization traits and types.
//! 
//! This module provides abstractions for serializing and deserializing
//! molecular models and their components:
//! - Format-specific serialization traits
//! - Version-aware serialization
//! - Error handling for serialization
//! - Format metadata and validation

use std::io::{Read, Write};
use serde::{Serialize, Deserialize};
use chrono::{DateTime, Utc};
use crate::core::{Model, Result};
use std::fmt;

/// Represents a specific format version
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct FormatVersion {
    pub major: u32,
    pub minor: u32,
    pub patch: u32,
}

impl fmt::Display for FormatVersion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}.{}.{}", self.major, self.minor, self.patch)
    }
}

impl FormatVersion {
    pub fn new(major: u32, minor: u32, patch: u32) -> Self {
        Self { major, minor, patch }
    }

    pub fn is_compatible_with(&self, other: &FormatVersion) -> bool {
        self.major == other.major && self.minor >= other.minor
    }
}

/// Metadata about a serialized model
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializationMetadata {
    pub format: String,
    pub version: FormatVersion,
    pub timestamp: DateTime<Utc>,
    pub model_type: String,
    pub capabilities: Vec<String>,
}

/// Trait for serializing models to a specific format
pub trait ModelSerializer<M: Model> {
    /// The format name (e.g., "molden", "xyz")
    fn format_name(&self) -> &'static str;
    
    /// The format version
    fn format_version(&self) -> FormatVersion;
    
    /// Serialize a model to a writer
    fn serialize<W: Write>(&self, model: &M, writer: W) -> Result<()>;
    
    /// Deserialize a model from a reader
    fn deserialize<R: Read>(&self, reader: R) -> Result<M>;
}

/// Trait for models that can be serialized
pub trait SerializableModel: Model + Sized {
    /// Get the model's serialization metadata
    fn serialization_metadata(&self) -> SerializationMetadata;
    
    /// Serialize the model to a specific format
    fn serialize_to<W: Write, S: ModelSerializer<Self>>(
        &self,
        serializer: &S,
        writer: W,
    ) -> Result<()> {
        serializer.serialize(self, writer)
    }
    
    /// Deserialize the model from a specific format
    fn deserialize_from<R: Read, S: ModelSerializer<Self>>(
        serializer: &S,
        reader: R,
    ) -> Result<Self> {
        serializer.deserialize(reader)
    }
}

/// Trait for models that can be created from their data
pub trait FromModelData: Model {
    /// Create a new model from its data
    fn from_data(data: Self::Data) -> Self;
}

/// Trait for format-specific serialization implementations
pub trait FormatSerializer<M: Model> {
    /// The format name
    fn format_name(&self) -> &'static str;
    
    /// The format version
    fn format_version(&self) -> FormatVersion;
    
    /// Serialize a model to a writer
    fn serialize<W: Write>(&self, model: &M, writer: W) -> Result<()>;
    
    /// Deserialize a model from a reader
    fn deserialize<R: Read>(&self, reader: R) -> Result<M>;
    
    /// Validate the format of the input
    fn validate<R: Read>(&self, reader: R) -> Result<()>;
} 