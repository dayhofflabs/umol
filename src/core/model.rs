//! Core model traits and types.
//! 
//! This module defines the fundamental abstractions for molecular models:
//! - Model trait for representing molecular systems
//! - Basic model operations and validations

use std::collections::HashSet;
use serde::{Serialize, Deserialize};
use crate::core::{Result, Capability};

/// A trait for molecular models
pub trait Model: Serialize + for<'de> Deserialize<'de> {
    /// The type of data stored in this model
    type Data: Serialize + for<'de> Deserialize<'de>;
    
    /// Get a reference to the model's data
    fn data(&self) -> &Self::Data;
    
    /// Get the capabilities provided by this model
    fn capabilities(&self) -> HashSet<Capability>;
    
    /// Check if this model provides a specific capability
    fn has_capability(&self, capability: &Capability) -> bool {
        self.capabilities().contains(capability)
    }

    /// Validate the model
    fn validate(&self) -> Result<()> {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json;

    /// A simple model that just counts atoms
    #[derive(Debug, Clone, Serialize, Deserialize)]
    pub struct AtomCount {
        data: AtomCountData,
    }

    #[derive(Debug, Clone, Serialize, Deserialize)]
    pub struct AtomCountData {
        count: usize,
    }

    impl Model for AtomCount {
        type Data = AtomCountData;
        
        fn data(&self) -> &Self::Data {
            &self.data
        }
        
        fn capabilities(&self) -> HashSet<Capability> {
            let mut caps = HashSet::new();
            caps.insert(Capability::local("atom_count", 1));
            caps
        }
    }

    impl AtomCount {
        pub fn new(count: usize) -> Self {
            Self {
                data: AtomCountData { count }
            }
        }
    }

    #[test]
    fn test_atom_count_serialization() {
        let model = AtomCount::new(42);
        
        // Serialize to JSON
        let json = serde_json::to_string(&model).unwrap();
        assert_eq!(json, r#"{"data":{"count":42}}"#);
        
        // Deserialize from JSON
        let deserialized: AtomCount = serde_json::from_str(&json).unwrap();
        assert_eq!(deserialized.data.count, 42);
    }

    #[test]
    fn test_atom_count_capabilities() {
        let model = AtomCount::new(42);
        let caps = model.capabilities();
        assert!(caps.contains(&Capability::local("atom_count", 1)));
    }
}
