//! Capability system for describing model features.
//!
//! This module defines the `Capability` type which is used to describe
//! what features a model can provide. A capability consists of:
//! - An optional namespace
//! - A name
//! - A version number

use std::fmt;

use serde::{Deserialize, Serialize};

/// A capability that a model can provide
#[derive(Debug, Clone, Hash, Eq, PartialEq, Serialize, Deserialize)]
pub struct Capability {
    /// The namespace of the capability
    pub namespace: Option<String>,
    /// The name of the capability
    pub name: String,
    /// The version of the capability
    pub version: u32,
}

impl Capability {
    /// Create a new capability with a namespace
    pub fn new(namespace: impl Into<String>, name: impl Into<String>, version: u32) -> Self {
        Self {
            namespace: Some(namespace.into()),
            name: name.into(),
            version,
        }
    }

    /// Create a new capability without a namespace
    pub fn local(name: impl Into<String>, version: u32) -> Self {
        Self {
            namespace: None,
            name: name.into(),
            version,
        }
    }
}

impl fmt::Display for Capability {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(ref ns) = self.namespace {
            write!(f, "{}:{}:{}", ns, self.name, self.version)
        } else {
            write!(f, "{}:{}", self.name, self.version)
        }
    }
}

#[cfg(test)]
mod tests {
    use serde_json;

    use super::*;

    #[test]
    fn test_capability_creation() {
        let cap = Capability::new("test", "feature", 1);
        assert_eq!(cap.namespace, Some("test".to_string()));
        assert_eq!(cap.name, "feature");
        assert_eq!(cap.version, 1);

        let local = Capability::local("feature", 1);
        assert_eq!(local.namespace, None);
        assert_eq!(local.name, "feature");
        assert_eq!(local.version, 1);
    }

    #[test]
    fn test_capability_display() {
        let cap = Capability::new("test", "feature", 1);
        assert_eq!(cap.to_string(), "test:feature:1");

        let local = Capability::local("feature", 1);
        assert_eq!(local.to_string(), "feature:1");
    }

    #[test]
    fn test_capability_serialization() {
        let cap = Capability::new("test", "feature", 1);
        let json = serde_json::to_string(&cap).unwrap();
        assert_eq!(json, r#"{"namespace":"test","name":"feature","version":1}"#);

        let deserialized: Capability = serde_json::from_str(&json).unwrap();
        assert_eq!(deserialized, cap);
    }
}
