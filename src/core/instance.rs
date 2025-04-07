//! Instance types and operations.
//!
//! Instances combine entities with their model representations.

use crate::core::error::Result;
use crate::core::{Entity, Model};

/// Core functionality for instances in the chemical domain
pub trait Instance {
    /// The entity type associated with this instance
    type Entity: Entity;
    /// The model type used for representation
    type Model: Model;

    /// Get a reference to the entity
    fn entity(&self) -> &Self::Entity;
    /// Get a reference to the model
    fn model(&self) -> &Self::Model;

    /// Create a new instance from its components
    fn from_components(entity: Self::Entity, model: Self::Model) -> Result<Self>
    where
        Self: Sized;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::Capability;
    use serde::{Deserialize, Serialize};
    use serde_json;
    use std::collections::HashSet;

    // Test entity
    #[derive(Debug, Clone, Serialize, Deserialize)]
    struct SimpleEntity {
        namespace: Option<String>,
        id: String,
        label: String,
    }

    impl SimpleEntity {
        fn local(id: impl Into<String>, label: impl Into<String>) -> Self {
            Self {
                namespace: None,
                id: id.into(),
                label: label.into(),
            }
        }
    }

    impl Entity for SimpleEntity {
        fn namespace(&self) -> Option<&str> {
            self.namespace.as_deref()
        }

        fn id(&self) -> &str {
            &self.id
        }

        fn label(&self) -> &str {
            &self.label
        }
    }

    // Test model
    #[derive(Debug, Clone, Serialize, Deserialize)]
    struct AtomCount {
        data: AtomCountData,
    }

    #[derive(Debug, Clone, Serialize, Deserialize)]
    struct AtomCountData {
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
        fn new(count: usize) -> Self {
            Self {
                data: AtomCountData { count }
            }
        }
    }

    /// A simple instance combining SimpleEntity with AtomCount
    #[derive(Debug, Clone, Serialize, Deserialize)]
    struct SimpleInstance {
        entity: SimpleEntity,
        model: AtomCount,
    }

    impl Instance for SimpleInstance {
        type Entity = SimpleEntity;
        type Model = AtomCount;

        fn entity(&self) -> &Self::Entity {
            &self.entity
        }

        fn model(&self) -> &Self::Model {
            &self.model
        }

        fn from_components(entity: Self::Entity, model: Self::Model) -> Result<Self> {
            Ok(Self { entity, model })
        }
    }

    #[test]
    fn test_instance_serialization() {
        let entity = SimpleEntity::local("C1", "carbon");
        let model = AtomCount::new(1);
        let instance = SimpleInstance::from_components(entity, model).unwrap();

        // Serialize to JSON
        let json = serde_json::to_string(&instance).unwrap();
        assert_eq!(
            json,
            r#"{"entity":{"namespace":null,"id":"C1","label":"carbon"},"model":{"data":{"count":1}}}"#
        );

        // Deserialize from JSON
        let deserialized: SimpleInstance = serde_json::from_str(&json).unwrap();
        assert_eq!(deserialized.entity.id, "C1");
        assert_eq!(deserialized.entity.label, "carbon");
        assert_eq!(deserialized.model.data.count, 1);
    }
}
