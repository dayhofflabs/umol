//! Stuff type.
//!
//! Stuff combines entities with their model representations.

use crate::{AsEntity, AsModel, Entity, Model, Result};

/// Stuff combines entities with their model representations.
/// TODO: Rename to Item, Object, or Instance.
pub trait Stuff {
    /// Associated entity type
    type Entity: Entity;
    /// Associated model type
    type Model: Model;

    /// Get entity reference
    fn entity(&self) -> &Self::Entity;
    /// Get model reference
    fn model(&self) -> &Self::Model;

    /// Create new stuff from its components
    fn from_components(entity: Self::Entity, model: Self::Model) -> Result<Self>
    where
        Self: Sized;
}

// Stuff with associated entity E can act as E
impl<E: Entity, I: Stuff<Entity = E>> AsEntity<E> for I {
    fn as_entity(&self) -> &E {
        self.entity()
    }
}

// Stuff with associated model M can act as M
impl<M: Model, I: Stuff<Model = M>> AsModel<M> for I {
    fn as_model(&self) -> &M {
        self.model()
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use serde::{Deserialize, Serialize};
    use serde_json;

    use super::*;
    use crate::Capability;

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
                data: AtomCountData { count },
            }
        }
    }

    /// A simple stuff combining SimpleEntity with AtomCount
    #[derive(Debug, Clone, Serialize, Deserialize)]
    struct SimpleStuff {
        entity: SimpleEntity,
        model: AtomCount,
    }

    impl Stuff for SimpleStuff {
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
    fn test_stuff_serialization() {
        let entity = SimpleEntity::local("C1", "carbon");
        let model = AtomCount::new(1);
        let stuff = SimpleStuff::from_components(entity, model).unwrap();

        // Serialize to JSON
        let json = serde_json::to_string(&stuff).unwrap();
        assert_eq!(
            json,
            r#"{"entity":{"namespace":null,"id":"C1","label":"carbon"},"model":{"data":{"count":1}}}"#
        );

        // Deserialize from JSON
        let deserialized: SimpleStuff = serde_json::from_str(&json).unwrap();
        assert_eq!(deserialized.entity.id, "C1");
        assert_eq!(deserialized.entity.label, "carbon");
        assert_eq!(deserialized.model.data.count, 1);
    }
}
