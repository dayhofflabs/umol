//! Entity and relation types.
//!
//! Entities represent the semantic objects in molecular modeling:
//! - Molecules
//! - Conformers
//! - Resonance structures
//! - Reactions
//! 
//! Relations are typed associations between entities:
//! - Transformations
//! - Generalizations
//! - Specializations

use serde::{Deserialize, Serialize};

/// Entities in the chemical domain
pub trait Entity: Serialize + for<'de> Deserialize<'de> {
    /// Get the namespace of this entity
    fn namespace(&self) -> Option<&str>;
    /// Get the unique identifier of this entity
    fn id(&self) -> &str;
    /// Get the human-readable label of this entity
    fn label(&self) -> &str;
}

pub trait AsEntity<E: Entity> {
    /// Get a reference to the underlying entity
    fn as_entity(&self) -> &E;
}

/// Relationships between entities
pub trait Relation: Serialize + for<'de> Deserialize<'de> {
    /// Entity type associated with this relation
    type Side: Entity;
    /// Get the namespace of this relation
    fn namespace(&self) -> Option<&str>;
    /// Get the unique identifier of this relation
    fn id(&self) -> &str;
    /// Get the human-readable label of this relation
    fn label(&self) -> &str;
    /// Get the source entity
    fn source(&self) -> &Self::Side;
    /// Get the target entity
    fn target(&self) -> &Self::Side;
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json;

    /// A simple concrete implementation of the Entity trait
    #[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
    pub struct SimpleEntity {
        /// The namespace this entity belongs to
        namespace: Option<String>,
        /// A unique identifier for this entity
        id: String,
        /// A human-readable label
        label: String,
    }

    impl SimpleEntity {
        /// Create a new entity with a namespace
        pub fn new(
            namespace: impl Into<String>,
            id: impl Into<String>,
            label: impl Into<String>,
        ) -> Self {
            Self {
                namespace: Some(namespace.into()),
                id: id.into(),
                label: label.into(),
            }
        }

        /// Create a new local entity (without namespace)
        pub fn local(id: impl Into<String>, label: impl Into<String>) -> Self {
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

    impl std::fmt::Display for SimpleEntity {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            if let Some(ns) = &self.namespace {
                write!(f, "{}::{}", ns, self.label)
            } else {
                write!(f, "{}", self.label)
            }
        }
    }

    /// A simple concrete implementation of the Relation trait
    #[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
    pub struct SimpleRelation {
        /// The namespace this relation belongs to
        namespace: Option<String>,
        /// A unique identifier for this relation
        id: String,
        /// A human-readable label
        label: String,
        /// The source entity
        source: SimpleEntity,
        /// The target entity
        target: SimpleEntity,
    }

    impl SimpleRelation {
        /// Create a new relation with a namespace
        pub fn new(
            namespace: impl Into<String>,
            id: impl Into<String>,
            label: impl Into<String>,
            source: SimpleEntity,
            target: SimpleEntity,
        ) -> Self {
            Self {
                namespace: Some(namespace.into()),
                id: id.into(),
                label: label.into(),
                source,
                target,
            }
        }
    }

    impl Relation for SimpleRelation {
        type Side = SimpleEntity;

        fn namespace(&self) -> Option<&str> {
            self.namespace.as_deref()
        }

        fn id(&self) -> &str {
            &self.id
        }

        fn label(&self) -> &str {
            &self.label
        }

        fn source(&self) -> &Self::Side {
            &self.source
        }

        fn target(&self) -> &Self::Side {
            &self.target
        }
    }

    impl std::fmt::Display for SimpleRelation {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(f, "{} {} {}", self.source, self.label, self.target)
        }
    }

    #[test]
    fn test_entity_display() {
        let entity = SimpleEntity::new("chebi", "CHEBI:24431", "chemical entity");
        assert_eq!(format!("{}", entity), "chebi::chemical entity");

        let entity = SimpleEntity::local("C1", "carbon");
        assert_eq!(format!("{}", entity), "carbon");
    }

    #[test]
    fn test_entity_serialization() {
        let entity = SimpleEntity::new("chebi", "CHEBI:24431", "chemical entity");

        // Serialize to JSON
        let json = serde_json::to_string(&entity).unwrap();
        assert_eq!(
            json,
            r#"{"namespace":"chebi","id":"CHEBI:24431","label":"chemical entity"}"#
        );

        // Deserialize from JSON
        let deserialized: SimpleEntity = serde_json::from_str(&json).unwrap();
        assert_eq!(deserialized.namespace, Some("chebi".to_string()));
        assert_eq!(deserialized.id, "CHEBI:24431");
        assert_eq!(deserialized.label, "chemical entity");

        // Test local entity
        let local = SimpleEntity::local("C1", "carbon");
        let json = serde_json::to_string(&local).unwrap();
        assert_eq!(json, r#"{"namespace":null,"id":"C1","label":"carbon"}"#);
    }

    #[test]
    fn test_relation_display() {
        let source = SimpleEntity::new("chebi", "CHEBI:24431", "chemical entity");
        let target = SimpleEntity::new("chebi", "CHEBI:23367", "molecular entity");
        let relation = SimpleRelation::new("chebi", "CHEBI:23367", "is_a", source, target);
        assert_eq!(
            format!("{}", relation),
            "chebi::chemical entity is_a chebi::molecular entity"
        );
    }

    #[test]
    fn test_relation_serialization() {
        let source = SimpleEntity::new("chebi", "CHEBI:24431", "chemical entity");
        let target = SimpleEntity::new("chebi", "CHEBI:23367", "molecular entity");
        let relation = SimpleRelation::new("chebi", "CHEBI:23367", "is_a", source, target);

        // Serialize to JSON
        let json = serde_json::to_string(&relation).unwrap();
        assert_eq!(
            json,
            r#"{"namespace":"chebi","id":"CHEBI:23367","label":"is_a","source":{"namespace":"chebi","id":"CHEBI:24431","label":"chemical entity"},"target":{"namespace":"chebi","id":"CHEBI:23367","label":"molecular entity"}}"#
        );

        // Deserialize from JSON
        let deserialized: SimpleRelation = serde_json::from_str(&json).unwrap();
        assert_eq!(deserialized.namespace, Some("chebi".to_string()));
        assert_eq!(deserialized.id, "CHEBI:23367");
        assert_eq!(deserialized.label, "is_a");
        assert_eq!(deserialized.source.id, "CHEBI:24431");
        assert_eq!(deserialized.target.id, "CHEBI:23367");
    }
}
