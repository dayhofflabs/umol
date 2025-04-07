//! Entity and relation types.
//!
//! Entities represent the semantic objects in molecular modeling:
//! - Molecules
//! - Conformers
//! - Resonance structures
//! - Reactions
//! Relations are typed associations between entities:
//! - Transformations
//! - Generalizations
//! - Specializations

use std::fmt;

/// Represents entities in the chemical domain (structures, conformers, etc.)
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Entity {
    /// The namespace this entity belongs to
    pub namespace: Option<String>,
    /// A unique identifier for this entity
    pub id: String,
    /// A human-readable label
    pub label: String,
}

impl Entity {
    /// Create a new entity
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

    pub fn local(id: impl Into<String>, label: impl Into<String>) -> Self {
        Self {
            namespace: None,
            id: id.into(),
            label: label.into(),
        }
    }
}

impl fmt::Display for Entity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(ns) = &self.namespace {
            write!(f, "{}::{}", ns, self.label)
        } else {
            write!(f, "{}", self.label)
        }
    }
}

/// Represent relationships between entities (transformations, reactions, etc.)
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Relation {
    /// The namespace this relation belongs to
    pub namespace: Option<String>,
    /// A unique identifier for this relation
    pub id: String,
    /// A human-readable label
    pub label: String,
    /// The source entity
    pub source: Entity,
    /// The target entity
    pub target: Entity,
}

impl Relation {
    /// Create a new relation
    pub fn new(
        namespace: impl Into<String>,
        id: impl Into<String>,
        label: impl Into<String>,
        source: Entity,
        target: Entity,
    ) -> Self {
        Self {
            namespace: Some(namespace.into()),
            id: id.into(),
            label: label.into(),
            source,
            target,
        }
    }

    pub fn local(
        id: impl Into<String>,
        label: impl Into<String>,
        source: Entity,
        target: Entity,
    ) -> Self {
        Self {
            namespace: None,
            id: id.into(),
            label: label.into(),
            source,
            target,
        }
    }
}

impl fmt::Display for Relation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {} {}", self.source, self.label, self.target)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_entity_display() {
        let entity = Entity::new("chebi", "CHEBI:24431", "chemical entity");
        assert_eq!(format!("{}", entity), "chebi::chemical entity");

        let entity = Entity::local("C1", "carbon");
        assert_eq!(format!("{}", entity), "carbon");
    }

    #[test]
    fn test_relation_display() {
        let source = Entity::new("chebi", "CHEBI:24431", "chemical entity");

        let target = Entity::new("chebi", "CHEBI:23367", "molecular entity");

        let relation = Relation::new("chebi", "CHEBI:23367", "is_a", source, target);
        assert_eq!(
            format!("{}", relation),
            "chebi::chemical entity is_a chebi::molecular entity"
        );
    }
}
