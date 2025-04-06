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
    /// A unique identifier for this entity
    pub id: String,
    /// A human-readable label
    pub label: String,
    /// The namespace this entity belongs to
    pub namespace: Option<String>,
}

impl Entity {
    /// Create a new entity
    pub fn new(
        id: impl Into<String>,
        label: impl Into<String>,
        namespace: Option<String>,
    ) -> Self {
        Self {
            id: id.into(),
            label: label.into(),
            namespace,
        }
    }

    /// Set the namespace
    pub fn with_namespace(mut self, namespace: impl Into<String>) -> Self {
        self.namespace = Some(namespace.into());
        self
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
    /// A unique identifier for this relation
    pub id: String,
    /// A human-readable label
    pub label: String,
    /// The namespace this relation belongs to
    pub namespace: Option<String>,
    /// The source entity
    pub source: Entity,
    /// The target entity
    pub target: Entity,
}

impl Relation {
    /// Create a new relation
    pub fn new(
        id: impl Into<String>,
        label: impl Into<String>,
        namespace: Option<String>,
        source: Entity,
        target: Entity,
    ) -> Self {
        Self {
            id: id.into(),
            label: label.into(),
            namespace,
            source,
            target,
        }
    }

    /// Set the namespace
    pub fn with_namespace(mut self, namespace: impl Into<String>) -> Self {
        self.namespace = Some(namespace.into());
        self
    }
}

impl fmt::Display for Relation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} {} {}",
            self.source,
            self.label,
            self.target
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_entity_display() {
        let entity = Entity::new(
            "CHEBI:24431",
            "chemical entity",
            Some("chebi".into()),
        );
        assert_eq!(format!("{}", entity), "chebi::chemical entity");

        let entity = Entity::new(
            "C1",
            "carbon",
            None,
        );
        assert_eq!(format!("{}", entity), "carbon");
    }

    #[test]
    fn test_relation_display() {
        let source = Entity::new(
            "CHEBI:24431",
            "chemical entity",
            Some("chebi".into()),
        );

        let target = Entity::new(
            "CHEBI:23367",
            "molecular entity",
            Some("chebi".into()),
        );
        
        let relation = Relation::new(
            "CHEBI:23367",
            "is_a",
            Some("chebi".into()),
            source,
            target,
        );
        assert_eq!(
            format!("{}", relation),
            "chebi::chemical entity is_a chebi::molecular entity"
        );
    }
}
