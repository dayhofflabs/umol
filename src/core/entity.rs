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

use std::collections::HashMap;
use std::fmt;
use serde::{Serialize, Deserialize};

/// A value that can be stored as an attribute
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum AttributeValue {
    String(String),
    Integer(i64),
    Float(f64),
    Boolean(bool),
    List(Vec<AttributeValue>),
    Map(HashMap<String, AttributeValue>),
}

/// Trait for entity types
pub trait EntityType: Send + Sync {
    /// Get the name of this entity type
    fn name(&self) -> &str;
    
    /// Get the namespace of this entity type
    fn namespace(&self) -> Option<&str>;
    
    /// Get the allowed attributes for this entity type
    fn allowed_attributes(&self) -> &[&str];
}

/// Trait for relation types
pub trait RelationType: Send + Sync {
    /// Get the name of this relation type
    fn name(&self) -> &str;
    
    /// Get the namespace of this relation type
    fn namespace(&self) -> Option<&str>;
    
    /// Get the source entity type
    fn source_type(&self) -> &str;
    
    /// Get the target entity type
    fn target_type(&self) -> &str;
    
    /// Get the allowed attributes for this relation type
    fn allowed_attributes(&self) -> &[&str];
}

/// Represents entities in the chemical domain (structures, conformers, etc.)
#[derive(Debug, Clone)]
pub struct Entity {
    /// A unique identifier for this entity
    pub label: String,
    /// The namespace this entity belongs to
    pub namespace: Option<String>,
    /// The type of this entity
    entity_type: Box<dyn EntityType>,
    /// The attributes of this entity
    attributes: HashMap<String, AttributeValue>,
}

impl Entity {
    /// Create a new entity
    pub fn new(
        label: impl Into<String>,
        namespace: Option<String>,
        entity_type: impl EntityType + 'static,
    ) -> Self {
        Self {
            label: label.into(),
            namespace,
            entity_type: Box::new(entity_type),
            attributes: HashMap::new(),
        }
    }
    
    /// Get the entity type
    pub fn entity_type(&self) -> &dyn EntityType {
        &*self.entity_type
    }
    
    /// Get an attribute value
    pub fn get_attribute(&self, name: &str) -> Option<&AttributeValue> {
        self.attributes.get(name)
    }
    
    /// Set an attribute value
    pub fn set_attribute(&mut self, name: impl Into<String>, value: AttributeValue) {
        self.attributes.insert(name.into(), value);
    }
    
    /// Validate the entity against its type requirements
    pub fn validate(&self) -> Result<(), String> {
        // Check that all attributes are allowed
        for attr in self.attributes.keys() {
            if !self.entity_type.allowed_attributes().contains(&attr.as_str()) {
                return Err(format!("Unknown attribute: {}", attr));
            }
        }
        
        Ok(())
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
pub struct Relation {
    source: Entity,
    target: Entity,
    relation_type: Box<dyn RelationType>,
    attributes: HashMap<String, AttributeValue>,
}

impl Relation {
    /// Create a new relation
    pub fn new(
        source: Entity,
        target: Entity,
        relation_type: impl RelationType + 'static,
    ) -> Self {
        Self {
            source,
            target,
            relation_type: Box::new(relation_type),
            attributes: HashMap::new(),
        }
    }
    
    /// Get the relation type
    pub fn relation_type(&self) -> &dyn RelationType {
        &*self.relation_type
    }
    
    /// Get an attribute value
    pub fn get_attribute(&self, name: &str) -> Option<&AttributeValue> {
        self.attributes.get(name)
    }
    
    /// Set an attribute value
    pub fn set_attribute(&mut self, name: impl Into<String>, value: AttributeValue) {
        self.attributes.insert(name.into(), value);
    }
    
    /// Validate the relation against its type requirements
    pub fn validate(&self) -> Result<(), String> {
        // Check source and target types
        if self.source.entity_type().name() != self.relation_type.source_type() {
            return Err(format!(
                "Invalid source type: expected {}, got {}",
                self.relation_type.source_type(),
                self.source.entity_type().name()
            ));
        }
        
        if self.target.entity_type().name() != self.relation_type.target_type() {
            return Err(format!(
                "Invalid target type: expected {}, got {}",
                self.relation_type.target_type(),
                self.target.entity_type().name()
            ));
        }
        
        // Check that all attributes are allowed
        for attr in self.attributes.keys() {
            if !self.relation_type.allowed_attributes().contains(&attr.as_str()) {
                return Err(format!("Unknown attribute: {}", attr));
            }
        }
        
        Ok(())
    }
}

impl fmt::Display for Relation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} {} {}",
            self.source,
            self.relation_type.name(),
            self.target
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Simple test implementations of EntityType and RelationType
    #[derive(Debug)]
    struct TestEntityType {
        name: String,
        namespace: Option<String>,
        allowed_attrs: Vec<&'static str>,
    }

    impl EntityType for TestEntityType {
        fn name(&self) -> &str {
            &self.name
        }
        
        fn namespace(&self) -> Option<&str> {
            self.namespace.as_deref()
        }
        
        fn allowed_attributes(&self) -> &[&str] {
            &self.allowed_attrs
        }
    }

    #[derive(Debug)]
    struct TestRelationType {
        name: String,
        namespace: Option<String>,
        source_type: String,
        target_type: String,
        allowed_attrs: Vec<&'static str>,
    }

    impl RelationType for TestRelationType {
        fn name(&self) -> &str {
            &self.name
        }
        
        fn namespace(&self) -> Option<&str> {
            self.namespace.as_deref()
        }
        
        fn source_type(&self) -> &str {
            &self.source_type
        }
        
        fn target_type(&self) -> &str {
            &self.target_type
        }
        
        fn allowed_attributes(&self) -> &[&str] {
            &self.allowed_attrs
        }
    }

    #[test]
    fn test_entity_display() {
        let entity_type = TestEntityType {
            name: "molecule".into(),
            namespace: Some("chemical".into()),
            allowed_attrs: vec!["smiles", "inchi"],
        };

        let entity = Entity::new("benzene", Some("molecule".into()), entity_type);
        assert_eq!(format!("{}", entity), "molecule::benzene");

        let entity_type = TestEntityType {
            name: "atom".into(),
            namespace: None,
            allowed_attrs: vec!["element", "charge"],
        };

        let entity = Entity::new("C1", None, entity_type);
        assert_eq!(format!("{}", entity), "C1");
    }

    #[test]
    fn test_entity_attributes() {
        let entity_type = TestEntityType {
            name: "molecule".into(),
            namespace: None,
            allowed_attrs: vec!["smiles", "inchi"],
        };

        let mut entity = Entity::new("benzene", None, entity_type);
        
        // Test setting and getting attributes
        entity.set_attribute("smiles", AttributeValue::String("c1ccccc1".into()));
        assert!(matches!(
            entity.get_attribute("smiles"),
            Some(AttributeValue::String(s)) if s == "c1ccccc1"
        ));

        // Test validation
        assert!(entity.validate().is_ok());
        
        // Test invalid attribute
        entity.set_attribute("invalid", AttributeValue::String("value".into()));
        assert!(entity.validate().is_err());
    }

    #[test]
    fn test_relation_display() {
        let source_type = TestEntityType {
            name: "molecule".into(),
            namespace: Some("chemical".into()),
            allowed_attrs: vec![],
        };

        let target_type = TestEntityType {
            name: "molecule".into(),
            namespace: Some("chemical".into()),
            allowed_attrs: vec![],
        };

        let relation_type = TestRelationType {
            name: "transforms_to".into(),
            namespace: Some("reaction".into()),
            source_type: "molecule".into(),
            target_type: "molecule".into(),
            allowed_attrs: vec!["energy", "barrier"],
        };

        let source = Entity::new("benzene", Some("molecule".into()), source_type);
        let target = Entity::new("toluene", Some("molecule".into()), target_type);
        
        let relation = Relation::new(source, target, relation_type);
        assert_eq!(
            format!("{}", relation),
            "molecule::benzene transforms_to molecule::toluene"
        );
    }

    #[test]
    fn test_relation_validation() {
        let source_type = TestEntityType {
            name: "molecule".into(),
            namespace: None,
            allowed_attrs: vec![],
        };

        let target_type = TestEntityType {
            name: "molecule".into(),
            namespace: None,
            allowed_attrs: vec![],
        };

        let relation_type = TestRelationType {
            name: "transforms_to".into(),
            namespace: None,
            source_type: "molecule".into(),
            target_type: "molecule".into(),
            allowed_attrs: vec!["energy"],
        };

        let source = Entity::new("benzene", None, source_type);
        let target = Entity::new("toluene", None, target_type);
        
        let mut relation = Relation::new(source, target, relation_type);
        
        // Test valid relation
        assert!(relation.validate().is_ok());
        
        // Test valid attribute
        relation.set_attribute("energy", AttributeValue::Float(42.0));
        assert!(relation.validate().is_ok());
        
        // Test invalid attribute
        relation.set_attribute("invalid", AttributeValue::String("value".into()));
        assert!(relation.validate().is_err());
    }
}
