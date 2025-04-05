use std::collections::HashSet;
use crate::core::{
    Capability, Entity, Model, Property,
    Instance, Result, Error,
    ConvertTo, ConvertToWithMetadata, ConversionMetadata,
};

/// Trait for testing model implementations
pub trait ModelTest {
    type E: Entity;
    type M: Model;
    
    /// Create a new instance for testing
    fn create_test_instance() -> Result<Instance<Self::E, Self::M>>;
    
    /// Test that the model correctly reports its capabilities
    fn test_capabilities() -> Result<()>;
    
    /// Test basic model operations
    fn test_model_operations() -> Result<()>;
    
    /// Test property calculations
    fn test_property_calculations() -> Result<()>;
}

/// Trait for testing property implementations
pub trait PropertyTest {
    type P: Property;
    
    /// Test that the property correctly reports its requirements
    fn test_requirements() -> Result<()>;
    
    /// Test property calculation on a simple case
    fn test_simple_calculation() -> Result<()>;
    
    /// Test property calculation on edge cases
    fn test_edge_cases() -> Result<()>;
}

/// Helper function to verify capability requirements
pub fn verify_capabilities<M: Model>(model: &M, required: &[Capability]) -> Result<()> {
    let required: HashSet<_> = required.iter().cloned().collect();
    let available = model.capabilities();
    
    if !required.is_subset(&available) {
        let missing: Vec<_> = required.difference(&available).collect();
        return Err(crate::core::Error::Model(
            crate::core::ModelError::MissingCapability(
                missing.first().cloned().unwrap()
            )
        ));
    }
    
    Ok(())
}

/// Helper function to verify property calculation
pub fn verify_property_calculation<P: Property, E: Entity, M: Model>(
    instance: &Instance<E, M>
) -> Result<()> {
    // Verify the model has required capabilities
    verify_capabilities(instance.model(), &P::required_capabilities().into_iter().collect::<Vec<_>>())?;
    
    // Try to compute the property
    P::compute(instance)?;
    
    Ok(())
}

// Mock implementations for testing
#[derive(Debug, Clone, PartialEq)]
pub struct MockEntity {
    id: String,
}

impl MockEntity {
    pub fn new(id: &str) -> Self {
        Self { id: id.to_string() }
    }
}

impl Entity for MockEntity {
    fn generalizes(&self, other: &Self) -> bool {
        self.id.len() <= other.id.len()
    }

    fn specializes(&self, other: &Self) -> bool {
        self.id.len() >= other.id.len()
    }
}

#[derive(Debug, Clone)]
pub struct MockModel {
    capabilities: HashSet<Capability>,
    data: String,
}

impl MockModel {
    pub fn new(data: &str, capabilities: &[Capability]) -> Self {
        Self {
            data: data.to_string(),
            capabilities: capabilities.iter().cloned().collect(),
        }
    }
}

impl Model for MockModel {
    fn capabilities(&self) -> HashSet<Capability> {
        self.capabilities.clone()
    }
}

#[derive(Debug, Clone)]
pub struct MockModelAdvanced {
    capabilities: HashSet<Capability>,
    data: String,
}

impl MockModelAdvanced {
    pub fn new(data: &str, capabilities: &[Capability]) -> Self {
        Self {
            data: data.to_string(),
            capabilities: capabilities.iter().cloned().collect(),
        }
    }
}

impl Model for MockModelAdvanced {
    fn capabilities(&self) -> HashSet<Capability> {
        self.capabilities.clone()
    }
}

// Implement conversion between mock models
#[derive(Default)]
pub struct MockConversionParams {
    pub preserve_capabilities: bool,
}

impl ConvertToWithMetadata<MockModelAdvanced> for MockModel {
    type Params = MockConversionParams;

    fn convert_to_with_metadata(
        &self,
        params: &Self::Params
    ) -> Result<(MockModelAdvanced, ConversionMetadata)> {
        let capabilities = if params.preserve_capabilities {
            self.capabilities.clone()
        } else {
            HashSet::new()
        };
        
        let mut metadata = ConversionMetadata::default();
        metadata.attributes.insert("source".to_string(), "MockModel".to_string());
        
        Ok((MockModelAdvanced::new(&self.data, &capabilities.iter().collect::<Vec<_>>()), metadata))
    }
}

// Mock property for testing
pub struct MockProperty;

impl Property for MockProperty {
    type Value = f64;

    fn name() -> &'static str {
        "Mock Property"
    }

    fn description() -> &'static str {
        "A mock property for testing"
    }

    fn units() -> Option<&'static str> {
        Some("mock_units")
    }

    fn required_capabilities() -> HashSet<Capability> {
        let mut caps = HashSet::new();
        caps.insert(Capability::HasAtoms);
        caps
    }

    fn compute<E: Entity, M: Model>(instance: &Instance<E, M>) -> Result<Self::Value> {
        if instance.model().has_capability(&Capability::HasAtoms) {
            Ok(42.0)
        } else {
            Err(Error::Property(crate::core::PropertyError::MissingCapability(Capability::HasAtoms)))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mock_entity_relations() {
        let general = MockEntity::new("a");
        let specific = MockEntity::new("abc");

        assert!(general.generalizes(&specific));
        assert!(specific.specializes(&general));
        assert!(!general.specializes(&specific));
        assert!(!specific.generalizes(&general));
    }

    #[test]
    fn test_mock_model_capabilities() {
        let caps = vec![Capability::HasAtoms, Capability::HasBonds];
        let model = MockModel::new("test", &caps);

        assert!(model.has_capability(&Capability::HasAtoms));
        assert!(model.has_capability(&Capability::HasBonds));
        assert!(!model.has_capability(&Capability::HasCoordinates3D));
    }

    #[test]
    fn test_mock_model_conversion() {
        let caps = vec![Capability::HasAtoms];
        let model = MockModel::new("test", &caps);
        
        // Test basic conversion
        let advanced = model.convert_to().unwrap();
        assert_eq!(model.capabilities(), advanced.capabilities());

        // Test parameterized conversion
        let params = MockConversionParams { preserve_capabilities: false };
        let (advanced, metadata) = model.convert_to_with_metadata(&params).unwrap();
        assert!(advanced.capabilities().is_empty());
        assert_eq!(metadata.attributes.get("source"), Some(&"MockModel".to_string()));
    }

    #[test]
    fn test_mock_property() {
        let entity = MockEntity::new("test");
        let model = MockModel::new("test", &[Capability::HasAtoms]);
        let instance = Instance::new(entity, model);

        let value = MockProperty::compute(&instance).unwrap();
        assert_eq!(value, 42.0);

        let model_without_atoms = MockModel::new("test", &[]);
        let instance_without_atoms = Instance::new(MockEntity::new("test"), model_without_atoms);
        assert!(MockProperty::compute(&instance_without_atoms).is_err());
    }

    #[test]
    fn test_instance_conversion() {
        let entity = MockEntity::new("test");
        let model = MockModel::new("test", &[Capability::HasAtoms]);
        let instance = Instance::new(entity.clone(), model);

        let converted = instance.convert_to().unwrap();
        assert_eq!(converted.entity(), &entity);
        assert!(converted.model().has_capability(&Capability::HasAtoms));
    }
} 