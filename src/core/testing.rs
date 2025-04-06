//! Testing utilities and traits.
//! 
//! This module provides testing infrastructure for umol:
//! - Model testing traits
//! - Property testing traits
//! - Test utilities and helpers

use std::collections::HashSet;
use crate::core::{
    Model, Property, Capability, Result, Error,
    Entity, Instance,
    error::{ModelError, PropertyError},
    ConvertToWithMetadata, ConversionMetadata,
};

/// Trait for testing model implementations
pub trait ModelTest {
    /// The entity type being tested
    type E: Entity;
    /// The model type being tested
    type M: Model;
    
    /// Create a test instance for testing
    fn create_test_instance() -> Result<Instance<Self::E, Self::M>>;
    
    /// Test that the model has the required capabilities
    fn test_capabilities() -> Result<()>;
    
    /// Test that model operations work correctly
    fn test_model_operations() -> Result<()>;
    
    /// Test that property calculations work correctly
    fn test_property_calculations() -> Result<()>;
}

/// Trait for testing property implementations
pub trait PropertyTest {
    /// The property type being tested
    type P: Property;
    
    /// Test that the property requirements are met
    fn test_requirements() -> Result<()>;
    
    /// Test that simple calculations work correctly
    fn test_simple_calculation() -> Result<()>;
    
    /// Test that edge cases are handled correctly
    fn test_edge_cases() -> Result<()>;
}

/// Helper functions for testing
pub mod helpers {
    use super::*;
    
    /// Verify that a model has all required capabilities
    pub fn verify_capabilities<M: Model>(
        model: &M,
        required: &[Capability],
    ) -> Result<()> {
        for cap in required {
            if !model.has_capability(cap) {
                return Err(Error::Model(ModelError::MissingCapability(cap.clone())));
            }
        }
        Ok(())
    }
    
    /// Verify that a property can be calculated on a model
    pub fn verify_property_calculation<P: Property, E: Entity, M: Model>(
        instance: &Instance<E, M>,
    ) -> Result<()> {
        P::compute(instance)?;
        Ok(())
    }
}

/// Default implementation of ModelTest
pub struct DefaultModelTest<E: Entity, M: Model> {
    _phantom: std::marker::PhantomData<(E, M)>,
}

impl<E: Entity, M: Model> DefaultModelTest<E, M> {
    /// Create a new default model test
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<E: Entity, M: Model> ModelTest for DefaultModelTest<E, M> {
    type E = E;
    type M = M;
    
    fn create_test_instance() -> Result<Instance<Self::E, Self::M>> {
        unimplemented!("Default implementation does not provide test instances")
    }
    
    fn test_capabilities() -> Result<()> {
        Ok(())
    }
    
    fn test_model_operations() -> Result<()> {
        Ok(())
    }
    
    fn test_property_calculations() -> Result<()> {
        Ok(())
    }
}

/// Default implementation of PropertyTest
pub struct DefaultPropertyTest<P: Property> {
    _phantom: std::marker::PhantomData<P>,
}

impl<P: Property> DefaultPropertyTest<P> {
    /// Create a new default property test
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<P: Property> PropertyTest for DefaultPropertyTest<P> {
    type P = P;
    
    fn test_requirements() -> Result<()> {
        Ok(())
    }
    
    fn test_simple_calculation() -> Result<()> {
        Ok(())
    }
    
    fn test_edge_cases() -> Result<()> {
        Ok(())
    }
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

/// A mock model for testing
pub struct MockModel {
    capabilities: HashSet<Capability>,
}

impl MockModel {
    pub fn new(capabilities: &[Capability]) -> Self {
        Self {
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
}

impl MockModelAdvanced {
    pub fn new(capabilities: &[Capability]) -> Self {
        Self {
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
        
        let caps: Vec<Capability> = capabilities.into_iter().collect();
        Ok((MockModelAdvanced::new(&caps), metadata))
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
        caps.insert(Capability {
            name: "has_atoms".to_string(),
            version: "1.0".to_string(),
        });
        caps
    }

    fn compute<E: Entity, M: Model>(instance: &Instance<E, M>) -> Result<Self::Value> {
        let has_atoms = Capability {
            name: "has_atoms".to_string(),
            version: "1.0".to_string(),
        };
        if instance.model.has_capability(&has_atoms) {
            Ok(42.0)
        } else {
            Err(Error::Property(PropertyError::MissingCapability(has_atoms)))
        }
    }
}

/// Helper function to verify property calculation
pub fn verify_property_calculation<P: Property>(
    instance: &Instance<impl Entity, impl Model>,
    expected: P::Value,
) -> Result<()> 
where P::Value: PartialEq {
    let result = P::compute(instance)?;
    if result == expected {
        Ok(())
    } else {
        Err(Error::Property(PropertyError::CalculationFailed(
            "Property calculation result mismatch".to_string(),
        )))
    }
}

/// Helper function to verify model capabilities
pub fn verify_capabilities(model: &impl Model, required: &[Capability]) -> Result<()> {
    for cap in required {
        if !model.has_capability(cap) {
            return Err(Error::Model(ModelError::MissingCapability(cap.clone())));
        }
    }
    Ok(())
}

/// Test that verifies a model has the required capabilities
pub fn test_model_capabilities(model: &impl Model) -> Result<()> {
    let mut caps = HashSet::new();
    caps.insert(Capability {
        name: "has_atoms".to_string(),
        version: "1.0".to_string(),
    });
    
    for cap in &caps {
        if !model.has_capability(cap) {
            return Err(Error::Model(ModelError::MissingCapability(cap.clone())));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::ConvertTo;

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
    fn test_model_capabilities() {
        let has_atoms = Capability {
            name: "has_atoms".to_string(),
            version: "1.0".to_string(),
        };
        let has_bonds = Capability {
            name: "has_bonds".to_string(),
            version: "1.0".to_string(),
        };
        let has_coords_3d = Capability {
            name: "has_coordinates_3d".to_string(),
            version: "1.0".to_string(),
        };

        let caps = vec![has_atoms.clone(), has_bonds.clone()];
        let model = MockModel::new(&caps);

        assert!(model.has_capability(&has_atoms));
        assert!(model.has_capability(&has_bonds));
        assert!(!model.has_capability(&has_coords_3d));
    }

    #[test]
    fn test_model_capability_intersection() {
        let has_atoms = Capability {
            name: "has_atoms".to_string(),
            version: "1.0".to_string(),
        };
        let has_bonds = Capability {
            name: "has_bonds".to_string(),
            version: "1.0".to_string(),
        };
        let has_coords_3d = Capability {
            name: "has_coordinates_3d".to_string(),
            version: "1.0".to_string(),
        };

        let model1 = MockModel::new(&[has_atoms.clone(), has_bonds.clone()]);
        let model2 = MockModel::new(&[has_atoms.clone(), has_coords_3d.clone()]);

        let common = model1.capabilities().intersection(&model2.capabilities()).cloned().collect::<HashSet<_>>();
        assert_eq!(common.len(), 1);
        assert!(common.contains(&has_atoms));
    }

    #[test]
    fn test_model_conversion() {
        let has_atoms = Capability {
            name: "has_atoms".to_string(),
            version: "1.0".to_string(),
        };
        let caps = vec![has_atoms];
        let model = MockModel::new(&caps);
        let advanced = model.convert_to().unwrap();

        assert_eq!(model.capabilities(), advanced.capabilities());
    }

    #[test]
    fn test_instance_conversion() {
        let has_atoms = Capability {
            name: "has_atoms".to_string(),
            version: "1.0".to_string(),
        };
        let entity = MockEntity::new("test");
        let model = MockModel::new(&[has_atoms.clone()]);
        let instance = Instance::new(entity.clone(), model);

        let converted = instance.convert_to().unwrap();
        assert_eq!(converted.entity(), &entity);
        assert!(converted.model().has_capability(&has_atoms));
    }

    #[test]
    fn test_property_metadata() {
        assert_eq!(MockProperty::name(), "Mock Property");
        assert_eq!(MockProperty::description(), "A mock property for testing");
        assert_eq!(MockProperty::units(), Some("mock_units"));
    }

    #[test]
    fn test_property_computation() {
        let has_atoms = Capability {
            name: "has_atoms".to_string(),
            version: "1.0".to_string(),
        };
        let entity = MockEntity::new("test");
        let model = MockModel::new(&[has_atoms.clone()]);
        let instance = Instance::new(entity, model);

        let value = MockProperty::compute(&instance).unwrap();
        assert_eq!(value, 42.0);
    }

    #[test]
    fn test_property_missing_capability() {
        let has_atoms = Capability {
            name: "has_atoms".to_string(),
            version: "1.0".to_string(),
        };
        let entity = MockEntity::new("test");
        let model = MockModel::new(&[]);
        let instance = Instance::new(entity, model);

        let result = MockProperty::compute(&instance);
        assert!(matches!(result, 
            Err(Error::Property(PropertyError::MissingCapability(cap))) if cap == has_atoms
        ));
    }

    #[test]
    fn test_property_required_capabilities() {
        let has_atoms = Capability {
            name: "has_atoms".to_string(),
            version: "1.0".to_string(),
        };
        let caps = MockProperty::required_capabilities();
        assert_eq!(caps.len(), 1);
        assert!(caps.contains(&has_atoms));
    }
} 