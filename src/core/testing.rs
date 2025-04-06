//! Testing utilities and traits.
//!
//! This module provides testing infrastructure for umol:
//! - Model testing traits
//! - Property testing traits
//! - Test utilities and helpers

use crate::core::{
    error::{ModelError, PropertyError},
    Capability, ConversionMetadata, ConvertToWithMetadata, Error, Instance, Model, Property,
    Result,
};
use std::collections::HashSet;

/// Trait for testing model implementations
pub trait ModelTest {
    /// The model type being tested
    type M: Model;

    /// Create a test instance for testing
    fn create_test_instance() -> Result<Instance<Self::M>>;

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
    pub fn verify_capabilities<M: Model>(model: &M, required: &[Capability]) -> Result<()> {
        for cap in required {
            if !model.has_capability(cap) {
                return Err(ModelError::MissingCapability(cap.clone()).into());
            }
        }
        Ok(())
    }

    /// Verify that a property can be computed on a model
    pub fn verify_property_computation<P, M>(instance: &Instance<M>) -> Result<()>
    where
        P: Property,
        M: Model,
    {
        P::compute(instance)?;
        Ok(())
    }
}

/// Default implementation of ModelTest
pub struct DefaultModelTest<M: Model> {
    _phantom: std::marker::PhantomData<M>,
}

impl<M: Model> DefaultModelTest<M> {
    /// Create a new default model test
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<M: Model> ModelTest for DefaultModelTest<M> {
    type M = M;

    fn create_test_instance() -> Result<Instance<Self::M>> {
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

/// A mock model for testing
pub struct MockModel {
    capabilities: HashSet<Capability>,
    data: (), // Empty data for testing
}

impl MockModel {
    pub fn new(capabilities: &[Capability]) -> Self {
        Self {
            capabilities: capabilities.iter().cloned().collect(),
            data: (),
        }
    }
}

impl Model for MockModel {
    type Data = ();

    fn data(&self) -> &Self::Data {
        &self.data
    }

    fn capabilities(&self) -> HashSet<Capability> {
        self.capabilities.clone()
    }
}

#[derive(Debug, Clone)]
pub struct MockModelAdvanced {
    capabilities: HashSet<Capability>,
    data: (), // Empty data for testing
}

impl MockModelAdvanced {
    pub fn new(capabilities: &[Capability]) -> Self {
        Self {
            capabilities: capabilities.iter().cloned().collect(),
            data: (),
        }
    }
}

impl Model for MockModelAdvanced {
    type Data = ();

    fn data(&self) -> &Self::Data {
        &self.data
    }

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
        params: &Self::Params,
    ) -> Result<(MockModelAdvanced, ConversionMetadata)> {
        let capabilities = if params.preserve_capabilities {
            self.capabilities.clone()
        } else {
            HashSet::new()
        };

        let mut metadata = ConversionMetadata::default();
        metadata
            .attributes
            .insert("source".to_string(), "MockModel".to_string());

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
        caps.insert(Capability::new("core", "has_atoms", 1));
        caps
    }

    fn compute<M: Model>(instance: &Instance<M>) -> Result<Self::Value> {
        let has_atoms = Capability::new("core", "has_atoms", 1);
        if instance.model().has_capability(&has_atoms) {
            Ok(42.0)
        } else {
            Err(Error::Property(PropertyError::MissingCapability(has_atoms)))
        }
    }
}

/// Helper function to verify property computation
pub fn verify_property_calculation<P: Property>(
    instance: &Instance<impl Model>,
    expected: P::Value,
) -> Result<()>
where
    P::Value: PartialEq + std::fmt::Display,
{
    let result = P::compute(instance)?;
    if result == expected {
        Ok(())
    } else {
        Err(PropertyError::CalculationFailed(format!(
            "Property computation result mismatch: expected {}, got {}",
            expected, result
        )).into())
    }
}

/// Test that verifies a model has the required capabilities
pub fn test_model_capabilities(model: &impl Model) -> Result<()> {
    let mut caps = HashSet::new();
    caps.insert(Capability::new("core", "has_atoms", 1));

    for cap in &caps {
        if !model.has_capability(cap) {
            return Err(ModelError::MissingCapability(cap.clone()).into());
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::Entity;

    #[test]
    fn test_model_capabilities() {
        let has_atoms = Capability::new("core", "has_atoms", 1);
        let has_bonds = Capability::new("core", "has_bonds", 1);
        let has_coords_3d = Capability::new("core", "has_coordinates_3d", 1);

        let caps = vec![has_atoms.clone(), has_bonds.clone()];
        let model = MockModel::new(&caps);

        assert!(model.has_capability(&has_atoms));
        assert!(model.has_capability(&has_bonds));
        assert!(!model.has_capability(&has_coords_3d));
    }

    #[test]
    fn test_model_capability_intersection() {
        let has_atoms = Capability::new("core", "has_atoms", 1);
        let has_bonds = Capability::new("core", "has_bonds", 1);
        let has_coords_3d = Capability::new("core", "has_coordinates_3d", 1);

        let model1 = MockModel::new(&[has_atoms.clone(), has_bonds.clone()]);
        let model2 = MockModel::new(&[has_atoms.clone(), has_coords_3d.clone()]);

        let common = model1
            .capabilities()
            .intersection(&model2.capabilities())
            .cloned()
            .collect::<HashSet<_>>();
        assert_eq!(common.len(), 1);
        assert!(common.contains(&has_atoms));
    }

    #[test]
    fn test_model_conversion() {
        let has_atoms = Capability::new("core", "has_atoms", 1);
        let caps = vec![has_atoms];
        let model = MockModel::new(&caps);
        let params = MockConversionParams {
            preserve_capabilities: true,
        };
        let (advanced, _) = model.convert_to_with_metadata(&params).unwrap();

        assert_eq!(model.capabilities(), advanced.capabilities());
    }

    #[test]
    fn test_instance_creation() {
        let has_atoms = Capability::new("core", "has_atoms", 1);
        let entity = Entity::new("test", "test", None);
        let model = MockModel::new(&[has_atoms.clone()]);
        let instance = Instance::new(entity, model).unwrap();

        assert!(instance.model().has_capability(&has_atoms));
        assert_eq!(instance.entity().id, "test");
    }

    #[test]
    fn test_instance_validation() {
        let has_atoms = Capability::new("core", "has_atoms", 1);
        let entity = Entity::new("test", "test", None);
        let model = MockModel::new(&[has_atoms.clone()]);
        let instance = Instance::new(entity, model).unwrap();

        assert!(instance.model().has_capability(&has_atoms));
    }

    #[test]
    fn test_property_computation() {
        let has_atoms = Capability::new("core", "has_atoms", 1);
        let entity = Entity::new("test", "test", None);
        let model = MockModel::new(&[has_atoms.clone()]);
        let instance = Instance::new(entity, model).unwrap();

        let result = MockProperty::compute(&instance);
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), 42.0);
    }

    #[test]
    fn test_model_conversion_with_metadata() {
        let has_atoms = Capability::new("core", "has_atoms", 1);
        let model = MockModel::new(&[has_atoms.clone()]);

        let params = MockConversionParams {
            preserve_capabilities: true,
        };

        let (advanced, metadata) = model.convert_to_with_metadata(&params).unwrap();

        assert_eq!(advanced.capabilities(), model.capabilities());
        assert_eq!(metadata.attributes.get("source").unwrap(), "MockModel");
    }

    #[test]
    fn test_model_conversion_without_metadata() {
        let has_atoms = Capability::new("core", "has_atoms", 1);
        let model = MockModel::new(&[has_atoms.clone()]);

        let params = MockConversionParams {
            preserve_capabilities: false,
        };

        let (advanced, _) = model.convert_to_with_metadata(&params).unwrap();

        assert!(advanced.capabilities().is_empty());
    }

    #[test]
    fn test_property_required_capabilities_validation() {
        // let has_atoms = Capability::new("core", "has_atoms", 1);
        let entity = Entity::new("test", "test", None);
        let model = MockModel::new(&[]); // No capabilities
        let instance = Instance::new(entity, model).unwrap();

        let result = MockProperty::compute(&instance);
        assert!(matches!(
            result,
            Err(Error::Property(PropertyError::MissingCapability(_)))
        ));
    }

    #[test]
    fn test_model_test_trait() {
        struct TestModelTest;
        impl ModelTest for TestModelTest {
            type M = MockModel;

            fn create_test_instance() -> Result<Instance<Self::M>> {
                let entity = Entity::new("test", "test", None);
                let model = MockModel::new(&[]);
                Instance::new(entity, model)
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

        let instance = TestModelTest::create_test_instance().unwrap();
        assert_eq!(instance.model().capabilities().len(), 0);
    }

    #[test]
    fn test_property_test_trait() {
        struct TestPropertyTest;
        impl PropertyTest for TestPropertyTest {
            type P = MockProperty;

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

        assert!(TestPropertyTest::test_requirements().is_ok());
    }
}
