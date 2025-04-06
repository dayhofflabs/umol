use crate::core::{
    error::{ModelError, PropertyError}, 
    Capability, ConvertTo, ConvertToWithMetadata, Error, Instance, Model, Property, Result,
    conversion::ConversionMetadata,
};
use std::collections::HashSet;

// Mock model implementations
#[derive(Debug, Clone)]
struct MockModel {
    capabilities: HashSet<Capability>,
    data: String,
}

impl MockModel {
    fn new(data: &str, capabilities: &[Capability]) -> Self {
        Self {
            data: data.to_string(),
            capabilities: capabilities.iter().cloned().collect(),
        }
    }
}

impl Model for MockModel {
    type Data = String;

    fn capabilities(&self) -> HashSet<Capability> {
        self.capabilities.clone()
    }

    fn data(&self) -> &Self::Data {
        &self.data
    }
}

#[derive(Debug, Clone)]
struct MockModelAdvanced {
    capabilities: HashSet<Capability>,
    data: String,
}

impl MockModelAdvanced {
    fn new(data: &str, capabilities: &[Capability]) -> Self {
        Self {
            data: data.to_string(),
            capabilities: capabilities.iter().cloned().collect(),
        }
    }
}

impl Model for MockModelAdvanced {
    type Data = String;

    fn capabilities(&self) -> HashSet<Capability> {
        self.capabilities.clone()
    }

    fn data(&self) -> &Self::Data {
        &self.data
    }
}

// Mock property implementation
struct MockProperty;

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

struct MockConversionParams {
    preserve_capabilities: bool,
}

impl Default for MockConversionParams {
    fn default() -> Self {
        Self {
            preserve_capabilities: true,
        }
    }
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
        metadata.attributes.insert("source".to_string(), "MockModel".to_string());

        let caps: Vec<Capability> = capabilities.into_iter().collect();
        Ok((MockModelAdvanced::new(&self.data, &caps), metadata))
    }
}

#[test]
fn test_model_capabilities() {
    let has_atoms = Capability::new("core", "has_atoms", 1);
    let has_bonds = Capability::new("core", "has_bonds", 1);
    let has_coords = Capability::new("core", "has_coordinates_3d", 1);

    let caps = vec![has_atoms.clone(), has_bonds.clone()];
    let model = MockModel::new("test", &caps);

    assert!(model.has_capability(&has_atoms));
    assert!(model.has_capability(&has_bonds));
    assert!(!model.has_capability(&has_coords));
}

#[test]
fn test_model_capability_intersection() {
    let has_atoms = Capability::new("core", "has_atoms", 1);
    let has_bonds = Capability::new("core", "has_bonds", 1);
    let has_coords = Capability::new("core", "has_coordinates_3d", 1);

    let model1 = MockModel::new("test", &[has_atoms.clone(), has_bonds.clone()]);
    let model2 = MockModel::new("test", &[has_atoms.clone(), has_coords.clone()]);

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
    let caps = vec![has_atoms.clone()];
    let model = MockModel::new("test", &caps);
    let advanced = model.convert_to().unwrap();

    assert_eq!(model.capabilities(), advanced.capabilities());
}

#[test]
fn test_instance_conversion() {
    use crate::core::Entity;
    let has_atoms = Capability::new("core", "has_atoms", 1);
    let entity = Entity::new("test", "Test Entity", None);
    let model = MockModel::new("test", &[has_atoms.clone()]);
    let instance: Instance<MockModel> = Instance::new(entity.clone(), model).unwrap();

    let converted: Instance<MockModelAdvanced> = instance.convert_to().unwrap();
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
    use crate::core::Entity;
    let has_atoms = Capability::new("core", "has_atoms", 1);
    let entity = Entity::new("test", "Test Entity", None);
    let model = MockModel::new("test", &[has_atoms.clone()]);
    let instance = Instance::new(entity, model).unwrap();

    let value = MockProperty::compute(&instance).unwrap();
    assert_eq!(value, 42.0);
}

#[test]
fn test_property_missing_capability() {
    use crate::core::Entity;
    let has_atoms = Capability::new("core", "has_atoms", 1);
    let entity = Entity::new("test", "Test Entity", None);
    let model = MockModel::new("test", &[]);
    let instance = Instance::new(entity, model).unwrap();

    let result = MockProperty::compute(&instance);
    assert!(matches!(result,
        Err(Error::Property(PropertyError::MissingCapability(cap))) if cap == has_atoms
    ));
}

#[test]
fn test_property_required_capabilities() {
    let has_atoms = Capability::new("core", "has_atoms", 1);
    let caps = MockProperty::required_capabilities();
    assert_eq!(caps.len(), 1);
    assert!(caps.contains(&has_atoms));
}

#[test]
fn test_model_validation() {
    struct ValidatingModel {
        capabilities: HashSet<Capability>,
        data: String,
        should_fail: bool,
    }

    impl ValidatingModel {
        fn new(data: &str, capabilities: &[Capability], should_fail: bool) -> Self {
            Self {
                data: data.to_string(),
                capabilities: capabilities.iter().cloned().collect(),
                should_fail,
            }
        }
    }

    impl Model for ValidatingModel {
        type Data = String;

        fn capabilities(&self) -> HashSet<Capability> {
            self.capabilities.clone()
        }

        fn data(&self) -> &Self::Data {
            &self.data
        }

        fn validate(&self) -> Result<()> {
            if self.should_fail {
                Err(Error::Model(ModelError::NotFound("Test validation error".into())))
            } else {
                Ok(())
            }
        }
    }

    // Test successful validation
    let model = ValidatingModel::new("test", &[], false);
    assert!(model.validate().is_ok());

    // Test failed validation
    let model = ValidatingModel::new("test", &[], true);
    let result = model.validate();
    assert!(result.is_err());
    if let Err(Error::Model(ModelError::NotFound(_))) = result {
        // Expected error
    } else {
        panic!("Unexpected error: {:?}", result);
    }
}

#[test]
fn test_complex_capability_requirements() {
    // Define a property that requires multiple capabilities
    struct ComplexProperty;

    impl Property for ComplexProperty {
        type Value = f64;

        fn name() -> &'static str {
            "Complex Property"
        }

        fn description() -> &'static str {
            "A property requiring multiple capabilities"
        }

        fn units() -> Option<&'static str> {
            Some("complex_units")
        }

        fn required_capabilities() -> HashSet<Capability> {
            let mut caps = HashSet::new();
            caps.insert(Capability::new("core", "has_atoms", 1));
            caps.insert(Capability::new("core", "has_bonds", 1));
            caps.insert(Capability::new("core", "has_coordinates_3d", 1));
            caps
        }

        fn compute<M: Model>(instance: &Instance<M>) -> Result<Self::Value> {
            let required = Self::required_capabilities();
            
            // Check all required capabilities
            for cap in required.iter() {
                if !instance.model().has_capability(cap) {
                    return Err(Error::Property(PropertyError::MissingCapability(cap.clone())));
                }
            }
            
            Ok(42.0)
        }
    }

    use crate::core::Entity;

    // Test with all required capabilities
    let all_caps = vec![
        Capability::new("core", "has_atoms", 1),
        Capability::new("core", "has_bonds", 1),
        Capability::new("core", "has_coordinates_3d", 1),
    ];
    let model = MockModel::new("test", &all_caps);
    let instance = Instance::new(Entity::new("test", "test", None), model).unwrap();
    assert!(ComplexProperty::compute(&instance).is_ok());

    // Test with missing capability
    let partial_caps = vec![
        Capability::new("core", "has_atoms", 1),
        Capability::new("core", "has_bonds", 1),
    ];
    let model = MockModel::new("test", &partial_caps);
    let instance = Instance::new(Entity::new("test", "test", None), model).unwrap();
    let result = ComplexProperty::compute(&instance);
    assert!(result.is_err());
    if let Err(Error::Property(PropertyError::MissingCapability(cap))) = result {
        assert_eq!(cap.name, "has_coordinates_3d");
    } else {
        panic!("Unexpected error: {:?}", result);
    }
}

#[test]
fn test_conversion_metadata() {
    use crate::core::conversion::ConversionMetadata;

    // Test basic conversion with metadata
    let has_atoms = Capability::new("core", "has_atoms", 1);
    let caps = vec![has_atoms.clone()];
    let model = MockModel::new("test", &caps);

    let mut metadata = ConversionMetadata::default();
    metadata.attributes.insert("source".to_string(), "test".to_string());
    metadata.attributes.insert("target".to_string(), "advanced".to_string());

    let params = MockConversionParams {
        preserve_capabilities: true,
    };

    let (advanced, result_metadata) = model.convert_to_with_metadata(&params).unwrap();

    // Verify model conversion was successful
    assert_eq!(model.capabilities(), advanced.capabilities());

    // Verify metadata was preserved
    assert!(result_metadata.attributes.contains_key("source"));
    assert_eq!(result_metadata.attributes.get("source").unwrap(), "MockModel");
}
