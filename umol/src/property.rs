//! Property definitions and calculations.
//!
//! Properties are calculations that can be performed on models.
//! PropertySpec allows for type-safe property definitions.

use crate::{Capability, Model, Result};
use serde::{Deserialize, Serialize};
use std::collections::HashSet;

/// Base trait for all properties
pub trait Property<M: Model> {
    /// The type of value this property computes
    type Value: Serialize + for<'de> Deserialize<'de>;
    /// The type of arguments this property accepts
    type Args;

    /// Get the name of this property
    fn name(&self) -> String;

    /// Get the description of this property
    fn description(&self) -> String {
        self.name()
    }

    /// Get the units of this property, if applicable
    fn units(&self) -> Option<String> {
        None
    }

    /// Get the capabilities required to compute this property
    fn required_capabilities(&self) -> HashSet<Capability> {
        HashSet::new()
    }

    /// Compute the property for a given model and arguments
    fn compute(&self, model: &M, args: Self::Args) -> Result<Self::Value>;
}

/// Defines the input/output types and computation contract for a property
pub trait PropertySpec<M: Model> {
    /// The type of value this property computes
    type Value: Serialize + for<'de> Deserialize<'de>;
    /// The type of arguments this property accepts
    type Args;

    /// Compute the property value for a given model and arguments
    fn compute_spec(&self, model: &M, args: Self::Args) -> Result<Self::Value>;
}

// Implement PropertySpec for any Property
impl<M: Model, P: Property<M>> PropertySpec<M> for P {
    type Value = P::Value;
    type Args = P::Args;

    fn compute_spec(&self, model: &M, args: Self::Args) -> Result<Self::Value> {
        self.compute(model, args)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{AsModel, Entity, Instance, Model};
    use map_macro::hash_set;
    use serde::{Deserialize, Serialize};
    use umol_macros::property;

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
    struct SimpleModel {
        data: SimpleModelData,
    }

    #[derive(Debug, Clone, Serialize, Deserialize)]
    struct SimpleModelData {
        pub value: i32,
    }

    impl Model for SimpleModel {
        type Data = SimpleModelData;

        fn data(&self) -> &Self::Data {
            &self.data
        }

        fn capabilities(&self) -> HashSet<Capability> {
            hash_set! {
                Capability::local("simple_model", 1),
            }
        }
    }

    impl SimpleModel {
        fn new(value: i32) -> Self {
            Self {
                data: SimpleModelData { value },
            }
        }
    }

    // Test instance
    #[derive(Debug, Clone)]
    struct SimpleInstance {
        entity: SimpleEntity,
        model: SimpleModel,
    }

    impl Instance for SimpleInstance {
        type Entity = SimpleEntity;
        type Model = SimpleModel;

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

    // Test property
    struct SimpleProperty {
        name: String,
        description: String,
        units: Option<String>,
        required_capabilities: HashSet<Capability>,
        arg: f64,
    }

    impl SimpleProperty {
        fn new(arg: f64) -> Self {
            Self {
                name: "simple_property".to_string(),
                description: "Simple property".to_string(),
                units: Some("unit".to_string()),
                required_capabilities: hash_set! {
                    Capability::local("simple_model", 1),
                },
                arg,
            }
        }
    }

    impl Property<SimpleModel> for SimpleProperty {
        type Value = f64;
        type Args = ();

        fn name(&self) -> String {
            self.name.clone()
        }

        fn description(&self) -> String {
            self.description.clone()
        }

        fn units(&self) -> Option<String> {
            self.units.clone()
        }

        fn required_capabilities(&self) -> HashSet<Capability> {
            self.required_capabilities.clone()
        }

        fn compute(&self, model: &SimpleModel, _args: Self::Args) -> Result<Self::Value> {
            let data = model.data();
            Ok(data.value as f64 * self.arg)
        }
    }

    #[test]
    fn test_property_computation() {
        let entity = SimpleEntity::local("test", "test");
        let model = SimpleModel::new(1);
        let property = SimpleProperty::new(4.0);

        // Test computation with a model
        let value = property.compute(&model, ()).unwrap();
        assert!((value - 4.0).abs() < 1e-6);

        // Test computation with an instance
        let model = SimpleModel::new(2);
        let instance = SimpleInstance::from_components(entity.clone(), model).unwrap();
        let value = property.compute(instance.as_model(), ()).unwrap();
        assert!((value - 8.0).abs() < 1e-6);
    }

    struct SimplePropertyWithArgs;

    #[property(method = "property_with_args")]
    impl Property<SimpleModel> for SimplePropertyWithArgs {
        type Value = i32;
        type Args = i32;

        fn name(&self) -> String {
            "simple_property_with_args".to_string()
        }

        fn description(&self) -> String {
            "Simple property with arguments".to_string()
        }

        fn units(&self) -> Option<String> {
            Some("unit".to_string())
        }

        fn compute(&self, model: &SimpleModel, args: Self::Args) -> Result<Self::Value> {
            Ok(model.data().value * args)
        }
    }

    struct SimplePropertyNoArgs;

    #[property]
    impl Property<SimpleModel> for SimplePropertyNoArgs {
        type Value = i32;
        type Args = ();

        fn name(&self) -> String {
            "simple_property_no_args".to_string()
        }

        fn compute(&self, model: &SimpleModel, _args: Self::Args) -> Result<Self::Value> {
            Ok(model.data().value)
        }
    }

    #[test]
    fn test_property_macro_with_args() {
        let model = SimpleModel::new(2);
        let property = SimplePropertyWithArgs;

        assert_eq!(property.name(), "simple_property_with_args");
        assert_eq!(property.description(), "Simple property with arguments");
        assert_eq!(property.units(), Some("unit".to_string()));
        assert!(property.required_capabilities().is_empty());
        assert_eq!(property.compute(&model, 3).unwrap(), 6);
        assert_eq!(model.property_with_args(3).unwrap(), 6);
    }

    #[test]
    fn test_property_macro_no_args() {
        let model = SimpleModel::new(2);
        let property = SimplePropertyNoArgs;

        assert_eq!(property.name(), "simple_property_no_args");
        assert_eq!(property.description(), "simple_property_no_args");
        assert_eq!(property.units(), None);
        assert!(property.required_capabilities().is_empty());
        assert_eq!(property.compute(&model, ()).unwrap(), 2);
        assert_eq!(model.simple_property_no_args().unwrap(), 2);
    }
}
