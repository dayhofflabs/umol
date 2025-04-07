//! Property definitions and calculations.
//!
//! Properties are calculations that can be performed on models.
//! PropertySpec allows for type-safe property definitions.

use crate::core::{Capability, Model, Result};
use serde::{Deserialize, Serialize};
use std::collections::HashSet;

/// Base trait for all properties
pub trait Property<M: Model>: Serialize {
    /// The type of value this property computes
    type Value: Serialize + for<'de> Deserialize<'de>;
    /// The type of arguments this property accepts
    type Args;

    /// Get the name of this property
    fn name(&self) -> String;


    /// Get the description of this property
    fn description(&self) -> String;

    /// Get the units of this property, if applicable
    fn units(&self) -> Option<String>;

    /// Get the capabilities required to compute this property
    fn required_capabilities(&self) -> HashSet<Capability>;

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
    use crate::core::{AsModel, Entity, Model, Instance};
    use serde::{Deserialize, Serialize};
    use serde_json;

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
        pub count: usize,
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

    // Test instance
    #[derive(Debug, Clone)]
    struct SimpleInstance {
        entity: SimpleEntity,
        model: AtomCount,
    }

    impl Instance for SimpleInstance {
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

    /// A property that computes molecular mass for AtomCount
    #[derive(Debug, Clone, Serialize, Deserialize)]
    struct MolecularMass {
        name: String,
        description: String,
        units: Option<String>,
        required_capabilities: HashSet<Capability>,
        atomic_mass: f64,
    }

    impl MolecularMass {
        /// Create a new molecular mass property
        fn new(atomic_mass: f64) -> Self {
            let mut caps = HashSet::new();
            caps.insert(Capability::local("atom_count", 1));

            Self {
                name: "molecular_mass".to_string(),
                description: "Molecular mass in amu".to_string(),
                units: Some("amu".to_string()),
                required_capabilities: caps,
                atomic_mass,
            }
        }
    }

    impl Property<AtomCount> for MolecularMass {
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

        fn compute(&self, model: &AtomCount, _args: Self::Args) -> Result<Self::Value> {
            let data = model.data();
            Ok(data.count as f64 * self.atomic_mass)
        }
    }

    #[test]
    fn test_molecular_mass_computation() {
        // Create a test instance with 1 carbon atom
        let entity = SimpleEntity::local("C1", "carbon");
        let model = AtomCount::new(1);

        // Create the property with carbon's atomic mass
        let property = MolecularMass::new(12.011);

        // Test computation for 1 carbon atom
        let mass = property.compute(&model, ()).unwrap();
        assert!((mass - 12.011).abs() < 1e-6);  // Should be exactly 12.011 g/mol

        // Test computation for 2 carbon atoms
        let model = AtomCount::new(2);
        let instance = SimpleInstance::from_components(entity.clone(), model).unwrap();
        let mass = property.compute(instance.as_model(), ()).unwrap();
        assert!((mass - 24.022).abs() < 1e-6);  // Should be exactly 24.022 g/mol
    }

    #[test]
    fn test_molecular_mass_serialization() {
        // Create the property
        let property = MolecularMass::new(12.011);

        // Test serialization
        let json = serde_json::to_string(&property).unwrap();
        let expected_json = r#"{"name":"molecular_mass","description":"Molecular mass in amu","units":"amu","required_capabilities":[{"namespace":null,"name":"atom_count","version":1}],"atomic_mass":12.011}"#;
        assert_eq!(json, expected_json);

        // Test deserialization
        let deserialized: MolecularMass = serde_json::from_str(&json).unwrap();
        assert_eq!(deserialized.name, "molecular_mass");
        assert_eq!(deserialized.description, "Molecular mass in amu");
        assert_eq!(deserialized.units, Some("amu".to_string()));
        assert!(deserialized.required_capabilities.contains(&Capability::local("atom_count", 1)));
        assert!((deserialized.atomic_mass - 12.011).abs() < 1e-6);
    }
}
