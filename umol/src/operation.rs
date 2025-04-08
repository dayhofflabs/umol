//! Operation traits and types.
//!
//! This module provides traits for operations on stuff:
//! - Basic operations
//! - Conversion operations
//! - Operation composition

use crate::{ConvertTo, Stuff, Model, Result};

/// A trait for operations that transform stuff
pub trait Operation {
    /// The input model type
    type Input: Model;
    /// The output model type
    type Output: Model;

    /// Apply operation to stuff
    fn apply<
        SI: Stuff<Model = Self::Input>,
        SO: Stuff<Model = Self::Output, Entity = SI::Entity>,
    >(
        &self,
        stuff: &SI,
    ) -> Result<SO>
    where
        SI::Entity: Clone;
}

/// A conversion operation lifted to work on stuff
pub struct ConversionOperation<M1: Model, M2: Model>
where
    M1: ConvertTo<M2>,
{
    _phantom: std::marker::PhantomData<(M1, M2)>,
}

impl<M1: Model, M2: Model> ConversionOperation<M1, M2>
where
    M1: ConvertTo<M2>,
{
    /// Create a new conversion operation
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<M1: Model, M2: Model> Operation for ConversionOperation<M1, M2>
where
    M1: ConvertTo<M2>,
{
    type Input = M1;
    type Output = M2;

    fn apply<
        SI: Stuff<Model = Self::Input>,
        SO: Stuff<Model = Self::Output, Entity = SI::Entity>,
    >(
        &self,
        stuff: &SI,
    ) -> Result<SO>
    where
        SI::Entity: Clone,
    {
        let new_model = stuff.model().convert_to()?;
        let entity = stuff.entity().clone();
        SO::from_components(entity, new_model)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Capability, Entity};
    use serde::{Deserialize, Serialize};
    use std::collections::HashSet;

    // Test models
    #[derive(Debug, Clone, Serialize, Deserialize)]
    struct SourceModel {
        data: SourceData,
    }

    #[derive(Debug, Clone, Serialize, Deserialize)]
    struct SourceData {
        value: i32,
    }

    impl Model for SourceModel {
        type Data = SourceData;

        fn data(&self) -> &Self::Data {
            &self.data
        }

        fn capabilities(&self) -> HashSet<Capability> {
            let mut caps = HashSet::new();
            caps.insert(Capability::local("source", 1));
            caps
        }
    }

    #[derive(Debug, Clone, Serialize, Deserialize)]
    struct TargetModel {
        data: TargetData,
    }

    #[derive(Debug, Clone, Serialize, Deserialize)]
    struct TargetData {
        value: i32,
        processed: bool,
    }

    impl Model for TargetModel {
        type Data = TargetData;

        fn data(&self) -> &Self::Data {
            &self.data
        }

        fn capabilities(&self) -> HashSet<Capability> {
            let mut caps = HashSet::new();
            caps.insert(Capability::local("target", 1));
            caps
        }
    }

    // Test entity
    #[derive(Debug, Clone, Serialize, Deserialize)]
    struct TestEntity {
        id: String,
    }

    impl Entity for TestEntity {
        fn namespace(&self) -> Option<&str> {
            None
        }

        fn id(&self) -> &str {
            &self.id
        }

        fn label(&self) -> &str {
            &self.id
        }
    }

    // Test stuff
    #[derive(Debug, Clone)]
    struct SourceStuff {
        entity: TestEntity,
        model: SourceModel,
    }

    impl Stuff for SourceStuff {
        type Model = SourceModel;
        type Entity = TestEntity;

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

    #[derive(Debug, Clone)]
    struct TargetStuff {
        entity: TestEntity,
        model: TargetModel,
    }

    impl Stuff for TargetStuff {
        type Model = TargetModel;
        type Entity = TestEntity;

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

    // Implement conversion between models
    impl ConvertTo<TargetModel> for SourceModel {
        fn convert_to(&self) -> Result<TargetModel> {
            Ok(TargetModel {
                data: TargetData {
                    value: self.data.value,
                    processed: false,
                },
            })
        }
    }

    // Test operation
    struct TestOperation;

    impl Operation for TestOperation {
        type Input = SourceModel;
        type Output = TargetModel;

        fn apply<
            SI: Stuff<Model = Self::Input>,
            SO: Stuff<Model = Self::Output, Entity = SI::Entity>,
        >(
            &self,
            stuff: &SI,
        ) -> Result<SO>
        where
            SI::Entity: Clone,
        {
            let new_model = stuff.model().convert_to()?;
            let entity = stuff.entity().clone();
            SO::from_components(entity, new_model)
        }
    }

    #[test]
    fn test_conversion_operation() {
        let source_model = SourceModel {
            data: SourceData { value: 2 },
        };

        let entity = TestEntity {
            id: "test".to_string(),
        };

        let source_stuff = SourceStuff {
            entity: entity.clone(),
            model: source_model,
        };

        let operation = ConversionOperation::<SourceModel, TargetModel>::new();
        let target_stuff = operation
            .apply::<SourceStuff, TargetStuff>(&source_stuff)
            .unwrap();

        assert_eq!(target_stuff.entity().id(), "test");
        assert_eq!(target_stuff.model().data.value, 2);
        assert!(!target_stuff.model().data.processed);
    }

    #[test]
    fn test_conversion_operation_preserves_entity() {
        let source_model = SourceModel {
            data: SourceData { value: 42 },
        };

        let entity = TestEntity {
            id: "test".to_string(),
        };

        let source_stuff = SourceStuff {
            entity: entity.clone(),
            model: source_model,
        };

        let operation = ConversionOperation::<SourceModel, TargetModel>::new();
        let target_stuff = operation
            .apply::<SourceStuff, TargetStuff>(&source_stuff)
            .unwrap();

        // Verify the entity is preserved
        assert_eq!(target_stuff.entity().id(), source_stuff.entity().id());
        assert_eq!(target_stuff.entity().label(), source_stuff.entity().label());
    }

    #[test]
    fn test_conversion_operation_model_capabilities() {
        let source_model = SourceModel {
            data: SourceData { value: 3 },
        };

        let entity = TestEntity {
            id: "test".to_string(),
        };

        let source_stuff = SourceStuff {
            entity: entity.clone(),
            model: source_model,
        };

        let operation = ConversionOperation::<SourceModel, TargetModel>::new();
        let target_stuff = operation
            .apply::<SourceStuff, TargetStuff>(&source_stuff)
            .unwrap();

        // Verify target model has correct capabilities
        let caps = target_stuff.model().capabilities();
        assert!(caps.contains(&Capability::local("target", 1)));
        assert!(!caps.contains(&Capability::local("source", 1)));
    }

    #[test]
    fn test_custom_operation() {
        let source_model = SourceModel {
            data: SourceData { value: 3 },
        };

        let entity = TestEntity {
            id: "test".to_string(),
        };

        let source_stuff = SourceStuff {
            entity: entity.clone(),
            model: source_model,
        };

        let operation = TestOperation;
        let target_stuff = operation
            .apply::<SourceStuff, TargetStuff>(&source_stuff)
            .unwrap();

        assert_eq!(target_stuff.entity().id(), "test");
        assert_eq!(target_stuff.model().data.value, 3);
        assert!(!target_stuff.model().data.processed);
    }
}
