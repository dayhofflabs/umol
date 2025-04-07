//! Operation traits and types.
//! 
//! This module provides traits for operations on instances:
//! - Basic operations
//! - Conversion operations
//! - Operation composition

use crate::core::{Model, Instance, Result};
use crate::core::conversion::ConvertTo;

/// A trait for operations that transform instances
pub trait Operation {
    /// The input model type
    type Input: Model;
    /// The output model type
    type Output: Model;

    /// Apply the operation to an instance
    fn apply<I: Instance<Model = Self::Input>, J: Instance<Model = Self::Output, Entity = I::Entity>>(
        &self,
        instance: &I,
    ) -> Result<J>
    where I::Entity: Clone;
}

/// A conversion operation that can be applied to instances
pub struct ConversionOperation<M1: Model, M2: Model> 
where M1: ConvertTo<M2> {
    _phantom: std::marker::PhantomData<(M1, M2)>,
}

impl<M1: Model, M2: Model> ConversionOperation<M1, M2> 
where M1: ConvertTo<M2> {
    /// Create a new conversion operation
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<M1: Model, M2: Model> Operation for ConversionOperation<M1, M2> 
where M1: ConvertTo<M2> {
    type Input = M1;
    type Output = M2;

    fn apply<I: Instance<Model = Self::Input>, J: Instance<Model = Self::Output, Entity = I::Entity>>(
        &self,
        instance: &I,
    ) -> Result<J>
    where I::Entity: Clone {
        let new_model = instance.model().convert_to()?;
        let entity = instance.entity().clone();
        J::from_components(entity, new_model)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde::{Serialize, Deserialize};
    use std::collections::HashSet;
    use crate::core::{Entity, Capability};

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

    // Test instances
    #[derive(Debug, Clone)]
    struct SourceInstance {
        entity: TestEntity,
        model: SourceModel,
    }

    impl Instance for SourceInstance {
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
    struct TargetInstance {
        entity: TestEntity,
        model: TargetModel,
    }

    impl Instance for TargetInstance {
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

        fn apply<I: Instance<Model = Self::Input>, J: Instance<Model = Self::Output, Entity = I::Entity>>(
            &self,
            instance: &I,
        ) -> Result<J>
        where I::Entity: Clone {
            let new_model = instance.model().convert_to()?;
            let entity = instance.entity().clone();
            J::from_components(entity, new_model)
        }
    }

    #[test]
    fn test_conversion_operation() {
        let source_model = SourceModel {
            data: SourceData { value: 42 },
        };
        
        let entity = TestEntity {
            id: "test".to_string(),
        };
        
        let source_instance = SourceInstance {
            entity: entity.clone(),
            model: source_model,
        };

        let operation = ConversionOperation::<SourceModel, TargetModel>::new();
        let target_instance = operation.apply::<SourceInstance, TargetInstance>(&source_instance).unwrap();

        assert_eq!(target_instance.entity().id(), "test");
        assert_eq!(target_instance.model().data.value, 42);
        assert!(!target_instance.model().data.processed);
    }

    #[test]
    fn test_conversion_operation_preserves_entity() {
        let source_model = SourceModel {
            data: SourceData { value: 42 },
        };
        
        let entity = TestEntity {
            id: "test".to_string(),
        };
        
        let source_instance = SourceInstance {
            entity: entity.clone(),
            model: source_model,
        };

        let operation = ConversionOperation::<SourceModel, TargetModel>::new();
        let target_instance = operation.apply::<SourceInstance, TargetInstance>(&source_instance).unwrap();

        // Verify the entity is preserved
        assert_eq!(target_instance.entity().id(), source_instance.entity().id());
        assert_eq!(target_instance.entity().label(), source_instance.entity().label());
    }

    #[test]
    fn test_conversion_operation_model_capabilities() {
        let source_model = SourceModel {
            data: SourceData { value: 42 },
        };
        
        let entity = TestEntity {
            id: "test".to_string(),
        };
        
        let source_instance = SourceInstance {
            entity: entity.clone(),
            model: source_model,
        };

        let operation = ConversionOperation::<SourceModel, TargetModel>::new();
        let target_instance = operation.apply::<SourceInstance, TargetInstance>(&source_instance).unwrap();

        // Verify target model has correct capabilities
        let caps = target_instance.model().capabilities();
        assert!(caps.contains(&Capability::local("target", 1)));
        assert!(!caps.contains(&Capability::local("source", 1)));
    }

    #[test]
    fn test_custom_operation() {
        let source_model = SourceModel {
            data: SourceData { value: 42 },
        };
        
        let entity = TestEntity {
            id: "test".to_string(),
        };
        
        let source_instance = SourceInstance {
            entity: entity.clone(),
            model: source_model,
        };

        let operation = TestOperation;
        let target_instance = operation.apply::<SourceInstance, TargetInstance>(&source_instance).unwrap();

        assert_eq!(target_instance.entity().id(), "test");
        assert_eq!(target_instance.model().data.value, 42);
        assert!(!target_instance.model().data.processed);
    }
} 