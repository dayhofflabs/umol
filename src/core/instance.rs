//! Instance types and operations.
//! 
//! Instances combine entities with their model representations:
//! - Instance creation and validation
//! - Model-specific operations
//! - Instance relationships and transformations
//! - Operation history tracking

use crate::core::{Entity, Model, ConvertTo, Result};

/// An instance pairs an entity with its representation in a specific model
pub struct Instance<E: Entity, M: Model> {
    pub entity: E,
    pub model: M,
}

impl<E: Entity + Clone, M: Model> Instance<E, M> {
    pub fn new(entity: E, model: M) -> Self {
        Self { entity, model }
    }

    pub fn entity(&self) -> &E {
        &self.entity
    }

    pub fn model(&self) -> &M {
        &self.model
    }

    /// Convert this instance to use a different model
    pub fn convert_to<M2: Model>(&self) -> Result<Instance<E, M2>> 
    where M: ConvertTo<M2> {
        let new_model = self.model.convert_to()?;
        Ok(Instance::new(self.entity.clone(), new_model))
    }
}

/// Operations are applied to instances
pub trait Operation<E: Entity, M1: Model, M2: Model> {
    fn apply(&self, instance: &Instance<E, M1>) -> Result<Instance<E, M2>>;
} 