//! Instance types and operations.
//! 
//! Instances combine entities with their model representations:
//! - Instance creation and validation
//! - Model-specific operations
//! - Instance relationships and transformations
//! - Operation history tracking

use crate::core::{Model, Result, Entity};
use crate::core::conversion::ConvertTo;
use crate::core::operation::Operation;

/// An instance pairs an entity with its representation in a specific model
pub struct Instance<M: Model> {
    /// The entity being represented
    pub entity: Entity,
    /// The model used for representation
    pub model: M,
}

impl<M: Model> Instance<M> {
    /// Create a new instance
    pub fn new(entity: Entity, model: M) -> Self {
        Self { entity, model }
    }

    /// Get a reference to the entity
    pub fn entity(&self) -> &Entity {
        &self.entity
    }

    /// Get a reference to the model
    pub fn model(&self) -> &M {
        &self.model
    }

    /// Convert this instance to use a different model
    pub fn convert_to<M2: Model>(&self) -> Result<Instance<M2>> 
    where M: ConvertTo<M2> {
        let new_model = self.model.convert_to()?;
        Ok(Instance::new(self.entity.clone(), new_model))
    }

    /// Apply an operation to this instance
    pub fn apply<M2: Model, O: Operation<M, M2>>(
        &self,
        op: &O
    ) -> Result<Instance<M2>> {
        op.apply(self)
    }

    /// Validate that the entity type is compatible with the model
    pub fn validate(&self) -> Result<()> {
        // Validate the entity itself
        self.entity.validate()?;

        Ok(())
    }
} 