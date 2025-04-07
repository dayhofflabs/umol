//! Instance types and operations.
//!
//! Instances combine entities with their model representations.

use crate::core::error::{OperationError, Result};
use crate::core::{ConvertTo, Entity, Model, Operation};

/// An instance pairs an entity with its representation in a specific model
pub struct Instance<M: Model> {
    /// The entity being represented
    entity: Entity,
    /// The model used for representation
    model: M,
}

impl<M: Model> Instance<M> {
    /// Create a new instance
    pub fn new(entity: Entity, model: M) -> Result<Self> {
        Ok(Self { entity, model })
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
    where
        M: ConvertTo<M2>,
    {
        let new_model = self.model.convert_to()?;
        Instance::new(self.entity.clone(), new_model)
    }

    /// Apply an operation to this instance
    pub fn apply<M2: Model, O: Operation<M, M2>>(&self, op: &O) -> Result<Instance<M2>> {
        // Validate that the operation is valid for this instance
        if !op.is_valid_for(self)? {
            return Err(OperationError::InvalidParameters(
                "Operation is not valid for this instance".into(),
            )
            .into());
        }
        op.apply(self)
    }
}
