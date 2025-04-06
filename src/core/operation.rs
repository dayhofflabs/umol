//! Operation traits and types.
//! 
//! This module provides traits for operations on instances:
//! - Basic operations
//! - Conversion operations
//! - Operation composition

use crate::core::{Model, Instance, Result};
use crate::core::conversion::ConvertTo;

/// A trait for operations that transform instances
pub trait Operation<M1: Model, M2: Model> {
    /// Check if this operation is valid for the given instance
    fn is_valid_for(&self, instance: &Instance<M1>) -> Result<bool>;

    /// Apply the operation to an instance
    fn apply(&self, instance: &Instance<M1>) -> Result<Instance<M2>>;
}

/// A conversion operation that can be applied to instances
pub struct ConversionOperation<M2: Model> {
    _phantom: std::marker::PhantomData<M2>,
}

impl<M2: Model> ConversionOperation<M2> {
    /// Create a new conversion operation
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<M1: Model, M2: Model> Operation<M1, M2> for ConversionOperation<M2>
where M1: ConvertTo<M2> {
    fn is_valid_for(&self, _instance: &Instance<M1>) -> Result<bool> {
        // A conversion operation is valid if the model can be converted
        Ok(true)
    }

    fn apply(&self, instance: &Instance<M1>) -> Result<Instance<M2>> {
        let new_model = instance.model().convert_to()?;
        Instance::new(instance.entity().clone(), new_model)
    }
} 