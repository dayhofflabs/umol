use std::collections::HashSet;
use crate::core::{
    Capability, Entity, Model, Property,
    Instance, Operation, Result,
};

/// Trait for testing model implementations
pub trait ModelTest {
    type E: Entity;
    type M: Model;
    
    /// Create a new instance for testing
    fn create_test_instance() -> Result<Instance<Self::E, Self::M>>;
    
    /// Test that the model correctly reports its capabilities
    fn test_capabilities() -> Result<()>;
    
    /// Test basic model operations
    fn test_model_operations() -> Result<()>;
    
    /// Test property calculations
    fn test_property_calculations() -> Result<()>;
}

/// Trait for testing property implementations
pub trait PropertyTest {
    type P: Property;
    
    /// Test that the property correctly reports its requirements
    fn test_requirements() -> Result<()>;
    
    /// Test property calculation on a simple case
    fn test_simple_calculation() -> Result<()>;
    
    /// Test property calculation on edge cases
    fn test_edge_cases() -> Result<()>;
}

/// Helper function to verify capability requirements
pub fn verify_capabilities<M: Model>(model: &M, required: &[Capability]) -> Result<()> {
    let required: HashSet<_> = required.iter().cloned().collect();
    let available = model.capabilities();
    
    if !required.is_subset(&available) {
        let missing: Vec<_> = required.difference(&available).collect();
        return Err(crate::core::Error::Model(
            crate::core::ModelError::MissingCapability(
                missing.first().cloned().unwrap()
            )
        ));
    }
    
    Ok(())
}

/// Helper function to verify property calculation
pub fn verify_property_calculation<P: Property, E: Entity, M: Model>(
    instance: &Instance<E, M>
) -> Result<()> {
    // Verify the model has required capabilities
    verify_capabilities(instance.model(), &P::required_capabilities().into_iter().collect::<Vec<_>>())?;
    
    // Try to compute the property
    P::compute(instance)?;
    
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    
    // Add tests for the testing utilities themselves
    #[test]
    fn test_verify_capabilities() {
        // TODO: Add tests
    }
    
    #[test]
    fn test_verify_property_calculation() {
        // TODO: Add tests
    }
} 