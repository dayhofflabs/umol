//! Property definitions and calculations.
//! 
//! Properties are calculations that can be performed on models:
//! - Property definitions with metadata
//! - Capability requirements for calculations
//! - Computation results and error handling
//! - Property relationships and dependencies

use std::collections::HashSet;
use crate::core::{Capability, Instance, Model, Result};

/// Base trait for all properties
pub trait Property {
    /// The type of value this property computes
    type Value;
    
    /// Get the name of this property
    fn name() -> &'static str where Self: Sized;
    
    /// Get a description of this property
    fn description() -> &'static str where Self: Sized;
    
    /// Get the units of this property, if applicable
    fn units() -> Option<&'static str> where Self: Sized;
    
    /// Get the capabilities required to compute this property
    fn required_capabilities() -> HashSet<Capability> where Self: Sized;
    
    /// Compute the property for a given instance
    fn compute<M: Model>(instance: &Instance<M>) -> Result<Self::Value>;
}