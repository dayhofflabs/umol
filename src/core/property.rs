use std::collections::HashSet;
use crate::core::{Capability, Entity, Instance, Model, Result};

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
    fn compute<E: Entity, M: Model>(instance: &Instance<E, M>) -> Result<Self::Value>;
}

/// Properties that can be computed for molecules
pub trait MolecularProperty: Property {
    /// The scope of this property (atomic, molecular, etc.)
    fn scope() -> PropertyScope;
    
    /// Whether this property is intensive (independent of system size)
    fn is_intensive() -> bool;
    
    /// Whether this property is extensive (scales with system size)
    fn is_extensive() -> bool {
        !Self::is_intensive()
    }
}

/// Properties related to energy
pub trait EnergyProperty: MolecularProperty {
    /// The type of energy this property computes
    fn energy_type() -> EnergyType;
    
    /// Whether this is a relative or absolute energy
    fn is_relative() -> bool;
}

/// Properties related to structure
pub trait StructuralProperty: MolecularProperty {
    /// The structural features this property describes
    fn features() -> HashSet<StructuralFeature>;
}

/// The scope of a molecular property
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PropertyScope {
    Atomic,      // Property of individual atoms
    Bond,        // Property of bonds
    Fragment,    // Property of molecular fragments
    Molecular,   // Property of entire molecules
    Ensemble,    // Property of molecular ensembles
    System,      // Property of the entire chemical system
}

/// Types of energy properties
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnergyType {
    Electronic,
    Nuclear,
    Total,
    Kinetic,
    Potential,
    ZeroPoint,
    Thermal,
    Enthalpy,
    Entropy,
    Gibbs,
}

/// Structural features that properties can describe
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum StructuralFeature {
    Geometry,
    Topology,
    Electronics,
    Bonding,
    Aromaticity,
    Stereochemistry,
    Symmetry,
} 