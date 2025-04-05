//! Entity traits and relationships.
//! 
//! Entities represent the semantic objects in molecular modeling:
//! - Abstract entities (atoms, bonds, molecules)
//! - Entity relationships (generalizes, specializes)
//! - Identity and equality semantics

use crate::core::Error;

/// Represents entities in the chemical domain (structures, conformers, etc.)
/// Forms a partial ordering through generalization/specialization relationships
pub trait Entity {
    fn generalizes(&self, other: &Self) -> bool;
    fn specializes(&self, other: &Self) -> bool;
}

/// Represents relationships between entities (transformations, reactions, etc.)
pub trait Relation {
    type Source: Entity;
    type Target: Entity;
    
    fn apply(&self, source: &Self::Source) -> Result<Self::Target, Error>;
    fn compose(&self, other: &Self) -> Option<Self> where Self: Sized;
} 