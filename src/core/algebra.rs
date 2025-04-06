//! Algebraic operations on molecular objects.
//! 
//! This module provides algebraic structures for molecular modeling:
//! - Ensembles of molecular objects
//! - Aggregation operations
//! - Mathematical operations on models
//! - Transformations and symmetries

use crate::core::Model;

/// Represents an ensemble of models with weights (e.g., conformers, resonance structures)
pub trait Ensemble<M: Model> {
    /// Returns a slice of model-weight pairs
    fn components(&self) -> &[(M, f64)];
    
    /// Returns a vector of weights
    fn weights(&self) -> Vec<f64> {
        self.components().iter().map(|(_, w)| *w).collect()
    }
    
    /// Returns a vector of references to models
    fn models(&self) -> Vec<&M> {
        self.components().iter().map(|(m, _)| m).collect()
    }
}

/// Represents an aggregate of models with coefficients (e.g., reaction systems)
pub trait Aggregate<M: Model> {
    /// Returns a slice of model-coefficient pairs
    fn components(&self) -> &[(M, f64)];
    
    /// Returns a vector of coefficients
    fn coefficients(&self) -> Vec<f64> {
        self.components().iter().map(|(_, c)| *c).collect()
    }
    
    /// Returns a vector of references to models
    fn models(&self) -> Vec<&M> {
        self.components().iter().map(|(m, _)| m).collect()
    }
} 