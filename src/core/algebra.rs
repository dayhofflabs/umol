//! Model algebra.
//!
//! Algebraic operations on models:
//! - Ensembles of chemical models (e.g., conformers, resonance structures,
//!   positive fractional coefficients)
//! - Aggregates of chemical models (e.g., mixtures, fixed-stoichiometry
//!   complexes, positive integer coefficients)
//! - Processes (e.g., reactions, positive and negative integer coefficients)

use crate::core::Model;

/// Represents an ensemble of models with weights (e.g., conformers, resonance
/// structures, positive fractional coefficients)
pub trait Ensemble<M: Model> {
    /// Returns a slice of model-weight pairs
    fn components(&self) -> &[(M, Option<f64>)];

    /// Returns a vector of weights
    fn weights(&self) -> Option<Vec<f64>> {
        self.components().iter().map(|(_, w)| *w).collect()
    }

    /// Returns a vector of references to models
    fn models(&self) -> Vec<&M> {
        self.components().iter().map(|(m, _)| m).collect()
    }
}

/// Represents an aggregate of models with coefficients (e.g., mixtures,
/// fixed-stoichiometry complexes, positive integer coefficients)
pub trait Aggregate<M: Model> {
    /// Returns a slice of model-coefficient pairs
    fn components(&self) -> &[(M, u32)];

    /// Returns a vector of coefficients
    fn coefficients(&self) -> Vec<u32> {
        self.components().iter().map(|(_, c)| *c).collect()
    }

    /// Returns a vector of references to models
    fn models(&self) -> Vec<&M> {
        self.components().iter().map(|(m, _)| m).collect()
    }
}

/// Represents a process of chemical models (e.g., reactions, positive and
/// negative integer coefficients)
pub trait Process<M: Model> {
    /// Returns a slice of model-coefficient pairs
    fn components(&self) -> &[(M, i32)];

    /// Returns a vector of coefficients
    fn coefficients(&self) -> Vec<i32> {
        self.components().iter().map(|(_, c)| *c).collect()
    }

    /// Returns a vector of references to models
    fn models(&self) -> Vec<&M> {
        self.components().iter().map(|(m, _)| m).collect()
    }
}
