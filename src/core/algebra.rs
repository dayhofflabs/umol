use crate::core::Model;

/// Represents an ensemble of models with weights (e.g., conformers, resonance structures)
pub trait Ensemble<M: Model> {
    fn components(&self) -> Vec<(M, f64)>;  // Models with weights
    
    fn weights(&self) -> Vec<f64> {
        self.components().into_iter().map(|(_, w)| w).collect()
    }
    
    fn models(&self) -> Vec<M> {
        self.components().into_iter().map(|(m, _)| m).collect()
    }
}

/// Represents an aggregate of models with coefficients (e.g., reaction systems)
pub trait Aggregate<M: Model> {
    fn components(&self) -> Vec<(M, f64)>;  // Models with coefficients
    
    fn coefficients(&self) -> Vec<f64> {
        self.components().into_iter().map(|(_, c)| c).collect()
    }
    
    fn models(&self) -> Vec<M> {
        self.components().into_iter().map(|(m, _)| m).collect()
    }
} 