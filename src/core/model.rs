use std::collections::HashSet;

/// Capabilities that models can provide
#[derive(Clone, Hash, Eq, PartialEq, Debug)]
pub enum Capability {
    // Structural capabilities
    HasAtoms,
    HasBonds,
    HasAromaticity,
    HasStereochemistry,
    HasCharges,
    HasRadicals,
    HasIsotopes,
    HasResonance,
    
    // Geometric capabilities
    HasCoordinates2D,
    HasCoordinates3D,
    HasSymmetry,
    HasConformers,
    
    // Electronic capabilities
    HasElectronDensity,
    HasChargeDistribution,
    HasOrbitalEnergies,
    HasWavefunction,
    HasSpinDensity,
    
    // Property calculation capabilities
    CanComputeEnergy,
    CanComputeGradient,
    CanOptimizeGeometry,
    CanComputeCharges,
    CanComputeVibrationalModes,
    
    // Ensemble capabilities
    CanHandleEnsembles,
    CanHandleAggregates,
    
    // Transformation capabilities
    CanTransformStructure,
    CanPerformReactions,
}

/// Base trait for all models
pub trait Model {
    /// Get the set of capabilities this model provides
    fn capabilities(&self) -> HashSet<Capability>;
    
    /// Check if the model has a specific capability
    fn has_capability(&self, cap: &Capability) -> bool {
        self.capabilities().contains(cap)
    }
    
    /// Check if the model has all the required capabilities
    fn has_capabilities(&self, caps: &[Capability]) -> bool {
        caps.iter().all(|cap| self.has_capability(cap))
    }
} 