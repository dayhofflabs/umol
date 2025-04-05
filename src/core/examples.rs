//! Example implementations of the plugin system
//! 
//! This module demonstrates how to implement a plugin for umol using
//! a quantum chemistry plugin as an example. It shows:
//! 
//! - How to create a plugin
//! - How to specify plugin dependencies
//! - How to register capabilities
//! - How to implement and register properties
//! - How to implement and register models
//! - How to implement model conversions
//! - How to implement file format handlers

use std::collections::HashSet;
use semver::Version;
use crate::core::{
    Plugin, Registry, Result, Error,
    Capability, PluginRequirements,
    PropertyDefinition, ModelDefinition,
    ConversionDefinition, FormatHandler,
};

/// Example quantum chemistry plugin that provides quantum mechanical
/// calculations and models.
/// 
/// This plugin demonstrates:
/// - Dependency management (requires core plugin)
/// - Capability registration
/// - Property implementation (electronic energy)
/// - Model implementation (wavefunction)
/// - Model conversion (wavefunction to density)
/// - File format support (Molden format)
pub struct QuantumPlugin;

impl Plugin for QuantumPlugin {
    fn name(&self) -> &str { "quantum" }
    fn version(&self) -> Version { "1.0.0".parse().unwrap() }
    
    fn requires(&self) -> PluginRequirements {
        PluginRequirements {
            plugins: [
                ("core".to_string(), "1.0.0".parse().unwrap())
            ].into_iter().collect(),
            capabilities: [
                Capability::new("core", "has_atoms"),
                Capability::new("core", "has_coordinates")
            ].into_iter().collect(),
        }
    }
    
    fn register(&self, registry: &mut Registry) {
        // Register quantum-specific capabilities
        registry.register_capability(Capability::new("quantum", "has_wavefunction"));
        registry.register_capability(Capability::new("quantum", "has_density_matrix"));
        
        // Register property for calculating electronic energy
        registry.register_property(
            "electronic_energy".to_string(),
            || Ok(Box::new(ElectronicEnergyProperty))
        );
        
        // Register wavefunction model
        registry.register_model(
            "wavefunction".to_string(),
            || Ok(Box::new(WavefunctionModel::new()))
        );
        
        // Register conversion from wavefunction to density
        registry.register_conversion(
            "wavefunction".to_string(),
            "density".to_string(),
            || Ok(Box::new(WavefunctionToDensityConversion))
        );
        
        // Register Molden file format support
        registry.register_format(
            "molden".to_string(),
            || Ok(Box::new(MoldenFormatHandler))
        );
    }
}

/// Property that calculates the electronic energy of a quantum system.
/// 
/// This demonstrates how to implement a property that:
/// - Requires specific capabilities (has_wavefunction)
/// - Performs a calculation on a model
/// - Returns a typed result
struct ElectronicEnergyProperty;

impl PropertyDefinition for ElectronicEnergyProperty {
    // Implementation details...
}

/// Model representing a quantum mechanical wavefunction.
/// 
/// This demonstrates how to implement a model that:
/// - Provides specific capabilities
/// - Stores quantum mechanical data
/// - Supports property calculations
struct WavefunctionModel {
    // Implementation details...
}

impl WavefunctionModel {
    fn new() -> Self {
        Self { /* ... */ }
    }
}

impl ModelDefinition for WavefunctionModel {
    // Implementation details...
}

/// Conversion from wavefunction to density matrix representation.
/// 
/// This demonstrates how to implement a model conversion that:
/// - Preserves essential information
/// - Handles conversion parameters
/// - Provides conversion metadata
struct WavefunctionToDensityConversion;

impl ConversionDefinition for WavefunctionToDensityConversion {
    // Implementation details...
}

/// Handler for the Molden file format.
/// 
/// This demonstrates how to implement a file format handler that:
/// - Reads molecular data from files
/// - Writes molecular data to files
/// - Validates file contents
struct MoldenFormatHandler;

impl FormatHandler for MoldenFormatHandler {
    // Implementation details...
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Test that demonstrates the plugin registration process:
    /// 1. Register core plugin (provides basic capabilities)
    /// 2. Register quantum plugin (depends on core)
    /// 3. Verify capabilities are available
    /// 4. Verify property is registered but not initialized
    #[test]
    fn test_plugin_registration() {
        let mut registry = Registry::new();
        
        // Register core plugin first
        let core = CorePlugin;
        registry.register_plugin(Box::new(core)).unwrap();
        
        // Register quantum plugin
        let quantum = QuantumPlugin;
        registry.register_plugin(Box::new(quantum)).unwrap();
        
        // Verify capabilities were registered
        assert!(registry.has_capability(&Capability::new("quantum", "has_wavefunction")));
        assert!(registry.has_capability(&Capability::new("quantum", "has_density_matrix")));
        
        // Verify property is available but not initialized
        let property = registry.get_property("electronic_energy");
        assert!(property.is_ok());
    }
    
    /// Test that demonstrates dependency checking:
    /// - Attempting to register quantum plugin without core plugin
    /// - Verifying appropriate error is returned
    #[test]
    fn test_missing_dependency() {
        let mut registry = Registry::new();
        
        // Try to register quantum plugin without core plugin
        let quantum = QuantumPlugin;
        let result = registry.register_plugin(Box::new(quantum));
        
        assert!(matches!(result, 
            Err(Error::MissingPlugin(p)) if p == "core"
        ));
    }
} 