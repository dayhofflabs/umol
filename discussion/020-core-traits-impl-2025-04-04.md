# Core Traits Implementation

## Overview
Today we implemented the foundational traits for the molecular modeling framework based on the semantic model discussed in [15-semantic-model-2025-03-10.md](15-semantic-model-2025-03-10.md). The implementation focuses on establishing a clear separation of concerns between entities, models, properties, and their interactions.

## Major Changes

### 1. Core Module Restructuring
- Moved from `core/mod.rs` to `core.rs` for better organization
- Removed `AtomLink` and `AtomSite` traits as they didn't fit the new semantic model
- Created separate modules for each core concept:
  - `entity.rs`: Entity and Relation traits
  - `model.rs`: Model trait and capabilities
  - `algebra.rs`: Ensemble and Aggregate traits
  - `instance.rs`: Instance struct and Operation trait
  - `property.rs`: Property traits and types
  - `conversion.rs`: Conversion traits
  - `error.rs`: Error types
  - `testing.rs`: Testing infrastructure

### 2. Error Handling
Implemented a comprehensive error handling system using `thiserror`:

```rust
pub enum Error {
    Entity(EntityError),
    Model(ModelError),
    Conversion(ConversionError),
    Operation(OperationError),
    Property(PropertyError),
    ValidationError(String),
    Multiple(Vec<Error>),
    Other(Box<dyn std::error::Error + Send + Sync>),
}
```

Each domain has its own specific error type with detailed variants.

### 3. Model Capabilities
Defined an extensive set of capabilities that models can provide:

```rust
pub enum Capability {
    // Structural capabilities
    HasAtoms, HasBonds, HasAromaticity, HasStereochemistry,
    HasCharges, HasRadicals, HasIsotopes, HasResonance,
    
    // Geometric capabilities
    HasCoordinates2D, HasCoordinates3D,
    HasSymmetry, HasConformers,
    
    // Electronic capabilities
    HasElectronDensity, HasChargeDistribution,
    HasOrbitalEnergies, HasWavefunction, HasSpinDensity,
    
    // Property calculation capabilities
    CanComputeEnergy, CanComputeGradient,
    CanOptimizeGeometry, CanComputeCharges,
    CanComputeVibrationalModes,
    
    // Ensemble capabilities
    CanHandleEnsembles, CanHandleAggregates,
    
    // Transformation capabilities
    CanTransformStructure, CanPerformReactions,
}
```

### 4. Property Framework
Implemented a hierarchical property system:

```rust
pub trait Property {
    type Value;
    fn name() -> &'static str;
    fn description() -> &'static str;
    fn units() -> Option<&'static str>;
    fn required_capabilities() -> HashSet<Capability>;
    fn compute<E: Entity, M: Model>(instance: &Instance<E, M>) -> Result<Self::Value>;
}

pub trait MolecularProperty: Property {
    fn scope() -> PropertyScope;
    fn is_intensive() -> bool;
    fn is_extensive() -> bool;
}

pub trait EnergyProperty: MolecularProperty {
    fn energy_type() -> EnergyType;
    fn is_relative() -> bool;
}

pub trait StructuralProperty: MolecularProperty {
    fn features() -> HashSet<StructuralFeature>;
}
```

### 5. Testing Infrastructure
Created a testing framework for models and properties:

```rust
pub trait ModelTest {
    type E: Entity;
    type M: Model;
    
    fn create_test_instance() -> Result<Instance<Self::E, Self::M>>;
    fn test_capabilities() -> Result<()>;
    fn test_model_operations() -> Result<()>;
    fn test_property_calculations() -> Result<()>;
}

pub trait PropertyTest {
    type P: Property;
    
    fn test_requirements() -> Result<()>;
    fn test_simple_calculation() -> Result<()>;
    fn test_edge_cases() -> Result<()>;
}
```

## Next Steps

1. Fix remaining linter errors in the testing module
2. Implement the graph model using this foundation
3. Add tests for the core traits
4. Consider adding more specific capabilities or error types as needed
5. Enhance the property framework with more specific traits if required

## Open Questions

1. Should we add more specific error types for each domain?
2. Do we need additional capabilities for the graph model?
3. Should we expand the testing framework with more specific test cases?
4. Do we need to add more property types or features? 