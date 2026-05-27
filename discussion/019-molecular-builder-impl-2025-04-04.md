# Implementing the Molecular Builder

Date: April 4, 2025

## Overview

This document summarizes the implementation of the fluent molecular builder API and related components for the umol project. These changes enable intuitive, type-safe construction of molecular graphs while enforcing chemical rules.

## Core Components Implemented

### 1. Trait Hierarchy for Molecular Entities

We established a hierarchy of traits for different types of molecular entities:

```rust
/// Trait for any structure that can be added to a molecule
pub trait MolecularEntity {
    /// Add this entity to a builder, potentially using an attachment point
    fn add_to_builder(&self, builder: MoleculeBuilder, attachment: Option<AtomIndex>) 
        -> Result<MoleculeBuilder, Error>;
    
    /// Create a standalone molecule from this entity
    fn to_molecule(&self) -> Result<GraphMolecule, Error>;
}

/// Trait for fragments that can be connected to existing structures
pub trait Fragment: MolecularEntity {
    /// Get the primary attachment point for this fragment
    fn primary_attachment(&self) -> Option<AtomIndex>;
    
    /// Get all available attachment points 
    fn attachment_points(&self) -> Vec<(AtomIndex, String)>;
}

/// Trait for substructure queries
pub trait SubstructureQuery {
    /// Check if this query matches a given molecule
    fn matches(&self, molecule: &GraphMolecule) -> bool;
    
    /// Find all matches of this query in a molecule
    fn find_matches(&self, molecule: &GraphMolecule) -> Vec<Vec<AtomIndex>>;
}

/// Trait for molecular templates (parameterized structures)
pub trait MolecularTemplate {
    type Params;
    
    /// Generate a concrete molecular entity from template parameters
    fn instantiate(&self, params: Self::Params) -> Result<Box<dyn MolecularEntity>, Error>;
}
```

### 2. Fragment Registry

We implemented a registry system for molecular fragments to avoid hard-coding:

```rust
/// Registry for named fragments
pub struct FragmentRegistry {
    fragments: std::collections::HashMap<String, Box<dyn Fragment>>,
}

impl FragmentRegistry {
    pub fn new() -> Self { /* ... */ }
    pub fn register(&mut self, name: &str, fragment: Box<dyn Fragment>) { /* ... */ }
    pub fn get(&self, name: &str) -> Option<&dyn Fragment> { /* ... */ }
    pub fn names(&self) -> Vec<String> { /* ... */ }
    
    /// Parse a SMILES string into a fragment - placeholder until SMILES parser is implemented
    pub fn from_smiles(&self, smiles: &str) -> Result<Box<dyn Fragment>, Error> { /* ... */ }
}
```

### 3. Fluent Builder API

We implemented a fluent API for molecular construction using context objects to guide users through valid operations:

```rust
// Main builder
pub struct MoleculeBuilder {
    graph: StableGraph<GraphAtom, GraphBond>,
    active_atom: Option<AtomIndex>,
    operations: Vec<BuildOperation>,
    validation_set: ValidationSet,
}

// Example usage
let molecule = MoleculeBuilder::new()
    .atom(Element::C)
    .attach(Element::C, BondOrder::Single)?
    .attach(Element::O, BondOrder::Single)?
    .with_hydrogens(1)
    .done()
    .build()?;
```

Context types to support the fluent API:
- `AtomContext`: Operations after adding an atom
- `BondContext`: Operations for creating bonds 
- `BondResultContext`: Operations after creating a bond
- `SelectionContext`: Operations on selected atoms

### 4. Validation Framework

We implemented a comprehensive validation system:

```rust
/// Trait for any validation rule
pub trait ValidationRule {
    /// Check if the molecule satisfies this rule
    fn validate(&self, molecule: &GraphMolecule) -> Result<(), Vec<Error>>;
    
    /// Get a description of this rule
    fn description(&self) -> &str;
}

// Specific rules
pub struct ValenceRule;
pub struct ConnectivityRule;
pub struct ChemicalReasonablenessRule;

/// Collection of validation rules
pub struct ValidationSet {
    rules: Vec<Box<dyn ValidationRule>>,
}
```

### 5. Error Handling

We enhanced the Error enum to support validation errors and structural issues:

```rust
#[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
pub enum Error {
    // Existing errors
    
    #[error("Invalid operation: {0}")]
    InvalidOperation(String),
    
    #[error("Invalid valence for {element:?} at index {atom_idx:?}, expected {expected} bonds but has {actual}")]
    InvalidValence {
        element: Element,
        atom_idx: AtomIndex,
        expected: u8,
        actual: u8,
    },
    
    #[error("Molecule contains disconnected fragments")]
    DisconnectedMolecule,
    
    #[error("Chemically unreasonable structure: {0}")]
    UnreasonableStructure(String),
    
    #[error("Validation failed with multiple errors")]
    ValidationFailed(Vec<String>),
}
```

## Design Decisions

### Modal vs. Fluent API

We decided to implement a fluent API with implicit modes rather than an explicit state machine approach:
- Each context type (like `AtomContext`) represents a mode implicitly
- Method chaining guides the user through valid operations
- Type safety ensures correct operation sequencing
- The active atom concept maintains context between operations

### Fragment Handling

Instead of hard-coding fragments, we created a registry system:
- Fragments are defined dynamically rather than hard-coded in Rust
- SMILES parsing (to be implemented next) will provide a compact way to define fragments
- The registry enables lookup by name for common fragments

### Validation Approach

We chose a modular approach to validation:
- Validation rules are composable
- Rules can be applied selectively or as a standard set
- Validation can happen during building or at finalization
- All errors are collected and reported, not just the first one found

## Next Steps

1. **SMILES Parsing Implementation:**
   - Implement a SMILES parser to enable creating fragments from SMILES strings
   - Define common fragments via SMILES strings

2. **Fragment Library:**
   - Create a standard library of fragments (rings, functional groups) using SMILES
   - Implement a resource-based approach for loading these definitions

3. **Advanced Builder Features:**
   - Implement additional building operations (ring closures, etc.)
   - Add convenience methods for common patterns

4. **Testing:**
   - Create comprehensive tests for builder and validation
   - Test with complex molecular structures

5. **Documentation:**
   - Document the API with examples
   - Create usage guides for the builder

## Conclusion

The implemented fluent molecular builder API provides a intuitive, type-safe way to construct molecular graphs. The trait-based approach ensures extensibility, while the validation system enforces chemical correctness. The next priority is implementing SMILES parsing to fully realize the fragment-based approach to molecular construction.