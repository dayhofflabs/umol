// Molecular fragments and substructures

use crate::error::Error;
use crate::graph::{GraphMolecule, AtomIndex};
use crate::graph::builder::MoleculeBuilder;
use crate::core::types::{AtomIndex, BondIndex};
use super::{Atom, Bond, Molecule};

/// Trait for any structure that can be added to a molecule
pub trait MolecularEntity {
    /// Add this entity to a molecule builder, potentially using an attachment point
    fn add_to_builder(&self, builder: MoleculeBuilder, attachment: Option<AtomIndex>) -> Result<MoleculeBuilder, Error>;
    
    /// Create a standalone molecule from this entity
    fn to_molecule(&self) -> Result<GraphMolecule, Error> {
        let builder = MoleculeBuilder::new();
        self.add_to_builder(builder, None)?.build()
    }
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

/// Registry for named fragments
pub struct FragmentRegistry {
    fragments: std::collections::HashMap<String, Box<dyn Fragment>>,
}

impl FragmentRegistry {
    pub fn new() -> Self {
        Self {
            fragments: std::collections::HashMap::new(),
        }
    }
    
    pub fn register(&mut self, name: &str, fragment: Box<dyn Fragment>) {
        self.fragments.insert(name.to_string(), fragment);
    }
    
    pub fn get(&self, name: &str) -> Option<&dyn Fragment> {
        self.fragments.get(name).map(|b| b.as_ref())
    }
    
    pub fn names(&self) -> Vec<String> {
        self.fragments.keys().cloned().collect()
    }
    
    /// Parse a SMILES string into a fragment - placeholder until SMILES parser is implemented
    pub fn from_smiles(&self, _smiles: &str) -> Result<Box<dyn Fragment>, Error> {
        // Placeholder until we implement SMILES parsing
        Err(Error::InvalidOperation("SMILES parsing not yet implemented".to_string()))
    }
}

/// Implementation of MolecularEntity for GraphMolecule
impl MolecularEntity for GraphMolecule {
    fn add_to_builder(&self, mut builder: MoleculeBuilder, _attachment: Option<AtomIndex>) -> Result<MoleculeBuilder, Error> {
        // For now, simply create a new molecule - in the future we could support merging
        // with proper atom/bond mapping
        
        // This is just a placeholder implementation
        Ok(builder)
    }
}

/// A molecular fragment that can be used as a query or template
pub struct Fragment {
    molecule: Molecule,
}

impl Fragment {
    /// Create a new fragment from a molecule
    pub fn new(molecule: Molecule) -> Self {
        Self { molecule }
    }

    /// Get the underlying molecule
    pub fn molecule(&self) -> &Molecule {
        &self.molecule
    }
}

/// A query fragment for matching patterns in molecules
pub struct Query {
    fragment: Fragment,
}

impl Query {
    /// Create a new query from a fragment
    pub fn new(fragment: Fragment) -> Self {
        Self { fragment }
    }

    /// Get the underlying fragment
    pub fn fragment(&self) -> &Fragment {
        &self.fragment
    }
}

/// A template fragment for generating molecules
pub struct Template {
    fragment: Fragment,
}

impl Template {
    /// Create a new template from a fragment
    pub fn new(fragment: Fragment) -> Self {
        Self { fragment }
    }

    /// Get the underlying fragment
    pub fn fragment(&self) -> &Fragment {
        &self.fragment
    }
}

// NOTE: Specific fragment implementations (rings, functional groups, etc.)
// will be defined once we have SMILES parsing and will be based on SMILES strings
// rather than being hard-coded.