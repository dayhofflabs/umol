//! Molecular graph model.

use crate::atom::AtomStandard;
use crate::bond::Bond;
use crate::conformer::Conformer;
use crate::sgroup::SGroup;
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::stable_graph::StableGraph;
use petgraph::Undirected;
use std::collections::HashMap;
use umol::error::DataError;
use umol::Result;

/// Type aliases for the node and edge indices
pub type AtomIndex = NodeIndex<usize>;
pub type BondIndex = EdgeIndex<usize>;

/// Graph-based molecule representation with MOL file semantics
#[derive(Debug, Clone)]
pub struct Molecule {
    pub graph: StableGraph<AtomStandard, Bond, Undirected, usize>,
    pub conformers: Vec<Conformer>,
    pub properties: HashMap<String, String>,
    pub sgroups: Vec<SGroup>,
}

impl Molecule {
    /// Create empty molecule
    pub fn new() -> Self {
        Self {
            graph: StableGraph::<AtomStandard, Bond, Undirected, usize>::default(),
            conformers: Vec::new(),
            properties: HashMap::new(),
            sgroups: Vec::new(),
        }
    }

    /// Parse molecule from MOL string (strict - only standard atoms)
    pub fn from_mol_str(input: &str) -> Result<Self> {
        Self::from_mol_bytes(input.as_bytes())
    }

    /// Parse molecule from MOL bytes (strict - only standard atoms)
    pub fn from_mol_bytes(input: &[u8]) -> Result<Self> {
        // This is a simplified implementation - in practice would need full MOL parsing
        // For now, just return an empty molecule
        Ok(Self::new())
    }

    /// Get number of atoms
    pub fn atom_count(&self) -> usize {
        self.graph.node_count()
    }

    /// Get number of bonds
    pub fn bond_count(&self) -> usize {
        self.graph.edge_count()
    }

    /// Get molecule-level property by key
    pub fn property(&self, key: &str) -> Option<&String> {
        self.properties.get(key)
    }

    /// Set molecule-level property by key
    pub fn set_property(&mut self, key: String, value: String) {
        self.properties.insert(key, value);
    }

    /// Get molecule-level properties as hashmap
    pub fn properties(&self) -> &HashMap<String, String> {
        &self.properties
    }

    /// Get mutable reference to molecule-level properties map
    pub fn properties_mut(&mut self) -> &mut HashMap<String, String> {
        &mut self.properties
    }

    /// Add atom to the molecule and update index mappings
    ///
    /// - `atom`: Atom to add (Molecule takes ownership)
    ///
    /// Return index of added atom.
    pub fn add_atom(&mut self, atom: AtomStandard) -> usize {
        self.graph.add_node(atom).index()
    }

    /// Add bond between two atoms specified by external/MOL indices
    ///
    /// - `idx1`, `idx2`: Atom indices
    /// - `bond`: Bond to add (Molecule takes ownership)
    ///
    /// Return index of added bond.
    pub fn add_bond(&mut self, idx1: usize, idx2: usize, bond: Bond) -> usize {
        self.graph
            .add_edge(AtomIndex::new(idx1), AtomIndex::new(idx2), bond)
            .index()
    }

    /// Get immutable reference to atom by index
    ///
    /// - `idx`: Atom index
    ///
    /// Return immutable reference to atom.
    pub fn atom(&self, idx: usize) -> Option<&AtomStandard> {
        self.graph.node_weight(AtomIndex::new(idx))
    }

    /// Get mutable reference to atom by index
    ///
    /// - `idx`: Atom index
    ///
    /// Return mutable reference to atom.
    pub fn atom_mut(&mut self, idx: usize) -> Option<&mut AtomStandard> {
        self.graph.node_weight_mut(AtomIndex::new(idx))
    }

    /// Get immutable reference to bond by index
    ///
    /// - `idx`: Bond index
    ///
    /// Return immutable reference to bond.
    pub fn bond(&self, idx: usize) -> Option<&Bond> {
        self.graph.edge_weight(BondIndex::new(idx))
    }

    /// Get mutable reference to bond by index
    ///
    /// - `idx`: Bond index
    ///
    /// Return mutable reference to bond.
    pub fn bond_mut(&mut self, idx: usize) -> Option<&mut Bond> {
        self.graph.edge_weight_mut(BondIndex::new(idx))
    }

    /// Get iterator over neighbor atom indices for atom index
    pub fn neighbors(&self, idx: usize) -> impl Iterator<Item = usize> + '_ {
        self.graph.neighbors(AtomIndex::new(idx)).map(|i| i.index())
    }

    /// Get immutable slice of all conformers
    pub fn conformers(&self) -> &[Conformer] {
        &self.conformers
    }

    /// Get mutable reference to vector of conformers
    pub fn conformers_mut(&mut self) -> &mut Vec<Conformer> {
        &mut self.conformers
    }

    /// Add conformer to the molecule
    ///
    /// - `conformer`: Conformer to add
    ///
    /// Return error if the number of positions in the conformer does not match
    /// the number of atoms in the molecule.
    pub fn add_conformer(&mut self, conformer: Conformer) -> Result<()> {
        let num_atoms = self.atom_count();
        let num_positions = conformer.positions.len();

        if num_atoms != num_positions {
            return Err(DataError::InvalidConformer(format!(
                "Expected {} positions, found {}",
                num_atoms, num_positions
            ))
            .into());
        }

        self.conformers.push(conformer);
        Ok(())
    }
}
