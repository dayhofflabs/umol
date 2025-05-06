//! Molecular graph model.

use crate::atom::Atom;
use crate::bond::Bond;
use crate::conformer::Conformer;
use crate::sgroup::SGroup;
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::stable_graph::StableGraph;
use petgraph::Undirected;
use std::collections::HashMap;
use umol::error::DataError;
use umol::{Error, Result};

/// Type alias for the node index type used in the molecular graph.
pub type AtomIndex = NodeIndex<usize>;

/// Type alias for the edge index type used in the molecular graph.
pub type BondIndex = EdgeIndex<usize>;

/// Molecule type represented by a graph of atoms and bonds.
#[derive(Debug, Clone)]
pub struct Molecule {
    /// The underlying graph storing atoms (nodes) and bonds (edges).
    pub graph: StableGraph<Atom, Bond, Undirected, usize>,
    /// A collection of conformers (3D coordinate sets) for the molecule.
    pub conformers: Vec<Conformer>,
    /// Molecule-level properties stored as key-value pairs.
    pub properties: HashMap<String, String>,
    /// Map from external/MOL 1-based indices to internal graph atom indices.
    /// Populated during parsing.
    pub external_indices: HashMap<usize, AtomIndex>,
    /// Map from internal graph atom indices back to external/MOL 1-based indices.
    /// Populated during parsing, used for writing.
    pub internal_indices: HashMap<AtomIndex, usize>,
    /// Collection of SGroups defined in the file.
    pub sgroups: Vec<SGroup>,
}

impl Molecule {
    /// Creates a new, empty Molecule.
    pub fn new() -> Self {
        Self {
            graph: StableGraph::<Atom, Bond, Undirected, usize>::default(),
            conformers: Vec::new(),
            properties: HashMap::new(),
            external_indices: HashMap::new(),
            internal_indices: HashMap::new(),
            sgroups: Vec::new(),
        }
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
    pub fn get_prop(&self, key: &str) -> Option<&String> {
        self.properties.get(key)
    }

    /// Set molecule-level property
    pub fn set_prop(&mut self, key: String, value: String) {
        self.properties.insert(key, value);
    }

    /// Get molecule-level properties as a map
    pub fn properties(&self) -> &HashMap<String, String> {
        &self.properties
    }

    /// Get mutable reference to molecule-level properties map
    pub fn properties_mut(&mut self) -> &mut HashMap<String, String> {
        &mut self.properties
    }

    /// Add atom to the molecule and update index mappings
    ///
    /// - `idx`: The external index (e.g., 1-based from MOL file).
    /// - `atom`: The Atom object to add.
    ///
    /// Returns the internal graph index (`AtomIndex`) of the added atom.
    pub fn add_atom(&mut self, idx: usize, atom: Atom) -> AtomIndex {
        let graph_index = self.graph.add_node(atom);
        self.external_indices.insert(idx, graph_index);
        self.internal_indices.insert(graph_index, idx);
        graph_index
    }

    /// Add bond between two atoms specified by external/MOL indices
    ///
    /// - `idx1`, `idx2`: External indices of the atoms to connect.
    /// - `bond`: The Bond object to add.
    ///
    /// Returns the internal graph index (`BondIndex`) of the added bond, or an error
    /// if either atom index is not found.
    pub fn add_bond(&mut self, idx1: usize, idx2: usize, bond: Bond) -> Result<BondIndex> {
        let graph_idx1 = *self
            .external_indices
            .get(&idx1)
            .ok_or_else::<Error, _>(|| DataError::MissingAtomIndex(idx1).into())?;
        let graph_idx2 = *self
            .external_indices
            .get(&idx2)
            .ok_or_else::<Error, _>(|| DataError::MissingAtomIndex(idx2).into())?;

        Ok(self.graph.add_edge(graph_idx1, graph_idx2, bond))
    }

    /// Get immutable reference to atom by internal graph index
    pub fn atom(&self, idx: AtomIndex) -> Option<&Atom> {
        self.graph.node_weight(idx)
    }

    /// Get mutable reference to atom by internal graph index
    pub fn atom_mut(&mut self, idx: AtomIndex) -> Option<&mut Atom> {
        self.graph.node_weight_mut(idx)
    }

    /// Get immutable reference to bond by internal graph index
    pub fn bond(&self, idx: BondIndex) -> Option<&Bond> {
        self.graph.edge_weight(idx)
    }

    /// Get mutable reference to bond by internal graph index
    pub fn bond_mut(&mut self, idx: BondIndex) -> Option<&mut Bond> {
        self.graph.edge_weight_mut(idx)
    }

    /// Get iterator over graph indices of neighbor atoms for a given atom index
    pub fn neighbors(&self, idx: AtomIndex) -> impl Iterator<Item = AtomIndex> + '_ {
        self.graph.neighbors(idx)
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
    /// Check if the conformer has the correct number of atomic positions.
    ///
    /// Returns an error if the number of positions in the conformer does not match
    /// the number of atoms in the molecule.
    pub fn add_conformer(&mut self, conformer: Conformer) -> Result<()> {
        let num_atoms = self.atom_count();
        let num_positions = conformer.positions.len();

        if num_atoms != num_positions {
            return Err(DataError::InvalidConformerDefinition(format!(
                "Expected {} positions, found {}",
                num_atoms, num_positions
            ))
            .into());
        }

        self.conformers.push(conformer);
        Ok(())
    }
}
