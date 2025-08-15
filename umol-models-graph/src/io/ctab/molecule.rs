//! Molecule type for CTab format.

use crate::io::ctab::atom::{Atom, AtomStandard};
use crate::io::ctab::bond::{Bond, BondStandard};
use crate::io::ctab::sgroup::SGroup;

use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::stable_graph::StableGraph;
use petgraph::Undirected;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};

/// Type aliases for the node and edge indices
pub type AtomIndex = NodeIndex<usize>;
pub type BondIndex = EdgeIndex<usize>;

/// Graph-based molecule representation with full MOL file semantics (including queries)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Molecule {
    pub graph: StableGraph<Atom, Bond, Undirected, usize>,
    pub sgroups: BTreeMap<usize, SGroup>,
    pub properties: HashMap<String, String>,
}

/// Graph-based molecule representation for standard (non-query) molecules only
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MoleculeStandard {
    pub graph: StableGraph<AtomStandard, BondStandard, Undirected, usize>,
    pub sgroups: BTreeMap<usize, SGroup>,
    pub properties: HashMap<String, String>,
}



impl Molecule {
    /// Create empty molecule
    pub fn new() -> Self {
        Self {
            graph: StableGraph::<Atom, Bond, Undirected, usize>::default(),
            sgroups: BTreeMap::new(),
            properties: HashMap::new(),
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
    pub fn add_atom(&mut self, atom: Atom) -> usize {
        self.graph.add_node(atom).index()
    }

    /// Add bond between two atoms specified by external/MOL indices
    pub fn add_bond(&mut self, idx1: usize, idx2: usize, bond: Bond) -> usize {
        self.graph
            .add_edge(AtomIndex::new(idx1), AtomIndex::new(idx2), bond)
            .index()
    }

    /// Get iterator over atom indices
    pub fn atom_indices(&self) -> impl Iterator<Item = usize> + '_ {
        self.graph.node_indices().map(|i| i.index())
    }

    /// Get iterator over atoms
    pub fn atoms(&self) -> impl Iterator<Item = &Atom> + '_ {
        self.graph.node_weights()
    }

    /// Get immutable reference to atom by index
    pub fn atom(&self, idx: usize) -> Option<&Atom> {
        self.graph.node_weight(AtomIndex::new(idx))
    }

    /// Get mutable reference to atom by index
    pub fn atom_mut(&mut self, idx: usize) -> Option<&mut Atom> {
        self.graph.node_weight_mut(AtomIndex::new(idx))
    }

    /// Get iterator over bond indices
    pub fn bond_indices(&self) -> impl Iterator<Item = usize> + '_ {
        self.graph.edge_indices().map(|i| i.index())
    }

    /// Get iterator over bonds
    pub fn bonds(&self) -> impl Iterator<Item = &Bond> + '_ {
        self.graph.edge_weights()
    }

    /// Get immutable reference to bond by index
    pub fn bond(&self, idx: usize) -> Option<&Bond> {
        self.graph.edge_weight(BondIndex::new(idx))
    }

    /// Get mutable reference to bond by index
    pub fn bond_mut(&mut self, idx: usize) -> Option<&mut Bond> {
        self.graph.edge_weight_mut(BondIndex::new(idx))
    }

    /// Get iterator over neighbor atom indices for atom index
    pub fn neighbors(&self, idx: usize) -> impl Iterator<Item = usize> + '_ {
        self.graph.neighbors(AtomIndex::new(idx)).map(|i| i.index())
    }

    /// Get iterator over sgroup indices
    pub fn sgroup_indices(&self) -> impl Iterator<Item = usize> + '_ {
        self.sgroups.keys().copied()
    }

    /// Get iterator over sgroups
    pub fn sgroups(&self) -> impl Iterator<Item = &SGroup> + '_ {
        self.sgroups.values()
    }

    /// Get immutable reference to sgroup by index
    pub fn sgroup(&self, idx: usize) -> Option<&SGroup> {
        self.sgroups.get(&idx)
    }

    /// Get mutable reference to sgroup by index
    pub fn sgroup_mut(&mut self, idx: usize) -> Option<&mut SGroup> {
        self.sgroups.get_mut(&idx)
    }

    /// Add sgroup to the molecule
    pub fn add_sgroup(&mut self, sgroup_index: usize, sgroup: SGroup) {
        self.sgroups.insert(sgroup_index, sgroup);
    }
}

impl MoleculeStandard {
    /// Create empty standard molecule
    pub fn new() -> Self {
        Self {
            graph: StableGraph::<AtomStandard, BondStandard, Undirected, usize>::default(),
            sgroups: BTreeMap::new(),
            properties: HashMap::new(),
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

    /// Add atom to the molecule and update index mappings
    pub fn add_atom(&mut self, atom: AtomStandard) -> usize {
        self.graph.add_node(atom).index()
    }

    /// Add bond between two atoms specified by external/MOL indices
    pub fn add_bond(&mut self, idx1: usize, idx2: usize, bond: BondStandard) -> usize {
        self.graph
            .add_edge(AtomIndex::new(idx1), AtomIndex::new(idx2), bond)
            .index()
    }

    /// Get immutable reference to atom by index
    pub fn atom(&self, idx: usize) -> Option<&AtomStandard> {
        self.graph.node_weight(AtomIndex::new(idx))
    }

    /// Get mutable reference to atom by index
    pub fn atom_mut(&mut self, idx: usize) -> Option<&mut AtomStandard> {
        self.graph.node_weight_mut(AtomIndex::new(idx))
    }

    /// Get immutable reference to bond by index
    pub fn bond(&self, idx: usize) -> Option<&BondStandard> {
        self.graph.edge_weight(BondIndex::new(idx))
    }

    /// Get mutable reference to bond by index
    pub fn bond_mut(&mut self, idx: usize) -> Option<&mut BondStandard> {
        self.graph.edge_weight_mut(BondIndex::new(idx))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_molecule_standard_serialize() {
        let graph =
            StableGraph::<AtomStandard, BondStandard, Undirected, usize>::with_capacity(0, 0);
        let sgroups = BTreeMap::new();
        let properties = HashMap::new();

        let molecule = MoleculeStandard {
            graph,
            sgroups,
            properties,
        };

        let yaml =
            serde_yaml::to_string(&molecule).expect("Failed to serialize MoleculeStandard to YAML");
        let deserialized: MoleculeStandard =
            serde_yaml::from_str(&yaml).expect("Failed to deserialize MoleculeStandard from YAML");
        assert_eq!(molecule.properties, deserialized.properties);
    }

    #[test]
    fn test_molecule_serialize() {
        let graph = StableGraph::<Atom, Bond, Undirected, usize>::with_capacity(0, 0);
        let sgroups = BTreeMap::new();
        let properties = HashMap::new();
        let molecule = Molecule {
            graph,
            sgroups,
            properties,
        };

        let yaml = serde_yaml::to_string(&molecule).expect("Failed to serialize Molecule to YAML");
        let deserialized: Molecule =
            serde_yaml::from_str(&yaml).expect("Failed to deserialize Molecule from YAML");
        assert_eq!(molecule.properties, deserialized.properties);
    }
}
