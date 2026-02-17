//! GraphIR molecule representation built on typed atoms and bonds.

use petgraph::graph::NodeIndex;
use petgraph::prelude::*;
use petgraph::stable_graph::StableGraph;
use petgraph::visit::{EdgeRef, NodeIndexable};

use super::aromatic::AromaticSystem;
use super::atom::{Atom, AtomBuilder};
use super::bond::Bond;
use super::config::ResolveConfig;
use super::error::ResolutionError;
use super::multicenter::MulticenterBond;

pub type AtomIndex = NodeIndex<u32>;
pub type BondIndex = EdgeIndex<u32>;

/// Resolved molecule in GraphIR. All atoms and bonds are fully validated.
#[derive(Debug, Clone)]
pub struct Molecule {
    graph: StableGraph<Atom, Bond, Undirected, u32>,
    aromatic_systems: Vec<AromaticSystem>,
    multicenter_bonds: Vec<MulticenterBond>,
}

impl Molecule {
    pub fn atom_count(&self) -> usize {
        self.graph.node_count()
    }

    pub fn atoms<'graph>(&'graph self) -> impl Iterator<Item = &'graph Atom> + 'graph {
        self.graph.node_weights()
    }

    pub fn atom(&self, index: AtomIndex) -> Option<&Atom> {
        self.graph.node_weight(index)
    }

    pub fn bond_count(&self) -> usize {
        self.graph.edge_count()
    }

    pub fn bonds<'graph>(&'graph self) -> impl Iterator<Item = &'graph Bond> + 'graph {
        self.graph.edge_weights()
    }

    pub fn bond(&self, index: BondIndex) -> Option<&Bond> {
        self.graph.edge_weight(index)
    }

    pub fn atom_indices(&self) -> impl Iterator<Item = AtomIndex> + '_ {
        self.graph.node_indices()
    }

    pub fn bond_indices(&self) -> impl Iterator<Item = BondIndex> + '_ {
        self.graph.edge_indices()
    }

    pub fn bond_atoms(&self, index: BondIndex) -> Option<(&Atom, &Atom)> {
        self.graph.edge_endpoints(index).map(|(a, b)| {
            (
                self.graph.node_weight(a).unwrap(),
                self.graph.node_weight(b).unwrap(),
            )
        })
    }

    pub fn bond_atom_indices(&self, index: BondIndex) -> Option<(AtomIndex, AtomIndex)> {
        self.graph.edge_endpoints(index)
    }

    pub fn bonds_between<'graph>(
        &'graph self,
        a: AtomIndex,
        b: AtomIndex,
    ) -> impl Iterator<Item = BondIndex> + 'graph {
        self.graph.edges_connecting(a, b).map(|edge| edge.id())
    }

    pub fn atom_bonds<'graph>(
        &'graph self,
        index: AtomIndex,
    ) -> impl Iterator<Item = &'graph Bond> + 'graph {
        self.graph.edges(index).map(|e| e.weight())
    }

    pub fn atom_bond_indices<'graph>(
        &'graph self,
        index: AtomIndex,
    ) -> impl Iterator<Item = BondIndex> + 'graph {
        self.graph.edges(index).map(|e| e.id())
    }

    pub fn atom_neighbors<'graph>(
        &'graph self,
        index: AtomIndex,
    ) -> impl Iterator<Item = &'graph Atom> + 'graph {
        self.graph
            .neighbors(index)
            .map(|n| self.graph.node_weight(n).unwrap())
    }

    pub fn atom_neighbor_indices<'graph>(
        &'graph self,
        index: AtomIndex,
    ) -> impl Iterator<Item = AtomIndex> + 'graph {
        self.graph.neighbors(index)
    }

    pub fn aromatic_systems(&self) -> &[AromaticSystem] {
        &self.aromatic_systems
    }

    pub fn multicenter_bonds(&self) -> &[MulticenterBond] {
        &self.multicenter_bonds
    }

    pub fn add_atom(&mut self, atom: Atom) -> AtomIndex {
        self.graph.add_node(atom)
    }

    pub fn remove_atom(&mut self, index: AtomIndex) -> Option<Atom> {
        self.graph.remove_node(index)
    }

    pub fn replace_atom(&mut self, index: AtomIndex, atom: Atom) -> Option<Atom> {
        if let Some(slot) = self.graph.node_weight_mut(index) {
            let old = std::mem::replace(slot, atom);
            Some(old)
        } else {
            None
        }
    }

    pub fn add_bond(&mut self, a: AtomIndex, b: AtomIndex, bond: Bond) -> Option<BondIndex> {
        if !self.graph.contains_node(a) || !self.graph.contains_node(b) {
            return None;
        }
        Some(self.graph.add_edge(a, b, bond))
    }

    pub fn remove_bond(&mut self, index: BondIndex) -> Option<Bond> {
        self.graph.remove_edge(index)
    }

    pub fn replace_bond(&mut self, index: BondIndex, bond: Bond) -> Option<Bond> {
        if let Some(slot) = self.graph.edge_weight_mut(index) {
            let old = std::mem::replace(slot, bond);
            Some(old)
        } else {
            None
        }
    }
}

/// Builder for constructing a `Molecule`. Carries `AtomBuilder` nodes during
/// resolution phases; `build()` finalizes each atom and produces a `Molecule`.
///
/// Used both by the resolution pipeline (from TableIR) and for manual
/// molecule construction.
#[derive(Debug, Clone)]
pub struct MoleculeBuilder {
    graph: StableGraph<AtomBuilder, Bond, Undirected, u32>,
    aromatic_systems: Vec<AromaticSystem>,
    multicenter_bonds: Vec<MulticenterBond>,
}

impl MoleculeBuilder {
    pub fn new() -> Self {
        Self {
            graph: StableGraph::default(),
            aromatic_systems: Vec::new(),
            multicenter_bonds: Vec::new(),
        }
    }

    pub fn with_capacity(atom_capacity: usize, bond_capacity: usize) -> Self {
        Self {
            graph: StableGraph::with_capacity(atom_capacity, bond_capacity),
            aromatic_systems: Vec::new(),
            multicenter_bonds: Vec::new(),
        }
    }

    pub fn atom_count(&self) -> usize {
        self.graph.node_count()
    }

    pub fn bond_count(&self) -> usize {
        self.graph.edge_count()
    }

    pub fn atom_builder(&self, index: AtomIndex) -> Option<&AtomBuilder> {
        self.graph.node_weight(index)
    }

    pub fn atom_builder_mut(&mut self, index: AtomIndex) -> Option<&mut AtomBuilder> {
        self.graph.node_weight_mut(index)
    }

    pub fn bond(&self, index: BondIndex) -> Option<&Bond> {
        self.graph.edge_weight(index)
    }

    pub fn atom_indices(&self) -> impl Iterator<Item = AtomIndex> + '_ {
        self.graph.node_indices()
    }

    pub fn bond_indices(&self) -> impl Iterator<Item = BondIndex> + '_ {
        self.graph.edge_indices()
    }

    pub fn bond_atom_indices(&self, index: BondIndex) -> Option<(AtomIndex, AtomIndex)> {
        self.graph.edge_endpoints(index)
    }

    pub fn atom_bond_indices<'graph>(
        &'graph self,
        index: AtomIndex,
    ) -> impl Iterator<Item = BondIndex> + 'graph {
        self.graph.edges(index).map(|e| e.id())
    }

    pub fn atom_neighbor_indices<'graph>(
        &'graph self,
        index: AtomIndex,
    ) -> impl Iterator<Item = AtomIndex> + 'graph {
        self.graph.neighbors(index)
    }

    pub fn add_atom(&mut self, atom: AtomBuilder) -> AtomIndex {
        self.graph.add_node(atom)
    }

    pub fn add_bond(
        &mut self,
        a: AtomIndex,
        b: AtomIndex,
        bond: Bond,
    ) -> Option<BondIndex> {
        if !self.graph.contains_node(a) || !self.graph.contains_node(b) {
            return None;
        }
        Some(self.graph.add_edge(a, b, bond))
    }

    pub fn remove_atom(&mut self, index: AtomIndex) -> Option<AtomBuilder> {
        self.graph.remove_node(index)
    }

    pub fn remove_bond(&mut self, index: BondIndex) -> Option<Bond> {
        self.graph.remove_edge(index)
    }

    /// Build the final `Molecule` by finalizing each `AtomBuilder` into an `Atom`.
    ///
    /// Requires all atom builders to have exactly one valence candidate
    /// remaining (i.e., resolution phases must have been run).
    pub fn build(self, _config: &ResolveConfig) -> Result<Molecule, ResolutionError> {
        let mut resolved_graph =
            StableGraph::with_capacity(self.graph.node_count(), self.graph.edge_count());

        let mut index_map = Vec::with_capacity(self.graph.node_bound());
        index_map.resize(self.graph.node_bound(), None);

        for old_idx in self.graph.node_indices() {
            let builder = self.graph.node_weight(old_idx).unwrap();
            let atom = builder.build()?;
            let new_idx = resolved_graph.add_node(atom);
            index_map[old_idx.index()] = Some(new_idx);
        }

        for old_edge in self.graph.edge_indices() {
            let (a, b) = self.graph.edge_endpoints(old_edge).unwrap();
            let bond = self.graph.edge_weight(old_edge).unwrap().clone();
            let new_a = index_map[a.index()].unwrap();
            let new_b = index_map[b.index()].unwrap();
            resolved_graph.add_edge(new_a, new_b, bond);
        }

        Ok(Molecule {
            graph: resolved_graph,
            aromatic_systems: self.aromatic_systems,
            multicenter_bonds: self.multicenter_bonds,
        })
    }
}
