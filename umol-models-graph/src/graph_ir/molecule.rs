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
pub type AromaticSystemIndex = u32;
pub type MulticenterBondIndex = u32;

/// Resolved molecule in GraphIR. All atoms and bonds are fully validated.
#[derive(Debug, Clone)]
pub struct Molecule {
    graph: StableGraph<Atom, Bond, Undirected, u32>,
    aromatic_systems: Vec<AromaticSystem>,
    multicenter_bonds: Vec<MulticenterBond>,
}

impl Molecule {
    // Atoms
    pub fn atom_count(&self) -> usize {
        self.graph.node_count()
    }

    pub fn atom_indices(&self) -> impl Iterator<Item = AtomIndex> + '_ {
        self.graph.node_indices()
    }

    pub fn atoms(&self) -> impl Iterator<Item = &Atom> + '_ {
        self.graph.node_weights()
    }

    pub fn atom(&self, index: AtomIndex) -> Option<&Atom> {
        self.graph.node_weight(index)
    }

    // Bonds
    pub fn bond_count(&self) -> usize {
        self.graph.edge_count()
    }

    pub fn bond_indices(&self) -> impl Iterator<Item = BondIndex> + '_ {
        self.graph.edge_indices()
    }

    pub fn bonds(&self) -> impl Iterator<Item = &Bond> + '_ {
        self.graph.edge_weights()
    }

    pub fn bond(&self, index: BondIndex) -> Option<&Bond> {
        self.graph.edge_weight(index)
    }

    // Aromatic systems
    pub fn aromatic_system_count(&self) -> usize {
        self.aromatic_systems.len()
    }

    pub fn aromatic_systems_indices(&self) -> impl Iterator<Item = AromaticSystemIndex> + '_ {
        (0..self.aromatic_system_count()).map(|i| i as AromaticSystemIndex)
    }

    pub fn aromatic_systems(&self) -> impl Iterator<Item = &AromaticSystem> + '_ {
        self.aromatic_systems.iter()
    }

    pub fn aromatic_system(&self, index: AromaticSystemIndex) -> Option<&AromaticSystem> {
        self.aromatic_systems.get(index as usize)
    }

    // Multicenter bonds
    pub fn multicenter_bond_count(&self) -> usize {
        self.multicenter_bonds.len()
    }

    pub fn multicenter_bonds_indices(&self) -> impl Iterator<Item = MulticenterBondIndex> + '_ {
        (0..self.multicenter_bond_count()).map(|i| i as MulticenterBondIndex)
    }

    pub fn multicenter_bonds(&self) -> impl Iterator<Item = &MulticenterBond> + '_ {
        self.multicenter_bonds.iter()
    }

    pub fn multicenter_bond(&self, index: MulticenterBondIndex) -> Option<&MulticenterBond> {
        self.multicenter_bonds.get(index as usize)
    }

    // Atom-atom relationships
    pub fn atom_neighbor_indices(&self, index: AtomIndex) -> impl Iterator<Item = AtomIndex> + '_ {
        self.graph.neighbors(index)
    }

    pub fn atom_neighbors(&self, index: AtomIndex) -> impl Iterator<Item = &Atom> + '_ {
        self.graph
            .neighbors(index)
            .map(|n| self.graph.node_weight(n).unwrap())
    }

    // Atom-bond relationships
    pub fn atom_bond_indices(&self, index: AtomIndex) -> impl Iterator<Item = BondIndex> + '_ {
        self.graph.edges(index).map(|e| e.id())
    }

    pub fn atom_bonds(&self, index: AtomIndex) -> impl Iterator<Item = &Bond> + '_ {
        self.graph.edges(index).map(|e| e.weight())
    }

    pub fn connecting_bond_indices(
        &self,
        a: AtomIndex,
        b: AtomIndex,
    ) -> impl Iterator<Item = BondIndex> + '_ {
        self.graph.edges_connecting(a, b).map(|edge| edge.id())
    }

    pub fn connecting_bond(&self, a: AtomIndex, b: AtomIndex) -> Option<&Bond> {
        self.graph.edges_connecting(a, b).next().map(|e| e.weight())
    }

    pub fn bond_atom_indices(&self, index: BondIndex) -> Option<(AtomIndex, AtomIndex)> {
        self.graph.edge_endpoints(index)
    }

    pub fn bond_atoms(&self, index: BondIndex) -> Option<(&Atom, &Atom)> {
        self.graph.edge_endpoints(index).map(|(a, b)| {
            (
                self.graph.node_weight(a).unwrap(),
                self.graph.node_weight(b).unwrap(),
            )
        })
    }

    // Atom-aromatic system relationships
    pub fn atom_aromatic_systems_indices(
        &self,
        index: AtomIndex,
    ) -> impl Iterator<Item = AromaticSystemIndex> + '_ {
        self.aromatic_systems_indices()
            .filter(move |&i| self.aromatic_system(i).unwrap().contains_atom(index))
    }

    pub fn atom_aromatic_systems(
        &self,
        index: AtomIndex,
    ) -> impl Iterator<Item = AromaticSystem> + '_ {
        self.aromatic_systems()
            .filter(move |s| s.contains_atom(index))
            .map(|s| s.clone())
    }

    // TODO: Consider adding connecting_aromatic_system_indices and connecting_aromatic_system methods
    // that take a IntoIterator<Item = AtomIndex>

    // Atom-multicenter bond relationships
    pub fn atom_multicenter_bonds_indices(
        &self,
        index: AtomIndex,
    ) -> impl Iterator<Item = MulticenterBondIndex> + '_ {
        self.multicenter_bonds_indices()
            .filter(move |&i| self.multicenter_bond(i).unwrap().contains_atom(index))
    }

    pub fn atom_multicenter_bonds(
        &self,
        index: AtomIndex,
    ) -> impl Iterator<Item = MulticenterBond> + '_ {
        self.multicenter_bonds()
            .filter(move |b| b.contains_atom(index))
            .map(|b| b.clone())
    }

    // TODO: Consider adding connecting_multicenter_bond_indices and connecting_multicenter_bond methods
    // that take a IntoIterator<Item = AtomIndex>
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

    // Atoms
    pub fn atom_count(&self) -> usize {
        self.graph.node_count()
    }

    pub fn atom_indices(&self) -> impl Iterator<Item = AtomIndex> + '_ {
        self.graph.node_indices()
    }

    pub fn atoms(&self) -> impl Iterator<Item = &AtomBuilder> + '_ {
        self.graph.node_weights()
    }

    pub fn atom(&self, index: AtomIndex) -> Option<&AtomBuilder> {
        self.graph.node_weight(index)
    }

    pub fn atom_mut(&mut self, index: AtomIndex) -> Option<&mut AtomBuilder> {
        self.graph.node_weight_mut(index)
    }

    pub fn add_atom(&mut self, atom: AtomBuilder) -> AtomIndex {
        self.graph.add_node(atom)
    }

    pub fn remove_atom(&mut self, index: AtomIndex) -> Option<AtomBuilder> {
        self.graph.remove_node(index)
    }

    pub fn replace_atom(&mut self, index: AtomIndex, atom: AtomBuilder) -> Option<AtomBuilder> {
        self.graph
            .node_weight_mut(index)
            .map(|old| std::mem::replace(old, atom))
    }

    // Bonds
    pub fn bond_count(&self) -> usize {
        self.graph.edge_count()
    }

    pub fn bond_indices(&self) -> impl Iterator<Item = BondIndex> + '_ {
        self.graph.edge_indices()
    }

    pub fn bond(&self, index: BondIndex) -> Option<&Bond> {
        self.graph.edge_weight(index)
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
        self.graph
            .edge_weight_mut(index)
            .map(|old| std::mem::replace(old, bond))
    }

    // Aromatic systems
    pub fn aromatic_system_count(&self) -> usize {
        self.aromatic_systems.len()
    }

    pub fn aromatic_system_indices(&self) -> impl Iterator<Item = AromaticSystemIndex> + '_ {
        (0..self.aromatic_system_count()).map(|i| i as AromaticSystemIndex)
    }

    pub fn aromatic_systems(&self) -> impl Iterator<Item = &AromaticSystem> + '_ {
        self.aromatic_systems.iter()
    }

    pub fn aromatic_system(&self, index: AromaticSystemIndex) -> Option<&AromaticSystem> {
        self.aromatic_systems.get(index as usize)
    }

    pub fn aromatic_system_mut(&mut self, index: AromaticSystemIndex) -> Option<&mut AromaticSystem> {
        self.aromatic_systems.get_mut(index as usize)
    }

    pub fn add_aromatic_system(&mut self, system: AromaticSystem) {
        self.aromatic_systems.push(system);
    }

    pub fn remove_aromatic_system(&mut self, index: AromaticSystemIndex) -> Option<AromaticSystem> {
        let index = index as usize;
        if index >= self.aromatic_systems.len() {
            return None;
        }
        Some(self.aromatic_systems.remove(index))
    }

    pub fn replace_aromatic_system(
        &mut self,
        index: AromaticSystemIndex,
        system: AromaticSystem,
    ) -> Option<AromaticSystem> {
        let index = index as usize;
        if index >= self.aromatic_systems.len() {
            return None;
        }
        self.aromatic_systems
            .get_mut(index)
            .map(|s| std::mem::replace(s, system))
    }

    // Multicenter bonds
    pub fn multicenter_bond_count(&self) -> usize {
        self.multicenter_bonds.len()
    }

    pub fn multicenter_bond_indices(&self) -> impl Iterator<Item = MulticenterBondIndex> + '_ {
        (0..self.multicenter_bond_count()).map(|i| i as MulticenterBondIndex)
    }

    pub fn multicenter_bonds(&self) -> impl Iterator<Item = &MulticenterBond> + '_ {
        self.multicenter_bonds.iter()
    }

    pub fn multicenter_bond(&self, index: MulticenterBondIndex) -> Option<&MulticenterBond> {
        self.multicenter_bonds.get(index as usize)
    }

    pub fn multicenter_bond_mut(&mut self, index: MulticenterBondIndex) -> Option<&mut MulticenterBond> {
        self.multicenter_bonds.get_mut(index as usize)
    }

    pub fn add_multicenter_bond(&mut self, bond: MulticenterBond) {
        self.multicenter_bonds.push(bond);
    }

    pub fn remove_multicenter_bond(
        &mut self,
        index: MulticenterBondIndex,
    ) -> Option<MulticenterBond> {
        let index = index as usize;
        if index >= self.multicenter_bonds.len() {
            return None;
        }
        Some(self.multicenter_bonds.remove(index))
    }

    pub fn replace_multicenter_bond(
        &mut self,
        index: MulticenterBondIndex,
        bond: MulticenterBond,
    ) -> Option<MulticenterBond> {
        let index = index as usize;
        if index >= self.multicenter_bonds.len() {
            return None;
        }
        self.multicenter_bonds
            .get_mut(index)
            .map(|b| std::mem::replace(b, bond))
    }

    // Atom-atom relationships
    pub fn atom_neighbor_indices(&self, index: AtomIndex) -> impl Iterator<Item = AtomIndex> + '_ {
        self.graph.neighbors(index)
    }

    pub fn atom_neighbors(&self, index: AtomIndex) -> impl Iterator<Item = &AtomBuilder> + '_ {
        self.graph
            .neighbors(index)
            .map(|n| self.graph.node_weight(n).unwrap())
    }

    // Atom-bond relationships
    pub fn atom_bond_indices(&self, index: AtomIndex) -> impl Iterator<Item = BondIndex> + '_ {
        self.graph.edges(index).map(|e| e.id())
    }

    pub fn atom_bonds(&self, index: AtomIndex) -> impl Iterator<Item = &Bond> + '_ {
        self.graph.edges(index).map(|e| e.weight())
    }

    pub fn connecting_bond_indices(
        &self,
        a: AtomIndex,
        b: AtomIndex,
    ) -> impl Iterator<Item = BondIndex> + '_ {
        self.graph.edges_connecting(a, b).map(|e| e.id())
    }

    pub fn connecting_bond(&self, a: AtomIndex, b: AtomIndex) -> Option<&Bond> {
        self.graph.edges_connecting(a, b).next().map(|e| e.weight())
    }

    pub fn bond_atom_indices(&self, index: BondIndex) -> Option<(AtomIndex, AtomIndex)> {
        self.graph.edge_endpoints(index)
    }

    pub fn bond_atoms(&self, index: BondIndex) -> Option<(&AtomBuilder, &AtomBuilder)> {
        self.graph.edge_endpoints(index).map(|(a, b)| {
            (
                self.graph.node_weight(a).unwrap(),
                self.graph.node_weight(b).unwrap(),
            )
        })
    }

    // Atom-aromatic system relationships
    pub fn atom_aromatic_systems_indices(
        &self,
        index: AtomIndex,
    ) -> impl Iterator<Item = AromaticSystemIndex> + '_ {
        self.aromatic_system_indices()
            .filter(move |&i| self.aromatic_system(i).unwrap().contains_atom(index))
    }

    pub fn atom_aromatic_systems(
        &self,
        index: AtomIndex,
    ) -> impl Iterator<Item = AromaticSystem> + '_ {
        self.aromatic_systems()
            .filter(move |s| s.contains_atom(index))
            .map(|s| s.clone())
    }

    // TODO: Consider adding connecting_aromatic_system_indices and connecting_aromatic_system methods
    // that take a IntoIterator<Item = AtomIndex>

    // Atom-multicenter bond relationships
    pub fn atom_multicenter_bonds_indices(
        &self,
        index: AtomIndex,
    ) -> impl Iterator<Item = MulticenterBondIndex> + '_ {
        self.multicenter_bond_indices()
            .filter(move |&i| self.multicenter_bond(i).unwrap().contains_atom(index))
    }

    pub fn atom_multicenter_bonds(
        &self,
        index: AtomIndex,
    ) -> impl Iterator<Item = MulticenterBond> + '_ {
        self.multicenter_bonds()
            .filter(move |b| b.contains_atom(index))
            .map(|b| b.clone())
    }

    // TODO: Consider adding connecting_multicenter_bond_indices and connecting_multicenter_bond methods
    // that take a IntoIterator<Item = AtomIndex>

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

impl Default for MoleculeBuilder {
    fn default() -> Self {
        Self::new()
    }
}
