//! GraphIR molecule builder.

use std::collections::{HashMap, HashSet};

use petgraph::prelude::*;
use petgraph::stable_graph::StableGraph;
use petgraph::visit::{EdgeRef, NodeIndexable};
use umol_data::SpinState;

use super::super::aromatic::AromaticSystem;
use super::super::atom::AtomBuilder;
use super::super::bond::BondBuilder;
use super::super::config::ResolveConfig;
use super::super::dative::DativeBond;
use super::super::error::ResolutionError;
use super::super::multicenter::MulticenterBond;
use super::super::noncovalent::NoncovalentBond;
use super::super::symmetry::compute_symmetry;
use super::utils::biconnected_components;
use super::{
    AromaticSystemIndex, AtomIndex, BondIndex, DativeBondIndex, Molecule, MulticenterBondIndex,
    NoncovalentBondIndex, TopologyExportError, TopologyGraph, TopologyProjection,
};
/// Builder for constructing a `Molecule`. Carries `AtomBuilder` nodes during
/// resolution phases; `build()` finalizes each atom and produces a `Molecule`.
///
/// Used both by the resolution pipeline (from TableIR) and for manual
/// molecule construction.
#[derive(Debug, Clone)]
pub struct MoleculeBuilder {
    graph: StableGraph<AtomBuilder, BondBuilder, Undirected, u32>,
    dative_bonds: Vec<DativeBond>,
    aromatic_systems: Vec<AromaticSystem>,
    multicenter_bonds: Vec<MulticenterBond>,
    noncovalent_bonds: Vec<NoncovalentBond>,
    charge: Option<i32>,
    spin: Option<SpinState>,
}

impl MoleculeBuilder {
    pub fn new() -> Self {
        Self {
            graph: StableGraph::default(),
            dative_bonds: Vec::new(),
            aromatic_systems: Vec::new(),
            multicenter_bonds: Vec::new(),
            noncovalent_bonds: Vec::new(),
            charge: None,
            spin: None,
        }
    }

    pub fn with_capacity(atom_capacity: usize, bond_capacity: usize) -> Self {
        Self {
            graph: StableGraph::with_capacity(atom_capacity, bond_capacity),
            dative_bonds: Vec::new(),
            aromatic_systems: Vec::new(),
            multicenter_bonds: Vec::new(),
            noncovalent_bonds: Vec::new(),
            charge: None,
            spin: None,
        }
    }

    pub fn topology_graph(&self, projection: TopologyProjection) -> TopologyGraph {
        TopologyGraph::from_builder(self, projection)
    }

    pub fn topology_canonical_bfs(&self, projection: TopologyProjection) -> Vec<NodeIndex> {
        let canonical_atoms = compute_symmetry(self).canonical_order();
        let mut atom_rank = HashMap::<AtomIndex, usize>::new();
        for (rank, atom) in canonical_atoms.into_iter().enumerate() {
            atom_rank.insert(atom, rank);
        }
        let graph = self.topology_graph(projection);
        graph.canonical_bfs_with_rank(|node_ref| match node_ref {
            super::TopologyNodeRef::Atom(ai) => atom_rank
                .get(&ai)
                .copied()
                .unwrap_or(usize::MAX / 4 + ai.index()),
            super::TopologyNodeRef::Bond(i) => usize::MAX / 2 + i.index(),
            super::TopologyNodeRef::DativeBond(i) => usize::MAX / 2 + 1_000_000 + i.index(),
            super::TopologyNodeRef::NoncovalentBond(i) => usize::MAX / 2 + 2_000_000 + i.index(),
            super::TopologyNodeRef::MulticenterBond(i) => usize::MAX / 2 + 3_000_000 + i.index(),
        })
    }

    pub fn topology_graph6_canonical(
        &self,
        projection: TopologyProjection,
    ) -> Result<String, TopologyExportError> {
        let canonical_atoms = compute_symmetry(self).canonical_order();
        let mut atom_rank = HashMap::<AtomIndex, usize>::new();
        for (rank, atom) in canonical_atoms.into_iter().enumerate() {
            atom_rank.insert(atom, rank);
        }
        let graph = self.topology_graph(projection);
        let (g6, _order) = graph.to_graph6_canonical_with_rank(|node_ref| match node_ref {
            super::TopologyNodeRef::Atom(ai) => atom_rank
                .get(&ai)
                .copied()
                .unwrap_or(usize::MAX / 4 + ai.index()),
            super::TopologyNodeRef::Bond(i) => usize::MAX / 2 + i.index(),
            super::TopologyNodeRef::DativeBond(i) => usize::MAX / 2 + 1_000_000 + i.index(),
            super::TopologyNodeRef::NoncovalentBond(i) => usize::MAX / 2 + 2_000_000 + i.index(),
            super::TopologyNodeRef::MulticenterBond(i) => usize::MAX / 2 + 3_000_000 + i.index(),
        })?;
        Ok(g6)
    }

    pub fn biconnected_components(&self) -> Vec<Vec<AtomIndex>> {
        biconnected_components(self.atom_indices(), self.adjacency_list())
    }

    pub(crate) fn topology_nodes(&self) -> impl Iterator<Item = AtomIndex> + '_ {
        self.atom_indices()
    }

    pub(crate) fn topology_edges(
        &self,
        projection: TopologyProjection,
    ) -> std::vec::IntoIter<super::topology::TopologyEdge> {
        let mut edges = Vec::new();

        for i in self.bond_indices() {
            if let Some((a, b)) = self.bond_atom_indices(i) {
                edges.push(super::topology::TopologyEdge::Edge {
                    node_ref: super::TopologyNodeRef::Bond(i),
                    a,
                    b,
                });
            }
        }

        if projection.dative == super::DativeProjection::Undirected {
            for i in self.dative_bond_indices() {
                if let Some(b) = self.dative_bond(i) {
                    edges.push(super::topology::TopologyEdge::Edge {
                        node_ref: super::TopologyNodeRef::DativeBond(i),
                        a: b.donor(),
                        b: b.acceptor(),
                    });
                }
            }
        }

        if projection.noncovalent == super::NoncovalentProjection::Undirected {
            for i in self.noncovalent_bond_indices() {
                if let Some(b) = self.noncovalent_bond(i) {
                    edges.push(super::topology::TopologyEdge::Edge {
                        node_ref: super::TopologyNodeRef::NoncovalentBond(i),
                        a: b.a(),
                        b: b.b(),
                    });
                }
            }
        }

        match projection.multicenter {
            super::MulticenterProjection::Skip => {}
            super::MulticenterProjection::CliqueExpansion => {
                for i in self.multicenter_bond_indices() {
                    if let Some(mc) = self.multicenter_bond(i) {
                        let mut seen = HashSet::new();
                        let atoms: Vec<AtomIndex> = mc
                            .all_atoms()
                            .into_iter()
                            .filter(|a| seen.insert(*a))
                            .collect();
                        for x in 0..atoms.len() {
                            for y in (x + 1)..atoms.len() {
                                edges.push(super::topology::TopologyEdge::Edge {
                                    node_ref: super::TopologyNodeRef::MulticenterBond(i),
                                    a: atoms[x],
                                    b: atoms[y],
                                });
                            }
                        }
                    }
                }
            }
            super::MulticenterProjection::IncidenceNode => {
                for i in self.multicenter_bond_indices() {
                    if let Some(mc) = self.multicenter_bond(i) {
                        let mut seen = HashSet::new();
                        let atoms: Vec<AtomIndex> = mc
                            .all_atoms()
                            .into_iter()
                            .filter(|a| seen.insert(*a))
                            .collect();
                        if !atoms.is_empty() {
                            edges.push(super::topology::TopologyEdge::Hyperedge {
                                node_ref: super::TopologyNodeRef::MulticenterBond(i),
                                atoms,
                            });
                        }
                    }
                }
            }
        }

        edges.into_iter()
    }

    // Atoms
    pub fn atom_count(&self) -> usize {
        self.graph.node_count()
    }

    pub fn atom_indices(&self) -> impl Iterator<Item = AtomIndex> + '_ {
        self.graph.node_indices()
    }

    /// Atoms that have at least one candidate with a non-None aromatic valence.
    pub fn aromatic_candidate_atoms(&self) -> impl Iterator<Item = AtomIndex> + '_ {
        self.atom_indices().filter(|&atom| {
            self.atom(atom)
                .map(|a| {
                    a.candidates()
                        .iter()
                        .any(|c| c.aromatic_valence().is_aromatic())
                })
                .unwrap_or(false)
        })
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

    pub fn bond(&self, index: BondIndex) -> Option<&BondBuilder> {
        self.graph.edge_weight(index)
    }

    pub fn bond_mut(&mut self, index: BondIndex) -> Option<&mut BondBuilder> {
        self.graph.edge_weight_mut(index)
    }

    pub fn add_bond(&mut self, a: AtomIndex, b: AtomIndex, bond: BondBuilder) -> Option<BondIndex> {
        if !self.graph.contains_node(a) || !self.graph.contains_node(b) {
            return None;
        }
        Some(self.graph.add_edge(a, b, bond))
    }

    pub fn add_bond_unchecked(
        &mut self,
        a: AtomIndex,
        b: AtomIndex,
        bond: BondBuilder,
    ) -> BondIndex {
        debug_assert!(
            self.graph.contains_node(a),
            "atom index {:?} not in builder",
            a
        );
        debug_assert!(
            self.graph.contains_node(b),
            "atom index {:?} not in builder",
            b
        );
        self.graph.add_edge(a, b, bond)
    }

    pub fn remove_bond(&mut self, index: BondIndex) -> Option<BondBuilder> {
        self.graph.remove_edge(index)
    }

    pub fn replace_bond(&mut self, index: BondIndex, bond: BondBuilder) -> Option<BondBuilder> {
        self.graph
            .edge_weight_mut(index)
            .map(|old| std::mem::replace(old, bond))
    }

    // Dative bonds
    pub fn dative_bond_count(&self) -> usize {
        self.dative_bonds.len()
    }

    pub fn dative_bond_indices(&self) -> impl Iterator<Item = DativeBondIndex> + '_ {
        (0..self.dative_bond_count()).map(|i| DativeBondIndex(i as u32))
    }

    pub fn dative_bonds(&self) -> impl Iterator<Item = &DativeBond> + '_ {
        self.dative_bonds.iter()
    }

    pub fn dative_bond(&self, index: DativeBondIndex) -> Option<&DativeBond> {
        self.dative_bonds.get(index.index())
    }

    pub fn dative_bond_mut(&mut self, index: DativeBondIndex) -> Option<&mut DativeBond> {
        self.dative_bonds.get_mut(index.index())
    }

    pub fn add_dative_bond(&mut self, bond: DativeBond) {
        self.dative_bonds.push(bond);
    }

    pub fn remove_dative_bond(&mut self, index: DativeBondIndex) -> Option<DativeBond> {
        let i = index.index();
        if i >= self.dative_bonds.len() {
            return None;
        }
        Some(self.dative_bonds.remove(i))
    }

    pub fn replace_dative_bond(
        &mut self,
        index: DativeBondIndex,
        bond: DativeBond,
    ) -> Option<DativeBond> {
        self.dative_bonds
            .get_mut(index.index())
            .map(|b| std::mem::replace(b, bond))
    }

    // Aromatic systems
    pub fn aromatic_system_count(&self) -> usize {
        self.aromatic_systems.len()
    }

    pub fn aromatic_system_indices(&self) -> impl Iterator<Item = AromaticSystemIndex> + '_ {
        (0..self.aromatic_system_count()).map(|i| AromaticSystemIndex(i as u32))
    }

    pub fn aromatic_systems(&self) -> impl Iterator<Item = &AromaticSystem> + '_ {
        self.aromatic_systems.iter()
    }

    pub fn aromatic_system(&self, index: AromaticSystemIndex) -> Option<&AromaticSystem> {
        self.aromatic_systems.get(index.index())
    }

    pub fn aromatic_system_mut(
        &mut self,
        index: AromaticSystemIndex,
    ) -> Option<&mut AromaticSystem> {
        self.aromatic_systems.get_mut(index.index())
    }

    pub fn add_aromatic_system(&mut self, system: AromaticSystem) {
        self.aromatic_systems.push(system);
    }

    pub fn remove_aromatic_system(&mut self, index: AromaticSystemIndex) -> Option<AromaticSystem> {
        let i = index.index();
        if i >= self.aromatic_systems.len() {
            return None;
        }
        Some(self.aromatic_systems.remove(i))
    }

    pub fn replace_aromatic_system(
        &mut self,
        index: AromaticSystemIndex,
        system: AromaticSystem,
    ) -> Option<AromaticSystem> {
        self.aromatic_systems
            .get_mut(index.index())
            .map(|s| std::mem::replace(s, system))
    }

    // Multicenter bonds
    pub fn multicenter_bond_count(&self) -> usize {
        self.multicenter_bonds.len()
    }

    pub fn multicenter_bond_indices(&self) -> impl Iterator<Item = MulticenterBondIndex> + '_ {
        (0..self.multicenter_bond_count()).map(|i| MulticenterBondIndex(i as u32))
    }

    pub fn multicenter_bonds(&self) -> impl Iterator<Item = &MulticenterBond> + '_ {
        self.multicenter_bonds.iter()
    }

    pub fn multicenter_bond(&self, index: MulticenterBondIndex) -> Option<&MulticenterBond> {
        self.multicenter_bonds.get(index.index())
    }

    pub fn multicenter_bond_mut(
        &mut self,
        index: MulticenterBondIndex,
    ) -> Option<&mut MulticenterBond> {
        self.multicenter_bonds.get_mut(index.index())
    }

    pub fn add_multicenter_bond(&mut self, bond: MulticenterBond) {
        self.multicenter_bonds.push(bond);
    }

    pub fn remove_multicenter_bond(
        &mut self,
        index: MulticenterBondIndex,
    ) -> Option<MulticenterBond> {
        let i = index.index();
        if i >= self.multicenter_bonds.len() {
            return None;
        }
        Some(self.multicenter_bonds.remove(i))
    }

    pub fn replace_multicenter_bond(
        &mut self,
        index: MulticenterBondIndex,
        bond: MulticenterBond,
    ) -> Option<MulticenterBond> {
        self.multicenter_bonds
            .get_mut(index.index())
            .map(|b| std::mem::replace(b, bond))
    }

    // Non-covalent bonds
    pub fn noncovalent_bond_count(&self) -> usize {
        self.noncovalent_bonds.len()
    }

    pub fn noncovalent_bond_indices(&self) -> impl Iterator<Item = NoncovalentBondIndex> + '_ {
        (0..self.noncovalent_bond_count()).map(|i| NoncovalentBondIndex(i as u32))
    }

    pub fn noncovalent_bonds(&self) -> impl Iterator<Item = &NoncovalentBond> + '_ {
        self.noncovalent_bonds.iter()
    }

    pub fn noncovalent_bond(&self, index: NoncovalentBondIndex) -> Option<&NoncovalentBond> {
        self.noncovalent_bonds.get(index.index())
    }

    pub fn noncovalent_bond_mut(
        &mut self,
        index: NoncovalentBondIndex,
    ) -> Option<&mut NoncovalentBond> {
        self.noncovalent_bonds.get_mut(index.index())
    }

    pub fn add_noncovalent_bond(&mut self, bond: NoncovalentBond) {
        self.noncovalent_bonds.push(bond);
    }

    pub fn remove_noncovalent_bond(
        &mut self,
        index: NoncovalentBondIndex,
    ) -> Option<NoncovalentBond> {
        let i = index.index();
        if i >= self.noncovalent_bonds.len() {
            return None;
        }
        Some(self.noncovalent_bonds.remove(i))
    }

    pub fn replace_noncovalent_bond(
        &mut self,
        index: NoncovalentBondIndex,
        bond: NoncovalentBond,
    ) -> Option<NoncovalentBond> {
        self.noncovalent_bonds
            .get_mut(index.index())
            .map(|b| std::mem::replace(b, bond))
    }

    // Molecular charge and spin
    pub fn charge(&self) -> Option<i32> {
        self.charge
    }

    pub fn spin(&self) -> Option<SpinState> {
        self.spin
    }

    pub fn set_charge(&mut self, charge: i32) {
        self.charge = Some(charge);
    }

    pub fn set_spin(&mut self, spin: SpinState) {
        self.spin = Some(spin);
    }

    pub fn update_charge(&mut self, f: impl FnOnce(Option<i32>) -> Option<i32>) {
        self.charge = f(self.charge);
    }

    pub fn update_spin(&mut self, f: impl FnOnce(Option<SpinState>) -> Option<SpinState>) {
        self.spin = f(self.spin);
    }

    // Atom-atom relationships
    pub fn adjacency_list(&self) -> HashMap<AtomIndex, Vec<AtomIndex>> {
        let mut adj = HashMap::with_capacity(self.graph.node_count());
        for atom in self.graph.node_indices() {
            adj.insert(atom, Vec::new());
        }
        for bond in self.graph.edge_indices() {
            let (a, b) = self.graph.edge_endpoints(bond).unwrap();
            adj.get_mut(&a).unwrap().push(b);
            adj.get_mut(&b).unwrap().push(a);
        }
        adj
    }

    pub fn atom_neighbor_indices(&self, index: AtomIndex) -> impl Iterator<Item = AtomIndex> + '_ {
        self.graph.neighbors(index)
    }

    pub fn atom_neighbors(&self, index: AtomIndex) -> impl Iterator<Item = &AtomBuilder> + '_ {
        self.graph
            .neighbors(index)
            .map(|n| self.graph.node_weight(n).unwrap())
    }

    // TODO: Add dative and noncovalent neighbors (+indices)
    // TODO: Add aromatic system and multicenter system partners (+indices)

    // Atom-bond relationships
    pub fn atom_bond_count(&self, index: AtomIndex) -> usize {
        self.graph.edges(index).count()
    }

    pub fn atom_bond_indices(&self, index: AtomIndex) -> impl Iterator<Item = BondIndex> + '_ {
        self.graph.edges(index).map(|e| e.id())
    }

    pub fn atom_bonds(&self, index: AtomIndex) -> impl Iterator<Item = &BondBuilder> + '_ {
        self.graph.edges(index).map(|e| e.weight())
    }

    pub fn atom_bond_order_sum(&self, index: AtomIndex) -> u8 {
        self.graph.edges(index).map(|e| e.weight().order()).sum()
    }

    pub fn atom_bond_order_sum_with_aromatic(&self, index: AtomIndex, aromatic_weight: f32) -> f32 {
        self.graph
            .edges(index)
            .map(|e| {
                let b = e.weight();
                if b.aromatic_hint() == Some(true) {
                    aromatic_weight
                } else {
                    b.order() as f32
                }
            })
            .sum()
    }

    pub fn connecting_bond_index(&self, a: AtomIndex, b: AtomIndex) -> Option<BondIndex> {
        self.graph.edges_connecting(a, b).next().map(|e| e.id())
    }

    pub fn connecting_bond(&self, a: AtomIndex, b: AtomIndex) -> Option<&BondBuilder> {
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

    // Atom-dative bond relationships
    pub fn atom_has_dative_bonds(&self, index: AtomIndex) -> bool {
        self.dative_bonds.iter().any(|b| b.contains_atom(index))
    }

    pub fn atom_dative_bond_counts(&self, index: AtomIndex) -> (usize, usize) {
        assert!(
            self.graph.contains_node(index),
            "atom index {:?} not in builder",
            index
        );
        let mut donated = 0;
        let mut accepted = 0;
        for db in &self.dative_bonds {
            if db.donor() == index {
                donated += 1;
            } else if db.acceptor() == index {
                accepted += 1;
            }
        }
        (donated, accepted)
    }

    pub fn atom_dative_bond_indices(
        &self,
        index: AtomIndex,
    ) -> impl Iterator<Item = DativeBondIndex> + '_ {
        self.dative_bond_indices()
            .filter(move |&i| self.dative_bond(i).unwrap().contains_atom(index))
    }

    pub fn atom_dative_bonds(&self, index: AtomIndex) -> impl Iterator<Item = &DativeBond> + '_ {
        self.dative_bonds().filter(move |b| b.contains_atom(index))
    }

    pub fn atom_dative_bond_order_sums(&self, index: AtomIndex) -> (u8, u8) {
        debug_assert!(
            self.graph.contains_node(index),
            "atom index {:?} not in builder",
            index
        );

        let mut donated = 0;
        let mut accepted = 0;
        for db in &self.dative_bonds {
            if db.donor() == index {
                donated += db.order();
            } else if db.acceptor() == index {
                accepted += db.order();
            }
        }
        (donated, accepted)
    }

    // Atom-aromatic system relationships
    pub fn atom_has_aromatic_systems(&self, index: AtomIndex) -> bool {
        debug_assert!(
            self.graph.contains_node(index),
            "atom index {:?} not in builder",
            index
        );
        self.atom_aromatic_systems_indices(index).next().is_some()
    }

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

    // Atom-multicenter bond relationships
    pub fn atom_has_multicenter_bonds(&self, index: AtomIndex) -> bool {
        debug_assert!(
            self.graph.contains_node(index),
            "atom index {:?} not in builder",
            index
        );
        self.atom_multicenter_bonds_indices(index).next().is_some()
    }

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

    // Atom-noncovalent bond relationships
    pub fn atom_has_noncovalent_bonds(&self, index: AtomIndex) -> bool {
        self.noncovalent_bonds
            .iter()
            .any(|b| b.contains_atom(index))
    }

    pub fn atom_noncovalent_bond_indices(
        &self,
        index: AtomIndex,
    ) -> impl Iterator<Item = NoncovalentBondIndex> + '_ {
        self.noncovalent_bond_indices()
            .filter(move |&i| self.noncovalent_bond(i).unwrap().contains_atom(index))
    }

    pub fn atom_noncovalent_bonds(
        &self,
        index: AtomIndex,
    ) -> impl Iterator<Item = &NoncovalentBond> + '_ {
        self.noncovalent_bonds()
            .filter(move |b| b.contains_atom(index))
    }

    /// Build the final `Molecule` by finalizing each `AtomBuilder` into an `Atom`.
    ///
    /// Requires all atom builders to have exactly one valence candidate
    /// remaining (i.e., resolution phases must have been run).
    pub fn build(self, _config: &ResolveConfig) -> Result<Molecule, ResolutionError> {
        let mut graph =
            StableGraph::with_capacity(self.graph.node_count(), self.graph.edge_count());

        let mut index_map = Vec::with_capacity(self.graph.node_bound());
        index_map.resize(self.graph.node_bound(), None);

        for old_idx in self.graph.node_indices() {
            let builder = self.graph.node_weight(old_idx).unwrap();
            let atom = builder.build()?;
            let new_idx = graph.add_node(atom);
            index_map[old_idx.index()] = Some(new_idx);
        }

        for old_edge in self.graph.edge_indices() {
            let (a, b) = self.graph.edge_endpoints(old_edge).unwrap();
            let bond_builder = self.graph.edge_weight(old_edge).unwrap();
            let new_a = index_map[a.index()].unwrap();
            let new_b = index_map[b.index()].unwrap();
            graph.add_edge(new_a, new_b, bond_builder.build());
        }

        // TODO: add aromatic system charge
        let charge: i32 = graph.node_weights().map(|a| a.charge() as i32).sum();

        if let Some(explicit) = self.charge {
            if explicit != charge {
                return Err(ResolutionError::MolecularChargeMismatch {
                    explicit,
                    atom_sum: charge,
                });
            }
        }

        // TODO: add aromatic system spin
        let atom_spins: Vec<SpinState> = graph.node_weights().map(|a| a.spin()).collect();

        let spin = match self.spin {
            Some(explicit) => {
                if !explicit.is_compatible(&atom_spins) {
                    let atom_unpaired_sum: u16 = atom_spins
                        .iter()
                        .map(|s| s.unpaired_electrons() as u16)
                        .sum();
                    return Err(ResolutionError::MolecularSpinIncompatible {
                        explicit_unpaired: explicit.unpaired_electrons(),
                        explicit_multiplicity: explicit.multiplicity().multiplicity(),
                        atom_unpaired_sum,
                    });
                }
                explicit
            }
            None => SpinState::high_spin(&atom_spins).ok_or_else(|| {
                ResolutionError::ValenceViolation(
                    graph.node_weights().next().unwrap().element(),
                    "molecular spin exceeds maximum representable".to_string(),
                )
            })?,
        };

        Ok(Molecule {
            graph,
            aromatic_systems: self.aromatic_systems,
            multicenter_bonds: self.multicenter_bonds,
            dative_bonds: self.dative_bonds,
            noncovalent_bonds: self.noncovalent_bonds,
            charge,
            spin,
        })
    }
}

impl Default for MoleculeBuilder {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use smallvec::SmallVec;
    use umol_data::Element;

    use super::*;
    use crate::graph_ir::atom::AtomBuilder;
    use crate::graph_ir::bond::BondBuilder;
    use crate::graph_ir::config::ResolveConfig;
    use crate::graph_ir::molecule::Molecule;
    use crate::spec;

    #[fixture]
    fn empty_builder() -> MoleculeBuilder {
        MoleculeBuilder::new()
    }

    #[fixture]
    fn single_atom_builder() -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomBuilder::new(Element::C));
        builder
    }

    #[fixture]
    fn ring_builder(#[default(6)] n: usize) -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..n)
            .map(|_| builder.add_atom(AtomBuilder::new(Element::C)))
            .collect();
        for i in 0..n {
            builder.add_bond_unchecked(atoms[i], atoms[(i + 1) % n], BondBuilder::new(1, None));
        }
        builder
    }

    #[fixture]
    fn chain_builder(#[default(5)] n: usize) -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..n)
            .map(|_| builder.add_atom(AtomBuilder::new(Element::C)))
            .collect();
        for i in 0..n - 1 {
            builder.add_bond_unchecked(atoms[i], atoms[i + 1], BondBuilder::new(1, None));
        }
        builder
    }

    #[fixture]
    fn naphthalene_builder() -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..10)
            .map(|_| builder.add_atom(AtomBuilder::new(Element::C)))
            .collect();
        let ring1_edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];
        for (a, b) in ring1_edges {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondBuilder::new(1, None));
        }
        let ring2_edges = [(3, 6), (6, 7), (7, 8), (8, 9), (9, 4)];
        for (a, b) in ring2_edges {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondBuilder::new(1, None));
        }
        builder
    }

    #[fixture]
    fn cubane_builder() -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..8)
            .map(|_| builder.add_atom(AtomBuilder::new(Element::C)))
            .collect();
        let edges = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ];
        for (a, b) in edges {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondBuilder::new(1, None));
        }
        builder
    }

    #[fixture]
    fn spiro_builder() -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..5)
            .map(|_| builder.add_atom(AtomBuilder::new(Element::C)))
            .collect();
        let edges = [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)];
        for (a, b) in edges {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondBuilder::new(1, None));
        }
        builder
    }

    #[fixture]
    fn bridged_builder() -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..6)
            .map(|_| builder.add_atom(AtomBuilder::new(Element::C)))
            .collect();
        let ring1_edges = [(0, 2), (2, 1), (1, 3), (3, 0)];
        for (a, b) in ring1_edges {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondBuilder::new(1, None));
        }
        let ring2_edges = [(0, 4), (4, 1), (1, 5), (5, 0)];
        for (a, b) in ring2_edges {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondBuilder::new(1, None));
        }
        builder
    }

    #[fixture]
    fn naphthalene_molecule(mut naphthalene_builder: MoleculeBuilder) -> Molecule {
        let carbon = spec!("{Cv4}");
        for atom in naphthalene_builder.atom_indices().collect::<Vec<_>>() {
            naphthalene_builder
                .atom_mut(atom)
                .unwrap()
                .set_candidates(SmallVec::from_elem(carbon, 1));
        }
        naphthalene_builder
            .build(&ResolveConfig::default())
            .expect("test molecule should build")
    }

    #[rstest]
    #[case::empty(empty_builder(), vec![])]
    #[case::single_atom(single_atom_builder(), vec![])]
    #[case::chain(chain_builder(5), vec![])]
    #[case::single_ring(ring_builder(6), vec![6])]
    #[case::naphthalene(naphthalene_builder(), vec![10])]
    #[case::spiro(spiro_builder(), vec![3, 3])]
    #[case::cubane(cubane_builder(), vec![8])]
    fn test_biconnected_components(
        #[case] builder: MoleculeBuilder,
        #[case] expected_sizes: Vec<usize>,
    ) {
        let mut actual_sizes: Vec<usize> = builder
            .biconnected_components()
            .iter()
            .map(|c| c.len())
            .collect();
        actual_sizes.sort_unstable();
        assert_eq!(actual_sizes, expected_sizes);
    }
}
