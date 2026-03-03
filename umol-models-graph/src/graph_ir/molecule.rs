//! GraphIR molecule representation using typed atoms and bonds

use std::collections::HashSet;

use petgraph::graph::NodeIndex;
use petgraph::prelude::*;
use petgraph::stable_graph::StableGraph;
use petgraph::visit::EdgeRef;
use umol_data::SpinState;

use super::aromatic::AromaticSystem;
use super::atom::Atom;
use super::bond::Bond;
use super::dative::DativeBond;
use super::multicenter::MulticenterBond;
use super::noncovalent::NoncovalentBond;

pub mod builder;
pub mod topology;
pub use builder::*;
pub use topology::*;

pub type AtomIndex = NodeIndex<u32>;
pub type BondIndex = EdgeIndex<u32>;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct DativeBondIndex(pub u32);

impl DativeBondIndex {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct AromaticSystemIndex(pub u32);

impl AromaticSystemIndex {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MulticenterBondIndex(pub u32);

impl MulticenterBondIndex {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct NoncovalentBondIndex(pub u32);

impl NoncovalentBondIndex {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

/// Resolved molecule in GraphIR. All atoms and bonds are fully validated.
#[derive(Debug, Clone)]
pub struct Molecule {
    graph: StableGraph<Atom, Bond, Undirected, u32>,
    dative_bonds: Vec<DativeBond>,
    aromatic_systems: Vec<AromaticSystem>,
    multicenter_bonds: Vec<MulticenterBond>,
    noncovalent_bonds: Vec<NoncovalentBond>,
    charge: i32,
    spin: SpinState,
}

impl Molecule {
    // Topology
    pub fn topology_graph(&self, projection: TopologyProjection) -> TopologyGraph {
        TopologyGraph::from_molecule(self, projection)
    }

    pub fn topology_canonical_bfs(&self, projection: TopologyProjection) -> Vec<NodeIndex> {
        self.topology_graph(projection).canonical_bfs()
    }

    pub(crate) fn topology_nodes(&self) -> impl Iterator<Item = AtomIndex> + '_ {
        self.atom_indices()
    }

    pub(crate) fn topology_edges(
        &self,
        projection: TopologyProjection,
    ) -> impl Iterator<Item = TopologyEdge> + '_ {
        let mut edges = Vec::new();

        for i in self.bond_indices() {
            if let Some((a, b)) = self.bond_atom_indices(i) {
                edges.push(topology::TopologyEdge::Edge {
                    node_ref: TopologyNodeRef::Bond(i),
                    a,
                    b,
                });
            }
        }

        if projection.dative == DativeProjection::Undirected {
            for i in self.dative_bond_indices() {
                if let Some(b) = self.dative_bond(i) {
                    edges.push(topology::TopologyEdge::Edge {
                        node_ref: TopologyNodeRef::DativeBond(i),
                        a: b.donor(),
                        b: b.acceptor(),
                    });
                }
            }
        }

        if projection.noncovalent == NoncovalentProjection::Undirected {
            for i in self.noncovalent_bond_indices() {
                if let Some(b) = self.noncovalent_bond(i) {
                    edges.push(topology::TopologyEdge::Edge {
                        node_ref: TopologyNodeRef::NoncovalentBond(i),
                        a: b.a(),
                        b: b.b(),
                    });
                }
            }
        }

        match projection.multicenter {
            MulticenterProjection::Skip => {}
            MulticenterProjection::CliqueExpansion => {
                for i in self.multicenter_bonds_indices() {
                    if let Some(mc) = self.multicenter_bond(i) {
                        let mut seen = HashSet::new();
                        let atoms: Vec<AtomIndex> = mc
                            .all_atoms()
                            .into_iter()
                            .filter(|a| seen.insert(*a))
                            .collect();
                        for x in 0..atoms.len() {
                            for y in (x + 1)..atoms.len() {
                                edges.push(topology::TopologyEdge::Edge {
                                    node_ref: TopologyNodeRef::MulticenterBond(i),
                                    a: atoms[x],
                                    b: atoms[y],
                                });
                            }
                        }
                    }
                }
            }
            MulticenterProjection::IncidenceNode => {
                for i in self.multicenter_bonds_indices() {
                    if let Some(mc) = self.multicenter_bond(i) {
                        let mut seen = HashSet::new();
                        let atoms: Vec<AtomIndex> = mc
                            .all_atoms()
                            .into_iter()
                            .filter(|a| seen.insert(*a))
                            .collect();
                        if !atoms.is_empty() {
                            edges.push(topology::TopologyEdge::Hyperedge {
                                node_ref: TopologyNodeRef::MulticenterBond(i),
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

    // Aromatic systems
    pub fn aromatic_system_count(&self) -> usize {
        self.aromatic_systems.len()
    }

    pub fn aromatic_systems_indices(&self) -> impl Iterator<Item = AromaticSystemIndex> + '_ {
        (0..self.aromatic_system_count()).map(|i| AromaticSystemIndex(i as u32))
    }

    pub fn aromatic_systems(&self) -> impl Iterator<Item = &AromaticSystem> + '_ {
        self.aromatic_systems.iter()
    }

    pub fn aromatic_system(&self, index: AromaticSystemIndex) -> Option<&AromaticSystem> {
        self.aromatic_systems.get(index.index())
    }

    // Multicenter bonds
    pub fn multicenter_bond_count(&self) -> usize {
        self.multicenter_bonds.len()
    }

    pub fn multicenter_bonds_indices(&self) -> impl Iterator<Item = MulticenterBondIndex> + '_ {
        (0..self.multicenter_bond_count()).map(|i| MulticenterBondIndex(i as u32))
    }

    pub fn multicenter_bonds(&self) -> impl Iterator<Item = &MulticenterBond> + '_ {
        self.multicenter_bonds.iter()
    }

    pub fn multicenter_bond(&self, index: MulticenterBondIndex) -> Option<&MulticenterBond> {
        self.multicenter_bonds.get(index.index())
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

    // Charge
    pub fn charge(&self) -> i32 {
        self.charge
    }

    pub fn spin(&self) -> SpinState {
        self.spin
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

    // TODO: Add dative and noncovalent neighbors (+indices)
    // TODO: Add aromatic system and multicenter system partners (+indices)

    // Atom-bond relationships
    pub fn bonded_atoms(&self) -> impl Iterator<Item = (AtomIndex, AtomIndex)> + '_ {
        self.graph
            .edge_indices()
            .map(|e| self.graph.edge_endpoints(e).unwrap())
    }

    pub fn atom_bond_count(&self, index: AtomIndex) -> usize {
        self.graph.edges(index).count()
    }

    pub fn atom_bond_indices(&self, index: AtomIndex) -> impl Iterator<Item = BondIndex> + '_ {
        self.graph.edges(index).map(|e| e.id())
    }

    pub fn atom_bonds(&self, index: AtomIndex) -> impl Iterator<Item = &Bond> + '_ {
        self.graph.edges(index).map(|e| e.weight())
    }

    pub fn atom_bond_order_sum(&self, index: AtomIndex) -> u8 {
        self.graph.edges(index).map(|e| e.weight().order()).sum()
    }

    pub fn connecting_bond_index(&self, a: AtomIndex, b: AtomIndex) -> Option<BondIndex> {
        self.graph
            .edges_connecting(a, b)
            .next()
            .map(|edge| edge.id())
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

    // Atom-dative bond relationships
    pub fn dative_bonded_atoms(&self) -> impl Iterator<Item = (AtomIndex, AtomIndex)> + '_ {
        self.dative_bonds.iter().map(|b| (b.donor(), b.acceptor()))
    }

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
    pub fn aromatic_bonded_atoms(&self) -> impl Iterator<Item = Vec<AtomIndex>> + '_ {
        self.aromatic_systems.iter().map(|s| s.atoms().collect())
    }

    pub fn atom_has_aromatic_systems(&self, index: AtomIndex) -> bool {
        self.aromatic_systems.iter().any(|s| s.contains_atom(index))
    }

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

    // Atom-multicenter bond relationships
    pub fn multicenter_bonded_atoms(&self) -> impl Iterator<Item = Vec<AtomIndex>> + '_ {
        self.multicenter_bonds.iter().map(|b| b.all_atoms())
    }

    pub fn atom_has_multicenter_bonds(&self, index: AtomIndex) -> bool {
        self.multicenter_bonds
            .iter()
            .any(|b| b.contains_atom(index))
    }

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

    // Atom-noncovalent bond relationships
    pub fn noncovalent_bonded_atoms(&self) -> impl Iterator<Item = (AtomIndex, AtomIndex)> + '_ {
        self.noncovalent_bonds.iter().map(|b| (b.a(), b.b()))
    }

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
}
