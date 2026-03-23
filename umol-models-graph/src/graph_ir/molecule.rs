//! GraphIR molecule representation using typed atoms and bonds

use std::collections::{HashMap, HashSet};

use petgraph::graph::NodeIndex;
use petgraph::prelude::*;
use petgraph::stable_graph::StableGraph;
use petgraph::visit::EdgeRef;
use umol_data::SpinState;

use crate::algorithms::biconnected_components;
use crate::atom::AromaticValence;
use crate::graph_ir::aromaticity::AromaticSystem;
use crate::graph_ir::atom::Atom;
use crate::graph_ir::bond::Bond;
use crate::graph_ir::dative::DativeBond;
use crate::graph_ir::multicenter::MulticenterBond;
use crate::graph_ir::noncovalent::NoncovalentBond;

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

    // Atom properties
    pub fn atom_valence(&self, index: AtomIndex) -> u8 {
        self.atom(index).map(|a| a.valence()).unwrap_or(0)
    }

    pub fn atom_aromatic_valence(&self, index: AtomIndex) -> u8 {
        self.atom(index)
            .map(|a| a.aromatic_valence())
            .map(|v| match v {
                AromaticValence::Valence(n) => n,
                AromaticValence::None => 0,
            })
            .unwrap_or(0)
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

    pub fn biconnected_components(&self) -> Vec<Vec<AtomIndex>> {
        let mut atoms: Vec<AtomIndex> = self.atom_indices().collect();
        atoms.sort_unstable();
        if atoms.is_empty() {
            return Vec::new();
        }

        let atom_to_id: HashMap<AtomIndex, usize> = atoms
            .iter()
            .copied()
            .enumerate()
            .map(|(i, a)| (a, i))
            .collect();
        let adj = self.adjacency_list();
        let mut adj_int: Vec<Vec<usize>> = vec![Vec::new(); atoms.len()];
        for &atom in &atoms {
            let mut neighbors = adj
                .get(&atom)
                .map(|ns| {
                    ns.iter()
                        .filter_map(|&n| atom_to_id.get(&n).copied())
                        .collect::<Vec<_>>()
                })
                .unwrap_or_default();
            neighbors.sort_unstable();
            neighbors.dedup();
            let u = atom_to_id[&atom];
            adj_int[u] = neighbors;
        }

        biconnected_components(atoms.len(), &adj_int)
            .into_iter()
            .map(|component| component.into_iter().map(|i| atoms[i]).collect())
            .collect()
    }

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
    // Iterator over all bonded atom pairs
    pub fn bonded_atom_pairs(&self) -> impl Iterator<Item = (AtomIndex, AtomIndex)> + '_ {
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

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
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
    fn naphthalene_molecule() -> Molecule {
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
        let carbon = spec!("{Cv4}");
        for atom in builder.atom_indices().collect::<Vec<_>>() {
            builder
                .atom_mut(atom)
                .unwrap()
                .set_candidates(SmallVec::from_elem(carbon, 1));
        }
        builder
            .build(&ResolveConfig::default())
            .expect("test molecule should build")
    }

    #[rstest]
    #[case::naphthalene(naphthalene_molecule(), vec![10])]
    fn test_biconnected_components(#[case] molecule: Molecule, #[case] expected_sizes: Vec<usize>) {
        let mut actual_sizes: Vec<usize> = molecule
            .biconnected_components()
            .iter()
            .map(|c| c.len())
            .collect();
        actual_sizes.sort_unstable();
        assert_eq!(actual_sizes, expected_sizes);
    }

    #[test]
    fn test_atom_aromatic_valence_resolved_semantics() {
        let mut aromatic_builder = MoleculeBuilder::new();
        let aromatic_atom = aromatic_builder.add_atom(AtomBuilder::new(Element::C));
        aromatic_builder
            .atom_mut(aromatic_atom)
            .expect("atom should exist")
            .set_candidates(SmallVec::from_elem(spec!("{Cv2a1H}"), 1));
        let aromatic = aromatic_builder
            .build(&ResolveConfig::default())
            .expect("aromatic molecule should build");
        assert_eq!(aromatic.atom_aromatic_valence(aromatic_atom), 1);

        let mut non_aromatic_builder = MoleculeBuilder::new();
        let non_aromatic_atom = non_aromatic_builder.add_atom(AtomBuilder::new(Element::C));
        non_aromatic_builder
            .atom_mut(non_aromatic_atom)
            .expect("atom should exist")
            .set_candidates(SmallVec::from_elem(spec!("{Cv4}"), 1));
        let non_aromatic = non_aromatic_builder
            .build(&ResolveConfig::default())
            .expect("non-aromatic molecule should build");
        assert_eq!(non_aromatic.atom_aromatic_valence(non_aromatic_atom), 0);
        assert_eq!(non_aromatic.atom_aromatic_valence(AtomIndex::new(999)), 0);
    }
}
