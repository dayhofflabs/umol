//! GraphIR molecule representation using typed atoms and bonds

use std::collections::HashMap;

use petgraph::graph::NodeIndex;
use petgraph::prelude::*;
use petgraph::stable_graph::StableGraph;
use petgraph::visit::EdgeRef;
use serde::{Deserialize, Serialize};
use umol_data::SpinState;

use super::aromaticity::AromaticSystem;
use super::atom::Atom;
use super::bond::Bond;
use super::dative::DativeBond;
use super::multicenter::MulticenterBond;
use super::noncovalent::NoncovalentBond;
use crate::algorithms::biconnected_components;
use crate::atom::AromaticValence;
use crate::dsl::ast::ToAst;
use crate::dsl::config::MoleculeDslConfig;
use crate::dsl::bond::BondAst;
use crate::dsl::molecule::{
    AromaticSystem as AromaticSystemAst, AtomRef, Atoms,
    LocalizedBond as LocalizedBondAst, DativeBond as DativeBondAst, MoleculeAst,
    MulticenterBond as MulticenterBondAst, NoncovalentBond as NoncovalentBondAst,
};

pub type AtomIndex = NodeIndex<u32>;
pub type BondIndex = EdgeIndex<u32>;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct DativeBondIndex(pub u32);

impl DativeBondIndex {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct AromaticSystemIndex(pub u32);

impl AromaticSystemIndex {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct MulticenterBondIndex(pub u32);

impl MulticenterBondIndex {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct NoncovalentBondIndex(pub u32);

impl NoncovalentBondIndex {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

/// Resolved molecule in GraphIR. All atoms and bonds are fully validated.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Molecule {
    graph: StableGraph<Atom, Bond, Undirected, u32>,
    dative_bonds: Vec<DativeBond>,
    aromatic_systems: Vec<AromaticSystem>,
    multicenter_bonds: Vec<MulticenterBond>,
    noncovalent_bonds: Vec<NoncovalentBond>,
    charge: i8,
    spin: SpinState,
}

impl Molecule {
    pub(crate) fn from_parts(
        graph: StableGraph<Atom, Bond, Undirected, u32>,
        dative_bonds: Vec<DativeBond>,
        aromatic_systems: Vec<AromaticSystem>,
        multicenter_bonds: Vec<MulticenterBond>,
        noncovalent_bonds: Vec<NoncovalentBond>,
        charge: i8,
        spin: SpinState,
    ) -> Self {
        Self {
            graph,
            dative_bonds,
            aromatic_systems,
            multicenter_bonds,
            noncovalent_bonds,
            charge,
            spin,
        }
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

    // Atom properties
    pub fn atom_valence(&self, index: AtomIndex) -> u8 {
        self.atom(index).map(|a| a.valence()).unwrap_or(0)
    }

    pub fn atom_aromatic_valence(&self, index: AtomIndex) -> u8 {
        self.atom(index)
            .map(|a| a.aromatic_valence())
            .map(|v| match v {
                AromaticValence::Valence(n) => n,
                AromaticValence::NotAromatic => 0,
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
    pub fn charge(&self) -> i8 {
        self.charge
    }

    pub fn spin(&self) -> SpinState {
        self.spin
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

impl ToAst<MoleculeAst> for Molecule {
    fn to_ast(&self, cfg: &MoleculeDslConfig) -> MoleculeAst {
        let mut atom_vec = Vec::new();

        for idx in self.atom_indices() {
            let atom_ast = self.atom(idx).unwrap().to_ast(&cfg.atom);
            atom_vec.push(atom_ast);
        }

        let label = |idx: AtomIndex| -> AtomRef {
            AtomRef::Index(idx.index())
        };

        let bonds: Vec<LocalizedBondAst> = self
            .bond_indices()
            .map(|bi| {
                let (a, b) = self.bond_atom_indices(bi).unwrap();
                LocalizedBondAst {
                    id: None,
                    a: label(a),
                    b: label(b),
                    bond: self.bond(bi).unwrap().to_ast(&cfg.bond),
                }
            })
            .collect();

        let dative_bonds: Vec<DativeBondAst> = self
            .dative_bonds()
            .map(|db| DativeBondAst {
                id: None,
                donor: label(db.donor()),
                acceptor: label(db.acceptor()),
                bond: BondAst::from_order(db.order()),
            })
            .collect();

        let aromatic_systems: Vec<AromaticSystemAst> = self
            .aromatic_systems()
            .map(|sys| AromaticSystemAst {
                id: None,
                atoms: sys.atoms().map(&label).collect(),
            })
            .collect();

        let multicenter_bonds: Vec<MulticenterBondAst> = self
            .multicenter_bonds()
            .map(|mc| MulticenterBondAst {
                id: None,
                atoms: mc.all_atoms().into_iter().map(&label).collect(),
            })
            .collect();

        let noncovalent_bonds: Vec<NoncovalentBondAst> = self
            .noncovalent_bonds()
            .map(|nc| NoncovalentBondAst {
                id: None,
                a: label(nc.a()),
                b: label(nc.b()),
                bond: BondAst::from_order(1),
            })
            .collect();

        MoleculeAst {
            atoms: Atoms::indexed(atom_vec),
            bonds,
            dative_bonds,
            aromatic_systems,
            multicenter_bonds,
            noncovalent_bonds,
            charge: Some(self.charge() as i64),
            spin: Some(self.spin()),
        }
    }
}

#[cfg(test)]
mod tests;
