//! GraphIR molecule representation built on typed atoms and bonds.

use std::collections::{HashMap, HashSet};
use std::fmt;

use indexmap::IndexMap;
use petgraph::graph::NodeIndex;
use petgraph::prelude::*;
use petgraph::stable_graph::StableGraph;
use petgraph::visit::EdgeRef;
use umol::error::DataError;
use umol::Result;

use super::atom::{Atom, AtomBuilder};
use super::atom_matcher::{AtomMatcher, STRICT_ATOM_MATCHER};
use super::atom_validator::{AtomValidator, STRICT_ATOM_VALIDATOR};
use super::bond::{Bond, BondBuilder};
use super::bond_matcher::{BondMatcher, STRICT_BOND_MATCHER};

pub type AtomIndex = NodeIndex<usize>;
pub type BondIndex = EdgeIndex<usize>;

#[derive(Debug, Clone)]
pub struct Molecule {
    data: StableGraph<Atom, Bond, Undirected, usize>,
}

impl Molecule {
    pub fn atom_count(&self) -> usize {
        self.data.node_count()
    }

    pub fn atoms<'graph>(&'graph self) -> impl Iterator<Item = &'graph Atom> + 'graph {
        self.data.node_weights()
    }

    pub fn atom(&self, index: AtomIndex) -> Option<&Atom> {
        self.data.node_weight(index)
    }

    pub fn bond_count(&self) -> usize {
        self.data.edge_count()
    }

    pub fn bonds<'graph>(&'graph self) -> impl Iterator<Item = &'graph Bond> + 'graph {
        self.data.edge_weights()
    }

    pub fn bond(&self, index: BondIndex) -> Option<&Bond> {
        self.data.edge_weight(index)
    }

    pub fn atom_indices(&self) -> impl Iterator<Item = AtomIndex> + '_ {
        self.data.node_indices()
    }

    pub fn bond_indices(&self) -> impl Iterator<Item = BondIndex> + '_ {
        self.data.edge_indices()
    }

    pub fn bond_atoms(&self, index: BondIndex) -> Option<(&Atom, &Atom)> {
        self.data.edge_endpoints(index).map(|(a, b)| {
            (
                self.data.node_weight(a).unwrap(),
                self.data.node_weight(b).unwrap(),
            )
        })
    }

    pub fn bond_atom_indices(&self, index: BondIndex) -> Option<(AtomIndex, AtomIndex)> {
        self.data.edge_endpoints(index)
    }

    pub fn bonds_between<'graph>(
        &'graph self,
        a: AtomIndex,
        b: AtomIndex,
    ) -> impl Iterator<Item = BondIndex> + 'graph {
        self.data.edges_connecting(a, b).map(|edge| edge.id())
    }

    pub fn atom_bonds<'graph>(
        &'graph self,
        index: AtomIndex,
    ) -> impl Iterator<Item = &'graph Bond> + 'graph {
        self.data.edges(index).map(|e| e.weight())
    }

    pub fn atom_bond_indices<'graph>(
        &'graph self,
        index: AtomIndex,
    ) -> impl Iterator<Item = BondIndex> + 'graph {
        self.data.edges(index).map(|e| e.id())
    }

    pub fn atom_neighbors<'graph>(
        &'graph self,
        index: AtomIndex,
    ) -> impl Iterator<Item = &'graph Atom> + 'graph {
        self.data
            .neighbors(index)
            .map(|n| self.data.node_weight(n).unwrap())
    }

    pub fn atom_neighbor_indices<'graph>(
        &'graph self,
        index: AtomIndex,
    ) -> impl Iterator<Item = AtomIndex> + 'graph {
        self.data.neighbors(index)
    }

    pub fn add_atom(&mut self, atom: Atom) -> AtomIndex {
        self.data.add_node(atom)
    }

    pub fn remove_atom(&mut self, index: AtomIndex) -> Option<Atom> {
        self.data.remove_node(index)
    }

    pub fn replace_atom(&mut self, index: AtomIndex, atom: Atom) -> Option<Atom> {
        if let Some(slot) = self.data.node_weight_mut(index) {
            let old = std::mem::replace(slot, atom);
            Some(old)
        } else {
            None
        }
    }

    pub fn add_bond(&mut self, a: AtomIndex, b: AtomIndex, bond: Bond) -> Option<BondIndex> {
        if !self.data.contains_node(a) || !self.data.contains_node(b) {
            return None;
        }
        Some(self.data.add_edge(a, b, bond))
    }

    pub fn remove_bond(&mut self, index: BondIndex) -> Option<Bond> {
        self.data.remove_edge(index)
    }

    pub fn replace_bond(&mut self, index: BondIndex, bond: Bond) -> Option<Bond> {
        if let Some(slot) = self.data.edge_weight_mut(index) {
            let old = std::mem::replace(slot, bond);
            Some(old)
        } else {
            None
        }
    }
}

impl fmt::Display for Molecule {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(
            f,
            "Molecule with {} atoms and {} bonds:",
            self.atom_count(),
            self.bond_count()
        )?;
        for (i, atom) in self.atoms().enumerate() {
            writeln!(f, "  Atom {}: {:?}", i, atom)?;
        }
        for (i, bond) in self.bonds().enumerate() {
            if let Some((a, b)) = self.bond_atoms(BondIndex::new(i)) {
                writeln!(
                    f,
                    "  Bond {}: {} between atoms {:?} and {:?}",
                    i, bond, a, b
                )?;
            }
        }
        Ok(())
    }
}

pub struct MoleculeBuilder {
    atom_builders: HashMap<usize, AtomBuilder>,
    bond_builders: HashMap<(usize, usize), BondBuilder>,
}

impl Default for MoleculeBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl MoleculeBuilder {
    pub fn new() -> Self {
        Self {
            atom_builders: HashMap::new(),
            bond_builders: HashMap::new(),
        }
    }

    pub fn with_capacity(atom_capacity: usize, bond_capacity: usize) -> Self {
        Self {
            atom_builders: HashMap::with_capacity(atom_capacity),
            bond_builders: HashMap::with_capacity(bond_capacity),
        }
    }

    pub fn create_atom<A: Into<AtomBuilder>>(&mut self, atom: A) -> (usize, &mut AtomBuilder) {
        let builder = atom.into();
        let idx = self.atom_builders.len();
        self.atom_builders.insert(idx, builder);
        (idx, self.atom_builders.get_mut(&idx).unwrap())
    }

    pub fn create_atoms<A: Into<AtomBuilder>>(
        &mut self,
        atoms: impl IntoIterator<Item = A>,
    ) -> impl Iterator<Item = usize> {
        let builders_iter = atoms.into_iter().map(|atom| atom.into());
        let (lbound, _) = builders_iter.size_hint();
        let offset = self.atom_builders.len();
        let indices = builders_iter.enumerate().fold(
            Vec::with_capacity(lbound),
            |mut acc, (idx, builder)| {
                self.atom_builders.insert(offset + idx, builder);
                acc.push(offset + idx);
                acc
            },
        );
        indices.into_iter()
    }

    pub fn add_atom<A: Into<AtomBuilder>>(
        &mut self,
        idx: usize,
        atom: A,
    ) -> Result<(usize, &mut AtomBuilder)> {
        let builder = atom.into();
        if self.atom_builders.contains_key(&idx) {
            return Err(DataError::DuplicateAtomIndex(idx).into());
        }
        self.atom_builders.insert(idx, builder);
        Ok((idx, self.atom_builders.get_mut(&idx).unwrap()))
    }

    pub fn add_atoms<A: Into<AtomBuilder>>(
        &mut self,
        atoms: impl IntoIterator<Item = (usize, A)>,
    ) -> Result<impl Iterator<Item = usize>> {
        let staged_atoms: Vec<(usize, AtomBuilder)> = atoms
            .into_iter()
            .map(|(idx, atom)| (idx, atom.into()))
            .collect();
        let mut seen_indices = HashSet::with_capacity(staged_atoms.len());
        for (idx, _) in &staged_atoms {
            if !seen_indices.insert(*idx) {
                return Err(DataError::DuplicateAtomIndex(*idx).into());
            }
            if self.atom_builders.contains_key(idx) {
                return Err(DataError::DuplicateAtomIndex(*idx).into());
            }
        }

        let mut indices = Vec::with_capacity(staged_atoms.len());
        for (idx, builder) in staged_atoms {
            self.atom_builders.insert(idx, builder);
            indices.push(idx);
        }

        Ok(indices.into_iter())
    }

    pub fn add_bond<B: Into<BondBuilder>>(
        &mut self,
        idx1: usize,
        idx2: usize,
        bond: B,
    ) -> Result<(usize, usize, &mut BondBuilder)> {
        if idx1 == idx2 {
            return Err(DataError::LoopBond(idx1).into());
        }
        if !self.atom_builders.contains_key(&idx1) {
            return Err(DataError::MissingAtomIndex(idx1).into());
        }
        if !self.atom_builders.contains_key(&idx2) {
            return Err(DataError::MissingAtomIndex(idx2).into());
        }
        let builder = bond.into();
        if self.bond_builders.contains_key(&(idx1, idx2)) {
            return Err(DataError::DuplicateBondIndex(idx1, idx2).into());
        }
        self.bond_builders.insert((idx1, idx2), builder);
        Ok((
            idx1,
            idx2,
            self.bond_builders.get_mut(&(idx1, idx2)).unwrap(),
        ))
    }

    pub fn add_bonds<B: Into<BondBuilder>>(
        &mut self,
        bonds: impl IntoIterator<Item = (usize, usize, B)>,
    ) -> Result<impl Iterator<Item = (usize, usize)>> {
        let canonical_bond_key = |idx1: usize, idx2: usize| (idx1.min(idx2), idx1.max(idx2));

        let staged_bonds: Vec<(usize, usize, BondBuilder)> = bonds
            .into_iter()
            .map(|(idx1, idx2, bond)| (idx1, idx2, bond.into()))
            .collect();
        let mut seen_keys = HashSet::with_capacity(staged_bonds.len());
        for (idx1, idx2, _) in &staged_bonds {
            if idx1 == idx2 {
                return Err(DataError::LoopBond(*idx1).into());
            }
            if !self.atom_builders.contains_key(idx1) {
                return Err(DataError::MissingAtomIndex(*idx1).into());
            }
            if !self.atom_builders.contains_key(idx2) {
                return Err(DataError::MissingAtomIndex(*idx2).into());
            }

            let key = canonical_bond_key(*idx1, *idx2);
            if !seen_keys.insert(key) {
                return Err(DataError::DuplicateBondIndex(*idx1, *idx2).into());
            }
            if self.bond_builders.contains_key(&key) {
                return Err(DataError::DuplicateBondIndex(*idx1, *idx2).into());
            }
        }

        let mut indices = Vec::with_capacity(staged_bonds.len());
        for (idx1, idx2, builder) in staged_bonds {
            let key = canonical_bond_key(idx1, idx2);
            self.bond_builders.insert(key, builder);
            indices.push((idx1, idx2));
        }

        Ok(indices.into_iter())
    }

    pub fn build(self) -> Result<Molecule> {
        self.build_with(
            &STRICT_ATOM_VALIDATOR,
            &STRICT_ATOM_MATCHER,
            &STRICT_BOND_MATCHER,
        )
    }

    pub fn build_with(
        self,
        atom_validator: &AtomValidator,
        atom_matcher: &AtomMatcher,
        bond_matcher: &BondMatcher,
    ) -> Result<Molecule> {
        let mut atom_builders = self.atom_builders;
        let bond_builders = self.bond_builders;
        let mut built_bonds = HashMap::with_capacity(bond_builders.len());
        let mut observed_valence: HashMap<usize, u32> = HashMap::with_capacity(atom_builders.len());

        for (key @ (idx1, idx2), bond_builder) in bond_builders {
            let bond = bond_builder.build_with(bond_matcher)?;
            let valence = u32::from(bond.order().value());
            built_bonds.insert(key, bond);
            observed_valence
                .entry(idx1)
                .and_modify(|v| *v = v.saturating_add(valence))
                .or_insert(valence);
            observed_valence
                .entry(idx2)
                .and_modify(|v| *v = v.saturating_add(valence))
                .or_insert(valence);
        }

        for (idx, observed) in observed_valence {
            if let Some(builder) = atom_builders.get_mut(&idx) {
                if builder.valence().is_none() {
                    builder.set_valence(observed);
                }
            }
        }

        let mut built_atoms = IndexMap::with_capacity(atom_builders.len());
        for (idx, atom_builder) in atom_builders {
            let atom = atom_builder.build_with(atom_validator, atom_matcher)?;
            built_atoms.insert(idx, atom);
        }

        let mut graph = StableGraph::with_capacity(built_atoms.len(), built_bonds.len());
        let mut atom_indices = HashMap::with_capacity(built_atoms.len());
        for (idx, atom) in built_atoms {
            let node_index = graph.add_node(atom);
            atom_indices.insert(idx, node_index);
        }

        for ((idx1, idx2), bond) in built_bonds {
            let node1 = *atom_indices
                .get(&idx1)
                .expect("Node index map missing mapping for idx1");
            let node2 = *atom_indices
                .get(&idx2)
                .expect("Node index map missing mapping for idx2");
            graph.add_edge(node1, node2, bond);
        }

        Ok(Molecule { data: graph })
    }
}
