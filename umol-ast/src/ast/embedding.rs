//! Embedding of a molecular substructure into a parent `MoleculeAst`:
//! per-entity-type local→parent index maps borrowing the parent. Produced by
//! `MoleculeAst::induced_subgraph` and (in the future) by subgraph isomorphism
//! matching.

use std::collections::{HashMap, HashSet};

use super::edit::{AtomRef, BondRef, Edit};
use super::idx::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};
use super::molecule::MoleculeAst;

#[derive(Clone, Debug)]
pub struct MoleculeEmbedding<'a> {
    parent_atoms: Vec<AtomId>,
    parent_bonds: Vec<BondId>,
    parent_dative_bonds: Vec<DativeBondId>,
    parent_aromatic_systems: Vec<AromaticSystemId>,
    parent_multicenter_bonds: Vec<MulticenterBondId>,
    parent_noncovalent_bonds: Vec<NoncovalentBondId>,
    inverse_atom: HashMap<AtomId, u32>,
    ast: &'a MoleculeAst,
}

impl<'a> MoleculeEmbedding<'a> {
    pub(in crate::ast) fn new(
        parent_atoms: Vec<AtomId>,
        parent_bonds: Vec<BondId>,
        parent_dative_bonds: Vec<DativeBondId>,
        parent_aromatic_systems: Vec<AromaticSystemId>,
        parent_multicenter_bonds: Vec<MulticenterBondId>,
        parent_noncovalent_bonds: Vec<NoncovalentBondId>,
        inverse_atom: HashMap<AtomId, u32>,
        ast: &'a MoleculeAst,
    ) -> Self {
        Self {
            parent_atoms,
            parent_bonds,
            parent_dative_bonds,
            parent_aromatic_systems,
            parent_multicenter_bonds,
            parent_noncovalent_bonds,
            inverse_atom,
            ast,
        }
    }

    pub fn ast(&self) -> &'a MoleculeAst {
        self.ast
    }

    pub fn parent_atom(&self, local: AtomId) -> AtomId {
        self.parent_atoms[local.index()]
    }

    pub fn parent_atoms(&self) -> &[AtomId] {
        &self.parent_atoms
    }

    pub fn local_atom(&self, parent: AtomId) -> Option<AtomId> {
        self.inverse_atom.get(&parent).copied().map(AtomId)
    }

    pub fn parent_bond(&self, local: BondId) -> BondId {
        self.parent_bonds[local.index()]
    }

    pub fn parent_bonds(&self) -> &[BondId] {
        &self.parent_bonds
    }

    pub fn parent_dative_bond(&self, local: DativeBondId) -> DativeBondId {
        self.parent_dative_bonds[local.index()]
    }

    pub fn parent_dative_bonds(&self) -> &[DativeBondId] {
        &self.parent_dative_bonds
    }

    pub fn parent_aromatic_system(&self, local: AromaticSystemId) -> AromaticSystemId {
        self.parent_aromatic_systems[local.index()]
    }

    pub fn parent_aromatic_systems(&self) -> &[AromaticSystemId] {
        &self.parent_aromatic_systems
    }

    pub fn parent_multicenter_bond(&self, local: MulticenterBondId) -> MulticenterBondId {
        self.parent_multicenter_bonds[local.index()]
    }

    pub fn parent_multicenter_bonds(&self) -> &[MulticenterBondId] {
        &self.parent_multicenter_bonds
    }

    pub fn parent_noncovalent_bond(&self, local: NoncovalentBondId) -> NoncovalentBondId {
        self.parent_noncovalent_bonds[local.index()]
    }

    pub fn parent_noncovalent_bonds(&self) -> &[NoncovalentBondId] {
        &self.parent_noncovalent_bonds
    }

    /// Materialize the embedded substructure as a standalone `MoleculeAst`.
    /// Atom/bond ordering follows the parent's id order, with gaps closed by
    /// dense compaction (the order produced by `MoleculeBuilder::remove`).
    pub fn extract(&self) -> MoleculeAst {
        let atom_set: HashSet<AtomId> = self.parent_atoms.iter().copied().collect();
        let remove_atoms: Vec<AtomId> = (0..self.ast.atoms().count())
            .map(AtomId::from)
            .filter(|a| !atom_set.contains(a))
            .collect();
        let remove_bonds: Vec<BondId> = self
            .ast
            .bonds()
            .iter()
            .filter(|b| {
                let [a, b_end] = b.atom_ids();
                !atom_set.contains(&a) || !atom_set.contains(&b_end)
            })
            .map(|b| b.id)
            .collect();
        let mut builder = self.ast.edit();
        builder.remove(&remove_atoms, &remove_bonds);
        builder.build()
    }

    /// Edits that transform `parent` into the extracted substructure: a single
    /// `RemoveTopology` over atoms and bonds not present in the embedding.
    /// Overlay drops cascade automatically per phase 8i.
    pub fn edits(&self) -> Vec<Edit> {
        let surviving_atoms: HashSet<AtomId> = self.parent_atoms.iter().copied().collect();
        let surviving_bonds: HashSet<BondId> = self.parent_bonds.iter().copied().collect();
        let removed_atoms: Vec<AtomRef> = (0..self.ast.atoms().count())
            .map(AtomId::from)
            .filter(|a| !surviving_atoms.contains(a))
            .map(AtomRef::Id)
            .collect();
        let removed_bonds: Vec<BondRef> = (0..self.ast.bonds().count())
            .map(BondId::from)
            .filter(|b| !surviving_bonds.contains(b))
            .map(BondRef::Id)
            .collect();
        if removed_atoms.is_empty() && removed_bonds.is_empty() {
            return Vec::new();
        }
        vec![Edit::RemoveTopology {
            atoms: removed_atoms,
            bonds: removed_bonds,
        }]
    }
}
