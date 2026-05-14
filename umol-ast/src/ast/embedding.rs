//! Embedding of one molecule into another: per-entity-type new→old index
//! maps plus the extracted sub-AST. Produced by `MoleculeAst::induced_subgraph`
//! and (in the future) by subgraph isomorphism matching.

use std::collections::HashSet;

use super::edit::{AtomRef, BondRef, Edit};
use super::idx::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};
use super::molecule::MoleculeAst;

#[derive(Clone, Debug)]
pub struct Embedding {
    ast: MoleculeAst,
    atom_map: Vec<AtomId>,
    bond_map: Vec<BondId>,
    dative_bond_map: Vec<DativeBondId>,
    aromatic_system_map: Vec<AromaticSystemId>,
    multicenter_bond_map: Vec<MulticenterBondId>,
    noncovalent_bond_map: Vec<NoncovalentBondId>,
}

impl Embedding {
    pub(in crate::ast) fn new(
        ast: MoleculeAst,
        atom_map: Vec<AtomId>,
        bond_map: Vec<BondId>,
        dative_bond_map: Vec<DativeBondId>,
        aromatic_system_map: Vec<AromaticSystemId>,
        multicenter_bond_map: Vec<MulticenterBondId>,
        noncovalent_bond_map: Vec<NoncovalentBondId>,
    ) -> Self {
        Self {
            ast,
            atom_map,
            bond_map,
            dative_bond_map,
            aromatic_system_map,
            multicenter_bond_map,
            noncovalent_bond_map,
        }
    }

    pub fn ast(&self) -> &MoleculeAst {
        &self.ast
    }

    pub fn into_ast(self) -> MoleculeAst {
        self.ast
    }

    pub fn parent_atom(&self, local: AtomId) -> AtomId {
        self.atom_map[local.index()]
    }

    pub fn parent_atoms(&self) -> &[AtomId] {
        &self.atom_map
    }

    pub fn parent_bond(&self, local: BondId) -> BondId {
        self.bond_map[local.index()]
    }

    pub fn parent_bonds(&self) -> &[BondId] {
        &self.bond_map
    }

    pub fn parent_dative_bond(&self, local: DativeBondId) -> DativeBondId {
        self.dative_bond_map[local.index()]
    }

    pub fn parent_dative_bonds(&self) -> &[DativeBondId] {
        &self.dative_bond_map
    }

    pub fn parent_aromatic_system(&self, local: AromaticSystemId) -> AromaticSystemId {
        self.aromatic_system_map[local.index()]
    }

    pub fn parent_aromatic_systems(&self) -> &[AromaticSystemId] {
        &self.aromatic_system_map
    }

    pub fn parent_multicenter_bond(&self, local: MulticenterBondId) -> MulticenterBondId {
        self.multicenter_bond_map[local.index()]
    }

    pub fn parent_multicenter_bonds(&self) -> &[MulticenterBondId] {
        &self.multicenter_bond_map
    }

    pub fn parent_noncovalent_bond(&self, local: NoncovalentBondId) -> NoncovalentBondId {
        self.noncovalent_bond_map[local.index()]
    }

    pub fn parent_noncovalent_bonds(&self) -> &[NoncovalentBondId] {
        &self.noncovalent_bond_map
    }

    /// Edits that transform `parent` into [`Self::ast`]: a single
    /// `RemoveTopology` over atoms and bonds not present in the embedding.
    /// Overlay drops cascade automatically per phase 8i.
    pub fn edits(&self, parent: &MoleculeAst) -> Vec<Edit> {
        let surviving_atoms: HashSet<AtomId> = self.atom_map.iter().copied().collect();
        let surviving_bonds: HashSet<BondId> = self.bond_map.iter().copied().collect();
        let removed_atoms: Vec<AtomRef> = (0..parent.atoms().count())
            .map(AtomId::from)
            .filter(|a| !surviving_atoms.contains(a))
            .map(AtomRef::Id)
            .collect();
        let removed_bonds: Vec<BondRef> = (0..parent.bonds().count())
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
