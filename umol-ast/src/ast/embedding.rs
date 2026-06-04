//! Embedding of a molecular substructure into a host `MoleculeAst`:
//! per-entity-type sub→host index maps borrowing the host. Produced by
//! `MoleculeAst::induced_subgraph` and (in the future) by subgraph isomorphism
//! matching.

use std::collections::{HashMap, HashSet};

use super::edit::{AtomRef, BondRef, Edit};
use super::ids::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::molecule::MoleculeAst;

#[derive(Clone, Debug)]
pub struct MoleculeEmbedding<'a> {
    host_atoms: Vec<AtomId>,
    host_bonds: Vec<BondId>,
    host_dative_bonds: Vec<DativeBondId>,
    host_aromatic_systems: Vec<AromaticSystemId>,
    host_multicenter_bonds: Vec<MulticenterBondId>,
    host_noncovalent_bonds: Vec<NoncovalentBondId>,
    host_stereo_atoms: Vec<StereoAtomId>,
    host_stereo_bonds: Vec<StereoBondId>,
    sub_atoms: HashMap<AtomId, AtomId>,
    ast: &'a MoleculeAst,
}

impl<'a> MoleculeEmbedding<'a> {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        host_atoms: Vec<AtomId>,
        host_bonds: Vec<BondId>,
        host_dative_bonds: Vec<DativeBondId>,
        host_aromatic_systems: Vec<AromaticSystemId>,
        host_multicenter_bonds: Vec<MulticenterBondId>,
        host_noncovalent_bonds: Vec<NoncovalentBondId>,
        host_stereo_atoms: Vec<StereoAtomId>,
        host_stereo_bonds: Vec<StereoBondId>,
        sub_atoms: HashMap<AtomId, AtomId>,
        ast: &'a MoleculeAst,
    ) -> Self {
        Self {
            host_atoms,
            host_bonds,
            host_dative_bonds,
            host_aromatic_systems,
            host_multicenter_bonds,
            host_noncovalent_bonds,
            host_stereo_atoms,
            host_stereo_bonds,
            sub_atoms,
            ast,
        }
    }

    pub fn ast(&self) -> &'a MoleculeAst {
        self.ast
    }

    pub fn host_atom(&self, sub: AtomId) -> AtomId {
        self.host_atoms[sub.index()]
    }

    pub fn host_atoms(&self) -> &[AtomId] {
        &self.host_atoms
    }

    pub fn sub_atom(&self, host: AtomId) -> Option<AtomId> {
        self.sub_atoms.get(&host).copied()
    }

    pub fn host_bond(&self, sub: BondId) -> BondId {
        self.host_bonds[sub.index()]
    }

    pub fn host_bonds(&self) -> &[BondId] {
        &self.host_bonds
    }

    pub fn host_dative_bond(&self, sub: DativeBondId) -> DativeBondId {
        self.host_dative_bonds[sub.index()]
    }

    pub fn host_dative_bonds(&self) -> &[DativeBondId] {
        &self.host_dative_bonds
    }

    pub fn host_aromatic_system(&self, sub: AromaticSystemId) -> AromaticSystemId {
        self.host_aromatic_systems[sub.index()]
    }

    pub fn host_aromatic_systems(&self) -> &[AromaticSystemId] {
        &self.host_aromatic_systems
    }

    pub fn host_multicenter_bond(&self, sub: MulticenterBondId) -> MulticenterBondId {
        self.host_multicenter_bonds[sub.index()]
    }

    pub fn host_multicenter_bonds(&self) -> &[MulticenterBondId] {
        &self.host_multicenter_bonds
    }

    pub fn host_noncovalent_bond(&self, sub: NoncovalentBondId) -> NoncovalentBondId {
        self.host_noncovalent_bonds[sub.index()]
    }

    pub fn host_noncovalent_bonds(&self) -> &[NoncovalentBondId] {
        &self.host_noncovalent_bonds
    }

    pub fn host_stereo_atom(&self, sub: StereoAtomId) -> StereoAtomId {
        self.host_stereo_atoms[sub.index()]
    }

    pub fn host_stereo_atoms(&self) -> &[StereoAtomId] {
        &self.host_stereo_atoms
    }

    pub fn host_stereo_bond(&self, sub: StereoBondId) -> StereoBondId {
        self.host_stereo_bonds[sub.index()]
    }

    pub fn host_stereo_bonds(&self) -> &[StereoBondId] {
        &self.host_stereo_bonds
    }

    /// Materialize the embedded substructure as a standalone `MoleculeAst`.
    /// Atom/bond ordering follows the host's id order, with gaps closed by
    /// dense compaction (the order produced by `MoleculeBuilder::remove`).
    pub fn extract(&self) -> MoleculeAst {
        let atom_set: HashSet<AtomId> = self.host_atoms.iter().copied().collect();
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

    /// Edits that transform the host into the extracted substructure: a single
    /// `RemoveTopology` over atoms and bonds not present in the embedding.
    /// Overlay drops cascade automatically per phase 8i.
    pub fn edits(&self) -> Vec<Edit> {
        let surviving_atoms: HashSet<AtomId> = self.host_atoms.iter().copied().collect();
        let surviving_bonds: HashSet<BondId> = self.host_bonds.iter().copied().collect();
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
