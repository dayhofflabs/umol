//! Induced molecule subgraph with index maps.

use super::idx::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};
use super::molecule::MoleculeAst;

/// Result of `MoleculeAst::induced_subgraph`. The maps translate new
/// (compacted) indices back to the original molecule's indices.
#[derive(Clone, Debug)]
pub struct MoleculeSubgraph {
    pub ast: MoleculeAst,
    pub atom_map: Vec<AtomId>,
    pub bond_map: Vec<BondId>,
    pub dative_bond_map: Vec<DativeBondId>,
    pub aromatic_system_map: Vec<AromaticSystemId>,
    pub multicenter_bond_map: Vec<MulticenterBondId>,
    pub noncovalent_bond_map: Vec<NoncovalentBondId>,
}
