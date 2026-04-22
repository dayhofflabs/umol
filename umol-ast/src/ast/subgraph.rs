//! Induced molecule subgraph with index maps.

use super::idx::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};
use super::molecule::MoleculeAst;

/// Result of `MoleculeAst::induced_subgraph`. The maps translate new
/// (compacted) indices back to the original molecule's indices.
#[derive(Clone, Debug)]
pub struct MoleculeSubgraph {
    pub ast: MoleculeAst,
    pub atom_map: Vec<AtomIdx>,
    pub bond_map: Vec<BondIdx>,
    pub dative_bond_map: Vec<DativeBondIdx>,
    pub aromatic_system_map: Vec<AromaticSystemIdx>,
    pub multicenter_bond_map: Vec<MulticenterBondIdx>,
    pub noncovalent_bond_map: Vec<NoncovalentBondIdx>,
}
