//! AST-level index remapping produced by `MoleculeBuilder::remove`.
//!
//! Wraps `umol_graph_core::Remapping` for node/edge (atom/bond) and carries
//! sorted removed-id lists for the four relation kinds. Storage is O(removed)
//! per kind; lookups are binary search + partition-point shift.

use umol_graph_core::{EdgeId, NodeId, Remapping};

use super::idx::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};

/// Index remapping produced by `MoleculeBuilder::remove`. Translates
/// pre-removal `AtomIdx` / `BondIdx` / relation indices to post-removal
/// indices, or signals that an entity was removed. Used to rewrite stale
/// index references against the new `MoleculeAst` layout.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IdxRemapping {
    graph: Remapping,
    removed_dative_bonds: Vec<u32>,
    removed_aromatic_systems: Vec<u32>,
    removed_multicenter_bonds: Vec<u32>,
    removed_noncovalent_bonds: Vec<u32>,
}

impl IdxRemapping {
    pub fn new(
        graph: Remapping,
        removed_dative_bonds: Vec<u32>,
        removed_aromatic_systems: Vec<u32>,
        removed_multicenter_bonds: Vec<u32>,
        removed_noncovalent_bonds: Vec<u32>,
    ) -> Self {
        Self {
            graph,
            removed_dative_bonds,
            removed_aromatic_systems,
            removed_multicenter_bonds,
            removed_noncovalent_bonds,
        }
    }

    pub fn atom(&self, idx: AtomIdx) -> Option<AtomIdx> {
        self.graph.node(NodeId::from(idx)).map(AtomIdx::from)
    }

    pub fn bond(&self, idx: BondIdx) -> Option<BondIdx> {
        self.graph.edge(EdgeId::from(idx)).map(BondIdx::from)
    }

    pub fn dative_bond(&self, idx: DativeBondIdx) -> Option<DativeBondIdx> {
        remap_relation(&self.removed_dative_bonds, idx.0).map(DativeBondIdx)
    }

    pub fn aromatic_system(&self, idx: AromaticSystemIdx) -> Option<AromaticSystemIdx> {
        remap_relation(&self.removed_aromatic_systems, idx.0).map(AromaticSystemIdx)
    }

    pub fn multicenter_bond(&self, idx: MulticenterBondIdx) -> Option<MulticenterBondIdx> {
        remap_relation(&self.removed_multicenter_bonds, idx.0).map(MulticenterBondIdx)
    }

    pub fn noncovalent_bond(&self, idx: NoncovalentBondIdx) -> Option<NoncovalentBondIdx> {
        remap_relation(&self.removed_noncovalent_bonds, idx.0).map(NoncovalentBondIdx)
    }

    pub fn graph(&self) -> &Remapping {
        &self.graph
    }
}

fn remap_relation(removed: &[u32], old: u32) -> Option<u32> {
    if removed.binary_search(&old).is_ok() {
        return None;
    }
    let shift = removed.partition_point(|&r| r < old);
    Some(old - shift as u32)
}
