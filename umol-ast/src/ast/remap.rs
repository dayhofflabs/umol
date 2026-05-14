//! AST-level index remapping produced by `MoleculeBuilder::remove`.
//!
//! Wraps `umol_graph_core::Remapping` for node/edge (atom/bond) and carries
//! sorted removed-id lists for the four relation kinds. Storage is O(removed)
//! per kind; lookups are binary search + partition-point shift.

use umol_graph_core::{EdgeId, NodeId, Remapping};

use super::idx::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};

/// Index remapping produced by `MoleculeBuilder::remove`. Translates
/// pre-removal `AtomId` / `BondId` / relation indices to post-removal
/// indices, or signals that an entity was removed. Used to rewrite stale
/// index references against the new `MoleculeAst` layout.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IdRemapping {
    graph: Remapping,
    removed_dative_bonds: Vec<u32>,
    removed_aromatic_systems: Vec<u32>,
    removed_multicenter_bonds: Vec<u32>,
    removed_noncovalent_bonds: Vec<u32>,
}

/// Inverse view of an [`IdRemapping`] for rollback. Translates surviving
/// post-removal ids back into the pre-removal coordinate system; removed ids
/// are restored from the explicit `Undo` payloads.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UndoRemapping {
    forward: IdRemapping,
}

impl IdRemapping {
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

    pub fn atom(&self, idx: AtomId) -> Option<AtomId> {
        self.graph.node(NodeId::from(idx)).map(AtomId::from)
    }

    pub fn bond(&self, idx: BondId) -> Option<BondId> {
        self.graph.edge(EdgeId::from(idx)).map(BondId::from)
    }

    pub fn dative_bond(&self, idx: DativeBondId) -> Option<DativeBondId> {
        remap_relation(&self.removed_dative_bonds, idx.0).map(DativeBondId)
    }

    pub fn aromatic_system(&self, idx: AromaticSystemId) -> Option<AromaticSystemId> {
        remap_relation(&self.removed_aromatic_systems, idx.0).map(AromaticSystemId)
    }

    pub fn multicenter_bond(&self, idx: MulticenterBondId) -> Option<MulticenterBondId> {
        remap_relation(&self.removed_multicenter_bonds, idx.0).map(MulticenterBondId)
    }

    pub fn noncovalent_bond(&self, idx: NoncovalentBondId) -> Option<NoncovalentBondId> {
        remap_relation(&self.removed_noncovalent_bonds, idx.0).map(NoncovalentBondId)
    }

    pub fn graph(&self) -> &Remapping {
        &self.graph
    }

    pub fn undo_remapping(&self) -> UndoRemapping {
        UndoRemapping {
            forward: self.clone(),
        }
    }
}

impl UndoRemapping {
    pub fn forward(&self) -> &IdRemapping {
        &self.forward
    }

    pub fn atom(&self, idx: AtomId) -> AtomId {
        AtomId(unmap_dense(&self.forward.graph.removed_nodes, idx.0))
    }

    pub fn bond(&self, idx: BondId) -> BondId {
        BondId(unmap_dense(&self.forward.graph.removed_edges, idx.0))
    }

    pub fn dative_bond(&self, idx: DativeBondId) -> DativeBondId {
        DativeBondId(unmap_dense(&self.forward.removed_dative_bonds, idx.0))
    }

    pub fn aromatic_system(&self, idx: AromaticSystemId) -> AromaticSystemId {
        AromaticSystemId(unmap_dense(&self.forward.removed_aromatic_systems, idx.0))
    }

    pub fn multicenter_bond(&self, idx: MulticenterBondId) -> MulticenterBondId {
        MulticenterBondId(unmap_dense(&self.forward.removed_multicenter_bonds, idx.0))
    }

    pub fn noncovalent_bond(&self, idx: NoncovalentBondId) -> NoncovalentBondId {
        NoncovalentBondId(unmap_dense(&self.forward.removed_noncovalent_bonds, idx.0))
    }
}

fn remap_relation(removed: &[u32], old: u32) -> Option<u32> {
    if removed.binary_search(&old).is_ok() {
        return None;
    }
    let shift = removed.partition_point(|&r| r < old);
    Some(old - shift as u32)
}

fn unmap_dense(removed: &[u32], post: u32) -> u32 {
    let mut old = post;
    loop {
        let next = post + removed.partition_point(|&r| r <= old) as u32;
        if next == old {
            return old;
        }
        old = next;
    }
}
