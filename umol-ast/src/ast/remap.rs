//! AST-level id remapping produced by `MoleculeBuilder::remove`.
//!
//! Wraps `umol_graph_core::RemovalRemapping` for node/edge (atom/bond) and carries
//! sorted removed-id lists for the six relation kinds. Storage is O(removed)
//! per kind; lookups are binary search + partition-point shift.

use umol_graph_core::{EdgeId, NodeId, RelationId, RemovalRemapping};

use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};

/// Id remapping produced by `MoleculeBuilder::remove`. Translates
/// pre-removal `AtomId` / `BondId` / relation ids to post-removal
/// ids, or signals that an entity was removed. Used to rewrite stale
/// id references against the new `MoleculeAst` layout.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IdRemapping {
    graph: RemovalRemapping,
    removed_dative_bonds: Vec<RelationId>,
    removed_aromatic_systems: Vec<RelationId>,
    removed_multicenter_bonds: Vec<RelationId>,
    removed_noncovalent_bonds: Vec<RelationId>,
    removed_stereo_atoms: Vec<RelationId>,
    removed_stereo_bonds: Vec<RelationId>,
}

/// Inverse view of an [`IdRemapping`] for rollback. Translates surviving
/// post-removal ids back into the pre-removal coordinate system; removed ids
/// are restored from the explicit `Undo` payloads.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UndoRemapping {
    forward: IdRemapping,
}

impl IdRemapping {
    pub fn empty() -> Self {
        Self::new(
            RemovalRemapping::new(Vec::new(), Vec::new()),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn relations(
        removed_dative_bonds: Vec<RelationId>,
        removed_aromatic_systems: Vec<RelationId>,
        removed_multicenter_bonds: Vec<RelationId>,
        removed_noncovalent_bonds: Vec<RelationId>,
        removed_stereo_atoms: Vec<RelationId>,
        removed_stereo_bonds: Vec<RelationId>,
    ) -> Self {
        Self::new(
            RemovalRemapping::new(Vec::new(), Vec::new()),
            removed_dative_bonds,
            removed_aromatic_systems,
            removed_multicenter_bonds,
            removed_noncovalent_bonds,
            removed_stereo_atoms,
            removed_stereo_bonds,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn new(
        graph: RemovalRemapping,
        mut removed_dative_bonds: Vec<RelationId>,
        mut removed_aromatic_systems: Vec<RelationId>,
        mut removed_multicenter_bonds: Vec<RelationId>,
        mut removed_noncovalent_bonds: Vec<RelationId>,
        mut removed_stereo_atoms: Vec<RelationId>,
        mut removed_stereo_bonds: Vec<RelationId>,
    ) -> Self {
        normalize_removed(&mut removed_dative_bonds);
        normalize_removed(&mut removed_aromatic_systems);
        normalize_removed(&mut removed_multicenter_bonds);
        normalize_removed(&mut removed_noncovalent_bonds);
        normalize_removed(&mut removed_stereo_atoms);
        normalize_removed(&mut removed_stereo_bonds);
        Self {
            graph,
            removed_dative_bonds,
            removed_aromatic_systems,
            removed_multicenter_bonds,
            removed_noncovalent_bonds,
            removed_stereo_atoms,
            removed_stereo_bonds,
        }
    }

    pub fn atom(&self, id: AtomId) -> Option<AtomId> {
        self.graph.map_node(NodeId::from(id)).map(AtomId::from)
    }

    pub fn bond(&self, id: BondId) -> Option<BondId> {
        self.graph.map_edge(EdgeId::from(id)).map(BondId::from)
    }

    pub fn dative_bond(&self, id: DativeBondId) -> Option<DativeBondId> {
        remap_relation(&self.removed_dative_bonds, id.into()).map(DativeBondId::from)
    }

    pub fn aromatic_system(&self, id: AromaticSystemId) -> Option<AromaticSystemId> {
        remap_relation(&self.removed_aromatic_systems, id.into()).map(AromaticSystemId::from)
    }

    pub fn multicenter_bond(&self, id: MulticenterBondId) -> Option<MulticenterBondId> {
        remap_relation(&self.removed_multicenter_bonds, id.into()).map(MulticenterBondId::from)
    }

    pub fn noncovalent_bond(&self, id: NoncovalentBondId) -> Option<NoncovalentBondId> {
        remap_relation(&self.removed_noncovalent_bonds, id.into()).map(NoncovalentBondId::from)
    }

    pub fn stereo_atom(&self, id: StereoAtomId) -> Option<StereoAtomId> {
        remap_relation(&self.removed_stereo_atoms, id.into()).map(StereoAtomId::from)
    }

    pub fn stereo_bond(&self, id: StereoBondId) -> Option<StereoBondId> {
        remap_relation(&self.removed_stereo_bonds, id.into()).map(StereoBondId::from)
    }

    pub fn graph(&self) -> &RemovalRemapping {
        &self.graph
    }

    pub fn undo_remapping(&self) -> UndoRemapping {
        UndoRemapping::from(self)
    }
}

impl From<&IdRemapping> for UndoRemapping {
    fn from(value: &IdRemapping) -> Self {
        Self {
            forward: value.clone(),
        }
    }
}

impl UndoRemapping {
    pub fn forward(&self) -> &IdRemapping {
        &self.forward
    }

    pub fn atom(&self, id: AtomId) -> AtomId {
        AtomId::from(self.forward.graph.unmap_node(NodeId::from(id)))
    }

    pub fn bond(&self, id: BondId) -> BondId {
        BondId::from(self.forward.graph.unmap_edge(EdgeId::from(id)))
    }

    pub fn dative_bond(&self, id: DativeBondId) -> DativeBondId {
        DativeBondId::from(unmap_dense(&self.forward.removed_dative_bonds, id.into()))
    }

    pub fn aromatic_system(&self, id: AromaticSystemId) -> AromaticSystemId {
        AromaticSystemId::from(unmap_dense(
            &self.forward.removed_aromatic_systems,
            id.into(),
        ))
    }

    pub fn multicenter_bond(&self, id: MulticenterBondId) -> MulticenterBondId {
        MulticenterBondId::from(unmap_dense(
            &self.forward.removed_multicenter_bonds,
            id.into(),
        ))
    }

    pub fn noncovalent_bond(&self, id: NoncovalentBondId) -> NoncovalentBondId {
        NoncovalentBondId::from(unmap_dense(
            &self.forward.removed_noncovalent_bonds,
            id.into(),
        ))
    }

    pub fn stereo_atom(&self, id: StereoAtomId) -> StereoAtomId {
        StereoAtomId::from(unmap_dense(&self.forward.removed_stereo_atoms, id.into()))
    }

    pub fn stereo_bond(&self, id: StereoBondId) -> StereoBondId {
        StereoBondId::from(unmap_dense(&self.forward.removed_stereo_bonds, id.into()))
    }
}

fn remap_relation(removed: &[RelationId], old: RelationId) -> Option<RelationId> {
    if removed.binary_search(&old).is_ok() {
        return None;
    }
    let shift = removed.partition_point(|&r| r < old);
    Some(RelationId(old.0 - shift as u32))
}

fn unmap_dense(removed: &[RelationId], post: RelationId) -> RelationId {
    let mut old = post;
    loop {
        let next = RelationId(post.0 + removed.partition_point(|&r| r <= old) as u32);
        if next == old {
            return old;
        }
        old = next;
    }
}

fn normalize_removed(removed: &mut Vec<RelationId>) {
    removed.sort_unstable();
    removed.dedup();
}

#[cfg(test)]
mod tests {
    use std::fmt::Debug;

    use rstest::*;
    use umol_graph_core::RemovalRemapping;

    use super::*;

    #[fixture]
    fn remapping() -> IdRemapping {
        IdRemapping::new(
            RemovalRemapping::new(vec![1, 3], vec![0, 2]),
            vec![RelationId(2), RelationId(0), RelationId(2)],
            vec![RelationId(1)],
            vec![RelationId(3), RelationId(0)],
            vec![RelationId(2)],
            vec![RelationId(1)],
            vec![RelationId(2)],
        )
    }

    #[rstest]
    #[case::before_removed(AtomId(0), Some(AtomId(0)))]
    #[case::removed_first(AtomId(1), None)]
    #[case::between_removed(AtomId(2), Some(AtomId(1)))]
    #[case::removed_second(AtomId(3), None)]
    #[case::after_removed(AtomId(4), Some(AtomId(2)))]
    fn test_id_remapping_atom(
        remapping: IdRemapping,
        #[case] input: AtomId,
        #[case] expected: Option<AtomId>,
    ) {
        assert_eq!(remapping.atom(input), expected);
    }

    #[rstest]
    #[case::before_removed(BondId(0), None)]
    #[case::between_removed(BondId(1), Some(BondId(0)))]
    #[case::removed_second(BondId(2), None)]
    #[case::after_removed(BondId(3), Some(BondId(1)))]
    fn test_id_remapping_bond(
        remapping: IdRemapping,
        #[case] input: BondId,
        #[case] expected: Option<BondId>,
    ) {
        assert_eq!(remapping.bond(input), expected);
    }

    #[rstest]
    #[case::dative_removed(DativeBondId(0), None)]
    #[case::dative_shifted(DativeBondId(1), Some(DativeBondId(0)))]
    #[case::dative_removed_duplicate_input(DativeBondId(2), None)]
    #[case::aromatic_removed(AromaticSystemId(1), None)]
    #[case::multicenter_shifted(MulticenterBondId(2), Some(MulticenterBondId(1)))]
    #[case::noncovalent_shifted(NoncovalentBondId(3), Some(NoncovalentBondId(2)))]
    #[case::stereo_atom_removed(StereoAtomId(1), None)]
    #[case::stereo_atom_shifted(StereoAtomId(2), Some(StereoAtomId(1)))]
    #[case::stereo_bond_removed(StereoBondId(2), None)]
    #[case::stereo_bond_shifted(StereoBondId(3), Some(StereoBondId(2)))]
    fn test_id_remapping_relations<T>(
        remapping: IdRemapping,
        #[case] input: T,
        #[case] expected: Option<T>,
    ) where
        T: RelationCase,
    {
        assert_eq!(input.map(&remapping), expected);
    }

    #[rstest]
    #[case::before_gap(AtomId(0), AtomId(0))]
    #[case::after_first_gap(AtomId(1), AtomId(2))]
    #[case::after_second_gap(AtomId(2), AtomId(4))]
    fn test_undo_remapping_atom(
        remapping: IdRemapping,
        #[case] input: AtomId,
        #[case] expected: AtomId,
    ) {
        assert_eq!(remapping.undo_remapping().atom(input), expected);
    }

    #[rstest]
    #[case::after_first_gap(BondId(0), BondId(1))]
    #[case::after_second_gap(BondId(1), BondId(3))]
    fn test_undo_remapping_bond(
        remapping: IdRemapping,
        #[case] input: BondId,
        #[case] expected: BondId,
    ) {
        assert_eq!(UndoRemapping::from(&remapping).bond(input), expected);
    }

    #[rstest]
    #[case::dative_after_two_gaps(DativeBondId(0), DativeBondId(1))]
    #[case::aromatic_after_gap(AromaticSystemId(0), AromaticSystemId(0))]
    #[case::aromatic_shifted(AromaticSystemId(1), AromaticSystemId(2))]
    #[case::multicenter_after_two_gaps(MulticenterBondId(1), MulticenterBondId(2))]
    #[case::noncovalent_after_gap(NoncovalentBondId(2), NoncovalentBondId(3))]
    #[case::stereo_atom_after_gap(StereoAtomId(1), StereoAtomId(2))]
    #[case::stereo_bond_after_gap(StereoBondId(2), StereoBondId(3))]
    fn test_undo_remapping_relations<T>(
        remapping: IdRemapping,
        #[case] input: T,
        #[case] expected: T,
    ) where
        T: RelationCase,
    {
        assert_eq!(input.unmap(&remapping.undo_remapping()), expected);
    }

    trait RelationCase: Copy + PartialEq + Debug {
        fn map(self, remapping: &IdRemapping) -> Option<Self>;
        fn unmap(self, remapping: &UndoRemapping) -> Self;
    }

    impl RelationCase for DativeBondId {
        fn map(self, remapping: &IdRemapping) -> Option<Self> {
            remapping.dative_bond(self)
        }

        fn unmap(self, remapping: &UndoRemapping) -> Self {
            remapping.dative_bond(self)
        }
    }

    impl RelationCase for AromaticSystemId {
        fn map(self, remapping: &IdRemapping) -> Option<Self> {
            remapping.aromatic_system(self)
        }

        fn unmap(self, remapping: &UndoRemapping) -> Self {
            remapping.aromatic_system(self)
        }
    }

    impl RelationCase for MulticenterBondId {
        fn map(self, remapping: &IdRemapping) -> Option<Self> {
            remapping.multicenter_bond(self)
        }

        fn unmap(self, remapping: &UndoRemapping) -> Self {
            remapping.multicenter_bond(self)
        }
    }

    impl RelationCase for NoncovalentBondId {
        fn map(self, remapping: &IdRemapping) -> Option<Self> {
            remapping.noncovalent_bond(self)
        }

        fn unmap(self, remapping: &UndoRemapping) -> Self {
            remapping.noncovalent_bond(self)
        }
    }

    impl RelationCase for StereoAtomId {
        fn map(self, remapping: &IdRemapping) -> Option<Self> {
            remapping.stereo_atom(self)
        }

        fn unmap(self, remapping: &UndoRemapping) -> Self {
            remapping.stereo_atom(self)
        }
    }

    impl RelationCase for StereoBondId {
        fn map(self, remapping: &IdRemapping) -> Option<Self> {
            remapping.stereo_bond(self)
        }

        fn unmap(self, remapping: &UndoRemapping) -> Self {
            remapping.stereo_bond(self)
        }
    }
}
