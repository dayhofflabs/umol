//! AST-level id mappings between `MoleculeAst` id spaces.
//!
//! [`IdCompaction`] is the removal compaction produced by
//! `MoleculeBuilder::remove` (wraps `umol_graph_core::Compaction` for atom/bond
//! and carries sorted removed-id lists for the six relation kinds; lookups are
//! binary search + partition-point shift). [`IdRemapping`] is the general total
//! relabeling used to move `Delta`s between id spaces (`reverse`, `compose`).

use std::collections::HashMap;

use umol_graph_core::{Compaction, EdgeId, NodeId, RelationId};

use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};

/// Id compaction produced by `MoleculeBuilder::remove`. Translates
/// pre-removal `AtomId` / `BondId` / relation ids to post-removal
/// ids, or signals that an entity was removed. Used to rewrite stale
/// id references against the new `MoleculeAst` layout.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IdCompaction {
    graph: Compaction,
    removed_dative_bonds: Vec<RelationId>,
    removed_aromatic_systems: Vec<RelationId>,
    removed_multicenter_bonds: Vec<RelationId>,
    removed_noncovalent_bonds: Vec<RelationId>,
    removed_stereo_atoms: Vec<RelationId>,
    removed_stereo_bonds: Vec<RelationId>,
}

/// Inverse view of an [`IdCompaction`] for rollback. Translates surviving
/// post-removal ids back into the pre-removal coordinate system; removed ids
/// are restored from the explicit `Undo` payloads.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UndoCompaction {
    forward: IdCompaction,
}

impl IdCompaction {
    pub fn empty() -> Self {
        Self::new(
            Compaction::new(Vec::new(), Vec::new()),
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
            Compaction::new(Vec::new(), Vec::new()),
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
        graph: Compaction,
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

    pub fn compact_atom(&self, id: AtomId) -> Option<AtomId> {
        self.graph.compact_node(NodeId::from(id)).map(AtomId::from)
    }

    pub fn compact_bond(&self, id: BondId) -> Option<BondId> {
        self.graph.compact_edge(EdgeId::from(id)).map(BondId::from)
    }

    pub fn compact_dative_bond(&self, id: DativeBondId) -> Option<DativeBondId> {
        compact_relation(&self.removed_dative_bonds, id.into()).map(DativeBondId::from)
    }

    pub fn compact_aromatic_system(&self, id: AromaticSystemId) -> Option<AromaticSystemId> {
        compact_relation(&self.removed_aromatic_systems, id.into()).map(AromaticSystemId::from)
    }

    pub fn compact_multicenter_bond(&self, id: MulticenterBondId) -> Option<MulticenterBondId> {
        compact_relation(&self.removed_multicenter_bonds, id.into()).map(MulticenterBondId::from)
    }

    pub fn compact_noncovalent_bond(&self, id: NoncovalentBondId) -> Option<NoncovalentBondId> {
        compact_relation(&self.removed_noncovalent_bonds, id.into()).map(NoncovalentBondId::from)
    }

    pub fn compact_stereo_atom(&self, id: StereoAtomId) -> Option<StereoAtomId> {
        compact_relation(&self.removed_stereo_atoms, id.into()).map(StereoAtomId::from)
    }

    pub fn compact_stereo_bond(&self, id: StereoBondId) -> Option<StereoBondId> {
        compact_relation(&self.removed_stereo_bonds, id.into()).map(StereoBondId::from)
    }

    pub fn graph(&self) -> &Compaction {
        &self.graph
    }

    pub fn undo_compaction(&self) -> UndoCompaction {
        UndoCompaction::from(self)
    }
}

impl From<&IdCompaction> for UndoCompaction {
    fn from(value: &IdCompaction) -> Self {
        Self {
            forward: value.clone(),
        }
    }
}

impl UndoCompaction {
    pub fn forward(&self) -> &IdCompaction {
        &self.forward
    }

    pub fn uncompact_atom(&self, id: AtomId) -> AtomId {
        AtomId::from(self.forward.graph.uncompact_node(NodeId::from(id)))
    }

    pub fn uncompact_bond(&self, id: BondId) -> BondId {
        BondId::from(self.forward.graph.uncompact_edge(EdgeId::from(id)))
    }

    pub fn uncompact_dative_bond(&self, id: DativeBondId) -> DativeBondId {
        DativeBondId::from(uncompact_dense(
            &self.forward.removed_dative_bonds,
            id.into(),
        ))
    }

    pub fn uncompact_aromatic_system(&self, id: AromaticSystemId) -> AromaticSystemId {
        AromaticSystemId::from(uncompact_dense(
            &self.forward.removed_aromatic_systems,
            id.into(),
        ))
    }

    pub fn uncompact_multicenter_bond(&self, id: MulticenterBondId) -> MulticenterBondId {
        MulticenterBondId::from(uncompact_dense(
            &self.forward.removed_multicenter_bonds,
            id.into(),
        ))
    }

    pub fn uncompact_noncovalent_bond(&self, id: NoncovalentBondId) -> NoncovalentBondId {
        NoncovalentBondId::from(uncompact_dense(
            &self.forward.removed_noncovalent_bonds,
            id.into(),
        ))
    }

    pub fn uncompact_stereo_atom(&self, id: StereoAtomId) -> StereoAtomId {
        StereoAtomId::from(uncompact_dense(
            &self.forward.removed_stereo_atoms,
            id.into(),
        ))
    }

    pub fn uncompact_stereo_bond(&self, id: StereoBondId) -> StereoBondId {
        StereoBondId::from(uncompact_dense(
            &self.forward.removed_stereo_bonds,
            id.into(),
        ))
    }
}

fn compact_relation(removed: &[RelationId], old: RelationId) -> Option<RelationId> {
    if removed.binary_search(&old).is_ok() {
        return None;
    }
    let shift = removed.partition_point(|&r| r < old);
    Some(RelationId(old.0 - shift as u32))
}

fn uncompact_dense(removed: &[RelationId], post: RelationId) -> RelationId {
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

/// Total id relabeling between two `MoleculeAst` id spaces — the general
/// counterpart to [`IdCompaction`]. Maps every referenced atom / bond / overlay
/// id to its image in the target id space; used to move `Delta`s between id spaces
/// (`reverse`, `compose`). Every id a moved delta references must be present.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct IdRemapping {
    atom: HashMap<AtomId, AtomId>,
    bond: HashMap<BondId, BondId>,
    dative: HashMap<DativeBondId, DativeBondId>,
    aromatic: HashMap<AromaticSystemId, AromaticSystemId>,
    multicenter: HashMap<MulticenterBondId, MulticenterBondId>,
    noncovalent: HashMap<NoncovalentBondId, NoncovalentBondId>,
    stereo_atom: HashMap<StereoAtomId, StereoAtomId>,
    stereo_bond: HashMap<StereoBondId, StereoBondId>,
}

impl IdRemapping {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        atom: HashMap<AtomId, AtomId>,
        bond: HashMap<BondId, BondId>,
        dative: HashMap<DativeBondId, DativeBondId>,
        aromatic: HashMap<AromaticSystemId, AromaticSystemId>,
        multicenter: HashMap<MulticenterBondId, MulticenterBondId>,
        noncovalent: HashMap<NoncovalentBondId, NoncovalentBondId>,
        stereo_atom: HashMap<StereoAtomId, StereoAtomId>,
        stereo_bond: HashMap<StereoBondId, StereoBondId>,
    ) -> Self {
        Self {
            atom,
            bond,
            dative,
            aromatic,
            multicenter,
            noncovalent,
            stereo_atom,
            stereo_bond,
        }
    }

    pub fn map_atom(&self, id: AtomId) -> AtomId {
        self.atom[&id]
    }

    pub fn map_bond(&self, id: BondId) -> BondId {
        self.bond[&id]
    }

    pub fn map_dative(&self, id: DativeBondId) -> DativeBondId {
        self.dative[&id]
    }

    pub fn map_aromatic(&self, id: AromaticSystemId) -> AromaticSystemId {
        self.aromatic[&id]
    }

    pub fn map_multicenter(&self, id: MulticenterBondId) -> MulticenterBondId {
        self.multicenter[&id]
    }

    pub fn map_noncovalent(&self, id: NoncovalentBondId) -> NoncovalentBondId {
        self.noncovalent[&id]
    }

    pub fn map_stereo_atom(&self, id: StereoAtomId) -> StereoAtomId {
        self.stereo_atom[&id]
    }

    pub fn map_stereo_bond(&self, id: StereoBondId) -> StereoBondId {
        self.stereo_bond[&id]
    }
}

#[cfg(test)]
mod tests {
    use std::fmt::Debug;

    use rstest::*;
    use umol_graph_core::Compaction;

    use super::*;

    #[fixture]
    fn compaction() -> IdCompaction {
        IdCompaction::new(
            Compaction::new(vec![1, 3], vec![0, 2]),
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
    fn test_id_compaction_atom(
        compaction: IdCompaction,
        #[case] input: AtomId,
        #[case] expected: Option<AtomId>,
    ) {
        assert_eq!(compaction.compact_atom(input), expected);
    }

    #[rstest]
    #[case::before_removed(BondId(0), None)]
    #[case::between_removed(BondId(1), Some(BondId(0)))]
    #[case::removed_second(BondId(2), None)]
    #[case::after_removed(BondId(3), Some(BondId(1)))]
    fn test_id_compaction_bond(
        compaction: IdCompaction,
        #[case] input: BondId,
        #[case] expected: Option<BondId>,
    ) {
        assert_eq!(compaction.compact_bond(input), expected);
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
    fn test_id_compaction_relations<T>(
        compaction: IdCompaction,
        #[case] input: T,
        #[case] expected: Option<T>,
    ) where
        T: RelationCase,
    {
        assert_eq!(input.compact(&compaction), expected);
    }

    #[rstest]
    #[case::before_gap(AtomId(0), AtomId(0))]
    #[case::after_first_gap(AtomId(1), AtomId(2))]
    #[case::after_second_gap(AtomId(2), AtomId(4))]
    fn test_undo_compaction_atom(
        compaction: IdCompaction,
        #[case] input: AtomId,
        #[case] expected: AtomId,
    ) {
        assert_eq!(compaction.undo_compaction().uncompact_atom(input), expected);
    }

    #[rstest]
    #[case::after_first_gap(BondId(0), BondId(1))]
    #[case::after_second_gap(BondId(1), BondId(3))]
    fn test_undo_compaction_bond(
        compaction: IdCompaction,
        #[case] input: BondId,
        #[case] expected: BondId,
    ) {
        assert_eq!(
            UndoCompaction::from(&compaction).uncompact_bond(input),
            expected
        );
    }

    #[rstest]
    #[case::dative_after_two_gaps(DativeBondId(0), DativeBondId(1))]
    #[case::aromatic_after_gap(AromaticSystemId(0), AromaticSystemId(0))]
    #[case::aromatic_shifted(AromaticSystemId(1), AromaticSystemId(2))]
    #[case::multicenter_after_two_gaps(MulticenterBondId(1), MulticenterBondId(2))]
    #[case::noncovalent_after_gap(NoncovalentBondId(2), NoncovalentBondId(3))]
    #[case::stereo_atom_after_gap(StereoAtomId(1), StereoAtomId(2))]
    #[case::stereo_bond_after_gap(StereoBondId(2), StereoBondId(3))]
    fn test_undo_compaction_relations<T>(
        compaction: IdCompaction,
        #[case] input: T,
        #[case] expected: T,
    ) where
        T: RelationCase,
    {
        assert_eq!(input.uncompact(&compaction.undo_compaction()), expected);
    }

    trait RelationCase: Copy + PartialEq + Debug {
        fn compact(self, compaction: &IdCompaction) -> Option<Self>;
        fn uncompact(self, compaction: &UndoCompaction) -> Self;
    }

    impl RelationCase for DativeBondId {
        fn compact(self, compaction: &IdCompaction) -> Option<Self> {
            compaction.compact_dative_bond(self)
        }

        fn uncompact(self, compaction: &UndoCompaction) -> Self {
            compaction.uncompact_dative_bond(self)
        }
    }

    impl RelationCase for AromaticSystemId {
        fn compact(self, compaction: &IdCompaction) -> Option<Self> {
            compaction.compact_aromatic_system(self)
        }

        fn uncompact(self, compaction: &UndoCompaction) -> Self {
            compaction.uncompact_aromatic_system(self)
        }
    }

    impl RelationCase for MulticenterBondId {
        fn compact(self, compaction: &IdCompaction) -> Option<Self> {
            compaction.compact_multicenter_bond(self)
        }

        fn uncompact(self, compaction: &UndoCompaction) -> Self {
            compaction.uncompact_multicenter_bond(self)
        }
    }

    impl RelationCase for NoncovalentBondId {
        fn compact(self, compaction: &IdCompaction) -> Option<Self> {
            compaction.compact_noncovalent_bond(self)
        }

        fn uncompact(self, compaction: &UndoCompaction) -> Self {
            compaction.uncompact_noncovalent_bond(self)
        }
    }

    impl RelationCase for StereoAtomId {
        fn compact(self, compaction: &IdCompaction) -> Option<Self> {
            compaction.compact_stereo_atom(self)
        }

        fn uncompact(self, compaction: &UndoCompaction) -> Self {
            compaction.uncompact_stereo_atom(self)
        }
    }

    impl RelationCase for StereoBondId {
        fn compact(self, compaction: &IdCompaction) -> Option<Self> {
            compaction.compact_stereo_bond(self)
        }

        fn uncompact(self, compaction: &UndoCompaction) -> Self {
            compaction.uncompact_stereo_bond(self)
        }
    }
}
