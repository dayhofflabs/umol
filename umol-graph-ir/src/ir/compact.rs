//! Graph-IR id compaction across all `Molecule` entity spaces.
//!
use umol_graph_core::{Compaction, EdgeId, GraphCompaction, NodeId};

use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};

/// Molecule-level compaction produced by `MoleculeEditor::tracked_remove`. Translates a pre-removal id in
/// any of the eight entity kinds to its post-removal id, or reports that the entity was removed.
///
/// Holds one [`GraphCompaction`] for atoms and bonds and one [`Compaction`] per relation set, so
/// every entity kind is renumbered by the same operation over its own id type.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MoleculeCompaction {
    graph: GraphCompaction,
    dative_bonds: Compaction<DativeBondId>,
    aromatic_systems: Compaction<AromaticSystemId>,
    multicenter_bonds: Compaction<MulticenterBondId>,
    noncovalent_bonds: Compaction<NoncovalentBondId>,
    stereo_atoms: Compaction<StereoAtomId>,
    stereo_bonds: Compaction<StereoBondId>,
}

/// Inverse view of a [`MoleculeCompaction`] for rollback. Translates surviving
/// post-removal ids back into the pre-removal coordinate system; removed ids
/// are restored from the explicit `Undo` payloads.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UndoCompaction {
    forward: MoleculeCompaction,
}

impl MoleculeCompaction {
    /// Empty source and result domains for all eight entity kinds.
    pub const fn empty() -> Self {
        Self {
            graph: GraphCompaction::empty(),
            dative_bonds: Compaction::empty(),
            aromatic_systems: Compaction::empty(),
            multicenter_bonds: Compaction::empty(),
            noncovalent_bonds: Compaction::empty(),
            stereo_atoms: Compaction::empty(),
            stereo_bonds: Compaction::empty(),
        }
    }

    /// Assemble count-bearing compactions for all eight entity kinds.
    ///
    /// Components are already validated; agreement with a molecule is contextual.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        graph: GraphCompaction,
        dative_bonds: Compaction<DativeBondId>,
        aromatic_systems: Compaction<AromaticSystemId>,
        multicenter_bonds: Compaction<MulticenterBondId>,
        noncovalent_bonds: Compaction<NoncovalentBondId>,
        stereo_atoms: Compaction<StereoAtomId>,
        stereo_bonds: Compaction<StereoBondId>,
    ) -> Self {
        Self {
            graph,
            dative_bonds,
            aromatic_systems,
            multicenter_bonds,
            noncovalent_bonds,
            stereo_atoms,
            stereo_bonds,
        }
    }

    pub fn compact_atom(&self, id: AtomId) -> Option<AtomId> {
        self.graph.compact_node(NodeId::from(id)).map(AtomId::from)
    }

    pub fn compact_bond(&self, id: BondId) -> Option<BondId> {
        self.graph.compact_edge(EdgeId::from(id)).map(BondId::from)
    }

    pub fn compact_dative_bond(&self, id: DativeBondId) -> Option<DativeBondId> {
        self.dative_bonds.compact(id)
    }

    pub fn compact_aromatic_system(&self, id: AromaticSystemId) -> Option<AromaticSystemId> {
        self.aromatic_systems.compact(id)
    }

    pub fn compact_multicenter_bond(&self, id: MulticenterBondId) -> Option<MulticenterBondId> {
        self.multicenter_bonds.compact(id)
    }

    pub fn compact_noncovalent_bond(&self, id: NoncovalentBondId) -> Option<NoncovalentBondId> {
        self.noncovalent_bonds.compact(id)
    }

    pub fn compact_stereo_atom(&self, id: StereoAtomId) -> Option<StereoAtomId> {
        self.stereo_atoms.compact(id)
    }

    pub fn compact_stereo_bond(&self, id: StereoBondId) -> Option<StereoBondId> {
        self.stereo_bonds.compact(id)
    }

    pub fn graph(&self) -> &GraphCompaction {
        &self.graph
    }

    pub fn dative_bonds(&self) -> &Compaction<DativeBondId> {
        &self.dative_bonds
    }

    pub fn aromatic_systems(&self) -> &Compaction<AromaticSystemId> {
        &self.aromatic_systems
    }

    pub fn multicenter_bonds(&self) -> &Compaction<MulticenterBondId> {
        &self.multicenter_bonds
    }

    pub fn noncovalent_bonds(&self) -> &Compaction<NoncovalentBondId> {
        &self.noncovalent_bonds
    }

    pub fn stereo_atoms(&self) -> &Compaction<StereoAtomId> {
        &self.stereo_atoms
    }

    pub fn stereo_bonds(&self) -> &Compaction<StereoBondId> {
        &self.stereo_bonds
    }

    pub fn undo_compaction(&self) -> UndoCompaction {
        UndoCompaction::from(self)
    }
}

impl From<&MoleculeCompaction> for UndoCompaction {
    fn from(value: &MoleculeCompaction) -> Self {
        Self {
            forward: value.clone(),
        }
    }
}

impl UndoCompaction {
    pub fn forward(&self) -> &MoleculeCompaction {
        &self.forward
    }

    pub fn uncompact_atom(&self, id: AtomId) -> AtomId {
        AtomId::from(self.forward.graph.uncompact_node(NodeId::from(id)))
    }

    pub fn uncompact_bond(&self, id: BondId) -> BondId {
        BondId::from(self.forward.graph.uncompact_edge(EdgeId::from(id)))
    }

    pub fn uncompact_dative_bond(&self, id: DativeBondId) -> DativeBondId {
        self.forward.dative_bonds.uncompact(id)
    }

    pub fn uncompact_aromatic_system(&self, id: AromaticSystemId) -> AromaticSystemId {
        self.forward.aromatic_systems.uncompact(id)
    }

    pub fn uncompact_multicenter_bond(&self, id: MulticenterBondId) -> MulticenterBondId {
        self.forward.multicenter_bonds.uncompact(id)
    }

    pub fn uncompact_noncovalent_bond(&self, id: NoncovalentBondId) -> NoncovalentBondId {
        self.forward.noncovalent_bonds.uncompact(id)
    }

    pub fn uncompact_stereo_atom(&self, id: StereoAtomId) -> StereoAtomId {
        self.forward.stereo_atoms.uncompact(id)
    }

    pub fn uncompact_stereo_bond(&self, id: StereoBondId) -> StereoBondId {
        self.forward.stereo_bonds.uncompact(id)
    }
}

#[cfg(test)]
mod tests {
    use std::fmt::Debug;

    use rstest::*;
    use umol_graph_core::{Correspondence, GraphCompaction};

    use super::*;

    #[rstest]
    fn test_molecule_compaction_empty() {
        let actual = MoleculeCompaction::empty();
        assert_eq!(actual.graph(), &GraphCompaction::empty());
        assert_eq!(actual.dative_bonds(), &Compaction::empty());
        assert_eq!(actual.aromatic_systems(), &Compaction::empty());
        assert_eq!(actual.multicenter_bonds(), &Compaction::empty());
        assert_eq!(actual.noncovalent_bonds(), &Compaction::empty());
        assert_eq!(actual.stereo_atoms(), &Compaction::empty());
        assert_eq!(actual.stereo_bonds(), &Compaction::empty());
    }
    use crate::ir::correspondence::MoleculeCorrespondence;

    #[fixture]
    fn compaction() -> MoleculeCompaction {
        MoleculeCompaction::new(
            GraphCompaction::new(
                Compaction::new(5, vec![NodeId(1), NodeId(3)]).unwrap(),
                Compaction::new(4, vec![EdgeId(0), EdgeId(2)]).unwrap(),
            ),
            Compaction::new(3, vec![DativeBondId(2), DativeBondId(0), DativeBondId(2)]).unwrap(),
            Compaction::new(3, vec![AromaticSystemId(1)]).unwrap(),
            Compaction::new(4, vec![MulticenterBondId(3), MulticenterBondId(0)]).unwrap(),
            Compaction::new(4, vec![NoncovalentBondId(2)]).unwrap(),
            Compaction::new(3, vec![StereoAtomId(1)]).unwrap(),
            Compaction::new(4, vec![StereoBondId(2)]).unwrap(),
        )
    }

    #[rstest]
    fn test_molecule_correspondence_from_compaction(compaction: MoleculeCompaction) {
        let expected = MoleculeCorrespondence::new(
            Correspondence::new(
                vec![
                    (AtomId(0), AtomId(0)),
                    (AtomId(2), AtomId(1)),
                    (AtomId(4), AtomId(2)),
                ],
                5,
                3,
            )
            .unwrap(),
            Correspondence::new(vec![(BondId(1), BondId(0)), (BondId(3), BondId(1))], 4, 2)
                .unwrap(),
            Correspondence::new(vec![(DativeBondId(1), DativeBondId(0))], 3, 1).unwrap(),
            Correspondence::new(
                vec![
                    (AromaticSystemId(0), AromaticSystemId(0)),
                    (AromaticSystemId(2), AromaticSystemId(1)),
                ],
                3,
                2,
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (MulticenterBondId(1), MulticenterBondId(0)),
                    (MulticenterBondId(2), MulticenterBondId(1)),
                ],
                4,
                2,
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (NoncovalentBondId(0), NoncovalentBondId(0)),
                    (NoncovalentBondId(1), NoncovalentBondId(1)),
                    (NoncovalentBondId(3), NoncovalentBondId(2)),
                ],
                4,
                3,
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (StereoAtomId(0), StereoAtomId(0)),
                    (StereoAtomId(2), StereoAtomId(1)),
                ],
                3,
                2,
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (StereoBondId(0), StereoBondId(0)),
                    (StereoBondId(1), StereoBondId(1)),
                    (StereoBondId(3), StereoBondId(2)),
                ],
                4,
                3,
            )
            .unwrap(),
        );
        assert_eq!(MoleculeCorrespondence::from(&compaction), expected);
    }

    #[rstest]
    fn test_undo_compaction_from(compaction: MoleculeCompaction) {
        let undo = UndoCompaction::from(&compaction);
        assert_eq!(undo.forward(), &compaction);
        assert_eq!(
            [
                (
                    undo.forward().graph().nodes().source_count(),
                    undo.forward().graph().nodes().result_count()
                ),
                (
                    undo.forward().graph().edges().source_count(),
                    undo.forward().graph().edges().result_count()
                ),
                (
                    undo.forward().dative_bonds().source_count(),
                    undo.forward().dative_bonds().result_count()
                ),
                (
                    undo.forward().aromatic_systems().source_count(),
                    undo.forward().aromatic_systems().result_count()
                ),
                (
                    undo.forward().multicenter_bonds().source_count(),
                    undo.forward().multicenter_bonds().result_count()
                ),
                (
                    undo.forward().noncovalent_bonds().source_count(),
                    undo.forward().noncovalent_bonds().result_count()
                ),
                (
                    undo.forward().stereo_atoms().source_count(),
                    undo.forward().stereo_atoms().result_count()
                ),
                (
                    undo.forward().stereo_bonds().source_count(),
                    undo.forward().stereo_bonds().result_count()
                ),
            ],
            [
                (5, 3),
                (4, 2),
                (3, 1),
                (3, 2),
                (4, 2),
                (4, 3),
                (3, 2),
                (4, 3)
            ],
        );
    }

    #[rstest]
    #[case::before_removed(AtomId(0), Some(AtomId(0)))]
    #[case::removed_first(AtomId(1), None)]
    #[case::between_removed(AtomId(2), Some(AtomId(1)))]
    #[case::removed_second(AtomId(3), None)]
    #[case::after_removed(AtomId(4), Some(AtomId(2)))]
    fn test_molecule_compaction_compact_atom(
        compaction: MoleculeCompaction,
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
    fn test_molecule_compaction_compact_bond(
        compaction: MoleculeCompaction,
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
    fn test_molecule_compaction_compact_relations<T>(
        compaction: MoleculeCompaction,
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
    fn test_undo_compaction_uncompact_atom(
        compaction: MoleculeCompaction,
        #[case] input: AtomId,
        #[case] expected: AtomId,
    ) {
        assert_eq!(compaction.undo_compaction().uncompact_atom(input), expected);
    }

    #[rstest]
    #[case::after_first_gap(BondId(0), BondId(1))]
    #[case::after_second_gap(BondId(1), BondId(3))]
    fn test_undo_compaction_uncompact_bond(
        compaction: MoleculeCompaction,
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
    fn test_undo_compaction_uncompact_relations<T>(
        compaction: MoleculeCompaction,
        #[case] input: T,
        #[case] expected: T,
    ) where
        T: RelationCase,
    {
        assert_eq!(input.uncompact(&compaction.undo_compaction()), expected);
    }

    trait RelationCase: Copy + PartialEq + Debug {
        fn compact(self, compaction: &MoleculeCompaction) -> Option<Self>;
        fn uncompact(self, compaction: &UndoCompaction) -> Self;
    }

    impl RelationCase for DativeBondId {
        fn compact(self, compaction: &MoleculeCompaction) -> Option<Self> {
            compaction.compact_dative_bond(self)
        }

        fn uncompact(self, compaction: &UndoCompaction) -> Self {
            compaction.uncompact_dative_bond(self)
        }
    }

    impl RelationCase for AromaticSystemId {
        fn compact(self, compaction: &MoleculeCompaction) -> Option<Self> {
            compaction.compact_aromatic_system(self)
        }

        fn uncompact(self, compaction: &UndoCompaction) -> Self {
            compaction.uncompact_aromatic_system(self)
        }
    }

    impl RelationCase for MulticenterBondId {
        fn compact(self, compaction: &MoleculeCompaction) -> Option<Self> {
            compaction.compact_multicenter_bond(self)
        }

        fn uncompact(self, compaction: &UndoCompaction) -> Self {
            compaction.uncompact_multicenter_bond(self)
        }
    }

    impl RelationCase for NoncovalentBondId {
        fn compact(self, compaction: &MoleculeCompaction) -> Option<Self> {
            compaction.compact_noncovalent_bond(self)
        }

        fn uncompact(self, compaction: &UndoCompaction) -> Self {
            compaction.uncompact_noncovalent_bond(self)
        }
    }

    impl RelationCase for StereoAtomId {
        fn compact(self, compaction: &MoleculeCompaction) -> Option<Self> {
            compaction.compact_stereo_atom(self)
        }

        fn uncompact(self, compaction: &UndoCompaction) -> Self {
            compaction.uncompact_stereo_atom(self)
        }
    }

    impl RelationCase for StereoBondId {
        fn compact(self, compaction: &MoleculeCompaction) -> Option<Self> {
            compaction.compact_stereo_bond(self)
        }

        fn uncompact(self, compaction: &UndoCompaction) -> Self {
            compaction.uncompact_stereo_bond(self)
        }
    }
}
