//! Graph-IR id remappings between `Molecule` id spaces.
//!
//! [`MoleculeRemapping`] is the bijective renumbering over all eight molecule entity id spaces.

use umol_graph_core::{EdgeId, GraphRemapping, NodeId, Remapping};

use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};

/// Independent bijective renumberings of all eight molecule entity id spaces.
/// Each component declares equal source and target counts. Agreement with a particular molecule
/// is a contextual requirement of the operation consuming this carrier.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MoleculeRemapping {
    graph: GraphRemapping,
    dative_bonds: Remapping<DativeBondId>,
    aromatic_systems: Remapping<AromaticSystemId>,
    multicenter_bonds: Remapping<MulticenterBondId>,
    noncovalent_bonds: Remapping<NoncovalentBondId>,
    stereo_atoms: Remapping<StereoAtomId>,
    stereo_bonds: Remapping<StereoBondId>,
}

impl MoleculeRemapping {
    /// Empty source and result domains for all eight entity kinds.
    pub const fn empty() -> Self {
        Self {
            graph: GraphRemapping::empty(),
            dative_bonds: Remapping::empty(),
            aromatic_systems: Remapping::empty(),
            multicenter_bonds: Remapping::empty(),
            noncovalent_bonds: Remapping::empty(),
            stereo_atoms: Remapping::empty(),
            stereo_bonds: Remapping::empty(),
        }
    }

    /// Assemble already-valid graph and overlay permutations.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        graph: GraphRemapping,
        dative_bonds: Remapping<DativeBondId>,
        aromatic_systems: Remapping<AromaticSystemId>,
        multicenter_bonds: Remapping<MulticenterBondId>,
        noncovalent_bonds: Remapping<NoncovalentBondId>,
        stereo_atoms: Remapping<StereoAtomId>,
        stereo_bonds: Remapping<StereoBondId>,
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

    /// Return the node and edge remappings.
    pub fn graph(&self) -> &GraphRemapping {
        &self.graph
    }

    /// Return the dative-bond remapping.
    pub fn dative_bonds(&self) -> &Remapping<DativeBondId> {
        &self.dative_bonds
    }

    /// Return the aromatic-system remapping.
    pub fn aromatic_systems(&self) -> &Remapping<AromaticSystemId> {
        &self.aromatic_systems
    }

    /// Return the multicenter-bond remapping.
    pub fn multicenter_bonds(&self) -> &Remapping<MulticenterBondId> {
        &self.multicenter_bonds
    }

    /// Return the noncovalent-bond remapping.
    pub fn noncovalent_bonds(&self) -> &Remapping<NoncovalentBondId> {
        &self.noncovalent_bonds
    }

    /// Return the stereo-atom remapping.
    pub fn stereo_atoms(&self) -> &Remapping<StereoAtomId> {
        &self.stereo_atoms
    }

    /// Return the stereo-bond remapping.
    pub fn stereo_bonds(&self) -> &Remapping<StereoBondId> {
        &self.stereo_bonds
    }

    /// Return the image of an atom id, or `None` when it lies outside the source domain.
    pub fn try_map_atom(&self, id: AtomId) -> Option<AtomId> {
        self.graph.try_map_node(NodeId::from(id)).map(AtomId::from)
    }

    /// Return the image of a bond id, or `None` when it lies outside the source domain.
    pub fn try_map_bond(&self, id: BondId) -> Option<BondId> {
        self.graph.try_map_edge(EdgeId::from(id)).map(BondId::from)
    }

    /// Return the image of a dative-bond id, or `None` when it lies outside the source domain.
    pub fn try_map_dative_bond(&self, id: DativeBondId) -> Option<DativeBondId> {
        self.dative_bonds.try_map(id)
    }

    /// Return the image of an aromatic-system id, or `None` when it lies outside the source domain.
    pub fn try_map_aromatic_system(&self, id: AromaticSystemId) -> Option<AromaticSystemId> {
        self.aromatic_systems.try_map(id)
    }

    /// Return the image of a multicenter-bond id, or `None` when it lies outside the source domain.
    pub fn try_map_multicenter_bond(&self, id: MulticenterBondId) -> Option<MulticenterBondId> {
        self.multicenter_bonds.try_map(id)
    }

    /// Return the image of a noncovalent-bond id, or `None` when it lies outside the source domain.
    pub fn try_map_noncovalent_bond(&self, id: NoncovalentBondId) -> Option<NoncovalentBondId> {
        self.noncovalent_bonds.try_map(id)
    }

    /// Return the image of a stereo-atom id, or `None` when it lies outside the source domain.
    pub fn try_map_stereo_atom(&self, id: StereoAtomId) -> Option<StereoAtomId> {
        self.stereo_atoms.try_map(id)
    }

    /// Return the image of a stereo-bond id, or `None` when it lies outside the source domain.
    pub fn try_map_stereo_bond(&self, id: StereoBondId) -> Option<StereoBondId> {
        self.stereo_bonds.try_map(id)
    }

    /// Return the image of an atom id.
    ///
    /// # Panics
    ///
    /// Panics when `id` lies outside the atom source domain.
    pub fn map_atom(&self, id: AtomId) -> AtomId {
        self.try_map_atom(id)
            .expect("atom id outside remapping source domain")
    }

    /// Return the image of a bond id.
    ///
    /// # Panics
    ///
    /// Panics when `id` lies outside the bond source domain.
    pub fn map_bond(&self, id: BondId) -> BondId {
        self.try_map_bond(id)
            .expect("bond id outside remapping source domain")
    }

    /// Return the image of a dative-bond id.
    ///
    /// # Panics
    ///
    /// Panics when `id` lies outside the dative-bond source domain.
    pub fn map_dative_bond(&self, id: DativeBondId) -> DativeBondId {
        self.try_map_dative_bond(id)
            .expect("dative-bond id outside remapping source domain")
    }

    /// Return the image of an aromatic-system id.
    ///
    /// # Panics
    ///
    /// Panics when `id` lies outside the aromatic-system source domain.
    pub fn map_aromatic_system(&self, id: AromaticSystemId) -> AromaticSystemId {
        self.try_map_aromatic_system(id)
            .expect("aromatic-system id outside remapping source domain")
    }

    /// Return the image of a multicenter-bond id.
    ///
    /// # Panics
    ///
    /// Panics when `id` lies outside the multicenter-bond source domain.
    pub fn map_multicenter_bond(&self, id: MulticenterBondId) -> MulticenterBondId {
        self.try_map_multicenter_bond(id)
            .expect("multicenter-bond id outside remapping source domain")
    }

    /// Return the image of a noncovalent-bond id.
    ///
    /// # Panics
    ///
    /// Panics when `id` lies outside the noncovalent-bond source domain.
    pub fn map_noncovalent_bond(&self, id: NoncovalentBondId) -> NoncovalentBondId {
        self.try_map_noncovalent_bond(id)
            .expect("noncovalent-bond id outside remapping source domain")
    }

    /// Return the image of a stereo-atom id.
    ///
    /// # Panics
    ///
    /// Panics when `id` lies outside the stereo-atom source domain.
    pub fn map_stereo_atom(&self, id: StereoAtomId) -> StereoAtomId {
        self.try_map_stereo_atom(id)
            .expect("stereo-atom id outside remapping source domain")
    }

    /// Return the image of a stereo-bond id.
    ///
    /// # Panics
    ///
    /// Panics when `id` lies outside the stereo-bond source domain.
    pub fn map_stereo_bond(&self, id: StereoBondId) -> StereoBondId {
        self.try_map_stereo_bond(id)
            .expect("stereo-bond id outside remapping source domain")
    }
}

impl Default for MoleculeRemapping {
    fn default() -> Self {
        Self::empty()
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::Correspondence;

    use super::*;

    #[rstest]
    fn test_molecule_remapping_empty() {
        let actual = MoleculeRemapping::empty();
        assert_eq!(actual.graph(), &GraphRemapping::empty());
        assert_eq!(actual.dative_bonds(), &Remapping::empty());
        assert_eq!(actual.aromatic_systems(), &Remapping::empty());
        assert_eq!(actual.multicenter_bonds(), &Remapping::empty());
        assert_eq!(actual.noncovalent_bonds(), &Remapping::empty());
        assert_eq!(actual.stereo_atoms(), &Remapping::empty());
        assert_eq!(actual.stereo_bonds(), &Remapping::empty());
    }
    use crate::ir::MoleculeCorrespondence;

    #[fixture]
    fn molecule_remapping() -> MoleculeRemapping {
        MoleculeRemapping::new(
            GraphRemapping::new(
                Remapping::new(vec![NodeId(1), NodeId(0)]).unwrap(),
                Remapping::new(vec![EdgeId(1), EdgeId(0)]).unwrap(),
            ),
            Remapping::new(vec![DativeBondId(1), DativeBondId(0)]).unwrap(),
            Remapping::new(vec![AromaticSystemId(0)]).unwrap(),
            Remapping::new(vec![MulticenterBondId(1), MulticenterBondId(0)]).unwrap(),
            Remapping::new(vec![NoncovalentBondId(0)]).unwrap(),
            Remapping::new(vec![StereoAtomId(1), StereoAtomId(0)]).unwrap(),
            Remapping::new(vec![StereoBondId(0)]).unwrap(),
        )
    }

    #[rstest]
    fn test_molecule_correspondence_from_remapping(molecule_remapping: MoleculeRemapping) {
        let expected = MoleculeCorrespondence::new(
            Correspondence::new(vec![(AtomId(0), AtomId(1)), (AtomId(1), AtomId(0))], 2, 2)
                .unwrap(),
            Correspondence::new(vec![(BondId(0), BondId(1)), (BondId(1), BondId(0))], 2, 2)
                .unwrap(),
            Correspondence::new(
                vec![
                    (DativeBondId(0), DativeBondId(1)),
                    (DativeBondId(1), DativeBondId(0)),
                ],
                2,
                2,
            )
            .unwrap(),
            Correspondence::new(vec![(AromaticSystemId(0), AromaticSystemId(0))], 1, 1).unwrap(),
            Correspondence::new(
                vec![
                    (MulticenterBondId(0), MulticenterBondId(1)),
                    (MulticenterBondId(1), MulticenterBondId(0)),
                ],
                2,
                2,
            )
            .unwrap(),
            Correspondence::new(vec![(NoncovalentBondId(0), NoncovalentBondId(0))], 1, 1).unwrap(),
            Correspondence::new(
                vec![
                    (StereoAtomId(0), StereoAtomId(1)),
                    (StereoAtomId(1), StereoAtomId(0)),
                ],
                2,
                2,
            )
            .unwrap(),
            Correspondence::new(vec![(StereoBondId(0), StereoBondId(0))], 1, 1).unwrap(),
        );
        assert_eq!(MoleculeCorrespondence::from(&molecule_remapping), expected);
    }

    #[rstest]
    fn test_molecule_correspondence_from_remapping_empty() {
        assert_eq!(
            MoleculeCorrespondence::from(&MoleculeRemapping::default()),
            MoleculeCorrespondence::empty()
        );
    }

    #[rstest]
    fn test_molecule_remapping_accessors(molecule_remapping: MoleculeRemapping) {
        assert_eq!(
            molecule_remapping.graph(),
            &GraphRemapping::new(
                Remapping::new(vec![NodeId(1), NodeId(0)]).unwrap(),
                Remapping::new(vec![EdgeId(1), EdgeId(0)]).unwrap()
            )
        );
        assert_eq!(
            molecule_remapping.dative_bonds(),
            &Remapping::new(vec![DativeBondId(1), DativeBondId(0)]).unwrap()
        );
        assert_eq!(
            molecule_remapping.aromatic_systems(),
            &Remapping::new(vec![AromaticSystemId(0)]).unwrap()
        );
        assert_eq!(
            molecule_remapping.multicenter_bonds(),
            &Remapping::new(vec![MulticenterBondId(1), MulticenterBondId(0)]).unwrap()
        );
        assert_eq!(
            molecule_remapping.noncovalent_bonds(),
            &Remapping::new(vec![NoncovalentBondId(0)]).unwrap()
        );
        assert_eq!(
            molecule_remapping.stereo_atoms(),
            &Remapping::new(vec![StereoAtomId(1), StereoAtomId(0)]).unwrap()
        );
        assert_eq!(
            molecule_remapping.stereo_bonds(),
            &Remapping::new(vec![StereoBondId(0)]).unwrap()
        );
    }

    #[rstest]
    #[case::first(AtomId(0), Some(AtomId(1)))]
    #[case::last(AtomId(1), Some(AtomId(0)))]
    #[case::uncovered(AtomId(2), None)]
    fn test_molecule_remapping_try_map_atom(
        molecule_remapping: MoleculeRemapping,
        #[case] id: AtomId,
        #[case] expected: Option<AtomId>,
    ) {
        assert_eq!(molecule_remapping.try_map_atom(id), expected);
    }

    #[rstest]
    #[case::first(BondId(0), Some(BondId(1)))]
    #[case::last(BondId(1), Some(BondId(0)))]
    #[case::uncovered(BondId(2), None)]
    fn test_molecule_remapping_try_map_bond(
        molecule_remapping: MoleculeRemapping,
        #[case] id: BondId,
        #[case] expected: Option<BondId>,
    ) {
        assert_eq!(molecule_remapping.try_map_bond(id), expected);
    }

    #[rstest]
    #[case::first(DativeBondId(0), Some(DativeBondId(1)))]
    #[case::last(DativeBondId(1), Some(DativeBondId(0)))]
    #[case::uncovered(DativeBondId(2), None)]
    fn test_molecule_remapping_try_map_dative_bond(
        molecule_remapping: MoleculeRemapping,
        #[case] id: DativeBondId,
        #[case] expected: Option<DativeBondId>,
    ) {
        assert_eq!(molecule_remapping.try_map_dative_bond(id), expected);
    }

    #[rstest]
    #[case::identity(AromaticSystemId(0), Some(AromaticSystemId(0)))]
    #[case::uncovered(AromaticSystemId(1), None)]
    fn test_molecule_remapping_try_map_aromatic_system(
        molecule_remapping: MoleculeRemapping,
        #[case] id: AromaticSystemId,
        #[case] expected: Option<AromaticSystemId>,
    ) {
        assert_eq!(molecule_remapping.try_map_aromatic_system(id), expected);
    }

    #[rstest]
    #[case::first(MulticenterBondId(0), Some(MulticenterBondId(1)))]
    #[case::last(MulticenterBondId(1), Some(MulticenterBondId(0)))]
    #[case::uncovered(MulticenterBondId(2), None)]
    fn test_molecule_remapping_try_map_multicenter_bond(
        molecule_remapping: MoleculeRemapping,
        #[case] id: MulticenterBondId,
        #[case] expected: Option<MulticenterBondId>,
    ) {
        assert_eq!(molecule_remapping.try_map_multicenter_bond(id), expected);
    }

    #[rstest]
    #[case::first(NoncovalentBondId(0), Some(NoncovalentBondId(0)))]
    #[case::uncovered(NoncovalentBondId(1), None)]
    fn test_molecule_remapping_try_map_noncovalent_bond(
        molecule_remapping: MoleculeRemapping,
        #[case] id: NoncovalentBondId,
        #[case] expected: Option<NoncovalentBondId>,
    ) {
        assert_eq!(molecule_remapping.try_map_noncovalent_bond(id), expected);
    }

    #[rstest]
    #[case::first(StereoAtomId(0), Some(StereoAtomId(1)))]
    #[case::last(StereoAtomId(1), Some(StereoAtomId(0)))]
    #[case::uncovered(StereoAtomId(2), None)]
    fn test_molecule_remapping_try_map_stereo_atom(
        molecule_remapping: MoleculeRemapping,
        #[case] id: StereoAtomId,
        #[case] expected: Option<StereoAtomId>,
    ) {
        assert_eq!(molecule_remapping.try_map_stereo_atom(id), expected);
    }

    #[rstest]
    #[case::first(StereoBondId(0), Some(StereoBondId(0)))]
    #[case::uncovered(StereoBondId(1), None)]
    fn test_molecule_remapping_try_map_stereo_bond(
        molecule_remapping: MoleculeRemapping,
        #[case] id: StereoBondId,
        #[case] expected: Option<StereoBondId>,
    ) {
        assert_eq!(molecule_remapping.try_map_stereo_bond(id), expected);
    }

    #[rstest]
    fn test_molecule_remapping_map(molecule_remapping: MoleculeRemapping) {
        assert_eq!(molecule_remapping.map_atom(AtomId(0)), AtomId(1));
        assert_eq!(molecule_remapping.map_bond(BondId(1)), BondId(0));
        assert_eq!(
            molecule_remapping.map_dative_bond(DativeBondId(0)),
            DativeBondId(1)
        );
        assert_eq!(
            molecule_remapping.map_aromatic_system(AromaticSystemId(0)),
            AromaticSystemId(0)
        );
        assert_eq!(
            molecule_remapping.map_multicenter_bond(MulticenterBondId(1)),
            MulticenterBondId(0)
        );
        assert_eq!(
            molecule_remapping.map_noncovalent_bond(NoncovalentBondId(0)),
            NoncovalentBondId(0)
        );
        assert_eq!(
            molecule_remapping.map_stereo_atom(StereoAtomId(1)),
            StereoAtomId(0)
        );
        assert_eq!(
            molecule_remapping.map_stereo_bond(StereoBondId(0)),
            StereoBondId(0)
        );
    }

    #[rstest]
    #[should_panic(expected = "stereo-bond id outside remapping source domain")]
    fn test_molecule_remapping_map_error(molecule_remapping: MoleculeRemapping) {
        molecule_remapping.map_stereo_bond(StereoBondId(1));
    }

    #[rstest]
    fn test_molecule_remapping_default() {
        assert_eq!(MoleculeRemapping::default(), MoleculeRemapping::empty());
    }
}
