//! Graph-IR id remappings between `Molecule` id spaces.
//!
//! [`MoleculeRemapping`] is the dense total relabeling over all eight molecule entity id spaces.
//! [`IdRemapping`] supplies sparse lookup tables for current reference-transport operations.

use std::collections::HashMap;

use umol_graph_core::{EdgeId, GraphRemapping, NodeId, Remapping};

use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};

/// Total relabeling of all eight `Molecule` entity id spaces.
///
/// Each component's image vector defines its dense source domain. Every id in that domain has an
/// image, while images may repeat or occupy only part of a larger target id space. Consumers that
/// require injective, surjective, or dense-target mappings establish those contextual properties.
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
    /// Construct a remapping from the complete image vector for each source id space.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        graph: GraphRemapping,
        dative_bonds: Vec<DativeBondId>,
        aromatic_systems: Vec<AromaticSystemId>,
        multicenter_bonds: Vec<MulticenterBondId>,
        noncovalent_bonds: Vec<NoncovalentBondId>,
        stereo_atoms: Vec<StereoAtomId>,
        stereo_bonds: Vec<StereoBondId>,
    ) -> Self {
        Self {
            graph,
            dative_bonds: Remapping::new(dative_bonds),
            aromatic_systems: Remapping::new(aromatic_systems),
            multicenter_bonds: Remapping::new(multicenter_bonds),
            noncovalent_bonds: Remapping::new(noncovalent_bonds),
            stereo_atoms: Remapping::new(stereo_atoms),
            stereo_bonds: Remapping::new(stereo_bonds),
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
        Self::new(
            GraphRemapping::new(Vec::new(), Vec::new()),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
    }
}

/// Sparse lookup tables used to transport references between `Molecule` id spaces.
///
/// Every id read by a consuming operation must be present. Unlike [`MoleculeRemapping`], this type
/// does not declare a dense source domain.
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
    use rstest::*;

    use super::*;

    #[fixture]
    fn molecule_remapping() -> MoleculeRemapping {
        MoleculeRemapping::new(
            GraphRemapping::new(vec![NodeId(4), NodeId(1)], vec![EdgeId(5), EdgeId(5)]),
            vec![DativeBondId(2), DativeBondId(0)],
            vec![AromaticSystemId(8)],
            vec![MulticenterBondId(1), MulticenterBondId(4)],
            vec![NoncovalentBondId(0)],
            vec![StereoAtomId(3), StereoAtomId(3)],
            vec![StereoBondId(6)],
        )
    }

    #[rstest]
    fn test_molecule_remapping_accessors(molecule_remapping: MoleculeRemapping) {
        assert_eq!(
            molecule_remapping.graph(),
            &GraphRemapping::new(vec![NodeId(4), NodeId(1)], vec![EdgeId(5), EdgeId(5)])
        );
        assert_eq!(
            molecule_remapping.dative_bonds(),
            &Remapping::new(vec![DativeBondId(2), DativeBondId(0)])
        );
        assert_eq!(
            molecule_remapping.aromatic_systems(),
            &Remapping::new(vec![AromaticSystemId(8)])
        );
        assert_eq!(
            molecule_remapping.multicenter_bonds(),
            &Remapping::new(vec![MulticenterBondId(1), MulticenterBondId(4)])
        );
        assert_eq!(
            molecule_remapping.noncovalent_bonds(),
            &Remapping::new(vec![NoncovalentBondId(0)])
        );
        assert_eq!(
            molecule_remapping.stereo_atoms(),
            &Remapping::new(vec![StereoAtomId(3), StereoAtomId(3)])
        );
        assert_eq!(
            molecule_remapping.stereo_bonds(),
            &Remapping::new(vec![StereoBondId(6)])
        );
    }

    #[rstest]
    #[case::first(AtomId(0), Some(AtomId(4)))]
    #[case::last(AtomId(1), Some(AtomId(1)))]
    #[case::uncovered(AtomId(2), None)]
    fn test_molecule_remapping_try_map_atom(
        molecule_remapping: MoleculeRemapping,
        #[case] id: AtomId,
        #[case] expected: Option<AtomId>,
    ) {
        assert_eq!(molecule_remapping.try_map_atom(id), expected);
    }

    #[rstest]
    #[case::first(BondId(0), Some(BondId(5)))]
    #[case::repeated_target(BondId(1), Some(BondId(5)))]
    #[case::uncovered(BondId(2), None)]
    fn test_molecule_remapping_try_map_bond(
        molecule_remapping: MoleculeRemapping,
        #[case] id: BondId,
        #[case] expected: Option<BondId>,
    ) {
        assert_eq!(molecule_remapping.try_map_bond(id), expected);
    }

    #[rstest]
    #[case::first(DativeBondId(0), Some(DativeBondId(2)))]
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
    #[case::sparse_target(AromaticSystemId(0), Some(AromaticSystemId(8)))]
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
    #[case::last(MulticenterBondId(1), Some(MulticenterBondId(4)))]
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
    #[case::first(StereoAtomId(0), Some(StereoAtomId(3)))]
    #[case::repeated_target(StereoAtomId(1), Some(StereoAtomId(3)))]
    #[case::uncovered(StereoAtomId(2), None)]
    fn test_molecule_remapping_try_map_stereo_atom(
        molecule_remapping: MoleculeRemapping,
        #[case] id: StereoAtomId,
        #[case] expected: Option<StereoAtomId>,
    ) {
        assert_eq!(molecule_remapping.try_map_stereo_atom(id), expected);
    }

    #[rstest]
    #[case::first(StereoBondId(0), Some(StereoBondId(6)))]
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
        assert_eq!(molecule_remapping.map_atom(AtomId(0)), AtomId(4));
        assert_eq!(molecule_remapping.map_bond(BondId(1)), BondId(5));
        assert_eq!(
            molecule_remapping.map_dative_bond(DativeBondId(0)),
            DativeBondId(2)
        );
        assert_eq!(
            molecule_remapping.map_aromatic_system(AromaticSystemId(0)),
            AromaticSystemId(8)
        );
        assert_eq!(
            molecule_remapping.map_multicenter_bond(MulticenterBondId(1)),
            MulticenterBondId(4)
        );
        assert_eq!(
            molecule_remapping.map_noncovalent_bond(NoncovalentBondId(0)),
            NoncovalentBondId(0)
        );
        assert_eq!(
            molecule_remapping.map_stereo_atom(StereoAtomId(1)),
            StereoAtomId(3)
        );
        assert_eq!(
            molecule_remapping.map_stereo_bond(StereoBondId(0)),
            StereoBondId(6)
        );
    }

    #[rstest]
    #[should_panic(expected = "stereo-bond id outside remapping source domain")]
    fn test_molecule_remapping_map_error(molecule_remapping: MoleculeRemapping) {
        molecule_remapping.map_stereo_bond(StereoBondId(1));
    }

    #[rstest]
    fn test_molecule_remapping_default() {
        assert_eq!(
            MoleculeRemapping::default(),
            MoleculeRemapping::new(
                GraphRemapping::new(vec![], vec![]),
                vec![],
                vec![],
                vec![],
                vec![],
                vec![],
                vec![],
            )
        );
    }
}
