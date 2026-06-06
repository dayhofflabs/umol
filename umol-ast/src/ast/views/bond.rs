//! Bond views: `BondViews` namespace, `BondView` / `BondViewMut` AST bundles,
//! `BondBuilderView` / `BondBuilderViewMut` builder bundles.

use std::ops::Index;

use umol_graph_core::{EdgeId, NodeId};

use super::super::bond::BondAst;
use super::super::constraint::BondConstraints;
use super::super::ids::{AtomId, BondId, StereoBondId};
use super::super::molecule::MoleculeAst;
use super::super::rings::{RingSet, RingView};
use super::super::spin::SpinStateAst;
use super::super::stereo::StereoKind;
use super::super::traits::Lattice;
use super::super::value::ValueAst;
use super::aromatic::AromaticSystemView;
use super::atom::AtomView;
use super::stereo::StereoBondView;

/// Namespace accessor for bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct BondViews<'a> {
    molecule: &'a MoleculeAst,
    bonds: &'a [BondAst],
}

impl<'a> BondViews<'a> {
    pub(crate) fn new(molecule: &'a MoleculeAst, bonds: &'a [BondAst]) -> Self {
        Self { molecule, bonds }
    }

    pub fn count(&self) -> usize {
        self.molecule.raw_graph().edge_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = BondId> {
        self.molecule.raw_graph().edge_ids().map(BondId::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = BondView<'a>> {
        let molecule = self.molecule;
        let bonds = self.bonds;
        let graph = molecule.raw_graph();
        graph.edge_ids().map(move |id| {
            let [s, t] = graph.edge_endpoints(id);
            BondView {
                id: BondId::from(id),
                atoms: [s, t],
                ast: &bonds[id.index()],
                molecule,
            }
        })
    }

    pub fn contains(&self, id: BondId) -> bool {
        self.molecule.raw_graph().contains_edge(EdgeId::from(id))
    }

    pub fn get(&self, id: BondId) -> Option<BondView<'a>> {
        if !self.contains(id) {
            return None;
        }
        let [s, t] = self.molecule.raw_graph().edge_endpoints(EdgeId::from(id));
        Some(BondView {
            id,
            atoms: [s, t],
            ast: &self.bonds[id.index()],
            molecule: self.molecule,
        })
    }

    /// Id of the bond between `a` and `b`, if any.
    pub fn connecting_id(&self, a: AtomId, b: AtomId) -> Option<BondId> {
        self.molecule
            .raw_graph()
            .find_edge(NodeId::from(a), NodeId::from(b))
            .map(BondId::from)
    }

    /// View of the bond between `a` and `b`, if any.
    pub fn connecting(&self, a: AtomId, b: AtomId) -> Option<BondView<'a>> {
        self.connecting_id(a, b).map(|id| {
            self.get(id)
                .expect("bond id from graph must refer to a bond in this molecule")
        })
    }

    /// Ids of bonds whose both endpoints lie in `atoms`.
    pub fn induced_ids(&self, atoms: &[AtomId]) -> Vec<BondId> {
        let mut nodes: Vec<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        nodes.sort_unstable();
        self.molecule
            .raw_graph()
            .induced_edges(&nodes)
            .map(BondId::from)
            .collect()
    }

    /// Views of bonds whose both endpoints lie in `atoms`.
    pub fn induced(&self, atoms: &[AtomId]) -> Vec<BondView<'a>> {
        self.induced_ids(atoms)
            .into_iter()
            .map(|id| {
                self.get(id)
                    .expect("bond id from graph must refer to a bond in this molecule")
            })
            .collect()
    }
}

impl<'a> Index<BondId> for BondViews<'a> {
    type Output = BondAst;
    fn index(&self, id: BondId) -> &BondAst {
        &self.bonds[id.index()]
    }
}

/// Borrowed view of a bond: its index, the two participating atoms, and data.
#[derive(Clone, Copy, Debug)]
pub struct BondView<'a> {
    pub id: BondId,
    atoms: [NodeId; 2],
    pub ast: &'a BondAst,
    molecule: &'a MoleculeAst,
}

impl<'a> BondView<'a> {
    #[inline]
    pub fn order(&self) -> &'a ValueAst {
        &self.ast.order
    }

    #[inline]
    pub fn charge(&self) -> &'a ValueAst {
        &self.ast.charge
    }

    #[inline]
    pub fn spin(&self) -> &'a SpinStateAst {
        &self.ast.spin
    }

    #[inline]
    pub fn constraints(&self) -> &'a BondConstraints {
        &self.ast.constraints
    }

    /// The two atom indices incident to this bond.
    pub fn atom_ids(&self) -> [AtomId; 2] {
        self.atoms.map(AtomId::from)
    }

    /// Views of the two atoms incident to this bond.
    pub fn atoms(&self) -> impl Iterator<Item = AtomView<'a>> + 'a {
        let molecule = self.molecule;
        self.atoms
            .into_iter()
            .map(move |id| molecule.atom(AtomId::from(id)))
    }

    /// The aromatic system this bond participates in, if any. A bond is in
    /// an aromatic system iff both endpoints belong to that system.
    pub fn aromatic_system(&self) -> Option<AromaticSystemView<'a>> {
        let [a, b] = self.atom_ids();
        self.molecule
            .aromatic_systems()
            .incident(a)
            .find(|sys| sys.atom_ids().any(|x| x == b))
    }

    pub fn is_in_aromatic_system(&self) -> bool {
        self.aromatic_system().is_some()
    }

    pub fn is_in_cis_trans_stereo(&self) -> bool {
        self.cis_trans_stereo().is_some()
    }

    pub fn cis_trans_stereo_id(&self) -> Option<StereoBondId> {
        self.cis_trans_stereo().map(|s| s.id)
    }

    /// The cis/trans stereo bond sited on this bond, if any. A bond is the
    /// site of at most one stereo bond; the kind filter selects the cis/trans
    /// case from any other bond-centered geometries that share the relation.
    pub fn cis_trans_stereo(&self) -> Option<StereoBondView<'a>> {
        self.molecule
            .stereo_bonds()
            .coincident(self.id)
            .filter(|s| s.kind() == StereoKind::CisTrans)
    }

    /// True if this bond belongs to any ring in the molecule's canonical
    /// ring set (Vismara relevant cycles, max ring size 22). Uses the
    /// molecule's cached canonical `RingSet`.
    pub fn is_in_ring(&self) -> bool {
        self.molecule.rings().contains_bond(self.id)
    }

    /// True if this bond appears in any ring of the supplied set.
    pub fn is_in_ring_from(&self, rings: &RingSet) -> bool {
        rings.contains_bond(self.id)
    }

    /// Rings containing this bond drawn from the molecule's canonical
    /// `RingSet` (Vismara relevant cycles, max ring size 22).
    pub fn rings(&self) -> impl Iterator<Item = RingView<'a>> + 'a {
        let id = self.id;
        self.molecule
            .rings()
            .iter()
            .filter(move |v| v.bonds().contains(&id))
    }

    /// Rings from the supplied set that contain this bond.
    pub fn rings_from<'r>(&self, rings: &'r RingSet) -> impl Iterator<Item = RingView<'r>> + 'r {
        let id = self.id;
        rings.iter().filter(move |v| v.bonds().contains(&id))
    }

    /// Count of canonical rings (Vismara / max ring size 22) containing
    /// this bond. Always `Lit`.
    pub fn ring_count(&self) -> ValueAst {
        ValueAst::Lit(self.rings().count() as i64)
    }

    /// Sizes of canonical rings containing this bond, in iteration order.
    /// Multi-valued: a bond shared between fused rings yields one size per
    /// ring.
    pub fn ring_size(&self) -> impl Iterator<Item = usize> + 'a {
        self.rings().map(|r| r.len())
    }

    /// Is bond ground
    pub fn is_ground(&self) -> bool {
        self.ast.is_ground()
    }

    /// Is bond undetermined
    pub fn is_undetermined(&self) -> bool {
        self.ast.is_undetermined()
    }
}

/// Mutable borrowed view of a bond.
#[derive(Debug)]
pub struct BondViewMut<'a> {
    pub id: BondId,
    atoms: [AtomId; 2],
    pub ast: &'a mut BondAst,
}

impl<'a> BondViewMut<'a> {
    pub(crate) fn new(id: BondId, atoms: [AtomId; 2], ast: &'a mut BondAst) -> Self {
        Self { id, atoms, ast }
    }

    /// The two atoms incident to this bond.
    pub fn atoms(&self) -> [AtomId; 2] {
        self.atoms
    }
}

// Builder-scope view bundles for bonds.

pub struct BondBuilderView<'a> {
    pub id: BondId,
    pub ast: &'a BondAst,
    pub atoms: [AtomId; 2],
}

pub struct BondBuilderViewMut<'a> {
    pub id: BondId,
    pub ast: &'a mut BondAst,
    pub atoms: [AtomId; 2],
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_shared::element::Element;

    use crate::ast::aromatic::AromaticSystemAst;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::constraint::Constraints;
    use crate::ast::dative::DativeBondAst;
    use crate::ast::ids::{AromaticSystemId, AtomId, BondId, StereoBondId};
    use crate::ast::ligand::{StereoLigand, StereoLigandKind};
    use crate::ast::molecule::MoleculeAst;
    use crate::ast::multicenter::MulticenterBondAst;
    use crate::ast::noncovalent::{NoncovalentBondAst, NoncovalentBondKind};
    use crate::ast::rings::RingFamily;
    use crate::ast::stereo::{StereoBondAst, StereoCosetAst, StereoKind};
    use crate::ast::value::ValueAst;

    #[fixture]
    fn molecule() -> MoleculeAst {
        MoleculeAst::from_parts(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::O),
            ],
            vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(2)),
                (AtomId(2), AtomId(3), BondAst::from_order(1)),
            ],
            vec![(vec![AtomId(2)], AtomId(3), DativeBondAst::from_order(1))],
            vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemAst::default(),
            )],
            vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                MulticenterBondAst::default(),
            )],
            vec![(
                AtomId(0),
                AtomId(3),
                NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            Vec::new(),
            Vec::new(),
            Constraints::default(),
        )
    }

    #[fixture]
    fn ring_with_chain() -> MoleculeAst {
        MoleculeAst::from_atoms_and_bonds(
            vec![AtomAst::from_element(Element::C); 7],
            vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(1)),
                (AtomId(2), AtomId(3), BondAst::from_order(1)),
                (AtomId(3), AtomId(4), BondAst::from_order(1)),
                (AtomId(4), AtomId(5), BondAst::from_order(1)),
                (AtomId(5), AtomId(0), BondAst::from_order(1)),
                (AtomId(0), AtomId(6), BondAst::from_order(1)),
            ],
        )
    }

    #[rstest]
    fn test_bond_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.bonds().count(), 3);
    }

    #[rstest]
    fn test_bond_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.bonds().ids().collect::<Vec<_>>(),
            vec![BondId(0), BondId(1), BondId(2)],
        );
    }

    #[rstest]
    fn test_bond_views_iter(molecule: MoleculeAst) {
        let collected: Vec<(BondId, [AtomId; 2], BondAst)> = molecule
            .bonds()
            .iter()
            .map(|v| (v.id, v.atom_ids(), v.ast.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![
                (BondId(0), [AtomId(0), AtomId(1)], BondAst::from_order(1)),
                (BondId(1), [AtomId(1), AtomId(2)], BondAst::from_order(2)),
                (BondId(2), [AtomId(2), AtomId(3)], BondAst::from_order(1)),
            ],
        );
    }

    #[rstest]
    #[case::present(BondId(1), true)]
    #[case::absent(BondId(99), false)]
    fn test_bond_views_contains(molecule: MoleculeAst, #[case] id: BondId, #[case] expected: bool) {
        assert_eq!(molecule.bonds().contains(id), expected);
    }

    #[rstest]
    fn test_bond_views_get(molecule: MoleculeAst) {
        let res = molecule.bonds().get(BondId(1));
        assert!(res.is_some());
        let view = res.unwrap();
        assert_eq!(view.id, BondId(1));
        assert_eq!(view.atom_ids(), [AtomId(1), AtomId(2)]);
        assert_eq!(*view.ast, BondAst::from_order(2));
    }

    #[rstest]
    fn test_bond_views_get_none(molecule: MoleculeAst) {
        let res = molecule.bonds().get(BondId(99));
        assert!(res.is_none());
    }

    #[rstest]
    fn test_bond_views_index(molecule: MoleculeAst) {
        let bond: &BondAst = &molecule.bonds()[BondId(1)];
        assert_eq!(*bond, BondAst::from_order(2));
    }

    #[rstest]
    fn test_bond_view_atom_ids(molecule: MoleculeAst) {
        assert_eq!(molecule.bond(BondId(1)).atom_ids(), [AtomId(1), AtomId(2)]);
    }

    #[rstest]
    fn test_bond_view_atoms(molecule: MoleculeAst) {
        let ids: Vec<AtomId> = molecule.bond(BondId(1)).atoms().map(|a| a.id).collect();
        assert_eq!(ids, vec![AtomId(1), AtomId(2)]);
    }

    #[rstest]
    #[case::both_endpoints_aromatic(BondId(0), Some(AromaticSystemId(0)))]
    #[case::both_endpoints_aromatic_alt(BondId(1), Some(AromaticSystemId(0)))]
    #[case::one_endpoint_outside(BondId(2), None)]
    fn test_bond_view_aromatic_system(
        molecule: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected: Option<AromaticSystemId>,
    ) {
        let id = molecule.bond(bond).aromatic_system().map(|v| v.id);
        assert_eq!(id, expected);
    }

    #[rstest]
    #[case::both_endpoints_aromatic(BondId(0), true)]
    #[case::both_endpoints_aromatic_alt(BondId(1), true)]
    #[case::one_endpoint_outside(BondId(2), false)]
    fn test_bond_view_is_in_aromatic_system(
        molecule: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected: bool,
    ) {
        assert_eq!(molecule.bond(bond).is_in_aromatic_system(), expected);
    }

    #[fixture]
    fn stereo_molecule() -> MoleculeAst {
        MoleculeAst::from_parts(
            vec![AtomAst::from_element(Element::C); 4],
            vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(2)),
                (AtomId(2), AtomId(3), BondAst::from_order(1)),
            ],
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            vec![(
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                ],
                StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Lit(1)),
            )],
            Constraints::default(),
        )
    }

    #[rstest]
    #[case::site(BondId(1), true)]
    #[case::non_site(BondId(0), false)]
    fn test_bond_view_is_in_cis_trans_stereo(
        stereo_molecule: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected: bool,
    ) {
        assert_eq!(
            stereo_molecule.bond(bond).is_in_cis_trans_stereo(),
            expected
        );
    }

    #[rstest]
    #[case::site(BondId(1), Some(StereoBondId(0)))]
    #[case::non_site(BondId(0), None)]
    fn test_bond_view_cis_trans_stereo_id(
        stereo_molecule: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected: Option<StereoBondId>,
    ) {
        assert_eq!(stereo_molecule.bond(bond).cis_trans_stereo_id(), expected);
    }

    #[rstest]
    fn test_bond_view_cis_trans_stereo(stereo_molecule: MoleculeAst) {
        let view = stereo_molecule.bond(BondId(1)).cis_trans_stereo().unwrap();
        assert_eq!(view.id, StereoBondId(0));
        assert_eq!(view.kind(), StereoKind::CisTrans);
        assert!(stereo_molecule.bond(BondId(0)).cis_trans_stereo().is_none());
    }

    #[rstest]
    #[case::ring_bond_0_1(BondId(0), true)]
    #[case::ring_bond_5_0(BondId(5), true)]
    #[case::chain_bond_0_6(BondId(6), false)]
    fn test_bond_view_is_in_ring(
        ring_with_chain: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected: bool,
    ) {
        assert_eq!(ring_with_chain.bond(bond).is_in_ring(), expected);
    }

    #[rstest]
    #[case::ring_bond(BondId(0), true)]
    #[case::chain_bond(BondId(6), false)]
    fn test_bond_view_is_in_ring_from(
        ring_with_chain: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected: bool,
    ) {
        let rings = ring_with_chain.rings_with(RingFamily::Relevant, 22, |_| true);
        assert_eq!(ring_with_chain.bond(bond).is_in_ring_from(&rings), expected);
    }

    #[rstest]
    #[case::ring_bond(BondId(0), 1)]
    #[case::chain_bond(BondId(6), 0)]
    fn test_bond_view_rings_from(
        ring_with_chain: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected_count: usize,
    ) {
        let rings = ring_with_chain.rings_with(RingFamily::Relevant, 22, |_| true);
        let count = ring_with_chain.bond(bond).rings_from(&rings).count();
        assert_eq!(count, expected_count);
    }

    #[rstest]
    #[case::ring_bond(BondId(0), ValueAst::Lit(1))]
    #[case::chain_bond(BondId(6), ValueAst::Lit(0))]
    fn test_bond_view_ring_count(
        ring_with_chain: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(ring_with_chain.bond(bond).ring_count(), expected);
    }

    #[rstest]
    #[case::ring_bond(BondId(0), vec![6])]
    #[case::chain_bond(BondId(6), vec![])]
    fn test_bond_view_ring_size(
        ring_with_chain: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected: Vec<usize>,
    ) {
        let sizes: Vec<_> = ring_with_chain.bond(bond).ring_size().collect();
        assert_eq!(sizes, expected);
    }
}
