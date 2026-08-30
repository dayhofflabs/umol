//! Bond views: `BondViews` namespace, `BondView` / `BondViewMut` attribute bundles,
//! `BondEditorView` / `BondEditorViewMut` builder bundles.

use umol_graph_core::{EdgeId, NodeId};

use super::super::bond::BondForm;
use super::super::boolean::BooleanForm;
use super::super::constraint::{BondConstraintForm, BondConstraintKey, BondConstraintsForm};
use super::super::id::{AtomId, BondId, StereoBondId};
use super::super::molecule::Molecule;
use super::super::num::NumForm;
use super::super::ring::RingSet;
use super::super::spin::UnpairedElectronsForm;
use super::super::stereo::{CisTransStereoForm, StereoKind};
use super::super::traits::Lattice;
use super::aromatic::AromaticSystemView;
use super::atom::AtomView;
use super::constraints::BondConstraintsView;
use super::ring::bond_ring_membership;
use super::stereo::StereoBondView;

/// Namespace accessor for bond views on a `Molecule`.
#[derive(Clone, Copy)]
pub struct BondViews<'a> {
    molecule: &'a Molecule,
    bonds: &'a [BondForm],
}

impl<'a> BondViews<'a> {
    pub(crate) fn new(molecule: &'a Molecule, bonds: &'a [BondForm]) -> Self {
        Self { molecule, bonds }
    }

    pub fn count(&self) -> usize {
        self.molecule.raw_graph().edge_count()
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = BondId> {
        self.molecule.raw_graph().edge_ids().map(BondId::from)
    }

    pub fn iter(&self) -> impl ExactSizeIterator<Item = BondView<'a>> {
        let molecule = self.molecule;
        let bonds = self.bonds;
        let graph = molecule.raw_graph();
        graph.edge_ids().map(move |id| {
            let [s, t] = graph.edge_endpoints(id);
            BondView {
                id: BondId::from(id),
                atoms: [s, t],
                attributes: &bonds[id.index()],
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
            attributes: &self.bonds[id.index()],
            molecule: self.molecule,
        })
    }

    /// Id of the bond between `first` and `second`, if any.
    pub fn of_id(&self, first: AtomId, second: AtomId) -> Option<BondId> {
        self.molecule
            .raw_graph()
            .find_edge(NodeId::from(first), NodeId::from(second))
            .map(BondId::from)
    }

    /// View of the bond between `first` and `second`, if any.
    pub fn of(&self, first: AtomId, second: AtomId) -> Option<BondView<'a>> {
        self.of_id(first, second).map(|id| {
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

/// Borrowed view of a bond: its index, the two participating atoms, and data.
#[derive(Clone, Copy, Debug)]
pub struct BondView<'a> {
    pub id: BondId,
    atoms: [NodeId; 2],
    pub attributes: &'a BondForm,
    molecule: &'a Molecule,
}

impl<'a> BondView<'a> {
    #[inline]
    pub fn order(&self) -> &'a NumForm {
        &self.attributes.order
    }

    #[inline]
    pub fn charge(&self) -> &'a NumForm {
        &self.attributes.charge
    }

    #[inline]
    pub fn unpaired_electrons(&self) -> &'a UnpairedElectronsForm {
        &self.attributes.unpaired_electrons
    }

    /// Constraint reading of this bond: the container's read API (asserted
    /// side, meanings intact) plus the keyed `asserted`/`derived`/
    /// `derived_complete` accessors. Mutation stays on the stored container.
    #[inline]
    pub fn constraints(&self) -> BondConstraintsView<'a> {
        BondConstraintsView::new(self.molecule, self.id)
    }

    /// The two atom indices incident to this bond.
    pub fn atom_ids(&self) -> [AtomId; 2] {
        self.atoms.map(AtomId::from)
    }

    /// Views of the two atoms incident to this bond.
    pub fn atoms(&self) -> impl ExactSizeIterator<Item = AtomView<'a>> + 'a {
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

    pub fn is_stereo_bond(&self) -> bool {
        self.stereo_bond().is_some()
    }

    pub fn stereo_bond_id(&self) -> Option<StereoBondId> {
        self.stereo_bond().map(|s| s.id)
    }

    /// The stereo bond sited on this bond, if any — any bond-centered geometry. A
    /// bond is the site of at most one stereo bond.
    pub fn stereo_bond(&self) -> Option<StereoBondView<'a>> {
        self.molecule.stereo_bonds().at(self.id)
    }

    /// Is bond ground
    pub fn is_ground(&self) -> bool {
        self.attributes.is_ground()
    }

    /// Is bond undetermined
    pub fn is_undetermined(&self) -> bool {
        self.attributes.is_undetermined()
    }
}

// Derivation layer beneath the bond facades: functions of the molecule and
// bond id, presented by `BondView` (typed) and `BondConstraintsView` (keyed).

/// Stored constraint container of `bond`.
pub(crate) fn bond_asserted_constraints(molecule: &Molecule, bond: BondId) -> &BondConstraintsForm {
    &molecule.bond(bond).attributes.constraints
}

/// Asserted side of one bond constraint key under resolution's closed-world
/// claim: the stored assertion, else the absence cell closed to its definite
/// negative. Never reads relations — a bond inside a stored aromatic system
/// without its own `#a` assertion reads `Lit(false)`.
pub(crate) fn bond_asserted_complete_constraint(
    molecule: &Molecule,
    bond: BondId,
    key: BondConstraintKey,
) -> Option<BondConstraintForm> {
    if let Some(asserted) = bond_asserted_constraints(molecule, bond).get(key) {
        return Some(asserted.clone());
    }
    match key {
        BondConstraintKey::Aromatic => Some(BondConstraintForm::aromatic(BooleanForm::Lit(false))),
        BondConstraintKey::CisTransStereo => Some(BondConstraintForm::cis_trans_stereo(
            CisTransStereoForm::NotStereo,
        )),
        BondConstraintKey::RingMembership(_) => None,
    }
}

/// Derived side of one bond constraint key, read from the molecule's
/// relations.
///
/// `complete` selects the closure reading: absence of a resolution-written
/// overlay yields its definite negative (`Aromatic(false)` / `NotStereo`)
/// instead of no value. A bond is in an aromatic system iff both endpoints
/// belong to that system. Ring keys require `rings` and panic without it —
/// the caller scanning keys decides whether to build the ring set.
pub(crate) fn bond_derived_constraint(
    molecule: &Molecule,
    bond: BondId,
    rings: Option<&RingSet>,
    key: BondConstraintKey,
    complete: bool,
) -> Option<BondConstraintForm> {
    match key {
        BondConstraintKey::Aromatic => {
            if molecule.bond(bond).is_in_aromatic_system() {
                Some(BondConstraintForm::aromatic(BooleanForm::Lit(true)))
            } else if complete {
                Some(BondConstraintForm::aromatic(BooleanForm::Lit(false)))
            } else {
                None
            }
        }
        BondConstraintKey::CisTransStereo => {
            // Total over staging entities: an undetermined kind or coset is
            // an open claim, not the absence of one.
            if let Some(stereo) = molecule.stereo_bonds().at(bond) {
                let form = match (
                    stereo.attributes.configuration.kind(),
                    stereo.attributes.configuration.coset(),
                ) {
                    (Some(StereoKind::CisTrans), Some(coset)) => {
                        CisTransStereoForm::stereo(coset.clone())
                    }
                    (Some(StereoKind::CisTrans), None) | (None, _) => {
                        CisTransStereoForm::Undetermined
                    }
                    (Some(_), _) => CisTransStereoForm::NotStereo,
                };
                Some(BondConstraintForm::cis_trans_stereo(form))
            } else if complete {
                Some(BondConstraintForm::cis_trans_stereo(
                    CisTransStereoForm::NotStereo,
                ))
            } else {
                None
            }
        }
        BondConstraintKey::RingMembership(scope) => {
            let rings = rings.expect("ring constraint key requires ring context (with_rings)");
            Some(BondConstraintForm::ring_membership(
                scope,
                bond_ring_membership(rings, bond, scope),
            ))
        }
    }
}

/// Mutable borrowed view of a bond.
#[derive(Debug)]
pub struct BondViewMut<'a> {
    pub id: BondId,
    pub atoms: [AtomId; 2],
    pub attributes: &'a mut BondForm,
}

// Editor-scope view bundles for bonds.

pub struct BondEditorView<'a> {
    pub id: BondId,
    pub atoms: [AtomId; 2],
    pub attributes: &'a BondForm,
}

pub struct BondEditorViewMut<'a> {
    pub id: BondId,
    pub atoms: [AtomId; 2],
    pub attributes: &'a mut BondForm,
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;

    use super::super::assert_exact_size_by;
    use crate::ir::aromatic::AromaticSystemForm;
    use crate::ir::atom::AtomForm;
    use crate::ir::bond::BondForm;
    use crate::ir::dative::DativeBondForm;
    use crate::ir::id::{AromaticSystemId, AtomId, BondId, StereoBondId};
    use crate::ir::ligand::{StereoLigand, StereoLigandKind};
    use crate::ir::molecule::{Molecule, MoleculeEntries};
    use crate::ir::multicenter::MulticenterBondForm;
    use crate::ir::noncovalent::{NoncovalentBondForm, NoncovalentBondKind};
    use crate::ir::stereo::{StereoBondForm, StereoCoset, StereoKind};

    #[fixture]
    fn molecule() -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(2)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
            ],
            dative: vec![(vec![AtomId(2)], AtomId(3), DativeBondForm::from_order(1))],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemForm::default(),
            )],
            multicenter: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                MulticenterBondForm::default(),
            )],
            noncovalent: vec![(
                [AtomId(0), AtomId(3)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            ..Default::default()
        })
    }

    #[rstest]
    fn test_bond_views_count(molecule: Molecule) {
        assert_eq!(molecule.bonds().count(), 3);
    }

    #[rstest]
    fn test_bond_views_ids(molecule: Molecule) {
        assert_exact_size_by(Molecule::default().bonds().ids(), vec![], |id| id);
        assert_exact_size_by(
            molecule.bonds().ids(),
            vec![BondId(0), BondId(1), BondId(2)],
            |id| id,
        );
    }

    #[rstest]
    fn test_bond_views_iter(molecule: Molecule) {
        assert_exact_size_by(Molecule::default().bonds().iter(), vec![], |view| {
            (view.id, view.atom_ids(), view.attributes.clone())
        });
        assert_exact_size_by(
            molecule.bonds().iter(),
            vec![
                (BondId(0), [AtomId(0), AtomId(1)], BondForm::from_order(1)),
                (BondId(1), [AtomId(1), AtomId(2)], BondForm::from_order(2)),
                (BondId(2), [AtomId(2), AtomId(3)], BondForm::from_order(1)),
            ],
            |view| (view.id, view.atom_ids(), view.attributes.clone()),
        );
    }

    #[rstest]
    #[case::present(BondId(1), true)]
    #[case::absent(BondId(99), false)]
    fn test_bond_views_contains(molecule: Molecule, #[case] id: BondId, #[case] expected: bool) {
        assert_eq!(molecule.bonds().contains(id), expected);
    }

    #[rstest]
    fn test_bond_views_get(molecule: Molecule) {
        let res = molecule.bonds().get(BondId(1));
        assert!(res.is_some());
        let view = res.unwrap();
        assert_eq!(view.id, BondId(1));
        assert_eq!(view.atom_ids(), [AtomId(1), AtomId(2)]);
        assert_eq!(*view.attributes, BondForm::from_order(2));
    }

    #[rstest]
    fn test_bond_views_get_none(molecule: Molecule) {
        let res = molecule.bonds().get(BondId(99));
        assert!(res.is_none());
    }

    #[rstest]
    fn test_bond_view_atom_ids(molecule: Molecule) {
        assert_eq!(molecule.bond(BondId(1)).atom_ids(), [AtomId(1), AtomId(2)]);
    }

    #[rstest]
    fn test_bond_view_atoms(molecule: Molecule) {
        assert_exact_size_by(
            molecule.bond(BondId(1)).atoms(),
            vec![AtomId(1), AtomId(2)],
            |atom| atom.id,
        );
    }

    #[rstest]
    #[case::both_endpoints_aromatic(BondId(0), Some(AromaticSystemId(0)))]
    #[case::both_endpoints_aromatic_alt(BondId(1), Some(AromaticSystemId(0)))]
    #[case::one_endpoint_outside(BondId(2), None)]
    fn test_bond_view_aromatic_system(
        molecule: Molecule,
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
        molecule: Molecule,
        #[case] bond: BondId,
        #[case] expected: bool,
    ) {
        assert_eq!(molecule.bond(bond).is_in_aromatic_system(), expected);
    }

    #[fixture]
    fn stereo_molecule() -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 4],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(2)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
            ],
            stereo_bonds: vec![(
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
            ..Default::default()
        })
    }

    #[rstest]
    #[case::site(BondId(1), true)]
    #[case::non_site(BondId(0), false)]
    fn test_bond_view_is_stereo_bond(
        stereo_molecule: Molecule,
        #[case] bond: BondId,
        #[case] expected: bool,
    ) {
        assert_eq!(stereo_molecule.bond(bond).is_stereo_bond(), expected);
    }

    #[rstest]
    #[case::site(BondId(1), Some(StereoBondId(0)))]
    #[case::non_site(BondId(0), None)]
    fn test_bond_view_stereo_bond_id(
        stereo_molecule: Molecule,
        #[case] bond: BondId,
        #[case] expected: Option<StereoBondId>,
    ) {
        assert_eq!(stereo_molecule.bond(bond).stereo_bond_id(), expected);
    }

    #[rstest]
    fn test_bond_view_stereo_bond(stereo_molecule: Molecule) {
        let view = stereo_molecule.bond(BondId(1)).stereo_bond().unwrap();
        assert_eq!(view.id, StereoBondId(0));
        assert_eq!(view.kind(), StereoKind::CisTrans);
        assert!(stereo_molecule.bond(BondId(0)).stereo_bond().is_none());
    }
}
