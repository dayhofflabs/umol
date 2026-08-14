//! Dative bond views.

use std::collections::HashSet;

use umol_graph_core::{FixedVarBirelationSet, NodeId, Ordered, RelationId, Unordered};

use super::super::constraint::{
    DativeBondConstraintForm, DativeBondConstraintKey, DativeBondConstraintsForm,
};
use super::super::dative::DativeBondForm;
use super::super::id::{AtomId, DativeBondId};
use super::super::molecule::Molecule;
use super::super::num::NumForm;
use super::super::traits::Lattice;
use super::atom::AtomView;
use super::constraints::DativeBondConstraintsView;

/// Namespace accessor for dative-bond views on a `Molecule`.
#[derive(Clone, Copy)]
pub struct DativeBondViews<'a> {
    molecule: &'a Molecule,
    dative_bonds: &'a FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, DativeBondForm>,
}

impl<'a> DativeBondViews<'a> {
    pub(crate) fn new(
        molecule: &'a Molecule,
        dative_bonds: &'a FixedVarBirelationSet<
            NodeId,
            Ordered,
            1,
            NodeId,
            Unordered,
            DativeBondForm,
        >,
    ) -> Self {
        Self {
            molecule,
            dative_bonds,
        }
    }

    pub fn count(&self) -> usize {
        self.dative_bonds.count()
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = DativeBondId> {
        self.dative_bonds.relation_ids().map(DativeBondId::from)
    }

    pub fn iter(&self) -> impl ExactSizeIterator<Item = DativeBondView<'a>> {
        let molecule = self.molecule;
        let set = self.dative_bonds;
        set.relation_ids().map(move |rid| DativeBondView {
            id: DativeBondId::from(rid),
            attributes: set.data(rid),
            acceptor_id: set.participants_1(rid)[0],
            donors: set.participants_2(rid),
            molecule,
        })
    }

    pub fn contains(&self, id: DativeBondId) -> bool {
        self.dative_bonds.contains(RelationId::from(id))
    }

    pub fn get(&self, id: DativeBondId) -> Option<DativeBondView<'a>> {
        if !self.contains(id) {
            return None;
        }
        let rid = RelationId::from(id);
        Some(DativeBondView {
            id,
            attributes: self.dative_bonds.data(rid),
            acceptor_id: self.dative_bonds.participants_1(rid)[0],
            donors: self.dative_bonds.participants_2(rid),
            molecule: self.molecule,
        })
    }

    /// Ids of dative bonds incident on `atom`.
    pub fn incident_ids(&self, atom: AtomId) -> impl ExactSizeIterator<Item = DativeBondId> + 'a {
        self.dative_bonds
            .incident(NodeId::from(atom))
            .iter()
            .map(|&rid| DativeBondId::from(rid))
    }

    /// Whether any dative bond is incident on `atom`.
    pub fn has_incident(&self, atom: AtomId) -> bool {
        self.dative_bonds.has_incident(NodeId::from(atom))
    }

    /// Views of dative bonds incident on `atom`.
    pub fn incident(&self, atom: AtomId) -> impl ExactSizeIterator<Item = DativeBondView<'a>> + 'a {
        let molecule = self.molecule;
        let set = self.dative_bonds;
        self.incident_ids(atom).map(move |id| {
            let rid = RelationId::from(id);
            DativeBondView {
                id,
                attributes: set.data(rid),
                acceptor_id: set.participants_1(rid)[0],
                donors: set.participants_2(rid),
                molecule,
            }
        })
    }

    /// Id of the dative bond with exactly this acceptor and donor set, if any. Per-factor: the
    /// donor/acceptor roles are matched, not the merged atom set.
    pub fn of_id(&self, acceptor: AtomId, donors: &[AtomId]) -> Option<DativeBondId> {
        let donor_nodes: Vec<NodeId> = donors.iter().map(|&a| NodeId::from(a)).collect();
        self.dative_bonds
            .find_by_participants(&[NodeId::from(acceptor)], &donor_nodes)
            .map(DativeBondId::from)
    }

    /// View of the dative bond with exactly this acceptor and donor set, if any.
    pub fn of(&self, acceptor: AtomId, donors: &[AtomId]) -> Option<DativeBondView<'a>> {
        self.of_id(acceptor, donors).map(|id| {
            self.get(id).expect(
                "dative bond id from relation set must refer to a dative bond in this molecule",
            )
        })
    }

    /// Ids of dative bonds whose participants all lie in `atoms`.
    pub fn induced_ids(&self, atoms: &[AtomId]) -> Vec<DativeBondId> {
        let set: HashSet<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        self.dative_bonds
            .relation_ids()
            .filter(|&rid| {
                self.dative_bonds
                    .participants_1(rid)
                    .iter()
                    .chain(self.dative_bonds.participants_2(rid))
                    .all(|p| set.contains(p))
            })
            .map(DativeBondId::from)
            .collect()
    }

    /// Views of dative bonds whose participants all lie in `atoms`.
    pub fn induced(&self, atoms: &[AtomId]) -> Vec<DativeBondView<'a>> {
        self.induced_ids(atoms)
            .into_iter()
            .map(|id| {
                self.get(id).expect(
                    "dative bond id from relation set must refer to a dative bond in this molecule",
                )
            })
            .collect()
    }
}

/// Borrowed view of a dative bond: index, the designated acceptor atom,
/// and underlying `DativeBondForm`. Donor atoms via `donors()` / `donor_ids()`;
/// the full participant set (donors then acceptor) via `atoms()` / `atom_ids()`.
#[derive(Clone, Copy, Debug)]
pub struct DativeBondView<'a> {
    pub id: DativeBondId,
    acceptor_id: NodeId,
    donors: &'a [NodeId],
    pub attributes: &'a DativeBondForm,
    molecule: &'a Molecule,
}

impl<'a> DativeBondView<'a> {
    #[inline]
    pub fn order(&self) -> &'a NumForm {
        &self.attributes.order
    }

    /// Constraint reading of this dative bond: the container's read API
    /// (asserted side, meanings intact) plus the keyed accessors. Mutation
    /// stays on the stored container.
    #[inline]
    pub fn constraints(&self) -> DativeBondConstraintsView<'a> {
        DativeBondConstraintsView::new(self.molecule, self.id)
    }

    /// Donor atom ids.
    pub fn donor_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + 'a {
        self.donors.iter().map(|&n| AtomId::from(n))
    }

    pub fn acceptor_id(&self) -> AtomId {
        AtomId::from(self.acceptor_id)
    }

    /// All atoms in this dative bond: the donors followed by the acceptor.
    pub fn atom_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + 'a {
        let donors = self.donors;
        let acceptor = self.acceptor_id;
        (0..donors.len() + 1).map(move |index| {
            AtomId::from(if index < donors.len() {
                donors[index]
            } else {
                acceptor
            })
        })
    }

    /// Donor atom views.
    pub fn donors(&self) -> impl ExactSizeIterator<Item = AtomView<'a>> + 'a {
        let molecule = self.molecule;
        self.donor_ids().map(move |id| molecule.atom(id))
    }

    /// View of the acceptor atom.
    pub fn acceptor(&self) -> AtomView<'a> {
        self.molecule.atom(self.acceptor_id())
    }

    /// Views of all atoms in this dative bond (donors then acceptor).
    pub fn atoms(&self) -> impl ExactSizeIterator<Item = AtomView<'a>> + 'a {
        let molecule = self.molecule;
        self.atom_ids().map(move |id| molecule.atom(id))
    }

    pub fn donor_count(&self) -> usize {
        self.donors.len()
    }

    pub fn atom_count(&self) -> usize {
        self.donor_count() + 1
    }

    /// Is dative bond ground
    pub fn is_ground(&self) -> bool {
        self.attributes.is_ground()
    }

    /// Is dative bond undetermined
    pub fn is_undetermined(&self) -> bool {
        self.attributes.is_undetermined()
    }
}

// Derivation layer beneath the dative-bond facades.

/// Stored constraint container of `bond`.
pub(crate) fn dative_bond_asserted_constraints(
    molecule: &Molecule,
    bond: DativeBondId,
) -> &DativeBondConstraintsForm {
    &molecule.dative_bond(bond).attributes.constraints
}

/// Derived side of one dative-bond constraint key. Aromatic incidence is
/// defined only for a binary dative bond (doc 117 stub for multi-donor
/// entries): the donor and acceptor share an aromatic system. The ring key
/// has no projection; both read vacuous under either mode where undefined.
pub(crate) fn dative_bond_derived_constraint(
    molecule: &Molecule,
    bond: DativeBondId,
    key: DativeBondConstraintKey,
    complete: bool,
) -> Option<DativeBondConstraintForm> {
    match key {
        DativeBondConstraintKey::Aromatic => {
            let view = molecule.dative_bond(bond);
            if view.donor_count() != 1 {
                return None;
            }
            let donor_system = view.donors().next().and_then(|d| d.aromatic_system_id());
            let shared =
                donor_system.is_some() && donor_system == view.acceptor().aromatic_system_id();
            if shared {
                Some(DativeBondConstraintForm::aromatic(true))
            } else if complete {
                Some(DativeBondConstraintForm::aromatic(false))
            } else {
                None
            }
        }
        DativeBondConstraintKey::RingMembership(_) => None,
    }
}

/// Mutable borrowed view of a dative bond: its id, participants (donors +
/// acceptor, owned) and mutable data. Molecule-scope peer of `DativeBondView`.
#[derive(Debug)]
pub struct DativeBondViewMut<'a> {
    pub id: DativeBondId,
    pub donors: Vec<AtomId>,
    pub acceptor: AtomId,
    pub attributes: &'a mut DativeBondForm,
}

// Editor-scope view bundles for dative bonds.

pub struct DativeBondEditorView<'a> {
    pub id: DativeBondId,
    donors: &'a [NodeId],
    acceptor: AtomId,
    pub attributes: &'a DativeBondForm,
}

impl<'a> DativeBondEditorView<'a> {
    pub(crate) fn new(
        id: DativeBondId,
        donors: &'a [NodeId],
        acceptor: AtomId,
        attributes: &'a DativeBondForm,
    ) -> Self {
        Self {
            id,
            donors,
            acceptor,
            attributes,
        }
    }

    /// All atoms: donors followed by the acceptor.
    pub fn atom_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + 'a {
        let donors = self.donors;
        let acceptor = self.acceptor;
        (0..donors.len() + 1).map(move |index| {
            if index < donors.len() {
                AtomId::from(donors[index])
            } else {
                acceptor
            }
        })
    }
}

pub struct DativeBondEditorViewMut<'a> {
    pub id: DativeBondId,
    donors: &'a [NodeId],
    acceptor: AtomId,
    pub attributes: &'a mut DativeBondForm,
}

impl<'a> DativeBondEditorViewMut<'a> {
    pub(crate) fn new(
        id: DativeBondId,
        donors: &'a [NodeId],
        acceptor: AtomId,
        attributes: &'a mut DativeBondForm,
    ) -> Self {
        Self {
            id,
            donors,
            acceptor,
            attributes,
        }
    }

    /// All atoms: donors followed by the acceptor.
    pub fn atom_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + '_ {
        let donors = self.donors;
        let acceptor = self.acceptor;
        (0..donors.len() + 1).map(move |index| {
            if index < donors.len() {
                AtomId::from(donors[index])
            } else {
                acceptor
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;
    use umol_graph_core::NodeId;

    use super::super::assert_exact_size_by;
    use super::{DativeBondEditorView, DativeBondEditorViewMut};
    use crate::ir::aromatic::AromaticSystemForm;
    use crate::ir::atom::AtomForm;
    use crate::ir::bond::BondForm;
    use crate::ir::dative::DativeBondForm;
    use crate::ir::id::{AtomId, DativeBondId};
    use crate::ir::molecule::{Molecule, MoleculeEntries};
    use crate::ir::multicenter::MulticenterBondForm;
    use crate::ir::noncovalent::{NoncovalentBondForm, NoncovalentBondKind};

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
                AtomId(0),
                AtomId(3),
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            ..Default::default()
        })
    }

    #[rstest]
    fn test_dative_bond_views_count(molecule: Molecule) {
        assert_eq!(molecule.dative_bonds().count(), 1);
    }

    #[rstest]
    fn test_dative_bond_views_ids(molecule: Molecule) {
        assert_exact_size_by(Molecule::default().dative_bonds().ids(), vec![], |id| id);
        assert_exact_size_by(molecule.dative_bonds().ids(), vec![DativeBondId(0)], |id| {
            id
        });
    }

    #[rstest]
    fn test_dative_bond_views_iter(molecule: Molecule) {
        assert_exact_size_by(Molecule::default().dative_bonds().iter(), vec![], |view| {
            (view.id, view.acceptor_id(), view.attributes.clone())
        });
        assert_exact_size_by(
            molecule.dative_bonds().iter(),
            vec![(DativeBondId(0), AtomId(3), DativeBondForm::from_order(1))],
            |view| (view.id, view.acceptor_id(), view.attributes.clone()),
        );
    }

    #[rstest]
    #[case::participant(AtomId(2), vec![DativeBondId(0)])]
    #[case::uninvolved(AtomId(0), vec![])]
    fn test_dative_bond_views_incident(
        molecule: Molecule,
        #[case] atom: AtomId,
        #[case] expected: Vec<DativeBondId>,
    ) {
        assert_exact_size_by(
            molecule.dative_bonds().incident_ids(atom),
            expected.clone(),
            |id| id,
        );
        assert_exact_size_by(molecule.dative_bonds().incident(atom), expected, |view| {
            view.id
        });
    }

    #[rstest]
    #[case::present(DativeBondId(0), true)]
    #[case::absent(DativeBondId(99), false)]
    fn test_dative_bond_views_contains(
        molecule: Molecule,
        #[case] id: DativeBondId,
        #[case] expected: bool,
    ) {
        assert_eq!(molecule.dative_bonds().contains(id), expected);
    }

    #[rstest]
    fn test_dative_bond_views_get(molecule: Molecule) {
        let res = molecule.dative_bonds().get(DativeBondId(0));
        assert!(res.is_some());
        let view = res.unwrap();
        assert_eq!(view.id, DativeBondId(0));
        assert_eq!(view.acceptor_id(), AtomId(3));
    }

    #[rstest]
    fn test_dative_bond_views_get_none(molecule: Molecule) {
        let res = molecule.dative_bonds().get(DativeBondId(99));
        assert!(res.is_none());
    }

    #[rstest]
    fn test_dative_bond_view_atom_ids(molecule: Molecule) {
        assert_exact_size_by(
            molecule.dative_bond(DativeBondId(0)).atom_ids(),
            vec![AtomId(2), AtomId(3)],
            |id| id,
        );
    }

    #[rstest]
    fn test_dative_bond_view_donor_ids(molecule: Molecule) {
        assert_exact_size_by(
            molecule.dative_bond(DativeBondId(0)).donor_ids(),
            vec![AtomId(2)],
            |id| id,
        );
    }

    #[rstest]
    fn test_dative_bond_view_acceptor_id(molecule: Molecule) {
        assert_eq!(
            molecule.dative_bond(DativeBondId(0)).acceptor_id(),
            AtomId(3)
        );
    }

    #[rstest]
    fn test_dative_bond_view_atoms(molecule: Molecule) {
        assert_exact_size_by(
            molecule.dative_bond(DativeBondId(0)).atoms(),
            vec![AtomId(2), AtomId(3)],
            |atom| atom.id,
        );
    }

    #[rstest]
    fn test_dative_bond_view_donors(molecule: Molecule) {
        assert_exact_size_by(
            molecule.dative_bond(DativeBondId(0)).donors(),
            vec![AtomId(2)],
            |atom| atom.id,
        );
    }

    #[rstest]
    fn test_dative_bond_view_acceptor(molecule: Molecule) {
        assert_eq!(
            molecule.dative_bond(DativeBondId(0)).acceptor().id,
            AtomId(3),
        );
    }

    #[rstest]
    fn test_dative_bond_view_atom_count(molecule: Molecule) {
        assert_eq!(molecule.dative_bond(DativeBondId(0)).atom_count(), 2);
    }

    #[rstest]
    fn test_dative_bond_editor_view_atom_ids() {
        let donors = [NodeId(1), NodeId(2)];
        let attributes = DativeBondForm::from_order(1);
        let view = DativeBondEditorView::new(DativeBondId(0), &donors, AtomId(3), &attributes);
        assert_exact_size_by(
            view.atom_ids(),
            vec![AtomId(1), AtomId(2), AtomId(3)],
            |id| id,
        );
    }

    #[rstest]
    fn test_dative_bond_editor_view_mut_atom_ids() {
        let donors = [NodeId(1), NodeId(2)];
        let mut attributes = DativeBondForm::from_order(1);
        let view =
            DativeBondEditorViewMut::new(DativeBondId(0), &donors, AtomId(3), &mut attributes);
        assert_exact_size_by(
            view.atom_ids(),
            vec![AtomId(1), AtomId(2), AtomId(3)],
            |id| id,
        );
    }
}
