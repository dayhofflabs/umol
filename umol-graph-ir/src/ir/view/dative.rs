//! Dative bond views.

use std::collections::HashSet;
use std::iter;

use umol_graph_core::NodeId;

use super::super::constraint::{
    DativeBondConstraintForm, DativeBondConstraintKey, DativeBondConstraintsForm,
};
use super::super::dative::{DativeBondForm, DativeBonds};
use super::super::id::{AtomId, AtomPosition, DativeBondId};
use super::super::molecule::Molecule;
use super::super::num::NumForm;
use super::super::traits::Lattice;
use super::atom::AtomView;
use super::constraints::DativeBondConstraintsView;

/// Namespace accessor for dative-bond views on a `Molecule`.
#[derive(Clone, Copy)]
pub struct DativeBondViews<'a> {
    molecule: &'a Molecule,
}

impl<'a> DativeBondViews<'a> {
    pub(crate) fn new(molecule: &'a Molecule) -> Self {
        Self { molecule }
    }

    pub fn count(&self) -> usize {
        self.molecule.raw_dative_bonds().count()
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = DativeBondId> {
        self.molecule.raw_dative_bonds().ids()
    }

    pub fn iter(&self) -> impl ExactSizeIterator<Item = DativeBondView<'a>> {
        let molecule = self.molecule;
        self.molecule
            .raw_dative_bonds()
            .ids()
            .map(move |id| DativeBondView { molecule, id })
    }

    pub fn contains(&self, id: DativeBondId) -> bool {
        self.molecule.raw_dative_bonds().contains(id)
    }

    pub fn get(&self, id: DativeBondId) -> Option<DativeBondView<'a>> {
        if !self.contains(id) {
            return None;
        }
        Some(DativeBondView {
            molecule: self.molecule,
            id,
        })
    }

    /// Ids of dative bonds incident on `atom`.
    pub fn incident_ids(&self, atom: AtomId) -> impl ExactSizeIterator<Item = DativeBondId> + 'a {
        self.molecule.raw_dative_bonds().incident_ids(atom)
    }

    /// Whether any dative bond is incident on `atom`.
    pub fn has_incident(&self, atom: AtomId) -> bool {
        self.molecule.raw_dative_bonds().has_incident(atom)
    }

    /// Views of dative bonds incident on `atom`.
    pub fn incident(&self, atom: AtomId) -> impl ExactSizeIterator<Item = DativeBondView<'a>> + 'a {
        let molecule = self.molecule;
        self.incident_ids(atom)
            .map(move |id| DativeBondView { molecule, id })
    }

    /// Id of the dative bond with exactly this acceptor and donor set, if any. Per-factor: the
    /// donor/acceptor roles are matched, not the merged atom set.
    pub fn of_id(&self, acceptor: AtomId, donors: &[AtomId]) -> Option<DativeBondId> {
        self.molecule
            .raw_dative_bonds()
            .coincident_id(acceptor, donors)
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
        let dative_bonds = self.molecule.raw_dative_bonds();
        dative_bonds
            .ids()
            .filter(|&id| {
                iter::once(&dative_bonds.acceptor_node(id))
                    .chain(dative_bonds.donor_nodes(id))
                    .all(|p| set.contains(p))
            })
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
    molecule: &'a Molecule,
    id: DativeBondId,
}

impl<'a> DativeBondView<'a> {
    #[inline]
    pub fn id(&self) -> DativeBondId {
        self.id
    }

    #[inline]
    pub fn attributes(&self) -> &'a DativeBondForm {
        self.molecule.raw_dative_bonds().attributes(self.id)
    }

    #[inline]
    pub fn order(&self) -> &'a NumForm {
        &self.attributes().order
    }

    /// Constraint reading of this dative bond: the container's read API
    /// (asserted side, meanings intact) plus the keyed accessors. Mutation
    /// stays on the stored container.
    #[inline]
    pub fn constraints(&self) -> DativeBondConstraintsView<'a> {
        DativeBondConstraintsView::new(self.molecule, self.id)
    }

    /// Donor atom ids.
    #[inline]
    pub fn donor_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + 'a {
        self.molecule.raw_dative_bonds().donors(self.id)
    }

    #[inline]
    pub fn acceptor_id(&self) -> AtomId {
        self.molecule.raw_dative_bonds().acceptor(self.id)
    }

    /// All atoms in this dative bond: the donors followed by the acceptor.
    #[inline]
    pub fn atom_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + 'a {
        dative_bond_atom_ids(self.molecule.raw_dative_bonds(), self.id)
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

    #[inline]
    pub fn donor_count(&self) -> usize {
        self.donor_ids().len()
    }

    #[inline]
    pub fn atom_count(&self) -> usize {
        self.donor_count() + 1
    }

    /// Is dative bond ground
    pub fn is_ground(&self) -> bool {
        self.attributes().is_ground()
    }

    /// Is dative bond undetermined
    pub fn is_undetermined(&self) -> bool {
        self.attributes().is_undetermined()
    }
}

/// Read-only editor access to a dative bond.
pub struct DativeBondEditorView<'a> {
    dative_bonds: &'a DativeBonds,
    id: DativeBondId,
}

impl<'a> DativeBondEditorView<'a> {
    pub(crate) fn new(dative_bonds: &'a DativeBonds, id: DativeBondId) -> Self {
        Self { dative_bonds, id }
    }

    #[inline]
    pub fn id(&self) -> DativeBondId {
        self.id
    }

    #[inline]
    pub fn attributes(&self) -> &'a DativeBondForm {
        self.dative_bonds.attributes(self.id)
    }

    #[inline]
    pub fn order(&self) -> &'a NumForm {
        &self.attributes().order
    }

    /// Donor atom ids.
    #[inline]
    pub fn donor_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + 'a {
        self.dative_bonds.donors(self.id)
    }

    #[inline]
    pub fn acceptor_id(&self) -> AtomId {
        self.dative_bonds.acceptor(self.id)
    }

    /// All atoms: donors followed by the acceptor.
    #[inline]
    pub fn atom_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + 'a {
        dative_bond_atom_ids(self.dative_bonds, self.id)
    }

    #[inline]
    pub fn donor_count(&self) -> usize {
        self.donor_ids().len()
    }

    #[inline]
    pub fn atom_count(&self) -> usize {
        self.donor_count() + 1
    }
}

/// Mutable attribute access to a dative bond.
#[derive(Debug)]
pub struct DativeBondViewMut<'a> {
    dative_bonds: &'a mut DativeBonds,
    id: DativeBondId,
}

impl<'a> DativeBondViewMut<'a> {
    pub(crate) fn new(dative_bonds: &'a mut DativeBonds, id: DativeBondId) -> Self {
        Self { dative_bonds, id }
    }

    #[inline]
    pub fn id(&self) -> DativeBondId {
        self.id
    }

    #[inline]
    pub fn attributes(&self) -> &DativeBondForm {
        self.dative_bonds.attributes(self.id)
    }

    #[inline]
    pub fn attributes_mut(&mut self) -> &mut DativeBondForm {
        self.dative_bonds.attributes_mut(self.id)
    }

    #[inline]
    pub fn order(&self) -> &NumForm {
        &self.attributes().order
    }

    #[inline]
    pub fn constraints(&self) -> &DativeBondConstraintsForm {
        &self.attributes().constraints
    }

    /// Donor atom ids.
    #[inline]
    pub fn donor_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + '_ {
        self.dative_bonds.donors(self.id)
    }

    #[inline]
    pub fn acceptor_id(&self) -> AtomId {
        self.dative_bonds.acceptor(self.id)
    }

    /// All atoms: donors followed by the acceptor.
    #[inline]
    pub fn atom_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + '_ {
        dative_bond_atom_ids(self.dative_bonds, self.id)
    }

    #[inline]
    pub fn donor_count(&self) -> usize {
        self.donor_ids().len()
    }

    #[inline]
    pub fn atom_count(&self) -> usize {
        self.donor_count() + 1
    }
}

/// Mutable editor access to a dative bond.
///
/// Structural mutations update incidence and preserve attributes and constraints.
/// Donor mutations preserve the acceptor; acceptor replacement preserves donors.
/// Atom existence and distinctness are checked at molecule publication.
#[derive(Debug)]
pub struct DativeBondEditorViewMut<'a> {
    dative_bonds: &'a mut DativeBonds,
    id: DativeBondId,
}

impl<'a> DativeBondEditorViewMut<'a> {
    pub(crate) fn new(dative_bonds: &'a mut DativeBonds, id: DativeBondId) -> Self {
        Self { dative_bonds, id }
    }

    #[inline]
    pub fn id(&self) -> DativeBondId {
        self.id
    }

    #[inline]
    pub fn attributes(&self) -> &DativeBondForm {
        self.dative_bonds.attributes(self.id)
    }

    #[inline]
    pub fn attributes_mut(&mut self) -> &mut DativeBondForm {
        self.dative_bonds.attributes_mut(self.id)
    }

    #[inline]
    pub fn order(&self) -> &NumForm {
        &self.attributes().order
    }

    #[inline]
    pub fn constraints(&self) -> &DativeBondConstraintsForm {
        &self.attributes().constraints
    }

    /// Donor atom ids.
    #[inline]
    pub fn donor_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + '_ {
        self.dative_bonds.donors(self.id)
    }

    #[inline]
    pub fn acceptor_id(&self) -> AtomId {
        self.dative_bonds.acceptor(self.id)
    }

    /// All atoms: donors followed by the acceptor.
    #[inline]
    pub fn atom_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + '_ {
        dative_bond_atom_ids(self.dative_bonds, self.id)
    }

    #[inline]
    pub fn donor_count(&self) -> usize {
        self.donor_ids().len()
    }

    #[inline]
    pub fn atom_count(&self) -> usize {
        self.donor_count() + 1
    }

    /// Replace the donor list, preserving the supplied order.
    pub fn replace_donors(&mut self, donors: &[AtomId]) {
        self.dative_bonds.replace_donors(self.id, donors);
    }

    /// Replace the acceptor atom.
    pub fn replace_acceptor(&mut self, acceptor: AtomId) {
        self.dative_bonds.replace_acceptor(self.id, acceptor);
    }

    /// Replace the donor at `position` without changing the donor count.
    ///
    /// # Panics
    ///
    /// Panics when `position` is outside the donor list.
    pub fn replace_donor(&mut self, position: AtomPosition, donor: AtomId) {
        self.dative_bonds
            .replace_donor(self.id, position.index(), donor);
    }

    /// Insert a donor at `position`, preserving the order of the existing donors.
    ///
    /// # Panics
    ///
    /// Panics when `position` exceeds the donor count. Insertion at the end is allowed.
    pub fn insert_donor(&mut self, position: AtomPosition, donor: AtomId) {
        self.dative_bonds
            .insert_donor(self.id, position.index(), donor);
    }

    /// Remove the donor at `position`, preserving the order of the remaining donors.
    ///
    /// # Panics
    ///
    /// Panics when `position` is outside the donor list.
    pub fn remove_donor(&mut self, position: AtomPosition) {
        self.dative_bonds.remove_donor(self.id, position.index());
    }
}

#[inline]
fn dative_bond_atom_ids(
    dative_bonds: &DativeBonds,
    id: DativeBondId,
) -> impl ExactSizeIterator<Item = AtomId> + '_ {
    let donors = dative_bonds.donor_nodes(id);
    let acceptor = dative_bonds.acceptor(id);
    (0..donors.len() + 1).map(move |index| {
        if index < donors.len() {
            AtomId::from(donors[index])
        } else {
            acceptor
        }
    })
}

/// Stored constraint container of `bond`.
pub(crate) fn dative_bond_asserted_constraints(
    molecule: &Molecule,
    id: DativeBondId,
) -> &DativeBondConstraintsForm {
    &molecule.dative_bond(id).attributes().constraints
}

/// Asserted side of one dative-bond constraint key under resolution's
/// closed-world claim: the stored assertion, else the absence cell closed to
/// its definite negative. Never reads relations.
pub(crate) fn dative_bond_asserted_complete_constraint(
    molecule: &Molecule,
    id: DativeBondId,
    key: DativeBondConstraintKey,
) -> Option<DativeBondConstraintForm> {
    if let Some(asserted) = dative_bond_asserted_constraints(molecule, id).get(key) {
        return Some(asserted.clone());
    }
    match key {
        DativeBondConstraintKey::Aromatic => Some(DativeBondConstraintForm::aromatic(false)),
        DativeBondConstraintKey::RingMembership(_) => None,
    }
}

/// Derived side of one dative-bond constraint key. Aromatic incidence is
/// defined only for a binary dative bond (doc 117 stub for multi-donor
/// entries): the donor and acceptor share an aromatic system. The ring key
/// has no projection; both read vacuous under either mode where undefined.
pub(crate) fn dative_bond_derived_constraint(
    molecule: &Molecule,
    id: DativeBondId,
    key: DativeBondConstraintKey,
    complete: bool,
) -> Option<DativeBondConstraintForm> {
    match key {
        DativeBondConstraintKey::Aromatic => {
            let view = molecule.dative_bond(id);
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

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;

    use super::super::assert_exact_size_by;
    use crate::ir::aromatic::AromaticSystemForm;
    use crate::ir::atom::AtomForm;
    use crate::ir::bond::BondForm;
    use crate::ir::constraint::DativeBondConstraintForm;
    use crate::ir::dative::DativeBondForm;
    use crate::ir::entity::Entity;
    use crate::ir::id::{AtomId, DativeBondId};
    use crate::ir::molecule::{Molecule, MoleculeEntries};
    use crate::ir::multicenter::MulticenterBondForm;
    use crate::ir::noncovalent::{NoncovalentBondForm, NoncovalentBondKind};
    use crate::ir::num::NumForm;
    use crate::ir::{AtomPosition, MoleculeIntegrityError};

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

    #[fixture]
    fn dative_molecule() -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::default(); 4],
            dative: vec![
                (vec![AtomId(1)], AtomId(0), DativeBondForm::from_order(1)),
                (
                    vec![AtomId(3), AtomId(1)],
                    AtomId(2),
                    DativeBondForm::from_order(2),
                ),
                (vec![], AtomId(3), DativeBondForm::default()),
            ],
            ..Default::default()
        })
    }

    #[fixture]
    fn dative_entries() -> MoleculeEntries {
        MoleculeEntries {
            atoms: vec![AtomForm::default(); 4],
            dative: vec![(
                vec![AtomId(1), AtomId(2)],
                AtomId(0),
                DativeBondForm::from_order(2)
                    .with_constraint(DativeBondConstraintForm::aromatic(true)),
            )],
            ..Default::default()
        }
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
            (view.id(), view.acceptor_id(), view.attributes().clone())
        });
        assert_exact_size_by(
            molecule.dative_bonds().iter(),
            vec![(DativeBondId(0), AtomId(3), DativeBondForm::from_order(1))],
            |view| (view.id(), view.acceptor_id(), view.attributes().clone()),
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
            view.id()
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
        assert_eq!(view.id(), DativeBondId(0));
        assert_eq!(view.acceptor_id(), AtomId(3));
    }

    #[rstest]
    fn test_dative_bond_views_get_none(molecule: Molecule) {
        let res = molecule.dative_bonds().get(DativeBondId(99));
        assert!(res.is_none());
    }

    #[rstest]
    #[case(DativeBondId(0))]
    #[case(DativeBondId(1))]
    #[case(DativeBondId(2))]
    fn test_dative_bond_view_id(
        #[from(dative_molecule)] molecule: Molecule,
        #[case] id: DativeBondId,
    ) {
        assert_eq!(molecule.dative_bond(id).id(), id);
    }

    #[rstest]
    #[case(DativeBondId(0), DativeBondForm { order: NumForm::Lit(1), ..Default::default() })]
    #[case(DativeBondId(1), DativeBondForm { order: NumForm::Lit(2), ..Default::default() })]
    #[case(DativeBondId(2), DativeBondForm::default())]
    fn test_dative_bond_view_attributes(
        #[from(dative_molecule)] molecule: Molecule,
        #[case] id: DativeBondId,
        #[case] expected: DativeBondForm,
    ) {
        let attributes = {
            let view = molecule.dative_bond(id);
            view.attributes()
        };
        assert_eq!(attributes, &expected);
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
            |atom| atom.id(),
        );
    }

    #[rstest]
    fn test_dative_bond_view_donors(molecule: Molecule) {
        assert_exact_size_by(
            molecule.dative_bond(DativeBondId(0)).donors(),
            vec![AtomId(2)],
            |atom| atom.id(),
        );
    }

    #[rstest]
    fn test_dative_bond_view_acceptor(molecule: Molecule) {
        assert_eq!(
            molecule.dative_bond(DativeBondId(0)).acceptor().id(),
            AtomId(3),
        );
    }

    #[rstest]
    fn test_dative_bond_view_atom_count(molecule: Molecule) {
        assert_eq!(molecule.dative_bond(DativeBondId(0)).atom_count(), 2);
    }

    #[rstest]
    #[case(DativeBondId(0))]
    #[case(DativeBondId(1))]
    #[case(DativeBondId(2))]
    fn test_dative_bond_editor_view_id(
        #[from(dative_molecule)] molecule: Molecule,
        #[case] id: DativeBondId,
    ) {
        let editor = molecule.edit();
        assert_eq!(editor.dative_bond(id).id(), id);
    }

    #[rstest]
    #[case(DativeBondId(0), DativeBondForm { order: NumForm::Lit(1), ..Default::default() })]
    #[case(DativeBondId(1), DativeBondForm { order: NumForm::Lit(2), ..Default::default() })]
    #[case(DativeBondId(2), DativeBondForm::default())]
    fn test_dative_bond_editor_view_attributes(
        #[from(dative_molecule)] molecule: Molecule,
        #[case] id: DativeBondId,
        #[case] expected: DativeBondForm,
    ) {
        let editor = molecule.edit();
        let attributes = {
            let view = editor.dative_bond(id);
            view.attributes()
        };
        assert_eq!(attributes, &expected);
    }

    #[rstest]
    #[case(DativeBondId(0), NumForm::Lit(1))]
    #[case(DativeBondId(1), NumForm::Lit(2))]
    #[case(DativeBondId(2), NumForm::Undetermined)]
    fn test_dative_bond_editor_view_order(
        #[from(dative_molecule)] molecule: Molecule,
        #[case] id: DativeBondId,
        #[case] expected: NumForm,
    ) {
        let editor = molecule.edit();
        let order = {
            let view = editor.dative_bond(id);
            view.order()
        };
        assert_eq!(order, &expected);
    }

    #[rstest]
    #[case::single(DativeBondId(0), vec![AtomId(1)])]
    #[case::multiple(DativeBondId(1), vec![AtomId(3), AtomId(1)])]
    #[case::empty(DativeBondId(2), vec![])]
    fn test_dative_bond_editor_view_donor_ids(
        #[from(dative_molecule)] molecule: Molecule,
        #[case] id: DativeBondId,
        #[case] expected: Vec<AtomId>,
    ) {
        let editor = molecule.edit();
        let view = editor.dative_bond(id);
        assert_eq!(view.donor_count(), expected.len());
        assert_eq!(view.atom_count(), expected.len() + 1);
        assert_exact_size_by(view.donor_ids(), expected, |id| id);
    }

    #[rstest]
    #[case(DativeBondId(0), AtomId(0))]
    #[case(DativeBondId(1), AtomId(2))]
    #[case(DativeBondId(2), AtomId(3))]
    fn test_dative_bond_editor_view_acceptor_id(
        #[from(dative_molecule)] molecule: Molecule,
        #[case] id: DativeBondId,
        #[case] expected: AtomId,
    ) {
        let editor = molecule.edit();
        let view = editor.dative_bond(id);
        assert_eq!(view.acceptor_id(), expected);
    }

    #[rstest]
    fn test_dative_bond_editor_view_atom_ids() {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::default(); 4],
            dative: vec![(
                vec![AtomId(1), AtomId(2)],
                AtomId(3),
                DativeBondForm::from_order(1),
            )],
            ..Default::default()
        });
        let editor = molecule.edit();
        let view = editor.dative_bond(DativeBondId(0));
        assert_exact_size_by(
            view.atom_ids(),
            vec![AtomId(1), AtomId(2), AtomId(3)],
            |id| id,
        );
    }

    #[rstest]
    #[case(DativeBondId(0))]
    #[case(DativeBondId(1))]
    fn test_dative_bond_view_mut_attributes_mut(
        #[from(dative_molecule)] mut molecule: Molecule,
        #[case] id: DativeBondId,
    ) {
        let expected = DativeBondForm {
            order: NumForm::Lit(3),
            ..Default::default()
        };
        {
            let mut view = molecule.dative_bond_mut(id);
            assert_eq!(view.id(), id);
            *view.attributes_mut() = expected.clone();
            assert_eq!(view.attributes(), &expected);
            assert_eq!(view.order(), &NumForm::Lit(3));
        }
        assert_eq!(molecule.dative_bond(id).attributes(), &expected);
    }

    #[rstest]
    #[case::single(DativeBondId(0), vec![AtomId(1)])]
    #[case::multiple(DativeBondId(1), vec![AtomId(3), AtomId(1)])]
    #[case::empty(DativeBondId(2), vec![])]
    fn test_dative_bond_view_mut_donor_ids(
        #[from(dative_molecule)] mut molecule: Molecule,
        #[case] id: DativeBondId,
        #[case] expected: Vec<AtomId>,
    ) {
        let view = molecule.dative_bond_mut(id);
        assert_eq!(view.donor_count(), expected.len());
        assert_eq!(view.atom_count(), expected.len() + 1);
        assert_exact_size_by(view.donor_ids(), expected, |id| id);
    }

    #[rstest]
    #[case(DativeBondId(0), AtomId(0))]
    #[case(DativeBondId(1), AtomId(2))]
    #[case(DativeBondId(2), AtomId(3))]
    fn test_dative_bond_view_mut_acceptor_id(
        #[from(dative_molecule)] mut molecule: Molecule,
        #[case] id: DativeBondId,
        #[case] expected: AtomId,
    ) {
        let view = molecule.dative_bond_mut(id);
        assert_eq!(view.acceptor_id(), expected);
    }

    #[rstest]
    #[case(DativeBondId(0))]
    #[case(DativeBondId(1))]
    fn test_dative_bond_editor_view_mut_attributes_mut(
        #[from(dative_molecule)] molecule: Molecule,
        #[case] id: DativeBondId,
    ) {
        let mut editor = molecule.edit();
        let expected = DativeBondForm {
            order: NumForm::Lit(3),
            ..Default::default()
        };
        {
            let mut view = editor.dative_bond_mut(id);
            assert_eq!(view.id(), id);
            *view.attributes_mut() = expected.clone();
            assert_eq!(view.attributes(), &expected);
            assert_eq!(view.order(), &NumForm::Lit(3));
        }
        assert_eq!(editor.dative_bond(id).attributes(), &expected);
    }

    #[rstest]
    #[case::single(DativeBondId(0), vec![AtomId(1)])]
    #[case::multiple(DativeBondId(1), vec![AtomId(3), AtomId(1)])]
    #[case::empty(DativeBondId(2), vec![])]
    fn test_dative_bond_editor_view_mut_donor_ids(
        #[from(dative_molecule)] molecule: Molecule,
        #[case] id: DativeBondId,
        #[case] expected: Vec<AtomId>,
    ) {
        let mut editor = molecule.edit();
        let view = editor.dative_bond_mut(id);
        assert_eq!(view.donor_count(), expected.len());
        assert_eq!(view.atom_count(), expected.len() + 1);
        assert_exact_size_by(view.donor_ids(), expected, |id| id);
    }

    #[rstest]
    #[case(DativeBondId(0), AtomId(0))]
    #[case(DativeBondId(1), AtomId(2))]
    #[case(DativeBondId(2), AtomId(3))]
    fn test_dative_bond_editor_view_mut_acceptor_id(
        #[from(dative_molecule)] molecule: Molecule,
        #[case] id: DativeBondId,
        #[case] expected: AtomId,
    ) {
        let mut editor = molecule.edit();
        let view = editor.dative_bond_mut(id);
        assert_eq!(view.acceptor_id(), expected);
    }

    #[rstest]
    fn test_dative_bond_editor_view_mut_atom_ids() {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::default(); 4],
            dative: vec![(
                vec![AtomId(1), AtomId(2)],
                AtomId(3),
                DativeBondForm::from_order(1),
            )],
            ..Default::default()
        });
        let mut editor = molecule.edit();
        let view = editor.dative_bond_mut(DativeBondId(0));
        assert_exact_size_by(
            view.atom_ids(),
            vec![AtomId(1), AtomId(2), AtomId(3)],
            |id| id,
        );
    }

    #[rstest]
    #[case::empty(vec![])]
    #[case::single(vec![AtomId(3)])]
    #[case::reordered(vec![AtomId(2), AtomId(1)])]
    #[case::expanded(vec![AtomId(3), AtomId(2), AtomId(1)])]
    fn test_dative_bond_editor_view_mut_replace_donors(
        mut dative_entries: MoleculeEntries,
        #[case] donors: Vec<AtomId>,
    ) {
        let mut editor = Molecule::from_entries(dative_entries.clone()).edit();
        editor
            .dative_bond_mut(DativeBondId(0))
            .replace_donors(&donors);
        assert_eq!(
            editor
                .dative_bond(DativeBondId(0))
                .donor_ids()
                .collect::<Vec<_>>(),
            donors
        );
        dative_entries.dative[0].0 = donors;
        assert_eq!(
            editor.try_build(),
            Ok(Molecule::from_entries(dative_entries))
        );
    }

    #[rstest]
    fn test_dative_bond_editor_view_mut_replace_acceptor(mut dative_entries: MoleculeEntries) {
        let mut editor = Molecule::from_entries(dative_entries.clone()).edit();
        editor
            .dative_bond_mut(DativeBondId(0))
            .replace_acceptor(AtomId(3));
        dative_entries.dative[0].1 = AtomId(3);
        let molecule = editor.try_build().unwrap();
        assert_eq!(molecule, Molecule::from_entries(dative_entries));
        assert_eq!(
            molecule
                .dative_bonds()
                .incident_ids(AtomId(0))
                .collect::<Vec<_>>(),
            vec![]
        );
        assert_eq!(
            molecule
                .dative_bonds()
                .incident_ids(AtomId(3))
                .collect::<Vec<_>>(),
            vec![DativeBondId(0)]
        );
    }

    #[rstest]
    #[case::first(AtomPosition(0), vec![AtomId(3), AtomId(2)])]
    #[case::last(AtomPosition(1), vec![AtomId(1), AtomId(3)])]
    fn test_dative_bond_editor_view_mut_replace_donor(
        mut dative_entries: MoleculeEntries,
        #[case] position: AtomPosition,
        #[case] expected: Vec<AtomId>,
    ) {
        let mut editor = Molecule::from_entries(dative_entries.clone()).edit();
        editor
            .dative_bond_mut(DativeBondId(0))
            .replace_donor(position, AtomId(3));
        assert_eq!(
            editor
                .dative_bond(DativeBondId(0))
                .donor_ids()
                .collect::<Vec<_>>(),
            expected
        );
        dative_entries.dative[0].0 = expected;
        assert_eq!(
            editor.try_build(),
            Ok(Molecule::from_entries(dative_entries))
        );
    }

    #[rstest]
    #[case::first(vec![AtomId(1), AtomId(2)], AtomPosition(0), vec![AtomId(3), AtomId(1), AtomId(2)])]
    #[case::middle(vec![AtomId(1), AtomId(2)], AtomPosition(1), vec![AtomId(1), AtomId(3), AtomId(2)])]
    #[case::end(vec![AtomId(1), AtomId(2)], AtomPosition(2), vec![AtomId(1), AtomId(2), AtomId(3)])]
    #[case::empty(vec![], AtomPosition(0), vec![AtomId(3)])]
    fn test_dative_bond_editor_view_mut_insert_donor(
        mut dative_entries: MoleculeEntries,
        #[case] donors: Vec<AtomId>,
        #[case] position: AtomPosition,
        #[case] expected: Vec<AtomId>,
    ) {
        dative_entries.dative[0].0 = donors;
        let mut editor = Molecule::from_entries(dative_entries.clone()).edit();
        editor
            .dative_bond_mut(DativeBondId(0))
            .insert_donor(position, AtomId(3));
        assert_eq!(
            editor
                .dative_bond(DativeBondId(0))
                .donor_ids()
                .collect::<Vec<_>>(),
            expected
        );
        dative_entries.dative[0].0 = expected;
        assert_eq!(
            editor.try_build(),
            Ok(Molecule::from_entries(dative_entries))
        );
    }

    #[rstest]
    #[case::first(vec![AtomId(1), AtomId(2)], AtomPosition(0), vec![AtomId(2)])]
    #[case::last(vec![AtomId(1), AtomId(2)], AtomPosition(1), vec![AtomId(1)])]
    #[case::only(vec![AtomId(1)], AtomPosition(0), vec![])]
    fn test_dative_bond_editor_view_mut_remove_donor(
        mut dative_entries: MoleculeEntries,
        #[case] donors: Vec<AtomId>,
        #[case] position: AtomPosition,
        #[case] expected: Vec<AtomId>,
    ) {
        dative_entries.dative[0].0 = donors;
        let mut editor = Molecule::from_entries(dative_entries.clone()).edit();
        editor
            .dative_bond_mut(DativeBondId(0))
            .remove_donor(position);
        assert_eq!(
            editor
                .dative_bond(DativeBondId(0))
                .donor_ids()
                .collect::<Vec<_>>(),
            expected
        );
        dative_entries.dative[0].0 = expected;
        assert_eq!(
            editor.try_build(),
            Ok(Molecule::from_entries(dative_entries))
        );
    }

    #[rstest]
    fn test_dative_bond_editor_view_mut_replace_donors_incidence(dative_entries: MoleculeEntries) {
        let mut editor = Molecule::from_entries(dative_entries).edit();
        {
            let mut view = editor.dative_bond_mut(DativeBondId(0));
            view.replace_donors(&[AtomId(3), AtomId(1)]);
            view.replace_donor(AtomPosition(1), AtomId(2));
            view.insert_donor(AtomPosition(2), AtomId(1));
            view.remove_donor(AtomPosition(0));
            view.replace_acceptor(AtomId(3));
        }
        let molecule = editor.try_build().unwrap();
        let incidence: Vec<Vec<_>> = molecule
            .atoms()
            .ids()
            .map(|atom| molecule.dative_bonds().incident_ids(atom).collect())
            .collect();
        assert_eq!(
            incidence,
            vec![
                vec![],
                vec![DativeBondId(0)],
                vec![DativeBondId(0)],
                vec![DativeBondId(0)]
            ]
        );
    }

    #[rstest]
    #[case::empty(vec![], AtomPosition(0))]
    #[case::end(vec![AtomId(1)], AtomPosition(1))]
    #[should_panic]
    fn test_dative_bond_editor_view_mut_replace_donor_error(
        mut dative_entries: MoleculeEntries,
        #[case] donors: Vec<AtomId>,
        #[case] position: AtomPosition,
    ) {
        dative_entries.dative[0].0 = donors;
        let mut editor = Molecule::from_entries(dative_entries).edit();
        editor
            .dative_bond_mut(DativeBondId(0))
            .replace_donor(position, AtomId(3));
    }

    #[rstest]
    #[case::empty(vec![], AtomPosition(1))]
    #[case::beyond_end(vec![AtomId(1)], AtomPosition(2))]
    #[should_panic]
    fn test_dative_bond_editor_view_mut_insert_donor_error(
        mut dative_entries: MoleculeEntries,
        #[case] donors: Vec<AtomId>,
        #[case] position: AtomPosition,
    ) {
        dative_entries.dative[0].0 = donors;
        let mut editor = Molecule::from_entries(dative_entries).edit();
        editor
            .dative_bond_mut(DativeBondId(0))
            .insert_donor(position, AtomId(3));
    }

    #[rstest]
    #[case::empty(vec![], AtomPosition(0))]
    #[case::end(vec![AtomId(1)], AtomPosition(1))]
    #[should_panic]
    fn test_dative_bond_editor_view_mut_remove_donor_error(
        mut dative_entries: MoleculeEntries,
        #[case] donors: Vec<AtomId>,
        #[case] position: AtomPosition,
    ) {
        dative_entries.dative[0].0 = donors;
        let mut editor = Molecule::from_entries(dative_entries).edit();
        editor
            .dative_bond_mut(DativeBondId(0))
            .remove_donor(position);
    }

    #[rstest]
    #[case::missing_atom(vec![AtomId(4)], MoleculeIntegrityError::InvalidReference { entity: Entity::Atom(AtomId(4)) })]
    #[case::repeated_donor(vec![AtomId(1), AtomId(1)], MoleculeIntegrityError::DuplicateAtom { entity: Entity::DativeBond(DativeBondId(0)), atom: AtomId(1) })]
    #[case::acceptor_as_donor(vec![AtomId(0)], MoleculeIntegrityError::DuplicateAtom { entity: Entity::DativeBond(DativeBondId(0)), atom: AtomId(0) })]
    fn test_dative_bond_editor_view_mut_replace_donors_publication(
        dative_entries: MoleculeEntries,
        #[case] donors: Vec<AtomId>,
        #[case] expected: MoleculeIntegrityError,
    ) {
        let mut editor = Molecule::from_entries(dative_entries).edit();
        editor
            .dative_bond_mut(DativeBondId(0))
            .replace_donors(&donors);
        assert_eq!(editor.try_build(), Err(expected));
    }
}
