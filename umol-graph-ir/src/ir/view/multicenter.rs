//! Multicenter bond views.

use std::collections::HashSet;

use umol_graph_core::NodeId;

use super::super::constraint::{
    MulticenterBondConstraintForm, MulticenterBondConstraintKey, MulticenterBondConstraintsForm,
};
use super::super::electrons::ElectronCountsForm;
use super::super::id::{AtomId, AtomPosition, MulticenterBondId};
use super::super::molecule::Molecule;
use super::super::multicenter::{MulticenterBondForm, MulticenterBonds};
use super::super::num::NumForm;
use super::super::spin::UnpairedElectronsForm;
use super::super::traits::Lattice;
use super::atom::AtomView;
use super::constraints::MulticenterBondConstraintsView;

/// Namespace accessor for multicenter-bond views on a `Molecule`.
#[derive(Clone, Copy)]
pub struct MulticenterBondViews<'a> {
    molecule: &'a Molecule,
}

impl<'a> MulticenterBondViews<'a> {
    pub(crate) fn new(molecule: &'a Molecule) -> Self {
        Self { molecule }
    }

    pub fn count(&self) -> usize {
        self.molecule.raw_multicenter_bonds().count()
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = MulticenterBondId> {
        self.molecule.raw_multicenter_bonds().ids()
    }

    pub fn iter(&self) -> impl ExactSizeIterator<Item = MulticenterBondView<'a>> {
        let molecule = self.molecule;
        self.molecule
            .raw_multicenter_bonds()
            .ids()
            .map(move |id| MulticenterBondView { molecule, id })
    }

    pub fn contains(&self, id: MulticenterBondId) -> bool {
        self.molecule.raw_multicenter_bonds().contains(id)
    }

    pub fn get(&self, id: MulticenterBondId) -> Option<MulticenterBondView<'a>> {
        if !self.contains(id) {
            return None;
        }
        Some(MulticenterBondView {
            molecule: self.molecule,
            id,
        })
    }

    /// Ids of multicenter bonds incident on `atom`.
    pub fn incident_ids(
        &self,
        atom: AtomId,
    ) -> impl ExactSizeIterator<Item = MulticenterBondId> + 'a {
        self.molecule.raw_multicenter_bonds().incident_ids(atom)
    }

    /// Whether any multicenter bond is incident on `atom`.
    pub fn has_incident(&self, atom: AtomId) -> bool {
        self.molecule.raw_multicenter_bonds().has_incident(atom)
    }

    /// Views of multicenter bonds incident on `atom`.
    pub fn incident(
        &self,
        atom: AtomId,
    ) -> impl ExactSizeIterator<Item = MulticenterBondView<'a>> + 'a {
        let molecule = self.molecule;
        self.incident_ids(atom)
            .map(move |id| MulticenterBondView { molecule, id })
    }

    /// ID of the multicenter bond whose participant set equals `atoms`, if any.
    pub fn of_id(&self, atoms: impl IntoIterator<Item = AtomId>) -> Option<MulticenterBondId> {
        let atoms: Vec<AtomId> = atoms.into_iter().collect();
        self.molecule.raw_multicenter_bonds().coincident_id(&atoms)
    }

    /// View of the multicenter bond whose participant set equals `atoms`, if any.
    pub fn of(&self, atoms: impl IntoIterator<Item = AtomId>) -> Option<MulticenterBondView<'a>> {
        self.of_id(atoms).map(|id| {
            self.get(id).expect(
                "multicenter bond id from relation set must refer to a multicenter bond in this molecule",
            )
        })
    }

    /// Ids of multicenter bonds whose participants all lie in `atoms`.
    pub fn induced_ids(&self, atoms: &[AtomId]) -> Vec<MulticenterBondId> {
        let set: HashSet<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        let multicenter_bonds = self.molecule.raw_multicenter_bonds();
        multicenter_bonds
            .ids()
            .filter(|&id| {
                multicenter_bonds
                    .atom_nodes(id)
                    .iter()
                    .all(|p| set.contains(p))
            })
            .collect()
    }

    /// Views of multicenter bonds whose participants all lie in `atoms`.
    pub fn induced(&self, atoms: &[AtomId]) -> Vec<MulticenterBondView<'a>> {
        self.induced_ids(atoms)
            .into_iter()
            .map(|id| {
                self.get(id).expect(
                    "multicenter bond id from relation set must refer to a multicenter bond in this molecule",
                )
            })
            .collect()
    }
}

/// Borrowed view of a multicenter bond: its index, member atoms via
/// `atoms()`, and underlying `MulticenterBondForm`.
#[derive(Clone, Copy, Debug)]
pub struct MulticenterBondView<'a> {
    molecule: &'a Molecule,
    id: MulticenterBondId,
}

impl<'a> MulticenterBondView<'a> {
    #[inline]
    pub fn id(&self) -> MulticenterBondId {
        self.id
    }

    #[inline]
    pub fn attributes(&self) -> &'a MulticenterBondForm {
        self.molecule.raw_multicenter_bonds().attributes(self.id)
    }

    #[inline]
    pub fn electrons(&self) -> &'a ElectronCountsForm {
        &self.attributes().electrons
    }

    #[inline]
    pub fn charge(&self) -> &'a NumForm {
        &self.attributes().charge
    }

    #[inline]
    pub fn unpaired_electrons(&self) -> &'a UnpairedElectronsForm {
        &self.attributes().unpaired_electrons
    }

    /// Constraint reading of this multicenter bond: the container's read API
    /// (asserted side, meanings intact) plus the keyed accessors. Mutation
    /// stays on the stored container.
    #[inline]
    pub fn constraints(&self) -> MulticenterBondConstraintsView<'a> {
        MulticenterBondConstraintsView::new(self.molecule, self.id)
    }

    #[inline]
    pub fn atom_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + 'a {
        self.molecule.raw_multicenter_bonds().atoms(self.id)
    }

    pub fn atoms(&self) -> impl ExactSizeIterator<Item = AtomView<'a>> + 'a {
        let molecule = self.molecule;
        self.atom_ids().map(move |id| molecule.atom(id))
    }

    /// Sum of per-atom electron contributions on this multicenter bond.
    /// `Lit(n)` when the counts are concrete; `Undetermined` otherwise.
    /// Includes every stored count, including counts beyond the atom list.
    pub fn electron_count(&self) -> NumForm {
        match &self.attributes().electrons {
            ElectronCountsForm::Lit(counts) => NumForm::Lit(counts.iter().sum()),
            ElectronCountsForm::Undetermined => NumForm::Undetermined,
        }
    }

    #[inline]
    pub fn atom_count(&self) -> usize {
        self.atom_ids().len()
    }

    /// Atom views for atoms in this multicenter bond that also appear in `subset`.
    pub fn overlapping_atoms<'s>(
        &self,
        subset: &'s [AtomId],
    ) -> impl Iterator<Item = AtomView<'a>> + 's
    where
        'a: 's,
    {
        let molecule = self.molecule;
        self.atom_ids()
            .filter(move |a| subset.contains(a))
            .map(move |id| molecule.atom(id))
    }

    /// Is multicenter bond ground
    pub fn is_ground(&self) -> bool {
        self.attributes().is_ground()
    }

    /// Is multicenter bond undetermined
    pub fn is_undetermined(&self) -> bool {
        self.attributes().is_undetermined()
    }
}

/// Read-only editor access to a multicenter bond.
pub struct MulticenterBondEditorView<'a> {
    multicenter_bonds: &'a MulticenterBonds,
    id: MulticenterBondId,
}

impl<'a> MulticenterBondEditorView<'a> {
    pub(crate) fn new(multicenter_bonds: &'a MulticenterBonds, id: MulticenterBondId) -> Self {
        Self {
            multicenter_bonds,
            id,
        }
    }

    #[inline]
    pub fn id(&self) -> MulticenterBondId {
        self.id
    }

    #[inline]
    pub fn attributes(&self) -> &'a MulticenterBondForm {
        self.multicenter_bonds.attributes(self.id)
    }

    #[inline]
    pub fn electrons(&self) -> &'a ElectronCountsForm {
        &self.attributes().electrons
    }

    #[inline]
    pub fn charge(&self) -> &'a NumForm {
        &self.attributes().charge
    }

    #[inline]
    pub fn unpaired_electrons(&self) -> &'a UnpairedElectronsForm {
        &self.attributes().unpaired_electrons
    }

    #[inline]
    pub fn atom_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + 'a {
        self.multicenter_bonds.atoms(self.id)
    }

    /// Sum of all stored literal electron contributions, or `Undetermined`.
    pub fn electron_count(&self) -> NumForm {
        match self.electrons() {
            ElectronCountsForm::Lit(counts) => NumForm::Lit(counts.iter().sum()),
            ElectronCountsForm::Undetermined => NumForm::Undetermined,
        }
    }

    #[inline]
    pub fn atom_count(&self) -> usize {
        self.atom_ids().len()
    }
}

/// Mutable attribute access to a multicenter bond.
#[derive(Debug)]
pub struct MulticenterBondViewMut<'a> {
    multicenter_bonds: &'a mut MulticenterBonds,
    id: MulticenterBondId,
}

impl<'a> MulticenterBondViewMut<'a> {
    pub(crate) fn new(multicenter_bonds: &'a mut MulticenterBonds, id: MulticenterBondId) -> Self {
        Self {
            multicenter_bonds,
            id,
        }
    }

    #[inline]
    pub fn id(&self) -> MulticenterBondId {
        self.id
    }

    #[inline]
    pub fn attributes(&self) -> &MulticenterBondForm {
        self.multicenter_bonds.attributes(self.id)
    }

    #[inline]
    pub fn attributes_mut(&mut self) -> &mut MulticenterBondForm {
        self.multicenter_bonds.attributes_mut(self.id)
    }

    #[inline]
    pub fn electrons(&self) -> &ElectronCountsForm {
        &self.attributes().electrons
    }

    #[inline]
    pub fn charge(&self) -> &NumForm {
        &self.attributes().charge
    }

    #[inline]
    pub fn unpaired_electrons(&self) -> &UnpairedElectronsForm {
        &self.attributes().unpaired_electrons
    }

    #[inline]
    pub fn constraints(&self) -> &MulticenterBondConstraintsForm {
        &self.attributes().constraints
    }

    #[inline]
    pub fn atom_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + '_ {
        self.multicenter_bonds.atoms(self.id)
    }

    /// Sum of all stored literal electron contributions, or `Undetermined`.
    pub fn electron_count(&self) -> NumForm {
        match self.electrons() {
            ElectronCountsForm::Lit(counts) => NumForm::Lit(counts.iter().sum()),
            ElectronCountsForm::Undetermined => NumForm::Undetermined,
        }
    }

    #[inline]
    pub fn atom_count(&self) -> usize {
        self.atom_ids().len()
    }
}

/// Mutable editor access to a multicenter bond.
///
/// Structural mutations update incidence and preserve attributes and constraints.
/// Atom existence, distinctness, and uniqueness of bond atom sets are checked at publication.
#[derive(Debug)]
pub struct MulticenterBondEditorViewMut<'a> {
    multicenter_bonds: &'a mut MulticenterBonds,
    id: MulticenterBondId,
}

impl<'a> MulticenterBondEditorViewMut<'a> {
    pub(crate) fn new(multicenter_bonds: &'a mut MulticenterBonds, id: MulticenterBondId) -> Self {
        Self {
            multicenter_bonds,
            id,
        }
    }

    #[inline]
    pub fn id(&self) -> MulticenterBondId {
        self.id
    }

    #[inline]
    pub fn attributes(&self) -> &MulticenterBondForm {
        self.multicenter_bonds.attributes(self.id)
    }

    #[inline]
    pub fn attributes_mut(&mut self) -> &mut MulticenterBondForm {
        self.multicenter_bonds.attributes_mut(self.id)
    }

    #[inline]
    pub fn electrons(&self) -> &ElectronCountsForm {
        &self.attributes().electrons
    }

    #[inline]
    pub fn charge(&self) -> &NumForm {
        &self.attributes().charge
    }

    #[inline]
    pub fn unpaired_electrons(&self) -> &UnpairedElectronsForm {
        &self.attributes().unpaired_electrons
    }

    #[inline]
    pub fn constraints(&self) -> &MulticenterBondConstraintsForm {
        &self.attributes().constraints
    }

    #[inline]
    pub fn atom_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + '_ {
        self.multicenter_bonds.atoms(self.id)
    }

    /// Sum of all stored literal electron contributions, or `Undetermined`.
    pub fn electron_count(&self) -> NumForm {
        match self.electrons() {
            ElectronCountsForm::Lit(counts) => NumForm::Lit(counts.iter().sum()),
            ElectronCountsForm::Undetermined => NumForm::Undetermined,
        }
    }

    #[inline]
    pub fn atom_count(&self) -> usize {
        self.atom_ids().len()
    }

    /// Replace the atom list, preserving the supplied order.
    pub fn replace_atoms(&mut self, atoms: &[AtomId]) {
        self.multicenter_bonds.replace_atoms(self.id, atoms);
    }

    /// Replace the atom at `position` without changing the atom count.
    ///
    /// # Panics
    ///
    /// Panics when `position` is outside the atom list.
    pub fn replace_atom(&mut self, position: AtomPosition, atom: AtomId) {
        self.multicenter_bonds
            .replace_atom(self.id, position.index(), atom);
    }

    /// Insert an atom at `position`, preserving the order of the existing atoms.
    ///
    /// # Panics
    ///
    /// Panics when `position` exceeds the atom count. Insertion at the end is allowed.
    pub fn insert_atom(&mut self, position: AtomPosition, atom: AtomId) {
        self.multicenter_bonds
            .insert_atom(self.id, position.index(), atom);
    }

    /// Remove the atom at `position`, preserving the order of the remaining atoms.
    ///
    /// # Panics
    ///
    /// Panics when `position` is outside the atom list.
    pub fn remove_atom(&mut self, position: AtomPosition) {
        self.multicenter_bonds
            .remove_atom(self.id, position.index());
    }
}

/// Stored constraint container of `bond`.
pub(crate) fn multicenter_bond_asserted_constraints(
    molecule: &Molecule,
    id: MulticenterBondId,
) -> &MulticenterBondConstraintsForm {
    &molecule.multicenter_bond(id).attributes().constraints
}

/// Derived side of one multicenter-bond constraint key: the electron count is
/// the bond's own per-atom contribution sum — a self-projection with no
/// absence cell, so both modes agree.
pub(crate) fn multicenter_bond_derived_constraint(
    molecule: &Molecule,
    id: MulticenterBondId,
    key: MulticenterBondConstraintKey,
    _complete: bool,
) -> Option<MulticenterBondConstraintForm> {
    match key {
        MulticenterBondConstraintKey::ElectronCount => {
            Some(MulticenterBondConstraintForm::electron_count(
                molecule.multicenter_bond(id).electron_count(),
            ))
        }
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
    use crate::ir::constraint::MulticenterBondConstraintForm;
    use crate::ir::dative::DativeBondForm;
    use crate::ir::electrons::ElectronCountsForm;
    use crate::ir::entity::Entity;
    use crate::ir::id::{AtomId, AtomPosition, MulticenterBondId};
    use crate::ir::molecule::{Molecule, MoleculeEntries, MoleculeIntegrityError};
    use crate::ir::multicenter::MulticenterBondForm;
    use crate::ir::noncovalent::{NoncovalentBondForm, NoncovalentBondKind};
    use crate::ir::num::NumForm;
    use crate::ir::spin::UnpairedElectronsForm;

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
    fn multicenter_molecule() -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::default(); 6],
            multicenter: vec![
                (
                    vec![AtomId(2), AtomId(0), AtomId(1)],
                    MulticenterBondForm {
                        charge: NumForm::Lit(1),
                        ..Default::default()
                    },
                ),
                (
                    vec![AtomId(5), AtomId(3), AtomId(4)],
                    MulticenterBondForm {
                        charge: NumForm::Lit(-1),
                        ..Default::default()
                    },
                ),
            ],
            ..Default::default()
        })
    }

    #[fixture]
    fn multicenter_entries() -> MoleculeEntries {
        MoleculeEntries {
            atoms: vec![AtomForm::default(); 5],
            multicenter: vec![
                (
                    vec![AtomId(0), AtomId(2), AtomId(1)],
                    MulticenterBondForm::from_electrons(vec![2, 1, 1])
                        .with_charge(-1_i64)
                        .with_constraint(MulticenterBondConstraintForm::electron_count(4)),
                ),
                (vec![AtomId(4)], MulticenterBondForm::default()),
            ],
            ..Default::default()
        }
    }

    #[rstest]
    fn test_multicenter_bond_views_count(molecule: Molecule) {
        assert_eq!(molecule.multicenter_bonds().count(), 1);
    }

    #[rstest]
    fn test_multicenter_bond_views_ids(molecule: Molecule) {
        assert_exact_size_by(
            Molecule::default().multicenter_bonds().ids(),
            vec![],
            |id| id,
        );
        assert_exact_size_by(
            molecule.multicenter_bonds().ids(),
            vec![MulticenterBondId(0)],
            |id| id,
        );
    }

    #[rstest]
    fn test_multicenter_bond_views_iter(molecule: Molecule) {
        assert_exact_size_by(
            Molecule::default().multicenter_bonds().iter(),
            vec![],
            |view| (view.id(), view.atom_ids().collect::<Vec<_>>()),
        );
        assert_exact_size_by(
            molecule.multicenter_bonds().iter(),
            vec![(MulticenterBondId(0), vec![AtomId(0), AtomId(1), AtomId(2)])],
            |view| (view.id(), view.atom_ids().collect::<Vec<_>>()),
        );
    }

    #[rstest]
    #[case::participant(AtomId(0), vec![MulticenterBondId(0)])]
    #[case::uninvolved(AtomId(3), vec![])]
    fn test_multicenter_bond_views_incident(
        molecule: Molecule,
        #[case] atom: AtomId,
        #[case] expected: Vec<MulticenterBondId>,
    ) {
        assert_exact_size_by(
            molecule.multicenter_bonds().incident_ids(atom),
            expected.clone(),
            |id| id,
        );
        assert_exact_size_by(
            molecule.multicenter_bonds().incident(atom),
            expected,
            |view| view.id(),
        );
    }

    #[rstest]
    #[case::present(MulticenterBondId(0), true)]
    #[case::absent(MulticenterBondId(99), false)]
    fn test_multicenter_bond_views_contains(
        molecule: Molecule,
        #[case] id: MulticenterBondId,
        #[case] expected: bool,
    ) {
        assert_eq!(molecule.multicenter_bonds().contains(id), expected);
    }

    #[rstest]
    fn test_multicenter_bond_views_get(molecule: Molecule) {
        let res = molecule.multicenter_bonds().get(MulticenterBondId(0));
        assert!(res.is_some());
        let view = res.unwrap();
        assert_eq!(view.id(), MulticenterBondId(0));
        assert_eq!(
            view.atom_ids().collect::<Vec<_>>(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
        );
    }

    #[rstest]
    fn test_multicenter_bond_views_get_none(molecule: Molecule) {
        let res = molecule.multicenter_bonds().get(MulticenterBondId(99));
        assert!(res.is_none());
    }

    #[rstest]
    #[case(MulticenterBondId(0))]
    #[case(MulticenterBondId(1))]
    fn test_multicenter_bond_view_id(
        #[from(multicenter_molecule)] molecule: Molecule,
        #[case] id: MulticenterBondId,
    ) {
        assert_eq!(molecule.multicenter_bond(id).id(), id);
    }

    #[rstest]
    #[case(MulticenterBondId(0), MulticenterBondForm { charge: NumForm::Lit(1), ..Default::default() })]
    #[case(MulticenterBondId(1), MulticenterBondForm { charge: NumForm::Lit(-1), ..Default::default() })]
    fn test_multicenter_bond_view_attributes(
        #[from(multicenter_molecule)] molecule: Molecule,
        #[case] id: MulticenterBondId,
        #[case] expected: MulticenterBondForm,
    ) {
        let attributes = {
            let view = molecule.multicenter_bond(id);
            view.attributes()
        };
        assert_eq!(attributes, &expected);
    }

    #[rstest]
    #[case(MulticenterBondId(0), vec![AtomId(2), AtomId(0), AtomId(1)])]
    #[case(MulticenterBondId(1), vec![AtomId(5), AtomId(3), AtomId(4)])]
    fn test_multicenter_bond_view_atom_ids_order(
        #[from(multicenter_molecule)] molecule: Molecule,
        #[case] id: MulticenterBondId,
        #[case] expected: Vec<AtomId>,
    ) {
        let atom_ids = {
            let view = molecule.multicenter_bond(id);
            view.atom_ids()
        };
        assert_exact_size_by(atom_ids, expected, |id| id);
    }

    #[rstest]
    fn test_multicenter_bond_view_atom_ids(molecule: Molecule) {
        assert_exact_size_by(
            molecule.multicenter_bond(MulticenterBondId(0)).atom_ids(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
            |id| id,
        );
    }

    #[rstest]
    fn test_multicenter_bond_view_atoms(molecule: Molecule) {
        assert_exact_size_by(
            molecule.multicenter_bond(MulticenterBondId(0)).atoms(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
            |atom| atom.id(),
        );
    }

    #[rstest]
    fn test_multicenter_bond_view_electron_count(molecule: Molecule) {
        assert_eq!(
            molecule
                .multicenter_bond(MulticenterBondId(0))
                .electron_count(),
            NumForm::Undetermined,
        );
    }

    #[rstest]
    fn test_multicenter_bond_view_atom_count(molecule: Molecule) {
        assert_eq!(
            molecule.multicenter_bond(MulticenterBondId(0)).atom_count(),
            3,
        );
    }

    #[rstest]
    #[case::two_in(vec![AtomId(0), AtomId(1)], vec![AtomId(0), AtomId(1)])]
    #[case::all_in(vec![AtomId(0), AtomId(1), AtomId(2)], vec![AtomId(0), AtomId(1), AtomId(2)])]
    #[case::disjoint(vec![AtomId(3)], vec![])]
    fn test_multicenter_bond_view_overlapping_atoms(
        molecule: Molecule,
        #[case] subset: Vec<AtomId>,
        #[case] expected: Vec<AtomId>,
    ) {
        let ids: Vec<AtomId> = molecule
            .multicenter_bond(MulticenterBondId(0))
            .overlapping_atoms(&subset)
            .map(|v| v.id())
            .collect();
        assert_eq!(ids, expected);
    }

    #[rstest]
    #[case(MulticenterBondId(0))]
    #[case(MulticenterBondId(1))]
    fn test_multicenter_bond_editor_view_id(
        #[from(multicenter_molecule)] molecule: Molecule,
        #[case] id: MulticenterBondId,
    ) {
        let editor = molecule.edit();
        assert_eq!(editor.multicenter_bond(id).id(), id);
    }

    #[rstest]
    #[case(MulticenterBondId(0), MulticenterBondForm { charge: NumForm::Lit(1), ..Default::default() })]
    #[case(MulticenterBondId(1), MulticenterBondForm { charge: NumForm::Lit(-1), ..Default::default() })]
    fn test_multicenter_bond_editor_view_attributes(
        #[from(multicenter_molecule)] molecule: Molecule,
        #[case] id: MulticenterBondId,
        #[case] expected: MulticenterBondForm,
    ) {
        let editor = molecule.edit();
        let attributes = {
            let view = editor.multicenter_bond(id);
            view.attributes()
        };
        assert_eq!(attributes, &expected);
    }

    #[rstest]
    #[case::literal(ElectronCountsForm::Lit(vec![2, 1, 1]), NumForm::Lit(4))]
    #[case::undetermined(ElectronCountsForm::Undetermined, NumForm::Undetermined)]
    fn test_multicenter_bond_editor_view_electron_count(
        mut multicenter_entries: MoleculeEntries,
        #[case] electrons: ElectronCountsForm,
        #[case] expected: NumForm,
    ) {
        multicenter_entries.multicenter[0].1.electrons = electrons.clone();
        let attributes = multicenter_entries.multicenter[0].1.clone();
        let editor = Molecule::from_entries(multicenter_entries).edit();
        let fields = {
            let view = editor.multicenter_bond(MulticenterBondId(0));
            assert_eq!(view.electron_count(), expected);
            assert_eq!(view.atom_count(), 3);
            (view.electrons(), view.charge(), view.unpaired_electrons())
        };
        assert_eq!(
            fields,
            (
                &electrons,
                &attributes.charge,
                &attributes.unpaired_electrons
            )
        );
    }

    #[rstest]
    #[case(MulticenterBondId(0), vec![AtomId(2), AtomId(0), AtomId(1)])]
    #[case(MulticenterBondId(1), vec![AtomId(5), AtomId(3), AtomId(4)])]
    fn test_multicenter_bond_editor_view_atom_ids_order(
        #[from(multicenter_molecule)] molecule: Molecule,
        #[case] id: MulticenterBondId,
        #[case] expected: Vec<AtomId>,
    ) {
        let editor = molecule.edit();
        let atom_ids = {
            let view = editor.multicenter_bond(id);
            view.atom_ids()
        };
        assert_exact_size_by(atom_ids, expected, |id| id);
    }

    #[rstest]
    fn test_multicenter_bond_editor_view_atom_ids(molecule: Molecule) {
        let editor = molecule.edit();
        let view = editor.multicenter_bond(MulticenterBondId(0));
        assert_exact_size_by(
            view.atom_ids(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
            |id| id,
        );
    }

    #[rstest]
    #[case(MulticenterBondId(0))]
    #[case(MulticenterBondId(1))]
    fn test_multicenter_bond_view_mut_attributes_mut(
        #[from(multicenter_molecule)] mut molecule: Molecule,
        #[case] id: MulticenterBondId,
    ) {
        let expected = MulticenterBondForm {
            electrons: ElectronCountsForm::Lit(vec![2, 1, 1]),
            charge: NumForm::Lit(-2),
            unpaired_electrons: UnpairedElectronsForm {
                count: NumForm::Lit(1),
                multiplicity: NumForm::Lit(2),
            },
            ..Default::default()
        };
        {
            let mut view = molecule.multicenter_bond_mut(id);
            assert_eq!(view.id(), id);
            assert_eq!(view.electron_count(), NumForm::Undetermined);
            assert_eq!(view.atom_count(), 3);
            *view.attributes_mut() = expected.clone();
            assert_eq!(view.attributes(), &expected);
            assert_eq!(view.electrons(), &expected.electrons);
            assert_eq!(view.charge(), &expected.charge);
            assert_eq!(view.unpaired_electrons(), &expected.unpaired_electrons);
            assert_eq!(view.electron_count(), NumForm::Lit(4));
        }
        assert_eq!(molecule.multicenter_bond(id).attributes(), &expected);
    }

    #[rstest]
    #[case(MulticenterBondId(0), vec![AtomId(2), AtomId(0), AtomId(1)])]
    #[case(MulticenterBondId(1), vec![AtomId(5), AtomId(3), AtomId(4)])]
    fn test_multicenter_bond_view_mut_atom_ids_order(
        #[from(multicenter_molecule)] mut molecule: Molecule,
        #[case] id: MulticenterBondId,
        #[case] expected: Vec<AtomId>,
    ) {
        let view = molecule.multicenter_bond_mut(id);
        assert_exact_size_by(view.atom_ids(), expected, |id| id);
    }

    #[rstest]
    #[case(MulticenterBondId(0))]
    #[case(MulticenterBondId(1))]
    fn test_multicenter_bond_editor_view_mut_attributes_mut(
        #[from(multicenter_molecule)] molecule: Molecule,
        #[case] id: MulticenterBondId,
    ) {
        let mut editor = molecule.edit();
        let expected = MulticenterBondForm {
            electrons: ElectronCountsForm::Lit(vec![2, 1, 1]),
            charge: NumForm::Lit(-2),
            unpaired_electrons: UnpairedElectronsForm {
                count: NumForm::Lit(1),
                multiplicity: NumForm::Lit(2),
            },
            ..Default::default()
        };
        {
            let mut view = editor.multicenter_bond_mut(id);
            assert_eq!(view.id(), id);
            assert_eq!(view.electron_count(), NumForm::Undetermined);
            assert_eq!(view.atom_count(), 3);
            *view.attributes_mut() = expected.clone();
            assert_eq!(view.attributes(), &expected);
            assert_eq!(view.electrons(), &expected.electrons);
            assert_eq!(view.charge(), &expected.charge);
            assert_eq!(view.unpaired_electrons(), &expected.unpaired_electrons);
            assert_eq!(view.electron_count(), NumForm::Lit(4));
        }
        assert_eq!(editor.multicenter_bond(id).attributes(), &expected);
    }

    #[rstest]
    #[case(MulticenterBondId(0), vec![AtomId(2), AtomId(0), AtomId(1)])]
    #[case(MulticenterBondId(1), vec![AtomId(5), AtomId(3), AtomId(4)])]
    fn test_multicenter_bond_editor_view_mut_atom_ids_order(
        #[from(multicenter_molecule)] molecule: Molecule,
        #[case] id: MulticenterBondId,
        #[case] expected: Vec<AtomId>,
    ) {
        let mut editor = molecule.edit();
        let view = editor.multicenter_bond_mut(id);
        assert_exact_size_by(view.atom_ids(), expected, |id| id);
    }

    #[rstest]
    fn test_multicenter_bond_editor_view_mut_atom_ids() {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::default(); 3],
            multicenter: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                MulticenterBondForm::default(),
            )],
            ..Default::default()
        });
        let mut editor = molecule.edit();
        let view = editor.multicenter_bond_mut(MulticenterBondId(0));
        assert_exact_size_by(
            view.atom_ids(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
            |id| id,
        );
    }

    #[rstest]
    #[case::empty(vec![])]
    #[case::single(vec![AtomId(3)])]
    #[case::overlapping(vec![AtomId(4), AtomId(0)])]
    #[case::reordered(vec![AtomId(1), AtomId(0), AtomId(2)])]
    #[case::expanded(vec![AtomId(3), AtomId(0), AtomId(2), AtomId(1)])]
    fn test_multicenter_bond_editor_view_mut_replace_atoms(
        mut multicenter_entries: MoleculeEntries,
        #[case] atoms: Vec<AtomId>,
    ) {
        let mut editor = Molecule::from_entries(multicenter_entries.clone()).edit();
        {
            let mut view = editor.multicenter_bond_mut(MulticenterBondId(0));
            view.replace_atoms(&atoms);
            assert_eq!(view.atom_ids().collect::<Vec<_>>(), atoms);
            assert_eq!(view.atom_count(), atoms.len());
            assert_eq!(view.electron_count(), NumForm::Lit(4));
        }
        multicenter_entries.multicenter[0].0 = atoms;
        assert_eq!(
            editor.finish(),
            Ok(Molecule::from_entries(multicenter_entries))
        );
    }

    #[rstest]
    #[case::missing_atom(vec![AtomId(5)], MoleculeIntegrityError::InvalidReference { entity: Entity::Atom(AtomId(5)) })]
    #[case::repeated_atom(vec![AtomId(0), AtomId(0)], MoleculeIntegrityError::DuplicateAtom { entity: Entity::MulticenterBond(MulticenterBondId(0)), atom: AtomId(0) })]
    #[case::identical(vec![AtomId(4)], MoleculeIntegrityError::IdenticalMulticenterBonds { atoms: vec![AtomId(4)] })]
    fn test_multicenter_bond_editor_view_mut_replace_atoms_publication(
        multicenter_entries: MoleculeEntries,
        #[case] atoms: Vec<AtomId>,
        #[case] expected: MoleculeIntegrityError,
    ) {
        let mut editor = Molecule::from_entries(multicenter_entries).edit();
        editor
            .multicenter_bond_mut(MulticenterBondId(0))
            .replace_atoms(&atoms);
        assert_eq!(editor.finish(), Err(expected));
    }

    #[rstest]
    #[case::attributes_first(true)]
    #[case::atoms_first(false)]
    fn test_multicenter_bond_editor_view_mut_replace_atoms_attributes(
        mut multicenter_entries: MoleculeEntries,
        #[case] attributes_first: bool,
    ) {
        let atoms = [AtomId(3), AtomId(0)];
        let electrons = ElectronCountsForm::Lit(vec![1, 2]);
        let mut editor = Molecule::from_entries(multicenter_entries.clone()).edit();
        {
            let mut view = editor.multicenter_bond_mut(MulticenterBondId(0));
            if attributes_first {
                view.attributes_mut().electrons = electrons.clone();
                view.replace_atoms(&atoms);
            } else {
                view.replace_atoms(&atoms);
                view.attributes_mut().electrons = electrons.clone();
            }
            assert_eq!(view.atom_ids().collect::<Vec<_>>(), atoms);
            assert_eq!(view.electron_count(), NumForm::Lit(3));
        }
        multicenter_entries.multicenter[0].0 = atoms.to_vec();
        multicenter_entries.multicenter[0].1.electrons = electrons;
        assert_eq!(
            editor.finish(),
            Ok(Molecule::from_entries(multicenter_entries))
        );
    }

    #[rstest]
    #[case::first(AtomPosition(0), vec![AtomId(3), AtomId(2), AtomId(1)])]
    #[case::last(AtomPosition(2), vec![AtomId(0), AtomId(2), AtomId(3)])]
    fn test_multicenter_bond_editor_view_mut_replace_atom(
        mut multicenter_entries: MoleculeEntries,
        #[case] position: AtomPosition,
        #[case] expected: Vec<AtomId>,
    ) {
        let mut editor = Molecule::from_entries(multicenter_entries.clone()).edit();
        editor
            .multicenter_bond_mut(MulticenterBondId(0))
            .replace_atom(position, AtomId(3));
        assert_eq!(
            editor
                .multicenter_bond(MulticenterBondId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            expected
        );
        multicenter_entries.multicenter[0].0 = expected;
        assert_eq!(
            editor.finish(),
            Ok(Molecule::from_entries(multicenter_entries))
        );
    }

    #[rstest]
    #[case::empty(vec![], AtomPosition(0))]
    #[case::end(vec![AtomId(0)], AtomPosition(1))]
    #[should_panic]
    fn test_multicenter_bond_editor_view_mut_replace_atom_error(
        mut multicenter_entries: MoleculeEntries,
        #[case] atoms: Vec<AtomId>,
        #[case] position: AtomPosition,
    ) {
        multicenter_entries.multicenter[0].0 = atoms;
        let mut editor = Molecule::from_entries(multicenter_entries).edit();
        editor
            .multicenter_bond_mut(MulticenterBondId(0))
            .replace_atom(position, AtomId(3));
    }

    #[rstest]
    #[case::first(vec![AtomId(0), AtomId(2)], AtomPosition(0), vec![AtomId(3), AtomId(0), AtomId(2)])]
    #[case::middle(vec![AtomId(0), AtomId(2)], AtomPosition(1), vec![AtomId(0), AtomId(3), AtomId(2)])]
    #[case::end(vec![AtomId(0), AtomId(2)], AtomPosition(2), vec![AtomId(0), AtomId(2), AtomId(3)])]
    #[case::empty(vec![], AtomPosition(0), vec![AtomId(3)])]
    fn test_multicenter_bond_editor_view_mut_insert_atom(
        mut multicenter_entries: MoleculeEntries,
        #[case] atoms: Vec<AtomId>,
        #[case] position: AtomPosition,
        #[case] expected: Vec<AtomId>,
    ) {
        multicenter_entries.multicenter[0].0 = atoms;
        let mut editor = Molecule::from_entries(multicenter_entries.clone()).edit();
        editor
            .multicenter_bond_mut(MulticenterBondId(0))
            .insert_atom(position, AtomId(3));
        assert_eq!(
            editor
                .multicenter_bond(MulticenterBondId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            expected
        );
        multicenter_entries.multicenter[0].0 = expected;
        assert_eq!(
            editor.finish(),
            Ok(Molecule::from_entries(multicenter_entries))
        );
    }

    #[rstest]
    #[case::empty(vec![], AtomPosition(1))]
    #[case::beyond_end(vec![AtomId(0)], AtomPosition(2))]
    #[should_panic]
    fn test_multicenter_bond_editor_view_mut_insert_atom_error(
        mut multicenter_entries: MoleculeEntries,
        #[case] atoms: Vec<AtomId>,
        #[case] position: AtomPosition,
    ) {
        multicenter_entries.multicenter[0].0 = atoms;
        let mut editor = Molecule::from_entries(multicenter_entries).edit();
        editor
            .multicenter_bond_mut(MulticenterBondId(0))
            .insert_atom(position, AtomId(3));
    }

    #[rstest]
    #[case::first(vec![AtomId(0), AtomId(2)], AtomPosition(0), vec![AtomId(2)])]
    #[case::last(vec![AtomId(0), AtomId(2)], AtomPosition(1), vec![AtomId(0)])]
    #[case::only(vec![AtomId(0)], AtomPosition(0), vec![])]
    fn test_multicenter_bond_editor_view_mut_remove_atom(
        mut multicenter_entries: MoleculeEntries,
        #[case] atoms: Vec<AtomId>,
        #[case] position: AtomPosition,
        #[case] expected: Vec<AtomId>,
    ) {
        multicenter_entries.multicenter[0].0 = atoms;
        let mut editor = Molecule::from_entries(multicenter_entries.clone()).edit();
        editor
            .multicenter_bond_mut(MulticenterBondId(0))
            .remove_atom(position);
        assert_eq!(
            editor
                .multicenter_bond(MulticenterBondId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            expected
        );
        multicenter_entries.multicenter[0].0 = expected;
        assert_eq!(
            editor.finish(),
            Ok(Molecule::from_entries(multicenter_entries))
        );
    }

    #[rstest]
    #[case::empty(vec![], AtomPosition(0))]
    #[case::end(vec![AtomId(0)], AtomPosition(1))]
    #[should_panic]
    fn test_multicenter_bond_editor_view_mut_remove_atom_error(
        mut multicenter_entries: MoleculeEntries,
        #[case] atoms: Vec<AtomId>,
        #[case] position: AtomPosition,
    ) {
        multicenter_entries.multicenter[0].0 = atoms;
        let mut editor = Molecule::from_entries(multicenter_entries).edit();
        editor
            .multicenter_bond_mut(MulticenterBondId(0))
            .remove_atom(position);
    }

    #[rstest]
    fn test_multicenter_bond_editor_view_mut_replace_atoms_incidence(
        multicenter_entries: MoleculeEntries,
    ) {
        let mut editor = Molecule::from_entries(multicenter_entries).edit();
        {
            let mut view = editor.multicenter_bond_mut(MulticenterBondId(0));
            view.replace_atoms(&[AtomId(3), AtomId(0)]);
            view.replace_atom(AtomPosition(1), AtomId(2));
            view.insert_atom(AtomPosition(2), AtomId(1));
            view.remove_atom(AtomPosition(0));
        }
        let molecule = editor.finish().unwrap();
        let incidence: Vec<Vec<_>> = molecule
            .atoms()
            .ids()
            .map(|atom| molecule.multicenter_bonds().incident_ids(atom).collect())
            .collect();
        assert_eq!(
            incidence,
            vec![
                vec![],
                vec![MulticenterBondId(0)],
                vec![MulticenterBondId(0)],
                vec![],
                vec![MulticenterBondId(1)]
            ]
        );
    }
}
