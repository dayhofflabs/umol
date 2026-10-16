//! Aromatic system views.

use std::collections::HashSet;

use umol_graph_core::NodeId;

use super::super::aromatic::{AromaticSystemForm, AromaticSystems};
use super::super::constraint::{
    AromaticSystemConstraintForm, AromaticSystemConstraintKey, AromaticSystemConstraintsForm,
};
use super::super::correspondence::MoleculeCorrespondence;
use super::super::electrons::ElectronCountsForm;
use super::super::id::{AromaticSystemId, AtomId, AtomPosition, BondId};
use super::super::molecule::Molecule;
use super::super::num::NumForm;
use super::super::spin::UnpairedElectronsForm;
use super::super::traits::Lattice;
use super::atom::AtomView;
use super::bond::BondView;
use super::constraints::AromaticSystemConstraintsView;

/// Namespace accessor for aromatic-system views on a `Molecule`.
#[derive(Clone, Copy)]
pub struct AromaticSystemViews<'a> {
    molecule: &'a Molecule,
}

impl<'a> AromaticSystemViews<'a> {
    pub(crate) fn new(molecule: &'a Molecule) -> Self {
        Self { molecule }
    }

    pub fn count(&self) -> usize {
        self.molecule.raw_aromatic_systems().count()
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = AromaticSystemId> {
        self.molecule.raw_aromatic_systems().ids()
    }

    pub fn iter(&self) -> impl ExactSizeIterator<Item = AromaticSystemView<'a>> {
        let molecule = self.molecule;
        self.molecule
            .raw_aromatic_systems()
            .ids()
            .map(move |id| AromaticSystemView { molecule, id })
    }

    pub fn contains(&self, id: AromaticSystemId) -> bool {
        self.molecule.raw_aromatic_systems().contains(id)
    }

    pub fn get(&self, id: AromaticSystemId) -> Option<AromaticSystemView<'a>> {
        if !self.contains(id) {
            return None;
        }
        Some(AromaticSystemView {
            molecule: self.molecule,
            id,
        })
    }

    /// Ids of aromatic systems incident on `atom`.
    pub fn incident_ids(
        &self,
        atom: AtomId,
    ) -> impl ExactSizeIterator<Item = AromaticSystemId> + 'a {
        self.molecule.raw_aromatic_systems().incident_ids(atom)
    }

    /// Whether any aromatic system is incident on `atom`.
    pub fn has_incident(&self, atom: AtomId) -> bool {
        self.molecule.raw_aromatic_systems().has_incident(atom)
    }

    /// Views of aromatic systems incident on `atom`.
    pub fn incident(
        &self,
        atom: AtomId,
    ) -> impl ExactSizeIterator<Item = AromaticSystemView<'a>> + 'a {
        let molecule = self.molecule;
        self.incident_ids(atom)
            .map(move |id| AromaticSystemView { molecule, id })
    }

    /// Id of the aromatic system whose atom set equals `atoms`, if any.
    pub fn of_id(&self, atoms: impl IntoIterator<Item = AtomId>) -> Option<AromaticSystemId> {
        let atoms: Vec<AtomId> = atoms.into_iter().collect();
        self.molecule.raw_aromatic_systems().coincident_id(&atoms)
    }

    /// View of the aromatic system whose atom set equals `atoms`, if any.
    pub fn of(&self, atoms: impl IntoIterator<Item = AtomId>) -> Option<AromaticSystemView<'a>> {
        self.of_id(atoms).map(|id| {
            self.get(id).expect(
                "aromatic system id from relation set must refer to an aromatic system in this molecule",
            )
        })
    }

    /// Ids of aromatic systems whose atoms all lie in `atoms`.
    pub fn induced_ids(&self, atoms: &[AtomId]) -> Vec<AromaticSystemId> {
        let set: HashSet<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        let aromatic_systems = self.molecule.raw_aromatic_systems();
        aromatic_systems
            .ids()
            .filter(|&id| {
                aromatic_systems
                    .atom_nodes(id)
                    .iter()
                    .all(|p| set.contains(p))
            })
            .collect()
    }

    /// Views of aromatic systems whose atoms all lie in `atoms`.
    pub fn induced(&self, atoms: &[AtomId]) -> Vec<AromaticSystemView<'a>> {
        self.induced_ids(atoms)
            .into_iter()
            .map(|id| {
                self.get(id).expect(
                    "aromatic system id from relation set must refer to an aromatic system in this molecule",
                )
            })
            .collect()
    }
}

/// Borrowed view of an aromatic system: its index, the `AromaticSystemForm`,
/// and accessors for member atoms and induced ring bonds via `atoms()` and
/// `bonds()`.
#[derive(Clone, Copy, Debug)]
pub struct AromaticSystemView<'a> {
    molecule: &'a Molecule,
    id: AromaticSystemId,
}

impl<'a> AromaticSystemView<'a> {
    #[inline]
    pub fn id(&self) -> AromaticSystemId {
        self.id
    }

    #[inline]
    pub fn attributes(&self) -> &'a AromaticSystemForm {
        self.molecule.raw_aromatic_systems().attributes(self.id)
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

    /// Constraint reading of this aromatic system: the container's read API
    /// (asserted side, meanings intact) plus the keyed accessors. Mutation
    /// stays on the stored container.
    #[inline]
    pub fn constraints(&self) -> AromaticSystemConstraintsView<'a> {
        AromaticSystemConstraintsView::new(self.molecule, self.id)
    }

    #[inline]
    pub fn atom_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + 'a {
        self.molecule.raw_aromatic_systems().atoms(self.id)
    }

    pub fn atoms(&self) -> impl ExactSizeIterator<Item = AtomView<'a>> + 'a {
        let molecule = self.molecule;
        self.atom_ids().map(move |id| molecule.atom(id))
    }

    pub fn bond_ids(&self) -> impl Iterator<Item = BondId> + 'a {
        self.molecule
            .raw_graph()
            .induced_edges(self.molecule.raw_aromatic_systems().atom_nodes(self.id))
            .map(BondId::from)
    }

    pub fn bonds(&self) -> impl Iterator<Item = BondView<'a>> + 'a {
        let molecule = self.molecule;
        self.molecule
            .raw_graph()
            .induced_edges(self.molecule.raw_aromatic_systems().atom_nodes(self.id))
            .map(move |edge| molecule.bond(BondId::from(edge)))
    }

    /// The molecule subgraph induced by this system's atoms, as a sub-to-host correspondence.
    pub fn induced_subgraph(&self) -> MoleculeCorrespondence {
        self.molecule
            .induced_subgraph(&self.atom_ids().collect::<Vec<_>>())
    }

    /// Sum of per-atom electron contributions on this aromatic system.
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

    pub fn bond_count(&self) -> usize {
        self.bond_ids().count()
    }

    /// Atom views for atoms in this system that also appear in `subset`.
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

    /// Bond views for bonds in this system that also appear in `subset`.
    pub fn overlapping_bonds<'s>(
        &self,
        subset: &'s [BondId],
    ) -> impl Iterator<Item = BondView<'a>> + 's
    where
        'a: 's,
    {
        let molecule = self.molecule;
        self.molecule
            .raw_graph()
            .induced_edges(self.molecule.raw_aromatic_systems().atom_nodes(self.id))
            .map(BondId::from)
            .filter(move |b| subset.contains(b))
            .map(move |id| molecule.bond(id))
    }

    /// Is aromatic system ground
    pub fn is_ground(&self) -> bool {
        self.attributes().is_ground()
    }

    /// Is aromatic system undetermined
    pub fn is_undetermined(&self) -> bool {
        self.attributes().is_undetermined()
    }
}

/// Read-only editor access to an aromatic system.
pub struct AromaticSystemEditorView<'a> {
    aromatic_systems: &'a AromaticSystems,
    id: AromaticSystemId,
}

impl<'a> AromaticSystemEditorView<'a> {
    pub(crate) fn new(aromatic_systems: &'a AromaticSystems, id: AromaticSystemId) -> Self {
        Self {
            aromatic_systems,
            id,
        }
    }

    #[inline]
    pub fn id(&self) -> AromaticSystemId {
        self.id
    }

    #[inline]
    pub fn attributes(&self) -> &'a AromaticSystemForm {
        self.aromatic_systems.attributes(self.id)
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
        self.aromatic_systems.atoms(self.id)
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

/// Mutable attribute access to an aromatic system.
#[derive(Debug)]
pub struct AromaticSystemViewMut<'a> {
    aromatic_systems: &'a mut AromaticSystems,
    id: AromaticSystemId,
}

impl<'a> AromaticSystemViewMut<'a> {
    pub(crate) fn new(aromatic_systems: &'a mut AromaticSystems, id: AromaticSystemId) -> Self {
        Self {
            aromatic_systems,
            id,
        }
    }

    #[inline]
    pub fn id(&self) -> AromaticSystemId {
        self.id
    }

    #[inline]
    pub fn attributes(&self) -> &AromaticSystemForm {
        self.aromatic_systems.attributes(self.id)
    }

    #[inline]
    pub fn attributes_mut(&mut self) -> &mut AromaticSystemForm {
        self.aromatic_systems.attributes_mut(self.id)
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
    pub fn constraints(&self) -> &AromaticSystemConstraintsForm {
        &self.attributes().constraints
    }

    #[inline]
    pub fn atom_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + '_ {
        self.aromatic_systems.atoms(self.id)
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

/// Mutable editor access to an aromatic system.
///
/// Structural mutations update incidence and preserve attributes and constraints.
/// Atom existence, distinctness, and overlap between systems are checked at publication.
#[derive(Debug)]
pub struct AromaticSystemEditorViewMut<'a> {
    aromatic_systems: &'a mut AromaticSystems,
    id: AromaticSystemId,
}

impl<'a> AromaticSystemEditorViewMut<'a> {
    pub(crate) fn new(aromatic_systems: &'a mut AromaticSystems, id: AromaticSystemId) -> Self {
        Self {
            aromatic_systems,
            id,
        }
    }

    #[inline]
    pub fn id(&self) -> AromaticSystemId {
        self.id
    }

    #[inline]
    pub fn attributes(&self) -> &AromaticSystemForm {
        self.aromatic_systems.attributes(self.id)
    }

    #[inline]
    pub fn attributes_mut(&mut self) -> &mut AromaticSystemForm {
        self.aromatic_systems.attributes_mut(self.id)
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
    pub fn constraints(&self) -> &AromaticSystemConstraintsForm {
        &self.attributes().constraints
    }

    #[inline]
    pub fn atom_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + '_ {
        self.aromatic_systems.atoms(self.id)
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
        self.aromatic_systems.replace_atoms(self.id, atoms);
    }

    /// Replace the atom at `position` without changing the atom count.
    ///
    /// # Panics
    ///
    /// Panics when `position` is outside the atom list.
    pub fn replace_atom(&mut self, position: AtomPosition, atom: AtomId) {
        self.aromatic_systems
            .replace_atom(self.id, position.index(), atom);
    }

    /// Insert an atom at `position`, preserving the order of the existing atoms.
    ///
    /// # Panics
    ///
    /// Panics when `position` exceeds the atom count. Insertion at the end is allowed.
    pub fn insert_atom(&mut self, position: AtomPosition, atom: AtomId) {
        self.aromatic_systems
            .insert_atom(self.id, position.index(), atom);
    }

    /// Remove the atom at `position`, preserving the order of the remaining atoms.
    ///
    /// # Panics
    ///
    /// Panics when `position` is outside the atom list.
    pub fn remove_atom(&mut self, position: AtomPosition) {
        self.aromatic_systems.remove_atom(self.id, position.index());
    }
}

/// Stored constraint container of aromatic system `id`.
pub(crate) fn aromatic_system_asserted_constraints(
    molecule: &Molecule,
    id: AromaticSystemId,
) -> &AromaticSystemConstraintsForm {
    &molecule.aromatic_system(id).attributes().constraints
}

/// Derived side of one aromatic-system constraint key: the electron count is
/// the system's own per-atom contribution sum — a self-projection with no
/// absence cell, so both modes agree.
pub(crate) fn aromatic_system_derived_constraint(
    molecule: &Molecule,
    id: AromaticSystemId,
    key: AromaticSystemConstraintKey,
    _complete: bool,
) -> Option<AromaticSystemConstraintForm> {
    match key {
        AromaticSystemConstraintKey::ElectronCount => {
            Some(AromaticSystemConstraintForm::electron_count(
                molecule.aromatic_system(id).electron_count(),
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
    use crate::ir::constraint::AromaticSystemConstraintForm;
    use crate::ir::dative::DativeBondForm;
    use crate::ir::electrons::ElectronCountsForm;
    use crate::ir::entity::Entity;
    use crate::ir::id::{AromaticSystemId, AtomId, BondId};
    use crate::ir::molecule::{Molecule, MoleculeEntries};
    use crate::ir::multicenter::MulticenterBondForm;
    use crate::ir::noncovalent::{NoncovalentBondForm, NoncovalentBondKind};
    use crate::ir::num::NumForm;
    use crate::ir::spin::UnpairedElectronsForm;
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
    fn aromatic_molecule() -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::default(); 6],
            aromatic: vec![
                (
                    vec![AtomId(2), AtomId(0), AtomId(1)],
                    AromaticSystemForm {
                        charge: NumForm::Lit(1),
                        ..Default::default()
                    },
                ),
                (
                    vec![AtomId(5), AtomId(3), AtomId(4)],
                    AromaticSystemForm {
                        charge: NumForm::Lit(-1),
                        ..Default::default()
                    },
                ),
            ],
            ..Default::default()
        })
    }

    #[fixture]
    fn aromatic_entries() -> MoleculeEntries {
        MoleculeEntries {
            atoms: vec![AtomForm::default(); 5],
            aromatic: vec![
                (
                    vec![AtomId(0), AtomId(2), AtomId(1)],
                    AromaticSystemForm::from_electrons(vec![2, 1, 1])
                        .with_charge(-1_i64)
                        .with_constraint(AromaticSystemConstraintForm::electron_count(4)),
                ),
                (vec![AtomId(4)], AromaticSystemForm::default()),
            ],
            ..Default::default()
        }
    }

    #[rstest]
    fn test_aromatic_system_views_count(molecule: Molecule) {
        assert_eq!(molecule.aromatic_systems().count(), 1);
    }

    #[rstest]
    fn test_aromatic_system_views_ids(molecule: Molecule) {
        assert_exact_size_by(Molecule::default().aromatic_systems().ids(), vec![], |id| {
            id
        });
        assert_exact_size_by(
            molecule.aromatic_systems().ids(),
            vec![AromaticSystemId(0)],
            |id| id,
        );
    }

    #[rstest]
    fn test_aromatic_system_views_iter(molecule: Molecule) {
        assert_exact_size_by(
            Molecule::default().aromatic_systems().iter(),
            vec![],
            |view| (view.id, view.atom_ids().collect::<Vec<_>>()),
        );
        assert_exact_size_by(
            molecule.aromatic_systems().iter(),
            vec![(AromaticSystemId(0), vec![AtomId(0), AtomId(1), AtomId(2)])],
            |view| (view.id, view.atom_ids().collect::<Vec<_>>()),
        );
    }

    #[rstest]
    #[case::participant(AtomId(0), vec![AromaticSystemId(0)])]
    #[case::uninvolved(AtomId(3), vec![])]
    fn test_aromatic_system_views_incident(
        molecule: Molecule,
        #[case] atom: AtomId,
        #[case] expected: Vec<AromaticSystemId>,
    ) {
        assert_exact_size_by(
            molecule.aromatic_systems().incident_ids(atom),
            expected.clone(),
            |id| id,
        );
        assert_exact_size_by(
            molecule.aromatic_systems().incident(atom),
            expected,
            |view| view.id,
        );
    }

    #[rstest]
    #[case::present(AromaticSystemId(0), true)]
    #[case::absent(AromaticSystemId(99), false)]
    fn test_aromatic_system_views_contains(
        molecule: Molecule,
        #[case] id: AromaticSystemId,
        #[case] expected: bool,
    ) {
        assert_eq!(molecule.aromatic_systems().contains(id), expected);
    }

    #[rstest]
    fn test_aromatic_system_views_get(molecule: Molecule) {
        let res = molecule.aromatic_systems().get(AromaticSystemId(0));
        assert!(res.is_some());
        let view = res.unwrap();
        assert_eq!(view.id, AromaticSystemId(0));
        assert_eq!(
            view.atom_ids().collect::<Vec<_>>(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
        );
    }

    #[rstest]
    fn test_aromatic_system_views_get_none(molecule: Molecule) {
        let res = molecule.aromatic_systems().get(AromaticSystemId(99));
        assert!(res.is_none());
    }

    #[rstest]
    #[case(AromaticSystemId(0))]
    #[case(AromaticSystemId(1))]
    fn test_aromatic_system_view_id(
        #[from(aromatic_molecule)] molecule: Molecule,
        #[case] id: AromaticSystemId,
    ) {
        assert_eq!(molecule.aromatic_system(id).id(), id);
    }

    #[rstest]
    #[case(AromaticSystemId(0), AromaticSystemForm { charge: NumForm::Lit(1), ..Default::default() })]
    #[case(AromaticSystemId(1), AromaticSystemForm { charge: NumForm::Lit(-1), ..Default::default() })]
    fn test_aromatic_system_view_attributes(
        #[from(aromatic_molecule)] molecule: Molecule,
        #[case] id: AromaticSystemId,
        #[case] expected: AromaticSystemForm,
    ) {
        let attributes = {
            let view = molecule.aromatic_system(id);
            view.attributes()
        };
        assert_eq!(attributes, &expected);
    }

    #[rstest]
    #[case(AromaticSystemId(0), vec![AtomId(2), AtomId(0), AtomId(1)])]
    #[case(AromaticSystemId(1), vec![AtomId(5), AtomId(3), AtomId(4)])]
    fn test_aromatic_system_view_atom_ids_order(
        #[from(aromatic_molecule)] molecule: Molecule,
        #[case] id: AromaticSystemId,
        #[case] expected: Vec<AtomId>,
    ) {
        let atom_ids = {
            let view = molecule.aromatic_system(id);
            view.atom_ids()
        };
        assert_exact_size_by(atom_ids, expected, |id| id);
    }

    #[rstest]
    fn test_aromatic_system_view_atom_ids(molecule: Molecule) {
        assert_exact_size_by(
            molecule.aromatic_system(AromaticSystemId(0)).atom_ids(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
            |id| id,
        );
    }

    #[rstest]
    fn test_aromatic_system_view_atoms(molecule: Molecule) {
        assert_exact_size_by(
            molecule.aromatic_system(AromaticSystemId(0)).atoms(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
            |atom| atom.id(),
        );
    }

    #[rstest]
    fn test_aromatic_system_view_bond_ids(molecule: Molecule) {
        assert_eq!(
            molecule
                .aromatic_system(AromaticSystemId(0))
                .bond_ids()
                .collect::<Vec<_>>(),
            vec![BondId(0), BondId(1)],
        );
    }

    #[rstest]
    fn test_aromatic_system_view_bonds(molecule: Molecule) {
        let ids: Vec<BondId> = molecule
            .aromatic_system(AromaticSystemId(0))
            .bonds()
            .map(|v| v.id())
            .collect();
        assert_eq!(ids, vec![BondId(0), BondId(1)]);
    }

    #[rstest]
    fn test_aromatic_system_view_induced_subgraph(molecule: Molecule) {
        let correspondence = molecule
            .aromatic_system(AromaticSystemId(0))
            .induced_subgraph();
        assert_eq!(
            correspondence.atoms().matched_pairs(),
            &[
                (AtomId(0), AtomId(0)),
                (AtomId(1), AtomId(1)),
                (AtomId(2), AtomId(2)),
            ],
        );
        assert_eq!(
            correspondence.bonds().matched_pairs(),
            &[(BondId(0), BondId(0)), (BondId(1), BondId(1))],
        );
    }

    #[rstest]
    fn test_aromatic_system_view_electron_count(molecule: Molecule) {
        assert_eq!(
            molecule
                .aromatic_system(AromaticSystemId(0))
                .electron_count(),
            NumForm::Undetermined,
        );
    }

    #[rstest]
    fn test_aromatic_system_view_atom_count(molecule: Molecule) {
        assert_eq!(
            molecule.aromatic_system(AromaticSystemId(0)).atom_count(),
            3
        );
    }

    #[rstest]
    fn test_aromatic_system_view_bond_count(molecule: Molecule) {
        assert_eq!(
            molecule.aromatic_system(AromaticSystemId(0)).bond_count(),
            2
        );
    }

    #[rstest]
    #[case::two_in(vec![AtomId(0), AtomId(1)], vec![AtomId(0), AtomId(1)])]
    #[case::all_in(vec![AtomId(0), AtomId(1), AtomId(2)], vec![AtomId(0), AtomId(1), AtomId(2)])]
    #[case::disjoint(vec![AtomId(3)], vec![])]
    fn test_aromatic_system_view_overlapping_atoms(
        molecule: Molecule,
        #[case] subset: Vec<AtomId>,
        #[case] expected: Vec<AtomId>,
    ) {
        let ids: Vec<AtomId> = molecule
            .aromatic_system(AromaticSystemId(0))
            .overlapping_atoms(&subset)
            .map(|v| v.id())
            .collect();
        assert_eq!(ids, expected);
    }

    #[rstest]
    #[case::one(vec![BondId(0)], vec![BondId(0)])]
    #[case::both(vec![BondId(0), BondId(1)], vec![BondId(0), BondId(1)])]
    #[case::other(vec![BondId(2)], vec![])]
    fn test_aromatic_system_view_overlapping_bonds(
        molecule: Molecule,
        #[case] subset: Vec<BondId>,
        #[case] expected: Vec<BondId>,
    ) {
        let ids: Vec<BondId> = molecule
            .aromatic_system(AromaticSystemId(0))
            .overlapping_bonds(&subset)
            .map(|v| v.id())
            .collect();
        assert_eq!(ids, expected);
    }

    #[rstest]
    #[case(AromaticSystemId(0))]
    #[case(AromaticSystemId(1))]
    fn test_aromatic_system_editor_view_id(
        #[from(aromatic_molecule)] molecule: Molecule,
        #[case] id: AromaticSystemId,
    ) {
        let editor = molecule.edit();
        assert_eq!(editor.aromatic_system(id).id(), id);
    }

    #[rstest]
    #[case(AromaticSystemId(0), AromaticSystemForm { charge: NumForm::Lit(1), ..Default::default() })]
    #[case(AromaticSystemId(1), AromaticSystemForm { charge: NumForm::Lit(-1), ..Default::default() })]
    fn test_aromatic_system_editor_view_attributes(
        #[from(aromatic_molecule)] molecule: Molecule,
        #[case] id: AromaticSystemId,
        #[case] expected: AromaticSystemForm,
    ) {
        let editor = molecule.edit();
        let attributes = {
            let view = editor.aromatic_system(id);
            view.attributes()
        };
        assert_eq!(attributes, &expected);
    }

    #[rstest]
    #[case::literal(ElectronCountsForm::Lit(vec![2, 1, 1]), NumForm::Lit(4))]
    #[case::undetermined(ElectronCountsForm::Undetermined, NumForm::Undetermined)]
    fn test_aromatic_system_editor_view_electron_count(
        mut aromatic_entries: MoleculeEntries,
        #[case] electrons: ElectronCountsForm,
        #[case] expected: NumForm,
    ) {
        aromatic_entries.aromatic[0].1.electrons = electrons.clone();
        let attributes = aromatic_entries.aromatic[0].1.clone();
        let editor = Molecule::from_entries(aromatic_entries).edit();
        let fields = {
            let view = editor.aromatic_system(AromaticSystemId(0));
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
    #[case(AromaticSystemId(0), vec![AtomId(2), AtomId(0), AtomId(1)])]
    #[case(AromaticSystemId(1), vec![AtomId(5), AtomId(3), AtomId(4)])]
    fn test_aromatic_system_editor_view_atom_ids_order(
        #[from(aromatic_molecule)] molecule: Molecule,
        #[case] id: AromaticSystemId,
        #[case] expected: Vec<AtomId>,
    ) {
        let editor = molecule.edit();
        let atom_ids = {
            let view = editor.aromatic_system(id);
            view.atom_ids()
        };
        assert_exact_size_by(atom_ids, expected, |id| id);
    }

    #[rstest]
    fn test_aromatic_system_editor_view_atom_ids(molecule: Molecule) {
        let editor = molecule.edit();
        let view = editor.aromatic_system(AromaticSystemId(0));
        assert_exact_size_by(
            view.atom_ids(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
            |id| id,
        );
    }

    #[rstest]
    #[case(AromaticSystemId(0))]
    #[case(AromaticSystemId(1))]
    fn test_aromatic_system_view_mut_attributes_mut(
        #[from(aromatic_molecule)] mut molecule: Molecule,
        #[case] id: AromaticSystemId,
    ) {
        let expected = AromaticSystemForm {
            electrons: ElectronCountsForm::Lit(vec![2, 1, 1]),
            charge: NumForm::Lit(-2),
            unpaired_electrons: UnpairedElectronsForm {
                count: NumForm::Lit(1),
                multiplicity: NumForm::Lit(2),
            },
            ..Default::default()
        };
        {
            let mut view = molecule.aromatic_system_mut(id);
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
        assert_eq!(molecule.aromatic_system(id).attributes(), &expected);
    }

    #[rstest]
    #[case(AromaticSystemId(0), vec![AtomId(2), AtomId(0), AtomId(1)])]
    #[case(AromaticSystemId(1), vec![AtomId(5), AtomId(3), AtomId(4)])]
    fn test_aromatic_system_view_mut_atom_ids_order(
        #[from(aromatic_molecule)] mut molecule: Molecule,
        #[case] id: AromaticSystemId,
        #[case] expected: Vec<AtomId>,
    ) {
        let view = molecule.aromatic_system_mut(id);
        assert_exact_size_by(view.atom_ids(), expected, |id| id);
    }

    #[rstest]
    #[case(AromaticSystemId(0))]
    #[case(AromaticSystemId(1))]
    fn test_aromatic_system_editor_view_mut_attributes_mut(
        #[from(aromatic_molecule)] molecule: Molecule,
        #[case] id: AromaticSystemId,
    ) {
        let mut editor = molecule.edit();
        let expected = AromaticSystemForm {
            electrons: ElectronCountsForm::Lit(vec![2, 1, 1]),
            charge: NumForm::Lit(-2),
            unpaired_electrons: UnpairedElectronsForm {
                count: NumForm::Lit(1),
                multiplicity: NumForm::Lit(2),
            },
            ..Default::default()
        };
        {
            let mut view = editor.aromatic_system_mut(id);
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
        assert_eq!(editor.aromatic_system(id).attributes(), &expected);
    }

    #[rstest]
    #[case(AromaticSystemId(0), vec![AtomId(2), AtomId(0), AtomId(1)])]
    #[case(AromaticSystemId(1), vec![AtomId(5), AtomId(3), AtomId(4)])]
    fn test_aromatic_system_editor_view_mut_atom_ids_order(
        #[from(aromatic_molecule)] molecule: Molecule,
        #[case] id: AromaticSystemId,
        #[case] expected: Vec<AtomId>,
    ) {
        let mut editor = molecule.edit();
        let view = editor.aromatic_system_mut(id);
        assert_exact_size_by(view.atom_ids(), expected, |id| id);
    }

    #[rstest]
    fn test_aromatic_system_editor_view_mut_atom_ids() {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::default(); 3],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemForm::default(),
            )],
            ..Default::default()
        });
        let mut editor = molecule.edit();
        let view = editor.aromatic_system_mut(AromaticSystemId(0));
        assert_exact_size_by(
            view.atom_ids(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
            |id| id,
        );
    }

    #[rstest]
    #[case::empty(vec![])]
    #[case::single(vec![AtomId(3)])]
    #[case::reordered(vec![AtomId(1), AtomId(0), AtomId(2)])]
    #[case::expanded(vec![AtomId(3), AtomId(0), AtomId(2), AtomId(1)])]
    fn test_aromatic_system_editor_view_mut_replace_atoms(
        mut aromatic_entries: MoleculeEntries,
        #[case] atoms: Vec<AtomId>,
    ) {
        let mut editor = Molecule::from_entries(aromatic_entries.clone()).edit();
        {
            let mut view = editor.aromatic_system_mut(AromaticSystemId(0));
            view.replace_atoms(&atoms);
            assert_eq!(view.atom_ids().collect::<Vec<_>>(), atoms);
            assert_eq!(view.atom_count(), atoms.len());
            assert_eq!(view.electron_count(), NumForm::Lit(4));
        }
        aromatic_entries.aromatic[0].0 = atoms;
        assert_eq!(
            editor.finish(),
            Ok(Molecule::from_entries(aromatic_entries))
        );
    }

    #[rstest]
    #[case::missing_atom(vec![AtomId(5)], MoleculeIntegrityError::InvalidReference { entity: Entity::Atom(AtomId(5)) })]
    #[case::repeated_atom(vec![AtomId(0), AtomId(0)], MoleculeIntegrityError::DuplicateAtom { entity: Entity::AromaticSystem(AromaticSystemId(0)), atom: AtomId(0) })]
    #[case::overlap(vec![AtomId(4)], MoleculeIntegrityError::AromaticSystemsOverlap { atom: AtomId(4) })]
    fn test_aromatic_system_editor_view_mut_replace_atoms_publication(
        aromatic_entries: MoleculeEntries,
        #[case] atoms: Vec<AtomId>,
        #[case] expected: MoleculeIntegrityError,
    ) {
        let mut editor = Molecule::from_entries(aromatic_entries).edit();
        editor
            .aromatic_system_mut(AromaticSystemId(0))
            .replace_atoms(&atoms);
        assert_eq!(editor.finish(), Err(expected));
    }

    #[rstest]
    #[case::attributes_first(true)]
    #[case::atoms_first(false)]
    fn test_aromatic_system_editor_view_mut_replace_atoms_attributes(
        mut aromatic_entries: MoleculeEntries,
        #[case] attributes_first: bool,
    ) {
        let atoms = [AtomId(3), AtomId(0)];
        let electrons = ElectronCountsForm::Lit(vec![1, 2]);
        let mut editor = Molecule::from_entries(aromatic_entries.clone()).edit();
        {
            let mut view = editor.aromatic_system_mut(AromaticSystemId(0));
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
        aromatic_entries.aromatic[0].0 = atoms.to_vec();
        aromatic_entries.aromatic[0].1.electrons = electrons;
        assert_eq!(
            editor.finish(),
            Ok(Molecule::from_entries(aromatic_entries))
        );
    }

    #[rstest]
    #[case::first(AtomPosition(0), vec![AtomId(3), AtomId(2), AtomId(1)])]
    #[case::last(AtomPosition(2), vec![AtomId(0), AtomId(2), AtomId(3)])]
    fn test_aromatic_system_editor_view_mut_replace_atom(
        mut aromatic_entries: MoleculeEntries,
        #[case] position: AtomPosition,
        #[case] expected: Vec<AtomId>,
    ) {
        let mut editor = Molecule::from_entries(aromatic_entries.clone()).edit();
        editor
            .aromatic_system_mut(AromaticSystemId(0))
            .replace_atom(position, AtomId(3));
        assert_eq!(
            editor
                .aromatic_system(AromaticSystemId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            expected
        );
        aromatic_entries.aromatic[0].0 = expected;
        assert_eq!(
            editor.finish(),
            Ok(Molecule::from_entries(aromatic_entries))
        );
    }

    #[rstest]
    #[case::empty(vec![], AtomPosition(0))]
    #[case::end(vec![AtomId(0)], AtomPosition(1))]
    #[should_panic]
    fn test_aromatic_system_editor_view_mut_replace_atom_error(
        mut aromatic_entries: MoleculeEntries,
        #[case] atoms: Vec<AtomId>,
        #[case] position: AtomPosition,
    ) {
        aromatic_entries.aromatic[0].0 = atoms;
        let mut editor = Molecule::from_entries(aromatic_entries).edit();
        editor
            .aromatic_system_mut(AromaticSystemId(0))
            .replace_atom(position, AtomId(3));
    }

    #[rstest]
    #[case::first(vec![AtomId(0), AtomId(2)], AtomPosition(0), vec![AtomId(3), AtomId(0), AtomId(2)])]
    #[case::middle(vec![AtomId(0), AtomId(2)], AtomPosition(1), vec![AtomId(0), AtomId(3), AtomId(2)])]
    #[case::end(vec![AtomId(0), AtomId(2)], AtomPosition(2), vec![AtomId(0), AtomId(2), AtomId(3)])]
    #[case::empty(vec![], AtomPosition(0), vec![AtomId(3)])]
    fn test_aromatic_system_editor_view_mut_insert_atom(
        mut aromatic_entries: MoleculeEntries,
        #[case] atoms: Vec<AtomId>,
        #[case] position: AtomPosition,
        #[case] expected: Vec<AtomId>,
    ) {
        aromatic_entries.aromatic[0].0 = atoms;
        let mut editor = Molecule::from_entries(aromatic_entries.clone()).edit();
        editor
            .aromatic_system_mut(AromaticSystemId(0))
            .insert_atom(position, AtomId(3));
        assert_eq!(
            editor
                .aromatic_system(AromaticSystemId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            expected
        );
        aromatic_entries.aromatic[0].0 = expected;
        assert_eq!(
            editor.finish(),
            Ok(Molecule::from_entries(aromatic_entries))
        );
    }

    #[rstest]
    #[case::empty(vec![], AtomPosition(1))]
    #[case::beyond_end(vec![AtomId(0)], AtomPosition(2))]
    #[should_panic]
    fn test_aromatic_system_editor_view_mut_insert_atom_error(
        mut aromatic_entries: MoleculeEntries,
        #[case] atoms: Vec<AtomId>,
        #[case] position: AtomPosition,
    ) {
        aromatic_entries.aromatic[0].0 = atoms;
        let mut editor = Molecule::from_entries(aromatic_entries).edit();
        editor
            .aromatic_system_mut(AromaticSystemId(0))
            .insert_atom(position, AtomId(3));
    }

    #[rstest]
    #[case::first(vec![AtomId(0), AtomId(2)], AtomPosition(0), vec![AtomId(2)])]
    #[case::last(vec![AtomId(0), AtomId(2)], AtomPosition(1), vec![AtomId(0)])]
    #[case::only(vec![AtomId(0)], AtomPosition(0), vec![])]
    fn test_aromatic_system_editor_view_mut_remove_atom(
        mut aromatic_entries: MoleculeEntries,
        #[case] atoms: Vec<AtomId>,
        #[case] position: AtomPosition,
        #[case] expected: Vec<AtomId>,
    ) {
        aromatic_entries.aromatic[0].0 = atoms;
        let mut editor = Molecule::from_entries(aromatic_entries.clone()).edit();
        editor
            .aromatic_system_mut(AromaticSystemId(0))
            .remove_atom(position);
        assert_eq!(
            editor
                .aromatic_system(AromaticSystemId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            expected
        );
        aromatic_entries.aromatic[0].0 = expected;
        assert_eq!(
            editor.finish(),
            Ok(Molecule::from_entries(aromatic_entries))
        );
    }

    #[rstest]
    #[case::empty(vec![], AtomPosition(0))]
    #[case::end(vec![AtomId(0)], AtomPosition(1))]
    #[should_panic]
    fn test_aromatic_system_editor_view_mut_remove_atom_error(
        mut aromatic_entries: MoleculeEntries,
        #[case] atoms: Vec<AtomId>,
        #[case] position: AtomPosition,
    ) {
        aromatic_entries.aromatic[0].0 = atoms;
        let mut editor = Molecule::from_entries(aromatic_entries).edit();
        editor
            .aromatic_system_mut(AromaticSystemId(0))
            .remove_atom(position);
    }

    #[rstest]
    fn test_aromatic_system_editor_view_mut_replace_atoms_incidence(
        aromatic_entries: MoleculeEntries,
    ) {
        let mut editor = Molecule::from_entries(aromatic_entries).edit();
        {
            let mut view = editor.aromatic_system_mut(AromaticSystemId(0));
            view.replace_atoms(&[AtomId(3), AtomId(0)]);
            view.replace_atom(AtomPosition(1), AtomId(2));
            view.insert_atom(AtomPosition(2), AtomId(1));
            view.remove_atom(AtomPosition(0));
        }
        let molecule = editor.finish().unwrap();
        let incidence: Vec<Vec<_>> = molecule
            .atoms()
            .ids()
            .map(|atom| molecule.aromatic_systems().incident_ids(atom).collect())
            .collect();
        assert_eq!(
            incidence,
            vec![
                vec![],
                vec![AromaticSystemId(0)],
                vec![AromaticSystemId(0)],
                vec![],
                vec![AromaticSystemId(1)]
            ]
        );
    }
}
