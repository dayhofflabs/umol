//! Multicenter bond views.

use std::collections::HashSet;

use umol_graph_core::NodeId;

use super::super::constraint::{
    MulticenterBondConstraintForm, MulticenterBondConstraintKey, MulticenterBondConstraintsForm,
};
use super::super::electrons::ElectronCountsForm;
use super::super::id::{AtomId, MulticenterBondId};
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
    multicenter_bonds: &'a MulticenterBonds,
}

impl<'a> MulticenterBondViews<'a> {
    pub(crate) fn new(molecule: &'a Molecule, multicenter_bonds: &'a MulticenterBonds) -> Self {
        Self {
            molecule,
            multicenter_bonds,
        }
    }

    pub fn count(&self) -> usize {
        self.multicenter_bonds.count()
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = MulticenterBondId> {
        self.multicenter_bonds.ids()
    }

    pub fn iter(&self) -> impl ExactSizeIterator<Item = MulticenterBondView<'a>> {
        let molecule = self.molecule;
        let set = self.multicenter_bonds;
        set.ids().map(move |id| MulticenterBondView {
            id,
            attributes: set.attributes(id),
            atoms: set.atom_nodes(id),
            molecule,
        })
    }

    pub fn contains(&self, id: MulticenterBondId) -> bool {
        self.multicenter_bonds.contains(id)
    }

    pub fn get(&self, id: MulticenterBondId) -> Option<MulticenterBondView<'a>> {
        if !self.contains(id) {
            return None;
        }
        Some(MulticenterBondView {
            id,
            attributes: self.multicenter_bonds.attributes(id),
            atoms: self.multicenter_bonds.atom_nodes(id),
            molecule: self.molecule,
        })
    }

    /// Ids of multicenter bonds incident on `atom`.
    pub fn incident_ids(
        &self,
        atom: AtomId,
    ) -> impl ExactSizeIterator<Item = MulticenterBondId> + 'a {
        self.multicenter_bonds.incident_ids(atom)
    }

    /// Whether any multicenter bond is incident on `atom`.
    pub fn has_incident(&self, atom: AtomId) -> bool {
        self.multicenter_bonds.has_incident(atom)
    }

    /// Views of multicenter bonds incident on `atom`.
    pub fn incident(
        &self,
        atom: AtomId,
    ) -> impl ExactSizeIterator<Item = MulticenterBondView<'a>> + 'a {
        let molecule = self.molecule;
        let set = self.multicenter_bonds;
        self.incident_ids(atom).map(move |id| MulticenterBondView {
            id,
            attributes: set.attributes(id),
            atoms: set.atom_nodes(id),
            molecule,
        })
    }

    /// ID of the multicenter bond whose participant set equals `atoms`, if any.
    pub fn of_id(&self, atoms: impl IntoIterator<Item = AtomId>) -> Option<MulticenterBondId> {
        let atoms: Vec<AtomId> = atoms.into_iter().collect();
        self.multicenter_bonds.coincident_id(&atoms)
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
        self.multicenter_bonds
            .ids()
            .filter(|&id| {
                self.multicenter_bonds
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
    pub id: MulticenterBondId,
    atoms: &'a [NodeId],
    pub attributes: &'a MulticenterBondForm,
    molecule: &'a Molecule,
}

impl<'a> MulticenterBondView<'a> {
    #[inline]
    pub fn electrons(&self) -> &'a ElectronCountsForm {
        &self.attributes.electrons
    }

    #[inline]
    pub fn charge(&self) -> &'a NumForm {
        &self.attributes.charge
    }

    #[inline]
    pub fn unpaired_electrons(&self) -> &'a UnpairedElectronsForm {
        &self.attributes.unpaired_electrons
    }

    /// Constraint reading of this multicenter bond: the container's read API
    /// (asserted side, meanings intact) plus the keyed accessors. Mutation
    /// stays on the stored container.
    #[inline]
    pub fn constraints(&self) -> MulticenterBondConstraintsView<'a> {
        MulticenterBondConstraintsView::new(self.molecule, self.id)
    }

    pub fn atom_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + 'a {
        self.atoms.iter().map(|&n| AtomId::from(n))
    }

    pub fn atoms(&self) -> impl ExactSizeIterator<Item = AtomView<'a>> + 'a {
        let molecule = self.molecule;
        self.atoms
            .iter()
            .map(move |&n| molecule.atom(AtomId::from(n)))
    }

    /// Sum of per-atom electron contributions on this multicenter bond.
    /// `Lit(n)` when the counts are concrete; `Undetermined` otherwise.
    pub fn electron_count(&self) -> NumForm {
        match &self.attributes.electrons {
            ElectronCountsForm::Lit(counts) => NumForm::Lit(counts.iter().sum()),
            ElectronCountsForm::Undetermined => NumForm::Undetermined,
        }
    }

    pub fn atom_count(&self) -> usize {
        self.atoms.len()
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
        self.atoms
            .iter()
            .map(|&n| AtomId::from(n))
            .filter(move |a| subset.contains(a))
            .map(move |id| molecule.atom(id))
    }

    /// Is multicenter bond ground
    pub fn is_ground(&self) -> bool {
        self.attributes.is_ground()
    }

    /// Is multicenter bond undetermined
    pub fn is_undetermined(&self) -> bool {
        self.attributes.is_undetermined()
    }
}

// Derivation layer beneath the multicenter-bond facades.

/// Stored constraint container of `bond`.
pub(crate) fn multicenter_bond_asserted_constraints(
    molecule: &Molecule,
    bond: MulticenterBondId,
) -> &MulticenterBondConstraintsForm {
    &molecule.multicenter_bond(bond).attributes.constraints
}

/// Derived side of one multicenter-bond constraint key: the electron count is
/// the bond's own per-atom contribution sum — a self-projection with no
/// absence cell, so both modes agree.
pub(crate) fn multicenter_bond_derived_constraint(
    molecule: &Molecule,
    bond: MulticenterBondId,
    key: MulticenterBondConstraintKey,
    _complete: bool,
) -> Option<MulticenterBondConstraintForm> {
    match key {
        MulticenterBondConstraintKey::ElectronCount => {
            Some(MulticenterBondConstraintForm::electron_count(
                molecule.multicenter_bond(bond).electron_count(),
            ))
        }
    }
}

/// Mutable borrowed view of a multicenter bond: its id, member atoms (owned)
/// and mutable data. Molecule-scope peer of `MulticenterBondView`.
#[derive(Debug)]
pub struct MulticenterBondViewMut<'a> {
    pub id: MulticenterBondId,
    pub atoms: Vec<AtomId>,
    pub attributes: &'a mut MulticenterBondForm,
}

// Builder-scope view bundles for multicenter bonds.

pub struct MulticenterBondEditorView<'a> {
    pub id: MulticenterBondId,
    atoms: &'a [NodeId],
    pub attributes: &'a MulticenterBondForm,
}

impl<'a> MulticenterBondEditorView<'a> {
    pub(crate) fn new(
        id: MulticenterBondId,
        atoms: &'a [NodeId],
        attributes: &'a MulticenterBondForm,
    ) -> Self {
        Self {
            id,
            atoms,
            attributes,
        }
    }

    pub fn atom_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + 'a {
        self.atoms.iter().map(|&n| AtomId::from(n))
    }
}

pub struct MulticenterBondEditorViewMut<'a> {
    pub id: MulticenterBondId,
    atoms: &'a [NodeId],
    pub attributes: &'a mut MulticenterBondForm,
}

impl<'a> MulticenterBondEditorViewMut<'a> {
    pub(crate) fn new(
        id: MulticenterBondId,
        atoms: &'a [NodeId],
        attributes: &'a mut MulticenterBondForm,
    ) -> Self {
        Self {
            id,
            atoms,
            attributes,
        }
    }

    pub fn atom_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + '_ {
        self.atoms.iter().map(|&n| AtomId::from(n))
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;
    use umol_graph_core::NodeId;

    use super::super::assert_exact_size_by;
    use super::{MulticenterBondEditorView, MulticenterBondEditorViewMut};
    use crate::ir::aromatic::AromaticSystemForm;
    use crate::ir::atom::AtomForm;
    use crate::ir::bond::BondForm;
    use crate::ir::dative::DativeBondForm;
    use crate::ir::id::{AtomId, MulticenterBondId};
    use crate::ir::molecule::{Molecule, MoleculeEntries};
    use crate::ir::multicenter::MulticenterBondForm;
    use crate::ir::noncovalent::{NoncovalentBondForm, NoncovalentBondKind};
    use crate::ir::num::NumForm;

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
            |view| (view.id, view.atom_ids().collect::<Vec<_>>()),
        );
        assert_exact_size_by(
            molecule.multicenter_bonds().iter(),
            vec![(MulticenterBondId(0), vec![AtomId(0), AtomId(1), AtomId(2)])],
            |view| (view.id, view.atom_ids().collect::<Vec<_>>()),
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
            |view| view.id,
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
        assert_eq!(view.id, MulticenterBondId(0));
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
            |atom| atom.id,
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
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, expected);
    }

    #[rstest]
    fn test_multicenter_bond_editor_view_atom_ids() {
        let atoms = [NodeId(0), NodeId(1), NodeId(2)];
        let attributes = MulticenterBondForm::default();
        let view = MulticenterBondEditorView::new(MulticenterBondId(0), &atoms, &attributes);
        assert_exact_size_by(
            view.atom_ids(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
            |id| id,
        );
    }

    #[rstest]
    fn test_multicenter_bond_editor_view_mut_atom_ids() {
        let atoms = [NodeId(0), NodeId(1), NodeId(2)];
        let mut attributes = MulticenterBondForm::default();
        let view = MulticenterBondEditorViewMut::new(MulticenterBondId(0), &atoms, &mut attributes);
        assert_exact_size_by(
            view.atom_ids(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
            |id| id,
        );
    }
}
