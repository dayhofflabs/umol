//! Noncovalent bond views.

use std::collections::HashSet;

use super::super::constraint::{
    NoncovalentBondConstraintForm, NoncovalentBondConstraintKey, NoncovalentBondConstraintsForm,
};
use super::super::id::{AtomId, AtomPosition, NoncovalentBondId};
use super::super::molecule::Molecule;
use super::super::noncovalent::{NoncovalentBondForm, NoncovalentBondKindForm, NoncovalentBonds};
use super::super::traits::Lattice;
use super::atom::AtomView;
use super::constraints::NoncovalentBondConstraintsView;

/// Namespace accessor for noncovalent-bond views on a `Molecule`.
#[derive(Clone, Copy)]
pub struct NoncovalentBondViews<'a> {
    molecule: &'a Molecule,
}

impl<'a> NoncovalentBondViews<'a> {
    pub(crate) fn new(molecule: &'a Molecule) -> Self {
        Self { molecule }
    }

    pub fn count(&self) -> usize {
        self.molecule.raw_noncovalent_bonds().count()
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = NoncovalentBondId> {
        self.molecule.raw_noncovalent_bonds().ids()
    }

    pub fn iter(&self) -> impl ExactSizeIterator<Item = NoncovalentBondView<'a>> {
        let molecule = self.molecule;
        self.molecule
            .raw_noncovalent_bonds()
            .ids()
            .map(move |id| NoncovalentBondView { molecule, id })
    }

    pub fn contains(&self, id: NoncovalentBondId) -> bool {
        self.molecule.raw_noncovalent_bonds().contains(id)
    }

    pub fn get(&self, id: NoncovalentBondId) -> Option<NoncovalentBondView<'a>> {
        if !self.contains(id) {
            return None;
        }
        Some(NoncovalentBondView {
            molecule: self.molecule,
            id,
        })
    }

    /// Ids of noncovalent bonds incident on `atom`.
    pub fn incident_ids(
        &self,
        atom: AtomId,
    ) -> impl ExactSizeIterator<Item = NoncovalentBondId> + 'a {
        self.molecule.raw_noncovalent_bonds().incident_ids(atom)
    }

    /// Whether any noncovalent bond is incident on `atom`.
    pub fn has_incident(&self, atom: AtomId) -> bool {
        self.molecule.raw_noncovalent_bonds().has_incident(atom)
    }

    /// Views of noncovalent bonds incident on `atom`.
    pub fn incident(
        &self,
        atom: AtomId,
    ) -> impl ExactSizeIterator<Item = NoncovalentBondView<'a>> + 'a {
        let molecule = self.molecule;
        self.incident_ids(atom)
            .map(move |id| NoncovalentBondView { molecule, id })
    }

    /// Id of the noncovalent bond between `a` and `b`, if any.
    pub fn of_id(&self, first: AtomId, second: AtomId) -> Option<NoncovalentBondId> {
        self.molecule
            .raw_noncovalent_bonds()
            .coincident_id(first, second)
    }

    /// View of the noncovalent bond between `first` and `second`, if any.
    pub fn of(&self, first: AtomId, second: AtomId) -> Option<NoncovalentBondView<'a>> {
        self.of_id(first, second).map(|id| {
            self.get(id).expect(
                "noncovalent bond id from relation set must refer to a noncovalent bond in this molecule",
            )
        })
    }

    /// Ids of noncovalent bonds whose endpoints both lie in `atoms`.
    pub fn induced_ids(&self, atoms: &[AtomId]) -> Vec<NoncovalentBondId> {
        let set: HashSet<AtomId> = atoms.iter().copied().collect();
        let noncovalent_bonds = self.molecule.raw_noncovalent_bonds();
        noncovalent_bonds
            .ids()
            .filter(|&id| {
                noncovalent_bonds
                    .atoms(id)
                    .iter()
                    .all(|atom| set.contains(atom))
            })
            .collect()
    }

    /// Views of noncovalent bonds whose endpoints both lie in `atoms`.
    pub fn induced(&self, atoms: &[AtomId]) -> Vec<NoncovalentBondView<'a>> {
        self.induced_ids(atoms)
            .into_iter()
            .map(|id| {
                self.get(id).expect(
                    "noncovalent bond id from relation set must refer to a noncovalent bond in this molecule",
                )
            })
            .collect()
    }
}

/// Borrowed view of a noncovalent bond: the two participating atoms plus data.
#[derive(Clone, Copy, Debug)]
pub struct NoncovalentBondView<'a> {
    molecule: &'a Molecule,
    id: NoncovalentBondId,
}

impl<'a> NoncovalentBondView<'a> {
    #[inline]
    pub fn id(&self) -> NoncovalentBondId {
        self.id
    }

    #[inline]
    pub fn attributes(&self) -> &'a NoncovalentBondForm {
        self.molecule.raw_noncovalent_bonds().attributes(self.id)
    }

    #[inline]
    pub fn kind(&self) -> &'a NoncovalentBondKindForm {
        &self.attributes().kind
    }

    /// Constraint reading of this noncovalent bond: the container's read API
    /// (asserted side, meanings intact) plus the keyed accessors. Mutation
    /// stays on the stored container.
    #[inline]
    pub fn constraints(&self) -> NoncovalentBondConstraintsView<'a> {
        NoncovalentBondConstraintsView::new(self.molecule, self.id)
    }

    /// The two atom ids in this noncovalent interaction.
    #[inline]
    pub fn atom_ids(&self) -> [AtomId; 2] {
        self.molecule.raw_noncovalent_bonds().atoms(self.id)
    }

    /// Views of the two atoms in this noncovalent interaction.
    pub fn atoms(&self) -> [AtomView<'a>; 2] {
        let [a, b] = self.atom_ids();
        [self.molecule.atom(a), self.molecule.atom(b)]
    }

    /// Is noncovalent bond ground
    pub fn is_ground(&self) -> bool {
        self.attributes().is_ground()
    }

    /// Is noncovalent bond undetermined
    pub fn is_undetermined(&self) -> bool {
        self.attributes().is_undetermined()
    }
}

/// Read-only editor access to a noncovalent bond.
pub struct NoncovalentBondEditorView<'a> {
    noncovalent_bonds: &'a NoncovalentBonds,
    id: NoncovalentBondId,
}

impl<'a> NoncovalentBondEditorView<'a> {
    pub(crate) fn new(noncovalent_bonds: &'a NoncovalentBonds, id: NoncovalentBondId) -> Self {
        Self {
            noncovalent_bonds,
            id,
        }
    }

    #[inline]
    pub fn id(&self) -> NoncovalentBondId {
        self.id
    }

    #[inline]
    pub fn attributes(&self) -> &'a NoncovalentBondForm {
        self.noncovalent_bonds.attributes(self.id)
    }

    #[inline]
    pub fn kind(&self) -> &'a NoncovalentBondKindForm {
        &self.attributes().kind
    }

    #[inline]
    pub fn atom_ids(&self) -> [AtomId; 2] {
        self.noncovalent_bonds.atoms(self.id)
    }
}

/// Mutable attribute access to a noncovalent bond.
#[derive(Debug)]
pub struct NoncovalentBondViewMut<'a> {
    noncovalent_bonds: &'a mut NoncovalentBonds,
    id: NoncovalentBondId,
}

impl<'a> NoncovalentBondViewMut<'a> {
    pub(crate) fn new(noncovalent_bonds: &'a mut NoncovalentBonds, id: NoncovalentBondId) -> Self {
        Self {
            noncovalent_bonds,
            id,
        }
    }

    #[inline]
    pub fn id(&self) -> NoncovalentBondId {
        self.id
    }

    #[inline]
    pub fn attributes(&self) -> &NoncovalentBondForm {
        self.noncovalent_bonds.attributes(self.id)
    }

    #[inline]
    pub fn attributes_mut(&mut self) -> &mut NoncovalentBondForm {
        self.noncovalent_bonds.attributes_mut(self.id)
    }

    #[inline]
    pub fn kind(&self) -> &NoncovalentBondKindForm {
        &self.attributes().kind
    }

    #[inline]
    pub fn constraints(&self) -> &NoncovalentBondConstraintsForm {
        &self.attributes().constraints
    }

    #[inline]
    pub fn atom_ids(&self) -> [AtomId; 2] {
        self.noncovalent_bonds.atoms(self.id)
    }
}

/// Mutable editor access to a noncovalent bond.
///
/// Structural mutations update incidence and preserve attributes and constraints.
/// Atom existence, distinctness, and uniqueness of endpoint pairs are checked at publication.
#[derive(Debug)]
pub struct NoncovalentBondEditorViewMut<'a> {
    noncovalent_bonds: &'a mut NoncovalentBonds,
    id: NoncovalentBondId,
}

impl<'a> NoncovalentBondEditorViewMut<'a> {
    pub(crate) fn new(noncovalent_bonds: &'a mut NoncovalentBonds, id: NoncovalentBondId) -> Self {
        Self {
            noncovalent_bonds,
            id,
        }
    }

    #[inline]
    pub fn id(&self) -> NoncovalentBondId {
        self.id
    }

    #[inline]
    pub fn attributes(&self) -> &NoncovalentBondForm {
        self.noncovalent_bonds.attributes(self.id)
    }

    #[inline]
    pub fn attributes_mut(&mut self) -> &mut NoncovalentBondForm {
        self.noncovalent_bonds.attributes_mut(self.id)
    }

    #[inline]
    pub fn kind(&self) -> &NoncovalentBondKindForm {
        &self.attributes().kind
    }

    #[inline]
    pub fn constraints(&self) -> &NoncovalentBondConstraintsForm {
        &self.attributes().constraints
    }

    #[inline]
    pub fn atom_ids(&self) -> [AtomId; 2] {
        self.noncovalent_bonds.atoms(self.id)
    }

    /// Replace both endpoints, preserving the supplied order.
    pub fn replace_atoms(&mut self, atoms: [AtomId; 2]) {
        self.noncovalent_bonds.replace_atoms(self.id, atoms);
    }

    /// Replace the endpoint at `position`.
    ///
    /// # Panics
    ///
    /// Panics when `position` is not 0 or 1.
    pub fn replace_atom(&mut self, position: AtomPosition, atom: AtomId) {
        self.noncovalent_bonds
            .replace_atom(self.id, position.index(), atom);
    }
}

/// Stored constraint container of `bond`.
pub(crate) fn noncovalent_bond_asserted_constraints(
    molecule: &Molecule,
    id: NoncovalentBondId,
) -> &NoncovalentBondConstraintsForm {
    &molecule.noncovalent_bond(id).attributes().constraints
}

/// Derived side of one noncovalent-bond constraint key: intramolecularity is
/// whether the two endpoints share a localized-bond component — derived from
/// the topology with no absence cell, so both modes agree.
pub(crate) fn noncovalent_bond_derived_constraint(
    molecule: &Molecule,
    id: NoncovalentBondId,
    key: NoncovalentBondConstraintKey,
    _complete: bool,
) -> Option<NoncovalentBondConstraintForm> {
    match key {
        NoncovalentBondConstraintKey::Intramolecular => {
            let [a, b] = molecule.noncovalent_bond(id).atom_ids();
            Some(NoncovalentBondConstraintForm::intramolecular(
                same_bond_component(molecule, a, b),
            ))
        }
    }
}

/// Whether `a` and `b` lie in one localized-bond component, by breadth-first
/// reachability over localized bonds.
fn same_bond_component(molecule: &Molecule, a: AtomId, b: AtomId) -> bool {
    if a == b {
        return true;
    }
    let mut visited = vec![false; molecule.atoms().count()];
    visited[a.index()] = true;
    let mut queue = vec![a];
    while let Some(current) = queue.pop() {
        for neighbor in molecule.neighbors(current) {
            let next = neighbor.atom_id();
            if next == b {
                return true;
            }
            if !visited[next.index()] {
                visited[next.index()] = true;
                queue.push(next);
            }
        }
    }
    false
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
    use crate::ir::constraint::NoncovalentBondConstraintForm;
    use crate::ir::dative::DativeBondForm;
    use crate::ir::entity::Entity;
    use crate::ir::id::{AtomId, AtomPosition, NoncovalentBondId};
    use crate::ir::molecule::{Molecule, MoleculeEntries, MoleculeIntegrityError};
    use crate::ir::multicenter::MulticenterBondForm;
    use crate::ir::noncovalent::{
        NoncovalentBondForm, NoncovalentBondKind, NoncovalentBondKindForm,
    };

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
    fn noncovalent_molecule() -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::default(); 4],
            noncovalent: vec![
                (
                    [AtomId(2), AtomId(0)],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                ),
                (
                    [AtomId(3), AtomId(1)],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HalogenBond),
                ),
            ],
            ..Default::default()
        })
    }

    #[fixture]
    fn noncovalent_entries() -> MoleculeEntries {
        MoleculeEntries {
            atoms: vec![AtomForm::default(); 5],
            noncovalent: vec![
                (
                    [AtomId(2), AtomId(0)],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond)
                        .with_constraint(NoncovalentBondConstraintForm::intramolecular(false)),
                ),
                (
                    [AtomId(3), AtomId(1)],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HalogenBond),
                ),
            ],
            ..Default::default()
        }
    }

    #[rstest]
    fn test_noncovalent_bond_views_count(molecule: Molecule) {
        assert_eq!(molecule.noncovalent_bonds().count(), 1);
    }

    #[rstest]
    fn test_noncovalent_bond_views_ids(molecule: Molecule) {
        assert_exact_size_by(
            Molecule::default().noncovalent_bonds().ids(),
            vec![],
            |id| id,
        );
        assert_exact_size_by(
            molecule.noncovalent_bonds().ids(),
            vec![NoncovalentBondId(0)],
            |id| id,
        );
    }

    #[rstest]
    fn test_noncovalent_bond_views_iter(molecule: Molecule) {
        assert_exact_size_by(
            Molecule::default().noncovalent_bonds().iter(),
            vec![],
            |view| (view.id(), view.atom_ids(), view.attributes().clone()),
        );
        assert_exact_size_by(
            molecule.noncovalent_bonds().iter(),
            vec![(
                NoncovalentBondId(0),
                [AtomId(0), AtomId(3)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            |view| (view.id(), view.atom_ids(), view.attributes().clone()),
        );
    }

    #[rstest]
    #[case::participant(AtomId(0), vec![NoncovalentBondId(0)])]
    #[case::uninvolved(AtomId(1), vec![])]
    fn test_noncovalent_bond_views_incident(
        molecule: Molecule,
        #[case] atom: AtomId,
        #[case] expected: Vec<NoncovalentBondId>,
    ) {
        assert_exact_size_by(
            molecule.noncovalent_bonds().incident_ids(atom),
            expected.clone(),
            |id| id,
        );
        assert_exact_size_by(
            molecule.noncovalent_bonds().incident(atom),
            expected,
            |view| view.id(),
        );
    }

    #[rstest]
    #[case::present(NoncovalentBondId(0), true)]
    #[case::absent(NoncovalentBondId(99), false)]
    fn test_noncovalent_bond_views_contains(
        molecule: Molecule,
        #[case] id: NoncovalentBondId,
        #[case] expected: bool,
    ) {
        assert_eq!(molecule.noncovalent_bonds().contains(id), expected);
    }

    #[rstest]
    fn test_noncovalent_bond_views_get(molecule: Molecule) {
        let res = molecule.noncovalent_bonds().get(NoncovalentBondId(0));
        assert!(res.is_some());
        let view = res.unwrap();
        assert_eq!(view.id(), NoncovalentBondId(0));
        assert_eq!(view.atom_ids(), [AtomId(0), AtomId(3)]);
    }

    #[rstest]
    fn test_noncovalent_bond_views_get_none(molecule: Molecule) {
        let res = molecule.noncovalent_bonds().get(NoncovalentBondId(99));
        assert!(res.is_none());
    }

    #[rstest]
    #[case(NoncovalentBondId(0))]
    #[case(NoncovalentBondId(1))]
    fn test_noncovalent_bond_view_id(
        #[from(noncovalent_molecule)] molecule: Molecule,
        #[case] id: NoncovalentBondId,
    ) {
        assert_eq!(molecule.noncovalent_bond(id).id(), id);
    }

    #[rstest]
    #[case(NoncovalentBondId(0), NoncovalentBondForm { kind: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), ..Default::default() })]
    #[case(NoncovalentBondId(1), NoncovalentBondForm { kind: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HalogenBond), ..Default::default() })]
    fn test_noncovalent_bond_view_attributes(
        #[from(noncovalent_molecule)] molecule: Molecule,
        #[case] id: NoncovalentBondId,
        #[case] expected: NoncovalentBondForm,
    ) {
        let attributes = {
            let view = molecule.noncovalent_bond(id);
            view.attributes()
        };
        assert_eq!(attributes, &expected);
    }

    #[rstest]
    #[case(NoncovalentBondId(0), [AtomId(2), AtomId(0)])]
    #[case(NoncovalentBondId(1), [AtomId(3), AtomId(1)])]
    fn test_noncovalent_bond_view_atom_ids_order(
        #[from(noncovalent_molecule)] molecule: Molecule,
        #[case] id: NoncovalentBondId,
        #[case] expected: [AtomId; 2],
    ) {
        let atom_ids = {
            let view = molecule.noncovalent_bond(id);
            view.atom_ids()
        };
        assert_eq!(atom_ids, expected);
    }

    #[rstest]
    fn test_noncovalent_bond_view_atom_ids(molecule: Molecule) {
        assert_eq!(
            molecule.noncovalent_bond(NoncovalentBondId(0)).atom_ids(),
            [AtomId(0), AtomId(3)],
        );
    }

    #[rstest]
    fn test_noncovalent_bond_view_atoms(molecule: Molecule) {
        let ids = molecule
            .noncovalent_bond(NoncovalentBondId(0))
            .atoms()
            .map(|v| v.id());
        assert_eq!(ids, [AtomId(0), AtomId(3)]);
    }

    #[rstest]
    #[case(NoncovalentBondId(0))]
    #[case(NoncovalentBondId(1))]
    fn test_noncovalent_bond_editor_view_id(
        #[from(noncovalent_molecule)] molecule: Molecule,
        #[case] id: NoncovalentBondId,
    ) {
        let editor = molecule.edit();
        assert_eq!(editor.noncovalent_bond(id).id(), id);
    }

    #[rstest]
    #[case(NoncovalentBondId(0), NoncovalentBondForm { kind: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), ..Default::default() })]
    #[case(NoncovalentBondId(1), NoncovalentBondForm { kind: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HalogenBond), ..Default::default() })]
    fn test_noncovalent_bond_editor_view_attributes(
        #[from(noncovalent_molecule)] molecule: Molecule,
        #[case] id: NoncovalentBondId,
        #[case] expected: NoncovalentBondForm,
    ) {
        let editor = molecule.edit();
        let attributes = {
            let view = editor.noncovalent_bond(id);
            view.attributes()
        };
        assert_eq!(attributes, &expected);
    }

    #[rstest]
    #[case::literal(NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond))]
    #[case::undetermined(NoncovalentBondKindForm::Undetermined)]
    fn test_noncovalent_bond_editor_view_kind(
        mut noncovalent_entries: MoleculeEntries,
        #[case] expected: NoncovalentBondKindForm,
    ) {
        noncovalent_entries.noncovalent[0].1.kind = expected.clone();
        let editor = Molecule::from_entries(noncovalent_entries).edit();
        let kind = {
            let view = editor.noncovalent_bond(NoncovalentBondId(0));
            view.kind()
        };
        assert_eq!(kind, &expected);
    }

    #[rstest]
    #[case(NoncovalentBondId(0), [AtomId(2), AtomId(0)])]
    #[case(NoncovalentBondId(1), [AtomId(3), AtomId(1)])]
    fn test_noncovalent_bond_editor_view_atom_ids_order(
        #[from(noncovalent_molecule)] molecule: Molecule,
        #[case] id: NoncovalentBondId,
        #[case] expected: [AtomId; 2],
    ) {
        let editor = molecule.edit();
        let atom_ids = {
            let view = editor.noncovalent_bond(id);
            view.atom_ids()
        };
        assert_eq!(atom_ids, expected);
    }

    #[rstest]
    #[case(NoncovalentBondId(0))]
    #[case(NoncovalentBondId(1))]
    fn test_noncovalent_bond_view_mut_attributes_mut(
        #[from(noncovalent_molecule)] mut molecule: Molecule,
        #[case] id: NoncovalentBondId,
    ) {
        let expected = NoncovalentBondForm {
            kind: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
            ..Default::default()
        };
        {
            let mut view = molecule.noncovalent_bond_mut(id);
            assert_eq!(view.id(), id);
            *view.attributes_mut() = expected.clone();
            assert_eq!(view.attributes(), &expected);
            assert_eq!(view.kind(), &expected.kind);
            view.attributes_mut().kind = NoncovalentBondKindForm::Undetermined;
            assert_eq!(view.kind(), &NoncovalentBondKindForm::Undetermined);
            view.attributes_mut().kind = expected.kind.clone();
        }
        assert_eq!(molecule.noncovalent_bond(id).attributes(), &expected);
    }

    #[rstest]
    #[case(NoncovalentBondId(0), [AtomId(2), AtomId(0)])]
    #[case(NoncovalentBondId(1), [AtomId(3), AtomId(1)])]
    fn test_noncovalent_bond_view_mut_atom_ids_order(
        #[from(noncovalent_molecule)] mut molecule: Molecule,
        #[case] id: NoncovalentBondId,
        #[case] expected: [AtomId; 2],
    ) {
        let view = molecule.noncovalent_bond_mut(id);
        assert_eq!(view.atom_ids(), expected);
    }

    #[rstest]
    #[case(NoncovalentBondId(0))]
    #[case(NoncovalentBondId(1))]
    fn test_noncovalent_bond_editor_view_mut_attributes_mut(
        #[from(noncovalent_molecule)] molecule: Molecule,
        #[case] id: NoncovalentBondId,
    ) {
        let mut editor = molecule.edit();
        let expected = NoncovalentBondForm {
            kind: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
            ..Default::default()
        };
        {
            let mut view = editor.noncovalent_bond_mut(id);
            assert_eq!(view.id(), id);
            *view.attributes_mut() = expected.clone();
            assert_eq!(view.attributes(), &expected);
            assert_eq!(view.kind(), &expected.kind);
            view.attributes_mut().kind = NoncovalentBondKindForm::Undetermined;
            assert_eq!(view.kind(), &NoncovalentBondKindForm::Undetermined);
            view.attributes_mut().kind = expected.kind.clone();
        }
        assert_eq!(editor.noncovalent_bond(id).attributes(), &expected);
    }

    #[rstest]
    #[case(NoncovalentBondId(0), [AtomId(2), AtomId(0)])]
    #[case(NoncovalentBondId(1), [AtomId(3), AtomId(1)])]
    fn test_noncovalent_bond_editor_view_mut_atom_ids_order(
        #[from(noncovalent_molecule)] molecule: Molecule,
        #[case] id: NoncovalentBondId,
        #[case] expected: [AtomId; 2],
    ) {
        let mut editor = molecule.edit();
        let view = editor.noncovalent_bond_mut(id);
        assert_eq!(view.atom_ids(), expected);
    }

    #[rstest]
    #[case::reversed([AtomId(0), AtomId(2)])]
    #[case::shared_endpoint([AtomId(3), AtomId(0)])]
    #[case::disjoint([AtomId(4), AtomId(1)])]
    fn test_noncovalent_bond_editor_view_mut_replace_atoms(
        mut noncovalent_entries: MoleculeEntries,
        #[case] atoms: [AtomId; 2],
    ) {
        let mut editor = Molecule::from_entries(noncovalent_entries.clone()).edit();
        {
            let mut view = editor.noncovalent_bond_mut(NoncovalentBondId(0));
            view.replace_atoms(atoms);
            assert_eq!(view.atom_ids(), atoms);
            assert_eq!(view.attributes(), &noncovalent_entries.noncovalent[0].1);
        }
        noncovalent_entries.noncovalent[0].0 = atoms;
        assert_eq!(
            editor.finish(),
            Ok(Molecule::from_entries(noncovalent_entries))
        );
    }

    #[rstest]
    #[case::missing_atom([AtomId(5), AtomId(0)], MoleculeIntegrityError::InvalidReference { entity: Entity::Atom(AtomId(5)) })]
    #[case::repeated_atom([AtomId(0), AtomId(0)], MoleculeIntegrityError::DuplicateAtom { entity: Entity::NoncovalentBond(NoncovalentBondId(0)), atom: AtomId(0) })]
    #[case::parallel([AtomId(1), AtomId(3)], MoleculeIntegrityError::ParallelNoncovalentBonds { atoms: [AtomId(1), AtomId(3)] })]
    fn test_noncovalent_bond_editor_view_mut_replace_atoms_publication(
        noncovalent_entries: MoleculeEntries,
        #[case] atoms: [AtomId; 2],
        #[case] expected: MoleculeIntegrityError,
    ) {
        let mut editor = Molecule::from_entries(noncovalent_entries).edit();
        editor
            .noncovalent_bond_mut(NoncovalentBondId(0))
            .replace_atoms(atoms);
        assert_eq!(editor.finish(), Err(expected));
    }

    #[rstest]
    #[case::first(AtomPosition(0), [AtomId(4), AtomId(0)])]
    #[case::second(AtomPosition(1), [AtomId(2), AtomId(4)])]
    fn test_noncovalent_bond_editor_view_mut_replace_atom(
        mut noncovalent_entries: MoleculeEntries,
        #[case] position: AtomPosition,
        #[case] expected: [AtomId; 2],
    ) {
        let mut editor = Molecule::from_entries(noncovalent_entries.clone()).edit();
        {
            let mut view = editor.noncovalent_bond_mut(NoncovalentBondId(0));
            view.replace_atom(position, AtomId(4));
            assert_eq!(view.atom_ids(), expected);
            assert_eq!(view.attributes(), &noncovalent_entries.noncovalent[0].1);
        }
        noncovalent_entries.noncovalent[0].0 = expected;
        assert_eq!(
            editor.finish(),
            Ok(Molecule::from_entries(noncovalent_entries))
        );
    }

    #[rstest]
    #[case::end(AtomPosition(2))]
    #[case::beyond_end(AtomPosition(3))]
    #[should_panic]
    fn test_noncovalent_bond_editor_view_mut_replace_atom_error(
        noncovalent_entries: MoleculeEntries,
        #[case] position: AtomPosition,
    ) {
        let mut editor = Molecule::from_entries(noncovalent_entries).edit();
        editor
            .noncovalent_bond_mut(NoncovalentBondId(0))
            .replace_atom(position, AtomId(4));
    }

    #[rstest]
    fn test_noncovalent_bond_editor_view_mut_replace_atoms_incidence(
        noncovalent_entries: MoleculeEntries,
    ) {
        let mut editor = Molecule::from_entries(noncovalent_entries).edit();
        {
            let mut view = editor.noncovalent_bond_mut(NoncovalentBondId(0));
            view.replace_atoms([AtomId(4), AtomId(1)]);
            view.replace_atom(AtomPosition(1), AtomId(3));
        }
        let molecule = editor.finish().unwrap();
        let incidence: Vec<Vec<_>> = molecule
            .atoms()
            .ids()
            .map(|atom| molecule.noncovalent_bonds().incident_ids(atom).collect())
            .collect();
        assert_eq!(
            incidence,
            vec![
                vec![],
                vec![NoncovalentBondId(1)],
                vec![],
                vec![NoncovalentBondId(0), NoncovalentBondId(1)],
                vec![NoncovalentBondId(0)],
            ]
        );
    }
}
