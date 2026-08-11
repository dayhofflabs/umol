//! Aromatic system views.

use std::collections::HashSet;

use umol_graph_core::{NodeId, RelationId, Unordered, VarRelationSet};

use super::super::aromatic::AromaticSystemForm;
use super::super::constraint::AromaticSystemConstraintsForm;
use super::super::correspondence::MoleculeCorrespondence;
use super::super::electrons::ElectronCountsForm;
use super::super::id::{AromaticSystemId, AtomId, BondId};
use super::super::molecule::Molecule;
use super::super::num::NumForm;
use super::super::spin::UnpairedElectronsForm;
use super::super::traits::Lattice;
use super::atom::AtomView;
use super::bond::BondView;

/// Namespace accessor for aromatic-system views on a `Molecule`.
#[derive(Clone, Copy)]
pub struct AromaticSystemViews<'a> {
    molecule: &'a Molecule,
    aromatic_systems: &'a VarRelationSet<NodeId, Unordered, AromaticSystemForm>,
}

impl<'a> AromaticSystemViews<'a> {
    pub(crate) fn new(
        molecule: &'a Molecule,
        aromatic_systems: &'a VarRelationSet<NodeId, Unordered, AromaticSystemForm>,
    ) -> Self {
        Self {
            molecule,
            aromatic_systems,
        }
    }

    pub fn count(&self) -> usize {
        self.aromatic_systems.count()
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = AromaticSystemId> {
        self.aromatic_systems
            .relation_ids()
            .map(AromaticSystemId::from)
    }

    pub fn iter(&self) -> impl ExactSizeIterator<Item = AromaticSystemView<'a>> {
        let molecule = self.molecule;
        let set = self.aromatic_systems;
        set.relation_ids().map(move |rid| AromaticSystemView {
            id: AromaticSystemId::from(rid),
            attributes: set.data(rid),
            atoms: set.participants(rid),
            molecule,
        })
    }

    pub fn contains(&self, id: AromaticSystemId) -> bool {
        self.aromatic_systems.contains(RelationId::from(id))
    }

    pub fn get(&self, id: AromaticSystemId) -> Option<AromaticSystemView<'a>> {
        if !self.contains(id) {
            return None;
        }
        let rid = RelationId::from(id);
        Some(AromaticSystemView {
            id,
            attributes: self.aromatic_systems.data(rid),
            atoms: self.aromatic_systems.participants(rid),
            molecule: self.molecule,
        })
    }

    /// Ids of aromatic systems incident on `atom`.
    pub fn incident_ids(
        &self,
        atom: AtomId,
    ) -> impl ExactSizeIterator<Item = AromaticSystemId> + 'a {
        self.aromatic_systems
            .incident(NodeId::from(atom))
            .iter()
            .map(|&rid| AromaticSystemId::from(rid))
    }

    /// Whether any aromatic system is incident on `atom`.
    pub fn has_incident(&self, atom: AtomId) -> bool {
        self.aromatic_systems.has_incident(NodeId::from(atom))
    }

    /// Views of aromatic systems incident on `atom`.
    pub fn incident(
        &self,
        atom: AtomId,
    ) -> impl ExactSizeIterator<Item = AromaticSystemView<'a>> + 'a {
        let molecule = self.molecule;
        let set = self.aromatic_systems;
        self.incident_ids(atom).map(move |id| {
            let rid = RelationId::from(id);
            AromaticSystemView {
                id,
                attributes: set.data(rid),
                atoms: set.participants(rid),
                molecule,
            }
        })
    }

    /// Id of the aromatic system whose atom set equals `atoms`, if any.
    pub fn of_id(&self, atoms: impl IntoIterator<Item = AtomId>) -> Option<AromaticSystemId> {
        let nodes: Vec<NodeId> = atoms.into_iter().map(NodeId::from).collect();
        self.aromatic_systems
            .find_by_participants(&nodes)
            .map(AromaticSystemId::from)
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
        self.aromatic_systems
            .relation_ids()
            .filter(|&rid| {
                self.aromatic_systems
                    .participants(rid)
                    .iter()
                    .all(|p| set.contains(p))
            })
            .map(AromaticSystemId::from)
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
    pub id: AromaticSystemId,
    atoms: &'a [NodeId],
    pub attributes: &'a AromaticSystemForm,
    molecule: &'a Molecule,
}

impl<'a> AromaticSystemView<'a> {
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

    #[inline]
    pub fn constraints(&self) -> &'a AromaticSystemConstraintsForm {
        &self.attributes.constraints
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

    pub fn bond_ids(&self) -> impl Iterator<Item = BondId> + 'a {
        self.molecule
            .raw_graph()
            .induced_edges(self.atoms)
            .map(BondId::from)
    }

    pub fn bonds(&self) -> impl Iterator<Item = BondView<'a>> + 'a {
        let molecule = self.molecule;
        self.molecule
            .raw_graph()
            .induced_edges(self.atoms)
            .map(move |edge| molecule.bond(BondId::from(edge)))
    }

    /// The molecule subgraph induced by this system's atoms, as a sub-to-host correspondence.
    pub fn induced_subgraph(&self) -> MoleculeCorrespondence {
        self.molecule
            .induced_subgraph(&self.atom_ids().collect::<Vec<_>>())
    }

    /// Sum of per-atom electron contributions on this aromatic system.
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
        self.atoms
            .iter()
            .map(|&n| AtomId::from(n))
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
            .induced_edges(self.atoms)
            .map(BondId::from)
            .filter(move |b| subset.contains(b))
            .map(move |id| molecule.bond(id))
    }

    /// Is aromatic system ground
    pub fn is_ground(&self) -> bool {
        self.attributes.is_ground()
    }

    /// Is aromatic system undetermined
    pub fn is_undetermined(&self) -> bool {
        self.attributes.is_undetermined()
    }
}

/// Mutable borrowed view of an aromatic system: its id, member atoms (owned)
/// and mutable data. Molecule-scope peer of `AromaticSystemView`.
#[derive(Debug)]
pub struct AromaticSystemViewMut<'a> {
    pub id: AromaticSystemId,
    pub atoms: Vec<AtomId>,
    pub attributes: &'a mut AromaticSystemForm,
}

// Editor-scope view bundles for aromatic systems.

pub struct AromaticSystemEditorView<'a> {
    pub id: AromaticSystemId,
    atoms: &'a [NodeId],
    pub attributes: &'a AromaticSystemForm,
}

impl<'a> AromaticSystemEditorView<'a> {
    pub(crate) fn new(
        id: AromaticSystemId,
        atoms: &'a [NodeId],
        attributes: &'a AromaticSystemForm,
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

pub struct AromaticSystemEditorViewMut<'a> {
    pub id: AromaticSystemId,
    atoms: &'a [NodeId],
    pub attributes: &'a mut AromaticSystemForm,
}

impl<'a> AromaticSystemEditorViewMut<'a> {
    pub(crate) fn new(
        id: AromaticSystemId,
        atoms: &'a [NodeId],
        attributes: &'a mut AromaticSystemForm,
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
    use super::{AromaticSystemEditorView, AromaticSystemEditorViewMut};
    use crate::ir::aromatic::AromaticSystemForm;
    use crate::ir::atom::AtomForm;
    use crate::ir::bond::BondForm;
    use crate::ir::dative::DativeBondForm;
    use crate::ir::id::{AromaticSystemId, AtomId, BondId};
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
                AtomId(0),
                AtomId(3),
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            ..Default::default()
        })
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
            |atom| atom.id,
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
            .map(|v| v.id)
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
            .map(|v| v.id)
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
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, expected);
    }

    #[rstest]
    fn test_aromatic_system_editor_view_atom_ids() {
        let atoms = [NodeId(0), NodeId(1), NodeId(2)];
        let attributes = AromaticSystemForm::default();
        let view = AromaticSystemEditorView::new(AromaticSystemId(0), &atoms, &attributes);
        assert_exact_size_by(
            view.atom_ids(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
            |id| id,
        );
    }

    #[rstest]
    fn test_aromatic_system_editor_view_mut_atom_ids() {
        let atoms = [NodeId(0), NodeId(1), NodeId(2)];
        let mut attributes = AromaticSystemForm::default();
        let view = AromaticSystemEditorViewMut::new(AromaticSystemId(0), &atoms, &mut attributes);
        assert_exact_size_by(
            view.atom_ids(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
            |id| id,
        );
    }
}
