//! Aromatic system views.

use std::collections::HashSet;

use umol_graph_core::{NodeId, RelationId, Unordered, VarRelationSet};

use super::super::aromatic::AromaticSystemForm;
use super::super::constraint::AromaticSystemConstraintsForm;
use super::super::correspondence::MoleculeCorrespondence;
use super::super::electrons::ElectronCountsForm;
use super::super::id::{AromaticSystemId, AtomId, BondId};
use super::super::molecule::MoleculeAst;
use super::super::spin::UnpairedElectronsForm;
use super::super::traits::Lattice;
use super::super::value::NumForm;
use super::atom::AtomView;
use super::bond::BondView;

/// Namespace accessor for aromatic-system views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct AromaticSystemViews<'a> {
    molecule: &'a MoleculeAst,
    aromatic_systems: &'a VarRelationSet<NodeId, Unordered, AromaticSystemForm>,
}

impl<'a> AromaticSystemViews<'a> {
    pub(crate) fn new(
        molecule: &'a MoleculeAst,
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
            ast: set.data(rid),
            atoms: set.participants(rid),
            molecule,
        })
    }

    /// Whether a system repeats a participant or two systems share an atom — the "aromatic systems are
    /// vertex-disjoint" structural conflict; the aromatic peer of [`StereoAtomViews::has_conflict`].
    pub fn has_conflict(&self) -> bool {
        let mut seen: HashSet<AtomId> = HashSet::new();
        for view in self.iter() {
            for atom in view.atom_ids() {
                if !seen.insert(atom) {
                    return true;
                }
            }
        }
        false
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
            ast: self.aromatic_systems.data(rid),
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
                ast: set.data(rid),
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
    pub ast: &'a AromaticSystemForm,
    molecule: &'a MoleculeAst,
}

impl<'a> AromaticSystemView<'a> {
    #[inline]
    pub fn electrons(&self) -> &'a ElectronCountsForm {
        &self.ast.electrons
    }

    #[inline]
    pub fn charge(&self) -> &'a NumForm {
        &self.ast.charge
    }

    #[inline]
    pub fn unpaired_electrons(&self) -> &'a UnpairedElectronsForm {
        &self.ast.unpaired_electrons
    }

    #[inline]
    pub fn constraints(&self) -> &'a AromaticSystemConstraintsForm {
        &self.ast.constraints
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
        match &self.ast.electrons {
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
        self.ast.is_ground()
    }

    /// Is aromatic system undetermined
    pub fn is_undetermined(&self) -> bool {
        self.ast.is_undetermined()
    }
}

/// Mutable borrowed view of an aromatic system: its id, member atoms (owned)
/// and mutable data. Molecule-scope peer of `AromaticSystemView`.
#[derive(Debug)]
pub struct AromaticSystemViewMut<'a> {
    pub id: AromaticSystemId,
    pub atoms: Vec<AtomId>,
    pub ast: &'a mut AromaticSystemForm,
}

// Editor-scope view bundles for aromatic systems.

pub struct AromaticSystemEditorView<'a> {
    pub id: AromaticSystemId,
    atoms: &'a [NodeId],
    pub ast: &'a AromaticSystemForm,
}

impl<'a> AromaticSystemEditorView<'a> {
    pub(crate) fn new(
        id: AromaticSystemId,
        atoms: &'a [NodeId],
        ast: &'a AromaticSystemForm,
    ) -> Self {
        Self { id, atoms, ast }
    }

    pub fn atom_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + 'a {
        self.atoms.iter().map(|&n| AtomId::from(n))
    }
}

pub struct AromaticSystemEditorViewMut<'a> {
    pub id: AromaticSystemId,
    atoms: &'a [NodeId],
    pub ast: &'a mut AromaticSystemForm,
}

impl<'a> AromaticSystemEditorViewMut<'a> {
    pub(crate) fn new(
        id: AromaticSystemId,
        atoms: &'a [NodeId],
        ast: &'a mut AromaticSystemForm,
    ) -> Self {
        Self { id, atoms, ast }
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
    use crate::ir::molecule::{MoleculeAst, MoleculeEntries};
    use crate::ir::multicenter::MulticenterBondForm;
    use crate::ir::noncovalent::{NoncovalentBondForm, NoncovalentBondKind};
    use crate::ir::value::NumForm;

    #[fixture]
    fn molecule() -> MoleculeAst {
        MoleculeAst::from_entries(MoleculeEntries {
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
    fn test_aromatic_system_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.aromatic_systems().count(), 1);
    }

    #[rstest]
    fn test_aromatic_system_views_ids(molecule: MoleculeAst) {
        assert_exact_size_by(
            MoleculeAst::default().aromatic_systems().ids(),
            vec![],
            |id| id,
        );
        assert_exact_size_by(
            molecule.aromatic_systems().ids(),
            vec![AromaticSystemId(0)],
            |id| id,
        );
    }

    #[rstest]
    fn test_aromatic_system_views_iter(molecule: MoleculeAst) {
        assert_exact_size_by(
            MoleculeAst::default().aromatic_systems().iter(),
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
        molecule: MoleculeAst,
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
        molecule: MoleculeAst,
        #[case] id: AromaticSystemId,
        #[case] expected: bool,
    ) {
        assert_eq!(molecule.aromatic_systems().contains(id), expected);
    }

    #[rstest]
    fn test_aromatic_system_views_get(molecule: MoleculeAst) {
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
    fn test_aromatic_system_views_get_none(molecule: MoleculeAst) {
        let res = molecule.aromatic_systems().get(AromaticSystemId(99));
        assert!(res.is_none());
    }

    #[rstest]
    fn test_aromatic_system_view_atom_ids(molecule: MoleculeAst) {
        assert_exact_size_by(
            molecule.aromatic_system(AromaticSystemId(0)).atom_ids(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
            |id| id,
        );
    }

    #[rstest]
    fn test_aromatic_system_view_atoms(molecule: MoleculeAst) {
        assert_exact_size_by(
            molecule.aromatic_system(AromaticSystemId(0)).atoms(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
            |atom| atom.id,
        );
    }

    #[rstest]
    fn test_aromatic_system_view_bond_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .aromatic_system(AromaticSystemId(0))
                .bond_ids()
                .collect::<Vec<_>>(),
            vec![BondId(0), BondId(1)],
        );
    }

    #[rstest]
    fn test_aromatic_system_view_bonds(molecule: MoleculeAst) {
        let ids: Vec<BondId> = molecule
            .aromatic_system(AromaticSystemId(0))
            .bonds()
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, vec![BondId(0), BondId(1)]);
    }

    #[rstest]
    fn test_aromatic_system_view_induced_subgraph(molecule: MoleculeAst) {
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
    fn test_aromatic_system_view_electron_count(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .aromatic_system(AromaticSystemId(0))
                .electron_count(),
            NumForm::Undetermined,
        );
    }

    #[rstest]
    fn test_aromatic_system_view_atom_count(molecule: MoleculeAst) {
        assert_eq!(
            molecule.aromatic_system(AromaticSystemId(0)).atom_count(),
            3
        );
    }

    #[rstest]
    fn test_aromatic_system_view_bond_count(molecule: MoleculeAst) {
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
        molecule: MoleculeAst,
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
        molecule: MoleculeAst,
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
        let ast = AromaticSystemForm::default();
        let view = AromaticSystemEditorView::new(AromaticSystemId(0), &atoms, &ast);
        assert_exact_size_by(
            view.atom_ids(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
            |id| id,
        );
    }

    #[rstest]
    fn test_aromatic_system_editor_view_mut_atom_ids() {
        let atoms = [NodeId(0), NodeId(1), NodeId(2)];
        let mut ast = AromaticSystemForm::default();
        let view = AromaticSystemEditorViewMut::new(AromaticSystemId(0), &atoms, &mut ast);
        assert_exact_size_by(
            view.atom_ids(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
            |id| id,
        );
    }
}
