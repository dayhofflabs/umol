//! Multicenter bond views.

use std::collections::{BTreeSet, HashSet};

use umol_graph_core::{NodeId, RelationId, Unordered, VarRelationSet};

use super::super::constraint::MulticenterBondConstraintsForm;
use super::super::electrons::ElectronCountsForm;
use super::super::id::{AtomId, MulticenterBondId};
use super::super::molecule::MoleculeAst;
use super::super::multicenter::MulticenterBondForm;
use super::super::spin::UnpairedElectronsForm;
use super::super::traits::Lattice;
use super::super::value::NumForm;
use super::atom::AtomView;

/// Namespace accessor for multicenter-bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct MulticenterBondViews<'a> {
    molecule: &'a MoleculeAst,
    multicenter_bonds: &'a VarRelationSet<NodeId, Unordered, MulticenterBondForm>,
}

impl<'a> MulticenterBondViews<'a> {
    pub(crate) fn new(
        molecule: &'a MoleculeAst,
        multicenter_bonds: &'a VarRelationSet<NodeId, Unordered, MulticenterBondForm>,
    ) -> Self {
        Self {
            molecule,
            multicenter_bonds,
        }
    }

    pub fn count(&self) -> usize {
        self.multicenter_bonds.count()
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = MulticenterBondId> {
        self.multicenter_bonds
            .relation_ids()
            .map(MulticenterBondId::from)
    }

    pub fn iter(&self) -> impl ExactSizeIterator<Item = MulticenterBondView<'a>> {
        let molecule = self.molecule;
        let set = self.multicenter_bonds;
        set.relation_ids().map(move |rid| MulticenterBondView {
            id: MulticenterBondId::from(rid),
            ast: set.data(rid),
            atoms: set.participants(rid),
            molecule,
        })
    }

    /// Whether a bond repeats a participant, or two bonds have an identical participant set (partial
    /// overlap between distinct sets is allowed). Emit-compliance peer of
    /// [`StereoAtomViews::has_conflict`].
    pub fn has_conflict(&self) -> bool {
        let mut seen_sets: HashSet<BTreeSet<AtomId>> = HashSet::new();
        for view in self.iter() {
            let mut set: BTreeSet<AtomId> = BTreeSet::new();
            for atom in view.atom_ids() {
                if !set.insert(atom) {
                    return true;
                }
            }
            if !seen_sets.insert(set) {
                return true;
            }
        }
        false
    }

    pub fn contains(&self, id: MulticenterBondId) -> bool {
        self.multicenter_bonds.contains(RelationId::from(id))
    }

    pub fn get(&self, id: MulticenterBondId) -> Option<MulticenterBondView<'a>> {
        if !self.contains(id) {
            return None;
        }
        let rid = RelationId::from(id);
        Some(MulticenterBondView {
            id,
            ast: self.multicenter_bonds.data(rid),
            atoms: self.multicenter_bonds.participants(rid),
            molecule: self.molecule,
        })
    }

    /// Ids of multicenter bonds incident on `atom`.
    pub fn incident_ids(
        &self,
        atom: AtomId,
    ) -> impl ExactSizeIterator<Item = MulticenterBondId> + 'a {
        self.multicenter_bonds
            .incident(NodeId::from(atom))
            .iter()
            .map(|&rid| MulticenterBondId::from(rid))
    }

    /// Whether any multicenter bond is incident on `atom`.
    pub fn has_incident(&self, atom: AtomId) -> bool {
        self.multicenter_bonds.has_incident(NodeId::from(atom))
    }

    /// Views of multicenter bonds incident on `atom`.
    pub fn incident(
        &self,
        atom: AtomId,
    ) -> impl ExactSizeIterator<Item = MulticenterBondView<'a>> + 'a {
        let molecule = self.molecule;
        let set = self.multicenter_bonds;
        self.incident_ids(atom).map(move |id| {
            let rid = RelationId::from(id);
            MulticenterBondView {
                id,
                ast: set.data(rid),
                atoms: set.participants(rid),
                molecule,
            }
        })
    }

    /// ID of the multicenter bond whose participant set equals `atoms`, if any.
    pub fn of_id(&self, atoms: impl IntoIterator<Item = AtomId>) -> Option<MulticenterBondId> {
        let nodes: Vec<NodeId> = atoms.into_iter().map(NodeId::from).collect();
        self.multicenter_bonds
            .find_by_participants(&nodes)
            .map(MulticenterBondId::from)
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
            .relation_ids()
            .filter(|&rid| {
                self.multicenter_bonds
                    .participants(rid)
                    .iter()
                    .all(|p| set.contains(p))
            })
            .map(MulticenterBondId::from)
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
    pub ast: &'a MulticenterBondForm,
    molecule: &'a MoleculeAst,
}

impl<'a> MulticenterBondView<'a> {
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
    pub fn constraints(&self) -> &'a MulticenterBondConstraintsForm {
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

    /// Sum of per-atom electron contributions on this multicenter bond.
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
        self.ast.is_ground()
    }

    /// Is multicenter bond undetermined
    pub fn is_undetermined(&self) -> bool {
        self.ast.is_undetermined()
    }
}

/// Mutable borrowed view of a multicenter bond: its id, member atoms (owned)
/// and mutable data. Molecule-scope peer of `MulticenterBondView`.
#[derive(Debug)]
pub struct MulticenterBondViewMut<'a> {
    pub id: MulticenterBondId,
    pub atoms: Vec<AtomId>,
    pub ast: &'a mut MulticenterBondForm,
}

// Builder-scope view bundles for multicenter bonds.

pub struct MulticenterBondEditorView<'a> {
    pub id: MulticenterBondId,
    atoms: &'a [NodeId],
    pub ast: &'a MulticenterBondForm,
}

impl<'a> MulticenterBondEditorView<'a> {
    pub(crate) fn new(
        id: MulticenterBondId,
        atoms: &'a [NodeId],
        ast: &'a MulticenterBondForm,
    ) -> Self {
        Self { id, atoms, ast }
    }

    pub fn atom_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + 'a {
        self.atoms.iter().map(|&n| AtomId::from(n))
    }
}

pub struct MulticenterBondEditorViewMut<'a> {
    pub id: MulticenterBondId,
    atoms: &'a [NodeId],
    pub ast: &'a mut MulticenterBondForm,
}

impl<'a> MulticenterBondEditorViewMut<'a> {
    pub(crate) fn new(
        id: MulticenterBondId,
        atoms: &'a [NodeId],
        ast: &'a mut MulticenterBondForm,
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
    use super::{MulticenterBondEditorView, MulticenterBondEditorViewMut};
    use crate::ir::aromatic::AromaticSystemForm;
    use crate::ir::atom::AtomForm;
    use crate::ir::bond::BondForm;
    use crate::ir::dative::DativeBondForm;
    use crate::ir::id::{AtomId, MulticenterBondId};
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
    fn test_multicenter_bond_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.multicenter_bonds().count(), 1);
    }

    #[rstest]
    fn test_multicenter_bond_views_ids(molecule: MoleculeAst) {
        assert_exact_size_by(
            MoleculeAst::default().multicenter_bonds().ids(),
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
    fn test_multicenter_bond_views_iter(molecule: MoleculeAst) {
        assert_exact_size_by(
            MoleculeAst::default().multicenter_bonds().iter(),
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
        molecule: MoleculeAst,
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
        molecule: MoleculeAst,
        #[case] id: MulticenterBondId,
        #[case] expected: bool,
    ) {
        assert_eq!(molecule.multicenter_bonds().contains(id), expected);
    }

    #[rstest]
    fn test_multicenter_bond_views_get(molecule: MoleculeAst) {
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
    fn test_multicenter_bond_views_get_none(molecule: MoleculeAst) {
        let res = molecule.multicenter_bonds().get(MulticenterBondId(99));
        assert!(res.is_none());
    }

    #[rstest]
    fn test_multicenter_bond_view_atom_ids(molecule: MoleculeAst) {
        assert_exact_size_by(
            molecule.multicenter_bond(MulticenterBondId(0)).atom_ids(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
            |id| id,
        );
    }

    #[rstest]
    fn test_multicenter_bond_view_atoms(molecule: MoleculeAst) {
        assert_exact_size_by(
            molecule.multicenter_bond(MulticenterBondId(0)).atoms(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
            |atom| atom.id,
        );
    }

    #[rstest]
    fn test_multicenter_bond_view_electron_count(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .multicenter_bond(MulticenterBondId(0))
                .electron_count(),
            NumForm::Undetermined,
        );
    }

    #[rstest]
    fn test_multicenter_bond_view_atom_count(molecule: MoleculeAst) {
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
        molecule: MoleculeAst,
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
        let ast = MulticenterBondForm::default();
        let view = MulticenterBondEditorView::new(MulticenterBondId(0), &atoms, &ast);
        assert_exact_size_by(
            view.atom_ids(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
            |id| id,
        );
    }

    #[rstest]
    fn test_multicenter_bond_editor_view_mut_atom_ids() {
        let atoms = [NodeId(0), NodeId(1), NodeId(2)];
        let mut ast = MulticenterBondForm::default();
        let view = MulticenterBondEditorViewMut::new(MulticenterBondId(0), &atoms, &mut ast);
        assert_exact_size_by(
            view.atom_ids(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
            |id| id,
        );
    }
}
