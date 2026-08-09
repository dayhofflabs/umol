//! Noncovalent bond views.

use std::collections::HashSet;

use umol_graph_core::{FixedRelationSet, NodeId, RelationId, Unordered};

use super::super::constraint::NoncovalentBondConstraintsForm;
use super::super::id::{AtomId, NoncovalentBondId};
use super::super::molecule::Molecule;
use super::super::noncovalent::{NoncovalentBondForm, NoncovalentBondKindForm};
use super::super::traits::Lattice;
use super::atom::AtomView;

/// Namespace accessor for noncovalent-bond views on a `Molecule`.
#[derive(Clone, Copy)]
pub struct NoncovalentBondViews<'a> {
    molecule: &'a Molecule,
    noncovalent_bonds: &'a FixedRelationSet<NodeId, Unordered, NoncovalentBondForm, 2>,
}

impl<'a> NoncovalentBondViews<'a> {
    pub(crate) fn new(
        molecule: &'a Molecule,
        noncovalent_bonds: &'a FixedRelationSet<NodeId, Unordered, NoncovalentBondForm, 2>,
    ) -> Self {
        Self {
            molecule,
            noncovalent_bonds,
        }
    }

    pub fn count(&self) -> usize {
        self.noncovalent_bonds.count()
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = NoncovalentBondId> {
        self.noncovalent_bonds
            .relation_ids()
            .map(NoncovalentBondId::from)
    }

    pub fn iter(&self) -> impl ExactSizeIterator<Item = NoncovalentBondView<'a>> {
        let molecule = self.molecule;
        let set = self.noncovalent_bonds;
        set.relation_ids().map(move |rid| NoncovalentBondView {
            id: NoncovalentBondId::from(rid),
            ast: set.data(rid),
            atoms: {
                let parts = set.participants(rid);
                [parts[0], parts[1]]
            },
            molecule,
        })
    }

    /// Whether a noncovalent bond is a self-loop or two share an unordered atom pair (at most one
    /// interaction per pair). Emit-compliance peer of [`StereoAtomViews::has_conflict`].
    pub fn has_conflict(&self) -> bool {
        let mut seen: HashSet<[AtomId; 2]> = HashSet::new();
        self.iter().any(|view| {
            let [a, b] = view.atom_ids();
            a == b || !seen.insert(if a <= b { [a, b] } else { [b, a] })
        })
    }

    pub fn contains(&self, id: NoncovalentBondId) -> bool {
        self.noncovalent_bonds.contains(RelationId::from(id))
    }

    pub fn get(&self, id: NoncovalentBondId) -> Option<NoncovalentBondView<'a>> {
        if !self.contains(id) {
            return None;
        }
        let rid = RelationId::from(id);
        let parts = self.noncovalent_bonds.participants(rid);
        Some(NoncovalentBondView {
            id,
            ast: self.noncovalent_bonds.data(rid),
            atoms: [parts[0], parts[1]],
            molecule: self.molecule,
        })
    }

    /// Ids of noncovalent bonds incident on `atom`.
    pub fn incident_ids(
        &self,
        atom: AtomId,
    ) -> impl ExactSizeIterator<Item = NoncovalentBondId> + 'a {
        self.noncovalent_bonds
            .incident(NodeId::from(atom))
            .iter()
            .map(|&rid| NoncovalentBondId::from(rid))
    }

    /// Whether any noncovalent bond is incident on `atom`.
    pub fn has_incident(&self, atom: AtomId) -> bool {
        self.noncovalent_bonds.has_incident(NodeId::from(atom))
    }

    /// Views of noncovalent bonds incident on `atom`.
    pub fn incident(
        &self,
        atom: AtomId,
    ) -> impl ExactSizeIterator<Item = NoncovalentBondView<'a>> + 'a {
        let molecule = self.molecule;
        let set = self.noncovalent_bonds;
        self.incident_ids(atom).map(move |id| {
            let rid = RelationId::from(id);
            let parts = set.participants(rid);
            NoncovalentBondView {
                id,
                ast: set.data(rid),
                atoms: [parts[0], parts[1]],
                molecule,
            }
        })
    }

    /// Id of the noncovalent bond between `a` and `b`, if any.
    pub fn of_id(&self, first: AtomId, second: AtomId) -> Option<NoncovalentBondId> {
        self.noncovalent_bonds
            .find_by_participants(&[NodeId::from(first), NodeId::from(second)])
            .map(NoncovalentBondId::from)
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
        let set: HashSet<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        self.noncovalent_bonds
            .relation_ids()
            .filter(|&rid| {
                self.noncovalent_bonds
                    .participants(rid)
                    .iter()
                    .all(|p| set.contains(p))
            })
            .map(NoncovalentBondId::from)
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
    pub id: NoncovalentBondId,
    atoms: [NodeId; 2],
    pub ast: &'a NoncovalentBondForm,
    molecule: &'a Molecule,
}

impl<'a> NoncovalentBondView<'a> {
    #[inline]
    pub fn kind(&self) -> &'a NoncovalentBondKindForm {
        &self.ast.kind
    }

    #[inline]
    pub fn constraints(&self) -> &'a NoncovalentBondConstraintsForm {
        &self.ast.constraints
    }

    /// The two atom ids in this noncovalent interaction.
    pub fn atom_ids(&self) -> [AtomId; 2] {
        self.atoms.map(AtomId::from)
    }

    /// Views of the two atoms in this noncovalent interaction.
    pub fn atoms(&self) -> [AtomView<'a>; 2] {
        let [a, b] = self.atom_ids();
        [self.molecule.atom(a), self.molecule.atom(b)]
    }

    /// Is noncovalent bond ground
    pub fn is_ground(&self) -> bool {
        self.ast.is_ground()
    }

    /// Is noncovalent bond undetermined
    pub fn is_undetermined(&self) -> bool {
        self.ast.is_undetermined()
    }
}

/// Mutable borrowed view of a noncovalent bond: its id, the two incident atoms
/// (owned) and mutable data. Molecule-scope peer of `NoncovalentBondView`.
#[derive(Debug)]
pub struct NoncovalentBondViewMut<'a> {
    pub id: NoncovalentBondId,
    pub atoms: [AtomId; 2],
    pub ast: &'a mut NoncovalentBondForm,
}

// Builder-scope view bundles for noncovalent bonds.

pub struct NoncovalentBondEditorView<'a> {
    pub id: NoncovalentBondId,
    pub atoms: [AtomId; 2],
    pub ast: &'a NoncovalentBondForm,
}

pub struct NoncovalentBondEditorViewMut<'a> {
    pub id: NoncovalentBondId,
    pub atoms: [AtomId; 2],
    pub ast: &'a mut NoncovalentBondForm,
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
    use crate::ir::id::{AtomId, NoncovalentBondId};
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
            |view| (view.id, view.atom_ids(), view.ast.clone()),
        );
        assert_exact_size_by(
            molecule.noncovalent_bonds().iter(),
            vec![(
                NoncovalentBondId(0),
                [AtomId(0), AtomId(3)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            |view| (view.id, view.atom_ids(), view.ast.clone()),
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
            |view| view.id,
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
        assert_eq!(view.id, NoncovalentBondId(0));
        assert_eq!(view.atom_ids(), [AtomId(0), AtomId(3)]);
    }

    #[rstest]
    fn test_noncovalent_bond_views_get_none(molecule: Molecule) {
        let res = molecule.noncovalent_bonds().get(NoncovalentBondId(99));
        assert!(res.is_none());
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
            .map(|v| v.id);
        assert_eq!(ids, [AtomId(0), AtomId(3)]);
    }
}
