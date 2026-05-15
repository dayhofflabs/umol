//! Dative bond views.

use std::collections::HashSet;
use std::ops::Index;

use umol_graph_core::{NodeId, RelationId, VarRelationSet};

use super::super::constraint::DativeBondConstraints;
use super::super::dative::DativeBondAst;
use super::super::idx::{AtomId, DativeBondId};
use super::super::molecule::MoleculeAst;
use super::super::value::ValueAst;
use super::atom::AtomView;

/// Namespace accessor for dative-bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct DativeBondViews<'a> {
    molecule: &'a MoleculeAst,
    set: &'a VarRelationSet<DativeBondAst>,
}

impl<'a> DativeBondViews<'a> {
    pub(crate) fn new(molecule: &'a MoleculeAst, set: &'a VarRelationSet<DativeBondAst>) -> Self {
        Self { molecule, set }
    }

    pub fn count(&self) -> usize {
        self.set.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = DativeBondId> {
        self.set.relation_ids().map(DativeBondId::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = DativeBondView<'a>> {
        let molecule = self.molecule;
        let set = self.set;
        set.relation_ids().map(move |rid| {
            let atoms = set.participants(rid);
            let ast = set.data(rid);
            let acceptor_id = AtomId::from(atoms[ast.acceptor_slot as usize]);
            DativeBondView {
                id: DativeBondId::from(rid),
                ast,
                acceptor_id,
                atoms,
                molecule,
            }
        })
    }

    pub fn get(&self, id: DativeBondId) -> DativeBondView<'a> {
        let rid = RelationId::from(id);
        let atoms = self.set.participants(rid);
        let ast = self.set.data(rid);
        let acceptor_id = AtomId::from(atoms[ast.acceptor_slot as usize]);
        DativeBondView {
            id,
            ast,
            acceptor_id,
            atoms,
            molecule: self.molecule,
        }
    }

    /// IDs of dative bonds incident on `atom`.
    pub fn incident_ids(&self, atom: AtomId) -> impl Iterator<Item = DativeBondId> + 'a {
        self.set
            .incident(NodeId::from(atom))
            .iter()
            .map(|&rid| DativeBondId::from(rid))
    }

    /// Views of dative bonds incident on `atom`.
    pub fn incident(&self, atom: AtomId) -> impl Iterator<Item = DativeBondView<'a>> + 'a {
        let molecule = self.molecule;
        let set = self.set;
        self.incident_ids(atom).map(move |id| {
            let rid = RelationId::from(id);
            let atoms = set.participants(rid);
            let ast = set.data(rid);
            let acceptor_id = AtomId::from(atoms[ast.acceptor_slot as usize]);
            DativeBondView {
                id,
                ast,
                acceptor_id,
                atoms,
                molecule,
            }
        })
    }

    /// ID of the dative bond whose participant set equals `atoms`, if any.
    pub fn connecting_id(&self, atoms: impl IntoIterator<Item = AtomId>) -> Option<DativeBondId> {
        let target: HashSet<AtomId> = atoms.into_iter().collect();
        let &first = target.iter().next()?;
        self.incident_ids(first).find(|&id| {
            let parts: HashSet<AtomId> = self
                .set
                .participants(RelationId::from(id))
                .iter()
                .map(|&n| AtomId::from(n))
                .collect();
            parts == target
        })
    }

    /// View of the dative bond whose participant set equals `atoms`, if any.
    pub fn connecting(
        &self,
        atoms: impl IntoIterator<Item = AtomId>,
    ) -> Option<DativeBondView<'a>> {
        self.connecting_id(atoms).map(|id| self.get(id))
    }

    /// IDs of dative bonds whose participants all lie in `atoms`.
    pub fn induced_ids(&self, atoms: &[AtomId]) -> Vec<DativeBondId> {
        let set: HashSet<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        self.set
            .relation_ids()
            .filter(|&rid| self.set.participants(rid).iter().all(|p| set.contains(p)))
            .map(DativeBondId::from)
            .collect()
    }

    /// Views of dative bonds whose participants all lie in `atoms`.
    pub fn induced(&self, atoms: &[AtomId]) -> Vec<DativeBondView<'a>> {
        self.induced_ids(atoms)
            .into_iter()
            .map(|id| self.get(id))
            .collect()
    }
}

impl<'a> Index<DativeBondId> for DativeBondViews<'a> {
    type Output = DativeBondAst;
    fn index(&self, id: DativeBondId) -> &DativeBondAst {
        self.set.data(RelationId::from(id))
    }
}

/// Borrowed view of a dative bond: index, the designated acceptor atom,
/// and underlying `DativeBondAst`. Donor atoms and the full participant
/// set are reachable through `donors()` and `atoms()`.
#[derive(Clone, Copy, Debug)]
pub struct DativeBondView<'a> {
    pub id: DativeBondId,
    pub acceptor_id: AtomId,
    atoms: &'a [NodeId],
    pub ast: &'a DativeBondAst,
    molecule: &'a MoleculeAst,
}

impl<'a> DativeBondView<'a> {
    #[inline]
    pub fn acceptor_slot(&self) -> u8 {
        self.ast.acceptor_slot
    }

    #[inline]
    pub fn order(&self) -> &'a ValueAst {
        &self.ast.order
    }

    #[inline]
    pub fn constraints(&self) -> &'a DativeBondConstraints {
        &self.ast.constraints
    }

    /// All atoms in this dative bond (donors + acceptor), sorted by `AtomId`.
    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.atoms.iter().map(|&n| AtomId::from(n))
    }

    /// Views of all atoms in this dative bond (donors + acceptor).
    pub fn atoms(&self) -> impl Iterator<Item = AtomView<'a>> + 'a {
        let molecule = self.molecule;
        self.atoms
            .iter()
            .map(move |&n| molecule.atom(AtomId::from(n)))
    }

    /// Donor atom ids (participants minus the acceptor slot).
    pub fn donor_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        let acceptor_slot = self.ast.acceptor_slot as usize;
        self.atoms
            .iter()
            .enumerate()
            .filter(move |(i, _)| *i != acceptor_slot)
            .map(|(_, &n)| AtomId::from(n))
    }

    /// Donor atom views (participants minus the acceptor slot).
    pub fn donors(&self) -> impl Iterator<Item = AtomView<'a>> + 'a {
        let molecule = self.molecule;
        self.donor_ids().map(move |id| molecule.atom(id))
    }

    /// View of the acceptor atom.
    pub fn acceptor(&self) -> AtomView<'a> {
        self.molecule.atom(self.acceptor_id)
    }

    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }
}

// Builder-scope view bundles for dative bonds.

pub struct DativeBondBuilderView<'a> {
    pub id: DativeBondId,
    pub ast: &'a DativeBondAst,
    pub(crate) atoms: &'a [NodeId],
    pub acceptor_id: AtomId,
}

impl<'a> DativeBondBuilderView<'a> {
    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.atoms.iter().map(|&n| AtomId::from(n))
    }
}

pub struct DativeBondBuilderViewMut<'a> {
    pub id: DativeBondId,
    pub ast: &'a mut DativeBondAst,
    pub(crate) atoms: &'a [NodeId],
    pub acceptor_id: AtomId,
}

impl<'a> DativeBondBuilderViewMut<'a> {
    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + '_ {
        self.atoms.iter().map(|&n| AtomId::from(n))
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_shared::element::Element;

    use crate::ast::aromatic::AromaticSystemAst;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::constraint::Constraints;
    use crate::ast::dative::DativeBondAst;
    use crate::ast::idx::{AtomId, DativeBondId};
    use crate::ast::molecule::MoleculeAst;
    use crate::ast::multicenter::MulticenterBondAst;
    use crate::ast::noncovalent::{NoncovalentBondAst, NoncovalentBondKind};
    use crate::ast::value::ValueAst;

    #[fixture]
    fn molecule() -> MoleculeAst {
        MoleculeAst::from_parts(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::O),
            ],
            vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(2)),
                (AtomId(2), AtomId(3), BondAst::from_order(1)),
            ],
            vec![(vec![AtomId(2)], AtomId(3), DativeBondAst::from_order(1))],
            vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemAst::default(),
            )],
            vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                MulticenterBondAst::default(),
            )],
            vec![(
                AtomId(0),
                AtomId(3),
                NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            Constraints::default(),
        )
    }

    #[rstest]
    fn test_dative_bond_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.dative_bonds().count(), 1);
    }

    #[rstest]
    fn test_dative_bond_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.dative_bonds().ids().collect::<Vec<_>>(),
            vec![DativeBondId(0)],
        );
    }

    #[rstest]
    fn test_dative_bond_views_iter(molecule: MoleculeAst) {
        let collected: Vec<(DativeBondId, AtomId, DativeBondAst)> = molecule
            .dative_bonds()
            .iter()
            .map(|v| (v.id, v.acceptor_id, v.ast.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![(
                DativeBondId(0),
                AtomId(3),
                DativeBondAst::from_order(1).with_acceptor_slot(1),
            )],
        );
    }

    #[rstest]
    fn test_dative_bond_views_get(molecule: MoleculeAst) {
        let view = molecule.dative_bonds().get(DativeBondId(0));
        assert_eq!(view.id, DativeBondId(0));
        assert_eq!(view.acceptor_id, AtomId(3));
    }

    #[rstest]
    fn test_dative_bond_views_index(molecule: MoleculeAst) {
        let dative: &DativeBondAst = &molecule.dative_bonds()[DativeBondId(0)];
        assert_eq!(dative.order, ValueAst::Lit(1));
    }

    #[rstest]
    fn test_dative_bond_view_atom_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .dative_bond(DativeBondId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(2), AtomId(3)],
        );
    }

    #[rstest]
    fn test_dative_bond_view_donor_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .dative_bond(DativeBondId(0))
                .donor_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(2)],
        );
    }

    #[rstest]
    fn test_dative_bond_view_acceptor_id(molecule: MoleculeAst) {
        assert_eq!(molecule.dative_bond(DativeBondId(0)).acceptor_id, AtomId(3));
    }

    #[rstest]
    fn test_dative_bond_view_atoms(molecule: MoleculeAst) {
        let ids: Vec<AtomId> = molecule
            .dative_bond(DativeBondId(0))
            .atoms()
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, vec![AtomId(2), AtomId(3)]);
    }

    #[rstest]
    fn test_dative_bond_view_donors(molecule: MoleculeAst) {
        let ids: Vec<AtomId> = molecule
            .dative_bond(DativeBondId(0))
            .donors()
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, vec![AtomId(2)]);
    }

    #[rstest]
    fn test_dative_bond_view_acceptor(molecule: MoleculeAst) {
        assert_eq!(
            molecule.dative_bond(DativeBondId(0)).acceptor().id,
            AtomId(3),
        );
    }

    #[rstest]
    fn test_dative_bond_view_atom_count(molecule: MoleculeAst) {
        assert_eq!(molecule.dative_bond(DativeBondId(0)).atom_count(), 2);
    }
}
