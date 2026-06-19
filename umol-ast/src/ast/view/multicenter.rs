//! Multicenter bond views.

use std::collections::HashSet;
use std::ops::Index;

use umol_graph_core::{NodeId, RelationId, Unordered, VarRelationSet};

use super::super::constraint::MulticenterBondConstraints;
use super::super::electrons::ElectronCountsAst;
use super::super::id::{AtomId, MulticenterBondId};
use super::super::molecule::MoleculeAst;
use super::super::multicenter::MulticenterBondAst;
use super::super::spin::SpinStateAst;
use super::super::traits::Lattice;
use super::super::value::ValueAst;
use super::atom::AtomView;

/// Namespace accessor for multicenter-bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct MulticenterBondViews<'a> {
    molecule: &'a MoleculeAst,
    multicenter_bonds: &'a VarRelationSet<NodeId, Unordered, MulticenterBondAst>,
}

impl<'a> MulticenterBondViews<'a> {
    pub(crate) fn new(
        molecule: &'a MoleculeAst,
        multicenter_bonds: &'a VarRelationSet<NodeId, Unordered, MulticenterBondAst>,
    ) -> Self {
        Self {
            molecule,
            multicenter_bonds,
        }
    }

    pub fn count(&self) -> usize {
        self.multicenter_bonds.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = MulticenterBondId> {
        self.multicenter_bonds
            .relation_ids()
            .map(MulticenterBondId::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = MulticenterBondView<'a>> {
        let molecule = self.molecule;
        let set = self.multicenter_bonds;
        set.relation_ids().map(move |rid| MulticenterBondView {
            id: MulticenterBondId::from(rid),
            ast: set.data(rid),
            atoms: set.participants(rid),
            molecule,
        })
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
    pub fn incident_ids(&self, atom: AtomId) -> impl Iterator<Item = MulticenterBondId> + 'a {
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
    pub fn incident(&self, atom: AtomId) -> impl Iterator<Item = MulticenterBondView<'a>> + 'a {
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
    pub fn connecting_id(
        &self,
        atoms: impl IntoIterator<Item = AtomId>,
    ) -> Option<MulticenterBondId> {
        let target: HashSet<AtomId> = atoms.into_iter().collect();
        let &first = target.iter().next()?;
        self.incident_ids(first).find(|&id| {
            let parts: HashSet<AtomId> = self
                .multicenter_bonds
                .participants(RelationId::from(id))
                .iter()
                .map(|&n| AtomId::from(n))
                .collect();
            parts == target
        })
    }

    /// View of the multicenter bond whose participant set equals `atoms`, if any.
    pub fn connecting(
        &self,
        atoms: impl IntoIterator<Item = AtomId>,
    ) -> Option<MulticenterBondView<'a>> {
        self.connecting_id(atoms).map(|id| {
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

impl<'a> Index<MulticenterBondId> for MulticenterBondViews<'a> {
    type Output = MulticenterBondAst;
    fn index(&self, id: MulticenterBondId) -> &MulticenterBondAst {
        self.multicenter_bonds.data(RelationId::from(id))
    }
}

/// Borrowed view of a multicenter bond: its index, member atoms via
/// `atoms()`, and underlying `MulticenterBondAst`.
#[derive(Clone, Copy, Debug)]
pub struct MulticenterBondView<'a> {
    pub id: MulticenterBondId,
    atoms: &'a [NodeId],
    pub ast: &'a MulticenterBondAst,
    molecule: &'a MoleculeAst,
}

impl<'a> MulticenterBondView<'a> {
    #[inline]
    pub fn electrons(&self) -> &'a ElectronCountsAst {
        &self.ast.electrons
    }

    #[inline]
    pub fn charge(&self) -> &'a ValueAst {
        &self.ast.charge
    }

    #[inline]
    pub fn spin(&self) -> &'a SpinStateAst {
        &self.ast.spin
    }

    #[inline]
    pub fn constraints(&self) -> &'a MulticenterBondConstraints {
        &self.ast.constraints
    }

    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.atoms.iter().map(|&n| AtomId::from(n))
    }

    pub fn atoms(&self) -> impl Iterator<Item = AtomView<'a>> + 'a {
        let molecule = self.molecule;
        self.atoms
            .iter()
            .map(move |&n| molecule.atom(AtomId::from(n)))
    }

    /// Sum of per-atom electron contributions on this multicenter bond.
    /// `Lit(n)` when the counts are concrete; `Undetermined` otherwise.
    pub fn electron_count(&self) -> ValueAst {
        match &self.ast.electrons {
            ElectronCountsAst::Lit(counts) => ValueAst::Lit(counts.iter().sum()),
            ElectronCountsAst::Undetermined => ValueAst::Undetermined,
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

// Builder-scope view bundles for multicenter bonds.

pub struct MulticenterBondBuilderView<'a> {
    pub id: MulticenterBondId,
    pub ast: &'a MulticenterBondAst,
    pub(crate) atoms: &'a [NodeId],
}

impl<'a> MulticenterBondBuilderView<'a> {
    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.atoms.iter().map(|&n| AtomId::from(n))
    }
}

pub struct MulticenterBondBuilderViewMut<'a> {
    pub id: MulticenterBondId,
    pub ast: &'a mut MulticenterBondAst,
    pub(crate) atoms: &'a [NodeId],
}

impl<'a> MulticenterBondBuilderViewMut<'a> {
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
    use crate::ast::id::{AtomId, MulticenterBondId};
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
            Vec::new(),
            Vec::new(),
            Constraints::default(),
        )
    }

    #[rstest]
    fn test_multicenter_bond_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.multicenter_bonds().count(), 1);
    }

    #[rstest]
    fn test_multicenter_bond_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.multicenter_bonds().ids().collect::<Vec<_>>(),
            vec![MulticenterBondId(0)],
        );
    }

    #[rstest]
    fn test_multicenter_bond_views_iter(molecule: MoleculeAst) {
        let collected: Vec<(MulticenterBondId, Vec<AtomId>)> = molecule
            .multicenter_bonds()
            .iter()
            .map(|v| (v.id, v.atom_ids().collect()))
            .collect();
        assert_eq!(
            collected,
            vec![(MulticenterBondId(0), vec![AtomId(0), AtomId(1), AtomId(2)],)],
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
    fn test_multicenter_bond_views_index(molecule: MoleculeAst) {
        let _: &MulticenterBondAst = &molecule.multicenter_bonds()[MulticenterBondId(0)];
    }

    #[rstest]
    fn test_multicenter_bond_view_atom_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .multicenter_bond(MulticenterBondId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
        );
    }

    #[rstest]
    fn test_multicenter_bond_view_atoms(molecule: MoleculeAst) {
        let ids: Vec<AtomId> = molecule
            .multicenter_bond(MulticenterBondId(0))
            .atoms()
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, vec![AtomId(0), AtomId(1), AtomId(2)]);
    }

    #[rstest]
    fn test_multicenter_bond_view_electron_count(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .multicenter_bond(MulticenterBondId(0))
                .electron_count(),
            ValueAst::Undetermined,
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
}
