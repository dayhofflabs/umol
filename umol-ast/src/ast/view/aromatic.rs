//! Aromatic system views.

use std::collections::HashSet;
use std::ops::Index;

use umol_graph_core::{NodeId, RelationId, Unordered, VarRelationSet};

use super::super::aromatic::AromaticSystemAst;
use super::super::constraint::AromaticSystemConstraints;
use super::super::electrons::ElectronCountsAst;
use super::super::id::{AromaticSystemId, AtomId, BondId};
use super::super::molecule::MoleculeAst;
use super::super::ring::RingView;
use super::super::spin::SpinStateAst;
use super::super::traits::Lattice;
use super::super::value::ValueAst;
use super::atom::AtomView;
use super::bond::BondView;

/// Namespace accessor for aromatic-system views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct AromaticSystemViews<'a> {
    molecule: &'a MoleculeAst,
    aromatic_systems: &'a VarRelationSet<NodeId, Unordered, AromaticSystemAst>,
}

impl<'a> AromaticSystemViews<'a> {
    pub(crate) fn new(
        molecule: &'a MoleculeAst,
        aromatic_systems: &'a VarRelationSet<NodeId, Unordered, AromaticSystemAst>,
    ) -> Self {
        Self {
            molecule,
            aromatic_systems,
        }
    }

    pub fn count(&self) -> usize {
        self.aromatic_systems.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = AromaticSystemId> {
        self.aromatic_systems
            .relation_ids()
            .map(AromaticSystemId::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = AromaticSystemView<'a>> {
        let molecule = self.molecule;
        let set = self.aromatic_systems;
        set.relation_ids().map(move |rid| AromaticSystemView {
            id: AromaticSystemId::from(rid),
            ast: set.data(rid),
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
            ast: self.aromatic_systems.data(rid),
            atoms: self.aromatic_systems.participants(rid),
            molecule: self.molecule,
        })
    }

    /// Ids of aromatic systems incident on `atom`.
    pub fn incident_ids(&self, atom: AtomId) -> impl Iterator<Item = AromaticSystemId> + 'a {
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
    pub fn incident(&self, atom: AtomId) -> impl Iterator<Item = AromaticSystemView<'a>> + 'a {
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
    pub fn connecting_id(
        &self,
        atoms: impl IntoIterator<Item = AtomId>,
    ) -> Option<AromaticSystemId> {
        let target: HashSet<AtomId> = atoms.into_iter().collect();
        let &first = target.iter().next()?;
        self.incident_ids(first).find(|&id| {
            let parts: HashSet<AtomId> = self
                .aromatic_systems
                .participants(RelationId::from(id))
                .iter()
                .map(|&n| AtomId::from(n))
                .collect();
            parts == target
        })
    }

    /// View of the aromatic system whose atom set equals `atoms`, if any.
    pub fn connecting(
        &self,
        atoms: impl IntoIterator<Item = AtomId>,
    ) -> Option<AromaticSystemView<'a>> {
        self.connecting_id(atoms).map(|id| {
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

impl<'a> Index<AromaticSystemId> for AromaticSystemViews<'a> {
    type Output = AromaticSystemAst;
    fn index(&self, id: AromaticSystemId) -> &AromaticSystemAst {
        self.aromatic_systems.data(RelationId::from(id))
    }
}

/// Borrowed view of an aromatic system: its index, the `AromaticSystemAst`,
/// and accessors for member atoms and induced ring bonds via `atoms()` and
/// `bonds()`.
#[derive(Clone, Copy, Debug)]
pub struct AromaticSystemView<'a> {
    pub id: AromaticSystemId,
    atoms: &'a [NodeId],
    pub ast: &'a AromaticSystemAst,
    molecule: &'a MoleculeAst,
}

impl<'a> AromaticSystemView<'a> {
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
    pub fn constraints(&self) -> &'a AromaticSystemConstraints {
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

    /// Sum of per-atom electron contributions on this aromatic system.
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

    /// Rings from the molecule's canonical `RingSet` that share at least
    /// one atom with this aromatic system.
    pub fn overlapping_rings(&self) -> impl Iterator<Item = RingView<'a>> + 'a {
        let atoms: Vec<AtomId> = self.atoms.iter().map(|&n| AtomId::from(n)).collect();
        self.molecule
            .rings()
            .iter()
            .filter(move |r| r.atoms().iter().any(|a| atoms.contains(a)))
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

// Builder-scope view bundles for aromatic systems.

pub struct AromaticSystemBuilderView<'a> {
    pub id: AromaticSystemId,
    pub ast: &'a AromaticSystemAst,
    pub(crate) atoms: &'a [NodeId],
}

impl<'a> AromaticSystemBuilderView<'a> {
    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.atoms.iter().map(|&n| AtomId::from(n))
    }
}

pub struct AromaticSystemBuilderViewMut<'a> {
    pub id: AromaticSystemId,
    pub ast: &'a mut AromaticSystemAst,
    pub(crate) atoms: &'a [NodeId],
}

impl<'a> AromaticSystemBuilderViewMut<'a> {
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
    use crate::ast::id::{AromaticSystemId, AtomId, BondId};
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
    fn test_aromatic_system_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.aromatic_systems().count(), 1);
    }

    #[rstest]
    fn test_aromatic_system_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.aromatic_systems().ids().collect::<Vec<_>>(),
            vec![AromaticSystemId(0)],
        );
    }

    #[rstest]
    fn test_aromatic_system_views_iter(molecule: MoleculeAst) {
        let collected: Vec<(AromaticSystemId, Vec<AtomId>)> = molecule
            .aromatic_systems()
            .iter()
            .map(|v| (v.id, v.atom_ids().collect()))
            .collect();
        assert_eq!(
            collected,
            vec![(AromaticSystemId(0), vec![AtomId(0), AtomId(1), AtomId(2)])],
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
    fn test_aromatic_system_views_index(molecule: MoleculeAst) {
        let _: &AromaticSystemAst = &molecule.aromatic_systems()[AromaticSystemId(0)];
    }

    #[rstest]
    fn test_aromatic_system_view_atom_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .aromatic_system(AromaticSystemId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
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
    fn test_aromatic_system_view_atoms(molecule: MoleculeAst) {
        let ids: Vec<AtomId> = molecule
            .aromatic_system(AromaticSystemId(0))
            .atoms()
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, vec![AtomId(0), AtomId(1), AtomId(2)]);
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
    fn test_aromatic_system_view_electron_count(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .aromatic_system(AromaticSystemId(0))
                .electron_count(),
            ValueAst::Undetermined,
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
    fn test_aromatic_system_view_overlapping_rings(molecule: MoleculeAst) {
        let ids: Vec<usize> = molecule
            .aromatic_system(AromaticSystemId(0))
            .overlapping_rings()
            .map(|r| r.len())
            .collect();
        assert_eq!(ids, Vec::<usize>::new());
    }
}
