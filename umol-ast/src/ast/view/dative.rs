//! Dative bond views.

use std::collections::{HashMap, HashSet};
use std::iter;
use std::ops::Index;

use umol_graph_core::{FixedVarBirelationSet, NodeId, Ordered, RelationId, Unordered};

use super::super::constraint::DativeBondConstraintsAst;
use super::super::dative::DativeBondAst;
use super::super::id::{AtomId, DativeBondId};
use super::super::molecule::MoleculeAst;
use super::super::traits::Lattice;
use super::super::value::ValueAst;
use super::atom::AtomView;

/// Namespace accessor for dative-bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct DativeBondViews<'a> {
    molecule: &'a MoleculeAst,
    dative_bonds: &'a FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, DativeBondAst>,
}

impl<'a> DativeBondViews<'a> {
    pub(crate) fn new(
        molecule: &'a MoleculeAst,
        dative_bonds: &'a FixedVarBirelationSet<
            NodeId,
            Ordered,
            1,
            NodeId,
            Unordered,
            DativeBondAst,
        >,
    ) -> Self {
        Self {
            molecule,
            dative_bonds,
        }
    }

    pub fn count(&self) -> usize {
        self.dative_bonds.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = DativeBondId> {
        self.dative_bonds.relation_ids().map(DativeBondId::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = DativeBondView<'a>> {
        let molecule = self.molecule;
        let set = self.dative_bonds;
        set.relation_ids().map(move |rid| DativeBondView {
            id: DativeBondId::from(rid),
            ast: set.data(rid),
            acceptor_id: set.participants_1(rid)[0],
            donors: set.participants_2(rid),
            molecule,
        })
    }

    /// Whether a dative bond names its acceptor as a donor, repeats a donor, or two bonds share an
    /// acceptor and a donor (parallel dative). Emit-compliance peer of
    /// [`StereoAtomViews::has_conflict`].
    pub fn has_conflict(&self) -> bool {
        let mut donors_by_acceptor: HashMap<AtomId, HashSet<AtomId>> = HashMap::new();
        for view in self.iter() {
            let acceptor = view.acceptor_id();
            let mut donors: HashSet<AtomId> = HashSet::new();
            for donor in view.donor_ids() {
                if donor == acceptor || !donors.insert(donor) {
                    return true;
                }
            }
            let accumulated = donors_by_acceptor.entry(acceptor).or_default();
            for donor in donors {
                if !accumulated.insert(donor) {
                    return true;
                }
            }
        }
        false
    }

    pub fn contains(&self, id: DativeBondId) -> bool {
        self.dative_bonds.contains(RelationId::from(id))
    }

    pub fn get(&self, id: DativeBondId) -> Option<DativeBondView<'a>> {
        if !self.contains(id) {
            return None;
        }
        let rid = RelationId::from(id);
        Some(DativeBondView {
            id,
            ast: self.dative_bonds.data(rid),
            acceptor_id: self.dative_bonds.participants_1(rid)[0],
            donors: self.dative_bonds.participants_2(rid),
            molecule: self.molecule,
        })
    }

    /// Ids of dative bonds incident on `atom`.
    pub fn incident_ids(&self, atom: AtomId) -> impl Iterator<Item = DativeBondId> + 'a {
        self.dative_bonds
            .incident(NodeId::from(atom))
            .iter()
            .map(|&rid| DativeBondId::from(rid))
    }

    /// Whether any dative bond is incident on `atom`.
    pub fn has_incident(&self, atom: AtomId) -> bool {
        self.dative_bonds.has_incident(NodeId::from(atom))
    }

    /// Views of dative bonds incident on `atom`.
    pub fn incident(&self, atom: AtomId) -> impl Iterator<Item = DativeBondView<'a>> + 'a {
        let molecule = self.molecule;
        let set = self.dative_bonds;
        self.incident_ids(atom).map(move |id| {
            let rid = RelationId::from(id);
            DativeBondView {
                id,
                ast: set.data(rid),
                acceptor_id: set.participants_1(rid)[0],
                donors: set.participants_2(rid),
                molecule,
            }
        })
    }

    /// Id of the dative bond with exactly this acceptor and donor set, if any. Per-factor: the
    /// donor/acceptor roles are matched, not the merged atom set.
    pub fn connecting_id(&self, acceptor: AtomId, donors: &[AtomId]) -> Option<DativeBondId> {
        let donor_nodes: Vec<NodeId> = donors.iter().map(|&a| NodeId::from(a)).collect();
        self.dative_bonds
            .find_by_participants(&[NodeId::from(acceptor)], &donor_nodes)
            .map(DativeBondId::from)
    }

    /// View of the dative bond with exactly this acceptor and donor set, if any.
    pub fn connecting(&self, acceptor: AtomId, donors: &[AtomId]) -> Option<DativeBondView<'a>> {
        self.connecting_id(acceptor, donors).map(|id| {
            self.get(id).expect(
                "dative bond id from relation set must refer to a dative bond in this molecule",
            )
        })
    }

    /// Ids of dative bonds whose participants all lie in `atoms`.
    pub fn induced_ids(&self, atoms: &[AtomId]) -> Vec<DativeBondId> {
        let set: HashSet<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        self.dative_bonds
            .relation_ids()
            .filter(|&rid| {
                self.dative_bonds
                    .participants_1(rid)
                    .iter()
                    .chain(self.dative_bonds.participants_2(rid))
                    .all(|p| set.contains(p))
            })
            .map(DativeBondId::from)
            .collect()
    }

    /// Views of dative bonds whose participants all lie in `atoms`.
    pub fn induced(&self, atoms: &[AtomId]) -> Vec<DativeBondView<'a>> {
        self.induced_ids(atoms)
            .into_iter()
            .map(|id| {
                self.get(id).expect(
                    "dative bond id from relation set must refer to a dative bond in this molecule",
                )
            })
            .collect()
    }
}

impl<'a> Index<DativeBondId> for DativeBondViews<'a> {
    type Output = DativeBondAst;
    fn index(&self, id: DativeBondId) -> &DativeBondAst {
        self.dative_bonds.data(RelationId::from(id))
    }
}

/// Borrowed view of a dative bond: index, the designated acceptor atom,
/// and underlying `DativeBondAst`. Donor atoms via `donors()` / `donor_ids()`;
/// the full participant set (donors then acceptor) via `atoms()` / `atom_ids()`.
#[derive(Clone, Copy, Debug)]
pub struct DativeBondView<'a> {
    pub id: DativeBondId,
    acceptor_id: NodeId,
    donors: &'a [NodeId],
    pub ast: &'a DativeBondAst,
    molecule: &'a MoleculeAst,
}

impl<'a> DativeBondView<'a> {
    #[inline]
    pub fn order(&self) -> &'a ValueAst {
        &self.ast.order
    }

    #[inline]
    pub fn constraints(&self) -> &'a DativeBondConstraintsAst {
        &self.ast.constraints
    }

    /// Donor atom ids.
    pub fn donor_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.donors.iter().map(|&n| AtomId::from(n))
    }

    pub fn acceptor_id(&self) -> AtomId {
        AtomId::from(self.acceptor_id)
    }

    /// All atoms in this dative bond: the donors followed by the acceptor.
    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.donors
            .iter()
            .copied()
            .chain(iter::once(self.acceptor_id))
            .map(AtomId::from)
    }

    /// Donor atom views.
    pub fn donors(&self) -> impl Iterator<Item = AtomView<'a>> + 'a {
        let molecule = self.molecule;
        self.donor_ids().map(move |id| molecule.atom(id))
    }

    /// View of the acceptor atom.
    pub fn acceptor(&self) -> AtomView<'a> {
        self.molecule.atom(self.acceptor_id())
    }

    /// Views of all atoms in this dative bond (donors then acceptor).
    pub fn atoms(&self) -> impl Iterator<Item = AtomView<'a>> + 'a {
        let molecule = self.molecule;
        self.atom_ids().map(move |id| molecule.atom(id))
    }

    pub fn donor_count(&self) -> usize {
        self.donors.len()
    }

    pub fn atom_count(&self) -> usize {
        self.donor_count() + 1
    }

    /// Is dative bond ground
    pub fn is_ground(&self) -> bool {
        self.ast.is_ground()
    }

    /// Is dative bond undetermined
    pub fn is_undetermined(&self) -> bool {
        self.ast.is_undetermined()
    }
}

/// Mutable borrowed view of a dative bond: its id, participants (donors +
/// acceptor, owned) and mutable data. Molecule-scope peer of `DativeBondView`.
#[derive(Debug)]
pub struct DativeBondViewMut<'a> {
    pub id: DativeBondId,
    acceptor: AtomId,
    donors: Vec<AtomId>,
    pub ast: &'a mut DativeBondAst,
}

impl<'a> DativeBondViewMut<'a> {
    /// The acceptor atom id.
    pub fn acceptor_id(&self) -> AtomId {
        self.acceptor
    }

    /// The donor atom ids.
    pub fn donor_ids(&self) -> impl Iterator<Item = AtomId> + '_ {
        self.donors.iter().copied()
    }

    /// All atoms: the donors followed by the acceptor.
    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + '_ {
        self.donors.iter().copied().chain(iter::once(self.acceptor))
    }
}

// Builder-scope view bundles for dative bonds.

pub struct DativeBondEditorView<'a> {
    pub id: DativeBondId,
    pub ast: &'a DativeBondAst,
    pub(crate) donors: &'a [NodeId],
    pub acceptor_id: AtomId,
}

impl<'a> DativeBondEditorView<'a> {
    /// All atoms: donors followed by the acceptor.
    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        let acceptor = self.acceptor_id;
        self.donors
            .iter()
            .map(|&n| AtomId::from(n))
            .chain(iter::once(acceptor))
    }
}

pub struct DativeBondEditorViewMut<'a> {
    pub id: DativeBondId,
    pub ast: &'a mut DativeBondAst,
    pub(crate) donors: &'a [NodeId],
    pub acceptor_id: AtomId,
}

impl<'a> DativeBondEditorViewMut<'a> {
    /// All atoms: donors followed by the acceptor.
    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + '_ {
        let acceptor = self.acceptor_id;
        self.donors
            .iter()
            .map(|&n| AtomId::from(n))
            .chain(iter::once(acceptor))
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;

    use crate::ast::aromatic::AromaticSystemAst;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::dative::DativeBondAst;
    use crate::ast::id::{AtomId, DativeBondId};
    use crate::ast::molecule::{MoleculeAst, MoleculeParts};
    use crate::ast::multicenter::MulticenterBondAst;
    use crate::ast::noncovalent::{NoncovalentBondAst, NoncovalentBondKind};
    use crate::ast::value::ValueAst;

    #[fixture]
    fn molecule() -> MoleculeAst {
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(2)),
                (AtomId(2), AtomId(3), BondAst::from_order(1)),
            ],
            dative: vec![(vec![AtomId(2)], AtomId(3), DativeBondAst::from_order(1))],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemAst::default(),
            )],
            multicenter: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                MulticenterBondAst::default(),
            )],
            noncovalent: vec![(
                AtomId(0),
                AtomId(3),
                NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            ..Default::default()
        })
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
            .map(|v| (v.id, v.acceptor_id(), v.ast.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![(DativeBondId(0), AtomId(3), DativeBondAst::from_order(1))],
        );
    }

    #[rstest]
    #[case::present(DativeBondId(0), true)]
    #[case::absent(DativeBondId(99), false)]
    fn test_dative_bond_views_contains(
        molecule: MoleculeAst,
        #[case] id: DativeBondId,
        #[case] expected: bool,
    ) {
        assert_eq!(molecule.dative_bonds().contains(id), expected);
    }

    #[rstest]
    fn test_dative_bond_views_get(molecule: MoleculeAst) {
        let res = molecule.dative_bonds().get(DativeBondId(0));
        assert!(res.is_some());
        let view = res.unwrap();
        assert_eq!(view.id, DativeBondId(0));
        assert_eq!(view.acceptor_id(), AtomId(3));
    }

    #[rstest]
    fn test_dative_bond_views_get_none(molecule: MoleculeAst) {
        let res = molecule.dative_bonds().get(DativeBondId(99));
        assert!(res.is_none());
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
        assert_eq!(
            molecule.dative_bond(DativeBondId(0)).acceptor_id(),
            AtomId(3)
        );
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
