//! Ring views over a molecule.

use std::borrow::Cow;

use super::super::constraint::RingScope;
use super::super::id::{AtomId, BondId};
use super::super::molecule::Molecule;
use super::super::num::NumForm;
use super::super::ring::{RingId, RingSet, RingSetKind};
use super::read_num;

/// Molecule ring views over owned or borrowed ring data.
///
/// The molecule is borrowed. The ring set must describe its covalent topology
/// and atom/bond ids for meaningful contextual readings. Mismatched context
/// does not panic but has no specified semantic result.
///
/// # Semantic properties
///
/// Owned and borrowed views over the same ring data and compatible molecule
/// give identical readings. Reconstructing a view after attribute-only mutation
/// preserves ring topology and reads the molecule's current attributes.
#[derive(Debug)]
pub struct RingViews<'a> {
    molecule: &'a Molecule,
    rings: Cow<'a, RingSet>,
}

impl<'a> RingViews<'a> {
    pub(crate) fn new(molecule: &'a Molecule, rings: RingSet) -> Self {
        Self {
            molecule,
            rings: Cow::Owned(rings),
        }
    }

    /// Borrow existing ring data and its molecule without copying or perceiving rings.
    ///
    /// The set must describe the molecule's covalent topology and atom/bond ids
    /// for meaningful readings. Reconstruction does not validate that relationship;
    /// mismatched context does not panic but has no specified semantic result.
    pub fn from_ring_set(rings: &'a RingSet, molecule: &'a Molecule) -> Self {
        Self {
            molecule,
            rings: Cow::Borrowed(rings),
        }
    }

    pub fn count(&self) -> usize {
        self.rings.count()
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = RingId> + '_ {
        self.rings.ids()
    }

    pub fn iter(&self) -> impl ExactSizeIterator<Item = RingView<'_>> + '_ {
        self.rings.iter()
    }

    pub fn get(&self, id: RingId) -> Option<RingView<'_>> {
        self.rings.get(id)
    }

    pub fn contains(&self, id: RingId) -> bool {
        self.rings.get(id).is_some()
    }

    /// The ring atom view of atom `id`.
    pub fn atom(&self, id: AtomId) -> RingAtomView<'_> {
        RingAtomView {
            rings: &self.rings,
            molecule: self.molecule,
            id,
        }
    }

    /// The ring bond view of bond `id`.
    pub fn bond(&self, id: BondId) -> RingBondView<'_> {
        RingBondView {
            rings: &self.rings,
            id,
        }
    }

    pub fn kind(&self) -> RingSetKind {
        self.rings.kind()
    }

    pub fn max_ring_size(&self) -> usize {
        self.rings.max_ring_size()
    }

    /// Consume this view, moving owned ring data or cloning borrowed data.
    ///
    /// The returned set preserves all ring ids, members, and relations and is
    /// independent of the molecule and any borrowed source set.
    pub fn into_ring_set(self) -> RingSet {
        self.rings.into_owned()
    }
}

/// Borrowed view of a single ring.
#[derive(Debug, Clone, Copy)]
pub struct RingView<'a> {
    pub id: RingId,
    atoms: &'a [AtomId],
    bonds: &'a [BondId],
}

impl<'a> RingView<'a> {
    pub(crate) fn new(id: RingId, atoms: &'a [AtomId], bonds: &'a [BondId]) -> Self {
        Self { id, atoms, bonds }
    }

    pub fn atoms(&self) -> &'a [AtomId] {
        self.atoms
    }

    pub fn bonds(&self) -> &'a [BondId] {
        self.bonds
    }

    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }
}

/// Ring atom data with reference to ring set and molecule.
#[derive(Debug, Clone, Copy)]
pub struct RingAtomView<'a> {
    rings: &'a RingSet,
    molecule: &'a Molecule,
    id: AtomId,
}

impl<'a> RingAtomView<'a> {
    /// Whether this atom lies in any ring.
    pub fn is_in_ring(&self) -> bool {
        self.rings.contains_atom(self.id)
    }

    /// The rings containing this atom.
    pub fn rings(&self) -> impl Iterator<Item = RingView<'a>> + 'a {
        let id = self.id;
        self.rings.iter().filter(move |v| v.atoms().contains(&id))
    }

    /// Count of rings in this view containing the atom and matching `scope` (`All` = any,
    /// `Size(s)` = size `s`). Constraint matching constructs the view with the fixed Relevant
    /// projection through size 22. Always `Lit`.
    pub fn ring_membership(&self, scope: RingScope) -> NumForm {
        atom_ring_membership(self.rings, self.id, scope)
    }

    pub fn ring_count(&self) -> NumForm {
        self.ring_membership(RingScope::All)
    }

    pub fn ring_size_count(&self, size: u8) -> NumForm {
        self.ring_membership(RingScope::Size(size))
    }

    /// Smallest containing ring size, or `None` if this atom is in no ring.
    pub fn smallest_ring_size(&self) -> Option<usize> {
        self.rings.atom_smallest_ring_size(self.id)
    }

    /// Count of incident bonds that lie in a ring. Always `Lit`.
    pub fn ring_degree(&self) -> NumForm {
        atom_ring_degree(self.molecule, self.rings, self.id)
    }

    /// Sum of bond orders of incident bonds that lie in a ring. `Undetermined`
    /// if a normalized order remains non-`Lit`; contradiction returns `None`.
    pub fn ring_valence(&self) -> Option<NumForm> {
        atom_ring_valence(self.molecule, self.rings, self.id)
    }
}

/// Ring bond data with reference to ring set.
#[derive(Debug, Clone, Copy)]
pub struct RingBondView<'a> {
    rings: &'a RingSet,
    id: BondId,
}

impl<'a> RingBondView<'a> {
    /// Whether this bond lies in any ring.
    pub fn is_in_ring(&self) -> bool {
        self.rings.contains_bond(self.id)
    }

    /// The rings containing this bond.
    pub fn rings(&self) -> impl Iterator<Item = RingView<'a>> + 'a {
        let id = self.id;
        self.rings.iter().filter(move |v| v.bonds().contains(&id))
    }

    /// Count of rings in this view containing the bond and matching `scope` (`All` = any,
    /// `Size(s)` = size `s`). Constraint matching constructs the view with the fixed Relevant
    /// projection through size 22. Always `Lit`.
    pub fn ring_membership(&self, scope: RingScope) -> NumForm {
        bond_ring_membership(self.rings, self.id, scope)
    }

    pub fn ring_count(&self) -> NumForm {
        self.ring_membership(RingScope::All)
    }

    pub fn ring_size_count(&self, s: u8) -> NumForm {
        self.ring_membership(RingScope::Size(s))
    }

    /// Smallest containing ring size, or `None` if this bond is in no ring.
    pub fn smallest_ring_size(&self) -> Option<usize> {
        self.rings.bond_smallest_ring_size(self.id)
    }
}

/// Count of rings containing atom `id` and matching `scope`. Always `Lit`.
pub(crate) fn atom_ring_membership(rings: &RingSet, id: AtomId, scope: RingScope) -> NumForm {
    let containing = rings.iter().filter(|v| v.atoms().contains(&id));
    let count = match scope {
        RingScope::All => containing.count(),
        RingScope::Size(s) => containing.filter(|r| r.len() == s as usize).count(),
    };
    NumForm::Lit(count as i64)
}

/// Count of incident bonds of atom `id` that lie in a ring. Always `Lit`.
pub(crate) fn atom_ring_degree(molecule: &Molecule, rings: &RingSet, id: AtomId) -> NumForm {
    if !molecule.raw_graph().contains_node(id.into()) {
        return NumForm::Lit(0);
    }
    let count = molecule
        .neighbors(id)
        .filter(|n| rings.contains_bond(n.bond().id()))
        .count();
    NumForm::Lit(count as i64)
}

/// Sum of bond orders of incident bonds of atom `id` that lie in a ring.
/// `Undetermined` if a normalized order remains non-`Lit`; contradiction returns `None`.
pub(crate) fn atom_ring_valence(
    molecule: &Molecule,
    rings: &RingSet,
    id: AtomId,
) -> Option<NumForm> {
    if !molecule.raw_graph().contains_node(id.into()) {
        return Some(NumForm::Lit(0));
    }
    molecule
        .neighbors(id)
        .filter(|n| rings.contains_bond(n.bond().id()))
        .try_fold(NumForm::Lit(0), |sum, n| {
            Some(sum + read_num(n.bond().order())?)
        })
}

/// Count of rings containing bond `id` and matching `scope`. Always `Lit`.
pub(crate) fn bond_ring_membership(rings: &RingSet, id: BondId, scope: RingScope) -> NumForm {
    let containing = rings.iter().filter(|v| v.bonds().contains(&id));
    let count = match scope {
        RingScope::All => containing.count(),
        RingScope::Size(s) => containing.filter(|r| r.len() == s as usize).count(),
    };
    NumForm::Lit(count as i64)
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;

    use super::super::assert_exact_size_by;
    use super::*;
    use crate::ir::atom::AtomForm;
    use crate::ir::bond::BondForm;
    use crate::ir::molecule::MoleculeEntries;
    use crate::ir::ring::{RingConfig, RingConnection, RingModel, RingRelation};

    // A six-membered ring (atoms 0-5, bonds 0-5) with a pendant chain atom 6
    // (bond 6, atom 0 to atom 6).
    #[fixture]
    fn ring_with_chain() -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 7],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
                (AtomId(3), AtomId(4), BondForm::from_order(1)),
                (AtomId(4), AtomId(5), BondForm::from_order(1)),
                (AtomId(5), AtomId(0), BondForm::from_order(1)),
                (AtomId(0), AtomId(6), BondForm::from_order(1)),
            ],
            ..Default::default()
        })
    }

    #[fixture]
    fn fused_rings() -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 4],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
                (AtomId(2), AtomId(0), BondForm::from_order(1)),
                (AtomId(1), AtomId(3), BondForm::from_order(1)),
                (AtomId(3), AtomId(2), BondForm::from_order(1)),
            ],
            ..Default::default()
        })
    }

    #[rstest]
    #[case::owned(false)]
    #[case::borrowed(true)]
    fn test_ring_views_from_ring_set(fused_rings: Molecule, #[case] borrowed: bool) {
        let ring_set = fused_rings
            .rings(RingModel::default(), RingConfig::default())
            .into_ring_set();
        let views = if borrowed {
            RingViews::from_ring_set(&ring_set, &fused_rings)
        } else {
            fused_rings.rings(RingModel::default(), RingConfig::default())
        };
        assert_eq!(views.kind(), RingSetKind::Relevant);
        assert_eq!(views.max_ring_size(), 22);
        assert_eq!(views.count(), 2);
        assert_eq!(views.ids().collect::<Vec<_>>(), vec![RingId(0), RingId(1)]);
        assert!(views.contains(RingId(0)));
        assert!(!views.contains(RingId(2)));
        assert_eq!(views.get(RingId(2)).map(|ring| ring.id), None);
        let members: Vec<_> = views
            .iter()
            .map(|ring| {
                let mut atoms = ring.atoms().to_vec();
                let mut bonds = ring.bonds().to_vec();
                atoms.sort_unstable();
                bonds.sort_unstable();
                (ring.id, atoms, bonds)
            })
            .collect();
        assert_eq!(
            members,
            vec![
                (
                    RingId(0),
                    vec![AtomId(0), AtomId(1), AtomId(2)],
                    vec![BondId(0), BondId(1), BondId(2)]
                ),
                (
                    RingId(1),
                    vec![AtomId(1), AtomId(2), AtomId(3)],
                    vec![BondId(1), BondId(3), BondId(4)]
                ),
            ]
        );
        assert_eq!(views.atom(AtomId(1)).ring_count(), NumForm::Lit(2));
        assert_eq!(views.atom(AtomId(1)).smallest_ring_size(), Some(3));
        assert_eq!(views.atom(AtomId(1)).ring_degree(), NumForm::Lit(3));
        assert_eq!(views.atom(AtomId(1)).ring_valence(), Some(NumForm::Lit(3)));
        assert_eq!(views.bond(BondId(1)).ring_count(), NumForm::Lit(2));
        assert_eq!(views.bond(BondId(1)).smallest_ring_size(), Some(3));
        if borrowed {
            assert_eq!(
                views.get(RingId(0)).unwrap().atoms().as_ptr(),
                ring_set.get(RingId(0)).unwrap().atoms().as_ptr()
            );
            assert_eq!(
                views.get(RingId(0)).unwrap().bonds().as_ptr(),
                ring_set.get(RingId(0)).unwrap().bonds().as_ptr()
            );
        }
    }

    #[rstest]
    #[case::literal(NumForm::Lit(2), Some(NumForm::Lit(3)))]
    #[case::unknown(NumForm::Undetermined, Some(NumForm::Undetermined))]
    #[case::singleton(NumForm::lit_set([2]), Some(NumForm::Lit(3)))]
    #[case::bottom(NumForm::lit_set([]), None)]
    fn test_ring_views_from_ring_set_mutation(
        mut ring_with_chain: Molecule,
        #[case] order: NumForm,
        #[case] expected: Option<NumForm>,
    ) {
        let ring_set = ring_with_chain
            .rings(RingModel::default(), RingConfig::default())
            .into_ring_set();
        assert_eq!(
            RingViews::from_ring_set(&ring_set, &ring_with_chain)
                .atom(AtomId(0))
                .ring_valence(),
            Some(NumForm::Lit(2))
        );
        ring_with_chain.bond_mut(BondId(0)).attributes_mut().order = order.clone();
        let views = RingViews::from_ring_set(&ring_set, &ring_with_chain);
        let atom = views.atom(AtomId(0));
        assert_eq!(atom.ring_valence(), expected);
        assert_eq!(atom.ring_degree(), NumForm::Lit(2));
        assert_eq!(atom.ring_count(), NumForm::Lit(1));
        assert_eq!(atom.ring_size_count(6), NumForm::Lit(1));
        assert_eq!(atom.smallest_ring_size(), Some(6));
        assert_eq!(
            atom.rings().map(|ring| ring.id).collect::<Vec<_>>(),
            vec![RingId(0)]
        );
        assert_eq!(views.bond(BondId(0)).ring_count(), NumForm::Lit(1));
        assert_eq!(ring_with_chain.bond(BondId(0)).order(), &order);
    }

    #[rstest]
    #[case::empty(0, vec![])]
    #[case::smaller(2, vec![(0, 1)])]
    #[case::different_topology(7, vec![(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6)])]
    #[case::larger(8, vec![(0, 7)])]
    fn test_ring_views_from_ring_set_context(
        ring_with_chain: Molecule,
        #[case] atom_count: usize,
        #[case] bonds: Vec<(u32, u32)>,
    ) {
        let ring_set = ring_with_chain
            .rings(RingModel::default(), RingConfig::default())
            .into_ring_set();
        let context = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); atom_count],
            bonds: bonds
                .into_iter()
                .map(|(a, b)| (AtomId(a), AtomId(b), BondForm::from_order(1)))
                .collect(),
            ..Default::default()
        });
        let views = RingViews::from_ring_set(&ring_set, &context);
        // Mismatched contextual readings promise only panic freedom.
        for id in ring_with_chain
            .atoms()
            .ids()
            .chain([AtomId(atom_count as u32)])
        {
            let atom = views.atom(id);
            let _ = (atom.ring_degree(), atom.ring_valence());
        }
        let mut atoms = views.get(RingId(0)).unwrap().atoms().to_vec();
        atoms.sort_unstable();
        assert_eq!(
            atoms,
            vec![
                AtomId(0),
                AtomId(1),
                AtomId(2),
                AtomId(3),
                AtomId(4),
                AtomId(5)
            ]
        );
    }

    #[rstest]
    fn test_ring_views_count(ring_with_chain: Molecule) {
        assert_eq!(
            ring_with_chain
                .rings(RingModel::default(), RingConfig::default())
                .count(),
            1
        );
    }

    #[rstest]
    fn test_ring_views_ids(ring_with_chain: Molecule) {
        let empty = Molecule::default();
        assert_exact_size_by(
            empty
                .rings(RingModel::default(), RingConfig::default())
                .ids(),
            vec![],
            |id| id,
        );
        assert_exact_size_by(
            ring_with_chain
                .rings(RingModel::default(), RingConfig::default())
                .ids(),
            vec![RingId(0)],
            |id| id,
        );
    }

    #[rstest]
    fn test_ring_views_iter(ring_with_chain: Molecule) {
        let empty = Molecule::default();
        let empty_rings = empty.rings(RingModel::default(), RingConfig::default());
        assert_exact_size_by(empty_rings.iter(), vec![], |ring| (ring.id, ring.len()));

        let rings = ring_with_chain.rings(RingModel::default(), RingConfig::default());
        assert_exact_size_by(rings.iter(), vec![(RingId(0), 6)], |ring| {
            (ring.id, ring.len())
        });
    }

    #[rstest]
    #[case::present(RingId(0), Some(6))]
    #[case::absent(RingId(1), None)]
    fn test_ring_views_get(
        ring_with_chain: Molecule,
        #[case] ring: RingId,
        #[case] expected_len: Option<usize>,
    ) {
        assert_eq!(
            ring_with_chain
                .rings(RingModel::default(), RingConfig::default())
                .get(ring)
                .map(|r| r.len()),
            expected_len
        );
    }

    #[rstest]
    #[case::present(RingId(0), true)]
    #[case::absent(RingId(1), false)]
    fn test_ring_views_contains(
        ring_with_chain: Molecule,
        #[case] ring: RingId,
        #[case] expected: bool,
    ) {
        assert_eq!(
            ring_with_chain
                .rings(RingModel::default(), RingConfig::default())
                .contains(ring),
            expected
        );
    }

    #[rstest]
    fn test_ring_views_kind(ring_with_chain: Molecule) {
        assert_eq!(
            ring_with_chain
                .rings(RingModel::default(), RingConfig::default())
                .kind(),
            RingSetKind::Relevant
        );
    }

    #[rstest]
    fn test_ring_views_max_ring_size(ring_with_chain: Molecule) {
        assert_eq!(
            ring_with_chain
                .rings(RingModel::default(), RingConfig::default())
                .max_ring_size(),
            22
        );
    }

    #[rstest]
    fn test_ring_views_into_ring_set(ring_with_chain: Molecule) {
        let views = ring_with_chain.rings(RingModel::default(), RingConfig::default());
        let atoms = views.get(RingId(0)).unwrap().atoms().as_ptr();
        let bonds = views.get(RingId(0)).unwrap().bonds().as_ptr();
        let ring_set = views.into_ring_set();
        drop(ring_with_chain);
        assert_eq!(ring_set.count(), 1);
        assert!(ring_set.contains_atom(AtomId(0)));
        assert!(!ring_set.contains_atom(AtomId(6)));
        assert_eq!(ring_set.get(RingId(0)).unwrap().atoms().as_ptr(), atoms);
        assert_eq!(ring_set.get(RingId(0)).unwrap().bonds().as_ptr(), bonds);
        assert_eq!(ring_set.atom_smallest_ring_size(AtomId(0)), Some(6));
        assert_eq!(ring_set.bond_smallest_ring_size(BondId(0)), Some(6));
        assert_eq!(ring_set.graph().edges(), &[]);
    }

    #[rstest]
    fn test_ring_views_into_ring_set_borrowed(fused_rings: Molecule) {
        let source = fused_rings
            .rings(RingModel::default(), RingConfig::default())
            .into_ring_set();
        let ring_set = RingViews::from_ring_set(&source, &fused_rings).into_ring_set();
        assert_ne!(
            ring_set.get(RingId(0)).unwrap().atoms().as_ptr(),
            source.get(RingId(0)).unwrap().atoms().as_ptr()
        );
        assert_ne!(
            ring_set.get(RingId(0)).unwrap().bonds().as_ptr(),
            source.get(RingId(0)).unwrap().bonds().as_ptr()
        );
        drop(source);
        drop(fused_rings);
        assert_eq!(ring_set.kind(), RingSetKind::Relevant);
        assert_eq!(ring_set.max_ring_size(), 22);
        assert_eq!(ring_set.count(), 2);
        assert_eq!(
            ring_set.shared_atoms(RingId(0), RingId(1)),
            vec![AtomId(1), AtomId(2)]
        );
        assert_eq!(ring_set.shared_bonds(RingId(0), RingId(1)), vec![BondId(1)]);
        assert_eq!(ring_set.relation(RingId(0), RingId(1)), RingRelation::Fused);
        assert_eq!(ring_set.atom_smallest_ring_size(AtomId(1)), Some(3));
        assert_eq!(ring_set.bond_smallest_ring_size(BondId(1)), Some(3));
        assert_eq!(
            ring_set.fused_components(),
            vec![vec![RingId(0), RingId(1)]]
        );
        assert_eq!(
            ring_set.graph().edges(),
            &[RingConnection {
                first: RingId(0),
                second: RingId(1),
                relation: RingRelation::Fused,
            }]
        );
    }

    #[rstest]
    #[case::ring_atom(AtomId(0), true)]
    #[case::chain_atom(AtomId(6), false)]
    fn test_ring_atom_view_is_in_ring(
        ring_with_chain: Molecule,
        #[case] atom: AtomId,
        #[case] expected: bool,
    ) {
        assert_eq!(
            ring_with_chain
                .rings(RingModel::default(), RingConfig::default())
                .atom(atom)
                .is_in_ring(),
            expected
        );
    }

    #[rstest]
    #[case::ring_atom(AtomId(0), vec![RingId(0)])]
    #[case::chain_atom(AtomId(6), vec![])]
    fn test_ring_atom_view_rings(
        ring_with_chain: Molecule,
        #[case] atom: AtomId,
        #[case] expected: Vec<RingId>,
    ) {
        let ids: Vec<RingId> = ring_with_chain
            .rings(RingModel::default(), RingConfig::default())
            .atom(atom)
            .rings()
            .map(|r| r.id)
            .collect();
        assert_eq!(ids, expected);
    }

    #[rstest]
    #[case::all_ring_atom(AtomId(0), RingScope::All, NumForm::Lit(1))]
    #[case::size_match(AtomId(0), RingScope::Size(6), NumForm::Lit(1))]
    #[case::size_no_match(AtomId(0), RingScope::Size(5), NumForm::Lit(0))]
    #[case::chain_atom(AtomId(6), RingScope::All, NumForm::Lit(0))]
    fn test_ring_atom_view_ring_membership(
        ring_with_chain: Molecule,
        #[case] atom: AtomId,
        #[case] scope: RingScope,
        #[case] expected: NumForm,
    ) {
        assert_eq!(
            ring_with_chain
                .rings(RingModel::default(), RingConfig::default())
                .atom(atom)
                .ring_membership(scope),
            expected
        );
    }

    #[rstest]
    #[case::ring_atom(AtomId(0), NumForm::Lit(1))]
    #[case::chain_atom(AtomId(6), NumForm::Lit(0))]
    fn test_ring_atom_view_ring_count(
        ring_with_chain: Molecule,
        #[case] atom: AtomId,
        #[case] expected: NumForm,
    ) {
        assert_eq!(
            ring_with_chain
                .rings(RingModel::default(), RingConfig::default())
                .atom(atom)
                .ring_count(),
            expected
        );
    }

    #[rstest]
    #[case::size_match(AtomId(0), 6, NumForm::Lit(1))]
    #[case::size_no_match(AtomId(0), 5, NumForm::Lit(0))]
    #[case::chain_atom(AtomId(6), 6, NumForm::Lit(0))]
    fn test_ring_atom_view_ring_size_count(
        ring_with_chain: Molecule,
        #[case] atom: AtomId,
        #[case] size: u8,
        #[case] expected: NumForm,
    ) {
        assert_eq!(
            ring_with_chain
                .rings(RingModel::default(), RingConfig::default())
                .atom(atom)
                .ring_size_count(size),
            expected
        );
    }

    #[rstest]
    #[case::ring_atom(AtomId(0), Some(6))]
    #[case::chain_atom(AtomId(6), None)]
    fn test_ring_atom_view_smallest_ring_size(
        ring_with_chain: Molecule,
        #[case] atom: AtomId,
        #[case] expected: Option<usize>,
    ) {
        assert_eq!(
            ring_with_chain
                .rings(RingModel::default(), RingConfig::default())
                .atom(atom)
                .smallest_ring_size(),
            expected
        );
    }

    #[rstest]
    #[case::ring_atom(AtomId(0), NumForm::Lit(2))]
    #[case::chain_atom(AtomId(6), NumForm::Lit(0))]
    fn test_ring_atom_view_ring_degree(
        ring_with_chain: Molecule,
        #[case] atom: AtomId,
        #[case] expected: NumForm,
    ) {
        assert_eq!(
            ring_with_chain
                .rings(RingModel::default(), RingConfig::default())
                .atom(atom)
                .ring_degree(),
            expected
        );
    }

    #[rstest]
    #[case::ring_atom(AtomId(0), NumForm::Lit(2))]
    #[case::chain_atom(AtomId(6), NumForm::Lit(0))]
    fn test_ring_atom_view_ring_valence(
        ring_with_chain: Molecule,
        #[case] atom: AtomId,
        #[case] expected: NumForm,
    ) {
        assert_eq!(
            ring_with_chain
                .rings(RingModel::default(), RingConfig::default())
                .atom(atom)
                .ring_valence(),
            Some(expected)
        );
    }

    #[rstest]
    #[case::ring_bond(BondId(0), true)]
    #[case::chain_bond(BondId(6), false)]
    fn test_ring_bond_view_is_in_ring(
        ring_with_chain: Molecule,
        #[case] bond: BondId,
        #[case] expected: bool,
    ) {
        assert_eq!(
            ring_with_chain
                .rings(RingModel::default(), RingConfig::default())
                .bond(bond)
                .is_in_ring(),
            expected
        );
    }

    #[rstest]
    #[case::ring_bond(BondId(0), vec![RingId(0)])]
    #[case::chain_bond(BondId(6), vec![])]
    fn test_ring_bond_view_rings(
        ring_with_chain: Molecule,
        #[case] bond: BondId,
        #[case] expected: Vec<RingId>,
    ) {
        let ids: Vec<RingId> = ring_with_chain
            .rings(RingModel::default(), RingConfig::default())
            .bond(bond)
            .rings()
            .map(|r| r.id)
            .collect();
        assert_eq!(ids, expected);
    }

    #[rstest]
    #[case::all_ring_bond(BondId(0), RingScope::All, NumForm::Lit(1))]
    #[case::size_match(BondId(0), RingScope::Size(6), NumForm::Lit(1))]
    #[case::size_no_match(BondId(0), RingScope::Size(5), NumForm::Lit(0))]
    #[case::chain_bond(BondId(6), RingScope::All, NumForm::Lit(0))]
    fn test_ring_bond_view_ring_membership(
        ring_with_chain: Molecule,
        #[case] bond: BondId,
        #[case] scope: RingScope,
        #[case] expected: NumForm,
    ) {
        assert_eq!(
            ring_with_chain
                .rings(RingModel::default(), RingConfig::default())
                .bond(bond)
                .ring_membership(scope),
            expected
        );
    }

    #[rstest]
    #[case::ring_bond(BondId(0), NumForm::Lit(1))]
    #[case::chain_bond(BondId(6), NumForm::Lit(0))]
    fn test_ring_bond_view_ring_count(
        ring_with_chain: Molecule,
        #[case] bond: BondId,
        #[case] expected: NumForm,
    ) {
        assert_eq!(
            ring_with_chain
                .rings(RingModel::default(), RingConfig::default())
                .bond(bond)
                .ring_count(),
            expected
        );
    }

    #[rstest]
    #[case::size_match(BondId(0), 6, NumForm::Lit(1))]
    #[case::size_no_match(BondId(0), 5, NumForm::Lit(0))]
    #[case::chain_bond(BondId(6), 6, NumForm::Lit(0))]
    fn test_ring_bond_view_ring_size_count(
        ring_with_chain: Molecule,
        #[case] bond: BondId,
        #[case] size: u8,
        #[case] expected: NumForm,
    ) {
        assert_eq!(
            ring_with_chain
                .rings(RingModel::default(), RingConfig::default())
                .bond(bond)
                .ring_size_count(size),
            expected
        );
    }

    #[rstest]
    #[case::ring_bond(BondId(0), Some(6))]
    #[case::chain_bond(BondId(6), None)]
    fn test_ring_bond_view_smallest_ring_size(
        ring_with_chain: Molecule,
        #[case] bond: BondId,
        #[case] expected: Option<usize>,
    ) {
        assert_eq!(
            ring_with_chain
                .rings(RingModel::default(), RingConfig::default())
                .bond(bond)
                .smallest_ring_size(),
            expected
        );
    }
}
