//! Ring views over a molecule.

use super::super::constraint::RingScope;
use super::super::id::{AtomId, BondId};
use super::super::molecule::Molecule;
use super::super::num::NumForm;
use super::super::ring::{RingId, RingSet, RingSetKind};

/// Molecule ring views: owned ring set plus borrow of molecule.
#[derive(Debug)]
pub struct RingViews<'a> {
    molecule: &'a Molecule,
    rings: RingSet,
}

impl<'a> RingViews<'a> {
    pub(crate) fn new(molecule: &'a Molecule, rings: RingSet) -> Self {
        Self { molecule, rings }
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

    /// The ring situation of `atom`.
    pub fn atom(&self, id: AtomId) -> RingAtomView<'_> {
        RingAtomView {
            rings: &self.rings,
            molecule: self.molecule,
            id,
        }
    }

    /// The ring situation of `bond`.
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

    /// Consume this view, returning its owned `RingSet`.
    pub fn into_ring_set(self) -> RingSet {
        self.rings
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

    pub fn ring_size_count(&self, s: u8) -> NumForm {
        self.ring_membership(RingScope::Size(s))
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
    /// if any contributing bond's order is non-`Lit`.
    pub fn ring_valence(&self) -> NumForm {
        atom_ring_valence(self.molecule, self.rings, self.id)
    }
}

// Derivation layer beneath the ring-atom facades: per-quantity functions of the
// ring set, molecule, and atom id, presented by `RingAtomView` (typed) and
// `AtomConstraintsView` (keyed).

/// Count of rings containing `atom` and matching `scope`. Always `Lit`.
pub(crate) fn atom_ring_membership(rings: &RingSet, atom: AtomId, scope: RingScope) -> NumForm {
    let containing = rings.iter().filter(|v| v.atoms().contains(&atom));
    let count = match scope {
        RingScope::All => containing.count(),
        RingScope::Size(s) => containing.filter(|r| r.len() == s as usize).count(),
    };
    NumForm::Lit(count as i64)
}

/// Count of incident bonds of `atom` that lie in a ring. Always `Lit`.
pub(crate) fn atom_ring_degree(molecule: &Molecule, rings: &RingSet, atom: AtomId) -> NumForm {
    let count = molecule
        .neighbors(atom)
        .filter(|n| rings.contains_bond(n.bond().id))
        .count();
    NumForm::Lit(count as i64)
}

/// Sum of bond orders of incident bonds of `atom` that lie in a ring.
/// `Undetermined` if any contributing bond's order is non-`Lit`.
pub(crate) fn atom_ring_valence(molecule: &Molecule, rings: &RingSet, atom: AtomId) -> NumForm {
    molecule
        .neighbors(atom)
        .filter(|n| rings.contains_bond(n.bond().id))
        .map(|n| n.bond().order().clone())
        .fold(NumForm::Lit(0), |acc, order| acc + order)
}

/// Count of rings containing `bond` and matching `scope`. Always `Lit`.
pub(crate) fn bond_ring_membership(rings: &RingSet, bond: BondId, scope: RingScope) -> NumForm {
    let containing = rings.iter().filter(|v| v.bonds().contains(&bond));
    let count = match scope {
        RingScope::All => containing.count(),
        RingScope::Size(s) => containing.filter(|r| r.len() == s as usize).count(),
    };
    NumForm::Lit(count as i64)
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

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;

    use super::super::assert_exact_size_by;
    use super::*;
    use crate::ir::atom::AtomForm;
    use crate::ir::bond::BondForm;
    use crate::ir::molecule::MoleculeEntries;
    use crate::ir::ring::{RingConfig, RingModel};

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
        let ring_set = ring_with_chain
            .rings(RingModel::default(), RingConfig::default())
            .into_ring_set();
        assert_eq!(ring_set.count(), 1);
        assert!(ring_set.contains_atom(AtomId(0)));
        assert!(!ring_set.contains_atom(AtomId(6)));
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
            expected
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
