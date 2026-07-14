//! Ring views over a molecule.
//!
//! [`RingsView`] is a molecule's rings — an owned canonical `RingSet` plus a
//! borrow of the molecule it was enumerated from — with the ring-collection
//! surface (`count` / `ids` / `iter` / `get` / `contains`) mirroring the entity
//! namespace accessors, and per-atom / per-bond ring sub-views carrying the
//! topology-derived ring queries. Built by [`MoleculeAst::rings_view`]. The
//! borrow means a ring view cannot be held across a structural mutation of the
//! molecule.

use super::super::constraint::RingScope;
use super::super::id::{AtomId, BondId};
use super::super::molecule::MoleculeAst;
use super::super::ring::{RingFamily, RingId, RingSet, RingView};
use super::super::value::ValueAst;

/// The rings of a molecule: an owned canonical `RingSet` plus a borrow of the
/// molecule. Shape mirrors `AtomsView` / `BondsView`.
#[derive(Debug)]
pub struct RingsView<'a> {
    molecule: &'a MoleculeAst,
    rings: RingSet,
}

impl<'a> RingsView<'a> {
    pub(crate) fn new(molecule: &'a MoleculeAst, rings: RingSet) -> Self {
        Self { molecule, rings }
    }

    pub fn count(&self) -> usize {
        self.rings.count()
    }

    pub fn ids(&self) -> impl Iterator<Item = RingId> + '_ {
        self.rings.ids()
    }

    pub fn iter(&self) -> impl Iterator<Item = RingView<'_>> + '_ {
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

    pub fn family(&self) -> RingFamily {
        self.rings.family()
    }

    pub fn max_ring_size(&self) -> usize {
        self.rings.max_ring_size()
    }
}

/// One atom's ring situation within a [`RingsView`].
#[derive(Debug, Clone, Copy)]
pub struct RingAtomView<'a> {
    rings: &'a RingSet,
    molecule: &'a MoleculeAst,
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

    /// Count of rings containing this atom matching `scope` (`All` = any,
    /// `Size(s)` = size `s`). Always `Lit`.
    pub fn ring_membership(&self, scope: RingScope) -> ValueAst {
        let count = match scope {
            RingScope::All => self.rings().count(),
            RingScope::Size(s) => self.rings().filter(|r| r.len() == s as usize).count(),
        };
        ValueAst::Lit(count as i64)
    }

    pub fn ring_count(&self) -> ValueAst {
        self.ring_membership(RingScope::All)
    }

    pub fn ring_size_count(&self, s: u8) -> ValueAst {
        self.ring_membership(RingScope::Size(s))
    }

    /// Smallest containing ring size, or `None` if this atom is in no ring.
    pub fn smallest_ring_size(&self) -> Option<usize> {
        self.rings.atom_smallest_ring_size(self.id)
    }

    /// Count of incident bonds that lie in a ring. Always `Lit`.
    pub fn ring_degree(&self) -> ValueAst {
        let count = self
            .molecule
            .atom(self.id)
            .neighbors()
            .filter(|n| self.rings.contains_bond(n.bond().id))
            .count();
        ValueAst::Lit(count as i64)
    }

    /// Sum of bond orders of incident bonds that lie in a ring. `Undetermined`
    /// if any contributing bond's order is non-`Lit`.
    pub fn ring_valence(&self) -> ValueAst {
        self.molecule
            .atom(self.id)
            .neighbors()
            .filter(|n| self.rings.contains_bond(n.bond().id))
            .map(|n| n.bond().order().clone())
            .fold(ValueAst::Lit(0), |acc, order| acc + order)
    }
}

/// One bond's ring situation within a [`RingsView`].
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

    /// Count of rings containing this bond matching `scope` (`All` = any,
    /// `Size(s)` = size `s`). Always `Lit`.
    pub fn ring_membership(&self, scope: RingScope) -> ValueAst {
        let count = match scope {
            RingScope::All => self.rings().count(),
            RingScope::Size(s) => self.rings().filter(|r| r.len() == s as usize).count(),
        };
        ValueAst::Lit(count as i64)
    }

    pub fn ring_count(&self) -> ValueAst {
        self.ring_membership(RingScope::All)
    }

    pub fn ring_size_count(&self, s: u8) -> ValueAst {
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

    use super::*;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::molecule::MoleculeParts;

    // A six-membered ring (atoms 0-5, bonds 0-5) with a pendant chain atom 6
    // (bond 6, atom 0 to atom 6).
    #[fixture]
    fn ring_with_chain() -> MoleculeAst {
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C); 7],
            bonds: vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(1)),
                (AtomId(2), AtomId(3), BondAst::from_order(1)),
                (AtomId(3), AtomId(4), BondAst::from_order(1)),
                (AtomId(4), AtomId(5), BondAst::from_order(1)),
                (AtomId(5), AtomId(0), BondAst::from_order(1)),
                (AtomId(0), AtomId(6), BondAst::from_order(1)),
            ],
            ..Default::default()
        })
    }

    #[rstest]
    fn test_rings_view_count(ring_with_chain: MoleculeAst) {
        let rings = ring_with_chain.rings_view();
        assert_eq!(rings.count(), 1);
        assert_eq!(rings.family(), RingFamily::Relevant);
        assert_eq!(rings.max_ring_size(), 22);
    }

    #[rstest]
    fn test_rings_view_get_contains(ring_with_chain: MoleculeAst) {
        let rings = ring_with_chain.rings_view();
        assert!(rings.contains(RingId(0)));
        assert!(!rings.contains(RingId(1)));
        assert_eq!(rings.get(RingId(0)).unwrap().len(), 6);
        assert!(rings.get(RingId(1)).is_none());
        assert_eq!(rings.ids().collect::<Vec<_>>(), vec![RingId(0)]);
        assert_eq!(rings.iter().count(), 1);
    }

    #[rstest]
    #[case::ring_atom(AtomId(0), true)]
    #[case::chain_atom(AtomId(6), false)]
    fn test_ring_atom_view_is_in_ring(
        ring_with_chain: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: bool,
    ) {
        assert_eq!(
            ring_with_chain.rings_view().atom(atom).is_in_ring(),
            expected
        );
    }

    #[rstest]
    fn test_ring_atom_view_counts(ring_with_chain: MoleculeAst) {
        let rings = ring_with_chain.rings_view();
        let atom = rings.atom(AtomId(0));
        assert_eq!(atom.ring_count(), ValueAst::Lit(1));
        assert_eq!(atom.ring_size_count(6), ValueAst::Lit(1));
        assert_eq!(atom.ring_size_count(5), ValueAst::Lit(0));
        assert_eq!(atom.smallest_ring_size(), Some(6));
        // atom 0 has two incident ring bonds (0-1, 0-5) and one chain bond (0-6)
        assert_eq!(atom.ring_degree(), ValueAst::Lit(2));
        assert_eq!(atom.ring_valence(), ValueAst::Lit(2));
        assert_eq!(rings.atom(AtomId(6)).smallest_ring_size(), None);
        assert_eq!(rings.atom(AtomId(6)).ring_count(), ValueAst::Lit(0));
    }

    #[rstest]
    #[case::ring_bond(BondId(0), true)]
    #[case::chain_bond(BondId(6), false)]
    fn test_ring_bond_view_is_in_ring(
        ring_with_chain: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected: bool,
    ) {
        assert_eq!(
            ring_with_chain.rings_view().bond(bond).is_in_ring(),
            expected
        );
    }

    #[rstest]
    fn test_ring_bond_view_counts(ring_with_chain: MoleculeAst) {
        let rings = ring_with_chain.rings_view();
        assert_eq!(rings.bond(BondId(0)).ring_count(), ValueAst::Lit(1));
        assert_eq!(rings.bond(BondId(0)).ring_size_count(6), ValueAst::Lit(1));
        assert_eq!(rings.bond(BondId(0)).smallest_ring_size(), Some(6));
        assert_eq!(rings.bond(BondId(6)).ring_count(), ValueAst::Lit(0));
        assert_eq!(rings.bond(BondId(6)).smallest_ring_size(), None);
    }
}
