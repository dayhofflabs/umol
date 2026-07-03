//! One firing of a reaction — `apply`'s codomain.
//!
//! The two concrete molecule sides of a single rule application plus the correspondence between them:
//! `lhs` is the molecule the rule was matched into (borrowed), `rhs` the molecule produced (owned),
//! and `comap` maps `lhs` → `rhs` (preserved entities mated, deleted `lhs` entities left-exposed,
//! created entities right-exposed). It is the *instance* of a `ReactionAst` (rule : derivation ∷
//! function : one evaluation) and carries the ground-truth atom map — `apply` created the atoms, so
//! no post-hoc diff is needed to recover it; `to_reaction` abstracts back to the rule layer.
//!
//! Only `lhs` borrows; `rhs` / `comap` are owned, so `reverse` / `to_reaction` are self-contained
//! while the borrow lives, and persisting past it is a drop to the owned `(rhs, comap)`.

use umol_graph_core::{Correspondence, NodeId};

use super::correspondence::MoleculeCorrespondence;
use super::molecule::MoleculeAst;
use super::reaction::ReactionAst;

/// A reaction fired once at a match: its two concrete molecule sides (`lhs` ⇒ `rhs`) plus the
/// correspondence between them.
#[derive(Clone, Debug)]
pub struct ReactionDerivation<'a> {
    lhs: &'a MoleculeAst,
    rhs: MoleculeAst,
    comap: MoleculeCorrespondence,
}

impl<'a> ReactionDerivation<'a> {
    pub(crate) fn new(
        lhs: &'a MoleculeAst,
        rhs: MoleculeAst,
        comap: MoleculeCorrespondence,
    ) -> Self {
        Self { lhs, rhs, comap }
    }

    /// The molecule the rule was matched into.
    pub fn lhs(&self) -> &MoleculeAst {
        self.lhs
    }

    /// The molecule produced by the firing.
    pub fn rhs(&self) -> &MoleculeAst {
        &self.rhs
    }

    /// The `lhs`↔`rhs` correspondence: preserved entities mated, deleted `lhs` entities left-exposed,
    /// created entities right-exposed.
    pub fn comap(&self) -> &MoleculeCorrespondence {
        &self.comap
    }

    /// The atom-level slice of the comap — the per-step atom map.
    pub fn atom_map(&self) -> &Correspondence<NodeId> {
        self.comap.atoms()
    }

    /// Abstract back to the rule layer: `lhs` as the reaction's `lhs` plus the deltas taking it to
    /// `rhs` under the known comap. Inverse of `ReactionAst::apply`, up to delta normal form.
    pub fn to_reaction(&self) -> ReactionAst {
        ReactionAst::new(
            self.lhs.clone(),
            self.lhs.difference_to(&self.rhs, &self.comap),
        )
    }

    /// The reverse derivation `rhs ⇒ lhs`: sides swapped, comap inverted. Its `lhs` is this
    /// derivation's `rhs`, so the returned derivation borrows `self` (a shorter lifetime than the
    /// forward one).
    pub fn reverse(&self) -> ReactionDerivation<'_> {
        ReactionDerivation {
            lhs: &self.rhs,
            rhs: self.lhs.clone(),
            comap: self.comap.reverse(),
        }
    }

    /// Chain onto a following derivation `next` (which fires on this one's `rhs`): the composite
    /// `lhs ⇒ next.rhs`, with the comaps composed (pathway atom-map propagation).
    pub fn chain(&self, next: &ReactionDerivation) -> ReactionDerivation<'a> {
        ReactionDerivation {
            lhs: self.lhs,
            rhs: next.rhs.clone(),
            comap: self.comap.compose(&next.comap),
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;
    use umol_graph_core::NodeId;

    use super::super::atom::AtomAst;
    use super::super::bond::BondAst;
    use super::super::delta::{BondDelta, Delta, Deltas};
    use super::super::edit::BondFieldChange;
    use super::super::id::{AtomId, BondId};
    use super::super::value::ValueAst;
    use super::*;

    fn bond_order_molecule(order: u8) -> MoleculeAst {
        MoleculeAst::from_atoms_and_bonds(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(order))],
        )
    }

    fn total_atoms() -> Correspondence<NodeId> {
        Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))], 2, 2)
    }

    fn bond_order_modify(old: i64, new: i64) -> Deltas {
        Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
            id: BondId(0),
            change: BondFieldChange::Order {
                old: ValueAst::Lit(old),
                new: ValueAst::Lit(new),
            },
        })])
    }

    /// A `lhs ⇒ rhs` derivation over C-C, bond order 1 → 2, total atom map.
    #[fixture]
    fn derivation_parts() -> (MoleculeAst, MoleculeAst, MoleculeCorrespondence) {
        let lhs = bond_order_molecule(1);
        let rhs = bond_order_molecule(2);
        let comap = MoleculeCorrespondence::induce(&lhs, &rhs, total_atoms());
        (lhs, rhs, comap)
    }

    #[rstest]
    fn test_reaction_derivation_lhs(
        derivation_parts: (MoleculeAst, MoleculeAst, MoleculeCorrespondence),
    ) {
        let (lhs, rhs, comap) = derivation_parts;
        let derivation = ReactionDerivation::new(&lhs, rhs, comap);
        assert_eq!(derivation.lhs(), &bond_order_molecule(1));
    }

    #[rstest]
    fn test_reaction_derivation_rhs(
        derivation_parts: (MoleculeAst, MoleculeAst, MoleculeCorrespondence),
    ) {
        let (lhs, rhs, comap) = derivation_parts;
        let derivation = ReactionDerivation::new(&lhs, rhs, comap);
        assert_eq!(derivation.rhs(), &bond_order_molecule(2));
    }

    #[rstest]
    fn test_reaction_derivation_atom_map(
        derivation_parts: (MoleculeAst, MoleculeAst, MoleculeCorrespondence),
    ) {
        let (lhs, rhs, comap) = derivation_parts;
        let derivation = ReactionDerivation::new(&lhs, rhs, comap);
        assert_eq!(
            derivation.atom_map().mates(),
            &[(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))]
        );
    }

    #[rstest]
    fn test_reaction_derivation_to_reaction(
        derivation_parts: (MoleculeAst, MoleculeAst, MoleculeCorrespondence),
    ) {
        let (lhs, rhs, comap) = derivation_parts;
        let derivation = ReactionDerivation::new(&lhs, rhs, comap);
        assert_eq!(
            derivation.to_reaction(),
            ReactionAst::new(bond_order_molecule(1), bond_order_modify(1, 2))
        );
    }

    #[rstest]
    fn test_reaction_derivation_reverse(
        derivation_parts: (MoleculeAst, MoleculeAst, MoleculeCorrespondence),
    ) {
        let (lhs, rhs, comap) = derivation_parts;
        let derivation = ReactionDerivation::new(&lhs, rhs, comap);
        let reversed = derivation.reverse();
        assert_eq!(reversed.rhs(), &bond_order_molecule(1));
        assert_eq!(
            reversed.to_reaction(),
            ReactionAst::new(bond_order_molecule(2), bond_order_modify(2, 1))
        );
    }

    #[rstest]
    fn test_reaction_derivation_chain() {
        // lhs C-C(1) ⇒ mid C-C(2) ⇒ rhs C-C(3): the chain is lhs ⇒ rhs, a single order 1→3 modify.
        let lhs = bond_order_molecule(1);
        let mid = bond_order_molecule(2);
        let rhs = bond_order_molecule(3);
        let first = ReactionDerivation::new(
            &lhs,
            mid.clone(),
            MoleculeCorrespondence::induce(&lhs, &mid, total_atoms()),
        );
        let second = ReactionDerivation::new(
            &mid,
            rhs,
            MoleculeCorrespondence::induce(&mid, &bond_order_molecule(3), total_atoms()),
        );
        let chained = first.chain(&second);
        assert_eq!(chained.rhs(), &bond_order_molecule(3));
        assert_eq!(
            chained.to_reaction(),
            ReactionAst::new(bond_order_molecule(1), bond_order_modify(1, 3))
        );
    }
}
