//! One firing of a reaction — `apply`'s codomain.
//!
//! The DPO direct derivation `host ⇒ product` reduced to its externally useful data: the produced
//! molecule plus the host↔product correspondence `apply` builds while firing (preserved entities
//! mated, deleted host entities left-exposed, created entities right-exposed). It is the *instance*
//! of a `ReactionAst` (rule : derivation ∷ function : one evaluation) and carries the ground-truth
//! atom map — `apply` created the atoms, so no post-hoc diff is needed to recover it; `to_reaction`
//! abstracts back to the rule layer.
//!
//! Only the host borrows; `product` / `comap` are owned, so `reverse` / `to_reaction` are
//! self-contained while the borrow lives, and persisting past the host is a drop to the owned
//! `(product, comap)`.

use umol_graph_core::{Correspondence, NodeId};

use super::correspondence::MoleculeCorrespondence;
use super::molecule::MoleculeAst;
use super::reaction::ReactionAst;

/// A reaction fired once at a host: `host ⇒ product` with the correspondence between the two.
#[derive(Clone, Debug)]
pub struct ReactionDerivation<'a> {
    host: &'a MoleculeAst,
    product: MoleculeAst,
    comap: MoleculeCorrespondence,
}

impl<'a> ReactionDerivation<'a> {
    pub(crate) fn new(
        host: &'a MoleculeAst,
        product: MoleculeAst,
        comap: MoleculeCorrespondence,
    ) -> Self {
        Self {
            host,
            product,
            comap,
        }
    }

    /// The molecule produced by the firing.
    pub fn product(&self) -> &MoleculeAst {
        &self.product
    }

    /// The host↔product correspondence: preserved entities mated, deleted host entities left-exposed,
    /// created entities right-exposed.
    pub fn comap(&self) -> &MoleculeCorrespondence {
        &self.comap
    }

    /// The atom-level slice of the comap — the per-step atom map.
    pub fn atom_map(&self) -> &Correspondence<NodeId> {
        self.comap.atoms()
    }

    /// Abstract back to the rule layer: the host as `lhs` plus the deltas taking it to `product`
    /// under the known comap. Inverse of `ReactionAst::apply`, up to delta normal form.
    pub fn to_reaction(&self) -> ReactionAst {
        ReactionAst::new(
            self.host.clone(),
            self.host.difference_to(&self.product, &self.comap),
        )
    }

    /// The reverse derivation `product ⇒ host`: sides swapped, comap inverted. Its host is this
    /// derivation's product, so the returned derivation borrows `self` (a shorter lifetime than the
    /// forward one).
    pub fn reverse(&self) -> ReactionDerivation<'_> {
        ReactionDerivation {
            host: &self.product,
            product: self.host.clone(),
            comap: self.comap.reverse(),
        }
    }

    /// Chain onto a following derivation `next` (which fires on this one's product): the composite
    /// `host ⇒ next.product`, with the comaps composed (pathway atom-map propagation).
    pub fn chain(&self, next: &ReactionDerivation) -> ReactionDerivation<'a> {
        ReactionDerivation {
            host: self.host,
            product: next.product.clone(),
            comap: self.comap.compose(&next.comap),
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;

    use super::super::atom::AtomAst;
    use super::super::bond::BondAst;
    use super::super::delta::{BondDelta, Delta, Deltas};
    use super::super::edit::BondFieldChange;
    use super::super::id::{AtomId, BondId};
    use super::super::value::ValueAst;
    use super::*;
    use umol_graph_core::NodeId;

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

    /// A `host ⇒ product` derivation over C-C, bond order `left` → `right`, total atom map.
    #[fixture]
    fn derivation_parts() -> (MoleculeAst, MoleculeAst, MoleculeCorrespondence) {
        let host = bond_order_molecule(1);
        let product = bond_order_molecule(2);
        let comap = MoleculeCorrespondence::induce(&host, &product, total_atoms());
        (host, product, comap)
    }

    #[rstest]
    fn test_reaction_derivation_product(
        derivation_parts: (MoleculeAst, MoleculeAst, MoleculeCorrespondence),
    ) {
        let (host, product, comap) = derivation_parts;
        let derivation = ReactionDerivation::new(&host, product, comap);
        assert_eq!(derivation.product(), &bond_order_molecule(2));
    }

    #[rstest]
    fn test_reaction_derivation_atom_map(
        derivation_parts: (MoleculeAst, MoleculeAst, MoleculeCorrespondence),
    ) {
        let (host, product, comap) = derivation_parts;
        let derivation = ReactionDerivation::new(&host, product, comap);
        assert_eq!(
            derivation.atom_map().mates(),
            &[(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))]
        );
    }

    #[rstest]
    fn test_reaction_derivation_to_reaction(
        derivation_parts: (MoleculeAst, MoleculeAst, MoleculeCorrespondence),
    ) {
        let (host, product, comap) = derivation_parts;
        let derivation = ReactionDerivation::new(&host, product, comap);
        assert_eq!(
            derivation.to_reaction(),
            ReactionAst::new(bond_order_molecule(1), bond_order_modify(1, 2))
        );
    }

    #[rstest]
    fn test_reaction_derivation_reverse(
        derivation_parts: (MoleculeAst, MoleculeAst, MoleculeCorrespondence),
    ) {
        let (host, product, comap) = derivation_parts;
        let derivation = ReactionDerivation::new(&host, product, comap);
        let reversed = derivation.reverse();
        assert_eq!(reversed.product(), &bond_order_molecule(1));
        assert_eq!(
            reversed.to_reaction(),
            ReactionAst::new(bond_order_molecule(2), bond_order_modify(2, 1))
        );
    }

    #[rstest]
    fn test_reaction_derivation_chain() {
        // host C-C(1) ⇒ mid C-C(2) ⇒ end C-C(3): the chain is host ⇒ end, a single order 1→3 modify.
        let host = bond_order_molecule(1);
        let mid = bond_order_molecule(2);
        let end = bond_order_molecule(3);
        let first = ReactionDerivation::new(
            &host,
            mid.clone(),
            MoleculeCorrespondence::induce(&host, &mid, total_atoms()),
        );
        let second = ReactionDerivation::new(
            &mid,
            end,
            MoleculeCorrespondence::induce(&mid, &bond_order_molecule(3), total_atoms()),
        );
        let chained = first.chain(&second);
        assert_eq!(chained.product(), &bond_order_molecule(3));
        assert_eq!(
            chained.to_reaction(),
            ReactionAst::new(bond_order_molecule(1), bond_order_modify(1, 3))
        );
    }
}
