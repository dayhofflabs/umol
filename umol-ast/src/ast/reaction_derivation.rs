//! One firing of a reaction — `apply`'s codomain.
//!
//! The two concrete molecule sides of a single rule application plus the correspondence between them:
//! `lhs` is an owned snapshot of the molecule the rule matched, `rhs` is the molecule produced, and
//! `comap` maps `lhs` → `rhs` (preserved entities mated, deleted `lhs` entities left-exposed, created
//! entities right-exposed). It is the *instance* of a `ReactionAst` (rule : derivation ∷ function :
//! one evaluation) and carries the ground-truth atom map — `apply` created the atoms, so no post-hoc
//! diff is needed to recover it; `to_reaction` abstracts back to the rule layer.

use umol_graph_core::{Correspondence, NodeId};

use super::correspondence::MoleculeCorrespondence;
use super::molecule::MoleculeAst;
#[cfg(test)]
use super::molecule::MoleculeParts;
use super::reaction::ReactionAst;

/// A reaction fired once at a match: its two concrete molecule sides (`lhs` ⇒ `rhs`) plus the
/// correspondence between them.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReactionDerivation {
    lhs: MoleculeAst,
    rhs: MoleculeAst,
    comap: MoleculeCorrespondence,
}

impl ReactionDerivation {
    pub(crate) fn new(lhs: MoleculeAst, rhs: MoleculeAst, comap: MoleculeCorrespondence) -> Self {
        Self { lhs, rhs, comap }
    }

    /// The molecule the rule was matched into.
    pub fn lhs(&self) -> &MoleculeAst {
        &self.lhs
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

    /// The reverse derivation `rhs ⇒ lhs`: sides swapped, comap inverted.
    pub fn reverse(&self) -> ReactionDerivation {
        ReactionDerivation {
            lhs: self.rhs.clone(),
            rhs: self.lhs.clone(),
            comap: self.comap.reverse(),
        }
    }

    /// Chain onto a following derivation `next` (which fires on this one's `rhs`): the composite
    /// `lhs ⇒ next.rhs`, with the comaps composed (pathway atom-map propagation).
    pub fn chain(&self, next: &ReactionDerivation) -> ReactionDerivation {
        ReactionDerivation {
            lhs: self.lhs.clone(),
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

    /// A `lhs ⇒ rhs` derivation over C-C, bond order 1 → 2, total atom map.
    #[fixture]
    fn derivation_parts() -> (MoleculeAst, MoleculeAst, MoleculeCorrespondence) {
        let lhs = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ..Default::default()
        });
        let rhs = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
            ..Default::default()
        });
        let comap = MoleculeCorrespondence::induce(
            &lhs,
            &rhs,
            Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))], 2, 2),
        );
        (lhs, rhs, comap)
    }

    #[rstest]
    fn test_reaction_derivation_new(
        derivation_parts: (MoleculeAst, MoleculeAst, MoleculeCorrespondence),
    ) {
        let (lhs, rhs, comap) = derivation_parts;
        let expected = ReactionDerivation {
            lhs: lhs.clone(),
            rhs: rhs.clone(),
            comap: comap.clone(),
        };
        assert_eq!(ReactionDerivation::new(lhs, rhs, comap), expected);
    }

    #[rstest]
    fn test_reaction_derivation_new_independence(
        derivation_parts: (MoleculeAst, MoleculeAst, MoleculeCorrespondence),
    ) {
        let (mut lhs, rhs, comap) = derivation_parts;
        let expected = lhs.clone();
        let derivation = ReactionDerivation::new(lhs.clone(), rhs, comap);
        *lhs.atom_mut(AtomId(0)).ast = AtomAst::from_element(Element::N);
        assert_eq!(derivation.lhs(), &expected);
    }

    #[rstest]
    fn test_reaction_derivation_lhs(
        derivation_parts: (MoleculeAst, MoleculeAst, MoleculeCorrespondence),
    ) {
        let (lhs, rhs, comap) = derivation_parts;
        let expected = lhs.clone();
        let derivation = ReactionDerivation::new(lhs, rhs, comap);
        assert_eq!(derivation.lhs(), &expected);
    }

    #[rstest]
    fn test_reaction_derivation_rhs(
        derivation_parts: (MoleculeAst, MoleculeAst, MoleculeCorrespondence),
    ) {
        let (lhs, rhs, comap) = derivation_parts;
        let expected = rhs.clone();
        let derivation = ReactionDerivation::new(lhs, rhs, comap);
        assert_eq!(derivation.rhs(), &expected);
    }

    #[rstest]
    fn test_reaction_derivation_comap(
        derivation_parts: (MoleculeAst, MoleculeAst, MoleculeCorrespondence),
    ) {
        let (lhs, rhs, comap) = derivation_parts;
        let expected = comap.clone();
        let derivation = ReactionDerivation::new(lhs, rhs, comap);
        assert_eq!(derivation.comap(), &expected);
    }

    #[rstest]
    fn test_reaction_derivation_atom_map(
        derivation_parts: (MoleculeAst, MoleculeAst, MoleculeCorrespondence),
    ) {
        let (lhs, rhs, comap) = derivation_parts;
        let derivation = ReactionDerivation::new(lhs, rhs, comap);
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
        let expected_lhs = lhs.clone();
        let derivation = ReactionDerivation::new(lhs, rhs, comap);
        assert_eq!(
            derivation.to_reaction(),
            ReactionAst::new(
                expected_lhs,
                Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order {
                        old: ValueAst::Lit(1),
                        new: ValueAst::Lit(2),
                    },
                })]),
            )
        );
    }

    #[rstest]
    fn test_reaction_derivation_reverse(
        derivation_parts: (MoleculeAst, MoleculeAst, MoleculeCorrespondence),
    ) {
        let (lhs, rhs, comap) = derivation_parts;
        let expected = ReactionDerivation {
            lhs: rhs.clone(),
            rhs: lhs.clone(),
            comap: comap.reverse(),
        };
        assert_eq!(ReactionDerivation::new(lhs, rhs, comap).reverse(), expected);
    }

    #[rstest]
    fn test_reaction_derivation_chain() {
        let lhs = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ..Default::default()
        });
        let mid = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
            ..Default::default()
        });
        let rhs = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(3))],
            ..Default::default()
        });
        let first_comap = MoleculeCorrespondence::induce(
            &lhs,
            &mid,
            Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))], 2, 2),
        );
        let second_comap = MoleculeCorrespondence::induce(
            &mid,
            &rhs,
            Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))], 2, 2),
        );
        let first = ReactionDerivation::new(lhs.clone(), mid.clone(), first_comap.clone());
        let second = ReactionDerivation::new(mid, rhs.clone(), second_comap.clone());
        assert_eq!(
            first.chain(&second),
            ReactionDerivation {
                lhs,
                rhs,
                comap: first_comap.compose(&second_comap),
            }
        );
    }
}
