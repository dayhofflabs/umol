//! Reaction fingerprints: a molecular featurizer applied to each side of a
//! reaction, then combined across roles. The product side is derived from the
//! reactant (`lhs`) and the `deltas` via `to_reaction_span().right()`. `Difference`
//! computes counts difference (product minus reactant); `DisjointUnion` side-tags
//! each feature and unions both sides.

use umol_ast::ast::ReactionAst;

use super::feature_set::{FeatureSet, SignedFeatureSet};
use super::featurizer::{Featurizer, FingerprintError};

/// Which side of a reaction a feature came from.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Side {
    Reactant,
    Product,
}

/// How the two sides' feature sets are combined into a reaction fingerprint.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReactionCombinator {
    Difference,
    DisjointUnion,
}

/// A reaction fingerprint, shaped by the [`ReactionCombinator`] that produced it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ReactionFingerprint {
    Difference(SignedFeatureSet<u64>),
    DisjointUnion(FeatureSet<(Side, u64)>),
}

/// Featurize `reaction` by applying `featurizer` to the reactant (`lhs`) and the derived
/// product (`to_reaction_span().right()`), then combining per `combinator`.
/// `Inconsistent` if the deltas cannot be resolved to a product.
pub fn featurize_reaction(
    reaction: &ReactionAst,
    featurizer: &Featurizer,
    combinator: ReactionCombinator,
) -> Result<ReactionFingerprint, FingerprintError> {
    let product = reaction
        .to_reaction_span()
        .map_err(|_| FingerprintError::Inconsistent)?
        .right();
    Ok(match combinator {
        ReactionCombinator::Difference => {
            let reactants = featurizer.featurize_counted(&reaction.lhs)?;
            let products = featurizer.featurize_counted(&product)?;
            ReactionFingerprint::Difference(SignedFeatureSet::difference(&products, &reactants))
        }
        ReactionCombinator::DisjointUnion => {
            let reactants = featurizer.featurize(&reaction.lhs)?;
            let products = featurizer.featurize(&product)?;
            let tagged = reactants
                .ids()
                .iter()
                .map(|&id| (Side::Reactant, id))
                .chain(products.ids().iter().map(|&id| (Side::Product, id)));
            ReactionFingerprint::DisjointUnion(FeatureSet::from_features(tagged))
        }
    })
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::{AtomDelta, AtomId, BondDelta, BondId, Delta, Deltas, ReactionAst};

    use super::*;
    use crate::fingerprint::{Featurizer, MorganFeaturizer};
    use crate::parse::parse_smiles;

    /// Morgan radius-0 oxygen invariant of ethanol — present in `CCO`, absent in `CC`.
    const ETHANOL_OXYGEN: u64 = 864662311;

    // Identity reaction (no deltas): every feature count cancels, so the difference is empty.
    #[rstest]
    fn test_featurize_reaction_difference_identity() {
        let reaction = ReactionAst::new(parse_smiles("CCO").unwrap(), Deltas::new());
        let featurizer = Featurizer::Morgan(MorganFeaturizer::new(1));
        let fingerprint =
            featurize_reaction(&reaction, &featurizer, ReactionCombinator::Difference).unwrap();
        match fingerprint {
            ReactionFingerprint::Difference(difference) => assert!(difference.is_empty()),
            other => panic!("expected Difference, got {other:?}"),
        }
    }

    // CCO with the oxygen (atom 2) and its bond removed: the oxygen feature is
    // products(0) − reactants(1) = −1.
    #[rstest]
    fn test_featurize_reaction_difference() {
        let lhs = parse_smiles("CCO").unwrap();
        let oxygen = lhs.atom(AtomId(2)).ast.clone();
        let bond = lhs.bond(BondId(1)).ast.clone();
        let reaction = ReactionAst::new(
            lhs,
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Remove {
                    id: AtomId(2),
                    ast: oxygen,
                }),
                Delta::Bond(BondDelta::Remove {
                    id: BondId(1),
                    atoms: [AtomId(1), AtomId(2)],
                    ast: bond,
                }),
            ]),
        );
        let featurizer = Featurizer::Morgan(MorganFeaturizer::new(1));
        let fingerprint =
            featurize_reaction(&reaction, &featurizer, ReactionCombinator::Difference).unwrap();
        match fingerprint {
            ReactionFingerprint::Difference(difference) => {
                assert_eq!(difference.count(&ETHANOL_OXYGEN), -1);
            }
            other => panic!("expected Difference, got {other:?}"),
        }
    }

    // DisjointUnion tags each side and never cancels: an identity reaction's set is
    // exactly the molecule's features once as Reactant and once as Product.
    #[rstest]
    fn test_featurize_reaction_disjoint_union() {
        let reaction = ReactionAst::new(parse_smiles("CCO").unwrap(), Deltas::new());
        let featurizer = Featurizer::Morgan(MorganFeaturizer::new(1));
        let single = featurizer.featurize(&parse_smiles("CCO").unwrap()).unwrap();
        let fingerprint =
            featurize_reaction(&reaction, &featurizer, ReactionCombinator::DisjointUnion).unwrap();
        match fingerprint {
            ReactionFingerprint::DisjointUnion(union) => {
                assert_eq!(union.len(), 2 * single.len());
                assert!(union.ids().contains(&(Side::Reactant, ETHANOL_OXYGEN)));
                assert!(union.ids().contains(&(Side::Product, ETHANOL_OXYGEN)));
            }
            other => panic!("expected DisjointUnion, got {other:?}"),
        }
    }
}
