//! Reaction fingerprints: a molecular featurizer applied to each side of a
//! reaction, then combined across roles. `Difference` is the RDKit reaction
//! difference fingerprint (product counts minus reactant counts); `DisjointUnion`
//! side-tags each feature and unions both sides (no cancellation). Neither uses
//! the atom map. Both sides must be ground.

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
    /// Signed product-minus-reactant counts (RDKit reaction difference FP).
    Difference,
    /// Side-tagged union of both sides' binary features; no cancellation, so an
    /// unchanged scaffold is retained on both sides.
    DisjointUnion,
}

/// A reaction fingerprint, shaped by the [`ReactionCombinator`] that produced it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ReactionFingerprint {
    Difference(SignedFeatureSet<u64>),
    DisjointUnion(FeatureSet<(Side, u64)>),
}

/// Featurize `reaction` by applying `featurizer` to each side (`lhs` = reactants,
/// `rhs` = products) and combining per `combinator`. Both sides must be ground.
pub fn featurize_reaction(
    reaction: &ReactionAst,
    featurizer: &Featurizer,
    combinator: ReactionCombinator,
) -> Result<ReactionFingerprint, FingerprintError> {
    Ok(match combinator {
        ReactionCombinator::Difference => {
            let reactants = featurizer.featurize_counted(&reaction.lhs)?;
            let products = featurizer.featurize_counted(&reaction.rhs)?;
            ReactionFingerprint::Difference(SignedFeatureSet::difference(&products, &reactants))
        }
        ReactionCombinator::DisjointUnion => {
            let reactants = featurizer.featurize(&reaction.lhs)?;
            let products = featurizer.featurize(&reaction.rhs)?;
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
    use umol_ast::ast::ReactionAst;

    use super::*;
    use crate::fingerprint::{Featurizer, MorganFeaturizer};
    use crate::parse::parse_smiles;

    /// Morgan radius-0 oxygen invariant of ethanol — present in `CCO`, absent in `CC`.
    const ETHANOL_OXYGEN: u64 = 864662311;

    // Identity reaction: every feature count cancels, so the difference is empty.
    #[rstest]
    fn test_featurize_reaction_difference_identity() {
        let reaction = ReactionAst {
            lhs: parse_smiles("CCO").unwrap(),
            rhs: parse_smiles("CCO").unwrap(),
            atom_map: vec![],
            stereo_atoms: vec![],
            stereo_bonds: vec![],
        };
        let featurizer = Featurizer::Morgan(MorganFeaturizer::new(1));
        let fingerprint =
            featurize_reaction(&reaction, &featurizer, ReactionCombinator::Difference).unwrap();
        match fingerprint {
            ReactionFingerprint::Difference(difference) => assert!(difference.is_empty()),
            other => panic!("expected Difference, got {other:?}"),
        }
    }

    // CCO → CC removes the oxygen: its feature is products(0) − reactants(1) = −1.
    #[rstest]
    fn test_featurize_reaction_difference() {
        let reaction = ReactionAst {
            lhs: parse_smiles("CCO").unwrap(),
            rhs: parse_smiles("CC").unwrap(),
            atom_map: vec![],
            stereo_atoms: vec![],
            stereo_bonds: vec![],
        };
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
        let reaction = ReactionAst {
            lhs: parse_smiles("CCO").unwrap(),
            rhs: parse_smiles("CCO").unwrap(),
            atom_map: vec![],
            stereo_atoms: vec![],
            stereo_bonds: vec![],
        };
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
