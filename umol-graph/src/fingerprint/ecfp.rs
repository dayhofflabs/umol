//! ECFP featurizer (Rogers & Hahn 2010): circular refinement over Daylight atom
//! invariants, with structural (bond-set) duplicate removal.
//!
//! The molecule must be ground; the caller ([`super::Featurizer::featurize`])
//! guarantees it, so the seed reads concrete literals directly. The hash is
//! `xxh3_64` with [`ECFP_SEED`]; the paper leaves the hash unspecified, so this is
//! a frozen choice (placeholder identity, like the WL schemes).

use umol_ast::ast::{AsLit, AtomId, BondId, MoleculeAst};
use umol_graph_core::CircularRefinementAlgorithm;

use crate::hash::{RogersHahn, ECFP_SEED};
use super::feature_set::{CountedFeatureSet, FeatureSet};

/// ECFP fingerprint of `radius` iterations (diameter `2 * radius`, i.e. ECFP_{2r}).
#[derive(Clone, Copy, Debug)]
pub struct EcfpFeaturizer {
    pub radius: u32,
    pub seed: u64,
}

impl EcfpFeaturizer {
    pub fn new(radius: u32) -> Self {
        Self {
            radius,
            seed: ECFP_SEED,
        }
    }

    /// `mol` must be ground. Returns the deduplicated set of feature identifiers.
    pub fn featurize(&self, mol: &MoleculeAst) -> FeatureSet<u64> {
        FeatureSet::from_features(self.identifiers(mol))
    }

    /// `mol` must be ground. Like [`Self::featurize`] but keeps per-identifier counts.
    pub fn featurize_counted(&self, mol: &MoleculeAst) -> CountedFeatureSet<u64> {
        CountedFeatureSet::from_features(self.identifiers(mol))
    }

    /// The circular-refinement identifier multiset (one per surviving environment).
    fn identifiers(&self, mol: &MoleculeAst) -> Vec<u64> {
        mol.raw_graph().circular_refine(
            |node| atom_components(mol, AtomId::from(node)),
            |edge| bond_label(mol, BondId::from(edge)),
            CircularRefinementAlgorithm::Ec {
                radius: self.radius,
                scheme: RogersHahn { seed: self.seed },
            },
        )
    }
}

/// Daylight initial invariant (Rogers & Hahn 2010): heavy-atom degree, heavy
/// valence, atomic number, isotope mass (natural → 0), formal charge, attached
/// hydrogens, and ring membership.
fn atom_components(mol: &MoleculeAst, id: AtomId) -> Vec<u32> {
    let atom = mol.atom(id);
    let element = atom.element().as_lit().expect("ground atom");
    vec![
        u32::from(element.atomic_number()),
        atom.heavy_atom_degree().as_lit().expect("ground atom") as u32,
        atom.heavy_atom_valence().as_lit().expect("ground atom") as u32,
        atom.total_hydrogens().as_lit().expect("ground atom") as u32,
        atom.charge().as_lit().expect("ground atom") as i32 as u32,
        atom.isotope_mass().as_lit().unwrap_or(0),
        u32::from(atom.is_in_ring()),
    ]
}

/// Bond label: bond order with the aromatic-system flag.
fn bond_label(mol: &MoleculeAst, id: BondId) -> u32 {
    let bond = mol.bond(id);
    let order = bond.order().as_lit().expect("ground bond") as u32;
    order | ((bond.is_in_aromatic_system() as u32) << 16)
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::mol_ground;

    use super::*;

    const BUTYRAMIDE: &str = r#"{
        :atoms ["C #h3" "C #h2" "C #h2" "C #h0" "O #h0" "N #h2"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "2"] [3 5 "1"]]
    }"#;

    // Rogers & Hahn 2010, Figure 8: butyramide feature counts per diameter.
    #[rstest]
    #[case::diameter_0(0, 5)]
    #[case::diameter_2(1, 11)]
    #[case::diameter_4(2, 14)]
    #[case::diameter_6(3, 14)]
    fn test_ecfp_featurizer_featurize_butyramide(#[case] radius: u32, #[case] expected: usize) {
        let fingerprint = EcfpFeaturizer::new(radius).featurize(&mol_ground!(BUTYRAMIDE));
        assert_eq!(fingerprint.len(), expected);
    }

    #[rstest]
    #[case::relabeled_propane(
        r#"{:atoms ["C #h3" "C #h2" "C #h3"] :bonds [[0 1 "1"] [1 2 "1"]]}"#,
        r#"{:atoms ["C #h2" "C #h3" "C #h3"] :bonds [[0 1 "1"] [0 2 "1"]]}"#
    )]
    fn test_ecfp_featurizer_featurize_order_independent(#[case] a: &str, #[case] b: &str) {
        let featurizer = EcfpFeaturizer::new(2);
        assert_eq!(
            featurizer.featurize(&mol_ground!(a)),
            featurizer.featurize(&mol_ground!(b))
        );
    }
}
