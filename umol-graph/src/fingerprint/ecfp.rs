//! ECFP featurizer (Rogers & Hahn 2010): circular refinement over Daylight atom
//! invariants, with structural (bond-set) duplicate removal.
//!
//! Public entry points reject non-ground molecules before seed extraction. The
//! hash is [`EcfpHashScheme::Xxh3Width64V1`]; the paper leaves the hash
//! unspecified, so this is the frozen umol ECFP identity.

use umol_graph_core::CircularRefinementAlgorithm;
use umol_graph_ir::ir::{
    AsLit, AtomId, BondId, IsotopeMass, Molecule, RingConfig, RingModel, RingSet,
};

use super::feature_set::{CountedFeatureSet, FeatureSet};
use super::featurizer::FingerprintError;
use crate::hash::EcfpHashScheme;

/// ECFP fingerprint of `radius` iterations (diameter `2 * radius`, i.e. ECFP_{2r}).
#[derive(Clone, Copy, Debug)]
pub struct EcfpFeaturizer {
    pub radius: u32,
    pub hashing_scheme: EcfpHashScheme,
    pub ring_config: RingConfig,
}

impl EcfpFeaturizer {
    pub fn new(radius: u32) -> Self {
        Self {
            radius,
            hashing_scheme: EcfpHashScheme::default(),
            ring_config: RingConfig::default(),
        }
    }

    /// Returns the deduplicated set of feature identifiers.
    pub fn featurize(&self, mol: &Molecule) -> Result<FeatureSet<u64>, FingerprintError> {
        if !mol.is_ground() {
            return Err(FingerprintError::NotGround);
        }
        Ok(FeatureSet::from_features(self.identifiers(mol)))
    }

    /// Like [`Self::featurize`] but keeps per-identifier counts.
    pub fn featurize_counted(
        &self,
        mol: &Molecule,
    ) -> Result<CountedFeatureSet<u64>, FingerprintError> {
        if !mol.is_ground() {
            return Err(FingerprintError::NotGround);
        }
        Ok(CountedFeatureSet::from_features(self.identifiers(mol)))
    }

    /// The circular-refinement identifier multiset (one per surviving environment).
    fn identifiers(&self, mol: &Molecule) -> Vec<u64> {
        let rings = mol
            .rings(RingModel::default(), self.ring_config)
            .into_ring_set();
        mol.raw_graph().circular_refine(
            |node| atom_components(mol, &rings, AtomId::from(node)),
            |edge| bond_label(mol, BondId::from(edge)),
            CircularRefinementAlgorithm::Ec {
                radius: self.radius,
                scheme: self.hashing_scheme.recipe(),
            },
        )
    }
}

/// Daylight initial invariant (Rogers & Hahn 2010): heavy-atom degree, heavy
/// valence, atomic number, isotope mass (natural → 0), formal charge, attached
/// hydrogens, and ring membership.
fn atom_components(mol: &Molecule, rings: &RingSet, id: AtomId) -> Vec<u32> {
    let atom = mol.atom(id);
    let element = atom.element().as_lit().expect("ground atom");
    let isotope_mass = match atom.isotope_mass().as_lit().expect("ground atom") {
        IsotopeMass::Natural => 0,
        IsotopeMass::MassNumber(mass) => mass,
    };
    vec![
        u32::from(element.atomic_number()),
        atom.heavy_atom_degree().as_lit().expect("ground atom") as u32,
        atom.heavy_atom_valence().as_lit().expect("ground atom") as u32,
        atom.total_hydrogens().as_lit().expect("ground atom") as u32,
        atom.charge().as_lit().expect("ground atom") as i32 as u32,
        isotope_mass,
        u32::from(rings.contains_atom(id)),
    ]
}

/// Bond label: bond order with the aromatic-system flag.
fn bond_label(mol: &Molecule, id: BondId) -> u32 {
    let bond = mol.bond(id);
    let order = bond.order().as_lit().expect("ground bond") as u32;
    order | ((bond.is_in_aromatic_system() as u32) << 16)
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_ir::mol_dsl_ground;

    use super::*;

    const BUTYRAMIDE: &str = r#"{
        :atoms ["C #h3" "C #h2" "C #h2" "C #h0" "O #h0" "N #h2"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "2"] [3 5 "1"]]
    }"#;

    #[rstest]
    fn test_ecfp_featurizer_new() {
        let featurizer = EcfpFeaturizer::new(2);

        assert_eq!(featurizer.radius, 2);
        assert_eq!(featurizer.hashing_scheme, EcfpHashScheme::default());
        assert_eq!(featurizer.ring_config, RingConfig::default());
    }

    // Rogers & Hahn 2010, Figure 8 fixes the count per diameter; the exact ids
    // pin the frozen umol hash recipe for those same environments.
    #[rstest]
    #[case::diameter_0(
        0,
        &[
            1189585227353469813,
            1343896606611716210,
            6816650886737406922,
            9398025501618298006,
            16149328945726899460,
        ]
    )]
    #[case::diameter_2(
        1,
        &[
            686136971914186761,
            1189585227353469813,
            1343896606611716210,
            1674899844642375346,
            5686907935783274670,
            6158447595325937241,
            6816650886737406922,
            9398025501618298006,
            13652293261850732425,
            14550739996647717087,
            16149328945726899460,
        ]
    )]
    #[case::diameter_4(
        2,
        &[
            686136971914186761,
            1189585227353469813,
            1343896606611716210,
            1674899844642375346,
            5686907935783274670,
            6158447595325937241,
            6816650886737406922,
            9129806645566723864,
            9398025501618298006,
            13652293261850732425,
            14550739996647717087,
            16149328945726899460,
            16204012715323123438,
            16670450973526877804,
        ]
    )]
    #[case::diameter_6(
        3,
        &[
            686136971914186761,
            1189585227353469813,
            1343896606611716210,
            1674899844642375346,
            5686907935783274670,
            6158447595325937241,
            6816650886737406922,
            9129806645566723864,
            9398025501618298006,
            13652293261850732425,
            14550739996647717087,
            16149328945726899460,
            16204012715323123438,
            16670450973526877804,
        ]
    )]
    fn test_ecfp_featurizer_featurize_butyramide(#[case] radius: u32, #[case] expected: &[u64]) {
        let fingerprint = EcfpFeaturizer::new(radius)
            .featurize(&mol_dsl_ground!(BUTYRAMIDE))
            .unwrap();
        assert_eq!(fingerprint.ids(), expected);
    }

    #[rstest]
    #[case::relabeled_propane(
        r#"{:atoms ["C #h3" "C #h2" "C #h3"] :bonds [[0 1 "1"] [1 2 "1"]]}"#,
        r#"{:atoms ["C #h2" "C #h3" "C #h3"] :bonds [[0 1 "1"] [0 2 "1"]]}"#
    )]
    fn test_ecfp_featurizer_featurize_order_independent(#[case] a: &str, #[case] b: &str) {
        let featurizer = EcfpFeaturizer::new(2);
        assert_eq!(
            featurizer.featurize(&mol_dsl_ground!(a)).unwrap(),
            featurizer.featurize(&mol_dsl_ground!(b)).unwrap()
        );
    }
}
