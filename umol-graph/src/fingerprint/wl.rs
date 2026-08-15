//! Weisfeiler–Lehman featurizer: frozen color refinement over the atom graph.

use umol_graph_core::{EdgeId, Refinement, RefinementAlgorithm, RefinementRounds};
use umol_graph_ir::ir::{AsLit, AtomId, BondId, Molecule};

use super::feature_set::{CountedFeatureSet, FeatureSet};
use super::featurizer::FingerprintError;
use crate::hash::WlHashScheme;

/// Weisfeiler–Lehman color-refinement fingerprint over the atom graph, hashed
/// through a frozen `hashing_scheme` for `rounds` rounds.
#[derive(Clone, Copy, Debug)]
pub struct WlFeaturizer {
    pub rounds: RefinementRounds,
    pub hashing_scheme: WlHashScheme,
}

impl WlFeaturizer {
    pub fn new(rounds: RefinementRounds) -> Self {
        Self {
            rounds,
            hashing_scheme: WlHashScheme::default(),
        }
    }

    /// Returns the deduplicated set of feature identifiers.
    pub fn featurize(&self, mol: &Molecule) -> Result<FeatureSet<u64>, FingerprintError> {
        if !mol.is_concrete() {
            return Err(FingerprintError::NotConcrete);
        }
        Ok(FeatureSet::from_sorted_unique(
            self.refinement(mol).features(),
        ))
    }

    /// Like [`Self::featurize`] but keeps per-identifier counts.
    pub fn featurize_counted(
        &self,
        mol: &Molecule,
    ) -> Result<CountedFeatureSet<u64>, FingerprintError> {
        if !mol.is_concrete() {
            return Err(FingerprintError::NotConcrete);
        }
        Ok(CountedFeatureSet::from_counts(
            self.refinement(mol).counts(),
        ))
    }

    fn refinement(&self, mol: &Molecule) -> Refinement<u64> {
        // Pre-extract seeds so the refine closures stay simple index lookups.
        let atom_seeds: Vec<u64> = (0..mol.atoms().count())
            .map(|i| atom_seed(mol, AtomId(i as u32)))
            .collect();
        let bond_seeds: Vec<u64> = (0..mol.bonds().count())
            .map(|e| bond_seed(mol, BondId(e as u32)))
            .collect();

        mol.raw_graph().refine(
            |node| atom_seeds[AtomId::from(node).index()],
            |edge: EdgeId| bond_seeds[edge.index()],
            RefinementAlgorithm::WeisfeilerLehman {
                rounds: self.rounds,
                scheme: self.hashing_scheme.refinement_scheme(),
            },
        )
    }
}

/// Seed an atom from (atomic number, formal charge, implicit hydrogens), each in
/// its own byte range so distinct tuples stay distinct before the scheme rehashes.
fn atom_seed(mol: &Molecule, id: AtomId) -> u64 {
    let atom = mol.atom(id);
    let atomic_number = atom
        .element()
        .as_lit()
        .expect("ground atom")
        .atomic_number();
    let charge = atom.charge().as_lit().expect("ground atom");
    let implicit_hydrogens = atom.implicit_hydrogens().as_lit().expect("ground atom");
    (atomic_number as u64)
        | (((charge as u16) as u64) << 8)
        | (((implicit_hydrogens as u16) as u64) << 24)
}

/// Seed a bond from (bond order, aromatic-system membership).
fn bond_seed(mol: &Molecule, id: BondId) -> u64 {
    let bond = mol.bond(id);
    let order = bond.order().as_lit().expect("ground bond");
    (order as u16 as u64) | ((bond.is_in_aromatic_system() as u64) << 16)
}

#[cfg(test)]
mod tests {
    use rstest::{fixture, rstest};
    use umol_graph_core::RefinementRounds;
    use umol_graph_ir::mol_dsl_concrete;

    use super::*;

    #[fixture]
    fn featurizer() -> WlFeaturizer {
        WlFeaturizer::new(RefinementRounds::Fixed(3))
    }

    #[rstest]
    #[case::ethane(
        r#"{:atoms ["C #h3" "C #h3"] :bonds [[0 1 "1"]]}"#,
        vec![
            2659163409134283895,
            7542810387455301591,
            9541344068636876323,
            12512207080905326651
        ]
    )]
    fn test_wl_featurizer_featurize(
        featurizer: WlFeaturizer,
        #[case] edn: &str,
        #[case] expected: Vec<u64>,
    ) {
        assert_eq!(
            featurizer.featurize(&mol_dsl_concrete!(edn)).unwrap().ids(),
            expected.as_slice()
        );
    }

    #[rstest]
    #[case::relabeled_propane(
        r#"{:atoms ["C #h3" "C #h2" "C #h3"] :bonds [[0 1 "1"] [1 2 "1"]]}"#,
        r#"{:atoms ["C #h2" "C #h3" "C #h3"] :bonds [[0 1 "1"] [0 2 "1"]]}"#,
        true
    )]
    #[case::ethane_vs_propane(
        r#"{:atoms ["C #h3" "C #h3"] :bonds [[0 1 "1"]]}"#,
        r#"{:atoms ["C #h3" "C #h2" "C #h3"] :bonds [[0 1 "1"] [1 2 "1"]]}"#,
        false
    )]
    #[case::propane_vs_isopropyl_cation(
        r#"{:atoms ["C #h3" "C #h2" "C #h3"] :bonds [[0 1 "1"] [1 2 "1"]]}"#,
        r#"{:atoms ["C #h3" "C #h1 #c+" "C #h3"] :bonds [[0 1 "1"] [1 2 "1"]]}"#,
        false
    )]
    fn test_wl_featurizer_featurize_equivalent(
        featurizer: WlFeaturizer,
        #[case] a: &str,
        #[case] b: &str,
        #[case] equal: bool,
    ) {
        let same = featurizer.featurize(&mol_dsl_concrete!(a)).unwrap()
            == featurizer.featurize(&mol_dsl_concrete!(b)).unwrap();
        assert_eq!(same, equal);
    }
}
