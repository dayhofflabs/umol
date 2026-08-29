//! Representation-integrity preservation by chemistry-layer publishers.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph::ingest::ingest_smiles;
use umol_graph::ops::model::{AromaticityModel, ChemistryModel, ValenceModel};
use umol_graph::ops::resolve::Resolver;
use umol_graph::ops::transform::{
    Aromatizer, DelocalizeCharge, KekulizeConfig, Kekulizer, MaximumMatchingAlgorithm, Transformer,
};
use umol_graph_ir::ir::{AtomId, Molecule};

fn molecule_from_smiles(source: &str) -> Molecule {
    ingest_smiles(source).expect("publication-property SMILES fixture resolves")
}

fn kekule_molecule_strategy() -> impl Strategy<Value = Molecule> {
    prop::sample::select(vec!["C1=CC=CC=C1", "C1=CC=CN=C1"]).prop_map(molecule_from_smiles)
}

fn aromatic_molecule_strategy() -> impl Strategy<Value = Molecule> {
    prop::sample::select(vec!["c1ccccc1", "c1ccncc1"]).prop_map(molecule_from_smiles)
}

fn resolved_molecule_strategy() -> impl Strategy<Value = Molecule> {
    prop::sample::select(vec!["C", "CC", "C=C", "N", "O", "c1ccccc1"])
        .prop_map(molecule_from_smiles)
}

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_aromatizer_integrity_preservation(source in kekule_molecule_strategy()) {
        let published = Aromatizer::new(&AromaticityModel::daylight())
            .transform(&source)
            .map_err(|error| TestCaseError::fail(format!("aromatization failed: {error}")))?;

        prop_assert_eq!(published.edit().try_build(), Ok(published));
    }

    #[test]
    fn test_delocalize_charge_integrity_preservation(source in aromatic_molecule_strategy()) {
        let published = DelocalizeCharge
            .transform(&source)
            .expect("delocalized-charge transformation is infallible");

        prop_assert_eq!(published.edit().try_build(), Ok(published));
    }

    #[test]
    fn test_kekulizer_integrity_preservation(
        source in aromatic_molecule_strategy(),
        use_bipartite in any::<bool>(),
    ) {
        let algorithm = if use_bipartite {
            MaximumMatchingAlgorithm::HopcroftKarp
        } else {
            MaximumMatchingAlgorithm::Edmonds
        };
        let node_order = source.atoms().ids().collect::<Vec<AtomId>>();
        let published = Kekulizer::new(KekulizeConfig::new(algorithm), node_order)
            .transform(&source)
            .map_err(|error| TestCaseError::fail(format!("kekulization failed: {error}")))?;

        prop_assert_eq!(published.edit().try_build(), Ok(published));
    }

    #[test]
    fn test_resolver_integrity_preservation(mut source in resolved_molecule_strategy()) {
        let model = ChemistryModel {
            valence: ValenceModel::smiles(),
            ..ChemistryModel::default()
        };
        Resolver::new(&model)
            .resolve(&mut source)
            .map_err(|error| TestCaseError::fail(format!("resolution failed: {error}")))?;

        prop_assert_eq!(source.edit().try_build(), Ok(source));
    }
}
