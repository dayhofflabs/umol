use rstest::rstest;
use umol_graph_core::AutomorphismAlgorithm as GraphCoreAutomorphismAlgorithm;
use umol_graph_ir::mol_dsl_concrete;

use super::*;

#[rstest]
fn test_graph_symmetry_copy() {
    let molecule = mol_dsl_concrete!(
        r#"{:atoms ["C" "H" "F" "Cl" "Br"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#
    );
    let symmetry = GraphSymmetry::from_rust(
        molecule
            .graph_symmetry(&GraphIrGraphSymmetryConfig::default())
            .into_graph_symmetry(),
    );
    let copied = symmetry.copy();
    assert!(molecule.is_chiral(symmetry.to_rust()));
    assert!(molecule.is_chiral(copied.to_rust()));
    assert_eq!(
        molecule.proper_orbit_of(copied.to_rust(), GraphIrAtomId(0)),
        vec![GraphIrAtomId(0)]
    );
}

#[rstest]
fn test_graph_symmetry_config_default() {
    assert_eq!(
        GraphSymmetryConfig::default().to_rust(),
        GraphIrGraphSymmetryConfig {
            automorphism_algorithm: GraphCoreAutomorphismAlgorithm::Nauty,
            max_classification_trials: 4096,
            include_constraints: false,
        }
    );
}

#[rstest]
#[case::without_constraints(GraphIrGraphSymmetryConfig {
    automorphism_algorithm: GraphCoreAutomorphismAlgorithm::Nauty,
    max_classification_trials: 4096,
    include_constraints: false,
})]
#[case::with_constraints(GraphIrGraphSymmetryConfig {
    automorphism_algorithm: GraphCoreAutomorphismAlgorithm::Nauty,
    max_classification_trials: 16,
    include_constraints: true,
})]
fn test_graph_symmetry_config_roundtrip(#[case] config: GraphIrGraphSymmetryConfig) {
    assert_eq!(GraphSymmetryConfig::from_rust(config).to_rust(), config);
}
