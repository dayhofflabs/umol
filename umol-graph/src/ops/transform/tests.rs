pub(super) mod stereo;

use pretty_assertions::assert_eq;
use rstest::{fixture, rstest};
use umol_graph_ir::ir::{AromaticSystemId, AtomId, Molecule};
use umol_graph_ir::mol_dsl_concrete;

use crate::ops::aromaticity::{AromaticityConfig, AromaticityContradiction};
use crate::ops::model::{AromaticityModel, AromaticityRule};
use crate::ops::transform::{
    AromatizeError, FoldHydrogensError, KekulizeConfig, KekulizeError, MaximumMatchingAlgorithm,
    Transform, UnfoldHydrogensError,
};

#[fixture]
fn localized_benzene() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [0 5 "1"]]
    }"#
    )
}

#[fixture]
fn aromatic_benzene() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
    }"#
    )
}

#[rstest]
fn test_molecule_aromatize(localized_benzene: Molecule, aromatic_benzene: Molecule) {
    let model = AromaticityModel::daylight();
    let config = AromaticityConfig::default();
    assert_eq!(
        localized_benzene.clone().aromatize(&model, config),
        Ok(aromatic_benzene.clone())
    );
    let mut molecule = localized_benzene;
    assert_eq!(molecule.aromatize_into(&model, config), Ok(()));
    assert_eq!(molecule, aromatic_benzene);
}

#[rstest]
#[case::clar_heterocycle(
    mol_dsl_concrete!(r#"{
        :atoms ["N#n" "C#h" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [0 5 "1"]]
    }"#),
    AromaticityModel { rule: AromaticityRule::Clar, ..AromaticityModel::daylight() },
    AromatizeError::Contradiction(AromaticityContradiction::ClarNonBenzenoid(
        "Clar model requires benzenoid input but non-carbon aromatic atoms are present".into()
    )),
)]
fn test_molecule_aromatize_error(
    #[case] input: Molecule,
    #[case] model: AromaticityModel,
    #[case] expected: AromatizeError,
) {
    assert_eq!(
        input
            .clone()
            .aromatize(&model, AromaticityConfig::default()),
        Err(expected.clone())
    );
    let mut molecule = input.clone();
    assert_eq!(
        molecule.aromatize_into(&model, AromaticityConfig::default()),
        Err(expected)
    );
    assert_eq!(molecule, input);
}

#[rstest]
fn test_molecule_kekulize(aromatic_benzene: Molecule, localized_benzene: Molecule) {
    assert_eq!(
        aromatic_benzene
            .clone()
            .kekulize((0..6).map(AtomId).collect(), KekulizeConfig::default()),
        Ok(localized_benzene.clone())
    );
    let mut molecule = aromatic_benzene;
    assert_eq!(
        molecule.kekulize_into((0..6).map(AtomId).collect(), KekulizeConfig::default()),
        Ok(())
    );
    assert_eq!(molecule, localized_benzene);
}

#[rstest]
#[case::missing_atom(
    aromatic_benzene(),
    vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4)],
    KekulizeConfig::default(),
    KekulizeError::InvalidNodeOrder {
        aromatic_system: AromaticSystemId(0),
        missing: vec![AtomId(5)],
        duplicates: vec![],
    },
)]
#[case::non_bipartite(
    mol_dsl_concrete!(r#"{
        :atoms ["C#a" "C#a" "C#a"] :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"]]
        :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,1]"}]
    }"#),
    vec![AtomId(0), AtomId(1), AtomId(2)],
    KekulizeConfig::new(MaximumMatchingAlgorithm::HopcroftKarp),
    KekulizeError::NonBipartiteMatching(AromaticSystemId(0)),
)]
fn test_molecule_kekulize_error(
    #[case] input: Molecule,
    #[case] node_order: Vec<AtomId>,
    #[case] config: KekulizeConfig,
    #[case] expected: KekulizeError,
) {
    assert_eq!(
        input.clone().kekulize(node_order.clone(), config),
        Err(expected.clone())
    );
    let mut molecule = input.clone();
    assert_eq!(molecule.kekulize_into(node_order, config), Err(expected));
    assert_eq!(molecule, input);
}

#[rstest]
fn test_molecule_kekulize_into_sequence(localized_benzene: Molecule, aromatic_benzene: Molecule) {
    let mut molecule = localized_benzene;
    assert_eq!(
        molecule.aromatize_into(&AromaticityModel::daylight(), AromaticityConfig::default()),
        Ok(())
    );
    assert_eq!(molecule, aromatic_benzene);
    assert_eq!(
        molecule.kekulize_into(vec![], KekulizeConfig::default()),
        Err(KekulizeError::InvalidNodeOrder {
            aromatic_system: AromaticSystemId(0),
            missing: (0..6).map(AtomId).collect(),
            duplicates: vec![],
        })
    );
    assert_eq!(molecule, aromatic_benzene);
}

#[rstest]
#[case::cyclopropenyl(
    mol_dsl_concrete!(r#"{
        :atoms ["C#c+#h#a0" "C#h#a" "C#h#a"]
        :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"]]
        :aromatic-systems [{:atoms [0 1 2] :attrs "[0,1,1]"}]
    }"#),
    mol_dsl_concrete!(r#"{
        :atoms ["C#h#a" "C#h#a" "C#h#a"]
        :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"]]
        :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,1]#c+"}]
    }"#),
)]
fn test_molecule_delocalize_charge(#[case] input: Molecule, #[case] expected: Molecule) {
    assert_eq!(input.clone().delocalize_charge(), Ok(expected.clone()));
    let mut molecule = input;
    assert_eq!(molecule.delocalize_charge_into(), Ok(()));
    assert_eq!(molecule, expected);
}

#[rstest]
#[case::methyl(
    mol_dsl_concrete!(r#"{:atoms ["H" "C#h2"] :bonds [[0 1 "1"]]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C#h3"]}"#),
)]
fn test_molecule_fold_hydrogens(#[case] input: Molecule, #[case] expected: Molecule) {
    assert_eq!(input.clone().fold_hydrogens(), Ok(expected.clone()));
    let mut molecule = input;
    assert_eq!(molecule.fold_hydrogens_into(), Ok(()));
    assert_eq!(molecule, expected);
}

#[rstest]
#[case::negative_count(
    mol_dsl_concrete!(r#"{:atoms ["C#h-1" "H"] :bonds [[0 1 "1"]]}"#),
    FoldHydrogensError::NegativeImplicitHydrogenCount { atom: AtomId(0), count: -1 },
)]
fn test_molecule_fold_hydrogens_error(
    #[case] input: Molecule,
    #[case] expected: FoldHydrogensError,
) {
    assert_eq!(input.clone().fold_hydrogens(), Err(expected.clone()));
    let mut molecule = input.clone();
    assert_eq!(molecule.fold_hydrogens_into(), Err(expected));
    assert_eq!(molecule, input);
}

#[rstest]
#[case::methyl(
    mol_dsl_concrete!(r#"{:atoms ["C#h3"]}"#),
    mol_dsl_concrete!(r#"{
        :atoms ["C" "H" "H" "H"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
    }"#),
)]
fn test_molecule_unfold_hydrogens(#[case] input: Molecule, #[case] expected: Molecule) {
    assert_eq!(input.clone().unfold_hydrogens(), Ok(expected.clone()));
    let mut molecule = input;
    assert_eq!(molecule.unfold_hydrogens_into(), Ok(()));
    assert_eq!(molecule, expected);
}

#[rstest]
#[case::negative_after_applicable(
    mol_dsl_concrete!(r#"{:atoms ["C#h3" "O#h-1"]}"#),
    UnfoldHydrogensError::NegativeImplicitHydrogenCount { atom: AtomId(1), count: -1 },
)]
fn test_molecule_unfold_hydrogens_error(
    #[case] input: Molecule,
    #[case] expected: UnfoldHydrogensError,
) {
    assert_eq!(input.clone().unfold_hydrogens(), Err(expected.clone()));
    let mut molecule = input.clone();
    assert_eq!(molecule.unfold_hydrogens_into(), Err(expected));
    assert_eq!(molecule, input);
}
