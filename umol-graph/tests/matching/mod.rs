use rstest::rstest;
use serde::Deserialize;
use umol_chem::element::Element;
use umol_graph::ops::aromaticity::{AromaticityConfig, AromaticityPerceiver};
use umol_graph::ops::model::{AromaticityModel, AromaticityRule};
use umol_graph_core::MaximumIndependentSetAlgorithm;
use umol_graph_ir::ir::{
    AromaticSystemForm, AromaticValence, AromaticValenceForm, AtomConstraintForm, AtomForm, AtomId,
    BondForm, Molecule, MoleculeEntries, NumForm, UnpairedElectronsForm,
};
use umol_utils::solution::Solution;

#[derive(Deserialize)]
struct GraphFixture {
    node_count: usize,
    edges: Vec<[u32; 2]>,
}

#[rstest]
fn test_aromaticity_perceiver_find_systems_clar() {
    let fixture: GraphFixture = toml::from_str(include_str!("data/coronene_planar.toml")).unwrap();
    let atoms: Vec<_> = (0..fixture.node_count)
        .map(|_| {
            let mut atom = AtomForm::from_element(Element::C);
            atom.constraints.set(AtomConstraintForm::AromaticValence(
                AromaticValenceForm::Aromatic(NumForm::Lit(1)),
            ));
            atom
        })
        .collect();
    let bonds = fixture
        .edges
        .iter()
        .map(|&[a, b]| (AtomId(a), AtomId(b), BondForm::from_order(1)))
        .collect();
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms,
        bonds,
        ..Default::default()
    });
    let perception = AromaticityPerceiver::new(&AromaticityModel {
        rule: AromaticityRule::Clar,
        ..AromaticityModel::daylight()
    });
    assert_eq!(
        perception.find_systems(
            &molecule,
            AromaticityConfig {
                maximum_independent_set_algorithm: MaximumIndependentSetAlgorithm::BranchAndBound,
                ..Default::default()
            },
            |_| AromaticValence::Aromatic(1),
        ),
        Ok(Solution::Determined(vec![(
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15, 18, 19, 20, 21]
                .into_iter()
                .map(AtomId)
                .collect(),
            AromaticSystemForm::from_electrons(vec![1; 18])
                .with_charge(0)
                .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
        )]))
    );
}
