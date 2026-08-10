use rstest::rstest;
use serde::Deserialize;
use umol_chem::element::Element;
use umol_graph::ops::aromaticity::ClarAromaticity;
use umol_graph_core::MaximumIndependentSetAlgorithm;
use umol_graph_ir::ir::{
    AromaticValenceForm, AtomConstraintForm, AtomForm, AtomId, BondForm, ElementForm, Molecule,
    MoleculeEntries, NumForm, RingConfig, RingModel, RingSetKind,
};

#[derive(Deserialize)]
struct GraphFixture {
    node_count: usize,
    edges: Vec<[u32; 2]>,
}

#[rstest]
fn test_clar_aromaticity_find_from_rings() {
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
    let rings = molecule
        .rings(
            RingModel {
                kind: RingSetKind::Relevant,
                max_ring_size: 6,
            },
            RingConfig::default(),
        )
        .into_ring_set();
    let systems = ClarAromaticity
        .find_from_rings(
            &molecule,
            &rings,
            MaximumIndependentSetAlgorithm::BranchAndBound,
            &|view| match &view.attributes.element {
                ElementForm::Lit(Element::C) => Some(1),
                _ => None,
            },
        )
        .unwrap();

    assert_eq!(systems.len(), 1);
    assert_eq!(systems[0].0.len(), 18);
}
