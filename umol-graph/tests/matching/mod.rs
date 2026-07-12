#[path = "../../../umol-graph-core/tests/matching/fixture.rs"]
#[allow(dead_code)]
mod fixture;

use rstest::rstest;
use umol_ast::ast::{
    AromaticValenceAst, AtomAst, AtomConstraintAst, AtomId, BondAst, ElementAst, MoleculeAst,
    MoleculeParts, RingFamily, ValueAst,
};
use umol_chem::element::Element;
use umol_graph::ops::aromaticity::ClarAromaticity;

#[rstest]
fn test_clar_aromaticity_find_from_rings() {
    let fixture = fixture::parse(fixture::CORONENE);
    let atoms: Vec<_> = (0..fixture.node_count)
        .map(|_| {
            let mut atom = AtomAst::from_element(Element::C);
            atom.constraints.set(AtomConstraintAst::AromaticValence(
                AromaticValenceAst::Aromatic(ValueAst::Lit(1)),
            ));
            atom
        })
        .collect();
    let bonds = fixture
        .edges
        .iter()
        .map(|&[a, b]| (AtomId(a), AtomId(b), BondAst::from_order(1)))
        .collect();
    let ast = MoleculeAst::from_parts(MoleculeParts {
        atoms,
        bonds,
        ..Default::default()
    });
    let rings = ast.rings_with(RingFamily::Simple, 6, |_| true);
    let systems = ClarAromaticity
        .find_from_rings(&ast, &rings, &|view| match &view.ast.element {
            ElementAst::Lit(Element::C) => Some(1),
            _ => None,
        })
        .unwrap();

    assert_eq!(systems.len(), 1);
    assert_eq!(systems[0].0.len(), 18);
}
