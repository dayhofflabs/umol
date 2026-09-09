use rstest::rstest;

use super::*;

#[rstest]
#[case::unmarked(None, vec![(0, 1), (1, 2), (2, 2), (2, 1)], vec![])]
#[case::closing_order(Some(Chirality::Clockwise), vec![(0, 1), (1, 2), (2, 2), (2, 1)], vec![(2, 1, 2), (2, 0, 2)])]
#[case::boundaries(Some(Chirality::CounterClockwise), vec![(0, 1), (2, 1), (1, 2), (2, 2)], vec![(2, 0, 1), (2, 1, 2)])]
#[case::opening_only(Some(Chirality::Tetrahedral { arr: 1 }), vec![(2, 1), (0, 1)], vec![])]
#[case::th2(Some(Chirality::Tetrahedral { arr: 2 }), vec![(0, 1), (2, 1)], vec![(2, 0, 1)])]
#[case::unsupported(Some(Chirality::SquarePlanar { arr: 1 }), vec![(0, 1), (2, 1)], vec![])]
fn test_molecule_editor_on_ring_bond(
    #[case] chirality: Option<Chirality>,
    #[case] digits: Vec<(usize, usize)>,
    #[case] expected: Vec<(u32, usize, usize)>,
) {
    let mut builder = MoleculeEditor::with_capacity(3, 2, false);
    builder.atoms = vec![Atom::from_element(Element::C); 3];
    builder.atoms[2].chirality = chirality;
    for (pos, (atom, digit)) in digits.into_iter().enumerate() {
        builder
            .on_ring_bond(atom, digit, None, None, None, pos, pos + 1, 0)
            .unwrap();
    }
    assert_eq!(builder.stereo_closures, expected);
    assert_eq!(builder.stereo_roots, vec![]);
    builder.on_component_end();
    assert_eq!(builder.stereo_closures, vec![]);
}

#[rstest]
fn test_molecule_editor_on_component_end() {
    let mut builder = MoleculeEditor::with_capacity(3, 0, false);
    builder.atoms = vec![Atom::from_element(Element::C); 3];
    builder.on_stereo_root(0);
    builder.on_stereo_root(2);
    assert_eq!(builder.stereo_roots, vec![0, 2]);
    builder.on_component_end();
    assert_eq!(builder.stereo_roots, vec![]);
}

#[rstest]
#[case::unmarked(None, vec![(0, 1), (1, 2), (2, 2), (2, 1)], vec![])]
#[case::closing_order(Some(Chirality::Clockwise), vec![(0, 1), (1, 2), (2, 2), (2, 1)], vec![(2, 1, 2), (2, 0, 2)])]
#[case::boundaries(Some(Chirality::CounterClockwise), vec![(0, 1), (2, 1), (1, 2), (2, 2)], vec![(2, 0, 1), (2, 1, 2)])]
#[case::opening_only(Some(Chirality::Tetrahedral { arr: 1 }), vec![(2, 1), (0, 1)], vec![])]
#[case::th2(Some(Chirality::Tetrahedral { arr: 2 }), vec![(0, 1), (2, 1)], vec![(2, 0, 1)])]
#[case::unsupported(Some(Chirality::SquarePlanar { arr: 1 }), vec![(0, 1), (2, 1)], vec![])]
fn test_extended_molecule_builder_on_ring_bond(
    #[case] chirality: Option<Chirality>,
    #[case] digits: Vec<(usize, usize)>,
    #[case] expected: Vec<(u32, usize, usize)>,
) {
    let mut builder = ExtendedMoleculeBuilder::with_capacity(3, 2, false);
    builder.atoms = vec![ExtendedAtom::from_element(Element::C); 3];
    builder.atoms[2].chirality = chirality;
    for (pos, (atom, digit)) in digits.into_iter().enumerate() {
        builder
            .on_ring_bond(atom, digit, None, None, None, pos, pos + 1, 0)
            .unwrap();
    }
    assert_eq!(builder.stereo_closures, expected);
    assert_eq!(builder.stereo_roots, vec![]);
    builder.on_component_end();
    assert_eq!(builder.stereo_closures, vec![]);
}

#[rstest]
fn test_extended_molecule_builder_on_component_end() {
    let mut builder = ExtendedMoleculeBuilder::with_capacity(3, 0, false);
    builder.atoms = vec![ExtendedAtom::from_element(Element::C); 3];
    builder.on_stereo_root(0);
    builder.on_stereo_root(2);
    assert_eq!(builder.stereo_roots, vec![0, 2]);
    builder.on_component_end();
    assert_eq!(builder.stereo_roots, vec![]);
}
