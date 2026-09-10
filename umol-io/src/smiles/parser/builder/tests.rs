use rstest::rstest;

use super::*;

#[rstest]
#[case::unmarked(None, vec![(0, 1), (1, 2), (2, 2), (2, 1)], None)]
#[case::closing_order(Some(Chirality::Clockwise), vec![(0, 1), (1, 2), (2, 2), (2, 1)], Some(vec![1, 0]))]
#[case::boundaries(Some(Chirality::CounterClockwise), vec![(0, 1), (2, 1), (1, 2), (2, 2)], Some(vec![0, 1]))]
#[case::opening_only(Some(Chirality::Tetrahedral { arr: 1 }), vec![(2, 1), (0, 1)], Some(vec![0]))]
#[case::th2(Some(Chirality::Tetrahedral { arr: 2 }), vec![(0, 1), (2, 1)], Some(vec![0]))]
#[case::unsupported(Some(Chirality::SquarePlanar { arr: 1 }), vec![(0, 1), (2, 1)], None)]
fn test_molecule_editor_on_ring_bond(
    #[case] chirality: Option<Chirality>,
    #[case] digits: Vec<(usize, usize)>,
    #[case] expected: Option<Vec<usize>>,
) {
    let mut builder = MoleculeEditor::with_capacity(3, 2, false, None);
    let mut cursors = Vec::new();
    for atom in 0..3 {
        builder.current = None;
        builder.on_atom(AtomData {
            element: Some(Element::C),
            isotope: None,
            charge: None,
            implicit_hydrogens: None,
            class: None,
            aromatic: Some(false),
            chirality: if atom == 2 { chirality } else { None },
            span: None,
        });
        cursors.push(builder.current.unwrap());
    }
    for (pos, (atom, digit)) in digits.into_iter().enumerate() {
        builder.current = Some(cursors[atom]);
        builder.on_ring_bond(digit, pos, pos + 1, 0).unwrap();
    }
    assert_eq!(
        builder
            .stereo
            .iter()
            .map(|frame| frame.bonds.to_vec())
            .collect::<Vec<_>>(),
        expected.into_iter().collect::<Vec<_>>()
    );
    assert_eq!(builder.open_rings, 0);
    assert_eq!(builder.ring_bonds, vec![]);
}

#[rstest]
fn test_molecule_editor_finish() {
    let mut builder = MoleculeEditor::with_capacity(3, 0, false, None);
    for chirality in [
        Some(Chirality::CounterClockwise),
        None,
        Some(Chirality::Clockwise),
    ] {
        builder.current = None;
        builder.on_atom(AtomData {
            element: Some(Element::C),
            isotope: None,
            charge: None,
            implicit_hydrogens: None,
            class: None,
            aromatic: Some(false),
            chirality,
            span: None,
        });
    }
    let (molecule, rings) = builder.finish(0).unwrap();
    assert_eq!(
        molecule.stereo_atoms,
        vec![
            StereoAtom {
                atom: 0,
                winding: Winding::CounterClockwise,
                ligands: vec![]
            },
            StereoAtom {
                atom: 2,
                winding: Winding::Clockwise,
                ligands: vec![]
            },
        ]
    );
    assert_eq!(rings, vec![]);
}

#[rstest]
#[case::unmarked(None, vec![(0, 1), (1, 2), (2, 2), (2, 1)], None)]
#[case::closing_order(Some(Chirality::Clockwise), vec![(0, 1), (1, 2), (2, 2), (2, 1)], Some(vec![1, 0]))]
#[case::boundaries(Some(Chirality::CounterClockwise), vec![(0, 1), (2, 1), (1, 2), (2, 2)], Some(vec![0, 1]))]
#[case::opening_only(Some(Chirality::Tetrahedral { arr: 1 }), vec![(2, 1), (0, 1)], Some(vec![0]))]
#[case::th2(Some(Chirality::Tetrahedral { arr: 2 }), vec![(0, 1), (2, 1)], Some(vec![0]))]
#[case::unsupported(Some(Chirality::SquarePlanar { arr: 1 }), vec![(0, 1), (2, 1)], None)]
fn test_extended_molecule_builder_on_ring_bond(
    #[case] chirality: Option<Chirality>,
    #[case] digits: Vec<(usize, usize)>,
    #[case] expected: Option<Vec<usize>>,
) {
    let mut builder = ExtendedMoleculeBuilder::with_capacity(3, 2, false, None);
    let mut cursors = Vec::new();
    for atom in 0..3 {
        builder.current = None;
        builder.on_atom(ExtendedAtomData {
            symbol: AtomSymbol::Element(Element::C),
            isotope: None,
            charge: None,
            implicit_hydrogens: None,
            class: None,
            aromatic: false,
            chirality: if atom == 2 { chirality } else { None },
            span: None,
        });
        cursors.push(builder.current.unwrap());
    }
    for (pos, (atom, digit)) in digits.into_iter().enumerate() {
        builder.current = Some(cursors[atom]);
        builder.on_ring_bond(digit, pos, pos + 1, 0).unwrap();
    }
    assert_eq!(
        builder
            .stereo
            .iter()
            .map(|frame| frame.bonds.to_vec())
            .collect::<Vec<_>>(),
        expected.into_iter().collect::<Vec<_>>()
    );
    assert_eq!(builder.open_rings, 0);
    assert_eq!(builder.ring_bonds, vec![]);
}

#[rstest]
fn test_extended_molecule_builder_finish() {
    let mut builder = ExtendedMoleculeBuilder::with_capacity(3, 0, false, None);
    for chirality in [
        Some(Chirality::CounterClockwise),
        None,
        Some(Chirality::Clockwise),
    ] {
        builder.current = None;
        builder.on_atom(ExtendedAtomData {
            symbol: AtomSymbol::Element(Element::C),
            isotope: None,
            charge: None,
            implicit_hydrogens: None,
            class: None,
            aromatic: false,
            chirality,
            span: None,
        });
    }
    let (molecule, rings) = builder.finish(0).unwrap();
    assert_eq!(
        molecule.stereo_atoms,
        vec![
            StereoAtom {
                atom: 0,
                winding: Winding::CounterClockwise,
                ligands: vec![]
            },
            StereoAtom {
                atom: 2,
                winding: Winding::Clockwise,
                ligands: vec![]
            },
        ]
    );
    assert_eq!(rings, vec![]);
}
