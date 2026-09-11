use rstest::rstest;
use umol_perm::{ClassKey, Permutation};

use super::{AtomVisit, RingVisit, Traversal};
use crate::smiles::Smiles;
use crate::table_ir::{Atom, Bond, BondOrder, Molecule, Neighbor, StereoLigand};

#[rustfmt::skip]
#[rstest]
#[case::empty(0, vec![], vec![], vec![])]
#[case::isolated(3, vec![], vec![
    AtomVisit { atom: 0, parent: None, subtree_end: 1, rings: 0..0 },
    AtomVisit { atom: 1, parent: None, subtree_end: 2, rings: 0..0 },
    AtomVisit { atom: 2, parent: None, subtree_end: 3, rings: 0..0 },
], vec![])]
#[case::branches(6, vec![(0,1), (1,2), (0,3), (3,4), (0,5)], vec![
    AtomVisit { atom: 0, parent: None, subtree_end: 6, rings: 0..0 },
    AtomVisit { atom: 1, parent: Some(Neighbor { atom: 0, bond: 0 }), subtree_end: 3, rings: 0..0 },
    AtomVisit { atom: 2, parent: Some(Neighbor { atom: 1, bond: 1 }), subtree_end: 3, rings: 0..0 },
    AtomVisit { atom: 3, parent: Some(Neighbor { atom: 0, bond: 2 }), subtree_end: 5, rings: 0..0 },
    AtomVisit { atom: 4, parent: Some(Neighbor { atom: 3, bond: 3 }), subtree_end: 5, rings: 0..0 },
    AtomVisit { atom: 5, parent: Some(Neighbor { atom: 0, bond: 4 }), subtree_end: 6, rings: 0..0 },
], vec![])]
#[case::table_order(4, vec![(2,3), (0,2), (0,1)], vec![
    AtomVisit { atom: 0, parent: None, subtree_end: 4, rings: 0..0 },
    AtomVisit { atom: 2, parent: Some(Neighbor { atom: 0, bond: 1 }), subtree_end: 3, rings: 0..0 },
    AtomVisit { atom: 3, parent: Some(Neighbor { atom: 2, bond: 0 }), subtree_end: 3, rings: 0..0 },
    AtomVisit { atom: 1, parent: Some(Neighbor { atom: 0, bond: 2 }), subtree_end: 4, rings: 0..0 },
], vec![])]
#[case::ring(3, vec![(0,1), (1,2), (0,2)], vec![
    AtomVisit { atom: 0, parent: None, subtree_end: 3, rings: 0..1 },
    AtomVisit { atom: 1, parent: Some(Neighbor { atom: 0, bond: 0 }), subtree_end: 3, rings: 1..1 },
    AtomVisit { atom: 2, parent: Some(Neighbor { atom: 1, bond: 1 }), subtree_end: 3, rings: 1..2 },
], vec![
    RingVisit { neighbor: Neighbor { atom: 2, bond: 2 }, label: 1, opening: true },
    RingVisit { neighbor: Neighbor { atom: 0, bond: 2 }, label: 1, opening: false },
])]
#[case::multiple_closures(4, vec![(0,1), (1,2), (2,3), (0,2), (0,3), (1,3)], vec![
    AtomVisit { atom: 0, parent: None, subtree_end: 4, rings: 0..2 },
    AtomVisit { atom: 1, parent: Some(Neighbor { atom: 0, bond: 0 }), subtree_end: 4, rings: 2..3 },
    AtomVisit { atom: 2, parent: Some(Neighbor { atom: 1, bond: 1 }), subtree_end: 4, rings: 3..4 },
    AtomVisit { atom: 3, parent: Some(Neighbor { atom: 2, bond: 2 }), subtree_end: 4, rings: 4..6 },
], vec![
    RingVisit { neighbor: Neighbor { atom: 2, bond: 3 }, label: 1, opening: true },
    RingVisit { neighbor: Neighbor { atom: 3, bond: 4 }, label: 2, opening: true },
    RingVisit { neighbor: Neighbor { atom: 3, bond: 5 }, label: 3, opening: true },
    RingVisit { neighbor: Neighbor { atom: 0, bond: 3 }, label: 1, opening: false },
    RingVisit { neighbor: Neighbor { atom: 0, bond: 4 }, label: 2, opening: false },
    RingVisit { neighbor: Neighbor { atom: 1, bond: 5 }, label: 3, opening: false },
])]
#[case::close_and_open(5, vec![(0,1), (1,2), (0,2), (2,3), (3,4), (2,4)], vec![
    AtomVisit { atom: 0, parent: None, subtree_end: 5, rings: 0..1 },
    AtomVisit { atom: 1, parent: Some(Neighbor { atom: 0, bond: 0 }), subtree_end: 5, rings: 1..1 },
    AtomVisit { atom: 2, parent: Some(Neighbor { atom: 1, bond: 1 }), subtree_end: 5, rings: 1..3 },
    AtomVisit { atom: 3, parent: Some(Neighbor { atom: 2, bond: 3 }), subtree_end: 5, rings: 3..3 },
    AtomVisit { atom: 4, parent: Some(Neighbor { atom: 3, bond: 4 }), subtree_end: 5, rings: 3..4 },
], vec![
    RingVisit { neighbor: Neighbor { atom: 2, bond: 2 }, label: 1, opening: true },
    RingVisit { neighbor: Neighbor { atom: 0, bond: 2 }, label: 1, opening: false },
    RingVisit { neighbor: Neighbor { atom: 4, bond: 5 }, label: 2, opening: true },
    RingVisit { neighbor: Neighbor { atom: 2, bond: 5 }, label: 2, opening: false },
])]
#[case::reuse_across_components(6, vec![(0,1), (1,2), (0,2), (3,4), (4,5), (3,5)], vec![
    AtomVisit { atom: 0, parent: None, subtree_end: 3, rings: 0..1 },
    AtomVisit { atom: 1, parent: Some(Neighbor { atom: 0, bond: 0 }), subtree_end: 3, rings: 1..1 },
    AtomVisit { atom: 2, parent: Some(Neighbor { atom: 1, bond: 1 }), subtree_end: 3, rings: 1..2 },
    AtomVisit { atom: 3, parent: None, subtree_end: 6, rings: 2..3 },
    AtomVisit { atom: 4, parent: Some(Neighbor { atom: 3, bond: 3 }), subtree_end: 6, rings: 3..3 },
    AtomVisit { atom: 5, parent: Some(Neighbor { atom: 4, bond: 4 }), subtree_end: 6, rings: 3..4 },
], vec![
    RingVisit { neighbor: Neighbor { atom: 2, bond: 2 }, label: 1, opening: true },
    RingVisit { neighbor: Neighbor { atom: 0, bond: 2 }, label: 1, opening: false },
    RingVisit { neighbor: Neighbor { atom: 5, bond: 5 }, label: 1, opening: true },
    RingVisit { neighbor: Neighbor { atom: 3, bond: 5 }, label: 1, opening: false },
])]
#[case::loop_and_parallel(2, vec![(0,0), (0,1), (0,1)], vec![
    AtomVisit { atom: 0, parent: None, subtree_end: 2, rings: 0..3 },
    AtomVisit { atom: 1, parent: Some(Neighbor { atom: 0, bond: 1 }), subtree_end: 2, rings: 3..4 },
], vec![
    RingVisit { neighbor: Neighbor { atom: 0, bond: 0 }, label: 1, opening: true },
    RingVisit { neighbor: Neighbor { atom: 0, bond: 0 }, label: 1, opening: false },
    RingVisit { neighbor: Neighbor { atom: 1, bond: 2 }, label: 2, opening: true },
    RingVisit { neighbor: Neighbor { atom: 0, bond: 2 }, label: 2, opening: false },
])]
#[case::invalid_endpoints(2, vec![(0,5), (0,1), (u32::MAX,1)], vec![
    AtomVisit { atom: 0, parent: None, subtree_end: 2, rings: 0..0 },
    AtomVisit { atom: 1, parent: Some(Neighbor { atom: 0, bond: 1 }), subtree_end: 2, rings: 0..0 },
], vec![])]
fn test_traversal_new(
    #[case] count: usize,
    #[case] bonds: Vec<(u32, u32)>,
    #[case] atoms: Vec<AtomVisit>,
    #[case] rings: Vec<RingVisit>,
) {
    let molecule = Molecule {
        atoms: vec![Atom::wildcard(); count],
        bonds: bonds.into_iter().map(|(a, b)| Bond::new(a, b, BondOrder::Single)).collect(),
        ..Molecule::empty()
    };
    let expected = Traversal { atoms, rings };
    assert_eq!(Traversal::new(&molecule), expected);
}

#[rstest]
#[case::branches("CC(C)C(C)C", 1, vec![2, 3])]
#[case::nested("CC(C)C(C)C", 3, vec![4, 5])]
#[case::leaf("CC(C)C(C)C", 2, vec![])]
#[case::root("C.C(C)(C)C", 1, vec![2, 3, 4])]
#[case::absent("", 0, vec![])]
fn test_traversal_children(
    #[case] input: &str,
    #[case] position: usize,
    #[case] expected: Vec<u32>,
) {
    let molecule = Smiles::parse(input).unwrap().into_table_ir();
    let traversal = Traversal::new(&molecule);
    assert_eq!(
        traversal
            .children(position)
            .map(|atom| atom.atom)
            .collect::<Vec<_>>(),
        expected
    );
}

#[rstest]
#[case::root("[C@](F)(Cl)(Br)I", 0, vec![1, 2, 3, 4], 0)]
#[case::branch_before_closure("O1CCC[C@](F)1Cl", 4, vec![0, 3, 5, 6], 0)]
#[case::ring_root("[C@]12(CCC1)CCC2", 0, vec![1, 4, 3, 6], 0)]
#[case::mixed_closures("O1CCC[C@]21CCCC2", 4, vec![0, 5, 3, 8], 0)]
#[case::branch_return("[C@](CCC1)(F)(Cl)1", 0, vec![3, 1, 4, 5], 1)]
#[case::explicit_hydrogen("C[C@]1([H])CCCCO1", 1, vec![0, 3, 7, 2], 0)]
fn test_traversal_neighbors(
    #[case] input: &str,
    #[case] site: u32,
    #[case] expected: Vec<u32>,
    #[case] expected_coset: u32,
) {
    let molecule = Smiles::parse(input).unwrap().into_table_ir();
    let original = molecule.clone();
    let traversal = Traversal::new(&molecule);
    let position = traversal
        .atoms
        .iter()
        .position(|atom| atom.atom == site)
        .unwrap();
    let actual: Vec<_> = traversal
        .neighbors(position)
        .map(|neighbor| neighbor.atom)
        .collect();
    assert_eq!(actual, expected);
    let frame = &molecule.stereo_atoms[0].ligands;
    let emitted: Vec<_> = actual.into_iter().map(StereoLigand::Atom).collect();
    let action = Permutation::between(frame, &emitted).unwrap();
    assert_eq!(
        ClassKey::Tetrahedral.space().reindex(0, action),
        Some(expected_coset)
    );
    assert_eq!(molecule, original);
}

#[rstest]
#[case::chain(3, vec![(2,0), (1,2)], vec![(0, [0,2]), (1, [2,1])])]
#[case::ring(3, vec![(0,1), (1,2), (0,2)], vec![(2, [0,2]), (0, [0,1]), (1, [1,2])])]
#[case::fused(4, vec![(0,1), (1,2), (2,3), (0,2), (0,3), (1,3)],
    vec![(3, [0,2]), (4, [0,3]), (0, [0,1]), (5, [1,3]), (1, [1,2]), (2, [2,3])])]
fn test_traversal_bonds(
    #[case] count: usize,
    #[case] bonds: Vec<(u32, u32)>,
    #[case] expected: Vec<(u32, [u32; 2])>,
) {
    let molecule = Molecule {
        atoms: vec![Atom::wildcard(); count],
        bonds: bonds
            .into_iter()
            .map(|(a, b)| Bond::new(a, b, BondOrder::Single))
            .collect(),
        ..Molecule::empty()
    };
    assert_eq!(
        Traversal::new(&molecule).bonds().collect::<Vec<_>>(),
        expected
    );
}
