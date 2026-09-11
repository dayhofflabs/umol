use rstest::rstest;
use smallvec::{smallvec, SmallVec};

use super::{
    assign_markers, derive_marker_components, select_markers, MarkerAssignmentError,
    MarkerComponent, MarkerComponentError, MarkerSite, SelectionSite,
};
use crate::smiles::render::traversal::Traversal;
use crate::smiles::Smiles;
use crate::table_ir::BondDirection::{Falling, Rising};
use crate::table_ir::{
    Atom, AtomPair, Bond, BondConfiguration, BondDirection, BondDonation, BondNoncovalent,
    BondOrder, BondRelation, Molecule, StereoBond,
};

#[rustfmt::skip]
#[rstest]
#[case::empty("", vec![], vec![])]
#[case::unasserted("CC=CC=CC", vec![], vec![])]
#[case::either_only("CC=CC", vec![StereoBond { bond: 1, configuration: BondConfiguration::Either }], vec![])]
#[case::four_substituents("FC(Cl)=C(Br)I", vec![StereoBond {
    bond: 2, configuration: BondConfiguration::Framed { references: [2,5], relation: BondRelation::SameSide },
}], vec![MarkerComponent {
    candidates: vec![0,1,3,4], sites: vec![MarkerSite {
        bond: 2, candidates: [smallvec![0,1], smallvec![3,4]],
        configuration: Some(BondConfiguration::Framed { references: [2,5], relation: BondRelation::SameSide }),
    }],
}])]
#[case::shared_chain("CC=CC=CC", vec![StereoBond {
    bond: 1, configuration: BondConfiguration::Framed { references: [0,3], relation: BondRelation::OppositeSide },
}], vec![MarkerComponent {
    candidates: vec![0,2,4], sites: vec![
        MarkerSite { bond: 1, candidates: [smallvec![0], smallvec![2]], configuration: Some(BondConfiguration::Framed { references: [0,3], relation: BondRelation::OppositeSide }) },
        MarkerSite { bond: 3, candidates: [smallvec![2], smallvec![4]], configuration: None },
    ],
}])]
#[case::unasserted_coupling("CC=CC=CC=CC", vec![
    StereoBond { bond: 1, configuration: BondConfiguration::Framed { references: [0,3], relation: BondRelation::OppositeSide } },
    StereoBond { bond: 5, configuration: BondConfiguration::Framed { references: [4,7], relation: BondRelation::SameSide } },
], vec![MarkerComponent {
    candidates: vec![0,2,4,6], sites: vec![
        MarkerSite { bond: 1, candidates: [smallvec![0], smallvec![2]], configuration: Some(BondConfiguration::Framed { references: [0,3], relation: BondRelation::OppositeSide }) },
        MarkerSite { bond: 3, candidates: [smallvec![2], smallvec![4]], configuration: None },
        MarkerSite { bond: 5, candidates: [smallvec![4], smallvec![6]], configuration: Some(BondConfiguration::Framed { references: [4,7], relation: BondRelation::SameSide }) },
    ],
}])]
#[case::either_coupling("CC=CC=CC", vec![
    StereoBond { bond: 3, configuration: BondConfiguration::Either },
    StereoBond { bond: 1, configuration: BondConfiguration::Framed { references: [0,3], relation: BondRelation::OppositeSide } },
], vec![MarkerComponent {
    candidates: vec![0,2,4], sites: vec![
        MarkerSite { bond: 1, candidates: [smallvec![0], smallvec![2]], configuration: Some(BondConfiguration::Framed { references: [0,3], relation: BondRelation::OppositeSide }) },
        MarkerSite { bond: 3, candidates: [smallvec![2], smallvec![4]], configuration: Some(BondConfiguration::Either) },
    ],
}])]
#[case::localized_link("CC=CCCC=CC", vec![
    StereoBond { bond: 5, configuration: BondConfiguration::Framed { references: [4,7], relation: BondRelation::SameSide } },
    StereoBond { bond: 1, configuration: BondConfiguration::Framed { references: [0,3], relation: BondRelation::OppositeSide } },
], vec![
    MarkerComponent { candidates: vec![0,2], sites: vec![MarkerSite { bond: 1, candidates: [smallvec![0], smallvec![2]], configuration: Some(BondConfiguration::Framed { references: [0,3], relation: BondRelation::OppositeSide }) }] },
    MarkerComponent { candidates: vec![4,6], sites: vec![MarkerSite { bond: 5, candidates: [smallvec![4], smallvec![6]], configuration: Some(BondConfiguration::Framed { references: [4,7], relation: BondRelation::SameSide }) }] },
])]
#[case::discard_unasserted_component("CC=CC.CC=CC", vec![StereoBond {
    bond: 1, configuration: BondConfiguration::Framed { references: [0,3], relation: BondRelation::OppositeSide },
}], vec![MarkerComponent { candidates: vec![0,2], sites: vec![MarkerSite {
    bond: 1, candidates: [smallvec![0], smallvec![2]],
    configuration: Some(BondConfiguration::Framed { references: [0,3], relation: BondRelation::OppositeSide }),
}] }])]
#[case::branched("CC=C(C=CC)C=CC", vec![StereoBond {
    bond: 1, configuration: BondConfiguration::Framed { references: [0,3], relation: BondRelation::OppositeSide },
}], vec![MarkerComponent {
    candidates: vec![0,2,4,5,7], sites: vec![
        MarkerSite { bond: 1, candidates: [smallvec![0], smallvec![2,5]], configuration: Some(BondConfiguration::Framed { references: [0,3], relation: BondRelation::OppositeSide }) },
        MarkerSite { bond: 3, candidates: [smallvec![2], smallvec![4]], configuration: None },
        MarkerSite { bond: 6, candidates: [smallvec![5], smallvec![7]], configuration: None },
    ],
}])]
#[case::cycle("C1=CC=CC=C1", vec![StereoBond {
    bond: 1, configuration: BondConfiguration::Framed { references: [5,2], relation: BondRelation::SameSide },
}], vec![MarkerComponent {
    candidates: vec![0,2,4], sites: vec![
        MarkerSite { bond: 1, candidates: [smallvec![0], smallvec![2]], configuration: Some(BondConfiguration::Framed { references: [5,2], relation: BondRelation::SameSide }) },
        MarkerSite { bond: 3, candidates: [smallvec![2], smallvec![4]], configuration: None },
        MarkerSite { bond: 5, candidates: [smallvec![4], smallvec![0]], configuration: None },
    ],
}])]
#[case::uncoupled_endpoint("CC=CC=C#C", vec![StereoBond {
    bond: 1, configuration: BondConfiguration::Framed { references: [0,3], relation: BondRelation::SameSide },
}], vec![MarkerComponent { candidates: vec![0,2], sites: vec![MarkerSite {
    bond: 1, candidates: [smallvec![0], smallvec![2]],
    configuration: Some(BondConfiguration::Framed { references: [0,3], relation: BondRelation::SameSide }),
}] }])]
fn test_derive_marker_components(
    #[case] input: &str,
    #[case] frames: Vec<StereoBond>,
    #[case] expected: Vec<MarkerComponent>,
) {
    let mut molecule = Smiles::parse(input).unwrap().into_table_ir();
    molecule.stereo_bonds = frames;
    let original = molecule.clone();
    assert_eq!(derive_marker_components(&molecule), Ok(expected));
    assert_eq!(molecule, original);
}

#[rstest]
#[case::aromatic(Bond::new(1, 2, BondOrder::Aromatic))]
#[case::triple(Bond::new(1, 2, BondOrder::Triple))]
#[case::dative(Bond::new_dative(1, 2, BondOrder::Single, BondDonation::Donating))]
#[case::shared(Bond::new_dative(1, 2, BondOrder::Single, BondDonation::Shared))]
#[case::noncovalent(Bond::new_noncovalent(1, 2, BondNoncovalent::Hydrogen))]
fn test_derive_marker_components_candidates(#[case] excluded: Bond) {
    let configuration = BondConfiguration::Framed {
        references: [0, 4],
        relation: BondRelation::SameSide,
    };
    let molecule = Molecule {
        atoms: vec![Atom::wildcard(); 6],
        bonds: vec![
            Bond::new(0, 1, BondOrder::Single),
            excluded,
            Bond::new(1, 3, BondOrder::Double),
            Bond::new(3, 4, BondOrder::Single),
            Bond::new(3, 5, BondOrder::Single),
        ],
        stereo_bonds: vec![StereoBond {
            bond: 2,
            configuration,
        }],
        ..Molecule::empty()
    };
    assert_eq!(
        derive_marker_components(&molecule),
        Ok(vec![MarkerComponent {
            candidates: vec![0, 3, 4],
            sites: vec![MarkerSite {
                bond: 2,
                candidates: [smallvec![0], smallvec![3, 4]],
                configuration: Some(configuration),
            }],
        }])
    );
}

#[rstest]
#[case::terminal("C=CC", 0, MarkerComponentError::UnsupportedSite { bond: 0 })]
#[case::cumulene("CC=C=CC", 1, MarkerComponentError::UnsupportedSite { bond: 1 })]
#[case::overfull("CC(C)(C)=CC", 3, MarkerComponentError::UnsupportedSite { bond: 3 })]
#[case::shared_substituent("C1=CC1", 1, MarkerComponentError::UnsupportedSite { bond: 1 })]
#[case::missing_left("C#C=CC", 1, MarkerComponentError::MissingCandidate { bond: 1, atom: 1 })]
#[case::missing_right("CC=C#C", 1, MarkerComponentError::MissingCandidate { bond: 1, atom: 2 })]
#[case::single_site("CCCC", 1, MarkerComponentError::UnsupportedSite { bond: 1 })]
#[case::missing_site("CC=CC", u32::MAX, MarkerComponentError::BondIndexOutOfBounds { bond: u32::MAX })]
fn test_derive_marker_components_error(
    #[case] input: &str,
    #[case] bond: u32,
    #[case] expected: MarkerComponentError,
) {
    let mut molecule = Smiles::parse(input).unwrap().into_table_ir();
    molecule.stereo_bonds = vec![StereoBond {
        bond,
        configuration: BondConfiguration::Framed {
            references: [0, 3],
            relation: BondRelation::SameSide,
        },
    }];
    assert_eq!(derive_marker_components(&molecule), Err(expected));
}

#[rstest]
#[case::out_of_range(AtomPair::new(1, u32::MAX), MarkerComponentError::AtomIndexOutOfBounds { atom: u32::MAX })]
#[case::self_loop(AtomPair::new(1, 1), MarkerComponentError::UnsupportedSite { bond: 1 })]
fn test_derive_marker_components_endpoints(
    #[case] endpoints: AtomPair,
    #[case] expected: MarkerComponentError,
) {
    let mut molecule = Smiles::parse("C/C=C/C").unwrap().into_table_ir();
    molecule.bonds[1].atoms = endpoints;
    assert_eq!(derive_marker_components(&molecule), Err(expected));
}

#[rstest]
fn test_derive_marker_components_duplicates() {
    let mut molecule = Smiles::parse("C/C=C/C").unwrap().into_table_ir();
    molecule.stereo_bonds.push(molecule.stereo_bonds[0].clone());
    assert_eq!(
        derive_marker_components(&molecule),
        Err(MarkerComponentError::DuplicateStereoBond { bond: 1 })
    );
}

#[rstest]
#[case::empty("", vec![], Ok(vec![]))]
#[case::unasserted("FC=CF", vec![], Ok(vec![]))]
#[case::trans("FC=CF", vec![(1, [0,3], BondRelation::OppositeSide)], Ok(vec![(0,Rising),(2,Rising)]))]
#[case::cis("FC=CF", vec![(1, [0,3], BondRelation::SameSide)], Ok(vec![(0,Rising),(2,Falling)]))]
#[case::four_substituents("FC(Cl)=C(Br)I", vec![(2,[0,4],BondRelation::SameSide)], Ok(vec![(1,Rising),(4,Rising)]))]
#[case::reference_swap("FC(Cl)=C(Br)I", vec![(2,[2,4],BondRelation::SameSide)], Ok(vec![(1,Rising),(4,Falling)]))]
#[case::both_references_swapped("FC(Cl)=C(Br)I", vec![(2,[2,5],BondRelation::SameSide)], Ok(vec![(1,Rising),(4,Rising)]))]
#[case::noncandidate_reference("C(F)(#N)=C(Br)I", vec![(2,[2,4],BondRelation::SameSide)], Ok(vec![(0,Rising),(4,Rising)]))]
#[case::shared_chain("CC=CC=CC", vec![(1,[0,3],BondRelation::OppositeSide),(3,[2,5],BondRelation::OppositeSide)], Ok(vec![(0,Rising),(2,Rising),(4,Rising)]))]
#[case::both_side_candidates("CC=C(C=CC)C=CC", vec![(1,[0,3],BondRelation::OppositeSide),(3,[2,5],BondRelation::OppositeSide),(6,[2,8],BondRelation::OppositeSide)], Ok(vec![(0,Rising),(2,Rising),(4,Rising),(5,Falling),(7,Falling)]))]
#[case::partial_chain("CC=CC=CC", vec![(1,[0,3],BondRelation::OppositeSide)], Ok(vec![(0,Rising),(2,Rising)]))]
#[case::coverage_conflict("CC=CC=CC=CC", vec![(1,[0,3],BondRelation::OppositeSide),(5,[4,7],BondRelation::SameSide)], Err(MarkerAssignmentError::NoAssignment { bond: 1 }))]
#[case::branch_alternative("CC=C(C)C=CC=CC", vec![(1,[0,3],BondRelation::SameSide),(6,[5,8],BondRelation::SameSide)], Ok(vec![(0,Rising),(2,Falling),(5,Rising),(7,Falling)]))]
#[case::separate_groups("CC=CCCC=CC", vec![(1,[0,3],BondRelation::OppositeSide),(5,[4,7],BondRelation::SameSide)], Ok(vec![(0,Rising),(2,Rising),(4,Rising),(6,Falling)]))]
#[case::cycle("C1=CC=CC=C1", vec![(1,[5,2],BondRelation::SameSide),(3,[1,4],BondRelation::SameSide),(5,[3,0],BondRelation::OppositeSide)], Ok(vec![(0,Rising),(2,Rising),(4,Falling)]))]
#[case::cycle_conflict("C1=CC=CC=C1", vec![(1,[5,2],BondRelation::SameSide),(3,[1,4],BondRelation::SameSide),(5,[3,0],BondRelation::SameSide)], Err(MarkerAssignmentError::NoAssignment { bond: 1 }))]
#[case::cycle_alternative("C1=CC=C(C)C(C)=C1", vec![(1,[7,2],BondRelation::SameSide),(3,[1,5],BondRelation::SameSide),(7,[3,0],BondRelation::SameSide)], Ok(vec![(0,Rising),(2,Rising),(4,Rising),(6,Rising)]))]
#[case::invalid_reference("FC=CF", vec![(1,[0,8],BondRelation::SameSide)], Err(MarkerAssignmentError::InvalidReference { bond: 1, atom: 8 }))]
#[case::endpoint_reference("FC=CF", vec![(1,[1,3],BondRelation::SameSide)], Err(MarkerAssignmentError::InvalidReference { bond: 1, atom: 1 }))]
#[case::opposite_endpoint_reference("FC=CF", vec![(1,[2,3],BondRelation::SameSide)], Err(MarkerAssignmentError::InvalidReference { bond: 1, atom: 2 }))]
fn test_assign_markers(
    #[case] input: &str,
    #[case] frames: Vec<(u32, [u32; 2], BondRelation)>,
    #[case] expected: Result<Vec<(u32, BondDirection)>, MarkerAssignmentError>,
) {
    let mut molecule = Smiles::parse(input).unwrap().into_table_ir();
    molecule.stereo_bonds = frames
        .into_iter()
        .map(|(bond, references, relation)| StereoBond {
            bond,
            configuration: BondConfiguration::Framed {
                references,
                relation,
            },
        })
        .collect();
    assert_eq!(
        assign_markers(&molecule, &Traversal::new(&molecule)),
        expected
    );
}

#[rstest]
#[case::same(BondRelation::SameSide, vec![(0,Falling),(2,Rising)])]
#[case::opposite(BondRelation::OppositeSide, vec![(0,Falling),(2,Falling)])]
fn test_assign_markers_orientation(
    #[case] relation: BondRelation,
    #[case] expected: Vec<(u32, BondDirection)>,
) {
    let molecule = Molecule {
        atoms: vec![Atom::wildcard(); 4],
        bonds: vec![
            Bond::new(1, 2, BondOrder::Single),
            Bond::new(0, 2, BondOrder::Double),
            Bond::new(0, 3, BondOrder::Single),
        ],
        stereo_bonds: vec![StereoBond {
            bond: 1,
            configuration: BondConfiguration::Framed {
                references: [3, 1],
                relation,
            },
        }],
        ..Molecule::empty()
    };
    assert_eq!(
        assign_markers(&molecule, &Traversal::new(&molecule)),
        Ok(expected)
    );
}

#[rstest]
#[case::either()]
fn test_assign_markers_either() {
    let mut molecule = Smiles::parse("CC=CC=CC").unwrap().into_table_ir();
    molecule.stereo_bonds = vec![
        StereoBond {
            bond: 1,
            configuration: BondConfiguration::Framed {
                references: [0, 3],
                relation: BondRelation::SameSide,
            },
        },
        StereoBond {
            bond: 3,
            configuration: BondConfiguration::Either,
        },
    ];
    assert_eq!(
        assign_markers(&molecule, &Traversal::new(&molecule)),
        Ok(vec![(0, Rising), (2, Falling)])
    );
}

#[rstest]
#[case::parity_backtrack(
    vec![SelectionSite { candidates: [smallvec![0,1],smallvec![2]], definite: true }],
    vec![smallvec![(2,false)],smallvec![(2,false),(2,true)],smallvec![(0,false),(1,false),(1,true)]],
    vec![false,false,false], Some(vec![Some(false),None,Some(false)]),
)]
#[case::coverage_backtrack(
    vec![
        SelectionSite { candidates: [smallvec![0,1],smallvec![2]], definite: true },
        SelectionSite { candidates: [smallvec![1],smallvec![2]], definite: false },
    ],
    vec![smallvec![(2,false)],smallvec![],smallvec![(0,false)]],
    vec![false,false,false], Some(vec![Some(false),None,Some(false)]),
)]
#[case::split_groups(
    vec![
        SelectionSite { candidates: [smallvec![0],smallvec![1]], definite: true },
        SelectionSite { candidates: [smallvec![2],smallvec![3]], definite: true },
    ],
    vec![smallvec![(1,true)],smallvec![(0,true)],smallvec![(3,false)],smallvec![(2,false)]],
    vec![true,false,false,true], Some(vec![Some(true),Some(false),Some(false),Some(false)]),
)]
#[case::exhausted(
    vec![SelectionSite { candidates: [smallvec![0,1],smallvec![2]], definite: true }],
    vec![smallvec![(2,false),(2,true)],smallvec![(2,false),(2,true)],smallvec![(0,false),(0,true),(1,false),(1,true)]],
    vec![false,false,false], None,
)]
fn test_select_markers(
    #[case] sites: Vec<SelectionSite>,
    #[case] constraints: Vec<SmallVec<[(usize, bool); 4]>>,
    #[case] seeds: Vec<bool>,
    #[case] expected: Option<Vec<Option<bool>>>,
) {
    assert_eq!(select_markers(&sites, &constraints, &seeds), expected);
}
