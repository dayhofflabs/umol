use rstest::rstest;
use smallvec::smallvec;

use super::{derive_marker_components, MarkerComponent, MarkerComponentError, MarkerSite};
use crate::smiles::Smiles;
use crate::table_ir::{
    Atom, AtomPair, Bond, BondConfiguration, BondDonation, BondNoncovalent, BondOrder,
    BondRelation, Molecule, StereoBond,
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
