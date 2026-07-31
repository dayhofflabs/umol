use std::collections::BTreeMap;

use indexmap::IndexMap;
use pretty_assertions::assert_eq;
use rstest::rstest;
use umol_chem::element::Element;
use umol_chem::spin::SpinMultiplicity;

use super::*;
use crate::table_ir::{
    Atom, AtomStereoCare, AtomSymbol, Bond, BondOrder, Chirality, ConversionError, CtfileData,
    CxAnnotationData, ExtendedAtom, ExtendedBond, MulticenterBond, MulticenterSet, RGroup,
    RGroupOccurrence, SGroup, SGroupType, SourceFormat, Span, WildcardAtom,
};

#[rstest]
fn test_molecule_empty() {
    assert_eq!(
        Molecule::empty(),
        Molecule {
            atoms: vec![],
            bonds: vec![],
            positions: None,
            multicenter_bonds: vec![],
            configuration_scope: None,
            chirality_frame: None,
            comments: vec![],
            properties: IndexMap::new(),
            source_format: SourceFormat::UNKNOWN,
        }
    );
}

#[rstest]
#[case::three(vec![
    Atom::aliphatic_atom(Element::C),
    Atom::aromatic_atom(Element::N),
    Atom::aliphatic_atom(Element::O),
], 3)]
fn test_molecule_atom_count(#[case] atoms: Vec<Atom>, #[case] expected: usize) {
    assert_eq!(
        Molecule {
            atoms,
            ..Molecule::empty()
        }
        .atom_count(),
        expected
    );
}

#[rstest]
#[case::two(vec![
    Bond::new(0, 1, BondOrder::Single),
    Bond::new(1, 2, BondOrder::Double),
], 2)]
fn test_molecule_bond_count(#[case] bonds: Vec<Bond>, #[case] expected: usize) {
    assert_eq!(
        Molecule {
            bonds,
            ..Molecule::empty()
        }
        .bond_count(),
        expected
    );
}

#[rstest]
#[case::one(
    vec![MulticenterBond::new(vec![MulticenterSet::new(vec![0, 1, 2])])],
    1
)]
fn test_molecule_multicenter_bond_count(
    #[case] multicenter_bonds: Vec<MulticenterBond>,
    #[case] expected: usize,
) {
    assert_eq!(
        Molecule {
            multicenter_bonds,
            ..Molecule::empty()
        }
        .multicenter_bond_count(),
        expected
    );
}

#[rstest]
#[case::two(
    IndexMap::from([
        ("name".to_owned(), "methane".to_owned()),
        ("source".to_owned(), "test".to_owned()),
    ]),
    2
)]
fn test_molecule_property_count(
    #[case] properties: IndexMap<String, String>,
    #[case] expected: usize,
) {
    assert_eq!(
        Molecule {
            properties,
            ..Molecule::empty()
        }
        .property_count(),
        expected
    );
}

#[rstest]
#[case::carbon(vec![Atom::from_element(Element::C)], "C")]
#[case::methane(vec![
    Atom::from_element(Element::C),
    Atom::from_element(Element::H),
    Atom::from_element(Element::H),
    Atom::from_element(Element::H),
    Atom::from_element(Element::H),
], "CH4")]
#[case::charged(vec![Atom {
    charge: Some(1),
    ..Atom::from_element(Element::C)
}], "C+")]
#[case::hill_order(vec![
    Atom::from_element(Element::C),
    Atom::from_element(Element::C),
    Atom::from_element(Element::O),
    Atom::from_element(Element::N),
], "C2NO")]
#[case::single_wildcard(vec![
    Atom::wildcard(),
    Atom::from_element(Element::C),
    Atom::from_element(Element::O),
    Atom::from_element(Element::O),
    Atom::from_element(Element::C),
    Atom::from_element(Element::C),
], "C3O2*")]
#[case::multiple_wildcards(vec![
    Atom::wildcard(),
    Atom::from_element(Element::C),
    Atom::wildcard(),
], "C*2")]
fn test_molecule_sum_formula(#[case] atoms: Vec<Atom>, #[case] expected: &str) {
    assert_eq!(
        Molecule {
            atoms,
            ..Molecule::empty()
        }
        .sum_formula(),
        expected
    );
}

#[rstest]
fn test_extended_molecule_empty() {
    assert_eq!(
        ExtendedMolecule::empty(),
        ExtendedMolecule {
            atoms: vec![],
            bonds: vec![],
            positions: None,
            multicenter_bonds: vec![],
            configuration_scope: None,
            chirality_frame: None,
            comments: vec![],
            properties: IndexMap::new(),
            ctfile_data: None,
            cx_data: None,
            source_format: SourceFormat::UNKNOWN,
        }
    );
}

#[rstest]
#[case::two(vec![
    ExtendedAtom::from_element(Element::C),
    ExtendedAtom::from_element(Element::N),
], 2)]
fn test_extended_molecule_atom_count(#[case] atoms: Vec<ExtendedAtom>, #[case] expected: usize) {
    assert_eq!(
        ExtendedMolecule {
            atoms,
            ..ExtendedMolecule::empty()
        }
        .atom_count(),
        expected
    );
}

#[rstest]
#[case::two(vec![
    ExtendedBond::new(0, 1, BondOrder::Single),
    ExtendedBond::new(1, 2, BondOrder::Double),
], 2)]
fn test_extended_molecule_bond_count(#[case] bonds: Vec<ExtendedBond>, #[case] expected: usize) {
    assert_eq!(
        ExtendedMolecule {
            bonds,
            ..ExtendedMolecule::empty()
        }
        .bond_count(),
        expected
    );
}

#[rstest]
#[case::two(
    vec![
        MulticenterBond::new(vec![MulticenterSet::new(vec![0, 1, 2])]),
        MulticenterBond::new(vec![MulticenterSet::new(vec![3, 4, 5])]),
    ],
    2
)]
fn test_extended_molecule_multicenter_bond_count(
    #[case] multicenter_bonds: Vec<MulticenterBond>,
    #[case] expected: usize,
) {
    assert_eq!(
        ExtendedMolecule {
            multicenter_bonds,
            ..ExtendedMolecule::empty()
        }
        .multicenter_bond_count(),
        expected
    );
}

#[rstest]
fn test_extended_molecule_sgroups() {
    let ctfile_group = SGroup::new(SGroupType::Superatom);
    let cx_group = SGroup::new(SGroupType::Data);
    let molecule = ExtendedMolecule {
        ctfile_data: Some(CtfileData {
            sgroups: BTreeMap::from([(1, ctfile_group.clone())]),
            ..Default::default()
        }),
        cx_data: Some(CxAnnotationData {
            sgroups: BTreeMap::from([(1, cx_group.clone())]),
            ..Default::default()
        }),
        ..ExtendedMolecule::empty()
    };

    let view = molecule.sgroups();
    assert!(!view.is_empty());
    assert_eq!(view.len(), 2);
    assert_eq!(view.get(&1), Some(&ctfile_group));
    assert!(view.contains_key(&1));
    assert!(!view.contains_key(&2));
    assert_eq!(
        view.iter()
            .map(|(key, group)| (key, group.clone()))
            .collect::<Vec<_>>(),
        vec![(1, ctfile_group), (1, cx_group)]
    );
}

#[rstest]
fn test_extended_molecule_sgroups_mut() {
    let group = SGroup::new(SGroupType::Superatom);
    let mut molecule = ExtendedMolecule::empty();
    molecule.sgroups_mut().insert(1, group.clone());
    assert_eq!(
        molecule.cx_data,
        Some(CxAnnotationData {
            sgroups: BTreeMap::from([(1, group)]),
            ..Default::default()
        })
    );
}

#[rstest]
fn test_extended_molecule_rgroups() {
    let ctfile_group = RGroup::new(Some(1));
    let cx_group = RGroup::new(Some(2));
    let molecule = ExtendedMolecule {
        ctfile_data: Some(CtfileData {
            rgroups: BTreeMap::from([(1, ctfile_group.clone())]),
            ..Default::default()
        }),
        cx_data: Some(CxAnnotationData {
            rgroups: BTreeMap::from([(1, cx_group.clone())]),
            ..Default::default()
        }),
        ..ExtendedMolecule::empty()
    };

    let view = molecule.rgroups();
    assert!(!view.is_empty());
    assert_eq!(view.len(), 2);
    assert_eq!(view.get(&1), Some(&ctfile_group));
    assert!(view.contains_key(&1));
    assert!(!view.contains_key(&2));
    assert_eq!(
        view.iter()
            .map(|(key, group)| (key, group.clone()))
            .collect::<Vec<_>>(),
        vec![(1, ctfile_group), (1, cx_group)]
    );
}

#[rstest]
fn test_extended_molecule_rgroups_mut() {
    let group = RGroup {
        label: Some(1),
        dependent_label: None,
        rgroup_or_h: false,
        occurrence: vec![RGroupOccurrence::Range(1, 3)],
    };
    let mut molecule = ExtendedMolecule::empty();
    molecule.rgroups_mut().insert(1, group.clone());
    assert_eq!(
        molecule.cx_data,
        Some(CxAnnotationData {
            rgroups: BTreeMap::from([(1, group)]),
            ..Default::default()
        })
    );
}

#[rstest]
#[case::two(
    IndexMap::from([
        ("name".to_owned(), "methane".to_owned()),
        ("source".to_owned(), "test".to_owned()),
    ]),
    2
)]
fn test_extended_molecule_property_count(
    #[case] properties: IndexMap<String, String>,
    #[case] expected: usize,
) {
    assert_eq!(
        ExtendedMolecule {
            properties,
            ..ExtendedMolecule::empty()
        }
        .property_count(),
        expected
    );
}

#[rstest]
#[case::methane(vec![
    ExtendedAtom::from_element(Element::C),
    ExtendedAtom::from_element(Element::H),
    ExtendedAtom::from_element(Element::H),
    ExtendedAtom::from_element(Element::H),
    ExtendedAtom::from_element(Element::H),
], "CH4")]
#[case::single_wildcard(vec![
    ExtendedAtom::from_atom_symbol(AtomSymbol::WildcardAtom(WildcardAtom::Any)),
    ExtendedAtom::from_element(Element::C),
    ExtendedAtom::from_element(Element::O),
    ExtendedAtom::from_element(Element::O),
    ExtendedAtom::from_element(Element::C),
    ExtendedAtom::from_element(Element::C),
], "C3O2*")]
#[case::multiple_wildcards(vec![
    ExtendedAtom::from_atom_symbol(AtomSymbol::WildcardAtom(WildcardAtom::Any)),
    ExtendedAtom::from_element(Element::C),
    ExtendedAtom::from_atom_symbol(AtomSymbol::WildcardAtom(WildcardAtom::Any)),
], "C*2")]
fn test_extended_molecule_sum_formula(#[case] atoms: Vec<ExtendedAtom>, #[case] expected: &str) {
    assert_eq!(
        ExtendedMolecule {
            atoms,
            ..ExtendedMolecule::empty()
        }
        .sum_formula(),
        expected
    );
}

#[rstest]
fn test_extended_molecule_from() {
    let atoms = vec![
        Atom::aliphatic_atom(Element::C),
        Atom::aromatic_atom(Element::N),
    ];
    let bonds = vec![Bond::new(0, 1, BondOrder::Double)];
    let molecule = Molecule {
        atoms: atoms.clone(),
        bonds: bonds.clone(),
        source_format: SourceFormat::MOL,
        ..Molecule::empty()
    };
    assert_eq!(
        ExtendedMolecule::from(molecule),
        ExtendedMolecule {
            atoms: atoms.into_iter().map(ExtendedAtom::from).collect(),
            bonds: bonds.into_iter().map(ExtendedBond::from).collect(),
            source_format: SourceFormat::MOL,
            ..ExtendedMolecule::empty()
        }
    );
}

#[rstest]
fn test_molecule_try_from() {
    let extended = ExtendedMolecule {
        atoms: vec![
            ExtendedAtom::from_element(Element::C),
            ExtendedAtom::from_element(Element::O),
        ],
        bonds: vec![ExtendedBond::new(0, 1, BondOrder::Single)],
        source_format: SourceFormat::SMILES,
        ..ExtendedMolecule::empty()
    };
    assert_eq!(
        Molecule::try_from(extended),
        Ok(Molecule {
            atoms: vec![
                Atom::from_element(Element::C),
                Atom::from_element(Element::O),
            ],
            bonds: vec![Bond::new(0, 1, BondOrder::Single)],
            source_format: SourceFormat::SMILES,
            ..Molecule::empty()
        })
    );
}

#[rstest]
fn test_molecule_try_from_error() {
    let extended = ExtendedMolecule {
        atoms: vec![ExtendedAtom {
            stereo_care: Some(AtomStereoCare::Care),
            ..ExtendedAtom::from_element(Element::C)
        }],
        ..ExtendedMolecule::empty()
    };
    assert_eq!(
        Molecule::try_from(extended),
        Err(ConversionError::HasExtendedFeatures)
    );
}

#[rstest]
fn test_molecule_try_from_roundtrip() {
    let molecule = Molecule {
        atoms: vec![
            Atom::from_element(Element::C),
            Atom {
                element: None,
                isotope_mass: Some(13),
                charge: Some(-1),
                implicit_hydrogens: Some(2),
                valence: Some(4),
                lone_pairs: Some(1),
                unpaired_electrons: Some(2),
                multiplicity: Some(SpinMultiplicity::SINGLET),
                aromatic: None,
                chirality: Some(Chirality::Unspecified),
                class: Some(7),
                label: Some("wildcard".to_owned()),
                value: Some("value".to_owned()),
                span: Some(Span::bytes(1, 8)),
            },
        ],
        bonds: vec![Bond::new(0, 1, BondOrder::Single)],
        source_format: SourceFormat::SMILES,
        ..Molecule::empty()
    };

    let extended = ExtendedMolecule::from(molecule.clone());
    assert_eq!(Molecule::try_from(extended), Ok(molecule));
}
