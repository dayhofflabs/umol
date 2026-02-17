use indexmap::IndexMap;
use pretty_assertions::assert_eq;
use umol_data::Element;

use super::*;
use crate::bond::BondDonation;
use crate::position::Point3D;
use crate::table_ir::{
    Atom, AtomStereoCare, BicycloStereo, BicycloStereoData, Bond, BondOrder, Chirality,
    ConversionError, CtfileData, CxAnnotationData, ExtendedAtom, ExtendedBond, JoinError,
    LegacyGroupAbbreviation, LocalParityCenter, MulticenterBond, MulticenterSet, RGroup,
    RGroupOccurrence, SGroup, SGroupType, SourceFormat, StereoSet, StereoSetMode,
};

#[test]
fn test_molecule_empty() {
    let mol = Molecule::empty();
    assert!(mol.atoms.is_empty());
    assert!(mol.bonds.is_empty());
    assert!(mol.rings.is_empty());
    assert_eq!(mol.source_format, SourceFormat::UNKNOWN);
}

#[test]
fn test_molecule_with_atoms_and_bonds() {
    let mol = Molecule {
        atoms: vec![
            Atom::aliphatic_atom(Element::C),
            Atom::aliphatic_atom(Element::O),
        ],
        bonds: vec![Bond::new(0, 1, BondOrder::Single)],
        rings: vec![],
        positions: None,
        multicenter_bonds: vec![],
        stereo_interpretation: None,
        comments: vec![],
        properties: IndexMap::new(),
        source_format: SourceFormat::SMILES,
    };
    assert_eq!(mol.atom_count(), 2);
    assert_eq!(mol.bond_count(), 1);
    assert_eq!(mol.bonds[0].start_atom(), 0);
    assert_eq!(mol.bonds[0].end_atom(), 1);
}

#[test]
fn test_molecule_atom_count() {
    let mol = Molecule {
        atoms: vec![
            Atom::aliphatic_atom(Element::C),
            Atom::aromatic_atom(Element::N),
            Atom::aliphatic_atom(Element::O),
        ],
        bonds: vec![],
        rings: vec![],
        positions: None,
        multicenter_bonds: vec![],
        stereo_interpretation: None,
        comments: vec![],
        properties: IndexMap::new(),
        source_format: SourceFormat::SMILES,
    };
    assert_eq!(mol.atom_count(), 3);
}

#[test]
fn test_molecule_bond_count() {
    let mol = Molecule {
        atoms: vec![
            Atom::aliphatic_atom(Element::C),
            Atom::aliphatic_atom(Element::C),
            Atom::aliphatic_atom(Element::C),
        ],
        bonds: vec![
            Bond::new(0, 1, BondOrder::Single),
            Bond::new(1, 2, BondOrder::Double),
        ],
        rings: vec![],
        positions: None,
        multicenter_bonds: vec![],
        stereo_interpretation: None,
        comments: vec![],
        properties: IndexMap::new(),
        source_format: SourceFormat::SMILES,
    };
    assert_eq!(mol.bond_count(), 2);
}

#[test]
fn test_molecule_multicenter_bond_count() {
    let mol = Molecule {
        atoms: vec![
            Atom::aliphatic_atom(Element::C),
            Atom::aliphatic_atom(Element::C),
            Atom::aliphatic_atom(Element::C),
        ],
        bonds: vec![],
        rings: vec![],
        positions: None,
        multicenter_bonds: vec![MulticenterBond::new(vec![MulticenterSet::new(
            vec![0, 1, 2],
            3,
        )])],
        stereo_interpretation: None,
        comments: vec![],
        properties: IndexMap::new(),
        source_format: SourceFormat::SMILES,
    };
    assert_eq!(mol.multicenter_bond_count(), 1);
}

#[test]
fn test_molecule_component_count() {
    let mol = Molecule {
        atoms: vec![
            Atom::aliphatic_atom(Element::C),
            Atom::aliphatic_atom(Element::C),
            Atom::aliphatic_atom(Element::C),
        ],
        bonds: vec![Bond::new(0, 1, BondOrder::Single)],
        rings: vec![],
        positions: None,
        multicenter_bonds: vec![],
        stereo_interpretation: None,
        comments: vec![],
        properties: IndexMap::new(),
        source_format: SourceFormat::SMILES,
    };
    assert_eq!(mol.component_count(), 2);
    assert_eq!(mol.component_atom_indices(), vec![vec![0, 1], vec![2]]);
}

#[test]
fn test_molecule_component_count_multicenter_bonds() {
    let mol = Molecule {
        atoms: vec![
            Atom::aliphatic_atom(Element::C),
            Atom::aliphatic_atom(Element::C),
            Atom::aliphatic_atom(Element::C),
        ],
        bonds: vec![],
        rings: vec![],
        positions: None,
        multicenter_bonds: vec![MulticenterBond::new(vec![MulticenterSet::new(
            vec![0, 1, 2],
            3,
        )])],
        stereo_interpretation: None,
        comments: vec![],
        properties: IndexMap::new(),
        source_format: SourceFormat::SMILES,
    };
    assert_eq!(mol.component_count(), 1);
    assert_eq!(mol.component_atom_indices(), vec![vec![0, 1, 2]]);
}

#[test]
fn test_molecule_split_components_remaps_indices() {
    let mut dative = Bond::new_dative(0, 1, BondOrder::Single, BondDonation::Donating);
    dative.ring = Some(7);
    let mol = Molecule {
        atoms: vec![
            Atom::aliphatic_atom(Element::C),
            Atom::aliphatic_atom(Element::N),
            Atom::aliphatic_atom(Element::O),
            Atom::aliphatic_atom(Element::H),
        ],
        bonds: vec![dative, Bond::new(2, 3, BondOrder::Double)],
        rings: vec![Ring {
            ring_idx: 7,
            start_atom: Some(0),
            end_atom: Some(1),
            open_span: None,
            close_span: None,
        }],
        positions: Some(vec![
            Point3D::new(10.0, 0.0, 0.0),
            Point3D::new(11.0, 0.0, 0.0),
            Point3D::new(20.0, 0.0, 0.0),
            Point3D::new(21.0, 0.0, 0.0),
        ]),
        multicenter_bonds: vec![MulticenterBond::new(vec![MulticenterSet::new(
            vec![2, 3],
            2,
        )])],
        stereo_interpretation: None,
        comments: vec!["molecule comment".to_string()],
        properties: IndexMap::from_iter([("k".to_string(), "v".to_string())]),
        source_format: SourceFormat::SMILES,
    };

    let components = mol.split_components();
    assert_eq!(components.len(), 2);

    let c0 = &components[0];
    assert_eq!(c0.atom_count(), 2);
    assert_eq!(c0.bond_count(), 1);
    assert_eq!(c0.bonds[0].atoms.as_tuple(), (0, 1));
    assert_eq!(c0.bonds[0].donation, Some(BondDonation::Donating));
    assert_eq!(c0.rings.len(), 1);
    assert_eq!(c0.rings[0].start_atom, Some(0));
    assert_eq!(c0.rings[0].end_atom, Some(1));
    assert_eq!(
        c0.positions.as_ref().unwrap()[0],
        Point3D::new(10.0, 0.0, 0.0)
    );
    assert_eq!(
        c0.positions.as_ref().unwrap()[1],
        Point3D::new(11.0, 0.0, 0.0)
    );
    assert_eq!(c0.multicenter_bond_count(), 0);
    assert_eq!(c0.comments, mol.comments);
    assert_eq!(c0.properties, mol.properties);
    assert_eq!(c0.source_format, mol.source_format);

    let c1 = &components[1];
    assert_eq!(c1.atom_count(), 2);
    assert_eq!(c1.bond_count(), 1);
    assert_eq!(c1.bonds[0].atoms.as_tuple(), (0, 1));
    assert_eq!(c1.bonds[0].order, BondOrder::Double);
    assert!(c1.rings.is_empty());
    assert_eq!(
        c1.positions.as_ref().unwrap()[0],
        Point3D::new(20.0, 0.0, 0.0)
    );
    assert_eq!(
        c1.positions.as_ref().unwrap()[1],
        Point3D::new(21.0, 0.0, 0.0)
    );
    assert_eq!(c1.multicenter_bond_count(), 1);
    assert_eq!(
        c1.multicenter_bonds[0].contributions()[0],
        MulticenterSet::new(vec![0, 1], 2)
    );
}

#[test]
fn test_molecule_join_components_empty() {
    let combined = Molecule::join_components(&[]);
    assert_eq!(combined, Molecule::empty());
}

#[test]
fn test_molecule_split_components_roundtrip() {
    let mut dative = Bond::new_dative(0, 1, BondOrder::Single, BondDonation::Donating);
    dative.ring = Some(7);
    let mol = Molecule {
        atoms: vec![
            Atom::aliphatic_atom(Element::C),
            Atom::aliphatic_atom(Element::N),
            Atom::aliphatic_atom(Element::O),
            Atom::aliphatic_atom(Element::H),
        ],
        bonds: vec![dative, Bond::new(2, 3, BondOrder::Double)],
        rings: vec![Ring {
            ring_idx: 7,
            start_atom: Some(0),
            end_atom: Some(1),
            open_span: None,
            close_span: None,
        }],
        positions: Some(vec![
            Point3D::new(10.0, 0.0, 0.0),
            Point3D::new(11.0, 0.0, 0.0),
            Point3D::new(20.0, 0.0, 0.0),
            Point3D::new(21.0, 0.0, 0.0),
        ]),
        multicenter_bonds: vec![MulticenterBond::new(vec![MulticenterSet::new(
            vec![2, 3],
            2,
        )])],
        stereo_interpretation: None,
        comments: vec!["molecule comment".to_string()],
        properties: IndexMap::from_iter([("k".to_string(), "v".to_string())]),
        source_format: SourceFormat::SMILES,
    };

    let split = mol.split_components();
    let combined = Molecule::join_components(&split);

    assert_eq!(combined.atoms, mol.atoms);
    assert_eq!(combined.bonds, mol.bonds);
    assert_eq!(combined.rings, mol.rings);
    assert_eq!(combined.positions, mol.positions);
    assert_eq!(combined.multicenter_bonds, mol.multicenter_bonds);
    assert_eq!(combined.source_format, mol.source_format);
    assert_eq!(combined.stereo_interpretation, mol.stereo_interpretation);
    assert_eq!(combined.properties, mol.properties);
    assert_eq!(
        combined.comments,
        vec!["molecule comment", "molecule comment"]
    );
}

#[test]
fn test_molecule_join_components_metadata() {
    let first = Molecule {
        atoms: vec![Atom::aliphatic_atom(Element::C)],
        bonds: vec![],
        rings: vec![],
        positions: None,
        multicenter_bonds: vec![],
        stereo_interpretation: Some(StereoInterpretation::Absolute),
        comments: vec!["first".to_string()],
        properties: IndexMap::from_iter([
            ("k".to_string(), "v1".to_string()),
            ("a".to_string(), "1".to_string()),
        ]),
        source_format: SourceFormat::MOL,
    };
    let second = Molecule {
        atoms: vec![Atom::aliphatic_atom(Element::O)],
        bonds: vec![],
        rings: vec![],
        positions: None,
        multicenter_bonds: vec![],
        stereo_interpretation: None,
        comments: vec!["second".to_string()],
        properties: IndexMap::from_iter([
            ("k".to_string(), "v2".to_string()),
            ("b".to_string(), "2".to_string()),
        ]),
        source_format: SourceFormat::SMILES,
    };

    let combined = Molecule::join_components(&[first, second]);
    assert_eq!(combined.comments, vec!["first", "second"]);
    assert_eq!(combined.properties.get("k"), Some(&"v2".to_string()));
    assert_eq!(combined.properties.get("a"), Some(&"1".to_string()));
    assert_eq!(combined.properties.get("b"), Some(&"2".to_string()));
    assert_eq!(combined.source_format, SourceFormat::MOL);
    assert_eq!(
        combined.stereo_interpretation,
        Some(StereoInterpretation::Absolute)
    );
}

#[test]
fn test_molecule_join_components_partial_positions() {
    let first = Molecule {
        atoms: vec![Atom::aliphatic_atom(Element::C)],
        bonds: vec![],
        rings: vec![],
        positions: Some(vec![Point3D::new(0.0, 0.0, 0.0)]),
        multicenter_bonds: vec![],
        stereo_interpretation: None,
        comments: vec![],
        properties: IndexMap::new(),
        source_format: SourceFormat::SMILES,
    };
    let second = Molecule {
        atoms: vec![Atom::aliphatic_atom(Element::N)],
        bonds: vec![],
        rings: vec![],
        positions: None,
        multicenter_bonds: vec![],
        stereo_interpretation: None,
        comments: vec![],
        properties: IndexMap::new(),
        source_format: SourceFormat::SMILES,
    };

    let combined = Molecule::join_components(&[first, second]);
    assert!(combined.positions.is_none());
}

#[test]
fn test_molecule_sum_formula() {
    let mol = Molecule {
        atoms: vec![Atom::aliphatic_atom(Element::C)],
        bonds: vec![],
        rings: vec![],
        positions: None,
        multicenter_bonds: vec![],
        stereo_interpretation: None,
        comments: vec![],
        properties: IndexMap::new(),
        source_format: SourceFormat::SMILES,
    };
    assert_eq!(mol.sum_formula(), "C");
}

#[test]
fn test_molecule_sum_formula_hydrogen() {
    let mol = Molecule {
        atoms: vec![
            Atom::aliphatic_atom(Element::C),
            Atom::aliphatic_atom(Element::H),
            Atom::aliphatic_atom(Element::H),
            Atom::aliphatic_atom(Element::H),
            Atom::aliphatic_atom(Element::H),
        ],
        bonds: vec![],
        rings: vec![],
        positions: None,
        multicenter_bonds: vec![],
        stereo_interpretation: None,
        comments: vec![],
        properties: IndexMap::new(),
        source_format: SourceFormat::SMILES,
    };
    assert_eq!(mol.sum_formula(), "CH4");
}

#[test]
fn test_molecule_sum_formula_charge() {
    let mut atom = Atom::aliphatic_atom(Element::C);
    atom.charge = Some(1);
    let mol = Molecule {
        atoms: vec![atom],
        bonds: vec![],
        rings: vec![],
        positions: None,
        multicenter_bonds: vec![],
        stereo_interpretation: None,
        comments: vec![],
        properties: IndexMap::new(),
        source_format: SourceFormat::SMILES,
    };
    let formula = mol.sum_formula();
    assert_eq!(formula, "C+");
}

#[test]
fn test_molecule_sum_formula_multiple_elements() {
    let mol = Molecule {
        atoms: vec![
            Atom::aliphatic_atom(Element::C),
            Atom::aliphatic_atom(Element::C),
            Atom::aliphatic_atom(Element::O),
            Atom::aliphatic_atom(Element::N),
        ],
        bonds: vec![],
        rings: vec![],
        positions: None,
        multicenter_bonds: vec![],
        stereo_interpretation: None,
        comments: vec![],
        properties: IndexMap::new(),
        source_format: SourceFormat::SMILES,
    };
    let formula = mol.sum_formula();
    assert_eq!(formula, "C2NO");
}

#[test]
fn test_extended_molecule_empty() {
    let ext = ExtendedMolecule::empty();
    assert!(ext.atoms.is_empty());
    assert!(ext.bonds.is_empty());
    assert!(ext.rings.is_empty());
    assert!(ext.positions.is_none());
    assert!(ext.multicenter_bonds.is_empty());
    assert!(ext.stereo_interpretation.is_none());
    assert!(ext.properties.is_empty());
    assert!(ext.comments.is_empty());
    assert_eq!(ext.source_format, SourceFormat::UNKNOWN);
    assert!(ext.ctfile_data.is_none());
}

#[test]
fn test_extended_molecule() {
    let ext = ExtendedMolecule {
        atoms: vec![ExtendedAtom::from_element(Element::C)],
        bonds: vec![ExtendedBond::new(0, 0, BondOrder::Single)],
        rings: vec![],
        positions: None,
        multicenter_bonds: vec![],
        stereo_interpretation: None,
        comments: vec!["test".to_string()],
        properties: IndexMap::new(),
        ctfile_data: None,
        cx_data: None,
        source_format: SourceFormat::SMILES,
    };

    assert_eq!(ext.atom_count(), 1);
    assert_eq!(ext.bond_count(), 1);
    assert_eq!(ext.comments.len(), 1);
}

#[test]
fn test_extended_molecule_from_molecule() {
    let mol = Molecule {
        atoms: vec![
            Atom::aliphatic_atom(Element::C),
            Atom::aromatic_atom(Element::N),
        ],
        bonds: vec![Bond::new(0, 1, BondOrder::Double)],
        rings: vec![],
        positions: None,
        multicenter_bonds: vec![],
        stereo_interpretation: None,
        comments: vec![],
        properties: IndexMap::new(),
        source_format: SourceFormat::MOL,
    };

    let ext = ExtendedMolecule::from(mol);
    assert_eq!(ext.atoms.len(), 2);
    assert_eq!(ext.bonds.len(), 1);
    assert_eq!(ext.source_format, SourceFormat::MOL);
    assert!(ext.sgroups().is_empty());
    assert!(ext.rgroups().is_empty());
    assert!(ext.ctfile_data.is_none());
}

#[test]
fn test_extended_molecule_atom_count() {
    let mut ext = ExtendedMolecule::empty();
    ext.atoms.push(ExtendedAtom::from_element(Element::C));
    ext.atoms.push(ExtendedAtom::from_element(Element::N));
    assert_eq!(ext.atom_count(), 2);
}

#[test]
fn test_extended_molecule_bond_count() {
    let mut ext = ExtendedMolecule::empty();
    ext.bonds.push(ExtendedBond::new(0, 1, BondOrder::Single));
    ext.bonds.push(ExtendedBond::new(1, 2, BondOrder::Double));
    assert_eq!(ext.bond_count(), 2);
}

#[test]
fn test_extended_molecule_multicenter_bond_count() {
    let mut ext = ExtendedMolecule::empty();
    ext.multicenter_bonds
        .push(MulticenterBond::new(vec![MulticenterSet::new(
            vec![0, 1, 2],
            3,
        )]));
    ext.multicenter_bonds
        .push(MulticenterBond::new(vec![MulticenterSet::new(
            vec![3, 4, 5],
            2,
        )]));
    assert_eq!(ext.multicenter_bond_count(), 2);
}

#[test]
fn test_extended_molecule_split_join_components_core() {
    let mut ext = ExtendedMolecule::empty();
    ext.atoms = vec![
        ExtendedAtom::from_element(Element::C),
        ExtendedAtom::from_element(Element::N),
        ExtendedAtom::from_element(Element::O),
        ExtendedAtom::from_element(Element::H),
    ];
    ext.bonds = vec![
        ExtendedBond::new_dative(0, 1, BondOrder::Single, BondDonation::Donating),
        ExtendedBond::new(2, 3, BondOrder::Double),
    ];
    ext.rings = vec![Ring {
        ring_idx: 5,
        start_atom: Some(2),
        end_atom: Some(3),
        open_span: None,
        close_span: None,
    }];
    ext.positions = Some(vec![
        Point3D::new(0.0, 0.0, 0.0),
        Point3D::new(1.0, 0.0, 0.0),
        Point3D::new(2.0, 0.0, 0.0),
        Point3D::new(3.0, 0.0, 0.0),
    ]);
    ext.multicenter_bonds = vec![MulticenterBond::new(vec![MulticenterSet::new(
        vec![2, 3],
        2,
    )])];

    let split = ext.split_components();
    assert_eq!(split.len(), 2);
    let joined = ExtendedMolecule::join_components(&split);
    assert_eq!(joined.atoms, ext.atoms);
    assert_eq!(joined.bonds, ext.bonds);
    assert_eq!(joined.rings, ext.rings);
    assert_eq!(joined.positions, ext.positions);
    assert_eq!(joined.multicenter_bonds, ext.multicenter_bonds);
}

#[test]
fn test_extended_molecule_split_join_ctfile_data() {
    let mut ext = ExtendedMolecule::empty();
    ext.atoms = vec![
        ExtendedAtom::from_element(Element::C),
        ExtendedAtom::from_element(Element::N),
        ExtendedAtom::from_element(Element::O),
        ExtendedAtom::from_element(Element::H),
    ];
    ext.bonds = vec![
        ExtendedBond::new(0, 1, BondOrder::Single),
        ExtendedBond::new(2, 3, BondOrder::Single),
    ];
    let mut sgroup_a = SGroup::new(SGroupType::Component);
    sgroup_a.atom_indices = vec![0, 1];
    sgroup_a.bond_indices = vec![0];
    let mut sgroup_b = SGroup::new(SGroupType::Component);
    sgroup_b.atom_indices = vec![2, 3];
    sgroup_b.bond_indices = vec![1];
    ext.ctfile_data = Some(CtfileData {
        sgroups: vec![(10u32, sgroup_a), (20u32, sgroup_b)]
            .into_iter()
            .collect(),
        rgroups: vec![(1u32, RGroup::new(Some(1)))].into_iter().collect(),
        legacy_group_abbreviations: vec![
            LegacyGroupAbbreviation {
                atom_index1: 0,
                atom_index2: 1,
                label: "A".to_string(),
            },
            LegacyGroupAbbreviation {
                atom_index1: 2,
                atom_index2: 3,
                label: "B".to_string(),
            },
        ],
    });

    let split = ext.split_components();
    assert_eq!(split.len(), 2);

    let joined = ExtendedMolecule::join_components(&split);
    assert!(joined.ctfile_data.is_some());
    assert_eq!(joined.ctfile_data.as_ref().unwrap().sgroups.len(), 2);
}

#[test]
fn test_extended_molecule_split_join_cx_data() {
    let mut ext = ExtendedMolecule::empty();
    ext.atoms = vec![
        ExtendedAtom::from_element(Element::C),
        ExtendedAtom::from_element(Element::N),
        ExtendedAtom::from_element(Element::O),
        ExtendedAtom::from_element(Element::H),
    ];
    ext.bonds = vec![
        ExtendedBond::new(0, 1, BondOrder::Single),
        ExtendedBond::new(2, 3, BondOrder::Single),
    ];
    let mut sgroup = SGroup::new(SGroupType::Component);
    sgroup.atom_indices = vec![2, 3];
    sgroup.bond_indices = vec![1];
    ext.cx_data = Some(CxAnnotationData {
        stereo_groups: vec![(
            1u32,
            StereoSet {
                atoms: vec![2, 3],
                mode: StereoSetMode::Correlated,
            },
        )]
        .into_iter()
        .collect(),
        components: Some(vec![vec![2, 3]]),
        sgroups: vec![(7u32, sgroup)].into_iter().collect(),
        rgroups: vec![(9u32, RGroup::new(Some(9)))].into_iter().collect(),
        rgroup_members: vec![(9u32, vec!["[*:1]C".to_string()])]
            .into_iter()
            .collect(),
        local_parity: Some(vec![LocalParityCenter {
            center: 2,
            substituents: vec![3],
            chirality: Chirality::Clockwise,
        }]),
        bicyclo_stereo: Some(vec![BicycloStereo::TowardsEitherBridge(
            BicycloStereoData {
                ligand_atom: 2,
                connection_atom: 3,
                lower_bridge_atoms: vec![2],
                higher_bridge_atoms: vec![3],
            },
        )]),
    });

    let split = ext.split_components();
    assert_eq!(
        split[1].cx_data.as_ref().unwrap().components,
        Some(vec![vec![0, 1]])
    );
    assert_eq!(
        split[1]
            .cx_data
            .as_ref()
            .unwrap()
            .stereo_groups
            .get(&1)
            .unwrap()
            .atoms,
        vec![0, 1]
    );

    let joined = ExtendedMolecule::join_components(&split);
    assert!(joined.cx_data.is_some());
}

#[test]
fn test_extended_molecule_try_join_components_collision() {
    let mut a = ExtendedMolecule::empty();
    a.atoms = vec![ExtendedAtom::from_element(Element::C)];
    a.ctfile_data = Some(CtfileData {
        sgroups: vec![(1u32, SGroup::new(SGroupType::Component))]
            .into_iter()
            .collect(),
        rgroups: Default::default(),
        legacy_group_abbreviations: vec![],
    });

    let mut b = ExtendedMolecule::empty();
    b.atoms = vec![ExtendedAtom::from_element(Element::N)];
    b.ctfile_data = Some(CtfileData {
        sgroups: vec![(1u32, SGroup::new(SGroupType::Component))]
            .into_iter()
            .collect(),
        rgroups: Default::default(),
        legacy_group_abbreviations: vec![],
    });

    let err = ExtendedMolecule::try_join_components(&[a, b]).unwrap_err();
    assert_eq!(err, JoinError::CtfileSgroupCollision { label: 1 });
}

#[test]
fn test_extended_molecule_sum_formula() {
    let mut ext = ExtendedMolecule::empty();
    ext.atoms.push(ExtendedAtom::from_element(Element::C));
    ext.atoms.push(ExtendedAtom::from_element(Element::H));
    ext.atoms.push(ExtendedAtom::from_element(Element::H));
    ext.atoms.push(ExtendedAtom::from_element(Element::H));
    ext.atoms.push(ExtendedAtom::from_element(Element::H));
    assert_eq!(ext.sum_formula(), "CH4");
}

#[test]
fn test_extended_molecule_sum_formula_with_wildcards() {
    use super::super::atom::{AtomSymbol, WildcardAtom};

    let mut ext = ExtendedMolecule::empty();
    // *C(=O)OCC - ethyl ester with wildcard
    ext.atoms
        .push(ExtendedAtom::from_atom_symbol(AtomSymbol::WildcardAtom(
            WildcardAtom::Any,
        )));
    ext.atoms.push(ExtendedAtom::from_element(Element::C));
    ext.atoms.push(ExtendedAtom::from_element(Element::O));
    ext.atoms.push(ExtendedAtom::from_element(Element::O));
    ext.atoms.push(ExtendedAtom::from_element(Element::C));
    ext.atoms.push(ExtendedAtom::from_element(Element::C));
    assert_eq!(ext.sum_formula(), "C3O2*");

    // Two wildcards
    let mut ext2 = ExtendedMolecule::empty();
    ext2.atoms
        .push(ExtendedAtom::from_atom_symbol(AtomSymbol::WildcardAtom(
            WildcardAtom::Any,
        )));
    ext2.atoms.push(ExtendedAtom::from_element(Element::C));
    ext2.atoms
        .push(ExtendedAtom::from_atom_symbol(AtomSymbol::WildcardAtom(
            WildcardAtom::Any,
        )));
    assert_eq!(ext2.sum_formula(), "C*2");
}

#[test]
fn test_extended_molecule_sgroups() {
    let mut ext = ExtendedMolecule::empty();
    assert!(ext.sgroups().is_empty());

    ext.sgroups_mut()
        .insert(0, SGroup::new(SGroupType::Superatom));
    assert_eq!(ext.sgroups().len(), 1);
    assert!(ext.sgroups().contains_key(&0));
}

#[test]
fn test_extended_molecule_rgroups() {
    let mut ext = ExtendedMolecule::empty();
    assert!(ext.rgroups().is_empty());

    let rgroup = RGroup {
        label: Some(1),
        dependent_label: None,
        rgroup_or_h: false,
        occurrence: vec![RGroupOccurrence::Range(1, 3)],
    };
    ext.rgroups_mut().insert(1, rgroup);
    assert_eq!(ext.rgroups().len(), 1);
    assert!(ext.rgroups().contains_key(&1));
}

#[test]
fn test_extended_molecule_to_molecule() {
    let mut ext = ExtendedMolecule::empty();
    ext.atoms.push(ExtendedAtom::from_element(Element::C));
    ext.atoms.push(ExtendedAtom::from_element(Element::O));
    ext.bonds.push(ExtendedBond::new(0, 1, BondOrder::Single));
    ext.source_format = SourceFormat::SMILES;

    let mol = ext.to_molecule().unwrap();
    assert_eq!(mol.atoms.len(), 2);
    assert_eq!(mol.bonds.len(), 1);
    assert_eq!(mol.source_format, SourceFormat::SMILES);
}

#[test]
fn test_extended_molecule_to_molecule_error() {
    let mut ext = ExtendedMolecule::empty();
    let mut atom = ExtendedAtom::from_element(Element::C);
    atom.stereo_care = Some(AtomStereoCare::Care);
    ext.atoms.push(atom);

    let result = ext.to_molecule();
    assert!(result.is_err());
    assert_eq!(result.unwrap_err(), ConversionError::HasExtendedFeatures);
}
