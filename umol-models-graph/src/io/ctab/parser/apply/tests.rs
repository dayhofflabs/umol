use super::*;
use crate::io::ctab::atom::{Atom, AtomRadical, AtomSymbol};
use crate::io::ctab::bond::{Bond, BondType};
use crate::io::ctab::sgroup::{SGroup, SGroupBracketStyle, SGroupType};
use pretty_assertions::assert_eq;
use rstest::*;
use umol::error::{DataError, Error, ValidationError};
use umol_data::{e, Element};

#[fixture]
fn basic_molecule() -> Molecule {
    let mut molecule = Molecule::new();
    for _ in 0..3 {
        let atom = Atom::new(AtomSymbol::Element(e!(C)));
        molecule.add_atom(atom);
    }
    molecule.add_bond(0, 1, Bond::new(BondType::Single));
    molecule.add_bond(1, 2, Bond::new(BondType::Single));
    molecule
}

#[fixture]
fn molecule_with_properties(basic_molecule: Molecule) -> Molecule {
    let mut molecule = basic_molecule;

    // Set some initial properties to test conflicts
    if let Some(atom) = molecule.atom_mut(0) {
        atom.charge = 1;
        atom.radical = Some(AtomRadical::Doublet);
        atom.isotope_mass = Some(13);
        atom.properties
            .insert("molFileAlias".to_string(), "existing".to_string());
        atom.properties
            .insert("molFileValue".to_string(), "existing".to_string());
    }

    // Add an SGroup with existing properties
    let mut sgroup = SGroup::new(SGroupType::Superatom);
    sgroup.atom_indices = vec![0, 1];
    sgroup.bond_indices = vec![0];
    sgroup.label = None;
    sgroup.bracket_style = Some(SGroupBracketStyle::Default);
    molecule.sgroups.insert(0, sgroup);

    molecule
}

#[fixture]
fn molecule_with_rgroup(basic_molecule: Molecule) -> Molecule {
    let mut molecule = basic_molecule;
    molecule.atom_mut(0).unwrap().symbol = AtomSymbol::RGroup(RGroup::new(None));
    molecule
}

#[fixture]
fn molecule_with_labeled_rgroup(molecule_with_rgroup: Molecule) -> Molecule {
    let mut molecule = molecule_with_rgroup;
    molecule.atom_mut(0).unwrap().symbol = AtomSymbol::RGroup(RGroup::new(Some(1)));
    molecule
}

#[fixture]
fn molecule_with_sgroup(basic_molecule: Molecule) -> Molecule {
    let mut molecule = basic_molecule;
    let mut sgroup = SGroup::new(SGroupType::Generic);
    sgroup.atom_indices = vec![0, 1];
    sgroup.bond_indices = vec![0];
    molecule.add_sgroup(sgroup);
    molecule
}

#[fixture]
fn molecule_with_labeled_sgroup(molecule_with_sgroup: Molecule) -> Molecule {
    let mut molecule = molecule_with_sgroup;
    molecule.sgroups.get_mut(&0).unwrap().label = Some(15);
    let sgroup = SGroup::new(SGroupType::Superatom);
    molecule.add_sgroup(sgroup);
    molecule
}

#[rstest]
fn test_apply_atom_alias_entry(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let entry = AtomAliasEntry {
        atom_index: 0,
        alias: "C1".to_string(),
    };

    entry.apply(&mut molecule).unwrap();

    let alias = molecule
        .atom(0)
        .unwrap()
        .properties
        .get("molFileAlias")
        .unwrap();
    assert_eq!(alias, "C1");
}

#[rstest]
fn test_atom_alias_entry_invalid(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let entry = AtomAliasEntry {
        atom_index: 5,
        alias: "C1".to_string(),
    };

    let result = entry.apply(&mut molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Data(DataError::MissingAtomIndex(idx)) => assert_eq!(idx, 5),
        _ => panic!("Expected MissingAtomIndex error"),
    }
}

#[rstest]
fn test_atom_alias_entry_conflict(molecule_with_properties: Molecule) {
    let mut molecule = molecule_with_properties;
    let entry = AtomAliasEntry {
        atom_index: 0,
        alias: "new".to_string(),
    };

    let result = entry.apply(&mut molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Validation(ValidationError::InvalidComponent(msg)) => {
            assert!(msg.contains("Atom alias conflict"));
            assert!(msg.contains("existing 'existing' vs new 'new'"));
        }
        _ => panic!("Expected InvalidComponent error"),
    }
}

#[rstest]
fn test_apply_atom_value_entry(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let entry = AtomValueEntry {
        atom_index: 0,
        value: "*".to_string(),
    };

    entry.apply(&mut molecule).unwrap();

    let value = molecule
        .atom(0)
        .unwrap()
        .properties
        .get("molFileValue")
        .unwrap();
    assert_eq!(value, "*");
}

#[rstest]
fn test_apply_atom_value_entry_invalid(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let entry = AtomValueEntry {
        atom_index: 5,
        value: "*".to_string(),
    };

    let result = entry.apply(&mut molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Data(DataError::MissingAtomIndex(idx)) => assert_eq!(idx, 5),
        _ => panic!("Expected MissingAtomIndex error"),
    }
}

#[rstest]
fn test_apply_atom_value_conflict(molecule_with_properties: Molecule) {
    let mut molecule = molecule_with_properties;
    let entry = AtomValueEntry {
        atom_index: 0,
        value: "new".to_string(),
    };

    let result = entry.apply(&mut molecule);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_charge_entries(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let entries = vec![ChargeEntry {
        atom_index: 0,
        charge: -1,
    }];

    entries.apply(&mut molecule).unwrap();
    assert_eq!(molecule.atom(0).unwrap().charge, -1);
}

#[rstest]
fn test_apply_charge_entries_multiple(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let entries = vec![
        ChargeEntry {
            atom_index: 0,
            charge: -1,
        },
        ChargeEntry {
            atom_index: 2,
            charge: 1,
        },
    ];

    let result = entries.apply(&mut molecule);
    assert!(result.is_ok());

    assert_eq!(molecule.atom(0).unwrap().charge, -1);
    assert_eq!(molecule.atom(1).unwrap().charge, 0); // unchanged
    assert_eq!(molecule.atom(2).unwrap().charge, 1);
}

#[rstest]
fn test_apply_charge_entries_empty(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let entries: Vec<ChargeEntry> = vec![];

    // Should succeed with no changes
    entries.apply(&mut molecule).unwrap();

    assert_eq!(molecule.atom(0).unwrap().charge, 0);
    assert_eq!(molecule.atom(1).unwrap().charge, 0);
}

#[rstest]
fn test_apply_charge_entries_out_of_bounds(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let entries = vec![ChargeEntry {
        atom_index: 5,
        charge: -1,
    }];

    let result = entries.apply(&mut molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Data(DataError::MissingAtomIndex(idx)) => assert_eq!(idx, 5),
        _ => panic!("Expected MissingAtomIndex error"),
    }
}

#[rstest]
fn test_apply_charge_entries_overwrite(molecule_with_properties: Molecule) {
    let mut molecule = molecule_with_properties;
    let entries = vec![ChargeEntry {
        atom_index: 0,
        charge: -1,
    }]; // overwrites existing charge 1

    let result = entries.apply(&mut molecule);
    assert!(result.is_ok());
    let atom = molecule.atom(0).unwrap();
    assert_eq!(atom.charge, -1);
    assert_eq!(atom.radical, None);
}

#[rstest]
#[case(0, None)]
#[case(1, Some(AtomRadical::Singlet))]
#[case(2, Some(AtomRadical::Doublet))]
#[case(3, Some(AtomRadical::Triplet))]
fn test_apply_radical_entries(
    basic_molecule: Molecule,
    #[case] radical_type: u8,
    #[case] expected: Option<AtomRadical>,
) {
    let mut molecule = basic_molecule;
    let entries = vec![RadicalEntry {
        atom_index: 0,
        radical_type,
    }];

    let result = entries.apply(&mut molecule);
    assert!(result.is_ok());
    assert_eq!(molecule.atom(0).unwrap().radical, expected);
}

#[rstest]
fn test_apply_radical_entries_overwrite(molecule_with_properties: Molecule) {
    let mut molecule = molecule_with_properties;
    let entries = vec![RadicalEntry {
        atom_index: 0,
        radical_type: 1,
    }]; // overwrites existing radical doublet

    let result = entries.apply(&mut molecule);
    assert!(result.is_ok());

    let atom = molecule.atom(0).unwrap();
    assert_eq!(atom.radical, Some(AtomRadical::Singlet));
    assert_eq!(atom.charge, 0);
}

#[rstest]
fn test_apply_radical_entries_invalid(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let entries = vec![RadicalEntry {
        atom_index: 0,
        radical_type: 4,
    }]; // invalid

    let result = entries.apply(&mut molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Validation(ValidationError::InvalidComponent(msg)) => {
            assert!(msg.contains("Invalid radical type"));
        }
        _ => panic!("Expected InvalidComponent error"),
    }
}

#[rstest]
fn test_apply_isotope_entries(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let entries = vec![IsotopeEntry {
        atom_index: 0,
        mass: 14,
    }];

    let result = entries.apply(&mut molecule);
    assert!(result.is_ok());

    let atom = molecule.atom(0).unwrap();
    assert_eq!(atom.isotope_mass, Some(14));
}

#[rstest]
fn test_apply_isotope_entries_invalid(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let entries = vec![IsotopeEntry {
        atom_index: 0,
        mass: 40,
    }];

    let result = entries.apply(&mut molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Data(DataError::InvalidIsotope(msg)) => {
            assert!(msg.contains("Invalid isotope mass number 40 for element C"));
        }
        _ => panic!("Expected InvalidIsotope error"),
    }
}

#[rstest]
fn test_apply_isotope_entries_conflict(molecule_with_properties: Molecule) {
    let mut molecule = molecule_with_properties;
    let entries = vec![IsotopeEntry {
        atom_index: 0,
        mass: 14,
    }]; // conflicts with existing 13

    let result = entries.apply(&mut molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Validation(ValidationError::InvalidComponent(msg)) => {
            assert!(msg.contains("Isotope conflict for atom 0"));
        }
        _ => panic!("Expected InvalidComponent error"),
    }

    let entries = vec![IsotopeEntry {
        atom_index: 0,
        mass: 0,
    }];

    let result = entries.apply(&mut molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Validation(ValidationError::InvalidComponent(msg)) => {
            assert!(msg.contains("Isotope conflict for atom 0"));
        }
        _ => panic!("Expected InvalidComponent error"),
    }
}

#[rstest]
fn test_apply_ring_bond_count_entries(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let entries = vec![RingBondCountEntry {
        atom_index: 0,
        ring_bond_count: 2,
    }];

    entries.apply(&mut molecule).unwrap();
    assert_eq!(molecule.atom(0).unwrap().ring_bond_count, Some(2));
}

#[rstest]
fn test_apply_ring_bond_count_entries_conflict(molecule_with_properties: Molecule) {
    let mut molecule = molecule_with_properties;
    molecule.atom_mut(0).unwrap().ring_bond_count = Some(3);
    let entries = vec![RingBondCountEntry {
        atom_index: 0,
        ring_bond_count: 2,
    }];

    let result = entries.apply(&mut molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Validation(ValidationError::InvalidComponent(msg)) => {
            assert!(msg.contains("Ring bond count conflict"));
            assert!(msg.contains("existing 1 vs new 2"));
        }
        _ => panic!("Expected InvalidComponent error"),
    }
}

#[rstest]
fn test_apply_substitution_count_entries(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let entries = vec![SubstitutionCountEntry {
        atom_index: 0,
        substitution_count: 2,
    }];

    entries.apply(&mut molecule).unwrap();
    assert_eq!(molecule.atom(0).unwrap().substitution_count, Some(2));
}

#[rstest]
fn test_apply_substitution_count_entries_invalid(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let entries = vec![SubstitutionCountEntry {
        atom_index: 0,
        substitution_count: 7,
    }];

    let result = entries.apply(&mut molecule);
    assert!(result.is_err());
    match result.unwrap_err() {
        Error::Validation(ValidationError::InvalidComponent(msg)) => {
            assert!(msg.contains("Invalid substitution count"));
        }
        _ => panic!("Expected InvalidComponent error"),
    }
}

#[rstest]
fn test_apply_substitution_count_entries_conflict(molecule_with_properties: Molecule) {
    let mut molecule = molecule_with_properties;
    molecule.atom_mut(0).unwrap().substitution_count = Some(1);
    let entries = vec![SubstitutionCountEntry {
        atom_index: 0,
        substitution_count: 2,
    }];

    let result = entries.apply(&mut molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Validation(ValidationError::InvalidComponent(msg)) => {
            assert!(msg.contains("Substitution count conflict"));
            assert!(msg.contains("existing 1 vs new 2"));
        }
        _ => panic!("Expected InvalidComponent error"),
    }
}

#[rstest]
fn test_apply_unsaturated_atom_entries(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let entries = vec![UnsaturatedAtomEntry {
        atom_index: 0,
        unsaturated: 1,
    }];

    entries.apply(&mut molecule).unwrap();
    assert_eq!(molecule.atom(0).unwrap().unsaturated, Some(true));
}

#[rstest]
fn test_apply_unsaturated_atom_entries_conflict(molecule_with_properties: Molecule) {
    let mut molecule = molecule_with_properties;
    molecule.atom_mut(0).unwrap().unsaturated = Some(true);
    let entries = vec![UnsaturatedAtomEntry {
        atom_index: 0,
        unsaturated: 0,
    }];

    let result = entries.apply(&mut molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Validation(ValidationError::InvalidComponent(msg)) => {
            assert!(msg.contains("Unsaturated conflict"));
        }
        _ => panic!("Expected InvalidComponent error"),
    }
}

#[rstest]
fn test_apply_unsaturated_atom_entries_invalid(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let entries = vec![UnsaturatedAtomEntry {
        atom_index: 0,
        unsaturated: 2,
    }];

    let result = entries.apply(&mut molecule);
    assert!(result.is_err());
    match result.unwrap_err() {
        Error::Validation(ValidationError::InvalidComponent(msg)) => {
            assert!(msg.contains("Invalid unsaturated atom"));
        }
        _ => panic!("Expected InvalidComponent error"),
    }
}

#[rstest]
fn test_apply_atom_list_entry(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let entry = AtomListEntry {
        atom_index: 0,
        exclusion: false,
        elements: vec![e!(C), e!(Si)],
    };

    entry.apply(&mut molecule).unwrap();
    assert_eq!(
        molecule.atom(0).unwrap().symbol,
        AtomSymbol::AtomList(AtomList {
            elements: vec![e!(C), e!(Si)],
        })
    );
}

#[rstest]
fn test_apply_atom_list_entry_invalid(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let entry = AtomListEntry {
        atom_index: 5,
        exclusion: false,
        elements: vec![e!(C), e!(Si)],
    };

    let result = entry.apply(&mut molecule);
    assert!(result.is_err());
    match result.unwrap_err() {
        Error::Data(DataError::MissingAtomIndex(idx)) => assert_eq!(idx, 5),
        _ => panic!("Expected MissingAtomIndex error"),
    }
}

#[rstest]
fn test_apply_atom_list_entry_conflict(molecule_with_properties: Molecule) {
    let mut molecule = molecule_with_properties;
    molecule.atom_mut(0).unwrap().symbol = AtomSymbol::AtomList(AtomList {
        elements: vec![e!(C), e!(Si)],
    });

    let entry = AtomListEntry {
        atom_index: 0,
        exclusion: false,
        elements: vec![e!(C), e!(Pb)],
    };

    let result = entry.apply(&mut molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Validation(ValidationError::InvalidComponent(msg)) => {
            assert!(msg.contains("Atom list conflict"));
        }
        _ => panic!("Expected InvalidComponent error"),
    }
}

#[rstest]
fn test_apply_rgroup_label_entries(molecule_with_rgroup: Molecule) {
    let mut molecule = molecule_with_rgroup;
    let entries = vec![RGroupLabelEntry {
        atom_index: 0,
        label: 1,
    }];

    entries.apply(&mut molecule).unwrap();

    assert_eq!(
        molecule.atom(0).unwrap().symbol,
        AtomSymbol::RGroup(RGroup::new(Some(1)))
    );
}

#[rstest]
fn test_apply_rgroup_label_entries_replace(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let entries = vec![RGroupLabelEntry {
        atom_index: 0,
        label: 1,
    }];

    entries.apply(&mut molecule).unwrap();

    assert_eq!(
        molecule.atom(0).unwrap().symbol,
        AtomSymbol::RGroup(RGroup::new(Some(1)))
    );
}

#[rstest]
fn test_apply_rgroup_label_entries_invalid(molecule_with_labeled_rgroup: Molecule) {
    let mut molecule = molecule_with_labeled_rgroup;
    let entries = vec![RGroupLabelEntry {
        atom_index: 5,
        label: 1,
    }];

    let result = entries.apply(&mut molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Data(DataError::MissingAtomIndex(idx)) => assert_eq!(idx, 5),
        _ => panic!("Expected MissingAtomIndex error"),
    }
}

#[rstest]
fn test_apply_rgroup_label_entries_conflict(molecule_with_labeled_rgroup: Molecule) {
    let mut molecule = molecule_with_labeled_rgroup;
    molecule.atom_mut(0).unwrap().symbol = AtomSymbol::RGroup(RGroup::new(Some(2)));
    let entries = vec![RGroupLabelEntry {
        atom_index: 0,
        label: 1,
    }];

    let result = entries.apply(&mut molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Validation(ValidationError::InvalidComponent(msg)) => {
            assert!(msg.contains("RGroup label conflict"));
        }
        _ => panic!("Expected InvalidComponent error"),
    }
}

#[rstest]
fn test_apply_rgroup_label_entries_duplicate(molecule_with_labeled_rgroup: Molecule) {
    let mut molecule = molecule_with_labeled_rgroup;
    let entries = vec![
        RGroupLabelEntry {
            atom_index: 0,
            label: 1,
        },
        RGroupLabelEntry {
            atom_index: 0,
            label: 2,
        },
    ];

    let result = entries.apply(&mut molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Validation(ValidationError::InvalidComponent(msg)) => {
            assert!(msg.contains("RGroup"));
        }
        _ => panic!("Expected InvalidComponent error"),
    }
}

#[rstest]
fn test_apply_sgroup_type_entries(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let entries = vec![SGroupTypeEntry {
        sgroup_index: 1,
        sgroup_type: SGroupType::Superatom,
    }];

    entries.apply(&mut molecule).unwrap();

    assert_eq!(molecule.sgroups.len(), 1);
    assert_eq!(molecule.sgroups[&1].group_type, SGroupType::Superatom);
}

#[rstest]
fn test_apply_sgroup_type_entries_multiple(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let entries = vec![
        SGroupTypeEntry {
            sgroup_index: 0,
            sgroup_type: SGroupType::Data,
        },
        SGroupTypeEntry {
            sgroup_index: 1,
            sgroup_type: SGroupType::Superatom,
        },
        SGroupTypeEntry {
            sgroup_index: 2,
            sgroup_type: SGroupType::RepeatingUnit,
        },
    ];

    entries.apply(&mut molecule).unwrap();

    assert_eq!(molecule.sgroups.len(), 3);
    assert_eq!(molecule.sgroups[&0].group_type, SGroupType::Data);
    assert_eq!(molecule.sgroups[&1].group_type, SGroupType::Superatom);
    assert_eq!(molecule.sgroups[&2].group_type, SGroupType::RepeatingUnit);
}

#[rstest]
fn test_apply_sgroup_type_entries_conflict(molecule_with_sgroup: Molecule) {
    let mut molecule = molecule_with_sgroup;

    // Try to create another SGroup with same index
    let entries = vec![SGroupTypeEntry {
        sgroup_index: 0,
        sgroup_type: SGroupType::Data,
    }];
    let result = entries.apply(&mut molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Validation(ValidationError::InvalidComponent(msg)) => {
            assert!(msg.contains("SGroup index conflict"));
        }
        _ => panic!("Expected InvalidComponent error"),
    }
}

#[rstest]
fn test_apply_sgroup_label_entries(molecule_with_sgroup: Molecule) {
    let mut molecule = molecule_with_sgroup;
    let entries = vec![SGroupLabelEntry {
        sgroup_index: 0,
        label: 19,
    }];

    entries.apply(&mut molecule).unwrap();

    assert_eq!(molecule.sgroups.len(), 1);
    assert_eq!(molecule.sgroups[&0].label, Some(19));
}

#[rstest]
fn test_apply_sgroup_label_entries_invalid(molecule_with_sgroup: Molecule) {
    let mut molecule = molecule_with_sgroup;
    let entries = vec![SGroupLabelEntry {
        sgroup_index: 1,
        label: 15,
    }];

    let result = entries.apply(&mut molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Validation(ValidationError::InvalidComponent(msg)) => {
            assert!(msg.contains("Invalid SGroup index"));
        }
        _ => panic!("Expected InvalidComponent error"),
    }
}

#[rstest]
fn test_apply_sgroup_label_entries_conflict(molecule_with_sgroup: Molecule) {
    let mut molecule = molecule_with_sgroup;
    molecule.sgroups.get_mut(&0).unwrap().label = Some(19);
    let entries = vec![SGroupLabelEntry {
        sgroup_index: 0,
        label: 20,
    }];

    let result = entries.apply(&mut molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Validation(ValidationError::InvalidComponent(msg)) => {
            assert!(msg.contains("SGroup label conflict"));
        }
        _ => panic!("Expected InvalidComponent error"),
    }
}

#[rstest]
fn test_apply_sgroup_label_entries_duplicate(molecule_with_labeled_sgroup: Molecule) {
    let mut molecule = molecule_with_labeled_sgroup;
    let entries = vec![SGroupLabelEntry {
        sgroup_index: 1,
        label: 15,
    }];

    let result = entries.apply(&mut molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Validation(ValidationError::InvalidComponent(msg)) => {
            assert!(msg.contains("SGroup label conflict: duplicate label '15'"));
        }
        _ => panic!("Expected InvalidComponent error"),
    }
}

#[rstest]
fn test_apply_sgroup_atom_list_entry(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let sgroup = SGroup::new(SGroupType::Superatom);
    molecule.sgroups.insert(0, sgroup);

    let entry = SGroupAtomListEntry {
        sgroup_index: 0,
        atom_indices: vec![0, 1, 2],
    };

    entry.apply(&mut molecule).unwrap();

    assert_eq!(molecule.sgroups.len(), 1);
    assert_eq!(molecule.sgroups[&0].atom_indices, vec![0, 1, 2]);
}

#[rstest]
fn test_apply_sgroup_atom_list_entry_invalid(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let sgroup = SGroup::new(SGroupType::Superatom);
    molecule.sgroups.insert(0, sgroup);
    let entry = SGroupAtomListEntry {
        sgroup_index: 0,
        atom_indices: vec![0, 5], // atom 5 doesn't exist
    };

    let result = entry.apply(&mut molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Data(DataError::MissingAtomIndex(idx)) => assert_eq!(idx, 5),
        _ => panic!("Expected MissingAtomIndex error"),
    }
}

#[rstest]
fn test_apply_sgroup_bond_list_entry(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let sgroup = SGroup::new(SGroupType::Superatom);
    molecule.sgroups.insert(0, sgroup);
    let entry = SGroupBondListEntry {
        sgroup_index: 0,
        bond_indices: vec![0, 1],
    };

    entry.apply(&mut molecule).unwrap();

    assert_eq!(molecule.sgroups.len(), 1);
    assert_eq!(molecule.sgroups[&0].bond_indices, vec![0, 1]);
}

#[rstest]
fn test_apply_sgroup_bond_list_entry_invalid(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let sgroup = SGroup::new(SGroupType::Superatom);
    molecule.sgroups.insert(0, sgroup);
    let entry = SGroupBondListEntry {
        sgroup_index: 0,
        bond_indices: vec![0, 5], // bond 5 doesn't exist
    };

    let result = entry.apply(&mut molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Data(DataError::MissingBondIndex(idx)) => assert_eq!(idx, 5),
        _ => panic!("Expected MissingBondIndex error"),
    }
}

#[rstest]
fn test_apply_sgroup_atom_bond_list_entries_shared(basic_molecule: Molecule) {
    let mut molecule = basic_molecule;
    let sgroup = SGroup::new(SGroupType::Superatom);
    molecule.sgroups.insert(0, sgroup);
    let sgroup = SGroup::new(SGroupType::Superatom);
    molecule.sgroups.insert(1, sgroup);

    // Assign overlapping atoms to both SGroups
    let entry1 = SGroupAtomListEntry {
        sgroup_index: 0,
        atom_indices: vec![0, 1],
    };
    entry1.apply(&mut molecule).unwrap();

    let entry2 = SGroupAtomListEntry {
        sgroup_index: 1,
        atom_indices: vec![1, 2],
    };
    entry2.apply(&mut molecule).unwrap();

    // Assign same bond to both SGroups
    let entry3 = SGroupBondListEntry {
        sgroup_index: 0,
        bond_indices: vec![0],
    };
    entry3.apply(&mut molecule).unwrap();

    let entry4 = SGroupBondListEntry {
        sgroup_index: 1,
        bond_indices: vec![0],
    };
    entry4.apply(&mut molecule).unwrap();

    // Verify overlapping assignments
    assert_eq!(molecule.sgroups[&0].atom_indices, vec![0, 1]);
    assert_eq!(molecule.sgroups[&1].atom_indices, vec![1, 2]);
    assert_eq!(molecule.sgroups[&0].bond_indices, vec![0]);
    assert_eq!(molecule.sgroups[&1].bond_indices, vec![0]);
}
