use super::*;
use crate::atom::{Atom, AtomSymbol};
use crate::bond::{Bond, BondType};
use crate::sgroup::{SGroup, SGroupType};
use rstest::rstest;
use umol::error::{DataError, Error, ValidationError};
use umol_data::{e, Element};

/// Helper to create a test molecule with specified number of atoms and bonds
fn create_test_molecule(atom_count: usize, bond_count: usize) -> Molecule {
    let mut molecule = Molecule::new();

    // Add atoms
    for _i in 0..atom_count {
        let atom = Atom::new(AtomSymbol::Element(e!(C)));
        molecule.add_atom(atom);
    }

    // Add bonds (connect consecutive atoms)
    for i in 0..bond_count.min(atom_count.saturating_sub(1)) {
        molecule.add_bond(i, i + 1, Bond::new(BondType::Single));
    }

    molecule
}

/// Helper to create a test molecule with pre-set properties for conflict testing
fn create_molecule_with_properties() -> Molecule {
    let mut molecule = create_test_molecule(3, 2);

    // Set some initial properties to test conflicts
    if let Some(atom) = molecule.atom_mut(0) {
        atom.charge = 1;
        atom.radical = Some(2);
        atom.isotope_mass = Some(13);
        atom.properties
            .insert("molFileAlias".to_string(), "existing".to_string());
        atom.properties
            .insert("molFileValue".to_string(), "existing".to_string());
    }

    // Add an SGroup with existing properties
    let mut sgroup = SGroup::new(0, SGroupType::Superatom);
    sgroup.atom_indices = vec![0, 1];
    sgroup.bond_indices = vec![0];
    sgroup.label = Some("existing".to_string());
    molecule.sgroups.push(sgroup);

    molecule
}

#[rstest]
#[case(vec![ChargeEntry { atom_index: 0, charge: -1 }], 0, -1)]
#[case(vec![ChargeEntry { atom_index: 1, charge: 2 }], 1, 2)]
fn test_apply_charge_entries(
    #[case] entries: Vec<ChargeEntry>,
    #[case] expected_atom_index: usize,
    #[case] expected_charge: i8,
) {
    let mut molecule = create_test_molecule(3, 0);

    // Apply charges
    entries.apply(&mut molecule).unwrap();

    // Verify charge was applied
    assert_eq!(
        molecule.atom(expected_atom_index).unwrap().charge,
        expected_charge
    );
}

#[test]
fn test_apply_charge_entries_multiple() {
    let mut molecule = create_test_molecule(3, 0);
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

    entries.apply(&mut molecule).unwrap();

    assert_eq!(molecule.atom(0).unwrap().charge, -1);
    assert_eq!(molecule.atom(1).unwrap().charge, 0); // unchanged
    assert_eq!(molecule.atom(2).unwrap().charge, 1);
}

#[test]
fn test_apply_charge_entries_empty() {
    let mut molecule = create_test_molecule(2, 0);
    let entries: Vec<ChargeEntry> = vec![];

    // Should succeed with no changes
    entries.apply(&mut molecule).unwrap();

    assert_eq!(molecule.atom(0).unwrap().charge, 0);
    assert_eq!(molecule.atom(1).unwrap().charge, 0);
}

#[test]
fn test_charge_entry_out_of_bounds() {
    let mut molecule = create_test_molecule(2, 0);
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

#[test]
fn test_charge_entry_conflict() {
    let mut molecule = create_molecule_with_properties();
    let entries = vec![ChargeEntry {
        atom_index: 0,
        charge: -1,
    }]; // conflicts with existing charge 1

    let result = entries.apply(&mut molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Validation(ValidationError::InvalidComponent(msg)) => {
            assert!(msg.contains("Charge conflict"));
            assert!(msg.contains("existing 1 vs new -1"));
        }
        _ => panic!("Expected InvalidComponent error"),
    }
}

#[rstest]
#[case(0, None)]
#[case(1, Some(1))]
#[case(2, Some(2))]
#[case(3, Some(3))]
fn test_apply_radical_entry(#[case] radical_type: i8, #[case] expected: Option<u8>) {
    let mut molecule = create_test_molecule(2, 0);
    let entries = vec![RadicalEntry {
        atom_index: 0,
        radical_type,
    }];

    entries.apply(&mut molecule).unwrap();

    assert_eq!(molecule.atom(0).unwrap().radical, expected);
}

#[test]
fn test_radical_entry_invalid() {
    let mut molecule = create_test_molecule(2, 0);
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

#[test]
fn test_apply_isotope_entry() {
    let mut molecule = create_test_molecule(2, 0);
    let entries = vec![IsotopeEntry {
        atom_index: 0,
        mass: 14,
    }];

    entries.apply(&mut molecule).unwrap();

    assert_eq!(molecule.atom(0).unwrap().isotope_mass, Some(14));
}

#[test]
fn test_isotope_entry_conflict() {
    let mut molecule = create_molecule_with_properties();
    let entries = vec![IsotopeEntry {
        atom_index: 0,
        mass: 14,
    }]; // conflicts with existing 13

    let result = entries.apply(&mut molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Validation(ValidationError::InvalidComponent(msg)) => {
            assert!(msg.contains("Isotope conflict"));
            assert!(msg.contains("existing 13 vs new 14"));
        }
        _ => panic!("Expected InvalidComponent error"),
    }
}

#[test]
fn test_apply_sgroup_type_entry() {
    let mut molecule = create_test_molecule(2, 0);
    let entries = vec![SGroupTypeEntry {
        sgroup_index: 0,
        sgroup_type: "SUP".to_string(),
    }];

    entries.apply(&mut molecule).unwrap();

    assert_eq!(molecule.sgroups.len(), 1);
    assert_eq!(molecule.sgroups[0].group_type, SGroupType::Superatom);
}

#[test]
fn test_apply_sgroup_type_entry_multiple() {
    let mut molecule = create_test_molecule(2, 0);
    let entries = vec![SGroupTypeEntry {
        sgroup_index: 2,
        sgroup_type: "DAT".to_string(),
    }];

    entries.apply(&mut molecule).unwrap();

    // Should create SGroups 0, 1, 2
    assert_eq!(molecule.sgroups.len(), 3);
    assert_eq!(molecule.sgroups[0].group_type, SGroupType::Generic);
    assert_eq!(molecule.sgroups[1].group_type, SGroupType::Generic);
    assert_eq!(molecule.sgroups[2].group_type, SGroupType::Data);
}

#[test]
fn test_apply_sgroup_label_entry() {
    let mut molecule = create_test_molecule(2, 0);
    let entries = vec![SGroupLabelEntry {
        sgroup_index: 0,
        label: "Ph".to_string(),
    }];

    entries.apply(&mut molecule).unwrap();

    assert_eq!(molecule.sgroups.len(), 1);
    assert_eq!(molecule.sgroups[0].label, Some("Ph".to_string()));
}

#[test]
fn test_apply_sgroup_atom_list_entry() {
    let mut molecule = create_test_molecule(3, 0);
    let entry = SGroupAtomListEntry {
        sgroup_index: 0,
        atom_indices: vec![0, 1, 2],
    };

    entry.apply(&mut molecule).unwrap();

    assert_eq!(molecule.sgroups.len(), 1);
    assert_eq!(molecule.sgroups[0].atom_indices, vec![0, 1, 2]);
}

#[test]
fn test_sgroup_atom_list_entry_invalid_atom_index() {
    let mut molecule = create_test_molecule(2, 0);
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

#[test]
fn test_apply_sgroup_bond_list_entry() {
    let mut molecule = create_test_molecule(3, 2);
    let entry = SGroupBondListEntry {
        sgroup_index: 0,
        bond_indices: vec![0, 1],
    };

    entry.apply(&mut molecule).unwrap();

    assert_eq!(molecule.sgroups.len(), 1);
    assert_eq!(molecule.sgroups[0].bond_indices, vec![0, 1]);
}

#[test]
fn test_sgroup_bond_list_entry_invalid_bond_index() {
    let mut molecule = create_test_molecule(3, 1);
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

#[test]
fn test_apply_atom_alias_entry() {
    let mut molecule = create_test_molecule(2, 0);
    let entry = AtomAliasEntry {
        atom_index: 0,
        alias: "CF3".to_string(),
    };

    entry.apply(&mut molecule).unwrap();

    let alias = molecule
        .atom(0)
        .unwrap()
        .properties
        .get("molFileAlias")
        .unwrap();
    assert_eq!(alias, "CF3");
}

#[test]
fn test_atom_alias_entry_conflict() {
    let mut molecule = create_molecule_with_properties();
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

#[test]
fn test_apply_atom_value_entry() {
    let mut molecule = create_test_molecule(2, 0);
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
