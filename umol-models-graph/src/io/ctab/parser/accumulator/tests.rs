use super::*;
use crate::io::ctab::{
    atom::{
        Atom, AtomList, AtomRadical, AtomStandard, AtomSymbol, AttachmentPointType, LinkAtom,
        UnsaturatedAtom,
    },
    bond::{Bond, BondStandard, BondType},
    molecule::MoleculeStandard,
    parser::properties::{
        AtomAliasEntry, AtomAttachmentOrderEntry, AtomHydrogenCountEntry, AtomListEntry,
        AtomValueEntry, AttachmentPointEntry, ChargeEntry, IsotopeEntry, LinkAtomEntry,
        PropertyEntries, RGroupLabelEntry, RGroupLogicEntry, RadicalEntry, RingBondCountEntry,
        SGroupAtomListEntry, SGroupBondListEntry, SGroupComponentEntry, SGroupConnectingBondEntry,
        SGroupConnectivityEntry, SGroupCorrespondenceEntry, SGroupDataDescriptionEntry,
        SGroupDataDisplayEntry, SGroupDataEntry, SGroupDisplayInfoEntry, SGroupExpansionEntry,
        SGroupHierarchyEntry, SGroupLabelEntry, SGroupParentAtomEntry, SGroupSubscriptData,
        SGroupSubscriptEntry, SGroupSubtypeEntry, SGroupTypeEntry, SubstitutionCountEntry,
        UnsaturatedAtomEntry, ZeroAtomChargeEntry, ZeroBondOrderEntry,
    },
    rgroup::{RGroup, RGroupOccurrence},
    sgroup::{
        SGroupBracketCoords, SGroupConnectingBond, SGroupConnectivity, SGroupDataDisplayChars,
        SGroupDataDisplayPlacement, SGroupDataDisplayType, SGroupDataDisplayUnits, SGroupDataType,
        SGroupSubtype, SGroupType,
    },
};
use pretty_assertions::assert_eq;
use rstest::*;
use umol::error::{DataError, Error, ValidationError};
use umol_data::{e, Element, NamedIsotope};

#[fixture]
fn basic_molecule() -> Molecule {
    let mut molecule = Molecule::new();
    molecule.add_atom(Atom::new(AtomSymbol::Element(e!(C))));
    molecule
}

#[fixture]
fn molecule_with_properties(mut basic_molecule: Molecule) -> Molecule {
    if let Some(atom) = basic_molecule.atom_mut(0) {
        atom.properties
            .insert("molFileAlias".to_string(), "existing".to_string());
        atom.properties
            .insert("molFileValue".to_string(), "existing".to_string());
        atom.charge = 1;
        atom.radical = Some(AtomRadical::Doublet);
        atom.isotope_mass = Some(13);
        atom.ring_bond_count = Some(RingBondCount::R2);
        atom.substitution_count = Some(SubstitutionCount::S2);
        atom.unsaturated = Some(UnsaturatedAtom);
    }
    basic_molecule
}

#[fixture]
fn molecule_with_rgroup(mut basic_molecule: Molecule) -> Molecule {
    basic_molecule.atom_mut(0).unwrap().symbol = AtomSymbol::RGroup(RGroup::new(Some(1)));
    basic_molecule
}

#[fixture]
fn molecule_with_bond() -> Molecule {
    let mut molecule = Molecule::new();
    molecule.add_atom(Atom::new(AtomSymbol::Element(e!(C))));
    molecule.add_atom(Atom::new(AtomSymbol::Element(e!(C))));
    molecule.add_bond(0, 1, Bond::new(BondType::Single));
    molecule
}

#[fixture]
fn basic_molecule_standard() -> MoleculeStandard {
    let mut molecule = MoleculeStandard::new();
    molecule.add_atom(AtomStandard::new(e!(C)));
    molecule
}

#[fixture]
fn molecule_with_bond_standard() -> MoleculeStandard {
    let mut molecule = MoleculeStandard::new();
    molecule.add_atom(AtomStandard::new(e!(C)));
    molecule.add_atom(AtomStandard::new(e!(C)));
    molecule.add_bond(0, 1, BondStandard::new(BondType::Single));
    molecule
}

#[fixture]
fn acc_with_superatom_sgroup() -> MoleculeProperties {
    let mut acc = MoleculeProperties::new();
    let type_entry = PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry {
        sgroup_index: 0,
        sgroup_type: SGroupType::Superatom,
    }]);
    acc.add_entry(type_entry).unwrap();
    acc
}

#[fixture]
fn acc_with_data_sgroup() -> MoleculeProperties {
    let mut acc = MoleculeProperties::new();
    let type_entry = PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry {
        sgroup_index: 0,
        sgroup_type: SGroupType::Data,
    }]);
    acc.add_entry(type_entry).unwrap();
    let data_description_entry =
        PropertyEntries::SGroupDataDescriptionEntry(SGroupDataDescriptionEntry {
            sgroup_index: 0,
            field_name: "test".to_string(),
            field_type: SGroupDataType::Text,
            field_units: None,
            query_identifier: None,
            data_query_operator: None,
        });
    acc.add_entry(data_description_entry).unwrap();
    acc
}

#[fixture]
fn acc_with_multiple_sgroup() -> MoleculeProperties {
    let mut acc = MoleculeProperties::new();
    let type_entry = PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry {
        sgroup_index: 0,
        sgroup_type: SGroupType::MultipleGroup,
    }]);
    acc.add_entry(type_entry).unwrap();
    acc
}

#[fixture]
fn acc_with_copolymer_sgroup() -> MoleculeProperties {
    let mut acc = MoleculeProperties::new();
    let type_entry = PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry {
        sgroup_index: 0,
        sgroup_type: SGroupType::Copolymer,
    }]);
    acc.add_entry(type_entry).unwrap();
    acc
}

#[rstest]
fn test_apply_atom_alias_entry(mut basic_molecule: Molecule) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::AtomAliasEntry(AtomAliasEntry {
        atom_index: 0,
        alias: "C1".to_string(),
    });
    acc.add_entry(entry).unwrap();
    acc.apply(&mut basic_molecule).unwrap();

    let alias = basic_molecule
        .atom(0)
        .unwrap()
        .properties
        .get("molFileAlias")
        .unwrap();
    assert_eq!(alias, "C1");
}

#[rstest]
fn test_atom_alias_entry_invalid(mut basic_molecule: Molecule) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::AtomAliasEntry(AtomAliasEntry {
        atom_index: 5,
        alias: "C1".to_string(),
    });
    acc.add_entry(entry).unwrap();

    let result = acc.apply(&mut basic_molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Data(DataError::MissingAtomIndex(idx)) => assert_eq!(idx, 5),
        _ => panic!("Expected MissingAtomIndex error"),
    }
}

#[rstest]
fn test_atom_alias_entry_conflict(mut molecule_with_properties: Molecule) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::AtomAliasEntry(AtomAliasEntry {
        atom_index: 0,
        alias: "new".to_string(),
    });
    acc.add_entry(entry).unwrap();

    let result = acc.apply(&mut molecule_with_properties);
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
fn test_apply_atom_value_entry(mut basic_molecule: Molecule) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::AtomValueEntry(AtomValueEntry {
        atom_index: 0,
        value: "*".to_string(),
    });
    acc.add_entry(entry).unwrap();
    acc.apply(&mut basic_molecule).unwrap();

    let value = basic_molecule
        .atom(0)
        .unwrap()
        .properties
        .get("molFileValue")
        .unwrap();
    assert_eq!(value, "*");
}

#[rstest]
fn test_apply_atom_value_entry_invalid(mut basic_molecule: Molecule) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::AtomValueEntry(AtomValueEntry {
        atom_index: 5,
        value: "*".to_string(),
    });
    acc.add_entry(entry).unwrap();

    let result = acc.apply(&mut basic_molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Data(DataError::MissingAtomIndex(idx)) => assert_eq!(idx, 5),
        _ => panic!("Expected MissingAtomIndex error"),
    }
}

#[rstest]
fn test_apply_atom_value_conflict(mut molecule_with_properties: Molecule) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::AtomValueEntry(AtomValueEntry {
        atom_index: 0,
        value: "new".to_string(),
    });
    acc.add_entry(entry).unwrap();

    let result = acc.apply(&mut molecule_with_properties);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Validation(ValidationError::InvalidComponent(msg)) => {
            assert!(msg.contains("Atom value conflict"));
            assert!(msg.contains("existing 'existing' vs new 'new'"));
        }
        _ => panic!("Expected InvalidComponent error"),
    }
}

#[rstest]
fn test_apply_charge_entries(mut basic_molecule: Molecule) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::ChargeEntries(vec![ChargeEntry {
        atom_index: 0,
        charge: -1,
    }]);
    acc.add_entry(entry).unwrap();
    acc.apply(&mut basic_molecule).unwrap();

    assert_eq!(basic_molecule.atom(0).unwrap().charge, -1);
}

#[rstest]
fn test_apply_charge_entries_multiple(mut basic_molecule: Molecule) {
    let mut acc = MoleculeProperties::new();
    basic_molecule.add_atom(Atom::new(AtomSymbol::Element(e!(C))));
    basic_molecule.add_atom(Atom::new(AtomSymbol::Element(e!(C))));

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
    acc.add_entry(PropertyEntries::ChargeEntries(entries))
        .unwrap();
    acc.apply(&mut basic_molecule).unwrap();

    assert_eq!(basic_molecule.atom(0).unwrap().charge, -1);
    assert_eq!(basic_molecule.atom(1).unwrap().charge, 0); // unchanged
    assert_eq!(basic_molecule.atom(2).unwrap().charge, 1);
}

#[rstest]
fn test_apply_charge_entries_invalid(mut basic_molecule: Molecule) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::ChargeEntries(vec![ChargeEntry {
        atom_index: 5,
        charge: -1,
    }]);
    acc.add_entry(entry).unwrap();
    let result = acc.apply(&mut basic_molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Data(DataError::MissingAtomIndex(idx)) => assert_eq!(idx, 5),
        _ => panic!("Expected MissingAtomIndex error"),
    }
}

#[rstest]
fn test_apply_charge_entries_overwrite(mut molecule_with_properties: Molecule) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::ChargeEntries(vec![ChargeEntry {
        atom_index: 0,
        charge: -2,
    }]);
    acc.add_entry(entry).unwrap();
    acc.apply(&mut molecule_with_properties).unwrap();

    let atom = molecule_with_properties.atom(0).unwrap();
    assert_eq!(atom.charge, -2);
    assert_eq!(atom.radical, None);
}

#[rstest]
#[case(0, None)]
#[case(1, Some(AtomRadical::Singlet))]
#[case(2, Some(AtomRadical::Doublet))]
#[case(3, Some(AtomRadical::Triplet))]
fn test_apply_radical_entries(
    mut basic_molecule: Molecule,
    #[case] radical_type: u8,
    #[case] expected: Option<AtomRadical>,
) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::RadicalEntries(vec![RadicalEntry {
        atom_index: 0,
        radical_type,
    }]);
    acc.add_entry(entry).unwrap();
    acc.apply(&mut basic_molecule).unwrap();
    assert_eq!(basic_molecule.atom(0).unwrap().radical, expected);
}

#[rstest]
fn test_apply_radical_entries_invalid() {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::RadicalEntries(vec![RadicalEntry {
        atom_index: 0,
        radical_type: 4,
    }]);
    let result = acc.add_entry(entry);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Validation(ValidationError::InvalidComponent(msg)) => {
            assert!(msg.contains("Invalid radical type"));
        }
        _ => panic!("Expected InvalidComponent error"),
    }
}

#[rstest]
fn test_apply_radical_entries_overwrite(mut molecule_with_properties: Molecule) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::RadicalEntries(vec![RadicalEntry {
        atom_index: 0,
        radical_type: 1,
    }]);
    acc.add_entry(entry).unwrap();
    acc.apply(&mut molecule_with_properties).unwrap();

    let atom = molecule_with_properties.atom(0).unwrap();
    assert_eq!(atom.radical, Some(AtomRadical::Singlet));
    assert_eq!(atom.charge, 0);
}

#[rstest]
fn test_apply_isotope_entries(mut basic_molecule: Molecule) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::IsotopeEntries(vec![IsotopeEntry {
        atom_index: 0,
        mass: 14,
    }]);
    acc.add_entry(entry).unwrap();
    acc.apply(&mut basic_molecule).unwrap();

    let atom = basic_molecule.atom(0).unwrap();
    assert_eq!(atom.isotope_mass, Some(14));
}

#[rstest]
fn test_apply_isotope_entries_conflict(mut molecule_with_properties: Molecule) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::IsotopeEntries(vec![IsotopeEntry {
        atom_index: 0,
        mass: 14,
    }]);
    acc.add_entry(entry).unwrap();
    let result = acc.apply(&mut molecule_with_properties);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Validation(ValidationError::InvalidComponent(msg)) => {
            assert!(msg.contains("Isotope conflict: existing"));
        }
        _ => panic!("Expected InvalidComponent error"),
    }
}

#[rstest]
fn test_apply_isotope_entries_invalid(mut basic_molecule: Molecule) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::IsotopeEntries(vec![IsotopeEntry {
        atom_index: 0,
        mass: 40,
    }]);
    acc.add_entry(entry).unwrap();
    let result = acc.apply(&mut basic_molecule);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Data(DataError::InvalidIsotope(msg)) => {
            assert!(msg.contains("Invalid isotope mass number 40 for element C"));
        }
        _ => panic!("Expected InvalidIsotope error"),
    }
}

#[rstest]
fn test_apply_isotope_on_named_isotope(mut basic_molecule: Molecule) {
    basic_molecule.atom_mut(0).unwrap().symbol = AtomSymbol::NamedIsotope(NamedIsotope::D);
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::IsotopeEntries(vec![IsotopeEntry {
        atom_index: 0,
        mass: 3,
    }]);
    acc.add_entry(entry).unwrap();
    acc.apply(&mut basic_molecule).unwrap();

    let atom = basic_molecule.atom(0).unwrap();
    assert_eq!(atom.isotope_mass, Some(3));
}

#[rstest]
#[case(0, None)]
#[case(2, Some(RingBondCount::R2))]
#[case(4, Some(RingBondCount::R4Plus))]
#[case(-1, Some(RingBondCount::NoRingBonds))]
fn test_apply_ring_bond_count_entries(
    mut basic_molecule: Molecule,
    #[case] code: i8,
    #[case] expected: Option<RingBondCount>,
) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::RingBondCountEntries(vec![RingBondCountEntry {
        atom_index: 0,
        ring_bond_count: code,
    }]);
    acc.add_entry(entry).unwrap();
    acc.apply(&mut basic_molecule).unwrap();
    assert_eq!(basic_molecule.atom(0).unwrap().ring_bond_count, expected);
}

#[rstest]
fn test_apply_ring_bond_count_entries_conflict(mut molecule_with_properties: Molecule) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::RingBondCountEntries(vec![RingBondCountEntry {
        atom_index: 0,
        ring_bond_count: 3,
    }]);
    acc.add_entry(entry).unwrap();
    let result = acc.apply(&mut molecule_with_properties);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Validation(ValidationError::InvalidComponent(msg)) => {
            assert!(msg.contains("Ring bond count conflict: existing"));
        }
        _ => panic!("Expected InvalidComponent error"),
    }
}

#[rstest]
#[case(1)]
#[case(5)]
fn test_apply_ring_bond_count_entries_invalid(#[case] code: i8) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::RingBondCountEntries(vec![RingBondCountEntry {
        atom_index: 0,
        ring_bond_count: code,
    }]);
    let result = acc.add_entry(entry);
    assert!(result.is_err());
}

#[rstest]
#[case(0, None)]
#[case(2, Some(SubstitutionCount::S2))]
#[case(4, Some(SubstitutionCount::S4))]
#[case(-1, Some(SubstitutionCount::NoSubstitution))]
fn test_apply_substitution_count_entries(
    mut basic_molecule: Molecule,
    #[case] code: i8,
    #[case] expected: Option<SubstitutionCount>,
) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::SubstitutionCountEntries(vec![SubstitutionCountEntry {
        atom_index: 0,
        substitution_count: code,
    }]);
    acc.add_entry(entry).unwrap();
    acc.apply(&mut basic_molecule).unwrap();
    assert_eq!(basic_molecule.atom(0).unwrap().substitution_count, expected);
}

#[rstest]
fn test_apply_substitution_count_entries_conflict(mut molecule_with_properties: Molecule) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::SubstitutionCountEntries(vec![SubstitutionCountEntry {
        atom_index: 0,
        substitution_count: 3,
    }]);
    acc.add_entry(entry).unwrap();
    let result = acc.apply(&mut molecule_with_properties);
    assert!(result.is_err());

    match result.unwrap_err() {
        Error::Validation(ValidationError::InvalidComponent(msg)) => {
            assert!(msg.contains("Substitution count conflict: existing"));
        }
        _ => panic!("Expected InvalidComponent error"),
    }
}

#[rstest]
#[case(-3)]
#[case(7)]
fn test_apply_substitution_count_entries_invalid(#[case] code: i8) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::SubstitutionCountEntries(vec![SubstitutionCountEntry {
        atom_index: 0,
        substitution_count: code,
    }]);
    let result = acc.add_entry(entry);
    assert!(result.is_err());
}

#[rstest]
#[case(0, None)]
#[case(1, Some(UnsaturatedAtom))]
fn test_apply_unsaturated_atom_entries(
    mut basic_molecule: Molecule,
    #[case] code: u8,
    #[case] expected: Option<UnsaturatedAtom>,
) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::UnsaturatedAtomEntries(vec![UnsaturatedAtomEntry {
        atom_index: 0,
        unsaturated: code,
    }]);
    acc.add_entry(entry).unwrap();
    acc.apply(&mut basic_molecule).unwrap();
    assert_eq!(basic_molecule.atom(0).unwrap().unsaturated, expected);
}

#[rstest]
#[case(2)]
fn test_apply_unsaturated_atom_entries_invalid(#[case] code: u8) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::UnsaturatedAtomEntries(vec![UnsaturatedAtomEntry {
        atom_index: 0,
        unsaturated: code,
    }]);
    let result = acc.add_entry(entry);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_link_atom_entries(mut basic_molecule: Molecule) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::LinkAtomEntries(vec![LinkAtomEntry {
        atom_index: 0,
        repeat_count: 2,
        subs_index1: 1,
        subs_index2: None,
    }]);
    acc.add_entry(entry).unwrap();
    acc.apply(&mut basic_molecule).unwrap();
    assert_eq!(
        basic_molecule.atom(0).unwrap().link_atom,
        Some(LinkAtom {
            repeat_count: 2,
            subs_index1: 1,
            subs_index2: None
        })
    );
}

#[rstest]
fn test_apply_link_atom_entries_conflict() {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::LinkAtomEntries(vec![
        LinkAtomEntry {
            atom_index: 0,
            repeat_count: 2,
            subs_index1: 1,
            subs_index2: None,
        },
        LinkAtomEntry {
            atom_index: 0,
            repeat_count: 3,
            subs_index1: 2,
            subs_index2: Some(3),
        },
    ]);
    let result = acc.add_entry(entry);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_atom_list_entry(mut basic_molecule: Molecule) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::AtomListEntry(AtomListEntry {
        atom_index: 0,
        elements: vec![e!(N), e!(O)],
        exclusion: false,
    });
    acc.add_entry(entry).unwrap();
    acc.apply(&mut basic_molecule).unwrap();
    assert_eq!(
        basic_molecule.atom(0).unwrap().symbol,
        AtomSymbol::AtomList(AtomList {
            elements: vec![e!(N), e!(O)],
            exclusion: false
        })
    );
}

#[rstest]
fn test_apply_atom_list_entry_conflict(mut basic_molecule: Molecule) {
    let mut acc = MoleculeProperties::new();
    basic_molecule.atom_mut(0).unwrap().symbol =
        AtomSymbol::RGroup(crate::io::ctab::rgroup::RGroup::new(Some(1)));
    let entry = PropertyEntries::AtomListEntry(AtomListEntry {
        atom_index: 0,
        elements: vec![e!(N), e!(O)],
        exclusion: false,
    });
    acc.add_entry(entry).unwrap();
    let result = acc.apply(&mut basic_molecule);
    assert!(result.is_err());
}

#[rstest]
#[case(1, Some(AttachmentPointType::First))]
#[case(2, Some(AttachmentPointType::Second))]
#[case(3, Some(AttachmentPointType::Both))]
#[case(0, None)]
fn test_apply_attachment_point_entries(
    mut basic_molecule: Molecule,
    #[case] code: u8,
    #[case] expected: Option<AttachmentPointType>,
) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::AttachmentPointEntries(vec![AttachmentPointEntry {
        atom_index: 0,
        attachment_type: code,
    }]);
    acc.add_entry(entry).unwrap();
    acc.apply(&mut basic_molecule).unwrap();
    assert_eq!(basic_molecule.atom(0).unwrap().attachment_point, expected);
}

#[rstest]
fn test_apply_attachment_point_entries_conflict() {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::AttachmentPointEntries(vec![
        AttachmentPointEntry {
            atom_index: 0,
            attachment_type: 1,
        },
        AttachmentPointEntry {
            atom_index: 0,
            attachment_type: 2,
        },
    ]);
    let result = acc.add_entry(entry);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_attachment_point_entries_invalid() {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::AttachmentPointEntries(vec![AttachmentPointEntry {
        atom_index: 0,
        attachment_type: 4,
    }]);
    let result = acc.add_entry(entry);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_atom_attachment_order_entry(mut basic_molecule: Molecule) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::AtomAttachmentOrderEntry(AtomAttachmentOrderEntry {
        atom_index: 0,
        attachments: vec![(1, 2), (2, 1)],
    });
    acc.add_entry(entry).unwrap();
    acc.apply(&mut basic_molecule).unwrap();
    assert_eq!(
        basic_molecule
            .atom(0)
            .unwrap()
            .attachment_order
            .as_ref()
            .unwrap(),
        &vec![(1, 2), (2, 1)]
    );
}

#[rstest]
fn test_apply_atom_attachment_order_entry_conflict(mut basic_molecule: Molecule) {
    let mut acc = MoleculeProperties::new();
    basic_molecule.atom_mut(0).unwrap().attachment_order = Some(vec![(1, 1)]);
    let entry = PropertyEntries::AtomAttachmentOrderEntry(AtomAttachmentOrderEntry {
        atom_index: 0,
        attachments: vec![(2, 2)],
    });
    acc.add_entry(entry).unwrap();
    let result = acc.apply(&mut basic_molecule);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_rgroup_label_entries(mut basic_molecule: Molecule) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::RGroupLabelEntries(vec![RGroupLabelEntry {
        atom_index: 0,
        label: 1,
    }]);
    acc.add_entry(entry).unwrap();
    acc.apply(&mut basic_molecule).unwrap();
    assert_eq!(
        basic_molecule.atom(0).unwrap().symbol,
        AtomSymbol::RGroup(RGroup::new(Some(1)))
    );
}

#[rstest]
fn test_apply_rgroup_label_entries_conflict(mut molecule_with_rgroup: Molecule) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::RGroupLabelEntries(vec![RGroupLabelEntry {
        atom_index: 0,
        label: 2,
    }]);
    acc.add_entry(entry).unwrap();
    let result = acc.apply(&mut molecule_with_rgroup);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_rgroup_logic_entry(mut molecule_with_rgroup: Molecule) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::RGroupLogicEntry(RGroupLogicEntry {
        label: 1,
        dependent_label: Some(2),
        rgroup_or_h: true,
        occurrence: vec![RGroupOccurrence::Exactly(1)],
    });
    acc.add_entry(entry).unwrap();
    acc.apply(&mut molecule_with_rgroup).unwrap();

    let atom = molecule_with_rgroup.atom(0).unwrap();
    if let AtomSymbol::RGroup(rgroup) = &atom.symbol {
        assert_eq!(rgroup.dependent_label, Some(2));
        assert!(rgroup.rgroup_or_h);
        assert_eq!(rgroup.occurrence.len(), 1);
        assert_eq!(rgroup.occurrence[0], RGroupOccurrence::Exactly(1));
    } else {
        panic!("Expected RGroup symbol");
    }
}

#[rstest]
fn test_apply_rgroup_logic_entry_multiple_occurrences(mut molecule_with_rgroup: Molecule) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::RGroupLogicEntry(RGroupLogicEntry {
        label: 1,
        dependent_label: Some(2),
        rgroup_or_h: true,
        occurrence: vec![
            RGroupOccurrence::Exactly(1),
            RGroupOccurrence::GreaterThan(5),
        ],
    });
    acc.add_entry(entry).unwrap();
    acc.apply(&mut molecule_with_rgroup).unwrap();

    let atom = molecule_with_rgroup.atom(0).unwrap();
    if let AtomSymbol::RGroup(rgroup) = &atom.symbol {
        assert_eq!(rgroup.occurrence.len(), 2);
        assert_eq!(rgroup.occurrence[0], RGroupOccurrence::Exactly(1));
        assert_eq!(rgroup.occurrence[1], RGroupOccurrence::GreaterThan(5));
    } else {
        panic!("Expected RGroup symbol");
    }
}

#[rstest]
fn test_apply_sgroup_type_entries(mut basic_molecule: Molecule) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry {
        sgroup_index: 0,
        sgroup_type: SGroupType::Superatom,
    }]);
    acc.add_entry(entry).unwrap();
    acc.apply(&mut basic_molecule).unwrap();

    assert_eq!(basic_molecule.sgroups().count(), 1);
    let sgroup = basic_molecule.sgroup(0).unwrap();
    assert_eq!(sgroup.group_type, SGroupType::Superatom);
}

#[rstest]
fn test_apply_sgroup_type_entries_conflict() {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::SGroupTypeEntries(vec![
        SGroupTypeEntry {
            sgroup_index: 0,
            sgroup_type: SGroupType::Superatom,
        },
        SGroupTypeEntry {
            sgroup_index: 0,
            sgroup_type: SGroupType::Data,
        },
    ]);
    let result = acc.add_entry(entry);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_sgroup_subtype_entries(
    mut basic_molecule: Molecule,
    mut acc_with_copolymer_sgroup: MoleculeProperties,
) {
    let subtype_entry = PropertyEntries::SGroupSubtypeEntries(vec![SGroupSubtypeEntry {
        sgroup_index: 0,
        sgroup_subtype: SGroupSubtype::Alternating,
    }]);
    acc_with_copolymer_sgroup.add_entry(subtype_entry).unwrap();
    acc_with_copolymer_sgroup
        .apply(&mut basic_molecule)
        .unwrap();

    assert_eq!(basic_molecule.sgroups().count(), 1);
    let sgroup = basic_molecule.sgroup(0).unwrap();
    assert_eq!(sgroup.group_subtype, Some(SGroupSubtype::Alternating));
}

#[rstest]
fn test_apply_sgroup_subtype_entries_no_type() {
    let mut acc = MoleculeProperties::new();
    let subtype_entry = PropertyEntries::SGroupSubtypeEntries(vec![SGroupSubtypeEntry {
        sgroup_index: 0,
        sgroup_subtype: SGroupSubtype::Alternating,
    }]);
    let result = acc.add_entry(subtype_entry);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_sgroup_label_entries(
    mut basic_molecule: Molecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
) {
    let label_entry = PropertyEntries::SGroupLabelEntries(vec![SGroupLabelEntry {
        sgroup_index: 0,
        label: 123,
    }]);
    acc_with_superatom_sgroup.add_entry(label_entry).unwrap();
    acc_with_superatom_sgroup
        .apply(&mut basic_molecule)
        .unwrap();

    assert_eq!(basic_molecule.sgroups().count(), 1);
    let sgroup = basic_molecule.sgroup(0).unwrap();
    assert_eq!(sgroup.label, Some(123));
}

#[rstest]
fn test_apply_sgroup_label_entries_no_type() {
    let mut acc = MoleculeProperties::new();
    let label_entry = PropertyEntries::SGroupLabelEntries(vec![SGroupLabelEntry {
        sgroup_index: 0,
        label: 123,
    }]);
    let result = acc.add_entry(label_entry);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_sgroup_connectivity_entries(
    mut basic_molecule: Molecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
) {
    let connectivity_entry =
        PropertyEntries::SGroupConnectivityEntries(vec![SGroupConnectivityEntry {
            sgroup_index: 0,
            connectivity: SGroupConnectivity::HeadToTail,
        }]);
    acc_with_superatom_sgroup
        .add_entry(connectivity_entry)
        .unwrap();
    acc_with_superatom_sgroup
        .apply(&mut basic_molecule)
        .unwrap();

    assert_eq!(basic_molecule.sgroups().count(), 1);
    let sgroup = basic_molecule.sgroup(0).unwrap();
    assert_eq!(sgroup.connectivity, Some(SGroupConnectivity::HeadToTail));
}

#[rstest]
fn test_apply_sgroup_connectivity_entries_no_type() {
    let mut acc = MoleculeProperties::new();
    let connectivity_entry =
        PropertyEntries::SGroupConnectivityEntries(vec![SGroupConnectivityEntry {
            sgroup_index: 0,
            connectivity: SGroupConnectivity::HeadToTail,
        }]);
    let result = acc.add_entry(connectivity_entry);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_sgroup_expansion_entries(
    mut basic_molecule: Molecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
) {
    let expansion_entry =
        PropertyEntries::SGroupExpansionEntries(vec![SGroupExpansionEntry { sgroup_index: 0 }]);
    acc_with_superatom_sgroup
        .add_entry(expansion_entry)
        .unwrap();
    acc_with_superatom_sgroup
        .apply(&mut basic_molecule)
        .unwrap();

    assert_eq!(basic_molecule.sgroups().count(), 1);
    let sgroup = basic_molecule.sgroup(0).unwrap();
    assert!(sgroup.expansion);
}

#[rstest]
fn test_apply_sgroup_expansion_entries_no_type() {
    let mut acc = MoleculeProperties::new();
    let expansion_entry =
        PropertyEntries::SGroupExpansionEntries(vec![SGroupExpansionEntry { sgroup_index: 0 }]);
    let result = acc.add_entry(expansion_entry);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_sgroup_atom_list_entry(
    mut basic_molecule: Molecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
) {
    let atom_list_entry = PropertyEntries::SGroupAtomListEntry(SGroupAtomListEntry {
        sgroup_index: 0,
        atom_indices: vec![0, 1],
    });
    acc_with_superatom_sgroup
        .add_entry(atom_list_entry)
        .unwrap();
    acc_with_superatom_sgroup
        .apply(&mut basic_molecule)
        .unwrap();

    assert_eq!(basic_molecule.sgroups().count(), 1);
    let sgroup = basic_molecule.sgroup(0).unwrap();
    assert_eq!(sgroup.atom_indices, vec![0, 1]);
}

#[rstest]
fn test_apply_sgroup_atom_list_entry_no_type() {
    let mut acc = MoleculeProperties::new();
    let atom_list_entry = PropertyEntries::SGroupAtomListEntry(SGroupAtomListEntry {
        sgroup_index: 0,
        atom_indices: vec![0, 1],
    });
    let result = acc.add_entry(atom_list_entry);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_sgroup_bond_list_entry(
    mut basic_molecule: Molecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
) {
    let bond_list_entry = PropertyEntries::SGroupBondListEntry(SGroupBondListEntry {
        sgroup_index: 0,
        bond_indices: vec![0, 1],
    });
    acc_with_superatom_sgroup
        .add_entry(bond_list_entry)
        .unwrap();
    acc_with_superatom_sgroup
        .apply(&mut basic_molecule)
        .unwrap();

    assert_eq!(basic_molecule.sgroups().count(), 1);
    let sgroup = basic_molecule.sgroup(0).unwrap();
    assert_eq!(sgroup.bond_indices, vec![0, 1]);
}

#[rstest]
fn test_apply_sgroup_bond_list_entry_no_type() {
    let mut acc = MoleculeProperties::new();
    let bond_list_entry = PropertyEntries::SGroupBondListEntry(SGroupBondListEntry {
        sgroup_index: 0,
        bond_indices: vec![0, 1],
    });
    let result = acc.add_entry(bond_list_entry);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_sgroup_parent_atom_entry(
    mut basic_molecule: Molecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
) {
    let parent_atom_entry = PropertyEntries::SGroupParentAtomEntry(SGroupParentAtomEntry {
        sgroup_index: 0,
        atom_indices: vec![0, 1],
    });
    acc_with_superatom_sgroup
        .add_entry(parent_atom_entry)
        .unwrap();
    acc_with_superatom_sgroup
        .apply(&mut basic_molecule)
        .unwrap();

    assert_eq!(basic_molecule.sgroups().count(), 1);
    let sgroup = basic_molecule.sgroup(0).unwrap();
    assert_eq!(sgroup.parent_atom_indices, Some(vec![0, 1]));
}

#[rstest]
fn test_apply_sgroup_parent_atom_entry_no_type() {
    let mut acc = MoleculeProperties::new();
    let parent_atom_entry = PropertyEntries::SGroupParentAtomEntry(SGroupParentAtomEntry {
        sgroup_index: 0,
        atom_indices: vec![0, 1],
    });
    let result = acc.add_entry(parent_atom_entry);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_sgroup_subscript_entry_subscript(
    mut basic_molecule: Molecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
) {
    let subscript_entry = PropertyEntries::SGroupSubscriptEntry(SGroupSubscriptEntry {
        sgroup_index: 0,
        data: SGroupSubscriptData::Subscript("Ph".to_string()),
    });
    acc_with_superatom_sgroup
        .add_entry(subscript_entry)
        .unwrap();
    acc_with_superatom_sgroup
        .apply(&mut basic_molecule)
        .unwrap();

    assert_eq!(basic_molecule.sgroups().count(), 1);
    let sgroup = basic_molecule.sgroup(0).unwrap();
    assert_eq!(sgroup.subscript, Some("Ph".to_string()));
    assert_eq!(sgroup.multiplier, None);
}

#[rstest]
fn test_apply_sgroup_subscript_entry_multiplier(
    mut basic_molecule: Molecule,
    mut acc_with_multiple_sgroup: MoleculeProperties,
) {
    let subscript_entry = PropertyEntries::SGroupSubscriptEntry(SGroupSubscriptEntry {
        sgroup_index: 0,
        data: SGroupSubscriptData::Multiplier(SGroupMultiplier::N),
    });
    acc_with_multiple_sgroup
        .add_entry(subscript_entry)
        .unwrap();
    acc_with_multiple_sgroup
        .apply(&mut basic_molecule)
        .unwrap();

    assert_eq!(basic_molecule.sgroups().count(), 1);
    let sgroup = basic_molecule.sgroup(0).unwrap();
    assert_eq!(sgroup.subscript, None);
    assert_eq!(sgroup.multiplier, Some(SGroupMultiplier::N));
}

#[rstest]
fn test_apply_sgroup_subscript_entry_no_type() {
    let mut acc = MoleculeProperties::new();
    let subscript_entry = PropertyEntries::SGroupSubscriptEntry(SGroupSubscriptEntry {
        sgroup_index: 0,
        data: SGroupSubscriptData::Subscript("n".to_string()),
    });
    let result = acc.add_entry(subscript_entry);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_sgroup_correspondence_entry(
    mut basic_molecule: Molecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
) {
    let correspondence_entry =
        PropertyEntries::SGroupCorrespondenceEntry(SGroupCorrespondenceEntry {
            sgroup_index: 0,
            bond_indices: vec![0, 1],
        });
    acc_with_superatom_sgroup
        .add_entry(correspondence_entry)
        .unwrap();
    acc_with_superatom_sgroup
        .apply(&mut basic_molecule)
        .unwrap();

    assert_eq!(basic_molecule.sgroups().count(), 1);
    let sgroup = basic_molecule.sgroup(0).unwrap();
    assert_eq!(sgroup.correspondence, Some(vec![0, 1]));
}

#[rstest]
fn test_apply_sgroup_correspondence_entry_no_type() {
    let mut acc = MoleculeProperties::new();
    let correspondence_entry =
        PropertyEntries::SGroupCorrespondenceEntry(SGroupCorrespondenceEntry {
            sgroup_index: 0,
            bond_indices: vec![0, 1],
        });
    let result = acc.add_entry(correspondence_entry);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_sgroup_display_info_entry(
    mut basic_molecule: Molecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
) {
    let display_info_entry = PropertyEntries::SGroupDisplayInfoEntry(SGroupDisplayInfoEntry {
        sgroup_index: 0,
        bracket_coords: vec![1.0, 2.0, 3.0, 4.0],
    });
    acc_with_superatom_sgroup
        .add_entry(display_info_entry)
        .unwrap();
    acc_with_superatom_sgroup
        .apply(&mut basic_molecule)
        .unwrap();

    assert_eq!(basic_molecule.sgroups().count(), 1);
    let sgroup = basic_molecule.sgroup(0).unwrap();
    assert_eq!(
        sgroup.bracket_coords,
        Some(SGroupBracketCoords {
            bracket1: (1.0, 2.0),
            bracket2: (3.0, 4.0)
        })
    );
}

#[rstest]
fn test_apply_sgroup_display_info_entry_no_type() {
    let mut acc = MoleculeProperties::new();
    let display_info_entry = PropertyEntries::SGroupDisplayInfoEntry(SGroupDisplayInfoEntry {
        sgroup_index: 0,
        bracket_coords: vec![1.0, 2.0, 3.0, 4.0],
    });
    let result = acc.add_entry(display_info_entry);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_sgroup_connecting_bond_entry(
    mut basic_molecule: Molecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
) {
    let connecting_bond_entry =
        PropertyEntries::SGroupConnectingBondEntry(SGroupConnectingBondEntry {
            sgroup_index: 0,
            bond_index: 0,
            bond_vector: (1.0, 2.0),
        });
    acc_with_superatom_sgroup
        .add_entry(connecting_bond_entry)
        .unwrap();
    acc_with_superatom_sgroup
        .apply(&mut basic_molecule)
        .unwrap();

    assert_eq!(basic_molecule.sgroups().count(), 1);
    let sgroup = basic_molecule.sgroup(0).unwrap();
    assert_eq!(
        sgroup.connecting_bond,
        Some(SGroupConnectingBond {
            bond_index: 0,
            bond_vector: (1.0, 2.0)
        })
    );
}

#[rstest]
fn test_apply_sgroup_connecting_bond_entry_no_type() {
    let mut acc = MoleculeProperties::new();
    let connecting_bond_entry =
        PropertyEntries::SGroupConnectingBondEntry(SGroupConnectingBondEntry {
            sgroup_index: 0,
            bond_index: 0,
            bond_vector: (1.0, 2.0),
        });
    let result = acc.add_entry(connecting_bond_entry);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_sgroup_data_description_entry(
    mut basic_molecule: Molecule,
    mut acc_with_data_sgroup: MoleculeProperties,
) {
    acc_with_data_sgroup.apply(&mut basic_molecule).unwrap();

    assert_eq!(basic_molecule.sgroups().count(), 1);
    let sgroup = basic_molecule.sgroup(0).unwrap();
    let data = sgroup.data.get("test").unwrap();
    assert_eq!(data.field_type, SGroupDataType::Text);
    assert_eq!(data.field_units, None);
    assert_eq!(data.query_identifier, None);
    assert_eq!(data.data_query_operator, None);
}

#[rstest]
fn test_apply_sgroup_data_description_entry_no_type() {
    let mut acc = MoleculeProperties::new();
    let data_description_entry =
        PropertyEntries::SGroupDataDescriptionEntry(SGroupDataDescriptionEntry {
            sgroup_index: 0,
            field_name: "test".to_string(),
            field_type: SGroupDataType::Text,
            field_units: Some("unit".to_string()),
            query_identifier: Some("Q".to_string()),
            data_query_operator: Some("=".to_string()),
        });
    let result = acc.add_entry(data_description_entry);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_sgroup_data_display_entry(
    mut basic_molecule: Molecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
) {
    let display_entry = PropertyEntries::SGroupDataDisplayEntry(SGroupDataDisplayEntry {
        sgroup_index: 0,
        coords: (10.0, 20.0),
        display_type: SGroupDataDisplayType::Detached,
        display_placement: SGroupDataDisplayPlacement::Relative,
        display_units: SGroupDataDisplayUnits::DisplayUnits,
        display_chars: SGroupDataDisplayChars::Number(5),
        display_tag: Some(5),
        display_position: 3,
    });
    acc_with_superatom_sgroup.add_entry(display_entry).unwrap();
    acc_with_superatom_sgroup
        .apply(&mut basic_molecule)
        .unwrap();

    assert_eq!(basic_molecule.sgroups().count(), 1);
    let sgroup = basic_molecule.sgroup(0).unwrap();
    let display = sgroup.display.as_ref().unwrap();
    assert_eq!(display.coords, (10.0, 20.0));
    assert_eq!(display.display_type, SGroupDataDisplayType::Detached);
    assert_eq!(
        display.display_placement,
        SGroupDataDisplayPlacement::Relative
    );
    assert_eq!(display.display_units, SGroupDataDisplayUnits::DisplayUnits);
    assert_eq!(display.display_chars, SGroupDataDisplayChars::Number(5));
}

#[rstest]
fn test_apply_sgroup_data_display_entry_no_type() {
    let mut acc = MoleculeProperties::new();
    let display_entry = PropertyEntries::SGroupDataDisplayEntry(SGroupDataDisplayEntry {
        sgroup_index: 0,
        coords: (10.0, 20.0),
        display_type: SGroupDataDisplayType::Detached,
        display_placement: SGroupDataDisplayPlacement::Relative,
        display_units: SGroupDataDisplayUnits::DisplayUnits,
        display_chars: SGroupDataDisplayChars::All,
        display_tag: None,
        display_position: 0,
    });
    let result = acc.add_entry(display_entry);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_sgroup_data_entry(
    mut basic_molecule: Molecule,
    mut acc_with_data_sgroup: MoleculeProperties,
) {
    let continuation_entry = PropertyEntries::SGroupDataEntry(SGroupDataEntry::Continuation {
        sgroup_index: 0,
        data_content: "content".to_string(),
    });
    acc_with_data_sgroup.add_entry(continuation_entry).unwrap();

    let data_entry =
        PropertyEntries::SGroupDataEntry(SGroupDataEntry::EndBlank { sgroup_index: 0 });
    acc_with_data_sgroup.add_entry(data_entry).unwrap();
    acc_with_data_sgroup.apply(&mut basic_molecule).unwrap();

    assert_eq!(basic_molecule.sgroups().count(), 1);
    let sgroup = basic_molecule.sgroup(0).unwrap();
    let data = sgroup.data.get("test").unwrap();
    assert_eq!(data.data_content, Some(vec!["content".to_string()]));
}

#[rstest]
fn test_apply_sgroup_data_entry_no_type() {
    let mut acc = MoleculeProperties::new();
    let data_entry =
        PropertyEntries::SGroupDataEntry(SGroupDataEntry::EndBlank { sgroup_index: 0 });
    let result = acc.add_entry(data_entry);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_sgroup_data_entry_no_description() {
    let mut acc = MoleculeProperties::new();
    let type_entry = PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry {
        sgroup_index: 0,
        sgroup_type: SGroupType::Data,
    }]);
    acc.add_entry(type_entry).unwrap();
    let data_entry =
        PropertyEntries::SGroupDataEntry(SGroupDataEntry::EndBlank { sgroup_index: 0 });
    let result = acc.add_entry(data_entry);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_sgroup_data_entry_auto_finalization(
    mut basic_molecule: Molecule,
    mut acc_with_data_sgroup: MoleculeProperties,
) {
    let continuation_entry = PropertyEntries::SGroupDataEntry(SGroupDataEntry::Continuation {
        sgroup_index: 0,
        data_content: "incomplete".to_string(),
    });
    acc_with_data_sgroup.add_entry(continuation_entry).unwrap();

    // Apply without SED - should auto-finalize
    acc_with_data_sgroup.apply(&mut basic_molecule).unwrap();

    assert_eq!(basic_molecule.sgroups().count(), 1);
    let sgroup = basic_molecule.sgroup(0).unwrap();
    let data = sgroup.data.get("test").unwrap();
    assert_eq!(data.data_content, Some(vec!["incomplete".to_string()]));
}

#[rstest]
fn test_apply_sgroup_data_entry_multi_continuation() {
    let mut acc = MoleculeProperties::new();
    let type_entry = PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry {
        sgroup_index: 0,
        sgroup_type: SGroupType::Data,
    }]);
    acc.add_entry(type_entry).unwrap();

    let data_description_entry =
        PropertyEntries::SGroupDataDescriptionEntry(SGroupDataDescriptionEntry {
            sgroup_index: 0,
            field_name: "test".to_string(),
            field_type: SGroupDataType::Text,
            field_units: None,
            query_identifier: None,
            data_query_operator: None,
        });
    acc.add_entry(data_description_entry).unwrap();

    let continuation1 = PropertyEntries::SGroupDataEntry(SGroupDataEntry::Continuation {
        sgroup_index: 0,
        data_content: "part1".to_string(),
    });
    acc.add_entry(continuation1).unwrap();

    let continuation2 = PropertyEntries::SGroupDataEntry(SGroupDataEntry::Continuation {
        sgroup_index: 0,
        data_content: "part2".to_string(),
    });
    acc.add_entry(continuation2).unwrap();

    let end_entry = PropertyEntries::SGroupDataEntry(SGroupDataEntry::EndWithData {
        sgroup_index: 0,
        data_content: "final".to_string(),
    });
    acc.add_entry(end_entry).unwrap();

    let mut molecule = Molecule::new();
    acc.apply(&mut molecule).unwrap();

    let sgroup = molecule.sgroup(0).unwrap();
    let data = sgroup.data.get("test").unwrap();
    assert_eq!(data.data_content, Some(vec!["part1part2final".to_string()]));
}

#[rstest]
fn test_apply_sgroup_hierarchy_entries(
    mut basic_molecule: Molecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
) {
    // Add a second SGroup as parent
    let parent_type_entry = PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry {
        sgroup_index: 1,
        sgroup_type: SGroupType::Superatom,
    }]);
    acc_with_superatom_sgroup
        .add_entry(parent_type_entry)
        .unwrap();

    let hierarchy_entry = PropertyEntries::SGroupHierarchyEntries(vec![SGroupHierarchyEntry {
        sgroup_index: 0,
        parent_sgroup_index: 1,
    }]);
    acc_with_superatom_sgroup
        .add_entry(hierarchy_entry)
        .unwrap();
    acc_with_superatom_sgroup
        .apply(&mut basic_molecule)
        .unwrap();

    assert_eq!(basic_molecule.sgroups().count(), 2);
    let child_sgroup = basic_molecule.sgroup(0).unwrap();
    assert_eq!(child_sgroup.hierarchy_parent, Some(1));
}

#[rstest]
fn test_apply_sgroup_hierarchy_entries_no_type() {
    let mut acc = MoleculeProperties::new();
    let hierarchy_entry = PropertyEntries::SGroupHierarchyEntries(vec![SGroupHierarchyEntry {
        sgroup_index: 0,
        parent_sgroup_index: 1,
    }]);
    let result = acc.add_entry(hierarchy_entry);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_sgroup_component_entries(
    mut basic_molecule: Molecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
) {
    let component_entry = PropertyEntries::SGroupComponentEntries(vec![SGroupComponentEntry {
        sgroup_index: 0,
        component_number: 42,
    }]);
    acc_with_superatom_sgroup
        .add_entry(component_entry)
        .unwrap();
    acc_with_superatom_sgroup
        .apply(&mut basic_molecule)
        .unwrap();

    assert_eq!(basic_molecule.sgroups().count(), 1);
    let sgroup = basic_molecule.sgroup(0).unwrap();
    assert_eq!(sgroup.component_number, Some(42));
}

#[rstest]
fn test_apply_sgroup_component_entries_no_type() {
    let mut acc = MoleculeProperties::new();
    let component_entry = PropertyEntries::SGroupComponentEntries(vec![SGroupComponentEntry {
        sgroup_index: 0,
        component_number: 42,
    }]);
    let result = acc.add_entry(component_entry);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_zero_order_bond_entries(mut molecule_with_bond: Molecule) {
    let mut acc = MoleculeProperties::new();

    let zero_order_entry = PropertyEntries::ZeroBondOrderEntries(vec![ZeroBondOrderEntry {
        bond_index: 0,
        bond_order: 0,
    }]);
    acc.add_entry(zero_order_entry).unwrap();
    acc.apply(&mut molecule_with_bond).unwrap();

    let bond = molecule_with_bond.bond(0).unwrap();
    assert_eq!(bond.bond_type, BondType::Zero);
}

#[rstest]
fn test_apply_zero_order_bond_entries_invalid(basic_molecule: Molecule) {
    let mut acc = MoleculeProperties::new();
    let mut molecule = basic_molecule;

    let zero_order_entry = PropertyEntries::ZeroBondOrderEntries(vec![ZeroBondOrderEntry {
        bond_index: 5,
        bond_order: 0,
    }]);
    acc.add_entry(zero_order_entry).unwrap();
    let result = acc.apply(&mut molecule);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_zero_atom_charge_entries(mut basic_molecule: Molecule) {
    let mut acc = MoleculeProperties::new();

    let zero_atom_charge_entry =
        PropertyEntries::ZeroAtomChargeEntries(vec![ZeroAtomChargeEntry {
            atom_index: 0,
            charge: -1,
        }]);
    acc.add_entry(zero_atom_charge_entry).unwrap();
    acc.apply(&mut basic_molecule).unwrap();

    let atom = basic_molecule.atom(0).unwrap();
    assert_eq!(atom.charge, -1);
}

#[rstest]
fn test_apply_zero_atom_charge_entries_invalid(basic_molecule: Molecule) {
    let mut acc = MoleculeProperties::new();
    let mut molecule = basic_molecule;

    let zero_atom_charge_entry =
        PropertyEntries::ZeroAtomChargeEntries(vec![ZeroAtomChargeEntry {
            atom_index: 5,
            charge: -1,
        }]);
    acc.add_entry(zero_atom_charge_entry).unwrap();
    let result = acc.apply(&mut molecule);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_atom_hydrogen_count_entries(mut basic_molecule: Molecule) {
    let mut acc = MoleculeProperties::new();

    let atom_hydrogen_count_entry =
        PropertyEntries::AtomHydrogenCountEntries(vec![AtomHydrogenCountEntry {
            atom_index: 0,
            hydrogen_count: 1,
        }]);
    acc.add_entry(atom_hydrogen_count_entry).unwrap();
    acc.apply(&mut basic_molecule).unwrap();

    let atom = basic_molecule.atom(0).unwrap();
    assert_eq!(atom.hydrogen_count, Some(1));
}

#[rstest]
fn test_apply_atom_hydrogen_count_entries_invalid(basic_molecule: Molecule) {
    let mut acc = MoleculeProperties::new();
    let mut molecule = basic_molecule;

    let atom_hydrogen_count_entry =
        PropertyEntries::AtomHydrogenCountEntries(vec![AtomHydrogenCountEntry {
            atom_index: 5,
            hydrogen_count: 1,
        }]);
    acc.add_entry(atom_hydrogen_count_entry).unwrap();
    let result = acc.apply(&mut molecule);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_atom_alias_standard(mut basic_molecule_standard: MoleculeStandard) {
    let mut acc = MoleculeProperties::new();

    let alias_entry = PropertyEntries::AtomAliasEntry(AtomAliasEntry {
        atom_index: 0,
        alias: "CF3".to_string(),
    });
    acc.add_entry(alias_entry).unwrap();
    acc.apply_standard(&mut basic_molecule_standard).unwrap();

    let atom = basic_molecule_standard.atom(0).unwrap();
    assert_eq!(
        atom.properties.get("molFileAlias"),
        Some(&"CF3".to_string())
    );
}

#[rstest]
fn test_apply_atom_alias_standard_invalid(basic_molecule_standard: MoleculeStandard) {
    let mut acc = MoleculeProperties::new();
    let mut molecule = basic_molecule_standard;

    let alias_entry = PropertyEntries::AtomAliasEntry(AtomAliasEntry {
        atom_index: 5,
        alias: "CF3".to_string(),
    });
    acc.add_entry(alias_entry).unwrap();
    let result = acc.apply_standard(&mut molecule);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_atom_value_standard(mut basic_molecule_standard: MoleculeStandard) {
    let mut acc = MoleculeProperties::new();

    let value_entry = PropertyEntries::AtomValueEntry(AtomValueEntry {
        atom_index: 0,
        value: "*".to_string(),
    });
    acc.add_entry(value_entry).unwrap();
    acc.apply_standard(&mut basic_molecule_standard).unwrap();

    let atom = basic_molecule_standard.atom(0).unwrap();
    assert_eq!(atom.properties.get("molFileValue"), Some(&"*".to_string()));
}

#[rstest]
fn test_apply_atom_value_standard_invalid(basic_molecule_standard: MoleculeStandard) {
    let mut acc = MoleculeProperties::new();
    let mut molecule = basic_molecule_standard;

    let value_entry = PropertyEntries::AtomValueEntry(AtomValueEntry {
        atom_index: 5,
        value: "*".to_string(),
    });
    acc.add_entry(value_entry).unwrap();
    let result = acc.apply_standard(&mut molecule);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_charge_standard(mut basic_molecule_standard: MoleculeStandard) {
    let mut acc = MoleculeProperties::new();

    let charge_entry = PropertyEntries::ChargeEntries(vec![ChargeEntry {
        atom_index: 0,
        charge: -1,
    }]);
    acc.add_entry(charge_entry).unwrap();
    acc.apply_standard(&mut basic_molecule_standard).unwrap();

    let atom = basic_molecule_standard.atom(0).unwrap();
    assert_eq!(atom.charge, -1);
    assert_eq!(atom.radical, None);
}

#[rstest]
fn test_apply_charge_standard_invalid(basic_molecule_standard: MoleculeStandard) {
    let mut acc = MoleculeProperties::new();
    let mut molecule = basic_molecule_standard;

    let charge_entry = PropertyEntries::ChargeEntries(vec![ChargeEntry {
        atom_index: 5,
        charge: -1,
    }]);
    acc.add_entry(charge_entry).unwrap();
    let result = acc.apply_standard(&mut molecule);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_radical_standard(mut basic_molecule_standard: MoleculeStandard) {
    let mut acc = MoleculeProperties::new();

    let radical_entry = PropertyEntries::RadicalEntries(vec![RadicalEntry {
        atom_index: 0,
        radical_type: 2, // Doublet
    }]);
    acc.add_entry(radical_entry).unwrap();
    acc.apply_standard(&mut basic_molecule_standard).unwrap();

    let atom = basic_molecule_standard.atom(0).unwrap();
    assert_eq!(atom.radical, Some(AtomRadical::Doublet));
    assert_eq!(atom.charge, 0);
}

#[rstest]
fn test_apply_radical_standard_invalid(basic_molecule_standard: MoleculeStandard) {
    let mut acc = MoleculeProperties::new();
    let mut molecule = basic_molecule_standard;

    let radical_entry = PropertyEntries::RadicalEntries(vec![RadicalEntry {
        atom_index: 5,
        radical_type: 2,
    }]);
    acc.add_entry(radical_entry).unwrap();
    let result = acc.apply_standard(&mut molecule);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_isotope_standard(mut basic_molecule_standard: MoleculeStandard) {
    let mut acc = MoleculeProperties::new();

    let isotope_entry = PropertyEntries::IsotopeEntries(vec![IsotopeEntry {
        atom_index: 0,
        mass: 13, // Carbon-13
    }]);
    acc.add_entry(isotope_entry).unwrap();
    acc.apply_standard(&mut basic_molecule_standard).unwrap();

    let atom = basic_molecule_standard.atom(0).unwrap();
    assert_eq!(atom.isotope_mass, Some(13));
}

#[rstest]
fn test_apply_isotope_standard_invalid(basic_molecule_standard: MoleculeStandard) {
    let mut acc = MoleculeProperties::new();
    let mut molecule = basic_molecule_standard;

    let isotope_entry = PropertyEntries::IsotopeEntries(vec![IsotopeEntry {
        atom_index: 5,
        mass: 13,
    }]);
    acc.add_entry(isotope_entry).unwrap();
    let result = acc.apply_standard(&mut molecule);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_zero_order_bond_standard(mut molecule_with_bond_standard: MoleculeStandard) {
    let mut acc = MoleculeProperties::new();

    let zero_order_entry = PropertyEntries::ZeroBondOrderEntries(vec![ZeroBondOrderEntry {
        bond_index: 0,
        bond_order: 0,
    }]);
    acc.add_entry(zero_order_entry).unwrap();
    acc.apply_standard(&mut molecule_with_bond_standard)
        .unwrap();

    let bond = molecule_with_bond_standard.bond(0).unwrap();
    assert_eq!(bond.bond_type, BondType::Zero);
}

#[rstest]
fn test_apply_zero_order_bond_standard_invalid(basic_molecule_standard: MoleculeStandard) {
    let mut acc = MoleculeProperties::new();
    let mut molecule = basic_molecule_standard;

    let zero_order_entry = PropertyEntries::ZeroBondOrderEntries(vec![ZeroBondOrderEntry {
        bond_index: 5,
        bond_order: 0,
    }]);
    acc.add_entry(zero_order_entry).unwrap();
    let result = acc.apply_standard(&mut molecule);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_zero_atom_charge_entries_standard(mut basic_molecule_standard: MoleculeStandard) {
    let mut acc = MoleculeProperties::new();

    let zero_atom_charge_entry =
        PropertyEntries::ZeroAtomChargeEntries(vec![ZeroAtomChargeEntry {
            atom_index: 0,
            charge: -1,
        }]);
    acc.add_entry(zero_atom_charge_entry).unwrap();
    acc.apply_standard(&mut basic_molecule_standard).unwrap();

    let atom = basic_molecule_standard.atom(0).unwrap();
    assert_eq!(atom.charge, -1);
}

#[rstest]
fn test_apply_zero_atom_charge_entries_standard_invalid(basic_molecule_standard: MoleculeStandard) {
    let mut acc = MoleculeProperties::new();
    let mut molecule = basic_molecule_standard;

    let zero_atom_charge_entry =
        PropertyEntries::ZeroAtomChargeEntries(vec![ZeroAtomChargeEntry {
            atom_index: 5,
            charge: -1,
        }]);
    acc.add_entry(zero_atom_charge_entry).unwrap();
    let result = acc.apply_standard(&mut molecule);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_atom_hydrogen_count_entries_standard(mut basic_molecule_standard: MoleculeStandard) {
    let mut acc = MoleculeProperties::new();

    let atom_hydrogen_count_entry =
        PropertyEntries::AtomHydrogenCountEntries(vec![AtomHydrogenCountEntry {
            atom_index: 0,
            hydrogen_count: 1,
        }]);

    acc.add_entry(atom_hydrogen_count_entry).unwrap();
    acc.apply_standard(&mut basic_molecule_standard).unwrap();

    let atom = basic_molecule_standard.atom(0).unwrap();
    assert_eq!(atom.hydrogen_count, Some(1));
}

#[rstest]
fn test_apply_atom_hydrogen_count_entries_standard_invalid(
    basic_molecule_standard: MoleculeStandard,
) {
    let mut acc = MoleculeProperties::new();
    let mut molecule = basic_molecule_standard;

    let atom_hydrogen_count_entry =
        PropertyEntries::AtomHydrogenCountEntries(vec![AtomHydrogenCountEntry {
            atom_index: 5,
            hydrogen_count: 1,
        }]);
    acc.add_entry(atom_hydrogen_count_entry).unwrap();
    let result = acc.apply_standard(&mut molecule);
    assert!(result.is_err());
}
