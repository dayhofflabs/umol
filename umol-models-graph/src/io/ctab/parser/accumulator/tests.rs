//! Unit tests for accumulator applying properties to ExtendedMolecule

use pretty_assertions::assert_eq;
use rstest::*;
use umol_data::{e, Element, NamedIsotope};

use super::*;
use crate::io::ctab::config::CtabParseFlags;
use crate::io::ctab::parser::properties::{
    AtomAliasEntry, AtomAttachmentOrderEntry, AtomHydrogenCountEntry, AtomListEntry,
    AtomValueEntry, AttachmentPointEntry, ChargeEntry, IsotopeEntry, LinkAtomEntry,
    PropertyEntries, RGroupLabelEntry, RGroupLogicEntry, RadicalEntry, RingBondCountEntry,
    SGroupAtomListEntry, SGroupBondListEntry, SGroupComponentEntry, SGroupConnectingBondEntry,
    SGroupConnectivityEntry, SGroupCorrespondenceEntry, SGroupDataDescriptionEntry,
    SGroupDataDisplayEntry, SGroupDataEntry, SGroupDisplayInfoEntry, SGroupExpansionEntry,
    SGroupHierarchyEntry, SGroupLabelEntry, SGroupParentAtomEntry, SGroupSubscriptEntry,
    SGroupSubtypeEntry, SGroupTypeEntry, SubstitutionCountEntry, UnsaturatedAtomEntry,
    ZeroAtomChargeEntry, ZeroBondOrderEntry,
};
use crate::simple_ir::{
    AtomList, AtomRadical, AtomSymbol, AttachmentPointType, BondOrder, ExtendedAtom, ExtendedBond,
    ExtendedMolecule, LinkAtom, RGroup, RGroupOccurrence, RingBondCount, SGroupBracketCoords,
    SGroupConnectingBond, SGroupConnectivity, SGroupDataDisplayChars, SGroupDataDisplayPlacement,
    SGroupDataDisplayType, SGroupDataDisplayUnits, SGroupDataType, SGroupMultiplier,
    SGroupMultiplierTerm, SGroupSubtype, SGroupType, SubstitutionCount, UnsaturatedAtom,
};

#[fixture]
fn single_atom() -> ExtendedMolecule {
    let mut molecule = ExtendedMolecule::default();
    molecule.atoms.push(ExtendedAtom {
        symbol: AtomSymbol::Element(e!(C)),
        ..Default::default()
    });
    molecule
}

#[fixture]
fn with_properties(mut single_atom: ExtendedMolecule) -> ExtendedMolecule {
    if let Some(atom) = single_atom.atoms.get_mut(0) {
        atom.properties
            .insert("molFileAlias".to_string(), "existing".to_string());
        atom.properties
            .insert("molFileValue".to_string(), "existing".to_string());
        atom.charge = Some(1);
        atom.radical = Some(AtomRadical::Doublet);
        atom.isotope = Some(13);
        atom.ring_bond_count = Some(RingBondCount::R2);
        atom.substitution_count = Some(SubstitutionCount::S2);
        atom.unsaturated = Some(UnsaturatedAtom);
    }
    single_atom
}

#[fixture]
fn with_rgroup(mut single_atom: ExtendedMolecule) -> ExtendedMolecule {
    single_atom.atoms[0].symbol = AtomSymbol::RGroup(RGroup::new(Some(1)));
    single_atom
}

#[fixture]
fn with_unlabeled_rgroup(mut single_atom: ExtendedMolecule) -> ExtendedMolecule {
    single_atom.atoms[0].symbol = AtomSymbol::RGroup(RGroup::new(None));
    single_atom
}

#[fixture]
fn with_bond() -> ExtendedMolecule {
    let mut molecule = ExtendedMolecule::default();
    molecule.atoms.push(ExtendedAtom {
        symbol: AtomSymbol::Element(e!(C)),
        ..Default::default()
    });
    molecule.atoms.push(ExtendedAtom {
        symbol: AtomSymbol::Element(e!(C)),
        ..Default::default()
    });
    molecule.bonds.push(ExtendedBond {
        start_atom: 0,
        end_atom: 1,
        order: BondOrder::Single,
        ..Default::default()
    });
    molecule
}

#[fixture]
fn acc_with_superatom_sgroup(flags_lenient: CtabParseFlags) -> MoleculeProperties {
    let mut acc = MoleculeProperties::new();
    let type_entry = PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry {
        sgroup_index: 0,
        sgroup_type: SGroupType::Superatom,
    }]);
    acc.add_entry(type_entry, flags_lenient).unwrap();
    acc
}

#[fixture]
fn acc_with_data_sgroup(flags_lenient: CtabParseFlags) -> MoleculeProperties {
    let mut acc = MoleculeProperties::new();
    let type_entry = PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry {
        sgroup_index: 0,
        sgroup_type: SGroupType::Data,
    }]);
    acc.add_entry(type_entry, flags_lenient).unwrap();
    let data_description_entry =
        PropertyEntries::SGroupDataDescriptionEntry(SGroupDataDescriptionEntry {
            sgroup_index: 0,
            field_name: "test".to_string(),
            field_type: SGroupDataType::Text,
            field_units: None,
            query_identifier: None,
            data_query_operator: None,
        });
    acc.add_entry(data_description_entry, flags_lenient).unwrap();
    acc
}

#[fixture]
fn acc_with_multiple_sgroup(flags_lenient: CtabParseFlags) -> MoleculeProperties {
    let mut acc = MoleculeProperties::new();
    let type_entry = PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry {
        sgroup_index: 0,
        sgroup_type: SGroupType::MultipleGroup,
    }]);
    acc.add_entry(type_entry, flags_lenient).unwrap();
    acc
}

#[fixture]
fn acc_with_copolymer_sgroup(flags_lenient: CtabParseFlags) -> MoleculeProperties {
    let mut acc = MoleculeProperties::new();
    let type_entry = PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry {
        sgroup_index: 0,
        sgroup_type: SGroupType::Copolymer,
    }]);
    acc.add_entry(type_entry, flags_lenient).unwrap();
    acc
}

#[fixture]
fn flags_lenient() -> CtabParseFlags {
    CtabParseFlags::LENIENT
}

#[fixture]
fn flags_strict() -> CtabParseFlags {
    CtabParseFlags::STRICT
}

#[fixture]
fn flags_extended() -> CtabParseFlags {
    CtabParseFlags::EXTENDED
}

#[rstest]
fn test_apply_atom_alias(mut single_atom: ExtendedMolecule, flags_lenient: CtabParseFlags) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::AtomAliasEntry(AtomAliasEntry {
        atom_index: 0,
        alias: "CF3".to_string(),
    });
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    let atom = &single_atom.atoms[0];
    assert_eq!(
        atom.properties.get("molFileAlias"),
        Some(&"CF3".to_string())
    );
}

#[rstest]
fn test_apply_atom_alias_invalid_index(
    mut single_atom: ExtendedMolecule,
    flags_strict: CtabParseFlags,
) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::AtomAliasEntry(AtomAliasEntry {
        atom_index: 5,
        alias: "CF3".to_string(),
    });
    acc.add_entry(entry, flags_strict).unwrap();
    let result = acc.update_extended_molecule(&mut single_atom, flags_strict);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_atom_value(mut single_atom: ExtendedMolecule, flags_lenient: CtabParseFlags) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::AtomValueEntry(AtomValueEntry {
        atom_index: 0,
        value: "*".to_string(),
    });
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    let atom = &single_atom.atoms[0];
    assert_eq!(atom.properties.get("molFileValue"), Some(&"*".to_string()));
}

#[rstest]
fn test_apply_atom_value_invalid_index(
    mut single_atom: ExtendedMolecule,
    flags_strict: CtabParseFlags,
) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::AtomValueEntry(AtomValueEntry {
        atom_index: 5,
        value: "*".to_string(),
    });
    acc.add_entry(entry, flags_strict).unwrap();
    let result = acc.update_extended_molecule(&mut single_atom, flags_strict);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_charge(mut single_atom: ExtendedMolecule, flags_lenient: CtabParseFlags) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::ChargeEntries(vec![ChargeEntry {
        atom_index: 0,
        charge: -1,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    let atom = &single_atom.atoms[0];
    assert_eq!(atom.charge, Some(-1));
    assert_eq!(atom.radical, None);
}

#[rstest]
fn test_apply_charge_multiple(mut single_atom: ExtendedMolecule, flags_lenient: CtabParseFlags) {
    let mut acc = MoleculeProperties::new();
    single_atom.atoms.push(ExtendedAtom {
        symbol: AtomSymbol::Element(e!(C)),
        ..Default::default()
    });
    single_atom.atoms.push(ExtendedAtom {
        symbol: AtomSymbol::Element(e!(C)),
        ..Default::default()
    });

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
    acc.add_entry(PropertyEntries::ChargeEntries(entries), flags_lenient)
        .unwrap();
    acc.update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    assert_eq!(single_atom.atoms[0].charge, Some(-1));
    assert_eq!(single_atom.atoms[1].charge, None); // unchanged
    assert_eq!(single_atom.atoms[2].charge, Some(1));
}

#[rstest]
fn test_apply_charge_invalid_index(
    mut single_atom: ExtendedMolecule,
    flags_strict: CtabParseFlags,
) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::ChargeEntries(vec![ChargeEntry {
        atom_index: 5,
        charge: -1,
    }]);
    acc.add_entry(entry, flags_strict).unwrap();
    let result = acc.update_extended_molecule(&mut single_atom, flags_strict);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_charge_overwrite(mut with_properties: ExtendedMolecule, flags_lenient: CtabParseFlags) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::ChargeEntries(vec![ChargeEntry {
        atom_index: 0,
        charge: -2,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut with_properties, flags_lenient)
        .unwrap();

    let atom = &with_properties.atoms[0];
    assert_eq!(atom.charge, Some(-2));
    assert_eq!(atom.radical, None);
}

#[rstest]
#[case(0, None)]
#[case(1, Some(AtomRadical::Singlet))]
#[case(2, Some(AtomRadical::Doublet))]
#[case(3, Some(AtomRadical::Triplet))]
fn test_apply_radical(
    mut single_atom: ExtendedMolecule,
    #[case] radical_type: u8,
    #[case] expected: Option<AtomRadical>,
    flags_lenient: CtabParseFlags,
) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::RadicalEntries(vec![RadicalEntry {
        atom_index: 0,
        radical_type,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();
    assert_eq!(single_atom.atoms[0].radical, expected);
}

#[rstest]
fn test_apply_radical_invalid_code(flags_strict: CtabParseFlags) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::RadicalEntries(vec![RadicalEntry {
        atom_index: 0,
        radical_type: 4,
    }]);
    let result = acc.add_entry(entry, flags_strict);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_radical_overwrite(
    mut with_properties: ExtendedMolecule,
    flags_lenient: CtabParseFlags,
) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::RadicalEntries(vec![RadicalEntry {
        atom_index: 0,
        radical_type: 1,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut with_properties, flags_lenient)
        .unwrap();

    let atom = &with_properties.atoms[0];
    assert_eq!(atom.radical, Some(AtomRadical::Singlet));
    assert_eq!(atom.charge, None);
}

#[rstest]
fn test_apply_isotope(mut single_atom: ExtendedMolecule, flags_lenient: CtabParseFlags) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::IsotopeEntries(vec![IsotopeEntry {
        atom_index: 0,
        mass: 14,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    let atom = &single_atom.atoms[0];
    assert_eq!(atom.isotope, Some(14));
}

#[rstest]
fn test_apply_isotope_extended(mut single_atom: ExtendedMolecule, flags_extended: CtabParseFlags) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::IsotopeEntries(vec![IsotopeEntry {
        atom_index: 0,
        mass: 40,
    }]);
    acc.add_entry(entry, flags_extended).unwrap();
    acc.update_extended_molecule(&mut single_atom, flags_extended)
        .unwrap();

    let atom = &single_atom.atoms[0];
    assert_eq!(atom.isotope, Some(40));
}

#[rstest]
fn test_apply_isotope_named_isotope(
    mut single_atom: ExtendedMolecule,
    flags_lenient: CtabParseFlags,
) {
    single_atom.atoms[0].symbol = AtomSymbol::NamedIsotope(NamedIsotope::D);
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::IsotopeEntries(vec![IsotopeEntry {
        atom_index: 0,
        mass: 3,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    let atom = &single_atom.atoms[0];
    assert_eq!(atom.isotope, Some(3));
}

#[rstest]
#[case(0, None)]
#[case(2, Some(RingBondCount::R2))]
#[case(4, Some(RingBondCount::R4Plus))]
#[case(-1, Some(RingBondCount::NoRingBonds))]
fn test_apply_ring_bond_count(
    mut single_atom: ExtendedMolecule,
    #[case] code: i8,
    #[case] expected: Option<RingBondCount>,
    flags_lenient: CtabParseFlags,
) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::RingBondCountEntries(vec![RingBondCountEntry {
        atom_index: 0,
        ring_bond_count: code,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();
    assert_eq!(single_atom.atoms[0].ring_bond_count, expected);
}

#[rstest]
#[case(1)]
#[case(5)]
fn test_apply_ring_bond_count_invalid(#[case] code: i8, flags_strict: CtabParseFlags) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::RingBondCountEntries(vec![RingBondCountEntry {
        atom_index: 0,
        ring_bond_count: code,
    }]);
    let result = acc.add_entry(entry, flags_strict);
    assert!(result.is_err());
}

#[rstest]
#[case(0, None)]
#[case(2, Some(SubstitutionCount::S2))]
#[case(4, Some(SubstitutionCount::S4))]
#[case(-1, Some(SubstitutionCount::NoSubstitution))]
fn test_apply_substitution_count(
    mut single_atom: ExtendedMolecule,
    #[case] code: i8,
    #[case] expected: Option<SubstitutionCount>,
    flags_lenient: CtabParseFlags,
) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::SubstitutionCountEntries(vec![SubstitutionCountEntry {
        atom_index: 0,
        substitution_count: code,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();
    assert_eq!(single_atom.atoms[0].substitution_count, expected);
}

#[rstest]
#[case(-3)]
#[case(7)]
fn test_apply_substitution_count_invalid(#[case] code: i8, flags_strict: CtabParseFlags) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::SubstitutionCountEntries(vec![SubstitutionCountEntry {
        atom_index: 0,
        substitution_count: code,
    }]);
    let result = acc.add_entry(entry, flags_strict);
    assert!(result.is_err());
}

#[rstest]
#[case(0, None)]
#[case(1, Some(UnsaturatedAtom))]
fn test_apply_unsaturated(
    mut single_atom: ExtendedMolecule,
    #[case] code: u8,
    #[case] expected: Option<UnsaturatedAtom>,
    flags_lenient: CtabParseFlags,
) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::UnsaturatedAtomEntries(vec![UnsaturatedAtomEntry {
        atom_index: 0,
        unsaturated: code,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();
    assert_eq!(single_atom.atoms[0].unsaturated, expected);
}

#[rstest]
#[case(2)]
fn test_apply_unsaturated_invalid(#[case] code: u8, flags_strict: CtabParseFlags) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::UnsaturatedAtomEntries(vec![UnsaturatedAtomEntry {
        atom_index: 0,
        unsaturated: code,
    }]);
    let result = acc.add_entry(entry, flags_strict);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_link_atom(mut single_atom: ExtendedMolecule, flags_lenient: CtabParseFlags) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::LinkAtomEntries(vec![LinkAtomEntry {
        atom_index: 0,
        repeat_count: 2,
        subs_index1: 1,
        subs_index2: None,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();
    assert_eq!(
        single_atom.atoms[0].link_atom,
        Some(LinkAtom {
            repeat_count: 2,
            subs_index1: 1,
            subs_index2: None
        })
    );
}

#[rstest]
fn test_apply_link_atom_conflict(flags_lenient: CtabParseFlags) {
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
    let result = acc.add_entry(entry, flags_lenient);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_atom_list(mut single_atom: ExtendedMolecule, flags_lenient: CtabParseFlags) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::AtomListEntry(AtomListEntry {
        atom_index: 0,
        elements: vec![e!(N), e!(O)],
        exclusion: false,
    });
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();
    assert_eq!(
        single_atom.atoms[0].symbol,
        AtomSymbol::AtomList(AtomList {
            elements: vec![e!(N), e!(O)],
            exclusion: false
        })
    );
}

#[rstest]
#[case(1, Some(AttachmentPointType::First))]
#[case(2, Some(AttachmentPointType::Second))]
#[case(3, Some(AttachmentPointType::Both))]
#[case(0, None)]
fn test_apply_attachment_point(
    mut single_atom: ExtendedMolecule,
    #[case] code: u8,
    #[case] expected: Option<AttachmentPointType>,
    flags_lenient: CtabParseFlags,
) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::AttachmentPointEntries(vec![AttachmentPointEntry {
        atom_index: 0,
        attachment_type: code,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();
    assert_eq!(single_atom.atoms[0].attachment_point, expected);
}

#[rstest]
fn test_apply_attachment_point_conflict(flags_lenient: CtabParseFlags) {
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
    let result = acc.add_entry(entry, flags_lenient);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_attachment_point_invalid(flags_strict: CtabParseFlags) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::AttachmentPointEntries(vec![AttachmentPointEntry {
        atom_index: 0,
        attachment_type: 4,
    }]);
    let result = acc.add_entry(entry, flags_strict);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_attachment_order(mut single_atom: ExtendedMolecule, flags_lenient: CtabParseFlags) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::AtomAttachmentOrderEntry(AtomAttachmentOrderEntry {
        atom_index: 0,
        attachments: vec![(1, 2), (2, 1)],
    });
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();
    assert_eq!(
        single_atom.atoms[0].attachment_order.as_ref().unwrap(),
        &vec![(1, 2), (2, 1)]
    );
}

#[rstest]
fn test_apply_rgroup_label(mut single_atom: ExtendedMolecule, flags_lenient: CtabParseFlags) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::RGroupLabelEntries(vec![RGroupLabelEntry {
        atom_index: 0,
        label: 1,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();
    // Accumulator creates minimal RGroup with empty occurrence
    if let AtomSymbol::RGroup(rgroup) = &single_atom.atoms[0].symbol {
        assert_eq!(rgroup.label, Some(1));
    } else {
        panic!("Expected RGroup symbol");
    }
}

#[rstest]
fn test_apply_rgroup_label_keep(mut with_rgroup: ExtendedMolecule, flags_lenient: CtabParseFlags) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::RGroupLabelEntries(vec![RGroupLabelEntry {
        atom_index: 0,
        label: 1,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut with_rgroup, flags_lenient)
        .unwrap();

    // Keeps original label (doesn't overwrite same label)
    if let AtomSymbol::RGroup(rgroup) = &with_rgroup.atoms[0].symbol {
        assert_eq!(rgroup.label, Some(1));
    } else {
        panic!("Expected RGroup symbol");
    }
}

#[rstest]
fn test_apply_rgroup_label_overwrite(
    mut with_unlabeled_rgroup: ExtendedMolecule,
    flags_lenient: CtabParseFlags,
) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::RGroupLabelEntries(vec![RGroupLabelEntry {
        atom_index: 0,
        label: 3,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut with_unlabeled_rgroup, flags_lenient)
        .unwrap();

    // Label is overwritten
    if let AtomSymbol::RGroup(rgroup) = &with_unlabeled_rgroup.atoms[0].symbol {
        assert_eq!(rgroup.label, Some(3));
    } else {
        panic!("Expected RGroup symbol");
    }
}

#[rstest]
fn test_apply_rgroup_logic(mut with_rgroup: ExtendedMolecule, flags_lenient: CtabParseFlags) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::RGroupLogicEntry(RGroupLogicEntry {
        label: 1,
        dependent_label: Some(2),
        rgroup_or_h: true,
        occurrence: vec![RGroupOccurrence::Exactly(1)],
    });
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut with_rgroup, flags_lenient)
        .unwrap();

    let rgroup = with_rgroup.rgroups.get(&1).unwrap();
    assert_eq!(rgroup.dependent_label, Some(2));
    assert!(rgroup.rgroup_or_h);
    assert_eq!(rgroup.occurrence.len(), 1);
    assert_eq!(rgroup.occurrence[0], RGroupOccurrence::Exactly(1));
}

#[rstest]
fn test_apply_rgroup_logic_multiple_occurrences(
    mut with_rgroup: ExtendedMolecule,
    flags_lenient: CtabParseFlags,
) {
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
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut with_rgroup, flags_lenient)
        .unwrap();

    let rgroup = with_rgroup.rgroups.get(&1).unwrap();
    assert_eq!(rgroup.occurrence.len(), 2);
    assert_eq!(rgroup.occurrence[0], RGroupOccurrence::Exactly(1));
    assert_eq!(rgroup.occurrence[1], RGroupOccurrence::GreaterThan(5));
}

#[rstest]
fn test_apply_zero_order_bond(mut with_bond: ExtendedMolecule, flags_lenient: CtabParseFlags) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::ZeroBondOrderEntries(vec![ZeroBondOrderEntry {
        bond_index: 0,
        bond_order: 0,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut with_bond, flags_lenient)
        .unwrap();

    assert_eq!(with_bond.bonds[0].order, BondOrder::Zero);
}

#[rstest]
fn test_apply_zero_atom_charge(mut single_atom: ExtendedMolecule, flags_lenient: CtabParseFlags) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::ZeroAtomChargeEntries(vec![ZeroAtomChargeEntry {
        atom_index: 0,
        charge: -1,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    assert_eq!(single_atom.atoms[0].charge, Some(-1));
}

#[rstest]
fn test_apply_hydrogen_count(mut single_atom: ExtendedMolecule, flags_lenient: CtabParseFlags) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::AtomHydrogenCountEntries(vec![AtomHydrogenCountEntry {
        atom_index: 0,
        hydrogen_count: 1,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    assert_eq!(single_atom.atoms[0].hydrogens, Some(1));
}

#[rstest]
fn test_apply_sgroup_type(mut single_atom: ExtendedMolecule, flags_lenient: CtabParseFlags) {
    let mut acc = MoleculeProperties::new();
    let entry = PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry {
        sgroup_index: 0,
        sgroup_type: SGroupType::Superatom,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    assert_eq!(single_atom.sgroups.len(), 1);
    let sgroup = single_atom.sgroups.get(&0).unwrap();
    assert_eq!(sgroup.group_type, SGroupType::Superatom);
}

#[rstest]
fn test_apply_sgroup_type_conflict(flags_lenient: CtabParseFlags) {
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
    let result = acc.add_entry(entry, flags_lenient);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_sgroup_subtype(
    mut single_atom: ExtendedMolecule,
    mut acc_with_copolymer_sgroup: MoleculeProperties,
    flags_lenient: CtabParseFlags,
) {
    let subtype_entry = PropertyEntries::SGroupSubtypeEntries(vec![SGroupSubtypeEntry {
        sgroup_index: 0,
        sgroup_subtype: SGroupSubtype::Alternating,
    }]);
    acc_with_copolymer_sgroup
        .add_entry(subtype_entry, flags_lenient)
        .unwrap();
    acc_with_copolymer_sgroup
        .update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    let sgroup = single_atom.sgroups.get(&0).unwrap();
    assert_eq!(sgroup.group_subtype, Some(SGroupSubtype::Alternating));
}

#[rstest]
fn test_apply_sgroup_subtype_missing(flags_strict: CtabParseFlags) {
    let mut acc = MoleculeProperties::new();
    let subtype_entry = PropertyEntries::SGroupSubtypeEntries(vec![SGroupSubtypeEntry {
        sgroup_index: 0,
        sgroup_subtype: SGroupSubtype::Alternating,
    }]);
    let result = acc.add_entry(subtype_entry, flags_strict);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_sgroup_label(
    mut single_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
    flags_lenient: CtabParseFlags,
) {
    let label_entry = PropertyEntries::SGroupLabelEntries(vec![SGroupLabelEntry {
        sgroup_index: 0,
        label: 123,
    }]);
    acc_with_superatom_sgroup
        .add_entry(label_entry, flags_lenient)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    let sgroup = single_atom.sgroups.get(&0).unwrap();
    assert_eq!(sgroup.label, Some(123));
}

#[rstest]
fn test_apply_sgroup_connectivity(
    mut single_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
    flags_lenient: CtabParseFlags,
) {
    let entry = PropertyEntries::SGroupConnectivityEntries(vec![SGroupConnectivityEntry {
        sgroup_index: 0,
        connectivity: SGroupConnectivity::HeadToTail,
    }]);
    acc_with_superatom_sgroup
        .add_entry(entry, flags_lenient)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    let sgroup = single_atom.sgroups.get(&0).unwrap();
    assert_eq!(sgroup.connectivity, Some(SGroupConnectivity::HeadToTail));
}

#[rstest]
fn test_apply_sgroup_expansion(
    mut single_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
    flags_lenient: CtabParseFlags,
) {
    let entry =
        PropertyEntries::SGroupExpansionEntries(vec![SGroupExpansionEntry { sgroup_index: 0 }]);
    acc_with_superatom_sgroup
        .add_entry(entry, flags_lenient)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    let sgroup = single_atom.sgroups.get(&0).unwrap();
    assert!(sgroup.expansion);
}

#[rstest]
fn test_apply_sgroup_atom_list(
    mut single_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
    flags_lenient: CtabParseFlags,
) {
    let entry = PropertyEntries::SGroupAtomListEntry(SGroupAtomListEntry {
        sgroup_index: 0,
        atom_indices: vec![0, 1],
    });
    acc_with_superatom_sgroup
        .add_entry(entry, flags_lenient)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    let sgroup = single_atom.sgroups.get(&0).unwrap();
    assert_eq!(sgroup.atom_indices, vec![0, 1]);
}

#[rstest]
fn test_apply_sgroup_bond_list(
    mut single_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
    flags_lenient: CtabParseFlags,
) {
    let entry = PropertyEntries::SGroupBondListEntry(SGroupBondListEntry {
        sgroup_index: 0,
        bond_indices: vec![0, 1],
    });
    acc_with_superatom_sgroup
        .add_entry(entry, flags_lenient)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    let sgroup = single_atom.sgroups.get(&0).unwrap();
    assert_eq!(sgroup.bond_indices, vec![0, 1]);
}

#[rstest]
fn test_apply_sgroup_parent_atom(
    mut single_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
    flags_lenient: CtabParseFlags,
) {
    let entry = PropertyEntries::SGroupParentAtomEntry(SGroupParentAtomEntry {
        sgroup_index: 0,
        atom_indices: vec![0, 1],
    });
    acc_with_superatom_sgroup
        .add_entry(entry, flags_lenient)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    let sgroup = single_atom.sgroups.get(&0).unwrap();
    assert_eq!(sgroup.parent_atom_indices, Some(vec![0, 1]));
}

#[rstest]
fn test_apply_sgroup_subscript(
    mut single_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
    flags_lenient: CtabParseFlags,
) {
    let entry = PropertyEntries::SGroupSubscriptEntry(SGroupSubscriptEntry {
        sgroup_index: 0,
        multiplier: None,
        subscript: Some("Ph".to_string()),
    });
    acc_with_superatom_sgroup
        .add_entry(entry, flags_lenient)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    let sgroup = single_atom.sgroups.get(&0).unwrap();
    assert_eq!(sgroup.subscript, Some("Ph".to_string()));
    assert_eq!(sgroup.multiplier, None);
}

#[rstest]
fn test_apply_sgroup_multiplier(
    mut single_atom: ExtendedMolecule,
    mut acc_with_multiple_sgroup: MoleculeProperties,
    flags_lenient: CtabParseFlags,
) {
    let entry = PropertyEntries::SGroupSubscriptEntry(SGroupSubscriptEntry {
        sgroup_index: 0,
        multiplier: Some(SGroupMultiplier::Single(SGroupMultiplierTerm::Variable(
            'n',
        ))),
        subscript: None,
    });
    acc_with_multiple_sgroup
        .add_entry(entry, flags_lenient)
        .unwrap();
    acc_with_multiple_sgroup
        .update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    let sgroup = single_atom.sgroups.get(&0).unwrap();
    assert_eq!(sgroup.subscript, None);
    assert_eq!(
        sgroup.multiplier,
        Some(SGroupMultiplier::Single(SGroupMultiplierTerm::Variable(
            'n'
        )))
    );
}

#[rstest]
fn test_apply_sgroup_correspondence(
    mut single_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
    flags_lenient: CtabParseFlags,
) {
    let entry = PropertyEntries::SGroupCorrespondenceEntry(SGroupCorrespondenceEntry {
        sgroup_index: 0,
        bond_indices: vec![0, 1],
    });
    acc_with_superatom_sgroup
        .add_entry(entry, flags_lenient)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    let sgroup = single_atom.sgroups.get(&0).unwrap();
    assert_eq!(sgroup.correspondence, Some(vec![0, 1]));
}

#[rstest]
fn test_apply_sgroup_display_info(
    mut single_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
    flags_lenient: CtabParseFlags,
) {
    let entry = PropertyEntries::SGroupDisplayInfoEntry(SGroupDisplayInfoEntry {
        sgroup_index: 0,
        bracket_coords: vec![1.0, 2.0, 3.0, 4.0],
    });
    acc_with_superatom_sgroup
        .add_entry(entry, flags_lenient)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    let sgroup = single_atom.sgroups.get(&0).unwrap();
    assert_eq!(
        sgroup.bracket_coords,
        Some(SGroupBracketCoords {
            bracket1: (1.0, 2.0),
            bracket2: (3.0, 4.0)
        })
    );
}

#[rstest]
fn test_apply_sgroup_connecting_bond(
    mut single_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
    flags_lenient: CtabParseFlags,
) {
    let entry = PropertyEntries::SGroupConnectingBondEntry(SGroupConnectingBondEntry {
        sgroup_index: 0,
        bond_index: 0,
        bond_vector: (1.0, 2.0),
    });
    acc_with_superatom_sgroup
        .add_entry(entry, flags_lenient)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    let sgroup = single_atom.sgroups.get(&0).unwrap();
    assert_eq!(
        sgroup.connecting_bond,
        Some(SGroupConnectingBond {
            bond_index: 0,
            bond_vector: (1.0, 2.0)
        })
    );
}

#[rstest]
fn test_apply_sgroup_data_description(
    mut single_atom: ExtendedMolecule,
    mut acc_with_data_sgroup: MoleculeProperties,
    flags_lenient: CtabParseFlags,
) {
    acc_with_data_sgroup
        .update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    let sgroup = single_atom.sgroups.get(&0).unwrap();
    let data = sgroup.data.get("test").unwrap();
    assert_eq!(data.field_type, SGroupDataType::Text);
    assert_eq!(data.field_units, None);
    assert_eq!(data.query_identifier, None);
    assert_eq!(data.data_query_operator, None);
}

#[rstest]
fn test_apply_sgroup_data_entry(
    mut single_atom: ExtendedMolecule,
    mut acc_with_data_sgroup: MoleculeProperties,
    flags_lenient: CtabParseFlags,
) {
    let continuation_entry = PropertyEntries::SGroupDataEntry(SGroupDataEntry::Continuation {
        sgroup_index: 0,
        data_content: "content".to_string(),
    });
    acc_with_data_sgroup
        .add_entry(continuation_entry, flags_lenient)
        .unwrap();

    let data_entry =
        PropertyEntries::SGroupDataEntry(SGroupDataEntry::EndBlank { sgroup_index: 0 });
    acc_with_data_sgroup
        .add_entry(data_entry, flags_lenient)
        .unwrap();
    acc_with_data_sgroup
        .update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    let sgroup = single_atom.sgroups.get(&0).unwrap();
    let data = sgroup.data.get("test").unwrap();
    assert_eq!(data.data_content, Some(vec!["content".to_string()]));
}

#[rstest]
fn test_apply_sgroup_data_display(
    mut single_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
    flags_lenient: CtabParseFlags,
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
    acc_with_superatom_sgroup
        .add_entry(display_entry, flags_lenient)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    let sgroup = single_atom.sgroups.get(&0).unwrap();
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
fn test_apply_sgroup_hierarchy(
    mut single_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
    flags_lenient: CtabParseFlags,
) {
    // Add a second SGroup as parent
    let parent_type_entry = PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry {
        sgroup_index: 1,
        sgroup_type: SGroupType::Superatom,
    }]);
    acc_with_superatom_sgroup
        .add_entry(parent_type_entry, flags_lenient)
        .unwrap();

    let hierarchy_entry = PropertyEntries::SGroupHierarchyEntries(vec![SGroupHierarchyEntry {
        sgroup_index: 0,
        parent_sgroup_index: 1,
    }]);
    acc_with_superatom_sgroup
        .add_entry(hierarchy_entry, flags_lenient)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    assert_eq!(single_atom.sgroups.len(), 2);
    let child_sgroup = single_atom.sgroups.get(&0).unwrap();
    assert_eq!(child_sgroup.hierarchy_parent, Some(1));
}

#[rstest]
fn test_apply_sgroup_component(
    mut single_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: MoleculeProperties,
    flags_lenient: CtabParseFlags,
) {
    let component_entry = PropertyEntries::SGroupComponentEntries(vec![SGroupComponentEntry {
        sgroup_index: 0,
        component_number: 42,
    }]);
    acc_with_superatom_sgroup
        .add_entry(component_entry, flags_lenient)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    let sgroup = single_atom.sgroups.get(&0).unwrap();
    assert_eq!(sgroup.component_number, Some(42));
}

