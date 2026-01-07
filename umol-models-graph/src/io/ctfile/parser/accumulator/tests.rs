//! Unit tests for accumulator applying properties to Molecule and ExtendedMolecule

use pretty_assertions::assert_eq;
use rstest::*;
use umol_data::{e, Element, NamedIsotope};

use super::*;
use crate::io::ctfile::config::CtabParseFlags;
use crate::io::ctfile::parser::properties::{
    AtomAliasEntry, AtomAttachmentOrderEntry, AtomChargeOverrideEntry, AtomHydrogenCountEntry,
    AtomListEntry, AtomValueEntry, AttachmentPointEntry, BondOrderOverrideEntry, ChargeEntry,
    IsotopeEntry, LegacyGroupAbbreviationEntry, LinkAtomEntry, MoleculeChiralFlagEntry,
    PropertyEntries, RGroupLabelEntry, RGroupLogicEntry, RadicalEntry, RingBondCountEntry,
    SGroupAtomListEntry, SGroupBondListEntry, SGroupComponentEntry, SGroupConnectingBondEntry,
    SGroupConnectivityEntry, SGroupCorrespondenceEntry, SGroupDataDescriptionEntry,
    SGroupDataDisplayEntry, SGroupDataEntry, SGroupDisplayInfoEntry, SGroupExpansionEntry,
    SGroupHierarchyEntry, SGroupLabelEntry, SGroupParentAtomEntry, SGroupSubscriptEntry,
    SGroupSubtypeEntry, SGroupTypeEntry, SubstitutionCountEntry, UnsaturatedAtomEntry,
};
use crate::table_ir::{
    Atom, AtomList, AtomSymbol, AttachmentPointType, Bond, BondOrder, ExtendedAtom, ExtendedBond,
    ExtendedMolecule, LinkAtom, Molecule, RGroup, RGroupOccurrence, RingBondCount,
    SGroupDataDisplayChars, SGroupDataDisplayPlacement, SGroupDataDisplayType,
    SGroupDataDisplayUnits, SGroupDataType, SGroupMultiplier, SGroupMultiplierTerm, SGroupType,
    SubstitutionCount, UnsaturatedAtom,
};

#[fixture]
fn flags_basic() -> CtabParseFlags {
    CtabParseFlags::BASIC
}

#[fixture]
fn flags_extended() -> CtabParseFlags {
    CtabParseFlags::EXTENDED
}

#[fixture]
fn flags_strict() -> CtabParseFlags {
    CtabParseFlags::STRICT
}

#[fixture]
fn flags_lenient() -> CtabParseFlags {
    CtabParseFlags::LENIENT
}

#[fixture]
fn single_atom() -> Molecule {
    let mut molecule = Molecule::empty();
    molecule.atoms.push(Atom::from_element(e!(C)));
    molecule
}

#[fixture]
fn triatomic_molecule() -> Molecule {
    let mut molecule = Molecule::empty();
    for _ in 0..3 {
        molecule.atoms.push(Atom::from_element(e!(C)));
    }
    molecule
}

#[fixture]
fn with_properties(mut single_atom: Molecule) -> Molecule {
    if let Some(atom) = single_atom.atoms.get_mut(0) {
        atom.charge = Some(1);
        atom.unpaired_e = Some(1); // Doublet: 1 unpaired electron
        atom.isotope_mass = Some(13);
        atom.alias = Some("existing".to_string());
        atom.value = Some("existing".to_string());
    }
    single_atom
}

#[fixture]
fn with_bond(mut triatomic_molecule: Molecule) -> Molecule {
    triatomic_molecule
        .bonds
        .push(Bond::new(0, 1, BondOrder::Single));
    triatomic_molecule
}

#[fixture]
fn single_extended_atom() -> ExtendedMolecule {
    let mut molecule = ExtendedMolecule::empty();
    molecule.atoms.push(ExtendedAtom::from_element(e!(C)));
    molecule
}

#[fixture]
fn triatomic_extended_molecule() -> ExtendedMolecule {
    let mut molecule = ExtendedMolecule::empty();
    for _ in 0..3 {
        molecule.atoms.push(ExtendedAtom::from_element(e!(C)));
    }
    molecule
}

#[fixture]
fn with_extended_properties(mut single_extended_atom: ExtendedMolecule) -> ExtendedMolecule {
    if let Some(atom) = single_extended_atom.atoms.get_mut(0) {
        atom.properties
            .insert("molFileAlias".to_string(), "existing".to_string());
        atom.properties
            .insert("molFileValue".to_string(), "existing".to_string());
        atom.charge = Some(1);
        atom.unpaired_e = Some(1); // Doublet: 1 unpaired electron
        atom.isotope_mass = Some(13);
        atom.ring_bond_count = Some(RingBondCount::R2);
        atom.substitution_count = Some(SubstitutionCount::S2);
        atom.unsaturated = Some(UnsaturatedAtom);
    }
    single_extended_atom
}

#[fixture]
fn with_rgroup(mut single_extended_atom: ExtendedMolecule) -> ExtendedMolecule {
    single_extended_atom.atoms[0].symbol = AtomSymbol::RGroup(RGroup::new(Some(1)));
    single_extended_atom
}

#[fixture]
fn with_unlabeled_rgroup(mut single_extended_atom: ExtendedMolecule) -> ExtendedMolecule {
    single_extended_atom.atoms[0].symbol = AtomSymbol::RGroup(RGroup::new(None));
    single_extended_atom
}

#[fixture]
fn with_extended_bond(mut single_extended_atom: ExtendedMolecule) -> ExtendedMolecule {
    single_extended_atom
        .atoms
        .push(ExtendedAtom::from_element(e!(C)));
    single_extended_atom
        .atoms
        .push(ExtendedAtom::from_element(e!(C)));
    single_extended_atom
        .bonds
        .push(ExtendedBond::new(0, 1, BondOrder::Single));
    single_extended_atom
}

#[fixture]
fn acc_with_superatom_sgroup(flags_extended: CtabParseFlags) -> PropertyAccumulator {
    let mut acc = PropertyAccumulator::new();
    let type_entry = PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry {
        sgroup_index: 0,
        sgroup_type: SGroupType::Superatom,
    }]);
    acc.add_entry(type_entry, flags_extended).unwrap();
    acc
}

#[fixture]
fn acc_with_data_sgroup(flags_extended: CtabParseFlags) -> PropertyAccumulator {
    let mut acc = PropertyAccumulator::new();
    let type_entry = PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry {
        sgroup_index: 0,
        sgroup_type: SGroupType::Data,
    }]);
    acc.add_entry(type_entry, flags_extended).unwrap();
    let data_description_entry =
        PropertyEntries::SGroupDataDescriptionEntry(SGroupDataDescriptionEntry {
            sgroup_index: 0,
            field_name: "test".to_string(),
            field_type: SGroupDataType::Text,
            field_units: None,
            query_identifier: None,
            data_query_operator: None,
        });
    acc.add_entry(data_description_entry, flags_extended)
        .unwrap();
    acc
}

#[fixture]
fn acc_with_multiple_sgroup(flags_extended: CtabParseFlags) -> PropertyAccumulator {
    let mut acc = PropertyAccumulator::new();
    let type_entry = PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry {
        sgroup_index: 0,
        sgroup_type: SGroupType::MultipleGroup,
    }]);
    acc.add_entry(type_entry, flags_extended).unwrap();
    acc
}

#[fixture]
fn acc_with_copolymer_sgroup(flags_extended: CtabParseFlags) -> PropertyAccumulator {
    let mut acc = PropertyAccumulator::new();
    let type_entry = PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry {
        sgroup_index: 0,
        sgroup_type: SGroupType::Copolymer,
    }]);
    acc.add_entry(type_entry, flags_extended).unwrap();
    acc
}

#[rstest]
fn test_apply_molecule_chiral_flag(mut single_atom: Molecule, flags_basic: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    let entry =
        PropertyEntries::MoleculeChiralFlagEntry(MoleculeChiralFlagEntry { chiral_flag: true });
    acc.add_entry(entry, flags_basic).unwrap();
    acc.update_molecule(&mut single_atom, flags_basic).unwrap();

    assert_eq!(
        single_atom.properties.get("chiral_flag"),
        Some(&"true".to_string())
    );
}

#[rstest]
fn test_apply_molecule_no_chiral_flag(mut single_atom: Molecule, flags_basic: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    acc.update_molecule(&mut single_atom, flags_basic).unwrap();

    assert!(single_atom.properties.get("chiral_flag").is_none());
}

#[rstest]
fn test_apply_atom_alias(mut single_atom: Molecule, flags_basic: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::AtomAliasEntry(AtomAliasEntry {
        atom_index: 0,
        alias: "CF3".to_string(),
    });
    acc.add_entry(entry, flags_basic).unwrap();
    acc.update_molecule(&mut single_atom, flags_basic).unwrap();

    let atom = &single_atom.atoms[0];
    assert_eq!(atom.alias, Some("CF3".to_string()));
}

#[rstest]
fn test_apply_atom_alias_invalid_index(mut single_atom: Molecule, flags_basic: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::AtomAliasEntry(AtomAliasEntry {
        atom_index: 5,
        alias: "CF3".to_string(),
    });
    acc.add_entry(entry, flags_basic).unwrap();
    let result = acc.update_molecule(&mut single_atom, flags_basic);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_atom_value(mut single_atom: Molecule, flags_basic: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::AtomValueEntry(AtomValueEntry {
        atom_index: 0,
        value: "*".to_string(),
    });
    acc.add_entry(entry, flags_basic).unwrap();
    acc.update_molecule(&mut single_atom, flags_basic).unwrap();

    let atom = &single_atom.atoms[0];
    assert_eq!(atom.value, Some("*".to_string()));
}

#[rstest]
fn test_apply_atom_value_invalid_index(mut single_atom: Molecule, flags_basic: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::AtomValueEntry(AtomValueEntry {
        atom_index: 5,
        value: "*".to_string(),
    });
    acc.add_entry(entry, flags_basic).unwrap();
    let result = acc.update_molecule(&mut single_atom, flags_basic);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_charge(mut single_atom: Molecule, flags_basic: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::ChargeEntries(vec![ChargeEntry {
        atom_index: 0,
        charge: -1,
    }]);
    acc.add_entry(entry, flags_basic).unwrap();
    acc.update_molecule(&mut single_atom, flags_basic).unwrap();

    let atom = &single_atom.atoms[0];
    assert_eq!(atom.charge, Some(-1));
    assert_eq!(atom.unpaired_e, None);
}

#[rstest]
fn test_apply_charge_multiple(mut triatomic_molecule: Molecule, flags_basic: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();

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
    acc.add_entry(PropertyEntries::ChargeEntries(entries), flags_basic)
        .unwrap();
    acc.update_molecule(&mut triatomic_molecule, flags_basic)
        .unwrap();

    assert_eq!(triatomic_molecule.atoms[0].charge, Some(-1));
    assert_eq!(triatomic_molecule.atoms[1].charge, None);
    assert_eq!(triatomic_molecule.atoms[2].charge, Some(1));
}

#[rstest]
fn test_apply_charge_invalid_index(mut single_atom: Molecule, flags_basic: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::ChargeEntries(vec![ChargeEntry {
        atom_index: 5,
        charge: -1,
    }]);
    acc.add_entry(entry, flags_basic).unwrap();
    let result = acc.update_molecule(&mut single_atom, flags_basic);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_charge_overwrite(mut with_properties: Molecule, flags_basic: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::ChargeEntries(vec![ChargeEntry {
        atom_index: 0,
        charge: -2,
    }]);
    acc.add_entry(entry, flags_basic).unwrap();
    acc.update_molecule(&mut with_properties, flags_basic)
        .unwrap();

    let atom = &with_properties.atoms[0];
    assert_eq!(atom.charge, Some(-2));
    assert_eq!(atom.unpaired_e, None);
}

#[rstest]
fn test_apply_radical(mut single_atom: Molecule, flags_basic: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::RadicalEntries(vec![RadicalEntry {
        atom_index: 0,
        radical_type: 2, // Doublet: 1 unpaired electron
    }]);
    acc.add_entry(entry, flags_basic).unwrap();
    acc.update_molecule(&mut single_atom, flags_basic).unwrap();

    let atom = &single_atom.atoms[0];
    assert_eq!(atom.unpaired_e, Some(1));
    assert_eq!(atom.charge, None);
}

#[rstest]
fn test_apply_radical_invalid_index(mut single_atom: Molecule, flags_basic: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::RadicalEntries(vec![RadicalEntry {
        atom_index: 5,
        radical_type: 2,
    }]);
    acc.add_entry(entry, flags_basic).unwrap();
    let result = acc.update_molecule(&mut single_atom, flags_basic);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_radical_invalid_code(flags_basic: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::RadicalEntries(vec![RadicalEntry {
        atom_index: 0,
        radical_type: 4,
    }]);
    let result = acc.add_entry(entry, flags_basic);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_radical_overwrite(mut with_properties: Molecule, flags_basic: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::RadicalEntries(vec![RadicalEntry {
        atom_index: 0,
        radical_type: 1,
    }]);
    acc.add_entry(entry, flags_basic).unwrap();
    acc.update_molecule(&mut with_properties, flags_basic)
        .unwrap();

    let atom = &with_properties.atoms[0];
    assert_eq!(atom.unpaired_e, Some(0)); // Singlet: 0 unpaired electrons
    assert_eq!(atom.charge, None);
}

#[rstest]
fn test_apply_isotope(mut single_atom: Molecule, flags_basic: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::IsotopeEntries(vec![IsotopeEntry {
        atom_index: 0,
        mass: 14,
    }]);
    acc.add_entry(entry, flags_basic).unwrap();
    acc.update_molecule(&mut single_atom, flags_basic).unwrap();

    let atom = &single_atom.atoms[0];
    assert_eq!(atom.isotope_mass, Some(14));
}

#[rstest]
fn test_apply_isotope_invalid_index(mut single_atom: Molecule, flags_basic: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::IsotopeEntries(vec![IsotopeEntry {
        atom_index: 5,
        mass: 13,
    }]);
    acc.add_entry(entry, flags_basic).unwrap();
    let result = acc.update_molecule(&mut single_atom, flags_basic);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_isotope_lenient(mut single_atom: Molecule, flags_lenient: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::IsotopeEntries(vec![IsotopeEntry {
        atom_index: 0,
        mass: 40,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    let atom = &single_atom.atoms[0];
    assert_eq!(atom.isotope_mass, Some(40));
}

#[rstest]
fn test_apply_bond_order_override(mut with_bond: Molecule, flags_lenient: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::BondOrderOverrideEntries(vec![BondOrderOverrideEntry {
        bond_index: 0,
        bond_order: BondOrder::Zero,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_molecule(&mut with_bond, flags_lenient).unwrap();

    let bond = &with_bond.bonds[0];
    assert_eq!(bond.order, BondOrder::Zero);
}

#[rstest]
fn test_apply_bond_order_override_invalid(
    mut single_atom: Molecule,
    flags_lenient: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::BondOrderOverrideEntries(vec![BondOrderOverrideEntry {
        bond_index: 5,
        bond_order: BondOrder::Zero,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    let result = acc.update_molecule(&mut single_atom, flags_lenient);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_atom_charge_override(mut single_atom: Molecule, flags_lenient: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::AtomChargeOverrideEntries(vec![AtomChargeOverrideEntry {
        atom_index: 0,
        charge: -1,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    let atom = &single_atom.atoms[0];
    assert_eq!(atom.charge, Some(-1));
}

#[rstest]
fn test_apply_atom_charge_override_invalid(
    mut single_atom: Molecule,
    flags_lenient: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::AtomChargeOverrideEntries(vec![AtomChargeOverrideEntry {
        atom_index: 5,
        charge: -1,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    let result = acc.update_molecule(&mut single_atom, flags_lenient);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_atom_hydrogen_count(mut single_atom: Molecule, flags_lenient: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::AtomHydrogenCountEntries(vec![AtomHydrogenCountEntry {
        atom_index: 0,
        hydrogen_count: Some(1),
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_molecule(&mut single_atom, flags_lenient)
        .unwrap();

    let atom = &single_atom.atoms[0];
    assert_eq!(atom.hydrogens, Some(1));
}

#[rstest]
fn test_apply_atom_hydrogen_count_invalid(
    mut single_atom: Molecule,
    flags_lenient: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::AtomHydrogenCountEntries(vec![AtomHydrogenCountEntry {
        atom_index: 5,
        hydrogen_count: Some(1),
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    let result = acc.update_molecule(&mut single_atom, flags_lenient);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_extended_molecule_chiral_flag(
    mut single_extended_atom: ExtendedMolecule,
    flags_extended: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry =
        PropertyEntries::MoleculeChiralFlagEntry(MoleculeChiralFlagEntry { chiral_flag: true });
    acc.add_entry(entry, flags_extended).unwrap();
    acc.update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();

    assert_eq!(
        single_extended_atom.properties.get("chiral_flag"),
        Some(&"true".to_string())
    );
}

#[rstest]
fn test_apply_extended_molecule_no_chiral_flag(
    mut single_extended_atom: ExtendedMolecule,
    flags_extended: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    acc.update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();

    assert!(single_extended_atom.properties.get("chiral_flag").is_none());
}

#[rstest]
fn test_apply_extended_atom_alias(
    mut single_extended_atom: ExtendedMolecule,
    flags_basic: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::AtomAliasEntry(AtomAliasEntry {
        atom_index: 0,
        alias: "CF3".to_string(),
    });
    acc.add_entry(entry, flags_basic).unwrap();
    acc.update_extended_molecule(&mut single_extended_atom, flags_basic)
        .unwrap();

    let atom = &single_extended_atom.atoms[0];
    assert_eq!(atom.alias, Some("CF3".to_string()));
}

#[rstest]
fn test_apply_extended_atom_alias_invalid_index(
    mut single_extended_atom: ExtendedMolecule,
    flags_basic: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::AtomAliasEntry(AtomAliasEntry {
        atom_index: 5,
        alias: "CF3".to_string(),
    });
    acc.add_entry(entry, flags_basic).unwrap();
    let result = acc.update_extended_molecule(&mut single_extended_atom, flags_basic);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_extended_legacy_group_abbreviation(
    mut triatomic_extended_molecule: ExtendedMolecule,
    flags_basic: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::LegacyGroupAbbreviationEntry(LegacyGroupAbbreviationEntry {
        atom_index1: 0,
        atom_index2: 1,
        label: "A".to_string(),
    });
    acc.add_entry(entry, flags_basic).unwrap();
    acc.update_extended_molecule(&mut triatomic_extended_molecule, flags_basic)
        .unwrap();

    let ctfile_data = triatomic_extended_molecule.ctfile_data.as_ref().unwrap();
    assert_eq!(ctfile_data.legacy_group_abbreviations.len(), 1);
    let abbreviation = &ctfile_data.legacy_group_abbreviations[0];
    assert_eq!(abbreviation.atom_index1, 0);
    assert_eq!(abbreviation.atom_index2, 1);
    assert_eq!(abbreviation.label, "A");
}

#[rstest]
fn test_apply_extended_atom_value(
    mut single_extended_atom: ExtendedMolecule,
    flags_basic: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::AtomValueEntry(AtomValueEntry {
        atom_index: 0,
        value: "*".to_string(),
    });
    acc.add_entry(entry, flags_basic).unwrap();
    acc.update_extended_molecule(&mut single_extended_atom, flags_basic)
        .unwrap();

    let atom = &single_extended_atom.atoms[0];
    assert_eq!(atom.value, Some("*".to_string()));
}

#[rstest]
fn test_apply_extended_atom_value_invalid_index(
    mut single_extended_atom: ExtendedMolecule,
    flags_basic: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::AtomValueEntry(AtomValueEntry {
        atom_index: 5,
        value: "*".to_string(),
    });
    acc.add_entry(entry, flags_basic).unwrap();
    let result = acc.update_extended_molecule(&mut single_extended_atom, flags_basic);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_extended_charge(
    mut single_extended_atom: ExtendedMolecule,
    flags_basic: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::ChargeEntries(vec![ChargeEntry {
        atom_index: 0,
        charge: -1,
    }]);
    acc.add_entry(entry, flags_basic).unwrap();
    acc.update_extended_molecule(&mut single_extended_atom, flags_basic)
        .unwrap();

    let atom = &single_extended_atom.atoms[0];
    assert_eq!(atom.charge, Some(-1));
    assert_eq!(atom.unpaired_e, None);
}

#[rstest]
fn test_apply_extended_charge_multiple(
    mut triatomic_extended_molecule: ExtendedMolecule,
    flags_basic: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();

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
    acc.add_entry(PropertyEntries::ChargeEntries(entries), flags_basic)
        .unwrap();
    acc.update_extended_molecule(&mut triatomic_extended_molecule, flags_basic)
        .unwrap();

    assert_eq!(triatomic_extended_molecule.atoms[0].charge, Some(-1));
    assert_eq!(triatomic_extended_molecule.atoms[1].charge, None);
    assert_eq!(triatomic_extended_molecule.atoms[2].charge, Some(1));
}

#[rstest]
fn test_apply_extended_charge_invalid_index(
    mut single_extended_atom: ExtendedMolecule,
    flags_basic: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::ChargeEntries(vec![ChargeEntry {
        atom_index: 5,
        charge: -1,
    }]);
    acc.add_entry(entry, flags_basic).unwrap();
    let result = acc.update_extended_molecule(&mut single_extended_atom, flags_basic);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_extended_charge_overwrite(
    mut with_extended_properties: ExtendedMolecule,
    flags_basic: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::ChargeEntries(vec![ChargeEntry {
        atom_index: 0,
        charge: -2,
    }]);
    acc.add_entry(entry, flags_basic).unwrap();
    acc.update_extended_molecule(&mut with_extended_properties, flags_basic)
        .unwrap();

    let atom = &with_extended_properties.atoms[0];
    assert_eq!(atom.charge, Some(-2));
    assert_eq!(atom.unpaired_e, None);
}

#[rstest]
fn test_apply_extended_radical(
    mut single_extended_atom: ExtendedMolecule,
    flags_basic: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::RadicalEntries(vec![RadicalEntry {
        atom_index: 0,
        radical_type: 2, // Doublet: 1 unpaired electron
    }]);
    acc.add_entry(entry, flags_basic).unwrap();
    acc.update_extended_molecule(&mut single_extended_atom, flags_basic)
        .unwrap();

    let atom = &single_extended_atom.atoms[0];
    assert_eq!(atom.unpaired_e, Some(1));
    assert_eq!(atom.charge, None);
}

#[rstest]
fn test_apply_extended_radical_invalid_index(
    mut single_extended_atom: ExtendedMolecule,
    flags_basic: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::RadicalEntries(vec![RadicalEntry {
        atom_index: 5,
        radical_type: 2,
    }]);
    acc.add_entry(entry, flags_basic).unwrap();
    let result = acc.update_extended_molecule(&mut single_extended_atom, flags_basic);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_extended_radical_invalid_code(flags_basic: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::RadicalEntries(vec![RadicalEntry {
        atom_index: 0,
        radical_type: 4,
    }]);
    let result = acc.add_entry(entry, flags_basic);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_extended_radical_overwrite(
    mut with_extended_properties: ExtendedMolecule,
    flags_basic: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::RadicalEntries(vec![RadicalEntry {
        atom_index: 0,
        radical_type: 1,
    }]);
    acc.add_entry(entry, flags_basic).unwrap();
    acc.update_extended_molecule(&mut with_extended_properties, flags_basic)
        .unwrap();

    let atom = &with_extended_properties.atoms[0];
    assert_eq!(atom.unpaired_e, Some(0)); // Singlet: 0 unpaired electrons
    assert_eq!(atom.charge, None);
}

#[rstest]
fn test_apply_extended_isotope(
    mut single_extended_atom: ExtendedMolecule,
    flags_basic: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::IsotopeEntries(vec![IsotopeEntry {
        atom_index: 0,
        mass: 14,
    }]);
    acc.add_entry(entry, flags_basic).unwrap();
    acc.update_extended_molecule(&mut single_extended_atom, flags_basic)
        .unwrap();

    let atom = &single_extended_atom.atoms[0];
    assert_eq!(atom.isotope_mass, Some(14));
}

#[rstest]
fn test_apply_extended_isotope_invalid_index(
    mut single_extended_atom: ExtendedMolecule,
    flags_basic: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::IsotopeEntries(vec![IsotopeEntry {
        atom_index: 5,
        mass: 13,
    }]);
    acc.add_entry(entry, flags_basic).unwrap();
    let result = acc.update_extended_molecule(&mut single_extended_atom, flags_basic);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_extended_isotope_named_isotope(
    mut single_extended_atom: ExtendedMolecule,
    flags_extended: CtabParseFlags,
) {
    single_extended_atom.atoms[0].symbol = AtomSymbol::NamedIsotope(NamedIsotope::D);
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::IsotopeEntries(vec![IsotopeEntry {
        atom_index: 0,
        mass: 3,
    }]);
    acc.add_entry(entry, flags_extended).unwrap();
    acc.update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();

    let atom = &single_extended_atom.atoms[0];
    assert_eq!(atom.isotope_mass, Some(3));
}

#[rstest]
fn test_apply_extended_isotope_lenient(
    mut single_extended_atom: ExtendedMolecule,
    flags_lenient: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::IsotopeEntries(vec![IsotopeEntry {
        atom_index: 0,
        mass: 40,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut single_extended_atom, flags_lenient)
        .unwrap();

    let atom = &single_extended_atom.atoms[0];
    assert_eq!(atom.isotope_mass, Some(40));
}

#[rstest]
fn test_apply_extended_bond_order_override(
    mut with_extended_bond: ExtendedMolecule,
    flags_lenient: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::BondOrderOverrideEntries(vec![BondOrderOverrideEntry {
        bond_index: 0,
        bond_order: BondOrder::Zero,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut with_extended_bond, flags_lenient)
        .unwrap();

    assert_eq!(with_extended_bond.bonds[0].order, BondOrder::Zero);
}

#[rstest]
fn test_apply_extended_bond_order_override_invalid(
    mut single_extended_atom: ExtendedMolecule,
    flags_lenient: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::BondOrderOverrideEntries(vec![BondOrderOverrideEntry {
        bond_index: 5,
        bond_order: BondOrder::Zero,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    let result = acc.update_extended_molecule(&mut single_extended_atom, flags_lenient);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_extended_atom_charge_override(
    mut single_extended_atom: ExtendedMolecule,
    flags_lenient: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::AtomChargeOverrideEntries(vec![AtomChargeOverrideEntry {
        atom_index: 0,
        charge: -1,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut single_extended_atom, flags_lenient)
        .unwrap();

    assert_eq!(single_extended_atom.atoms[0].charge, Some(-1));
}

#[rstest]
fn test_apply_extended_atom_charge_override_invalid(
    mut single_extended_atom: ExtendedMolecule,
    flags_lenient: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::AtomChargeOverrideEntries(vec![AtomChargeOverrideEntry {
        atom_index: 5,
        charge: -1,
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    let result = acc.update_extended_molecule(&mut single_extended_atom, flags_lenient);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_extended_atom_hydrogen_count(
    mut single_extended_atom: ExtendedMolecule,
    flags_lenient: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::AtomHydrogenCountEntries(vec![AtomHydrogenCountEntry {
        atom_index: 0,
        hydrogen_count: Some(1),
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    acc.update_extended_molecule(&mut single_extended_atom, flags_lenient)
        .unwrap();

    assert_eq!(single_extended_atom.atoms[0].hydrogens, Some(1));
}

#[rstest]
fn test_apply_extended_atom_hydrogen_count_invalid(
    mut single_extended_atom: ExtendedMolecule,
    flags_lenient: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::AtomHydrogenCountEntries(vec![AtomHydrogenCountEntry {
        atom_index: 5,
        hydrogen_count: Some(1),
    }]);
    acc.add_entry(entry, flags_lenient).unwrap();
    let result = acc.update_extended_molecule(&mut single_extended_atom, flags_lenient);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_extended_sgroup_type(
    mut single_extended_atom: ExtendedMolecule,
    flags_extended: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry {
        sgroup_index: 0,
        sgroup_type: SGroupType::Superatom,
    }]);
    acc.add_entry(entry, flags_extended).unwrap();
    acc.update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();

    assert_eq!(single_extended_atom.sgroups().len(), 1);
    let sgroup = single_extended_atom.sgroups().get(&0).unwrap();
    assert_eq!(sgroup.group_type, SGroupType::Superatom);
}

#[rstest]
fn test_apply_extended_sgroup_type_conflict(flags_extended: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
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
    let result = acc.add_entry(entry, flags_extended);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_extended_sgroup_subtype(
    mut single_extended_atom: ExtendedMolecule,
    mut acc_with_copolymer_sgroup: PropertyAccumulator,
    flags_extended: CtabParseFlags,
) {
    let subtype_entry = PropertyEntries::SGroupSubtypeEntries(vec![SGroupSubtypeEntry {
        sgroup_index: 0,
        sgroup_subtype: SGroupSubtype::Alternating,
    }]);
    acc_with_copolymer_sgroup
        .add_entry(subtype_entry, flags_extended)
        .unwrap();
    acc_with_copolymer_sgroup
        .update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();

    let sgroup = single_extended_atom.sgroups().get(&0).unwrap();
    assert_eq!(sgroup.group_subtype, Some(SGroupSubtype::Alternating));
}

#[rstest]
fn test_apply_extended_sgroup_subtype_missing(flags_strict: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    let subtype_entry = PropertyEntries::SGroupSubtypeEntries(vec![SGroupSubtypeEntry {
        sgroup_index: 0,
        sgroup_subtype: SGroupSubtype::Alternating,
    }]);
    let result = acc.add_entry(subtype_entry, flags_strict);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_extended_sgroup_label(
    mut single_extended_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: PropertyAccumulator,
    flags_extended: CtabParseFlags,
) {
    let label_entry = PropertyEntries::SGroupLabelEntries(vec![SGroupLabelEntry {
        sgroup_index: 0,
        label: 123,
    }]);
    acc_with_superatom_sgroup
        .add_entry(label_entry, flags_extended)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();

    let sgroup = single_extended_atom.sgroups().get(&0).unwrap();
    assert_eq!(sgroup.label, Some(123));
}

#[rstest]
fn test_apply_extended_sgroup_connectivity(
    mut single_extended_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: PropertyAccumulator,
    flags_extended: CtabParseFlags,
) {
    let entry = PropertyEntries::SGroupConnectivityEntries(vec![SGroupConnectivityEntry {
        sgroup_index: 0,
        connectivity: SGroupConnectivity::HeadToTail,
    }]);
    acc_with_superatom_sgroup
        .add_entry(entry, flags_extended)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();

    let sgroup = single_extended_atom.sgroups().get(&0).unwrap();
    assert_eq!(sgroup.connectivity, Some(SGroupConnectivity::HeadToTail));
}

#[rstest]
fn test_apply_extended_sgroup_expansion(
    mut single_extended_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: PropertyAccumulator,
    flags_extended: CtabParseFlags,
) {
    let entry =
        PropertyEntries::SGroupExpansionEntries(vec![SGroupExpansionEntry { sgroup_index: 0 }]);
    acc_with_superatom_sgroup
        .add_entry(entry, flags_extended)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();

    let sgroup = single_extended_atom.sgroups().get(&0).unwrap();
    assert!(sgroup.expansion);
}

#[rstest]
fn test_apply_extended_sgroup_atom_list(
    mut single_extended_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: PropertyAccumulator,
    flags_extended: CtabParseFlags,
) {
    let entry = PropertyEntries::SGroupAtomListEntry(SGroupAtomListEntry {
        sgroup_index: 0,
        atom_indices: vec![0, 1],
    });
    acc_with_superatom_sgroup
        .add_entry(entry, flags_extended)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();

    let sgroup = single_extended_atom.sgroups().get(&0).unwrap();
    assert_eq!(sgroup.atom_indices, vec![0, 1]);
}

#[rstest]
fn test_apply_extended_sgroup_bond_list(
    mut single_extended_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: PropertyAccumulator,
    flags_extended: CtabParseFlags,
) {
    let entry = PropertyEntries::SGroupBondListEntry(SGroupBondListEntry {
        sgroup_index: 0,
        bond_indices: vec![0, 1],
    });
    acc_with_superatom_sgroup
        .add_entry(entry, flags_extended)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();

    let sgroup = single_extended_atom.sgroups().get(&0).unwrap();
    assert_eq!(sgroup.bond_indices, vec![0, 1]);
}

#[rstest]
fn test_apply_extended_sgroup_parent_atom(
    mut single_extended_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: PropertyAccumulator,
    flags_extended: CtabParseFlags,
) {
    let entry = PropertyEntries::SGroupParentAtomEntry(SGroupParentAtomEntry {
        sgroup_index: 0,
        atom_indices: vec![0, 1],
    });
    acc_with_superatom_sgroup
        .add_entry(entry, flags_extended)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();

    let sgroup = single_extended_atom.sgroups().get(&0).unwrap();
    assert_eq!(sgroup.parent_atom_indices, Some(vec![0, 1]));
}

#[rstest]
fn test_apply_extended_sgroup_subscript(
    mut single_extended_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: PropertyAccumulator,
    flags_extended: CtabParseFlags,
) {
    let entry = PropertyEntries::SGroupSubscriptEntry(SGroupSubscriptEntry {
        sgroup_index: 0,
        multiplier: None,
        subscript: Some("Ph".to_string()),
    });
    acc_with_superatom_sgroup
        .add_entry(entry, flags_extended)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();

    let sgroup = single_extended_atom.sgroups().get(&0).unwrap();
    assert_eq!(sgroup.subscript, Some("Ph".to_string()));
    assert_eq!(sgroup.multiplier, None);
}

#[rstest]
fn test_apply_extended_sgroup_multiplier(
    mut single_extended_atom: ExtendedMolecule,
    mut acc_with_multiple_sgroup: PropertyAccumulator,
    flags_extended: CtabParseFlags,
) {
    let entry = PropertyEntries::SGroupSubscriptEntry(SGroupSubscriptEntry {
        sgroup_index: 0,
        multiplier: Some(SGroupMultiplier::Single(SGroupMultiplierTerm::Variable(
            'n',
        ))),
        subscript: None,
    });
    acc_with_multiple_sgroup
        .add_entry(entry, flags_extended)
        .unwrap();
    acc_with_multiple_sgroup
        .update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();

    let sgroup = single_extended_atom.sgroups().get(&0).unwrap();
    assert_eq!(sgroup.subscript, None);
    assert_eq!(
        sgroup.multiplier,
        Some(SGroupMultiplier::Single(SGroupMultiplierTerm::Variable(
            'n'
        )))
    );
}

#[rstest]
fn test_apply_extended_sgroup_correspondence(
    mut single_extended_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: PropertyAccumulator,
    flags_extended: CtabParseFlags,
) {
    let entry = PropertyEntries::SGroupCorrespondenceEntry(SGroupCorrespondenceEntry {
        sgroup_index: 0,
        bond_indices: vec![0, 1],
    });
    acc_with_superatom_sgroup
        .add_entry(entry, flags_extended)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();

    let sgroup = single_extended_atom.sgroups().get(&0).unwrap();
    assert_eq!(sgroup.correspondence, Some(vec![0, 1]));
}

#[rstest]
fn test_apply_extended_sgroup_display_info(
    mut single_extended_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: PropertyAccumulator,
    flags_extended: CtabParseFlags,
) {
    let entry = PropertyEntries::SGroupDisplayInfoEntry(SGroupDisplayInfoEntry {
        sgroup_index: 0,
        bracket_coords: vec![1.0, 2.0, 3.0, 4.0],
    });
    acc_with_superatom_sgroup
        .add_entry(entry, flags_extended)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();

    let sgroup = single_extended_atom.sgroups().get(&0).unwrap();
    assert_eq!(
        sgroup.bracket_coords,
        Some(SGroupBracketCoords {
            bracket1: (1.0, 2.0),
            bracket2: (3.0, 4.0)
        })
    );
}

#[rstest]
fn test_apply_extended_sgroup_connecting_bond(
    mut single_extended_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: PropertyAccumulator,
    flags_extended: CtabParseFlags,
) {
    let entry = PropertyEntries::SGroupConnectingBondEntry(SGroupConnectingBondEntry {
        sgroup_index: 0,
        bond_index: 0,
        bond_vector: (1.0, 2.0),
    });
    acc_with_superatom_sgroup
        .add_entry(entry, flags_extended)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();

    let sgroup = single_extended_atom.sgroups().get(&0).unwrap();
    assert_eq!(
        sgroup.connecting_bond,
        Some(SGroupConnectingBond {
            bond_index: 0,
            bond_vector: (1.0, 2.0)
        })
    );
}

#[rstest]
fn test_apply_extended_sgroup_data_description(
    mut single_extended_atom: ExtendedMolecule,
    mut acc_with_data_sgroup: PropertyAccumulator,
    flags_extended: CtabParseFlags,
) {
    acc_with_data_sgroup
        .update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();

    let sgroup = single_extended_atom.sgroups().get(&0).unwrap();
    let data = sgroup.data.get("test").unwrap();
    assert_eq!(data.field_type, SGroupDataType::Text);
    assert_eq!(data.field_units, None);
    assert_eq!(data.query_identifier, None);
    assert_eq!(data.data_query_operator, None);
}

#[rstest]
fn test_apply_extended_sgroup_data_entry(
    mut single_extended_atom: ExtendedMolecule,
    mut acc_with_data_sgroup: PropertyAccumulator,
    flags_extended: CtabParseFlags,
) {
    let continuation_entry = PropertyEntries::SGroupDataEntry(SGroupDataEntry::Continuation {
        sgroup_index: 0,
        data_content: "content".to_string(),
    });
    acc_with_data_sgroup
        .add_entry(continuation_entry, flags_extended)
        .unwrap();

    let data_entry =
        PropertyEntries::SGroupDataEntry(SGroupDataEntry::EndBlank { sgroup_index: 0 });
    acc_with_data_sgroup
        .add_entry(data_entry, flags_extended)
        .unwrap();
    acc_with_data_sgroup
        .update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();

    let sgroup = single_extended_atom.sgroups().get(&0).unwrap();
    let data = sgroup.data.get("test").unwrap();
    assert_eq!(data.data_content, Some(vec!["content".to_string()]));
}

#[rstest]
fn test_apply_extended_sgroup_data_display(
    mut single_extended_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: PropertyAccumulator,
    flags_extended: CtabParseFlags,
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
        .add_entry(display_entry, flags_extended)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();

    let sgroup = single_extended_atom.sgroups().get(&0).unwrap();
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
fn test_apply_extended_sgroup_hierarchy(
    mut single_extended_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: PropertyAccumulator,
    flags_extended: CtabParseFlags,
) {
    // Add a second SGroup as parent
    let parent_type_entry = PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry {
        sgroup_index: 1,
        sgroup_type: SGroupType::Superatom,
    }]);
    acc_with_superatom_sgroup
        .add_entry(parent_type_entry, flags_extended)
        .unwrap();

    let hierarchy_entry = PropertyEntries::SGroupHierarchyEntries(vec![SGroupHierarchyEntry {
        sgroup_index: 0,
        parent_sgroup_index: 1,
    }]);
    acc_with_superatom_sgroup
        .add_entry(hierarchy_entry, flags_extended)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();

    assert_eq!(single_extended_atom.sgroups().len(), 2);
    let child_sgroup = single_extended_atom.sgroups().get(&0).unwrap();
    assert_eq!(child_sgroup.hierarchy_parent, Some(1));
}

#[rstest]
fn test_apply_extended_sgroup_component(
    mut single_extended_atom: ExtendedMolecule,
    mut acc_with_superatom_sgroup: PropertyAccumulator,
    flags_extended: CtabParseFlags,
) {
    let component_entry = PropertyEntries::SGroupComponentEntries(vec![SGroupComponentEntry {
        sgroup_index: 0,
        component_number: 12,
    }]);
    acc_with_superatom_sgroup
        .add_entry(component_entry, flags_extended)
        .unwrap();
    acc_with_superatom_sgroup
        .update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();

    let sgroup = single_extended_atom.sgroups().get(&0).unwrap();
    assert_eq!(sgroup.component_number, Some(12));
}

#[rstest]
#[case::none(0, None)]
#[case::r2(2, Some(RingBondCount::R2))]
#[case::r4_plus(4, Some(RingBondCount::R4Plus))]
#[case::no_ring_bonds(-1, Some(RingBondCount::NoRingBonds))]
fn test_apply_extended_ring_bond_count(
    mut single_extended_atom: ExtendedMolecule,
    #[case] code: i8,
    #[case] expected: Option<RingBondCount>,
    flags_extended: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::RingBondCountEntries(vec![RingBondCountEntry {
        atom_index: 0,
        ring_bond_count: code,
    }]);
    acc.add_entry(entry, flags_extended).unwrap();
    acc.update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();
    assert_eq!(single_extended_atom.atoms[0].ring_bond_count, expected);
}

#[rstest]
fn test_apply_extended_ring_bond_count_conflict(
    mut with_extended_properties: ExtendedMolecule,
    flags_extended: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::RingBondCountEntries(vec![RingBondCountEntry {
        atom_index: 0,
        ring_bond_count: 3,
    }]);
    acc.add_entry(entry, flags_extended).unwrap();
    let result = acc.update_extended_molecule(&mut with_extended_properties, flags_extended);
    assert!(result.is_err());
}

#[rstest]
#[case::out_of_range_low(1)]
#[case::out_of_range_high(5)]
fn test_apply_extended_ring_bond_count_invalid(#[case] code: i8, flags_strict: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::RingBondCountEntries(vec![RingBondCountEntry {
        atom_index: 0,
        ring_bond_count: code,
    }]);
    let result = acc.add_entry(entry, flags_strict);
    assert!(result.is_err());
}

#[rstest]
#[case::none(0, None)]
#[case::s2(2, Some(SubstitutionCount::S2))]
#[case::s4(4, Some(SubstitutionCount::S4))]
#[case::no_substitution(-1, Some(SubstitutionCount::NoSubstitution))]
fn test_apply_extended_substitution_count(
    mut single_extended_atom: ExtendedMolecule,
    #[case] code: i8,
    #[case] expected: Option<SubstitutionCount>,
    flags_extended: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::SubstitutionCountEntries(vec![SubstitutionCountEntry {
        atom_index: 0,
        substitution_count: code,
    }]);
    acc.add_entry(entry, flags_extended).unwrap();
    acc.update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();
    assert_eq!(single_extended_atom.atoms[0].substitution_count, expected);
}

#[rstest]
fn test_apply_extended_substitution_count_conflict(
    mut with_extended_properties: ExtendedMolecule,
    flags_extended: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::SubstitutionCountEntries(vec![SubstitutionCountEntry {
        atom_index: 0,
        substitution_count: 3,
    }]);
    acc.add_entry(entry, flags_extended).unwrap();
    let result = acc.update_extended_molecule(&mut with_extended_properties, flags_extended);
    assert!(result.is_err());
}

#[rstest]
#[case::out_of_range_low(-3)]
#[case::out_of_range_high(7)]
fn test_apply_extended_substitution_count_invalid(#[case] code: i8, flags_strict: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::SubstitutionCountEntries(vec![SubstitutionCountEntry {
        atom_index: 0,
        substitution_count: code,
    }]);
    let result = acc.add_entry(entry, flags_strict);
    assert!(result.is_err());
}

#[rstest]
#[case::none(0, None)]
#[case::unsaturated(1, Some(UnsaturatedAtom))]
fn test_apply_extended_unsaturated(
    mut single_extended_atom: ExtendedMolecule,
    #[case] code: u8,
    #[case] expected: Option<UnsaturatedAtom>,
    flags_extended: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::UnsaturatedAtomEntries(vec![UnsaturatedAtomEntry {
        atom_index: 0,
        unsaturated: code,
    }]);
    acc.add_entry(entry, flags_extended).unwrap();
    acc.update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();
    assert_eq!(single_extended_atom.atoms[0].unsaturated, expected);
}

#[rstest]
#[case::out_of_range_high(2)]
fn test_apply_extended_unsaturated_invalid(#[case] code: u8, flags_strict: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::UnsaturatedAtomEntries(vec![UnsaturatedAtomEntry {
        atom_index: 0,
        unsaturated: code,
    }]);
    let result = acc.add_entry(entry, flags_strict);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_extended_link_atom(
    mut single_extended_atom: ExtendedMolecule,
    flags_extended: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::LinkAtomEntries(vec![LinkAtomEntry {
        atom_index: 0,
        repeat_count: 2,
        subs_index1: 1,
        subs_index2: None,
    }]);
    acc.add_entry(entry, flags_extended).unwrap();
    acc.update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();
    assert_eq!(
        single_extended_atom.atoms[0].link_atom,
        Some(LinkAtom {
            repeat_count: 2,
            subs_index1: 1,
            subs_index2: None
        })
    );
}

#[rstest]
fn test_apply_extended_link_atom_conflict(flags_extended: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
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
    let result = acc.add_entry(entry, flags_extended);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_extended_atom_list(
    mut single_extended_atom: ExtendedMolecule,
    flags_extended: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::AtomListEntry(AtomListEntry {
        atom_index: 0,
        elements: vec![e!(N), e!(O)],
        exclusion: false,
    });
    acc.add_entry(entry, flags_extended).unwrap();
    acc.update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();
    assert_eq!(
        single_extended_atom.atoms[0].symbol,
        AtomSymbol::AtomList(AtomList {
            elements: vec![e!(N), e!(O)],
            exclusion: false
        })
    );
}

#[rstest]
fn test_apply_extended_atom_list_conflict(
    mut with_rgroup: ExtendedMolecule,
    flags_extended: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::AtomListEntry(AtomListEntry {
        atom_index: 0,
        elements: vec![e!(N), e!(O)],
        exclusion: false,
    });
    acc.add_entry(entry, flags_extended).unwrap();
    let result = acc.update_extended_molecule(&mut with_rgroup, flags_extended);
    assert!(result.is_err());
}

#[rstest]
#[case::first(1, Some(AttachmentPointType::First))]
#[case::second(2, Some(AttachmentPointType::Second))]
#[case::both(3, Some(AttachmentPointType::Both))]
#[case::none(0, None)]
fn test_apply_extended_attachment_point(
    mut single_extended_atom: ExtendedMolecule,
    #[case] code: u8,
    #[case] expected: Option<AttachmentPointType>,
    flags_extended: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::AttachmentPointEntries(vec![AttachmentPointEntry {
        atom_index: 0,
        attachment_type: code,
    }]);
    acc.add_entry(entry, flags_extended).unwrap();
    acc.update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();
    assert_eq!(single_extended_atom.atoms[0].attachment_point, expected);
}

#[rstest]
fn test_apply_extended_attachment_point_conflict(flags_extended: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
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
    let result = acc.add_entry(entry, flags_extended);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_extended_attachment_point_invalid(flags_strict: CtabParseFlags) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::AttachmentPointEntries(vec![AttachmentPointEntry {
        atom_index: 0,
        attachment_type: 4,
    }]);
    let result = acc.add_entry(entry, flags_strict);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_extended_attachment_order(
    mut single_extended_atom: ExtendedMolecule,
    flags_extended: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::AtomAttachmentOrderEntry(AtomAttachmentOrderEntry {
        atom_index: 0,
        attachments: vec![(1, 2), (2, 1)],
    });
    acc.add_entry(entry, flags_extended).unwrap();
    acc.update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();
    assert_eq!(
        single_extended_atom.atoms[0]
            .attachment_order
            .as_ref()
            .unwrap(),
        &vec![(1, 2), (2, 1)]
    );
}

#[rstest]
fn test_apply_extended_rgroup_label(
    mut single_extended_atom: ExtendedMolecule,
    flags_extended: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::RGroupLabelEntries(vec![RGroupLabelEntry {
        atom_index: 0,
        label: 1,
    }]);
    acc.add_entry(entry, flags_extended).unwrap();
    acc.update_extended_molecule(&mut single_extended_atom, flags_extended)
        .unwrap();
    // Accumulator creates minimal RGroup with empty occurrence
    if let AtomSymbol::RGroup(rgroup) = &single_extended_atom.atoms[0].symbol {
        assert_eq!(rgroup.label, Some(1));
    } else {
        panic!("Expected RGroup symbol");
    }
}

#[rstest]
fn test_apply_extended_rgroup_label_keep(
    mut with_rgroup: ExtendedMolecule,
    flags_extended: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::RGroupLabelEntries(vec![RGroupLabelEntry {
        atom_index: 0,
        label: 1,
    }]);
    acc.add_entry(entry, flags_extended).unwrap();
    acc.update_extended_molecule(&mut with_rgroup, flags_extended)
        .unwrap();

    // Keeps original label (doesn't overwrite same label)
    if let AtomSymbol::RGroup(rgroup) = &with_rgroup.atoms[0].symbol {
        assert_eq!(rgroup.label, Some(1));
    } else {
        panic!("Expected RGroup symbol");
    }
}

#[rstest]
fn test_apply_extended_rgroup_label_overwrite(
    mut with_unlabeled_rgroup: ExtendedMolecule,
    flags_extended: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::RGroupLabelEntries(vec![RGroupLabelEntry {
        atom_index: 0,
        label: 3,
    }]);
    acc.add_entry(entry, flags_extended).unwrap();
    acc.update_extended_molecule(&mut with_unlabeled_rgroup, flags_extended)
        .unwrap();

    // Label is overwritten
    if let AtomSymbol::RGroup(rgroup) = &with_unlabeled_rgroup.atoms[0].symbol {
        assert_eq!(rgroup.label, Some(3));
    } else {
        panic!("Expected RGroup symbol");
    }
}

#[rstest]
fn test_apply_extended_rgroup_label_conflict(
    mut with_rgroup: ExtendedMolecule,
    flags_extended: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::RGroupLabelEntries(vec![RGroupLabelEntry {
        atom_index: 0,
        label: 2,
    }]);
    acc.add_entry(entry, flags_extended).unwrap();
    let result = acc.update_extended_molecule(&mut with_rgroup, flags_extended);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_extended_rgroup_label_invalid(
    mut single_extended_atom: ExtendedMolecule,
    flags_strict: CtabParseFlags,
) {
    single_extended_atom.atoms[0].symbol = AtomSymbol::AtomList(AtomList {
        elements: vec![e!(N), e!(O)],
        exclusion: false,
    });
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::RGroupLabelEntries(vec![RGroupLabelEntry {
        atom_index: 0,
        label: 1,
    }]);
    acc.add_entry(entry, flags_strict).unwrap();
    let result = acc.update_extended_molecule(&mut single_extended_atom, flags_strict);
    assert!(result.is_err());
}

#[rstest]
fn test_apply_extended_rgroup_logic(
    mut with_rgroup: ExtendedMolecule,
    flags_extended: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::RGroupLogicEntry(RGroupLogicEntry {
        label: 1,
        dependent_label: Some(2),
        rgroup_or_h: true,
        occurrence: vec![RGroupOccurrence::Exactly(1)],
    });
    acc.add_entry(entry, flags_extended).unwrap();
    acc.update_extended_molecule(&mut with_rgroup, flags_extended)
        .unwrap();

    let rgroup = with_rgroup.rgroups().get(&1).unwrap();
    assert_eq!(rgroup.dependent_label, Some(2));
    assert!(rgroup.rgroup_or_h);
    assert_eq!(rgroup.occurrence.len(), 1);
    assert_eq!(rgroup.occurrence[0], RGroupOccurrence::Exactly(1));
}

#[rstest]
fn test_apply_extended_rgroup_logic_multiple_occurrences(
    mut with_rgroup: ExtendedMolecule,
    flags_extended: CtabParseFlags,
) {
    let mut acc = PropertyAccumulator::new();
    let entry = PropertyEntries::RGroupLogicEntry(RGroupLogicEntry {
        label: 1,
        dependent_label: Some(2),
        rgroup_or_h: true,
        occurrence: vec![
            RGroupOccurrence::Exactly(1),
            RGroupOccurrence::GreaterThan(5),
        ],
    });
    acc.add_entry(entry, flags_extended).unwrap();
    acc.update_extended_molecule(&mut with_rgroup, flags_extended)
        .unwrap();

    let rgroup = with_rgroup.rgroups().get(&1).unwrap();
    assert_eq!(rgroup.occurrence.len(), 2);
    assert_eq!(rgroup.occurrence[0], RGroupOccurrence::Exactly(1));
    assert_eq!(rgroup.occurrence[1], RGroupOccurrence::GreaterThan(5));
}
