use bstr::ByteSlice;
use float_cmp::*;
use nom::error::ErrorKind as NomErrorKind;
use nom::Err;
use pretty_assertions::assert_eq;
use rstest::*;
use umol_data::Element;

use super::*;
use crate::io::ctfile::config::CtabParseFlags;
use crate::table_ir::{SGroupMultiplierOp, SGroupMultiplierTerm};

#[test]
fn test_properties_block() {
    let ctab_data = b"M  CHG  1   2  -1\nM  END";
    let flags = CtabParseFlags::BASIC;
    let result = properties_block(0, flags).parse(ctab_data);
    assert!(
        result.is_ok(),
        "Properties block should parse successfully in BASIC mode"
    );

    let (remaining, (property_entries, _)) = result.unwrap();
    assert!(remaining.is_empty(), "All input should be consumed");
    assert_eq!(property_entries.len(), 1, "Should have 1 property entry");
    assert!(matches!(
        property_entries[0],
        PropertyEntries::ChargeEntries(_)
    ));
}

#[test]
fn test_extended_properties_block() {
    let ctab_data = b"M  ALS   1  2 F Cl  Br\nM  END";
    let flags = CtabParseFlags::EXTENDED;
    let result = extended_properties_block(0, flags).parse(ctab_data);
    assert!(
        result.is_ok(),
        "Extended properties block should parse successfully in EXTENDED mode"
    );
    let (remaining, _) = result.unwrap();
    assert!(remaining.is_empty(), "All input should be consumed");
}

#[test]
fn test_properties_block_m_end_only() {
    let ctab_data = b"M  END";
    let flags = CtabParseFlags::BASIC;
    let result = properties_block(0, flags).parse(ctab_data);
    assert!(
        result.is_ok(),
        "Properties block should parse with just M  END"
    );

    let (remaining, (property_entries, _)) = result.unwrap();
    assert!(remaining.is_empty(), "All input should be consumed");
    assert_eq!(property_entries.len(), 0, "Should have 0 property entries");
}

#[rstest]
#[case::v_atom_value(b"V    1 *", PropertyEntries::AtomValueEntry(AtomValueEntry { atom_index: 0, value: "*".to_string() }))]
#[case::chg_atom(b"M  CHG  1   1  -1", PropertyEntries::ChargeEntries(vec![ChargeEntry { atom_index: 0, charge: -1 }]))]
#[case::rad_atom(b"M  RAD  1   1   2", PropertyEntries::RadicalEntries(vec![RadicalEntry { atom_index: 0, radical_type: 2 }]))]
#[case::iso_atom(b"M  ISO  1   1  13", PropertyEntries::IsotopeEntries(vec![IsotopeEntry { atom_index: 0, mass: 13 }]))]
fn test_property_input(#[case] input: &[u8], #[case] expected: PropertyEntries) {
    let (remaining, result) = all_consuming(property_input(CtabParseFlags::BASIC))
        .parse(input)
        .unwrap();
    let input_str = input.to_str_lossy();
    assert!(
        remaining.is_empty(),
        "remaining should be empty for {:?}",
        input_str
    );
    assert_eq!(result, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::rbc_query(b"M  RBC  1   1   2")]
#[case::sub_query(b"M  SUB  1   1   3")]
#[case::uns_query(b"M  UNS  1   1   1")]
#[case::lin_query(b"M  LIN  1   1   2   5   7")]
#[case::als_query(b"M  ALS  1  3FC   N   O   ")]
#[case::apo_query(b"M  APO  1   1   1")]
#[case::aal_query(b"M  AAL  1 1   2   1")]
#[case::rgp_rgroup(b"M  RGP   1   1   2")]
#[case::log_rgroup(b"M  LOG   1   1   0   0  >2")]
#[case::sty_sgroup(b"M  STY  1   1 SUP")]
#[case::sst_sgroup(b"M  SST  1   1 ALT")]
#[case::slb_sgroup(b"M  SLB  1   1  19")]
#[case::scn_sgroup(b"M  SCN  1   1 HH ")]
#[case::sal_sgroup(b"M  SAL   1  1   5")]
#[case::sbl_sgroup(b"M  SBL   1  1   3")]
#[case::smt_sgroup(b"M  SMT   1 n")]
#[case::sds_sgroup(b"M  SDS EXP  1   1")]
#[case::crs_sgroup(b"M  CRS   1  3  10   9   4")]
#[case::sdi_sgroup(b"M  SDI   3  4    4.4700   -3.1700    4.4700   -5.7500")]
#[case::sbv_sgroup(b"M  SBV   1  11    0.6400    0.9700")]
#[case::sdt_sgroup(b"M  SDT   1 pH   ")]
#[case::sdd_sgroup(b"M  SDD   1     0.0000    0.0000    DR    ALL  1       6")]
#[case::scd_sgroup(b"M  SCD   1   1   0")]
#[case::sed_sgroup(b"M  SED   1   1   0")]
#[case::spl_sgroup(b"M  SPL   1   1   0")]
#[case::snc_sgroup(b"M  SNC   1   1   0")]
#[case::sty_sgroup(b"M  STY  1   1 SUP")]
#[case::sst_sgroup(b"M  SST  1   1 ALT")]
#[case::slb_sgroup(b"M  SLB  1   1  19")]
#[case::sal_sgroup(b"M  SAL   1  1   5")]
#[case::sbl_sgroup(b"M  SBL   1  1   3")]
#[case::smt_sgroup(b"M  SMT   1 n")]
#[case::zbo_clark_extensions(b"M  ZBO  1   1   0")]
#[case::zch_clark_extensions(b"M  ZCH  1   1  -1")]
#[case::hyd_clark_extensions(b"M  HYD  1   1   1")]
fn test_property_input_invalid(#[case] input: &[u8]) {
    let result = all_consuming(property_input(CtabParseFlags::BASIC)).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?}", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == NomErrorKind::Tag),
        "Expected Tag error for {:?}, got {:?}",
        input_str,
        result
    );
}

#[rustfmt::skip]
#[rstest]
#[case::bond_order_override(b"M  ZBO  1   1   0", PropertyEntries::BondOrderOverrideEntries(vec![BondOrderOverrideEntry { bond_index: 0, bond_order: BondOrder::Zero }]))]
#[case::atom_charge_override(b"M  ZCH  1   1  -1", PropertyEntries::AtomChargeOverrideEntries(vec![AtomChargeOverrideEntry { atom_index: 0, charge: -1 }]))]
#[case::atom_hydrogen_count(b"M  HYD  1   1   1", PropertyEntries::AtomHydrogenCountEntries(vec![AtomHydrogenCountEntry { atom_index: 0, hydrogen_count: Some(1) }]))]
fn test_property_input_lenient(#[case] input: &[u8], #[case] expected: PropertyEntries) {
    let (remaining, result) = all_consuming(property_input(CtabParseFlags::LENIENT))
        .parse(input)
        .unwrap();
    let input_str = input.to_str_lossy();
    assert!(
        remaining.is_empty(),
        "remaining should be empty for {:?}",
        input_str
    );
    assert_eq!(result, expected);
}

#[rstest]
#[case::v_atom_value(b"V    1 *", PropertyEntries::AtomValueEntry(AtomValueEntry { atom_index: 0, value: "*".to_string() }))]
#[case::chg(b"M  CHG  1   1  -1", PropertyEntries::ChargeEntries(vec![ChargeEntry { atom_index: 0, charge: -1 }]))]
#[case::rad(b"M  RAD  1   1   2", PropertyEntries::RadicalEntries(vec![RadicalEntry { atom_index: 0, radical_type: 2 }]))]
#[case::iso(b"M  ISO  1   1  13", PropertyEntries::IsotopeEntries(vec![IsotopeEntry { atom_index: 0, mass: 13 }]))]
#[case::rbc_query(b"M  RBC  1   1   2", PropertyEntries::RingBondCountEntries(vec![RingBondCountEntry { atom_index: 0, ring_bond_count: 2 }]))]
#[case::sub_query(b"M  SUB  1   1   3", PropertyEntries::SubstitutionCountEntries(vec![SubstitutionCountEntry { atom_index: 0, substitution_count: 3 }]))]
#[case::uns_query(b"M  UNS  1   1   1", PropertyEntries::UnsaturatedAtomEntries(vec![UnsaturatedAtomEntry { atom_index: 0, unsaturated: 1 }]))]
#[case::lin_query(b"M  LIN  1   1   2   5   7", PropertyEntries::LinkAtomEntries(vec![LinkAtomEntry { atom_index: 0, repeat_count: 2, subs_index1: 4, subs_index2: Some(6) }]))]
#[case::als_normal_f_query(b"M  ALS   1  3   F   N   O   ", PropertyEntries::AtomListEntry(AtomListEntry { atom_index: 0, exclusion: false, elements: vec![Element::F, Element::N, Element::O] }))]
#[case::als_exclusion_t_query(b"M  ALS   1  3 T F   N   O   ", PropertyEntries::AtomListEntry(AtomListEntry { atom_index: 0, exclusion: true, elements: vec![Element::F, Element::N, Element::O] }))]
#[case::apo_query(b"M  APO  1   1   1", PropertyEntries::AttachmentPointEntries(vec![AttachmentPointEntry { atom_index: 0, attachment_type: 1 }]))]
#[case::aal_query(b"M  AAL   4  2  14   1   9   2", PropertyEntries::AtomAttachmentOrderEntry(AtomAttachmentOrderEntry {
        atom_index: 3, attachments: vec![(13, 1), (8, 2)] }))]
#[case::rgp_rgroup(b"M  RGP  1   1   2", PropertyEntries::RGroupLabelEntries(vec![RGroupLabelEntry { atom_index: 0, label: 2 }]))]
#[case::log_rgroup(b"M  LOG  1   1   0   0  >2", PropertyEntries::RGroupLogicEntry(RGroupLogicEntry { label: 1, dependent_label: None, rgroup_or_h: false, occurrence: vec![RGroupOccurrence::GreaterThan(2)] }))]
#[case::sty_sgroup(b"M  STY  1   1 SUP", PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry { sgroup_index: 0, sgroup_type: SGroupType::Superatom }]))]
#[case::sst_sgroup(b"M  SST  1   1 ALT", PropertyEntries::SGroupSubtypeEntries(vec![SGroupSubtypeEntry { sgroup_index: 0, sgroup_subtype: SGroupSubtype::Alternating }]))]
#[case::slb_sgroup(b"M  SLB  1   1  19", PropertyEntries::SGroupLabelEntries(vec![SGroupLabelEntry { sgroup_index: 0, label: 19 }]))]
#[case::scn_sgroup(b"M  SCN  1   1 HH ", PropertyEntries::SGroupConnectivityEntries(vec![SGroupConnectivityEntry { sgroup_index: 0, connectivity: SGroupConnectivity::HeadToHead }]))]
#[case::sds_sgroup(b"M  SDS EXP  1   1", PropertyEntries::SGroupExpansionEntries(vec![SGroupExpansionEntry { sgroup_index: 0 }]))]
#[case::sal_sgroup(b"M  SAL   1  1   5", PropertyEntries::SGroupAtomListEntry(SGroupAtomListEntry { sgroup_index: 0, atom_indices: vec![4] }))]
#[case::sbl_sgroup(b"M  SBL   1  1   3", PropertyEntries::SGroupBondListEntry(SGroupBondListEntry { sgroup_index: 0, bond_indices: vec![2] }))]
#[case::spa_sgroup(b"M  SPA   1 12   3   4   5   6   9  10  11  12  13  14  15  16",
       PropertyEntries::SGroupParentAtomEntry(SGroupParentAtomEntry { sgroup_index: 0, atom_indices: vec![2, 3, 4, 5, 8, 9, 10, 11, 12, 13, 14, 15] }))]
#[case::smt_sgroup(b"M  SMT   1 n", PropertyEntries::SGroupSubscriptEntry(SGroupSubscriptEntry { sgroup_index: 0,
       multiplier: Some(SGroupMultiplier::Single(SGroupMultiplierTerm::Variable('n'))), subscript: Some("n".to_string()) }))]
#[case::crs_sgroup(b"M  CRS   1  3  10   9   4", PropertyEntries::SGroupCorrespondenceEntry(SGroupCorrespondenceEntry { sgroup_index: 0, bond_indices: vec![9, 8, 3] }))]
#[case::sdi_sgroup(b"M  SDI   3  4    4.4700   -3.1700    4.4700   -5.7500", PropertyEntries::SGroupDisplayInfoEntry(SGroupDisplayInfoEntry { sgroup_index: 2, bracket_coords: vec![4.4700, -3.1700, 4.4700, -5.7500] }))]
#[case::sbv_sgroup(b"M  SBV   1  11    0.6400    0.9700", PropertyEntries::SGroupConnectingBondEntry(SGroupConnectingBondEntry { sgroup_index: 0, bond_index: 10, bond_vector: (0.6400, 0.9700) }))]
#[case::sdt_sgroup(b"M  SDT   1 pH   ", PropertyEntries::SGroupDataDescriptionEntry(SGroupDataDescriptionEntry { sgroup_index: 0, field_name: "pH".to_string(),
    field_type: SGroupDataType::Text, field_units: None, query_identifier: None, data_query_operator: None }))]
#[case::sdd_sgroup(b"M  SDD   1     0.0000    0.0000    DR    ALL  1       6", PropertyEntries::SGroupDataDisplayEntry(SGroupDataDisplayEntry { sgroup_index: 0,
    coords: (0.0000, 0.0000), display_type: SGroupDataDisplayType::Detached, display_placement: SGroupDataDisplayPlacement::Relative, display_units: SGroupDataDisplayUnits::None,
    display_chars: SGroupDataDisplayChars::All, display_tag: None, display_position: 6 }))]
#[case::scd_sgroup(b"M  SCD   1 4.6", PropertyEntries::SGroupDataEntry(SGroupDataEntry::Continuation { sgroup_index: 0, data_content: "4.6".to_string() }))]
#[case::sed_sgroup(b"M  SED   2 E/Z unknown", PropertyEntries::SGroupDataEntry(SGroupDataEntry::EndWithData { sgroup_index: 1, data_content: "E/Z unknown".to_string() }))]
#[case::sed_sgroup_empty(b"M  SED   1", PropertyEntries::SGroupDataEntry(SGroupDataEntry::EndBlank { sgroup_index: 0 }))]
#[case::spl_sgroup(b"M  SPL  1   1   2", PropertyEntries::SGroupHierarchyEntries(vec![SGroupHierarchyEntry { sgroup_index: 0, parent_sgroup_index: 1 }]))]
#[case::snc_sgroup(b"M  SNC  2   1   1   2   2", PropertyEntries::SGroupComponentEntries(vec![SGroupComponentEntry { sgroup_index: 0, component_number: 1 }, SGroupComponentEntry { sgroup_index: 1, component_number: 2 }]))]
fn test_extended_property_input(#[case] input: &[u8], #[case] expected: PropertyEntries) {
    let result = all_consuming(extended_property_input(CtabParseFlags::EXTENDED)).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "remaining should be empty for {:?}",
        input_str
    );
    assert_eq!(result, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::trailing_chars(b"M  CHG  1   1  -1  a", NomErrorKind::Eof)]
#[case::count_mismatch(b"M  CHG  2   1  -1", NomErrorKind::Tag)]
#[case::count_zero(b"M  CHG  0", NomErrorKind::Verify)]
#[case::atom_index_zero(b"M  CHG  1   0 -10", NomErrorKind::Verify)]
#[case::invalid_property_tag(b"M  XXX  1   1  -1", NomErrorKind::Tag)]
#[case::bond_order_override(b"M  ZBO  1   1   0", NomErrorKind::Tag)]
#[case::atom_charge_override(b"M  ZCH  1   1  -1", NomErrorKind::Tag)]
#[case::atom_hydrogen_count(b"M  HYD  1   1   1", NomErrorKind::Tag)]
fn test_extended_property_input_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let result = all_consuming(extended_property_input(CtabParseFlags::EXTENDED)).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.as_ref(), Err(Err::Error(e)) if e.code == expected_kind),
        "Expected {:?} error for {:?}, got {:?}",
        expected_kind,
        input_str,
        result
    );
}

#[rustfmt::skip]
#[rstest]
#[case::bond_order_override(b"M  ZBO  1   1   0", PropertyEntries::BondOrderOverrideEntries(vec![BondOrderOverrideEntry { bond_index: 0, bond_order: BondOrder::Zero }]))]
#[case::atom_charge_override(b"M  ZCH  1   1  -1", PropertyEntries::AtomChargeOverrideEntries(vec![AtomChargeOverrideEntry { atom_index: 0, charge: -1 }]))]
#[case::atom_hydrogen_count(b"M  HYD  1   1   1", PropertyEntries::AtomHydrogenCountEntries(vec![AtomHydrogenCountEntry { atom_index: 0, hydrogen_count: Some(1) }]))]
fn test_extended_property_input_lenient(#[case] input: &[u8], #[case] expected: PropertyEntries) {
    let result = all_consuming(extended_property_input(CtabParseFlags::LENIENT)).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "remaining should be empty for {:?}",
        input_str
    );
    assert_eq!(result, expected);
}

#[rstest]
#[case::no_space(b"A    1", b"CF3", AtomAliasEntry { atom_index: 0, alias: "CF3".to_string() })]
#[case::leading_space(b"A   15", b"  Et", AtomAliasEntry { atom_index: 14, alias: "  Et".to_string() })]
fn test_parse_atom_alias_input(
    #[case] first_line: &[u8],
    #[case] second_line: &[u8],
    #[case] expected: AtomAliasEntry,
) {
    let result = parse_atom_alias_input(first_line, second_line);
    assert!(result.is_ok(), "{:?} should have succeeded", first_line);
    assert_eq!(result.unwrap(), PropertyEntries::AtomAliasEntry(expected));
}

#[rstest]
#[case::atom_index_is_zero(b"A    0", b"Et", NomErrorKind::Digit)]
#[case::too_short(b"A   1", b"Et", NomErrorKind::Digit)]
fn test_parse_atom_alias_input_invalid(
    #[case] first_line: &[u8],
    #[case] second_line: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = parse_atom_alias_input(first_line, second_line);
    assert!(result.is_err(), "{:?} should have failed", first_line);
    assert!(
        matches!(result.as_ref(), Err(e) if e.code == expected_kind),
        "Expected {:?} error, got {:?}",
        expected_kind,
        result
    );
}

#[rstest]
#[case::basic(b"G   12 11", b"SH", LegacyGroupAbbreviationEntry { atom_index1: 11, atom_index2: 10, label: "SH".to_string() })]
#[case::single_digit(b"G    5  1", b"NH2", LegacyGroupAbbreviationEntry { atom_index1: 4, atom_index2: 0, label: "NH2".to_string() })]
#[case::double_digit(b"G   17 16", b"COOH", LegacyGroupAbbreviationEntry { atom_index1: 16, atom_index2: 15, label: "COOH".to_string() })]
#[case::no_padding(b"G  123456", b"X", LegacyGroupAbbreviationEntry { atom_index1: 122, atom_index2: 455, label: "X".to_string() })]
fn test_parse_legacy_group_abbreviation_input(
    #[case] first_line: &[u8],
    #[case] second_line: &[u8],
    #[case] expected: LegacyGroupAbbreviationEntry,
) {
    let result = parse_legacy_group_abbreviation_input(first_line, second_line);
    assert!(result.is_ok(), "{:?} should have succeeded", first_line);
    assert_eq!(
        result.unwrap(),
        PropertyEntries::LegacyGroupAbbreviationEntry(expected)
    );
}

#[rstest]
#[case::from_atom_is_zero(b"G    0  1", b"Et", NomErrorKind::Digit)]
#[case::to_atom_is_zero(b"G    1  0", b"Et", NomErrorKind::Digit)]
#[case::too_short(b"G   12", b"Et", NomErrorKind::Digit)]
fn test_parse_legacy_group_abbreviation_input_invalid(
    #[case] first_line: &[u8],
    #[case] second_line: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = parse_legacy_group_abbreviation_input(first_line, second_line);
    assert!(result.is_err(), "{:?} should have failed", first_line);
    assert!(
        matches!(result.as_ref(), Err(e) if e.code == expected_kind),
        "Expected {:?} error, got {:?}",
        expected_kind,
        result
    );
}

#[rstest]
#[case::asterisk(b"  1 *", AtomValueEntry { atom_index: 0, value: "*".to_string() })]
#[case::text(b" 15 query", AtomValueEntry { atom_index: 14, value: "query".to_string() })]
fn test_atom_value_entry(#[case] input: &[u8], #[case] expected: AtomValueEntry) {
    let result = all_consuming(atom_value_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case::atom_index_is_zero(b"  0 *", NomErrorKind::Verify)]
fn test_atom_value_entry_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let result = all_consuming(atom_value_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?}", input_str);
    assert!(
        matches!(result.as_ref(), Err(Err::Error(e)) if e.code == expected_kind),
        "Expected {:?} error for {:?}, got {:?}",
        expected_kind,
        input_str,
        result
    );
}

#[rustfmt::skip]
#[rstest]
#[case::single_entry(b"  1   1  -1", vec![ChargeEntry { atom_index: 0, charge: -1 }])]
#[case::two_entries(b"  2   1  -1   4   1", vec![ChargeEntry { atom_index: 0, charge: -1 }, ChargeEntry { atom_index: 3, charge: 1 }])]
#[case::max_entries(b"  8   1   1   2   2   3   3   4   4   5   5   6   6   7   7   8   8",
       vec![ChargeEntry { atom_index: 0, charge: 1 }, ChargeEntry { atom_index: 1, charge: 2 },
            ChargeEntry { atom_index: 2, charge: 3 }, ChargeEntry { atom_index: 3, charge: 4 },
            ChargeEntry { atom_index: 4, charge: 5 }, ChargeEntry { atom_index: 5, charge: 6 },
            ChargeEntry { atom_index: 6, charge: 7 }, ChargeEntry { atom_index: 7, charge: 8 }])]
#[case::max_charge(b"  1  25  15", vec![ChargeEntry { atom_index: 24, charge: 15 }])]
fn test_charge_entries(
    #[case] input: &[u8],
    #[case] expected: Vec<ChargeEntry>,
) {
    let result = all_consuming(charge_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::trailing_characters(b"  1   1  -1  a", NomErrorKind::Eof)]
#[case::count_exceeds_item_list_length(b"  2   1  -1", NomErrorKind::Tag)]
#[case::item_list_exceeds_count(b"  2   1  -1   4   1   5   6", NomErrorKind::Eof)]
#[case::count_is_zero(b"  0", NomErrorKind::Verify)]
#[case::atom_index_is_zero(b"  1   0 -10", NomErrorKind::Verify)]
fn test_charge_entries_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let result = all_consuming(charge_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rustfmt::skip]
#[rstest]
#[case::single_entry(b"  1   1   2", vec![RadicalEntry { atom_index: 0, radical_type: 2 }])]
#[case::two_entries(b"  2   1   1   4   3", vec![RadicalEntry { atom_index: 0, radical_type: 1 }, RadicalEntry { atom_index: 3, radical_type: 3 }])]
fn test_radical_entries(
    #[case] input: &[u8],
    #[case] expected: Vec<RadicalEntry>,
) {
    let result = all_consuming(radical_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::value_out_of_range(b"  1   1   4", NomErrorKind::Verify)]
#[case::value_is_negative(b"  1   1  -1", NomErrorKind::Digit)]
#[case::count_exceeds_item_list_length(b"  2   1   1", NomErrorKind::Tag)]
#[case::item_list_exceeds_count(b"  2   1   1   4   1   5   1", NomErrorKind::Eof)]
#[case::trailing_characters(b"  1   1   2 a", NomErrorKind::Eof)]
#[case::count_is_zero(b"  0", NomErrorKind::Verify)]
#[case::atom_index_is_zero(b"  1   0   2", NomErrorKind::Verify)]
fn test_radical_entries_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = all_consuming(radical_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case::single_entry(b"  1   1  13", vec![IsotopeEntry { atom_index: 0, mass: 13 }])]
#[case::two_entries(b"  2  12   2  15  14", vec![IsotopeEntry { atom_index: 11, mass: 2 }, IsotopeEntry { atom_index: 14, mass: 14 }])]
fn test_isotope_entries(#[case] input: &[u8], #[case] expected: Vec<IsotopeEntry>) {
    let result = all_consuming(isotope_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::value_is_negative(b"  1   1  -1", NomErrorKind::Digit)]
#[case::count_exceeds_item_list_length(b"  2   1  10", NomErrorKind::Tag)]
#[case::item_list_exceeds_count(b"  2   1  10   4   1   5   1", NomErrorKind::Eof)]
#[case::trailing_characters(b"  1   1  12 a", NomErrorKind::Eof)]
#[case::count_is_zero(b"  0", NomErrorKind::Verify)]
#[case::atom_index_is_zero(b"  1   0  12", NomErrorKind::Verify)]
fn test_isotope_entries_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = all_consuming(isotope_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case::single_entry(b"  1   1   2", vec![RingBondCountEntry { atom_index: 0, ring_bond_count: 2 }])]
#[case::two_entries(b"  2   1  -1   4  -2", vec![RingBondCountEntry { atom_index: 0, ring_bond_count: -1 }, RingBondCountEntry { atom_index: 3, ring_bond_count: -2 }])]
#[case::zero_value(b"  1   3   0", vec![RingBondCountEntry { atom_index: 2, ring_bond_count: 0 }])]
#[case::max_value(b"  1  10   4", vec![RingBondCountEntry { atom_index: 9, ring_bond_count: 4 }])]
fn test_ring_bond_count_entries(#[case] input: &[u8], #[case] expected: Vec<RingBondCountEntry>) {
    let result = all_consuming(ring_bond_count_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case::value_out_of_range(b"  1   1   5", NomErrorKind::Verify)]
#[case::value_out_of_range_negative(b"  1   1  -3", NomErrorKind::Verify)]
#[case::trailing_characters(b"  1   1   2 a", NomErrorKind::Eof)]
#[case::count_is_zero(b"  0", NomErrorKind::Verify)]
#[case::atom_index_is_zero(b"  1   0   2", NomErrorKind::Verify)]
fn test_ring_bond_count_entries_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let result = all_consuming(ring_bond_count_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case::single_entry(b"  1   1   3", vec![SubstitutionCountEntry { atom_index: 0, substitution_count: 3 }])]
#[case::two_entries(b"  2   1  -1   4   6", vec![SubstitutionCountEntry { atom_index: 0, substitution_count: -1 }, SubstitutionCountEntry { atom_index: 3, substitution_count: 6 }])]
#[case::negative_value(b"  1   5  -2", vec![SubstitutionCountEntry { atom_index: 4, substitution_count: -2 }])]
fn test_substitution_count_entries(
    #[case] input: &[u8],
    #[case] expected: Vec<SubstitutionCountEntry>,
) {
    let result = all_consuming(substitution_count_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case::value_out_of_range(b"  1   1  16", NomErrorKind::Verify)]
#[case::value_out_of_range_negative(b"  1   1  -3", NomErrorKind::Verify)]
#[case::trailing_characters(b"  1   1   3 a", NomErrorKind::Eof)]
#[case::count_is_zero(b"  0", NomErrorKind::Verify)]
#[case::atom_index_is_zero(b"  1   0   3", NomErrorKind::Verify)]
fn test_substitution_count_entries_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = all_consuming(substitution_count_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case::single_entry(b"  1   1   1", vec![UnsaturatedAtomEntry { atom_index: 0, unsaturated: 1 }])]
#[case::two_entries(b"  2   1   0   3   1", vec![UnsaturatedAtomEntry { atom_index: 0, unsaturated: 0 }, UnsaturatedAtomEntry { atom_index: 2, unsaturated: 1 }])]
#[case::zero_value(b"  1  10   0", vec![UnsaturatedAtomEntry { atom_index: 9, unsaturated: 0 }])]
fn test_unsaturated_atom_entries(
    #[case] input: &[u8],
    #[case] expected: Vec<UnsaturatedAtomEntry>,
) {
    let result = all_consuming(unsaturated_atom_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case::value_out_of_range(b"  1   1   2", NomErrorKind::Verify)]
#[case::unsigned_value_is_negative(b"  1   1  -1", NomErrorKind::Digit)]
#[case::trailing_characters(b"  1   1   1 a", NomErrorKind::Eof)]
#[case::count_is_zero(b"  0", NomErrorKind::Verify)]
#[case::atom_index_is_zero(b"  1   0   1", NomErrorKind::Verify)]
fn test_unsaturated_atom_entries_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = all_consuming(unsaturated_atom_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case::single_entry(b"  1   1   2   5   7", vec![LinkAtomEntry { atom_index: 0, repeat_count: 2, subs_index1: 4, subs_index2: Some(6) }])]
#[case::two_entries(b"  2   3   3   1   3   8   4   5   6",
       vec![LinkAtomEntry { atom_index: 2, repeat_count: 3, subs_index1: 0, subs_index2: Some(2) },
            LinkAtomEntry { atom_index: 7, repeat_count: 4, subs_index1: 4, subs_index2: Some(5) }])]
fn test_link_atom_entries(#[case] input: &[u8], #[case] expected: Vec<LinkAtomEntry>) {
    let result = all_consuming(link_atom_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::repeat_count_less_than_2(b"  1   1   1   5   7", NomErrorKind::Verify)]
#[case::count_exceeds_4(b"  5   1   2   5   7", NomErrorKind::Verify)]
#[case::trailing_characters(b"  1   1   2   5   7 a", NomErrorKind::Eof)]
#[case::count_is_zero(b"  0", NomErrorKind::Verify)]
#[case::atom_index_is_zero(b"  1   0   2   5   7", NomErrorKind::Verify)]
fn test_link_atom_entries_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = all_consuming(link_atom_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {:?}",
        input_str,
        expected_kind,
        result.unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case::exclusion_flag_false(b"   1  3 F C   N   O   ",
       AtomListEntry { atom_index: 0, exclusion: false, elements: vec![Element::C, Element::N, Element::O] })]
#[case::no_right_padding(b"   1  3 F C   N   O",
       AtomListEntry { atom_index: 0, exclusion: false, elements: vec![Element::C, Element::N, Element::O] })]
#[case::exclusion_flag_true(b"   5  2 T Cl  Br  ",
       AtomListEntry { atom_index: 4, exclusion: true, elements: vec![Element::Cl, Element::Br] })]
#[case::no_exclusion_flag(b"  10  1   H   ",
       AtomListEntry { atom_index: 9, exclusion: false, elements: vec![Element::H] })]
fn test_atom_list_entry(#[case] input: &[u8], #[case] expected: AtomListEntry) {
    let result = all_consuming(atom_list_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case::count_is_zero(b"   1  0 F C   ", NomErrorKind::Verify)]
#[case::count_exceeds_16(b"   1 17 F C   N   ", NomErrorKind::Verify)]
#[case::invalid_exclusion_flag(b"   1  1 X C   ", NomErrorKind::Tag)]
#[case::invalid_element_symbol(b"   1  1 F XX  ", NomErrorKind::MapOpt)]
fn test_atom_list_entry_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let result = all_consuming(atom_list_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {:?}",
        input_str,
        expected_kind,
        result.unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case::single_attachment_point(b"  1   1   2", vec![AttachmentPointEntry { atom_index: 0, attachment_type: 2 }])]
#[case::two_attachment_points(b"  2   1   1   2   3",
       vec![AttachmentPointEntry { atom_index: 0, attachment_type: 1 },
            AttachmentPointEntry { atom_index: 1, attachment_type: 3 }])]
fn test_attachment_point_entries(
    #[case] input: &[u8],
    #[case] expected: Vec<AttachmentPointEntry>,
) {
    let result = all_consuming(attachment_point_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::attachment_type_out_of_range(b"  1   1   4", NomErrorKind::Verify)]
#[case::count_is_zero(b"  0   1", NomErrorKind::Verify)]
#[case::atom_index_is_zero(b"  1   0   1", NomErrorKind::Verify)]
#[case::count_exceeds_2(b"  3   1   1   2   2   3   3", NomErrorKind::Verify)]
#[case::trailing_characters(b"  1   1   1 a", NomErrorKind::Eof)]
fn test_attachment_point_entries_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = all_consuming(attachment_point_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {:?}",
        input_str,
        expected_kind,
        result.unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case::two_attachments(b"   4  2  14   1   9   2", AtomAttachmentOrderEntry { atom_index: 3, attachments: vec![(13, 1), (8, 2)] })]
fn test_atom_attachment_order_entry(
    #[case] input: &[u8],
    #[case] expected: AtomAttachmentOrderEntry,
) {
    let result = all_consuming(atom_attachment_order_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::count_is_zero(b"   1   0", NomErrorKind::Verify)]
#[case::count_exceeds_2(b"   1   3   1   2", NomErrorKind::Verify)]
#[case::atom_index_is_zero(b"   0   1   1   2", NomErrorKind::Verify)]
#[case::attachment_type_is_zero(b"   1   1   1   2", NomErrorKind::Verify)]
#[case::attachment_type_out_of_range(b"   1   1   1   3", NomErrorKind::Verify)]
#[case::trailing_characters(b"   1   1   1   2 a", NomErrorKind::Verify)]
fn test_atom_attachment_order_entry_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = all_consuming(atom_attachment_order_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {:?}",
        input_str,
        expected_kind,
        result.unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case::single_entry(b"  1   1   2", vec![RGroupLabelEntry { atom_index: 0, label: 2 }])]
#[case::two_entries(b"  2   1   1   2   2",
       vec![RGroupLabelEntry { atom_index: 0, label: 1 }, RGroupLabelEntry { atom_index: 1, label: 2 }])]
fn test_rgroup_label_entries(#[case] input: &[u8], #[case] expected: Vec<RGroupLabelEntry>) {
    let result = all_consuming(rgroup_label_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case::label_is_zero(b"  1   0", NomErrorKind::Verify)]
#[case::count_exceeds_8(b"  9   1   2", NomErrorKind::Verify)]
#[case::atom_index_is_zero(b"  1   0   2", NomErrorKind::Verify)]
#[case::trailing_characters(b"  1   1   2 a", NomErrorKind::Eof)]
fn test_rgroup_label_entries_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let result = all_consuming(rgroup_label_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?}", input_str);
    assert!(
        matches!(result.as_ref(), Err(Err::Error(e)) if e.code == expected_kind),
        "Expected {:?} error for {:?}, got {:?}",
        expected_kind,
        input_str,
        result
    );
}

#[rstest]
#[case::greater_than(b"  1   1   0   0  >2",
       RGroupLogicEntry { label: 1, dependent_label: None, rgroup_or_h: false, occurrence: vec![RGroupOccurrence::GreaterThan(2)] })]
#[case::exactly_and_greater_than(b"  1   1   0   0  0,>0",
       RGroupLogicEntry { label: 1, dependent_label: None, rgroup_or_h: false, occurrence: vec![RGroupOccurrence::Exactly(0), RGroupOccurrence::GreaterThan(0)] })]
#[case::dependent_label(b"  1   1   2   0",
       RGroupLogicEntry { label: 1, dependent_label: Some(2), rgroup_or_h: false, occurrence: vec![RGroupOccurrence::GreaterThan(0)] })]
#[case::rgroup_or_h(b"  1   1   0   1",
       RGroupLogicEntry { label: 1, dependent_label: None, rgroup_or_h: true, occurrence: vec![RGroupOccurrence::GreaterThan(0)] })]
#[case::no_occurrence(b"  1   1   2",
       RGroupLogicEntry { label: 1, dependent_label: Some(2), rgroup_or_h: false, occurrence: vec![RGroupOccurrence::GreaterThan(0)] })]
fn test_rgroup_logic_entry(#[case] input: &[u8], #[case] expected: RGroupLogicEntry) {
    let result = all_consuming(rgroup_logic_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::count_is_zero(b"  0   1   0", NomErrorKind::Verify)]
#[case::count_exceeds_1(b"  2   1   0", NomErrorKind::Verify)]
#[case::label_is_zero(b"  1   0   0", NomErrorKind::Verify)]
#[case::rgroup_or_h_out_of_range(b"  1   1   0   2", NomErrorKind::Verify)]
fn test_rgroup_logic_entry_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = all_consuming(rgroup_logic_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case::single_entry(b"  1   1 SUP", vec![SGroupTypeEntry { sgroup_index: 0, sgroup_type: SGroupType::Superatom }])]
#[case::two_entries(b"  2   1 SUP   2 DAT", vec![
    SGroupTypeEntry { sgroup_index: 0, sgroup_type: SGroupType::Superatom },
    SGroupTypeEntry { sgroup_index: 1, sgroup_type: SGroupType::Data }
])]
fn test_sgroup_type_entries(#[case] input: &[u8], #[case] expected: Vec<SGroupTypeEntry>) {
    let result = all_consuming(sgroup_type_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case::invalid_sgroup_type(b"  1   1 FOO", NomErrorKind::MapRes)]
#[case::trailing_characters(b"  1   1 SUP a", NomErrorKind::Eof)]
fn test_sgroup_type_entries_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let result = all_consuming(sgroup_type_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case::single_entry(b"  1   1 ALT", vec![SGroupSubtypeEntry { sgroup_index: 0, sgroup_subtype: SGroupSubtype::Alternating }])]
#[case::two_entries(b"  2   1 RAN   2 BLO", vec![
    SGroupSubtypeEntry { sgroup_index: 0, sgroup_subtype: SGroupSubtype::Random },
    SGroupSubtypeEntry { sgroup_index: 1, sgroup_subtype: SGroupSubtype::Block }
])]
fn test_sgroup_subtype_entries(#[case] input: &[u8], #[case] expected: Vec<SGroupSubtypeEntry>) {
    let result = all_consuming(sgroup_subtype_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case::invalid_sgroup_subtype(b"  1   1 FOO", NomErrorKind::MapRes)]
#[case::trailing_characters(b"  1   1 ALT a", NomErrorKind::Eof)]
fn test_sgroup_subtype_entries_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let result = all_consuming(sgroup_subtype_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case::single_entry(b"  1   1   1", vec![SGroupLabelEntry { sgroup_index: 0, label: 1 }])]
#[case::two_entries(b"  2   1  14   2  15",
       vec![SGroupLabelEntry { sgroup_index: 0, label: 14 }, SGroupLabelEntry { sgroup_index: 1, label: 15 }])]
fn test_sgroup_label_entries(#[case] input: &[u8], #[case] expected: Vec<SGroupLabelEntry>) {
    let result = all_consuming(sgroup_label_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case::label_out_of_range(b"  1   1   0", NomErrorKind::Verify)]
#[case::label_out_of_range_high(b"  1   1 513", NomErrorKind::Verify)]
#[case::trailing_characters(b"  1   1   1 a", NomErrorKind::Eof)]
fn test_sgroup_label_entries_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let result = all_consuming(sgroup_label_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
//       nn8 sss ttt sss ttt sss ttt
#[case::two_entries(b"  2   1 HT    2 HT ", vec![
    SGroupConnectivityEntry { sgroup_index: 0, connectivity: SGroupConnectivity::HeadToTail },
    SGroupConnectivityEntry { sgroup_index: 1, connectivity: SGroupConnectivity::HeadToTail },
])]
#[case::partial_last_entry(b"  3   1 HT    2 HT    3 HT", vec![
    SGroupConnectivityEntry { sgroup_index: 0, connectivity: SGroupConnectivity::HeadToTail },
    SGroupConnectivityEntry { sgroup_index: 1, connectivity: SGroupConnectivity::HeadToTail },
    SGroupConnectivityEntry { sgroup_index: 2, connectivity: SGroupConnectivity::HeadToTail }])]
fn test_sgroup_connectivity_entries(
    #[case] input: &[u8],
    #[case] expected: Vec<SGroupConnectivityEntry>,
) {
    let result = all_consuming(sgroup_connectivity_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case::count_is_zero(b"  0", NomErrorKind::Verify)]
#[case::invalid_connectivity(b"  1   1 FOO", NomErrorKind::MapRes)]
#[case::trailing_characters(b"  1   1 HT a", NomErrorKind::Eof)]
fn test_sgroup_connectivity_entries_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = all_consuming(sgroup_connectivity_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case::single_entry(b" EXP  1   1", vec![SGroupExpansionEntry { sgroup_index: 0 }])]
#[case::two_entries(b" EXP  2   1   2", vec![
    SGroupExpansionEntry { sgroup_index: 0 },
    SGroupExpansionEntry { sgroup_index: 1 }
])]
fn test_sgroup_expansion_entries(
    #[case] input: &[u8],
    #[case] expected: Vec<SGroupExpansionEntry>,
) {
    let result = all_consuming(sgroup_expansion_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case::count_is_zero(b" EXP  0   1", NomErrorKind::Verify)]
#[case::trailing_characters(b" EXP  1   1 a", NomErrorKind::Eof)]
fn test_sgroup_expansion_entries_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = all_consuming(sgroup_expansion_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case::two_entries(b"   1  2   1   2", SGroupAtomListEntry { sgroup_index: 0, atom_indices: vec![0, 1] })]
#[case::single_entry(b"   3  1  15", SGroupAtomListEntry { sgroup_index: 2, atom_indices: vec![14] })]
fn test_sgroup_atom_list_entry(#[case] input: &[u8], #[case] expected: SGroupAtomListEntry) {
    let result = all_consuming(sgroup_atom_list_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case::count_exceeds_15(b"   1 16   1", NomErrorKind::Verify)]
#[case::trailing_characters(b"   1  1   1 a", NomErrorKind::Eof)]
fn test_sgroup_atom_list_entry_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let result = all_consuming(sgroup_atom_list_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case::two_entries(b"   1  2   1   2", SGroupBondListEntry { sgroup_index: 0, bond_indices: vec![0, 1] })]
#[case::single_entry(b"   3  1  15", SGroupBondListEntry { sgroup_index: 2, bond_indices: vec![14] })]
fn test_sgroup_bond_list_entry(#[case] input: &[u8], #[case] expected: SGroupBondListEntry) {
    let result = all_consuming(sgroup_bond_list_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case::count_exceeds_15(b"   1 16   1", NomErrorKind::Verify)]
#[case::trailing_characters(b"   1  1   1 a", NomErrorKind::Eof)]
fn test_sgroup_bond_list_entry_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let result = all_consuming(sgroup_bond_list_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case::four_entries(b"   1  4   3   4   5   6", SGroupParentAtomEntry { sgroup_index: 0, atom_indices: vec![2, 3, 4, 5] })]
fn test_sgroup_parent_atom_entries(#[case] input: &[u8], #[case] expected: SGroupParentAtomEntry) {
    let result = all_consuming(sgroup_parent_atom_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case::count_exceeds_15(b"   1 16   3", NomErrorKind::Verify)]
#[case::trailing_characters(b"   1  4   3   4   5   6 a", NomErrorKind::Eof)]
fn test_sgroup_parent_atom_entries_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = all_consuming(sgroup_parent_atom_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case::monomer(b"   1 1", SGroupSubscriptEntry { sgroup_index: 0, multiplier: Some(SGroupMultiplier::Single(SGroupMultiplierTerm::Integer(1))), subscript: Some("1".to_string()) })]
#[case::n_mer(b"   1 n", SGroupSubscriptEntry { sgroup_index: 0, multiplier: Some(SGroupMultiplier::Single(SGroupMultiplierTerm::Variable('n'))), subscript: Some("n".to_string()) })]
#[case::ph_subscript(b"   1 Ph", SGroupSubscriptEntry { sgroup_index: 0, multiplier: Some(SGroupMultiplier::Expression
     { left: SGroupMultiplierTerm::Variable('P'), op: SGroupMultiplierOp::Mul, right: SGroupMultiplierTerm::Variable('h') }), subscript: Some("Ph".to_string()) })]
fn test_sgroup_subscript_entry(#[case] input: &[u8], #[case] expected: SGroupSubscriptEntry) {
    let result = all_consuming(sgroup_subscript_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case::sgroup_index_is_zero(b"   0 1", NomErrorKind::Verify)]
fn test_sgroup_subscript_entry_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let result = all_consuming(sgroup_subscript_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case::three_entries(b"   3  3  10   9   4", SGroupCorrespondenceEntry { sgroup_index: 2, bond_indices: vec![9, 8, 3] })]
fn test_sgroup_correspondence_entry(
    #[case] input: &[u8],
    #[case] expected: SGroupCorrespondenceEntry,
) {
    let result = all_consuming(sgroup_correspondence_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case::count_is_zero(b"   3  0", NomErrorKind::Verify)]
#[case::trailing_characters(b"   3  3  10   9   4 a", NomErrorKind::Eof)]
fn test_sgroup_correspondence_entry_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = all_consuming(sgroup_correspondence_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case::all_coordinates(b"   1  4  -13.0153    4.4289  -13.0153    8.2211",
       SGroupDisplayInfoEntry { sgroup_index: 0,  bracket_coords: vec![-13.0153, 4.4289, -13.0153, 8.2211]})]
fn test_sgroup_display_info_entry(#[case] input: &[u8], #[case] expected: SGroupDisplayInfoEntry) {
    let result = all_consuming(sgroup_display_info_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result.sgroup_index, expected.sgroup_index);
    for i in 0..result.bracket_coords.len() {
        assert!(approx_eq!(
            f64,
            result.bracket_coords[i],
            expected.bracket_coords[i]
        ));
    }
}

#[rustfmt::skip]
#[rstest]
#[case::sgroup_index_is_zero(b"   0  4    4.4700   -3.1700    4.4700   -5.7500", NomErrorKind::Verify)]
#[case::count_is_zero(b"   1  0", NomErrorKind::Verify)]
fn test_sgroup_display_info_entry_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = all_consuming(sgroup_display_info_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"   1   6   -0.7200   -0.4200", SGroupConnectingBondEntry { sgroup_index: 0, bond_index: 5, bond_vector: (-0.7200, -0.4200) })]
fn test_sgroup_connecting_bond_entry(
    #[case] input: &[u8],
    #[case] expected: SGroupConnectingBondEntry,
) {
    let (remaining, result) = all_consuming(sgroup_connecting_bond_entry())
        .parse(input)
        .unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result.sgroup_index, expected.sgroup_index);
    assert_eq!(result.bond_index, expected.bond_index);
    assert!(approx_eq!(
        f64,
        result.bond_vector.0,
        expected.bond_vector.0
    ));
    assert!(approx_eq!(
        f64,
        result.bond_vector.1,
        expected.bond_vector.1
    ));
}

#[rustfmt::skip]
#[rstest]
#[case::sgroup_index_is_zero(b"   0   1   -0.7200   -0.4200", NomErrorKind::Verify)]
#[case::bond_index_is_zero(b"   1   0   -0.7200   -0.4200", NomErrorKind::Verify)]
#[case::trailing_characters(b"   1   1   -0.7200   -0.4200 a", NomErrorKind::Eof)]
fn test_sgroup_connecting_bond_entry_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = all_consuming(sgroup_connecting_bond_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case::text_field(b"   1 pH   ",
       SGroupDataDescriptionEntry { sgroup_index: 0, field_name: "pH".to_string(), field_type: SGroupDataType::Text, field_units: None, query_identifier: None, data_query_operator: None })]
#[case::marvin_extensions(b"   3 MRV_COORDINATE_BOND_TYPE                              ",
       SGroupDataDescriptionEntry { sgroup_index: 2, field_name: "MRV_COORDINATE_BOND_TYPE".to_string(), field_type: SGroupDataType::Text, field_units: None, query_identifier: None, data_query_operator: None })]
#[case::numerical_field(b"   3 WEIGHT_PERCENT                N %",
       SGroupDataDescriptionEntry { sgroup_index: 2, field_name: "WEIGHT_PERCENT".to_string(), field_type: SGroupDataType::Numeric, field_units: Some("%".to_string()), query_identifier: None, data_query_operator: None })]
fn test_sgroup_data_description_entry(
    #[case] input: &[u8],
    #[case] expected: SGroupDataDescriptionEntry,
) {
    let result = all_consuming(sgroup_data_description_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::sgroup_index_is_zero(b"   0 pH   ", NomErrorKind::Verify)]
#[case::invalid_field_type(b"   1 pH                            X", NomErrorKind::MapRes)]
fn test_sgroup_data_description_entry_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = all_consuming(sgroup_data_description_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rstest]
#[case::detached(b"   1     0.0000    0.0000    DR    ALL  0       0",
        SGroupDataDisplayEntry { sgroup_index: 0, coords: (0.0000, 0.0000), display_type: SGroupDataDisplayType::Detached,
        display_placement: SGroupDataDisplayPlacement::Relative, display_units: SGroupDataDisplayUnits::None,
        display_chars: SGroupDataDisplayChars::All, display_tag: None, display_position: 0 })]
#[case::relative_position(b"   2     0.0000    0.0000    DR    ALL  1       6",
        SGroupDataDisplayEntry { sgroup_index: 1, coords: (0.0000, 0.0000), display_type: SGroupDataDisplayType::Detached,
        display_placement: SGroupDataDisplayPlacement::Relative, display_units: SGroupDataDisplayUnits::None,
        display_chars: SGroupDataDisplayChars::All, display_tag: None, display_position: 6 })]
#[case::number(b"   2     0.0000    0.0000    DR    1    1       6",
        SGroupDataDisplayEntry { sgroup_index: 1, coords: (0.0000, 0.0000), display_type: SGroupDataDisplayType::Detached,
        display_placement: SGroupDataDisplayPlacement::Relative, display_units: SGroupDataDisplayUnits::None,
        display_chars: SGroupDataDisplayChars::Number(1), display_tag: None, display_position: 6 })]
#[case::absolute_position(b"   3     0.0000    0.0000    DR    ALL  0       0 ",
        SGroupDataDisplayEntry { sgroup_index: 2, coords: (0.0000, 0.0000), display_type: SGroupDataDisplayType::Detached,
        display_placement: SGroupDataDisplayPlacement::Relative, display_units: SGroupDataDisplayUnits::None,
        display_chars: SGroupDataDisplayChars::All, display_tag: None, display_position: 0 })]
fn test_sgroup_data_display_entry(#[case] input: &[u8], #[case] expected: SGroupDataDisplayEntry) {
    let result = all_consuming(terminated(sgroup_data_display_entry(true), space0)).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::sgroup_index_is_zero(b"   0     0.0000    0.0000    DR    ALL  0       0", NomErrorKind::Verify)]
#[case::invalid_display_type(b"   1     0.0000    0.0000    XR    ALL  0       0", NomErrorKind::MapRes)]
#[case::invalid_placement_type(b"   1     0.0000    0.0000    DX    ALL  0       0", NomErrorKind::MapRes)]
#[case::invalid_chars_type(b"   1     0.0000    0.0000    DR    NON  0       0", NomErrorKind::Digit)]
fn test_sgroup_data_display_entry_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = all_consuming(sgroup_data_display_entry(true)).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected: {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rstest]
#[case::numerical(b"   1 4.6", SGroupDataEntry::Continuation { sgroup_index: 0, data_content: "4.6".to_string() })]
#[case::text(b"   2 E/Z unknown", SGroupDataEntry::Continuation { sgroup_index: 1, data_content: "E/Z unknown".to_string() })]
#[case::empty(b"   1", SGroupDataEntry::Continuation { sgroup_index: 0, data_content: "".to_string() })]
fn tests_sgroup_data_continuation_entry(#[case] input: &[u8], #[case] expected: SGroupDataEntry) {
    let result = all_consuming(sgroup_data_continuation_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case::sgroup_index_is_zero(b"   0 4.6", NomErrorKind::Verify)]
fn test_sgroup_data_continuation_entry_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = all_consuming(sgroup_data_continuation_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case::numerical(b"   1 4.6", SGroupDataEntry::EndWithData { sgroup_index: 0, data_content: "4.6".to_string() })]
#[case::text(b"   2 E/Z unknown", SGroupDataEntry::EndWithData { sgroup_index: 1, data_content: "E/Z unknown".to_string() })]
#[case::empty(b"   1", SGroupDataEntry::EndBlank { sgroup_index: 0 })]
fn tests_sgroup_data_end_entry(#[case] input: &[u8], #[case] expected: SGroupDataEntry) {
    let result = all_consuming(sgroup_data_end_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case::sgroup_index_is_zero(b"   0  1", NomErrorKind::Verify)]
fn test_sgroup_data_end_entry_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let result = all_consuming(sgroup_data_end_entry()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rstest]
#[case::one_entry(b"  1   1   2", vec![SGroupHierarchyEntry { sgroup_index: 0, parent_sgroup_index: 1 }])]
#[case::multiple_entries(b"  3   1   4   2   4   3   2", vec![
    SGroupHierarchyEntry { sgroup_index: 0, parent_sgroup_index: 3 },
    SGroupHierarchyEntry { sgroup_index: 1, parent_sgroup_index: 3 },
    SGroupHierarchyEntry { sgroup_index: 2, parent_sgroup_index: 1 }
])]
fn test_sgroup_hierarchy_entries(
    #[case] input: &[u8],
    #[case] expected: Vec<SGroupHierarchyEntry>,
) {
    let result = all_consuming(sgroup_hierarchy_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::count_is_zero(b"  0   1   2", NomErrorKind::Verify)]
#[case::sgroup_index_is_zero(b"  1   0   2", NomErrorKind::Verify)]
#[case::parent_sgroup_index_is_zero(b"  1   1   0", NomErrorKind::Verify)]
#[case::trailing_characters(b"  1   1   2 a", NomErrorKind::Eof)]
fn test_sgroup_hierarchy_entries_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = all_consuming(sgroup_hierarchy_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case::multiple_entries(b"  2   1   1   2   2", vec![
    SGroupComponentEntry { sgroup_index: 0, component_number: 1 },
    SGroupComponentEntry { sgroup_index: 1, component_number: 2 }
])]
fn test_sgroup_component_entries(
    #[case] input: &[u8],
    #[case] expected: Vec<SGroupComponentEntry>,
) {
    let result = all_consuming(sgroup_component_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case::count_is_zero(b"  0   1   2", NomErrorKind::Verify)]
#[case::component_number_is_zero(b"  1   0   2", NomErrorKind::Verify)]
#[case::trailing_characters(b"  1   1   0 a", NomErrorKind::Eof)]
fn test_sgroup_component_entries_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = all_consuming(sgroup_component_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rstest]
#[case::one_entry(b"  1   1   0", vec![BondOrderOverrideEntry { bond_index: 0, bond_order: BondOrder::Zero }])]
#[case::multiple_entries(b"  2   1   2   3   4", vec![
    BondOrderOverrideEntry { bond_index: 0, bond_order: BondOrder::Double },
    BondOrderOverrideEntry { bond_index: 2, bond_order: BondOrder::Quadruple },
])]
fn test_bond_order_override_entries(
    #[case] input: &[u8],
    #[case] expected: Vec<BondOrderOverrideEntry>,
) {
    let result = all_consuming(bond_order_override_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case::count_is_zero(b"  0   1   0", NomErrorKind::Verify)]
#[case::bond_index_is_zero(b"  1   0   0", NomErrorKind::Verify)]
#[case::bond_order_out_of_range(b"  1   1   7", NomErrorKind::Verify)]
#[case::trailing_characters(b"  1   1   0 a", NomErrorKind::Eof)]
fn test_bond_order_override_entries_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let result = all_consuming(bond_order_override_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rstest]
#[case::one_entry(b"  1   1   0", vec![AtomChargeOverrideEntry { atom_index: 0, charge: 0 }])]
#[case::multiple_entries(b"  2   1   2   3  -1", vec![
    AtomChargeOverrideEntry { atom_index: 0, charge: 2 },
    AtomChargeOverrideEntry { atom_index: 2, charge: -1 },
])]
fn test_atom_charge_override_entries(
    #[case] input: &[u8],
    #[case] expected: Vec<AtomChargeOverrideEntry>,
) {
    let result = all_consuming(atom_charge_overrides_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case::count_is_zero(b"  0   1   0", NomErrorKind::Verify)]
#[case::atom_index_is_zero(b"  1   0   0", NomErrorKind::Verify)]
#[case::hydrogen_count_out_of_range(b"  1   1   9", NomErrorKind::Verify)]
#[case::trailing_characters(b"  1   1   0 a", NomErrorKind::Eof)]
fn test_atom_charge_override_entries_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = all_consuming(atom_charge_overrides_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rstest]
#[case::one_entry(b"  1   1   0", vec![AtomHydrogenCountEntry { atom_index: 0, hydrogen_count: Some(0) }])]
#[case::no_override(b"  1   1  -1", vec![AtomHydrogenCountEntry { atom_index: 0, hydrogen_count: None }])]
#[case::multiple_entries(b"  2   1   2   3   4", vec![
    AtomHydrogenCountEntry { atom_index: 0, hydrogen_count: Some(2) },
    AtomHydrogenCountEntry { atom_index: 2, hydrogen_count: Some(4) },
])]
fn test_atom_hydrogen_count_entries(
    #[case] input: &[u8],
    #[case] expected: Vec<AtomHydrogenCountEntry>,
) {
    let result = all_consuming(atom_hydrogen_count_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case::count_is_zero(b"  0   1   0", NomErrorKind::Eof)]
#[case::atom_index_is_zero(b"  1   0   0", NomErrorKind::Verify)]
#[case::hydrogen_count_out_of_range(b"  1   1   9", NomErrorKind::Verify)]
#[case::trailing_characters(b"  1   1   0 a", NomErrorKind::Eof)]
fn test_atom_hydrogen_count_entries_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = all_consuming(atom_hydrogen_count_entries()).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {:?}, expected {:?}, got {}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}
