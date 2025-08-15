use super::*;
use pretty_assertions::assert_eq;
use nom::{error, Err};
use rstest::*;
use umol_data::Element;
use float_cmp::approx_eq;

#[rstest]
#[case(b"  1 F    3   9   7   8  ", PropertyEntries::AtomListEntry(AtomListEntry { atom_index: 0, exclusion: false, elements: vec![Element::F, Element::N, Element::O] }))]   
#[case(b"  1 T    3   9   7   8  ", PropertyEntries::AtomListEntry(AtomListEntry { atom_index: 0, exclusion: true, elements: vec![Element::F, Element::N, Element::O] }))]   
#[case(b"  1      3   9   7   8  ", PropertyEntries::AtomListEntry(AtomListEntry { atom_index: 0, exclusion: false, elements: vec![Element::F, Element::N, Element::O] }))]   
#[case(b"  1 F    3   9   7   8", PropertyEntries::AtomListEntry(AtomListEntry { atom_index: 0, exclusion: false, elements: vec![Element::F, Element::N, Element::O] }))]   
#[case(b"  1 T    3   9   7   8", PropertyEntries::AtomListEntry(AtomListEntry { atom_index: 0, exclusion: true, elements: vec![Element::F, Element::N, Element::O] }))]   
fn test_legacy_atom_list_input(#[case] input: &[u8], #[case] expected: PropertyEntries) {
    let (remaining, result) = legacy_atom_list_input().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  1 F    0   9  ", "count is zero", error::ErrorKind::Verify)]
#[case(b"  1 F    6    9  ", "count exceeds 5", error::ErrorKind::Verify)]
#[case(b"  1 X    1   9  ", "invalid exclusion flag", error::ErrorKind::MapRes)]
#[case(b"  1 F    1   0  ", "invalid element atomic number", error::ErrorKind::MapOpt)]
fn test_legacy_atom_list_input_invalid(#[case] input: &[u8], #[case] desc: &str, #[case] expected_kind: error::ErrorKind) {
    let result = legacy_atom_list_input().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.as_ref(), Err(Err::Error(e)) if e.code == expected_kind),
        "Expected {:?} error for {}, got {:?}",
        expected_kind,
        desc,
        result
    );
}

#[rstest]
#[case(b"A    1\nCF3", "A atom alias property", PropertyEntries::AtomAliasEntry(AtomAliasEntry { atom_index: 0, alias: "CF3".to_string() }))]
#[case(b"V    1 *", "V atom value property", PropertyEntries::AtomValueEntry(AtomValueEntry { atom_index: 0, value: "*".to_string() }))]
#[case(b"M  CHG  1   1  -1", "CHG standard property", PropertyEntries::ChargeEntries(vec![ChargeEntry { atom_index: 0, charge: -1 }]))]
#[case(b"M  RAD  1   1   2", "RAD standard property", PropertyEntries::RadicalEntries(vec![RadicalEntry { atom_index: 0, radical_type: 2 }]))]
#[case(b"M  ISO  1   1  13", "ISO standard property", PropertyEntries::IsotopeEntries(vec![IsotopeEntry { atom_index: 0, mass: 13 }]))]
#[case(b"M  RBC  1   1   2", "RBC query property", PropertyEntries::RingBondCountEntries(vec![RingBondCountEntry { atom_index: 0, ring_bond_count: 2 }]))]
#[case(b"M  SUB  1   1   3", "SUB query property", PropertyEntries::SubstitutionCountEntries(vec![SubstitutionCountEntry { atom_index: 0, substitution_count: 3 }]))]
#[case(b"M  UNS  1   1   1", "UNS query property", PropertyEntries::UnsaturatedAtomEntries(vec![UnsaturatedAtomEntry { atom_index: 0, unsaturated: 1 }]))]
#[case(b"M  LIN  1   1   2   5   7", "LIN query property", PropertyEntries::LinkAtomEntries(vec![LinkAtomEntry { atom_index: 0, repeat_count: 2, subs_index1: 4, subs_index2: Some(6) }]))]
#[case(b"M  ALS   1  3   F   N   O   ", "ALS normal (F) query property", PropertyEntries::AtomListEntry(AtomListEntry { atom_index: 0, exclusion: false, elements: vec![Element::F, Element::N, Element::O] }))]
#[case(b"M  ALS   1  3 T F   N   O   ", "ALS exclusion (T) query property", PropertyEntries::AtomListEntry(AtomListEntry { atom_index: 0, exclusion: true, elements: vec![Element::F, Element::N, Element::O] }))]
#[case(b"M  APO  1   1   1", "APO query property", PropertyEntries::AttachmentPointEntries(vec![AttachmentPointEntry { atom_index: 0, attachment_type: 1 }]))]
#[case(b"M  AAL   4  2  14   1   9   2", "AAL query property", PropertyEntries::AtomAttachmentOrderEntry(AtomAttachmentOrderEntry { atom_index: 3, attachments: vec![(13, 1), (8, 2)] }))]
#[case(b"M  RGP  1   1   2", "RGP RGroup property", PropertyEntries::RGroupLabelEntries(vec![RGroupLabelEntry { atom_index: 0, label: 2 }]))]
#[case(b"M  LOG  1   1   0   0  >2", "LOG RGroup property", PropertyEntries::RGroupLogicEntry(RGroupLogicEntry { label: 1, dependent_label: None, rgroup_or_h: false, occurrence: vec![RGroupOccurrence::GreaterThan(2)] }))]
#[case(b"M  STY  1   1 SUP", "STY SGroup property", PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry { sgroup_index: 0, sgroup_type: SGroupType::Superatom }]))]
#[case(b"M  SST  1   1 ALT", "SST SGroup property", PropertyEntries::SGroupSubtypeEntries(vec![SGroupSubtypeEntry { sgroup_index: 0, sgroup_subtype: SGroupSubtype::Alternating }]))]
#[case(b"M  SLB  1   1  19", "SLB SGroup property", PropertyEntries::SGroupLabelEntries(vec![SGroupLabelEntry { sgroup_index: 0, label: 19 }]))]
#[case(b"M  SCN  1   1 HH ", "SCN SGroup property", PropertyEntries::SGroupConnectivityEntries(vec![SGroupConnectivityEntry { sgroup_index: 0, connectivity: SGroupConnectivity::HeadToHead }]))]
#[case(b"M  SDS EXP  1   1", "SDS SGroup property", PropertyEntries::SGroupExpansionEntries(vec![SGroupExpansionEntry { sgroup_index: 0 }]))]
#[case(b"M  SAL   1  1   5", "SAL SGroup property", PropertyEntries::SGroupAtomListEntry(SGroupAtomListEntry { sgroup_index: 0, atom_indices: vec![4] }))]
#[case(b"M  SBL   1  1   3", "SBL SGroup property", PropertyEntries::SGroupBondListEntry(SGroupBondListEntry { sgroup_index: 0, bond_indices: vec![2] }))]
#[case(b"M  SPA   1 12   3   4   5   6   9  10  11  12  13  14  15  16", "SPA SGroup property", PropertyEntries::SGroupParentAtomEntry(SGroupParentAtomEntry { sgroup_index: 0,
    atom_indices: vec![2, 3, 4, 5, 8, 9, 10, 11, 12, 13, 14, 15] }))]
#[case(b"M  SMT   1 n", "SMT SGroup property", PropertyEntries::SGroupSubscriptEntry(SGroupSubscriptEntry { sgroup_index: 0, data: SGroupSubscriptData::Multiplier(SGroupMultiplier::N) }))]
#[case(b"M  CRS   1  3  10   9   4", "CRS SGroup property", PropertyEntries::SGroupCorrespondenceEntry(SGroupCorrespondenceEntry { sgroup_index: 0, bond_indices: vec![9, 8, 3] }))]
#[case(b"M  SDI   3  4    4.4700   -3.1700    4.4700   -5.7500", "SDI SGroup property", PropertyEntries::SGroupDisplayInfoEntry(SGroupDisplayInfoEntry { sgroup_index: 2, bracket_coords: vec![4.4700, -3.1700, 4.4700, -5.7500] }))]
#[case(b"M  SBV   1  11    0.6400    0.9700", "SBV SGroup property", PropertyEntries::SGroupConnectingBondEntry(SGroupConnectingBondEntry { sgroup_index: 0, bond_index: 10, bond_vector: (0.6400, 0.9700) }))]
#[case(b"M  SDT   1 pH   ", "SDT SGroup property", PropertyEntries::SGroupDataDescriptionEntry(SGroupDataDescriptionEntry { sgroup_index: 0, field_name: "pH".to_string(), 
    field_type: SGroupDataType::Text, field_units: None, query_identifier: None, data_query_operator: None }))]
#[case(b"M  SDD   1     0.0000    0.0000    DR    ALL  1       6", "SDD SGroup property", PropertyEntries::SGroupDataDisplayEntry(SGroupDataDisplayEntry { sgroup_index: 0,
    coords: (0.0000, 0.0000), display_type: SGroupDataDisplayType::Detached, display_placement: SGroupDataDisplayPlacement::Relative, display_units: SGroupDataDisplayUnits::None,
    display_chars: SGroupDataDisplayChars::All, display_tag: None, display_position: 6 }))]
#[case(b"M  SCD   1 4.6", "SCD SGroup property", PropertyEntries::SGroupDataEntry(SGroupDataEntry::Continuation { sgroup_index: 0, data_content: "4.6".to_string() }))]
#[case(b"M  SED   2 E/Z unknown", "SED SGroup property", PropertyEntries::SGroupDataEntry(SGroupDataEntry::EndWithData { sgroup_index: 1, data_content: "E/Z unknown".to_string() }))]
#[case(b"M  SED   1", "SED SGroup property", PropertyEntries::SGroupDataEntry(SGroupDataEntry::EndBlank { sgroup_index: 0 }))]
#[case(b"M  SPL  1   1   2", "SPL SGroup property", PropertyEntries::SGroupHierarchyEntries(vec![SGroupHierarchyEntry { sgroup_index: 0, parent_sgroup_index: 1 }]))]
#[case(b"M  SNC  2   1   1   2   2", "SNC SGroup property", PropertyEntries::SGroupComponentEntries(vec![SGroupComponentEntry { sgroup_index: 0, component_number: 1 }, SGroupComponentEntry { sgroup_index: 1, component_number: 2 }]))]
#[case(b"M  ZBO  1   1   0", "ZBO bond property", PropertyEntries::ZeroBondOrderEntries(vec![ZeroBondOrderEntry { bond_index: 0, bond_order: 0 }]))]
#[case(b"M  ZCH  1   1  -1", "ZCH atom property", PropertyEntries::ZeroAtomChargeEntries(vec![ZeroAtomChargeEntry { atom_index: 0, charge: -1 }]))]
#[case(b"M  HYD  1   1   1", "HYD atom property", PropertyEntries::AtomHydrogenCountEntries(vec![AtomHydrogenCountEntry { atom_index: 0, hydrogen_count: 1 }]))]
fn test_property_input(#[case] input: &[u8], #[case] desc: &str, #[case] expected: PropertyEntries) {
    let (remaining, result) = property_input().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty for {}", desc);
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"M  CHG  1   1  -1  a", "trailing chars", error::ErrorKind::Eof)]
#[case(b"M  CHG  2   1  -1", "count does not match item list", error::ErrorKind::Eof)]
#[case(b"M  CHG  0", "count is zero", error::ErrorKind::Verify)]
#[case(b"M  CHG  1   0 -10", "atom index is zero", error::ErrorKind::Verify)]
#[case(b"M  XXX  1   1  -1", "invalid property tag", error::ErrorKind::Tag)]
fn test_property_input_invalid(#[case] input: &[u8], #[case] desc: &str, #[case] expected_kind: error::ErrorKind) {
    let result = property_input().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.as_ref(), Err(Err::Error(e)) if e.code == expected_kind),
        "Expected {:?} error for {}, got {:?}",
        expected_kind,
        desc,
        result
    );
}

#[rstest]
#[case(b"A    1\nCF3", "A atom alias property", PropertyEntries::AtomAliasEntry(AtomAliasEntry { atom_index: 0, alias: "CF3".to_string() }))]
#[case(b"V    1 *", "V atom value property", PropertyEntries::AtomValueEntry(AtomValueEntry { atom_index: 0, value: "*".to_string() }))]
#[case(b"M  CHG  1   1  -1", "CHG atom property", PropertyEntries::ChargeEntries(vec![ChargeEntry { atom_index: 0, charge: -1 }]))]
#[case(b"M  RAD  1   1   2", "RAD atom property", PropertyEntries::RadicalEntries(vec![RadicalEntry { atom_index: 0, radical_type: 2 }]))]
#[case(b"M  ISO  1   1  13", "ISO atom property", PropertyEntries::IsotopeEntries(vec![IsotopeEntry { atom_index: 0, mass: 13 }]))]
#[case(b"M  STY  1   1 SUP", "STY SGroup property", PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry { sgroup_index: 0, sgroup_type: SGroupType::Superatom }]))]
#[case(b"M  SST  1   1 ALT", "SST SGroup property", PropertyEntries::SGroupSubtypeEntries(vec![SGroupSubtypeEntry { sgroup_index: 0, sgroup_subtype: SGroupSubtype::Alternating }]))]
#[case(b"M  SLB  1   1  19", "SLB SGroup property", PropertyEntries::SGroupLabelEntries(vec![SGroupLabelEntry { sgroup_index: 0, label: 19 }]))]
#[case(b"M  SAL   1  1   5", "SAL SGroup property", PropertyEntries::SGroupAtomListEntry(SGroupAtomListEntry { sgroup_index: 0, atom_indices: vec![4] }))]
#[case(b"M  SBL   1  1   3", "SBL SGroup property", PropertyEntries::SGroupBondListEntry(SGroupBondListEntry { sgroup_index: 0, bond_indices: vec![2] }))]
#[case(b"M  SMT   1 n", "SMT SGroup property", PropertyEntries::SGroupSubscriptEntry(SGroupSubscriptEntry { sgroup_index: 0, data: SGroupSubscriptData::Multiplier(SGroupMultiplier::N) }))]
#[case(b"M  ZBO  1   1   0", "ZBO bond property", PropertyEntries::ZeroBondOrderEntries(vec![ZeroBondOrderEntry { bond_index: 0, bond_order: 0 }]))]
#[case(b"M  ZCH  1   1  -1", "ZCH atom property", PropertyEntries::ZeroAtomChargeEntries(vec![ZeroAtomChargeEntry { atom_index: 0, charge: -1 }]))]
#[case(b"M  HYD  1   1   1", "HYD atom property", PropertyEntries::AtomHydrogenCountEntries(vec![AtomHydrogenCountEntry { atom_index: 0, hydrogen_count: 1 }]))]
fn test_property_input_standard(#[case] input: &[u8], #[case] desc: &str, #[case] expected: PropertyEntries) {
    let (remaining, result) = property_input_standard().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty for {}", desc);
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"M  RBC  1   1   2", "RBC query property not supported in standard parser")]
#[case(b"M  SUB  1   1   3", "SUB query property not supported in standard parser")]
#[case(b"M  UNS  1   1   1", "UNS query property not supported in standard parser")]
#[case(b"M  LIN  1   1   2   5   7", "LIN query property not supported in standard parser")]
#[case(b"M  ALS  1  3FC   N   O   ", "ALS query property not supported in standard parser")]
#[case(b"M  APO  1   1   1", "APO query property not supported in standard parser")]
#[case(b"M  AAL  1 1   2   1", "AAL query property not supported in standard parser")]
#[case(b"M  RGP   1   1   2", "RGP RGroup property not supported in standard parser")]
#[case(b"M  LOG   1   1   0   0  >2", "LOG RGroup property not supported in standard parser")]
#[case(b"M  SCN  1   1 HH ", "SCN SGroup property not supported in standard parser")]
#[case(b"M  SDS EXP  1   1", "SDS SGroup property not supported in standard parser")]
#[case(b"M  CRS   1  3  10   9   4", "CRS SGroup property not supported in standard parser")]
#[case(b"M  SDI   3  4    4.4700   -3.1700    4.4700   -5.7500", "SDI SGroup property not supported in standard parser")]
#[case(b"M  SBV   1  11    0.6400    0.9700", "SBV SGroup property not supported in standard parser")]
#[case(b"M  SDT   1 pH   ", "SDT SGroup property not supported in standard parser")]
#[case(b"M  SDD   1     0.0000    0.0000    DR    ALL  1       6", "SDD SGroup property not supported in standard parser")]
#[case(b"M  SCD   1   1   0", "SCD SGroup property not supported in standard parser")]
#[case(b"M  SED   1   1   0", "SED SGroup property not supported in standard parser")]
#[case(b"M  SPL   1   1   0", "SPL SGroup property not supported in standard parser")]
#[case(b"M  SNC   1   1   0", "SNC SGroup property not supported in standard parser")]
fn test_property_input_standard_invalid(#[case] input: &[u8], #[case] desc: &str) {
    let result = property_input_standard().parse(input);
    assert!(result.is_err(), "{}", desc);
    assert!(
        matches!(result.as_ref(), Err(Err::Error(e)) if e.code == error::ErrorKind::Tag),
        "Expected Tag error for {}, got {:?}",
        desc,
        result
    );
}

#[rstest]
#[case(b"  1\nCF3", AtomAliasEntry { atom_index: 0, alias: "CF3".to_string() })]
#[case(b" 15\n  Et", AtomAliasEntry { atom_index: 14, alias: "Et".to_string() })]
fn test_atom_alias_entry(#[case] input: &[u8], #[case] expected: AtomAliasEntry) {
    let (remaining, result) = atom_alias_entry().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  0\n  Et", "atom index is zero", error::ErrorKind::Verify)]
fn test_atom_alias_entry_invalid(#[case] input: &[u8], #[case] desc: &str, #[case] expected_kind: error::ErrorKind) {
    let result = atom_alias_entry().parse(input);
    assert!(result.is_err(), "{}", desc);
    assert!(
        matches!(result.as_ref(), Err(Err::Error(e)) if e.code == expected_kind),
        "Expected {:?} error for {}, got {:?}",
        expected_kind,
        desc,
        result
    );
}

#[rstest]
#[case(b"  1 *", AtomValueEntry { atom_index: 0, value: "*".to_string() })]
#[case(b" 15 query", AtomValueEntry { atom_index: 14, value: "query".to_string() })]
fn test_atom_value_entry(#[case] input: &[u8], #[case] expected: AtomValueEntry) {
    let (remaining, result) = atom_value_entry().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  0 *", "atom index is zero", error::ErrorKind::Verify)]
fn test_atom_value_entry_invalid(#[case] input: &[u8], #[case] desc: &str, #[case] expected_kind: error::ErrorKind) {
    let result = atom_value_entry().parse(input);
    assert!(result.is_err(), "{}", desc);
    assert!(
        matches!(result.as_ref(), Err(Err::Error(e)) if e.code == expected_kind),
        "Expected {:?} error for {}, got {:?}",
        expected_kind,
        desc,
        result
    );
}

#[rstest]
#[case(b"  1   1  -1", vec![ChargeEntry { atom_index: 0, charge: -1 }])]
#[case(b"  2   1  -1   4   1", vec![ChargeEntry { atom_index: 0, charge: -1 }, ChargeEntry { atom_index: 3, charge: 1 }])]
#[case(b"  8   1   1   2   2   3   3   4   4   5   5   6   6   7   7   8   8",
   vec![
        ChargeEntry { atom_index: 0, charge: 1 },
        ChargeEntry { atom_index: 1, charge: 2 },
        ChargeEntry { atom_index: 2, charge: 3 },
        ChargeEntry { atom_index: 3, charge: 4 },
        ChargeEntry { atom_index: 4, charge: 5 },
        ChargeEntry { atom_index: 5, charge: 6 },
        ChargeEntry { atom_index: 6, charge: 7 },
        ChargeEntry { atom_index: 7, charge: 8 },
    ])]
#[case(b"  1  25  15", vec![ChargeEntry { atom_index: 24, charge: 15 }])]
fn test_charge_entries(#[case] input: &[u8], #[case] expected: Vec<ChargeEntry>) {
    let (remaining, result) = charge_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  1   1  -1  a", "trailing characters", error::ErrorKind::Eof)]
#[case(b"  2   1  -1", "count exceeds item list length", error::ErrorKind::Eof)]
#[case(b"  2   1  -1   4   1   5   6", "item list exceeds count", error::ErrorKind::Eof)]
#[case(b"  0", "count is zero", error::ErrorKind::Verify)]
#[case(b"  1   0 -10", "atom index is zero", error::ErrorKind::Verify)]
fn test_charge_entries_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = charge_entries().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"  1   1   2", vec![RadicalEntry { atom_index: 0, radical_type: 2 }])]
#[case(b"  2   1   1   4   3", vec![RadicalEntry { atom_index: 0, radical_type: 1 }, RadicalEntry { atom_index: 3, radical_type: 3 }])]
fn test_radical_entries(#[case] input: &[u8], #[case] expected: Vec<RadicalEntry>) {
    let (remaining, result) = radical_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  1   1   4", "value out of range", error::ErrorKind::Verify)]
#[case(b"  1   1  -1", "value is negative", error::ErrorKind::Digit)]
#[case(b"  2   1   1", "count exceeds item list length", error::ErrorKind::Eof)]
#[case(b"  2   1   1   4   1   5   1", "item list exceeds count", error::ErrorKind::Eof)]
#[case(b"  1   1   2 a", "trailing characters", error::ErrorKind::Eof)]
#[case(b"  0", "count is zero", error::ErrorKind::Verify)]
#[case(b"  1   0   2", "atom index is zero", error::ErrorKind::Verify)]
fn test_radical_entries_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = radical_entries().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"  1   1  13", vec![IsotopeEntry { atom_index: 0, mass: 13 }])]
#[case(b"  2  12   2  15  14", vec![IsotopeEntry { atom_index: 11, mass: 2 }, IsotopeEntry { atom_index: 14, mass: 14 }])]
fn test_isotope_entries(#[case] input: &[u8], #[case] expected: Vec<IsotopeEntry>) {
    let (remaining, result) = isotope_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  1   1  -1", "value is negative", error::ErrorKind::Digit)]
#[case(b"  2   1  10", "count exceeds item list length", error::ErrorKind::Eof)]
#[case(b"  2   1  10   4   1   5   1", "item list exceeds count", error::ErrorKind::Eof)]
#[case(b"  1   1  12 a", "trailing characters", error::ErrorKind::Eof)]
#[case(b"  0", "count is zero", error::ErrorKind::Verify)]
#[case(b"  1   0  12", "atom index is zero", error::ErrorKind::Verify)]
fn test_isotope_entries_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = isotope_entries().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"  1   1   2", vec![RingBondCountEntry { atom_index: 0, ring_bond_count: 2 }])]
#[case(b"  2   1  -1   4  -2", vec![RingBondCountEntry { atom_index: 0, ring_bond_count: -1 }, RingBondCountEntry { atom_index: 3, ring_bond_count: -2 }])]
#[case(b"  1   3   0", vec![RingBondCountEntry { atom_index: 2, ring_bond_count: 0 }])]
#[case(b"  1  10   4", vec![RingBondCountEntry { atom_index: 9, ring_bond_count: 4 }])]
fn test_ring_bond_count_entries(#[case] input: &[u8], #[case] expected: Vec<RingBondCountEntry>) {
    let (remaining, result) = ring_bond_count_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  1   1   5", "value out of range", error::ErrorKind::Verify)]
#[case(b"  1   1  -3", "value out of range", error::ErrorKind::Verify)]
#[case(b"  1   1   2 a", "trailing characters", error::ErrorKind::Eof)]
fn test_ring_bond_count_entries_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = ring_bond_count_entries().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}


#[rstest]
#[case(b"  1   1   3", vec![SubstitutionCountEntry { atom_index: 0, substitution_count: 3 }])]
#[case(b"  2   1  -1   4   6", vec![SubstitutionCountEntry { atom_index: 0, substitution_count: -1 }, SubstitutionCountEntry { atom_index: 3, substitution_count: 6 }])]
#[case(b"  1   5  -2", vec![SubstitutionCountEntry { atom_index: 4, substitution_count: -2 }])]
fn test_substitution_count_entries(#[case] input: &[u8], #[case] expected: Vec<SubstitutionCountEntry>) {
    let (remaining, result) = substitution_count_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  1   1  16", "value out of range", error::ErrorKind::Verify)]
#[case(b"  1   1  -3", "value out of range", error::ErrorKind::Verify)]
#[case(b"  1   1   3 a", "trailing characters", error::ErrorKind::Eof)]
fn test_substitution_count_entries_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = substitution_count_entries().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"  1   1   1", vec![UnsaturatedAtomEntry { atom_index: 0, unsaturated: 1 }])]
#[case(b"  2   1   0   3   1", vec![UnsaturatedAtomEntry { atom_index: 0, unsaturated: 0 }, UnsaturatedAtomEntry { atom_index: 2, unsaturated: 1 }])]
#[case(b"  1  10   0", vec![UnsaturatedAtomEntry { atom_index: 9, unsaturated: 0 }])]
fn test_unsaturated_atom_entries(#[case] input: &[u8], #[case] expected: Vec<UnsaturatedAtomEntry>) {
    let (remaining, result) = unsaturated_atom_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  1   1   2", "value out of range", error::ErrorKind::Verify)]
#[case(b"  1   1  -1", "unsigned value is negative", error::ErrorKind::Digit)]
#[case(b"  1   1   1 a", "trailing characters", error::ErrorKind::Eof)]
fn test_unsaturated_atom_entries_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = unsaturated_atom_entries().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"  1   1   2   5   7", vec![LinkAtomEntry { atom_index: 0, repeat_count: 2, subs_index1: 4, subs_index2: Some(6) }])]
#[case(b"  2   3   3   1   3   8   4   5   6", vec![
    LinkAtomEntry { atom_index: 2, repeat_count: 3, subs_index1: 0, subs_index2: Some(2) },
    LinkAtomEntry { atom_index: 7, repeat_count: 4, subs_index1: 4, subs_index2: Some(5) }
])]
fn test_link_atom_entries(#[case] input: &[u8], #[case] expected: Vec<LinkAtomEntry>) {
    let (remaining, result) = link_atom_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  1   1   1   5   7", "repeat count less than 2", error::ErrorKind::Verify)]
#[case(b"  5   1   2   5   7", "count exceeds 4", error::ErrorKind::Verify)]
#[case(b"  1   1   2   5   7 a", "trailing characters", error::ErrorKind::Eof)]
fn test_link_atom_entries_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = link_atom_entries().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {:?}",
        desc,
        expected_kind,
        result.unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"   1  3 F C   N   O   ", AtomListEntry { atom_index: 0, exclusion: false, elements: vec![Element::C, Element::N, Element::O] })]
#[case(b"   1  3 F C   N   O", AtomListEntry { atom_index: 0, exclusion: false, elements: vec![Element::C, Element::N, Element::O] })]
#[case(b"   5  2 T Cl  Br  ", AtomListEntry { atom_index: 4, exclusion: true, elements: vec![Element::Cl, Element::Br] })]
#[case(b"  10  1   H   ", AtomListEntry { atom_index: 9, exclusion: false, elements: vec![Element::H] })]
fn test_atom_list_entry(#[case] input: &[u8], #[case] expected: AtomListEntry) {
    let (remaining, result) = atom_list_entry().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"   1  0 F C   ", "count is zero", error::ErrorKind::Verify)]
#[case(b"   1 17 F C   N   ", "count exceeds 16", error::ErrorKind::Verify)]
#[case(b"   1  1 X C   ", "invalid exclusion flag", error::ErrorKind::Tag)]
#[case(b"   1  1 F XX  ", "invalid element symbol", error::ErrorKind::MapOpt)]
fn test_atom_list_entry_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = atom_list_entry().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {:?}",
        desc,
        expected_kind,
        result.unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"  1   1   2", vec![AttachmentPointEntry { atom_index: 0, attachment_type: 2 }])]
#[case(b"  2   1   1   2   3", vec![AttachmentPointEntry { atom_index: 0, attachment_type: 1 }, AttachmentPointEntry { atom_index: 1, attachment_type: 3 }])]
fn test_attachment_point_entries(#[case] input: &[u8], #[case] expected: Vec<AttachmentPointEntry>) {
    let (remaining, result) = attachment_point_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  1   1   4", "attachment type out of range", error::ErrorKind::Verify)]
#[case(b"  3   1   1   2   2   3   3", "count exceeds 2", error::ErrorKind::Verify)]
#[case(b"  1   1   1 a", "trailing characters", error::ErrorKind::Eof)]
fn test_attachment_point_entries_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = attachment_point_entries().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {:?}",
        desc,
        expected_kind,
        result.unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"   4  2  14   1   9   2", AtomAttachmentOrderEntry { atom_index: 3, attachments: vec![(13, 1), (8, 2)] })]
fn test_atom_attachment_order_entry(#[case] input: &[u8], #[case] expected: AtomAttachmentOrderEntry) {
    let (remaining, result) = atom_attachment_order_entry().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"   1   0", "count is zero", error::ErrorKind::Verify)]
#[case(b"   1   3   1   2", "count exceeds 2", error::ErrorKind::Verify)]
#[case(b"   0   1   1   2", "atom index is zero", error::ErrorKind::Verify)]
fn test_atom_attachment_order_entry_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = atom_attachment_order_entry().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {:?}",
        desc,
        expected_kind,
        result.unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"  1   1   2", vec![RGroupLabelEntry { atom_index: 0, label: 2 }])]
#[case(b"  2   1   1   2   2", vec![RGroupLabelEntry { atom_index: 0, label: 1 }, RGroupLabelEntry { atom_index: 1, label: 2 }])]
fn test_rgroup_label_entries(#[case] input: &[u8], #[case] expected: Vec<RGroupLabelEntry>) {
    let (remaining, result) = rgroup_label_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  1   0", "label is zero", error::ErrorKind::Verify)]
#[case(b"  9   1   2", "count exceeds 8", error::ErrorKind::Verify)]
#[case(b"  1   0   2", "atom index is zero", error::ErrorKind::Verify)]
#[case(b"  1   1   2 a", "trailing characters", error::ErrorKind::Eof)]
fn test_rgroup_label_entries_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = rgroup_label_entries().parse(input);
    assert!(result.is_err(), "{}", desc);
    assert!(
        matches!(result.as_ref(), Err(Err::Error(e)) if e.code == expected_kind),
        "Expected {:?} error for {}, got {:?}",
        expected_kind,
        desc,
        result
    );
}

#[rstest]
#[case(b"  1   1   0   0  >2", RGroupLogicEntry { label: 1, dependent_label: None, rgroup_or_h: false, occurrence: vec![RGroupOccurrence::GreaterThan(2)] })]
#[case(b"  1   1   0   0  0,>0", RGroupLogicEntry { label: 1, dependent_label: None, rgroup_or_h: false, occurrence: vec![RGroupOccurrence::Exactly(0), RGroupOccurrence::GreaterThan(0)] })]
#[case(b"  1   1   2   0", RGroupLogicEntry { label: 1, dependent_label: Some(2), rgroup_or_h: false, occurrence: vec![RGroupOccurrence::GreaterThan(0)] })]
#[case(b"  1   1   0   1", RGroupLogicEntry { label: 1, dependent_label: None, rgroup_or_h: true, occurrence: vec![RGroupOccurrence::GreaterThan(0)] })]
#[case(b"  1   1   2", RGroupLogicEntry { label: 1, dependent_label: Some(2), rgroup_or_h: false, occurrence: vec![RGroupOccurrence::GreaterThan(0)] })]
fn test_rgroup_logic_entry(#[case] input: &[u8], #[case] expected: RGroupLogicEntry) {
    let (remaining, result) = rgroup_logic_entry().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  0   1   0", "count is zero", error::ErrorKind::Verify)]
#[case(b"  2   1   0", "count exceeds 1", error::ErrorKind::Verify)]
#[case(b"  1   0   0", "label is zero", error::ErrorKind::Verify)]
#[case(b"  1   1   0   2", "rgroup_or_h out of range", error::ErrorKind::Verify)]
fn test_rgroup_logic_entry_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = rgroup_logic_entry().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"  1   1 SUP", vec![SGroupTypeEntry { sgroup_index: 0, sgroup_type: SGroupType::Superatom }])]
#[case(b"  2   1 SUP   2 DAT", vec![
    SGroupTypeEntry { sgroup_index: 0, sgroup_type: SGroupType::Superatom }, 
    SGroupTypeEntry { sgroup_index: 1, sgroup_type: SGroupType::Data }
])]
fn test_sgroup_type_entries(#[case] input: &[u8], #[case] expected: Vec<SGroupTypeEntry>) {
    let (remaining, result) = sgroup_type_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  1   1 FOO", "invalid sgroup type", error::ErrorKind::MapRes)]
#[case(b"  1   1 SUP a", "trailing characters", error::ErrorKind::Eof)]
fn test_sgroup_type_entries_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = sgroup_type_entries().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"  1   1 ALT", vec![SGroupSubtypeEntry { sgroup_index: 0, sgroup_subtype: SGroupSubtype::Alternating }])]
#[case(b"  2   1 RAN   2 BLO", vec![
    SGroupSubtypeEntry { sgroup_index: 0, sgroup_subtype: SGroupSubtype::Random }, 
    SGroupSubtypeEntry { sgroup_index: 1, sgroup_subtype: SGroupSubtype::Block }
])]
fn test_sgroup_subtype_entries(#[case] input: &[u8], #[case] expected: Vec<SGroupSubtypeEntry>) {
    let (remaining, result) = sgroup_subtype_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  1   1 FOO", "invalid sgroup subtype", error::ErrorKind::MapRes)]
#[case(b"  1   1 ALT a", "trailing characters", error::ErrorKind::Eof)]
fn test_sgroup_subtype_entries_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = sgroup_subtype_entries().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}   


#[rstest]
#[case(b"  1   1   1", vec![SGroupLabelEntry { sgroup_index: 0, label: 1 }])]
#[case(b"  2   1  14   2  15", vec![
    SGroupLabelEntry { sgroup_index: 0, label: 14 }, 
    SGroupLabelEntry { sgroup_index: 1, label: 15 }
])]
fn test_sgroup_label_entries(#[case] input: &[u8], #[case] expected: Vec<SGroupLabelEntry>) {
    let (remaining, result) = sgroup_label_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  1   1   0", "label out of range", error::ErrorKind::Verify)]
#[case(b"  1   1 513", "label out of range", error::ErrorKind::Verify)]
#[case(b"  1   1   1 a", "trailing characters", error::ErrorKind::Eof)]
fn test_sgroup_label_entries_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = sgroup_label_entries().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"  3   1 HT    2 HT    3 HT ", vec![
    SGroupConnectivityEntry { sgroup_index: 0, connectivity: SGroupConnectivity::HeadToTail },
    SGroupConnectivityEntry { sgroup_index: 1, connectivity: SGroupConnectivity::HeadToTail },
    SGroupConnectivityEntry { sgroup_index: 2, connectivity: SGroupConnectivity::HeadToTail }
])]
#[case(b"  2   1 HT    2 EU ", vec![
    SGroupConnectivityEntry { sgroup_index: 0, connectivity: SGroupConnectivity::HeadToTail },
    SGroupConnectivityEntry { sgroup_index: 1, connectivity: SGroupConnectivity::EitherUnknown }
])]
fn test_sgroup_connectivity_entries(#[case] input: &[u8], #[case] expected: Vec<SGroupConnectivityEntry>) {
    let (remaining, result) = sgroup_connectivity_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  0", "count is zero", error::ErrorKind::Verify)]
#[case(b"  1   1 FOO", "invalid connectivity", error::ErrorKind::MapRes)]
#[case(b"  1   1 HT a", "trailing characters", error::ErrorKind::Eof)]
fn test_sgroup_connectivity_entries_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = sgroup_connectivity_entries().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b" EXP  1   1", vec![SGroupExpansionEntry { sgroup_index: 0 }])]
#[case(b" EXP  2   1   2", vec![
    SGroupExpansionEntry { sgroup_index: 0 },
    SGroupExpansionEntry { sgroup_index: 1 }
])]
fn test_sgroup_expansion_entries(#[case] input: &[u8], #[case] expected: Vec<SGroupExpansionEntry>) {
    let (remaining, result) = sgroup_expansion_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b" EXP  0   1", "count is zero", error::ErrorKind::Verify)]
#[case(b" EXP  1   1 a", "trailing characters", error::ErrorKind::Eof)]
fn test_sgroup_expansion_entries_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = sgroup_expansion_entries().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"   1  2   1   2", SGroupAtomListEntry { sgroup_index: 0, atom_indices: vec![0, 1] })]
#[case(b"   3  1  15", SGroupAtomListEntry { sgroup_index: 2, atom_indices: vec![14] })]
fn test_sgroup_atom_list_entry(#[case] input: &[u8], #[case] expected: SGroupAtomListEntry) {
    let (remaining, result) = sgroup_atom_list_entry().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"   1 16   1", "count exceeds 15", error::ErrorKind::Verify)]
#[case(b"   1  1   1 a", "trailing characters", error::ErrorKind::Eof)]
fn test_sgroup_atom_list_entry_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = sgroup_atom_list_entry().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"   1  2   1   2", SGroupBondListEntry { sgroup_index: 0, bond_indices: vec![0, 1] })]
#[case(b"   3  1  15", SGroupBondListEntry { sgroup_index: 2, bond_indices: vec![14] })]
fn test_sgroup_bond_list_entry(#[case] input: &[u8], #[case] expected: SGroupBondListEntry) {
    let (remaining, result) = sgroup_bond_list_entry().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"   1 16   1", "count exceeds 15", error::ErrorKind::Verify)]
#[case(b"   1  1   1 a", "trailing characters", error::ErrorKind::Eof)]
fn test_sgroup_bond_list_entry_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = sgroup_bond_list_entry().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"   1  4   3   4   5   6", SGroupParentAtomEntry { sgroup_index: 0, atom_indices: vec![2, 3, 4, 5] })]
fn test_sgroup_parent_atom_entries(#[case] input: &[u8], #[case] expected: SGroupParentAtomEntry) {
    let (remaining, result) = sgroup_parent_atom_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"   1 16   3", "count exceeds 15", error::ErrorKind::Verify)]
#[case(b"   1  4   3   4   5   6 a", "trailing characters", error::ErrorKind::Eof)]
fn test_sgroup_parent_atom_entries_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = sgroup_parent_atom_entries().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"   1 1", SGroupSubscriptEntry { sgroup_index: 0, data: SGroupSubscriptData::Multiplier(SGroupMultiplier::Count(1)) })]
#[case(b"   1 n", SGroupSubscriptEntry { sgroup_index: 0, data: SGroupSubscriptData::Multiplier(SGroupMultiplier::N) })]
#[case(b"   1 Ph", SGroupSubscriptEntry { sgroup_index: 0, data: SGroupSubscriptData::Subscript("Ph".to_string()) })]
fn test_sgroup_subscript_entry(#[case] input: &[u8], #[case] expected: SGroupSubscriptEntry) {
    let (remaining, result) = sgroup_subscript_entry().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"   0 1", "sgroup index is zero", error::ErrorKind::Verify)]
fn test_sgroup_subscript_entry_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = sgroup_subscript_entry().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"   3  3  10   9   4", SGroupCorrespondenceEntry { sgroup_index: 2, bond_indices: vec![9, 8, 3] })]
fn test_sgroup_correspondence_entry(#[case] input: &[u8], #[case] expected: SGroupCorrespondenceEntry) {
    let (remaining, result) = sgroup_correspondence_entry().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"   3  0", "count is zero", error::ErrorKind::Verify)]
#[case(b"   3  3  10   9   4 a", "trailing characters", error::ErrorKind::Eof)]
fn test_sgroup_correspondence_entry_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = sgroup_correspondence_entry().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"   1  4  -13.0153    4.4289  -13.0153    8.2211", SGroupDisplayInfoEntry
       { sgroup_index: 0,  bracket_coords: vec![-13.0153, 4.4289, -13.0153, 8.2211]})]
fn test_sgroup_display_info_entry(#[case] input: &[u8], #[case] expected: SGroupDisplayInfoEntry) {

    let (remaining, result) = sgroup_display_info_entry().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");

    assert_eq!(result.sgroup_index, expected.sgroup_index);
    for i in 0..result.bracket_coords.len() {
        assert!(approx_eq!(f64, result.bracket_coords[i], expected.bracket_coords[i]));
    }
}

#[rstest]
#[case(b"   0  4    4.4700   -3.1700    4.4700   -5.7500", "sgroup index is zero", error::ErrorKind::Verify)]
#[case(b"   1  0", "count is zero", error::ErrorKind::Verify)]
fn test_sgroup_display_info_entry_invalid(#[case] input: &[u8], #[case] desc: &str, #[case] expected_kind: error::ErrorKind) {
    let result = sgroup_display_info_entry().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"   1   6   -0.7200   -0.4200", SGroupConnectingBondEntry { sgroup_index: 0, bond_index: 5, bond_vector: (-0.7200, -0.4200) })]
fn test_sgroup_connecting_bond_entry(#[case] input: &[u8], #[case] expected: SGroupConnectingBondEntry) {
    let (remaining, result) = sgroup_connecting_bond_entry().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result.sgroup_index, expected.sgroup_index);
    assert_eq!(result.bond_index, expected.bond_index);
    assert!(approx_eq!(f64, result.bond_vector.0, expected.bond_vector.0));
    assert!(approx_eq!(f64, result.bond_vector.1, expected.bond_vector.1));
}

#[rstest]
#[case(b"   0   1   -0.7200   -0.4200", "sgroup index is zero", error::ErrorKind::Verify)]
#[case(b"   1   0   -0.7200   -0.4200", "bond index is zero", error::ErrorKind::Verify)]
#[case(b"   1   1   -0.7200   -0.4200 a", "trailing characters", error::ErrorKind::Eof)]
fn test_sgroup_connecting_bond_entry_invalid(#[case] input: &[u8], #[case] desc: &str, #[case] expected_kind: error::ErrorKind) {
    let result = sgroup_connecting_bond_entry().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"   1 pH   ", SGroupDataDescriptionEntry { sgroup_index: 0, field_name: "pH".to_string(), field_type: SGroupDataType::Text, field_units: None, query_identifier: None, data_query_operator: None })]
#[case(b"   3 WEIGHT_PERCENT                N %", SGroupDataDescriptionEntry { sgroup_index: 2, field_name: "WEIGHT_PERCENT".to_string(), field_type: SGroupDataType::Numeric, field_units: Some("%".to_string()), query_identifier: None, data_query_operator: None })]
fn test_sgroup_data_description_entry(#[case] input: &[u8], #[case] expected: SGroupDataDescriptionEntry) {
    let (remaining, result) = sgroup_data_description_entry().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"   0 pH   ", "sgroup index is zero", error::ErrorKind::Verify)]
#[case(b"   1 pH                            X", "invalid field type", error::ErrorKind::MapRes)]
fn test_sgroup_data_description_entry_invalid(#[case] input: &[u8], #[case] desc: &str, #[case] expected_kind: error::ErrorKind) {
    let result = sgroup_data_description_entry().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind), "Mismatched error kind for {}, expected {:?}, got {}", desc, expected_kind, result.clone().unwrap_err().map(|e| e.code));
}

#[rstest]
#[case(b"   1     0.0000    0.0000    DR    ALL  0       0",
    SGroupDataDisplayEntry { sgroup_index: 0, coords: (0.0000, 0.0000), display_type: SGroupDataDisplayType::Detached,
        display_placement: SGroupDataDisplayPlacement::Relative, display_units: SGroupDataDisplayUnits::None,
        display_chars: SGroupDataDisplayChars::All, display_tag: None, display_position: 0 })]
#[case(b"   2     0.0000    0.0000    DR    ALL  1       6",
    SGroupDataDisplayEntry { sgroup_index: 1, coords: (0.0000, 0.0000), display_type: SGroupDataDisplayType::Detached,
        display_placement: SGroupDataDisplayPlacement::Relative, display_units: SGroupDataDisplayUnits::None,
        display_chars: SGroupDataDisplayChars::All, display_tag: None, display_position: 6 })]
#[case(b"   2     0.0000    0.0000    DR    1    1       6",
    SGroupDataDisplayEntry { sgroup_index: 1, coords: (0.0000, 0.0000), display_type: SGroupDataDisplayType::Detached,
        display_placement: SGroupDataDisplayPlacement::Relative, display_units: SGroupDataDisplayUnits::None,
        display_chars: SGroupDataDisplayChars::Number(1), display_tag: None, display_position: 6 })]
fn test_sgroup_data_display_entry(#[case] input: &[u8], #[case] expected: SGroupDataDisplayEntry) {
    let (remaining, result) = sgroup_data_display_entry().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"   0     0.0000    0.0000    DR    ALL  0       0", "sgroup index is zero", error::ErrorKind::Verify)]
#[case(b"   1     0.0000    0.0000    XR    ALL  0       0", "invalid display type", error::ErrorKind::MapRes)]
#[case(b"   1     0.0000    0.0000    DX    ALL  0       0", "invalid placement type", error::ErrorKind::MapRes)]
#[case(b"   1     0.0000    0.0000    DR    NON  0       0", "invalid chars type", error::ErrorKind::MapRes)]
fn test_sgroup_data_display_entry_invalid(
    #[case] input: &[u8], #[case] desc: &str, #[case] expected_kind: error::ErrorKind
) {
    let result = sgroup_data_display_entry().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected: {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rstest]
#[case(b"   1 4.6", SGroupDataEntry::Continuation { sgroup_index: 0, data_content: "4.6".to_string() })]
#[case(b"   2 E/Z unknown", SGroupDataEntry::Continuation { sgroup_index: 1, data_content: "E/Z unknown".to_string() })]
#[case(b"   1", SGroupDataEntry::Continuation { sgroup_index: 0, data_content: "".to_string() })]
fn tests_sgroup_data_continuation_entry(#[case] input: &[u8], #[case] expected: SGroupDataEntry) {
    let (remaining, result) = sgroup_data_continuation_entry().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"   0 4.6", "sgroup index is zero", error::ErrorKind::Verify)]
fn test_sgroup_data_continuation_entry_invalid(#[case] input: &[u8], #[case] desc: &str, #[case] expected_kind: error::ErrorKind) {
    let result = sgroup_data_continuation_entry().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"   1 4.6", SGroupDataEntry::EndWithData { sgroup_index: 0, data_content: "4.6".to_string() })]
#[case(b"   2 E/Z unknown", SGroupDataEntry::EndWithData { sgroup_index: 1, data_content: "E/Z unknown".to_string() })]
#[case(b"   1", SGroupDataEntry::EndBlank { sgroup_index: 0 })]
#[case(b"   2", SGroupDataEntry::EndBlank { sgroup_index: 1 })]
fn tests_sgroup_data_end_entry(#[case] input: &[u8], #[case] expected: SGroupDataEntry) {
    let (remaining, result) = sgroup_data_end_entry().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"   0  1", "sgroup index is zero", error::ErrorKind::Verify)]
fn test_sgroup_data_end_entry_invalid(#[case] input: &[u8], #[case] desc: &str, #[case] expected_kind: error::ErrorKind) {
    let result = sgroup_data_end_entry().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rstest]
#[case(b"  1   1   2", vec![SGroupHierarchyEntry { sgroup_index: 0, parent_sgroup_index: 1 }])]
#[case(b"  3   1   4   2   4   3   2", vec![
    SGroupHierarchyEntry { sgroup_index: 0, parent_sgroup_index: 3 },
    SGroupHierarchyEntry { sgroup_index: 1, parent_sgroup_index: 3 },
    SGroupHierarchyEntry { sgroup_index: 2, parent_sgroup_index: 1 }
])]
fn test_sgroup_hierarchy_entries(#[case] input: &[u8], #[case] expected: Vec<SGroupHierarchyEntry>) {
    let (remaining, result) = sgroup_hierarchy_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  0   1   2", "count is zero", error::ErrorKind::Verify)]
#[case(b"  1   0   2", "sgroup index is zero", error::ErrorKind::Verify)]
#[case(b"  1   1   0", "parent sgroup index is zero", error::ErrorKind::Verify)]
#[case(b"  1   1   2 a", "trailing characters", error::ErrorKind::Eof)]
fn test_sgroup_hierarchy_entries_invalid(#[case] input: &[u8], #[case] desc: &str, #[case] expected_kind: error::ErrorKind) {
    let result = sgroup_hierarchy_entries().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"  2   1   1   2   2", vec![
    SGroupComponentEntry { sgroup_index: 0, component_number: 1 },
    SGroupComponentEntry { sgroup_index: 1, component_number: 2 }
])]
fn test_sgroup_component_entries(#[case] input: &[u8], #[case] expected: Vec<SGroupComponentEntry>) {
    let (remaining, result) = sgroup_component_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  0   1   2", "count is zero", error::ErrorKind::Verify)]
#[case(b"  1   0   2", "component number is zero", error::ErrorKind::Verify)]
#[case(b"  1   1   0 a", "trailing characters", error::ErrorKind::Eof)]
fn test_sgroup_component_entries_invalid(#[case] input: &[u8], #[case] desc: &str, #[case] expected_kind: error::ErrorKind) {
    let result = sgroup_component_entries().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind), "Mismatched error kind for {}, expected {:?}, got {}", desc, expected_kind, result.clone().unwrap_err().map(|e| e.code));
}

#[rstest]
#[case(b"  1   1   0", vec![ZeroBondOrderEntry { bond_index: 0, bond_order: 0 }])]
#[case(b"  2   1   2   3   4", vec![
    ZeroBondOrderEntry { bond_index: 0, bond_order: 2 },
    ZeroBondOrderEntry { bond_index: 2, bond_order: 4 },
])]
fn test_zero_order_bond_entries(#[case] input: &[u8], #[case] expected: Vec<ZeroBondOrderEntry>) {
    let (remaining, result) = zero_bond_order_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  0   1   0", "count is zero", error::ErrorKind::Verify)]
#[case(b"  1   0   0", "bond index is zero", error::ErrorKind::Verify)]
#[case(b"  1   1   0 a", "trailing characters", error::ErrorKind::Eof)]
fn test_zero_bond_order_entries_invalid(#[case] input: &[u8], #[case] desc: &str, #[case] expected_kind: error::ErrorKind) {
    let result = zero_bond_order_entries().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind), "Mismatched error kind for {}, expected {:?}, got {}", desc, expected_kind, result.clone().unwrap_err().map(|e| e.code));
}

#[rstest]
#[case(b"  1   1   0", vec![ZeroAtomChargeEntry { atom_index: 0, charge: 0 }])]
#[case(b"  2   1   2   3  -1", vec![
    ZeroAtomChargeEntry { atom_index: 0, charge: 2 },
    ZeroAtomChargeEntry { atom_index: 2, charge: -1 },
])]
fn test_zero_atom_charge_entries(#[case] input: &[u8], #[case] expected: Vec<ZeroAtomChargeEntry>) {
    let (remaining, result) = zero_atom_charge_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  0   1   0", "count is zero", error::ErrorKind::Verify)]
#[case(b"  1   0   0", "atom index is zero", error::ErrorKind::Verify)]
#[case(b"  1   1   0 a", "trailing characters", error::ErrorKind::Eof)]
fn test_zero_atom_charge_entries_invalid(#[case] input: &[u8], #[case] desc: &str, #[case] expected_kind: error::ErrorKind) {
    let result = zero_atom_charge_entries().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind), "Mismatched error kind for {}, expected {:?}, got {}", desc, expected_kind, result.clone().unwrap_err().map(|e| e.code));
}

#[rstest]
#[case(b"  1   1   0", vec![AtomHydrogenCountEntry { atom_index: 0, hydrogen_count: 0 }])]
#[case(b"  2   1   2   3   4", vec![
    AtomHydrogenCountEntry { atom_index: 0, hydrogen_count: 2 },
    AtomHydrogenCountEntry { atom_index: 2, hydrogen_count: 4 },
])]
fn test_atom_hydrogen_count_entries(#[case] input: &[u8], #[case] expected: Vec<AtomHydrogenCountEntry>) {
    let (remaining, result) = atom_hydrogen_count_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  0   1   0", "count is zero", error::ErrorKind::Verify)]
#[case(b"  1   0   0", "atom index is zero", error::ErrorKind::Verify)]
#[case(b"  1   1   0 a", "trailing characters", error::ErrorKind::Eof)]
fn test_atom_hydrogen_count_entries_invalid(#[case] input: &[u8], #[case] desc: &str, #[case] expected_kind: error::ErrorKind) {
    let result = atom_hydrogen_count_entries().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind), "Mismatched error kind for {}, expected {:?}, got {}", desc, expected_kind, result.clone().unwrap_err().map(|e| e.code));
}