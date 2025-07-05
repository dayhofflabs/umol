use super::*;
use pretty_assertions::assert_eq;
use nom::{error::ErrorKind, Err};
use rstest::rstest;
use umol_data::Element;

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
#[case(b"M  CHG  1   1  -1  a", "trailing chars", ErrorKind::Digit)]
#[case(b"M  CHG  2   1  -1", "count does not match item list", ErrorKind::Digit)]
#[case(
    b"M  CHG  1   1  -1   4   1",
    "item list longer than count",
    ErrorKind::Digit
)]
#[case(b"M  CHG  0", "count is zero", ErrorKind::Digit)]
#[case(b"M  CHG  1   0 -10", "atom index is zero", ErrorKind::Digit)]
#[case(b"M  XXX  1   1  -1", "invalid property tag", ErrorKind::Digit)]
#[case(b"X  CHG  1   1  -1", "invalid prefix", ErrorKind::Digit)]
fn test_charge_entries_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: ErrorKind,
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
#[case(b"M  RAD  1   1   4", "value out of range", ErrorKind::Digit)]
#[case(b"M  RAD  1   1  -1", "value out of range", ErrorKind::Digit)]
fn test_radical_entries_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: ErrorKind,
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
#[case(b"M  ISO  1   1  13", "value out of range", ErrorKind::Digit)]
#[case(b"M  ISO  1   1  -1", "value out of range", ErrorKind::Digit)]
fn test_isotope_entries_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: ErrorKind,
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
#[case(b"  1  2   1   2", SGroupAtomListEntry { sgroup_index: 0, atom_indices: vec![0, 1] })]
#[case(b"  3  1  15", SGroupAtomListEntry { sgroup_index: 2, atom_indices: vec![14] })]
fn test_sgroup_atom_list_entry(#[case] input: &[u8], #[case] expected: SGroupAtomListEntry) {
    let (remaining, result) = sgroup_atom_list_entry().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  1  2   1   2", SGroupBondListEntry { sgroup_index: 0, bond_indices: vec![0, 1] })]
#[case(b"  3  1  15", SGroupBondListEntry { sgroup_index: 2, bond_indices: vec![14] })]
fn test_sgroup_bond_list_entry(#[case] input: &[u8], #[case] expected: SGroupBondListEntry) {
    let (remaining, result) = sgroup_bond_list_entry().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  1 CF3", AtomAliasEntry { atom_index: 0, alias: "CF3".to_string() })]
#[case(b" 15 Et", AtomAliasEntry { atom_index: 14, alias: "Et".to_string() })]
fn test_atom_alias_entry(#[case] input: &[u8], #[case] expected: AtomAliasEntry) {
    let (remaining, result) = atom_alias_entry().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
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
#[case(b"M  CHG  1   1  -1", "CHG standard property", PropertyEntries::ChargeEntries(vec![ChargeEntry { atom_index: 0, charge: -1 }]))]
#[case(b"M  RAD  1   1   2", "RAD standard property", PropertyEntries::RadicalEntries(vec![RadicalEntry { atom_index: 0, radical_type: 2 }]))]
#[case(b"M  ISO  1   1  13", "ISO standard property", PropertyEntries::IsotopeEntries(vec![IsotopeEntry { atom_index: 0, mass: 13 }]))]
#[case(b"M  STY  1   1 SUP", "STY SGroup property", PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry { sgroup_index: 0, sgroup_type: SGroupType::Superatom }]))]
#[case(b"M  SLB  1   1  19", "SLB SGroup property", PropertyEntries::SGroupLabelEntries(vec![SGroupLabelEntry { sgroup_index: 0, label: 19 }]))]
#[case(b"M  SAL  1  1   5", "SAL SGroup property", PropertyEntries::SGroupAtomListEntry(SGroupAtomListEntry { sgroup_index: 0, atom_indices: vec![4] }))]
#[case(b"M  SBL  1  1   3", "SBL SGroup property", PropertyEntries::SGroupBondListEntry(SGroupBondListEntry { sgroup_index: 0, bond_indices: vec![2] }))]
#[case(b"A    1 CF3", "A atom alias property", PropertyEntries::AtomAliasEntry(AtomAliasEntry { atom_index: 0, alias: "CF3".to_string() }))]
#[case(b"V    1 *", "V atom value property", PropertyEntries::AtomValueEntry(AtomValueEntry { atom_index: 0, value: "*".to_string() }))]
#[case(b"M  ALS  1  3FC   N   O   ", "ALS query property", PropertyEntries::AtomListEntry(AtomListEntry { atom_index: 0, exclusion: false, elements: vec![Element::C, Element::N, Element::O] }))]
#[case(b"M  APO  1   1   1", "APO query property", PropertyEntries::AttachmentPointEntries(vec![AttachmentPointEntry { atom_index: 0, attachment_type: 1 }]))]
#[case(b"M  AAL  1 1   2   1", "AAL query property", PropertyEntries::AtomAttachmentOrderEntry(AtomAttachmentOrderEntry { atom_index: 0, attachments: vec![(1, 1)] }))]
#[case(b"M  RBC  1   1   2", "RBC query property", PropertyEntries::RingBondCountEntries(vec![RingBondCountEntry { atom_index: 0, ring_bond_count: 2 }]))]
#[case(b"M  SUB  1   1   3", "SUB query property", PropertyEntries::SubstitutionCountEntries(vec![SubstitutionCountEntry { atom_index: 0, substitution_count: 3 }]))]
#[case(b"M  UNS  1   1   1", "UNS query property", PropertyEntries::UnsaturatedAtomEntries(vec![UnsaturatedAtomEntry { atom_index: 0, unsaturated: 1 }]))]
#[case(b"M  LIN  1   1   2   5   7", "LIN query property", PropertyEntries::LinkAtomEntries(vec![LinkAtomEntry { atom_index: 0, repeat_count: 2, bond1: 5, bond2: 7 }]))]
fn test_property_input(#[case] input: &[u8], #[case] desc: &str, #[case] expected: PropertyEntries) {
    let (remaining, result) = property_input(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty for {}", desc);
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"M  CHG  1   1  -1", "CHG standard property", PropertyEntries::ChargeEntries(vec![ChargeEntry { atom_index: 0, charge: -1 }]))]
#[case(b"M  RAD  1   1   2", "RAD standard property", PropertyEntries::RadicalEntries(vec![RadicalEntry { atom_index: 0, radical_type: 2 }]))]
#[case(b"M  ISO  1   1  13", "ISO standard property", PropertyEntries::IsotopeEntries(vec![IsotopeEntry { atom_index: 0, mass: 13 }]))]
#[case(b"M  STY  1   1 SUP", "STY SGroup property", PropertyEntries::SGroupTypeEntries(vec![SGroupTypeEntry { sgroup_index: 0, sgroup_type: SGroupType::Superatom }]))]
#[case(b"M  SLB  1   1  19", "SLB SGroup property", PropertyEntries::SGroupLabelEntries(vec![SGroupLabelEntry { sgroup_index: 0, label: 19 }]))]
#[case(b"M  SAL  1  1   5", "SAL SGroup property", PropertyEntries::SGroupAtomListEntry(SGroupAtomListEntry { sgroup_index: 0, atom_indices: vec![4] }))]
#[case(b"M  SBL  1  1   3", "SBL SGroup property", PropertyEntries::SGroupBondListEntry(SGroupBondListEntry { sgroup_index: 0, bond_indices: vec![2] }))]
#[case(b"A    1 CF3", "A atom alias property", PropertyEntries::AtomAliasEntry(AtomAliasEntry { atom_index: 0, alias: "CF3".to_string() }))]
#[case(b"V    1 *", "V atom value property", PropertyEntries::AtomValueEntry(AtomValueEntry { atom_index: 0, value: "*".to_string() }))]
fn test_property_input_standard(#[case] input: &[u8], #[case] desc: &str, #[case] expected: PropertyEntries) {
    let (remaining, result) = property_input_standard(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty for {}", desc);
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"M  ALS  1  3FC   N   O   ", "ALS query property not supported in standard parser")]
#[case(b"M  APO  1   1   1", "APO query property not supported in standard parser")]
#[case(b"M  AAL  1 1   2   1", "AAL query property not supported in standard parser")]
#[case(b"M  RBC  1   1   2", "RBC query property not supported in standard parser")]
#[case(b"M  SUB  1   1   3", "SUB query property not supported in standard parser")]
#[case(b"M  UNS  1   1   1", "UNS query property not supported in standard parser")]
#[case(b"M  LIN  1   1   2   5   7", "LIN query property not supported in standard parser")]
fn test_property_input_standard_invalid(#[case] input: &[u8], #[case] desc: &str) {
    let result = property_input_standard(input);
    assert!(result.is_err(), "{}", desc);
    assert!(
        matches!(result.as_ref(), Err(nom::Err::Error(e)) if e.code == ErrorKind::Tag),
        "Expected Tag error for {}, got {:?}",
        desc,
        result
    );
}

// Individual parser tests for new atom properties

#[rstest]
#[case(b"  1  3FC   N   O   ", AtomListEntry { atom_index: 0, exclusion: false, elements: vec![Element::C, Element::N, Element::O] })]
#[case(b"  5  2TCl  Br  ", AtomListEntry { atom_index: 4, exclusion: true, elements: vec![Element::Cl, Element::Br] })]
#[case(b" 10  1FH   ", AtomListEntry { atom_index: 9, exclusion: false, elements: vec![Element::H] })]
fn test_atom_list_entry(#[case] input: &[u8], #[case] expected: AtomListEntry) {
    let (remaining, result) = atom_list_entry().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  1   1   1", vec![AttachmentPointEntry { atom_index: 0, attachment_type: 1 }])]
#[case(b"  2   1   2   3   3", vec![AttachmentPointEntry { atom_index: 0, attachment_type: 2 }, AttachmentPointEntry { atom_index: 2, attachment_type: 3 }])]
#[case(b"  1   5   0", vec![AttachmentPointEntry { atom_index: 4, attachment_type: 0 }])]
fn test_attachment_point_entries(#[case] input: &[u8], #[case] expected: Vec<AttachmentPointEntry>) {
    let (remaining, result) = attachment_point_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"  1 1   2   1", AtomAttachmentOrderEntry { atom_index: 0, attachments: vec![(1, 1)] })]
#[case(b"  5 2   3   1   7   2", AtomAttachmentOrderEntry { atom_index: 4, attachments: vec![(2, 1), (6, 2)] })]
#[case(b" 10 1   4   2", AtomAttachmentOrderEntry { atom_index: 9, attachments: vec![(3, 2)] })]
fn test_atom_attachment_order_entry(#[case] input: &[u8], #[case] expected: AtomAttachmentOrderEntry) {
    let (remaining, result) = atom_attachment_order_entry().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
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
#[case(b"  1   1   3", vec![SubstitutionCountEntry { atom_index: 0, substitution_count: 3 }])]
#[case(b"  2   1  -1   4   6", vec![SubstitutionCountEntry { atom_index: 0, substitution_count: -1 }, SubstitutionCountEntry { atom_index: 3, substitution_count: 6 }])]
#[case(b"  1   5  -2", vec![SubstitutionCountEntry { atom_index: 4, substitution_count: -2 }])]
fn test_substitution_count_entries(#[case] input: &[u8], #[case] expected: Vec<SubstitutionCountEntry>) {
    let (remaining, result) = substitution_count_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
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
#[case(b"  1   1   2   5   7", vec![LinkAtomEntry { atom_index: 0, repeat_count: 2, bond1: 5, bond2: 7 }])]
#[case(b"  2   3   3   0   0   8   4   1   2", vec![
    LinkAtomEntry { atom_index: 2, repeat_count: 3, bond1: 0, bond2: 0 },
    LinkAtomEntry { atom_index: 7, repeat_count: 4, bond1: 1, bond2: 2 }
])]
fn test_link_atom_entries(#[case] input: &[u8], #[case] expected: Vec<LinkAtomEntry>) {
    let (remaining, result) = link_atom_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

// Error handling tests for new properties

#[rstest]
#[case(b"  1  0FC   ", "count is zero", ErrorKind::Verify)]
#[case(b"  1 17FC   N   ", "count exceeds 16", ErrorKind::Verify)]
#[case(b"  1  1XC   ", "invalid exclusion flag", ErrorKind::Tag)]
#[case(b"  1  1FXX  ", "invalid element symbol", ErrorKind::MapRes)]
fn test_atom_list_entry_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: ErrorKind,
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
#[case(b"  0   1   1", "count is zero", ErrorKind::Verify)]
#[case(b"  3   1   1", "count exceeds 2", ErrorKind::Verify)]
#[case(b"  1   1   4", "attachment type out of range", ErrorKind::Verify)]
fn test_attachment_point_entries_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: ErrorKind,
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
#[case(b"  1 0   2   1", "pair count is zero", ErrorKind::Verify)]
#[case(b"  1 3   2   1", "pair count exceeds 2", ErrorKind::Verify)]
#[case(b"  1 1   2   0", "order is zero", ErrorKind::Verify)]
#[case(b"  1 1   2   3", "order exceeds 2", ErrorKind::Verify)]
fn test_atom_attachment_order_entry_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: ErrorKind,
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
#[case(b"  1   1   1", "repeat count less than 2", ErrorKind::Verify)]
#[case(b"  5   1   2", "count exceeds 4", ErrorKind::Verify)]
fn test_link_atom_entries_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: ErrorKind,
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
