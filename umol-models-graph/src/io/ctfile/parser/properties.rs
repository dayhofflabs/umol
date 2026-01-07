//! Parsers for CTab property lines.

use bstr::ByteSlice;
use nom::branch::alt;
use nom::bytes::complete::{tag, take};
use nom::character::complete::{space0, u8 as nom_u8};
use nom::combinator::{all_consuming, cond, map, map_opt, map_parser, map_res, opt, rest};
use nom::error::{Error as NomError, ErrorKind as NomErrorKind};
use nom::multi::{count as nom_count, length_count};
use nom::sequence::{delimited, preceded, terminated};
use nom::{Err, Parser};
use umol_data::Element;

use super::sgroup::{sgroup_connectivity, sgroup_subtype, sgroup_type};
use super::utils::{
    fixed_width_element_partial, fixed_width_float, fixed_width_int, fixed_width_int_in_range,
    fixed_width_int_minus1, fixed_width_str_partial, fixed_width_unused, rgroup_occurrences,
    LinesWithOffsetExt,
};
use crate::io::ctfile::config::CtabParseFlags;
use crate::io::ctfile::error::ParseError;
use crate::io::ctfile::parser::sgroup::{
    sgroup_data_display_chars, sgroup_data_display_placement, sgroup_data_display_type,
    sgroup_data_display_units, sgroup_data_type, sgroup_multiplier, sgroup_subscript,
};
use crate::table_ir::{
    BondOrder, RGroupOccurrence, SGroupConnectivity, SGroupDataDisplayChars,
    SGroupDataDisplayPlacement, SGroupDataDisplayType, SGroupDataDisplayUnits, SGroupDataType,
    SGroupMultiplier, SGroupSubtype, SGroupType,
};

#[derive(Debug, Clone, PartialEq)]
pub struct MoleculeChiralFlagEntry {
    pub chiral_flag: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtomAliasEntry {
    pub atom_index: usize,
    pub alias: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct LegacyGroupAbbreviationEntry {
    pub atom_index1: usize, // Atoms on this side are abbreviated
    pub atom_index2: usize, // Attachment point to main structure
    pub label: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtomValueEntry {
    pub atom_index: usize,
    pub value: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ChargeEntry {
    pub atom_index: usize,
    pub charge: i8,
}

#[derive(Debug, Clone, PartialEq)]
pub struct RadicalEntry {
    pub atom_index: usize,
    pub radical_type: u8,
}

#[derive(Debug, Clone, PartialEq)]
pub struct IsotopeEntry {
    pub atom_index: usize,
    pub mass: u32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct RingBondCountEntry {
    pub atom_index: usize,
    pub ring_bond_count: i8, // -2=r*, -1=r0, 0=off, 2=r2, 3=r3, 4+=r4
}

#[derive(Debug, Clone, PartialEq)]
pub struct SubstitutionCountEntry {
    pub atom_index: usize,
    pub substitution_count: i8, // -2=s*, -1=s0, 0=off, 1-5=s1-s5, 6+=s6
}

#[derive(Debug, Clone, PartialEq)]
pub struct UnsaturatedAtomEntry {
    pub atom_index: usize,
    pub unsaturated: u8, // 0=off, 1=on
}

#[derive(Debug, Clone, PartialEq)]
pub struct LinkAtomEntry {
    pub atom_index: usize,
    pub repeat_count: u8,           // vvv >= 2
    pub subs_index1: usize,         // bbb
    pub subs_index2: Option<usize>, // ccc (optional)
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtomListEntry {
    pub atom_index: usize,
    pub exclusion: bool,        // T = NOT list, F = normal list
    pub elements: Vec<Element>, // Converted from 4-char symbols
}

#[derive(Debug, Clone, PartialEq)]
pub struct AttachmentPointEntry {
    pub atom_index: usize,
    pub attachment_type: u8, // 0=none, 1=first, 2=second, 3=both
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtomAttachmentOrderEntry {
    pub atom_index: usize,
    pub attachments: Vec<(usize, u8)>, // (neighbor index, attachment order)
}

#[derive(Debug, Clone, PartialEq)]
pub struct RGroupLabelEntry {
    pub atom_index: usize,
    pub label: u32, // 0 = no label, 1-32 in CTab spec, 1-999 in RDKit
}

#[derive(Debug, Clone, PartialEq)]
pub struct RGroupLogicEntry {
    pub label: u32,                   // 0 = no label, 1-32 in CTab spec
    pub dependent_label: Option<u32>, // None = no dependent label
    pub rgroup_or_h: bool,            // false = off, true = on: RGroup or H atom
    pub occurrence: Vec<RGroupOccurrence>,
    // n=exactly n, n-m=from n through m (inclusive), >n=greater n,
    // <n=fewer than n, blank (default): > 0.
    // Any non-contradictory combination of the preceding values is allowed,
    //separated by commas.
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupTypeEntry {
    pub sgroup_index: usize,
    pub sgroup_type: SGroupType,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupSubtypeEntry {
    pub sgroup_index: usize,
    pub sgroup_subtype: SGroupSubtype,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupLabelEntry {
    pub sgroup_index: usize,
    pub label: u32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupConnectivityEntry {
    pub sgroup_index: usize,
    pub connectivity: SGroupConnectivity,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupExpansionEntry {
    pub sgroup_index: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupAtomListEntry {
    pub sgroup_index: usize,
    pub atom_indices: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupBondListEntry {
    pub sgroup_index: usize,
    pub bond_indices: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupParentAtomEntry {
    pub sgroup_index: usize,
    pub atom_indices: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupSubscriptEntry {
    pub sgroup_index: usize,
    pub multiplier: Option<SGroupMultiplier>,
    pub subscript: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupCorrespondenceEntry {
    pub sgroup_index: usize,
    pub bond_indices: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupDisplayInfoEntry {
    pub sgroup_index: usize,
    pub bracket_coords: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupConnectingBondEntry {
    pub sgroup_index: usize,
    pub bond_index: usize,
    pub bond_vector: (f64, f64),
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupDataDescriptionEntry {
    pub sgroup_index: usize,
    pub field_name: String,
    pub field_type: SGroupDataType,
    pub field_units: Option<String>,
    pub query_identifier: Option<String>,
    pub data_query_operator: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum SGroupDataEntry {
    /// SCD - Continuation data (exactly 69 characters)
    Continuation {
        sgroup_index: usize,
        data_content: String, // Should be exactly 69 chars
    },
    /// SED - End with data (<= 69 characters)
    EndWithData {
        sgroup_index: usize,
        data_content: String, // 0-69 chars
    },
    /// SED - End without data (blank, processes buffered data)
    EndBlank { sgroup_index: usize },
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupDataDisplayEntry {
    pub sgroup_index: usize,
    pub coords: (f64, f64),
    pub display_type: SGroupDataDisplayType,
    pub display_placement: SGroupDataDisplayPlacement,
    pub display_units: SGroupDataDisplayUnits,
    pub display_chars: SGroupDataDisplayChars,
    pub display_tag: Option<u8>, // 0 = no tag, 1-9 = tag
    pub display_position: u8,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupComponentEntry {
    pub sgroup_index: usize,
    pub component_number: u32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupHierarchyEntry {
    pub sgroup_index: usize,
    pub parent_sgroup_index: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BondOrderOverrideEntry {
    pub bond_index: usize,
    pub bond_order: BondOrder,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtomChargeOverrideEntry {
    pub atom_index: usize,
    pub charge: i8,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtomHydrogenCountEntry {
    pub atom_index: usize,
    pub hydrogen_count: Option<u8>,
}

/// Parsed property entries
#[derive(Debug, Clone, PartialEq)]
pub enum PropertyEntries {
    MoleculeChiralFlagEntry(MoleculeChiralFlagEntry),
    AtomAliasEntry(AtomAliasEntry),
    LegacyGroupAbbreviationEntry(LegacyGroupAbbreviationEntry),
    AtomValueEntry(AtomValueEntry),
    ChargeEntries(Vec<ChargeEntry>),
    RadicalEntries(Vec<RadicalEntry>),
    IsotopeEntries(Vec<IsotopeEntry>),
    RingBondCountEntries(Vec<RingBondCountEntry>),
    SubstitutionCountEntries(Vec<SubstitutionCountEntry>),
    UnsaturatedAtomEntries(Vec<UnsaturatedAtomEntry>),
    LinkAtomEntries(Vec<LinkAtomEntry>),
    AtomListEntry(AtomListEntry),
    AttachmentPointEntries(Vec<AttachmentPointEntry>),
    AtomAttachmentOrderEntry(AtomAttachmentOrderEntry),
    RGroupLabelEntries(Vec<RGroupLabelEntry>),
    RGroupLogicEntry(RGroupLogicEntry),
    SGroupTypeEntries(Vec<SGroupTypeEntry>),
    SGroupSubtypeEntries(Vec<SGroupSubtypeEntry>),
    SGroupLabelEntries(Vec<SGroupLabelEntry>),
    SGroupConnectivityEntries(Vec<SGroupConnectivityEntry>),
    SGroupExpansionEntries(Vec<SGroupExpansionEntry>),
    SGroupAtomListEntry(SGroupAtomListEntry),
    SGroupBondListEntry(SGroupBondListEntry),
    SGroupParentAtomEntry(SGroupParentAtomEntry),
    SGroupSubscriptEntry(SGroupSubscriptEntry),
    SGroupCorrespondenceEntry(SGroupCorrespondenceEntry),
    SGroupDisplayInfoEntry(SGroupDisplayInfoEntry),
    SGroupConnectingBondEntry(SGroupConnectingBondEntry),
    SGroupDataDescriptionEntry(SGroupDataDescriptionEntry),
    SGroupDataEntry(SGroupDataEntry),
    SGroupDataDisplayEntry(SGroupDataDisplayEntry),
    SGroupHierarchyEntries(Vec<SGroupHierarchyEntry>),
    SGroupComponentEntries(Vec<SGroupComponentEntry>),
    BondOrderOverrideEntries(Vec<BondOrderOverrideEntry>),
    AtomChargeOverrideEntries(Vec<AtomChargeOverrideEntry>),
    AtomHydrogenCountEntries(Vec<AtomHydrogenCountEntry>),
}

/// Parse properties block (basic properties only)
pub(super) fn properties_block<'inp>(
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (Vec<PropertyEntries>, u32), Error = ParseError> + use<'inp> {
    let no_v2000_end_tags = flags.contains(CtabParseFlags::NO_V2000_END_TAGS);

    move |input: &'inp [u8]| {
        let mut properties = Vec::new();
        let mut lines_iter = input.lines_with_offset();
        let mut byte_offset = 0;
        let mut line_index = 0;
        let mut end_found = false;

        while let Some((line, byte_len)) = lines_iter.next() {
            if line.starts_with(b"M  END") {
                byte_offset += byte_len;
                line_index += 1;
                end_found = true;
                break;
            }

            // Handle atom alias (two-line property)
            if line.starts_with(b"A  ") {
                let (next_line, next_byte_len) = lines_iter.next().ok_or_else(|| {
                    Err::Error(ParseError::UnexpectedEof {
                        line: line_offset + line_index + 1,
                        block: "atom alias",
                    })
                })?;
                let property = parse_atom_alias_input(line, next_line).map_err(|_| {
                    Err::Error(ParseError::InvalidPropertyLine {
                        line: line_offset + line_index,
                        col: 0,
                    })
                })?;
                properties.push(property);
                byte_offset += byte_len + next_byte_len;
                line_index += 2;
            } else {
                let (_, property) = all_consuming(terminated(property_input(flags), space0))
                    .parse(line)
                    .map_err(|_| {
                        Err::Error(ParseError::InvalidPropertyLine {
                            line: line_offset + line_index,
                            col: 0,
                        })
                    })?;
                properties.push(property);
                byte_offset += byte_len;
                line_index += 1;
            }
        }

        if !end_found && !no_v2000_end_tags {
            return Err(Err::Error(ParseError::MissingMEndTag {
                line: line_offset + line_index,
            }));
        }

        let remaining = &input[byte_offset..];
        Ok((remaining, (properties, line_offset + line_index)))
    }
}

/// Parse extended properties block
pub(super) fn extended_properties_block<'inp>(
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (Vec<PropertyEntries>, u32), Error = ParseError> + use<'inp> {
    let no_v2000_end_tags = flags.contains(CtabParseFlags::NO_V2000_END_TAGS);

    move |input: &'inp [u8]| {
        let mut properties = Vec::new();
        let mut lines_iter = input.lines_with_offset();
        let mut byte_offset = 0;
        let mut line_index = 0;
        let mut end_found = false;

        while let Some((line, byte_len)) = lines_iter.next() {
            if line.starts_with(b"M  END") {
                byte_offset += byte_len;
                line_index += 1;
                end_found = true;
                break;
            }

            // Handle atom alias (two-line property)
            if line.starts_with(b"A  ") {
                let (next_line, next_byte_len) = lines_iter.next().ok_or_else(|| {
                    Err::Error(ParseError::UnexpectedEof {
                        line: line_offset + line_index + 1,
                        block: "atom alias",
                    })
                })?;
                let property = parse_atom_alias_input(line, next_line).map_err(|_| {
                    Err::Error(ParseError::InvalidPropertyLine {
                        line: line_offset + line_index,
                        col: 0,
                    })
                })?;
                properties.push(property);
                byte_offset += byte_len + next_byte_len;
                line_index += 2;
            // Handle legacy group abbreviation (two-line property)
            } else if line.starts_with(b"G  ") {
                let (next_line, next_byte_len) = lines_iter.next().ok_or_else(|| {
                    Err::Error(ParseError::UnexpectedEof {
                        line: line_offset + line_index + 1,
                        block: "legacy group abbreviation",
                    })
                })?;
                let property =
                    parse_legacy_group_abbreviation_input(line, next_line).map_err(|_| {
                        Err::Error(ParseError::InvalidPropertyLine {
                            line: line_offset + line_index,
                            col: 0,
                        })
                    })?;
                properties.push(property);
                byte_offset += byte_len + next_byte_len;
                line_index += 2;
            } else {
                let (_, property) =
                    all_consuming(terminated(extended_property_input(flags), space0))
                        .parse(line)
                        .map_err(|_| {
                            Err::Error(ParseError::InvalidPropertyLine {
                                line: line_offset + line_index,
                                col: 0,
                            })
                        })?;
                properties.push(property);
                byte_offset += byte_len;
                line_index += 1;
            }
        }

        if !end_found && !no_v2000_end_tags {
            return Err(Err::Error(ParseError::MissingMEndTag {
                line: line_offset + line_index,
            }));
        }

        let remaining = &input[byte_offset..];
        Ok((remaining, (properties, line_offset + line_index)))
    }
}

/// Parse property line (basic properties only)
/// Note: A  (atom alias) and M  END are handled at the block level, not here.
pub fn property_input<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = PropertyEntries, Error = NomError<&'inp [u8]>> + use<'inp> {
    let allow_clark_extensions = flags.contains(CtabParseFlags::CLARK_EXTENSIONS);
    move |input: &'inp [u8]| {
        if input.len() < 3 {
            return Err(Err::Error(NomError::new(input, NomErrorKind::Eof)));
        }

        // A  (atom alias) and M  END are handled at block level
        debug_assert!(
            &input[0..3] != b"A  ",
            "A   should be handled at block level"
        );

        // Handle V lines
        if &input[0..3] == b"V  " {
            return atom_value_entry()
                .parse(&input[3..])
                .map(|(i, o)| (i, PropertyEntries::AtomValueEntry(o)));
        }

        // Handle M lines
        if input.len() < 6 {
            return Err(Err::Error(NomError::new(input, NomErrorKind::Eof)));
        }

        debug_assert!(
            &input[0..6] != b"M  END",
            "M  END should be handled at block level"
        );
        let (remaining, tag_bytes) = take(6u8)(input)?;

        // General M properties
        let (remaining, result) = match tag_bytes {
            b"M  CHG" => charge_entries()
                .parse(remaining)
                .map(|(i, o)| (i, PropertyEntries::ChargeEntries(o))),
            b"M  RAD" => radical_entries()
                .parse(remaining)
                .map(|(i, o)| (i, PropertyEntries::RadicalEntries(o))),
            b"M  ISO" => isotope_entries()
                .parse(remaining)
                .map(|(i, o)| (i, PropertyEntries::IsotopeEntries(o))),
            // Clark extensions
            tag @ (b"M  ZBO" | b"M  ZCH" | b"M  HYD") => {
                if allow_clark_extensions {
                    match tag {
                        b"M  ZBO" => bond_order_override_entries()
                            .parse(remaining)
                            .map(|(i, o)| (i, PropertyEntries::BondOrderOverrideEntries(o))),
                        b"M  ZCH" => atom_charge_overrides_entries()
                            .parse(remaining)
                            .map(|(i, o)| (i, PropertyEntries::AtomChargeOverrideEntries(o))),
                        b"M  HYD" => atom_hydrogen_count_entries()
                            .parse(remaining)
                            .map(|(i, o)| (i, PropertyEntries::AtomHydrogenCountEntries(o))),
                        _ => unreachable!(),
                    }
                } else {
                    Err(Err::Error(NomError::new(input, NomErrorKind::Tag)))
                }
            }
            _ => Err(Err::Error(NomError::new(input, NomErrorKind::Tag))),
        }?;
        Ok((remaining, result))
    }
}

/// Parse extended property line
/// Note: A  (atom alias) and M  END are handled at the block level, not here.
pub fn extended_property_input<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = PropertyEntries, Error = NomError<&'inp [u8]>> + use<'inp> {
    let allow_queries = flags.contains(CtabParseFlags::WILDCARDS);
    let allow_rgroups = flags.contains(CtabParseFlags::RGROUPS);
    let allow_sgroups = flags.contains(CtabParseFlags::SGROUPS);
    let allow_clark_extensions = flags.contains(CtabParseFlags::CLARK_EXTENSIONS);
    let skip_unused_fields = flags.contains(CtabParseFlags::SKIP_UNUSED_FIELDS);
    move |input: &'inp [u8]| {
        if input.len() < 3 {
            return Err(Err::Error(NomError::new(input, NomErrorKind::Eof)));
        }

        // A  (atom alias), G  (legacy group abbreviation), and M  END are handled at block level
        debug_assert!(
            &input[0..3] != b"A  " && &input[0..3] != b"G  ",
            "A or G should be handled at block level"
        );

        // Handle V lines
        if &input[0..3] == b"V  " {
            return atom_value_entry()
                .parse(&input[3..])
                .map(|(i, o)| (i, PropertyEntries::AtomValueEntry(o)));
        }

        // Handle M lines
        if input.len() < 6 {
            return Err(Err::Error(NomError::new(input, NomErrorKind::Eof)));
        }

        debug_assert!(
            &input[0..6] != b"M  END",
            "M  END should be handled at block level"
        );
        let (remaining, tag_bytes) = take(6u8)(input)?;

        let (remaining, result) = match tag_bytes {
            // General M properties
            b"M  CHG" => charge_entries()
                .parse(remaining)
                .map(|(i, o)| (i, PropertyEntries::ChargeEntries(o))),
            b"M  RAD" => radical_entries()
                .parse(remaining)
                .map(|(i, o)| (i, PropertyEntries::RadicalEntries(o))),
            b"M  ISO" => isotope_entries()
                .parse(remaining)
                .map(|(i, o)| (i, PropertyEntries::IsotopeEntries(o))),
            // Query properties
            tag @ (b"M  RBC" | b"M  SUB" | b"M  UNS" | b"M  LIN" | b"M  ALS") if allow_queries => {
                match tag {
                    b"M  RBC" => ring_bond_count_entries()
                        .parse(remaining)
                        .map(|(i, o)| (i, PropertyEntries::RingBondCountEntries(o))),
                    b"M  SUB" => substitution_count_entries()
                        .parse(remaining)
                        .map(|(i, o)| (i, PropertyEntries::SubstitutionCountEntries(o))),
                    b"M  UNS" => unsaturated_atom_entries()
                        .parse(remaining)
                        .map(|(i, o)| (i, PropertyEntries::UnsaturatedAtomEntries(o))),
                    b"M  LIN" => link_atom_entries()
                        .parse(remaining)
                        .map(|(i, o)| (i, PropertyEntries::LinkAtomEntries(o))),
                    b"M  ALS" => atom_list_entry()
                        .parse(remaining)
                        .map(|(i, o)| (i, PropertyEntries::AtomListEntry(o))),
                    _ => unreachable!(),
                }
            }
            // RGroup properties
            tag @ (b"M  APO" | b"M  AAL" | b"M  RGP" | b"M  LOG") if allow_rgroups => match tag {
                b"M  APO" => attachment_point_entries()
                    .parse(remaining)
                    .map(|(i, o)| (i, PropertyEntries::AttachmentPointEntries(o))),
                b"M  AAL" => atom_attachment_order_entry()
                    .parse(remaining)
                    .map(|(i, o)| (i, PropertyEntries::AtomAttachmentOrderEntry(o))),
                b"M  RGP" => rgroup_label_entries()
                    .parse(remaining)
                    .map(|(i, o)| (i, PropertyEntries::RGroupLabelEntries(o))),
                b"M  LOG" => rgroup_logic_entry()
                    .parse(remaining)
                    .map(|(i, o)| (i, PropertyEntries::RGroupLogicEntry(o))),
                _ => unreachable!(),
            },
            // SGroup properties
            tag @ (b"M  STY" | b"M  SST" | b"M  SLB" | b"M  SCN" | b"M  SAL" | b"M  SBL"
            | b"M  SMT" | b"M  SDS" | b"M  SPA" | b"M  CRS" | b"M  SDI" | b"M  SBV"
            | b"M  SDT" | b"M  SDD" | b"M  SCD" | b"M  SED" | b"M  SPL" | b"M  SNC")
                if allow_sgroups =>
            {
                match tag {
                    b"M  STY" => sgroup_type_entries()
                        .parse(remaining)
                        .map(|(i, o)| (i, PropertyEntries::SGroupTypeEntries(o))),
                    b"M  SST" => sgroup_subtype_entries()
                        .parse(remaining)
                        .map(|(i, o)| (i, PropertyEntries::SGroupSubtypeEntries(o))),
                    b"M  SLB" => sgroup_label_entries()
                        .parse(remaining)
                        .map(|(i, o)| (i, PropertyEntries::SGroupLabelEntries(o))),
                    b"M  SCN" => sgroup_connectivity_entries()
                        .parse(remaining)
                        .map(|(i, o)| (i, PropertyEntries::SGroupConnectivityEntries(o))),
                    b"M  SAL" => sgroup_atom_list_entry()
                        .parse(remaining)
                        .map(|(i, o)| (i, PropertyEntries::SGroupAtomListEntry(o))),
                    b"M  SBL" => sgroup_bond_list_entry()
                        .parse(remaining)
                        .map(|(i, o)| (i, PropertyEntries::SGroupBondListEntry(o))),
                    b"M  SMT" => sgroup_subscript_entry()
                        .parse(remaining)
                        .map(|(i, o)| (i, PropertyEntries::SGroupSubscriptEntry(o))),
                    b"M  SDS" => sgroup_expansion_entries()
                        .parse(remaining)
                        .map(|(i, o)| (i, PropertyEntries::SGroupExpansionEntries(o))),
                    b"M  SPA" => sgroup_parent_atom_entries()
                        .parse(remaining)
                        .map(|(i, o)| (i, PropertyEntries::SGroupParentAtomEntry(o))),
                    b"M  CRS" => sgroup_correspondence_entry()
                        .parse(remaining)
                        .map(|(i, o)| (i, PropertyEntries::SGroupCorrespondenceEntry(o))),
                    b"M  SDI" => sgroup_display_info_entry()
                        .parse(remaining)
                        .map(|(i, o)| (i, PropertyEntries::SGroupDisplayInfoEntry(o))),
                    b"M  SBV" => sgroup_connecting_bond_entry()
                        .parse(remaining)
                        .map(|(i, o)| (i, PropertyEntries::SGroupConnectingBondEntry(o))),
                    b"M  SDT" => sgroup_data_description_entry()
                        .parse(remaining)
                        .map(|(i, o)| (i, PropertyEntries::SGroupDataDescriptionEntry(o))),
                    b"M  SDD" => sgroup_data_display_entry(skip_unused_fields)
                        .parse(remaining)
                        .map(|(i, o)| (i, PropertyEntries::SGroupDataDisplayEntry(o))),
                    b"M  SCD" => sgroup_data_continuation_entry()
                        .parse(remaining)
                        .map(|(i, o)| (i, PropertyEntries::SGroupDataEntry(o))),
                    b"M  SED" => sgroup_data_end_entry()
                        .parse(remaining)
                        .map(|(i, o)| (i, PropertyEntries::SGroupDataEntry(o))),
                    b"M  SPL" => sgroup_hierarchy_entries()
                        .parse(remaining)
                        .map(|(i, o)| (i, PropertyEntries::SGroupHierarchyEntries(o))),
                    b"M  SNC" => sgroup_component_entries()
                        .parse(remaining)
                        .map(|(i, o)| (i, PropertyEntries::SGroupComponentEntries(o))),
                    _ => unreachable!(),
                }
            }
            // Clark extensions
            tag @ (b"M  ZBO" | b"M  ZCH" | b"M  HYD") if allow_clark_extensions => match tag {
                b"M  ZBO" => bond_order_override_entries()
                    .parse(remaining)
                    .map(|(i, o)| (i, PropertyEntries::BondOrderOverrideEntries(o))),
                b"M  ZCH" => atom_charge_overrides_entries()
                    .parse(remaining)
                    .map(|(i, o)| (i, PropertyEntries::AtomChargeOverrideEntries(o))),
                b"M  HYD" => atom_hydrogen_count_entries()
                    .parse(remaining)
                    .map(|(i, o)| (i, PropertyEntries::AtomHydrogenCountEntries(o))),
                _ => unreachable!(),
            },
            _ => Err(Err::Error(NomError::new(input, NomErrorKind::Tag))),
        }?;
        Ok((remaining, result))
    }
}

/// Create atom alias property from two lines.
/// A  aaa
/// x...
/// aaa: atom index, x...: alias text
pub(super) fn parse_atom_alias_input<'a>(
    first_line: &'a [u8],
    second_line: &[u8],
) -> Result<PropertyEntries, NomError<&'a [u8]>> {
    let (_, atom_index) = preceded(tag("A  "), fixed_width_int_minus1::<usize>(3))
        .parse(first_line)
        .map_err(|_: Err<NomError<_>>| NomError::new(first_line, NomErrorKind::Digit))?;

    Ok(PropertyEntries::AtomAliasEntry(AtomAliasEntry {
        atom_index,
        alias: second_line.to_str_lossy().into_owned(),
    }))
}

/// Create legacy group abbreviation property from two lines.
/// G  aaappp
/// x...
/// aaa: atom index1, ppp: atom index2, x...: abbreviation label
/// The atoms on the side of atom index1 are abbreviated
pub(super) fn parse_legacy_group_abbreviation_input<'a>(
    first_line: &'a [u8],
    second_line: &[u8],
) -> Result<PropertyEntries, NomError<&'a [u8]>> {
    let (_, (atom_index1, atom_index2)) = preceded(
        tag("G  "),
        (
            fixed_width_int_minus1::<usize>(3),
            fixed_width_int_minus1::<usize>(3),
        ),
    )
    .parse(first_line)
    .map_err(|_: Err<NomError<_>>| NomError::new(first_line, NomErrorKind::Digit))?;

    Ok(PropertyEntries::LegacyGroupAbbreviationEntry(
        LegacyGroupAbbreviationEntry {
            atom_index1,
            atom_index2,
            label: second_line.to_str_lossy().into_owned(),
        },
    ))
}

/// Parse atom value entry.
/// V  aaa v..
/// aaa: atom index, v..: value string (can contain spaces)
fn atom_value_entry<'inp>(
) -> impl Parser<&'inp [u8], Output = AtomValueEntry, Error = NomError<&'inp [u8]>> {
    map_res(
        (fixed_width_int_minus1::<usize>(3), preceded(tag(" "), rest)),
        move |(atom_index, value)| -> Result<AtomValueEntry, NomError<&[u8]>> {
            Ok(AtomValueEntry {
                atom_index,
                value: value.to_str_lossy().into_owned(),
            })
        },
    )
}

/// Parse charge property entries.
/// nn8 aaa vvv ...
/// vvv: -15..= 15.
fn charge_entries<'inp>(
) -> impl Parser<&'inp [u8], Output = Vec<ChargeEntry>, Error = NomError<&'inp [u8]>> {
    length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                preceded(tag(" "), fixed_width_int_in_range::<i8, _>(3, -15..=15)),
            ),
            |(atom_index, charge)| ChargeEntry { atom_index, charge },
        ),
    )
}

/// Parse radical property entries.
/// nn8 aaa vvv ...
/// vvv: 0..= 3: 0 = no radical, 1 = singlet (:), 2 = doublet (. or ^), 3 = triplet (^^).
fn radical_entries<'inp>(
) -> impl Parser<&'inp [u8], Output = Vec<RadicalEntry>, Error = NomError<&'inp [u8]>> {
    length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                preceded(tag(" "), fixed_width_int_in_range::<u8, _>(3, 0..=3)),
            ),
            |(atom_index, radical_type)| RadicalEntry {
                atom_index,
                radical_type,
            },
        ),
    )
}

/// Parse isotope property entries.
/// nn8 aaa vvv ...
/// vvv: isotope mass number (not difference)
/// Difference between the isotope mass number and reference isotope mass number
/// should be in the range -18..=12.
fn isotope_entries<'inp>(
) -> impl Parser<&'inp [u8], Output = Vec<IsotopeEntry>, Error = NomError<&'inp [u8]>> {
    length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                preceded(tag(" "), fixed_width_int::<u32>(3)),
            ),
            |(atom_index, mass)| IsotopeEntry { atom_index, mass },
        ),
    )
}

/// Parse ring bond count property entries.
/// M  RBCnn8 aaa vvv ...
/// vvv: Ring bond count (-2 = as drawn (r*), -1 = no ring bonds (r0), 0 = off, 2 = r2, 3 = r3, 4 = r4+)
fn ring_bond_count_entries<'inp>(
) -> impl Parser<&'inp [u8], Output = Vec<RingBondCountEntry>, Error = NomError<&'inp [u8]>> {
    move |input: &'inp [u8]| {
        length_count(
            fixed_width_int_in_range::<u8, _>(3, 1..=8),
            map(
                (
                    preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                    preceded(tag(" "), fixed_width_int_in_range::<i8, _>(3, -2..=4)),
                ),
                |(atom_index, ring_bond_count)| RingBondCountEntry {
                    atom_index,
                    ring_bond_count,
                },
            ),
        )
        .parse(input)
    }
}

/// Parse substitution count property entries.
/// M  SUBnn8 aaa vvv ...
/// vvv: Substitution count (-2 = as drawn (s*), -1 = no substitution (s0), 0 = off, 1-5 = s1-s5,
/// 6 = s6+)
fn substitution_count_entries<'inp>(
) -> impl Parser<&'inp [u8], Output = Vec<SubstitutionCountEntry>, Error = NomError<&'inp [u8]>> {
    length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                preceded(tag(" "), fixed_width_int_in_range::<i8, _>(3, -2..=15)),
            ),
            |(atom_index, substitution_count)| SubstitutionCountEntry {
                atom_index,
                substitution_count,
            },
        ),
    )
}

/// Parse unsaturated atom property entries.
/// M  UNSnn8 aaa vvv ...
/// vvv: Unsaturated flag (0 = off, 1 = on)
fn unsaturated_atom_entries<'inp>(
) -> impl Parser<&'inp [u8], Output = Vec<UnsaturatedAtomEntry>, Error = NomError<&'inp [u8]>> {
    length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                preceded(tag(" "), fixed_width_int_in_range::<u8, _>(3, 0..=1)),
            ),
            |(atom_index, unsaturated)| UnsaturatedAtomEntry {
                atom_index,
                unsaturated,
            },
        ),
    )
}

/// Parse link atom property entries.
/// M  LINnn4 aaa vvv bbb ccc
/// nn4: Count (max 4), aaa: Atom index, vvv: Upper repeat count (>= 2, lower repeat count is 1),
/// bbb/ccc: Substituent indices (can be 0)
fn link_atom_entries<'inp>(
) -> impl Parser<&'inp [u8], Output = Vec<LinkAtomEntry>, Error = NomError<&'inp [u8]>> {
    length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=4),
        map(
            (
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                preceded(tag(" "), fixed_width_int_in_range::<u8, _>(3, 2..=255)),
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                opt(preceded(tag(" "), fixed_width_int_minus1::<usize>(3))),
            ),
            |(atom_index, repeat_count, subs_index1, subs_index2)| LinkAtomEntry {
                atom_index,
                repeat_count,
                subs_index1,
                subs_index2,
            },
        ),
    )
}

/// Parse atom list property entry.
/// M  ALS aaannn e 11112222333344445555...
/// aaa: Atom number, nnn: Number of entries (16 max), e: Exclusion (T/F), 1111: 4-char symbols
fn atom_list_entry<'inp>(
) -> impl Parser<&'inp [u8], Output = AtomListEntry, Error = NomError<&'inp [u8]>> {
    move |input: &'inp [u8]| {
        let (remaining, atom_index) =
            preceded(tag(" "), fixed_width_int_minus1::<usize>(3)).parse(input)?;

        let (remaining, count) =
            fixed_width_int_in_range::<usize, _>(3, 1..=16).parse(remaining)?;

        let (remaining, exclusion_byte) =
            delimited(tag(" "), take(1usize), tag(" ")).parse(remaining)?;
        let exclusion = match exclusion_byte {
            b"T" => true,
            b"F" | b" " => false,
            _ => return Err(Err::Error(NomError::new(input, NomErrorKind::Tag))),
        };

        let (remaining, elements) =
            nom_count(fixed_width_element_partial(4), count).parse(remaining)?;
        let elements = elements
            .iter()
            .copied()
            .collect::<Option<Vec<Element>>>()
            .ok_or_else(|| Err::Error(NomError::new(remaining, NomErrorKind::Verify)))?;

        Ok((
            remaining,
            AtomListEntry {
                atom_index,
                exclusion,
                elements,
            },
        ))
    }
}

/// Parse attachment point property entries.
/// Atom aaa is typically on ordinary atom, does not have to be a RGroup (opposite of AAL)
/// Attachment point does not appear in the atom list
/// M  APOnn2 aaa vvv ...
/// nn2: count (max 2), aaa: atom index, vvv: attachment type (0-3)
/// 0 = no attachment, 1 = first attachment point, 2 = second attachment point, 3 = both attachment points
fn attachment_point_entries<'inp>(
) -> impl Parser<&'inp [u8], Output = Vec<AttachmentPointEntry>, Error = NomError<&'inp [u8]>> {
    length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=2),
        map(
            (
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                preceded(tag(" "), fixed_width_int_in_range::<u8, _>(3, 0..=3)),
            ),
            |(atom_index, attachment_type)| AttachmentPointEntry {
                atom_index,
                attachment_type,
            },
        ),
    )
}

/// Parse atom attachment order entry.
/// M  AAL aaann2 111 v1v 222 v2v ...
/// Atom aaa refers to an RGroup, atoms 111, 222 are ordinary atoms (opposite of APO)
/// aaa: atom index, n2: pair count (max 2), 111/222: neighbor indices, v1v/v2v: attachment orders
fn atom_attachment_order_entry<'inp>(
) -> impl Parser<&'inp [u8], Output = AtomAttachmentOrderEntry, Error = NomError<&'inp [u8]>> {
    map(
        (
            preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
            length_count(
                fixed_width_int_in_range::<u8, _>(3, 1..=2),
                (
                    preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                    preceded(tag(" "), fixed_width_int_in_range::<u8, _>(3, 1..=2)),
                ),
            ),
        ),
        |(atom_index, attachments)| AtomAttachmentOrderEntry {
            atom_index,
            attachments,
        },
    )
}

/// Parse RGroup label property entries.
/// M  RGPnn8 aaa rrr ...
/// aaa: atom index, rrr: RGroup label (0 = no label, 1-32 in CTab spec, 1-999 in RDKit)
fn rgroup_label_entries<'inp>(
) -> impl Parser<&'inp [u8], Output = Vec<RGroupLabelEntry>, Error = NomError<&'inp [u8]>> {
    length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                preceded(tag(" "), fixed_width_int::<u32>(3)),
            ),
            |(atom_index, label)| RGroupLabelEntry { atom_index, label },
        ),
    )
}

// Parse RGroup logic entry.
/// M  LOGnn1 rrr iii hhh ooo...
/// nn1: count (max 1), rrr: RGroup label (1-32 in CTab spec, 1-999 in RDKit)
/// iii: Number of dependent Rgroup (IF rrr THEN iii)
/// hhh: REstH property of rrr (0 (default)=off, 1=on): RGroup or H atom
/// ooo...: Range of RGroup occurrence required: n=exactly n, n-m=from n through m (inclusive),
/// >n=greater n, <n=fewer than n, blank (default): > 0.
/// > Any non-contradictory combination of the preceding values is allowed, separated by commas.
fn rgroup_logic_entry<'inp>(
) -> impl Parser<&'inp [u8], Output = RGroupLogicEntry, Error = NomError<&'inp [u8]>> {
    move |input: &'inp [u8]| {
        let count = fixed_width_int_in_range::<u8, _>(3, 1..=1);
        let label = preceded(tag(" "), fixed_width_int_in_range::<u32, _>(3, 1..=999));
        let dependent_label = preceded(
            tag(" "),
            map(fixed_width_int::<u32>(3), |label| {
                if label == 0 {
                    None
                } else {
                    Some(label)
                }
            }),
        );
        let (i, (count, label, dependent_label)) = (count, label, dependent_label).parse(input)?;
        debug_assert!(count == 1, "count should be 1");

        let (i, rgroup_or_h) = map(
            cond(
                i.len() >= 4,
                preceded(tag(" "), fixed_width_int_in_range::<u8, _>(3, 0..=1)),
            ),
            |r| r.is_some_and(|r| r != 0),
        )
        .parse(i)?;

        let (i, occurrence) = map(cond(!i.is_empty(), rgroup_occurrences()), |o| {
            o.unwrap_or(vec![RGroupOccurrence::GreaterThan(0)])
        })
        .parse(i)?;

        Ok((
            i,
            RGroupLogicEntry {
                label,
                dependent_label,
                rgroup_or_h,
                occurrence,
            },
        ))
    }
}

/// Parse SGroup type entries.
/// M  STYnn8 sss ttt ...
/// sss: SGroup index, ttt: SGroup type (3-character string)
fn sgroup_type_entries<'inp>(
) -> impl Parser<&'inp [u8], Output = Vec<SGroupTypeEntry>, Error = NomError<&'inp [u8]>> {
    length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                preceded(tag(" "), sgroup_type()),
            ),
            |(sgroup_index, sgroup_type)| SGroupTypeEntry {
                sgroup_index,
                sgroup_type,
            },
        ),
    )
}

/// Parse SGroup subtype entries.
/// M  SSTnn8 sss ttt ...
/// sss: SGroup index, ttt: SGroup subtype (3-character string)
fn sgroup_subtype_entries<'inp>(
) -> impl Parser<&'inp [u8], Output = Vec<SGroupSubtypeEntry>, Error = NomError<&'inp [u8]>> {
    all_consuming(length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                preceded(tag(" "), sgroup_subtype()),
            ),
            |(sgroup_index, sgroup_subtype)| SGroupSubtypeEntry {
                sgroup_index,
                sgroup_subtype,
            },
        ),
    ))
}

/// Parse SGroup label entries.
/// M  SLBnn8 sss vvv ...
/// sss: SGroup index, vvv: integer label is from 1-512
fn sgroup_label_entries<'inp>(
) -> impl Parser<&'inp [u8], Output = Vec<SGroupLabelEntry>, Error = NomError<&'inp [u8]>> {
    length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                preceded(tag(" "), fixed_width_int_in_range::<u32, _>(3, 1..=512)),
            ),
            |(sgroup_index, label)| SGroupLabelEntry {
                sgroup_index,
                label,
            },
        ),
    )
}

/// Parse SGroup connectivity entries.
/// M  SCNnn8 sss ttt ...
/// sss: SGroup index, ttt: SGroup connectivity (2-character string), left-justified
fn sgroup_connectivity_entries<'inp>(
) -> impl Parser<&'inp [u8], Output = Vec<SGroupConnectivityEntry>, Error = NomError<&'inp [u8]>> {
    length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                preceded(tag(" "), sgroup_connectivity()),
            ),
            |(sgroup_index, connectivity)| SGroupConnectivityEntry {
                sgroup_index,
                connectivity,
            },
        ),
    )
}

/// Parse SGroup expansion entries.
/// M SDS EXPn15 sss ...
/// sss: SGroup index, n15: count (max 15)
fn sgroup_expansion_entries<'inp>(
) -> impl Parser<&'inp [u8], Output = Vec<SGroupExpansionEntry>, Error = NomError<&'inp [u8]>> {
    preceded(
        tag(" EXP"),
        length_count(
            fixed_width_int_in_range::<u8, _>(3, 1..=15),
            map(
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                |sgroup_index| SGroupExpansionEntry { sgroup_index },
            ),
        ),
    )
}

/// Parse SGroup atom list entry.
/// M  SAL sssn15 aaa ...
/// sss: SGroup index, n15: count (max 15), aaa: atom indices (each 4 chars: " aaa")
fn sgroup_atom_list_entry<'inp>(
) -> impl Parser<&'inp [u8], Output = SGroupAtomListEntry, Error = NomError<&'inp [u8]>> {
    map(
        (
            preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
            length_count(
                fixed_width_int_in_range::<u8, _>(3, 1..=15),
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
            ),
        ),
        |(sgroup_index, atom_indices)| SGroupAtomListEntry {
            sgroup_index,
            atom_indices,
        },
    )
}

/// Parse SGroup bond list entry.
/// M  SBL sssn15 bbb ...
/// sss: SGroup index (3 chars), n: count (3 chars), bbb: bond indices (each 4 chars: " bbb")
fn sgroup_bond_list_entry<'inp>(
) -> impl Parser<&'inp [u8], Output = SGroupBondListEntry, Error = NomError<&'inp [u8]>> {
    map(
        (
            preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
            length_count(
                fixed_width_int_in_range::<u8, _>(3, 1..=15),
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
            ),
        ),
        |(sgroup_index, bond_indices)| SGroupBondListEntry {
            sgroup_index,
            bond_indices,
        },
    )
}

/// Parse SGroup parent atom entries.
/// M  SPA sssn15 aaa ...
/// sss: SGroup index, n15: count (max 15), aaa: atom indices (each 4 chars: " aaa")
fn sgroup_parent_atom_entries<'inp>(
) -> impl Parser<&'inp [u8], Output = SGroupParentAtomEntry, Error = NomError<&'inp [u8]>> {
    map(
        (
            preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
            length_count(
                fixed_width_int_in_range::<u8, _>(3, 1..=15),
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
            ),
        ),
        |(sgroup_index, atom_indices)| SGroupParentAtomEntry {
            sgroup_index,
            atom_indices,
        },
    )
}

/// Parse SGroup subscript entry.
/// M  SMT sss m...
/// sss: SGroup index, m: subscript text
/// For multiple groups, m... is the text representation of the multiple group multiplier.For superatoms,
/// m... is the text of the superatom label.)
fn sgroup_subscript_entry<'inp>(
) -> impl Parser<&'inp [u8], Output = SGroupSubscriptEntry, Error = NomError<&'inp [u8]>> {
    map(
        (
            preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
            preceded(
                tag(" "),
                map_res(rest, move |s: &[u8]| {
                    let multiplier = sgroup_multiplier()
                        .parse(s)
                        .ok()
                        .map(|(_, multiplier)| multiplier);
                    let subscript = sgroup_subscript()
                        .parse(s)
                        .ok()
                        .map(|(_, subscript)| subscript);
                    if multiplier.is_none() && subscript.is_none() {
                        Err(NomError::new(s, NomErrorKind::MapRes))
                    } else {
                        Ok((multiplier, subscript))
                    }
                }),
            ),
        ),
        |(sgroup_index, (multiplier, subscript))| SGroupSubscriptEntry {
            sgroup_index,
            multiplier,
            subscript,
        },
    )
}

/// Parse SGroup correspondence entry.
/// M  CRS sssnn6 bb1 bb2 bb3
/// sss: SGroup index, nn6: count (max 6), bb1-bb3: bond indices
/// bb1, bb2: Crossing bonds that share a common bracket
/// bb3: Crossing bond in repeating unit that connects to bond bb1
fn sgroup_correspondence_entry<'inp>(
) -> impl Parser<&'inp [u8], Output = SGroupCorrespondenceEntry, Error = NomError<&'inp [u8]>> {
    map(
        (
            preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
            length_count(
                fixed_width_int_in_range::<usize, _>(3, 1..=6),
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
            ),
        ),
        |(sgroup_index, bond_indices)| SGroupCorrespondenceEntry {
            sgroup_index,
            bond_indices,
        },
    )
}

/// Parse SGroup display info entry.
/// M  SDI sssnn4 x1 y1 x2 y2
/// sss: SGroup index, nn4: count (max 4), x1, y1: opening bracket, x2, y2: closing bracket (f10.4)
fn sgroup_display_info_entry<'inp>(
) -> impl Parser<&'inp [u8], Output = SGroupDisplayInfoEntry, Error = NomError<&'inp [u8]>> {
    map(
        (
            preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
            length_count(
                fixed_width_int_in_range::<usize, _>(3, 1..=4),
                fixed_width_float::<f64>(10, 4),
            ),
        ),
        |(sgroup_index, bracket_coords)| SGroupDisplayInfoEntry {
            sgroup_index,
            bracket_coords,
        },
    )
}

/// Parse SGroup connecting bond entry.
/// M  SBV sss bb1 x1 y1
/// sss: SGroup index, bb1: bond index, x1, y1: bond vector (f10.4)
fn sgroup_connecting_bond_entry<'inp>(
) -> impl Parser<&'inp [u8], Output = SGroupConnectingBondEntry, Error = NomError<&'inp [u8]>> {
    map(
        (
            preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
            preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
            fixed_width_float::<f64>(10, 4),
            fixed_width_float::<f64>(10, 4),
        ),
        |(sgroup_index, bond_index, x1, y1)| SGroupConnectingBondEntry {
            sgroup_index,
            bond_index,
            bond_vector: (x1, y1),
        },
    )
}

/// Parse SGroup data field description entry.
/// M  SDT sss fff..fffgghh...hhhiijjj..,
/// sss: SGroup index, fff..fff: field name (30 chars), gg: field type (2 chars, F=formatted, N=numeric, T=text)
/// hh: field units or format (20 chars), ii: Query identifier (MQ: MACCS-II, IQ: ISIS, PQ: program code)
/// jjj: data query operator
fn sgroup_data_description_entry<'inp>(
) -> impl Parser<&'inp [u8], Output = SGroupDataDescriptionEntry, Error = NomError<&'inp [u8]>> {
    map(
        (
            preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
            preceded(tag(" "), fixed_width_str_partial(30)),
            sgroup_data_type(),
            fixed_width_str_partial(20),
            fixed_width_str_partial(2),
            fixed_width_str_partial(1),
        ),
        |(
            sgroup_index,
            field_name,
            field_type,
            field_units,
            query_identifier,
            data_query_operator,
        )| {
            SGroupDataDescriptionEntry {
                sgroup_index,
                field_name: field_name.unwrap_or_default(),
                field_type,
                field_units,
                query_identifier,
                data_query_operator,
            }
        },
    )
}

/// Parse SGroup data display entry.
/// M  SDD sss xxxxx.xxxxyyyyy.yyyy eeefgh i jjjkkk ll m  noo
/// sss: SGroup index, x, y: coordinates (f10.4), eee: skipped, f: data display (A=attached, D=detached)
/// g: absolute, relative placement (A=absolute, R=relative), h: display units (" "=none, U=display units)
/// i: skipped, jjj: number of characters to display (1-999 or ALL), kkk: number of lines to display (unused, always 1)
/// ll: skipped, m: tag character (if non-blank), n: Data display DASP position (1-9), oo: skipped
fn sgroup_data_display_entry<'inp>(
    skip_unused_fields: bool,
) -> impl Parser<&'inp [u8], Output = SGroupDataDisplayEntry, Error = NomError<&'inp [u8]>> {
    map(
        (
            preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
            preceded(tag(" "), fixed_width_float::<f64>(10, 4)),
            terminated(fixed_width_float::<f64>(10, 4), tag(" ")),
            preceded(
                fixed_width_unused(3, skip_unused_fields),
                sgroup_data_display_type(),
            ),
            sgroup_data_display_placement(),
            sgroup_data_display_units(),
            preceded(
                fixed_width_unused(3, skip_unused_fields),
                sgroup_data_display_chars(),
            ),
            preceded(
                fixed_width_unused(7, skip_unused_fields),
                map_parser(take(1usize), opt(nom_u8)),
            ),
            preceded(tag("  "), fixed_width_int::<u8>(1)),
        ),
        |(
            sgroup_index,
            x,
            y,
            display_type,
            display_placement,
            display_units,
            display_chars,
            display_tag,
            display_position,
        )| SGroupDataDisplayEntry {
            sgroup_index,
            coords: (x, y),
            display_type,
            display_placement,
            display_units,
            display_chars,
            display_tag,
            display_position,
        },
    )
}

/// Parse SGroup data continuation entry.
/// M  SCD sss d...
/// sss: SGroup index, d...: data content (max 69 chars)
fn sgroup_data_continuation_entry<'inp>(
) -> impl Parser<&'inp [u8], Output = SGroupDataEntry, Error = NomError<&'inp [u8]>> {
    map_res(
        (
            preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
            opt(preceded(tag(" "), rest)),
        ),
        move |(sgroup_index, data_content)| -> Result<SGroupDataEntry, NomError<&[u8]>> {
            let data_content = data_content.unwrap_or_default();
            Ok(SGroupDataEntry::Continuation {
                sgroup_index,
                data_content: data_content[..data_content.len().min(69)]
                    .to_str_lossy()
                    .into_owned(),
            })
        },
    )
}

/// Parse SGroup data end entry.
/// M  SED sss d... OR
/// M  SED
/// sss: SGroup index, d...: data content (max 69 chars)
/// The data content is the same as the data content in the SCD entry.
/// The second form is used to indicate the end of the data content, in which case it should be empty.
fn sgroup_data_end_entry<'inp>(
) -> impl Parser<&'inp [u8], Output = SGroupDataEntry, Error = NomError<&'inp [u8]>> {
    alt((
        map_res(
            (
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                preceded(tag(" "), rest),
            ),
            move |(sgroup_index, data_content)| -> Result<SGroupDataEntry, NomError<&[u8]>> {
                Ok(SGroupDataEntry::EndWithData {
                    sgroup_index,
                    data_content: data_content[..data_content.len().min(69)]
                        .to_str_lossy()
                        .into_owned(),
                })
            },
        ),
        map(
            preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
            |sgroup_index| SGroupDataEntry::EndBlank { sgroup_index },
        ),
    ))
}

/// Parse SGroup hierarchy entries.
/// M  SPLnn8 sss ppp ...
/// sss: SGroup index, ppp: Parent SGroup index
fn sgroup_hierarchy_entries<'inp>(
) -> impl Parser<&'inp [u8], Output = Vec<SGroupHierarchyEntry>, Error = NomError<&'inp [u8]>> {
    all_consuming(length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
            ),
            |(sgroup_index, parent_sgroup_index)| SGroupHierarchyEntry {
                sgroup_index,
                parent_sgroup_index,
            },
        ),
    ))
}

/// Parse SGroup component number entries.
/// M  SNCnn8 sss ccc ...
/// sss: SGroup index, ccc: Component number
fn sgroup_component_entries<'inp>(
) -> impl Parser<&'inp [u8], Output = Vec<SGroupComponentEntry>, Error = NomError<&'inp [u8]>> {
    length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                preceded(tag(" "), fixed_width_int::<u32>(3)),
            ),
            |(sgroup_index, component_number)| SGroupComponentEntry {
                sgroup_index,
                component_number,
            },
        ),
    )
}

/// Parse zero bond order entries.
/// M  ZBOnn8 bbb vvv ...
/// bbb: bond index, vvv: bond vector override (>= 0, limited to 6 here)
fn bond_order_override_entries<'inp>(
) -> impl Parser<&'inp [u8], Output = Vec<BondOrderOverrideEntry>, Error = NomError<&'inp [u8]>> {
    length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map_opt(
            (
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                preceded(tag(" "), fixed_width_int_in_range::<u8, _>(3, 0..=6)),
            ),
            |(bond_index, bond_order)| {
                BondOrder::from_value(bond_order).map(|bo| BondOrderOverrideEntry {
                    bond_index,
                    bond_order: bo,
                })
            },
        ),
    )
}

/// Parse zero atom charge entries.
/// M  ZCHnn8 aaa ccc ...
/// aaa: atom index, ccc: atom charge override (-8..=8)
fn atom_charge_overrides_entries<'inp>(
) -> impl Parser<&'inp [u8], Output = Vec<AtomChargeOverrideEntry>, Error = NomError<&'inp [u8]>> {
    length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                preceded(tag(" "), fixed_width_int_in_range::<i8, _>(3, -8..=8)),
            ),
            |(atom_index, charge)| AtomChargeOverrideEntry { atom_index, charge },
        ),
    )
}

/// Parse atom explicit hydrogen count entries.
/// M  HCTnn8 aaa hhh ...
/// aaa: atom index, hhh: atom explicit hydrogen count (>= 0, limited to 8 here, -1 = no override)
fn atom_hydrogen_count_entries<'inp>(
) -> impl Parser<&'inp [u8], Output = Vec<AtomHydrogenCountEntry>, Error = NomError<&'inp [u8]>> {
    length_count(
        fixed_width_int_in_range::<u8, _>(3, 0..=8),
        map(
            (
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                preceded(tag(" "), fixed_width_int_in_range::<i8, _>(3, -1..=8)),
            ),
            |(atom_index, hydrogen_count)| AtomHydrogenCountEntry {
                atom_index,
                hydrogen_count: if hydrogen_count == -1 {
                    None
                } else {
                    Some(hydrogen_count as u8)
                },
            },
        ),
    )
}

#[cfg(test)]
mod tests;
