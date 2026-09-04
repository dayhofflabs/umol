//! Parsers for CTab property lines.

use bstr::ByteSlice;
use umol_chem::element::Element;
use winnow::combinator::{cond, delimited, opt, preceded, terminated};
use winnow::error::ErrMode;
use winnow::stream::Location;
use winnow::token::{rest, take};
use winnow::Parser;

use super::sgroup::{sgroup_connectivity, sgroup_subtype, sgroup_type};
use super::utils::{
    finish_line, fixed_width_element_partial, fixed_width_float_f10_4, fixed_width_int,
    fixed_width_int_in_range, fixed_width_int_minus1, fixed_width_partial, fixed_width_str_partial,
    fixed_width_unused, input_error_column, next_line, rgroup_occurrences, Input, InputError,
    IntParser, PResult,
};
use crate::ctfile::config::CtabParseFlags;
use crate::ctfile::error::ParseError;
use crate::ctfile::parser::sgroup::{
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
    pub atom_index: u32,
    pub alias: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct LegacyGroupAbbreviationEntry {
    pub atom_index1: u32, // Atoms on this side are abbreviated
    pub atom_index2: u32, // Attachment point to main structure
    pub label: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtomValueEntry {
    pub atom_index: u32,
    pub value: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ChargeEntry {
    pub atom_index: u32,
    pub charge: i8,
}

#[derive(Debug, Clone, PartialEq)]
pub struct RadicalEntry {
    pub atom_index: u32,
    pub radical_type: u8,
}

#[derive(Debug, Clone, PartialEq)]
pub struct IsotopeEntry {
    pub atom_index: u32,
    pub mass: u32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct RingBondCountEntry {
    pub atom_index: u32,
    pub ring_bond_count: i8, // -2=r*, -1=r0, 0=off, 2=r2, 3=r3, 4+=r4
}

#[derive(Debug, Clone, PartialEq)]
pub struct SubstitutionCountEntry {
    pub atom_index: u32,
    pub substitution_count: i8, // -2=s*, -1=s0, 0=off, 1-5=s1-s5, 6+=s6
}

#[derive(Debug, Clone, PartialEq)]
pub struct UnsaturatedAtomEntry {
    pub atom_index: u32,
    pub unsaturated: u8, // 0=off, 1=on
}

#[derive(Debug, Clone, PartialEq)]
pub struct LinkAtomEntry {
    pub atom_index: u32,
    pub repeat_count: u8,         // vvv >= 2
    pub subs_index1: u32,         // bbb
    pub subs_index2: Option<u32>, // ccc (optional)
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtomListEntry {
    pub atom_index: u32,
    pub exclusion: bool,        // T = NOT list, F = normal list
    pub elements: Vec<Element>, // Converted from 4-char symbols
}

#[derive(Debug, Clone, PartialEq)]
pub struct AttachmentPointEntry {
    pub atom_index: u32,
    pub attachment_type: u8, // 0=none, 1=first, 2=second, 3=both
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtomAttachmentOrderEntry {
    pub atom_index: u32,
    pub attachments: Vec<(u32, u8)>, // (neighbor index, attachment order)
}

#[derive(Debug, Clone, PartialEq)]
pub struct RGroupLabelEntry {
    pub atom_index: u32,
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
    pub sgroup_index: u32,
    pub sgroup_type: SGroupType,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupSubtypeEntry {
    pub sgroup_index: u32,
    pub sgroup_subtype: SGroupSubtype,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupLabelEntry {
    pub sgroup_index: u32,
    pub label: u32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupConnectivityEntry {
    pub sgroup_index: u32,
    pub connectivity: SGroupConnectivity,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupExpansionEntry {
    pub sgroup_index: u32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupAtomListEntry {
    pub sgroup_index: u32,
    pub atom_indices: Vec<u32>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupBondListEntry {
    pub sgroup_index: u32,
    pub bond_indices: Vec<u32>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupParentAtomEntry {
    pub sgroup_index: u32,
    pub atom_indices: Vec<u32>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupSubscriptEntry {
    pub sgroup_index: u32,
    pub multiplier: Option<SGroupMultiplier>,
    pub subscript: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupCorrespondenceEntry {
    pub sgroup_index: u32,
    pub bond_indices: Vec<u32>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupDisplayInfoEntry {
    pub sgroup_index: u32,
    pub bracket_coords: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupConnectingBondEntry {
    pub sgroup_index: u32,
    pub bond_index: u32,
    pub bond_vector: (f64, f64),
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupDataDescriptionEntry {
    pub sgroup_index: u32,
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
        sgroup_index: u32,
        data_content: String, // Should be exactly 69 chars
    },
    /// SED - End with data (<= 69 characters)
    EndWithData {
        sgroup_index: u32,
        data_content: String, // 0-69 chars
    },
    /// SED - End without data (blank, processes buffered data)
    EndBlank { sgroup_index: u32 },
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupDataDisplayEntry {
    pub sgroup_index: u32,
    pub coords: (f64, f64),
    pub display_type: SGroupDataDisplayType,
    pub display_placement: SGroupDataDisplayPlacement,
    pub display_units: SGroupDataDisplayUnits,
    pub display_chars: SGroupDataDisplayChars,
    pub display_tag: Option<u8>, // 0 = no tag, 1-9 = tag
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupComponentEntry {
    pub sgroup_index: u32,
    pub component_number: u32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupHierarchyEntry {
    pub sgroup_index: u32,
    pub parent_sgroup_index: u32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BondOrderOverrideEntry {
    pub bond_index: u32,
    pub bond_order: BondOrder,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtomChargeOverrideEntry {
    pub atom_index: u32,
    pub charge: i8,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtomHydrogenCountEntry {
    pub atom_index: u32,
    pub hydrogen_count: Option<u8>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ChemSketchLabelEntry {
    pub atom_index: u32,
    pub label: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MarvinSmartsPatternEntry {
    pub atom_index: u32,
    pub smarts_pattern: String,
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
    ChemSketchLabelEntry(ChemSketchLabelEntry),
    MarvinSmartsPatternEntry(MarvinSmartsPatternEntry),
}

/// Parse properties block (basic properties only)
pub(super) fn properties_block<'inp>(
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], (Vec<PropertyEntries>, u32), ErrMode<ParseError>> + use<'inp> {
    let no_v2000_end_tags = flags.contains(CtabParseFlags::NO_V2000_END_TAGS);

    move |input: &mut &'inp [u8]| {
        let mut properties = Vec::new();
        let mut line_index = 0;
        let mut end_found = false;

        while !input.is_empty() {
            let mut line = next_line(input).expect("non-empty input contains a physical line");
            if line.as_ref().starts_with(b"M  END") {
                line_index += 1;
                end_found = true;
                break;
            }

            if line.as_ref().starts_with(b"A  ") {
                let alias = next_line(input).map_err(|_| {
                    ErrMode::Cut(ParseError::UnexpectedEof {
                        line: line_offset + line_index + 1,
                        block: "atom alias",
                    })
                })?;
                let property =
                    parse_atom_alias_input(line.as_ref(), alias.as_ref()).map_err(|error| {
                        ErrMode::Cut(ParseError::InvalidPropertyLine {
                            line: line_offset + line_index,
                            col: input_error_column(error, &line),
                        })
                    })?;
                properties.push(property);
                line_index += 2;
            } else {
                let property = property_input(flags)
                    .parse_next(&mut line)
                    .and_then(|value| {
                        finish_line(&mut line)?;
                        Ok(value)
                    })
                    .map_err(|error| {
                        ErrMode::Cut(ParseError::InvalidPropertyLine {
                            line: line_offset + line_index,
                            col: input_error_column(error, &line),
                        })
                    })?;
                properties.push(property);
                line_index += 1;
            }
        }

        if !end_found && !no_v2000_end_tags {
            return Err(ErrMode::Cut(ParseError::MissingMEndTag {
                line: line_offset + line_index,
            }));
        }

        Ok((properties, line_offset + line_index))
    }
}

/// Parse extended properties block.
pub(super) fn extended_properties_block<'inp>(
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], (Vec<PropertyEntries>, u32), ErrMode<ParseError>> + use<'inp> {
    let no_v2000_end_tags = flags.contains(CtabParseFlags::NO_V2000_END_TAGS);

    move |input: &mut &'inp [u8]| {
        let mut properties = Vec::new();
        let mut line_index = 0;
        let mut end_found = false;

        while !input.is_empty() {
            let mut line = next_line(input).expect("non-empty input contains a physical line");
            if line.as_ref().starts_with(b"M  END") {
                line_index += 1;
                end_found = true;
                break;
            }

            if line.as_ref().starts_with(b"A  ") {
                let alias = next_line(input).map_err(|_| {
                    ErrMode::Cut(ParseError::UnexpectedEof {
                        line: line_offset + line_index + 1,
                        block: "atom alias",
                    })
                })?;
                let property =
                    parse_atom_alias_input(line.as_ref(), alias.as_ref()).map_err(|error| {
                        ErrMode::Cut(ParseError::InvalidPropertyLine {
                            line: line_offset + line_index,
                            col: input_error_column(error, &line),
                        })
                    })?;
                properties.push(property);
                line_index += 2;
            } else if line.as_ref().starts_with(b"G  ") {
                let label = next_line(input).map_err(|_| {
                    ErrMode::Cut(ParseError::UnexpectedEof {
                        line: line_offset + line_index + 1,
                        block: "legacy group abbreviation",
                    })
                })?;
                let property = parse_legacy_group_abbreviation_input(line.as_ref(), label.as_ref())
                    .map_err(|error| {
                        ErrMode::Cut(ParseError::InvalidPropertyLine {
                            line: line_offset + line_index,
                            col: input_error_column(error, &line),
                        })
                    })?;
                properties.push(property);
                line_index += 2;
            } else {
                let property = extended_property_input(flags)
                    .parse_next(&mut line)
                    .and_then(|value| {
                        finish_line(&mut line)?;
                        Ok(value)
                    })
                    .map_err(|error| {
                        ErrMode::Cut(ParseError::InvalidPropertyLine {
                            line: line_offset + line_index,
                            col: input_error_column(error, &line),
                        })
                    })?;
                properties.push(property);
                line_index += 1;
            }
        }

        if !end_found && !no_v2000_end_tags {
            return Err(ErrMode::Cut(ParseError::MissingMEndTag {
                line: line_offset + line_index,
            }));
        }

        Ok((properties, line_offset + line_index))
    }
}

/// Parse property line (basic properties only)
/// Note: A  (atom alias) and M  END are handled at the block level, not here.
fn property_input<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<Input<'inp>, PropertyEntries, ErrMode<InputError>> + use<'inp> {
    let allow_clark_extensions = flags.contains(CtabParseFlags::CLARK_EXTENSIONS);
    let allow_editor_extensions = flags.contains(CtabParseFlags::EDITOR_EXTENSIONS);
    move |input: &mut Input<'inp>| {
        let column = input.current_token_start();
        let bytes: &[u8] = input.as_ref();
        if bytes.len() < 3 {
            return Err(ErrMode::Backtrack(InputError::at_column(column)));
        }

        debug_assert!(&bytes[0..3] != b"A  ", "A should be handled at block level");
        if &bytes[0..3] == b"V  " {
            let _: &[u8] = take(3usize).parse_next(input)?;
            return atom_value_entry
                .map(PropertyEntries::AtomValueEntry)
                .parse_next(input)
                .map_err(ErrMode::cut);
        }

        if bytes.len() < 6 {
            return Err(ErrMode::Backtrack(InputError::at_column(column)));
        }

        debug_assert!(
            &bytes[0..6] != b"M  END",
            "M  END should be handled at block level"
        );
        let tag = &bytes[0..6];
        let known = matches!(tag, b"M  CHG" | b"M  RAD" | b"M  ISO")
            || (allow_clark_extensions && matches!(tag, b"M  ZBO" | b"M  ZCH" | b"M  HYD"))
            || (allow_editor_extensions && tag == b"M  ZZC");
        if !known {
            return Err(ErrMode::Backtrack(InputError::at_column(column)));
        }
        let _: &[u8] = take(6usize).parse_next(input)?;

        let result = match tag {
            b"M  CHG" => charge_entries
                .map(PropertyEntries::ChargeEntries)
                .parse_next(input),
            b"M  RAD" => radical_entries
                .map(PropertyEntries::RadicalEntries)
                .parse_next(input),
            b"M  ISO" => isotope_entries
                .map(PropertyEntries::IsotopeEntries)
                .parse_next(input),
            b"M  ZBO" => bond_order_override_entries
                .map(PropertyEntries::BondOrderOverrideEntries)
                .parse_next(input),
            b"M  ZCH" => atom_charge_overrides_entries
                .map(PropertyEntries::AtomChargeOverrideEntries)
                .parse_next(input),
            b"M  HYD" => atom_hydrogen_count_entries
                .map(PropertyEntries::AtomHydrogenCountEntries)
                .parse_next(input),
            b"M  ZZC" => chemsketch_label_entry
                .map(PropertyEntries::ChemSketchLabelEntry)
                .parse_next(input),
            _ => unreachable!(),
        };
        result.map_err(ErrMode::cut)
    }
}

/// Parse extended property line
/// Note: A  (atom alias) and M  END are handled at the block level, not here.
fn extended_property_input<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<Input<'inp>, PropertyEntries, ErrMode<InputError>> + use<'inp> {
    let allow_queries = flags.contains(CtabParseFlags::WILDCARDS);
    let allow_rgroups = flags.contains(CtabParseFlags::RGROUPS);
    let allow_sgroups = flags.contains(CtabParseFlags::SGROUPS);
    let allow_clark_extensions = flags.contains(CtabParseFlags::CLARK_EXTENSIONS);
    let allow_editor_extensions = flags.contains(CtabParseFlags::EDITOR_EXTENSIONS);
    let skip_unused_fields = flags.contains(CtabParseFlags::SKIP_UNUSED_FIELDS);
    move |input: &mut Input<'inp>| {
        let column = input.current_token_start();
        let bytes: &[u8] = input.as_ref();
        if bytes.len() < 3 {
            return Err(ErrMode::Backtrack(InputError::at_column(column)));
        }

        debug_assert!(
            &bytes[0..3] != b"A  " && &bytes[0..3] != b"G  ",
            "A or G should be handled at block level"
        );

        if &bytes[0..3] == b"V  " {
            let _: &[u8] = take(3usize).parse_next(input)?;
            return atom_value_entry
                .map(PropertyEntries::AtomValueEntry)
                .parse_next(input)
                .map_err(ErrMode::cut);
        }

        if bytes.len() < 6 {
            return Err(ErrMode::Backtrack(InputError::at_column(column)));
        }

        debug_assert!(
            &bytes[0..6] != b"M  END",
            "M  END should be handled at block level"
        );
        let tag = &bytes[0..6];
        let known = matches!(tag, b"M  CHG" | b"M  RAD" | b"M  ISO")
            || (allow_queries
                && matches!(
                    tag,
                    b"M  RBC" | b"M  SUB" | b"M  UNS" | b"M  LIN" | b"M  ALS"
                ))
            || (allow_rgroups && matches!(tag, b"M  APO" | b"M  AAL" | b"M  RGP" | b"M  LOG"))
            || (allow_sgroups
                && matches!(
                    tag,
                    b"M  STY"
                        | b"M  SST"
                        | b"M  SLB"
                        | b"M  SCN"
                        | b"M  SAL"
                        | b"M  SBL"
                        | b"M  SMT"
                        | b"M  SDS"
                        | b"M  SPA"
                        | b"M  CRS"
                        | b"M  SDI"
                        | b"M  SBV"
                        | b"M  SDT"
                        | b"M  SDD"
                        | b"M  SCD"
                        | b"M  SED"
                        | b"M  SPL"
                        | b"M  SNC"
                ))
            || (allow_clark_extensions && matches!(tag, b"M  ZBO" | b"M  ZCH" | b"M  HYD"))
            || (allow_editor_extensions
                && (tag == b"M  ZZC" || (allow_queries && tag == b"M  MRV")));
        if !known {
            return Err(ErrMode::Backtrack(InputError::at_column(column)));
        }
        let _: &[u8] = take(6usize).parse_next(input)?;

        let result = match tag {
            b"M  CHG" => charge_entries
                .map(PropertyEntries::ChargeEntries)
                .parse_next(input),
            b"M  RAD" => radical_entries
                .map(PropertyEntries::RadicalEntries)
                .parse_next(input),
            b"M  ISO" => isotope_entries
                .map(PropertyEntries::IsotopeEntries)
                .parse_next(input),
            tag @ (b"M  RBC" | b"M  SUB" | b"M  UNS" | b"M  LIN" | b"M  ALS") if allow_queries => {
                match tag {
                    b"M  RBC" => ring_bond_count_entries
                        .map(PropertyEntries::RingBondCountEntries)
                        .parse_next(input),
                    b"M  SUB" => substitution_count_entries
                        .map(PropertyEntries::SubstitutionCountEntries)
                        .parse_next(input),
                    b"M  UNS" => unsaturated_atom_entries
                        .map(PropertyEntries::UnsaturatedAtomEntries)
                        .parse_next(input),
                    b"M  LIN" => link_atom_entries
                        .map(PropertyEntries::LinkAtomEntries)
                        .parse_next(input),
                    b"M  ALS" => atom_list_entry
                        .map(PropertyEntries::AtomListEntry)
                        .parse_next(input),
                    _ => unreachable!(),
                }
            }
            tag @ (b"M  APO" | b"M  AAL" | b"M  RGP" | b"M  LOG") if allow_rgroups => match tag {
                b"M  APO" => attachment_point_entries
                    .map(PropertyEntries::AttachmentPointEntries)
                    .parse_next(input),
                b"M  AAL" => atom_attachment_order_entry
                    .map(PropertyEntries::AtomAttachmentOrderEntry)
                    .parse_next(input),
                b"M  RGP" => rgroup_label_entries
                    .map(PropertyEntries::RGroupLabelEntries)
                    .parse_next(input),
                b"M  LOG" => rgroup_logic_entry
                    .map(PropertyEntries::RGroupLogicEntry)
                    .parse_next(input),
                _ => unreachable!(),
            },
            tag @ (b"M  STY" | b"M  SST" | b"M  SLB" | b"M  SCN" | b"M  SAL" | b"M  SBL"
            | b"M  SMT" | b"M  SDS" | b"M  SPA" | b"M  CRS" | b"M  SDI" | b"M  SBV"
            | b"M  SDT" | b"M  SDD" | b"M  SCD" | b"M  SED" | b"M  SPL" | b"M  SNC")
                if allow_sgroups =>
            {
                match tag {
                    b"M  STY" => sgroup_type_entries
                        .map(PropertyEntries::SGroupTypeEntries)
                        .parse_next(input),
                    b"M  SST" => sgroup_subtype_entries
                        .map(PropertyEntries::SGroupSubtypeEntries)
                        .parse_next(input),
                    b"M  SLB" => sgroup_label_entries
                        .map(PropertyEntries::SGroupLabelEntries)
                        .parse_next(input),
                    b"M  SCN" => sgroup_connectivity_entries
                        .map(PropertyEntries::SGroupConnectivityEntries)
                        .parse_next(input),
                    b"M  SAL" => sgroup_atom_list_entry
                        .map(PropertyEntries::SGroupAtomListEntry)
                        .parse_next(input),
                    b"M  SBL" => sgroup_bond_list_entry
                        .map(PropertyEntries::SGroupBondListEntry)
                        .parse_next(input),
                    b"M  SMT" => sgroup_subscript_entry
                        .map(PropertyEntries::SGroupSubscriptEntry)
                        .parse_next(input),
                    b"M  SDS" => sgroup_expansion_entries
                        .map(PropertyEntries::SGroupExpansionEntries)
                        .parse_next(input),
                    b"M  SPA" => sgroup_parent_atom_entries
                        .map(PropertyEntries::SGroupParentAtomEntry)
                        .parse_next(input),
                    b"M  CRS" => sgroup_correspondence_entry
                        .map(PropertyEntries::SGroupCorrespondenceEntry)
                        .parse_next(input),
                    b"M  SDI" => sgroup_display_info_entry
                        .map(PropertyEntries::SGroupDisplayInfoEntry)
                        .parse_next(input),
                    b"M  SBV" => sgroup_connecting_bond_entry
                        .map(PropertyEntries::SGroupConnectingBondEntry)
                        .parse_next(input),
                    b"M  SDT" => sgroup_data_description_entry
                        .map(PropertyEntries::SGroupDataDescriptionEntry)
                        .parse_next(input),
                    b"M  SDD" => sgroup_data_display_entry(skip_unused_fields)
                        .map(PropertyEntries::SGroupDataDisplayEntry)
                        .parse_next(input),
                    b"M  SCD" => sgroup_data_continuation_entry
                        .map(PropertyEntries::SGroupDataEntry)
                        .parse_next(input),
                    b"M  SED" => sgroup_data_end_entry
                        .map(PropertyEntries::SGroupDataEntry)
                        .parse_next(input),
                    b"M  SPL" => sgroup_hierarchy_entries
                        .map(PropertyEntries::SGroupHierarchyEntries)
                        .parse_next(input),
                    b"M  SNC" => sgroup_component_entries
                        .map(PropertyEntries::SGroupComponentEntries)
                        .parse_next(input),
                    _ => unreachable!(),
                }
            }
            tag @ (b"M  ZBO" | b"M  ZCH" | b"M  HYD") if allow_clark_extensions => match tag {
                b"M  ZBO" => bond_order_override_entries
                    .map(PropertyEntries::BondOrderOverrideEntries)
                    .parse_next(input),
                b"M  ZCH" => atom_charge_overrides_entries
                    .map(PropertyEntries::AtomChargeOverrideEntries)
                    .parse_next(input),
                b"M  HYD" => atom_hydrogen_count_entries
                    .map(PropertyEntries::AtomHydrogenCountEntries)
                    .parse_next(input),
                _ => unreachable!(),
            },
            tag @ (b"M  ZZC" | b"M  MRV") if allow_editor_extensions => match tag {
                b"M  ZZC" => chemsketch_label_entry
                    .map(PropertyEntries::ChemSketchLabelEntry)
                    .parse_next(input),
                b"M  MRV" if allow_queries => marvin_smarts_pattern_entry
                    .map(PropertyEntries::MarvinSmartsPatternEntry)
                    .parse_next(input),
                _ => unreachable!(),
            },
            _ => unreachable!(),
        };
        result.map_err(ErrMode::cut)
    }
}

/// Create atom alias property from two lines.
/// A  aaa
/// x...
/// aaa: atom index, x...: alias text
pub(super) fn parse_atom_alias_input(
    first_line: &[u8],
    second_line: &[u8],
) -> PResult<PropertyEntries> {
    let mut input = Input::new(first_line);
    let atom_index = preceded(b"A  ", fixed_width_int_minus1::<u32>(3)).parse_next(&mut input)?;

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
pub(super) fn parse_legacy_group_abbreviation_input(
    first_line: &[u8],
    second_line: &[u8],
) -> PResult<PropertyEntries> {
    let mut input = Input::new(first_line);
    let (atom_index1, atom_index2) = preceded(
        b"G  ",
        (
            fixed_width_int_minus1::<u32>(3),
            fixed_width_int_minus1::<u32>(3),
        ),
    )
    .parse_next(&mut input)?;

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
fn atom_value_entry(input: &mut Input<'_>) -> PResult<AtomValueEntry> {
    (fixed_width_int_minus1::<u32>(3), preceded(b' ', rest))
        .map(|(atom_index, value): (u32, &[u8])| AtomValueEntry {
            atom_index,
            value: value.to_str_lossy().into_owned(),
        })
        .parse_next(input)
}

fn counted<'inp, N, O, Count, Item>(
    mut count: Count,
    mut item: Item,
) -> impl Parser<Input<'inp>, Vec<O>, ErrMode<InputError>>
where
    N: Into<usize>,
    Count: Parser<Input<'inp>, N, ErrMode<InputError>>,
    Item: Parser<Input<'inp>, O, ErrMode<InputError>>,
{
    move |input: &mut Input<'inp>| {
        let count = count.parse_next(input)?.into();
        let mut values = Vec::with_capacity(count);
        for _ in 0..count {
            values.push(item.parse_next(input)?);
        }
        Ok(values)
    }
}

/// Parse charge property entries.
/// nn8 aaa vvv ...
/// vvv: -15..= 15.
fn charge_entries(input: &mut Input<'_>) -> PResult<Vec<ChargeEntry>> {
    counted(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        (
            preceded(b' ', fixed_width_int_minus1::<u32>(3)),
            preceded(b' ', fixed_width_int_in_range::<i8, _>(3, -15..=15)),
        )
            .map(|(atom_index, charge)| ChargeEntry { atom_index, charge }),
    )
    .parse_next(input)
}

/// Parse radical property entries.
/// nn8 aaa vvv ...
/// vvv: 0..= 3: 0 = no radical, 1 = singlet (:), 2 = doublet (. or ^), 3 = triplet (^^).
fn radical_entries(input: &mut Input<'_>) -> PResult<Vec<RadicalEntry>> {
    counted(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        (
            preceded(b' ', fixed_width_int_minus1::<u32>(3)),
            preceded(b' ', fixed_width_int_in_range::<u8, _>(3, 0..=3)),
        )
            .map(|(atom_index, radical_type)| RadicalEntry {
                atom_index,
                radical_type,
            }),
    )
    .parse_next(input)
}

/// Parse isotope property entries.
/// nn8 aaa vvv ...
/// vvv: isotope mass number (not difference)
/// Difference between the isotope mass number and reference isotope mass number
/// should be in the range -18..=12.
fn isotope_entries(input: &mut Input<'_>) -> PResult<Vec<IsotopeEntry>> {
    counted(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        (
            preceded(b' ', fixed_width_int_minus1::<u32>(3)),
            preceded(b' ', fixed_width_int::<u32>(3)),
        )
            .map(|(atom_index, mass)| IsotopeEntry { atom_index, mass }),
    )
    .parse_next(input)
}

/// Parse ring bond count property entries.
/// M  RBCnn8 aaa vvv ...
/// vvv: Ring bond count (-2 = as drawn (r*), -1 = no ring bonds (r0), 0 = off, 2 = r2, 3 = r3, 4 = r4+)
fn ring_bond_count_entries(input: &mut Input<'_>) -> PResult<Vec<RingBondCountEntry>> {
    counted(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        (
            preceded(b' ', fixed_width_int_minus1::<u32>(3)),
            preceded(b' ', fixed_width_int_in_range::<i8, _>(3, -2..=4)),
        )
            .map(|(atom_index, ring_bond_count)| RingBondCountEntry {
                atom_index,
                ring_bond_count,
            }),
    )
    .parse_next(input)
}

/// Parse substitution count property entries.
/// M  SUBnn8 aaa vvv ...
/// vvv: Substitution count (-2 = as drawn (s*), -1 = no substitution (s0), 0 = off, 1-5 = s1-s5,
/// 6 = s6+)
fn substitution_count_entries(input: &mut Input<'_>) -> PResult<Vec<SubstitutionCountEntry>> {
    counted(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        (
            preceded(b' ', fixed_width_int_minus1::<u32>(3)),
            preceded(b' ', fixed_width_int_in_range::<i8, _>(3, -2..=15)),
        )
            .map(|(atom_index, substitution_count)| SubstitutionCountEntry {
                atom_index,
                substitution_count,
            }),
    )
    .parse_next(input)
}

/// Parse unsaturated atom property entries.
/// M  UNSnn8 aaa vvv ...
/// vvv: Unsaturated flag (0 = off, 1 = on)
fn unsaturated_atom_entries(input: &mut Input<'_>) -> PResult<Vec<UnsaturatedAtomEntry>> {
    counted(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        (
            preceded(b' ', fixed_width_int_minus1::<u32>(3)),
            preceded(b' ', fixed_width_int_in_range::<u8, _>(3, 0..=1)),
        )
            .map(|(atom_index, unsaturated)| UnsaturatedAtomEntry {
                atom_index,
                unsaturated,
            }),
    )
    .parse_next(input)
}

/// Parse link atom property entries.
/// M  LINnn4 aaa vvv bbb ccc
/// nn4: Count (max 4), aaa: Atom index, vvv: Upper repeat count (>= 2, lower repeat count is 1),
/// bbb/ccc: Substituent indices (can be 0)
fn link_atom_entries(input: &mut Input<'_>) -> PResult<Vec<LinkAtomEntry>> {
    counted(
        fixed_width_int_in_range::<u8, _>(3, 1..=4),
        (
            preceded(b' ', fixed_width_int_minus1::<u32>(3)),
            preceded(b' ', fixed_width_int_in_range::<u8, _>(3, 2..=255)),
            preceded(b' ', fixed_width_int_minus1::<u32>(3)),
            opt(preceded(b' ', fixed_width_int_minus1::<u32>(3))),
        )
            .map(
                |(atom_index, repeat_count, subs_index1, subs_index2)| LinkAtomEntry {
                    atom_index,
                    repeat_count,
                    subs_index1,
                    subs_index2,
                },
            ),
    )
    .parse_next(input)
}

/// Parse atom list property entry.
/// M  ALS aaannn e 11112222333344445555...
/// aaa: Atom number, nnn: Number of entries (16 max), e: Exclusion (T/F), 1111: 4-char symbols
fn atom_list_entry(input: &mut Input<'_>) -> PResult<AtomListEntry> {
    let atom_index = preceded(b' ', fixed_width_int_minus1::<u32>(3)).parse_next(input)?;
    let count = fixed_width_int_in_range::<usize, _>(3, 1..=16).parse_next(input)?;
    let column = input.current_token_start() + 1;
    let exclusion_byte: &[u8] = delimited(b' ', take(1usize), b' ').parse_next(input)?;
    let exclusion = match exclusion_byte {
        b"T" => true,
        b"F" | b" " => false,
        _ => return Err(ErrMode::Backtrack(InputError::at_column(column))),
    };

    let mut elements = Vec::with_capacity(count);
    for _ in 0..count {
        let column = input.current_token_start();
        let element = fixed_width_element_partial(4)
            .parse_next(input)?
            .ok_or(ErrMode::Backtrack(InputError::at_column(column)))?;
        elements.push(element);
    }

    Ok(AtomListEntry {
        atom_index,
        exclusion,
        elements,
    })
}

/// Parse attachment point property entries.
/// Atom aaa is typically on ordinary atom, does not have to be a RGroup (opposite of AAL)
/// Attachment point does not appear in the atom list
/// M  APOnn2 aaa vvv ...
/// nn2: count (max 2), aaa: atom index, vvv: attachment type (0-3)
/// 0 = no attachment, 1 = first attachment point, 2 = second attachment point, 3 = both attachment points
fn attachment_point_entries(input: &mut Input<'_>) -> PResult<Vec<AttachmentPointEntry>> {
    counted(
        fixed_width_int_in_range::<u8, _>(3, 1..=2),
        (
            preceded(b' ', fixed_width_int_minus1::<u32>(3)),
            preceded(b' ', fixed_width_int_in_range::<u8, _>(3, 0..=3)),
        )
            .map(|(atom_index, attachment_type)| AttachmentPointEntry {
                atom_index,
                attachment_type,
            }),
    )
    .parse_next(input)
}

/// Parse atom attachment order entry.
/// M  AAL aaann2 111 v1v 222 v2v ...
/// Atom aaa refers to an RGroup, atoms 111, 222 are ordinary atoms (opposite of APO)
/// aaa: atom index, n2: pair count (max 2), 111/222: neighbor indices, v1v/v2v: attachment orders
fn atom_attachment_order_entry(input: &mut Input<'_>) -> PResult<AtomAttachmentOrderEntry> {
    (
        preceded(b' ', fixed_width_int_minus1::<u32>(3)),
        counted(
            fixed_width_int_in_range::<u8, _>(3, 1..=2),
            (
                preceded(b' ', fixed_width_int_minus1::<u32>(3)),
                preceded(b' ', fixed_width_int_in_range::<u8, _>(3, 1..=2)),
            ),
        ),
    )
        .map(|(atom_index, attachments)| AtomAttachmentOrderEntry {
            atom_index,
            attachments,
        })
        .parse_next(input)
}

/// Parse RGroup label property entries.
/// M  RGPnn8 aaa rrr ...
/// aaa: atom index, rrr: RGroup label (0 = no label, 1-32 in CTab spec, 1-999 in RDKit)
fn rgroup_label_entries(input: &mut Input<'_>) -> PResult<Vec<RGroupLabelEntry>> {
    counted(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        (
            preceded(b' ', fixed_width_int_minus1::<u32>(3)),
            preceded(b' ', fixed_width_int::<u32>(3)),
        )
            .map(|(atom_index, label)| RGroupLabelEntry { atom_index, label }),
    )
    .parse_next(input)
}

// Parse RGroup logic entry.
/// M  LOGnn1 rrr iii hhh ooo...
/// nn1: count (max 1), rrr: RGroup label (1-32 in CTab spec, 1-999 in RDKit)
/// iii: Number of dependent Rgroup (IF rrr THEN iii)
/// hhh: REstH property of rrr (0 (default)=off, 1=on): RGroup or H atom
/// ooo...: Range of RGroup occurrence required: n=exactly n, n-m=from n through m (inclusive),
/// >n=greater n, <n=fewer than n, blank (default): > 0.
/// > Any non-contradictory combination of the preceding values is allowed, separated by commas.
fn rgroup_logic_entry(input: &mut Input<'_>) -> PResult<RGroupLogicEntry> {
    let _count = fixed_width_int_in_range::<u8, _>(3, 1..=1).parse_next(input)?;
    let label = preceded(b' ', fixed_width_int_in_range::<u32, _>(3, 1..=999)).parse_next(input)?;
    let dependent_label = preceded(b' ', fixed_width_int::<u32>(3))
        .map(|label| (label != 0).then_some(label))
        .parse_next(input)?;
    let rgroup_or_h = cond(
        input.len() >= 4,
        preceded(b' ', fixed_width_int_in_range::<u8, _>(3, 0..=1)),
    )
    .map(|value| value.is_some_and(|value| value != 0))
    .parse_next(input)?;
    let occurrence = cond(!input.is_empty(), rgroup_occurrences())
        .map(|value| value.unwrap_or_else(|| vec![RGroupOccurrence::GreaterThan(0)]))
        .parse_next(input)?;

    Ok(RGroupLogicEntry {
        label,
        dependent_label,
        rgroup_or_h,
        occurrence,
    })
}

/// Parse SGroup type entries.
/// M  STYnn8 sss ttt ...
/// sss: SGroup index, ttt: SGroup type (3-character string)
fn sgroup_type_entries(input: &mut Input<'_>) -> PResult<Vec<SGroupTypeEntry>> {
    counted(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        (
            preceded(b' ', fixed_width_int_minus1::<u32>(3)),
            preceded(b' ', sgroup_type),
        )
            .map(|(sgroup_index, sgroup_type)| SGroupTypeEntry {
                sgroup_index,
                sgroup_type,
            }),
    )
    .parse_next(input)
}

/// Parse SGroup subtype entries.
/// M  SSTnn8 sss ttt ...
/// sss: SGroup index, ttt: SGroup subtype (3-character string)
fn sgroup_subtype_entries(input: &mut Input<'_>) -> PResult<Vec<SGroupSubtypeEntry>> {
    counted(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        (
            preceded(b' ', fixed_width_int_minus1::<u32>(3)),
            preceded(b' ', sgroup_subtype),
        )
            .map(|(sgroup_index, sgroup_subtype)| SGroupSubtypeEntry {
                sgroup_index,
                sgroup_subtype,
            }),
    )
    .parse_next(input)
}

/// Parse SGroup label entries.
/// M  SLBnn8 sss vvv ...
/// sss: SGroup index, vvv: integer label is from 1-512
fn sgroup_label_entries(input: &mut Input<'_>) -> PResult<Vec<SGroupLabelEntry>> {
    counted(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        (
            preceded(b' ', fixed_width_int_minus1::<u32>(3)),
            preceded(b' ', fixed_width_int_in_range::<u32, _>(3, 1..=512)),
        )
            .map(|(sgroup_index, label)| SGroupLabelEntry {
                sgroup_index,
                label,
            }),
    )
    .parse_next(input)
}

/// Parse SGroup connectivity entries.
/// M  SCNnn8 sss ttt ...
/// sss: SGroup index, ttt: SGroup connectivity (2-character string), left-justified
fn sgroup_connectivity_entries(input: &mut Input<'_>) -> PResult<Vec<SGroupConnectivityEntry>> {
    counted(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        (
            preceded(b' ', fixed_width_int_minus1::<u32>(3)),
            preceded(b' ', sgroup_connectivity),
        )
            .map(|(sgroup_index, connectivity)| SGroupConnectivityEntry {
                sgroup_index,
                connectivity,
            }),
    )
    .parse_next(input)
}

/// Parse SGroup expansion entries.
/// M SDS EXPn15 sss ...
/// sss: SGroup index, n15: count (max 15)
fn sgroup_expansion_entries(input: &mut Input<'_>) -> PResult<Vec<SGroupExpansionEntry>> {
    preceded(
        b" EXP".as_slice(),
        counted(
            fixed_width_int_in_range::<u8, _>(3, 1..=15),
            preceded(b' ', fixed_width_int_minus1::<u32>(3))
                .map(|sgroup_index| SGroupExpansionEntry { sgroup_index }),
        ),
    )
    .parse_next(input)
}

/// Parse SGroup atom list entry.
/// M  SAL sssn15 aaa ...
/// sss: SGroup index, n15: count (max 15), aaa: atom indices (each 4 chars: " aaa")
fn sgroup_atom_list_entry(input: &mut Input<'_>) -> PResult<SGroupAtomListEntry> {
    (
        preceded(b' ', fixed_width_int_minus1::<u32>(3)),
        counted(
            fixed_width_int_in_range::<u8, _>(3, 1..=15),
            preceded(b' ', fixed_width_int_minus1::<u32>(3)),
        ),
    )
        .map(|(sgroup_index, atom_indices)| SGroupAtomListEntry {
            sgroup_index,
            atom_indices,
        })
        .parse_next(input)
}

/// Parse SGroup bond list entry.
/// M  SBL sssn15 bbb ...
/// sss: SGroup index (3 chars), n: count (3 chars), bbb: bond indices (each 4 chars: " bbb")
fn sgroup_bond_list_entry(input: &mut Input<'_>) -> PResult<SGroupBondListEntry> {
    (
        preceded(b' ', fixed_width_int_minus1::<u32>(3)),
        counted(
            fixed_width_int_in_range::<u8, _>(3, 1..=15),
            preceded(b' ', fixed_width_int_minus1::<u32>(3)),
        ),
    )
        .map(|(sgroup_index, bond_indices)| SGroupBondListEntry {
            sgroup_index,
            bond_indices,
        })
        .parse_next(input)
}

/// Parse SGroup parent atom entries.
/// M  SPA sssn15 aaa ...
/// sss: SGroup index, n15: count (max 15), aaa: atom indices (each 4 chars: " aaa")
fn sgroup_parent_atom_entries(input: &mut Input<'_>) -> PResult<SGroupParentAtomEntry> {
    (
        preceded(b' ', fixed_width_int_minus1::<u32>(3)),
        counted(
            fixed_width_int_in_range::<u8, _>(3, 1..=15),
            preceded(b' ', fixed_width_int_minus1::<u32>(3)),
        ),
    )
        .map(|(sgroup_index, atom_indices)| SGroupParentAtomEntry {
            sgroup_index,
            atom_indices,
        })
        .parse_next(input)
}

/// Parse SGroup subscript entry.
/// M  SMT sss m...
/// sss: SGroup index, m: subscript text
/// For multiple groups, m... is the text representation of the multiple group multiplier.For superatoms,
/// m... is the text of the superatom label.)
fn sgroup_subscript_entry(input: &mut Input<'_>) -> PResult<SGroupSubscriptEntry> {
    let sgroup_index = preceded(b' ', fixed_width_int_minus1::<u32>(3)).parse_next(input)?;
    let value: &[u8] = preceded(b' ', rest).parse_next(input)?;
    let mut multiplier_input = Input::new(value);
    let multiplier = sgroup_multiplier.parse_next(&mut multiplier_input).ok();
    let mut subscript_input = Input::new(value);
    let subscript = sgroup_subscript.parse_next(&mut subscript_input).ok();

    Ok(SGroupSubscriptEntry {
        sgroup_index,
        multiplier,
        subscript,
    })
}

/// Parse SGroup correspondence entry.
/// M  CRS sssnn6 bb1 bb2 bb3
/// sss: SGroup index, nn6: count (max 6), bb1-bb3: bond indices
/// bb1, bb2: Crossing bonds that share a common bracket
/// bb3: Crossing bond in repeating unit that connects to bond bb1
fn sgroup_correspondence_entry(input: &mut Input<'_>) -> PResult<SGroupCorrespondenceEntry> {
    (
        preceded(b' ', fixed_width_int_minus1::<u32>(3)),
        counted(
            fixed_width_int_in_range::<usize, _>(3, 1..=6),
            preceded(b' ', fixed_width_int_minus1::<u32>(3)),
        ),
    )
        .map(|(sgroup_index, bond_indices)| SGroupCorrespondenceEntry {
            sgroup_index,
            bond_indices,
        })
        .parse_next(input)
}

/// Parse SGroup display info entry.
/// M  SDI sssnn4 x1 y1 x2 y2
/// sss: SGroup index, nn4: count (max 4), x1, y1: opening bracket, x2, y2: closing bracket (f10.4)
fn sgroup_display_info_entry(input: &mut Input<'_>) -> PResult<SGroupDisplayInfoEntry> {
    (
        preceded(b' ', fixed_width_int_minus1::<u32>(3)),
        counted(
            fixed_width_int_in_range::<usize, _>(3, 1..=4),
            fixed_width_float_f10_4::<f64>(),
        ),
    )
        .map(|(sgroup_index, bracket_coords)| SGroupDisplayInfoEntry {
            sgroup_index,
            bracket_coords,
        })
        .parse_next(input)
}

/// Parse SGroup connecting bond entry.
/// M  SBV sss bb1 x1 y1
/// sss: SGroup index, bb1: bond index, x1, y1: bond vector (f10.4)
fn sgroup_connecting_bond_entry(input: &mut Input<'_>) -> PResult<SGroupConnectingBondEntry> {
    (
        preceded(b' ', fixed_width_int_minus1::<u32>(3)),
        preceded(b' ', fixed_width_int_minus1::<u32>(3)),
        fixed_width_float_f10_4::<f64>(),
        fixed_width_float_f10_4::<f64>(),
    )
        .map(
            |(sgroup_index, bond_index, x1, y1)| SGroupConnectingBondEntry {
                sgroup_index,
                bond_index,
                bond_vector: (x1, y1),
            },
        )
        .parse_next(input)
}

/// Parse SGroup data field description entry.
/// M  SDT sss fff..fffgghh...hhhiijjj..,
/// sss: SGroup index, fff..fff: field name (30 chars), gg: field type (2 chars, F=formatted, N=numeric, T=text)
/// hh: field units or format (20 chars), ii: Query identifier (MQ: MACCS-II, IQ: ISIS, PQ: program code)
/// jjj: data query operator
fn sgroup_data_description_entry(input: &mut Input<'_>) -> PResult<SGroupDataDescriptionEntry> {
    (
        preceded(b' ', fixed_width_int_minus1::<u32>(3)),
        preceded(b' ', fixed_width_str_partial(30)),
        sgroup_data_type,
        fixed_width_str_partial(20),
        fixed_width_str_partial(2),
        fixed_width_str_partial(1),
    )
        .map(
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
        .parse_next(input)
}

/// Parse SGroup data display entry.
/// M  SDD sss xxxxx.xxxxyyyyy.yyyy eeefgh i jjjkkk ll m  noo
/// sss: SGroup index, x, y: coordinates (f10.4), eee: skipped, f: data display (A=attached, D=detached)
/// g: absolute, relative placement (A=absolute, R=relative), h: display units (" "=none, U=display units)
/// i: skipped, jjj: number of characters to display (1-999 or ALL), kkk: number of lines to display (unused, always 1)
/// ll: skipped, m: tag character (if non-blank),
/// nn, oo: skipped (nn is MACCS-II only in the CTFile spec and refers to the second position, however example files
/// show numerical data in the first position as well).
fn sgroup_data_display_entry<'inp>(
    skip_unused_fields: bool,
) -> impl Parser<Input<'inp>, SGroupDataDisplayEntry, ErrMode<InputError>> {
    (
        preceded(b' ', fixed_width_int_minus1::<u32>(3)),
        preceded(b' ', fixed_width_float_f10_4::<f64>()),
        terminated(fixed_width_float_f10_4::<f64>(), b' '),
        preceded(
            fixed_width_unused(3, skip_unused_fields),
            sgroup_data_display_type,
        ),
        sgroup_data_display_placement,
        sgroup_data_display_units,
        preceded(
            fixed_width_unused(3, skip_unused_fields),
            sgroup_data_display_chars,
        ),
        delimited(
            fixed_width_unused(7, skip_unused_fields),
            fixed_width_partial(1, <u8 as IntParser>::parse, true),
            opt((b' ', fixed_width_unused(2, skip_unused_fields))),
        ),
    )
        .map(
            |(
                sgroup_index,
                x,
                y,
                display_type,
                display_placement,
                display_units,
                display_chars,
                display_tag,
            )| SGroupDataDisplayEntry {
                sgroup_index,
                coords: (x, y),
                display_type,
                display_placement,
                display_units,
                display_chars,
                display_tag,
            },
        )
}

/// Parse SGroup data continuation entry.
/// M  SCD sss d...
/// sss: SGroup index, d...: data content (max 69 chars)
fn sgroup_data_continuation_entry(input: &mut Input<'_>) -> PResult<SGroupDataEntry> {
    (
        preceded(b' ', fixed_width_int_minus1::<u32>(3)),
        opt(preceded(b' ', rest)),
    )
        .map(|(sgroup_index, data_content): (u32, Option<&[u8]>)| {
            let data_content = data_content.unwrap_or_default();
            SGroupDataEntry::Continuation {
                sgroup_index,
                data_content: data_content[..data_content.len().min(69)]
                    .to_str_lossy()
                    .into_owned(),
            }
        })
        .parse_next(input)
}

/// Parse SGroup data end entry.
/// M  SED sss d... OR
/// M  SED
/// sss: SGroup index, d...: data content (max 69 chars)
/// The data content is the same as the data content in the SCD entry.
/// The second form is used to indicate the end of the data content, in which case it should be empty.
fn sgroup_data_end_entry(input: &mut Input<'_>) -> PResult<SGroupDataEntry> {
    let sgroup_index = preceded(b' ', fixed_width_int_minus1::<u32>(3)).parse_next(input)?;
    match opt(preceded(b' ', rest)).parse_next(input)? {
        Some(data_content) => Ok(SGroupDataEntry::EndWithData {
            sgroup_index,
            data_content: data_content[..data_content.len().min(69)]
                .to_str_lossy()
                .into_owned(),
        }),
        None => Ok(SGroupDataEntry::EndBlank { sgroup_index }),
    }
}

/// Parse SGroup hierarchy entries.
/// M  SPLnn8 sss ppp ...
/// sss: SGroup index, ppp: Parent SGroup index
fn sgroup_hierarchy_entries(input: &mut Input<'_>) -> PResult<Vec<SGroupHierarchyEntry>> {
    counted(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        (
            preceded(b' ', fixed_width_int_minus1::<u32>(3)),
            preceded(b' ', fixed_width_int_minus1::<u32>(3)),
        )
            .map(|(sgroup_index, parent_sgroup_index)| SGroupHierarchyEntry {
                sgroup_index,
                parent_sgroup_index,
            }),
    )
    .parse_next(input)
}

/// Parse SGroup component number entries.
/// M  SNCnn8 sss ccc ...
/// sss: SGroup index, ccc: Component number
fn sgroup_component_entries(input: &mut Input<'_>) -> PResult<Vec<SGroupComponentEntry>> {
    counted(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        (
            preceded(b' ', fixed_width_int_minus1::<u32>(3)),
            preceded(b' ', fixed_width_int::<u32>(3)),
        )
            .map(|(sgroup_index, component_number)| SGroupComponentEntry {
                sgroup_index,
                component_number,
            }),
    )
    .parse_next(input)
}

/// Parse zero bond order entries.
/// M  ZBOnn8 bbb vvv ...
/// bbb: bond index, vvv: bond vector override (>= 0, limited to 6 here)
fn bond_order_override_entries(input: &mut Input<'_>) -> PResult<Vec<BondOrderOverrideEntry>> {
    counted(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        (
            preceded(b' ', fixed_width_int_minus1::<u32>(3)),
            preceded(b' ', fixed_width_int_in_range::<u8, _>(3, 0..=6)),
        )
            .verify_map(|(bond_index, bond_order)| {
                BondOrder::from_value(bond_order).map(|bo| BondOrderOverrideEntry {
                    bond_index,
                    bond_order: bo,
                })
            }),
    )
    .parse_next(input)
}

/// Parse zero atom charge entries.
/// M  ZCHnn8 aaa ccc ...
/// aaa: atom index, ccc: atom charge override (-8..=8)
fn atom_charge_overrides_entries(input: &mut Input<'_>) -> PResult<Vec<AtomChargeOverrideEntry>> {
    counted(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        (
            preceded(b' ', fixed_width_int_minus1::<u32>(3)),
            preceded(b' ', fixed_width_int_in_range::<i8, _>(3, -8..=8)),
        )
            .map(|(atom_index, charge)| AtomChargeOverrideEntry { atom_index, charge }),
    )
    .parse_next(input)
}

/// Parse atom explicit hydrogen count entries.
/// M  HYDnn8 aaa hhh ...
/// aaa: atom index, hhh: atom explicit hydrogen count (>= 0, limited to 8 here, -1 = no override)
fn atom_hydrogen_count_entries(input: &mut Input<'_>) -> PResult<Vec<AtomHydrogenCountEntry>> {
    counted(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        (
            preceded(b' ', fixed_width_int_minus1::<u32>(3)),
            preceded(b' ', fixed_width_int_in_range::<i8, _>(3, -1..=8)),
        )
            .map(|(atom_index, hydrogen_count)| AtomHydrogenCountEntry {
                atom_index,
                hydrogen_count: if hydrogen_count == -1 {
                    None
                } else {
                    Some(hydrogen_count as u8)
                },
            }),
    )
    .parse_next(input)
}

/// Parse ACD/ChemSketch label entry
/// M  ZZC aaa x...
/// aaa: atom index, x: text label
fn chemsketch_label_entry(input: &mut Input<'_>) -> PResult<ChemSketchLabelEntry> {
    (
        preceded(b' ', fixed_width_int_minus1::<u32>(3)),
        preceded(b' ', rest),
    )
        .map(|(atom_index, label): (u32, &[u8])| ChemSketchLabelEntry {
            atom_index,
            label: label.to_str_lossy().into_owned(),
        })
        .parse_next(input)
}

/// Parse Marvin SMARTS pattern entry
/// M  MRV SMA aaa p...
/// aaa: atom index, p: SMARTS pattern
fn marvin_smarts_pattern_entry(input: &mut Input<'_>) -> PResult<MarvinSmartsPatternEntry> {
    (
        preceded(b" SMA ".as_slice(), fixed_width_int_minus1::<u32>(3)),
        preceded(b' ', rest),
    )
        .map(
            |(atom_index, smarts_pattern): (u32, &[u8])| MarvinSmartsPatternEntry {
                atom_index,
                smarts_pattern: smarts_pattern.to_str_lossy().into_owned(),
            },
        )
        .parse_next(input)
}

#[cfg(test)]
mod tests;
