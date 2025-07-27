//! Parsers for CTab property lines.

use nom::bytes::complete::{tag, take};
use nom::character::complete::{line_ending, not_line_ending};
use nom::combinator::{all_consuming, map, map_opt, map_parser, map_res, opt, success};
use nom::error;
use nom::multi::{length_count, many_m_n};
use nom::sequence::{delimited, preceded, terminated};
use nom::Parser;

use crate::io::ctab::rgroup::RGroupOccurrence;
use crate::io::ctab::sgroup::{SGroup, SGroupConnectivity, SGroupSubtype, SGroupType};

use super::utils::{
    fixed_width_element_partial, fixed_width_float, fixed_width_int, fixed_width_int_in_range,
    fixed_width_int_minus1, fixed_width_int_partial, rgroup_occurrences,
};
use umol_data::Element;

#[derive(Debug, Clone, PartialEq)]
pub struct AtomAliasEntry {
    pub atom_index: usize,
    pub alias: String,
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
    pub exclusion: bool,                   // T = NOT list, F = normal list
    pub elements: Vec<umol_data::Element>, // Converted from 4-char symbols
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
    pub subscript: String,
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

/// An enum representing a parsed property modification, containing the raw data.
/// This avoids allocating a new Vec for every single property line in a file.
#[derive(Debug, Clone, PartialEq)]
pub enum PropertyEntries {
    AtomAliasEntry(AtomAliasEntry),
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
    // Superatom bond and vector [SBV]
    // Data sgroup fields [SDT]
    // Data sgroup display information [SDD]
    // Data sgroup data [SCD] / [SED]
    // Sgroup hierarchy [SPL]
    // Sgroup component numbers [SNC]
    End,
}

/// Parse a legacy atom list entry (e.g., "  1  2 C  N   ")
pub fn legacy_atom_list_input<'a>(
) -> impl Parser<&'a [u8], Output = PropertyEntries, Error = error::Error<&'a [u8]>> {
    |input: &'a [u8]| {
        // Parse atom index (3 chars)
        let (remaining, atom_index) = fixed_width_int_minus1::<usize>(3).parse(input)?;

        // Parse exclusion flag (1 char)
        let (remaining, exclusion_byte) =
            delimited(tag(" "), take(1usize), tag("    ")).parse(remaining)?;
        let exclusion = match exclusion_byte {
            b"T" => true,
            b"F" | b" " => false,
            _ => {
                return Err(nom::Err::Error(error::Error::new(
                    input,
                    error::ErrorKind::Tag,
                )))
            }
        };

        // Parse atom list (at most 5, each 3 chars, delimited by spaces)
        let (remaining, elements) = length_count(
            fixed_width_int_in_range::<u8, _>(3, 1..=5),
            preceded(
                tag(" "),
                map_opt(
                    fixed_width_int_partial::<u8>(3),
                    Element::from_atomic_number,
                ),
            ),
        )
        .parse(remaining)?;

        Ok((
            remaining,
            PropertyEntries::AtomListEntry(AtomListEntry {
                atom_index,
                exclusion,
                elements,
            }),
        ))
    }
}

/// Parse property line (standard properties only - no queries)
pub fn property_input_standard<'a>(
    input: &'a [u8],
) -> nom::IResult<&'a [u8], PropertyEntries, error::Error<&'a [u8]>> {
    if input.len() < 3 {
        return Err(nom::Err::Error(error::Error::new(
            input,
            error::ErrorKind::Eof,
        )));
    }

    // Handle A and V lines (different format from M lines)
    match &input[0..3] {
        b"A  " => atom_alias_entry()
            .parse(&input[3..])
            .map(|(i, o)| (i, PropertyEntries::AtomAliasEntry(o))),
        b"V  " => atom_value_entry()
            .parse(&input[3..])
            .map(|(i, o)| (i, PropertyEntries::AtomValueEntry(o))),
        _ => {
            // Handle M lines
            if input.len() < 6 {
                return Err(nom::Err::Error(error::Error::new(
                    input,
                    error::ErrorKind::Eof,
                )));
            }
            let (rest, tag) = take(6u8)(input)?;
            match tag {
                b"M  CHG" => charge_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::ChargeEntries(o))),
                b"M  RAD" => radical_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::RadicalEntries(o))),
                b"M  ISO" => isotope_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::IsotopeEntries(o))),
                b"M  STY" => sgroup_type_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SGroupTypeEntries(o))),
                b"M  SST" => sgroup_subtype_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SGroupSubtypeEntries(o))),
                b"M  SLB" => sgroup_label_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SGroupLabelEntries(o))),
                b"M  SAL" => sgroup_atom_list_entry()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SGroupAtomListEntry(o))),
                b"M  SBL" => sgroup_bond_list_entry()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SGroupBondListEntry(o))),
                b"M  SMT" => sgroup_subscript_entry()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SGroupSubscriptEntry(o))),
                b"M  END" => success(PropertyEntries::End).parse(rest),
                _ => Err(nom::Err::Error(error::Error::new(
                    input,
                    error::ErrorKind::Tag,
                ))),
            }
        }
    }
}

/// Parse property line (all properties including queries)
pub fn property_input<'a>(
    input: &'a [u8],
) -> nom::IResult<&'a [u8], PropertyEntries, error::Error<&'a [u8]>> {
    if input.len() < 3 {
        return Err(nom::Err::Error(error::Error::new(
            input,
            error::ErrorKind::Eof,
        )));
    }

    // Handle A and V lines (different format from M lines)
    match &input[0..3] {
        b"A  " => atom_alias_entry()
            .parse(&input[3..])
            .map(|(i, o)| (i, PropertyEntries::AtomAliasEntry(o))),
        b"V  " => atom_value_entry()
            .parse(&input[3..])
            .map(|(i, o)| (i, PropertyEntries::AtomValueEntry(o))),
        _ => {
            // Handle M lines
            if input.len() < 6 {
                return Err(nom::Err::Error(error::Error::new(
                    input,
                    error::ErrorKind::Eof,
                )));
            }
            let (rest, tag) = take(6u8)(input)?;
            match tag {
                b"M  CHG" => charge_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::ChargeEntries(o))),
                b"M  RAD" => radical_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::RadicalEntries(o))),
                b"M  ISO" => isotope_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::IsotopeEntries(o))),
                b"M  RBC" => ring_bond_count_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::RingBondCountEntries(o))),
                b"M  SUB" => substitution_count_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SubstitutionCountEntries(o))),
                b"M  UNS" => unsaturated_atom_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::UnsaturatedAtomEntries(o))),
                b"M  LIN" => link_atom_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::LinkAtomEntries(o))),
                b"M  ALS" => atom_list_entry()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::AtomListEntry(o))),
                b"M  APO" => attachment_point_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::AttachmentPointEntries(o))),
                b"M  AAL" => atom_attachment_order_entry()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::AtomAttachmentOrderEntry(o))),
                b"M  RGP" => rgroup_label_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::RGroupLabelEntries(o))),
                b"M  LOG" => rgroup_logic_entry()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::RGroupLogicEntry(o))),
                b"M  STY" => sgroup_type_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SGroupTypeEntries(o))),
                b"M  SST" => sgroup_subtype_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SGroupSubtypeEntries(o))),
                b"M  SLB" => sgroup_label_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SGroupLabelEntries(o))),
                b"M  SCN" => sgroup_connectivity_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SGroupConnectivityEntries(o))),
                b"M  SDS" => sgroup_expansion_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SGroupExpansionEntries(o))),
                b"M  SAL" => sgroup_atom_list_entry()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SGroupAtomListEntry(o))),
                b"M  SBL" => sgroup_bond_list_entry()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SGroupBondListEntry(o))),
                b"M  SPA" => sgroup_parent_atom_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SGroupParentAtomEntry(o))),
                b"M  SMT" => sgroup_subscript_entry()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SGroupSubscriptEntry(o))),
                b"M  CRS" => sgroup_correspondence_entry()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SGroupCorrespondenceEntry(o))),
                b"M  SDI" => sgroup_display_info_entry()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SGroupDisplayInfoEntry(o))),
                b"M  END" => success(PropertyEntries::End).parse(rest),
                _ => Err(nom::Err::Error(error::Error::new(
                    input,
                    error::ErrorKind::Tag,
                ))),
            }
        }
    }
}

/// Parse atom alias entry.
/// A  aaa
/// x..
/// aaa: atom index, x..: alias string (can contain spaces)
fn atom_alias_entry<'a>(
) -> impl Parser<&'a [u8], Output = AtomAliasEntry, Error = error::Error<&'a [u8]>> {
    map(
        (
            terminated(
                map_parser(not_line_ending, fixed_width_int_minus1::<usize>(3)),
                line_ending,
            ),
            not_line_ending,
        ),
        |(atom_index, alias_bytes)| AtomAliasEntry {
            atom_index,
            alias: String::from_utf8_lossy(alias_bytes).trim().to_string(),
        },
    )
}

/// Parse atom value entry.
/// V  aaa v..
/// aaa: atom index, v..: value string (can contain spaces)
fn atom_value_entry<'a>(
) -> impl Parser<&'a [u8], Output = AtomValueEntry, Error = error::Error<&'a [u8]>> {
    map(
        (
            fixed_width_int_minus1::<usize>(3),
            preceded(tag(" "), not_line_ending),
        ),
        |(atom_index, value_bytes)| AtomValueEntry {
            atom_index,
            value: String::from_utf8_lossy(value_bytes).trim().to_string(),
        },
    )
}

/// Parse charge property entries.
/// nn8 aaa vvv ...
/// vvv: -15..= 15.
fn charge_entries<'a>(
) -> impl Parser<&'a [u8], Output = Vec<ChargeEntry>, Error = error::Error<&'a [u8]>> {
    all_consuming(length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                ),
                map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_in_range::<i8, _>(3, -15..=15)),
                ),
            ),
            |(atom_index, charge)| ChargeEntry { atom_index, charge },
        ),
    ))
}

/// Parse radical property entries.
/// nn8 aaa vvv ...
/// vvv: 0..= 3: 0 = no radical, 1 = singlet (:), 2 = doublet (. or ^), 3 = triplet (^^).
fn radical_entries<'a>(
) -> impl Parser<&'a [u8], Output = Vec<RadicalEntry>, Error = error::Error<&'a [u8]>> {
    all_consuming(length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                ),
                map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_in_range::<u8, _>(3, 0..=3)),
                ),
            ),
            |(atom_index, radical_type)| RadicalEntry {
                atom_index,
                radical_type,
            },
        ),
    ))
}

/// Parse isotope property entries.
/// nn8 aaa vvv ...
/// vvv: isotope mass number (not difference)
/// Difference between the isotope mass number and reference isotope mass number
/// should be in the range -18..=12.
fn isotope_entries<'a>(
) -> impl Parser<&'a [u8], Output = Vec<IsotopeEntry>, Error = error::Error<&'a [u8]>> {
    all_consuming(length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                ),
                map_parser(take(4usize), preceded(tag(" "), fixed_width_int::<u32>(3))),
            ),
            |(atom_index, mass)| IsotopeEntry { atom_index, mass },
        ),
    ))
}

/// Parse ring bond count property entries.
/// M  RBCnn8 aaa vvv ...
/// vvv: Ring bond count (-2 = as drawn (r*), -1 = no ring bonds (r0), 0 = off, 2 = r2, 3 = r3, 4 = r4+)
fn ring_bond_count_entries<'a>(
) -> impl Parser<&'a [u8], Output = Vec<RingBondCountEntry>, Error = error::Error<&'a [u8]>> {
    all_consuming(length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                ),
                map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_in_range::<i8, _>(3, -2..=4)),
                ),
            ),
            |(atom_index, ring_bond_count)| RingBondCountEntry {
                atom_index,
                ring_bond_count,
            },
        ),
    ))
}

/// Parse substitution count property entries.
/// M  SUBnn8 aaa vvv ...
/// vvv: Substitution count (-2 = as drawn (s*), -1 = no substitution (s0), 0 = off, 1-5 = s1-s5,
/// 6 = s6+)
fn substitution_count_entries<'a>(
) -> impl Parser<&'a [u8], Output = Vec<SubstitutionCountEntry>, Error = error::Error<&'a [u8]>> {
    all_consuming(length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                ),
                map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_in_range::<i8, _>(3, -2..=15)),
                ),
            ),
            |(atom_index, substitution_count)| SubstitutionCountEntry {
                atom_index,
                substitution_count,
            },
        ),
    ))
}

/// Parse unsaturated atom property entries.
/// M  UNSnn8 aaa vvv ...
/// vvv: Unsaturated flag (0 = off, 1 = on)
fn unsaturated_atom_entries<'a>(
) -> impl Parser<&'a [u8], Output = Vec<UnsaturatedAtomEntry>, Error = error::Error<&'a [u8]>> {
    all_consuming(length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                ),
                map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_in_range::<u8, _>(3, 0..=1)),
                ),
            ),
            |(atom_index, unsaturated)| UnsaturatedAtomEntry {
                atom_index,
                unsaturated,
            },
        ),
    ))
}

/// Parse link atom property entries.
/// M  LINnn4 aaa vvv bbb ccc
/// nn4: Count (max 4), aaa: Atom index, vvv: Upper repeat count (>= 2, lower repeat count is 1),
/// bbb/ccc: Substituent indices (can be 0)
fn link_atom_entries<'a>(
) -> impl Parser<&'a [u8], Output = Vec<LinkAtomEntry>, Error = error::Error<&'a [u8]>> {
    all_consuming(length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=4),
        map(
            (
                map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                ),
                map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_in_range::<u8, _>(3, 2..=255)),
                ),
                map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                ),
                opt(map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                )),
            ),
            |(atom_index, repeat_count, subs_index1, subs_index2)| LinkAtomEntry {
                atom_index,
                repeat_count,
                subs_index1,
                subs_index2,
            },
        ),
    ))
}

/// Parse atom list property entry.
/// M  ALS aaannn e 11112222333344445555...
/// aaa: Atom number, nnn: Number of entries (16 max), e: Exclusion (T/F), 1111: 4-char symbols
fn atom_list_entry<'a>(
) -> impl Parser<&'a [u8], Output = AtomListEntry, Error = error::Error<&'a [u8]>> {
    |input: &'a [u8]| {
        // Parse atom index (3 chars)
        let (remaining, atom_index) =
            preceded(tag(" "), fixed_width_int_minus1::<usize>(3)).parse(input)?;

        // Parse count (3 chars, max 16)
        let (remaining, count) =
            fixed_width_int_in_range::<usize, _>(3, 1..=16).parse(remaining)?;

        // Parse exclusion flag (1 char)
        let (remaining, exclusion_byte) =
            delimited(tag(" "), take(1usize), tag(" ")).parse(remaining)?;
        let exclusion = match exclusion_byte {
            b"T" => true,
            b"F" | b" " => false,
            _ => {
                return Err(nom::Err::Error(error::Error::new(
                    input,
                    error::ErrorKind::Tag,
                )))
            }
        };

        // Parse 4-char atom symbols
        let (remaining, elements) =
            many_m_n(count, count, fixed_width_element_partial(4)).parse(remaining)?;

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
fn attachment_point_entries<'a>(
) -> impl Parser<&'a [u8], Output = Vec<AttachmentPointEntry>, Error = error::Error<&'a [u8]>> {
    all_consuming(length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=2),
        map(
            (
                map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                ),
                map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_in_range::<u8, _>(3, 0..=3)),
                ),
            ),
            |(atom_index, attachment_type)| AttachmentPointEntry {
                atom_index,
                attachment_type,
            },
        ),
    ))
}

/// Parse atom attachment order entry.
/// M  AAL aaan2 111 v1v 222 v2v ...
/// Atom aaa refers to an RGroup, atoms 111, 222 are ordinary atoms (opposite of APO)
/// aaa: atom index, n2: pair count (max 2), 111/222: neighbor indices, v1v/v2v: attachment orders
fn atom_attachment_order_entry<'a>(
) -> impl Parser<&'a [u8], Output = AtomAttachmentOrderEntry, Error = error::Error<&'a [u8]>> {
    map(
        (
            preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
            length_count(
                fixed_width_int_in_range::<u8, _>(3, 1..=2),
                (
                    map_parser(
                        take(4usize),
                        preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                    ),
                    map_parser(
                        take(4usize),
                        preceded(tag(" "), fixed_width_int_in_range::<u8, _>(3, 1..=2)),
                    ),
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
fn rgroup_label_entries<'a>(
) -> impl Parser<&'a [u8], Output = Vec<RGroupLabelEntry>, Error = error::Error<&'a [u8]>> {
    all_consuming(length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                ),
                map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_in_range::<u32, _>(3, 0..=999)),
                ),
            ),
            |(atom_index, label)| RGroupLabelEntry { atom_index, label },
        ),
    ))
}

// Parse RGroup logic entry.
/// M  LOGnn1 rrr iii hhh ooo...
/// nn1: count (max 1), rrr: RGroup label (1-32 in CTab spec, 1-999 in RDKit)
/// iii: Number of dependent Rgroup (IF rrr THEN iii)
/// hhh: REstH property of rrr (0 (default)=off, 1=on): RGroup or H atom
/// ooo...: Range of RGroup occurrence required: n=exactly n, n-m=from n through m (inclusive),
/// >n=greater n, <n=fewer than n, blank (default): > 0.
/// Any non-contradictory combination of the preceding values is allowed, separated by commas.
fn rgroup_logic_entry<'a>(
) -> impl Parser<&'a [u8], Output = RGroupLogicEntry, Error = error::Error<&'a [u8]>> {
    |input: &'a [u8]| {
        let (remaining, _) = fixed_width_int_in_range::<u8, _>(3, 1..=1).parse(input)?;

        let (remaining, label) = map_parser(
            take(4usize),
            preceded(tag(" "), fixed_width_int_in_range::<u32, _>(3, 1..=999)),
        )
        .parse(remaining)?;

        let (remaining, dependent_label) =
            map_parser(take(4usize), preceded(tag(" "), fixed_width_int::<u32>(3)))
                .parse(remaining)?;

        let (remaining, rgroup_or_h) = if remaining.len() >= 4 {
            let (rest, value) = map_parser(
                take(4usize),
                preceded(tag(" "), fixed_width_int_in_range::<u8, _>(3, 0..=1)),
            )
            .parse(remaining)?;
            (rest, value)
        } else {
            (remaining, 0) // default value
        };

        let (remaining, occurrence) = if !remaining.is_empty() {
            let (rest, occurrence_bytes) = preceded(tag(" "), not_line_ending).parse(remaining)?;
            let (_, occurrence) = rgroup_occurrences().parse(occurrence_bytes)?;
            (rest, occurrence)
        } else {
            (remaining, vec![RGroupOccurrence::GreaterThan(0)])
        };

        Ok((
            remaining,
            RGroupLogicEntry {
                label,
                dependent_label: if dependent_label == 0 {
                    None
                } else {
                    Some(dependent_label)
                },
                rgroup_or_h: rgroup_or_h != 0,
                occurrence,
            },
        ))
    }
}

/// Parse SGroup type entries.
/// M  STYnn8 sss ttt ...
/// sss: SGroup index, ttt: SGroup type (3-character string)
fn sgroup_type_entries<'a>(
) -> impl Parser<&'a [u8], Output = Vec<SGroupTypeEntry>, Error = error::Error<&'a [u8]>> {
    all_consuming(length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                ),
                map_parser(
                    take(4usize),
                    preceded(
                        tag(" "),
                        map_res(take(3usize), |bytes: &[u8]| {
                            let s = std::str::from_utf8(bytes)
                                .map_err(|_| error::Error::new(bytes, error::ErrorKind::MapRes))?;
                            SGroup::get_type(s)
                                .map_err(|_| error::Error::new(bytes, error::ErrorKind::MapRes))
                        }),
                    ),
                ),
            ),
            |(sgroup_index, sgroup_type)| SGroupTypeEntry {
                sgroup_index,
                sgroup_type,
            },
        ),
    ))
}

/// Parse SGroup subtype entries.
/// M  SSTnn8 sss ttt ...
/// sss: SGroup index, ttt: SGroup subtype (3-character string)
fn sgroup_subtype_entries<'a>(
) -> impl Parser<&'a [u8], Output = Vec<SGroupSubtypeEntry>, Error = error::Error<&'a [u8]>> {
    all_consuming(length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                ),
                map_parser(
                    take(4usize),
                    preceded(
                        tag(" "),
                        map_res(take(3usize), |bytes: &[u8]| {
                            let s = std::str::from_utf8(bytes)
                                .map_err(|_| error::Error::new(bytes, error::ErrorKind::MapRes))?;
                            SGroup::get_subtype(s)
                                .map_err(|_| error::Error::new(bytes, error::ErrorKind::MapRes))
                        }),
                    ),
                ),
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
fn sgroup_label_entries<'a>(
) -> impl Parser<&'a [u8], Output = Vec<SGroupLabelEntry>, Error = error::Error<&'a [u8]>> {
    all_consuming(length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                ),
                map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_in_range::<u32, _>(3, 1..=512)),
                ),
            ),
            |(sgroup_index, label)| SGroupLabelEntry {
                sgroup_index,
                label,
            },
        ),
    ))
}

/// Parse SGroup connectivity entries.
/// M  SCNnn8 sss ttt ...
/// sss: SGroup index, ttt: SGroup connectivity (2-character string), left-justified
fn sgroup_connectivity_entries<'a>(
) -> impl Parser<&'a [u8], Output = Vec<SGroupConnectivityEntry>, Error = error::Error<&'a [u8]>> {
    all_consuming(length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                ),
                map_parser(
                    take(4usize),
                    preceded(
                        tag(" "),
                        map_res(take(2usize), |bytes: &[u8]| {
                            let s = std::str::from_utf8(bytes)
                                .map_err(|_| error::Error::new(bytes, error::ErrorKind::MapRes))?;
                            SGroup::get_connectivity(s)
                                .map_err(|_| error::Error::new(bytes, error::ErrorKind::MapRes))
                        }),
                    ),
                ),
            ),
            |(sgroup_index, connectivity)| SGroupConnectivityEntry {
                sgroup_index,
                connectivity,
            },
        ),
    ))
}

/// Parse SGroup expansion entries.
/// M SDS EXPn15 sss ...
/// sss: SGroup index, n15: count (max 15)
fn sgroup_expansion_entries<'a>(
) -> impl Parser<&'a [u8], Output = Vec<SGroupExpansionEntry>, Error = error::Error<&'a [u8]>> {
    all_consuming(preceded(
        tag(" EXP"),
        length_count(
            fixed_width_int_in_range::<u8, _>(3, 1..=15),
            map(
                map_parser(
                    take(4usize),
                    preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                ),
                |sgroup_index| SGroupExpansionEntry { sgroup_index },
            ),
        ),
    ))
}

/// Parse SGroup atom list entry.
/// M  SAL sssn15 aaa ...
/// sss: SGroup index, n15: count (max 15), aaa: atom indices (each 4 chars: " aaa")
fn sgroup_atom_list_entry<'a>(
) -> impl Parser<&'a [u8], Output = SGroupAtomListEntry, Error = error::Error<&'a [u8]>> {
    all_consuming(preceded(
        tag(" "),
        map(
            (
                fixed_width_int_minus1::<usize>(3),
                length_count(
                    fixed_width_int_in_range::<u8, _>(3, 1..=15),
                    map_parser(
                        take(4usize),
                        preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                    ),
                ),
            ),
            |(sgroup_index, atom_indices)| SGroupAtomListEntry {
                sgroup_index,
                atom_indices,
            },
        ),
    ))
}

/// Parse SGroup bond list entry.
/// M  SBL sssn15 bbb ...
/// sss: SGroup index (3 chars), n: count (3 chars), bbb: bond indices (each 4 chars: " bbb")
fn sgroup_bond_list_entry<'a>(
) -> impl Parser<&'a [u8], Output = SGroupBondListEntry, Error = error::Error<&'a [u8]>> {
    all_consuming(preceded(
        tag(" "),
        map(
            (
                fixed_width_int_minus1::<usize>(3),
                length_count(
                    fixed_width_int_in_range::<u8, _>(3, 1..=15),
                    map_parser(
                        take(4usize),
                        preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                    ),
                ),
            ),
            |(sgroup_index, bond_indices)| SGroupBondListEntry {
                sgroup_index,
                bond_indices,
            },
        ),
    ))
}

/// Parse SGroup parent atom entries.
/// M  SPA sssn15 aaa ...
/// sss: SGroup index, n15: count (max 15), aaa: atom indices (each 4 chars: " aaa")
fn sgroup_parent_atom_entries<'a>(
) -> impl Parser<&'a [u8], Output = SGroupParentAtomEntry, Error = error::Error<&'a [u8]>> {
    all_consuming(preceded(
        tag(" "),
        map(
            (
                fixed_width_int_minus1::<usize>(3),
                length_count(
                    fixed_width_int_in_range::<u8, _>(3, 1..=15),
                    map_parser(
                        take(4usize),
                        preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                    ),
                ),
            ),
            |(sgroup_index, atom_indices)| SGroupParentAtomEntry {
                sgroup_index,
                atom_indices,
            },
        ),
    ))
}

/// Parse SGroup subscript entry.
/// M  SMT sss m...
/// sss: SGroup index, m: subscript text
/// For multiple groups, m... is the text representation of the multiple group multiplier.For superatoms,
/// m... is the text of the superatom label.)
fn sgroup_subscript_entry<'a>(
) -> impl Parser<&'a [u8], Output = SGroupSubscriptEntry, Error = error::Error<&'a [u8]>> {
    all_consuming(map(
        (
            preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
            preceded(tag(" "), not_line_ending),
        ),
        |(sgroup_index, subscript)| SGroupSubscriptEntry {
            sgroup_index,
            subscript: String::from_utf8_lossy(subscript).trim().to_string(),
        },
    ))
}

/// Parse SGroup correspondence entry.
/// M  CRS sssnn6 bb1 bb2 bb3
/// sss: SGroup index, nn6: count (max 6), bb1-bb3: bond indices
/// bb1, bb2: Crossing bonds that share a common bracket
/// bb3: Crossing bond in repeating unit that connects to bond bb1
fn sgroup_correspondence_entry<'a>(
) -> impl Parser<&'a [u8], Output = SGroupCorrespondenceEntry, Error = error::Error<&'a [u8]>> {
    all_consuming(preceded(
        tag(" "),
        map(
            (
                fixed_width_int_minus1::<usize>(3),
                length_count(
                    fixed_width_int_in_range::<usize, _>(3, 1..=6),
                    map_parser(
                        take(4usize),
                        preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                    ),
                ),
            ),
            |(sgroup_index, bond_indices)| SGroupCorrespondenceEntry {
                sgroup_index,
                bond_indices,
            },
        ),
    ))
}

/// Parse SGroup display info entry.
/// M SDI sssnn4 x1 y1 x2 y2
/// sss: SGroup index, nn4: count (max 4), x1, y1: coordinates of opening bracket, x2, y2: coordinates of closing bracket
/// (f10.4)
fn sgroup_display_info_entry<'a>(
) -> impl Parser<&'a [u8], Output = SGroupDisplayInfoEntry, Error = error::Error<&'a [u8]>> {
    all_consuming(preceded(
        tag(" "),
        map(
            (
                fixed_width_int_minus1::<usize>(3),
                length_count(
                    fixed_width_int_in_range::<usize, _>(3, 1..=4),
                    fixed_width_float::<f64>(10, 4),
                ),
            ),
            |(sgroup_index, bracket_coords)| SGroupDisplayInfoEntry {
                sgroup_index,
                bracket_coords,
            },
        ),
    ))
}

#[cfg(test)]
mod tests;
