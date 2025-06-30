//! Properties block parser for CTab files.

use nom::bytes::complete::{tag, take};
use nom::character::complete::space0;
use nom::combinator::{map, map_parser, rest};
use nom::error;
use nom::multi::length_count;
use nom::sequence::preceded;
use nom::Parser;

use super::utils::{fixed_width_int, fixed_width_int_in_range, fixed_width_int_minus1};
use umol_data::Element;

#[derive(Debug, Clone, PartialEq)]
pub struct ChargeEntry {
    pub atom_index: usize,
    pub charge: i8,
}

#[derive(Debug, Clone, PartialEq)]
pub struct RadicalEntry {
    pub atom_index: usize,
    pub radical_type: i8,
}

#[derive(Debug, Clone, PartialEq)]
pub struct IsotopeEntry {
    pub atom_index: usize,
    pub mass: u32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupTypeEntry {
    pub sgroup_index: usize,
    pub sgroup_type: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SGroupLabelEntry {
    pub sgroup_index: usize,
    pub label: String,
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
    pub attachments: Vec<(usize, u8)>, // (neighbor_index, order)
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
    pub repeat_count: u8, // vvv >= 2
    pub bond1: usize,     // bbb (can be 0)
    pub bond2: usize,     // ccc (can be 0)
}

/// An enum representing a parsed property modification, containing the raw data.
/// This avoids allocating a new Vec for every single property line in a file.
#[derive(Debug, Clone, PartialEq)]
pub enum PropertyEntries {
    ChargeEntries(Vec<ChargeEntry>),
    RadicalEntries(Vec<RadicalEntry>),
    IsotopeEntries(Vec<IsotopeEntry>),
    SGroupTypeEntries(Vec<SGroupTypeEntry>),
    SGroupLabelEntries(Vec<SGroupLabelEntry>),
    SGroupAtomListEntry(SGroupAtomListEntry),
    SGroupBondListEntry(SGroupBondListEntry),
    AtomAliasEntry(AtomAliasEntry),
    AtomValueEntry(AtomValueEntry),
    AtomListEntry(AtomListEntry),
    AttachmentPointEntries(Vec<AttachmentPointEntry>),
    AtomAttachmentOrderEntry(AtomAttachmentOrderEntry),
    RingBondCountEntries(Vec<RingBondCountEntry>),
    SubstitutionCountEntries(Vec<SubstitutionCountEntry>),
    UnsaturatedAtomEntries(Vec<UnsaturatedAtomEntry>),
    LinkAtomEntries(Vec<LinkAtomEntry>),
}

/// Parse charge property entries.
/// nn8 aaa vvv ...
/// vvv: -15..= 15.
fn charge_entries<'a>(
) -> impl Parser<&'a [u8], Output = Vec<ChargeEntry>, Error = error::Error<&'a [u8]>> {
    length_count(
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
    )
}

/// Parse radical property entries.
/// nn8 aaa vvv ...
/// vvv: 0..= 3: 0 = no radical, 1 = singlet (:), 2 = doublet (. or ^), 3 = triplet (^^).
fn radical_entries<'a>(
) -> impl Parser<&'a [u8], Output = Vec<RadicalEntry>, Error = error::Error<&'a [u8]>> {
    length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                preceded(tag(" "), fixed_width_int_in_range::<i8, _>(3, 0..=3)),
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
fn isotope_entries<'a>(
) -> impl Parser<&'a [u8], Output = Vec<IsotopeEntry>, Error = error::Error<&'a [u8]>> {
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

/// Parse SGroup type entries.
/// M  STYnn8 sss ttt ...
/// sss: SGroup index, ttt: SGroup type (3-character string)
fn sgroup_type_entries<'a>(
) -> impl Parser<&'a [u8], Output = Vec<SGroupTypeEntry>, Error = error::Error<&'a [u8]>> {
    length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                map_parser(take(4usize), preceded(tag(" "), take(3usize))),
            ),
            |(sgroup_index, type_bytes)| SGroupTypeEntry {
                sgroup_index,
                sgroup_type: String::from_utf8_lossy(type_bytes).trim().to_string(),
            },
        ),
    )
}

/// Parse SGroup label entries.
/// M  SLBnn8 sss vvv ...
/// sss: SGroup index, vvv: label (3-character string, can be longer in practice)
fn sgroup_label_entries<'a>(
) -> impl Parser<&'a [u8], Output = Vec<SGroupLabelEntry>, Error = error::Error<&'a [u8]>> {
    length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                map_parser(take(4usize), preceded(tag(" "), take(3usize))),
            ),
            |(sgroup_index, label_bytes)| SGroupLabelEntry {
                sgroup_index,
                label: String::from_utf8_lossy(label_bytes).trim().to_string(),
            },
        ),
    )
}

/// Parse SGroup atom list entry.
/// M  SAL sssn15 aaa ...
/// sss: SGroup index (3 chars), n: count (3 chars), aaa: atom indices (each 4 chars: " aaa")
fn sgroup_atom_list_entry<'a>(
) -> impl Parser<&'a [u8], Output = SGroupAtomListEntry, Error = error::Error<&'a [u8]>> {
    |input: &'a [u8]| {
        // Parse sgroup index and count from first 6 characters
        let (remaining, (sgroup_index, count)) = map_parser(
            take(6usize),
            (
                fixed_width_int_minus1::<usize>(3),
                fixed_width_int_in_range::<u8, _>(3, 1..=15),
            ),
        )
        .parse(input)?;

        // Parse the atom indices - each is 4 chars with format " aaa"
        let mut indices = Vec::with_capacity(count as usize);
        let mut remaining = remaining;
        for _ in 0..count {
            let (rest, index) = map_parser(
                take(4usize),
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
            )
            .parse(remaining)?;
            indices.push(index);
            remaining = rest;
        }

        Ok((
            remaining,
            SGroupAtomListEntry {
                sgroup_index,
                atom_indices: indices,
            },
        ))
    }
}

/// Parse SGroup bond list entry.
/// M  SBL sssn15 bbb ...
/// sss: SGroup index (3 chars), n: count (3 chars), bbb: bond indices (each 4 chars: " bbb")
fn sgroup_bond_list_entry<'a>(
) -> impl Parser<&'a [u8], Output = SGroupBondListEntry, Error = error::Error<&'a [u8]>> {
    |input: &'a [u8]| {
        // Parse sgroup index and count from first 6 characters
        let (remaining, (sgroup_index, count)) = map_parser(
            take(6usize),
            (
                fixed_width_int_minus1::<usize>(3),
                fixed_width_int_in_range::<u8, _>(3, 1..=15),
            ),
        )
        .parse(input)?;

        // Parse the bond indices - each is 4 chars with format " bbb"
        let mut indices = Vec::with_capacity(count as usize);
        let mut remaining = remaining;
        for _ in 0..count {
            let (rest, index) = map_parser(
                take(4usize),
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
            )
            .parse(remaining)?;
            indices.push(index);
            remaining = rest;
        }

        Ok((
            remaining,
            SGroupBondListEntry {
                sgroup_index,
                bond_indices: indices,
            },
        ))
    }
}

/// Parse atom alias entry.
/// A  aaa alias_text
/// aaa: atom index, alias_text: alias string (can contain spaces)
fn atom_alias_entry<'a>(
) -> impl Parser<&'a [u8], Output = AtomAliasEntry, Error = error::Error<&'a [u8]>> {
    map(
        (
            preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
            preceded(space0, rest),
        ),
        |(atom_index, alias_bytes)| AtomAliasEntry {
            atom_index,
            alias: String::from_utf8_lossy(alias_bytes).trim().to_string(),
        },
    )
}

/// Parse atom value entry.
/// V  aaa value_text
/// aaa: atom index, value_text: value string (can contain spaces)
fn atom_value_entry<'a>(
) -> impl Parser<&'a [u8], Output = AtomValueEntry, Error = error::Error<&'a [u8]>> {
    map(
        (
            preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
            preceded(space0, rest),
        ),
        |(atom_index, value_bytes)| AtomValueEntry {
            atom_index,
            value: String::from_utf8_lossy(value_bytes).trim().to_string(),
        },
    )
}

/// Parse atom list property entry.
/// M  ALS aaannn e 11112222333344445555...
/// aaa: Atom number, nnn: Number of entries (16 max), e: Exclusion (T/F), 1111: 4-char symbols
fn atom_list_entry<'a>(
) -> impl Parser<&'a [u8], Output = AtomListEntry, Error = error::Error<&'a [u8]>> {
    |input: &'a [u8]| {
        // Parse atom index (3 chars)
        let (remaining, atom_index) = fixed_width_int_minus1::<usize>(3).parse(input)?;

        // Parse count (3 chars, max 16)
        let (remaining, count) = fixed_width_int_in_range::<u8, _>(3, 1..=16).parse(remaining)?;

        // Parse exclusion flag (1 char)
        let (remaining, exclusion_byte) = take(1usize).parse(remaining)?;
        let exclusion = match exclusion_byte {
            b"T" => true,
            b"F" => false,
            _ => {
                return Err(nom::Err::Error(error::Error::new(
                    input,
                    error::ErrorKind::Tag,
                )))
            }
        };

        // Parse 4-char atom symbols
        let mut elements = Vec::with_capacity(count as usize);
        let mut remaining = remaining;
        for _ in 0..count {
            let (rest, symbol_bytes) = take(4usize).parse(remaining)?;
            let symbol_cow = String::from_utf8_lossy(symbol_bytes);
            let symbol_str = symbol_cow.trim();

            // Convert to Element
            let element = Element::from_symbol(symbol_str).ok_or_else(|| {
                nom::Err::Error(error::Error::new(input, error::ErrorKind::MapRes))
            })?;
            elements.push(element);
            remaining = rest;
        }

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
/// M  APOnn2 aaa vvv ...
/// nn2: Count (max 2), aaa: Atom index, vvv: Attachment type (0-3)
fn attachment_point_entries<'a>(
) -> impl Parser<&'a [u8], Output = Vec<AttachmentPointEntry>, Error = error::Error<&'a [u8]>> {
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
/// M  AAL aaan2 111 v1v 222 v2v ...
/// aaa: Atom index, n2: Pair count (max 2), 111/222: Neighbor indices, v1v/v2v: Orders
fn atom_attachment_order_entry<'a>(
) -> impl Parser<&'a [u8], Output = AtomAttachmentOrderEntry, Error = error::Error<&'a [u8]>> {
    |input: &'a [u8]| {
        // Parse atom index (3 chars)
        let (remaining, atom_index) = fixed_width_int_minus1::<usize>(3).parse(input)?;

        // Parse pair count (2 chars, max 2)
        let (remaining, count) = fixed_width_int_in_range::<u8, _>(2, 1..=2).parse(remaining)?;

        // Parse neighbor-order pairs
        let mut attachments = Vec::with_capacity(count as usize);
        let mut remaining = remaining;
        for _ in 0..count {
            let (rest, neighbor_index) =
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)).parse(remaining)?;
            let (rest, order) =
                preceded(tag(" "), fixed_width_int_in_range::<u8, _>(3, 1..=2)).parse(rest)?;
            attachments.push((neighbor_index, order));
            remaining = rest;
        }

        Ok((
            remaining,
            AtomAttachmentOrderEntry {
                atom_index,
                attachments,
            },
        ))
    }
}

/// Parse ring bond count property entries.
/// M  RBCnn8 aaa vvv ...
/// vvv: Ring bond count (-2=r*, -1=r0, 0=off, 2=r2, 3=r3, 4+=r4)
fn ring_bond_count_entries<'a>(
) -> impl Parser<&'a [u8], Output = Vec<RingBondCountEntry>, Error = error::Error<&'a [u8]>> {
    length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=8),
        map(
            (
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                preceded(tag(" "), fixed_width_int_in_range::<i8, _>(3, -2..=15)),
            ),
            |(atom_index, ring_bond_count)| RingBondCountEntry {
                atom_index,
                ring_bond_count,
            },
        ),
    )
}

/// Parse substitution count property entries.
/// M  SUBnn8 aaa vvv ...
/// vvv: Substitution count (-2=s*, -1=s0, 0=off, 1-5=s1-s5, 6+=s6)
fn substitution_count_entries<'a>(
) -> impl Parser<&'a [u8], Output = Vec<SubstitutionCountEntry>, Error = error::Error<&'a [u8]>> {
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
/// vvv: Unsaturated flag (0=off, 1=on)
fn unsaturated_atom_entries<'a>(
) -> impl Parser<&'a [u8], Output = Vec<UnsaturatedAtomEntry>, Error = error::Error<&'a [u8]>> {
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
/// nn4: Count (max 4), aaa: Atom index, vvv: Repeat count (>=2), bbb/ccc: Bond indices (can be 0)
fn link_atom_entries<'a>(
) -> impl Parser<&'a [u8], Output = Vec<LinkAtomEntry>, Error = error::Error<&'a [u8]>> {
    length_count(
        fixed_width_int_in_range::<u8, _>(3, 1..=4),
        map(
            (
                preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                preceded(tag(" "), fixed_width_int_in_range::<u8, _>(3, 2..=255)),
                preceded(tag(" "), fixed_width_int::<usize>(3)),
                preceded(tag(" "), fixed_width_int::<usize>(3)),
            ),
            |(atom_index, repeat_count, bond1, bond2)| LinkAtomEntry {
                atom_index,
                repeat_count,
                bond1,
                bond2,
            },
        ),
    )
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
                b"M  SLB" => sgroup_label_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SGroupLabelEntries(o))),
                b"M  SAL" => sgroup_atom_list_entry()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SGroupAtomListEntry(o))),
                b"M  SBL" => sgroup_bond_list_entry()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SGroupBondListEntry(o))),
                // Query properties are not supported in standard parser
                b"M  ALS" | b"M  APO" | b"M  AAL" | b"M  RBC" | b"M  SUB" | b"M  UNS"
                | b"M  LIN" => Err(nom::Err::Error(error::Error::new(
                    input,
                    error::ErrorKind::Tag,
                ))),
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
                b"M  STY" => sgroup_type_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SGroupTypeEntries(o))),
                b"M  SLB" => sgroup_label_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SGroupLabelEntries(o))),
                b"M  SAL" => sgroup_atom_list_entry()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SGroupAtomListEntry(o))),
                b"M  SBL" => sgroup_bond_list_entry()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::SGroupBondListEntry(o))),
                b"M  ALS" => atom_list_entry()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::AtomListEntry(o))),
                b"M  APO" => attachment_point_entries()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::AttachmentPointEntries(o))),
                b"M  AAL" => atom_attachment_order_entry()
                    .parse(rest)
                    .map(|(i, o)| (i, PropertyEntries::AtomAttachmentOrderEntry(o))),
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
                _ => Err(nom::Err::Error(error::Error::new(
                    input,
                    error::ErrorKind::Tag,
                ))),
            }
        }
    }
}

#[cfg(test)]
mod tests;
