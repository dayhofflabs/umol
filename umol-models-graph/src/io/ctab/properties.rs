//! Properties block parser for CTab files.

use nom::{
    bytes::complete::{tag, take},
    character::complete::{space0},
    combinator::{map, map_parser, rest},
    error,
    multi::length_count,
    sequence::preceded,
    Parser,
};

use super::utils::{fixed_width_int, fixed_width_int_in_range, fixed_width_int_minus1};

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

/// Parse property line
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
