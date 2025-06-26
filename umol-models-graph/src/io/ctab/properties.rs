//! Properties block parser for CTab files.

use nom::{
    bytes::complete::{tag, take},
    combinator::{map, map_parser},
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

/// An enum representing a parsed property modification, containing the raw data.
/// This avoids allocating a new Vec for every single property line in a file.
#[derive(Debug, Clone, PartialEq)]
pub enum PropertyEntries {
    ChargeEntries(Vec<ChargeEntry>),
    RadicalEntries(Vec<RadicalEntry>),
    IsotopeEntries(Vec<IsotopeEntry>),
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

/// Parse property line
pub fn property_input_standard<'a>(
    input: &'a [u8],
) -> nom::IResult<&'a [u8], PropertyEntries, error::Error<&'a [u8]>> {
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
        _ => Err(nom::Err::Error(error::Error::new(
            input,
            error::ErrorKind::Tag,
        ))),
    }
}

#[cfg(test)]
mod tests;
