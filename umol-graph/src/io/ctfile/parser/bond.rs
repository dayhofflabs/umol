//! Bond block parsers for CTab files.

use nom::character::complete::space0;
use nom::combinator::all_consuming;
use nom::error::{Error as NomError, ErrorKind as NomErrorKind};
use nom::sequence::terminated;
use nom::{Err, Parser};

use super::convert::{
    convert_bond_reacting_center_code, convert_bond_stereo_direction_code,
    convert_bond_topology_code, convert_bond_type_code, convert_extended_bond_type_code,
};
use super::utils::{parse_int_opt, validate_unused_n, LinesWithOffsetExt};
use crate::io::ctfile::config::CtabParseFlags;
use crate::io::ctfile::error::ParseError;
use crate::table_ir::bond::{Bond, BondOrder, ExtendedBond};

/// Parse bond block (basic bonds only)
pub(super) fn bond_block<'inp>(
    bond_count: u32,
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (Vec<Bond>, u32), Error = ParseError> + use<'inp> {
    move |input: &'inp [u8]| {
        let mut bonds = Vec::with_capacity(bond_count as usize);
        let mut lines_iter = input.lines_with_offset();
        let mut byte_offset = 0;

        for line_index in 0..bond_count {
            let (line, byte_len) = lines_iter.next().ok_or_else(|| {
                Err::Error(ParseError::UnexpectedEof {
                    line: line_offset + line_index,
                    block: "bond",
                })
            })?;

            let (_, bond) = all_consuming(terminated(bond_input(flags), space0))
                .parse(line)
                .map_err(|e| {
                    Err::Error(ParseError::bond_from_nom(e, line_offset + line_index, line))
                })?;
            bonds.push(bond);
            byte_offset += byte_len;
        }

        let remaining = &input[byte_offset..];
        Ok((remaining, (bonds, line_offset + bond_count)))
    }
}

/// Parse extended bond block
pub(super) fn extended_bond_block<'inp>(
    bond_count: u32,
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (Vec<ExtendedBond>, u32), Error = ParseError> + use<'inp> {
    move |input: &'inp [u8]| {
        let mut bonds = Vec::with_capacity(bond_count as usize);
        let mut lines_iter = input.lines_with_offset();
        let mut byte_offset = 0;

        for line_index in 0..bond_count {
            let (line, byte_len) = lines_iter.next().ok_or_else(|| {
                Err::Error(ParseError::UnexpectedEof {
                    line: line_offset + line_index,
                    block: "bond",
                })
            })?;

            let (_, bond) = all_consuming(terminated(extended_bond_input(flags), space0))
                .parse(line)
                .map_err(|e| {
                    Err::Error(ParseError::bond_from_nom(e, line_offset + line_index, line))
                })?;
            bonds.push(bond);
            byte_offset += byte_len;
        }

        let remaining = &input[byte_offset..];
        Ok((remaining, (bonds, line_offset + bond_count)))
    }
}

/// Parse bond input (optimized for performance)
/// Fails immediately on query bond properties. For parsing all bond types, see extended_bond_input.
///
/// "111222tttsssxxxrrrccc" (21 characters wide)
///
/// *Values in the bond block*
/// -------------------------------------------------------------------------
/// | Field | Position | Meaning              | Values     | Notes          |
/// |-------|----------|----------------------|------------|----------------|
/// | 111   | 1-3      | first atom           | 1..=aaa    | Generic        |
/// | 222   | 4-6      | second atom          | 1..=aaa    | Generic        |
/// | ttt   | 7-9      | bond type            | 1..=8      | Generic,Query  |
/// | sss   | 10-12    | bond stereo          | 0..=6      | Generic        |
/// | rrr   | 16-18    | bond topology        | 0..=2      | Query          |
/// | ccc   | 19-21    | bond reacting center | 0..=15     | Reaction,Query |
/// -------------------------------------------------------------------------
///
/// *Behavior in unused and extended fields*
/// ---------------------------------------------------------------------
/// | Field    | Basic      | Basic strict | Extended | Extended strict |
/// ---------------------------------------------------------------------
/// | Unused   | skip       | zero/blank   | skip     | zero/blank      |
/// | Extended | zero/blank | zero/blank   | validate | validate        |
/// ---------------------------------------------------------------------
/// skip: skip field, regardless of content
/// zero/blank: validate and accept only zero or blank values, reject any other content
/// validate: validate field according to its specification, reject any other content
///
/// NOTE: Basic parser should accept a strict subset of inputs accepted by extended parser
///       Increasing strictness: skip < validate < zero/blank
///
pub fn bond_input<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = Bond, Error = NomError<&'inp [u8]>> + use<'inp> {
    move |input: &'inp [u8]| {
        if input.len() < 9 {
            return Err(Err::Error(NomError::new(input, NomErrorKind::Eof)));
        }

        let mut offset;

        let extended_range = flags.contains(CtabParseFlags::EXTENDED_RANGE);
        let skip_unused_fields = flags.contains(CtabParseFlags::SKIP_UNUSED_FIELDS);

        // Atom indices (0-2, 3-5)
        let first_atom = parse_int_opt::<u32>(input, &input[0..3])?
            .and_then(|val| (val >= 1).then_some((val - 1) as usize))
            .ok_or_else(|| Err::Error(NomError::new(input, NomErrorKind::Verify)))?;
        let second_atom = parse_int_opt::<u32>(input, &input[3..6])?
            .and_then(|val| (val >= 1).then_some((val - 1) as usize))
            .ok_or_else(|| Err::Error(NomError::new(input, NomErrorKind::Verify)))?;

        // Bond type (6-8)
        let order_code = parse_int_opt::<u8>(input, &input[6..9])?.unwrap_or(0);
        let order = convert_bond_type_code(order_code, extended_range)
            .map_err(|_| Err::Error(NomError::new(input, NomErrorKind::MapRes)))?;

        offset = 9;

        // Stereo/direction (9-11)
        let (stereo, wedge) = if input.len() >= 12 {
            offset = 12;
            let stereo_code = parse_int_opt::<u8>(input, &input[9..12])?.unwrap_or(0);
            convert_bond_stereo_direction_code(stereo_code)
        } else {
            (None, None)
        };

        // Ignored field xxx (12-14)
        if input.len() >= 15 {
            offset = 15;
            validate_unused_n(input, &input[12..15], 1, 3, skip_unused_fields)?;
        }

        // Bond topology and reacting center (15-20) - extended
        let count = ((input.len().saturating_sub(15)) / 3).min(2);
        if count > 0 {
            offset = 15 + count * 3;
            validate_unused_n(input, &input[15..15 + count * 3], count, 3, false)?;
        }

        let mut bond = Bond::new(first_atom as u32, second_atom as u32, order);
        // The stereo/direction field applies only to single and double bonds.
        if matches!(order, BondOrder::Single | BondOrder::Double) {
            bond.stereo = stereo;
            bond.wedge = wedge;
        }

        Ok((&input[offset..], bond))
    }
}

/// Parse extended bond input
/// Allows all bond types. For faster parsing of basic bonds, see bond_input.
///
/// "111222tttsssxxxrrrccc" (21 characters wide)
///
/// *Values in the bond block*
/// -------------------------------------------------------------------------
/// | Field | Position | Meaning              | Values     | Notes          |
/// |-------|----------|----------------------|------------|----------------|
/// | 111   | 1-3      | first atom           | 1..=aaa    | Generic        |
/// | 222   | 4-6      | second atom          | 1..=aaa    | Generic        |
/// | ttt   | 7-9      | bond type            | 1..=8      | Generic, Query |
/// | sss   | 10-12    | bond stereo          | 0..=6      | Generic        |
/// | xxx   | 13-15    | ignored              |            |                |
/// | rrr   | 16-18    | bond topology        | 0..=2      | Query          |
/// | ccc   | 19-21    | bond reacting center | 0..=15     | Reaction,Query |
/// -------------------------------------------------------------------------
///
///
pub fn extended_bond_input<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = ExtendedBond, Error = NomError<&'inp [u8]>> + use<'inp> {
    let skip_unused_fields = flags.contains(CtabParseFlags::SKIP_UNUSED_FIELDS);
    let extended_range = flags.contains(CtabParseFlags::EXTENDED_RANGE);
    let allow_wildcards = flags.contains(CtabParseFlags::WILDCARDS);

    move |input: &'inp [u8]| {
        if input.len() < 9 {
            return Err(Err::Error(NomError::new(input, NomErrorKind::Eof)));
        }

        let mut offset = 9;

        // Atom indices (0-2, 3-5)
        let first_atom = parse_int_opt::<u32>(input, &input[0..3])?
            .and_then(|val| (val >= 1).then_some((val - 1) as usize))
            .ok_or_else(|| Err::Error(NomError::new(input, NomErrorKind::Verify)))?;
        let second_atom = parse_int_opt::<u32>(input, &input[3..6])?
            .and_then(|val| (val >= 1).then_some((val - 1) as usize))
            .ok_or_else(|| Err::Error(NomError::new(input, NomErrorKind::Verify)))?;

        // Bond type (6-8)
        let order_code = parse_int_opt::<u8>(input, &input[6..9])?.unwrap_or(0);
        let order = convert_extended_bond_type_code(order_code, extended_range, allow_wildcards)
            .map_err(|_| Err::Error(NomError::new(input, NomErrorKind::MapRes)))?;

        // Stereo/direction (9-11)
        let (stereo, wedge) = if input.len() >= 12 {
            offset = 12;
            let stereo_code = parse_int_opt::<u8>(input, &input[9..12])?.unwrap_or(0);
            convert_bond_stereo_direction_code(stereo_code)
        } else {
            (None, None)
        };

        // Ignored field xxx (12-14)
        if input.len() >= 15 {
            offset = 15;
            validate_unused_n(input, &input[12..15], 1, 3, skip_unused_fields)?;
        }

        // Topology rrr (15-17)
        let topology = if input.len() >= 18 {
            offset = 18;
            let val = parse_int_opt::<u8>(input, &input[15..18])?.unwrap_or(0);
            if val > 2 {
                return Err(Err::Error(NomError::new(input, NomErrorKind::Verify)));
            }
            convert_bond_topology_code(val)
        } else {
            None
        };

        // Reacting center ccc (18-20)
        let reacting_center = if input.len() >= 21 {
            offset = 21;
            let val = parse_int_opt::<i8>(input, &input[18..21])?.unwrap_or(0);
            convert_bond_reacting_center_code(val, extended_range)
                .map_err(|_| Err::Error(NomError::new(input, NomErrorKind::Verify)))?
        } else {
            None
        };

        // Validate any trailing bytes as whitespace
        if input.len() > offset && !input[offset..].trim_ascii().is_empty() {
            return Err(Err::Error(NomError::new(input, NomErrorKind::Verify)));
        }

        let mut bond = ExtendedBond::new(first_atom as u32, second_atom as u32, order);
        // The stereo/direction field applies only to single and double bonds.
        if matches!(order, BondOrder::Single | BondOrder::Double) {
            bond.stereo = stereo;
            bond.wedge = wedge;
        }
        bond.topology = topology;
        bond.reacting_center = reacting_center;

        Ok((&input[input.len()..], bond))
    }
}

#[cfg(test)]
mod tests;
