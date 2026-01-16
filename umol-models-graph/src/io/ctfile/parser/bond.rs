//! Bond block parsers for CTab files.

use nom::character::complete::space0;
use nom::combinator::{all_consuming, cond, map_res};
use nom::error::{Error as NomError, ErrorKind as NomErrorKind};
use nom::sequence::terminated;
use nom::{Err, Parser};

use super::convert::{
    convert_bond_reacting_center_code, convert_bond_stereo_direction_code,
    convert_bond_topology_code, convert_bond_type_code, convert_extended_bond_type_code,
};
use super::utils::{
    fixed_width_int, fixed_width_int_minus1, fixed_width_int_opt_unsigned, fixed_width_unused_n,
    LinesWithOffsetExt,
};
use crate::io::ctfile::config::CtabParseFlags;
use crate::io::ctfile::error::ParseError;
use crate::table_ir::bond::{Bond, BondOrder, ExtendedBond};

/// Parse bond block (basic bonds only)
pub(super) fn bond_block<'inp>(
    bond_count: u32,
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (Vec<(usize, usize, Bond)>, u32), Error = ParseError> + use<'inp>
{
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

            let (_, (atom1, atom2, bond)) = all_consuming(terminated(bond_input(flags), space0))
                .parse(line)
                .map_err(|e| {
                    Err::Error(ParseError::bond_from_nom(e, line_offset + line_index, line))
                })?;
            bonds.push((atom1, atom2, bond));
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
) -> impl Parser<&'inp [u8], Output = (Vec<(usize, usize, ExtendedBond)>, u32), Error = ParseError>
       + use<'inp> {
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

            let (_, (atom1, atom2, bond)) =
                all_consuming(terminated(extended_bond_input(flags), space0))
                    .parse(line)
                    .map_err(|e| {
                        Err::Error(ParseError::bond_from_nom(e, line_offset + line_index, line))
                    })?;
            bonds.push((atom1, atom2, bond));
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
) -> impl Parser<&'inp [u8], Output = (usize, usize, Bond), Error = NomError<&'inp [u8]>> + use<'inp>
{
    move |input: &'inp [u8]| {
        if input.len() < 9 {
            return Err(Err::Error(NomError::new(input, NomErrorKind::Eof)));
        }

        if flags == CtabParseFlags::BASIC {
            match input.len() {
                21 => return basic_bond_input21(input),
                18 => return basic_bond_input18(input),
                _ => {}
            }
        }

        let extended_range = flags.contains(CtabParseFlags::EXTENDED_RANGE);
        let skip_unused_fields = flags.contains(CtabParseFlags::SKIP_UNUSED_FIELDS);

        // Atom indices
        let (remaining, first_atom) = fixed_width_int_minus1::<usize>(3).parse(input)?;
        let (remaining, second_atom) = fixed_width_int_minus1::<usize>(3).parse(remaining)?;

        // Bond type
        let (remaining, order) = map_res(fixed_width_int::<u8>(3), move |code| {
            convert_bond_type_code(code, extended_range)
        })
        .parse(remaining)?;

        // Stereo/direction
        let (remaining, stereo_direction) = cond(
            remaining.len() >= 3,
            map_res(fixed_width_int::<u8>(3), convert_bond_stereo_direction_code),
        )
        .parse(remaining)?;

        // Ignored field
        let (remaining, _) = cond(
            !remaining.is_empty(),
            fixed_width_unused_n((remaining.len() / 3).min(1), 3, skip_unused_fields),
        )
        .parse(remaining)?;

        // Unused fields
        let (remaining, _) = cond(
            !remaining.is_empty(),
            fixed_width_unused_n((remaining.len() / 3).min(2), 3, false),
        )
        .parse(remaining)?;

        let mut bond = Bond::with_order(order);
        if let Some((stereo_val, direction)) = stereo_direction {
            match order {
                BondOrder::Single => {
                    bond.direction = direction;
                }
                BondOrder::Double => {
                    bond.stereo = stereo_val;
                }
                _ => (),
            }
        }

        Ok((remaining, (first_atom, second_atom, bond)))
    }
}

/// Fast-path basic bond input parser for 21-character lines.  Hard-codes
/// CtabParseFlags::BASIC = NAMED_ISOTOPES | ATOM_MAP_HCOUNT_FIELDS | SKIP_UNUSED_FIELDS behavior.
/// behavior
fn basic_bond_input21<'inp>(
    input: &'inp [u8],
) -> Result<(&'inp [u8], (usize, usize, Bond)), Err<NomError<&'inp [u8]>>> {
    let line = &input[..21];
    let remaining = &input[21..];

    let first_atom = fixed_width_int_opt_unsigned(input, &line[0..3])?
        .and_then(|val| (val >= 1).then_some((val - 1) as usize))
        .ok_or_else(|| Err::Error(NomError::new(input, NomErrorKind::Verify)))?;
    let second_atom = fixed_width_int_opt_unsigned(input, &line[3..6])?
        .and_then(|val| (val >= 1).then_some((val - 1) as usize))
        .ok_or_else(|| Err::Error(NomError::new(input, NomErrorKind::Verify)))?;
    let order_code = fixed_width_int_opt_unsigned(input, &line[6..9])?.unwrap_or(0) as u8;
    let order = convert_bond_type_code(order_code, false)
        .map_err(|_| Err::Error(NomError::new(input, NomErrorKind::MapRes)))?;

    let stereo_code = fixed_width_int_opt_unsigned(input, &line[9..12])?.unwrap_or(0) as u8;
    let stereo_direction = convert_bond_stereo_direction_code(stereo_code)
        .map_err(|_| Err::Error(NomError::new(input, NomErrorKind::Verify)))?;

    if !super::utils::is_all_whitespace_or_zeroes(&line[15..18]) {
        return Err(Err::Error(NomError::new(input, NomErrorKind::Verify)));
    }
    if !super::utils::is_all_whitespace_or_zeroes(&line[18..21]) {
        return Err(Err::Error(NomError::new(input, NomErrorKind::Verify)));
    }

    let mut bond = Bond::with_order(order);
    if let (Some(stereo), Some(direction)) = stereo_direction {
        match order {
            BondOrder::Single => bond.direction = Some(direction),
            BondOrder::Double => bond.stereo = Some(stereo),
            _ => {}
        }
    }
    Ok((remaining, (first_atom, second_atom, bond)))
}

/// Fast-path basic bond input parser for 18-character lines.  Hard-codes
/// CtabParseFlags::BASIC = NAMED_ISOTOPES | ATOM_MAP_HCOUNT_FIELDS | SKIP_UNUSED_FIELDS behavior.
///
fn basic_bond_input18<'inp>(
    input: &'inp [u8],
) -> Result<(&'inp [u8], (usize, usize, Bond)), Err<NomError<&'inp [u8]>>> {
    let line = &input[..18];
    let remaining = &input[18..];

    let first_atom = fixed_width_int_opt_unsigned(input, &line[0..3])?
        .and_then(|val| (val >= 1).then_some((val - 1) as usize))
        .ok_or_else(|| Err::Error(NomError::new(input, NomErrorKind::Verify)))?;
    let second_atom = fixed_width_int_opt_unsigned(input, &line[3..6])?
        .and_then(|val| (val >= 1).then_some((val - 1) as usize))
        .ok_or_else(|| Err::Error(NomError::new(input, NomErrorKind::Verify)))?;
    let order_code = fixed_width_int_opt_unsigned(input, &line[6..9])?.unwrap_or(0) as u8;
    let order = convert_bond_type_code(order_code, false)
        .map_err(|_| Err::Error(NomError::new(input, NomErrorKind::MapRes)))?;

    let stereo_code = fixed_width_int_opt_unsigned(input, &line[9..12])?.unwrap_or(0) as u8;
    let stereo_direction = convert_bond_stereo_direction_code(stereo_code)
        .map_err(|_| Err::Error(NomError::new(input, NomErrorKind::Verify)))?;

    if !super::utils::is_all_whitespace_or_zeroes(&line[15..18]) {
        return Err(Err::Error(NomError::new(input, NomErrorKind::Verify)));
    }

    let mut bond = Bond::with_order(order);
    if let (Some(stereo), Some(direction)) = stereo_direction {
        match order {
            BondOrder::Single => bond.direction = Some(direction),
            BondOrder::Double => bond.stereo = Some(stereo),
            _ => {}
        }
    }
    Ok((remaining, (first_atom, second_atom, bond)))
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
) -> impl Parser<&'inp [u8], Output = (usize, usize, ExtendedBond), Error = NomError<&'inp [u8]>>
       + use<'inp> {
    let skip_unused_fields = flags.contains(CtabParseFlags::SKIP_UNUSED_FIELDS);
    let extended_range = flags.contains(CtabParseFlags::EXTENDED_RANGE);
    let allow_wildcards = flags.contains(CtabParseFlags::WILDCARDS);
    move |input: &'inp [u8]| {
        if input.len() < 9 {
            return Err(Err::Error(NomError::new(input, NomErrorKind::Eof)));
        }

        // Atom indices
        let (i, first_atom) = fixed_width_int_minus1::<usize>(3).parse(input)?;
        let (i, second_atom) = fixed_width_int_minus1::<usize>(3).parse(i)?;

        let (i, order) = map_res(fixed_width_int::<u8>(3), move |code| {
            convert_extended_bond_type_code(code, extended_range, allow_wildcards)
        })
        .parse(i)?;

        // Stereo/dir
        let (i, stereo_direction) = cond(
            i.len() >= 3,
            map_res(fixed_width_int::<u8>(3), convert_bond_stereo_direction_code),
        )
        .parse(i)?;

        // Ignore xxx field
        let (i, _) = cond(
            !i.is_empty(),
            fixed_width_unused_n((i.len() / 3).min(1), 3, skip_unused_fields),
        )
        .parse(i)?;

        // Topology, reacting center
        let (i, topology) = cond(
            i.len() >= 3,
            map_res(fixed_width_int::<u8>(3), convert_bond_topology_code),
        )
        .parse(i)?;
        let (i, reacting_center) = cond(
            i.len() >= 3,
            map_res(fixed_width_int::<i8>(3), move |code| {
                convert_bond_reacting_center_code(code, extended_range)
            }),
        )
        .parse(i)?;

        let mut bond = ExtendedBond::with_order(order);
        if let Some((stereo_val, direction)) = stereo_direction {
            match order {
                BondOrder::Single => {
                    bond.direction = direction;
                }
                BondOrder::Double => {
                    bond.stereo = stereo_val;
                }
                _ => (),
            }
        }
        bond.topology = topology.flatten();
        bond.reacting_center = reacting_center.flatten();

        Ok((i, (first_atom, second_atom, bond)))
    }
}

#[cfg(test)]
mod tests;
