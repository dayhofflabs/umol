//! Bond block parsers for CTab files.

use winnow::error::ErrMode;
use winnow::token::take;
use winnow::{ModalResult, Parser};

use super::convert::{
    convert_bond_reacting_center_code, convert_bond_stereo_direction_code,
    convert_bond_topology_code, convert_bond_type_code, convert_extended_bond_type_code,
};
use super::utils::{
    finish_line, input_error_column, next_line, parse_int_opt, validate_unused_n, Input, InputError,
};
use crate::ctfile::config::CtabParseFlags;
use crate::ctfile::error::ParseError;
use crate::table_ir::bond::{Bond, BondOrder, ExtendedBond};

/// Parse bond block (basic bonds only)
pub(super) fn bond_block(
    input: &mut &[u8],
    bond_count: u32,
    line_offset: u32,
    flags: CtabParseFlags,
) -> ModalResult<(Vec<Bond>, u32), ParseError> {
    let mut bonds = Vec::with_capacity(bond_count as usize);
    for line_index in 0..bond_count {
        let physical_line = line_offset + line_index;
        let mut line = next_line(input).map_err(|_| {
            ErrMode::Cut(ParseError::UnexpectedEof {
                line: physical_line,
                block: "bond",
            })
        })?;
        let result = bond_input(flags).parse_next(&mut line).and_then(|value| {
            finish_line(&mut line)?;
            Ok(value)
        });
        let bond = result.map_err(|error| {
            ErrMode::Cut(ParseError::InvalidBondLine {
                line: physical_line,
                col: input_error_column(error, &line),
            })
        })?;
        bonds.push(bond);
    }
    Ok((bonds, line_offset + bond_count))
}

/// Parse extended bond block
pub(super) fn extended_bond_block(
    input: &mut &[u8],
    bond_count: u32,
    line_offset: u32,
    flags: CtabParseFlags,
) -> ModalResult<(Vec<ExtendedBond>, u32), ParseError> {
    let mut bonds = Vec::with_capacity(bond_count as usize);
    for line_index in 0..bond_count {
        let physical_line = line_offset + line_index;
        let mut line = next_line(input).map_err(|_| {
            ErrMode::Cut(ParseError::UnexpectedEof {
                line: physical_line,
                block: "bond",
            })
        })?;
        let result = extended_bond_input(flags)
            .parse_next(&mut line)
            .and_then(|value| {
                finish_line(&mut line)?;
                Ok(value)
            });
        let bond = result.map_err(|error| {
            ErrMode::Cut(ParseError::InvalidBondLine {
                line: physical_line,
                col: input_error_column(error, &line),
            })
        })?;
        bonds.push(bond);
    }
    Ok((bonds, line_offset + bond_count))
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
fn bond_input<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<Input<'inp>, Bond, ErrMode<InputError>> + use<'inp> {
    move |input: &mut Input<'inp>| {
        let bytes: &[u8] = input.as_ref();
        if bytes.len() < 9 {
            return Err(ErrMode::Backtrack(InputError {
                column: bytes.len() as u32,
            }));
        }

        let mut offset;

        let extended_range = flags.contains(CtabParseFlags::EXTENDED_RANGE);
        let skip_unused_fields = flags.contains(CtabParseFlags::SKIP_UNUSED_FIELDS);

        // Atom indices (0-2, 3-5)
        let first_atom = parse_int_opt::<u32>(&bytes[0..3], 0)?
            .and_then(|value| value.checked_sub(1))
            .ok_or(ErrMode::Backtrack(InputError { column: 0 }))?;
        let second_atom = parse_int_opt::<u32>(&bytes[3..6], 3)?
            .and_then(|value| value.checked_sub(1))
            .ok_or(ErrMode::Backtrack(InputError { column: 3 }))?;

        // Bond type (6-8)
        let order_code = parse_int_opt::<u8>(&bytes[6..9], 6)?.unwrap_or(0);
        let order = convert_bond_type_code(order_code, extended_range)
            .map_err(|_| ErrMode::Backtrack(InputError { column: 6 }))?;

        offset = 9;

        // Stereo/direction (9-11)
        let (stereo, wedge) = if bytes.len() >= 12 {
            offset = 12;
            let stereo_code = parse_int_opt::<u8>(&bytes[9..12], 9)?.unwrap_or(0);
            if !matches!(stereo_code, 0 | 1 | 3 | 4 | 6) {
                return Err(ErrMode::Backtrack(InputError { column: 9 }));
            }
            convert_bond_stereo_direction_code(stereo_code)
        } else {
            (None, None)
        };

        // Ignored field xxx (12-14)
        if bytes.len() >= 15 {
            offset = 15;
            validate_unused_n(&bytes[12..15], 1, 3, skip_unused_fields, 12)?;
        }

        // Bond topology and reacting center (15-20) - extended
        let count = ((bytes.len().saturating_sub(15)) / 3).min(2);
        if count > 0 {
            offset = 15 + count * 3;
            validate_unused_n(&bytes[15..15 + count * 3], count, 3, false, 15)?;
        }

        let mut bond = Bond::new(first_atom, second_atom, order);
        // The stereo/direction field applies only to single and double bonds.
        if matches!(order, BondOrder::Single | BondOrder::Double) {
            bond.stereo = stereo;
            bond.wedge = wedge;
        }

        let _: &[u8] = take(offset).parse_next(input)?;
        Ok(bond)
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
fn extended_bond_input<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<Input<'inp>, ExtendedBond, ErrMode<InputError>> + use<'inp> {
    let skip_unused_fields = flags.contains(CtabParseFlags::SKIP_UNUSED_FIELDS);
    let extended_range = flags.contains(CtabParseFlags::EXTENDED_RANGE);
    let allow_wildcards = flags.contains(CtabParseFlags::WILDCARDS);

    move |input: &mut Input<'inp>| {
        let bytes: &[u8] = input.as_ref();
        if bytes.len() < 9 {
            return Err(ErrMode::Backtrack(InputError {
                column: bytes.len() as u32,
            }));
        }

        let mut offset: usize = 9;

        // Atom indices (0-2, 3-5)
        let first_atom = parse_int_opt::<u32>(&bytes[0..3], 0)?
            .and_then(|value| value.checked_sub(1))
            .ok_or(ErrMode::Backtrack(InputError { column: 0 }))?;
        let second_atom = parse_int_opt::<u32>(&bytes[3..6], 3)?
            .and_then(|value| value.checked_sub(1))
            .ok_or(ErrMode::Backtrack(InputError { column: 3 }))?;

        // Bond type (6-8)
        let order_code = parse_int_opt::<u8>(&bytes[6..9], 6)?.unwrap_or(0);
        let order = convert_extended_bond_type_code(order_code, extended_range, allow_wildcards)
            .map_err(|_| ErrMode::Backtrack(InputError { column: 6 }))?;

        // Stereo/direction (9-11)
        let (stereo, wedge) = if bytes.len() >= 12 {
            offset = 12;
            let stereo_code = parse_int_opt::<u8>(&bytes[9..12], 9)?.unwrap_or(0);
            if !matches!(stereo_code, 0 | 1 | 3 | 4 | 6) {
                return Err(ErrMode::Backtrack(InputError { column: 9 }));
            }
            convert_bond_stereo_direction_code(stereo_code)
        } else {
            (None, None)
        };

        // Ignored field xxx (12-14)
        if bytes.len() >= 15 {
            offset = 15;
            validate_unused_n(&bytes[12..15], 1, 3, skip_unused_fields, 12)?;
        }

        // Topology rrr (15-17)
        let topology = if bytes.len() >= 18 {
            offset = 18;
            let val = parse_int_opt::<u8>(&bytes[15..18], 15)?.unwrap_or(0);
            if val > 2 {
                return Err(ErrMode::Backtrack(InputError { column: 15 }));
            }
            convert_bond_topology_code(val)
        } else {
            None
        };

        // Reacting center ccc (18-20)
        let reacting_center = if bytes.len() >= 21 {
            offset = 21;
            let val = parse_int_opt::<i8>(&bytes[18..21], 18)?.unwrap_or(0);
            convert_bond_reacting_center_code(val, extended_range)
                .map_err(|_| ErrMode::Backtrack(InputError { column: 18 }))?
        } else {
            None
        };

        let mut bond = ExtendedBond::new(first_atom, second_atom, order);
        // The stereo/direction field applies only to single and double bonds.
        if matches!(order, BondOrder::Single | BondOrder::Double) {
            bond.stereo = stereo;
            bond.wedge = wedge;
        }
        bond.topology = topology;
        bond.reacting_center = reacting_center;

        let _: &[u8] = take(offset).parse_next(input)?;
        Ok(bond)
    }
}

#[cfg(test)]
mod tests;
