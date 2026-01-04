//! Bond block parser for CTab files.
//!
//! Parses bond blocks from CTFile format and produces TableIR types directly.

use bstr::ByteSlice;
use nom::character::complete::space0;
use nom::combinator::{all_consuming, cond, map, map_res};
use nom::error::{Error as NomError, ErrorKind as NomErrorKind};
use nom::sequence::terminated;
use nom::{Err as NomErr, IResult, Parser};

use super::convert::{
    convert_bond_reacting_center_code, convert_bond_stereo_direction_code,
    convert_bond_topology_code, convert_bond_type_code, convert_extended_bond_type_code,
};
use super::utils::{fixed_width_int, fixed_width_int_minus1, fixed_width_unused_n};
use crate::io::ctfile::config::CtabParseFlags;
use crate::io::ctfile::error::ParseError;
use crate::table_ir::bond::{Bond, BondOrder, ExtendedBond};

/// Parse bond inputs with 12-21 characters (s. `bond_input` for more details).
/// Lacks trailing stereo/dir fields (substituted by defaults).
fn bond_input12<'inp>(
    input: &'inp [u8],
    flags: CtabParseFlags,
) -> IResult<&'inp [u8], (usize, usize, Bond)> {
    let extended_range = flags.contains(CtabParseFlags::EXTENDED_RANGE);
    let skip_unused_fields = flags.contains(CtabParseFlags::SKIP_UNUSED_FIELDS);
    let first_atom = fixed_width_int_minus1::<usize>(3);
    let second_atom = fixed_width_int_minus1::<usize>(3);
    let bond_type = map_res(fixed_width_int::<u8>(3), move |code| {
        convert_bond_type_code(code, extended_range)
    });
    let stereo_direction = map_res(fixed_width_int::<u8>(3), |code| {
        convert_bond_stereo_direction_code(code, false)
    });
    let n = input.len().saturating_sub(12) / 3;
    let padding1 = fixed_width_unused_n(n, 3, skip_unused_fields);

    map(
        (
            first_atom,
            second_atom,
            bond_type,
            terminated(stereo_direction, padding1),
        ),
        |(first_atom, second_atom, order, (stereo, dir))| {
            let mut bond = Bond::with_order(order);
            match order {
                BondOrder::Single => bond.direction = dir,
                BondOrder::Double => bond.stereo = stereo,
                _ => (),
            }
            (first_atom, second_atom, bond)
        },
    )
    .parse(input)
}

/// Parse  bond inputs with 9 characters (s. `bond_input` for more details).
/// Lacks trailing stereo/dir fields (substituted by defaults).
fn bond_input9<'inp>(
    input: &'inp [u8],
    flags: CtabParseFlags,
) -> IResult<&'inp [u8], (usize, usize, Bond)> {
    let extended_range = flags.contains(CtabParseFlags::EXTENDED_RANGE);
    let first_atom = fixed_width_int_minus1::<usize>(3);
    let second_atom = fixed_width_int_minus1::<usize>(3);
    let bond_type = map_res(fixed_width_int::<u8>(3), move |code| {
        convert_bond_type_code(code, extended_range)
    });
    map(
        (first_atom, second_atom, bond_type),
        |(first_atom, second_atom, order)| (first_atom, second_atom, Bond::with_order(order)),
    )
    .parse(input)
}

/// Parse bond input (optimized for performance)
/// Fails immediately on query bond properties. For parsing all bond types, see extended_bond_input.
///
/// "111222tttsssxxxrrrccc" (21 characters wide)
///
/// *Values in the bond block*
/// --------------------------------------------------------------
/// | Field | Meaning              | Values     | Notes          |
/// |-------|----------------------|------------|----------------|
/// | 111   | first atom           | 1..=aaa    | Generic        |
/// | 222   | second atom          | 1..=aaa    | Generic        |
/// | ttt   | bond type            | 1..=8      | Generic,Query  |
/// | sss   | bond stereo          | 0..=6      | Generic        |
/// | xxx   | ignored              |            |                |
/// | rrr   | bond topology        | 0..=2      | Query          |
/// | ccc   | bond reacting center | 0..=3      | Reaction,Query |
/// --------------------------------------------------------------
///
pub fn bond_input<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (usize, usize, Bond), Error = NomError<&'inp [u8]>> + use<'inp>
{
    move |input: &'inp [u8]| {
        let input = input.trim_end_with(|c| c == '\r' || c == '\n');
        let len = input.len();
        let parser = match len {
            10.. => bond_input12,
            9 => bond_input9,
            _ => return Err(NomErr::Error(NomError::new(input, NomErrorKind::Eof))),
        };

        terminated(move |input| parser(input, flags), space0).parse(input)
    }
}

fn extended_bond_input_inner<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (usize, usize, ExtendedBond), Error = NomError<&'inp [u8]>>
       + use<'inp> {
    let skip_unused_fields = flags.contains(CtabParseFlags::SKIP_UNUSED_FIELDS);
    let extended_range = flags.contains(CtabParseFlags::EXTENDED_RANGE);
    let allow_wildcards = flags.contains(CtabParseFlags::WILDCARDS);
    move |input: &'inp [u8]| {
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
            map_res(fixed_width_int::<u8>(3), |code| {
                convert_bond_stereo_direction_code(code, true)
            }),
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
            map_res(fixed_width_int::<u8>(3), |code| {
                convert_bond_topology_code(code, true)
            }),
        )
        .parse(i)?;
        let (i, reacting_center) = cond(
            i.len() >= 3,
            map_res(fixed_width_int::<i8>(3), |code| {
                convert_bond_reacting_center_code(code, true, extended_range)
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

pub fn extended_bond_input<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (usize, usize, ExtendedBond), Error = NomError<&'inp [u8]>>
       + use<'inp> {
    move |input: &'inp [u8]| {
        let input = input.trim_end_with(|c| c == '\r' || c == '\n');
        terminated(extended_bond_input_inner(flags), space0).parse(input)
    }
}

/// Parse bond block (basic bonds only)
pub(super) fn bond_block<'inp>(
    bond_count: u32,
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (Vec<(usize, usize, Bond)>, u32), Error = ParseError> + use<'inp>
{
    move |input: &'inp [u8]| {
        let mut bonds = Vec::with_capacity(bond_count as usize);
        let mut lines_iter = input.lines_with_terminator();
        let mut offset = 0;

        for i in 0..bond_count {
            let line = lines_iter.next().ok_or_else(|| {
                NomErr::Error(ParseError::UnexpectedEof {
                    line: line_offset + i,
                    block: "bond",
                })
            })?;

            let (_, (atom1, atom2, bond)) = all_consuming(bond_input(flags))
                .parse(line)
                .map_err(|e| NomErr::Error(ParseError::bond_from_nom(e, line_offset + i, line)))?;
            bonds.push((atom1, atom2, bond));
            offset += line.len();
        }

        let remaining = &input[offset..];
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
        let mut lines_iter = input.lines_with_terminator();
        let mut offset = 0;

        for i in 0..bond_count {
            let line = lines_iter.next().ok_or_else(|| {
                NomErr::Error(ParseError::UnexpectedEof {
                    line: line_offset + i,
                    block: "bond",
                })
            })?;

            let (_, (atom1, atom2, bond)) = all_consuming(extended_bond_input(flags))
                .parse(line)
                .map_err(|e| {
                NomErr::Error(ParseError::bond_from_nom(e, line_offset + i, line))
            })?;
            bonds.push((atom1, atom2, bond));
            offset += line.len();
        }

        let remaining = &input[offset..];
        Ok((remaining, (bonds, line_offset + bond_count)))
    }
}

#[cfg(test)]
mod tests;
