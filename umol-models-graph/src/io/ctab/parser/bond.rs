//! Bond block parser for CTab files.

use nom::bytes::complete::take;
use nom::character::complete::space0;
use nom::combinator::{cond, map, map_res};
use nom::error;
use nom::sequence::terminated;
use nom::{Err, IResult, Parser};

use crate::io::ctab::bond::{Bond, BondStandard, BondType};

use super::convert::{
    convert_bond_reacting_center_code, convert_bond_stereo_dir_code, convert_bond_topology_code,
    convert_bond_type_code, convert_bond_type_code_standard,
};
use super::utils::{fixed_width_int, fixed_width_int_minus1};

fn bond_input_standard21(input: &[u8]) -> IResult<&[u8], (usize, usize, BondStandard)> {
    terminated(
        map(
            (
                fixed_width_int_minus1::<usize>(3),
                fixed_width_int_minus1::<usize>(3),
                map_res(fixed_width_int::<u8>(3), |code| {
                    convert_bond_type_code_standard(code)
                }),
                map_res(fixed_width_int::<u8>(3), |code| {
                    convert_bond_stereo_dir_code(code, false)
                }),
            ),
            |(first_atom, second_atom, bond_type, (stereo, dir))| {
                let mut bond = BondStandard::new(bond_type);
                match bond.bond_type {
                    BondType::Single => bond.dir = dir,
                    BondType::Double => bond.stereo = stereo,
                    _ => (),
                }
                (first_atom, second_atom, bond)
            },
        ),
        take(9usize), // ignore rrrccc fields
    )
    .parse(input)
}

/// Parse standard bond inputs with 12 characters (s. `bond_input` for more details).
/// Lacks trailing stereo/dir fields (substituted by defaults).
fn bond_input_standard12(input: &[u8]) -> IResult<&[u8], (usize, usize, BondStandard)> {
    map(
        (
            fixed_width_int_minus1::<usize>(3),
            fixed_width_int_minus1::<usize>(3),
            map_res(fixed_width_int::<u8>(3), |code| {
                convert_bond_type_code_standard(code)
            }),
            map_res(fixed_width_int::<u8>(3), |code| {
                convert_bond_stereo_dir_code(code, false)
            }),
        ),
        |(first_atom, second_atom, bond_type, (stereo, dir))| {
            let mut bond = BondStandard::new(bond_type);
            match bond.bond_type {
                BondType::Single => bond.dir = dir,
                BondType::Double => bond.stereo = stereo,
                _ => (),
            }
            (first_atom, second_atom, bond)
        },
    )
    .parse(input)
}

/// Parse standard bond inputs with 9 characters (s. `bond_input` for more details).
/// Lacks trailing stereo/dir fields (substituted by defaults).
fn bond_input_standard9(input: &[u8]) -> IResult<&[u8], (usize, usize, BondStandard)> {
    map(
        (
            fixed_width_int_minus1::<usize>(3),
            fixed_width_int_minus1::<usize>(3),
            map_res(fixed_width_int::<u8>(3), |code| {
                convert_bond_type_code_standard(code)
            }),
        ),
        |(first_atom, second_atom, bond_type)| {
            (first_atom, second_atom, BondStandard::new(bond_type))
        },
    )
    .parse(input)
}

/// Parse bond input (optimized for performance)
/// Fails immediately on non-standard bond properties. For parsing all bond types, see bond_like_input.
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
/// | rrr   | bond topology        | 0..=2      | Query          |
/// | ccc   | bond reacting center | 0..=3      | Reaction,Query |
/// --------------------------------------------------------------
pub fn bond_input_standard<'a>(
) -> impl Parser<&'a [u8], Output = (usize, usize, BondStandard), Error = error::Error<&'a [u8]>> {
    move |input: &'a [u8]| {
        let len = input.len();
        let parser = match len {
            21.. => bond_input_standard21,
            10..=20 => bond_input_standard12,
            9 => bond_input_standard9,
            _ => {
                return Err(Err::Error(error::Error::new(
                    input,
                    error::ErrorKind::Eof,
                )))
            }
        };
        terminated(parser, space0).parse(input)
    }
}

fn bond_input_inner<'a>(
) -> impl Parser<&'a [u8], Output = (usize, usize, Bond), Error = error::Error<&'a [u8]>> {
    move |input: &'a [u8]| {
        let (i, first_atom) = fixed_width_int_minus1::<usize>(3).parse(input)?;
        let (i, second_atom) = fixed_width_int_minus1::<usize>(3).parse(i)?;
        let (i, bond_type) = map_res(fixed_width_int::<u8>(3), |code| {
            convert_bond_type_code(code, false)
        })
        .parse(i)?;

        let (i, stereo_dir) = cond(
            i.len() >= 3,
            map_res(fixed_width_int::<u8>(3), |code| {
                convert_bond_stereo_dir_code(code, true)
            }),
        )
        .parse(i)?;

        let (i, _) = cond(i.len() >= 3, take(3usize)).parse(i)?; // xxx field

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
                convert_bond_reacting_center_code(code, true)
            }),
        )
        .parse(i)?;

        let mut bond = Bond::new(bond_type);
        if let Some((stereo_val, dir_val)) = stereo_dir {
            match bond.bond_type {
                BondType::Single => {
                    bond.dir = dir_val;
                }
                BondType::Double => {
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

pub fn bond_input<'a>(
) -> impl Parser<&'a [u8], Output = (usize, usize, Bond), Error = error::Error<&'a [u8]>> {
    move |input: &'a [u8]| {
        terminated(bond_input_inner(), space0).parse(input)
    }
}

#[cfg(test)]
mod tests;
