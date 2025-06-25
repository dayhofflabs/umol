//! Bond block parser for CTab files.

use nom::bytes::complete::take;
use nom::character::complete::multispace0;
use nom::combinator::{all_consuming, map, map_res};
use nom::error;
use nom::sequence::terminated;
use nom::{IResult, Parser};

use crate::bond::{Bond, BondDir, BondStereo, BondType};

use super::convert::{convert_bond_stereo_dir_code, convert_bond_type_code_standard};
use super::utils::{fixed_width_int, fixed_width_int_minus1};

fn bond_input9(input: &[u8]) -> IResult<&[u8], (usize, usize, Bond)> {
    map(
        (
            fixed_width_int_minus1::<usize>(3),
            fixed_width_int_minus1::<usize>(3),
            map_res(fixed_width_int::<u8>(3), |code| {
                convert_bond_type_code_standard(code)
            }),
        ),
        |(first_atom, second_atom, bond_type)| (first_atom, second_atom, Bond::new(bond_type)),
    )
    .parse(input)
}

fn bond_input12(input: &[u8]) -> IResult<&[u8], (usize, usize, Bond)> {
    map(
        (
            fixed_width_int_minus1::<usize>(3),
            fixed_width_int_minus1::<usize>(3),
            map_res(fixed_width_int::<u8>(3), |code| {
                convert_bond_type_code_standard(code)
            }),
            map_res(fixed_width_int::<u8>(3), |code| {
                convert_bond_stereo_dir_code(code)
            }),
        ),
        |(first_atom, second_atom, bond_type, (stereo, dir))| {
            let mut bond = Bond::new(bond_type);
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

fn bond_input21(input: &[u8]) -> IResult<&[u8], (usize, usize, Bond)> {
    terminated(
        map(
            (
                fixed_width_int_minus1::<usize>(3),
                fixed_width_int_minus1::<usize>(3),
                map_res(fixed_width_int::<u8>(3), |code| {
                    convert_bond_type_code_standard(code)
                }),
                map_res(fixed_width_int::<u8>(3), |code| {
                    convert_bond_stereo_dir_code(code)
                }),
            ),
            |(first_atom, second_atom, bond_type, (stereo, dir))| {
                let mut bond = Bond::new(bond_type);
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
/// | ttt   | bond type            | 1..=8      | Query          |
/// | sss   | bond stereo          | 0..=6      | Generic        |
/// | rrr   | bond topology        | 0..=2      | Query          |
/// | ccc   | bond reacting center | 0..=3      | Reaction,Query |
/// --------------------------------------------------------------
pub fn bond_input<'a>(
) -> impl Parser<&'a [u8], Output = (usize, usize, Bond), Error = error::Error<&'a [u8]>> {
    move |input: &'a [u8]| {
        let len = input.len();
        let parser = match len {
            21.. => bond_input21,
            10..=20 => bond_input12,
            9 => bond_input9,
            _ => {
                return Err(nom::Err::Error(error::Error::new(
                    input,
                    error::ErrorKind::Eof,
                )))
            }
        };
        all_consuming(terminated(parser, multispace0)).parse(input)
    }
}

#[cfg(test)]
mod tests;
