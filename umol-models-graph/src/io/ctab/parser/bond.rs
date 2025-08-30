//! Bond block parser for CTab files.

use nom::character::complete::space0;
use nom::combinator::{cond, map, map_res};
use nom::error;
use nom::sequence::terminated;
use nom::{Err, IResult, Parser};

use crate::io::config::ParseFlags;
use crate::io::ctab::bond::{Bond, BondLike, BondType};

use super::convert::{
    convert_bond_reacting_center_code, convert_bond_stereo_dir_code, convert_bond_topology_code,
    convert_bond_type_code, convert_bondlike_type_code,
};
use super::utils::{fixed_width_int, fixed_width_int_minus1, fixed_width_padding_n};

/// Parse bond inputs with 12-21 characters (s. `bond_input` for more details).
/// Lacks trailing stereo/dir fields (substituted by defaults).
fn bond_input12(input: &[u8], flags: ParseFlags) -> IResult<&[u8], (usize, usize, Bond)> {
    let allow_unicode = flags.contains(ParseFlags::UNICODE);
    let strict_padding = flags.contains(ParseFlags::STRICT_PADDING);
    let first_atom = fixed_width_int_minus1::<usize>(3, allow_unicode);
    let second_atom = fixed_width_int_minus1::<usize>(3, allow_unicode);
    let bond_type = map_res(
        fixed_width_int::<u8>(3, allow_unicode),
        convert_bond_type_code,
    );
    let stereo_dir = map_res(fixed_width_int::<u8>(3, allow_unicode), |code| {
        convert_bond_stereo_dir_code(code, false)
    });
    let n = input.len().saturating_sub(12) / 3;
    let padding1 = fixed_width_padding_n(n, 3, allow_unicode, strict_padding);

    map(
        (
            first_atom,
            second_atom,
            bond_type,
            terminated(stereo_dir, padding1),
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

/// Parse  bond inputs with 9 characters (s. `bond_input` for more details).
/// Lacks trailing stereo/dir fields (substituted by defaults).
fn bond_input9(input: &[u8], flags: ParseFlags) -> IResult<&[u8], (usize, usize, Bond)> {
    let allow_unicode = flags.contains(ParseFlags::UNICODE);
    map(
        (
            fixed_width_int_minus1::<usize>(3, allow_unicode),
            fixed_width_int_minus1::<usize>(3, allow_unicode),
            map_res(fixed_width_int::<u8>(3, allow_unicode), |code| {
                convert_bond_type_code(code)
            }),
        ),
        |(first_atom, second_atom, bond_type)| (first_atom, second_atom, Bond::new(bond_type)),
    )
    .parse(input)
}

/// Parse bond input (optimized for performance)
/// Fails immediately on query bond properties. For parsing all bond types, see bondlike_input.
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
pub fn bond_input<'a>(
    flags: ParseFlags,
) -> impl Parser<&'a [u8], Output = (usize, usize, Bond), Error = error::Error<&'a [u8]>> {
    move |input: &'a [u8]| {
        let len = input.len();
        let parser = match len {
            10.. => bond_input12,
            9 => bond_input9,
            _ => return Err(Err::Error(error::Error::new(input, error::ErrorKind::Eof))),
        };
        terminated(move |input| parser(input, flags), space0).parse(input)
    }
}

fn bondlike_input_inner<'a>(
    flags: ParseFlags,
) -> impl Parser<&'a [u8], Output = (usize, usize, BondLike), Error = error::Error<&'a [u8]>> {
    let allow_unicode = flags.contains(ParseFlags::UNICODE);
    let strict_padding = flags.contains(ParseFlags::STRICT_PADDING);
    move |input: &'a [u8]| {
        // Atom indices
        let (i, first_atom) = fixed_width_int_minus1::<usize>(3, allow_unicode).parse(input)?;
        let (i, second_atom) = fixed_width_int_minus1::<usize>(3, allow_unicode).parse(i)?;

        // Bond type
        let (i, bond_type) = map_res(fixed_width_int::<u8>(3, allow_unicode), |code| {
            convert_bondlike_type_code(code, false)
        })
        .parse(i)?;

        // Stereo/dir
        let (i, stereo_dir) = cond(
            i.len() >= 3,
            map_res(fixed_width_int::<u8>(3, allow_unicode), |code| {
                convert_bond_stereo_dir_code(code, true)
            }),
        )
        .parse(i)?;

        // Ignore xxx field
        let (i, _) = cond(
            i.len() > 0,
            fixed_width_padding_n((i.len() / 3).min(1), 3, allow_unicode, strict_padding),
        )
        .parse(i)?;

        // Topology, reacting center
        let (i, topology) = cond(
            i.len() >= 3,
            map_res(fixed_width_int::<u8>(3, allow_unicode), |code| {
                convert_bond_topology_code(code, true)
            }),
        )
        .parse(i)?;
        let (i, reacting_center) = cond(
            i.len() >= 3,
            map_res(fixed_width_int::<i8>(3, allow_unicode), |code| {
                convert_bond_reacting_center_code(code, true)
            }),
        )
        .parse(i)?;

        let mut bond = BondLike::new(bond_type);
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

pub fn bondlike_input<'a>(
    flags: ParseFlags,
) -> impl Parser<&'a [u8], Output = (usize, usize, BondLike), Error = error::Error<&'a [u8]>> {
    terminated(bondlike_input_inner(flags), space0)
}

#[cfg(test)]
mod tests;
