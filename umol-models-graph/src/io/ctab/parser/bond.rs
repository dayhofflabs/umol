//! Bond block parser for CTab files.

use bstr::ByteSlice;
use nom::character::complete::space0;
use nom::combinator::{cond, map, map_res};
use nom::sequence::terminated;
use nom::{error, Err, IResult, Parser};

use super::convert::{
    convert_bond_reacting_center_code, convert_bond_stereo_dir_code, convert_bond_topology_code,
    convert_bond_type_code, convert_bondlike_type_code,
};
use super::utils::{fixed_width_int, fixed_width_int_minus1, fixed_width_padding_n};
use crate::io::ctab::bond::{Bond, BondLike, BondType};
use crate::io::ctab::config::CtabParseFlags;

/// Parse bond inputs with 12-21 characters (s. `bond_input` for more details).
/// Lacks trailing stereo/dir fields (substituted by defaults).
fn bond_input12<'inp, 'fl>(
    input: &'inp [u8],
    flags: &'fl CtabParseFlags,
) -> IResult<&'inp [u8], (usize, usize, Bond)> {
    let strict_padding = flags.contains(CtabParseFlags::STRICT_PADDING);
    let first_atom = fixed_width_int_minus1::<usize>(3);
    let second_atom = fixed_width_int_minus1::<usize>(3);
    let bond_type = map_res(fixed_width_int::<u8>(3), move |code| {
        convert_bond_type_code(code)
    });
    let stereo_dir = map_res(fixed_width_int::<u8>(3), |code| {
        convert_bond_stereo_dir_code(code, false)
    });
    let n = input.len().saturating_sub(12) / 3;
    let padding1 = fixed_width_padding_n(n, 3, strict_padding);

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
fn bond_input9<'inp, 'fl>(
    input: &'inp [u8],
    _flags: &'fl CtabParseFlags,
) -> IResult<&'inp [u8], (usize, usize, Bond)> {
    map(
        (
            fixed_width_int_minus1::<usize>(3),
            fixed_width_int_minus1::<usize>(3),
            map_res(fixed_width_int::<u8>(3), convert_bond_type_code),
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
pub fn bond_input<'inp, 'fl>(
    flags: &'fl CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (usize, usize, Bond), Error = error::Error<&'inp [u8]>>
       + use<'inp, 'fl> {
    move |input: &'inp [u8]| {
        let input = input.trim_end_with(|c| c == '\r' || c == '\n');
        let len = input.len();
        let parser = match len {
            10.. => bond_input12,
            9 => bond_input9,
            _ => return Err(Err::Error(error::Error::new(input, error::ErrorKind::Eof))),
        };

        terminated(move |input| parser(input, flags), space0).parse(input)
    }
}

fn bondlike_input_inner<'inp, 'fl>(
    flags: &'fl CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (usize, usize, BondLike), Error = error::Error<&'inp [u8]>>
       + use<'inp, 'fl> {
    let strict_padding = flags.contains(CtabParseFlags::STRICT_PADDING);
    // Bond type - allow zero-order/high-order bonds and queries based on flags
    let extended_range = flags.contains(CtabParseFlags::EXTENDED_RANGE);
    let allow_queries = flags.contains(CtabParseFlags::QUERIES);
    move |input: &'inp [u8]| {
        // Atom indices
        let (i, first_atom) = fixed_width_int_minus1::<usize>(3).parse(input)?;
        let (i, second_atom) = fixed_width_int_minus1::<usize>(3).parse(i)?;

        let (i, bond_type) = map_res(fixed_width_int::<u8>(3), move |code| {
            convert_bondlike_type_code(code, extended_range, allow_queries)
        })
        .parse(i)?;

        // Stereo/dir
        let (i, stereo_dir) = cond(
            i.len() >= 3,
            map_res(fixed_width_int::<u8>(3), |code| {
                convert_bond_stereo_dir_code(code, true)
            }),
        )
        .parse(i)?;

        // Ignore xxx field
        let (i, _) = cond(
            !i.is_empty(),
            fixed_width_padding_n((i.len() / 3).min(1), 3, strict_padding),
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

pub fn bondlike_input<'inp, 'fl>(
    flags: &'fl CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (usize, usize, BondLike), Error = error::Error<&'inp [u8]>>
       + use<'inp, 'fl> {
    move |input: &'inp [u8]| {
        let input = input.trim_end_with(|c| c == '\r' || c == '\n');
        terminated(bondlike_input_inner(flags), space0).parse(input)
    }
}

#[cfg(test)]
mod tests;
