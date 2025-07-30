//! Atom block parser for CTab files.

use nom::bytes::complete::take;
use nom::character::complete::space0;
use nom::combinator::{complete, cond, map, map_res, verify};
use nom::sequence::{delimited, preceded, terminated};
use nom::{error, IResult, Parser};
use umol_data::{Element, NamedIsotope};

use super::utils::is_blanks_or_zeros;
use crate::io::ctab::atom::{Atom, AtomList, AtomStandard, AtomSymbol, Point3D};
use crate::io::ctab::parser::utils::trim_whitespace_3char;
use crate::io::ctab::query::QueryAtom;
use crate::io::ctab::rgroup::RGroup;

use super::convert::{
    convert_atom_charge_code, convert_atom_exact_change_flag_code,
    convert_atom_hydrogen_count_code, convert_atom_inversion_flag_code,
    convert_atom_mass_diff_code, convert_atom_stereo_care_code, convert_atom_stereo_parity_code,
    convert_atom_symbol_mass_diff, convert_atom_valence_code,
};
use super::utils::{
    fixed_width_float, fixed_width_int, fixed_width_int_in_range, fixed_width_int_in_range_opt,
    repeat,
};

/// Parse atom symbol (standard atoms, Element and NamedIsotope only)
/// Returns error for non-standard atom symbols (L, A, Q, *, LP, R#).
/// See `atom_symbol` for more details.
fn atom_symbol_standard<'a>(
) -> impl Parser<&'a [u8], Output = AtomSymbol, Error = error::Error<&'a [u8]>> {
    map_res(take(3usize), |s| {
        let s = trim_whitespace_3char(s);
        if let Some(element) = Element::from_symbol_bytes(s) {
            Ok(AtomSymbol::Element(element))
        } else if let Some(isotope) = NamedIsotope::from_symbol_bytes(s) {
            Ok(AtomSymbol::NamedIsotope(isotope))
        } else {
            Err(error::Error::new(s, error::ErrorKind::MapRes))
        }
    })
}

/// Parse atom symbol (all atom types allowed in MOL specification).
/// --------------------------------------------------------------------------------
/// | Symbol      | Type          | Parser* | Notes                                |
/// --------------------------------------------------------------------------------
/// | H-Og        | Element       | S, A    |                                      |
/// | D, T        | Named Isotope | S, A    | Heavy H isotopes as extension        |
/// | L           | Atom List     | A       | Query molecules                      |
/// | *,A,Q,X,M   | Query Atom    | A       | Query molecules, rarely in oligomers |
/// | AH,QH,XH,MH | Query Atom    | A       | Query molecules, CXSMILES extension  |
/// | LP          | Lone Pair     | A       | Rarely used                          |
/// | R#          | R Group       | A       | Query molecules                      |
/// --------------------------------------------------------------------------------
/// | *Parsers: S: standard, A: all                                                |
/// --------------------------------------------------------------------------------
///
fn atom_symbol<'a>() -> impl Parser<&'a [u8], Output = AtomSymbol, Error = error::Error<&'a [u8]>> {
    map_res(take(3usize), |s| {
        let s = trim_whitespace_3char(s);
        if let Some(element) = Element::from_symbol_bytes(s) {
            Ok(AtomSymbol::Element(element))
        } else if let Some(isotope) = NamedIsotope::from_symbol_bytes(s) {
            Ok(AtomSymbol::NamedIsotope(isotope))
        } else if let Some(rgroup) = RGroup::from_symbol_bytes(s) {
            Ok(AtomSymbol::RGroup(rgroup))
        } else if let Some(query) = QueryAtom::from_symbol_bytes(s) {
            Ok(AtomSymbol::Query(query))
        } else if s == b"L" {
            Ok(AtomSymbol::AtomList(AtomList {
                elements: vec![],
                exclusion: false,
            }))
        } else if s == b"LP" {
            Ok(AtomSymbol::LonePair)
        } else {
            Err(error::Error::new(s, error::ErrorKind::MapRes))
        }
    })
}

/// Parse standard atom inputs with 52-69 characters (s. `atom_input` for more details).
/// Includes atom mapping number. Ignores whitespace padding.
fn atom_input_standard69(input: &[u8]) -> IResult<&[u8], AtomStandard, error::Error<&[u8]>> {
    let x = fixed_width_float::<f64>(10, 4);
    let y = fixed_width_float::<f64>(10, 4);
    let z = fixed_width_float::<f64>(10, 4);
    let symbol = atom_symbol_standard();
    let mass_diff = map(fixed_width_int_in_range_opt::<i8, _>(2, -3..=4), |opt| {
        convert_atom_mass_diff_code(opt.unwrap_or(0))
    });
    let charge_radical = map(fixed_width_int_in_range_opt::<u8, _>(3, 0..=7), |opt| {
        convert_atom_charge_code(opt.unwrap_or(0))
    });
    let valence = map_res(
        fixed_width_int_in_range::<u8, _>(3, 0..=15),
        convert_atom_valence_code,
    );
    let atom_map_num = fixed_width_int_in_range_opt::<u32, _>(3, 1..=999);

    map(
        (
            x,
            y,
            z,
            preceded(take(1usize), symbol),
            complete(mass_diff),
            complete(charge_radical),
            preceded(
                repeat(3usize, verify(take(3usize), is_blanks_or_zeros)),
                complete(valence),
            ),
            delimited(
                take(9usize),
                complete(atom_map_num),
                repeat(2usize, verify(take(3usize), is_blanks_or_zeros)),
            ),
        ),
        |(x, y, z, symbol, mass_diff, charge_radical, valence, atom_map_num)| {
            let (element, isotope_mass) = convert_atom_symbol_mass_diff(symbol, mass_diff);
            let (charge, radical) = charge_radical;
            AtomStandard {
                element,
                charge,
                radical,
                isotope_mass,
                stereo_parity: None,
                hydrogen_count: None,
                valence,
                atom_map_num,
                position: Some(Point3D::new(x, y, z)),
                properties: std::collections::HashMap::new(),
            }
        },
    )
    .parse(input)
}

/// Parse standard atom inputs with 49-51 characters (s. `atom_input` for more details).
/// Includes mass difference, charge/radical and valence fields.
fn atom_input_standard51(input: &[u8]) -> IResult<&[u8], AtomStandard, error::Error<&[u8]>> {
    let x = fixed_width_float::<f64>(10, 4);
    let y = fixed_width_float::<f64>(10, 4);
    let z = fixed_width_float::<f64>(10, 4);
    let symbol = atom_symbol_standard();
    let mass_diff = map(fixed_width_int_in_range_opt::<i8, _>(2, -3..=4), |opt| {
        convert_atom_mass_diff_code(opt.unwrap_or(0))
    });
    let charge_radical = map(fixed_width_int_in_range_opt::<u8, _>(3, 0..=7), |opt| {
        convert_atom_charge_code(opt.unwrap_or(0))
    });
    let stereo_parity = map_res(
        fixed_width_int_in_range::<u8, _>(3, 0..=3),
        convert_atom_stereo_parity_code,
    );
    let valence = map_res(
        fixed_width_int_in_range::<u8, _>(3, 0..=15),
        convert_atom_valence_code,
    );

    map(
        (
            x,
            y,
            z,
            preceded(take(1usize), symbol),
            complete(mass_diff),
            complete(charge_radical),
            complete(stereo_parity),
            preceded(
                repeat(2usize, verify(take(3usize), is_blanks_or_zeros)),
                complete(valence),
            ),
        ),
        |(x, y, z, symbol, mass_diff, charge_radical, stereo_parity, valence)| {
            let (element, isotope_mass) = convert_atom_symbol_mass_diff(symbol, mass_diff);
            let (charge, radical) = charge_radical;
            AtomStandard {
                element,
                charge,
                radical,
                isotope_mass,
                stereo_parity,
                hydrogen_count: None,
                valence,
                atom_map_num: None,
                position: Some(Point3D::new(x, y, z)),
                properties: std::collections::HashMap::new(),
            }
        },
    )
    .parse(input)
}

/// Parse standard atom inputs with 42-48 characters including up to 6 characters of ignored data
/// (s. `atom_input` for more details). Lacks trailing valence and atom mapping fields
/// (substituted by defaults).
fn atom_input_standard42(input: &[u8]) -> IResult<&[u8], AtomStandard, error::Error<&[u8]>> {
    let x = fixed_width_float::<f64>(10, 4);
    let y = fixed_width_float::<f64>(10, 4);
    let z = fixed_width_float::<f64>(10, 4);
    let symbol = atom_symbol_standard();
    let mass_diff = map(fixed_width_int_in_range_opt::<i8, _>(2, -3..=4), |opt| {
        convert_atom_mass_diff_code(opt.unwrap_or(0))
    });
    let charge_radical = map(fixed_width_int_in_range_opt::<u8, _>(3, 0..=7), |opt| {
        convert_atom_charge_code(opt.unwrap_or(0))
    });
    let stereo_parity = map_res(
        fixed_width_int_in_range::<u8, _>(3, 0..=3),
        convert_atom_stereo_parity_code,
    );

    let (input, (x, y, z, symbol, mass_diff, charge_radical, stereo_parity)) = (
        x,
        y,
        z,
        preceded(take(1usize), symbol),
        complete(mass_diff),
        complete(charge_radical),
        complete(stereo_parity),
    )
        .parse(input)?;

    // Verify that ignored block has up to 2 fields of width 3 containing blanks or zeros
    let n_to_take = 6.min(input.len());
    let (remaining, _) =
        repeat(n_to_take / 3, verify(take(3usize), is_blanks_or_zeros)).parse(input)?;

    let (element, isotope_mass) = convert_atom_symbol_mass_diff(symbol, mass_diff);
    let (charge, radical) = charge_radical;
    Ok((
        remaining,
        AtomStandard {
            element,
            charge,
            radical,
            isotope_mass,
            stereo_parity,
            hydrogen_count: None,
            valence: None,
            atom_map_num: None,
            position: Some(Point3D::new(x, y, z)),
            properties: std::collections::HashMap::new(),
        },
    ))
}

/// Parse standard atom inputs with 37-41 characters (s. `atom_input` for more details).
/// Lacks trailing stereo parity, valence and atom mapping fields (substituted by defaults).
fn atom_input_standard39(input: &[u8]) -> IResult<&[u8], AtomStandard, error::Error<&[u8]>> {
    let x = fixed_width_float::<f64>(10, 4);
    let y = fixed_width_float::<f64>(10, 4);
    let z = fixed_width_float::<f64>(10, 4);
    let symbol = atom_symbol_standard();
    let mass_diff = map(fixed_width_int_in_range_opt::<i8, _>(2, -3..=4), |opt| {
        convert_atom_mass_diff_code(opt.unwrap_or(0))
    });
    let charge_radical = map(fixed_width_int_in_range_opt::<u8, _>(3, 0..=7), |opt| {
        convert_atom_charge_code(opt.unwrap_or(0))
    });

    map(
        (
            x,
            y,
            z,
            preceded(take(1usize), symbol),
            complete(mass_diff),
            complete(charge_radical),
        ),
        |(x, y, z, symbol, mass_diff, charge_radical)| {
            let (element, isotope_mass) = convert_atom_symbol_mass_diff(symbol, mass_diff);
            let (charge, radical) = charge_radical;
            AtomStandard {
                element,
                charge,
                radical,
                isotope_mass,
                stereo_parity: None,
                hydrogen_count: None,
                valence: None,
                atom_map_num: None,
                position: Some(Point3D::new(x, y, z)),
                properties: std::collections::HashMap::new(),
            }
        },
    )
    .parse(input)
}

/// Parse standard atom inputs with 35-36 characters (s. `atom_input` for more details).
/// Lacks trailing charge/radical, valence and atom mapping fields (substituted by defaults).
fn atom_input_standard36(input: &[u8]) -> IResult<&[u8], AtomStandard, error::Error<&[u8]>> {
    let x = fixed_width_float::<f64>(10, 4);
    let y = fixed_width_float::<f64>(10, 4);
    let z = fixed_width_float::<f64>(10, 4);
    let symbol = atom_symbol_standard();
    let mass_diff = map(fixed_width_int_in_range_opt::<i8, _>(2, -3..=4), |opt| {
        convert_atom_mass_diff_code(opt.unwrap_or(0))
    });

    map(
        (x, y, z, preceded(take(1usize), symbol), complete(mass_diff)),
        |(x, y, z, symbol, mass_diff)| {
            let (element, isotope_mass) = convert_atom_symbol_mass_diff(symbol, mass_diff);
            AtomStandard {
                element,
                charge: 0,
                radical: None,
                isotope_mass,
                stereo_parity: None,
                hydrogen_count: None,
                valence: None,
                atom_map_num: None,
                position: Some(Point3D::new(x, y, z)),
                properties: std::collections::HashMap::new(),
            }
        },
    )
    .parse(input)
}

/// Parse standard atom inputs with 34 characters (s. `atom_input` for more details).
/// Lacks trailing mass difference, charge/radical, valence and atom mapping fields
/// (substituted by defaults).
fn atom_input_standard34(input: &[u8]) -> IResult<&[u8], AtomStandard, error::Error<&[u8]>> {
    let x = fixed_width_float::<f64>(10, 4);
    let y = fixed_width_float::<f64>(10, 4);
    let z = fixed_width_float::<f64>(10, 4);
    let symbol = atom_symbol_standard();

    map(
        (x, y, z, preceded(take(1usize), symbol)),
        |(x, y, z, symbol)| {
            let (element, isotope_mass) = convert_atom_symbol_mass_diff(symbol, None);
            AtomStandard {
                element,
                charge: 0,
                radical: None,
                isotope_mass,
                stereo_parity: None,
                hydrogen_count: None,
                valence: None,
                atom_map_num: None,
                position: Some(Point3D::new(x, y, z)),
                properties: std::collections::HashMap::new(),
            }
        },
    )
    .parse(input)
}

/// Parse standard atom input (optimized for performance)
/// Fails immediately on non-standard atom symbols. For parsing all atom types, see atom_like_input.
///
/// xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee (69 characters wide)
///
/// *Values in the atom block*
/// -------------------------------------------------------------------------------
/// | Field | Meaning            | Values       | Notes                           |
/// |-------|--------------------|--------------|----------------------------------
/// | x,y,z | atom coordinates   |              | Generic, F10.4 format           |
/// | aaa   | atom symbol        | s. above     | Generic, Query, 3D, RGroup      |
/// | dd    | mass difference    | -3..=4       | Generic, s. also M  ISO         |
/// | ccc   | charge code        | 0..=7        | Generic, s. also M  CHG/M  RAD  |
/// | sss   | stereo parity      | 0..=3        | Generic, used in practice       |
/// | vvv   | valence code       | 0..=15       | Generic                         |
/// | mmm   | atom mapping       | 1..=#atoms   | Reaction, accepted as extension |
/// -------------------------------------------------------------------------------
///
pub fn atom_input_standard<'a>(
) -> impl Parser<&'a [u8], Output = AtomStandard, Error = error::Error<&'a [u8]>> {
    move |input: &'a [u8]| {
        let len = input.len();
        let parser = match len {
            67.. => atom_input_standard69,
            49..=66 => atom_input_standard51,
            42..=48 => atom_input_standard42,
            37..=41 => atom_input_standard39,
            35..=36 => atom_input_standard36,
            34 => atom_input_standard34,
            _ => {
                return Err(nom::Err::Error(error::Error::new(
                    input,
                    error::ErrorKind::Eof,
                )))
            }
        };
        terminated(parser, space0).parse(input)
    }
}

/// Parse atom and atom-like input
/// Allows all atom types. For faster parsing of standard molecules, see atom_input.
///
/// xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee (69 characters wide)
///
/// *Values in the atom block*
/// -----------------------------------------------------------------------------------------
/// | Field | Meaning            | Values       | Notes                                     |
/// |-------|--------------------|--------------|-------------------------------------------|
/// | x,y,z | atom coordinates   |              | Generic, F10.4 format                     |
/// | aaa   | atom symbol        | s. above     | Generic, Query, 3D, RGroup                |
/// | dd    | mass difference    | -3..=4       | Generic, s. also M  ISO                   |
/// | ccc   | charge code        | 0..=7        | Generic, s. also M  CHG/M  RAD            |
/// | sss   | stereo parity      | 0..=3        | Generic, ignored when read                |
/// | hhh   | hydrogen code      | 0..=5        | Query, H0 means no implicit Hs            |
/// |       |                    |              | Hn means >=n implicit Hs                  |
/// | bbb   | stereo care        | 0, 1         | Query, consider double bond stereo        |
/// |       |                    |              | when stereo care is 1 for both bond atoms |
/// | vvv   | valence code       | 0..=15       | Generic                                   |
/// | mmm   | atom mapping       | 1..=#atoms   | Reaction, accepted as extension           |
/// | nnn   | inversion          | 0..=2        | Reaction                                  |
/// | eee   | exact change       | 0, 1         | Reaction                                  |
/// -----------------------------------------------------------------------------------------
///
fn atom_input_inner(input: &[u8]) -> IResult<&[u8], Atom, error::Error<&[u8]>> {
    let (i, x) = fixed_width_float::<f64>(10, 4).parse(input)?;
    let (i, y) = fixed_width_float::<f64>(10, 4).parse(i)?;
    let (i, z) = fixed_width_float::<f64>(10, 4).parse(i)?;
    let (i, _) = take(1usize).parse(i)?;
    let (i, atom_symbol) = atom_symbol().parse(i)?;
    let (i, mass_diff) = map(fixed_width_int(2), convert_atom_mass_diff_code).parse(i)?;
    let (i, (charge, radical)) = map(fixed_width_int(3), convert_atom_charge_code).parse(i)?;

    let (i, stereo_parity) = cond(
        i.len() >= 3,
        map_res(fixed_width_int(3), convert_atom_stereo_parity_code),
    )
    .parse(i)?;

    let (i, hydrogen_count) = cond(
        i.len() >= 3,
        map_res(fixed_width_int(3), convert_atom_hydrogen_count_code),
    )
    .parse(i)?;

    let (i, stereo_care) = cond(
        i.len() >= 3,
        map_res(fixed_width_int(3), convert_atom_stereo_care_code),
    )
    .parse(i)?;

    let (i, valence) = cond(
        i.len() >= 3,
        map_res(fixed_width_int(3), convert_atom_valence_code),
    )
    .parse(i)?;

    // Skip unused HHH, rrr, iii fields
    let (i, _) = cond(i.len() >= 9, take(9usize)).parse(i)?;

    let (i, atom_map_num) = cond(i.len() >= 3, fixed_width_int::<u32>(3)).parse(i)?;

    let (i, inversion_flag) = cond(
        i.len() >= 3,
        map_res(fixed_width_int(3), convert_atom_inversion_flag_code),
    )
    .parse(i)?;

    let (i, exact_change_flag) = cond(
        i.len() >= 3,
        map_res(fixed_width_int(3), convert_atom_exact_change_flag_code),
    )
    .parse(i)?;

    let isotope_mass = match &atom_symbol {
        AtomSymbol::Element(e) => {
            mass_diff.map(|diff| (e.reference_mass_number() as i8 + diff) as u32)
        }
        AtomSymbol::NamedIsotope(i) => Some(i.mass_number()),
        _ => mass_diff.map(|diff| diff.unsigned_abs() as u32),
    };

    Ok((
        i,
        Atom {
            symbol: atom_symbol,
            charge,
            radical,
            isotope_mass,
            stereo_parity: stereo_parity.flatten(),
            hydrogen_count: hydrogen_count.flatten(),
            stereo_care: stereo_care.flatten(),
            valence: valence.flatten(),
            atom_map_num,
            inversion_retention: inversion_flag.flatten(),
            exact_change: exact_change_flag.flatten(),
            attachment_point: None,
            attachment_order: None,
            ring_bond_count: None,
            substitution_count: None,
            unsaturated: None,
            link_atom: None,
            position: Some(Point3D::new(x, y, z)),
            properties: std::collections::HashMap::new(),
        },
    ))
}

pub fn atom_input<'a>() -> impl Parser<&'a [u8], Output = Atom, Error = error::Error<&'a [u8]>> {
    terminated(atom_input_inner, space0)
}

#[cfg(test)]
mod tests;
