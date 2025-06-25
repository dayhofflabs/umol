//! Atom block parser for CTab files.

use nom::branch::alt;
use nom::bytes::complete::{tag, take};
use nom::character::complete::{alpha1, multispace0, space0};
use nom::combinator::{all_consuming, complete, map, map_parser, map_res, value, verify};
use nom::sequence::{delimited, preceded, terminated};
use nom::{error, IResult, Parser};
use umol_data::{Element, NamedIsotope};

use crate::atom::{Atom, AtomLike, AtomList, AtomSymbol};
use crate::conformer::Point3D;

use super::convert::{
    convert_atom_charge_code, convert_atom_hydrogen_count_code, convert_atom_mass_diff_code,
    convert_atom_stereo_parity_code, convert_atom_symbol_mass_diff, convert_atom_valence_code,
};
use super::utils::{
    fixed_width_float, fixed_width_int, fixed_width_int_in_range, fixed_width_int_in_range_minus1,
    fixed_width_int_in_range_opt,
};

/// Parse atom symbol (standard atoms, Element and NamedIsotope only)
/// Returns error for non-standard atom symbols (L, A, Q, *, LP, R#).
/// See `atom_symbol` for more details.
fn atom_symbol_standard<'a>(
) -> impl Parser<&'a [u8], Output = AtomSymbol, Error = error::Error<&'a [u8]>> {
    |input| {
        map_parser(
            take(3usize),
            all_consuming(alt((
                map(
                    map_res(delimited(space0, alpha1, space0), |s| {
                        Element::from_symbol_bytes(s)
                            .ok_or_else(|| error::Error::new(s, error::ErrorKind::MapRes))
                    }),
                    |element| AtomSymbol::Element(element),
                ),
                map(
                    map_res(delimited(space0, alpha1, space0), |s| {
                        NamedIsotope::from_symbol_bytes(s)
                            .ok_or_else(|| error::Error::new(s, error::ErrorKind::MapRes))
                    }),
                    |isotope| AtomSymbol::NamedIsotope(isotope),
                ),
            ))),
        )
        .parse(input)
    }
}

/// Parse atom symbol (all atom types allowed in MOL specification).
/// -----------------------------------------------------------------------------
/// | Symbol   | Type          | Parser* | Notes                                |
/// -----------------------------------------------------------------------------
/// | H-Og     | Element       | S, A    |                                      |
/// | D, T     | Named Isotope | S, A    | Heavy H isotopes are an extension    |
/// | L        | Atom List     | A       | Query molecules                      |
/// | A, Q, *  | Unspecified   | A       | Query molecules, rarely in oligomers |
/// | LP       | Lone Pair     | A       | Rarely used                          |
/// | R#       | R Group       | A       | Query molecules                      |
/// -----------------------------------------------------------------------------
/// | *Parsers: S: standard, A: all                                             |
/// |----------------------------------------------------------------------------
///
fn atom_symbol<'a>() -> impl Parser<&'a [u8], Output = AtomSymbol, Error = error::Error<&'a [u8]>> {
    |input| {
        map_parser(
            take(3usize),
            all_consuming(alt((
                value(AtomSymbol::LonePair, delimited(space0, tag("LP"), space0)),
                value(
                    AtomSymbol::AtomList(AtomList { elements: vec![] }),
                    delimited(space0, tag("L"), space0),
                ),
                value(
                    AtomSymbol::Unspecified('A'),
                    delimited(space0, tag("A"), space0),
                ),
                value(
                    AtomSymbol::Unspecified('Q'),
                    delimited(space0, tag("Q"), space0),
                ),
                value(
                    AtomSymbol::Unspecified('*'),
                    delimited(space0, tag("*"), space0),
                ),
                map(
                    delimited(
                        space0,
                        (
                            tag("R"),
                            fixed_width_int_in_range_minus1::<usize, _>(2, 1..=31),
                        ),
                        space0,
                    ),
                    |(_, idx)| AtomSymbol::RGroup(idx),
                ),
                map(
                    map_res(delimited(space0, alpha1, space0), |s| {
                        Element::from_symbol_bytes(s)
                            .ok_or_else(|| error::Error::new(s, error::ErrorKind::MapRes))
                    }),
                    |element| AtomSymbol::Element(element),
                ),
                map(
                    map_res(delimited(space0, alpha1, space0), |s| {
                        NamedIsotope::from_symbol_bytes(s)
                            .ok_or_else(|| error::Error::new(s, error::ErrorKind::MapRes))
                    }),
                    |isotope| AtomSymbol::NamedIsotope(isotope),
                ),
            ))),
        )
        .parse(input)
    }
}

/// Parse standard atom inputs with 52-69 characters (s. `atom_input` for more details).
/// Includes atom mapping number. Ignores whitespace padding.
fn atom_input69<'a>(
    input: &'a [u8],
) -> IResult<&'a [u8], (Atom, Point3D), error::Error<&'a [u8]>> {
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
            preceded(take(9usize), complete(valence)),
            delimited(take(9usize), complete(atom_map_num), take(6usize)),
        ),
        |(x, y, z, symbol, mass_diff, charge_radical, valence, atom_map_num)| {
            let (element, isotope_mass) = convert_atom_symbol_mass_diff(symbol, mass_diff);
            let (charge, _radical) = charge_radical;
            let atom = Atom {
                element,
                charge,
                isotope_mass,
                stereo_parity: None,
                hydrogen_count: None,
                valence,
                atom_map_num,
                radical: None,
                properties: std::collections::HashMap::new(),
            };

            (atom, Point3D::new(x, y, z))
        },
    )
    .parse(input)
}

/// Parse standard atom inputs with 49-51 characters (s. `atom_input` for more details).
/// Includes mass difference, charge/radical and valence fields.
fn atom_input51<'a>(
    input: &'a [u8],
) -> IResult<&'a [u8], (Atom, Point3D), error::Error<&'a [u8]>> {
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

    map(
        (
            x,
            y,
            z,
            preceded(take(1usize), symbol),
            complete(mass_diff),
            complete(charge_radical),
            preceded(take(9usize), complete(valence)),
        ),
        |(x, y, z, symbol, mass_diff, charge_radical, valence)| {
            let (element, isotope_mass) = convert_atom_symbol_mass_diff(symbol, mass_diff);
            let (charge, _radical) = charge_radical;
            let atom = Atom {
                element,
                charge,
                isotope_mass,
                stereo_parity: None,
                hydrogen_count: None,
                valence,
                atom_map_num: None,
                radical: None,
                properties: std::collections::HashMap::new(),
            };

            (atom, Point3D::new(x, y, z))
        },
    )
    .parse(input)
}

/// Parse standard atom inputs with 37-48 characters, including up to 9 characters of ignored data
/// (s. `atom_input` for more details). Lacks trailing valence and atom mapping fields
/// (substituted by defaults).
fn atom_input39<'a>(
    input: &'a [u8],
) -> IResult<&'a [u8], (Atom, Point3D), error::Error<&'a [u8]>> {
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

    let (input, (x, y, z, symbol, mass_diff, charge_radical)) = (
        x,
        y,
        z,
        preceded(take(1usize), symbol),
        complete(mass_diff),
        complete(charge_radical),
    )
        .parse(input)?;

    // Consume and discard up to 9 bytes from the remainder (ignored block).
    let n_to_take = 9.min(input.len());
    let (remaining, _) = take(n_to_take)(input)?;

    let (element, isotope_mass) = convert_atom_symbol_mass_diff(symbol, mass_diff);
    let (charge, _radical) = charge_radical;
    let atom = Atom {
        element,
        charge,
        isotope_mass,
        stereo_parity: None,
        hydrogen_count: None,
        valence: None,
        atom_map_num: None,
        radical: None,
        properties: std::collections::HashMap::new(),
    };

    Ok((remaining, (atom, Point3D::new(x, y, z))))
}

/// Parse standard atom inputs with 35-36 characters (s. `atom_input` for more details).
/// Lacks trailing charge/radical, valence and atom mapping fields (substituted by defaults).
fn atom_input36<'a>(
    input: &'a [u8],
) -> IResult<&'a [u8], (Atom, Point3D), error::Error<&'a [u8]>> {
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
            let atom = Atom {
                element,
                charge: 0,
                isotope_mass,
                stereo_parity: None,
                hydrogen_count: None,
                valence: None,
                atom_map_num: None,
                radical: None,
                properties: std::collections::HashMap::new(),
            };

            (atom, Point3D::new(x, y, z))
        },
    )
    .parse(input)
}

/// Parse standard atom inputs with 34 characters (s. `atom_input` for more details).
/// Lacks trailing mass difference, charge/radical, valence and atom mapping fields
/// (substituted by defaults).
fn atom_input34<'a>(
    input: &'a [u8],
) -> IResult<&'a [u8], (Atom, Point3D), error::Error<&'a [u8]>> {
    let x = fixed_width_float::<f64>(10, 4);
    let y = fixed_width_float::<f64>(10, 4);
    let z = fixed_width_float::<f64>(10, 4);
    let symbol = atom_symbol_standard();

    map(
        (x, y, z, preceded(take(1usize), symbol)),
        |(x, y, z, symbol)| {
            let (element, isotope_mass) = convert_atom_symbol_mass_diff(symbol, None);
            let atom = Atom {
                element,
                charge: 0,
                isotope_mass,
                stereo_parity: None,
                hydrogen_count: None,
                valence: None,
                atom_map_num: None,
                radical: None,
                properties: std::collections::HashMap::new(),
            };

            (atom, Point3D::new(x, y, z))
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
/// ------------------------------------------------------------------------------
/// | Field | Meaning            | Values       | Notes                          |
/// |-------|--------------------|--------------|---------------------------------
/// | x,y,z | atom coordinates   |              | Generic, F10.4 format          |
/// | aaa   | atom symbol        | s. above     | Generic, Query, 3D, RGroup     |
/// | dd    | mass difference    | -3..=4       | Generic, s. also M  ISO        |
/// | ccc   | charge code        | 0..=7        | Generic, s. also M  CHG/M  RAD |
/// | vvv   | valence code       | 0..=15       | Generic                        |
/// | mmm   | atom mapping       | 1..=#atoms   | Reaction, accepted as extnsion |
/// ------------------------------------------------------------------------------
///
pub fn atom_input<'a>(
) -> impl Parser<&'a [u8], Output = (Atom, Point3D), Error = error::Error<&'a [u8]>> {
    move |input: &'a [u8]| {
        let len = input.len();
        let parser = match len {
            67.. => atom_input69,
            49..=66 => atom_input51,
            37..=48 => atom_input39,
            35..=36 => atom_input36,
            34 => atom_input34,
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
/// | aaa   | atom symbol        | s. above     | Generic, Query, 3D, RGroup]               |
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
pub fn atom_like_input<'a>(
) -> impl Parser<&'a [u8], Output = (AtomLike, Point3D), Error = error::Error<&'a [u8]>> {
    let x = fixed_width_float::<f64>(10, 4);
    let y = fixed_width_float::<f64>(10, 4);
    let z = fixed_width_float::<f64>(10, 4);
    let symbol = atom_symbol();
    let mass_diff = map(fixed_width_int_in_range_opt::<i8, _>(2, -3..=4), |opt| {
        convert_atom_mass_diff_code(opt.unwrap_or(0))
    });
    let charge = map(fixed_width_int_in_range_opt::<u8, _>(3, 0..=7), |opt| {
        convert_atom_charge_code(opt.unwrap_or(0))
    });
    let stereo_parity = map_res(
        fixed_width_int_in_range::<u8, _>(3, 0..=3),
        convert_atom_stereo_parity_code,
    );
    let hydrogen_count = map_res(
        fixed_width_int_in_range::<u8, _>(3, 0..=5),
        convert_atom_hydrogen_count_code,
    );
    let stereo_care = fixed_width_int_in_range::<u8, _>(3, 0..=1);
    let valence = map_res(
        fixed_width_int_in_range::<u8, _>(3, 0..=15),
        convert_atom_valence_code,
    );
    let ignore_gap = |input: &'a [u8]| take(9.min(input.len()))(input);
    let atom_mapping = fixed_width_int::<u8>(3);
    let inversion = fixed_width_int_in_range::<u8, _>(3, 0..=2);
    let exact_change = fixed_width_int_in_range::<u8, _>(3, 0..=1);

    terminated(
        map(
            (
                x,
                y,
                z,
                take(1usize),
                symbol,
                mass_diff,
                charge,
                stereo_parity,
                hydrogen_count,
                stereo_care,
                valence,
                ignore_gap,
                atom_mapping,
                inversion,
                exact_change,
            ),
            |(
                x,
                y,
                z,
                _,
                symbol,
                mass_diff,
                charge_radical,
                stereo_parity,
                hydrogen_count,
                _stereo_care,
                valence,
                _ignored_gap,
                atom_mapping,
                _inversion,
                _exact_change,
            )| {
                // Calculate isotope mass based on symbol type
                let isotope_mass = match &symbol {
                    AtomSymbol::Element(e) => mass_diff.and_then(|diff| {
                        if diff != 0 {
                            Some((e.reference_mass_number() as i8 + diff) as u32)
                        } else {
                            None
                        }
                    }),
                    AtomSymbol::NamedIsotope(i) => {
                        // For named isotopes, use the isotope's mass, ignore mass_diff
                        Some(i.mass_number())
                    }
                    _ => {
                        // For other atom types, mass_diff might still be meaningful
                        mass_diff.and_then(|diff| {
                            if diff != 0 {
                                Some(diff.abs() as u32)
                            } else {
                                None
                            }
                        })
                    }
                };

                let (charge, _radical) = charge_radical;

                let atom_like = AtomLike {
                    symbol,
                    charge: charge,
                    isotope_mass,
                    stereo_parity,
                    hydrogen_count,
                    valence,
                    atom_map_num: Some(atom_mapping as u32),
                };

                (atom_like, Point3D::new(x, y, z))
            },
        ),
        multispace0,
    )
}

#[cfg(test)]
mod tests;