//! Atom block parser for CTab files.

use nom::branch::alt;
use nom::bytes::complete::{tag, take};
use nom::character::complete::{alpha1, multispace0, space0};
use nom::combinator::{all_consuming, complete, map, map_parser, map_res, value};
use nom::sequence::{delimited, preceded, terminated, tuple};
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
pub(crate) fn atom_input69<'a>(
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
//"    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0"
//"                                      ^         vvv"
//"                                      3           ^         mmm"
//"                                      9           5           ^"
//"                                                  1           6     ^"
//"                                                              3     6"
//"                                                                    9"
    eprintln!(
        "Calling <atom_input69> with input: {:?}",
        String::from_utf8_lossy(input)
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
            delimited(take(9usize), complete(atom_map_num), take(6usize)),
        ),
        |(x, y, z, symbol, mass_diff, charge_radical, valence, atom_map_num)| {
            let (element, isotope_mass) = convert_atom_symbol_mass_diff(symbol, mass_diff);
            let (charge, radical) = charge_radical;
            let atom = Atom {
                element,
                charge,
                isotope_mass,
                stereo_parity: None,
                hydrogen_count: None,
                valence,
                atom_map_num,
                radical,
                properties: std::collections::HashMap::new(),
            };

            (atom, Point3D::new(x, y, z))
        },
    )
    .parse(input)
}

/// Parse standard atom inputs with 49-51 characters (s. `atom_input` for more details).
/// Includes mass difference, charge/radical and valence fields.
pub(crate) fn atom_input51<'a>(
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

    eprintln!(
        "Calling <atom_input51> with input: {:?}",
        String::from_utf8_lossy(input)
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
            let (charge, radical) = charge_radical;
            let atom = Atom {
                element,
                charge,
                isotope_mass,
                stereo_parity: None,
                hydrogen_count: None,
                valence,
                atom_map_num: None,
                radical,
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
pub(crate) fn atom_input39<'a>(
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

    eprintln!(
        "Calling <atom_input39> with input: {:?}",
        String::from_utf8_lossy(input)
    );

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
    let (charge, radical) = charge_radical;
    let atom = Atom {
        element,
        charge,
        isotope_mass,
        stereo_parity: None,
        hydrogen_count: None,
        valence: None,
        atom_map_num: None,
        radical,
        properties: std::collections::HashMap::new(),
    };

    Ok((remaining, (atom, Point3D::new(x, y, z))))
}

/// Parse standard atom inputs with 35-36 characters (s. `atom_input` for more details).
/// Lacks trailing charge/radical, valence and atom mapping fields (substituted by defaults).
pub(crate) fn atom_input36<'a>(
    input: &'a [u8],
) -> IResult<&'a [u8], (Atom, Point3D), error::Error<&'a [u8]>> {
    let x = fixed_width_float::<f64>(10, 4);
    let y = fixed_width_float::<f64>(10, 4);
    let z = fixed_width_float::<f64>(10, 4);
    let symbol = atom_symbol_standard();
    let mass_diff = map(fixed_width_int_in_range_opt::<i8, _>(2, -3..=4), |opt| {
        convert_atom_mass_diff_code(opt.unwrap_or(0))
    });

    eprintln!(
        "Calling <atom_input36> with input: {:?}",
        String::from_utf8_lossy(input)
    );

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
pub(crate) fn atom_input34<'a>(
    input: &'a [u8],
) -> IResult<&'a [u8], (Atom, Point3D), error::Error<&'a [u8]>> {
    let x = fixed_width_float::<f64>(10, 4);
    let y = fixed_width_float::<f64>(10, 4);
    let z = fixed_width_float::<f64>(10, 4);
    let symbol = atom_symbol_standard();

    eprintln!(
        "Calling <atom_input34> with input: {:?}",
        String::from_utf8_lossy(input)
    );

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
pub(crate) fn atom_input<'a>(
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
pub(crate) fn atom_like_input<'a>(
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

                let (charge, radical) = charge_radical;

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
mod tests {
    use super::*;
    use float_cmp::approx_eq;
    use nom::{error::ErrorKind, Err};
    use rstest::rstest;

    #[rstest]
    #[case(b"A  ", AtomSymbol::Unspecified('A'))]
    #[case(b"Q  ", AtomSymbol::Unspecified('Q'))]
    #[case(b"*  ", AtomSymbol::Unspecified('*'))]
    #[case(b"L  ", AtomSymbol::AtomList(AtomList { elements: vec![] }))]
    #[case(b"LP ", AtomSymbol::LonePair)]
    #[case(b"R1 ", AtomSymbol::RGroup(0))]
    #[case(b"R3 ", AtomSymbol::RGroup(2))]
    #[case(b"H  ", AtomSymbol::Element(Element::H))]
    #[case(b"C  ", AtomSymbol::Element(Element::C))]
    #[case(b" C ", AtomSymbol::Element(Element::C))]
    #[case(b"  C", AtomSymbol::Element(Element::C))]
    #[case(b"Cu ", AtomSymbol::Element(Element::Cu))]
    #[case(b"cu ", AtomSymbol::Element(Element::Cu))]
    #[case(b"CU ", AtomSymbol::Element(Element::Cu))]
    #[case(b"D  ", AtomSymbol::NamedIsotope(NamedIsotope::D))]
    #[case(b"d  ", AtomSymbol::NamedIsotope(NamedIsotope::D))]
    #[case(b"T  ", AtomSymbol::NamedIsotope(NamedIsotope::T))]
    fn test_atom_symbol(#[case] input: &[u8], #[case] expected: AtomSymbol) {
        let (remaining, symbol) = atom_symbol().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(symbol, expected);
    }

    #[rstest]
    #[case(b"   ", "empty field", ErrorKind::Alpha)]
    #[case(b"H", "too short", ErrorKind::Eof)]
    #[case(b"R  ", "R group index missing", ErrorKind::MapRes)]
    #[case(b"R0 ", "R group index must be between 1 and 31", ErrorKind::MapRes)]
    #[case(b"R32", "R group index must be between 1 and 31", ErrorKind::MapRes)]
    #[case(b"Xx ", "Unknown atom symbol", ErrorKind::MapRes)]
    #[case(b"LQ ", "Unknown atom symbol", ErrorKind::Eof)]
    fn test_atom_symbol_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let result = atom_symbol().parse(input);
        assert!(result.is_err(), "{} should have failed", desc);
        assert!(
            matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "Mismatched error kind for {}, expected {:?}, got {}",
            desc,
            expected_kind,
            result.clone().unwrap_err().map(|e| e.code),
        );
    }

    #[rstest]
    #[case(b"H  ", AtomSymbol::Element(Element::H))]
    #[case(b"C  ", AtomSymbol::Element(Element::C))]
    #[case(b"Cu ", AtomSymbol::Element(Element::Cu))]
    #[case(b"D  ", AtomSymbol::NamedIsotope(NamedIsotope::D))]
    fn test_atom_symbol_standard(#[case] input: &[u8], #[case] expected: AtomSymbol) {
        let (remaining, symbol) = atom_symbol_standard().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(symbol, expected);
    }

    #[rstest]
    #[case(b"A  ", "unspecified atom", ErrorKind::MapRes)]
    #[case(b"L  ", "atom list", ErrorKind::MapRes)]
    #[case(b"LP ", "lone pair", ErrorKind::MapRes)]
    #[case(b"R1 ", "R group", ErrorKind::MapRes)]
    fn test_atom_symbol_standard_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let result = atom_symbol_standard().parse(input);
        assert!(
            result.is_err(),
            "{} should be rejected by standard parser",
            desc
        );
        assert!(
            matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "Mismatched error kind for {}, expected {:?}, got {}",
            desc,
            expected_kind,
            result.clone().unwrap_err().map(|e| e.code),
        );
    }

    #[rstest]
    #[case(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0  0  0", "standard valid", Element::C, Some(10), 1, Some(4), None)]
    #[case(b"    1.2345    2.3456    3.4567 C  -3  3  0  0  0  4  0  0  0  0  0  0", "mass diff lower bound", Element::C, Some(9), 1, Some(4), None)]
    #[case(b"    1.2345    2.3456    3.4567 C   4  3  0  0  0  4  0  0  0  0  0  0", "mass diff upper bound", Element::C, Some(16), 1, Some(4), None)]
    #[case(b"    1.2345    2.3456    3.4567 C  -4  3  0  0  0  4  0  0  0  0  0  0", "mass diff out-of-range low", Element::C, None, 1, Some(4), None)]
    #[case(b"    1.2345    2.3456    3.4567 C   5  3  0  0  0  4  0  0  0  0  0  0", "mass diff out-of-range high", Element::C, None, 1, Some(4), None)]
    #[case(b"    1.2345    2.3456    3.4567 C  -2  8  0  0  0  4  0  0  0  0  0  0", "charge out-of-range high", Element::C, Some(10), 0, Some(4), None)]
    #[case(b"    1.2345    2.3456    3.4567 C  -2  0  0  0  0  4  0  0  0  1  0  0", "atom map num non-zero", Element::C, Some(10), 0, Some(4), Some(1))]
    #[case(b"    1.2345    2.3456    3.4567 C  -2  3  0  a  0  4  0  0  0  0  0  0", "ignore block 1", Element::C, Some(10), 1, Some(4), None)]
    #[case(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  a  0  0  0  0  0", "ignore block 2", Element::C, Some(10), 1, Some(4), None)]
    #[case(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0  a  0", "ignore block 3", Element::C, Some(10), 1, Some(4), None)]
    fn test_atom_input69(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_element: Element,
        #[case] expected_isotope_mass: Option<u32>,
        #[case] expected_charge: i8,
        #[case] expected_valence: Option<u8>,
        #[case] expected_atom_map_num: Option<u32>,
    ) {
        let result = atom_input69(input);
        assert!(
            result.is_ok(),
            "Parser failed for case '{}': {:?}",
            desc,
            result
        );
        let (remaining, (atom, pos)) = result.unwrap();
        assert!(
            remaining.is_empty(),
            "Non-empty input for case '{}'",
            desc
        );
        assert_eq!(
            atom.element, expected_element,
            "Mismatched element for '{}'",
            desc
        );
        assert_eq!(atom.isotope_mass, expected_isotope_mass, "Mismatched isotope mass for '{}'", desc);
        assert_eq!(atom.charge, expected_charge, "Mismatched charge for '{}'", desc);
        assert_eq!(atom.valence, expected_valence, "Mismatched valence for '{}'", desc);
        assert_eq!(atom.atom_map_num, expected_atom_map_num, "Mismatched atom map num for '{}'", desc);
    }

    #[rstest]
    #[case(b"    1.234a    2.3456    3.4567 C   0  0  0  0  0  0  0  0  0  0  0  0", "non-numeric coordinate", ErrorKind::Eof)]
    #[case(b"    1.2345    2.3456    3.4567 C   0  0  0  0  0  0  0  0  0  a  0  0", "non-numeric atom map number", ErrorKind::Digit)]
    #[case(b"    1.2345    2.3456    3.4567 L   0  0  0  0  0  0  0  0  0  0  0  0", "non-standard atom symbol", ErrorKind::MapRes)]
    fn test_atom_input69_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let result = atom_input69(input);
        assert!(result.is_err(), "Parser should have failed for '{}'", desc);
        let err = result.unwrap_err();
        if let Err::Error(e) = err {
            assert_eq!(
                e.code, expected_kind,
                "Mismatched error kind for '{}'",
                desc
            );
        }
    }

    #[rstest]
    #[case(
        b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4",
        "standard valid",
        Element::C,
        Some(10),
        1,
        Some(4)
    )]
    #[case(
        b"    1.2345    2.3456    3.4567 C  -3  3  0  0  0  4",
        "mass diff lower bound",
        Element::C,
        Some(9),
        1,
        Some(4)
    )]
    #[case(
        b"    1.2345    2.3456    3.4567 C   4  3  0  0  0  4",
        "mass diff upper bound",
        Element::C,
        Some(16),
        1,
        Some(4)
    )]
    #[case(
        b"    1.2345    2.3456    3.4567 C  -4  3  0  0  0  4",
        "mass diff out-of-range low",
        Element::C,
        None,
        1,
        Some(4)
    )]
    #[case(
        b"    1.2345    2.3456    3.4567 C   5  3  0  0  0  4",
        "mass diff out-of-range high",
        Element::C,
        None,
        1,
        Some(4)
    )]
    #[case(
        b"    1.2345    2.3456    3.4567 C  -2  8  0  0  0  4",
        "charge out-of-range high",
        Element::C,
        Some(10),
        0,
        Some(4)
    )]
    fn test_atom_input51(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_element: Element,
        #[case] expected_isotope_mass: Option<u32>,
        #[case] expected_charge: i8,
        #[case] expected_valence: Option<u8>,
    ) {
        let result = atom_input51(input);
        assert!(
            result.is_ok(),
            "Parser failed for case '{}': {:?}",
            desc,
            result
        );
        let (remaining, (atom, pos)) = result.unwrap();
        assert!(
            remaining.is_empty(),
            "Remaining non-empty for case '{}': '{}'",
            desc,
            String::from_utf8_lossy(remaining)
        );
        assert_eq!(
            atom.element, expected_element,
            "Mismatched element for '{}'",
            desc
        );
        assert!(
            approx_eq!(f64, pos.x, 1.2345),
            "Mismatched x for '{}'",
            desc
        );
        assert!(
            approx_eq!(f64, pos.y, 2.3456),
            "Mismatched y for '{}'",
            desc
        );
        assert!(
            approx_eq!(f64, pos.z, 3.4567),
            "Mismatched z for '{}'",
            desc
        );
        assert_eq!(
            atom.isotope_mass, expected_isotope_mass,
            "Mismatched isotope mass for '{}'",
            desc
        );
        assert_eq!(
            atom.charge, expected_charge,
            "Mismatched charge for '{}'",
            desc
        );
        assert_eq!(
            atom.valence, expected_valence,
            "Mismatched valence for '{}'",
            desc
        );
    }

    #[rstest]
    #[case(
        b"    1.234a    2.3456    3.4567 C  -2  3  0  0  0  4",
        "non-numeric coordinate",
        ErrorKind::Eof
    )]
    #[case(
        b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  a",
        "non-numeric valence",
        ErrorKind::Digit
    )]
    #[case(
        b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0 16",
        "out-of-range valence",
        ErrorKind::Verify
    )]
    #[case(
        b"    1.2345    2.3456    3.4567 L  -2  3  0  0  0  4",
        "invalid atom symbol",
        ErrorKind::MapRes
    )]
    fn test_atom_input51_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let result = atom_input51(input);
        assert!(result.is_err(), "Parser should have failed for '{}'", desc);
        let err = result.unwrap_err();
        assert!(
            matches!(err, Err::Error(_)),
            "Error should be a nom::Err::Error for '{}'",
            desc
        );
        if let Err::Error(e) = err {
            assert_eq!(
                e.code, expected_kind,
                "Mismatched error kind for '{}'",
                desc
            );
        }
    }

    #[rstest]
    #[case(
        b"    1.2345    2.3456    3.4567 C      3  0  0  0  4",
        "blank mass diff",
        None,
        1,
        Some(4)
    )]
    #[case(
        b"    1.2345    2.3456    3.4567 C  -2     0  0  0  4",
        "blank charge",
        Some(10),
        0,
        Some(4)
    )]
    #[case(
        b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0   ",
        "blank valence",
        Some(10),
        1,
        None
    )]
    fn test_atom_input51_empty_fields(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_isotope_mass: Option<u32>,
        #[case] expected_charge: i8,
        #[case] expected_valence: Option<u8>,
    ) {
        let result = atom_input51(input);
        assert!(
            result.is_ok(),
            "Parser failed for case '{}': {:?}",
            desc,
            result
        );
        let (remaining, (atom, _)) = result.unwrap();
        assert!(
            remaining.is_empty(),
            "Remaining non-empty for case '{}': '{}'",
            desc,
            String::from_utf8_lossy(remaining)
        );
        assert_eq!(
            atom.isotope_mass, expected_isotope_mass,
            "Mismatched isotope mass for '{}'",
            desc
        );
        assert_eq!(
            atom.charge, expected_charge,
            "Mismatched charge for '{}'",
            desc
        );
        assert_eq!(
            atom.valence, expected_valence,
            "Mismatched valence for '{}'",
            desc
        );
    }

    #[rstest]
    #[case(
        b"    1.2345    2.3456    3.4567 C  -2  3",
        "standard valid",
        Element::C,
        Some(10),
        1
    )]
    #[case(
        b"    1.2345    2.3456    3.4567 C  -4  3",
        "mass diff out-of-range low",
        Element::C,
        None,
        1
    )]
    #[case(
        b"    1.2345    2.3456    3.4567 C  -2  8",
        "charge out-of-range high",
        Element::C,
        Some(10),
        0
    )]
    fn test_atom_input39(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_element: Element,
        #[case] expected_isotope_mass: Option<u32>,
        #[case] expected_charge: i8,
    ) {
        let result = atom_input39(input);
        assert!(
            result.is_ok(),
            "Parser failed for case '{}': {:?}",
            desc,
            result
        );
        let (remaining, (atom, pos)) = result.unwrap();
        assert!(
            remaining.is_empty(),
            "Remaining non-empty for case '{}': '{}'",
            desc,
            String::from_utf8_lossy(remaining)
        );
        assert_eq!(
            atom.element, expected_element,
            "Mismatched element for '{}'",
            desc
        );
        assert_eq!(
            atom.isotope_mass, expected_isotope_mass,
            "Mismatched isotope mass for '{}'",
            desc
        );
        assert_eq!(
            atom.charge, expected_charge,
            "Mismatched charge for '{}'",
            desc
        );
        assert_eq!(atom.valence, None, "Valence is not None for '{}'", desc);
        assert!(
            approx_eq!(f64, pos.x, 1.2345),
            "Mismatched x for '{}'",
            desc
        );
    }

    #[rstest]
    #[case(
        b"    1.234a    2.3456    3.4567 C  -2  3",
        "non-numeric coordinate",
        ErrorKind::Eof
    )]
    #[case(
        b"    1.2345    2.3456    3.4567 C  -a  3",
        "non-numeric mass diff",
        ErrorKind::Digit
    )]
    #[case(
        b"    1.2345    2.3456    3.4567 L  -2  3",
        "invalid atom symbol",
        ErrorKind::MapRes
    )]
    fn test_atom_input39_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let result = atom_input39(input);
        assert!(result.is_err(), "Parser should have failed for '{}'", desc);
        if let Err(Err::Error(e)) = result {
            assert_eq!(
                e.code, expected_kind,
                "Mismatched error kind for '{}'",
                desc
            );
        } else {
            panic!(
                "Expected a nom::Err::Error for '{}', got {:?}",
                desc, result
            );
        }
    }

    #[rstest]
    #[case(b"    1.2345    2.3456    3.4567 C      3", "blank mass diff", None, 1)]
    #[case(
        b"    1.2345    2.3456    3.4567 C  -2   ",
        "blank charge",
        Some(10),
        0
    )]
    fn test_atom_input39_empty_fields(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_isotope_mass: Option<u32>,
        #[case] expected_charge: i8,
    ) {
        let result = atom_input39(input);
        assert!(
            result.is_ok(),
            "Parser failed for case '{}': {:?}",
            desc,
            result
        );
        let (_, (atom, _)) = result.unwrap();
        assert_eq!(
            atom.isotope_mass, expected_isotope_mass,
            "Mismatched isotope mass for '{}'",
            desc
        );
        assert_eq!(
            atom.charge, expected_charge,
            "Mismatched charge for '{}'",
            desc
        );
    }

    #[rstest]
    #[case(b"    1.2345    2.3456    3.4567 C  -2  3\n", "trailing newline")]
    #[case(b"    1.2345    2.3456    3.4567 C  -2  3   ", "trailing spaces")]
    #[case(b"    1.2345    2.3456    3.4567 C  -2  3\t\t", "trailing tabs")]
    fn test_atom_input39_whitespace_padded(#[case] input: &[u8], #[case] desc: &str) {
        let result = all_consuming(terminated(atom_input39, multispace0)).parse(input);
        assert!(
            result.is_ok(),
            "Parser failed for case '{}': {:?}",
            desc,
            result
        );
        let (remaining, (atom, _)) = result.unwrap();
        assert!(
            remaining.is_empty(),
            "Remaining non-empty for case '{}': '{}'",
            desc,
            String::from_utf8_lossy(remaining)
        );
        assert_eq!(atom.charge, 1);
    }

    #[rstest]
    #[case(
        b"    1.2345    2.3456    3.4567 C  -2  3abc",
        "valid numeric data in ignored gap"
    )]
    #[case(
        b"    1.2345    2.3456    3.4567 C  -2  3         ",
        "whitespace in ignored gap"
    )]
    fn test_atom_input39_ignored_gap(#[case] input: &[u8], #[case] desc: &str) {
        let result = atom_input39(input);
        assert!(
            result.is_ok(),
            "Parser failed for case '{}': {:?}",
            desc,
            result
        );
        let (remaining, (atom, _)) = result.unwrap();
        assert!(
            remaining.is_empty(),
            "Remaining input should be empty for case '{}', but was: '{}'",
            desc,
            String::from_utf8_lossy(remaining)
        );
        assert_eq!(atom.charge, 1);
    }

    #[rstest]
    #[case(
        b"    1.2345    2.3456    3.4567 C  -2",
        "standard valid",
        Element::C,
        Some(10)
    )]
    #[case(
        b"    1.2345    2.3456    3.4567 C  -4",
        "mass diff out-of-range low",
        Element::C,
        None
    )]
    fn test_atom_input36(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_element: Element,
        #[case] expected_isotope_mass: Option<u32>,
    ) {
        let result = atom_input36(input);
        assert!(
            result.is_ok(),
            "Parser failed for case '{}': {:?}",
            desc,
            result
        );
        let (remaining, (atom, _)) = result.unwrap();
        assert!(
            remaining.is_empty(),
            "Remaining non-empty for case '{}': '{}'",
            desc,
            String::from_utf8_lossy(remaining)
        );
        assert_eq!(
            atom.element, expected_element,
            "Mismatched element for '{}'",
            desc
        );
        assert_eq!(
            atom.isotope_mass, expected_isotope_mass,
            "Mismatched isotope mass for '{}'",
            desc
        );
        assert_eq!(atom.charge, 0, "Charge should be 0 for '{}'", desc);
        assert_eq!(atom.valence, None, "Valence should be None for '{}'", desc);
    }

    #[rstest]
    #[case(
        b"    1.234a    2.3456    3.4567 C  -2",
        "non-numeric coordinate",
        ErrorKind::Eof
    )]
    #[case(
        b"    1.2345    2.3456    3.4567 L  -2",
        "invalid atom symbol",
        ErrorKind::MapRes
    )]
    fn test_atom_input36_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let result = atom_input36(input);
        assert!(result.is_err(), "Parser should have failed for '{}'", desc);
        if let Err(Err::Error(e)) = result {
            assert_eq!(
                e.code, expected_kind,
                "Mismatched error kind for '{}'",
                desc
            );
        } else {
            panic!(
                "Expected a nom::Err::Error for '{}', got {:?}",
                desc, result
            );
        }
    }

    #[rstest]
    #[case(b"    1.2345    2.3456    3.4567 C    ", "blank mass diff")]
    fn test_atom_input36_empty_fields(#[case] input: &[u8], #[case] desc: &str) {
        let result = atom_input36(input);
        assert!(
            result.is_ok(),
            "Parser failed for case '{}': {:?}",
            desc,
            result
        );
        let (_, (atom, _)) = result.unwrap();
        assert_eq!(
            atom.isotope_mass, None,
            "Mismatched isotope mass for '{}'",
            desc
        );
    }

    #[rstest]
    #[case(
        b"    1.2345    2.3456    3.4567 C  -2  \n",
        "trailing whitespace and newline"
    )]
    fn test_atom_input36_whitespace_padded(#[case] input: &[u8], #[case] _desc: &str) {
        let result = all_consuming(terminated(atom_input36, multispace0)).parse(input);
        assert!(
            result.is_ok(),
            "Parser failed for case '{}': {:?}",
            _desc,
            result
        );
        let (remaining, (atom, _)) = result.unwrap();
        assert!(
            remaining.is_empty(),
            "Remaining input should be empty for case '{}'",
            _desc
        );
        assert_eq!(atom.isotope_mass, Some(10));
    }

    #[rstest]
    #[case(b"    1.2345    2.3456    3.4567 C  ", "standard valid", Element::C)]
    fn test_atom_input34(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_element: Element,
    ) {
        let result = atom_input34(input);
        assert!(
            result.is_ok(),
            "Parser failed for case '{}': {:?}",
            desc,
            result
        );
        let (remaining, (atom, _)) = result.unwrap();
        assert!(
            remaining.is_empty(),
            "Remaining non-empty for case '{}': '{}'",
            desc,
            String::from_utf8_lossy(remaining)
        );
        assert_eq!(
            atom.element, expected_element,
            "Mismatched element for '{}'",
            desc
        );
        assert_eq!(
            atom.isotope_mass, None,
            "Isotope mass should be None for '{}'",
            desc
        );
        assert_eq!(atom.charge, 0, "Charge should be 0 for '{}'", desc);
        assert_eq!(atom.valence, None, "Valence should be None for '{}'", desc);
    }

    #[rstest]
    #[case(
        b"    1.234a    2.3456    3.4567 C  ",
        "non-numeric coordinate",
        ErrorKind::Eof
    )]
    #[case(
        b"    1.2345    2.3456    3.4567 L  ",
        "invalid atom symbol",
        ErrorKind::MapRes
    )]
    fn test_atom_input34_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let result = atom_input34(input);
        assert!(result.is_err(), "Parser should have failed for '{}'", desc);
        if let Err(Err::Error(e)) = result {
            assert_eq!(
                e.code, expected_kind,
                "Mismatched error kind for '{}'",
                desc
            );
        } else {
            panic!(
                "Expected a nom::Err::Error for '{}', got {:?}",
                desc, result
            );
        }
    }

    #[rstest]
    #[case(
        b"    1.2345    2.3456    3.4567 C    \n",
        "trailing whitespace and newline"
    )]
    fn test_atom_input34_whitespace_padded(#[case] input: &[u8], #[case] _desc: &str) {
        let result = all_consuming(terminated(atom_input34, multispace0)).parse(input);
        assert!(
            result.is_ok(),
            "Parser failed for case '{}': {:?}",
            _desc,
            result
        );
        let (remaining, (atom, _)) = result.unwrap();
        assert!(
            remaining.is_empty(),
            "Remaining input should be empty for case '{}'",
            _desc
        );
        assert_eq!(atom.element, Element::C);
    }

    #[rstest]
    #[case(
        b"    1.0000    2.0000    3.0000 C  ",
        "len 34",
        Element::C,
        None,
        0,
        None
    )]
    #[case(
        b"    1.0000    2.0000    3.0000 C   ",
        "len 35 padded",
        Element::C,
        None,
        0,
        None
    )]
    #[case(
        b"    1.0000    2.0000    3.0000 C  -2",
        "len 36",
        Element::C,
        Some(10),
        0,
        None
    )]
    #[case(
        b"    1.0000    2.0000    3.0000 C  -2  ",
        "len 38 padded",
        Element::C,
        Some(10),
        0,
        None
    )]
    #[case(
        b"    1.0000    2.0000    3.0000 C  -2  3",
        "len 39",
        Element::C,
        Some(10),
        1,
        None
    )]
    #[case(
        b"    1.0000    2.0000    3.0000 C  -2  3abc",
        "len 42 with gap",
        Element::C,
        Some(10),
        1,
        None
    )]
    #[case(
        b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4",
        "len 51",
        Element::C,
        Some(10),
        1,
        Some(4)
    )]
    #[case(
        b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4 ",
        "len 52 padded",
        Element::C,
        Some(10),
        1,
        Some(4)
    )]
    #[case(
        b"    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
        "len 69 zeros",
        Element::C,
        None,
        0,
        None
    )]
    fn test_atom_input(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_element: Element,
        #[case] expected_isotope_mass: Option<u32>,
        #[case] expected_charge: i8,
        #[case] expected_valence: Option<u8>,
    ) {
        let mut parser = atom_input();
        let result = parser.parse(input);
        assert!(
            result.is_ok(),
            "Parser failed for case '{}': {:?}",
            desc,
            result
        );
        let (remaining, (atom, _)) = result.unwrap();
        assert!(
            remaining.is_empty(),
            "Remaining non-empty for case '{}': '{}'",
            desc,
            String::from_utf8_lossy(remaining)
        );
        assert_eq!(
            atom.element, expected_element,
            "Mismatched element for '{}'",
            desc
        );
        assert_eq!(
            atom.isotope_mass, expected_isotope_mass,
            "Mismatched isotope mass for '{}'",
            desc
        );
        assert_eq!(
            atom.charge, expected_charge,
            "Mismatched charge for '{}'",
            desc
        );
        assert_eq!(
            atom.valence, expected_valence,
            "Mismatched valence for '{}'",
            desc
        );
    }

    #[rstest]
    #[case(
        b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  a",
        "non-numeric valence",
        ErrorKind::Digit
    )]
    #[case(
        b"    1.0000    2.0000    3.0000 L  -2  3  0  0  0  4",
        "invalid element",
        ErrorKind::MapRes
    )]
    fn test_atom_input_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let mut parser = atom_input();
        let result = parser.parse(input);
        assert!(result.is_err(), "Parser should have failed for '{}'", desc);
        if let Err(Err::Error(e)) = result {
            assert_eq!(
                e.code, expected_kind,
                "Mismatched error kind for '{}'",
                desc
            );
        } else {
            panic!(
                "Expected a nom::Err::Error for '{}', got {:?}",
                desc, result
            );
        }
    }

    #[test]
    fn test_atom_input_partial_fields() {
        let input = b"    1.0000    2.0000    3.0000 C  -2 3"; // len 38
        let mut parser = atom_input();
        let result = parser.parse(input);
        assert!(
            result.is_err(),
            "Parser should have failed for partial field"
        );
        if let Err(Err::Error(e)) = result {
            // The charge field is 3 chars, we provided 2, fixed_width_opt should see a partial non-whitespace field and fail with Eof
            assert_eq!(
                e.code,
                ErrorKind::Eof,
                "Mismatched error kind for partial field"
            );
        } else {
            panic!(
                "Expected a nom::Err::Error for partial field, got {:?}",
                result
            );
        }
    }

    #[rstest]
    #[case(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4   \t", "len 55")]
    #[case(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0  0  0  0           ", "len 80")]
    fn test_atom_input_whitespace_padded(#[case] input: &[u8], #[case] desc: &str) {
        let mut parser = atom_input();
        let result = parser.parse(input);
        assert!(
            result.is_ok(),
            "Parser failed for whitespace padded input: {:?}",
            result
        );
        let (remaining, (atom, _)) = result.unwrap();
        assert!(
            remaining.is_empty(),
            "Non-empty input for case '{}'",
            desc
        );
        assert_eq!(atom.charge, 1);
        assert_eq!(atom.valence, Some(4));
    }
}
