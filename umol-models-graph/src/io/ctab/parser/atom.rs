//! Atom block parser for CTab files.

use bstr::ByteSlice;
use nom::bytes::complete::tag;
use nom::character::complete::space0;
use nom::combinator::{cond, map, map_res};
use nom::sequence::{preceded, terminated};
use nom::{error, Err, IResult, Parser};
use umol_data::{Element, NamedIsotope};

use super::convert::{
    convert_atom_charge_code, convert_atom_exact_change_flag_code,
    convert_atom_hydrogen_count_code, convert_atom_inversion_flag_code,
    convert_atom_mass_diff_code, convert_atom_stereo_care_code, convert_atom_stereo_parity_code,
    convert_atom_symbol_mass_diff, convert_atom_valence_code,
};
use super::utils::{
    fixed_width_float, fixed_width_int, fixed_width_int_in_range, fixed_width_int_in_range_opt,
    fixed_width_padding_n, to_string,
};
use crate::io::ctab::atom::{Atom, AtomLike, AtomList, AtomSymbol};
use crate::io::ctab::config::CtabParseFlags;
use crate::io::ctab::parser::utils::fixed_width_partial;
use crate::io::ctab::query::QueryAtom;
use crate::io::ctab::rgroup::RGroup;
use crate::position::Point3D;

/// Parse atom symbol (Element and NamedIsotope only).
///
/// Returns error for atom-like symbols (L, A, Q, *, LP, R#).
/// If `allow_named_isotopes` is true, allow named isotopes (D, T).
fn atom_symbol<'inp, 'fl>(
    flags: &'fl CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = AtomSymbol, Error = error::Error<&'inp [u8]>> {
    let allow_named_isotopes = flags.contains(CtabParseFlags::NAMED_ISOTOPES);
    move |input: &'inp [u8]| {
        fixed_width_partial(
            3,
            move |s: &'inp [u8]| {
                let s = s.trim_ascii();
                if let Some(element) = Element::from_symbol_bytes(s) {
                    Ok((&b""[..], AtomSymbol::Element(element)))
                } else if allow_named_isotopes {
                    if let Some(isotope) = NamedIsotope::from_symbol_bytes(s) {
                        Ok((&b""[..], AtomSymbol::NamedIsotope(isotope)))
                    } else {
                        Err(Err::Error(error::Error::new(s, error::ErrorKind::MapRes)))
                    }
                } else {
                    Err(Err::Error(error::Error::new(s, error::ErrorKind::MapRes)))
                }
            },
            true,
        )
        .parse(input)
        .and_then(|(remaining, symbol)| match symbol {
            Some(symbol) => Ok((remaining, symbol)),
            None => Err(Err::Error(error::Error::new(input, error::ErrorKind::Eof))),
        })
    }
}

/// Parse atom-like symbol (all atom types allowed in MOL specification).
///
/// --------------------------------------------------------------------------------
/// | Symbol      | Type          | Parser* | Notes                                |
/// --------------------------------------------------------------------------------
/// | H-Og        | Element       | B, A    |                                      |
/// | D, T        | Named Isotope | B, A    | Heavy H isotopes as extension        |
/// | L           | Atom List     | A       | Query molecules                      |
/// | *,A,Q,X,M   | Query Atom    | A       | Query molecules, rarely in oligomers |
/// | AH,QH,XH,MH | Query Atom    | A       | Query molecules, CXSMILES extension  |
/// | LP          | Lone Pair     | A       | Rarely used                          |
/// | R#          | R Group       | A       | Query molecules                      |
/// --------------------------------------------------------------------------------
/// | *Parsers: B: basic, A: all                                                   |
/// --------------------------------------------------------------------------------
///
/// If `allow_named_isotopes` is true, allow named isotopes (D, T).
/// If `allow_rgroups` is true, allow rgroups (R#).
/// If `allow_queries` is true, allow queries (A, Q, *, LP, R#).
/// If `allow_extended_queries` is true, allow extended queries (AH, QH, XH, MH).
/// If `allow_subatoms` is true, allow subatoms (LP).
///
fn atomlike_symbol<'inp, 'fl>(
    flags: &'fl CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = AtomSymbol, Error = error::Error<&'inp [u8]>> {
    let allow_named_isotopes = flags.contains(CtabParseFlags::NAMED_ISOTOPES);
    let allow_rgroups = flags.contains(CtabParseFlags::RGROUPS);
    let allow_queries = flags.contains(CtabParseFlags::QUERIES);
    let allow_extended_queries = flags.contains(CtabParseFlags::EXTENDED_QUERIES);
    let allow_electrons = flags.contains(CtabParseFlags::ELECTRONS);
    let allow_pseudoatoms = flags.contains(CtabParseFlags::PSEUDOATOMS);
    move |input: &'inp [u8]| {
        fixed_width_partial(
            3,
            move |s: &'inp [u8]| {
                let s = s.trim_ascii();
                if let Some(element) = Element::from_symbol_bytes(s) {
                    return Ok((&b""[..], AtomSymbol::Element(element)));
                }
                if allow_named_isotopes {
                    if let Some(isotope) = NamedIsotope::from_symbol_bytes(s) {
                        return Ok((&b""[..], AtomSymbol::NamedIsotope(isotope)));
                    }
                }
                if allow_rgroups {
                    if let Some(rgroup) = RGroup::from_symbol_bytes(s) {
                        return Ok((&b""[..], AtomSymbol::RGroup(rgroup)));
                    }
                }
                if allow_electrons && s == b"LP" {
                    return Ok((&b""[..], AtomSymbol::LonePair));
                }
                if allow_queries {
                    match s {
                        b"A" | b"Q" | b"*" | b"X" | b"M" => {
                            if let Some(query) = QueryAtom::from_symbol_bytes(s) {
                                return Ok((&b""[..], AtomSymbol::Query(query)));
                            }
                        }
                        b"AH" | b"QH" | b"XH" | b"MH" => {
                            if allow_extended_queries {
                                if let Some(query) = QueryAtom::from_symbol_bytes(s) {
                                    return Ok((&b""[..], AtomSymbol::Query(query)));
                                }
                            }
                        }
                        b"L" => return Ok((&b""[..], AtomSymbol::AtomList(AtomList::default()))),
                        _ => {}
                    }
                }
                if allow_pseudoatoms {
                    if let Ok(s) = to_string(s) {
                        return Ok((&b""[..], AtomSymbol::Pseudoatom(s)));
                    }
                }
                Err(Err::Error(error::Error::new(s, error::ErrorKind::MapRes)))
            },
            true,
        )
        .parse(input)
        .and_then(|(remaining, symbol)| match symbol {
            Some(symbol) => Ok((remaining, symbol)),
            None => Err(Err::Error(error::Error::new(input, error::ErrorKind::Eof))),
        })
    }
}

/// Parse atom inputs with 52-69 characters (s. `atom_input` for more details).
/// Includes atom mapping number.
fn atom_input69<'inp, 'fl>(
    input: &'inp [u8],
    flags: &'fl CtabParseFlags,
) -> IResult<&'inp [u8], Atom, error::Error<&'inp [u8]>> {
    let strict_padding = flags.contains(CtabParseFlags::STRICT_PADDING);
    let x = fixed_width_float::<f64>(10, 4);
    let y = fixed_width_float::<f64>(10, 4);
    let z = fixed_width_float::<f64>(10, 4);
    let symbol = atom_symbol(flags);
    let mass_diff = map(fixed_width_int_in_range_opt::<i8, _>(2, -3..=4), |opt| {
        convert_atom_mass_diff_code(opt.unwrap_or(0))
    });
    let charge_radical = map(fixed_width_int_in_range_opt::<u8, _>(3, 0..=7), |opt| {
        convert_atom_charge_code(opt.unwrap_or(0))
    });
    let stereo_parity = map_res(fixed_width_int_in_range::<u8, _>(3, 0..=3), |code| {
        convert_atom_stereo_parity_code(code, false)
    });
    let padding1 = fixed_width_padding_n(2, 3, strict_padding);
    let valence = map_res(
        fixed_width_int_in_range::<u8, _>(3, 0..=15),
        convert_atom_valence_code,
    );
    let padding2 = fixed_width_padding_n(3, 3, strict_padding);
    let atom_map_num = fixed_width_int_in_range_opt::<u32, _>(3, 1..=999);
    let padding3 = fixed_width_padding_n(2, 3, strict_padding);

    map(
        (
            x,
            y,
            z,
            preceded(tag(" "), symbol),
            mass_diff,
            charge_radical,
            terminated(stereo_parity, padding1),
            terminated(valence, padding2),
            terminated(atom_map_num, padding3),
        ),
        |(x, y, z, symbol, mass_diff, charge_radical, stereo_parity, valence, atom_map_num)| {
            let (element, isotope_mass) = convert_atom_symbol_mass_diff(symbol, mass_diff);
            let (charge, radical) = charge_radical;
            Atom {
                element,
                charge,
                radical,
                isotope_mass,
                stereo_parity,
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

/// Parse atom inputs with 49-51 characters (s. `atom_input` for more details).
/// Includes mass difference, charge/radical and valence fields.
///
/// If `allow_named_isotopes` is true, allow named isotopes (D, T).
/// If `strict_padding` is true, require strict padding.
///
fn atom_input51<'inp, 'fl>(
    input: &'inp [u8],
    flags: &'fl CtabParseFlags,
) -> IResult<&'inp [u8], Atom, error::Error<&'inp [u8]>> {
    let strict_padding = flags.contains(CtabParseFlags::STRICT_PADDING);
    let x = fixed_width_float::<f64>(10, 4);
    let y = fixed_width_float::<f64>(10, 4);
    let z = fixed_width_float::<f64>(10, 4);
    let symbol = atom_symbol(flags);
    let mass_diff = map(fixed_width_int_in_range_opt::<i8, _>(2, -3..=4), |opt| {
        convert_atom_mass_diff_code(opt.unwrap_or(0))
    });
    let charge_radical = map(fixed_width_int_in_range_opt::<u8, _>(3, 0..=7), |opt| {
        convert_atom_charge_code(opt.unwrap_or(0))
    });
    let padding1 = fixed_width_padding_n(2, 3, strict_padding);
    let stereo_parity = map_res(fixed_width_int_in_range::<u8, _>(3, 0..=3), |code| {
        convert_atom_stereo_parity_code(code, false)
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
            preceded(tag(" "), symbol),
            mass_diff,
            charge_radical,
            terminated(stereo_parity, padding1),
            valence,
        ),
        |(x, y, z, symbol, mass_diff, charge_radical, stereo_parity, valence)| {
            let (element, isotope_mass) = convert_atom_symbol_mass_diff(symbol, mass_diff);
            let (charge, radical) = charge_radical;
            Atom {
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

/// Parse atom inputs with 42-48 characters including up to 6 characters of ignored data
/// (s. `atom_input` for more details). Lacks trailing valence and atom mapping fields
/// (substituted by defaults).
///
/// If `allow_named_isotopes` is true, allow named isotopes (D, T).
/// If `strict_padding` is true, require strict padding.
///
fn atom_input42<'inp, 'fl>(
    input: &'inp [u8],
    flags: &'fl CtabParseFlags,
) -> IResult<&'inp [u8], Atom, error::Error<&'inp [u8]>> {
    let strict_padding = flags.contains(CtabParseFlags::STRICT_PADDING);
    let x = fixed_width_float::<f64>(10, 4);
    let y = fixed_width_float::<f64>(10, 4);
    let z = fixed_width_float::<f64>(10, 4);
    let symbol = atom_symbol(flags);
    let mass_diff = map(
        fixed_width_int_in_range::<i8, _>(2, -3..=4),
        convert_atom_mass_diff_code,
    );
    let charge_radical = map(
        fixed_width_int_in_range::<u8, _>(3, 0..=7),
        convert_atom_charge_code,
    );
    let stereo_parity = map_res(fixed_width_int_in_range::<u8, _>(3, 0..=3), |code| {
        convert_atom_stereo_parity_code(code, false)
    });
    let n = input.len().saturating_sub(42) / 3;
    let padding1 = fixed_width_padding_n(n, 3, strict_padding);

    map(
        (
            x,
            y,
            z,
            preceded(tag(" "), symbol),
            mass_diff,
            charge_radical,
            terminated(stereo_parity, padding1),
        ),
        |(x, y, z, symbol, mass_diff, charge_radical, stereo_parity)| {
            let (element, isotope_mass) = convert_atom_symbol_mass_diff(symbol, mass_diff);
            let (charge, radical) = charge_radical;
            Atom {
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
            }
        },
    )
    .parse(input)
}

/// Parse atom inputs with 37-41 characters (s. `atom_input` for more details).
/// Lacks trailing stereo parity, valence and atom mapping fields (substituted by defaults).
///
/// If `allow_named_isotopes` is true, allow named isotopes (D, T).
///
fn atom_input39<'inp, 'fl>(
    input: &'inp [u8],
    flags: &'fl CtabParseFlags,
) -> IResult<&'inp [u8], Atom, error::Error<&'inp [u8]>> {
    let x = fixed_width_float::<f64>(10, 4);
    let y = fixed_width_float::<f64>(10, 4);
    let z = fixed_width_float::<f64>(10, 4);
    let symbol = atom_symbol(flags);
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
            preceded(tag(" "), symbol),
            mass_diff,
            charge_radical,
        ),
        |(x, y, z, symbol, mass_diff, charge_radical)| {
            let (element, isotope_mass) = convert_atom_symbol_mass_diff(symbol, mass_diff);
            let (charge, radical) = charge_radical;
            Atom {
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

/// Parse atom inputs with 35-36 characters (s. `atom_input` for more details).
/// Lacks trailing charge/radical, valence and atom mapping fields (substituted by defaults).
///
/// If `allow_named_isotopes` is true, allow named isotopes (D, T).
/// If `strict_padding` is true, require strict padding.
///
fn atom_input36<'inp, 'fl>(
    input: &'inp [u8],
    flags: &'fl CtabParseFlags,
) -> IResult<&'inp [u8], Atom, error::Error<&'inp [u8]>> {
    let x = fixed_width_float::<f64>(10, 4);
    let y = fixed_width_float::<f64>(10, 4);
    let z = fixed_width_float::<f64>(10, 4);
    let symbol = atom_symbol(flags);
    let mass_diff = map(fixed_width_int_in_range_opt::<i8, _>(2, -3..=4), |opt| {
        convert_atom_mass_diff_code(opt.unwrap_or(0))
    });

    map(
        (x, y, z, preceded(tag(" "), symbol), mass_diff),
        |(x, y, z, symbol, mass_diff)| {
            let (element, isotope_mass) = convert_atom_symbol_mass_diff(symbol, mass_diff);
            Atom {
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

/// Parse atom inputs with 34 characters (s. `atom_input` for more details).
/// Lacks trailing mass difference, charge/radical, valence and atom mapping fields
/// (substituted by defaults).
///
/// If `allow_named_isotopes` is true, allow named isotopes (D, T).
/// If `strict_padding` is true, require strict padding.
///
fn atom_input34<'inp, 'fl>(
    input: &'inp [u8],
    flags: &'fl CtabParseFlags,
) -> IResult<&'inp [u8], Atom, error::Error<&'inp [u8]>> {
    let x = fixed_width_float::<f64>(10, 4);
    let y = fixed_width_float::<f64>(10, 4);
    let z = fixed_width_float::<f64>(10, 4);
    let symbol = atom_symbol(flags);

    map(
        (x, y, z, preceded(tag(" "), symbol)),
        |(x, y, z, symbol)| {
            let (element, isotope_mass) = convert_atom_symbol_mass_diff(symbol, None);
            Atom {
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

/// Parse atom input (optimized for performance)
/// Fails immediately on atomlike symbols. For parsing all atom types, see atomlike_input.
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
pub fn atom_input<'inp, 'fl>(
    flags: &'fl CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = Atom, Error = error::Error<&'inp [u8]>> + use<'inp, 'fl> {
    move |input: &'inp [u8]| {
        let input = input.trim_end_with(|c| c == '\r' || c == '\n');
        let len = input.len();

        let parser = match len {
            67.. => atom_input69,
            49..=66 => atom_input51,
            42..=48 => atom_input42,
            37..=41 => atom_input39,
            35..=36 => atom_input36,
            32..=34 => atom_input34,
            _ => return Err(Err::Error(error::Error::new(input, error::ErrorKind::Eof))),
        };
        terminated(move |input| parser(input, flags), space0).parse(input)
    }
}

/// Parse atom and atom-like input
/// Allows all atom types. For faster parsing of basic molecules, see atom_input.
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
/// | HHH   | ignored            |              |                                           |   
/// | rrr   | ignored            |              |                                           |   
/// | iii   | ignored            |              |                                           |   
/// | mmm   | atom mapping       | 1..=#atoms   | Reaction, accepted as extension           |
/// | nnn   | inversion          | 0..=2        | Reaction                                  |
/// | eee   | exact change       | 0, 1         | Reaction                                  |
/// -----------------------------------------------------------------------------------------
pub fn atomlike_input<'inp, 'fl>(
    flags: &'fl CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = AtomLike, Error = error::Error<&'inp [u8]>> + use<'inp, 'fl> {
    move |input: &'inp [u8]| {
        let input = input.trim_end_with(|c| c == '\r' || c == '\n');
        terminated(atomlike_input_inner(flags), space0).parse(input)
    }
}

// Internal parser for atomlike_input
fn atomlike_input_inner<'inp, 'fl>(
    flags: &'fl CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = AtomLike, Error = error::Error<&'inp [u8]>> + use<'inp, 'fl> {
    let strict_padding = flags.contains(CtabParseFlags::STRICT_PADDING);
    move |input: &'inp [u8]| {
        // x, y, z coordinates
        let (i, x) = fixed_width_float::<f64>(10, 4).parse(input)?;
        let (i, y) = fixed_width_float::<f64>(10, 4).parse(i)?;
        let (i, z) = fixed_width_float::<f64>(10, 4).parse(i)?;

        // Atom symbol
        let (i, atom_symbol) = preceded(tag(" "), atomlike_symbol(flags)).parse(i)?;

        // Mass difference
        let (i, mass_diff) = map(fixed_width_int_in_range_opt::<i8, _>(2, -3..=4), |opt| {
            convert_atom_mass_diff_code(opt.unwrap_or(0))
        })
        .parse(i)?;

        // Charge/radical
        let (i, (charge, radical)) = map(fixed_width_int_in_range_opt::<u8, _>(3, 0..=7), |opt| {
            convert_atom_charge_code(opt.unwrap_or(0))
        })
        .parse(i)?;

        // Stereo parity
        let (i, stereo_parity) = cond(
            i.len() >= 3,
            map_res(fixed_width_int_in_range::<u8, _>(3, 0..=3), |code| {
                convert_atom_stereo_parity_code(code, true)
            }),
        )
        .parse(i)?;

        // Hydrogen count
        let extended_range = flags.contains(CtabParseFlags::EXTENDED_RANGE);
        let max_hydrogen = if extended_range { 255 } else { 5 };
        let (i, hydrogen_count) = cond(
            i.len() >= 3,
            map_res(
                fixed_width_int_in_range::<u8, _>(3, 0..=max_hydrogen),
                move |code| convert_atom_hydrogen_count_code(code, extended_range),
            ),
        )
        .parse(i)?;

        // Stereo care
        let (i, stereo_care) = cond(
            i.len() >= 3,
            map_res(fixed_width_int(3), convert_atom_stereo_care_code),
        )
        .parse(i)?;

        // Valence
        let (i, valence) = cond(
            i.len() >= 3,
            map_res(fixed_width_int(3), convert_atom_valence_code),
        )
        .parse(i)?;

        // Ignored fields
        let (i, _) = cond(
            !i.is_empty(),
            fixed_width_padding_n((i.len() / 3).min(3), 3, strict_padding),
        )
        .parse(i)?;

        // Atom mapping number
        let (i, atom_map_num) = cond(
            i.len() >= 3,
            fixed_width_int_in_range_opt::<u32, _>(3, 1..=999),
        )
        .parse(i)?;

        // Inversion flag
        let (i, inversion_flag) = cond(
            i.len() >= 3,
            map_res(
                fixed_width_int_in_range::<u8, _>(3, 0..=2),
                convert_atom_inversion_flag_code,
            ),
        )
        .parse(i)?;

        // Exact change flag
        let (i, exact_change_flag) = cond(
            i.len() >= 3,
            map_res(fixed_width_int(3), convert_atom_exact_change_flag_code),
        )
        .parse(i)?;

        // Combine atom mass difference and named isotope information
        let isotope_mass = match &atom_symbol {
            AtomSymbol::Element(e) => {
                mass_diff.map(|diff| (e.reference_mass_number() as i8 + diff) as u32)
            }
            AtomSymbol::NamedIsotope(i) => Some(i.mass_number()),
            _ => mass_diff.map(|diff| diff.unsigned_abs() as u32),
        };

        Ok((
            i,
            AtomLike {
                symbol: atom_symbol,
                charge,
                radical,
                isotope_mass,
                stereo_parity: stereo_parity.flatten(),
                hydrogen_count: hydrogen_count.flatten(),
                stereo_care: stereo_care.flatten(),
                valence: valence.flatten(),
                atom_map_num: atom_map_num.flatten(),
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
}

#[cfg(test)]
mod tests;
