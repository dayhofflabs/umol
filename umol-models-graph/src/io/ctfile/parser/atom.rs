//! Atom block parsers for CTab files.

use bstr::ByteSlice;
use nom::bytes::complete::tag;
use nom::character::complete::space0;
use nom::combinator::{all_consuming, cond, map, map_res};
use nom::error::{Error as NomError, ErrorKind as NomErrorKind};
use nom::sequence::{preceded, terminated};
use nom::{Err, IResult, Parser};
use umol_data::{Element, NamedIsotope};

use super::convert::{
    convert_atom_charge_code, convert_atom_exact_change_flag_code,
    convert_atom_hydrogen_count_code, convert_atom_inversion_flag_code,
    convert_atom_mass_diff_code, convert_atom_stereo_care_code, convert_atom_stereo_parity_code,
    convert_atom_symbol_mass_diff, convert_atom_valence_code,
};
use super::utils::{
    fixed_width_int, fixed_width_int_in_range, fixed_width_int_in_range_opt, fixed_width_partial,
    fixed_width_position, fixed_width_unused_n, is_reserved_atom_symbol, LinesWithOffsetExt,
};
use crate::io::ctfile::config::CtabParseFlags;
use crate::io::ctfile::error::ParseError;
use crate::io::ctfile::parser::rgroup::rgroup_symbol;
use crate::position::{all_zero, Point3D};
use crate::table_ir::{Atom, AtomList, AtomSymbol, ExtendedAtom, WildcardAtom};

/// Parse atom block (basic atoms only)
pub(super) fn atom_block<'inp>(
    atom_count: u32,
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (Vec<Atom>, Option<Vec<Point3D>>, u32), Error = ParseError>
       + use<'inp> {
    move |input: &'inp [u8]| {
        let ignore_positions = flags.contains(CtabParseFlags::IGNORE_POSITIONS);
        let mut atoms = Vec::with_capacity(atom_count as usize);
        let mut positions = Vec::with_capacity(if ignore_positions {
            0
        } else {
            atom_count as usize
        });
        let mut lines_iter = input.lines_with_offset();
        let mut byte_offset = 0;

        for line_index in 0..atom_count {
            let (line, byte_len) = lines_iter.next().ok_or_else(|| {
                Err::Error(ParseError::UnexpectedEof {
                    line: line_offset + line_index,
                    block: "atom",
                })
            })?;

            let (_, (atom, pos)) = all_consuming(terminated(atom_input(flags), space0))
                .parse(line)
                .map_err(|e| {
                    Err::Error(ParseError::atom_from_nom(e, line_offset + line_index, line))
                })?;
            atoms.push(atom);
            if !ignore_positions {
                positions.push(pos);
            }
            byte_offset += byte_len;
        }

        let remaining = &input[byte_offset..];
        if ignore_positions || (atom_count > 1 && all_zero(&positions)) {
            Ok((remaining, (atoms, None, line_offset + atom_count)))
        } else {
            Ok((
                remaining,
                (atoms, Some(positions), line_offset + atom_count),
            ))
        }
    }
}

/// Parse extended atom block
pub(super) fn extended_atom_block<'inp>(
    atom_count: u32,
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<
    &'inp [u8],
    Output = (Vec<ExtendedAtom>, Option<Vec<Point3D>>, u32),
    Error = ParseError,
> + use<'inp> {
    move |input: &'inp [u8]| {
        let ignore_positions = flags.contains(CtabParseFlags::IGNORE_POSITIONS);
        let mut atoms = Vec::with_capacity(atom_count as usize);
        let mut positions = Vec::with_capacity(if ignore_positions {
            0
        } else {
            atom_count as usize
        });
        let mut lines_iter = input.lines_with_offset();
        let mut byte_offset = 0;

        for line_index in 0..atom_count {
            let (line, byte_len) = lines_iter.next().ok_or_else(|| {
                Err::Error(ParseError::UnexpectedEof {
                    line: line_offset + line_index,
                    block: "atom",
                })
            })?;

            let (_, (atom, pos)) = all_consuming(terminated(extended_atom_input(flags), space0))
                .parse(line)
                .map_err(|e| {
                    Err::Error(ParseError::atom_from_nom(e, line_offset + line_index, line))
                })?;
            atoms.push(atom);
            if !ignore_positions {
                positions.push(pos);
            }
            byte_offset += byte_len;
        }

        let remaining = &input[byte_offset..];
        if ignore_positions || (atom_count > 1 && all_zero(&positions)) {
            Ok((remaining, (atoms, None, line_offset + atom_count)))
        } else {
            Ok((
                remaining,
                (atoms, Some(positions), line_offset + atom_count),
            ))
        }
    }
}

/// Parse basic atom input (optimized for performance)
/// Fails immediately on extended atom symbols. For parsing all atom types, see extended_atom_input.
///
/// xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee (69 characters wide)
///
/// *Values in the atom block*
/// ----------------------------------------------------------------------------------------------
/// | Field | Position | Meaning            | Values       | Notes                               |
/// |-------|----------|--------------------|--------------|--------------------------------------
/// | x,y,z |  1-30    | atom coordinates   |              | Generic, F10.4 format               |
/// |       | 31       | blank              |              |                                     |
/// | aaa   | 32-34    | atom symbol        | s. below     | Generic, Query, 3D, RGroup          |
/// | dd    | 35-36    | mass difference    | -3..=4       | Generic, s. also M  ISO             |
/// | ccc   | 37-39    | charge code        | 0..=7        | Generic, s. also M  CHG/M  RAD      |
/// | sss   | 40-42    | stereo parity      | 0..=3        | Generic, ignored when read          |
/// |       |          |                    |              | according to docs, used in practice |
/// | hhh   | 43-45    | hydrogen code      | 0..=5        | Query, H0 means no implicit Hs      |
/// |       |          |                    |              | Hn means >=n implicit Hs            |
/// | vvv   | 49-51    | valence code       | 0..=15       | Generic                             |
/// | mmm   | 61-63    | atom mapping       | 1..=#atoms   | Reaction, accepted as extension     |
/// ----------------------------------------------------------------------------------------------
///
/// *Behavior in unused and extended fields*
/// ---------------------------------------------------------------------
/// | Field    | Basic      | Basic strict | Extended | Extended strict |
/// |-------------------------------------------------------------------|
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
pub fn atom_input<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (Atom, Point3D), Error = NomError<&'inp [u8]>> + use<'inp> {
    move |input: &'inp [u8]| {
        let len = input.len();

        let parser = match len {
            61.. => atom_input69,
            49..=60 => atom_input60,
            40..=48 => atom_input48,
            37..=39 => atom_input39,
            35..=36 => atom_input36,
            32..=34 => atom_input34,
            _ => return Err(Err::Error(NomError::new(input, NomErrorKind::Eof))),
        };
        (move |input| parser(input, flags)).parse(input)
    }
}

/// Parse extended atom input
/// Allows all atom types. For faster parsing of basic molecules, see atom_input.
///
/// xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee (69 characters wide)
///
/// *Values in the atom block*
/// --------------------------------------------------------------------------------------------------
/// | Field | Position | Meaning          | Values       | Notes                                     |
/// |-------|----------|------------------|--------------|-------------------------------------------|
/// | x,y,z |  1-30    | atom coordinates |              | Generic, F10.4 format                     |
/// |       | 31       | blank            |              |                                           |
/// | aaa   | 32-34    | atom symbol      | s. above     | Generic, Query, 3D, RGroup                |
/// | dd    | 35-36    | mass difference  | -3..=4       | Generic, s. also M  ISO                   |
/// | ccc   | 37-39    | charge code      | 0..=7        | Generic, s. also M  CHG/M  RAD            |
/// | sss   | 40-42    | stereo parity    | 0..=3        | Generic, ignored when read according to   |
/// |       |          |                  |              | docs, used in practice                    |
/// | hhh   | 43-45    | hydrogen code    | 0..=5        | Query, H0 means no implicit Hs            |
/// |       |          |                  |              | Hn means >=n implicit Hs                  |
/// | bbb   | 46-48    | stereo care      | 0, 1         | Query, consider double bond stereo        |
/// |       |          |                  |              | when stereo care is 1 for both bond atoms |
/// | vvv   | 49-51    | valence code     | 0..=15       | Generic                                   |
/// | HHH   | 52-54    | ignored          |              |                                           |   
/// | rrr   | 55-57    | ignored          |              |                                           |   
/// | iii   | 58-60    | ignored          |              |                                           |   
/// | mmm   | 61-63    | atom mapping     | 1..=#atoms   | Reaction, accepted as extension           |
/// | nnn   | 64-66    | inversion        | 0..=2        | Reaction                                  |
/// | eee   | 67-69    | exact change     | 0, 1         | Reaction                                  |
/// --------------------------------------------------------------------------------------------------
///
pub fn extended_atom_input<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (ExtendedAtom, Point3D), Error = NomError<&'inp [u8]>> + use<'inp>
{
    extended_atom_input_inner(flags)
}

/// Parse atom symbol (Element and NamedIsotope only).
///
/// Returns error for extended atom symbols (L, A, Q, *, LP, R#, pseudoatoms).
/// If `allow_named_isotopes` is true, allow named isotopes (D, T).
///
fn atom_symbol<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = AtomSymbol, Error = NomError<&'inp [u8]>> {
    let allow_named_isotopes = flags.contains(CtabParseFlags::NAMED_ISOTOPES);
    move |input: &'inp [u8]| {
        map_res(
            fixed_width_partial(
                3,
                move |s: &'inp [u8]| {
                    let s = s.trim_ascii();
                    Element::from_symbol_bytes(s)
                        .map(AtomSymbol::Element)
                        .or_else(|| {
                            allow_named_isotopes
                                .then(|| NamedIsotope::from_symbol_bytes(s))
                                .flatten()
                                .map(AtomSymbol::NamedIsotope)
                        })
                        .ok_or(Err::Error(NomError::new(s, NomErrorKind::MapRes)))
                        .map(|symbol| (&b""[..], symbol))
                },
                true,
            ),
            move |symbol| symbol.ok_or(NomError::new(input, NomErrorKind::MapRes)),
        )
        .parse(input)
    }
}

/// Parse extended atom symbol (all atom types allowed in MOL specification).
///
/// --------------------------------------------------------------------------------
/// | Symbol      | Type          | Parser* | Notes                                |
/// --------------------------------------------------------------------------------
/// | H-Og        | Element       | B, E    |                                      |
/// | D, T        | Named Isotope | B, E    | Heavy H isotopes as extension        |
/// | L           | Atom List     | E       | Query molecules                      |
/// | *,A,Q,X,M   | Wildcard Atom | E       | Query molecules, rarely in oligomers |
/// | AH,QH,XH,MH | Wildcard Atom | E       | Query molecules, CXSMILES extension  |
/// | LP          | Lone Pair     | E       | Rarely used                          |
/// | R#          | R Group       | E       | Query molecules                      |
/// | any string  | Pseudoatom    | E       | Any string as pseudoatom              |
/// --------------------------------------------------------------------------------
/// | *Parsers: B: basic, E: extended                                               |
/// --------------------------------------------------------------------------------
///
/// If `allow_named_isotopes` is true, allow named isotopes (D, T).
/// If `allow_rgroups` is true, allow rgroups (R#).
/// If `allow_wildcards` is true, allow wildcard atoms (A, Q, *, X, M) and atom lists (L).
/// If `allow_chemaxon_wildcards` is true, allow CXSMILES wildcard atoms (AH, QH, XH, MH).
/// If `allow_electrons` is true, allow electrons (LP).
/// If `allow_pseudoatoms` is true, allow pseudoatoms (any string).
///
fn extended_atom_symbol<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = AtomSymbol, Error = NomError<&'inp [u8]>> {
    let allow_wildcards = flags.contains(CtabParseFlags::WILDCARDS);
    let allow_chemaxon_wildcards = flags.contains(CtabParseFlags::CHEMAXON_WILDCARDS);
    let allow_electrons = flags.contains(CtabParseFlags::ELECTRONS);
    let allow_rgroups = flags.contains(CtabParseFlags::RGROUPS);
    let allow_named_isotopes = flags.contains(CtabParseFlags::NAMED_ISOTOPES);
    let allow_pseudoatoms = flags.contains(CtabParseFlags::PSEUDOATOMS);
    move |input: &'inp [u8]| {
        map_res(
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
                    if allow_wildcards {
                        match s {
                            b"A" | b"Q" | b"*" | b"X" | b"M" => {
                                if let Some(wildcard) = WildcardAtom::from_symbol_bytes(s) {
                                    return Ok((&b""[..], AtomSymbol::WildcardAtom(wildcard)));
                                }
                            }
                            b"AH" | b"QH" | b"XH" | b"MH" => {
                                if allow_chemaxon_wildcards {
                                    if let Some(wildcard) = WildcardAtom::from_symbol_bytes(s) {
                                        return Ok((&b""[..], AtomSymbol::WildcardAtom(wildcard)));
                                    }
                                } else {
                                    return Err(Err::Error(NomError::new(s, NomErrorKind::MapRes)));
                                }
                            }
                            b"L" => return Ok((&b""[..], AtomSymbol::AtomList(AtomList::empty()))),
                            _ => {} // Fall through
                        }
                    }
                    if allow_rgroups {
                        if let Ok((_, rgroup)) = rgroup_symbol(s) {
                            return Ok((&b""[..], AtomSymbol::RGroup(rgroup)));
                        }
                    }
                    if allow_electrons && s == b"LP" {
                        return Ok((&b""[..], AtomSymbol::LonePair));
                    }

                    // Reject reserved atom symbols if corresponding flag is not set
                    if is_reserved_atom_symbol(
                        s,
                        allow_named_isotopes,
                        allow_wildcards,
                        allow_chemaxon_wildcards,
                        allow_electrons,
                        allow_rgroups,
                    ) {
                        return Err(Err::Error(NomError::new(s, NomErrorKind::MapRes)));
                    }

                    if allow_pseudoatoms && s.is_ascii() {
                        let s = s.to_str_lossy().into_owned();
                        return Ok((&b""[..], AtomSymbol::Pseudoatom(s)));
                    }
                    Err(Err::Error(NomError::new(s, NomErrorKind::MapRes)))
                },
                true,
            ),
            move |symbol| symbol.ok_or(NomError::new(input, NomErrorKind::MapRes)),
        )
        .parse(input)
    }
}

/// Parse basic atom input with 61-69 characters (s. `atom_input` for more details).
/// Includes atom mapping number (stored as class, see below).
///
/// Includes extended inversion/retention and exact change fields.
///
/// If `allow_named_isotopes` is true, allow named isotopes (D, T).
/// If `ignore_positions` is true, set atom positions to zero.
/// If `skip_unused_fields` is true, do no validate unused fields.
/// If `atom_map_hcount_fields` is true, set class to atom mapping number.
/// If `extended_range` is true, allow extended hydrogen count range.
///
fn atom_input69<'inp>(
    input: &'inp [u8],
    flags: CtabParseFlags,
) -> IResult<&'inp [u8], (Atom, Point3D), NomError<&'inp [u8]>> {
    let skip_unused_fields = flags.contains(CtabParseFlags::SKIP_UNUSED_FIELDS);
    let ignore_positions = flags.contains(CtabParseFlags::IGNORE_POSITIONS);
    let atom_map_hcount_fields = flags.contains(CtabParseFlags::ATOM_MAP_HCOUNT_FIELDS);
    let extended_range = flags.contains(CtabParseFlags::EXTENDED_RANGE);
    let position = fixed_width_position(ignore_positions);
    let symbol = atom_symbol(flags);
    let mass_diff = map(fixed_width_int_in_range_opt::<i8, _>(2, -3..=4), |opt| {
        convert_atom_mass_diff_code(opt.unwrap_or(0))
    });
    let charge_radical = map(fixed_width_int_in_range_opt::<u8, _>(3, 0..=7), |opt| {
        convert_atom_charge_code(opt.unwrap_or(0))
    });
    let stereo_parity = map_res(fixed_width_int_in_range::<u8, _>(3, 0..=3), convert_atom_stereo_parity_code);
    let max_hydrogen_count = if extended_range { 13 } else { 5 };
    let hydrogen_count = cond(
        atom_map_hcount_fields,
        map_res(
            fixed_width_int_in_range::<u8, _>(3, 0..=max_hydrogen_count),
            move |code| convert_atom_hydrogen_count_code(code, extended_range),
        ),
    );
    let count1 = if atom_map_hcount_fields { 1 } else { 2 };
    let extended1 = fixed_width_unused_n(count1, 3, false);
    let valence = map_res(
        fixed_width_int_in_range::<u8, _>(3, 0..=15),
        convert_atom_valence_code,
    );
    let unused2 = fixed_width_unused_n(3, 3, skip_unused_fields);
    let atom_map_num = cond(
        atom_map_hcount_fields,
        fixed_width_int_in_range_opt::<u32, _>(3, 1..=999),
    );
    let remaining_len = input
        .len()
        .min(69)
        .saturating_sub(if atom_map_hcount_fields { 63 } else { 60 });
    let count3 = remaining_len / 3;
    let extended3 = fixed_width_unused_n(count3, 3, false);

    map(
        (
            position,
            preceded(tag(" "), symbol),
            mass_diff,
            charge_radical,
            stereo_parity,
            terminated(hydrogen_count, extended1),
            terminated(valence, unused2),
            terminated(atom_map_num, extended3),
        ),
        |(
            position,
            symbol,
            mass_diff,
            charge_radical,
            chirality,
            hydrogen_count,
            valence,
            atom_map_num,
        )| {
            let (element, isotope) = convert_atom_symbol_mass_diff(symbol, mass_diff);
            let (charge, unpaired_e) = charge_radical;
            (
                Atom {
                    element,
                    charge,
                    isotope_mass: isotope,
                    hydrogens: hydrogen_count.flatten(),
                    implicit_h: false,
                    valence,
                    unpaired_e,
                    aromatic: None,
                    chirality,
                    class: atom_map_num.flatten(),
                    span: None,
                    label: None,
                    value: None,
                },
                position,
            )
        },
    )
    .parse(input)
}

/// Parse basic atom input with 49-60 characters (s. `atom_input` for more details).
/// Includes mass difference, charge/radical and valence fields. Lacks atom mapping field.
/// Includes extended hydrogen code and stereo care fields and three ignored fields.
///
/// If `allow_named_isotopes` is true, allow named isotopes (D, T).
/// If `ignore_positions` is true, set atom positions to zero.
/// If `skip_unused_fields` is true, do no validate unused fields.
/// If `atom_map_hcount_fields` is true, include hydrogen count field.
/// If `extended_range` is true, allow extended hydrogen count range.
///
fn atom_input60<'inp>(
    input: &'inp [u8],
    flags: CtabParseFlags,
) -> IResult<&'inp [u8], (Atom, Point3D), NomError<&'inp [u8]>> {
    let ignore_positions = flags.contains(CtabParseFlags::IGNORE_POSITIONS);
    let skip_unused_fields = flags.contains(CtabParseFlags::SKIP_UNUSED_FIELDS);
    let atom_map_hcount_fields = flags.contains(CtabParseFlags::ATOM_MAP_HCOUNT_FIELDS);
    let extended_range = flags.contains(CtabParseFlags::EXTENDED_RANGE);
    let position = fixed_width_position(ignore_positions);
    let symbol = atom_symbol(flags);
    let mass_diff = map(fixed_width_int_in_range_opt::<i8, _>(2, -3..=4), |opt| {
        convert_atom_mass_diff_code(opt.unwrap_or(0))
    });
    let charge_radical = map(fixed_width_int_in_range_opt::<u8, _>(3, 0..=7), |opt| {
        convert_atom_charge_code(opt.unwrap_or(0))
    });
    let stereo_parity = map_res(fixed_width_int_in_range::<u8, _>(3, 0..=3), convert_atom_stereo_parity_code);
    let max_hydrogen_count = if extended_range { 13 } else { 5 };
    let hydrogen_count = cond(
        atom_map_hcount_fields,
        map_res(
            fixed_width_int_in_range::<u8, _>(3, 0..=max_hydrogen_count),
            move |code| convert_atom_hydrogen_count_code(code, extended_range),
        ),
    );
    let count1 = if atom_map_hcount_fields { 1 } else { 2 };
    let extended1 = fixed_width_unused_n(count1, 3, false);
    let valence = map_res(
        fixed_width_int_in_range::<u8, _>(3, 0..=15),
        convert_atom_valence_code,
    );
    let remaining_len = input.len().saturating_sub(51);
    let count2 = remaining_len / 3;
    let unused2 = fixed_width_unused_n(count2, 3, skip_unused_fields);

    map(
        (
            position,
            preceded(tag(" "), symbol),
            mass_diff,
            charge_radical,
            stereo_parity,
            terminated(hydrogen_count, extended1),
            terminated(valence, unused2),
        ),
        |(position, symbol, mass_diff, charge_radical, chirality, hydrogen_count, valence)| {
            let (element, isotope) = convert_atom_symbol_mass_diff(symbol, mass_diff);
            let (charge, unpaired_e) = charge_radical;
            (
                Atom {
                    element,
                    charge,
                    isotope_mass: isotope,
                    hydrogens: hydrogen_count.flatten(),
                    implicit_h: false,
                    valence,
                    unpaired_e,
                    aromatic: None,
                    chirality,
                    class: None,
                    span: None,
                    label: None,
                    value: None,
                },
                position,
            )
        },
    )
    .parse(input)
}

/// Parse basic atom input with 40-48 characters including up to 6 characters of ignored data
/// (s. `atom_input` for more details). Lacks trailing valence and atom mapping fields.
/// Includes extended hydrogen code and stereo care fields.
///
/// If `allow_named_isotopes` is true, allow named isotopes (D, T).
/// If `ignore_positions` is true, set atom positions to zero.
/// If `atom_map_hcount_fields` is true, include hydrogen count field.
/// If `extended_range` is true, allow extended hydrogen count range.
///
fn atom_input48<'inp>(
    input: &'inp [u8],
    flags: CtabParseFlags,
) -> IResult<&'inp [u8], (Atom, Point3D), NomError<&'inp [u8]>> {
    let ignore_positions = flags.contains(CtabParseFlags::IGNORE_POSITIONS);
    let atom_map_hcount_fields = flags.contains(CtabParseFlags::ATOM_MAP_HCOUNT_FIELDS);
    let extended_range = flags.contains(CtabParseFlags::EXTENDED_RANGE);
    let position = fixed_width_position(ignore_positions);
    let symbol = atom_symbol(flags);
    let mass_diff = map(
        fixed_width_int_in_range::<i8, _>(2, -3..=4),
        convert_atom_mass_diff_code,
    );
    let charge_radical = map(
        fixed_width_int_in_range::<u8, _>(3, 0..=7),
        convert_atom_charge_code,
    );
    let stereo_parity = map_res(fixed_width_int_in_range::<u8, _>(3, 0..=3), convert_atom_stereo_parity_code);
    let max_hydrogen_count = if extended_range { 13 } else { 5 };
    let hydrogen_count = cond(
        atom_map_hcount_fields,
        map_res(
            fixed_width_int_in_range::<u8, _>(3, 0..=max_hydrogen_count),
            move |code| convert_atom_hydrogen_count_code(code, extended_range),
        ),
    );
    let count1 = input
        .len()
        .saturating_sub(if atom_map_hcount_fields { 45 } else { 42 })
        / 3;
    let extended1 = fixed_width_unused_n(count1, 3, false);

    map(
        (
            position,
            preceded(tag(" "), symbol),
            mass_diff,
            charge_radical,
            stereo_parity,
            terminated(hydrogen_count, extended1),
        ),
        |(position, symbol, mass_diff, charge_radical, chirality, hydrogen_count)| {
            let (element, isotope) = convert_atom_symbol_mass_diff(symbol, mass_diff);
            let (charge, unpaired_e) = charge_radical;
            (
                Atom {
                    element,
                    charge,
                    isotope_mass: isotope,
                    hydrogens: hydrogen_count.flatten(),
                    implicit_h: false,
                    valence: None,
                    unpaired_e,
                    aromatic: None,
                    chirality,
                    class: None,
                    span: None,
                    label: None,
                    value: None,
                },
                position,
            )
        },
    )
    .parse(input)
}

/// Parse basic atom input with 37-39 characters (s. `atom_input` for more details).
/// Lacks trailing stereo parity, valence and atom mapping fields.
///
/// If `allow_named_isotopes` is true, allow named isotopes (D, T).
/// If `ignore_positions` is true, set atom positions to zero.
///
fn atom_input39<'inp>(
    input: &'inp [u8],
    flags: CtabParseFlags,
) -> IResult<&'inp [u8], (Atom, Point3D), NomError<&'inp [u8]>> {
    let ignore_positions = flags.contains(CtabParseFlags::IGNORE_POSITIONS);
    let position = fixed_width_position(ignore_positions);
    let symbol = atom_symbol(flags);
    let mass_diff = map(fixed_width_int_in_range_opt::<i8, _>(2, -3..=4), |opt| {
        convert_atom_mass_diff_code(opt.unwrap_or(0))
    });
    let charge_radical = map(fixed_width_int_in_range_opt::<u8, _>(3, 0..=7), |opt| {
        convert_atom_charge_code(opt.unwrap_or(0))
    });

    map(
        (
            position,
            preceded(tag(" "), symbol),
            mass_diff,
            charge_radical,
        ),
        |(position, symbol, mass_diff, charge_radical)| {
            let (element, isotope) = convert_atom_symbol_mass_diff(symbol, mass_diff);
            let (charge, unpaired_e) = charge_radical;
            (
                Atom {
                    element,
                    charge,
                    isotope_mass: isotope,
                    hydrogens: None,
                    implicit_h: false,
                    valence: None,
                    unpaired_e,
                    aromatic: None,
                    chirality: None,
                    class: None,
                    span: None,
                    label: None,
                    value: None,
                },
                position,
            )
        },
    )
    .parse(input)
}

/// Parse basic atom input with 34-36 characters (s. `atom_input` for more details).
/// Lacks trailing charge/radical, valence and atom mapping fields.
///
/// If `allow_named_isotopes` is true, allow named isotopes (D, T).
/// If `ignore_positions` is true, set atom positions to zero.
///
fn atom_input36<'inp>(
    input: &'inp [u8],
    flags: CtabParseFlags,
) -> IResult<&'inp [u8], (Atom, Point3D), NomError<&'inp [u8]>> {
    let ignore_positions = flags.contains(CtabParseFlags::IGNORE_POSITIONS);
    let position = fixed_width_position(ignore_positions);
    let symbol = atom_symbol(flags);
    let mass_diff = map(fixed_width_int_in_range_opt::<i8, _>(2, -3..=4), |opt| {
        convert_atom_mass_diff_code(opt.unwrap_or(0))
    });

    map(
        (position, preceded(tag(" "), symbol), mass_diff),
        |(position, symbol, mass_diff)| {
            let (element, isotope) = convert_atom_symbol_mass_diff(symbol, mass_diff);
            (
                Atom {
                    element,
                    charge: None,
                    isotope_mass: isotope,
                    hydrogens: None,
                    implicit_h: false,
                    valence: None,
                    unpaired_e: None,
                    aromatic: None,
                    chirality: None,
                    class: None,
                    span: None,
                    label: None,
                    value: None,
                },
                position,
            )
        },
    )
    .parse(input)
}

/// Parse basic atom input with 31-33 characters (s. `atom_input` for more details).
/// Lacks trailing mass difference, charge/radical, valence and atom mapping fields.
///
/// If `allow_named_isotopes` is true, allow named isotopes (D, T).
/// If `ignore_positions` is true, set atom positions to zero.
///
fn atom_input34<'inp>(
    input: &'inp [u8],
    flags: CtabParseFlags,
) -> IResult<&'inp [u8], (Atom, Point3D), NomError<&'inp [u8]>> {
    let ignore_positions = flags.contains(CtabParseFlags::IGNORE_POSITIONS);
    let position = fixed_width_position(ignore_positions);
    let symbol = atom_symbol(flags);

    map(
        (position, preceded(tag(" "), symbol)),
        |(position, symbol)| {
            let (element, isotope) = convert_atom_symbol_mass_diff(symbol, None);
            (
                Atom {
                    element,
                    charge: None,
                    isotope_mass: isotope,
                    hydrogens: None,
                    implicit_h: false,
                    valence: None,
                    unpaired_e: None,
                    aromatic: None,
                    chirality: None,
                    class: None,
                    span: None,
                    label: None,
                    value: None,
                },
                position,
            )
        },
    )
    .parse(input)
}

/// Internal parser for extended_atom_input
fn extended_atom_input_inner<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (ExtendedAtom, Point3D), Error = NomError<&'inp [u8]>> + use<'inp>
{
    let skip_unused_fields = flags.contains(CtabParseFlags::SKIP_UNUSED_FIELDS);
    let ignore_positions = flags.contains(CtabParseFlags::IGNORE_POSITIONS);
    let extended_range = flags.contains(CtabParseFlags::EXTENDED_RANGE);
    move |input: &'inp [u8]| {
        // x, y, z coordinates
        let (i, position) = fixed_width_position(ignore_positions).parse(input)?;

        // Atom symbol
        let (i, atom_symbol) = preceded(tag(" "), extended_atom_symbol(flags)).parse(i)?;

        // Mass difference
        let (i, mass_diff) = map(fixed_width_int_in_range_opt::<i8, _>(2, -3..=4), |opt| {
            convert_atom_mass_diff_code(opt.unwrap_or(0))
        })
        .parse(i)?;

        // Charge/radical
        let (i, (charge, unpaired_e)) =
            map(fixed_width_int_in_range_opt::<u8, _>(3, 0..=7), |opt| {
                convert_atom_charge_code(opt.unwrap_or(0))
            })
            .parse(i)?;

        // Stereo parity
        let (i, stereo_parity) = cond(
            i.len() >= 3,
            map_res(fixed_width_int_in_range::<u8, _>(3, 0..=3), convert_atom_stereo_parity_code),
        )
        .parse(i)?;

        // Hydrogen count
        let max_hydrogen_count = if extended_range { 13 } else { 5 };
        let (i, hydrogen_count) = cond(
            i.len() >= 3,
            map_res(
                fixed_width_int_in_range::<u8, _>(3, 0..=max_hydrogen_count),
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
            fixed_width_unused_n((i.len() / 3).min(3), 3, skip_unused_fields),
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
            (
                ExtendedAtom {
                    symbol: atom_symbol,
                    charge,
                    isotope_mass,
                    unpaired_e,
                    hydrogens: hydrogen_count.flatten(),
                    stereo_care: stereo_care.flatten(),
                    valence: valence.flatten(),
                    inversion_retention: inversion_flag.flatten(),
                    exact_change: exact_change_flag.flatten(),
                    implicit_h: false,
                    aromatic: None,
                    chirality: stereo_parity.flatten(),
                    class: atom_map_num.flatten(),
                    span: None,
                    label: None,
                    value: None,
                    ring_bond_count: None,
                    substitution_count: None,
                    unsaturated: None,
                    link_atom: None,
                    attachment_point: None,
                    attachment_order: None,
                    properties: std::collections::HashMap::new(),
                },
                position,
            ),
        ))
    }
}

#[cfg(test)]
mod tests;
