//! Atom block parsers for CTab files.

use bstr::ByteSlice;
use nom::bytes::complete::tag;
use nom::character::complete::space0;
use nom::combinator::{all_consuming, cond, map, map_res};
use nom::error::{Error as NomError, ErrorKind as NomErrorKind};
use nom::sequence::{preceded, terminated};
use nom::{Err, Parser};
use umol_data::{Element, NamedIsotope};

use super::convert::{
    convert_atom_charge_code, convert_atom_exact_change_flag_code,
    convert_atom_hydrogen_count_code, convert_atom_inversion_flag_code,
    convert_atom_mass_diff_code, convert_atom_stereo_care_code, convert_atom_stereo_parity_code,
    convert_atom_symbol_mass_diff, convert_atom_valence_code,
};
use super::utils::{
    fixed_width_float_f10_4, fixed_width_int, fixed_width_int_in_range,
    fixed_width_int_in_range_opt, fixed_width_int_opt_signed, fixed_width_int_opt_unsigned,
    fixed_width_partial, fixed_width_position, fixed_width_unused_n, is_all_whitespace_or_zeroes,
    is_reserved_atom_symbol, LinesWithOffsetExt,
};
use crate::io::ctfile::config::CtabParseFlags;
use crate::io::ctfile::error::ParseError;
use crate::io::ctfile::parser::convert::convert_extended_atom_symbol_mass_diff;
use crate::io::ctfile::parser::rgroup::rgroup_symbol;
use crate::io::ctfile::parser::utils::fixed_width_unused;
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
        if input.len() < 32 {
            return Err(Err::Error(NomError::new(input, NomErrorKind::Eof)));
        }
        if flags == CtabParseFlags::BASIC && input.len() >= 69 {
            return basic_atom_input69().parse(input);
        }

        let skip_unused_fields = flags.contains(CtabParseFlags::SKIP_UNUSED_FIELDS);
        let ignore_positions = flags.contains(CtabParseFlags::IGNORE_POSITIONS);
        let atom_map_hcount_fields = flags.contains(CtabParseFlags::ATOM_MAP_HCOUNT_FIELDS);
        let extended_range = flags.contains(CtabParseFlags::EXTENDED_RANGE);

        // x, y, z coordinates
        let (remaining, position) =
            fixed_width_position(ignore_positions, skip_unused_fields).parse(input)?;

        // Atom symbol
        let (remaining, symbol) = preceded(tag(" "), atom_symbol(flags)).parse(remaining)?;

        // Mass difference
        let (remaining, mass_diff) = map(
            fixed_width_int_in_range_opt::<i8, _>(2, -3..=4),
            |opt| convert_atom_mass_diff_code(opt.unwrap_or(0)),
        )
        .parse(remaining)?;

        // Combine atom mass difference and named isotope information
        let (element, isotope_mass) = convert_atom_symbol_mass_diff(&symbol, mass_diff);

        // Charge/radical
        let (remaining, (charge, unpaired_e)) =
            map(fixed_width_int_in_range_opt::<u8, _>(3, 0..=7), |opt| {
                convert_atom_charge_code(opt.unwrap_or(0))
            })
            .parse(remaining)?;

        // Stereo parity
        let (remaining, stereo_parity) = cond(
            remaining.len() >= 3,
            map_res(
                fixed_width_int_in_range::<u8, _>(3, 0..=3),
                convert_atom_stereo_parity_code,
            ),
        )
        .parse(remaining)?;

        // Hydrogen count
        let (remaining, hydrogen_count) = cond(
            remaining.len() >= 3,
            map_res(fixed_width_int_in_range::<u8, _>(3, 0..=13), move |code| {
                convert_atom_hydrogen_count_code(code, extended_range)
            }),
        )
        .parse(remaining)?;

        // Unused field
        let (remaining, _) =
            cond(remaining.len() >= 3, fixed_width_unused(3, false)).parse(remaining)?;

        // Valence
        let (remaining, valence) = cond(
            remaining.len() >= 3,
            map_res(
                fixed_width_int_in_range::<u8, _>(3, 0..=15),
                convert_atom_valence_code,
            ),
        )
        .parse(remaining)?;

        // Ignored fields
        let (remaining, _) = cond(
            !remaining.is_empty(),
            fixed_width_unused_n((remaining.len() / 3).min(3), 3, skip_unused_fields),
        )
        .parse(remaining)?;

        // Atom mapping number
        let (remaining, atom_map_num) = cond(
            remaining.len() >= 3,
            fixed_width_int_in_range_opt::<u32, _>(3, 1..=999),
        )
        .parse(remaining)?;

        // Unused fields
        let (remaining, _) = cond(
            !remaining.is_empty(),
            fixed_width_unused_n((remaining.len() / 3).min(2), 3, false),
        )
        .parse(remaining)?;

        // Verify hydrogen count
        let hydrogen_count = hydrogen_count.flatten();
        let hydrogens = if atom_map_hcount_fields {
            Ok(hydrogen_count)
        } else if hydrogen_count.is_some() {
            Err(Err::Error(NomError::new(input, NomErrorKind::Verify)))
        } else {
            Ok(None)
        }?;

        // Verify atom mapping number
        let atom_map_num = atom_map_num.flatten();
        let class = if atom_map_hcount_fields {
            Ok(atom_map_num)
        } else if atom_map_num.is_some() {
            Err(Err::Error(NomError::new(input, NomErrorKind::Verify)))
        } else {
            Ok(None)
        }?;

        Ok((
            remaining,
            (
                Atom {
                    element,
                    charge,
                    isotope_mass,
                    hydrogens,
                    implicit_h: false,
                    valence: valence.flatten(),
                    unpaired_e,
                    aromatic: None,
                    chirality: stereo_parity.flatten(),
                    class,
                    span: None,
                    label: None,
                    value: None,
                },
                position,
            ),
        ))
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
    let skip_unused_fields = flags.contains(CtabParseFlags::SKIP_UNUSED_FIELDS);
    let ignore_positions = flags.contains(CtabParseFlags::IGNORE_POSITIONS);
    let extended_range = flags.contains(CtabParseFlags::EXTENDED_RANGE);
    move |input: &'inp [u8]| {
        if input.len() < 32 {
            return Err(Err::Error(NomError::new(input, NomErrorKind::Eof)));
        }

        // x, y, z coordinates
        let (remaining, position) =
            fixed_width_position(ignore_positions, skip_unused_fields).parse(input)?;

        // Atom symbol
        let (remaining, symbol) =
            preceded(tag(" "), extended_atom_symbol(flags)).parse(remaining)?;

        // Mass difference
        let (remaining, mass_diff) = map(
            fixed_width_int_in_range_opt::<i8, _>(2, -3..=4),
            |opt| convert_atom_mass_diff_code(opt.unwrap_or(0)),
        )
        .parse(remaining)?;

        // Combine atom mass difference and named isotope information
        let isotope_mass = convert_extended_atom_symbol_mass_diff(&symbol, mass_diff);

        // Charge/radical
        let (remaining, (charge, unpaired_e)) =
            map(fixed_width_int_in_range_opt::<u8, _>(3, 0..=7), |opt| {
                convert_atom_charge_code(opt.unwrap_or(0))
            })
            .parse(remaining)?;

        // Stereo parity
        let (remaining, stereo_parity) = cond(
            remaining.len() >= 3,
            map_res(
                fixed_width_int_in_range::<u8, _>(3, 0..=3),
                convert_atom_stereo_parity_code,
            ),
        )
        .parse(remaining)?;

        // Hydrogen count
        let (remaining, hydrogen_count) = cond(
            remaining.len() >= 3,
            map_res(
                fixed_width_int_in_range::<u8, _>(3, 0..=13),
                move |code| convert_atom_hydrogen_count_code(code, extended_range),
            ),
        )
        .parse(remaining)?;

        // Stereo care
        let (remaining, stereo_care) = cond(
            remaining.len() >= 3,
            map_res(fixed_width_int(3), convert_atom_stereo_care_code),
        )
        .parse(remaining)?;

        // Valence
        let (remaining, valence) = cond(
            remaining.len() >= 3,
            map_res(fixed_width_int(3), convert_atom_valence_code),
        )
        .parse(remaining)?;

        // Ignored fields
        let (remaining, _) = cond(
            !remaining.is_empty(),
            fixed_width_unused_n((remaining.len() / 3).min(3), 3, skip_unused_fields),
        )
        .parse(remaining)?;

        // Atom mapping number
        let (remaining, atom_map_num) = cond(
            remaining.len() >= 3,
            fixed_width_int_in_range_opt::<u32, _>(3, 1..=999),
        )
        .parse(remaining)?;

        // Inversion flag
        let (remaining, inversion_flag) = cond(
            remaining.len() >= 3,
            map_res(
                fixed_width_int_in_range::<u8, _>(3, 0..=2),
                convert_atom_inversion_flag_code,
            ),
        )
        .parse(remaining)?;

        // Exact change flag
        let (remaining, exact_change_flag) = cond(
            remaining.len() >= 3,
            map_res(fixed_width_int(3), convert_atom_exact_change_flag_code),
        )
        .parse(remaining)?;

        Ok((
            remaining,
            (
                ExtendedAtom {
                    symbol,
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
                    pattern: None,
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

/// Fast-path basic atom input parser for 69-character lines.
/// Hard-codes CtabParseFlags::BASIC = NAMED_ISOTOPES | ATOM_MAP_HCOUNT_FIELDS | SKIP_UNUSED_FIELDS
/// behavior.
///
pub fn basic_atom_input69<'inp>(
) -> impl Parser<&'inp [u8], Output = (Atom, Point3D), Error = NomError<&'inp [u8]>> + use<'inp> {
    move |input: &'inp [u8]| {
        if input.len() < 69 {
            return Err(Err::Error(NomError::new(input, NomErrorKind::Eof)));
        }

        let line = &input[..69];
        let remaining = &input[69..];

        let x = fixed_width_float_f10_4::<f64>().parse(&line[0..10])?.1;
        let y = fixed_width_float_f10_4::<f64>().parse(&line[10..20])?.1;
        let z = fixed_width_float_f10_4::<f64>().parse(&line[20..30])?.1;
        let position = Point3D::new(x, y, z);

        if line[30] != b' ' {
            return Err(Err::Error(NomError::new(input, NomErrorKind::Char)));
        }

        let symbol_field = &line[31..34];
        let symbol_trimmed = symbol_field.trim_ascii();
        if symbol_trimmed.is_empty() {
            return Err(Err::Error(NomError::new(input, NomErrorKind::MapRes)));
        }
        let symbol = Element::from_symbol_bytes(symbol_trimmed)
            .map(AtomSymbol::Element)
            .or_else(|| {
                NamedIsotope::from_symbol_bytes(symbol_trimmed).map(AtomSymbol::NamedIsotope)
            })
            .ok_or_else(|| Err::Error(NomError::new(input, NomErrorKind::MapRes)))?;

        let mass_diff = fixed_width_int_opt_signed(input, &line[34..36])?;
        let mass_diff = mass_diff.filter(|val| (-3..=4).contains(val)).unwrap_or(0);
        let (element, isotope_mass) =
            convert_atom_symbol_mass_diff(&symbol, convert_atom_mass_diff_code(mass_diff as i8));

        let charge_code = fixed_width_int_opt_signed(input, &line[36..39])?;
        let charge_code = charge_code.filter(|val| (0..=7).contains(val)).unwrap_or(0);
        let (charge, unpaired_e) = convert_atom_charge_code(charge_code as u8);

        let stereo_code = fixed_width_int_opt_signed(input, &line[39..42])?.unwrap_or(0);
        if !(0..=3).contains(&stereo_code) {
            return Err(Err::Error(NomError::new(input, NomErrorKind::Verify)));
        }
        let chirality = convert_atom_stereo_parity_code(stereo_code as u8)
            .map_err(|_| Err::Error(NomError::new(input, NomErrorKind::Verify)))?;

        let h_code = fixed_width_int_opt_signed(input, &line[42..45])?.unwrap_or(0);
        if !(0..=5).contains(&h_code) {
            return Err(Err::Error(NomError::new(input, NomErrorKind::Verify)));
        }
        let hydrogens = convert_atom_hydrogen_count_code(h_code as u8, false)
            .map_err(|_| Err::Error(NomError::new(input, NomErrorKind::Verify)))?;

        if !is_all_whitespace_or_zeroes(&line[45..48]) {
            return Err(Err::Error(NomError::new(input, NomErrorKind::Verify)));
        }

        let valence_code = fixed_width_int_opt_signed(input, &line[48..51])?.unwrap_or(0);
        if !(0..=15).contains(&valence_code) {
            return Err(Err::Error(NomError::new(input, NomErrorKind::Verify)));
        }
        let valence = convert_atom_valence_code(valence_code as u8)
            .map_err(|_| Err::Error(NomError::new(input, NomErrorKind::Verify)))?;

        let class = fixed_width_int_opt_unsigned(input, &line[60..63])?
            .and_then(|val| (1..=999).contains(&val).then_some(val));

        for i in 0..2 {
            let start = i * 3;
            let end = start + 3;
            if !is_all_whitespace_or_zeroes(&line[63 + start..63 + end]) {
                return Err(Err::Error(NomError::new(input, NomErrorKind::Verify)));
            }
        }

        Ok((
            remaining,
            (
                Atom {
                    element,
                    charge,
                    isotope_mass,
                    hydrogens,
                    implicit_h: false,
                    valence,
                    unpaired_e,
                    aromatic: None,
                    chirality,
                    class,
                    span: None,
                    label: None,
                    value: None,
                },
                position,
            ),
        ))
    }
}

#[cfg(test)]
mod tests;
