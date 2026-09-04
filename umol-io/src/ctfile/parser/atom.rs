//! Atom block parsers for CTab files.

use std::collections::HashMap;

use bstr::ByteSlice;
use umol_chem::element::Element;
use umol_chem::isotope::NamedIsotope;
use umol_geometric_core::Point3D;
use winnow::error::ErrMode;
use winnow::token::take;
use winnow::{ModalResult, Parser};

use super::convert::{
    convert_atom_charge_code, convert_atom_exact_change_flag_code,
    convert_atom_hydrogen_count_code, convert_atom_inversion_flag_code,
    convert_atom_stereo_care_code, convert_atom_stereo_parity_code, convert_atom_symbol_mass_diff,
    convert_atom_valence_code, convert_extended_atom_symbol_mass_diff,
};
use super::utils::{
    finish_line, input_error_column, is_reserved_atom_symbol, next_line, parse_float_f10_4,
    parse_int_opt, validate_unused_n, Input, InputError, PResult,
};
use crate::ctfile::config::CtabParseFlags;
use crate::ctfile::error::ParseError;
use crate::ctfile::parser::rgroup::rgroup_symbol;
use crate::table_ir::{Atom, AtomList, AtomSymbol, ExtendedAtom, WildcardAtom};

type AtomBlockOutput = (Vec<Atom>, Option<Vec<Point3D>>, u32);
type ExtendedAtomBlockOutput = (Vec<ExtendedAtom>, Option<Vec<Point3D>>, u32);

/// Parse atom block (basic atoms only)
pub(super) fn atom_block(
    input: &mut &[u8],
    atom_count: u32,
    line_offset: u32,
    flags: CtabParseFlags,
) -> ModalResult<AtomBlockOutput, ParseError> {
    let ignore_positions = flags.contains(CtabParseFlags::IGNORE_POSITIONS);
    let mut atoms = Vec::with_capacity(atom_count as usize);
    let mut positions = Vec::with_capacity(if ignore_positions {
        0
    } else {
        atom_count as usize
    });

    for line_index in 0..atom_count {
        let physical_line = line_offset + line_index;
        let mut line = next_line(input).map_err(|_| {
            ErrMode::Cut(ParseError::UnexpectedEof {
                line: physical_line,
                block: "atom",
            })
        })?;
        let result = atom_input(flags).parse_next(&mut line).and_then(|value| {
            finish_line(&mut line)?;
            Ok(value)
        });
        let (atom, position) = result.map_err(|error| {
            ErrMode::Cut(ParseError::InvalidAtomLine {
                line: physical_line,
                col: input_error_column(error, &line),
            })
        })?;
        atoms.push(atom);
        if !ignore_positions {
            positions.push(position);
        }
    }

    let positions = if ignore_positions || (atom_count > 1 && Point3D::all_zero(&positions)) {
        None
    } else {
        Some(positions)
    };
    Ok((atoms, positions, line_offset + atom_count))
}

/// Parse extended atom block
pub(super) fn extended_atom_block(
    input: &mut &[u8],
    atom_count: u32,
    line_offset: u32,
    flags: CtabParseFlags,
) -> ModalResult<ExtendedAtomBlockOutput, ParseError> {
    let ignore_positions = flags.contains(CtabParseFlags::IGNORE_POSITIONS);
    let mut atoms = Vec::with_capacity(atom_count as usize);
    let mut positions = Vec::with_capacity(if ignore_positions {
        0
    } else {
        atom_count as usize
    });

    for line_index in 0..atom_count {
        let physical_line = line_offset + line_index;
        let mut line = next_line(input).map_err(|_| {
            ErrMode::Cut(ParseError::UnexpectedEof {
                line: physical_line,
                block: "atom",
            })
        })?;
        let result = extended_atom_input(flags)
            .parse_next(&mut line)
            .and_then(|value| {
                finish_line(&mut line)?;
                Ok(value)
            });
        let (atom, position) = result.map_err(|error| {
            ErrMode::Cut(ParseError::InvalidAtomLine {
                line: physical_line,
                col: input_error_column(error, &line),
            })
        })?;
        atoms.push(atom);
        if !ignore_positions {
            positions.push(position);
        }
    }

    let positions = if ignore_positions || (atom_count > 1 && Point3D::all_zero(&positions)) {
        None
    } else {
        Some(positions)
    };
    Ok((atoms, positions, line_offset + atom_count))
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
fn atom_input<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<Input<'inp>, (Atom, Point3D), ErrMode<InputError>> + use<'inp> {
    move |input: &mut Input<'inp>| {
        let bytes: &[u8] = input.as_ref();
        if bytes.len() < 32 {
            return Err(ErrMode::Backtrack(InputError {
                column: bytes.len() as u32,
            }));
        }

        let mut offset;

        let skip_unused_fields = flags.contains(CtabParseFlags::SKIP_UNUSED_FIELDS);
        let ignore_positions = flags.contains(CtabParseFlags::IGNORE_POSITIONS);
        let atom_map_hcount_fields = flags.contains(CtabParseFlags::ATOM_MAP_HCOUNT_FIELDS);
        let allow_named_isotopes = flags.contains(CtabParseFlags::NAMED_ISOTOPES);
        let extended_range = flags.contains(CtabParseFlags::EXTENDED_RANGE);

        // x, y, z coordinates (0-29)
        let position = if ignore_positions && skip_unused_fields {
            Point3D::zero()
        } else {
            let x = parse_float_f10_4(&bytes[0..10], 0)?;
            let y = parse_float_f10_4(&bytes[10..20], 10)?;
            let z = parse_float_f10_4(&bytes[20..30], 20)?;
            if ignore_positions {
                Point3D::zero()
            } else {
                Point3D::new(x, y, z)
            }
        };

        // Blank (30)
        if bytes[30] != b' ' {
            return Err(ErrMode::Backtrack(InputError { column: 30 }));
        }

        // Atom symbol (31-33)
        let end = 34.min(bytes.len());
        let symbol = parse_atom_symbol(&bytes[31..end], allow_named_isotopes, 31)?;
        offset = end;

        // Mass difference (34-35)
        let mass_diff_code = if bytes.len() >= 36 {
            offset = 36;
            parse_int_opt::<i8>(&bytes[34..36], 34)?
                .filter(|val| (-3..=4).contains(val))
                .unwrap_or(0)
        } else {
            0
        };
        let (element, isotope_mass) = convert_atom_symbol_mass_diff(&symbol, mass_diff_code);

        // Charge/radical (36-38)
        let (charge, unpaired_count) = if bytes.len() >= 39 {
            offset = 39;
            let val = parse_int_opt::<u8>(&bytes[36..39], 36)?
                .filter(|val| (0..=7).contains(val))
                .unwrap_or(0);
            convert_atom_charge_code(val)
        } else {
            (None, None)
        };
        let unpaired_electrons = unpaired_count;
        let multiplicity = None;

        // Stereo parity (39-41)
        let chirality = if bytes.len() >= 42 {
            offset = 42;
            let val = parse_int_opt::<u8>(&bytes[39..42], 39)?.unwrap_or(0);
            if val > 3 {
                return Err(ErrMode::Backtrack(InputError { column: 39 }));
            }
            convert_atom_stereo_parity_code(val)
        } else {
            None
        };

        // Hydrogen count (42-44)
        let hydrogen_count = if bytes.len() >= 45 {
            offset = 45;
            let val = parse_int_opt::<u8>(&bytes[42..45], 42)?.unwrap_or(0);
            let max_val = if extended_range { 13 } else { 5 };
            if val > max_val {
                return Err(ErrMode::Backtrack(InputError { column: 42 }));
            }
            convert_atom_hydrogen_count_code(val)
        } else {
            None
        };
        let hydrogens = if atom_map_hcount_fields {
            Ok(hydrogen_count)
        } else if hydrogen_count.is_some() {
            Err(ErrMode::Backtrack(InputError { column: 42 }))
        } else {
            Ok(None)
        }?;

        // Stereo care box (45-47) - extended
        if bytes.len() >= 48 {
            offset = 48;
            validate_unused_n(&bytes[45..48], 1, 3, false, 45)?;
        }

        // Valence (48-50)
        let valence = if bytes.len() >= 51 {
            offset = 51;
            let val = parse_int_opt::<u8>(&bytes[48..51], 48)?.unwrap_or(0);
            if val > 15 {
                return Err(ErrMode::Backtrack(InputError { column: 48 }));
            }
            convert_atom_valence_code(val)
        } else {
            None
        };

        // Unused fields (51-59)
        let count = ((bytes.len().saturating_sub(51)) / 3).min(3);
        if count > 0 {
            offset = 51 + count * 3;
            validate_unused_n(&bytes[51..51 + count * 3], count, 3, skip_unused_fields, 51)?;
        }

        // Atom mapping number (60-62)
        let atom_map_num = if bytes.len() >= 63 {
            offset = 63;
            parse_int_opt::<u32>(&bytes[60..63], 60)?.filter(|val| (1..=999).contains(val))
        } else {
            None
        };
        // Verify atom mapping number
        let class = if atom_map_hcount_fields {
            Ok(atom_map_num)
        } else if atom_map_num.is_some() {
            Err(ErrMode::Backtrack(InputError { column: 60 }))
        } else {
            Ok(None)
        }?;

        // Inversion flag and exact change (63-68) - extended
        let count = ((bytes.len().saturating_sub(63)) / 3).min(2);
        if count > 0 {
            offset = 63 + count * 3;
            validate_unused_n(&bytes[63..63 + count * 3], count, 3, false, 63)?;
        }

        let _: &[u8] = take(offset).parse_next(input)?;
        Ok((
            Atom {
                element: Some(element),
                charge,
                isotope_mass,
                implicit_hydrogens: hydrogens,
                valence,
                lone_pairs: None,
                unpaired_electrons,
                multiplicity,
                aromatic: None,
                chirality,
                class,
                span: None,
                label: None,
                value: None,
            },
            position,
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
fn extended_atom_input<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<Input<'inp>, (ExtendedAtom, Point3D), ErrMode<InputError>> + use<'inp> {
    let skip_unused_fields = flags.contains(CtabParseFlags::SKIP_UNUSED_FIELDS);
    let ignore_positions = flags.contains(CtabParseFlags::IGNORE_POSITIONS);
    let extended_range = flags.contains(CtabParseFlags::EXTENDED_RANGE);
    let allow_named_isotopes = flags.contains(CtabParseFlags::NAMED_ISOTOPES);
    let allow_wildcards = flags.contains(CtabParseFlags::WILDCARDS);
    let allow_chemaxon_wildcards = flags.contains(CtabParseFlags::CHEMAXON_WILDCARDS);
    let allow_electrons = flags.contains(CtabParseFlags::ELECTRONS);
    let allow_rgroups = flags.contains(CtabParseFlags::RGROUPS);
    let allow_pseudoatoms = flags.contains(CtabParseFlags::PSEUDOATOMS);

    move |input: &mut Input<'inp>| {
        let bytes: &[u8] = input.as_ref();
        if bytes.len() < 32 {
            return Err(ErrMode::Backtrack(InputError {
                column: bytes.len() as u32,
            }));
        }

        let mut offset;

        // x, y, z coordinates (0-29)
        let position = if ignore_positions && skip_unused_fields {
            Point3D::zero()
        } else {
            let x = parse_float_f10_4(&bytes[0..10], 0)?;
            let y = parse_float_f10_4(&bytes[10..20], 10)?;
            let z = parse_float_f10_4(&bytes[20..30], 20)?;
            if ignore_positions {
                Point3D::zero()
            } else {
                Point3D::new(x, y, z)
            }
        };

        // Blank (30)
        if bytes[30] != b' ' {
            return Err(ErrMode::Backtrack(InputError { column: 30 }));
        }

        // Atom symbol (31-33)
        let end = 34.min(bytes.len());
        let symbol = parse_extended_atom_symbol(
            &bytes[31..end],
            allow_named_isotopes,
            allow_wildcards,
            allow_chemaxon_wildcards,
            allow_electrons,
            allow_rgroups,
            allow_pseudoatoms,
            31,
        )?;
        offset = end;

        // Mass difference (34-35)
        let mass_diff_code = if bytes.len() >= 36 {
            offset = 36;
            parse_int_opt::<i8>(&bytes[34..36], 34)?
                .filter(|val| (-3..=4).contains(val))
                .unwrap_or(0)
        } else {
            0
        };
        let isotope_mass = convert_extended_atom_symbol_mass_diff(&symbol, mass_diff_code);

        // Charge/radical (36-38)
        let (charge, unpaired_count) = if bytes.len() >= 39 {
            offset = 39;
            let val = parse_int_opt::<u8>(&bytes[36..39], 36)?
                .filter(|val| (0..=7).contains(val))
                .unwrap_or(0);
            convert_atom_charge_code(val)
        } else {
            (None, None)
        };
        let unpaired_electrons = unpaired_count;
        let multiplicity = None;

        // Stereo parity (39-41)
        let chirality = if bytes.len() >= 42 {
            offset = 42;
            let val = parse_int_opt::<u8>(&bytes[39..42], 39)?.unwrap_or(0);
            if val > 3 {
                return Err(ErrMode::Backtrack(InputError { column: 39 }));
            }
            convert_atom_stereo_parity_code(val)
        } else {
            None
        };

        // Hydrogen count (42-44)
        let hydrogens = if bytes.len() >= 45 {
            offset = 45;
            let val = parse_int_opt::<u8>(&bytes[42..45], 42)?.unwrap_or(0);
            let max_val = if extended_range { 13 } else { 5 };
            if val > max_val {
                return Err(ErrMode::Backtrack(InputError { column: 42 }));
            }
            convert_atom_hydrogen_count_code(val)
        } else {
            None
        };

        // Stereo care box (45-47)
        let stereo_care = if bytes.len() >= 48 {
            offset = 48;
            let val = parse_int_opt::<u8>(&bytes[45..48], 45)?.unwrap_or(0);
            if val > 1 {
                return Err(ErrMode::Backtrack(InputError { column: 45 }));
            }
            convert_atom_stereo_care_code(val)
        } else {
            None
        };

        // Valence (48-50)
        let valence = if bytes.len() >= 51 {
            offset = 51;
            let val = parse_int_opt::<u8>(&bytes[48..51], 48)?
                .filter(|val| (0..=15).contains(val))
                .unwrap_or(0);
            if val > 15 {
                return Err(ErrMode::Backtrack(InputError { column: 48 }));
            }
            convert_atom_valence_code(val)
        } else {
            None
        };

        // Unused fields (51-59)
        let count = ((bytes.len().saturating_sub(51)) / 3).min(3);
        if count > 0 {
            offset = 51 + count * 3;
            validate_unused_n(&bytes[51..51 + count * 3], count, 3, skip_unused_fields, 51)?;
        }

        // Atom mapping number (60-62)
        let class = if bytes.len() >= 63 {
            offset = 63;
            parse_int_opt::<u32>(&bytes[60..63], 60)?.filter(|val| (1..=999).contains(val))
        } else {
            None
        };

        // Inversion flag (63-65)
        let inversion_retention = if bytes.len() >= 66 {
            offset = 66;
            let val = parse_int_opt::<u8>(&bytes[63..66], 63)?.unwrap_or(0);
            if val > 2 {
                return Err(ErrMode::Backtrack(InputError { column: 63 }));
            }
            convert_atom_inversion_flag_code(val)
        } else {
            None
        };

        // Exact change flag (66-68)
        let exact_change = if bytes.len() >= 69 {
            offset = 69;
            let val = parse_int_opt::<u8>(&bytes[66..69], 66)?.unwrap_or(0);
            if val > 1 {
                return Err(ErrMode::Backtrack(InputError { column: 66 }));
            }
            convert_atom_exact_change_flag_code(val)
        } else {
            None
        };

        let _: &[u8] = take(offset).parse_next(input)?;
        Ok((
            ExtendedAtom {
                symbol,
                charge,
                isotope_mass,
                implicit_hydrogens: hydrogens,
                stereo_care,
                valence,
                lone_pairs: None,
                unpaired_electrons,
                multiplicity,
                inversion_retention,
                exact_change,
                aromatic: None,
                chirality,
                class,
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
                ligand_order: None,
                properties: HashMap::new(),
            },
            position,
        ))
    }
}

/// Parse atom symbol (Element and NamedIsotope only).
///
/// Returns error for extended atom symbols (L, A, Q, *, LP, R#, pseudoatoms).
/// If `allow_named_isotopes` is true, allow named isotopes (D, T).
///
#[inline(always)]
fn parse_atom_symbol(field: &[u8], allow_named_isotopes: bool, column: u32) -> PResult<AtomSymbol> {
    let trimmed = field.trim_ascii();
    if trimmed.is_empty() {
        return Err(ErrMode::Backtrack(InputError { column }));
    }
    Element::from_symbol_bytes(trimmed)
        .map(AtomSymbol::Element)
        .or_else(|| {
            allow_named_isotopes
                .then(|| NamedIsotope::from_symbol_bytes(trimmed))
                .flatten()
                .map(AtomSymbol::NamedIsotope)
        })
        .ok_or(ErrMode::Backtrack(InputError { column }))
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
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn parse_extended_atom_symbol(
    field: &[u8],
    allow_named_isotopes: bool,
    allow_wildcards: bool,
    allow_chemaxon_wildcards: bool,
    allow_electrons: bool,
    allow_rgroups: bool,
    allow_pseudoatoms: bool,
    column: u32,
) -> PResult<AtomSymbol> {
    let s = field.trim_ascii();
    if s.is_empty() {
        return Err(ErrMode::Backtrack(InputError { column }));
    }

    if let Some(element) = Element::from_symbol_bytes(s) {
        return Ok(AtomSymbol::Element(element));
    }
    if let Some(isotope) = allow_named_isotopes
        .then(|| NamedIsotope::from_symbol_bytes(s))
        .flatten()
    {
        return Ok(AtomSymbol::NamedIsotope(isotope));
    }
    if allow_wildcards {
        match s {
            b"A" | b"Q" | b"*" | b"X" | b"M" => {
                if let Some(wildcard) = WildcardAtom::from_symbol_bytes(s) {
                    return Ok(AtomSymbol::WildcardAtom(wildcard));
                }
            }
            b"AH" | b"QH" | b"XH" | b"MH" => {
                if allow_chemaxon_wildcards {
                    if let Some(wildcard) = WildcardAtom::from_symbol_bytes(s) {
                        return Ok(AtomSymbol::WildcardAtom(wildcard));
                    }
                } else {
                    return Err(ErrMode::Backtrack(InputError { column }));
                }
            }
            b"L" => return Ok(AtomSymbol::AtomList(AtomList::empty())),
            _ => {} // Fall through
        }
    }
    if allow_rgroups {
        let mut input = Input::new(s);
        let rgroup = rgroup_symbol(&mut input).ok().filter(|_| input.is_empty());
        if let Some(rgroup) = rgroup {
            return Ok(AtomSymbol::RGroup(rgroup));
        }
    }
    if allow_electrons && s == b"LP" {
        return Ok(AtomSymbol::LonePair);
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
        return Err(ErrMode::Backtrack(InputError { column }));
    }

    if allow_pseudoatoms && s.is_ascii() {
        let s = s.to_str_lossy().into_owned();
        return Ok(AtomSymbol::Pseudoatom(s));
    }
    Err(ErrMode::Backtrack(InputError { column }))
}

#[cfg(test)]
mod tests;
