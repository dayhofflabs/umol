//! Utilities for SMILES parser.

use std::borrow::Cow;
use std::str::from_utf8;

use umol_chem::element::Element;

use super::super::config::SmilesSyntaxFlags;
use super::super::error::ParseError;
use super::builder::{BondData, ExtendedMoleculeBuilder, MoleculeEditor};
use crate::table_ir::atom::Chirality;
use crate::table_ir::{
    AtomSymbol, Bond, BondDirection, BondDonation, BondOrder, ExtendedBond, Span, WildcardAtom,
};

#[derive(Debug, Clone, Copy)]
pub(super) enum Frame {
    Branch {
        base: usize,
        had_atom: bool,
        open_pos: usize,
    },
    Group {
        had_atom: bool,
        open_pos: usize,
    },
}

#[inline]
pub(super) fn parse_ring_index(
    input: &[u8],
    i: usize,
    offset: usize,
) -> Result<Option<(usize, usize, bool)>, ParseError> {
    let n = input.len();
    if i >= n {
        return Ok(None);
    }
    let b0 = input[i];
    if b0.is_ascii_digit() {
        let idx = (b0 - b'0') as usize;
        return Ok(Some((idx, i + 1, false)));
    }
    if b0 == b'%' {
        if i + 2 >= n || !input[i + 1].is_ascii_digit() || !input[i + 2].is_ascii_digit() {
            return Err(ParseError::InvalidRingIndex { pos: offset + i });
        }
        let idx = ((input[i + 1] - b'0') as usize) * 10 + (input[i + 2] - b'0') as usize;
        return Ok(Some((idx, i + 3, true)));
    }
    Ok(None)
}

#[inline]
pub(super) fn invalid_ring_context(pstack: &[Frame]) -> bool {
    matches!(
        pstack.last(),
        Some(
            Frame::Branch {
                had_atom: false,
                ..
            } | Frame::Group {
                had_atom: false,
                ..
            }
        )
    )
}

#[inline]
pub(super) fn parse_organic_aliphatic_element(input: &[u8], i: usize) -> Option<(Element, usize)> {
    let n = input.len();
    if i >= n {
        return None;
    }
    match input[i] {
        b'B' => {
            if i + 1 < n && input[i + 1] == b'r' {
                Some((Element::Br, 2))
            } else {
                Some((Element::B, 1))
            }
        }
        b'C' => {
            if i + 1 < n && input[i + 1] == b'l' {
                Some((Element::Cl, 2))
            } else {
                Some((Element::C, 1))
            }
        }
        b'N' => Some((Element::N, 1)),
        b'O' => Some((Element::O, 1)),
        b'S' => Some((Element::S, 1)),
        b'P' => Some((Element::P, 1)),
        b'F' => Some((Element::F, 1)),
        b'I' => Some((Element::I, 1)),
        _ => None,
    }
}

#[inline]
pub(super) fn parse_organic_aromatic_element(input: &[u8], i: usize) -> Option<(Element, usize)> {
    if i >= input.len() {
        return None;
    }
    match input[i] {
        b'b' => Some((Element::B, 1)),
        b'c' => Some((Element::C, 1)),
        b'n' => Some((Element::N, 1)),
        b'o' => Some((Element::O, 1)),
        b'p' => Some((Element::P, 1)),
        b's' => Some((Element::S, 1)),
        _ => None,
    }
}

#[inline]
pub(super) fn parse_bracket_aliphatic_element(input: &[u8], i: usize) -> Option<(Element, usize)> {
    // Only allow uppercase-starting symbols for aliphatic branch
    let n = input.len();
    if i >= n || !input[i].is_ascii_uppercase() {
        return None;
    }
    if i + 1 < n && input[i + 1].is_ascii_lowercase() {
        if let Some(e) = Element::from_symbol_bytes(&input[i..i + 2]) {
            return Some((e, 2));
        }
    }
    if let Some(e) = Element::from_symbol_bytes(&input[i..i + 1]) {
        return Some((e, 1));
    }
    None
}

#[inline]
pub(super) fn parse_bracket_aromatic_element(
    input: &[u8],
    i: usize,
    extended_aromatics: bool,
) -> Option<(Element, usize)> {
    let n = input.len();
    if i >= n {
        return None;
    }
    match input[i] {
        b'b' => Some((Element::B, 1)),
        b'c' => Some((Element::C, 1)),
        b'n' => Some((Element::N, 1)),
        b'o' => Some((Element::O, 1)),
        b'p' => Some((Element::P, 1)),
        b's' => {
            if i + 1 < n && input[i + 1] == b'e' {
                Some((Element::Se, 2))
            } else if extended_aromatics && i + 1 < n && input[i + 1] == b'i' {
                Some((Element::Si, 2))
            } else {
                Some((Element::S, 1))
            }
        }
        b't' => {
            if extended_aromatics && i + 1 < n && input[i + 1] == b'e' {
                Some((Element::Te, 2))
            } else {
                None
            }
        }
        b'a' => {
            if i + 1 < n && input[i + 1] == b's' {
                Some((Element::As, 2))
            } else {
                None
            }
        }
        _ => None,
    }
}

#[inline]
pub(super) fn parse_u32(input: &[u8], mut i: usize, max_digits: usize) -> (u32, usize, usize) {
    let mut v: u32 = 0;
    let mut cnt = 0usize;
    while i < input.len() && input[i].is_ascii_digit() && cnt < max_digits {
        v = v
            .saturating_mul(10)
            .saturating_add((input[i] - b'0') as u32);
        i += 1;
        cnt += 1;
    }
    (v, i, cnt)
}

#[inline]
pub(super) fn parse_charge(input: &[u8], i: usize, sign_char: u8) -> (i8, usize) {
    let n = input.len();
    let sign: i8 = if sign_char == b'+' { 1 } else { -1 };
    let j = i + 1;
    if j < n && input[j] == sign_char {
        return (2 * sign, j + 1);
    }
    let (mut val, j2, cnt) = parse_u32(input, j, 2);
    if cnt == 0 {
        val = 1;
    }
    (val as i8 * sign, j2)
}

#[inline]
pub(super) fn parse_class_index(
    input: &[u8],
    i: usize,
    pos_base: usize,
) -> Result<(u32, usize), ParseError> {
    let n = input.len();
    if i + 1 >= n || !input[i + 1].is_ascii_digit() {
        return Err(ParseError::MissingClassIndex {
            pos: pos_in_bracket(pos_base, i),
        });
    }
    let (v, j, _) = parse_u32(input, i + 1, 10);
    Ok((v, j))
}

#[inline]
pub(super) fn parse_chirality(
    input: &[u8],
    i: usize,
    pos_base: usize,
) -> Result<(Option<Chirality>, usize), ParseError> {
    let n = input.len();
    let k = i;
    if k + 1 < n && input[k + 1] == b'@' {
        return Ok((Some(Chirality::Clockwise), k + 2));
    }
    if k + 2 < n && input[k + 1] == b'T' && input[k + 2] == b'H' {
        if k + 3 >= n || !input[k + 3].is_ascii_digit() {
            return Err(ParseError::MissingChiralityIndex {
                pos: pos_in_bracket(pos_base, k),
            });
        }
        let v = (input[k + 3] - b'0') as u32;
        if v == 1 || v == 2 {
            return Ok((Some(Chirality::Tetrahedral { arr: v }), k + 4));
        }
        return Err(ParseError::ChiralityOutOfRange {
            pos: pos_in_bracket(pos_base, k),
        });
    }
    if k + 2 < n && input[k + 1] == b'A' && input[k + 2] == b'L' {
        if k + 3 >= n || !input[k + 3].is_ascii_digit() {
            return Err(ParseError::MissingChiralityIndex {
                pos: pos_in_bracket(pos_base, k),
            });
        }
        let v = (input[k + 3] - b'0') as u32;
        if v == 1 || v == 2 {
            return Ok((Some(Chirality::Allenal { arr: v }), k + 4));
        }
        return Err(ParseError::ChiralityOutOfRange {
            pos: pos_in_bracket(pos_base, k),
        });
    }
    if k + 2 < n && input[k + 1] == b'S' && input[k + 2] == b'P' {
        if k + 3 >= n || !input[k + 3].is_ascii_digit() {
            return Err(ParseError::MissingChiralityIndex {
                pos: pos_in_bracket(pos_base, k),
            });
        }
        let v = (input[k + 3] - b'0') as u32;
        if (1..=3).contains(&v) {
            return Ok((Some(Chirality::SquarePlanar { arr: v }), k + 4));
        }
        return Err(ParseError::ChiralityOutOfRange {
            pos: pos_in_bracket(pos_base, k),
        });
    }
    if k + 2 < n && input[k + 1] == b'T' && input[k + 2] == b'B' {
        let (v, j, cnt) = parse_u32(input, k + 3, 2);
        if cnt == 0 {
            return Err(ParseError::MissingChiralityIndex {
                pos: pos_in_bracket(pos_base, k),
            });
        }
        if (1..=20).contains(&v) {
            return Ok((Some(Chirality::TrigonalBipyramidal { arr: v }), j));
        }
        return Err(ParseError::ChiralityOutOfRange {
            pos: pos_in_bracket(pos_base, k),
        });
    }
    if k + 2 < n && input[k + 1] == b'O' && input[k + 2] == b'H' {
        let (v, j, cnt) = parse_u32(input, k + 3, 2);
        if cnt == 0 {
            return Err(ParseError::MissingChiralityIndex {
                pos: pos_in_bracket(pos_base, k),
            });
        }
        if (1..=30).contains(&v) {
            return Ok((Some(Chirality::Octahedral { arr: v }), j));
        }
        return Err(ParseError::ChiralityOutOfRange {
            pos: pos_in_bracket(pos_base, k),
        });
    }
    Ok((Some(Chirality::CounterClockwise), k + 1))
}

#[inline]
pub(super) fn pos_in_bracket(base: usize, local: usize) -> usize {
    base + 1 + local
}

#[inline]
pub(super) fn make_bond(start: usize, end: usize, b: BondData) -> Bond {
    let mut bond = Bond::new(start as u32, end as u32, b.order);
    // AtomPair normalization sorts the atoms; a swap reverses the start-atom
    // viewpoint, so flip both the direction and the donation.
    bond.direction = if start > end {
        b.direction.map(|w| w.flip())
    } else {
        b.direction
    };
    bond.donation = if start > end {
        b.donation.map(|d| d.flip())
    } else {
        b.donation
    };
    bond.span = b.span;
    bond
}

#[inline]
pub(super) fn make_extended_bond(start: usize, end: usize, b: BondData) -> ExtendedBond {
    let mut bond = ExtendedBond::new(start as u32, end as u32, b.order);
    // AtomPair normalization sorts the atoms; a swap reverses the start-atom
    // viewpoint, so flip both the direction and the donation.
    bond.direction = if start > end {
        b.direction.map(|w| w.flip())
    } else {
        b.direction
    };
    bond.donation = if start > end {
        b.donation.map(|d| d.flip())
    } else {
        b.donation
    };
    bond.span = b.span;
    bond
}

#[inline]
#[allow(clippy::too_many_arguments)]
pub(super) fn attach_atom(
    builder: &mut MoleculeEditor,
    last_atom_idx: Option<usize>,
    curr_atom_idx: usize,
    pending_bond: &mut Option<(
        BondOrder,
        Option<BondDirection>,
        Option<BondDonation>,
        usize,
    )>,
    curr_aromatic: bool,
    curr_atom_start: u32,
    curr_atom_end: u32,
) {
    if let Some(last) = last_atom_idx {
        if let Some((order, direction, donation, pos)) = pending_bond.take() {
            builder.on_bond(
                last,
                curr_atom_idx,
                BondData {
                    order,
                    direction,
                    donation,
                    span: Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
                },
            );
        } else if builder.is_aromatic(last) && curr_aromatic {
            builder.on_bond(
                last,
                curr_atom_idx,
                BondData {
                    order: BondOrder::Aromatic,
                    direction: None,
                    donation: None,
                    span: Span::from_bytes_opt(Some(curr_atom_start), Some(curr_atom_end)),
                },
            );
        } else {
            builder.on_bond_single_fast(
                last,
                curr_atom_idx,
                Some(curr_atom_start),
                Some(curr_atom_end),
            );
        };
    }
}

#[inline]
pub(super) fn parse_bond(b: u8) -> (BondOrder, Option<BondDirection>) {
    match b {
        b'-' => (BondOrder::Single, None),
        b'=' => (BondOrder::Double, None),
        b'#' => (BondOrder::Triple, None),
        b'$' => (BondOrder::Quadruple, None),
        b':' => (BondOrder::Aromatic, None),
        b'/' => (BondOrder::Single, Some(BondDirection::Rising)),
        b'\\' => (BondOrder::Single, Some(BondDirection::Falling)),
        b'~' => (BondOrder::Any, None),
        _ => (BondOrder::Single, None),
    }
}

/// Parse extended bond tokens including dative bonds (-> and <-).
/// Returns (order, direction, donation, bytes_consumed).
#[inline]
pub(super) fn parse_extended_bond(
    input: &[u8],
    pos: usize,
) -> (
    BondOrder,
    Option<BondDirection>,
    Option<BondDonation>,
    usize,
) {
    let b0 = input[pos];
    let next = if pos + 1 < input.len() {
        Some(input[pos + 1])
    } else {
        None
    };

    match (b0, next) {
        (b'-', Some(b'>')) => (BondOrder::Single, None, Some(BondDonation::Donating), 2),
        (b'<', Some(b'-')) => (BondOrder::Single, None, Some(BondDonation::Accepting), 2),
        (b'~', _) => (BondOrder::Any, None, None, 1),
        _ => {
            let (order, direction) = parse_bond(b0);
            (order, direction, None, 1)
        }
    }
}

#[inline]
#[allow(clippy::type_complexity)]
pub(super) fn parse_bracket(
    input: &[u8],
    pos_offset: usize,
    flags: SmilesSyntaxFlags,
) -> Result<
    (
        Option<Element>,
        Option<bool>,
        Option<u32>,
        Option<i8>,
        Option<u32>,
        Option<u8>,
        Option<Chirality>,
    ),
    ParseError,
> {
    let n = input.len();
    let mut i = 0usize;

    // 1) Optional isotope (one or more digits)
    let mut isotope: Option<u32> = None;
    let start_digits = i;
    while i < n && input[i].is_ascii_digit() {
        i += 1;
    }
    if i > start_digits {
        let mut v: u32 = 0;
        for &b in &input[start_digits..i] {
            v = v.saturating_mul(10).saturating_add((b - b'0') as u32);
        }
        isotope = Some(v);
    }

    // 2) Element symbol
    let element: Option<Element>;
    let aromatic: Option<bool>;
    if i < n && input[i] == b'*' {
        element = None;
        aromatic = None;
        i += 1;
    } else if i < n && input[i].is_ascii_alphabetic() {
        if let Some((e, consumed)) = parse_bracket_aliphatic_element(input, i) {
            element = Some(e);
            i += consumed;
            aromatic = Some(false);
        } else if let Some((e, consumed)) = parse_bracket_aromatic_element(
            input,
            i,
            flags.contains(SmilesSyntaxFlags::EXTENDED_AROMATICS),
        ) {
            element = Some(e);
            i += consumed;
            aromatic = Some(true);
        } else {
            return Err(ParseError::InvalidBracket {
                pos: pos_offset + 1 + i,
            });
        }
    } else {
        return Err(ParseError::InvalidBracket {
            pos: pos_offset + 1 + i,
        });
    }

    // 3) Tail fields in any order
    let mut charge: Option<i8> = None;
    let mut class: Option<u32> = None;
    let mut hydrogens: Option<u8> = None;
    let mut chirality: Option<Chirality> = None;

    while i < n {
        let b0 = input[i];
        match b0 {
            b'H' => {
                if element == Some(Element::H) {
                    return Err(ParseError::BracketHwithHcount {
                        pos: pos_offset + 1 + i,
                    });
                }
                if hydrogens.is_some() {
                    return Err(ParseError::DuplicateBracketField {
                        pos: pos_offset + 1 + i,
                    });
                }
                let mut val: u8 = 1; // default H
                if i + 1 < n && input[i + 1].is_ascii_digit() {
                    val = input[i + 1] - b'0';
                    i += 1;
                }
                hydrogens = Some(val);
                i += 1;
            }
            b'+' | b'-' => {
                if charge.is_some() {
                    return Err(ParseError::DuplicateBracketField {
                        pos: pos_offset + 1 + i,
                    });
                }
                let (val, j2) = parse_charge(input, i, b0);
                charge = Some(val);
                i = j2;
            }
            b':' => {
                if class.is_some() {
                    return Err(ParseError::DuplicateBracketField {
                        pos: pos_offset + 1 + i,
                    });
                }
                let (v, j2) = parse_class_index(input, i, pos_offset)?;
                class = Some(v);
                i = j2;
            }
            b'@' => {
                if chirality.is_some() {
                    return Err(ParseError::DuplicateBracketField {
                        pos: pos_offset + 1 + i,
                    });
                }
                let (chir_opt, j2) = parse_chirality(input, i, pos_offset)?;
                chirality = chir_opt;
                i = j2;
            }
            _ => {
                return Err(ParseError::InvalidBracket {
                    pos: pos_offset + 1 + i,
                });
            }
        }
    }

    Ok((
        element, aromatic, isotope, charge, class, hydrogens, chirality,
    ))
}

#[inline]
#[allow(clippy::too_many_arguments)]
pub(super) fn attach_extended_atom(
    builder: &mut ExtendedMoleculeBuilder,
    last_atom_idx: Option<usize>,
    curr_atom_idx: usize,
    pending_bond: &mut Option<(
        BondOrder,
        Option<BondDirection>,
        Option<BondDonation>,
        usize,
    )>,
    curr_aromatic: bool,
    curr_atom_start: u32,
    curr_atom_end: u32,
) {
    if let Some(last) = last_atom_idx {
        if let Some((order, direction, donation, pos)) = pending_bond.take() {
            builder.on_bond(
                last,
                curr_atom_idx,
                BondData {
                    order,
                    direction,
                    donation,
                    span: Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
                },
            );
        } else if builder.is_aromatic(last) && curr_aromatic {
            builder.on_bond(
                last,
                curr_atom_idx,
                BondData {
                    order: BondOrder::Aromatic,
                    direction: None,
                    donation: None,
                    span: Span::from_bytes_opt(Some(curr_atom_start), Some(curr_atom_end)),
                },
            );
        } else {
            builder.on_bond_single_fast(
                last,
                curr_atom_idx,
                Some(curr_atom_start),
                Some(curr_atom_end),
            );
        };
    }
}

#[inline]
#[allow(clippy::type_complexity)]
pub(super) fn parse_extended_bracket(
    input: &[u8],
    pos_offset: usize,
    flags: SmilesSyntaxFlags,
) -> Result<
    (
        AtomSymbol,
        bool,
        Option<u32>,
        Option<i8>,
        Option<u32>,
        Option<u8>,
        Option<Chirality>,
    ),
    ParseError,
> {
    let n = input.len();
    let mut i = 0usize;

    let mut isotope: Option<u32> = None;
    let start_digits = i;
    while i < n && input[i].is_ascii_digit() {
        i += 1;
    }
    if i > start_digits {
        let mut v: u32 = 0;
        for &b in &input[start_digits..i] {
            v = v.saturating_mul(10).saturating_add((b - b'0') as u32);
        }
        isotope = Some(v);
    }

    let symbol: AtomSymbol;
    let aromatic: bool;
    if i < n && input[i] == b'*' {
        symbol = AtomSymbol::WildcardAtom(WildcardAtom::Any);
        aromatic = false;
        i += 1;
    } else if i < n && input[i].is_ascii_alphabetic() {
        if let Some((e, consumed)) = parse_bracket_aliphatic_element(input, i) {
            symbol = AtomSymbol::Element(e);
            i += consumed;
            aromatic = false;
        } else if let Some((e, consumed)) = parse_bracket_aromatic_element(
            input,
            i,
            flags.contains(SmilesSyntaxFlags::EXTENDED_AROMATICS),
        ) {
            symbol = AtomSymbol::Element(e);
            i += consumed;
            aromatic = true;
        } else {
            return Err(ParseError::InvalidBracket {
                pos: pos_offset + 1 + i,
            });
        }
    } else {
        return Err(ParseError::InvalidBracket {
            pos: pos_offset + 1 + i,
        });
    }

    let mut charge: Option<i8> = None;
    let mut class: Option<u32> = None;
    let mut hydrogens: Option<u8> = None;
    let mut chirality: Option<Chirality> = None;

    while i < n {
        let b0 = input[i];
        match b0 {
            b'H' => {
                if symbol == AtomSymbol::Element(Element::H) {
                    return Err(ParseError::BracketHwithHcount {
                        pos: pos_offset + 1 + i,
                    });
                }
                if hydrogens.is_some() {
                    return Err(ParseError::DuplicateBracketField {
                        pos: pos_offset + 1 + i,
                    });
                }
                let mut val: u8 = 1;
                if i + 1 < n && input[i + 1].is_ascii_digit() {
                    val = input[i + 1] - b'0';
                    i += 1;
                }
                hydrogens = Some(val);
                i += 1;
            }
            b'+' | b'-' => {
                if charge.is_some() {
                    return Err(ParseError::DuplicateBracketField {
                        pos: pos_offset + 1 + i,
                    });
                }
                let (val, j2) = parse_charge(input, i, b0);
                charge = Some(val);
                i = j2;
            }
            b':' => {
                if class.is_some() {
                    return Err(ParseError::DuplicateBracketField {
                        pos: pos_offset + 1 + i,
                    });
                }
                let (v, j2) = parse_class_index(input, i, pos_offset)?;
                class = Some(v);
                i = j2;
            }
            b'@' => {
                if chirality.is_some() {
                    return Err(ParseError::DuplicateBracketField {
                        pos: pos_offset + 1 + i,
                    });
                }
                let (chir_opt, j2) = parse_chirality(input, i, pos_offset)?;
                chirality = chir_opt;
                i = j2;
            }
            _ => {
                return Err(ParseError::InvalidBracket {
                    pos: pos_offset + 1 + i,
                });
            }
        }
    }

    Ok((
        symbol, aromatic, isotope, charge, class, hydrogens, chirality,
    ))
}

/// Unescape `&#code;` sequences in a byte string.
///
/// CXSMILES uses HTML-style numeric character references to escape special characters.
/// Returns borrowed data when no escapes are present.
pub(crate) fn unescape_html_entities(input: &[u8]) -> Cow<'_, [u8]> {
    if !input.windows(2).any(|w| w == b"&#") {
        return Cow::Borrowed(input);
    }

    let mut result = Vec::with_capacity(input.len());
    let mut i = 0;

    while i < input.len() {
        if i + 2 < input.len() && input[i] == b'&' && input[i + 1] == b'#' {
            let mut j = i + 2;
            while j < input.len() && input[j] != b';' && input[j].is_ascii_digit() {
                j += 1;
            }
            if j < input.len() && input[j] == b';' {
                if let Ok(s) = from_utf8(&input[i + 2..j]) {
                    if let Ok(code) = s.parse::<u8>() {
                        result.push(code);
                        i = j + 1;
                        continue;
                    }
                }
            }
        }
        result.push(input[i]);
        i += 1;
    }

    Cow::Owned(result)
}

/// Split on semicolons while respecting `&#n;` escape sequences.
///
/// In CXSMILES labels, semicolons separate entries, but `&#59;` represents a literal semicolon.
/// Returns slices into the original input.
pub(crate) fn split_escaped_semicolons(input: &[u8]) -> Vec<&[u8]> {
    let mut result = Vec::new();
    let mut start = 0;
    let mut i = 0;

    while i < input.len() {
        if i + 2 < input.len() && input[i] == b'&' && input[i + 1] == b'#' {
            let mut j = i + 2;
            while j < input.len() && input[j].is_ascii_digit() {
                j += 1;
            }
            if j < input.len() && input[j] == b';' {
                i = j + 1;
                continue;
            }
        }

        if input[i] == b';' {
            result.push(&input[start..i]);
            i += 1;
            start = i;
        } else {
            i += 1;
        }
    }

    result.push(&input[start..]);
    result
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::plain(b"hello", b"hello")]
    #[case::semicolon(b"&#59;", b";")]
    #[case::comma(b"&#44;", b",")]
    #[case::mixed(b"a&#59;b", b"a;b")]
    #[case::multiple(b"&#59;&#59;", b";;")]
    #[case::incomplete(b"&#", b"&#")]
    #[case::no_semicolon(b"&#65x", b"&#65x")]
    fn test_unescape_html_entities(#[case] input: &[u8], #[case] expected: &[u8]) {
        let result = unescape_html_entities(input);
        assert_eq!(&*result, expected);
    }

    #[rstest]
    #[case::single(b"abc", vec![&b"abc"[..]])]
    #[case::two(b"a;b", vec![&b"a"[..], &b"b"[..]])]
    #[case::empty_parts(b";", vec![&b""[..], &b""[..]])]
    #[case::escaped(b"a&#59;b", vec![&b"a&#59;b"[..]])]
    #[case::mixed(b"a;b&#59;c;d", vec![&b"a"[..], &b"b&#59;c"[..], &b"d"[..]])]
    fn test_split_escaped_semicolons(#[case] input: &[u8], #[case] expected: Vec<&[u8]>) {
        let result = split_escaped_semicolons(input);
        assert_eq!(result, expected);
    }
}
