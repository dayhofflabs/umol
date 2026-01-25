//! Utilities for SMILES parser.

use umol_data::Element;

use super::super::error::ParseError;
use super::builder::{BondData, ExtendedMoleculeBuilder, MoleculeBuilder};
use crate::span::Span;
use crate::table_ir::{AtomSymbol, BondDirection, BondOrder, Chirality, WildcardAtom};

#[derive(Debug, Clone, Copy)]
pub(super) struct OpenRing {
    pub(super) atom_id: u32,
    pub(super) order: Option<BondOrder>,
    pub(super) direction: Option<BondDirection>,
    pub(super) open_pos: usize,
    pub(super) open_end: usize,
    pub(super) open_aromatic: bool,
}

#[derive(Debug, Clone, Copy)]
pub(super) enum Frame {
    Branch {
        base: u32,
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
            return Err(ParseError::InvalidRingIndex { pos: i });
        }
        let idx = ((input[i + 1] - b'0') as usize) * 10 + (input[i + 2] - b'0') as usize;
        return Ok(Some((idx, i + 3, true)));
    }
    Ok(None)
}

#[inline]
#[allow(clippy::too_many_arguments)]
pub(super) fn process_ring_closure(
    ring_table: &mut Vec<Option<OpenRing>>,
    builder: &mut MoleculeBuilder,
    last_aromatic: bool,
    last_atom_idx: u32,
    idx: usize,
    order_opt: Option<BondOrder>,
    dir_opt: Option<BondDirection>,
    pos: usize,
    token_end: usize,
) -> Result<(), ParseError> {
    if ring_table.len() <= idx {
        ring_table.resize_with(idx + 1, || None);
    }
    let entry = &mut ring_table[idx];
    match entry.take() {
        None => {
            *entry = Some(OpenRing {
                atom_id: last_atom_idx,
                order: order_opt,
                direction: dir_opt,
                open_pos: pos,
                open_end: token_end,
                open_aromatic: last_aromatic,
            });
            builder.on_ring_open(
                idx as u32,
                Some(pos as u32),
                Some(token_end as u32),
                Some(last_atom_idx),
            );
        }
        Some(open) => {
            if let (Some(d1), Some(d2)) = (open.direction, dir_opt) {
                if d1 != d2 {
                    return Err(ParseError::MismatchedRingBondDirs {
                        pos,
                        open_pos: open.open_pos,
                    });
                }
            }
            if let (Some(o1), Some(o2)) = (open.order, order_opt) {
                if o1 != o2 {
                    return Err(ParseError::MismatchedRingBondOrders {
                        pos,
                        open_pos: open.open_pos,
                    });
                }
            }
            if open.direction.is_some() || dir_opt.is_some() {
                let ord = open.order.or(order_opt).unwrap_or(BondOrder::Single);
                if ord != BondOrder::Single {
                    return Err(ParseError::MismatchedRingBondOrders {
                        pos,
                        open_pos: open.open_pos,
                    });
                }
            }
            let mut final_order = match (open.order, order_opt) {
                (Some(o1), Some(o2)) => {
                    if o1 == o2 {
                        o1
                    } else {
                        o2
                    }
                }
                (Some(o), None) | (None, Some(o)) => o,
                (None, None) => BondOrder::Single,
            };
            let final_dir = open.direction.or(dir_opt);
            let a = open.atom_id;
            let b = last_atom_idx;
            if final_order == BondOrder::Single && open.open_aromatic && last_aromatic {
                final_order = BondOrder::Aromatic;
            }
            builder.on_bond(
                a,
                b,
                BondData {
                    order: final_order,
                    direction: final_dir,
                    span: Span::from_bytes_opt(
                        Some(open.open_pos as u32),
                        Some(open.open_end as u32),
                    ),
                },
            );
            builder.on_ring_close(
                idx as u32,
                Some(pos as u32),
                Some(token_end as u32),
                Some(b),
            );
        }
    }
    Ok(())
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
pub(super) fn parse_bracket_aromatic_element(input: &[u8], i: usize) -> Option<(Element, usize)> {
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
            } else {
                Some((Element::S, 1))
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
        return Ok((Some(Chirality::CounterClockwise), k + 2));
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
    Ok((Some(Chirality::Clockwise), k + 1))
}

#[inline]
pub(super) fn pos_in_bracket(base: usize, local: usize) -> usize {
    base + 1 + local
}

#[inline]
#[allow(clippy::too_many_arguments)]
pub(super) fn attach_atom(
    builder: &mut MoleculeBuilder,
    last_atom_idx: Option<u32>,
    curr_atom_idx: u32,
    pending_bond: &mut Option<(BondOrder, Option<BondDirection>, usize)>,
    last_aromatic: bool,
    curr_aromatic: bool,
    curr_atom_start: u32,
    curr_atom_end: u32,
) {
    if let Some(last) = last_atom_idx {
        if let Some((order, bond_dir, pos)) = pending_bond.take() {
            builder.on_bond(
                last,
                curr_atom_idx,
                BondData {
                    order,
                    direction: bond_dir,
                    span: Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
                },
            );
        } else if last_aromatic && curr_aromatic {
            builder.on_bond(
                last,
                curr_atom_idx,
                BondData {
                    order: BondOrder::Aromatic,
                    direction: None,
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
        b'/' => (BondOrder::Single, Some(BondDirection::Up)),
        b'\\' => (BondOrder::Single, Some(BondDirection::Down)),
        _ => (BondOrder::Single, None),
    }
}

#[inline]
#[allow(clippy::type_complexity)]
pub(super) fn parse_bracket(
    input: &[u8],
    pos_offset: usize,
) -> Result<
    (
        Option<Element>,
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
    // Wildcards not supported in basic SMILES parser
    let element: Option<Element>;
    let aromatic: bool;
    if i < n && input[i] == b'*' {
        return Err(ParseError::InvalidBracket {
            pos: pos_offset + 1 + i,
        });
    } else if i < n && input[i].is_ascii_alphabetic() {
        if let Some((e, consumed)) = parse_bracket_aliphatic_element(input, i) {
            element = Some(e);
            i += consumed;
            aromatic = false;
        } else if let Some((e, consumed)) = parse_bracket_aromatic_element(input, i) {
            element = Some(e);
            i += consumed;
            aromatic = true;
        } else {
            return Err(ParseError::InvalidBracket {
                pos: pos_offset + 1 + i,
            });
        }
    } else {
        // Neither '*' nor element
        return Err(ParseError::InvalidBracket {
            pos: pos_offset + 1 + i,
        });
    }

    // 3) Tail fields in any order
    let mut charge: Option<i8> = None;
    let mut class_num: Option<u32> = None;
    let mut hcount: Option<u8> = None;
    let mut chir: Option<Chirality> = None;

    while i < n {
        let b0 = input[i];
        match b0 {
            b'H' => {
                if element == Some(Element::H) {
                    return Err(ParseError::BracketHwithHcount {
                        pos: pos_offset + 1 + i,
                    });
                }
                if hcount.is_some() {
                    return Err(ParseError::DuplicateBracketField {
                        pos: pos_offset + 1 + i,
                    });
                }
                let mut val: u8 = 1; // default H
                if i + 1 < n && input[i + 1].is_ascii_digit() {
                    val = input[i + 1] - b'0';
                    i += 1;
                }
                hcount = Some(val);
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
                if class_num.is_some() {
                    return Err(ParseError::DuplicateBracketField {
                        pos: pos_offset + 1 + i,
                    });
                }
                let (v, j2) = parse_class_index(input, i, pos_offset)?;
                class_num = Some(v);
                i = j2;
            }
            b'@' => {
                if chir.is_some() {
                    return Err(ParseError::DuplicateBracketField {
                        pos: pos_offset + 1 + i,
                    });
                }
                let (chir_opt, j2) = parse_chirality(input, i, pos_offset)?;
                chir = chir_opt;
                i = j2;
            }
            _ => {
                return Err(ParseError::InvalidBracket {
                    pos: pos_offset + 1 + i,
                });
            }
        }
    }

    Ok((element, aromatic, isotope, charge, class_num, hcount, chir))
}

#[inline]
#[allow(clippy::too_many_arguments)]
pub(super) fn process_ring_closure_extended(
    ring_table: &mut Vec<Option<OpenRing>>,
    builder: &mut ExtendedMoleculeBuilder,
    last_aromatic: bool,
    last_atom_idx: u32,
    idx: usize,
    order_opt: Option<BondOrder>,
    dir_opt: Option<BondDirection>,
    pos: usize,
    token_end: usize,
) -> Result<(), ParseError> {
    if ring_table.len() <= idx {
        ring_table.resize_with(idx + 1, || None);
    }
    let entry = &mut ring_table[idx];
    match entry.take() {
        None => {
            *entry = Some(OpenRing {
                atom_id: last_atom_idx,
                order: order_opt,
                direction: dir_opt,
                open_pos: pos,
                open_end: token_end,
                open_aromatic: last_aromatic,
            });
            builder.on_ring_open(
                idx as u32,
                Some(pos as u32),
                Some(token_end as u32),
                Some(last_atom_idx),
            );
        }
        Some(open) => {
            if let (Some(d1), Some(d2)) = (open.direction, dir_opt) {
                if d1 != d2 {
                    return Err(ParseError::MismatchedRingBondDirs {
                        pos,
                        open_pos: open.open_pos,
                    });
                }
            }
            if let (Some(o1), Some(o2)) = (open.order, order_opt) {
                if o1 != o2 {
                    return Err(ParseError::MismatchedRingBondOrders {
                        pos,
                        open_pos: open.open_pos,
                    });
                }
            }
            if open.direction.is_some() || dir_opt.is_some() {
                let ord = open.order.or(order_opt).unwrap_or(BondOrder::Single);
                if ord != BondOrder::Single {
                    return Err(ParseError::MismatchedRingBondOrders {
                        pos,
                        open_pos: open.open_pos,
                    });
                }
            }
            let mut final_order = match (open.order, order_opt) {
                (Some(o1), Some(o2)) => {
                    if o1 == o2 {
                        o1
                    } else {
                        o2
                    }
                }
                (Some(o), None) | (None, Some(o)) => o,
                (None, None) => BondOrder::Single,
            };
            let final_dir = open.direction.or(dir_opt);
            let a = open.atom_id;
            let b = last_atom_idx;
            if final_order == BondOrder::Single && open.open_aromatic && last_aromatic {
                final_order = BondOrder::Aromatic;
            }
            builder.on_bond(
                a,
                b,
                BondData {
                    order: final_order,
                    direction: final_dir,
                    span: Span::from_bytes_opt(
                        Some(open.open_pos as u32),
                        Some(open.open_end as u32),
                    ),
                },
            );
            builder.on_ring_close(
                idx as u32,
                Some(pos as u32),
                Some(token_end as u32),
                Some(b),
            );
        }
    }
    Ok(())
}

#[inline]
#[allow(clippy::too_many_arguments)]
pub(super) fn attach_atom_extended(
    builder: &mut ExtendedMoleculeBuilder,
    last_atom_idx: Option<u32>,
    curr_atom_idx: u32,
    pending_bond: &mut Option<(BondOrder, Option<BondDirection>, usize)>,
    last_aromatic: bool,
    curr_aromatic: bool,
    curr_atom_start: u32,
    curr_atom_end: u32,
) {
    if let Some(last) = last_atom_idx {
        if let Some((order, bond_dir, pos)) = pending_bond.take() {
            builder.on_bond(
                last,
                curr_atom_idx,
                BondData {
                    order,
                    direction: bond_dir,
                    span: Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
                },
            );
        } else if last_aromatic && curr_aromatic {
            builder.on_bond(
                last,
                curr_atom_idx,
                BondData {
                    order: BondOrder::Aromatic,
                    direction: None,
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
pub(super) fn parse_bracket_extended(
    input: &[u8],
    pos_offset: usize,
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
        } else if let Some((e, consumed)) = parse_bracket_aromatic_element(input, i) {
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
    let mut class_num: Option<u32> = None;
    let mut hcount: Option<u8> = None;
    let mut chir: Option<Chirality> = None;

    while i < n {
        let b0 = input[i];
        match b0 {
            b'H' => {
                if symbol == AtomSymbol::Element(Element::H) {
                    return Err(ParseError::BracketHwithHcount {
                        pos: pos_offset + 1 + i,
                    });
                }
                if hcount.is_some() {
                    return Err(ParseError::DuplicateBracketField {
                        pos: pos_offset + 1 + i,
                    });
                }
                let mut val: u8 = 1;
                if i + 1 < n && input[i + 1].is_ascii_digit() {
                    val = input[i + 1] - b'0';
                    i += 1;
                }
                hcount = Some(val);
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
                if class_num.is_some() {
                    return Err(ParseError::DuplicateBracketField {
                        pos: pos_offset + 1 + i,
                    });
                }
                let (v, j2) = parse_class_index(input, i, pos_offset)?;
                class_num = Some(v);
                i = j2;
            }
            b'@' => {
                if chir.is_some() {
                    return Err(ParseError::DuplicateBracketField {
                        pos: pos_offset + 1 + i,
                    });
                }
                let (chir_opt, j2) = parse_chirality(input, i, pos_offset)?;
                chir = chir_opt;
                i = j2;
            }
            _ => {
                return Err(ParseError::InvalidBracket {
                    pos: pos_offset + 1 + i,
                });
            }
        }
    }

    Ok((symbol, aromatic, isotope, charge, class_num, hcount, chir))
}
