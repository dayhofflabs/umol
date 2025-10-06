//! SMILES parser

use strum::{AsRefStr, EnumDiscriminants, EnumIter, IntoEnumIterator};
use umol_data::Element;

use super::api::ParseMeta;
use crate::io::config::SmilesParseFlags;
use crate::io::ir::builder::{AtomData, BondData, MoleculeBuilder};
use crate::io::ir::{BondDir, BondOrder, Chirality, Molecule};
use crate::io::smiles::diagnostics::DiagnosticCode;

#[derive(Debug, Clone, PartialEq, EnumIter, AsRefStr, EnumDiscriminants)]
#[strum(serialize_all = "SCREAMING_SNAKE_CASE")]
pub enum ParseError {
    InvalidWhitespace { pos: usize },
    InvalidComment { pos: usize },
    UnterminatedBlockComment { pos: usize },
    InvalidElement { pos: usize },
    InvalidToken { pos: usize },

    UnbalancedOpenParen { pos: usize },
    UnbalancedCloseParen { pos: usize },
    EmptyBranch { pos: usize },
    EmptyGroup { pos: usize },
    NonfinalGroup { pos: usize },

    LeadingBond { pos: usize },
    TrailingBond { pos: usize },
    ConsecutiveBonds { pos: usize },

    LeadingRing { pos: usize },
    UnbalancedRingIndex { open_pos: usize },
    InvalidRingIndex { pos: usize },
    MismatchedRingBondDirs { pos: usize, open_pos: usize },
    MismatchedRingBondOrders { pos: usize, open_pos: usize },

    LeadingDot { pos: usize },
    TrailingDot { pos: usize },
    ConsecutiveDots { pos: usize },
    DotBeforeRing { pos: usize },

    EmptyBracket { pos: usize },
    UnbalancedOpenBracket { pos: usize },
    UnbalancedCloseBracket { pos: usize },
    StrayBracketField { pos: usize },
    DuplicateBracketField { pos: usize },
    MissingClassIndex { pos: usize },
    MissingChiralityIndex { pos: usize },
    ChiralityOutOfRange { pos: usize },
    BracketHwithHcount { pos: usize },
    InvalidBracket { pos: usize },
}

impl ParseError {
    pub fn all() -> impl Iterator<Item = ParseError> {
        ParseError::iter()
    }
    pub fn as_str(&self) -> &str {
        self.as_ref()
    }
}

impl From<ParseErrorDiscriminants> for DiagnosticCode {
    fn from(discriminant: ParseErrorDiscriminants) -> Self {
        match discriminant {
            ParseErrorDiscriminants::InvalidWhitespace => DiagnosticCode::InvalidWhitespace,
            ParseErrorDiscriminants::InvalidComment => DiagnosticCode::InvalidComment,
            ParseErrorDiscriminants::UnterminatedBlockComment => {
                DiagnosticCode::UnterminatedBlockComment
            }
            ParseErrorDiscriminants::InvalidElement => DiagnosticCode::InvalidElement,
            ParseErrorDiscriminants::InvalidToken => DiagnosticCode::InvalidToken,
            ParseErrorDiscriminants::UnbalancedOpenParen => DiagnosticCode::UnbalancedOpenParen,
            ParseErrorDiscriminants::UnbalancedCloseParen => DiagnosticCode::UnbalancedCloseParen,
            ParseErrorDiscriminants::EmptyBranch => DiagnosticCode::EmptyBranch,
            ParseErrorDiscriminants::EmptyGroup => DiagnosticCode::EmptyGroup,
            ParseErrorDiscriminants::NonfinalGroup => DiagnosticCode::NonfinalGroup,
            ParseErrorDiscriminants::LeadingBond => DiagnosticCode::LeadingBond,
            ParseErrorDiscriminants::TrailingBond => DiagnosticCode::TrailingBond,
            ParseErrorDiscriminants::ConsecutiveBonds => DiagnosticCode::ConsecutiveBonds,
            ParseErrorDiscriminants::LeadingRing => DiagnosticCode::LeadingRing,
            ParseErrorDiscriminants::UnbalancedRingIndex => DiagnosticCode::UnbalancedRingIndex,
            ParseErrorDiscriminants::InvalidRingIndex => DiagnosticCode::InvalidRingIndex,
            ParseErrorDiscriminants::MismatchedRingBondDirs => {
                DiagnosticCode::MismatchedRingBondDirs
            }
            ParseErrorDiscriminants::MismatchedRingBondOrders => {
                DiagnosticCode::MismatchedRingBondOrders
            }
            ParseErrorDiscriminants::LeadingDot => DiagnosticCode::LeadingDot,
            ParseErrorDiscriminants::TrailingDot => DiagnosticCode::TrailingDot,
            ParseErrorDiscriminants::ConsecutiveDots => DiagnosticCode::ConsecutiveDots,
            ParseErrorDiscriminants::DotBeforeRing => DiagnosticCode::DotBeforeRing,
            ParseErrorDiscriminants::EmptyBracket => DiagnosticCode::EmptyBracket,
            ParseErrorDiscriminants::UnbalancedOpenBracket => DiagnosticCode::UnbalancedOpenBracket,
            ParseErrorDiscriminants::UnbalancedCloseBracket => {
                DiagnosticCode::UnbalancedCloseBracket
            }
            ParseErrorDiscriminants::StrayBracketField => DiagnosticCode::StrayBracketField,
            ParseErrorDiscriminants::DuplicateBracketField => DiagnosticCode::DuplicateBracketField,
            ParseErrorDiscriminants::MissingClassIndex => DiagnosticCode::MissingClassIndex,
            ParseErrorDiscriminants::MissingChiralityIndex => DiagnosticCode::MissingChiralityIndex,
            ParseErrorDiscriminants::ChiralityOutOfRange => DiagnosticCode::ChiralityOutOfRange,
            ParseErrorDiscriminants::BracketHwithHcount => DiagnosticCode::BracketHwithHcount,
            ParseErrorDiscriminants::InvalidBracket => DiagnosticCode::InvalidBracket,
        }
    }
}

// Public entrypoint: strict OpenSMILES
pub fn parse_smiles(input: &[u8]) -> Result<Molecule, ParseError> {
    let flags = SmilesParseFlags::STRICT_OPENSMILES;
    parse_smiles_with(input, flags)
}

// Flags-aware inner parser
pub fn parse_smiles_with(input: &[u8], flags: SmilesParseFlags) -> Result<Molecule, ParseError> {
    let allow_ws = flags.contains(SmilesParseFlags::INTERTOKEN_WS);
    let allow_comments = flags.contains(SmilesParseFlags::COMMENTS);
    let use_eoi = flags.contains(SmilesParseFlags::EXPLICIT_EOI);

    let input = if use_eoi {
        let cut = truncate_at_eoi(input, allow_comments);
        &input[..cut]
    } else {
        input
    };

    if !allow_ws && !allow_comments {
        let mut start = 0usize;
        while start < input.len() && matches!(input[start], b' ' | b'\t' | b'\n' | b'\r') {
            start += 1;
        }
        if start == input.len() {
            return Ok(Molecule::default());
        }
        if start > 0 {
            return Err(ParseError::InvalidWhitespace { pos: 0 });
        }
        let mut end = input.len();
        while end > 0 && matches!(input[end - 1], b' ' | b'\t' | b'\n' | b'\r') {
            end -= 1;
        }
        for (k, b) in input[start..end].iter().enumerate() {
            if matches!(*b, b' ' | b'\t' | b'\n' | b'\r') {
                return Err(ParseError::InvalidWhitespace { pos: start + k });
            }
        }
        return parse_smiles_inner(&input[start..end], flags).map(|(m, _)| m);
    }

    parse_smiles_inner(input, flags).map(|(m, _)| m)
}

// Fixed-size resources for parser internals
const RING_TABLE_LEN: usize = 100; // OpenSMILES ring indices: 0..9 and %00..%99
const BRANCH_STACK_DEPTH: usize = 16; // Branch stack depth
const RING_SEQUENCE_CAPACITY: usize = 8; // Ring sequence capacity

#[derive(Debug, Clone, Copy)]
struct OpenRing {
    atom_id: u32,
    order: Option<BondOrder>,
    dir: Option<BondDir>,
    open_pos: usize,
    open_aromatic: bool,
}

#[derive(Debug, Clone, Copy)]
enum Frame {
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
fn is_digit(b: u8) -> bool {
    (b'0'..=b'9').contains(&b)
}

#[inline]
fn skip_line_comment(input: &[u8], mut i: usize) -> usize {
    while i < input.len() && input[i] != b'\n' && input[i] != b'\r' {
        i += 1;
    }
    i
}

#[inline]
fn skip_block_comment(input: &[u8], mut i: usize, start_pos: usize) -> Result<usize, usize> {
    while i + 1 < input.len() {
        if input[i] == b'*' && input[i + 1] == b'/' {
            return Ok(i + 2);
        }
        i += 1;
    }
    Err(start_pos)
}

#[inline]
fn parse_ring_index(input: &[u8], i: usize) -> Result<Option<(usize, usize, bool)>, ParseError> {
    let n = input.len();
    if i >= n {
        return Ok(None);
    }
    let b0 = input[i];
    if is_digit(b0) {
        let idx = (b0 - b'0') as usize;
        return Ok(Some((idx, i + 1, false)));
    }
    if b0 == b'%' {
        if i + 2 >= n || !is_digit(input[i + 1]) || !is_digit(input[i + 2]) {
            return Err(ParseError::InvalidRingIndex { pos: i });
        }
        let idx = ((input[i + 1] - b'0') as usize) * 10 + (input[i + 2] - b'0') as usize;
        return Ok(Some((idx, i + 3, true)));
    }
    Ok(None)
}

#[inline]
fn process_ring_closure(
    ring_table: &mut [Option<OpenRing>; RING_TABLE_LEN],
    builder: &mut MoleculeBuilder,
    last_aromatic: bool,
    last_atom_idx: u32,
    idx: usize,
    order_opt: Option<BondOrder>,
    dir_opt: Option<BondDir>,
    pos: usize,
) -> Result<(), ParseError> {
    let entry = &mut ring_table[idx];
    match entry.take() {
        None => {
            *entry = Some(OpenRing {
                atom_id: last_atom_idx,
                order: order_opt,
                dir: dir_opt,
                open_pos: pos,
                open_aromatic: last_aromatic,
            });
        }
        Some(open) => {
            if let (Some(d1), Some(d2)) = (open.dir, dir_opt) {
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
            if open.dir.is_some() || dir_opt.is_some() {
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
            let final_dir = open.dir.or(dir_opt);
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
                    dir: final_dir,
                },
            );
        }
    }
    Ok(())
}

#[inline]
fn invalid_ring_context(pstack: &Vec<Frame>) -> bool {
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
fn parse_organic_aliphatic_element(input: &[u8], i: usize) -> Option<(Element, usize)> {
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
fn parse_organic_aromatic_element(input: &[u8], i: usize) -> Option<(Element, usize)> {
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
fn parse_bracket_aliphatic_element(inner: &[u8], i: usize) -> Option<(Element, usize)> {
    // Only allow uppercase-starting symbols for aliphatic branch
    let n = inner.len();
    if i >= n || !inner[i].is_ascii_uppercase() {
        return None;
    }
    if i + 1 < n && inner[i + 1].is_ascii_lowercase() {
        if let Some(e) = Element::from_symbol_bytes(&inner[i..i + 2]) {
            return Some((e, 2));
        }
    }
    if let Some(e) = Element::from_symbol_bytes(&inner[i..i + 1]) {
        return Some((e, 1));
    }
    None
}

#[inline]
fn parse_bracket_aromatic_element(inner: &[u8], i: usize) -> Option<(Element, usize)> {
    let n = inner.len();
    if i >= n {
        return None;
    }
    match inner[i] {
        b'b' => Some((Element::B, 1)),
        b'c' => Some((Element::C, 1)),
        b'n' => Some((Element::N, 1)),
        b'o' => Some((Element::O, 1)),
        b'p' => Some((Element::P, 1)),
        b's' => {
            if i + 1 < n && inner[i + 1] == b'e' {
                Some((Element::Se, 2))
            } else {
                Some((Element::S, 1))
            }
        }
        b'a' => {
            if i + 1 < n && inner[i + 1] == b's' {
                Some((Element::As, 2))
            } else {
                None
            }
        }
        _ => None,
    }
}

#[inline]
fn parse_u32(input: &[u8], mut i: usize, max_digits: usize) -> (u32, usize, usize) {
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
fn parse_charge(inner: &[u8], i: usize, sign_char: u8) -> (i32, usize) {
    let n = inner.len();
    let sign = if sign_char == b'+' { 1 } else { -1 };
    let j = i + 1;
    if j < n && inner[j] == sign_char {
        return (2 * sign, j + 1);
    }
    let (mut val, j2, cnt) = parse_u32(inner, j, 2);
    if cnt == 0 {
        val = 1;
    }
    (val as i32 * sign, j2)
}

#[inline]
fn parse_class_index(inner: &[u8], i: usize, pos_base: usize) -> Result<(u32, usize), ParseError> {
    let n = inner.len();
    if i + 1 >= n || !inner[i + 1].is_ascii_digit() {
        return Err(ParseError::MissingClassIndex {
            pos: pos_in_bracket(pos_base, i),
        });
    }
    let (v, j, _) = parse_u32(inner, i + 1, 10);
    Ok((v, j))
}

#[inline]
fn parse_chirality(
    inner: &[u8],
    i: usize,
    pos_base: usize,
) -> Result<(Option<Chirality>, usize), ParseError> {
    let n = inner.len();
    let k = i;
    if k + 1 < n && inner[k + 1] == b'@' {
        return Ok((Some(Chirality::CounterClockwise), k + 2));
    }
    if k + 2 < n && inner[k + 1] == b'T' && inner[k + 2] == b'H' {
        if k + 3 >= n || !inner[k + 3].is_ascii_digit() {
            return Err(ParseError::MissingChiralityIndex {
                pos: pos_in_bracket(pos_base, k),
            });
        }
        let v = (inner[k + 3] - b'0') as u32;
        if v == 1 || v == 2 {
            return Ok((Some(Chirality::Tetrahedral { arr: v }), k + 4));
        }
        return Err(ParseError::ChiralityOutOfRange {
            pos: pos_in_bracket(pos_base, k),
        });
    }
    if k + 2 < n && inner[k + 1] == b'A' && inner[k + 2] == b'L' {
        if k + 3 >= n || !inner[k + 3].is_ascii_digit() {
            return Err(ParseError::MissingChiralityIndex {
                pos: pos_in_bracket(pos_base, k),
            });
        }
        let v = (inner[k + 3] - b'0') as u32;
        if v == 1 || v == 2 {
            return Ok((Some(Chirality::Allenal { arr: v }), k + 4));
        }
        return Err(ParseError::ChiralityOutOfRange {
            pos: pos_in_bracket(pos_base, k),
        });
    }
    if k + 2 < n && inner[k + 1] == b'S' && inner[k + 2] == b'P' {
        if k + 3 >= n || !inner[k + 3].is_ascii_digit() {
            return Err(ParseError::MissingChiralityIndex {
                pos: pos_in_bracket(pos_base, k),
            });
        }
        let v = (inner[k + 3] - b'0') as u32;
        if (1..=3).contains(&v) {
            return Ok((Some(Chirality::SquarePlanar { arr: v }), k + 4));
        }
        return Err(ParseError::ChiralityOutOfRange {
            pos: pos_in_bracket(pos_base, k),
        });
    }
    if k + 2 < n && inner[k + 1] == b'T' && inner[k + 2] == b'B' {
        let (v, j, cnt) = parse_u32(inner, k + 3, 2);
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
    if k + 2 < n && inner[k + 1] == b'O' && inner[k + 2] == b'H' {
        let (v, j, cnt) = parse_u32(inner, k + 3, 2);
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
fn pos_in_bracket(base: usize, local: usize) -> usize {
    base + 1 + local
}

#[inline]
fn attach_atom(
    builder: &mut MoleculeBuilder,
    last_atom_idx: Option<u32>,
    curr_atom_idx: u32,
    pending_bond: &mut Option<(BondOrder, Option<BondDir>, usize)>,
    last_aromatic: bool,
    curr_aromatic: bool,
) {
    if let Some(last) = last_atom_idx {
        if let Some((order, dir, _)) = pending_bond.take() {
            builder.on_bond(last, curr_atom_idx, BondData { order, dir });
        } else if last_aromatic && curr_aromatic {
            builder.on_bond(
                last,
                curr_atom_idx,
                BondData {
                    order: BondOrder::Aromatic,
                    dir: None,
                },
            );
        } else {
            builder.on_bond_single_fast(last, curr_atom_idx);
        }
    }
}

#[inline]
fn parse_bond(b: u8) -> (BondOrder, Option<BondDir>) {
    match b {
        b'-' => (BondOrder::Single, None),
        b'=' => (BondOrder::Double, None),
        b'#' => (BondOrder::Triple, None),
        b'$' => (BondOrder::Quadruple, None),
        b':' => (BondOrder::Aromatic, None),
        b'/' => (BondOrder::Single, Some(BondDir::Up)),
        b'\\' => (BondOrder::Single, Some(BondDir::Down)),
        _ => (BondOrder::Single, None),
    }
}

#[inline]
fn parse_bracket(
    inner: &[u8],
    pos_offset: usize,
) -> Result<
    (
        Option<Element>,
        bool,
        Option<u32>,
        Option<i32>,
        Option<u32>,
        Option<u32>,
        Option<Chirality>,
        bool,
    ),
    ParseError,
> {
    let n = inner.len();
    let mut i = 0usize;

    // 1) Optional isotope (one or more digits)
    let mut isotope: Option<u32> = None;
    let start_digits = i;
    while i < n && inner[i].is_ascii_digit() {
        i += 1;
    }
    if i > start_digits {
        let mut v: u32 = 0;
        for &b in &inner[start_digits..i] {
            v = v.saturating_mul(10).saturating_add((b - b'0') as u32);
        }
        isotope = Some(v);
    }

    // 2) Element symbol or wildcard '*'
    let element: Option<Element>;
    let mut aromatic = false;
    let mut unknown_symbol = false;
    if i < n && inner[i] == b'*' {
        element = None;
        unknown_symbol = true;
        i += 1;
    } else if i < n && inner[i].is_ascii_alphabetic() {
        if let Some((e, consumed)) = parse_bracket_aliphatic_element(inner, i) {
            element = Some(e);
            i += consumed;
            aromatic = false;
        } else if let Some((e, consumed)) = parse_bracket_aromatic_element(inner, i) {
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
    let mut charge: Option<i32> = None;
    let mut class_num: Option<u32> = None;
    let mut hcount: Option<u32> = None;
    let mut chir: Option<Chirality> = None;

    while i < n {
        let b0 = inner[i];
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
                let mut val: u32 = 1; // default H
                if i + 1 < n && inner[i + 1].is_ascii_digit() {
                    val = (inner[i + 1] - b'0') as u32;
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
                let (val, j2) = parse_charge(inner, i, b0);
                charge = Some(val);
                i = j2;
            }
            b':' => {
                if class_num.is_some() {
                    return Err(ParseError::DuplicateBracketField {
                        pos: pos_offset + 1 + i,
                    });
                }
                let (v, j2) = parse_class_index(inner, i, pos_offset)?;
                class_num = Some(v);
                i = j2;
            }
            b'@' => {
                if chir.is_some() {
                    return Err(ParseError::DuplicateBracketField {
                        pos: pos_offset + 1 + i,
                    });
                }
                let (chir_opt, j2) = parse_chirality(inner, i, pos_offset)?;
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

    Ok((
        element,
        aromatic,
        isotope,
        charge,
        class_num,
        hcount,
        chir,
        unknown_symbol,
    ))
}

#[inline]
fn truncate_at_eoi(input: &[u8], allow_comments: bool) -> usize {
    let n = input.len();
    let mut i = 0usize;
    let mut line_start = 0usize;
    let mut had_content = false;
    while i < n {
        let b0 = input[i];
        if b0 == b' ' || b0 == b'\t' {
            i += 1;
            continue;
        }
        if allow_comments && b0 == b'/' && i + 1 < n && input[i + 1] == b'/' {
            i = skip_line_comment(input, i + 2);
        }
        if allow_comments && b0 == b'/' && i + 1 < n && input[i + 1] == b'*' {
            match skip_block_comment(input, i + 2, i) {
                Ok(next) => {
                    i = next;
                }
                Err(_) => {
                    i = n;
                }
            }
            continue;
        }
        if b0 == b'\r' {
            if !had_content {
                return line_start;
            }
            i += 1;
            if i < n && input[i] == b'\n' {
                i += 1;
            }
            line_start = i;
            had_content = false;
            continue;
        }
        if b0 == b'\n' {
            if !had_content {
                return line_start;
            }
            i += 1;
            line_start = i;
            had_content = false;
            continue;
        }
        had_content = true;
        i += 1;
    }
    n
}

fn parse_smiles_inner(
    input: &[u8],
    flags: SmilesParseFlags,
) -> Result<(Molecule, Option<ParseMeta>), ParseError> {
    let allow_ws = flags.contains(SmilesParseFlags::INTERTOKEN_WS);
    let allow_comments = flags.contains(SmilesParseFlags::COMMENTS);
    let no_lints = flags.contains(SmilesParseFlags::NO_LINTS);

    let mut i = 0usize;
    let n = input.len();
    let mut builder = MoleculeBuilder::with_capacity(n.max(1), n.max(1).saturating_sub(1));
    let mut branch_stack: Vec<Frame> = Vec::with_capacity(BRANCH_STACK_DEPTH);
    let mut ring_table: [Option<OpenRing>; RING_TABLE_LEN] = [None; RING_TABLE_LEN];
    let mut ring_sequence: Option<Vec<(u32, usize)>> = if !no_lints {
        Some(Vec::with_capacity(RING_SEQUENCE_CAPACITY))
    } else {
        None
    };
    let mut last_atom_idx: Option<u32> = None;
    let mut pending_bond: Option<(BondOrder, Option<BondDir>, usize)> = None;
    let mut last_aromatic: bool = false;
    let mut just_closed_group: bool = false;
    let mut first_ring_digit: Option<u32> = None;

    while i < n {
        let b0 = input[i];
        if b0 == b'/' && i + 1 < n {
            let b1 = input[i + 1];
            if allow_comments {
                if b1 == b'/' {
                    i += 2;
                    while i < n && input[i] != b'\n' {
                        i += 1;
                    }
                    continue;
                }
                if b1 == b'*' {
                    let start_pos = i;
                    i += 2;
                    let mut closed = false;
                    while i + 1 < n {
                        if input[i] == b'*' && input[i + 1] == b'/' {
                            i += 2;
                            closed = true;
                            break;
                        }
                        i += 1;
                    }
                    if !closed {
                        return Err(ParseError::UnterminatedBlockComment { pos: start_pos });
                    }
                    continue;
                }
            } else {
                if b1 == b'/' || b1 == b'*' {
                    return Err(ParseError::InvalidComment { pos: i });
                }
            }
        }
        if matches!(b0, b' ' | b'\t' | b'\n' | b'\r') {
            if allow_ws {
                i += 1;
                continue;
            }
            return Err(ParseError::InvalidWhitespace { pos: i });
        }
        if b0 != b'(' {
            just_closed_group = false;
        }
        if b0 == b'(' {
            if let Some((_, _, pos)) = pending_bond {
                return Err(ParseError::TrailingBond { pos });
            }
            if just_closed_group {
                last_atom_idx = None;
                branch_stack.push(Frame::Group {
                    had_atom: false,
                    open_pos: i,
                });
                just_closed_group = false;
            } else {
                match last_atom_idx {
                    Some(idx) => branch_stack.push(Frame::Branch {
                        base: idx,
                        had_atom: false,
                        open_pos: i,
                    }),
                    None => branch_stack.push(Frame::Group {
                        had_atom: false,
                        open_pos: i,
                    }),
                }
            }
            i += 1;
            continue;
        }
        if b0 == b')' {
            if let Some((_, _, pos)) = pending_bond {
                return Err(ParseError::TrailingBond { pos });
            }
            let Some(frame) = branch_stack.pop() else {
                return Err(ParseError::UnbalancedCloseParen { pos: i });
            };
            match frame {
                Frame::Branch { base, had_atom, .. } => {
                    if !had_atom {
                        return Err(ParseError::EmptyBranch { pos: i });
                    }
                    last_atom_idx = Some(base);
                }
                Frame::Group {
                    had_atom, open_pos, ..
                } => {
                    if !had_atom {
                        if i + 1 != n {
                            return Err(ParseError::EmptyGroup { pos: i });
                        }
                        if i > 0 && input[i - 1] == b'.' {
                            return Err(ParseError::LeadingDot { pos: i - 1 });
                        }
                        if open_pos != 0 {
                            return Err(ParseError::EmptyGroup { pos: i });
                        }
                        last_atom_idx = None;
                        just_closed_group = false;
                    } else {
                        just_closed_group = true;
                        if branch_stack.is_empty() && i + 1 != n {
                            let next = input[i + 1];
                            if next != b'.' {
                                return Err(ParseError::NonfinalGroup { pos: i });
                            }
                        }
                    }
                }
            }
            i += 1;
            continue;
        }
        if b0 == b'.' {
            if let Some((_, _, pos)) = pending_bond {
                return Err(ParseError::TrailingBond { pos });
            }
            if i == 0 {
                return Err(ParseError::LeadingDot { pos: i });
            }
            if let Some(Frame::Group {
                had_atom: false, ..
            }) = branch_stack.last()
            {
                return Err(ParseError::LeadingDot { pos: i });
            }
            if i + 1 == n {
                return Err(ParseError::TrailingDot { pos: i });
            }
            if input[i + 1] == b'.' {
                return Err(ParseError::ConsecutiveDots { pos: i });
            }
            // Detect dot before ring (single digit ring index)
            if input[i + 1].is_ascii_digit() {
                return Err(ParseError::DotBeforeRing { pos: i });
            }
            // Detect dot before percent ring index
            if input[i + 1] == b'%' {
                return Err(ParseError::DotBeforeRing { pos: i });
            }
            last_atom_idx = None;
            last_aromatic = false;
            i += 1;
            continue;
        }
        match parse_ring_index(input, i) {
            Ok(Some((idx, next_i, percent))) => {
                if last_atom_idx.is_none() {
                    return Err(ParseError::LeadingRing { pos: i });
                }
                if invalid_ring_context(&branch_stack) {
                    return Err(ParseError::LeadingRing { pos: 0 });
                }
                if let Some(seq) = ring_sequence.as_mut() {
                    let d = if percent {
                        ((input[i + 1] - b'0') as u32) * 10 + (input[i + 2] - b'0') as u32
                    } else {
                        (input[i] - b'0') as u32
                    };
                    if first_ring_digit.is_none() {
                        first_ring_digit = Some(d);
                    }
                    seq.push((d, i));
                }
                let bond = pending_bond.take();
                let (order_opt, dir_opt) = bond.map_or((None, None), |(o, d, _)| (Some(o), d));
                process_ring_closure(
                    &mut ring_table,
                    &mut builder,
                    last_aromatic,
                    last_atom_idx.unwrap(),
                    idx,
                    order_opt,
                    dir_opt,
                    i,
                )?;
                i = next_i;
                continue;
            }
            Err(e) => return Err(e),
            Ok(None) => {}
        }
        // percent branch is handled by parse_ring_index above
        if matches!(b0, b'-' | b'=' | b'#' | b'$' | b':' | b'/' | b'\\') {
            if pending_bond.is_some() {
                return Err(ParseError::ConsecutiveBonds { pos: i });
            }
            if last_atom_idx.is_none() {
                if let Some(Frame::Group {
                    had_atom: false, ..
                }) = branch_stack.last()
                {
                    return Err(ParseError::LeadingBond { pos: i });
                }
                return Err(ParseError::LeadingBond { pos: i });
            }
            let (order, dir) = parse_bond(b0);
            pending_bond = Some((order, dir, i));
            i += 1;
            continue;
        }
        if b0 == b'[' {
            let start = i + 1;
            let mut j = start;
            while j < n && input[j] != b']' {
                j += 1;
            }
            if j >= n {
                return Err(ParseError::UnbalancedOpenBracket { pos: i });
            }
            // Empty bracket []
            if j == start {
                return Err(ParseError::EmptyBracket { pos: i });
            }
            let inner = &input[start..j];
            let (elem_opt, aromatic, iso_opt, charge_opt, class_opt, h_opt, chir_opt, unknown) =
                parse_bracket(inner, i)?;
            let (element, aromatic) = match elem_opt {
                Some(e) => (e, aromatic),
                None => (Element::C, false),
            };
            let atom = AtomData {
                element,
                isotope: iso_opt,
                charge: charge_opt,
                hydrogen_count: h_opt,
                class: class_opt,
                aromatic,
                implicit_h: false,
                chirality: chir_opt,
                unknown_symbol: unknown,
            };
            let curr = builder.on_atom(atom);
            attach_atom(
                &mut builder,
                last_atom_idx,
                curr,
                &mut pending_bond,
                last_aromatic,
                aromatic,
            );
            last_atom_idx = Some(curr);
            last_aromatic = aromatic;
            if let Some(top) = branch_stack.last_mut() {
                match top {
                    Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                        *had_atom = true
                    }
                }
            }
            i = j + 1;
            continue;
        }
        if b0 == b'C' {
            if i + 1 < n && input[i + 1] == b'l' {
                let curr = builder.on_atom_fast(Element::Cl, true, false);
                attach_atom(
                    &mut builder,
                    last_atom_idx,
                    curr,
                    &mut pending_bond,
                    last_aromatic,
                    false,
                );
                last_atom_idx = Some(curr);
                last_aromatic = false;
                if let Some(top) = branch_stack.last_mut() {
                    match top {
                        Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                            *had_atom = true
                        }
                    }
                }
                i += 2;
                continue;
            }
            let curr = builder.on_atom_fast(Element::C, true, false);
            attach_atom(
                &mut builder,
                last_atom_idx,
                curr,
                &mut pending_bond,
                last_aromatic,
                false,
            );
            last_atom_idx = Some(curr);
            last_aromatic = false;
            if let Some(top) = branch_stack.last_mut() {
                match top {
                    Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                        *had_atom = true
                    }
                }
            }
            i += 1;
            continue;
        }
        if b0 == b'B' {
            if i + 1 < n && input[i + 1] == b'r' {
                let curr = builder.on_atom_fast(Element::Br, true, false);
                attach_atom(
                    &mut builder,
                    last_atom_idx,
                    curr,
                    &mut pending_bond,
                    last_aromatic,
                    false,
                );
                last_atom_idx = Some(curr);
                last_aromatic = false;
                if let Some(top) = branch_stack.last_mut() {
                    match top {
                        Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                            *had_atom = true
                        }
                    }
                }
                i += 2;
                continue;
            }
            let curr = builder.on_atom_fast(Element::B, true, false);
            attach_atom(
                &mut builder,
                last_atom_idx,
                curr,
                &mut pending_bond,
                last_aromatic,
                false,
            );
            last_atom_idx = Some(curr);
            last_aromatic = false;
            if let Some(top) = branch_stack.last_mut() {
                match top {
                    Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                        *had_atom = true
                    }
                }
            }
            i += 1;
            continue;
        }
        // Elements
        if b0.is_ascii_alphabetic() {
            if let Some((element, consumed)) = parse_organic_aliphatic_element(input, i) {
                let curr = builder.on_atom_fast(element, true, false);
                attach_atom(
                    &mut builder,
                    last_atom_idx,
                    curr,
                    &mut pending_bond,
                    last_aromatic,
                    false,
                );
                last_atom_idx = Some(curr);
                last_aromatic = false;
                if let Some(top) = branch_stack.last_mut() {
                    match top {
                        Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                            *had_atom = true
                        }
                    }
                }
                i += consumed;
                continue;
            }
            if let Some((element, consumed)) = parse_organic_aromatic_element(input, i) {
                let curr = builder.on_atom_fast(element, true, true);
                attach_atom(
                    &mut builder,
                    last_atom_idx,
                    curr,
                    &mut pending_bond,
                    last_aromatic,
                    true,
                );
                last_atom_idx = Some(curr);
                last_aromatic = true;
                if let Some(top) = branch_stack.last_mut() {
                    match top {
                        Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                            *had_atom = true
                        }
                    }
                }
                i += consumed;
                continue;
            }
            return Err(ParseError::InvalidElement { pos: i });
        }
        if b0 == b'*' {
            let atom = AtomData {
                element: Element::C,
                isotope: Some(0),
                charge: Some(0),
                hydrogen_count: Some(0),
                class: None,
                aromatic: false,
                implicit_h: false,
                chirality: None,
                unknown_symbol: true,
            };
            let curr = builder.on_atom(atom);
            attach_atom(
                &mut builder,
                last_atom_idx,
                curr,
                &mut pending_bond,
                last_aromatic,
                false,
            );
            last_atom_idx = Some(curr);
            last_aromatic = false;
            if let Some(top) = branch_stack.last_mut() {
                match top {
                    Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                        *had_atom = true
                    }
                }
            }
            i += 1;
            continue;
        }
        if b0 == b']' {
            return Err(ParseError::UnbalancedCloseBracket { pos: i });
        }
        // Bracket-only fields outside bracket
        if b0 == b'@' || b0 == b'+' {
            return Err(ParseError::StrayBracketField { pos: i });
        }
        return Err(ParseError::InvalidToken { pos: i });
    }

    if pending_bond.is_some() {
        let (_, _, pos) = pending_bond.unwrap();
        return Err(ParseError::TrailingBond { pos });
    }
    if !branch_stack.is_empty() {
        let pos = match branch_stack.last().unwrap() {
            Frame::Branch { open_pos, .. } | Frame::Group { open_pos, .. } => *open_pos,
        };
        return Err(ParseError::UnbalancedOpenParen { pos });
    }
    let mut last_open: Option<usize> = None;
    for entry in ring_table.iter().flatten() {
        match last_open {
            None => last_open = Some(entry.open_pos),
            Some(p) => {
                if entry.open_pos > p {
                    last_open = Some(entry.open_pos)
                }
            }
        }
    }
    if let Some(pos_open) = last_open {
        return Err(ParseError::UnbalancedRingIndex { open_pos: pos_open });
    }
    let meta_opt = if let Some(seq) = ring_sequence {
        Some(ParseMeta {
            token_spans: Vec::new(),
            ring_events: seq.into_iter().map(|(d, _)| d).collect(),
        })
    } else {
        None
    };
    let mut mols = builder.finish();
    let mol = mols.pop().unwrap_or_default();
    Ok((mol, meta_opt))
}

#[cfg(test)]
mod tests;
