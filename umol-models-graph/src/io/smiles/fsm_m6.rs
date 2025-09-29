//! SMILES M6 parser (lex/syntax only, flags-ready)

use smallvec::SmallVec;
use umol_data::Element;

use crate::io::config::SmilesParseFlags;
use crate::io::ir::builder::{AtomData, BondData, MoleculeBuilder};
use crate::io::ir::{BondDir, BondOrder, Molecule};
use crate::io::smiles::parser::utils::BracketField;

#[derive(Debug, Clone, PartialEq)]
pub enum M6Error {
    // Reserved for future WS/comment handling in the M6 loop
    InvalidWhitespace { pos: usize },
    InvalidComment { pos: usize },
    UnterminatedBlockComment { pos: usize },
    UnsupportedToken { pos: usize },
    UnbalancedBranchOpen { pos: usize },
    UnbalancedBranchClose { pos: usize },
    EmptyBranch { pos: usize },
    EmptyGroup { pos: usize },
    TopLevelGroupTrailing { pos: usize },
    TrailingBond { pos: usize },
    ConsecutiveBond { pos: usize },
    LeadingBond { pos: usize },
    RingIndexInvalid { pos: usize },
    LeadingRing { pos: usize },
    RingBondDirConflict { pos: usize, open_pos: usize },
    RingBondOrderConflict { pos: usize, open_pos: usize },
    RingSelfLoop { pos: usize },
    RingTwoMember { pos: usize },
    RingUnclosed { open_pos: usize },
    LeadingDot { pos: usize },
    TrailingDot { pos: usize },
    ConsecutiveDot { pos: usize },
    UnbalancedOpenBracket { pos: usize },
    UnbalancedCloseBracket { pos: usize },
    InvalidBracket { pos: usize },
    BracketHCountTwoDigits { pos: usize },
    BracketEmptyClass { pos: usize },
}

// Public M6 entrypoint: set strict flags and call inner
pub fn parse_smiles_m6(input: &[u8]) -> Result<Molecule, M6Error> {
    let flags = SmilesParseFlags::STRICT_OPENSMILES;
    parse_smiles_inner(input, flags)
}

// Flags-aware inner parser. For now it forwards to M5 without allocation or preprocessing.
pub fn parse_smiles_inner(input: &[u8], flags: SmilesParseFlags) -> Result<Molecule, M6Error> {
    let allow_ws = flags.contains(SmilesParseFlags::INTERTOKEN_WS);
    let allow_comments = flags.contains(SmilesParseFlags::COMMENTS);
    let use_eoi = flags.contains(SmilesParseFlags::EXPLICIT_EOI);

    // If explicit EOI is enabled, truncate at the first blank-line boundary.
    // A blank line contains only ASCII whitespace and, if COMMENTS is enabled,
    // C-style comments; otherwise comments are treated as content.
    let input = if use_eoi {
        let cut = truncate_at_eoi(input, allow_comments);
        &input[..cut]
    } else {
        input
    };

    // In strict mode (no inter-token WS/comments), accept only leading/trailing terminators
    if !allow_ws && !allow_comments {
        // Find first non-terminator
        let mut start = 0usize;
        while start < input.len() && matches!(input[start], b' ' | b'\t' | b'\n' | b'\r') {
            start += 1;
        }
        if start == input.len() {
            return Ok(Molecule::default());
        }
        if start > 0 {
            return Err(M6Error::InvalidWhitespace { pos: 0 });
        }
        // Find last non-terminator
        let mut end = input.len();
        while end > 0 && matches!(input[end - 1], b' ' | b'\t' | b'\n' | b'\r') {
            end -= 1;
        }
        // Ensure no whitespace inside
        for (k, b) in input[start..end].iter().enumerate() {
            if matches!(*b, b' ' | b'\t' | b'\n' | b'\r') {
                return Err(M6Error::InvalidWhitespace { pos: start + k });
            }
        }
        return m6_parse_core(&input[start..end], flags);
    }

    // Extended mode: parse in-loop with WS/comments skipping
    m6_parse_core(input, flags)
}

fn is_digit(b: u8) -> bool {
    (b'0'..=b'9').contains(&b)
}

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

fn map_bond(b: u8) -> (BondOrder, Option<BondDir>) {
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

// Local bracket parser that emits M6Error directly (reusing utils)
fn parse_bracket_m6(
    inner: &[u8],
    pos_offset: usize,
) -> Result<(Option<Element>, Option<u32>, SmallVec<[BracketField; 4]>), M6Error> {
    let bytes = inner;
    if bytes.last() == Some(&b':') {
        return Err(M6Error::BracketEmptyClass { pos: pos_offset });
    }
    if bytes
        .windows(3)
        .any(|w| w[0] == b'H' && w[1].is_ascii_digit() && w[2].is_ascii_digit())
    {
        return Err(M6Error::BracketHCountTwoDigits { pos: pos_offset });
    }
    let s = match std::str::from_utf8(inner) {
        Ok(v) => v,
        Err(_) => return Err(M6Error::InvalidBracket { pos: pos_offset }),
    };
    if !crate::io::smiles::parser::utils::is_valid_bracket_inner(s) {
        return Err(M6Error::InvalidBracket { pos: pos_offset });
    }
    let (elem_opt, iso_opt, tails) = crate::io::smiles::parser::utils::parse_bracket(s);
    Ok((elem_opt, iso_opt, tails))
}

// Zero-allocation scan for the first blank-line boundary. Returns index to cut.
fn truncate_at_eoi(input: &[u8], allow_comments: bool) -> usize {
    let n = input.len();
    let mut i = 0usize;
    let mut line_start = 0usize;
    let mut had_content = false;
    while i < n {
        let b0 = input[i];

        // Horizontal whitespace
        if b0 == b' ' || b0 == b'\t' {
            i += 1;
            continue;
        }

        // Line comment // ... to EOL (only if comments are allowed)
        if allow_comments && b0 == b'/' && i + 1 < n && input[i + 1] == b'/' {
            i += 2;
            while i < n && input[i] != b'\n' && input[i] != b'\r' {
                i += 1;
            }
            // fallthrough to newline handling
        }

        // Block comment /* ... */ (only if comments are allowed)
        if allow_comments && b0 == b'/' && i + 1 < n && input[i + 1] == b'*' {
            i += 2;
            while i + 1 < n {
                if input[i] == b'*' && input[i + 1] == b'/' {
                    i += 2;
                    break;
                }
                i += 1;
            }
            continue; // comments do not add content to the line
        }

        // Newline(s)
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

        // Any other byte is content
        had_content = true;
        i += 1;
    }
    // No blank line found
    n
}

fn m6_parse_core(input: &[u8], flags: SmilesParseFlags) -> Result<Molecule, M6Error> {
    let allow_ws = flags.contains(SmilesParseFlags::INTERTOKEN_WS);
    let allow_comments = flags.contains(SmilesParseFlags::COMMENTS);

    let mut i = 0usize;
    let n = input.len();

    let mut builder = MoleculeBuilder::with_capacity(n.max(1), n.max(1).saturating_sub(1));
    let mut last_atom_idx: Option<u32> = None;
    let mut prev_atom_idx: Option<u32> = None;
    let mut pending_bond: Option<(BondOrder, Option<BondDir>, usize)> = None;
    let mut last_aromatic: bool = false;

    let mut pstack: Vec<Frame> = Vec::new();
    let mut ring_table: [Option<OpenRing>; 100] = [None; 100];
    let mut just_closed_group: bool = false;

    while i < n {
        let b0 = input[i];

        // Comments
        if allow_comments && b0 == b'/' && i + 1 < n {
            let b1 = input[i + 1];
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
                    return Err(M6Error::UnterminatedBlockComment { pos: start_pos });
                }
                continue;
            }
        }

        // Whitespace
        if matches!(b0, b' ' | b'\t' | b'\n' | b'\r') {
            if allow_ws {
                i += 1;
                continue;
            }
            return Err(M6Error::InvalidWhitespace { pos: i });
        }

        if b0 != b'(' {
            just_closed_group = false;
        }

        // Parentheses open
        if b0 == b'(' {
            if let Some((_, _, pos)) = pending_bond {
                return Err(M6Error::TrailingBond { pos });
            }
            if just_closed_group {
                last_atom_idx = None;
                prev_atom_idx = None;
                pstack.push(Frame::Group {
                    had_atom: false,
                    open_pos: i,
                });
                just_closed_group = false;
            } else {
                match last_atom_idx {
                    Some(idx) => pstack.push(Frame::Branch {
                        base: idx,
                        had_atom: false,
                        open_pos: i,
                    }),
                    None => pstack.push(Frame::Group {
                        had_atom: false,
                        open_pos: i,
                    }),
                }
            }
            i += 1;
            continue;
        }

        // Parentheses close
        if b0 == b')' {
            if let Some((_, _, pos)) = pending_bond {
                return Err(M6Error::TrailingBond { pos });
            }
            let Some(frame) = pstack.pop() else {
                return Err(M6Error::UnbalancedBranchClose { pos: i });
            };
            match frame {
                Frame::Branch { base, had_atom, .. } => {
                    if !had_atom {
                        return Err(M6Error::EmptyBranch { pos: i });
                    }
                    last_atom_idx = Some(base);
                    prev_atom_idx = None;
                }
                Frame::Group {
                    had_atom, open_pos, ..
                } => {
                    if !had_atom {
                        if i + 1 != n {
                            return Err(M6Error::EmptyGroup { pos: i });
                        }
                        if i > 0 && input[i - 1] == b'.' {
                            return Err(M6Error::LeadingDot { pos: i });
                        }
                        if open_pos != 0 {
                            return Err(M6Error::EmptyGroup { pos: i });
                        }
                        last_atom_idx = None;
                        prev_atom_idx = None;
                        just_closed_group = false;
                    } else {
                        just_closed_group = true;
                        if pstack.is_empty() && i + 1 != n {
                            let next = input[i + 1];
                            if next != b'.' {
                                return Err(M6Error::TopLevelGroupTrailing { pos: i });
                            }
                        }
                    }
                }
            }
            i += 1;
            continue;
        }

        // Component separator '.'
        if b0 == b'.' {
            if let Some((_, _, pos)) = pending_bond {
                return Err(M6Error::TrailingBond { pos });
            }
            if i == 0 {
                return Err(M6Error::LeadingDot { pos: i });
            }
            if i + 1 == n {
                return Err(M6Error::TrailingDot { pos: i });
            }
            if input[i + 1] == b'.' {
                return Err(M6Error::ConsecutiveDot { pos: i });
            }
            last_atom_idx = None;
            prev_atom_idx = None;
            last_aromatic = false;
            i += 1;
            continue;
        }

        // Ring tokens: digit and %DD
        if is_digit(b0) {
            if last_atom_idx.is_none() {
                return Err(M6Error::LeadingRing { pos: i });
            }
            let idx: usize = (b0 - b'0') as usize;
            let bond = pending_bond.take();
            let (order_opt, dir_opt) = bond.map_or((None, None), |(o, d, _)| (Some(o), d));
            let entry = &mut ring_table[idx];
            match entry.take() {
                None => {
                    *entry = Some(OpenRing {
                        atom_id: last_atom_idx.unwrap(),
                        order: order_opt,
                        dir: dir_opt,
                        open_pos: i,
                        open_aromatic: last_aromatic,
                    });
                }
                Some(open) => {
                    let b = last_atom_idx.unwrap();
                    if open.atom_id == b {
                        return Err(M6Error::RingSelfLoop { pos: i });
                    }
                    if prev_atom_idx == Some(open.atom_id) {
                        return Err(M6Error::RingTwoMember { pos: i });
                    }
                    if let (Some(d1), Some(d2)) = (open.dir, dir_opt) {
                        if d1 != d2 {
                            return Err(M6Error::RingBondDirConflict {
                                pos: i,
                                open_pos: open.open_pos,
                            });
                        }
                    }
                    if let (Some(o1), Some(o2)) = (open.order, order_opt) {
                        if o1 != o2 {
                            return Err(M6Error::RingBondOrderConflict {
                                pos: i,
                                open_pos: open.open_pos,
                            });
                        }
                    }
                    if open.dir.is_some() || dir_opt.is_some() {
                        let ord = open.order.or(order_opt).unwrap_or(BondOrder::Single);
                        if ord != BondOrder::Single {
                            return Err(M6Error::RingBondOrderConflict {
                                pos: i,
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
                    let b = last_atom_idx.unwrap();
                    if final_order == BondOrder::Single && open.open_aromatic && last_aromatic {
                        final_order = BondOrder::Aromatic;
                    }
                    if a != b {
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
            }
            i += 1;
            continue;
        }
        if b0 == b'%' {
            if i + 2 >= n || !is_digit(input[i + 1]) || !is_digit(input[i + 2]) {
                return Err(M6Error::RingIndexInvalid { pos: i });
            }
            if last_atom_idx.is_none() {
                return Err(M6Error::LeadingRing { pos: i });
            }
            let idx: usize = ((input[i + 1] - b'0') as usize) * 10 + (input[i + 2] - b'0') as usize;
            let bond = pending_bond.take();
            let (order_opt, dir_opt) = bond.map_or((None, None), |(o, d, _)| (Some(o), d));
            let entry = &mut ring_table[idx];
            match entry.take() {
                None => {
                    *entry = Some(OpenRing {
                        atom_id: last_atom_idx.unwrap(),
                        order: order_opt,
                        dir: dir_opt,
                        open_pos: i,
                        open_aromatic: last_aromatic,
                    });
                }
                Some(open) => {
                    let b = last_atom_idx.unwrap();
                    if open.atom_id == b {
                        return Err(M6Error::RingSelfLoop { pos: i });
                    }
                    if prev_atom_idx == Some(open.atom_id) {
                        return Err(M6Error::RingTwoMember { pos: i });
                    }
                    if let (Some(d1), Some(d2)) = (open.dir, dir_opt) {
                        if d1 != d2 {
                            return Err(M6Error::RingBondDirConflict {
                                pos: i,
                                open_pos: open.open_pos,
                            });
                        }
                    }
                    if let (Some(o1), Some(o2)) = (open.order, order_opt) {
                        if o1 != o2 {
                            return Err(M6Error::RingBondOrderConflict {
                                pos: i,
                                open_pos: open.open_pos,
                            });
                        }
                    }
                    if open.dir.is_some() || dir_opt.is_some() {
                        let ord = open.order.or(order_opt).unwrap_or(BondOrder::Single);
                        if ord != BondOrder::Single {
                            return Err(M6Error::RingBondOrderConflict {
                                pos: i,
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
                    let b = last_atom_idx.unwrap();
                    if final_order == BondOrder::Single && open.open_aromatic && last_aromatic {
                        final_order = BondOrder::Aromatic;
                    }
                    if a != b {
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
            }
            i += 3;
            continue;
        }

        // Bonds
        if matches!(b0, b'-' | b'=' | b'#' | b'$' | b':' | b'/' | b'\\') {
            if pending_bond.is_some() {
                return Err(M6Error::ConsecutiveBond { pos: i });
            }
            if last_atom_idx.is_none() {
                return Err(M6Error::LeadingBond { pos: i });
            }
            let (order, dir) = map_bond(b0);
            pending_bond = Some((order, dir, i));
            i += 1;
            continue;
        }

        // Brackets
        if b0 == b'[' {
            let start = i + 1;
            let mut j = start;
            while j < n && input[j] != b']' {
                j += 1;
            }
            if j >= n {
                return Err(M6Error::UnbalancedOpenBracket { pos: i });
            }
            let inner = &input[start..j];
            let (elem_opt, iso_opt, fields) = parse_bracket_m6(inner, i)?;
            let (element, aromatic) = match elem_opt {
                Some(e) => {
                    let first = inner.first().copied().unwrap_or_default();
                    (e, first.is_ascii_lowercase())
                }
                None => (Element::C, false),
            };
            let mut atom = AtomData {
                element,
                isotope: iso_opt,
                charge: None,
                hydrogen_count: None,
                class: None,
                aromatic,
                implicit_h: false,
                chirality: None,
                unknown_symbol: elem_opt.is_none(),
            };
            for f in fields {
                match f {
                    BracketField::Chiral(ch) => atom.chirality = Some(ch),
                    BracketField::HydrogenCount(h) => {
                        atom.hydrogen_count = Some((h as u8).min(u8::MAX))
                    }
                    BracketField::Charge(q) => atom.charge = Some(q),
                    BracketField::Class(c) => atom.class = Some(c),
                }
            }
            let curr = builder.on_atom(atom);
            if let Some(last) = last_atom_idx {
                if let Some((order, dir, _)) = pending_bond.take() {
                    builder.on_bond(last, curr, BondData { order, dir });
                } else {
                    if last_aromatic && aromatic {
                        builder.on_bond(
                            last,
                            curr,
                            BondData {
                                order: BondOrder::Aromatic,
                                dir: None,
                            },
                        );
                    } else {
                        builder.on_bond_single_fast(last, curr);
                    }
                }
                prev_atom_idx = Some(last);
            }
            last_atom_idx = Some(curr);
            last_aromatic = aromatic;
            if let Some(top) = pstack.last_mut() {
                match top {
                    Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                        *had_atom = true
                    }
                }
            }
            i = j + 1;
            continue;
        }

        // Two-letter halogens first
        if b0 == b'C' {
            if i + 1 < n && input[i + 1] == b'l' {
                let curr = builder.on_atom_fast(Element::Cl, true, false);
                if let Some(last) = last_atom_idx {
                    if let Some((order, dir, _)) = pending_bond.take() {
                        builder.on_bond(last, curr, BondData { order, dir });
                    } else {
                        builder.on_bond_single_fast(last, curr);
                    }
                }
                last_atom_idx = Some(curr);
                last_aromatic = false;
                if let Some(top) = pstack.last_mut() {
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
            if let Some(last) = last_atom_idx {
                if let Some((order, dir, _)) = pending_bond.take() {
                    builder.on_bond(last, curr, BondData { order, dir });
                } else {
                    builder.on_bond_single_fast(last, curr);
                }
                prev_atom_idx = Some(last);
            }
            last_atom_idx = Some(curr);
            last_aromatic = false;
            if let Some(top) = pstack.last_mut() {
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
                if let Some(last) = last_atom_idx {
                    if let Some((order, dir, _)) = pending_bond.take() {
                        builder.on_bond(last, curr, BondData { order, dir });
                    } else {
                        builder.on_bond_single_fast(last, curr);
                    }
                }
                last_atom_idx = Some(curr);
                last_aromatic = false;
                if let Some(top) = pstack.last_mut() {
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
            if let Some(last) = last_atom_idx {
                if let Some((order, dir, _)) = pending_bond.take() {
                    builder.on_bond(last, curr, BondData { order, dir });
                } else {
                    builder.on_bond_single_fast(last, curr);
                }
                prev_atom_idx = Some(last);
            }
            last_atom_idx = Some(curr);
            last_aromatic = false;
            if let Some(top) = pstack.last_mut() {
                match top {
                    Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                        *had_atom = true
                    }
                }
            }
            i += 1;
            continue;
        }

        // Single-letter organics
        let elem = match b0 {
            b'N' => Some(Element::N),
            b'O' => Some(Element::O),
            b'P' => Some(Element::P),
            b'S' => Some(Element::S),
            b'F' => Some(Element::F),
            b'I' => Some(Element::I),
            _ => None,
        };
        if let Some(element) = elem {
            let curr = builder.on_atom_fast(element, true, false);
            if let Some(last) = last_atom_idx {
                if let Some((order, dir, _)) = pending_bond.take() {
                    builder.on_bond(last, curr, BondData { order, dir });
                } else {
                    builder.on_bond_single_fast(last, curr);
                }
                prev_atom_idx = Some(last);
            }
            last_atom_idx = Some(curr);
            last_aromatic = false;
            if let Some(top) = pstack.last_mut() {
                match top {
                    Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                        *had_atom = true
                    }
                }
            }
            i += 1;
            continue;
        }

        // Aromatic bare atoms
        if matches!(b0, b'b' | b'c' | b'n' | b'o' | b'p' | b's') {
            let element = match b0 {
                b'b' => Element::B,
                b'c' => Element::C,
                b'n' => Element::N,
                b'o' => Element::O,
                b'p' => Element::P,
                _ => Element::S,
            };
            let curr = builder.on_atom_fast(element, true, true);
            if let Some(last) = last_atom_idx {
                if let Some((order, dir, _)) = pending_bond.take() {
                    builder.on_bond(last, curr, BondData { order, dir });
                } else {
                    if last_aromatic {
                        builder.on_bond(
                            last,
                            curr,
                            BondData {
                                order: BondOrder::Aromatic,
                                dir: None,
                            },
                        );
                    } else {
                        builder.on_bond_single_fast(last, curr);
                    }
                }
                prev_atom_idx = Some(last);
            }
            last_atom_idx = Some(curr);
            last_aromatic = true;
            if let Some(top) = pstack.last_mut() {
                match top {
                    Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                        *had_atom = true
                    }
                }
            }
            i += 1;
            continue;
        }

        // Bare wildcard '*'
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
            if let Some(last) = last_atom_idx {
                if let Some((order, dir, _)) = pending_bond.take() {
                    builder.on_bond(last, curr, BondData { order, dir });
                } else {
                    builder.on_bond_single_fast(last, curr);
                }
                prev_atom_idx = Some(last);
            }
            last_atom_idx = Some(curr);
            last_aromatic = false;
            if let Some(top) = pstack.last_mut() {
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
            return Err(M6Error::UnbalancedCloseBracket { pos: i });
        }
        return Err(M6Error::UnsupportedToken { pos: i });
    }

    if pending_bond.is_some() {
        let (_, _, pos) = pending_bond.unwrap();
        return Err(M6Error::TrailingBond { pos });
    }
    if !pstack.is_empty() {
        let pos = match pstack.last().unwrap() {
            Frame::Branch { open_pos, .. } | Frame::Group { open_pos, .. } => *open_pos,
        };
        return Err(M6Error::UnbalancedBranchOpen { pos });
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
        return Err(M6Error::RingUnclosed { open_pos: pos_open });
    }
    let mut mols = builder.finish();
    Ok(mols.pop().unwrap_or_default())
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_data::Element;

    use super::*;
    use crate::io::ir::{AtomSymbol, BondDir, BondSymbol, Chirality};
    use crate::io::smiles::test_support::build_from_graph;

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(b"", Molecule::default())]
    #[case::chain_c_1(b"C", build_from_graph("C |"))]
    #[case::chain_c_5(b"CCCCC", build_from_graph("C C C C C | 0-1 1-2 2-3 3-4"))]
    #[case::aromatic_c_6(b"cccccc", build_from_graph("C* C* C* C* C* C* | 0-1: 1-2: 2-3: 3-4: 4-5:"))]
    #[case::chain_mixed_5(b"CClOBrN", build_from_graph("C Cl O Br N | 0-1 1-2 2-3 3-4"))]
    fn m6_chain(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m6(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty_group(b"()", Molecule::default())]
    #[case::group_c_1(b"(C)", build_from_graph("C |"))]
    #[case::group_c_1_aromatic(b"(c)", build_from_graph("C* |"))]
    #[case::group_c_4(b"(CCCC)", build_from_graph("C C C C | 0-1 1-2 2-3"))]
    #[case::group_nested(b"((CC))", build_from_graph("C C | 0-1"))]
    #[case::branch_c_111(b"C(C)(C)", build_from_graph("C C C | 0-1 0-2"))]
    #[case::branch_c_211(b"CC(C)C", build_from_graph("C C C C | 0-1 1-2 1-3"))]
    #[case::branch_c_222_aromatic(b"cc(cc)cc", build_from_graph("C* C* C* C* C* C* | 0-1: 1-2: 2-3: 1-4: 4-5:"))]
    #[case::trailing_branch(b"C(CC)", build_from_graph("C C C | 0-1 1-2"))]
    fn m6_tree(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m6(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rstest]
    #[case::unbalanced_closing_paren_1(b")C", M6Error::UnbalancedBranchClose { pos: 0 })]
    #[case::unbalanced_closing_paren_2(b"C)C", M6Error::UnbalancedBranchClose { pos: 1 })]
    #[case::unclosed_group(b"(C", M6Error::UnbalancedBranchOpen { pos: 0 })]
    #[case::unclosed_branch(b"C(C", M6Error::UnbalancedBranchOpen { pos: 1 })]
    #[case::empty_branch(b"C()", M6Error::EmptyBranch { pos: 2 })]
    #[case::empty_group_before_atom(b"()C", M6Error::EmptyGroup { pos: 1 })]
    #[case::two_top_level_groups(b"(C)(C)", M6Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::three_top_level_groups(b"(C)(C)(C)", M6Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::three_top_level_groups_aromatic(b"(c)(c)(c)", M6Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::two_top_level_groups_rings(b"(C1CC1)(C2CC2)", M6Error::TopLevelGroupTrailing { pos: 6 })]
    #[case::group_before_atom(b"(C)C", M6Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::group_before_atom_aromatic(b"(c)c", M6Error::TopLevelGroupTrailing { pos: 2 })]
    fn m6_tree_invalid(#[case] input: &[u8], #[case] expected: M6Error) {
        let err = parse_smiles_m6(input);
        assert!(err.is_err(), "{:?} should have failed", input);
        let err = err.unwrap_err();
        assert_eq!(err, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ring_c_3(b"C1CC1", build_from_graph("C C C | 0-1 1-2 0-2"))]
    #[case::ring_c_10(b"C1CCCCCCCCC1", build_from_graph("C C C C C C C C C C | 0-1 1-2 2-3 3-4 4-5 5-6 6-7 7-8 8-9 0-9"))]
    #[case::ring_aromatic_c_6(b"c1ccccc1", build_from_graph("C* C* C* C* C* C* | 0-1: 1-2: 2-3: 3-4: 4-5: 0-5:"))]
    #[case::ring_index_0(b"C0CC0", build_from_graph("C C C | 0-1 1-2 0-2"))]
    #[case::ring_index_percent(b"C%12CC%12", build_from_graph("C C C | 0-1 1-2 0-2"))]
    #[case::ring_index_numeric_equiv_1(b"C1CC%01", build_from_graph("C C C | 0-1 1-2 0-2"))]
    #[case::ring_index_numeric_equiv_0(b"C0CC%00", build_from_graph("C C C | 0-1 1-2 0-2"))]
    #[case::ring_index_numeric_equiv_9(b"C9CC%09", build_from_graph("C C C | 0-1 1-2 0-2"))]
    #[case::ring_index_max_percent(b"C%99CC%99", build_from_graph("C C C | 0-1 1-2 0-2"))]
    #[case::two_rings_bonded_0(b"C1CC1C2CC2", build_from_graph("C C C C C C | 0-1 1-2 0-2 2-3 3-4 4-5 3-5"))]
    #[case::two_rings_bonded_0_aromatic(b"c1cc1c2cc2", build_from_graph("C* C* C* C* C* C* | 0-1: 1-2: 0-2: 2-3: 3-4: 4-5: 3-5:"))]
    #[case::two_rings_bonded_0_aromatic_1(b"c1cc1C2CC2", build_from_graph("C* C* C* C C C | 0-1: 1-2: 0-2: 2-3 3-4 4-5 3-5"))]
    #[case::two_rings_index_reused(b"C1CC1C1CC1", build_from_graph("C C C C C C | 0-1 1-2 0-2 2-3 3-4 4-5 3-5"))]
    #[case::two_rings_bonded_2(b"C1CC1CCC2CC2", build_from_graph("C C C C C C C C | 0-1 1-2 0-2 2-3 3-4 4-5 5-6 6-7 5-7"))]
    #[case::two_rings_spiro(b"C1CC12CC2", build_from_graph("C C C C C | 0-1 1-2 0-2 2-3 3-4 2-4"))]
    #[case::two_rings_fused(b"C12CC1C2", build_from_graph("C C C C | 0-1 1-2 0-2 2-3 0-3"))]
    #[case::two_rings_bridged(b"C12CC(C2)C1", build_from_graph("C C C C C | 0-1 1-2 2-3 0-3 2-4 0-4"))]
    #[case::two_rings_fused_aromatic(b"c12ccccc1cccc2", build_from_graph("C* C* C* C* C* C* C* C* C* C* | 0-1: 1-2: 2-3: 3-4: 4-5: 0-5: 5-6: 6-7: 7-8: 8-9: 0-9:"))]
    #[case::three_rings_fused(b"C12C3C1C32", build_from_graph("C C C C | 0-1 1-2 0-2 2-3 1-3 0-3"))]
    #[case::ring_group(b"(C1CC1)", build_from_graph("C C C | 0-1 1-2 0-2"))]
    #[case::ring_branch_1(b"CC(C1)(C1)", build_from_graph("C C C C | 0-1 1-2 1-3 2-3"))]
    #[case::ring_branch_2(b"C(C1)CC1", build_from_graph("C C C C | 0-1 0-2 2-3 1-3 "))]
    #[case::substituted_ring_1(b"CC1CC1", build_from_graph("C C C C | 0-1 1-2 2-3 1-3"))]
    #[case::substituted_ring_2(b"C1(C)CC1", build_from_graph("C C C C | 0-1 0-2 2-3 0-3"))]
    #[case::substituted_ring_3(b"C(C)1CC1", build_from_graph("C C C C | 0-1 0-2 2-3 0-3"))]
    #[case::substituted_ring_4(b"C1C(C)C1", build_from_graph("C C C C | 0-1 1-2 1-3 0-3"))]
    #[case::substituted_ring_5(b"C1CC(C)1", build_from_graph("C C C C | 0-1 1-2 2-3 0-2"))]
    #[case::substituted_ring_6(b"C1CC1C", build_from_graph("C C C C | 0-1 1-2 0-2 2-3"))]
    #[case::substituted_ring_7(b"C1CC1(C)", build_from_graph("C C C C | 0-1 1-2 0-2 2-3"))]
    #[case::substituted_ring_aromatic(b"c1c(c)c1", build_from_graph("C* C* C* C* | 0-1: 1-2: 1-3: 0-3:"))]
    #[case::substituted_ring_branch_1(b"C1C(C(C)C)C1", build_from_graph("C C C C C C | 0-1 1-2 2-3 2-4 1-5 0-5"))]
    fn m6_ring(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m6(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rstest]
    #[case::leading_ring(b"1C", M6Error::LeadingRing { pos: 0 })]
    #[case::bad_percent_short(b"C%1", M6Error::RingIndexInvalid { pos: 1 })]
    #[case::bad_percent_char(b"C%1a", M6Error::RingIndexInvalid { pos: 1 })]
    #[case::bad_percent_eoi(b"C%", M6Error::RingIndexInvalid { pos: 1 })]
    #[case::bad_percent_zero(b"C%0", M6Error::RingIndexInvalid { pos: 1 })]
    #[case::ring_self_loop(b"C11", M6Error::RingSelfLoop { pos: 2 })]
    #[case::ring_two_member(b"C1C1", M6Error::RingTwoMember { pos: 3 })]
    #[case::ring_bond_order_conflict_3(b"C=1CC#1", M6Error::RingBondOrderConflict { pos: 6, open_pos: 2 })]
    #[case::ring_bond_order_conflict_4(b"C/1CC=1", M6Error::RingBondOrderConflict { pos: 6, open_pos: 2 })]
    #[case::ring_bond_order_conflict_5(b"C\\1CC=1", M6Error::RingBondOrderConflict { pos: 6, open_pos: 2 })]
    #[case::ring_bond_order_conflict_6(b"C=1CC/1", M6Error::RingBondOrderConflict { pos: 6, open_pos: 2 })]
    #[case::ring_bond_order_conflict_7(b"C=1CC\\1", M6Error::RingBondOrderConflict { pos: 6, open_pos: 2 })]
    #[case::ring_bond_order_conflict_8(b"C=%10CC#%10", M6Error::RingBondOrderConflict { pos: 8, open_pos: 2 })]
    #[case::ring_bond_dir_conflict_1(b"C/1CC\\1", M6Error::RingBondDirConflict { pos: 6, open_pos: 2 })]
    #[case::ring_bond_dir_conflict_2(b"C\\1CC/1", M6Error::RingBondDirConflict { pos: 6, open_pos: 2 })]
    #[case::ring_bond_dir_conflict_3(b"C/%12CC\\%12", M6Error::RingBondDirConflict { pos: 8, open_pos: 2 })]
    #[case::ring_bond_dir_conflict_4(b"C\\%12CC/%12", M6Error::RingBondDirConflict { pos: 8, open_pos: 2 })]
    #[case::ring_unclosed_1(b"C1CC", M6Error::RingUnclosed { open_pos: 1 })]
    #[case::ring_unclosed_2(b"C1CC1C1", M6Error::RingUnclosed { open_pos: 6 })]
    fn m6_ring_invalid(#[case] input: &[u8], #[case] expected: M6Error) {
        let err = parse_smiles_m6(input);
        assert!(err.is_err(), "{:?} should have failed", input);
        let err = err.unwrap_err();
        assert_eq!(err, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single_bond(b"C-C", build_from_graph("C C | 0-1:-"))]
    #[case::double_bond(b"C=C", build_from_graph("C C | 0-1:="))]
    #[case::triple_bond(b"C#C", build_from_graph("C C | 0-1:#"))]
    #[case::quadruple_bond(b"C$C", build_from_graph("C C | 0-1:$"))]
    #[case::aromatic_bond(b"C:C", build_from_graph("C C | 0-1::"))]
    #[case::up_bond(b"C/C", build_from_graph("C C | 0-1:/"))]
    #[case::down_bond(b"C\\C", build_from_graph("C C | 0-1:\\"))]
    #[case::single_bond_aromatic(b"c-c", build_from_graph("C* C* | 0-1:-"))]
    #[case::double_bond_aromatic(b"c=c", build_from_graph("C* C* | 0-1:="))]
    #[case::triple_bond_aromatic(b"c#c", build_from_graph("C* C* | 0-1:#"))]
    #[case::quadruple_bond_aromatic(b"c$c", build_from_graph("C* C* | 0-1:$"))]
    #[case::aromatic_bond_aromatic(b"c:c", build_from_graph("C* C* | 0-1::"))]
    #[case::up_bond_aromatic(b"c/c", build_from_graph("C* C* | 0-1:/"))]
    #[case::down_bond_aromatic(b"c\\c", build_from_graph("C* C* | 0-1:\\"))]
    #[case::cumulated_bonds(b"C=C=C", build_from_graph("C C C | 0-1:= 1-2:="))]
    #[case::conjugated_bonds(b"C=CC=C", build_from_graph("C C C C | 0-1:= 1-2:- 2-3:="))]
    #[case::cumulated_bonds_aromatic(b"c=c=c", build_from_graph("C* C* C* | 0-1:= 1-2:="))]
    #[case::conjugated_bonds_aromatic(b"c=c-c=c", build_from_graph("C* C* C* C* | 0-1:= 1-2:- 2-3:="))]
    #[case::branch_leading_bond(b"CC(-C)C", build_from_graph("C C C C | 0-1 1-2 1-3"))]
    #[case::branch_leading_double_bond(b"CC(=C)C", build_from_graph("C C C C | 0-1 1-2:= 1-3"))]
    #[case::branch_internal_bond(b"CC(C-C)C", build_from_graph("C C C C C | 0-1 1-2 2-3 1-4"))]
    #[case::branch_internal_double_bond(b"CC(C=C)C", build_from_graph("C C C C C | 0-1 1-2 2-3:= 1-4"))]
    #[case::branch_followed_by_bond(b"CC(C)-C", build_from_graph("C C C C | 0-1 1-2 1-3"))]
    #[case::branch_followed_by_double_bond(b"CC(C)=C", build_from_graph("C C C C | 0-1 1-2 1-3:="))]
    #[case::branch_leading_bond_aromatic(b"cc(:c)c", build_from_graph("C* C* C* C* | 0-1: 1-2: 1-3:"))]
    #[case::branch_internal_bond_aromatic(b"cc(c:c)c", build_from_graph("C* C* C* C* C* | 0-1: 1-2: 2-3: 1-4:"))]
    #[case::branch_followed_by_bond_aromatic(b"cc(c):c", build_from_graph("C* C* C* C* | 0-1: 1-2: 1-3:"))]
    #[case::branch_trans_double_bond_1(b"C/C=C/C", build_from_graph("C C C C | 0-1:/ 1-2:= 2-3:/"))]
    #[case::branch_trans_double_bond_2(b"C\\C=C\\C", build_from_graph("C C C C | 0-1:\\ 1-2:= 2-3:\\"))]
    #[case::branch_cis_double_bond_1(b"C\\C=C/C", build_from_graph("C C C C | 0-1:\\ 1-2:= 2-3:/"))]
    #[case::branch_cis_double_bond_2(b"C/C=C\\C", build_from_graph("C C C C | 0-1:/ 1-2:= 2-3:\\"))]
    #[case::ring_single_bond(b"C-1-C-C-1", build_from_graph("C C C | 0-1 1-2 0-2"))]
    #[case::ring_double_bond_1(b"C1-C=C1", build_from_graph("C C C | 0-1 1-2:= 0-2"))]
    #[case::ring_double_bond_2(b"C1-CC=1", build_from_graph("C C C | 0-1 1-2 0-2:="))]
    #[case::ring_double_bond_3(b"C=1-CC1", build_from_graph("C C C | 0-1 1-2 0-2:="))]
    #[case::ring_double_bond_4(b"C=1-C-C=1", build_from_graph("C C C | 0-1 1-2 0-2:="))]
    #[case::ring_double_bond_unilateral_close(b"C1CC=1", build_from_graph("C C C | 0-1 1-2 0-2:="))]
    #[case::ring_double_bond_unilateral_open(b"C=1CC1", build_from_graph("C C C | 0-1 1-2 0-2:="))]
    #[case::ring_triple_bond(b"C1-C-C#1", build_from_graph("C C C | 0-1 1-2 0-2:#"))]
    #[case::ring_quadruple_bond(b"C1-C-C$1", build_from_graph("C C C | 0-1 1-2 0-2:$"))]
    #[case::ring_aromatic_bond(b"c1:c:c:1", build_from_graph("C* C* C* | 0-1: 1-2: 0-2:"))]
    #[case::ring_up_bond_1(b"C1CC/1", build_from_graph("C C C | 0-1 1-2 0-2:/"))]
    #[case::ring_up_bond_2(b"C/1CC1", build_from_graph("C C C | 0-1 1-2 0-2:/"))]
    #[case::ring_up_bond_3(b"C/1CC/1", build_from_graph("C C C | 0-1 1-2 0-2:/"))]
    #[case::ring_down_bond(b"C1CC\\1", build_from_graph("C C C | 0-1 1-2 0-2:\\"))]
    #[case::ring_down_bond_both(b"C\\1CC\\1", build_from_graph("C C C | 0-1 1-2 0-2:\\"))]
    #[case::ring_up_bond_percent_open(b"C/%12CC%12", build_from_graph("C C C | 0-1 1-2 0-2:/"))]
    #[case::ring_up_bond_percent_close(b"C%12CC/%12", build_from_graph("C C C | 0-1 1-2 0-2:/"))]
    #[case::ring_down_bond_percent_both(b"C\\%12CC\\%12", build_from_graph("C C C | 0-1 1-2 0-2:\\"))]
    #[case::ring_between_bonds(b"C1CC-1-C", build_from_graph("C C C C | 0-1 1-2 0-2 2-3"))]
    fn m6_bonds(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m6(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rstest]
    #[case::trailing_bond_1(b"C-", M6Error::TrailingBond { pos: 1 })]
    #[case::trailing_bond_2(b"C=", M6Error::TrailingBond { pos: 1 })]
    #[case::trailing_bond_3(b"C#", M6Error::TrailingBond { pos: 1 })]
    #[case::trailing_bond_4(b"C$", M6Error::TrailingBond { pos: 1 })]
    #[case::trailing_stereo_bond_1(b"C/", M6Error::TrailingBond { pos: 1 })]
    #[case::trailing_stereo_bond_2(b"C\\", M6Error::TrailingBond { pos: 1 })]
    #[case::trailing_aromatic_bond(b"C:", M6Error::TrailingBond { pos: 1 })]
    #[case::branch_trailing_bond_1(b"C(C-)C", M6Error::TrailingBond { pos: 3 })]
    #[case::branch_trailing_bond_2(b"C(C=)C", M6Error::TrailingBond { pos: 3 })]
    #[case::branch_trailing_stereo_bond(b"CC(C/)CC", M6Error::TrailingBond { pos: 4 })]
    #[case::group_trailing_bond_1(b"(C-)", M6Error::TrailingBond { pos: 2 })]
    #[case::group_trailing_bond_2(b"(C=)", M6Error::TrailingBond { pos: 2 })]
    #[case::group_trailing_stereo_bond(b"(C/)", M6Error::TrailingBond { pos: 2 })]
    #[case::group_trailing_aromatic_bond(b"(C:)", M6Error::TrailingBond { pos: 2 })]
    #[case::trailing_bond_before_dot_1(b"C-.C", M6Error::TrailingBond { pos: 1 })]
    #[case::trailing_bond_before_dot_2(b"C=.C", M6Error::TrailingBond { pos: 1 })]
    #[case::trailing_bond_before_dot_aromatic(b"C:.C", M6Error::TrailingBond { pos: 1 })]
    #[case::trailing_stereo_bond_before_dot_up(b"C/.C", M6Error::TrailingBond { pos: 1 })]
    #[case::trailing_stereo_bond_before_dot_down(b"C\\.C", M6Error::TrailingBond { pos: 1 })]
    #[case::bond_after_group_1(b"(C)-", M6Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::bond_after_group_2(b"(C)=", M6Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::group_after_group_1(b"(C)(C)", M6Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::group_after_group_2(b"(c)(c)", M6Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::ring_after_group(b"(C1CCC)1", M6Error::TopLevelGroupTrailing { pos : 6})]
    #[case::consecutive_bonds_1(b"C--C", M6Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bonds_2(b"C-=C", M6Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bonds_3(b"C-#C", M6Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bonds_4(b"C-$C", M6Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bonds_5(b"C-:C", M6Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_stereo_bonds_1(b"C//C", M6Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_stereo_bonds_2(b"C\\\\C", M6Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bond_and_stereo_bond_1(b"C-/C", M6Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bond_and_stereo_bond_2(b"C=\\C", M6Error::ConsecutiveBond { pos: 2 })]
    #[case::leading_bond_1(b"-C", M6Error::LeadingBond { pos: 0 })]
    #[case::leading_bond_2(b"=C", M6Error::LeadingBond { pos: 0 })]
    #[case::leading_bond_3(b"#C", M6Error::LeadingBond { pos: 0 })]
    #[case::leading_bond_4(b"$C", M6Error::LeadingBond { pos: 0 })]
    #[case::leading_aromatic_bond(b":C", M6Error::LeadingBond { pos: 0 })]
    #[case::leading_sterebond_1(b"/C", M6Error::LeadingBond { pos: 0 })]
    #[case::leading_sterebond_2(b"\\C", M6Error::LeadingBond { pos: 0 })]
    #[case::group_leading_bond_1(b"(-C)C", M6Error::LeadingBond { pos: 1 })]
    #[case::group_leading_bond_2(b"(=C)C", M6Error::LeadingBond { pos: 1 })]
    #[case::group_leading_bond_3(b"(#C)C", M6Error::LeadingBond { pos: 1 })]
    #[case::group_leading_bond_4(b"($C)C", M6Error::LeadingBond { pos: 1 })]
    #[case::group_leading_sterebond_1(b"(/C)C", M6Error::LeadingBond { pos: 1 })]
    #[case::group_leading_sterebond_2(b"(\\C)C", M6Error::LeadingBond { pos: 1 })]
    #[case::group_leading_aromatic_bond(b"(:C)C", M6Error::LeadingBond { pos: 1 })]
    fn m6_bonds_invalid(#[case] input: &[u8], #[case] expected: M6Error) {
        let err = parse_smiles_m6(input);
        assert!(err.is_err(), "{:?} should have failed", input);
        let err = err.unwrap_err();
        assert_eq!(err, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::components_2(b"CC.CC", build_from_graph("C C C C | 0-1 2-3"))]
    #[case::components_5(b"C.C.C.C.C", build_from_graph("C C C C C | "))]
    #[case::ring_components_1(b"C1.CC1", build_from_graph("C C C | 1-2 0-2"))]
    #[case::ring_components_2(b"C%12.CC%12", build_from_graph("C C C | 1-2 0-2"))]
    #[case::ring_components_aromatic(b"c1.ccccc1", build_from_graph("C* C* C* C* C* C* | 1-2: 2-3: 3-4: 4-5: 0-5:"))]
    #[case::branch_components(b"C(C.C)", build_from_graph("C C C | 0-1"))]
    #[case::leading_dot_in_branch_1(b"C(.C)", build_from_graph("C C | "))]
    #[case::leading_dot_in_branch_2(b"C(.C)(C)", build_from_graph("C C C | 0-2"))]
    #[case::leading_dot_in_branch_3(b"C(C)(.C)", build_from_graph("C C C | 0-1"))]
    #[case::trailing_dot_in_branch_1(b"C(C.)", build_from_graph("C C | 0-1"))]
    #[case::trailing_dot_in_branch_2(b"C(C.)C", build_from_graph("C C C | 0-1 0-2"))]
    #[case::trailing_dot_in_branch_3(b"C(C.)(C)", build_from_graph("C C C | 0-1 0-2"))]
    #[case::group_components_1(b"(C.CC.C)", build_from_graph("C C C C | 1-2"))]
    #[case::group_components_2(b"(CC).(CC)", build_from_graph("C C C C | 0-1 2-3"))]
    #[case::group_components_3(b"(C.C).C", build_from_graph("C C C |"))]
    #[case::group_components_4(b"C.(C).C", build_from_graph("C C C |"))]
    #[case::group_components_5(b"C.C.(C)", build_from_graph("C C C |"))]
    #[case::leading_dot_in_group_1(b"(.CC)", build_from_graph("C C | 0-1"))]
    #[case::leading_dot_in_group_2(b"(.CC).(CC)", build_from_graph("C C C C | 0-1 2-3"))]
    #[case::leading_dot_in_group_3(b"(CC).(.CC)", build_from_graph("C C C C | 0-1 2-3"))]
    #[case::leading_dot_in_group_4(b"C(.C.C)", build_from_graph("C C C |"))]
    #[case::trailing_dot_in_group_1(b"(CC.)", build_from_graph("C C | 0-1"))]
    #[case::trailing_dot_in_group_2(b"(CC.).CC", build_from_graph("C C C C | 0-1 2-3"))]
    #[case::trailing_dot_in_group_3(b"(CC).(CC.)", build_from_graph("C C C C | 0-1 2-3"))]
    #[case::group_ring_components_1(b"(CC1.C1)", build_from_graph("C C C | 0-1 1-2"))]
    #[case::group_ring_components_2(b"C1.(C).CC1", build_from_graph("C C C C | 2-3 0-3 "))]
    #[case::group_ring_components_3(b"C%12.(C).CC%12", build_from_graph("C C C C | 2-3 0-3 "))]
    #[case::rings_across_multiple_dots_digit(b"C1.C.CC1", build_from_graph("C C C C | 2-3 0-3"))]
    #[case::rings_across_multiple_dots_percent(b"C%12.C.CC%12", build_from_graph("C C C C | 2-3 0-3"))]
    #[case::ring_double_unilateral_open(b"C=1.CC1", build_from_graph("C C C | 1-2 0-2:="))]
    #[case::ring_double_unilateral_close(b"C1.CC=1", build_from_graph("C C C | 1-2 0-2:="))]
    #[case::ring_dir_up_both(b"C/1.CC/1", build_from_graph("C C C | 1-2 0-2:/"))]
    #[case::ring_dir_down_both(b"C\\1.CC\\1", build_from_graph("C C C | 1-2 0-2:\\"))]
    #[case::ring_dir_up_both_percent(b"C/%12.CC/%12", build_from_graph("C C C | 1-2 0-2:/"))]
    #[case::ring_dir_down_both_percent(b"C\\%12.CC\\%12", build_from_graph("C C C | 1-2 0-2:\\"))]
    #[case::branch_multiple_components(b"C(.C.C)", build_from_graph("C C C |"))]
    #[case::groups_leading_dot_both(b"(.C).(.C)", build_from_graph("C C |"))]
    #[case::group_leading_dot_middle(b"C.(.C).C", build_from_graph("C C C |"))]
    fn m6_components(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m6(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rstest]
    #[case::leading_dot_1(b".", M6Error::LeadingDot { pos: 0 })]
    #[case::leading_dot_2(b".C", M6Error::LeadingDot { pos: 0 })]
    #[case::leading_dot_3(b"..C", M6Error::LeadingDot { pos: 0 })]
    #[case::leading_dot_4(b".C.", M6Error::LeadingDot { pos: 0 })]
    #[case::trailing_dot_1(b"C.", M6Error::TrailingDot { pos: 1 })]
    #[case::trailing_dot_2(b"C..", M6Error::ConsecutiveDot { pos: 1 })]
    #[case::double_dot(b"C..C", M6Error::ConsecutiveDot { pos: 1 })]
    #[case::dot_before_ring_digit(b"C.1", M6Error::LeadingRing { pos: 2 })]
    #[case::dot_before_ring_percent(b"C.%12", M6Error::LeadingRing { pos: 2 })]
    #[case::dot_in_group_1(b"(.)", M6Error::LeadingDot { pos: 2 })]
    #[case::dot_in_group_2(b"(.)C", M6Error::EmptyGroup { pos: 2 })]
    #[case::dot_in_group_3(b"(.).C", M6Error::EmptyGroup { pos: 2 })]
    #[case::dot_before_group(b"C.(C)C", M6Error::TopLevelGroupTrailing { pos: 4 })]
    #[case::dot_in_branch_1(b"C(.)", M6Error::EmptyBranch { pos: 3 })]
    #[case::dot_in_branch_2(b"C(.)C", M6Error::EmptyBranch { pos: 3 })]
    #[case::dot_in_branch_3(b"C(.)(C)", M6Error::EmptyBranch { pos: 3 })]
    #[case::dot_in_component_1(b"().C", M6Error::EmptyGroup { pos: 1 })]
    #[case::dot_in_component_2(b"(.).C", M6Error::EmptyGroup { pos: 2 })]
    #[case::dot_in_component_3(b"(.).(C)", M6Error::EmptyGroup { pos: 2 })]
    #[case::dot_in_component_4(b"C.()", M6Error::EmptyGroup { pos: 3 })]
    #[case::dot_in_component_5(b"C.(.)", M6Error::LeadingDot { pos: 4 })]
    #[case::dot_in_component_6(b"(C).(.)", M6Error::LeadingDot { pos: 6 })]
    #[case::dot_unclosed_ring_1(b"C1.C", M6Error::RingUnclosed { open_pos: 1 })]
    #[case::dot_unclosed_ring_2(b"C.C1", M6Error::RingUnclosed { open_pos: 3 })]
    #[case::dot_unclosed_ring_before_group(b"C1.(C)(C)C1", M6Error::TopLevelGroupTrailing { pos: 5 })]
    #[case::ring_order_conflict_digit(b"C=1.CC#1", M6Error::RingBondOrderConflict { pos: 7, open_pos: 2 })]
    #[case::ring_order_conflict_percent(b"C=%12.CC#%12", M6Error::RingBondOrderConflict { pos: 9, open_pos: 2 })]
    #[case::ring_dir_conflict_digit(b"C/1.CC\\1", M6Error::RingBondDirConflict { pos: 7, open_pos: 2 })]
    #[case::ring_dir_conflict_percent(b"C/%12.CC\\%12", M6Error::RingBondDirConflict { pos: 9, open_pos: 2 })]
    #[case::ring_dir_conflict_aromatic(b"c/1.cc\\1", M6Error::RingBondDirConflict { pos: 7, open_pos: 2 })]
    #[case::group_dot_before_ring_digit(b"(.1)", M6Error::LeadingRing { pos: 2 })]
    #[case::group_dot_before_ring_percent(b"(.%12)", M6Error::LeadingRing { pos: 2 })]
    #[case::branch_dot_before_ring_digit(b"C(.1)", M6Error::LeadingRing { pos: 3 })]
    #[case::branch_dot_before_ring_percent(b"C(.%12)", M6Error::LeadingRing { pos: 3 })]
    #[case::branch_dot_before_bond_1(b"(.-C)", M6Error::LeadingBond { pos: 2 })]
    #[case::branch_dot_before_bond_2(b"C(.-C)", M6Error::LeadingBond { pos: 3 })]
    #[case::leading_bond_after_dot_1(b"C.-C", M6Error::LeadingBond { pos: 2 })]
    #[case::leading_bond_after_dot_2(b"C.=-C", M6Error::LeadingBond { pos: 2 })]
    #[case::leading_sterebond_after_dot_up(b"C./C", M6Error::LeadingBond { pos: 2 })]
    #[case::leading_sterebond_after_dot_down(b"C.\\C", M6Error::LeadingBond { pos: 2 })]
    #[case::trailing_bond_dot_aromatic(b"C:.", M6Error::TrailingBond { pos: 1 })]
    #[case::group_trailing_bond_dot(b"(C-.)", M6Error::TrailingBond { pos: 2 })]
    #[case::branch_trailing_bond_dot(b"C(C-.)", M6Error::TrailingBond { pos: 3 })]
    fn m6_components_invalid(#[case] input: &[u8], #[case] expected: M6Error) {
        let err = parse_smiles_m6(input);
        assert!(err.is_err(), "{:?} should have failed", input);
        let err = err.unwrap_err();
        assert_eq!(err, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::basic_c(b"[C]", Some(Element::C), false, None, None, None, None, None)]
    #[case::basic_aromatic_c(b"[c]", Some(Element::C), true, None, None, None, None, None)]
    #[case::wildcard(b"[*]", None, false, None, None, None, None, None)]
    #[case::isotope_elem(b"[13C]", Some(Element::C), false, Some(13), None, None, None, None)]
    #[case::isotope_zero(b"[0C]", Some(Element::C), false, Some(0), None, None, None, None)]
    #[case::isotope_wild(b"[13*]", None, false, Some(13), None, None, None, None)]
    #[case::chirality_cw(b"[C@]", Some(Element::C), false, None, Some(Chirality::Clockwise), None, None, None)]
    #[case::chirality_ccw(b"[C@@]", Some(Element::C), false, None, Some(Chirality::CounterClockwise), None, None, None)]
    #[case::chirality_th2(b"[C@TH2]", Some(Element::C), false, None, Some(Chirality::Tetrahedral { arr: 2 }), None, None, None)]
    #[case::chirality_al1(b"[C@AL1]", Some(Element::C), false, None, Some(Chirality::Allenal { arr: 1 }), None, None, None)]
    #[case::chirality_sp3(b"[C@SP3]", Some(Element::C), false, None, Some(Chirality::SquarePlanar { arr: 3 }), None, None, None)]
    #[case::chirality_tb5(b"[C@TB5]", Some(Element::C), false, None, Some(Chirality::TrigonalBipyramidal { arr: 5 }), None, None, None)]
    #[case::chirality_oh7(b"[C@OH7]", Some(Element::C), false, None, Some(Chirality::Octahedral { arr: 7 }), None, None, None)]
    #[case::hcount(b"[CH]", Some(Element::C), false, None, None, Some(1), None, None)]
    #[case::hcount_1(b"[CH1]", Some(Element::C), false, None, None, Some(1), None, None)]
    #[case::hcount_0(b"[CH0]", Some(Element::C), false, None,None, Some(0), None, None)]
    #[case::hcount_3(b"[CH3]", Some(Element::C), false, None, None, Some(3), None, None)]
    #[case::hcount_aromatic(b"[cH]", Some(Element::C), true, None, None, Some(1), None, None)]
    #[case::wildcard_h1(b"[*H]", None, false, None, None, Some(1), None, None)]
    #[case::wildcard_h2(b"[*H2]", None, false, None, None, Some(2), None, None)]
    #[case::wildcard_h0(b"[*H0]", None, false, None, None, Some(0), None, None)]
    #[case::chirality_cw_hydrogen(b"[C@H]", Some(Element::C), false, None, Some(Chirality::Clockwise), Some(1), None, None)]
    #[case::chirality_ccw_hydrogen(b"[C@@H]", Some(Element::C), false, None, Some(Chirality::CounterClockwise), Some(1), None, None)]
    #[case::charge_plus(b"[C+]", Some(Element::C), false, None, None, None, Some(1), None)]
    #[case::charge_minus(b"[C-]", Some(Element::C), false, None, None, None, Some(-1), None)]
    #[case::charge_pp(b"[C++]", Some(Element::C), false, None, None, None, Some(2), None)]
    #[case::charge_mm(b"[C--]", Some(Element::C), false, None, None, None, Some(-2), None)]
    #[case::zero_charge_pos(b"[C+0]", Some(Element::C), false, None, None, None, Some(0), None)]
    #[case::zero_charge_neg(b"[C-0]", Some(Element::C), false, None, None, None, Some(0), None)]
    #[case::charge_plus_10(b"[C+10]", Some(Element::C), false, None, None, None, Some(10), None)]
    #[case::charge_minus_10(b"[C-10]", Some(Element::C), false, None, None, None, Some(-10), None)]
    #[case::charge_plus_hcount(b"[C+H]", Some(Element::C), false, None, None, Some(1), Some(1), None)]
    #[case::charge_plus_1_hcount(b"[C+1H]", Some(Element::C), false, None, None, Some(1), Some(1), None)]
    #[case::charge_minus_hcount(b"[C-H]", Some(Element::C), false, None, None, Some(1), Some(-1), None)]
    #[case::charge_minus_1_hcount(b"[C-1H]", Some(Element::C), false, None, None, Some(1), Some(-1), None)]
    #[case::class_elem(b"[C:12]", Some(Element::C), false, None, None, None, None, Some(12))]
    #[case::class_wild(b"[*:5]", None, false, None, None, None, None, Some(5))]
    fn m6_bracket(
        #[case] input: &[u8],
        #[case] elem: Option<Element>,
        #[case] aromatic: bool,
        #[case] isotope: Option<u32>,
        #[case] chirality: Option<Chirality>,
        #[case] hcount: Option<u32>,
        #[case] charge: Option<i32>,
        #[case] class_: Option<u32>,
    ) {
        let res = parse_smiles_m6(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol.atoms.len(), 1, "expected single atom");
        let a = &mol.atoms[0];
        match elem {
            Some(e) => match &a.symbol {
                AtomSymbol::Element(el) => assert_eq!(*el, e),
                other => panic!("expected element {:?}, got {:?}", e, other),
            },
            None => assert!(matches!(a.symbol, AtomSymbol::Unknown)),
        }
        assert_eq!(a.aromatic, Some(aromatic));
        assert_eq!(a.isotope, isotope);
        assert_eq!(a.chirality, chirality);
        assert_eq!(a.hydrogen_count, hcount);
        assert_eq!(a.charge, charge);
        assert_eq!(a.class, class_);
    }

    #[rstest]
    #[case::aliphatic_before(b"C[C]", BondOrder::Single, None)]
    #[case::aliphatic_before_single(b"C-[C]", BondOrder::Single, None)]
    #[case::aliphatic_before_double(b"C=[C]", BondOrder::Double, None)]
    #[case::aliphatic_before_triple(b"C#[C]", BondOrder::Triple, None)]
    #[case::aliphatic_before_quadruple(b"C$[C]", BondOrder::Quadruple, None)]
    #[case::aliphatic_before_aromatic(b"C:[C]", BondOrder::Aromatic, None)]
    #[case::aliphatic_before_up(b"C/[C]", BondOrder::Single, Some(BondDir::Up))]
    #[case::aliphatic_before_down(b"C\\[C]", BondOrder::Single, Some(BondDir::Down))]
    #[case::aliphatic_after(b"[C]C", BondOrder::Single, None)]
    #[case::aliphatic_after_single(b"[C]-C", BondOrder::Single, None)]
    #[case::aliphatic_after_double(b"[C]=C", BondOrder::Double, None)]
    #[case::aliphatic_after_triple(b"[C]#C", BondOrder::Triple, None)]
    #[case::aliphatic_after_quadruple(b"[C]$C", BondOrder::Quadruple, None)]
    #[case::aliphatic_after_aromatic(b"[C]:C", BondOrder::Aromatic, None)]
    #[case::aliphatic_after_up(b"[C]/C", BondOrder::Single, Some(BondDir::Up))]
    #[case::aliphatic_after_down(b"[C]\\C", BondOrder::Single, Some(BondDir::Down))]
    #[case::aromatic_before(b"c[c]", BondOrder::Aromatic, None)]
    #[case::aromatic_before_single(b"c-[c]", BondOrder::Single, None)]
    #[case::aromatic_before_aromatic(b"c:[c]", BondOrder::Aromatic, None)]
    #[case::aromatic_after(b"[c]c", BondOrder::Aromatic, None)]
    #[case::aromatic_after_single(b"[c]-c", BondOrder::Single, None)]
    #[case::aromatic_after_aromatic(b"[c]:c", BondOrder::Aromatic, None)]
    #[case::aliphatic_before_aromatic(b"C[c]", BondOrder::Single, None)]
    #[case::aliphatic_single_before_aromatic(b"C-[c]", BondOrder::Single, None)]
    #[case::aliphatic_aromatic_before_aromatic(b"C:[c]", BondOrder::Aromatic, None)]
    #[case::aliphatic_after_aromatic(b"[c]C", BondOrder::Single, None)]
    #[case::aromatic_after_aliphatic(b"[C]c", BondOrder::Single, None)]
    #[case::aromatic_after_aliphatic_single(b"[C]-c", BondOrder::Single, None)]
    #[case::aromatic_after_aliphatic_aromatic(b"[c]:c", BondOrder::Aromatic, None)]
    #[case::aromatic_after_aliphatic_up(b"[C]/c", BondOrder::Single, Some(BondDir::Up))]
    #[case::aromatic_after_aliphatic_down(b"[C]\\c", BondOrder::Single, Some(BondDir::Down))]
    #[case::bracket_branch_1(b"[C](C)", BondOrder::Single, None)]
    #[case::bracket_branch_2(b"C([C])", BondOrder::Single, None)]
    #[case::bracket_branch_single(b"C(-[C])", BondOrder::Single, None)]
    #[case::bracket_branch_double(b"C(=[C])", BondOrder::Double, None)]
    #[case::bracket_branch_triple(b"C(#[C])", BondOrder::Triple, None)]
    #[case::bracket_branch_quadruple(b"C($[C])", BondOrder::Quadruple, None)]
    #[case::bracket_branch_aromatic(b"C(:[C])", BondOrder::Aromatic, None)]
    #[case::bracket_branch_up(b"C(/[C])", BondOrder::Single, Some(BondDir::Up))]
    #[case::bracket_branch_down(b"C(\\[C])", BondOrder::Single, Some(BondDir::Down))]
    #[case::bracket_branch_down(b"C(\\[C])", BondOrder::Single, Some(BondDir::Down))]
    #[case::bracket_group_1(b"([C]C)", BondOrder::Single, None)]
    #[case::bracket_group_1(b"(C[C])", BondOrder::Single, None)]
    #[case::bracket_ring_1(b"[C]1CC1", BondOrder::Single, None)]
    #[case::bracket_ring_2(b"[C]1cc1", BondOrder::Single, None)]
    #[case::bracket_ring_double_1(b"[C]1=cc1", BondOrder::Double, None)]
    #[case::bracket_ring_double_2(b"[C]=1cc1", BondOrder::Single, None)]
    #[case::bracket_aromatic_ring(b"[c]1cc1", BondOrder::Aromatic, None)]
    fn m6_bracket_bonds(
        #[case] input: &[u8],
        #[case] expected: BondOrder,
        #[case] dir: Option<BondDir>,
    ) {
        let res = parse_smiles_m6(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol.bonds[0].symbol, BondSymbol::Bond(expected));
        assert_eq!(mol.bonds[0].direction, dir);
    }

    #[rstest]
    #[case::empty_bracket(b"[]", M6Error::InvalidBracket { pos: 0 })]
    #[case::bracket_in_chain_empty(b"C[]", M6Error::InvalidBracket { pos: 1 })]
    #[case::bracket_in_group_empty(b"(C[])", M6Error::InvalidBracket { pos: 2 })]
    #[case::bracket_in_branch_empty(b"C([])C", M6Error::InvalidBracket { pos: 2 })]
    #[case::bracket_in_component_empty(b"[].C", M6Error::InvalidBracket { pos: 0 })]
    #[case::invalid_element_1(b"[X]", M6Error::InvalidBracket { pos: 0 })]
    #[case::invalid_element_2(b"[Z]", M6Error::InvalidBracket { pos: 0 })]
    #[case::invalid_element_3(b"[Aq]", M6Error::InvalidBracket { pos: 0 })]
    #[case::invalid_element_4(b"[Sh]", M6Error::InvalidBracket { pos: 0 })]
    #[case::two_elements_1(b"[CF]", M6Error::InvalidBracket { pos: 0 })]
    #[case::two_elements_2(b"[AsF]", M6Error::InvalidBracket { pos: 0 })]
    #[case::two_elements_3(b"[FAs]", M6Error::InvalidBracket { pos: 0 })]
    #[case::two_elements_4(b"[AsBr]", M6Error::InvalidBracket { pos: 0 })]
    #[case::two_elements_wildcard_1(b"[*C]", M6Error::InvalidBracket { pos: 0 })]
    #[case::two_elements_wildcard_2(b"[C*]", M6Error::InvalidBracket { pos: 0 })]
    #[case::wildcard_invalid_element_1(b"[*X]", M6Error::InvalidBracket { pos: 0 })]
    #[case::wildcard_invalid_element_2(b"[X*]", M6Error::InvalidBracket { pos: 0 })]
    #[case::double_wildcard(b"[**]", M6Error::InvalidBracket { pos: 0 })]
    #[case::zero_charge_no_sign(b"[C0]", M6Error::InvalidBracket { pos: 0 })]
    #[case::pos_charge_no_sign(b"[C1]", M6Error::InvalidBracket { pos: 0 })]
    #[case::charge_no_element_1(b"[+]", M6Error::InvalidBracket { pos: 0 })]
    #[case::charge_no_element_2(b"[-]", M6Error::InvalidBracket { pos: 0 })]
    #[case::charge_no_element_3(b"[+0]", M6Error::InvalidBracket { pos: 0 })]
    #[case::charge_no_element_4(b"[-0]", M6Error::InvalidBracket { pos: 0 })]
    #[case::charge_no_element_5(b"[+1]", M6Error::InvalidBracket { pos: 0 })]
    #[case::charge_no_element_6(b"[-1]", M6Error::InvalidBracket { pos: 0 })]
    #[case::zero_isotope_no_element(b"[0]", M6Error::InvalidBracket { pos: 0 })]
    #[case::isotope_no_element(b"[13]", M6Error::InvalidBracket { pos: 0 })]
    #[case::chirality_no_element_1(b"[@]", M6Error::InvalidBracket { pos: 0 })]
    #[case::chirality_no_element_2(b"[@@]", M6Error::InvalidBracket { pos: 0 })]
    #[case::chirality_no_element_4(b"[@@TH1]", M6Error::InvalidBracket { pos: 0 })]
    #[case::class_no_element(b"[:12]", M6Error::InvalidBracket { pos: 0 })]
    #[case::hcount_two_digits(b"[CH10]", M6Error::BracketHCountTwoDigits { pos: 0 })]
    #[case::colon_no_class(b"[C:]", M6Error::BracketEmptyClass { pos: 0 })]
    #[case::unbalanced_open_bracket_1(b"[", M6Error::UnbalancedOpenBracket { pos: 0 })]
    #[case::unbalanced_open_bracket_2(b"C[", M6Error::UnbalancedOpenBracket { pos: 1 })]
    #[case::unbalanced_open_bracket_3(b"[C", M6Error::UnbalancedOpenBracket { pos: 0 })]
    #[case::unbalanced_open_bracket_4(b"[*", M6Error::UnbalancedOpenBracket { pos: 0 })]
    #[case::unbalanced_open_bracket_5(b"[)", M6Error::UnbalancedOpenBracket { pos: 0 })]
    #[case::unbalanced_open_bracket_6(b"[[", M6Error::UnbalancedOpenBracket { pos: 0 })]
    #[case::unbalanced_open_bracket_7(b"[.", M6Error::UnbalancedOpenBracket { pos: 0 })]
    #[case::unbalanced_open_bracket_8(b"C[", M6Error::UnbalancedOpenBracket { pos: 1 })]
    #[case::unbalanced_open_bracket_9(b"[C)", M6Error::UnbalancedOpenBracket { pos: 0 })]
    #[case::unbalanced_open_bracket_10(b"[.C", M6Error::UnbalancedOpenBracket { pos: 0 })]
    #[case::unbalanced_open_bracket_11(b"C.[", M6Error::UnbalancedOpenBracket { pos: 2 })]
    #[case::dot_in_bracket(b"[.]", M6Error::InvalidBracket { pos: 0 })]
    #[case::branch_open_in_bracket(b"[(]", M6Error::InvalidBracket { pos: 0 })]
    #[case::branch_close_in_bracket(b"[)]", M6Error::InvalidBracket { pos: 0 })]
    #[case::bracket_in_bracket_1(b"[[]", M6Error::InvalidBracket { pos: 0 })]
    #[case::bracket_in_bracket_2(b"[]]", M6Error::InvalidBracket { pos: 0 })]
    #[case::open_bracket_in_branch(b"C([)", M6Error::UnbalancedOpenBracket { pos: 2 })]
    #[case::close_bracket_in_branch(b"C(])", M6Error::UnbalancedCloseBracket { pos: 2 })]
    #[case::unbalanced_close_bracket_1(b"]", M6Error::UnbalancedCloseBracket { pos: 0 })]
    #[case::unbalanced_close_bracket_2(b"]C", M6Error::UnbalancedCloseBracket { pos: 0 })]
    #[case::unbalanced_close_bracket_3(b"C]", M6Error::UnbalancedCloseBracket { pos: 1 })]
    #[case::unbalanced_close_bracket_4(b"*]", M6Error::UnbalancedCloseBracket { pos: 1 })]
    #[case::unbalanced_close_bracket_5(b"C.]", M6Error::UnbalancedCloseBracket { pos: 2 })]
    #[case::unbalanced_close_bracket_6(b"].", M6Error::UnbalancedCloseBracket { pos: 0 })]
    #[case::unbalanced_close_bracket_7(b"].C", M6Error::UnbalancedCloseBracket { pos: 0 })]
    #[case::unbalanced_close_bracket_8(b"(]", M6Error::UnbalancedCloseBracket { pos: 1 })]
    #[case::unbalanced_close_bracket_9(b"(C]", M6Error::UnbalancedCloseBracket { pos: 2 })]
    fn m6_bracket_invalid(#[case] input: &[u8], #[case] expected: M6Error) {
        let err = parse_smiles_m6(input);
        assert!(err.is_err(), "{:?} should have failed", input);
        let err = err.unwrap_err();
        assert_eq!(err, expected);
    }

    #[rstest]
    #[case::wildcard(b"*", 0, false)]
    #[case::two_wildcards(b"**", 0, true)]
    #[case::wildcard_after_c(b"C*", 1, true)]
    #[case::wildcard_before_c(b"*C", 0, true)]
    #[case::wildcard_bond_single(b"C-*", 1, true)]
    #[case::wildcard_bond_single_rev(b"*-C", 0, true)]
    #[case::wildcard_branch_1(b"*(C)", 0, true)]
    #[case::wildcard_branch_2(b"C(*)", 1, true)]
    #[case::wildcard_branch_3(b"C(*C)", 1, true)]
    #[case::wildcard_group_1(b"(*)", 0, false)]
    #[case::wildcard_group_2(b"(*C)", 0, true)]
    #[case::wildcard_group_3(b"(C*)", 1, true)]
    #[case::wildcard_ring_1(b"*1CC1", 0, true)]
    #[case::wildcard_ring_2(b"C1*C1", 1, true)]
    #[case::wildcard_ring_3(b"C1C*1", 2, true)]
    #[case::wildcard_component_1(b"*.C", 0, false)]
    #[case::wildcard_component_2(b"C.*", 1, false)]
    #[case::wildcard_dot_bond_1(b"*1.C1", 0, true)]
    #[case::wildcard_dot_bond_2(b"C1.*1", 1, true)]
    fn m6_wildcard(#[case] input: &[u8], #[case] star_idx: usize, #[case] has_bonds: bool) {
        let res = parse_smiles_m6(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert!(star_idx < mol.atoms.len());
        let a = &mol.atoms[star_idx];
        assert!(matches!(a.symbol, AtomSymbol::Unknown));
        assert_eq!(a.isotope, Some(0));
        assert_eq!(a.charge, Some(0));
        assert_eq!(a.hydrogen_count, Some(0));
        assert_eq!(a.aromatic, Some(false));
        assert_eq!(a.implicit_h, false);
        if has_bonds {
            assert!(mol.bonds.len() > 0);
        }
    }

    #[rstest]
    #[case::wildcard_after_group(b"(C)*", M6Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::wildcard_unclosed_ring(b"*1", M6Error::RingUnclosed { open_pos: 1 })]
    #[case::wildcard_unclosed_branch(b"C(*", M6Error::UnbalancedBranchOpen { pos: 1 })]
    #[case::wildcard_unclosed_group(b"(C*", M6Error::UnbalancedBranchOpen { pos: 0 })]
    #[case::wildcard_unclosed_bracket(b"[*", M6Error::UnbalancedOpenBracket { pos: 0 })]
    #[case::wildcard_trailing_bond(b"*-", M6Error::TrailingBond { pos: 1 })]
    #[case::wildcard_trailing_dot(b"*.", M6Error::TrailingDot { pos: 1 })]
    fn m6_wildcard_invalid(#[case] input: &[u8], #[case] expected: M6Error) {
        let res = parse_smiles_m6(input);
        assert!(res.is_err(), "{:?} should have failed", input);
        let err = res.unwrap_err();
        assert_eq!(err, expected);
    }

    #[rstest]
    #[case::terminator_space(b"CC ", build_from_graph("C C | 0-1"))]
    #[case::terminator_tab(b"CC\t", build_from_graph("C C | 0-1"))]
    #[case::terminator_newline(b"CC\n", build_from_graph("C C | 0-1"))]
    #[case::terminator_cr(b"CC\r", build_from_graph("C C | 0-1"))]
    #[case::terminator_crlf(b"CC\r\n", build_from_graph("C C | 0-1"))]
    fn m6_whitespace_strict(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m6(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rstest]
    #[case::leading_space(b" CC", M6Error::InvalidWhitespace { pos: 0 })]
    #[case::leading_tab(b"\tCC", M6Error::InvalidWhitespace { pos: 0 })]
    #[case::leading_newline(b"\nCC", M6Error::InvalidWhitespace { pos: 0 })]
    #[case::leading_cr(b"\rCC", M6Error::InvalidWhitespace { pos: 0 })]
    #[case::leading_crlf(b"\r\nCC", M6Error::InvalidWhitespace { pos: 0 })]
    #[case::terminator_space_trailing_structure(b"CC CC", M6Error::InvalidWhitespace { pos: 2 })]
    #[case::terminator_tab_trailing_structure(b"CC\tCC", M6Error::InvalidWhitespace { pos: 2 })]
    #[case::terminator_cr_trailing_structure(b"CC\rCC", M6Error::InvalidWhitespace { pos: 2 })]
    #[case::terminator_newline_trailing_structure(b"CC\nCC", M6Error::InvalidWhitespace { pos: 2 })]
    #[case::terminator_crlf_trailing_structure(b"CC\r\nCC", M6Error::InvalidWhitespace { pos: 2 })]
    fn m6_whitespace_strict_invalid(#[case] input: &[u8], #[case] expected: M6Error) {
        let res = parse_smiles_m6(input);
        assert!(res.is_err(), "{:?} should have failed", input);
        let err = res.unwrap_err();
        assert_eq!(err, expected);
    }

    #[rstest]
    #[case::ws_intertoken_spaces_flags(b"C C", build_from_graph("C C | 0-1"))]
    #[case::ws_intertoken_tabs_flags(b"C\tC", build_from_graph("C C | 0-1"))]
    #[case::ws_newlines_flags(b"C\nC", build_from_graph("C C | 0-1"))]
    #[case::line_comment_flags(b"C// x\nC", build_from_graph("C C | 0-1"))]
    #[case::block_comment_flags(b"C/* x */C", build_from_graph("C C | 0-1"))]
    #[case::block_comment_multiline_flags(b"C/* x\n y */C", build_from_graph("C C | 0-1"))]
    #[case::eoi_blank_line(b"C\n\nC", build_from_graph("C |"))]
    #[case::eoi_blank_line_crlf(b"C\r\n\r\nC", build_from_graph("C |"))]
    #[case::eoi_blank_line_with_comment(b"C\n/* comment */\n\nC", build_from_graph("C |"))]
    fn m6_whitespace_lenient(#[case] input: &[u8], #[case] expected: Molecule) {
        let flags = SmilesParseFlags::INTERTOKEN_WS
            | SmilesParseFlags::COMMENTS
            | SmilesParseFlags::EXPLICIT_EOI;
        let res = parse_smiles_inner(input, flags);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rstest]
    #[case::split_halogen_ws(b"C l", M6Error::UnsupportedToken { pos: 2 })]
    #[case::percent_ring_ws_split(b"C% 12", M6Error::RingIndexInvalid { pos: 1 })]
    #[case::percent_ring_nl_split(b"C%\n12", M6Error::RingIndexInvalid { pos: 1 })]
    #[case::unterminated_block_comment(b"C/* x", M6Error::UnterminatedBlockComment { pos: 1 })]
    #[case::bracket_inner_ws(b"[ C ]", M6Error::InvalidBracket { pos: 0 })]
    fn m6_whitespace_lenient_invalid(#[case] input: &[u8], #[case] expected: M6Error) {
        let flags = SmilesParseFlags::INTERTOKEN_WS | SmilesParseFlags::COMMENTS | SmilesParseFlags::EXPLICIT_EOI;
        let res = parse_smiles_inner(input, flags);
        assert!(res.is_err(), "{:?} should have failed", input);
        let err = res.unwrap_err();
        assert_eq!(err, expected);
    }
}
