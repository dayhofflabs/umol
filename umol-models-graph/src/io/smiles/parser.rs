//! SMILES parser (FSM-based)

use umol_data::Element;

use crate::io::config::SmilesParseFlags;
use crate::io::ir::builder::{AtomData, BondData, MoleculeBuilder};
use crate::io::ir::{BondDir, BondOrder, Molecule, Chirality};

// Bracket parsing is implemented inline in this module to avoid re-parsing and allocations.

#[derive(Debug, Clone, PartialEq)]
pub enum ParseError {
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
    FieldOutsideBracket { pos: usize },
    BracketDuplicateField { pos: usize },
    BracketHOnH { pos: usize },
    GroupLeadingConnector { pos: usize },
}

// Public entrypoint: strict OpenSMILES
pub fn parse_smiles(input: &[u8]) -> Result<Molecule, ParseError> {
    let flags = SmilesParseFlags::STRICT_OPENSMILES;
    parse_smiles_inner(input, flags)
}

// Flags-aware inner parser
pub fn parse_smiles_inner(input: &[u8], flags: SmilesParseFlags) -> Result<Molecule, ParseError> {
    let allow_ws = flags.contains(SmilesParseFlags::INTERTOKEN_WS);
    let allow_comments = flags.contains(SmilesParseFlags::COMMENTS);
    let use_eoi = flags.contains(SmilesParseFlags::EXPLICIT_EOI);
    let _record_lint = flags.contains(SmilesParseFlags::LINT_SIDECHANNEL);

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
        return parse_core(&input[start..end], flags);
    }

    parse_core(input, flags)
}

fn is_digit(b: u8) -> bool { (b'0'..=b'9').contains(&b) }

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
    Branch { base: u32, had_atom: bool, open_pos: usize },
    Group { had_atom: bool, open_pos: usize },
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

fn parse_bracket_inner_bytes(inner: &[u8], pos_offset: usize)
    -> Result<(Option<Element>, bool, Option<u32>, Option<i32>, Option<u32>, Option<u8>, Option<Chirality>, bool), ParseError>
{
    if inner.last() == Some(&b':') { return Err(ParseError::BracketEmptyClass { pos: pos_offset }); }
    // Detect H followed by two digits
    if inner.windows(3).any(|w| w[0] == b'H' && w[1].is_ascii_digit() && w[2].is_ascii_digit()) {
        return Err(ParseError::BracketHCountTwoDigits { pos: pos_offset });
    }

    let n = inner.len();
    let mut i = 0usize;

    // 1) Optional isotope (one or more digits)
    let mut isotope: Option<u32> = None;
    let start_digits = i;
    while i < n && inner[i].is_ascii_digit() { i += 1; }
    if i > start_digits {
        let mut v: u32 = 0;
        for &b in &inner[start_digits..i] { v = v.saturating_mul(10).saturating_add((b - b'0') as u32); }
        isotope = Some(v);
    }

    // 2) Element symbol or wildcard '*'
    let mut element: Option<Element> = None;
    let mut aromatic = false;
    let mut unknown_symbol = false;
    if i < n && inner[i] == b'*' {
        element = None;
        unknown_symbol = true;
        i += 1;
    } else if i < n && inner[i].is_ascii_alphabetic() {
        // Hybrid fast path + fallback
        let b0 = inner[i];
        if b0.is_ascii_uppercase() {
            // Fast path: common aliphatic
            match b0 {
                b'C' => {
                    if i + 1 < n && inner[i + 1].is_ascii_lowercase() {
                        if inner[i + 1] == b'l' { element = Some(Element::Cl); i += 2; }
                        else if let Some(e) = umol_data::Element::from_symbol_bytes(&inner[i..i + 2]) { element = Some(e); i += 2; }
                        else { element = Some(Element::C); i += 1; }
                    } else { element = Some(Element::C); i += 1; }
                }
                b'B' => {
                    if i + 1 < n && inner[i + 1].is_ascii_lowercase() {
                        if inner[i + 1] == b'r' { element = Some(Element::Br); i += 2; }
                        else if let Some(e) = umol_data::Element::from_symbol_bytes(&inner[i..i + 2]) { element = Some(e); i += 2; }
                        else { element = Some(Element::B); i += 1; }
                    } else { element = Some(Element::B); i += 1; }
                }
                b'N' => {
                    if i + 1 < n && inner[i + 1].is_ascii_lowercase() {
                        if let Some(e) = umol_data::Element::from_symbol_bytes(&inner[i..i + 2]) { element = Some(e); i += 2; }
                        else { element = Some(Element::N); i += 1; }
                    } else { element = Some(Element::N); i += 1; }
                }
                b'O' => {
                    if i + 1 < n && inner[i + 1].is_ascii_lowercase() {
                        if let Some(e) = umol_data::Element::from_symbol_bytes(&inner[i..i + 2]) { element = Some(e); i += 2; }
                        else { element = Some(Element::O); i += 1; }
                    } else { element = Some(Element::O); i += 1; }
                }
                b'S' => {
                    if i + 1 < n && inner[i + 1].is_ascii_lowercase() {
                        if let Some(e) = umol_data::Element::from_symbol_bytes(&inner[i..i + 2]) { element = Some(e); i += 2; }
                        else { element = Some(Element::S); i += 1; }
                    } else { element = Some(Element::S); i += 1; }
                }
                b'P' => {
                    if i + 1 < n && inner[i + 1].is_ascii_lowercase() {
                        if let Some(e) = umol_data::Element::from_symbol_bytes(&inner[i..i + 2]) { element = Some(e); i += 2; }
                        else { element = Some(Element::P); i += 1; }
                    } else { element = Some(Element::P); i += 1; }
                }
                b'F' => {
                    if i + 1 < n && inner[i + 1].is_ascii_lowercase() {
                        if let Some(e) = umol_data::Element::from_symbol_bytes(&inner[i..i + 2]) { element = Some(e); i += 2; }
                        else { element = Some(Element::F); i += 1; }
                    } else { element = Some(Element::F); i += 1; }
                }
                b'I' => {
                    if i + 1 < n && inner[i + 1].is_ascii_lowercase() {
                        if let Some(e) = umol_data::Element::from_symbol_bytes(&inner[i..i + 2]) { element = Some(e); i += 2; }
                        else { element = Some(Element::I); i += 1; }
                    } else { element = Some(Element::I); i += 1; }
                }
                _ => {
                    // Fallback: full periodic table (aliphatic only)
                    let mut consumed = 0usize;
                    if i + 1 < n && inner[i + 1].is_ascii_alphabetic() {
                        if let Some(e) = umol_data::Element::from_symbol_bytes(&inner[i..i + 2]) {
                            element = Some(e);
                            consumed = 2;
                        }
                    }
                    if consumed == 0 {
                        if let Some(e) = umol_data::Element::from_symbol_bytes(&inner[i..i + 1]) {
                            element = Some(e);
                            consumed = 1;
                        } else {
                            return Err(ParseError::InvalidBracket { pos: pos_offset });
                        }
                    }
                    i += consumed;
                }
            }
            aromatic = false;
        } else {
            // Lowercase start: aromatic only if in allowed set
            match b0 {
                b'b' => { element = Some(Element::B); aromatic = true; i += 1; }
                b'c' => { element = Some(Element::C); aromatic = true; i += 1; }
                b'n' => { element = Some(Element::N); aromatic = true; i += 1; }
                b'o' => { element = Some(Element::O); aromatic = true; i += 1; }
                b'p' => { element = Some(Element::P); aromatic = true; i += 1; }
                b's' => {
                    // 's' (S) or 'se' (Se) aromatic
                    if i + 1 < n && inner[i + 1] == b'e' { element = Some(Element::Se); aromatic = true; i += 2; }
                    else { element = Some(Element::S); aromatic = true; i += 1; }
                }
                b'a' => {
                    // 'as' (As) aromatic
                    if i + 1 < n && inner[i + 1] == b's' { element = Some(Element::As); aromatic = true; i += 2; }
                    else { return Err(ParseError::InvalidBracket { pos: pos_offset }); }
                }
                _ => { return Err(ParseError::InvalidBracket { pos: pos_offset }); }
            }
        }
    } else {
        // Neither '*' nor element
        return Err(ParseError::InvalidBracket { pos: pos_offset });
    }

    // 3) Tail fields in any order
    let mut charge: Option<i32> = None;
    let mut class_num: Option<u32> = None;
    let mut hcount: Option<u8> = None;
    let mut chir: Option<Chirality> = None;

    while i < n {
        let b0 = inner[i];
        match b0 {
            b'H' => {
                if element == Some(Element::H) { return Err(ParseError::BracketHOnH { pos: pos_offset }); }
                if hcount.is_some() { return Err(ParseError::BracketDuplicateField { pos: pos_offset }); }
                let mut val: u32 = 1; // default H
                if i + 1 < n && inner[i + 1].is_ascii_digit() {
                    val = (inner[i + 1] - b'0') as u32; i += 1;
                }
                hcount = Some((val as u8).min(u8::MAX));
                i += 1;
            }
            b'+' | b'-' => {
                if charge.is_some() { return Err(ParseError::BracketDuplicateField { pos: pos_offset }); }
                let sign = if b0 == b'+' { 1 } else { -1 };
                let mut j = i + 1; let mut val: i32 = 0; let mut cnt = 0;
                if j < n && inner[j] == b'+' && b0 == b'+' { charge = Some(2); i = j + 1; continue; }
                if j < n && inner[j] == b'-' && b0 == b'-' { charge = Some(-2); i = j + 1; continue; }
                while j < n && inner[j].is_ascii_digit() && cnt < 2 { val = val.saturating_mul(10) + (inner[j] - b'0') as i32; j += 1; cnt += 1; }
                if cnt == 0 { val = 1; }
                charge = Some(val * sign);
                i = j;
            }
            b':' => {
                if class_num.is_some() { return Err(ParseError::BracketDuplicateField { pos: pos_offset }); }
                if i + 1 >= n || !inner[i + 1].is_ascii_digit() { return Err(ParseError::BracketEmptyClass { pos: pos_offset }); }
                let mut j = i + 1; let mut v: u32 = 0;
                while j < n && inner[j].is_ascii_digit() { v = v.saturating_mul(10).saturating_add((inner[j]-b'0') as u32); j += 1; }
                class_num = Some(v);
                i = j;
            }
            b'@' => {
                if chir.is_some() { return Err(ParseError::BracketDuplicateField { pos: pos_offset }); }
                // @@
                if i + 1 < n && inner[i + 1] == b'@' { chir = Some(Chirality::CounterClockwise); i += 2; continue; }
                // @THn
                if i + 3 < n && inner[i + 1] == b'T' && inner[i + 2] == b'H' && inner[i + 3].is_ascii_digit() {
                    chir = Some(Chirality::Tetrahedral { arr: (inner[i + 3] - b'0') as u32 }); i += 4; continue;
                }
                // @AL[12]
                if i + 3 < n && inner[i + 1] == b'A' && inner[i + 2] == b'L' && (inner[i + 3] == b'1' || inner[i + 3] == b'2') {
                    chir = Some(Chirality::Allenal { arr: (inner[i + 3] - b'0') as u32 }); i += 4; continue;
                }
                // @SP[123]
                if i + 3 < n && inner[i + 1] == b'S' && inner[i + 2] == b'P' && (inner[i + 3] == b'1' || inner[i + 3] == b'2' || inner[i + 3] == b'3') {
                    chir = Some(Chirality::SquarePlanar { arr: (inner[i + 3] - b'0') as u32 }); i += 4; continue;
                }
                // @TBn (n: first digit only)
                if i + 3 <= n && i + 3 - 1 < n && inner[i + 1] == b'T' && inner[i + 2] == b'B' && i + 3 < n && inner[i + 3].is_ascii_digit() {
                    chir = Some(Chirality::TrigonalBipyramidal { arr: (inner[i + 3] - b'0') as u32 }); i += 4; continue;
                }
                // @OHn (n: first digit only)
                if i + 3 <= n && i + 3 - 1 < n && inner[i + 1] == b'O' && inner[i + 2] == b'H' && i + 3 < n && inner[i + 3].is_ascii_digit() {
                    chir = Some(Chirality::Octahedral { arr: (inner[i + 3] - b'0') as u32 }); i += 4; continue;
                }
                // '@' alone
                chir = Some(Chirality::Clockwise);
                i += 1;
            }
            _ => { return Err(ParseError::InvalidBracket { pos: pos_offset }); }
        }
    }

    Ok((element, aromatic, isotope, charge, class_num, hcount, chir, unknown_symbol))
}

fn truncate_at_eoi(input: &[u8], allow_comments: bool) -> usize {
    let n = input.len();
    let mut i = 0usize;
    let mut line_start = 0usize;
    let mut had_content = false;
    while i < n {
        let b0 = input[i];
        if b0 == b' ' || b0 == b'\t' { i += 1; continue; }
        if allow_comments && b0 == b'/' && i + 1 < n && input[i + 1] == b'/' {
            i += 2; while i < n && input[i] != b'\n' && input[i] != b'\r' { i += 1; }
        }
        if allow_comments && b0 == b'/' && i + 1 < n && input[i + 1] == b'*' {
            i += 2; while i + 1 < n { if input[i] == b'*' && input[i + 1] == b'/' { i += 2; break; } i += 1; } continue;
        }
        if b0 == b'\r' {
            if !had_content { return line_start; }
            i += 1; if i < n && input[i] == b'\n' { i += 1; }
            line_start = i; had_content = false; continue;
        }
        if b0 == b'\n' { if !had_content { return line_start; } i += 1; line_start = i; had_content = false; continue; }
        had_content = true; i += 1;
    }
    n
}

fn parse_core(input: &[u8], flags: SmilesParseFlags) -> Result<Molecule, ParseError> {
    let allow_ws = flags.contains(SmilesParseFlags::INTERTOKEN_WS);
    let allow_comments = flags.contains(SmilesParseFlags::COMMENTS);
    let _record_lint = flags.contains(SmilesParseFlags::LINT_SIDECHANNEL);

    let mut i = 0usize; let n = input.len();
    let mut builder = MoleculeBuilder::with_capacity(n.max(1), n.max(1).saturating_sub(1));
    let mut last_atom_idx: Option<u32> = None; let mut prev_atom_idx: Option<u32> = None;
    let mut pending_bond: Option<(BondOrder, Option<BondDir>, usize)> = None; let mut last_aromatic: bool = false;
    let mut pstack: Vec<Frame> = Vec::new(); let mut ring_table: [Option<OpenRing>; 100] = [None; 100];
    let mut just_closed_group: bool = false;

    // Gated side-channel (not exposed yet): ring sequence and percent padding
    let mut _ring_sequence: Option<Vec<(u32, usize)>> = if _record_lint { Some(Vec::with_capacity(8)) } else { None };
    let mut _first_ring_digit: Option<u32> = None;
    let mut _percent_padded: bool = false;

    while i < n {
        let b0 = input[i];
        if allow_comments && b0 == b'/' && i + 1 < n {
            let b1 = input[i + 1];
            if b1 == b'/' { i += 2; while i < n && input[i] != b'\n' { i += 1; } continue; }
            if b1 == b'*' { let start_pos = i; i += 2; let mut closed = false; while i + 1 < n { if input[i] == b'*' && input[i + 1] == b'/' { i += 2; closed = true; break; } i += 1; } if !closed { return Err(ParseError::UnterminatedBlockComment { pos: start_pos }); } continue; }
        }
        if matches!(b0, b' ' | b'\t' | b'\n' | b'\r') { if allow_ws { i += 1; continue; } return Err(ParseError::InvalidWhitespace { pos: i }); }
        if b0 != b'(' { just_closed_group = false; }
        if b0 == b'(' { if let Some((_, _, pos)) = pending_bond { return Err(ParseError::TrailingBond { pos }); } if just_closed_group { last_atom_idx=None; prev_atom_idx=None; pstack.push(Frame::Group{had_atom:false,open_pos:i}); just_closed_group=false; } else { match last_atom_idx { Some(idx)=>pstack.push(Frame::Branch{base:idx,had_atom:false,open_pos:i}), None=>pstack.push(Frame::Group{had_atom:false,open_pos:i}), } } i+=1; continue; }
        if b0 == b')' { if let Some((_, _, pos)) = pending_bond { return Err(ParseError::TrailingBond { pos }); } let Some(frame)=pstack.pop() else { return Err(ParseError::UnbalancedBranchClose { pos: i }); }; match frame { Frame::Branch { base, had_atom, .. } => { if !had_atom { return Err(ParseError::EmptyBranch { pos: i }); } last_atom_idx=Some(base); prev_atom_idx=None; } Frame::Group { had_atom, open_pos, .. } => { if !had_atom { if i + 1 != n { return Err(ParseError::EmptyGroup { pos: i }); } if i>0 && input[i-1]==b'.' { return Err(ParseError::LeadingDot { pos: i }); } if open_pos != 0 { return Err(ParseError::EmptyGroup { pos: i }); } last_atom_idx=None; prev_atom_idx=None; just_closed_group=false; } else { just_closed_group=true; if pstack.is_empty() && i + 1 != n { let next = input[i + 1]; if next != b'.' { return Err(ParseError::TopLevelGroupTrailing { pos: i }); } } } } } i+=1; continue; }
        if b0 == b'.' { if let Some((_, _, pos)) = pending_bond { return Err(ParseError::TrailingBond { pos }); } if i==0 { return Err(ParseError::LeadingDot { pos: i }); } if let Some(Frame::Group { had_atom: false, .. }) = pstack.last() { return Err(ParseError::GroupLeadingConnector { pos: i }); } if i + 1 == n { return Err(ParseError::TrailingDot { pos: i }); } if input[i + 1] == b'.' { return Err(ParseError::ConsecutiveDot { pos: i }); } last_atom_idx=None; prev_atom_idx=None; last_aromatic=false; i+=1; continue; }
        if is_digit(b0) { if last_atom_idx.is_none() { return Err(ParseError::LeadingRing { pos: i }); } let idx: usize = (b0 - b'0') as usize; if let Some(seq) = _ring_sequence.as_mut() { let d = idx as u32; if _first_ring_digit.is_none() { _first_ring_digit = Some(d); } seq.push((d, i)); } let bond = pending_bond.take(); let (order_opt, dir_opt) = bond.map_or((None,None), |(o,d,_)|(Some(o),d)); let entry=&mut ring_table[idx]; match entry.take() { None=>{ *entry=Some(OpenRing{ atom_id:last_atom_idx.unwrap(), order:order_opt, dir:dir_opt, open_pos:i, open_aromatic:last_aromatic }); } Some(open)=>{ let b=last_atom_idx.unwrap(); if open.atom_id==b { return Err(ParseError::RingSelfLoop { pos: i }); } if prev_atom_idx==Some(open.atom_id) { return Err(ParseError::RingTwoMember { pos: i }); } if let (Some(d1),Some(d2))=(open.dir,dir_opt) { if d1!=d2 { return Err(ParseError::RingBondDirConflict{pos:i, open_pos: open.open_pos}); } } if let (Some(o1),Some(o2))=(open.order,order_opt) { if o1!=o2 { return Err(ParseError::RingBondOrderConflict{pos:i, open_pos: open.open_pos}); } } if open.dir.is_some() || dir_opt.is_some() { let ord = open.order.or(order_opt).unwrap_or(BondOrder::Single); if ord != BondOrder::Single { return Err(ParseError::RingBondOrderConflict{pos:i, open_pos: open.open_pos}); } }
            let mut final_order = match (open.order, order_opt) { (Some(o1), Some(o2)) => { if o1==o2 { o1 } else { o2 } }, (Some(o), None) | (None, Some(o)) => o, (None, None) => BondOrder::Single };
            let final_dir = open.dir.or(dir_opt); let a=open.atom_id; let b=last_atom_idx.unwrap(); if final_order==BondOrder::Single && open.open_aromatic && last_aromatic { final_order=BondOrder::Aromatic; } if a != b { builder.on_bond(a,b, BondData{order:final_order,dir:final_dir}); } } } i+=1; continue; }
        if b0 == b'%' { if i + 2 >= n || !is_digit(input[i+1]) || !is_digit(input[i+2]) { return Err(ParseError::RingIndexInvalid { pos: i }); } if last_atom_idx.is_none() { return Err(ParseError::LeadingRing { pos: i }); } if let Some(seq) = _ring_sequence.as_mut() { if input[i+1] == b'0' && (input[i+2] >= b'1' && input[i+2] <= b'9') { _percent_padded = true; } let d = ((input[i+1]-b'0') as u32)*10 + (input[i+2]-b'0') as u32; if _first_ring_digit.is_none() { _first_ring_digit = Some(d); } seq.push((d, i)); }
            let idx: usize = ((input[i+1]-b'0') as usize)*10 + (input[i+2]-b'0') as usize; let bond=pending_bond.take(); let (order_opt, dir_opt) = bond.map_or((None,None), |(o,d,_)|(Some(o),d)); let entry=&mut ring_table[idx]; match entry.take() { None=>{ *entry=Some(OpenRing{ atom_id:last_atom_idx.unwrap(), order:order_opt, dir:dir_opt, open_pos:i, open_aromatic:last_aromatic }); } Some(open)=>{ let b=last_atom_idx.unwrap(); if open.atom_id==b { return Err(ParseError::RingSelfLoop { pos: i }); } if prev_atom_idx==Some(open.atom_id) { return Err(ParseError::RingTwoMember { pos: i }); } if let (Some(d1),Some(d2))=(open.dir,dir_opt) { if d1!=d2 { return Err(ParseError::RingBondDirConflict{pos:i, open_pos: open.open_pos}); } } if let (Some(o1),Some(o2))=(open.order,order_opt) { if o1!=o2 { return Err(ParseError::RingBondOrderConflict{pos:i, open_pos: open.open_pos}); } } if open.dir.is_some() || dir_opt.is_some() { let ord = open.order.or(order_opt).unwrap_or(BondOrder::Single); if ord != BondOrder::Single { return Err(ParseError::RingBondOrderConflict{pos:i, open_pos: open.open_pos}); } }
            let mut final_order = match (open.order, order_opt) { (Some(o1), Some(o2)) => { if o1==o2 { o1 } else { o2 } }, (Some(o), None) | (None, Some(o)) => o, (None, None) => BondOrder::Single };
            let final_dir = open.dir.or(dir_opt); let a=open.atom_id; let b=last_atom_idx.unwrap(); if final_order==BondOrder::Single && open.open_aromatic && last_aromatic { final_order=BondOrder::Aromatic; } if a != b { builder.on_bond(a,b, BondData{order:final_order,dir:final_dir}); } } } i+=3; continue; }
        if matches!(b0, b'-' | b'=' | b'#' | b'$' | b':' | b'/' | b'\\') { if pending_bond.is_some() { return Err(ParseError::ConsecutiveBond { pos: i }); } if last_atom_idx.is_none() { if let Some(Frame::Group { had_atom: false, .. }) = pstack.last() { return Err(ParseError::GroupLeadingConnector { pos: i }); } return Err(ParseError::LeadingBond { pos: i }); } let (order,dir)=map_bond(b0); pending_bond=Some((order,dir,i)); i+=1; continue; }
        if b0 == b'[' { let start=i+1; let mut j=start; while j < n && input[j] != b']' { j+=1; } if j >= n { return Err(ParseError::UnbalancedOpenBracket { pos: i }); } let inner=&input[start..j]; let (elem_opt, aromatic, iso_opt, charge_opt, class_opt, h_opt, chir_opt, unknown) = parse_bracket_inner_bytes(inner, i)?; let (element, aromatic) = match elem_opt { Some(e) => (e, aromatic), None => (Element::C, false) }; let atom = AtomData { element, isotope: iso_opt, charge: charge_opt, hydrogen_count: h_opt.map(|h| (h as u8).min(u8::MAX)), class: class_opt, aromatic, implicit_h: false, chirality: chir_opt, unknown_symbol: unknown, };
            let curr = builder.on_atom(atom); if let Some(last)=last_atom_idx { if let Some((order,dir,_))=pending_bond.take() { builder.on_bond(last,curr,BondData{order,dir}); } else { if last_aromatic && aromatic { builder.on_bond(last,curr,BondData{order:BondOrder::Aromatic, dir:None}); } else { builder.on_bond_single_fast(last,curr); } } prev_atom_idx = Some(last); } last_atom_idx = Some(curr); last_aromatic = aromatic; if let Some(top)=pstack.last_mut() { match top { Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => { *had_atom = true } } } i = j + 1; continue; }
        if b0 == b'C' { if i + 1 < n && input[i+1] == b'l' { let curr = builder.on_atom_fast(Element::Cl, true, false); if let Some(last)=last_atom_idx { if let Some((order,dir,_))=pending_bond.take() { builder.on_bond(last,curr,BondData{order,dir}); } else { builder.on_bond_single_fast(last,curr); } } last_atom_idx=Some(curr); last_aromatic=false; if let Some(top)=pstack.last_mut() { match top { Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => { *had_atom = true } } } i += 2; continue; } let curr = builder.on_atom_fast(Element::C, true, false); if let Some(last)=last_atom_idx { if let Some((order,dir,_))=pending_bond.take() { builder.on_bond(last,curr,BondData{order,dir}); } else { builder.on_bond_single_fast(last,curr); } prev_atom_idx = Some(last); } last_atom_idx=Some(curr); last_aromatic=false; if let Some(top)=pstack.last_mut() { match top { Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => { *had_atom = true } } } i += 1; continue; }
        if b0 == b'B' { if i + 1 < n && input[i+1] == b'r' { let curr = builder.on_atom_fast(Element::Br, true, false); if let Some(last)=last_atom_idx { if let Some((order,dir,_))=pending_bond.take() { builder.on_bond(last,curr,BondData{order,dir}); } else { builder.on_bond_single_fast(last,curr); } } last_atom_idx=Some(curr); last_aromatic=false; if let Some(top)=pstack.last_mut() { match top { Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => { *had_atom = true } } } i += 2; continue; } let curr = builder.on_atom_fast(Element::B, true, false); if let Some(last)=last_atom_idx { if let Some((order,dir,_))=pending_bond.take() { builder.on_bond(last,curr,BondData{order,dir}); } else { builder.on_bond_single_fast(last,curr); } prev_atom_idx = Some(last); } last_atom_idx=Some(curr); last_aromatic=false; if let Some(top)=pstack.last_mut() { match top { Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => { *had_atom = true } } } i += 1; continue; }
        let elem = match b0 { b'N' => Some(Element::N), b'O' => Some(Element::O), b'P' => Some(Element::P), b'S' => Some(Element::S), b'F' => Some(Element::F), b'I' => Some(Element::I), _ => None };
        if let Some(element) = elem { let curr = builder.on_atom_fast(element, true, false); if let Some(last)=last_atom_idx { if let Some((order,dir,_))=pending_bond.take() { builder.on_bond(last,curr,BondData{order,dir}); } else { builder.on_bond_single_fast(last,curr); } prev_atom_idx = Some(last); } last_atom_idx=Some(curr); last_aromatic=false; if let Some(top)=pstack.last_mut() { match top { Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => { *had_atom = true } } } i += 1; continue; }
        if matches!(b0, b'b' | b'c' | b'n' | b'o' | b'p' | b's') { let element = match b0 { b'b' => Element::B, b'c' => Element::C, b'n' => Element::N, b'o' => Element::O, b'p' => Element::P, _ => Element::S }; let curr = builder.on_atom_fast(element, true, true); if let Some(last)=last_atom_idx { if let Some((order,dir,_))=pending_bond.take() { builder.on_bond(last,curr,BondData{order,dir}); } else { if last_aromatic { builder.on_bond(last,curr,BondData{order:BondOrder::Aromatic, dir:None}); } else { builder.on_bond_single_fast(last,curr); } } prev_atom_idx = Some(last); } last_atom_idx=Some(curr); last_aromatic=true; if let Some(top)=pstack.last_mut() { match top { Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => { *had_atom = true } } } i += 1; continue; }
        if b0 == b'*' { let atom = AtomData { element: Element::C, isotope: Some(0), charge: Some(0), hydrogen_count: Some(0), class: None, aromatic: false, implicit_h: false, chirality: None, unknown_symbol: true }; let curr = builder.on_atom(atom); if let Some(last)=last_atom_idx { if let Some((order,dir,_))=pending_bond.take() { builder.on_bond(last,curr,BondData{order,dir}); } else { builder.on_bond_single_fast(last,curr); } prev_atom_idx = Some(last); } last_atom_idx=Some(curr); last_aromatic=false; if let Some(top)=pstack.last_mut() { match top { Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => { *had_atom = true } } } i += 1; continue; }
        if b0 == b']' { return Err(ParseError::UnbalancedCloseBracket { pos: i }); }
        // Bracket-only fields outside bracket
        if b0 == b'@' || b0 == b'+' {
            return Err(ParseError::FieldOutsideBracket { pos: i });
        }
        return Err(ParseError::UnsupportedToken { pos: i });
    }

    if pending_bond.is_some() { let (_,_,pos)=pending_bond.unwrap(); return Err(ParseError::TrailingBond { pos }); }
    if !pstack.is_empty() { let pos = match pstack.last().unwrap() { Frame::Branch{open_pos,..} | Frame::Group{open_pos,..} => *open_pos }; return Err(ParseError::UnbalancedBranchOpen { pos }); }
    let mut last_open: Option<usize> = None; for entry in ring_table.iter().flatten() { match last_open { None => last_open = Some(entry.open_pos), Some(p) => { if entry.open_pos > p { last_open = Some(entry.open_pos) } } } }
    if let Some(pos_open) = last_open { return Err(ParseError::RingUnclosed { open_pos: pos_open }); }
    let mut mols = builder.finish(); Ok(mols.pop().unwrap_or_default())
}

#[cfg(test)]
mod tests;