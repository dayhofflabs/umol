//! SMILES parser (FSM-based)

use smallvec::SmallVec;
use umol_data::Element;

use crate::io::config::SmilesParseFlags;
use crate::io::ir::builder::{AtomData, BondData, MoleculeBuilder};
use crate::io::ir::{BondDir, BondOrder, Molecule};
use crate::io::smiles::parser::utils::BracketField;

pub mod utils;

#[derive(Debug, Clone, PartialEq)]
pub enum M6Error {
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

// Public entrypoint: strict OpenSMILES
pub fn parse_smiles(input: &[u8]) -> Result<Molecule, M6Error> {
    let flags = SmilesParseFlags::STRICT_OPENSMILES;
    parse_smiles_inner(input, flags)
}

// Flags-aware inner parser
pub fn parse_smiles_inner(input: &[u8], flags: SmilesParseFlags) -> Result<Molecule, M6Error> {
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
            return Err(M6Error::InvalidWhitespace { pos: 0 });
        }
        let mut end = input.len();
        while end > 0 && matches!(input[end - 1], b' ' | b'\t' | b'\n' | b'\r') {
            end -= 1;
        }
        for (k, b) in input[start..end].iter().enumerate() {
            if matches!(*b, b' ' | b'\t' | b'\n' | b'\r') {
                return Err(M6Error::InvalidWhitespace { pos: start + k });
            }
        }
        return m6_parse_core(&input[start..end], flags);
    }

    m6_parse_core(input, flags)
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

fn parse_bracket_m6(
    inner: &[u8],
    pos_offset: usize,
) -> Result<(Option<Element>, Option<u32>, SmallVec<[BracketField; 4]>), M6Error> {
    let bytes = inner;
    if bytes.last() == Some(&b':') { return Err(M6Error::BracketEmptyClass { pos: pos_offset }); }
    if bytes.windows(3).any(|w| w[0] == b'H' && w[1].is_ascii_digit() && w[2].is_ascii_digit()) {
        return Err(M6Error::BracketHCountTwoDigits { pos: pos_offset });
    }
    let s = match std::str::from_utf8(inner) { Ok(v) => v, Err(_) => return Err(M6Error::InvalidBracket { pos: pos_offset }) };
    if !crate::io::smiles::parser::utils::is_valid_bracket_inner(s) {
        return Err(M6Error::InvalidBracket { pos: pos_offset });
    }
    let (elem_opt, iso_opt, tails) = crate::io::smiles::parser::utils::parse_bracket(s);
    Ok((elem_opt, iso_opt, tails))
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

fn m6_parse_core(input: &[u8], flags: SmilesParseFlags) -> Result<Molecule, M6Error> {
    let allow_ws = flags.contains(SmilesParseFlags::INTERTOKEN_WS);
    let allow_comments = flags.contains(SmilesParseFlags::COMMENTS);

    let mut i = 0usize; let n = input.len();
    let mut builder = MoleculeBuilder::with_capacity(n.max(1), n.max(1).saturating_sub(1));
    let mut last_atom_idx: Option<u32> = None; let mut prev_atom_idx: Option<u32> = None;
    let mut pending_bond: Option<(BondOrder, Option<BondDir>, usize)> = None; let mut last_aromatic: bool = false;
    let mut pstack: Vec<Frame> = Vec::new(); let mut ring_table: [Option<OpenRing>; 100] = [None; 100];
    let mut just_closed_group: bool = false;

    while i < n {
        let b0 = input[i];
        if allow_comments && b0 == b'/' && i + 1 < n {
            let b1 = input[i + 1];
            if b1 == b'/' { i += 2; while i < n && input[i] != b'\n' { i += 1; } continue; }
            if b1 == b'*' { let start_pos = i; i += 2; let mut closed = false; while i + 1 < n { if input[i] == b'*' && input[i + 1] == b'/' { i += 2; closed = true; break; } i += 1; } if !closed { return Err(M6Error::UnterminatedBlockComment { pos: start_pos }); } continue; }
        }
        if matches!(b0, b' ' | b'\t' | b'\n' | b'\r') { if allow_ws { i += 1; continue; } return Err(M6Error::InvalidWhitespace { pos: i }); }
        if b0 != b'(' { just_closed_group = false; }
        if b0 == b'(' { if let Some((_, _, pos)) = pending_bond { return Err(M6Error::TrailingBond { pos }); } if just_closed_group { last_atom_idx=None; prev_atom_idx=None; pstack.push(Frame::Group{had_atom:false,open_pos:i}); just_closed_group=false; } else { match last_atom_idx { Some(idx)=>pstack.push(Frame::Branch{base:idx,had_atom:false,open_pos:i}), None=>pstack.push(Frame::Group{had_atom:false,open_pos:i}), } } i+=1; continue; }
        if b0 == b')' { if let Some((_, _, pos)) = pending_bond { return Err(M6Error::TrailingBond { pos }); } let Some(frame)=pstack.pop() else { return Err(M6Error::UnbalancedBranchClose { pos: i }); }; match frame { Frame::Branch { base, had_atom, .. } => { if !had_atom { return Err(M6Error::EmptyBranch { pos: i }); } last_atom_idx=Some(base); prev_atom_idx=None; } Frame::Group { had_atom, open_pos, .. } => { if !had_atom { if i + 1 != n { return Err(M6Error::EmptyGroup { pos: i }); } if i>0 && input[i-1]==b'.' { return Err(M6Error::LeadingDot { pos: i }); } if open_pos != 0 { return Err(M6Error::EmptyGroup { pos: i }); } last_atom_idx=None; prev_atom_idx=None; just_closed_group=false; } else { just_closed_group=true; if pstack.is_empty() && i + 1 != n { let next = input[i + 1]; if next != b'.' { return Err(M6Error::TopLevelGroupTrailing { pos: i }); } } } } } i+=1; continue; }
        if b0 == b'.' { if let Some((_, _, pos)) = pending_bond { return Err(M6Error::TrailingBond { pos }); } if i==0 { return Err(M6Error::LeadingDot { pos: i }); } if i + 1 == n { return Err(M6Error::TrailingDot { pos: i }); } if input[i + 1] == b'.' { return Err(M6Error::ConsecutiveDot { pos: i }); } last_atom_idx=None; prev_atom_idx=None; last_aromatic=false; i+=1; continue; }
        if is_digit(b0) { if last_atom_idx.is_none() { return Err(M6Error::LeadingRing { pos: i }); } let idx: usize = (b0 - b'0') as usize; let bond = pending_bond.take(); let (order_opt, dir_opt) = bond.map_or((None,None), |(o,d,_)|(Some(o),d)); let entry=&mut ring_table[idx]; match entry.take() { None=>{ *entry=Some(OpenRing{ atom_id:last_atom_idx.unwrap(), order:order_opt, dir:dir_opt, open_pos:i, open_aromatic:last_aromatic }); } Some(open)=>{ let b=last_atom_idx.unwrap(); if open.atom_id==b { return Err(M6Error::RingSelfLoop { pos: i }); } if prev_atom_idx==Some(open.atom_id) { return Err(M6Error::RingTwoMember { pos: i }); } if let (Some(d1),Some(d2))=(open.dir,dir_opt) { if d1!=d2 { return Err(M6Error::RingBondDirConflict{pos:i, open_pos: open.open_pos}); } } if let (Some(o1),Some(o2))=(open.order,order_opt) { if o1!=o2 { return Err(M6Error::RingBondOrderConflict{pos:i, open_pos: open.open_pos}); } } if open.dir.is_some() || dir_opt.is_some() { let ord = open.order.or(order_opt).unwrap_or(BondOrder::Single); if ord != BondOrder::Single { return Err(M6Error::RingBondOrderConflict{pos:i, open_pos: open.open_pos}); } }
            let mut final_order = match (open.order, order_opt) { (Some(o1), Some(o2)) => { if o1==o2 { o1 } else { o2 } }, (Some(o), None) | (None, Some(o)) => o, (None, None) => BondOrder::Single };
            let final_dir = open.dir.or(dir_opt); let a=open.atom_id; let b=last_atom_idx.unwrap(); if final_order==BondOrder::Single && open.open_aromatic && last_aromatic { final_order=BondOrder::Aromatic; } if a != b { builder.on_bond(a,b, BondData{order:final_order,dir:final_dir}); } } } i+=1; continue; }
        if b0 == b'%' { if i + 2 >= n || !is_digit(input[i+1]) || !is_digit(input[i+2]) { return Err(M6Error::RingIndexInvalid { pos: i }); } if last_atom_idx.is_none() { return Err(M6Error::LeadingRing { pos: i }); } let idx: usize = ((input[i+1]-b'0') as usize)*10 + (input[i+2]-b'0') as usize; let bond=pending_bond.take(); let (order_opt, dir_opt) = bond.map_or((None,None), |(o,d,_)|(Some(o),d)); let entry=&mut ring_table[idx]; match entry.take() { None=>{ *entry=Some(OpenRing{ atom_id:last_atom_idx.unwrap(), order:order_opt, dir:dir_opt, open_pos:i, open_aromatic:last_aromatic }); } Some(open)=>{ let b=last_atom_idx.unwrap(); if open.atom_id==b { return Err(M6Error::RingSelfLoop { pos: i }); } if prev_atom_idx==Some(open.atom_id) { return Err(M6Error::RingTwoMember { pos: i }); } if let (Some(d1),Some(d2))=(open.dir,dir_opt) { if d1!=d2 { return Err(M6Error::RingBondDirConflict{pos:i, open_pos: open.open_pos}); } } if let (Some(o1),Some(o2))=(open.order,order_opt) { if o1!=o2 { return Err(M6Error::RingBondOrderConflict{pos:i, open_pos: open.open_pos}); } } if open.dir.is_some() || dir_opt.is_some() { let ord = open.order.or(order_opt).unwrap_or(BondOrder::Single); if ord != BondOrder::Single { return Err(M6Error::RingBondOrderConflict{pos:i, open_pos: open.open_pos}); } }
            let mut final_order = match (open.order, order_opt) { (Some(o1), Some(o2)) => { if o1==o2 { o1 } else { o2 } }, (Some(o), None) | (None, Some(o)) => o, (None, None) => BondOrder::Single };
            let final_dir = open.dir.or(dir_opt); let a=open.atom_id; let b=last_atom_idx.unwrap(); if final_order==BondOrder::Single && open.open_aromatic && last_aromatic { final_order=BondOrder::Aromatic; } if a != b { builder.on_bond(a,b, BondData{order:final_order,dir:final_dir}); } } } i+=3; continue; }
        if matches!(b0, b'-' | b'=' | b'#' | b'$' | b':' | b'/' | b'\\') { if pending_bond.is_some() { return Err(M6Error::ConsecutiveBond { pos: i }); } if last_atom_idx.is_none() { return Err(M6Error::LeadingBond { pos: i }); } let (order,dir)=map_bond(b0); pending_bond=Some((order,dir,i)); i+=1; continue; }
        if b0 == b'[' { let start=i+1; let mut j=start; while j < n && input[j] != b']' { j+=1; } if j >= n { return Err(M6Error::UnbalancedOpenBracket { pos: i }); } let inner=&input[start..j]; let (elem_opt, iso_opt, fields) = parse_bracket_m6(inner, i)?; let (element, aromatic) = match elem_opt { Some(e) => { let first = inner.first().copied().unwrap_or_default(); (e, first.is_ascii_lowercase()) } None => (Element::C, false) }; let mut atom = AtomData { element, isotope: iso_opt, charge: None, hydrogen_count: None, class: None, aromatic, implicit_h: false, chirality: None, unknown_symbol: elem_opt.is_none(), };
            for f in fields { match f { BracketField::Chiral(ch) => atom.chirality = Some(ch), BracketField::HydrogenCount(h) => { atom.hydrogen_count = Some((h as u8).min(u8::MAX)) }, BracketField::Charge(q) => atom.charge = Some(q), BracketField::Class(c) => atom.class = Some(c), } }
            let curr = builder.on_atom(atom); if let Some(last)=last_atom_idx { if let Some((order,dir,_))=pending_bond.take() { builder.on_bond(last,curr,BondData{order,dir}); } else { if last_aromatic && aromatic { builder.on_bond(last,curr,BondData{order:BondOrder::Aromatic, dir:None}); } else { builder.on_bond_single_fast(last,curr); } } prev_atom_idx = Some(last); } last_atom_idx = Some(curr); last_aromatic = aromatic; if let Some(top)=pstack.last_mut() { match top { Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => { *had_atom = true } } } i = j + 1; continue; }
        if b0 == b'C' { if i + 1 < n && input[i+1] == b'l' { let curr = builder.on_atom_fast(Element::Cl, true, false); if let Some(last)=last_atom_idx { if let Some((order,dir,_))=pending_bond.take() { builder.on_bond(last,curr,BondData{order,dir}); } else { builder.on_bond_single_fast(last,curr); } } last_atom_idx=Some(curr); last_aromatic=false; if let Some(top)=pstack.last_mut() { match top { Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => { *had_atom = true } } } i += 2; continue; } let curr = builder.on_atom_fast(Element::C, true, false); if let Some(last)=last_atom_idx { if let Some((order,dir,_))=pending_bond.take() { builder.on_bond(last,curr,BondData{order,dir}); } else { builder.on_bond_single_fast(last,curr); } prev_atom_idx = Some(last); } last_atom_idx=Some(curr); last_aromatic=false; if let Some(top)=pstack.last_mut() { match top { Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => { *had_atom = true } } } i += 1; continue; }
        if b0 == b'B' { if i + 1 < n && input[i+1] == b'r' { let curr = builder.on_atom_fast(Element::Br, true, false); if let Some(last)=last_atom_idx { if let Some((order,dir,_))=pending_bond.take() { builder.on_bond(last,curr,BondData{order,dir}); } else { builder.on_bond_single_fast(last,curr); } } last_atom_idx=Some(curr); last_aromatic=false; if let Some(top)=pstack.last_mut() { match top { Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => { *had_atom = true } } } i += 2; continue; } let curr = builder.on_atom_fast(Element::B, true, false); if let Some(last)=last_atom_idx { if let Some((order,dir,_))=pending_bond.take() { builder.on_bond(last,curr,BondData{order,dir}); } else { builder.on_bond_single_fast(last,curr); } prev_atom_idx = Some(last); } last_atom_idx=Some(curr); last_aromatic=false; if let Some(top)=pstack.last_mut() { match top { Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => { *had_atom = true } } } i += 1; continue; }
        let elem = match b0 { b'N' => Some(Element::N), b'O' => Some(Element::O), b'P' => Some(Element::P), b'S' => Some(Element::S), b'F' => Some(Element::F), b'I' => Some(Element::I), _ => None };
        if let Some(element) = elem { let curr = builder.on_atom_fast(element, true, false); if let Some(last)=last_atom_idx { if let Some((order,dir,_))=pending_bond.take() { builder.on_bond(last,curr,BondData{order,dir}); } else { builder.on_bond_single_fast(last,curr); } prev_atom_idx = Some(last); } last_atom_idx=Some(curr); last_aromatic=false; if let Some(top)=pstack.last_mut() { match top { Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => { *had_atom = true } } } i += 1; continue; }
        if matches!(b0, b'b' | b'c' | b'n' | b'o' | b'p' | b's') { let element = match b0 { b'b' => Element::B, b'c' => Element::C, b'n' => Element::N, b'o' => Element::O, b'p' => Element::P, _ => Element::S }; let curr = builder.on_atom_fast(element, true, true); if let Some(last)=last_atom_idx { if let Some((order,dir,_))=pending_bond.take() { builder.on_bond(last,curr,BondData{order,dir}); } else { if last_aromatic { builder.on_bond(last,curr,BondData{order:BondOrder::Aromatic, dir:None}); } else { builder.on_bond_single_fast(last,curr); } } prev_atom_idx = Some(last); } last_atom_idx=Some(curr); last_aromatic=true; if let Some(top)=pstack.last_mut() { match top { Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => { *had_atom = true } } } i += 1; continue; }
        if b0 == b'*' { let atom = AtomData { element: Element::C, isotope: Some(0), charge: Some(0), hydrogen_count: Some(0), class: None, aromatic: false, implicit_h: false, chirality: None, unknown_symbol: true }; let curr = builder.on_atom(atom); if let Some(last)=last_atom_idx { if let Some((order,dir,_))=pending_bond.take() { builder.on_bond(last,curr,BondData{order,dir}); } else { builder.on_bond_single_fast(last,curr); } prev_atom_idx = Some(last); } last_atom_idx=Some(curr); last_aromatic=false; if let Some(top)=pstack.last_mut() { match top { Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => { *had_atom = true } } } i += 1; continue; }
        if b0 == b']' { return Err(M6Error::UnbalancedCloseBracket { pos: i }); }
        return Err(M6Error::UnsupportedToken { pos: i });
    }

    if pending_bond.is_some() { let (_,_,pos)=pending_bond.unwrap(); return Err(M6Error::TrailingBond { pos }); }
    if !pstack.is_empty() { let pos = match pstack.last().unwrap() { Frame::Branch{open_pos,..} | Frame::Group{open_pos,..} => *open_pos }; return Err(M6Error::UnbalancedBranchOpen { pos }); }
    let mut last_open: Option<usize> = None; for entry in ring_table.iter().flatten() { match last_open { None => last_open = Some(entry.open_pos), Some(p) => { if entry.open_pos > p { last_open = Some(entry.open_pos) } } } }
    if let Some(pos_open) = last_open { return Err(M6Error::RingUnclosed { open_pos: pos_open }); }
    let mut mols = builder.finish(); Ok(mols.pop().unwrap_or_default())
}