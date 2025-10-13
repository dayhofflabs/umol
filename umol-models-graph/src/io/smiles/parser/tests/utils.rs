//! Utilities for SMILES parser tests.

use regex::Regex;
use umol_data::Element;

use crate::io::ir::builder::{AtomData, BondData, MoleculeBuilder};
use crate::io::ir::{BondDir, BondOrder, Molecule, Ring};

fn parse_atom_token(tok: &str) -> (Element, bool, Option<u32>, Option<u32>) {
    // Asterisk denotes aromatic variant of the organic subset: C*, N*, O*, P*, S*, B*, ...
    // Two-letter atoms (Cl, Br) are recognized as-is; aromaticity via trailing '*'.
    // Optional '@<start>' or '@<start>..<end>' sets span positions.
    let pat = Regex::new(r"^([A-Z][a-z]?)(\*)?(?:@(\d+)(?:\.\.(\d+))?)?$").unwrap();
    let caps = pat.captures(tok).unwrap();
    let el = Element::from_symbol(caps.get(1).unwrap().as_str()).unwrap();
    let aromatic = caps.get(2).is_some();
    let start = caps
        .get(3)
        .map(|m| m.as_str().parse::<u32>().expect("valid u32 span"));
    let end = caps
        .get(4)
        .map(|m| m.as_str().parse::<u32>().expect("valid u32 span"));
    (el, aromatic, start, end)
}

fn parse_bond_token(tok: &str) -> (usize, usize, BondOrder, Option<BondDir>, Option<u32>, Option<u32>) {
    // Token forms:
    //   i-j
    //   i-j:<spec> where spec in -,=,#,$,:,/,\
    // Optional '@<start>' or '@<start>..<end>' sets span positions.
    let (core, pos_opt) = tok
        .split_once('@')
        .map_or((tok, None), |(c, p)| (c, Some(p)));
    let (span_start, span_end) = pos_opt.map_or((None, None), |p| {
        if let Some((s, e)) = p.split_once("..") {
            let ss = s.parse::<u32>().expect("valid u32 span");
            let ee = e.parse::<u32>().expect("valid u32 span");
            (Some(ss), Some(ee))
        } else {
            (Some(p.parse::<u32>().expect("valid u32 span")), None)
        }
    });
    let (ends, spec) = core.split_once(':').unwrap_or((core, "-"));
    let (lhs, rhs) = ends.split_once('-').expect("edge must be i-j");
    let i = lhs.parse::<usize>().expect("left index");
    let j = rhs.parse::<usize>().expect("right index");
    let spec_norm = if spec.is_empty() { ":" } else { spec };
    let (order, dir) = match spec_norm {
        "-" => (BondOrder::Single, None),
        "=" => (BondOrder::Double, None),
        "#" => (BondOrder::Triple, None),
        "$" => (BondOrder::Quadruple, None),
        ":" => (BondOrder::Aromatic, None),
        "/" => (BondOrder::Single, Some(BondDir::Up)),
        "\\" => (BondOrder::Single, Some(BondDir::Down)),
        other => panic!("unknown bond spec: {}", other),
    };
    (i, j, order, dir, span_start, span_end)
}

pub fn build_from_graph(spec: &str) -> Molecule {
    // Format: "atoms... | edges... [| rings...]"
    // atoms: tokens like "C", "Cl", optional aromatic '*' and optional span "@<pos>": e.g. "C*@5"
    // edges: tokens like "i-j" or with type/dir: "i-j:=" "/" "\\" etc., optional span "@<pos>"
    // rings (optional third section): tokens encoding ring events with optional positions/atoms.
    // Grammar (examples):
    //   full:      idx@open-close:a-b   e.g. "1@2-7:0-5"
    //   open-only: idx@open:a           e.g. "3@10:2"
    //   close-only:idx-close:-b         e.g. "9-21:-4"
    let parts: Vec<_> = spec.split('|').map(|s| s.trim()).collect();
    assert!(parts.len() >= 2, "spec must have at least atoms | edges");
    let atoms_s = parts[0];
    let edges_s = parts[1];
    let rings_s = if parts.len() >= 3 {
        Some(parts[2])
    } else {
        None
    };

    let atoms: Vec<_> = atoms_s.split_whitespace().collect();
    let edges: Vec<_> = edges_s.split_whitespace().collect();

    let mut b = MoleculeBuilder::with_capacity(atoms.len(), edges.len());
    // map of insertion index to atom id is identity by construction; still collect ids
    let mut ids: Vec<u32> = Vec::with_capacity(atoms.len());
    // Keep a map from atom span_start -> span_end for bond defaulting
    use std::collections::HashMap;
    let mut atom_span_map: HashMap<u32, u32> = HashMap::new();

    for tok in atoms {
        let (el, arom, start, mut end) = parse_atom_token(tok);
        if end.is_none() {
            if let Some(s) = start {
                // Default width: aromatic tokens are 1 byte; aliphatic Cl/Br are 2; others 1
                let w: u32 = if arom { 1 } else {
                    match el {
                        Element::Cl | Element::Br => 2,
                        _ => 1,
                    }
                };
                end = Some(s + w);
            }
        }
        let id = b.on_atom(AtomData {
            element: el,
            aromatic: arom,
            implicit_h: true,
            isotope: None,
            charge: None,
            hydrogen_count: None,
            class: None,
            chirality: None,
            unknown_symbol: false,
            span_start: start,
            span_end: end,
        });
        if let (Some(s), Some(e)) = (start, end) { atom_span_map.insert(s, e); }
        ids.push(id);
    }
    for etok in edges {
        let (i, j, order, dir, span_start, mut span_end) = parse_bond_token(etok);
        if span_end.is_none() {
            if let Some(s) = span_start {
                span_end = atom_span_map.get(&s).copied().or(Some(s + 1));
            }
        }
        b.on_bond(
            ids[i],
            ids[j],
            BondData {
                order,
                dir,
                span_start,
                span_end,
            },
        );
    }
    let mut mols = b.finish();
    let mut mol = mols.pop().unwrap_or_default();

    if let Some(rings_src) = rings_s {
        let ring_tokens: Vec<_> = rings_src.split_whitespace().collect();
        if !ring_tokens.is_empty() {
            let mut rings: Vec<Ring> = Vec::with_capacity(ring_tokens.len());
            for rtok in ring_tokens {
                if rtok.is_empty() {
                    continue;
                }
                // Parse left/right around ':'
                let (left, right) = rtok.split_once(':').map_or((rtok, ""), |(l, r)| (l, r));
                // Parse left: idx[@open][-close]
                let mut chars = left.chars().peekable();
                let mut idx_str = String::new();
                while let Some(c) = chars.peek() {
                    if c.is_ascii_digit() {
                        idx_str.push(*c);
                        chars.next();
                    } else {
                        break;
                    }
                }
                let ring_idx: u32 = idx_str.parse().expect("ring idx");
                let mut open_start: Option<u32> = None;
                let mut close_pos: Option<u32> = None;
                let mut open_end: Option<u32> = None;
                let mut close_end: Option<u32> = None;
                if chars.peek() == Some(&'@') {
                    chars.next();
                    // parse start
                    let mut num = String::new();
                    while let Some(c) = chars.peek() {
                        if c.is_ascii_digit() { num.push(*c); chars.next(); } else { break; }
                    }
                    if !num.is_empty() { open_start = Some(num.parse::<u32>().expect("open pos")); }
                    // optional ..end
                    if chars.peek() == Some(&'.') { chars.next(); if chars.peek() == Some(&'.') { chars.next(); let mut num2 = String::new(); while let Some(c) = chars.peek() { if c.is_ascii_digit() { num2.push(*c); chars.next(); } else { break; } } if !num2.is_empty() { open_end = Some(num2.parse::<u32>().expect("open end")); } } }
                }
                if chars.peek() == Some(&'-') {
                    chars.next();
                    let mut num = String::new();
                    while let Some(c) = chars.peek() { if c.is_ascii_digit() { num.push(*c); chars.next(); } else { break; } }
                    if !num.is_empty() { close_pos = Some(num.parse::<u32>().expect("close pos")); }
                    if chars.peek() == Some(&'.') { chars.next(); if chars.peek() == Some(&'.') { chars.next(); let mut num2 = String::new(); while let Some(c) = chars.peek() { if c.is_ascii_digit() { num2.push(*c); chars.next(); } else { break; } } if !num2.is_empty() { close_end = Some(num2.parse::<u32>().expect("close end")); } } }
                }
                // Parse right: a-b (either side may be empty)
                let (a_opt, b_opt) = if right.is_empty() {
                    (None, None)
                } else {
                    let mut split = right.splitn(2, '-');
                    let a_str = split.next().unwrap_or("");
                    let b_str = split.next().unwrap_or("");
                    let a = if a_str.is_empty() {
                        None
                    } else {
                        Some(a_str.parse::<usize>().expect("a"))
                    };
                    let b = if b_str.is_empty() {
                        None
                    } else {
                        Some(b_str.parse::<usize>().expect("b"))
                    };
                    (a, b)
                };
                let atom_a = a_opt.map(|ai| ids[ai]);
                let atom_b = b_opt.map(|bi| ids[bi]);
                // Default ends if positions are provided without ends
                if open_end.is_none() { if let Some(s) = open_start { open_end = Some(s + 1); } }
                if close_end.is_none() { if let Some(s) = close_pos { close_end = Some(s + 1); } }
                rings.push(Ring {
                    ring_idx,
                    open_start,
                    close_start: close_pos,
                    atom_a,
                    atom_b,
                    open_end,
                    close_end,
                });
            }
            mol.ring_events = rings;
        }
    }

    mol
}
