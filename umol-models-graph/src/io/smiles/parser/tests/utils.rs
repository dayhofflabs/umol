//! Utilities for SMILES parser tests.

use regex::Regex;
use umol_data::Element;

use crate::io::ir::builder::{AtomData, BondData, MoleculeBuilder};
use crate::io::ir::{BondDir, BondOrder, Molecule, Ring};

fn parse_atom_token(tok: &str) -> (Element, bool, Option<u32>) {
    // Asterisk denotes aromatic variant of the organic subset: C*, N*, O*, P*, S*, B*, ...
    // Two-letter atoms (Cl, Br) are recognized as-is; aromaticity via trailing '*'.
    // Optional '@<pos>' sets the atom span start.
    let pat = Regex::new(r"^([A-Z][a-z]?)(\*)?(?:@(\d+))?$").unwrap();
    let caps = pat.captures(tok).unwrap();
    let el = Element::from_symbol(caps.get(1).unwrap().as_str()).unwrap();
    let aromatic = caps.get(2).is_some();
    let start = caps
        .get(3)
        .map(|m| m.as_str().parse::<u32>().expect("valid u32 span"));
    (el, aromatic, start)
}

fn parse_bond_token(tok: &str) -> (usize, usize, BondOrder, Option<BondDir>, Option<u32>) {
    // Token forms:
    //   i-j
    //   i-j:<spec> where spec in -,=,#,$,:,/,\
    // Optional '@<pos>' suffix sets the bond span start.
    let (core, pos_opt) = tok
        .split_once('@')
        .map_or((tok, None), |(c, p)| (c, Some(p)));
    let span_start = pos_opt.map(|p| p.parse::<u32>().expect("valid u32 span"));
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
    (i, j, order, dir, span_start)
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
    for tok in atoms {
        let (el, arom, start) = parse_atom_token(tok);
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
        });
        ids.push(id);
    }
    for etok in edges {
        let (i, j, order, dir, span_start) = parse_bond_token(etok);
        b.on_bond(
            ids[i],
            ids[j],
            BondData {
                order,
                dir,
                span_start,
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
                let mut open_pos: Option<u32> = None;
                let mut close_pos: Option<u32> = None;
                if chars.peek() == Some(&'@') {
                    chars.next();
                    let mut num = String::new();
                    while let Some(c) = chars.peek() {
                        if c.is_ascii_digit() {
                            num.push(*c);
                            chars.next();
                        } else {
                            break;
                        }
                    }
                    if !num.is_empty() {
                        open_pos = Some(num.parse::<u32>().expect("open pos"));
                    }
                }
                if chars.peek() == Some(&'-') {
                    chars.next();
                    let mut num = String::new();
                    while let Some(c) = chars.peek() {
                        if c.is_ascii_digit() {
                            num.push(*c);
                            chars.next();
                        } else {
                            break;
                        }
                    }
                    if !num.is_empty() {
                        close_pos = Some(num.parse::<u32>().expect("close pos"));
                    }
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
                rings.push(Ring {
                    ring_idx,
                    open_pos,
                    close_pos,
                    atom_a,
                    atom_b,
                });
            }
            mol.ring_events = rings;
        }
    }

    mol
}
