//! Utilities for SMILES parser tests.

#![allow(clippy::type_complexity)]

use std::collections::HashMap;

use regex::Regex;
use umol_shared::element::Element;

use super::super::builder::{
    AtomData, BondData, ExtendedAtomData, ExtendedMoleculeBuilder, MoleculeBuilder,
};
use crate::span::Span;
use crate::table_ir::atom::Chirality;
use crate::table_ir::{
    AtomSymbol, BondDonation, BondOrder, BondWedge, ChiralityFrame, ExtendedMolecule, Molecule,
    SourceFormat, WildcardAtom,
};

/// Returns the sorted list of neighbor atom indices for a given atom in a Molecule.
pub fn get_atom_neighbors(mol: &Molecule, atom_idx: u32) -> Vec<u32> {
    let mut result = Vec::new();
    for bond in &mol.bonds {
        let (a, b) = bond.atoms.as_tuple();
        if a == atom_idx {
            result.push(b);
        } else if b == atom_idx {
            result.push(a);
        }
    }
    result.sort();
    result
}

/// Finds the first chiral atom in a Molecule.
/// Returns (atom_index, element, chirality, sorted_neighbors) or None if no chiral atom found.
pub fn find_chiral_center(mol: &Molecule) -> Option<(usize, Element, Chirality, Vec<u32>)> {
    for (idx, atom) in mol.atoms.iter().enumerate() {
        if let Some(chir) = atom.chirality {
            let neighbors = get_atom_neighbors(mol, idx as u32);
            return Some((idx, atom.element, chir, neighbors));
        }
    }
    None
}

/// Returns the sorted list of neighbor atom indices for a given atom in an ExtendedMolecule.
pub fn get_extended_atom_neighbors(mol: &ExtendedMolecule, atom_idx: u32) -> Vec<u32> {
    let mut result = Vec::new();
    for bond in &mol.bonds {
        let (a, b) = bond.atoms.as_tuple();
        if a == atom_idx {
            result.push(b);
        } else if b == atom_idx {
            result.push(a);
        }
    }
    result.sort();
    result
}

/// Finds the first chiral atom in an ExtendedMolecule.
/// Returns (atom_index, element, chirality, sorted_neighbors) or None if no chiral atom found.
pub fn find_extended_chiral_center(
    mol: &ExtendedMolecule,
) -> Option<(usize, Element, Chirality, Vec<u32>)> {
    for (idx, atom) in mol.atoms.iter().enumerate() {
        if let Some(chir) = atom.chirality {
            if let AtomSymbol::Element(el) = atom.symbol {
                let neighbors = get_extended_atom_neighbors(mol, idx as u32);
                return Some((idx, el, chir, neighbors));
            }
        }
    }
    None
}

/// Finds the first bond with stereo direction in a Molecule.
/// Returns (atom1, atom2, wedge) or None if no stereo bond found.
pub fn find_stereo_bond(mol: &Molecule) -> Option<(u32, u32, BondWedge)> {
    for bond in &mol.bonds {
        if let Some(wedge) = bond.wedge {
            let (a, b) = bond.atoms.as_tuple();
            return Some((a, b, wedge));
        }
    }
    None
}

/// Finds the first bond with stereo wedge in an ExtendedMolecule.
/// Returns (atom1, atom2, wedge) or None if no stereo bond found.
pub fn find_extended_stereo_bond(mol: &ExtendedMolecule) -> Option<(u32, u32, BondWedge)> {
    for bond in &mol.bonds {
        if let Some(wedge) = bond.wedge {
            let (a, b) = bond.atoms.as_tuple();
            return Some((a, b, wedge));
        }
    }
    None
}

fn parse_atom_token(tok: &str) -> (Element, bool, Option<u32>, Option<u32>) {
    // Underscore denotes aromatic variant of the organic subset: C_, N_, O_, P_, S_, B_, ...
    // Optional '@<start>' or '@<start>..<end>' sets span positions.
    let pat = Regex::new(r"^([A-Z][a-z]?)(_)?(?:@(\d+)(?:\.\.(\d+))?)?$").unwrap();
    let caps = pat.captures(tok).expect("valid atom token");
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

enum ExtendedAtomSymbol {
    Element(Element),
    Wildcard,
}

fn parse_extended_atom_token(tok: &str) -> (ExtendedAtomSymbol, bool, Option<u32>, Option<u32>) {
    // Supports element symbols (C, Cl, etc.), aromaticity via '_', and wildcard '*'.
    // Optional '@<start>' or '@<start>..<end>' sets span positions.
    let pat = Regex::new(r"^(\*|[A-Z][a-z]?)(_)?(?:@(\d+)(?:\.\.(\d+))?)?$").unwrap();
    let caps = pat.captures(tok).expect("valid extended atom token");
    let sym_str = caps.get(1).unwrap().as_str();
    let symbol = if sym_str == "*" {
        ExtendedAtomSymbol::Wildcard
    } else {
        ExtendedAtomSymbol::Element(Element::from_symbol(sym_str).unwrap())
    };
    let aromatic = caps.get(2).is_some();
    let start = caps
        .get(3)
        .map(|m| m.as_str().parse::<u32>().expect("valid u32 span"));
    let end = caps
        .get(4)
        .map(|m| m.as_str().parse::<u32>().expect("valid u32 span"));
    (symbol, aromatic, start, end)
}

fn parse_bond_token(
    tok: &str,
) -> (
    usize,
    usize,
    BondOrder,
    Option<BondWedge>,
    Option<BondDonation>,
    Option<u32>,
    Option<u32>,
) {
    // Token forms:
    //   i-j                     (single bond)
    //   i-j:<spec>              (with explicit spec)
    //   i-j<spec>               (spec without colon for ~, ->, <-)
    // Optional '@<start>' or '@<start>..<end>' sets span positions.
    // Use regex to extract: <idx>-<idx><rest>
    let re = Regex::new(r"^(\d+)-(\d+)(.*)$").unwrap();
    let caps = re.captures(tok).expect("edge must match i-j pattern");
    let i: usize = caps[1].parse().expect("left index");
    let j: usize = caps[2].parse().expect("right index");
    let rest = &caps[3];

    // Parse rest: optional ':'<spec> or just <spec>, then optional @<span>
    let (spec_part, span_part) = if let Some(at_pos) = rest.find('@') {
        (&rest[..at_pos], Some(&rest[at_pos + 1..]))
    } else {
        (rest, None)
    };

    let (span_start, span_end) = span_part.map_or((None, None), |p| {
        if let Some((s, e)) = p.split_once("..") {
            let ss = s.parse::<u32>().expect("valid u32 span");
            let ee = e.parse::<u32>().expect("valid u32 span");
            (Some(ss), Some(ee))
        } else {
            (Some(p.parse::<u32>().expect("valid u32 span")), None)
        }
    });

    // Parse spec: may have leading ':' as separator, or ':' alone means aromatic
    let spec_norm = if spec_part.is_empty() {
        "-" // default to single bond
    } else if spec_part == ":" {
        ":" // just ':' means aromatic
    } else if let Some(s) = spec_part.strip_prefix(':') {
        if s.is_empty() {
            ":"
        } else {
            s
        }
    } else {
        spec_part // no colon prefix
    };

    let (order, wedge, donation) = match spec_norm {
        "-" => (BondOrder::Single, None, None),
        "=" => (BondOrder::Double, None, None),
        "#" => (BondOrder::Triple, None, None),
        "$" => (BondOrder::Quadruple, None, None),
        ":" => (BondOrder::Aromatic, None, None),
        "/" => (BondOrder::Single, Some(BondWedge::Up), None),
        "\\" => (BondOrder::Single, Some(BondWedge::Down), None),
        "~" => (BondOrder::Any, None, None),
        "->" => (BondOrder::Single, None, Some(BondDonation::Donating)),
        "<-" => (BondOrder::Single, None, Some(BondDonation::Accepting)),
        other => panic!("unknown bond spec: {}", other),
    };
    (i, j, order, wedge, donation, span_start, span_end)
}

pub fn build_from_graph(spec: &str) -> Molecule {
    // Format: "atoms... | edges... [| rings...]"
    // atoms: tokens like "C", "Cl", optional aromatic '_' and optional span "@<pos>": e.g. "C_@5"
    // edges: tokens like "i-j" or with type/wedge: "i-j:=" "/" "\\" etc., optional span "@<pos>"
    // rings (optional third section): tokens encoding ring events with optional positions/atoms.
    // Grammar (examples):
    //   full:      idx@open-close:a-b   e.g. "1@2-7:0-5"
    //   open-only: idx@open:a           e.g. "3@10:2"
    //   close-only:idx-close:-b         e.g. "9-21:-4"
    let parts: Vec<_> = spec.split('|').map(|s| s.trim()).collect();
    assert!(parts.len() >= 2, "spec must have at least atoms | edges");
    let atoms_s = parts[0];
    let edges_s = parts[1];

    let atoms: Vec<_> = atoms_s.split_whitespace().collect();
    let edges: Vec<_> = edges_s.split_whitespace().collect();

    let mut b = MoleculeBuilder::with_capacity(atoms.len(), edges.len());
    // map of insertion index to atom id is identity by construction; still collect ids
    let mut ids: Vec<u32> = Vec::with_capacity(atoms.len());
    // Keep a map from atom span_start -> span_end for bond defaulting
    let mut atom_span_map: HashMap<u32, u32> = HashMap::new();

    for tok in atoms {
        let (el, arom, start, mut end) = parse_atom_token(tok);
        if end.is_none() {
            if let Some(s) = start {
                // Default width: aromatic tokens are 1 byte; aliphatic Cl/Br are 2; others 1
                let w: u32 = if arom {
                    1
                } else {
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
            isotope: None,
            charge: None,
            implicit_hydrogens: None,
            class: None,
            chirality: None,
            span: Span::from_bytes_opt(start, end),
        });
        if let (Some(s), Some(e)) = (start, end) {
            atom_span_map.insert(s, e);
        }
        ids.push(id);
    }
    for etok in edges {
        let (i, j, order, wedge, donation, span_start, mut span_end) = parse_bond_token(etok);
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
                wedge,
                donation,
                span: Span::from_bytes_opt(span_start, span_end),
            },
        );
    }
    let mut mols = b.finish();
    let mut mol = mols.pop().unwrap_or_else(Molecule::empty);

    mol.source_format = SourceFormat::SMILES;
    mol.chirality_frame = Some(ChiralityFrame::FirstNeighborToward);
    mol
}

pub fn build_extended_from_graph(spec: &str) -> ExtendedMolecule {
    // Format: "atoms... | edges... [| rings...]"
    // atoms: tokens like "C", "Cl", "*" (wildcard), optional aromatic '_' and span "@<pos>"
    // edges: tokens like "i-j" or with type/wedge: "i-j:=" "/" "\\" etc., optional span "@<pos>"
    // rings (optional third section): tokens encoding ring events with optional positions/atoms.
    let parts: Vec<_> = spec.split('|').map(|s| s.trim()).collect();
    assert!(parts.len() >= 2, "spec must have at least atoms | edges");
    let atoms_s = parts[0];
    let edges_s = parts[1];

    let atoms: Vec<_> = atoms_s.split_whitespace().collect();
    let edges: Vec<_> = edges_s.split_whitespace().collect();

    let mut b = ExtendedMoleculeBuilder::with_capacity(atoms.len(), edges.len());
    let mut ids: Vec<u32> = Vec::with_capacity(atoms.len());
    let mut atom_span_map: HashMap<u32, u32> = HashMap::new();

    for tok in atoms {
        let (sym, arom, start, mut end) = parse_extended_atom_token(tok);
        if end.is_none() {
            if let Some(s) = start {
                let w: u32 = match &sym {
                    ExtendedAtomSymbol::Wildcard => 1,
                    ExtendedAtomSymbol::Element(_) if arom => 1,
                    ExtendedAtomSymbol::Element(el) => match el {
                        Element::Cl | Element::Br => 2,
                        _ => 1,
                    },
                };
                end = Some(s + w);
            }
        }
        let id = match sym {
            ExtendedAtomSymbol::Wildcard => b.on_wildcard(WildcardAtom::Any, None, start, end),
            ExtendedAtomSymbol::Element(el) => b.on_atom(ExtendedAtomData {
                symbol: AtomSymbol::Element(el),
                aromatic: arom,
                isotope: None,
                charge: None,
                implicit_hydrogens: None,
                class: None,
                chirality: None,
                span: Span::from_bytes_opt(start, end),
            }),
        };
        if let (Some(s), Some(e)) = (start, end) {
            atom_span_map.insert(s, e);
        }
        ids.push(id);
    }
    for etok in edges {
        let (i, j, order, wedge, donation, span_start, mut span_end) = parse_bond_token(etok);
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
                wedge,
                donation,
                span: Span::from_bytes_opt(span_start, span_end),
            },
        );
    }
    let mut mols = b.finish();
    let mut mol = mols.pop().unwrap_or_else(ExtendedMolecule::empty);

    mol.source_format = SourceFormat::SMILES;
    mol.chirality_frame = Some(ChiralityFrame::FirstNeighborToward);
    mol
}
