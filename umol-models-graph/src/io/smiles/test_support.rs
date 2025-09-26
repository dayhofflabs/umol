//! Test support: compact graph builder for expected Molecule IR

#[cfg(test)]
use crate::io::ir::builder::{BondData, MoleculeBuilder};
#[cfg(test)]
use crate::io::ir::{BondDir, BondOrder, Molecule};
#[cfg(test)]
use umol_data::Element;

#[cfg(test)]
fn parse_atom_token(tok: &str) -> (Element, bool) {
    // Asterisk denotes aromatic variant of the organic subset: C*, N*, O*, P*, S*, B*, c*, n*, ...
    // Two-letter atoms (Cl, Br) are recognized as-is; aromaticity via trailing '*'.
    let aromatic = tok.ends_with('*');
    let base = if aromatic { &tok[..tok.len() - 1] } else { tok };
    let el = match base {
        "C" | "c" => Element::C,
        "N" | "n" => Element::N,
        "O" | "o" => Element::O,
        "P" | "p" => Element::P,
        "S" | "s" => Element::S,
        "B" | "b" => Element::B,
        "F" | "f" => Element::F,
        "I" | "i" => Element::I,
        "Cl" | "cl" => Element::Cl,
        "Br" | "br" => Element::Br,
        other => panic!("unknown element token: {}", other),
    };
    (el, aromatic)
}

#[cfg(test)]
fn parse_bond_token(tok: &str) -> (usize, usize, BondOrder, Option<BondDir>) {
    // token forms: i-j, i-j:-, i-j:=, i-j:#, i-j:$, i-j::, i-j:/, i-j:\
    // indices are base-10, 0-based; direction optional only for single bonds
    let (ends, spec) = tok.split_once(':').unwrap_or((tok, "-"));
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
    (i, j, order, dir)
}

#[cfg(test)]
pub fn build_from_graph(spec: &str) -> Molecule {
    // Format: "atoms... | edges..."; atoms are space- or comma-separated tokens.
    // '*' suffix means aromatic atom.
    // edges are tokens like "i-j" or "i-j:=" etc., separated by spaces or commas.
    let (atoms_s, edges_s) = spec.split_once('|').expect("spec must have atoms | edges");
    let atoms: Vec<_> = atoms_s
        .split(|c: char| c == ' ' || c == ',')
        .filter(|t| !t.is_empty())
        .collect();
    let edges: Vec<_> = edges_s
        .split(|c: char| c == ' ' || c == ',')
        .filter(|t| !t.is_empty())
        .collect();

    let mut b = MoleculeBuilder::with_capacity(atoms.len(), edges.len());
    // map of insertion index to atom id is identity by construction; still collect ids
    let mut ids: Vec<u32> = Vec::with_capacity(atoms.len());
    for tok in atoms {
        let (el, arom) = parse_atom_token(tok);
        let id = b.on_atom_fast(el, true, arom);
        ids.push(id);
    }
    for etok in edges {
        let (i, j, order, dir) = parse_bond_token(etok);
        b.on_bond(ids[i], ids[j], BondData { order, dir });
    }
    let mut mols = b.finish();
    mols.pop().unwrap_or_default()
}


