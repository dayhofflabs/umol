use crate::io::ir::Molecule;
use crate::io::ir::builder::MoleculeBuilder;
use umol_data::Element;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum M0Error {
    UnsupportedToken { pos: usize },
}

pub fn parse_smiles_m0(input: &[u8]) -> Result<Molecule, M0Error> {
    let mut i = 0usize;
    let n = input.len();

    let mut builder = MoleculeBuilder::with_capacity(n.max(1), n.max(1).saturating_sub(1));
    let mut last_atom_idx: Option<u32> = None;

    while i < n {
        let b0 = input[i];

        // Recognize two-letter halogens first: Cl, Br
        if b0 == b'C' {
            if i + 1 < n && input[i + 1] == b'l' {
                let curr = builder.on_atom_fast(Element::Cl, true, false);
                if let Some(last) = last_atom_idx { builder.on_bond_single_fast(last, curr); }
                last_atom_idx = Some(curr);
                i += 2;
                continue;
            }
            // Single C
            let curr = builder.on_atom_fast(Element::C, true, false);
            if let Some(last) = last_atom_idx { builder.on_bond_single_fast(last, curr); }
            last_atom_idx = Some(curr);
            i += 1;
            continue;
        }
        if b0 == b'B' {
            if i + 1 < n && input[i + 1] == b'r' {
                let curr = builder.on_atom_fast(Element::Br, true, false);
                if let Some(last) = last_atom_idx { builder.on_bond_single_fast(last, curr); }
                last_atom_idx = Some(curr);
                i += 2;
                continue;
            }
            // Single B
            let curr = builder.on_atom_fast(Element::B, true, false);
            if let Some(last) = last_atom_idx { builder.on_bond_single_fast(last, curr); }
            last_atom_idx = Some(curr);
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
            if let Some(last) = last_atom_idx { builder.on_bond_single_fast(last, curr); }
            last_atom_idx = Some(curr);
            i += 1;
            continue;
        }

        return Err(M0Error::UnsupportedToken { pos: i });
    }

    let mut mols = builder.finish();
    Ok(mols.pop().unwrap_or_default())
}


