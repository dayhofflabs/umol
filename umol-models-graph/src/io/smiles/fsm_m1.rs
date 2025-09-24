use crate::io::ir::Molecule;
use crate::io::ir::builder::MoleculeBuilder;
use umol_data::Element;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum M1Error {
    UnsupportedToken { pos: usize },
    UnbalancedBranchOpen,
    UnbalancedBranchClose { pos: usize },
}

// M1: chains + branches (no rings, no aromatic, no charged/bracketed atoms)
pub fn parse_smiles_m1(input: &[u8]) -> Result<Molecule, M1Error> {
    let mut i = 0usize;
    let n = input.len();

    let mut builder = MoleculeBuilder::with_capacity(n.max(1), n.max(1).saturating_sub(1));
    let mut last_atom_idx: Option<u32> = None;

    // Stack of branch attach points (indices of base atom)
    let mut branch_stack: Vec<u32> = Vec::new();

    while i < n {
        let b0 = input[i];

        // Branch start: '(' — requires an attach point
        if b0 == b'(' {
            let base = match last_atom_idx {
                Some(idx) => idx,
                None => return Err(M1Error::UnsupportedToken { pos: i }),
            };
            branch_stack.push(base);
            i += 1;
            continue;
        }

        // Branch end: ')' — pop attach point and restore as current
        if b0 == b')' {
            let Some(base) = branch_stack.pop() else { return Err(M1Error::UnbalancedBranchClose { pos: i }); };
            last_atom_idx = Some(base);
            i += 1;
            continue;
        }

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

        return Err(M1Error::UnsupportedToken { pos: i });
    }

    if !branch_stack.is_empty() { return Err(M1Error::UnbalancedBranchOpen); }

    let mut mols = builder.finish();
    Ok(mols.pop().unwrap_or_default())
}


