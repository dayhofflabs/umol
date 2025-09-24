use umol_data::Element;

use crate::io::ir::builder::MoleculeBuilder;
use crate::io::ir::Molecule;

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
                if let Some(last) = last_atom_idx {
                    builder.on_bond_single_fast(last, curr);
                }
                last_atom_idx = Some(curr);
                i += 2;
                continue;
            }
            // Single C
            let curr = builder.on_atom_fast(Element::C, true, false);
            if let Some(last) = last_atom_idx {
                builder.on_bond_single_fast(last, curr);
            }
            last_atom_idx = Some(curr);
            i += 1;
            continue;
        }
        if b0 == b'B' {
            if i + 1 < n && input[i + 1] == b'r' {
                let curr = builder.on_atom_fast(Element::Br, true, false);
                if let Some(last) = last_atom_idx {
                    builder.on_bond_single_fast(last, curr);
                }
                last_atom_idx = Some(curr);
                i += 2;
                continue;
            }
            // Single B
            let curr = builder.on_atom_fast(Element::B, true, false);
            if let Some(last) = last_atom_idx {
                builder.on_bond_single_fast(last, curr);
            }
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
            if let Some(last) = last_atom_idx {
                builder.on_bond_single_fast(last, curr);
            }
            last_atom_idx = Some(curr);
            i += 1;
            continue;
        }

        return Err(M0Error::UnsupportedToken { pos: i });
    }

    let mut mols = builder.finish();
    Ok(mols.pop().unwrap_or_default())
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::empty(b"", (0, 0))]
    #[case::c_chain_1_atom(b"C", (1, 0))]
    #[case::c_chain_5_atoms(b"CCCCC", (5, 4))]
    #[case::mixed_chain_5_atoms(b"CClCBrC", (5, 4))]
    fn m0_chain(#[case] input: &[u8], #[case] (exp_atoms, exp_bonds): (usize, usize)) {
        let res = parse_smiles_m0(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol.atoms.len(), exp_atoms);
        assert_eq!(mol.bonds.len(), exp_bonds);
    }

    #[rstest]
    #[case::bond_order(b"C-C", M0Error::UnsupportedToken { pos: 1 })]
    #[case::bracket(b"[C]", M0Error::UnsupportedToken { pos: 0 })]
    #[case::group(b"(C)", M0Error::UnsupportedToken { pos: 0 })]
    #[case::branch(b"CC(C)C", M0Error::UnsupportedToken { pos: 2 })]
    #[case::ring(b"C1CC1", M0Error::UnsupportedToken { pos: 1 })]
    #[case::component(b"CC.CC", M0Error::UnsupportedToken { pos: 2 })]
    fn m0_chain_invalid(#[case] input: &[u8], #[case] expected: M0Error) {
        let res = parse_smiles_m0(input);
        assert!(res.is_err(), "{:?} should have failed", input);
        let err = res.unwrap_err();
        assert_eq!(err, expected);
    }
}
