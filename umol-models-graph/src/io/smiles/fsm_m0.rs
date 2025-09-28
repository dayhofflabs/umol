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
    use crate::io::smiles::test_support::build_from_graph;

    #[rstest]
    #[case::empty(b"", Molecule::default())]
    #[case::chain_c_1(b"C", build_from_graph("C |"))]
    #[case::chain_c_5(b"CCCCC", build_from_graph("C C C C C | 0-1 1-2 2-3 3-4"))]
    #[case::chain_mixed_5(b"CClOBrN", build_from_graph("C Cl O Br N | 0-1 1-2 2-3 3-4"))]
    fn m0_chain(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m0(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rstest]
    #[case::non_ascii(b"\xf0\x9f\x9c\x8d", M0Error::UnsupportedToken { pos: 0 })]
    #[case::comma(b",", M0Error::UnsupportedToken { pos: 0 })]
    #[case::semicolon(b";", M0Error::UnsupportedToken { pos: 0 })]
    #[case::question_mark(b"?", M0Error::UnsupportedToken { pos: 0 })]
    #[case::caret(b"^", M0Error::UnsupportedToken { pos: 0 })]
    #[case::pipe(b"|", M0Error::UnsupportedToken { pos: 0 })]
    #[case::open_angle_bracket(b"<", M0Error::UnsupportedToken { pos: 0 })]
    #[case::close_angle_bracket(b"<", M0Error::UnsupportedToken { pos: 0 })]
    #[case::open_brace(b"{", M0Error::UnsupportedToken { pos: 0 })]
    #[case::close_brace(b"}", M0Error::UnsupportedToken { pos: 0 })]
    #[case::single_quote(b"'", M0Error::UnsupportedToken { pos: 0 })]
    #[case::double_quote(b"\"", M0Error::UnsupportedToken { pos: 0 })]
    #[case::backtick(b"`", M0Error::UnsupportedToken { pos: 0 })]
    #[case::tilde(b"~", M0Error::UnsupportedToken { pos: 0 })]
    #[case::exclamation_mark(b"!", M0Error::UnsupportedToken { pos: 0 })]
    #[case::ampersand(b"&", M0Error::UnsupportedToken { pos: 0 })]
    #[case::underscore(b"_", M0Error::UnsupportedToken { pos: 0 })]
    #[case::bare_chirality(b"C@", M0Error::UnsupportedToken { pos: 1 })]
    #[case::bare_charge_pos(b"C+", M0Error::UnsupportedToken { pos: 1 })]
    #[case::bare_charge_neg(b"C-", M0Error::UnsupportedToken { pos: 1 })]
    #[case::bare_hcount(b"CH", M0Error::UnsupportedToken { pos: 1 })]
    #[case::bare_digit(b"1", M0Error::UnsupportedToken { pos: 0 })]
    fn m0_tokens_invalid(#[case] input: &[u8], #[case] expected: M0Error) {
        let res = parse_smiles_m0(input);
        assert!(res.is_err(), "{:?} should have failed", input);
        let err = res.unwrap_err();
        assert_eq!(err, expected);
    }

    #[rstest]
    #[case::invalid_element_1(b"X", M0Error::UnsupportedToken { pos: 0 })]
    #[case::invalid_element_2(b"Z", M0Error::UnsupportedToken { pos: 0 })]
    #[case::invalid_element_3(b"Aq", M0Error::UnsupportedToken { pos: 0 })]
    #[case::invalid_element_4(b"Sh", M0Error::UnsupportedToken { pos: 1 })]
    fn m0_chain_invalid(#[case] input: &[u8], #[case] expected: M0Error) {
        let res = parse_smiles_m0(input);
        assert!(res.is_err(), "{:?} should have failed", input);
        let err = res.unwrap_err();
        assert_eq!(err, expected);
    }

    #[rstest]
    #[case::wildcard(b"*", M0Error::UnsupportedToken { pos: 0 })]
    #[case::aromatic(b"c", M0Error::UnsupportedToken { pos: 0 })]
    #[case::bond_order(b"C-C", M0Error::UnsupportedToken { pos: 1 })]
    #[case::bracket(b"[C]", M0Error::UnsupportedToken { pos: 0 })]
    #[case::group(b"(C)", M0Error::UnsupportedToken { pos: 0 })]
    #[case::branch(b"CC(C)C", M0Error::UnsupportedToken { pos: 2 })]
    #[case::ring(b"C1CC1", M0Error::UnsupportedToken { pos: 1 })]
    #[case::ring_percent(b"C%12CC1%2", M0Error::UnsupportedToken { pos: 1 })]
    #[case::component(b"CC.CC", M0Error::UnsupportedToken { pos: 2 })]
    #[case::whitespace_1(b"C ", M0Error::UnsupportedToken { pos: 1 })]
    #[case::whitespace_2(b"C\t", M0Error::UnsupportedToken { pos: 1 })]
    #[case::whitespace_3(b"C\n", M0Error::UnsupportedToken { pos: 1 })]
    #[case::whitespace_4(b"C\r", M0Error::UnsupportedToken { pos: 1 })]
    #[case::whitespace_5(b"C\r\n", M0Error::UnsupportedToken { pos: 1 })]
    fn m0_unimplemented(#[case] input: &[u8], #[case] expected: M0Error) {
        let res = parse_smiles_m0(input);
        assert!(res.is_err(), "{:?} should have failed", input);
        let err = res.unwrap_err();
        assert_eq!(err, expected);
    }
}
