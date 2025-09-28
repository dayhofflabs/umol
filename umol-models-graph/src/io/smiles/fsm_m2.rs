use umol_data::Element;

use crate::io::ir::builder::{BondData, MoleculeBuilder};
use crate::io::ir::{BondDir, BondOrder, Molecule};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum M2Error {
    UnsupportedToken { pos: usize },
    UnbalancedBranchOpen { pos: usize },
    UnbalancedBranchClose { pos: usize },
    EmptyBranch { pos: usize },
    EmptyGroup { pos: usize },
    TopLevelGroupTrailing { pos: usize },
    TrailingBond { pos: usize },
    ConsecutiveBond { pos: usize },
    LeadingBond { pos: usize },
}

#[derive(Debug, Clone, Copy)]
enum Frame {
    Branch {
        base: u32,
        had_atom: bool,
        open_pos: usize,
    },
    Group {
        had_atom: bool,
        open_pos: usize,
    },
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

// M2 (bonds): chains + groups/branches + explicit bond tokens; no rings
pub fn parse_smiles_m2(input: &[u8]) -> Result<Molecule, M2Error> {
    let mut i = 0usize;
    let n = input.len();

    let mut builder = MoleculeBuilder::with_capacity(n.max(1), n.max(1).saturating_sub(1));
    let mut last_atom_idx: Option<u32> = None;
    let mut pending_bond: Option<(BondOrder, Option<BondDir>, usize)> = None;

    let mut pstack: Vec<Frame> = Vec::new();

    while i < n {
        let b0 = input[i];

        // Parentheses
        if b0 == b'(' {
            if let Some((_, _, pos)) = pending_bond {
                return Err(M2Error::TrailingBond { pos });
            }
            match last_atom_idx {
                Some(idx) => pstack.push(Frame::Branch {
                    base: idx,
                    had_atom: false,
                    open_pos: i,
                }),
                None => pstack.push(Frame::Group {
                    had_atom: false,
                    open_pos: i,
                }),
            }
            i += 1;
            continue;
        }
        if b0 == b')' {
            if let Some((_, _, pos)) = pending_bond {
                return Err(M2Error::TrailingBond { pos });
            }
            let Some(frame) = pstack.pop() else {
                return Err(M2Error::UnbalancedBranchClose { pos: i });
            };
            match frame {
                Frame::Branch { base, had_atom, .. } => {
                    if !had_atom {
                        return Err(M2Error::EmptyBranch { pos: i });
                    }
                    last_atom_idx = Some(base);
                }
                Frame::Group { had_atom, .. } => {
                    if !had_atom {
                        if i + 1 != n {
                            return Err(M2Error::EmptyGroup { pos: i });
                        }
                    } else {
                        // Disallow anything following a non-empty TOP-LEVEL group
                        if pstack.is_empty() && i + 1 != n {
                            return Err(M2Error::TopLevelGroupTrailing { pos: i });
                        }
                    }
                }
            }
            i += 1;
            continue;
        }

        // Bond tokens
        if matches!(b0, b'-' | b'=' | b'#' | b'$' | b':' | b'/' | b'\\') {
            if pending_bond.is_some() {
                return Err(M2Error::ConsecutiveBond { pos: i });
            }
            if last_atom_idx.is_none() {
                return Err(M2Error::LeadingBond { pos: i });
            }
            let (order, dir) = map_bond(b0);
            pending_bond = Some((order, dir, i));
            i += 1;
            continue;
        }

        // Two-letter halogens first: Cl, Br
        if b0 == b'C' {
            if i + 1 < n && input[i + 1] == b'l' {
                let curr = builder.on_atom_fast(Element::Cl, true, false);
                if let Some(last) = last_atom_idx {
                    if let Some((order, dir, _pos)) = pending_bond.take() {
                        builder.on_bond(last, curr, BondData { order, dir });
                    } else {
                        builder.on_bond_single_fast(last, curr);
                    }
                }
                last_atom_idx = Some(curr);
                if let Some(top) = pstack.last_mut() {
                    match top {
                        Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                            *had_atom = true
                        }
                    }
                }
                i += 2;
                continue;
            }
            // Single C
            let curr = builder.on_atom_fast(Element::C, true, false);
            if let Some(last) = last_atom_idx {
                if let Some((order, dir, _pos)) = pending_bond.take() {
                    builder.on_bond(last, curr, BondData { order, dir });
                } else {
                    builder.on_bond_single_fast(last, curr);
                }
            }
            last_atom_idx = Some(curr);
            if let Some(top) = pstack.last_mut() {
                match top {
                    Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                        *had_atom = true
                    }
                }
            }
            i += 1;
            continue;
        }
        if b0 == b'B' {
            if i + 1 < n && input[i + 1] == b'r' {
                let curr = builder.on_atom_fast(Element::Br, true, false);
                if let Some(last) = last_atom_idx {
                    if let Some((order, dir, _pos)) = pending_bond.take() {
                        builder.on_bond(last, curr, BondData { order, dir });
                    } else {
                        builder.on_bond_single_fast(last, curr);
                    }
                }
                last_atom_idx = Some(curr);
                if let Some(top) = pstack.last_mut() {
                    match top {
                        Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                            *had_atom = true
                        }
                    }
                }
                i += 2;
                continue;
            }
            // Single B
            let curr = builder.on_atom_fast(Element::B, true, false);
            if let Some(last) = last_atom_idx {
                if let Some((order, dir, _pos)) = pending_bond.take() {
                    builder.on_bond(last, curr, BondData { order, dir });
                } else {
                    builder.on_bond_single_fast(last, curr);
                }
            }
            last_atom_idx = Some(curr);
            if let Some(top) = pstack.last_mut() {
                match top {
                    Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                        *had_atom = true
                    }
                }
            }
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
                if let Some((order, dir, _pos)) = pending_bond.take() {
                    builder.on_bond(last, curr, BondData { order, dir });
                } else {
                    builder.on_bond_single_fast(last, curr);
                }
            }
            last_atom_idx = Some(curr);
            if let Some(top) = pstack.last_mut() {
                match top {
                    Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                        *had_atom = true
                    }
                }
            }
            i += 1;
            continue;
        }

        return Err(M2Error::UnsupportedToken { pos: i });
    }

    if pending_bond.is_some() {
        let (_, _, pos) = pending_bond.unwrap();
        return Err(M2Error::TrailingBond { pos });
    }

    if !pstack.is_empty() {
        let pos = match pstack.last().unwrap() {
            Frame::Branch { open_pos, .. } | Frame::Group { open_pos, .. } => *open_pos,
        };
        return Err(M2Error::UnbalancedBranchOpen { pos });
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
    fn m2_chain(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m2(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rstest]
    #[case::non_ascii(b"\xf0\x9f\x9c\x8d", M2Error::UnsupportedToken { pos: 0 })]
    #[case::comma(b",", M2Error::UnsupportedToken { pos: 0 })]
    #[case::semicolon(b";", M2Error::UnsupportedToken { pos: 0 })]
    #[case::question_mark(b"?", M2Error::UnsupportedToken { pos: 0 })]
    #[case::caret(b"^", M2Error::UnsupportedToken { pos: 0 })]
    #[case::pipe(b"|", M2Error::UnsupportedToken { pos: 0 })]
    #[case::open_angle_bracket(b"<", M2Error::UnsupportedToken { pos: 0 })]
    #[case::close_angle_bracket(b"<", M2Error::UnsupportedToken { pos: 0 })]
    #[case::open_brace(b"{", M2Error::UnsupportedToken { pos: 0 })]
    #[case::close_brace(b"}", M2Error::UnsupportedToken { pos: 0 })]
    #[case::single_quote(b"'", M2Error::UnsupportedToken { pos: 0 })]
    #[case::double_quote(b"\"", M2Error::UnsupportedToken { pos: 0 })]
    #[case::backtick(b"`", M2Error::UnsupportedToken { pos: 0 })]
    #[case::tilde(b"~", M2Error::UnsupportedToken { pos: 0 })]
    #[case::exclamation_mark(b"!", M2Error::UnsupportedToken { pos: 0 })]
    #[case::ampersand(b"&", M2Error::UnsupportedToken { pos: 0 })]
    #[case::underscore(b"_", M2Error::UnsupportedToken { pos: 0 })]
    #[case::bare_chirality(b"C@", M2Error::UnsupportedToken { pos: 1 })]
    #[case::bare_charge_pos(b"C+", M2Error::UnsupportedToken { pos: 1 })]
    #[case::bare_charge_neg(b"C-", M2Error::TrailingBond { pos: 1 })]
    #[case::bare_hcount(b"CH", M2Error::UnsupportedToken { pos: 1 })]
    #[case::bare_digit(b"1", M2Error::UnsupportedToken { pos: 0 })]
    fn m2_tokens_invalid(#[case] input: &[u8], #[case] expected: M2Error) {
        let res = parse_smiles_m2(input);
        assert!(res.is_err(), "{:?} should have failed", input);
        let err = res.unwrap_err();
        assert_eq!(err, expected);
    }

    #[rstest]
    #[case::invalid_element_1(b"X", M2Error::UnsupportedToken { pos: 0 })]
    #[case::invalid_element_2(b"Z", M2Error::UnsupportedToken { pos: 0 })]
    #[case::invalid_element_3(b"Aq", M2Error::UnsupportedToken { pos: 0 })]
    #[case::invalid_element_4(b"Sh", M2Error::UnsupportedToken { pos: 1 })]
    fn m2_chain_invalid(#[case] input: &[u8], #[case] expected: M2Error) {
        let res = parse_smiles_m2(input);
        assert!(res.is_err(), "{:?} should have failed", input);
        let err = res.unwrap_err();
        assert_eq!(err, expected);
    }

    #[rstest]
    #[case::empty_group(b"()", Molecule::default())]
    #[case::group_c_1(b"(C)", build_from_graph("C |"))]
    #[case::group_c_4(b"(CCCC)", build_from_graph("C C C C | 0-1 1-2 2-3"))]
    #[case::group_nested(b"((CC))", build_from_graph("C C | 0-1"))]
    #[case::branch_c_111(b"C(C)(C)", build_from_graph("C C C | 0-1 0-2"))]
    #[case::branch_c_211(b"CC(C)C", build_from_graph("C C C C | 0-1 1-2 1-3"))]
    #[case::trailing_branch(b"C(CC)", build_from_graph("C C C | 0-1 1-2"))]
    #[case::group_branched_c_3(b"C(C)(C)", build_from_graph("C C C | 0-1 0-2"))]
    fn m2_tree(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m2(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rstest]
    #[case::unbalanced_closing_paren_1(b")C", M2Error::UnbalancedBranchClose { pos: 0 })]
    #[case::unbalanced_closing_paren_2(b"C)C", M2Error::UnbalancedBranchClose { pos: 1 })]
    #[case::unclosed_group(b"(C", M2Error::UnbalancedBranchOpen { pos: 0 })]
    #[case::unclosed_branch(b"C(C", M2Error::UnbalancedBranchOpen { pos: 1 })]
    #[case::empty_branch(b"C()", M2Error::EmptyBranch { pos: 2 })]
    #[case::empty_group_before_atom(b"()C", M2Error::EmptyGroup { pos: 1 })]
    #[case::group_before_atom(b"(C)C", M2Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::two_top_level_groups(b"(C)(C)", M2Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::three_top_level_groups(b"(C)(C)(C)", M2Error::TopLevelGroupTrailing { pos: 2 })]
    fn m2_tree_invalid(#[case] input: &[u8], #[case] expected: M2Error) {
        let err = parse_smiles_m2(input);
        assert!(err.is_err(), "{:?} should have failed", input);
        let err = err.unwrap_err();
        assert_eq!(err, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single_bond(b"C-C", build_from_graph("C C | 0-1:-"))]
    #[case::double_bond(b"C=C", build_from_graph("C C | 0-1:="))]
    #[case::triple_bond(b"C#C", build_from_graph("C C | 0-1:#"))]
    #[case::quadruple_bond(b"C$C", build_from_graph("C C | 0-1:$"))]
    #[case::aromatic_bond(b"C:C", build_from_graph("C C | 0-1::"))]
    #[case::up_bond(b"C/C", build_from_graph("C C | 0-1:/"))]
    #[case::down_bond(b"C\\C", build_from_graph("C C | 0-1:\\"))]
    #[case::cumulated_bonds(b"C=C=C", build_from_graph("C C C | 0-1:= 1-2:="))]
    #[case::conjugated_bonds(b"C=CC=C", build_from_graph("C C C C | 0-1:= 1-2:- 2-3:="))]
    #[case::branch_leading_bond(b"CC(-C)C", build_from_graph("C C C C | 0-1 1-2 1-3"))]
    #[case::branch_leading_double_bond(b"CC(=C)C", build_from_graph("C C C C | 0-1 1-2:= 1-3"))]
    #[case::branch_internal_bond(b"CC(C-C)C", build_from_graph("C C C C C | 0-1 1-2 2-3 1-4"))]
    #[case::branch_internal_double_bond(b"CC(C=C)C", build_from_graph("C C C C C | 0-1 1-2 2-3:= 1-4"))]
    #[case::branch_followed_by_bond(b"CC(C)-C", build_from_graph("C C C C | 0-1 1-2 1-3"))]
    #[case::branch_followed_by_double_bond(b"CC(C)=C", build_from_graph("C C C C | 0-1 1-2 1-3:="))]
    #[case::branch_trans_double_bond_1(b"C/C=C/C", build_from_graph("C C C C | 0-1:/ 1-2:= 2-3:/"))]
    #[case::branch_trans_double_bond_2(b"C\\C=C\\C", build_from_graph("C C C C | 0-1:\\ 1-2:= 2-3:\\"))]
    #[case::branch_cis_double_bond_1(b"C\\C=C/C", build_from_graph("C C C C | 0-1:\\ 1-2:= 2-3:/"))]
    #[case::branch_cis_double_bond_2(b"C/C=C\\C", build_from_graph("C C C C | 0-1:/ 1-2:= 2-3:\\"))]
    fn m2_bonds(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m2(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rstest]
    #[case::trailing_bond_1(b"C-", M2Error::TrailingBond { pos: 1 })]
    #[case::trailing_bond_2(b"C=", M2Error::TrailingBond { pos: 1 })]
    #[case::trailing_bond_3(b"C#", M2Error::TrailingBond { pos: 1 })]
    #[case::trailing_bond_4(b"C$", M2Error::TrailingBond { pos: 1 })]
    #[case::trailing_stereo_bond_1(b"C/", M2Error::TrailingBond { pos: 1 })]
    #[case::trailing_stereo_bond_2(b"C\\", M2Error::TrailingBond { pos: 1 })]
    #[case::trailing_aromatic_bond(b"C:", M2Error::TrailingBond { pos: 1 })]
    #[case::branch_trailing_bond_1(b"C(C-)C", M2Error::TrailingBond { pos: 3 })]
    #[case::branch_trailing_bond_2(b"C(C=)C", M2Error::TrailingBond { pos: 3 })]
    #[case::branch_trailing_stereo_bond(b"CC(C/)CC", M2Error::TrailingBond { pos: 4 })]
    #[case::branch_trailing_aromatic_bond(b"C(C:)", M2Error::TrailingBond { pos: 3 })]
    #[case::group_trailing_bond_1(b"(C-)", M2Error::TrailingBond { pos: 2 })]
    #[case::group_trailing_bond_2(b"(C=)", M2Error::TrailingBond { pos: 2 })]
    #[case::group_trailing_stereo_bond(b"(C/)", M2Error::TrailingBond { pos: 2 })]
    #[case::group_trailing_aromatic_bond(b"(C:)", M2Error::TrailingBond { pos: 2 })]
    #[case::bond_after_group_1(b"(C)-", M2Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::bond_after_group_2(b"(C)=", M2Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::consecutive_bonds_1(b"C--C", M2Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bonds_2(b"C-=C", M2Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bonds_3(b"C-#C", M2Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bonds_4(b"C-$C", M2Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bonds_5(b"C-:C", M2Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_stereo_bonds_1(b"C//C", M2Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_stereo_bonds_2(b"C\\\\C", M2Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bond_and_stereo_bond_1(b"C-/C", M2Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bond_and_stereo_bond_2(b"C=\\C", M2Error::ConsecutiveBond { pos: 2 })]
    #[case::leading_bond_1(b"-C", M2Error::LeadingBond { pos: 0 })]
    #[case::leading_bond_2(b"=C", M2Error::LeadingBond { pos: 0 })]
    #[case::leading_bond_3(b"#C", M2Error::LeadingBond { pos: 0 })]
    #[case::leading_bond_4(b"$C", M2Error::LeadingBond { pos: 0 })]
    #[case::leading_aromatic_bond(b":C", M2Error::LeadingBond { pos: 0 })]
    #[case::leading_sterebond_1(b"/C", M2Error::LeadingBond { pos: 0 })]
    #[case::leading_sterebond_2(b"\\C", M2Error::LeadingBond { pos: 0 })]
    #[case::group_leading_bond_1(b"(-C)C", M2Error::LeadingBond { pos: 1 })]
    #[case::group_leading_bond_2(b"(=C)C", M2Error::LeadingBond { pos: 1 })]
    #[case::group_leading_bond_3(b"(#C)C", M2Error::LeadingBond { pos: 1 })]
    #[case::group_leading_bond_4(b"($C)C", M2Error::LeadingBond { pos: 1 })]
    #[case::group_leading_sterebond_1(b"(/C)C", M2Error::LeadingBond { pos: 1 })]
    #[case::group_leading_sterebond_2(b"(\\C)C", M2Error::LeadingBond { pos: 1 })]
    #[case::group_leading_aromatic_bond(b"(:C)C", M2Error::LeadingBond { pos: 1 })]
    fn m2_bonds_invalid(#[case] input: &[u8], #[case] expected: M2Error) {
        let err = parse_smiles_m2(input);
        assert!(err.is_err(), "{:?} should have failed", input);
        let err = err.unwrap_err();
        assert_eq!(err, expected);
    }

    #[rstest]
    #[case::wildcard(b"*", M2Error::UnsupportedToken { pos: 0 })]
    #[case::aromatic(b"c", M2Error::UnsupportedToken { pos: 0 })]
    #[case::bracket(b"[C]", M2Error::UnsupportedToken { pos: 0 })]
    #[case::ring(b"C1CC1", M2Error::UnsupportedToken { pos: 1 })]
    #[case::ring_percent(b"C%12CC1%2", M2Error::UnsupportedToken { pos: 1 })]
    #[case::component(b"CC.CC", M2Error::UnsupportedToken { pos: 2 })]
    #[case::whitespace_1(b"C ", M2Error::UnsupportedToken { pos: 1 })]
    #[case::whitespace_2(b"C\t", M2Error::UnsupportedToken { pos: 1 })]
    #[case::whitespace_3(b"C\n", M2Error::UnsupportedToken { pos: 1 })]
    #[case::whitespace_4(b"C\r", M2Error::UnsupportedToken { pos: 1 })]
    #[case::whitespace_5(b"C\r\n", M2Error::UnsupportedToken { pos: 1 })]
    fn m2_unimplemented(#[case] input: &[u8], #[case] expected: M2Error) {
        let err = parse_smiles_m2(input);
        assert!(err.is_err(), "{:?} should have failed", input);
        let err = err.unwrap_err();
        assert_eq!(err, expected);
    }
}
