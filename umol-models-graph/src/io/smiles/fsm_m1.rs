use umol_data::Element;

use crate::io::ir::builder::MoleculeBuilder;
use crate::io::ir::Molecule;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum M1Error {
    UnsupportedToken { pos: usize },
    UnbalancedBranchOpen { pos: usize },
    UnbalancedBranchClose { pos: usize },
    EmptyBranch { pos: usize },
    EmptyGroup { pos: usize },
    TopLevelGroupTrailing { pos: usize },
}

// Frames for parentheses: Branch attaches to a base atom; Group is top-level grouping
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

// M1: chains + branches (no rings, no aromatic, no charged/bracketed atoms)
pub fn parse_smiles_m1(input: &[u8]) -> Result<Molecule, M1Error> {
    let mut i = 0usize;
    let n = input.len();

    let mut builder = MoleculeBuilder::with_capacity(n.max(1), n.max(1).saturating_sub(1));
    let mut last_atom_idx: Option<u32> = None;

    let mut branch_stack: Vec<Frame> = Vec::new();

    while i < n {
        let b0 = input[i];

        // Parenthesis open: '(' — branch if we have an attach point, else top-level group
        if b0 == b'(' {
            match last_atom_idx {
                Some(idx) => branch_stack.push(Frame::Branch {
                    base: idx,
                    had_atom: false,
                    open_pos: i,
                }),
                None => branch_stack.push(Frame::Group {
                    had_atom: false,
                    open_pos: i,
                }),
            }
            i += 1;
            continue;
        }

        // Parenthesis close: ')' — close branch or group
        if b0 == b')' {
            let Some(frame) = branch_stack.pop() else {
                return Err(M1Error::UnbalancedBranchClose { pos: i });
            };
            match frame {
                Frame::Branch { base, had_atom, .. } => {
                    if !had_atom {
                        return Err(M1Error::EmptyBranch { pos: i });
                    }
                    last_atom_idx = Some(base);
                }
                Frame::Group { had_atom, .. } => {
                    if !had_atom {
                        // Permit top-level empty group only when it is the entire input ("()")
                        if i + 1 != n {
                            return Err(M1Error::EmptyGroup { pos: i });
                        }
                    } else {
                        // Disallow anything following a non-empty TOP-LEVEL group only
                        // i.e., when we just closed the outermost group (stack now empty)
                        if branch_stack.is_empty() && i + 1 != n {
                            return Err(M1Error::TopLevelGroupTrailing { pos: i });
                        }
                    }
                    // Do not change last_atom_idx for groups
                }
            }
            i += 1;
            continue;
        }

        // Recognize two-letter halogens first: Cl, Br
        if b0 == b'C' {
            if i + 1 < n && input[i + 1] == b'l' {
                let curr = builder.on_atom_fast(Element::Cl, true, false);
                if let Some(last) = last_atom_idx {
                    builder.on_bond_single_fast(last, curr);
                }
                last_atom_idx = Some(curr);
                if !branch_stack.is_empty() {
                    for f in branch_stack.iter_mut() {
                        match f {
                            Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                                *had_atom = true
                            }
                        }
                    }
                }
                i += 2;
                continue;
            }
            // Single C
            let curr = builder.on_atom_fast(Element::C, true, false);
            if let Some(last) = last_atom_idx {
                builder.on_bond_single_fast(last, curr);
            }
            last_atom_idx = Some(curr);
            if !branch_stack.is_empty() {
                for f in branch_stack.iter_mut() {
                    match f {
                        Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                            *had_atom = true
                        }
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
                    builder.on_bond_single_fast(last, curr);
                }
                last_atom_idx = Some(curr);
                if !branch_stack.is_empty() {
                    for f in branch_stack.iter_mut() {
                        match f {
                            Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                                *had_atom = true
                            }
                        }
                    }
                }
                i += 2;
                continue;
            }
            // Single B
            let curr = builder.on_atom_fast(Element::B, true, false);
            if let Some(last) = last_atom_idx {
                builder.on_bond_single_fast(last, curr);
            }
            last_atom_idx = Some(curr);
            if !branch_stack.is_empty() {
                for f in branch_stack.iter_mut() {
                    match f {
                        Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                            *had_atom = true
                        }
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
                builder.on_bond_single_fast(last, curr);
            }
            last_atom_idx = Some(curr);
            if !branch_stack.is_empty() {
                for f in branch_stack.iter_mut() {
                    match f {
                        Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                            *had_atom = true
                        }
                    }
                }
            }
            i += 1;
            continue;
        }

        return Err(M1Error::UnsupportedToken { pos: i });
    }

    if !branch_stack.is_empty() {
        // report last unmatched '('
        let pos = match branch_stack.last().unwrap() {
            Frame::Branch { open_pos, .. } | Frame::Group { open_pos, .. } => *open_pos,
        };
        return Err(M1Error::UnbalancedBranchOpen { pos });
    }

    let mut mols = builder.finish();
    Ok(mols.pop().unwrap_or_default())
}

#[cfg(test)]
mod tests {
    use rstest::*;

    use super::*;
    use crate::io::smiles::test_support::build_from_graph;

    #[rstest]
    #[case::empty(b"", Molecule::default())]
    #[case::chain_c_1(b"C", build_from_graph("C |"))]
    #[case::chain_c_5(b"CCCCC", build_from_graph("C C C C C | 0-1 1-2 2-3 3-4"))]
    #[case::chain_mixed_5(b"CClOBrN", build_from_graph("C Cl O Br N | 0-1 1-2 2-3 3-4"))]
    fn m1_chain(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m1(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
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
    fn m1_tree(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m1(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rstest]
    #[case::aromatic(b"c", M1Error::UnsupportedToken { pos: 0 })]
    #[case::bond_order(b"C-C", M1Error::UnsupportedToken { pos: 1 })]
    #[case::bracket(b"[C]", M1Error::UnsupportedToken { pos: 0 })]
    #[case::ring(b"C1CC1", M1Error::UnsupportedToken { pos: 1 })]
    #[case::component(b"CC.CC", M1Error::UnsupportedToken { pos: 2 })]
    #[case::stray_closing_paren(b")C", M1Error::UnbalancedBranchClose { pos: 0 })]
    #[case::unclosed_group(b"(C", M1Error::UnbalancedBranchOpen { pos: 0 })]
    #[case::unclosed_branch(b"C(C", M1Error::UnbalancedBranchOpen { pos: 1 })]
    #[case::empty_branch(b"C()", M1Error::EmptyBranch { pos: 2 })]
    #[case::empty_group_before_atom(b"()C", M1Error::EmptyGroup { pos: 1 })]
    #[case::two_top_level_groups(b"(C)(C)", M1Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::group_before_atom(b"(C)C", M1Error::TopLevelGroupTrailing { pos: 2 })]
    fn m1_invalid(#[case] input: &[u8], #[case] expected: M1Error) {
        let err = parse_smiles_m1(input);
        assert!(err.is_err(), "{:?} should have failed", input);
        let err = err.unwrap_err();
        assert_eq!(err, expected);
    }
}
