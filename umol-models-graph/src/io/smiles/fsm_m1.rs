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
}

// M1: chains + branches (no rings, no aromatic, no charged/bracketed atoms)
pub fn parse_smiles_m1(input: &[u8]) -> Result<Molecule, M1Error> {
    let mut i = 0usize;
    let n = input.len();

    let mut builder = MoleculeBuilder::with_capacity(n.max(1), n.max(1).saturating_sub(1));
    let mut last_atom_idx: Option<u32> = None;

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
    use rstest::rstest;
    use umol_data::Element;

    use super::*;
    use crate::io::ir::builder::MoleculeBuilder;

    fn build_chain_c(n: usize) -> Molecule {
        let mut b = MoleculeBuilder::with_capacity(n, n.saturating_sub(1));
        let mut last: Option<u32> = None;
        for _ in 0..n {
            let curr = b.on_atom_fast(Element::C, true, false);
            if let Some(s) = last {
                b.on_bond_single_fast(s, curr);
            }
            last = Some(curr);
        }
        let mut mols = b.finish();
        mols.pop().unwrap_or_default()
    }

    fn build_branch_c(n_initial: usize, n_branch1: usize, n_branch2: usize) -> Molecule {
        let total_atoms = n_initial + n_branch1 + n_branch2;
        let mut b = MoleculeBuilder::with_capacity(total_atoms, total_atoms.saturating_sub(1));
        if n_initial == 0 {
            return Molecule::default();
        }

        // Trunk up to and including the branching point
        let mut last: Option<u32> = None;
        for _ in 0..n_initial {
            let curr = b.on_atom_fast(Element::C, true, false);
            if let Some(prev) = last {
                b.on_bond_single_fast(prev, curr);
            }
            last = Some(curr);
        }
        let base = last.expect("at least one initial atom required");

        // Branch 1 from base
        let mut prev = base;
        for _ in 0..n_branch1 {
            let curr = b.on_atom_fast(Element::C, true, false);
            b.on_bond_single_fast(prev, curr);
            prev = curr;
        }

        // Branch 2 from base
        let mut prev2 = base;
        for _ in 0..n_branch2 {
            let curr = b.on_atom_fast(Element::C, true, false);
            b.on_bond_single_fast(prev2, curr);
            prev2 = curr;
        }

        let mut mols = b.finish();
        mols.pop().unwrap_or_default()
    }

    #[rstest]
    #[case::empty_group(b"()", Molecule::default())]
    #[case::chain(b"CCC", build_chain_c(3))]
    fn m1_chain(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m1(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rstest]
    #[case::branch(b"CC(C)C", build_branch_c(2, 1, 1))]
    #[case::trailing_branch(b"C(CC)", build_chain_c(3))]
    #[case::nested_group(b"((CC))", build_chain_c(2))]
    #[case::group_around_cccc(b"(CCCC)", build_chain_c(4))]
    fn m1_tree(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m1(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rstest]
    #[case::bond_order(b"C-C", M1Error::UnsupportedToken { pos: 1 })]
    #[case::bracket(b"[C]", M1Error::UnsupportedToken { pos: 0 })]
    #[case::ring(b"C1CC1", M1Error::UnsupportedToken { pos: 1 })]
    #[case::component(b"CC.CC", M1Error::UnsupportedToken { pos: 2 })]
    #[case::stray_closing_paren(b")C", M1Error::UnbalancedBranchClose { pos: 0 })]
    #[case::unclosed_group(b"(C", M1Error::UnbalancedBranchOpen { pos: 0 })]
    #[case::unclosed_branch(b"C(C", M1Error::UnbalancedBranchOpen { pos: 1 })]
    #[case::trailing_bond_in_branch(b"C(C-)C", M1Error::UnsupportedToken { pos: 3 })]
    #[case::empty_branch(b"C()", M1Error::EmptyBranch { pos: 2 })]
    #[case::empty_group_before_atom(b"()C", M1Error::EmptyGroup { pos: 1 })]
    fn invalid_groups(#[case] input: &[u8], #[case] expected: M1Error) {
        let err = parse_smiles_m1(input);
        assert!(err.is_err(), "{:?} should have failed", input);
        let err = err.unwrap_err();
        assert_eq!(err, expected);
    }
}
