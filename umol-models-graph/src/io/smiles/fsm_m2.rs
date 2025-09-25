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

        let mut last: Option<u32> = None;
        for _ in 0..n_initial {
            let curr = b.on_atom_fast(Element::C, true, false);
            if let Some(prev) = last {
                b.on_bond_single_fast(prev, curr);
            }
            last = Some(curr);
        }
        let base = last.expect("at least one initial atom required");

        let mut prev = base;
        for _ in 0..n_branch1 {
            let curr = b.on_atom_fast(Element::C, true, false);
            b.on_bond_single_fast(prev, curr);
            prev = curr;
        }

        let mut prev2 = base;
        for _ in 0..n_branch2 {
            let curr = b.on_atom_fast(Element::C, true, false);
            b.on_bond_single_fast(prev2, curr);
            prev2 = curr;
        }

        let mut mols = b.finish();
        mols.pop().unwrap_or_default()
    }

    fn build_bonds_c(n: usize, orders: &[BondOrder], dirs: &[Option<BondDir>]) -> Molecule {
        assert_eq!(orders.len(), n.saturating_sub(1));
        assert_eq!(dirs.len(), n.saturating_sub(1));

        let mut b = MoleculeBuilder::with_capacity(n, n.saturating_sub(1));
        let mut last: Option<u32> = None;
        for i in 0..n {
            let curr = b.on_atom_fast(Element::C, true, false);
            if let Some(s) = last {
                b.on_bond(
                    s,
                    curr,
                    BondData {
                        order: orders[i - 1],
                        dir: dirs[i - 1],
                    },
                );
            }
            last = Some(curr);
        }
        let mut mols = b.finish();
        mols.pop().unwrap_or_default()
    }

    fn build_stereo_double_bond(dir1: BondDir, dir2: BondDir) -> Molecule {
        let mut b = MoleculeBuilder::with_capacity(4, 3);
        let idx1 = b.on_atom_fast(Element::C, true, false);
        let idx2 = b.on_atom_fast(Element::C, true, false);
        let idx3 = b.on_atom_fast(Element::C, true, false);
        let idx4 = b.on_atom_fast(Element::C, true, false);
        b.on_bond(
            idx1,
            idx2,
            BondData {
                order: BondOrder::Single,
                dir: Some(dir1),
            },
        );
        b.on_bond(
            idx2,
            idx3,
            BondData {
                order: BondOrder::Double,
                dir: None,
            },
        );
        b.on_bond(
            idx3,
            idx4,
            BondData {
                order: BondOrder::Single,
                dir: Some(dir2),
            },
        );
        let mut mols = b.finish();
        mols.pop().unwrap_or_default()
    }

    #[rstest]
    #[case::empty_group(b"()", Molecule::default())]
    #[case::chain(b"CCC", build_chain_c(3))]
    fn m2_chain(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m2(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rstest]
    #[case::branch(b"CC(C)C", build_branch_c(2, 1, 1))]
    #[case::trailing_branch(b"C(CC)", build_chain_c(3))]
    #[case::top_level_group(b"(CCCC)", build_chain_c(4))]
    #[case::nested_group(b"((CC))", build_chain_c(2))]
    fn m2_tree(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m2(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single_bond(b"C-C", build_bonds_c(2, &[BondOrder::Single], &[None]))]
    #[case::double_bond(b"C=C", build_bonds_c(2, &[BondOrder::Double], &[None]))]
    #[case::triple_bond(b"C#C", build_bonds_c(2, &[BondOrder::Triple], &[None]))]
    #[case::quadruple_bond(b"C$C", build_bonds_c(2, &[BondOrder::Quadruple], &[None]))]
    #[case::aromatic_bond(b"C:C", build_bonds_c(2, &[BondOrder::Aromatic], &[None]))]
    #[case::up_bond(b"C/C", build_bonds_c(2, &[BondOrder::Single], &[Some(BondDir::Up)]))]
    #[case::down_bond(b"C\\C", build_bonds_c(2, &[BondOrder::Single], &[Some(BondDir::Down)]))]
    #[case::branch_leading_bond(b"CC(-C)C", build_branch_c(2, 1, 1))]
    #[case::branch_internal_bond(b"CC(C-C)C", build_branch_c(2, 2, 1))]
    #[case::branch_followed_by_bond(b"CC(C)-C", build_branch_c(2, 1, 1))]
    #[case::branch_cis_double_bond(b"CC(C)-C", build_branch_c(2, 1, 1))]
    #[case::branch_trans_double_bond_1(b"C/C=C/C", build_stereo_double_bond(BondDir::Up, BondDir::Up))]
    #[case::branch_trans_double_bond_2(b"C\\C=C\\C", build_stereo_double_bond(BondDir::Down, BondDir::Down))]
    #[case::branch_cis_double_bond_1(b"C\\C=C/C", build_stereo_double_bond(BondDir::Down, BondDir::Up))]
    #[case::branch_cis_double_bond_2(b"C/C=C\\C", build_stereo_double_bond(BondDir::Up, BondDir::Down))]
    #[case::cumulated_bonds(b"C=C=C", build_bonds_c(3, &[BondOrder::Double, BondOrder::Double], &[None, None]))]
    #[case::conjugated_bonds(b"C=CC=C", build_bonds_c(4, &[BondOrder::Double, BondOrder::Single, BondOrder::Double],
                             &[None, None, None]))]
    fn m2_bonds(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m2(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rstest]
    #[case::ring(b"C1CC1", M2Error::UnsupportedToken { pos: 1 })]
    #[case::component(b"CC.CC", M2Error::UnsupportedToken { pos: 2 })]
    #[case::unbalanced_closing_paren_1(b")C", M2Error::UnbalancedBranchClose { pos: 0 })]
    #[case::unbalanced_closing_paren_2(b"C)C", M2Error::UnbalancedBranchClose { pos: 1 })]
    #[case::unclosed_group(b"(C", M2Error::UnbalancedBranchOpen { pos: 0 })]
    #[case::unclosed_branch(b"C(C", M2Error::UnbalancedBranchOpen { pos: 1 })]
    #[case::empty_branch(b"C()", M2Error::EmptyBranch { pos: 2 })]
    #[case::empty_group_before_atom(b"()C", M2Error::EmptyGroup { pos: 1 })]
    #[case::trailing_bond_1(b"C-", M2Error::TrailingBond { pos: 1 })]
    #[case::trailing_bond_2(b"C=", M2Error::TrailingBond { pos: 1 })]
    #[case::trailing_bond_3(b"C#", M2Error::TrailingBond { pos: 1 })]
    #[case::trailing_bond_4(b"C$", M2Error::TrailingBond { pos: 1 })]
    #[case::trailing_stereo_bond_1(b"C/", M2Error::TrailingBond { pos: 1 })]
    #[case::trailing_stereo_bond_2(b"C\\", M2Error::TrailingBond { pos: 1 })]
    #[case::trailing_aromatic_bond(b"C:", M2Error::TrailingBond { pos: 1 })]
    #[case::trailing_bond(b"C(C-)C", M2Error::TrailingBond { pos: 3 })]
    #[case::trailing_stereo_bond(b"CC(C/)CC", M2Error::TrailingBond { pos: 4 })]
    #[case::group_trailing_bond(b"(C-)", M2Error::TrailingBond { pos: 2 })]
    #[case::group_trailing_stereo_bond(b"(C/)", M2Error::TrailingBond { pos: 2 })]
    #[case::bond_after_group(b"(C)-", M2Error::TrailingBond { pos: 3 })]
    #[case::consecutive_bonds_1(b"C--C", M2Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bonds_2(b"C-=C", M2Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bonds_3(b"C-#C", M2Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bonds_4(b"C-$C", M2Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bonds_5(b"C-:C", M2Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_sterebonds_1(b"C//C", M2Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_sterebonds_2(b"C\\\\C", M2Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bond_and_sterebond_1(b"C-/C", M2Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bond_and_sterebond_2(b"C=\\C", M2Error::ConsecutiveBond { pos: 2 })]
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
    fn m2_invalid(#[case] input: &[u8], #[case] expected: M2Error) {
        let err = parse_smiles_m2(input);
        assert!(err.is_err(), "{:?} should have failed", input);
        let err = err.unwrap_err();
        assert_eq!(err, expected);
    }
}
