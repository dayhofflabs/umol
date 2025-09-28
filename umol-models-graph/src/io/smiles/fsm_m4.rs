use umol_data::Element;

use crate::io::ir::builder::{BondData, MoleculeBuilder};
use crate::io::ir::{BondDir, BondOrder, Molecule};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum M4Error {
    UnsupportedToken { pos: usize },
    UnbalancedBranchOpen { pos: usize },
    UnbalancedBranchClose { pos: usize },
    EmptyBranch { pos: usize },
    EmptyGroup { pos: usize },
    TopLevelGroupTrailing { pos: usize },
    TrailingBond { pos: usize },
    ConsecutiveBond { pos: usize },
    LeadingBond { pos: usize },
    RingIndexInvalid { pos: usize },
    LeadingRing { pos: usize },
    RingBondDirConflict { pos: usize, open_pos: usize },
    RingBondOrderConflict { pos: usize, open_pos: usize },
    RingSelfLoop { pos: usize },
    RingTwoMember { pos: usize },
    RingUnclosed { open_pos: usize },
    LeadingDot { pos: usize },
    TrailingDot { pos: usize },
    ConsecutiveDot { pos: usize },
}

fn is_digit(b: u8) -> bool {
    (b'0'..=b'9').contains(&b)
}

#[derive(Debug, Clone, Copy)]
struct OpenRing {
    atom_id: u32,
    order: Option<BondOrder>,
    dir: Option<BondDir>,
    open_pos: usize,
    open_aromatic: bool,
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

pub fn parse_smiles_m4(input: &[u8]) -> Result<Molecule, M4Error> {
    let mut i = 0usize;
    let n = input.len();

    let mut builder = MoleculeBuilder::with_capacity(n.max(1), n.max(1).saturating_sub(1));
    let mut last_atom_idx: Option<u32> = None;
    let mut prev_atom_idx: Option<u32> = None; // immediate predecessor along current path
    let mut pending_bond: Option<(BondOrder, Option<BondDir>, usize)> = None;
    let mut last_aromatic: bool = false;

    let mut pstack: Vec<Frame> = Vec::new();
    let mut ring_table: [Option<OpenRing>; 100] = [None; 100]; // indices 0..99
    let mut just_closed_group: bool = false;

    while i < n {
        let b0 = input[i];
        if b0 != b'(' {
            just_closed_group = false;
        }

        // Parentheses
        if b0 == b'(' {
            if let Some((_, _, pos)) = pending_bond {
                return Err(M4Error::TrailingBond { pos });
            }
            if just_closed_group {
                // Starting a new top-level group after closing a non-empty group:
                // isolate components by clearing adjacency so atoms inside won't bond
                last_atom_idx = None;
                prev_atom_idx = None;
                pstack.push(Frame::Group {
                    had_atom: false,
                    open_pos: i,
                });
                just_closed_group = false;
            } else {
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
            }
            i += 1;
            continue;
        }
        if b0 == b')' {
            if let Some((_, _, pos)) = pending_bond {
                return Err(M4Error::TrailingBond { pos });
            }
            let Some(frame) = pstack.pop() else {
                return Err(M4Error::UnbalancedBranchClose { pos: i });
            };
            match frame {
                Frame::Branch { base, had_atom, .. } => {
                    if !had_atom {
                        return Err(M4Error::EmptyBranch { pos: i });
                    }
                    last_atom_idx = Some(base);
                    prev_atom_idx = None;
                }
                Frame::Group {
                    had_atom, open_pos, ..
                } => {
                    if !had_atom {
                        if i + 1 != n {
                            return Err(M4Error::EmptyGroup { pos: i });
                        }
                        // End-of-input: if group ends with a dot, treat as leading dot in component
                        if i > 0 && input[i - 1] == b'.' {
                            return Err(M4Error::LeadingDot { pos: i });
                        }
                        // Allow exactly "()" at EOI (group starting at 0). Otherwise, error.
                        if open_pos != 0 {
                            return Err(M4Error::EmptyGroup { pos: i });
                        }
                        // For exactly "()" keep no adjacency
                        last_atom_idx = None;
                        prev_atom_idx = None;
                        just_closed_group = false;
                    } else {
                        // Mark that we just closed a non-empty top-level group
                        just_closed_group = true;
                        // If this was the OUTERMOST group and there are trailing tokens, error
                        if pstack.is_empty() && i + 1 != n {
                            // Allow a connector (dot) to follow a top-level group in M4
                            let next = input[i + 1];
                            if next != b'.' {
                                return Err(M4Error::TopLevelGroupTrailing { pos: i });
                            }
                        }
                    }
                    // If group had atoms, leave last_atom_idx as-is to allow adjacency
                }
            }
            i += 1;
            continue;
        }

        // Component separator '.' resets adjacency but preserves ring_table and state
        if b0 == b'.' {
            if let Some((_, _, pos)) = pending_bond {
                return Err(M4Error::TrailingBond { pos });
            }
            // Invalid dot forms: leading, trailing, or consecutive dots
            if i == 0 {
                return Err(M4Error::LeadingDot { pos: i });
            }
            if i + 1 == n {
                return Err(M4Error::TrailingDot { pos: i });
            }
            if input[i + 1] == b'.' {
                return Err(M4Error::ConsecutiveDot { pos: i });
            }
            last_atom_idx = None;
            prev_atom_idx = None;
            last_aromatic = false;
            i += 1;
            continue;
        }

        // Ring tokens: single digit 0..9 and %DD
        if is_digit(b0) {
            if last_atom_idx.is_none() {
                return Err(M4Error::LeadingRing { pos: i });
            }
            let idx: usize = (b0 - b'0') as usize;
            let bond = pending_bond.take();
            let (order_opt, dir_opt) = bond.map_or((None, None), |(o, d, _)| (Some(o), d));
            let entry = &mut ring_table[idx];
            match entry.take() {
                None => {
                    *entry = Some(OpenRing {
                        atom_id: last_atom_idx.unwrap(),
                        order: order_opt,
                        dir: dir_opt,
                        open_pos: i,
                        open_aromatic: last_aromatic,
                    });
                }
                Some(open) => {
                    let b = last_atom_idx.unwrap();
                    if open.atom_id == b {
                        return Err(M4Error::RingSelfLoop { pos: i });
                    }
                    if prev_atom_idx == Some(open.atom_id) {
                        return Err(M4Error::RingTwoMember { pos: i });
                    }
                    if let (Some(d1), Some(d2)) = (open.dir, dir_opt) {
                        if d1 != d2 {
                            return Err(M4Error::RingBondDirConflict {
                                pos: i,
                                open_pos: open.open_pos,
                            });
                        }
                    }
                    // If both sides specified an order and they differ -> conflict
                    if let (Some(o1), Some(o2)) = (open.order, order_opt) {
                        if o1 != o2 {
                            return Err(M4Error::RingBondOrderConflict {
                                pos: i,
                                open_pos: open.open_pos,
                            });
                        }
                    }
                    // Direction combined with non-single order on the ring bond is unsupported
                    if open.dir.is_some() || dir_opt.is_some() {
                        let ord = open.order.or(order_opt).unwrap_or(BondOrder::Single);
                        if ord != BondOrder::Single {
                            return Err(M4Error::RingBondOrderConflict {
                                pos: i,
                                open_pos: open.open_pos,
                            });
                        }
                    }
                    let mut final_order = match (open.order, order_opt) {
                        (Some(o1), Some(o2)) => {
                            if o1 == o2 {
                                o1
                            } else {
                                o2
                            }
                        }
                        (Some(o), None) | (None, Some(o)) => o,
                        (None, None) => BondOrder::Single,
                    };
                    let final_dir = open.dir.or(dir_opt);
                    let a = open.atom_id;
                    let b = last_atom_idx.unwrap();
                    if final_order == BondOrder::Single && open.open_aromatic && last_aromatic {
                        final_order = BondOrder::Aromatic;
                    }
                    if a != b {
                        builder.on_bond(
                            a,
                            b,
                            BondData {
                                order: final_order,
                                dir: final_dir,
                            },
                        );
                    }
                }
            }
            i += 1;
            continue;
        }
        if b0 == b'%' {
            if i + 2 >= n || !is_digit(input[i + 1]) || !is_digit(input[i + 2]) {
                return Err(M4Error::RingIndexInvalid { pos: i });
            }
            if last_atom_idx.is_none() {
                return Err(M4Error::LeadingRing { pos: i });
            }
            let idx: usize = ((input[i + 1] - b'0') as usize) * 10 + (input[i + 2] - b'0') as usize;
            let bond = pending_bond.take();
            let (order_opt, dir_opt) = bond.map_or((None, None), |(o, d, _)| (Some(o), d));
            let entry = &mut ring_table[idx];
            match entry.take() {
                None => {
                    *entry = Some(OpenRing {
                        atom_id: last_atom_idx.unwrap(),
                        order: order_opt,
                        dir: dir_opt,
                        open_pos: i,
                        open_aromatic: last_aromatic,
                    });
                }
                Some(open) => {
                    let b = last_atom_idx.unwrap();
                    if open.atom_id == b {
                        return Err(M4Error::RingSelfLoop { pos: i });
                    }
                    if prev_atom_idx == Some(open.atom_id) {
                        return Err(M4Error::RingTwoMember { pos: i });
                    }
                    if let (Some(d1), Some(d2)) = (open.dir, dir_opt) {
                        if d1 != d2 {
                            return Err(M4Error::RingBondDirConflict {
                                pos: i,
                                open_pos: open.open_pos,
                            });
                        }
                    }
                    if let (Some(o1), Some(o2)) = (open.order, order_opt) {
                        if o1 != o2 {
                            return Err(M4Error::RingBondOrderConflict {
                                pos: i,
                                open_pos: open.open_pos,
                            });
                        }
                    }
                    if open.dir.is_some() || dir_opt.is_some() {
                        let ord = open.order.or(order_opt).unwrap_or(BondOrder::Single);
                        if ord != BondOrder::Single {
                            return Err(M4Error::RingBondOrderConflict {
                                pos: i,
                                open_pos: open.open_pos,
                            });
                        }
                    }
                    let mut final_order = match (open.order, order_opt) {
                        (Some(o1), Some(o2)) => {
                            if o1 == o2 {
                                o1
                            } else {
                                o2
                            }
                        }
                        (Some(o), None) | (None, Some(o)) => o,
                        (None, None) => BondOrder::Single,
                    };
                    let final_dir = open.dir.or(dir_opt);
                    let a = open.atom_id;
                    let b = last_atom_idx.unwrap();
                    if final_order == BondOrder::Single && open.open_aromatic && last_aromatic {
                        final_order = BondOrder::Aromatic;
                    }
                    if a != b {
                        builder.on_bond(
                            a,
                            b,
                            BondData {
                                order: final_order,
                                dir: final_dir,
                            },
                        );
                    }
                }
            }
            i += 3; // %DD consumed; %DDD will naturally see D next
            continue;
        }

        // Bond tokens
        if matches!(b0, b'-' | b'=' | b'#' | b'$' | b':' | b'/' | b'\\') {
            if pending_bond.is_some() {
                return Err(M4Error::ConsecutiveBond { pos: i });
            }
            if last_atom_idx.is_none() {
                return Err(M4Error::LeadingBond { pos: i });
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
                        // Implicit bond to non-aromatic atom is single regardless of last_aromatic
                        builder.on_bond_single_fast(last, curr);
                    }
                }
                last_atom_idx = Some(curr);
                last_aromatic = false;
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
                    // Implicit bond to non-aromatic atom is single regardless of last_aromatic
                    builder.on_bond_single_fast(last, curr);
                }
                prev_atom_idx = Some(last);
            }
            last_atom_idx = Some(curr);
            last_aromatic = false;
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
                        // Implicit bond to non-aromatic atom is single regardless of last_aromatic
                        builder.on_bond_single_fast(last, curr);
                    }
                }
                last_atom_idx = Some(curr);
                last_aromatic = false;
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
                    // Implicit bond to non-aromatic atom is single regardless of last_aromatic
                    builder.on_bond_single_fast(last, curr);
                }
                prev_atom_idx = Some(last);
            }
            last_atom_idx = Some(curr);
            last_aromatic = false;
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
                    // Implicit bond to non-aromatic atom is single regardless of last_aromatic
                    builder.on_bond_single_fast(last, curr);
                }
                prev_atom_idx = Some(last);
            }
            last_atom_idx = Some(curr);
            last_aromatic = false;
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

        // Aromatic bare atoms: b c n o p s
        if matches!(b0, b'b' | b'c' | b'n' | b'o' | b'p' | b's') {
            let element = match b0 {
                b'b' => Element::B,
                b'c' => Element::C,
                b'n' => Element::N,
                b'o' => Element::O,
                b'p' => Element::P,
                _ => Element::S,
            };
            let curr = builder.on_atom_fast(element, true, true);
            if let Some(last) = last_atom_idx {
                if let Some((order, dir, _pos)) = pending_bond.take() {
                    builder.on_bond(last, curr, BondData { order, dir });
                } else {
                    if last_aromatic {
                        builder.on_bond(
                            last,
                            curr,
                            BondData {
                                order: BondOrder::Aromatic,
                                dir: None,
                            },
                        );
                    } else {
                        builder.on_bond_single_fast(last, curr);
                    }
                }
                prev_atom_idx = Some(last);
            }
            last_atom_idx = Some(curr);
            last_aromatic = true;
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

        return Err(M4Error::UnsupportedToken { pos: i });
    }

    if pending_bond.is_some() {
        let (_, _, pos) = pending_bond.unwrap();
        return Err(M4Error::TrailingBond { pos });
    }

    if !pstack.is_empty() {
        let pos = match pstack.last().unwrap() {
            Frame::Branch { open_pos, .. } | Frame::Group { open_pos, .. } => *open_pos,
        };
        return Err(M4Error::UnbalancedBranchOpen { pos });
    }

    // Unclosed rings: report the unmatched with the latest pos_open
    let mut last_open: Option<usize> = None;
    for entry in ring_table.iter().flatten() {
        match last_open {
            None => last_open = Some(entry.open_pos),
            Some(p) => {
                if entry.open_pos > p {
                    last_open = Some(entry.open_pos)
                }
            }
        }
    }
    if let Some(pos_open) = last_open {
        return Err(M4Error::RingUnclosed { open_pos: pos_open });
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

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(b"", Molecule::default())]
    #[case::chain_c_1(b"C", build_from_graph("C |"))]
    #[case::chain_c_5(b"CCCCC", build_from_graph("C C C C C | 0-1 1-2 2-3 3-4"))]
    #[case::aromatic_c_6(b"cccccc", build_from_graph("C* C* C* C* C* C* | 0-1: 1-2: 2-3: 3-4: 4-5:"))]
    #[case::chain_mixed_5(b"CClOBrN", build_from_graph("C Cl O Br N | 0-1 1-2 2-3 3-4"))]
    fn m4_chain(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m4(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rstest]
    #[case::non_ascii(b"\xf0\x9f\x9c\x8d", M4Error::UnsupportedToken { pos: 0 })]
    #[case::comma(b",", M4Error::UnsupportedToken { pos: 0 })]
    #[case::semicolon(b";", M4Error::UnsupportedToken { pos: 0 })]
    #[case::question_mark(b"?", M4Error::UnsupportedToken { pos: 0 })]
    #[case::caret(b"^", M4Error::UnsupportedToken { pos: 0 })]
    #[case::pipe(b"|", M4Error::UnsupportedToken { pos: 0 })]
    #[case::open_angle_bracket(b"<", M4Error::UnsupportedToken { pos: 0 })]
    #[case::close_angle_bracket(b"<", M4Error::UnsupportedToken { pos: 0 })]
    #[case::open_brace(b"{", M4Error::UnsupportedToken { pos: 0 })]
    #[case::close_brace(b"}", M4Error::UnsupportedToken { pos: 0 })]
    #[case::single_quote(b"'", M4Error::UnsupportedToken { pos: 0 })]
    #[case::double_quote(b"\"", M4Error::UnsupportedToken { pos: 0 })]
    #[case::backtick(b"`", M4Error::UnsupportedToken { pos: 0 })]
    #[case::tilde(b"~", M4Error::UnsupportedToken { pos: 0 })]
    #[case::exclamation_mark(b"!", M4Error::UnsupportedToken { pos: 0 })]
    #[case::ampersand(b"&", M4Error::UnsupportedToken { pos: 0 })]
    #[case::underscore(b"_", M4Error::UnsupportedToken { pos: 0 })]
    #[case::bare_chirality(b"C@", M4Error::UnsupportedToken { pos: 1 })]
    #[case::bare_hcount(b"CH", M4Error::UnsupportedToken { pos: 1 })]
    fn m4_tokens_invalid(#[case] input: &[u8], #[case] expected: M4Error) {
        let res = parse_smiles_m4(input);
        assert!(res.is_err(), "{:?} should have failed", input);
        let err = res.unwrap_err();
        assert_eq!(err, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty_group(b"()", Molecule::default())]
    #[case::group_c_1(b"(C)", build_from_graph("C |"))]
    #[case::group_c_1_aromatic(b"(c)", build_from_graph("C* |"))]
    #[case::group_c_4(b"(CCCC)", build_from_graph("C C C C | 0-1 1-2 2-3"))]
    #[case::group_nested(b"((CC))", build_from_graph("C C | 0-1"))]
    #[case::branch_c_111(b"C(C)(C)", build_from_graph("C C C | 0-1 0-2"))]
    #[case::branch_c_211(b"CC(C)C", build_from_graph("C C C C | 0-1 1-2 1-3"))]
    #[case::branch_c_222_aromatic(b"cc(cc)cc", build_from_graph("C* C* C* C* C* C* | 0-1: 1-2: 2-3: 1-4: 4-5:"))]
    #[case::trailing_branch(b"C(CC)", build_from_graph("C C C | 0-1 1-2"))]
    fn m4_tree(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m4(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rstest]
    #[case::unbalanced_closing_paren_1(b")C", M4Error::UnbalancedBranchClose { pos: 0 })]
    #[case::unbalanced_closing_paren_2(b"C)C", M4Error::UnbalancedBranchClose { pos: 1 })]
    #[case::unclosed_group(b"(C", M4Error::UnbalancedBranchOpen { pos: 0 })]
    #[case::unclosed_branch(b"C(C", M4Error::UnbalancedBranchOpen { pos: 1 })]
    #[case::empty_branch(b"C()", M4Error::EmptyBranch { pos: 2 })]
    #[case::empty_group_before_atom(b"()C", M4Error::EmptyGroup { pos: 1 })]
    #[case::group_before_atom(b"(C)C", M4Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::group_before_atom_aromatic(b"(c)c", M4Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::two_top_level_groups(b"(C)(C)", M4Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::two_top_level_groups_aromatic(b"(c)(c)", M4Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::three_top_level_groups(b"(C)(C)(C)", M4Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::three_top_level_groups_aromatic(b"(c)(c)(c)", M4Error::TopLevelGroupTrailing { pos: 2 })]
    fn m4_tree_invalid(#[case] input: &[u8], #[case] expected: M4Error) {
        let err = parse_smiles_m4(input);
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
    #[case::single_bond_aromatic(b"c-c", build_from_graph("C* C* | 0-1:-"))]
    #[case::double_bond_aromatic(b"c=c", build_from_graph("C* C* | 0-1:="))]
    #[case::triple_bond_aromatic(b"c#c", build_from_graph("C* C* | 0-1:#"))]
    #[case::quadruple_bond_aromatic(b"c$c", build_from_graph("C* C* | 0-1:$"))]
    #[case::aromatic_bond_aromatic(b"c:c", build_from_graph("C* C* | 0-1::"))]
    #[case::up_bond_aromatic(b"c/c", build_from_graph("C* C* | 0-1:/"))]
    #[case::down_bond_aromatic(b"c\\c", build_from_graph("C* C* | 0-1:\\"))]
    #[case::cumulated_bonds(b"C=C=C", build_from_graph("C C C | 0-1:= 1-2:="))]
    #[case::conjugated_bonds(b"C=CC=C", build_from_graph("C C C C | 0-1:= 1-2:- 2-3:="))]
    #[case::cumulated_bonds_aromatic(b"c=c=c", build_from_graph("C* C* C* | 0-1:= 1-2:="))]
    #[case::conjugated_bonds_aromatic(b"c=c-c=c", build_from_graph("C* C* C* C* | 0-1:= 1-2:- 2-3:="))]
    #[case::branch_leading_bond(b"CC(-C)C", build_from_graph("C C C C | 0-1 1-2 1-3"))]
    #[case::branch_leading_double_bond(b"CC(=C)C", build_from_graph("C C C C | 0-1 1-2:= 1-3"))]
    #[case::branch_internal_bond(b"CC(C-C)C", build_from_graph("C C C C C | 0-1 1-2 2-3 1-4"))]
    #[case::branch_internal_double_bond(b"CC(C=C)C", build_from_graph("C C C C C | 0-1 1-2 2-3:= 1-4"))]
    #[case::branch_followed_by_bond(b"CC(C)-C", build_from_graph("C C C C | 0-1 1-2 1-3"))]
    #[case::branch_followed_by_double_bond(b"CC(C)=C", build_from_graph("C C C C | 0-1 1-2 1-3:="))]
    #[case::branch_leading_bond_aromatic(b"cc(:c)c", build_from_graph("C* C* C* C* | 0-1: 1-2: 1-3:"))]
    #[case::branch_internal_bond_aromatic(b"cc(c:c)c", build_from_graph("C* C* C* C* C* | 0-1: 1-2: 2-3: 1-4:"))]
    #[case::branch_followed_by_bond_aromatic(b"cc(c):c", build_from_graph("C* C* C* C* | 0-1: 1-2: 1-3:"))]
    #[case::branch_trans_double_bond_1(b"C/C=C/C", build_from_graph("C C C C | 0-1:/ 1-2:= 2-3:/"))]
    #[case::branch_trans_double_bond_2(b"C\\C=C\\C", build_from_graph("C C C C | 0-1:\\ 1-2:= 2-3:\\"))]
    #[case::branch_cis_double_bond_1(b"C\\C=C/C", build_from_graph("C C C C | 0-1:\\ 1-2:= 2-3:/"))]
    #[case::branch_cis_double_bond_2(b"C/C=C\\C", build_from_graph("C C C C | 0-1:/ 1-2:= 2-3:\\"))]
    #[case::ring_single_bond(b"C-1-C-C-1", build_from_graph("C C C | 0-1 1-2 0-2"))]
    #[case::ring_double_bond_1(b"C1-C=C1", build_from_graph("C C C | 0-1 1-2:= 0-2"))]
    #[case::ring_double_bond_2(b"C1-CC=1", build_from_graph("C C C | 0-1 1-2 0-2:="))]
    #[case::ring_double_bond_3(b"C=1-CC1", build_from_graph("C C C | 0-1 1-2 0-2:="))]
    #[case::ring_double_bond_4(b"C=1-C-C=1", build_from_graph("C C C | 0-1 1-2 0-2:="))]
    #[case::ring_double_bond_unilateral_close(b"C1CC=1", build_from_graph("C C C | 0-1 1-2 0-2:="))]
    #[case::ring_double_bond_unilateral_open(b"C=1CC1", build_from_graph("C C C | 0-1 1-2 0-2:="))]
    #[case::ring_triple_bond(b"C1-C-C#1", build_from_graph("C C C | 0-1 1-2 0-2:#"))]
    #[case::ring_quadruple_bond(b"C1-C-C$1", build_from_graph("C C C | 0-1 1-2 0-2:$"))]
    #[case::ring_aromatic_bond(b"c1:c:c:1", build_from_graph("C* C* C* | 0-1: 1-2: 0-2:"))]
    #[case::ring_up_bond_1(b"C1CC/1", build_from_graph("C C C | 0-1 1-2 0-2:/"))]
    #[case::ring_up_bond_2(b"C/1CC1", build_from_graph("C C C | 0-1 1-2 0-2:/"))]
    #[case::ring_up_bond_3(b"C/1CC/1", build_from_graph("C C C | 0-1 1-2 0-2:/"))]
    #[case::ring_down_bond(b"C1CC\\1", build_from_graph("C C C | 0-1 1-2 0-2:\\"))]
    #[case::ring_down_bond_both(b"C\\1CC\\1", build_from_graph("C C C | 0-1 1-2 0-2:\\"))]
    #[case::ring_up_bond_percent_open(b"C/%12CC%12", build_from_graph("C C C | 0-1 1-2 0-2:/"))]
    #[case::ring_up_bond_percent_close(b"C%12CC/%12", build_from_graph("C C C | 0-1 1-2 0-2:/"))]
    #[case::ring_down_bond_percent_both(b"C\\%12CC\\%12", build_from_graph("C C C | 0-1 1-2 0-2:\\"))]
    #[case::ring_between_bonds(b"C1CC-1-C", build_from_graph("C C C C | 0-1 1-2 0-2 2-3"))]
    fn m4_bonds(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m4(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rstest]
    #[case::trailing_bond_1(b"C-", M4Error::TrailingBond { pos: 1 })]
    #[case::trailing_bond_2(b"C=", M4Error::TrailingBond { pos: 1 })]
    #[case::trailing_bond_3(b"C#", M4Error::TrailingBond { pos: 1 })]
    #[case::trailing_bond_4(b"C$", M4Error::TrailingBond { pos: 1 })]
    #[case::trailing_stereo_bond_1(b"C/", M4Error::TrailingBond { pos: 1 })]
    #[case::trailing_stereo_bond_2(b"C\\", M4Error::TrailingBond { pos: 1 })]
    #[case::trailing_aromatic_bond(b"C:", M4Error::TrailingBond { pos: 1 })]
    #[case::trailing_single_bond_aromatic(b"c-", M4Error::TrailingBond { pos: 1 })]
    #[case::trailing_double_bond_aromatic(b"c=", M4Error::TrailingBond { pos: 1 })]
    #[case::trailing_triple_bond_aromatic(b"c#", M4Error::TrailingBond { pos: 1 })]
    #[case::trailing_quadruple_bond_aromatic(b"c$", M4Error::TrailingBond { pos: 1 })]
    #[case::trailing_up_bond_aromatic(b"c/", M4Error::TrailingBond { pos: 1 })]
    #[case::trailing_down_bond_aromatic(b"c\\", M4Error::TrailingBond { pos: 1 })]
    #[case::trailing_aromatic_bond_aromatic(b"c:", M4Error::TrailingBond { pos: 1 })]
    #[case::branch_trailing_bond_1(b"C(C-)C", M4Error::TrailingBond { pos: 3 })]
    #[case::branch_trailing_bond_2(b"C(C=)C", M4Error::TrailingBond { pos: 3 })]
    #[case::branch_trailing_stereo_bond(b"CC(C/)CC", M4Error::TrailingBond { pos: 4 })]
    #[case::branch_trailing_aromatic_bond(b"C(C:)", M4Error::TrailingBond { pos: 3 })]
    #[case::group_trailing_bond_1(b"(C-)", M4Error::TrailingBond { pos: 2 })]
    #[case::group_trailing_bond_2(b"(C=)", M4Error::TrailingBond { pos: 2 })]
    #[case::group_trailing_stereo_bond(b"(C/)", M4Error::TrailingBond { pos: 2 })]
    #[case::group_trailing_aromatic_bond(b"(C:)", M4Error::TrailingBond { pos: 2 })]
    #[case::bond_after_group_1(b"(C)-", M4Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::bond_after_group_2(b"(C)=", M4Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::consecutive_bonds_1(b"C--C", M4Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bonds_2(b"C-=C", M4Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bonds_3(b"C-#C", M4Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bonds_4(b"C-$C", M4Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bonds_5(b"C-:C", M4Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_stereo_bonds_1(b"C//C", M4Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_stereo_bonds_2(b"C\\\\C", M4Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bond_and_stereo_bond_1(b"C-/C", M4Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bond_and_stereo_bond_2(b"C=\\C", M4Error::ConsecutiveBond { pos: 2 })]
    #[case::leading_bond_1(b"-C", M4Error::LeadingBond { pos: 0 })]
    #[case::leading_bond_2(b"=C", M4Error::LeadingBond { pos: 0 })]
    #[case::leading_bond_3(b"#C", M4Error::LeadingBond { pos: 0 })]
    #[case::leading_bond_4(b"$C", M4Error::LeadingBond { pos: 0 })]
    #[case::leading_aromatic_bond(b":C", M4Error::LeadingBond { pos: 0 })]
    #[case::leading_sterebond_1(b"/C", M4Error::LeadingBond { pos: 0 })]
    #[case::leading_sterebond_2(b"\\C", M4Error::LeadingBond { pos: 0 })]
    #[case::leading_single_bond_aromatic(b"-c", M4Error::LeadingBond { pos: 0 })]
    #[case::leading_double_bond_aromatic(b"=c", M4Error::LeadingBond { pos: 0 })]
    #[case::leading_triple_bond_aromatic(b"#c", M4Error::LeadingBond { pos: 0 })]
    #[case::leading_quadruple_bond_aromatic(b"$c", M4Error::LeadingBond { pos: 0 })]
    #[case::leading_up_bond_aromatic(b"/c", M4Error::LeadingBond { pos: 0 })]
    #[case::leading_down_bond_aromatic(b"\\c", M4Error::LeadingBond { pos: 0 })]
    #[case::leading_aromatic_bond_aromatic(b":c", M4Error::LeadingBond { pos: 0 })]
    #[case::group_leading_bond_1(b"(-C)C", M4Error::LeadingBond { pos: 1 })]
    #[case::group_leading_bond_2(b"(=C)C", M4Error::LeadingBond { pos: 1 })]
    #[case::group_leading_bond_3(b"(#C)C", M4Error::LeadingBond { pos: 1 })]
    #[case::group_leading_bond_4(b"($C)C", M4Error::LeadingBond { pos: 1 })]
    #[case::group_leading_sterebond_1(b"(/C)C", M4Error::LeadingBond { pos: 1 })]
    #[case::group_leading_sterebond_2(b"(\\C)C", M4Error::LeadingBond { pos: 1 })]
    #[case::group_leading_aromatic_bond(b"(:C)C", M4Error::LeadingBond { pos: 1 })]
    fn m4_bonds_invalid(#[case] input: &[u8], #[case] expected: M4Error) {
        let err = parse_smiles_m4(input);
        assert!(err.is_err(), "{:?} should have failed", input);
        let err = err.unwrap_err();
        assert_eq!(err, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ring_c_3(b"C1CC1", build_from_graph("C C C | 0-1 1-2 0-2"))]
    #[case::ring_c_10(b"C1CCCCCCCCC1", build_from_graph("C C C C C C C C C C | 0-1 1-2 2-3 3-4 4-5 5-6 6-7 7-8 8-9 0-9"))]
    #[case::ring_aromatic_c_6(b"c1ccccc1", build_from_graph("C* C* C* C* C* C* | 0-1: 1-2: 2-3: 3-4: 4-5: 0-5:"))]
    #[case::ring_index_0(b"C0CC0", build_from_graph("C C C | 0-1 1-2 0-2"))]
    #[case::ring_index_percent(b"C%12CC%12", build_from_graph("C C C | 0-1 1-2 0-2"))]
    #[case::ring_index_numeric_equiv_1(b"C1CC%01", build_from_graph("C C C | 0-1 1-2 0-2"))]
    #[case::ring_index_numeric_equiv_0(b"C0CC%00", build_from_graph("C C C | 0-1 1-2 0-2"))]
    #[case::ring_index_numeric_equiv_9(b"C9CC%09", build_from_graph("C C C | 0-1 1-2 0-2"))]
    #[case::ring_index_max_percent(b"C%99CC%99", build_from_graph("C C C | 0-1 1-2 0-2"))]
    #[case::two_rings_bonded_0(b"C1CC1C2CC2", build_from_graph("C C C C C C | 0-1 1-2 0-2 2-3 3-4 4-5 3-5"))]
    #[case::two_rings_bonded_0_aromatic(b"c1cc1c2cc2", build_from_graph("C* C* C* C* C* C* | 0-1: 1-2: 0-2: 2-3: 3-4: 4-5: 3-5:"))]
    #[case::two_rings_bonded_0_aromatic_1(b"c1cc1C2CC2", build_from_graph("C* C* C* C C C | 0-1: 1-2: 0-2: 2-3 3-4 4-5 3-5"))]
    #[case::two_rings_index_reused(b"C1CC1C1CC1", build_from_graph("C C C C C C | 0-1 1-2 0-2 2-3 3-4 4-5 3-5"))]
    #[case::two_rings_bonded_2(b"C1CC1CCC2CC2", build_from_graph("C C C C C C C C | 0-1 1-2 0-2 2-3 3-4 4-5 5-6 6-7 5-7"))]
    #[case::two_rings_spiro(b"C1CC12CC2", build_from_graph("C C C C C | 0-1 1-2 0-2 2-3 3-4 2-4"))]
    #[case::two_rings_fused(b"C12CC1C2", build_from_graph("C C C C | 0-1 1-2 0-2 2-3 0-3"))]
    #[case::two_rings_bridged(b"C12CC(C2)C1", build_from_graph("C C C C C | 0-1 1-2 2-3 0-3 2-4 0-4"))]
    #[case::two_rings_fused_aromatic(b"c12ccccc1cccc2", build_from_graph("C* C* C* C* C* C* C* C* C* C* | 0-1: 1-2: 2-3: 3-4: 4-5: 0-5: 5-6: 6-7: 7-8: 8-9: 0-9:"))]
    #[case::three_rings_fused(b"C12C3C1C32", build_from_graph("C C C C | 0-1 1-2 0-2 2-3 1-3 0-3"))]
    #[case::ring_group(b"(C1CC1)", build_from_graph("C C C | 0-1 1-2 0-2"))]
    #[case::ring_branch_1(b"CC(C1)(C1)", build_from_graph("C C C C | 0-1 1-2 1-3 2-3"))]
    #[case::ring_branch_2(b"C(C1)CC1", build_from_graph("C C C C | 0-1 0-2 2-3 1-3 "))]
    #[case::substituted_ring_1(b"CC1CC1", build_from_graph("C C C C | 0-1 1-2 2-3 1-3"))]
    #[case::substituted_ring_2(b"C1(C)CC1", build_from_graph("C C C C | 0-1 0-2 2-3 0-3"))]
    #[case::substituted_ring_3(b"C(C)1CC1", build_from_graph("C C C C | 0-1 0-2 2-3 0-3"))]
    #[case::substituted_ring_4(b"C1C(C)C1", build_from_graph("C C C C | 0-1 1-2 1-3 0-3"))]
    #[case::substituted_ring_5(b"C1CC(C)1", build_from_graph("C C C C | 0-1 1-2 2-3 0-2"))]
    #[case::substituted_ring_6(b"C1CC1C", build_from_graph("C C C C | 0-1 1-2 0-2 2-3"))]
    #[case::substituted_ring_7(b"C1CC1(C)", build_from_graph("C C C C | 0-1 1-2 0-2 2-3"))]
    #[case::substituted_ring_aromatic(b"c1c(c)c1", build_from_graph("C* C* C* C* | 0-1: 1-2: 1-3: 0-3:"))]
    #[case::substituted_ring_branch_1(b"C1C(C(C)C)C1", build_from_graph("C C C C C C | 0-1 1-2 2-3 2-4 1-5 0-5"))]
    fn m4_ring(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m4(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rstest]
    #[case::leading_ring(b"1C", M4Error::LeadingRing { pos: 0 })]
    #[case::bad_percent_short(b"C%1", M4Error::RingIndexInvalid { pos: 1 })]
    #[case::bad_percent_char(b"C%1a", M4Error::RingIndexInvalid { pos: 1 })]
    #[case::bad_percent_eoi(b"C%", M4Error::RingIndexInvalid { pos: 1 })]
    #[case::bad_percent_zero(b"C%0", M4Error::RingIndexInvalid { pos: 1 })]
    #[case::ring_self_loop(b"C11", M4Error::RingSelfLoop { pos: 2 })]
    #[case::ring_two_member(b"C1C1", M4Error::RingTwoMember { pos: 3 })]
    #[case::ring_bond_order_conflict_3(b"C=1CC#1", M4Error::RingBondOrderConflict { pos: 6, open_pos: 2 })]
    #[case::ring_bond_order_conflict_4(b"C/1CC=1", M4Error::RingBondOrderConflict { pos: 6, open_pos: 2 })]
    #[case::ring_bond_order_conflict_5(b"C\\1CC=1", M4Error::RingBondOrderConflict { pos: 6, open_pos: 2 })]
    #[case::ring_bond_order_conflict_6(b"C=1CC/1", M4Error::RingBondOrderConflict { pos: 6, open_pos: 2 })]
    #[case::ring_bond_order_conflict_7(b"C=1CC\\1", M4Error::RingBondOrderConflict { pos: 6, open_pos: 2 })]
    #[case::ring_bond_order_conflict_8(b"C=%10CC#%10", M4Error::RingBondOrderConflict { pos: 8, open_pos: 2 })]
    #[case::ring_bond_dir_conflict_1(b"C/1CC\\1", M4Error::RingBondDirConflict { pos: 6, open_pos: 2 })]
    #[case::ring_bond_dir_conflict_2(b"C\\1CC/1", M4Error::RingBondDirConflict { pos: 6, open_pos: 2 })]
    #[case::ring_bond_dir_conflict_3(b"C/%12CC\\%12", M4Error::RingBondDirConflict { pos: 8, open_pos: 2 })]
    #[case::ring_bond_dir_conflict_4(b"C\\%12CC/%12", M4Error::RingBondDirConflict { pos: 8, open_pos: 2 })]
    #[case::ring_unclosed_1(b"C1CC", M4Error::RingUnclosed { open_pos: 1 })]
    #[case::ring_unclosed_2(b"C1CC1C1", M4Error::RingUnclosed { open_pos: 6 })]
    #[case::ring_after_group(b"(C1CCC)1", M4Error::TopLevelGroupTrailing { pos : 6})]
    #[case::two_top_level_groups_rings(b"(C1CC1)(C2CC2)", M4Error::TopLevelGroupTrailing { pos: 6 })]
    fn m4_ring_invalid(#[case] input: &[u8], #[case] expected: M4Error) {
        let err = parse_smiles_m4(input);
        assert!(err.is_err(), "{:?} should have failed", input);
        let err = err.unwrap_err();
        assert_eq!(err, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::components_2(b"CC.CC", build_from_graph("C C C C | 0-1 2-3"))]
    #[case::components_5(b"C.C.C.C.C", build_from_graph("C C C C C | "))]
    #[case::ring_components_1(b"C1.CC1", build_from_graph("C C C | 1-2 0-2"))]
    #[case::ring_components_2(b"C%12.CC%12", build_from_graph("C C C | 1-2 0-2"))]
    #[case::ring_components_aromatic(b"c1.ccccc1", build_from_graph("C* C* C* C* C* C* | 1-2: 2-3: 3-4: 4-5: 0-5:"))]
    #[case::branch_components(b"C(C.C)", build_from_graph("C C C | 0-1"))]
    #[case::leading_dot_in_branch_1(b"C(.C)", build_from_graph("C C | "))]
    #[case::leading_dot_in_branch_2(b"C(.C)(C)", build_from_graph("C C C | 0-2"))]
    #[case::leading_dot_in_branch_3(b"C(C)(.C)", build_from_graph("C C C | 0-1"))]
    #[case::trailing_dot_in_branch_1(b"C(C.)", build_from_graph("C C | 0-1"))]
    #[case::trailing_dot_in_branch_2(b"C(C.)C", build_from_graph("C C C | 0-1 0-2"))]
    #[case::trailing_dot_in_branch_3(b"C(C.)(C)", build_from_graph("C C C | 0-1 0-2"))]
    #[case::group_components_1(b"(C.CC.C)", build_from_graph("C C C C | 1-2"))]
    #[case::group_components_2(b"(CC).(CC)", build_from_graph("C C C C | 0-1 2-3"))]
    #[case::group_components_3(b"(C.C).C", build_from_graph("C C C |"))]
    #[case::group_components_4(b"C.(C).C", build_from_graph("C C C |"))]
    #[case::group_components_5(b"C.C.(C)", build_from_graph("C C C |"))]
    #[case::leading_dot_in_group_1(b"(.CC)", build_from_graph("C C | 0-1"))]
    #[case::leading_dot_in_group_2(b"(.CC).(CC)", build_from_graph("C C C C | 0-1 2-3"))]
    #[case::leading_dot_in_group_3(b"(CC).(.CC)", build_from_graph("C C C C | 0-1 2-3"))]
    #[case::leading_dot_in_group_4(b"C(.C.C)", build_from_graph("C C C |"))]
    #[case::trailing_dot_in_group_1(b"(CC.)", build_from_graph("C C | 0-1"))]
    #[case::trailing_dot_in_group_2(b"(CC.).CC", build_from_graph("C C C C | 0-1 2-3"))]
    #[case::trailing_dot_in_group_3(b"(CC).(CC.)", build_from_graph("C C C C | 0-1 2-3"))]
    #[case::group_ring_components_1(b"(CC1.C1)", build_from_graph("C C C | 0-1 1-2"))]
    #[case::group_ring_components_2(b"C1.(C).CC1", build_from_graph("C C C C | 2-3 0-3 "))]
    #[case::group_ring_components_3(b"C%12.(C).CC%12", build_from_graph("C C C C | 2-3 0-3 "))]
    #[case::rings_across_multiple_dots_digit(b"C1.C.CC1", build_from_graph("C C C C | 2-3 0-3"))]
    #[case::rings_across_multiple_dots_percent(b"C%12.C.CC%12", build_from_graph("C C C C | 2-3 0-3"))]
    #[case::ring_double_unilateral_open(b"C=1.CC1", build_from_graph("C C C | 1-2 0-2:="))]
    #[case::ring_double_unilateral_close(b"C1.CC=1", build_from_graph("C C C | 1-2 0-2:="))]
    #[case::ring_dir_up_both(b"C/1.CC/1", build_from_graph("C C C | 1-2 0-2:/"))]
    #[case::ring_dir_down_both(b"C\\1.CC\\1", build_from_graph("C C C | 1-2 0-2:\\"))]
    #[case::ring_dir_up_both_percent(b"C/%12.CC/%12", build_from_graph("C C C | 1-2 0-2:/"))]
    #[case::ring_dir_down_both_percent(b"C\\%12.CC\\%12", build_from_graph("C C C | 1-2 0-2:\\"))]
    #[case::branch_multiple_components(b"C(.C.C)", build_from_graph("C C C |"))]
    #[case::groups_leading_dot_both(b"(.C).(.C)", build_from_graph("C C |"))]
    #[case::group_leading_dot_middle(b"C.(.C).C", build_from_graph("C C C |"))]
    fn m4_components(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m4(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rstest]
    #[case::leading_dot_1(b".", M4Error::LeadingDot { pos: 0 })]
    #[case::leading_dot_2(b".C", M4Error::LeadingDot { pos: 0 })]
    #[case::leading_dot_3(b"..C", M4Error::LeadingDot { pos: 0 })]
    #[case::leading_dot_4(b".C.", M4Error::LeadingDot { pos: 0 })]
    #[case::trailing_dot_1(b"C.", M4Error::TrailingDot { pos: 1 })]
    #[case::trailing_dot_2(b"C..", M4Error::ConsecutiveDot { pos: 1 })]
    #[case::double_dot(b"C..C", M4Error::ConsecutiveDot { pos: 1 })]
    #[case::dot_in_group_1(b"(.)", M4Error::LeadingDot { pos: 2 })]
    #[case::dot_in_group_2(b"(.)C", M4Error::EmptyGroup { pos: 2 })]
    #[case::dot_in_group_3(b"(.).C", M4Error::EmptyGroup { pos: 2 })]
    #[case::dot_in_branch_1(b"C(.)", M4Error::EmptyBranch { pos: 3 })]
    #[case::dot_in_branch_2(b"C(.)C", M4Error::EmptyBranch { pos: 3 })]
    #[case::dot_in_branch_3(b"C(.)(C)", M4Error::EmptyBranch { pos: 3 })]
    #[case::dot_in_component_1(b"().C", M4Error::EmptyGroup { pos: 1 })]
    #[case::dot_in_component_2(b"(.).C", M4Error::EmptyGroup { pos: 2 })]
    #[case::dot_in_component_3(b"(.).(C)", M4Error::EmptyGroup { pos: 2 })]
    #[case::dot_in_component_4(b"C.()", M4Error::EmptyGroup { pos: 3 })]
    #[case::dot_in_component_5(b"C.(.)", M4Error::LeadingDot { pos: 4 })]
    #[case::dot_in_component_6(b"(C).(.)", M4Error::LeadingDot { pos: 6 })]
    #[case::trailing_bond_before_dot_1(b"C-.C", M4Error::TrailingBond { pos: 1 })]
    #[case::trailing_bond_before_dot_2(b"C=.C", M4Error::TrailingBond { pos: 1 })]
    #[case::trailing_bond_before_dot_aromatic(b"C:.C", M4Error::TrailingBond { pos: 1 })]
    #[case::trailing_stereo_bond_before_dot_up(b"C/.C", M4Error::TrailingBond { pos: 1 })]
    #[case::trailing_stereo_bond_before_dot_down(b"C\\.C", M4Error::TrailingBond { pos: 1 })]
    #[case::dot_before_ring_digit(b"C.1", M4Error::LeadingRing { pos: 2 })]
    #[case::dot_before_ring_percent(b"C.%12", M4Error::LeadingRing { pos: 2 })]
    #[case::dot_before_branch_1(b"C1.(C)(C)C1", M4Error::TopLevelGroupTrailing { pos: 5 })]
    #[case::dot_unclosed_ring_1(b"C1.C", M4Error::RingUnclosed { open_pos: 1 })]
    #[case::dot_unclosed_ring_2(b"C.C1", M4Error::RingUnclosed { open_pos: 3 })]
    #[case::ring_order_conflict_digit(b"C=1.CC#1", M4Error::RingBondOrderConflict { pos: 7, open_pos: 2 })]
    #[case::ring_order_conflict_percent(b"C=%12.CC#%12", M4Error::RingBondOrderConflict { pos: 9, open_pos: 2 })]
    #[case::ring_dir_conflict_digit(b"C/1.CC\\1", M4Error::RingBondDirConflict { pos: 7, open_pos: 2 })]
    #[case::ring_dir_conflict_percent(b"C/%12.CC\\%12", M4Error::RingBondDirConflict { pos: 9, open_pos: 2 })]
    #[case::ring_dir_conflict_aromatic(b"c/1.cc\\1", M4Error::RingBondDirConflict { pos: 7, open_pos: 2 })]
    #[case::group_dot_before_ring_digit(b"(.1)", M4Error::LeadingRing { pos: 2 })]
    #[case::group_dot_before_ring_percent(b"(.%12)", M4Error::LeadingRing { pos: 2 })]
    #[case::branch_dot_before_ring_digit(b"C(.1)", M4Error::LeadingRing { pos: 3 })]
    #[case::branch_dot_before_ring_percent(b"C(.%12)", M4Error::LeadingRing { pos: 3 })]
    #[case::branch_dot_before_bond_1(b"(.-C)", M4Error::LeadingBond { pos: 2 })]
    #[case::branch_dot_before_bond_2(b"C(.-C)", M4Error::LeadingBond { pos: 3 })]
    #[case::leading_bond_after_dot_1(b"C.-C", M4Error::LeadingBond { pos: 2 })]
    #[case::leading_bond_after_dot_2(b"C.=-C", M4Error::LeadingBond { pos: 2 })]
    #[case::leading_sterebond_after_dot_up(b"C./C", M4Error::LeadingBond { pos: 2 })]
    #[case::leading_sterebond_after_dot_down(b"C.\\C", M4Error::LeadingBond { pos: 2 })]
    #[case::trailing_bond_dot_aromatic(b"C:.", M4Error::TrailingBond { pos: 1 })]
    #[case::group_trailing_bond_dot(b"(C-.)", M4Error::TrailingBond { pos: 2 })]
    #[case::branch_trailing_bond_dot(b"C(C-.)", M4Error::TrailingBond { pos: 3 })]
    #[case::group_before_atom_dot(b"C.(C)C", M4Error::TopLevelGroupTrailing { pos: 4 })]
    fn m4_components_invalid(#[case] input: &[u8], #[case] expected: M4Error) {
        let err = parse_smiles_m4(input);
        assert!(err.is_err(), "{:?} should have failed", input);
        let err = err.unwrap_err();
        assert_eq!(err, expected);
    }

    #[rstest]
    #[case::wildcard(b"*", M4Error::UnsupportedToken { pos: 0 })]
    #[case::bracket(b"[C]", M4Error::UnsupportedToken { pos: 0 })]
    #[case::whitespace_1(b"C ", M4Error::UnsupportedToken { pos: 1 })]
    #[case::whitespace_2(b"C\t", M4Error::UnsupportedToken { pos: 1 })]
    #[case::whitespace_3(b"C\n", M4Error::UnsupportedToken { pos: 1 })]
    #[case::whitespace_4(b"C\r", M4Error::UnsupportedToken { pos: 1 })]
    #[case::whitespace_5(b"C\r\n", M4Error::UnsupportedToken { pos: 1 })]
    fn m4_unimplemented(#[case] input: &[u8], #[case] expected: M4Error) {
        let err = parse_smiles_m4(input);
        assert!(err.is_err(), "{:?} should have failed", input);
        let err = err.unwrap_err();
        assert_eq!(err, expected);
    }
}
