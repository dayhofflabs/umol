use umol_data::Element;

use crate::io::ir::builder::{BondData, MoleculeBuilder};
use crate::io::ir::{BondDir, BondOrder, Molecule};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum M3Error {
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

pub fn parse_smiles_m3(input: &[u8]) -> Result<Molecule, M3Error> {
    let mut i = 0usize;
    let n = input.len();

    let mut builder = MoleculeBuilder::with_capacity(n.max(1), n.max(1).saturating_sub(1));
    let mut last_atom_idx: Option<u32> = None;
    let mut prev_atom_idx: Option<u32> = None; // immediate predecessor along current path
    let mut pending_bond: Option<(BondOrder, Option<BondDir>, usize)> = None;
    let mut last_aromatic: bool = false;

    let mut pstack: Vec<Frame> = Vec::new();
    let mut ring_table: Vec<Option<OpenRing>> = vec![None; 100]; // indices 0..99
    let mut just_closed_group: bool = false;

    while i < n {
        let b0 = input[i];
        if b0 != b'(' {
            just_closed_group = false;
        }

        // Parentheses
        if b0 == b'(' {
            if let Some((_, _, pos)) = pending_bond {
                return Err(M3Error::TrailingBond { pos });
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
                return Err(M3Error::TrailingBond { pos });
            }
            let Some(frame) = pstack.pop() else {
                return Err(M3Error::UnbalancedBranchClose { pos: i });
            };
            match frame {
                Frame::Branch { base, had_atom, .. } => {
                    if !had_atom {
                        return Err(M3Error::EmptyBranch { pos: i });
                    }
                    last_atom_idx = Some(base);
                    prev_atom_idx = None;
                }
                Frame::Group { had_atom, .. } => {
                    if !had_atom {
                        if i + 1 != n {
                            return Err(M3Error::EmptyGroup { pos: i });
                        }
                        // For empty top-level group at end-of-input keep no adjacency
                        last_atom_idx = None;
                        prev_atom_idx = None;
                        just_closed_group = false;
                    } else {
                        // Mark that we just closed a non-empty top-level group
                        just_closed_group = true;
                        // If this was the OUTERMOST group and there are trailing tokens, error
                        if pstack.is_empty() && i + 1 != n {
                            return Err(M3Error::TopLevelGroupTrailing { pos: i });
                        }
                    }
                    // If group had atoms, leave last_atom_idx as-is to allow adjacency
                }
            }
            i += 1;
            continue;
        }

        // Ring tokens: single digit 0..9 and %DD
        if is_digit(b0) {
            if last_atom_idx.is_none() {
                return Err(M3Error::LeadingRing { pos: i });
            }
            let idx: usize = (b0 - b'0') as usize;
            if idx >= ring_table.len() {
                // Should not happen, but guard
                return Err(M3Error::RingIndexInvalid { pos: i });
            }
            let bond = pending_bond.take();
            let (order_opt, dir_opt) = bond.map_or((None, None), |(o, d, _)| (Some(o), d));
            match ring_table[idx] {
                None => {
                    ring_table[idx] = Some(OpenRing {
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
                        return Err(M3Error::RingSelfLoop { pos: i });
                    }
                    if prev_atom_idx == Some(open.atom_id) {
                        return Err(M3Error::RingTwoMember { pos: i });
                    }
                    if let (Some(d1), Some(d2)) = (open.dir, dir_opt) {
                        if d1 != d2 {
                            return Err(M3Error::RingBondDirConflict {
                                pos: i,
                                open_pos: open.open_pos,
                            });
                        }
                    }
                    // Bond order conflict handling for ring closure
                    // If both sides specified an order and they differ -> conflict
                    if let (Some(o1), Some(o2)) = (open.order, order_opt) {
                        if o1 != o2 {
                            return Err(M3Error::RingBondOrderConflict {
                                pos: i,
                                open_pos: open.open_pos,
                            });
                        }
                    }
                    // Direction in combination with non-single order on the ring bond is unsupported
                    if open.dir.is_some() || dir_opt.is_some() {
                        let ord = open.order.or(order_opt).unwrap_or(BondOrder::Single);
                        if ord != BondOrder::Single {
                            return Err(M3Error::RingBondOrderConflict {
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
                    ring_table[idx] = None;
                }
            }
            i += 1;
            continue;
        }
        if b0 == b'%' {
            if i + 2 >= n || !is_digit(input[i + 1]) || !is_digit(input[i + 2]) {
                return Err(M3Error::RingIndexInvalid { pos: i });
            }
            if last_atom_idx.is_none() {
                return Err(M3Error::LeadingRing { pos: i });
            }
            let idx: usize = ((input[i + 1] - b'0') as usize) * 10 + (input[i + 2] - b'0') as usize;
            if idx >= ring_table.len() {
                return Err(M3Error::RingIndexInvalid { pos: i });
            }
            let bond = pending_bond.take();
            let (order_opt, dir_opt) = bond.map_or((None, None), |(o, d, _)| (Some(o), d));
            match ring_table[idx] {
                None => {
                    ring_table[idx] = Some(OpenRing {
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
                        return Err(M3Error::RingSelfLoop { pos: i });
                    }
                    if prev_atom_idx == Some(open.atom_id) {
                        return Err(M3Error::RingTwoMember { pos: i });
                    }
                    if let (Some(d1), Some(d2)) = (open.dir, dir_opt) {
                        if d1 != d2 {
                            return Err(M3Error::RingBondDirConflict {
                                pos: i,
                                open_pos: open.open_pos,
                            });
                        }
                    }
                    // Bond order conflict handling for ring closure (%NN form)
                    // If both sides specified an order and they differ -> conflict
                    if let (Some(o1), Some(o2)) = (open.order, order_opt) {
                        if o1 != o2 {
                            return Err(M3Error::RingBondOrderConflict {
                                pos: i,
                                open_pos: open.open_pos,
                            });
                        }
                    }
                    // Direction in combination with non-single order on the ring bond is unsupported
                    if open.dir.is_some() || dir_opt.is_some() {
                        let ord = open.order.or(order_opt).unwrap_or(BondOrder::Single);
                        if ord != BondOrder::Single {
                            return Err(M3Error::RingBondOrderConflict {
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
                    ring_table[idx] = None;
                }
            }
            i += 3; // %DD consumed; %DDD will naturally see D next
            continue;
        }

        // Bond tokens
        if matches!(b0, b'-' | b'=' | b'#' | b'$' | b':' | b'/' | b'\\') {
            if pending_bond.is_some() {
                return Err(M3Error::ConsecutiveBond { pos: i });
            }
            if last_atom_idx.is_none() {
                return Err(M3Error::LeadingBond { pos: i });
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

        return Err(M3Error::UnsupportedToken { pos: i });
    }

    if pending_bond.is_some() {
        let (_, _, pos) = pending_bond.unwrap();
        return Err(M3Error::TrailingBond { pos });
    }

    if !pstack.is_empty() {
        let pos = match pstack.last().unwrap() {
            Frame::Branch { open_pos, .. } | Frame::Group { open_pos, .. } => *open_pos,
        };
        return Err(M3Error::UnbalancedBranchOpen { pos });
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
        return Err(M3Error::RingUnclosed { open_pos: pos_open });
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
    fn m3_chain(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m3(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
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
    fn m3_tree(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m3(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ring_c_3(b"C1CC1", build_from_graph("C C C | 0-1 1-2 0-2"))]
    #[case::ring_c_10(b"C1CCCCCCCCC1", build_from_graph("C C C C C C C C C C | 0-1 1-2 2-3 3-4 4-5 5-6 6-7 7-8 8-9 0-9"))]
    #[case::ring_aromatic_c_6(b"c1ccccc1", build_from_graph("C* C* C* C* C* C* | 0-1: 1-2: 2-3: 3-4: 4-5: 0-5:"))]
    #[case::ring_index_0(b"C0CC0", build_from_graph("C C C | 0-1 1-2 0-2"))]
    #[case::ring_index_percent(b"C%12CC%12", build_from_graph("C C C | 0-1 1-2 0-2"))]
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
    fn m3_ring(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m3(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
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
    fn m3_bonds(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m3(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rstest]
    #[case::leading_ring(b"1C", M3Error::LeadingRing { pos: 0 })]
    #[case::bad_percent_short(b"C%1", M3Error::RingIndexInvalid { pos: 1 })]
    #[case::bad_percent_char(b"C%1a", M3Error::RingIndexInvalid { pos: 1 })]
    #[case::bad_percent_eoi(b"C%", M3Error::RingIndexInvalid { pos: 1 })]
    #[case::ring_self_loop(b"C11", M3Error::RingSelfLoop { pos: 2 })]
    #[case::ring_two_member(b"C1C1", M3Error::RingTwoMember { pos: 3 })]
    #[case::ring_bond_order_conflict_3(b"C=1CC#1", M3Error::RingBondOrderConflict { pos: 6, open_pos: 2 })]
    #[case::ring_bond_order_conflict_4(b"C/1CC=1", M3Error::RingBondOrderConflict { pos: 6, open_pos: 2 })]
    #[case::ring_bond_order_conflict_5(b"C\\1CC=1", M3Error::RingBondOrderConflict { pos: 6, open_pos: 2 })]
    #[case::ring_bond_order_conflict_6(b"C=1CC/1", M3Error::RingBondOrderConflict { pos: 6, open_pos: 2 })]
    #[case::ring_bond_order_conflict_7(b"C=1CC\\1", M3Error::RingBondOrderConflict { pos: 6, open_pos: 2 })]
    #[case::ring_bond_order_conflict_8(b"C=%10CC#%10", M3Error::RingBondOrderConflict { pos: 8, open_pos: 2 })]
    #[case::ring_bond_dir_conflict_1(b"C/1CC\\1", M3Error::RingBondDirConflict { pos: 6, open_pos: 2 })]
    #[case::ring_bond_dir_conflict_2(b"C\\1CC/1", M3Error::RingBondDirConflict { pos: 6, open_pos: 2 })]
    #[case::ring_bond_dir_conflict_3(b"C/%12CC\\%12", M3Error::RingBondDirConflict { pos: 8, open_pos: 2 })]
    #[case::ring_bond_dir_conflict_4(b"C\\%12CC/%12", M3Error::RingBondDirConflict { pos: 8, open_pos: 2 })]
    #[case::ring_unclosed_1(b"C1CC", M3Error::RingUnclosed { open_pos: 1 })]
    #[case::ring_unclosed_2(b"C1CC1C1", M3Error::RingUnclosed { open_pos: 6 })]
    #[case::component(b"CC.CC", M3Error::UnsupportedToken { pos: 2 })]
    #[case::unbalanced_closing_paren_1(b")C", M3Error::UnbalancedBranchClose { pos: 0 })]
    #[case::unbalanced_closing_paren_2(b"C)C", M3Error::UnbalancedBranchClose { pos: 1 })]
    #[case::unclosed_group(b"(C", M3Error::UnbalancedBranchOpen { pos: 0 })]
    #[case::unclosed_branch(b"C(C", M3Error::UnbalancedBranchOpen { pos: 1 })]
    #[case::empty_branch(b"C()", M3Error::EmptyBranch { pos: 2 })]
    #[case::empty_group_before_atom(b"()C", M3Error::EmptyGroup { pos: 1 })]
    #[case::two_top_level_groups(b"(C)(C)", M3Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::three_top_level_groups(b"(C)(C)(C)", M3Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::three_top_level_groups_aromatic(b"(c)(c)(c)", M3Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::two_top_level_groups_rings(b"(C1CC1)(C2CC2)", M3Error::TopLevelGroupTrailing { pos: 6 })]
    #[case::group_before_atom(b"(C)C", M3Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::group_before_atom_aromatic(b"(c)c", M3Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::trailing_bond_1(b"C-", M3Error::TrailingBond { pos: 1 })]
    #[case::trailing_bond_2(b"C=", M3Error::TrailingBond { pos: 1 })]
    #[case::trailing_bond_3(b"C#", M3Error::TrailingBond { pos: 1 })]
    #[case::trailing_bond_4(b"C$", M3Error::TrailingBond { pos: 1 })]
    #[case::trailing_stereo_bond_1(b"C/", M3Error::TrailingBond { pos: 1 })]
    #[case::trailing_stereo_bond_2(b"C\\", M3Error::TrailingBond { pos: 1 })]
    #[case::trailing_aromatic_bond(b"C:", M3Error::TrailingBond { pos: 1 })]
    #[case::branch_trailing_bond_1(b"C(C-)C", M3Error::TrailingBond { pos: 3 })]
    #[case::branch_trailing_bond_2(b"C(C=)C", M3Error::TrailingBond { pos: 3 })]
    #[case::branch_trailing_stereo_bond(b"CC(C/)CC", M3Error::TrailingBond { pos: 4 })]
    #[case::group_trailing_bond_1(b"(C-)", M3Error::TrailingBond { pos: 2 })]
    #[case::group_trailing_bond_2(b"(C=)", M3Error::TrailingBond { pos: 2 })]
    #[case::group_trailing_stereo_bond(b"(C/)", M3Error::TrailingBond { pos: 2 })]
    #[case::bond_after_group_1(b"(C)-", M3Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::bond_after_group_2(b"(C)=", M3Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::group_after_group_1(b"(C)(C)", M3Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::group_after_group_2(b"(c)(c)", M3Error::TopLevelGroupTrailing { pos: 2 })]
    #[case::ring_after_group(b"(C1CCC)1", M3Error::TopLevelGroupTrailing { pos : 6})]
    #[case::consecutive_bonds_1(b"C--C", M3Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bonds_2(b"C-=C", M3Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bonds_3(b"C-#C", M3Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bonds_4(b"C-$C", M3Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bonds_5(b"C-:C", M3Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_stereo_bonds_1(b"C//C", M3Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_stereo_bonds_2(b"C\\\\C", M3Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bond_and_stereo_bond_1(b"C-/C", M3Error::ConsecutiveBond { pos: 2 })]
    #[case::consecutive_bond_and_stereo_bond_2(b"C=\\C", M3Error::ConsecutiveBond { pos: 2 })]
    #[case::leading_bond_1(b"-C", M3Error::LeadingBond { pos: 0 })]
    #[case::leading_bond_2(b"=C", M3Error::LeadingBond { pos: 0 })]
    #[case::leading_bond_3(b"#C", M3Error::LeadingBond { pos: 0 })]
    #[case::leading_bond_4(b"$C", M3Error::LeadingBond { pos: 0 })]
    #[case::leading_aromatic_bond(b":C", M3Error::LeadingBond { pos: 0 })]
    #[case::leading_sterebond_1(b"/C", M3Error::LeadingBond { pos: 0 })]
    #[case::leading_sterebond_2(b"\\C", M3Error::LeadingBond { pos: 0 })]
    #[case::group_leading_bond_1(b"(-C)C", M3Error::LeadingBond { pos: 1 })]
    #[case::group_leading_bond_2(b"(=C)C", M3Error::LeadingBond { pos: 1 })]
    #[case::group_leading_bond_3(b"(#C)C", M3Error::LeadingBond { pos: 1 })]
    #[case::group_leading_bond_4(b"($C)C", M3Error::LeadingBond { pos: 1 })]
    #[case::group_leading_sterebond_1(b"(/C)C", M3Error::LeadingBond { pos: 1 })]
    #[case::group_leading_sterebond_2(b"(\\C)C", M3Error::LeadingBond { pos: 1 })]
    #[case::group_leading_aromatic_bond(b"(:C)C", M3Error::LeadingBond { pos: 1 })]
    fn m3_invalid(#[case] input: &[u8], #[case] expected: M3Error) {
        let err = parse_smiles_m3(input);
        assert!(err.is_err(), "{:?} should have failed", input);
        let err = err.unwrap_err();
        assert_eq!(err, expected);
    }
}
