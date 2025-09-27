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
    // M4 component-specific
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

// M4: components ('.') as a strict superset of M3
// Streaming FSM (copied from M3) with an added '.' branch that resets adjacency.
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
                Frame::Group { had_atom, .. } => {
                    if !had_atom {
                        if i + 1 != n {
                            return Err(M4Error::EmptyGroup { pos: i });
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
            if i == 0 { return Err(M4Error::LeadingDot { pos: i }); }
            if i + 1 == n { return Err(M4Error::TrailingDot { pos: i }); }
            if input[i + 1] == b'.' { return Err(M4Error::ConsecutiveDot { pos: i }); }
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
    #[case::two_components(b"CC.CC", build_from_graph("C C C C | 0-1 2-3"))]
    #[case::ring_across_dot_digit(b"C1.CC1", build_from_graph("C C C | 1-2 0-2"))]
    #[case::ring_across_dot_percent(b"C%12.CC%12", build_from_graph("C C C | 1-2 0-2"))]
    #[case::branch_local_components(b"C(C.C)", build_from_graph("C C C | 0-1"))]
    #[case::groups_and_dots(b"(CC).(CC)", build_from_graph("C C C C | 0-1 2-3"))]
    #[case::group_components_and_atom(b"(C.C).C", build_from_graph("C C C |"))]
    #[case::dots_around_groups_1(b"C.(C).C", build_from_graph("C C C |"))]
    #[case::dots_around_groups_2(b"C.C.(C)", build_from_graph("C C C |"))]
    #[case::rings_across_multiple_dots_digit(b"C1.C.CC1", build_from_graph("C C C C | 2-3 0-3"))]
    #[case::rings_across_multiple_dots_percent(b"C%12.C.CC%12", build_from_graph("C C C C | 2-3 0-3"))]
    fn m4_components_valid(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m4(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rstest]
    #[case::leading_dot(b".C")]
    #[case::trailing_dot(b"C.")]
    #[case::double_dot(b"C..C")]
    #[case::dot_before_ring_digit(b"C.1")]
    #[case::dot_before_ring_percent(b"C.%12")]
    fn m4_components_invalid(#[case] input: &[u8]) {
        let res = parse_smiles_m4(input);
        assert!(res.is_err(), "{:?} should have failed", input);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ring_double_unilateral_open(b"C=1.CC1", build_from_graph("C C C | 1-2 0-2:="))]
    #[case::ring_double_unilateral_close(b"C1.CC=1", build_from_graph("C C C | 1-2 0-2:="))]
    #[case::ring_dir_up_both(b"C/1.CC/1", build_from_graph("C C C | 1-2 0-2:/"))]
    #[case::ring_dir_down_both(b"C\\1.CC\\1", build_from_graph("C C C | 1-2 0-2:\\"))]
    #[case::ring_dir_up_both_percent(b"C/%12.CC/%12", build_from_graph("C C C | 1-2 0-2:/"))]
    #[case::ring_dir_down_both_percent(b"C\\%12.CC\\%12", build_from_graph("C C C | 1-2 0-2:\\"))]
    fn m4_rings_components_valid(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m4(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rstest]
    #[case::ring_order_conflict_digit(b"C=1.CC#1", M4Error::RingBondOrderConflict { pos: 7, open_pos: 2 })]
    #[case::ring_order_conflict_percent(b"C=%12.CC#%12", M4Error::RingBondOrderConflict { pos: 9, open_pos: 2 })]
    #[case::ring_dir_conflict_digit(b"C/1.CC\\1", M4Error::RingBondDirConflict { pos: 7, open_pos: 2 })]
    #[case::ring_dir_conflict_percent(b"C/%12.CC\\%12", M4Error::RingBondDirConflict { pos: 9, open_pos: 2 })]
    fn m4_rings_components_invalid(#[case] input: &[u8], #[case] expected: M4Error) {
        let err = parse_smiles_m4(input);
        assert!(err.is_err(), "{:?} should have failed", input);
        let err = err.unwrap_err();
        assert_eq!(err, expected);
    }

  
}
