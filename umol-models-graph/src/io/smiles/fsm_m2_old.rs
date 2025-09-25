pub use crate::io::smiles::fsm_m2_bonds::*;

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[derive(Debug, Clone, Copy)]
    enum Connection {
        Bonded(usize),  // number of linear atoms between rings (0 => direct bond)
        Bridged(usize), // number of shared atoms (>=1)
    }

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

    fn build_ring_c(n: usize) -> Molecule {
        if n == 0 {
            return Molecule::default();
        }

        let mut b = MoleculeBuilder::with_capacity(n, n);

        // Create first atom
        let first = b.on_atom_fast(Element::C, true, false);
        let mut prev = first;

        // Create remaining atoms and linear bonds
        for _ in 1..n {
            let curr = b.on_atom_fast(Element::C, true, false);
            b.on_bond_single_fast(prev, curr);
            prev = curr;
        }

        // Close the ring to the first atom (match parser orientation: first -> last)
        b.on_bond_single_fast(first, prev);

        let mut mols = b.finish();
        mols.pop().unwrap_or_default()
    }

    fn build_two_rings_c(n_ring1: usize, n_ring2: usize, connection: Connection) -> Molecule {
        match connection {
            Connection::Bonded(n_between) => {
                let total = n_ring1 + n_between + n_ring2;
                let mut b = MoleculeBuilder::with_capacity(total, total + 2);

                // Ring 1
                let r1_first = b.on_atom_fast(Element::C, true, false);
                let mut prev = r1_first;
                for _ in 1..n_ring1 {
                    let curr = b.on_atom_fast(Element::C, true, false);
                    b.on_bond_single_fast(prev, curr);
                    prev = curr;
                }
                b.on_bond_single_fast(r1_first, prev);
                let r1_last = prev;

                // Connector chain (possibly zero length)
                let mut tail = r1_last;
                for _ in 0..n_between {
                    let curr = b.on_atom_fast(Element::C, true, false);
                    b.on_bond_single_fast(tail, curr);
                    tail = curr;
                }

                // Ring 2, anchor at tail
                let r2_first = b.on_atom_fast(Element::C, true, false);
                b.on_bond_single_fast(tail, r2_first);
                let mut prev2 = r2_first;
                for _ in 1..n_ring2 {
                    let curr = b.on_atom_fast(Element::C, true, false);
                    b.on_bond_single_fast(prev2, curr);
                    prev2 = curr;
                }
                b.on_bond_single_fast(r2_first, prev2);

                let mut mols = b.finish();
                mols.pop().unwrap_or_default()
            }
            Connection::Bridged(shared) => {
                assert!(shared >= 1, "shared atoms must be >= 1");
                let add1 = n_ring1.saturating_sub(shared);
                let add2 = n_ring2.saturating_sub(shared);
                let total = shared + add1 + add2;
                let mut b = MoleculeBuilder::with_capacity(total, total + 2);

                // Shared path S0..S{shared-1}
                let s_first = b.on_atom_fast(Element::C, true, false);
                let mut prev = s_first;
                for _ in 1..shared {
                    let curr = b.on_atom_fast(Element::C, true, false);
                    b.on_bond_single_fast(prev, curr);
                    prev = curr;
                }
                let s_last = prev;

                // Ring 1 extension off s_last, then close to s_first
                let mut tail1 = s_last;
                for _ in 0..add1 {
                    let curr = b.on_atom_fast(Element::C, true, false);
                    b.on_bond_single_fast(tail1, curr);
                    tail1 = curr;
                }
                b.on_bond_single_fast(s_first, tail1);

                // Ring 2 extension also off s_last, then close to s_first (mirror orientation to match parser DFS)
                let mut tail2 = s_last;
                for _ in 0..add2 {
                    let curr = b.on_atom_fast(Element::C, true, false);
                    b.on_bond_single_fast(tail2, curr);
                    tail2 = curr;
                }
                b.on_bond_single_fast(s_first, tail2);

                let mut mols = b.finish();
                mols.pop().unwrap_or_default()
            }
        }
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

    #[rstest]
    #[case::ring3(b"C1CC1", build_ring_c(3))]
    #[case::ring10(b"C1CCCCCCCCC1", build_ring_c(10))]
    #[case::ring_index0(b"C0CC0", build_ring_c(3))]
    #[case::ring_percent(b"C%12CC%12", build_ring_c(3))]
    #[case::two_rings_bonded(b"C1CC1C2CC2", build_two_rings_c(3, 3, Connection::Bonded(0)))]
    #[case::two_rings_bonded_2(b"C1CC1CCC2CC2", build_two_rings_c(3, 3, Connection::Bonded(2)))]
    #[case::two_rings_spiro(b"C1CC12CC2", build_two_rings_c(3, 3, Connection::Bridged(1)))]
    #[case::two_rings_fused(b"C", build_two_rings_c(3, 3, Connection::Bridged(2)))]
    #[case::two_rings_bridged(b"C(CC(C1)(C2))12", build_two_rings_c(4, 4, Connection::Bridged(3)))]
    fn m2_ring(#[case] input: &[u8], #[case] expected: Molecule) {
        let res = parse_smiles_m2(input);
        assert!(res.is_ok(), "{:?} should have succeeded", input);
        let mol = res.unwrap();
        assert_eq!(mol, expected);
    }

    #[rstest]
    #[case::bond_order(b"C-C", M2Error::UnsupportedToken { pos: 1 })]
    #[case::bracket(b"[C]", M2Error::UnsupportedToken { pos: 0 })]
    #[case::component(b"CC.CC", M2Error::UnsupportedToken { pos: 2 })]
    #[case::stray_closing_paren(b")C", M2Error::UnbalancedBranchClose { pos: 0 })]
    #[case::unclosed_group(b"(C", M2Error::UnbalancedBranchOpen { pos: 0 })]
    #[case::unclosed_branch(b"C(C", M2Error::UnbalancedBranchOpen { pos: 1 })]
    #[case::trailing_bond_in_branch(b"C(C-)C", M2Error::UnsupportedToken { pos: 3 })]
    #[case::empty_branch(b"C()", M2Error::EmptyBranch { pos: 2 })]
    #[case::empty_group_before_atom(b"()C", M2Error::EmptyGroup { pos: 1 })]
    fn m2_invalid(#[case] input: &[u8], #[case] expected: M2Error) {
        let err = parse_smiles_m2(input);
        assert!(err.is_err(), "{:?} should have failed", input);
        let err = err.unwrap_err();
        assert_eq!(err, expected);
    }
}
