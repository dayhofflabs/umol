use pretty_assertions::assert_eq;
use rstest::*;
use umol_data::Element;

use self::utils::build_from_graph;
use super::*;
use crate::io::ir::{AtomSymbol, BondDir, BondSymbol, Chirality};

#[rstest]
#[case::organic_c(b"C", build_from_graph("C |"))]
#[case::organic_b(b"B", build_from_graph("B |"))]
#[case::organic_n(b"N", build_from_graph("N |"))]
#[case::organic_o(b"O", build_from_graph("O |"))]
#[case::organic_s(b"S", build_from_graph("S |"))]
#[case::organic_p(b"P", build_from_graph("P |"))]
#[case::organic_f(b"F", build_from_graph("F |"))]
#[case::organic_cl(b"Cl", build_from_graph("Cl |"))]
#[case::organic_br(b"Br", build_from_graph("Br |"))]
#[case::organic_i(b"I", build_from_graph("I |"))]
#[case::organic_b_aromatic(b"b", build_from_graph("B* |"))]
#[case::organic_c_aromatic(b"c", build_from_graph("C* |"))]
#[case::organic_n_aromatic(b"n", build_from_graph("N* |"))]
#[case::organic_o_aromatic(b"o", build_from_graph("O* |"))]
#[case::organic_s_aromatic(b"s", build_from_graph("S* |"))]
#[case::organic_p_aromatic(b"p", build_from_graph("P* |"))]
fn element(#[case] input: &[u8], #[case] expected: Molecule) {
    let res = parse_smiles(input);
    assert!(res.is_ok(), "{:?} should have succeeded", input);
    let mol = res.unwrap();
    assert_eq!(mol, expected);
}

#[rstest]
// TODO: Should be InvalidElement
#[case::element_h(b"H", ParseError::UnsupportedToken { pos: 0 })]
#[case::element_he(b"He", ParseError::UnsupportedToken { pos: 0 })]
#[case::element_q(b"Q", ParseError::UnsupportedToken { pos: 0 })]
#[case::element_f_aromatic(b"f", ParseError::UnsupportedToken { pos: 0 })]
fn element_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let res = parse_smiles(input);
    assert!(res.is_err(), "{:?} should have failed", input);
    let mol = res.unwrap_err();
    assert_eq!(mol, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::empty(b"", Molecule::default())]
#[case::chain_c_1(b"C", build_from_graph("C |"))]
#[case::chain_c_5(b"CCCCC", build_from_graph("C C C C C | 0-1 1-2 2-3 3-4"))]
#[case::aromatic_c_6(b"cccccc", build_from_graph("C* C* C* C* C* C* | 0-1: 1-2: 2-3: 3-4: 4-5:"))]
#[case::chain_mixed_5(b"CClOBrN", build_from_graph("C Cl O Br N | 0-1 1-2 2-3 3-4"))]
fn chain(#[case] input: &[u8], #[case] expected: Molecule) {
    let res = parse_smiles(input);
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
#[case::branch_c_211(b"CC(C)C", build_from_graph("C C C C | 0-1 1-2 1-3"))]
#[case::branch_c_222_aromatic(b"cc(cc)cc", build_from_graph("C* C* C* C* C* C* | 0-1: 1-2: 2-3: 1-4: 4-5:"))]
#[case::branch_trailing(b"C(CC)", build_from_graph("C C C | 0-1 1-2"))]
#[case::branch_multiple(b"CC(C)(C)C", build_from_graph("C C C C C | 0-1 1-2 1-3 1-4"))]
#[case::branch_multiple_trailing(b"C(C)(C)", build_from_graph("C C C | 0-1 0-2"))]
#[case::branch_nested(b"C(C(C)C)C", build_from_graph("C C C C C | 0-1 1-2 1-3 0-4"))]
fn tree(#[case] input: &[u8], #[case] expected: Molecule) {
    let res = parse_smiles(input);
    assert!(res.is_ok(), "{:?} should have succeeded", input);
    let mol = res.unwrap();
    assert_eq!(mol, expected);
}

#[rstest]
#[case::unbalanced_closing_paren_1(b")C", ParseError::UnbalancedBranchClose { pos: 0 })]
#[case::unbalanced_closing_paren_2(b"C)C", ParseError::UnbalancedBranchClose { pos: 1 })]
#[case::unclosed_group(b"(C", ParseError::UnbalancedBranchOpen { pos: 0 })]
#[case::unclosed_branch(b"C(C", ParseError::UnbalancedBranchOpen { pos: 1 })]
#[case::empty_branch(b"C()", ParseError::EmptyBranch { pos: 2 })]
#[case::empty_group_before_atom(b"()C", ParseError::EmptyGroup { pos: 1 })]
#[case::two_top_level_groups(b"(C)(C)", ParseError::TopLevelGroupTrailing { pos: 2 })]
#[case::three_top_level_groups(b"(C)(C)(C)", ParseError::TopLevelGroupTrailing { pos: 2 })]
#[case::three_top_level_groups_aromatic(b"(c)(c)(c)", ParseError::TopLevelGroupTrailing { pos: 2 })]
#[case::two_top_level_groups_rings(b"(C1CC1)(C2CC2)", ParseError::TopLevelGroupTrailing { pos: 6 })]
#[case::group_before_atom(b"(C)C", ParseError::TopLevelGroupTrailing { pos: 2 })]
#[case::group_before_atom_aromatic(b"(c)c", ParseError::TopLevelGroupTrailing { pos: 2 })]
fn tree_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let err = parse_smiles(input);
    assert!(err.is_err(), "{:?} should have failed", input);
    let err = err.unwrap_err();
    assert_eq!(err, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::ring_c_3(b"C1CC1", build_from_graph("C C C | 0-1 1-2 0-2"))]
#[case::ring_c_6(b"C1CCCCC1", build_from_graph("C C C C C C | 0-1 1-2 2-3 3-4 4-5 0-5"))]
#[case::ring_c_10(b"C1CCCCCCCCC1", build_from_graph("C C C C C C C C C C | 0-1 1-2 2-3 3-4 4-5 5-6 6-7 7-8 8-9 0-9"))]
#[case::ring_aromatic_c_6(b"c1ccccc1", build_from_graph("C* C* C* C* C* C* | 0-1: 1-2: 2-3: 3-4: 4-5: 0-5:"))]
#[case::ring_index_0(b"C0CC0", build_from_graph("C C C | 0-1 1-2 0-2"))]
#[case::ring_index_percent(b"C%12CC%12", build_from_graph("C C C | 0-1 1-2 0-2"))]
#[case::ring_index_zero_prefix_1(b"C1CC%01", build_from_graph("C C C | 0-1 1-2 0-2"))]
#[case::ring_index_zero_prefix_2(b"C0CC%00", build_from_graph("C C C | 0-1 1-2 0-2"))]
#[case::ring_index_zero_prefix_3(b"C9CC%09", build_from_graph("C C C | 0-1 1-2 0-2"))]
#[case::ring_index_max_99(b"C%99CC%99", build_from_graph("C C C | 0-1 1-2 0-2"))]
#[case::ring_indices_single_percent(b"C%123CCC%12CC3", build_from_graph("C C C C C C | 0-1 1-2 2-3 0-3 3-4 4-5 0-5"))]
#[case::two_rings_bonded_0(b"C1CC1C2CC2", build_from_graph("C C C C C C | 0-1 1-2 0-2 2-3 3-4 4-5 3-5"))]
#[case::two_rings_bonded_0_aromatic_1(b"c1cc1c2cc2", build_from_graph("C* C* C* C* C* C* | 0-1: 1-2: 0-2: 2-3: 3-4: 4-5: 3-5:"))]
#[case::two_rings_bonded_0_aromatic_2(b"c1cc1C2CC2", build_from_graph("C* C* C* C C C | 0-1: 1-2: 0-2: 2-3 3-4 4-5 3-5"))]
#[case::two_rings_index_reused(b"C1CC1C1CC1", build_from_graph("C C C C C C | 0-1 1-2 0-2 2-3 3-4 4-5 3-5"))]
#[case::two_rings_bonded_2(b"C1CC1CCC2CC2", build_from_graph("C C C C C C C C | 0-1 1-2 0-2 2-3 3-4 4-5 5-6 6-7 5-7"))]
#[case::two_rings_spiro(b"C1CC12CC2", build_from_graph("C C C C C | 0-1 1-2 0-2 2-3 3-4 2-4"))]
#[case::two_rings_fused(b"C12CC1C2", build_from_graph("C C C C | 0-1 1-2 0-2 2-3 0-3"))]
#[case::two_rings_bridged(b"C12CC(C2)C1", build_from_graph("C C C C C | 0-1 1-2 2-3 0-3 2-4 0-4"))]
#[case::two_rings_fused_aromatic(b"c12ccccc1cccc2", build_from_graph("C* C* C* C* C* C* C* C* C* C* | 0-1: 1-2: 2-3: 3-4: 4-5: 0-5: 5-6: 6-7: 7-8: 8-9: 0-9:"))]
#[case::two_rings_interleaved_indices(b"N1CC2CCCCC2CC1", build_from_graph("N C C C C C C C C C | 0-1 1-2 2-3 3-4 4-5 5-6 6-7 2-7 7-8 8-9 0-9"))]
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
#[case::substituted_ring_branch(b"C1C(C(C)C)C1", build_from_graph("C C C C C C | 0-1 1-2 2-3 2-4 1-5 0-5"))]
fn ring(#[case] input: &[u8], #[case] expected: Molecule) {
    let res = parse_smiles(input);
    assert!(res.is_ok(), "{:?} should have succeeded", input);
    let mol = res.unwrap();
    assert_eq!(mol, expected);
}

#[rstest]
#[case::leading_ring_0(b"0C", ParseError::LeadingRing { pos: 0 })]
#[case::leading_ring_1(b"1C", ParseError::LeadingRing { pos: 0 })]
#[case::leading_ring_percent(b"%12C", ParseError::LeadingRing { pos: 0 })]
#[case::leading_ring_group(b"(1CCC)", ParseError::LeadingRing { pos: 1 })]
#[case::leading_ring_branch(b"C(1CCC)", ParseError::LeadingRing { pos: 0 })]
#[case::ring_unclosed_1(b"C1CC", ParseError::RingUnclosed { open_pos: 1 })]
#[case::ring_unclosed_2(b"C1CC1C1", ParseError::RingUnclosed { open_pos: 6 })]
#[case::ring_unclosed_3(b"C1CC2C", ParseError::RingUnclosed { open_pos: 4 })]
#[case::ring_unclosed_self_loop(b"C111", ParseError::RingUnclosed { open_pos: 3 })]
#[case::ring_unclosed_percent(b"C%12CC", ParseError::RingUnclosed { open_pos: 1 })]
#[case::bad_percent_no_index_1(b"C%", ParseError::RingIndexInvalid { pos: 1 })]
#[case::bad_percent_no_index_2(b"C%C", ParseError::RingIndexInvalid { pos: 1 })]
#[case::bad_percent_single_digit_0(b"C%0", ParseError::RingIndexInvalid { pos: 1 })]
#[case::bad_percent_single_digit_1(b"C%1", ParseError::RingIndexInvalid { pos: 1 })]
#[case::bad_percent_char(b"C%1a", ParseError::RingIndexInvalid { pos: 1 })]
fn ring_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let err = parse_smiles(input);
    assert!(err.is_err(), "{:?} should have failed", input);
    let err = err.unwrap_err();
    assert_eq!(err, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::ring_self_loop(b"C11", build_from_graph("C | 0-0"))]
#[case::ring_self_loop_percent(b"C%11%11", build_from_graph("C | 0-0"))]
#[case::ring_two_member(b"C1C1", build_from_graph("C C | 0-1 0-1"))]
#[case::ring_two_member_multiple(b"C12C12", build_from_graph("C C | 0-1 0-1 0-1"))]
#[case::ring_two_member_percent(b"C%12C%12", build_from_graph("C C | 0-1 0-1"))]
#[case::ring_two_member_single_percent(b"C%123CCC%123", build_from_graph("C C C C | 0-1 1-2 2-3 0-3 0-3"))]
#[case::ring_multiple_rings(b"C12CCCCC12", build_from_graph("C C C C C C | 0-1 1-2 2-3 3-4 4-5 0-5 0-5"))]
#[case::ring_multiple_rings_triple(b"C123CCCCC123", build_from_graph("C C C C C C | 0-1 1-2 2-3 3-4 4-5 0-5 0-5 0-5"))]
#[case::ring_multiple_rings_percent(b"C%12%13CCCCC%12%13", build_from_graph("C C C C C C | 0-1 1-2 2-3 3-4 4-5 0-5 0-5"))]
fn ring_invalid_topology(#[case] input: &[u8], #[case] expected: Molecule) {
    // Expected to pass here, fail in post-parse topology check
    let res = parse_smiles(input);
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
#[case::branch_leading_single_bond(b"CC(-C)C", build_from_graph("C C C C | 0-1 1-2 1-3"))]
#[case::branch_leading_single_bond_multiple(b"CC(-C)(-C)C", build_from_graph("C C C C C | 0-1 1-2 1-3 1-4"))]
#[case::branch_leading_double_bond(b"CC(=C)C", build_from_graph("C C C C | 0-1 1-2:= 1-3"))]
#[case::branch_leading_double_bond_multiple(b"OS(=O)(=O)O", build_from_graph("O S O O O | 0-1 1-2:= 1-3:= 1-4"))]
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
#[case::ring_single_bond_percent(b"C-%12-C-C-%12", build_from_graph("C C C | 0-1 1-2 0-2"))]
#[case::ring_double_bond_1(b"C1-C=C1", build_from_graph("C C C | 0-1 1-2:= 0-2"))]
#[case::ring_double_bond_2(b"C1-CC=1", build_from_graph("C C C | 0-1 1-2 0-2:="))]
#[case::ring_double_bond_3(b"C=1-CC1", build_from_graph("C C C | 0-1 1-2 0-2:="))]
#[case::ring_double_bond_4(b"C=1-C-C=1", build_from_graph("C C C | 0-1 1-2 0-2:="))]
#[case::ring_double_bond_5(b"C=1CCCCC=1", build_from_graph("C C C C C C | 0-1 1-2 2-3 3-4 4-5 0-5:="))]
#[case::ring_double_bond_unilateral_close_1(b"C1CC=1", build_from_graph("C C C | 0-1 1-2 0-2:="))]
#[case::ring_double_bond_unilateral_close_2(b"C1CCCCC=1", build_from_graph("C C C C C C | 0-1 1-2 2-3 3-4 4-5 0-5:="))]
#[case::ring_double_bond_unilateral_open_1(b"C=1CC1", build_from_graph("C C C | 0-1 1-2 0-2:="))]
#[case::ring_double_bond_unilateral_open_2(b"C=1CCCCC1", build_from_graph("C C C C C C | 0-1 1-2 2-3 3-4 4-5 0-5:="))]
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
fn bonds(#[case] input: &[u8], #[case] expected: Molecule) {
    let res = parse_smiles(input);
    assert!(res.is_ok(), "{:?} should have succeeded", input);
    let mol = res.unwrap();
    assert_eq!(mol, expected);
}

#[rstest]
#[case::trailing_bond_1(b"C-", ParseError::TrailingBond { pos: 1 })]
#[case::trailing_bond_2(b"C=", ParseError::TrailingBond { pos: 1 })]
#[case::trailing_bond_3(b"C#", ParseError::TrailingBond { pos: 1 })]
#[case::trailing_bond_4(b"C$", ParseError::TrailingBond { pos: 1 })]
#[case::trailing_stereo_bond_1(b"C/", ParseError::TrailingBond { pos: 1 })]
#[case::trailing_stereo_bond_2(b"C\\", ParseError::TrailingBond { pos: 1 })]
#[case::trailing_aromatic_bond(b"C:", ParseError::TrailingBond { pos: 1 })]
#[case::branch_trailing_bond_1(b"C(C-)C", ParseError::TrailingBond { pos: 3 })]
#[case::branch_trailing_bond_2(b"C(C=)C", ParseError::TrailingBond { pos: 3 })]
#[case::branch_trailing_stereo_bond(b"CC(C/)CC", ParseError::TrailingBond { pos: 4 })]
#[case::group_trailing_bond_1(b"(C-)", ParseError::TrailingBond { pos: 2 })]
#[case::group_trailing_bond_2(b"(C=)", ParseError::TrailingBond { pos: 2 })]
#[case::group_trailing_stereo_bond(b"(C/)", ParseError::TrailingBond { pos: 2 })]
#[case::group_trailing_aromatic_bond(b"(C:)", ParseError::TrailingBond { pos: 2 })]
#[case::trailing_bond_before_dot_1(b"C-.C", ParseError::TrailingBond { pos: 1 })]
#[case::trailing_bond_before_dot_2(b"C=.C", ParseError::TrailingBond { pos: 1 })]
#[case::trailing_bond_before_dot_aromatic(b"C:.C", ParseError::TrailingBond { pos: 1 })]
#[case::trailing_stereo_bond_before_dot_up(b"C/.C", ParseError::TrailingBond { pos: 1 })]
#[case::trailing_stereo_bond_before_dot_down(b"C\\.C", ParseError::TrailingBond { pos: 1 })]
#[case::bond_after_group_1(b"(C)-", ParseError::TopLevelGroupTrailing { pos: 2 })]
#[case::bond_after_group_2(b"(C)=", ParseError::TopLevelGroupTrailing { pos: 2 })]
#[case::group_after_group_1(b"(C)(C)", ParseError::TopLevelGroupTrailing { pos: 2 })]
#[case::group_after_group_2(b"(c)(c)", ParseError::TopLevelGroupTrailing { pos: 2 })]
#[case::ring_after_group(b"(C1CCC)1", ParseError::TopLevelGroupTrailing { pos : 6})]
#[case::consecutive_bonds_1(b"C--C", ParseError::ConsecutiveBond { pos: 2 })]
#[case::consecutive_bonds_2(b"C-=C", ParseError::ConsecutiveBond { pos: 2 })]
#[case::consecutive_bonds_3(b"C-#C", ParseError::ConsecutiveBond { pos: 2 })]
#[case::consecutive_bonds_4(b"C-$C", ParseError::ConsecutiveBond { pos: 2 })]
#[case::consecutive_bonds_5(b"C-:C", ParseError::ConsecutiveBond { pos: 2 })]
#[case::consecutive_stereo_bonds_1(b"C//C", ParseError::ConsecutiveBond { pos: 2 })]
#[case::consecutive_stereo_bonds_2(b"C\\\\C", ParseError::ConsecutiveBond { pos: 2 })]
#[case::consecutive_bond_and_stereo_bond_1(b"C-/C", ParseError::ConsecutiveBond { pos: 2 })]
#[case::consecutive_bond_and_stereo_bond_2(b"C=\\C", ParseError::ConsecutiveBond { pos: 2 })]
#[case::leading_bond_1(b"-C", ParseError::LeadingBond { pos: 0 })]
#[case::leading_bond_2(b"=C", ParseError::LeadingBond { pos: 0 })]
#[case::leading_bond_3(b"#C", ParseError::LeadingBond { pos: 0 })]
#[case::leading_bond_4(b"$C", ParseError::LeadingBond { pos: 0 })]
#[case::leading_aromatic_bond(b":C", ParseError::LeadingBond { pos: 0 })]
#[case::leading_sterebond_1(b"/C", ParseError::LeadingBond { pos: 0 })]
#[case::leading_sterebond_2(b"\\C", ParseError::LeadingBond { pos: 0 })]
#[case::group_leading_bond_1(b"(-C)C", ParseError::GroupLeadingConnector { pos: 1 })]
#[case::group_leading_bond_2(b"(=C)C", ParseError::GroupLeadingConnector { pos: 1 })]
#[case::group_leading_bond_3(b"(#C)C", ParseError::GroupLeadingConnector { pos: 1 })]
#[case::group_leading_bond_4(b"($C)C", ParseError::GroupLeadingConnector { pos: 1 })]
#[case::group_leading_sterebond_1(b"(/C)C", ParseError::GroupLeadingConnector { pos: 1 })]
#[case::group_leading_sterebond_2(b"(\\C)C", ParseError::GroupLeadingConnector { pos: 1 })]
#[case::group_leading_aromatic_bond(b"(:C)C", ParseError::GroupLeadingConnector { pos: 1 })]
#[case::ring_bond_order_conflict_1(b"C-1CCCCC=1", ParseError::RingBondOrderConflict { pos: 9, open_pos: 2 })]
#[case::ring_bond_order_conflict_2(b"C=1CCCCC-1", ParseError::RingBondOrderConflict { pos: 9, open_pos: 2 })]
#[case::ring_bond_order_conflict_3(b"C=1CC#1", ParseError::RingBondOrderConflict { pos: 6, open_pos: 2 })]
#[case::ring_bond_order_conflict_4(b"C/1CC=1", ParseError::RingBondOrderConflict { pos: 6, open_pos: 2 })]
#[case::ring_bond_order_conflict_5(b"C\\1CC=1", ParseError::RingBondOrderConflict { pos: 6, open_pos: 2 })]
#[case::ring_bond_order_conflict_6(b"C=1CC/1", ParseError::RingBondOrderConflict { pos: 6, open_pos: 2 })]
#[case::ring_bond_order_conflict_7(b"C=1CC\\1", ParseError::RingBondOrderConflict { pos: 6, open_pos: 2 })]
#[case::ring_bond_order_conflict_8(b"C=%10CC#%10", ParseError::RingBondOrderConflict { pos: 8, open_pos: 2 })]
#[case::ring_bond_dir_conflict_1(b"C/1CC\\1", ParseError::RingBondDirConflict { pos: 6, open_pos: 2 })]
#[case::ring_bond_dir_conflict_2(b"C\\1CC/1", ParseError::RingBondDirConflict { pos: 6, open_pos: 2 })]
#[case::ring_bond_dir_conflict_3(b"C/%12CC\\%12", ParseError::RingBondDirConflict { pos: 8, open_pos: 2 })]
#[case::ring_bond_dir_conflict_4(b"C\\%12CC/%12", ParseError::RingBondDirConflict { pos: 8, open_pos: 2 })]
fn bonds_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let err = parse_smiles(input);
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
#[case::leading_dot_in_branch_3(b"C(.C.C)", build_from_graph("C C C |"))]
#[case::leading_dot_in_branch_4(b"C(C)(.C)", build_from_graph("C C C | 0-1"))]
#[case::trailing_dot_in_branch_1(b"C(C.)", build_from_graph("C C | 0-1"))]
#[case::trailing_dot_in_branch_2(b"C(C.)C", build_from_graph("C C C | 0-1 0-2"))]
#[case::trailing_dot_in_branch_3(b"C(C.)(C)", build_from_graph("C C C | 0-1 0-2"))]
#[case::group_components_1(b"(C.CC.C)", build_from_graph("C C C C | 1-2"))]
#[case::group_components_2(b"(CC).(CC)", build_from_graph("C C C C | 0-1 2-3"))]
#[case::group_components_3(b"(C.C).C", build_from_graph("C C C |"))]
#[case::group_components_4(b"C.(C).C", build_from_graph("C C C |"))]
#[case::group_components_5(b"C.C.(C)", build_from_graph("C C C |"))]
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
fn components(#[case] input: &[u8], #[case] expected: Molecule) {
    let res = parse_smiles(input);
    assert!(res.is_ok(), "{:?} should have succeeded", input);
    let mol = res.unwrap();
    assert_eq!(mol, expected);
}

#[rstest]
#[case::leading_dot_1(b".", ParseError::LeadingDot { pos: 0 })]
#[case::leading_dot_2(b".C", ParseError::LeadingDot { pos: 0 })]
#[case::leading_dot_3(b"..C", ParseError::LeadingDot { pos: 0 })]
#[case::leading_dot_4(b".C.", ParseError::LeadingDot { pos: 0 })]
#[case::trailing_dot_1(b"C.", ParseError::TrailingDot { pos: 1 })]
#[case::trailing_dot_2(b"C..", ParseError::ConsecutiveDot { pos: 1 })]
#[case::double_dot(b"C..C", ParseError::ConsecutiveDot { pos: 1 })]
#[case::dot_before_ring_digit(b"C.1", ParseError::LeadingRing { pos: 2 })]
#[case::dot_before_ring_percent(b"C.%12", ParseError::LeadingRing { pos: 2 })]
// TODO: Replace GroupLeadingConnector by LeadingDot / LeadingBond?
#[case::dot_in_group_1(b"(.)", ParseError::GroupLeadingConnector { pos: 1 })]
#[case::dot_in_group_2(b"(.)C", ParseError::GroupLeadingConnector { pos: 1 })]
#[case::dot_in_group_3(b"(.).C", ParseError::GroupLeadingConnector { pos: 1 })]
#[case::leading_dot_in_group_1(b"(.CC)", ParseError::GroupLeadingConnector { pos: 1 })]
#[case::leading_dot_in_group_2(b"(.CC).(CC)", ParseError::GroupLeadingConnector { pos: 1 })]
#[case::leading_dot_in_group_3(b"(CC).(.CC)", ParseError::GroupLeadingConnector { pos: 6 })]
#[case::leading_dot_in_group_4(b"(.C).(.C)", ParseError::GroupLeadingConnector { pos: 1 })]
#[case::leading_dot_in_group_5(b"C.(.C).C", ParseError::GroupLeadingConnector { pos: 3 })]
#[case::dot_before_group(b"C.(C)C", ParseError::TopLevelGroupTrailing { pos: 4 })]
#[case::dot_in_branch_1(b"C(.)", ParseError::EmptyBranch { pos: 3 })]
#[case::dot_in_branch_2(b"C(.)C", ParseError::EmptyBranch { pos: 3 })]
#[case::dot_in_branch_3(b"C(.)(C)", ParseError::EmptyBranch { pos: 3 })]
#[case::dot_in_component_1(b"().C", ParseError::EmptyGroup { pos: 1 })]
#[case::dot_in_component_2(b"(.).C", ParseError::GroupLeadingConnector { pos: 1 })]
#[case::dot_in_component_3(b"(.).(C)", ParseError::GroupLeadingConnector { pos: 1})]
#[case::dot_in_component_4(b"C.()", ParseError::EmptyGroup { pos: 3 })]
#[case::dot_in_component_5(b"C.(.)", ParseError::GroupLeadingConnector { pos: 3 })]
#[case::dot_in_component_6(b"(C).(.)", ParseError::GroupLeadingConnector { pos: 5 })]
#[case::dot_unclosed_ring_1(b"C1.C", ParseError::RingUnclosed { open_pos: 1 })]
#[case::dot_unclosed_ring_2(b"C.C1", ParseError::RingUnclosed { open_pos: 3 })]
#[case::dot_unclosed_ring_before_group(b"C1.(C)(C)C1", ParseError::TopLevelGroupTrailing { pos: 5 })]
#[case::ring_order_conflict_digit(b"C=1.CC#1", ParseError::RingBondOrderConflict { pos: 7, open_pos: 2 })]
#[case::ring_order_conflict_percent(b"C=%12.CC#%12", ParseError::RingBondOrderConflict { pos: 9, open_pos: 2 })]
#[case::ring_dir_conflict_digit(b"C/1.CC\\1", ParseError::RingBondDirConflict { pos: 7, open_pos: 2 })]
#[case::ring_dir_conflict_percent(b"C/%12.CC\\%12", ParseError::RingBondDirConflict { pos: 9, open_pos: 2 })]
#[case::ring_dir_conflict_aromatic(b"c/1.cc\\1", ParseError::RingBondDirConflict { pos: 7, open_pos: 2 })]
#[case::group_dot_before_ring_digit(b"(.1)", ParseError::GroupLeadingConnector { pos: 1 })]
#[case::group_dot_before_ring_percent(b"(.%12)", ParseError::GroupLeadingConnector { pos: 1 })]
#[case::branch_dot_before_ring_digit(b"C(.1)", ParseError::LeadingRing { pos: 3 })]
#[case::branch_dot_before_ring_percent(b"C(.%12)", ParseError::LeadingRing { pos: 3 })]
#[case::group_dot_before_bond(b"(.-C)", ParseError::GroupLeadingConnector { pos: 1 })]
#[case::branch_dot_before_bond(b"C(.-C)", ParseError::LeadingBond { pos: 3 })]
#[case::leading_bond_after_dot_1(b"C.-C", ParseError::LeadingBond { pos: 2 })]
#[case::leading_bond_after_dot_2(b"C.=-C", ParseError::LeadingBond { pos: 2 })]
#[case::leading_stereobond_after_dot_up(b"C./C", ParseError::LeadingBond { pos: 2 })]
#[case::leading_stereobond_after_dot_down(b"C.\\C", ParseError::LeadingBond { pos: 2 })]
#[case::trailing_bond_dot_aromatic(b"C:.", ParseError::TrailingBond { pos: 1 })]
#[case::group_trailing_bond_dot(b"(C-.)", ParseError::TrailingBond { pos: 2 })]
#[case::branch_trailing_bond_dot(b"C(C-.)", ParseError::TrailingBond { pos: 3 })]
fn components_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let err = parse_smiles(input);
    assert!(err.is_err(), "{:?} should have failed", input);
    let err = err.unwrap_err();
    assert_eq!(err, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::atom_c(b"[C]", Some(Element::C), false, None, None, None, None, None)]
#[case::atom_h(b"[H]", Some(Element::H), false, None, None, None, None, None)]
#[case::atom_zn(b"[Zn]", Some(Element::Zn), false, None, None, None, None, None)]
#[case::atom_og(b"[Og]", Some(Element::Og), false, None, None, None, None, None)]
#[case::atom_aromatic_c(b"[c]", Some(Element::C), true, None, None, None, None, None)]
#[case::atom_aromatic_se(b"[se]", Some(Element::Se), true, None, None, None, None, None)]
#[case::wildcard(b"[*]", None, false, None, None, None, None, None)]
#[case::isotope_element(b"[13C]", Some(Element::C), false, Some(13), None, None, None, None)]
#[case::isotope_zero(b"[0C]", Some(Element::C), false, Some(0), None, None, None, None)]
#[case::isotope_wildcard(b"[13*]", None, false, Some(13), None, None, None, None)]
#[case::isotope_zero_prefix_1(b"[02H]", Some(Element::H), false, Some(2), None, None, None, None)]
#[case::isotope_zero_prefix_2(b"[002H]", Some(Element::H), false, Some(2), None, None, None, None)]
#[case::isotope_three_digits_1(b"[238U]", Some(Element::U), false, Some(238), None, None, None, None)]
#[case::isotope_three_digits_2(b"[208Pb]", Some(Element::Pb), false, Some(208), None, None, None, None)]
#[case::isotope_unstable(b"[36Cl]", Some(Element::Cl), false, Some(36), None, None, None, None)]
#[case::isotope_max_999(b"[999Og]", Some(Element::Og), false, Some(999), None, None, None, None)]
#[case::isotope_hcount(b"[13CH4]", Some(Element::C), false, Some(13), None, Some(4), None, None)]
#[case::isotope_charge(b"[2H+]", Some(Element::H), false, Some(2), None, None, Some(1), None)]
#[case::chirality_cw(b"[C@]", Some(Element::C), false, None, Some(Chirality::Clockwise), None, None, None)]
#[case::chirality_ccw(b"[C@@]", Some(Element::C), false, None, Some(Chirality::CounterClockwise), None, None, None)]
#[case::chirality_th2(b"[C@TH2]", Some(Element::C), false, None, Some(Chirality::Tetrahedral { arr: 2 }), None, None, None)]
#[case::chirality_al1(b"[C@AL1]", Some(Element::C), false, None, Some(Chirality::Allenal { arr: 1 }), None, None, None)]
#[case::chirality_sp3(b"[C@SP3]", Some(Element::C), false, None, Some(Chirality::SquarePlanar { arr: 3 }), None, None, None)]
#[case::chirality_tb5(b"[C@TB5]", Some(Element::C), false, None, Some(Chirality::TrigonalBipyramidal { arr: 5 }), None, None, None)]
#[case::chirality_oh7(b"[C@OH7]", Some(Element::C), false, None, Some(Chirality::Octahedral { arr: 7 }), None, None, None)]
#[case::hcount(b"[CH]", Some(Element::C), false, None, None, Some(1), None, None)]
#[case::hcount_1(b"[CH1]", Some(Element::C), false, None, None, Some(1), None, None)]
#[case::hcount_0(b"[CH0]", Some(Element::C), false, None,None, Some(0), None, None)]
#[case::hcount_3(b"[CH3]", Some(Element::C), false, None, None, Some(3), None, None)]
#[case::hcount_aromatic(b"[cH]", Some(Element::C), true, None, None, Some(1), None, None)]
#[case::hcount_two_characters_1(b"[ClH]", Some(Element::Cl), false, None, None, Some(1), None, None)]
#[case::hcount_two_character_2(b"[ClH1]", Some(Element::Cl), false, None, None, Some(1), None, None)]
#[case::wildcard_h1(b"[*H]", None, false, None, None, Some(1), None, None)]
#[case::wildcard_h2(b"[*H2]", None, false, None, None, Some(2), None, None)]
#[case::wildcard_h0(b"[*H0]", None, false, None, None, Some(0), None, None)]
#[case::chirality_cw_hydrogen(b"[C@H]", Some(Element::C), false, None, Some(Chirality::Clockwise), Some(1), None, None)]
#[case::chirality_ccw_hydrogen(b"[C@@H]", Some(Element::C), false, None, Some(Chirality::CounterClockwise), Some(1), None, None)]
#[case::charge_plus(b"[C+]", Some(Element::C), false, None, None, None, Some(1), None)]
#[case::charge_minus(b"[C-]", Some(Element::C), false, None, None, None, Some(-1), None)]
#[case::charge_pp(b"[C++]", Some(Element::C), false, None, None, None, Some(2), None)]
#[case::charge_mm(b"[C--]", Some(Element::C), false, None, None, None, Some(-2), None)]
#[case::zero_charge_pos(b"[C+0]", Some(Element::C), false, None, None, None, Some(0), None)]
#[case::zero_charge_neg(b"[C-0]", Some(Element::C), false, None, None, None, Some(0), None)]
#[case::charge_plus_15(b"[C+15]", Some(Element::C), false, None, None, None, Some(15), None)]
#[case::charge_minus_15(b"[C-15]", Some(Element::C), false, None, None, None, Some(-15), None)]
#[case::charge_two_characters_plus_1(b"[Na+]", Some(Element::Na), false, None, None, None, Some(1), None)]
#[case::charge_two_characters_plus_2(b"[Ca+2]", Some(Element::Ca), false, None, None, None, Some(2), None)]
#[case::charge_two_characters_pp(b"[Ca++]", Some(Element::Ca), false, None, None, None, Some(2), None)]
#[case::charge_two_characters_minus_1(b"[Cl-]", Some(Element::Cl), false, None, None, None, Some(-1), None)]
#[case::charge_two_characters_minus_2(b"[Se-2]", Some(Element::Se), false, None, None, None, Some(-2), None)]
#[case::charge_two_characters_mm(b"[Se--]", Some(Element::Se), false, None, None, None, Some(-2), None)]
#[case::charge_plus_hcount(b"[C+H]", Some(Element::C), false, None, None, Some(1), Some(1), None)]
#[case::charge_plus_1_hcount(b"[C+1H]", Some(Element::C), false, None, None, Some(1), Some(1), None)]
#[case::charge_minus_hcount(b"[C-H]", Some(Element::C), false, None, None, Some(1), Some(-1), None)]
#[case::charge_minus_1_hcount(b"[C-1H]", Some(Element::C), false, None, None, Some(1), Some(-1), None)]
#[case::hcount_charge_pos_1(b"[NH+]", Some(Element::N), false, None, None, Some(1), Some(1), None)]
#[case::hcount_charge_pos_2(b"[NH+1]", Some(Element::N), false, None, None, Some(1), Some(1), None)]
#[case::hcount_charge_pos_two_characters_1(b"[NaH+]", Some(Element::Na), false, None, None, Some(1), Some(1), None)]
#[case::hcount_charge_pos_two_characters_2(b"[AlH+2]", Some(Element::Al), false, None, None, Some(1), Some(2), None)]
#[case::hcount_charge_pos_two_characters_pp(b"[AlH++]", Some(Element::Al), false, None, None, Some(1), Some(2), None)]
#[case::hcount_charge_neg_1(b"[NH-]", Some(Element::N), false, None, None, Some(1), Some(-1), None)]
#[case::hcount_charge_neg_2(b"[NH-1]", Some(Element::N), false, None, None, Some(1), Some(-1), None)]
#[case::hcount_charge_neg_3(b"[N-H1]", Some(Element::N), false, None, None, Some(1), Some(-1), None)]
#[case::hcount_charge_neg_two_characters_1(b"[AsH-]", Some(Element::As), false, None, None, Some(1), Some(-1), None)]
#[case::hcount_charge_neg_two_characters_2(b"[AsH-2]", Some(Element::As), false, None, None, Some(1), Some(-2), None)]
#[case::hcount_charge_neg_two_characters_mm(b"[AsH--]", Some(Element::As), false, None, None, Some(1), Some(-2), None)]
#[case::class_elem(b"[C:12]", Some(Element::C), false, None, None, None, None, Some(12))]
#[case::class_wildcard(b"[*:5]", None, false, None, None, None, None, Some(5))]
#[case::class_zero(b"[C:0]", Some(Element::C), false, None, None, None, None, Some(0))]
#[case::class_zero_prefix_1(b"[C:03]", Some(Element::C), false, None, None, None, None, Some(3))]
#[case::class_zero_prefix_2(b"[C:003]", Some(Element::C), false, None, None, None, None, Some(3))]
#[case::class_max_9999(b"[C:9999]", Some(Element::C), false, None, None, None, None, Some(9999))]
#[case::ordering_1(b"[C@H+1:2]", Some(Element::C), false, None, Some(Chirality::Clockwise), Some(1), Some(1), Some(2))]
#[case::ordering_2(b"[CH@+1:2]", Some(Element::C), false, None, Some(Chirality::Clockwise), Some(1), Some(1), Some(2))]
#[case::ordering_3(b"[CH+1@:2]", Some(Element::C), false, None, Some(Chirality::Clockwise), Some(1), Some(1), Some(2))]
#[case::ordering_4(b"[CH+1:2@]", Some(Element::C), false, None, Some(Chirality::Clockwise), Some(1), Some(1), Some(2))]
#[case::ordering_5(b"[C+1@H:2]", Some(Element::C), false, None, Some(Chirality::Clockwise), Some(1), Some(1), Some(2))]
#[case::ordering_6(b"[C:2@H+1]", Some(Element::C), false, None, Some(Chirality::Clockwise), Some(1), Some(1), Some(2))]
fn bracket(
    #[case] input: &[u8],
    #[case] elem: Option<Element>,
    #[case] aromatic: bool,
    #[case] isotope: Option<u32>,
    #[case] chirality: Option<Chirality>,
    #[case] hcount: Option<u32>,
    #[case] charge: Option<i32>,
    #[case] class_: Option<u32>,
) {
    let res = parse_smiles(input);
    assert!(res.is_ok(), "{:?} should have succeeded", input);
    let mol = res.unwrap();
    assert_eq!(mol.atoms.len(), 1, "expected single atom");
    let a = &mol.atoms[0];
    match elem {
        Some(e) => match &a.symbol {
            AtomSymbol::Element(el) => assert_eq!(*el, e),
            other => panic!("expected element {:?}, got {:?}", e, other),
        },
        None => assert!(matches!(a.symbol, AtomSymbol::Unknown)),
    }
    assert_eq!(a.aromatic, Some(aromatic));
    assert_eq!(a.isotope, isotope);
    assert_eq!(a.chirality, chirality);
    assert_eq!(a.hydrogen_count, hcount);
    assert_eq!(a.charge, charge);
    assert_eq!(a.class, class_);
}

#[rstest]
#[case::aliphatic_before(b"C[C]", BondOrder::Single, None)]
#[case::aliphatic_before_single(b"C-[C]", BondOrder::Single, None)]
#[case::aliphatic_before_double(b"C=[C]", BondOrder::Double, None)]
#[case::aliphatic_before_triple(b"C#[C]", BondOrder::Triple, None)]
#[case::aliphatic_before_quadruple(b"C$[C]", BondOrder::Quadruple, None)]
#[case::aliphatic_before_aromatic(b"C:[C]", BondOrder::Aromatic, None)]
#[case::aliphatic_before_up(b"C/[C]", BondOrder::Single, Some(BondDir::Up))]
#[case::aliphatic_before_down(b"C\\[C]", BondOrder::Single, Some(BondDir::Down))]
#[case::aliphatic_after(b"[C]C", BondOrder::Single, None)]
#[case::aliphatic_after_single(b"[C]-C", BondOrder::Single, None)]
#[case::aliphatic_after_double(b"[C]=C", BondOrder::Double, None)]
#[case::aliphatic_after_triple(b"[C]#C", BondOrder::Triple, None)]
#[case::aliphatic_after_quadruple(b"[C]$C", BondOrder::Quadruple, None)]
#[case::aliphatic_after_aromatic(b"[C]:C", BondOrder::Aromatic, None)]
#[case::aliphatic_after_up(b"[C]/C", BondOrder::Single, Some(BondDir::Up))]
#[case::aliphatic_after_down(b"[C]\\C", BondOrder::Single, Some(BondDir::Down))]
#[case::aromatic_before(b"c[c]", BondOrder::Aromatic, None)]
#[case::aromatic_before_single(b"c-[c]", BondOrder::Single, None)]
#[case::aromatic_before_aromatic(b"c:[c]", BondOrder::Aromatic, None)]
#[case::aromatic_after(b"[c]c", BondOrder::Aromatic, None)]
#[case::aromatic_after_single(b"[c]-c", BondOrder::Single, None)]
#[case::aromatic_after_aromatic(b"[c]:c", BondOrder::Aromatic, None)]
#[case::aliphatic_before_aromatic(b"C[c]", BondOrder::Single, None)]
#[case::aliphatic_single_before_aromatic(b"C-[c]", BondOrder::Single, None)]
#[case::aliphatic_aromatic_before_aromatic(b"C:[c]", BondOrder::Aromatic, None)]
#[case::aliphatic_after_aromatic(b"[c]C", BondOrder::Single, None)]
#[case::aromatic_after_aliphatic(b"[C]c", BondOrder::Single, None)]
#[case::aromatic_after_aliphatic_single(b"[C]-c", BondOrder::Single, None)]
#[case::aromatic_after_aliphatic_aromatic(b"[c]:c", BondOrder::Aromatic, None)]
#[case::aromatic_after_aliphatic_up(b"[C]/c", BondOrder::Single, Some(BondDir::Up))]
#[case::aromatic_after_aliphatic_down(b"[C]\\c", BondOrder::Single, Some(BondDir::Down))]
#[case::bracket_branch_1(b"[C](C)", BondOrder::Single, None)]
#[case::bracket_branch_2(b"C([C])", BondOrder::Single, None)]
#[case::bracket_branch_single(b"C(-[C])", BondOrder::Single, None)]
#[case::bracket_branch_double(b"C(=[C])", BondOrder::Double, None)]
#[case::bracket_branch_triple(b"C(#[C])", BondOrder::Triple, None)]
#[case::bracket_branch_quadruple(b"C($[C])", BondOrder::Quadruple, None)]
#[case::bracket_branch_aromatic(b"C(:[C])", BondOrder::Aromatic, None)]
#[case::bracket_branch_up(b"C(/[C])", BondOrder::Single, Some(BondDir::Up))]
#[case::bracket_branch_down(b"C(\\[C])", BondOrder::Single, Some(BondDir::Down))]
#[case::bracket_branch_down(b"C(\\[C])", BondOrder::Single, Some(BondDir::Down))]
#[case::bracket_group_1(b"([C]C)", BondOrder::Single, None)]
#[case::bracket_group_1(b"(C[C])", BondOrder::Single, None)]
#[case::bracket_ring_1(b"[C]1CC1", BondOrder::Single, None)]
#[case::bracket_ring_2(b"[C]1cc1", BondOrder::Single, None)]
#[case::bracket_ring_double_1(b"[C]1=cc1", BondOrder::Double, None)]
#[case::bracket_ring_double_2(b"[C]=1cc1", BondOrder::Single, None)]
#[case::bracket_aromatic_ring(b"[c]1cc1", BondOrder::Aromatic, None)]
fn bracket_bonds(#[case] input: &[u8], #[case] expected: BondOrder, #[case] dir: Option<BondDir>) {
    let res = parse_smiles(input);
    assert!(res.is_ok(), "{:?} should have succeeded", input);
    let mol = res.unwrap();
    assert_eq!(mol.bonds[0].symbol, BondSymbol::Bond(expected));
    assert_eq!(mol.bonds[0].direction, dir);
}

#[rstest]
#[case::empty_bracket(b"[]", ParseError::InvalidBracket { pos: 0 })]
#[case::bracket_in_chain_empty(b"C[]", ParseError::InvalidBracket { pos: 1 })]
#[case::bracket_in_group_empty(b"(C[])", ParseError::InvalidBracket { pos: 2 })]
#[case::bracket_in_branch_empty(b"C([])C", ParseError::InvalidBracket { pos: 2 })]
#[case::bracket_in_component_empty(b"[].C", ParseError::InvalidBracket { pos: 0 })]
#[case::bracket_in_ring_empty(b"C1[]C1", ParseError::InvalidBracket { pos: 2 })]
#[case::double_bracket(b"[[C]]", ParseError::InvalidBracket { pos: 0 })]
#[case::invalid_element_1(b"[X]", ParseError::InvalidBracket { pos: 0 })]
#[case::invalid_element_2(b"[Z]", ParseError::InvalidBracket { pos: 0 })]
#[case::invalid_element_3(b"[Aq]", ParseError::InvalidBracket { pos: 0 })]
#[case::invalid_element_4(b"[Sh]", ParseError::InvalidBracket { pos: 0 })]
#[case::invalid_aromatic_element_1(b"[f]", ParseError::InvalidBracket { pos: 0 })]
#[case::invalid_aromatic_element_2(b"[ca]", ParseError::InvalidBracket { pos: 0 })]
#[case::two_elements_1(b"[CF]", ParseError::InvalidBracket { pos: 0 })]
#[case::two_elements_2(b"[AsF]", ParseError::InvalidBracket { pos: 0 })]
#[case::two_elements_3(b"[FAs]", ParseError::InvalidBracket { pos: 0 })]
#[case::two_elements_4(b"[AsBr]", ParseError::InvalidBracket { pos: 0 })]
#[case::two_elements_wildcard_1(b"[*C]", ParseError::InvalidBracket { pos: 0 })]
#[case::two_elements_wildcard_2(b"[C*]", ParseError::InvalidBracket { pos: 0 })]
#[case::wildcard_invalid_element_1(b"[*X]", ParseError::InvalidBracket { pos: 0 })]
#[case::wildcard_invalid_element_2(b"[X*]", ParseError::InvalidBracket { pos: 0 })]
#[case::double_wildcard(b"[**]", ParseError::InvalidBracket { pos: 0 })]
#[case::zero_charge_no_sign(b"[C0]", ParseError::InvalidBracket { pos: 0 })]
#[case::pos_charge_no_sign(b"[C1]", ParseError::InvalidBracket { pos: 0 })]
#[case::charge_no_element_1(b"[+]", ParseError::InvalidBracket { pos: 0 })]
#[case::charge_no_element_2(b"[-]", ParseError::InvalidBracket { pos: 0 })]
#[case::charge_no_element_3(b"[+0]", ParseError::InvalidBracket { pos: 0 })]
#[case::charge_no_element_4(b"[-0]", ParseError::InvalidBracket { pos: 0 })]
#[case::charge_no_element_5(b"[+1]", ParseError::InvalidBracket { pos: 0 })]
#[case::charge_no_element_6(b"[-1]", ParseError::InvalidBracket { pos: 0 })]
#[case::zero_isotope_no_element(b"[0]", ParseError::InvalidBracket { pos: 0 })]
#[case::isotope_no_element(b"[13]", ParseError::InvalidBracket { pos: 0 })]
#[case::chirality_no_element_1(b"[@]", ParseError::InvalidBracket { pos: 0 })]
#[case::chirality_no_element_2(b"[@@]", ParseError::InvalidBracket { pos: 0 })]
#[case::chirality_no_element_4(b"[@@TH1]", ParseError::InvalidBracket { pos: 0 })]
#[case::class_no_element(b"[:12]", ParseError::InvalidBracket { pos: 0 })]
#[case::hcount_two_digits_1(b"[CH10]", ParseError::InvalidBracket { pos: 0 })]
#[case::hcount_two_digits_2(b"[SeH10]", ParseError::InvalidBracket { pos: 0 })]
#[case::colon_no_class(b"[C:]", ParseError::BracketEmptyClass { pos: 0 })]
#[case::unbalanced_open_bracket_1(b"[", ParseError::UnbalancedOpenBracket { pos: 0 })]
#[case::unbalanced_open_bracket_2(b"C[", ParseError::UnbalancedOpenBracket { pos: 1 })]
#[case::unbalanced_open_bracket_3(b"[C", ParseError::UnbalancedOpenBracket { pos: 0 })]
#[case::unbalanced_open_bracket_4(b"[*", ParseError::UnbalancedOpenBracket { pos: 0 })]
#[case::unbalanced_open_bracket_5(b"[)", ParseError::UnbalancedOpenBracket { pos: 0 })]
#[case::unbalanced_open_bracket_6(b"[[", ParseError::UnbalancedOpenBracket { pos: 0 })]
#[case::unbalanced_open_bracket_7(b"[.", ParseError::UnbalancedOpenBracket { pos: 0 })]
#[case::unbalanced_open_bracket_8(b"C[", ParseError::UnbalancedOpenBracket { pos: 1 })]
#[case::unbalanced_open_bracket_9(b"[C)", ParseError::UnbalancedOpenBracket { pos: 0 })]
#[case::unbalanced_open_bracket_10(b"[.C", ParseError::UnbalancedOpenBracket { pos: 0 })]
#[case::unbalanced_open_bracket_11(b"C.[", ParseError::UnbalancedOpenBracket { pos: 2 })]
#[case::dot_in_bracket(b"[.]", ParseError::InvalidBracket { pos: 0 })]
#[case::branch_open_in_bracket(b"[(]", ParseError::InvalidBracket { pos: 0 })]
#[case::branch_close_in_bracket(b"[)]", ParseError::InvalidBracket { pos: 0 })]
#[case::bracket_in_bracket_1(b"[[]", ParseError::InvalidBracket { pos: 0 })]
#[case::bracket_in_bracket_2(b"[]]", ParseError::InvalidBracket { pos: 0 })]
#[case::open_bracket_in_branch(b"C([)", ParseError::UnbalancedOpenBracket { pos: 2 })]
#[case::close_bracket_in_branch(b"C(])", ParseError::UnbalancedCloseBracket { pos: 2 })]
#[case::unbalanced_close_bracket_1(b"]", ParseError::UnbalancedCloseBracket { pos: 0 })]
#[case::unbalanced_close_bracket_2(b"]C", ParseError::UnbalancedCloseBracket { pos: 0 })]
#[case::unbalanced_close_bracket_3(b"C]", ParseError::UnbalancedCloseBracket { pos: 1 })]
#[case::unbalanced_close_bracket_4(b"*]", ParseError::UnbalancedCloseBracket { pos: 1 })]
#[case::unbalanced_close_bracket_5(b"C.]", ParseError::UnbalancedCloseBracket { pos: 2 })]
#[case::unbalanced_close_bracket_6(b"].", ParseError::UnbalancedCloseBracket { pos: 0 })]
#[case::unbalanced_close_bracket_7(b"].C", ParseError::UnbalancedCloseBracket { pos: 0 })]
#[case::unbalanced_close_bracket_8(b"(]", ParseError::UnbalancedCloseBracket { pos: 1 })]
#[case::unbalanced_close_bracket_9(b"(C]", ParseError::UnbalancedCloseBracket { pos: 2 })]
// TODO: Allow this?
#[case::h_on_h_1(b"[HH]", ParseError::BracketHOnH { pos: 0 })]
#[case::h_on_h_2(b"[HH1]", ParseError::BracketHOnH { pos: 0 })]
#[case::h_on_h_3(b"[HH0]", ParseError::BracketHOnH { pos: 0 })]
#[case::duplicate_hcount_1(b"[CHH]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::duplicate_hcount_2(b"[CHH1]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::duplicate_hcount_3(b"[CH1H1]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::duplicate_hcount_4(b"[CH1H]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::duplicate_hcount_5(b"[CH+H]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::duplicate_charge_pos_1(b"[C++1]", ParseError::InvalidBracket { pos: 0 })]
#[case::duplicate_charge_pos_2(b"[C+1+1]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::duplicate_charge_pos_3(b"[C+1+]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::duplicate_charge_pos_4(b"[C+-]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::duplicate_charge_pos_5(b"[C+-1]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::duplicate_charge_pos_6(b"[C+1-1]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::duplicate_charge_pos_7(b"[C+1-]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::duplicate_charge_pos_8(b"[C+H+]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::duplicate_charge_neg_1(b"[C--1]", ParseError::InvalidBracket { pos: 0 })]
#[case::duplicate_charge_neg_2(b"[C-1-1]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::duplicate_charge_neg_3(b"[C-1-]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::duplicate_charge_neg_4(b"[C-+]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::duplicate_charge_neg_5(b"[C-+1]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::duplicate_charge_neg_6(b"[C-1+1]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::duplicate_charge_neg_7(b"[C-1+]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::duplicate_charge_neg_8(b"[C-H-]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::duplicate_chirality_1(b"[C@@@]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::duplicate_chirality_2(b"[C@TH1@@]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::duplicate_chirality_3(b"[C@H@]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::empty_class(b"[C:]", ParseError::BracketEmptyClass { pos: 0 })]
#[case::empty_class_two_characters(b"[Cl:]", ParseError::BracketEmptyClass { pos: 0 })]
#[case::empty_class_hcount(b"[Cl:H]", ParseError::BracketEmptyClass { pos: 0 })]
#[case::empty_class_charge_pos(b"[Na:+]", ParseError::BracketEmptyClass { pos: 0 })]
#[case::empty_class_charge_neg(b"[Cl:-]", ParseError::BracketEmptyClass { pos: 0 })]
#[case::empty_class_chirality_cw(b"[C:@]", ParseError::BracketEmptyClass { pos: 0 })]
#[case::empty_class_chirality_ccw(b"[C:@@]", ParseError::BracketEmptyClass { pos: 0 })]
#[case::empty_class_double_colon(b"[C::]", ParseError::BracketEmptyClass { pos: 0 })]
#[case::duplicate_class_1(b"[C:1:1]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::duplicate_class_2(b"[C:12:1]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::duplicate_class_3(b"[C:12:12]", ParseError::BracketDuplicateField { pos: 0 })]
#[case::duplicate_class_4(b"[C:1:12]", ParseError::BracketDuplicateField { pos: 0 })]
fn bracket_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let err = parse_smiles(input);
    assert!(err.is_err(), "{:?} should have failed", input);
    let err = err.unwrap_err();
    assert_eq!(err, expected);
}

#[rstest]
#[case::chirality_outside_bracket(b"C@C", ParseError::FieldOutsideBracket { pos: 1 })]
#[case::charge_outside_bracket(b"C+C", ParseError::FieldOutsideBracket { pos: 1 })]
// TODO: Should be InvalidElement
#[case::hcount_outside_bracket(b"CHC", ParseError::UnsupportedToken { pos: 1 })]
#[case::class_outside_bracket(b"C:1C", ParseError::RingUnclosed { open_pos: 2 })]
fn bracket_fields_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let err = parse_smiles(input);
    assert!(err.is_err(), "{:?} should have failed", input);
    let err = err.unwrap_err();
    assert_eq!(err, expected);
}

#[rstest]
#[case::wildcard(b"*", 0, false, false)]
#[case::two_wildcards(b"**", 0, true, false)]
#[case::wildcard_after_c(b"C*", 1, true, false)]
#[case::wildcard_before_c(b"*C", 0, true, false)]
#[case::wildcard_bond_single(b"C-*", 1, true, false)]
#[case::wildcard_bond_single_rev(b"*-C", 0, true, false)]
#[case::wildcard_branch_1(b"*(C)", 0, true, false)]
#[case::wildcard_branch_2(b"C(*)", 1, true, false)]
#[case::wildcard_branch_3(b"C(*C)", 1, true, false)]
#[case::wildcard_group_1(b"(*)", 0, false, false)]
#[case::wildcard_group_2(b"(*C)", 0, true, false)]
#[case::wildcard_group_3(b"(C*)", 1, true, false)]
#[case::wildcard_ring_1(b"*1CC1", 0, true, false)]
#[case::wildcard_ring_2(b"C1*C1", 1, true, false)]
#[case::wildcard_ring_3(b"C1C*1", 2, true, false)]
// TODO: Should be aromatic? Not important for the atom but for incident bonds
#[case::wildcard_ring_aromatic(b"c1*c1", 1, true, false)]
#[case::wildcard_component_1(b"*.C", 0, false, false)]
#[case::wildcard_component_2(b"C.*", 1, false, false)]
#[case::wildcard_dot_bond_1(b"*1.C1", 0, true, false)]
#[case::wildcard_dot_bond_2(b"C1.*1", 1, true, false)]
fn wildcard(
    #[case] input: &[u8],
    #[case] star_idx: usize,
    #[case] has_bonds: bool,
    #[case] aromatic: bool,
) {
    let res = parse_smiles(input);
    assert!(res.is_ok(), "{:?} should have succeeded", input);
    let mol = res.unwrap();
    assert!(star_idx < mol.atoms.len());
    let a = &mol.atoms[star_idx];
    assert!(matches!(a.symbol, AtomSymbol::Unknown));
    assert_eq!(a.isotope, Some(0));
    assert_eq!(a.charge, Some(0));
    assert_eq!(a.hydrogen_count, Some(0));
    assert_eq!(a.aromatic, Some(aromatic));
    assert_eq!(a.implicit_h, false);
    if has_bonds {
        assert!(mol.bonds.len() > 0);
    }
}

#[rstest]
#[case::wildcard_after_group(b"(C)*", ParseError::TopLevelGroupTrailing { pos: 2 })]
#[case::wildcard_unclosed_ring(b"*1", ParseError::RingUnclosed { open_pos: 1 })]
#[case::wildcard_unclosed_branch(b"C(*", ParseError::UnbalancedBranchOpen { pos: 1 })]
#[case::wildcard_unclosed_group(b"(C*", ParseError::UnbalancedBranchOpen { pos: 0 })]
#[case::wildcard_unclosed_bracket(b"[*", ParseError::UnbalancedOpenBracket { pos: 0 })]
#[case::wildcard_trailing_bond(b"*-", ParseError::TrailingBond { pos: 1 })]
#[case::wildcard_trailing_dot(b"*.", ParseError::TrailingDot { pos: 1 })]
fn wildcard_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let res = parse_smiles(input);
    assert!(res.is_err(), "{:?} should have failed", input);
    let err = res.unwrap_err();
    assert_eq!(err, expected);
}

#[rstest]
#[case::space(b" ", Molecule::default())]
#[case::tab(b"\t", Molecule::default())]
#[case::newline(b"\n", Molecule::default())]
#[case::cr(b"\r", Molecule::default())]
#[case::crlf(b"\r\n", Molecule::default())]
#[case::terminator_space(b"CC ", build_from_graph("C C | 0-1"))]
#[case::terminator_tab(b"CC\t", build_from_graph("C C | 0-1"))]
#[case::terminator_newline(b"CC\n", build_from_graph("C C | 0-1"))]
#[case::terminator_cr(b"CC\r", build_from_graph("C C | 0-1"))]
#[case::terminator_crlf(b"CC\r\n", build_from_graph("C C | 0-1"))]
fn whitespace_strict(#[case] input: &[u8], #[case] expected: Molecule) {
    let res = parse_smiles(input);
    assert!(res.is_ok(), "{:?} should have succeeded", input);
    let mol = res.unwrap();
    assert_eq!(mol, expected);
}

#[rstest]
#[case::leading_space(b" CC", ParseError::InvalidWhitespace { pos: 0 })]
#[case::leading_tab(b"\tCC", ParseError::InvalidWhitespace { pos: 0 })]
#[case::leading_newline(b"\nCC", ParseError::InvalidWhitespace { pos: 0 })]
#[case::leading_cr(b"\rCC", ParseError::InvalidWhitespace { pos: 0 })]
#[case::leading_crlf(b"\r\nCC", ParseError::InvalidWhitespace { pos: 0 })]
#[case::terminator_space_trailing_structure(b"CC CC", ParseError::InvalidWhitespace { pos: 2 })]
#[case::terminator_tab_trailing_structure(b"CC\tCC", ParseError::InvalidWhitespace { pos: 2 })]
#[case::terminator_cr_trailing_structure(b"CC\rCC", ParseError::InvalidWhitespace { pos: 2 })]
#[case::terminator_newline_trailing_structure(b"CC\nCC", ParseError::InvalidWhitespace { pos: 2 })]
#[case::terminator_crlf_trailing_structure(b"CC\r\nCC", ParseError::InvalidWhitespace { pos: 2 })]
fn whitespace_strict_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let res = parse_smiles(input);
    assert!(res.is_err(), "{:?} should have failed", input);
    let err = res.unwrap_err();
    assert_eq!(err, expected);
}

#[rstest]
#[case::ws_intertoken_spaces_flags(b"C C", build_from_graph("C C | 0-1"))]
#[case::ws_intertoken_tabs_flags(b"C\tC", build_from_graph("C C | 0-1"))]
#[case::ws_newlines_flags(b"C\nC", build_from_graph("C C | 0-1"))]
#[case::line_comment_flags(b"C// x\nC", build_from_graph("C C | 0-1"))]
#[case::block_comment_flags(b"C/* x */C", build_from_graph("C C | 0-1"))]
#[case::block_comment_multiline_flags(b"C/* x\n y */C", build_from_graph("C C | 0-1"))]
#[case::eoi_blank_line(b"C\n\nC", build_from_graph("C |"))]
#[case::eoi_blank_line_crlf(b"C\r\n\r\nC", build_from_graph("C |"))]
#[case::eoi_blank_line_with_comment(b"C\n/* comment */\n\nC", build_from_graph("C |"))]
fn whitespace_lenient(#[case] input: &[u8], #[case] expected: Molecule) {
    let flags = SmilesParseFlags::INTERTOKEN_WS
        | SmilesParseFlags::COMMENTS
        | SmilesParseFlags::EXPLICIT_EOI;
    let res = parse_smiles_inner(input, flags);
    assert!(res.is_ok(), "{:?} should have succeeded", input);
    let mol = res.unwrap();
    assert_eq!(mol, expected);
}

#[rstest]
#[case::split_halogen_ws(b"C l", ParseError::UnsupportedToken { pos: 2 })]
#[case::percent_ring_ws_split(b"C% 12", ParseError::RingIndexInvalid { pos: 1 })]
#[case::percent_ring_nl_split(b"C%\n12", ParseError::RingIndexInvalid { pos: 1 })]
#[case::unterminated_block_comment(b"C/* x", ParseError::UnterminatedBlockComment { pos: 1 })]
#[case::bracket_inner_ws(b"[ C ]", ParseError::InvalidBracket { pos: 0 })]
fn whitespace_lenient_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let flags = SmilesParseFlags::INTERTOKEN_WS
        | SmilesParseFlags::COMMENTS
        | SmilesParseFlags::EXPLICIT_EOI;
    let res = parse_smiles_inner(input, flags);
    assert!(res.is_err(), "{:?} should have failed", input);
    let err = res.unwrap_err();
    assert_eq!(err, expected);
}

mod utils;
