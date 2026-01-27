#![allow(clippy::too_many_arguments)]

use bstr::ByteSlice;
use pretty_assertions::assert_eq;
use rstest::*;
use umol_data::Element;

use super::super::*;
use super::utils::{
    build_extended_from_graph, find_extended_chiral_center, find_extended_stereo_bond,
};
use crate::table_ir::{AtomSymbol, BondOrder, BondWedge, Chirality, ExtendedMolecule};

#[rstest]
#[case::organic_c(b"C", build_extended_from_graph("C@0 |"))]
#[case::organic_b(b"B", build_extended_from_graph("B@0 |"))]
#[case::organic_n(b"N", build_extended_from_graph("N@0 |"))]
#[case::organic_o(b"O", build_extended_from_graph("O@0 |"))]
#[case::organic_s(b"S", build_extended_from_graph("S@0 |"))]
#[case::organic_p(b"P", build_extended_from_graph("P@0 |"))]
#[case::organic_f(b"F", build_extended_from_graph("F@0 |"))]
#[case::organic_cl(b"Cl", build_extended_from_graph("Cl@0 |"))]
#[case::organic_br(b"Br", build_extended_from_graph("Br@0 |"))]
#[case::organic_i(b"I", build_extended_from_graph("I@0 |"))]
#[case::organic_b_aromatic(b"b", build_extended_from_graph("B_@0 |"))]
#[case::organic_c_aromatic(b"c", build_extended_from_graph("C_@0 |"))]
#[case::organic_n_aromatic(b"n", build_extended_from_graph("N_@0 |"))]
#[case::organic_o_aromatic(b"o", build_extended_from_graph("O_@0 |"))]
#[case::organic_s_aromatic(b"s", build_extended_from_graph("S_@0 |"))]
#[case::organic_p_aromatic(b"p", build_extended_from_graph("P_@0 |"))]
fn element(#[case] input: &[u8], #[case] expected: ExtendedMolecule) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_ok(), "{:?} should have succeeded", input_str);
    let mol = res.unwrap();
    assert_eq!(mol, expected);
}

#[rstest]
#[case::element_h(b"H", ParseError::InvalidElement { pos: 0 })]
#[case::element_he(b"He", ParseError::InvalidElement { pos: 0 })]
#[case::element_al(b"Al", ParseError::InvalidElement { pos: 0 })]
#[case::element_q(b"Q", ParseError::InvalidElement { pos: 0 })]
#[case::element_f_aromatic(b"f", ParseError::InvalidElement { pos: 0 })]
fn element_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_err(), "{:?} should have failed", input_str);
    let err = res.unwrap_err();
    assert_eq!(err, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::empty(b"", ExtendedMolecule::empty())]
#[case::chain_c_1(b"C", build_extended_from_graph("C@0..1 |"))]
#[case::chain_c_5(b"CCCCC", build_extended_from_graph("C@0..1 C@1..2 C@2..3 C@3..4 C@4..5 | 0-1@1..2 1-2@2..3 2-3@3..4 3-4@4..5"))]
#[case::aromatic_c_6(b"cccccc", build_extended_from_graph("C_@0..1 C_@1..2 C_@2..3 C_@3..4 C_@4..5 C_@5..6 | 0-1:@1..2 1-2:@2..3 2-3:@3..4 3-4:@4..5 4-5:@5..6"))]
#[case::chain_mixed_5(b"CClOBrN", build_extended_from_graph("C@0..1 Cl@1..3 O@3..4 Br@4..6 N@6..7 | 0-1@1..3 1-2@3..4 2-3@4..6 3-4@6..7"))]
fn chain(#[case] input: &[u8], #[case] expected: ExtendedMolecule) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_ok(), "{:?} should have succeeded", input_str);
    let mol = res.unwrap();
    assert_eq!(mol, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::empty_group(b"()", ExtendedMolecule::empty())]
#[case::group_c_1(b"(C)", build_extended_from_graph("C@1 |"))]
#[case::group_c_1_aromatic(b"(c)", build_extended_from_graph("C_@1 |"))]
#[case::group_c_4(b"(CCCC)", build_extended_from_graph("C@1 C@2 C@3 C@4 | 0-1@2 1-2@3 2-3@4"))]
#[case::group_nested(b"((CC))", build_extended_from_graph("C@2 C@3 | 0-1@3"))]
#[case::branch_c_211(b"CC(C)C", build_extended_from_graph("C@0 C@1 C@3 C@5 | 0-1@1 1-2@3 1-3@5"))]
#[case::branch_c_222_aromatic(b"cc(cc)cc", build_extended_from_graph("C_@0 C_@1 C_@3 C_@4 C_@6 C_@7 | 0-1:@1 1-2:@3 2-3:@4 1-4:@6 4-5:@7"))]
#[case::branch_trailing(b"C(CC)", build_extended_from_graph("C@0 C@2 C@3 | 0-1@2 1-2@3"))]
#[case::branch_multiple(b"CC(C)(C)C", build_extended_from_graph("C@0 C@1 C@3 C@6 C@8 | 0-1@1 1-2@3 1-3@6 1-4@8"))]
#[case::branch_multiple_trailing(b"C(C)(C)", build_extended_from_graph("C@0 C@2 C@5 | 0-1@2 0-2@5"))]
#[case::branch_nested(b"C(C(C)C)C", build_extended_from_graph("C@0 C@2 C@4 C@6 C@8 | 0-1@2 1-2@4 1-3@6 0-4@8"))]
fn tree(#[case] input: &[u8], #[case] expected: ExtendedMolecule) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_ok(), "{:?} should have succeeded", input_str);
    let mol = res.unwrap();
    assert_eq!(mol, expected);
}

#[rstest]
#[case::unbalanced_closing_paren_1(b")C", ParseError::UnbalancedCloseParen { pos: 0 })]
#[case::unbalanced_closing_paren_2(b"C)C", ParseError::UnbalancedCloseParen { pos: 1 })]
#[case::unclosed_group(b"(C", ParseError::UnbalancedOpenParen { pos: 0 })]
#[case::unclosed_branch(b"C(C", ParseError::UnbalancedOpenParen { pos: 1 })]
#[case::empty_branch(b"C()", ParseError::EmptyBranch { pos: 2 })]
#[case::empty_group_before_atom(b"()C", ParseError::EmptyGroup { pos: 1 })]
#[case::two_top_level_groups(b"(C)(C)", ParseError::NonfinalGroup { pos: 2 })]
#[case::three_top_level_groups(b"(C)(C)(C)", ParseError::NonfinalGroup { pos: 2 })]
#[case::three_top_level_groups_aromatic(b"(c)(c)(c)", ParseError::NonfinalGroup { pos: 2 })]
#[case::two_top_level_groups_rings(b"(C1CC1)(C2CC2)", ParseError::NonfinalGroup { pos: 6 })]
#[case::group_before_atom(b"(C)C", ParseError::NonfinalGroup { pos: 2 })]
#[case::group_before_atom_aromatic(b"(c)c", ParseError::NonfinalGroup { pos: 2 })]
fn tree_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_err(), "{:?} should have failed", input_str);
    let err = res.unwrap_err();
    assert_eq!(err, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::ring_3(b"C1CC1", build_extended_from_graph("C@0 C@2 C@3 | 0-1@2 1-2@3 0-2@1 | 1@1-4:0-2"))]
#[case::ring_6(b"C1CCCCC1", build_extended_from_graph("C@0 C@2 C@3 C@4 C@5 C@6 | 0-1@2 1-2@3 2-3@4 3-4@5 4-5@6 0-5@1 | 1@1-7:0-5"))]
#[case::ring_10(b"C1CCCCCCCCC1", build_extended_from_graph("C@0 C@2 C@3 C@4 C@5 C@6 C@7 C@8 C@9 C@10 | 0-1@2 1-2@3 2-3@4 3-4@5 4-5@6 5-6@7 6-7@8 7-8@9 8-9@10 0-9@1 | 1@1-11:0-9"))]
#[case::ring_aromatic(b"c1ccccc1", build_extended_from_graph("C_@0 C_@2 C_@3 C_@4 C_@5 C_@6 | 0-1:@2 1-2:@3 2-3:@4 3-4:@5 4-5:@6 0-5:@1 | 1@1-7:0-5"))]
#[case::ring_aromatic_anti(b"c1ccc1", build_extended_from_graph("C_@0 C_@2 C_@3 C_@4 | 0-1:@2 1-2:@3 2-3:@4 0-3:@1 | 1@1-5:0-3"))]
#[case::ring_aromatic_heteroatom(b"c1occc1", build_extended_from_graph("C_@0 O_@2 C_@3 C_@4 C_@5 | 0-1:@2 1-2:@3 2-3:@4 3-4:@5 0-4:@1 | 1@1-6:0-4"))]
#[case::ring_index_0(b"C0CC0", build_extended_from_graph("C@0 C@2 C@3 | 0-1@2 1-2@3 0-2@1 | 0@1-4:0-2"))]
#[case::ring_index_percent(b"C%12CC%12", build_extended_from_graph("C@0 C@4 C@5 | 0-1@4 1-2@5 0-2@1 | 12@1-6:0-2"))]
#[case::ring_index_percent_zero(b"C%00CC%00", build_extended_from_graph("C@0 C@4 C@5 | 0-1@4 1-2@5 0-2@1 | 0@1-6:0-2"))]
#[case::ring_index_percent_zero_prefix(b"C%01CC%01", build_extended_from_graph("C@0 C@4 C@5 | 0-1@4 1-2@5 0-2@1 | 1@1-6:0-2"))]
#[case::ring_index_zero_prefix_1(b"C1CC%01", build_extended_from_graph("C@0 C@2 C@3 | 0-1@2 1-2@3 0-2@1 | 1@1-4:0-2"))]
#[case::ring_index_zero_prefix_2(b"C0CC%00", build_extended_from_graph("C@0 C@2 C@3 | 0-1@2 1-2@3 0-2@1 | 0@1-4:0-2"))]
#[case::ring_index_zero_prefix_3(b"C9CC%09", build_extended_from_graph("C@0 C@2 C@3 | 0-1@2 1-2@3 0-2@1 | 9@1-4:0-2"))]
#[case::ring_index_max_99(b"C%99CC%99", build_extended_from_graph("C@0 C@4 C@5 | 0-1@4 1-2@5 0-2@1 | 99@1-6:0-2"))]
#[case::ring_indices_single_percent_1(b"C%123CCC%12CC3", build_extended_from_graph("C@0 C@5 C@6 C@7 C@11 C@12 | 0-1@5 1-2@6 2-3@7 0-3@1 3-4@11 4-5@12 0-5@4 | 12@1-8:0-3 3@4-13:0-5"))]
#[case::ring_indices_single_percent_2(b"C3%12CCC%12CC3", build_extended_from_graph("C@0 C@5 C@6 C@7 C@11 C@12 | 0-1@5 1-2@6 2-3@7 0-3@2 3-4@11 4-5@12 0-5@1 | 3@1-13:0-5 12@2-8:0-3"))]
#[case::two_rings_bonded_0(b"C1CC1C2CC2", build_extended_from_graph("C@0 C@2 C@3 C@5 C@7 C@8 | 0-1@2 1-2@3 0-2@1 2-3@5 3-4@7 4-5@8 3-5@6 | 1@1-4:0-2 2@6-9:3-5"))]
#[case::two_rings_bonded_0_aromatic_1(b"c1cc1c2cc2", build_extended_from_graph("C_@0 C_@2 C_@3 C_@5 C_@7 C_@8 | 0-1:@2 1-2:@3 0-2:@1 2-3:@5 3-4:@7 4-5:@8 3-5:@6 | 1@1-4:0-2 2@6-9:3-5"))]
#[case::two_rings_bonded_0_aromatic_2(b"c1cc1C2CC2", build_extended_from_graph("C_@0 C_@2 C_@3 C@5 C@7 C@8 | 0-1:@2 1-2:@3 0-2:@1 2-3@5 3-4@7 4-5@8 3-5@6 | 1@1-4:0-2 2@6-9:3-5"))]
#[case::two_rings_index_reused(b"C1CC1C1CC1", build_extended_from_graph("C@0 C@2 C@3 C@5 C@7 C@8 | 0-1@2 1-2@3 0-2@1 2-3@5 3-4@7 4-5@8 3-5@6 | 1@1-4:0-2 1@6-9:3-5"))]
#[case::two_rings_bonded_2(b"C1CC1CCC2CC2", build_extended_from_graph("C@0 C@2 C@3 C@5 C@6 C@7 C@9 C@10 | 0-1@2 1-2@3 0-2@1 2-3@5 3-4@6 4-5@7 5-6@9 6-7@10 5-7@8 | 1@1-4:0-2 2@8-11:5-7"))]
#[case::two_rings_spiro(b"C1CC12CC2", build_extended_from_graph("C@0 C@2 C@3 C@6 C@7 | 0-1@2 1-2@3 0-2@1 2-3@6 3-4@7 2-4@5 | 1@1-4:0-2 2@5-8:2-4"))]
#[case::two_rings_spiro_branch(b"C12(CCCCC1)CCCCC2", build_extended_from_graph("C@0 C@4 C@5 C@6 C@7 C@8 C@11 C@12 C@13 C@14 C@15 | 0-1@4 1-2@5 2-3@6 3-4@7 4-5@8 0-5@1 0-6@11 6-7@12 7-8@13 8-9@14 9-10@15 0-10@2 | 1@1-9:0-5 2@2-16:0-10"))]
#[case::two_rings_fused(b"C12CC1C2", build_extended_from_graph("C@0 C@3 C@4 C@6 | 0-1@3 1-2@4 0-2@1 2-3@6 0-3@2 | 1@1-5:0-2 2@2-7:0-3"))]
#[case::two_rings_bridged(b"C12CC(C2)C1", build_extended_from_graph("C@0 C@3 C@4 C@6 C@9 | 0-1@3 1-2@4 2-3@6 0-3@2 2-4@9 0-4@1 | 1@1-10:0-4 2@2-7:0-3"))]
#[case::two_rings_fused_aromatic(b"c12ccccc1cccc2", build_extended_from_graph("C_@0 C_@3 C_@4 C_@5 C_@6 C_@7 C_@9 C_@10 C_@11 C_@12 | 0-1:@3 1-2:@4 2-3:@5 3-4:@6 4-5:@7 0-5:@1 5-6:@9 6-7:@10 7-8:@11 8-9:@12 0-9:@2 | 1@1-8:0-5 2@2-13:0-9"))]
#[case::two_rings_fused_aromatic_aliphatic(b"c1ccc2CCCc2c1", build_extended_from_graph("C_@0 C_@2 C_@3 C_@4 C@6 C@7 C@8 C_@9 C_@11 | 0-1:@2 1-2:@3 2-3:@4 3-4@6 4-5@7 5-6@8 6-7@9 3-7:@5 7-8:@11 0-8:@1 | 1@1-12:0-8 2@5-10:3-7"))]
#[case::two_rings_interleaved_indices(b"N1CC2CCCCC2CC1", build_extended_from_graph("N@0 C@2 C@3 C@5 C@6 C@7 C@8 C@9 C@11 C@12 | 0-1@2 1-2@3 2-3@5 3-4@6 4-5@7 5-6@8 6-7@9 2-7@4 7-8@11 8-9@12 0-9@1 | 1@1-13:0-9 2@4-10:2-7"))]
#[case::three_rings_fused(b"C12C3C1C32", build_extended_from_graph("C@0 C@3 C@5 C@7 | 0-1@3 1-2@5 0-2@1 2-3@7 1-3@4 0-3@2 | 1@1-6:0-2 2@2-9:0-3 3@4-8:1-3"))]
#[case::ring_group(b"(C1CC1)", build_extended_from_graph("C@1 C@3 C@4 | 0-1@3 1-2@4 0-2@2 | 1@2-5:0-2"))]
#[case::ring_branch_1(b"CC(C1)(C1)", build_extended_from_graph("C@0 C@1 C@3 C@7 | 0-1@1 1-2@3 1-3@7 2-3@4 | 1@4-8:2-3"))]
#[case::ring_branch_2(b"C(C1)CC1", build_extended_from_graph("C@0 C@2 C@5 C@6 | 0-1@2 0-2@5 2-3@6 1-3@3 | 1@3-7:1-3"))]
#[case::substituted_ring_1(b"CC1CC1", build_extended_from_graph("C@0 C@1 C@3 C@4 | 0-1@1 1-2@3 2-3@4 1-3@2 | 1@2-5:1-3"))]
#[case::substituted_ring_2(b"C1(C)CC1", build_extended_from_graph("C@0 C@3 C@5 C@6 | 0-1@3 0-2@5 2-3@6 0-3@1 | 1@1-7:0-3"))]
#[case::substituted_ring_3(b"C(C)1CC1", build_extended_from_graph("C@0 C@2 C@5 C@6 | 0-1@2 0-2@5 2-3@6 0-3@4 | 1@4-7:0-3"))]
#[case::substituted_ring_4(b"C1C(C)C1", build_extended_from_graph("C@0 C@2 C@4 C@6 | 0-1@2 1-2@4 1-3@6 0-3@1 | 1@1-7:0-3"))]
#[case::substituted_ring_5(b"C1CC(C)1", build_extended_from_graph("C@0 C@2 C@3 C@5 | 0-1@2 1-2@3 2-3@5 0-2@1 | 1@1-7:0-2"))]
#[case::substituted_ring_6(b"C1CC1C", build_extended_from_graph("C@0 C@2 C@3 C@5 | 0-1@2 1-2@3 0-2@1 2-3@5 | 1@1-4:0-2"))]
#[case::substituted_ring_7(b"C1CC1(C)", build_extended_from_graph("C@0 C@2 C@3 C@6 | 0-1@2 1-2@3 0-2@1 2-3@6 | 1@1-4:0-2"))]
#[case::substituted_ring_aromatic(b"c1c(c)c1", build_extended_from_graph("C_@0 C_@2 C_@4 C_@6 | 0-1:@2 1-2:@4 1-3:@6 0-3:@1 | 1@1-7:0-3"))]
#[case::substituted_ring_branch(b"C1C(C(C)C)C1", build_extended_from_graph("C@0 C@2 C@4 C@6 C@8 C@10 | 0-1@2 1-2@4 2-3@6 2-4@8 1-5@10 0-5@1 | 1@1-11:0-5"))]
fn ring(#[case] input: &[u8], #[case] expected: ExtendedMolecule) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_ok(), "{:?} should have succeeded", input_str);
    let mol = res.unwrap();
    assert_eq!(mol, expected);
}

#[rstest]
#[case::leading_ring_0(b"0C", ParseError::LeadingRing { pos: 0 })]
#[case::leading_ring_1(b"1C", ParseError::LeadingRing { pos: 0 })]
#[case::leading_ring_percent(b"%12C", ParseError::LeadingRing { pos: 0 })]
#[case::leading_ring_group(b"(1CCC)", ParseError::LeadingRing { pos: 1 })]
#[case::leading_ring_branch(b"C(1CCC)", ParseError::LeadingRing { pos: 0 })]
#[case::ring_unclosed_1(b"C1CC", ParseError::UnbalancedRingIndex { open_pos: 1 })]
#[case::ring_unclosed_2(b"C1CC1C1", ParseError::UnbalancedRingIndex { open_pos: 6 })]
#[case::ring_unclosed_3(b"C1CC2C", ParseError::UnbalancedRingIndex { open_pos: 4 })]
#[case::ring_unclosed_self_loop(b"C111", ParseError::UnbalancedRingIndex { open_pos: 3 })]
#[case::ring_unclosed_percent(b"C%12CC", ParseError::UnbalancedRingIndex { open_pos: 1 })]
#[case::bad_percent_no_index_1(b"C%", ParseError::InvalidRingIndex { pos: 1 })]
#[case::bad_percent_no_index_2(b"C%C", ParseError::InvalidRingIndex { pos: 1 })]
#[case::bad_percent_single_digit_0(b"C%0", ParseError::InvalidRingIndex { pos: 1 })]
#[case::bad_percent_single_digit_1(b"C%1", ParseError::InvalidRingIndex { pos: 1 })]
#[case::bad_percent_char(b"C%1a", ParseError::InvalidRingIndex { pos: 1 })]
fn ring_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_err(), "{:?} should have failed", input_str);
    let err = res.unwrap_err();
    assert_eq!(err, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::ring_self_loop(b"C11", build_extended_from_graph("C@0 | 0-0@1 | 1@1-2:0-0"))]
#[case::ring_self_loop_percent(b"C%11%11", build_extended_from_graph("C@0 | 0-0@1 | 11@1-4:0-0"))]
#[case::ring_two_member(b"C1C1", build_extended_from_graph("C@0 C@2 | 0-1@2 0-1@1 | 1@1-3:0-1"))]
#[case::ring_two_member_multiple(b"C12C12", build_extended_from_graph("C@0 C@3 | 0-1@3 0-1@1 0-1@2 | 1@1-4:0-1 2@2-5:0-1"))]
#[case::ring_two_member_percent(b"C%12C%12", build_extended_from_graph("C@0 C@4 | 0-1@4 0-1@1 | 12@1-5:0-1"))]
#[case::ring_two_member_single_percent(b"C%123CCC%123", build_extended_from_graph("C@0 C@5 C@6 C@7 | 0-1@5 1-2@6 2-3@7 0-3@1 0-3@4 | 12@1-8:0-3 3@4-11:0-3"))]
#[case::ring_multiple_rings(b"C12CCCCC12", build_extended_from_graph("C@0 C@3 C@4 C@5 C@6 C@7 | 0-1@3 1-2@4 2-3@5 3-4@6 4-5@7 0-5@1 0-5@2 | 1@1-8:0-5 2@2-9:0-5"))]
#[case::ring_multiple_rings_triple(b"C123CCCCC123", build_extended_from_graph("C@0 C@4 C@5 C@6 C@7 C@8 | 0-1@4 1-2@5 2-3@6 3-4@7 4-5@8 0-5@1 0-5@2 0-5@3 | 1@1-9:0-5 2@2-10:0-5 3@3-11:0-5"))]
#[case::ring_multiple_rings_percent(b"C%12%13CCCCC%12%13", build_extended_from_graph("C@0 C@7 C@8 C@9 C@10 C@11 | 0-1@7 1-2@8 2-3@9 3-4@10 4-5@11 0-5@1 0-5@4 | 12@1-12:0-5 13@4-15:0-5"))]
fn ring_invalid_topology(#[case] input: &[u8], #[case] expected: ExtendedMolecule) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_ok(), "{:?} should have succeeded", input_str);
    let mol = res.unwrap();
    assert_eq!(mol, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::single_bond(b"C-C", build_extended_from_graph("C@0 C@2 | 0-1:-@1"))]
#[case::double_bond(b"C=C", build_extended_from_graph("C@0 C@2 | 0-1:=@1"))]
#[case::triple_bond(b"C#C", build_extended_from_graph("C@0 C@2 | 0-1:#@1"))]
#[case::quadruple_bond(b"C$C", build_extended_from_graph("C@0 C@2 | 0-1:$@1"))]
#[case::aromatic_bond(b"C:C", build_extended_from_graph("C@0 C@2 | 0-1::@1"))]
#[case::up_bond(b"C/C", build_extended_from_graph("C@0 C@2 | 0-1:/@1"))]
#[case::down_bond(b"C\\C", build_extended_from_graph("C@0 C@2 | 0-1:\\@1"))]
#[case::single_bond_aromatic(b"c-c", build_extended_from_graph("C_@0 C_@2 | 0-1:-@1"))]
#[case::double_bond_aromatic(b"c=c", build_extended_from_graph("C_@0 C_@2 | 0-1:=@1"))]
#[case::triple_bond_aromatic(b"c#c", build_extended_from_graph("C_@0 C_@2 | 0-1:#@1"))]
#[case::quadruple_bond_aromatic(b"c$c", build_extended_from_graph("C_@0 C_@2 | 0-1:$@1"))]
#[case::aromatic_bond_aromatic(b"c:c", build_extended_from_graph("C_@0 C_@2 | 0-1::@1"))]
#[case::up_bond_aromatic(b"c/c", build_extended_from_graph("C_@0 C_@2 | 0-1:/@1"))]
#[case::down_bond_aromatic(b"c\\c", build_extended_from_graph("C_@0 C_@2 | 0-1:\\@1"))]
#[case::allene_bonds(b"C=C=C", build_extended_from_graph("C@0 C@2 C@4 | 0-1:=@1 1-2:=@3"))]
#[case::conjugated_bonds(b"C=CC=C", build_extended_from_graph("C@0 C@2 C@3 C@5 | 0-1:=@1 1-2:-@3 2-3:=@4"))]
#[case::cumulene_bonds(b"C=C=C=C", build_extended_from_graph("C@0 C@2 C@4 C@6 | 0-1:=@1 1-2:=@3 2-3:=@5"))]
#[case::allene_bonds_aromatic(b"c=c=c", build_extended_from_graph("C_@0 C_@2 C_@4 | 0-1:=@1 1-2:=@3"))]
#[case::trans_bonds_1(b"C/C=C/C", build_extended_from_graph("C@0 C@2 C@4 C@6 | 0-1:/@1 1-2:=@3 2-3:/@5"))]
#[case::trans_bonds_2(b"C\\C=C\\C", build_extended_from_graph("C@0 C@2 C@4 C@6 | 0-1:\\@1 1-2:=@3 2-3:\\@5"))]
#[case::cis_bonds_1(b"C\\C=C/C", build_extended_from_graph("C@0 C@2 C@4 C@6 | 0-1:\\@1 1-2:=@3 2-3:/@5"))]
#[case::cis_bonds_2(b"C/C=C\\C", build_extended_from_graph("C@0 C@2 C@4 C@6 | 0-1:/@1 1-2:=@3 2-3:\\@5"))]
#[case::trans_cumulene_bonds(b"F/C=C=C=C/F", build_extended_from_graph("F@0 C@2 C@4 C@6 C@8 F@10 | 0-1:/@1 1-2:=@3 2-3:=@5 3-4:=@7 4-5:/@9"))]
#[case::cis_cumulene_bonds(b"F/C=C=C=C\\F", build_extended_from_graph("F@0 C@2 C@4 C@6 C@8 F@10 | 0-1:/@1 1-2:=@3 2-3:=@5 3-4:=@7 4-5:\\@9"))]
#[case::conjugated_bonds_aromatic(b"c=c-c=c", build_extended_from_graph("C_@0 C_@2 C_@4 C_@6 | 0-1:=@1 1-2:-@3 2-3:=@5"))]
#[case::branch_leading_single_bond(b"CC(-C)C", build_extended_from_graph("C@0 C@1 C@4 C@6 | 0-1@1 1-2@3 1-3@6"))]
#[case::branch_leading_single_bond_multiple(b"CC(-C)(-C)C", build_extended_from_graph("C@0 C@1 C@4 C@8 C@10 | 0-1@1 1-2@3 1-3@7 1-4@10"))]
#[case::branch_leading_double_bond(b"CC(=C)C", build_extended_from_graph("C@0 C@1 C@4 C@6 | 0-1@1 1-2:=@3 1-3@6"))]
#[case::branch_leading_double_bond_multiple(b"OS(=O)(=O)O", build_extended_from_graph("O@0 S@1 O@4 O@8 O@10 | 0-1@1 1-2:=@3 1-3:=@7 1-4@10"))]
#[case::branch_internal_bond(b"CC(C-C)C", build_extended_from_graph("C@0 C@1 C@3 C@5 C@7 | 0-1@1 1-2@3 2-3@4 1-4@7"))]
#[case::branch_internal_double_bond(b"CC(C=C)C", build_extended_from_graph("C@0 C@1 C@3 C@5 C@7 | 0-1@1 1-2@3 2-3:=@4 1-4@7"))]
#[case::branch_followed_by_bond(b"CC(C)-C", build_extended_from_graph("C@0 C@1 C@3 C@6 | 0-1@1 1-2@3 1-3@5"))]
#[case::branch_followed_by_double_bond(b"CC(C)=C", build_extended_from_graph("C@0 C@1 C@3 C@6 | 0-1@1 1-2@3 1-3:=@5"))]
#[case::branch_leading_bond_aromatic(b"cc(:c)c", build_extended_from_graph("C_@0 C_@1 C_@4 C_@6 | 0-1:@1 1-2:@3 1-3:@6"))]
#[case::branch_internal_bond_aromatic(b"cc(c:c)c", build_extended_from_graph("C_@0 C_@1 C_@3 C_@5 C_@7 | 0-1:@1 1-2:@3 2-3:@4 1-4:@7"))]
#[case::branch_followed_by_bond_aromatic(b"cc(c):c", build_extended_from_graph("C_@0 C_@1 C_@3 C_@6 | 0-1:@1 1-2:@3 1-3:@5"))]
#[case::branch_trans_double_bond_1(b"C/C=C/C", build_extended_from_graph("C@0 C@2 C@4 C@6 | 0-1:/@1 1-2:=@3 2-3:/@5"))]
#[case::branch_trans_double_bond_2(b"C\\C=C\\C", build_extended_from_graph("C@0 C@2 C@4 C@6 | 0-1:\\@1 1-2:=@3 2-3:\\@5"))]
#[case::branch_cis_double_bond_1(b"C\\C=C/C", build_extended_from_graph("C@0 C@2 C@4 C@6 | 0-1:\\@1 1-2:=@3 2-3:/@5"))]
#[case::branch_cis_double_bond_2(b"C/C=C\\C", build_extended_from_graph("C@0 C@2 C@4 C@6 | 0-1:/@1 1-2:=@3 2-3:\\@5"))]
#[case::ring_single_bond(b"C-1-C-C-1", build_extended_from_graph("C@0 C@4 C@6 | 0-1@3 1-2@5 0-2@2 | 1@2-8:0-2"))]
#[case::ring_single_bond_percent(b"C-%12-C-C-%12", build_extended_from_graph("C@0 C@6 C@8 | 0-1@5 1-2@7 0-2@2 | 12@2-10:0-2"))]
#[case::ring_double_bond_1(b"C1-C=C1", build_extended_from_graph("C@0 C@3 C@5 | 0-1@2 1-2:=@4 0-2@1 | 1@1-6:0-2"))]
#[case::ring_double_bond_2(b"C1-CC=1", build_extended_from_graph("C@0 C@3 C@4 | 0-1@2 1-2@4 0-2:=@1 | 1@1-6:0-2"))]
#[case::ring_double_bond_3(b"C=1-CC1", build_extended_from_graph("C@0 C@4 C@5 | 0-1@3 1-2@5 0-2:=@2 | 1@2-6:0-2"))]
#[case::ring_double_bond_4(b"C=1-C-C=1", build_extended_from_graph("C@0 C@4 C@6 | 0-1@3 1-2@5 0-2:=@2 | 1@2-8:0-2"))]
#[case::ring_double_bond_5(b"C=1CCCCC=1", build_extended_from_graph("C@0 C@3 C@4 C@5 C@6 C@7 | 0-1@3 1-2@4 2-3@5 3-4@6 4-5@7 0-5:=@2 | 1@2-9:0-5"))]
#[case::ring_double_bond_unilateral_close_1(b"C1CC=1", build_extended_from_graph("C@0 C@2 C@3 | 0-1@2 1-2@3 0-2:=@1 | 1@1-5:0-2"))]
#[case::ring_double_bond_unilateral_close_2(b"C1CCCCC=1", build_extended_from_graph("C@0 C@2 C@3 C@4 C@5 C@6 | 0-1@2 1-2@3 2-3@4 3-4@5 4-5@6 0-5:=@1 | 1@1-8:0-5"))]
#[case::ring_double_bond_unilateral_open_1(b"C=1CC1", build_extended_from_graph("C@0 C@3 C@4 | 0-1@3 1-2@4 0-2:=@2 | 1@2-5:0-2"))]
#[case::ring_double_bond_unilateral_open_2(b"C=1CCCCC1", build_extended_from_graph("C@0 C@3 C@4 C@5 C@6 C@7 | 0-1@3 1-2@4 2-3@5 3-4@6 4-5@7 0-5:=@2 | 1@2-8:0-5"))]
#[case::ring_triple_bond(b"C1-C-C#1", build_extended_from_graph("C@0 C@3 C@5 | 0-1@2 1-2@4 0-2:#@1 | 1@1-7:0-2"))]
#[case::ring_quadruple_bond(b"C1-C-C$1", build_extended_from_graph("C@0 C@3 C@5 | 0-1@2 1-2@4 0-2:$@1 | 1@1-7:0-2"))]
#[case::ring_aromatic_bond(b"c1:c:c:1", build_extended_from_graph("C_@0 C_@3 C_@5 | 0-1:@2 1-2:@4 0-2:@1 | 1@1-7:0-2"))]
#[case::ring_aromatic_single_bond(b"c1ccccc1-c2ccccc2", build_extended_from_graph("C_@0 C_@2 C_@3 C_@4 C_@5 C_@6 C_@9 C_@11 C_@12 C_@13 C_@14 C_@15 | 0-1:@2 1-2:@3 2-3:@4 3-4:@5 4-5:@6 0-5:@1 5-6@8 6-7:@11 7-8:@12 8-9:@13 9-10:@14 10-11:@15 6-11:@10 | 1@1-7:0-5 2@10-16:6-11"))]
#[case::ring_up_bond_1(b"C1CC/1", build_extended_from_graph("C@0 C@2 C@3 | 0-1@2 1-2@3 0-2:/@1 | 1@1-5:0-2"))]
#[case::ring_up_bond_2(b"C/1CC1", build_extended_from_graph("C@0 C@3 C@4 | 0-1@3 1-2@4 0-2:/@2 | 1@2-5:0-2"))]
#[case::ring_up_bond_3(b"C/1CC/1", build_extended_from_graph("C@0 C@3 C@4 | 0-1@3 1-2@4 0-2:/@2 | 1@2-6:0-2"))]
#[case::ring_down_bond(b"C1CC\\1", build_extended_from_graph("C@0 C@2 C@3 | 0-1@2 1-2@3 0-2:\\@1 | 1@1-5:0-2"))]
#[case::ring_down_bond_both(b"C\\1CC\\1", build_extended_from_graph("C@0 C@3 C@4 | 0-1@3 1-2@4 0-2:\\@2 | 1@2-6:0-2"))]
#[case::ring_up_bond_percent_open(b"C/%12CC%12", build_extended_from_graph("C@0 C@5 C@6 | 0-1@5 1-2@6 0-2:/@2 | 12@2-7:0-2"))]
#[case::ring_up_bond_percent_close(b"C%12CC/%12", build_extended_from_graph("C@0 C@4 C@5 | 0-1@4 1-2@5 0-2:/@1 | 12@1-7:0-2"))]
#[case::ring_down_bond_percent_both(b"C\\%12CC\\%12", build_extended_from_graph("C@0 C@5 C@6 | 0-1@5 1-2@6 0-2:\\@2 | 12@2-8:0-2"))]
#[case::ring_between_bonds(b"C1CC-1-C", build_extended_from_graph("C@0 C@2 C@3 C@7 | 0-1@2 1-2@3 0-2@1 2-3@6 | 1@1-5:0-2"))]
fn bonds(#[case] input: &[u8], #[case] expected: ExtendedMolecule) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_ok(), "{:?} should have succeeded", input_str);
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
#[case::trailing_stereo_bond_before_dot_down(b"C.\\C", ParseError::LeadingBond { pos: 2 })]
#[case::bond_after_group_1(b"(C)-", ParseError::NonfinalGroup { pos: 2 })]
#[case::bond_after_group_2(b"(C)=", ParseError::NonfinalGroup { pos: 2 })]
#[case::group_after_group_1(b"(C)(C)", ParseError::NonfinalGroup { pos: 2 })]
#[case::group_after_group_2(b"(c)(c)", ParseError::NonfinalGroup { pos: 2 })]
#[case::ring_after_group(b"(C1CCC)1", ParseError::NonfinalGroup { pos : 6})]
#[case::consecutive_bonds_1(b"C--C", ParseError::ConsecutiveBonds { pos: 2 })]
#[case::consecutive_bonds_2(b"C-=C", ParseError::ConsecutiveBonds { pos: 2 })]
#[case::consecutive_bonds_3(b"C-#C", ParseError::ConsecutiveBonds { pos: 2 })]
#[case::consecutive_bonds_4(b"C-$C", ParseError::ConsecutiveBonds { pos: 2 })]
#[case::consecutive_bonds_5(b"C-:C", ParseError::ConsecutiveBonds { pos: 2 })]
#[case::consecutive_stereo_bonds_1(b"C//C", ParseError::ConsecutiveBonds { pos: 2 })]
#[case::consecutive_stereo_bonds_2(b"C\\\\C", ParseError::ConsecutiveBonds { pos: 2 })]
#[case::consecutive_bond_and_stereo_bond_1(b"C-/C", ParseError::ConsecutiveBonds { pos: 2 })]
#[case::consecutive_bond_and_stereo_bond_2(b"C=\\C", ParseError::ConsecutiveBonds { pos: 2 })]
#[case::leading_bond_1(b"-C", ParseError::LeadingBond { pos: 0 })]
#[case::leading_bond_2(b"=C", ParseError::LeadingBond { pos: 0 })]
#[case::leading_bond_3(b"#C", ParseError::LeadingBond { pos: 0 })]
#[case::leading_bond_4(b"$C", ParseError::LeadingBond { pos: 0 })]
#[case::leading_aromatic_bond(b":C", ParseError::LeadingBond { pos: 0 })]
#[case::leading_stereobond_1(b"/C", ParseError::LeadingBond { pos: 0 })]
#[case::leading_stereobond_2(b"\\C", ParseError::LeadingBond { pos: 0 })]
#[case::group_leading_bond_1(b"(-C)C", ParseError::LeadingBond { pos: 1 })]
#[case::group_leading_bond_2(b"(=C)C", ParseError::LeadingBond { pos: 1 })]
#[case::group_leading_bond_3(b"(#C)C", ParseError::LeadingBond { pos: 1 })]
#[case::group_leading_bond_4(b"($C)C", ParseError::LeadingBond { pos: 1 })]
#[case::group_leading_stereobond_1(b"(/C)C", ParseError::LeadingBond { pos: 1 })]
#[case::group_leading_stereobond_2(b"(\\C)C", ParseError::LeadingBond { pos: 1 })]
#[case::group_leading_aromatic_bond(b"(:C)C", ParseError::LeadingBond { pos: 1 })]
#[case::ring_bond_order_conflict_1(b"C-1CCCCC=1", ParseError::MismatchedRingBondOrders { pos: 9, open_pos: 2 })]
#[case::ring_bond_order_conflict_2(b"C=1CCCCC-1", ParseError::MismatchedRingBondOrders { pos: 9, open_pos: 2 })]
#[case::ring_bond_order_conflict_3(b"C=1CC#1", ParseError::MismatchedRingBondOrders { pos: 6, open_pos: 2 })]
#[case::ring_bond_order_conflict_4(b"C/1CC=1", ParseError::MismatchedRingBondOrders { pos: 6, open_pos: 2 })]
#[case::ring_bond_order_conflict_5(b"C\\1CC=1", ParseError::MismatchedRingBondOrders { pos: 6, open_pos: 2 })]
#[case::ring_bond_order_conflict_6(b"C=1CC/1", ParseError::MismatchedRingBondOrders { pos: 6, open_pos: 2 })]
#[case::ring_bond_order_conflict_7(b"C=1CC\\1", ParseError::MismatchedRingBondOrders { pos: 6, open_pos: 2 })]
#[case::ring_bond_order_conflict_8(b"C=%10CC#%10", ParseError::MismatchedRingBondOrders { pos: 8, open_pos: 2 })]
#[case::ring_bond_dir_conflict_1(b"C/1CC\\1", ParseError::MismatchedRingBondDirs { pos: 6, open_pos: 2 })]
#[case::ring_bond_dir_conflict_2(b"C\\1CC/1", ParseError::MismatchedRingBondDirs { pos: 6, open_pos: 2 })]
#[case::ring_bond_dir_conflict_3(b"C/%12CC\\%12", ParseError::MismatchedRingBondDirs { pos: 8, open_pos: 2 })]
#[case::ring_bond_dir_conflict_4(b"C\\%12CC/%12", ParseError::MismatchedRingBondDirs { pos: 8, open_pos: 2 })]
#[case::extended_bonds_1(b"C~C", ParseError::InvalidToken { pos: 1 })]
#[case::extended_bonds_2(b"C->N", ParseError::InvalidToken { pos: 2 })]
#[case::extended_bonds_3(b"C<-N", ParseError::InvalidToken { pos: 1 })]
#[case::extended_bonds_consecutive(b"C~~C", ParseError::InvalidToken { pos: 1 })]
fn bonds_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_err(), "{:?} should have failed", input_str);
    let err = res.unwrap_err();
    assert_eq!(err, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::any_bond(b"C~C", build_extended_from_graph("C@0 C@2 | 0-1~@1"))]
#[case::any_bond_in_ring(b"C1~CC~C1", build_extended_from_graph("C@0 C@3 C@4 C@6 | 0-1~@2 1-2@4 2-3~@5 0-3@1 | 1@1-7:0-3"))]
#[case::any_ring_bond_1(b"C1CC~1", build_extended_from_graph("C@0 C@2 C@3 | 0-1@2 1-2@3 0-2~@1 | 1@1-5:0-2"))]
#[case::any_ring_bond_2(b"C~1CC1", build_extended_from_graph("C@0 C@3 C@4 | 0-1@3 1-2@4 0-2~@2 | 1@2-5:0-2"))]
#[case::any_ring_bond_3(b"C~1CC~1", build_extended_from_graph("C@0 C@3 C@4 | 0-1@3 1-2@4 0-2~@2 | 1@2-6:0-2"))]
#[case::dative_accepting_1(b"C<-N", build_extended_from_graph("C@0 N@3 | 0-1<-@1"))]
#[case::dative_accepting_2(b"N<-C", build_extended_from_graph("N@0 C@3 | 0-1<-@1"))]
#[case::dative_accepting_multiple(b"C<-N<-O", build_extended_from_graph("C@0 N@3 O@6 | 0-1<-@1 1-2<-@4"))]
#[case::dative_ring_bond_1(b"C<-1CC1", build_extended_from_graph("C@0 C@4 C@5 | 0-1@4 1-2@5 0-2<-@3 | 1@3-6:0-2"))]
#[case::dative_ring_bond_2(b"C1CC->1", build_extended_from_graph("C@0 C@2 C@3 | 0-1@2 1-2@3 0-2<-@1 | 1@1-6:0-2"))]
#[case::dative_ring_bond_3(b"C<-1CC->1", build_extended_from_graph("C@0 C@4 C@5 | 0-1@4 1-2@5 0-2<-@3 | 1@3-8:0-2"))]
#[case::dative_donating_1(b"C->N", build_extended_from_graph("C@0 N@3 | 0-1->@1"))]
#[case::dative_donating_2(b"N->C", build_extended_from_graph("N@0 C@3 | 0-1->@1"))]
#[case::dative_donating_multiple(b"C->N->O", build_extended_from_graph("C@0 N@3 O@6 | 0-1->@1 1-2->@4"))]
fn bonds_lenient(#[case] input: &[u8], #[case] expected: ExtendedMolecule) {
    let res = parse_extended_smiles_bytes_with(input, &SmilesIoConfig::lenient());
    assert!(res.is_ok(), "{:?} should have succeeded", input);
    let mol = res.unwrap();
    assert_eq!(mol, expected);
}

#[rstest]
#[case::any_bond_consecutive(b"C~~C", ParseError::ConsecutiveBonds { pos: 2 })]
#[case::dative_bond_consecutive_1(b"C->->C", ParseError::ConsecutiveBonds { pos: 3 })]
#[case::dative_bond_consecutive_2(b"C-><-C", ParseError::ConsecutiveBonds { pos: 3 })]
#[case::dative_bond_consecutive_3(b"C<-<-C", ParseError::ConsecutiveBonds { pos: 3 })]
#[case::dative_bond_consecutive_4(b"C<-->C", ParseError::ConsecutiveBonds { pos: 3 })]
#[case::any_dative_bond_consecutive_1(b"C~->C", ParseError::ConsecutiveBonds { pos: 2 })]
#[case::any_dative_bond_consecutive_2(b"C<-~C", ParseError::ConsecutiveBonds { pos: 3 })]
#[case::any_ring_bond_order_conflict_1(b"C~1CC-1", ParseError::MismatchedRingBondOrders { pos: 6, open_pos: 2 })]
#[case::any_ring_bond_order_conflict_2(b"C-1CC~1", ParseError::MismatchedRingBondOrders { pos: 6, open_pos: 2 })]
#[case::dative_ring_bond_donation_conflict_1(b"C->1CC->1", ParseError::MismatchedRingBondDonations { pos: 8, open_pos: 3 })]
#[case::dative_ring_bond_donation_conflict_2(b"C<-1CC<-1", ParseError::MismatchedRingBondDonations { pos: 8, open_pos: 3 })]
fn bonds_lenient_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let res = parse_extended_smiles_bytes_with(input, &SmilesIoConfig::lenient());
    let input_str = input.to_str_lossy();
    assert!(res.is_err(), "{:?} should have failed", input_str);
    let err = res.unwrap_err();
    assert_eq!(err, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::components_2(b"CC.CC", build_extended_from_graph("C@0 C@1 C@3 C@4 | 0-1@1 2-3@4"))]
#[case::components_5(b"C.C.C.C.C", build_extended_from_graph("C@0 C@2 C@4 C@6 C@8 | "))]
#[case::ring_components_1(b"C1.CC1", build_extended_from_graph("C@0 C@3 C@4 | 1-2@4 0-2@1 | 1@1-5:0-2"))]
#[case::ring_components_2(b"C%12.CC%12", build_extended_from_graph("C@0 C@5 C@6 | 1-2@6 0-2@1 | 12@1-7:0-2"))]
#[case::ring_components_3(b"C1.C12.C2", build_extended_from_graph("C@0 C@3 C@7 | 0-1@1 1-2@5 | 1@1-4:0-1 2@5-8:1-2"))]
#[case::ring_components_aromatic_1(b"c1.ccccc1", build_extended_from_graph("C_@0 C_@3 C_@4 C_@5 C_@6 C_@7 | 1-2:@4 2-3:@5 3-4:@6 4-5:@7 0-5:@1 | 1@1-8:0-5"))]
#[case::ring_components_aromatic_2(b"c1c2c3c4cc1.Br2.Cl3.Cl4", build_extended_from_graph("C_@0 C_@2 C_@4 C_@6 C_@8 C_@9 Br@12 Cl@16 Cl@20 | 0-1:@2 1-2:@4 2-3:@6 3-4:@8 4-5:@9 0-5:@1 1-6@3 2-7@5 3-8@7 | 1@1-10:0-5 2@3-14:1-6 3@5-18:2-7 4@7-22:3-8"))]
#[case::branch_components(b"C(C.C)", build_extended_from_graph("C@0 C@2 C@4 | 0-1@2"))]
#[case::branch_leading_dot_1(b"C(.C)", build_extended_from_graph("C@0 C@3 | "))]
#[case::branch_leading_dot_2(b"C(.C)(C)", build_extended_from_graph("C@0 C@3 C@6 | 0-2@6"))]
#[case::branch_leading_dot_3(b"C(.C.C)", build_extended_from_graph("C@0 C@3 C@5 |"))]
#[case::branch_leading_dot_4(b"C(C)(.C)", build_extended_from_graph("C@0 C@2 C@6 | 0-1@2"))]
#[case::branch_trailing_dot_1(b"C(C.)", build_extended_from_graph("C@0 C@2 | 0-1@2"))]
#[case::branch_trailing_dot_2(b"C(C.)C", build_extended_from_graph("C@0 C@2 C@5 | 0-1@2 0-2@5"))]
#[case::branch_trailing_dot_3(b"C(C.)(C)", build_extended_from_graph("C@0 C@2 C@6 | 0-1@2 0-2@6"))]
#[case::branch_inner_dot(b"C(C.C)C", build_extended_from_graph("C@0 C@2 C@4 C@6 | 0-1@2 0-3@6"))]
#[case::group_components_1(b"(C.CC.C)", build_extended_from_graph("C@1 C@3 C@4 C@6 | 1-2@4"))]
#[case::group_components_2(b"(CC).(CC)", build_extended_from_graph("C@1 C@2 C@6 C@7 | 0-1@2 2-3@7"))]
#[case::group_components_3(b"(C.C).C", build_extended_from_graph("C@1 C@3 C@6 |"))]
#[case::group_components_4(b"C.(C).C", build_extended_from_graph("C@0 C@3 C@6 |"))]
#[case::group_components_5(b"C.C.(C)", build_extended_from_graph("C@0 C@2 C@5 |"))]
#[case::group_trailing_dot_1(b"(CC.)", build_extended_from_graph("C@1 C@2 | 0-1@2"))]
#[case::group_trailing_dot_2(b"(CC.).CC", build_extended_from_graph("C@1 C@2 C@6 C@7 | 0-1@2 2-3@7"))]
#[case::group_trailing_dot_3(b"(CC).(CC.)", build_extended_from_graph("C@1 C@2 C@6 C@7 | 0-1@2 2-3@7"))]
#[case::branch_ring_components_1(b"C1(C.C)CC1", build_extended_from_graph("C@0 C@3 C@5 C@7 C@8 | 0-1@3 0-3@7 3-4@8 0-4@1 | 1@1-9:0-4"))]
#[case::branch_ring_components_2(b"C1(C.C1)CC", build_extended_from_graph("C@0 C@3 C@5 C@8 C@9 | 0-1@3 0-2@1 0-3@8 3-4@9 | 1@1-6:0-2"))]
#[case::branch_ring_components_3(b"C(C1.C)CC1", build_extended_from_graph("C@0 C@2 C@5 C@7 C@8 | 0-1@2 0-3@7 3-4@8 1-4@3 | 1@3-9:1-4"))]
#[case::group_ring_components_1(b"(CC1.C1)", build_extended_from_graph("C@1 C@2 C@5 | 0-1@2 1-2@3 | 1@3-6:1-2"))]
#[case::group_ring_components_2(b"C1.(C).CC1", build_extended_from_graph("C@0 C@4 C@7 C@8 | 2-3@8 0-3@1 | 1@1-9:0-3"))]
#[case::group_ring_components_3(b"C%12.(C).CC%12", build_extended_from_graph("C@0 C@6 C@9 C@10 | 2-3@10 0-3@1 | 12@1-11:0-3"))]
#[case::rings_across_multiple_dots_digit(b"C1.C.CC1", build_extended_from_graph("C@0 C@3 C@5 C@6 | 2-3@6 0-3@1 | 1@1-7:0-3"))]
#[case::rings_across_multiple_dots_percent(b"C%12.C.CC%12", build_extended_from_graph("C@0 C@5 C@7 C@8 | 2-3@8 0-3@1 | 12@1-9:0-3"))]
#[case::ring_double_unilateral_open(b"C=1.CC1", build_extended_from_graph("C@0 C@4 C@5 | 1-2@5 0-2:=@2 | 1@2-6:0-2"))]
#[case::ring_double_unilateral_close(b"C1.CC=1", build_extended_from_graph("C@0 C@3 C@4 | 1-2@4 0-2:=@1 | 1@1-6:0-2"))]
#[case::ring_dir_up_both(b"C/1.CC/1", build_extended_from_graph("C@0 C@4 C@5 | 1-2@5 0-2:/@2 | 1@2-7:0-2"))]
#[case::ring_dir_down_both(b"C\\1.CC\\1", build_extended_from_graph("C@0 C@4 C@5 | 1-2@5 0-2:\\@2 | 1@2-7:0-2"))]
#[case::ring_dir_up_both_percent(b"C/%12.CC/%12", build_extended_from_graph("C@0 C@6 C@7 | 1-2@7 0-2:/@2 | 12@2-9:0-2"))]
#[case::ring_dir_down_both_percent(b"C\\%12CC\\%12", build_extended_from_graph("C@0 C@5 C@6 | 0-1@5 1-2@6 0-2:\\@2 | 12@2-8:0-2"))]
#[case::branch_multiple_components(b"C(.C.C)", build_extended_from_graph("C@0 C@3 C@5 |"))]
fn components(#[case] input: &[u8], #[case] expected: ExtendedMolecule) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_ok(), "{:?} should have succeeded", input_str);
    let mol = res.unwrap();
    assert_eq!(mol, expected);
}

#[rstest]
#[case::leading_dot_1(b".", ParseError::LeadingDot { pos: 0 })]
#[case::leading_dot_2(b".C", ParseError::LeadingDot { pos: 0 })]
#[case::leading_dot_3(b"..C", ParseError::LeadingDot { pos: 0 })]
#[case::leading_dot_4(b".C.", ParseError::LeadingDot { pos: 0 })]
#[case::trailing_dot_1(b"C.", ParseError::TrailingDot { pos: 1 })]
#[case::trailing_dot_2(b"C..", ParseError::ConsecutiveDots { pos: 1 })]
#[case::double_dot(b"C..C", ParseError::ConsecutiveDots { pos: 1 })]
#[case::dot_before_ring_digit(b"C.1", ParseError::DotBeforeRing { pos: 1 })]
#[case::dot_before_ring_digit_close(b"C1CCCCC.1", ParseError::DotBeforeRing { pos: 7 })]
#[case::dot_before_ring_digit_both(b"C.1CCCCC.1", ParseError::DotBeforeRing { pos: 1 })]
#[case::dot_before_ring_percent(b"C.%12", ParseError::DotBeforeRing { pos: 1 })]
#[case::dot_in_group_1(b"(.)", ParseError::LeadingDot { pos: 1 })]
#[case::dot_in_group_2(b"(.)C", ParseError::LeadingDot { pos: 1 })]
#[case::dot_in_group_3(b"(.).C", ParseError::LeadingDot { pos: 1 })]
#[case::leading_dot_in_group_1(b"(.CC)", ParseError::LeadingDot { pos: 1 })]
#[case::leading_dot_in_group_2(b"(.CC).(CC)", ParseError::LeadingDot { pos: 1 })]
#[case::leading_dot_in_group_3(b"(CC).(.CC)", ParseError::LeadingDot { pos: 6 })]
#[case::leading_dot_in_group_4(b"(.C).(.C)", ParseError::LeadingDot { pos: 1 })]
#[case::leading_dot_in_group_5(b"C.(.C).C", ParseError::LeadingDot { pos: 3 })]
#[case::dot_before_group(b"C.(C)C", ParseError::NonfinalGroup { pos: 4 })]
#[case::dot_in_branch_1(b"C(.)", ParseError::EmptyBranch { pos: 3 })]
#[case::dot_in_branch_2(b"C(.)C", ParseError::EmptyBranch { pos: 3 })]
#[case::dot_in_branch_3(b"C(.)(C)", ParseError::EmptyBranch { pos: 3 })]
#[case::dot_in_component_1(b"().C", ParseError::EmptyGroup { pos: 1 })]
#[case::dot_in_component_2(b"(.).C", ParseError::LeadingDot { pos: 1 })]
#[case::dot_in_component_3(b"(.).(C)", ParseError::LeadingDot { pos: 1})]
#[case::dot_in_component_4(b"C.()", ParseError::EmptyGroup { pos: 3 })]
#[case::dot_in_component_5(b"C.(.)", ParseError::LeadingDot { pos: 3 })]
#[case::dot_in_component_6(b"(C).(.)", ParseError::LeadingDot { pos: 5 })]
#[case::dot_unclosed_ring_1(b"C1.C", ParseError::UnbalancedRingIndex { open_pos: 1 })]
#[case::dot_unclosed_ring_2(b"C.C1", ParseError::UnbalancedRingIndex { open_pos: 3 })]
#[case::dot_unclosed_ring_before_group(b"C1.(C)(C)C1", ParseError::NonfinalGroup { pos: 5 })]
#[case::ring_order_conflict_digit(b"C=1.CC#1", ParseError::MismatchedRingBondOrders { pos: 7, open_pos: 2 })]
#[case::ring_order_conflict_percent(b"C=%12.CC#%12", ParseError::MismatchedRingBondOrders { pos: 9, open_pos: 2 })]
#[case::ring_dir_conflict_digit(b"C/1.CC\\1", ParseError::MismatchedRingBondDirs { pos: 7, open_pos: 2 })]
#[case::ring_dir_conflict_percent(b"C/%12.CC\\%12", ParseError::MismatchedRingBondDirs { pos: 9, open_pos: 2 })]
#[case::ring_dir_conflict_aromatic(b"c/1.cc\\1", ParseError::MismatchedRingBondDirs { pos: 7, open_pos: 2 })]
#[case::group_dot_before_ring_digit(b"(.1)", ParseError::LeadingDot { pos: 1 })]
#[case::group_dot_before_ring_percent(b"(.%12)", ParseError::LeadingDot { pos: 1 })]
#[case::branch_dot_before_ring_digit(b"C(.1)", ParseError::DotBeforeRing { pos: 2 })]
#[case::branch_dot_before_ring_percent(b"C(.%12)", ParseError::DotBeforeRing { pos: 2 })]
#[case::group_dot_before_bond(b"(.-C)", ParseError::LeadingDot { pos: 1 })]
#[case::branch_dot_before_bond(b"C(.-C)", ParseError::LeadingBond { pos: 3 })]
#[case::leading_bond_after_dot_1(b"C.-C", ParseError::LeadingBond { pos: 2 })]
#[case::leading_bond_after_dot_2(b"C.=-C", ParseError::LeadingBond { pos: 2 })]
#[case::leading_stereobond_after_dot_up(b"C./C", ParseError::LeadingBond { pos: 2 })]
#[case::leading_stereobond_after_dot_down(b"C.\\C", ParseError::LeadingBond { pos: 2 })]
#[case::trailing_bond_dot_aromatic(b"C:.", ParseError::TrailingBond { pos: 1 })]
#[case::group_trailing_bond_dot(b"(C-.)", ParseError::TrailingBond { pos: 2 })]
#[case::branch_trailing_bond_dot(b"C(C-.)", ParseError::TrailingBond { pos: 3 })]
fn components_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_err(), "{:?} should have failed", input_str);
    let err = res.unwrap_err();
    assert_eq!(err, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::atom_c(b"[C]", Element::C, false, None, None, None, None, None)]
#[case::atom_h(b"[H]", Element::H, false, None, None, None, None, None)]
#[case::atom_zn(b"[Zn]", Element::Zn, false, None, None, None, None, None)]
#[case::atom_og(b"[Og]", Element::Og, false, None, None, None, None, None)]
#[case::atom_aromatic_c(b"[c]", Element::C, true, None, None, None, None, None)]
#[case::atom_aromatic_se(b"[se]", Element::Se, true, None, None, None, None, None)]
#[case::atom_aromatic_as(b"[as]", Element::As, true, None, None, None, None, None)]
#[case::isotope_element(b"[13C]", Element::C, false, Some(13), None, None, None, None)]
#[case::isotope_zero(b"[0C]", Element::C, false, Some(0), None, None, None, None)]
#[case::isotope_hydrogen_1h(b"[1H]", Element::H, false, Some(1), None, None, None, None)]
#[case::isotope_hydrogen_2h(b"[2H]", Element::H, false, Some(2), None, None, None, None)]
#[case::isotope_hydrogen_3h(b"[3H]", Element::H, false, Some(3), None, None, None, None)]
#[case::isotope_zero_prefix_1(b"[02H]", Element::H, false, Some(2), None, None, None, None)]
#[case::isotope_zero_prefix_2(b"[002H]", Element::H, false, Some(2), None, None, None, None)]
#[case::isotope_three_digits_1(b"[238U]", Element::U, false, Some(238), None, None, None, None)]
#[case::isotope_three_digits_2(b"[208Pb]", Element::Pb, false, Some(208), None, None, None, None)]
#[case::isotope_unstable(b"[36Cl]", Element::Cl, false, Some(36), None, None, None, None)]
#[case::isotope_max_999(b"[999Og]", Element::Og, false, Some(999), None, None, None, None)]
#[case::isotope_hcount(b"[13CH4]", Element::C, false, Some(13), None, Some(4), None, None)]
#[case::isotope_charge(b"[2H+]", Element::H, false, Some(2), None, None, Some(1), None)]
#[case::chirality_cw(b"[C@]", Element::C, false, None, Some(Chirality::Clockwise), None, None, None)]
#[case::chirality_ccw(b"[C@@]", Element::C, false, None, Some(Chirality::CounterClockwise), None, None, None)]
#[case::chirality_th2(b"[C@TH2]", Element::C, false, None, Some(Chirality::Tetrahedral { arr: 2 }), None, None, None)]
#[case::chirality_al1(b"[C@AL1]", Element::C, false, None, Some(Chirality::Allenal { arr: 1 }), None, None, None)]
#[case::chirality_sp3(b"[C@SP3]", Element::C, false, None, Some(Chirality::SquarePlanar { arr: 3 }), None, None, None)]
#[case::chirality_tb5(b"[C@TB5]", Element::C, false, None, Some(Chirality::TrigonalBipyramidal { arr: 5 }), None, None, None)]
#[case::chirality_oh7(b"[C@OH7]", Element::C, false, None, Some(Chirality::Octahedral { arr: 7 }), None, None, None)]
#[case::hcount_default(b"[CH]", Element::C, false, None, None, Some(1), None, None)]
#[case::hcount_h1(b"[CH1]", Element::C, false, None, None, Some(1), None, None)]
#[case::hcount_h0(b"[CH0]", Element::C, false, None, None, Some(0), None, None)]
#[case::hcount_h3(b"[CH3]", Element::C, false, None, None, Some(3), None, None)]
#[case::hcount_h4(b"[CH4]", Element::C, false, None, None, Some(4), None, None)]
#[case::hcount_aromatic(b"[cH]", Element::C, true, None, None, Some(1), None, None)]
#[case::hcount_two_characters_1(b"[ClH]", Element::Cl, false, None, None, Some(1), None, None)]
#[case::hcount_two_character_2(b"[ClH1]", Element::Cl, false, None, None, Some(1), None, None)]
#[case::chirality_cw_hydrogen(b"[C@H]", Element::C, false, None, Some(Chirality::Clockwise), Some(1), None, None)]
#[case::chirality_ccw_hydrogen(b"[C@@H]", Element::C, false, None, Some(Chirality::CounterClockwise), Some(1), None, None)]
#[case::charge_plus(b"[C+]", Element::C, false, None, None, None, Some(1), None)]
#[case::charge_minus(b"[C-]", Element::C, false, None, None, None, Some(-1), None)]
#[case::charge_pp(b"[C++]", Element::C, false, None, None, None, Some(2), None)]
#[case::charge_mm(b"[C--]", Element::C, false, None, None, None, Some(-2), None)]
#[case::zero_charge_pos(b"[C+0]", Element::C, false, None, None, None, Some(0), None)]
#[case::zero_charge_neg(b"[C-0]", Element::C, false, None, None, None, Some(0), None)]
#[case::charge_plus_15(b"[C+15]", Element::C, false, None, None, None, Some(15), None)]
#[case::charge_minus_15(b"[C-15]", Element::C, false, None, None, None, Some(-15), None)]
#[case::charge_two_characters_plus_1(b"[Na+]", Element::Na, false, None, None, None, Some(1), None)]
#[case::charge_two_characters_plus_2(b"[Ca+2]", Element::Ca, false, None, None, None, Some(2), None)]
#[case::charge_two_characters_pp(b"[Ca++]", Element::Ca, false, None, None, None, Some(2), None)]
#[case::charge_two_characters_minus_1(b"[Cl-]", Element::Cl, false, None, None, None, Some(-1), None)]
#[case::charge_two_characters_minus_2(b"[Se-2]", Element::Se, false, None, None, None, Some(-2), None)]
#[case::charge_two_characters_mm(b"[Se--]", Element::Se, false, None, None, None, Some(-2), None)]
#[case::charge_plus_hcount(b"[C+H]", Element::C, false, None, None, Some(1), Some(1), None)]
#[case::charge_plus_1_hcount(b"[C+1H]", Element::C, false, None, None, Some(1), Some(1), None)]
#[case::charge_minus_hcount(b"[C-H]", Element::C, false, None, None, Some(1), Some(-1), None)]
#[case::charge_minus_1_hcount(b"[C-1H]", Element::C, false, None, None, Some(1), Some(-1), None)]
#[case::hcount_charge_pos_1(b"[NH+]", Element::N, false, None, None, Some(1), Some(1), None)]
#[case::hcount_charge_pos_2(b"[NH+1]", Element::N, false, None, None, Some(1), Some(1), None)]
#[case::hcount_charge_pos_two_characters_1(b"[NaH+]", Element::Na, false, None, None, Some(1), Some(1), None)]
#[case::hcount_charge_pos_two_characters_2(b"[AlH+2]", Element::Al, false, None, None, Some(1), Some(2), None)]
#[case::hcount_charge_pos_two_characters_pp(b"[AlH++]", Element::Al, false, None, None, Some(1), Some(2), None)]
#[case::hcount_charge_neg_1(b"[NH-]", Element::N, false, None, None, Some(1), Some(-1), None)]
#[case::hcount_charge_neg_2(b"[NH-1]", Element::N, false, None, None, Some(1), Some(-1), None)]
#[case::hcount_charge_neg_3(b"[N-H1]", Element::N, false, None, None, Some(1), Some(-1), None)]
#[case::hcount_charge_neg_two_characters_1(b"[AsH-]", Element::As, false, None, None, Some(1), Some(-1), None)]
#[case::hcount_charge_neg_two_characters_2(b"[AsH-2]", Element::As, false, None, None, Some(1), Some(-2), None)]
#[case::hcount_charge_neg_two_characters_mm(b"[AsH--]", Element::As, false, None, None, Some(1), Some(-2), None)]
#[case::class_elem(b"[C:12]", Element::C, false, None, None, None, None, Some(12))]
#[case::class_zero(b"[C:0]", Element::C, false, None, None, None, None, Some(0))]
#[case::class_zero_prefix_1(b"[C:03]", Element::C, false, None, None, None, None, Some(3))]
#[case::class_zero_prefix_2(b"[C:003]", Element::C, false, None, None, None, None, Some(3))]
#[case::class_max_9999(b"[C:9999]", Element::C, false, None, None, None, None, Some(9999))]
#[case::ordering_1(b"[C@H+1:2]", Element::C, false, None, Some(Chirality::Clockwise), Some(1), Some(1), Some(2))]
#[case::ordering_2(b"[CH@+1:2]", Element::C, false, None, Some(Chirality::Clockwise), Some(1), Some(1), Some(2))]
#[case::ordering_3(b"[CH+1@:2]", Element::C, false, None, Some(Chirality::Clockwise), Some(1), Some(1), Some(2))]
#[case::ordering_4(b"[CH+1:2@]", Element::C, false, None, Some(Chirality::Clockwise), Some(1), Some(1), Some(2))]
#[case::ordering_5(b"[C+1@H:2]", Element::C, false, None, Some(Chirality::Clockwise), Some(1), Some(1), Some(2))]
#[case::ordering_6(b"[C:2@H+1]", Element::C, false, None, Some(Chirality::Clockwise), Some(1), Some(1), Some(2))]
fn bracket(
    #[case] input: &[u8],
    #[case] elem: Element,
    #[case] aromatic: bool,
    #[case] isotope: Option<u32>,
    #[case] chirality: Option<Chirality>,
    #[case] hcount: Option<u8>,
    #[case] charge: Option<i8>,
    #[case] class_: Option<u32>,
) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_ok(), "{:?} should have succeeded: {:?}", input_str, res);
    let mol = res.unwrap();
    assert_eq!(mol.atoms.len(), 1, "expected single atom");
    let a = &mol.atoms[0];
    assert!(matches!(a.symbol, AtomSymbol::Element(e) if e == elem));
    assert_eq!(a.aromatic, Some(aromatic));
    assert_eq!(a.isotope_mass, isotope);
    assert_eq!(a.chirality, chirality);
    assert_eq!(a.hydrogens, hcount);
    assert_eq!(a.charge, charge);
    assert_eq!(a.class, class_);
}

#[rstest]
#[case::empty_bracket(b"[]", ParseError::EmptyBracket { pos: 0 })]
#[case::bracket_in_chain_empty(b"C[]", ParseError::EmptyBracket { pos: 1 })]
#[case::bracket_in_group_empty(b"(C[])", ParseError::EmptyBracket { pos: 2 })]
#[case::bracket_in_branch_empty(b"C([])C", ParseError::EmptyBracket { pos: 2 })]
#[case::bracket_in_component_empty(b"[].C", ParseError::EmptyBracket { pos: 0 })]
#[case::bracket_in_ring_empty(b"C1[]C1", ParseError::EmptyBracket { pos: 2 })]
#[case::double_bracket(b"[[C]]", ParseError::InvalidBracket { pos: 1 })]
#[case::invalid_element_1(b"[X]", ParseError::InvalidBracket { pos: 1 })]
#[case::invalid_element_2(b"[Z]", ParseError::InvalidBracket { pos: 1 })]
#[case::invalid_element_3(b"[Aq]", ParseError::InvalidBracket { pos: 1 })]
#[case::invalid_element_4(b"[Sh]", ParseError::InvalidBracket { pos: 2 })]
#[case::invalid_aromatic_element_1(b"[f]", ParseError::InvalidBracket { pos: 1 })]
#[case::invalid_aromatic_element_2(b"[ca]", ParseError::InvalidBracket { pos: 2 })]
#[case::two_elements_1(b"[CF]", ParseError::InvalidBracket { pos: 2 })]
#[case::two_elements_2(b"[AsF]", ParseError::InvalidBracket { pos: 3 })]
#[case::two_elements_3(b"[FAs]", ParseError::InvalidBracket { pos: 2 })]
#[case::two_elements_4(b"[AsBr]", ParseError::InvalidBracket { pos: 3 })]
#[case::zero_charge_no_sign(b"[C0]", ParseError::InvalidBracket { pos: 2 })]
#[case::pos_charge_no_sign(b"[C1]", ParseError::InvalidBracket { pos: 2 })]
#[case::charge_no_element_1(b"[+]", ParseError::InvalidBracket { pos: 1 })]
#[case::charge_no_element_2(b"[-]", ParseError::InvalidBracket { pos: 1 })]
#[case::charge_no_element_3(b"[+0]", ParseError::InvalidBracket { pos: 1 })]
#[case::charge_no_element_4(b"[-0]", ParseError::InvalidBracket { pos: 1 })]
#[case::charge_no_element_5(b"[+1]", ParseError::InvalidBracket { pos: 1 })]
#[case::charge_no_element_6(b"[-1]", ParseError::InvalidBracket { pos: 1 })]
#[case::zero_isotope_no_element(b"[0]", ParseError::InvalidBracket { pos: 2 })]
#[case::isotope_no_element(b"[13]", ParseError::InvalidBracket { pos: 3 })]
#[case::class_no_element(b"[:12]", ParseError::InvalidBracket { pos: 1 })]
#[case::hcount_two_digits_1(b"[CH10]", ParseError::InvalidBracket { pos: 4 })]
#[case::hcount_two_digits_2(b"[SeH10]", ParseError::InvalidBracket { pos: 5 })]
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
#[case::dot_in_bracket(b"[.]", ParseError::InvalidBracket { pos: 1 })]
#[case::branch_open_in_bracket(b"[(]", ParseError::InvalidBracket { pos: 1 })]
#[case::branch_close_in_bracket(b"[)]", ParseError::InvalidBracket { pos: 1 })]
#[case::bracket_in_bracket_1(b"[[]", ParseError::InvalidBracket { pos: 1 })]
#[case::bracket_in_bracket_2(b"[]]", ParseError::EmptyBracket { pos: 0 })]
#[case::open_bracket_in_branch(b"C([)", ParseError::UnbalancedOpenBracket { pos: 2 })]
#[case::close_bracket_in_branch(b"C(])", ParseError::UnbalancedCloseBracket { pos: 2 })]
#[case::unbalanced_close_bracket_1(b"]", ParseError::UnbalancedCloseBracket { pos: 0 })]
#[case::unbalanced_close_bracket_2(b"]C", ParseError::UnbalancedCloseBracket { pos: 0 })]
#[case::unbalanced_close_bracket_3(b"C]", ParseError::UnbalancedCloseBracket { pos: 1 })]
#[case::unbalanced_close_bracket_5(b"C.]", ParseError::UnbalancedCloseBracket { pos: 2 })]
#[case::unbalanced_close_bracket_6(b"].", ParseError::UnbalancedCloseBracket { pos: 0 })]
#[case::unbalanced_close_bracket_7(b"].C", ParseError::UnbalancedCloseBracket { pos: 0 })]
#[case::unbalanced_close_bracket_8(b"(]", ParseError::UnbalancedCloseBracket { pos: 1 })]
#[case::unbalanced_close_bracket_9(b"(C]", ParseError::UnbalancedCloseBracket { pos: 2 })]
#[case::bracket_h_with_hcount_1(b"[HH]", ParseError::BracketHwithHcount { pos: 2 })]
#[case::bracket_h_with_hcount_2(b"[HH1]", ParseError::BracketHwithHcount { pos: 2 })]
#[case::bracket_h_with_hcount_3(b"[HH0]", ParseError::BracketHwithHcount { pos: 2 })]
#[case::duplicate_hcount_1(b"[CHH]", ParseError::DuplicateBracketField { pos: 3 })]
#[case::duplicate_hcount_2(b"[CHH1]", ParseError::DuplicateBracketField { pos: 3 })]
#[case::duplicate_hcount_3(b"[CH1H1]", ParseError::DuplicateBracketField { pos: 4 })]
#[case::duplicate_hcount_4(b"[CH1H]", ParseError::DuplicateBracketField { pos: 4 })]
#[case::duplicate_hcount_5(b"[CH+H]", ParseError::DuplicateBracketField { pos: 4 })]
#[case::duplicate_charge_pos_1(b"[C++1]", ParseError::InvalidBracket { pos: 4 })]
#[case::duplicate_charge_pos_2(b"[C+1+1]", ParseError::DuplicateBracketField { pos: 4 })]
#[case::duplicate_charge_pos_3(b"[C+1+]", ParseError::DuplicateBracketField { pos: 4 })]
#[case::duplicate_charge_pos_4(b"[C+-]", ParseError::DuplicateBracketField { pos: 3 })]
#[case::duplicate_charge_pos_5(b"[C+-1]", ParseError::DuplicateBracketField { pos: 3 })]
#[case::duplicate_charge_pos_6(b"[C+1-1]", ParseError::DuplicateBracketField { pos: 4 })]
#[case::duplicate_charge_pos_7(b"[C+1-]", ParseError::DuplicateBracketField { pos: 4 })]
#[case::duplicate_charge_pos_8(b"[C+H+]", ParseError::DuplicateBracketField { pos: 4 })]
#[case::duplicate_charge_neg_1(b"[C--1]", ParseError::InvalidBracket { pos: 4 })]
#[case::duplicate_charge_neg_2(b"[C-1-1]", ParseError::DuplicateBracketField { pos: 4 })]
#[case::duplicate_charge_neg_3(b"[C-1-]", ParseError::DuplicateBracketField { pos: 4 })]
#[case::duplicate_charge_neg_4(b"[C-+]", ParseError::DuplicateBracketField { pos: 3 })]
#[case::duplicate_charge_neg_5(b"[C-+1]", ParseError::DuplicateBracketField { pos: 3 })]
#[case::duplicate_charge_neg_6(b"[C-1+1]", ParseError::DuplicateBracketField { pos: 4 })]
#[case::duplicate_charge_neg_7(b"[C-1+]", ParseError::DuplicateBracketField { pos: 4 })]
#[case::duplicate_charge_neg_8(b"[C-H-]", ParseError::DuplicateBracketField { pos: 4 })]
#[case::empty_class(b"[C:]", ParseError::MissingClassIndex { pos: 2 })]
#[case::empty_class_two_characters(b"[Cl:]", ParseError::MissingClassIndex { pos: 3 })]
#[case::empty_class_hcount(b"[Cl:H]", ParseError::MissingClassIndex { pos: 3 })]
#[case::empty_class_charge_pos(b"[Na:+]", ParseError::MissingClassIndex { pos: 3 })]
#[case::empty_class_charge_neg(b"[Cl:-]", ParseError::MissingClassIndex { pos: 3 })]
#[case::empty_class_chirality_cw(b"[C:@]", ParseError::MissingClassIndex { pos: 2 })]
#[case::empty_class_chirality_ccw(b"[C:@@]", ParseError::MissingClassIndex { pos: 2 })]
#[case::empty_class_double_colon(b"[C::]", ParseError::MissingClassIndex { pos: 2 })]
#[case::duplicate_class_1(b"[C:1:1]", ParseError::DuplicateBracketField { pos: 4 })]
#[case::duplicate_class_2(b"[C:12:1]", ParseError::DuplicateBracketField { pos: 5 })]
#[case::duplicate_class_3(b"[C:12:12]", ParseError::DuplicateBracketField { pos: 5 })]
#[case::duplicate_class_4(b"[C:1:12]", ParseError::DuplicateBracketField { pos: 4 })]
fn bracket_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_err(), "{:?} should have failed", input_str);
    let err = res.unwrap_err();
    assert_eq!(err, expected);
}

#[rstest]
#[case::stray_chirality_field(b"C@C", ParseError::StrayBracketField { pos: 1 })]
#[case::stray_charge_field(b"C+C", ParseError::StrayBracketField { pos: 1 })]
#[case::stray_hcount_field(b"CHC", ParseError::InvalidElement { pos: 1 })]
#[case::stray_class_field(b"C:1C", ParseError::UnbalancedRingIndex { open_pos: 2 })]
fn bracket_fields_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_err(), "{:?} should have failed", input_str);
    let err = res.unwrap_err();
    assert_eq!(err, expected);
}

#[rstest]
#[case::aromatic_te(b"[te]", Element::Te, true)]
#[case::aromatic_si(b"[si]", Element::Si, true)]
fn bracket_lenient(#[case] input: &[u8], #[case] elem: Element, #[case] aromatic: bool) {
    let res = parse_extended_smiles_bytes_with(input, &SmilesIoConfig::lenient());
    assert!(res.is_ok(), "{:?} should have succeeded", input);
    let mol = res.unwrap();
    assert_eq!(mol.atoms.len(), 1, "expected single atom");
    let a = &mol.atoms[0];
    assert_eq!(a.symbol, AtomSymbol::Element(elem));
    assert_eq!(a.aromatic, Some(aromatic));
}

#[rustfmt::skip]
#[rstest]
#[case::aliphatic_before(b"C[C]", Some(BondOrder::Single), None)]
#[case::aliphatic_before_single(b"C-[C]", Some(BondOrder::Single), None)]
#[case::aliphatic_before_double(b"C=[C]", Some(BondOrder::Double), None)]
#[case::aliphatic_before_triple(b"C#[C]", Some(BondOrder::Triple), None)]
#[case::aliphatic_before_quadruple(b"C$[C]", Some(BondOrder::Quadruple), None)]
#[case::aliphatic_before_aromatic(b"C:[C]", Some(BondOrder::Aromatic), None)]
#[case::aliphatic_before_up(b"C/[C]", Some(BondOrder::Single), Some(BondWedge::Up))]
#[case::aliphatic_before_down(b"C\\[C]", Some(BondOrder::Single), Some(BondWedge::Down))]
#[case::aliphatic_after(b"[C]C", Some(BondOrder::Single), None)]
#[case::aliphatic_after_single(b"[C]-C", Some(BondOrder::Single), None)]
#[case::aliphatic_after_double(b"[C]=C", Some(BondOrder::Double), None)]
#[case::aliphatic_after_triple(b"[C]#C", Some(BondOrder::Triple), None)]
#[case::aliphatic_after_quadruple(b"[C]$C", Some(BondOrder::Quadruple), None)]
#[case::aliphatic_after_aromatic(b"[C]:C", Some(BondOrder::Aromatic), None)]
#[case::aliphatic_after_up(b"[C]/C", Some(BondOrder::Single), Some(BondWedge::Up))]
#[case::aliphatic_after_down(b"[C]\\C", Some(BondOrder::Single), Some(BondWedge::Down))]
#[case::aromatic_before(b"c[c]", Some(BondOrder::Aromatic), None)]
#[case::aromatic_before_single(b"c-[c]", Some(BondOrder::Single), None)]
#[case::aromatic_before_aromatic(b"c:[c]", Some(BondOrder::Aromatic), None)]
#[case::aromatic_after(b"[c]c", Some(BondOrder::Aromatic), None)]
#[case::aromatic_after_single(b"[c]-c", Some(BondOrder::Single), None)]
#[case::aromatic_after_aromatic(b"[c]:c", Some(BondOrder::Aromatic), None)]
#[case::aliphatic_before_aromatic(b"C[c]", Some(BondOrder::Single), None)]
#[case::aliphatic_single_before_aromatic(b"C-[c]", Some(BondOrder::Single), None)]
#[case::aliphatic_aromatic_before_aromatic(b"C:[c]", Some(BondOrder::Aromatic), None)]
#[case::aliphatic_after_aromatic(b"[c]C", Some(BondOrder::Single), None)]
#[case::aromatic_after_aliphatic(b"[C]c", Some(BondOrder::Single), None)]
#[case::aromatic_after_aliphatic_single(b"[C]-c", Some(BondOrder::Single), None)]
#[case::aromatic_after_aliphatic_aromatic(b"[c]:c", Some(BondOrder::Aromatic), None)]
#[case::aromatic_after_aliphatic_up(b"[C]/c", Some(BondOrder::Single), Some(BondWedge::Up))]
#[case::aromatic_after_aliphatic_down(b"[C]\\c", Some(BondOrder::Single), Some(BondWedge::Down))]
#[case::bracket_branch_1(b"[C](C)", Some(BondOrder::Single), None)]
#[case::bracket_branch_2(b"C([C])", Some(BondOrder::Single), None)]
#[case::bracket_branch_single(b"C(-[C])", Some(BondOrder::Single), None)]
#[case::bracket_branch_double(b"C(=[C])", Some(BondOrder::Double), None)]
#[case::bracket_branch_triple(b"C(#[C])", Some(BondOrder::Triple), None)]
#[case::bracket_branch_quadruple(b"C($[C])", Some(BondOrder::Quadruple), None)]
#[case::bracket_branch_aromatic(b"C(:[C])", Some(BondOrder::Aromatic), None)]
#[case::bracket_branch_up(b"C(/[C])", Some(BondOrder::Single), Some(BondWedge::Up))]
#[case::bracket_branch_down(b"C(\\[C])", Some(BondOrder::Single), Some(BondWedge::Down))]
#[case::bracket_group_1(b"([C]C)", Some(BondOrder::Single), None)]
#[case::bracket_group_2(b"(C[C])", Some(BondOrder::Single), None)]
#[case::bracket_ring_1(b"[C]1CC1", Some(BondOrder::Single), None)]
#[case::bracket_ring_2(b"[C]1cc1", Some(BondOrder::Single), None)]
#[case::bracket_ring_double_1(b"[C]1=cc1", Some(BondOrder::Double), None)]
#[case::bracket_ring_double_2(b"[C]=1cc1", Some(BondOrder::Single), None)]
#[case::bracket_aromatic_ring(b"[c]1cc1", Some(BondOrder::Aromatic), None)]
#[case::two_brackets_h2(b"[H][H]", Some(BondOrder::Single), None)]
#[case::two_brackets_hcl(b"[Cl][H]", Some(BondOrder::Single), None)]
#[case::two_brackets_ch4(b"[CH3][H]", Some(BondOrder::Single), None)]
#[case::two_brackets_double_bond(b"[CH2]=[O]", Some(BondOrder::Double), None)]
#[case::two_brackets_triple_bond(b"[C-]#[O+]", Some(BondOrder::Triple), None)]
#[case::two_brackets_quadruple_bond(b"[C]$[C]", Some(BondOrder::Quadruple), None)]
#[case::two_brackets_aromatic_bond(b"[CH]:[CH]", Some(BondOrder::Aromatic), None)]
#[case::two_brackets_up_bond(b"[CH]/[OH]", Some(BondOrder::Single), Some(BondWedge::Up))]
#[case::two_brackets_down_bond(b"[CH]\\[OH]", Some(BondOrder::Single), Some(BondWedge::Down))]
#[case::bracket_before_dot(b"[Na+].[Cl-]", None, None)]
fn bracket_bonds(
    #[case] input: &[u8],
    #[case] expected_order: Option<BondOrder>,
    #[case] expected_dir: Option<BondWedge>,
) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_ok(), "{:?} should have succeeded: {:?}", input_str, res);
    let mol = res.unwrap();
    if let Some(bond1) = mol.bonds.first() {
        if let Some(expected_order) = expected_order {
            assert_eq!(bond1.order, expected_order);
        }
        assert_eq!(bond1.wedge, expected_dir);
    }
}

#[rustfmt::skip]
#[rstest]
// Tetrahedral chirality - equivalent representations
// All these represent the same stereochemistry, center is C at index 1
#[case::chirality_1_eq_1(b"N[C@](Br)(O)C", 1, Element::C, Chirality::Clockwise, vec![0, 2, 3, 4])]
#[case::chirality_1_eq_2(b"N[C@@](Br)(C)O", 1, Element::C, Chirality::CounterClockwise, vec![0, 2, 3, 4])]
#[case::chirality_1_eq_3(b"O[C@](Br)(C)N", 1, Element::C, Chirality::Clockwise, vec![0, 2, 3, 4])]
#[case::chirality_1_eq_4(b"C[C@](Br)(N)O", 1, Element::C, Chirality::Clockwise, vec![0, 2, 3, 4])]
#[case::chirality_1_eq_5(b"C[C@@](Br)(O)N", 1, Element::C, Chirality::CounterClockwise, vec![0, 2, 3, 4])]
#[case::chirality_1_eq_6(b"[C@@](C)(Br)(O)N", 0, Element::C, Chirality::CounterClockwise, vec![1, 2, 3, 4])]
#[case::chirality_1_eq_7(b"Br[C@](O)(N)C", 1, Element::C, Chirality::Clockwise, vec![0, 2, 3, 4])]
#[case::chirality_1_eq_8(b"Br[C@](C)(O)N", 1, Element::C, Chirality::Clockwise, vec![0, 2, 3, 4])]
#[case::chirality_1_eq_9(b"Br[C@](N)(C)O", 1, Element::C, Chirality::Clockwise, vec![0, 2, 3, 4])]
#[case::chiraliry_1_eq_10(b"Br[C@@](N)(O)C", 1, Element::C, Chirality::CounterClockwise, vec![0, 2, 3, 4])]
#[case::chirality_1_eq_11(b"[C@@](Br)(N)(O)C", 0, Element::C, Chirality::CounterClockwise, vec![1, 2, 3, 4])]
// Tetrahedral chirality with ring
#[case::chirality_2_eq_1(b"FC1C[C@](Br)(Cl)CCC1", 3, Element::C, Chirality::Clockwise, vec![2, 4, 5, 6])]
#[case::chirality_2_eq_2(b"[C@]1(Br)(Cl)CCCC(F)C1", 0, Element::C, Chirality::Clockwise, vec![1, 2, 3, 8])]
// Tetrahedral chirality with explicit hydrogen
#[case::chirality_3_eq_1(b"N[C@H](O)C", 1, Element::C, Chirality::Clockwise, vec![0, 2, 3])]
// Allene chirality (2 double bonds - even)
#[case::chirality_allene_1(b"NC(Br)=[C@]=C(O)C", 3, Element::C, Chirality::Clockwise, vec![1, 4])]
#[case::chirality_allene_2(b"NC(Br)=[C@AL1]=C(O)C", 3, Element::C, Chirality::Allenal { arr: 1 }, vec![1, 4])]
// Extended allene chirality (4 double bonds - even)
#[case::chirality_cumulene_4_trans(b"NC(Br)=C=[C@]=C=C(O)C", 4, Element::C, Chirality::Clockwise, vec![3, 5])]
#[case::chirality_cumulene_4_al1(b"NC(Br)=C=[C@AL1]=C=C(O)C", 4, Element::C, Chirality::Allenal { arr: 1 }, vec![3, 5])]
// Extended allene chirality (6 double bonds - even)
#[case::chirality_cumulene_6_trans(b"NC(Br)=C=C=[C@]=C=C=C(O)C", 5, Element::C, Chirality::Clockwise, vec![4, 6])]
#[case::chirality_cumulene_6_al1(b"NC(Br)=C=C=[C@AL1]=C=C=C(O)C", 5, Element::C, Chirality::Allenal { arr: 1 }, vec![4, 6])]
// Trigonal bipyramidal chirality
#[case::chirality_tb_1(b"S[As@TB1](F)(Cl)(Br)N", 1, Element::As, Chirality::TrigonalBipyramidal { arr: 1 }, vec![0, 2, 3, 4, 5])]
#[case::chirality_tb_2(b"S[As@TB5](F)(N)(Cl)Br", 1, Element::As, Chirality::TrigonalBipyramidal { arr: 5 }, vec![0, 2, 3, 4, 5])]
#[case::chirality_tb_3(b"F[As@TB15](Cl)(S)(Br)N", 1, Element::As, Chirality::TrigonalBipyramidal { arr: 15 }, vec![0, 2, 3, 4, 5])]
#[case::chirality_tb_4(b"S[As@TB2](Br)(Cl)(F)N", 1, Element::As, Chirality::TrigonalBipyramidal { arr: 2 }, vec![0, 2, 3, 4, 5])]
#[case::chirality_tb_5(b"F[As@TB10](S)(Cl)(N)Br", 1, Element::As, Chirality::TrigonalBipyramidal { arr: 10 }, vec![0, 2, 3, 4, 5])]
#[case::chirality_tb_6(b"Br[As@TB20](Cl)(S)(F)N", 1, Element::As, Chirality::TrigonalBipyramidal { arr: 20 }, vec![0, 2, 3, 4, 5])]
// Octahedral chirality
#[case::chirality_oh_1(b"C[Co@](F)(Cl)(Br)(I)S", 1, Element::Co, Chirality::Clockwise, vec![0, 2, 3, 4, 5, 6])]
#[case::chirality_oh_2(b"S[Co@OH5](F)(I)(Cl)(C)Br", 1, Element::Co, Chirality::Octahedral { arr: 5 }, vec![0, 2, 3, 4, 5, 6])]
#[case::chirality_oh_3(b"Br[Co@OH12](Cl)(I)(F)(S)C", 1, Element::Co, Chirality::Octahedral { arr: 12 }, vec![0, 2, 3, 4, 5, 6])]
#[case::chirality_oh_4(b"Cl[Co@OH19](C)(I)(F)(S)Br", 1, Element::Co, Chirality::Octahedral { arr: 19 }, vec![0, 2, 3, 4, 5, 6])]
#[case::chirality_oh_5(b"F[Co@@](S)(I)(C)(Cl)Br", 1, Element::Co, Chirality::CounterClockwise, vec![0, 2, 3, 4, 5, 6])]
#[case::chirality_oh_6(b"Br[Co@OH9](C)(S)(Cl)(F)I", 1, Element::Co, Chirality::Octahedral { arr: 9 }, vec![0, 2, 3, 4, 5, 6])]
#[case::chirality_oh_7(b"Cl[Co@OH15](C)(Br)(F)(I)S", 1, Element::Co, Chirality::Octahedral { arr: 15 }, vec![0, 2, 3, 4, 5, 6])]
#[case::chirality_oh_8(b"I[Co@OH27](Cl)(Br)(F)(S)C", 1, Element::Co, Chirality::Octahedral { arr: 27 }, vec![0, 2, 3, 4, 5, 6])]
// Partial stereo specification (first chiral center)
#[case::partial_chirality_1(b"N1[C@H](Cl)[C@@H](Cl)C(Cl)CC1", 1, Element::C, Chirality::Clockwise, vec![0, 2, 3])]
// Single-atom chirality markers (max values)
#[case::chirality_tetrahedral_max2(b"[Ni@TH2]", 0, Element::Ni, Chirality::Tetrahedral { arr: 2 }, vec![])]
#[case::chirality_allenal_max2(b"[C@AL2]", 0, Element::C, Chirality::Allenal { arr: 2 }, vec![])]
#[case::chirality_square_planar_max3(b"[Cu@SP3]", 0, Element::Cu, Chirality::SquarePlanar { arr: 3 }, vec![])]
#[case::chirality_trigonal_bipyramidal_max20(b"[P@TB20]", 0, Element::P, Chirality::TrigonalBipyramidal { arr: 20 }, vec![])]
#[case::chirality_trigonal_bipyramidal_zero_prefix(b"[P@TB02]", 0, Element::P, Chirality::TrigonalBipyramidal { arr: 2 }, vec![])]
#[case::chirality_octahedral_max30(b"[Co@OH30]", 0, Element::Co, Chirality::Octahedral { arr: 30 }, vec![])]
#[case::chirality_octahedral_zero_prefix(b"[Co@OH03]", 0, Element::Co, Chirality::Octahedral { arr: 3 }, vec![])]
fn stereo_chiral(
    #[case] input: &[u8],
    #[case] exp_idx: usize,
    #[case] exp_element: Element,
    #[case] exp_chirality: Chirality,
    #[case] exp_neighbors: Vec<u32>,
) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_ok(), "{:?} should have succeeded: {:?}", input_str, res);
    let mol = res.unwrap();
    let (idx, element, chirality, neighbors) = find_extended_chiral_center(&mol)
        .expect("expected a chiral center");
    assert_eq!(idx, exp_idx, "chiral center index mismatch for {:?}", input_str);
    assert_eq!(element, exp_element, "element mismatch for {:?}", input_str);
    assert_eq!(chirality, exp_chirality, "chirality mismatch for {:?}", input_str);
    assert_eq!(neighbors, exp_neighbors, "neighbors mismatch for {:?}", input_str);
}

#[rstest]
// Double bond stereo - trans
#[case::trans_1_eq_1(b"F/C=C/F", 0, 1, BondWedge::Up)]
#[case::trans_1_eq_2(b"F\\C=C\\F", 0, 1, BondWedge::Down)]
#[case::trans_1_eq_3(b"C(\\F)=C/F", 0, 1, BondWedge::Down)]
// Double bond stereo - cis
#[case::cis_1_eq_1(b"F\\C=C/F", 0, 1, BondWedge::Down)]
#[case::cis_1_eq_2(b"F/C=C\\F", 0, 1, BondWedge::Up)]
#[case::cis_1_eq_3(b"C(/F)=C/F", 0, 1, BondWedge::Up)]
// Cis with substituents
#[case::cis_2_eq_1(b"C/C(/F)=C(\\F)/C", 0, 1, BondWedge::Up)]
#[case::cis_2_eq_2(b"C/C(/F)=C(/C)\\F", 0, 1, BondWedge::Up)]
#[case::cis_2_eq_3(b"C/C(F)=C(/C)F", 0, 1, BondWedge::Up)]
#[case::cis_2_eq_4(b"CC(/F)=C(/C)F", 1, 2, BondWedge::Up)]
#[case::cis_2_eq_5(b"C/C(F)=C(C)\\F", 0, 1, BondWedge::Up)]
#[case::cis_2_eq_6(b"CC(/F)=C(C)\\F", 1, 2, BondWedge::Up)]
// Partial stereo specification
#[case::partial_cis_trans_1(b"F/C=C/C/C=C\\C", 0, 1, BondWedge::Up)]
#[case::partial_cis_trans_2(b"F/C=C/CC=CC", 0, 1, BondWedge::Up)]
// Extended cis/trans for cumulenes (3 double bonds - odd)
#[case::cumulene_3_trans(b"F/C=C=C=C/F", 0, 1, BondWedge::Up)]
#[case::cumulene_3_cis(b"F/C=C=C=C\\F", 0, 1, BondWedge::Up)]
// Extended cis/trans for cumulenes (5 double bonds - odd)
#[case::cumulene_5_trans(b"F/C=C=C=C=C=C/F", 0, 1, BondWedge::Up)]
#[case::cumulene_5_cis(b"F/C=C=C=C=C=C\\F", 0, 1, BondWedge::Up)]
fn stereo_bonds(
    #[case] input: &[u8],
    #[case] exp_a: u32,
    #[case] exp_b: u32,
    #[case] exp_dir: BondWedge,
) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    let (a, b, dir) = find_extended_stereo_bond(&mol).expect("expected a stereo bond");
    assert_eq!(a, exp_a, "atom1 mismatch for {:?}", input_str);
    assert_eq!(b, exp_b, "atom2 mismatch for {:?}", input_str);
    assert_eq!(dir, exp_dir, "direction mismatch for {:?}", input_str);
}

#[rstest]
#[case::chirality_no_element_1(b"[@]", ParseError::InvalidBracket { pos: 1 })]
#[case::chirality_no_element_2(b"[@@]", ParseError::InvalidBracket { pos: 1 })]
#[case::chirality_no_element_3(b"[@TH1]", ParseError::InvalidBracket { pos: 1 })]
#[case::missing_chirality_index_1(b"[C@TH]", ParseError::MissingChiralityIndex { pos: 2 })]
#[case::missing_chirality_index_2(b"[Fe@OH]", ParseError::MissingChiralityIndex { pos: 3 })]
#[case::duplicate_chirality_1(b"[C@@@]", ParseError::DuplicateBracketField { pos: 4 })]
#[case::duplicate_chirality_2(b"[C@TH1@@]", ParseError::DuplicateBracketField { pos: 6 })]
#[case::duplicate_chirality_3(b"[C@TH1@AL1]", ParseError::DuplicateBracketField { pos: 6 })]
#[case::duplicate_chirality_4(b"[C@H@]", ParseError::DuplicateBracketField { pos: 4 })]
#[case::tetrahedral_zero(b"[C@TH0]", ParseError::ChiralityOutOfRange { pos: 2 })]
#[case::tetrahedral_out_of_range(b"[C@TH3]", ParseError::ChiralityOutOfRange { pos: 2 })]
#[case::allenal_zero(b"[C@AL0]", ParseError::ChiralityOutOfRange { pos: 2 })]
#[case::allenal_out_of_range(b"[C@AL3]", ParseError::ChiralityOutOfRange { pos: 2 })]
#[case::square_planar_zero(b"[C@SP0]", ParseError::ChiralityOutOfRange { pos: 2 })]
#[case::square_planar_out_of_range(b"[C@SP4]", ParseError::ChiralityOutOfRange { pos: 2 })]
#[case::trigonal_bipyramidal_zero(b"[C@TB0]", ParseError::ChiralityOutOfRange { pos: 2 })]
#[case::trigonal_bipyramidal_out_of_range(b"[C@TB21]", ParseError::ChiralityOutOfRange { pos: 2 })]
#[case::octahedral_zero(b"[C@OH0]", ParseError::ChiralityOutOfRange { pos: 2 })]
#[case::octahedral_out_of_range(b"[C@OH31]", ParseError::ChiralityOutOfRange { pos: 2 })]
fn stereo_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_err(), "{:?} should have failed", input_str);
    let err = res.unwrap_err();
    assert_eq!(err, expected);
}

#[rstest]
#[case::conflicting_stereo_bonds(
    b"C/C(\\F)=C/FC",
    build_extended_from_graph(
        "C@0 C@2 F@5 C@8 F@10 C@11 | 0-1:/@1 1-2:\\@4 1-3:=@7 3-4:/@9 4-5@11"
    )
)]
fn stereo_invalid_semantics(#[case] input: &[u8], #[case] expected: ExtendedMolecule) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_ok(), "{:?} should have succeeded", input_str);
    let mol = res.unwrap();
    assert_eq!(mol, expected);
}

#[rstest]
#[case::exclamation(b"C!C", ParseError::InvalidToken { pos: 1 })]
#[case::tilde(b"C~C", ParseError::InvalidToken { pos: 1 })]
#[case::underscore(b"C_C", ParseError::InvalidToken { pos: 1 })]
#[case::comma(b"C,C", ParseError::InvalidToken { pos: 1 })]
#[case::semicolon(b"C;C", ParseError::InvalidToken { pos: 1 })]
#[case::less_than(b"C<C", ParseError::InvalidToken { pos: 1 })]
#[case::greater_than(b"C>C", ParseError::InvalidToken { pos: 1 })]
#[case::question(b"C?C", ParseError::InvalidToken { pos: 1 })]
#[case::caret(b"C^C", ParseError::InvalidToken { pos: 1 })]
#[case::backtick(b"C`C", ParseError::InvalidToken { pos: 1 })]
#[case::left_brace(b"C{C", ParseError::InvalidToken { pos: 1 })]
#[case::right_brace(b"C}C", ParseError::InvalidToken { pos: 1 })]
#[case::pipe(b"C|C", ParseError::InvalidToken { pos: 1 })]
#[case::double_quote(b"C\"C", ParseError::InvalidToken { pos: 1 })]
#[case::single_quote(b"C'C", ParseError::InvalidToken { pos: 1 })]
#[case::nul(b"C\x00C", ParseError::InvalidToken { pos: 1 })]
#[case::soh(b"C\x01C", ParseError::InvalidToken { pos: 1 })]
#[case::del(b"C\x7FC", ParseError::InvalidToken { pos: 1 })]
#[case::high_byte(b"C\x80C", ParseError::InvalidToken { pos: 1 })]
#[case::utf8_2byte(b"C\xC3\xA9C", ParseError::InvalidToken { pos: 1 })]
#[case::utf8_3byte(b"C\xE2\x82\xACCC", ParseError::InvalidToken { pos: 1 })]
#[case::en_dash(b"C\xE2\x80\x93C", ParseError::InvalidToken { pos: 1 })]
#[case::em_dash(b"C\xE2\x80\x94C", ParseError::InvalidToken { pos: 1 })]
#[case::cyrillic_es(b"\xD0\xA1C", ParseError::InvalidToken { pos: 0 })]
#[case::greek_omicron(b"C\xCE\xBFC", ParseError::InvalidToken { pos: 1 })]
#[case::greek_capital_omicron(b"C\xCE\x9FC", ParseError::InvalidToken { pos: 1 })]
fn token_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_err(), "{:?} should have failed", input_str);
    let err = res.unwrap_err();
    assert_eq!(err, expected);
}

#[rstest]
#[case::bracketed_organic_atoms(b"[CH3][CH3]", 2, 1)]
#[case::minus_one_charge(b"[CH3-1]", 1, 0)]
#[case::hcount_one(b"C[13CH1](C)C", 4, 3)]
#[case::order_1(b"[C-H3]", 1, 0)]
#[case::order_2(b"C[CH@](Br)Cl", 4, 3)]
#[case::h_bracket_atoms(b"[H][C-]([H])[H]", 4, 3)]
#[case::explicit_single_bond(b"C-C", 2, 1)]
#[case::explicit_aromatic_bond(b"c:1:c:c:c:c:c:1", 6, 6)]
#[case::reused_ring_indices(b"c1ccccc1C1CCCC1", 11, 12)]
#[case::first_ring_not_one(b"c0ccccc0C1CCCC1", 11, 12)]
#[case::non_single_ring_closure(b"CC=1CCCCC=1", 7, 7)]
#[case::adjacent_ring_closures(b"C12(CCCCC1)CCCCC2", 11, 12)]
#[case::zero_prefix_ring_index(b"C%01CCCCC%01", 6, 6)]
#[case::unnecessary_chiral_marker(b"Br[C@H](Br)C", 4, 3)]
#[case::unnecessary_stereo_marker(b"F/C(/F)=C/F", 5, 4)]
#[case::redundant_top_level_parens(b"(N1CCCC1)", 5, 5)]
#[case::aromatic_atoms_in_chain_1(b"CccccC", 6, 5)]
#[case::aromatic_atoms_in_chain_2(b"Ccc", 3, 2)]
#[case::incomplete_stereo_1(b"C/C=C", 3, 2)]
#[case::incomplete_stereo_2(b"C/C=CC", 4, 3)]
fn style_warnings(#[case] input: &[u8], #[case] atoms: usize, #[case] bonds: usize) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_ok(), "{:?} should have succeeded", input_str);
    let mol = res.unwrap();
    assert_eq!(mol.atoms.len(), atoms);
    assert_eq!(mol.bonds.len(), bonds);
}

#[rstest]
#[case::redundant_nested_parens(b"C((C))O", ParseError::EmptyBranch { pos: 5 })]
#[case::consecutive_dots(b"[Na+]..[Cl-]", ParseError::ConsecutiveDots { pos: 5 })]
#[case::leading_dot(b".CCO", ParseError::LeadingDot { pos: 0 })]
#[case::trailing_dot(b"CCO.", ParseError::TrailingDot { pos: 3 })]
#[case::unclosed_ring(b"C1CCC", ParseError::UnbalancedRingIndex { open_pos: 1 })]
#[case::stereo_double_bond(b"CC/=C/C", ParseError::ConsecutiveBonds { pos: 3 })]
#[case::named_isotope_d(b"D[CH3]", ParseError::InvalidElement { pos: 0 })]
#[case::named_isotope_t(b"T[CH3]", ParseError::InvalidElement { pos: 0 })]
fn style_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_err(), "{:?} should have failed", input_str);
    let err = res.unwrap_err();
    assert_eq!(err, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::space(b" ", ExtendedMolecule::empty())]
#[case::tab(b"\t", ExtendedMolecule::empty())]
#[case::newline(b"\n", ExtendedMolecule::empty())]
#[case::cr(b"\r", ExtendedMolecule::empty())]
#[case::crlf(b"\r\n", ExtendedMolecule::empty())]
#[case::terminator_space(b"CC ", build_extended_from_graph("C@0 C@1 | 0-1@1"))]
#[case::terminator_tab(b"CC\t", build_extended_from_graph("C@0 C@1 | 0-1@1"))]
#[case::terminator_newline(b"CC\n", build_extended_from_graph("C@0 C@1 | 0-1@1"))]
#[case::terminator_cr(b"CC\r", build_extended_from_graph("C@0 C@1 | 0-1@1"))]
#[case::terminator_crlf(b"CC\r\n", build_extended_from_graph("C@0 C@1 | 0-1@1"))]
fn whitespace(#[case] input: &[u8], #[case] expected: ExtendedMolecule) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_ok(), "{:?} should have succeeded", input_str);
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
fn whitespace_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_err(), "{:?} should have failed", input_str);
    let err = res.unwrap_err();
    assert_eq!(err, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::wildcard_bare(b"*", build_extended_from_graph("*@0 |"))]
#[case::wildcard_in_chain(b"C*C", build_extended_from_graph("C@0 *@1 C@2 | 0-1@1 1-2@2"))]
#[case::wildcard_branch(b"C(*)C", build_extended_from_graph("C@0 *@2 C@4 | 0-1@2 0-2@4"))]
#[case::wildcard_bonded(b"C-*", build_extended_from_graph("C@0 *@2 | 0-1:-@1"))]
#[case::multiple_wildcards(b"*.*", build_extended_from_graph("*@0 *@2 |"))]
fn wildcard(#[case] input: &[u8], #[case] expected: ExtendedMolecule) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(res.is_ok(), "{:?} should have succeeded: {:?}", input_str, res);
    let mol = res.unwrap();
    assert_eq!(mol, expected);
}

#[rstest]
#[case::wildcard_with_class(b"[*:1]", 1, 0, Some(1))]
#[case::wildcard_with_class_zero(b"[*:0]", 1, 0, Some(0))]
fn wildcard_bracket(
    #[case] input: &[u8],
    #[case] atoms: usize,
    #[case] bonds: usize,
    #[case] class: Option<u32>,
) {
    let res = parse_extended_smiles_bytes(input);
    let input_str = input.to_str_lossy();
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    assert_eq!(mol.atom_count(), atoms);
    assert_eq!(mol.bond_count(), bonds);
    if let Some(expected_class) = class {
        assert_eq!(mol.atoms[0].class, Some(expected_class));
    }
}
