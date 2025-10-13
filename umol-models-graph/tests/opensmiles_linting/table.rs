//! Tests for OpenSMILES (UMOL) linting codes

use rstest::*;
use umol_models_graph::io::smiles::diagnostics::DiagnosticList;
use umol_models_graph::io::smiles::linter::lint_smiles;

fn codes(report: &DiagnosticList) -> Vec<&'static str> {
    report.diagnostics.iter().map(|d| d.code.0).collect()
}

#[rstest]
#[case::lex_ring_index_invalid("C%0", &["LEX_RING_INDEX_INVALID"])]
#[case::lex_trailing_bond("C-", &["LEX_TRAILING_BOND"])]
#[case::lex_leading_dot("C.1", &["LEX_LEADING_DOT"])]
#[case::brkt_unbalanced_open("[C", &["BRKT_UNBALANCED_OPEN"])]
#[case::brkt_hcount_one("[CH1]", &[])]
#[case::brkt_charge_plus_one("[C+1]", &[])]
#[case::brkt_unnecessary_bracket("[C]", &[])]
#[case::style_unnecessary_percent_ring_index("C%01", &["STYLE_UNNECESSARY_PERCENT_RING_INDEX"])]
fn test_style_and_lex_table(#[case] input: &str, #[case] expected: &[&str]) {
    let r = lint_smiles(input);
    let mut got = codes(&r);
    got.sort();
    let mut exp = expected.to_vec();
    exp.sort();
    assert!(
        exp.iter().all(|e| got.contains(e)),
        "codes {:?} do not include all of {:?}",
        got,
        exp
    );
}

#[rstest]
#[case::brkt_hcount_one_hh("[HH]", &[])]
#[case::brkt_hcount_one_hh1("[HH1]", &[])]
#[case::brkt_invalid_class("[C:-1]", &[])]
#[case::lex_ring_index_invalid_percent("%", &["LEX_RING_INDEX_INVALID"])]
#[case::lex_ring_index_invalid_percent_one_x("%1x", &["LEX_RING_INDEX_INVALID"])]
#[case::brkt_unbalanced_close("]", &["BRKT_UNBALANCED_CLOSE"])]
#[case::brkt_unexpected_close("]", &["BRKT_UNBALANCED_CLOSE"])]
#[case::brkt_field_outside("+", &["BRKT_FIELD_OUTSIDE"])]
#[case::lex_leading_dot(".", &["LEX_LEADING_DOT"])]
#[case::lex_trailing_dot("C.", &["LEX_TRAILING_DOT"])]
#[case::lex_multiple_dots("C..C", &["LEX_MULTIPLE_DOTS"])]
#[case::brkt_field_outside_x("X", &[])]
#[case::brch_unexpected_close("C)", &["BRCH_UNEXPECTED_CLOSE"])]
#[case::brch_unclosed("(C", &["BRCH_UNCLOSED"])]
#[case::grp_leading_dot("(.C)", &["GRP_LEADING_DOT"])]
#[case::grp_leading_bond("(-C)", &["GRP_LEADING_BOND"])]
#[case::brch_empty_branch("C()", &["BRCH_EMPTY_BRANCH"])]
#[case::lex_trailing_bond("C(-)", &["LEX_TRAILING_BOND"])]
#[case::brch_unclosed("C(C.C", &["BRCH_UNCLOSED"])]
#[case::brch_empty_branch("C(())", &["BRCH_EMPTY_BRANCH"])]
fn test_error_table(#[case] input: &str, #[case] expected_any: &[&str]) {
    let r = lint_smiles(input);
    let mut got = codes(&r);
    got.sort();
    let mut exp = expected_any.to_vec();
    exp.sort();
    assert!(
        exp.iter().all(|e| got.contains(e)),
        "codes {:?} do not include all of {:?}",
        got,
        exp
    );
}

#[rstest]
#[case::ring_unclosed("C1C", &["RING_UNCLOSED"])]
#[case::ring_bond_dir_conflict("C/1CC\\1", &["RING_BOND_DIR_CONFLICT"])]
#[case::ring_bond_order_conflict("C=1CC#1", &["RING_BOND_ORDER_CONFLICT"])]
#[case::ring_unclosed("C1.C", &["RING_UNCLOSED"])]
fn test_ring_errors(#[case] input: &str, #[case] expected_any: &[&str]) {
    let r = lint_smiles(input);
    let mut got = codes(&r);
    got.sort();
    let mut exp = expected_any.to_vec();
    exp.sort();
    assert!(
        got.iter().any(|c| exp.contains(c)),
        "codes {:?} do not include any of {:?}",
        got,
        exp
    );
}
