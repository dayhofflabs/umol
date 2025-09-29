//! Tests for OpenSMILES (UMOL) linting codes

use rstest::*;
use umol_models_graph::diagnostics::DiagnosticsReport;
use umol_models_graph::io::smiles::linter::lint_smiles;

fn codes(report: &DiagnosticsReport) -> Vec<&'static str> {
    report.diagnostics.iter().map(|d| d.code.0).collect()
}

#[rstest]
#[case("C%0", &["LEX_RING_INDEX_INVALID"])]
#[case("C-", &["LEX_TRAILING_BOND"])]
#[case("C.1", &["LEX_DOT_BEFORE_RING"])]
#[case("[C", &["PARSER_UNBALANCED_OPEN_BRACKET"])]
#[case("[CH1]", &["STYLE_HCOUNT_ONE_SIMPLE"])]
#[case("[C+1]", &["STYLE_CHARGE_SIGN_SIMPLE"])]
#[case("[C]", &["STYLE_BRKT_ORGANIC"])]
#[case("C%01", &["STYLE_UNNECESSARY_PERCENT_RING_INDEX"])]
fn test_style_and_lex_table(#[case] input: &str, #[case] expected: &[&str]) {
    println!("input: {:?}", input);
    let r = lint_smiles(input);
    let mut got = codes(&r);
    got.sort();
    let mut exp = expected.to_vec();
    let mut expv: Vec<&str> = exp.drain(..).collect();
    expv.sort();
    println!("got: {:?}", got);
    println!("expv: {:?}", expv);
    assert!(expv.iter().all(|e| got.contains(e)), "codes {:?} do not include all of {:?}", got, expv);
}

#[rstest]
#[case("[HH]", &[])]
#[case("[HH1]", &[])]
#[case("[C:-1]", &[])]
#[case("%", &["LEX_RING_INDEX_INVALID"])]
#[case("%1x", &["LEX_RING_INDEX_INVALID"])]
#[case("]", &["PARSER_UNBALANCED_CLOSE_BRACKET"])]
#[case("@", &[])]
#[case("+", &[])]
#[case(".", &["LEX_LEADING_DOT"])]
#[case("C.", &["LEX_TRAILING_DOT"])]
#[case("C..C", &["LEX_MULTIPLE_DOTS"])]
#[case("X", &[])]
#[case("C)", &["BRCH_UNEXPECTED_CLOSE"])]
#[case("(C", &["BRCH_UNCLOSED"])]
#[case("C()", &["BRCH_EMPTY_BRANCH"])]
#[case("C(-)", &["BRCH_DANGLING_BOND"])]
#[case("C(C.C", &["BRCH_UNCLOSED"])]
#[case("C(())", &["BRCH_EMPTY_BRANCH"])]
fn test_error_table(#[case] input: &str, #[case] expected_any: &[&str]) {
    println!("input: {:?}", input);
    let r = lint_smiles(input);
    let mut got = codes(&r);
    got.sort();
    let mut exp = expected_any.to_vec();
    let mut expv: Vec<&str> = exp.drain(..).collect();
    expv.sort();
    println!("got: {:?}", got);
    println!("expv: {:?}", expv);
    assert!(expv.iter().all(|e| got.contains(e)), "codes {:?} do not include all of {:?}", got, expv);
}

#[rstest]
#[case("C1C", &["RING_UNCLOSED"])]
#[case("C11", &["RING_SELF_LOOP"])]
#[case("C/1CC\\1", &["RING_BOND_DIR_CONFLICT"])]
#[case("C=1CC#1", &["RING_BOND_ORDER_CONFLICT"])]
#[case("C1.C", &["RING_UNCLOSED"])]
fn test_ring_errors(#[case] input: &str, #[case] expected_any: &[&str]) {
    let r = lint_smiles(input);
    let mut got = codes(&r);
    got.sort();
    let mut exp = expected_any.to_vec();
    let mut expv: Vec<&str> = exp.drain(..).collect();
    expv.sort();
    assert!(
        got.iter().any(|c| expv.contains(c)),
        "codes {:?} do not include any of {:?}",
        got,
        expv
    );
}
