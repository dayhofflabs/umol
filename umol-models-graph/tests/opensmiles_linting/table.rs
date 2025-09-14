//! Tests for OpenSMILES (UMOL) linting codes

use rstest::*;
use umol_models_graph::diagnostics::DiagnosticsReport;
use umol_models_graph::io::smiles::linter::{lint_smiles, lint_smiles_parse};

fn codes(report: &DiagnosticsReport) -> Vec<&'static str> {
    report.diagnostics.iter().map(|d| d.code.0).collect()
}

#[rstest]
#[case("C%0", &["LEX_BAD_PERCENT_FORM"])]
#[case("C-", &["SYN_TRAILING_BOND"])]
#[case("C.1", &["SYN_DOT_BEFORE_RING"])]
#[case("[CH1]", &["STYLE_HCOUNT_ONE_SIMPLE"])]
#[case("[C+1]", &["STYLE_CHARGE_SIGN_SIMPLE"])]
#[case("[C]", &["STYLE_BARE_ORGANIC"])]
#[case("C%01", &["STYLE_SINGLE_DIGIT_RING"])]
fn test_style_and_lex_table(#[case] input: &str, #[case] expected: &[&str]) {
    let r = lint_smiles(input);
    let mut got = codes(&r);
    got.sort();
    let mut exp = expected.to_vec();
    let mut expv: Vec<&str> = exp.drain(..).collect();
    expv.sort();
    assert!(
        got.iter().any(|c| expv.contains(c)),
        "codes {:?} do not include any of {:?}",
        got,
        expv
    );
}

#[rstest]
#[case("C1C", &["RING_UNCLOSED"])]
#[case("C11", &["RING_SELF_LOOP"])]
fn test_ring_errors(#[case] input: &str, #[case] expected_any: &[&str]) {
    let r = lint_smiles_parse(input);
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
