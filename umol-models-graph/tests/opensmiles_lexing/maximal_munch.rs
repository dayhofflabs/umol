//! Maximal munch tests for OpenSMILES (UMOL)

use logos::Logos;
use pretty_assertions::assert_eq;
use rstest::rstest;
use umol_models_graph::io::smiles::lexer::Token;

use super::fixtures::toks;

#[rstest]
#[case("Cl", vec![Token::Cl])]
#[case("Br", vec![Token::Br])]
#[case("se", vec![Token::AromSe])]
#[case("as", vec![Token::AromAs])]
fn multi_char_atoms_take_precedence(
    toks: impl Fn(&str) -> Vec<Token>,
    #[case] input: &str,
    #[case] expected: Vec<Token>,
) {
    assert_eq!(toks(input), expected);
}

#[rstest]
#[case("C", vec![Token::C])]
#[case("B", vec![Token::B])]
#[case("c", vec![Token::AromC])]
fn single_char_when_no_longer_match(
    toks: impl Fn(&str) -> Vec<Token>,
    #[case] input: &str,
    #[case] expected: Vec<Token>,
) {
    assert_eq!(toks(input), expected);
}

#[rstest]
#[case("::", vec![Token::Colon, Token::Colon])]
#[case("//", vec![Token::Slash, Token::Slash])]
#[case("\\\\", vec![Token::Backslash, Token::Backslash])]
fn bond_runs_are_longest_tokens(
    toks: impl Fn(&str) -> Vec<Token>,
    #[case] input: &str,
    #[case] expected: Vec<Token>,
) {
    // Two-character runs that should be two tokens, not one unknown
    assert_eq!(toks(input), expected);
}

#[rstest]
fn percent_vs_digits_is_single_token() {
    // %12 must be one token
    let tokens = Token::lexer(b"%12")
        .map(|t| t.ok())
        .collect::<Option<Vec<_>>>();
    assert_eq!(tokens, Some(vec![Token::Percent(12)]));
}
