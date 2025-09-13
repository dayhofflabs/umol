//! Whitespace tests for OpenSMILES (UMOL)

use pretty_assertions::assert_eq;
use rstest::rstest;
use umol_models_graph::io::smiles::lexer::Token;

use super::fixtures::toks;

#[rstest]
#[case("C C", vec![Token::C, Token::Stop, Token::C])]
#[case("C\tC", vec![Token::C, Token::Stop, Token::C])]
fn inter_token_whitespace_is_stop_tokens(toks: impl Fn(&str) -> Vec<Token>, #[case] input: &str, #[case] expected: Vec<Token>) {
    assert_eq!(toks(input), expected);
}

#[rstest]
fn trailing_whitespace_yields_stop_and_eoi(toks: impl Fn(&str) -> Vec<Token>) {
    // Two trailing whitespace chars yield two Stop tokens
    assert_eq!(toks("C \n"), vec![Token::C, Token::Stop, Token::Stop]);
}

#[rstest]
fn leading_whitespace_is_just_stop(toks: impl Fn(&str) -> Vec<Token>) {
    assert_eq!(toks(" C"), vec![Token::Stop, Token::C]);
}
