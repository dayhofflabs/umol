//! Punctuation (bond) tokenization tests for OpenSMILES (UMOL)

use pretty_assertions::assert_eq;
use rstest::rstest;
use umol_models_graph::io::smiles::lexer::Token;

use super::fixtures::toks;

#[rstest]
#[case("-", vec![Token::Minus])]
#[case("=", vec![Token::Equal])]
#[case("#", vec![Token::Hash])]
#[case("$", vec![Token::Dollar])]
#[case(":", vec![Token::Colon])]
#[case("/", vec![Token::Slash])]
#[case("\\", vec![Token::Backslash])]
fn bond_symbols_tokenize(
    toks: impl Fn(&str) -> Vec<Token>,
    #[case] input: &str,
    #[case] expected: Vec<Token>,
) {
    assert_eq!(toks(input), expected);
}

#[rstest]
#[case("==", vec![Token::Equal, Token::Equal])]
#[case("##", vec![Token::Hash, Token::Hash])]
#[case("::", vec![Token::Colon, Token::Colon])]
fn repeated_bond_symbols_are_multiple_tokens(
    toks: impl Fn(&str) -> Vec<Token>,
    #[case] input: &str,
    #[case] expected: Vec<Token>,
) {
    assert_eq!(toks(input), expected);
}

#[rstest]
fn repeated_dots_are_multiple_tokens(toks: impl Fn(&str) -> Vec<Token>) {
    assert_eq!(toks(".."), vec![Token::Dot, Token::Dot]);
}
