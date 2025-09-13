//! Charge tokenization tests for OpenSMILES (UMOL)

use pretty_assertions::assert_eq;
use rstest::rstest;
use umol_models_graph::io::smiles::lexer::Token;

use super::fixtures::toks;

#[rstest]
#[case("+", vec![Token::Plus])]
#[case("-", vec![Token::Minus])]
#[case("++", vec![Token::PlusTwo])]
#[case("--", vec![Token::MinusTwo])]
fn simple_charges_tokenize(
    toks: impl Fn(&str) -> Vec<Token>,
    #[case] input: &str,
    #[case] expected: Vec<Token>,
) {
    assert_eq!(toks(input), expected);
}

#[rstest]
fn plus_or_minus_followed_by_digits(toks: impl Fn(&str) -> Vec<Token>) {
    // Lexer does not combine +/- with digits; it yields separate tokens
    assert_eq!(toks("+0"), vec![Token::Plus, Token::Digit(0)]);
    assert_eq!(
        toks("+07"),
        vec![Token::Plus, Token::Digit(0), Token::Digit(7)]
    );
    assert_eq!(
        toks("-99"),
        vec![Token::Minus, Token::Digit(9), Token::Digit(9)]
    );
}
