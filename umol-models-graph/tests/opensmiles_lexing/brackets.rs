//! Bracket tokenization tests for OpenSMILES (UMOL)

use pretty_assertions::assert_eq;
use rstest::rstest;
use umol_models_graph::io::smiles::lexer::Token;

use super::fixtures::toks;

#[rstest]
#[case("[C]", vec![Token::OpenBracket, Token::C, Token::CloseBracket])]
#[case("[CH3]", vec![Token::OpenBracket, Token::C, Token::H, Token::Digit(3), Token::CloseBracket])]
#[case("[Cl]", vec![Token::OpenBracket, Token::Cl, Token::CloseBracket])]
#[case("[se]", vec![Token::OpenBracket, Token::AromSe, Token::CloseBracket])]
#[case("[as]", vec![Token::OpenBracket, Token::AromAs, Token::CloseBracket])]
#[case("[13C]", vec![Token::OpenBracket, Token::Digit(1), Token::Digit(3), Token::C, Token::CloseBracket])]
#[case("[C:1]", vec![Token::OpenBracket, Token::C, Token::Colon, Token::Digit(1), Token::CloseBracket])]
#[case("[C:12]", vec![Token::OpenBracket, Token::C, Token::Colon, Token::Digit(1), Token::Digit(2), Token::CloseBracket])]
#[case("[C@@]", vec![Token::OpenBracket, Token::C, Token::CounterClockwise, Token::CloseBracket])]
#[case("[C@]", vec![Token::OpenBracket, Token::C, Token::Clockwise, Token::CloseBracket])]
#[case("[C@H+]", vec![Token::OpenBracket, Token::C, Token::Clockwise, Token::H, Token::Plus, Token::CloseBracket])]
#[case("[C+07]", vec![Token::OpenBracket, Token::C, Token::Plus, Token::Digit(0), Token::Digit(7), Token::CloseBracket])]
#[case("[CH10]", vec![Token::OpenBracket, Token::C, Token::H, Token::Digit(1), Token::Digit(0), Token::CloseBracket])]
#[case("[*]", vec![Token::OpenBracket, Token::Asterisk, Token::CloseBracket])]
#[case("[C H]", vec![Token::OpenBracket, Token::C, Token::Stop, Token::H, Token::CloseBracket])]
#[case("[C%12]", vec![Token::OpenBracket, Token::C, Token::Percent(12), Token::CloseBracket])]
fn bracket_various_tokenize(
    toks: impl Fn(&str) -> Vec<Token>,
    #[case] input: &str,
    #[case] expected: Vec<Token>,
) {
    assert_eq!(toks(input), expected);
}

#[rstest]
fn bracket_field_order_violations_emit_errors_at_lex_level_or_tokens(
    toks: impl Fn(&str) -> Vec<Token>,
) {
    // Lexer does not enforce order; invalid sequences still tokenize
    // '[HC]' -> '[', 'H', 'C', ']'
    assert_eq!(
        toks("[HC]"),
        vec![Token::OpenBracket, Token::H, Token::C, Token::CloseBracket]
    );
    // '[C:1:2]' -> '[', 'C', ':', '1', ':', '2', ']'
    assert_eq!(
        toks("[C:1:2]"),
        vec![
            Token::OpenBracket,
            Token::C,
            Token::Colon,
            Token::Digit(1),
            Token::Colon,
            Token::Digit(2),
            Token::CloseBracket,
        ]
    );
}
