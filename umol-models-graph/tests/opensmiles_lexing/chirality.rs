//! Chirality tokenization tests for OpenSMILES (UMOL)

use pretty_assertions::assert_eq;
use rstest::rstest;
use umol_models_graph::io::smiles::lexer::Token;

use super::fixtures::toks;

#[rstest]
#[case("@", vec![Token::Clockwise])]
#[case("@@", vec![Token::CounterClockwise])]
#[case("@TH", vec![Token::Tetrahedral])]
#[case("@AL", vec![Token::Allenal])]
#[case("@SP", vec![Token::SquarePlanar])]
#[case("@TB", vec![Token::TrigonalBipyramidal])]
#[case("@OH", vec![Token::Octahedral])]
fn basic_chirality_tokens(
    toks: impl Fn(&str) -> Vec<Token>,
    #[case] input: &str,
    #[case] expected: Vec<Token>,
) {
    assert_eq!(toks(input), expected);
}

#[rstest]
fn chirality_with_numbers_tokenize_as_separate_tokens(toks: impl Fn(&str) -> Vec<Token>) {
    // '@TH1' -> Tetrahedral, Digit(1)
    assert_eq!(toks("@TH1"), vec![Token::Tetrahedral, Token::Digit(1)]);
    assert_eq!(
        toks("@OH30"),
        vec![Token::Octahedral, Token::Digit(3), Token::Digit(0)]
    );
}
