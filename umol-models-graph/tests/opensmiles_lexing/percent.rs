//! Percent tests for OpenSMILES (UMOL)

use logos::Logos;
use pretty_assertions::assert_eq;
use rstest::rstest;
use umol_models_graph::io::smiles::lexer::Token;

use super::fixtures::toks;

#[rstest]
fn percent_two_digits_valid(toks: impl Fn(&str) -> Vec<Token>) {
    let got = toks("C%12C");
    assert_eq!(got, vec![Token::C, Token::Percent(12), Token::C]);
}

#[rstest]
fn percent_leading_zero_invalid_tokenized() {
    // Logos will treat "%01" as Error then Digit(1); accept that contract for lexing-level tests
    let mut it = Token::lexer("%01");
    assert!(it.next().unwrap().is_err());
}

#[rstest]
fn percent_zero_invalid_tokenized() {
    // "%0" is not matched by Percent rule; first '%' is Error
    let mut it = Token::lexer("%0");
    assert!(it.next().unwrap().is_err());
}

#[rstest]
fn percent_followed_by_space_splits_tokens() {
    // "%1 2" -> Error, Digit(1), Stop, Digit(2)
    let mut it = Token::lexer("%1 2");
    assert!(it.next().unwrap().is_err());
    let rest = it.map(|t| t.ok()).collect::<Option<Vec<_>>>().unwrap();
    assert_eq!(rest, vec![Token::Digit(1), Token::Stop, Token::Digit(2)]);
}
