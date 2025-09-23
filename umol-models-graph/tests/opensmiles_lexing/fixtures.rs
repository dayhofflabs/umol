//! Fixtures for OpenSMILES (UMOL) lexing tests

use logos::Logos;
use rstest::fixture;
use umol_models_graph::io::smiles::lexer::Token;

#[fixture]
pub fn toks() -> impl Fn(&str) -> Vec<Token> {
    |input: &str| Token::lexer(input.as_bytes())
        .map(|t| t.ok())
        .collect::<Option<Vec<_>>>()
        .unwrap_or_default()
}
