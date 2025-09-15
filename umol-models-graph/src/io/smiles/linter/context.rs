//! Context for SMILES linting.

use crate::io::smiles::lexer::Lexer;

pub struct LintContext<'a> {
    pub input: &'a str,
    // Lazily available resources as needed later
    pub lexer: Lexer<'a>,
}

impl<'a> LintContext<'a> {
    pub fn new(input: &'a str) -> Self {
        Self {
            input,
            lexer: Lexer::new(input),
        }
    }
}
