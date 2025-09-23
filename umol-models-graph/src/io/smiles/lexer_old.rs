//! Legacy string-backed lexer wrapper for quick decoupling.

use logos::{Logos, SpannedIter};

pub use super::lexer::{LexicalError, Token};

pub type Spanned<Tok, Loc, Error> = Result<(Loc, Tok, Loc), Error>;

pub struct Lexer<'input> { token_stream: SpannedIter<'input, Token> }

impl<'input> Lexer<'input> {
    pub fn new(input: &'input str) -> Self {
        Self { token_stream: Token::lexer(input).spanned() }
    }
}

impl<'input> Iterator for Lexer<'input> {
    type Item = Spanned<Token, usize, LexicalError>;
    fn next(&mut self) -> Option<Self::Item> {
        self.token_stream.next().map(|(tok, span)| match tok {
            Ok(token) => Ok((span.start, token, span.end)),
            Err(_) => Ok((span.start, Token::Error, span.end)),
        })
    }
}


