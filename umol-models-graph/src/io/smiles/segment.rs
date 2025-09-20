//! SMILES segmentation: iterator over structural segments (supertokens).
//!
//! This groups bracket atoms, compresses whitespace, and classifies tokens into
//! atoms, bonds, branches, ring closures, component separators, and errors.

use crate::diagnostics::Span;
use crate::io::smiles::lexer::{Lexer, Token};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondKind {
    Single,
    Double,
    Triple,
    Quadruple,
    Aromatic,
    Up,
    Down,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Segment<'a> {
    /// Organic/aromatic atom or `*` outside brackets.
    AtomSimple { span: Span, raw: &'a str },
    /// Bracket atom (well-formed): `[ ... ]` including the brackets.
    AtomBracket { span: Span, inner: &'a str },
    /// Malformed bracket atom (missing closing `]`).
    MalformedBracket { span: Span, inner: &'a str },
    /// Bond symbol outside brackets.
    Bond { span: Span, kind: BondKind },
    /// `(` outside brackets.
    BranchOpen { span: Span },
    /// `)` outside brackets.
    BranchClose { span: Span },
    /// Ring closure outside brackets: digit or percent-two-digit.
    RingClosure { span: Span, index: u32 },
    /// Component separator `.` outside brackets.
    ComponentSeparator { span: Span },
    /// One or more whitespace characters.
    WhitespaceBlock { span: Span },
    /// Lexical error token from the lexer.
    LexError { span: Span },
    /// Token that is only valid inside brackets, but appears outside (e.g., '@', '++').
    StrayBracketField { span: Span },
    /// Syntactically invalid token in this context (e.g., ']' or '@' outside brackets).
    Invalid { span: Span },
}

pub struct Segments<'input> {
    input: &'input str,
    lexer: Lexer<'input>,
    peeked: Option<(usize, Token, usize)>,
}

impl<'input> Segments<'input> {
    pub fn new(input: &'input str) -> Self {
        Self {
            input,
            lexer: Lexer::new(input),
            peeked: None,
        }
    }

    fn next_tok(&mut self) -> Option<(usize, Token, usize)> {
        if let Some(t) = self.peeked.take() {
            return Some(t);
        }
        self.lexer.next().and_then(|r| r.ok())
    }

    fn peek_tok(&mut self) -> Option<(usize, Token, usize)> {
        if self.peeked.is_none() {
            self.peeked = self.lexer.next().and_then(|r| r.ok());
        }
        self.peeked.clone()
    }
}

impl<'input> Iterator for Segments<'input> {
    type Item = Segment<'input>;

    fn next(&mut self) -> Option<Self::Item> {
        let (l, tok, r) = self.next_tok()?;

        // Compress whitespace blocks
        if matches!(tok, Token::Stop) {
            let start = l;
            let mut end = r;
            while let Some((_, nt, nr)) = self.peek_tok() {
                if matches!(nt, Token::Stop) {
                    // consume
                    let _ = self.next_tok();
                    end = nr;
                } else {
                    break;
                }
            }
            return Some(Segment::WhitespaceBlock {
                span: Span::new(start, end),
            });
        }

        match tok {
            Token::OpenBracket => {
                let start = l;
                let inner_start = r; // after '['
                let mut end = r;
                let mut closed = false;
                while let Some((_, nt, nr)) = self.next_tok() {
                    end = nr;
                    if matches!(nt, Token::CloseBracket) {
                        closed = true;
                        break;
                    }
                }
                if closed {
                    let span = Span::new(start, end);
                    // Safe slicing: inner range is within input and excludes ']'
                    let inner_end = end.saturating_sub(1);
                    let inner = if inner_start <= inner_end && inner_end <= self.input.len() {
                        &self.input[inner_start..inner_end]
                    } else {
                        ""
                    };
                    Some(Segment::AtomBracket { span, inner })
                } else {
                    let span = Span::new(start, self.input.len());
                    let inner = if inner_start <= self.input.len() {
                        &self.input[inner_start..]
                    } else {
                        ""
                    };
                    Some(Segment::MalformedBracket { span, inner })
                }
            }

            Token::CloseBracket => Some(Segment::Invalid {
                span: Span::new(l, r),
            }),

            Token::OpenParen => Some(Segment::BranchOpen {
                span: Span::new(l, r),
            }),
            Token::CloseParen => Some(Segment::BranchClose {
                span: Span::new(l, r),
            }),

            Token::Digit(v) => Some(Segment::RingClosure {
                span: Span::new(l, r),
                index: v,
            }),
            Token::Percent(v) => Some(Segment::RingClosure {
                span: Span::new(l, r),
                index: v,
            }),

            Token::Dot => Some(Segment::ComponentSeparator {
                span: Span::new(l, r),
            }),

            Token::Minus => Some(Segment::Bond {
                span: Span::new(l, r),
                kind: BondKind::Single,
            }),
            Token::Equal => Some(Segment::Bond {
                span: Span::new(l, r),
                kind: BondKind::Double,
            }),
            Token::Hash => Some(Segment::Bond {
                span: Span::new(l, r),
                kind: BondKind::Triple,
            }),
            Token::Dollar => Some(Segment::Bond {
                span: Span::new(l, r),
                kind: BondKind::Quadruple,
            }),
            Token::Colon => Some(Segment::Bond {
                span: Span::new(l, r),
                kind: BondKind::Aromatic,
            }),
            Token::Slash => Some(Segment::Bond {
                span: Span::new(l, r),
                kind: BondKind::Up,
            }),
            Token::Backslash => Some(Segment::Bond {
                span: Span::new(l, r),
                kind: BondKind::Down,
            }),

            // Outside brackets, any alphabetic token (elements, aromatic atoms) or '*' is a simple atom.
            other => {
                let raw = &self.input[l..r];
                let is_alpha = raw.as_bytes().iter().all(|b| b.is_ascii_alphabetic());
                if is_alpha || raw == "*" {
                    Some(Segment::AtomSimple {
                        span: Span::new(l, r),
                        raw,
                    })
                } else {
                    match other {
                        // Tokens that should only appear inside brackets => mark as stray bracket field
                        Token::Plus
                        | Token::PlusTwo
                        | Token::MinusTwo
                        | Token::Clockwise
                        | Token::CounterClockwise
                        | Token::Tetrahedral
                        | Token::Allenal
                        | Token::SquarePlanar
                        | Token::TrigonalBipyramidal
                        | Token::Octahedral => Some(Segment::StrayBracketField {
                            span: Span::new(l, r),
                        }),
                        Token::Error => Some(Segment::LexError {
                            span: Span::new(l, r),
                        }),
                        // Fallback: treat unknown as error to avoid misclassifying.
                        _ => Some(Segment::Invalid {
                            span: Span::new(l, r),
                        }),
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_segments() {
        let input = "C.1[CH3]C(=O)N";
        let segs = Segments::new(input).collect::<Vec<_>>();

        // Expect: C, ., 1, [CH3], C, (, =, O, ), N
        assert!(matches!(segs[0], Segment::AtomSimple { .. }));
        assert!(matches!(segs[1], Segment::ComponentSeparator { .. }));
        assert!(matches!(segs[2], Segment::RingClosure { index: 1, .. }));
        assert!(matches!(segs[3], Segment::AtomBracket { .. }));
        assert!(matches!(segs[4], Segment::AtomSimple { .. }));
        assert!(matches!(segs[5], Segment::BranchOpen { .. }));
        assert!(matches!(
            segs[6],
            Segment::Bond {
                kind: BondKind::Double,
                ..
            }
        ));
        assert!(matches!(segs[7], Segment::AtomSimple { .. }));
        assert!(matches!(segs[8], Segment::BranchClose { .. }));
        assert!(matches!(segs[9], Segment::AtomSimple { .. }));
    }

    #[test]
    fn test_malformed_bracket() {
        let input = "C[CH3"; // missing closing bracket
        let segs = Segments::new(input).collect::<Vec<_>>();
        assert!(matches!(segs[0], Segment::AtomSimple { .. }));
        assert!(matches!(segs[1], Segment::MalformedBracket { .. }));
    }

    #[test]
    fn test_stray_close_bracket() {
        let input = "]C";
        let segs = Segments::new(input).collect::<Vec<_>>();
        assert!(matches!(segs[0], Segment::Invalid { .. }));
        assert!(matches!(segs[1], Segment::AtomSimple { .. }));
    }

    #[test]
    fn test_percent_two_digit() {
        let input = "C%12C";
        let segs = Segments::new(input).collect::<Vec<_>>();
        assert!(matches!(segs[0], Segment::AtomSimple { .. }));
        assert!(matches!(segs[1], Segment::RingClosure { index: 12, .. }));
        assert!(matches!(segs[2], Segment::AtomSimple { .. }));
    }

    #[test]
    fn test_stray_chirality() {
        let input = "C@C";
        let segs = Segments::new(input).collect::<Vec<_>>();
        assert!(matches!(segs[0], Segment::AtomSimple { .. }));
        assert!(matches!(segs[1], Segment::StrayBracketField { .. }));
        assert!(matches!(segs[2], Segment::AtomSimple { .. }));
    }
}
