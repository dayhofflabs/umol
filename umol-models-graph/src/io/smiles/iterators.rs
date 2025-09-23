//! SMILES segmentation: iterator over structural segments (supertokens).
//!
//! This groups bracket atoms, compresses whitespace, and classifies tokens into
//! atoms, bonds, branches, ring closures, component separators, and errors.

use crate::diagnostics::Span;
use crate::io::smiles::lexer_old::{Lexer, Token};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
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
    // TODO: Review naming MalformedBracket -> UnclosedBracket/SoleBracketOpen,
    // TODO: BracketClose -> SoleBracketClose, FreeBracketField -> UnenclosedBracketField,
    // TODO: Verify the use of Invalid, the semantics are not clear. ] or @ outside brackets
    // TODO: have separate segments
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
    NewComponent { span: Span },
    /// One or more whitespace characters.
    WhitespaceBlock { span: Span },
    /// Lexical error token from the lexer.
    LexError { span: Span },
    /// Token that is only valid inside brackets, but appears outside (e.g., '@', '++').
    FreeBracketField { span: Span },
    /// Closing ']' that appears outside of a bracket atom.
    BracketClose { span: Span },
    /// Syntactically invalid token in this context
    Invalid { span: Span },
}

pub struct Segments<'input> {
    input: &'input str,
    lexer: Lexer<'input>,
    peeked: Option<(usize, Token, usize)>,
}

impl<'input> Segments<'input> {
    pub fn new(input: &'input str) -> Self {
        Self { input, lexer: Lexer::new(input), peeked: None }
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
        if matches!(tok, Token::Stop) {
            let start = l;
            let mut end = r;
            while let Some((_, nt, nr)) = self.peek_tok() {
                if matches!(nt, Token::Stop) {
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
                let inner_start = r;
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
            Token::CloseBracket => Some(Segment::BracketClose {
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
            Token::Dot => Some(Segment::NewComponent {
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
                        Token::Plus
                        | Token::PlusTwo
                        | Token::MinusTwo
                        | Token::Clockwise
                        | Token::CounterClockwise
                        | Token::Tetrahedral
                        | Token::Allenal
                        | Token::SquarePlanar
                        | Token::TrigonalBipyramidal
                        | Token::Octahedral => Some(Segment::FreeBracketField {
                            span: Span::new(l, r),
                        }),
                        Token::Error => Some(Segment::LexError {
                            span: Span::new(l, r),
                        }),
                        _ => Some(Segment::Invalid {
                            span: Span::new(l, r),
                        }),
                    }
                }
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BranchEventKind {
    Open,
    Close,
    NewComponent,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BranchEvent {
    pub kind: BranchEventKind,
    pub span: Span,
    pub depth_after: usize,
}

pub struct Branches<'a> {
    segs: &'a [Segment<'a>],
    idx: usize,
    depth: usize,
}

impl<'a> Branches<'a> {
    pub fn new(segs: &'a [Segment<'a>]) -> Self {
        Self {
            segs,
            idx: 0,
            depth: 0,
        }
    }
}

impl<'a> Iterator for Branches<'a> {
    type Item = BranchEvent;
    fn next(&mut self) -> Option<Self::Item> {
        while self.idx < self.segs.len() {
            let ev = match self.segs[self.idx] {
                Segment::BranchOpen { span } => {
                    self.depth = self.depth.saturating_add(1);
                    Some(BranchEvent {
                        kind: BranchEventKind::Open,
                        span,
                        depth_after: self.depth,
                    })
                }
                Segment::BranchClose { span } => {
                    let depth_now = self.depth;
                    self.depth = self.depth.saturating_sub(1);
                    Some(BranchEvent {
                        kind: BranchEventKind::Close,
                        span,
                        depth_after: depth_now.saturating_sub(0),
                    })
                }
                Segment::NewComponent { span } => {
                    if self.depth > 0 {
                        Some(BranchEvent {
                            kind: BranchEventKind::NewComponent,
                            span,
                            depth_after: self.depth,
                        })
                    } else {
                        None
                    }
                }
                _ => None,
            };
            self.idx += 1;
            if let Some(ev) = ev {
                return Some(ev);
            }
        }
        None
    }
}
