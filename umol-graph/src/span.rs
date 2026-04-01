//! Span type for TableIR.

// TODO: Move to table_ir::span.rs
use std::fmt;

/// Span type for TableIR.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Span {
    Bytes { start: u32, end: u32 },
    Line { line: u32 },
    None,
}

impl Span {
    pub fn bytes(start: u32, end: u32) -> Self {
        Span::Bytes { start, end }
    }

    pub fn line(line: u32) -> Self {
        Span::Line { line }
    }

    pub fn from_bytes_opt(start: Option<u32>, end: Option<u32>) -> Option<Self> {
        start.zip(end).map(|(s, e)| Span::bytes(s, e))
    }

    pub fn bytes_range(&self) -> Option<(u32, u32)> {
        match *self {
            Span::Bytes { start, end } => Some((start, end)),
            Span::Line { .. } => None,
            Span::None => None,
        }
    }
}

impl fmt::Display for Span {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Span::Bytes { start, end } => write!(f, "@{}..{}", start, end),
            Span::Line { line } => write!(f, "@line {}", line),
            Span::None => write!(f, "none"),
        }
    }
}
