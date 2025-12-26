use serde::{Deserialize, Serialize};

use std::fmt;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Span {
    Bytes {
        start: u32,
        end: u32,
    },
    Line {
        line: u32,
        col: u32,
        len: u32,
    },
}

impl Span {
    pub fn bytes(start: u32, end: u32) -> Self {
        Span::Bytes { start, end }
    }

    pub fn line(line: u32, col: u32, len: u32) -> Self {
        Span::Line { line, col, len }
    }

    pub fn from_bytes_opt(start: Option<u32>, end: Option<u32>) -> Option<Self> {
        start.zip(end).map(|(s, e)| Span::bytes(s, e))
    }

    pub fn bytes_range(&self) -> Option<(u32, u32)> {
        match *self {
            Span::Bytes { start, end } => Some((start, end)),
            Span::Line { .. } => None,
        }
    }

    pub fn with_start(self, start: u32) -> Self {
        match self {
            Span::Bytes { end, .. } => Span::bytes(start, end),
            Span::Line { line, len, .. } => Span::line(line, start, len),
        }
    }

    pub fn with_end(self, end: u32) -> Self {
        match self {
            Span::Bytes { start, .. } => Span::bytes(start, end),
            Span::Line { line, col, .. } => {
                let len = end.saturating_sub(col);
                Span::line(line, col, len)
            }
        }
    }
}

impl fmt::Display for Span {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Span::Bytes { start, end } => write!(f, "@{}..{}", start, end),
            Span::Line { line, col, len } => write!(f, "@{}:{}+{}", line, col, len),
        }
    }
}
