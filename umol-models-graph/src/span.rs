use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Span {
    Bytes { start: u32, end: u32 },
}

impl Span {
    pub fn bytes(start: u32, end: u32) -> Self {
        Span::Bytes { start, end }
    }

    pub fn from_bytes_opt(start: Option<u32>, end: Option<u32>) -> Option<Self> {
        start.zip(end).map(|(s, e)| Span::bytes(s, e))
    }

    pub fn bytes_range(&self) -> (u32, u32) {
        match *self {
            Span::Bytes { start, end } => (start, end),
        }
    }

    pub fn bytes_opt(&self) -> Option<(u32, u32)> {
        Some(self.bytes_range())
    }

    pub fn with_start(self, start: u32) -> Self {
        match self {
            Span::Bytes { end, .. } => Span::bytes(start, end),
        }
    }

    pub fn with_end(self, end: u32) -> Self {
        match self {
            Span::Bytes { start, .. } => Span::bytes(start, end),
        }
    }
}
