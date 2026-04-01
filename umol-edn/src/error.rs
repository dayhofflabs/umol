//! Error types for EDN parsing.

use std::fmt;

/// Byte-offset span in the input.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Span {
    pub start: usize,
    pub end: usize,
}

/// EDN parsing/formatting error.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum EdnError {
    UnexpectedToken { span: Span, found: char },
    UnexpectedEof { offset: usize },
    IntegerOverflow { span: Span },
    InvalidEscape { span: Span, sequence: String },
    InvalidCharLiteral { span: Span },
    InvalidNumber { span: Span },
    DuplicateMapKey { span: Span },
    TrailingContent { offset: usize },
    InvalidUtf8 { offset: usize },
    Custom(String),
}

impl fmt::Display for EdnError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            EdnError::UnexpectedToken { span, found } => {
                write!(f, "unexpected token '{found}' at byte {}", span.start)
            }
            EdnError::UnexpectedEof { offset } => {
                write!(f, "unexpected end of input at byte {offset}")
            }
            EdnError::IntegerOverflow { span } => {
                write!(f, "integer overflow at byte {}", span.start)
            }
            EdnError::InvalidEscape { span, sequence } => {
                write!(f, "invalid escape sequence '{sequence}' at byte {}", span.start)
            }
            EdnError::InvalidCharLiteral { span } => {
                write!(f, "invalid character literal at byte {}", span.start)
            }
            EdnError::InvalidNumber { span } => {
                write!(f, "invalid number at byte {}", span.start)
            }
            EdnError::DuplicateMapKey { span } => {
                write!(f, "duplicate map key at byte {}", span.start)
            }
            EdnError::TrailingContent { offset } => {
                write!(f, "trailing content at byte {offset}")
            }
            EdnError::InvalidUtf8 { offset } => {
                write!(f, "invalid UTF-8 at byte {offset}")
            }
            EdnError::Custom(msg) => write!(f, "{msg}"),
        }
    }
}

impl std::error::Error for EdnError {}
