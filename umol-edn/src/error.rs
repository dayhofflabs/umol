//! Error types for EDN parsing.

use std::error::Error as StdError;
use std::fmt;

#[cfg(feature = "serde")]
use serde::de::Error as SerdeDeError;
#[cfg(feature = "serde")]
use serde::ser::Error as SerdeSerError;
use winnow::error::{ErrMode, ParserError};
use winnow::stream::Location;
use winnow::LocatingSlice;

/// EDN parsing/formatting error.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum EdnError {
    UnexpectedEof {
        offset: usize,
    },
    UnexpectedToken {
        offset: usize,
        found: char,
    },
    InvalidNumber {
        offset: usize,
    },
    InvalidEscape {
        offset: usize,
    },
    InvalidCharLiteral {
        offset: usize,
    },
    InvalidSymbol {
        offset: usize,
    },
    DuplicateKey {
        offset: usize,
    },
    InvalidTag {
        offset: usize,
        tag: String,
    },
    InvalidInst {
        reason: String,
    },
    InvalidUuid {
        reason: String,
    },
    TrailingContent {
        offset: usize,
    },
    UnsupportedFeature {
        offset: usize,
        feature: &'static str,
    },
    /// Required map field absent. Raised by `FromEdn` impls and the
    /// `EdnMapHelper` builder.
    MissingField {
        key: String,
        path: Vec<String>,
    },
    /// Strict-mode finalize found a key the consumer did not consume.
    UnknownField {
        key: String,
        path: Vec<String>,
    },
    /// Wrong `Edn` variant for the requested type. `expected` and `got` are
    /// short discriminator strings (e.g. `"map"`, `"keyword"`, `"int"`).
    TypeMismatch {
        expected: &'static str,
        got: &'static str,
        path: Vec<String>,
    },
    /// Numeric value did not fit the target type.
    OutOfRange {
        value: String,
        target: &'static str,
        path: Vec<String>,
    },
    Custom(String),
}

/// Format a `path` slice as a dotted EDN access path for error messages.
/// Empty path renders as `<root>`.
fn fmt_path(path: &[String]) -> String {
    if path.is_empty() {
        "<root>".to_string()
    } else {
        path.join("")
    }
}

impl fmt::Display for EdnError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            EdnError::UnexpectedEof { offset } => {
                write!(f, "unexpected end of input at byte {offset}")
            }
            EdnError::UnexpectedToken { offset, found } => {
                write!(f, "unexpected token '{found}' at byte {offset}")
            }
            EdnError::InvalidNumber { offset } => {
                write!(f, "invalid number at byte {offset}")
            }
            EdnError::InvalidEscape { offset } => {
                write!(f, "invalid escape sequence at byte {offset}")
            }
            EdnError::InvalidCharLiteral { offset } => {
                write!(f, "invalid character literal at byte {offset}")
            }
            EdnError::InvalidSymbol { offset } => {
                write!(f, "invalid symbol at byte {offset}")
            }
            EdnError::DuplicateKey { offset } => {
                write!(f, "duplicate map key at byte {offset}")
            }
            EdnError::InvalidTag { offset, tag } => {
                write!(f, "invalid tag '{tag}' at byte {offset}")
            }
            EdnError::InvalidInst { reason } => {
                write!(f, "invalid #inst: {reason}")
            }
            EdnError::InvalidUuid { reason } => {
                write!(f, "invalid #uuid: {reason}")
            }
            EdnError::TrailingContent { offset } => {
                write!(f, "trailing content at byte {offset}")
            }
            EdnError::UnsupportedFeature { offset, feature } => {
                write!(f, "unsupported feature '{feature}' at byte {offset}")
            }
            EdnError::MissingField { key, path } => {
                write!(f, "missing required field '{key}' at {}", fmt_path(path))
            }
            EdnError::UnknownField { key, path } => {
                write!(f, "unknown field '{key}' at {}", fmt_path(path))
            }
            EdnError::TypeMismatch {
                expected,
                got,
                path,
            } => {
                write!(
                    f,
                    "type mismatch at {}: expected {expected}, got {got}",
                    fmt_path(path)
                )
            }
            EdnError::OutOfRange {
                value,
                target,
                path,
            } => {
                write!(
                    f,
                    "value '{value}' out of range for {target} at {}",
                    fmt_path(path)
                )
            }
            EdnError::Custom(msg) => write!(f, "{msg}"),
        }
    }
}

impl StdError for EdnError {}

impl<'a> ParserError<LocatingSlice<&'a str>> for EdnError {
    type Inner = Self;

    fn from_input(input: &LocatingSlice<&'a str>) -> Self {
        let offset = input.current_token_start();
        match input.as_ref().chars().next() {
            Some(c) => EdnError::UnexpectedToken { offset, found: c },
            None => EdnError::UnexpectedEof { offset },
        }
    }

    fn into_inner(self) -> Result<Self::Inner, Self> {
        Ok(self)
    }
}

/// Extract `EdnError` from `ErrMode<EdnError>`.
pub(crate) fn unwrap_err(e: ErrMode<EdnError>) -> EdnError {
    match e {
        ErrMode::Backtrack(e) | ErrMode::Cut(e) => e,
        // Incomplete means the parser needs more input — treat as EOF.
        // No precise offset is available from winnow's Incomplete variant;
        // report 0 (callers that know the input length can improve this).
        ErrMode::Incomplete(_) => EdnError::UnexpectedEof { offset: 0 },
    }
}

#[cfg(feature = "serde")]
impl SerdeDeError for EdnError {
    fn custom<T: fmt::Display>(msg: T) -> Self {
        EdnError::Custom(msg.to_string())
    }
}

#[cfg(feature = "serde")]
impl SerdeSerError for EdnError {
    fn custom<T: fmt::Display>(msg: T) -> Self {
        EdnError::Custom(msg.to_string())
    }
}
