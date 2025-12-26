use crate::diagnostics::{Diagnostic, DiagnosticKind, Severity};
use crate::span::Span;
use nom::error::{ErrorKind, ParseError as NomParseErrorTrait};
use nom::Err;
use thiserror::Error;

#[derive(Debug, Error, Clone, PartialEq)]
pub enum ParseError {
    #[error("Invalid header block at line {line}")]
    InvalidHeader { line: u32 },
    #[error("Invalid counts line at line {line}")]
    InvalidCountsLine { line: u32 },
    #[error("Invalid atom line at line {line}, col {col}")]
    InvalidAtomLine { line: u32, col: u32 },
    #[error("Invalid bond line at line {line}, col {col}")]
    InvalidBondLine { line: u32, col: u32 },
    #[error("Invalid property line at line {line}, col {col}")]
    InvalidPropertyLine { line: u32, col: u32 },
    #[error("Invalid Sgroup line at line {line}, col {col}")]
    InvalidSgroupLine { line: u32, col: u32 },
    #[error("Invalid SDF data header at line {line}")]
    InvalidSdfDataHeader { line: u32 },
    #[error("Invalid SDF data value at line {line}")]
    InvalidSdfDataValue { line: u32 },
    #[error("Missing record delimiter at line {line}")]
    MissingDelimiter { line: u32 },
    #[error("Parse error at line {line}, col {col}: {message}")]
    Generic {
        line: u32,
        col: u32,
        message: String,
    },
}

#[derive(Debug, Error, Clone, PartialEq)]
pub enum SemanticError {
    #[error("Invalid charge code: {0}")]
    InvalidChargeCode(u8),
    #[error("Invalid valence code: {0}")]
    InvalidValenceCode(u8),
    #[error("Invalid stereo parity code: {0}")]
    InvalidStereoParity(u8),
    #[error("Property mismatch: {0}")]
    PropertyMismatch(String),
    #[error("Inconsistent Sgroups: {0}")]
    InconsistentSgroups(String),
    #[error("Index out of bounds: {0}")]
    IndexOutOfBounds(usize),
    #[error("Incomplete structure: {0}")]
    IncompleteStructure(String),
    #[error("Duplicate property: {0}")]
    DuplicateProperty(String),
    #[error("Generic semantic error: {0}")]
    Generic(String),
}

impl ParseError {
    pub fn generic(line: u32, col: u32, message: impl Into<String>) -> Self {
        ParseError::Generic {
            line,
            col,
            message: message.into(),
        }
    }

    pub fn from_nom(e: Err<nom::error::Error<&[u8]>>, line: u32, line_input: &[u8]) -> Self {
        match e {
            Err::Error(err) | Err::Failure(err) => {
                let col = (err.input.as_ptr() as usize) - (line_input.as_ptr() as usize);
                ParseError::generic(line, col as u32, format!("nom error: {:?}", err.code))
            }
            Err::Incomplete(_) => ParseError::generic(line, 0, "Incomplete"),
        }
    }

    pub fn with_line(mut self, line: u32) -> Self {
        match &mut self {
            ParseError::InvalidHeader { line: l } => *l = line,
            ParseError::InvalidCountsLine { line: l } => *l = line,
            ParseError::InvalidAtomLine { line: l, .. } => *l = line,
            ParseError::InvalidBondLine { line: l, .. } => *l = line,
            ParseError::InvalidPropertyLine { line: l, .. } => *l = line,
            ParseError::InvalidSgroupLine { line: l, .. } => *l = line,
            ParseError::InvalidSdfDataHeader { line: l } => *l = line,
            ParseError::InvalidSdfDataValue { line: l } => *l = line,
            ParseError::MissingDelimiter { line: l } => *l = line,
            ParseError::Generic { line: l, .. } => *l = line,
        }
        self
    }
}

impl NomParseErrorTrait<&[u8]> for ParseError {
    fn from_error_kind(_input: &[u8], kind: ErrorKind) -> Self {
        ParseError::Generic {
            line: 0,
            col: 0,
            message: format!("nom error: {:?}", kind),
        }
    }

    fn append(_input: &[u8], _kind: ErrorKind, other: Self) -> Self {
        other
    }
}

impl From<SemanticError> for Diagnostic {
    fn from(error: SemanticError) -> Self {
        use DiagnosticKind::*;
        let (kind, details) = match error {
            SemanticError::InvalidChargeCode(_) => (CtfileInvalidPropertyLine, Some(error.to_string())),
            SemanticError::InvalidValenceCode(_) => (CtfileInvalidPropertyLine, Some(error.to_string())),
            SemanticError::InvalidStereoParity(_) => (CtfileInvalidAtomLine, Some(error.to_string())),
            SemanticError::PropertyMismatch(ref s) => (CtfileInvalidPropertyLine, Some(s.clone())),
            SemanticError::InconsistentSgroups(ref s) => (CtfileInvalidSgroupLine, Some(s.clone())),
            SemanticError::IndexOutOfBounds(_) => (Unknown, Some(error.to_string())),
            SemanticError::IncompleteStructure(ref s) => (CtfileInvalidCountsLine, Some(s.clone())),
            SemanticError::DuplicateProperty(ref s) => (CtfileInvalidPropertyLine, Some(s.clone())),
            SemanticError::Generic(ref s) => (Unknown, Some(s.clone())),
        };
        Diagnostic {
            kind,
            category: kind.category(),
            severity: Severity::Error,
            span: None,
            details,
        }
    }
}

impl From<SemanticError> for umol::error::ParseError {
    fn from(error: SemanticError) -> Self {
        umol::error::ParseError::Format(Box::new(error))
    }
}

impl From<SemanticError> for umol::Error {
    fn from(error: SemanticError) -> Self {
        umol::Error::Parse(error.into())
    }
}

impl From<ParseError> for Diagnostic {
    fn from(error: ParseError) -> Self {
        let (kind, span, details) = match error {
            ParseError::InvalidHeader { line } => (
                DiagnosticKind::CtfileInvalidHeader,
                Span::line(line, 0, 0),
                None,
            ),
            ParseError::InvalidCountsLine { line } => (
                DiagnosticKind::CtfileInvalidCountsLine,
                Span::line(line, 0, 0),
                None,
            ),
            ParseError::InvalidAtomLine { line, col } => (
                DiagnosticKind::CtfileInvalidAtomLine,
                Span::line(line, col, 1),
                None,
            ),
            ParseError::InvalidBondLine { line, col } => (
                DiagnosticKind::CtfileInvalidBondLine,
                Span::line(line, col, 1),
                None,
            ),
            ParseError::InvalidPropertyLine { line, col } => (
                DiagnosticKind::CtfileInvalidPropertyLine,
                Span::line(line, col, 1),
                None,
            ),
            ParseError::InvalidSgroupLine { line, col } => (
                DiagnosticKind::CtfileInvalidSgroupLine,
                Span::line(line, col, 1),
                None,
            ),
            ParseError::InvalidSdfDataHeader { line } => (
                DiagnosticKind::CtfileInvalidSdfDataHeader,
                Span::line(line, 0, 0),
                None,
            ),
            ParseError::InvalidSdfDataValue { line } => (
                DiagnosticKind::CtfileInvalidSdfDataValue,
                Span::line(line, 0, 0),
                None,
            ),
            ParseError::MissingDelimiter { line } => (
                DiagnosticKind::CtfileMissingDelimiter,
                Span::line(line, 0, 0),
                None,
            ),
            ParseError::Generic {
                line,
                col,
                message,
            } => (
                DiagnosticKind::Unknown,
                Span::line(line, col, 1),
                Some(message),
            ),
        };
        Diagnostic {
            kind,
            category: kind.category(),
            severity: Severity::Error,
            span: Some(span),
            details,
        }
    }
}

impl From<ParseError> for umol::error::ParseError {
    fn from(error: ParseError) -> Self {
        umol::error::ParseError::Format(Box::new(error))
    }
}

impl From<ParseError> for umol::Error {
    fn from(error: ParseError) -> Self {
        umol::Error::Parse(error.into())
    }
}

