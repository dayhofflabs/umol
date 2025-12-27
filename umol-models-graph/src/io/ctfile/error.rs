use nom::Err;
use thiserror::Error;

use crate::diagnostics::{Diagnostic, DiagnosticKind, Severity};
use crate::span::Span;

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
    #[error("Unexpected end of file in {block} block at line {line}")]
    UnexpectedEof { line: u32, block: &'static str },
    #[error("Incomplete input at line {line}")]
    Incomplete { line: u32 },
    #[error("Parse error at line {line}, col {col}: {message}")]
    Generic {
        line: u32,
        col: u32,
        message: String,
    },
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
            Err::Incomplete(_) => ParseError::Incomplete { line },
        }
    }

    /// Create InvalidCountsLine from a nom error
    pub fn counts_from_nom(e: Err<nom::error::Error<&[u8]>>, line: u32) -> Self {
        match e {
            Err::Incomplete(_) => ParseError::Incomplete { line },
            _ => ParseError::InvalidCountsLine { line },
        }
    }

    /// Create InvalidAtomLine from a nom error
    pub fn atom_from_nom(e: Err<nom::error::Error<&[u8]>>, line: u32, line_input: &[u8]) -> Self {
        match e {
            Err::Error(err) | Err::Failure(err) => {
                let col = (err.input.as_ptr() as usize) - (line_input.as_ptr() as usize);
                ParseError::InvalidAtomLine { line, col: col as u32 }
            }
            Err::Incomplete(_) => ParseError::Incomplete { line },
        }
    }

    /// Create InvalidBondLine from a nom error
    pub fn bond_from_nom(e: Err<nom::error::Error<&[u8]>>, line: u32, line_input: &[u8]) -> Self {
        match e {
            Err::Error(err) | Err::Failure(err) => {
                let col = (err.input.as_ptr() as usize) - (line_input.as_ptr() as usize);
                ParseError::InvalidBondLine { line, col: col as u32 }
            }
            Err::Incomplete(_) => ParseError::Incomplete { line },
        }
    }

    /// Create InvalidHeader from a nom error
    pub fn header_from_nom(e: Err<nom::error::Error<&[u8]>>, line: u32) -> Self {
        match e {
            Err::Incomplete(_) => ParseError::Incomplete { line },
            _ => ParseError::InvalidHeader { line },
        }
    }

    /// Create InvalidSdfDataHeader from a nom error
    pub fn sdf_data_from_nom(e: Err<nom::error::Error<&[u8]>>, line: u32) -> Self {
        match e {
            Err::Incomplete(_) => ParseError::Incomplete { line },
            _ => ParseError::InvalidSdfDataHeader { line },
        }
    }

    /// Create MissingDelimiter from a nom error
    pub fn delimiter_from_nom(e: Err<nom::error::Error<&[u8]>>, line: u32) -> Self {
        match e {
            Err::Incomplete(_) => ParseError::Incomplete { line },
            _ => ParseError::MissingDelimiter { line },
        }
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
            ParseError::UnexpectedEof { line, block } => (
                DiagnosticKind::CtfileUnexpectedEof,
                Span::line(line, 0, 0),
                Some(format!("in {} block", block)),
            ),
            ParseError::Incomplete { line } => (
                DiagnosticKind::CtfileIncomplete,
                Span::line(line, 0, 0),
                None,
            ),
            ParseError::Generic { line, col, message } => (
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
    #[error("Invalid {field} code: {value}")]
    InvalidCode { field: &'static str, value: i32 },
    #[error("Invalid isotope mass {mass} for element {element}")]
    InvalidIsotopeMass {
        mass: u32,
        element: umol_data::Element,
    },
    #[error("Undefined S-group {index}: {property}")]
    UndefinedSGroup {
        index: usize,
        property: &'static str,
    },
    #[error("S-group {0} has no type")]
    SGroupMissingType(usize),
    #[error("S-group type {sgroup_type:?}: {message}")]
    SGroupTypeConstraint {
        sgroup_type: crate::table_ir::SGroupType,
        message: &'static str,
    },
    #[error("Missing S-group data context")]
    MissingSGroupDataContext,
    #[error("S-group index mismatch: expected {expected}, got {actual}")]
    SGroupIndexMismatch { expected: usize, actual: usize },
}

impl From<SemanticError> for Diagnostic {
    fn from(error: SemanticError) -> Self {
        use DiagnosticKind::*;
        let (kind, details) = match &error {
            SemanticError::InvalidChargeCode(_) => {
                (CtfileInvalidPropertyLine, Some(error.to_string()))
            }
            SemanticError::InvalidValenceCode(_) => {
                (CtfileInvalidPropertyLine, Some(error.to_string()))
            }
            SemanticError::InvalidStereoParity(_) => {
                (CtfileInvalidAtomLine, Some(error.to_string()))
            }
            SemanticError::PropertyMismatch(s) => (CtfileInvalidPropertyLine, Some(s.clone())),
            SemanticError::InconsistentSgroups(s) => (CtfileInvalidSgroupLine, Some(s.clone())),
            SemanticError::IndexOutOfBounds(_) => {
                (CtfileInvalidPropertyLine, Some(error.to_string()))
            }
            SemanticError::IncompleteStructure(s) => (CtfileInvalidCountsLine, Some(s.clone())),
            SemanticError::DuplicateProperty(s) => (CtfileInvalidPropertyLine, Some(s.clone())),
            SemanticError::InvalidCode { .. } => {
                (CtfileInvalidPropertyLine, Some(error.to_string()))
            }
            SemanticError::InvalidIsotopeMass { .. } => {
                (CtfileInvalidPropertyLine, Some(error.to_string()))
            }
            SemanticError::UndefinedSGroup { .. } => {
                (CtfileInvalidSgroupLine, Some(error.to_string()))
            }
            SemanticError::SGroupMissingType(_) => {
                (CtfileInvalidSgroupLine, Some(error.to_string()))
            }
            SemanticError::SGroupTypeConstraint { .. } => {
                (CtfileInvalidSgroupLine, Some(error.to_string()))
            }
            SemanticError::MissingSGroupDataContext => {
                (CtfileInvalidSgroupLine, Some(error.to_string()))
            }
            SemanticError::SGroupIndexMismatch { .. } => {
                (CtfileInvalidSgroupLine, Some(error.to_string()))
            }
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
