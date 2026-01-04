//! Error types for CTFile parsing

use nom::error::{Error as NomError, ErrorKind as NomErrorKind, ParseError as NomParseError};
use nom::Err;
use thiserror::Error;
use umol_data::Element;

use crate::diagnostics::{Diagnostic, DiagnosticKind, Severity};
use crate::span::Span;
use crate::table_ir::SGroupType;

#[derive(Debug, Clone, PartialEq, Error)]
pub enum ParseError {
    #[error("Invalid header block at line {line}")]
    InvalidHeader { line: u32 },
    #[error("Invalid counts line at line {line}")]
    InvalidCountsLine { line: u32 },
    #[error("Invalid atom line at line {line}, col {col}")]
    InvalidAtomLine { line: u32, col: u32 },
    #[error("Invalid bond line at line {line}, col {col}")]
    InvalidBondLine { line: u32, col: u32 },
    #[error("Invalid legacy atom list line at line {line}, col {col}")]
    InvalidLegacyAtomListLine { line: u32, col: u32 },
    #[error("Unsupported legacy atom list at line {line}")]
    UnsupportedLegacyAtomList { line: u32 },
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
    #[error("Missing M  END tag at line {line}")]
    MissingMEndTag { line: u32 },
    #[error("Unexpected end of file in {block} block at line {line}")]
    UnexpectedEof { line: u32, block: &'static str },
    #[error("Incomplete input at line {line}")]
    Incomplete { line: u32 },
    #[error("Nom parser error: {0:?}")]
    NomError(NomErrorKind),
    #[error("Parse error at line {line}, col {col}: {message}")]
    Generic {
        line: u32,
        col: u32,
        message: String,
    },
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
    InvalidIsotopeMass { mass: u32, element: Element },
    #[error("Undefined S-group {index}: {property}")]
    UndefinedSGroup {
        index: usize,
        property: &'static str,
    },
    #[error("S-group {0} has no type")]
    SGroupMissingType(usize),
    #[error("S-group type {sgroup_type:?}: {message}")]
    SGroupTypeConstraint {
        sgroup_type: SGroupType,
        message: &'static str,
    },
    #[error("Missing S-group data context")]
    MissingSGroupDataContext,
    #[error("S-group index mismatch: expected {expected}, got {actual}")]
    SGroupIndexMismatch { expected: usize, actual: usize },
}

/// TODO: Use LineSpan wrapper to track column offset
impl<I> NomParseError<I> for ParseError {
    fn from_error_kind(_input: I, kind: NomErrorKind) -> Self {
        ParseError::NomError(kind)
    }

    fn append(_input: I, _kind: NomErrorKind, other: Self) -> Self {
        other
    }
}

impl From<Err<ParseError>> for ParseError {
    fn from(e: Err<ParseError>) -> Self {
        match e {
            Err::Error(inner) | Err::Failure(inner) => inner,
            Err::Incomplete(_) => ParseError::Incomplete { line: 0 },
        }
    }
}

impl ParseError {
    pub fn generic(line: u32, col: u32, message: impl Into<String>) -> Self {
        ParseError::Generic {
            line,
            col,
            message: message.into(),
        }
    }

    pub fn from_nom(e: Err<NomError<&[u8]>>, line: u32, line_input: &[u8]) -> Self {
        match e {
            Err::Error(err) | Err::Failure(err) => {
                let col = (err.input.as_ptr() as usize) - (line_input.as_ptr() as usize);
                ParseError::generic(line, col as u32, format!("nom error: {:?}", err.code))
            }
            Err::Incomplete(_) => ParseError::Incomplete { line },
        }
    }

    /// Create InvalidCountsLine from a nom error
    pub fn counts_from_nom(e: Err<NomError<&[u8]>>, line: u32) -> Self {
        match e {
            Err::Incomplete(_) => ParseError::Incomplete { line },
            _ => ParseError::InvalidCountsLine { line },
        }
    }

    /// Create InvalidAtomLine from a nom error
    pub fn atom_from_nom(e: Err<NomError<&[u8]>>, line: u32, line_input: &[u8]) -> Self {
        match e {
            Err::Error(err) | Err::Failure(err) => {
                let col = (err.input.as_ptr() as usize) - (line_input.as_ptr() as usize);
                ParseError::InvalidAtomLine {
                    line,
                    col: col as u32,
                }
            }
            Err::Incomplete(_) => ParseError::Incomplete { line },
        }
    }

    /// Create InvalidBondLine from a nom error
    pub fn bond_from_nom(e: Err<NomError<&[u8]>>, line: u32, line_input: &[u8]) -> Self {
        match e {
            Err::Error(err) | Err::Failure(err) => {
                let col = (err.input.as_ptr() as usize) - (line_input.as_ptr() as usize);
                ParseError::InvalidBondLine {
                    line,
                    col: col as u32,
                }
            }
            Err::Incomplete(_) => ParseError::Incomplete { line },
        }
    }

    /// Create InvalidLegacyAtomListLine from a nom error
    pub fn legacy_atom_list_from_nom(
        e: Err<NomError<&[u8]>>,
        line: u32,
        line_input: &[u8],
    ) -> Self {
        match e {
            Err::Error(err) | Err::Failure(err) => {
                let col = (err.input.as_ptr() as usize) - (line_input.as_ptr() as usize);
                ParseError::InvalidLegacyAtomListLine {
                    line,
                    col: col as u32,
                }
            }
            Err::Incomplete(_) => ParseError::Incomplete { line },
        }
    }

    /// Create InvalidHeader from a nom error
    pub fn header_from_nom(e: Err<NomError<&[u8]>>, line: u32) -> Self {
        match e {
            Err::Incomplete(_) => ParseError::Incomplete { line },
            _ => ParseError::InvalidHeader { line },
        }
    }

    /// Create InvalidSdfDataHeader from a nom error
    pub fn sdf_data_from_nom(e: Err<NomError<&[u8]>>, line: u32) -> Self {
        match e {
            Err::Incomplete(_) => ParseError::Incomplete { line },
            _ => ParseError::InvalidSdfDataHeader { line },
        }
    }

    /// Create MissingDelimiter from a nom error
    pub fn delimiter_from_nom(e: Err<NomError<&[u8]>>, line: u32) -> Self {
        match e {
            Err::Incomplete(_) => ParseError::Incomplete { line },
            _ => ParseError::MissingDelimiter { line },
        }
    }

    /// Create MissingMEndTag from a nom error
    pub fn m_end_from_nom(e: Err<NomError<&[u8]>>, line: u32) -> Self {
        match e {
            Err::Incomplete(_) => ParseError::Incomplete { line },
            _ => ParseError::MissingMEndTag { line },
        }
    }
}

impl From<ParseError> for Diagnostic {
    fn from(error: ParseError) -> Self {
        let (kind, span, details) = match error {
            ParseError::InvalidHeader { line } => (
                DiagnosticKind::CtfileInvalidHeader,
                Span::line(line),
                None,
            ),
            ParseError::InvalidCountsLine { line } => (
                DiagnosticKind::CtfileInvalidCountsLine,
                Span::line(line),
                None,
            ),
            ParseError::InvalidAtomLine { line, .. } => (
                DiagnosticKind::CtfileInvalidAtomLine,
                Span::line(line),
                None,
            ),
            ParseError::InvalidBondLine { line, .. } => (
                DiagnosticKind::CtfileInvalidBondLine,
                Span::line(line),
                None,
            ),
            ParseError::InvalidLegacyAtomListLine { line, .. } => (
                DiagnosticKind::CtfileInvalidLegacyAtomListLine,
                Span::line(line),
                None,
            ),
            ParseError::UnsupportedLegacyAtomList { line } => (
                DiagnosticKind::CtfileInvalidLegacyAtomListLine,
                Span::line(line),
                None,
            ),
            ParseError::InvalidPropertyLine { line, .. } => (
                DiagnosticKind::CtfileInvalidPropertyLine,
                Span::line(line),
                None,
            ),
            ParseError::InvalidSgroupLine { line, .. } => (
                DiagnosticKind::CtfileInvalidSgroupLine,
                Span::line(line),
                None,
            ),
            ParseError::InvalidSdfDataHeader { line } => (
                DiagnosticKind::CtfileInvalidSdfDataHeader,
                Span::line(line),
                None,
            ),
            ParseError::InvalidSdfDataValue { line } => (
                DiagnosticKind::CtfileInvalidSdfDataValue,
                Span::line(line),
                None,
            ),
            ParseError::MissingDelimiter { line } => (
                DiagnosticKind::CtfileMissingDelimiter,
                Span::line(line),
                None,
            ),
            ParseError::MissingMEndTag { line } => (
                DiagnosticKind::CtfileMissingMEndTag,
                Span::line(line),
                None,
            ),
            ParseError::UnexpectedEof { line, block } => (
                DiagnosticKind::CtfileUnexpectedEof,
                Span::line(line),
                Some(format!("in {} block", block)),
            ),
            ParseError::Incomplete { line } => (
                DiagnosticKind::CtfileIncomplete,
                Span::line(line),
                None,
            ),
            ParseError::NomError(kind) => (
                DiagnosticKind::CtfileParserError,
                Span::None,
                Some(format!("nom error: {:?}", kind)),
            ),
            ParseError::Generic { line, message, .. } => (
                DiagnosticKind::CtfileInvalidPropertyLine,
                Span::line(line),
                Some(message),
            ),
            ParseError::InvalidChargeCode(_) => (
                DiagnosticKind::CtfileInvalidPropertyLine,
                Span::None,
                Some(error.to_string()),
            ),
            ParseError::InvalidValenceCode(_) => (
                DiagnosticKind::CtfileInvalidPropertyLine,
                Span::None,
                Some(error.to_string()),
            ),
            ParseError::InvalidStereoParity(_) => (
                DiagnosticKind::CtfileInvalidAtomLine,
                Span::None,
                Some(error.to_string()),
            ),
            ParseError::PropertyMismatch(s) => (
                DiagnosticKind::CtfileInvalidPropertyLine,
                Span::None,
                Some(s.clone()),
            ),
            ParseError::InconsistentSgroups(s) => (
                DiagnosticKind::CtfileInvalidSgroupLine,
                Span::None,
                Some(s.clone()),
            ),
            ParseError::IndexOutOfBounds(_) => (
                DiagnosticKind::CtfileInvalidPropertyLine,
                Span::None,
                Some(error.to_string()),
            ),
            ParseError::IncompleteStructure(s) => (
                DiagnosticKind::CtfileInvalidCountsLine,
                Span::None,
                Some(s.clone()),
            ),
            ParseError::DuplicateProperty(s) => (
                DiagnosticKind::CtfileInvalidPropertyLine,
                Span::None,
                Some(s.clone()),
            ),
            ParseError::InvalidCode { .. } => (
                DiagnosticKind::CtfileInvalidPropertyLine,
                Span::None,
                Some(error.to_string()),
            ),
            ParseError::InvalidIsotopeMass { .. } => (
                DiagnosticKind::CtfileInvalidPropertyLine,
                Span::None,
                Some(error.to_string()),
            ),
            ParseError::UndefinedSGroup { .. } => (
                DiagnosticKind::CtfileInvalidSgroupLine,
                Span::None,
                Some(error.to_string()),
            ),
            ParseError::SGroupMissingType(_) => (
                DiagnosticKind::CtfileInvalidSgroupLine,
                Span::None,
                Some(error.to_string()),
            ),
            ParseError::SGroupTypeConstraint { .. } => (
                DiagnosticKind::CtfileInvalidSgroupLine,
                Span::None,
                Some(error.to_string()),
            ),
            ParseError::MissingSGroupDataContext => (
                DiagnosticKind::CtfileInvalidSgroupLine,
                Span::None,
                Some(error.to_string()),
            ),
            ParseError::SGroupIndexMismatch { .. } => (
                DiagnosticKind::CtfileInvalidSgroupLine,
                Span::None,
                Some(error.to_string()),
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
