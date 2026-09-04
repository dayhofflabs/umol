//! Error types for CTfile parsing.

use std::any::Any;

use thiserror::Error;
use umol_chem::element::Element;
use umol_utils::error::UmolError;

use crate::table_ir::SGroupType;

/// An error encountered while parsing a complete MOL or SDF input.
///
/// Syntax variants identify the CTfile construct that could not be parsed. Line and
/// column values are zero-based physical byte positions. Construction variants retain
/// the domain condition that prevented the parsed records from forming a TableIR value.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum ParseError {
    #[error("Invalid counts line at line {line}, col {col}")]
    InvalidCountsLine { line: u32, col: u32 },
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
    #[error("Invalid SDF data header at line {line}, col {col}")]
    InvalidSdfDataHeader { line: u32, col: u32 },
    #[error("Missing record delimiter at line {line}")]
    MissingDelimiter { line: u32 },
    #[error("Missing M  END tag at line {line}")]
    MissingMEndTag { line: u32 },
    #[error("Unexpected end of file in {block} block at line {line}")]
    UnexpectedEof { line: u32, block: &'static str },
    #[error("Duplicate property: {0}")]
    DuplicateProperty(String),
    #[error("Invalid {field} code: {value}")]
    InvalidCode { field: &'static str, value: i32 },
    #[error("Invalid isotope mass {mass} for element {element}")]
    InvalidIsotopeMass { mass: u32, element: Element },
    #[error("Atom index out of bounds: {0}")]
    AtomIndexOutOfBounds(u32),
    #[error("Bond index out of bounds: {0}")]
    BondIndexOutOfBounds(u32),
    #[error("Undefined S-group {index}: {property}")]
    UndefinedSGroup { index: u32, property: &'static str },
    #[error("S-group {0} has no type")]
    SGroupMissingType(u32),
    #[error("S-group type {sgroup_type:?}: {message}")]
    SGroupTypeConstraint {
        sgroup_type: SGroupType,
        message: &'static str,
    },
    #[error("Missing context for data SGroup {index} in {location}")]
    MissingSGroupDataContext { index: u32, location: &'static str },
    #[error("Unfinalized data SGroup {index}")]
    MissingSGroupDataEnd { index: u32 },
    #[error("S-group index mismatch: expected {expected}, got {actual}")]
    SGroupIndexMismatch { expected: u32, actual: u32 },
}

impl UmolError for ParseError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}
