//! Domain errors for DSL parsing

use nom::error::{ErrorKind as NomErrorKind, ParseError as NomParseError};
use thiserror::Error;
use umol_data::SpinStateError;
use umol_edn::{EdnError, ParseError as EdnParseError};

#[derive(Clone, Debug, PartialEq, Error)]
pub enum AtomDslError {
    #[error("Invalid atom element: {0}")]
    InvalidElement(String),
    #[error("Invalid isotope: {0}")]
    InvalidIsotope(String),
    #[error("Invalid charge: {0}")]
    InvalidCharge(String),
    #[error("Invalid implicit hydrogens: {0}")]
    InvalidImplicitHydrogens(String),
    #[error("Invalid aromatic valence: {0}")]
    InvalidAromaticValence(String),
    #[error("Invalid value: {0}")]
    InvalidValue(String),
    #[error("Unknown atom predicate: {0}")]
    UnknownAtomPredicate(String),
    #[error("Duplicate {0} atom predicate")]
    DuplicateAtomPredicate(String),
    #[error("Trailing input: {0:?}")]
    TrailingInput(String),
    #[error("Incomplete input")]
    Incomplete,
    #[error("Nom error: {0:?}")]
    NomError(NomErrorKind),
}

impl<I> NomParseError<I> for AtomDslError {
    fn from_error_kind(_input: I, kind: NomErrorKind) -> Self {
        AtomDslError::NomError(kind)
    }

    fn append(_input: I, _kind: NomErrorKind, other: Self) -> Self {
        other
    }
}

#[derive(Clone, Debug, PartialEq, Error)]
pub enum BondDslError {
    #[error("Invalid bond order: {0}")]
    InvalidBondOrder(String),
    #[error("Invalid charge: {0}")]
    InvalidCharge(String),
    #[error("Unknown bond predicate: {0}")]
    UnknownBondPredicate(String),
    #[error("Duplicate {0} bond predicate")]
    DuplicateBondPredicate(String),
    #[error("Invalid value: {0}")]
    InvalidValue(String),
    #[error("Trailing input: {0:?}")]
    TrailingInput(String),
    #[error("Incomplete input")]
    Incomplete,
    #[error("Nom error: {0:?}")]
    NomError(NomErrorKind),
}

impl<I> NomParseError<I> for BondDslError {
    fn from_error_kind(_input: I, kind: NomErrorKind) -> Self {
        BondDslError::NomError(kind)
    }

    fn append(_input: I, _kind: NomErrorKind, other: Self) -> Self {
        other
    }
}

#[derive(Clone, Debug, PartialEq, Error)]
pub enum ParseError {
    #[error("Invalid number")]
    InvalidNumber,
    #[error("Trailing input: {0:?}")]
    TrailingInput(String),
    #[error("Invalid value: {0}")]
    InvalidValue(String),
    #[error("Nom error: {0:?}")]
    NomError(NomErrorKind),
    #[error("EDN parse error: {0}")]
    EdnParse(String),
    #[error("Missing required key: {0}")]
    MissingKey(String),
    #[error("expected {expected} for :{field}")]
    WrongFieldType { field: String, expected: String },
    #[error("invalid atom DSL: {0}")]
    InvalidAtomSpec(String),
    #[error("invalid bond DSL: {0}")]
    InvalidBondSpec(String),
    #[error("invalid bond entry: expected map-based {{[:id keyword] :a :b :bond}} or vector-based [a b bond-spec]")]
    InvalidBond,
    #[error("Duplicate structural id: {0}")]
    DuplicateId(String),
    #[error("Unknown atom endpoint: {0}")]
    InvalidAtomIndex(String),
    #[error("Unknown alias: {0}")]
    UnknownAlias(String),
    #[error("Invalid spin state: {0}")]
    InvalidSpinState(#[from] SpinStateError),
    #[error("Incomplete input")]
    Incomplete,
}

impl From<AtomDslError> for ParseError {
    fn from(e: AtomDslError) -> Self {
        ParseError::InvalidAtomSpec(e.to_string())
    }
}

impl From<BondDslError> for ParseError {
    fn from(e: BondDslError) -> Self {
        ParseError::InvalidBondSpec(e.to_string())
    }
}

impl From<EdnError> for ParseError {
    fn from(e: EdnError) -> Self {
        match &e {
            EdnError::Parse(EdnParseError::UnexpectedEof { .. }) => ParseError::Incomplete,
            _ => ParseError::EdnParse(e.to_string()),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Error)]
pub enum EvaluationError {
    #[error("Unbound variable: {0}")]
    UnboundVariable(String),
    #[error("Division by zero")]
    DivisionByZero,
    #[error("Type mismatch")]
    TypeMismatch,
}

#[derive(Clone, Debug, PartialEq, Error)]
pub enum LoweringError {
    #[error("non-ground value for field '{field}'")]
    NonGround { field: &'static str },
    #[error("value {value} out of range for field '{field}'")]
    OutOfRange { field: &'static str, value: i64 },
    #[error("field '{field}' is required but not present")]
    MissingField { field: &'static str },
    #[error("invalid spin multiplicity: {0}")]
    InvalidMultiplicity(u8),
    #[error("incompatible spin state: {0}")]
    SpinState(#[from] SpinStateError),
    #[error("invalid atom spec: {0}")]
    Atom(String),
    #[error("unknown atom label: {0}")]
    UnknownLabel(String),
    #[error("invalid molecule spec: {0}")]
    Molecule(String),
}
