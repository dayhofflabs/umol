//! Domain errors for DSL parsing.

use thiserror::Error;
use umol_shared::error::SpinStateError;
use umol_shared::spin::SpinMultiplicity;
use winnow::error::{ErrMode, ParserError};
use winnow::stream::Stream;

pub(crate) type PResult<T> = Result<T, ErrMode<ParseError>>;

/// Error raised when a DSL input fails to parse (invalid syntax, unknown
/// predicate, duplicate predicate, unresolved ref, etc.).
#[rustfmt::skip]
#[derive(Clone, Debug, PartialEq, Error)]
pub enum ParseError {
    #[error("expected atom element")]
    ExpectedElement,
    #[error("expected predicate body")]
    ExpectedPredicateBody,
    #[error("unknown atom predicate: {0}")]
    UnknownAtomPredicate(String),
    #[error("duplicate atom predicate: {0}")]
    DuplicateAtomPredicate(String),
    #[error("unknown bond predicate: {0}")]
    UnknownBondPredicate(String),
    #[error("duplicate bond predicate: {0}")]
    DuplicateBondPredicate(String),
    #[error("unknown aromatic-system predicate: {0}")]
    UnknownAromaticSystemPredicate(String),
    #[error("duplicate aromatic-system predicate: {0}")]
    DuplicateAromaticSystemPredicate(String),
    #[error("unknown multicenter-bond predicate: {0}")]
    UnknownMulticenterBondPredicate(String),
    #[error("duplicate multicenter-bond predicate: {0}")]
    DuplicateMulticenterBondPredicate(String),
    #[error("unknown dative-bond predicate: {0}")]
    UnknownDativeBondPredicate(String),
    #[error("duplicate dative-bond predicate: {0}")]
    DuplicateDativeBondPredicate(String),
    #[error("expected noncovalent-bond kind")]
    ExpectedNoncovalentBondKind,
    #[error("trailing input: {0:?}")]
    TrailingInput(String),
    #[error("unpaired electrons {unpaired} out of range")]
    UnpairedElectronsOutOfRange { unpaired: u8 },
    #[error("multiplicity {multiplicity} out of range")]
    MultiplicityOutOfRange { multiplicity: u8 },
    #[error("{unpaired} unpaired electrons, {multiplicity} multiplicity incompatible")]
    IncompatibleSpin { unpaired: u8, multiplicity: SpinMultiplicity },
    #[error("raising error: {0}")]
    RaisingError(String),
    #[error("lowering error: {0}")]
    LoweringError(String),
    #[error("syntax error")]
    Syntax,
    #[error("EDN parse: {0}")]
    EdnParse(String),
    #[error("missing key: {0}")]
    MissingKey(String),
    #[error("duplicate id: {0}")]
    DuplicateId(String),
    #[error("invalid value: {0}")]
    InvalidValue(String),
    #[error("invalid {kind} ref: {value}")]
    InvalidRef { kind: &'static str, value: String },
    #[error("{field}: expected {expected}")]
    WrongFieldType { field: String, expected: String },
}

impl ParseError {
    pub(crate) fn from_spin_state_error(err: SpinStateError) -> Self {
        match err {
            SpinStateError::UnpairedElectronsOutOfRange { unpaired } => {
                ParseError::UnpairedElectronsOutOfRange { unpaired }
            }
            SpinStateError::MultiplicityOutOfRange { multiplicity } => {
                ParseError::MultiplicityOutOfRange { multiplicity }
            }
            SpinStateError::Incompatible {
                unpaired,
                multiplicity,
            } => ParseError::IncompatibleSpin {
                unpaired,
                multiplicity,
            },
            _ => ParseError::Syntax,
        }
    }
}

impl<I: Stream> ParserError<I> for ParseError {
    type Inner = Self;

    fn from_input(_input: &I) -> Self {
        ParseError::Syntax
    }

    fn into_inner(self) -> Result<Self::Inner, Self> {
        Ok(self)
    }
}
