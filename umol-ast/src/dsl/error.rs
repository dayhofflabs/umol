//! Domain errors for DSL parsing.

use thiserror::Error;
use umol_shared::error::SpinStateError;
use umol_shared::spin::SpinMultiplicity;
use winnow::error::{ErrMode, ParserError};
use winnow::stream::Stream;

pub type WinnowErrorMode = ErrMode<ParseError>;
pub type PResult<T> = Result<T, WinnowErrorMode>;


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
    #[error("unknown aromatic predicate: {0}")]
    UnknownAromaticPredicate(String),
    #[error("duplicate aromatic predicate: {0}")]
    DuplicateAromaticPredicate(String),
    #[error("trailing input: {0:?}")]
    TrailingInput(String),
    #[error("unpaired electrons {unpaired} out of range")]
    UnpairedElectronsOutOfRange { unpaired: u8 },
    #[error("multiplicity {multiplicity} out of range")]
    MultiplicityOutOfRange { multiplicity: u8 },
    #[error("{unpaired} unpaired electrons, {multiplicity} multiplicity incompatible")]
    IncompatibleSpin { unpaired: u8, multiplicity: SpinMultiplicity },
    #[error("syntax error")]
    Syntax,
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
