//! Error types for umol-shared.

use std::any::Any;
use std::error::Error as StdError;

use thiserror::Error;

use crate::spin::SpinMultiplicity;

/// Trait for all umol module-level error types.
///
/// Used as `Box<dyn UmolError>` at cross-module boundaries.
pub trait UmolError: StdError + Send + Sync + 'static {
    fn as_any(&self) -> &dyn Any;
}

#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum DataError {
    #[error(transparent)]
    Element(#[from] ElementError),
    #[error(transparent)]
    Isotope(#[from] IsotopeError),
    #[error(transparent)]
    Occupation(#[from] OccupationError),
    #[error(transparent)]
    SpinState(#[from] SpinStateError),
}

impl UmolError for DataError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum ElementError {
    #[error("invalid element symbol: {symbol}")]
    InvalidSymbol { symbol: String },
    #[error("invalid atomic number: {atomic_number}")]
    InvalidAtomicNumber { atomic_number: u8 },
}

impl UmolError for ElementError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum IsotopeError {
    #[error("invalid isotope: {symbol}")]
    InvalidSymbol { symbol: String },
}

impl UmolError for IsotopeError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum OccupationError {
    #[error("invalid occupation: {occupation}")]
    Invalid { occupation: String },
}

impl UmolError for OccupationError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum SpinStateError {
    #[error("unexpected token '{token}'")]
    UnexpectedToken { token: char },
    #[error("invalid tag '{tag}'")]
    InvalidTag { tag: String },
    #[error("duplicate tag '{tag}'")]
    DuplicateTag { tag: String },
    #[error("unpaired electrons {unpaired} out of range")]
    UnpairedElectronsOutOfRange { unpaired: u8 },
    #[error("multiplicity {multiplicity} out of range")]
    MultiplicityOutOfRange { multiplicity: u8 },
    #[error("spin state is underdetermined")]
    Underdetermined,
    #[error("{unpaired} unpaired electrons, {multiplicity} multiplicity incompatible")]
    Incompatible {
        unpaired: u8,
        multiplicity: SpinMultiplicity,
    },
}

impl UmolError for SpinStateError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}
