//! Error types for the solver.

use thiserror::Error;
use umol_shared::element::Element;
use umol_shared::error::SpinStateError;
use umol_shared::spin::SpinMultiplicity;

use crate::diagnostics::Diagnostic;
use crate::solver::aromaticity::AromaticityError;

impl From<ResolutionError> for Diagnostic {
    fn from(_error: ResolutionError) -> Self {
        todo!()
    }
}

#[derive(Debug, Error, Clone, PartialEq)]
pub enum ValidationError {
    #[error("non-ground value for field '{field}'")]
    NonGround { field: &'static str },
    #[error("invalid spin multiplicity: {0}")]
    InvalidMultiplicity(u8),
    #[error("field '{field}' out of range: {value} not in [{min}, {max}]")]
    OutOfRange {
        field: &'static str,
        value: i64,
        min: i64,
        max: i64,
    },
    #[error("charge {charge} out of bounds for {element}: expected [{min_charge}, {max_charge}]")]
    ChargeOutOfBounds {
        element: Element,
        charge: i8,
        min_charge: i8,
        max_charge: i8,
    },
    #[error(
        "electron invariant mismatch for {element}: inv_o={orbital_invariant}, inv_e={electron_invariant}"
    )]
    ElectronInvariantMismatch {
        element: Element,
        orbital_invariant: i16,
        electron_invariant: i16,
    },
    #[error("spin state is underdetermined")]
    SpinUnderdetermined,
    #[error("{unpaired_electrons} unpaired electrons, {multiplicity} multiplicity incompatible")]
    SpinIncompatible {
        unpaired_electrons: u8,
        multiplicity: SpinMultiplicity,
    },
}

impl From<SpinStateError> for ValidationError {
    fn from(value: SpinStateError) -> Self {
        match value {
            SpinStateError::Underdetermined => ValidationError::SpinUnderdetermined,
            SpinStateError::Incompatible {
                unpaired_electrons,
                multiplicity,
            } => ValidationError::SpinIncompatible {
                unpaired_electrons,
                multiplicity,
            },
            SpinStateError::UnpairedElectronsOutOfRange { unpaired_electrons } => {
                ValidationError::OutOfRange {
                    field: "unpaired_electrons",
                    value: unpaired_electrons as i64,
                    min: 0,
                    max: 254,
                }
            }
            SpinStateError::MultiplicityOutOfRange { multiplicity } => {
                ValidationError::OutOfRange {
                    field: "multiplicity",
                    value: multiplicity as i64,
                    min: 1,
                    max: 255,
                }
            }
            // Parse-related variants — shouldn't reach validation, but map them sensibly
            // TODO: Rethink
            SpinStateError::UnexpectedToken { token } => {
                ValidationError::InvalidMultiplicity(u8::try_from(token as u32).unwrap_or(0))
            }
            SpinStateError::InvalidTag { .. } => ValidationError::NonGround { field: "spin" },
            SpinStateError::DuplicateTag { .. } => ValidationError::NonGround { field: "spin" },
        }
    }
}

#[derive(Debug, Error, Clone, PartialEq)]
pub enum ConfigError {
    #[error("Invalid atom type registry: {0}")]
    InvalidAtomTypeRegistry(String),
    #[error("Invalid valence table: {0}")]
    InvalidValenceTable(String),
}

#[derive(Debug, Error, Clone, PartialEq)]
pub enum ResolutionError {
    #[error("resolution underdetermined")]
    Underdetermined,
    #[error("resolution contradictory")]
    Contradictory,
    #[error(transparent)]
    Aromaticity(#[from] AromaticityError),
}
