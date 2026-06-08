//! Error types for TableIR conversions.

use thiserror::Error;

/// Error type for conversions between TableIR types.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum ConversionError {
    #[error("Has extended features")]
    HasExtendedFeatures,
}
