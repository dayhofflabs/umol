//! Error types for TableIR conversions.

use thiserror::Error;

/// Error type for conversions between TableIR types.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum ConversionError {
    #[error("Has extended features")]
    HasExtendedFeatures,
}

/// Error type for lossy/invalid join operations.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum JoinError {
    #[error("CTFile SGroup label collision: {label}")]
    CtfileSgroupCollision { label: u32 },
    #[error("CTFile RGroup label collision: {label}")]
    CtfileRgroupCollision { label: u32 },
    #[error("CX SGroup label collision: {label}")]
    CxSgroupCollision { label: u32 },
    #[error("CX RGroup label collision: {label}")]
    CxRgroupCollision { label: u32 },
}
