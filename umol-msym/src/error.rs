use std::error::Error as StdError;
use std::ffi::CStr;
use std::fmt;

use umol_msym_sys as ffi;

#[derive(Debug, thiserror::Error)]
#[error("{message}")]
pub struct MsymError {
    pub code: i32,
    pub message: String,
}

impl MsymError {
    pub(crate) fn from_code(code: ffi::msym_error_t) -> Self {
        let message = unsafe {
            let ptr = ffi::msymErrorString(code);
            if ptr.is_null() {
                "unknown error".to_string()
            } else {
                CStr::from_ptr(ptr).to_string_lossy().into_owned()
            }
        };
        let detail = unsafe {
            let ptr = ffi::msymGetErrorDetails();
            if ptr.is_null() {
                None
            } else {
                let s = CStr::from_ptr(ptr).to_string_lossy();
                if s.is_empty() {
                    None
                } else {
                    Some(s.into_owned())
                }
            }
        };
        let message = match detail {
            Some(d) => format!("{message}: {d}"),
            None => message,
        };
        Self { code, message }
    }
}

pub(crate) fn check(code: ffi::msym_error_t) -> Result<(), MsymError> {
    if code == ffi::MSYM_SUCCESS {
        Ok(())
    } else {
        Err(MsymError::from_code(code))
    }
}

#[derive(Debug, Clone)]
pub struct ParseError(pub(crate) String);

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid symbol: '{}'", self.0)
    }
}

impl StdError for ParseError {}

#[derive(Debug, Clone)]
pub enum ReductionError {
    InfiniteGroup,
    DimensionMismatch { expected: usize, got: usize },
    NonIntegralMultiplicity { irrep: String, value: f64 },
}

impl fmt::Display for ReductionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InfiniteGroup => {
                write!(f, "character-based reduction undefined for infinite groups")
            }
            Self::DimensionMismatch { expected, got } => {
                write!(f, "expected {expected} class characters, got {got}")
            }
            Self::NonIntegralMultiplicity { irrep, value } => {
                write!(f, "non-integral multiplicity {value:.4} for irrep {irrep}")
            }
        }
    }
}

impl StdError for ReductionError {}
