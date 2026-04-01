//! Error types and handling.

use std::error::Error as StdError;

/// Trait for all umol module-level error types.
///
/// Used as `Box<dyn UmolError>` at cross-module boundaries.
pub trait Error: StdError + Send + Sync + 'static {
    fn as_any(&self) -> &dyn std::any::Any;
}
