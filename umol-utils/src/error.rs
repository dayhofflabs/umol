//! The workspace-wide error trait.

use std::any::Any;
use std::error::Error as StdError;

/// Trait for all umol module-level error types.
///
/// Used as `Box<dyn UmolError>` at cross-module boundaries.
pub trait UmolError: StdError + Send + Sync + 'static {
    fn as_any(&self) -> &dyn Any;
}

/// Box any `UmolError` at a cross-module boundary, so `?` promotes into
/// `Box<dyn UmolError>`
impl<E: UmolError> From<E> for Box<dyn UmolError> {
    fn from(error: E) -> Self {
        Box::new(error)
    }
}
