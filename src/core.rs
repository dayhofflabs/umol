// Core module exports

pub mod entity;
pub mod model;
pub mod algebra;
pub mod instance;
pub mod property;
pub mod conversion;
pub mod error;
pub mod io;

#[cfg(test)]
pub mod testing;

// Re-export commonly used types
pub use entity::{Entity, Relation};
pub use model::{Model, Capability};
pub use algebra::{Ensemble, Aggregate};
pub use instance::{Instance, Operation};
pub use property::Property;
pub use conversion::Conversion;
pub use error::{Error, Result}; 