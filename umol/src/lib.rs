//! umol - Warning: May have substance!

mod algebra;
mod capability;
mod conversion;
mod entity;
pub mod error;
pub mod logging;
mod model;
mod operation;
mod property;
mod stuff;

pub use algebra::{Aggregate, Ensemble, Process};
pub use capability::Capability;
pub use conversion::{ConversionMetadata, ConvertTo, ConvertToWithMetadata};
pub use entity::{AsEntity, Entity, Relation};
pub use error::{Error, Result};
pub use model::{AsModel, Model};
pub use operation::{ConversionOperation, Operation};
pub use property::{Property, PropertySpec};
pub use stuff::Stuff;

pub use umol_macros::property;
