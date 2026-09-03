//! External-format boundaries, two-dimensional layout, and depiction for umol.
//!
//! The `coordgen` feature enables explicit CoordGen layout through [`layout`]. The `depiction`
//! feature enables the high-level `depict` API and currently selects CoordGen as its default and
//! only layout backend.

pub mod ctfile;
#[cfg(feature = "depiction")]
pub mod depict;
pub mod layout;
pub mod smiles;
#[cfg(feature = "depiction")]
mod svg;
pub mod table_ir;
mod utils;
