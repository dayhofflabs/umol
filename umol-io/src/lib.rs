//! I/O for umol: molecular file-format parsers (SMILES, CTFile)
//! and the TableIR boundary type.

pub mod ctfile;
#[cfg(feature = "coordgen")]
pub mod depict;
pub mod layout;
pub mod smiles;
#[cfg(feature = "coordgen")]
mod svg;
pub mod table_ir;
mod utils;
