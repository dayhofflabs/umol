//! I/O for umol: molecular file-format parsers (SMILES, CTFile)
//! and the TableIR boundary type.

pub mod ctfile;
#[cfg(feature = "depiction")]
pub mod depict;
pub mod layout;
pub mod smiles;
#[cfg(feature = "depiction")]
mod svg;
pub mod table_ir;
mod utils;
