//! I/O for umol: molecular file-format parsers (SMILES, CTFile), the TableIR
//! boundary type and its raise to `MoleculeAst`, byte spans, and diagnostics.

pub mod diagnostics;
pub mod io;
pub mod span;
pub mod table_ir;
