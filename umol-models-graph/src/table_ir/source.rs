//! Input molecular formats for TableIR.

/// Input molecular format
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SourceFormat {
    MOL,
    SMILES,
    SMARTS,
    UNKNOWN,
}
