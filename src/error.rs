// Error types for molecule validation

#[derive(Debug, thiserror::Error)]
pub enum MoleculeError {
    #[error("Invalid atom index")]
    InvalidAtomIndex,
    #[error("Invalid element symbol: {0}")]
    InvalidElementSymbol(String),
    #[error("Invalid bond index")]
    InvalidBondIndex,
}
