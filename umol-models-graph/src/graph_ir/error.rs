use thiserror::Error;
use umol_data::{Element, SpinStateError};

use crate::diagnostics::Diagnostic;
use crate::graph_ir::aromaticity::AromaticityError;
use crate::graph_ir::atom_type::AtomError;
use crate::graph_ir::bond::BondError;
use crate::graph_ir::kekule::KekulizationError;
use crate::table_ir::bond::BondOrder;

#[derive(Debug, Error, Clone, PartialEq)]
pub enum ResolutionError {
    #[error("Invalid bond order: {0}")]
    InvalidBondOrder(BondOrder),
    #[error("Atom index out of range: {0}")]
    AtomIndexOutOfRange(u32),
    #[error("Invalid atom: {0}")]
    InvalidAtom(String),
    #[error("Invalid bond: {0}")]
    InvalidBond(String),
    #[error("Invalid atom type registry: {0}")]
    InvalidAtomTypeRegistry(String),
    #[error("Invalid valence table: {0}")]
    InvalidValenceTable(String),

    #[error("Molecule has more than one connected component")]
    TopologyDisconnected,
    #[error("Self-loop on bond {0}")]
    TopologySelfLoop(u32),
    #[error("Parallel bonds {0} and {1}")]
    TopologyParallelEdges(u32, u32),
    #[error("Valence violation for element {0:?}: {1}")]
    ValenceViolation(Element, String),
    #[error("Bond invariant violation: {0}")]
    BondInvariantViolation(String),
    #[error("No valence match for {0}")]
    ValenceNoMatch(String),
    #[error("Valence ambiguous for {0}")]
    ValenceAmbiguous(String),

    #[error(transparent)]
    SpinState(#[from] SpinStateError),

    #[error("Aromaticity inconsistent: {0}")]
    AromaticityInconsistent(String),
    #[error(transparent)]
    Aromaticity(#[from] AromaticityError),
    #[error(transparent)]
    Kekulization(#[from] KekulizationError),

    #[error("Molecular charge mismatch: explicit {explicit}, from atoms {atom_sum}")]
    MolecularChargeMismatch { explicit: i8, atom_sum: i8 },
    #[error(
        "Molecular spin incompatible: {explicit_unpaired} unpaired electrons (multiplicity {explicit_multiplicity}) \
         from atoms (total unpaired: {atom_unpaired_sum})"
    )]
    MolecularSpinIncompatible {
        explicit_unpaired: u8,
        explicit_multiplicity: u8,
        atom_unpaired_sum: u16,
    },
    #[error(
        "Molecular spin incomplete: explicit multiplicity is required (compatible multiplicities: {compatible_multiplicities:?}, total unpaired: {atom_unpaired_sum})"
    )]
    MolecularSpinIncomplete {
        atom_unpaired_sum: u16,
        compatible_multiplicities: Vec<u8>,
    },
}

#[derive(Debug, Error, Clone, PartialEq)]
pub enum ValidationError {
    #[error("non-ground value for field '{field}'")]
    NonGround { field: &'static str },
    #[error("invalid spin multiplicity: {0}")]
    InvalidMultiplicity(u8),
    #[error(transparent)]
    Atom(#[from] AtomError),
    #[error(transparent)]
    Bond(#[from] BondError),
}

impl From<ValidationError> for ResolutionError {
    fn from(value: ValidationError) -> Self {
        match value {
            ValidationError::NonGround { field } => {
                ResolutionError::InvalidBond(format!("non-ground field: {field}"))
            }
            ValidationError::InvalidMultiplicity(n) => {
                ResolutionError::InvalidBond(format!("invalid multiplicity: {n}"))
            }
            ValidationError::Atom(ae) => ResolutionError::InvalidAtom(ae.to_string()),
            ValidationError::Bond(be) => ResolutionError::InvalidBond(be.to_string()),
        }
    }
}

impl From<ResolutionError> for Diagnostic {
    fn from(_error: ResolutionError) -> Self {
        todo!()
    }
}
