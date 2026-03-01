use thiserror::Error;
use umol_data::Element;

use crate::diagnostics::Diagnostic;
use crate::table_ir::bond::BondOrder;

#[derive(Debug, Error, Clone, PartialEq)]
pub enum ResolutionError {
    #[error("Invalid bond order: {0}")]
    InvalidBondOrder(BondOrder),
    #[error("Atom index out of range: {0}")]
    AtomIndexOutOfRange(u32),
    #[error("Invalid atom specification: {0}")]
    InvalidAtomSpec(String),
    #[error("Invalid atom type registry: {0}")]
    InvalidAtomTypeRegistry(String),

    #[error("Molecule has more than one connected component")]
    TopologyDisconnected,
    #[error("Self-loop on bond {0}")]
    TopologySelfLoop(u32),
    #[error("Parallel bonds {0} and {1}")]
    TopologyParallelEdges(u32, u32),
    #[error("Valence violation for element {0:?}: {1}")]
    ValenceViolation(Element, String),
    #[error("No valence match for {0}")]
    ValenceNoMatch(String),
    #[error("Valence ambiguous for {0}")]
    ValenceAmbiguous(String),
}

impl From<ResolutionError> for Diagnostic {
    fn from(_error: ResolutionError) -> Self {
        todo!()
    }
}
