//! Error types for graph_ir module.

use thiserror::Error;
use umol_shared::{Element, SpinMultiplicity, SpinStateError};

use super::aromaticity::AromaticityError;
use super::kekule::KekulizationError;
use super::transform::TransformError;
use crate::diagnostics::Diagnostic;
use crate::table_ir::bond::BondOrder;

#[derive(Debug, Error, Clone, PartialEq)]
pub enum GraphIrError {
    #[error(transparent)]
    Validation(#[from] ValidationError),
    #[error(transparent)]
    Resolution(#[from] ResolutionError),
    #[error(transparent)]
    Aromaticity(#[from] AromaticityError),
    #[error(transparent)]
    Kekulization(#[from] KekulizationError),
    #[error(transparent)]
    Transform(#[from] TransformError),
}

impl From<ResolutionError> for Diagnostic {
    fn from(_error: ResolutionError) -> Self {
        todo!()
    }
}

#[derive(Debug, Error, Clone, PartialEq)]
pub enum ValidationError {
    #[error("non-ground value for field '{field}'")]
    NonGround { field: &'static str },
    #[error("invalid spin multiplicity: {0}")]
    InvalidMultiplicity(u8),
    #[error("field '{field}' out of range: {value} not in [{min}, {max}]")]
    OutOfRange {
        field: &'static str,
        value: i64,
        min: i64,
        max: i64,
    },
    #[error("charge {charge} out of bounds for {element}: expected [{min_charge}, {max_charge}]")]
    ChargeOutOfBounds {
        element: Element,
        charge: i8,
        min_charge: i8,
        max_charge: i8,
    },
    #[error(
        "electron invariant mismatch for {element}: inv_o={orbital_invariant}, inv_e={electron_invariant}"
    )]
    ElectronInvariantMismatch {
        element: Element,
        orbital_invariant: i16,
        electron_invariant: i16,
    },
    #[error("spin state is underdetermined")]
    SpinUnderdetermined,
    #[error("{unpaired_electrons} unpaired electrons, {multiplicity} multiplicity incompatible")]
    SpinIncompatible {
        unpaired_electrons: u8,
        multiplicity: SpinMultiplicity,
    },
}

impl From<SpinStateError> for ValidationError {
    fn from(value: SpinStateError) -> Self {
        match value {
            SpinStateError::Underdetermined => ValidationError::SpinUnderdetermined,
            SpinStateError::Incompatible {
                unpaired_electrons,
                multiplicity,
            } => ValidationError::SpinIncompatible {
                unpaired_electrons,
                multiplicity,
            },
            SpinStateError::UnpairedElectronsOutOfRange { unpaired_electrons } => {
                ValidationError::OutOfRange {
                    field: "unpaired_electrons",
                    value: unpaired_electrons as i64,
                    min: 0,
                    max: 254,
                }
            }
            SpinStateError::MultiplicityOutOfRange { multiplicity } => {
                ValidationError::OutOfRange {
                    field: "multiplicity",
                    value: multiplicity as i64,
                    min: 1,
                    max: 255,
                }
            }
            // Parse-related variants — shouldn't reach validation, but map them sensibly
            // TODO: Rethink
            SpinStateError::UnexpectedToken { token } => {
                ValidationError::InvalidMultiplicity(u8::try_from(token as u32).unwrap_or(0))
            }
            SpinStateError::InvalidTag { .. } => ValidationError::NonGround { field: "spin" },
            SpinStateError::DuplicateTag { .. } => ValidationError::NonGround { field: "spin" },
        }
    }
}

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
    #[error("Aromaticity inconsistent: {0}")]
    AromaticityInconsistent(String),
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
