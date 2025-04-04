// Error types

use crate::graph::{AtomIndex, BondIndex};
use crate::Element;
use thiserror;
use std::fmt;

#[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
pub enum Error {
    #[error("Invalid element symbol: {0}")]
    InvalidElementSymbol(String),
    
    #[error("Invalid charge {charge} for element {element:?} (allowed range: {min} to {max})")]
    InvalidCharge {
        element: Element,
        charge: i8,
        min: i8,
        max: i8,
    },
    
    #[error("Invalid number of unpaired electrons {unpaired} for element {element:?} (max allowed: {max})")]
    InvalidUnpairedElectrons {
        element: Element,
        unpaired: u8,
        max: u8,
    },
    
    #[error("Invalid atom index {0}")]
    InvalidAtomIndex(AtomIndex),
    
    #[error("Invalid bond order: {0}")]
    InvalidBondOrder(String),
    
    #[error("Invalid bond index {0}")]
    InvalidBondIndex(BondIndex),
    
    #[error("Invalid operation: {0}")]
    InvalidOperation(String),
    
    #[error("Invalid valence for {element:?} at index {atom_idx:?}, expected {expected} bonds but has {actual}")]
    InvalidValence {
        element: Element,
        atom_idx: AtomIndex,
        expected: u8,
        actual: u8,
    },
    
    #[error("Molecule contains disconnected fragments")]
    DisconnectedMolecule,
    
    #[error("Chemically unreasonable structure: {0}")]
    UnreasonableStructure(String),
    
    #[error("Validation failed with multiple errors")]
    ValidationFailed(Vec<String>),
}