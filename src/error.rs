// Error types for molecule validation

use thiserror::Error;
use crate::atom::Element;
use crate::graph::{AtomIndex, BondIndex};

#[derive(Error, Debug)]
pub enum MoleculeError {
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
}
