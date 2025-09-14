//! Utilities for SMILES parser.

use crate::io::ir::Chirality;

#[derive(Clone, Debug)]
pub enum BracketField {
    Chiral(Chirality),
    H(u32),
    Q(i32),
    Class(u32),
}
