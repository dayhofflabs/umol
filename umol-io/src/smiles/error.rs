use std::any::Any;

use thiserror::Error;
use umol_shared::error::UmolError;

#[derive(Debug, Clone, PartialEq, Error)]
pub enum ParseError {
    #[error("Leading whitespace")]
    LeadingWhitespace,
    #[error("Invalid element at position {pos}")]
    InvalidElement { pos: usize },
    #[error("Invalid token at position {pos}")]
    InvalidToken { pos: usize },
    #[error("Unbalanced open parenthesis at position {pos}")]
    UnbalancedOpenParen { pos: usize },
    #[error("Unbalanced close parenthesis at position {pos}")]
    UnbalancedCloseParen { pos: usize },
    #[error("Empty branch at position {pos}")]
    EmptyBranch { pos: usize },
    #[error("Empty group at position {pos}")]
    EmptyGroup { pos: usize },
    #[error("Nonfinal group at position {pos}")]
    NonfinalGroup { pos: usize },
    #[error("Leading bond at position {pos}")]
    LeadingBond { pos: usize },
    #[error("Trailing bond at position {pos}")]
    TrailingBond { pos: usize },
    #[error("Consecutive bonds at position {pos}")]
    ConsecutiveBonds { pos: usize },
    #[error("Leading ring at position {pos}")]
    LeadingRing { pos: usize },
    #[error("Unbalanced ring index opening at position {open_pos}")]
    UnbalancedRingIndex { open_pos: usize },
    #[error("Invalid ring index at position {pos}")]
    InvalidRingIndex { pos: usize },
    #[error("Mismatched ring bond orders at position {pos}")]
    MismatchedRingBondOrders { pos: usize, open_pos: usize },
    #[error("Mismatched ring bond directions at position {pos}")]
    MismatchedRingBondDirections { pos: usize, open_pos: usize },
    #[error("Mismatched ring bond donations at position {pos}")]
    MismatchedRingBondDonations { pos: usize, open_pos: usize },
    #[error("Leading dot at position {pos}")]
    LeadingDot { pos: usize },
    #[error("Trailing dot at position {pos}")]
    TrailingDot { pos: usize },
    #[error("Consecutive dots at position {pos}")]
    ConsecutiveDots { pos: usize },
    #[error("Dot before ring at position {pos}")]
    DotBeforeRing { pos: usize },
    #[error("Empty bracket at position {pos}")]
    EmptyBracket { pos: usize },
    #[error("Unbalanced open bracket at position {pos}")]
    UnbalancedOpenBracket { pos: usize },
    #[error("Unbalanced close bracket at position {pos}")]
    UnbalancedCloseBracket { pos: usize },
    #[error("Stray bracket field at position {pos}")]
    StrayBracketField { pos: usize },
    #[error("Duplicate bracket field at position {pos}")]
    DuplicateBracketField { pos: usize },
    #[error("Missing class index at position {pos}")]
    MissingClassIndex { pos: usize },
    #[error("Missing chirality index at position {pos}")]
    MissingChiralityIndex { pos: usize },
    #[error("Chirality out of range at position {pos}")]
    ChiralityOutOfRange { pos: usize },
    #[error("Bracket H with Hcount at position {pos}")]
    BracketHwithHcount { pos: usize },
    #[error("Invalid bracket at position {pos}")]
    InvalidBracket { pos: usize },
    #[error("Invalid CX tag at position {pos}")]
    InvalidCxTag { pos: usize },
    #[error("Missing reaction arrow at position {pos}")]
    MissingReactionArrow { pos: usize },
    #[error("Atom index out of bounds: {atom_idx}")]
    AtomIndexOutOfBounds { atom_idx: u32 },
    #[error("Bond index out of bounds: {bond_idx}")]
    BondIndexOutOfBounds { bond_idx: u32 },
    #[error("Mismatched atom/bond indices: atom {atom_idx} is not incident on bond {bond_idx}")]
    MismatchedAtomBondIndices { atom_idx: u32, bond_idx: u32 },
    #[error("S-group index out of bounds: {sgroup_idx}")]
    SgroupIndexOutOfBounds { sgroup_idx: u32 },
}

impl UmolError for ParseError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}
