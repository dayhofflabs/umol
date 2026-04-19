use thiserror::Error;

use crate::diagnostics::{Diagnostic, DiagnosticKind, Severity};
use crate::ops::error::ResolutionError;
use crate::span::Span;

#[derive(Debug, Clone, PartialEq, Error)]
pub enum SmilesError {
    #[error(transparent)]
    Parse(#[from] ParseError),
    #[error(transparent)]
    Resolve(#[from] ResolutionError),
    #[error("resolution contradictory")]
    ResolveContradictory,
    #[error("resolution underdetermined")]
    ResolveUnderdetermined,
}

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
    MismatchedRingBondDirs { pos: usize, open_pos: usize },
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

impl From<ParseError> for Diagnostic {
    fn from(error: ParseError) -> Self {
        use DiagnosticKind::*;
        use ParseError::*;

        let (kind, span) = match error {
            LeadingWhitespace => (
                SmilesLeadingWhitespace,
                Span::from_bytes_opt(Some(0), Some(1)),
            ),
            InvalidElement { pos } => (
                SmilesInvalidElement,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            InvalidToken { pos } => (
                SmilesInvalidToken,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            UnbalancedOpenParen { pos } => (
                SmilesUnbalancedOpenParen,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            UnbalancedCloseParen { pos } => (
                SmilesUnbalancedCloseParen,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            EmptyBranch { pos } => (
                SmilesEmptyBranch,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            EmptyGroup { pos } => (
                SmilesEmptyGroup,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            NonfinalGroup { pos } => (
                SmilesNonfinalGroup,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            LeadingBond { pos } => (
                SmilesLeadingBond,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            TrailingBond { pos } => (
                SmilesTrailingBond,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            ConsecutiveBonds { pos } => (
                SmilesConsecutiveBonds,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            LeadingRing { pos } => (
                SmilesLeadingRing,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            UnbalancedRingIndex { open_pos } => (
                SmilesUnbalancedRingIndex,
                Span::from_bytes_opt(Some(open_pos as u32), Some(open_pos as u32 + 1)),
            ),
            InvalidRingIndex { pos } => (
                SmilesInvalidRingIndex,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            MismatchedRingBondOrders { pos, .. } => (
                SmilesMismatchedRingBondOrders,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            MismatchedRingBondDirs { pos, .. } => (
                SmilesMismatchedRingBondDirs,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            MismatchedRingBondDonations { pos, .. } => (
                SmilesMismatchedRingBondDonations,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            LeadingDot { pos } => (
                SmilesLeadingDot,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            TrailingDot { pos } => (
                SmilesTrailingDot,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            ConsecutiveDots { pos } => (
                SmilesConsecutiveDots,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            DotBeforeRing { pos } => (
                SmilesDotBeforeRing,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            EmptyBracket { pos } => (
                SmilesEmptyBracket,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            UnbalancedOpenBracket { pos } => (
                SmilesUnbalancedOpenBracket,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            UnbalancedCloseBracket { pos } => (
                SmilesUnbalancedCloseBracket,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            StrayBracketField { pos } => (
                SmilesStrayBracketField,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            DuplicateBracketField { pos } => (
                SmilesDuplicateBracketField,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            MissingClassIndex { pos } => (
                SmilesMissingClassIndex,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            MissingChiralityIndex { pos } => (
                SmilesMissingChiralityIndex,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            ChiralityOutOfRange { pos } => (
                SmilesChiralityOutOfRange,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            BracketHwithHcount { pos } => (
                SmilesBracketHwithHcount,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            InvalidBracket { pos } => (
                SmilesInvalidBracket,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            InvalidCxTag { pos } => (
                SmilesInvalidCxTag,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            MissingReactionArrow { pos } => (
                SmilesMissingReactionArrow,
                Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
            ),
            AtomIndexOutOfBounds { .. } => (SmilesAtomIndexOutOfBounds, None),
            BondIndexOutOfBounds { .. } => (SmilesBondIndexOutOfBounds, None),
            MismatchedAtomBondIndices { .. } => (SmilesMismatchedAtomBondIndices, None),
            SgroupIndexOutOfBounds { .. } => (SmilesInvalidCxTag, None),
        };

        Diagnostic {
            kind,
            category: kind.category(),
            severity: Severity::Error,
            span,
            details: None,
        }
    }
}
