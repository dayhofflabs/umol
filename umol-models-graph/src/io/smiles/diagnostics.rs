//! Diagnostics types for graph-based models

use std::fmt;

use strum::{AsRefStr, EnumIter, IntoEnumIterator};
use umol_data::Element;

use crate::simple_ir::{BondDir, BondOrder, BondStereo, BondSymbol, Chirality};
use crate::span::Span;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, EnumIter, AsRefStr)]
#[strum(serialize_all = "SCREAMING_SNAKE_CASE")]
pub enum Code {
    // Lexical errors
    InvalidWhitespace,
    InvalidComment,
    UnterminatedBlockComment,
    InvalidElement,
    InvalidToken,

    // Syntactic errors
    UnbalancedOpenParen,
    UnbalancedCloseParen,
    EmptyBranch,
    EmptyGroup,
    NonfinalGroup,
    LeadingBond,
    TrailingBond,
    ConsecutiveBonds,
    LeadingRing,
    UnbalancedRingIndex,
    InvalidRingIndex,
    MismatchedRingBondDirs,
    MismatchedRingBondOrders,
    LeadingDot,
    TrailingDot,
    ConsecutiveDots,
    DotBeforeRing,
    EmptyBracket,
    UnbalancedOpenBracket,
    UnbalancedCloseBracket,
    StrayBracketField,
    DuplicateBracketField,
    MissingClassIndex,
    MissingChiralityIndex,
    BracketHwithHcount,
    InvalidBracket,

    // Topology errors
    SelfLoopRing,
    ParallelEdges,

    // Valence errors
    ValenceOutOfElementRange,
    HcountOutOfElementRange,
    ChargeOutOElementfRange,
    HcountMismatch,
    NoMatch,
    AmbiguousMatch,
    NoKnownValenceStates,
    ValenceUnknownBondOrder,
    MissingBracketH,

    // Aromaticity errors
    AromaticAtomNotInRing,
    AromaticBondNotInRing,
    NoMatchingAromaticAtomConfig,
    InvalidAromaticAtom,
    InvalidAromaticBondAtom,
    AromaticBondOrderMismatch,
    KekuleInconsistent,
    HuckelFail,

    // Aromaticity warnings
    AvoidMixedAromaticity,
    AvoidInconsistentAromaticity,
    HuckelInconsistent,

    // Stereochemistry errors
    DoubleConflict,
    DoubleInsufficient,

    // Stereochemistry warnings
    AvoidUnnecessaryStereoDescriptor,
    UnsupportedCentralChiralityElement,
    ChiralitySubstituentMismatch,
    NonChiralAnnotated,

    // Numeric errors
    Overflow,
    ClassOutOfRange,
    HcountOutOfRange,
    ChargeOutOfRange,
    IsotopeOutOfRange,
    ChiralityOutOfRange,

    // Numeric warnings
    IsotopeUncatalogued,

    // Internal errors
    InternalError,

    // Style warnings
    PreferBareOrganicAtom,
    PreferImplicitH,
    PreferBracketFieldOrder,
    PreferSimpleChargeSign,
    PreferSimpleHcount,
    AvoidExplicitSingleBond,
    AvoidExplicitAromaticBond,
    AvoidUnnecessaryGroup,
    AvoidRedundantNestedParens,
    PreferBranchesBeforeRingBonds,
    PreferFirstRingOne,
    PreferConsecutiveRingNumbering,
    AvoidReusedRingIndices,
    PreferSingleDigitRingIndex,
    PreferSingleRingClosure,
    AvoidAdjacentRingClosures,
    PreferBondSymbolAtRingOpen,
    AvoidRingClosureAcrossDot,
    PreferAromaticForm,
}

impl Code {
    // Iterate all enum variants
    pub fn all() -> impl Iterator<Item = Code> {
        Code::iter()
    }
    // Stable string name for this code (SCREAMING_SNAKE_CASE variant name)
    pub fn as_str(&self) -> &str {
        self.as_ref()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Severity {
    Error,
    Warning,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, AsRefStr, EnumIter)]
#[strum(serialize_all = "UPPERCASE")]
pub enum Category {
    Lexical,
    Syntactic,
    Topology,
    Valence,
    Aromaticity,
    Stereo,
    Style,
    Internal,
}

impl Category {
    pub fn all() -> impl Iterator<Item = Category> {
        Category::iter()
    }

    pub fn as_str(&self) -> &str {
        self.as_ref()
    }

    pub fn default_category(&self) -> Category {
        Category::Internal
    }
    pub fn default_severity(&self) -> Severity {
        Severity::Error
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Diagnostic {
    pub code: Code,
    pub severity: Severity,
    pub category: Category,
    pub span: Span,
    pub message: &'static str,
    pub details: Option<String>,
}

impl Diagnostic {}

#[derive(Debug, Default, Clone)]
pub struct DiagnosticList {
    pub diagnostics: Vec<Diagnostic>,
}

impl DiagnosticList {
    pub fn new() -> Self {
        Self {
            diagnostics: Vec::new(),
        }
    }
    pub fn push(&mut self, d: Diagnostic) {
        self.diagnostics.push(d);
    }
    pub fn extend<I: IntoIterator<Item = Diagnostic>>(&mut self, it: I) {
        self.diagnostics.extend(it);
    }
    pub fn append_list(&mut self, other: &mut DiagnosticList) {
        self.diagnostics.append(&mut other.diagnostics);
    }
    pub fn iter(&self) -> impl Iterator<Item = &Diagnostic> {
        self.diagnostics.iter()
    }

    pub fn into_vec(self) -> Vec<Diagnostic> {
        self.diagnostics
    }
    pub fn from(d: Diagnostic) -> Self {
        let mut list = DiagnosticList::new();
        list.push(d);
        list
    }
    pub fn has_errors(&self) -> bool {
        self.diagnostics
            .iter()
            .any(|d| d.severity == Severity::Error)
    }

    pub fn upgrade_warnings(&mut self) {
        for d in &mut self.diagnostics {
            if d.severity == Severity::Warning {
                d.severity = Severity::Error;
            }
        }
    }

    pub fn errors(&self) -> impl Iterator<Item = &Diagnostic> {
        self.diagnostics
            .iter()
            .filter(|d| d.severity == Severity::Error)
    }
    pub fn warnings(&self) -> impl Iterator<Item = &Diagnostic> {
        self.diagnostics
            .iter()
            .filter(|d| d.severity == Severity::Warning)
    }
}

impl IntoIterator for DiagnosticList {
    type Item = Diagnostic;
    type IntoIter = std::vec::IntoIter<Diagnostic>;
    fn into_iter(self) -> Self::IntoIter {
        self.diagnostics.into_iter()
    }
}

impl fmt::Display for Diagnostic {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (start, end) = self.span.bytes_range();
        write!(
            f,
            "{} [{}:{}] @{}..{}",
            self.message,
            self.category.as_ref(),
            self.code.as_str(),
            start,
            end,
        )
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum Edit {
    SetAtomCharge {
        atom: usize,
        charge: i32,
    },
    SetAtomExplicitHCount {
        atom: usize,
        count: u32,
    },
    SetAtomImplicitHCount {
        atom: usize,
        count: u32,
    },
    SetAtomImplicitH {
        atom: usize,
        implicit: bool,
    },
    SetAtomAromaticFlag {
        atom: usize,
        aromatic: Option<bool>,
    },
    SetAtomChirality {
        atom: usize,
        chirality: Option<Chirality>,
    },
    SetAtomClass {
        atom: usize,
        class: Option<u32>,
    },
    SetAtomUnpairedECount {
        atom: usize,
        count: u32,
    },
    SetBondOrder {
        bond: usize,
        order: BondOrder,
    },
    SetBondSymbol {
        bond: usize,
        symbol: BondSymbol,
    },
    SetBondDirection {
        bond: usize,
        direction: Option<BondDir>,
    },
    SetBondStereo {
        bond: usize,
        stereo: Option<BondStereo>,
    },
    AddAtom {
        atom: usize,
        element: Element,
    },
    RemoveAtom {
        atom: usize,
    },
    AddBond {
        bond: usize,
        atoms: (usize, usize),
        order: BondOrder,
    },
    RemoveBond {
        bond: usize,
    },
    RetargetBond {
        bond: usize,
        atoms: (usize, usize),
    },
    SetBondRing {
        bond: usize,
        ring: Option<u32>,
    },
}

#[derive(Debug, Default, Clone)]
pub struct EditList {
    pub edits: Vec<Edit>,
}

impl EditList {
    pub fn new() -> Self {
        Self { edits: Vec::new() }
    }

    pub fn push(&mut self, edit: Edit) {
        self.edits.push(edit);
    }

    pub fn extend<I: IntoIterator<Item = Edit>>(&mut self, edits: I) {
        self.edits.extend(edits);
    }

    pub fn append_list(&mut self, other: &mut EditList) {
        self.edits.append(&mut other.edits);
    }

    pub fn iter(&self) -> impl Iterator<Item = &Edit> {
        self.edits.iter()
    }

    pub fn into_vec(self) -> Vec<Edit> {
        self.edits
    }
}

impl IntoIterator for EditList {
    type Item = Edit;
    type IntoIter = std::vec::IntoIter<Edit>;

    fn into_iter(self) -> Self::IntoIter {
        self.edits.into_iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_diagnostic_display() {
        let d = Diagnostic {
            code: Code::InvalidToken,
            severity: Severity::Error,
            category: Category::Lexical,
            span: Span::bytes(0, 10),
            message: "Invalid token",
            details: None,
        };
        assert_eq!(
            d.to_string(),
            "Invalid token [LEXICAL:INVALID_TOKEN] @0..10"
        );
    }
}
