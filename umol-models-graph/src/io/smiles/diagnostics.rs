//! Diagnostics types for graph-based models

use std::fmt;

use strum::{AsRefStr, EnumIter, IntoEnumIterator};

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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Span {
    pub start: usize,
    pub end: usize,
}

impl Span {
    pub fn new(start: usize, end: usize) -> Self {
        Self { start, end }
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
        write!(
            f,
            "{} [{}:{}] @{}..{}",
            self.message,
            self.category.as_ref(),
            self.code.as_str(),
            self.span.start,
            self.span.end,
        )
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
            span: Span::new(0, 10),
            message: "Invalid token",
            details: None,
        };
        assert_eq!(d.to_string(), "Invalid token [LEXICAL:INVALID_TOKEN] @0..10");
    }
}
