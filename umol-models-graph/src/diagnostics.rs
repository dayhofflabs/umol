//! Diagnostics for atom/bond-based molecule models

use std::fmt::{self, Display};

use serde::{Deserialize, Serialize};
use strum::{Display, EnumIter, EnumMessage, IntoEnumIterator};

use crate::span::Span;

#[derive(
    Debug,
    Copy,
    Clone,
    Default,
    PartialEq,
    Eq,
    Display,
    EnumIter,
    EnumMessage,
    Serialize,
    Deserialize,
)]
#[strum(serialize_all = "UPPERCASE")]
pub enum Severity {
    #[default]
    Error,
    Warning,
}

#[derive(
    Debug,
    Copy,
    Clone,
    Default,
    PartialEq,
    Eq,
    Display,
    EnumIter,
    Serialize,
    Deserialize,
)]
#[strum(serialize_all = "SCREAMING_SNAKE_CASE")]
pub enum Category {
    #[default]
    Lexical,
    Syntactic,
    Semantic,
}

#[derive(
    Debug,
    Copy,
    Clone,
    PartialEq,
    Eq,
    Display,
    EnumIter,
    EnumMessage,
    Serialize,
    Deserialize,
    Default,
)]
#[strum(serialize_all = "SCREAMING_SNAKE_CASE")]
pub enum DiagnosticKind {
    #[strum(message = "Invalid whitespace")]
    SmilesInvalidWhitespace,
    #[strum(message = "Invalid comment")]
    SmilesInvalidComment,
    #[strum(message = "Unterminated block comment")]
    SmilesUnterminatedBlockComment,
    #[strum(message = "Invalid element")]
    SmilesInvalidElement,
    #[strum(message = "Invalid token")]
    SmilesInvalidToken,

    #[strum(message = "Unbalanced open parenthesis")]
    SmilesUnbalancedOpenParen,
    #[strum(message = "Unbalanced close parenthesis")]
    SmilesUnbalancedCloseParen,
    #[strum(message = "Empty branch")]
    SmilesEmptyBranch,
    #[strum(message = "Empty group")]
    SmilesEmptyGroup,
    #[strum(message = "Nonfinal group")]
    SmilesNonfinalGroup,
    #[strum(message = "Leading bond")]
    SmilesLeadingBond,
    #[strum(message = "Trailing bond")]
    SmilesTrailingBond,
    #[strum(message = "Consecutive bonds")]
    SmilesConsecutiveBonds,
    #[strum(message = "Leading ring")]
    SmilesLeadingRing,
    #[strum(message = "Unbalanced ring index")]
    SmilesUnbalancedRingIndex,
    #[strum(message = "Invalid ring index")]
    SmilesInvalidRingIndex,
    #[strum(message = "Mismatched ring bond directions")]
    SmilesMismatchedRingBondDirs,
    #[strum(message = "Mismatched ring bond orders")]
    SmilesMismatchedRingBondOrders,
    #[strum(message = "Leading dot")]
    SmilesLeadingDot,
    #[strum(message = "Trailing dot")]
    SmilesTrailingDot,
    #[strum(message = "Consecutive dots")]
    SmilesConsecutiveDots,
    #[strum(message = "Dot before ring")]
    SmilesDotBeforeRing,
    #[strum(message = "Empty bracket")]
    SmilesEmptyBracket,
    #[strum(message = "Unbalanced open bracket")]
    SmilesUnbalancedOpenBracket,
    #[strum(message = "Unbalanced close bracket")]
    SmilesUnbalancedCloseBracket,
    #[strum(message = "Stray bracket field")]
    SmilesStrayBracketField,
    #[strum(message = "Duplicate bracket field")]
    SmilesDuplicateBracketField,
    #[strum(message = "Missing class index")]
    SmilesMissingClassIndex,
    #[strum(message = "Missing chirality index")]
    SmilesMissingChiralityIndex,
    #[strum(message = "Chirality out of range")]
    SmilesChiralityOutOfRange,
    #[strum(message = "Bracket H with Hcount")]
    SmilesBracketHwithHcount,
    #[strum(message = "Invalid bracket")]
    SmilesInvalidBracket,

    #[strum(message = "Class out of range")]
    SmilesClassOutOfRange,
    #[strum(message = "Hcount out of range")]
    SmilesHcountOutOfRange,
    #[strum(message = "Charge out of range")]
    SmilesChargeOutOfRange,
    #[strum(message = "Isotope out of range")]
    SmilesIsotopeOutOfRange,

    #[strum(message = "Isotope uncatalogued")]
    SmilesIsotopeUncatalogued,

    #[strum(message = "Prefer bare organic atom")]
    SmilesPreferBareOrganicAtom,
    #[strum(message = "Prefer implicit H")]
    SmilesPreferImplicitH,
    #[strum(message = "Prefer bracket field order")]
    SmilesPreferBracketFieldOrder,
    #[strum(message = "Prefer simple charge sign")]
    SmilesPreferSimpleChargeSign,
    #[strum(message = "Prefer simple H count")]
    SmilesPreferSimpleHcount,
    #[strum(message = "Avoid explicit single bond")]
    SmilesAvoidExplicitSingleBond,
    #[strum(message = "Avoid explicit aromatic bond")]
    SmilesAvoidExplicitAromaticBond,
    #[strum(message = "Avoid unnecessary group")]
    SmilesAvoidUnnecessaryGroup,
    #[strum(message = "Avoid redundant nested parentheses")]
    SmilesAvoidRedundantNestedParens,
    #[strum(message = "Prefer branches before ring bonds")]
    SmilesPreferBranchesBeforeRingBonds,
    #[strum(message = "Prefer start ring numbering with one")]
    SmilesPreferFirstRingOne,
    #[strum(message = "Prefer consecutive ring numbering")]
    SmilesPreferConsecutiveRingNumbering,
    #[strum(message = "Avoid reused ring indices")]
    SmilesAvoidReusedRingIndices,
    #[strum(message = "Prefer single digit ring index")]
    SmilesPreferSingleDigitRingIndex,
    #[strum(message = "Prefer single ring closure")]
    SmilesPreferSingleRingClosure,
    #[strum(message = "Avoid adjacent ring closures")]
    SmilesAvoidAdjacentRingClosures,
    #[strum(message = "Prefer bond symbol at ring open")]
    SmilesPreferBondSymbolAtRingOpen,
    #[strum(message = "Avoid ring closure across dot")]
    SmilesAvoidRingClosureAcrossDot,
    #[strum(message = "Prefer aromatic form")]
    SmilesPreferAromaticForm,

    #[strum(message = "Invalid counts line")]
    CtfileInvalidCountsLine,
    #[strum(message = "Invalid atom line")]
    CtfileInvalidAtomLine,
    #[strum(message = "Invalid bond line")]
    CtfileInvalidBondLine,
    #[strum(message = "Invalid property line")]
    CtfileInvalidPropertyLine,
    #[strum(message = "Invalid Sgroup line")]
    CtfileInvalidSgroupLine,
    #[strum(message = "Invalid header block")]
    CtfileInvalidHeader,
    #[strum(message = "Invalid SDF data header")]
    CtfileInvalidSdfDataHeader,
    #[strum(message = "Invalid SDF data value")]
    CtfileInvalidSdfDataValue,
    #[strum(message = "Missing record delimiter")]
    CtfileMissingDelimiter,

    #[strum(message = "GraphIR conversion failed")]
    GraphConversionUnknown,

    #[strum(message = "Self-loop ring")]
    GraphTopologySelfLoopRing,
    #[strum(message = "Parallel edges")]
    GraphTopologyParallelEdges,
    #[strum(message = "Unknown topology error")]
    GraphTopologyUnknown,

    #[strum(message = "Out of element range")]
    GraphValenceOutOfElementRange,
    #[strum(message = "H count out of element range")]
    GraphValenceHcountOutOfElementRange,
    #[strum(message = "Charge out of element range")]
    GraphValenceChargeOutOfElementRange,
    #[strum(message = "H count mismatch")]
    GraphValenceHcountMismatch,
    #[strum(message = "No match")]
    GraphValenceNoMatch,
    #[strum(message = "Ambiguous match")]
    GraphValenceAmbiguousMatch,
    #[strum(message = "No known valence states")]
    GraphValenceNoKnownValenceStates,
    #[strum(message = "Valence unknown bond order")]
    GraphValenceUnknownBondOrder,
    #[strum(message = "Missing bracket H")]
    GraphValenceMissingBracketH,
    #[strum(message = "Unknown valence error")]
    GraphValenceUnknown,

    #[strum(message = "Aromatic atom not in ring")]
    GraphAromaticityAromaticAtomNotInRing,
    #[strum(message = "Aromatic bond not in ring")]
    GraphAromaticityAromaticBondNotInRing,
    #[strum(message = "No matching aromatic atom config")]
    GraphAromaticityNoMatchingAromaticAtomConfig,
    #[strum(message = "Invalid aromatic atom")]
    GraphAromaticityInvalidAromaticAtom,
    #[strum(message = "Invalid aromatic bond atom")]
    GraphAromaticityInvalidAromaticBondAtom,
    #[strum(message = "Aromatic bond order mismatch")]
    GraphAromaticityAromaticBondOrderMismatch,
    #[strum(message = "Kekule inconsistent")]
    GraphAromaticityKekuleInconsistent,
    #[strum(message = "Huckel fail")]
    GraphAromaticityHuckelFail,
    #[strum(message = "Unknown aromaticity error")]
    GraphAromaticityUnknown,

    #[strum(message = "Avoid mixed aromaticity")]
    GraphAromaticityWarningAvoidMixedAromaticity,
    #[strum(message = "Avoid inconsistent aromaticity")]
    GraphAromaticityWarningAvoidInconsistentAromaticity,
    #[strum(message = "Huckel inconsistent")]
    GraphAromaticityWarningHuckelInconsistent,
    #[strum(message = "Unknown aromaticity warning")]
    GraphAromaticityWarningUnknown,

    #[strum(message = "Double conflict")]
    GraphStereoDoubleConflict,
    #[strum(message = "Double insufficient")]
    GraphStereoDoubleInsufficient,
    #[strum(message = "Unknown stereo error")]
    GraphStereoUnknown,

    #[strum(message = "Avoid unnecessary stereo descriptor")]
    GraphStereoWarningAvoidUnnecessaryStereoDescriptor,
    #[strum(message = "Unsupported central chirality element")]
    GraphStereoWarningUnsupportedCentralChiralityElement,
    #[strum(message = "Chirality substituent mismatch")]
    GraphStereoWarningChiralitySubstituentMismatch,
    #[strum(message = "Non chiral annotated")]
    GraphStereoWarningNonChiralAnnotated,
    #[strum(message = "Unknown stereo warning")]
    GraphStereoWarningUnknown,

    #[default]
    #[strum(message = "Unknown diagnostic")]
    Unknown,
}

impl DiagnosticKind {
    pub fn all() -> impl Iterator<Item = DiagnosticKind> {
        DiagnosticKind::iter()
    }

    pub fn category(&self) -> Category {
        use Category::*;
        use DiagnosticKind::*;

        match self {
            SmilesInvalidWhitespace
            | SmilesInvalidComment
            | SmilesUnterminatedBlockComment
            | SmilesInvalidElement
            | SmilesInvalidToken => Lexical,

            SmilesUnbalancedOpenParen
            | SmilesUnbalancedCloseParen
            | SmilesEmptyBranch
            | SmilesEmptyGroup
            | SmilesNonfinalGroup
            | SmilesLeadingBond
            | SmilesTrailingBond
            | SmilesConsecutiveBonds
            | SmilesLeadingRing
            | SmilesUnbalancedRingIndex
            | SmilesInvalidRingIndex
            | SmilesMismatchedRingBondDirs
            | SmilesMismatchedRingBondOrders
            | SmilesLeadingDot
            | SmilesTrailingDot
            | SmilesConsecutiveDots
            | SmilesDotBeforeRing
            | SmilesEmptyBracket
            | SmilesUnbalancedOpenBracket
            | SmilesUnbalancedCloseBracket
            | SmilesStrayBracketField
            | SmilesDuplicateBracketField
            | SmilesMissingClassIndex
            | SmilesMissingChiralityIndex
            | SmilesChiralityOutOfRange
            | SmilesBracketHwithHcount
            | SmilesInvalidBracket
            | SmilesClassOutOfRange
            | SmilesHcountOutOfRange
            | SmilesChargeOutOfRange
            | SmilesIsotopeOutOfRange
            | SmilesIsotopeUncatalogued
            | SmilesPreferBareOrganicAtom
            | SmilesPreferImplicitH
            | SmilesPreferBracketFieldOrder
            | SmilesPreferSimpleChargeSign
            | SmilesPreferSimpleHcount
            | SmilesAvoidExplicitSingleBond
            | SmilesAvoidExplicitAromaticBond
            | SmilesAvoidUnnecessaryGroup
            | SmilesAvoidRedundantNestedParens
            | SmilesPreferBranchesBeforeRingBonds
            | SmilesPreferFirstRingOne
            | SmilesPreferConsecutiveRingNumbering
            | SmilesAvoidReusedRingIndices
            | SmilesPreferSingleDigitRingIndex
            | SmilesPreferSingleRingClosure
            | SmilesAvoidAdjacentRingClosures
            | SmilesPreferBondSymbolAtRingOpen
            | SmilesAvoidRingClosureAcrossDot
            | SmilesPreferAromaticForm
            | CtfileInvalidCountsLine
            | CtfileInvalidAtomLine
            | CtfileInvalidBondLine
            | CtfileInvalidPropertyLine
            | CtfileInvalidSgroupLine
            | CtfileInvalidHeader
            | CtfileInvalidSdfDataHeader
            | CtfileInvalidSdfDataValue
            | CtfileMissingDelimiter => Syntactic,

            GraphConversionUnknown
            | GraphTopologySelfLoopRing
            | GraphTopologyParallelEdges
            | GraphTopologyUnknown
            | GraphValenceOutOfElementRange
            | GraphValenceHcountOutOfElementRange
            | GraphValenceChargeOutOfElementRange
            | GraphValenceHcountMismatch
            | GraphValenceNoMatch
            | GraphValenceAmbiguousMatch
            | GraphValenceNoKnownValenceStates
            | GraphValenceUnknownBondOrder
            | GraphValenceMissingBracketH
            | GraphValenceUnknown
            | GraphAromaticityAromaticAtomNotInRing
            | GraphAromaticityAromaticBondNotInRing
            | GraphAromaticityNoMatchingAromaticAtomConfig
            | GraphAromaticityInvalidAromaticAtom
            | GraphAromaticityInvalidAromaticBondAtom
            | GraphAromaticityAromaticBondOrderMismatch
            | GraphAromaticityKekuleInconsistent
            | GraphAromaticityHuckelFail
            | GraphAromaticityUnknown
            | GraphAromaticityWarningAvoidMixedAromaticity
            | GraphAromaticityWarningAvoidInconsistentAromaticity
            | GraphAromaticityWarningHuckelInconsistent
            | GraphAromaticityWarningUnknown
            | GraphStereoDoubleConflict
            | GraphStereoDoubleInsufficient
            | GraphStereoUnknown
            | GraphStereoWarningAvoidUnnecessaryStereoDescriptor
            | GraphStereoWarningUnsupportedCentralChiralityElement
            | GraphStereoWarningChiralitySubstituentMismatch
            | GraphStereoWarningNonChiralAnnotated
            | GraphStereoWarningUnknown
            | Unknown => Semantic,
        }
    }

    pub fn message(&self) -> &'static str {
        self.get_message().unwrap_or("Unknown diagnostic")
    }
}

#[derive(Debug, Default, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Diagnostic {
    pub kind: DiagnosticKind,
    pub category: Category,
    pub severity: Severity,
    pub span: Option<Span>,
    pub details: Option<String>,
}

impl Diagnostic {
    pub fn new(
        kind: DiagnosticKind,
        severity: Severity,
        span: Option<Span>,
        details: Option<String>,
    ) -> Self {
        Self {
            kind,
            category: kind.category(),
            severity,
            span,
            details,
        }
    }

    pub fn from_kind(kind: DiagnosticKind) -> Self {
        Self {
            kind,
            category: kind.category(),
            ..Default::default()
        }
    }

    pub fn severity(&self) -> Severity {
        self.severity
    }

    pub fn message(&self) -> &'static str {
        self.kind.message()
    }

    pub fn details(&self) -> Option<&String> {
        self.details.as_ref()
    }
}

impl Display for Diagnostic {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let span_str = match self.span {
            Some(span) => format!(" {}", span),
            None => "".to_string(),
        };

        write!(
            f,
            "{} [{}:{}]{}",
            self.message(),
            self.severity(),
            self.category,
            span_str,
        )
    }
}

#[derive(Debug, Default, Clone, PartialEq, Eq, Serialize, Deserialize)]
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
            .any(|d| d.severity() == Severity::Error)
    }

    pub fn errors(&self) -> impl Iterator<Item = &Diagnostic> {
        self.diagnostics
            .iter()
            .filter(|d| d.severity() == Severity::Error)
    }
    pub fn warnings(&self) -> impl Iterator<Item = &Diagnostic> {
        self.diagnostics
            .iter()
            .filter(|d| d.severity() == Severity::Warning)
    }
}

impl IntoIterator for DiagnosticList {
    type Item = Diagnostic;
    type IntoIter = std::vec::IntoIter<Diagnostic>;
    fn into_iter(self) -> Self::IntoIter {
        self.diagnostics.into_iter()
    }
}

impl Display for DiagnosticList {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for diagnostic in &self.diagnostics {
            writeln!(f, "- {}", diagnostic)?;
        }
        Ok(())
    }
}
