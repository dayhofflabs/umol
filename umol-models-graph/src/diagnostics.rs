//! Diagnostics for atom/bond-based molecule models

use std::fmt;

use strum::{Display, EnumIter, EnumMessage, IntoEnumIterator};

use crate::span::Span;

#[derive(Copy, Clone, Debug, Display, PartialEq, EnumIter, EnumMessage)]
#[strum(serialize_all = "UPPERCASE")]
pub enum Severity {
    Error,
    Warning,
}

#[derive(Copy, Clone, Debug, Display, PartialEq, EnumIter)]
#[strum(serialize_all = "SCREAMING_SNAKE_CASE")]
pub enum Category {
    Lexical,
    Syntactic,
    Semantic,
}

#[derive(Copy, Clone, Debug, Display, PartialEq, EnumIter, EnumMessage)]
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
    #[strum(message = "Unexpected end of file")]
    CtfileUnexpectedEof,
    #[strum(message = "Incomplete input")]
    CtfileIncomplete,

    #[strum(message = "GraphIR conversion failed")]
    GraphConversionFailed,

    #[strum(message = "Self-loop ring")]
    GraphSelfLoopRing,
    #[strum(message = "Parallel edges")]
    GraphParallelEdges,
    #[strum(message = "Unknown topology error")]
    GraphTopologyError,

    #[strum(message = "Out of element range")]
    GraphOutOfElementRange,
    #[strum(message = "H count out of element range")]
    GraphHcountOutOfElementRange,
    #[strum(message = "Charge out of element range")]
    GraphChargeOutOfElementRange,
    #[strum(message = "H count mismatch")]
    GraphHcountMismatch,
    #[strum(message = "No match")]
    GraphNoMatch,
    #[strum(message = "Ambiguous match")]
    GraphAmbiguousMatch,
    #[strum(message = "No known valence states")]
    GraphNoKnownValenceStates,
    #[strum(message = "Valence unknown bond order")]
    GraphUnknownBondOrder,
    #[strum(message = "Missing bracket H")]
    GraphMissingBracketH,
    #[strum(message = "Unknown valence error")]
    GraphValenceError,

    #[strum(message = "Aromatic atom not in ring")]
    GraphAromaticAtomNotInRing,
    #[strum(message = "Aromatic bond not in ring")]
    GraphAromaticBondNotInRing,
    #[strum(message = "No matching aromatic atom config")]
    GraphNoMatchingAromaticAtomConfig,
    #[strum(message = "Invalid aromatic atom")]
    GraphInvalidAromaticAtom,
    #[strum(message = "Invalid aromatic bond atom")]
    GraphInvalidAromaticBondAtom,
    #[strum(message = "Aromatic bond order mismatch")]
    GraphAromaticBondOrderMismatch,
    #[strum(message = "Kekule inconsistent")]
    GraphKekuleInconsistent,
    #[strum(message = "Huckel failed")]
    GraphHuckelFailed,
    #[strum(message = "Unknown aromaticity error")]
    GraphAromaticityError,

    #[strum(message = "Avoid mixed aromaticity")]
    GraphAvoidMixedAromaticity,
    #[strum(message = "Avoid inconsistent aromaticity")]
    GraphAvoidInconsistentAromaticity,
    #[strum(message = "Huckel inconsistent")]
    GraphHuckelInconsistent,
    #[strum(message = "Unknown aromaticity warning")]
    GraphAromaticityWarning,

    #[strum(message = "Double conflict")]
    GraphStereoDoubleConflict,
    #[strum(message = "Double insufficient")]
    GraphStereoDoubleInsufficient,
    #[strum(message = "Unknown stereo error")]
    GraphStereoError,

    #[strum(message = "Avoid unnecessary stereo descriptor")]
    GraphAvoidUnnecessaryStereoDescriptor,
    #[strum(message = "Unsupported central chirality element")]
    GraphUnsupportedCentralChiralityElement,
    #[strum(message = "Chirality substituent mismatch")]
    GraphChiralitySubstituentMismatch,
    #[strum(message = "Non chiral annotated")]
    GraphNonChiralAnnotated,
    #[strum(message = "Unknown stereo warning")]
    GraphStereoWarning,
}

impl DiagnosticKind {
    pub fn all() -> impl Iterator<Item = DiagnosticKind> {
        DiagnosticKind::iter()
    }

    pub fn category(&self) -> Category {
        use DiagnosticKind::*;

        match self {
            SmilesInvalidWhitespace
            | SmilesInvalidComment
            | SmilesUnterminatedBlockComment
            | SmilesInvalidElement
            | SmilesInvalidToken => Category::Lexical,

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
            | CtfileMissingDelimiter
            | CtfileUnexpectedEof
            | CtfileIncomplete => Category::Syntactic,

            GraphConversionFailed
            | GraphSelfLoopRing
            | GraphParallelEdges
            | GraphTopologyError
            | GraphOutOfElementRange
            | GraphHcountOutOfElementRange
            | GraphChargeOutOfElementRange
            | GraphHcountMismatch
            | GraphNoMatch
            | GraphAmbiguousMatch
            | GraphNoKnownValenceStates
            | GraphUnknownBondOrder
            | GraphMissingBracketH
            | GraphValenceError
            | GraphAromaticAtomNotInRing
            | GraphAromaticBondNotInRing
            | GraphNoMatchingAromaticAtomConfig
            | GraphInvalidAromaticAtom
            | GraphInvalidAromaticBondAtom
            | GraphAromaticBondOrderMismatch
            | GraphKekuleInconsistent
            | GraphHuckelFailed
            | GraphAromaticityError
            | GraphAvoidMixedAromaticity
            | GraphAvoidInconsistentAromaticity
            | GraphHuckelInconsistent
            | GraphAromaticityWarning
            | GraphStereoDoubleConflict
            | GraphStereoDoubleInsufficient
            | GraphStereoError
            | GraphAvoidUnnecessaryStereoDescriptor
            | GraphUnsupportedCentralChiralityElement
            | GraphChiralitySubstituentMismatch
            | GraphNonChiralAnnotated
            | GraphStereoWarning => Category::Semantic,
        }
    }

    pub fn default_severity(&self) -> Severity {
        use DiagnosticKind::*;

        match self {
            SmilesInvalidWhitespace
            | SmilesInvalidComment
            | SmilesUnterminatedBlockComment
            | SmilesInvalidElement
            | SmilesInvalidToken
            | SmilesUnbalancedOpenParen
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
            | CtfileInvalidCountsLine
            | CtfileInvalidAtomLine
            | CtfileInvalidBondLine
            | CtfileInvalidPropertyLine
            | CtfileInvalidSgroupLine
            | CtfileInvalidHeader
            | CtfileInvalidSdfDataHeader
            | CtfileInvalidSdfDataValue
            | CtfileMissingDelimiter
            | CtfileUnexpectedEof
            | CtfileIncomplete
            | GraphConversionFailed
            | GraphSelfLoopRing
            | GraphParallelEdges
            | GraphTopologyError
            | GraphOutOfElementRange
            | GraphHcountOutOfElementRange
            | GraphChargeOutOfElementRange
            | GraphHcountMismatch
            | GraphNoMatch
            | GraphAmbiguousMatch
            | GraphNoKnownValenceStates
            | GraphUnknownBondOrder
            | GraphMissingBracketH
            | GraphValenceError
            | GraphAromaticAtomNotInRing
            | GraphAromaticBondNotInRing
            | GraphNoMatchingAromaticAtomConfig
            | GraphInvalidAromaticAtom
            | GraphInvalidAromaticBondAtom
            | GraphAromaticBondOrderMismatch
            | GraphKekuleInconsistent
            | GraphHuckelFailed
            | GraphAromaticityError
            | GraphStereoDoubleConflict
            | GraphStereoDoubleInsufficient
            | GraphStereoError => Severity::Error,

            SmilesPreferBareOrganicAtom
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
            | GraphAvoidMixedAromaticity
            | GraphAvoidInconsistentAromaticity
            | GraphHuckelInconsistent
            | GraphAromaticityWarning
            | GraphAvoidUnnecessaryStereoDescriptor
            | GraphUnsupportedCentralChiralityElement
            | GraphChiralitySubstituentMismatch
            | GraphNonChiralAnnotated
            | GraphStereoWarning => Severity::Warning,
        }
    }

    pub fn message(&self) -> &'static str {
        self.get_message().unwrap_or("Unknown diagnostic")
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Diagnostic {
    pub kind: DiagnosticKind,
    pub category: Category,
    pub severity: Severity,
    pub span: Option<Span>,
    pub details: Option<String>,
}

impl Diagnostic {
    pub fn from_kind(kind: DiagnosticKind) -> Self {
        Self {
            kind,
            category: kind.category(),
            severity: kind.default_severity(),
            span: None,
            details: None,
        }
    }
    pub fn from_kind_and_severity(
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

    pub fn message(&self) -> &'static str {
        self.kind.message()
    }

    pub fn details(&self) -> Option<&String> {
        self.details.as_ref()
    }
}

impl fmt::Display for Diagnostic {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let span_str = match self.span {
            Some(span) => format!(" {}", span),
            None => "".to_string(),
        };

        write!(
            f,
            "{} [{}:{}]{}",
            self.message(),
            self.severity,
            self.category,
            span_str,
        )
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct DiagnosticList {
    pub diagnostics: Vec<Diagnostic>,
}

impl DiagnosticList {
    pub fn new() -> Self {
        Self::default()
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

impl fmt::Display for DiagnosticList {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for diagnostic in &self.diagnostics {
            writeln!(f, "- {}", diagnostic)?;
        }
        Ok(())
    }
}
