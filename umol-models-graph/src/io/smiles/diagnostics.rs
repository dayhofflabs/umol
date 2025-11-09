//! Diagnostics types for SMILES lexing and parsing.

use serde::{Deserialize, Serialize};
use strum::{AsRefStr, EnumIter, EnumMessage};

#[derive(
    Debug,
    Default,
    Clone,
    Copy,
    PartialEq,
    Eq,
    Hash,
    EnumIter,
    EnumMessage,
    AsRefStr,
    Serialize,
    Deserialize,
)]
#[strum(serialize_all = "SCREAMING_SNAKE_CASE")]
pub enum Lexical {
    #[strum(message = "Invalid whitespace")]
    InvalidWhitespace,
    #[strum(message = "Invalid comment")]
    InvalidComment,
    #[strum(message = "Unterminated block comment")]
    UnterminatedBlockComment,
    #[strum(message = "Invalid element")]
    InvalidElement,
    #[strum(message = "Invalid token")]
    InvalidToken,
    #[default]
    #[strum(message = "Unknown lexical error")]
    Unknown,
}

#[derive(
    Debug,
    Default,
    Clone,
    Copy,
    PartialEq,
    Eq,
    Hash,
    AsRefStr,
    EnumIter,
    EnumMessage,
    Serialize,
    Deserialize,
)]
#[strum(serialize_all = "SCREAMING_SNAKE_CASE")]
pub enum Syntactic {
    #[strum(message = "Unbalanced open parenthesis")]
    UnbalancedOpenParen,
    #[strum(message = "Unbalanced close parenthesis")]
    UnbalancedCloseParen,
    #[strum(message = "Empty branch")]
    EmptyBranch,
    #[strum(message = "Empty group")]
    EmptyGroup,
    #[strum(message = "Nonfinal group")]
    NonfinalGroup,
    #[strum(message = "Leading bond")]
    LeadingBond,
    #[strum(message = "Trailing bond")]
    TrailingBond,
    #[strum(message = "Consecutive bonds")]
    ConsecutiveBonds,
    LeadingRing,
    #[strum(message = "Unbalanced ring index")]
    UnbalancedRingIndex,
    #[strum(message = "Invalid ring index")]
    InvalidRingIndex,
    #[strum(message = "Mismatched ring bond directions")]
    MismatchedRingBondDirs,
    #[strum(message = "Mismatched ring bond orders")]
    MismatchedRingBondOrders,
    LeadingDot,
    #[strum(message = "Trailing dot")]
    TrailingDot,
    #[strum(message = "Consecutive dots")]
    ConsecutiveDots,
    #[strum(message = "Dot before ring")]
    DotBeforeRing,
    #[strum(message = "Empty bracket")]
    EmptyBracket,
    #[strum(message = "Unbalanced open bracket")]
    UnbalancedOpenBracket,
    #[strum(message = "Unbalanced close bracket")]
    UnbalancedCloseBracket,
    #[strum(message = "Stray bracket field")]
    StrayBracketField,
    #[strum(message = "Duplicate bracket field")]
    DuplicateBracketField,
    #[strum(message = "Missing class index")]
    MissingClassIndex,
    #[strum(message = "Missing chirality index")]
    MissingChiralityIndex,
    #[strum(message = "Bracket H with Hcount")]
    BracketHwithHcount,
    #[strum(message = "Invalid bracket")]
    InvalidBracket,
    #[default]
    #[strum(message = "Unknown syntactic error")]
    Unknown,
}

#[derive(
    Debug,
    Default,
    Clone,
    Copy,
    PartialEq,
    Eq,
    Hash,
    AsRefStr,
    EnumIter,
    EnumMessage,
    Serialize,
    Deserialize,
)]
#[strum(serialize_all = "SCREAMING_SNAKE_CASE")]
pub enum Numeric {
    #[strum(message = "Overflow")]
    Overflow,
    #[strum(message = "Class out of range")]
    ClassOutOfRange,
    #[strum(message = "Hcount out of range")]
    HcountOutOfRange,
    #[strum(message = "Charge out of range")]
    ChargeOutOfRange,
    #[strum(message = "Isotope out of range")]
    IsotopeOutOfRange,
    #[strum(message = "Chirality out of range")]
    ChiralityOutOfRange,
    #[default]
    #[strum(message = "Unknown numeric error")]
    Unknown,
}

#[derive(
    Debug,
    Default,
    Clone,
    Copy,
    PartialEq,
    Eq,
    Hash,
    AsRefStr,
    EnumIter,
    EnumMessage,
    Serialize,
    Deserialize,
)]
#[strum(serialize_all = "SCREAMING_SNAKE_CASE")]
pub enum NumericWarning {
    #[strum(message = "Isotope uncatalogued")]
    IsotopeUncatalogued,
    #[default]
    #[strum(message = "Unknown numeric warning")]
    Unknown,
}

#[derive(
    Debug,
    Default,
    Clone,
    Copy,
    PartialEq,
    Eq,
    Hash,
    AsRefStr,
    EnumIter,
    EnumMessage,
    Serialize,
    Deserialize,
)]
#[strum(serialize_all = "SCREAMING_SNAKE_CASE")]
pub enum StyleWarning {
    #[strum(message = "Prefer bare organic atom")]
    PreferBareOrganicAtom,
    #[strum(message = "Prefer implicit H")]
    PreferImplicitH,
    #[strum(message = "Prefer bracket field order")]
    PreferBracketFieldOrder,
    #[strum(message = "Prefer simple charge sign")]
    PreferSimpleChargeSign,
    #[strum(message = "Prefer simple H count")]
    PreferSimpleHcount,
    #[strum(message = "Avoid explicit single bond")]
    AvoidExplicitSingleBond,
    #[strum(message = "Avoid explicit aromatic bond")]
    AvoidExplicitAromaticBond,
    #[strum(message = "Avoid unnecessary group")]
    AvoidUnnecessaryGroup,
    #[strum(message = "Avoid redundant nested parentheses")]
    AvoidRedundantNestedParens,
    #[strum(message = "Prefer branches before ring bonds")]
    PreferBranchesBeforeRingBonds,
    #[strum(message = "Prefer start ring numbering with one")]
    PreferFirstRingOne,
    #[strum(message = "Prefer consecutive ring numbering")]
    PreferConsecutiveRingNumbering,
    #[strum(message = "Avoid reused ring indices")]
    AvoidReusedRingIndices,
    #[strum(message = "Prefer single digit ring index")]
    PreferSingleDigitRingIndex,
    #[strum(message = "Prefer single ring closure")]
    PreferSingleRingClosure,
    #[strum(message = "Avoid adjacent ring closures")]
    AvoidAdjacentRingClosures,
    #[strum(message = "Prefer bond symbol at ring open")]
    PreferBondSymbolAtRingOpen,
    #[strum(message = "Avoid ring closure across dot")]
    AvoidRingClosureAcrossDot,
    #[strum(message = "Prefer aromatic form")]
    PreferAromaticForm,
    #[default]
    #[strum(message = "Unknown style warning")]
    Unknown,
}
