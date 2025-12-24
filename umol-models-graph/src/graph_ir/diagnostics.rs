//! Diagnostics for GraphIR.

use serde::{Deserialize, Serialize};
use strum::{Display, EnumIter, EnumMessage};

#[derive(
    Debug,
    Default,
    Clone,
    Copy,
    PartialEq,
    Eq,
    Hash,
    Display,
    EnumIter,
    EnumMessage,
    Serialize,
    Deserialize,
)]
#[strum(serialize_all = "SCREAMING_SNAKE_CASE")]
pub enum Conversion {
    #[default]
    #[strum(message = "GraphIR conversion failed")]
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
    Display,
    EnumIter,
    EnumMessage,
    Serialize,
    Deserialize,
)]
#[strum(serialize_all = "SCREAMING_SNAKE_CASE")]
pub enum Topology {
    #[strum(message = "Self-loop ring")]
    SelfLoopRing,
    #[strum(message = "Parallel edges")]
    ParallelEdges,
    #[default]
    #[strum(message = "Unknown topology error")]
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
    Display,
    EnumIter,
    EnumMessage,
    Serialize,
    Deserialize,
)]
#[strum(serialize_all = "SCREAMING_SNAKE_CASE")]
pub enum Valence {
    #[strum(message = "Out of element range")]
    OutOfElementRange,
    #[strum(message = "H count out of element range")]
    HcountOutOfElementRange,
    #[strum(message = "Charge out of element range")]
    ChargeOutOfElementRange,
    #[strum(message = "H count mismatch")]
    HcountMismatch,
    #[strum(message = "No match")]
    NoMatch,
    #[strum(message = "Ambiguous match")]
    AmbiguousMatch,
    #[strum(message = "No known valence states")]
    NoKnownValenceStates,
    #[strum(message = "Valence unknown bond order")]
    ValenceUnknownBondOrder,
    #[strum(message = "Missing bracket H")]
    MissingBracketH,
    #[default]
    #[strum(message = "Unknown valence error")]
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
    Display,
    EnumIter,
    EnumMessage,
    Serialize,
    Deserialize,
)]
#[strum(serialize_all = "SCREAMING_SNAKE_CASE")]
pub enum Aromaticity {
    #[strum(message = "Aromatic atom not in ring")]
    AromaticAtomNotInRing,
    #[strum(message = "Aromatic bond not in ring")]
    AromaticBondNotInRing,
    #[strum(message = "No matching aromatic atom config")]
    NoMatchingAromaticAtomConfig,
    #[strum(message = "Invalid aromatic atom")]
    InvalidAromaticAtom,
    #[strum(message = "Invalid aromatic bond atom")]
    InvalidAromaticBondAtom,
    #[strum(message = "Aromatic bond order mismatch")]
    AromaticBondOrderMismatch,
    #[strum(message = "Kekule inconsistent")]
    KekuleInconsistent,
    #[strum(message = "Huckel fail")]
    HuckelFail,
    #[default]
    #[strum(message = "Unknown aromaticity error")]
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
    Display,
    EnumIter,
    EnumMessage,
    Serialize,
    Deserialize,
)]
#[strum(serialize_all = "SCREAMING_SNAKE_CASE")]
pub enum AromaticityWarning {
    #[strum(message = "Avoid mixed aromaticity")]
    AvoidMixedAromaticity,
    #[strum(message = "Avoid inconsistent aromaticity")]
    AvoidInconsistentAromaticity,
    #[strum(message = "Huckel inconsistent")]
    HuckelInconsistent,
    #[default]
    #[strum(message = "Unknown aromaticity warning")]
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
    Display,
    EnumIter,
    EnumMessage,
    Serialize,
    Deserialize,
)]
#[strum(serialize_all = "SCREAMING_SNAKE_CASE")]
pub enum Stereo {
    #[strum(message = "Double conflict")]
    DoubleConflict,
    #[strum(message = "Double insufficient")]
    DoubleInsufficient,
    #[default]
    #[strum(message = "Unknown stereo error")]
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
    Display,
    EnumIter,
    EnumMessage,
    Serialize,
    Deserialize,
)]
#[strum(serialize_all = "SCREAMING_SNAKE_CASE")]
pub enum StereoWarning {
    #[strum(message = "Avoid unnecessary stereo descriptor")]
    AvoidUnnecessaryStereoDescriptor,
    #[strum(message = "Unsupported central chirality element")]
    UnsupportedCentralChiralityElement,
    #[strum(message = "Chirality substituent mismatch")]
    ChiralitySubstituentMismatch,
    #[strum(message = "Non chiral annotated")]
    NonChiralAnnotated,
    #[default]
    #[strum(message = "Unknown stereo warning")]
    Unknown,
}
