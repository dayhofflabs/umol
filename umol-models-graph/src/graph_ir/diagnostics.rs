//! Diagnostics for GraphIR.

use thiserror::Error;

#[derive(Debug, Clone, PartialEq, Error)]
pub enum ResolutionError {
    #[error("Self-loop ring")]
    SelfLoopRing,
    #[error("Parallel edges")]
    ParallelEdges,
    #[error("Out of element range")]
    OutOfElementRange,
    #[error("H count out of element range")]
    HcountOutOfElementRange,
    #[error("Charge out of element range")]
    ChargeOutOfElementRange,
    #[error("H count mismatch")]
    HcountMismatch,
    #[error("No match")]
    NoMatch,
    #[error("Ambiguous match")]
    AmbiguousMatch,
    #[error("No known valence states")]
    NoKnownValenceStates,
    #[error("Valence unknown bond order")]
    ValenceUnknownBondOrder,
    #[error("Missing bracket H")]
    MissingBracketH,
    #[error("Aromatic atom not in ring")]
    AtomNotInRing,
    #[error("Aromatic bond not in ring")]
    BondNotInRing,
    #[error("No matching aromatic atom config")]
    NoMatchingAtomConfig,
    #[error("Invalid aromatic atom")]
    InvalidAtom,
    #[error("Invalid aromatic bond atom")]
    InvalidBondAtom,
    #[error("Bond order mismatch")]
    BondOrderMismatch,
    #[error("Kekule inconsistent")]
    KekuleInconsistent,
    #[error("Huckel failed")]
    HuckelFailed,
    #[error("Double conflict")]
    DoubleConflict,
    #[error("Double insufficient")]
    DoubleInsufficient,
    #[error("Unsupported central element")]
    UnsupportedCentralElement,
    #[error("Chirality substituent mismatch")]
    SubstituentMismatch,
    #[error("Not chiral annotated")]
    NotChiralAnnotated,
}

pub enum ResolutionWarning {
    AvoidMixedAromaticity,
    AvoidInconsistentAromaticity,
    HuckelInconsistent,
    AvoidUnnecessaryStereoDescriptor,
    UnsupportedCentralChiralityElement,
    ChiralitySubstituentMismatch,
    NonChiralAnnotated,
}
