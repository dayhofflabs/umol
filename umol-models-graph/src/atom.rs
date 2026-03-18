//! Shared atomic value types used across IR layers.

use std::fmt::{self, Display};
use std::str::FromStr;

use umol_data::SpinMultiplicity;

/// Implicit hydrogens (Normal - inferred from normal valences)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ImplicitHydrogens {
    Hydrogens(u8),
    Normal,
}

/// Unpaired electron configuration.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct UnpairedElectrons {
    pub count: u8,
    pub multiplicity: Option<SpinMultiplicity>,
}

impl UnpairedElectrons {
    pub fn new(count: u8, multiplicity: Option<SpinMultiplicity>) -> Self {
        Self {
            count,
            multiplicity,
        }
    }

    pub fn from_count(count: u8) -> Self {
        Self {
            count,
            multiplicity: None,
        }
    }
}

/// Chirality.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Chirality {
    Clockwise,
    CounterClockwise,
    Unspecified,
    Tetrahedral { arr: u32 },
    Allenal { arr: u32 },
    SquarePlanar { arr: u32 },
    TrigonalBipyramidal { arr: u32 },
    Octahedral { arr: u32 },
}

/// Aromatic valence of an atom: none (non-aromatic) or contributing valence (n >= 0)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum AromaticValence {
    None,
    Valence(u8),
}

impl AromaticValence {
    pub fn valence(&self) -> u8 {
        match self {
            AromaticValence::None => 0,
            AromaticValence::Valence(n) => *n,
        }
    }

    /// Atom is aromatic if it contributes valence (n >= 0) to an aromatic system.
    pub fn is_aromatic(&self) -> bool {
        matches!(self, AromaticValence::Valence(_))
    }
}

impl Display for AromaticValence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AromaticValence::None => Ok(()),
            AromaticValence::Valence(n) => write!(f, "a{}", n),
        }
    }
}

impl FromStr for AromaticValence {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some(rest) = s.strip_prefix('a') {
            let n: u8 = rest
                .parse()
                .map_err(|_| format!("invalid aromatic valence: {}", s))?;
            Ok(AromaticValence::Valence(n))
        } else {
            Err(format!("expected 'a' prefix: {}", s))
        }
    }
}
