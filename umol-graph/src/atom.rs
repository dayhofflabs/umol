//! Shared atomic value types used across IR layers.

/// Isotope mass ground term.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IsotopeMass {
    Natural,
    MassNumber(u32),
}

impl IsotopeMass {
    pub fn mass_number(&self) -> Option<u32> {
        match self {
            IsotopeMass::Natural => None,
            IsotopeMass::MassNumber(mass) => Some(*mass),
        }
    }
}

/// Aromatic valence ground term.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum AromaticValence {
    NotAromatic,
    Valence(u8),
}

impl AromaticValence {
    pub fn valence(&self) -> u8 {
        match self {
            AromaticValence::NotAromatic => 0,
            AromaticValence::Valence(n) => *n,
        }
    }

    /// Atom is aromatic if it contributes valence (n >= 0) to an aromatic system.
    pub fn is_aromatic(&self) -> bool {
        matches!(self, AromaticValence::Valence(_))
    }

    /// Extra count toward total valence for atoms in aromatic systems
    /// valence = 1 -> increment = 1 (benzene C)
    pub fn valence_increment(&self) -> u8 {
        match self {
            AromaticValence::NotAromatic => 0,
            AromaticValence::Valence(0) => 0,
            AromaticValence::Valence(1) => 1,
            AromaticValence::Valence(2) => 0,
            AromaticValence::Valence(_) => 0,
        }
    }
}

/// Chirality
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
