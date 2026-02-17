//! Common bond data types for graph-based molecular models.

/// Ordered pair of atom indices, `first <= second`
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct AtomPair {
    first: u32,
    second: u32,
}

impl AtomPair {
    /// Create a new AtomPair, normalizing to ensure first <= second.
    pub fn new(a: u32, b: u32) -> Self {
        if a <= b {
            Self {
                first: a,
                second: b,
            }
        } else {
            Self {
                first: b,
                second: a,
            }
        }
    }

    /// Get the first (smaller) atom index.
    pub fn first(&self) -> u32 {
        self.first
    }

    /// Get the second (larger) atom index.
    pub fn second(&self) -> u32 {
        self.second
    }

    /// Get both atom indices as a tuple (first, second).
    pub fn as_tuple(&self) -> (u32, u32) {
        (self.first, self.second)
    }

    /// Check if this bond contains the given atom index.
    pub fn contains(&self, index: u32) -> bool {
        self.first == index || self.second == index
    }

    /// Get the other atom index.
    pub fn other(&self, index: u32) -> Option<u32> {
        if self.first == index {
            Some(self.second)
        } else if self.second == index {
            Some(self.first)
        } else {
            None
        }
    }
}

/// Electron pair donation for dative/coordinate bonds.
/// Direction is defined from the perspective of the first (smaller-indexed) atom.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum BondDonation {
    Shared,    // Normal covalent bond (both atoms contribute)
    Donating,  // First atom donates electron pair to second
    Accepting, // First atom accepts electron pair from second
}

impl BondDonation {
    /// Flip the donation
    pub fn flip(self) -> Self {
        match self {
            Self::Shared => Self::Shared,
            Self::Donating => Self::Accepting,
            Self::Accepting => Self::Donating,
        }
    }
}

/// Non-covalent interaction type for weak bonds (H-bonds, halogen bonds, etc.)
/// These bonds do not contribute to valence calculations.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BondNoncovalent {
    Hydrogen,
}
