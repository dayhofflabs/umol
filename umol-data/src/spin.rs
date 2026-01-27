//! Spin multiplicity data

use std::fmt::{self, Display};
use std::str::FromStr;

use umol::error::DataError;
use umol::{Error, Result};

/// Spin multiplicity descriptor (2S+1)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SpinMultiplicity {
    Singlet = 0,
    Doublet = 1,
    Triplet = 2,
    Quartet = 3,
    Quintet = 4,
    Sextet = 5,
    Septet = 6,
    Octet = 7,
    Nonet = 8,
    Decet = 9,
}

/// Highest spin multiplicity
pub const HIGHEST_SPIN_MULTIPLICITY: SpinMultiplicity = SpinMultiplicity::Decet;
/// Maximum number of unpaired electrons
pub const MAX_UNPAIRED_ELECTRONS: u8 = 9;

impl SpinMultiplicity {
    /// Create a spin multiplicity from the number of unpaired electrons
    pub fn from_unpaired_electrons(unpaired_electrons: u8) -> Option<Self> {
        match unpaired_electrons {
            0 => Some(SpinMultiplicity::Singlet),
            1 => Some(SpinMultiplicity::Doublet),
            2 => Some(SpinMultiplicity::Triplet),
            3 => Some(SpinMultiplicity::Quartet),
            4 => Some(SpinMultiplicity::Quintet),
            5 => Some(SpinMultiplicity::Sextet),
            6 => Some(SpinMultiplicity::Septet),
            7 => Some(SpinMultiplicity::Octet),
            8 => Some(SpinMultiplicity::Nonet),
            9 => Some(SpinMultiplicity::Decet),
            _ => None,
        }
    }

    /// Create spin multiplicity from the multiplet name
    pub fn from_multiplet_name(name: &str) -> Option<Self> {
        match name.to_lowercase().as_str() {
            "singlet" => Some(SpinMultiplicity::Singlet),
            "doublet" => Some(SpinMultiplicity::Doublet),
            "triplet" => Some(SpinMultiplicity::Triplet),
            "quartet" => Some(SpinMultiplicity::Quartet),
            "quintet" => Some(SpinMultiplicity::Quintet),
            "sextet" => Some(SpinMultiplicity::Sextet),
            "septet" => Some(SpinMultiplicity::Septet),
            "octet" => Some(SpinMultiplicity::Octet),
            "nonet" => Some(SpinMultiplicity::Nonet),
            "decet" => Some(SpinMultiplicity::Decet),
            _ => None,
        }
    }

    /// Create spin multiplicity from multiplicity value
    pub fn from_multiplicity(multiplicity: u8) -> Option<Self> {
        Self::from_unpaired_electrons(multiplicity - 1)
    }

    /// Get number of unpaired electrons
    pub fn unpaired_electrons(&self) -> u8 {
        *self as u8
    }

    /// Get multiplet name
    pub fn name(&self) -> &str {
        match self {
            SpinMultiplicity::Singlet => "singlet",
            SpinMultiplicity::Doublet => "doublet",
            SpinMultiplicity::Triplet => "triplet",
            SpinMultiplicity::Quartet => "quartet",
            SpinMultiplicity::Quintet => "quintet",
            SpinMultiplicity::Sextet => "sextet",
            SpinMultiplicity::Septet => "septet",
            SpinMultiplicity::Octet => "octet",
            SpinMultiplicity::Nonet => "nonet",
            SpinMultiplicity::Decet => "decet",
        }
    }

    /// Get multiplicity value (2S+1)
    pub fn multiplicity(&self) -> u8 {
        self.unpaired_electrons() + 1
    }
}

impl Display for SpinMultiplicity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.name())
    }
}

impl TryFrom<&str> for SpinMultiplicity {
    type Error = Error;

    fn try_from(s: &str) -> Result<Self> {
        Self::from_multiplet_name(s.to_lowercase().as_str())
            .ok_or_else(|| DataError::InvalidSpinMultiplicity(s.to_string()).into())
    }
}

impl FromStr for SpinMultiplicity {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        Self::try_from(s)
    }
}

/// Shorthand macro for spin multiplicity
/// Allows using multiplicity names directly without quotes
#[macro_export]
macro_rules! mult {
    ($state:ident) => {
        SpinMultiplicity::$state
    };
}

// TODO: Implement Serialize, Deserialize for SpinMultiplicity

#[cfg(test)]
mod tests {
    use rstest::*;

    use super::*;

    #[rstest]
    #[case(0, SpinMultiplicity::Singlet)]
    #[case(1, SpinMultiplicity::Doublet)]
    #[case(2, SpinMultiplicity::Triplet)]
    fn test_from_unpaired_electrons(#[case] unpaired_electrons: u8, #[case] expected: SpinMultiplicity) {
        assert_eq!(
            SpinMultiplicity::from_unpaired_electrons(unpaired_electrons).unwrap(),
            expected
        );
    }

    #[test]
    fn test_from_unpaired_electrons_error() {
        assert!(SpinMultiplicity::from_unpaired_electrons(10).is_none());
    }

    #[rstest]
    #[case("singlet", SpinMultiplicity::Singlet)]
    #[case("doublet", SpinMultiplicity::Doublet)]
    #[case("triplet", SpinMultiplicity::Triplet)]
    #[case("Singlet", SpinMultiplicity::Singlet)]
    fn test_from_multiplet_name(#[case] multiplet_name: &str, #[case] expected: SpinMultiplicity) {
        assert_eq!(
            SpinMultiplicity::from_multiplet_name(multiplet_name).unwrap(),
            expected
        );
    }

    #[test]
    fn test_from_multiplet_name_error() {
        assert!(SpinMultiplicity::from_multiplet_name("nosuchplet").is_none());
    }

    #[rstest]
    #[case(1, SpinMultiplicity::Singlet)]
    #[case(2, SpinMultiplicity::Doublet)]
    #[case(3, SpinMultiplicity::Triplet)]
    fn test_from_multiplicity(#[case] multiplicity: u8, #[case] expected: SpinMultiplicity) {
        assert_eq!(
            SpinMultiplicity::from_multiplicity(multiplicity).unwrap(),
            expected
        );
    }

    #[test]
    fn test_from_multiplicity_error() {
        assert!(SpinMultiplicity::from_multiplicity(11).is_none());
    }

    #[rstest]
    #[case(SpinMultiplicity::Singlet, 0)]
    #[case(SpinMultiplicity::Doublet, 1)]
    #[case(SpinMultiplicity::Triplet, 2)]
    fn test_unpaired_electrons(#[case] spin: SpinMultiplicity, #[case] expected: u8) {
        assert_eq!(spin.unpaired_electrons(), expected);
    }

    #[rstest]
    #[case(SpinMultiplicity::Singlet, "singlet")]
    #[case(SpinMultiplicity::Doublet, "doublet")]
    #[case(SpinMultiplicity::Triplet, "triplet")]
    fn test_name(#[case] spin: SpinMultiplicity, #[case] expected: &str) {
        assert_eq!(spin.name(), expected);
    }

    #[rstest]
    #[case(SpinMultiplicity::Singlet, 1)]
    #[case(SpinMultiplicity::Doublet, 2)]
    #[case(SpinMultiplicity::Triplet, 3)]
    fn test_multiplicity(#[case] spin: SpinMultiplicity, #[case] expected: u8) {
        assert_eq!(spin.multiplicity(), expected);
    }

    #[rstest]
    #[case(SpinMultiplicity::Singlet, "singlet")]
    #[case(SpinMultiplicity::Doublet, "doublet")]
    #[case(SpinMultiplicity::Triplet, "triplet")]
    fn test_display(#[case] spin: SpinMultiplicity, #[case] expected: &str) {
        assert_eq!(format!("{}", spin), expected);
    }

    #[test]
    fn test_mult_macro() {
        assert_eq!(mult!(Singlet), SpinMultiplicity::Singlet);
        assert_eq!(mult!(Doublet), SpinMultiplicity::Doublet);
        assert_eq!(mult!(Triplet), SpinMultiplicity::Triplet);
    }
}
