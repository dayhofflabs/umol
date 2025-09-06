//! Spin state (multiplet) data

use std::fmt::{self, Display};
use std::str::FromStr;

use umol::error::DataError;
use umol::{Error, Result};

/// Spin state descriptor
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SpinState {
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

/// Highest spin state
pub const HIGHEST_SPIN_STATE: SpinState = SpinState::Decet;
/// Maximum number of unpaired electrons
pub const MAX_UNPAIRED_ELECTRONS: u8 = 9;

impl SpinState {
    /// Create a spin state from the number of unpaired electrons
    pub fn from_unpaired_electrons(unpaired_electrons: u8) -> Option<Self> {
        match unpaired_electrons {
            0 => Some(SpinState::Singlet),
            1 => Some(SpinState::Doublet),
            2 => Some(SpinState::Triplet),
            3 => Some(SpinState::Quartet),
            4 => Some(SpinState::Quintet),
            5 => Some(SpinState::Sextet),
            6 => Some(SpinState::Septet),
            7 => Some(SpinState::Octet),
            8 => Some(SpinState::Nonet),
            9 => Some(SpinState::Decet),
            _ => None,
        }
    }

    /// Create spin state from the multiplet name
    pub fn from_multiplet_name(name: &str) -> Option<Self> {
        match name.to_lowercase().as_str() {
            "singlet" => Some(SpinState::Singlet),
            "doublet" => Some(SpinState::Doublet),
            "triplet" => Some(SpinState::Triplet),
            "quartet" => Some(SpinState::Quartet),
            "quintet" => Some(SpinState::Quintet),
            "sextet" => Some(SpinState::Sextet),
            "septet" => Some(SpinState::Septet),
            "octet" => Some(SpinState::Octet),
            "nonet" => Some(SpinState::Nonet),
            "decet" => Some(SpinState::Decet),
            _ => None,
        }
    }

    /// Create spin state from multiplicity
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
            SpinState::Singlet => "singlet",
            SpinState::Doublet => "doublet",
            SpinState::Triplet => "triplet",
            SpinState::Quartet => "quartet",
            SpinState::Quintet => "quintet",
            SpinState::Sextet => "sextet",
            SpinState::Septet => "septet",
            SpinState::Octet => "octet",
            SpinState::Nonet => "nonet",
            SpinState::Decet => "decet",
        }
    }

    /// Get multiplicity
    pub fn multiplicity(&self) -> u8 {
        self.unpaired_electrons() + 1
    }
}

impl Display for SpinState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.name())
    }
}

impl TryFrom<&str> for SpinState {
    type Error = Error;

    fn try_from(s: &str) -> Result<Self> {
        Self::from_multiplet_name(s.to_lowercase().as_str())
            .ok_or_else(|| DataError::InvalidSpinState(s.to_string()).into())
    }
}

impl FromStr for SpinState {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        Self::try_from(s)
    }
}

/// Shorthand macro for spin states
/// Allows using spin state names directly without quotes
#[macro_export]
macro_rules! spin {
    ($state:ident) => {
        SpinState::$state
    };
}

// TODO: Implement Serialize, Deserialize for SpinState

#[cfg(test)]
mod tests {
    use rstest::*;

    use super::*;

    #[rstest]
    #[case(0, SpinState::Singlet)]
    #[case(1, SpinState::Doublet)]
    #[case(2, SpinState::Triplet)]
    fn test_spin_state_from_unpaired_electrons(
        #[case] unpaired_electrons: u8,
        #[case] expected: SpinState,
    ) {
        assert_eq!(
            SpinState::from_unpaired_electrons(unpaired_electrons).unwrap(),
            expected
        );
    }

    #[test]
    fn test_spin_state_from_unpaired_electrons_error() {
        assert!(SpinState::from_unpaired_electrons(10).is_none());
    }

    #[rstest]
    #[case("singlet", SpinState::Singlet)]
    #[case("doublet", SpinState::Doublet)]
    #[case("triplet", SpinState::Triplet)]
    #[case("Singlet", SpinState::Singlet)]
    fn test_spin_state_from_multiplet_name(
        #[case] multiplet_name: &str,
        #[case] expected: SpinState,
    ) {
        assert_eq!(
            SpinState::from_multiplet_name(multiplet_name).unwrap(),
            expected
        );
    }

    #[test]
    fn test_spin_state_from_multiplet_name_error() {
        assert!(SpinState::from_multiplet_name("nosuchplet").is_none());
    }

    #[rstest]
    #[case(1, SpinState::Singlet)]
    #[case(2, SpinState::Doublet)]
    #[case(3, SpinState::Triplet)]
    fn test_spin_state_from_multiplicity(#[case] multiplicity: u8, #[case] expected: SpinState) {
        assert_eq!(
            SpinState::from_multiplicity(multiplicity).unwrap(),
            expected
        );
    }

    #[test]
    fn test_spin_state_from_multiplicity_error() {
        assert!(SpinState::from_multiplicity(11).is_none());
    }

    #[rstest]
    #[case(SpinState::Singlet, 0)]
    #[case(SpinState::Doublet, 1)]
    #[case(SpinState::Triplet, 2)]
    fn test_spin_state_unpaired_electrons(#[case] spin_state: SpinState, #[case] expected: u8) {
        assert_eq!(spin_state.unpaired_electrons(), expected);
    }

    #[rstest]
    #[case(SpinState::Singlet, "singlet")]
    #[case(SpinState::Doublet, "doublet")]
    #[case(SpinState::Triplet, "triplet")]
    fn test_spin_state_name(#[case] spin_state: SpinState, #[case] expected: &str) {
        assert_eq!(spin_state.name(), expected);
    }

    #[rstest]
    #[case(SpinState::Singlet, 1)]
    #[case(SpinState::Doublet, 2)]
    #[case(SpinState::Triplet, 3)]
    fn test_spin_state_multiplicity(#[case] spin_state: SpinState, #[case] expected: u8) {
        assert_eq!(spin_state.multiplicity(), expected);
    }

    #[rstest]
    #[case(SpinState::Singlet, "singlet")]
    #[case(SpinState::Doublet, "doublet")]
    #[case(SpinState::Triplet, "triplet")]
    fn test_spin_state_display(#[case] spin_state: SpinState, #[case] expected: &str) {
        assert_eq!(format!("{}", spin_state), expected);
    }

    #[test]
    fn test_spin_state_macro() {
        assert_eq!(spin!(Singlet), SpinState::Singlet);
        assert_eq!(spin!(Doublet), SpinState::Doublet);
        assert_eq!(spin!(Triplet), SpinState::Triplet);
    }
}
