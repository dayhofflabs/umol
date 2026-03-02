//! Spin multiplicity and spin state data

use std::fmt;
use std::str::FromStr;

use serde::de::{Deserializer, Error as SerdeError};
use serde::ser::Serializer;
use serde::{Deserialize, Serialize};
use strum::{Display, EnumString};
use umol::error::DataError;

/// Spin multiplicity descriptor (2S+1)
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Display, EnumString)]
#[strum(serialize_all = "snake_case", ascii_case_insensitive)]
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

/// Maximum number of unpaired electrons representable by a `SpinState`.
pub const MAX_UNPAIRED_ELECTRONS: u8 = 9;

impl SpinMultiplicity {
    /// Create spin multiplicity from the numeric multiplicity value (2S+1).
    pub fn from_multiplicity(multiplicity: u8) -> Option<Self> {
        match multiplicity {
            1 => Some(SpinMultiplicity::Singlet),
            2 => Some(SpinMultiplicity::Doublet),
            3 => Some(SpinMultiplicity::Triplet),
            4 => Some(SpinMultiplicity::Quartet),
            5 => Some(SpinMultiplicity::Quintet),
            6 => Some(SpinMultiplicity::Sextet),
            7 => Some(SpinMultiplicity::Septet),
            8 => Some(SpinMultiplicity::Octet),
            9 => Some(SpinMultiplicity::Nonet),
            10 => Some(SpinMultiplicity::Decet),
            _ => None,
        }
    }

    /// Get multiplicity value (2S+1).
    pub fn multiplicity(&self) -> u8 {
        *self as u8 + 1
    }
}

/// Shorthand macro for spin multiplicity.
/// Allows using multiplicity names directly without quotes.
#[macro_export]
macro_rules! mult {
    ($state:ident) => {
        SpinMultiplicity::$state
    };
}

/// Check whether a (unpaired_electrons, multiplicity) pair is physically valid.
///
/// Valid multiplicities for `n` unpaired electrons are `n%2+1, n%2+3, ..., n+1`.
pub fn is_valid_spin_state(unpaired_electrons: u8, multiplicity: SpinMultiplicity) -> bool {
    let m = multiplicity.multiplicity();
    let n = unpaired_electrons;
    m <= n + 1 && m % 2 == (n + 1) % 2
}

/// Validated (unpaired_electrons, multiplicity) pair.
///
/// Invariant: `m <= n+1` and `m` has the same parity as `n+1`, where
/// `m = multiplicity.multiplicity()` and `n = unpaired_electrons`.
///
/// String format: `"{n},{multiplicity_name}"`, e.g. `"2,singlet"`, `"0,singlet"`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct SpinState {
    unpaired_electrons: u8,
    multiplicity: SpinMultiplicity,
}

impl SpinState {
    /// Create a spin state, panicking on invalid input.
    pub fn new(unpaired_electrons: u8, multiplicity: SpinMultiplicity) -> Self {
        Self::try_new(unpaired_electrons, multiplicity).unwrap_or_else(|| {
            panic!(
                "invalid spin state: {} unpaired electrons, {} multiplicity",
                unpaired_electrons,
                multiplicity.multiplicity()
            )
        })
    }

    /// Create a spin state, returning `None` if the combination is invalid.
    pub fn try_new(unpaired_electrons: u8, multiplicity: SpinMultiplicity) -> Option<Self> {
        if is_valid_spin_state(unpaired_electrons, multiplicity) {
            Some(Self {
                unpaired_electrons,
                multiplicity,
            })
        } else {
            None
        }
    }

    /// Create a spin state assuming maximum multiplicity (Hund's rule: m = n+1).
    pub fn max_multiplicity(unpaired_electrons: u8) -> Option<Self> {
        let m = SpinMultiplicity::from_multiplicity(unpaired_electrons + 1)?;
        Some(Self {
            unpaired_electrons,
            multiplicity: m,
        })
    }

    /// Closed-shell singlet: 0 unpaired electrons, singlet multiplicity.
    pub fn closed_shell() -> Self {
        Self {
            unpaired_electrons: 0,
            multiplicity: SpinMultiplicity::Singlet,
        }
    }

    pub fn unpaired_electrons(&self) -> u8 {
        self.unpaired_electrons
    }

    pub fn multiplicity(&self) -> SpinMultiplicity {
        self.multiplicity
    }
}

impl fmt::Display for SpinState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{},{}", self.unpaired_electrons, self.multiplicity)
    }
}

impl FromStr for SpinState {
    type Err = DataError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (n_str, m_str) = s.split_once(',').ok_or_else(|| {
            DataError::InvalidSpinState(format!("expected 'n,multiplicity', got '{}'", s))
        })?;
        let n: u8 = n_str.trim().parse().map_err(|_| {
            DataError::InvalidSpinState(format!("invalid unpaired electrons: '{}'", n_str))
        })?;
        let m: SpinMultiplicity = m_str.trim().parse().map_err(|_| {
            DataError::InvalidSpinState(format!("invalid multiplicity: '{}'", m_str))
        })?;
        Self::try_new(n, m).ok_or_else(|| {
            DataError::InvalidSpinState(format!(
                "inconsistent spin state: {} unpaired electrons, {}",
                n, m
            ))
        })
    }
}

impl Serialize for SpinState {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(&self.to_string())
    }
}

impl<'de> Deserialize<'de> for SpinState {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let s = String::deserialize(deserializer)?;
        s.parse().map_err(SerdeError::custom)
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;

    use super::*;

    #[rstest]
    #[case(1, SpinMultiplicity::Singlet)]
    #[case(2, SpinMultiplicity::Doublet)]
    #[case(3, SpinMultiplicity::Triplet)]
    #[case(10, SpinMultiplicity::Decet)]
    fn test_spin_multiplicity_from_multiplicity(
        #[case] multiplicity: u8,
        #[case] expected: SpinMultiplicity,
    ) {
        assert_eq!(
            SpinMultiplicity::from_multiplicity(multiplicity).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case(0)]
    #[case(11)]
    fn test_spin_multiplicity_from_multiplicity_invalid(#[case] multiplicity: u8) {
        assert!(SpinMultiplicity::from_multiplicity(multiplicity).is_none());
    }

    #[rstest]
    #[case(SpinMultiplicity::Singlet, 1)]
    #[case(SpinMultiplicity::Doublet, 2)]
    #[case(SpinMultiplicity::Triplet, 3)]
    #[case(SpinMultiplicity::Decet, 10)]
    fn test_spin_multiplicity_multiplicity(#[case] spin: SpinMultiplicity, #[case] expected: u8) {
        assert_eq!(spin.multiplicity(), expected);
    }

    #[rstest]
    #[case(SpinMultiplicity::Singlet, "singlet")]
    #[case(SpinMultiplicity::Doublet, "doublet")]
    #[case(SpinMultiplicity::Triplet, "triplet")]
    fn test_spin_multiplicity_to_string(#[case] spin: SpinMultiplicity, #[case] expected: &str) {
        assert_eq!(spin.to_string(), expected);
    }

    #[rstest]
    #[case("singlet", SpinMultiplicity::Singlet)]
    #[case("doublet", SpinMultiplicity::Doublet)]
    #[case("triplet", SpinMultiplicity::Triplet)]
    #[case("Singlet", SpinMultiplicity::Singlet)]
    #[case("DOUBLET", SpinMultiplicity::Doublet)]
    fn test_spin_multiplicity_parse(#[case] input: &str, #[case] expected: SpinMultiplicity) {
        assert_eq!(input.parse::<SpinMultiplicity>().unwrap(), expected);
    }

    #[test]
    fn test_spin_multiplicity_parse_error() {
        assert!("nosuchplet".parse::<SpinMultiplicity>().is_err());
    }

    #[test]
    fn test_mult_macro() {
        assert_eq!(mult!(Singlet), SpinMultiplicity::Singlet);
        assert_eq!(mult!(Doublet), SpinMultiplicity::Doublet);
        assert_eq!(mult!(Triplet), SpinMultiplicity::Triplet);
    }

    #[rstest]
    #[case(0, SpinMultiplicity::Singlet, true)]
    #[case(1, SpinMultiplicity::Doublet, true)]
    #[case(2, SpinMultiplicity::Singlet, true)]
    #[case(2, SpinMultiplicity::Triplet, true)]
    #[case(3, SpinMultiplicity::Doublet, true)]
    #[case(3, SpinMultiplicity::Quartet, true)]
    #[case(4, SpinMultiplicity::Singlet, true)]
    #[case(4, SpinMultiplicity::Triplet, true)]
    #[case(4, SpinMultiplicity::Quintet, true)]
    #[case(0, SpinMultiplicity::Doublet, false)]
    #[case(0, SpinMultiplicity::Triplet, false)]
    #[case(1, SpinMultiplicity::Singlet, false)]
    #[case(1, SpinMultiplicity::Triplet, false)]
    #[case(2, SpinMultiplicity::Doublet, false)]
    #[case(2, SpinMultiplicity::Quartet, false)]
    #[case(3, SpinMultiplicity::Singlet, false)]
    #[case(3, SpinMultiplicity::Triplet, false)]
    fn test_is_valid_spin_state(
        #[case] n: u8,
        #[case] m: SpinMultiplicity,
        #[case] expected: bool,
    ) {
        assert_eq!(is_valid_spin_state(n, m), expected);
    }

    #[rstest]
    #[case(0, SpinMultiplicity::Singlet)]
    #[case(1, SpinMultiplicity::Doublet)]
    #[case(2, SpinMultiplicity::Singlet)]
    #[case(2, SpinMultiplicity::Triplet)]
    #[case(3, SpinMultiplicity::Doublet)]
    #[case(3, SpinMultiplicity::Quartet)]
    fn test_spin_state_new(#[case] n: u8, #[case] m: SpinMultiplicity) {
        let state = SpinState::new(n, m);
        assert_eq!(state.unpaired_electrons(), n);
        assert_eq!(state.multiplicity(), m);
    }

    #[rstest]
    #[case(0, SpinMultiplicity::Singlet)]
    #[case(1, SpinMultiplicity::Doublet)]
    #[case(2, SpinMultiplicity::Singlet)]
    fn test_spin_state_try_new(#[case] n: u8, #[case] m: SpinMultiplicity) {
        assert!(SpinState::try_new(n, m).is_some());
    }

    #[rstest]
    #[case(0, SpinMultiplicity::Triplet)]
    #[case(1, SpinMultiplicity::Singlet)]
    #[case(2, SpinMultiplicity::Doublet)]
    fn test_spin_state_try_new_invalid(#[case] n: u8, #[case] m: SpinMultiplicity) {
        assert!(SpinState::try_new(n, m).is_none());
    }

    #[test]
    #[should_panic]
    fn test_spin_state_new_panics() {
        SpinState::new(0, SpinMultiplicity::Triplet);
    }

    #[rstest]
    #[case(0, SpinMultiplicity::Singlet)]
    #[case(1, SpinMultiplicity::Doublet)]
    #[case(2, SpinMultiplicity::Triplet)]
    #[case(9, SpinMultiplicity::Decet)]
    fn test_spin_state_max_multiplicity(#[case] n: u8, #[case] expected_m: SpinMultiplicity) {
        let state = SpinState::max_multiplicity(n).unwrap();
        assert_eq!(state.unpaired_electrons(), n);
        assert_eq!(state.multiplicity(), expected_m);
    }

    #[test]
    fn test_spin_state_max_multiplicity_overflow() {
        assert!(SpinState::max_multiplicity(10).is_none());
    }

    #[test]
    fn test_spin_state_closed_shell() {
        let state = SpinState::closed_shell();
        assert_eq!(state.unpaired_electrons(), 0);
        assert_eq!(state.multiplicity(), SpinMultiplicity::Singlet);
    }

    #[rstest]
    #[case(0, SpinMultiplicity::Singlet, "0,singlet")]
    #[case(1, SpinMultiplicity::Doublet, "1,doublet")]
    #[case(2, SpinMultiplicity::Singlet, "2,singlet")]
    #[case(2, SpinMultiplicity::Triplet, "2,triplet")]
    #[case(3, SpinMultiplicity::Quartet, "3,quartet")]
    fn test_spin_state_to_string(
        #[case] n: u8,
        #[case] m: SpinMultiplicity,
        #[case] expected: &str,
    ) {
        assert_eq!(SpinState::new(n, m).to_string(), expected);
    }

    #[rstest]
    #[case("0,singlet", 0, SpinMultiplicity::Singlet)]
    #[case("1,doublet", 1, SpinMultiplicity::Doublet)]
    #[case("2,singlet", 2, SpinMultiplicity::Singlet)]
    #[case("2,triplet", 2, SpinMultiplicity::Triplet)]
    #[case("3,quartet", 3, SpinMultiplicity::Quartet)]
    #[case(" 2 , triplet ", 2, SpinMultiplicity::Triplet)]
    fn test_spin_state_parse(
        #[case] input: &str,
        #[case] expected_n: u8,
        #[case] expected_m: SpinMultiplicity,
    ) {
        let state: SpinState = input.parse().unwrap();
        assert_eq!(state.unpaired_electrons(), expected_n);
        assert_eq!(state.multiplicity(), expected_m);
    }

    #[rstest]
    #[case("")]
    #[case("singlet")]
    #[case("0")]
    #[case("0,nosuchplet")]
    #[case("x,singlet")]
    #[case("1,singlet")] // invalid combination
    fn test_spin_state_parse_error(#[case] input: &str) {
        assert!(input.parse::<SpinState>().is_err());
    }

    #[rstest]
    #[case(0, SpinMultiplicity::Singlet)]
    #[case(2, SpinMultiplicity::Triplet)]
    #[case(3, SpinMultiplicity::Quartet)]
    fn test_spin_state_roundtrip(#[case] n: u8, #[case] m: SpinMultiplicity) {
        let state = SpinState::new(n, m);
        let s = state.to_string();
        let parsed: SpinState = s.parse().unwrap();
        assert_eq!(state, parsed);
    }

    #[rstest]
    #[case(0, SpinMultiplicity::Singlet)]
    #[case(2, SpinMultiplicity::Triplet)]
    fn test_spin_state_serde_roundtrip(#[case] n: u8, #[case] m: SpinMultiplicity) {
        let state = SpinState::new(n, m);
        let json = serde_json::to_string(&state).unwrap();
        let parsed: SpinState = serde_json::from_str(&json).unwrap();
        assert_eq!(state, parsed);
    }
}
