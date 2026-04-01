//! Spin multiplicity and spin state data

use std::fmt;
use std::str::FromStr;

use serde::de::{Deserializer, Error as SerdeError};
use serde::ser::Serializer;
use serde::{Deserialize, Serialize};
use strum::{Display, EnumString};
use thiserror::Error;

/// Spin multiplicity descriptor (2S+1)
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Display, EnumString,
         Serialize, Deserialize)]
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

/// Shorthand macro for spin-state literals parsed via `SpinState::from_str`.
///
/// Syntax: `#u<u> #m<m>` (e.g. `#u2 #m3`).
#[macro_export]
macro_rules! spin {
    ($s:expr) => {{
        $crate::SpinState::from_str($s).expect("invalid spin state")
    }};
}

/// Error for invalid spin-state values and parsing.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum SpinStateError {
    #[error("unexpected token '{token}'")]
    UnexpectedToken { token: char },

    #[error("invalid tag '{tag}'")]
    InvalidTag { tag: String },

    #[error("duplicate tag '{tag}'")]
    DuplicateTag { tag: String },

    #[error("unpaired electrons {unpaired_electrons} out of range")]
    UnpairedElectronsOutOfRange { unpaired_electrons: u8 },

    #[error("multiplicity {multiplicity} out of range")]
    MultiplicityOutOfRange { multiplicity: u8 },

    #[error("spin state is underdetermined")]
    Underdetermined,

    #[error("{unpaired_electrons} unpaired electrons, {multiplicity} multiplicity incompatible")]
    Incompatible {
        unpaired_electrons: u8,
        multiplicity: SpinMultiplicity,
    },
}

/// Validated (unpaired_electrons, multiplicity) pair.
///
/// Invariant: `m <= u + 1` and `m` has the same parity as `u+1`, where
/// `m = multiplicity.multiplicity()` and `u = unpaired_electrons`.
///
/// String format (canonical): `"#u<u> #m<m>"`, e.g. `"#u2 #m3"`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct SpinState {
    unpaired_electrons: u8,
    multiplicity: SpinMultiplicity,
}

impl SpinState {
    /// Check whether `(unpaired_electrons, multiplicity)` is physically valid.
    ///
    /// Valid multiplicities for `u` unpaired electrons are `u%2+1, u%2+3, ..., u+1`.
    pub fn are_compatible(unpaired_electrons: u8, multiplicity: SpinMultiplicity) -> bool {
        let u = unpaired_electrons;
        let m = multiplicity.multiplicity();
        m <= u + 1 && m % 2 == (u + 1) % 2
    }

    /// Create a spin state, panicking on invalid input.
    pub fn new(unpaired_electrons: u8, multiplicity: SpinMultiplicity) -> Self {
        Self::try_new(unpaired_electrons, multiplicity)
            .unwrap_or_else(|e| panic!("invalid spin state: {}", e))
    }

    /// Create a spin state, returning a domain error if the combination is invalid.
    pub fn try_new(
        unpaired_electrons: u8,
        multiplicity: SpinMultiplicity,
    ) -> Result<Self, SpinStateError> {
        if Self::are_compatible(unpaired_electrons, multiplicity) {
            Ok(Self {
                unpaired_electrons,
                multiplicity,
            })
        } else {
            Err(SpinStateError::Incompatible {
                unpaired_electrons,
                multiplicity,
            })
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

    /// High-spin molecular state from a collection of atomic spin states.
    ///
    /// Unpaired electrons = sum of atomic unpaired electrons.
    /// Multiplicity = max coupled multiplicity (all spins parallel).
    /// Returns `None` if the result exceeds `MAX_UNPAIRED_ELECTRONS`.
    pub fn high_spin(states: &[SpinState]) -> Option<Self> {
        let unpaired: u32 = states.iter().map(|s| s.unpaired_electrons as u32).sum();
        Self::max_multiplicity(unpaired as u8)
    }

    /// Check if molecular spin state is compatible with electron count.
    pub fn is_compatible_with(&self, electrons: u8) -> bool {
        self.unpaired_electrons <= electrons && (electrons - self.unpaired_electrons) % 2 == 0
    }

    /// Check whether this molecular spin state is achievable by coupling
    /// the given atomic spin states.
    ///
    /// Uses sequential angular momentum coupling: for spins S1, S2, the
    /// total S ranges from |S1-S2| to S1+S2 in integer steps. The set of
    /// allowed S values is order-independent.
    pub fn is_constructible_from(&self, states: &[SpinState]) -> bool {
        let unpaired: u32 = states.iter().map(|s| s.unpaired_electrons as u32).sum();
        if self.unpaired_electrons as u32 != unpaired {
            return false;
        }

        let target_two_s = (self.multiplicity.multiplicity() - 1) as u32;

        if states.is_empty() {
            return target_two_s == 0;
        }

        let mut possible: Vec<u32> = vec![(states[0].multiplicity.multiplicity() - 1) as u32];

        for state in &states[1..] {
            let two_s = (state.multiplicity.multiplicity() - 1) as u32;
            let mut next = Vec::new();
            for &prev in &possible {
                let lo = prev.abs_diff(two_s);
                let hi = prev + two_s;
                let mut s = lo;
                while s <= hi {
                    next.push(s);
                    s += 2;
                }
            }
            next.sort_unstable();
            next.dedup();
            possible = next;
        }

        possible.contains(&target_two_s)
    }
}

impl fmt::Display for SpinState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "#u{} #s{}",
            self.unpaired_electrons,
            self.multiplicity.multiplicity()
        )
    }
}

/// Parses a DSL ground spin literal from one or both of `#u` and `#s` tags, in any order,
/// separated by optional whitespace. Omitting the decimal after a tag implies 1.
///
/// - `#u` alone: multiplicity = maximum for given unpaired electrons (Hund's rule).
/// - `#s` alone: unpaired electrons = `m - 1` (minimum for that multiplicity).
/// Each tag must not appear more than once.
/// - Both: validated as a pair.
impl FromStr for SpinState {
    type Err = SpinStateError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut rest = s.trim();
        let mut unpaired: Option<u8> = None;
        let mut multiplicity: Option<SpinMultiplicity> = None;

        while !rest.is_empty() {
            if let Some(digits_str) = rest.strip_prefix("#u") {
                if unpaired.is_some() {
                    return Err(SpinStateError::DuplicateTag {
                        tag: "#u".to_string(),
                    });
                }
                let mut value = 0u8;
                let mut empty = true;
                let mut digits_len = 0;
                for digit in digits_str.chars() {
                    if !digit.is_ascii_digit() {
                        break;
                    }
                    value = value * 10 + (digit.to_digit(10).unwrap() as u8);
                    digits_len += 1;
                    empty = false;
                }
                if empty {
                    value = 1;
                }
                unpaired = Some(value);
                rest = &digits_str[digits_len..].trim_start();
            } else if let Some(digits_str) = rest.strip_prefix("#s") {
                if multiplicity.is_some() {
                    return Err(SpinStateError::DuplicateTag {
                        tag: "#s".to_string(),
                    });
                }
                let mut value = 0u8;
                let mut empty = true;
                let mut digits_len = 0;
                for digit in digits_str.chars() {
                    if !digit.is_ascii_digit() {
                        break;
                    }
                    value = value * 10 + (digit.to_digit(10).unwrap() as u8);
                    digits_len += 1;
                    empty = false;
                }
                if empty {
                    value = 1;
                }
                let mult = SpinMultiplicity::from_multiplicity(value).ok_or(
                    SpinStateError::MultiplicityOutOfRange {
                        multiplicity: value,
                    },
                )?;
                multiplicity = Some(mult);
                rest = &digits_str[digits_len..].trim_start();
            } else if rest.starts_with("#") {
                return Err(SpinStateError::InvalidTag {
                    tag: rest.to_string(),
                });
            } else {
                return Err(SpinStateError::UnexpectedToken {
                    token: rest.chars().next().expect("non-empty"),
                });
            }
        }

        match (unpaired, multiplicity) {
            (Some(u), Some(m)) => SpinState::try_new(u, m),
            (Some(u), None) => {
                SpinState::max_multiplicity(u).ok_or(SpinStateError::UnpairedElectronsOutOfRange {
                    unpaired_electrons: u,
                })
            }
            (None, Some(m)) => SpinState::try_new(m.multiplicity() - 1, m),
            (None, None) => Err(SpinStateError::Underdetermined),
        }
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
    use pretty_assertions::assert_eq;
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
    fn test_spin_macro() {
        assert_eq!(spin!("#u2"), SpinState::new(2, SpinMultiplicity::Triplet));
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
    fn test_are_compatible(#[case] n: u8, #[case] m: SpinMultiplicity, #[case] expected: bool) {
        assert_eq!(SpinState::are_compatible(n, m), expected);
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
        assert!(SpinState::try_new(n, m).is_ok());
    }

    #[rstest]
    #[case(0, SpinMultiplicity::Triplet)]
    #[case(1, SpinMultiplicity::Singlet)]
    #[case(2, SpinMultiplicity::Doublet)]
    fn test_spin_state_try_new_invalid(#[case] n: u8, #[case] m: SpinMultiplicity) {
        assert!(SpinState::try_new(n, m).is_err());
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
    fn test_spin_state_max_multiplicity_out_of_range() {
        assert!(SpinState::max_multiplicity(10).is_none());
    }

    #[test]
    fn test_spin_state_closed_shell() {
        let state = SpinState::closed_shell();
        assert_eq!(state.unpaired_electrons(), 0);
        assert_eq!(state.multiplicity(), SpinMultiplicity::Singlet);
    }

    #[rstest]
    #[case(0, SpinMultiplicity::Singlet, "#u0 #s1")]
    #[case(1, SpinMultiplicity::Doublet, "#u1 #s2")]
    #[case(2, SpinMultiplicity::Singlet, "#u2 #s1")]
    #[case(2, SpinMultiplicity::Triplet, "#u2 #s3")]
    #[case(3, SpinMultiplicity::Quartet, "#u3 #s4")]
    fn test_spin_state_to_string(
        #[case] n: u8,
        #[case] m: SpinMultiplicity,
        #[case] expected: &str,
    ) {
        assert_eq!(SpinState::new(n, m).to_string(), expected);
    }

    #[rstest]
    #[case("#u0 #s1", 0, SpinMultiplicity::Singlet)]
    #[case("#u1 #s2", 1, SpinMultiplicity::Doublet)]
    #[case("#u2 #s1", 2, SpinMultiplicity::Singlet)]
    #[case("#u2 #s3", 2, SpinMultiplicity::Triplet)]
    #[case("#u3 #s4", 3, SpinMultiplicity::Quartet)]
    #[case("#u2 #s3 ", 2, SpinMultiplicity::Triplet)]
    #[case("#u0", 0, SpinMultiplicity::Singlet)]
    #[case("#s1", 0, SpinMultiplicity::Singlet)]
    #[case("#s", 0, SpinMultiplicity::Singlet)]
    #[case("#u1", 1, SpinMultiplicity::Doublet)]
    #[case("#u", 1, SpinMultiplicity::Doublet)]
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
    #[case("", SpinStateError::Underdetermined)]
    #[case("singlet", SpinStateError::UnexpectedToken { token: 's' })]
    #[case("0", SpinStateError::UnexpectedToken { token: '0' })]
    #[case("x3", SpinStateError::UnexpectedToken { token: 'x' })]
    #[case("#", SpinStateError::InvalidTag { tag: "#".to_string() })]
    #[case("#x3", SpinStateError::InvalidTag { tag: "#x3".to_string() })]
    #[case("#u20", SpinStateError::UnpairedElectronsOutOfRange { unpaired_electrons: 20 })]
    #[case("#s20", SpinStateError::MultiplicityOutOfRange { multiplicity: 20 })]
    #[case("#s2#u2", SpinStateError::Incompatible { unpaired_electrons: 2, multiplicity: SpinMultiplicity::Doublet })]
    fn test_spin_state_parse_error(#[case] input: &str, #[case] expected: SpinStateError) {
        let result = input.parse::<SpinState>();
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), expected);
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

    #[rustfmt::skip]
    #[rstest]
    #[case(vec![], spin!("#s"))]
    #[case(vec![spin!("#u")], spin!("#u"))]
    #[case(vec![spin!("#u"), spin!("#u")], spin!("#u2"))]
    #[case(vec![spin!("#u2"), spin!("#u2"), spin!("#u2")], spin!("#u6"))]
    #[case(vec![spin!("#u2#s1"), spin!("#u")], spin!("#u3"))]
    fn test_high_spin(#[case] states: Vec<SpinState>, #[case] expected: SpinState) {
        let hs = SpinState::high_spin(&states).unwrap();
        assert_eq!(hs, expected);
    }

    #[rstest]
    #[case(spin!("#s"), 0, true)]
    #[case(spin!("#u"), 1, true)]
    #[case(spin!("#s"), 1, false)]
    #[case(spin!("#u2"), 0, false)]
    fn test_is_compatible_with(
        #[case] state: SpinState,
        #[case] electrons: u8,
        #[case] expected: bool,
    ) {
        assert_eq!(state.is_compatible_with(electrons), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case(vec![], spin!("#u0"), true)]
    #[case(vec![], spin!("#s2"), false)]
    #[case(vec![spin!("#s2"), spin!("#s2")], spin!("#u0"), false)]
    #[case(vec![spin!("#s2"), spin!("#s2")], spin!("#u2"), true)]
    #[case(vec![spin!("#s2"), spin!("#s2")], spin!("#u2 #s3"), true)]
    #[case(vec![spin!("#s2"), spin!("#s2")], spin!("#u4 #s5"), false)]
    #[case(vec![spin!("#u2"), spin!("#u2"), spin!("#u2")], spin!("#u6"), true)]
    #[case(vec![spin!("#u2"), spin!("#u2"), spin!("#u2")], spin!("#u6 #s3"), true)]
    #[case(vec![spin!("#u2"), spin!("#u2"), spin!("#u2")], spin!("#u6 #s5"), true)]
    #[case(vec![spin!("#u2"), spin!("#u2"), spin!("#u2")], spin!("#u6 #s7"), true)]
    #[case(vec![spin!("#u2"), spin!("#u2"), spin!("#u2")], spin!("#s2"), false)]
    fn test_is_constructible_from(
        #[case] states: Vec<SpinState>,
        #[case] target: SpinState,
        #[case] expected: bool,
    ) {
        assert_eq!(target.is_constructible_from(&states), expected);
    }
}
