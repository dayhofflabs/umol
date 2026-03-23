//! Spin multiplicity and spin state data

use std::fmt;
use std::str::FromStr;

use serde::de::{Deserializer, Error as SerdeError};
use serde::ser::Serializer;
use serde::{Deserialize, Serialize};
use strum::{Display, EnumString};
use thiserror::Error;

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

/// Shorthand macro for spin-state literals parsed via `SpinState::from_str`.
///
// Syntax: `s{^nxm}` (e.g. `s{^2x3}`)
#[macro_export]
macro_rules! spin {
    ($s:expr) => {{
        $crate::SpinState::from_str($s).expect("invalid spin state")
    }};
}

/// Error for invalid spin-state values and parsing.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Error)]
pub enum SpinStateError {
    #[error("expected spin literal format 's{{^nxm}}'")]
    InvalidFormat,

    #[error("expected digits after '^'")]
    MissingUnpairedDigits,

    #[error("invalid unpaired electron count")]
    InvalidUnpairedElectrons,

    #[error("expected digits after 'x'")]
    MissingMultiplicityDigits,

    #[error("invalid multiplicity value")]
    InvalidMultiplicity,

    #[error("unexpected token '{token}', expected '^' or 'x'")]
    UnexpectedToken { token: char },

    #[error("missing '^n' in spin literal")]
    MissingUnpairedField,

    #[error("missing 'xm' in spin literal")]
    MissingMultiplicityField,

    #[error(
        "unpaired electron count exceeds maximum: {unpaired_electrons} > {MAX_UNPAIRED_ELECTRONS}"
    )]
    UnpairedElectronsExceedMax { unpaired_electrons: u8 },

    #[error("spin state is underdetermined")]
    Underdetermined,

    #[error("incompatible spin state combination: {unpaired_electrons} unpaired electrons, {multiplicity}")]
    Incompatible {
        unpaired_electrons: u8,
        multiplicity: SpinMultiplicity,
    },
}

/// Validated (unpaired_electrons, multiplicity) pair.
///
/// Invariant: `m <= n+1` and `m` has the same parity as `n+1`, where
/// `m = multiplicity.multiplicity()` and `n = unpaired_electrons`.
///
/// String format (canonical): `"s{^nxm}"`, e.g. `"s{^2x1}"`, `"s{^0x1}"`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct SpinState {
    unpaired_electrons: u8,
    multiplicity: SpinMultiplicity,
}

impl SpinState {
    /// Check whether `(unpaired_electrons, multiplicity)` is physically valid.
    ///
    /// Valid multiplicities for `n` unpaired electrons are `n%2+1, n%2+3, ..., n+1`.
    pub fn are_compatible(unpaired_electrons: u8, multiplicity: SpinMultiplicity) -> bool {
        let m = multiplicity.multiplicity();
        let n = unpaired_electrons;
        m <= n + 1 && m % 2 == (n + 1) % 2
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
            "s{{^{}x{}}}",
            self.unpaired_electrons,
            self.multiplicity.multiplicity()
        )
    }
}

impl FromStr for SpinState {
    type Err = SpinStateError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let trimmed = s.trim();
        let body = trimmed
            .strip_prefix("s{")
            .and_then(|rest| rest.strip_suffix('}'))
            .ok_or(SpinStateError::InvalidFormat)?;
        parse_spin_literal(body)
    }
}

fn parse_spin_literal(body: &str) -> Result<SpinState, SpinStateError> {
    let mut rest = body.trim();
    let mut unpaired: Option<u8> = None;
    let mut multiplicity: Option<SpinMultiplicity> = None;

    while !rest.is_empty() {
        rest = rest.trim_start();
        let Some(head) = rest.chars().next() else {
            break;
        };
        match head {
            '^' => {
                let digits_len = rest[1..]
                    .chars()
                    .take_while(|c| c.is_ascii_digit())
                    .map(char::len_utf8)
                    .sum::<usize>();
                if digits_len == 0 {
                    return Err(SpinStateError::MissingUnpairedDigits);
                }
                let n: u8 = rest[1..1 + digits_len]
                    .parse()
                    .map_err(|_| SpinStateError::InvalidUnpairedElectrons)?;
                unpaired = Some(n);
                rest = &rest[1 + digits_len..];
            }
            'x' => {
                let digits_len = rest[1..]
                    .chars()
                    .take_while(|c| c.is_ascii_digit())
                    .map(char::len_utf8)
                    .sum::<usize>();
                if digits_len == 0 {
                    return Err(SpinStateError::MissingMultiplicityDigits);
                }
                let m_u8: u8 = rest[1..1 + digits_len]
                    .parse()
                    .map_err(|_| SpinStateError::InvalidMultiplicity)?;
                let m = SpinMultiplicity::from_multiplicity(m_u8)
                    .ok_or(SpinStateError::InvalidMultiplicity)?;
                multiplicity = Some(m);
                rest = &rest[1 + digits_len..];
            }
            _ => {
                return Err(SpinStateError::UnexpectedToken { token: head });
            }
        }
    }

    let n = unpaired.ok_or(SpinStateError::MissingUnpairedField)?;
    let m = multiplicity.ok_or(SpinStateError::MissingMultiplicityField)?;
    SpinState::try_new(n, m)
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
    use pretty_assertions::assert_eq;

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
        assert_eq!(
            spin!("s{^2x3}"),
            SpinState::new(2, SpinMultiplicity::Triplet)
        );
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
    #[case(0, SpinMultiplicity::Singlet, "s{^0x1}")]
    #[case(1, SpinMultiplicity::Doublet, "s{^1x2}")]
    #[case(2, SpinMultiplicity::Singlet, "s{^2x1}")]
    #[case(2, SpinMultiplicity::Triplet, "s{^2x3}")]
    #[case(3, SpinMultiplicity::Quartet, "s{^3x4}")]
    fn test_spin_state_to_string(
        #[case] n: u8,
        #[case] m: SpinMultiplicity,
        #[case] expected: &str,
    ) {
        assert_eq!(SpinState::new(n, m).to_string(), expected);
    }

    #[rstest]
    #[case("s{^0x1}", 0, SpinMultiplicity::Singlet)]
    #[case("s{^1x2}", 1, SpinMultiplicity::Doublet)]
    #[case("s{^2x1}", 2, SpinMultiplicity::Singlet)]
    #[case("s{^2x3}", 2, SpinMultiplicity::Triplet)]
    #[case("s{^3x4}", 3, SpinMultiplicity::Quartet)]
    #[case("s{ ^2 x3 }", 2, SpinMultiplicity::Triplet)]
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
    #[case("s{}")]
    #[case("s{^2}")]
    #[case("s{x3}")]
    #[case("s{^2x2}")] // invalid parity for n=2
    #[case("s{^x3}")]
    #[case("s{^2x}")]
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

    #[rustfmt::skip]
    #[rstest]
    #[case(vec![], spin!("s{^0x1}"))]
    #[case(vec![spin!("s{^1x2}")], spin!("s{^1x2}"))]
    #[case(vec![spin!("s{^1x2}"), spin!("s{^1x2}")], spin!("s{^2x3}"))]
    #[case(vec![spin!("s{^2x3}"), spin!("s{^2x3}"), spin!("s{^2x3}")], spin!("s{^6x7}"))]
    #[case(vec![spin!("s{^2x1}"), spin!("s{^1x2}")], spin!("s{^3x4}"))]
    fn test_high_spin(#[case] states: Vec<SpinState>, #[case] expected: SpinState) {
        let hs = SpinState::high_spin(&states).unwrap();
        assert_eq!(hs, expected);
    }

    #[rstest]
    #[case(spin!("s{^0x1}"), 0, true)]
    #[case(spin!("s{^1x2}"), 1, true)]
    #[case(spin!("s{^0x1}"), 1, false)]
    #[case(spin!("s{^2x3}"), 0, false)]
    fn test_is_compatible_with(
        #[case] state: SpinState,
        #[case] electrons: u8,
        #[case] expected: bool,
    ) {
        assert_eq!(state.is_compatible_with(electrons), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case(vec![], spin!("s{^0x1}"), true)]
    #[case(vec![], spin!("s{^1x2}"), false)]
    #[case(vec![spin!("s{^1x2}"), spin!("s{^1x2}")], spin!("s{^0x1}"), false)]
    #[case(vec![spin!("s{^1x2}"), spin!("s{^1x2}")], spin!("s{^2x1}"), true)]
    #[case(vec![spin!("s{^1x2}"), spin!("s{^1x2}")], spin!("s{^2x3}"), true)]
    #[case(vec![spin!("s{^1x2}"), spin!("s{^1x2}")], spin!("s{^4x5}"), false)]
    #[case(vec![spin!("s{^2x3}"), spin!("s{^2x3}"), spin!("s{^2x3}")], spin!("s{^6x1}"), true)]
    #[case(vec![spin!("s{^2x3}"), spin!("s{^2x3}"), spin!("s{^2x3}")], spin!("s{^6x3}"), true)]
    #[case(vec![spin!("s{^2x3}"), spin!("s{^2x3}"), spin!("s{^2x3}")], spin!("s{^6x5}"), true)]
    #[case(vec![spin!("s{^2x3}"), spin!("s{^2x3}"), spin!("s{^2x3}")], spin!("s{^6x7}"), true)]
    #[case(vec![spin!("s{^2x3}"), spin!("s{^2x3}"), spin!("s{^2x3}")], spin!("s{^1x2}"), false)]
    fn test_is_constructible_from(
        #[case] states: Vec<SpinState>,
        #[case] target: SpinState,
        #[case] expected: bool,
    ) {
        assert_eq!(target.is_constructible_from(&states), expected);
    }
}
