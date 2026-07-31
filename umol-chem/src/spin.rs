//! Spin multiplicity and spin state data

use std::fmt;
use std::num::NonZeroU8;
use std::str::FromStr;

use serde::de::{Deserializer, Error as SerdeError};
use serde::ser::Serializer;
use serde::{Deserialize, Serialize};

use crate::error::SpinStateError;

/// Exact unpaired-electron count and spin multiplicity.
///
/// This pair preserves structurally complete values without imposing the
/// physical invariants enforced by [`SpinState`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct UnpairedElectrons {
    pub count: i64,
    pub multiplicity: i64,
}

/// Positive spin multiplicity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(transparent)]
pub struct SpinMultiplicity(NonZeroU8);

impl SpinMultiplicity {
    pub const SINGLET: Self = Self(NonZeroU8::new(1).unwrap());
    pub const DOUBLET: Self = Self(NonZeroU8::new(2).unwrap());
    pub const TRIPLET: Self = Self(NonZeroU8::new(3).unwrap());
    pub const QUARTET: Self = Self(NonZeroU8::new(4).unwrap());
    pub const QUINTET: Self = Self(NonZeroU8::new(5).unwrap());
    pub const SEXTET: Self = Self(NonZeroU8::new(6).unwrap());
    pub const SEPTET: Self = Self(NonZeroU8::new(7).unwrap());
    pub const OCTET: Self = Self(NonZeroU8::new(8).unwrap());
    pub const NONET: Self = Self(NonZeroU8::new(9).unwrap());
    pub const DECET: Self = Self(NonZeroU8::new(10).unwrap());

    /// Construct a positive spin multiplicity.
    pub const fn new(value: u8) -> Option<Self> {
        match NonZeroU8::new(value) {
            Some(value) => Some(Self(value)),
            None => None,
        }
    }

    /// Return the conventional name for multiplicities one through ten.
    pub const fn name(self) -> Option<&'static str> {
        match self.0.get() {
            1 => Some("singlet"),
            2 => Some("doublet"),
            3 => Some("triplet"),
            4 => Some("quartet"),
            5 => Some("quintet"),
            6 => Some("sextet"),
            7 => Some("septet"),
            8 => Some("octet"),
            9 => Some("nonet"),
            10 => Some("decet"),
            _ => None,
        }
    }

    /// Construct from a conventional multiplicity name.
    pub fn from_name(name: &str) -> Option<Self> {
        if name.eq_ignore_ascii_case("singlet") {
            Some(Self::SINGLET)
        } else if name.eq_ignore_ascii_case("doublet") {
            Some(Self::DOUBLET)
        } else if name.eq_ignore_ascii_case("triplet") {
            Some(Self::TRIPLET)
        } else if name.eq_ignore_ascii_case("quartet") {
            Some(Self::QUARTET)
        } else if name.eq_ignore_ascii_case("quintet") {
            Some(Self::QUINTET)
        } else if name.eq_ignore_ascii_case("sextet") {
            Some(Self::SEXTET)
        } else if name.eq_ignore_ascii_case("septet") {
            Some(Self::SEPTET)
        } else if name.eq_ignore_ascii_case("octet") {
            Some(Self::OCTET)
        } else if name.eq_ignore_ascii_case("nonet") {
            Some(Self::NONET)
        } else if name.eq_ignore_ascii_case("decet") {
            Some(Self::DECET)
        } else {
            None
        }
    }
}

impl From<SpinMultiplicity> for u8 {
    fn from(multiplicity: SpinMultiplicity) -> Self {
        multiplicity.0.get()
    }
}

impl fmt::Display for SpinMultiplicity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

impl Serialize for SpinMultiplicity {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_u8(self.0.get())
    }
}

impl<'de> Deserialize<'de> for SpinMultiplicity {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let value = u8::deserialize(deserializer)?;
        Self::new(value).ok_or_else(|| SerdeError::custom("spin multiplicity must be nonzero"))
    }
}

/// Shorthand macro for spin-state literals parsed via `SpinState::from_str`.
///
/// Syntax: `#u<u>#s<s>` (e.g. `#u2#s3`).
#[macro_export]
macro_rules! spin {
    ($s:expr) => {{
        $crate::spin::SpinState::from_str($s).expect("invalid spin state")
    }};
}

/// Validated unpaired-electron count and multiplicity pair.
///
/// Invariant: `m <= u + 1` and `m` has the same parity as `u+1`, where
/// `m = u8::from(multiplicity)` and `u = unpaired_electrons`.
///
/// String format (canonical): `"#u<u>#s<s>"`, e.g. `"#u2#s3"`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct SpinState {
    unpaired_electrons: u8,
    multiplicity: SpinMultiplicity,
}

impl SpinState {
    /// Create a physically valid spin state.
    pub fn new(
        unpaired_electrons: u8,
        multiplicity: SpinMultiplicity,
    ) -> Result<Self, SpinStateError> {
        let max_multiplicity = u16::from(unpaired_electrons) + 1;
        let multiplicity_value = u16::from(u8::from(multiplicity));
        if multiplicity_value <= max_multiplicity && multiplicity_value % 2 == max_multiplicity % 2
        {
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

    /// Closed-shell singlet: 0 unpaired electrons, singlet multiplicity.
    pub fn closed_shell() -> Self {
        Self {
            unpaired_electrons: 0,
            multiplicity: SpinMultiplicity::SINGLET,
        }
    }

    pub fn unpaired_electrons(&self) -> u8 {
        self.unpaired_electrons
    }

    pub fn multiplicity(&self) -> SpinMultiplicity {
        self.multiplicity
    }
}

impl TryFrom<UnpairedElectrons> for SpinState {
    type Error = SpinStateError;

    fn try_from(unpaired_electrons: UnpairedElectrons) -> Result<Self, Self::Error> {
        let count = u8::try_from(unpaired_electrons.count).map_err(|_| {
            SpinStateError::UnpairedElectronsOutOfRange {
                count: unpaired_electrons.count,
            }
        })?;
        let multiplicity_value = u8::try_from(unpaired_electrons.multiplicity).map_err(|_| {
            SpinStateError::MultiplicityOutOfRange {
                multiplicity: unpaired_electrons.multiplicity,
            }
        })?;
        let multiplicity = SpinMultiplicity::new(multiplicity_value).ok_or(
            SpinStateError::MultiplicityOutOfRange {
                multiplicity: unpaired_electrons.multiplicity,
            },
        )?;
        Self::new(count, multiplicity)
    }
}

impl From<SpinState> for UnpairedElectrons {
    fn from(spin_state: SpinState) -> Self {
        Self {
            count: i64::from(spin_state.unpaired_electrons),
            multiplicity: i64::from(u8::from(spin_state.multiplicity)),
        }
    }
}

impl fmt::Display for SpinState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "#u{}#s{}",
            self.unpaired_electrons,
            u8::from(self.multiplicity)
        )
    }
}

/// Parses a DSL ground spin literal from one or both of `#u` and `#s` tags, in any order,
/// separated by optional whitespace. Omitting the decimal after a tag implies 1.
///
/// - `#u` alone: multiplicity = maximum for given unpaired electrons (Hund's rule).
/// - `#s` alone: unpaired electrons = `m - 1` (minimum for that multiplicity).
///
/// Each tag must not appear more than once.
/// - Both: validated as a pair.
impl FromStr for SpinState {
    type Err = SpinStateError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut rest = s.trim();
        let mut unpaired_electrons: Option<i64> = None;
        let mut multiplicity: Option<i64> = None;

        while !rest.is_empty() {
            if let Some(digits_str) = rest.strip_prefix("#u") {
                if unpaired_electrons.is_some() {
                    return Err(SpinStateError::DuplicateTag {
                        tag: "#u".to_string(),
                    });
                }
                let digits_len = digits_str.bytes().take_while(u8::is_ascii_digit).count();
                let value = if digits_len == 0 {
                    1
                } else {
                    digits_str[..digits_len].parse::<i64>().map_err(|_| {
                        SpinStateError::UnpairedElectronsOutOfRange { count: i64::MAX }
                    })?
                };
                unpaired_electrons = Some(value);
                rest = digits_str[digits_len..].trim_start();
            } else if let Some(digits_str) = rest.strip_prefix("#s") {
                if multiplicity.is_some() {
                    return Err(SpinStateError::DuplicateTag {
                        tag: "#s".to_string(),
                    });
                }
                let digits_len = digits_str.bytes().take_while(u8::is_ascii_digit).count();
                let value = if digits_len == 0 {
                    1
                } else {
                    digits_str[..digits_len].parse::<i64>().map_err(|_| {
                        SpinStateError::MultiplicityOutOfRange {
                            multiplicity: i64::MAX,
                        }
                    })?
                };
                multiplicity = Some(value);
                rest = digits_str[digits_len..].trim_start();
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

        match (unpaired_electrons, multiplicity) {
            (Some(count), Some(multiplicity)) => Self::try_from(UnpairedElectrons {
                count,
                multiplicity,
            }),
            (Some(count), None) => {
                let count_value = u8::try_from(count)
                    .map_err(|_| SpinStateError::UnpairedElectronsOutOfRange { count })?;
                let multiplicity = i64::from(u16::from(count_value) + 1);
                Self::try_from(UnpairedElectrons {
                    count,
                    multiplicity,
                })
            }
            (None, Some(multiplicity)) => {
                let multiplicity_value = u8::try_from(multiplicity)
                    .ok()
                    .and_then(SpinMultiplicity::new)
                    .ok_or(SpinStateError::MultiplicityOutOfRange { multiplicity })?;
                Self::new(u8::from(multiplicity_value) - 1, multiplicity_value)
            }
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
    use std::cmp::Ordering;

    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    fn test_unpaired_electrons_fields() {
        let unpaired_electrons = UnpairedElectrons {
            count: 2,
            multiplicity: 3,
        };

        assert_eq!(unpaired_electrons.count, 2);
        assert_eq!(unpaired_electrons.multiplicity, 3);
    }

    #[rstest]
    fn test_unpaired_electrons_eq() {
        let left = UnpairedElectrons {
            count: -1,
            multiplicity: 0,
        };
        let right = UnpairedElectrons {
            count: -1,
            multiplicity: 0,
        };

        assert_eq!(left, right);
    }

    #[rstest]
    #[case::count(
        UnpairedElectrons { count: 1, multiplicity: 5 },
        UnpairedElectrons { count: 2, multiplicity: 1 },
        Ordering::Less,
    )]
    #[case::multiplicity(
        UnpairedElectrons { count: 2, multiplicity: 1 },
        UnpairedElectrons { count: 2, multiplicity: 3 },
        Ordering::Less,
    )]
    #[case::equal(
        UnpairedElectrons { count: 2, multiplicity: 3 },
        UnpairedElectrons { count: 2, multiplicity: 3 },
        Ordering::Equal,
    )]
    fn test_unpaired_electrons_cmp(
        #[case] left: UnpairedElectrons,
        #[case] right: UnpairedElectrons,
        #[case] expected: Ordering,
    ) {
        assert_eq!(left.cmp(&right), expected);
    }

    #[rstest]
    #[case::physical(
        UnpairedElectrons { count: 2, multiplicity: 3 },
        r#"{"count":2,"multiplicity":3}"#,
    )]
    #[case::unvalidated(
        UnpairedElectrons { count: -1, multiplicity: 0 },
        r#"{"count":-1,"multiplicity":0}"#,
    )]
    fn test_unpaired_electrons_serde(
        #[case] unpaired_electrons: UnpairedElectrons,
        #[case] expected: &str,
    ) {
        let serialized = serde_json::to_string(&unpaired_electrons).unwrap();
        assert_eq!(serialized, expected);
        assert_eq!(
            serde_json::from_str::<UnpairedElectrons>(&serialized).unwrap(),
            unpaired_electrons,
        );
    }

    #[rstest]
    fn test_spin_multiplicity_new() {
        assert_eq!(SpinMultiplicity::new(0), None);
        for value in 1..=u8::MAX {
            assert_eq!(u8::from(SpinMultiplicity::new(value).unwrap()), value);
        }
    }

    #[rstest]
    #[case::singlet("singlet", "singlet", SpinMultiplicity::SINGLET)]
    #[case::doublet("doublet", "doublet", SpinMultiplicity::DOUBLET)]
    #[case::triplet("triplet", "triplet", SpinMultiplicity::TRIPLET)]
    #[case::quartet("quartet", "quartet", SpinMultiplicity::QUARTET)]
    #[case::quintet("quintet", "quintet", SpinMultiplicity::QUINTET)]
    #[case::sextet("sextet", "sextet", SpinMultiplicity::SEXTET)]
    #[case::septet("septet", "septet", SpinMultiplicity::SEPTET)]
    #[case::octet("octet", "octet", SpinMultiplicity::OCTET)]
    #[case::nonet("nonet", "nonet", SpinMultiplicity::NONET)]
    #[case::decet("decet", "decet", SpinMultiplicity::DECET)]
    #[case::case_insensitive("TrIpLeT", "triplet", SpinMultiplicity::TRIPLET)]
    fn test_spin_multiplicity_name(
        #[case] input: &str,
        #[case] expected_name: &str,
        #[case] expected: SpinMultiplicity,
    ) {
        assert_eq!(SpinMultiplicity::from_name(input), Some(expected));
        assert_eq!(expected.name(), Some(expected_name));
    }

    #[rstest]
    fn test_spin_multiplicity_name_nonconventional() {
        let multiplicity = SpinMultiplicity::new(11).unwrap();

        assert_eq!(multiplicity.name(), None);
        assert_eq!(SpinMultiplicity::from_name("undecet"), None);
    }

    #[rstest]
    #[case::singlet(SpinMultiplicity::SINGLET, "1")]
    #[case::triplet(SpinMultiplicity::TRIPLET, "3")]
    #[case::maximum(SpinMultiplicity::new(255).unwrap(), "255")]
    fn test_spin_multiplicity_fmt(#[case] multiplicity: SpinMultiplicity, #[case] expected: &str) {
        assert_eq!(multiplicity.to_string(), expected);
    }

    #[rstest]
    #[case::singlet(SpinMultiplicity::SINGLET, "1")]
    #[case::triplet(SpinMultiplicity::TRIPLET, "3")]
    #[case::maximum(SpinMultiplicity::new(255).unwrap(), "255")]
    fn test_spin_multiplicity_serde(
        #[case] multiplicity: SpinMultiplicity,
        #[case] expected: &str,
    ) {
        let serialized = serde_json::to_string(&multiplicity).unwrap();
        assert_eq!(serialized, expected);
        assert_eq!(
            serde_json::from_str::<SpinMultiplicity>(&serialized).unwrap(),
            multiplicity,
        );
    }

    #[rstest]
    fn test_spin_multiplicity_deserialize_zero() {
        assert!(serde_json::from_str::<SpinMultiplicity>("0").is_err());
    }

    #[rstest]
    #[case::closed_shell(0, SpinMultiplicity::SINGLET)]
    #[case::doublet(1, SpinMultiplicity::DOUBLET)]
    #[case::open_shell_singlet(2, SpinMultiplicity::SINGLET)]
    #[case::triplet(2, SpinMultiplicity::TRIPLET)]
    #[case::lower_spin(3, SpinMultiplicity::DOUBLET)]
    #[case::maximum_count(255, SpinMultiplicity::new(254).unwrap())]
    #[case::maximum_multiplicity(254, SpinMultiplicity::new(255).unwrap())]
    fn test_spin_state_new(#[case] unpaired_electrons: u8, #[case] multiplicity: SpinMultiplicity) {
        let spin_state = SpinState::new(unpaired_electrons, multiplicity).unwrap();

        assert_eq!(spin_state.unpaired_electrons(), unpaired_electrons);
        assert_eq!(spin_state.multiplicity(), multiplicity);
    }

    #[rstest]
    #[case::multiplicity_too_large(0, SpinMultiplicity::TRIPLET)]
    #[case::wrong_parity(1, SpinMultiplicity::SINGLET)]
    #[case::interior_wrong_parity(2, SpinMultiplicity::DOUBLET)]
    #[case::upper_boundary(253, SpinMultiplicity::new(255).unwrap())]
    fn test_spin_state_new_error(
        #[case] unpaired_electrons: u8,
        #[case] multiplicity: SpinMultiplicity,
    ) {
        assert_eq!(
            SpinState::new(unpaired_electrons, multiplicity),
            Err(SpinStateError::Incompatible {
                unpaired_electrons,
                multiplicity,
            }),
        );
    }

    #[rstest]
    fn test_spin_state_closed_shell() {
        assert_eq!(
            SpinState::closed_shell(),
            SpinState::new(0, SpinMultiplicity::SINGLET).unwrap(),
        );
    }

    #[rstest]
    #[case::closed_shell(
        UnpairedElectrons { count: 0, multiplicity: 1 },
        SpinState::closed_shell(),
    )]
    #[case::triplet(
        UnpairedElectrons { count: 2, multiplicity: 3 },
        SpinState::new(2, SpinMultiplicity::TRIPLET).unwrap(),
    )]
    #[case::upper_boundary(
        UnpairedElectrons { count: 254, multiplicity: 255 },
        SpinState::new(254, SpinMultiplicity::new(255).unwrap()).unwrap(),
    )]
    fn test_spin_state_try_from_unpaired_electrons(
        #[case] unpaired_electrons: UnpairedElectrons,
        #[case] expected: SpinState,
    ) {
        assert_eq!(SpinState::try_from(unpaired_electrons).unwrap(), expected);
    }

    #[rstest]
    #[case::negative_count(
        UnpairedElectrons { count: -1, multiplicity: 1 },
        SpinStateError::UnpairedElectronsOutOfRange { count: -1 },
    )]
    #[case::large_count(
        UnpairedElectrons { count: 256, multiplicity: 1 },
        SpinStateError::UnpairedElectronsOutOfRange { count: 256 },
    )]
    #[case::zero_multiplicity(
        UnpairedElectrons { count: 0, multiplicity: 0 },
        SpinStateError::MultiplicityOutOfRange { multiplicity: 0 },
    )]
    #[case::negative_multiplicity(
        UnpairedElectrons { count: 0, multiplicity: -1 },
        SpinStateError::MultiplicityOutOfRange { multiplicity: -1 },
    )]
    #[case::large_multiplicity(
        UnpairedElectrons { count: 0, multiplicity: 256 },
        SpinStateError::MultiplicityOutOfRange { multiplicity: 256 },
    )]
    #[case::incompatible(
        UnpairedElectrons { count: 2, multiplicity: 2 },
        SpinStateError::Incompatible {
            unpaired_electrons: 2,
            multiplicity: SpinMultiplicity::DOUBLET,
        },
    )]
    fn test_spin_state_try_from_unpaired_electrons_error(
        #[case] unpaired_electrons: UnpairedElectrons,
        #[case] expected: SpinStateError,
    ) {
        assert_eq!(SpinState::try_from(unpaired_electrons), Err(expected));
    }

    #[rstest]
    #[case::closed_shell(SpinState::closed_shell())]
    #[case::triplet(SpinState::new(2, SpinMultiplicity::TRIPLET).unwrap())]
    #[case::upper_boundary(
        SpinState::new(254, SpinMultiplicity::new(255).unwrap()).unwrap(),
    )]
    fn test_unpaired_electrons_from_spin_state(#[case] spin_state: SpinState) {
        let unpaired_electrons = UnpairedElectrons::from(spin_state);

        assert_eq!(SpinState::try_from(unpaired_electrons).unwrap(), spin_state);
    }

    #[rstest]
    #[case::closed_shell(SpinState::closed_shell(), "#u0#s1")]
    #[case::open_shell_singlet(
        SpinState::new(2, SpinMultiplicity::SINGLET).unwrap(),
        "#u2#s1",
    )]
    #[case::maximum_multiplicity(
        SpinState::new(254, SpinMultiplicity::new(255).unwrap()).unwrap(),
        "#u254#s255",
    )]
    fn test_spin_state_fmt(#[case] spin_state: SpinState, #[case] expected: &str) {
        assert_eq!(spin_state.to_string(), expected);
    }

    #[rstest]
    #[case::both("#u0#s1", 0, SpinMultiplicity::SINGLET)]
    #[case::reverse_order("#s3#u2", 2, SpinMultiplicity::TRIPLET)]
    #[case::whitespace(" #u2 #s1 ", 2, SpinMultiplicity::SINGLET)]
    #[case::count_only("#u2", 2, SpinMultiplicity::TRIPLET)]
    #[case::count_tag_only("#u", 1, SpinMultiplicity::DOUBLET)]
    #[case::multiplicity_only("#s4", 3, SpinMultiplicity::QUARTET)]
    #[case::multiplicity_tag_only("#s", 0, SpinMultiplicity::SINGLET)]
    #[case::above_conventional_range("#u10#s11", 10, SpinMultiplicity::new(11).unwrap())]
    #[case::maximum_multiplicity("#s255", 254, SpinMultiplicity::new(255).unwrap())]
    fn test_spin_state_from_str(
        #[case] input: &str,
        #[case] expected_count: u8,
        #[case] expected_multiplicity: SpinMultiplicity,
    ) {
        let spin_state = input.parse::<SpinState>().unwrap();

        assert_eq!(spin_state.unpaired_electrons(), expected_count);
        assert_eq!(spin_state.multiplicity(), expected_multiplicity);
    }

    #[rstest]
    #[case::empty("", SpinStateError::Underdetermined)]
    #[case::word("singlet", SpinStateError::UnexpectedToken { token: 's' })]
    #[case::number("0", SpinStateError::UnexpectedToken { token: '0' })]
    #[case::unknown_tag("#x3", SpinStateError::InvalidTag { tag: "#x3".to_string() })]
    #[case::duplicate_count("#u1#u2", SpinStateError::DuplicateTag { tag: "#u".to_string() })]
    #[case::duplicate_multiplicity("#s1#s2", SpinStateError::DuplicateTag { tag: "#s".to_string() })]
    #[case::count_out_of_range(
        "#u256#s1",
        SpinStateError::UnpairedElectronsOutOfRange { count: 256 },
    )]
    #[case::multiplicity_zero(
        "#u0#s0",
        SpinStateError::MultiplicityOutOfRange { multiplicity: 0 },
    )]
    #[case::multiplicity_out_of_range(
        "#u0#s256",
        SpinStateError::MultiplicityOutOfRange { multiplicity: 256 },
    )]
    #[case::derived_multiplicity_out_of_range(
        "#u255",
        SpinStateError::MultiplicityOutOfRange { multiplicity: 256 },
    )]
    #[case::incompatible(
        "#s2#u2",
        SpinStateError::Incompatible {
            unpaired_electrons: 2,
            multiplicity: SpinMultiplicity::DOUBLET,
        },
    )]
    fn test_spin_state_from_str_error(#[case] input: &str, #[case] expected: SpinStateError) {
        assert_eq!(input.parse::<SpinState>(), Err(expected));
    }

    #[rstest]
    #[case::closed_shell(SpinState::closed_shell(), r##""#u0#s1""##)]
    #[case::triplet(
        SpinState::new(2, SpinMultiplicity::TRIPLET).unwrap(),
        r##""#u2#s3""##,
    )]
    fn test_spin_state_serde(#[case] spin_state: SpinState, #[case] expected: &str) {
        let serialized = serde_json::to_string(&spin_state).unwrap();
        assert_eq!(serialized, expected);
        assert_eq!(
            serde_json::from_str::<SpinState>(&serialized).unwrap(),
            spin_state,
        );
    }

    #[rstest]
    fn test_spin_macro() {
        assert_eq!(
            spin!("#u2#s3"),
            SpinState::new(2, SpinMultiplicity::TRIPLET).unwrap(),
        );
    }
}
