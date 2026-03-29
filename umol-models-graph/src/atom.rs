//! Shared atomic value types used across IR layers.

use std::fmt::{self, Display};
use std::str::FromStr;

/// Isotope mass
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

impl Display for IsotopeMass {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            IsotopeMass::Natural => write!(f, "="),
            IsotopeMass::MassNumber(mass) => write!(f, "{}", mass),
        }
    }
}

impl FromStr for IsotopeMass {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s == "=" {
            return Ok(IsotopeMass::Natural);
        }
        if let Some(rest) = s.strip_prefix('i') {
            return rest
                .parse::<u32>()
                .map(IsotopeMass::MassNumber)
                .map_err(|_| format!("invalid isotope mass: {}", s));
        }
        s.parse::<u32>()
            .map(IsotopeMass::MassNumber)
            .map_err(|_| format!("invalid isotope mass: {}", s))
    }
}

/// Implicit hydrogens (Normal - inferred from normal valences)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ImplicitHydrogens {
    Hydrogens(u8),
    Normal,
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
        if s == "!" {
            return Ok(AromaticValence::None);
        }
        if let Some(rest) = s.strip_prefix('a') {
            return rest
                .parse::<u8>()
                .map(AromaticValence::Valence)
                .map_err(|_| format!("invalid aromatic valence: {}", s));
        }
        s.parse::<u8>()
            .map(AromaticValence::Valence)
            .map_err(|_| format!("invalid aromatic valence: {}", s))
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

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use pretty_assertions::assert_eq;
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::natural("=", IsotopeMass::Natural)]
    #[case::mass_number("12", IsotopeMass::MassNumber(12))]
    #[case::prefixed_mass("i13", IsotopeMass::MassNumber(13))]
    fn test_isotope_mass_parse(#[case] input: &str, #[case] expected: IsotopeMass) {
        let parsed = IsotopeMass::from_str(input);
        assert!(parsed.is_ok(), "isotope parse failed for {}", input);
        assert_eq!(parsed.unwrap(), expected);
    }

    #[rstest]
    #[case::empty("")]
    #[case::bad_prefix("x12")]
    #[case::alpha("foo")]
    #[case::negative("-1")]
    fn test_parse_isotope_mass_invalid(#[case] input: &str) {
        let parsed = IsotopeMass::from_str(input);
        assert!(
            parsed.is_err(),
            "expected isotope parse to fail for {}",
            input
        );
    }

    #[rstest]
    #[case::natural(IsotopeMass::Natural, "=")]
    #[case::mass_number(IsotopeMass::MassNumber(12), "12")]
    fn test_isotope_mass_display(#[case] value: IsotopeMass, #[case] expected: &str) {
        assert_eq!(value.to_string(), expected);
    }

    #[rstest]
    #[case::none_symbol("!", AromaticValence::None)]
    #[case::bare_number("0", AromaticValence::Valence(0))]
    #[case::bare_number_one("1", AromaticValence::Valence(1))]
    #[case::prefixed_number("a2", AromaticValence::Valence(2))]
    fn test_aromatic_valence_parse(#[case] input: &str, #[case] expected: AromaticValence) {
        let parsed = AromaticValence::from_str(input);
        assert!(parsed.is_ok(), "aromatic parse failed for {}", input);
        assert_eq!(parsed.unwrap(), expected);
    }

    #[rstest]
    #[case::empty("")]
    #[case::bad_prefix("b1")]
    #[case::alpha("foo")]
    #[case::negative("-1")]
    fn test_parse_aromatic_valence_invalid(#[case] input: &str) {
        let parsed = AromaticValence::from_str(input);
        assert!(
            parsed.is_err(),
            "expected aromatic parse to fail for {}",
            input
        );
    }

    #[rstest]
    #[case::none(AromaticValence::None, "")]
    #[case::valence_zero(AromaticValence::Valence(0), "a0")]
    #[case::valence_one(AromaticValence::Valence(1), "a1")]
    fn test_aromatic_valence_display(#[case] value: AromaticValence, #[case] expected: &str) {
        assert_eq!(value.to_string(), expected);
    }
}
