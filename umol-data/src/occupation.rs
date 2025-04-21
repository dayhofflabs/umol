//! Atomic occupations

use regex::Regex;
use std::fmt::{self, Display};
use std::str::FromStr;
use umol::error::DataError;
use umol::{Error, Result};

/// Atomic occupation of s, p, d, f orbitals.
/// Instances should typically be constructed via Configuration methods.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Occupation {
    /// s orbital occupation
    s: u8,
    /// p orbital occupation
    p: u8,
    /// d orbital occupation
    d: u8,
    /// f orbital occupation
    f: u8,
}

impl Occupation {
    /// Create a new occupation from spdf counts
    pub fn new(s: u8, p: u8, d: u8, f: u8) -> Self {
        Self { s, p, d, f }
    }

    /// Total number of electrons
    pub fn electron_count(&self) -> u8 {
        self.s + self.p + self.d + self.f
    }

    /// Number of s electrons
    pub fn s(&self) -> u8 {
        self.s
    }

    /// Number of p electrons
    pub fn p(&self) -> u8 {
        self.p
    }

    /// Number of d electrons
    pub fn d(&self) -> u8 {
        self.d
    }

    /// Number of f electrons
    pub fn f(&self) -> u8 {
        self.f
    }
}

impl Display for Occupation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut result = String::new();

        if self.s > 0 {
            result.push_str(&format!("s{}", self.s));
        }

        if self.p > 0 {
            result.push_str(&format!("p{}", self.p));
        }

        if self.d > 0 {
            result.push_str(&format!("d{}", self.d));
        }

        if self.f > 0 {
            result.push_str(&format!("f{}", self.f));
        }

        // Handle empty case (all zeros)
        if result.is_empty() {
            result = "".to_string();
        }

        write!(f, "{}", result)
    }
}

impl TryFrom<&str> for Occupation {
    type Error = Error;

    fn try_from(s: &str) -> Result<Self> {
        if s.is_empty() {
            return Err(DataError::InvalidOccupation(s.to_string()).into());
        }

        // Check if occupation string is valid
        let valid_occ_pattern = Regex::new(r"^([spdf](\d+))+$").unwrap();
        if !valid_occ_pattern.is_match(s) {
            return Err(DataError::InvalidOccupation(s.to_string()).into());
        }

        // Occupation regex: s<num>p<num>d<num>f<num>, in any order
        let occ_block_pattern = Regex::new(r"([spdf])(\d+)").unwrap();

        // Validate the string only contains valid orbital fragments
        if !occ_block_pattern.is_match(s) {
            return Err(DataError::InvalidOccupation(s.to_string()).into());
        }

        // Initialize occupation values
        let mut s_occ = 0;
        let mut p_occ = 0;
        let mut d_occ = 0;
        let mut f_occ = 0;

        // Process each capture
        for cap in occ_block_pattern.captures_iter(s) {
            let orbital_type = &cap[1];
            let count: u8 = cap[2].parse().unwrap_or(0);

            match orbital_type {
                "s" => s_occ = count,
                "p" => p_occ = count,
                "d" => d_occ = count,
                "f" => f_occ = count,
                _ => unreachable!(),
            }
        }

        Ok(Self {
            s: s_occ,
            p: p_occ,
            d: d_occ,
            f: f_occ,
        })
    }
}

impl FromStr for Occupation {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        Self::try_from(s)
    }
}

// TODO: Implement Serialize, Deserialize for Occupation

/// Shorthand macro for occupations
/// Allows using occupation strings directly without quotes
#[macro_export]
macro_rules! occ {
    ($occ:ident) => {
        stringify!($occ).parse::<Occupation>().unwrap()
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;

    #[rstest]
    #[case("s1p1d1f1", Occupation::new(1, 1, 1, 1))]
    #[case("s0p0d0f0", Occupation::new(0, 0, 0, 0))]
    #[case("s1", Occupation::new(1, 0, 0, 0))]
    #[case("p1", Occupation::new(0, 1, 0, 0))]
    #[case("d1", Occupation::new(0, 0, 1, 0))]
    #[case("f1", Occupation::new(0, 0, 0, 1))]
    #[case("s1d1", Occupation::new(1, 0, 1, 0))]
    #[case("p2d2", Occupation::new(0, 2, 2, 0))]
    #[case("d2p2", Occupation::new(0, 2, 2, 0))]
    #[case("s10p6", Occupation::new(10, 6, 0, 0))]
    fn test_from_str(#[case] s: &str, #[case] expected: Occupation) {
        assert_eq!(Occupation::from_str(s).unwrap(), expected);
    }

    #[rstest]
    #[case("")]
    #[case("s1p1d1f1x")]
    fn test_from_str_err(#[case] s: &str) {
        assert!(Occupation::from_str(s).is_err());
    }

    #[rstest]
    #[case(Occupation::new(1, 0, 0, 0), "s1")]
    #[case(Occupation::new(0, 1, 0, 0), "p1")]
    #[case(Occupation::new(0, 0, 1, 0), "d1")]
    #[case(Occupation::new(0, 0, 0, 1), "f1")]
    #[case(Occupation::new(1, 1, 1, 1), "s1p1d1f1")]
    #[case(Occupation::new(0, 0, 0, 0), "")]
    fn test_display(#[case] occupation: Occupation, #[case] expected: &str) {
        assert_eq!(format!("{}", occupation), expected);
    }

    #[test]
    fn test_macro() {
        assert_eq!(occ!(s1p1d1f1), Occupation::new(1, 1, 1, 1));
        assert_eq!(occ!(s0p0d0f0), Occupation::new(0, 0, 0, 0));
        assert_eq!(occ!(s1), Occupation::new(1, 0, 0, 0));
    }
}
