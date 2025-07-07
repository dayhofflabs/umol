//! Isotope definitions and data

use map_macro::hash_map;
use once_cell::sync::Lazy;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fmt::{self, Display};
use std::str::FromStr;
use umol::error::DataError;
use umol::{Error, Result};

use crate::Element;
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, PartialOrd)]
pub enum NamedIsotope {
    D,
    T,
}

/// Named isotope data:
///
/// 0. element
/// 1. isotope mass number
/// 2. isotope mass (in amu)
/// 3. isotope symbol
static NAMED_ISOTOPE_DATA: Lazy<HashMap<NamedIsotope, (Element, u32, &'static str)>> =
    Lazy::new(|| {
        hash_map! {
            NamedIsotope::D => (Element::H, 2, "D"),
            NamedIsotope::T => (Element::H, 3, "T"),
        }
    });

static SYMBOL_TO_NAMED_ISOTOPE: Lazy<HashMap<&'static str, NamedIsotope>> = Lazy::new(|| {
    NAMED_ISOTOPE_DATA
        .iter()
        .map(|(isotope, data)| (data.2, *isotope))
        .collect()
});

static NAMED_ISOTOPE_SYMBOLS: Lazy<HashSet<&'static str>> = Lazy::new(|| {
    NAMED_ISOTOPE_DATA
        .iter()
        .map(|(_, data)| data.2) // data.2 is the &'static str symbol (e.g., "D", "T")
        .collect()
});

impl NamedIsotope {
    // Get named isotope from symbol bytestring (allocation-free)
    pub fn from_symbol_bytes(symbol: &[u8]) -> Option<Self> {
        match symbol.len() {
            1 => {
                if !symbol[0].is_ascii_alphabetic() {
                    return None;
                }
                let upper_b = symbol[0].to_ascii_uppercase();
                let mut buf = [0u8; 1];
                buf[0] = upper_b;
                let key_str = core::str::from_utf8(&buf).unwrap();
                SYMBOL_TO_NAMED_ISOTOPE.get(key_str).copied()
            }
            _ => None,
        }
    }

    // Get named isotope from symbol string (allocation-free)
    pub fn from_symbol(symbol: &str) -> Option<Self> {
        Self::from_symbol_bytes(symbol.as_bytes())
    }

    /// Only hydrogen isotope are named
    pub fn element(&self) -> Element {
        NAMED_ISOTOPE_DATA.get(self).unwrap().0
    }

    /// Get the isotope mass number
    pub fn mass_number(&self) -> u32 {
        NAMED_ISOTOPE_DATA.get(self).unwrap().1
    }

    /// Get the isotope mass (in amu)
    pub fn mass(&self) -> f64 {
        ISOTOPE_DATA
            .get(&(self.element(), self.mass_number()))
            .unwrap()
            .0
    }

    /// Get the isotope half-life (in s)
    pub fn half_life(&self) -> Option<f64> {
        ISOTOPE_DATA
            .get(&(self.element(), self.mass_number()))
            .unwrap()
            .1
    }

    /// Get the isotope symbol
    pub fn symbol(&self) -> &'static str {
        NAMED_ISOTOPE_DATA.get(self).unwrap().2
    }

    /// Check if bytestring contains valid named isotope
    pub fn is_named_isotope_bytes(symbol: &[u8]) -> bool {
        match symbol.len() {
            1 => {
                if !symbol[0].is_ascii_alphabetic() {
                    return false;
                }
                let upper_b = symbol[0].to_ascii_uppercase();
                let mut key_buf = [0u8; 1];
                key_buf[0] = upper_b;
                // This unwrap is safe: a single ASCII char is valid UTF-8.
                let lookup_key_str = core::str::from_utf8(&key_buf).unwrap();
                NAMED_ISOTOPE_SYMBOLS.contains(lookup_key_str)
            }
            _ => false, // Named isotopes are single characters like D or T
        }
    }

    /// Check if string contains valid named isotope
    pub fn is_named_isotope(symbol: &str) -> bool {
        Self::is_named_isotope_bytes(symbol.as_bytes())
    }
}

impl TryFrom<&str> for NamedIsotope {
    type Error = Error;

    fn try_from(s: &str) -> Result<Self> {
        Self::from_symbol(s).ok_or_else(|| DataError::InvalidIsotope(s.to_string()).into())
    }
}

impl FromStr for NamedIsotope {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        Self::try_from(s)
    }
}

impl Display for NamedIsotope {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.symbol())
    }
}

impl From<NamedIsotope> for Element {
    fn from(isotope: NamedIsotope) -> Self {
        isotope.element()
    }
}

impl From<NamedIsotope> for Isotope {
    fn from(isotope: NamedIsotope) -> Self {
        Isotope {
            element: isotope.element(),
            mass_number: isotope.mass_number(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, PartialOrd)]
pub struct Isotope {
    element: Element,
    mass_number: u32,
}

/// Isotope data:
///
/// 0. isotope mass (in amu)
/// 1. half-life (in s)
static ISOTOPE_DATA: Lazy<HashMap<(Element, u32), (f64, Option<f64>)>> = Lazy::new(|| {
    hash_map! {
        // Hydrogen
        (Element::H, 1) => (1.00782503223, None), // Protium
        (Element::H, 2) => (2.01410177812, None), // Deuterium
        (Element::H, 3) => (3.0160492779, Some(3.8878e8)),  // Tritium (12.33 years)
        // Helium
        (Element::He, 3) => (3.0160293201, None),
        (Element::He, 4) => (4.00260325413, None),
        // Lithium
        (Element::Li, 6) => (6.0151228874, None),
        (Element::Li, 7) => (7.0160034366, None),
        // Beryllium
        (Element::Be, 9) => (9.012183065, None),
        // Boron
        (Element::B, 10) => (10.01293695, None),
        (Element::B, 11) => (11.00930536, None),
        // Carbon
        (Element::C, 12) => (12.0000000, None),
        (Element::C, 13) => (13.00335483507, None),
        // Nitrogen
        (Element::N, 14) => (14.00307400443, None),
        (Element::N, 15) => (15.00010889888, None),
        // Oxygen
        (Element::O, 16) => (15.99491461957, None),
        (Element::O, 17) => (16.99913175650, None),
        (Element::O, 18) => (17.99915961286, None),
        // Fluorine
        (Element::F, 19) => (18.99840316273, None),
        // Neon
        (Element::Ne, 20) => (19.9924401762, None),
        (Element::Ne, 21) => (20.993846685, None),
        (Element::Ne, 22) => (21.991385114, None),
    }
});

impl Isotope {
    /// Create a new isotope if it exists in ISOTOPE_DATA.  
    pub fn checked_new(element: Element, mass_number: u32) -> Option<Self> {
        if ISOTOPE_DATA.contains_key(&(element, mass_number)) {
            Some(Isotope {
                element,
                mass_number,
            })
        } else {
            None
        }
    }

    /// Create a new isotope from symbol bytestring
    /// The symbol must be in the format of "AZ": A: mass number, Z: element symbol
    /// Example: "1H", "2H", "3H", "4He", "13C", "226Ra"
    pub fn from_symbol_bytes(symbol: &[u8]) -> Option<Self> {
        let symbol_len = symbol.len();

        // Smallest is "1H" (len 2), largest could be e.g., "294Og" (len 5)
        if !(2..=5).contains(&symbol_len) {
            return None;
        }

        let mut idx: usize = 0;
        while idx < symbol_len && symbol[idx].is_ascii_digit() {
            idx += 1;
        }

        // Need >= 1 digit for mass number, >= 1 character for element symbol
        if idx == 0 || idx == symbol_len {
            return None;
        }

        let mass_number_bytes = &symbol[..idx];
        let element_symbol_bytes = &symbol[idx..];

        // Element symbol part must be 1 or 2 chars
        if !(1..=2).contains(&element_symbol_bytes.len()) {
            return None;
        }

        let mass_number_str = core::str::from_utf8(mass_number_bytes).ok()?;
        let mass_number = mass_number_str.parse::<u32>().ok()?;
        if mass_number == 0 {
            return None;
        }

        let element = Element::from_symbol_bytes(element_symbol_bytes)?;

        Self::checked_new(element, mass_number)
    }

    /// Create a new isotope from symbol string
    /// see `from_symbol_bytes` for format details
    pub fn from_symbol(symbol: &str) -> Option<Self> {
        Self::from_symbol_bytes(symbol.as_bytes())
    }

    /// Get isotope mass (in amu)
    pub fn mass(&self) -> f64 {
        ISOTOPE_DATA
            .get(&(self.element, self.mass_number))
            .unwrap()
            .0
    }

    /// Get isotope half-life (in s)
    pub fn half_life(&self) -> Option<f64> {
        ISOTOPE_DATA
            .get(&(self.element, self.mass_number))
            .unwrap()
            .1
    }

    /// Get isotope symbol (AZ notation)
    pub fn symbol(&self) -> String {
        format!("{}{}", self.mass_number, self.element.symbol())
    }
}

impl TryFrom<&str> for Isotope {
    type Error = Error;

    fn try_from(s: &str) -> Result<Self> {
        Self::from_symbol(s).ok_or_else(|| DataError::InvalidIsotope(s.to_string()).into())
    }
}

impl FromStr for Isotope {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        Self::try_from(s)
    }
}

impl Display for Isotope {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.symbol())
    }
}

impl From<Isotope> for Element {
    fn from(isotope: Isotope) -> Self {
        isotope.element
    }
}

#[macro_export]
macro_rules! iso {
    ($isotope:expr) => {
        Isotope::from_symbol($isotope).unwrap()
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use float_cmp::*;
    use rstest::*;
    use serde_json;

    #[test]
    fn test_named_isotope_from_symbol_bytes() {
        assert_eq!(NamedIsotope::from_symbol_bytes(b"D"), Some(NamedIsotope::D));
        assert_eq!(NamedIsotope::from_symbol_bytes(b"d"), Some(NamedIsotope::D));
        assert_eq!(NamedIsotope::from_symbol_bytes(b"T"), Some(NamedIsotope::T));
        assert_eq!(NamedIsotope::from_symbol_bytes(b"H"), None);
        assert_eq!(NamedIsotope::from_symbol_bytes(b"X"), None);
        assert_eq!(NamedIsotope::from_symbol_bytes(b""), None);
    }

    #[test]
    fn test_named_isotope_from_symbol() {
        assert_eq!(NamedIsotope::from_symbol("D"), Some(NamedIsotope::D));
        assert_eq!(NamedIsotope::from_symbol("d"), Some(NamedIsotope::D));
        assert_eq!(NamedIsotope::from_symbol("T"), Some(NamedIsotope::T));
        assert_eq!(NamedIsotope::from_symbol("H"), None);
        assert_eq!(NamedIsotope::from_symbol("X"), None);
        assert_eq!(NamedIsotope::from_symbol(""), None);
    }

    #[rstest]
    #[case(NamedIsotope::D, Element::H, 2, 2.01410177812, None, "D")]
    #[case(NamedIsotope::T, Element::H, 3, 3.0160492779, Some(3.8878e8), "T")]
    fn test_named_isotope_properties(
        #[case] named_isotope: NamedIsotope,
        #[case] expected_element: Element,
        #[case] expected_mass_number: u32,
        #[case] expected_mass: f64,
        #[case] expected_half_life: Option<f64>,
        #[case] expected_symbol: &str,
    ) {
        assert_eq!(named_isotope.element(), expected_element);
        assert_eq!(named_isotope.mass_number(), expected_mass_number);
        assert!(approx_eq!(
            f64,
            named_isotope.mass(),
            expected_mass,
            ulps = 4
        ));
        assert_eq!(named_isotope.half_life(), expected_half_life);
        assert_eq!(named_isotope.symbol(), expected_symbol);
    }

    #[test]
    fn test_named_isotope_is_named_isotope_bytes() {
        assert!(NamedIsotope::is_named_isotope_bytes(b"D"));
        assert!(NamedIsotope::is_named_isotope_bytes(b"d"));
        assert!(NamedIsotope::is_named_isotope_bytes(b"T"));
        assert!(NamedIsotope::is_named_isotope_bytes(b"t"));
        assert!(!NamedIsotope::is_named_isotope_bytes(b"H"));
        assert!(!NamedIsotope::is_named_isotope_bytes(b"X"));
        assert!(!NamedIsotope::is_named_isotope_bytes(b""));
    }

    #[test]
    fn test_named_isotope_is_named_isotope() {
        assert!(NamedIsotope::is_named_isotope("D"));
        assert!(NamedIsotope::is_named_isotope("d"));
        assert!(NamedIsotope::is_named_isotope("T"));
        assert!(NamedIsotope::is_named_isotope("t"));
        assert!(!NamedIsotope::is_named_isotope("H"));
        assert!(!NamedIsotope::is_named_isotope("X"));
        assert!(!NamedIsotope::is_named_isotope(""));
    }

    #[test]
    fn test_named_isotope_display() {
        assert_eq!(NamedIsotope::D.to_string(), "D");
        assert_eq!(NamedIsotope::T.to_string(), "T");
    }

    #[test]
    fn test_named_isotope_from_str() {
        assert_eq!("D".parse::<NamedIsotope>().unwrap(), NamedIsotope::D);
        assert_eq!("T".parse::<NamedIsotope>().unwrap(), NamedIsotope::T);
        assert!("X".parse::<NamedIsotope>().is_err());
    }

    #[test]
    fn test_named_isotope_try_from_str() {
        assert_eq!(NamedIsotope::try_from("D").unwrap(), NamedIsotope::D);
        assert_eq!(NamedIsotope::try_from("T").unwrap(), NamedIsotope::T);
        assert!(NamedIsotope::try_from("X").is_err());
    }

    #[test]
    fn test_named_isotope_to_element_conversion() {
        assert_eq!(Element::from(NamedIsotope::D), Element::H);
        assert_eq!(Element::from(NamedIsotope::T), Element::H);
    }

    #[test]
    fn test_named_isotope_to_isotope_conversion() {
        let isotope_d = Isotope::from(NamedIsotope::D);
        assert_eq!(isotope_d.element, Element::H);
        assert_eq!(isotope_d.mass_number, 2);

        let isotope_t = Isotope::from(NamedIsotope::T);
        assert_eq!(isotope_t.element, Element::H);
        assert_eq!(isotope_t.mass_number, 3);
    }

    #[test]
    fn test_named_isotope_serialization() {
        assert_eq!(serde_json::to_string(&NamedIsotope::D).unwrap(), r#""D""#);
        assert_eq!(serde_json::to_string(&NamedIsotope::T).unwrap(), r#""T""#);
        let isotopes = vec![NamedIsotope::D, NamedIsotope::T];
        assert_eq!(serde_json::to_string(&isotopes).unwrap(), r#"["D","T"]"#);
    }

    #[test]
    fn test_named_isotope_deserialization() {
        assert_eq!(
            serde_json::from_str::<NamedIsotope>(r#""D""#).unwrap(),
            NamedIsotope::D
        );
        assert_eq!(
            serde_json::from_str::<NamedIsotope>(r#""T""#).unwrap(),
            NamedIsotope::T
        );
        let isotopes: Vec<NamedIsotope> = serde_json::from_str(r#"["D","T"]"#).unwrap();
        assert_eq!(isotopes, vec![NamedIsotope::D, NamedIsotope::T]);
        assert!(serde_json::from_str::<NamedIsotope>(r#""X""#).is_err());
    }

    #[test]
    fn test_named_isotope_roundtrip() {
        let isotopes = vec![NamedIsotope::D, NamedIsotope::T];
        let serialized = serde_json::to_string(&isotopes).unwrap();
        let deserialized: Vec<NamedIsotope> = serde_json::from_str(&serialized).unwrap();
        assert_eq!(isotopes, deserialized);
    }

    // Tests for Isotope struct and its implementations
    #[test]
    fn test_isotope_checked_new() {
        assert!(Isotope::checked_new(Element::H, 1).is_some());
        assert!(Isotope::checked_new(Element::C, 12).is_some());
        assert_eq!(
            Isotope::checked_new(Element::H, 1).unwrap().element,
            Element::H
        );
        assert_eq!(
            Isotope::checked_new(Element::C, 12).unwrap().mass_number,
            12
        );
        assert!(Isotope::checked_new(Element::U, 235).is_none()); // Assuming U-235 is not in our small ISOTOPE_DATA table
        assert!(Isotope::checked_new(Element::H, 0).is_none()); // Invalid mass number
    }

    #[test]
    fn test_isotope_from_symbol_bytes() {
        assert_eq!(
            Isotope::from_symbol_bytes(b"1H"),
            Isotope::checked_new(Element::H, 1)
        );
        assert_eq!(
            Isotope::from_symbol_bytes(b"1h"),
            Isotope::checked_new(Element::H, 1)
        );
        assert_eq!(
            Isotope::from_symbol_bytes(b"12C"),
            Isotope::checked_new(Element::C, 12)
        );
        assert_eq!(
            Isotope::from_symbol_bytes(b"4He"),
            Isotope::checked_new(Element::He, 4)
        );
        assert_eq!(
            Isotope::from_symbol_bytes(b"22Ne"),
            Isotope::checked_new(Element::Ne, 22)
        );

        assert!(Isotope::from_symbol_bytes(b"H1").is_none()); // Invalid format
        assert!(Isotope::from_symbol_bytes(b"H").is_none()); // Element only, not an isotope symbol
        assert!(Isotope::from_symbol_bytes(b"1").is_none()); // Mass number only
        assert!(Isotope::from_symbol_bytes(b"").is_none()); // Empty string
        assert!(Isotope::from_symbol_bytes(b"1X").is_none()); // Invalid element
        assert!(Isotope::from_symbol_bytes(b"235U").is_none()); // Valid format, but U-235 not in ISOTOPE_DATA
        assert!(Isotope::from_symbol_bytes(b"0H").is_none()); // Invalid mass number
    }

    #[test]
    fn test_isotope_from_symbol() {
        assert_eq!(
            Isotope::from_symbol("1H"),
            Isotope::checked_new(Element::H, 1)
        );
        assert_eq!(
            Isotope::from_symbol("1h"),
            Isotope::checked_new(Element::H, 1)
        );
        assert_eq!(
            Isotope::from_symbol("12C"),
            Isotope::checked_new(Element::C, 12)
        );
        assert_eq!(
            Isotope::from_symbol("4He"),
            Isotope::checked_new(Element::He, 4)
        );
        assert_eq!(
            Isotope::from_symbol("22Ne"),
            Isotope::checked_new(Element::Ne, 22)
        );

        assert!(Isotope::from_symbol("H1").is_none()); // Invalid format
        assert!(Isotope::from_symbol("H").is_none()); // Element only, not an isotope symbol
        assert!(Isotope::from_symbol("1").is_none()); // Mass number only
        assert!(Isotope::from_symbol("").is_none()); // Empty string
        assert!(Isotope::from_symbol("1X").is_none()); // Invalid element
        assert!(Isotope::from_symbol("235U").is_none()); // Valid format, but U-235 not in ISOTOPE_DATA
        assert!(Isotope::from_symbol("0H").is_none()); // Invalid mass number
    }

    #[rstest]
    #[case("1H", 1.00782503223, None, "1H")]
    #[case("2H", 2.01410177812, None, "2H")]
    #[case("3H", 3.0160492779, Some(3.8878e8), "3H")]
    #[case("4He", 4.00260325413, None, "4He")]
    #[case("12C", 12.0000000, None, "12C")]
    #[case("13C", 13.00335483507, None, "13C")]
    #[case("20Ne", 19.9924401762, None, "20Ne")]
    fn test_isotope_properties(
        #[case] sym: &str,
        #[case] mass: f64,
        #[case] half_life: Option<f64>,
        #[case] expected_sym: &str,
    ) {
        let isotope = Isotope::from_symbol(sym).unwrap();
        assert!(approx_eq!(f64, isotope.mass(), mass, ulps = 4));
        match (isotope.half_life(), half_life) {
            (Some(v1), Some(v2)) => assert!(approx_eq!(f64, v1, v2, ulps = 4)),
            (None, None) => (),
            _ => panic!("Half-life mismatch"),
        }
        assert_eq!(isotope.symbol(), expected_sym);
    }

    #[test]
    fn test_isotope_display() {
        let isotope_h1 = Isotope::from_symbol("1H").unwrap();
        assert_eq!(isotope_h1.to_string(), "1H");
        let isotope_c12 = Isotope::from_symbol("12C").unwrap();
        assert_eq!(isotope_c12.to_string(), "12C");
    }

    #[test]
    fn test_isotope_from_str() {
        assert_eq!(
            "1H".parse::<Isotope>().unwrap(),
            Isotope::from_symbol("1H").unwrap()
        );
        assert_eq!(
            "12C".parse::<Isotope>().unwrap(),
            Isotope::from_symbol("12C").unwrap()
        );
        assert!("Invalid".parse::<Isotope>().is_err());
        assert!("235U".parse::<Isotope>().is_err()); // Not in ISOTOPE_DATA
    }

    #[test]
    fn test_isotope_try_from_str() {
        assert_eq!(
            Isotope::try_from("1H").unwrap(),
            Isotope::from_symbol("1H").unwrap()
        );
        assert_eq!(
            Isotope::try_from("12C").unwrap(),
            Isotope::from_symbol("12C").unwrap()
        );
        assert!(Isotope::try_from("Invalid").is_err());
        assert!(Isotope::try_from("OH").is_err());
    }

    #[test]
    fn test_isotope_to_element_conversion() {
        let isotope_h1 = Isotope::from_symbol("1H").unwrap();
        assert_eq!(Element::from(isotope_h1), Element::H);
        let isotope_c12 = Isotope::from_symbol("12C").unwrap();
        assert_eq!(Element::from(isotope_c12), Element::C);
    }

    #[test]
    fn test_iso_macro() {
        assert_eq!(iso!("1H"), Isotope::from_symbol("1H").unwrap());
        assert_eq!(iso!("12C"), Isotope::from_symbol("12C").unwrap());
    }
}

/// Represents a unit of time for half-life values.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum TimeUnit {
    Yoctoseconds, // 1e-24 s
    Zeptoseconds, // 1e-21 s
    Attoseconds,  // 1e-18 s
    Femtoseconds, // 1e-15 s
    Picoseconds,  // 1e-12 s
    Nanoseconds,  // 1e-9 s
    Microseconds, // 1e-6 s
    Milliseconds, // 1e-3 s
    Seconds,
    Minutes,
    Hours,
    Days,
    Years,
    KiloYears,    // 1e3 y
    MegaYears,    // 1e6 y
    GigaYears,    // 1e9 y
    TeraYears,    // 1e12 y
    PetaYears,    // 1e15 y
    ExaYears,     // 1e18 y
    ZettaYears,   // 1e21 y
    YottaYears,   // 1e24 y
    ElectronVolts,
    KiloElectronVolts,
    MegaElectronVolts,
}

impl TimeUnit {
    /// Returns the conversion factor to seconds.
    pub fn to_seconds_factor(&self) -> f64 {
        const S_PER_YEAR: f64 = 3.155_695_2e7; // seconds per mean tropical year
        const H_BAR: f64 = 6.582_119_569e-16; // eV·s

        match *self {
            TimeUnit::Yoctoseconds => 1e-24,
            TimeUnit::Zeptoseconds => 1e-21,
            TimeUnit::Attoseconds => 1e-18,
            TimeUnit::Femtoseconds => 1e-15,
            TimeUnit::Picoseconds => 1e-12,
            TimeUnit::Nanoseconds => 1e-9,
            TimeUnit::Microseconds => 1e-6,
            TimeUnit::Milliseconds => 1e-3,
            TimeUnit::Seconds => 1.0,
            TimeUnit::Minutes => 60.0,
            TimeUnit::Hours => 3600.0,
            TimeUnit::Days => 86400.0,
            TimeUnit::Years => S_PER_YEAR,
            TimeUnit::KiloYears => S_PER_YEAR * 1e3,
            TimeUnit::MegaYears => S_PER_YEAR * 1e6,
            TimeUnit::GigaYears => S_PER_YEAR * 1e9,
            TimeUnit::TeraYears => S_PER_YEAR * 1e12,
            TimeUnit::PetaYears => S_PER_YEAR * 1e15,
            TimeUnit::ExaYears => S_PER_YEAR * 1e18,
            TimeUnit::ZettaYears => S_PER_YEAR * 1e21,
            TimeUnit::YottaYears => S_PER_YEAR * 1e24,
            TimeUnit::ElectronVolts => H_BAR, // τ = ħ/Γ
            TimeUnit::KiloElectronVolts => H_BAR / 1e3,
            TimeUnit::MegaElectronVolts => H_BAR / 1e6,
        }
    }
}

/// Represents the half-life of an isotope.
#[derive(Debug, Copy, Clone, PartialEq, Serialize, Deserialize)]
pub struct HalfLife {
    pub value: f64,
    pub unit: TimeUnit,
}

impl HalfLife {
    /// Calculates the half-life in seconds.
    pub fn to_seconds(&self) -> f64 {
        self.value * self.unit.to_seconds_factor()
    }
}

/// Stores properties of a specific isotope.
#[derive(Debug, Copy, Clone, PartialEq, Serialize, Deserialize)]
pub struct IsotopeProperties {
    pub mass: f64,
    pub half_life: Option<HalfLife>,
}
