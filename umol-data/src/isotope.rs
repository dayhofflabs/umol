//! Isotope definitions and data

use crate::half_life::HalfLife;
use crate::isotope_data::{ISOTOPE_DATA, LIGHT_ISOTOPE_MAP};
use crate::Element;
use serde::{Deserialize, Serialize};
use std::fmt::{self, Display};
use std::str::FromStr;
use umol::error::DataError;
use umol::{Error, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, PartialOrd)]
pub enum NamedIsotope {
    D,
    T,
}

impl NamedIsotope {
    // Get named isotope from symbol bytestring (allocation-free)
    pub fn from_symbol_bytes(symbol: &[u8]) -> Option<Self> {
        match symbol {
            b"D" | b"d" => Some(NamedIsotope::D),
            b"T" | b"t" => Some(NamedIsotope::T),
            _ => None,
        }
    }

    // Get named isotope from symbol string (allocation-free)
    pub fn from_symbol(symbol: &str) -> Option<Self> {
        Self::from_symbol_bytes(symbol.as_bytes())
    }

    /// Only hydrogen isotope are named
    pub fn element(&self) -> Element {
        Element::H
    }

    /// Get the isotope mass number
    pub fn mass_number(&self) -> u32 {
        match self {
            NamedIsotope::D => 2,
            NamedIsotope::T => 3,
        }
    }

    /// Get the isotope key (element, mass number) as a 32-bit integer
    pub fn key(&self) -> u32 {
        (Element::H.atomic_number() as u32) << 16 | self.mass_number()
    }

    /// Get the isotope mass (in amu)
    pub fn mass(&self) -> f64 {
        ISOTOPE_DATA.get(&self.key()).unwrap().0
    }

    /// Get the isotope half-life (in s)
    pub fn half_life(&self) -> Option<HalfLife> {
        ISOTOPE_DATA.get(&self.key()).unwrap().1
    }

    /// Get the isotope symbol
    pub fn symbol(&self) -> &'static str {
        match self {
            NamedIsotope::D => "D",
            NamedIsotope::T => "T",
        }
    }

    /// Check if bytestring contains valid named isotope
    pub fn is_named_isotope_bytes(symbol: &[u8]) -> bool {
        symbol == b"D" || symbol == b"d" || symbol == b"T" || symbol == b"t"
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

/// Isotope data is auto-generated from AME2020 and NUBASE2020 and stored in isotope_data.rs.
///
/// 0. isotope mass (in amu)
/// 1. half-life (value, unit)
impl Isotope {
    /// Create a new isotope if it exists in ISOTOPE_DATA.  
    pub fn checked_new(element: Element, mass_number: u32) -> Option<Self> {
        if !Self::is_catalogued(element, mass_number) {
            return None;
        }
        Some(Isotope {
            element,
            mass_number,
        })
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

        let mass_number = atoi::atoi::<u32>(mass_number_bytes)?;
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

    /// Get element
    pub const fn element(&self) -> Element {
        self.element
    }

    /// Get mass number
    pub const fn mass_number(&self) -> u32 {
        self.mass_number
    }

    /// Get key (element, mass number) as a 32-bit integer
    pub const fn key(&self) -> u32 {
        (self.element.atomic_number() as u32) << 16 | self.mass_number
    }

    /// Get isotope mass (in amu)
    pub fn mass(&self) -> f64 {
        ISOTOPE_DATA.get(&self.key()).unwrap().0
    }

    /// Get isotope half-life (in s)
    pub fn half_life(&self) -> Option<HalfLife> {
        ISOTOPE_DATA.get(&self.key()).unwrap().1
    }

    /// Get isotope symbol (AZ notation)
    pub fn symbol(&self) -> String {
        format!("{}{}", self.mass_number, self.element.symbol())
    }

    /// Check if isotope is catalogued (stable or has >= 1 ns half-life)
    pub fn is_catalogued(element: Element, mass_number: u32) -> bool {
        let atomic_number = element.atomic_number() as usize;
        if atomic_number > 0 && atomic_number <= 20 && mass_number < 60 {
            return LIGHT_ISOTOPE_MAP[atomic_number - 1][mass_number as usize];
        }
        ISOTOPE_DATA.contains_key(&((element.atomic_number() as u32) << 16 | mass_number))
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
    use crate::half_life::TimeUnit;
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
    #[case(NamedIsotope::D, Element::H, 2, 2.014102, None, "D")]
    #[case(NamedIsotope::T, Element::H, 3, 3.016049, Some(HalfLife { value: 12.32, unit: TimeUnit::Years }), "T")]
    fn test_named_isotope_properties(
        #[case] named_isotope: NamedIsotope,
        #[case] expected_element: Element,
        #[case] expected_mass_number: u32,
        #[case] expected_mass: f64,
        #[case] expected_half_life: Option<HalfLife>,
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
        assert!(Isotope::checked_new(Element::C, 40).is_none()); // 40C not in ISOTOPE_DATA
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
        assert!(Isotope::from_symbol_bytes(b"40C").is_none()); // Valid format, but 40C not in ISOTOPE_DATA
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
        assert!(Isotope::from_symbol("40C").is_none()); // Valid format, but 40C not in ISOTOPE_DATA
        assert!(Isotope::from_symbol("0H").is_none()); // Invalid mass number
    }

    #[rstest]
    #[case("1H", 1.007825, None, "1H")]
    #[case("2H", 2.014102, None, "2H")]
    #[case("3H", 3.016049, Some(HalfLife { value: 12.32, unit: TimeUnit::Years }), "3H")]
    #[case("4He", 4.002603, None, "4He")]
    #[case("12C", 12.000000, None, "12C")]
    #[case("13C", 13.003355, None, "13C")]
    #[case("20Ne", 19.992440, None, "20Ne")]
    fn test_isotope_properties(
        #[case] sym: &str,
        #[case] mass: f64,
        #[case] half_life: Option<HalfLife>,
        #[case] expected_sym: &str,
    ) {
        let isotope = Isotope::from_symbol(sym).unwrap();
        assert!(approx_eq!(f64, isotope.mass(), mass, ulps = 4));
        match (isotope.half_life(), half_life) {
            (Some(v1), Some(v2)) => assert!(approx_eq!(f64, v1.value, v2.value, ulps = 4)),
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
        assert!("40C".parse::<Isotope>().is_err()); // 40C not in ISOTOPE_DATA
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

    #[rstest]
    #[case(Element::H, 1, true)]
    #[case(Element::H, 2, true)]
    #[case(Element::H, 3, true)]
    #[case(Element::H, 4, false)]
    #[case(Element::C, 12, true)]
    #[case(Element::C, 40, false)]
    fn test_isotope_is_catalogued(
        #[case] element: Element,
        #[case] mass_number: u32,
        #[case] expected: bool,
    ) {
        assert_eq!(Isotope::is_catalogued(element, mass_number), expected);
    }
}
