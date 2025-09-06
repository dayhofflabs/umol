//! R-group type for CTab format.

use std::fmt::{self, Display};

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum RGroupOccurrence {
    Exactly(u8),
    Range(u8, u8),   // Inclusive
    GreaterThan(u8), // Default is > 0
    FewerThan(u8),
}

impl RGroupOccurrence {
    pub fn contains(&self, count: u8) -> bool {
        match self {
            RGroupOccurrence::Exactly(n) => *n == count,
            RGroupOccurrence::Range(n, m) => count >= *n && count <= *m,
            RGroupOccurrence::GreaterThan(n) => count > *n,
            RGroupOccurrence::FewerThan(n) => count < *n,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RGroup {
    pub label: Option<u32>,
    pub dependent_label: Option<u32>,
    pub rgroup_or_h: bool,
    pub occurrence: Vec<RGroupOccurrence>,
}

impl RGroup {
    pub fn new(label: Option<u32>) -> Self {
        Self {
            label,
            dependent_label: None,
            rgroup_or_h: false,
            occurrence: vec![RGroupOccurrence::GreaterThan(0)],
        }
    }

    pub fn from_symbol_bytes(input: &[u8]) -> Option<Self> {
        debug_assert!(input.len() <= 3, "R-group symbol must be 1-3 characters");

        if input.is_empty() || input[0] != b'R' {
            None
        } else if input.len() == 1 || input.len() == 2 && input[1] == b'#' {
            Some(Self::new(None))
        } else {
            let num_str = &input[1..];
            if num_str.len() == 1 {
                if num_str[0] < b'0' || num_str[0] > b'9' {
                    None
                } else {
                    let label = (num_str[0] - b'0') as u32;
                    if label == 0 {
                        Some(Self::new(None))
                    } else {
                        Some(Self::new(Some(label)))
                    }
                }
            } else if num_str.len() == 2 {
                if num_str[0] < b'0' || num_str[0] > b'9' || num_str[1] < b'0' || num_str[1] > b'9'
                {
                    None
                } else {
                    let label = ((num_str[0] - b'0') * 10 + (num_str[1] - b'0')) as u32;
                    if label == 0 {
                        Some(Self::new(None))
                    } else {
                        Some(Self::new(Some(label)))
                    }
                }
            } else {
                None
            }
        }
    }

    pub fn from_symbol_str(input: &str) -> Option<Self> {
        Self::from_symbol_bytes(input.as_bytes())
    }
}

impl Display for RGroup {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.label.is_some() {
            write!(f, "R{}", self.label.unwrap_or(0))
        } else {
            write!(f, "R#")
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case(b"R", RGroup::new(None))]
    #[case(b"R#", RGroup::new(None))]
    #[case(b"R0", RGroup::new(None))]
    #[case(b"R1", RGroup::new(Some(1)))]
    #[case(b"R12", RGroup::new(Some(12)))]
    fn test_from_symbol_bytes(#[case] input: &[u8], #[case] expected: RGroup) {
        let symbol = RGroup::from_symbol_bytes(input);
        assert!(symbol.is_some());
        assert_eq!(symbol.unwrap(), expected);
    }

    #[rstest]
    #[case("R", RGroup::new(None))]
    #[case("R#", RGroup::new(None))]
    #[case("R0", RGroup::new(None))]
    #[case("R1", RGroup::new(Some(1)))]
    #[case("R12", RGroup::new(Some(12)))]
    fn test_from_symbol_str(#[case] input: &str, #[case] expected: RGroup) {
        let symbol = RGroup::from_symbol_str(input);
        assert!(symbol.is_some());
        assert_eq!(symbol.unwrap(), expected);
    }

    #[rstest]
    #[case(b"Q")]
    #[case(b"R-1")]
    #[case(b"R#1")]
    fn test_from_symbol_bytes_invalid(#[case] input: &[u8]) {
        let symbol = RGroup::from_symbol_bytes(input);
        assert!(symbol.is_none());
    }

    #[test]
    fn test_rgroup_serialize() {
        let rgroup = RGroup::new(Some(1));
        let yaml = serde_yaml::to_string(&rgroup).expect("Failed to serialize RGroup to YAML");
        let deserialized: RGroup =
            serde_yaml::from_str(&yaml).expect("Failed to deserialize RGroup from YAML");
        assert_eq!(rgroup, deserialized);
    }
}
