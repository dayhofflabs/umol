//! Query objects for CTab format.

use serde::{Deserialize, Serialize};
use std::fmt::{self, Display};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum QueryAtom {
    Any,           // * = any atom
    Heavy,         // A = all except H
    Heteroatom,    // Q = any heteroatom (all except H, C)
    Halogen,       // X = F, Cl, Br, I
    Metal, // M = any metal (all except H, He, B, C, N, O, F, Ne, Si, P, S, Cl, Ar, As, Br, Kr, Te, I, Xe, At, Rn)
    HeavyOrH, // AH = any atom (CXSMILES extension, semantically equivalent to Any)
    HeteroatomOrH, // QH = Q or H (CXSMILES extension)
    HalogenOrH, // XH = X or H (CXSMILES extension)
    MetalOrH, // MH = M or H (CXSMILES extension)
}

impl QueryAtom {
    pub fn symbol(&self) -> &str {
        match self {
            QueryAtom::Any => "*",
            QueryAtom::Heavy => "A",
            QueryAtom::Heteroatom => "Q",
            QueryAtom::Halogen => "X",
            QueryAtom::Metal => "M",
            QueryAtom::HeavyOrH => "AH",
            QueryAtom::HeteroatomOrH => "QH",
            QueryAtom::HalogenOrH => "XH",
            QueryAtom::MetalOrH => "MH",
        }
    }

    pub fn from_symbol_bytes(s: &[u8]) -> Option<QueryAtom> {
        match s {
            b"*" => Some(QueryAtom::Any),
            b"A" => Some(QueryAtom::Heavy),
            b"Q" => Some(QueryAtom::Heteroatom),
            b"X" => Some(QueryAtom::Halogen),
            b"M" => Some(QueryAtom::Metal),
            b"AH" => Some(QueryAtom::HeavyOrH),
            b"QH" => Some(QueryAtom::HeteroatomOrH),
            b"XH" => Some(QueryAtom::HalogenOrH),
            b"MH" => Some(QueryAtom::MetalOrH),
            _ => None,
        }
    }

    pub fn from_symbol_str(s: &str) -> Option<QueryAtom> {
        Self::from_symbol_bytes(s.as_bytes())
    }
}

impl Display for QueryAtom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.symbol())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[rstest]
    #[case(b"*", Some(QueryAtom::Any))]
    #[case(b"A", Some(QueryAtom::Heavy))]
    #[case(b"Q", Some(QueryAtom::Heteroatom))]
    #[case(b"X", Some(QueryAtom::Halogen))]
    #[case(b"M", Some(QueryAtom::Metal))]
    #[case(b"AH", Some(QueryAtom::HeavyOrH))]
    #[case(b"QH", Some(QueryAtom::HeteroatomOrH))]
    #[case(b"XH", Some(QueryAtom::HalogenOrH))]
    #[case(b"MH", Some(QueryAtom::MetalOrH))]
    #[case(b"", None)]
    #[case(b" ", None)]
    #[case(b"  ", None)]
    #[case(b"   ", None)]
    fn test_from_symbol_bytes(#[case] input: &[u8], #[case] expected: Option<QueryAtom>) {
        assert_eq!(QueryAtom::from_symbol_bytes(input), expected);
    }

    #[rstest]
    #[case("*", Some(QueryAtom::Any))]
    #[case("A", Some(QueryAtom::Heavy))]
    #[case("Q", Some(QueryAtom::Heteroatom))]
    #[case("X", Some(QueryAtom::Halogen))]
    #[case("M", Some(QueryAtom::Metal))]
    #[case("AH", Some(QueryAtom::HeavyOrH))]
    #[case("QH", Some(QueryAtom::HeteroatomOrH))]
    #[case("XH", Some(QueryAtom::HalogenOrH))]
    #[case("MH", Some(QueryAtom::MetalOrH))]
    #[case("", None)]
    #[case(" ", None)]
    #[case("  ", None)]
    #[case("   ", None)]
    fn test_from_symbol_str(#[case] input: &str, #[case] expected: Option<QueryAtom>) {
        assert_eq!(QueryAtom::from_symbol_str(input), expected);
    }

    #[test]
    fn test_query_atom_serialize() {
        let query_atom = QueryAtom::Any;
        let yaml =
            serde_yaml::to_string(&query_atom).expect("Failed to serialize QueryAtom to YAML");
        let deserialized: QueryAtom =
            serde_yaml::from_str(&yaml).expect("Failed to deserialize QueryAtom from YAML");
        assert_eq!(query_atom, deserialized);
    }
}
