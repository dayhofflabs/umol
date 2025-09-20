//! Shared SMILES linter utilities: bracket helpers and style helpers.

use bstr::ByteSlice;
use memchr::{memchr, memchr2_iter, memchr_iter};
use regex::Regex;

use crate::io::smiles::parser::utils::{
    parse_bracket as parse_bracket_parser, BracketField, BracketFields,
};

/// Returns the index of the first ']' at or after `idx`.
pub fn find_closing_bracket(bytes: &[u8], idx: usize) -> Option<usize> {
    bytes
        .get(idx..)
        .and_then(|r| memchr(b']', r))
        .map(|o| idx + o)
}

pub fn is_bare_organic(s: &str) -> bool {
    matches!(
        s,
        "B" | "C" | "N" | "O" | "S" | "P" | "F" | "Cl" | "Br" | "I"
    )
}

pub fn contains_h1(s: &str) -> bool {
    s.as_bytes().contains_str(b"H1")
}

pub fn find_subslice(hay: &str, needle: &str) -> Option<(usize, usize)> {
    hay.find(needle).map(|s| (s, s + needle.len()))
}

pub fn find_charge_plus_minus_one(s: &str) -> Option<(usize, usize)> {
    let bytes = s.as_bytes();
    for idx in memchr2_iter(b'+', b'-', bytes) {
        if idx + 1 < bytes.len() && bytes[idx + 1] == b'1' {
            let end = idx + 2;
            if end >= bytes.len() || !bytes[end].is_ascii_digit() {
                return Some((idx, end));
            }
        }
    }
    None
}

/// Returns true if bracket tail fields are out of canonical order.
/// Canonical order: Chirality -> H count -> Charge -> Class.
/// Uses the shared bracket parser once and checks the order of parsed fields.
pub fn bracket_order_misordered(s: &str) -> bool {
    let (_elem, _iso, tails) = parse_bracket_parser(s);
    let mut prev: Option<usize> = None;
    for f in tails.iter() {
        let id = match f {
            BracketField::Chiral(_) => 0,
            BracketField::HydrogenCount(_) => 1,
            BracketField::Charge(_) => 2,
            BracketField::Class(_) => 3,
        };
        if let Some(p) = prev {
            if id < p {
                return true;
            }
        }
        prev = Some(id);
    }
    false
}

pub fn find_h_two_digits(s: &str) -> Option<(usize, usize)> {
    let re = Regex::new(r"H\d{2,}").unwrap();
    re.find(s).map(|m| (m.start(), m.end()))
}

/// Detects invalid class index after ':' in a bracket tail.
/// Returns (start, end, is_negative) where:
/// - is_negative = true if immediately followed by '-'
/// - otherwise true if non-digit or end-of-input after ':'
pub fn invalid_class_index(s: &str) -> Option<(usize, usize, bool)> {
    let bytes = s.as_bytes();
    for idx in memchr_iter(b':', bytes) {
        let next = idx + 1;
        if next >= bytes.len() {
            return Some((idx, idx + 1, false));
        }
        match bytes[next] {
            b'-' => return Some((idx, next + 1, true)),
            d if !d.is_ascii_digit() => return Some((idx, next + 1, false)),
            _ => {}
        }
    }
    None
}

pub fn parse_bracket(inner: &str) -> BracketFields {
    let (elem, iso, tails) = parse_bracket_parser(inner);
    let mut parsed = BracketFields::default();
    parsed.element = elem;
    parsed.isotope = iso;
    for f in tails {
        match f {
            BracketField::Chiral(_) => {}
            BracketField::HydrogenCount(v) => parsed.hcount = Some(v),
            BracketField::Charge(q) => parsed.charge = Some(q),
            BracketField::Class(c) => parsed.class = Some(c),
        }
    }
    parsed
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case(b"[CH3]".as_slice(), 1, Some(4))]
    #[case(b"[CH3".as_slice(), 1, None)]
    fn test_find_closing_bracket(
        #[case] s: &[u8],
        #[case] start: usize,
        #[case] expect: Option<usize>,
    ) {
        assert_eq!(find_closing_bracket(s, start), expect);
    }

    #[rstest]
    #[case("+1", Some((0, 2)))]
    #[case("-1]", Some((0, 2)))]
    #[case("+10", None)]
    #[case("+2", None)]
    fn test_find_charge_plus_minus_one(#[case] s: &str, #[case] expect: Option<(usize, usize)>) {
        assert_eq!(find_charge_plus_minus_one(s), expect);
    }

    #[rstest]
    #[case(":", Some((0, 1, false)))]
    #[case(":x", Some((0, 2, false)))]
    #[case(":-1", Some((0, 2, true)))]
    #[case(":12", None)]
    fn test_invalid_class_index(#[case] s: &str, #[case] expect: Option<(usize, usize, bool)>) {
        assert_eq!(invalid_class_index(s), expect);
    }

    #[rstest]
    // Misordered examples
    #[case("*H@", true)]
    #[case("CH@", true)]
    #[case("*+1H", true)]
    #[case("C:1+1", true)]
    // Ordered examples
    #[case("*@H", false)]
    #[case("C@H", false)]
    #[case("C@H+1:1", false)]
    #[case("CH+1:2", false)]
    fn test_bracket_order_misordered(#[case] inner: &str, #[case] misordered: bool) {
        assert_eq!(bracket_order_misordered(inner), misordered);
    }
}
