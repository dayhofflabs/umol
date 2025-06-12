//! Properties block parser for CTab files.

use nom::{
    branch::alt,
    bytes::complete::tag,
    character::complete::space0,
    combinator::{all_consuming, map, peek},
    error,
    multi::length_count,
    sequence::{delimited, preceded},
    Parser,
};

use super::utils::{fixed_width_int, fixed_width_int_in_range, fixed_width_int_minus1};

#[derive(Debug, Clone, PartialEq)]
pub struct ChargeEntry {
    pub atom_index: usize,
    pub charge: i8,
}

#[derive(Debug, Clone, PartialEq)]
pub struct RadicalEntry {
    pub atom_index: usize,
    pub radical_type: i8,
}

#[derive(Debug, Clone, PartialEq)]
pub struct IsotopeEntry {
    pub atom_index: usize,
    pub mass: u32,
}

/// An enum representing a parsed property modification, containing the raw data.
/// This avoids allocating a new Vec for every single property line in a file.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum PropertyEntries {
    ChargeEntries(Vec<ChargeEntry>),
    RadicalEntries(Vec<RadicalEntry>),
    IsotopeEntries(Vec<IsotopeEntry>),
}

/// Parse charge property line (*[Generic]*).
/// M  CHGnn8 aaa vvv ...
/// vvv: -15..= 15.
fn charge_line<'a>() -> impl Parser<&'a [u8], Output = PropertyEntries, Error = error::Error<&'a [u8]>> {
    all_consuming(map(
        delimited(
            tag("M  CHG"),
            length_count(
                fixed_width_int_in_range::<u8, _>(3, 1..=8),
                map(
                    (
                        preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                        preceded(tag(" "), fixed_width_int_in_range::<i8, _>(3, -15..=15)),
                    ),
                    |(atom_index, charge)| ChargeEntry { atom_index, charge },
                ),
            ),
            space0,
        ),
        PropertyEntries::ChargeEntries,
    ))
}

/// Parse radical property line (*[Generic]*).
/// M  RADnn8 aaa vvv ...
/// vvv: 0..= 3: 0 = no radical, 1 = singlet (:), 2 = doublet (. or ^), 3 = triplet (^^).
fn radical_line<'a>() -> impl Parser<&'a [u8], Output = PropertyEntries, Error = error::Error<&'a [u8]>>
{
    all_consuming(map(
        delimited(
            tag("M  RAD"),
            length_count(
                fixed_width_int_in_range::<u8, _>(3, 1..=8),
                map(
                    (
                        preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                        preceded(tag(" "), fixed_width_int_in_range::<i8, _>(3, 0..=3)),
                    ),
                    |(atom_index, radical_type)| RadicalEntry {
                        atom_index,
                        radical_type,
                    },
                ),
            ),
            space0,
        ),
        PropertyEntries::RadicalEntries,
    ))
}

/// Parse isotope property line (*[Generic]*).
/// M  ISOnn8 aaa vvv ...
/// vvv: isotope mass number (not difference)
/// Difference between the isotope mass number and reference isotope mass number
/// should be in the range -18..=12.
fn isotope_line<'a>() -> impl Parser<&'a [u8], Output = PropertyEntries, Error = error::Error<&'a [u8]>>
{
    all_consuming(map(
        delimited(
            tag("M  ISO"),
            length_count(
                fixed_width_int_in_range::<u8, _>(3, 1..=8),
                map(
                    (
                        preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                        preceded(tag(" "), fixed_width_int::<u32>(3)),
                    ),
                    |(atom_index, mass)| IsotopeEntry { atom_index, mass },
                ),
            ),
            space0,
        ),
        PropertyEntries::IsotopeEntries,
    ))
}

/// Parse property line
pub(crate) fn property_line<'a>(
) -> impl Parser<&'a [u8], Output = PropertyEntries, Error = error::Error<&'a [u8]>> {
    alt((
        preceded(peek(tag("M  CHG")), charge_line()),
        preceded(peek(tag("M  RAD")), radical_line()),
        preceded(peek(tag("M  ISO")), isotope_line()),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use nom::{error::ErrorKind, Err, Parser};
    use rstest::rstest;

    #[rstest]
    #[case(b"M  CHG  1   1  -1", PropertyEntries::ChargeEntries(vec![ChargeEntry { atom_index: 0, charge: -1 }]))]
    #[case(b"M  CHG  2   1  -1   4   1", PropertyEntries::ChargeEntries(vec![ChargeEntry { atom_index: 0, charge: -1 }, ChargeEntry { atom_index: 3, charge: 1 }]))]
    #[case(b"M  CHG  8   1   1   2   2   3   3   4   4   5   5   6   6   7   7   8   8",
        PropertyEntries::ChargeEntries(vec![
            ChargeEntry { atom_index: 0, charge: 1 },
            ChargeEntry { atom_index: 1, charge: 2 },
            ChargeEntry { atom_index: 2, charge: 3 },
            ChargeEntry { atom_index: 3, charge: 4 },
            ChargeEntry { atom_index: 4, charge: 5 },
            ChargeEntry { atom_index: 5, charge: 6 },
            ChargeEntry { atom_index: 6, charge: 7 },
            ChargeEntry { atom_index: 7, charge: 8 },
        ])
    )]
    #[case(b"M  CHG  1  25  15", PropertyEntries::ChargeEntries(vec![ChargeEntry { atom_index: 24, charge: 15 }]))]
    fn test_charge_line(#[case] input: &[u8], #[case] expected: PropertyEntries) {
        let (remaining, result) = charge_line().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case(b"M  CHG  1   1  -1  a", "trailing chars", ErrorKind::Eof)]
    #[case(b"M  CHG  2   1  -1", "count does not match item list", ErrorKind::Tag)]
    #[case(
        b"M  CHG  1   1  -1   4   1",
        "item list longer than count",
        ErrorKind::Eof
    )]
    #[case(b"M  CHG  0", "count is zero", ErrorKind::Verify)]
    #[case(b"M  CHG  1   0 -10", "atom index is zero", ErrorKind::Verify)]
    #[case(b"M  XXX  1   1  -1", "invalid property tag", ErrorKind::Tag)]
    #[case(b"X  CHG  1   1  -1", "invalid prefix", ErrorKind::Tag)]
    fn test_charge_line_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let result = charge_line().parse(input);
        assert!(result.is_err(), "{} should have failed", desc);
        assert!(
            matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "Mismatched error kind for {}, expected {:?}, got {}",
            desc,
            expected_kind,
            result.clone().unwrap_err().map(|e| e.code),
        );
    }

    #[rstest]
    #[case(b"M  RAD  1   1   2", PropertyEntries::RadicalEntries(vec![RadicalEntry { atom_index: 0, radical_type: 2 }]))]
    #[case(b"M  RAD  2   1   1   4   3", PropertyEntries::RadicalEntries(vec![RadicalEntry { atom_index: 0, radical_type: 1 }, RadicalEntry { atom_index: 3, radical_type: 3 }]))]
    fn test_radical_line(#[case] input: &[u8], #[case] expected: PropertyEntries) {
        let (remaining, result) = radical_line().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case(b"M  RAD  1   1   4", "value out of range", ErrorKind::Verify)]
    #[case(b"M  RAD  1   1  -1", "value out of range", ErrorKind::Verify)]
    fn test_radical_line_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let result = radical_line().parse(input);
        assert!(result.is_err(), "{} should have failed", desc);
        assert!(
            matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "Mismatched error kind for {}, expected {:?}, got {}",
            desc,
            expected_kind,
            result.clone().unwrap_err().map(|e| e.code),
        );
    }

    #[rstest]
    #[case(b"M  ISO  1   1  13", PropertyEntries::IsotopeEntries(vec![IsotopeEntry { atom_index: 0, mass: 13 }]))]
    #[case(b"M  ISO  2  12   2  15  14", PropertyEntries::IsotopeEntries(vec![IsotopeEntry { atom_index: 11, mass: 2 }, IsotopeEntry { atom_index: 14, mass: 14 }]))]
    fn test_isotope_line(#[case] input: &[u8], #[case] expected: PropertyEntries) {
        let (remaining, result) = isotope_line().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case(b"M  CHG  1   1  -1", PropertyEntries::ChargeEntries(vec![ChargeEntry { atom_index: 0, charge: -1 }]))]
    #[case(b"M  RAD  1   1   2", PropertyEntries::RadicalEntries(vec![RadicalEntry { atom_index: 0, radical_type: 2 }]))]
    #[case(b"M  ISO  1   1  13", PropertyEntries::IsotopeEntries(vec![IsotopeEntry { atom_index: 0, mass: 13 }]))]
    fn test_property_line(#[case] input: &[u8], #[case] expected: PropertyEntries) {
        let (remaining, result) = property_line().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }
}
