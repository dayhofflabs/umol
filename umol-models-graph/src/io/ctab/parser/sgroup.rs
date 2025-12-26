//! Auxiliary parsers for SGroup properties.

use nom::branch::alt;
use nom::bytes::complete::{tag, tag_no_case, take, take_while_m_n};
use nom::character::complete::{space0, u32 as nom_u32};
use nom::combinator::{map, map_parser, map_res, rest, value};
use nom::sequence::separated_pair;
use nom::{error, AsChar, Err, Parser};

use crate::io::ctab::parser::utils::{fixed_width_partial, to_string};
use crate::simple_ir::{
    SGroupConnectivity, SGroupDataDisplayChars, SGroupDataDisplayPlacement, SGroupDataDisplayType,
    SGroupDataDisplayUnits, SGroupDataType, SGroupMultiplier, SGroupMultiplierOp,
    SGroupMultiplierTerm, SGroupSubtype, SGroupType,
};

/// Parse SGroup type string
pub fn sgroup_type<'a>(
) -> impl Parser<&'a [u8], Output = SGroupType, Error = error::Error<&'a [u8]>> {
    map_res(take(3usize), move |s: &[u8]| match s {
        b"SUP" => Ok(SGroupType::Superatom),
        b"MUL" => Ok(SGroupType::MultipleGroup),
        b"SRU" => Ok(SGroupType::RepeatingUnit),
        b"MON" => Ok(SGroupType::Monomer),
        b"MER" => Ok(SGroupType::Mer),
        b"COP" => Ok(SGroupType::Copolymer),
        b"CRO" => Ok(SGroupType::Crosslink),
        b"MOD" => Ok(SGroupType::Modification),
        b"GRA" => Ok(SGroupType::Graft),
        b"COM" => Ok(SGroupType::Component),
        b"MIX" => Ok(SGroupType::Mixture),
        b"FOR" => Ok(SGroupType::Formulation),
        b"DAT" => Ok(SGroupType::Data),
        b"ANY" => Ok(SGroupType::AnyPolymer),
        b"GEN" => Ok(SGroupType::Generic),
        _ => Err(error::Error::new(s, error::ErrorKind::MapRes)),
    })
}

/// Parse SGroup subtype string
pub fn sgroup_subtype<'a>(
) -> impl Parser<&'a [u8], Output = SGroupSubtype, Error = error::Error<&'a [u8]>> {
    map_res(take(3usize), move |s: &[u8]| match s {
        b"ALT" => Ok(SGroupSubtype::Alternating),
        b"RAN" => Ok(SGroupSubtype::Random),
        b"BLO" => Ok(SGroupSubtype::Block),
        _ => Err(error::Error::new(s, error::ErrorKind::MapRes)),
    })
}

/// Parse SGroup connectivity string
pub fn sgroup_connectivity<'a>(
) -> impl Parser<&'a [u8], Output = SGroupConnectivity, Error = error::Error<&'a [u8]>> {
    move |input: &'a [u8]| {
        fixed_width_partial(
            3usize,
            map_res(rest, move |s: &[u8]| {
                let s = s.trim_ascii();
                match s {
                    b"HH" => Ok(SGroupConnectivity::HeadToHead),
                    b"HT" => Ok(SGroupConnectivity::HeadToTail),
                    b"EU" => Ok(SGroupConnectivity::EitherUnknown),
                    _ => Err(error::Error::new(s, error::ErrorKind::MapRes)),
                }
            }),
            true,
        )
        .parse(input)
        .and_then(|(remaining, opt)| {
            Ok((
                remaining,
                opt.ok_or_else(|| Err::Error(error::Error::new(input, error::ErrorKind::MapRes)))?,
            ))
        })
    }
}

/// Parse SGroup multiplier string
pub fn sgroup_multiplier<'a>(
) -> impl Parser<&'a [u8], Output = SGroupMultiplier, Error = error::Error<&'a [u8]>> {
    alt((
        map(
            separated_pair(
                sgroup_multiplier_term(),
                space0,
                sgroup_multiplier_variable(),
            ),
            |(left, right)| SGroupMultiplier::Expression {
                left,
                op: SGroupMultiplierOp::Mul,
                right,
            },
        ),
        map(
            (
                sgroup_multiplier_term(),
                space0,
                sgroup_multiplier_op(),
                space0,
                sgroup_multiplier_term(),
            ),
            |(left, _, op, _, right)| SGroupMultiplier::Expression { left, op, right },
        ),
        map(sgroup_multiplier_term(), SGroupMultiplier::Single),
    ))
}

/// Parse an integer multiplier
fn sgroup_multiplier_integer<'a>(
) -> impl Parser<&'a [u8], Output = SGroupMultiplierTerm, Error = error::Error<&'a [u8]>> {
    map(nom_u32, SGroupMultiplierTerm::Integer)
}

/// Parse a single-character variable multiplier
fn sgroup_multiplier_variable<'a>(
) -> impl Parser<&'a [u8], Output = SGroupMultiplierTerm, Error = error::Error<&'a [u8]>> {
    map(take_while_m_n(1, 1, AsChar::is_alpha), |s: &[u8]| {
        SGroupMultiplierTerm::Variable(s[0] as char)
    })
}

/// Parse a single multiplier term (variable or integer)
fn sgroup_multiplier_term<'a>(
) -> impl Parser<&'a [u8], Output = SGroupMultiplierTerm, Error = error::Error<&'a [u8]>> {
    alt((sgroup_multiplier_integer(), sgroup_multiplier_variable()))
}

/// Parse arithmetic operator
fn sgroup_multiplier_op<'a>(
) -> impl Parser<&'a [u8], Output = SGroupMultiplierOp, Error = error::Error<&'a [u8]>> {
    alt((
        value(SGroupMultiplierOp::Add, tag("+")),
        value(SGroupMultiplierOp::Sub, tag("-")),
        value(SGroupMultiplierOp::Mul, tag("*")),
        value(SGroupMultiplierOp::Div, tag("/")),
    ))
}

/// Parse SGroup subscript string
pub fn sgroup_subscript<'a>(
) -> impl Parser<&'a [u8], Output = String, Error = error::Error<&'a [u8]>> {
    map_res(rest, move |s: &[u8]| to_string(s))
}

// Parse SGroup data type string
pub fn sgroup_data_type<'a>(
) -> impl Parser<&'a [u8], Output = SGroupDataType, Error = error::Error<&'a [u8]>> {
    map_res(fixed_width_partial(2usize, rest, true), move |s| {
        if let Some(s) = s {
            let s = s.trim_ascii();
            match s {
                b"T" | b"" => Ok(SGroupDataType::Text),
                b"F" => Ok(SGroupDataType::Formatted),
                b"N" => Ok(SGroupDataType::Numeric),
                _ => Err(error::Error::new(s, error::ErrorKind::MapRes)),
            }
        } else {
            Ok(SGroupDataType::Text)
        }
    })
}

// Parse SGroup data display type string
pub fn sgroup_data_display_type<'a>(
) -> impl Parser<&'a [u8], Output = SGroupDataDisplayType, Error = error::Error<&'a [u8]>> {
    map_res(take(1usize), move |s: &[u8]| match s {
        b"A" => Ok(SGroupDataDisplayType::Attached),
        b"D" => Ok(SGroupDataDisplayType::Detached),
        _ => Err(error::Error::new(s, error::ErrorKind::MapRes)),
    })
}

// Parse SGroup data display placement string
pub fn sgroup_data_display_placement<'a>(
) -> impl Parser<&'a [u8], Output = SGroupDataDisplayPlacement, Error = error::Error<&'a [u8]>> {
    map_res(take(1usize), move |s: &[u8]| match s {
        b"A" => Ok(SGroupDataDisplayPlacement::Absolute),
        b"R" => Ok(SGroupDataDisplayPlacement::Relative),
        _ => Err(error::Error::new(s, error::ErrorKind::MapRes)),
    })
}

// Parse SGroup data display units string
pub fn sgroup_data_display_units<'a>(
) -> impl Parser<&'a [u8], Output = SGroupDataDisplayUnits, Error = error::Error<&'a [u8]>> {
    map_res(take(1usize), move |s: &[u8]| match s {
        b" " => Ok(SGroupDataDisplayUnits::None),
        b"U" => Ok(SGroupDataDisplayUnits::DisplayUnits),
        _ => Err(error::Error::new(s, error::ErrorKind::MapRes)),
    })
}

// Parse SGroup data display chars string
pub fn sgroup_data_display_chars<'a>(
) -> impl Parser<&'a [u8], Output = SGroupDataDisplayChars, Error = error::Error<&'a [u8]>> {
    map_parser(take(3usize), move |s: &'a [u8]| {
        let s = s.trim_ascii();
        alt((
            value(SGroupDataDisplayChars::All, tag_no_case("ALL")),
            map(nom_u32, SGroupDataDisplayChars::Number),
        ))
        .parse(s)
    })
}

#[cfg(test)]
mod tests {
    use nom::Err;
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case(b"SUP", SGroupType::Superatom)]
    #[case(b"MUL", SGroupType::MultipleGroup)]
    #[case(b"SRU", SGroupType::RepeatingUnit)]
    #[case(b"MON", SGroupType::Monomer)]
    #[case(b"MER", SGroupType::Mer)]
    #[case(b"COP", SGroupType::Copolymer)]
    #[case(b"CRO", SGroupType::Crosslink)]
    #[case(b"MOD", SGroupType::Modification)]
    #[case(b"GRA", SGroupType::Graft)]
    #[case(b"COM", SGroupType::Component)]
    #[case(b"MIX", SGroupType::Mixture)]
    #[case(b"FOR", SGroupType::Formulation)]
    #[case(b"DAT", SGroupType::Data)]
    #[case(b"ANY", SGroupType::AnyPolymer)]
    #[case(b"GEN", SGroupType::Generic)]
    fn test_sgroup_type(#[case] input: &[u8], #[case] expected: SGroupType) {
        let (remaining, result) = sgroup_type().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case(b"XYZ", "unknown type", error::ErrorKind::MapRes)]
    #[case(b"SU", "too short", error::ErrorKind::Eof)]
    fn test_sgroup_type_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: error::ErrorKind,
    ) {
        let result = sgroup_type().parse(input);
        assert!(result.is_err());
        assert!(
            matches!(result.as_ref(), Err(Err::Error(e)) if e.code == expected_kind),
            "Expected {:?} error for {}, got {:?}",
            expected_kind,
            desc,
            result
        );
    }

    #[rstest]
    #[case(b"ALT", SGroupSubtype::Alternating)]
    #[case(b"RAN", SGroupSubtype::Random)]
    #[case(b"BLO", SGroupSubtype::Block)]
    fn test_sgroup_subtype(#[case] input: &[u8], #[case] expected: SGroupSubtype) {
        let (remaining, result) = sgroup_subtype().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case(b"XYZ", "unknown subtype", error::ErrorKind::MapRes)]
    #[case(b"SU", "too short", error::ErrorKind::Eof)]
    fn test_sgroup_subtype_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: error::ErrorKind,
    ) {
        let result = sgroup_subtype().parse(input);
        assert!(result.is_err());
        assert!(
            matches!(result.as_ref(), Err(Err::Error(e)) if e.code == expected_kind),
            "Expected {:?} error for {}, got {:?}",
            expected_kind,
            desc,
            result
        );
    }

    #[rstest]
    #[case(b"HH ", SGroupConnectivity::HeadToHead)]
    #[case(b"HT ", SGroupConnectivity::HeadToTail)]
    #[case(b"EU ", SGroupConnectivity::EitherUnknown)]
    #[case(b" HH", SGroupConnectivity::HeadToHead)]
    fn test_sgroup_connectivity(#[case] input: &[u8], #[case] expected: SGroupConnectivity) {
        let (remaining, result) = sgroup_connectivity().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case(b"XY ", "unknown connectivity", error::ErrorKind::MapRes)]
    fn test_sgroup_connectivity_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: error::ErrorKind,
    ) {
        let result = sgroup_connectivity().parse(input);
        assert!(result.is_err());
        assert!(
            matches!(result.as_ref(), Err(Err::Error(e)) if e.code == expected_kind),
            "Expected {:?} error for {}, got {:?}",
            expected_kind,
            desc,
            result
        );
    }

    #[rstest]
    #[case(b"N", SGroupMultiplier::Single(SGroupMultiplierTerm::Variable('N')))]
    #[case(b"n", SGroupMultiplier::Single(SGroupMultiplierTerm::Variable('n')))]
    #[case(b"m", SGroupMultiplier::Single(SGroupMultiplierTerm::Variable('m')))]
    #[case(b"1", SGroupMultiplier::Single(SGroupMultiplierTerm::Integer(1)))]
    #[case(b"2", SGroupMultiplier::Single(SGroupMultiplierTerm::Integer(2)))]
    #[case(b"2n", SGroupMultiplier::Expression {left: SGroupMultiplierTerm::Integer(2), op: SGroupMultiplierOp::Mul, right: SGroupMultiplierTerm::Variable('n')})]
    #[case(b"n+1", SGroupMultiplier::Expression {left: SGroupMultiplierTerm::Variable('n'), op: SGroupMultiplierOp::Add, right: SGroupMultiplierTerm::Integer(1)})]
    #[case(b"n*m", SGroupMultiplier::Expression {left: SGroupMultiplierTerm::Variable('n'), op: SGroupMultiplierOp::Mul, right: SGroupMultiplierTerm::Variable('m')})]
    #[case(b"n m", SGroupMultiplier::Expression {left: SGroupMultiplierTerm::Variable('n'), op: SGroupMultiplierOp::Mul, right: SGroupMultiplierTerm::Variable('m')})]
    #[case(b"2 m", SGroupMultiplier::Expression {left: SGroupMultiplierTerm::Integer(2), op: SGroupMultiplierOp::Mul, right: SGroupMultiplierTerm::Variable('m')})]
    #[case(b"X", SGroupMultiplier::Single(SGroupMultiplierTerm::Variable('X')))]
    fn test_sgroup_multiplier(#[case] input: &[u8], #[case] expected: SGroupMultiplier) {
        let (remaining, result) = sgroup_multiplier().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case(b"@", "invalid symbol", error::ErrorKind::TakeWhileMN)]
    fn test_sgroup_multiplier_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: error::ErrorKind,
    ) {
        let result = sgroup_multiplier().parse(input);
        assert!(result.is_err());
        assert!(
            matches!(result.as_ref(), Err(Err::Error(e)) if e.code == expected_kind),
            "Expected {:?} error for {}, got {:?}",
            expected_kind,
            desc,
            result
        );
    }

    #[rstest]
    #[case(b"F ", SGroupDataType::Formatted)]
    #[case(b"N ", SGroupDataType::Numeric)]
    #[case(b"T ", SGroupDataType::Text)]
    #[case(b"  ", SGroupDataType::Text)]
    fn test_sgroup_data_type(#[case] input: &[u8], #[case] expected: SGroupDataType) {
        let (remaining, result) = sgroup_data_type().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case(b"X", "unknown data type", error::ErrorKind::MapRes)]
    fn test_sgroup_data_type_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: error::ErrorKind,
    ) {
        let result = sgroup_data_type().parse(input);
        assert!(result.is_err());
        assert!(
            matches!(result.as_ref(), Err(Err::Error(e)) if e.code == expected_kind),
            "Expected {:?} error for {}, got {:?}",
            expected_kind,
            desc,
            result
        );
    }

    #[rstest]
    #[case(b"A", SGroupDataDisplayType::Attached)]
    #[case(b"D", SGroupDataDisplayType::Detached)]
    fn test_sgroup_data_display_type(
        #[case] input: &[u8],
        #[case] expected: SGroupDataDisplayType,
    ) {
        let (remaining, result) = sgroup_data_display_type().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case(b"X", "unknown display type", error::ErrorKind::MapRes)]
    fn test_sgroup_data_display_type_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: error::ErrorKind,
    ) {
        let result = sgroup_data_display_type().parse(input);
        assert!(result.is_err());
        assert!(
            matches!(result.as_ref(), Err(Err::Error(e)) if e.code == expected_kind),
            "Expected {:?} error for {}, got {:?}",
            expected_kind,
            desc,
            result
        );
    }

    #[rstest]
    #[case(b"A", SGroupDataDisplayPlacement::Absolute)]
    #[case(b"R", SGroupDataDisplayPlacement::Relative)]
    fn test_sgroup_data_display_placement(
        #[case] input: &[u8],
        #[case] expected: SGroupDataDisplayPlacement,
    ) {
        let (remaining, result) = sgroup_data_display_placement().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case(b"X", "unknown display placement", error::ErrorKind::MapRes)]
    fn test_sgroup_data_display_placement_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: error::ErrorKind,
    ) {
        let result = sgroup_data_display_placement().parse(input);
        assert!(result.is_err());
        assert!(
            matches!(result.as_ref(), Err(Err::Error(e)) if e.code == expected_kind),
            "Expected {:?} error for {}, got {:?}",
            expected_kind,
            desc,
            result
        );
    }

    #[rstest]
    #[case(b" ", SGroupDataDisplayUnits::None)]
    #[case(b"U", SGroupDataDisplayUnits::DisplayUnits)]
    fn test_sgroup_data_display_units(
        #[case] input: &[u8],
        #[case] expected: SGroupDataDisplayUnits,
    ) {
        let (remaining, result) = sgroup_data_display_units().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case(b"X", "unknown display units", error::ErrorKind::MapRes)]
    fn test_sgroup_data_display_units_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: error::ErrorKind,
    ) {
        let result = sgroup_data_display_units().parse(input);
        assert!(result.is_err());
        assert!(
            matches!(result.as_ref(), Err(Err::Error(e)) if e.code == expected_kind),
            "Expected {:?} error for {}, got {:?}",
            expected_kind,
            desc,
            result
        );
    }

    #[rstest]
    #[case(b"ALL", SGroupDataDisplayChars::All)]
    #[case(b"  1", SGroupDataDisplayChars::Number(1))]
    fn test_sgroup_data_display_chars(
        #[case] input: &[u8],
        #[case] expected: SGroupDataDisplayChars,
    ) {
        let (remaining, result) = sgroup_data_display_chars().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case(b"X", "unknown display chars", error::ErrorKind::Eof)]
    fn test_sgroup_data_display_chars_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: error::ErrorKind,
    ) {
        let result = sgroup_data_display_chars().parse(input);
        assert!(result.is_err());
        assert!(
            matches!(result.as_ref(), Err(Err::Error(e)) if e.code == expected_kind),
            "Expected {:?} error for {}, got {:?}",
            expected_kind,
            desc,
            result
        );
    }
}
