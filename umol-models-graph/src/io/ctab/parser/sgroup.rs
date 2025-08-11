//! Auxiliary parsers for SGroup properties.

use crate::io::ctab::parser::utils::{trim_whitespace_2char, trim_whitespace_3char};
use crate::io::ctab::sgroup::{SGroupConnectivity, SGroupMultiplier, SGroupSubtype, SGroupType};
use nom::branch::alt;
use nom::bytes::complete::{tag, take};
use nom::combinator::{map, map_parser, map_res};
use nom::number::complete::u8 as nom_u8;
use nom::{error, Parser};

/// Parse SGroup type string
pub fn sgroup_type<'a>(
) -> impl Parser<&'a [u8], Output = SGroupType, Error = error::Error<&'a [u8]>> {
    map_res(take(3usize), |s| {
        let s = trim_whitespace_3char(s);
        match s {
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
        }
    })
}

/// Parse SGroup subtype string
pub fn sgroup_subtype<'a>(
) -> impl Parser<&'a [u8], Output = SGroupSubtype, Error = error::Error<&'a [u8]>> {
    map_res(take(3usize), |s| {
        let s = trim_whitespace_3char(s);
        match s {
            b"ALT" => Ok(SGroupSubtype::Alternating),
            b"RAN" => Ok(SGroupSubtype::Random),
            b"BLO" => Ok(SGroupSubtype::Block),
            _ => Err(error::Error::new(s, error::ErrorKind::MapRes)),
        }
    })
}

/// Parse SGroup connectivity string
pub fn sgroup_connectivity<'a>(
) -> impl Parser<&'a [u8], Output = SGroupConnectivity, Error = error::Error<&'a [u8]>> {
    map_res(take(2usize), |s| {
        let s = trim_whitespace_2char(s);
        match s {
            b"HH" => Ok(SGroupConnectivity::HeadToHead),
            b"HT" => Ok(SGroupConnectivity::HeadToTail),
            b"EU" => Ok(SGroupConnectivity::EitherUnknown),
            _ => Err(error::Error::new(s, error::ErrorKind::MapRes)),
        }
    })
}

/// Parse SGroup multiplier string
pub fn sgroup_multiplier<'a>(
) -> impl Parser<&'a [u8], Output = SGroupMultiplier, Error = error::Error<&'a [u8]>> {
    map_parser(
        take(1usize),
        alt((
            map(tag("N"), |_| SGroupMultiplier::N),
            map(tag("n"), |_| SGroupMultiplier::N),
            map(nom_u8, |count| SGroupMultiplier::Count(count as u32)),
        )),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use nom::error::ErrorKind;
    use rstest::*;

    #[rstest]
    #[case("SUP", SGroupType::Superatom)]
    #[case("MUL", SGroupType::MultipleGroup)]
    #[case("SRU", SGroupType::RepeatingUnit)]
    #[case("MON", SGroupType::Monomer)]
    #[case("MER", SGroupType::Mer)]
    #[case("COP", SGroupType::Copolymer)]
    #[case("CRO", SGroupType::Crosslink)]
    #[case("MOD", SGroupType::Modification)]
    #[case("GRA", SGroupType::Graft)]
    #[case("COM", SGroupType::Component)]
    #[case("MIX", SGroupType::Mixture)]
    #[case("FOR", SGroupType::Formulation)]
    #[case("DAT", SGroupType::Data)]
    #[case("ANY", SGroupType::AnyPolymer)]
    #[case("GEN", SGroupType::Generic)]
    fn test_sgroup_type(#[case] input: &str, #[case] expected: SGroupType) {
        let (remaining, result) = sgroup_type().parse(input.as_bytes()).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case("XYZ", "unknown type", ErrorKind::MapRes)]
    #[case("SU", "too short", ErrorKind::Eof)]
    fn test_sgroup_type_invalid(
        #[case] input: &str,
        #[case] desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let result = sgroup_type().parse(input.as_bytes());
        assert!(result.is_err());
        assert!(
            matches!(result.as_ref(), Err(nom::Err::Error(e)) if e.code == expected_kind),
            "Expected {:?} error for {}, got {:?}",
            expected_kind,
            desc,
            result
        );
    }

    #[rstest]
    #[case("ALT", SGroupSubtype::Alternating)]
    #[case("RAN", SGroupSubtype::Random)]
    #[case("BLO", SGroupSubtype::Block)]
    fn test_sgroup_subtype(#[case] input: &str, #[case] expected: SGroupSubtype) {
        let (remaining, result) = sgroup_subtype().parse(input.as_bytes()).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case("XYZ", "unknown subtype", ErrorKind::MapRes)]
    #[case("SU", "too short", ErrorKind::Eof)]
    fn test_sgroup_subtype_invalid(
        #[case] input: &str,
        #[case] desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let result = sgroup_subtype().parse(input.as_bytes());
        assert!(result.is_err());
        assert!(
            matches!(result.as_ref(), Err(nom::Err::Error(e)) if e.code == expected_kind),
            "Expected {:?} error for {}, got {:?}",
            expected_kind,
            desc,
            result
        );
    }

    #[rstest]
    #[case("HH", SGroupConnectivity::HeadToHead)]
    #[case("HT", SGroupConnectivity::HeadToTail)]
    #[case("EU", SGroupConnectivity::EitherUnknown)]
    fn test_sgroup_connectivity(#[case] input: &str, #[case] expected: SGroupConnectivity) {
        let (remaining, result) = sgroup_connectivity().parse(input.as_bytes()).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case("XY", "unknown connectivity", ErrorKind::MapRes)]
    #[case("H", "too short", ErrorKind::Eof)]
    fn test_sgroup_connectivity_invalid(
        #[case] input: &str,
        #[case] desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let result = sgroup_connectivity().parse(input.as_bytes());
        assert!(result.is_err());
        assert!(
            matches!(result.as_ref(), Err(nom::Err::Error(e)) if e.code == expected_kind),
            "Expected {:?} error for {}, got {:?}",
            expected_kind,
            desc,
            result
        );
    }

    #[rstest]
    #[case("N", SGroupMultiplier::N)]
    #[case("1", SGroupMultiplier::Count(1))]
    #[case("2", SGroupMultiplier::Count(2))]
    fn test_sgroup_multiplier(#[case] input: &str, #[case] expected: SGroupMultiplier) {
        let (remaining, result) = sgroup_multiplier().parse(input.as_bytes()).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case("X", "unknown multiplier", ErrorKind::MapRes)]
    fn test_sgroup_multiplier_invalid(#[case] input: &str, #[case] desc: &str, #[case] expected_kind: ErrorKind) {
        let result = sgroup_multiplier().parse(input.as_bytes());
        assert!(result.is_err());
        assert!(
            matches!(result.as_ref(), Err(nom::Err::Error(e)) if e.code == expected_kind),
            "Expected {:?} error for {}, got {:?}",
            expected_kind,
            desc,
            result
        );
    }
}
