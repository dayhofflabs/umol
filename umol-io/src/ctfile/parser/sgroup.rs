//! Auxiliary parsers for SGroup properties.

use bstr::ByteSlice;
use winnow::ascii::space0;
use winnow::combinator::alt;
use winnow::stream::{Location, Stream};
use winnow::token::{rest, take, take_while};
use winnow::Parser;

use crate::ctfile::parser::utils::{fixed_width_partial, Input, InputError, IntParser, PResult};
use crate::table_ir::{
    SGroupConnectivity, SGroupDataDisplayChars, SGroupDataDisplayPlacement, SGroupDataDisplayType,
    SGroupDataDisplayUnits, SGroupDataType, SGroupMultiplier, SGroupMultiplierOp,
    SGroupMultiplierTerm, SGroupSubtype, SGroupType,
};

/// Parse SGroup type string.
pub(super) fn sgroup_type(input: &mut Input<'_>) -> PResult<SGroupType> {
    take(3usize)
        .verify_map(|s: &[u8]| match s {
            b"SUP" => Some(SGroupType::Superatom),
            b"MUL" => Some(SGroupType::MultipleGroup),
            b"SRU" => Some(SGroupType::RepeatingUnit),
            b"MON" => Some(SGroupType::Monomer),
            b"MER" => Some(SGroupType::Mer),
            b"COP" => Some(SGroupType::Copolymer),
            b"CRO" => Some(SGroupType::Crosslink),
            b"MOD" => Some(SGroupType::Modification),
            b"GRA" => Some(SGroupType::Graft),
            b"COM" => Some(SGroupType::Component),
            b"MIX" => Some(SGroupType::Mixture),
            b"FOR" => Some(SGroupType::Formulation),
            b"DAT" => Some(SGroupType::Data),
            b"ANY" => Some(SGroupType::AnyPolymer),
            b"GEN" => Some(SGroupType::Generic),
            _ => None,
        })
        .parse_next(input)
}

/// Parse SGroup subtype string.
pub(super) fn sgroup_subtype(input: &mut Input<'_>) -> PResult<SGroupSubtype> {
    take(3usize)
        .verify_map(|s: &[u8]| match s {
            b"ALT" => Some(SGroupSubtype::Alternating),
            b"RAN" => Some(SGroupSubtype::Random),
            b"BLO" => Some(SGroupSubtype::Block),
            _ => None,
        })
        .parse_next(input)
}

/// Parse SGroup connectivity string.
pub(super) fn sgroup_connectivity(input: &mut Input<'_>) -> PResult<SGroupConnectivity> {
    fixed_width_partial(
        3usize,
        rest.verify_map(|s: &[u8]| match s.trim_ascii() {
            b"HH" => Some(SGroupConnectivity::HeadToHead),
            b"HT" => Some(SGroupConnectivity::HeadToTail),
            b"EU" => Some(SGroupConnectivity::EitherUnknown),
            _ => None,
        }),
        true,
    )
    .verify_map(|value| value)
    .parse_next(input)
}

/// Parse SGroup multiplier string.
pub(super) fn sgroup_multiplier(input: &mut Input<'_>) -> PResult<SGroupMultiplier> {
    alt((
        (sgroup_multiplier_term, space0, sgroup_multiplier_variable).map(|(left, _, right)| {
            SGroupMultiplier::Expression {
                left,
                op: SGroupMultiplierOp::Mul,
                right,
            }
        }),
        (
            sgroup_multiplier_term,
            space0,
            sgroup_multiplier_op,
            space0,
            sgroup_multiplier_term,
        )
            .map(|(left, _, op, _, right)| SGroupMultiplier::Expression {
                left,
                op,
                right,
            }),
        sgroup_multiplier_term.map(SGroupMultiplier::Single),
    ))
    .parse_next(input)
}

/// Parse an integer multiplier.
fn sgroup_multiplier_integer(input: &mut Input<'_>) -> PResult<SGroupMultiplierTerm> {
    <u32 as IntParser>::parse
        .map(SGroupMultiplierTerm::Integer)
        .parse_next(input)
}

/// Parse a single-character variable multiplier.
fn sgroup_multiplier_variable(input: &mut Input<'_>) -> PResult<SGroupMultiplierTerm> {
    take_while(1..=1, |byte: u8| byte.is_ascii_alphabetic())
        .map(|s: &[u8]| SGroupMultiplierTerm::Variable(s[0] as char))
        .parse_next(input)
}

/// Parse a single multiplier term (variable or integer).
fn sgroup_multiplier_term(input: &mut Input<'_>) -> PResult<SGroupMultiplierTerm> {
    alt((sgroup_multiplier_integer, sgroup_multiplier_variable)).parse_next(input)
}

/// Parse arithmetic operator.
fn sgroup_multiplier_op(input: &mut Input<'_>) -> PResult<SGroupMultiplierOp> {
    alt((
        b'+'.value(SGroupMultiplierOp::Add),
        b'-'.value(SGroupMultiplierOp::Sub),
        b'*'.value(SGroupMultiplierOp::Mul),
        b'/'.value(SGroupMultiplierOp::Div),
    ))
    .parse_next(input)
}

/// Parse SGroup subscript string.
pub(super) fn sgroup_subscript(input: &mut Input<'_>) -> PResult<String> {
    rest.map(|s: &[u8]| s.to_str_lossy().into_owned())
        .parse_next(input)
}

// Parse SGroup data type string.
pub(super) fn sgroup_data_type(input: &mut Input<'_>) -> PResult<SGroupDataType> {
    fixed_width_partial(2usize, rest, true)
        .verify_map(|value: Option<&[u8]>| match value.map(<[u8]>::trim_ascii) {
            None | Some(b"") | Some(b"T") => Some(SGroupDataType::Text),
            Some(b"F") => Some(SGroupDataType::Formatted),
            Some(b"N") => Some(SGroupDataType::Numeric),
            _ => None,
        })
        .parse_next(input)
}

// Parse SGroup data display type string.
pub(super) fn sgroup_data_display_type(input: &mut Input<'_>) -> PResult<SGroupDataDisplayType> {
    take(1usize)
        .verify_map(|s: &[u8]| match s {
            b"A" => Some(SGroupDataDisplayType::Attached),
            b"D" => Some(SGroupDataDisplayType::Detached),
            _ => None,
        })
        .parse_next(input)
}

// Parse SGroup data display placement string.
pub(super) fn sgroup_data_display_placement(
    input: &mut Input<'_>,
) -> PResult<SGroupDataDisplayPlacement> {
    take(1usize)
        .verify_map(|s: &[u8]| match s {
            b"A" => Some(SGroupDataDisplayPlacement::Absolute),
            b"R" => Some(SGroupDataDisplayPlacement::Relative),
            _ => None,
        })
        .parse_next(input)
}

// Parse SGroup data display units string.
pub(super) fn sgroup_data_display_units(input: &mut Input<'_>) -> PResult<SGroupDataDisplayUnits> {
    take(1usize)
        .verify_map(|s: &[u8]| match s {
            b" " => Some(SGroupDataDisplayUnits::None),
            b"U" => Some(SGroupDataDisplayUnits::DisplayUnits),
            _ => None,
        })
        .parse_next(input)
}

// Parse SGroup data display chars string.
pub(super) fn sgroup_data_display_chars(input: &mut Input<'_>) -> PResult<SGroupDataDisplayChars> {
    let start = input.checkpoint();
    let column = input.current_token_start();
    let field: &[u8] = take(3usize).parse_next(input)?;
    let trimmed = field.trim_ascii();
    if trimmed.eq_ignore_ascii_case(b"ALL") {
        return Ok(SGroupDataDisplayChars::All);
    }

    let mut field_input = Input::new(trimmed);
    match <u32 as IntParser>::parse(&mut field_input) {
        Ok(value) if field_input.is_empty() => Ok(SGroupDataDisplayChars::Number(value)),
        Ok(_) | Err(_) => {
            input.reset(&start);
            Err(winnow::error::ErrMode::Backtrack(InputError::at_column(
                column,
            )))
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use winnow::error::ErrMode;

    use super::*;

    fn expected_error(column: u32) -> ErrMode<InputError> {
        ErrMode::Backtrack(InputError { column })
    }

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
        assert_eq!(sgroup_type.parse(Input::new(input)), Ok(expected));
    }

    #[rstest]
    #[case::unknown(b"XYZ")]
    #[case::short(b"SU")]
    fn test_sgroup_type_error(#[case] input: &[u8]) {
        assert_eq!(
            sgroup_type.parse_peek(Input::new(input)),
            Err(expected_error(0))
        );
    }

    #[rstest]
    #[case(b"ALT", SGroupSubtype::Alternating)]
    #[case(b"RAN", SGroupSubtype::Random)]
    #[case(b"BLO", SGroupSubtype::Block)]
    fn test_sgroup_subtype(#[case] input: &[u8], #[case] expected: SGroupSubtype) {
        assert_eq!(sgroup_subtype.parse(Input::new(input)), Ok(expected));
    }

    #[rstest]
    #[case::unknown(b"XYZ")]
    #[case::short(b"SU")]
    fn test_sgroup_subtype_error(#[case] input: &[u8]) {
        assert_eq!(
            sgroup_subtype.parse_peek(Input::new(input)),
            Err(expected_error(0))
        );
    }

    #[rstest]
    #[case(b"HH ", SGroupConnectivity::HeadToHead)]
    #[case(b"HT ", SGroupConnectivity::HeadToTail)]
    #[case(b"EU ", SGroupConnectivity::EitherUnknown)]
    #[case(b" HH", SGroupConnectivity::HeadToHead)]
    fn test_sgroup_connectivity(#[case] input: &[u8], #[case] expected: SGroupConnectivity) {
        assert_eq!(sgroup_connectivity.parse(Input::new(input)), Ok(expected));
    }

    #[rstest]
    #[case(b"XY ")]
    fn test_sgroup_connectivity_error(#[case] input: &[u8]) {
        assert_eq!(
            sgroup_connectivity.parse_peek(Input::new(input)),
            Err(expected_error(0))
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
        assert_eq!(sgroup_multiplier.parse(Input::new(input)), Ok(expected));
    }

    #[rstest]
    #[case(b"@")]
    fn test_sgroup_multiplier_error(#[case] input: &[u8]) {
        assert_eq!(
            sgroup_multiplier.parse_peek(Input::new(input)),
            Err(expected_error(0))
        );
    }

    #[rstest]
    #[case(b"F ", SGroupDataType::Formatted)]
    #[case(b"N ", SGroupDataType::Numeric)]
    #[case(b"T ", SGroupDataType::Text)]
    #[case(b"  ", SGroupDataType::Text)]
    fn test_sgroup_data_type(#[case] input: &[u8], #[case] expected: SGroupDataType) {
        assert_eq!(sgroup_data_type.parse(Input::new(input)), Ok(expected));
    }

    #[rstest]
    #[case(b"X")]
    fn test_sgroup_data_type_error(#[case] input: &[u8]) {
        assert_eq!(
            sgroup_data_type.parse_peek(Input::new(input)),
            Err(expected_error(0))
        );
    }

    #[rstest]
    #[case(b"A", SGroupDataDisplayType::Attached)]
    #[case(b"D", SGroupDataDisplayType::Detached)]
    fn test_sgroup_data_display_type(
        #[case] input: &[u8],
        #[case] expected: SGroupDataDisplayType,
    ) {
        assert_eq!(
            sgroup_data_display_type.parse(Input::new(input)),
            Ok(expected)
        );
    }

    #[rstest]
    #[case(b"X")]
    fn test_sgroup_data_display_type_error(#[case] input: &[u8]) {
        assert_eq!(
            sgroup_data_display_type.parse_peek(Input::new(input)),
            Err(expected_error(0))
        );
    }

    #[rstest]
    #[case(b"A", SGroupDataDisplayPlacement::Absolute)]
    #[case(b"R", SGroupDataDisplayPlacement::Relative)]
    fn test_sgroup_data_display_placement(
        #[case] input: &[u8],
        #[case] expected: SGroupDataDisplayPlacement,
    ) {
        assert_eq!(
            sgroup_data_display_placement.parse(Input::new(input)),
            Ok(expected)
        );
    }

    #[rstest]
    #[case(b"X")]
    fn test_sgroup_data_display_placement_error(#[case] input: &[u8]) {
        assert_eq!(
            sgroup_data_display_placement.parse_peek(Input::new(input)),
            Err(expected_error(0))
        );
    }

    #[rstest]
    #[case(b" ", SGroupDataDisplayUnits::None)]
    #[case(b"U", SGroupDataDisplayUnits::DisplayUnits)]
    fn test_sgroup_data_display_units(
        #[case] input: &[u8],
        #[case] expected: SGroupDataDisplayUnits,
    ) {
        assert_eq!(
            sgroup_data_display_units.parse(Input::new(input)),
            Ok(expected)
        );
    }

    #[rstest]
    #[case(b"X")]
    fn test_sgroup_data_display_units_error(#[case] input: &[u8]) {
        assert_eq!(
            sgroup_data_display_units.parse_peek(Input::new(input)),
            Err(expected_error(0))
        );
    }

    #[rstest]
    #[case(b"ALL", SGroupDataDisplayChars::All)]
    #[case(b"  1", SGroupDataDisplayChars::Number(1))]
    fn test_sgroup_data_display_chars(
        #[case] input: &[u8],
        #[case] expected: SGroupDataDisplayChars,
    ) {
        assert_eq!(
            sgroup_data_display_chars.parse(Input::new(input)),
            Ok(expected)
        );
    }

    #[rstest]
    #[case(b"X")]
    fn test_sgroup_data_display_chars_error(#[case] input: &[u8]) {
        assert_eq!(
            sgroup_data_display_chars.parse_peek(Input::new(input)),
            Err(expected_error(0))
        );
    }
}
