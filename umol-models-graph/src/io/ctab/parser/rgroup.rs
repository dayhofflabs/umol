//! RGroup parser for CTab files.

use nom::branch::alt;
use nom::bytes::complete::tag;
use nom::character::complete::u32 as nom_u32;
use nom::combinator::{map, opt, value};
use nom::{error, IResult, Parser};

use crate::table_ir::RGroup;

/// Parse RGroup symbol from byte slice.
///
/// Handles:
/// - "R" / "R#" / "R0" -> No label
/// - "R1", "R2", etc. -> Label n (n > 0)
pub(super) fn rgroup_symbol(input: &[u8]) -> IResult<&[u8], RGroup, error::Error<&[u8]>> {
    let (remaining, _) = tag("R").parse(input)?;
    let (remaining, label) = opt(alt((
        value(None, tag("#")),
        map(nom_u32, |n| if n == 0 { None } else { Some(n) }),
    )))
    .parse(remaining)?;
    Ok((remaining, RGroup::new(label.flatten())))
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
    fn test_rgroup_from_symbol_bytes(#[case] input: &[u8], #[case] expected: RGroup) {
        let result = rgroup_symbol(input);
        assert!(result.is_ok());
        let (remaining, symbol) = result.unwrap();
        assert!(remaining.is_empty());
        assert_eq!(symbol, expected);
    }
}
