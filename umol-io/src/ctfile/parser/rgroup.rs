//! RGroup parser for CTab files.

use winnow::combinator::{alt, opt};
use winnow::Parser;

use super::utils::{Input, IntParser, PResult};
use crate::table_ir::RGroup;

/// Parse RGroup symbol from byte slice.
///
/// Handles:
/// - "R" / "R#" / "R0" -> No label
/// - "R1", "R2", etc. -> Label n (n > 0)
pub(super) fn rgroup_symbol(input: &mut Input<'_>) -> PResult<RGroup> {
    b'R'.parse_next(input)?;
    let label = opt(alt((
        b'#'.value(None),
        <u32 as IntParser>::parse.map(|value| (value != 0).then_some(value)),
    )))
    .parse_next(input)?;
    Ok(RGroup::new(label.flatten()))
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use winnow::Parser;

    use super::*;

    #[rstest]
    #[case::unlabelled(b"R", RGroup::new(None))]
    #[case::hash(b"R#", RGroup::new(None))]
    #[case::zero(b"R0", RGroup::new(None))]
    #[case::label(b"R1", RGroup::new(Some(1)))]
    #[case::multiple_digits(b"R12", RGroup::new(Some(12)))]
    #[case::zero_padded(b"R012", RGroup::new(Some(12)))]
    fn test_rgroup_symbol(#[case] input: &[u8], #[case] expected: RGroup) {
        assert_eq!(rgroup_symbol.parse(Input::new(input)), Ok(expected));
    }
}
