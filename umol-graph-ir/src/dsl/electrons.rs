//! Electron-counts head DSL: the mandatory per-atom electron-count head at the
//! start of an aromatic-system or multicenter-bond string, shared by both.

use std::fmt;

use winnow::ascii::{dec_int, multispace0};
use winnow::combinator::{delimited, opt, separated};
use winnow::error::ErrMode;
use winnow::Parser;

use super::error::{PResult, ParseError};
use crate::ast::electrons::ElectronCountsAst;

/// The mandatory per-atom electron counts at the head of an aromatic-system /
/// multicenter-bond string: `*` (undetermined) or a non-empty `[n,n,…]` vector
/// (whitespace ignored). Leading and unprefixed, before any `#` predicate.
pub(crate) fn electron_counts(i: &mut &str) -> PResult<ElectronCountsAst> {
    if opt('*').parse_next(i)?.is_some() {
        return Ok(ElectronCountsAst::Undetermined);
    }
    if i.starts_with('[') {
        let start = *i;
        return electron_counts_vector(i)
            .map(ElectronCountsAst::Lit)
            .map_err(|_| ErrMode::Cut(ParseError::MalformedElectronCounts(start.to_string())));
    }
    Err(ErrMode::Cut(ParseError::ExpectedElectronCounts))
}

fn electron_counts_vector(i: &mut &str) -> PResult<Vec<i64>> {
    delimited(
        '[',
        delimited(
            multispace0,
            separated(
                1..,
                dec_int::<_, i64, _>,
                delimited(multispace0, ',', multispace0),
            ),
            multispace0,
        ),
        ']',
    )
    .parse_next(i)
}

pub(crate) fn fmt_electron_counts(
    f: &mut fmt::Formatter<'_>,
    electrons: &ElectronCountsAst,
) -> fmt::Result {
    match electrons {
        ElectronCountsAst::Undetermined => write!(f, "*"),
        ElectronCountsAst::Lit(counts) => {
            write!(f, "[")?;
            for (idx, n) in counts.iter().enumerate() {
                if idx > 0 {
                    write!(f, ",")?;
                }
                write!(f, "{}", n)?;
            }
            write!(f, "]")
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::undetermined("*", ElectronCountsAst::Undetermined)]
    #[case::single("[1]", ElectronCountsAst::Lit(vec![1]))]
    #[case::triple("[1,1,1]", ElectronCountsAst::Lit(vec![1, 1, 1]))]
    #[case::mixed("[2,0,2]", ElectronCountsAst::Lit(vec![2, 0, 2]))]
    #[case::whitespace("[ 1 , 0 , 2 ]", ElectronCountsAst::Lit(vec![1, 0, 2]))]
    #[case::negative("[-1,2]", ElectronCountsAst::Lit(vec![-1, 2]))]
    fn test_electron_counts(#[case] input: &str, #[case] expected: ElectronCountsAst) {
        assert_eq!(
            electron_counts.parse(input).map_err(|e| e.into_inner()),
            Ok(expected)
        );
    }

    #[rstest]
    #[case::missing("#c0", ParseError::ExpectedElectronCounts)]
    #[case::empty_input("", ParseError::ExpectedElectronCounts)]
    #[case::no_bracket("1,2", ParseError::ExpectedElectronCounts)]
    #[case::empty_vector("[]", ParseError::MalformedElectronCounts("[]".to_string()))]
    #[case::unmatched_open("[1,2", ParseError::MalformedElectronCounts("[1,2".to_string()))]
    #[case::non_numeric("[a]", ParseError::MalformedElectronCounts("[a]".to_string()))]
    #[case::trailing_comma("[1,2,]", ParseError::MalformedElectronCounts("[1,2,]".to_string()))]
    fn test_electron_counts_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(
            electron_counts.parse(input).map_err(|e| e.into_inner()),
            Err(expected)
        );
    }
}
