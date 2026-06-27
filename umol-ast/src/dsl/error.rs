//! Domain errors for DSL parsing.

use thiserror::Error;
use winnow::error::{ErrMode, ParserError};
use winnow::stream::Stream;

pub(crate) type PResult<T> = Result<T, ErrMode<ParseError>>;

/// Error raised when a DSL input fails to parse (invalid syntax, unknown
/// predicate, duplicate predicate, unresolved ref, etc.).
#[rustfmt::skip]
#[derive(Clone, Debug, PartialEq, Error)]
pub enum ParseError {
    #[error("expected atom element")]
    ExpectedElement,
    #[error("expected predicate body")]
    ExpectedPredicateBody,
    #[error("unknown atom predicate: {0}")]
    UnknownAtomPredicate(String),
    #[error("duplicate atom predicate: {0}")]
    DuplicateAtomPredicate(String),
    #[error("unknown bond predicate: {0}")]
    UnknownBondPredicate(String),
    #[error("duplicate bond predicate: {0}")]
    DuplicateBondPredicate(String),
    #[error("unknown aromatic-system predicate: {0}")]
    UnknownAromaticSystemPredicate(String),
    #[error("duplicate aromatic-system predicate: {0}")]
    DuplicateAromaticSystemPredicate(String),
    #[error("unknown multicenter-bond predicate: {0}")]
    UnknownMulticenterBondPredicate(String),
    #[error("duplicate multicenter-bond predicate: {0}")]
    DuplicateMulticenterBondPredicate(String),
    #[error("unknown dative-bond predicate: {0}")]
    UnknownDativeBondPredicate(String),
    #[error("duplicate dative-bond predicate: {0}")]
    DuplicateDativeBondPredicate(String),
    #[error("unknown stereo predicate: {0}")]
    UnknownStereoPredicate(String),
    #[error("duplicate stereo predicate: {0}")]
    DuplicateStereoPredicate(String),
    #[error("expected noncovalent-bond kind")]
    ExpectedNoncovalentBondKind,
    #[error("expected electron counts ('*' or '[...]')")]
    ExpectedElectronCounts,
    #[error("malformed electron counts: {0:?}")]
    MalformedElectronCounts(String),
    #[error("trailing input: {0:?}")]
    TrailingInput(String),
    #[error("raising error: {0}")]
    RaisingError(String),
    #[error("lowering error: {0}")]
    LoweringError(String),
    #[error("syntax error")]
    Syntax,
    #[error("EDN parse: {0}")]
    EdnParse(String),
    #[error("missing key: {0}")]
    MissingKey(String),
    #[error("duplicate id: {0}")]
    DuplicateId(String),
    #[error("invalid value: {0}")]
    InvalidValue(String),
    #[error("invalid {kind} ref: {value}")]
    InvalidRef { kind: &'static str, value: String },
    #[error("{field}: expected {expected}")]
    WrongFieldType { field: String, expected: String },
}

impl<I: Stream> ParserError<I> for ParseError {
    type Inner = Self;

    fn from_input(_input: &I) -> Self {
        ParseError::Syntax
    }

    fn into_inner(self) -> Result<Self::Inner, Self> {
        Ok(self)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    fn test_parse_error_parser_error_from_input_is_syntax() {
        let input: &[u8] = b"bogus";
        let err = <ParseError as ParserError<&[u8]>>::from_input(&input);
        assert_eq!(err, ParseError::Syntax);
    }

    #[rstest]
    #[case::trailing(ParseError::TrailingInput("rest".into()))]
    #[case::syntax(ParseError::Syntax)]
    #[case::invalid_value(ParseError::InvalidValue("x".into()))]
    fn test_parse_error_into_inner_is_identity(#[case] err: ParseError) {
        let cloned = err.clone();
        let result = <ParseError as ParserError<&[u8]>>::into_inner(err);
        assert_eq!(result, Ok(cloned));
    }

    #[rstest]
    #[case::expected_element(ParseError::ExpectedElement, "expected atom element")]
    #[case::unknown_atom_predicate(
        ParseError::UnknownAtomPredicate("foo".into()),
        "unknown atom predicate: foo",
    )]
    #[case::invalid_ref(
        ParseError::InvalidRef { kind: "atom", value: "7".into() },
        "invalid atom ref: 7",
    )]
    #[case::wrong_field_type(
        ParseError::WrongFieldType { field: "charge".into(), expected: "int".into() },
        "charge: expected int",
    )]
    fn test_parse_error_display(#[case] err: ParseError, #[case] expected: &str) {
        assert_eq!(err.to_string(), expected);
    }
}
