//! Utilities for molecular file input/output.

use std::borrow::Cow;
use std::sync::LazyLock;

use regex::bytes::Regex;

/// Regex to match Unicode whitespace characters
static WHITESPACE_REGEX: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"[\p{White-Space}--\r\n]").unwrap());

/// Replace Unicode whitespace characters with ASCII spaces
pub(crate) fn normalize_whitespace(input: &[u8]) -> Cow<'_, [u8]> {
    WHITESPACE_REGEX.replace_all(input, &b" "[..])
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::blank(b" ", b" ")]
    #[case::newline(b"\n", b"\n")]
    #[case::crlf(b"\r\n", b"\r\n")]
    #[case::tab(b"\t", b" ")]
    #[case::letters(b"abcd", b"abcd")]
    #[case::unicode_letters(b"\xce\xb1", b"\xce\xb1")]
    #[case::nbsp(b"\xc2\xa0", b" ")]
    fn test_normalize_whitespace(#[case] input: &[u8], #[case] expected: &[u8]) {
        let normalized = normalize_whitespace(input).into_owned();
        assert_eq!(normalized, expected);
    }
}
