//! Utilities for molecular file input/output.

use std::borrow::Cow;
use std::sync::LazyLock;

use regex::bytes::Regex;

/// Regex to match Unicode whitespace characters
static WHITESPACE_REGEX: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"[\p{White_Space}--[\r\n]]").unwrap());

/// Replace Unicode whitespace characters with ASCII spaces
pub(crate) fn normalize_whitespace(input: &[u8]) -> Cow<'_, [u8]> {
    WHITESPACE_REGEX.replace_all(input, &b" "[..])
}
