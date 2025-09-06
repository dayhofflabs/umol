//! Utilities for molecular file input/output.

use std::borrow::Cow;

use once_cell::sync::Lazy;
use regex::bytes::Regex;

/// Regex to match Unicode whitespace characters
static WHITESPACE_REGEX: Lazy<Regex> =
    Lazy::new(|| Regex::new(r"[\p{White_Space}--[\r\n]]").unwrap());

/// Replace Unicode whitespace characters with ASCII spaces
pub(crate) fn normalize_whitespace(input: &[u8]) -> Cow<'_, [u8]> {
    WHITESPACE_REGEX.replace_all(input, &b" "[..])
}
