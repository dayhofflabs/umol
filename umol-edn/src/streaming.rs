//! Streaming EDN reader.
//!
//! `EdnStreamDeserializer` exposes byte-level primitives for fused
//! `FromEdn::from_edn_str` overrides. It is a low-level primitive used by
//! the `umol-edn-macros` derive output and by hand-written hot-path
//! parsers; it is not a serde Deserializer.

use std::borrow::Cow;
use std::fmt;
use std::str::{from_utf8, FromStr};

use crate::error::{DeError, EdnError, ParseError};
use crate::parser::{is_symbol_char, is_symbol_start, validate_symbol};

pub struct EdnStreamDeserializer<'de> {
    input: &'de str,
    pos: usize,
    scratch: Vec<u8>,
}

impl<'de> EdnStreamDeserializer<'de> {
    pub fn new(input: &'de str) -> Self {
        Self {
            input,
            pos: 0,
            scratch: Vec::new(),
        }
    }

    pub fn expect_eof(&mut self) -> Result<(), ParseError> {
        self.skip_ws()?;
        if self.pos < self.input.len() {
            Err(ParseError::TrailingContent { offset: self.pos })
        } else {
            Ok(())
        }
    }

    #[inline]
    fn peek(&self) -> Option<u8> {
        self.input.as_bytes().get(self.pos).copied()
    }

    fn next_byte(&mut self) -> Result<u8, ParseError> {
        match self.input.as_bytes().get(self.pos) {
            Some(&b) => {
                self.pos += 1;
                Ok(b)
            }
            None => Err(ParseError::UnexpectedEof { offset: self.pos }),
        }
    }

    #[inline]
    fn skip_ws(&mut self) -> Result<(), ParseError> {
        let bytes = self.input.as_bytes();
        while self.pos < bytes.len() {
            match bytes[self.pos] {
                b' ' | b'\t' | b'\n' | b'\r' | b',' => self.pos += 1,
                b';' => {
                    self.pos += 1;
                    while self.pos < bytes.len() && bytes[self.pos] != b'\n' {
                        self.pos += 1;
                    }
                }
                b'#' if bytes.get(self.pos + 1) == Some(&b'_') => {
                    return self.skip_ws_discard();
                }
                _ => break,
            }
        }
        Ok(())
    }

    #[cold]
    fn skip_ws_discard(&mut self) -> Result<(), ParseError> {
        self.pos += 2;
        self.skip_ws()?;
        self.skip_value()?;
        self.skip_ws()
    }

    #[inline]
    fn parse_keyword_name(&mut self) -> Result<Cow<'de, str>, ParseError> {
        debug_assert_eq!(self.input.as_bytes()[self.pos], b':');
        self.pos += 1;

        // :/ and :/foo are not legal keywords per the EDN spec.
        if self.input.as_bytes().get(self.pos) == Some(&b'/') {
            return Err(ParseError::UnexpectedToken {
                offset: self.pos,
                found: '/',
            });
        }
        self.parse_symbol_str().map(Cow::Borrowed)
    }

    #[inline]
    fn parse_symbol_str(&mut self) -> Result<&'de str, ParseError> {
        let start = self.pos;
        let bytes = self.input.as_bytes();
        if self.pos >= bytes.len() {
            return Err(ParseError::UnexpectedEof { offset: start });
        }
        let first = bytes[self.pos] as char;
        if !is_symbol_start(first) && !first.is_ascii_digit() {
            return Err(ParseError::UnexpectedToken {
                offset: start,
                found: first,
            });
        }
        self.pos += 1;
        while self.pos < bytes.len() && is_symbol_char(bytes[self.pos] as char) {
            self.pos += 1;
        }
        let s = &self.input[start..self.pos];
        validate_symbol(s, start)?;
        Ok(s)
    }

    fn parse_number_i64(&mut self) -> Result<i64, ParseError> {
        let s = self.scan_number_str()?;
        s.parse::<i64>().map_err(|_| ParseError::InvalidNumber {
            offset: self.pos - s.len(),
        })
    }

    fn scan_number_str(&mut self) -> Result<&'de str, ParseError> {
        let start = self.pos;
        let bytes = self.input.as_bytes();
        // optional sign
        if self.pos < bytes.len() && matches!(bytes[self.pos], b'+' | b'-') {
            self.pos += 1;
        }
        // digits
        if self.pos >= bytes.len() {
            return Err(ParseError::UnexpectedEof { offset: start });
        }
        if !bytes[self.pos].is_ascii_digit() {
            return Err(ParseError::InvalidNumber { offset: start });
        }
        let digit_start = self.pos;
        while self.pos < bytes.len() && bytes[self.pos].is_ascii_digit() {
            self.pos += 1;
        }
        let digit_len = self.pos - digit_start;
        let has_dot = self.pos < bytes.len() && bytes[self.pos] == b'.';
        if digit_len > 1 && bytes[digit_start] == b'0' && !has_dot {
            return Err(ParseError::InvalidNumber { offset: start });
        }
        // optional .digits
        if has_dot {
            self.pos += 1;
            while self.pos < bytes.len() && bytes[self.pos].is_ascii_digit() {
                self.pos += 1;
            }
        }
        // optional exponent
        if self.pos < bytes.len() && matches!(bytes[self.pos], b'e' | b'E') {
            self.pos += 1;
            if self.pos < bytes.len() && matches!(bytes[self.pos], b'+' | b'-') {
                self.pos += 1;
            }
            while self.pos < bytes.len() && bytes[self.pos].is_ascii_digit() {
                self.pos += 1;
            }
        }
        // bignum suffix — reject without feature, leave for caller with feature
        if self.pos < bytes.len() && matches!(bytes[self.pos], b'N' | b'M') {
            #[cfg(not(feature = "bignum"))]
            return Err(ParseError::UnsupportedFeature {
                offset: start,
                feature: "bignum",
            });
        }
        Ok(&self.input[start..self.pos])
    }

    fn parse_string(&mut self) -> Result<Cow<'de, str>, ParseError> {
        debug_assert_eq!(self.input.as_bytes()[self.pos], b'"');
        self.pos += 1;
        let start = self.pos;
        let bytes = self.input.as_bytes();

        // Fast path: scan for closing quote with no escapes.
        if let Some(end) = memchr::memchr2(b'"', b'\\', &bytes[self.pos..]) {
            if bytes[self.pos + end] == b'"' {
                let s = &self.input[self.pos..self.pos + end];
                self.pos += end + 1;
                return Ok(Cow::Borrowed(s));
            }
        }

        // Slow path: has escapes.
        self.scratch.clear();
        // Copy bytes before the first backslash.
        if let Some(bs) = memchr::memchr(b'\\', &bytes[self.pos..]) {
            self.scratch
                .extend_from_slice(&bytes[self.pos..self.pos + bs]);
            self.pos += bs;
        }
        loop {
            if self.pos >= bytes.len() {
                return Err(ParseError::UnexpectedEof { offset: start });
            }
            match bytes[self.pos] {
                b'"' => {
                    self.pos += 1;
                    let s = from_utf8(&self.scratch)
                        .map_err(|_| ParseError::InvalidUtf8 { offset: start })?;
                    return Ok(Cow::Owned(s.to_string()));
                }
                b'\\' => {
                    let esc_offset = self.pos;
                    self.pos += 1;
                    let esc = self
                        .next_byte()
                        .map_err(|_| ParseError::InvalidEscape { offset: esc_offset })?;
                    match esc {
                        b't' => self.scratch.push(b'\t'),
                        b'r' => self.scratch.push(b'\r'),
                        b'n' => self.scratch.push(b'\n'),
                        b'\\' => self.scratch.push(b'\\'),
                        b'"' => self.scratch.push(b'"'),
                        b'u' => {
                            if self.pos + 4 > bytes.len() {
                                return Err(ParseError::InvalidEscape { offset: esc_offset });
                            }
                            let hex = from_utf8(&bytes[self.pos..self.pos + 4])
                                .map_err(|_| ParseError::InvalidEscape { offset: esc_offset })?;
                            let cp = u32::from_str_radix(hex, 16)
                                .map_err(|_| ParseError::InvalidEscape { offset: esc_offset })?;
                            let ch = char::from_u32(cp)
                                .ok_or(ParseError::InvalidEscape { offset: esc_offset })?;
                            self.pos += 4;
                            let mut buf = [0u8; 4];
                            let encoded = ch.encode_utf8(&mut buf);
                            self.scratch.extend_from_slice(encoded.as_bytes());
                        }
                        _ => return Err(ParseError::InvalidEscape { offset: esc_offset }),
                    }
                    // Batch copy until next " or \.
                    if let Some(span) = memchr::memchr2(b'"', b'\\', &bytes[self.pos..]) {
                        self.scratch
                            .extend_from_slice(&bytes[self.pos..self.pos + span]);
                        self.pos += span;
                    } else {
                        self.scratch.extend_from_slice(&bytes[self.pos..]);
                        self.pos = bytes.len();
                    }
                }
                other => {
                    self.scratch.push(other);
                    self.pos += 1;
                }
            }
        }
    }

    fn skip_value(&mut self) -> Result<(), ParseError> {
        self.skip_ws()?;
        let b = self
            .peek()
            .ok_or(ParseError::UnexpectedEof { offset: self.pos })?;
        match b {
            b'(' => {
                self.pos += 1;
                self.skip_delimited(b')')
            }
            b'[' => {
                self.pos += 1;
                self.skip_delimited(b']')
            }
            b'{' => {
                self.pos += 1;
                self.skip_delimited(b'}')
            }
            b'#' => {
                self.pos += 1;
                match self.peek() {
                    Some(b'{') => {
                        self.pos += 1;
                        self.skip_delimited(b'}')
                    }
                    Some(b'_') => {
                        self.pos += 1;
                        self.skip_ws()?;
                        self.skip_value()?;
                        self.skip_value()
                    }
                    _ => {
                        self.skip_atom()?;
                        self.skip_ws()?;
                        self.skip_value()
                    }
                }
            }
            b'"' => self.skip_string(),
            b'\\' => {
                self.pos += 1;
                self.skip_atom()
            }
            b':' => {
                self.pos += 1;
                self.skip_atom()
            }
            _ => self.skip_atom(),
        }
    }

    fn skip_delimited(&mut self, close: u8) -> Result<(), ParseError> {
        loop {
            self.skip_ws()?;
            match self.peek() {
                None => return Err(ParseError::UnexpectedEof { offset: self.pos }),
                Some(b) if b == close => {
                    self.pos += 1;
                    return Ok(());
                }
                _ => self.skip_value()?,
            }
        }
    }

    fn skip_string(&mut self) -> Result<(), ParseError> {
        debug_assert_eq!(self.input.as_bytes()[self.pos], b'"');
        self.pos += 1;
        let bytes = self.input.as_bytes();
        loop {
            match memchr::memchr2(b'"', b'\\', &bytes[self.pos..]) {
                None => return Err(ParseError::UnexpectedEof { offset: self.pos }),
                Some(i) => {
                    self.pos += i;
                    if bytes[self.pos] == b'"' {
                        self.pos += 1;
                        return Ok(());
                    }
                    // Skip escape: backslash + escaped content
                    self.pos += 1; // skip backslash
                    if self.pos >= bytes.len() {
                        return Err(ParseError::UnexpectedEof { offset: self.pos });
                    }
                    if bytes[self.pos] == b'u' {
                        // \uXXXX — skip u + 4 hex digits
                        self.pos = (self.pos + 5).min(bytes.len());
                    } else {
                        self.pos += 1;
                    }
                }
            }
        }
    }

    fn skip_atom(&mut self) -> Result<(), ParseError> {
        let bytes = self.input.as_bytes();
        if self.pos >= bytes.len() {
            return Err(ParseError::UnexpectedEof { offset: self.pos });
        }
        let start = self.pos;
        while self.pos < bytes.len() && is_symbol_char(bytes[self.pos] as char) {
            self.pos += 1;
        }
        if self.pos == start {
            return Err(ParseError::UnexpectedToken {
                offset: self.pos,
                found: bytes[self.pos] as char,
            });
        }
        Ok(())
    }
}

impl<'de> EdnStreamDeserializer<'de> {
    /// Current byte offset into the input.
    #[inline]
    pub fn position(&self) -> usize {
        self.pos
    }

    /// Skip whitespace and return the next byte without consuming.
    #[inline]
    pub fn peek_byte(&mut self) -> Result<Option<u8>, ParseError> {
        self.skip_ws()?;
        Ok(self.peek())
    }

    /// Skip whitespace and consume `expected`. Errors if it does not match.
    pub fn consume_byte(&mut self, expected: u8) -> Result<(), ParseError> {
        self.skip_ws()?;
        match self.peek() {
            Some(b) if b == expected => {
                self.pos += 1;
                Ok(())
            }
            other => Err(unexpected_or_eof(self.pos, other)),
        }
    }

    /// Skip whitespace; consume `expected` if present, returning whether it matched.
    pub fn try_consume_byte(&mut self, expected: u8) -> Result<bool, ParseError> {
        self.skip_ws()?;
        if self.peek() == Some(expected) {
            self.pos += 1;
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// Read a keyword (`:foo`) at the current position; returns its name without `:`.
    pub fn read_keyword_name(&mut self) -> Result<Cow<'de, str>, ParseError> {
        self.skip_ws()?;
        if self.peek() != Some(b':') {
            return Err(unexpected_or_eof(self.pos, self.peek()));
        }
        self.parse_keyword_name()
    }

    /// Read a string literal at the current position.
    pub fn read_string(&mut self) -> Result<Cow<'de, str>, ParseError> {
        self.skip_ws()?;
        if self.peek() != Some(b'"') {
            return Err(unexpected_or_eof(self.pos, self.peek()));
        }
        self.parse_string()
    }

    /// Read either a string literal or a keyword name (DSL convention for atom names).
    pub fn read_string_or_keyword(&mut self) -> Result<Cow<'de, str>, ParseError> {
        self.skip_ws()?;
        match self.peek() {
            Some(b'"') => self.parse_string(),
            Some(b':') => self.parse_keyword_name(),
            other => Err(unexpected_or_eof(self.pos, other)),
        }
    }

    /// Read a signed 64-bit integer.
    pub fn read_i64(&mut self) -> Result<i64, ParseError> {
        self.skip_ws()?;
        self.parse_number_i64()
    }

    /// Skip an arbitrary EDN value at the current position.
    pub fn read_skip_value(&mut self) -> Result<(), ParseError> {
        self.skip_value()
    }

    /// Return the source slice spanning one EDN value at the current position.
    /// Whitespace is skipped before measuring; the returned slice contains the
    /// raw EDN form of one value (atom, list, vector, map, set, or tagged).
    pub fn read_value_slice(&mut self) -> Result<&'de str, ParseError> {
        self.skip_ws()?;
        let start = self.pos;
        self.skip_value()?;
        Ok(&self.input[start..self.pos])
    }

    /// Read a string token and parse its contents via `FromStr`.
    ///
    /// On failure, wraps the subgrammar error in `DeError::Subgrammar` with
    /// the byte offset of the string token in the outer EDN source.
    pub fn read_subgrammar<T: FromStr>(&mut self, grammar: &'static str) -> Result<T, EdnError>
    where
        T::Err: fmt::Display,
    {
        let offset = self.position();
        let s = self.read_string()?;
        s.parse::<T>().map_err(|e| {
            DeError::Subgrammar {
                grammar,
                message: e.to_string(),
                path: vec![format!("@{offset}")],
            }
            .into()
        })
    }

    /// Like [`read_subgrammar`], but requires the string to be the only
    /// value in the input: reads the string token, parses its contents via
    /// `FromStr`, and then requires EOF. Use this to implement
    /// `FromEdn::from_edn_str` on a wrapper type whose EDN form is a single
    /// string literal handed to another grammar.
    ///
    /// [`read_subgrammar`]: Self::read_subgrammar
    pub fn read_subgrammar_all<T: FromStr>(&mut self, grammar: &'static str) -> Result<T, EdnError>
    where
        T::Err: fmt::Display,
    {
        let out = self.read_subgrammar::<T>(grammar)?;
        self.expect_eof()?;
        Ok(out)
    }
}

#[inline]
fn unexpected_or_eof(offset: usize, b: Option<u8>) -> ParseError {
    match b {
        Some(b) => ParseError::UnexpectedToken {
            offset,
            found: b as char,
        },
        None => ParseError::UnexpectedEof { offset },
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[test]
    fn test_edn_stream_deserializer_position() {
        let mut d = EdnStreamDeserializer::new("abc");
        assert_eq!(d.position(), 0);
        d.consume_byte(b'a').unwrap();
        assert_eq!(d.position(), 1);
    }

    #[rstest]
    #[case::empty("", true)]
    #[case::whitespace_only("  \t\n", true)]
    #[case::comment_only("; comment\n", true)]
    #[case::discard_form("#_ 123", true)]
    #[case::trailing_atom("x", false)]
    #[case::trailing_after_ws("  x", false)]
    fn test_edn_stream_deserializer_expect_eof(#[case] input: &str, #[case] ok: bool) {
        let mut d = EdnStreamDeserializer::new(input);
        assert_eq!(d.expect_eof().is_ok(), ok);
    }

    #[rstest]
    #[case::first_byte("abc", Some(b'a'))]
    #[case::skips_spaces("  x", Some(b'x'))]
    #[case::skips_commas(",, y", Some(b'y'))]
    #[case::skips_comment("; comment\nz", Some(b'z'))]
    #[case::empty_none("", None)]
    #[case::ws_only_none("  ", None)]
    fn test_edn_stream_deserializer_peek_byte(#[case] input: &str, #[case] expected: Option<u8>) {
        let mut d = EdnStreamDeserializer::new(input);
        assert_eq!(d.peek_byte().unwrap(), expected);
    }

    #[test]
    fn test_edn_stream_deserializer_consume_byte() {
        let mut d = EdnStreamDeserializer::new("  [1]");
        d.consume_byte(b'[').unwrap();
        assert_eq!(d.position(), 3);
    }

    #[rstest]
    #[case::wrong_byte("x")]
    #[case::eof("")]
    fn test_edn_stream_deserializer_consume_byte_error(#[case] input: &str) {
        let mut d = EdnStreamDeserializer::new(input);
        assert!(d.consume_byte(b'[').is_err());
    }

    #[rstest]
    #[case::matches("[x", b'[', true, 1)]
    #[case::no_match("x", b'[', false, 0)]
    #[case::skips_ws("  ]", b']', true, 3)]
    fn test_edn_stream_deserializer_try_consume_byte(
        #[case] input: &str,
        #[case] expected: u8,
        #[case] matched: bool,
        #[case] pos: usize,
    ) {
        let mut d = EdnStreamDeserializer::new(input);
        assert_eq!(d.try_consume_byte(expected).unwrap(), matched);
        assert_eq!(d.position(), pos);
    }

    #[rstest]
    #[case::simple(":foo", "foo")]
    #[case::with_ws("  :bar", "bar")]
    #[case::namespaced(":ns/name", "ns/name")]
    fn test_edn_stream_deserializer_read_keyword_name(#[case] input: &str, #[case] expected: &str) {
        let mut d = EdnStreamDeserializer::new(input);
        assert_eq!(d.read_keyword_name().unwrap().as_ref(), expected);
    }

    #[rstest]
    #[case::not_keyword("foo")]
    #[case::leading_slash(":/foo")]
    fn test_edn_stream_deserializer_read_keyword_name_error(#[case] input: &str) {
        let mut d = EdnStreamDeserializer::new(input);
        assert!(d.read_keyword_name().is_err());
    }

    #[rstest]
    #[case::simple(r#""hello""#, "hello")]
    #[case::with_ws(r#"  "world""#, "world")]
    #[case::empty(r#""""#, "")]
    fn test_edn_stream_deserializer_read_string(#[case] input: &str, #[case] expected: &str) {
        let mut d = EdnStreamDeserializer::new(input);
        assert_eq!(d.read_string().unwrap().as_ref(), expected);
    }

    #[rstest]
    #[case::tab(r#""a\tb""#, "a\tb")]
    #[case::newline(r#""a\nb""#, "a\nb")]
    #[case::carriage_return(r#""a\rb""#, "a\rb")]
    #[case::backslash(r#""a\\b""#, "a\\b")]
    #[case::quote(r#""a\"b""#, "a\"b")]
    #[case::unicode(r#""a\u0041b""#, "aAb")]
    fn test_edn_stream_deserializer_read_string_escapes(
        #[case] input: &str,
        #[case] expected: &str,
    ) {
        let mut d = EdnStreamDeserializer::new(input);
        assert_eq!(d.read_string().unwrap().as_ref(), expected);
    }

    #[rstest]
    #[case::not_string("foo")]
    #[case::unterminated(r#""abc"#)]
    #[case::bad_escape(r#""a\qb""#)]
    #[case::truncated_unicode(r#""a\u00""#)]
    #[case::escape_at_eof(r#""abc\"#)]
    fn test_edn_stream_deserializer_read_string_error(#[case] input: &str) {
        let mut d = EdnStreamDeserializer::new(input);
        assert!(d.read_string().is_err());
    }

    #[rstest]
    #[case::string(r#""hello""#, "hello")]
    #[case::keyword(":foo", "foo")]
    fn test_edn_stream_deserializer_read_string_or_keyword(
        #[case] input: &str,
        #[case] expected: &str,
    ) {
        let mut d = EdnStreamDeserializer::new(input);
        assert_eq!(d.read_string_or_keyword().unwrap().as_ref(), expected);
    }

    #[rstest]
    #[case::zero("0", 0)]
    #[case::positive("123", 123)]
    #[case::negative("-7", -7)]
    #[case::plus_sign("+5", 5)]
    #[case::with_ws("  99", 99)]
    fn test_edn_stream_deserializer_read_i64(#[case] input: &str, #[case] expected: i64) {
        let mut d = EdnStreamDeserializer::new(input);
        assert_eq!(d.read_i64().unwrap(), expected);
    }

    #[rstest]
    #[case::number("123")]
    #[case::boolean("true")]
    #[case::string(r#""hello""#)]
    #[case::vector("[1 2 3]")]
    #[case::list("(1 2)")]
    #[case::map("{:a 1}")]
    #[case::set("#{1 2}")]
    #[case::tagged("#tag value")]
    #[case::discard("#_ foo bar")]
    #[case::char(r#"\c"#)]
    #[case::keyword(":keyword")]
    fn test_edn_stream_deserializer_read_skip_value(#[case] input: &str) {
        let mut d = EdnStreamDeserializer::new(input);
        d.read_skip_value().unwrap();
    }

    #[test]
    fn test_edn_stream_deserializer_read_skip_value_nested() {
        let mut d = EdnStreamDeserializer::new("[{:a [1 2]} #{3}]");
        d.read_skip_value().unwrap();
        assert!(d.expect_eof().is_ok());
    }

    #[rstest]
    #[case::eof("")]
    #[case::unterminated_vector("[1 2")]
    #[case::unterminated_string(r#""abc"#)]
    #[case::unexpected_close(")")]
    fn test_edn_stream_deserializer_read_skip_value_error(#[case] input: &str) {
        let mut d = EdnStreamDeserializer::new(input);
        assert!(d.read_skip_value().is_err());
    }

    #[rstest]
    #[case::number("123 rest", "123")]
    #[case::vector("[1 2] rest", "[1 2]")]
    #[case::boolean("  true rest", "true")]
    #[case::string(r#""hi" rest"#, r#""hi""#)]
    fn test_edn_stream_deserializer_read_value_slice(#[case] input: &str, #[case] expected: &str) {
        let mut d = EdnStreamDeserializer::new(input);
        assert_eq!(d.read_value_slice().unwrap(), expected);
    }

    #[test]
    fn test_edn_stream_deserializer_read_subgrammar() {
        let mut d = EdnStreamDeserializer::new(r#""123""#);
        let v: i64 = d.read_subgrammar("test").unwrap();
        assert_eq!(v, 123);
    }

    #[test]
    fn test_edn_stream_deserializer_read_subgrammar_error() {
        let mut d = EdnStreamDeserializer::new(r#""not_a_number""#);
        let err = d.read_subgrammar::<i64>("test").unwrap_err();
        assert!(matches!(
            err,
            EdnError::De(DeError::Subgrammar {
                grammar: "test",
                ..
            })
        ));
    }

    #[test]
    fn test_edn_stream_deserializer_read_subgrammar_all() {
        let mut d = EdnStreamDeserializer::new(r#""3""#);
        let v: i64 = d.read_subgrammar_all("test").unwrap();
        assert_eq!(v, 3);
    }

    #[test]
    fn test_edn_stream_deserializer_read_subgrammar_all_trailing_content() {
        // Trailing content after the string literal must be rejected (expect_eof).
        let mut d = EdnStreamDeserializer::new(r#""5" extra"#);
        let err = d.read_subgrammar_all::<i64>("test").unwrap_err();
        assert!(matches!(
            err,
            EdnError::Parse(ParseError::TrailingContent { .. })
        ));
    }

    #[test]
    fn test_edn_stream_deserializer_read_subgrammar_all_subgrammar_error() {
        let mut d = EdnStreamDeserializer::new(r#""nope""#);
        let err = d.read_subgrammar_all::<i64>("test").unwrap_err();
        assert!(matches!(
            err,
            EdnError::De(DeError::Subgrammar {
                grammar: "test",
                ..
            })
        ));
    }

    #[test]
    fn test_edn_stream_deserializer_skip_comment() {
        let mut d = EdnStreamDeserializer::new("; comment\n123");
        assert_eq!(d.read_i64().unwrap(), 123);
    }

    #[test]
    fn test_edn_stream_deserializer_skip_discard() {
        let mut d = EdnStreamDeserializer::new("#_ [ignored] 7");
        assert_eq!(d.read_i64().unwrap(), 7);
    }

    #[test]
    fn test_edn_stream_deserializer_skip_nested_discard() {
        let mut d = EdnStreamDeserializer::new("#_ #_ a b 9");
        assert_eq!(d.read_i64().unwrap(), 9);
    }

    #[test]
    fn test_edn_stream_deserializer_skip_string_with_escapes() {
        let mut d = EdnStreamDeserializer::new(r#""a\"b\u0041c" 5"#);
        d.read_skip_value().unwrap();
        assert_eq!(d.read_i64().unwrap(), 5);
    }

    #[rstest]
    #[case::decimal("3.14")]
    #[case::exponent("1e2")]
    #[case::letters("abc")]
    #[case::bare_sign("+")]
    #[case::leading_zeros("00123")]
    fn test_edn_stream_deserializer_read_i64_not_i64(#[case] input: &str) {
        let mut d = EdnStreamDeserializer::new(input);
        assert!(d.read_i64().is_err());
    }

    #[test]
    fn test_edn_stream_deserializer_sequential_reads() {
        let mut d = EdnStreamDeserializer::new("{:name \"water\" :charge -1}");
        d.consume_byte(b'{').unwrap();
        assert_eq!(d.read_keyword_name().unwrap().as_ref(), "name");
        assert_eq!(d.read_string().unwrap().as_ref(), "water");
        assert_eq!(d.read_keyword_name().unwrap().as_ref(), "charge");
        assert_eq!(d.read_i64().unwrap(), -1);
        d.consume_byte(b'}').unwrap();
        d.expect_eof().unwrap();
    }

    #[rstest]
    #[case::eof("")]
    #[case::number("123")]
    fn test_edn_stream_deserializer_read_string_or_keyword_error(#[case] input: &str) {
        let mut d = EdnStreamDeserializer::new(input);
        assert!(d.read_string_or_keyword().is_err());
    }
}
