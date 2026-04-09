//! Streaming EDN reader.
//!
//! `EdnStreamDeserializer` exposes byte-level primitives for fused
//! `FromEdn::from_edn_str` overrides. It is a low-level primitive used by
//! the `umol-edn-macros` derive output and by hand-written hot-path
//! parsers; it is not a serde Deserializer.

use std::borrow::Cow;
use std::str::from_utf8;

use crate::error::ParseError;
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
