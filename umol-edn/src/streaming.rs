//! Streaming serde Deserializer that parses EDN directly into Rust types
//! without building an intermediate `Edn` value tree.

use std::borrow::Cow;
use std::str::from_utf8;

use serde::de::{
    self, DeserializeSeed, EnumAccess, MapAccess, SeqAccess, VariantAccess, Visitor,
};

use crate::config::{ParseConfig, TagReaders};
use crate::error::EdnError;
use crate::parser::{is_symbol_char, is_symbol_start};

const MAX_DEPTH: u16 = 128;

pub struct EdnStreamDeserializer<'de> {
    input: &'de str,
    pos: usize,
    scratch: Vec<u8>,
    depth: u16,
    tag_readers: TagReaders,
}

impl<'de> EdnStreamDeserializer<'de> {
    pub fn new(input: &'de str) -> Self {
        Self::with_config(input, &ParseConfig::default())
    }

    pub fn with_config(input: &'de str, config: &ParseConfig) -> Self {
        Self {
            input,
            pos: 0,
            scratch: Vec::new(),
            depth: MAX_DEPTH,
            tag_readers: config.tag_readers.clone(),
        }
    }

    pub fn expect_eof(&mut self) -> Result<(), EdnError> {
        self.skip_ws()?;
        if self.pos < self.input.len() {
            Err(EdnError::TrailingContent { offset: self.pos })
        } else {
            Ok(())
        }
    }

    // --- Low-level ---

    #[inline]
    fn peek(&self) -> Option<u8> {
        self.input.as_bytes().get(self.pos).copied()
    }

    fn next_byte(&mut self) -> Result<u8, EdnError> {
        match self.input.as_bytes().get(self.pos) {
            Some(&b) => {
                self.pos += 1;
                Ok(b)
            }
            None => Err(EdnError::UnexpectedEof { offset: self.pos }),
        }
    }

    #[inline]
    fn skip_ws(&mut self) -> Result<(), EdnError> {
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
    fn skip_ws_discard(&mut self) -> Result<(), EdnError> {
        self.pos += 2;
        self.skip_ws()?;
        self.skip_value()?;
        self.skip_ws()
    }

    fn enter_scope(&mut self) -> Result<(), EdnError> {
        self.depth = self.depth.checked_sub(1).ok_or_else(|| {
            EdnError::Custom("recursion limit exceeded".to_string())
        })?;
        Ok(())
    }

    fn leave_scope(&mut self) {
        self.depth += 1;
    }

    // --- Token parsing ---

    #[inline]
    fn parse_keyword_name(&mut self) -> Result<Cow<'de, str>, EdnError> {
        debug_assert_eq!(self.input.as_bytes()[self.pos], b':');
        self.pos += 1;

        // :/ and :/foo are not legal keywords per the EDN spec.
        if self.input.as_bytes().get(self.pos) == Some(&b'/') {
            return Err(EdnError::UnexpectedToken { offset: self.pos, found: '/' });
        }
        self.parse_symbol_str().map(Cow::Borrowed)
    }

    #[inline]
    fn parse_symbol_str(&mut self) -> Result<&'de str, EdnError> {
        let start = self.pos;
        let bytes = self.input.as_bytes();
        if self.pos >= bytes.len() {
            return Err(EdnError::UnexpectedEof { offset: start });
        }
        let first = bytes[self.pos] as char;
        if !is_symbol_start(first) && !first.is_ascii_digit() {
            return Err(EdnError::UnexpectedToken {
                offset: start,
                found: first,
            });
        }
        self.pos += 1;
        while self.pos < bytes.len() && is_symbol_char(bytes[self.pos] as char) {
            self.pos += 1;
        }
        let s = &self.input[start..self.pos];
        self.validate_symbol_str(s, start)?;
        Ok(s)
    }

    /// Validate symbol slash rules (mirrors parser::validate_symbol).
    #[inline]
    fn validate_symbol_str(&self, s: &str, offset: usize) -> Result<(), EdnError> {
        if s == "/" {
            return Ok(());
        }
        if let Some(slash_pos) = s.find('/') {
            let prefix = &s[..slash_pos];
            let name = &s[slash_pos + 1..];
            if prefix.is_empty() || name.is_empty() {
                return Err(EdnError::InvalidSymbol { offset });
            }
            let first_name_char = name.chars().next().unwrap();
            if first_name_char.is_ascii_digit() {
                return Err(EdnError::InvalidSymbol { offset });
            }
            if name.contains('/') {
                return Err(EdnError::InvalidSymbol { offset });
            }
            if !is_symbol_start(first_name_char) {
                return Err(EdnError::InvalidSymbol { offset });
            }
        }
        Ok(())
    }

    /// Skip a `#tag` prefix if present. Rejects bare (unqualified) tags unless registered.
    /// Recurses for nested tags like `#a #b value`.
    #[inline]
    fn skip_tag_if_present(&mut self) -> Result<(), EdnError> {
        let bytes = self.input.as_bytes();
        if self.pos < bytes.len() && bytes[self.pos] == b'#' {
            match bytes.get(self.pos + 1) {
                Some(b'{') | Some(b'_') => return Ok(()),
                _ => {}
            }
            self.pos += 1;
            let offset = self.pos;
            let tag = self.parse_symbol_str()?;
            if !tag.contains('/') && self.tag_readers.get(tag).is_none() {
                return Err(EdnError::InvalidTag {
                    offset,
                    tag: tag.to_string(),
                });
            }
            self.skip_ws()?;
            self.skip_tag_if_present()?;
        }
        Ok(())
    }

    fn parse_number_i64(&mut self) -> Result<i64, EdnError> {
        let s = self.scan_number_str()?;
        s.parse::<i64>()
            .map_err(|_| EdnError::InvalidNumber { offset: self.pos - s.len() })
    }

    fn parse_number_f64(&mut self) -> Result<f64, EdnError> {
        let s = self.scan_number_str()?;
        s.parse::<f64>()
            .map_err(|_| EdnError::InvalidNumber { offset: self.pos - s.len() })
    }

    fn scan_number_str(&mut self) -> Result<&'de str, EdnError> {
        let start = self.pos;
        let bytes = self.input.as_bytes();
        // optional sign
        if self.pos < bytes.len() && matches!(bytes[self.pos], b'+' | b'-') {
            self.pos += 1;
        }
        // digits
        if self.pos >= bytes.len() {
            return Err(EdnError::UnexpectedEof { offset: start });
        }
        if !bytes[self.pos].is_ascii_digit() {
            return Err(EdnError::InvalidNumber { offset: start });
        }
        while self.pos < bytes.len() && bytes[self.pos].is_ascii_digit() {
            self.pos += 1;
        }
        // optional .digits
        if self.pos < bytes.len() && bytes[self.pos] == b'.' {
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
        // reject bignum suffix
        if self.pos < bytes.len() && matches!(bytes[self.pos], b'N' | b'M') {
            return Err(EdnError::UnsupportedFeature {
                offset: start,
                feature: "bignum",
            });
        }
        Ok(&self.input[start..self.pos])
    }

    fn parse_number_any<V: Visitor<'de>>(
        &mut self,
        visitor: V,
    ) -> Result<V::Value, EdnError> {
        let s = self.scan_number_str()?;
        if memchr::memchr3(b'.', b'e', b'E', s.as_bytes()).is_some() {
            let f: f64 = s
                .parse()
                .map_err(|_| EdnError::InvalidNumber { offset: self.pos - s.len() })?;
            visitor.visit_f64(f)
        } else {
            let n: i64 = s
                .parse()
                .map_err(|_| EdnError::InvalidNumber { offset: self.pos - s.len() })?;
            visitor.visit_i64(n)
        }
    }

    fn parse_string(&mut self) -> Result<Cow<'de, str>, EdnError> {
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
            self.scratch.extend_from_slice(&bytes[self.pos..self.pos + bs]);
            self.pos += bs;
        }
        loop {
            if self.pos >= bytes.len() {
                return Err(EdnError::UnexpectedEof { offset: start });
            }
            match bytes[self.pos] {
                b'"' => {
                    self.pos += 1;
                    let s = from_utf8(&self.scratch)
                        .map_err(|_| EdnError::Custom("invalid UTF-8".to_string()))?;
                    return Ok(Cow::Owned(s.to_string()));
                }
                b'\\' => {
                    let esc_offset = self.pos;
                    self.pos += 1;
                    let esc = self.next_byte().map_err(|_| EdnError::InvalidEscape { offset: esc_offset })?;
                    match esc {
                        b't' => self.scratch.push(b'\t'),
                        b'r' => self.scratch.push(b'\r'),
                        b'n' => self.scratch.push(b'\n'),
                        b'\\' => self.scratch.push(b'\\'),
                        b'"' => self.scratch.push(b'"'),
                        b'u' => {
                            if self.pos + 4 > bytes.len() {
                                return Err(EdnError::InvalidEscape { offset: esc_offset });
                            }
                            let hex = &self.input[self.pos..self.pos + 4];
                            let cp = u32::from_str_radix(hex, 16)
                                .map_err(|_| EdnError::InvalidEscape { offset: esc_offset })?;
                            let ch = char::from_u32(cp)
                                .ok_or(EdnError::InvalidEscape { offset: esc_offset })?;
                            self.pos += 4;
                            let mut buf = [0u8; 4];
                            let encoded = ch.encode_utf8(&mut buf);
                            self.scratch.extend_from_slice(encoded.as_bytes());
                        }
                        _ => return Err(EdnError::InvalidEscape { offset: esc_offset }),
                    }
                    // Batch copy until next " or \.
                    if let Some(span) = memchr::memchr2(b'"', b'\\', &bytes[self.pos..]) {
                        self.scratch.extend_from_slice(&bytes[self.pos..self.pos + span]);
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

    /// Visit a string value, using borrowed str for no-escape case.
    fn parse_string_visitor<V: Visitor<'de>>(
        &mut self,
        visitor: V,
    ) -> Result<V::Value, EdnError> {
        let s = self.parse_string()?;
        match s {
            Cow::Borrowed(b) => visitor.visit_borrowed_str(b),
            Cow::Owned(o) => visitor.visit_string(o),
        }
    }

    // --- Value skipping (for ignored_any / #_ discard) ---

    fn skip_value(&mut self) -> Result<(), EdnError> {
        self.skip_ws()?;
        let b = self.peek().ok_or(EdnError::UnexpectedEof { offset: self.pos })?;
        match b {
            b'(' => { self.pos += 1; self.skip_delimited(b')') }
            b'[' => { self.pos += 1; self.skip_delimited(b']') }
            b'{' => { self.pos += 1; self.skip_delimited(b'}') }
            b'#' => {
                self.pos += 1;
                match self.peek() {
                    Some(b'{') => { self.pos += 1; self.skip_delimited(b'}') }
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
            b'\\' => { self.pos += 1; self.skip_atom() }
            b':' => { self.pos += 1; self.skip_atom() }
            _ => self.skip_atom(),
        }
    }

    fn skip_delimited(&mut self, close: u8) -> Result<(), EdnError> {
        loop {
            self.skip_ws()?;
            match self.peek() {
                None => return Err(EdnError::UnexpectedEof { offset: self.pos }),
                Some(b) if b == close => {
                    self.pos += 1;
                    return Ok(());
                }
                _ => self.skip_value()?,
            }
        }
    }

    fn skip_string(&mut self) -> Result<(), EdnError> {
        debug_assert_eq!(self.input.as_bytes()[self.pos], b'"');
        self.pos += 1;
        let bytes = self.input.as_bytes();
        loop {
            match memchr::memchr2(b'"', b'\\', &bytes[self.pos..]) {
                None => return Err(EdnError::UnexpectedEof { offset: self.pos }),
                Some(i) => {
                    self.pos += i;
                    if bytes[self.pos] == b'"' {
                        self.pos += 1;
                        return Ok(());
                    }
                    // Skip escape: backslash + one char
                    self.pos += 2;
                }
            }
        }
    }

    fn skip_atom(&mut self) -> Result<(), EdnError> {
        let bytes = self.input.as_bytes();
        if self.pos >= bytes.len() {
            return Err(EdnError::UnexpectedEof { offset: self.pos });
        }
        while self.pos < bytes.len() && is_symbol_char(bytes[self.pos] as char) {
            self.pos += 1;
        }
        Ok(())
    }
}

fn visit_cow_str<'de, V: Visitor<'de>>(visitor: V, s: Cow<'de, str>) -> Result<V::Value, EdnError> {
    match s {
        Cow::Borrowed(b) => visitor.visit_borrowed_str(b),
        Cow::Owned(o) => visitor.visit_string(o),
    }
}

// ---------------------------------------------------------------------------
// Deserializer trait
// ---------------------------------------------------------------------------

impl<'de, 'a> de::Deserializer<'de> for &'a mut EdnStreamDeserializer<'de> {
    type Error = EdnError;

    fn deserialize_any<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        self.skip_ws()?;
        let b = self.peek().ok_or(EdnError::UnexpectedEof { offset: self.pos })?;
        match b {
            b'"' => self.parse_string_visitor(visitor),
            b':' => {
                let name = self.parse_keyword_name()?;
                visit_cow_str(visitor, name)
            }
            b'(' | b'[' => self.deserialize_seq(visitor),
            b'{' => self.deserialize_map(visitor),
            b'#' => {
                let bytes = self.input.as_bytes();
                match bytes.get(self.pos + 1) {
                    Some(b'{') => {
                        // Set #{...} — deserialize as seq
                        self.enter_scope()?;
                        self.pos += 2; // skip #{
                        let result = visitor.visit_seq(SeqAccessor {
                            de: self,
                            closing: b'}',
                            first: true,
                        });
                        self.leave_scope();
                        result
                    }
                    Some(b'_') => {
                        self.pos += 2; // skip #_
                        self.skip_ws()?;
                        self.skip_value()?;
                        return self.deserialize_any(visitor);
                    }
                    _ => {
                        // Tagged literal: #tag value — unwrap tag, deserialize value.
                        self.pos += 1; // skip #
                        let offset = self.pos;
                        let tag = self.parse_symbol_str()?;
                        if !tag.contains('/') && self.tag_readers.get(tag).is_none() {
                            return Err(EdnError::InvalidTag {
                                offset,
                                tag: tag.to_string(),
                            });
                        }
                        self.skip_ws()?;
                        self.deserialize_any(visitor)
                    }
                }
            }
            b'\\' => {
                // Character literal
                self.pos += 1;
                let ch = self.parse_char_literal()?;
                visitor.visit_char(ch)
            }
            b'+' | b'-' => {
                // Could be number or symbol.
                let second = self.input.as_bytes().get(self.pos + 1).copied();
                if matches!(second, Some(b'0'..=b'9')) {
                    self.parse_number_any(visitor)
                } else {
                    let s = self.parse_symbol_str()?;
                    visitor.visit_borrowed_str(s)
                }
            }
            b'0'..=b'9' => self.parse_number_any(visitor),
            _ => {
                let s = self.parse_symbol_str()?;
                match s {
                    "nil" => visitor.visit_unit(),
                    "true" => visitor.visit_bool(true),
                    "false" => visitor.visit_bool(false),
                    _ => visitor.visit_borrowed_str(s),
                }
            }
        }
    }

    fn deserialize_bool<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        self.skip_ws()?;
        self.skip_tag_if_present()?;
        let s = self.parse_symbol_str()?;
        match s {
            "true" => visitor.visit_bool(true),
            "false" => visitor.visit_bool(false),
            _ => Err(EdnError::Custom(format!("expected bool, got {s}"))),
        }
    }

    fn deserialize_i8<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        self.skip_ws()?;
        self.skip_tag_if_present()?;
        let n = self.parse_number_i64()?;
        let v = i8::try_from(n).map_err(|_| EdnError::Custom(format!("{n} out of range for i8")))?;
        visitor.visit_i8(v)
    }

    fn deserialize_i16<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        self.skip_ws()?;
        self.skip_tag_if_present()?;
        let n = self.parse_number_i64()?;
        let v = i16::try_from(n).map_err(|_| EdnError::Custom(format!("{n} out of range for i16")))?;
        visitor.visit_i16(v)
    }

    fn deserialize_i32<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        self.skip_ws()?;
        self.skip_tag_if_present()?;
        let n = self.parse_number_i64()?;
        let v = i32::try_from(n).map_err(|_| EdnError::Custom(format!("{n} out of range for i32")))?;
        visitor.visit_i32(v)
    }

    fn deserialize_i64<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        self.skip_ws()?;
        self.skip_tag_if_present()?;
        visitor.visit_i64(self.parse_number_i64()?)
    }

    fn deserialize_u8<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        self.skip_ws()?;
        self.skip_tag_if_present()?;
        let n = self.parse_number_i64()?;
        let v = u8::try_from(n).map_err(|_| EdnError::Custom(format!("{n} out of range for u8")))?;
        visitor.visit_u8(v)
    }

    fn deserialize_u16<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        self.skip_ws()?;
        self.skip_tag_if_present()?;
        let n = self.parse_number_i64()?;
        let v = u16::try_from(n).map_err(|_| EdnError::Custom(format!("{n} out of range for u16")))?;
        visitor.visit_u16(v)
    }

    fn deserialize_u32<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        self.skip_ws()?;
        self.skip_tag_if_present()?;
        let n = self.parse_number_i64()?;
        let v = u32::try_from(n).map_err(|_| EdnError::Custom(format!("{n} out of range for u32")))?;
        visitor.visit_u32(v)
    }

    fn deserialize_u64<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        self.skip_ws()?;
        self.skip_tag_if_present()?;
        let n = self.parse_number_i64()?;
        let v = u64::try_from(n).map_err(|_| EdnError::Custom(format!("{n} out of range for u64")))?;
        visitor.visit_u64(v)
    }

    fn deserialize_f32<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        self.skip_ws()?;
        self.skip_tag_if_present()?;
        visitor.visit_f32(self.parse_number_f64()? as f32)
    }

    fn deserialize_f64<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        self.skip_ws()?;
        self.skip_tag_if_present()?;
        visitor.visit_f64(self.parse_number_f64()?)
    }

    fn deserialize_char<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        self.skip_ws()?;
        self.skip_tag_if_present()?;
        if self.peek() == Some(b'\\') {
            self.pos += 1;
            visitor.visit_char(self.parse_char_literal()?)
        } else {
            // Accept a single-char string as char.
            self.deserialize_any(visitor)
        }
    }

    fn deserialize_str<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        self.skip_ws()?;
        self.skip_tag_if_present()?;
        match self.peek() {
            Some(b'"') => self.parse_string_visitor(visitor),
            Some(b':') => {
                let name = self.parse_keyword_name()?;
                visit_cow_str(visitor, name)
            }
            _ => {
                let s = self.parse_symbol_str()?;
                visitor.visit_borrowed_str(s)
            }
        }
    }

    fn deserialize_string<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        self.deserialize_str(visitor)
    }

    fn deserialize_bytes<V: Visitor<'de>>(self, _visitor: V) -> Result<V::Value, Self::Error> {
        Err(EdnError::Custom("bytes not supported".to_string()))
    }

    fn deserialize_byte_buf<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        self.deserialize_bytes(visitor)
    }

    fn deserialize_option<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        self.skip_ws()?;
        self.skip_tag_if_present()?;
        // Check for nil: 'n' followed by 'i', 'l', then non-symbol-char.
        let bytes = self.input.as_bytes();
        if self.pos + 3 <= bytes.len()
            && bytes[self.pos] == b'n'
            && bytes[self.pos + 1] == b'i'
            && bytes[self.pos + 2] == b'l'
            && !bytes.get(self.pos + 3).map_or(false, |&b| is_symbol_char(b as char))
        {
            self.pos += 3;
            visitor.visit_none()
        } else {
            visitor.visit_some(self)
        }
    }

    fn deserialize_unit<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        self.skip_ws()?;
        self.skip_tag_if_present()?;
        let s = self.parse_symbol_str()?;
        if s == "nil" {
            visitor.visit_unit()
        } else {
            Err(EdnError::Custom(format!("expected nil, got {s}")))
        }
    }

    fn deserialize_unit_struct<V: Visitor<'de>>(
        self,
        _name: &'static str,
        visitor: V,
    ) -> Result<V::Value, Self::Error> {
        self.deserialize_unit(visitor)
    }

    fn deserialize_newtype_struct<V: Visitor<'de>>(
        self,
        _name: &'static str,
        visitor: V,
    ) -> Result<V::Value, Self::Error> {
        visitor.visit_newtype_struct(self)
    }

    fn deserialize_seq<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        self.skip_ws()?;
        self.skip_tag_if_present()?;
        let b = self.peek().ok_or(EdnError::UnexpectedEof { offset: self.pos })?;
        let closing = match b {
            b'[' => b']',
            b'(' => b')',
            _ => return Err(EdnError::UnexpectedToken { offset: self.pos, found: b as char }),
        };
        self.enter_scope()?;
        self.pos += 1;
        let result = visitor.visit_seq(SeqAccessor {
            de: self,
            closing,
            first: true,
        })?;
        // Fixed-size tuples don't ask for the terminal None element, so the
        // closing bracket may not have been consumed yet.
        self.skip_ws()?;
        if self.peek() == Some(closing) {
            self.pos += 1;
        }
        self.leave_scope();
        Ok(result)
    }

    fn deserialize_tuple<V: Visitor<'de>>(
        self,
        _len: usize,
        visitor: V,
    ) -> Result<V::Value, Self::Error> {
        self.deserialize_seq(visitor)
    }

    fn deserialize_tuple_struct<V: Visitor<'de>>(
        self,
        _name: &'static str,
        _len: usize,
        visitor: V,
    ) -> Result<V::Value, Self::Error> {
        self.deserialize_seq(visitor)
    }

    fn deserialize_map<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        self.skip_ws()?;
        self.skip_tag_if_present()?;
        match self.peek() {
            Some(b'{') => {
                self.enter_scope()?;
                self.pos += 1;
                let result = visitor.visit_map(MapAccessor { de: self, first: true });
                self.leave_scope();
                result
            }
            _ => {
                // nil as empty map (matching existing behavior)
                let bytes = self.input.as_bytes();
                if self.pos + 3 <= bytes.len()
                    && bytes[self.pos] == b'n'
                    && bytes[self.pos + 1] == b'i'
                    && bytes[self.pos + 2] == b'l'
                    && !bytes.get(self.pos + 3).map_or(false, |&b| is_symbol_char(b as char))
                {
                    self.pos += 3;
                    visitor.visit_map(EmptyMapAccessor)
                } else {
                    Err(EdnError::UnexpectedToken {
                        offset: self.pos,
                        found: self.peek().map(|b| b as char).unwrap_or('\0'),
                    })
                }
            }
        }
    }

    fn deserialize_struct<V: Visitor<'de>>(
        self,
        _name: &'static str,
        _fields: &'static [&'static str],
        visitor: V,
    ) -> Result<V::Value, Self::Error> {
        self.deserialize_map(visitor)
    }

    fn deserialize_enum<V: Visitor<'de>>(
        self,
        _name: &'static str,
        _variants: &'static [&'static str],
        visitor: V,
    ) -> Result<V::Value, Self::Error> {
        self.skip_ws()?;
        match self.peek() {
            Some(b':') => {
                let name = self.parse_keyword_name()?;
                match name {
                    Cow::Borrowed(b) => visitor.visit_enum(de::value::BorrowedStrDeserializer::new(b)),
                    Cow::Owned(o) => visitor.visit_enum(de::value::StrDeserializer::<EdnError>::new(&o)),
                }
            }
            Some(b'"') => {
                let s = self.parse_string()?;
                match s {
                    Cow::Borrowed(b) => {
                        visitor.visit_enum(de::value::BorrowedStrDeserializer::new(b))
                    }
                    Cow::Owned(o) => {
                        visitor.visit_enum(de::value::StrDeserializer::<EdnError>::new(&o))
                    }
                }
            }
            Some(b'#') => {
                let bytes = self.input.as_bytes();
                // Don't match #{, #_ — only tagged literals.
                if matches!(bytes.get(self.pos + 1), Some(b'{') | Some(b'_')) {
                    return self.deserialize_any(visitor);
                }
                self.pos += 1; // skip #
                let tag = self.parse_symbol_str()?;
                self.skip_ws()?;
                visitor.visit_enum(StreamingTaggedEnumAccess { de: self, tag })
            }
            _ => self.deserialize_any(visitor),
        }
    }

    fn deserialize_identifier<V: Visitor<'de>>(
        self,
        visitor: V,
    ) -> Result<V::Value, Self::Error> {
        self.skip_ws()?;
        match self.peek() {
            Some(b':') => {
                let name = self.parse_keyword_name()?;
                visit_cow_str(visitor, name)
            }
            Some(b'"') => self.parse_string_visitor(visitor),
            _ => {
                let s = self.parse_symbol_str()?;
                visitor.visit_borrowed_str(s)
            }
        }
    }

    fn deserialize_ignored_any<V: Visitor<'de>>(
        self,
        visitor: V,
    ) -> Result<V::Value, Self::Error> {
        self.skip_value()?;
        visitor.visit_unit()
    }
}

// ---------------------------------------------------------------------------
// Character literal parsing
// ---------------------------------------------------------------------------

impl<'de> EdnStreamDeserializer<'de> {
    fn parse_char_literal(&mut self) -> Result<char, EdnError> {
        let offset = self.pos;
        let bytes = self.input.as_bytes();
        if self.pos >= bytes.len() {
            return Err(EdnError::UnexpectedEof { offset });
        }

        // Fast path: single character followed by non-symbol-char or EOF.
        let first = bytes[self.pos];
        let second = bytes.get(self.pos + 1).copied();
        let is_single = second.is_none()
            || !is_symbol_char(second.unwrap() as char)
            || matches!(second, Some(b' ' | b'\t' | b'\n' | b'\r' | b','));
        if is_single && first != b'u' {
            let c = self.input[self.pos..].chars().next().unwrap();
            if c.is_whitespace() {
                return Err(EdnError::InvalidCharLiteral { offset });
            }
            self.pos += c.len_utf8();
            return Ok(c);
        }

        // Unicode escape \uNNNN
        if first == b'u' && self.pos + 5 <= bytes.len() {
            let hex = &self.input[self.pos + 1..self.pos + 5];
            if hex.bytes().all(|b| b.is_ascii_hexdigit()) {
                let cp = u32::from_str_radix(hex, 16)
                    .map_err(|_| EdnError::InvalidCharLiteral { offset })?;
                let ch = char::from_u32(cp)
                    .ok_or(EdnError::InvalidCharLiteral { offset })?;
                self.pos += 5;
                return Ok(ch);
            }
        }

        // Named characters
        let named: &[(&str, char)] = &[
            ("newline", '\n'), ("return", '\r'), ("space", ' '), ("tab", '\t'),
        ];
        let remaining = &self.input[self.pos..];
        for &(name, ch) in named {
            if remaining.starts_with(name) {
                let after = bytes.get(self.pos + name.len()).copied();
                let terminates = after.is_none()
                    || !is_symbol_char(after.unwrap() as char)
                    || matches!(after, Some(b' ' | b'\t' | b'\n' | b'\r' | b','));
                if terminates {
                    self.pos += name.len();
                    return Ok(ch);
                }
            }
        }

        Err(EdnError::InvalidCharLiteral { offset })
    }
}

// ---------------------------------------------------------------------------
// SeqAccessor
// ---------------------------------------------------------------------------

struct SeqAccessor<'a, 'de> {
    de: &'a mut EdnStreamDeserializer<'de>,
    closing: u8,
    first: bool,
}

impl<'a, 'de> SeqAccess<'de> for SeqAccessor<'a, 'de> {
    type Error = EdnError;

    fn next_element_seed<T: DeserializeSeed<'de>>(
        &mut self,
        seed: T,
    ) -> Result<Option<T::Value>, Self::Error> {
        self.de.skip_ws()?;
        match self.de.peek() {
            None => Err(EdnError::UnexpectedEof { offset: self.de.pos }),
            Some(b) if b == self.closing => {
                self.de.pos += 1;
                Ok(None)
            }
            _ => {
                self.first = false;
                seed.deserialize(&mut *self.de).map(Some)
            }
        }
    }
}

// ---------------------------------------------------------------------------
// MapAccessor
// ---------------------------------------------------------------------------

struct MapAccessor<'a, 'de> {
    de: &'a mut EdnStreamDeserializer<'de>,
    first: bool,
}

impl<'a, 'de> MapAccess<'de> for MapAccessor<'a, 'de> {
    type Error = EdnError;

    fn next_key_seed<K: DeserializeSeed<'de>>(
        &mut self,
        seed: K,
    ) -> Result<Option<K::Value>, Self::Error> {
        self.de.skip_ws()?;
        match self.de.peek() {
            None => Err(EdnError::UnexpectedEof { offset: self.de.pos }),
            Some(b'}') => {
                self.de.pos += 1;
                Ok(None)
            }
            _ => {
                self.first = false;
                seed.deserialize(&mut *self.de).map(Some)
            }
        }
    }

    fn next_value_seed<V: DeserializeSeed<'de>>(
        &mut self,
        seed: V,
    ) -> Result<V::Value, Self::Error> {
        seed.deserialize(&mut *self.de)
    }
}

// ---------------------------------------------------------------------------
// EmptyMapAccessor (for nil → empty map)
// ---------------------------------------------------------------------------

struct EmptyMapAccessor;

impl<'de> MapAccess<'de> for EmptyMapAccessor {
    type Error = EdnError;

    fn next_key_seed<K: DeserializeSeed<'de>>(
        &mut self,
        _seed: K,
    ) -> Result<Option<K::Value>, Self::Error> {
        Ok(None)
    }

    fn next_value_seed<V: DeserializeSeed<'de>>(
        &mut self,
        _seed: V,
    ) -> Result<V::Value, Self::Error> {
        Err(EdnError::Custom("empty map has no values".to_string()))
    }
}

// ---------------------------------------------------------------------------
// Tagged enum access (for #Variant value round-tripping)
// ---------------------------------------------------------------------------

struct StreamingTaggedEnumAccess<'a, 'de> {
    de: &'a mut EdnStreamDeserializer<'de>,
    tag: &'de str,
}

impl<'a, 'de> EnumAccess<'de> for StreamingTaggedEnumAccess<'a, 'de> {
    type Error = EdnError;
    type Variant = &'a mut EdnStreamDeserializer<'de>;

    fn variant_seed<V: DeserializeSeed<'de>>(
        self,
        seed: V,
    ) -> Result<(V::Value, Self::Variant), Self::Error> {
        let variant =
            seed.deserialize(de::value::BorrowedStrDeserializer::new(self.tag))?;
        Ok((variant, self.de))
    }
}

impl<'a, 'de> VariantAccess<'de> for &'a mut EdnStreamDeserializer<'de> {
    type Error = EdnError;

    fn unit_variant(self) -> Result<(), Self::Error> {
        Ok(())
    }

    fn newtype_variant_seed<T: DeserializeSeed<'de>>(
        self,
        seed: T,
    ) -> Result<T::Value, Self::Error> {
        seed.deserialize(self)
    }

    fn tuple_variant<V: Visitor<'de>>(
        self,
        _len: usize,
        visitor: V,
    ) -> Result<V::Value, Self::Error> {
        serde::Deserializer::deserialize_seq(self, visitor)
    }

    fn struct_variant<V: Visitor<'de>>(
        self,
        _fields: &'static [&'static str],
        visitor: V,
    ) -> Result<V::Value, Self::Error> {
        serde::Deserializer::deserialize_map(self, visitor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use serde::Deserialize;

    fn streaming_from_str<'a, T: Deserialize<'a>>(s: &'a str) -> Result<T, EdnError> {
        let mut de = EdnStreamDeserializer::new(s);
        let val = T::deserialize(&mut de)?;
        de.expect_eof()?;
        Ok(val)
    }

    // --- Primitives ---

    #[rstest]
    #[case("true", true)]
    #[case("false", false)]
    fn test_streaming_bool(#[case] input: &str, #[case] expected: bool) {
        assert_eq!(streaming_from_str::<bool>(input).unwrap(), expected);
    }

    #[rstest]
    #[case("0", 0i64)]
    #[case("12", 12)]
    #[case("-1", -1)]
    #[case("+5", 5)]
    fn test_streaming_i64(#[case] input: &str, #[case] expected: i64) {
        assert_eq!(streaming_from_str::<i64>(input).unwrap(), expected);
    }

    #[rstest]
    #[case("1.0", 1.0f64)]
    #[case("-3.14", -3.14)]
    #[case("1e10", 1e10)]
    fn test_streaming_f64(#[case] input: &str, #[case] expected: f64) {
        let val: f64 = streaming_from_str(input).unwrap();
        assert!((val - expected).abs() < 1e-10);
    }

    #[rstest]
    #[case(r#""hello""#, "hello")]
    #[case(r#""line\nbreak""#, "line\nbreak")]
    #[case(r#""tab\there""#, "tab\there")]
    fn test_streaming_string(#[case] input: &str, #[case] expected: &str) {
        assert_eq!(streaming_from_str::<String>(input).unwrap(), expected);
    }

    #[test]
    fn test_streaming_vec() {
        let v: Vec<i64> = streaming_from_str("[1 2 3]").unwrap();
        assert_eq!(v, vec![1, 2, 3]);
    }

    #[test]
    fn test_streaming_tuple() {
        let t: (String, String, String) = streaming_from_str("[:a :b :single]").unwrap();
        assert_eq!(t, ("a".to_string(), "b".to_string(), "single".to_string()));
    }

    #[test]
    fn test_streaming_option() {
        assert_eq!(streaming_from_str::<Option<i64>>("nil").unwrap(), None);
        assert_eq!(streaming_from_str::<Option<i64>>("5").unwrap(), Some(5));
    }

    // --- Struct ---

    #[derive(Deserialize, Debug, PartialEq)]
    struct MoleculeProxy {
        atoms: Vec<String>,
        bonds: Vec<(String, String, String)>,
    }

    #[test]
    fn test_streaming_struct() {
        let m: MoleculeProxy =
            streaming_from_str(r#"{:atoms [C O] :bonds [["0" "1" :single]]}"#).unwrap();
        assert_eq!(
            m,
            MoleculeProxy {
                atoms: vec!["C".to_string(), "O".to_string()],
                bonds: vec![("0".to_string(), "1".to_string(), "single".to_string())],
            }
        );
    }

    #[test]
    fn test_streaming_nested_struct() {
        #[derive(Deserialize, Debug, PartialEq)]
        struct Inner {
            x: i64,
        }
        #[derive(Deserialize, Debug, PartialEq)]
        struct Outer {
            inner: Inner,
        }
        let o: Outer = streaming_from_str("{:inner {:x 5}}").unwrap();
        assert_eq!(o, Outer { inner: Inner { x: 5 } });
    }

    // --- Enum ---

    #[derive(Deserialize, Debug, PartialEq)]
    enum Color {
        #[serde(rename = "red")]
        Red,
        #[serde(rename = "blue")]
        Blue,
    }

    #[rstest]
    #[case(":red", Color::Red)]
    #[case(":blue", Color::Blue)]
    fn test_streaming_enum(#[case] input: &str, #[case] expected: Color) {
        assert_eq!(streaming_from_str::<Color>(input).unwrap(), expected);
    }

    // --- HashMap ---

    #[test]
    fn test_streaming_hashmap() {
        use std::collections::HashMap;
        let m: HashMap<String, i64> = streaming_from_str("{:a 1 :b 2}").unwrap();
        assert_eq!(m.len(), 2);
        assert_eq!(m["a"], 1);
        assert_eq!(m["b"], 2);
    }

    // --- Error cases ---

    #[test]
    fn test_streaming_trailing_content() {
        let err = streaming_from_str::<i64>("12 43").unwrap_err();
        assert!(matches!(err, EdnError::TrailingContent { .. }));
    }

    #[test]
    fn test_streaming_empty() {
        let err = streaming_from_str::<i64>("").unwrap_err();
        assert!(matches!(err, EdnError::UnexpectedEof { .. }));
    }
}
