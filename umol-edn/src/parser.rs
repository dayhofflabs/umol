//! Winnow-based EDN parser.

use std::borrow::Cow;

use memchr::memchr2;

use winnow::ascii::digit1;
use winnow::combinator::opt;
use winnow::error::ErrMode;
use winnow::stream::Location;
use winnow::token::{any, one_of, take, take_while};
use winnow::{LocatingSlice, Parser};

use crate::config::{Dialect, DuplicateKeyPolicy, ParseConfig};
use crate::edn::{Edn, EdnMap, EdnSet, Keyword, Symbol};
use crate::error::{unwrap_err, EdnError};

type Input<'a> = LocatingSlice<&'a str>;
type E = ErrMode<EdnError>;
type PResult<T> = Result<T, E>;

/// Get the remaining input as a `&str`.
#[inline(always)]
fn rest<'a>(input: &Input<'a>) -> &'a str {
    *input.as_ref()
}

/// Peek at the next byte without consuming.
#[inline(always)]
fn peek_byte(input: &Input) -> Option<u8> {
    input.as_ref().as_bytes().first().copied()
}

/// Check if a byte is whitespace.
#[inline(always)]
fn is_ws_byte(b: u8) -> bool {
    matches!(b, b' ' | b'\t' | b'\n' | b'\r' | b',')
}

/// Parse a single EDN value from the input, returning the value and remaining input.
pub fn parse_value<'a>(
    input: &'a str,
    config: &ParseConfig,
) -> Result<(Edn<'a>, &'a str), EdnError> {
    let mut located = LocatingSlice::new(input);
    ws_and_comments(&mut located, config).map_err(unwrap_err)?;
    let val = edn_value(config)
        .parse_next(&mut located)
        .map_err(unwrap_err)?;
    let remainder = rest(&located);
    Ok((val, remainder))
}

/// Parse a single EDN value, rejecting trailing non-whitespace content.
pub fn parse_value_strict<'a>(
    input: &'a str,
    config: &ParseConfig,
) -> Result<Edn<'a>, EdnError> {
    let (val, remaining) = parse_value(input, config)?;
    let mut loc = LocatingSlice::new(remaining);
    ws_and_comments(&mut loc, config).ok();
    let after = rest(&loc);
    if !after.is_empty() {
        let trailing_offset = input.len() - after.len();
        return Err(EdnError::TrailingContent {
            offset: trailing_offset,
        });
    }
    Ok(val)
}

/// Parse all EDN values from the input.
pub fn parse_all<'a>(
    input: &'a str,
    config: &ParseConfig,
) -> Result<Vec<Edn<'a>>, EdnError> {
    let mut located = LocatingSlice::new(input);
    let mut values = Vec::new();
    loop {
        ws_and_comments(&mut located, config).ok();
        if peek_byte(&located).is_none() {
            break;
        }
        let val = edn_value(config)
            .parse_next(&mut located)
            .map_err(unwrap_err)?;
        values.push(val);
    }
    Ok(values)
}

// --- Whitespace and comments ---

fn is_ws(c: char) -> bool {
    matches!(c, ' ' | '\t' | '\n' | '\r' | ',')
}

fn ws_and_comments<'a>(input: &mut Input<'a>, config: &ParseConfig) -> PResult<()> {
    loop {
        take_while(0.., is_ws).parse_next(input)?;
        match peek_byte(input) {
            Some(b';') => {
                take_while(0.., |c: char| c != '\n').parse_next(input)?;
                opt('\n').parse_next(input)?;
            }
            Some(b'#') if config.dialect == Dialect::Clojure && rest(input).starts_with("#_") => {
                let _ = "#_".parse_next(input)?;
                ws_and_comments(input, config).ok();
                let _ = edn_value(config).parse_next(input)?;
            }
            _ => break,
        }
    }
    Ok(())
}

// --- Top-level value dispatch ---

fn edn_value<'a, 'b>(
    config: &'b ParseConfig,
) -> impl Parser<Input<'a>, Edn<'a>, E> + 'b
where
    'a: 'b,
{
    move |input: &mut Input<'a>| -> PResult<Edn<'a>> {
        ws_and_comments(input, config).ok();
        let offset = input.current_token_start();
        let b = peek_byte(input)
            .ok_or(ErrMode::Backtrack(EdnError::UnexpectedEof { offset }))?;
        match b {
            b'(' => edn_list(config).parse_next(input),
            b'[' => edn_vector(config).parse_next(input),
            b'{' => edn_map(config).parse_next(input),
            b'#' => edn_dispatch(config).parse_next(input),
            b':' => edn_keyword(input, config.dialect),
            b'"' => edn_string(config.dialect).parse_next(input),
            b'\\' => edn_char(config.dialect).parse_next(input),
            b'+' | b'-' => edn_number_or_symbol(config).parse_next(input),
            b'0'..=b'9' => edn_number(input),
            _ => edn_symbol_or_literal(input),
        }
    }
}

// --- Nil, booleans, symbols ---

pub(crate) fn is_symbol_start(c: char) -> bool {
    matches!(
        c,
        'a'..='z' | 'A'..='Z' | '.' | '*' | '!' | '_' | '?' | '$' | '%' | '&' | '=' | '<' | '>' | '/'
    )
}

pub(crate) fn is_symbol_char(c: char) -> bool {
    is_symbol_start(c) || matches!(c, '0'..='9' | '+' | '-' | '#' | ':' | '\'')
}

fn raw_symbol<'a>(input: &mut Input<'a>) -> PResult<&'a str> {
    (one_of(is_symbol_start), take_while(0.., is_symbol_char))
        .take()
        .parse_next(input)
}

fn edn_symbol_or_literal<'a>(input: &mut Input<'a>) -> PResult<Edn<'a>> {
    let s = raw_symbol(input)?;
    match s {
        "nil" => Ok(Edn::Nil),
        "true" => Ok(Edn::Bool(true)),
        "false" => Ok(Edn::Bool(false)),
        _ => Ok(Edn::Symbol(Symbol::new(s))),
    }
}

// --- Keywords ---

fn edn_keyword<'a>(input: &mut Input<'a>, dialect: Dialect) -> PResult<Edn<'a>> {
    let _ = ':'.parse_next(input)?;
    let s = match dialect {
        Dialect::Clojure => {
            // Clojure allows digit-starting keywords (e.g. :0, :1).
            (one_of(|c: char| is_symbol_start(c) || c.is_ascii_digit()), take_while(0.., is_symbol_char))
                .take()
                .parse_next(input)?
        }
        Dialect::Edn => raw_symbol(input)?,
    };
    Ok(Edn::Keyword(Keyword::new(s)))
}

// --- Numbers ---

fn edn_number<'a>(input: &mut Input<'a>) -> PResult<Edn<'a>> {
    let start = input.current_token_start();
    let num_str: &str = (
        opt(one_of(['+', '-'])),
        digit1,
        opt(('.', digit1)),
        opt((one_of(['e', 'E']), opt(one_of(['+', '-'])), digit1)),
    )
        .take()
        .parse_next(input)?;

    // Check for N/M suffix (bignum)
    if matches!(peek_byte(input), Some(b'N') | Some(b'M')) {
        let _ = any.parse_next(input)?;
        return Err(ErrMode::Cut(EdnError::UnsupportedFeature {
            offset: start,
            feature: "bignum",
        }));
    }

    if memchr::memchr3(b'.', b'e', b'E', num_str.as_bytes()).is_some() {
        let f: f64 = num_str
            .parse()
            .map_err(|_| ErrMode::Cut(EdnError::InvalidNumber { offset: start }))?;
        Ok(Edn::Float(f))
    } else {
        let n: i64 = num_str
            .parse()
            .map_err(|_| ErrMode::Cut(EdnError::InvalidNumber { offset: start }))?;
        Ok(Edn::Int(n))
    }
}

fn edn_number_or_symbol<'a, 'b>(
    _config: &'b ParseConfig,
) -> impl Parser<Input<'a>, Edn<'a>, E> + 'b
where
    'a: 'b,
{
    move |input: &mut Input<'a>| -> PResult<Edn<'a>> {
        // Peek at second byte: if digit, it's a number like +5 or -3.
        let second = input.as_ref().as_bytes().get(1).copied();
        if matches!(second, Some(b'0'..=b'9')) {
            edn_number(input)
        } else {
            edn_symbol_or_literal(input)
        }
    }
}

// --- Strings ---

fn edn_string<'a>(dialect: Dialect) -> impl Parser<Input<'a>, Edn<'a>, E> {
    move |input: &mut Input<'a>| -> PResult<Edn<'a>> {
        let _ = '"'.parse_next(input)?;

        // Fast path: scan for closing quote with no escapes.
        let s = rest(input);
        let bytes = s.as_bytes();
        if let Some(end) = memchr2(b'"', b'\\', bytes) {
            if bytes[end] == b'"' {
                let borrowed: &'a str = take(end).parse_next(input)?;
                let _ = '"'.parse_next(input)?;
                return Ok(Edn::Str(Cow::Borrowed(borrowed)));
            }
        }

        // Slow path: has escapes (or unterminated).
        // Copy text before the first backslash (found by the fast path scan above).
        let pre_escape = memchr::memchr(b'\\', bytes).unwrap_or(0);
        let mut result = String::new();
        if pre_escape > 0 {
            let span: &str = take(pre_escape).parse_next(input)?;
            result.push_str(span);
        }
        loop {
            let offset = input.current_token_start();
            let s = rest(input);
            let c = s
                .chars()
                .next()
                .ok_or(ErrMode::Cut(EdnError::UnexpectedEof { offset }))?;
            match c {
                '"' => {
                    let _ = any.parse_next(input)?;
                    return Ok(Edn::Str(Cow::Owned(result)));
                }
                '\\' => {
                    let esc_offset = input.current_token_start();
                    let _ = any.parse_next(input)?;
                    let s2 = rest(input);
                    let esc = s2
                        .chars()
                        .next()
                        .ok_or(ErrMode::Cut(EdnError::UnexpectedEof { offset: esc_offset }))?;
                    let _ = any.parse_next(input)?;
                    match esc {
                        't' => result.push('\t'),
                        'r' => result.push('\r'),
                        'n' => result.push('\n'),
                        '\\' => result.push('\\'),
                        '"' => result.push('"'),
                        'b' if dialect == Dialect::Clojure => result.push('\u{0008}'),
                        'f' if dialect == Dialect::Clojure => result.push('\u{000C}'),
                        'u' => {
                            let hex: &str =
                                take_while(4..=4, |c: char| c.is_ascii_hexdigit())
                                    .parse_next(input)
                                    .map_err(|_: E| ErrMode::Cut(EdnError::InvalidEscape { offset: esc_offset }))?;
                            let cp = u32::from_str_radix(hex, 16)
                                .map_err(|_| ErrMode::Cut(EdnError::InvalidEscape { offset: esc_offset }))?;
                            let ch = char::from_u32(cp)
                                .ok_or_else(|| ErrMode::Cut(EdnError::InvalidEscape { offset: esc_offset }))?;
                            result.push(ch);
                        }
                        '0'..='7' if dialect == Dialect::Clojure => {
                            let mut val = (esc as u32) - ('0' as u32);
                            for _ in 0..2 {
                                let s3 = rest(input);
                                match s3.chars().next() {
                                    Some(next) if ('0'..='7').contains(&next) => {
                                        val = val * 8 + (next as u32 - '0' as u32);
                                        let _ = any.parse_next(input)?;
                                    }
                                    _ => break,
                                }
                            }
                            if val > 0o377 {
                                return Err(ErrMode::Cut(EdnError::InvalidEscape { offset: esc_offset }));
                            }
                            let ch = char::from_u32(val)
                                .ok_or_else(|| ErrMode::Cut(EdnError::InvalidEscape { offset: esc_offset }))?;
                            result.push(ch);
                        }
                        _ => return Err(ErrMode::Cut(EdnError::InvalidEscape { offset: esc_offset })),
                    }
                    // After escape, batch copy until next " or \.
                    let s = rest(input);
                    let sb = s.as_bytes();
                    let span_end = memchr2(b'"', b'\\', sb).unwrap_or(sb.len());
                    if span_end > 0 {
                        let span: &str = take(span_end).parse_next(input)?;
                        result.push_str(span);
                    }
                }
                _ => {
                    result.push(c);
                    let _ = any.parse_next(input)?;
                }
            }
        }
    }
}

// --- Characters ---

fn edn_char<'a>(dialect: Dialect) -> impl Parser<Input<'a>, Edn<'a>, E> {
    move |input: &mut Input<'a>| -> PResult<Edn<'a>> {
        let _ = '\\'.parse_next(input)?;
        let char_offset = input.current_token_start();
        let s = rest(input);
        let bytes = s.as_bytes();

        // Fast path: single character followed by non-symbol-char or EOF.
        if let Some(&first) = bytes.first() {
            let second = bytes.get(1).copied();
            let single_char = second.is_none()
                || !is_symbol_char(second.unwrap() as char)
                || is_ws_byte(second.unwrap());
            if single_char && first != b'u' {
                let c = s.chars().next().unwrap();
                if c.is_whitespace() {
                    return Err(ErrMode::Cut(EdnError::InvalidCharLiteral { offset: char_offset }));
                }
                let _ = any.parse_next(input)?;
                return Ok(Edn::Char(c));
            }
        } else {
            return Err(ErrMode::Cut(EdnError::UnexpectedEof { offset: char_offset }));
        }

        // Unicode escape \uNNNN
        if bytes[0] == b'u' {
            let _ = 'u'.parse_next(input)?;
            let hex: &str = take_while(4..=4, |c: char| c.is_ascii_hexdigit())
                .parse_next(input)
                .map_err(|_: E| ErrMode::Cut(EdnError::InvalidCharLiteral { offset: char_offset }))?;
            let cp = u32::from_str_radix(hex, 16)
                .map_err(|_| ErrMode::Cut(EdnError::InvalidCharLiteral { offset: char_offset }))?;
            let ch =
                char::from_u32(cp).ok_or_else(|| ErrMode::Cut(EdnError::InvalidCharLiteral { offset: char_offset }))?;
            return Ok(Edn::Char(ch));
        }

        // Named characters (multi-byte sequences like "newline", "return", etc.)
        let named: &[(&str, char)] = if dialect == Dialect::Edn {
            &[
                ("newline", '\n'),
                ("return", '\r'),
                ("space", ' '),
                ("tab", '\t'),
            ]
        } else {
            &[
                ("newline", '\n'),
                ("return", '\r'),
                ("space", ' '),
                ("tab", '\t'),
                ("formfeed", '\u{000C}'),
                ("backspace", '\u{0008}'),
            ]
        };

        for &(name, ch) in named {
            if s.starts_with(name) {
                let after = &s[name.len()..];
                let terminates = after.is_empty()
                    || after.as_bytes().first().map_or(true, |&b| {
                        !is_symbol_char(b as char) || is_ws_byte(b)
                    });
                if terminates {
                    let _ = winnow::token::literal(name).parse_next(input)?;
                    return Ok(Edn::Char(ch));
                }
            }
        }

        // Multi-byte but not a recognized named char — error.
        // (Single chars were handled by the fast path above.)
        let c = s.chars().next().unwrap();
        if c.is_whitespace() {
            return Err(ErrMode::Cut(EdnError::InvalidCharLiteral { offset: char_offset }));
        }
        let _ = any.parse_next(input)?;
        if let Some(&next) = rest(input).as_bytes().first() {
            if is_symbol_char(next as char) && !is_ws_byte(next) {
                return Err(ErrMode::Cut(EdnError::InvalidCharLiteral { offset: char_offset }));
            }
        }
        Ok(Edn::Char(c))
    }
}

// --- Collections ---

fn edn_list<'a, 'b>(
    config: &'b ParseConfig,
) -> impl Parser<Input<'a>, Edn<'a>, E> + 'b
where
    'a: 'b,
{
    move |input: &mut Input<'a>| -> PResult<Edn<'a>> {
        let _ = '('.parse_next(input)?;
        let mut items = Vec::new();
        loop {
            ws_and_comments(input, config).ok();
            match peek_byte(input) {
                None => return Err(ErrMode::Cut(EdnError::UnexpectedEof { offset: input.current_token_start() })),
                Some(b')') => {
                    let _ = ')'.parse_next(input)?;
                    return Ok(Edn::List(items));
                }
                _ => items.push(edn_value(config).parse_next(input)?),
            }
        }
    }
}

fn edn_vector<'a, 'b>(
    config: &'b ParseConfig,
) -> impl Parser<Input<'a>, Edn<'a>, E> + 'b
where
    'a: 'b,
{
    move |input: &mut Input<'a>| -> PResult<Edn<'a>> {
        let _ = '['.parse_next(input)?;
        let mut items = Vec::new();
        loop {
            ws_and_comments(input, config).ok();
            match peek_byte(input) {
                None => return Err(ErrMode::Cut(EdnError::UnexpectedEof { offset: input.current_token_start() })),
                Some(b']') => {
                    let _ = ']'.parse_next(input)?;
                    return Ok(Edn::Vector(items));
                }
                _ => items.push(edn_value(config).parse_next(input)?),
            }
        }
    }
}

fn edn_map<'a, 'b>(
    config: &'b ParseConfig,
) -> impl Parser<Input<'a>, Edn<'a>, E> + 'b
where
    'a: 'b,
{
    move |input: &mut Input<'a>| -> PResult<Edn<'a>> {
        let _ = '{'.parse_next(input)?;
        let mut map = EdnMap::new();
        loop {
            ws_and_comments(input, config).ok();
            match peek_byte(input) {
                None => return Err(ErrMode::Cut(EdnError::UnexpectedEof { offset: input.current_token_start() })),
                Some(b'}') => {
                    let _ = '}'.parse_next(input)?;
                    return Ok(Edn::Map(map));
                }
                _ => {
                    let key_offset = input.current_token_start();
                    let key = edn_value(config).parse_next(input)?;
                    ws_and_comments(input, config).ok();
                    let val = edn_value(config).parse_next(input)?;
                    if config.duplicate_keys == DuplicateKeyPolicy::Error && map.contains_key(&key) {
                        return Err(ErrMode::Cut(EdnError::DuplicateKey { offset: key_offset }));
                    }
                    map.insert(key, val);
                }
            }
        }
    }
}

fn edn_set<'a, 'b>(
    config: &'b ParseConfig,
) -> impl Parser<Input<'a>, Edn<'a>, E> + 'b
where
    'a: 'b,
{
    move |input: &mut Input<'a>| -> PResult<Edn<'a>> {
        let _ = '{'.parse_next(input)?;
        let mut set = EdnSet::new();
        loop {
            ws_and_comments(input, config).ok();
            match peek_byte(input) {
                None => return Err(ErrMode::Cut(EdnError::UnexpectedEof { offset: input.current_token_start() })),
                Some(b'}') => {
                    let _ = '}'.parse_next(input)?;
                    return Ok(Edn::Set(set));
                }
                _ => { set.insert(edn_value(config).parse_next(input)?); }
            }
        }
    }
}

// --- # dispatch ---

fn edn_dispatch<'a, 'b>(
    config: &'b ParseConfig,
) -> impl Parser<Input<'a>, Edn<'a>, E> + 'b
where
    'a: 'b,
{
    move |input: &mut Input<'a>| -> PResult<Edn<'a>> {
        let _ = '#'.parse_next(input)?;
        let b = peek_byte(input)
            .ok_or(ErrMode::Cut(EdnError::UnexpectedEof { offset: input.current_token_start() }))?;
        match b {
            b'{' => edn_set(config).parse_next(input),
            b'_' if config.dialect == Dialect::Clojure => {
                let _ = '_'.parse_next(input)?;
                ws_and_comments(input, config).ok();
                let _ = edn_value(config).parse_next(input)?;
                edn_value(config).parse_next(input)
            }
            b'#' if config.dialect == Dialect::Clojure => {
                let _ = '#'.parse_next(input)?;
                match peek_byte(input) {
                    Some(b'N') => { let _ = "NaN".parse_next(input)?; Ok(Edn::Float(f64::NAN)) }
                    Some(b'-') => { let _ = "-Inf".parse_next(input)?; Ok(Edn::Float(f64::NEG_INFINITY)) }
                    Some(b'I') => { let _ = "Inf".parse_next(input)?; Ok(Edn::Float(f64::INFINITY)) }
                    _ => {
                        let offset = input.current_token_start();
                        let found = rest(input).chars().next().unwrap_or('\0');
                        Err(ErrMode::Cut(EdnError::UnexpectedToken { offset, found }))
                    }
                }
            }
            _ => {
                let tag = raw_symbol(input)?;
                ws_and_comments(input, config).ok();
                let val = edn_value(config).parse_next(input)?;
                Ok(Edn::Tagged(tag.to_string(), Box::new(val)))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;
    use std::borrow::Cow;
    use crate::edn::{Edn, EdnMap, EdnSet, Keyword};
    use crate::error::EdnError;
    use crate::reader::{read_all, read_string, read_string_with, Reader};
    use crate::config::{DuplicateKeyPolicy, ParseConfig, Dialect};

    // --- Primitives ---

    #[rstest]
    #[case("nil", Edn::Nil)]
    #[case("true", Edn::Bool(true))]
    #[case("false", Edn::Bool(false))]
    fn test_read_string_literals(#[case] input: &str, #[case] expected: Edn<'_>) {
        assert_eq!(read_string(input).unwrap(), expected);
    }

    #[rstest]
    #[case("0", 0)]
    #[case("12", 12)]
    #[case("-1", -1)]
    #[case("+5", 5)]
    #[case("9223372036854775807", i64::MAX)]
    #[case("-9223372036854775808", i64::MIN)]
    fn test_read_string_int(#[case] input: &str, #[case] expected: i64) {
        assert_eq!(read_string(input).unwrap(), Edn::Int(expected));
    }

    #[rstest]
    #[case("1.0", 1.0)]
    #[case("-3.14", -3.14)]
    #[case("1e10", 1e10)]
    #[case("1.5e-3", 1.5e-3)]
    #[case("1E10", 1e10)]
    fn test_read_string_float(#[case] input: &str, #[case] expected: f64) {
        assert_eq!(read_string(input).unwrap(), Edn::Float(expected));
    }

    #[rstest]
    #[case("##NaN")]
    fn test_read_string_nan(#[case] input: &str) {
        match read_string(input).unwrap() {
            Edn::Float(f) => assert!(f.is_nan()),
            other => panic!("expected Float(NaN), got {other:?}"),
        }
    }

    #[rstest]
    #[case("##Inf", f64::INFINITY)]
    #[case("##-Inf", f64::NEG_INFINITY)]
    fn test_read_string_special_float(#[case] input: &str, #[case] expected: f64) {
        assert_eq!(read_string(input).unwrap(), Edn::Float(expected));
    }

    // --- Strings ---

    #[rstest]
    #[case(r#""""#, "")]
    #[case(r#""hello""#, "hello")]
    #[case(r#""hello world""#, "hello world")]
    #[case(r#""line\nbreak""#, "line\nbreak")]
    #[case(r#""tab\there""#, "tab\there")]
    #[case(r#""quote\"here""#, "quote\"here")]
    #[case(r#""back\\slash""#, "back\\slash")]
    #[case(r#""cr\rhere""#, "cr\rhere")]
    fn test_read_string_str(#[case] input: &str, #[case] expected: &str) {
        assert_eq!(
            read_string(input).unwrap(),
            Edn::Str(Cow::Owned(expected.to_string()))
        );
    }

    #[rstest]
    #[case(r#""\u0041""#, "A")]
    #[case(r#""\u03BB""#, "\u{03BB}")]
    fn test_read_string_unicode_escape(#[case] input: &str, #[case] expected: &str) {
        assert_eq!(
            read_string(input).unwrap(),
            Edn::Str(Cow::Owned(expected.to_string()))
        );
    }

    // --- Characters ---

    #[rstest]
    #[case("\\a", 'a')]
    #[case("\\Z", 'Z')]
    #[case("\\newline", '\n')]
    #[case("\\return", '\r')]
    #[case("\\space", ' ')]
    #[case("\\tab", '\t')]
    #[case("\\u0041", 'A')]
    fn test_read_string_char(#[case] input: &str, #[case] expected: char) {
        assert_eq!(read_string(input).unwrap(), Edn::Char(expected));
    }

    // --- Keywords ---

    #[rstest]
    #[case(":foo", "foo")]
    #[case(":ns/name", "ns/name")]
    #[case(":a.b/c", "a.b/c")]
    fn test_read_string_keyword(#[case] input: &str, #[case] expected_name: &str) {
        let val = read_string(input).unwrap();
        match &val {
            Edn::Keyword(k) => assert_eq!(k.as_str(), expected_name),
            other => panic!("expected Keyword, got {other:?}"),
        }
    }

    #[test]
    fn test_keyword_namespace() {
        let k = Keyword::new("ns/name");
        assert_eq!(k.namespace(), Some("ns"));
        assert_eq!(k.name(), "name");

        let k2 = Keyword::new("bare");
        assert_eq!(k2.namespace(), None);
        assert_eq!(k2.name(), "bare");
    }

    // --- Symbols ---

    #[rstest]
    #[case("foo", "foo")]
    #[case("ns/name", "ns/name")]
    #[case("my.ns/sym", "my.ns/sym")]
    fn test_read_string_symbol(#[case] input: &str, #[case] expected_name: &str) {
        let val = read_string(input).unwrap();
        match &val {
            Edn::Symbol(s) => assert_eq!(s.as_str(), expected_name),
            other => panic!("expected Symbol, got {other:?}"),
        }
    }

    // --- Collections ---

    #[rstest]
    #[case("()", Edn::List(vec![]))]
    #[case("(1 2 3)", Edn::List(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)]))]
    #[case("(nil true)", Edn::List(vec![Edn::Nil, Edn::Bool(true)]))]
    fn test_read_string_list(#[case] input: &str, #[case] expected: Edn<'_>) {
        assert_eq!(read_string(input).unwrap(), expected);
    }

    #[rstest]
    #[case("[]", Edn::Vector(vec![]))]
    #[case("[1 2 3]", Edn::Vector(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)]))]
    fn test_read_string_vector(#[case] input: &str, #[case] expected: Edn<'_>) {
        assert_eq!(read_string(input).unwrap(), expected);
    }

    #[test]
    fn test_read_string_map() {
        let val = read_string("{:a 1 :b 2}").unwrap();
        let mut expected = EdnMap::new();
        expected.insert(Edn::Keyword(Keyword::new("a")), Edn::Int(1));
        expected.insert(Edn::Keyword(Keyword::new("b")), Edn::Int(2));
        assert_eq!(val, Edn::Map(expected));
    }

    #[test]
    fn test_read_string_empty_map() {
        assert_eq!(read_string("{}").unwrap(), Edn::Map(EdnMap::new()));
    }

    #[test]
    fn test_read_string_set() {
        let val = read_string("#{1 2 3}").unwrap();
        let mut expected = EdnSet::new();
        expected.insert(Edn::Int(1));
        expected.insert(Edn::Int(2));
        expected.insert(Edn::Int(3));
        assert_eq!(val, Edn::Set(expected));
    }

    #[test]
    fn test_read_string_empty_set() {
        assert_eq!(read_string("#{}").unwrap(), Edn::Set(EdnSet::new()));
    }

    // --- Nested ---

    #[test]
    fn test_read_string_nested() {
        let val = read_string("{:items [1 (2 3)] :flag true}").unwrap();
        let mut expected = EdnMap::new();
        expected.insert(
            Edn::Keyword(Keyword::new("items")),
            Edn::Vector(vec![
                Edn::Int(1),
                Edn::List(vec![Edn::Int(2), Edn::Int(3)]),
            ]),
        );
        expected.insert(Edn::Keyword(Keyword::new("flag")), Edn::Bool(true));
        assert_eq!(val, Edn::Map(expected));
    }

    // --- Tagged literals ---

    #[test]
    fn test_read_string_tagged() {
        let val = read_string("#myapp/Person {:name \"Alice\"}").unwrap();
        match val {
            Edn::Tagged(tag, inner) => {
                assert_eq!(tag, "myapp/Person");
                assert!(inner.is_map());
            }
            other => panic!("expected Tagged, got {other:?}"),
        }
    }

    // --- Comments ---

    #[rstest]
    #[case("; comment\n12", Edn::Int(12))]
    #[case("12 ; trailing", Edn::Int(12))]
    fn test_read_string_comment(#[case] input: &str, #[case] expected: Edn<'_>) {
        assert_eq!(read_string(input).unwrap(), expected);
    }

    // --- Discard ---

    #[rstest]
    #[case("#_ foo 12", Edn::Int(12))]
    #[case("[1 #_ 2 3]", Edn::Vector(vec![Edn::Int(1), Edn::Int(3)]))]
    fn test_read_string_discard(#[case] input: &str, #[case] expected: Edn<'_>) {
        assert_eq!(read_string(input).unwrap(), expected);
    }

    // --- Whitespace variants ---

    #[rstest]
    #[case("  12  ", Edn::Int(12))]
    #[case("\t12\n", Edn::Int(12))]
    #[case(",12,", Edn::Int(12))]
    #[case("[1,,2,,3]", Edn::Vector(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)]))]
    fn test_read_string_whitespace(#[case] input: &str, #[case] expected: Edn<'_>) {
        assert_eq!(read_string(input).unwrap(), expected);
    }

    // --- Error cases ---

    #[rstest]
    #[case("", EdnError::UnexpectedEof { offset: 0 })]
    #[case("   ", EdnError::UnexpectedEof { offset: 3 })]
    fn test_read_string_error_empty(#[case] input: &str, #[case] expected: EdnError) {
        assert_eq!(read_string(input).unwrap_err(), expected);
    }

    #[test]
    fn test_read_string_error_trailing() {
        assert_eq!(
            read_string("12 43").unwrap_err(),
            EdnError::TrailingContent { offset: 3 },
        );
    }

    #[test]
    fn test_read_string_error_duplicate_key() {
        assert_eq!(
            read_string("{:a 1 :a 2}").unwrap_err(),
            EdnError::DuplicateKey { offset: 6 },
        );
    }

    #[rstest]
    #[case(r#""\q""#, EdnError::InvalidEscape { offset: 1 })]
    #[case(r#""\u000G""#, EdnError::InvalidEscape { offset: 1 })]
    fn test_read_string_error_invalid_escape(#[case] input: &str, #[case] expected: EdnError) {
        assert_eq!(read_string(input).unwrap_err(), expected);
    }

    #[rstest]
    #[case(r"\abc", EdnError::InvalidCharLiteral { offset: 1 })]
    fn test_read_string_error_invalid_char_literal(#[case] input: &str, #[case] expected: EdnError) {
        assert_eq!(read_string(input).unwrap_err(), expected);
    }

    #[test]
    fn test_read_string_error_unexpected_token() {
        let err = read_string("##xyz").unwrap_err();
        assert!(
            matches!(err, EdnError::UnexpectedToken { offset: 2, found: 'x' }),
            "expected UnexpectedToken, got {err:?}"
        );
    }

    #[rstest]
    #[case("12N", EdnError::UnsupportedFeature { offset: 0, feature: "bignum" })]
    #[case("12M", EdnError::UnsupportedFeature { offset: 0, feature: "bignum" })]
    fn test_read_string_error_unsupported_feature(#[case] input: &str, #[case] expected: EdnError) {
        assert_eq!(read_string(input).unwrap_err(), expected);
    }

    #[test]
    fn test_read_string_error_unclosed_string() {
        let err = read_string(r#""hello"#).unwrap_err();
        assert!(
            matches!(err, EdnError::UnexpectedEof { .. }),
            "expected UnexpectedEof, got {err:?}"
        );
    }

    #[test]
    fn test_read_string_error_unclosed_vector() {
        let err = read_string("[1 2").unwrap_err();
        assert!(
            matches!(err, EdnError::UnexpectedEof { .. }),
            "expected UnexpectedEof, got {err:?}"
        );
    }

    #[test]
    fn test_read_string_duplicate_key_last_wins() {
        let config = ParseConfig {
            duplicate_keys: DuplicateKeyPolicy::LastWins,
            ..Default::default()
        };
        let val = read_string_with("{:a 1 :a 2}", &config).unwrap();
        let mut expected = EdnMap::new();
        expected.insert(Edn::Keyword(Keyword::new("a")), Edn::Int(2));
        assert_eq!(val, Edn::Map(expected));
    }

    // --- Integer overflow ---

    #[test]
    fn test_read_string_error_integer_overflow() {
        let err = read_string("99999999999999999999").unwrap_err();
        assert!(
            matches!(err, EdnError::InvalidNumber { offset: 0 }),
            "expected InvalidNumber, got {err:?}"
        );
    }

    // --- read_all ---

    #[test]
    fn test_read_all() {
        let values = read_all("1 2 3").unwrap();
        assert_eq!(values, vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)]);
    }

    #[test]
    fn test_read_all_empty() {
        let values = read_all("").unwrap();
        assert!(values.is_empty());
    }

    #[test]
    fn test_read_all_mixed() {
        let values = read_all(":a [1 2] nil").unwrap();
        assert_eq!(values.len(), 3);
        assert!(values[0].is_keyword());
        assert!(values[1].is_vector());
        assert!(values[2].is_nil());
    }

    // --- Reader iterator ---

    #[test]
    fn test_reader_iterator() {
        let reader = Reader::new("1 :foo [3]");
        let values: Result<Vec<_>, _> = reader.collect();
        let values = values.unwrap();
        assert_eq!(values.len(), 3);
        assert_eq!(values[0], Edn::Int(1));
        assert!(values[1].is_keyword());
        assert!(values[2].is_vector());
    }

    #[test]
    fn test_reader_empty() {
        let reader = Reader::new("");
        let values: Vec<_> = reader.collect();
        assert!(values.is_empty());
    }

    // --- Round-trip ---

    #[rstest]
    #[case("nil")]
    #[case("true")]
    #[case("false")]
    #[case("12")]
    #[case("-1")]
    #[case("1.5")]
    #[case(":keyword")]
    #[case(":ns/name")]
    #[case("symbol")]
    #[case("()")]
    #[case("(1 2 3)")]
    #[case("[]")]
    #[case("[1 2 3]")]
    #[case("{}")]
    #[case("\\a")]
    #[case("\\newline")]
    #[case("\\space")]
    fn test_roundtrip(#[case] input: &str) {
        let val = read_string(input).unwrap();
        let formatted = val.to_string();
        let reparsed = read_string(&formatted).unwrap();
        assert_eq!(val, reparsed);
    }

    #[test]
    fn test_roundtrip_string() {
        let val = read_string(r#""hello\nworld""#).unwrap();
        let formatted = val.to_string();
        let reparsed = read_string(&formatted).unwrap();
        assert_eq!(val, reparsed);
    }

    #[test]
    fn test_roundtrip_map() {
        let val = read_string("{:a 1, :b 2}").unwrap();
        let formatted = val.to_string();
        let reparsed = read_string(&formatted).unwrap();
        assert_eq!(val, reparsed);
    }

    #[test]
    fn test_roundtrip_set() {
        let val = read_string("#{1 2 3}").unwrap();
        let formatted = val.to_string();
        let reparsed = read_string(&formatted).unwrap();
        assert_eq!(val, reparsed);
    }

    #[test]
    fn test_roundtrip_special_floats() {
        for input in &["##Inf", "##-Inf"] {
            let val = read_string(input).unwrap();
            let formatted = val.to_string();
            let reparsed = read_string(&formatted).unwrap();
            assert_eq!(val, reparsed);
        }
    }

    // --- Edn accessors ---

    #[test]
    fn test_edn_get() {
        let val = read_string("{:name \"Alice\" :age 30}").unwrap();
        assert_eq!(val.get("name").unwrap().as_str(), Some("Alice"));
        assert_eq!(val.get("age").unwrap().as_i64(), Some(30));
        assert!(val.get("missing").is_none());
    }

    #[test]
    fn test_edn_numeric_narrowing() {
        let val = Edn::Int(12);
        assert_eq!(val.as_u8(), Some(12));
        assert_eq!(val.as_u16(), Some(12));
        assert_eq!(val.as_u32(), Some(12));
        assert_eq!(val.as_i32(), Some(12));

        let val_neg = Edn::Int(-1);
        assert_eq!(val_neg.as_u8(), None);
        assert_eq!(val_neg.as_i8(), Some(-1));
    }

    #[test]
    fn test_edn_iter() {
        let val = read_string("[1 2 3]").unwrap();
        let items: Vec<_> = val.iter().collect();
        assert_eq!(items.len(), 3);
    }

    // --- Clojure-compatible escapes (non-strict mode) ---

    #[test]
    fn test_read_string_backspace_formfeed() {
        let val = read_string(r#""\b\f""#).unwrap();
        assert_eq!(
            val,
            Edn::Str(Cow::Owned("\u{0008}\u{000C}".to_string()))
        );
    }

    #[test]
    fn test_read_string_char_formfeed_backspace() {
        assert_eq!(read_string("\\formfeed").unwrap(), Edn::Char('\u{000C}'));
        assert_eq!(read_string("\\backspace").unwrap(), Edn::Char('\u{0008}'));
    }

    // --- Edn strict dialect ---

    fn edn_strict() -> ParseConfig {
        ParseConfig {
            dialect: Dialect::Edn,
            ..Default::default()
        }
    }

    #[rstest]
    #[case("##NaN")]
    #[case("##Inf")]
    #[case("##-Inf")]
    fn test_edn_strict_rejects_special_floats(#[case] input: &str) {
        assert!(read_string_with(input, &edn_strict()).is_err());
    }

    #[test]
    fn test_edn_strict_no_discard() {
        // In Edn mode, #_ is not discard — it's a tagged literal with tag "_".
        // So "#_ foo 12" parses as Tagged("_", Symbol("foo")), with trailing "12".
        assert!(read_string_with("#_ foo 12", &edn_strict()).is_err());
        // Inside a vector, #_ parses as a tagged value.
        let val = read_string_with("[1 #_ 2 3]", &edn_strict()).unwrap();
        assert_eq!(
            val,
            Edn::Vector(vec![
                Edn::Int(1),
                Edn::Tagged("_".to_string(), Box::new(Edn::Int(2))),
                Edn::Int(3),
            ])
        );
    }

    #[rstest]
    #[case(":0")]
    #[case(":1")]
    #[case(":123")]
    fn test_edn_strict_rejects_digit_keywords(#[case] input: &str) {
        assert!(read_string_with(input, &edn_strict()).is_err());
    }

    #[rstest]
    #[case(":0", "0")]
    #[case(":1", "1")]
    #[case(":123abc", "123abc")]
    fn test_clojure_allows_digit_keywords(#[case] input: &str, #[case] expected: &str) {
        let val = read_string(input).unwrap();
        match &val {
            Edn::Keyword(k) => assert_eq!(k.as_str(), expected),
            other => panic!("expected Keyword, got {other:?}"),
        }
    }
}
