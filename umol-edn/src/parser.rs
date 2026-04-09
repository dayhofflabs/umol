//  Winnow-based EDN parser.

use std::borrow::Cow;
use std::cell::Cell;

use memchr::memchr2;
use winnow::ascii::digit1;
use winnow::combinator::opt;
use winnow::error::ErrMode;
use winnow::stream::{Location, Stream};
use winnow::token::{any, literal, one_of, take_while};
use winnow::{LocatingSlice, Parser};

use crate::collections::{EdnMap, EdnSeq, EdnSet};
use crate::config::{DuplicateKeyPolicy, ParseConfig};
#[cfg(feature = "bignum")]
use crate::edn::EdnBigDecimal;
use crate::edn::{Edn, EdnKeyword, EdnSymbol};
use crate::error::{unwrap_err, ParseError};

type Input<'a> = LocatingSlice<&'a str>;
type E = ErrMode<ParseError>;
type PResult<T> = Result<T, E>;

const MAX_DEPTH: u16 = 128;

struct ParseCtx<'c> {
    config: &'c ParseConfig,
    depth: Cell<u16>,
}

impl<'c> ParseCtx<'c> {
    fn new(config: &'c ParseConfig) -> Self {
        Self {
            config,
            depth: Cell::new(MAX_DEPTH),
        }
    }

    fn enter_scope(&self, offset: usize) -> PResult<()> {
        let d = self
            .depth
            .get()
            .checked_sub(1)
            .ok_or(ErrMode::Cut(ParseError::RecursionLimit { offset }))?;
        self.depth.set(d);
        Ok(())
    }

    fn leave_scope(&self) {
        self.depth.set(self.depth.get() + 1);
    }
}

/// Get the remaining input as a `&str`.
#[inline(always)]
fn rest<'a>(input: &Input<'a>) -> &'a str {
    input.as_ref()
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
) -> Result<(Edn<'a>, &'a str), ParseError> {
    let ctx = ParseCtx::new(config);
    let mut located = LocatingSlice::new(input);
    ws_and_comments(&mut located, &ctx);
    let val = edn_value(&ctx)
        .parse_next(&mut located)
        .map_err(unwrap_err)?;
    let remainder = rest(&located);
    Ok((val, remainder))
}

/// Parse a single EDN value, rejecting trailing non-whitespace content.
pub fn parse_value_strict<'a>(input: &'a str, config: &ParseConfig) -> Result<Edn<'a>, ParseError> {
    let (val, remaining) = parse_value(input, config)?;
    let ctx = ParseCtx::new(config);
    let mut loc = LocatingSlice::new(remaining);
    ws_and_comments(&mut loc, &ctx);
    let after = rest(&loc);
    if !after.is_empty() {
        let trailing_offset = input.len() - after.len();
        return Err(ParseError::TrailingContent {
            offset: trailing_offset,
        });
    }
    Ok(val)
}

/// Parse all EDN values from the input.
pub fn parse_all<'a>(input: &'a str, config: &ParseConfig) -> Result<Vec<Edn<'a>>, ParseError> {
    let ctx = ParseCtx::new(config);
    let mut located = LocatingSlice::new(input);
    let mut values = Vec::new();
    loop {
        ws_and_comments(&mut located, &ctx);
        if peek_byte(&located).is_none() {
            break;
        }
        let val = edn_value(&ctx)
            .parse_next(&mut located)
            .map_err(unwrap_err)?;
        values.push(val);
    }
    Ok(values)
}

fn is_ws(c: char) -> bool {
    matches!(c, ' ' | '\t' | '\n' | '\r' | ',')
}

fn ws_and_comments<'a>(input: &mut Input<'a>, ctx: &ParseCtx<'_>) {
    loop {
        let _: PResult<_> = take_while(0.., is_ws).parse_next(input);
        match peek_byte(input) {
            Some(b';') => {
                let _: PResult<_> = take_while(0.., |c: char| c != '\n').parse_next(input);
                let _: PResult<_> = opt('\n').parse_next(input);
            }
            Some(b'#') if rest(input).starts_with("#_") => {
                let _: PResult<_> = "#_".parse_next(input);
                ws_and_comments(input, ctx);
                if edn_value(ctx).parse_next(input).is_err() {
                    break;
                }
            }
            _ => break,
        }
    }
}

fn edn_value<'a, 'b>(ctx: &'b ParseCtx<'_>) -> impl Parser<Input<'a>, Edn<'a>, E> + 'b
where
    'a: 'b,
{
    move |input: &mut Input<'a>| -> PResult<Edn<'a>> {
        ws_and_comments(input, ctx);
        edn_value_dispatch(input, ctx)
    }
}

/// Dispatch to the correct value parser. Caller must skip whitespace first.
#[inline(always)]
fn edn_value_dispatch<'a>(input: &mut Input<'a>, ctx: &ParseCtx<'_>) -> PResult<Edn<'a>> {
    let offset = input.current_token_start();
    let b = peek_byte(input).ok_or(ErrMode::Backtrack(ParseError::UnexpectedEof { offset }))?;
    match b {
        b'(' => edn_list(ctx).parse_next(input),
        b'[' => edn_vector(ctx).parse_next(input),
        b'{' => edn_map(ctx).parse_next(input),
        b'#' => edn_dispatch(ctx).parse_next(input),
        b':' => edn_keyword(input),
        b'"' => edn_string().parse_next(input),
        b'\\' => edn_char().parse_next(input),
        b'+' | b'-' => edn_number_or_symbol().parse_next(input),
        b'0'..=b'9' => edn_number(input),
        _ => edn_symbol_or_literal(input),
    }
}

pub(crate) fn is_symbol_start(c: char) -> bool {
    matches!(
        c,
        'a'..='z' | 'A'..='Z' | '.' | '*' | '+' | '!' | '-' | '_' | '?' | '$' | '%' | '&' | '=' | '<' | '>' | '/'
    )
}

pub(crate) fn is_symbol_char(c: char) -> bool {
    is_symbol_start(c) || matches!(c, '0'..='9' | '#' | ':' | '\'')
}

#[inline]
fn raw_symbol<'a>(input: &mut Input<'a>) -> PResult<&'a str> {
    let s: &str = (one_of(is_symbol_start), take_while(0.., is_symbol_char))
        .take()
        .parse_next(input)?;
    let start = input.current_token_start() - s.len();
    validate_symbol(s, start).map_err(ErrMode::Cut)?;
    Ok(s)
}

#[inline]
pub(crate) fn validate_symbol(s: &str, offset: usize) -> Result<(), ParseError> {
    if s == "/" {
        return Ok(());
    }
    if let Some(slash_pos) = s.find('/') {
        let prefix = &s[..slash_pos];
        let name = &s[slash_pos + 1..];
        if prefix.is_empty() || name.is_empty() {
            return Err(ParseError::InvalidSymbol { offset });
        }
        let first_name_char = name.chars().next().unwrap();
        if first_name_char.is_ascii_digit()
            || name.contains('/')
            || !is_symbol_start(first_name_char)
        {
            return Err(ParseError::InvalidSymbol { offset });
        }
    }
    Ok(())
}

#[inline]
fn edn_symbol_or_literal<'a>(input: &mut Input<'a>) -> PResult<Edn<'a>> {
    let s = raw_symbol(input)?;
    match s {
        "nil" => Ok(Edn::Nil),
        "true" => Ok(Edn::Bool(true)),
        "false" => Ok(Edn::Bool(false)),
        _ => Ok(Edn::Symbol(EdnSymbol::new(s))),
    }
}

#[inline]
fn edn_keyword<'a>(input: &mut Input<'a>) -> PResult<Edn<'a>> {
    let _ = ':'.parse_next(input)?;
    let start = input.current_token_start() - 1;
    let s: &str = (
        one_of(|c: char| is_symbol_start(c) && c != '/'),
        take_while(0.., is_symbol_char),
    )
        .take()
        .parse_next(input)?;
    validate_symbol(s, start).map_err(ErrMode::Cut)?;
    Ok(Edn::Keyword(EdnKeyword::new(s)))
}

#[inline]
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
        let suffix = any.parse_next(input)?;

        #[cfg(not(feature = "bignum"))]
        {
            let _ = suffix;
            return Err(ErrMode::Cut(ParseError::UnsupportedFeature {
                offset: start,
                feature: "bignum",
            }));
        }

        #[cfg(feature = "bignum")]
        {
            let is_float = memchr::memchr3(b'.', b'e', b'E', num_str.as_bytes()).is_some();
            return if suffix == 'N' {
                if is_float {
                    Err(ErrMode::Cut(ParseError::InvalidNumber { offset: start }))
                } else {
                    let n: num_bigint::BigInt = num_str
                        .parse()
                        .map_err(|_| ErrMode::Cut(ParseError::InvalidNumber { offset: start }))?;
                    Ok(Edn::BigInt(n))
                }
            } else {
                let d: bigdecimal::BigDecimal = num_str
                    .parse()
                    .map_err(|_| ErrMode::Cut(ParseError::InvalidNumber { offset: start }))?;
                Ok(Edn::BigDecimal(EdnBigDecimal::new(d)))
            };
        }
    }

    // Reject leading zeros (007 is invalid, 0 is fine).
    let digits = num_str.trim_start_matches(['+', '-']);
    if digits.len() > 1 && digits.starts_with('0') && !digits.starts_with("0.") {
        return Err(ErrMode::Cut(ParseError::InvalidNumber { offset: start }));
    }

    if memchr::memchr3(b'.', b'e', b'E', num_str.as_bytes()).is_some() {
        let f: f64 = num_str
            .parse()
            .map_err(|_| ErrMode::Cut(ParseError::InvalidNumber { offset: start }))?;
        if !f.is_finite() {
            return Err(ErrMode::Cut(ParseError::InvalidNumber { offset: start }));
        }
        Ok(Edn::Float(f))
    } else {
        match num_str.parse::<i64>() {
            Ok(n) => Ok(Edn::Int(n)),
            #[cfg(feature = "bignum")]
            Err(_) => {
                let n: num_bigint::BigInt = num_str
                    .parse()
                    .map_err(|_| ErrMode::Cut(ParseError::InvalidNumber { offset: start }))?;
                Ok(Edn::BigInt(n))
            }
            #[cfg(not(feature = "bignum"))]
            Err(_) => Err(ErrMode::Cut(ParseError::InvalidNumber { offset: start })),
        }
    }
}

fn edn_number_or_symbol<'a>() -> impl Parser<Input<'a>, Edn<'a>, E> {
    move |input: &mut Input<'a>| -> PResult<Edn<'a>> {
        let second = input.as_ref().as_bytes().get(1).copied();
        if matches!(second, Some(b'0'..=b'9')) {
            edn_number(input)
        } else {
            edn_symbol_or_literal(input)
        }
    }
}

fn edn_string<'a>() -> impl Parser<Input<'a>, Edn<'a>, E> {
    move |input: &mut Input<'a>| -> PResult<Edn<'a>> {
        let _ = '"'.parse_next(input)?;

        // Fast path: scan for closing quote with no escapes.
        // memchr returns byte offsets; use next_slice (byte-based) instead of
        // take (character-based) to advance correctly on multi-byte UTF-8.
        let s = rest(input);
        let bytes = s.as_bytes();
        if let Some(end) = memchr2(b'"', b'\\', bytes) {
            if bytes[end] == b'"' {
                let borrowed: &'a str = input.next_slice(end);
                let _ = '"'.parse_next(input)?;
                return Ok(Edn::Str(Cow::Borrowed(borrowed)));
            }
        }

        // Slow path: has escapes (or unterminated).
        let pre_escape = memchr::memchr(b'\\', bytes).unwrap_or(0);
        let mut result = String::new();
        if pre_escape > 0 {
            let span: &str = input.next_slice(pre_escape);
            result.push_str(span);
        }
        loop {
            let offset = input.current_token_start();
            let s = rest(input);
            let c = s
                .chars()
                .next()
                .ok_or(ErrMode::Cut(ParseError::UnexpectedEof { offset }))?;
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
                        .ok_or(ErrMode::Cut(ParseError::UnexpectedEof {
                            offset: esc_offset,
                        }))?;
                    let _ = any.parse_next(input)?;
                    match esc {
                        't' => result.push('\t'),
                        'r' => result.push('\r'),
                        'n' => result.push('\n'),
                        '\\' => result.push('\\'),
                        '"' => result.push('"'),
                        'u' => {
                            let hex: &str = take_while(4..=4, |c: char| c.is_ascii_hexdigit())
                                .parse_next(input)
                                .map_err(|_: E| {
                                    ErrMode::Cut(ParseError::InvalidEscape { offset: esc_offset })
                                })?;
                            let cp = u32::from_str_radix(hex, 16).map_err(|_| {
                                ErrMode::Cut(ParseError::InvalidEscape { offset: esc_offset })
                            })?;
                            let ch = char::from_u32(cp).ok_or(ErrMode::Cut(
                                ParseError::InvalidEscape { offset: esc_offset },
                            ))?;
                            result.push(ch);
                        }
                        _ => {
                            return Err(ErrMode::Cut(ParseError::InvalidEscape {
                                offset: esc_offset,
                            }))
                        }
                    }
                    let s = rest(input);
                    let sb = s.as_bytes();
                    let span_end = memchr2(b'"', b'\\', sb).unwrap_or(sb.len());
                    if span_end > 0 {
                        let span: &str = input.next_slice(span_end);
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

fn edn_char<'a>() -> impl Parser<Input<'a>, Edn<'a>, E> {
    move |input: &mut Input<'a>| -> PResult<Edn<'a>> {
        let _ = '\\'.parse_next(input)?;
        let char_offset = input.current_token_start();
        let s = rest(input);
        let bytes = s.as_bytes();

        // Fast path: single character followed by non-symbol-char or EOF.
        if let Some(&_first) = bytes.first() {
            let second = bytes.get(1).copied();
            let single_char = second.is_none()
                || !is_symbol_char(second.unwrap() as char)
                || is_ws_byte(second.unwrap());
            if single_char {
                let c = s.chars().next().unwrap();
                if c.is_whitespace() {
                    return Err(ErrMode::Cut(ParseError::InvalidCharLiteral {
                        offset: char_offset,
                    }));
                }
                let _ = any.parse_next(input)?;
                return Ok(Edn::Char(c));
            }
        } else {
            return Err(ErrMode::Cut(ParseError::UnexpectedEof {
                offset: char_offset,
            }));
        }

        // Unicode escape \uNNNN — only reached when \u is followed by symbol chars.
        if bytes[0] == b'u' {
            let _ = 'u'.parse_next(input)?;
            let hex: &str = take_while(4..=4, |c: char| c.is_ascii_hexdigit())
                .parse_next(input)
                .map_err(|_: E| {
                    ErrMode::Cut(ParseError::InvalidCharLiteral {
                        offset: char_offset,
                    })
                })?;
            let cp = u32::from_str_radix(hex, 16).map_err(|_| {
                ErrMode::Cut(ParseError::InvalidCharLiteral {
                    offset: char_offset,
                })
            })?;
            let ch = char::from_u32(cp).ok_or(ErrMode::Cut(ParseError::InvalidCharLiteral {
                offset: char_offset,
            }))?;
            return Ok(Edn::Char(ch));
        }

        // Named characters
        const NAMED: &[(&str, char)] = &[
            ("newline", '\n'),
            ("return", '\r'),
            ("space", ' '),
            ("tab", '\t'),
        ];

        for &(name, ch) in NAMED {
            if let Some(after) = s.strip_prefix(name) {
                let terminates = after.is_empty()
                    || after
                        .as_bytes()
                        .first()
                        .is_none_or(|&b| !is_symbol_char(b as char) || is_ws_byte(b));
                if terminates {
                    let _ = literal(name).parse_next(input)?;
                    return Ok(Edn::Char(ch));
                }
            }
        }

        // Multi-byte but not a recognized named char — error.
        // (Single chars were handled by the fast path above.)
        let c = s.chars().next().unwrap();
        if c.is_whitespace() {
            return Err(ErrMode::Cut(ParseError::InvalidCharLiteral {
                offset: char_offset,
            }));
        }
        let _ = any.parse_next(input)?;
        if let Some(&next) = rest(input).as_bytes().first() {
            if is_symbol_char(next as char) && !is_ws_byte(next) {
                return Err(ErrMode::Cut(ParseError::InvalidCharLiteral {
                    offset: char_offset,
                }));
            }
        }
        Ok(Edn::Char(c))
    }
}

fn edn_list<'a, 'b>(ctx: &'b ParseCtx<'_>) -> impl Parser<Input<'a>, Edn<'a>, E> + 'b
where
    'a: 'b,
{
    move |input: &mut Input<'a>| -> PResult<Edn<'a>> {
        let _ = '('.parse_next(input)?;
        ctx.enter_scope(input.current_token_start())?;
        let mut items = EdnSeq::new();
        loop {
            ws_and_comments(input, ctx);
            match peek_byte(input) {
                None => {
                    return Err(ErrMode::Cut(ParseError::UnexpectedEof {
                        offset: input.current_token_start(),
                    }))
                }
                Some(b')') => {
                    let _ = ')'.parse_next(input)?;
                    ctx.leave_scope();
                    return Ok(Edn::List(items));
                }
                _ => items.push(edn_value_dispatch(input, ctx)?),
            }
        }
    }
}

fn edn_vector<'a, 'b>(ctx: &'b ParseCtx<'_>) -> impl Parser<Input<'a>, Edn<'a>, E> + 'b
where
    'a: 'b,
{
    move |input: &mut Input<'a>| -> PResult<Edn<'a>> {
        let _ = '['.parse_next(input)?;
        ctx.enter_scope(input.current_token_start())?;
        let mut items = EdnSeq::new();
        loop {
            ws_and_comments(input, ctx);
            match peek_byte(input) {
                None => {
                    return Err(ErrMode::Cut(ParseError::UnexpectedEof {
                        offset: input.current_token_start(),
                    }))
                }
                Some(b']') => {
                    let _ = ']'.parse_next(input)?;
                    ctx.leave_scope();
                    return Ok(Edn::Vector(items));
                }
                _ => items.push(edn_value_dispatch(input, ctx)?),
            }
        }
    }
}

fn edn_map<'a, 'b>(ctx: &'b ParseCtx<'_>) -> impl Parser<Input<'a>, Edn<'a>, E> + 'b
where
    'a: 'b,
{
    move |input: &mut Input<'a>| -> PResult<Edn<'a>> {
        let _ = '{'.parse_next(input)?;
        ctx.enter_scope(input.current_token_start())?;
        let mut map = EdnMap::new();
        loop {
            ws_and_comments(input, ctx);
            match peek_byte(input) {
                None => {
                    return Err(ErrMode::Cut(ParseError::UnexpectedEof {
                        offset: input.current_token_start(),
                    }))
                }
                Some(b'}') => {
                    let _ = '}'.parse_next(input)?;
                    ctx.leave_scope();
                    return Ok(Edn::Map(map));
                }
                _ => {
                    let key_offset = input.current_token_start();
                    let key = edn_value_dispatch(input, ctx)?;
                    ws_and_comments(input, ctx);
                    let val = edn_value_dispatch(input, ctx)?;
                    if ctx.config.duplicate_keys == DuplicateKeyPolicy::Error
                        && map.contains_key(&key)
                    {
                        return Err(ErrMode::Cut(ParseError::DuplicateKey {
                            offset: key_offset,
                        }));
                    }
                    map.insert(key, val);
                }
            }
        }
    }
}

fn edn_set<'a, 'b>(ctx: &'b ParseCtx<'_>) -> impl Parser<Input<'a>, Edn<'a>, E> + 'b
where
    'a: 'b,
{
    move |input: &mut Input<'a>| -> PResult<Edn<'a>> {
        let _ = '{'.parse_next(input)?;
        ctx.enter_scope(input.current_token_start())?;
        let mut set = EdnSet::new();
        loop {
            ws_and_comments(input, ctx);
            match peek_byte(input) {
                None => {
                    return Err(ErrMode::Cut(ParseError::UnexpectedEof {
                        offset: input.current_token_start(),
                    }))
                }
                Some(b'}') => {
                    let _ = '}'.parse_next(input)?;
                    ctx.leave_scope();
                    return Ok(Edn::Set(set));
                }
                _ => {
                    set.insert(edn_value_dispatch(input, ctx)?);
                }
            }
        }
    }
}

fn edn_dispatch<'a, 'b>(ctx: &'b ParseCtx<'_>) -> impl Parser<Input<'a>, Edn<'a>, E> + 'b
where
    'a: 'b,
{
    move |input: &mut Input<'a>| -> PResult<Edn<'a>> {
        let _ = '#'.parse_next(input)?;
        let b = peek_byte(input).ok_or(ErrMode::Cut(ParseError::UnexpectedEof {
            offset: input.current_token_start(),
        }))?;
        match b {
            b'{' => edn_set(ctx).parse_next(input),
            b'_' => {
                let _ = '_'.parse_next(input)?;
                ws_and_comments(input, ctx);
                let _ = edn_value(ctx).parse_next(input)?;
                edn_value(ctx).parse_next(input)
            }
            _ => {
                let offset = input.current_token_start();
                let tag = raw_symbol(input)?;

                // Built-in tags pass through as Tagged(...) even without features.
                const BUILTIN_TAGS: &[&str] = &["inst", "uuid"];

                // Reject bare (unqualified) tags unless registered or built-in.
                // The serde path opts out via `allow_unknown_tags`.
                if !ctx.config.allow_unknown_tags
                    && !tag.contains('/')
                    && ctx.config.tag_readers.get(tag).is_none()
                    && !BUILTIN_TAGS.contains(&tag)
                {
                    return Err(ErrMode::Cut(ParseError::InvalidTag {
                        offset,
                        tag: tag.to_string(),
                    }));
                }

                ws_and_comments(input, ctx);
                let val = edn_value(ctx).parse_next(input)?;

                match ctx.config.tag_readers.get(tag) {
                    Some(reader) => reader(val).map_err(ErrMode::Cut),
                    None => Ok(Edn::Tagged(Cow::Borrowed(tag), Box::new(val))),
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use pretty_assertions::assert_eq;
    use rstest::rstest;

    use crate::collections::{EdnMap, EdnSet};
    use crate::config::{DuplicateKeyPolicy, ParseConfig};
    use crate::edn::{Edn, EdnKeyword};
    use crate::error::{EdnError, ParseError};
    use crate::reader::{read_all, read_string, read_string_with, Reader};

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
    #[case("##Inf")]
    #[case("##-Inf")]
    #[case("8E1313")]
    #[case("1e999")]
    #[case("-1e999")]
    fn test_read_string_error_special_float(#[case] input: &str) {
        assert!(read_string(input).is_err());
    }

    #[rstest]
    #[case(r#""""#, "")]
    #[case(r#""hello""#, "hello")]
    #[case(r#""hello world""#, "hello world")]
    #[case(r#""line\nbreak""#, "line\nbreak")]
    #[case(r#""tab\there""#, "tab\there")]
    #[case(r#""quote\"here""#, "quote\"here")]
    #[case(r#""back\\slash""#, "back\\slash")]
    #[case(r#""cr\rhere""#, "cr\rhere")]
    #[case("\"é\"", "é")]
    #[case("\"α\"", "α")]
    #[case("\"世界\"", "世界")]
    #[case("\"hello é world\"", "hello é world")]
    fn test_read_string_str(#[case] input: &str, #[case] expected: &str) {
        assert_eq!(
            read_string(input).unwrap(),
            Edn::Str(Cow::Owned(expected.to_string()))
        );
    }

    #[test]
    fn test_read_string_unicode_escape() {
        assert_eq!(
            read_string(r#""\u0041""#).unwrap(),
            Edn::Str(Cow::Borrowed("A"))
        );
        assert_eq!(
            read_string(r#""\u00e9""#).unwrap(),
            Edn::Str(Cow::Owned("é".to_string()))
        );
    }

    #[rstest]
    #[case(r#""\b""#)]
    #[case(r#""\f""#)]
    fn test_read_string_error_invalid_string_escape(#[case] input: &str) {
        assert!(read_string(input).is_err());
    }

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

    #[rstest]
    #[case(":foo", "foo")]
    #[case(":ns/name", "ns/name")]
    #[case(":a.b/c", "a.b/c")]
    fn test_read_string_keyword(#[case] input: &str, #[case] expected_name: &str) {
        let val = read_string(input).unwrap();
        match &val {
            Edn::Keyword(k) => assert_eq!(k.as_str(), expected_name),
            other => panic!("expected EdnKeyword, got {other:?}"),
        }
    }

    #[test]
    fn test_keyword_namespace() {
        let k = EdnKeyword::new("ns/name");
        assert_eq!(k.namespace(), Some("ns"));
        assert_eq!(k.name(), "name");

        let k2 = EdnKeyword::new("bare");
        assert_eq!(k2.namespace(), None);
        assert_eq!(k2.name(), "bare");
    }

    #[rstest]
    #[case("foo", "foo")]
    #[case("ns/name", "ns/name")]
    #[case("my.ns/sym", "my.ns/sym")]
    fn test_read_string_symbol(#[case] input: &str, #[case] expected_name: &str) {
        let val = read_string(input).unwrap();
        match &val {
            Edn::Symbol(s) => assert_eq!(s.as_str(), expected_name),
            other => panic!("expected EdnSymbol, got {other:?}"),
        }
    }

    #[rstest]
    #[case("()", Edn::List(vec![].into()))]
    #[case("(1 2 3)", Edn::List(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)].into()))]
    #[case("(nil true)", Edn::List(vec![Edn::Nil, Edn::Bool(true)].into()))]
    fn test_read_string_list(#[case] input: &str, #[case] expected: Edn<'_>) {
        assert_eq!(read_string(input).unwrap(), expected);
    }

    #[rstest]
    #[case("[]", Edn::Vector(vec![].into()))]
    #[case("[1 2 3]", Edn::Vector(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)].into()))]
    fn test_read_string_vector(#[case] input: &str, #[case] expected: Edn<'_>) {
        assert_eq!(read_string(input).unwrap(), expected);
    }

    #[test]
    fn test_read_string_map() {
        let val = read_string("{:a 1 :b 2}").unwrap();
        let mut expected = EdnMap::new();
        expected.insert(Edn::Keyword(EdnKeyword::new("a")), Edn::Int(1));
        expected.insert(Edn::Keyword(EdnKeyword::new("b")), Edn::Int(2));
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

    #[test]
    fn test_read_string_nested() {
        let val = read_string("{:items [1 (2 3)] :flag true}").unwrap();
        let mut expected = EdnMap::new();
        expected.insert(
            Edn::Keyword(EdnKeyword::new("items")),
            Edn::Vector(
                vec![
                    Edn::Int(1),
                    Edn::List(vec![Edn::Int(2), Edn::Int(3)].into()),
                ]
                .into(),
            ),
        );
        expected.insert(Edn::Keyword(EdnKeyword::new("flag")), Edn::Bool(true));
        assert_eq!(val, Edn::Map(expected));
    }

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

    #[rstest]
    #[case("; comment\n12", Edn::Int(12))]
    #[case("12 ; trailing", Edn::Int(12))]
    fn test_read_string_comment(#[case] input: &str, #[case] expected: Edn<'_>) {
        assert_eq!(read_string(input).unwrap(), expected);
    }

    #[rstest]
    #[case("#_ foo 12", Edn::Int(12))]
    #[case("[1 #_ 2 3]", Edn::Vector(vec![Edn::Int(1), Edn::Int(3)].into()))]
    fn test_read_string_discard(#[case] input: &str, #[case] expected: Edn<'_>) {
        assert_eq!(read_string(input).unwrap(), expected);
    }

    #[rstest]
    #[case("  12  ", Edn::Int(12))]
    #[case("\t12\n", Edn::Int(12))]
    #[case(",12,", Edn::Int(12))]
    #[case("[1,,2,,3]", Edn::Vector(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)].into()))]
    fn test_read_string_whitespace(#[case] input: &str, #[case] expected: Edn<'_>) {
        assert_eq!(read_string(input).unwrap(), expected);
    }

    #[rstest]
    #[case("", EdnError::Parse(ParseError::UnexpectedEof { offset: 0 }))]
    #[case("   ", EdnError::Parse(ParseError::UnexpectedEof { offset: 3 }))]
    fn test_read_string_error_empty(#[case] input: &str, #[case] expected: EdnError) {
        assert_eq!(read_string(input).unwrap_err(), expected);
    }

    #[test]
    fn test_read_string_error_trailing() {
        assert_eq!(
            read_string("12 43").unwrap_err(),
            EdnError::Parse(ParseError::TrailingContent { offset: 3 }),
        );
    }

    #[test]
    fn test_read_string_error_duplicate_key() {
        assert_eq!(
            read_string("{:a 1 :a 2}").unwrap_err(),
            EdnError::Parse(ParseError::DuplicateKey { offset: 6 }),
        );
    }

    #[rstest]
    #[case(r#""\q""#, EdnError::Parse(ParseError::InvalidEscape { offset: 1 }))]
    #[case(r#""\u000G""#, EdnError::Parse(ParseError::InvalidEscape { offset: 1 }))]
    fn test_read_string_error_invalid_escape(#[case] input: &str, #[case] expected: EdnError) {
        assert_eq!(read_string(input).unwrap_err(), expected);
    }

    #[rstest]
    #[case(r"\abc", EdnError::Parse(ParseError::InvalidCharLiteral { offset: 1 }))]
    fn test_read_string_error_invalid_char_literal(
        #[case] input: &str,
        #[case] expected: EdnError,
    ) {
        assert_eq!(read_string(input).unwrap_err(), expected);
    }

    #[test]
    fn test_read_string_error_unexpected_token() {
        let err = read_string("##xyz").unwrap_err();
        assert!(
            matches!(
                err,
                EdnError::Parse(ParseError::UnexpectedToken { .. })
                    | EdnError::Parse(ParseError::InvalidSymbol { .. })
            ),
            "expected error for ##xyz, got {err:?}"
        );
    }

    #[cfg(not(feature = "bignum"))]
    #[rstest]
    #[case("12N", EdnError::Parse(ParseError::UnsupportedFeature { offset: 0, feature: "bignum" }))]
    #[case("12M", EdnError::Parse(ParseError::UnsupportedFeature { offset: 0, feature: "bignum" }))]
    fn test_read_string_error_unsupported_feature(#[case] input: &str, #[case] expected: EdnError) {
        assert_eq!(read_string(input).unwrap_err(), expected);
    }

    #[test]
    fn test_read_string_error_unclosed_string() {
        let err = read_string(r#""hello"#).unwrap_err();
        assert!(
            matches!(err, EdnError::Parse(ParseError::UnexpectedEof { .. })),
            "expected UnexpectedEof, got {err:?}"
        );
    }

    #[test]
    fn test_read_string_error_unclosed_vector() {
        let err = read_string("[1 2").unwrap_err();
        assert!(
            matches!(err, EdnError::Parse(ParseError::UnexpectedEof { .. })),
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
        expected.insert(Edn::Keyword(EdnKeyword::new("a")), Edn::Int(2));
        assert_eq!(val, Edn::Map(expected));
    }

    #[cfg(not(feature = "bignum"))]
    #[test]
    fn test_read_string_error_integer_overflow() {
        let err = read_string("99999999999999999999").unwrap_err();
        assert!(
            matches!(
                err,
                EdnError::Parse(ParseError::InvalidNumber { offset: 0 })
            ),
            "expected InvalidNumber, got {err:?}"
        );
    }

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
    fn test_edn_get() {
        let val = read_string("{:name \"Alice\" :age 30}").unwrap();
        assert_eq!(val.get_keyword("name").unwrap().as_str(), Some("Alice"));
        assert_eq!(val.get_keyword("age").unwrap().as_i64(), Some(30));
        assert!(val.get_keyword("missing").is_none());
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

    #[rstest]
    #[case(":0")]
    #[case(":1")]
    #[case(":123")]
    fn test_read_string_error_digit_keywords(#[case] input: &str) {
        assert!(read_string(input).is_err());
    }

    #[test]
    fn test_read_string_error_formfeed_backspace() {
        assert!(read_string("\\formfeed").is_err());
        assert!(read_string("\\backspace").is_err());
    }

    #[test]
    fn test_read_string_depth_limit() {
        let depth = super::MAX_DEPTH as usize;
        let open: String = "[".repeat(depth + 1);
        let close: String = "]".repeat(depth + 1);
        let input = format!("{open}1{close}");
        let err = read_string(&input).unwrap_err();
        assert!(err.to_string().contains("recursion limit"), "{err}");
    }

    #[test]
    fn test_read_string_at_depth_limit() {
        let depth = super::MAX_DEPTH as usize;
        let open: String = "[".repeat(depth);
        let close: String = "]".repeat(depth);
        let input = format!("{open}1{close}");
        assert!(read_string(&input).is_ok());
    }

    #[test]
    fn test_read_string_discard_depth_limit() {
        // Deeply nested collections inside a discard are still caught.
        let depth = super::MAX_DEPTH as usize;
        let open: String = "[".repeat(depth + 1);
        let close: String = "]".repeat(depth + 1);
        let input = format!("[#_ {open}1{close} 2]");
        assert!(read_string(&input).is_err());
    }
}
