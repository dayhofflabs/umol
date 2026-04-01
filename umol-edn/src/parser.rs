//! Winnow-based EDN parser.

use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet};

use winnow::ascii::digit1;
use winnow::combinator::opt;
use winnow::error::{ContextError, ErrMode};
use winnow::stream::{Location, Stream};
use winnow::token::{any, one_of, take_while};
use winnow::{LocatingSlice, ModalResult, Parser};

use crate::edn::{Edn, Keyword, Symbol};
use crate::error::EdnError;

/// Parser configuration.
#[derive(Clone, Debug)]
pub struct ParseConfig {
    pub strict: bool,
    pub duplicate_keys: DuplicateKeyPolicy,
}

impl Default for ParseConfig {
    fn default() -> Self {
        ParseConfig {
            strict: false,
            duplicate_keys: DuplicateKeyPolicy::Error,
        }
    }
}

/// Behavior on duplicate map keys.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DuplicateKeyPolicy {
    Error,
    LastWins,
}

type Input<'a> = LocatingSlice<&'a str>;
type E = ErrMode<ContextError>;

/// Get the remaining input as a `&str`.
fn rest<'a>(input: &Input<'a>) -> &'a str {
    *input.as_ref()
}

/// Parse a single EDN value from the input, returning the value and remaining input.
pub fn parse_value<'a>(
    input: &'a str,
    config: &ParseConfig,
) -> Result<(Edn<'a>, &'a str), EdnError> {
    let mut located = LocatingSlice::new(input);
    ws_and_comments(&mut located).map_err(|_| EdnError::UnexpectedEof { offset: 0 })?;
    let val = edn_value(config)
        .parse_next(&mut located)
        .map_err(|e| to_edn_error(e))?;
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
    ws_and_comments(&mut loc).ok();
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
        ws_and_comments(&mut located).ok();
        if rest(&located).is_empty() {
            break;
        }
        let val = edn_value(config)
            .parse_next(&mut located)
            .map_err(|e| to_edn_error(e))?;
        values.push(val);
    }
    Ok(values)
}

fn to_edn_error(e: E) -> EdnError {
    match e {
        ErrMode::Backtrack(ctx) | ErrMode::Cut(ctx) => {
            let msg = format!("{ctx}");
            if msg.is_empty() {
                EdnError::Custom("parse error".to_string())
            } else {
                EdnError::Custom(msg)
            }
        }
        ErrMode::Incomplete(_) => EdnError::UnexpectedEof { offset: 0 },
    }
}

// --- Whitespace and comments ---

fn is_ws(c: char) -> bool {
    matches!(c, ' ' | '\t' | '\n' | '\r' | ',')
}

fn ws_and_comments<'a>(input: &mut Input<'a>) -> ModalResult<()> {
    loop {
        let before = input.current_token_start();
        take_while(0.., is_ws).parse_next(input)?;
        let s = rest(input);
        if s.starts_with(';') {
            take_while(0.., |c: char| c != '\n').parse_next(input)?;
            opt('\n').parse_next(input)?;
        } else if s.starts_with("#_") {
            let _ = "#_".parse_next(input)?;
            ws_and_comments(input).ok();
            let discard_config = ParseConfig {
                strict: false,
                duplicate_keys: DuplicateKeyPolicy::LastWins,
            };
            let _ = edn_value(&discard_config).parse_next(input)?;
        } else if input.current_token_start() == before {
            break;
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
    move |input: &mut Input<'a>| -> ModalResult<Edn<'a>> {
        ws_and_comments(input).ok();
        let s = rest(input);
        let c = s
            .chars()
            .next()
            .ok_or(ErrMode::Backtrack(ContextError::new()))?;
        match c {
            '(' => edn_list(config).parse_next(input),
            '[' => edn_vector(config).parse_next(input),
            '{' => edn_map(config).parse_next(input),
            '#' => edn_dispatch(config).parse_next(input),
            ':' => edn_keyword(input),
            '"' => edn_string(config.strict).parse_next(input),
            '\\' => edn_char(config.strict).parse_next(input),
            '+' | '-' => edn_number_or_symbol(config).parse_next(input),
            c if c.is_ascii_digit() => edn_number(input),
            _ => edn_symbol_or_literal(input),
        }
    }
}

// --- Nil, booleans, symbols ---

fn is_symbol_start(c: char) -> bool {
    matches!(
        c,
        'a'..='z' | 'A'..='Z' | '.' | '*' | '!' | '_' | '?' | '$' | '%' | '&' | '=' | '<' | '>' | '/'
    )
}

fn is_symbol_char(c: char) -> bool {
    is_symbol_start(c) || matches!(c, '0'..='9' | '+' | '-' | '#' | ':' | '\'')
}

fn raw_symbol<'a>(input: &mut Input<'a>) -> ModalResult<&'a str> {
    (one_of(is_symbol_start), take_while(0.., is_symbol_char))
        .take()
        .parse_next(input)
}

fn edn_symbol_or_literal<'a>(input: &mut Input<'a>) -> ModalResult<Edn<'a>> {
    let s = raw_symbol(input)?;
    match s {
        "nil" => Ok(Edn::Nil),
        "true" => Ok(Edn::Bool(true)),
        "false" => Ok(Edn::Bool(false)),
        _ => Ok(Edn::Symbol(Symbol::new(s))),
    }
}

// --- Keywords ---

fn edn_keyword<'a>(input: &mut Input<'a>) -> ModalResult<Edn<'a>> {
    let _ = ':'.parse_next(input)?;
    let s = raw_symbol(input)?;
    Ok(Edn::Keyword(Keyword::new(s)))
}

// --- Numbers ---

fn edn_number<'a>(input: &mut Input<'a>) -> ModalResult<Edn<'a>> {
    let num_str: &str = (
        opt(one_of(['+', '-'])),
        digit1,
        opt(('.', digit1)),
        opt((one_of(['e', 'E']), opt(one_of(['+', '-'])), digit1)),
    )
        .take()
        .parse_next(input)?;

    // Check for N/M suffix (bignum)
    let s = rest(input);
    if s.starts_with('N') || s.starts_with('M') {
        let _ = any.parse_next(input)?;
        return Err(ErrMode::Cut(ContextError::new()));
    }

    if num_str.contains('.') || num_str.contains('e') || num_str.contains('E') {
        let f: f64 = num_str
            .parse()
            .map_err(|_| ErrMode::Cut(ContextError::new()))?;
        Ok(Edn::Float(f))
    } else {
        let n: i64 = num_str
            .parse()
            .map_err(|_| ErrMode::Cut(ContextError::new()))?;
        Ok(Edn::Int(n))
    }
}

fn edn_number_or_symbol<'a, 'b>(
    _config: &'b ParseConfig,
) -> impl Parser<Input<'a>, Edn<'a>, E> + 'b
where
    'a: 'b,
{
    move |input: &mut Input<'a>| -> ModalResult<Edn<'a>> {
        let checkpoint = input.checkpoint();
        if let Ok(val) = edn_number(input) {
            return Ok(val);
        }
        input.reset(&checkpoint);
        edn_symbol_or_literal(input)
    }
}

// --- Strings ---

fn edn_string<'a>(strict: bool) -> impl Parser<Input<'a>, Edn<'a>, E> {
    move |input: &mut Input<'a>| -> ModalResult<Edn<'a>> {
        let _ = '"'.parse_next(input)?;
        let mut result = String::new();

        loop {
            let s = rest(input);
            let c = s
                .chars()
                .next()
                .ok_or(ErrMode::Cut(ContextError::new()))?;
            match c {
                '"' => {
                    let _ = any.parse_next(input)?;
                    return Ok(Edn::Str(Cow::Owned(result)));
                }
                '\\' => {
                    let _ = any.parse_next(input)?;
                    let s2 = rest(input);
                    let esc = s2
                        .chars()
                        .next()
                        .ok_or(ErrMode::Cut(ContextError::new()))?;
                    let _ = any.parse_next(input)?;
                    match esc {
                        't' => result.push('\t'),
                        'r' => result.push('\r'),
                        'n' => result.push('\n'),
                        '\\' => result.push('\\'),
                        '"' => result.push('"'),
                        'b' if !strict => result.push('\u{0008}'),
                        'f' if !strict => result.push('\u{000C}'),
                        'u' => {
                            let hex: &str =
                                take_while(4..=4, |c: char| c.is_ascii_hexdigit())
                                    .parse_next(input)?;
                            let cp = u32::from_str_radix(hex, 16)
                                .map_err(|_| ErrMode::Cut(ContextError::new()))?;
                            let ch = char::from_u32(cp)
                                .ok_or_else(|| ErrMode::Cut(ContextError::new()))?;
                            result.push(ch);
                        }
                        '0'..='7' if !strict => {
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
                                return Err(ErrMode::Cut(ContextError::new()));
                            }
                            let ch = char::from_u32(val)
                                .ok_or_else(|| ErrMode::Cut(ContextError::new()))?;
                            result.push(ch);
                        }
                        _ => return Err(ErrMode::Cut(ContextError::new())),
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

fn edn_char<'a>(strict: bool) -> impl Parser<Input<'a>, Edn<'a>, E> {
    move |input: &mut Input<'a>| -> ModalResult<Edn<'a>> {
        let _ = '\\'.parse_next(input)?;
        let s = rest(input);

        let named: &[(&str, char)] = if strict {
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
                    || after
                        .chars()
                        .next()
                        .map_or(true, |c| !is_symbol_char(c) || is_ws(c));
                if terminates {
                    let _ = winnow::token::literal(name).parse_next(input)?;
                    return Ok(Edn::Char(ch));
                }
            }
        }

        // Unicode escape \uNNNN
        if s.starts_with('u') {
            let _ = 'u'.parse_next(input)?;
            let hex: &str = take_while(4..=4, |c: char| c.is_ascii_hexdigit())
                .parse_next(input)?;
            let cp = u32::from_str_radix(hex, 16)
                .map_err(|_| ErrMode::Cut(ContextError::new()))?;
            let ch =
                char::from_u32(cp).ok_or_else(|| ErrMode::Cut(ContextError::new()))?;
            return Ok(Edn::Char(ch));
        }

        // Single character
        let c = s
            .chars()
            .next()
            .ok_or(ErrMode::Cut(ContextError::new()))?;
        if c.is_whitespace() {
            return Err(ErrMode::Cut(ContextError::new()));
        }
        let _ = any.parse_next(input)?;
        if let Some(next) = rest(input).chars().next() {
            if is_symbol_char(next) && !is_ws(next) {
                return Err(ErrMode::Cut(ContextError::new()));
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
    move |input: &mut Input<'a>| -> ModalResult<Edn<'a>> {
        let _ = '('.parse_next(input)?;
        let mut items = Vec::new();
        loop {
            ws_and_comments(input).ok();
            let s = rest(input);
            if s.is_empty() {
                return Err(ErrMode::Cut(ContextError::new()));
            }
            if s.starts_with(')') {
                let _ = ')'.parse_next(input)?;
                return Ok(Edn::List(items));
            }
            items.push(edn_value(config).parse_next(input)?);
        }
    }
}

fn edn_vector<'a, 'b>(
    config: &'b ParseConfig,
) -> impl Parser<Input<'a>, Edn<'a>, E> + 'b
where
    'a: 'b,
{
    move |input: &mut Input<'a>| -> ModalResult<Edn<'a>> {
        let _ = '['.parse_next(input)?;
        let mut items = Vec::new();
        loop {
            ws_and_comments(input).ok();
            let s = rest(input);
            if s.is_empty() {
                return Err(ErrMode::Cut(ContextError::new()));
            }
            if s.starts_with(']') {
                let _ = ']'.parse_next(input)?;
                return Ok(Edn::Vector(items));
            }
            items.push(edn_value(config).parse_next(input)?);
        }
    }
}

fn edn_map<'a, 'b>(
    config: &'b ParseConfig,
) -> impl Parser<Input<'a>, Edn<'a>, E> + 'b
where
    'a: 'b,
{
    move |input: &mut Input<'a>| -> ModalResult<Edn<'a>> {
        let _ = '{'.parse_next(input)?;
        let mut map = BTreeMap::new();
        loop {
            ws_and_comments(input).ok();
            let s = rest(input);
            if s.is_empty() {
                return Err(ErrMode::Cut(ContextError::new()));
            }
            if s.starts_with('}') {
                let _ = '}'.parse_next(input)?;
                return Ok(Edn::Map(map));
            }
            let key = edn_value(config).parse_next(input)?;
            ws_and_comments(input).ok();
            let val = edn_value(config).parse_next(input)?;
            if config.duplicate_keys == DuplicateKeyPolicy::Error && map.contains_key(&key) {
                return Err(ErrMode::Cut(ContextError::new()));
            }
            map.insert(key, val);
        }
    }
}

fn edn_set<'a, 'b>(
    config: &'b ParseConfig,
) -> impl Parser<Input<'a>, Edn<'a>, E> + 'b
where
    'a: 'b,
{
    move |input: &mut Input<'a>| -> ModalResult<Edn<'a>> {
        let _ = '{'.parse_next(input)?;
        let mut set = BTreeSet::new();
        loop {
            ws_and_comments(input).ok();
            let s = rest(input);
            if s.is_empty() {
                return Err(ErrMode::Cut(ContextError::new()));
            }
            if s.starts_with('}') {
                let _ = '}'.parse_next(input)?;
                return Ok(Edn::Set(set));
            }
            set.insert(edn_value(config).parse_next(input)?);
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
    move |input: &mut Input<'a>| -> ModalResult<Edn<'a>> {
        let _ = '#'.parse_next(input)?;
        let s = rest(input);
        let c = s
            .chars()
            .next()
            .ok_or(ErrMode::Cut(ContextError::new()))?;
        match c {
            '{' => edn_set(config).parse_next(input),
            '_' => {
                let _ = '_'.parse_next(input)?;
                ws_and_comments(input).ok();
                let _ = edn_value(config).parse_next(input)?;
                edn_value(config).parse_next(input)
            }
            '#' => {
                let _ = '#'.parse_next(input)?;
                let s2 = rest(input);
                if s2.starts_with("NaN") {
                    let _ = "NaN".parse_next(input)?;
                    Ok(Edn::Float(f64::NAN))
                } else if s2.starts_with("-Inf") {
                    let _ = "-Inf".parse_next(input)?;
                    Ok(Edn::Float(f64::NEG_INFINITY))
                } else if s2.starts_with("Inf") {
                    let _ = "Inf".parse_next(input)?;
                    Ok(Edn::Float(f64::INFINITY))
                } else {
                    Err(ErrMode::Cut(ContextError::new()))
                }
            }
            _ => {
                let tag = raw_symbol(input)?;
                ws_and_comments(input).ok();
                let val = edn_value(config).parse_next(input)?;
                Ok(Edn::Tagged(tag.to_string(), Box::new(val)))
            }
        }
    }
}
