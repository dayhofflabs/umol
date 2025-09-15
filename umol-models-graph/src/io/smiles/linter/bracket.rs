//! Bracket microparser and helpers used by Bracket rules.

use indexmap::IndexMap;
use nom::branch::alt;
use nom::bytes::complete::tag;
use nom::character::complete::{digit1, one_of};
use nom::combinator::{map, map_opt, opt, recognize, value};
use nom::error::Error;
use nom::multi::fold_many0;
use nom::sequence::{pair, preceded};
use nom::Parser;
use smallvec::SmallVec;
use umol_data::Element;

use crate::diagnostics::{Category, Code, Diagnostic, DiagnosticsReport, Span};
use crate::io::ir::Chirality;
use crate::io::smiles::parser::utils::BracketField;

const DIGITS: &str = "0123456789";

#[derive(Default, Debug, Clone, Copy)]
pub struct BracketParsed {
    pub element: Option<Element>,
    pub isotope: Option<u32>,
    pub hcount: Option<u32>,
    pub charge: Option<i32>,
    pub class: Option<u32>,
}

pub fn find_closing_bracket(bytes: &[u8], mut idx: usize) -> Option<usize> {
    while idx < bytes.len() {
        if bytes[idx] == b']' {
            return Some(idx);
        }
        idx += 1;
    }
    None
}

pub fn is_bare_organic(inner: &str) -> bool {
    matches!(
        inner,
        "B" | "C" | "N" | "O" | "S" | "P" | "F" | "Cl" | "Br" | "I"
    )
}

pub fn inner_contains_h1(inner: &str) -> bool {
    inner.contains("H1")
}

pub fn find_subslice(hay: &str, needle: &str) -> Option<(usize, usize)> {
    hay.find(needle).map(|s| (s, s + needle.len()))
}

pub fn find_charge_plus_minus_one(inner: &str) -> Option<(usize, usize)> {
    let bytes = inner.as_bytes();
    let mut i = 0usize;
    while i + 1 < bytes.len() {
        if (bytes[i] == b'+' || bytes[i] == b'-') && bytes[i + 1] == b'1' {
            let end = i + 2;
            if end >= bytes.len() || !bytes[end].is_ascii_digit() {
                return Some((i, end));
            }
        }
        i += 1;
    }
    None
}

pub fn bracket_order_misordered(inner: &str) -> bool {
    let chiral_idx = find_first_chiral(inner);
    let h_idx = inner.find('H').map(|x| x as isize).unwrap_or(-1);
    let charge_idx = find_first_charge(inner);
    let class_idx = inner.find(':').map(|x| x as isize).unwrap_or(-1);
    let present = [chiral_idx, h_idx, charge_idx, class_idx]
        .into_iter()
        .filter(|&i| i >= 0)
        .count();
    if present < 2 {
        return false;
    }
    let mut last = -1isize;
    for idx in [chiral_idx, h_idx, charge_idx, class_idx] {
        if idx >= 0 {
            if idx < last {
                return true;
            }
            last = idx;
        }
    }
    false
}

pub fn find_first_chiral(inner: &str) -> isize {
    let patterns = ["@@", "@TH", "@AL", "@SP", "@TB", "@OH", "@"]; // '@@' before '@'
    let mut best: Option<usize> = None;
    for p in patterns.iter() {
        if let Some(i) = inner.find(p) {
            best = Some(best.map_or(i, |b| b.min(i)));
        }
    }
    best.map(|x| x as isize).unwrap_or(-1)
}

pub fn find_first_charge(inner: &str) -> isize {
    let mut best: Option<usize> = None;
    for p in ["++", "--", "+", "-"] {
        if let Some(i) = inner.find(p) {
            best = Some(best.map_or(i, |b| b.min(i)));
        }
    }
    best.map(|x| x as isize).unwrap_or(-1)
}

pub fn find_h_two_digits(inner: &str) -> Option<(usize, usize)> {
    let bytes = inner.as_bytes();
    let mut i = 0usize;
    while i + 2 < bytes.len() {
        if bytes[i] == b'H' && bytes[i + 1].is_ascii_digit() && bytes[i + 2].is_ascii_digit() {
            return Some((i, i + 3));
        }
        i += 1;
    }
    None
}

pub fn find_class_issues(inner: &str) -> Option<(usize, usize, bool)> {
    let bytes = inner.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() {
        if bytes[i] == b':' {
            if i + 1 >= bytes.len() {
                return Some((i, i + 1, false));
            }
            if bytes[i + 1] == b'-' {
                return Some((i, i + 2, true));
            }
            if !bytes[i + 1].is_ascii_digit() {
                return Some((i, i + 2, false));
            }
        }
        i += 1;
    }
    None
}

pub fn lint_style_percent_single_digit(input: &str, report: &mut DiagnosticsReport) {
    let bytes = input.as_bytes();
    let mut i = 0usize;
    while i + 2 < bytes.len() {
        if bytes[i] == b'%' && bytes[i + 1] == b'0' && (b'1'..=b'9').contains(&bytes[i + 2]) {
            report.push(Diagnostic::warning(
                Code("STYLE_UNNECESSARY_PERCENT_RING_INDEX"),
                Category::Style,
                Span::new(i, i + 3),
                "Prefer single-digit ring number for 1..9",
            ));
            i += 3;
            continue;
        }
        i += 1;
    }
}

pub fn lint_trailing_bond(input: &str, report: &mut DiagnosticsReport) {
    let bytes = input.as_bytes();
    let mut end = bytes.len();
    while end > 0 && matches!(bytes[end - 1], b' ' | b'\t' | b'\n' | b'\r') {
        end -= 1;
    }
    if end == 0 {
        return;
    }
    match bytes[end - 1] {
        b'-' | b'=' | b'#' | b'$' | b':' | b'/' | b'\\' => {
            report.push(Diagnostic::error(
                Code("SYN_TRAILING_BOND"),
                Category::Syn,
                Span::new(end - 1, end),
                "Trailing bond symbol",
            ));
        }
        _ => {}
    }
}

pub fn lint_dot_before_ring(input: &str, report: &mut DiagnosticsReport) {
    let bytes = input.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() {
        if bytes[i] == b'.' {
            if i + 1 < bytes.len() {
                if bytes[i + 1].is_ascii_digit() {
                    report.push(Diagnostic::error(
                        Code("SYN_DOT_BEFORE_RING"),
                        Category::Syn,
                        Span::new(i, i + 2),
                        "Dot before ring index is invalid",
                    ));
                } else if bytes[i + 1] == b'%' {
                    let mut end = i + 2;
                    if i + 3 < bytes.len()
                        && bytes[i + 2].is_ascii_digit()
                        && bytes[i + 3].is_ascii_digit()
                    {
                        end = i + 4;
                    }
                    report.push(Diagnostic::error(
                        Code("SYN_DOT_BEFORE_RING"),
                        Category::Syn,
                        Span::new(i, end),
                        "Dot before ring index is invalid",
                    ));
                }
            }
        }
        i += 1;
    }
}

pub fn lint_intertoken_whitespace(input: &str, report: &mut DiagnosticsReport) {
    let bytes = input.as_bytes();
    let mut last_non_ws: isize = -1;
    for (i, &b) in bytes.iter().enumerate() {
        if !matches!(b, b' ' | b'\t' | b'\n' | b'\r') {
            last_non_ws = i as isize;
        }
    }
    if last_non_ws < 0 {
        return;
    }
    let mut i = 0usize;
    while i < (last_non_ws as usize) {
        if matches!(bytes[i], b' ' | b'\t' | b'\n' | b'\r') {
            let start = i;
            while i < bytes.len() && matches!(bytes[i], b' ' | b'\t' | b'\n' | b'\r') {
                i += 1;
            }
            report.push(Diagnostic::error(
                Code("LEX_INTERTOKEN_WHITESPACE"),
                Category::Lex,
                Span::new(start, i),
                "Inter-token whitespace is not allowed",
            ));
            continue;
        }
        i += 1;
    }
}

pub fn lint_style_bonds(input: &str, report: &mut DiagnosticsReport) {
    let bytes = input.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() {
        match bytes[i] {
            b'[' => {
                let mut j = i + 1;
                while j < bytes.len() && bytes[j] != b']' {
                    j += 1;
                }
                i = j.saturating_add(1);
                continue;
            }
            b':' => {
                let prev = input[..i]
                    .as_bytes()
                    .iter()
                    .rfind(|&&b| b != b' ' && b != b'\t' && b != b'\n' && b != b'\r')
                    .copied();
                let next = input[i + 1..]
                    .as_bytes()
                    .iter()
                    .find(|&&b| b != b' ' && b != b'\t' && b != b'\n' && b != b'\r')
                    .copied();
                let is_arom = |b: u8| matches!(b, b'b' | b'c' | b'n' | b'o' | b'p' | b's');
                if prev.map(is_arom).unwrap_or(false) && next.map(is_arom).unwrap_or(false) {
                    report.push(Diagnostic::warning(
                        Code("STYLE_EXPLICIT_AROMATIC_BOND"),
                        Category::Style,
                        Span::new(i, i + 1),
                        "Avoid explicit ':' between aromatic atoms",
                    ));
                }
            }
            b'-' => {
                let prev = input[..i]
                    .as_bytes()
                    .iter()
                    .rfind(|&&b| !b.is_ascii_whitespace())
                    .copied();
                let next = input[i + 1..]
                    .as_bytes()
                    .iter()
                    .find(|&&b| !b.is_ascii_whitespace())
                    .copied();
                let is_arom = |b: u8| matches!(b, b'b' | b'c' | b'n' | b'o' | b'p' | b's');
                if !(prev.map(is_arom).unwrap_or(false) && next.map(is_arom).unwrap_or(false)) {
                    report.push(Diagnostic::warning(
                        Code("STYLE_EXPLICIT_SINGLE_BOND"),
                        Category::Style,
                        Span::new(i, i + 1),
                        "Avoid explicit '-' when default applies",
                    ));
                }
            }
            _ => {}
        }
        i += 1;
    }
}

pub fn lint_ring_style(input: &str, report: &mut DiagnosticsReport) {
    let mut i = 0usize;
    let bytes = input.as_bytes();
    let mut used: Vec<u32> = Vec::new();
    let mut counts: IndexMap<u32, u32> = IndexMap::new();
    while i < bytes.len() {
        match bytes[i] {
            b'[' => {
                let mut j = i + 1;
                while j < bytes.len() && bytes[j] != b']' {
                    j += 1;
                }
                i = j + 1;
                continue;
            }
            b'%' => {
                if i + 2 < bytes.len()
                    && bytes[i + 1].is_ascii_digit()
                    && bytes[i + 2].is_ascii_digit()
                {
                    let val = (bytes[i + 1] - b'0') as u32 * 10 + (bytes[i + 2] - b'0') as u32;
                    used.push(val);
                    *counts.entry(val).or_insert(0) += 1;
                    i += 3;
                    continue;
                }
            }
            b'0'..=b'9' => {
                let val = (bytes[i] - b'0') as u32;
                used.push(val);
                *counts.entry(val).or_insert(0) += 1;
                i += 1;
                continue;
            }
            _ => {}
        }
        i += 1;
    }
    if !used.is_empty() {
        let mut set = used.clone();
        set.sort_unstable();
        set.dedup();
        let first = set[0];
        if first != 1 {
            report.push(Diagnostic::warning(
                Code("STYLE_FIRST_RING_NOT_ONE"),
                Category::Style,
                Span::new(0, 0),
                "Prefer starting ring numbering at 1",
            ));
        }
        if set.len() >= 2 {
            let mut prev = set[0];
            for &v in &set[1..] {
                if v > prev + 1 {
                    report.push(Diagnostic::warning(
                        Code("STYLE_NONCONSECUTIVE_RING_NUMBERING"),
                        Category::Style,
                        Span::new(0, 0),
                        "Prefer consecutive ring numbering",
                    ));
                    break;
                }
                prev = v;
            }
        }
        for (_k, c) in counts.iter() {
            if *c > 2 {
                report.push(Diagnostic::warning(
                    Code("STYLE_REUSED_RING_INDICES"),
                    Category::Style,
                    Span::new(0, 0),
                    "Avoid reusing the same ring number",
                ));
                break;
            }
        }
    }
}

pub fn lint_dot_positions(input: &str, report: &mut DiagnosticsReport) {
    let bytes = input.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() && bytes[i].is_ascii_whitespace() {
        i += 1;
    }
    if i < bytes.len() && bytes[i] == b'.' {
        report.push(Diagnostic::error(
            Code("SYN_LEADING_DOT"),
            Category::Syn,
            Span::new(i, i + 1),
            "Leading dot",
        ));
    }
    let mut j = bytes.len();
    while j > 0 && bytes[j - 1].is_ascii_whitespace() {
        j -= 1;
    }
    if j > 0 && bytes[j - 1] == b'.' {
        report.push(Diagnostic::error(
            Code("SYN_TRAILING_DOT"),
            Category::Syn,
            Span::new(j - 1, j),
            "Trailing dot",
        ));
    }
    let mut k = 0usize;
    while k + 1 < bytes.len() {
        if bytes[k] == b'.' && bytes[k + 1] == b'.' {
            report.push(Diagnostic::error(
                Code("SYN_MULTIPLE_DOTS"),
                Category::Syn,
                Span::new(k, k + 2),
                "Multiple dots",
            ));
            break;
        }
        k += 1;
    }
}

pub fn parse_bracket_inner(inner: &str) -> BracketParsed {
    fn isotope<'a>() -> impl Parser<&'a str, Output = u32, Error = Error<&'a str>> {
        map(digit1, |digits: &str| digits.parse::<u32>().unwrap_or(0))
    }

    fn element_symbol<'a>() -> impl Parser<&'a str, Output = Option<Element>, Error = Error<&'a str>>
    {
        alt((
            value(None, tag("*")),
            map(
                map_opt(
                    recognize(pair(
                        one_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ"),
                        opt(one_of("abcdefghijklmnopqrstuvwxyz")),
                    )),
                    |s: &str| Element::from_symbol(s),
                ),
                Some,
            ),
        ))
    }

    fn chiral<'a>() -> impl Parser<&'a str, Output = BracketField, Error = Error<&'a str>> {
        alt((
            value(BracketField::Chiral(Chirality::CounterClockwise), tag("@@")),
            map(preceded(tag("@TH"), digit1), |d: &str| {
                let n = d.as_bytes()[0] - b'0';
                BracketField::Chiral(Chirality::Tetrahedral { arr: n as u32 })
            }),
            map(preceded(tag("@AL"), one_of("12")), |c: char| {
                BracketField::Chiral(Chirality::Allenal {
                    arr: c.to_digit(10).unwrap() as u32,
                })
            }),
            map(preceded(tag("@SP"), one_of("123")), |c: char| {
                BracketField::Chiral(Chirality::SquarePlanar {
                    arr: c.to_digit(10).unwrap() as u32,
                })
            }),
            map(preceded(tag("@TB"), digit1), |d: &str| {
                let n = d.as_bytes()[0] - b'0';
                BracketField::Chiral(Chirality::TrigonalBipyramidal { arr: n as u32 })
            }),
            map(preceded(tag("@OH"), digit1), |d: &str| {
                let n = d.as_bytes()[0] - b'0';
                BracketField::Chiral(Chirality::Octahedral { arr: n as u32 })
            }),
            value(BracketField::Chiral(Chirality::Clockwise), tag("@")),
        ))
    }

    fn d1<'a>() -> impl Parser<&'a str, Output = u32, Error = Error<&'a str>> {
        map(one_of(DIGITS), |c| c.to_digit(10).unwrap())
    }

    fn d1_to_2<'a>() -> impl Parser<&'a str, Output = u32, Error = Error<&'a str>> {
        map(pair(one_of(DIGITS), opt(one_of(DIGITS))), |(d1, d2)| {
            let mut v = d1.to_digit(10).unwrap();
            if let Some(c2) = d2 {
                v = v * 10 + c2.to_digit(10).unwrap();
            }
            v
        })
    }

    fn hcount<'a>() -> impl Parser<&'a str, Output = BracketField, Error = Error<&'a str>> {
        map(pair(tag("H"), opt(d1())), |(_, d): (&str, Option<u32>)| {
            BracketField::H(d.unwrap_or(1))
        })
    }

    fn charge<'a>() -> impl Parser<&'a str, Output = BracketField, Error = Error<&'a str>> {
        alt((
            value(BracketField::Q(2), tag("++")),
            value(BracketField::Q(-2), tag("--")),
            map(
                pair(tag("+"), opt(d1_to_2())),
                |(_, n): (&str, Option<u32>)| BracketField::Q(n.unwrap_or(1) as i32),
            ),
            map(
                pair(tag("-"), opt(d1_to_2())),
                |(_, n): (&str, Option<u32>)| BracketField::Q(-(n.unwrap_or(1) as i32)),
            ),
        ))
    }

    fn class<'a>() -> impl Parser<&'a str, Output = BracketField, Error = Error<&'a str>> {
        map(preceded(tag(":"), digit1), |digits: &str| {
            BracketField::Class(digits.parse::<u32>().unwrap_or(0))
        })
    }

    fn field_tail<'a>() -> impl Parser<&'a str, Output = BracketField, Error = Error<&'a str>> {
        alt((chiral(), hcount(), charge(), class()))
    }

    fn fields<'a>() -> impl Parser<&'a str, Output = BracketParsed, Error = Error<&'a str>> {
        map(
            (
                opt(isotope()),
                element_symbol(),
                fold_many0::<&str, Error<&str>, _, _, _, SmallVec<[BracketField; 4]>>(
                    field_tail(),
                    || SmallVec::<[BracketField; 4]>::new(),
                    |mut acc, f| {
                        acc.push(f);
                        acc
                    },
                ),
            ),
            |(iso_opt, elem_opt, tails): (
                Option<u32>,
                Option<Element>,
                SmallVec<[BracketField; 4]>,
            )| {
                let mut parsed = BracketParsed::default();
                if let Some(i) = iso_opt {
                    parsed.isotope = Some(i);
                }
                parsed.element = elem_opt;
                for f in tails {
                    match f {
                        BracketField::Chiral(_) => {}
                        BracketField::H(v) => parsed.hcount = Some(v),
                        BracketField::Q(q) => parsed.charge = Some(q),
                        BracketField::Class(c) => parsed.class = Some(c),
                    }
                }
                parsed
            },
        )
    }

    fields().parse(inner).map(|(_, f)| f).unwrap_or_default()
}
