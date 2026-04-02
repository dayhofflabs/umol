//! Configurable EDN formatter (pretty-printing).

use crate::edn::{Edn, EdnMap};
use crate::error::EdnError;
use crate::reader::read_string;
use crate::ser::to_string;

/// Configurable EDN formatter.
#[derive(Clone, Debug)]
pub struct EdnFormatter {
    /// Indent string per nesting level.
    pub indent: String,
    /// Newline string (`"\n"` or `"\r\n"`).
    pub newline: String,
    /// Target line width before wrapping collections. `None` = always expand.
    pub line_width: Option<usize>,
    /// Insert commas between map entries.
    pub commas: bool,
    /// Allow single-line maps when they fit within `line_width`.
    pub compact_maps: bool,
    /// Allow single-line lists/vectors when they fit within `line_width`.
    pub compact_seqs: bool,
    /// Allow single-line sets when they fit within `line_width`.
    pub compact_sets: bool,
    /// Sort map keys before output.
    pub sort_maps: bool,
    /// Sort set elements before output.
    pub sort_sets: bool,
}

impl Default for EdnFormatter {
    fn default() -> Self {
        Self {
            indent: "  ".to_string(),
            newline: "\n".to_string(),
            line_width: Some(80),
            commas: false,
            compact_maps: true,
            compact_seqs: true,
            compact_sets: true,
            sort_maps: true,
            sort_sets: true,
        }
    }
}

// TODO: Check if this should be moved to edn.rs
impl<'a> Edn<'a> {
    /// Format this value using the given formatter.
    pub fn to_string_with(&self, fmt: &EdnFormatter) -> String {
        let mut out = String::new();
        write_edn(&mut out, self, fmt, 0);
        out
    }
}

fn write_edn(out: &mut String, edn: &Edn<'_>, fmt: &EdnFormatter, depth: usize) {
    match edn {
        Edn::List(items) => write_seq(out, "(", ")", items, fmt, depth),
        Edn::Vector(items) => write_seq(out, "[", "]", items, fmt, depth),
        Edn::Map(m) => write_map(out, m, fmt, depth),
        Edn::Set(s) => {
            let mut items: Vec<&Edn<'_>> = s.iter().collect();
            if fmt.sort_sets {
                items.sort();
            }
            write_set(out, &items, fmt, depth);
        }
        Edn::Tagged(tag, inner) => {
            out.push('#');
            out.push_str(tag);
            out.push(' ');
            write_edn(out, inner, fmt, depth);
        }
        other => {
            out.push_str(&other.to_string());
        }
    }
}

/// Compute the compact (single-line) display length of an EDN value.
/// Returns `None` if it exceeds `limit` (early exit to avoid full traversal).
fn compact_len(edn: &Edn<'_>, limit: usize) -> Option<usize> {
    match edn {
        Edn::Nil => Some(3),
        Edn::Bool(true) => Some(4),
        Edn::Bool(false) => Some(5),
        Edn::Int(n) => Some(itoa_len(*n)),
        Edn::Float(v) => Some(format_float_len(*v)),
        Edn::Char(c) => Some(display_char_len(*c)),
        Edn::Str(s) => Some(display_string_len(s)),
        Edn::Keyword(k) => Some(1 + k.as_str().len()), // :name
        Edn::Symbol(s) => Some(s.as_str().len()),
        Edn::List(items) => seq_compact_len(items, 2, limit), // ()
        Edn::Vector(items) => seq_compact_len(items, 2, limit), // []
        Edn::Map(m) => map_compact_len(m, limit),
        Edn::Set(s) => {
            let items: Vec<&Edn<'_>> = s.iter().collect();
            set_compact_len(&items, limit)
        }
        Edn::Tagged(tag, inner) => {
            let prefix = 1 + tag.len() + 1; // # + tag + space
            compact_len(inner, limit.checked_sub(prefix)?).map(|n| prefix + n)
        }
    }
}

fn itoa_len(n: i64) -> usize {
    let mut buf = itoa::Buffer::new();
    buf.format(n).len()
}

fn format_float_len(v: f64) -> usize {
    if v.is_nan() {
        return 5; // ##NaN
    }
    if v == f64::INFINITY {
        return 5; // ##Inf
    }
    if v == f64::NEG_INFINITY {
        return 6; // ##-Inf
    }
    let s = format!("{v}");
    if s.contains('.') || s.contains('e') || s.contains('E') {
        s.len()
    } else {
        s.len() + 2 // append .0
    }
}

fn display_char_len(c: char) -> usize {
    match c {
        '\n' => 8,  // \newline
        '\r' => 7,  // \return
        ' ' => 6,   // \space
        '\t' => 4,  // \tab
        '\u{000C}' => 9,  // \formfeed
        '\u{0008}' => 10, // \backspace
        _ => 1 + c.len_utf8(), // \ + char
    }
}

fn display_string_len(s: &str) -> usize {
    let mut len = 2; // opening and closing quotes
    for c in s.chars() {
        len += match c {
            '"' | '\\' | '\n' | '\r' | '\t' | '\u{0008}' | '\u{000C}' => 2,
            _ => c.len_utf8(),
        };
    }
    len
}

fn seq_compact_len(items: &[Edn<'_>], overhead: usize, limit: usize) -> Option<usize> {
    if items.is_empty() {
        return Some(overhead);
    }
    let mut total = overhead; // delimiters
    for (i, item) in items.iter().enumerate() {
        if i > 0 {
            total += 1; // space
        }
        total += compact_len(item, limit.checked_sub(total)?)?;
        if total > limit {
            return None;
        }
    }
    Some(total)
}

fn map_compact_len(
    m: &EdnMap<'_>,
    limit: usize,
) -> Option<usize> {
    if m.is_empty() {
        return Some(2);
    }
    let mut total = 2usize; // { }
    for (i, (k, v)) in m.iter().enumerate() {
        if i > 0 {
            total += 1; // space between pairs
        }
        total += compact_len(k, limit.checked_sub(total)?)?;
        total += 1; // space between key and value
        total += compact_len(v, limit.checked_sub(total)?)?;
        if total > limit {
            return None;
        }
    }
    Some(total)
}

fn set_compact_len(items: &[&Edn<'_>], limit: usize) -> Option<usize> {
    if items.is_empty() {
        return Some(3); // #{}
    }
    let mut total = 3usize; // #{ }
    for (i, item) in items.iter().enumerate() {
        if i > 0 {
            total += 1; // space
        }
        total += compact_len(item, limit.checked_sub(total)?)?;
        if total > limit {
            return None;
        }
    }
    Some(total)
}

fn fits_compact(fmt: &EdnFormatter, depth: usize, edn_len: Option<usize>) -> bool {
    match (fmt.line_width, edn_len) {
        (Some(width), Some(len)) => depth * fmt.indent.len() + len <= width,
        _ => false,
    }
}

fn write_newline(out: &mut String, fmt: &EdnFormatter, depth: usize) {
    out.push_str(&fmt.newline);
    for _ in 0..depth {
        out.push_str(&fmt.indent);
    }
}

fn write_seq(
    out: &mut String,
    open: &str,
    close: &str,
    items: &[Edn<'_>],
    fmt: &EdnFormatter,
    depth: usize,
) {
    if items.is_empty() {
        out.push_str(open);
        out.push_str(close);
        return;
    }

    let limit = fmt.line_width.unwrap_or(0);
    if fmt.compact_seqs && fits_compact(fmt, depth, seq_compact_len(items, 2, limit)) {
        out.push_str(open);
        for (i, item) in items.iter().enumerate() {
            if i > 0 {
                out.push(' ');
            }
            out.push_str(&item.to_string());
        }
        out.push_str(close);
        return;
    }

    out.push_str(open);
    for item in items.iter() {
        write_newline(out, fmt, depth + 1);
        write_edn(out, item, fmt, depth + 1);
    }
    out.push_str(close);
}

fn write_map(
    out: &mut String,
    m: &EdnMap<'_>,
    fmt: &EdnFormatter,
    depth: usize,
) {
    if m.is_empty() {
        out.push_str("{}");
        return;
    }

    let mut entries: Vec<_> = m.iter().collect();
    if fmt.sort_maps {
        entries.sort_by(|a, b| a.0.cmp(b.0));
    }

    let separator = if fmt.commas { ", " } else { " " };
    let limit = fmt.line_width.unwrap_or(0);

    let comma_extra = if fmt.commas {
        m.len().saturating_sub(1)
    } else {
        0
    };
    if fmt.compact_maps
        && fits_compact(
            fmt,
            depth,
            map_compact_len(m, limit).map(|n| n + comma_extra),
        )
    {
        out.push('{');
        for (i, (k, v)) in entries.iter().enumerate() {
            if i > 0 {
                out.push_str(separator);
            }
            out.push_str(&k.to_string());
            out.push(' ');
            out.push_str(&v.to_string());
        }
        out.push('}');
        return;
    }

    out.push('{');
    for (i, (k, v)) in entries.iter().enumerate() {
        if i > 0 && fmt.commas {
            out.push(',');
        }
        write_newline(out, fmt, depth + 1);
        write_edn(out, k, fmt, depth + 1);
        out.push(' ');
        write_edn(out, v, fmt, depth + 1);
    }
    out.push('}');
}

fn write_set(out: &mut String, items: &[&Edn<'_>], fmt: &EdnFormatter, depth: usize) {
    if items.is_empty() {
        out.push_str("#{}");
        return;
    }

    let limit = fmt.line_width.unwrap_or(0);
    if fmt.compact_sets && fits_compact(fmt, depth, set_compact_len(items, limit)) {
        out.push_str("#{");
        for (i, item) in items.iter().enumerate() {
            if i > 0 {
                out.push(' ');
            }
            out.push_str(&item.to_string());
        }
        out.push('}');
        return;
    }

    out.push_str("#{");
    for item in items.iter() {
        write_newline(out, fmt, depth + 1);
        write_edn(out, item, fmt, depth + 1);
    }
    out.push('}');
}

/// Serialize a serde value to a pretty-printed EDN string.
#[cfg(feature = "serde")]
pub fn to_string_pretty<T: serde::Serialize>(value: &T) -> Result<String, EdnError> {
    to_string_pretty_with(value, &EdnFormatter::default())
}

/// Serialize a serde value to an EDN string with custom formatting.
#[cfg(feature = "serde")]
pub fn to_string_pretty_with<T: serde::Serialize>(
    value: &T,
    fmt: &EdnFormatter,
) -> Result<String, EdnError> {
    let compact = to_string(value)?;
    let edn = read_string(&compact)?;
    Ok(edn.to_string_with(fmt))
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;
    use crate::edn::EdnMap;

    use rstest::rstest;

    use crate::edn::{Edn, Keyword};
    use super::*;

    fn fmt_default() -> EdnFormatter {
        EdnFormatter::default()
    }

    fn fmt_narrow() -> EdnFormatter {
        EdnFormatter {
            line_width: Some(20),
            ..Default::default()
        }
    }

    #[rstest]
    #[case(Edn::Nil, "nil")]
    #[case(Edn::Int(12), "12")]
    #[case(Edn::Bool(true), "true")]
    #[case(Edn::Str(Cow::Borrowed("hello")), r#""hello""#)]
    fn test_formatter_atoms(#[case] edn: Edn<'_>, #[case] expected: &str) {
        assert_eq!(edn.to_string_with(&fmt_default()), expected);
    }

    #[test]
    fn test_formatter_short_vector_compact() {
        let v = Edn::Vector(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)]);
        assert_eq!(v.to_string_with(&fmt_default()), "[1 2 3]");
    }

    #[test]
    fn test_formatter_long_vector_multiline() {
        let items: Vec<Edn<'_>> = (0..20).map(|i| Edn::Str(Cow::Owned(format!("item-{i}")))).collect();
        let v = Edn::Vector(items);
        let result = v.to_string_with(&fmt_narrow());
        assert!(result.contains('\n'));
        assert!(result.starts_with('['));
        assert!(result.ends_with(']'));
    }

    #[test]
    fn test_formatter_short_map_compact() {
        let mut m = EdnMap::new();
        m.insert(Edn::keyword("a"), Edn::Int(1));
        let result = Edn::Map(m).to_string_with(&fmt_default());
        assert_eq!(result, "{:a 1}");
    }

    #[test]
    fn test_formatter_long_map_multiline() {
        let mut m = EdnMap::new();
        for i in 0..5 {
            m.insert(
                Edn::Keyword(Keyword::owned(format!("key-{i}"))),
                Edn::Str(Cow::Owned(format!("value-{i}"))),
            );
        }
        let result = Edn::Map(m).to_string_with(&fmt_narrow());
        assert!(result.contains('\n'));
        assert!(result.starts_with('{'));
        assert!(result.ends_with('}'));
        let lines: Vec<&str> = result.lines().collect();
        assert_eq!(lines.len(), 6); // { + 5 entries on separate lines
    }

    #[test]
    fn test_formatter_nested_map() {
        let mut inner = EdnMap::new();
        inner.insert(Edn::keyword("x"), Edn::Int(1));
        inner.insert(Edn::keyword("y"), Edn::Int(2));
        let mut outer = EdnMap::new();
        outer.insert(Edn::keyword("point"), Edn::Map(inner));
        outer.insert(Edn::keyword("label"), Edn::Str(Cow::Borrowed("origin")));
        let result = Edn::Map(outer).to_string_with(&fmt_narrow());
        assert!(result.contains('\n'));
        assert!(result.contains(":point"));
        assert!(result.contains(":label"));
    }

    #[test]
    fn test_formatter_empty_collections() {
        assert_eq!(Edn::Vector(vec![]).to_string_with(&fmt_default()), "[]");
        assert_eq!(Edn::List(vec![]).to_string_with(&fmt_default()), "()");
        assert_eq!(Edn::Map(EdnMap::new()).to_string_with(&fmt_default()), "{}");
    }

    #[test]
    fn test_formatter_no_line_width() {
        let fmt = EdnFormatter {
            line_width: None,
            ..Default::default()
        };
        let v = Edn::Vector(vec![Edn::Int(1), Edn::Int(2)]);
        let result = v.to_string_with(&fmt);
        assert!(result.contains('\n'));
    }

    #[test]
    fn test_formatter_custom_indent() {
        let fmt = EdnFormatter {
            indent: "    ".to_string(),
            line_width: Some(10),
            ..Default::default()
        };
        let v = Edn::Vector(vec![
            Edn::Int(100), Edn::Int(200), Edn::Int(300), Edn::Int(400),
        ]);
        let result = v.to_string_with(&fmt);
        assert!(result.contains("    "));
    }

    #[test]
    fn test_formatter_commas() {
        let mut m = EdnMap::new();
        m.insert(Edn::keyword("a"), Edn::Int(1));
        m.insert(Edn::keyword("b"), Edn::Int(2));
        let fmt = EdnFormatter {
            commas: true,
            ..Default::default()
        };
        let result = Edn::Map(m).to_string_with(&fmt);
        assert_eq!(result, "{:a 1, :b 2}");
    }

    #[test]
    fn test_formatter_commas_multiline() {
        let mut m = EdnMap::new();
        for i in 0..5 {
            m.insert(
                Edn::Keyword(Keyword::owned(format!("key-{i}"))),
                Edn::Str(Cow::Owned(format!("value-{i}"))),
            );
        }
        let fmt = EdnFormatter {
            commas: true,
            line_width: Some(20),
            ..Default::default()
        };
        let result = Edn::Map(m).to_string_with(&fmt);
        assert!(result.contains(",\n"));
    }

    #[test]
    fn test_formatter_no_compact_seqs() {
        let fmt = EdnFormatter {
            compact_seqs: false,
            ..Default::default()
        };
        let v = Edn::Vector(vec![Edn::Int(1), Edn::Int(2)]);
        let result = v.to_string_with(&fmt);
        assert!(result.contains('\n'));
    }

    #[test]
    fn test_formatter_no_compact_maps() {
        let fmt = EdnFormatter {
            compact_maps: false,
            ..Default::default()
        };
        let mut m = EdnMap::new();
        m.insert(Edn::keyword("a"), Edn::Int(1));
        let result = Edn::Map(m).to_string_with(&fmt);
        assert!(result.contains('\n'));
    }

    #[test]
    fn test_formatter_crlf() {
        let fmt = EdnFormatter {
            newline: "\r\n".to_string(),
            line_width: None,
            ..Default::default()
        };
        let v = Edn::Vector(vec![Edn::Int(1), Edn::Int(2)]);
        let result = v.to_string_with(&fmt);
        assert!(result.contains("\r\n"));
    }

    #[rstest]
    #[case(
        "{:atoms [\"C\" \"O\"] :bonds [[:0 :1 :single]]}",
        20,
        true
    )]
    #[case(
        "{:a 1}",
        80,
        false
    )]
    fn test_formatter_molecule_like(#[case] input: &str, #[case] width: usize, #[case] expect_multiline: bool) {
        let edn = read_string(input).unwrap();
        let fmt = EdnFormatter {
            line_width: Some(width),
            ..Default::default()
        };
        let result = edn.to_string_with(&fmt);
        assert_eq!(result.contains('\n'), expect_multiline);
    }

    #[cfg(feature = "serde")]
    mod serde_tests {
        use serde::Serialize;

        use super::*;

        #[derive(Serialize)]
        struct Point {
            x: f64,
            y: f64,
        }

        #[derive(Serialize)]
        struct Molecule {
            atoms: Vec<String>,
            bonds: Vec<(String, String, String)>,
        }

        #[test]
        fn test_to_string_pretty_struct() {
            let p = Point { x: 1.0, y: 2.0 };
            let result = to_string_pretty(&p).unwrap();
            assert_eq!(result, "{:x 1.0 :y 2.0}");
        }

        #[test]
        fn test_to_string_pretty_with_narrow() {
            let mol = Molecule {
                atoms: vec!["C".into(), "O".into(), "H".into()],
                bonds: vec![("0".into(), "1".into(), "single".into())],
            };
            let fmt = EdnFormatter {
                line_width: Some(30),
                ..Default::default()
            };
            let result = to_string_pretty_with(&mol, &fmt).unwrap();
            assert!(result.contains('\n'));
        }
    }
}
