//! Configurable EDN formatter (pretty-printing).

use crate::edn::Edn;

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
    /// Sort map keys before output. No-op for `BTreeMap` (already sorted).
    pub sort_maps: bool,
    /// Sort set elements before output. No-op for `BTreeSet` (already sorted).
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
            let items: Vec<&Edn<'_>> = s.iter().collect();
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

fn fits_on_line(fmt: &EdnFormatter, depth: usize, content_len: usize) -> bool {
    match fmt.line_width {
        Some(width) => {
            let indent_len = depth * fmt.indent.len();
            indent_len + content_len <= width
        }
        None => false,
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

    if fmt.compact_seqs {
        let compact_parts: Vec<String> = items.iter().map(|e| e.to_string()).collect();
        let compact_len = open.len()
            + close.len()
            + compact_parts.iter().map(|s| s.len()).sum::<usize>()
            + compact_parts.len().saturating_sub(1);
        if fits_on_line(fmt, depth, compact_len) {
            out.push_str(open);
            for (i, s) in compact_parts.iter().enumerate() {
                if i > 0 {
                    out.push(' ');
                }
                out.push_str(s);
            }
            out.push_str(close);
            return;
        }
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
    m: &std::collections::BTreeMap<Edn<'_>, Edn<'_>>,
    fmt: &EdnFormatter,
    depth: usize,
) {
    if m.is_empty() {
        out.push_str("{}");
        return;
    }

    let separator = if fmt.commas { ", " } else { " " };

    if fmt.compact_maps {
        let compact_parts: Vec<String> = m.iter().map(|(k, v)| format!("{k} {v}")).collect();
        let compact_len = 2 + compact_parts.iter().map(|s| s.len()).sum::<usize>()
            + separator.len() * compact_parts.len().saturating_sub(1);
        if fits_on_line(fmt, depth, compact_len) {
            out.push('{');
            for (i, s) in compact_parts.iter().enumerate() {
                if i > 0 {
                    out.push_str(separator);
                }
                out.push_str(s);
            }
            out.push('}');
            return;
        }
    }

    out.push('{');
    for (i, (k, v)) in m.iter().enumerate() {
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

    if fmt.compact_sets {
        let compact_parts: Vec<String> = items.iter().map(|e| e.to_string()).collect();
        let compact_len = 3 + compact_parts.iter().map(|s| s.len()).sum::<usize>()
            + compact_parts.len().saturating_sub(1); // #{ } + spaces
        if fits_on_line(fmt, depth, compact_len) {
            out.push_str("#{");
            for (i, s) in compact_parts.iter().enumerate() {
                if i > 0 {
                    out.push(' ');
                }
                out.push_str(s);
            }
            out.push('}');
            return;
        }
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
pub fn to_string_pretty<T: serde::Serialize>(value: &T) -> Result<String, crate::error::EdnError> {
    to_string_pretty_with(value, &EdnFormatter::default())
}

/// Serialize a serde value to an EDN string with custom formatting.
#[cfg(feature = "serde")]
pub fn to_string_pretty_with<T: serde::Serialize>(
    value: &T,
    fmt: &EdnFormatter,
) -> Result<String, crate::error::EdnError> {
    let compact = crate::ser::to_string(value)?;
    let edn = crate::reader::read_string(&compact)?;
    Ok(edn.to_string_with(fmt))
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;
    use std::collections::BTreeMap;

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
        let mut m = BTreeMap::new();
        m.insert(Edn::keyword("a"), Edn::Int(1));
        let result = Edn::Map(m).to_string_with(&fmt_default());
        assert_eq!(result, "{:a 1}");
    }

    #[test]
    fn test_formatter_long_map_multiline() {
        let mut m = BTreeMap::new();
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
        let mut inner = BTreeMap::new();
        inner.insert(Edn::keyword("x"), Edn::Int(1));
        inner.insert(Edn::keyword("y"), Edn::Int(2));
        let mut outer = BTreeMap::new();
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
        assert_eq!(Edn::Map(BTreeMap::new()).to_string_with(&fmt_default()), "{}");
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
        let mut m = BTreeMap::new();
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
        let mut m = BTreeMap::new();
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
        let mut m = BTreeMap::new();
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
        let edn = crate::read_string(input).unwrap();
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
