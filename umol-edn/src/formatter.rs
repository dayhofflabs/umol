//! Configurable EDN formatter (pretty-printing).

use crate::edn::Edn;

/// Configurable EDN formatter.
#[derive(Clone, Debug)]
pub struct EdnFormatter {
    /// Indent string per nesting level (default: two spaces).
    pub indent: String,
    /// Target line width before wrapping collections. `None` disables wrapping
    /// (always multi-line for non-empty collections).
    pub line_width: Option<usize>,
}

impl Default for EdnFormatter {
    fn default() -> Self {
        Self {
            indent: "  ".to_string(),
            line_width: Some(80),
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
            // Atoms: use Display (compact)
            out.push_str(&other.to_string());
        }
    }
}

fn fits_on_line(fmt: &EdnFormatter, depth: usize, content_len: usize, overhead: usize) -> bool {
    match fmt.line_width {
        Some(width) => {
            let indent_len = depth * fmt.indent.len();
            indent_len + overhead + content_len <= width
        }
        None => false,
    }
}

fn write_indent(out: &mut String, fmt: &EdnFormatter, depth: usize) {
    for _ in 0..depth {
        out.push_str(&fmt.indent);
    }
}

fn write_seq(out: &mut String, open: &str, close: &str, items: &[Edn<'_>], fmt: &EdnFormatter, depth: usize) {
    if items.is_empty() {
        out.push_str(open);
        out.push_str(close);
        return;
    }

    // Try compact representation
    let compact = format!("{}", Edn::Vector(items.to_vec()));
    // Replace the vector delimiters with the actual ones
    let compact = if open == "(" {
        format!("({})", &compact[1..compact.len() - 1])
    } else {
        compact
    };
    if fits_on_line(fmt, depth, compact.len(), 0) {
        out.push_str(&compact);
        return;
    }

    // Multi-line
    out.push_str(open);
    for (i, item) in items.iter().enumerate() {
        if i > 0 {
            out.push('\n');
            write_indent(out, fmt, depth + 1);
        }
        write_edn(out, item, fmt, depth + 1);
    }
    out.push_str(close);
}

fn write_map(out: &mut String, m: &std::collections::BTreeMap<Edn<'_>, Edn<'_>>, fmt: &EdnFormatter, depth: usize) {
    if m.is_empty() {
        out.push_str("{}");
        return;
    }

    // Try compact
    let edn = Edn::Map(m.clone());
    let compact = edn.to_string();
    if fits_on_line(fmt, depth, compact.len(), 0) {
        out.push_str(&compact);
        return;
    }

    // Multi-line: one key-value pair per line
    out.push('{');
    for (i, (k, v)) in m.iter().enumerate() {
        if i > 0 {
            out.push('\n');
            write_indent(out, fmt, depth + 1);
        }
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

    // Try compact
    let compact_parts: Vec<String> = items.iter().map(|e| e.to_string()).collect();
    let compact_len = 2 + compact_parts.iter().map(|s| s.len()).sum::<usize>()
        + compact_parts.len().saturating_sub(1); // #{...} + spaces
    if fits_on_line(fmt, depth, compact_len, 0) {
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

    // Multi-line
    out.push_str("#{");
    for (i, item) in items.iter().enumerate() {
        if i > 0 {
            out.push('\n');
            write_indent(out, fmt, depth + 1);
        }
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
    // Serialize to compact EDN, parse back as Edn, then pretty-print.
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
        // Each key-value pair on its own line
        let lines: Vec<&str> = result.lines().collect();
        assert_eq!(lines.len(), 5);
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
        // Inner map should also be formatted
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
        // With None line_width, always multi-line for non-empty
        assert!(result.contains('\n'));
    }

    #[test]
    fn test_formatter_custom_indent() {
        let fmt = EdnFormatter {
            indent: "    ".to_string(),
            line_width: Some(10),
        };
        let v = Edn::Vector(vec![
            Edn::Int(100), Edn::Int(200), Edn::Int(300), Edn::Int(400),
        ]);
        let result = v.to_string_with(&fmt);
        assert!(result.contains("    ")); // 4-space indent
    }

    #[rstest]
    #[case(
        "{:atoms [\"C\" \"O\"] :bonds [[:0 :1 :single]]}",
        20,
        // Should go multi-line at narrow width
        true
    )]
    #[case(
        "{:a 1}",
        80,
        // Short map, fits on one line
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
            // Small struct fits on one line
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
