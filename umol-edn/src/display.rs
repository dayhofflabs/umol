//! Display implementation for Edn values (compact EDN format).

use std::fmt;

use crate::edn::Edn;
use crate::parser::{is_symbol_char, is_symbol_start, validate_symbol};

impl fmt::Display for Edn<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Edn::Nil => write!(f, "nil"),
            Edn::Bool(true) => write!(f, "true"),
            Edn::Bool(false) => write!(f, "false"),
            Edn::Int(n) => write!(f, "{n}"),
            #[cfg(feature = "bignum")]
            Edn::BigInt(n) => write!(f, "{n}N"),
            Edn::Float(v) => format_float(f, *v),
            #[cfg(feature = "bignum")]
            Edn::BigDecimal(d) => write!(f, "{d}M"),
            Edn::Char(c) => format_char(f, *c),
            Edn::Str(s) => format_string(f, s),
            Edn::Keyword(k) => write!(f, "{k}"),
            Edn::Symbol(s) => write!(f, "{s}"),
            Edn::List(items) => format_seq(f, "(", ")", items),
            Edn::Vector(items) => format_seq(f, "[", "]", items),
            Edn::Map(m) => {
                write!(f, "{{")?;
                let mut entries: Vec<_> = m.iter().collect();
                entries.sort_by_key(|(k, _)| *k);
                for (i, (k, v)) in entries.iter().enumerate() {
                    if i > 0 {
                        write!(f, " ")?;
                    }
                    write!(f, "{k} {v}")?;
                }
                write!(f, "}}")
            }
            Edn::Set(s) => {
                write!(f, "#{{")?;
                let mut elems: Vec<_> = s.iter().collect();
                elems.sort();
                for (i, v) in elems.iter().enumerate() {
                    if i > 0 {
                        write!(f, " ")?;
                    }
                    write!(f, "{v}")?;
                }
                write!(f, "}}")
            }
            Edn::Tagged(tag, inner) => {
                if !is_valid_tag(tag) {
                    return Err(fmt::Error);
                }
                write!(f, "#{tag} {inner}")
            }
        }
    }
}

fn is_valid_tag(s: &str) -> bool {
    let mut chars = s.chars();
    match chars.next() {
        None => return false,
        Some(c) if !is_symbol_start(c) => return false,
        _ => {}
    }
    if !chars.all(is_symbol_char) {
        return false;
    }
    validate_symbol(s, 0).is_ok()
}

fn format_float(f: &mut fmt::Formatter<'_>, v: f64) -> fmt::Result {
    if v.is_nan() || !v.is_finite() {
        return Err(fmt::Error);
    }
    let mut buf = zmij::Buffer::new();
    f.write_str(buf.format_finite(v))
}

fn format_char(f: &mut fmt::Formatter<'_>, c: char) -> fmt::Result {
    match c {
        '\n' => write!(f, "\\newline"),
        '\r' => write!(f, "\\return"),
        ' ' => write!(f, "\\space"),
        '\t' => write!(f, "\\tab"),
        c if (c as u32) < 0x20 || c == '\u{7F}' => write!(f, "\\u{:04X}", c as u32),
        _ => write!(f, "\\{c}"),
    }
}

fn format_string(f: &mut fmt::Formatter<'_>, s: &str) -> fmt::Result {
    write!(f, "\"")?;
    for c in s.chars() {
        match c {
            '"' => write!(f, "\\\"")?,
            '\\' => write!(f, "\\\\")?,
            '\n' => write!(f, "\\n")?,
            '\r' => write!(f, "\\r")?,
            '\t' => write!(f, "\\t")?,
            _ => write!(f, "{c}")?,
        }
    }
    write!(f, "\"")
}

fn format_seq(
    f: &mut fmt::Formatter<'_>,
    open: &str,
    close: &str,
    items: &[Edn<'_>],
) -> fmt::Result {
    write!(f, "{open}")?;
    for (i, item) in items.iter().enumerate() {
        if i > 0 {
            write!(f, " ")?;
        }
        write!(f, "{item}")?;
    }
    write!(f, "{close}")
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use rstest::rstest;

    use crate::collections::{EdnMap, EdnSet};
    use crate::edn::{Edn, EdnKeyword, EdnSymbol};

    #[rstest]
    #[case(Edn::Nil, "nil")]
    #[case(Edn::Bool(true), "true")]
    #[case(Edn::Bool(false), "false")]
    #[case(Edn::Int(12), "12")]
    #[case(Edn::Int(-1), "-1")]
    #[case(Edn::Int(0), "0")]
    fn test_display_primitives(#[case] edn: Edn<'_>, #[case] expected: &str) {
        assert_eq!(edn.to_string(), expected);
    }

    #[rstest]
    #[case(3.14, "3.14")]
    #[case(1.0, "1.0")]
    #[case(12.0, "12.0")]
    #[case(-0.5, "-0.5")]
    fn test_display_float(#[case] v: f64, #[case] expected: &str) {
        assert_eq!(Edn::Float(v).to_string(), expected);
    }

    #[rstest]
    #[case('a', "\\a")]
    #[case('Z', "\\Z")]
    #[case('\n', "\\newline")]
    #[case('\r', "\\return")]
    #[case(' ', "\\space")]
    #[case('\t', "\\tab")]
    #[case('\u{000C}', "\\u000C")]
    #[case('\u{0008}', "\\u0008")]
    fn test_display_char(#[case] c: char, #[case] expected: &str) {
        assert_eq!(Edn::Char(c).to_string(), expected);
    }

    #[rstest]
    #[case("hello", r#""hello""#)]
    #[case("with \"quotes\"", r#""with \"quotes\"""#)]
    #[case("line\nbreak", r#""line\nbreak""#)]
    #[case("tab\there", r#""tab\there""#)]
    #[case("back\\slash", r#""back\\slash""#)]
    fn test_display_string(#[case] s: &str, #[case] expected: &str) {
        assert_eq!(Edn::Str(Cow::Borrowed(s)).to_string(), expected);
    }

    #[rstest]
    #[case(EdnKeyword::new("foo"), ":foo")]
    #[case(EdnKeyword::namespaced("ns", "bar"), ":ns/bar")]
    fn test_display_keyword(#[case] k: EdnKeyword<'_>, #[case] expected: &str) {
        assert_eq!(Edn::Keyword(k).to_string(), expected);
    }

    #[rstest]
    #[case(EdnSymbol::new("foo"), "foo")]
    #[case(EdnSymbol::namespaced("ns", "bar"), "ns/bar")]
    fn test_display_symbol(#[case] s: EdnSymbol<'_>, #[case] expected: &str) {
        assert_eq!(Edn::Symbol(s).to_string(), expected);
    }

    #[rstest]
    #[case(Edn::List(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)].into()), "(1 2 3)")]
    #[case(Edn::List(vec![].into()), "()")]
    #[case(Edn::Vector(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)].into()), "[1 2 3]")]
    #[case(Edn::Vector(vec![].into()), "[]")]
    fn test_display_seq(#[case] edn: Edn<'_>, #[case] expected: &str) {
        assert_eq!(edn.to_string(), expected);
    }

    #[test]
    fn test_display_map() {
        let mut m = EdnMap::new();
        m.insert(Edn::Keyword(EdnKeyword::new("a")), Edn::Int(1));
        assert_eq!(Edn::Map(m).to_string(), "{:a 1}");
    }

    #[test]
    fn test_display_map_multi() {
        let mut m = EdnMap::new();
        m.insert(Edn::Keyword(EdnKeyword::new("a")), Edn::Int(1));
        m.insert(Edn::Keyword(EdnKeyword::new("b")), Edn::Int(2));
        let s = Edn::Map(m).to_string();
        assert!(s == "{:a 1 :b 2}" || s == "{:b 2 :a 1}");
    }

    #[test]
    fn test_display_map_empty() {
        assert_eq!(Edn::Map(EdnMap::new()).to_string(), "{}");
    }

    #[test]
    fn test_display_set() {
        let mut s = EdnSet::new();
        s.insert(Edn::Int(1));
        s.insert(Edn::Int(2));
        let result = Edn::Set(s).to_string();
        assert!(result == "#{1 2}" || result == "#{2 1}");
    }

    #[test]
    fn test_display_set_empty() {
        assert_eq!(Edn::Set(EdnSet::new()).to_string(), "#{}");
    }

    #[test]
    fn test_display_tagged() {
        let tagged = Edn::Tagged(
            Cow::Borrowed("inst"),
            Box::new(Edn::Str(Cow::Borrowed("2023-01-01"))),
        );
        assert_eq!(tagged.to_string(), r#"#inst "2023-01-01""#);
    }

    #[test]
    fn test_display_nested() {
        let inner = Edn::Vector(vec![Edn::Int(1), Edn::Int(2)].into());
        let outer = Edn::List(vec![Edn::Keyword(EdnKeyword::new("data")), inner].into());
        assert_eq!(outer.to_string(), "(:data [1 2])");
    }

    #[rstest]
    #[case("my tag")]
    #[case("")]
    #[case("123")]
    #[case("/bad")]
    fn test_display_tagged_invalid(#[case] tag: &str) {
        use std::fmt::Write;
        let tagged = Edn::Tagged(Cow::Borrowed(tag), Box::new(Edn::Nil));
        let mut buf = String::new();
        assert!(write!(buf, "{tagged}").is_err());
    }
}
