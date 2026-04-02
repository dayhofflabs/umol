//! Display implementation for Edn values (compact EDN format).

use std::fmt;

use crate::edn::Edn;

impl fmt::Display for Edn<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Edn::Nil => write!(f, "nil"),
            Edn::Bool(true) => write!(f, "true"),
            Edn::Bool(false) => write!(f, "false"),
            Edn::Int(n) => write!(f, "{n}"),
            Edn::Float(v) => format_float(f, *v),
            Edn::Char(c) => format_char(f, *c),
            Edn::Str(s) => format_string(f, s),
            Edn::Keyword(k) => write!(f, "{k}"),
            Edn::Symbol(s) => write!(f, "{s}"),
            Edn::List(items) => format_seq(f, "(", ")", items),
            Edn::Vector(items) => format_seq(f, "[", "]", items),
            Edn::Map(m) => {
                write!(f, "{{")?;
                for (i, (k, v)) in m.iter().enumerate() {
                    if i > 0 {
                        write!(f, " ")?;
                    }
                    write!(f, "{k} {v}")?;
                }
                write!(f, "}}")
            }
            Edn::Set(s) => {
                write!(f, "#{{")?;
                for (i, v) in s.iter().enumerate() {
                    if i > 0 {
                        write!(f, " ")?;
                    }
                    write!(f, "{v}")?;
                }
                write!(f, "}}")
            }
            Edn::Tagged(tag, inner) => write!(f, "#{tag} {inner}"),
        }
    }
}

fn format_float(f: &mut fmt::Formatter<'_>, v: f64) -> fmt::Result {
    // EDN has no syntax for NaN/Inf — these values cannot be serialized.
    assert!(!v.is_nan() && v.is_finite(), "EDN cannot represent NaN or Infinity");
    let s = format!("{v}");
    if s.contains('.') || s.contains('e') || s.contains('E') {
        write!(f, "{s}")
    } else {
        write!(f, "{s}.0")
    }
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

fn format_seq(f: &mut fmt::Formatter<'_>, open: &str, close: &str, items: &[Edn<'_>]) -> fmt::Result {
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
    use crate::edn::{Edn, Keyword, Symbol};

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

    #[test]
    #[should_panic(expected = "EDN cannot represent NaN or Infinity")]
    fn test_display_nan_panics() {
        let _ = Edn::Float(f64::NAN).to_string();
    }

    #[test]
    #[should_panic(expected = "EDN cannot represent NaN or Infinity")]
    fn test_display_inf_panics() {
        let _ = Edn::Float(f64::INFINITY).to_string();
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
    #[case(Keyword::new("foo"), ":foo")]
    #[case(Keyword::namespaced("ns", "bar"), ":ns/bar")]
    fn test_display_keyword(#[case] k: Keyword<'_>, #[case] expected: &str) {
        assert_eq!(Edn::Keyword(k).to_string(), expected);
    }

    #[rstest]
    #[case(Symbol::new("foo"), "foo")]
    #[case(Symbol::namespaced("ns", "bar"), "ns/bar")]
    fn test_display_symbol(#[case] s: Symbol<'_>, #[case] expected: &str) {
        assert_eq!(Edn::Symbol(s).to_string(), expected);
    }

    #[test]
    fn test_display_list() {
        let v = Edn::List(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)].into());
        assert_eq!(v.to_string(), "(1 2 3)");
    }

    #[test]
    fn test_display_list_empty() {
        assert_eq!(Edn::List(vec![].into()).to_string(), "()");
    }

    #[test]
    fn test_display_vector() {
        let v = Edn::Vector(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)].into());
        assert_eq!(v.to_string(), "[1 2 3]");
    }

    #[test]
    fn test_display_vector_empty() {
        assert_eq!(Edn::Vector(vec![].into()).to_string(), "[]");
    }

    #[test]
    fn test_display_map() {
        let mut m = EdnMap::new();
        m.insert(Edn::Keyword(Keyword::new("a")), Edn::Int(1));
        assert_eq!(Edn::Map(m).to_string(), "{:a 1}");
    }

    #[test]
    fn test_display_map_multi() {
        let mut m = EdnMap::new();
        m.insert(Edn::Keyword(Keyword::new("a")), Edn::Int(1));
        m.insert(Edn::Keyword(Keyword::new("b")), Edn::Int(2));
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
            "inst".to_string(),
            Box::new(Edn::Str(Cow::Borrowed("2023-01-01"))),
        );
        assert_eq!(tagged.to_string(), r#"#inst "2023-01-01""#);
    }

    #[test]
    fn test_display_nested() {
        let inner = Edn::Vector(vec![Edn::Int(1), Edn::Int(2)].into());
        let outer = Edn::List(vec![Edn::Keyword(Keyword::new("data")), inner].into());
        assert_eq!(outer.to_string(), "(:data [1 2])");
    }
}
