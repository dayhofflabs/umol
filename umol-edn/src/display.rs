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
                        write!(f, ", ")?;
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
    if v.is_nan() {
        write!(f, "##NaN")
    } else if v == f64::INFINITY {
        write!(f, "##Inf")
    } else if v == f64::NEG_INFINITY {
        write!(f, "##-Inf")
    } else {
        let s = format!("{v}");
        // Ensure there's always a decimal point
        if s.contains('.') || s.contains('e') || s.contains('E') {
            write!(f, "{s}")
        } else {
            write!(f, "{s}.0")
        }
    }
}

fn format_char(f: &mut fmt::Formatter<'_>, c: char) -> fmt::Result {
    match c {
        '\n' => write!(f, "\\newline"),
        '\r' => write!(f, "\\return"),
        ' ' => write!(f, "\\space"),
        '\t' => write!(f, "\\tab"),
        '\u{000C}' => write!(f, "\\formfeed"),
        '\u{0008}' => write!(f, "\\backspace"),
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
            '\u{0008}' => write!(f, "\\b")?,
            '\u{000C}' => write!(f, "\\f")?,
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
