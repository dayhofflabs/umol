//! Predicate sub-parsers and apply helpers shared between atom and bond DSL.

use std::fmt;

use winnow::ascii::multispace0;
use winnow::combinator::{alt, empty, preceded};
use winnow::error::ErrMode;
use winnow::Parser;

use super::error::{PResult, ParseError};
use super::value::value;
use crate::ast::spin::SpinStateAst;
use crate::ast::value::{Expr, RelOp, ValueAst};

pub(crate) fn charge(i: &mut &str) -> PResult<ValueAst> {
    preceded(
        multispace0,
        alt((
            value,
            "+".value(ValueAst::Lit(1)),
            "-".value(ValueAst::Lit(-1)),
        )),
    )
    .parse_next(i)
    .map_err(|_: ErrMode<ParseError>| ErrMode::Backtrack(ParseError::ExpectedPredicateBody))
}

pub(crate) fn optional_value(i: &mut &str) -> PResult<ValueAst> {
    preceded(multispace0, alt((value, empty.value(ValueAst::Lit(1)))))
        .parse_next(i)
        .map_err(|_: ErrMode<ParseError>| ErrMode::Backtrack(ParseError::ExpectedPredicateBody))
}

pub(crate) fn ring_count(i: &mut &str) -> PResult<ValueAst> {
    preceded(
        multispace0,
        alt((
            value,
            "+".value(ValueAst::Expr(Expr::Rel(
                Box::new(Expr::Var("r".to_string())),
                RelOp::Ge,
                Box::new(Expr::Lit(1)),
            ))),
            empty.value(ValueAst::Lit(1)),
        )),
    )
    .parse_next(i)
    .map_err(|_: ErrMode<ParseError>| ErrMode::Backtrack(ParseError::ExpectedPredicateBody))
}

/// Recognize `Expr(Rel(Var(name), Ge, Lit(threshold)))` used as the `+` sugar
/// for aromatic valence (`a >= 0`), multicenter valence (`m >= 0`), and
/// ring count (`r >= 1`).
pub(crate) fn is_plus_sugar(v: &ValueAst, name: &str, threshold: i64) -> bool {
    match v {
        ValueAst::Expr(Expr::Rel(l, RelOp::Ge, r)) => {
            matches!(l.as_ref(), Expr::Var(n) if n == name)
                && matches!(r.as_ref(), Expr::Lit(n) if *n == threshold)
        }
        _ => false,
    }
}

pub(crate) fn fmt_value(f: &mut fmt::Formatter<'_>, v: &ValueAst) -> fmt::Result {
    match v {
        ValueAst::Undetermined => write!(f, "*"),
        ValueAst::Lit(n) => write!(f, "{}", n),
        ValueAst::LitSet(s) => {
            write!(f, "{{")?;
            for (i, n) in s.iter().enumerate() {
                if i > 0 {
                    write!(f, ",")?;
                }
                write!(f, "{}", n)?;
            }
            write!(f, "}}")
        }
        ValueAst::Expr(_) => write!(f, "<expr>"),
    }
}

pub(crate) fn fmt_charge(f: &mut fmt::Formatter<'_>, v: &ValueAst) -> fmt::Result {
    match v {
        ValueAst::Undetermined => Ok(()),
        ValueAst::Lit(0) => write!(f, "#c0"),
        ValueAst::Lit(1) => write!(f, "#c+"),
        ValueAst::Lit(-1) => write!(f, "#c-"),
        ValueAst::Lit(n) if *n > 0 => write!(f, "#c+{}", n),
        ValueAst::Lit(n) => write!(f, "#c{}", n),
        v => {
            write!(f, "#c")?;
            fmt_value(f, v)
        }
    }
}

pub(crate) fn fmt_spin_pair(f: &mut fmt::Formatter<'_>, spin: &SpinStateAst) -> fmt::Result {
    match &spin.unpaired {
        ValueAst::Undetermined => {}
        ValueAst::Lit(1) => write!(f, "#u")?,
        ValueAst::Lit(n) => write!(f, "#u{}", n)?,
        v => {
            write!(f, "#u")?;
            fmt_value(f, v)?;
        }
    }
    match &spin.multiplicity {
        ValueAst::Undetermined => {}
        ValueAst::Lit(1) => write!(f, "#s")?,
        ValueAst::Lit(n) => write!(f, "#s{}", n)?,
        v => {
            write!(f, "#s")?;
            fmt_value(f, v)?;
        }
    }
    Ok(())
}

pub(crate) fn fmt_ring_count(f: &mut fmt::Formatter<'_>, v: &ValueAst) -> fmt::Result {
    if is_plus_sugar(v, "r", 1) {
        return write!(f, "#R+");
    }
    match v {
        ValueAst::Undetermined => write!(f, "#R*"),
        ValueAst::Lit(1) => write!(f, "#R"),
        ValueAst::Lit(n) => write!(f, "#R{}", n),
        v => {
            write!(f, "#R")?;
            fmt_value(f, v)
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SpinPredicate {
    Unpaired(ValueAst),
    Multiplicity(ValueAst),
}

pub(crate) fn apply_spin_pair(
    spin: &mut SpinStateAst,
    pred: SpinPredicate,
    dup: fn(String) -> ParseError,
) -> Result<(), ParseError> {
    match pred {
        SpinPredicate::Unpaired(v) => {
            if !matches!(&spin.unpaired, ValueAst::Undetermined) {
                return Err(dup("#u".to_string()));
            }
            spin.unpaired = v;
        }
        SpinPredicate::Multiplicity(v) => {
            if !matches!(&spin.multiplicity, ValueAst::Undetermined) {
                return Err(dup("#s".to_string()));
            }
            spin.multiplicity = v;
        }
    }
    spin.validate().map_err(ParseError::from_spin_state_error)?;
    Ok(())
}
