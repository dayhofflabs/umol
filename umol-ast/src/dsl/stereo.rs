//! Stereo config-string DSL (Phase D1): the `class config` surface — a class head
//! (`Th`/`Ct`/`Sp`/`Tb`/`Oh`) plus the `* | ! | + | <coset-term>` config —
//! ↔ `(StereoKind, StereoConfigurationAst)`. One config parser/writer, reused by
//! the `#T`/`#C` constraint strings (D5) and the element `:type` head (D3).
//!
//! Whitespace is ignored. `~`/`^` build *symbolic* `SwapOp`/`ApplyOp` nodes; the
//! class involution and the `^` action are resolved later by
//! `StereoConfigurationAst::simplify(kind)`, never at parse. Unary `~` binds
//! tighter than postfix `^` (`~3^2134` = `(~3)^2134`), left-associative.

// The config parser and class head are exercised by the tests but not yet
// reached from the library build; D5 wires them into the atom/bond
// constraint-string parsers. (The writer is already used by `fmt_constraint`.)
#![allow(dead_code)]

use std::fmt::{self, Write};

use winnow::ascii::{digit1, multispace0};
use winnow::combinator::{alt, opt, preceded, separated, terminated};
use winnow::error::ErrMode;
use winnow::Parser;

use umol_perm::Permutation;

use super::error::{PResult, ParseError};
use super::value::{id, terminator};
use crate::ast::stereo::{StereoConfigurationAst, StereoExpr, StereoIndexAst, StereoKind};

/// `Th` / `Ct` / `Sp` / `Tb` / `Oh` head → `StereoKind`.
pub(crate) fn class(i: &mut &str) -> PResult<StereoKind> {
    alt((
        "Th".value(StereoKind::Tetrahedral),
        "Ct".value(StereoKind::CisTrans),
        "Sp".value(StereoKind::SquarePlanar),
        "Tb".value(StereoKind::TrigonalBipyramidal),
        "Oh".value(StereoKind::Octahedral),
    ))
    .parse_next(i)
}

/// Write the class head for `kind`.
pub(crate) fn fmt_class(f: &mut fmt::Formatter<'_>, kind: StereoKind) -> fmt::Result {
    match kind {
        StereoKind::Tetrahedral => f.write_str("Th"),
        StereoKind::CisTrans => f.write_str("Ct"),
        StereoKind::SquarePlanar => f.write_str("Sp"),
        StereoKind::TrigonalBipyramidal => f.write_str("Tb"),
        StereoKind::Octahedral => f.write_str("Oh"),
    }
}

/// Parse the `config` grammar into a `StereoConfigurationAst`.
pub(crate) fn stereo_config(i: &mut &str) -> PResult<StereoConfigurationAst> {
    preceded(
        multispace0,
        alt((
            '*'.value(StereoConfigurationAst::Undetermined),
            '!'.value(StereoConfigurationAst::NotStereo),
            '+'.value(StereoConfigurationAst::Stereo(StereoIndexAst::Undetermined)),
            terminated(nat, (multispace0, terminator))
                .map(|n| StereoConfigurationAst::Stereo(StereoIndexAst::Lit(n))),
            coset_term.map(|e| StereoConfigurationAst::Stereo(StereoIndexAst::Expr(Box::new(e)))),
        )),
    )
    .parse_next(i)
}

/// Write the `config` for `c` (no class head — the caller writes `#T`/`#C` or
/// the `Th`/`Ct` head).
pub(crate) fn fmt_stereo_config(
    f: &mut fmt::Formatter<'_>,
    c: &StereoConfigurationAst,
) -> fmt::Result {
    match c {
        StereoConfigurationAst::Undetermined => f.write_char('*'),
        StereoConfigurationAst::NotStereo => f.write_char('!'),
        StereoConfigurationAst::Stereo(StereoIndexAst::Undetermined) => f.write_char('+'),
        StereoConfigurationAst::Stereo(StereoIndexAst::Lit(n)) => write!(f, "{n}"),
        StereoConfigurationAst::Stereo(StereoIndexAst::Expr(e)) => fmt_expr(f, e),
    }
}

/// `coset-term`: a `~`-prefixed base carrying zero or more left-associative
/// `^image` postfixes.
fn coset_term(i: &mut &str) -> PResult<StereoExpr> {
    let mut e = prefix_term(i)?;
    loop {
        multispace0.parse_next(i)?;
        if opt('^').parse_next(i)?.is_some() {
            e = StereoExpr::apply(e, image(i)?);
        } else {
            return Ok(e);
        }
    }
}

/// `'~' prefix-term | base` — unary `~` binds tighter than `^`.
fn prefix_term(i: &mut &str) -> PResult<StereoExpr> {
    multispace0.parse_next(i)?;
    if opt('~').parse_next(i)?.is_some() {
        Ok(StereoExpr::swap(prefix_term(i)?))
    } else {
        base(i)
    }
}

/// `nat | '?' id ('::' set)? | set`.
fn base(i: &mut &str) -> PResult<StereoExpr> {
    preceded(
        multispace0,
        alt((
            nat_set.map(StereoExpr::LitSet),
            var_or_domain,
            nat.map(StereoExpr::Lit),
        )),
    )
    .parse_next(i)
}

fn var_or_domain(i: &mut &str) -> PResult<StereoExpr> {
    let name = preceded('?', id).parse_next(i)?;
    let domain = opt(preceded((multispace0, "::", multispace0), nat_set)).parse_next(i)?;
    Ok(match domain {
        Some(set) => StereoExpr::VarDomain(name, set),
        None => StereoExpr::Var(name),
    })
}

fn nat(i: &mut &str) -> PResult<u32> {
    let s: &str = digit1.parse_next(i)?;
    s.parse::<u32>()
        .map_err(|_| ErrMode::Backtrack(ParseError::Syntax))
}

/// `'{' nat (',' nat)* '}'`, whitespace-insensitive.
fn nat_set(i: &mut &str) -> PResult<Vec<u32>> {
    preceded(
        ('{', multispace0),
        terminated(
            separated(1.., nat, (multispace0, ',', multispace0)),
            (multispace0, '}'),
        ),
    )
    .parse_next(i)
}

/// The 1-indexed one-line image read as digits → `Permutation`. Kept literal
/// (never canonicalized); validated as a bijection of `1..=degree`.
fn image(i: &mut &str) -> PResult<Permutation> {
    let digits: &str = digit1.parse_next(i)?;
    let img: Vec<u8> = digits
        .bytes()
        .map(|b| b.checked_sub(b'1'))
        .collect::<Option<Vec<u8>>>()
        .ok_or(ErrMode::Cut(ParseError::Syntax))?;
    let degree = img.len();
    let mut seen = [false; 6];
    for &x in &img {
        let x = x as usize;
        if x >= degree || x >= 6 || seen[x] {
            return Err(ErrMode::Cut(ParseError::Syntax));
        }
        seen[x] = true;
    }
    Ok(Permutation::from_image(degree, &img))
}

fn fmt_expr(f: &mut fmt::Formatter<'_>, e: &StereoExpr) -> fmt::Result {
    match e {
        StereoExpr::Lit(n) => write!(f, "{n}"),
        StereoExpr::Var(name) => write!(f, "?{name}"),
        StereoExpr::SwapOp(inner) => {
            f.write_char('~')?;
            fmt_expr(f, inner)
        }
        StereoExpr::ApplyOp(inner, perm) => {
            fmt_expr(f, inner)?;
            f.write_char('^')?;
            fmt_image(f, *perm)
        }
        StereoExpr::LitSet(set) => fmt_nat_set(f, set),
        StereoExpr::VarDomain(name, set) => {
            write!(f, "?{name} :: ")?;
            fmt_nat_set(f, set)
        }
    }
}

fn fmt_image(f: &mut fmt::Formatter<'_>, perm: Permutation) -> fmt::Result {
    for i in 0..perm.degree() {
        write!(f, "{}", perm.apply(i) + 1)?;
    }
    Ok(())
}

fn fmt_nat_set(f: &mut fmt::Formatter<'_>, set: &[u32]) -> fmt::Result {
    f.write_char('{')?;
    for (i, n) in set.iter().enumerate() {
        if i > 0 {
            f.write_char(',')?;
        }
        write!(f, "{n}")?;
    }
    f.write_char('}')
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined("*", StereoConfigurationAst::Undetermined)]
    #[case::not_stereo("!", StereoConfigurationAst::NotStereo)]
    #[case::stereogenic("+", StereoConfigurationAst::Stereo(StereoIndexAst::Undetermined))]
    #[case::lit("1", StereoConfigurationAst::from(1_u32))]
    #[case::var("?o", StereoConfigurationAst::Stereo(StereoIndexAst::expr(StereoExpr::Var("o".to_string()))))]
    #[case::lit_set("{1,2}", StereoConfigurationAst::Stereo(StereoIndexAst::expr(StereoExpr::LitSet(vec![1, 2]))))]
    #[case::var_domain("?o :: {1,2}", StereoConfigurationAst::Stereo(StereoIndexAst::expr(StereoExpr::VarDomain("o".to_string(), vec![1, 2]))))]
    #[case::swap("~1", StereoConfigurationAst::Stereo(StereoIndexAst::expr(StereoExpr::swap(StereoExpr::Lit(1)))))]
    #[case::apply("1^2134", StereoConfigurationAst::Stereo(StereoIndexAst::expr(StereoExpr::apply(StereoExpr::Lit(1), Permutation::from_image(4, &[1, 0, 2, 3])))))]
    #[case::swap_binds_tighter_than_apply("~1^2134", StereoConfigurationAst::Stereo(StereoIndexAst::expr(StereoExpr::apply(StereoExpr::swap(StereoExpr::Lit(1)), Permutation::from_image(4, &[1, 0, 2, 3])))))]
    #[case::whitespace_ignored("  ?o :: { 1 , 2 }", StereoConfigurationAst::Stereo(StereoIndexAst::expr(StereoExpr::VarDomain("o".to_string(), vec![1, 2]))))]
    fn test_stereo_config(#[case] input: &str, #[case] expected: StereoConfigurationAst) {
        assert_eq!(stereo_config.parse(input).unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(StereoConfigurationAst::Undetermined, "*")]
    #[case::not_stereo(StereoConfigurationAst::NotStereo, "!")]
    #[case::stereogenic(StereoConfigurationAst::Stereo(StereoIndexAst::Undetermined), "+")]
    #[case::lit(StereoConfigurationAst::from(1_u32), "1")]
    #[case::var(StereoConfigurationAst::Stereo(StereoIndexAst::expr(StereoExpr::Var("o".to_string()))), "?o")]
    #[case::lit_set(StereoConfigurationAst::Stereo(StereoIndexAst::expr(StereoExpr::LitSet(vec![1, 2]))), "{1,2}")]
    #[case::var_domain(StereoConfigurationAst::Stereo(StereoIndexAst::expr(StereoExpr::VarDomain("o".to_string(), vec![1, 2]))), "?o :: {1,2}")]
    #[case::swap(StereoConfigurationAst::Stereo(StereoIndexAst::expr(StereoExpr::swap(StereoExpr::Lit(1)))), "~1")]
    #[case::apply(StereoConfigurationAst::Stereo(StereoIndexAst::expr(StereoExpr::apply(StereoExpr::Lit(1), Permutation::from_image(4, &[1, 0, 2, 3])))), "1^2134")]
    fn test_fmt_stereo_config(#[case] c: StereoConfigurationAst, #[case] expected: &str) {
        struct W<'a>(&'a StereoConfigurationAst);
        impl fmt::Display for W<'_> {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                fmt_stereo_config(f, self.0)
            }
        }
        assert_eq!(W(&c).to_string(), expected);
    }

    #[rstest]
    #[case::tetrahedral("Th", StereoKind::Tetrahedral)]
    #[case::cis_trans("Ct", StereoKind::CisTrans)]
    #[case::square_planar("Sp", StereoKind::SquarePlanar)]
    #[case::trigonal_bipyramidal("Tb", StereoKind::TrigonalBipyramidal)]
    #[case::octahedral("Oh", StereoKind::Octahedral)]
    fn test_class(#[case] input: &str, #[case] expected: StereoKind) {
        assert_eq!(class.parse(input).unwrap(), expected);
    }
}
