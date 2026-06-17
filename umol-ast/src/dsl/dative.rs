//! Dative-bond-string DSL.

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use umol_edn::{DeError, Edn, EdnError, EdnKeyword, EdnMap, EdnStreamDeserializer, FromEdn, ToEdn};
use winnow::ascii::multispace0;
use winnow::combinator::{preceded, repeat, terminated};
use winnow::error::ErrMode;
use winnow::token::take;
use winnow::Parser;

use super::config::DativeBondDefaults;
use super::error::{PResult, ParseError};
use super::predicates::{fmt_ring_count, optional_value, ring_count};
use super::value::{fmt_value, value, ValueDsl};
use crate::ast::constraint::DativeBondConstraint;
use crate::ast::dative::DativeBondAst;
use crate::ast::traits::{FromAst, IntoAst};
use crate::ast::value::ValueAst;

/// Surface DSL wrapper around `DativeBondAst`. The string form is the order
/// (number of donated electron pairs) followed by `#…` predicates,
/// paralleling `BondDsl`. Inline-capable constraints from
/// `DativeBondConstraint` are `RingCount` (`#R`) and `RingSize` (`#r`); the
/// remaining variants reference other entities and stay in the molecule
/// constraints container.
#[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DativeBondDsl(pub DativeBondAst);

impl DativeBondDsl {
    /// Zero-cost reference cast from `&DativeBondAst`. Relies on `repr(transparent)`.
    pub fn from_ref(ast: &DativeBondAst) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(ast as *const DativeBondAst as *const Self) }
    }
}

impl From<DativeBondAst> for DativeBondDsl {
    fn from(ast: DativeBondAst) -> Self {
        Self(ast)
    }
}

impl FromStr for DativeBondDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_dative_bond(s)
    }
}

impl Display for DativeBondDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_order(f, &self.0.order)?;
        for c in self.0.constraints.iter() {
            fmt_constraint(f, c)?;
        }
        Ok(())
    }
}

impl<'de> FromEdn<'de> for DativeBondDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s.parse().map_err(|e| DeError::subgrammar("dative", e)),
            Edn::Keyword(k) => {
                let s = expand_dative_keyword(k.name()).ok_or_else(|| {
                    DeError::Custom(format!("unknown dative keyword :{}", k.name()))
                })?;
                s.parse().map_err(|e| DeError::subgrammar("dative", e))
            }
            other => Err(DeError::TypeMismatch {
                expected: "string or dative-keyword",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        EdnStreamDeserializer::new(input).read_subgrammar_all("dative")
    }
}

/// Expand a dative-entry keyword shorthand to its equivalent dative-string
/// payload. Parallels [`super::bond::expand_bond_keyword`]:
///
/// - `:single` → `"1"`
/// - `:double` → `"2"`
/// - `:triple` → `"3"`
/// - `:quadruple` → `"4"`
///
/// Returns `None` for unrecognized keywords. Input sugar only — the AST
/// renders back to dative-string form.
pub(crate) fn expand_dative_keyword(name: &str) -> Option<&'static str> {
    match name {
        "single" => Some("1"),
        "double" => Some("2"),
        "triple" => Some("3"),
        "quadruple" => Some("4"),
        _ => None,
    }
}

impl ToEdn for DativeBondDsl {
    fn to_edn(&self) -> Edn<'static> {
        match dative_keyword_for(&self.0) {
            Some(kw) => Edn::Keyword(EdnKeyword::owned(kw.to_string())),
            None => Edn::Str(Cow::Owned(self.to_string())),
        }
    }
}

/// Return the dative-keyword shorthand for canonical dative shapes, or
/// `None` when the bond requires the full dative-string form. Inverse of
/// [`expand_dative_keyword`]: every shape this returns must round-trip.
///
/// Canonical means: no constraints and an integer order in 1..=4.
fn dative_keyword_for(ast: &DativeBondAst) -> Option<&'static str> {
    if !ast.constraints.is_empty() {
        return None;
    }
    match &ast.order {
        ValueAst::Lit(1) => Some("single"),
        ValueAst::Lit(2) => Some("double"),
        ValueAst::Lit(3) => Some("triple"),
        ValueAst::Lit(4) => Some("quadruple"),
        _ => None,
    }
}

impl FromAst<DativeBondAst> for DativeBondDsl {
    type Ctx = DativeBondDefaults;

    fn from_ast(ast: &DativeBondAst, _cfg: &Self::Ctx) -> Self {
        DativeBondDsl(ast.clone())
    }
}

impl IntoAst<DativeBondAst> for DativeBondDsl {
    type Ctx = DativeBondDefaults;

    fn into_ast(self, _cfg: &Self::Ctx) -> DativeBondAst {
        self.0
    }
}

impl FromStr for DativeBondAst {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(DativeBondDsl::from_str(s)?.into_ast(&DativeBondDefaults::default()))
    }
}

impl Display for DativeBondAst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        DativeBondDsl::from_ref(self).fmt(f)
    }
}

impl<'de> FromEdn<'de> for DativeBondAst {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        Ok(DativeBondDsl::from_edn(edn)?.into_ast(&DativeBondDefaults::default()))
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        Ok(DativeBondDsl::from_edn_str(input)?.into_ast(&DativeBondDefaults::default()))
    }
}

impl ToEdn for DativeBondAst {
    fn to_edn(&self) -> Edn<'static> {
        DativeBondDsl::from_ref(self).to_edn()
    }
}

/// Parse a complete dative-bond-string into a `DativeBondDsl`.
pub fn parse_dative_bond(input: &str) -> Result<DativeBondDsl, ParseError> {
    dative_bond.parse(input).map_err(|e| e.into_inner())
}

pub(crate) fn dative_bond(i: &mut &str) -> PResult<DativeBondDsl> {
    let order = preceded(multispace0, terminated(value, multispace0)).parse_next(i)?;
    let preds: Vec<DativeBondPredicate> =
        repeat(0.., terminated(dative_bond_predicate, multispace0)).parse_next(i)?;
    let mut form = DativeBondDsl(DativeBondAst::new(order));
    apply_predicates(&mut form, preds).map_err(ErrMode::Cut)?;
    Ok(form)
}

fn constraint_tag(c: &DativeBondConstraint) -> &'static str {
    match c {
        DativeBondConstraint::Aromatic => "#a",
        DativeBondConstraint::RingCount(_) => "#R",
        DativeBondConstraint::RingSize(_) => "#r",
    }
}

/// One predicate from a dative-bond-string; the parser yields a `Vec` of
/// these and the applier folds them into the `DativeBondAst`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DativeBondPredicate {
    Constraint(DativeBondConstraint),
}

fn dative_bond_predicate(i: &mut &str) -> PResult<DativeBondPredicate> {
    let start = *i;
    let prefix: &str = take(2usize).parse_next(i)?;
    match prefix {
        "#a" => Ok(DativeBondPredicate::Constraint(
            DativeBondConstraint::Aromatic,
        )),
        "#R" => ring_count
            .map(|v| DativeBondPredicate::Constraint(DativeBondConstraint::RingCount(v)))
            .parse_next(i),
        "#r" => optional_value
            .map(|v| DativeBondPredicate::Constraint(DativeBondConstraint::RingSize(v)))
            .parse_next(i),
        p if p.starts_with('#') => Err(ErrMode::Cut(ParseError::UnknownDativeBondPredicate(
            p.to_string(),
        ))),
        _ => Err(ErrMode::Cut(ParseError::TrailingInput(start.to_string()))),
    }
}

fn apply_predicates(
    form: &mut DativeBondDsl,
    preds: Vec<DativeBondPredicate>,
) -> Result<(), ParseError> {
    let ast = &mut form.0;
    for pred in preds {
        let DativeBondPredicate::Constraint(c) = pred;
        let tag = constraint_tag(&c);
        if ast
            .constraints
            .iter()
            .any(|existing| constraint_tag(existing) == tag)
        {
            return Err(ParseError::DuplicateDativeBondPredicate(tag.to_string()));
        }
        ast.constraints.add(c);
    }
    Ok(())
}

fn fmt_order(f: &mut fmt::Formatter<'_>, order: &ValueAst) -> fmt::Result {
    match order {
        ValueAst::Lit(n) => write!(f, "{}", n),
        ValueAst::Undetermined => write!(f, "*"),
        v => fmt_value(f, v),
    }
}

fn fmt_constraint(f: &mut fmt::Formatter<'_>, c: &DativeBondConstraint) -> fmt::Result {
    match c {
        DativeBondConstraint::Aromatic => write!(f, "#a"),
        DativeBondConstraint::RingCount(v) => fmt_ring_count(f, v),
        DativeBondConstraint::RingSize(v) => match v {
            ValueAst::Undetermined => Ok(()),
            ValueAst::Lit(1) => write!(f, "#r"),
            ValueAst::Lit(n) => write!(f, "#r{}", n),
            v => {
                write!(f, "#r")?;
                fmt_value(f, v)
            }
        },
    }
}

/// Surface DSL wrapper around the narrow `DativeBondConstraint`. EDN form is
/// either the bare keyword `:aromatic` (flag variant, no value) or a
/// single-key map keyed by the constraint kind: `{:ring-count <value>}` or
/// `{:ring-size <value>}`. Ref-bearing variants moved to
/// [`super::relational::RelationalConstraintDsl`].
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum DativeBondConstraintDsl {
    Aromatic,
    RingCount(ValueAst),
    RingSize(ValueAst),
}

impl<'de> FromEdn<'de> for DativeBondConstraintDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Keyword(k) if k.name() == "aromatic" => Ok(Self::Aromatic),
            Edn::Map(m) => {
                if m.len() != 1 {
                    return Err(DeError::Custom(format!(
                        "dative-bond-constraint must have exactly one key, got {}",
                        m.len()
                    )));
                }
                let (k, v) = m.iter().next().unwrap();
                let Edn::Keyword(key) = k else {
                    return Err(DeError::TypeMismatch {
                        expected: "keyword key",
                        got: k.kind(),
                        path: vec!["dative-bond-constraint".into()],
                    });
                };
                Ok(match key.name() {
                    "ring-count" => Self::RingCount(ValueDsl::from_edn(v)?.into_ast(&())),
                    "ring-size" => Self::RingSize(ValueDsl::from_edn(v)?.into_ast(&())),
                    other => {
                        return Err(DeError::UnknownField {
                            key: other.to_string(),
                            path: vec!["dative-bond-constraint".into()],
                        });
                    }
                })
            }
            other => Err(DeError::TypeMismatch {
                expected: ":aromatic / {:ring-count <value>} / {:ring-size <value>}",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for DativeBondConstraintDsl {
    fn to_edn(&self) -> Edn<'static> {
        match self {
            Self::Aromatic => Edn::Keyword(EdnKeyword::owned("aromatic".into())),
            Self::RingCount(v) | Self::RingSize(v) => {
                let key = match self {
                    Self::RingCount(_) => "ring-count",
                    Self::RingSize(_) => "ring-size",
                    Self::Aromatic => unreachable!(),
                };
                let mut m = EdnMap::with_capacity(1);
                m.insert(
                    Edn::Keyword(EdnKeyword::owned(key.into())),
                    ValueDsl::from_ast(v, &()).to_edn(),
                );
                Edn::Map(m)
            }
        }
    }
}

impl DativeBondConstraintDsl {
    /// Build from the narrow inline AST form.
    pub(crate) fn from_ast(c: &DativeBondConstraint) -> Self {
        match c {
            DativeBondConstraint::Aromatic => Self::Aromatic,
            DativeBondConstraint::RingCount(v) => Self::RingCount(v.clone()),
            DativeBondConstraint::RingSize(v) => Self::RingSize(v.clone()),
        }
    }

    /// Convert into the narrow inline AST form.
    pub(crate) fn into_ast(self) -> DativeBondConstraint {
        match self {
            Self::Aromatic => DativeBondConstraint::Aromatic,
            Self::RingCount(v) => DativeBondConstraint::RingCount(v),
            Self::RingSize(v) => DativeBondConstraint::RingSize(v),
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_edn::read_string;

    use super::*;
    use crate::ast::constraint::DativeBondConstraints;

    fn dative(order: ValueAst, constraints: DativeBondConstraints) -> DativeBondAst {
        DativeBondAst {
            acceptor_slot: 0,
            order,
            constraints,
        }
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single("1", DativeBondDsl(DativeBondAst::from_order(1)))]
    #[case::triple("3", DativeBondDsl(DativeBondAst::from_order(3)))]
    #[case::single_whitespace("  1  ", DativeBondDsl(DativeBondAst::from_order(1)))]
    #[case::undetermined_order("*", DativeBondDsl(DativeBondAst::default()))]
    #[case::aromatic("1#a", DativeBondDsl(dative(ValueAst::Lit(1), DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic]))))]
    #[case::aromatic_with_ring("1#a#r6", DativeBondDsl(dative(ValueAst::Lit(1), DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic, DativeBondConstraint::RingSize(ValueAst::Lit(6))]))))]
    #[case::ring_count("1#R2", DativeBondDsl(dative(ValueAst::Lit(1), DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Lit(2))]))))]
    #[case::ring_bare("1#R", DativeBondDsl(dative(ValueAst::Lit(1), DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Lit(1))]))))]
    #[case::ring_plus("1#R+", DativeBondDsl(dative(ValueAst::Lit(1), DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::var_at_least("r", 1))]))))]
    #[case::ring_undetermined("1#R*", DativeBondDsl(dative(ValueAst::Lit(1), DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Undetermined)]))))]
    #[case::ring_size("1#r6", DativeBondDsl(dative(ValueAst::Lit(1), DativeBondConstraints::from_iter([DativeBondConstraint::RingSize(ValueAst::Lit(6))]))))]
    #[case::ring_size_bare("1#r", DativeBondDsl(dative(ValueAst::Lit(1), DativeBondConstraints::from_iter([DativeBondConstraint::RingSize(ValueAst::Lit(1))]))))]
    #[case::ring_count_and_size("1#R2#r6", DativeBondDsl(dative(ValueAst::Lit(1), DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Lit(2)), DativeBondConstraint::RingSize(ValueAst::Lit(6))]))))]
    #[case::triple_with_constraint("3#R+", DativeBondDsl(dative(ValueAst::Lit(3), DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::var_at_least("r", 1))]))))]
    fn test_parse_dative(#[case] input: &str, #[case] expected: DativeBondDsl) {
        let result = dative_bond.parse(input);
        assert!(result.is_ok(), "{:?} should succeed, got {:?}", input, result.clone().unwrap_err());
        let form = result.unwrap();
        assert_eq!(form, expected);
    }

    /// Vacuous dative-bond constraints elide on rendering. `#R*` and `#r*`
    /// parse but the canonical form drops them.
    #[rstest]
    #[case::ring_count("1#R*", "1")]
    #[case::ring_size("1#r*", "1")]
    fn test_dative_render_elides_vacuous_constraints(
        #[case] input: &str,
        #[case] expected_canonical: &str,
    ) {
        let parsed: DativeBondDsl = dative_bond.parse(input).unwrap();
        assert_eq!(parsed.to_string(), expected_canonical);
        let reparsed: DativeBondDsl = dative_bond.parse(&parsed.to_string()).unwrap();
        assert!(
            reparsed.0.constraints.is_empty(),
            "vacuous constraint should be absent after render → reparse, got {:?}",
            reparsed.0.constraints,
        );
    }

    #[rstest]
    #[case::unknown("1#x", ParseError::UnknownDativeBondPredicate("#x".to_string()))]
    #[case::unknown_c("1#c", ParseError::UnknownDativeBondPredicate("#c".to_string()))]
    #[case::dup_ring("1#R1#R2", ParseError::DuplicateDativeBondPredicate("#R".to_string()))]
    #[case::dup_ring_size("1#r6#r5", ParseError::DuplicateDativeBondPredicate("#r".to_string()))]
    #[case::dup_aromatic("1#a#a", ParseError::DuplicateDativeBondPredicate("#a".to_string()))]
    #[case::trailing("1#R2 foo", ParseError::TrailingInput("foo".to_string()))]
    fn test_parse_dative_error(#[case] input: &str, #[case] expected: ParseError) {
        let result = dative_bond.parse(input);
        assert!(result.is_err(), "{:?} should fail", input);
        let err = result.unwrap_err().into_inner();
        assert_eq!(err, expected);
    }

    #[rstest]
    #[case::single("1")]
    #[case::triple("3")]
    #[case::undetermined("*")]
    #[case::ring_count("1#R2")]
    #[case::ring_size("1#r6")]
    #[case::both("1#R2#r6")]
    #[case::aromatic("1#a")]
    #[case::aromatic_with_ring("1#a#r6")]
    fn test_dative_roundtrip(#[case] input: &str) {
        let form: DativeBondDsl = input.parse().unwrap();
        let rendered = form.to_string();
        let reparsed: DativeBondDsl = rendered.parse().unwrap();
        assert_eq!(form, reparsed);
    }

    #[rstest]
    fn test_dative_dsl_to_ast_passthrough() {
        let dsl = DativeBondDsl(dative(
            ValueAst::Lit(1),
            DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Lit(2))]),
        ));
        let cfg = DativeBondDefaults::zeroed();
        let ast = dsl.into_ast(&cfg);
        assert_eq!(
            ast.constraints,
            DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Lit(2))])
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(":single", DativeBondAst::from_order(1))]
    #[case::double(":double", DativeBondAst::from_order(2))]
    #[case::triple(":triple", DativeBondAst::from_order(3))]
    #[case::quadruple(":quadruple", DativeBondAst::from_order(4))]
    fn test_dative_dsl_keyword_shorthand(#[case] input: &str, #[case] expected: DativeBondAst) {
        let edn = read_string(input).unwrap();
        let dsl = DativeBondDsl::from_edn(&edn).unwrap();
        assert_eq!(dsl.0, expected);
    }

    #[rstest]
    fn test_dative_dsl_keyword_shorthand_unknown_rejected() {
        let edn = read_string(":bogus").unwrap();
        let err = DativeBondDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::Custom(_)));
    }

    #[rstest]
    #[case::single(r##""1""##)]
    #[case::ring_count(r##""1#R2""##)]
    #[case::ring_count_and_size(r##""1#R2#r6""##)]
    fn test_dative_dsl_from_edn_str_matches_from_edn(#[case] input: &str) {
        let via_stream = DativeBondDsl::from_edn_str(input).unwrap();
        let tree = read_string(input).unwrap();
        let via_tree = DativeBondDsl::from_edn(&tree).unwrap();
        assert_eq!(via_stream, via_tree);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraint::Aromatic, ":aromatic")]
    #[case::ring_count(DativeBondConstraint::RingCount(ValueAst::Lit(2)), "{:ring-count 2}")]
    #[case::ring_size(DativeBondConstraint::RingSize(ValueAst::Lit(6)), "{:ring-size 6}")]
    fn test_dative_bond_constraint_dsl_roundtrip(
        #[case] input: DativeBondConstraint,
        #[case] edn_source: &str,
    ) {
        let dsl = DativeBondConstraintDsl::from_ast(&input);
        let edn = dsl.clone().to_edn();
        let expected = read_string(edn_source).unwrap();
        assert_eq!(edn, expected, "render mismatch");
        let parsed = DativeBondConstraintDsl::from_edn(&edn).unwrap();
        let back = parsed.into_ast();
        assert_eq!(back, input, "parse-back mismatch");
    }

    #[rstest]
    fn test_dative_bond_constraint_dsl_rejects_wrong_shape() {
        let err = DativeBondConstraintDsl::from_edn(&Edn::Int(3)).unwrap_err();
        assert!(matches!(err, DeError::TypeMismatch { .. }));
    }

    #[rstest]
    fn test_dative_bond_constraint_dsl_rejects_unknown_key() {
        let edn = read_string("{:bogus 1}").unwrap();
        let err = DativeBondConstraintDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::UnknownField { .. }));
    }

    #[rstest]
    #[case::single("1")]
    #[case::triple("3")]
    fn test_dative_bond_ast_from_str_to_string_roundtrip(#[case] s: &str) {
        let ast: DativeBondAst = s.parse().unwrap();
        assert_eq!(ast.to_string(), s);
    }
}
