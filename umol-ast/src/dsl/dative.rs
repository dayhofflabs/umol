//! Dative-bond-string DSL.

use std::borrow::Cow;
use std::convert::Infallible;
use std::fmt::{self, Display};
use std::str::FromStr;

use umol_edn::{DeError, Edn, EdnError, EdnStreamDeserializer, FromEdn, ToEdn};
use winnow::ascii::multispace0;
use winnow::combinator::{repeat, terminated};
use winnow::error::ErrMode;
use winnow::token::take;
use winnow::Parser;

use super::atom::AtomConstraintDsl;
use super::constraint::{AtomRef, BondRef, EntityCounts};
use super::molecule::Metadata;
use super::error::{PResult, ParseError};
use super::predicates::{fmt_ring_count, optional_value, ring_count};
use super::value::{fmt_value, ValueDsl};
use crate::ast::constraint::DativeBondConstraint;
use crate::ast::dative::DativeBondAst;
use crate::ast::traits::{FromAst, IntoAst};
use crate::ast::value::ValueAst;
use crate::dsl::config::DativeBondDefaults;

/// Surface DSL wrapper around `DativeBondAst`. No leading token; the string
/// form is a sequence of `#…` predicates. Inline-capable constraints from
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

impl FromStr for DativeBondDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_dative(s)
    }
}

impl Display for DativeBondDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
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
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        EdnStreamDeserializer::new(input).read_subgrammar_all("dative")
    }
}

impl ToEdn for DativeBondDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

impl FromAst<DativeBondAst> for DativeBondDsl {
    type Ctx = DativeBondDefaults;
    type Error = ParseError;

    fn from_ast(ast: &DativeBondAst, _cfg: &Self::Ctx) -> Result<Self, ParseError> {
        Ok(DativeBondDsl(ast.clone()))
    }
}

impl IntoAst<DativeBondAst> for DativeBondDsl {
    type Ctx = DativeBondDefaults;
    type Error = ParseError;

    fn into_ast(self, _cfg: &Self::Ctx) -> Result<DativeBondAst, ParseError> {
        Ok(self.0)
    }
}

// -- Parse --------------------

pub fn parse_dative(input: &str) -> Result<DativeBondDsl, ParseError> {
    dative.parse(input).map_err(|e| e.into_inner())
}

pub(crate) fn dative(i: &mut &str) -> PResult<DativeBondDsl> {
    multispace0.parse_next(i)?;
    let preds: Vec<DativePredicate> =
        repeat(0.., terminated(dative_predicate, multispace0)).parse_next(i)?;
    let mut form = DativeBondDsl::default();
    apply_predicates(&mut form, preds).map_err(ErrMode::Cut)?;
    Ok(form)
}

fn inline_constraint_tag(c: &DativeBondConstraint) -> Option<&'static str> {
    match c {
        DativeBondConstraint::RingCount(_) => Some("#R"),
        DativeBondConstraint::RingSize(_) => Some("#r"),
        _ => None,
    }
}

fn constraint_tag(c: &DativeBondConstraint) -> &'static str {
    inline_constraint_tag(c).expect("non-inline-capable dative constraint produced by parser")
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DativePredicate {
    Constraint(DativeBondConstraint),
}

fn dative_predicate(i: &mut &str) -> PResult<DativePredicate> {
    let start = *i;
    let prefix: &str = take(2usize).parse_next(i)?;
    match prefix {
        "#R" => ring_count
            .map(|v| DativePredicate::Constraint(DativeBondConstraint::RingCount(v)))
            .parse_next(i),
        "#r" => optional_value
            .map(|v| DativePredicate::Constraint(DativeBondConstraint::RingSize(v)))
            .parse_next(i),
        p if p.starts_with('#') => Err(ErrMode::Cut(ParseError::UnknownDativePredicate(
            p.to_string(),
        ))),
        _ => Err(ErrMode::Cut(ParseError::TrailingInput(start.to_string()))),
    }
}

fn apply_predicates(
    form: &mut DativeBondDsl,
    preds: Vec<DativePredicate>,
) -> Result<(), ParseError> {
    let ast = &mut form.0;
    for pred in preds {
        let DativePredicate::Constraint(c) = pred;
        let tag = constraint_tag(&c);
        if ast
            .constraints
            .iter()
            .any(|existing| inline_constraint_tag(existing) == Some(tag))
        {
            return Err(ParseError::DuplicateDativePredicate(tag.to_string()));
        }
        ast.constraints.add(c);
    }
    Ok(())
}

// -- Format --------------------

fn fmt_constraint(f: &mut fmt::Formatter<'_>, c: &DativeBondConstraint) -> fmt::Result {
    match c {
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
        _ => Ok(()),
    }
}

// -- Constraint DSL -------------------

/// Surface DSL wrapper around `DativeBondConstraint`. Mirrors the AST enum
/// with atom/bond refs in place of `AtomIdx` / `BondIdx`. EDN form is a
/// single-key map keyed by the constraint kind (e.g. `{:donor :n1}`,
/// `{:donor-satisfies {:valence 4}}`, `{:parallels 0}`).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum DativeBondConstraintDsl {
    RingCount(ValueAst),
    RingSize(ValueAst),
    Donor(AtomRef),
    Acceptor(AtomRef),
    DonorSatisfies(Box<AtomConstraintDsl>),
    AcceptorSatisfies(Box<AtomConstraintDsl>),
    Parallels(BondRef),
}

impl<'de> FromEdn<'de> for DativeBondConstraintDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let Edn::Map(m) = edn else {
            return Err(DeError::TypeMismatch {
                expected: "dative-bond-constraint single-key map",
                got: edn.kind(),
                path: Vec::new(),
            });
        };
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
            "ring-count" => Self::RingCount(ValueDsl::from_edn(v)?.into_ast(&()).unwrap()),
            "ring-size" => Self::RingSize(ValueDsl::from_edn(v)?.into_ast(&()).unwrap()),
            "donor" => Self::Donor(AtomRef::from_edn(v)?),
            "acceptor" => Self::Acceptor(AtomRef::from_edn(v)?),
            "donor-satisfies" => Self::DonorSatisfies(Box::new(AtomConstraintDsl::from_edn(v)?)),
            "acceptor-satisfies" => {
                Self::AcceptorSatisfies(Box::new(AtomConstraintDsl::from_edn(v)?))
            }
            "parallels" => Self::Parallels(BondRef::from_edn(v)?),
            other => {
                return Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["dative-bond-constraint".into()],
                });
            }
        })
    }
}

impl ToEdn for DativeBondConstraintDsl {
    fn to_edn(&self) -> Edn<'static> {
        let (key, value) = match self {
            Self::RingCount(v) => ("ring-count", ValueDsl::from_ast(v, &()).unwrap().to_edn()),
            Self::RingSize(v) => ("ring-size", ValueDsl::from_ast(v, &()).unwrap().to_edn()),
            Self::Donor(r) => ("donor", r.to_edn()),
            Self::Acceptor(r) => ("acceptor", r.to_edn()),
            Self::DonorSatisfies(c) => ("donor-satisfies", c.to_edn()),
            Self::AcceptorSatisfies(c) => ("acceptor-satisfies", c.to_edn()),
            Self::Parallels(r) => ("parallels", r.to_edn()),
        };
        let mut m = umol_edn::EdnMap::with_capacity(1);
        m.insert(Edn::Keyword(umol_edn::EdnKeyword::owned(key.into())), value);
        Edn::Map(m)
    }
}

impl DativeBondConstraintDsl {
    pub(crate) fn from_ast(
        c: &DativeBondConstraint,
        meta: &Metadata,
    ) -> Result<Self, Infallible> {
        Ok(match c {
            DativeBondConstraint::RingCount(v) => Self::RingCount(v.clone()),
            DativeBondConstraint::RingSize(v) => Self::RingSize(v.clone()),
            DativeBondConstraint::Donor(a) => Self::Donor(AtomRef::from_ast(*a, meta)),
            DativeBondConstraint::Acceptor(a) => Self::Acceptor(AtomRef::from_ast(*a, meta)),
            DativeBondConstraint::DonorSatisfies(c) => {
                Self::DonorSatisfies(Box::new(AtomConstraintDsl::from_ast(c, &()).unwrap()))
            }
            DativeBondConstraint::AcceptorSatisfies(c) => {
                Self::AcceptorSatisfies(Box::new(AtomConstraintDsl::from_ast(c, &()).unwrap()))
            }
            DativeBondConstraint::Parallels(b) => Self::Parallels(BondRef::from_ast(*b, meta)),
        })
    }

    pub(crate) fn into_ast(
        self,
        counts: &EntityCounts,
        meta: &Metadata,
    ) -> Result<DativeBondConstraint, ParseError> {
        Ok(match self {
            Self::RingCount(v) => DativeBondConstraint::RingCount(v),
            Self::RingSize(v) => DativeBondConstraint::RingSize(v),
            Self::Donor(r) => DativeBondConstraint::Donor(r.into_ast(counts.atom_count, meta)?),
            Self::Acceptor(r) => {
                DativeBondConstraint::Acceptor(r.into_ast(counts.atom_count, meta)?)
            }
            Self::DonorSatisfies(c) => {
                DativeBondConstraint::DonorSatisfies(Box::new(c.into_ast(&()).unwrap()))
            }
            Self::AcceptorSatisfies(c) => {
                DativeBondConstraint::AcceptorSatisfies(Box::new(c.into_ast(&()).unwrap()))
            }
            Self::Parallels(r) => {
                DativeBondConstraint::Parallels(r.into_ast(counts.bond_count, meta)?)
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::super::molecule::{Metadata, MetadataBuilder};
    use super::*;
    use crate::ast::constraint::{AtomConstraint, DativeBondConstraints};
    use crate::ast::dative::DativeDirection;
    use crate::ast::idx::{AtomIdx, BondIdx};
    use crate::ast::value::{Expr, RelOp};

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("", DativeBondDsl(DativeBondAst::default()))]
    #[case::whitespace("   ", DativeBondDsl(DativeBondAst::default()))]
    #[case::ring_count("#R2", DativeBondDsl(DativeBondAst { direction: DativeDirection::Forward, constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Lit(2))]) }))]
    #[case::ring_bare("#R", DativeBondDsl(DativeBondAst { direction: DativeDirection::Forward, constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Lit(1))]) }))]
    #[case::ring_plus("#R+", DativeBondDsl(DativeBondAst { direction: DativeDirection::Forward, constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("r".to_string())), RelOp::Ge, Box::new(Expr::Lit(1)))))]) }))]
    #[case::ring_undetermined("#R*", DativeBondDsl(DativeBondAst { direction: DativeDirection::Forward, constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Undetermined)]) }))]
    #[case::ring_size("#r6", DativeBondDsl(DativeBondAst { direction: DativeDirection::Forward, constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingSize(ValueAst::Lit(6))]) }))]
    #[case::ring_size_bare("#r", DativeBondDsl(DativeBondAst { direction: DativeDirection::Forward, constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingSize(ValueAst::Lit(1))]) }))]
    #[case::ring_count_and_size("#R2#r6", DativeBondDsl(DativeBondAst { direction: DativeDirection::Forward, constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Lit(2)), DativeBondConstraint::RingSize(ValueAst::Lit(6))]) }))]
    fn test_parse_dative(#[case] input: &str, #[case] expected: DativeBondDsl) {
        let result = dative.parse(input);
        assert!(result.is_ok(), "{:?} should succeed, got {:?}", input, result.clone().unwrap_err());
        let form = result.unwrap();
        assert_eq!(form, expected);
    }

    #[rstest]
    #[case::unknown("#x", ParseError::UnknownDativePredicate("#x".to_string()))]
    #[case::unknown_c("#c", ParseError::UnknownDativePredicate("#c".to_string()))]
    #[case::dup_ring("#R1#R2", ParseError::DuplicateDativePredicate("#R".to_string()))]
    #[case::dup_ring_size("#r6#r5", ParseError::DuplicateDativePredicate("#r".to_string()))]
    #[case::trailing("#R2 foo", ParseError::TrailingInput("foo".to_string()))]
    fn test_parse_dative_error(#[case] input: &str, #[case] expected: ParseError) {
        let result = dative.parse(input);
        assert!(result.is_err(), "{:?} should fail", input);
        let err = result.unwrap_err().into_inner();
        assert_eq!(err, expected);
    }

    #[rstest]
    #[case::empty("")]
    #[case::ring_count("#R2")]
    #[case::ring_size("#r6")]
    #[case::both("#R2#r6")]
    fn test_dative_roundtrip(#[case] input: &str) {
        let form: DativeBondDsl = input.parse().unwrap();
        let rendered = form.to_string();
        let reparsed: DativeBondDsl = rendered.parse().unwrap();
        assert_eq!(form, reparsed);
    }

    #[rstest]
    fn test_dative_dsl_to_ast_passthrough() {
        let dsl = DativeBondDsl(DativeBondAst {
            direction: DativeDirection::Forward,
            constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(
                ValueAst::Lit(2),
            )]),
        });
        let cfg = DativeBondDefaults::zeroed();
        let ast = dsl.into_ast(&cfg).unwrap();
        assert_eq!(
            ast.constraints,
            DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Lit(2))])
        );
    }

    #[rstest]
    #[case::empty(r##""""##)]
    #[case::ring_count(r##""#R2""##)]
    #[case::ring_count_and_size(r##""#R2#r6""##)]
    fn test_dative_dsl_from_edn_str_matches_from_edn(#[case] input: &str) {
        let via_stream = DativeBondDsl::from_edn_str(input).unwrap();
        let tree = umol_edn::read_string(input).unwrap();
        let via_tree = DativeBondDsl::from_edn(&tree).unwrap();
        assert_eq!(via_stream, via_tree);
    }

    // -- DativeBondConstraintDsl ----------------

    fn metadata_with_atom_and_bond_id() -> Metadata {
        let mut b = MetadataBuilder::default();
        b.set_atom_id(AtomIdx(1), "n1".to_string());
        b.set_bond_id(BondIdx(2), "b1".to_string());
        b.build()
    }

    fn counts_with(atom_count: usize, bond_count: usize) -> EntityCounts {
        EntityCounts {
            atom_count,
            bond_count,
            dative_bond_count: 0,
            aromatic_system_count: 0,
            multicenter_bond_count: 0,
            noncovalent_bond_count: 0,
        }
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ring_count(DativeBondConstraint::RingCount(ValueAst::Lit(2)), "{:ring-count 2}")]
    #[case::ring_size(DativeBondConstraint::RingSize(ValueAst::Lit(6)), "{:ring-size 6}")]
    #[case::donor_idx(DativeBondConstraint::Donor(AtomIdx(3)), "{:donor 3}")]
    #[case::acceptor_idx(DativeBondConstraint::Acceptor(AtomIdx(4)), "{:acceptor 4}")]
    #[case::parallels_idx(DativeBondConstraint::Parallels(BondIdx(5)), "{:parallels 5}")]
    #[case::donor_satisfies(DativeBondConstraint::DonorSatisfies(Box::new(AtomConstraint::Valence(ValueAst::Lit(3)))), "{:donor-satisfies {:valence 3}}")]
    #[case::acceptor_satisfies(DativeBondConstraint::AcceptorSatisfies(Box::new(AtomConstraint::Degree(ValueAst::Lit(2)))), "{:acceptor-satisfies {:degree 2}}")]
    fn test_dative_bond_constraint_dsl_roundtrip_indices(
        #[case] input: DativeBondConstraint,
        #[case] edn_source: &str,
    ) {
        let meta = Metadata::default();
        let dsl = DativeBondConstraintDsl::from_ast(&input, &meta).unwrap();
        let edn = dsl.clone().to_edn();
        let expected = umol_edn::read_string(edn_source).unwrap();
        assert_eq!(edn, expected, "render mismatch");
        let parsed = DativeBondConstraintDsl::from_edn(&edn).unwrap();
        let back = parsed.into_ast(&counts_with(10, 10), &meta).unwrap();
        assert_eq!(back, input, "parse-back mismatch");
    }

    #[rstest]
    fn test_dative_bond_constraint_dsl_uses_keyword_for_known_atom() {
        let meta = metadata_with_atom_and_bond_id();
        let dsl = DativeBondConstraintDsl::from_ast(
            &DativeBondConstraint::Donor(AtomIdx(1)),
            &meta,
        )
        .unwrap();
        let edn = dsl.to_edn();
        assert_eq!(edn, umol_edn::read_string("{:donor :n1}").unwrap());
    }

    #[rstest]
    fn test_dative_bond_constraint_dsl_resolves_keyword_id() {
        let meta = metadata_with_atom_and_bond_id();
        let edn = umol_edn::read_string("{:donor :n1}").unwrap();
        let dsl = DativeBondConstraintDsl::from_edn(&edn).unwrap();
        let back = dsl.into_ast(&counts_with(10, 10), &meta).unwrap();
        assert_eq!(back, DativeBondConstraint::Donor(AtomIdx(1)));
    }

    #[rstest]
    fn test_dative_bond_constraint_dsl_resolves_bond_keyword() {
        let meta = metadata_with_atom_and_bond_id();
        let edn = umol_edn::read_string("{:parallels :b1}").unwrap();
        let dsl = DativeBondConstraintDsl::from_edn(&edn).unwrap();
        let back = dsl.into_ast(&counts_with(10, 10), &meta).unwrap();
        assert_eq!(back, DativeBondConstraint::Parallels(BondIdx(2)));
    }

    #[rstest]
    fn test_dative_bond_constraint_dsl_rejects_out_of_range_index() {
        let meta = Metadata::default();
        let edn = umol_edn::read_string("{:donor 99}").unwrap();
        let dsl = DativeBondConstraintDsl::from_edn(&edn).unwrap();
        let err = dsl.into_ast(&counts_with(5, 5), &meta).unwrap_err();
        assert_eq!(
            err,
            ParseError::InvalidRef {
                kind: "atom",
                value: "99".into()
            }
        );
    }

    #[rstest]
    fn test_dative_bond_constraint_dsl_rejects_unknown_id() {
        let meta = Metadata::default();
        let edn = umol_edn::read_string("{:acceptor :nope}").unwrap();
        let dsl = DativeBondConstraintDsl::from_edn(&edn).unwrap();
        let err = dsl.into_ast(&counts_with(5, 5), &meta).unwrap_err();
        assert_eq!(
            err,
            ParseError::InvalidRef {
                kind: "atom",
                value: "nope".into()
            }
        );
    }

    #[rstest]
    fn test_dative_bond_constraint_dsl_rejects_wrong_shape() {
        let err = DativeBondConstraintDsl::from_edn(&Edn::Int(3)).unwrap_err();
        assert!(matches!(err, DeError::TypeMismatch { .. }));
    }

    #[rstest]
    fn test_dative_bond_constraint_dsl_rejects_unknown_key() {
        let edn = umol_edn::read_string("{:bogus 1}").unwrap();
        let err = DativeBondConstraintDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::UnknownField { .. }));
    }
}
