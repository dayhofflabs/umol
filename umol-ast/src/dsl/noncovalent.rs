//! Noncovalent-bond-string DSL.

use std::borrow::Cow;
use std::convert::Infallible;
use std::fmt::{self, Display};
use std::str::FromStr;

use umol_edn::{DeError, Edn, EdnError, EdnStreamDeserializer, FromEdn, ToEdn};
use winnow::ascii::multispace0;
use winnow::combinator::{alt, delimited, preceded, separated, terminated};
use winnow::error::{ErrMode, ParserError};
use winnow::token::one_of;
use winnow::Parser;

use super::atom::AtomConstraintDsl;
use super::constraint::{AtomRef, EntityCounts};
use super::molecule::Metadata;
use super::error::{PResult, ParseError};
use super::value::id;
use crate::ast::constraint::NoncovalentBondConstraint;
use crate::ast::noncovalent::{NoncovalentBondAst, NoncovalentKind, NoncovalentKindAst};
use crate::ast::traits::{FromAst, IntoAst};
use crate::dsl::config::NoncovalentBondDefaults;

/// Surface DSL wrapper around `NoncovalentBondAst`. String form is the
/// noncovalent-kind expression (three-letter literal, set, bind, ref, or `*`).
/// All `NoncovalentBondConstraint` variants are molecule-scope.
#[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct NoncovalentBondDsl(pub NoncovalentBondAst);

impl NoncovalentBondDsl {
    /// Zero-cost reference cast from `&NoncovalentBondAst`. Relies on `repr(transparent)`.
    pub fn from_ref(ast: &NoncovalentBondAst) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(ast as *const NoncovalentBondAst as *const Self) }
    }
}

impl FromStr for NoncovalentBondDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_noncovalent(s)
    }
}

impl Display for NoncovalentBondDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_noncovalent_ast(f, &self.0)
    }
}

impl<'de> FromEdn<'de> for NoncovalentBondDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s.parse().map_err(|e| DeError::subgrammar("noncovalent", e)),
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        EdnStreamDeserializer::new(input).read_subgrammar_all("noncovalent")
    }
}

impl ToEdn for NoncovalentBondDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

impl FromAst<NoncovalentBondAst> for NoncovalentBondDsl {
    type Ctx = NoncovalentBondDefaults;
    type Error = ParseError;

    fn from_ast(
        ast: &NoncovalentBondAst,
        _cfg: &Self::Ctx,
    ) -> Result<Self, ParseError> {
        Ok(NoncovalentBondDsl(ast.clone()))
    }
}

impl IntoAst<NoncovalentBondAst> for NoncovalentBondDsl {
    type Ctx = NoncovalentBondDefaults;
    type Error = ParseError;

    fn into_ast(self, _cfg: &Self::Ctx) -> Result<NoncovalentBondAst, ParseError> {
        Ok(self.0)
    }
}

// -- Parse --------------------

pub fn parse_noncovalent(input: &str) -> Result<NoncovalentBondDsl, ParseError> {
    noncovalent.parse(input).map_err(|e| e.into_inner())
}

pub(crate) fn noncovalent(i: &mut &str) -> PResult<NoncovalentBondDsl> {
    let kind = delimited(multispace0, kind_expr, multispace0).parse_next(i)?;
    Ok(NoncovalentBondDsl(NoncovalentBondAst::new(kind)))
}

fn kind_expr(i: &mut &str) -> PResult<NoncovalentKindAst> {
    alt((
        '*'.value(NoncovalentKindAst::Undetermined),
        kind_set.map(NoncovalentKindAst::Set),
        kind_bind.map(|(id, set)| NoncovalentKindAst::Bind { id, set }),
        kind_ref.map(NoncovalentKindAst::Ref),
        kind_literal.map(NoncovalentKindAst::Lit),
    ))
    .parse_next(i)
    .map_err(|_: ErrMode<ParseError>| ErrMode::Backtrack(ParseError::ExpectedNoncovalentKind))
}

fn kind_literal(i: &mut &str) -> PResult<NoncovalentKind> {
    let sym: &str = (
        one_of(|c: char| c.is_ascii_uppercase()),
        one_of(|c: char| c.is_ascii_lowercase()),
        one_of(|c: char| c.is_ascii_lowercase()),
    )
        .take()
        .parse_next(i)?;
    match kind_from_symbol(sym) {
        Some(k) => Ok(k),
        None => Err(ErrMode::Backtrack(ParseError::from_input(i))),
    }
}

fn kind_set(i: &mut &str) -> PResult<Vec<NoncovalentKind>> {
    delimited(
        '{',
        delimited(
            multispace0,
            separated(1.., kind_literal, delimited(multispace0, ',', multispace0)),
            multispace0,
        ),
        '}',
    )
    .parse_next(i)
}

fn kind_bind(i: &mut &str) -> PResult<(String, Vec<NoncovalentKind>)> {
    delimited(
        '(',
        (
            delimited(multispace0, preceded('?', id), multispace0),
            preceded(("::", multispace0), terminated(kind_set, multispace0)),
        ),
        ')',
    )
    .parse_next(i)
}

fn kind_ref(i: &mut &str) -> PResult<String> {
    delimited(
        '(',
        delimited(multispace0, preceded('?', id), multispace0),
        ')',
    )
    .parse_next(i)
}

fn kind_from_symbol(sym: &str) -> Option<NoncovalentKind> {
    match sym {
        "Hbd" => Some(NoncovalentKind::HydrogenBond),
        "Xbd" => Some(NoncovalentKind::HalogenBond),
        "Ybd" => Some(NoncovalentKind::ChalcogenBond),
        "Ion" => Some(NoncovalentKind::Ionic),
        "Vdw" => Some(NoncovalentKind::VanDerWaals),
        _ => None,
    }
}

fn kind_symbol(k: NoncovalentKind) -> &'static str {
    match k {
        NoncovalentKind::HydrogenBond => "Hbd",
        NoncovalentKind::HalogenBond => "Xbd",
        NoncovalentKind::ChalcogenBond => "Ybd",
        NoncovalentKind::Ionic => "Ion",
        NoncovalentKind::VanDerWaals => "Vdw",
    }
}

// -- Format --------------------

fn fmt_noncovalent_ast(f: &mut fmt::Formatter<'_>, ast: &NoncovalentBondAst) -> fmt::Result {
    fmt_kind(f, &ast.kind)
}

fn fmt_kind(f: &mut fmt::Formatter<'_>, kind: &NoncovalentKindAst) -> fmt::Result {
    match kind {
        NoncovalentKindAst::Lit(k) => write!(f, "{}", kind_symbol(*k)),
        NoncovalentKindAst::Undetermined => write!(f, "*"),
        NoncovalentKindAst::Set(ks) => {
            write!(f, "{{")?;
            for (i, k) in ks.iter().enumerate() {
                if i > 0 {
                    write!(f, ",")?;
                }
                write!(f, "{}", kind_symbol(*k))?;
            }
            write!(f, "}}")
        }
        NoncovalentKindAst::Bind { id, set } => {
            write!(f, "(?{} :: {{", id)?;
            for (i, k) in set.iter().enumerate() {
                if i > 0 {
                    write!(f, ",")?;
                }
                write!(f, "{}", kind_symbol(*k))?;
            }
            write!(f, "}})")
        }
        NoncovalentKindAst::Ref(id) => write!(f, "(?{})", id),
    }
}

// -- Constraint DSL -------------------

/// Surface DSL wrapper around `NoncovalentBondConstraint`. Mirrors the AST
/// enum with atom refs in place of `AtomIdx`. EDN form is a single-key map:
/// `{:ends [a b]}`, `{:contains ref}`, or `{:ends-satisfy [<atom-constraint>
/// <atom-constraint>]}`.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum NoncovalentBondConstraintDsl {
    Ends([AtomRef; 2]),
    Contains(AtomRef),
    EndsSatisfy([Box<AtomConstraintDsl>; 2]),
}

impl NoncovalentBondConstraintDsl {
    pub(crate) fn from_ast(
        c: &NoncovalentBondConstraint,
        meta: &Metadata,
    ) -> Result<Self, Infallible> {
        Ok(match c {
            NoncovalentBondConstraint::Ends([a, b]) => Self::Ends([
                AtomRef::from_ast(*a, meta),
                AtomRef::from_ast(*b, meta),
            ]),
            NoncovalentBondConstraint::Contains(a) => {
                Self::Contains(AtomRef::from_ast(*a, meta))
            }
            NoncovalentBondConstraint::EndsSatisfy([a, b]) => Self::EndsSatisfy([
                Box::new(AtomConstraintDsl::from_ast(a, &()).unwrap()),
                Box::new(AtomConstraintDsl::from_ast(b, &()).unwrap()),
            ]),
        })
    }

    pub(crate) fn into_ast(
        self,
        counts: &EntityCounts,
        meta: &Metadata,
    ) -> Result<NoncovalentBondConstraint, ParseError> {
        Ok(match self {
            Self::Ends([a, b]) => NoncovalentBondConstraint::Ends([
                a.into_ast(counts.atom_count, meta)?,
                b.into_ast(counts.atom_count, meta)?,
            ]),
            Self::Contains(r) => {
                NoncovalentBondConstraint::Contains(r.into_ast(counts.atom_count, meta)?)
            }
            Self::EndsSatisfy([a, b]) => NoncovalentBondConstraint::EndsSatisfy([
                Box::new(a.into_ast(&()).unwrap()),
                Box::new(b.into_ast(&()).unwrap()),
            ]),
        })
    }
}

impl<'de> FromEdn<'de> for NoncovalentBondConstraintDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let Edn::Map(m) = edn else {
            return Err(DeError::TypeMismatch {
                expected: "noncovalent-bond-constraint single-key map",
                got: edn.kind(),
                path: Vec::new(),
            });
        };
        if m.len() != 1 {
            return Err(DeError::Custom(format!(
                "noncovalent-bond-constraint must have exactly one key, got {}",
                m.len()
            )));
        }
        let (k, v) = m.iter().next().unwrap();
        let Edn::Keyword(key) = k else {
            return Err(DeError::TypeMismatch {
                expected: "keyword key",
                got: k.kind(),
                path: vec!["noncovalent-bond-constraint".into()],
            });
        };
        Ok(match key.name() {
            "ends" => Self::Ends(parse_pair::<AtomRef>(v, "ends")?),
            "contains" => Self::Contains(AtomRef::from_edn(v)?),
            "ends-satisfy" => Self::EndsSatisfy(parse_pair_boxed::<AtomConstraintDsl>(
                v,
                "ends-satisfy",
            )?),
            other => {
                return Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["noncovalent-bond-constraint".into()],
                });
            }
        })
    }
}

impl ToEdn for NoncovalentBondConstraintDsl {
    fn to_edn(&self) -> Edn<'static> {
        let (key, value) = match self {
            Self::Ends([a, b]) => (
                "ends",
                Edn::Vector(vec![a.to_edn(), b.to_edn()].into()),
            ),
            Self::Contains(r) => ("contains", r.to_edn()),
            Self::EndsSatisfy([a, b]) => (
                "ends-satisfy",
                Edn::Vector(vec![a.to_edn(), b.to_edn()].into()),
            ),
        };
        let mut m = umol_edn::EdnMap::with_capacity(1);
        m.insert(Edn::Keyword(umol_edn::EdnKeyword::owned(key.into())), value);
        Edn::Map(m)
    }
}

fn parse_pair<T>(edn: &Edn<'_>, context: &'static str) -> Result<[T; 2], DeError>
where
    T: for<'de> FromEdn<'de>,
{
    let Edn::Vector(v) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "2-element vector",
            got: edn.kind(),
            path: vec![context.into()],
        });
    };
    if v.len() != 2 {
        return Err(DeError::Custom(format!(
            "{}: expected 2 elements, got {}",
            context,
            v.len()
        )));
    }
    Ok([T::from_edn(&v[0])?, T::from_edn(&v[1])?])
}

fn parse_pair_boxed<T>(edn: &Edn<'_>, context: &'static str) -> Result<[Box<T>; 2], DeError>
where
    T: for<'de> FromEdn<'de>,
{
    let [a, b] = parse_pair::<T>(edn, context)?;
    Ok([Box::new(a), Box::new(b)])
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::super::molecule::Metadata;
    use super::*;
    use crate::ast::constraint::{AtomConstraint, NoncovalentBondConstraint};
    use crate::ast::idx::AtomIdx;
    use crate::ast::value::ValueAst;

    #[rustfmt::skip]
    #[rstest]
    #[case::hbond("Hbd", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentKind::HydrogenBond)))]
    #[case::xbond("Xbd", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentKind::HalogenBond)))]
    #[case::ybond("Ybd", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentKind::ChalcogenBond)))]
    #[case::ion("Ion", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentKind::Ionic)))]
    #[case::vdw("Vdw", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentKind::VanDerWaals)))]
    #[case::whitespace("  Hbd  ", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentKind::HydrogenBond)))]
    #[case::undetermined("*", NoncovalentBondDsl(NoncovalentBondAst::new(NoncovalentKindAst::Undetermined)))]
    #[case::set("{Hbd,Ion}", NoncovalentBondDsl(NoncovalentBondAst::new(NoncovalentKindAst::Set(vec![NoncovalentKind::HydrogenBond, NoncovalentKind::Ionic]))))]
    #[case::set_spaced("{ Hbd, Vdw }", NoncovalentBondDsl(NoncovalentBondAst::new(NoncovalentKindAst::Set(vec![NoncovalentKind::HydrogenBond, NoncovalentKind::VanDerWaals]))))]
    #[case::bind("(?k :: {Hbd,Ion})", NoncovalentBondDsl(NoncovalentBondAst::new(NoncovalentKindAst::Bind { id: "k".to_string(), set: vec![NoncovalentKind::HydrogenBond, NoncovalentKind::Ionic] })))]
    #[case::ref_("(?k)", NoncovalentBondDsl(NoncovalentBondAst::new(NoncovalentKindAst::Ref("k".to_string()))))]
    fn test_parse_noncovalent(#[case] input: &str, #[case] expected: NoncovalentBondDsl) {
        let result = noncovalent.parse(input);
        assert!(result.is_ok(), "{:?} should succeed, got {:?}", input, result.clone().unwrap_err());
        let form = result.unwrap();
        assert_eq!(form, expected);
    }

    #[rstest]
    #[case::empty("")]
    #[case::unknown_literal("Abc")]
    #[case::two_letter("Hb")]
    #[case::bare_paren("(")]
    fn test_parse_noncovalent_invalid(#[case] input: &str) {
        let result = noncovalent.parse(input);
        assert!(result.is_err(), "{:?} should fail", input);
    }

    #[rstest]
    #[case::hbond("Hbd")]
    #[case::ion("Ion")]
    #[case::undetermined("*")]
    #[case::set("{Hbd,Ion}")]
    #[case::bind("(?k :: {Hbd,Ion})")]
    #[case::ref_("(?k)")]
    fn test_noncovalent_roundtrip(#[case] input: &str) {
        let form: NoncovalentBondDsl = input.parse().unwrap();
        let rendered = form.to_string();
        let reparsed: NoncovalentBondDsl = rendered.parse().unwrap();
        assert_eq!(form, reparsed);
    }

    #[rstest]
    fn test_noncovalent_dsl_to_ast_passthrough() {
        let dsl = NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentKind::HydrogenBond));
        let cfg = NoncovalentBondDefaults::zeroed();
        let ast = dsl.into_ast(&cfg).unwrap();
        assert_eq!(
            ast.kind,
            NoncovalentKindAst::Lit(NoncovalentKind::HydrogenBond)
        );
    }

    #[rstest]
    #[case::single(r##""Hbd""##)]
    #[case::set(r##""{Hbd,Ion}""##)]
    #[case::undetermined(r##""*""##)]
    fn test_noncovalent_dsl_from_edn_str_matches_from_edn(#[case] input: &str) {
        let via_stream = NoncovalentBondDsl::from_edn_str(input).unwrap();
        let tree = umol_edn::read_string(input).unwrap();
        let via_tree = NoncovalentBondDsl::from_edn(&tree).unwrap();
        assert_eq!(via_stream, via_tree);
    }

    // -- NoncovalentBondConstraintDsl ----------------

    fn counts_with_atoms(atom_count: usize) -> EntityCounts {
        EntityCounts {
            atom_count,
            bond_count: 0,
            dative_bond_count: 0,
            aromatic_system_count: 0,
            multicenter_bond_count: 0,
            noncovalent_bond_count: 0,
        }
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ends(NoncovalentBondConstraint::Ends([AtomIdx(0), AtomIdx(1)]), "{:ends [0 1]}")]
    #[case::contains(NoncovalentBondConstraint::Contains(AtomIdx(3)), "{:contains 3}")]
    #[case::ends_satisfy(NoncovalentBondConstraint::EndsSatisfy([
        Box::new(AtomConstraint::Valence(ValueAst::Lit(4))),
        Box::new(AtomConstraint::Degree(ValueAst::Lit(3))),
    ]), "{:ends-satisfy [{:valence 4} {:degree 3}]}")]
    fn test_noncovalent_bond_constraint_dsl_roundtrip(
        #[case] input: NoncovalentBondConstraint,
        #[case] edn_source: &str,
    ) {
        let meta = Metadata::default();
        let dsl = NoncovalentBondConstraintDsl::from_ast(&input, &meta).unwrap();
        let edn = dsl.clone().to_edn();
        let expected = umol_edn::read_string(edn_source).unwrap();
        assert_eq!(edn, expected, "render mismatch");
        let parsed = NoncovalentBondConstraintDsl::from_edn(&edn).unwrap();
        let back = parsed.into_ast(&counts_with_atoms(10), &meta).unwrap();
        assert_eq!(back, input, "parse-back mismatch");
    }

    #[rstest]
    fn test_noncovalent_bond_constraint_dsl_rejects_out_of_range_atom() {
        let meta = Metadata::default();
        let edn = umol_edn::read_string("{:contains 99}").unwrap();
        let dsl = NoncovalentBondConstraintDsl::from_edn(&edn).unwrap();
        let err = dsl.into_ast(&counts_with_atoms(5), &meta).unwrap_err();
        assert_eq!(
            err,
            ParseError::InvalidRef {
                kind: "atom",
                value: "99".into()
            }
        );
    }

    #[rstest]
    fn test_noncovalent_bond_constraint_dsl_rejects_unknown_key() {
        let edn = umol_edn::read_string("{:bogus 1}").unwrap();
        let err = NoncovalentBondConstraintDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::UnknownField { .. }));
    }

    #[rstest]
    fn test_noncovalent_bond_constraint_dsl_rejects_wrong_pair_length() {
        let edn = umol_edn::read_string("{:ends [0]}").unwrap();
        let err = NoncovalentBondConstraintDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::Custom(_)));
    }
}
