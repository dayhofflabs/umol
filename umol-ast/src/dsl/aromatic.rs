//! Aromatic-system-string DSL.

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
use super::constraint::{AtomRef, ResolveContext};
use super::error::{PResult, ParseError};
use super::predicates::{
    apply_spin_pair, charge, fmt_charge, fmt_spin_pair, lower_spin, optional_value,
    raise_spin, SpinPredicate,
};
use super::value::fmt_value;
use crate::ast::aromatic::AromaticSystemAst;
use crate::dsl::config::{AromaticSystemDefaults, NumericDefault};
use crate::ast::constraint::AromaticSystemConstraint;
use crate::ast::traits::{FromAst, IntoAst};
use crate::ast::value::ValueAst;

/// Surface DSL wrapper around `AromaticSystemAst`. Parses and renders the
/// aromatic-system-string form. All `AromaticSystemConstraint` variants are
/// molecule-scope, so nothing from the constraint vec serializes inline.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AromaticSystemDsl(pub AromaticSystemAst);

impl FromStr for AromaticSystemDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_aromatic(s)
    }
}

impl Display for AromaticSystemDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_aromatic_ast(f, &self.0)
    }
}

impl<'de> FromEdn<'de> for AromaticSystemDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s.parse().map_err(|e| DeError::subgrammar("aromatic", e)),
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        EdnStreamDeserializer::new(input).read_subgrammar_all("aromatic")
    }
}

impl ToEdn for AromaticSystemDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

impl FromAst<AromaticSystemAst> for AromaticSystemDsl {
    type Ctx<'a> = AromaticSystemDefaults;
    type Error = ParseError;

    fn from_ast<'a>(
        ast: &AromaticSystemAst,
        cfg: &Self::Ctx<'a>,
    ) -> Result<Self, ParseError> {
        let mut out = ast.clone();
        lower_aromatic(&mut out, cfg);
        Ok(AromaticSystemDsl(out))
    }
}

impl IntoAst<AromaticSystemAst> for AromaticSystemDsl {
    type Ctx<'a> = AromaticSystemDefaults;
    type Error = ParseError;

    fn into_ast<'a>(
        mut self,
        cfg: &Self::Ctx<'a>,
    ) -> Result<AromaticSystemAst, ParseError> {
        raise_aromatic(&mut self.0, cfg);
        Ok(self.0)
    }
}

// -- Parse --------------------

pub fn parse_aromatic(input: &str) -> Result<AromaticSystemDsl, ParseError> {
    aromatic.parse(input).map_err(|e| e.into_inner())
}

pub(crate) fn aromatic(i: &mut &str) -> PResult<AromaticSystemDsl> {
    multispace0.parse_next(i)?;
    let preds: Vec<AromaticPredicate> =
        repeat(0.., terminated(aromatic_predicate, multispace0)).parse_next(i)?;
    let mut form = AromaticSystemDsl::default();
    apply_predicates(&mut form, preds).map_err(ErrMode::Cut)?;
    Ok(form)
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AromaticPredicate {
    Charge(ValueAst),
    Spin(SpinPredicate),
    Electrons(ValueAst),
}

fn aromatic_predicate(i: &mut &str) -> PResult<AromaticPredicate> {
    let start = *i;
    let prefix: &str = take(2usize).parse_next(i)?;
    match prefix {
        "#c" => charge.map(AromaticPredicate::Charge).parse_next(i),
        "#u" => optional_value
            .map(|v| AromaticPredicate::Spin(SpinPredicate::Unpaired(v)))
            .parse_next(i),
        "#s" => optional_value
            .map(|v| AromaticPredicate::Spin(SpinPredicate::Multiplicity(v)))
            .parse_next(i),
        "#e" => optional_value
            .map(AromaticPredicate::Electrons)
            .parse_next(i),
        p if p.starts_with('#') => Err(ErrMode::Cut(ParseError::UnknownAromaticPredicate(
            p.to_string(),
        ))),
        _ => Err(ErrMode::Cut(ParseError::TrailingInput(start.to_string()))),
    }
}

fn apply_predicates(
    form: &mut AromaticSystemDsl,
    preds: Vec<AromaticPredicate>,
) -> Result<(), ParseError> {
    let ast = &mut form.0;
    for pred in preds {
        match pred {
            AromaticPredicate::Charge(v) => {
                if !matches!(ast.charge, ValueAst::Undetermined) {
                    return Err(ParseError::DuplicateAromaticPredicate("#c".to_string()));
                }
                ast.charge = v;
            }
            AromaticPredicate::Spin(sp) => {
                apply_spin_pair(&mut ast.spin, sp, ParseError::DuplicateAromaticPredicate)?;
            }
            AromaticPredicate::Electrons(v) => {
                if !matches!(ast.electrons, ValueAst::Undetermined) {
                    return Err(ParseError::DuplicateAromaticPredicate("#e".to_string()));
                }
                ast.electrons = v;
            }
        }
    }
    Ok(())
}

// -- Format --------------------

fn fmt_aromatic_ast(f: &mut fmt::Formatter<'_>, ast: &AromaticSystemAst) -> fmt::Result {
    fmt_charge(f, &ast.charge)?;
    fmt_spin_pair(f, &ast.spin)?;
    fmt_electrons(f, &ast.electrons)
}

fn fmt_electrons(f: &mut fmt::Formatter<'_>, v: &ValueAst) -> fmt::Result {
    match v {
        ValueAst::Undetermined => Ok(()),
        ValueAst::Lit(1) => write!(f, "#e"),
        ValueAst::Lit(n) => write!(f, "#e{}", n),
        v => {
            write!(f, "#e")?;
            fmt_value(f, v)
        }
    }
}

// -- Raise --------------------

fn raise_aromatic(ast: &mut AromaticSystemAst, cfg: &AromaticSystemDefaults) {
    // Exhaustive destructure: adding a new AromaticSystemAst field is a
    // compile error here, forcing the author to decide how raising should
    // handle it.
    let AromaticSystemAst {
        charge,
        spin,
        electrons,
        constraints: _,
    } = ast;

    if matches!(*charge, ValueAst::Undetermined) {
        *charge = match cfg.charge{
            NumericDefault::Zero => ValueAst::Lit(0),
            NumericDefault::Required => ValueAst::Undetermined,
        };
    }
    if matches!(*electrons, ValueAst::Undetermined) {
        *electrons = match cfg.electrons{
            NumericDefault::Zero => ValueAst::Lit(0),
            NumericDefault::Required => ValueAst::Undetermined,
        };
    }
    raise_spin(spin, cfg.unpaired_electrons, cfg.multiplicity);
}

// -- Lower --------------------

fn lower_aromatic(ast: &mut AromaticSystemAst, cfg: &AromaticSystemDefaults) {
    // Exhaustive destructure: adding a new AromaticSystemAst field is a
    // compile error here, forcing the author to decide how lowering should
    // handle it.
    let AromaticSystemAst {
        charge,
        spin,
        electrons,
        constraints: _,
    } = ast;

    if matches!(
        (&cfg.charge, &*charge),
        (NumericDefault::Zero, ValueAst::Lit(0))
    ) {
        *charge = ValueAst::Undetermined;
    }
    if matches!(
        (&cfg.electrons, &*electrons),
        (NumericDefault::Zero, ValueAst::Lit(0))
    ) {
        *electrons = ValueAst::Undetermined;
    }
    lower_spin(spin, cfg.unpaired_electrons, cfg.multiplicity);
}

// -- Constraint DSL -------------------

/// Surface DSL wrapper around `AromaticSystemConstraint`. Mirrors the AST
/// enum with atom refs in place of `AtomIdx`. EDN form is a single-key map
/// keyed by the constraint kind (e.g. `{:atoms [0 1 2]}`, `{:contains :c1}`,
/// `{:all-atoms {:valence 4}}`).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum AromaticSystemConstraintDsl {
    Atoms(Vec<AtomRef>),
    Contains(AtomRef),
    ContainsAll(Vec<AtomRef>),
    AllAtoms(Box<AtomConstraintDsl>),
    AnyAtom(Box<AtomConstraintDsl>),
}

impl FromAst<AromaticSystemConstraint> for AromaticSystemConstraintDsl {
    type Ctx<'a> = ResolveContext<'a>;
    type Error = Infallible;

    fn from_ast<'a>(
        c: &AromaticSystemConstraint,
        ctx: &Self::Ctx<'a>,
    ) -> Result<Self, Infallible> {
        let meta = ctx.metadata;
        Ok(match c {
            AromaticSystemConstraint::Atoms(atoms) => {
                Self::Atoms(atoms.iter().map(|&a| AtomRef::from_ast(a, meta)).collect())
            }
            AromaticSystemConstraint::Contains(a) => {
                Self::Contains(AtomRef::from_ast(*a, meta))
            }
            AromaticSystemConstraint::ContainsAll(atoms) => Self::ContainsAll(
                atoms.iter().map(|&a| AtomRef::from_ast(a, meta)).collect(),
            ),
            AromaticSystemConstraint::AllAtoms(c) => {
                Self::AllAtoms(Box::new(AtomConstraintDsl::from_ast(c, &()).unwrap()))
            }
            AromaticSystemConstraint::AnyAtom(c) => {
                Self::AnyAtom(Box::new(AtomConstraintDsl::from_ast(c, &()).unwrap()))
            }
        })
    }
}

impl IntoAst<AromaticSystemConstraint> for AromaticSystemConstraintDsl {
    type Ctx<'a> = ResolveContext<'a>;
    type Error = ParseError;

    fn into_ast<'a>(
        self,
        ctx: &Self::Ctx<'a>,
    ) -> Result<AromaticSystemConstraint, ParseError> {
        let meta = ctx.metadata;
        let resolve_atoms = |refs: Vec<AtomRef>| -> Result<Vec<_>, ParseError> {
            refs.into_iter()
                .map(|r| r.into_ast(ctx.atom_count, meta))
                .collect()
        };
        Ok(match self {
            Self::Atoms(refs) => AromaticSystemConstraint::Atoms(resolve_atoms(refs)?),
            Self::Contains(r) => {
                AromaticSystemConstraint::Contains(r.into_ast(ctx.atom_count, meta)?)
            }
            Self::ContainsAll(refs) => {
                AromaticSystemConstraint::ContainsAll(resolve_atoms(refs)?)
            }
            Self::AllAtoms(c) => {
                AromaticSystemConstraint::AllAtoms(Box::new(c.into_ast(&()).unwrap()))
            }
            Self::AnyAtom(c) => {
                AromaticSystemConstraint::AnyAtom(Box::new(c.into_ast(&()).unwrap()))
            }
        })
    }
}

impl<'de> FromEdn<'de> for AromaticSystemConstraintDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let Edn::Map(m) = edn else {
            return Err(DeError::TypeMismatch {
                expected: "aromatic-system-constraint single-key map",
                got: edn.kind(),
                path: Vec::new(),
            });
        };
        if m.len() != 1 {
            return Err(DeError::Custom(format!(
                "aromatic-system-constraint must have exactly one key, got {}",
                m.len()
            )));
        }
        let (k, v) = m.iter().next().unwrap();
        let Edn::Keyword(key) = k else {
            return Err(DeError::TypeMismatch {
                expected: "keyword key",
                got: k.kind(),
                path: vec!["aromatic-system-constraint".into()],
            });
        };
        Ok(match key.name() {
            "atoms" => Self::Atoms(parse_atom_refs(v)?),
            "contains" => Self::Contains(AtomRef::from_edn(v)?),
            "contains-all" => Self::ContainsAll(parse_atom_refs(v)?),
            "all-atoms" => Self::AllAtoms(Box::new(AtomConstraintDsl::from_edn(v)?)),
            "any-atom" => Self::AnyAtom(Box::new(AtomConstraintDsl::from_edn(v)?)),
            other => {
                return Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["aromatic-system-constraint".into()],
                });
            }
        })
    }
}

impl ToEdn for AromaticSystemConstraintDsl {
    fn to_edn(&self) -> Edn<'static> {
        let (key, value) = match self {
            Self::Atoms(refs) => ("atoms", render_atom_refs(refs)),
            Self::Contains(r) => ("contains", r.to_edn()),
            Self::ContainsAll(refs) => ("contains-all", render_atom_refs(refs)),
            Self::AllAtoms(c) => ("all-atoms", c.to_edn()),
            Self::AnyAtom(c) => ("any-atom", c.to_edn()),
        };
        let mut m = umol_edn::EdnMap::with_capacity(1);
        m.insert(Edn::Keyword(umol_edn::EdnKeyword::owned(key.into())), value);
        Edn::Map(m)
    }
}

fn parse_atom_refs(edn: &Edn<'_>) -> Result<Vec<AtomRef>, DeError> {
    let Edn::Vector(v) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "vector of atom refs",
            got: edn.kind(),
            path: Vec::new(),
        });
    };
    v.iter().map(AtomRef::from_edn).collect()
}

fn render_atom_refs(refs: &[AtomRef]) -> Edn<'static> {
    Edn::Vector(refs.iter().map(AtomRef::to_edn).collect::<Vec<_>>().into())
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ast::constraint::AromaticSystemConstraints;
    use crate::ast::spin::SpinStateAst;

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("", AromaticSystemDsl(AromaticSystemAst::default()))]
    #[case::whitespace("   ", AromaticSystemDsl(AromaticSystemAst::default()))]
    #[case::charge_pos("#c+1", AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Lit(1), spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: AromaticSystemConstraints::new() }))]
    #[case::charge_neg("#c-2", AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Lit(-2), spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: AromaticSystemConstraints::new() }))]
    #[case::charge_plus_only("#c+", AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Lit(1), spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: AromaticSystemConstraints::new() }))]
    #[case::charge_minus_only("#c-", AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Lit(-1), spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: AromaticSystemConstraints::new() }))]
    #[case::electrons("#e6", AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Undetermined, spin: SpinStateAst::default(), electrons: ValueAst::Lit(6), constraints: AromaticSystemConstraints::new() }))]
    #[case::electrons_bare("#e", AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Undetermined, spin: SpinStateAst::default(), electrons: ValueAst::Lit(1), constraints: AromaticSystemConstraints::new() }))]
    #[case::electrons_wild("#e*", AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Undetermined, spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: AromaticSystemConstraints::new() }))]
    #[case::unpaired("#u1", AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Lit(1), multiplicity: ValueAst::Undetermined }, electrons: ValueAst::Undetermined, constraints: AromaticSystemConstraints::new() }))]
    #[case::mult("#s2", AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(2) }, electrons: ValueAst::Undetermined, constraints: AromaticSystemConstraints::new() }))]
    #[case::charge_electrons("#c+#e6", AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Lit(1), spin: SpinStateAst::default(), electrons: ValueAst::Lit(6), constraints: AromaticSystemConstraints::new() }))]
    #[case::full("#c0#u0#s1#e6", AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Lit(0), spin: SpinStateAst::new(0, 1), electrons: ValueAst::Lit(6), constraints: AromaticSystemConstraints::new() }))]
    fn test_parse_aromatic(#[case] input: &str, #[case] expected: AromaticSystemDsl) {
        let result = aromatic.parse(input);
        assert!(result.is_ok(), "{:?} should succeed, got {:?}", input, result.clone().unwrap_err());
        let form = result.unwrap();
        assert_eq!(form, expected);
    }

    #[rstest]
    #[case::unknown("#x", ParseError::UnknownAromaticPredicate("#x".to_string()))]
    #[case::unknown_a("#a", ParseError::UnknownAromaticPredicate("#a".to_string()))]
    #[case::dup_charge("#c+#c-", ParseError::DuplicateAromaticPredicate("#c".to_string()))]
    #[case::dup_electrons("#e6#e4", ParseError::DuplicateAromaticPredicate("#e".to_string()))]
    #[case::dup_unpaired("#u1#u2", ParseError::DuplicateAromaticPredicate("#u".to_string()))]
    #[case::dup_multiplicity("#s1#s2", ParseError::DuplicateAromaticPredicate("#s".to_string()))]
    #[case::trailing("#c+ foo", ParseError::TrailingInput("foo".to_string()))]
    fn test_parse_aromatic_error(#[case] input: &str, #[case] expected: ParseError) {
        let result = aromatic.parse(input);
        assert!(
            result.is_err(),
            "{:?} should fail, got {:?}",
            input,
            result.unwrap()
        );
        let err = result.unwrap_err().into_inner();
        assert_eq!(err, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(AromaticSystemDsl::default(), "")]
    #[case::charge_one(AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Lit(1), spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: AromaticSystemConstraints::new() }), "#c+")]
    #[case::charge_neg_two(AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Lit(-2), spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: AromaticSystemConstraints::new() }), "#c-2")]
    #[case::electrons_six(AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Undetermined, spin: SpinStateAst::default(), electrons: ValueAst::Lit(6), constraints: AromaticSystemConstraints::new() }), "#e6")]
    #[case::electrons_one(AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Undetermined, spin: SpinStateAst::default(), electrons: ValueAst::Lit(1), constraints: AromaticSystemConstraints::new() }), "#e")]
    #[case::full(AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Lit(0), spin: SpinStateAst::new(0, 1), electrons: ValueAst::Lit(6), constraints: AromaticSystemConstraints::new() }), "#c0#u0#s#e6")]
    fn test_display_aromatic(#[case] form: AromaticSystemDsl, #[case] expected: &str) {
        assert_eq!(form.to_string(), expected);
    }

    #[rstest]
    #[case::empty("")]
    #[case::charge("#c+1")]
    #[case::electrons("#e6")]
    #[case::unpaired("#u2")]
    #[case::explicit_mult("#s2")]
    fn test_aromatic_roundtrip(#[case] input: &str) {
        let form: AromaticSystemDsl = input.parse().unwrap();
        let rendered = form.to_string();
        let reparsed: AromaticSystemDsl = rendered.parse().unwrap();
        assert_eq!(form, reparsed);
    }

    #[rstest]
    fn test_aromatic_dsl_to_ast_fills_zero_defaults() {
        let dsl = AromaticSystemDsl::default();
        let cfg = AromaticSystemDefaults::zeroed();
        let ast = dsl.into_ast(&cfg).unwrap();
        assert_eq!(ast.charge, ValueAst::Lit(0));
        assert_eq!(ast.electrons, ValueAst::Lit(0));
        assert_eq!(ast.spin, SpinStateAst::new(0, 1));
    }

    #[rstest]
    fn test_aromatic_dsl_from_ast_strips_zero_defaults() {
        let ast = AromaticSystemAst {
            charge: ValueAst::Lit(0),
            spin: SpinStateAst::new(0, 1),
            electrons: ValueAst::Lit(0),
            constraints: AromaticSystemConstraints::new(),
        };
        let cfg = AromaticSystemDefaults::zeroed();
        let dsl = AromaticSystemDsl::from_ast(&ast, &cfg).unwrap();
        assert_eq!(dsl.0.charge, ValueAst::Undetermined);
        assert_eq!(dsl.0.electrons, ValueAst::Undetermined);
        assert_eq!(dsl.0.spin, SpinStateAst::default());
    }

    #[rstest]
    fn test_aromatic_dsl_roundtrip_zeroed() {
        let input = AromaticSystemDsl::default();
        let cfg = AromaticSystemDefaults::zeroed();
        let raised = input.clone().into_ast(&cfg).unwrap();
        let lowered = AromaticSystemDsl::from_ast(&raised, &cfg).unwrap();
        assert_eq!(input, lowered);
    }

    #[rstest]
    #[case::empty(r##""""##)]
    #[case::charge(r##""#c+""##)]
    #[case::full(r##""#c0#u0#s1#e6""##)]
    fn test_aromatic_dsl_from_edn_str_matches_from_edn(#[case] input: &str) {
        let via_stream = AromaticSystemDsl::from_edn_str(input).unwrap();
        let tree = umol_edn::read_string(input).unwrap();
        let via_tree = AromaticSystemDsl::from_edn(&tree).unwrap();
        assert_eq!(via_stream, via_tree);
    }

    // -- AromaticSystemConstraintDsl ----------------

    use super::super::molecule::Metadata;
    use crate::ast::constraint::{AromaticSystemConstraint, AtomConstraint};
    use crate::ast::idx::AtomIdx;
    use bimap::BiMap;
    use indexmap::IndexMap;

    fn empty_metadata() -> Metadata {
        Metadata {
            atom_ids: IndexMap::new(),
            atom_aliases: BiMap::new(),
            bond_ids: IndexMap::new(),
            dative_bond_ids: IndexMap::new(),
            aromatic_system_ids: IndexMap::new(),
            multicenter_bond_ids: IndexMap::new(),
            noncovalent_bond_ids: IndexMap::new(),
        }
    }

    fn ctx_with_atoms(atom_count: usize, meta: &Metadata) -> ResolveContext<'_> {
        ResolveContext {
            atom_count,
            bond_count: 0,
            dative_bond_count: 0,
            aromatic_system_count: 0,
            multicenter_bond_count: 0,
            noncovalent_bond_count: 0,
            metadata: meta,
        }
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::atoms(AromaticSystemConstraint::Atoms(vec![AtomIdx(0), AtomIdx(1)]), "{:atoms [0 1]}")]
    #[case::contains(AromaticSystemConstraint::Contains(AtomIdx(3)), "{:contains 3}")]
    #[case::contains_all(AromaticSystemConstraint::ContainsAll(vec![AtomIdx(2), AtomIdx(5)]), "{:contains-all [2 5]}")]
    #[case::all_atoms(AromaticSystemConstraint::AllAtoms(Box::new(AtomConstraint::Valence(ValueAst::Lit(4)))), "{:all-atoms {:valence 4}}")]
    #[case::any_atom(AromaticSystemConstraint::AnyAtom(Box::new(AtomConstraint::Degree(ValueAst::Lit(3)))), "{:any-atom {:degree 3}}")]
    fn test_aromatic_system_constraint_dsl_roundtrip(
        #[case] input: AromaticSystemConstraint,
        #[case] edn_source: &str,
    ) {
        let meta = empty_metadata();
        let render_ctx = ResolveContext::for_rendering(&meta);
        let dsl = AromaticSystemConstraintDsl::from_ast(&input, &render_ctx).unwrap();
        let edn = dsl.clone().to_edn();
        let expected = umol_edn::read_string(edn_source).unwrap();
        assert_eq!(edn, expected, "render mismatch");
        let parsed = AromaticSystemConstraintDsl::from_edn(&edn).unwrap();
        let back = parsed.into_ast(&ctx_with_atoms(10, &meta)).unwrap();
        assert_eq!(back, input, "parse-back mismatch");
    }

    #[rstest]
    fn test_aromatic_system_constraint_dsl_rejects_out_of_range_atom() {
        let meta = empty_metadata();
        let edn = umol_edn::read_string("{:contains 99}").unwrap();
        let dsl = AromaticSystemConstraintDsl::from_edn(&edn).unwrap();
        let err = dsl.into_ast(&ctx_with_atoms(5, &meta)).unwrap_err();
        assert_eq!(
            err,
            ParseError::InvalidRef {
                kind: "atom",
                value: "99".into()
            }
        );
    }

    #[rstest]
    fn test_aromatic_system_constraint_dsl_rejects_unknown_key() {
        let edn = umol_edn::read_string("{:bogus 1}").unwrap();
        let err = AromaticSystemConstraintDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::UnknownField { .. }));
    }
}
