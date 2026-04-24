//! Multicenter-bond-string DSL.

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
use super::constraint::{AtomRef, EntityCounts};
use super::molecule::Metadata;
use super::error::{PResult, ParseError};
use super::predicates::{
    apply_spin_pair, charge, fmt_charge, fmt_spin_pair, lower_spin, optional_value, raise_spin,
    SpinPredicate,
};
use super::value::fmt_value;
use crate::ast::constraint::MulticenterBondConstraint;
use crate::ast::multicenter::MulticenterBondAst;
use crate::ast::traits::{FromAst, IntoAst};
use crate::ast::value::ValueAst;
use crate::dsl::config::{MulticenterBondDefaults, NumericDefault};

/// Surface DSL wrapper around `MulticenterBondAst`. Parses and renders the
/// multicenter-bond-string form. All `MulticenterBondConstraint` variants are
/// molecule-scope, so nothing from the constraint vec serializes inline.
#[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MulticenterBondDsl(pub MulticenterBondAst);

impl MulticenterBondDsl {
    /// Zero-cost reference cast from `&MulticenterBondAst`. Relies on `repr(transparent)`.
    pub fn from_ref(ast: &MulticenterBondAst) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(ast as *const MulticenterBondAst as *const Self) }
    }
}

impl FromStr for MulticenterBondDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_multicenter(s)
    }
}

impl Display for MulticenterBondDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_multicenter_ast(f, &self.0)
    }
}

impl<'de> FromEdn<'de> for MulticenterBondDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s.parse().map_err(|e| DeError::subgrammar("multicenter", e)),
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        EdnStreamDeserializer::new(input).read_subgrammar_all("multicenter")
    }
}

impl ToEdn for MulticenterBondDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

impl FromAst<MulticenterBondAst> for MulticenterBondDsl {
    type Ctx = MulticenterBondDefaults;
    type Error = ParseError;

    fn from_ast(ast: &MulticenterBondAst, cfg: &Self::Ctx) -> Result<Self, ParseError> {
        let mut out = ast.clone();
        lower_multicenter(&mut out, cfg);
        Ok(MulticenterBondDsl(out))
    }
}

impl IntoAst<MulticenterBondAst> for MulticenterBondDsl {
    type Ctx = MulticenterBondDefaults;
    type Error = ParseError;

    fn into_ast(mut self, cfg: &Self::Ctx) -> Result<MulticenterBondAst, ParseError> {
        raise_multicenter(&mut self.0, cfg);
        Ok(self.0)
    }
}

// -- Parse --------------------

pub fn parse_multicenter(input: &str) -> Result<MulticenterBondDsl, ParseError> {
    multicenter.parse(input).map_err(|e| e.into_inner())
}

pub(crate) fn multicenter(i: &mut &str) -> PResult<MulticenterBondDsl> {
    multispace0.parse_next(i)?;
    let preds: Vec<MulticenterPredicate> =
        repeat(0.., terminated(multicenter_predicate, multispace0)).parse_next(i)?;
    let mut form = MulticenterBondDsl::default();
    apply_predicates(&mut form, preds).map_err(ErrMode::Cut)?;
    Ok(form)
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MulticenterPredicate {
    Charge(ValueAst),
    Spin(SpinPredicate),
    Electrons(ValueAst),
}

fn multicenter_predicate(i: &mut &str) -> PResult<MulticenterPredicate> {
    let start = *i;
    let prefix: &str = take(2usize).parse_next(i)?;
    match prefix {
        "#c" => charge.map(MulticenterPredicate::Charge).parse_next(i),
        "#u" => optional_value
            .map(|v| MulticenterPredicate::Spin(SpinPredicate::Unpaired(v)))
            .parse_next(i),
        "#s" => optional_value
            .map(|v| MulticenterPredicate::Spin(SpinPredicate::Multiplicity(v)))
            .parse_next(i),
        "#e" => optional_value
            .map(MulticenterPredicate::Electrons)
            .parse_next(i),
        p if p.starts_with('#') => Err(ErrMode::Cut(ParseError::UnknownMulticenterPredicate(
            p.to_string(),
        ))),
        _ => Err(ErrMode::Cut(ParseError::TrailingInput(start.to_string()))),
    }
}

fn apply_predicates(
    form: &mut MulticenterBondDsl,
    preds: Vec<MulticenterPredicate>,
) -> Result<(), ParseError> {
    let ast = &mut form.0;
    for pred in preds {
        match pred {
            MulticenterPredicate::Charge(v) => {
                if !matches!(ast.charge, ValueAst::Undetermined) {
                    return Err(ParseError::DuplicateMulticenterPredicate("#c".to_string()));
                }
                ast.charge = v;
            }
            MulticenterPredicate::Spin(sp) => {
                apply_spin_pair(&mut ast.spin, sp, ParseError::DuplicateMulticenterPredicate)?;
            }
            MulticenterPredicate::Electrons(v) => {
                if !matches!(ast.electrons, ValueAst::Undetermined) {
                    return Err(ParseError::DuplicateMulticenterPredicate("#e".to_string()));
                }
                ast.electrons = v;
            }
        }
    }
    Ok(())
}

// -- Format --------------------

fn fmt_multicenter_ast(f: &mut fmt::Formatter<'_>, ast: &MulticenterBondAst) -> fmt::Result {
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

fn raise_multicenter(ast: &mut MulticenterBondAst, cfg: &MulticenterBondDefaults) {
    // Exhaustive destructure: adding a new MulticenterBondAst field is a
    // compile error here, forcing the author to decide how raising should
    // handle it.
    let MulticenterBondAst {
        charge,
        spin,
        electrons,
        constraints: _,
    } = ast;

    if matches!(*charge, ValueAst::Undetermined) {
        *charge = match cfg.charge {
            NumericDefault::Zero => ValueAst::Lit(0),
            NumericDefault::Required => ValueAst::Undetermined,
        };
    }
    if matches!(*electrons, ValueAst::Undetermined) {
        *electrons = match cfg.electrons {
            NumericDefault::Zero => ValueAst::Lit(0),
            NumericDefault::Required => ValueAst::Undetermined,
        };
    }
    raise_spin(spin, cfg.unpaired_electrons, cfg.multiplicity);
}

// -- Format --------------------

fn lower_multicenter(ast: &mut MulticenterBondAst, cfg: &MulticenterBondDefaults) {
    // Exhaustive destructure: adding a new MulticenterBondAst field is a
    // compile error here, forcing the author to decide how lowering should
    // handle it.
    let MulticenterBondAst {
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

/// Surface DSL wrapper around `MulticenterBondConstraint`. Mirrors the AST
/// enum with atom refs in place of `AtomIdx`. EDN form is a single-key map
/// keyed by the constraint kind (`:atoms`, `:contains`, `:contains-all`,
/// `:all-atoms`, `:any-atom`).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum MulticenterBondConstraintDsl {
    Atoms(Vec<AtomRef>),
    Contains(AtomRef),
    ContainsAll(Vec<AtomRef>),
    AllAtoms(Box<AtomConstraintDsl>),
    AnyAtom(Box<AtomConstraintDsl>),
}

impl MulticenterBondConstraintDsl {
    pub(crate) fn from_ast(
        c: &MulticenterBondConstraint,
        meta: &Metadata,
    ) -> Result<Self, Infallible> {
        Ok(match c {
            MulticenterBondConstraint::Atoms(atoms) => {
                Self::Atoms(atoms.iter().map(|&a| AtomRef::from_ast(a, meta)).collect())
            }
            MulticenterBondConstraint::Contains(a) => Self::Contains(AtomRef::from_ast(*a, meta)),
            MulticenterBondConstraint::ContainsAll(atoms) => {
                Self::ContainsAll(atoms.iter().map(|&a| AtomRef::from_ast(a, meta)).collect())
            }
            MulticenterBondConstraint::AllAtoms(c) => {
                Self::AllAtoms(Box::new(AtomConstraintDsl::from_ast(c, &()).unwrap()))
            }
            MulticenterBondConstraint::AnyAtom(c) => {
                Self::AnyAtom(Box::new(AtomConstraintDsl::from_ast(c, &()).unwrap()))
            }
        })
    }

    pub(crate) fn into_ast(
        self,
        counts: &EntityCounts,
        meta: &Metadata,
    ) -> Result<MulticenterBondConstraint, ParseError> {
        let resolve_atoms = |refs: Vec<AtomRef>| -> Result<Vec<_>, ParseError> {
            refs.into_iter()
                .map(|r| r.into_ast(counts.atom_count, meta))
                .collect()
        };
        Ok(match self {
            Self::Atoms(refs) => MulticenterBondConstraint::Atoms(resolve_atoms(refs)?),
            Self::Contains(r) => {
                MulticenterBondConstraint::Contains(r.into_ast(counts.atom_count, meta)?)
            }
            Self::ContainsAll(refs) => MulticenterBondConstraint::ContainsAll(resolve_atoms(refs)?),
            Self::AllAtoms(c) => {
                MulticenterBondConstraint::AllAtoms(Box::new(c.into_ast(&()).unwrap()))
            }
            Self::AnyAtom(c) => {
                MulticenterBondConstraint::AnyAtom(Box::new(c.into_ast(&()).unwrap()))
            }
        })
    }
}

impl<'de> FromEdn<'de> for MulticenterBondConstraintDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let Edn::Map(m) = edn else {
            return Err(DeError::TypeMismatch {
                expected: "multicenter-bond-constraint single-key map",
                got: edn.kind(),
                path: Vec::new(),
            });
        };
        if m.len() != 1 {
            return Err(DeError::Custom(format!(
                "multicenter-bond-constraint must have exactly one key, got {}",
                m.len()
            )));
        }
        let (k, v) = m.iter().next().unwrap();
        let Edn::Keyword(key) = k else {
            return Err(DeError::TypeMismatch {
                expected: "keyword key",
                got: k.kind(),
                path: vec!["multicenter-bond-constraint".into()],
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
                    path: vec!["multicenter-bond-constraint".into()],
                });
            }
        })
    }
}

impl ToEdn for MulticenterBondConstraintDsl {
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

    use super::super::molecule::Metadata;
    use super::*;
    use crate::ast::constraint::{
        AtomConstraint, MulticenterBondConstraint, MulticenterBondConstraints,
    };
    use crate::ast::idx::AtomIdx;
    use crate::ast::spin::SpinStateAst;

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("", MulticenterBondDsl(MulticenterBondAst::default()))]
    #[case::whitespace("   ", MulticenterBondDsl(MulticenterBondAst::default()))]
    #[case::charge_pos("#c+1", MulticenterBondDsl(MulticenterBondAst { charge: ValueAst::Lit(1), spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: MulticenterBondConstraints::new() }))]
    #[case::charge_neg("#c-2", MulticenterBondDsl(MulticenterBondAst { charge: ValueAst::Lit(-2), spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: MulticenterBondConstraints::new() }))]
    #[case::electrons("#e6", MulticenterBondDsl(MulticenterBondAst { charge: ValueAst::Undetermined, spin: SpinStateAst::default(), electrons: ValueAst::Lit(6), constraints: MulticenterBondConstraints::new() }))]
    #[case::electrons_bare("#e", MulticenterBondDsl(MulticenterBondAst { charge: ValueAst::Undetermined, spin: SpinStateAst::default(), electrons: ValueAst::Lit(1), constraints: MulticenterBondConstraints::new() }))]
    #[case::unpaired("#u1", MulticenterBondDsl(MulticenterBondAst { charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Lit(1), multiplicity: ValueAst::Undetermined }, electrons: ValueAst::Undetermined, constraints: MulticenterBondConstraints::new() }))]
    #[case::mult("#s2", MulticenterBondDsl(MulticenterBondAst { charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(2) }, electrons: ValueAst::Undetermined, constraints: MulticenterBondConstraints::new() }))]
    #[case::charge_electrons("#c+#e2", MulticenterBondDsl(MulticenterBondAst { charge: ValueAst::Lit(1), spin: SpinStateAst::default(), electrons: ValueAst::Lit(2), constraints: MulticenterBondConstraints::new() }))]
    #[case::full("#c0#u0#s1#e2", MulticenterBondDsl(MulticenterBondAst { charge: ValueAst::Lit(0), spin: SpinStateAst::new(0, 1), electrons: ValueAst::Lit(2), constraints: MulticenterBondConstraints::new() }))]
    fn test_parse_multicenter(#[case] input: &str, #[case] expected: MulticenterBondDsl) {
        let result = multicenter.parse(input);
        assert!(result.is_ok(), "{:?} should succeed, got {:?}", input, result.clone().unwrap_err());
        let form = result.unwrap();
        assert_eq!(form, expected);
    }

    #[rstest]
    #[case::unknown("#x", ParseError::UnknownMulticenterPredicate("#x".to_string()))]
    #[case::unknown_a("#a", ParseError::UnknownMulticenterPredicate("#a".to_string()))]
    #[case::dup_charge("#c+#c-", ParseError::DuplicateMulticenterPredicate("#c".to_string()))]
    #[case::dup_electrons("#e2#e4", ParseError::DuplicateMulticenterPredicate("#e".to_string()))]
    #[case::dup_unpaired("#u1#u2", ParseError::DuplicateMulticenterPredicate("#u".to_string()))]
    #[case::dup_multiplicity("#s1#s2", ParseError::DuplicateMulticenterPredicate("#s".to_string()))]
    #[case::trailing("#c+ foo", ParseError::TrailingInput("foo".to_string()))]
    fn test_parse_multicenter_error(#[case] input: &str, #[case] expected: ParseError) {
        let result = multicenter.parse(input);
        assert!(result.is_err(), "{:?} should fail", input);
        let err = result.unwrap_err().into_inner();
        assert_eq!(err, expected);
    }

    #[rstest]
    #[case::empty("")]
    #[case::charge("#c+1")]
    #[case::electrons("#e6")]
    #[case::unpaired("#u2")]
    #[case::explicit_mult("#s2")]
    fn test_multicenter_roundtrip(#[case] input: &str) {
        let form: MulticenterBondDsl = input.parse().unwrap();
        let rendered = form.to_string();
        let reparsed: MulticenterBondDsl = rendered.parse().unwrap();
        assert_eq!(form, reparsed);
    }

    #[rstest]
    fn test_multicenter_dsl_to_ast_fills_zero_defaults() {
        let dsl = MulticenterBondDsl::default();
        let cfg = MulticenterBondDefaults::zeroed();
        let ast = dsl.into_ast(&cfg).unwrap();
        assert_eq!(ast.charge, ValueAst::Lit(0));
        assert_eq!(ast.electrons, ValueAst::Lit(0));
        assert_eq!(ast.spin, SpinStateAst::new(0, 1));
    }

    #[rstest]
    fn test_multicenter_dsl_from_ast_strips_zero_defaults() {
        let ast = MulticenterBondAst {
            charge: ValueAst::Lit(0),
            spin: SpinStateAst::new(0, 1),
            electrons: ValueAst::Lit(0),
            constraints: MulticenterBondConstraints::new(),
        };
        let cfg = MulticenterBondDefaults::zeroed();
        let dsl = MulticenterBondDsl::from_ast(&ast, &cfg).unwrap();
        assert_eq!(dsl.0.charge, ValueAst::Undetermined);
        assert_eq!(dsl.0.electrons, ValueAst::Undetermined);
        assert_eq!(dsl.0.spin, SpinStateAst::default());
    }

    #[rstest]
    #[case::empty(r##""""##)]
    #[case::charge(r##""#c+""##)]
    #[case::full(r##""#c0#u0#s1#e2""##)]
    fn test_multicenter_dsl_from_edn_str_matches_from_edn(#[case] input: &str) {
        let via_stream = MulticenterBondDsl::from_edn_str(input).unwrap();
        let tree = umol_edn::read_string(input).unwrap();
        let via_tree = MulticenterBondDsl::from_edn(&tree).unwrap();
        assert_eq!(via_stream, via_tree);
    }

    // -- MulticenterBondConstraintDsl ----------------

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
    #[case::atoms(MulticenterBondConstraint::Atoms(vec![AtomIdx(0), AtomIdx(1)]), "{:atoms [0 1]}")]
    #[case::contains(MulticenterBondConstraint::Contains(AtomIdx(3)), "{:contains 3}")]
    #[case::contains_all(MulticenterBondConstraint::ContainsAll(vec![AtomIdx(2), AtomIdx(5)]), "{:contains-all [2 5]}")]
    #[case::all_atoms(MulticenterBondConstraint::AllAtoms(Box::new(AtomConstraint::Valence(ValueAst::Lit(4)))), "{:all-atoms {:valence 4}}")]
    #[case::any_atom(MulticenterBondConstraint::AnyAtom(Box::new(AtomConstraint::Degree(ValueAst::Lit(3)))), "{:any-atom {:degree 3}}")]
    fn test_multicenter_bond_constraint_dsl_roundtrip(
        #[case] input: MulticenterBondConstraint,
        #[case] edn_source: &str,
    ) {
        let meta = Metadata::default();
        let dsl = MulticenterBondConstraintDsl::from_ast(&input, &meta).unwrap();
        let edn = dsl.clone().to_edn();
        let expected = umol_edn::read_string(edn_source).unwrap();
        assert_eq!(edn, expected, "render mismatch");
        let parsed = MulticenterBondConstraintDsl::from_edn(&edn).unwrap();
        let back = parsed.into_ast(&counts_with_atoms(10), &meta).unwrap();
        assert_eq!(back, input, "parse-back mismatch");
    }

    #[rstest]
    fn test_multicenter_bond_constraint_dsl_rejects_out_of_range_atom() {
        let meta = Metadata::default();
        let edn = umol_edn::read_string("{:contains 99}").unwrap();
        let dsl = MulticenterBondConstraintDsl::from_edn(&edn).unwrap();
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
    fn test_multicenter_bond_constraint_dsl_rejects_unknown_key() {
        let edn = umol_edn::read_string("{:bogus 1}").unwrap();
        let err = MulticenterBondConstraintDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::UnknownField { .. }));
    }
}
