//! Bond-string DSL.

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use umol_edn::{DeError, Edn, FromEdn, ToEdn};
use winnow::ascii::multispace0;
use winnow::combinator::{preceded, repeat, terminated};
use winnow::error::ErrMode;
use winnow::token::take;
use winnow::Parser;

use super::error::{PResult, ParseError};
use super::predicates::{
    apply_spin_pair, charge, fmt_charge, fmt_ring_count, fmt_spin_pair, fmt_value, optional_value,
    ring_count, SpinPredicate,
};
use super::value::value;
use crate::ast::bond::BondAst;
use crate::ast::constraint::BondConstraint;
use crate::ast::value::ValueAst;

/// AST-layer bundle pairing a `BondAst` with the bond-level constraints that
/// round-trip through the bond-string DSL.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BondTypeAst {
    pub ast: BondAst,
    pub constraints: Vec<BondConstraint>,
}

impl BondTypeAst {
    pub fn new(ast: BondAst) -> Self {
        Self {
            ast,
            constraints: Vec::new(),
        }
    }

    pub fn with_constraints(ast: BondAst, constraints: Vec<BondConstraint>) -> Self {
        Self { ast, constraints }
    }
}

impl FromStr for BondTypeAst {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_bond(s)
    }
}

impl Display for BondTypeAst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_bond_ast(f, &self.ast)?;
        for c in &self.constraints {
            fmt_constraint(f, c)?;
        }
        Ok(())
    }
}

impl<'de> FromEdn<'de> for BondTypeAst {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s.parse().map_err(|e| DeError::subgrammar("bond", e)),
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for BondTypeAst {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

/// Parse a complete bond-string into a `BondTypeAst` (base AST + lifted constraints).
pub fn parse_bond(input: &str) -> Result<BondTypeAst, ParseError> {
    bond.parse(input).map_err(|e| e.into_inner())
}

/// Bond-string parser (does not require consuming all input).
pub(crate) fn bond(i: &mut &str) -> PResult<BondTypeAst> {
    let order = preceded(multispace0, terminated(value, multispace0)).parse_next(i)?;
    let preds: Vec<BondPredicate> =
        repeat(0.., terminated(bond_predicate, multispace0)).parse_next(i)?;
    let mut form = BondTypeAst::new(BondAst::new(order));
    apply_predicates(&mut form, preds).map_err(ErrMode::Cut)?;
    Ok(form)
}

fn apply_predicates(form: &mut BondTypeAst, preds: Vec<BondPredicate>) -> Result<(), ParseError> {
    let ast = &mut form.ast;
    for pred in preds {
        match pred {
            BondPredicate::Charge(v) => {
                if !matches!(ast.charge, ValueAst::Undetermined) {
                    return Err(ParseError::DuplicateBondPredicate("#c".to_string()));
                }
                ast.charge = v;
            }
            BondPredicate::Spin(sp) => {
                apply_spin_pair(&mut ast.spin, sp, ParseError::DuplicateBondPredicate)?;
            }
            BondPredicate::Constraint(c) => {
                let tag = constraint_tag(&c);
                if form
                    .constraints
                    .iter()
                    .any(|existing| constraint_tag(existing) == tag)
                {
                    return Err(ParseError::DuplicateBondPredicate(tag.to_string()));
                }
                form.constraints.push(c);
            }
        }
    }
    Ok(())
}

fn constraint_tag(c: &BondConstraint) -> &'static str {
    match c {
        BondConstraint::Aromatic => "#a",
        BondConstraint::RingCount(_) => "#R",
        BondConstraint::RingSize(_) => "#r",
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BondPredicate {
    Charge(ValueAst),
    Spin(SpinPredicate),
    Constraint(BondConstraint),
}

fn bond_predicate(i: &mut &str) -> PResult<BondPredicate> {
    let start = *i;
    let prefix: &str = take(2usize).parse_next(i)?;
    match prefix {
        "#c" => charge.map(BondPredicate::Charge).parse_next(i),
        "#u" => optional_value
            .map(|v| BondPredicate::Spin(SpinPredicate::Unpaired(v)))
            .parse_next(i),
        "#s" => optional_value
            .map(|v| BondPredicate::Spin(SpinPredicate::Multiplicity(v)))
            .parse_next(i),
        "#a" => Ok(BondPredicate::Constraint(BondConstraint::Aromatic)),
        "#R" => ring_count
            .map(|v| BondPredicate::Constraint(BondConstraint::RingCount(v)))
            .parse_next(i),
        "#r" => optional_value
            .map(|v| BondPredicate::Constraint(BondConstraint::RingSize(v)))
            .parse_next(i),
        p if p.starts_with('#') => Err(ErrMode::Cut(ParseError::UnknownBondPredicate(
            p.to_string(),
        ))),
        _ => Err(ErrMode::Cut(ParseError::TrailingInput(start.to_string()))),
    }
}

fn fmt_bond_ast(f: &mut fmt::Formatter<'_>, ast: &BondAst) -> fmt::Result {
    match &ast.order {
        ValueAst::Lit(n) => write!(f, "{}", n)?,
        ValueAst::Undetermined => write!(f, "*")?,
        v => fmt_value(f, v)?,
    }

    fmt_charge(f, &ast.charge)?;
    fmt_spin_pair(f, &ast.spin)
}

fn fmt_constraint(f: &mut fmt::Formatter<'_>, c: &BondConstraint) -> fmt::Result {
    match c {
        BondConstraint::Aromatic => write!(f, "#a"),
        BondConstraint::RingCount(v) => fmt_ring_count(f, v),
        BondConstraint::RingSize(v) => match v {
            ValueAst::Undetermined => write!(f, "#r*"),
            ValueAst::Lit(1) => write!(f, "#r"),
            ValueAst::Lit(n) => write!(f, "#r{}", n),
            v => {
                write!(f, "#r")?;
                fmt_value(f, v)
            }
        },
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ast::{spin::SpinStateAst, value::{Expr, RelOp}};

    #[rustfmt::skip]
    #[rstest]
    #[case::single("1", BondTypeAst::new(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: Vec::new() }))]
    #[case::double("2", BondTypeAst::new(BondAst { order: ValueAst::Lit(2), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: Vec::new() }))]
    #[case::triple("3", BondTypeAst::new(BondAst { order: ValueAst::Lit(3), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: Vec::new() }))]
    #[case::quadruple("4", BondTypeAst::new(BondAst { order: ValueAst::Lit(4), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: Vec::new() }))]
    #[case::single_whitespace("  1  ", BondTypeAst::new(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: Vec::new() }))]
    #[case::single_pos_charge("1#c+2", BondTypeAst::new(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(2), spin: SpinStateAst::default(), constraints: Vec::new() }))]
    #[case::single_neg_charge("1#c-2", BondTypeAst::new(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(-2), spin: SpinStateAst::default(), constraints: Vec::new() }))]
    #[case::single_zero_charge("1#c0", BondTypeAst::new(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(0), spin: SpinStateAst::default(), constraints: Vec::new() }))]
    #[case::single_plus_only("1#c+", BondTypeAst::new(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: Vec::new() }))]
    #[case::single_minus_only("1#c-", BondTypeAst::new(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(-1), spin: SpinStateAst::default(), constraints: Vec::new() }))]
    #[case::double_unpaired("2#u3", BondTypeAst::new(BondAst { order: ValueAst::Lit(2), charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Lit(3), multiplicity: ValueAst::Undetermined }, constraints: Vec::new() }))]
    #[case::single_u_only("1#u", BondTypeAst::new(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Lit(1), multiplicity: ValueAst::Undetermined }, constraints: Vec::new() }))]
    #[case::single_mult("1#s2", BondTypeAst::new(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(2) }, constraints: Vec::new() }))]
    #[case::single_s_only("1#s", BondTypeAst::new(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(1) }, constraints: Vec::new() }))]
    #[case::double_charge_unpaired("2#c+#u2", BondTypeAst::new(BondAst { order: ValueAst::Lit(2), charge: ValueAst::Lit(1), spin: SpinStateAst { unpaired: ValueAst::Lit(2), multiplicity: ValueAst::Undetermined }, constraints: Vec::new() }))]
    #[case::double_charge_mult("2#c-1#s3", BondTypeAst::new(BondAst { order: ValueAst::Lit(2), charge: ValueAst::Lit(-1), spin: SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(3) }, constraints: Vec::new() }))]
    #[case::aromatic("1#a", BondTypeAst::with_constraints(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: Vec::new() }, vec![BondConstraint::Aromatic]))]
    #[case::charged_aromatic("1#c+#a", BondTypeAst::with_constraints(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: Vec::new() }, vec![BondConstraint::Aromatic]))]
    #[case::ring_count("1#R2", BondTypeAst::with_constraints(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: Vec::new() }, vec![BondConstraint::RingCount(ValueAst::Lit(2))]))]
    #[case::ring_bare("1#R", BondTypeAst::with_constraints(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: Vec::new() }, vec![BondConstraint::RingCount(ValueAst::Lit(1))]))]
    #[case::ring_plus("1#R+", BondTypeAst::with_constraints(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: Vec::new() }, vec![BondConstraint::RingCount(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("r".to_string())), RelOp::Ge, Box::new(Expr::Lit(1)))))]))]
    #[case::ring_undetermined("1#R*", BondTypeAst::with_constraints(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: Vec::new() }, vec![BondConstraint::RingCount(ValueAst::Undetermined)]))]
    #[case::ring_size("1#r6", BondTypeAst::with_constraints(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: Vec::new() }, vec![BondConstraint::RingSize(ValueAst::Lit(6))]))]
    fn test_parse_bond(#[case] input: &str, #[case] expected: BondTypeAst) {
        let result = bond.parse(input);
        assert!(result.is_ok(), "{:?} should succeed, got {:?}", input, result.clone().unwrap_err());
        let form = result.unwrap();
        assert_eq!(form, expected);
    }

    #[rstest]
    #[case::empty("", ParseError::Syntax)]
    #[case::unknown_pred("1#x", ParseError::UnknownBondPredicate("#x".to_string()))]
    #[case::dup_charge("1#c+#c-", ParseError::DuplicateBondPredicate("#c".to_string()))]
    #[case::dup_unpaired("1#u2#u3", ParseError::DuplicateBondPredicate("#u".to_string()))]
    #[case::dup_multiplicity("1#s1#s2", ParseError::DuplicateBondPredicate("#s".to_string()))]
    #[case::dup_aromatic("1#a#a", ParseError::DuplicateBondPredicate("#a".to_string()))]
    #[case::trailing("1#c+ foo", ParseError::TrailingInput("foo".to_string()))]
    fn test_parse_bond_error(#[case] input: &str, #[case] expected: ParseError) {
        let result = bond.parse(input);
        assert!(
            result.is_err(),
            "{:?} should fail, got {:?}",
            input,
            result.unwrap()
        );
        let err = result.unwrap_err().into_inner();
        assert_eq!(
            err, expected,
            "{:?} should fail with {:?}, got {:?}",
            input, expected, err
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::charge_pos("#c+2", BondPredicate::Charge(ValueAst::Lit(2)))]
    #[case::charge_neg("#c-2", BondPredicate::Charge(ValueAst::Lit(-2)))]
    #[case::charge_plus("#c+", BondPredicate::Charge(ValueAst::Lit(1)))]
    #[case::charge_minus("#c-", BondPredicate::Charge(ValueAst::Lit(-1)))]
    #[case::charge_zero("#c0", BondPredicate::Charge(ValueAst::Lit(0)))]
    #[case::charge_undetermined("#c*", BondPredicate::Charge(ValueAst::Undetermined))]
    #[case::unpaired("#u2", BondPredicate::Spin(SpinPredicate::Unpaired(ValueAst::Lit(2))))]
    #[case::unpaired_omit("#u", BondPredicate::Spin(SpinPredicate::Unpaired(ValueAst::Lit(1))))]
    #[case::unpaired_undetermined("#u*", BondPredicate::Spin(SpinPredicate::Unpaired(ValueAst::Undetermined)))]
    #[case::multiplicity("#s3", BondPredicate::Spin(SpinPredicate::Multiplicity(ValueAst::Lit(3))))]
    #[case::multiplicity_omit("#s", BondPredicate::Spin(SpinPredicate::Multiplicity(ValueAst::Lit(1))))]
    #[case::multiplicity_undetermined("#s*", BondPredicate::Spin(SpinPredicate::Multiplicity(ValueAst::Undetermined)))]
    #[case::aromatic("#a", BondPredicate::Constraint(BondConstraint::Aromatic))]
    fn test_bond_predicate(#[case] input: &str, #[case] expected: BondPredicate) {
        let result = bond_predicate.parse(input);
        assert!(result.is_ok(), "{input:?} should succeed, got {:?}", result.unwrap_err());
        let pred = result.unwrap();
        assert_eq!(pred, expected);
    }

    #[rstest]
    #[case::unknown("#x", ParseError::UnknownBondPredicate("#x".to_string()))]
    #[case::unknown_tag("#z", ParseError::UnknownBondPredicate("#z".to_string()))]
    #[case::trailing_no_hash("fo", ParseError::TrailingInput("fo".to_string()))]
    fn test_bond_predicate_error(#[case] input: &str, #[case] expected: ParseError) {
        let result = bond_predicate.parse(input);
        assert!(
            result.is_err(),
            "{input:?} should fail, got {:?}",
            result.unwrap()
        );
        let err = result.unwrap_err().into_inner();
        assert_eq!(err, expected);
    }
}
