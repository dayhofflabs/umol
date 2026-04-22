//! Bond-string DSL.

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::mem;
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
use crate::ast::config::{
    BondAstConfig, MultiplicityMode, NumericMode, UnpairedElectronsMode,
};
use crate::ast::constraint::BondConstraint;
use crate::ast::spin::SpinStateAst;
use crate::ast::traits::{FromAst, ToAst};
use crate::ast::value::ValueAst;

/// Surface DSL wrapper around `BondAst`. Parses and renders the bond-string
/// form (order plus `#…` predicates); inline constraints land in
/// `self.0.constraints`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BondDsl(pub BondAst);

impl FromStr for BondDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_bond(s)
    }
}

impl Display for BondDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_bond_ast(f, &self.0)?;
        for c in &self.0.constraints {
            fmt_constraint(f, c)?;
        }
        Ok(())
    }
}

impl<'de> FromEdn<'de> for BondDsl {
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

impl ToEdn for BondDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

impl FromAst<BondAst> for BondDsl {
    type Error = ParseError;

    fn from_ast(ast: &BondAst, cfg: &BondAstConfig) -> Result<Self, ParseError> {
        let mut out = ast.clone();
        lower_bond(&mut out, cfg);
        Ok(BondDsl(out))
    }
}

impl ToAst<BondAst> for BondDsl {
    type Error = ParseError;

    fn to_ast(&self, cfg: &BondAstConfig) -> Result<BondAst, ParseError> {
        let mut out = self.0.clone();
        raise_bond(&mut out, cfg);
        Ok(out)
    }
}

fn raise_bond(ast: &mut BondAst, cfg: &BondAstConfig) {
    if matches!(ast.charge, ValueAst::Undetermined) {
        ast.charge = match cfg.charge_mode {
            NumericMode::Zero => ValueAst::Lit(0),
            NumericMode::Required => ValueAst::Undetermined,
        };
    }
    raise_spin(&mut ast.spin, cfg);
}

fn raise_spin(spin: &mut SpinStateAst, cfg: &BondAstConfig) {
    let u = mem::replace(&mut spin.unpaired, ValueAst::Undetermined);
    let m = mem::replace(&mut spin.multiplicity, ValueAst::Undetermined);
    let resolved_u = if matches!(u, ValueAst::Undetermined) {
        match cfg.unpaired_electrons_mode {
            UnpairedElectronsMode::Zero => ValueAst::Lit(0),
            UnpairedElectronsMode::Required => ValueAst::Undetermined,
            UnpairedElectronsMode::Derived => match &m {
                ValueAst::Lit(mm) => ValueAst::Lit(mm - 1),
                _ => ValueAst::Undetermined,
            },
        }
    } else {
        u
    };
    let resolved_m = if matches!(m, ValueAst::Undetermined) {
        match cfg.multiplicity_mode {
            MultiplicityMode::Required => ValueAst::Undetermined,
            MultiplicityMode::Derived => match &resolved_u {
                ValueAst::Lit(uu) => ValueAst::Lit(uu + 1),
                _ => ValueAst::Undetermined,
            },
        }
    } else {
        m
    };
    spin.unpaired = resolved_u;
    spin.multiplicity = resolved_m;
}

fn lower_bond(ast: &mut BondAst, cfg: &BondAstConfig) {
    if matches!(
        (&cfg.charge_mode, &ast.charge),
        (NumericMode::Zero, ValueAst::Lit(0))
    ) {
        ast.charge = ValueAst::Undetermined;
    }
    lower_spin(&mut ast.spin, cfg);
}

fn lower_spin(spin: &mut SpinStateAst, cfg: &BondAstConfig) {
    if let (ValueAst::Lit(uu), ValueAst::Lit(mm)) = (&spin.unpaired, &spin.multiplicity) {
        let derived = *mm == uu + 1;
        let strip_u = match cfg.unpaired_electrons_mode {
            UnpairedElectronsMode::Zero => *uu == 0,
            UnpairedElectronsMode::Derived => {
                derived && matches!(cfg.multiplicity_mode, MultiplicityMode::Derived)
            }
            UnpairedElectronsMode::Required => false,
        };
        let strip_m = matches!(cfg.multiplicity_mode, MultiplicityMode::Derived) && derived;
        if strip_u {
            spin.unpaired = ValueAst::Undetermined;
        }
        if strip_m {
            spin.multiplicity = ValueAst::Undetermined;
        }
    }
}

/// Parse a complete bond-string into a `BondDsl`.
pub fn parse_bond(input: &str) -> Result<BondDsl, ParseError> {
    bond.parse(input).map_err(|e| e.into_inner())
}

/// Bond-string parser (does not require consuming all input).
pub(crate) fn bond(i: &mut &str) -> PResult<BondDsl> {
    let order = preceded(multispace0, terminated(value, multispace0)).parse_next(i)?;
    let preds: Vec<BondPredicate> =
        repeat(0.., terminated(bond_predicate, multispace0)).parse_next(i)?;
    let mut form = BondDsl(BondAst::new(order));
    apply_predicates(&mut form, preds).map_err(ErrMode::Cut)?;
    Ok(form)
}

fn apply_predicates(form: &mut BondDsl, preds: Vec<BondPredicate>) -> Result<(), ParseError> {
    let ast = &mut form.0;
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
                if ast
                    .constraints
                    .iter()
                    .any(|existing| constraint_tag(existing) == tag)
                {
                    return Err(ParseError::DuplicateBondPredicate(tag.to_string()));
                }
                ast.constraints.push(c);
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

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ast::{spin::SpinStateAst, value::{Expr, RelOp}};

    #[rustfmt::skip]
    #[rstest]
    #[case::single("1", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: Vec::new() }))]
    #[case::double("2", BondDsl(BondAst { order: ValueAst::Lit(2), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: Vec::new() }))]
    #[case::triple("3", BondDsl(BondAst { order: ValueAst::Lit(3), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: Vec::new() }))]
    #[case::quadruple("4", BondDsl(BondAst { order: ValueAst::Lit(4), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: Vec::new() }))]
    #[case::single_whitespace("  1  ", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: Vec::new() }))]
    #[case::single_pos_charge("1#c+2", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(2), spin: SpinStateAst::default(), constraints: Vec::new() }))]
    #[case::single_neg_charge("1#c-2", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(-2), spin: SpinStateAst::default(), constraints: Vec::new() }))]
    #[case::single_zero_charge("1#c0", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(0), spin: SpinStateAst::default(), constraints: Vec::new() }))]
    #[case::single_plus_only("1#c+", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: Vec::new() }))]
    #[case::single_minus_only("1#c-", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(-1), spin: SpinStateAst::default(), constraints: Vec::new() }))]
    #[case::double_unpaired("2#u3", BondDsl(BondAst { order: ValueAst::Lit(2), charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Lit(3), multiplicity: ValueAst::Undetermined }, constraints: Vec::new() }))]
    #[case::single_u_only("1#u", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Lit(1), multiplicity: ValueAst::Undetermined }, constraints: Vec::new() }))]
    #[case::single_mult("1#s2", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(2) }, constraints: Vec::new() }))]
    #[case::single_s_only("1#s", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(1) }, constraints: Vec::new() }))]
    #[case::double_charge_unpaired("2#c+#u2", BondDsl(BondAst { order: ValueAst::Lit(2), charge: ValueAst::Lit(1), spin: SpinStateAst { unpaired: ValueAst::Lit(2), multiplicity: ValueAst::Undetermined }, constraints: Vec::new() }))]
    #[case::double_charge_mult("2#c-1#s3", BondDsl(BondAst { order: ValueAst::Lit(2), charge: ValueAst::Lit(-1), spin: SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(3) }, constraints: Vec::new() }))]
    #[case::aromatic("1#a", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: vec![BondConstraint::Aromatic] }))]
    #[case::charged_aromatic("1#c+#a", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: vec![BondConstraint::Aromatic] }))]
    #[case::ring_count("1#R2", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: vec![BondConstraint::RingCount(ValueAst::Lit(2))] }))]
    #[case::ring_bare("1#R", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: vec![BondConstraint::RingCount(ValueAst::Lit(1))] }))]
    #[case::ring_plus("1#R+", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: vec![BondConstraint::RingCount(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("r".to_string())), RelOp::Ge, Box::new(Expr::Lit(1)))))] }))]
    #[case::ring_undetermined("1#R*", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: vec![BondConstraint::RingCount(ValueAst::Undetermined)] }))]
    #[case::ring_size("1#r6", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: vec![BondConstraint::RingSize(ValueAst::Lit(6))] }))]
    fn test_parse_bond(#[case] input: &str, #[case] expected: BondDsl) {
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

    #[rstest]
    fn test_bond_dsl_to_ast_fills_zero_defaults() {
        let dsl = BondDsl(BondAst::new(ValueAst::Lit(1)));
        let cfg = BondAstConfig::zeroed();
        let ast = dsl.to_ast(&cfg).unwrap();
        assert_eq!(ast.charge, ValueAst::Lit(0));
        assert_eq!(ast.spin, SpinStateAst::new(0, 1));
    }

    #[rstest]
    fn test_bond_dsl_from_ast_strips_zero_defaults() {
        let mut ast = BondAst::new(ValueAst::Lit(1));
        ast.charge = ValueAst::Lit(0);
        ast.spin = SpinStateAst::new(0, 1);
        let cfg = BondAstConfig::zeroed();
        let dsl = BondDsl::from_ast(&ast, &cfg).unwrap();
        assert_eq!(dsl.0.charge, ValueAst::Undetermined);
        assert_eq!(dsl.0.spin, SpinStateAst::default());
    }

    #[rstest]
    fn test_bond_dsl_roundtrip_zeroed() {
        let input = BondDsl(BondAst::new(ValueAst::Lit(2)));
        let cfg = BondAstConfig::zeroed();
        let raised = input.to_ast(&cfg).unwrap();
        let lowered = BondDsl::from_ast(&raised, &cfg).unwrap();
        assert_eq!(input, lowered);
    }
}
