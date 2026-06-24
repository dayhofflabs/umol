//! Predicate sub-parsers and apply helpers shared between atom and bond DSL.

use std::{fmt, mem};

use winnow::ascii::{dec_uint, multispace0};
use winnow::combinator::{alt, delimited, empty, opt, preceded};
use winnow::error::ErrMode;
use winnow::Parser;

use super::config::{MultiplicityDefault, UnpairedElectronsDefault};
use super::error::{PResult, ParseError};
use super::value::{fmt_value, value};
use crate::ast::constraint::{RingMembershipAst, RingScope};
use crate::ast::spin::SpinStateAst;
use crate::ast::value::ValueAst;

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
            "+".value(ValueAst::RangeFrom(1)),
            "!".value(ValueAst::Lit(0)),
            empty.value(ValueAst::Lit(1)),
        )),
    )
    .parse_next(i)
    .map_err(|_: ErrMode<ParseError>| ErrMode::Backtrack(ParseError::ExpectedPredicateBody))
}

/// Parse the body after `#R`: an optional `(size)` then a count, yielding the
/// `RingScope` (`All` when no size) and its count.
pub(crate) fn ring_membership(i: &mut &str) -> PResult<RingMembershipAst> {
    let size: Option<u8> = opt(delimited('(', dec_uint::<_, u8, _>, ')'))
        .parse_next(i)
        .map_err(|_: ErrMode<ParseError>| ErrMode::Backtrack(ParseError::ExpectedPredicateBody))?;
    let count = ring_count(i)?;
    Ok(RingMembershipAst::new(
        size.map_or(RingScope::All, RingScope::Size),
        count,
    ))
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
    Ok(())
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

// Canonical-rendering rules:
//
// - Vacuous constraints (any constraint whose payload is `ValueAst::Undetermined`)
//   are elided from the rendered surface form. The AST may carry them; the
//   canonical entity / molecule string does not show them.
// - Inherent fields whose payload is `ValueAst::Undetermined` are likewise
//   elided when they have a leading prefix (`#c`, `#u`, `#s`, `#e`, …).
// - Exception: leading **unprefixed** fields — bond order, atom element,
//   noncovalent bond type — cannot be elided because they fix the entity
//   string's start position. For these, `Undetermined` renders as `*`.
//
// Consequence for round-trip: rendering a vacuous constraint and reparsing
// produces an AST without that constraint, so AST equality across a
// render/parse cycle requires either generating only non-vacuous payloads
// or normalizing the input before comparing.

pub(crate) fn fmt_ring_membership(
    f: &mut fmt::Formatter<'_>,
    m: &RingMembershipAst,
) -> fmt::Result {
    let v = &m.count;
    if matches!(v, ValueAst::Undetermined) {
        return Ok(());
    }
    write!(f, "#R")?;
    if let RingScope::Size(s) = m.scope {
        write!(f, "({})", s)?;
    }
    if *v == ValueAst::RangeFrom(1) {
        return write!(f, "+");
    }
    match v {
        ValueAst::Lit(0) => write!(f, "!"),
        ValueAst::Lit(1) => Ok(()),
        ValueAst::Lit(n) => write!(f, "{}", n),
        v => fmt_value(f, v),
    }
}

/// Fill defaults on a `SpinStateAst` per the given modes. Shared across
/// atom/bond/aromatic-system/multicenter-bond DSL lowering (all entities that
/// carry a spin state except `NoncovalentBond`).
pub(crate) fn raise_spin(
    spin: &mut SpinStateAst,
    u_mode: UnpairedElectronsDefault,
    m_mode: MultiplicityDefault,
) {
    let u = mem::replace(&mut spin.unpaired, ValueAst::Undetermined);
    let m = mem::replace(&mut spin.multiplicity, ValueAst::Undetermined);
    let resolved_u = if matches!(u, ValueAst::Undetermined) {
        match u_mode {
            UnpairedElectronsDefault::Zero => ValueAst::Lit(0),
            UnpairedElectronsDefault::Required => ValueAst::Undetermined,
            UnpairedElectronsDefault::Derived => match &m {
                ValueAst::Lit(mm) => ValueAst::Lit(mm - 1),
                _ => ValueAst::Undetermined,
            },
        }
    } else {
        u
    };
    let resolved_m = if matches!(m, ValueAst::Undetermined) {
        match m_mode {
            MultiplicityDefault::Required => ValueAst::Undetermined,
            MultiplicityDefault::Derived => match &resolved_u {
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

/// Strip defaults from a `SpinStateAst` per the given modes. Compute `strip_m`
/// first so that under (Derived, Derived) the tie-break keeps `u` explicit:
/// `strip_u` under `Derived` backs off when `strip_m` has already fired, so at
/// most one of the two is stripped and re-raising recovers the original AST.
pub(crate) fn lower_spin(
    spin: &mut SpinStateAst,
    u_mode: UnpairedElectronsDefault,
    m_mode: MultiplicityDefault,
) {
    let uu = if let ValueAst::Lit(n) = spin.unpaired {
        Some(n)
    } else {
        None
    };
    let mm = if let ValueAst::Lit(n) = spin.multiplicity {
        Some(n)
    } else {
        None
    };
    let derived = matches!((uu, mm), (Some(u), Some(m)) if m == u + 1);

    let strip_m = matches!(m_mode, MultiplicityDefault::Derived) && derived;
    let strip_u = match u_mode {
        UnpairedElectronsDefault::Zero => uu == Some(0),
        UnpairedElectronsDefault::Derived => derived && mm.is_some() && !strip_m,
        UnpairedElectronsDefault::Required => false,
    };
    if strip_u {
        spin.unpaired = ValueAst::Undetermined;
    }
    if strip_m {
        spin.multiplicity = ValueAst::Undetermined;
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[fixture]
    fn spin(
        #[default(ValueAst::Undetermined)] unpaired: ValueAst,
        #[default(ValueAst::Undetermined)] multiplicity: ValueAst,
    ) -> SpinStateAst {
        SpinStateAst {
            unpaired,
            multiplicity,
        }
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::plus("+", ValueAst::Lit(1))]
    #[case::minus("-", ValueAst::Lit(-1))]
    #[case::zero("0", ValueAst::Lit(0))]
    #[case::pos_lit("+2", ValueAst::Lit(2))]
    #[case::neg_lit("-3", ValueAst::Lit(-3))]
    #[case::undetermined("*", ValueAst::Undetermined)]
    #[case::set("{1,2,3}", ValueAst::lit_set([1, 2, 3]))]
    fn test_charge(#[case] input: &str, #[case] expected: ValueAst) {
        let result = charge.parse(input).unwrap();
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case::letters("abc")]
    fn test_charge_error(#[case] input: &str) {
        assert_eq!(
            charge.parse(input).unwrap_err().into_inner(),
            ParseError::ExpectedPredicateBody
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("", ValueAst::Lit(1))]
    #[case::lit("3", ValueAst::Lit(3))]
    #[case::zero("0", ValueAst::Lit(0))]
    #[case::neg("-5", ValueAst::Lit(-5))]
    #[case::undetermined("*", ValueAst::Undetermined)]
    #[case::set("{1,2}", ValueAst::lit_set([1, 2]))]
    fn test_optional_value(#[case] input: &str, #[case] expected: ValueAst) {
        let result = optional_value.parse(input).unwrap();
        assert_eq!(result, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("", ValueAst::Lit(1))]
    #[case::lit("4", ValueAst::Lit(4))]
    #[case::undetermined("*", ValueAst::Undetermined)]
    #[case::plus_sugar("+", ValueAst::RangeFrom(1))]
    #[case::bang_sugar("!", ValueAst::Lit(0))]
    #[case::zero_numeric("0", ValueAst::Lit(0))]
    #[case::set("{2,3}", ValueAst::lit_set([2, 3]))]
    fn test_ring_count(#[case] input: &str, #[case] expected: ValueAst) {
        let result = ring_count.parse(input).unwrap();
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case::all_zero_renders_bang(RingScope::All, ValueAst::Lit(0), "#R!")]
    #[case::all_one_renders_bare(RingScope::All, ValueAst::Lit(1), "#R")]
    #[case::all_two(RingScope::All, ValueAst::Lit(2), "#R2")]
    #[case::all_plus(RingScope::All, ValueAst::RangeFrom(1), "#R+")]
    #[case::size_bare(RingScope::Size(6), ValueAst::Lit(1), "#R(6)")]
    #[case::size_plus(RingScope::Size(6), ValueAst::RangeFrom(1), "#R(6)+")]
    fn test_fmt_ring_membership(
        #[case] scope: RingScope,
        #[case] count: ValueAst,
        #[case] expected: &str,
    ) {
        struct W(RingMembershipAst);
        impl fmt::Display for W {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                fmt_ring_membership(f, &self.0)
            }
        }
        assert_eq!(W(RingMembershipAst { scope, count }).to_string(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::sets_unpaired(spin(ValueAst::Undetermined, ValueAst::Undetermined), SpinPredicate::Unpaired(ValueAst::Lit(1)),
        spin(ValueAst::Lit(1), ValueAst::Undetermined))]
    #[case::sets_multiplicity(spin(ValueAst::Undetermined, ValueAst::Undetermined), SpinPredicate::Multiplicity(ValueAst::Lit(2)),
        spin(ValueAst::Undetermined, ValueAst::Lit(2)))]
    #[case::sets_unpaired_over_existing_multiplicity(spin(ValueAst::Undetermined, ValueAst::Lit(2)), SpinPredicate::Unpaired(ValueAst::Lit(1)),
        spin(ValueAst::Lit(1), ValueAst::Lit(2)))]
    #[case::sets_multiplicity_over_existing_unpaired(spin(ValueAst::Lit(0), ValueAst::Undetermined), SpinPredicate::Multiplicity(ValueAst::Lit(1)),
        spin(ValueAst::Lit(0), ValueAst::Lit(1)))]
    fn test_apply_spin_pair(
        #[case] mut initial: SpinStateAst,
        #[case] pred: SpinPredicate,
        #[case] expected: SpinStateAst,
    ) {
        apply_spin_pair(&mut initial, pred, ParseError::DuplicateAtomPredicate).unwrap();
        assert_eq!(initial, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::duplicate_unpaired(spin(ValueAst::Lit(1), ValueAst::Undetermined), SpinPredicate::Unpaired(ValueAst::Lit(2)),
        ParseError::DuplicateAtomPredicate("#u".to_string()))]
    #[case::duplicate_multiplicity(spin(ValueAst::Undetermined, ValueAst::Lit(2)), SpinPredicate::Multiplicity(ValueAst::Lit(3)),
        ParseError::DuplicateAtomPredicate("#s".to_string()))]
    fn test_apply_spin_pair_error(
        #[case] mut initial: SpinStateAst,
        #[case] pred: SpinPredicate,
        #[case] expected: ParseError,
    ) {
        let err = apply_spin_pair(&mut initial, pred, ParseError::DuplicateAtomPredicate)
            .unwrap_err();
        assert_eq!(err, expected);
    }

    /// For every (initial u, m) pair reachable after parsing the seven canonical DSL fragments (`""`, `"#u0"`, `"#u1"`, `"#m1"`, `"#m2"`,
    /// `"#u1#m1"`, `"#u1#m2"`) under each of the six mode combinations, the raised `(u, m)` must match the expected value.
    #[rustfmt::skip]
    #[rstest]
    // u: Zero, m: Derived
    #[case::zd_empty(ValueAst::Undetermined, ValueAst::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived, ValueAst::Lit(0), ValueAst::Lit(1))]
    #[case::zd_u0(ValueAst::Lit(0), ValueAst::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived, ValueAst::Lit(0), ValueAst::Lit(1))]
    #[case::zd_u1(ValueAst::Lit(1), ValueAst::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived, ValueAst::Lit(1), ValueAst::Lit(2))]
    #[case::zd_m1(ValueAst::Undetermined, ValueAst::Lit(1), UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived, ValueAst::Lit(0), ValueAst::Lit(1))]
    #[case::zd_m2(ValueAst::Undetermined, ValueAst::Lit(2), UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived, ValueAst::Lit(0), ValueAst::Lit(2))]
    #[case::zd_u1m1(ValueAst::Lit(1), ValueAst::Lit(1), UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived, ValueAst::Lit(1), ValueAst::Lit(1))]
    #[case::zd_u1m2(ValueAst::Lit(1), ValueAst::Lit(2), UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived, ValueAst::Lit(1), ValueAst::Lit(2))]
    // u: Zero, m: Required
    #[case::zr_empty(ValueAst::Undetermined, ValueAst::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, ValueAst::Lit(0), ValueAst::Undetermined)]
    #[case::zr_u0(ValueAst::Lit(0), ValueAst::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, ValueAst::Lit(0), ValueAst::Undetermined)]
    #[case::zr_u1(ValueAst::Lit(1), ValueAst::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, ValueAst::Lit(1), ValueAst::Undetermined)]
    #[case::zr_m1(ValueAst::Undetermined, ValueAst::Lit(1), UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, ValueAst::Lit(0), ValueAst::Lit(1))]
    #[case::zr_m2(ValueAst::Undetermined, ValueAst::Lit(2), UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, ValueAst::Lit(0), ValueAst::Lit(2))]
    #[case::zr_u1m1(ValueAst::Lit(1), ValueAst::Lit(1), UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, ValueAst::Lit(1), ValueAst::Lit(1))]
    #[case::zr_u1m2(ValueAst::Lit(1), ValueAst::Lit(2), UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, ValueAst::Lit(1), ValueAst::Lit(2))]
    // u: Required, m: Derived
    #[case::rd_empty(ValueAst::Undetermined, ValueAst::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, ValueAst::Undetermined, ValueAst::Undetermined)]
    #[case::rd_u0(ValueAst::Lit(0), ValueAst::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, ValueAst::Lit(0), ValueAst::Lit(1))]
    #[case::rd_u1(ValueAst::Lit(1), ValueAst::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, ValueAst::Lit(1), ValueAst::Lit(2))]
    #[case::rd_m1(ValueAst::Undetermined, ValueAst::Lit(1), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, ValueAst::Undetermined, ValueAst::Lit(1))]
    #[case::rd_m2(ValueAst::Undetermined, ValueAst::Lit(2), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, ValueAst::Undetermined, ValueAst::Lit(2))]
    #[case::rd_u1m1(ValueAst::Lit(1), ValueAst::Lit(1), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, ValueAst::Lit(1), ValueAst::Lit(1))]
    #[case::rd_u1m2(ValueAst::Lit(1), ValueAst::Lit(2), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, ValueAst::Lit(1), ValueAst::Lit(2))]
    // u: Required, m: Required
    #[case::rr_empty(ValueAst::Undetermined, ValueAst::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Required, ValueAst::Undetermined, ValueAst::Undetermined)]
    #[case::rr_u0(ValueAst::Lit(0), ValueAst::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Required, ValueAst::Lit(0), ValueAst::Undetermined)]
    #[case::rr_u1(ValueAst::Lit(1), ValueAst::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Required, ValueAst::Lit(1), ValueAst::Undetermined)]
    #[case::rr_m1(ValueAst::Undetermined, ValueAst::Lit(1), UnpairedElectronsDefault::Required, MultiplicityDefault::Required, ValueAst::Undetermined, ValueAst::Lit(1))]
    #[case::rr_m2(ValueAst::Undetermined, ValueAst::Lit(2), UnpairedElectronsDefault::Required, MultiplicityDefault::Required, ValueAst::Undetermined, ValueAst::Lit(2))]
    #[case::rr_u1m1(ValueAst::Lit(1), ValueAst::Lit(1), UnpairedElectronsDefault::Required, MultiplicityDefault::Required, ValueAst::Lit(1), ValueAst::Lit(1))]
    #[case::rr_u1m2(ValueAst::Lit(1), ValueAst::Lit(2), UnpairedElectronsDefault::Required, MultiplicityDefault::Required, ValueAst::Lit(1), ValueAst::Lit(2))]
    // u: Derived, m: Required
    #[case::dr_empty(ValueAst::Undetermined, ValueAst::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, ValueAst::Undetermined, ValueAst::Undetermined)]
    #[case::dr_u0(ValueAst::Lit(0), ValueAst::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, ValueAst::Lit(0), ValueAst::Undetermined)]
    #[case::dr_u1(ValueAst::Lit(1), ValueAst::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, ValueAst::Lit(1), ValueAst::Undetermined)]
    #[case::dr_m1(ValueAst::Undetermined, ValueAst::Lit(1), UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, ValueAst::Lit(0), ValueAst::Lit(1))]
    #[case::dr_m2(ValueAst::Undetermined, ValueAst::Lit(2), UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, ValueAst::Lit(1), ValueAst::Lit(2))]
    #[case::dr_u1m1(ValueAst::Lit(1), ValueAst::Lit(1), UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, ValueAst::Lit(1), ValueAst::Lit(1))]
    #[case::dr_u1m2(ValueAst::Lit(1), ValueAst::Lit(2), UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, ValueAst::Lit(1), ValueAst::Lit(2))]
    // u: Derived, m: Derived
    #[case::dd_empty(ValueAst::Undetermined, ValueAst::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived, ValueAst::Undetermined, ValueAst::Undetermined)]
    #[case::dd_u0(ValueAst::Lit(0), ValueAst::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived, ValueAst::Lit(0), ValueAst::Lit(1))]
    #[case::dd_u1(ValueAst::Lit(1), ValueAst::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived, ValueAst::Lit(1), ValueAst::Lit(2))]
    #[case::dd_m1(ValueAst::Undetermined, ValueAst::Lit(1), UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived, ValueAst::Lit(0), ValueAst::Lit(1))]
    #[case::dd_m2(ValueAst::Undetermined, ValueAst::Lit(2), UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived, ValueAst::Lit(1), ValueAst::Lit(2))]
    #[case::dd_u1m1(ValueAst::Lit(1), ValueAst::Lit(1), UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived, ValueAst::Lit(1), ValueAst::Lit(1))]
    #[case::dd_u1m2(ValueAst::Lit(1), ValueAst::Lit(2), UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived, ValueAst::Lit(1), ValueAst::Lit(2))]
    fn test_raise_spin(
        #[case] init_u: ValueAst,
        #[case] init_m: ValueAst,
        #[case] u_mode: UnpairedElectronsDefault,
        #[case] m_mode: MultiplicityDefault,
        #[case] expected_u: ValueAst,
        #[case] expected_m: ValueAst,
    ) {
        let mut spin = SpinStateAst { unpaired: init_u, multiplicity: init_m };
        raise_spin(&mut spin, u_mode, m_mode);
        assert_eq!(spin.unpaired, expected_u);
        assert_eq!(spin.multiplicity, expected_m);
    }

    /// Per-mode lowering: covers the AST states reachable by raising the canonical DSL fragments, plus (U, U) where applicable.
    #[rustfmt::skip]
    #[rstest]
    // u: Zero, m: Derived
    #[case::zd_derived_zero(ValueAst::Lit(0), ValueAst::Lit(1), UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived, ValueAst::Undetermined, ValueAst::Undetermined)]
    #[case::zd_derived_nonzero(ValueAst::Lit(1), ValueAst::Lit(2), UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived, ValueAst::Lit(1), ValueAst::Undetermined)]
    #[case::zd_zero_nonderived(ValueAst::Lit(0), ValueAst::Lit(2), UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived, ValueAst::Undetermined, ValueAst::Lit(2))]
    #[case::zd_nonzero_nonderived(ValueAst::Lit(1), ValueAst::Lit(1), UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived, ValueAst::Lit(1), ValueAst::Lit(1))]
    // u: Zero, m: Required
    #[case::zr_zero_mundef(ValueAst::Lit(0), ValueAst::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, ValueAst::Undetermined, ValueAst::Undetermined)]
    #[case::zr_nonzero_mundef(ValueAst::Lit(1), ValueAst::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, ValueAst::Lit(1), ValueAst::Undetermined)]
    #[case::zr_zero_derived(ValueAst::Lit(0), ValueAst::Lit(1), UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, ValueAst::Undetermined, ValueAst::Lit(1))]
    #[case::zr_zero_doublet(ValueAst::Lit(0), ValueAst::Lit(2), UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, ValueAst::Undetermined, ValueAst::Lit(2))]
    #[case::zr_nonzero_singlet(ValueAst::Lit(1), ValueAst::Lit(1), UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, ValueAst::Lit(1), ValueAst::Lit(1))]
    #[case::zr_nonzero_doublet(ValueAst::Lit(1), ValueAst::Lit(2), UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, ValueAst::Lit(1), ValueAst::Lit(2))]
    // u: Required, m: Derived
    #[case::rd_both_undef(ValueAst::Undetermined, ValueAst::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, ValueAst::Undetermined, ValueAst::Undetermined)]
    #[case::rd_derived_zero(ValueAst::Lit(0), ValueAst::Lit(1), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, ValueAst::Lit(0), ValueAst::Undetermined)]
    #[case::rd_derived_nonzero(ValueAst::Lit(1), ValueAst::Lit(2), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, ValueAst::Lit(1), ValueAst::Undetermined)]
    #[case::rd_uundef_msinglet(ValueAst::Undetermined, ValueAst::Lit(1), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, ValueAst::Undetermined, ValueAst::Lit(1))]
    #[case::rd_uundef_mdoublet(ValueAst::Undetermined, ValueAst::Lit(2), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, ValueAst::Undetermined, ValueAst::Lit(2))]
    #[case::rd_nonderived(ValueAst::Lit(1), ValueAst::Lit(1), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, ValueAst::Lit(1), ValueAst::Lit(1))]
    // u: Required, m: Required — nothing strips.
    #[case::rr_both_undef(ValueAst::Undetermined, ValueAst::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Required, ValueAst::Undetermined, ValueAst::Undetermined)]
    #[case::rr_full(ValueAst::Lit(0), ValueAst::Lit(1), UnpairedElectronsDefault::Required, MultiplicityDefault::Required, ValueAst::Lit(0), ValueAst::Lit(1))]
    // u: Derived, m: Required
    #[case::dr_both_undef(ValueAst::Undetermined, ValueAst::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, ValueAst::Undetermined, ValueAst::Undetermined)]
    #[case::dr_zero_mundef(ValueAst::Lit(0), ValueAst::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, ValueAst::Lit(0), ValueAst::Undetermined)]
    #[case::dr_nonzero_mundef(ValueAst::Lit(1), ValueAst::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, ValueAst::Lit(1), ValueAst::Undetermined)]
    #[case::dr_derived_zero(ValueAst::Lit(0), ValueAst::Lit(1), UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, ValueAst::Undetermined, ValueAst::Lit(1))]
    #[case::dr_derived_nonzero(ValueAst::Lit(1), ValueAst::Lit(2), UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, ValueAst::Undetermined, ValueAst::Lit(2))]
    #[case::dr_nonderived(ValueAst::Lit(1), ValueAst::Lit(1), UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, ValueAst::Lit(1), ValueAst::Lit(1))]
    // u: Derived, m: Derived — tie-break keeps u explicit.
    #[case::dd_both_undef(ValueAst::Undetermined, ValueAst::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived, ValueAst::Undetermined, ValueAst::Undetermined)]
    #[case::dd_derived_zero(ValueAst::Lit(0), ValueAst::Lit(1), UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived, ValueAst::Lit(0), ValueAst::Undetermined)]
    #[case::dd_derived_nonzero(ValueAst::Lit(1), ValueAst::Lit(2), UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived, ValueAst::Lit(1), ValueAst::Undetermined)]
    #[case::dd_nonderived(ValueAst::Lit(1), ValueAst::Lit(1), UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived, ValueAst::Lit(1), ValueAst::Lit(1))]
    fn test_lower_spin(
        #[case] init_u: ValueAst,
        #[case] init_m: ValueAst,
        #[case] u_mode: UnpairedElectronsDefault,
        #[case] m_mode: MultiplicityDefault,
        #[case] expected_u: ValueAst,
        #[case] expected_m: ValueAst,
    ) {
        let mut spin = SpinStateAst { unpaired: init_u, multiplicity: init_m };
        lower_spin(&mut spin, u_mode, m_mode);
        assert_eq!(spin.unpaired, expected_u);
        assert_eq!(spin.multiplicity, expected_m);
    }

    /// AST preservation: the raised AST is a fixed point of `lower → raise`. Lowering strips default content; re-raising the result must recover the same AST.
    #[rustfmt::skip]
    #[rstest]
    #[case::zd_empty(ValueAst::Undetermined, ValueAst::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived)]
    #[case::zd_u0(ValueAst::Lit(0), ValueAst::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived)]
    #[case::zd_u1(ValueAst::Lit(1), ValueAst::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived)]
    #[case::zd_m1(ValueAst::Undetermined, ValueAst::Lit(1), UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived)]
    #[case::zd_m2(ValueAst::Undetermined, ValueAst::Lit(2), UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived)]
    #[case::zd_u1m1(ValueAst::Lit(1), ValueAst::Lit(1), UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived)]
    #[case::zd_u1m2(ValueAst::Lit(1), ValueAst::Lit(2), UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived)]
    #[case::zr_empty(ValueAst::Undetermined, ValueAst::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Required)]
    #[case::zr_u0(ValueAst::Lit(0), ValueAst::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Required)]
    #[case::zr_u1(ValueAst::Lit(1), ValueAst::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Required)]
    #[case::zr_m1(ValueAst::Undetermined, ValueAst::Lit(1), UnpairedElectronsDefault::Zero, MultiplicityDefault::Required)]
    #[case::zr_m2(ValueAst::Undetermined, ValueAst::Lit(2), UnpairedElectronsDefault::Zero, MultiplicityDefault::Required)]
    #[case::zr_u1m1(ValueAst::Lit(1), ValueAst::Lit(1), UnpairedElectronsDefault::Zero, MultiplicityDefault::Required)]
    #[case::zr_u1m2(ValueAst::Lit(1), ValueAst::Lit(2), UnpairedElectronsDefault::Zero, MultiplicityDefault::Required)]
    #[case::rd_empty(ValueAst::Undetermined, ValueAst::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Derived)]
    #[case::rd_u0(ValueAst::Lit(0), ValueAst::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Derived)]
    #[case::rd_u1(ValueAst::Lit(1), ValueAst::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Derived)]
    #[case::rd_m1(ValueAst::Undetermined, ValueAst::Lit(1), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived)]
    #[case::rd_m2(ValueAst::Undetermined, ValueAst::Lit(2), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived)]
    #[case::rd_u1m1(ValueAst::Lit(1), ValueAst::Lit(1), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived)]
    #[case::rd_u1m2(ValueAst::Lit(1), ValueAst::Lit(2), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived)]
    #[case::rr_empty(ValueAst::Undetermined, ValueAst::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Required)]
    #[case::rr_u0(ValueAst::Lit(0), ValueAst::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Required)]
    #[case::rr_u1(ValueAst::Lit(1), ValueAst::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Required)]
    #[case::rr_m1(ValueAst::Undetermined, ValueAst::Lit(1), UnpairedElectronsDefault::Required, MultiplicityDefault::Required)]
    #[case::rr_m2(ValueAst::Undetermined, ValueAst::Lit(2), UnpairedElectronsDefault::Required, MultiplicityDefault::Required)]
    #[case::rr_u1m1(ValueAst::Lit(1), ValueAst::Lit(1), UnpairedElectronsDefault::Required, MultiplicityDefault::Required)]
    #[case::rr_u1m2(ValueAst::Lit(1), ValueAst::Lit(2), UnpairedElectronsDefault::Required, MultiplicityDefault::Required)]
    #[case::dr_empty(ValueAst::Undetermined, ValueAst::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Required)]
    #[case::dr_u0(ValueAst::Lit(0), ValueAst::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Required)]
    #[case::dr_u1(ValueAst::Lit(1), ValueAst::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Required)]
    #[case::dr_m1(ValueAst::Undetermined, ValueAst::Lit(1), UnpairedElectronsDefault::Derived, MultiplicityDefault::Required)]
    #[case::dr_m2(ValueAst::Undetermined, ValueAst::Lit(2), UnpairedElectronsDefault::Derived, MultiplicityDefault::Required)]
    #[case::dr_u1m1(ValueAst::Lit(1), ValueAst::Lit(1), UnpairedElectronsDefault::Derived, MultiplicityDefault::Required)]
    #[case::dr_u1m2(ValueAst::Lit(1), ValueAst::Lit(2), UnpairedElectronsDefault::Derived, MultiplicityDefault::Required)]
    #[case::dd_empty(ValueAst::Undetermined, ValueAst::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived)]
    #[case::dd_u0(ValueAst::Lit(0), ValueAst::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived)]
    #[case::dd_u1(ValueAst::Lit(1), ValueAst::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived)]
    #[case::dd_m1(ValueAst::Undetermined, ValueAst::Lit(1), UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived)]
    #[case::dd_m2(ValueAst::Undetermined, ValueAst::Lit(2), UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived)]
    #[case::dd_u1m1(ValueAst::Lit(1), ValueAst::Lit(1), UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived)]
    #[case::dd_u1m2(ValueAst::Lit(1), ValueAst::Lit(2), UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived)]
    fn test_spin_roundtrip_preserves_ast(
        #[case] init_u: ValueAst,
        #[case] init_m: ValueAst,
        #[case] u_mode: UnpairedElectronsDefault,
        #[case] m_mode: MultiplicityDefault,
    ) {
        let mut raised = SpinStateAst { unpaired: init_u, multiplicity: init_m };
        raise_spin(&mut raised, u_mode, m_mode);

        let mut lowered_then_raised = raised.clone();
        lower_spin(&mut lowered_then_raised, u_mode, m_mode);
        raise_spin(&mut lowered_then_raised, u_mode, m_mode);

        assert_eq!(lowered_then_raised, raised);
    }
}
