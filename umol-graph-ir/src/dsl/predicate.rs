//! Predicate sub-parsers and apply helpers shared across the entity-string DSLs.

use std::{fmt, mem};

use winnow::ascii::{dec_uint, multispace0};
use winnow::combinator::{alt, delimited, empty, opt, preceded};
use winnow::error::ErrMode;
use winnow::Parser;

use super::config::{MultiplicityDefault, UnpairedElectronsDefault};
use super::error::{PResult, ParseError};
use super::value::{fmt_value, value};
use crate::ir::constraint::{RingMembershipAst, RingScope};
use crate::ir::spin::UnpairedElectronsForm;
use crate::ir::value::NumForm;

pub(crate) fn charge(i: &mut &str) -> PResult<NumForm> {
    preceded(
        multispace0,
        alt((
            value,
            "+".value(NumForm::Lit(1)),
            "-".value(NumForm::Lit(-1)),
        )),
    )
    .parse_next(i)
    .map_err(|_: ErrMode<ParseError>| ErrMode::Backtrack(ParseError::ExpectedPredicateBody))
}

pub(crate) fn optional_value(i: &mut &str) -> PResult<NumForm> {
    preceded(multispace0, alt((value, empty.value(NumForm::Lit(1)))))
        .parse_next(i)
        .map_err(|_: ErrMode<ParseError>| ErrMode::Backtrack(ParseError::ExpectedPredicateBody))
}

pub(crate) fn ring_count(i: &mut &str) -> PResult<NumForm> {
    preceded(
        multispace0,
        alt((
            value,
            "+".value(NumForm::RangeFrom(1)),
            "!".value(NumForm::Lit(0)),
            empty.value(NumForm::Lit(1)),
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
pub enum UnpairedElectronsPredicate {
    Count(NumForm),
    Multiplicity(NumForm),
}

pub(crate) fn apply_unpaired_electrons_predicate(
    unpaired_electrons: &mut UnpairedElectronsForm,
    predicate: UnpairedElectronsPredicate,
    duplicate_error: fn(String) -> ParseError,
) -> Result<(), ParseError> {
    match predicate {
        UnpairedElectronsPredicate::Count(v) => {
            if !matches!(&unpaired_electrons.count, NumForm::Undetermined) {
                return Err(duplicate_error("#u".to_string()));
            }
            unpaired_electrons.count = v;
        }
        UnpairedElectronsPredicate::Multiplicity(v) => {
            if !matches!(&unpaired_electrons.multiplicity, NumForm::Undetermined) {
                return Err(duplicate_error("#s".to_string()));
            }
            unpaired_electrons.multiplicity = v;
        }
    }
    Ok(())
}

pub(crate) fn fmt_charge(f: &mut fmt::Formatter<'_>, v: &NumForm) -> fmt::Result {
    match v {
        NumForm::Undetermined => Ok(()),
        NumForm::Lit(0) => write!(f, "#c0"),
        NumForm::Lit(1) => write!(f, "#c+"),
        NumForm::Lit(-1) => write!(f, "#c-"),
        NumForm::Lit(n) if *n > 0 => write!(f, "#c+{}", n),
        NumForm::Lit(n) => write!(f, "#c{}", n),
        v => {
            write!(f, "#c")?;
            fmt_value(f, v)
        }
    }
}

pub(crate) fn fmt_unpaired_electrons(
    f: &mut fmt::Formatter<'_>,
    unpaired_electrons: &UnpairedElectronsForm,
) -> fmt::Result {
    match &unpaired_electrons.count {
        NumForm::Undetermined => {}
        NumForm::Lit(1) => write!(f, "#u")?,
        NumForm::Lit(n) => write!(f, "#u{}", n)?,
        v => {
            write!(f, "#u")?;
            fmt_value(f, v)?;
        }
    }
    match &unpaired_electrons.multiplicity {
        NumForm::Undetermined => {}
        NumForm::Lit(1) => write!(f, "#s")?,
        NumForm::Lit(n) => write!(f, "#s{}", n)?,
        v => {
            write!(f, "#s")?;
            fmt_value(f, v)?;
        }
    }
    Ok(())
}

// Canonical-rendering rules:
//
// - Vacuous constraints (any constraint whose payload is `NumForm::Undetermined`)
//   are elided from the rendered surface form. The AST may carry them; the
//   canonical entity / molecule string does not show them.
// - Inherent fields whose payload is `NumForm::Undetermined` are likewise
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
    if matches!(v, NumForm::Undetermined) {
        return Ok(());
    }
    write!(f, "#R")?;
    if let RingScope::Size(s) = m.scope {
        write!(f, "({})", s)?;
    }
    if *v == NumForm::RangeFrom(1) {
        return write!(f, "+");
    }
    match v {
        NumForm::Lit(0) => write!(f, "!"),
        NumForm::Lit(1) => Ok(()),
        NumForm::Lit(n) => write!(f, "{}", n),
        v => fmt_value(f, v),
    }
}

/// Fill defaults on a `UnpairedElectronsForm` per the given modes. Shared across
/// atom/bond/aromatic-system/multicenter-bond DSL lowering (all entities that
/// carry an unpaired-electron state except `NoncovalentBond`).
pub(crate) fn raise_unpaired_electrons(
    unpaired_electrons: &mut UnpairedElectronsForm,
    count_default: UnpairedElectronsDefault,
    multiplicity_default: MultiplicityDefault,
) {
    let count = mem::replace(&mut unpaired_electrons.count, NumForm::Undetermined);
    let multiplicity = mem::replace(&mut unpaired_electrons.multiplicity, NumForm::Undetermined);
    let resolved_count = if matches!(count, NumForm::Undetermined) {
        match count_default {
            UnpairedElectronsDefault::Zero => NumForm::Lit(0),
            UnpairedElectronsDefault::Required => NumForm::Undetermined,
            UnpairedElectronsDefault::Derived => match &multiplicity {
                NumForm::Lit(value) => NumForm::Lit(value - 1),
                _ => NumForm::Undetermined,
            },
        }
    } else {
        count
    };
    let resolved_multiplicity = if matches!(multiplicity, NumForm::Undetermined) {
        match multiplicity_default {
            MultiplicityDefault::Required => NumForm::Undetermined,
            MultiplicityDefault::Derived => match &resolved_count {
                NumForm::Lit(value) => NumForm::Lit(value + 1),
                _ => NumForm::Undetermined,
            },
        }
    } else {
        multiplicity
    };
    unpaired_electrons.count = resolved_count;
    unpaired_electrons.multiplicity = resolved_multiplicity;
}

/// Strip defaults from a `UnpairedElectronsForm` per the given modes. Compute
/// `strip_multiplicity`
/// first so that under (Derived, Derived) the tie-break keeps `u` explicit:
/// `strip_count` under `Derived` backs off when `strip_multiplicity` has already fired, so at
/// most one of the two is stripped and re-raising recovers the original AST.
pub(crate) fn lower_unpaired_electrons(
    unpaired_electrons: &mut UnpairedElectronsForm,
    count_default: UnpairedElectronsDefault,
    multiplicity_default: MultiplicityDefault,
) {
    let literal_count = if let NumForm::Lit(n) = unpaired_electrons.count {
        Some(n)
    } else {
        None
    };
    let literal_multiplicity = if let NumForm::Lit(n) = unpaired_electrons.multiplicity {
        Some(n)
    } else {
        None
    };
    let derived = matches!(
        (literal_count, literal_multiplicity),
        (Some(count), Some(multiplicity)) if multiplicity == count + 1
    );

    let strip_multiplicity =
        matches!(multiplicity_default, MultiplicityDefault::Derived) && derived;
    let strip_count = match count_default {
        UnpairedElectronsDefault::Zero => literal_count == Some(0),
        UnpairedElectronsDefault::Derived => {
            derived && literal_multiplicity.is_some() && !strip_multiplicity
        }
        UnpairedElectronsDefault::Required => false,
    };
    if strip_count {
        unpaired_electrons.count = NumForm::Undetermined;
    }
    if strip_multiplicity {
        unpaired_electrons.multiplicity = NumForm::Undetermined;
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::plus("+", NumForm::Lit(1))]
    #[case::minus("-", NumForm::Lit(-1))]
    #[case::zero("0", NumForm::Lit(0))]
    #[case::pos_lit("+2", NumForm::Lit(2))]
    #[case::neg_lit("-3", NumForm::Lit(-3))]
    #[case::undetermined("*", NumForm::Undetermined)]
    #[case::set("{1,2,3}", NumForm::lit_set([1, 2, 3]))]
    fn test_charge(#[case] input: &str, #[case] expected: NumForm) {
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
    #[case::empty("", NumForm::Lit(1))]
    #[case::lit("3", NumForm::Lit(3))]
    #[case::zero("0", NumForm::Lit(0))]
    #[case::neg("-5", NumForm::Lit(-5))]
    #[case::undetermined("*", NumForm::Undetermined)]
    #[case::set("{1,2}", NumForm::lit_set([1, 2]))]
    fn test_optional_value(#[case] input: &str, #[case] expected: NumForm) {
        let result = optional_value.parse(input).unwrap();
        assert_eq!(result, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("", NumForm::Lit(1))]
    #[case::lit("4", NumForm::Lit(4))]
    #[case::undetermined("*", NumForm::Undetermined)]
    #[case::plus_sugar("+", NumForm::RangeFrom(1))]
    #[case::bang_sugar("!", NumForm::Lit(0))]
    #[case::zero_numeric("0", NumForm::Lit(0))]
    #[case::set("{2,3}", NumForm::lit_set([2, 3]))]
    fn test_ring_count(#[case] input: &str, #[case] expected: NumForm) {
        let result = ring_count.parse(input).unwrap();
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case::all_zero_renders_bang(RingScope::All, NumForm::Lit(0), "#R!")]
    #[case::all_one_renders_bare(RingScope::All, NumForm::Lit(1), "#R")]
    #[case::all_two(RingScope::All, NumForm::Lit(2), "#R2")]
    #[case::all_plus(RingScope::All, NumForm::RangeFrom(1), "#R+")]
    #[case::size_bare(RingScope::Size(6), NumForm::Lit(1), "#R(6)")]
    #[case::size_plus(RingScope::Size(6), NumForm::RangeFrom(1), "#R(6)+")]
    fn test_fmt_ring_membership(
        #[case] scope: RingScope,
        #[case] count: NumForm,
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
    #[case::count(UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Undetermined }, UnpairedElectronsPredicate::Count(NumForm::Lit(1)),
        UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Undetermined })]
    #[case::multiplicity(UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Undetermined }, UnpairedElectronsPredicate::Multiplicity(NumForm::Lit(2)),
        UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(2) })]
    #[case::count_with_multiplicity(UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(2) }, UnpairedElectronsPredicate::Count(NumForm::Lit(1)),
        UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) })]
    #[case::multiplicity_with_count(UnpairedElectronsForm { count: NumForm::Lit(0), multiplicity: NumForm::Undetermined }, UnpairedElectronsPredicate::Multiplicity(NumForm::Lit(1)),
        UnpairedElectronsForm { count: NumForm::Lit(0), multiplicity: NumForm::Lit(1) })]
    fn test_apply_unpaired_electrons_predicate(
        #[case] mut initial: UnpairedElectronsForm,
        #[case] predicate: UnpairedElectronsPredicate,
        #[case] expected: UnpairedElectronsForm,
    ) {
        apply_unpaired_electrons_predicate(
            &mut initial,
            predicate,
            ParseError::DuplicateAtomPredicate,
        )
        .unwrap();
        assert_eq!(initial, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::count(UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Undetermined }, UnpairedElectronsPredicate::Count(NumForm::Lit(2)),
        ParseError::DuplicateAtomPredicate("#u".to_string()))]
    #[case::multiplicity(UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(2) }, UnpairedElectronsPredicate::Multiplicity(NumForm::Lit(3)),
        ParseError::DuplicateAtomPredicate("#s".to_string()))]
    fn test_apply_unpaired_electrons_predicate_error(
        #[case] mut initial: UnpairedElectronsForm,
        #[case] predicate: UnpairedElectronsPredicate,
        #[case] expected: ParseError,
    ) {
        let err = apply_unpaired_electrons_predicate(
            &mut initial,
            predicate,
            ParseError::DuplicateAtomPredicate,
        )
        .unwrap_err();
        assert_eq!(err, expected);
    }

    #[rstest]
    #[case::undetermined(UnpairedElectronsForm::default(), "")]
    #[case::count_one(UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Undetermined }, "#u")]
    #[case::count(UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Undetermined }, "#u2")]
    #[case::multiplicity_one(UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(1) }, "#s")]
    #[case::multiplicity(UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(3) }, "#s3")]
    #[case::complete(UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Lit(1) }, "#u2#s")]
    fn test_fmt_unpaired_electrons(
        #[case] unpaired_electrons: UnpairedElectronsForm,
        #[case] expected: &str,
    ) {
        struct DisplayUnpairedElectrons(UnpairedElectronsForm);

        impl fmt::Display for DisplayUnpairedElectrons {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                fmt_unpaired_electrons(f, &self.0)
            }
        }

        assert_eq!(
            DisplayUnpairedElectrons(unpaired_electrons).to_string(),
            expected
        );
    }

    /// For every (initial u, s) pair reachable after parsing the seven canonical DSL fragments (`""`, `"#u0"`, `"#u1"`, `"#s1"`, `"#s2"`,
    /// `"#u1#s1"`, `"#u1#s2"`) under each of the six mode combinations, the raised `(u, s)` must match the expected value.
    #[rustfmt::skip]
    #[rstest]
    // u: Zero, m: Derived
    #[case::zd_empty(NumForm::Undetermined, NumForm::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived, NumForm::Lit(0), NumForm::Lit(1))]
    #[case::zd_u0(NumForm::Lit(0), NumForm::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived, NumForm::Lit(0), NumForm::Lit(1))]
    #[case::zd_u1(NumForm::Lit(1), NumForm::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived, NumForm::Lit(1), NumForm::Lit(2))]
    #[case::zd_m1(NumForm::Undetermined, NumForm::Lit(1), UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived, NumForm::Lit(0), NumForm::Lit(1))]
    #[case::zd_m2(NumForm::Undetermined, NumForm::Lit(2), UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived, NumForm::Lit(0), NumForm::Lit(2))]
    #[case::zd_u1m1(NumForm::Lit(1), NumForm::Lit(1), UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived, NumForm::Lit(1), NumForm::Lit(1))]
    #[case::zd_u1m2(NumForm::Lit(1), NumForm::Lit(2), UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived, NumForm::Lit(1), NumForm::Lit(2))]
    // u: Zero, m: Required
    #[case::zr_empty(NumForm::Undetermined, NumForm::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, NumForm::Lit(0), NumForm::Undetermined)]
    #[case::zr_u0(NumForm::Lit(0), NumForm::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, NumForm::Lit(0), NumForm::Undetermined)]
    #[case::zr_u1(NumForm::Lit(1), NumForm::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, NumForm::Lit(1), NumForm::Undetermined)]
    #[case::zr_m1(NumForm::Undetermined, NumForm::Lit(1), UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, NumForm::Lit(0), NumForm::Lit(1))]
    #[case::zr_m2(NumForm::Undetermined, NumForm::Lit(2), UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, NumForm::Lit(0), NumForm::Lit(2))]
    #[case::zr_u1m1(NumForm::Lit(1), NumForm::Lit(1), UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, NumForm::Lit(1), NumForm::Lit(1))]
    #[case::zr_u1m2(NumForm::Lit(1), NumForm::Lit(2), UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, NumForm::Lit(1), NumForm::Lit(2))]
    // u: Required, m: Derived
    #[case::rd_empty(NumForm::Undetermined, NumForm::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, NumForm::Undetermined, NumForm::Undetermined)]
    #[case::rd_u0(NumForm::Lit(0), NumForm::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, NumForm::Lit(0), NumForm::Lit(1))]
    #[case::rd_u1(NumForm::Lit(1), NumForm::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, NumForm::Lit(1), NumForm::Lit(2))]
    #[case::rd_m1(NumForm::Undetermined, NumForm::Lit(1), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, NumForm::Undetermined, NumForm::Lit(1))]
    #[case::rd_m2(NumForm::Undetermined, NumForm::Lit(2), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, NumForm::Undetermined, NumForm::Lit(2))]
    #[case::rd_u1m1(NumForm::Lit(1), NumForm::Lit(1), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, NumForm::Lit(1), NumForm::Lit(1))]
    #[case::rd_u1m2(NumForm::Lit(1), NumForm::Lit(2), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, NumForm::Lit(1), NumForm::Lit(2))]
    // u: Required, m: Required
    #[case::rr_empty(NumForm::Undetermined, NumForm::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Required, NumForm::Undetermined, NumForm::Undetermined)]
    #[case::rr_u0(NumForm::Lit(0), NumForm::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Required, NumForm::Lit(0), NumForm::Undetermined)]
    #[case::rr_u1(NumForm::Lit(1), NumForm::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Required, NumForm::Lit(1), NumForm::Undetermined)]
    #[case::rr_m1(NumForm::Undetermined, NumForm::Lit(1), UnpairedElectronsDefault::Required, MultiplicityDefault::Required, NumForm::Undetermined, NumForm::Lit(1))]
    #[case::rr_m2(NumForm::Undetermined, NumForm::Lit(2), UnpairedElectronsDefault::Required, MultiplicityDefault::Required, NumForm::Undetermined, NumForm::Lit(2))]
    #[case::rr_u1m1(NumForm::Lit(1), NumForm::Lit(1), UnpairedElectronsDefault::Required, MultiplicityDefault::Required, NumForm::Lit(1), NumForm::Lit(1))]
    #[case::rr_u1m2(NumForm::Lit(1), NumForm::Lit(2), UnpairedElectronsDefault::Required, MultiplicityDefault::Required, NumForm::Lit(1), NumForm::Lit(2))]
    // u: Derived, m: Required
    #[case::dr_empty(NumForm::Undetermined, NumForm::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, NumForm::Undetermined, NumForm::Undetermined)]
    #[case::dr_u0(NumForm::Lit(0), NumForm::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, NumForm::Lit(0), NumForm::Undetermined)]
    #[case::dr_u1(NumForm::Lit(1), NumForm::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, NumForm::Lit(1), NumForm::Undetermined)]
    #[case::dr_m1(NumForm::Undetermined, NumForm::Lit(1), UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, NumForm::Lit(0), NumForm::Lit(1))]
    #[case::dr_m2(NumForm::Undetermined, NumForm::Lit(2), UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, NumForm::Lit(1), NumForm::Lit(2))]
    #[case::dr_u1m1(NumForm::Lit(1), NumForm::Lit(1), UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, NumForm::Lit(1), NumForm::Lit(1))]
    #[case::dr_u1m2(NumForm::Lit(1), NumForm::Lit(2), UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, NumForm::Lit(1), NumForm::Lit(2))]
    // u: Derived, m: Derived
    #[case::dd_empty(NumForm::Undetermined, NumForm::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived, NumForm::Undetermined, NumForm::Undetermined)]
    #[case::dd_u0(NumForm::Lit(0), NumForm::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived, NumForm::Lit(0), NumForm::Lit(1))]
    #[case::dd_u1(NumForm::Lit(1), NumForm::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived, NumForm::Lit(1), NumForm::Lit(2))]
    #[case::dd_m1(NumForm::Undetermined, NumForm::Lit(1), UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived, NumForm::Lit(0), NumForm::Lit(1))]
    #[case::dd_m2(NumForm::Undetermined, NumForm::Lit(2), UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived, NumForm::Lit(1), NumForm::Lit(2))]
    #[case::dd_u1m1(NumForm::Lit(1), NumForm::Lit(1), UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived, NumForm::Lit(1), NumForm::Lit(1))]
    #[case::dd_u1m2(NumForm::Lit(1), NumForm::Lit(2), UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived, NumForm::Lit(1), NumForm::Lit(2))]
    fn test_raise_unpaired_electrons(
        #[case] initial_count: NumForm,
        #[case] initial_multiplicity: NumForm,
        #[case] count_default: UnpairedElectronsDefault,
        #[case] multiplicity_default: MultiplicityDefault,
        #[case] expected_count: NumForm,
        #[case] expected_multiplicity: NumForm,
    ) {
        let mut unpaired_electrons = UnpairedElectronsForm {
            count: initial_count,
            multiplicity: initial_multiplicity,
        };
        raise_unpaired_electrons(
            &mut unpaired_electrons,
            count_default,
            multiplicity_default,
        );
        assert_eq!(unpaired_electrons.count, expected_count);
        assert_eq!(unpaired_electrons.multiplicity, expected_multiplicity);
    }

    /// Per-mode lowering: covers the AST states reachable by raising the canonical DSL fragments, plus (U, U) where applicable.
    #[rustfmt::skip]
    #[rstest]
    // u: Zero, m: Derived
    #[case::zd_derived_zero(NumForm::Lit(0), NumForm::Lit(1), UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived, NumForm::Undetermined, NumForm::Undetermined)]
    #[case::zd_derived_nonzero(NumForm::Lit(1), NumForm::Lit(2), UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived, NumForm::Lit(1), NumForm::Undetermined)]
    #[case::zd_zero_nonderived(NumForm::Lit(0), NumForm::Lit(2), UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived, NumForm::Undetermined, NumForm::Lit(2))]
    #[case::zd_nonzero_nonderived(NumForm::Lit(1), NumForm::Lit(1), UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived, NumForm::Lit(1), NumForm::Lit(1))]
    // u: Zero, m: Required
    #[case::zr_zero_mundef(NumForm::Lit(0), NumForm::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, NumForm::Undetermined, NumForm::Undetermined)]
    #[case::zr_nonzero_mundef(NumForm::Lit(1), NumForm::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, NumForm::Lit(1), NumForm::Undetermined)]
    #[case::zr_zero_derived(NumForm::Lit(0), NumForm::Lit(1), UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, NumForm::Undetermined, NumForm::Lit(1))]
    #[case::zr_zero_doublet(NumForm::Lit(0), NumForm::Lit(2), UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, NumForm::Undetermined, NumForm::Lit(2))]
    #[case::zr_nonzero_singlet(NumForm::Lit(1), NumForm::Lit(1), UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, NumForm::Lit(1), NumForm::Lit(1))]
    #[case::zr_nonzero_doublet(NumForm::Lit(1), NumForm::Lit(2), UnpairedElectronsDefault::Zero, MultiplicityDefault::Required, NumForm::Lit(1), NumForm::Lit(2))]
    // u: Required, m: Derived
    #[case::rd_both_undef(NumForm::Undetermined, NumForm::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, NumForm::Undetermined, NumForm::Undetermined)]
    #[case::rd_derived_zero(NumForm::Lit(0), NumForm::Lit(1), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, NumForm::Lit(0), NumForm::Undetermined)]
    #[case::rd_derived_nonzero(NumForm::Lit(1), NumForm::Lit(2), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, NumForm::Lit(1), NumForm::Undetermined)]
    #[case::rd_uundef_msinglet(NumForm::Undetermined, NumForm::Lit(1), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, NumForm::Undetermined, NumForm::Lit(1))]
    #[case::rd_uundef_mdoublet(NumForm::Undetermined, NumForm::Lit(2), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, NumForm::Undetermined, NumForm::Lit(2))]
    #[case::rd_nonderived(NumForm::Lit(1), NumForm::Lit(1), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived, NumForm::Lit(1), NumForm::Lit(1))]
    // u: Required, m: Required — nothing strips.
    #[case::rr_both_undef(NumForm::Undetermined, NumForm::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Required, NumForm::Undetermined, NumForm::Undetermined)]
    #[case::rr_full(NumForm::Lit(0), NumForm::Lit(1), UnpairedElectronsDefault::Required, MultiplicityDefault::Required, NumForm::Lit(0), NumForm::Lit(1))]
    // u: Derived, m: Required
    #[case::dr_both_undef(NumForm::Undetermined, NumForm::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, NumForm::Undetermined, NumForm::Undetermined)]
    #[case::dr_zero_mundef(NumForm::Lit(0), NumForm::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, NumForm::Lit(0), NumForm::Undetermined)]
    #[case::dr_nonzero_mundef(NumForm::Lit(1), NumForm::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, NumForm::Lit(1), NumForm::Undetermined)]
    #[case::dr_derived_zero(NumForm::Lit(0), NumForm::Lit(1), UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, NumForm::Undetermined, NumForm::Lit(1))]
    #[case::dr_derived_nonzero(NumForm::Lit(1), NumForm::Lit(2), UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, NumForm::Undetermined, NumForm::Lit(2))]
    #[case::dr_nonderived(NumForm::Lit(1), NumForm::Lit(1), UnpairedElectronsDefault::Derived, MultiplicityDefault::Required, NumForm::Lit(1), NumForm::Lit(1))]
    // u: Derived, m: Derived — tie-break keeps u explicit.
    #[case::dd_both_undef(NumForm::Undetermined, NumForm::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived, NumForm::Undetermined, NumForm::Undetermined)]
    #[case::dd_derived_zero(NumForm::Lit(0), NumForm::Lit(1), UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived, NumForm::Lit(0), NumForm::Undetermined)]
    #[case::dd_derived_nonzero(NumForm::Lit(1), NumForm::Lit(2), UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived, NumForm::Lit(1), NumForm::Undetermined)]
    #[case::dd_nonderived(NumForm::Lit(1), NumForm::Lit(1), UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived, NumForm::Lit(1), NumForm::Lit(1))]
    fn test_lower_unpaired_electrons(
        #[case] initial_count: NumForm,
        #[case] initial_multiplicity: NumForm,
        #[case] count_default: UnpairedElectronsDefault,
        #[case] multiplicity_default: MultiplicityDefault,
        #[case] expected_count: NumForm,
        #[case] expected_multiplicity: NumForm,
    ) {
        let mut unpaired_electrons = UnpairedElectronsForm {
            count: initial_count,
            multiplicity: initial_multiplicity,
        };
        lower_unpaired_electrons(
            &mut unpaired_electrons,
            count_default,
            multiplicity_default,
        );
        assert_eq!(unpaired_electrons.count, expected_count);
        assert_eq!(unpaired_electrons.multiplicity, expected_multiplicity);
    }

    /// AST preservation: the raised AST is a fixed point of `lower → raise`. Lowering strips default content; re-raising the result must recover the same AST.
    #[rustfmt::skip]
    #[rstest]
    #[case::zd_empty(NumForm::Undetermined, NumForm::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived)]
    #[case::zd_u0(NumForm::Lit(0), NumForm::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived)]
    #[case::zd_u1(NumForm::Lit(1), NumForm::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived)]
    #[case::zd_m1(NumForm::Undetermined, NumForm::Lit(1), UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived)]
    #[case::zd_m2(NumForm::Undetermined, NumForm::Lit(2), UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived)]
    #[case::zd_u1m1(NumForm::Lit(1), NumForm::Lit(1), UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived)]
    #[case::zd_u1m2(NumForm::Lit(1), NumForm::Lit(2), UnpairedElectronsDefault::Zero, MultiplicityDefault::Derived)]
    #[case::zr_empty(NumForm::Undetermined, NumForm::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Required)]
    #[case::zr_u0(NumForm::Lit(0), NumForm::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Required)]
    #[case::zr_u1(NumForm::Lit(1), NumForm::Undetermined, UnpairedElectronsDefault::Zero, MultiplicityDefault::Required)]
    #[case::zr_m1(NumForm::Undetermined, NumForm::Lit(1), UnpairedElectronsDefault::Zero, MultiplicityDefault::Required)]
    #[case::zr_m2(NumForm::Undetermined, NumForm::Lit(2), UnpairedElectronsDefault::Zero, MultiplicityDefault::Required)]
    #[case::zr_u1m1(NumForm::Lit(1), NumForm::Lit(1), UnpairedElectronsDefault::Zero, MultiplicityDefault::Required)]
    #[case::zr_u1m2(NumForm::Lit(1), NumForm::Lit(2), UnpairedElectronsDefault::Zero, MultiplicityDefault::Required)]
    #[case::rd_empty(NumForm::Undetermined, NumForm::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Derived)]
    #[case::rd_u0(NumForm::Lit(0), NumForm::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Derived)]
    #[case::rd_u1(NumForm::Lit(1), NumForm::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Derived)]
    #[case::rd_m1(NumForm::Undetermined, NumForm::Lit(1), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived)]
    #[case::rd_m2(NumForm::Undetermined, NumForm::Lit(2), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived)]
    #[case::rd_u1m1(NumForm::Lit(1), NumForm::Lit(1), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived)]
    #[case::rd_u1m2(NumForm::Lit(1), NumForm::Lit(2), UnpairedElectronsDefault::Required, MultiplicityDefault::Derived)]
    #[case::rr_empty(NumForm::Undetermined, NumForm::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Required)]
    #[case::rr_u0(NumForm::Lit(0), NumForm::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Required)]
    #[case::rr_u1(NumForm::Lit(1), NumForm::Undetermined, UnpairedElectronsDefault::Required, MultiplicityDefault::Required)]
    #[case::rr_m1(NumForm::Undetermined, NumForm::Lit(1), UnpairedElectronsDefault::Required, MultiplicityDefault::Required)]
    #[case::rr_m2(NumForm::Undetermined, NumForm::Lit(2), UnpairedElectronsDefault::Required, MultiplicityDefault::Required)]
    #[case::rr_u1m1(NumForm::Lit(1), NumForm::Lit(1), UnpairedElectronsDefault::Required, MultiplicityDefault::Required)]
    #[case::rr_u1m2(NumForm::Lit(1), NumForm::Lit(2), UnpairedElectronsDefault::Required, MultiplicityDefault::Required)]
    #[case::dr_empty(NumForm::Undetermined, NumForm::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Required)]
    #[case::dr_u0(NumForm::Lit(0), NumForm::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Required)]
    #[case::dr_u1(NumForm::Lit(1), NumForm::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Required)]
    #[case::dr_m1(NumForm::Undetermined, NumForm::Lit(1), UnpairedElectronsDefault::Derived, MultiplicityDefault::Required)]
    #[case::dr_m2(NumForm::Undetermined, NumForm::Lit(2), UnpairedElectronsDefault::Derived, MultiplicityDefault::Required)]
    #[case::dr_u1m1(NumForm::Lit(1), NumForm::Lit(1), UnpairedElectronsDefault::Derived, MultiplicityDefault::Required)]
    #[case::dr_u1m2(NumForm::Lit(1), NumForm::Lit(2), UnpairedElectronsDefault::Derived, MultiplicityDefault::Required)]
    #[case::dd_empty(NumForm::Undetermined, NumForm::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived)]
    #[case::dd_u0(NumForm::Lit(0), NumForm::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived)]
    #[case::dd_u1(NumForm::Lit(1), NumForm::Undetermined, UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived)]
    #[case::dd_m1(NumForm::Undetermined, NumForm::Lit(1), UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived)]
    #[case::dd_m2(NumForm::Undetermined, NumForm::Lit(2), UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived)]
    #[case::dd_u1m1(NumForm::Lit(1), NumForm::Lit(1), UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived)]
    #[case::dd_u1m2(NumForm::Lit(1), NumForm::Lit(2), UnpairedElectronsDefault::Derived, MultiplicityDefault::Derived)]
    fn test_unpaired_electrons_defaults_roundtrip(
        #[case] initial_count: NumForm,
        #[case] initial_multiplicity: NumForm,
        #[case] count_default: UnpairedElectronsDefault,
        #[case] multiplicity_default: MultiplicityDefault,
    ) {
        let mut raised = UnpairedElectronsForm {
            count: initial_count,
            multiplicity: initial_multiplicity,
        };
        raise_unpaired_electrons(&mut raised, count_default, multiplicity_default);

        let mut lowered_then_raised = raised.clone();
        lower_unpaired_electrons(
            &mut lowered_then_raised,
            count_default,
            multiplicity_default,
        );
        raise_unpaired_electrons(
            &mut lowered_then_raised,
            count_default,
            multiplicity_default,
        );

        assert_eq!(lowered_then_raised, raised);
    }
}
