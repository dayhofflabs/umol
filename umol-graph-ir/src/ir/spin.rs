//! Unpaired-electron AST and component updates.

use umol_chem::spin::{SpinState, UnpairedElectrons};
use umol_graph_ir_macros::{Canonicalize, Lattice};

use super::traits::{AsLit, Canonicalize};
use super::value::ValueAst;

/// Unpaired-electron count and multiplicity as independent `ValueAst` fields.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Lattice, Canonicalize)]
pub struct UnpairedElectronsAst {
    pub count: ValueAst,
    pub multiplicity: ValueAst,
}

impl UnpairedElectronsAst {
    pub fn closed_shell() -> Self {
        UnpairedElectrons {
            count: 0,
            multiplicity: 1,
        }
        .into()
    }

    /// Apply an update independently to the unpaired-electron and multiplicity components.
    pub fn update(&self, update: &UnpairedElectronsUpdate) -> Self {
        Self {
            count: update.count.clone().unwrap_or_else(|| self.count.clone()),
            multiplicity: update
                .multiplicity
                .clone()
                .unwrap_or_else(|| self.multiplicity.clone()),
        }
    }

    /// Derive the minimal canonical component update carrying `self` to `other`.
    pub fn difference_to(&self, other: &Self) -> UnpairedElectronsUpdate {
        UnpairedElectronsUpdate {
            count: (!self.count.canonical_eq(&other.count)).then(|| other.count.clone()),
            multiplicity: (!self.multiplicity.canonical_eq(&other.multiplicity))
                .then(|| other.multiplicity.clone()),
        }
    }

    pub fn high_spin_complete(&mut self) {
        match (&self.count, &self.multiplicity) {
            (ValueAst::Undetermined, ValueAst::Lit(m)) => {
                let unpaired_electrons = *m - 1;
                self.count = ValueAst::Lit(unpaired_electrons);
            }
            (ValueAst::Lit(u), ValueAst::Undetermined) => {
                let multiplicity = u + 1;
                self.multiplicity = ValueAst::Lit(multiplicity);
            }
            _ => {}
        }
    }
}

impl Default for UnpairedElectronsAst {
    fn default() -> Self {
        Self {
            count: ValueAst::Undetermined,
            multiplicity: ValueAst::Undetermined,
        }
    }
}

impl From<(u8, u8)> for UnpairedElectronsAst {
    fn from((count, multiplicity): (u8, u8)) -> Self {
        Self {
            count: ValueAst::Lit(i64::from(count)),
            multiplicity: ValueAst::Lit(multiplicity as i64),
        }
    }
}

impl From<UnpairedElectrons> for UnpairedElectronsAst {
    fn from(unpaired_electrons: UnpairedElectrons) -> Self {
        Self {
            count: ValueAst::Lit(unpaired_electrons.count),
            multiplicity: ValueAst::Lit(unpaired_electrons.multiplicity),
        }
    }
}

impl From<SpinState> for UnpairedElectronsAst {
    fn from(state: SpinState) -> Self {
        UnpairedElectrons::from(state).into()
    }
}

impl AsLit for UnpairedElectronsAst {
    type Lit = UnpairedElectrons;

    /// Exact unpaired-electron components when both fields are literal.
    #[inline]
    fn as_lit(&self) -> Option<UnpairedElectrons> {
        Some(UnpairedElectrons {
            count: self.count.as_lit()?,
            multiplicity: self.multiplicity.as_lit()?,
        })
    }
}

/// Leaf-wise update for unpaired-electron components. `None` leaves that component unchanged;
/// `Some(value)` sets it exactly, including to `Undetermined`.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct UnpairedElectronsUpdate {
    pub count: Option<ValueAst>,
    pub multiplicity: Option<ValueAst>,
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::spin;

    use super::*;
    use crate::ir::error::Contradiction;
    use crate::ir::traits::{Canonicalize, Lattice};
    use crate::ir::value::ValueTerm;

    #[rstest]
    fn test_unpaired_electrons_ast_closed_shell() {
        assert_eq!(
            UnpairedElectronsAst::closed_shell(),
            UnpairedElectronsAst {
                count: ValueAst::Lit(0),
                multiplicity: ValueAst::Lit(1),
            }
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::count((2_u8, 3_u8).into(), UnpairedElectronsUpdate { count: Some(ValueAst::Lit(0)), multiplicity: None }, (0_u8, 3_u8).into())]
    #[case::multiplicity((2_u8, 3_u8).into(), UnpairedElectronsUpdate { count: None, multiplicity: Some(ValueAst::Lit(1)) }, (2_u8, 1_u8).into())]
    #[case::count_undetermined((2_u8, 3_u8).into(), UnpairedElectronsUpdate { count: Some(ValueAst::Undetermined), multiplicity: None }, UnpairedElectronsAst { count: ValueAst::Undetermined, multiplicity: ValueAst::Lit(3) })]
    #[case::multiplicity_undetermined((2_u8, 3_u8).into(), UnpairedElectronsUpdate { count: None, multiplicity: Some(ValueAst::Undetermined) }, UnpairedElectronsAst { count: ValueAst::Lit(2), multiplicity: ValueAst::Undetermined })]
    #[case::both((2_u8, 3_u8).into(), UnpairedElectronsUpdate { count: Some(ValueAst::Lit(0)), multiplicity: Some(ValueAst::Lit(1)) }, (0_u8, 1_u8).into())]
    fn test_unpaired_electrons_ast_update(
        #[case] unpaired_electrons: UnpairedElectronsAst,
        #[case] update: UnpairedElectronsUpdate,
        #[case] expected: UnpairedElectronsAst,
    ) {
        assert_eq!(unpaired_electrons.update(&update), expected);
    }

    #[rstest]
    #[case::empty((2_u8, 3_u8).into())]
    fn test_unpaired_electrons_ast_update_identity(
        #[case] unpaired_electrons: UnpairedElectronsAst,
    ) {
        assert_eq!(
            unpaired_electrons.update(&UnpairedElectronsUpdate::default()),
            unpaired_electrons,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::multiplicity((2_u8, 3_u8).into(), (2_u8, 1_u8).into(), UnpairedElectronsUpdate { count: None, multiplicity: Some(ValueAst::Lit(1)) })]
    #[case::count_undetermined((2_u8, 3_u8).into(), UnpairedElectronsAst { count: ValueAst::Undetermined, multiplicity: ValueAst::Lit(3) }, UnpairedElectronsUpdate { count: Some(ValueAst::Undetermined), multiplicity: None })]
    #[case::both((2_u8, 3_u8).into(), (0_u8, 1_u8).into(), UnpairedElectronsUpdate { count: Some(ValueAst::Lit(0)), multiplicity: Some(ValueAst::Lit(1)) })]
    fn test_unpaired_electrons_ast_difference_to(
        #[case] unpaired_electrons: UnpairedElectronsAst,
        #[case] other: UnpairedElectronsAst,
        #[case] expected: UnpairedElectronsUpdate,
    ) {
        assert_eq!(unpaired_electrons.difference_to(&other), expected);
    }

    #[rstest]
    #[case::canonical(
        UnpairedElectronsAst { count: ValueAst::Lit(2), multiplicity: ValueAst::Lit(1) },
        UnpairedElectronsAst { count: ValueAst::lit_set([2]), multiplicity: ValueAst::lit_set([1]) },
    )]
    fn test_unpaired_electrons_ast_difference_to_identity(
        #[case] unpaired_electrons: UnpairedElectronsAst,
        #[case] other: UnpairedElectronsAst,
    ) {
        assert_eq!(
            unpaired_electrons.difference_to(&other),
            UnpairedElectronsUpdate::default()
        );
    }

    #[rstest]
    #[case::count_undetermined(UnpairedElectronsAst { count: ValueAst::Undetermined, multiplicity: ValueAst::Lit(3) }, UnpairedElectronsAst { count: ValueAst::Lit(2), multiplicity: ValueAst::Lit(3) })]
    #[case::multiplicity_undetermined(UnpairedElectronsAst { count: ValueAst::Lit(2), multiplicity: ValueAst::Undetermined }, UnpairedElectronsAst { count: ValueAst::Lit(2), multiplicity: ValueAst::Lit(3) })]
    #[case::both_determined(UnpairedElectronsAst { count: ValueAst::Lit(2), multiplicity: ValueAst::Lit(3) }, UnpairedElectronsAst { count: ValueAst::Lit(2), multiplicity: ValueAst::Lit(3) })]
    fn test_unpaired_electrons_ast_high_spin_complete(
        #[case] input: UnpairedElectronsAst,
        #[case] expected: UnpairedElectronsAst,
    ) {
        let mut ast = input.clone();
        ast.high_spin_complete();
        assert_eq!(ast, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::folds_field(
        UnpairedElectronsAst { count: ValueAst::term(ValueTerm::Sum(vec![ValueTerm::Lit(1), ValueTerm::Lit(1)])), multiplicity: ValueAst::Lit(3) },
        Ok(UnpairedElectronsAst { count: ValueAst::Lit(2), multiplicity: ValueAst::Lit(3) }),
    )]
    #[case::parity_invalid_is_allowed(
        UnpairedElectronsAst { count: ValueAst::Lit(1), multiplicity: ValueAst::Lit(1) },
        Ok(UnpairedElectronsAst { count: ValueAst::Lit(1), multiplicity: ValueAst::Lit(1) }),
    )]
    fn test_unpaired_electrons_ast_canonicalize(
        #[case] input: UnpairedElectronsAst,
        #[case] expected: Result<UnpairedElectronsAst, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::both_undetermined(UnpairedElectronsAst::default())]
    #[case::valid_triplet(UnpairedElectronsAst { count: ValueAst::Lit(2), multiplicity: ValueAst::Lit(3) })]
    #[case::parity_invalid(UnpairedElectronsAst { count: ValueAst::Lit(1), multiplicity: ValueAst::Lit(1) })]
    fn test_unpaired_electrons_ast_canonicalize_identity(#[case] input: UnpairedElectronsAst) {
        assert_eq!(input.clone().canonicalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_undetermined(UnpairedElectronsAst::default(), false)]
    #[case::from_dsl(spin!("#u2").into(), true)]
    #[case::from_pair((2, 3).into(), true)]
    #[case::partial(UnpairedElectronsAst { count: ValueAst::Lit(2), multiplicity: ValueAst::Undetermined }, false)]
    fn test_unpaired_electrons_ast_is_ground(#[case] ast: UnpairedElectronsAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::both_undetermined(UnpairedElectronsAst::default(), true)]
    #[case::ground((2_u8, 3_u8).into(), false)]
    #[case::partial(UnpairedElectronsAst { count: ValueAst::Lit(2), multiplicity: ValueAst::Undetermined }, false)]
    fn test_unpaired_electrons_ast_is_undetermined(#[case] ast: UnpairedElectronsAst, #[case] expected: bool) {
        assert_eq!(ast.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_ground(UnpairedElectronsAst::default(), (2_u8, 3_u8).into(), Some((2_u8, 3_u8).into()))]
    #[case::exact((2_u8, 3_u8).into(), (2_u8, 3_u8).into(), Some((2_u8, 3_u8).into()))]
    #[case::count_conflict((2_u8, 3_u8).into(), (0_u8, 1_u8).into(), None)]
    #[case::field_wise(
        UnpairedElectronsAst { count: ValueAst::Undetermined, multiplicity: ValueAst::Lit(3) },
        UnpairedElectronsAst { count: ValueAst::Lit(2), multiplicity: ValueAst::Undetermined },
        Some((2_u8, 3_u8).into()),
    )]
    fn test_unpaired_electrons_ast_meet(
        #[case] a: UnpairedElectronsAst,
        #[case] b: UnpairedElectronsAst,
        #[case] expected: Option<UnpairedElectronsAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und((2_u8, 3_u8).into(), UnpairedElectronsAst::default(), UnpairedElectronsAst::default())]
    #[case::exact((2_u8, 3_u8).into(), (2_u8, 3_u8).into(), (2_u8, 3_u8).into())]
    #[case::field_wise_widen(
        (2_u8, 3_u8).into(),
        UnpairedElectronsAst { count: ValueAst::Lit(2), multiplicity: ValueAst::Lit(1) },
        UnpairedElectronsAst { count: ValueAst::Lit(2), multiplicity: ValueAst::lit_set([1, 3]) },
    )]
    fn test_unpaired_electrons_ast_join(
        #[case] a: UnpairedElectronsAst,
        #[case] b: UnpairedElectronsAst,
        #[case] expected: UnpairedElectronsAst,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::both_undetermined(UnpairedElectronsAst::default(), (2_u8, 3_u8).into(), true)]
    #[case::pattern_specific_target_undetermined((2_u8, 3_u8).into(), UnpairedElectronsAst::default(), false)]
    #[case::exact((2_u8, 3_u8).into(), (2_u8, 3_u8).into(), true)]
    #[case::count_mismatch((2_u8, 3_u8).into(), (0_u8, 3_u8).into(), false)]
    #[case::partial_pattern(
        UnpairedElectronsAst { count: ValueAst::Undetermined, multiplicity: ValueAst::Lit(3) },
        (2_u8, 3_u8).into(),
        true,
    )]
    fn test_unpaired_electrons_ast_matches(
        #[case] pattern: UnpairedElectronsAst,
        #[case] target: UnpairedElectronsAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    #[case::closed_shell(
        UnpairedElectrons { count: 0, multiplicity: 1 },
        UnpairedElectronsAst { count: ValueAst::Lit(0), multiplicity: ValueAst::Lit(1) }
    )]
    #[case::structurally_invalid(
        UnpairedElectrons { count: -1, multiplicity: 0 },
        UnpairedElectronsAst { count: ValueAst::Lit(-1), multiplicity: ValueAst::Lit(0) }
    )]
    fn test_unpaired_electrons_ast_from_unpaired_electrons(
        #[case] unpaired_electrons: UnpairedElectrons,
        #[case] expected: UnpairedElectronsAst,
    ) {
        assert_eq!(UnpairedElectronsAst::from(unpaired_electrons), expected);
    }

    #[rstest]
    #[case::closed_shell(
        SpinState::closed_shell(),
        UnpairedElectronsAst { count: ValueAst::Lit(0), multiplicity: ValueAst::Lit(1) }
    )]
    #[case::triplet(
        spin!("#u2#s3"),
        UnpairedElectronsAst { count: ValueAst::Lit(2), multiplicity: ValueAst::Lit(3) }
    )]
    fn test_unpaired_electrons_ast_from_spin_state(
        #[case] spin_state: SpinState,
        #[case] expected: UnpairedElectronsAst,
    ) {
        assert_eq!(UnpairedElectronsAst::from(spin_state), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::both_undetermined(UnpairedElectronsAst::default(), None)]
    #[case::count_lit_only(
        UnpairedElectronsAst { count: ValueAst::Lit(2), multiplicity: ValueAst::Undetermined },
        None,
    )]
    #[case::multiplicity_lit_only(
        UnpairedElectronsAst { count: ValueAst::Undetermined, multiplicity: ValueAst::Lit(3) },
        None,
    )]
    #[case::valid_triplet(
        (2_u8, 3_u8).into(),
        Some(UnpairedElectrons { count: 2, multiplicity: 3 })
    )]
    #[case::valid_closed_shell(
        (0_u8, 1_u8).into(),
        Some(UnpairedElectrons { count: 0, multiplicity: 1 })
    )]
    #[case::parity_invalid(
        (1_u8, 1_u8).into(),
        Some(UnpairedElectrons { count: 1, multiplicity: 1 })
    )]
    #[case::negative_count(
        UnpairedElectronsAst { count: ValueAst::Lit(-1), multiplicity: ValueAst::Lit(1) },
        Some(UnpairedElectrons { count: -1, multiplicity: 1 }),
    )]
    #[case::count_out_of_range(
        UnpairedElectronsAst { count: ValueAst::Lit(256), multiplicity: ValueAst::Lit(1) },
        Some(UnpairedElectrons { count: 256, multiplicity: 1 }),
    )]
    #[case::zero_multiplicity(
        UnpairedElectronsAst { count: ValueAst::Lit(0), multiplicity: ValueAst::Lit(0) },
        Some(UnpairedElectrons { count: 0, multiplicity: 0 }),
    )]
    #[case::multiplicity_out_of_range(
        UnpairedElectronsAst { count: ValueAst::Lit(0), multiplicity: ValueAst::Lit(256) },
        Some(UnpairedElectrons { count: 0, multiplicity: 256 }),
    )]
    fn test_unpaired_electrons_ast_as_lit(
        #[case] ast: UnpairedElectronsAst,
        #[case] expected: Option<UnpairedElectrons>,
    ) {
        assert_eq!(ast.as_lit(), expected);
    }
}
