//! Unpaired-electron form and component updates.

use umol_chem::spin::{SpinState, UnpairedElectrons};
use umol_graph_ir_macros::{Canonicalize, Lattice};

use super::num::NumForm;
use super::traits::{AsLit, Canonicalize};

/// Unpaired-electron count and multiplicity as independent `NumForm` fields.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Lattice, Canonicalize)]
pub struct UnpairedElectronsForm {
    pub count: NumForm,
    pub multiplicity: NumForm,
}

impl UnpairedElectronsForm {
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
            (NumForm::Undetermined, NumForm::Lit(m)) => {
                let unpaired_electrons = *m - 1;
                self.count = NumForm::Lit(unpaired_electrons);
            }
            (NumForm::Lit(u), NumForm::Undetermined) => {
                let multiplicity = u + 1;
                self.multiplicity = NumForm::Lit(multiplicity);
            }
            _ => {}
        }
    }
}

impl Default for UnpairedElectronsForm {
    fn default() -> Self {
        Self {
            count: NumForm::Undetermined,
            multiplicity: NumForm::Undetermined,
        }
    }
}

impl From<(u8, u8)> for UnpairedElectronsForm {
    fn from((count, multiplicity): (u8, u8)) -> Self {
        Self {
            count: NumForm::Lit(i64::from(count)),
            multiplicity: NumForm::Lit(multiplicity as i64),
        }
    }
}

impl From<UnpairedElectrons> for UnpairedElectronsForm {
    fn from(unpaired_electrons: UnpairedElectrons) -> Self {
        Self {
            count: NumForm::Lit(unpaired_electrons.count),
            multiplicity: NumForm::Lit(unpaired_electrons.multiplicity),
        }
    }
}

impl From<SpinState> for UnpairedElectronsForm {
    fn from(state: SpinState) -> Self {
        UnpairedElectrons::from(state).into()
    }
}

impl AsLit for UnpairedElectronsForm {
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
    pub count: Option<NumForm>,
    pub multiplicity: Option<NumForm>,
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::spin;

    use super::*;
    use crate::ir::error::Contradiction;
    use crate::ir::num::ArithExpr;
    use crate::ir::traits::{Canonicalize, Lattice};

    #[rstest]
    fn test_unpaired_electrons_form_closed_shell() {
        assert_eq!(
            UnpairedElectronsForm::closed_shell(),
            UnpairedElectronsForm {
                count: NumForm::Lit(0),
                multiplicity: NumForm::Lit(1),
            }
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::count((2_u8, 3_u8).into(), UnpairedElectronsUpdate { count: Some(NumForm::Lit(0)), multiplicity: None }, (0_u8, 3_u8).into())]
    #[case::multiplicity((2_u8, 3_u8).into(), UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Lit(1)) }, (2_u8, 1_u8).into())]
    #[case::count_undetermined((2_u8, 3_u8).into(), UnpairedElectronsUpdate { count: Some(NumForm::Undetermined), multiplicity: None }, UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(3) })]
    #[case::multiplicity_undetermined((2_u8, 3_u8).into(), UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Undetermined) }, UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Undetermined })]
    #[case::both((2_u8, 3_u8).into(), UnpairedElectronsUpdate { count: Some(NumForm::Lit(0)), multiplicity: Some(NumForm::Lit(1)) }, (0_u8, 1_u8).into())]
    fn test_unpaired_electrons_form_update(
        #[case] unpaired_electrons: UnpairedElectronsForm,
        #[case] update: UnpairedElectronsUpdate,
        #[case] expected: UnpairedElectronsForm,
    ) {
        assert_eq!(unpaired_electrons.update(&update), expected);
    }

    #[rstest]
    #[case::empty((2_u8, 3_u8).into())]
    fn test_unpaired_electrons_form_update_identity(
        #[case] unpaired_electrons: UnpairedElectronsForm,
    ) {
        assert_eq!(
            unpaired_electrons.update(&UnpairedElectronsUpdate::default()),
            unpaired_electrons,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::multiplicity((2_u8, 3_u8).into(), (2_u8, 1_u8).into(), UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Lit(1)) })]
    #[case::count_undetermined((2_u8, 3_u8).into(), UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(3) }, UnpairedElectronsUpdate { count: Some(NumForm::Undetermined), multiplicity: None })]
    #[case::both((2_u8, 3_u8).into(), (0_u8, 1_u8).into(), UnpairedElectronsUpdate { count: Some(NumForm::Lit(0)), multiplicity: Some(NumForm::Lit(1)) })]
    fn test_unpaired_electrons_form_difference_to(
        #[case] unpaired_electrons: UnpairedElectronsForm,
        #[case] other: UnpairedElectronsForm,
        #[case] expected: UnpairedElectronsUpdate,
    ) {
        assert_eq!(unpaired_electrons.difference_to(&other), expected);
    }

    #[rstest]
    #[case::canonical(
        UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Lit(1) },
        UnpairedElectronsForm { count: NumForm::lit_set([2]), multiplicity: NumForm::lit_set([1]) },
    )]
    fn test_unpaired_electrons_form_difference_to_identity(
        #[case] unpaired_electrons: UnpairedElectronsForm,
        #[case] other: UnpairedElectronsForm,
    ) {
        assert_eq!(
            unpaired_electrons.difference_to(&other),
            UnpairedElectronsUpdate::default()
        );
    }

    #[rstest]
    #[case::count_undetermined(UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(3) }, UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Lit(3) })]
    #[case::multiplicity_undetermined(UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Undetermined }, UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Lit(3) })]
    #[case::both_determined(UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Lit(3) }, UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Lit(3) })]
    fn test_unpaired_electrons_form_high_spin_complete(
        #[case] input: UnpairedElectronsForm,
        #[case] expected: UnpairedElectronsForm,
    ) {
        let mut ast = input.clone();
        ast.high_spin_complete();
        assert_eq!(ast, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::folds_field(
        UnpairedElectronsForm { count: NumForm::arith_expr(ArithExpr::Sum(vec![ArithExpr::Lit(1), ArithExpr::Lit(1)])), multiplicity: NumForm::Lit(3) },
        Ok(UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Lit(3) }),
    )]
    #[case::parity_invalid_is_allowed(
        UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(1) },
        Ok(UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(1) }),
    )]
    fn test_unpaired_electrons_form_canonicalize(
        #[case] input: UnpairedElectronsForm,
        #[case] expected: Result<UnpairedElectronsForm, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::both_undetermined(UnpairedElectronsForm::default())]
    #[case::valid_triplet(UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Lit(3) })]
    #[case::parity_invalid(UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(1) })]
    fn test_unpaired_electrons_form_canonicalize_identity(#[case] input: UnpairedElectronsForm) {
        assert_eq!(input.clone().canonicalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_undetermined(UnpairedElectronsForm::default(), false)]
    #[case::from_dsl(spin!("#u2").into(), true)]
    #[case::from_pair((2, 3).into(), true)]
    #[case::partial(UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Undetermined }, false)]
    fn test_unpaired_electrons_form_is_ground(#[case] ast: UnpairedElectronsForm, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::both_undetermined(UnpairedElectronsForm::default(), true)]
    #[case::ground((2_u8, 3_u8).into(), false)]
    #[case::partial(UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Undetermined }, false)]
    fn test_unpaired_electrons_form_is_undetermined(#[case] ast: UnpairedElectronsForm, #[case] expected: bool) {
        assert_eq!(ast.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_ground(UnpairedElectronsForm::default(), (2_u8, 3_u8).into(), Some((2_u8, 3_u8).into()))]
    #[case::exact((2_u8, 3_u8).into(), (2_u8, 3_u8).into(), Some((2_u8, 3_u8).into()))]
    #[case::count_conflict((2_u8, 3_u8).into(), (0_u8, 1_u8).into(), None)]
    #[case::field_wise(
        UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(3) },
        UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Undetermined },
        Some((2_u8, 3_u8).into()),
    )]
    fn test_unpaired_electrons_form_meet(
        #[case] a: UnpairedElectronsForm,
        #[case] b: UnpairedElectronsForm,
        #[case] expected: Option<UnpairedElectronsForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und((2_u8, 3_u8).into(), UnpairedElectronsForm::default(), UnpairedElectronsForm::default())]
    #[case::exact((2_u8, 3_u8).into(), (2_u8, 3_u8).into(), (2_u8, 3_u8).into())]
    #[case::field_wise_widen(
        (2_u8, 3_u8).into(),
        UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Lit(1) },
        UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::lit_set([1, 3]) },
    )]
    fn test_unpaired_electrons_form_join(
        #[case] a: UnpairedElectronsForm,
        #[case] b: UnpairedElectronsForm,
        #[case] expected: UnpairedElectronsForm,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::both_undetermined(UnpairedElectronsForm::default(), (2_u8, 3_u8).into(), true)]
    #[case::pattern_specific_target_undetermined((2_u8, 3_u8).into(), UnpairedElectronsForm::default(), false)]
    #[case::exact((2_u8, 3_u8).into(), (2_u8, 3_u8).into(), true)]
    #[case::count_mismatch((2_u8, 3_u8).into(), (0_u8, 3_u8).into(), false)]
    #[case::partial_pattern(
        UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(3) },
        (2_u8, 3_u8).into(),
        true,
    )]
    fn test_unpaired_electrons_form_matches(
        #[case] pattern: UnpairedElectronsForm,
        #[case] target: UnpairedElectronsForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    #[case::closed_shell(
        UnpairedElectrons { count: 0, multiplicity: 1 },
        UnpairedElectronsForm { count: NumForm::Lit(0), multiplicity: NumForm::Lit(1) }
    )]
    #[case::structurally_invalid(
        UnpairedElectrons { count: -1, multiplicity: 0 },
        UnpairedElectronsForm { count: NumForm::Lit(-1), multiplicity: NumForm::Lit(0) }
    )]
    fn test_unpaired_electrons_form_from_unpaired_electrons(
        #[case] unpaired_electrons: UnpairedElectrons,
        #[case] expected: UnpairedElectronsForm,
    ) {
        assert_eq!(UnpairedElectronsForm::from(unpaired_electrons), expected);
    }

    #[rstest]
    #[case::closed_shell(
        SpinState::closed_shell(),
        UnpairedElectronsForm { count: NumForm::Lit(0), multiplicity: NumForm::Lit(1) }
    )]
    #[case::triplet(
        spin!("#u2#s3"),
        UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Lit(3) }
    )]
    fn test_unpaired_electrons_form_from_spin_state(
        #[case] spin_state: SpinState,
        #[case] expected: UnpairedElectronsForm,
    ) {
        assert_eq!(UnpairedElectronsForm::from(spin_state), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::both_undetermined(UnpairedElectronsForm::default(), None)]
    #[case::count_lit_only(
        UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Undetermined },
        None,
    )]
    #[case::multiplicity_lit_only(
        UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(3) },
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
        UnpairedElectronsForm { count: NumForm::Lit(-1), multiplicity: NumForm::Lit(1) },
        Some(UnpairedElectrons { count: -1, multiplicity: 1 }),
    )]
    #[case::count_out_of_range(
        UnpairedElectronsForm { count: NumForm::Lit(256), multiplicity: NumForm::Lit(1) },
        Some(UnpairedElectrons { count: 256, multiplicity: 1 }),
    )]
    #[case::zero_multiplicity(
        UnpairedElectronsForm { count: NumForm::Lit(0), multiplicity: NumForm::Lit(0) },
        Some(UnpairedElectrons { count: 0, multiplicity: 0 }),
    )]
    #[case::multiplicity_out_of_range(
        UnpairedElectronsForm { count: NumForm::Lit(0), multiplicity: NumForm::Lit(256) },
        Some(UnpairedElectrons { count: 0, multiplicity: 256 }),
    )]
    fn test_unpaired_electrons_form_as_lit(
        #[case] ast: UnpairedElectronsForm,
        #[case] expected: Option<UnpairedElectrons>,
    ) {
        assert_eq!(ast.as_lit(), expected);
    }
}
