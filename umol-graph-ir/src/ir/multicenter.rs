//! Multicenter bond AST.

use umol_graph_core::{ParticipantPosition, RelationData};
use umol_graph_ir_macros::{Canonicalize, Lattice};

use super::constraint::{MulticenterBondConstraintAst, MulticenterBondConstraintsAst};
use super::electrons::ElectronCountsAst;
use super::spin::{UnpairedElectronsAst, UnpairedElectronsUpdate};
use super::traits::{Canonicalize, Lattice};
use super::value::NumForm;

#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Canonicalize, Lattice)]
pub struct MulticenterBondAst {
    pub electrons: ElectronCountsAst,
    pub charge: NumForm,
    pub unpaired_electrons: UnpairedElectronsAst,
    pub constraints: MulticenterBondConstraintsAst,
}

/// Attribute update for a multicenter bond. Ordinary fields are optional,
/// unpaired-electron components are updated independently, and undetermined constraints remove
/// their key.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MulticenterBondUpdate {
    pub electrons: Option<ElectronCountsAst>,
    pub charge: Option<NumForm>,
    pub unpaired_electrons: UnpairedElectronsUpdate,
    pub constraints: MulticenterBondConstraintsAst,
}

impl From<&str> for MulticenterBondAst {
    fn from(s: &str) -> Self {
        s.parse().expect("invalid multicenter bond string")
    }
}

impl RelationData for MulticenterBondAst {
    /// The per-member electron counts are positional, so they follow a participant reorder.
    fn on_permutation(&mut self, order: &[ParticipantPosition]) {
        self.electrons.permute(order);
    }

    fn is_permutation_invariant(&self) -> bool {
        self.electrons.is_undetermined()
    }
}

impl MulticenterBondAst {
    pub fn new(electrons: ElectronCountsAst) -> Self {
        Self {
            electrons,
            ..Default::default()
        }
    }

    pub fn from_electrons(electrons: Vec<i64>) -> Self {
        Self::new(ElectronCountsAst::Lit(electrons))
    }

    pub fn with_charge(mut self, charge: impl Into<NumForm>) -> Self {
        self.charge = charge.into();
        self
    }

    pub fn with_unpaired_electrons(
        mut self,
        unpaired_electrons: impl Into<UnpairedElectronsAst>,
    ) -> Self {
        self.unpaired_electrons = unpaired_electrons.into();
        self
    }

    /// Add a single constraint, replacing any existing entry of the same
    /// kind (last-wins per `MulticenterBondConstraintsAst::set`). Chainable.
    pub fn with_constraint(mut self, constraint: impl Into<MulticenterBondConstraintAst>) -> Self {
        self.constraints.set(constraint.into());
        self
    }

    /// Add each constraint from the iterator, replacing any existing entry
    /// of the same kind (last-wins per `MulticenterBondConstraintsAst::set`).
    /// Does not clear existing constraints; use `bond.constraints.clear()`
    /// or direct field assignment for wipe-and-replace.
    pub fn with_constraints<I>(mut self, constraints: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<MulticenterBondConstraintAst>,
    {
        for c in constraints {
            self.constraints.set(c.into());
        }
        self
    }

    /// Fill `Undetermined` value-bearing struct fields with zero defaults:
    /// charge → `Lit(0)`, unpaired electrons → closed-shell singlet `(0, 1)`. `electrons`
    /// and `constraints` are preserved. The result is ground iff `electrons`
    /// is already `Lit`.
    pub fn into_ground(mut self) -> Self {
        if self.charge.is_undetermined() {
            self.charge = NumForm::Lit(0);
        }
        if self.unpaired_electrons.is_undetermined() {
            self.unpaired_electrons = UnpairedElectronsAst::from((0_u8, 1_u8));
        }
        self
    }

    /// Apply an attribute update, leaving omitted leaves and constraint keys unchanged.
    pub fn update(&self, update: &MulticenterBondUpdate) -> MulticenterBondAst {
        let mut constraints = self.constraints.clone();
        constraints.update(&update.constraints);
        MulticenterBondAst {
            electrons: update
                .electrons
                .clone()
                .unwrap_or_else(|| self.electrons.clone()),
            charge: update.charge.clone().unwrap_or_else(|| self.charge.clone()),
            unpaired_electrons: self.unpaired_electrons.update(&update.unpaired_electrons),
            constraints,
        }
    }

    /// Derive the minimal canonical attribute update carrying `self` to `other`.
    pub fn difference_to(&self, other: &Self) -> MulticenterBondUpdate {
        let mut constraints = MulticenterBondConstraintsAst::new();
        for new in other.constraints.iter() {
            if self
                .constraints
                .get(new.key())
                .is_none_or(|old| !old.canonical_eq(new))
            {
                constraints.set(new.clone());
            }
        }
        for old in self.constraints.iter() {
            if other.constraints.get(old.key()).is_none() {
                constraints.set(old.as_undetermined());
            }
        }
        MulticenterBondUpdate {
            electrons: (!self.electrons.canonical_eq(&other.electrons))
                .then(|| other.electrons.clone()),
            charge: (!self.charge.canonical_eq(&other.charge)).then(|| other.charge.clone()),
            unpaired_electrons: self
                .unpaired_electrons
                .difference_to(&other.unpaired_electrons),
            constraints,
        }
    }

    /// Reorder the positional `electrons` by `order`, tracking a participant
    /// reordering; charge / unpaired electrons / constraints are positionless and unchanged.
    pub fn permute(&mut self, order: &[ParticipantPosition]) {
        self.electrons.permute(order);
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ir::error::Contradiction;
    use crate::ir::traits::Canonicalize;

    #[rustfmt::skip]
    #[rstest]
    #[case::literal(MulticenterBondAst::new(ElectronCountsAst::Lit(vec![1; 3])),
        MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1; 3]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsAst::default(),
            constraints: MulticenterBondConstraintsAst::new() })]
    fn test_multicenter_bond_ast_new(
        #[case] actual: MulticenterBondAst,
        #[case] expected: MulticenterBondAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::literal(MulticenterBondAst::from_electrons(vec![1, 1, 1]),
        MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1; 3]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsAst::default(),
            constraints: MulticenterBondConstraintsAst::new() })]
    fn test_multicenter_bond_ast_from_electrons(
        #[case] actual: MulticenterBondAst,
        #[case] expected: MulticenterBondAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::with_charge(MulticenterBondAst::from_electrons(vec![1, 1, 1]).with_charge(-1),
        MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]),
            charge: NumForm::Lit(-1), unpaired_electrons: UnpairedElectronsAst::default(),
            constraints: MulticenterBondConstraintsAst::new() })]
    #[case::with_unpaired_electrons(MulticenterBondAst::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((0_u8, 1_u8)),
        MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsAst::closed_shell(),
            constraints: MulticenterBondConstraintsAst::new() })]
    #[case::with_constraint(
        MulticenterBondAst::from_electrons(vec![1, 1, 1]).with_constraint(MulticenterBondConstraintAst::electron_count(2)),
        MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsAst::default(),
            constraints: MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(2)) })]
    #[case::with_constraints_extends(
        MulticenterBondAst::from_electrons(vec![1, 1, 1])
            .with_constraints([MulticenterBondConstraintAst::electron_count(2)]),
        MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsAst::default(),
            constraints: MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(2)) })]
    #[case::with_constraint_replaces_same_kind(
        MulticenterBondAst::from_electrons(vec![1, 1, 1])
            .with_constraint(MulticenterBondConstraintAst::electron_count(2))
            .with_constraint(MulticenterBondConstraintAst::electron_count(4)),
        MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsAst::default(),
            constraints: MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(4)) })]
    fn test_multicenter_bond_ast_with_methods(
        #[case] actual: MulticenterBondAst,
        #[case] expected: MulticenterBondAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::from_ground_electrons(
        MulticenterBondAst::from_electrons(vec![1; 3]).into_ground(),
        MulticenterBondAst {
            electrons: ElectronCountsAst::Lit(vec![1; 3]),
            charge: NumForm::Lit(0),
            unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)),
            constraints: MulticenterBondConstraintsAst::new(),
        },
    )]
    #[case::preserves_set_charge(
        MulticenterBondAst::from_electrons(vec![1; 3]).with_charge(1_i64).into_ground(),
        MulticenterBondAst {
            electrons: ElectronCountsAst::Lit(vec![1; 3]),
            charge: NumForm::Lit(1),
            unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)),
            constraints: MulticenterBondConstraintsAst::new(),
        },
    )]
    #[case::preserves_constraints(
        MulticenterBondAst::from_electrons(vec![1; 3])
            .with_constraint(MulticenterBondConstraintAst::electron_count(3))
            .into_ground(),
        MulticenterBondAst {
            electrons: ElectronCountsAst::Lit(vec![1; 3]),
            charge: NumForm::Lit(0),
            unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)),
            constraints: MulticenterBondConstraintsAst::from(
                MulticenterBondConstraintAst::electron_count(3),
            ),
        },
    )]
    fn test_multicenter_bond_ast_into_ground(
        #[case] actual: MulticenterBondAst,
        #[case] expected: MulticenterBondAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electrons(MulticenterBondAst::from_electrons(vec![1, 1, 1]), MulticenterBondUpdate { electrons: Some(ElectronCountsAst::Lit(vec![2, 2, 2])), ..Default::default() }, MulticenterBondAst::from_electrons(vec![2, 2, 2]))]
    #[case::electrons_undetermined(MulticenterBondAst::from_electrons(vec![1, 1, 1]), MulticenterBondUpdate { electrons: Some(ElectronCountsAst::Undetermined), ..Default::default() }, MulticenterBondAst::default())]
    #[case::charge(MulticenterBondAst::from_electrons(vec![1, 1, 1]).with_charge(0_i64), MulticenterBondUpdate { charge: Some(NumForm::Lit(-1)), ..Default::default() }, MulticenterBondAst::from_electrons(vec![1, 1, 1]).with_charge(-1_i64))]
    #[case::charge_undetermined(MulticenterBondAst::from_electrons(vec![1, 1, 1]).with_charge(-1_i64), MulticenterBondUpdate { charge: Some(NumForm::Undetermined), ..Default::default() }, MulticenterBondAst::from_electrons(vec![1, 1, 1]))]
    #[case::unpaired_electrons_count(MulticenterBondAst::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((2_u8, 3_u8)), MulticenterBondUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: Some(NumForm::Lit(0)), multiplicity: None }, ..Default::default() }, MulticenterBondAst::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((0_u8, 3_u8)))]
    #[case::unpaired_electrons_multiplicity(MulticenterBondAst::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((2_u8, 3_u8)), MulticenterBondUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Lit(1)) }, ..Default::default() }, MulticenterBondAst::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((2_u8, 1_u8)))]
    #[case::constraint_set(MulticenterBondAst::from_electrons(vec![1, 1, 1]), MulticenterBondUpdate { constraints: MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6_i64)), ..Default::default() }, MulticenterBondAst::from_electrons(vec![1, 1, 1]).with_constraint(MulticenterBondConstraintAst::electron_count(6_i64)))]
    #[case::constraint_replace(MulticenterBondAst::from_electrons(vec![1, 1, 1]).with_constraint(MulticenterBondConstraintAst::electron_count(6_i64)), MulticenterBondUpdate { constraints: MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(4_i64)), ..Default::default() }, MulticenterBondAst::from_electrons(vec![1, 1, 1]).with_constraint(MulticenterBondConstraintAst::electron_count(4_i64)))]
    #[case::constraint_remove(MulticenterBondAst::from_electrons(vec![1, 1, 1]).with_constraint(MulticenterBondConstraintAst::electron_count(6_i64)), MulticenterBondUpdate { constraints: MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(NumForm::Undetermined)), ..Default::default() }, MulticenterBondAst::from_electrons(vec![1, 1, 1]))]
    fn test_multicenter_bond_ast_update(
        #[case] bond: MulticenterBondAst,
        #[case] update: MulticenterBondUpdate,
        #[case] expected: MulticenterBondAst,
    ) {
        assert_eq!(bond.update(&update), expected);
    }

    #[rstest]
    #[case::empty(MulticenterBondAst::from_electrons(vec![1, 1, 1]).with_charge(-1_i64).with_unpaired_electrons((2_u8, 3_u8)).with_constraint(MulticenterBondConstraintAst::electron_count(6_i64)))]
    fn test_multicenter_bond_ast_update_identity(#[case] bond: MulticenterBondAst) {
        assert_eq!(bond.update(&MulticenterBondUpdate::default()), bond);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fields_and_constraints(
        MulticenterBondAst::from_electrons(vec![1, 1, 1]).with_charge(0_i64).with_unpaired_electrons((2_u8, 3_u8)).with_constraint(MulticenterBondConstraintAst::electron_count(6_i64)),
        MulticenterBondAst::from_electrons(vec![2, 2, 2]).with_unpaired_electrons((2_u8, 1_u8)),
        MulticenterBondUpdate {
            electrons: Some(ElectronCountsAst::Lit(vec![2, 2, 2])),
            charge: Some(NumForm::Undetermined),
            unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Lit(1)) },
            constraints: MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(NumForm::Undetermined)),
        },
    )]
    fn test_multicenter_bond_ast_difference_to(
        #[case] bond: MulticenterBondAst,
        #[case] other: MulticenterBondAst,
        #[case] expected: MulticenterBondUpdate,
    ) {
        assert_eq!(bond.difference_to(&other), expected);
    }

    #[rstest]
    #[case::canonical(
        MulticenterBondAst::from_electrons(vec![1, 1, 1]).with_charge(1_i64),
        MulticenterBondAst::from_electrons(vec![1, 1, 1]).with_charge(NumForm::lit_set([1])),
    )]
    fn test_multicenter_bond_ast_difference_to_identity(
        #[case] bond: MulticenterBondAst,
        #[case] other: MulticenterBondAst,
    ) {
        assert_eq!(bond.difference_to(&other), MulticenterBondUpdate::default());
    }

    #[rstest]
    #[case::three_members(
        MulticenterBondAst::from_electrons(vec![10, 20, 30]).with_charge(-1),
        vec![
            ParticipantPosition(2),
            ParticipantPosition(0),
            ParticipantPosition(1),
        ],
        MulticenterBondAst::from_electrons(vec![30, 10, 20]).with_charge(-1),
    )]
    fn test_multicenter_bond_ast_permute(
        #[case] mut input: MulticenterBondAst,
        #[case] order: Vec<ParticipantPosition>,
        #[case] expected: MulticenterBondAst,
    ) {
        input.permute(&order);
        assert_eq!(input, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_(MulticenterBondAst::default(), false)]
    #[case::charge_only(MulticenterBondAst::new(ElectronCountsAst::Undetermined).with_charge(0), false)]
    #[case::ground_no_atoms(MulticenterBondAst::new(ElectronCountsAst::Lit(Vec::new())).with_charge(0).with_unpaired_electrons((0, 1)), true)]
    #[case::all_ground_three(
        MulticenterBondAst::new(ElectronCountsAst::Lit(vec![1; 3])).with_charge(0).with_unpaired_electrons((0, 1)),
        true,
    )]
    #[case::ground_with_constraint(
        MulticenterBondAst::new(ElectronCountsAst::Lit(vec![1; 3]))
            .with_charge(0).with_unpaired_electrons((0, 1))
            .with_constraint(MulticenterBondConstraintAst::electron_count(3)),
        true,
    )]
    fn test_multicenter_bond_ast_is_ground(
        #[case] ast: MulticenterBondAst,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::folds_charge(
        MulticenterBondAst::default().with_charge(NumForm::lit_set([0])),
        Ok(MulticenterBondAst::default().with_charge(0)),
    )]
    #[case::charge_empty_litset_contradiction(
        MulticenterBondAst::default().with_charge(NumForm::lit_set(Vec::<i64>::new())),
        Err(Contradiction),
    )]
    fn test_multicenter_bond_ast_canonicalize(
        #[case] input: MulticenterBondAst,
        #[case] expected: Result<MulticenterBondAst, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_matches_default(MulticenterBondAst::default(), MulticenterBondAst::default(), true)]
    #[case::default_matches_ground(
        MulticenterBondAst::default(),
        MulticenterBondAst::new(ElectronCountsAst::Lit(Vec::new())).with_charge(0).with_unpaired_electrons((0, 1)),
        true,
    )]
    #[case::exact(
        MulticenterBondAst::new(ElectronCountsAst::Lit(vec![1; 3])).with_charge(0).with_unpaired_electrons((0, 1)),
        MulticenterBondAst::new(ElectronCountsAst::Lit(vec![1; 3])).with_charge(0).with_unpaired_electrons((0, 1)),
        true,
    )]
    #[case::electrons_length_mismatch(
        MulticenterBondAst::new(ElectronCountsAst::Lit(vec![1; 2])),
        MulticenterBondAst::new(ElectronCountsAst::Lit(vec![1; 3])).with_charge(0).with_unpaired_electrons((0, 1)),
        false,
    )]
    #[case::electrons_value_mismatch(
        MulticenterBondAst::new(ElectronCountsAst::Lit(vec![2; 3])),
        MulticenterBondAst::new(ElectronCountsAst::Lit(vec![1; 3])).with_charge(0).with_unpaired_electrons((0, 1)),
        false,
    )]
    #[case::charge_mismatch(
        MulticenterBondAst::new(ElectronCountsAst::Undetermined).with_charge(1),
        MulticenterBondAst::new(ElectronCountsAst::Undetermined).with_charge(0),
        false,
    )]
    #[case::unpaired_electrons_mismatch(
        MulticenterBondAst::new(ElectronCountsAst::Undetermined).with_unpaired_electrons((2_u8, 3_u8)),
        MulticenterBondAst::new(ElectronCountsAst::Undetermined).with_unpaired_electrons((0_u8, 1_u8)),
        false,
    )]
    #[case::constraint_required_present(
        MulticenterBondAst::new(ElectronCountsAst::Undetermined)
            .with_constraint(MulticenterBondConstraintAst::electron_count(3)),
        MulticenterBondAst::new(ElectronCountsAst::Undetermined)
            .with_constraint(MulticenterBondConstraintAst::electron_count(3)),
        true,
    )]
    #[case::constraint_required_absent(
        MulticenterBondAst::new(ElectronCountsAst::Undetermined)
            .with_constraint(MulticenterBondConstraintAst::electron_count(3)),
        MulticenterBondAst::new(ElectronCountsAst::Undetermined),
        false,
    )]
    fn test_multicenter_bond_ast_matches(
        #[case] pattern: MulticenterBondAst,
        #[case] target: MulticenterBondAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    #[case::both_default(
        MulticenterBondAst::default(),
        MulticenterBondAst::default(),
        Some(MulticenterBondAst::default())
    )]
    #[case::electrons_length_mismatch(
        MulticenterBondAst::new(ElectronCountsAst::Lit(vec![2; 3])),
        MulticenterBondAst::new(ElectronCountsAst::Lit(vec![2; 4])),
        None,
    )]
    #[case::narrows_electrons(
        MulticenterBondAst::new(ElectronCountsAst::Undetermined),
        MulticenterBondAst::from_electrons(vec![1, 2]),
        Some(MulticenterBondAst::from_electrons(vec![1, 2])),
    )]
    fn test_multicenter_bond_ast_meet(
        #[case] a: MulticenterBondAst,
        #[case] b: MulticenterBondAst,
        #[case] expected: Option<MulticenterBondAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::electrons_length_mismatch_widens_to_default(
        MulticenterBondAst::new(ElectronCountsAst::Lit(vec![2; 3])),
        MulticenterBondAst::new(ElectronCountsAst::Lit(vec![2; 4])),
        MulticenterBondAst::default(),
    )]
    fn test_multicenter_bond_ast_join(
        #[case] a: MulticenterBondAst,
        #[case] b: MulticenterBondAst,
        #[case] expected: MulticenterBondAst,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }
}
