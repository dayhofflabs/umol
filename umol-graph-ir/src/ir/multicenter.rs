//! Multicenter bonds: the molecule's collection and one bond's attribute form.

use std::sync::Arc;

use umol_graph_core::{
    EdgeId, GraphCompaction, NodeId, ParticipantPosition, RelationData, RelationId,
    RelationPushout, Remapping, Unordered, VarRelationSet,
};
use umol_graph_ir_macros::{Lattice, Normalize};

use super::constraint::{MulticenterBondConstraintForm, MulticenterBondConstraintsForm};
use super::electrons::ElectronCountsForm;
use super::id::AtomId;
use super::num::NumForm;
use super::spin::{UnpairedElectronsForm, UnpairedElectronsUpdate};
use super::traits::{Equiv, Lattice};

/// The molecule's multicenter bonds. The member atoms bear the frame, in parallel with the per-member electron
/// counts. There is no site.
///
/// Owns the frame structure its storage shape cannot state: which factor bears the participant
/// frame, and which, if any, is a site.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct MulticenterBonds(Arc<VarRelationSet<NodeId, Unordered, MulticenterBondForm>>);

impl From<Arc<VarRelationSet<NodeId, Unordered, MulticenterBondForm>>> for MulticenterBonds {
    fn from(set: Arc<VarRelationSet<NodeId, Unordered, MulticenterBondForm>>) -> Self {
        Self(set)
    }
}

impl MulticenterBonds {
    pub fn new(entries: Vec<(Vec<NodeId>, MulticenterBondForm)>) -> Self {
        Self(Arc::new(VarRelationSet::new(entries)))
    }

    pub fn participants(&self, id: RelationId) -> &[NodeId] {
        self.0.participants(id)
    }

    pub fn find_by_participants(&self, query: &[NodeId]) -> Option<RelationId> {
        self.0.find_by_participants(query)
    }

    pub fn participant_permutation(
        &self,
        id: RelationId,
        query: &[NodeId],
    ) -> Option<Vec<ParticipantPosition>> {
        self.0.participant_permutation(id, query)
    }

    pub fn pushout(
        &self,
        right: &Self,
        combine: impl FnMut(&MulticenterBondForm, &MulticenterBondForm) -> Option<MulticenterBondForm>,
    ) -> Option<RelationPushout<Self>> {
        self.0
            .pushout(&right.0, combine)
            .map(|pushout| RelationPushout {
                object: Self(Arc::new(pushout.object)),
                left: pushout.left,
                right: pushout.right,
            })
    }

    /// Glue `right`, relabelled into this molecule's id space, onto `self`: coinciding entries meet,
    /// non-coinciding entries are carried. `None` when a coincident meet is bottom.
    pub fn glue(
        &self,
        right: &Self,
        remapping: &Remapping,
    ) -> Option<Vec<(Vec<AtomId>, MulticenterBondForm)>> {
        let merged = self.pushout(&right.remap(remapping), |a, b| a.meet(b))?;
        Some(
            merged
                .object
                .relation_ids()
                .map(|id| {
                    (
                        merged
                            .object
                            .participants(id)
                            .iter()
                            .map(|&n| AtomId::from(n))
                            .collect(),
                        merged.object.data(id).clone(),
                    )
                })
                .collect(),
        )
    }

    pub fn count(&self) -> usize {
        self.0.count()
    }

    pub fn contains(&self, id: RelationId) -> bool {
        self.0.contains(id)
    }

    pub fn relation_ids(&self) -> impl ExactSizeIterator<Item = RelationId> {
        self.0.relation_ids()
    }

    pub fn data(&self, id: RelationId) -> &MulticenterBondForm {
        self.0.data(id)
    }

    pub fn data_mut(&mut self, id: RelationId) -> &mut MulticenterBondForm {
        Arc::make_mut(&mut self.0).data_mut(id)
    }

    pub fn data_iter_mut(&mut self) -> impl ExactSizeIterator<Item = &mut MulticenterBondForm> {
        Arc::make_mut(&mut self.0).data_iter_mut()
    }

    pub fn incident(&self, node: NodeId) -> &[RelationId] {
        self.0.incident(node)
    }

    pub fn incident_edge(&self, edge: EdgeId) -> &[RelationId] {
        self.0.incident_edge(edge)
    }

    pub fn has_incident(&self, node: NodeId) -> bool {
        self.0.has_incident(node)
    }

    pub fn has_incident_edge(&self, edge: EdgeId) -> bool {
        self.0.has_incident_edge(edge)
    }

    pub fn into_entries(self) -> Vec<(Vec<NodeId>, MulticenterBondForm)> {
        Arc::try_unwrap(self.0)
            .unwrap_or_else(|shared| (*shared).clone())
            .into_entries()
    }

    pub fn remap(&self, remapping: &Remapping) -> Self {
        Self(Arc::new(self.0.remap(remapping)))
    }

    pub fn compact(&self, compaction: &GraphCompaction) -> Self {
        Self(Arc::new(self.0.compact(compaction)))
    }

    pub fn into_arc(self) -> Arc<VarRelationSet<NodeId, Unordered, MulticenterBondForm>> {
        self.0
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Normalize, Lattice)]
pub struct MulticenterBondForm {
    pub electrons: ElectronCountsForm,
    pub charge: NumForm,
    pub unpaired_electrons: UnpairedElectronsForm,
    pub constraints: MulticenterBondConstraintsForm,
}

/// Attribute update for a multicenter bond. Ordinary fields are optional,
/// unpaired-electron components are updated independently, and undetermined constraints remove
/// their key.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MulticenterBondUpdate {
    pub electrons: Option<ElectronCountsForm>,
    pub charge: Option<NumForm>,
    pub unpaired_electrons: UnpairedElectronsUpdate,
    pub constraints: MulticenterBondConstraintsForm,
}

impl From<&str> for MulticenterBondForm {
    fn from(s: &str) -> Self {
        s.parse().expect("invalid multicenter bond string")
    }
}

impl RelationData for MulticenterBondForm {
    /// The per-member electron counts are positional, so they follow a participant reorder.
    fn on_permutation(&mut self, order: &[ParticipantPosition]) {
        self.electrons.permute(order);
    }

    fn is_permutation_invariant(&self) -> bool {
        self.electrons.is_undetermined()
    }
}

impl MulticenterBondForm {
    /// Concrete: every inherent field is ground; the constraint channel does
    /// not bear on concreteness.
    pub fn is_concrete(&self) -> bool {
        let MulticenterBondForm {
            electrons,
            charge,
            unpaired_electrons,
            constraints: _,
        } = self;
        electrons.is_ground() && charge.is_ground() && unpaired_electrons.is_ground()
    }
    pub fn new(electrons: ElectronCountsForm) -> Self {
        Self {
            electrons,
            ..Default::default()
        }
    }

    pub fn from_electrons(electrons: Vec<i64>) -> Self {
        Self::new(ElectronCountsForm::Lit(electrons))
    }

    pub fn with_charge(mut self, charge: impl Into<NumForm>) -> Self {
        self.charge = charge.into();
        self
    }

    pub fn with_unpaired_electrons(
        mut self,
        unpaired_electrons: impl Into<UnpairedElectronsForm>,
    ) -> Self {
        self.unpaired_electrons = unpaired_electrons.into();
        self
    }

    /// Add a single constraint, replacing any existing entry of the same
    /// kind (last-wins per `MulticenterBondConstraintsForm::set`). Chainable.
    pub fn with_constraint(mut self, constraint: impl Into<MulticenterBondConstraintForm>) -> Self {
        self.constraints.set(constraint.into());
        self
    }

    /// Add each constraint from the iterator, replacing any existing entry
    /// of the same kind (last-wins per `MulticenterBondConstraintsForm::set`).
    /// Does not clear existing constraints; use `bond.constraints.clear()`
    /// or direct field assignment for wipe-and-replace.
    pub fn with_constraints<I>(mut self, constraints: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<MulticenterBondConstraintForm>,
    {
        for c in constraints {
            self.constraints.set(c.into());
        }
        self
    }

    /// Fill `Undetermined` value-bearing struct fields with zero defaults:
    /// charge → `Lit(0)`, unpaired electrons → closed-shell singlet `(0, 1)`. `electrons`
    /// and `constraints` are preserved. The result is concrete iff `electrons`
    /// is already `Lit`.
    pub fn into_concrete(mut self) -> Self {
        if self.charge.is_undetermined() {
            self.charge = NumForm::Lit(0);
        }
        if self.unpaired_electrons.is_undetermined() {
            self.unpaired_electrons = UnpairedElectronsForm::from((0_u8, 1_u8));
        }
        self
    }

    /// Apply an attribute update, leaving omitted leaves and constraint keys unchanged.
    pub fn update(&self, update: &MulticenterBondUpdate) -> MulticenterBondForm {
        let mut constraints = self.constraints.clone();
        constraints.update(&update.constraints);
        MulticenterBondForm {
            electrons: update
                .electrons
                .clone()
                .unwrap_or_else(|| self.electrons.clone()),
            charge: update.charge.clone().unwrap_or_else(|| self.charge.clone()),
            unpaired_electrons: self.unpaired_electrons.update(&update.unpaired_electrons),
            constraints,
        }
    }

    /// Derive the minimal normalized attribute update carrying `self` to `other`.
    pub fn difference_to(&self, other: &Self) -> MulticenterBondUpdate {
        let mut constraints = MulticenterBondConstraintsForm::new();
        for new in other.constraints.iter() {
            if self
                .constraints
                .get(new.key())
                .is_none_or(|old| !old.equiv(new))
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
            electrons: (!self.electrons.equiv(&other.electrons)).then(|| other.electrons.clone()),
            charge: (!self.charge.equiv(&other.charge)).then(|| other.charge.clone()),
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
    use crate::ir::traits::Normalize;

    #[rustfmt::skip]
    #[rstest]
    #[case::literal(MulticenterBondForm::new(ElectronCountsForm::Lit(vec![1; 3])),
        MulticenterBondForm { electrons: ElectronCountsForm::Lit(vec![1; 3]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: MulticenterBondConstraintsForm::new() })]
    fn test_multicenter_bond_form_new(
        #[case] actual: MulticenterBondForm,
        #[case] expected: MulticenterBondForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::literal(MulticenterBondForm::from_electrons(vec![1, 1, 1]),
        MulticenterBondForm { electrons: ElectronCountsForm::Lit(vec![1; 3]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: MulticenterBondConstraintsForm::new() })]
    fn test_multicenter_bond_form_from_electrons(
        #[case] actual: MulticenterBondForm,
        #[case] expected: MulticenterBondForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::with_charge(MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_charge(-1),
        MulticenterBondForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]),
            charge: NumForm::Lit(-1), unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: MulticenterBondConstraintsForm::new() })]
    #[case::with_unpaired_electrons(MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((0_u8, 1_u8)),
        MulticenterBondForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::closed_shell(),
            constraints: MulticenterBondConstraintsForm::new() })]
    #[case::with_constraint(
        MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_constraint(MulticenterBondConstraintForm::electron_count(2)),
        MulticenterBondForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(2)) })]
    #[case::with_constraints_extends(
        MulticenterBondForm::from_electrons(vec![1, 1, 1])
            .with_constraints([MulticenterBondConstraintForm::electron_count(2)]),
        MulticenterBondForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(2)) })]
    #[case::with_constraint_replaces_same_kind(
        MulticenterBondForm::from_electrons(vec![1, 1, 1])
            .with_constraint(MulticenterBondConstraintForm::electron_count(2))
            .with_constraint(MulticenterBondConstraintForm::electron_count(4)),
        MulticenterBondForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(4)) })]
    fn test_multicenter_bond_form_with_methods(
        #[case] actual: MulticenterBondForm,
        #[case] expected: MulticenterBondForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::from_ground_electrons(
        MulticenterBondForm::from_electrons(vec![1; 3]).into_concrete(),
        MulticenterBondForm {
            electrons: ElectronCountsForm::Lit(vec![1; 3]),
            charge: NumForm::Lit(0),
            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
            constraints: MulticenterBondConstraintsForm::new(),
        },
    )]
    #[case::preserves_set_charge(
        MulticenterBondForm::from_electrons(vec![1; 3]).with_charge(1_i64).into_concrete(),
        MulticenterBondForm {
            electrons: ElectronCountsForm::Lit(vec![1; 3]),
            charge: NumForm::Lit(1),
            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
            constraints: MulticenterBondConstraintsForm::new(),
        },
    )]
    #[case::preserves_constraints(
        MulticenterBondForm::from_electrons(vec![1; 3])
            .with_constraint(MulticenterBondConstraintForm::electron_count(3))
            .into_concrete(),
        MulticenterBondForm {
            electrons: ElectronCountsForm::Lit(vec![1; 3]),
            charge: NumForm::Lit(0),
            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
            constraints: MulticenterBondConstraintsForm::from(
                MulticenterBondConstraintForm::electron_count(3),
            ),
        },
    )]
    fn test_multicenter_bond_form_into_concrete(
        #[case] actual: MulticenterBondForm,
        #[case] expected: MulticenterBondForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electrons(MulticenterBondForm::from_electrons(vec![1, 1, 1]), MulticenterBondUpdate { electrons: Some(ElectronCountsForm::Lit(vec![2, 2, 2])), ..Default::default() }, MulticenterBondForm::from_electrons(vec![2, 2, 2]))]
    #[case::electrons_undetermined(MulticenterBondForm::from_electrons(vec![1, 1, 1]), MulticenterBondUpdate { electrons: Some(ElectronCountsForm::Undetermined), ..Default::default() }, MulticenterBondForm::default())]
    #[case::charge(MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_charge(0_i64), MulticenterBondUpdate { charge: Some(NumForm::Lit(-1)), ..Default::default() }, MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_charge(-1_i64))]
    #[case::charge_undetermined(MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_charge(-1_i64), MulticenterBondUpdate { charge: Some(NumForm::Undetermined), ..Default::default() }, MulticenterBondForm::from_electrons(vec![1, 1, 1]))]
    #[case::unpaired_electrons_count(MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((2_u8, 3_u8)), MulticenterBondUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: Some(NumForm::Lit(0)), multiplicity: None }, ..Default::default() }, MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((0_u8, 3_u8)))]
    #[case::unpaired_electrons_multiplicity(MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((2_u8, 3_u8)), MulticenterBondUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Lit(1)) }, ..Default::default() }, MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((2_u8, 1_u8)))]
    #[case::constraint_set(MulticenterBondForm::from_electrons(vec![1, 1, 1]), MulticenterBondUpdate { constraints: MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6_i64)), ..Default::default() }, MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_constraint(MulticenterBondConstraintForm::electron_count(6_i64)))]
    #[case::constraint_replace(MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_constraint(MulticenterBondConstraintForm::electron_count(6_i64)), MulticenterBondUpdate { constraints: MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(4_i64)), ..Default::default() }, MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_constraint(MulticenterBondConstraintForm::electron_count(4_i64)))]
    #[case::constraint_remove(MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_constraint(MulticenterBondConstraintForm::electron_count(6_i64)), MulticenterBondUpdate { constraints: MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(NumForm::Undetermined)), ..Default::default() }, MulticenterBondForm::from_electrons(vec![1, 1, 1]))]
    fn test_multicenter_bond_form_update(
        #[case] bond: MulticenterBondForm,
        #[case] update: MulticenterBondUpdate,
        #[case] expected: MulticenterBondForm,
    ) {
        assert_eq!(bond.update(&update), expected);
    }

    #[rstest]
    #[case::empty(MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_charge(-1_i64).with_unpaired_electrons((2_u8, 3_u8)).with_constraint(MulticenterBondConstraintForm::electron_count(6_i64)))]
    fn test_multicenter_bond_form_update_identity(#[case] bond: MulticenterBondForm) {
        assert_eq!(bond.update(&MulticenterBondUpdate::default()), bond);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fields_and_constraints(
        MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_charge(0_i64).with_unpaired_electrons((2_u8, 3_u8)).with_constraint(MulticenterBondConstraintForm::electron_count(6_i64)),
        MulticenterBondForm::from_electrons(vec![2, 2, 2]).with_unpaired_electrons((2_u8, 1_u8)),
        MulticenterBondUpdate {
            electrons: Some(ElectronCountsForm::Lit(vec![2, 2, 2])),
            charge: Some(NumForm::Undetermined),
            unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Lit(1)) },
            constraints: MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(NumForm::Undetermined)),
        },
    )]
    fn test_multicenter_bond_form_difference_to(
        #[case] bond: MulticenterBondForm,
        #[case] other: MulticenterBondForm,
        #[case] expected: MulticenterBondUpdate,
    ) {
        assert_eq!(bond.difference_to(&other), expected);
    }

    #[rstest]
    #[case::normalized(
        MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_charge(1_i64),
        MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_charge(NumForm::lit_set([1])),
    )]
    fn test_multicenter_bond_form_difference_to_identity(
        #[case] bond: MulticenterBondForm,
        #[case] other: MulticenterBondForm,
    ) {
        assert_eq!(bond.difference_to(&other), MulticenterBondUpdate::default());
    }

    #[rstest]
    #[case::three_members(
        MulticenterBondForm::from_electrons(vec![10, 20, 30]).with_charge(-1),
        vec![
            ParticipantPosition(2),
            ParticipantPosition(0),
            ParticipantPosition(1),
        ],
        MulticenterBondForm::from_electrons(vec![30, 10, 20]).with_charge(-1),
    )]
    fn test_multicenter_bond_form_permute(
        #[case] mut input: MulticenterBondForm,
        #[case] order: Vec<ParticipantPosition>,
        #[case] expected: MulticenterBondForm,
    ) {
        input.permute(&order);
        assert_eq!(input, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_(MulticenterBondForm::default(), false)]
    #[case::charge_only(MulticenterBondForm::new(ElectronCountsForm::Undetermined).with_charge(0), false)]
    #[case::ground_no_atoms(MulticenterBondForm::new(ElectronCountsForm::Lit(Vec::new())).with_charge(0).with_unpaired_electrons((0, 1)), true)]
    #[case::all_ground_three(
        MulticenterBondForm::new(ElectronCountsForm::Lit(vec![1; 3])).with_charge(0).with_unpaired_electrons((0, 1)),
        true,
    )]
    #[case::ground_with_constraint(
        MulticenterBondForm::new(ElectronCountsForm::Lit(vec![1; 3]))
            .with_charge(0).with_unpaired_electrons((0, 1))
            .with_constraint(MulticenterBondConstraintForm::electron_count(3)),
        true,
    )]
    fn test_multicenter_bond_form_is_ground(
        #[case] form: MulticenterBondForm,
        #[case] expected: bool,
    ) {
        assert_eq!(form.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::folds_charge(
        MulticenterBondForm::default().with_charge(NumForm::lit_set([0])),
        Ok(MulticenterBondForm::default().with_charge(0)),
    )]
    #[case::charge_empty_litset_contradiction(
        MulticenterBondForm::default().with_charge(NumForm::lit_set(Vec::<i64>::new())),
        Err(Contradiction),
    )]
    fn test_multicenter_bond_form_normalize(
        #[case] input: MulticenterBondForm,
        #[case] expected: Result<MulticenterBondForm, Contradiction>,
    ) {
        assert_eq!(input.normalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_matches_default(MulticenterBondForm::default(), MulticenterBondForm::default(), true)]
    #[case::default_matches_ground(
        MulticenterBondForm::default(),
        MulticenterBondForm::new(ElectronCountsForm::Lit(Vec::new())).with_charge(0).with_unpaired_electrons((0, 1)),
        true,
    )]
    #[case::exact(
        MulticenterBondForm::new(ElectronCountsForm::Lit(vec![1; 3])).with_charge(0).with_unpaired_electrons((0, 1)),
        MulticenterBondForm::new(ElectronCountsForm::Lit(vec![1; 3])).with_charge(0).with_unpaired_electrons((0, 1)),
        true,
    )]
    #[case::electrons_length_mismatch(
        MulticenterBondForm::new(ElectronCountsForm::Lit(vec![1; 2])),
        MulticenterBondForm::new(ElectronCountsForm::Lit(vec![1; 3])).with_charge(0).with_unpaired_electrons((0, 1)),
        false,
    )]
    #[case::electrons_value_mismatch(
        MulticenterBondForm::new(ElectronCountsForm::Lit(vec![2; 3])),
        MulticenterBondForm::new(ElectronCountsForm::Lit(vec![1; 3])).with_charge(0).with_unpaired_electrons((0, 1)),
        false,
    )]
    #[case::charge_mismatch(
        MulticenterBondForm::new(ElectronCountsForm::Undetermined).with_charge(1),
        MulticenterBondForm::new(ElectronCountsForm::Undetermined).with_charge(0),
        false,
    )]
    #[case::unpaired_electrons_mismatch(
        MulticenterBondForm::new(ElectronCountsForm::Undetermined).with_unpaired_electrons((2_u8, 3_u8)),
        MulticenterBondForm::new(ElectronCountsForm::Undetermined).with_unpaired_electrons((0_u8, 1_u8)),
        false,
    )]
    #[case::constraint_required_present(
        MulticenterBondForm::new(ElectronCountsForm::Undetermined)
            .with_constraint(MulticenterBondConstraintForm::electron_count(3)),
        MulticenterBondForm::new(ElectronCountsForm::Undetermined)
            .with_constraint(MulticenterBondConstraintForm::electron_count(3)),
        true,
    )]
    #[case::constraint_required_absent(
        MulticenterBondForm::new(ElectronCountsForm::Undetermined)
            .with_constraint(MulticenterBondConstraintForm::electron_count(3)),
        MulticenterBondForm::new(ElectronCountsForm::Undetermined),
        false,
    )]
    fn test_multicenter_bond_form_matches(
        #[case] pattern: MulticenterBondForm,
        #[case] target: MulticenterBondForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    #[case::both_default(
        MulticenterBondForm::default(),
        MulticenterBondForm::default(),
        Some(MulticenterBondForm::default())
    )]
    #[case::electrons_length_mismatch(
        MulticenterBondForm::new(ElectronCountsForm::Lit(vec![2; 3])),
        MulticenterBondForm::new(ElectronCountsForm::Lit(vec![2; 4])),
        None,
    )]
    #[case::narrows_electrons(
        MulticenterBondForm::new(ElectronCountsForm::Undetermined),
        MulticenterBondForm::from_electrons(vec![1, 2]),
        Some(MulticenterBondForm::from_electrons(vec![1, 2])),
    )]
    fn test_multicenter_bond_form_meet(
        #[case] a: MulticenterBondForm,
        #[case] b: MulticenterBondForm,
        #[case] expected: Option<MulticenterBondForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::electrons_length_mismatch_widens_to_default(
        MulticenterBondForm::new(ElectronCountsForm::Lit(vec![2; 3])),
        MulticenterBondForm::new(ElectronCountsForm::Lit(vec![2; 4])),
        MulticenterBondForm::default(),
    )]
    fn test_multicenter_bond_form_join(
        #[case] a: MulticenterBondForm,
        #[case] b: MulticenterBondForm,
        #[case] expected: MulticenterBondForm,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }
}
