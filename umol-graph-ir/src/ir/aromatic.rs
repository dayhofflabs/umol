//! Aromatic systems: the molecule's collection and one system's attribute form.

use std::sync::Arc;

use umol_graph_core::{
    EdgeId, GraphCompaction, NodeId, ParticipantPosition, RelationData, RelationId,
    RelationPushout, Remapping, Unordered, VarRelationSet,
};
use umol_graph_ir_macros::{Lattice, Normalize};

use super::constraint::{AromaticSystemConstraintForm, AromaticSystemConstraintsForm};
use super::electrons::ElectronCountsForm;
use super::id::{AromaticSystemId, AtomId};
use super::num::NumForm;
use super::spin::{UnpairedElectronsForm, UnpairedElectronsUpdate};
use super::traits::{Equiv, Lattice};

/// The molecule's aromatic systems.
///
/// The atoms bear the participant frame: the per-member electron counts of
/// [`AromaticSystemForm`] are read against it, position by position.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct AromaticSystems(Arc<VarRelationSet<NodeId, Unordered, AromaticSystemForm>>);

impl AromaticSystems {
    pub fn new(entries: Vec<(Vec<AtomId>, AromaticSystemForm)>) -> Self {
        Self(Arc::new(VarRelationSet::new(
            entries
                .into_iter()
                .map(|(atoms, attributes)| {
                    (atoms.into_iter().map(NodeId::from).collect(), attributes)
                })
                .collect(),
        )))
    }

    pub fn count(&self) -> usize {
        self.0.count()
    }

    pub fn contains(&self, id: AromaticSystemId) -> bool {
        self.0.contains(RelationId::from(id))
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = AromaticSystemId> {
        self.0.relation_ids().map(AromaticSystemId::from)
    }

    /// The atoms of `id`, in their stored frame.
    pub fn atoms(&self, id: AromaticSystemId) -> impl ExactSizeIterator<Item = AtomId> + '_ {
        self.0
            .participants(RelationId::from(id))
            .iter()
            .map(|&atom| AtomId::from(atom))
    }

    pub fn attributes(&self, id: AromaticSystemId) -> &AromaticSystemForm {
        self.0.data(RelationId::from(id))
    }

    pub fn attributes_mut(&mut self, id: AromaticSystemId) -> &mut AromaticSystemForm {
        Arc::make_mut(&mut self.0).data_mut(RelationId::from(id))
    }

    /// Ids of the systems `atom` belongs to. Systems are atom-disjoint, so there is at most one.
    pub fn incident_ids(
        &self,
        atom: AtomId,
    ) -> impl ExactSizeIterator<Item = AromaticSystemId> + '_ {
        self.0
            .incident(NodeId::from(atom))
            .iter()
            .map(|&id| AromaticSystemId::from(id))
    }

    pub fn has_incident(&self, atom: AtomId) -> bool {
        self.0.has_incident(NodeId::from(atom))
    }

    pub fn into_entries(self) -> Vec<(Vec<AtomId>, AromaticSystemForm)> {
        Arc::try_unwrap(self.0)
            .unwrap_or_else(|shared| (*shared).clone())
            .into_entries()
            .into_iter()
            .map(|(atoms, attributes)| {
                (atoms.into_iter().map(AtomId::from).collect(), attributes)
            })
            .collect()
    }

    /// The atoms of `id` as graph nodes, for graph-core interop that is not yet typed in graph-IR
    /// ids. The public accessor is [`Self::atoms`].
    pub(crate) fn atom_nodes(&self, id: AromaticSystemId) -> &[NodeId] {
        self.0.participants(RelationId::from(id))
    }

    pub(crate) fn attributes_iter_mut(
        &mut self,
    ) -> impl ExactSizeIterator<Item = &mut AromaticSystemForm> {
        Arc::make_mut(&mut self.0).data_iter_mut()
    }

    pub(crate) fn remap(&self, remapping: &Remapping) -> Self {
        Self(Arc::new(self.0.remap(remapping)))
    }

    pub(crate) fn compact(&self, compaction: &GraphCompaction) -> Self {
        Self(Arc::new(self.0.compact(compaction)))
    }

    pub(crate) fn into_arc(self) -> Arc<VarRelationSet<NodeId, Unordered, AromaticSystemForm>> {
        self.0
    }

    /// Glue `right`, relabelled into this molecule's id space, onto `self`: coinciding systems meet,
    /// non-coinciding systems are carried. `None` when a coincident meet is bottom.
    pub(crate) fn glue(
        &self,
        right: &Self,
        remapping: &Remapping,
    ) -> Option<Vec<(Vec<AtomId>, AromaticSystemForm)>> {
        let merged = self
            .0
            .pushout(&right.remap(remapping).0, |a, b| a.meet(b))?;
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
                            .map(|&atom| AtomId::from(atom))
                            .collect(),
                        merged.object.data(id).clone(),
                    )
                })
                .collect(),
        )
    }

    /// Removed by the lookup relocation: replaced by `of_id` keyed on the uniqueness key, which for
    /// an aromatic system is any member atom.
    pub(crate) fn find_by_participants(&self, atoms: &[AtomId]) -> Option<AromaticSystemId> {
        let query: Vec<NodeId> = atoms.iter().map(|&atom| NodeId::from(atom)).collect();
        self.0
            .find_by_participants(&query)
            .map(AromaticSystemId::from)
    }

    /// Removed by S4e, which moves frame alignment onto `reframe_to`.
    pub(crate) fn participant_permutation(
        &self,
        id: AromaticSystemId,
        atoms: &[AtomId],
    ) -> Option<Vec<ParticipantPosition>> {
        let query: Vec<NodeId> = atoms.iter().map(|&atom| NodeId::from(atom)).collect();
        self.0
            .participant_permutation(RelationId::from(id), &query)
    }
}

impl From<VarRelationSet<NodeId, Unordered, AromaticSystemForm>> for AromaticSystems {
    fn from(set: VarRelationSet<NodeId, Unordered, AromaticSystemForm>) -> Self {
        Self(Arc::new(set))
    }
}

impl From<Arc<VarRelationSet<NodeId, Unordered, AromaticSystemForm>>> for AromaticSystems {
    fn from(set: Arc<VarRelationSet<NodeId, Unordered, AromaticSystemForm>>) -> Self {
        Self(set)
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Normalize, Lattice)]
pub struct AromaticSystemForm {
    pub electrons: ElectronCountsForm,
    pub charge: NumForm,
    pub unpaired_electrons: UnpairedElectronsForm,
    pub constraints: AromaticSystemConstraintsForm,
}

/// Attribute update for an aromatic system. Ordinary fields are optional,
/// unpaired-electron components are updated independently, and undetermined constraints remove
/// their key.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AromaticSystemUpdate {
    pub electrons: Option<ElectronCountsForm>,
    pub charge: Option<NumForm>,
    pub unpaired_electrons: UnpairedElectronsUpdate,
    pub constraints: AromaticSystemConstraintsForm,
}

impl From<&str> for AromaticSystemForm {
    fn from(s: &str) -> Self {
        s.parse().expect("invalid aromatic system string")
    }
}

impl RelationData for AromaticSystemForm {
    /// The per-member electron counts are positional, so they follow a participant reorder.
    fn on_permutation(&mut self, order: &[ParticipantPosition]) {
        self.electrons.permute(order);
    }

    fn is_permutation_invariant(&self) -> bool {
        self.electrons.is_undetermined()
    }
}

impl AromaticSystemForm {
    /// Concrete: every inherent field is ground; the constraint channel does
    /// not bear on concreteness.
    pub fn is_concrete(&self) -> bool {
        let AromaticSystemForm {
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
    /// kind (last-wins per `AromaticSystemConstraintsForm::set`). Chainable.
    pub fn with_constraint(mut self, constraint: impl Into<AromaticSystemConstraintForm>) -> Self {
        self.constraints.set(constraint.into());
        self
    }

    /// Add each constraint from the iterator, replacing any existing entry
    /// of the same kind (last-wins per `AromaticSystemConstraintsForm::set`).
    /// Does not clear existing constraints; use `system.constraints.clear()`
    /// or direct field assignment for wipe-and-replace.
    pub fn with_constraints<I>(mut self, constraints: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<AromaticSystemConstraintForm>,
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
    pub fn update(&self, update: &AromaticSystemUpdate) -> AromaticSystemForm {
        let mut constraints = self.constraints.clone();
        constraints.update(&update.constraints);
        AromaticSystemForm {
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
    pub fn difference_to(&self, other: &Self) -> AromaticSystemUpdate {
        let mut constraints = AromaticSystemConstraintsForm::new();
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
        AromaticSystemUpdate {
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
    #[case::literal(AromaticSystemForm::new(ElectronCountsForm::Lit(vec![1; 3])),
        AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1; 3]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: AromaticSystemConstraintsForm::new() })]
    fn test_aromatic_system_form_new(
        #[case] actual: AromaticSystemForm,
        #[case] expected: AromaticSystemForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::literal(AromaticSystemForm::from_electrons(vec![1, 1, 1]),
        AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1; 3]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: AromaticSystemConstraintsForm::new() })]
    fn test_aromatic_system_form_from_electrons(
        #[case] actual: AromaticSystemForm,
        #[case] expected: AromaticSystemForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::with_charge(AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_charge(-1),
        AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]),
            charge: NumForm::Lit(-1), unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: AromaticSystemConstraintsForm::new() })]
    #[case::with_unpaired_electrons(AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((0_u8, 1_u8)),
        AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::closed_shell(),
            constraints: AromaticSystemConstraintsForm::new() })]
    #[case::with_constraint(
        AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_constraint(AromaticSystemConstraintForm::electron_count(6)),
        AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)) })]
    #[case::with_constraint_replaces_same_kind(
        AromaticSystemForm::default()
            .with_constraint(AromaticSystemConstraintForm::electron_count(2))
            .with_constraint(AromaticSystemConstraintForm::electron_count(6)),
        AromaticSystemForm { electrons: ElectronCountsForm::Undetermined,
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)) })]
    fn test_aromatic_system_form_with_methods(
        #[case] actual: AromaticSystemForm,
        #[case] expected: AromaticSystemForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::from_ground_electrons(AromaticSystemForm::from_electrons(vec![1; 6]).into_concrete(),
        AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1; 6]), charge: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
        constraints: AromaticSystemConstraintsForm::new() })]
    #[case::preserves_set_charge(AromaticSystemForm::from_electrons(vec![1; 6]).with_charge(1_i64).into_concrete(),
        AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1; 6]), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
        constraints: AromaticSystemConstraintsForm::new() })]
    #[case::preserves_constraints(AromaticSystemForm::from_electrons(vec![1; 6]).with_constraint(AromaticSystemConstraintForm::electron_count(6)).into_concrete(),
        AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1; 6]), charge: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
        constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)) })]
    fn test_aromatic_system_form_into_concrete(
        #[case] actual: AromaticSystemForm,
        #[case] expected: AromaticSystemForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electrons(AromaticSystemForm::from_electrons(vec![1, 1, 1]), AromaticSystemUpdate { electrons: Some(ElectronCountsForm::Lit(vec![2, 2, 2])), ..Default::default() }, AromaticSystemForm::from_electrons(vec![2, 2, 2]))]
    #[case::electrons_undetermined(AromaticSystemForm::from_electrons(vec![1, 1, 1]), AromaticSystemUpdate { electrons: Some(ElectronCountsForm::Undetermined), ..Default::default() }, AromaticSystemForm::default())]
    #[case::charge(AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_charge(0_i64), AromaticSystemUpdate { charge: Some(NumForm::Lit(-1)), ..Default::default() }, AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_charge(-1_i64))]
    #[case::charge_undetermined(AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_charge(-1_i64), AromaticSystemUpdate { charge: Some(NumForm::Undetermined), ..Default::default() }, AromaticSystemForm::from_electrons(vec![1, 1, 1]))]
    #[case::unpaired_electrons_count(AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((2_u8, 3_u8)), AromaticSystemUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: Some(NumForm::Lit(0)), multiplicity: None }, ..Default::default() }, AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((0_u8, 3_u8)))]
    #[case::unpaired_electrons_multiplicity(AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((2_u8, 3_u8)), AromaticSystemUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Lit(1)) }, ..Default::default() }, AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((2_u8, 1_u8)))]
    #[case::constraint_set(AromaticSystemForm::from_electrons(vec![1, 1, 1]), AromaticSystemUpdate { constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6_i64)), ..Default::default() }, AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_constraint(AromaticSystemConstraintForm::electron_count(6_i64)))]
    #[case::constraint_replace(AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_constraint(AromaticSystemConstraintForm::electron_count(6_i64)), AromaticSystemUpdate { constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(4_i64)), ..Default::default() }, AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_constraint(AromaticSystemConstraintForm::electron_count(4_i64)))]
    #[case::constraint_remove(AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_constraint(AromaticSystemConstraintForm::electron_count(6_i64)), AromaticSystemUpdate { constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(NumForm::Undetermined)), ..Default::default() }, AromaticSystemForm::from_electrons(vec![1, 1, 1]))]
    fn test_aromatic_system_form_update(
        #[case] system: AromaticSystemForm,
        #[case] update: AromaticSystemUpdate,
        #[case] expected: AromaticSystemForm,
    ) {
        assert_eq!(system.update(&update), expected);
    }

    #[rstest]
    #[case::empty(AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_charge(-1_i64).with_unpaired_electrons((2_u8, 3_u8)).with_constraint(AromaticSystemConstraintForm::electron_count(6_i64)))]
    fn test_aromatic_system_form_update_identity(#[case] system: AromaticSystemForm) {
        assert_eq!(system.update(&AromaticSystemUpdate::default()), system);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fields_and_constraints(
        AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_charge(0_i64).with_unpaired_electrons((2_u8, 3_u8)).with_constraint(AromaticSystemConstraintForm::electron_count(6_i64)),
        AromaticSystemForm::from_electrons(vec![2, 2, 2]).with_unpaired_electrons((2_u8, 1_u8)),
        AromaticSystemUpdate {
            electrons: Some(ElectronCountsForm::Lit(vec![2, 2, 2])),
            charge: Some(NumForm::Undetermined),
            unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Lit(1)) },
            constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(NumForm::Undetermined)),
        },
    )]
    fn test_aromatic_system_form_difference_to(
        #[case] system: AromaticSystemForm,
        #[case] other: AromaticSystemForm,
        #[case] expected: AromaticSystemUpdate,
    ) {
        assert_eq!(system.difference_to(&other), expected);
    }

    #[rstest]
    #[case::normalized(
        AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_charge(1_i64),
        AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_charge(NumForm::lit_set([1])),
    )]
    fn test_aromatic_system_form_difference_to_identity(
        #[case] system: AromaticSystemForm,
        #[case] other: AromaticSystemForm,
    ) {
        assert_eq!(
            system.difference_to(&other),
            AromaticSystemUpdate::default()
        );
    }

    #[rstest]
    #[case::three_members(
        AromaticSystemForm::from_electrons(vec![10, 20, 30]).with_charge(-1),
        vec![
            ParticipantPosition(2),
            ParticipantPosition(0),
            ParticipantPosition(1),
        ],
        AromaticSystemForm::from_electrons(vec![30, 10, 20]).with_charge(-1),
    )]
    fn test_aromatic_system_form_permute(
        #[case] mut input: AromaticSystemForm,
        #[case] order: Vec<ParticipantPosition>,
        #[case] expected: AromaticSystemForm,
    ) {
        input.permute(&order);
        assert_eq!(input, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::all_undetermined(AromaticSystemForm::default(), false)]
    #[case::charge_only(AromaticSystemForm::new(ElectronCountsForm::Undetermined).with_charge(0), false)]
    #[case::ground_no_atoms(AromaticSystemForm::new(ElectronCountsForm::Lit(Vec::new())).with_charge(0).with_unpaired_electrons((0, 1)), true)]
    #[case::all_ground_six(AromaticSystemForm::new(ElectronCountsForm::Lit(vec![1; 6])).with_charge(0).with_unpaired_electrons((0, 1)), true)]
    #[case::ground_with_constraint(AromaticSystemForm::new(ElectronCountsForm::Lit(vec![1; 6])).with_charge(0).with_unpaired_electrons((0, 1)).with_constraint(AromaticSystemConstraintForm::electron_count(6)), true)]
    fn test_aromatic_system_form_is_ground(
        #[case] form: AromaticSystemForm,
        #[case] expected: bool,
    ) {
        assert_eq!(form.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::folds_charge(
        AromaticSystemForm::default().with_charge(NumForm::lit_set([0])),
        Ok(AromaticSystemForm::default().with_charge(0)),
    )]
    #[case::charge_empty_litset_contradiction(
        AromaticSystemForm::default().with_charge(NumForm::lit_set(Vec::<i64>::new())),
        Err(Contradiction),
    )]
    fn test_aromatic_system_form_normalize(
        #[case] input: AromaticSystemForm,
        #[case] expected: Result<AromaticSystemForm, Contradiction>,
    ) {
        assert_eq!(input.normalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_matches_default(AromaticSystemForm::default(), AromaticSystemForm::default(), true)]
    #[case::default_matches_ground(AromaticSystemForm::default(), AromaticSystemForm::new(ElectronCountsForm::Lit(Vec::new())).with_charge(0).with_unpaired_electrons((0, 1)), true)]
    #[case::exact(AromaticSystemForm::new(ElectronCountsForm::Lit(vec![1; 6])).with_charge(0).with_unpaired_electrons((0, 1)),
        AromaticSystemForm::new(ElectronCountsForm::Lit(vec![1; 6])).with_charge(0).with_unpaired_electrons((0, 1)), true)]
    #[case::electrons_length_mismatch(AromaticSystemForm::new(ElectronCountsForm::Lit(vec![1; 5])),
        AromaticSystemForm::new(ElectronCountsForm::Lit(vec![1; 6])).with_charge(0).with_unpaired_electrons((0, 1)), false)]
    #[case::electrons_value_mismatch(AromaticSystemForm::new(ElectronCountsForm::Lit(vec![2; 6])),
        AromaticSystemForm::new(ElectronCountsForm::Lit(vec![1; 6])).with_charge(0).with_unpaired_electrons((0, 1)), false)]
    #[case::pattern_undetermined_electron_matches_lit(AromaticSystemForm::new(ElectronCountsForm::Undetermined),
      AromaticSystemForm::new(ElectronCountsForm::Lit(vec![1; 6])).with_charge(0).with_unpaired_electrons((0, 1)), true)]
    #[case::charge_mismatch(AromaticSystemForm::new(ElectronCountsForm::Undetermined).with_charge(1),
        AromaticSystemForm::new(ElectronCountsForm::Undetermined).with_charge(0), false)]
    #[case::unpaired_electrons_mismatch(AromaticSystemForm::new(ElectronCountsForm::Undetermined).with_unpaired_electrons((2_u8, 3_u8)),
        AromaticSystemForm::new(ElectronCountsForm::Undetermined).with_unpaired_electrons((0_u8, 1_u8)), false)]
    #[case::constraint_required_present(
        AromaticSystemForm::new(ElectronCountsForm::Undetermined).with_constraint(AromaticSystemConstraintForm::electron_count(6)),
        AromaticSystemForm::new(ElectronCountsForm::Undetermined).with_constraint(AromaticSystemConstraintForm::electron_count(6)),
        true)]
    #[case::constraint_required_absent(
        AromaticSystemForm::new(ElectronCountsForm::Undetermined).with_constraint(AromaticSystemConstraintForm::electron_count(6)),
        AromaticSystemForm::new(ElectronCountsForm::Undetermined),
        false)]
    fn test_aromatic_system_form_matches(
        #[case] pattern: AromaticSystemForm,
        #[case] target: AromaticSystemForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    #[case::both_default(
        AromaticSystemForm::default(),
        AromaticSystemForm::default(),
        Some(AromaticSystemForm::default())
    )]
    #[case::electrons_length_mismatch(
        AromaticSystemForm::from_electrons(vec![1; 6]),
        AromaticSystemForm::from_electrons(vec![1; 5]),
        None,
    )]
    #[case::narrows_electrons(
        AromaticSystemForm::new(ElectronCountsForm::Undetermined),
        AromaticSystemForm::from_electrons(vec![1; 3]),
        Some(AromaticSystemForm::from_electrons(vec![1; 3])),
    )]
    fn test_aromatic_system_form_meet(
        #[case] a: AromaticSystemForm,
        #[case] b: AromaticSystemForm,
        #[case] expected: Option<AromaticSystemForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::electrons_length_mismatch_widens_to_default(
        AromaticSystemForm::from_electrons(vec![1; 6]),
        AromaticSystemForm::from_electrons(vec![1; 5]),
        AromaticSystemForm::default(),
    )]
    fn test_aromatic_system_form_join(
        #[case] a: AromaticSystemForm,
        #[case] b: AromaticSystemForm,
        #[case] expected: AromaticSystemForm,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }
}
