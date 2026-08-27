//! Dative bond form.

use std::sync::Arc;

use umol_graph_core::{
    BiRelationData, FixedVarBirelationSet, NodeId, Ordered, ParticipantPosition, RelationId,
    Remapping, Unordered,
};
use umol_graph_ir_macros::{Lattice, Normalize};

use super::constraint::{DativeBondConstraintForm, DativeBondConstraintsForm};
use super::id::{AtomId, DativeBondId};
use super::num::NumForm;
use super::traits::{Equiv, Lattice};

/// The site this entry is borne by.
/// The molecule's dative bonds.
///
/// The donors bear the participant frame; the acceptor is a single distinguished atom. Neither
/// factor is position-sensitive for [`DativeBondForm`], whose order is a scalar.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct DativeBonds(
    Arc<FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, DativeBondForm>>,
);

impl DativeBonds {
    pub fn new(entries: Vec<(Vec<AtomId>, AtomId, DativeBondForm)>) -> Self {
        Self(Arc::new(FixedVarBirelationSet::new(
            entries
                .into_iter()
                .map(|(donors, acceptor, attributes)| {
                    (
                        [NodeId::from(acceptor)],
                        donors.into_iter().map(NodeId::from).collect(),
                        attributes,
                    )
                })
                .collect(),
        )))
    }

    pub fn count(&self) -> usize {
        self.0.count()
    }

    pub fn contains(&self, id: DativeBondId) -> bool {
        self.0.contains(RelationId::from(id))
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = DativeBondId> {
        self.0.relation_ids().map(DativeBondId::from)
    }

    /// The atom accepting the donated pairs.
    pub fn acceptor(&self, id: DativeBondId) -> AtomId {
        AtomId::from(self.0.participants_1(RelationId::from(id))[0])
    }

    /// The donating atoms of `id`, in their stored frame.
    pub fn donors(&self, id: DativeBondId) -> impl ExactSizeIterator<Item = AtomId> + '_ {
        self.0
            .participants_2(RelationId::from(id))
            .iter()
            .map(|&atom| AtomId::from(atom))
    }

    pub fn attributes(&self, id: DativeBondId) -> &DativeBondForm {
        self.0.data(RelationId::from(id))
    }

    pub fn attributes_mut(&mut self, id: DativeBondId) -> &mut DativeBondForm {
        Arc::make_mut(&mut self.0).data_mut(RelationId::from(id))
    }

    /// Ids of the dative bonds `atom` takes part in, as acceptor or donor.
    pub fn incident_ids(&self, atom: AtomId) -> impl ExactSizeIterator<Item = DativeBondId> + '_ {
        self.0
            .incident(NodeId::from(atom))
            .iter()
            .map(|&id| DativeBondId::from(id))
    }

    pub fn has_incident(&self, atom: AtomId) -> bool {
        self.0.has_incident(NodeId::from(atom))
    }

    pub fn into_entries(self) -> Vec<(Vec<AtomId>, AtomId, DativeBondForm)> {
        Arc::try_unwrap(self.0)
            .unwrap_or_else(|shared| (*shared).clone())
            .into_entries()
            .into_iter()
            .map(|(acceptor, donors, attributes)| {
                (
                    donors.into_iter().map(AtomId::from).collect(),
                    AtomId::from(acceptor[0]),
                    attributes,
                )
            })
            .collect()
    }

    /// The acceptor of `id` as a graph node, for graph-core interop that is not yet typed in
    /// graph-IR ids. The public accessor is [`Self::acceptor`].
    pub(crate) fn acceptor_node(&self, id: DativeBondId) -> NodeId {
        self.0.participants_1(RelationId::from(id))[0]
    }

    /// The donors of `id` as graph nodes, for graph-core interop that is not yet typed in graph-IR
    /// ids. The public accessor is [`Self::donors`].
    pub(crate) fn donor_nodes(&self, id: DativeBondId) -> &[NodeId] {
        self.0.participants_2(RelationId::from(id))
    }

    pub(crate) fn attributes_iter_mut(
        &mut self,
    ) -> impl ExactSizeIterator<Item = &mut DativeBondForm> {
        Arc::make_mut(&mut self.0).data_iter_mut()
    }

    pub(crate) fn remap(&self, remapping: &Remapping) -> Self {
        Self(Arc::new(self.0.remap(remapping)))
    }

    pub(crate) fn into_arc(
        self,
    ) -> Arc<FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, DativeBondForm>> {
        self.0
    }

    /// Glue `right`, relabelled into this molecule's id space, onto `self`: coinciding bonds meet,
    /// non-coinciding bonds are carried. `None` when a coincident meet is bottom.
    pub(crate) fn glue(&self, right: &Self, remapping: &Remapping) -> Option<Self> {
        self.0
            .pushout(&right.remap(remapping).0, |a, b| a.meet(b))
            .map(|merged| Self(Arc::new(merged.object)))
    }

    /// Removed by the lookup relocation: replaced by `of_id` keyed on the uniqueness key, which for
    /// a dative bond is its acceptor and donor set.
    pub(crate) fn find_by_participants(
        &self,
        acceptor: AtomId,
        donors: &[AtomId],
    ) -> Option<DativeBondId> {
        let donors: Vec<NodeId> = donors.iter().map(|&atom| NodeId::from(atom)).collect();
        self.0
            .find_by_participants(&[NodeId::from(acceptor)], &donors)
            .map(DativeBondId::from)
    }

    /// Removed by S4e, which moves frame alignment onto `reframe_to`.
    #[allow(clippy::type_complexity)]
    pub(crate) fn participant_permutation(
        &self,
        id: DativeBondId,
        acceptor: &[AtomId],
        donors: &[AtomId],
    ) -> Option<(Vec<ParticipantPosition>, Vec<ParticipantPosition>)> {
        let acceptor: Vec<NodeId> = acceptor.iter().map(|&a| NodeId::from(a)).collect();
        let donors: Vec<NodeId> = donors.iter().map(|&a| NodeId::from(a)).collect();
        self.0
            .participant_permutation(RelationId::from(id), &acceptor, &donors)
    }
}

impl From<FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, DativeBondForm>>
    for DativeBonds
{
    fn from(
        set: FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, DativeBondForm>,
    ) -> Self {
        Self(Arc::new(set))
    }
}

impl From<Arc<FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, DativeBondForm>>>
    for DativeBonds
{
    fn from(
        set: Arc<FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, DativeBondForm>>,
    ) -> Self {
        Self(set)
    }
}

/// Dative bond data: bond order (number of electron pairs donated) and
/// constraints. The acceptor and donor atoms are the participants of the
/// owning molecule's dative relation; this value holds no participant refs.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Normalize, Lattice)]
pub struct DativeBondForm {
    /// Bond order — number of electron pairs donated.
    pub order: NumForm,
    pub constraints: DativeBondConstraintsForm,
}

/// Attribute update for a dative bond. The order is optional, and
/// undetermined constraints remove their key.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DativeBondUpdate {
    pub order: Option<NumForm>,
    pub constraints: DativeBondConstraintsForm,
}

impl DativeBondForm {
    /// Concrete: every inherent field is ground; the constraint channel does
    /// not bear on concreteness.
    pub fn is_concrete(&self) -> bool {
        let DativeBondForm {
            order,
            constraints: _,
        } = self;
        order.is_ground()
    }
    pub fn new(order: NumForm) -> Self {
        Self {
            order,
            constraints: DativeBondConstraintsForm::new(),
        }
    }

    pub fn from_order(order: u8) -> Self {
        Self::new(NumForm::Lit(order as i64))
    }

    pub fn with_order(mut self, order: impl Into<NumForm>) -> Self {
        self.order = order.into();
        self
    }

    /// Add a single constraint, replacing any existing entry of the same
    /// key (last-wins per `DativeBondConstraintsForm::set`). Chainable.
    pub fn with_constraint(mut self, constraint: impl Into<DativeBondConstraintForm>) -> Self {
        self.constraints.set(constraint.into());
        self
    }

    /// Add each constraint from the iterator, replacing any existing entry
    /// of the same key (last-wins per `DativeBondConstraintsForm::set`). Does
    /// not clear existing constraints; use `bond.constraints.clear()` or
    /// direct field assignment for wipe-and-replace.
    pub fn with_constraints<I>(mut self, constraints: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<DativeBondConstraintForm>,
    {
        self.constraints
            .extend(constraints.into_iter().map(Into::into));
        self
    }

    /// No-op: `DativeBondForm` has no value-bearing fields besides `order`,
    /// which is essential and never filled. Provided for API symmetry.
    pub fn into_concrete(self) -> Self {
        self
    }

    /// Apply an attribute update, leaving omitted leaves and constraint keys unchanged.
    pub fn update(&self, update: &DativeBondUpdate) -> DativeBondForm {
        let mut constraints = self.constraints.clone();
        constraints.update(&update.constraints);
        DativeBondForm {
            order: update.order.clone().unwrap_or_else(|| self.order.clone()),
            constraints,
        }
    }

    /// Derive the minimal normalized attribute update carrying `self` to `other`.
    /// Frame-invariant: no field is position-indexed, so a frame change carries the form
    /// unchanged.
    ///
    /// Destructured exhaustively on purpose: a new position-indexed field must fail to compile
    /// here rather than be silently left in the old frame.
    pub fn reframe_to(self, _from: &[AtomId], _to: &[AtomId]) -> Option<Self> {
        let Self { order, constraints } = self;
        Some(Self { order, constraints })
    }

    pub fn difference_to(&self, other: &Self) -> DativeBondUpdate {
        let mut constraints = DativeBondConstraintsForm::new();
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
        DativeBondUpdate {
            order: (!self.order.equiv(&other.order)).then(|| other.order.clone()),
            constraints,
        }
    }
}

impl From<&str> for DativeBondForm {
    fn from(s: &str) -> Self {
        s.parse().expect("invalid dative bond string")
    }
}

impl BiRelationData for DativeBondForm {
    /// `order` is a scalar; neither the acceptor nor the donor factor is position-indexed.
    fn on_permutation(
        &mut self,
        _order_1: &[ParticipantPosition],
        _order_2: &[ParticipantPosition],
    ) {
    }

    fn is_permutation_invariant(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ir::boolean::BooleanForm;
    use crate::ir::constraint::RingScope;
    use crate::ir::error::Contradiction;
    use crate::ir::traits::{Lattice, Normalize};

    #[rustfmt::skip]
    #[rstest]
    #[case::new(DativeBondForm::new(NumForm::Lit(2)),
        DativeBondForm { order: NumForm::Lit(2), constraints: DativeBondConstraintsForm::new() })]
    #[case::from_order(DativeBondForm::from_order(3),
        DativeBondForm { order: NumForm::Lit(3), constraints: DativeBondConstraintsForm::new() })]
    fn test_dative_bond_form_new(#[case] actual: DativeBondForm, #[case] expected: DativeBondForm) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::with_order(DativeBondForm::default().with_order(2),
        DativeBondForm { order: NumForm::Lit(2), constraints: DativeBondConstraintsForm::new() })]
    #[case::with_constraint(DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        DativeBondForm { order: NumForm::Lit(1),
            constraints: DativeBondConstraintsForm::from(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))) })]
    #[case::with_constraints_extends(
        DativeBondForm::from_order(1)
            .with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)))
            .with_constraints([DativeBondConstraintForm::ring_membership(RingScope::All, 1), DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1)]),
        DativeBondForm { order: NumForm::Lit(1),
            constraints: DativeBondConstraintsForm::from_iter([
                DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)),
                DativeBondConstraintForm::ring_membership(RingScope::All, 1),
                DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1),
            ]) })]
    #[case::with_constraint_appends_different_scopes(
        DativeBondForm::from_order(1)
            .with_constraint(DativeBondConstraintForm::ring_membership(RingScope::Size(5), 1))
            .with_constraint(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1)),
        DativeBondForm { order: NumForm::Lit(1),
            constraints: DativeBondConstraintsForm::from_iter([
                DativeBondConstraintForm::ring_membership(RingScope::Size(5), 1),
                DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1),
            ]) })]
    fn test_dative_bond_form_with_methods(#[case] actual: DativeBondForm, #[case] expected: DativeBondForm) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::from_order(DativeBondForm::from_order(1))]
    #[case::with_constraint(DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))))]
    fn test_dative_bond_form_into_concrete(#[case] bond: DativeBondForm) {
        assert_eq!(bond.clone().into_concrete(), bond);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::order(DativeBondForm::from_order(1), DativeBondUpdate { order: Some(NumForm::Lit(2)), ..Default::default() }, DativeBondForm::from_order(2))]
    #[case::order_undetermined(DativeBondForm::from_order(1), DativeBondUpdate { order: Some(NumForm::Undetermined), ..Default::default() }, DativeBondForm::default())]
    #[case::constraint_set(DativeBondForm::from_order(1), DativeBondUpdate { constraints: DativeBondConstraintsForm::from(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))), ..Default::default() }, DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))))]
    #[case::constraint_replace(DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1_i64)), DativeBondUpdate { constraints: DativeBondConstraintsForm::from(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 2_i64)), ..Default::default() }, DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 2_i64)))]
    #[case::constraint_remove(DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1_i64)), DativeBondUpdate { constraints: DativeBondConstraintsForm::from(DativeBondConstraintForm::ring_membership(RingScope::Size(6), NumForm::Undetermined)), ..Default::default() }, DativeBondForm::from_order(1))]
    fn test_dative_bond_form_update(
        #[case] bond: DativeBondForm,
        #[case] update: DativeBondUpdate,
        #[case] expected: DativeBondForm,
    ) {
        assert_eq!(bond.update(&update), expected);
    }

    #[rstest]
    #[case::empty(DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))))]
    fn test_dative_bond_form_update_identity(#[case] bond: DativeBondForm) {
        assert_eq!(bond.update(&DativeBondUpdate::default()), bond);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ground(DativeBondForm::from_order(1))]
    #[case::with_constraint(DativeBondForm::from_order(1)
        .with_constraint(DativeBondConstraintForm::aromatic(true)))]
    #[case::undetermined(DativeBondForm::default())]
    fn test_dative_bond_form_reframe_to(#[case] input: DativeBondForm) {
        assert_eq!(
            input.clone().reframe_to(&[AtomId(4), AtomId(7)], &[AtomId(7), AtomId(4)]),
            Some(input),
        );
    }

    #[rstest]
    #[case::fields_and_constraints(
        DativeBondForm::from_order(1).with_constraints([
            DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1_i64),
        ]),
        DativeBondForm::from_order(2).with_constraints([
            DativeBondConstraintForm::Aromatic(BooleanForm::Lit(false)),
            DativeBondConstraintForm::ring_membership(RingScope::All, 2_i64),
        ]),
        DativeBondUpdate {
            order: Some(NumForm::Lit(2)),
            constraints: DativeBondConstraintsForm::from_iter([
                DativeBondConstraintForm::Aromatic(BooleanForm::Lit(false)),
                DativeBondConstraintForm::ring_membership(RingScope::All, 2_i64),
                DativeBondConstraintForm::ring_membership(RingScope::Size(6), NumForm::Undetermined),
            ]),
        },
    )]
    fn test_dative_bond_form_difference_to(
        #[case] bond: DativeBondForm,
        #[case] other: DativeBondForm,
        #[case] expected: DativeBondUpdate,
    ) {
        assert_eq!(bond.difference_to(&other), expected);
    }

    #[rstest]
    #[case::normalized(
        DativeBondForm::from_order(1),
        DativeBondForm::new(NumForm::lit_set([1])),
    )]
    fn test_dative_bond_form_difference_to_identity(
        #[case] bond: DativeBondForm,
        #[case] other: DativeBondForm,
    ) {
        assert_eq!(bond.difference_to(&other), DativeBondUpdate::default());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_(DativeBondForm::default(), false)]
    #[case::order_lit(DativeBondForm::from_order(1), true)]
    #[case::order_undetermined(DativeBondForm::new(NumForm::Undetermined), false)]
    #[case::ground_with_constraint(DativeBondForm { order: NumForm::Lit(1),
        constraints: DativeBondConstraintsForm::from(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1)) }, true)]
    fn test_dative_bond_form_is_ground(#[case] form: DativeBondForm, #[case] expected: bool) {
        assert_eq!(form.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::folds_order(
        DativeBondForm::new(NumForm::lit_set([1])),
        Ok(DativeBondForm::from_order(1)),
    )]
    #[case::order_empty_litset_contradiction(
        DativeBondForm::new(NumForm::lit_set(Vec::<i64>::new())),
        Err(Contradiction),
    )]
    fn test_dative_bond_form_normalize(
        #[case] input: DativeBondForm,
        #[case] expected: Result<DativeBondForm, Contradiction>,
    ) {
        assert_eq!(input.normalize(), expected);
    }

    #[rstest]
    #[case::same_order(DativeBondForm::from_order(1), DativeBondForm::from_order(1), true)]
    #[case::order_mismatch(DativeBondForm::from_order(1), DativeBondForm::from_order(2), false)]
    #[case::pattern_wildcard_matches(
        DativeBondForm::new(NumForm::Undetermined),
        DativeBondForm::from_order(2),
        true
    )]
    #[case::pattern_more_specific(
        DativeBondForm::from_order(2),
        DativeBondForm::new(NumForm::Undetermined),
        false
    )]
    #[case::constraint_required_present(
        DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        true
    )]
    #[case::constraint_required_absent(
        DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        DativeBondForm::from_order(1),
        false
    )]
    fn test_dative_bond_form_matches(
        #[case] pattern: DativeBondForm,
        #[case] target: DativeBondForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    #[case::both_default(
        DativeBondForm::default(),
        DativeBondForm::default(),
        Some(DativeBondForm::default())
    )]
    #[case::order_mismatch(DativeBondForm::from_order(1), DativeBondForm::from_order(2), None)]
    #[case::narrows_order(
        DativeBondForm { order: NumForm::Undetermined, constraints: DativeBondConstraintsForm::new() },
        DativeBondForm { order: NumForm::Lit(1), constraints: DativeBondConstraintsForm::new() },
        Some(DativeBondForm { order: NumForm::Lit(1), constraints: DativeBondConstraintsForm::new() }),
    )]
    fn test_dative_bond_form_meet(
        #[case] a: DativeBondForm,
        #[case] b: DativeBondForm,
        #[case] expected: Option<DativeBondForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::order_widens_to_lit_set(
        DativeBondForm::from_order(1),
        DativeBondForm::from_order(2),
        DativeBondForm { order: NumForm::lit_set([1, 2]), constraints: DativeBondConstraintsForm::new() },
    )]
    fn test_dative_bond_form_join(
        #[case] a: DativeBondForm,
        #[case] b: DativeBondForm,
        #[case] expected: DativeBondForm,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }
}
