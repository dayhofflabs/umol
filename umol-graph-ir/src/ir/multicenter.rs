//! Multicenter bonds: the molecule's collection and one bond's attribute form.

use std::sync::Arc;

use umol_graph_core::{NodeId, ParticipantPosition, RelationId, Remapping, VarRelationSet};
use umol_graph_ir_macros::{Lattice, Normalize};
use umol_perm::DynPermutation;

use super::constraint::{MulticenterBondConstraintForm, MulticenterBondConstraintsForm};
use super::delta::EntitySpan;
use super::electrons::ElectronCountsForm;
use super::error::Contradiction;
use super::frame::MulticenterBondsFrameAction;
use super::id::{AtomId, MulticenterBondId};
use super::num::NumForm;
use super::spin::{UnpairedElectronsForm, UnpairedElectronsUpdate};
use super::traits::{FrameTransport, Lattice, Normalize, Reframe};

/// The molecule's multicenter bonds.
///
/// The atoms bear the participant frame: the per-member electron counts of
/// [`MulticenterBondForm`] are read against it, position by position. Values are issued by checked
/// molecule construction and trusted graph-IR transformations; raw assembly is not a public
/// construction path.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct MulticenterBonds(Arc<VarRelationSet<NodeId, MulticenterBondForm>>);

impl MulticenterBonds {
    pub(crate) fn new(entries: Vec<(Vec<AtomId>, MulticenterBondForm)>) -> Self {
        Self(Arc::new(VarRelationSet::new(
            entries
                .into_iter()
                .map(|(atoms, attributes)| {
                    (atoms.into_iter().map(NodeId::from).collect(), attributes)
                })
                .collect(),
        )))
    }

    pub(crate) fn from_arc(set: Arc<VarRelationSet<NodeId, MulticenterBondForm>>) -> Self {
        Self(set)
    }

    pub fn count(&self) -> usize {
        self.0.count()
    }

    pub fn contains(&self, id: MulticenterBondId) -> bool {
        self.0.contains(RelationId::from(id))
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = MulticenterBondId> {
        self.0.ids().map(MulticenterBondId::from)
    }

    /// The atoms of `id`, in their stored frame.
    pub fn atoms(&self, id: MulticenterBondId) -> impl ExactSizeIterator<Item = AtomId> + '_ {
        self.0
            .participants(RelationId::from(id))
            .iter()
            .map(|&atom| AtomId::from(atom))
    }

    pub fn attributes(&self, id: MulticenterBondId) -> &MulticenterBondForm {
        self.0.data(RelationId::from(id))
    }

    pub(crate) fn attributes_mut(&mut self, id: MulticenterBondId) -> &mut MulticenterBondForm {
        Arc::make_mut(&mut self.0).data_mut(RelationId::from(id))
    }

    /// Ids of the multicenter bonds `atom` belongs to. Unlike aromatic systems these may overlap,
    /// so an atom can belong to several; integrity rejects only identical atom sets.
    pub fn incident_ids(
        &self,
        atom: AtomId,
    ) -> impl ExactSizeIterator<Item = MulticenterBondId> + '_ {
        self.0
            .incident(NodeId::from(atom))
            .iter()
            .map(|&id| MulticenterBondId::from(id))
    }

    pub fn has_incident(&self, atom: AtomId) -> bool {
        self.0.has_incident(NodeId::from(atom))
    }

    pub(crate) fn into_entries(self) -> Vec<(Vec<AtomId>, MulticenterBondForm)> {
        Arc::try_unwrap(self.0)
            .unwrap_or_else(|shared| (*shared).clone())
            .into_entries()
            .into_iter()
            .map(|(atoms, attributes)| (atoms.into_iter().map(AtomId::from).collect(), attributes))
            .collect()
    }

    /// The atoms of `id` as graph nodes, for graph-core interop that is not yet typed in graph-IR
    /// ids. The public accessor is [`Self::atoms`].
    pub(crate) fn atom_nodes(&self, id: MulticenterBondId) -> &[NodeId] {
        self.0.participants(RelationId::from(id))
    }

    pub(crate) fn attributes_iter_mut(
        &mut self,
    ) -> impl ExactSizeIterator<Item = &mut MulticenterBondForm> {
        Arc::make_mut(&mut self.0)
            .iter_mut()
            .map(|(_, _, attributes)| attributes)
    }

    pub(crate) fn remap(&self, remapping: &Remapping) -> Self {
        Self(Arc::new(self.0.remap(remapping)))
    }

    pub(crate) fn into_arc(self) -> Arc<VarRelationSet<NodeId, MulticenterBondForm>> {
        self.0
    }

    /// Glue `right`, relabelled into this molecule's id space, onto `self`: coinciding bonds meet,
    /// non-coinciding bonds are carried. `None` when a coincident meet is bottom.
    pub(crate) fn glue(&self, right: &Self, remapping: &Remapping) -> Option<Self> {
        self.0
            .pushout(
                &right.remap(remapping).0,
                // Multicenter bonds anchor on their atoms: the node index.
                |set, atoms| atoms.first().and_then(|&node| set.coincident(node, atoms)),
                |(left_atoms, left), (right_atoms, right)| {
                    let left_atoms: Vec<AtomId> =
                        left_atoms.iter().map(|&atom| AtomId::from(atom)).collect();
                    let right_atoms: Vec<AtomId> =
                        right_atoms.iter().map(|&atom| AtomId::from(atom)).collect();
                    let action = DynPermutation::between(&right_atoms, &left_atoms)?;
                    right.clone().reframe_by(&action)?.meet(left)
                },
            )
            .map(|merged| Self(Arc::new(merged.object)))
    }

    /// Whether bond `id` is the one over `atoms` — the known-id sibling of
    /// [`coincident_id`](Self::coincident_id).
    pub fn is_coincident(&self, id: MulticenterBondId, atoms: &[AtomId]) -> bool {
        let query: Vec<NodeId> = atoms.iter().map(|&atom| NodeId::from(atom)).collect();
        self.0.is_coincident(RelationId::from(id), &query)
    }

    /// Id of the entity coinciding with these participants — the one whose participants equal
    /// them as a multiset. The identity question, distinct from lookup.
    pub fn coincident_id(&self, atoms: &[AtomId]) -> Option<MulticenterBondId> {
        // Multicenter bonds anchor on their atoms, so the node index is the one to scan.
        let query: Vec<NodeId> = atoms.iter().map(|&atom| NodeId::from(atom)).collect();
        let anchor = *query.first()?;
        self.0
            .coincident(anchor, &query)
            .map(MulticenterBondId::from)
    }
}

impl Normalize for MulticenterBonds {
    fn normalize(mut self) -> Result<Self, Contradiction> {
        for attributes in self.attributes_iter_mut() {
            *attributes = attributes.clone().normalize()?;
        }
        Ok(self)
    }
}

impl FrameTransport for MulticenterBonds {
    type Action = MulticenterBondsFrameAction;

    fn reframe_by(mut self, actions: &Self::Action) -> Option<Self> {
        let set = Arc::make_mut(&mut self.0);
        for relation_id in set.ids().collect::<Vec<_>>() {
            let action = actions.action(MulticenterBondId::from(relation_id))?;
            if action.degree() != set.participants(relation_id).len() {
                return None;
            }
            *set.data_mut(relation_id) = set.data(relation_id).clone().reframe_by(action)?;
            set.permute_with(relation_id, &participant_order(action));
        }
        Some(self)
    }
}

impl Reframe for MulticenterBonds {
    fn representative_action(&self) -> Self::Action {
        let actions = self
            .ids()
            .map(|id| multicenter_bond_representative_action(self.atoms(id).collect()))
            .collect();
        MulticenterBondsFrameAction::from_vec(actions)
            .expect("every dynamic permutation is a multicenter-bond action")
    }

    fn reframe(self) -> Result<Self, Contradiction> {
        reframe_multicenter_bonds_with(self, |_, _| {})
    }
}

pub(crate) fn reframe_multicenter_bonds_with(
    mut multicenter_bonds: MulticenterBonds,
    mut visit: impl FnMut(MulticenterBondId, &DynPermutation),
) -> Result<MulticenterBonds, Contradiction> {
    let set = Arc::make_mut(&mut multicenter_bonds.0);
    for relation_id in set.ids().collect::<Vec<_>>() {
        let id = MulticenterBondId::from(relation_id);
        let stored = set
            .participants(relation_id)
            .iter()
            .map(|&atom| AtomId::from(atom))
            .collect();
        let action = multicenter_bond_representative_action(stored);
        let attributes = set.data(relation_id).clone().normalize()?;
        *set.data_mut(relation_id) = attributes
            .reframe_by(&action)
            .ok_or(Contradiction)?
            .normalize()?;
        set.permute_with(relation_id, &participant_order(&action));
        visit(id, &action);
    }
    Ok(multicenter_bonds)
}

/// The reaction span's multicenter bonds, one [`EntitySpan`] per entity against a single participant frame.
///
/// The `Molecule` peer is [`MulticenterBonds`]. The surface is deliberately duplicated rather than shared
/// through a payload parameter: a type parameter on the molecule-level aggregates would complicate
/// the primary carrier to serve this one. Values are issued by
/// [`ReactionSpan`](super::reaction_span::ReactionSpan).
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct MulticenterBondSpans(VarRelationSet<NodeId, EntitySpan<MulticenterBondForm>>);

impl MulticenterBondSpans {
    pub(crate) fn into_entries(self) -> Vec<(Vec<AtomId>, EntitySpan<MulticenterBondForm>)> {
        self.0
            .into_entries()
            .into_iter()
            .map(|(atoms, span)| (atoms.into_iter().map(AtomId::from).collect(), span))
            .collect()
    }

    pub(crate) fn new(entries: Vec<(Vec<AtomId>, EntitySpan<MulticenterBondForm>)>) -> Self {
        Self(VarRelationSet::new(
            entries
                .into_iter()
                .map(|(atoms, span)| (atoms.into_iter().map(NodeId::from).collect(), span))
                .collect(),
        ))
    }

    pub fn count(&self) -> usize {
        self.0.count()
    }

    pub fn contains(&self, id: MulticenterBondId) -> bool {
        self.0.contains(RelationId::from(id))
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = MulticenterBondId> {
        self.0.ids().map(MulticenterBondId::from)
    }

    /// The atoms of `id`, in their stored frame.
    pub fn atoms(&self, id: MulticenterBondId) -> impl ExactSizeIterator<Item = AtomId> + '_ {
        self.0
            .participants(RelationId::from(id))
            .iter()
            .map(|&atom| AtomId::from(atom))
    }

    pub fn attributes(&self, id: MulticenterBondId) -> &EntitySpan<MulticenterBondForm> {
        self.0.data(RelationId::from(id))
    }

    pub(crate) fn remap(&self, remapping: &Remapping) -> Self {
        Self(self.0.remap(remapping))
    }
}

impl Normalize for MulticenterBondSpans {
    fn normalize(mut self) -> Result<Self, Contradiction> {
        for id in self.0.ids().collect::<Vec<_>>() {
            *self.0.data_mut(id) = self.0.data(id).clone().normalize()?;
        }
        Ok(self)
    }
}

impl FrameTransport for MulticenterBondSpans {
    type Action = MulticenterBondsFrameAction;

    fn reframe_by(mut self, actions: &Self::Action) -> Option<Self> {
        for relation_id in self.0.ids().collect::<Vec<_>>() {
            let action = actions.action(MulticenterBondId::from(relation_id))?;
            if action.degree() != self.0.participants(relation_id).len() {
                return None;
            }
            *self.0.data_mut(relation_id) = self.0.data(relation_id).clone().reframe_by(action)?;
            self.0.permute_with(relation_id, &participant_order(action));
        }
        Some(self)
    }
}

impl Reframe for MulticenterBondSpans {
    fn representative_action(&self) -> Self::Action {
        let actions = self
            .ids()
            .map(|id| multicenter_bond_representative_action(self.atoms(id).collect()))
            .collect();
        MulticenterBondsFrameAction::from_vec(actions)
            .expect("every dynamic permutation is a multicenter-bond action")
    }

    fn reframe(self) -> Result<Self, Contradiction> {
        reframe_multicenter_bond_spans_with(self, |_, _| {})
    }
}

pub(crate) fn reframe_multicenter_bond_spans_with(
    mut multicenter_bonds: MulticenterBondSpans,
    mut visit: impl FnMut(MulticenterBondId, &DynPermutation),
) -> Result<MulticenterBondSpans, Contradiction> {
    for relation_id in multicenter_bonds.0.ids().collect::<Vec<_>>() {
        let id = MulticenterBondId::from(relation_id);
        let stored = multicenter_bonds
            .0
            .participants(relation_id)
            .iter()
            .map(|&atom| AtomId::from(atom))
            .collect();
        let action = multicenter_bond_representative_action(stored);
        let span = multicenter_bonds.0.data(relation_id).clone().normalize()?;
        *multicenter_bonds.0.data_mut(relation_id) =
            span.reframe_by(&action).ok_or(Contradiction)?.normalize()?;
        multicenter_bonds
            .0
            .permute_with(relation_id, &participant_order(&action));
        visit(id, &action);
    }
    Ok(multicenter_bonds)
}

pub(crate) fn multicenter_bond_representative_action(frame: Vec<AtomId>) -> DynPermutation {
    let mut image: Vec<usize> = (0..frame.len()).collect();
    image.sort_unstable_by_key(|&position| frame[position]);
    DynPermutation::try_from(image).expect("sorted positions form a permutation")
}

fn participant_order(action: &DynPermutation) -> Vec<ParticipantPosition> {
    action
        .image()
        .iter()
        .map(|&position| ParticipantPosition(position as u32))
        .collect()
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
                .is_none_or(|old| !old.normalized_eq(new))
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
            electrons: (!self.electrons.normalized_eq(&other.electrons))
                .then(|| other.electrons.clone()),
            charge: (!self.charge.normalized_eq(&other.charge)).then(|| other.charge.clone()),
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

impl FrameTransport for MulticenterBondForm {
    type Action = DynPermutation;

    fn reframe_by(self, action: &Self::Action) -> Option<Self> {
        let Self {
            electrons,
            charge,
            unpaired_electrons,
            constraints,
        } = self;
        Some(Self {
            electrons: electrons.reframe_by(action)?,
            charge,
            unpaired_electrons,
            constraints: constraints.reframe_by(action)?,
        })
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ir::error::Contradiction;
    use crate::ir::traits::Normalize;

    /// Both sides of a `Modified` span are read against one participant list, so one action carries
    /// both. Selection never consults the payload here, so the two sides cannot disagree about it.
    #[rstest]
    fn test_multicenter_bond_spans_reframe() {
        let mut spans = MulticenterBondSpans::new(vec![(
            vec![AtomId(1), AtomId(4), AtomId(7)],
            EntitySpan::Modified {
                lhs: MulticenterBondForm::from_electrons(vec![10, 20, 30]),
                rhs: MulticenterBondForm::from_electrons(vec![11, 21, 31]),
            },
        )]);
        spans.0.permute_with(
            RelationId(0),
            &[
                ParticipantPosition(2),
                ParticipantPosition(0),
                ParticipantPosition(1),
            ],
        );

        let source = spans.clone();
        let (reframed, actions) = spans
            .reframe_with_action()
            .expect("the forms are satisfiable");

        assert_eq!(
            reframed.atoms(MulticenterBondId(0)).collect::<Vec<_>>(),
            vec![AtomId(1), AtomId(4), AtomId(7)],
        );
        assert_eq!(
            reframed.attributes(MulticenterBondId(0)),
            &EntitySpan::Modified {
                lhs: MulticenterBondForm::from_electrons(vec![20, 30, 10]),
                rhs: MulticenterBondForm::from_electrons(vec![21, 31, 11]),
            },
        );
        assert_eq!(
            actions.action(MulticenterBondId(0)),
            Some(&DynPermutation::try_from(vec![1, 2, 0]).expect("expected action is valid")),
        );
        assert_eq!(source.reframe_by(&actions), Some(reframed));
    }

    #[rstest]
    fn test_multicenter_bond_spans_normalize() {
        let spans = MulticenterBondSpans::new(vec![(
            vec![AtomId(1), AtomId(2)],
            EntitySpan::Modified {
                lhs: MulticenterBondForm::default().with_charge(NumForm::lit_set([0])),
                rhs: MulticenterBondForm::default().with_charge(0),
            },
        )]);

        let normalized = spans.normalize().expect("the forms are satisfiable");

        assert_eq!(
            normalized.attributes(MulticenterBondId(0)),
            &EntitySpan::Unchanged(MulticenterBondForm::default().with_charge(0)),
        );
    }

    #[rstest]
    fn test_multicenter_bond_spans_reframe_identity() {
        let spans = MulticenterBondSpans::new(vec![(
            vec![AtomId(1), AtomId(4)],
            EntitySpan::Modified {
                lhs: MulticenterBondForm::from_electrons(vec![10, 20]),
                rhs: MulticenterBondForm::from_electrons(vec![11, 21]),
            },
        )]);
        let once = spans.reframe().expect("the forms are satisfiable");
        let twice = once.clone().reframe().expect("the forms are satisfiable");
        assert_eq!(twice, once);
    }

    /// A side that declines the frame change takes the whole span with it: one action serves both,
    /// so there is no partial result to keep. Here the rhs electron vector disagrees in length with
    /// the participant frame.
    #[rstest]
    fn test_multicenter_bond_spans_reframe_error() {
        let spans = MulticenterBondSpans::new(vec![(
            vec![AtomId(1), AtomId(4), AtomId(7)],
            EntitySpan::Modified {
                lhs: MulticenterBondForm::from_electrons(vec![10, 20, 30]),
                rhs: MulticenterBondForm::from_electrons(vec![11, 21]),
            },
        )]);
        assert_eq!(spans.reframe(), Err(Contradiction));
    }

    #[rstest]
    fn test_multicenter_bond_spans_framed_eq() {
        let mut unsorted = MulticenterBondSpans::new(vec![(
            vec![AtomId(1), AtomId(4)],
            EntitySpan::Modified {
                lhs: MulticenterBondForm::from_electrons(vec![10, 20]),
                rhs: MulticenterBondForm::from_electrons(vec![11, 21]),
            },
        )]);
        unsorted.0.permute_with(
            RelationId(0),
            &[ParticipantPosition(1), ParticipantPosition(0)],
        );
        // The permute leaves the payload where it was, so the stored frame now reads
        // atom 4 -> 10, atom 1 -> 20; the selected presentation states the same fact sorted.
        let selected = MulticenterBondSpans::new(vec![(
            vec![AtomId(1), AtomId(4)],
            EntitySpan::Modified {
                lhs: MulticenterBondForm::from_electrons(vec![20, 10]),
                rhs: MulticenterBondForm::from_electrons(vec![21, 11]),
            },
        )]);

        assert!(unsorted.framed_eq(&selected));
        assert!(unsorted != selected);
    }

    /// Storage sorts an `Unordered` factor on construction, so the stored frame is permuted first
    /// to model the frame-preserving storage S5 introduces.
    #[fixture]
    fn unsorted_bond() -> MulticenterBonds {
        let mut systems = MulticenterBonds::new(vec![(
            vec![AtomId(1), AtomId(4), AtomId(7)],
            MulticenterBondForm::from_electrons(vec![10, 20, 30]),
        )]);
        Arc::make_mut(&mut systems.0).permute_with(
            RelationId(0),
            &[
                ParticipantPosition(2),
                ParticipantPosition(0),
                ParticipantPosition(1),
            ],
        );
        systems
    }

    #[rstest]
    fn test_multicenter_bonds_reframe(unsorted_bond: MulticenterBonds) {
        assert_eq!(
            unsorted_bond
                .atoms(MulticenterBondId(0))
                .collect::<Vec<_>>(),
            vec![AtomId(7), AtomId(1), AtomId(4)],
        );

        let reframed = unsorted_bond.reframe().expect("the form is satisfiable");

        assert_eq!(
            reframed.atoms(MulticenterBondId(0)).collect::<Vec<_>>(),
            vec![AtomId(1), AtomId(4), AtomId(7)],
        );
        assert_eq!(
            reframed.attributes(MulticenterBondId(0)),
            &MulticenterBondForm::from_electrons(vec![20, 30, 10]),
        );
    }

    #[rstest]
    fn test_multicenter_bonds_reframe_identity(unsorted_bond: MulticenterBonds) {
        let once = unsorted_bond.reframe().expect("the form is satisfiable");
        let twice = once.clone().reframe().expect("the form is satisfiable");
        assert_eq!(twice, once);
    }

    #[rstest]
    fn test_multicenter_bonds_reframe_with_action(unsorted_bond: MulticenterBonds) {
        let (reframed, actions) = unsorted_bond
            .clone()
            .reframe_with_action()
            .expect("the form is satisfiable");

        let action = actions
            .action(MulticenterBondId(0))
            .expect("the dense action covers the bond");
        assert_eq!(action.image(), [1, 2, 0]);
        assert_eq!(unsorted_bond.reframe_by(&actions), Some(reframed));
    }

    #[rstest]
    fn test_reframe_multicenter_bonds_with(unsorted_bond: MulticenterBonds) {
        let mut visited = None;
        let reframed = reframe_multicenter_bonds_with(unsorted_bond.clone(), |id, action| {
            visited = Some((id, action.clone()));
        })
        .expect("the form is satisfiable");

        assert_eq!(
            visited,
            Some((
                MulticenterBondId(0),
                DynPermutation::try_from(vec![1, 2, 0]).expect("the expected action is valid"),
            )),
        );
        assert_eq!(unsorted_bond.reframe(), Ok(reframed));
    }

    #[rstest]
    fn test_multicenter_bonds_normalize() {
        let bonds = MulticenterBonds::new(vec![(
            vec![AtomId(1), AtomId(2)],
            MulticenterBondForm::default().with_charge(NumForm::lit_set([0])),
        )]);

        let normalized = bonds.normalize().expect("the form is satisfiable");

        assert_eq!(
            normalized.attributes(MulticenterBondId(0)),
            &MulticenterBondForm::default().with_charge(0),
        );
    }

    #[rstest]
    fn test_multicenter_bonds_framed_eq(unsorted_bond: MulticenterBonds) {
        let selected = MulticenterBonds::new(vec![(
            vec![AtomId(1), AtomId(4), AtomId(7)],
            MulticenterBondForm::from_electrons(vec![20, 30, 10]),
        )]);
        assert!(unsorted_bond.framed_eq(&selected));
        assert!(!unsorted_bond.eq(&selected));

        let different = MulticenterBonds::new(vec![(
            vec![AtomId(1), AtomId(4), AtomId(7)],
            MulticenterBondForm::from_electrons(vec![10, 20, 30]),
        )]);
        assert!(!unsorted_bond.framed_eq(&different));
    }

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

    #[rustfmt::skip]
    #[rstest]
    #[case::positioned(
        MulticenterBondForm::from_electrons(vec![10, 20, 30]).with_charge(-1),
        vec![2, 0, 1],
        Some(MulticenterBondForm::from_electrons(vec![30, 10, 20]).with_charge(-1)),
    )]
    #[case::positioned_degree(
        MulticenterBondForm::from_electrons(vec![10, 20]),
        vec![2, 0, 1],
        None,
    )]
    #[case::dimensionless(
        MulticenterBondForm::default(),
        vec![3, 1, 0, 2],
        Some(MulticenterBondForm::default()),
    )]
    #[case::frame_invariant_constraint(
        MulticenterBondForm::default()
            .with_constraint(MulticenterBondConstraintForm::electron_count(2)),
        vec![1, 0],
        Some(
            MulticenterBondForm::default()
                .with_constraint(MulticenterBondConstraintForm::electron_count(2)),
        ),
    )]
    fn test_multicenter_bond_form_reframe_by(
        #[case] input: MulticenterBondForm,
        #[case] image: Vec<usize>,
        #[case] expected: Option<MulticenterBondForm>,
    ) {
        let action = DynPermutation::try_from(image).expect("case is a permutation");
        assert_eq!(input.reframe_by(&action), expected);
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
