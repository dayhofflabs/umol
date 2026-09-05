//! Stereochemistry forms: the configuration value and the operator-expression
//! tree over it.
//!
//! A configuration value is a dense coset index per stereo kind, corresponds to OpenSMILES
//! numbering for SP, TB, and OH.
//! `~` and `^` are group actions on the index; the owning configuration's
//! `normalize` folds closed operator-expressions against the coset algebra.

use std::borrow::Cow;
use std::collections::BTreeSet;
use std::sync::Arc;

use strum::VariantArray;
use umol_graph_core::{
    EdgeId, FixedVarBirelationSet, GraphCorrespondence, GraphRemapping, NodeId,
    ParticipantPosition, RelationId,
};
use umol_graph_ir_macros::{Lattice, Normalize};
use umol_perm::{ClassKey, Permutation};

use super::constraint::{
    FluxionalityForm, LigandPermutation, LigandSymmetryForm, OrientedLigandPermutation,
    StereoAtomConstraintForm, StereoAtomConstraintsForm, StereoBondConstraintForm,
    StereoBondConstraintsForm, StereoLigandPair, TopicityForm,
};
use super::delta::EntitySpan;
use super::error::{Contradiction, NoJoin};
use super::frame::{StereoAtomsFrameAction, StereoBondsFrameAction};
use super::id::{AtomId, BondId, StereoAtomId, StereoBondId};
use super::ligand::StereoLigand;
use super::traits::{AsLit, FrameTransport, Lattice, Normalize, Reframe};

/// The molecule's stereo atoms. The ligands bear the frame the configuration is read against; the
/// site is an atom.
///
/// Owns the frame structure its storage shape cannot state: which factor bears the participant
/// frame, and which is a site. Values are issued by checked molecule construction and trusted
/// graph-IR transformations; raw assembly is not a public construction path.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct StereoAtoms(Arc<FixedVarBirelationSet<NodeId, 1, StereoLigand, StereoAtomForm>>);

impl StereoAtoms {
    pub(crate) fn new(entries: Vec<(AtomId, Vec<StereoLigand>, StereoAtomForm)>) -> Self {
        Self(Arc::new(FixedVarBirelationSet::new(
            entries
                .into_iter()
                .map(|(site, ligands, attributes)| ([NodeId::from(site)], ligands, attributes))
                .collect(),
        )))
    }

    pub(crate) fn from_arc(
        set: Arc<FixedVarBirelationSet<NodeId, 1, StereoLigand, StereoAtomForm>>,
    ) -> Self {
        Self(set)
    }

    pub fn count(&self) -> usize {
        self.0.count()
    }

    pub fn contains(&self, id: StereoAtomId) -> bool {
        self.0.contains(RelationId::from(id))
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = StereoAtomId> {
        self.0.ids().map(StereoAtomId::from)
    }

    /// The atom `id` is borne by.
    pub fn site(&self, id: StereoAtomId) -> AtomId {
        AtomId::from(self.0.participants_1(RelationId::from(id))[0])
    }

    /// The ligands of `id`, in the frame its configuration is read against.
    pub fn ligands(&self, id: StereoAtomId) -> &[StereoLigand] {
        self.0.participants_2(RelationId::from(id))
    }

    pub fn attributes(&self, id: StereoAtomId) -> &StereoAtomForm {
        self.0.data(RelationId::from(id))
    }

    pub(crate) fn attributes_mut(&mut self, id: StereoAtomId) -> &mut StereoAtomForm {
        Arc::make_mut(&mut self.0).data_mut(RelationId::from(id))
    }

    /// Ids of the stereo atoms `atom` takes part in, as site or as ligand.
    pub fn incident_ids(&self, atom: AtomId) -> impl ExactSizeIterator<Item = StereoAtomId> + '_ {
        self.0
            .incident(NodeId::from(atom))
            .iter()
            .map(|&id| StereoAtomId::from(id))
    }

    pub fn has_incident(&self, atom: AtomId) -> bool {
        self.0.has_incident(NodeId::from(atom))
    }

    pub(crate) fn into_entries(self) -> Vec<(AtomId, Vec<StereoLigand>, StereoAtomForm)> {
        Arc::try_unwrap(self.0)
            .unwrap_or_else(|shared| (*shared).clone())
            .into_entries()
            .into_iter()
            .map(|(site, ligands, attributes)| (AtomId::from(site[0]), ligands, attributes))
            .collect()
    }

    pub(crate) fn attributes_iter_mut(
        &mut self,
    ) -> impl ExactSizeIterator<Item = &mut StereoAtomForm> {
        Arc::make_mut(&mut self.0)
            .iter_mut()
            .map(|(_, _, _, attributes)| attributes)
    }

    /// Map participant references, preserving entity ids, row order, attributes, and frames.
    ///
    /// # Panics
    /// Panics if a referenced node or edge has no image in `correspondence`.
    pub fn map(&self, correspondence: &GraphCorrespondence) -> Self {
        self.try_map(correspondence)
            .expect("correspondence must cover every participant reference")
    }

    /// Map participant references, or return `None` if any reference has no image.
    /// Unreferenced nodes and edges need not have images. No entity is dropped.
    pub fn try_map(&self, correspondence: &GraphCorrespondence) -> Option<Self> {
        Some(Self(Arc::new(self.0.try_map(correspondence)?)))
    }

    pub(crate) fn remap(&self, remapping: &GraphRemapping) -> Self {
        Self(Arc::new(self.0.remap(remapping)))
    }

    pub(crate) fn into_arc(
        self,
    ) -> Arc<FixedVarBirelationSet<NodeId, 1, StereoLigand, StereoAtomForm>> {
        self.0
    }

    /// Glue `right`, relabelled into this molecule's id space, onto `self`: a coinciding entry has
    /// its configuration reframed into the retained `self` ligand frame and then meets, a
    /// non-coinciding entry is carried in its own frame. `None` when a reframing is inadmissible or
    /// a coincident meet is bottom.
    pub(crate) fn glue(&self, right: &Self, remapping: &GraphRemapping) -> Option<Self> {
        self.0
            .pushout(
                &right.remap(remapping).0,
                // A stereo atom's site is an atom, unique by integrity: the sharpest node anchor.
                |set, site, ligands| {
                    site.first()
                        .and_then(|&node| set.coincident(node, site, ligands))
                },
                |(_, left_ligands, left), (_, right_ligands, right)| {
                    let action = Permutation::between(right_ligands, left_ligands)?;
                    right.clone().reframe_by(&action)?.meet(left)
                },
            )
            .map(|merged| Self(Arc::new(merged.object)))
    }

    /// Whether stereo atom `id` is the one on `site` over `ligands` — the known-id sibling of
    /// [`coincident_id`](Self::coincident_id).
    pub fn is_coincident(&self, id: StereoAtomId, site: AtomId, ligands: &[StereoLigand]) -> bool {
        self.0
            .is_coincident(RelationId::from(id), &[NodeId::from(site)], ligands)
    }

    /// Id of the entity coinciding with these participants — the one whose participants equal
    /// them as a multiset. The identity question, distinct from lookup.
    pub fn coincident_id(&self, site: AtomId, ligands: &[StereoLigand]) -> Option<StereoAtomId> {
        // A stereo atom's site is an atom, and integrity makes it unique, so it is the sharpest
        // node anchor available.
        self.0
            .coincident(NodeId::from(site), &[NodeId::from(site)], ligands)
            .map(StereoAtomId::from)
    }
}

pub(crate) fn stereo_atom_representative_action(frame: &[StereoLigand]) -> Option<Permutation> {
    let mut image: Vec<usize> = (0..frame.len()).collect();
    image.sort_by(|&left, &right| frame[left].cmp(&frame[right]).then(left.cmp(&right)));
    Permutation::try_from(image.as_slice()).ok()
}

pub(crate) fn stereo_bond_representative_action(frame: &[StereoLigand]) -> Option<Permutation> {
    ClassKey::CisTrans.space().normalizer(frame)
}

fn participant_order(action: Permutation) -> Vec<ParticipantPosition> {
    (0..action.degree())
        .map(|position| ParticipantPosition(action.apply(position) as u32))
        .collect()
}

impl Normalize for StereoAtoms {
    fn normalize(mut self) -> Result<Self, Contradiction> {
        for attributes in self.attributes_iter_mut() {
            *attributes = attributes.clone().normalize()?;
        }
        Ok(self)
    }
}

impl FrameTransport for StereoAtoms {
    type Action = StereoAtomsFrameAction;

    fn reframe_by(mut self, actions: &Self::Action) -> Option<Self> {
        let set = Arc::make_mut(&mut self.0);
        for relation_id in set.ids().collect::<Vec<_>>() {
            let action = actions.action(StereoAtomId::from(relation_id))?;
            if action.degree() != set.participants_2(relation_id).len() {
                return None;
            }
            *set.data_mut(relation_id) = set.data(relation_id).clone().reframe_by(action)?;
            set.permute_2_with(relation_id, &participant_order(*action));
        }
        Some(self)
    }
}

impl Reframe for StereoAtoms {
    fn representative_action(&self) -> Self::Action {
        let actions = self
            .ids()
            .map(|id| {
                stereo_atom_representative_action(self.ligands(id))
                    .expect("integrity-valid stereo-atom frames fit the bounded action")
            })
            .collect();
        StereoAtomsFrameAction::from_vec(actions)
            .expect("every bounded permutation is a stereo-atom action")
    }

    fn reframe(self) -> Result<Self, Contradiction> {
        reframe_stereo_atoms_with(self, |_, _| {})
    }
}

pub(crate) fn reframe_stereo_atoms_with(
    mut stereo_atoms: StereoAtoms,
    mut visit: impl FnMut(StereoAtomId, Permutation),
) -> Result<StereoAtoms, Contradiction> {
    let set = Arc::make_mut(&mut stereo_atoms.0);
    for relation_id in set.ids().collect::<Vec<_>>() {
        let id = StereoAtomId::from(relation_id);
        let action = stereo_atom_representative_action(set.participants_2(relation_id))
            .ok_or(Contradiction)?;
        let attributes = set.data(relation_id).clone().normalize()?;
        *set.data_mut(relation_id) = attributes
            .reframe_by(&action)
            .ok_or(Contradiction)?
            .normalize()?;
        set.permute_2_with(relation_id, &participant_order(action));
        visit(id, action);
    }
    Ok(stereo_atoms)
}

/// The molecule's stereo bonds. The ligands bear the frame the configuration is read against; the
/// site is a bond.
///
/// Owns the frame structure its storage shape cannot state: which factor bears the participant
/// frame, and which is a site. Values are issued by checked molecule construction and trusted
/// graph-IR transformations; raw assembly is not a public construction path.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct StereoBonds(Arc<FixedVarBirelationSet<EdgeId, 1, StereoLigand, StereoBondForm>>);

impl StereoBonds {
    pub(crate) fn new(entries: Vec<(BondId, Vec<StereoLigand>, StereoBondForm)>) -> Self {
        Self(Arc::new(FixedVarBirelationSet::new(
            entries
                .into_iter()
                .map(|(site, ligands, attributes)| ([EdgeId::from(site)], ligands, attributes))
                .collect(),
        )))
    }

    pub(crate) fn from_arc(
        set: Arc<FixedVarBirelationSet<EdgeId, 1, StereoLigand, StereoBondForm>>,
    ) -> Self {
        Self(set)
    }

    pub fn count(&self) -> usize {
        self.0.count()
    }

    pub fn contains(&self, id: StereoBondId) -> bool {
        self.0.contains(RelationId::from(id))
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = StereoBondId> {
        self.0.ids().map(StereoBondId::from)
    }

    /// The bond `id` is borne by.
    pub fn site(&self, id: StereoBondId) -> BondId {
        BondId::from(self.0.participants_1(RelationId::from(id))[0])
    }

    /// The ligands of `id`, in the frame its configuration is read against.
    pub fn ligands(&self, id: StereoBondId) -> &[StereoLigand] {
        self.0.participants_2(RelationId::from(id))
    }

    pub fn attributes(&self, id: StereoBondId) -> &StereoBondForm {
        self.0.data(RelationId::from(id))
    }

    pub(crate) fn attributes_mut(&mut self, id: StereoBondId) -> &mut StereoBondForm {
        Arc::make_mut(&mut self.0).data_mut(RelationId::from(id))
    }

    /// Ids of the stereo atoms `atom` takes part in, as site or as ligand.
    pub fn incident_ids(&self, atom: AtomId) -> impl ExactSizeIterator<Item = StereoBondId> + '_ {
        self.0
            .incident(NodeId::from(atom))
            .iter()
            .map(|&id| StereoBondId::from(id))
    }

    pub fn has_incident(&self, atom: AtomId) -> bool {
        self.0.has_incident(NodeId::from(atom))
    }

    /// Ids of the stereo bonds `bond` is the site of.
    pub fn incident_bond_ids(
        &self,
        bond: BondId,
    ) -> impl ExactSizeIterator<Item = StereoBondId> + '_ {
        self.0
            .incident_edge(EdgeId::from(bond))
            .iter()
            .map(|&id| StereoBondId::from(id))
    }

    pub fn has_incident_bond(&self, bond: BondId) -> bool {
        self.0.has_incident_edge(EdgeId::from(bond))
    }

    pub(crate) fn into_entries(self) -> Vec<(BondId, Vec<StereoLigand>, StereoBondForm)> {
        Arc::try_unwrap(self.0)
            .unwrap_or_else(|shared| (*shared).clone())
            .into_entries()
            .into_iter()
            .map(|(site, ligands, attributes)| (BondId::from(site[0]), ligands, attributes))
            .collect()
    }

    pub(crate) fn attributes_iter_mut(
        &mut self,
    ) -> impl ExactSizeIterator<Item = &mut StereoBondForm> {
        Arc::make_mut(&mut self.0)
            .iter_mut()
            .map(|(_, _, _, attributes)| attributes)
    }

    /// Map participant references, preserving entity ids, row order, attributes, and frames.
    ///
    /// # Panics
    /// Panics if a referenced node or edge has no image in `correspondence`.
    pub fn map(&self, correspondence: &GraphCorrespondence) -> Self {
        self.try_map(correspondence)
            .expect("correspondence must cover every participant reference")
    }

    /// Map participant references, or return `None` if any reference has no image.
    /// Unreferenced nodes and edges need not have images. No entity is dropped.
    pub fn try_map(&self, correspondence: &GraphCorrespondence) -> Option<Self> {
        Some(Self(Arc::new(self.0.try_map(correspondence)?)))
    }

    pub(crate) fn remap(&self, remapping: &GraphRemapping) -> Self {
        Self(Arc::new(self.0.remap(remapping)))
    }

    pub(crate) fn into_arc(
        self,
    ) -> Arc<FixedVarBirelationSet<EdgeId, 1, StereoLigand, StereoBondForm>> {
        self.0
    }

    /// Glue `right`, relabelled into this molecule's id space, onto `self`: a coinciding entry has
    /// its configuration reframed into the retained `self` ligand frame and then meets, a
    /// non-coinciding entry is carried in its own frame. `None` when a reframing is inadmissible or
    /// a coincident meet is bottom.
    pub(crate) fn glue(&self, right: &Self, remapping: &GraphRemapping) -> Option<Self> {
        self.0
            .pushout(
                &right.remap(remapping).0,
                // A stereo bond's site is a bond: the one entity kind scanning the edge index.
                |set, site, ligands| {
                    site.first()
                        .and_then(|&edge| set.coincident_edge(edge, site, ligands))
                },
                |(_, left_ligands, left), (_, right_ligands, right)| {
                    let action = Permutation::between(right_ligands, left_ligands)?;
                    right.clone().reframe_by(&action)?.meet(left)
                },
            )
            .map(|merged| Self(Arc::new(merged.object)))
    }

    /// Whether stereo bond `id` is the one on `site` over `ligands` — the known-id sibling of
    /// [`coincident_id`](Self::coincident_id).
    pub fn is_coincident(&self, id: StereoBondId, site: BondId, ligands: &[StereoLigand]) -> bool {
        self.0
            .is_coincident(RelationId::from(id), &[EdgeId::from(site)], ligands)
    }

    /// Id of the entity coinciding with these participants — the one whose participants equal
    /// them as a multiset. The identity question, distinct from lookup.
    pub fn coincident_id(&self, site: BondId, ligands: &[StereoLigand]) -> Option<StereoBondId> {
        // A stereo bond's site is a bond, so this is the one entity kind that scans the edge index.
        self.0
            .coincident_edge(EdgeId::from(site), &[EdgeId::from(site)], ligands)
            .map(StereoBondId::from)
    }
}

impl Normalize for StereoBonds {
    fn normalize(mut self) -> Result<Self, Contradiction> {
        for attributes in self.attributes_iter_mut() {
            *attributes = attributes.clone().normalize()?;
        }
        Ok(self)
    }
}

impl FrameTransport for StereoBonds {
    type Action = StereoBondsFrameAction;

    fn reframe_by(mut self, actions: &Self::Action) -> Option<Self> {
        let set = Arc::make_mut(&mut self.0);
        for relation_id in set.ids().collect::<Vec<_>>() {
            let action = actions.action(StereoBondId::from(relation_id))?;
            if action.degree() != set.participants_2(relation_id).len() {
                return None;
            }
            *set.data_mut(relation_id) = set.data(relation_id).clone().reframe_by(action)?;
            set.permute_2_with(relation_id, &participant_order(*action));
        }
        Some(self)
    }
}

impl Reframe for StereoBonds {
    fn representative_action(&self) -> Self::Action {
        let actions = self
            .ids()
            .map(|id| {
                stereo_bond_representative_action(self.ligands(id))
                    .expect("integrity-valid stereo-bond frames have four ligands")
            })
            .collect();
        StereoBondsFrameAction::from_vec(actions)
            .expect("every selected stereo-bond action is block-preserving")
    }

    fn reframe(self) -> Result<Self, Contradiction> {
        reframe_stereo_bonds_with(self, |_, _| {})
    }
}

pub(crate) fn reframe_stereo_bonds_with(
    mut stereo_bonds: StereoBonds,
    mut visit: impl FnMut(StereoBondId, Permutation),
) -> Result<StereoBonds, Contradiction> {
    let set = Arc::make_mut(&mut stereo_bonds.0);
    for relation_id in set.ids().collect::<Vec<_>>() {
        let id = StereoBondId::from(relation_id);
        let action = stereo_bond_representative_action(set.participants_2(relation_id))
            .ok_or(Contradiction)?;
        let attributes = set.data(relation_id).clone().normalize()?;
        *set.data_mut(relation_id) = attributes
            .reframe_by(&action)
            .ok_or(Contradiction)?
            .normalize()?;
        set.permute_2_with(relation_id, &participant_order(action));
        visit(id, action);
    }
    Ok(stereo_bonds)
}

/// The reaction span's stereo atoms, one [`EntitySpan`] per entity against a single ligand frame.
/// The `Molecule` peer is [`StereoAtoms`]. Values are issued by
/// [`ReactionSpan`](super::reaction_span::ReactionSpan).
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct StereoAtomSpans(
    FixedVarBirelationSet<NodeId, 1, StereoLigand, EntitySpan<StereoAtomForm>>,
);

impl StereoAtomSpans {
    pub(crate) fn into_entries(
        self,
    ) -> Vec<(AtomId, Vec<StereoLigand>, EntitySpan<StereoAtomForm>)> {
        self.0
            .into_entries()
            .into_iter()
            .map(|(site, ligands, span)| (AtomId::from(site[0]), ligands, span))
            .collect()
    }

    pub(crate) fn new(
        entries: Vec<(AtomId, Vec<StereoLigand>, EntitySpan<StereoAtomForm>)>,
    ) -> Self {
        Self(FixedVarBirelationSet::new(
            entries
                .into_iter()
                .map(|(site, ligands, span)| ([NodeId::from(site)], ligands, span))
                .collect(),
        ))
    }

    pub fn count(&self) -> usize {
        self.0.count()
    }

    pub fn contains(&self, id: StereoAtomId) -> bool {
        self.0.contains(RelationId::from(id))
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = StereoAtomId> {
        self.0.ids().map(StereoAtomId::from)
    }

    /// The atom `id` is borne by.
    pub fn site(&self, id: StereoAtomId) -> AtomId {
        AtomId::from(self.0.participants_1(RelationId::from(id))[0])
    }

    /// The ligands of `id`, in the frame both carried sides are read against.
    pub fn ligands(&self, id: StereoAtomId) -> &[StereoLigand] {
        self.0.participants_2(RelationId::from(id))
    }

    pub fn attributes(&self, id: StereoAtomId) -> &EntitySpan<StereoAtomForm> {
        self.0.data(RelationId::from(id))
    }

    /// Map participant references, preserving entity ids, row order, attributes, and frames.
    ///
    /// # Panics
    /// Panics if a referenced node or edge has no image in `correspondence`.
    pub fn map(&self, correspondence: &GraphCorrespondence) -> Self {
        self.try_map(correspondence)
            .expect("correspondence must cover every participant reference")
    }

    /// Map participant references, or return `None` if any reference has no image.
    /// Unreferenced nodes and edges need not have images. No entity is dropped.
    pub fn try_map(&self, correspondence: &GraphCorrespondence) -> Option<Self> {
        self.0.try_map(correspondence).map(Self)
    }

    pub(crate) fn remap(&self, remapping: &GraphRemapping) -> Self {
        Self(self.0.remap(remapping))
    }
}

impl Normalize for StereoAtomSpans {
    fn normalize(mut self) -> Result<Self, Contradiction> {
        for id in self.0.ids().collect::<Vec<_>>() {
            *self.0.data_mut(id) = self.0.data(id).clone().normalize()?;
        }
        Ok(self)
    }
}

impl FrameTransport for StereoAtomSpans {
    type Action = StereoAtomsFrameAction;

    fn reframe_by(mut self, actions: &Self::Action) -> Option<Self> {
        for relation_id in self.0.ids().collect::<Vec<_>>() {
            let action = actions.action(StereoAtomId::from(relation_id))?;
            if action.degree() != self.0.participants_2(relation_id).len() {
                return None;
            }
            *self.0.data_mut(relation_id) = self.0.data(relation_id).clone().reframe_by(action)?;
            self.0
                .permute_2_with(relation_id, &participant_order(*action));
        }
        Some(self)
    }
}

impl Reframe for StereoAtomSpans {
    fn representative_action(&self) -> Self::Action {
        let actions = self
            .ids()
            .map(|id| {
                stereo_atom_representative_action(self.ligands(id))
                    .expect("integrity-valid stereo-atom frames fit the bounded action")
            })
            .collect();
        StereoAtomsFrameAction::from_vec(actions)
            .expect("every bounded permutation is a stereo-atom action")
    }

    fn reframe(self) -> Result<Self, Contradiction> {
        reframe_stereo_atom_spans_with(self, |_, _| {})
    }
}

pub(crate) fn reframe_stereo_atom_spans_with(
    mut stereo_atoms: StereoAtomSpans,
    mut visit: impl FnMut(StereoAtomId, Permutation),
) -> Result<StereoAtomSpans, Contradiction> {
    for relation_id in stereo_atoms.0.ids().collect::<Vec<_>>() {
        let id = StereoAtomId::from(relation_id);
        let action = stereo_atom_representative_action(stereo_atoms.0.participants_2(relation_id))
            .ok_or(Contradiction)?;
        let span = stereo_atoms.0.data(relation_id).clone().normalize()?;
        *stereo_atoms.0.data_mut(relation_id) =
            span.reframe_by(&action).ok_or(Contradiction)?.normalize()?;
        stereo_atoms
            .0
            .permute_2_with(relation_id, &participant_order(action));
        visit(id, action);
    }
    Ok(stereo_atoms)
}

/// The reaction span's stereo bonds, one [`EntitySpan`] per entity against a single ligand frame.
/// The `Molecule` peer is [`StereoBonds`]. Values are issued by
/// [`ReactionSpan`](super::reaction_span::ReactionSpan).
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct StereoBondSpans(
    FixedVarBirelationSet<EdgeId, 1, StereoLigand, EntitySpan<StereoBondForm>>,
);

impl StereoBondSpans {
    pub(crate) fn into_entries(
        self,
    ) -> Vec<(BondId, Vec<StereoLigand>, EntitySpan<StereoBondForm>)> {
        self.0
            .into_entries()
            .into_iter()
            .map(|(site, ligands, span)| (BondId::from(site[0]), ligands, span))
            .collect()
    }

    pub(crate) fn new(
        entries: Vec<(BondId, Vec<StereoLigand>, EntitySpan<StereoBondForm>)>,
    ) -> Self {
        Self(FixedVarBirelationSet::new(
            entries
                .into_iter()
                .map(|(site, ligands, span)| ([EdgeId::from(site)], ligands, span))
                .collect(),
        ))
    }

    pub fn count(&self) -> usize {
        self.0.count()
    }

    pub fn contains(&self, id: StereoBondId) -> bool {
        self.0.contains(RelationId::from(id))
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = StereoBondId> {
        self.0.ids().map(StereoBondId::from)
    }

    /// The bond `id` is borne by.
    pub fn site(&self, id: StereoBondId) -> BondId {
        BondId::from(self.0.participants_1(RelationId::from(id))[0])
    }

    /// The ligands of `id`, in the frame both carried sides are read against.
    pub fn ligands(&self, id: StereoBondId) -> &[StereoLigand] {
        self.0.participants_2(RelationId::from(id))
    }

    pub fn attributes(&self, id: StereoBondId) -> &EntitySpan<StereoBondForm> {
        self.0.data(RelationId::from(id))
    }

    /// Map participant references, preserving entity ids, row order, attributes, and frames.
    ///
    /// # Panics
    /// Panics if a referenced node or edge has no image in `correspondence`.
    pub fn map(&self, correspondence: &GraphCorrespondence) -> Self {
        self.try_map(correspondence)
            .expect("correspondence must cover every participant reference")
    }

    /// Map participant references, or return `None` if any reference has no image.
    /// Unreferenced nodes and edges need not have images. No entity is dropped.
    pub fn try_map(&self, correspondence: &GraphCorrespondence) -> Option<Self> {
        self.0.try_map(correspondence).map(Self)
    }

    pub(crate) fn remap(&self, remapping: &GraphRemapping) -> Self {
        Self(self.0.remap(remapping))
    }
}

impl Normalize for StereoBondSpans {
    fn normalize(mut self) -> Result<Self, Contradiction> {
        for id in self.0.ids().collect::<Vec<_>>() {
            *self.0.data_mut(id) = self.0.data(id).clone().normalize()?;
        }
        Ok(self)
    }
}

impl FrameTransport for StereoBondSpans {
    type Action = StereoBondsFrameAction;

    fn reframe_by(mut self, actions: &Self::Action) -> Option<Self> {
        for relation_id in self.0.ids().collect::<Vec<_>>() {
            let action = actions.action(StereoBondId::from(relation_id))?;
            if action.degree() != self.0.participants_2(relation_id).len() {
                return None;
            }
            *self.0.data_mut(relation_id) = self.0.data(relation_id).clone().reframe_by(action)?;
            self.0
                .permute_2_with(relation_id, &participant_order(*action));
        }
        Some(self)
    }
}

impl Reframe for StereoBondSpans {
    fn representative_action(&self) -> Self::Action {
        let actions = self
            .ids()
            .map(|id| {
                stereo_bond_representative_action(self.ligands(id))
                    .expect("integrity-valid stereo-bond frames have four ligands")
            })
            .collect();
        StereoBondsFrameAction::from_vec(actions)
            .expect("every selected stereo-bond action is block-preserving")
    }

    fn reframe(self) -> Result<Self, Contradiction> {
        reframe_stereo_bond_spans_with(self, |_, _| {})
    }
}

pub(crate) fn reframe_stereo_bond_spans_with(
    mut stereo_bonds: StereoBondSpans,
    mut visit: impl FnMut(StereoBondId, Permutation),
) -> Result<StereoBondSpans, Contradiction> {
    for relation_id in stereo_bonds.0.ids().collect::<Vec<_>>() {
        let id = StereoBondId::from(relation_id);
        let action = stereo_bond_representative_action(stereo_bonds.0.participants_2(relation_id))
            .ok_or(Contradiction)?;
        let span = stereo_bonds.0.data(relation_id).clone().normalize()?;
        *stereo_bonds.0.data_mut(relation_id) =
            span.reframe_by(&action).ok_or(Contradiction)?.normalize()?;
        stereo_bonds
            .0
            .permute_2_with(relation_id, &participant_order(action));
        visit(id, action);
    }
    Ok(stereo_bonds)
}

/// Defines the stereo entity forms.
macro_rules! stereo_element {
    (
        $(#[doc = $doc:literal])+
        $name:ident, $constraints:ident, $constraint:ident, $allows_action:expr
    ) => {
        $(#[doc = $doc])+
        #[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Normalize, Lattice)]
        pub struct $name {
            pub configuration: StereoConfigurationForm,
            pub constraints: $constraints,
        }

        impl $name {
            /// Concrete: every inherent field is ground; the constraint
            /// channel does not bear on concreteness.
            pub fn is_concrete(&self) -> bool {
                let $name {
                    configuration,
                    constraints: _,
                } = self;
                configuration.is_ground()
            }
        }


        impl FrameTransport for $constraint {
            type Action = Permutation;

            fn reframe_by(self, action: &Self::Action) -> Option<Self> {
                if !($allows_action)(*action) {
                    return None;
                }
                let inverse = action.inverse();
                let reframe_permutation = |value: LigandPermutation| {
                    (value.0.degree() == action.degree()).then(|| {
                        LigandPermutation(inverse.compose(value.0).compose(*action))
                    })
                };
                Some(match self {
                    Self::LigandSymmetry(symmetry) => Self::LigandSymmetry(LigandSymmetryForm {
                        permutation: OrientedLigandPermutation {
                            permutation: reframe_permutation(symmetry.permutation.permutation)?,
                            orientation: symmetry.permutation.orientation,
                        },
                        invariant: symmetry.invariant,
                    }),
                    Self::Fluxionality(fluxionality) => Self::Fluxionality(FluxionalityForm {
                        permutation: reframe_permutation(fluxionality.permutation)?,
                        active: fluxionality.active,
                    }),
                    Self::Topicity(topicity) => {
                        let first = topicity.pair.first().index();
                        let second = topicity.pair.second().index();
                        if first >= action.degree() || second >= action.degree() {
                            return None;
                        }
                        Self::Topicity(TopicityForm {
                            pair: StereoLigandPair::new(
                                inverse.apply(first).into(),
                                inverse.apply(second).into(),
                            ),
                            relation: topicity.relation,
                        })
                    }
                    Self::Stereogenicity(stereogenicity) => {
                        Self::Stereogenicity(stereogenicity)
                    }
                })
            }
        }

        impl FrameTransport for $constraints {
            type Action = Permutation;

            fn reframe_by(self, action: &Self::Action) -> Option<Self> {
                if !($allows_action)(*action) {
                    return None;
                }
                self.into_iter()
                    .map(|constraint| {
                        constraint.reframe_by(action)
                    })
                    .collect()
            }
        }

        /// Configuration, permutation-valued constraints, and topicity positions move together.
        /// `None` when the permutation is not an action of the configured stereo kind or its degree
        /// is incompatible with a frame-relative constraint.
        impl FrameTransport for $name {
            type Action = Permutation;

            fn reframe_by(self, action: &Self::Action) -> Option<Self> {
                if !($allows_action)(*action) {
                    return None;
                }
                let Self {
                    configuration,
                    constraints,
                } = self;
                Some(Self {
                    configuration: configuration.reframe_by(action)?,
                    constraints: constraints.reframe_by(action)?,
                })
            }
        }

        impl From<&str> for $name {
            fn from(s: &str) -> Self {
                s.parse()
                    .expect(concat!("invalid ", stringify!($name), " string"))
            }
        }

        impl $name {
            pub fn new(kind: StereoKind, coset: impl Into<StereoCoset>) -> Self {
                Self {
                    configuration: StereoConfigurationForm::kinded(kind, coset),
                    constraints: $constraints::new(),
                }
            }

            /// Add a single constraint, replacing any existing entry with the same key.
            pub fn with_constraint(mut self, constraint: impl Into<$constraint>) -> Self {
                self.constraints.set(constraint.into());
                self
            }

            /// Add each constraint from iterator.
            pub fn with_constraints<I>(mut self, constraints: I) -> Self
            where
                I: IntoIterator,
                I::Item: Into<$constraint>,
            {
                self.constraints.extend(constraints.into_iter().map(Into::into));
                self
            }

            /// No-op. A stereo element is always stereogenic, so its coset has no
            /// zero default; it is ground iff its coset is ground.
            pub fn into_concrete(self) -> Self {
                self
            }

            /// Apply a ligand-position permutation to the configuration in its current frame.
            /// Frame-relative constraints are unchanged; use [`FrameTransport::reframe_by`] when the
            /// stored ligand frame itself is reordered.
            pub fn apply(&self, permutation: Permutation) -> Option<Self> {
                Some(Self {
                    configuration: self.configuration.apply(permutation)?,
                    constraints: self.constraints.clone(),
                })
            }

            /// The kind involution (`~`).
            pub fn swap(&self) -> Option<Self> {
                Some(Self {
                    configuration: self.configuration.swap()?,
                    constraints: self.constraints.clone(),
                })
            }

            /// The enantiomer / mirror (`'`).
            pub fn mirror(&self) -> Option<Self> {
                Some(Self {
                    configuration: self.configuration.mirror()?,
                    constraints: self.constraints.clone(),
                })
            }

        }

    };
}

stereo_element! {
    /// Stereo atom form with geometry class, configuration, and per-site constraints.
    StereoAtomForm, StereoAtomConstraintsForm, StereoAtomConstraintForm, |_: Permutation| true
}

stereo_element! {
    /// Stereo bond form with cis/trans configuration and per-site constraints.
    StereoBondForm, StereoBondConstraintsForm, StereoBondConstraintForm, |action: Permutation| {
        StereoKind::CisTrans.class_key().space().allows(action)
    }
}

/// Configuration portion of a stereo-element update.
///
/// `Unchanged` omits the field, `Undetermined` explicitly clears it, and
/// `Kinded` carries either an absolute coset or a kind-only relative update.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoConfigurationUpdate {
    #[default]
    Unchanged,
    Undetermined,
    Kinded {
        kind: StereoKind,
        coset: Option<StereoCoset>,
    },
}

impl StereoConfigurationUpdate {
    fn apply_to(&self, current: &StereoConfigurationForm) -> StereoConfigurationForm {
        match self {
            Self::Unchanged => current.clone(),
            Self::Undetermined => StereoConfigurationForm::Undetermined,
            Self::Kinded {
                kind,
                coset: Some(coset),
            } => StereoConfigurationForm::kinded(*kind, coset.clone()),
            Self::Kinded { kind, coset: None } => match current {
                StereoConfigurationForm::Kinded(current_kind, current_coset)
                    if current_kind == kind =>
                {
                    StereoConfigurationForm::kinded(*kind, current_coset.clone())
                }
                _ => StereoConfigurationForm::kinded(*kind, StereoCoset::Undetermined),
            },
        }
    }

    pub(crate) fn kind(&self) -> Option<StereoKind> {
        match self {
            Self::Kinded { kind, .. } => Some(*kind),
            Self::Unchanged | Self::Undetermined => None,
        }
    }
}

/// Attribute update for a stereo atom.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StereoAtomUpdate {
    pub configuration: StereoConfigurationUpdate,
    pub constraints: StereoAtomConstraintsForm,
}

/// Attribute update for a stereo bond.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StereoBondUpdate {
    pub configuration: StereoConfigurationUpdate,
    pub constraints: StereoBondConstraintsForm,
}

impl StereoAtomForm {
    /// Apply an attribute update.
    pub fn update(&self, update: &StereoAtomUpdate) -> Self {
        let mut constraints = self.constraints.clone();
        constraints.update(&update.constraints);
        Self {
            configuration: update.configuration.apply_to(&self.configuration),
            constraints,
        }
    }

    /// Derive the minimal normalized attribute update carrying `self` to `other`.
    pub fn difference_to(&self, other: &Self) -> StereoAtomUpdate {
        let mut constraints = StereoAtomConstraintsForm::new();
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
        let configuration = if self.configuration.normalized_eq(&other.configuration) {
            match self.configuration.kind() {
                Some(kind) if !constraints.is_empty() => {
                    StereoConfigurationUpdate::Kinded { kind, coset: None }
                }
                _ => StereoConfigurationUpdate::Unchanged,
            }
        } else {
            match &other.configuration {
                StereoConfigurationForm::Undetermined => StereoConfigurationUpdate::Undetermined,
                StereoConfigurationForm::Kinded(kind, coset) => StereoConfigurationUpdate::Kinded {
                    kind: *kind,
                    coset: Some(coset.clone()),
                },
            }
        };
        StereoAtomUpdate {
            configuration,
            constraints,
        }
    }
}

impl StereoBondForm {
    /// Apply an attribute update.
    pub fn update(&self, update: &StereoBondUpdate) -> Self {
        let mut constraints = self.constraints.clone();
        constraints.update(&update.constraints);
        Self {
            configuration: update.configuration.apply_to(&self.configuration),
            constraints,
        }
    }

    /// Derive the minimal normalized attribute update carrying `self` to `other`.
    pub fn difference_to(&self, other: &Self) -> StereoBondUpdate {
        let mut constraints = StereoBondConstraintsForm::new();
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
        let configuration = if self.configuration.normalized_eq(&other.configuration) {
            match self.configuration.kind() {
                Some(kind) if !constraints.is_empty() => {
                    StereoConfigurationUpdate::Kinded { kind, coset: None }
                }
                _ => StereoConfigurationUpdate::Unchanged,
            }
        } else {
            match &other.configuration {
                StereoConfigurationForm::Undetermined => StereoConfigurationUpdate::Undetermined,
                StereoConfigurationForm::Kinded(kind, coset) => StereoConfigurationUpdate::Kinded {
                    kind: *kind,
                    coset: Some(coset.clone()),
                },
            }
        };
        StereoBondUpdate {
            configuration,
            constraints,
        }
    }
}

/// Stereo kind: the atom-centered coordination geometries and the bond-centered cis/trans kind.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, strum::EnumCount)]
pub enum StereoKind {
    Tetrahedral,
    CisTrans,
    Axial,
    SquarePlanar,
    TrigonalBipyramidal,
    Octahedral,
}

impl StereoKind {
    /// The `umol-perm` class key for this stereo kind.
    pub fn class_key(self) -> ClassKey {
        match self {
            StereoKind::Tetrahedral => ClassKey::Tetrahedral,
            StereoKind::CisTrans => ClassKey::CisTrans,
            StereoKind::Axial => ClassKey::Axial,
            StereoKind::SquarePlanar => ClassKey::SquarePlanar,
            StereoKind::TrigonalBipyramidal => ClassKey::TrigonalBipyramidal,
            StereoKind::Octahedral => ClassKey::Octahedral,
        }
    }

    /// Number of ligand positions in this stereo kind.
    pub fn degree(self) -> usize {
        self.class_key().space().degree()
    }

    /// Number of cosets/configurations in this stereo kind.
    pub fn count(self) -> usize {
        self.class_key().space().count()
    }

    /// Whether this stereo kind can encode local handedness.
    pub fn is_chiral_class(self) -> bool {
        self.class_key().space().is_chiral()
    }

    /// Kind-specific `~` involution. Chiral kinds borrow the orientation-reversing
    /// generator from umol-perm; achiral kinds use a chosen ligand swap (no improper
    /// generator to borrow — theirs is the identity):
    /// - cis/trans: swap the two configurations
    /// - square-planar: swap the diagonal ligand pair
    pub fn involution(self) -> Permutation {
        let coset_space = self.class_key().space();
        if coset_space.is_chiral() {
            coset_space.improper()
        } else {
            match self {
                StereoKind::CisTrans => Permutation::from_image(&[1, 0, 2, 3]),
                StereoKind::SquarePlanar => Permutation::from_image(&[2, 1, 0, 3]),
                _ => unreachable!("only achiral kinds reach the chosen-swap branch"),
            }
        }
    }

    /// Act on coset index `index` by `permutation`, through the class's coset algebra.
    pub fn act(self, index: u32, permutation: Permutation) -> Option<u32> {
        self.class_key().space().reindex(index, permutation)
    }

    /// The mirror (improper, μ) generator as a permutation: chiral kinds use the
    /// orientation-reversing generator; achiral kinds act trivially on cosets.
    pub fn mirror_permutation(self) -> Permutation {
        if self.is_chiral_class() {
            self.class_key().space().improper()
        } else {
            Permutation::identity(self.degree())
        }
    }

    /// Whether `g` and `h` induce the same coset permutation for this kind.
    fn coset_action_eq(self, g: Permutation, h: Permutation) -> bool {
        let s = self.class_key().space();
        (0..s.count() as u32).all(|i| s.reindex(i, g) == s.reindex(i, h))
    }

    /// Normalize coset permutation, priority `Mirror > Swap > Apply`; `None`
    /// when it acts as the identity on cosets.
    pub fn canonicalize_permutation(self, g: Permutation) -> Option<CosetOp> {
        if self.coset_action_eq(g, Permutation::identity(self.degree())) {
            None
        } else if self.is_chiral_class() && self.coset_action_eq(g, self.mirror_permutation()) {
            Some(CosetOp::Mirror)
        } else if self.coset_action_eq(g, self.involution()) {
            Some(CosetOp::Swap)
        } else {
            Some(CosetOp::Apply(g))
        }
    }
}

/// Permutation in canonical priority form, `Mirror` > `Swap` > `Apply`, kind-dependent.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CosetOp {
    Swap,
    Mirror,
    Apply(Permutation),
}

/// Topicity of two ligand positions of a stereo carrier (derived ground value).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, VariantArray)]
pub enum Topicity {
    Homotopic,
    Enantiotopic,
    Diastereotopic,
}

/// Stereogenicity classification of a stereo carrier (derived ground value).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, VariantArray)]
pub enum Stereogenicity {
    Symmetric,
    Prochiral,
    Stereogenic,
}

/// Element-side stereo configuration: either undetermined (geometry not yet
/// known, so no coset) or `Kinded` — a concrete geometry bound to a coset that
/// may still be open. `*` (`Undetermined`) and `Th*` (`Kinded(Tetrahedral,
/// Undetermined)`) are distinct. `normalize` folds the coset under the kind;
/// no physical range-check (tier-2; the validator does it).
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoConfigurationForm {
    #[default]
    Undetermined,
    Kinded(StereoKind, StereoCoset),
}

impl StereoConfigurationForm {
    pub fn undetermined() -> Self {
        Self::Undetermined
    }

    pub fn kinded(kind: StereoKind, coset: impl Into<StereoCoset>) -> Self {
        Self::Kinded(kind, coset.into())
    }

    /// The coordination-geometry kind, or `None` when undetermined.
    pub fn kind(&self) -> Option<StereoKind> {
        match self {
            Self::Kinded(kind, _) => Some(*kind),
            Self::Undetermined => None,
        }
    }

    /// The coset, or `None` when undetermined.
    pub fn coset(&self) -> Option<&StereoCoset> {
        match self {
            Self::Kinded(_, coset) => Some(coset),
            Self::Undetermined => None,
        }
    }

    /// Mutable access to the coset, or `None` when undetermined.
    pub fn coset_mut(&mut self) -> Option<&mut StereoCoset> {
        match self {
            Self::Kinded(_, coset) => Some(coset),
            Self::Undetermined => None,
        }
    }

    /// Relabel the ligand positions (`^`); `Undetermined` is fixed.
    pub fn apply(&self, permutation: Permutation) -> Option<Self> {
        self.map_kinded(|kind, coset| coset.apply(kind, permutation))
    }

    /// The kind involution (`~`).
    pub fn swap(&self) -> Option<Self> {
        self.map_kinded(|kind, coset| coset.swap(kind))
    }

    /// The enantiomer / mirror (`'`).
    pub fn mirror(&self) -> Option<Self> {
        self.map_kinded(|kind, coset| coset.mirror(kind))
    }

    fn map_kinded(
        &self,
        f: impl FnOnce(StereoKind, &StereoCoset) -> Option<StereoCoset>,
    ) -> Option<Self> {
        Some(match self {
            Self::Undetermined => Self::Undetermined,
            Self::Kinded(kind, coset) => Self::Kinded(*kind, f(*kind, coset)?),
        })
    }

    /// Overwrite with `other`, field-wise: an `Undetermined` `other` keeps `self`; a same-kind
    /// `other` with an `Undetermined` coset keeps `self`'s coset (a partial that fixes only the
    /// kind); a differing kind or a determined coset overrides wholesale.
    pub fn update(&self, other: &Self) -> Self {
        match (self, other) {
            (_, Self::Undetermined) => self.clone(),
            (Self::Kinded(ks, cs), Self::Kinded(ko, StereoCoset::Undetermined)) if ks == ko => {
                Self::Kinded(*ks, cs.clone())
            }
            (_, Self::Kinded(..)) => other.clone(),
        }
    }
}

impl From<(StereoKind, u32)> for StereoConfigurationForm {
    fn from((kind, coset): (StereoKind, u32)) -> Self {
        Self::Kinded(kind, StereoCoset::Lit(coset))
    }
}

impl Normalize for StereoConfigurationForm {
    fn normalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::Kinded(kind, coset) => Self::Kinded(kind, canon_coset(coset, kind)?),
            Self::Undetermined => Self::Undetermined,
        })
    }

    fn normalized(&self) -> Result<Cow<'_, Self>, Contradiction> {
        match self {
            Self::Kinded(..) => Ok(Cow::Owned(self.clone().normalize()?)),
            Self::Undetermined => Ok(Cow::Borrowed(self)),
        }
    }
}

impl AsLit for StereoConfigurationForm {
    type Lit = StereoConfiguration;

    fn as_lit(&self) -> Option<StereoConfiguration> {
        match self {
            Self::Kinded(kind, coset) => coset
                .as_lit()
                .map(|coset| StereoConfiguration { kind: *kind, coset }),
            Self::Undetermined => None,
        }
    }
}

impl FrameTransport for StereoConfigurationForm {
    type Action = Permutation;

    /// The coset under the kind's action; an undetermined configuration has no kind and carries
    /// unchanged. `None` when the permutation is not an action of the configured kind.
    fn reframe_by(self, action: &Self::Action) -> Option<Self> {
        if self.kind().is_some_and(|kind| {
            kind.degree() != action.degree() || !kind.class_key().space().allows(*action)
        }) {
            return None;
        }
        self.apply(*action)
    }
}

impl Lattice for StereoConfigurationForm {
    fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }

    fn is_ground(&self) -> bool {
        matches!(self, Self::Kinded(_, StereoCoset::Lit(_)))
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        let a = self.normalized().ok()?;
        let b = other.normalized().ok()?;
        match (a.as_ref(), b.as_ref()) {
            (Self::Undetermined, x) | (x, Self::Undetermined) => Some(x.clone()),
            (Self::Kinded(k1, ca), Self::Kinded(k2, cb)) => {
                if k1 != k2 {
                    return None;
                }
                Some(Self::Kinded(*k1, coset_meet(ca, cb, *k1)?))
            }
        }
    }

    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        let a = self.normalized().unwrap_or(Cow::Owned(Self::Undetermined));
        let b = other.normalized().unwrap_or(Cow::Owned(Self::Undetermined));
        Ok(match (a.as_ref(), b.as_ref()) {
            (Self::Undetermined, _) | (_, Self::Undetermined) => Self::Undetermined,
            (Self::Kinded(k1, ca), Self::Kinded(k2, cb)) => {
                if k1 != k2 {
                    Self::Undetermined
                } else {
                    Self::Kinded(*k1, coset_join(ca, cb, *k1))
                }
            }
        })
    }
}

/// Generates a constraint-side stereo state for a fixed geometry (`#T`/`#C`):
/// undetermined, explicitly not-stereo, or a stereo center with a coset. The
/// geometry is the type's identity (`$kind`), so the coset folds/meets under that
/// constant kind — no kind field.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum TetrahedralStereo {
    NotStereo,
    Stereo(u32),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CisTransStereo {
    NotStereo,
    Stereo(u32),
}

/// Named tetrahedral configurations and their canonical coset indices.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum TetrahedralConfiguration {
    Ccw,
    Cw,
}

/// Named cis/trans configurations and their canonical coset indices.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CisTransConfiguration {
    Z,
    E,
}

macro_rules! stereo_site {
    ($name:ident, $lit:ident, $kind:expr) => {
        #[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub enum $name {
            #[default]
            Undetermined,
            NotStereo,
            Stereo(StereoCoset),
        }

        impl $name {
            pub fn undetermined() -> Self {
                Self::Undetermined
            }

            pub fn not_stereo() -> Self {
                Self::NotStereo
            }

            pub fn stereo(coset: impl Into<StereoCoset>) -> Self {
                Self::Stereo(coset.into())
            }

            pub fn is_stereo(&self) -> bool {
                matches!(self, Self::Stereo(_))
            }

            /// Matches literal coset index `value` under the type's kind.
            pub fn matches_value(&self, value: u32) -> bool {
                match self {
                    Self::Stereo(coset) => coset_matches(coset, &StereoCoset::Lit(value), $kind),
                    Self::NotStereo => false,
                    Self::Undetermined => true,
                }
            }
        }

        impl Normalize for $name {
            fn normalize(self) -> Result<Self, Contradiction> {
                Ok(match self {
                    Self::Stereo(coset) => Self::Stereo(canon_coset(coset, $kind)?),
                    other => other,
                })
            }

            fn normalized(&self) -> Result<Cow<'_, Self>, Contradiction> {
                match self {
                    Self::Stereo(_) => Ok(Cow::Owned(self.clone().normalize()?)),
                    _ => Ok(Cow::Borrowed(self)),
                }
            }
        }

        impl AsLit for $name {
            type Lit = $lit;

            /// The exact absence or stereo-coset value when ground.
            fn as_lit(&self) -> Option<$lit> {
                match self {
                    Self::NotStereo => Some($lit::NotStereo),
                    Self::Stereo(StereoCoset::Lit(coset)) => Some($lit::Stereo(*coset)),
                    _ => None,
                }
            }
        }

        impl Lattice for $name {
            fn is_undetermined(&self) -> bool {
                matches!(self, Self::Undetermined)
            }

            fn is_ground(&self) -> bool {
                match self {
                    Self::NotStereo => true,
                    Self::Stereo(coset) => matches!(coset, StereoCoset::Lit(_)),
                    Self::Undetermined => false,
                }
            }

            fn meet(&self, other: &Self) -> Option<Self> {
                let a = self.normalized().ok()?;
                let b = other.normalized().ok()?;
                match (a.as_ref(), b.as_ref()) {
                    (Self::Undetermined, x) | (x, Self::Undetermined) => Some(x.clone()),
                    (Self::NotStereo, Self::NotStereo) => Some(Self::NotStereo),
                    (Self::NotStereo, Self::Stereo(_)) | (Self::Stereo(_), Self::NotStereo) => None,
                    (Self::Stereo(ca), Self::Stereo(cb)) => {
                        Some(Self::Stereo(coset_meet(ca, cb, $kind)?))
                    }
                }
            }

            fn join(&self, other: &Self) -> Result<Self, NoJoin> {
                let a = self.normalized().unwrap_or(Cow::Owned(Self::Undetermined));
                let b = other.normalized().unwrap_or(Cow::Owned(Self::Undetermined));
                Ok(match (a.as_ref(), b.as_ref()) {
                    (Self::Undetermined, _) | (_, Self::Undetermined) => Self::Undetermined,
                    (Self::NotStereo, Self::NotStereo) => Self::NotStereo,
                    (Self::NotStereo, Self::Stereo(_)) | (Self::Stereo(_), Self::NotStereo) => {
                        Self::Undetermined
                    }
                    (Self::Stereo(ca), Self::Stereo(cb)) => Self::Stereo(coset_join(ca, cb, $kind)),
                })
            }
        }
    };
}

stereo_site! { TetrahedralStereoForm, TetrahedralStereo, StereoKind::Tetrahedral }
stereo_site! { CisTransStereoForm, CisTransStereo, StereoKind::CisTrans }

impl From<TetrahedralStereo> for TetrahedralStereoForm {
    fn from(stereo: TetrahedralStereo) -> Self {
        match stereo {
            TetrahedralStereo::NotStereo => Self::NotStereo,
            TetrahedralStereo::Stereo(coset) => Self::Stereo(StereoCoset::Lit(coset)),
        }
    }
}

impl From<CisTransStereo> for CisTransStereoForm {
    fn from(stereo: CisTransStereo) -> Self {
        match stereo {
            CisTransStereo::NotStereo => Self::NotStereo,
            CisTransStereo::Stereo(coset) => Self::Stereo(StereoCoset::Lit(coset)),
        }
    }
}

impl From<TetrahedralConfiguration> for TetrahedralStereo {
    fn from(configuration: TetrahedralConfiguration) -> Self {
        match configuration {
            TetrahedralConfiguration::Ccw => Self::Stereo(0),
            TetrahedralConfiguration::Cw => Self::Stereo(1),
        }
    }
}

impl From<CisTransConfiguration> for CisTransStereo {
    fn from(configuration: CisTransConfiguration) -> Self {
        match configuration {
            CisTransConfiguration::Z => Self::Stereo(0),
            CisTransConfiguration::E => Self::Stereo(1),
        }
    }
}

impl From<TetrahedralConfiguration> for TetrahedralStereoForm {
    fn from(configuration: TetrahedralConfiguration) -> Self {
        TetrahedralStereo::from(configuration).into()
    }
}

impl From<CisTransConfiguration> for CisTransStereoForm {
    fn from(configuration: CisTransConfiguration) -> Self {
        CisTransStereo::from(configuration).into()
    }
}

/// Operator-expression term: a `Var` (with optional finite domain), a literal
/// `Lit`/`LitSet` base, or one of these under the permutation-action operators
/// `~` (swap), `'` (mirror), `^` (apply). Kind-relative — **no
/// `Lattice`/`Normalize`** (structural `Eq` only); the owning configuration
/// normalizes it under its concrete kind. Normalization composes the operator
/// word into one net permutation: over a literal base it folds to a concrete
/// coset; over a `Var` it leaves at most one operator layer.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoTerm {
    Var(Box<(String, Option<BTreeSet<u32>>)>),
    Lit(u32),
    LitSet(BTreeSet<u32>),
    Swap(Box<StereoTerm>),
    Mirror(Box<StereoTerm>),
    Apply(Box<StereoTerm>, Permutation),
}

impl StereoTerm {
    /// A free coset variable `?name`.
    pub fn var(name: impl Into<String>) -> Self {
        Self::Var(Box::new((name.into(), None)))
    }

    /// A variable restricted to a finite coset domain `?name :: {…}`.
    pub fn var_in(name: impl Into<String>, domain: impl IntoIterator<Item = u32>) -> Self {
        Self::Var(Box::new((name.into(), Some(domain.into_iter().collect()))))
    }

    /// A finite literal-set base `{…}` (folds under the owner's kind).
    pub fn lit_set(values: impl IntoIterator<Item = u32>) -> Self {
        Self::LitSet(values.into_iter().collect())
    }

    /// `~inner` — the kind involution applied to `inner`.
    pub fn swap(inner: Self) -> Self {
        Self::Swap(Box::new(inner))
    }

    /// `'inner` — the enantiomer (mirror) of `inner`.
    pub fn mirror(inner: Self) -> Self {
        Self::Mirror(Box::new(inner))
    }

    /// `inner ^ permutation` — the group action of `permutation` on `inner`.
    pub fn apply(inner: Self, permutation: Permutation) -> Self {
        Self::Apply(Box::new(inner), permutation)
    }
}

/// Dense coset-index expression (0-indexed per stereo kind): undetermined, a
/// single index, a finite set, or an operator `Term` over a variable.
/// Kind-relative — no `Lattice` or `Normalize`; the owning configuration or
/// site normalizes it under its kind.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoCoset {
    #[default]
    Undetermined,
    Lit(u32),
    LitSet(BTreeSet<u32>),
    Term(Box<StereoTerm>),
}

impl StereoCoset {
    pub fn lit_set(values: impl IntoIterator<Item = u32>) -> Self {
        Self::LitSet(values.into_iter().collect())
    }

    pub fn term(term: StereoTerm) -> Self {
        Self::Term(Box::new(term))
    }

    /// Relabel the ligand positions (the `^` op): move each literal coset index through the kind's
    /// coset algebra, eager on `Lit`/`LitSet`; an open `Term` keeps the operator layer.
    fn apply(&self, kind: StereoKind, permutation: Permutation) -> Option<Self> {
        self.map_index(
            |c| kind.act(c, permutation),
            |t| StereoTerm::apply(t, permutation),
        )
    }

    /// The kind involution (the `~` op).
    fn swap(&self, kind: StereoKind) -> Option<Self> {
        self.map_index(|c| kind.act(c, kind.involution()), StereoTerm::swap)
    }

    /// The enantiomer / mirror (the `'` op).
    fn mirror(&self, kind: StereoKind) -> Option<Self> {
        self.map_index(
            |c| kind.act(c, kind.mirror_permutation()),
            StereoTerm::mirror,
        )
    }

    /// Map each literal index by `lit`; an open `Term` is wrapped by `term` (the only case that keeps
    /// an operator layer — a bare variable cannot be evaluated). `Undetermined` is fixed.
    fn map_index(
        &self,
        lit: impl Fn(u32) -> Option<u32>,
        term: impl FnOnce(StereoTerm) -> StereoTerm,
    ) -> Option<Self> {
        Some(match self {
            Self::Undetermined => Self::Undetermined,
            Self::Lit(c) => Self::Lit(lit(*c)?),
            Self::LitSet(s) => Self::LitSet(s.iter().map(|&c| lit(c)).collect::<Option<_>>()?),
            Self::Term(t) => Self::term(term((**t).clone())),
        })
    }
}

impl From<u32> for StereoCoset {
    fn from(index: u32) -> Self {
        Self::Lit(index)
    }
}

impl AsLit for StereoCoset {
    type Lit = u32;

    /// The single coset index, only when literal. Kind-independent — so it lives on
    /// the bare coset, unlike the kind-aware folding ops.
    #[inline]
    fn as_lit(&self) -> Option<u32> {
        match self {
            Self::Lit(i) => Some(*i),
            _ => None,
        }
    }
}

/// A ground stereo configuration: a concrete geometry plus its coset index. The
/// `AsLit` target of `StereoConfigurationForm` and the per-kind site types.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct StereoConfiguration {
    pub kind: StereoKind,
    pub coset: u32,
}

/// The literal coset-index set a positive coset denotes; `None` for the wildcard
/// `Undetermined` and the symbolic `Term`. Used by `coset_meet`/`coset_join`
/// after those two cases are handled.
fn coset_to_set(coset: &StereoCoset) -> Option<BTreeSet<u32>> {
    match coset {
        StereoCoset::Lit(i) => Some(BTreeSet::from([*i])),
        StereoCoset::LitSet(s) => Some(s.clone()),
        StereoCoset::Undetermined | StereoCoset::Term(_) => None,
    }
}

/// Walk a term's operator word into one net coset permutation (composed inner →
/// outer), returning the base leaf (`Var`/`Lit`/`LitSet`) and that permutation.
fn compose_term(term: &StereoTerm, kind: StereoKind) -> (&StereoTerm, Permutation) {
    match term {
        StereoTerm::Swap(inner) => {
            let (base, g) = compose_term(inner, kind);
            (base, g.compose(kind.involution()))
        }
        StereoTerm::Mirror(inner) => {
            let (base, g) = compose_term(inner, kind);
            (base, g.compose(kind.mirror_permutation()))
        }
        StereoTerm::Apply(inner, p) => {
            let (base, g) = compose_term(inner, kind);
            (base, g.compose(*p))
        }
        base => (base, Permutation::identity(kind.degree())),
    }
}

/// Normalize a coset under `kind`. A `Term` over a `Var` renders by priority
/// Mirror > Swap > Apply (canonicalizing the domain); every other form reduces to
/// a literal index set that folds: ∅ → `Err` (the bottom `meet` uses to signal
/// incompatible cosets), singleton → `Lit`, else `LitSet`. No universe folding
/// (`full → Undetermined`) and no range-check — both are tier-2 (the validator).
pub(crate) fn canon_coset(
    coset: StereoCoset,
    kind: StereoKind,
) -> Result<StereoCoset, Contradiction> {
    let s = kind.class_key().space();
    let set: BTreeSet<u32> = match &coset {
        StereoCoset::Undetermined => return Ok(StereoCoset::Undetermined),
        StereoCoset::Lit(i) => BTreeSet::from([*i]),
        StereoCoset::LitSet(values) => values.clone(),
        StereoCoset::Term(t) => {
            let (base, g) = compose_term(t, kind);
            match base {
                StereoTerm::Var(v) => {
                    let n = kind.count() as u32;
                    let domain = match &v.1 {
                        Some(set) if set.is_empty() => return Err(Contradiction),
                        Some(set) if set.len() as u32 == n => None,
                        Some(set) => Some(set.clone()),
                        None => None,
                    };
                    let var = StereoTerm::Var(Box::new((v.0.clone(), domain)));
                    let term = match kind.canonicalize_permutation(g) {
                        None => var,
                        Some(CosetOp::Mirror) => StereoTerm::Mirror(Box::new(var)),
                        Some(CosetOp::Swap) => StereoTerm::Swap(Box::new(var)),
                        Some(CosetOp::Apply(g)) => StereoTerm::Apply(Box::new(var), g),
                    };
                    return Ok(StereoCoset::term(term));
                }
                StereoTerm::Lit(i) => BTreeSet::from([s.reindex(*i, g).ok_or(Contradiction)?]),
                StereoTerm::LitSet(values) => values
                    .iter()
                    .map(|i| s.reindex(*i, g).ok_or(Contradiction))
                    .collect::<Result<_, _>>()?,
                StereoTerm::Swap(_) | StereoTerm::Mirror(_) | StereoTerm::Apply(..) => {
                    unreachable!("compose_term returns a base leaf")
                }
            }
        }
    };
    if set.is_empty() {
        Err(Contradiction)
    } else if set.len() == 1 {
        Ok(StereoCoset::Lit(set.into_iter().next().unwrap()))
    } else {
        Ok(StereoCoset::LitSet(set))
    }
}

/// Greatest lower bound of two cosets under `kind` (canonicalizing operands);
/// `Term` meets only an equal normalized `Term`.
pub(crate) fn coset_meet(
    a: &StereoCoset,
    b: &StereoCoset,
    kind: StereoKind,
) -> Option<StereoCoset> {
    let ca = canon_coset(a.clone(), kind).ok()?;
    let cb = canon_coset(b.clone(), kind).ok()?;
    use StereoCoset::{Term, Undetermined};
    match (&ca, &cb) {
        (Undetermined, _) => Some(cb),
        (_, Undetermined) => Some(ca),
        (Term(_), Term(_)) => (ca == cb).then_some(ca),
        (Term(_), _) | (_, Term(_)) => None,
        _ => {
            let sa = coset_to_set(&ca).unwrap();
            let sb = coset_to_set(&cb).unwrap();
            canon_coset(
                StereoCoset::LitSet(sa.intersection(&sb).copied().collect()),
                kind,
            )
            .ok()
        }
    }
}

/// Least upper bound of two cosets under `kind`.
pub(crate) fn coset_join(a: &StereoCoset, b: &StereoCoset, kind: StereoKind) -> StereoCoset {
    let ca = canon_coset(a.clone(), kind).unwrap_or(StereoCoset::Undetermined);
    let cb = canon_coset(b.clone(), kind).unwrap_or(StereoCoset::Undetermined);
    use StereoCoset::{Term, Undetermined};
    match (&ca, &cb) {
        (Undetermined, _) | (_, Undetermined) => StereoCoset::Undetermined,
        (Term(_), Term(_)) if ca == cb => ca,
        (Term(_), _) | (_, Term(_)) => StereoCoset::Undetermined,
        _ => {
            let sa = coset_to_set(&ca).unwrap();
            let sb = coset_to_set(&cb).unwrap();
            canon_coset(StereoCoset::LitSet(sa.union(&sb).copied().collect()), kind)
                .unwrap_or(StereoCoset::Undetermined)
        }
    }
}

/// `target` refines `pattern` under `kind` (meet-derived).
pub(crate) fn coset_matches(pattern: &StereoCoset, target: &StereoCoset, kind: StereoKind) -> bool {
    match (
        coset_meet(pattern, target, kind),
        canon_coset(target.clone(), kind),
    ) {
        (Some(m), Ok(ct)) => m == ct,
        _ => false,
    }
}

/// Apply a ligand-order permutation to a coset under `kind`.
pub(crate) fn coset_apply_permutation(
    coset: &StereoCoset,
    permutation: Permutation,
    kind: StereoKind,
) -> Option<StereoCoset> {
    let s = kind.class_key().space();
    match coset {
        StereoCoset::Undetermined => Some(StereoCoset::Undetermined),
        StereoCoset::Lit(i) => Some(StereoCoset::Lit(s.reindex(*i, permutation)?)),
        StereoCoset::LitSet(set) => Some(StereoCoset::LitSet(
            set.iter()
                .map(|i| s.reindex(*i, permutation))
                .collect::<Option<_>>()?,
        )),
        StereoCoset::Term(t) => Some(
            canon_coset(
                StereoCoset::term(StereoTerm::apply((**t).clone(), permutation)),
                kind,
            )
            .unwrap_or(StereoCoset::Undetermined),
        ),
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::Correspondence;
    use umol_perm::Orientation;

    use super::super::boolean::BooleanForm;
    use super::super::constraint::{StereogenicityForm, TopicityRelationForm};
    use super::super::id::{AtomId, StereoLigandPosition};
    use super::super::ligand::StereoLigandKind;
    use super::*;

    #[rstest]
    #[case::covered(None, None)]
    #[case::missing_atom_0(Some(NodeId(0)), None)]
    #[case::missing_atom_2(Some(NodeId(2)), None)]
    #[case::missing_atom_4(Some(NodeId(4)), None)]
    #[case::missing_atom_6(Some(NodeId(6)), None)]

    fn test_stereo_atoms_try_map(
        #[case] missing_node: Option<NodeId>,
        #[case] missing_edge: Option<EdgeId>,
    ) {
        let input = StereoAtoms::new(vec![
            (
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
            ),
            (
                AtomId(6),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(6), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(6), StereoLigandKind::LonePair),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, 2_u32),
            ),
        ]);
        let correspondence = GraphCorrespondence::new(
            Correspondence::new(
                vec![
                    (NodeId(0), NodeId(5)),
                    (NodeId(2), NodeId(1)),
                    (NodeId(4), NodeId(7)),
                    (NodeId(6), NodeId(3)),
                ]
                .into_iter()
                .filter(|(id, _)| Some(*id) != missing_node)
                .collect(),
                8,
                9,
            )
            .unwrap(),
            Correspondence::new(
                vec![(EdgeId(0), EdgeId(4)), (EdgeId(2), EdgeId(1))]
                    .into_iter()
                    .filter(|(id, _)| Some(*id) != missing_edge)
                    .collect(),
                4,
                6,
            )
            .unwrap(),
        );
        let expected = if missing_node.is_none() && missing_edge.is_none() {
            Some(StereoAtoms::new(vec![
                (
                    AtomId(5),
                    vec![
                        StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(5), StereoLigandKind::LonePair),
                    ],
                    StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
                ),
                (
                    AtomId(3),
                    vec![
                        StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(3), StereoLigandKind::LonePair),
                    ],
                    StereoAtomForm::new(StereoKind::Tetrahedral, 2_u32),
                ),
            ]))
        } else {
            None
        };
        assert_eq!(input.try_map(&correspondence), expected);
        if let Some(expected) = expected {
            assert_eq!(input.map(&correspondence), expected);
        }
    }

    #[rstest]
    #[should_panic(expected = "correspondence must cover every participant reference")]
    fn test_stereo_atoms_map_error() {
        let input = StereoAtoms::new(vec![
            (
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
            ),
            (
                AtomId(6),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(6), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(6), StereoLigandKind::LonePair),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, 2_u32),
            ),
        ]);
        input.map(&GraphCorrespondence::new(
            Correspondence::empty(),
            Correspondence::empty(),
        ));
    }

    #[rstest]
    #[case::covered(None, None)]
    #[case::missing_atom_0(Some(NodeId(0)), None)]
    #[case::missing_atom_2(Some(NodeId(2)), None)]
    #[case::missing_atom_4(Some(NodeId(4)), None)]
    #[case::missing_atom_6(Some(NodeId(6)), None)]

    fn test_stereo_atom_spans_try_map(
        #[case] missing_node: Option<NodeId>,
        #[case] missing_edge: Option<EdgeId>,
    ) {
        let input = StereoAtomSpans::new(vec![
            (
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                ],
                EntitySpan::Modified {
                    lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
                    rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 2_u32),
                },
            ),
            (
                AtomId(6),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(6), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(6), StereoLigandKind::LonePair),
                ],
                EntitySpan::Added(StereoAtomForm::new(StereoKind::Tetrahedral, 2_u32)),
            ),
        ]);
        let correspondence = GraphCorrespondence::new(
            Correspondence::new(
                vec![
                    (NodeId(0), NodeId(5)),
                    (NodeId(2), NodeId(1)),
                    (NodeId(4), NodeId(7)),
                    (NodeId(6), NodeId(3)),
                ]
                .into_iter()
                .filter(|(id, _)| Some(*id) != missing_node)
                .collect(),
                8,
                9,
            )
            .unwrap(),
            Correspondence::new(
                vec![(EdgeId(0), EdgeId(4)), (EdgeId(2), EdgeId(1))]
                    .into_iter()
                    .filter(|(id, _)| Some(*id) != missing_edge)
                    .collect(),
                4,
                6,
            )
            .unwrap(),
        );
        let expected = if missing_node.is_none() && missing_edge.is_none() {
            Some(StereoAtomSpans::new(vec![
                (
                    AtomId(5),
                    vec![
                        StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(5), StereoLigandKind::LonePair),
                    ],
                    EntitySpan::Modified {
                        lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
                        rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 2_u32),
                    },
                ),
                (
                    AtomId(3),
                    vec![
                        StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(3), StereoLigandKind::LonePair),
                    ],
                    EntitySpan::Added(StereoAtomForm::new(StereoKind::Tetrahedral, 2_u32)),
                ),
            ]))
        } else {
            None
        };
        assert_eq!(input.try_map(&correspondence), expected);
        if let Some(expected) = expected {
            assert_eq!(input.map(&correspondence), expected);
        }
    }

    #[rstest]
    #[should_panic(expected = "correspondence must cover every participant reference")]
    fn test_stereo_atom_spans_map_error() {
        let input = StereoAtomSpans::new(vec![
            (
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                ],
                EntitySpan::Modified {
                    lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
                    rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 2_u32),
                },
            ),
            (
                AtomId(6),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(6), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(6), StereoLigandKind::LonePair),
                ],
                EntitySpan::Added(StereoAtomForm::new(StereoKind::Tetrahedral, 2_u32)),
            ),
        ]);
        input.map(&GraphCorrespondence::new(
            Correspondence::empty(),
            Correspondence::empty(),
        ));
    }

    #[rstest]
    #[case::covered(None, None)]
    #[case::missing_atom_0(Some(NodeId(0)), None)]
    #[case::missing_atom_2(Some(NodeId(2)), None)]
    #[case::missing_atom_4(Some(NodeId(4)), None)]
    #[case::missing_atom_6(Some(NodeId(6)), None)]
    #[case::missing_bond_0(None, Some(EdgeId(0)))]
    #[case::missing_bond_2(None, Some(EdgeId(2)))]
    fn test_stereo_bonds_try_map(
        #[case] missing_node: Option<NodeId>,
        #[case] missing_edge: Option<EdgeId>,
    ) {
        let input = StereoBonds::new(vec![
            (
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(6), StereoLigandKind::LonePair),
                ],
                StereoBondForm::new(StereoKind::CisTrans, 1_u32),
            ),
            (
                BondId(2),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(6), StereoLigandKind::LonePair),
                ],
                StereoBondForm::new(StereoKind::CisTrans, 2_u32),
            ),
        ]);
        let correspondence = GraphCorrespondence::new(
            Correspondence::new(
                vec![
                    (NodeId(0), NodeId(5)),
                    (NodeId(2), NodeId(1)),
                    (NodeId(4), NodeId(7)),
                    (NodeId(6), NodeId(3)),
                ]
                .into_iter()
                .filter(|(id, _)| Some(*id) != missing_node)
                .collect(),
                8,
                9,
            )
            .unwrap(),
            Correspondence::new(
                vec![(EdgeId(0), EdgeId(4)), (EdgeId(2), EdgeId(1))]
                    .into_iter()
                    .filter(|(id, _)| Some(*id) != missing_edge)
                    .collect(),
                4,
                6,
            )
            .unwrap(),
        );
        let expected = if missing_node.is_none() && missing_edge.is_none() {
            Some(StereoBonds::new(vec![
                (
                    BondId(4),
                    vec![
                        StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(3), StereoLigandKind::LonePair),
                    ],
                    StereoBondForm::new(StereoKind::CisTrans, 1_u32),
                ),
                (
                    BondId(1),
                    vec![
                        StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(3), StereoLigandKind::LonePair),
                    ],
                    StereoBondForm::new(StereoKind::CisTrans, 2_u32),
                ),
            ]))
        } else {
            None
        };
        assert_eq!(input.try_map(&correspondence), expected);
        if let Some(expected) = expected {
            assert_eq!(input.map(&correspondence), expected);
        }
    }

    #[rstest]
    #[should_panic(expected = "correspondence must cover every participant reference")]
    fn test_stereo_bonds_map_error() {
        let input = StereoBonds::new(vec![
            (
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(6), StereoLigandKind::LonePair),
                ],
                StereoBondForm::new(StereoKind::CisTrans, 1_u32),
            ),
            (
                BondId(2),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(6), StereoLigandKind::LonePair),
                ],
                StereoBondForm::new(StereoKind::CisTrans, 2_u32),
            ),
        ]);
        input.map(&GraphCorrespondence::new(
            Correspondence::empty(),
            Correspondence::empty(),
        ));
    }

    #[rstest]
    #[case::covered(None, None)]
    #[case::missing_atom_0(Some(NodeId(0)), None)]
    #[case::missing_atom_2(Some(NodeId(2)), None)]
    #[case::missing_atom_4(Some(NodeId(4)), None)]
    #[case::missing_atom_6(Some(NodeId(6)), None)]
    #[case::missing_bond_0(None, Some(EdgeId(0)))]
    #[case::missing_bond_2(None, Some(EdgeId(2)))]
    fn test_stereo_bond_spans_try_map(
        #[case] missing_node: Option<NodeId>,
        #[case] missing_edge: Option<EdgeId>,
    ) {
        let input = StereoBondSpans::new(vec![
            (
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(6), StereoLigandKind::LonePair),
                ],
                EntitySpan::Modified {
                    lhs: StereoBondForm::new(StereoKind::CisTrans, 1_u32),
                    rhs: StereoBondForm::new(StereoKind::CisTrans, 2_u32),
                },
            ),
            (
                BondId(2),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(6), StereoLigandKind::LonePair),
                ],
                EntitySpan::Added(StereoBondForm::new(StereoKind::CisTrans, 2_u32)),
            ),
        ]);
        let correspondence = GraphCorrespondence::new(
            Correspondence::new(
                vec![
                    (NodeId(0), NodeId(5)),
                    (NodeId(2), NodeId(1)),
                    (NodeId(4), NodeId(7)),
                    (NodeId(6), NodeId(3)),
                ]
                .into_iter()
                .filter(|(id, _)| Some(*id) != missing_node)
                .collect(),
                8,
                9,
            )
            .unwrap(),
            Correspondence::new(
                vec![(EdgeId(0), EdgeId(4)), (EdgeId(2), EdgeId(1))]
                    .into_iter()
                    .filter(|(id, _)| Some(*id) != missing_edge)
                    .collect(),
                4,
                6,
            )
            .unwrap(),
        );
        let expected = if missing_node.is_none() && missing_edge.is_none() {
            Some(StereoBondSpans::new(vec![
                (
                    BondId(4),
                    vec![
                        StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(3), StereoLigandKind::LonePair),
                    ],
                    EntitySpan::Modified {
                        lhs: StereoBondForm::new(StereoKind::CisTrans, 1_u32),
                        rhs: StereoBondForm::new(StereoKind::CisTrans, 2_u32),
                    },
                ),
                (
                    BondId(1),
                    vec![
                        StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(3), StereoLigandKind::LonePair),
                    ],
                    EntitySpan::Added(StereoBondForm::new(StereoKind::CisTrans, 2_u32)),
                ),
            ]))
        } else {
            None
        };
        assert_eq!(input.try_map(&correspondence), expected);
        if let Some(expected) = expected {
            assert_eq!(input.map(&correspondence), expected);
        }
    }

    #[rstest]
    #[should_panic(expected = "correspondence must cover every participant reference")]
    fn test_stereo_bond_spans_map_error() {
        let input = StereoBondSpans::new(vec![
            (
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(6), StereoLigandKind::LonePair),
                ],
                EntitySpan::Modified {
                    lhs: StereoBondForm::new(StereoKind::CisTrans, 1_u32),
                    rhs: StereoBondForm::new(StereoKind::CisTrans, 2_u32),
                },
            ),
            (
                BondId(2),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(6), StereoLigandKind::LonePair),
                ],
                EntitySpan::Added(StereoBondForm::new(StereoKind::CisTrans, 2_u32)),
            ),
        ]);
        input.map(&GraphCorrespondence::new(
            Correspondence::empty(),
            Correspondence::empty(),
        ));
    }

    #[rstest]
    #[case::tetrahedral(StereoKind::Tetrahedral, ClassKey::Tetrahedral)]
    #[case::cis_trans(StereoKind::CisTrans, ClassKey::CisTrans)]
    #[case::axial(StereoKind::Axial, ClassKey::Axial)]
    #[case::square_planar(StereoKind::SquarePlanar, ClassKey::SquarePlanar)]
    #[case::trigonal_bipyramidal(StereoKind::TrigonalBipyramidal, ClassKey::TrigonalBipyramidal)]
    #[case::octahedral(StereoKind::Octahedral, ClassKey::Octahedral)]
    fn test_stereo_kind_class_key(#[case] kind: StereoKind, #[case] expected: ClassKey) {
        assert_eq!(kind.class_key(), expected);
    }

    #[rstest]
    #[case::tetrahedral(StereoKind::Tetrahedral, 4)]
    #[case::cis_trans(StereoKind::CisTrans, 4)]
    #[case::axial(StereoKind::Axial, 4)]
    #[case::square_planar(StereoKind::SquarePlanar, 4)]
    #[case::trigonal_bipyramidal(StereoKind::TrigonalBipyramidal, 5)]
    #[case::octahedral(StereoKind::Octahedral, 6)]
    fn test_stereo_kind_degree(#[case] kind: StereoKind, #[case] expected: usize) {
        assert_eq!(kind.degree(), expected);
    }

    #[rstest]
    #[case::tetrahedral(StereoKind::Tetrahedral, 2)]
    #[case::cis_trans(StereoKind::CisTrans, 2)]
    #[case::axial(StereoKind::Axial, 2)]
    #[case::square_planar(StereoKind::SquarePlanar, 3)]
    #[case::trigonal_bipyramidal(StereoKind::TrigonalBipyramidal, 20)]
    #[case::octahedral(StereoKind::Octahedral, 30)]
    fn test_stereo_kind_count(#[case] kind: StereoKind, #[case] expected: usize) {
        assert_eq!(kind.count(), expected);
    }

    #[rstest]
    #[case::tetrahedral(StereoKind::Tetrahedral, true)]
    #[case::cis_trans(StereoKind::CisTrans, false)]
    #[case::axial(StereoKind::Axial, true)]
    #[case::square_planar(StereoKind::SquarePlanar, false)]
    #[case::trigonal_bipyramidal(StereoKind::TrigonalBipyramidal, true)]
    #[case::octahedral(StereoKind::Octahedral, true)]
    fn test_stereo_kind_is_chiral_class(#[case] kind: StereoKind, #[case] expected: bool) {
        assert_eq!(kind.is_chiral_class(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::tetrahedral((StereoKind::Tetrahedral, 1), StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)))]
    #[case::octahedral((StereoKind::Octahedral, 5), StereoConfigurationForm::Kinded(StereoKind::Octahedral, StereoCoset::Lit(5)))]
    fn test_stereo_configuration_form_from(#[case] input: (StereoKind, u32), #[case] expected: StereoConfigurationForm) {
        assert_eq!(StereoConfigurationForm::from(input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::term_swap_folds_to_lit(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(0)))), StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)))]
    fn test_stereo_configuration_form_normalize(#[case] input: StereoConfigurationForm, #[case] expected: StereoConfigurationForm) {
        assert_eq!(input.normalize(), Ok(expected));
    }

    #[rstest]
    #[case::undetermined(StereoConfigurationForm::Undetermined)]
    #[case::kind_lit(StereoConfigurationForm::Kinded(
        StereoKind::Tetrahedral,
        StereoCoset::Lit(0)
    ))]
    #[case::kind_open(StereoConfigurationForm::Kinded(
        StereoKind::Tetrahedral,
        StereoCoset::Undetermined
    ))]
    // Multi-element / full coset sets are preserved (no complement or full→Undetermined fold).
    #[case::multi_element_set(StereoConfigurationForm::Kinded(StereoKind::SquarePlanar, StereoCoset::lit_set([0, 1])))]
    #[case::full_set(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::lit_set([0, 1])))]
    fn test_stereo_configuration_form_normalize_identity(#[case] input: StereoConfigurationForm) {
        assert_eq!(input.clone().normalize(), Ok(input));
    }

    #[rstest]
    #[case::empty_set(StereoConfigurationForm::Kinded(
        StereoKind::SquarePlanar,
        StereoCoset::LitSet(BTreeSet::new())
    ))]
    fn test_stereo_configuration_form_normalize_error(#[case] input: StereoConfigurationForm) {
        assert_eq!(input.normalize(), Err(Contradiction));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::kind_lit(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)), Some(StereoConfiguration { kind: StereoKind::Tetrahedral, coset: 1 }))]
    #[case::undetermined(StereoConfigurationForm::Undetermined, None)]
    #[case::kind_open(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined), None)]
    fn test_stereo_configuration_form_as_lit(#[case] config: StereoConfigurationForm, #[case] expected: Option<StereoConfiguration>) {
        assert_eq!(config.as_lit(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(StereoConfigurationForm::Undetermined, true)]
    #[case::kind_open(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined), false)]
    #[case::kind_lit(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)), false)]
    fn test_stereo_configuration_form_is_undetermined(#[case] config: StereoConfigurationForm, #[case] expected: bool) {
        assert_eq!(config.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::kind_lit(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)), true)]
    #[case::undetermined(StereoConfigurationForm::Undetermined, false)]
    #[case::kind_open(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined), false)]
    fn test_stereo_configuration_form_is_ground(#[case] config: StereoConfigurationForm, #[case] expected: bool) {
        assert_eq!(config.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_narrows(StereoConfigurationForm::Undetermined, StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), Some(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0))))]
    #[case::coset_same(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), Some(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0))))]
    #[case::open_narrows(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined), StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), Some(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0))))]
    #[case::coset_conflict(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), StereoConfigurationForm::from((StereoKind::Tetrahedral, 1)), None)]
    #[case::kind_conflict(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), StereoConfigurationForm::from((StereoKind::CisTrans, 0)), None)]
    fn test_stereo_configuration_form_meet(#[case] a: StereoConfigurationForm, #[case] b: StereoConfigurationForm, #[case] expected: Option<StereoConfigurationForm>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_absorbs(StereoConfigurationForm::Undetermined, StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), StereoConfigurationForm::Undetermined)]
    #[case::coset_same(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)))]
    #[case::coset_widens(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), StereoConfigurationForm::from((StereoKind::Tetrahedral, 1)), StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::lit_set([0, 1])))]
    #[case::kind_conflict(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), StereoConfigurationForm::from((StereoKind::CisTrans, 0)), StereoConfigurationForm::Undetermined)]
    fn test_stereo_configuration_form_join(#[case] a: StereoConfigurationForm, #[case] b: StereoConfigurationForm, #[case] expected: StereoConfigurationForm) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_matches_any(StereoConfigurationForm::Undetermined, StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), true)]
    #[case::open_matches_lit(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined), StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), true)]
    #[case::specific_vs_undetermined(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), StereoConfigurationForm::Undetermined, false)]
    #[case::coset_match(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), true)]
    #[case::coset_mismatch(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), StereoConfigurationForm::from((StereoKind::Tetrahedral, 1)), false)]
    #[case::kind_mismatch(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), StereoConfigurationForm::from((StereoKind::CisTrans, 0)), false)]
    fn test_stereo_configuration_form_matches(#[case] pattern: StereoConfigurationForm, #[case] target: StereoConfigurationForm, #[case] expected: bool) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    #[case::not_stereo(TetrahedralStereo::NotStereo, TetrahedralStereoForm::NotStereo)]
    #[case::stereo(
        TetrahedralStereo::Stereo(1),
        TetrahedralStereoForm::Stereo(StereoCoset::Lit(1))
    )]
    fn test_tetrahedral_stereo_form_from(
        #[case] stereo: TetrahedralStereo,
        #[case] expected: TetrahedralStereoForm,
    ) {
        assert_eq!(TetrahedralStereoForm::from(stereo), expected);
    }

    #[rstest]
    #[case::ccw(
        TetrahedralConfiguration::Ccw,
        TetrahedralStereoForm::Stereo(StereoCoset::Lit(0))
    )]
    #[case::cw(
        TetrahedralConfiguration::Cw,
        TetrahedralStereoForm::Stereo(StereoCoset::Lit(1))
    )]
    fn test_tetrahedral_stereo_form_from_configuration(
        #[case] configuration: TetrahedralConfiguration,
        #[case] expected: TetrahedralStereoForm,
    ) {
        assert_eq!(TetrahedralStereoForm::from(configuration), expected);
    }

    #[rstest]
    #[case::not_stereo(CisTransStereo::NotStereo, CisTransStereoForm::NotStereo)]
    #[case::stereo(
        CisTransStereo::Stereo(1),
        CisTransStereoForm::Stereo(StereoCoset::Lit(1))
    )]
    fn test_cis_trans_stereo_form_from(
        #[case] stereo: CisTransStereo,
        #[case] expected: CisTransStereoForm,
    ) {
        assert_eq!(CisTransStereoForm::from(stereo), expected);
    }

    #[rstest]
    #[case::z(
        CisTransConfiguration::Z,
        CisTransStereoForm::Stereo(StereoCoset::Lit(0))
    )]
    #[case::e(
        CisTransConfiguration::E,
        CisTransStereoForm::Stereo(StereoCoset::Lit(1))
    )]
    fn test_cis_trans_stereo_form_from_configuration(
        #[case] configuration: CisTransConfiguration,
        #[case] expected: CisTransStereoForm,
    ) {
        assert_eq!(CisTransStereoForm::from(configuration), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::term_swap_folds(TetrahedralStereoForm::Stereo(StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(0)))), TetrahedralStereoForm::Stereo(StereoCoset::Lit(1)))]
    fn test_tetrahedral_stereo_form_normalize(#[case] input: TetrahedralStereoForm, #[case] expected: TetrahedralStereoForm) {
        assert_eq!(input.normalize(), Ok(expected));
    }

    #[rstest]
    #[case::undetermined(TetrahedralStereoForm::Undetermined)]
    #[case::not_stereo(TetrahedralStereoForm::NotStereo)]
    #[case::stereo_lit(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)))]
    #[case::stereo_open(TetrahedralStereoForm::Stereo(StereoCoset::Undetermined))]
    fn test_tetrahedral_stereo_form_normalize_identity(#[case] input: TetrahedralStereoForm) {
        assert_eq!(input.clone().normalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stereo_lit(TetrahedralStereoForm::Stereo(StereoCoset::Lit(1)), Some(TetrahedralStereo::Stereo(1)))]
    #[case::stereo_zero(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), Some(TetrahedralStereo::Stereo(0)))]
    #[case::not_stereo(TetrahedralStereoForm::NotStereo, Some(TetrahedralStereo::NotStereo))]
    #[case::undetermined(TetrahedralStereoForm::Undetermined, None)]
    #[case::stereo_open(TetrahedralStereoForm::Stereo(StereoCoset::Undetermined), None)]
    fn test_tetrahedral_stereo_form_as_lit(#[case] site: TetrahedralStereoForm, #[case] expected: Option<TetrahedralStereo>) {
        assert_eq!(site.as_lit(), expected);
        assert_eq!(site.is_ground(), expected.is_some());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(TetrahedralStereoForm::Undetermined, true)]
    #[case::not_stereo(TetrahedralStereoForm::NotStereo, false)]
    #[case::stereo(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), false)]
    fn test_tetrahedral_stereo_form_is_undetermined(#[case] site: TetrahedralStereoForm, #[case] expected: bool) {
        assert_eq!(site.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::not_stereo(TetrahedralStereoForm::NotStereo, true)]
    #[case::stereo_lit(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), true)]
    #[case::undetermined(TetrahedralStereoForm::Undetermined, false)]
    #[case::stereo_open(TetrahedralStereoForm::Stereo(StereoCoset::Undetermined), false)]
    fn test_tetrahedral_stereo_form_is_ground(#[case] site: TetrahedralStereoForm, #[case] expected: bool) {
        assert_eq!(site.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(TetrahedralStereoForm::Undetermined, TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), Some(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0))))]
    #[case::not_stereo_same(TetrahedralStereoForm::NotStereo, TetrahedralStereoForm::NotStereo, Some(TetrahedralStereoForm::NotStereo))]
    #[case::not_stereo_vs_stereo(TetrahedralStereoForm::NotStereo, TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), None)]
    #[case::stereo_same(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), Some(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0))))]
    #[case::stereo_disjoint(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), TetrahedralStereoForm::Stereo(StereoCoset::Lit(1)), None)]
    #[case::open_narrows(TetrahedralStereoForm::Stereo(StereoCoset::Undetermined), TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), Some(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0))))]
    fn test_tetrahedral_stereo_form_meet(#[case] a: TetrahedralStereoForm, #[case] b: TetrahedralStereoForm, #[case] expected: Option<TetrahedralStereoForm>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(TetrahedralStereoForm::Undetermined, TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), TetrahedralStereoForm::Undetermined)]
    #[case::not_stereo_same(TetrahedralStereoForm::NotStereo, TetrahedralStereoForm::NotStereo, TetrahedralStereoForm::NotStereo)]
    #[case::not_stereo_vs_stereo(TetrahedralStereoForm::NotStereo, TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), TetrahedralStereoForm::Undetermined)]
    #[case::stereo_same(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)))]
    #[case::stereo_widens(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), TetrahedralStereoForm::Stereo(StereoCoset::Lit(1)), TetrahedralStereoForm::Stereo(StereoCoset::lit_set([0, 1])))]
    fn test_tetrahedral_stereo_form_join(#[case] a: TetrahedralStereoForm, #[case] b: TetrahedralStereoForm, #[case] expected: TetrahedralStereoForm) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(TetrahedralStereoForm::Undetermined, TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), true)]
    #[case::open_matches_lit(TetrahedralStereoForm::Stereo(StereoCoset::Undetermined), TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), true)]
    #[case::specific_vs_undetermined(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), TetrahedralStereoForm::Undetermined, false)]
    #[case::lit_match(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), true)]
    #[case::lit_mismatch(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), TetrahedralStereoForm::Stereo(StereoCoset::Lit(1)), false)]
    #[case::not_stereo_match(TetrahedralStereoForm::NotStereo, TetrahedralStereoForm::NotStereo, true)]
    #[case::not_stereo_vs_stereo(TetrahedralStereoForm::NotStereo, TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), false)]
    fn test_tetrahedral_stereo_form_matches(#[case] pattern: TetrahedralStereoForm, #[case] target: TetrahedralStereoForm, #[case] expected: bool) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(TetrahedralStereoForm::Undetermined, 0, true)]
    #[case::not_stereo(TetrahedralStereoForm::NotStereo, 0, false)]
    #[case::stereo_match(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), 0, true)]
    #[case::stereo_miss(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), 1, false)]
    fn test_tetrahedral_stereo_form_matches_value(#[case] site: TetrahedralStereoForm, #[case] value: u32, #[case] expected: bool) {
        assert_eq!(site.matches_value(value), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::term_swap_folds(CisTransStereoForm::Stereo(StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(0)))), CisTransStereoForm::Stereo(StereoCoset::Lit(1)))]
    fn test_cis_trans_stereo_form_normalize(#[case] input: CisTransStereoForm, #[case] expected: CisTransStereoForm) {
        assert_eq!(input.normalize(), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stereo_zero(CisTransStereoForm::Stereo(StereoCoset::Lit(0)), Some(CisTransStereo::Stereo(0)))]
    #[case::stereo_lit(CisTransStereoForm::Stereo(StereoCoset::Lit(1)), Some(CisTransStereo::Stereo(1)))]
    #[case::not_stereo(CisTransStereoForm::NotStereo, Some(CisTransStereo::NotStereo))]
    #[case::undetermined(CisTransStereoForm::Undetermined, None)]
    #[case::stereo_open(CisTransStereoForm::Stereo(StereoCoset::Undetermined), None)]
    fn test_cis_trans_stereo_form_as_lit(#[case] site: CisTransStereoForm, #[case] expected: Option<CisTransStereo>) {
        assert_eq!(site.as_lit(), expected);
        assert_eq!(site.is_ground(), expected.is_some());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(CisTransStereoForm::Undetermined, CisTransStereoForm::Stereo(StereoCoset::Lit(0)), Some(CisTransStereoForm::Stereo(StereoCoset::Lit(0))))]
    #[case::not_stereo_vs_stereo(CisTransStereoForm::NotStereo, CisTransStereoForm::Stereo(StereoCoset::Lit(0)), None)]
    #[case::stereo_disjoint(CisTransStereoForm::Stereo(StereoCoset::Lit(0)), CisTransStereoForm::Stereo(StereoCoset::Lit(1)), None)]
    fn test_cis_trans_stereo_form_meet(#[case] a: CisTransStereoForm, #[case] b: CisTransStereoForm, #[case] expected: Option<CisTransStereoForm>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(StereoCoset::Lit(2), Some(2))]
    #[case::undetermined(StereoCoset::Undetermined, None)]
    #[case::lit_set(StereoCoset::lit_set([1, 3]), None)]
    #[case::term(StereoCoset::term(StereoTerm::var("o")), None)]
    fn test_stereo_coset_as_lit(#[case] coset: StereoCoset, #[case] expected: Option<u32>) {
        assert_eq!(coset.as_lit(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_identity(StereoCoset::Lit(1), StereoKind::Tetrahedral, StereoCoset::Lit(1))]
    #[case::swap_lit_even(StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(0))), StereoKind::Tetrahedral, StereoCoset::Lit(1))]
    #[case::swap_lit_odd(StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(1))), StereoKind::Tetrahedral, StereoCoset::Lit(0))]
    #[case::mirror_chiral(StereoCoset::term(StereoTerm::mirror(StereoTerm::Lit(0))), StereoKind::Tetrahedral, StereoCoset::Lit(1))]
    #[case::mirror_achiral_noop(StereoCoset::term(StereoTerm::mirror(StereoTerm::Lit(0))), StereoKind::CisTrans, StereoCoset::Lit(0))]
    #[case::apply_lit(StereoCoset::term(StereoTerm::apply(StereoTerm::Lit(0), Permutation::from_image(&[1, 0, 2, 3]))), StereoKind::Tetrahedral, StereoCoset::Lit(1))]
    #[case::sp_swap_four(StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(1))), StereoKind::SquarePlanar, StereoCoset::Lit(2))]
    #[case::swap_var_chiral_to_mirror(StereoCoset::term(StereoTerm::swap(StereoTerm::var("o"))), StereoKind::Tetrahedral, StereoCoset::term(StereoTerm::mirror(StereoTerm::var("o"))))]
    #[case::swap_var_achiral_stays(StereoCoset::term(StereoTerm::swap(StereoTerm::var("o"))), StereoKind::CisTrans, StereoCoset::term(StereoTerm::swap(StereoTerm::var("o"))))]
    #[case::multi_element_set_preserved(StereoCoset::lit_set([0, 1]), StereoKind::SquarePlanar, StereoCoset::lit_set([0, 1]))]
    #[case::singleton_set_to_lit(StereoCoset::lit_set([1]), StereoKind::Octahedral, StereoCoset::Lit(1))]
    fn test_canon_coset(#[case] coset: StereoCoset, #[case] kind: StereoKind, #[case] expected: StereoCoset) {
        assert_eq!(canon_coset(coset, kind), Ok(expected));
    }

    #[rstest]
    #[case::empty_set(StereoCoset::LitSet(BTreeSet::new()), StereoKind::SquarePlanar)]
    fn test_canon_coset_error(#[case] coset: StereoCoset, #[case] kind: StereoKind) {
        assert_eq!(canon_coset(coset, kind), Err(Contradiction));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(StereoCoset::Undetermined, StereoCoset::Lit(1), StereoKind::Tetrahedral, Some(StereoCoset::Lit(1)))]
    #[case::lit_same(StereoCoset::Lit(1), StereoCoset::Lit(1), StereoKind::Tetrahedral, Some(StereoCoset::Lit(1)))]
    #[case::lit_disjoint(StereoCoset::Lit(0), StereoCoset::Lit(1), StereoKind::Tetrahedral, None)]
    #[case::set_intersect(StereoCoset::lit_set([1, 3]), StereoCoset::lit_set([3, 5]), StereoKind::Octahedral, Some(StereoCoset::Lit(3)))]
    #[case::term_equal(StereoCoset::term(StereoTerm::var("o")), StereoCoset::term(StereoTerm::var("o")), StereoKind::Tetrahedral, Some(StereoCoset::term(StereoTerm::var("o"))))]
    #[case::term_distinct(StereoCoset::term(StereoTerm::var("o")), StereoCoset::term(StereoTerm::var("p")), StereoKind::Tetrahedral, None)]
    fn test_coset_meet(#[case] a: StereoCoset, #[case] b: StereoCoset, #[case] kind: StereoKind, #[case] expected: Option<StereoCoset>) {
        assert_eq!(coset_meet(&a, &b, kind), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(StereoCoset::Undetermined, StereoCoset::Lit(1), StereoKind::Tetrahedral, StereoCoset::Undetermined)]
    #[case::lit_same(StereoCoset::Lit(1), StereoCoset::Lit(1), StereoKind::Tetrahedral, StereoCoset::Lit(1))]
    #[case::lit_union(StereoCoset::Lit(0), StereoCoset::Lit(1), StereoKind::Tetrahedral, StereoCoset::lit_set([0, 1]))]
    #[case::set_union(StereoCoset::lit_set([1, 3]), StereoCoset::lit_set([3, 5]), StereoKind::Octahedral, StereoCoset::lit_set([1, 3, 5]))]
    fn test_coset_join(#[case] a: StereoCoset, #[case] b: StereoCoset, #[case] kind: StereoKind, #[case] expected: StereoCoset) {
        assert_eq!(coset_join(&a, &b, kind), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(StereoCoset::Undetermined, StereoCoset::Lit(1), StereoKind::Tetrahedral, true)]
    #[case::lit_match(StereoCoset::Lit(1), StereoCoset::Lit(1), StereoKind::Tetrahedral, true)]
    #[case::lit_miss(StereoCoset::Lit(0), StereoCoset::Lit(1), StereoKind::Tetrahedral, false)]
    #[case::set_member(StereoCoset::lit_set([1, 3]), StereoCoset::Lit(3), StereoKind::Octahedral, true)]
    #[case::set_nonmember(StereoCoset::lit_set([1, 3]), StereoCoset::Lit(2), StereoKind::Octahedral, false)]
    #[case::specific_vs_wildcard(StereoCoset::Lit(0), StereoCoset::Undetermined, StereoKind::Tetrahedral, false)]
    fn test_coset_matches(#[case] pattern: StereoCoset, #[case] target: StereoCoset, #[case] kind: StereoKind, #[case] expected: bool) {
        assert_eq!(coset_matches(&pattern, &target, kind), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(StereoCoset::Lit(0), Permutation::from_image(&[1, 0, 2, 3]), StereoKind::Tetrahedral, StereoCoset::Lit(1))]
    #[case::undetermined(StereoCoset::Undetermined, Permutation::from_image(&[1, 0, 2, 3]), StereoKind::Tetrahedral, StereoCoset::Undetermined)]
    fn test_coset_apply_permutation(#[case] coset: StereoCoset, #[case] permutation: Permutation, #[case] kind: StereoKind, #[case] expected: StereoCoset) {
        assert_eq!(coset_apply_permutation(&coset, permutation, kind), Some(expected));
    }

    #[rstest]
    fn test_stereo_atom_form_new() {
        let stereo_atom = StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined);
        assert_eq!(
            stereo_atom.configuration,
            StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined)
        );
        assert_eq!(stereo_atom.constraints, StereoAtomConstraintsForm::new());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stereogenicity(
        StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
        StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        StereoAtomForm { configuration: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32), constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
    )]
    fn test_stereo_atom_form_with_methods(
        #[case] atom: StereoAtomForm,
        #[case] constraint: StereoAtomConstraintForm,
        #[case] expected: StereoAtomForm,
    ) {
        assert_eq!(atom.clone().with_constraint(constraint.clone()), expected);
        assert_eq!(atom.with_constraints([constraint]), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined), false)]
    #[case::ground(StereoAtomForm::new(StereoKind::Tetrahedral, 1u32), true)]
    fn test_stereo_atom_form_is_ground(#[case] atom: StereoAtomForm, #[case] expected: bool) {
        assert_eq!(atom.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::open_coset(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined))]
    #[case::ground(StereoAtomForm::new(StereoKind::Tetrahedral, 1u32))]
    fn test_stereo_atom_form_into_concrete(#[case] atom: StereoAtomForm) {
        assert_eq!(atom.clone().into_concrete(), atom);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::absolute(
        StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
        StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: Some(StereoCoset::Lit(1)) }, ..Default::default() },
        StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
    )]
    #[case::relative(
        StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
        StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: None }, ..Default::default() },
        StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
    )]
    #[case::undetermined(
        StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
        StereoAtomUpdate { configuration: StereoConfigurationUpdate::Undetermined, ..Default::default() },
        StereoAtomForm::default(),
    )]
    #[case::explicit_open(
        StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
        StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: Some(StereoCoset::Undetermined) }, ..Default::default() },
        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
    )]
    #[case::constraint_set(
        StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
        StereoAtomUpdate { constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))), ..Default::default() },
        StereoAtomForm { configuration: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32), constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
    )]
    #[case::constraint_remove(
        StereoAtomForm { configuration: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32), constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
        StereoAtomUpdate { constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)), ..Default::default() },
        StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
    )]
    fn test_stereo_atom_form_update(
        #[case] atom: StereoAtomForm,
        #[case] update: StereoAtomUpdate,
        #[case] expected: StereoAtomForm,
    ) {
        assert_eq!(atom.update(&update), expected);
    }

    #[rstest]
    #[case::empty(StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32))]
    fn test_stereo_atom_form_update_identity(#[case] atom: StereoAtomForm) {
        assert_eq!(atom.update(&StereoAtomUpdate::default()), atom);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::configuration_and_constraint(
        StereoAtomForm { configuration: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32), constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
        StereoAtomForm::default(),
        StereoAtomUpdate { configuration: StereoConfigurationUpdate::Undetermined, constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)) },
    )]
    #[case::constraint_context(
        StereoAtomForm { configuration: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32), constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
        StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
        StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: None }, constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)) },
    )]
    fn test_stereo_atom_form_difference_to(
        #[case] atom: StereoAtomForm,
        #[case] other: StereoAtomForm,
        #[case] expected: StereoAtomUpdate,
    ) {
        assert_eq!(atom.difference_to(&other), expected);
    }

    #[rstest]
    #[case::same(StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32))]
    fn test_stereo_atom_form_difference_to_identity(#[case] atom: StereoAtomForm) {
        assert_eq!(atom.difference_to(&atom), StereoAtomUpdate::default());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stereogenicity(
        StereoBondForm::new(StereoKind::CisTrans, 1_u32),
        StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        StereoBondForm { configuration: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1_u32), constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
    )]
    fn test_stereo_bond_form_with_methods(
        #[case] bond: StereoBondForm,
        #[case] constraint: StereoBondConstraintForm,
        #[case] expected: StereoBondForm,
    ) {
        assert_eq!(bond.clone().with_constraint(constraint.clone()), expected);
        assert_eq!(bond.with_constraints([constraint]), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::absolute(
        StereoBondForm::new(StereoKind::CisTrans, 0_u32),
        StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: Some(StereoCoset::Lit(1)) }, ..Default::default() },
        StereoBondForm::new(StereoKind::CisTrans, 1_u32),
    )]
    #[case::relative(
        StereoBondForm::new(StereoKind::CisTrans, 1_u32),
        StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: None }, ..Default::default() },
        StereoBondForm::new(StereoKind::CisTrans, 1_u32),
    )]
    #[case::undetermined(
        StereoBondForm::new(StereoKind::CisTrans, 1_u32),
        StereoBondUpdate { configuration: StereoConfigurationUpdate::Undetermined, ..Default::default() },
        StereoBondForm::default(),
    )]
    #[case::explicit_open(
        StereoBondForm::new(StereoKind::CisTrans, 1_u32),
        StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: Some(StereoCoset::Undetermined) }, ..Default::default() },
        StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined),
    )]
    #[case::constraint_set(
        StereoBondForm::new(StereoKind::CisTrans, 1_u32),
        StereoBondUpdate { constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))), ..Default::default() },
        StereoBondForm { configuration: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1_u32), constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
    )]
    #[case::constraint_remove(
        StereoBondForm { configuration: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1_u32), constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
        StereoBondUpdate { constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)), ..Default::default() },
        StereoBondForm::new(StereoKind::CisTrans, 1_u32),
    )]
    fn test_stereo_bond_form_update(
        #[case] bond: StereoBondForm,
        #[case] update: StereoBondUpdate,
        #[case] expected: StereoBondForm,
    ) {
        assert_eq!(bond.update(&update), expected);
    }

    #[rstest]
    #[case::empty(StereoBondForm::new(StereoKind::CisTrans, 1_u32))]
    fn test_stereo_bond_form_update_identity(#[case] bond: StereoBondForm) {
        assert_eq!(bond.update(&StereoBondUpdate::default()), bond);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::configuration_and_constraint(
        StereoBondForm { configuration: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1_u32), constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
        StereoBondForm::default(),
        StereoBondUpdate { configuration: StereoConfigurationUpdate::Undetermined, constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)) },
    )]
    #[case::constraint_context(
        StereoBondForm { configuration: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1_u32), constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
        StereoBondForm::new(StereoKind::CisTrans, 1_u32),
        StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: None }, constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)) },
    )]
    fn test_stereo_bond_form_difference_to(
        #[case] bond: StereoBondForm,
        #[case] other: StereoBondForm,
        #[case] expected: StereoBondUpdate,
    ) {
        assert_eq!(bond.difference_to(&other), expected);
    }

    #[rstest]
    #[case::same(StereoBondForm::new(StereoKind::CisTrans, 1_u32))]
    fn test_stereo_bond_form_difference_to_identity(#[case] bond: StereoBondForm) {
        assert_eq!(bond.difference_to(&bond), StereoBondUpdate::default());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_kind_narrows(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined), StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
        Some(StereoAtomForm::new(StereoKind::Tetrahedral, 1u32)))]
    #[case::different_kind(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
        StereoAtomForm::new(StereoKind::SquarePlanar, StereoCoset::Undetermined), None)]
    #[case::config_conflict(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32), StereoAtomForm::new(StereoKind::Tetrahedral, 1u32), None)]
    fn test_stereo_atom_form_meet(
        #[case] a: StereoAtomForm,
        #[case] b: StereoAtomForm,
        #[case] expected: Option<StereoAtomForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_coset(StereoAtomForm::new(StereoKind::Tetrahedral, 1u32), StereoAtomForm::new(StereoKind::Tetrahedral, 1u32), StereoAtomForm::new(StereoKind::Tetrahedral, 1u32))]
    #[case::distinct_cosets_widen(StereoAtomForm::new(StereoKind::Tetrahedral, 1u32), StereoAtomForm::new(StereoKind::Tetrahedral, 2u32), StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::lit_set([1, 2])))]
    fn test_stereo_atom_form_join(#[case] a: StereoAtomForm, #[case] b: StereoAtomForm, #[case] expected: StereoAtomForm) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_kind_match(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined), StereoAtomForm::new(StereoKind::Tetrahedral, 1u32), true)]
    #[case::different_kind(StereoAtomForm::new(StereoKind::Tetrahedral, 1u32), StereoAtomForm::new(StereoKind::SquarePlanar, 1u32), false)]
    fn test_stereo_atom_form_matches(
        #[case] pattern: StereoAtomForm,
        #[case] target: StereoAtomForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::folds_coset(
        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::lit_set([1])),
        Ok(StereoAtomForm::new(StereoKind::Tetrahedral, 1u32)),
    )]
    #[case::empty_coset_litset_contradiction(
        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::lit_set(Vec::<u32>::new())),
        Err(Contradiction),
    )]
    fn test_stereo_atom_form_normalize(
        #[case] input: StereoAtomForm,
        #[case] expected: Result<StereoAtomForm, Contradiction>,
    ) {
        assert_eq!(input.normalize(), expected);
    }

    #[rstest]
    fn test_stereo_bond_form_new() {
        let stereo_bond = StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined);
        assert_eq!(
            stereo_bond.configuration,
            StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Undetermined)
        );
        assert_eq!(stereo_bond.constraints, StereoBondConstraintsForm::new())
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_kind_narrows(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined), StereoBondForm::new(StereoKind::CisTrans, 1u32),
        Some(StereoBondForm::new(StereoKind::CisTrans, 1u32)))]
    #[case::config_conflict(StereoBondForm::new(StereoKind::CisTrans, 0u32), StereoBondForm::new(StereoKind::CisTrans, 1u32), None)]
    fn test_stereo_bond_form_meet(
        #[case] a: StereoBondForm,
        #[case] b: StereoBondForm,
        #[case] expected: Option<StereoBondForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::folds_coset(
        StereoBondForm::new(StereoKind::CisTrans, StereoCoset::lit_set([1])),
        Ok(StereoBondForm::new(StereoKind::CisTrans, 1u32)),
    )]
    #[case::empty_coset_litset_contradiction(
        StereoBondForm::new(StereoKind::CisTrans, StereoCoset::lit_set(Vec::<u32>::new())),
        Err(Contradiction),
    )]
    fn test_stereo_bond_form_normalize(
        #[case] input: StereoBondForm,
        #[case] expected: Result<StereoBondForm, Contradiction>,
    ) {
        assert_eq!(input.normalize(), expected);
    }

    #[rstest]
    #[case::identity(StereoKind::Tetrahedral, 0, Permutation::identity(4), 0)]
    #[case::involution(StereoKind::Tetrahedral, 0, StereoKind::Tetrahedral.involution(), 1)]
    #[case::involution_back(StereoKind::Tetrahedral, 1, StereoKind::Tetrahedral.involution(), 0)]
    fn test_stereo_kind_act(
        #[case] kind: StereoKind,
        #[case] index: u32,
        #[case] permutation: Permutation,
        #[case] expected: u32,
    ) {
        assert_eq!(kind.act(index, permutation), Some(expected));
    }

    #[rstest]
    #[case::coset_out_of_range(StereoKind::Tetrahedral, 2, Permutation::identity(4))]
    #[case::outside_parent_group(StereoKind::CisTrans, 0, Permutation::from_image(&[1, 2, 0, 3]))]
    #[case::wrong_degree(StereoKind::Tetrahedral, 0, Permutation::identity(3))]
    fn test_stereo_kind_act_error(
        #[case] kind: StereoKind,
        #[case] index: u32,
        #[case] permutation: Permutation,
    ) {
        assert_eq!(kind.act(index, permutation), None);
    }

    #[rstest]
    #[case::kinded(
        StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0_u32),
        Permutation::from_image(&[1, 0, 2, 3]),
        Some(StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32)),
    )]
    #[case::undetermined(
        StereoConfigurationForm::Undetermined,
        Permutation::identity(3),
        Some(StereoConfigurationForm::Undetermined)
    )]
    #[case::outside_parent_group(
        StereoConfigurationForm::kinded(StereoKind::CisTrans, 0_u32),
        Permutation::from_image(&[1, 2, 0, 3]),
        None,
    )]
    #[case::kinded_undetermined_outside_parent_group(
        StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Undetermined),
        Permutation::from_image(&[1, 2, 0, 3]),
        None,
    )]
    fn test_stereo_configuration_form_reframe_by(
        #[case] input: StereoConfigurationForm,
        #[case] action: Permutation,
        #[case] expected: Option<StereoConfigurationForm>,
    ) {
        assert_eq!(input.reframe_by(&action), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ligand_symmetry(
        StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm {
            permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[0, 2, 1, 3])), orientation: Orientation::Improper },
            invariant: BooleanForm::Lit(true),
        }),
        Permutation::from_image(&[1, 0, 2, 3]),
        StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm {
            permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[2, 1, 0, 3])), orientation: Orientation::Improper },
            invariant: BooleanForm::Lit(true),
        }),
    )]
    #[case::fluxionality(
        StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[0, 2, 1, 3])), active: BooleanForm::Lit(false) }),
        Permutation::from_image(&[1, 0, 2, 3]),
        StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[2, 1, 0, 3])), active: BooleanForm::Lit(false) }),
    )]
    #[case::topicity(
        StereoAtomConstraintForm::Topicity(TopicityForm {
            pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2)),
            relation: TopicityRelationForm::Lit(Topicity::Diastereotopic),
        }),
        Permutation::from_image(&[1, 0, 2, 3]),
        StereoAtomConstraintForm::Topicity(TopicityForm {
            pair: StereoLigandPair::new(StereoLigandPosition(1), StereoLigandPosition(2)),
            relation: TopicityRelationForm::Lit(Topicity::Diastereotopic),
        }),
    )]
    #[case::stereogenicity(
        StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Prochiral)),
        Permutation::identity(3),
        StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Prochiral)),
    )]
    fn test_stereo_atom_constraint_form_reframe_by(
        #[case] input: StereoAtomConstraintForm,
        #[case] action: Permutation,
        #[case] expected: StereoAtomConstraintForm,
    ) {
        assert_eq!(input.reframe_by(&action), Some(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::permutation_degree(StereoAtomConstraintForm::Fluxionality(FluxionalityForm {
        permutation: LigandPermutation(Permutation::identity(4)),
        active: BooleanForm::Lit(true),
    }))]
    #[case::topicity_position(StereoAtomConstraintForm::Topicity(TopicityForm {
        pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(3)),
        relation: TopicityRelationForm::Lit(Topicity::Homotopic),
    }))]
    fn test_stereo_atom_constraint_form_reframe_by_error(
        #[case] input: StereoAtomConstraintForm,
    ) {
        let action = Permutation::identity(3);
        assert_eq!(input.reframe_by(&action), None);
    }

    #[rstest]
    #[case::within_endpoint(
        Permutation::from_image(&[1, 0, 2, 3]),
        Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
            Stereogenicity::Stereogenic,
        ))),
    )]
    #[case::across_endpoints(
        Permutation::from_image(&[1, 2, 0, 3]),
        None,
    )]
    fn test_stereo_bond_constraint_frame_transport_domain(
        #[case] action: Permutation,
        #[case] expected: Option<StereoBondConstraintForm>,
    ) {
        assert_eq!(
            StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
                Stereogenicity::Stereogenic,
            ))
            .reframe_by(&action),
            expected,
        );
    }

    #[rstest]
    fn test_stereo_atom_form_reframe_by_constraints() {
        let action = Permutation::from_image(&[1, 0, 2, 3]);
        let input = StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32).with_constraint(
            StereoAtomConstraintForm::Fluxionality(FluxionalityForm {
                permutation: LigandPermutation(Permutation::from_image(&[0, 2, 1, 3])),
                active: BooleanForm::Lit(true),
            }),
        );
        let expected = StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32).with_constraint(
            StereoAtomConstraintForm::Fluxionality(FluxionalityForm {
                permutation: LigandPermutation(Permutation::from_image(&[2, 1, 0, 3])),
                active: BooleanForm::Lit(true),
            }),
        );

        assert_eq!(input.reframe_by(&action), Some(expected),);
    }

    #[rstest]
    #[case::undetermined(
        StereoCoset::Undetermined,
        StereoKind::Tetrahedral,
        Permutation::identity(4),
        StereoCoset::Undetermined
    )]
    #[case::lit_identity(
        StereoCoset::Lit(0),
        StereoKind::Tetrahedral,
        Permutation::identity(4),
        StereoCoset::Lit(0)
    )]
    #[case::lit_involution(StereoCoset::Lit(0), StereoKind::Tetrahedral, StereoKind::Tetrahedral.involution(), StereoCoset::Lit(1))]
    #[case::lit_set(StereoCoset::lit_set([0]), StereoKind::Tetrahedral, StereoKind::Tetrahedral.involution(), StereoCoset::lit_set([1]))]
    #[case::term_layers(
        StereoCoset::term(StereoTerm::var("x")),
        StereoKind::Tetrahedral,
        Permutation::identity(4),
        StereoCoset::term(StereoTerm::apply(StereoTerm::var("x"), Permutation::identity(4)))
    )]
    fn test_stereo_coset_apply(
        #[case] coset: StereoCoset,
        #[case] kind: StereoKind,
        #[case] permutation: Permutation,
        #[case] expected: StereoCoset,
    ) {
        assert_eq!(coset.apply(kind, permutation), Some(expected));
    }

    #[rstest]
    #[case::undetermined(
        StereoCoset::Undetermined,
        StereoKind::Tetrahedral,
        StereoCoset::Undetermined
    )]
    #[case::tetrahedral_0(StereoCoset::Lit(0), StereoKind::Tetrahedral, StereoCoset::Lit(1))]
    #[case::tetrahedral_1(StereoCoset::Lit(1), StereoKind::Tetrahedral, StereoCoset::Lit(0))]
    #[case::cis_trans(StereoCoset::Lit(0), StereoKind::CisTrans, StereoCoset::Lit(1))]
    #[case::term_layers(
        StereoCoset::term(StereoTerm::var("x")),
        StereoKind::Tetrahedral,
        StereoCoset::term(StereoTerm::swap(StereoTerm::var("x")))
    )]
    fn test_stereo_coset_swap(
        #[case] coset: StereoCoset,
        #[case] kind: StereoKind,
        #[case] expected: StereoCoset,
    ) {
        assert_eq!(coset.swap(kind), Some(expected));
    }

    #[rstest]
    #[case::undetermined(
        StereoCoset::Undetermined,
        StereoKind::Tetrahedral,
        StereoCoset::Undetermined
    )]
    #[case::chiral(StereoCoset::Lit(0), StereoKind::Tetrahedral, StereoCoset::Lit(1))]
    #[case::achiral_noop(StereoCoset::Lit(0), StereoKind::CisTrans, StereoCoset::Lit(0))]
    #[case::term_layers(
        StereoCoset::term(StereoTerm::var("x")),
        StereoKind::Tetrahedral,
        StereoCoset::term(StereoTerm::mirror(StereoTerm::var("x")))
    )]
    fn test_stereo_coset_mirror(
        #[case] coset: StereoCoset,
        #[case] kind: StereoKind,
        #[case] expected: StereoCoset,
    ) {
        assert_eq!(coset.mirror(kind), Some(expected));
    }

    #[rstest]
    #[case::undetermined(
        StereoConfigurationForm::Undetermined,
        Permutation::identity(4),
        StereoConfigurationForm::Undetermined
    )]
    #[case::kinded(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)), StereoKind::Tetrahedral.involution(), StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)))]
    fn test_stereo_configuration_form_apply(
        #[case] config: StereoConfigurationForm,
        #[case] permutation: Permutation,
        #[case] expected: StereoConfigurationForm,
    ) {
        assert_eq!(config.apply(permutation), Some(expected));
    }

    #[rstest]
    #[case::undetermined(
        StereoConfigurationForm::Undetermined,
        StereoConfigurationForm::Undetermined
    )]
    #[case::kinded(
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1))
    )]
    fn test_stereo_configuration_form_swap(
        #[case] config: StereoConfigurationForm,
        #[case] expected: StereoConfigurationForm,
    ) {
        assert_eq!(config.swap(), Some(expected));
    }

    #[rstest]
    #[case::undetermined(
        StereoConfigurationForm::Undetermined,
        StereoConfigurationForm::Undetermined
    )]
    #[case::chiral(
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1))
    )]
    #[case::achiral_noop(
        StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(0)),
        StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(0))
    )]
    fn test_stereo_configuration_form_mirror(
        #[case] config: StereoConfigurationForm,
        #[case] expected: StereoConfigurationForm,
    ) {
        assert_eq!(config.mirror(), Some(expected));
    }

    #[rstest]
    #[case::other_undetermined_keeps_self(
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        StereoConfigurationForm::Undetermined,
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1))
    )]
    #[case::same_kind_undetermined_coset_keeps_coset(
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined),
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1))
    )]
    #[case::same_kind_determined_coset_overrides(
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0))
    )]
    #[case::different_kind_overrides(
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Undetermined),
        StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Undetermined)
    )]
    #[case::self_undetermined_takes_other(
        StereoConfigurationForm::Undetermined,
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined),
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined)
    )]
    fn test_stereo_configuration_form_update(
        #[case] base: StereoConfigurationForm,
        #[case] other: StereoConfigurationForm,
        #[case] expected: StereoConfigurationForm,
    ) {
        assert_eq!(base.update(&other), expected);
    }

    #[rstest]
    #[case::apply(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32), StereoKind::Tetrahedral.involution(), StereoAtomForm::new(StereoKind::Tetrahedral, 1u32))]
    fn test_stereo_atom_form_apply(
        #[case] input: StereoAtomForm,
        #[case] permutation: Permutation,
        #[case] expected: StereoAtomForm,
    ) {
        assert_eq!(input.apply(permutation), Some(expected));
    }

    #[rstest]
    #[case::tetrahedral(
        StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        StereoAtomForm::new(StereoKind::Tetrahedral, 1u32)
    )]
    fn test_stereo_atom_form_swap(#[case] input: StereoAtomForm, #[case] expected: StereoAtomForm) {
        assert_eq!(input.swap(), Some(expected));
    }

    #[rstest]
    #[case::chiral(
        StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        StereoAtomForm::new(StereoKind::Tetrahedral, 1u32)
    )]
    fn test_stereo_atom_form_mirror(
        #[case] input: StereoAtomForm,
        #[case] expected: StereoAtomForm,
    ) {
        assert_eq!(input.mirror(), Some(expected));
    }

    #[rstest]
    #[case::cis_trans(
        StereoBondForm::new(StereoKind::CisTrans, 0u32),
        StereoBondForm::new(StereoKind::CisTrans, 1u32)
    )]
    fn test_stereo_bond_form_swap(#[case] input: StereoBondForm, #[case] expected: StereoBondForm) {
        assert_eq!(input.swap(), Some(expected));
    }

    #[rstest]
    #[case::identity(
        [StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)],
        [StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)],
        0,
    )]
    #[case::transposition(
        [StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)],
        [StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)],
        1,
    )]
    #[case::even_cycle(
        [StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)],
        [StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)],
        0,
    )]
    #[case::virtual_explicit_swap(
        [StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen)],
        [StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
        1,
    )]
    fn test_stereo_atom_form_reframe_by(
        #[case] before: [StereoLigand; 4],
        #[case] after: [StereoLigand; 4],
        #[case] expected_coset: u32,
    ) {
        let permutation = Permutation::between(&before, &after).expect("case has a unique frame");
        assert_eq!(
            StereoAtomForm::new(StereoKind::Tetrahedral, 0u32).reframe_by(&permutation),
            Some(StereoAtomForm::new(StereoKind::Tetrahedral, expected_coset)),
        );
    }

    /// A `Modified` span holds two configurations against one ligand frame, so one action carries
    /// both. Integrity guarantees the sides share a kind, hence a parent group, hence one candidate
    /// set.
    #[rstest]
    fn test_stereo_atom_spans_reframe() {
        let frame: Vec<StereoLigand> = [2, 1, 3, 4]
            .into_iter()
            .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
            .collect();
        let spans = StereoAtomSpans::new(vec![(
            AtomId(0),
            frame.clone(),
            EntitySpan::Modified {
                lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
                rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
            },
        )]);

        let source = spans.clone();
        let (reframed, actions) = spans
            .reframe_with_action()
            .expect("the forms are satisfiable");
        let action = actions
            .action(StereoAtomId(0))
            .expect("the dense action covers the stereo atom");

        assert_eq!(
            reframed.ligands(StereoAtomId(0)),
            [1, 2, 3, 4]
                .into_iter()
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                .collect::<Vec<_>>(),
        );
        // The same action carried both sides: neither was reframed independently.
        assert_eq!(
            reframed.attributes(StereoAtomId(0)),
            &EntitySpan::Modified {
                lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32)
                    .reframe_by(action)
                    .expect("a parent-group action"),
                rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 1u32)
                    .reframe_by(action)
                    .expect("a parent-group action"),
            },
        );
        assert_eq!(source.reframe_by(&actions), Some(reframed));
    }

    #[rstest]
    fn test_stereo_atom_spans_normalize() {
        let spans = StereoAtomSpans::new(vec![(
            AtomId(0),
            [1, 2, 3, 4]
                .into_iter()
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                .collect(),
            EntitySpan::Modified {
                lhs: StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::lit_set([1])),
                rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
            },
        )]);

        let normalized = spans.normalize().expect("the forms are satisfiable");

        assert_eq!(
            normalized.attributes(StereoAtomId(0)),
            &EntitySpan::Unchanged(StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32)),
        );
    }

    #[rstest]
    fn test_stereo_atom_spans_reframe_identity() {
        let spans = StereoAtomSpans::new(vec![(
            AtomId(0),
            [4, 2, 3, 1]
                .into_iter()
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                .collect(),
            EntitySpan::Modified {
                lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
                rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
            },
        )]);
        let once = spans.reframe().expect("the forms are satisfiable");
        let twice = once.clone().reframe().expect("the forms are satisfiable");
        assert_eq!(twice, once);
    }

    /// One side declining takes the whole span: a single action serves both, so there is no partial
    /// result to keep. Here the rhs asserts a coset outside its kind's coset space.
    #[rstest]
    fn test_stereo_atom_spans_reframe_error() {
        let spans = StereoAtomSpans::new(vec![(
            AtomId(0),
            [1, 2, 3, 4]
                .into_iter()
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                .collect(),
            EntitySpan::Modified {
                lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
                rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 5u32),
            },
        )]);
        assert_eq!(spans.reframe(), Err(Contradiction));
    }

    /// The ligands are an `Ordered` factor, so the supplied frame is the stored frame and no
    /// fixture surgery is needed to obtain an unselected one.
    #[rstest]
    fn test_stereo_atoms_reframe() {
        let atoms = StereoAtoms::new(vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        )]);

        let reframed = atoms.reframe().expect("the form is satisfiable");

        assert_eq!(
            reframed.ligands(StereoAtomId(0)),
            [
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
        );
        assert_eq!(reframed.site(StereoAtomId(0)), AtomId(0));
    }

    #[rstest]
    fn test_stereo_atoms_reframe_identity() {
        let atoms = StereoAtoms::new(vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
        )]);
        let once = atoms.reframe().expect("the form is satisfiable");
        let twice = once.clone().reframe().expect("the form is satisfiable");
        assert_eq!(twice, once);
    }

    #[rstest]
    fn test_stereo_atoms_reframe_with_action() {
        let atoms = StereoAtoms::new(vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        )]);

        let (reframed, actions) = atoms
            .clone()
            .reframe_with_action()
            .expect("the form is satisfiable");

        let action = actions
            .action(StereoAtomId(0))
            .expect("the dense action covers the stereo atom");
        assert_eq!(action, &Permutation::from_image(&[3, 2, 1, 0]));
        assert_eq!(atoms.reframe_by(&actions), Some(reframed));
    }

    #[rstest]
    fn test_stereo_atoms_normalize() {
        let atoms = StereoAtoms::new(vec![(
            AtomId(0),
            [1, 2, 3, 4]
                .into_iter()
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                .collect(),
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::lit_set([1])),
        )]);

        let normalized = atoms.normalize().expect("the form is satisfiable");

        assert_eq!(
            normalized.attributes(StereoAtomId(0)),
            &StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
        );
    }

    #[rstest]
    fn test_stereo_atoms_framed_eq_distinct_frames() {
        let frame = [2, 1, 3, 4]
            .into_iter()
            .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
            .collect::<Vec<_>>();
        let presentation = Permutation::from_image(&[1, 0, 3, 2]);
        let form = StereoAtomForm::new(StereoKind::Tetrahedral, 1u32);
        let stored = StereoAtoms::new(vec![(AtomId(0), frame.clone(), form.clone())]);
        let restated = StereoAtoms::new(vec![(
            AtomId(0),
            presentation.act(&frame),
            form.reframe_by(&presentation)
                .expect("a parent-group action"),
        )]);

        assert!(stored.framed_eq(&restated));

        let other_coset = StereoAtoms::new(vec![(
            AtomId(0),
            frame,
            StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        )]);
        assert!(!stored.framed_eq(&other_coset));
    }

    /// A cis/trans span keeps each endpoint's ligands in its own block, so the partitioned parent
    /// group sorts within blocks; one action still carries both sides.
    #[rstest]
    fn test_stereo_bond_spans_reframe() {
        let frame: Vec<StereoLigand> = [5, 3, 4, 2]
            .into_iter()
            .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
            .collect();
        let spans = StereoBondSpans::new(vec![(
            BondId(0),
            frame,
            EntitySpan::Modified {
                lhs: StereoBondForm::new(StereoKind::CisTrans, 0u32),
                rhs: StereoBondForm::new(StereoKind::CisTrans, 1u32),
            },
        )]);

        let source = spans.clone();
        let (reframed, actions) = spans
            .reframe_with_action()
            .expect("the forms are satisfiable");
        let action = actions
            .action(StereoBondId(0))
            .expect("the dense action covers the stereo bond");

        assert_eq!(
            reframed.ligands(StereoBondId(0)),
            [2, 4, 3, 5]
                .into_iter()
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                .collect::<Vec<_>>(),
        );
        assert_eq!(
            reframed.attributes(StereoBondId(0)),
            &EntitySpan::Modified {
                lhs: StereoBondForm::new(StereoKind::CisTrans, 0u32)
                    .reframe_by(action)
                    .expect("a parent-group action"),
                rhs: StereoBondForm::new(StereoKind::CisTrans, 1u32)
                    .reframe_by(action)
                    .expect("a parent-group action"),
            },
        );
        assert_eq!(source.reframe_by(&actions), Some(reframed));
    }

    #[rstest]
    fn test_stereo_bond_spans_normalize() {
        let spans = StereoBondSpans::new(vec![(
            BondId(0),
            [1, 2, 3, 4]
                .into_iter()
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                .collect(),
            EntitySpan::Modified {
                lhs: StereoBondForm::new(StereoKind::CisTrans, StereoCoset::lit_set([1])),
                rhs: StereoBondForm::new(StereoKind::CisTrans, 1_u32),
            },
        )]);

        let normalized = spans.normalize().expect("the forms are satisfiable");

        assert_eq!(
            normalized.attributes(StereoBondId(0)),
            &EntitySpan::Unchanged(StereoBondForm::new(StereoKind::CisTrans, 1_u32)),
        );
    }

    #[rstest]
    fn test_stereo_bond_spans_reframe_identity() {
        let spans = StereoBondSpans::new(vec![(
            BondId(0),
            [4, 2, 5, 3]
                .into_iter()
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                .collect(),
            EntitySpan::Modified {
                lhs: StereoBondForm::new(StereoKind::CisTrans, 0u32),
                rhs: StereoBondForm::new(StereoKind::CisTrans, 1u32),
            },
        )]);
        let once = spans.reframe().expect("the forms are satisfiable");
        let twice = once.clone().reframe().expect("the forms are satisfiable");
        assert_eq!(twice, once);
    }

    /// A cis/trans frame keeps each endpoint's two ligands in its own block, so the partitioned
    /// parent group sorts within blocks and orders the blocks rather than sorting outright.
    #[rstest]
    fn test_stereo_bonds_reframe() {
        let bonds = StereoBonds::new(vec![(
            BondId(0),
            vec![
                StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            ],
            StereoBondForm::new(StereoKind::CisTrans, 0u32),
        )]);

        let reframed = bonds.reframe().expect("the form is satisfiable");

        assert_eq!(
            reframed.ligands(StereoBondId(0)),
            [
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
            ],
        );
        assert_eq!(reframed.site(StereoBondId(0)), BondId(0));
    }

    #[rstest]
    fn test_stereo_bonds_reframe_identity() {
        let bonds = StereoBonds::new(vec![(
            BondId(0),
            vec![
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            ],
            StereoBondForm::new(StereoKind::CisTrans, 1u32),
        )]);
        let once = bonds.reframe().expect("the form is satisfiable");
        let twice = once.clone().reframe().expect("the form is satisfiable");
        assert_eq!(twice, once);
    }

    #[rstest]
    fn test_stereo_bonds_reframe_with_action() {
        let bonds = StereoBonds::new(vec![(
            BondId(0),
            vec![
                StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            ],
            StereoBondForm::new(StereoKind::CisTrans, 0u32),
        )]);

        let (reframed, actions) = bonds
            .clone()
            .reframe_with_action()
            .expect("the form is satisfiable");

        let action = actions
            .action(StereoBondId(0))
            .expect("the dense action covers the stereo bond");
        assert_eq!(action, &Permutation::from_image(&[3, 2, 1, 0]));
        assert_eq!(bonds.reframe_by(&actions), Some(reframed));
    }

    #[rstest]
    fn test_stereo_bonds_normalize() {
        let bonds = StereoBonds::new(vec![(
            BondId(0),
            [1, 2, 3, 4]
                .into_iter()
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                .collect(),
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::lit_set([1])),
        )]);

        let normalized = bonds.normalize().expect("the form is satisfiable");

        assert_eq!(
            normalized.attributes(StereoBondId(0)),
            &StereoBondForm::new(StereoKind::CisTrans, 1_u32),
        );
    }

    #[rstest]
    fn test_stereo_bonds_framed_eq_distinct_frames() {
        let frame = [2, 4, 3, 5]
            .into_iter()
            .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
            .collect::<Vec<_>>();
        let presentation = Permutation::from_image(&[1, 0, 3, 2]);
        let form = StereoBondForm::new(StereoKind::CisTrans, 1u32);
        let stored = StereoBonds::new(vec![(BondId(0), frame.clone(), form.clone())]);
        let restated = StereoBonds::new(vec![(
            BondId(0),
            presentation.act(&frame),
            form.reframe_by(&presentation)
                .expect("a parent-group action"),
        )]);

        assert!(stored.framed_eq(&restated));

        let other_coset = StereoBonds::new(vec![(
            BondId(0),
            frame,
            StereoBondForm::new(StereoKind::CisTrans, 0u32),
        )]);
        assert!(!stored.framed_eq(&other_coset));
    }

    /// A form's constructors are permissive, so a coset index outside its kind's coset space is
    /// representable until molecule publication checks it. The frame action is the first operation
    /// that requires the index to be in range, so it declines there.
    #[rstest]
    #[case::atom_coset_out_of_range(StereoKind::Tetrahedral, 2u32)]
    #[case::atom_coset_far_out_of_range(StereoKind::SquarePlanar, 40u32)]
    fn test_stereo_atom_form_reframe_by_error(#[case] kind: StereoKind, #[case] coset: u32) {
        let form = StereoAtomForm::new(kind, coset);
        assert_eq!(form.clone().reframe_by(&Permutation::identity(4)), None);
    }

    #[rstest]
    #[case::degree(
        StereoBondForm::new(StereoKind::CisTrans, 0u32),
        Permutation::identity(3)
    )]
    #[case::kinded_outside_parent(
        StereoBondForm::new(StereoKind::CisTrans, 0u32),
        Permutation::from_image(&[1, 2, 0, 3]),
    )]
    #[case::kindless_outside_parent(
        StereoBondForm::default(),
        Permutation::from_image(&[1, 2, 0, 3]),
    )]
    fn test_stereo_bond_form_frame_transport_domain(
        #[case] form: StereoBondForm,
        #[case] permutation: Permutation,
    ) {
        assert_eq!(form.reframe_by(&permutation), None,);
    }

    #[rstest]
    fn test_stereo_atom_form_reframe_by_roundtrip() {
        let before = [
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
        ];
        let after = [
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
        ];
        let permutation = Permutation::between(&before, &after).expect("a unique frame change");
        let atom = StereoAtomForm::new(StereoKind::Tetrahedral, 0u32);
        assert_eq!(
            atom.clone()
                .reframe_by(&permutation)
                .and_then(|reframed| reframed.reframe_by(&permutation.inverse())),
            Some(atom),
        );
    }

    /// A within-endpoint swap is an action of the cis/trans parent group and exchanges the two
    /// cosets.
    #[rstest]
    fn test_stereo_bond_form_reframe_by() {
        let before = [
            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
        ];
        let after = [
            StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
        ];
        let permutation = Permutation::between(&before, &after).expect("a unique frame change");
        assert_eq!(
            StereoBondForm::new(StereoKind::CisTrans, 0u32).reframe_by(&permutation),
            Some(StereoBondForm::new(StereoKind::CisTrans, 1u32)),
        );
    }
}
