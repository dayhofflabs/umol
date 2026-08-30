//! Aggregate canonicalization inputs and failures.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::hash::{DefaultHasher, Hash, Hasher};

use thiserror::Error;
use umol_graph_core::{
    AutomorphismAlgorithm, AutomorphismOutput, Correspondence, Graph, NodeId,
    SubdivisionNodeSource, UnionFind,
};
use umol_perm::{Orientation, Permutation};

use super::atom::{AtomForm, ElementForm, IsotopeMassForm};
use super::bond::BondForm;
use super::boolean::BooleanForm;
use super::constraint::{
    AromaticSystemConstraintForm, AromaticValenceForm, AtomConstraintForm, BondConstraintForm,
    Constraint, DativeBondConstraintForm, LigandPermutation, MoleculeConstraint,
    MulticenterBondConstraintForm, MulticenterValenceForm, NoncovalentBondConstraintForm,
    OrientedLigandPermutation, RelationalConstraint, RingMembershipForm, RingScope,
    StereoAtomConstraintForm, StereoBondConstraintForm, StereoLigandPair, StereogenicityForm,
    TopicityRelationForm,
};
use super::correspondence::MoleculeCorrespondence;
use super::delta::{
    AromaticSystemDelta, AtomDelta, BondDelta, ConstraintSpan, DativeBondDelta, Delta, EntitySpan,
    MulticenterBondDelta, NoncovalentBondDelta, StereoAtomDelta, StereoBondDelta,
};
use super::electrons::ElectronCountsForm;
use super::entity::{Entity, EntityKind};
use super::error::Contradiction;
use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::incidence::{Incidence, IncidenceGraph, IncidenceLevel};
use super::ligand::{StereoLigand, StereoLigandKind};
use super::molecule::Molecule;
use super::noncovalent::{NoncovalentBondKind, NoncovalentBondKindForm};
use super::num::{ArithExpr, NumForm, PredExpr};
use super::operators::{MemOp, RelOp};
use super::reaction::Reaction;
use super::reaction_span::ReactionSpan;
use super::spin::UnpairedElectronsForm;
use super::stereo::{
    CisTransStereoForm, StereoAtomForm, StereoBondForm, StereoConfigurationForm, StereoCoset,
    StereoKind, StereoTerm, Stereogenicity, TetrahedralStereoForm, Topicity,
};
use super::traits::{FrameTransport, Normalize, Reframe};

/// Semantic and operational inputs to aggregate canonicalization.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CanonicalizeContext {
    /// Whether stereo-sensitive refinement is iterated to a para-stereo fixpoint.
    pub para_stereo: bool,
    /// Graph automorphism algorithm used during canonical-frame search.
    pub automorphism_algorithm: AutomorphismAlgorithm,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum DescriptionLevel {
    Topology,
    Constitution,
    Structure,
    Full,
}

fn canonicalize_level_with_constraints(
    base: DescriptionLevel,
    has_inline_constraints: bool,
) -> DescriptionLevel {
    if has_inline_constraints {
        DescriptionLevel::Full
    } else {
        base
    }
}

fn entity_span_canonicalize_level<T>(
    span: &EntitySpan<T>,
    base: DescriptionLevel,
    has_inline_constraints: impl Fn(&T) -> bool,
) -> DescriptionLevel {
    let has_inline_constraints = match span {
        EntitySpan::Unchanged(attributes)
        | EntitySpan::Added(attributes)
        | EntitySpan::Removed(attributes) => has_inline_constraints(attributes),
        EntitySpan::Modified { lhs, rhs } => {
            has_inline_constraints(lhs) || has_inline_constraints(rhs)
        }
    };
    canonicalize_level_with_constraints(base, has_inline_constraints)
}

fn delta_canonicalize_level(delta: &Delta) -> DescriptionLevel {
    match delta {
        Delta::Atom(AtomDelta::Add { attributes, .. } | AtomDelta::Remove { attributes, .. }) => {
            canonicalize_level_with_constraints(
                DescriptionLevel::Topology,
                !attributes.constraints.is_empty(),
            )
        }
        Delta::Atom(AtomDelta::ModifyField { .. }) => DescriptionLevel::Topology,
        Delta::Atom(AtomDelta::ModifyConstraint { .. }) => DescriptionLevel::Full,
        Delta::Bond(BondDelta::Add { attributes, .. } | BondDelta::Remove { attributes, .. }) => {
            canonicalize_level_with_constraints(
                DescriptionLevel::Topology,
                !attributes.constraints.is_empty(),
            )
        }
        Delta::Bond(BondDelta::ModifyField { .. }) => DescriptionLevel::Topology,
        Delta::Bond(BondDelta::ModifyConstraint { .. }) => DescriptionLevel::Full,
        Delta::DativeBond(
            DativeBondDelta::Add { attributes, .. } | DativeBondDelta::Remove { attributes, .. },
        ) => canonicalize_level_with_constraints(
            DescriptionLevel::Constitution,
            !attributes.constraints.is_empty(),
        ),
        Delta::DativeBond(DativeBondDelta::ModifyField { .. }) => DescriptionLevel::Constitution,
        Delta::DativeBond(DativeBondDelta::ModifyConstraint { .. }) => DescriptionLevel::Full,
        Delta::AromaticSystem(
            AromaticSystemDelta::Add { attributes, .. }
            | AromaticSystemDelta::Remove { attributes, .. },
        ) => canonicalize_level_with_constraints(
            DescriptionLevel::Constitution,
            !attributes.constraints.is_empty(),
        ),
        Delta::AromaticSystem(AromaticSystemDelta::ModifyField { .. }) => {
            DescriptionLevel::Constitution
        }
        Delta::AromaticSystem(AromaticSystemDelta::ModifyConstraint { .. }) => {
            DescriptionLevel::Full
        }
        Delta::MulticenterBond(
            MulticenterBondDelta::Add { attributes, .. }
            | MulticenterBondDelta::Remove { attributes, .. },
        ) => canonicalize_level_with_constraints(
            DescriptionLevel::Constitution,
            !attributes.constraints.is_empty(),
        ),
        Delta::MulticenterBond(MulticenterBondDelta::ModifyField { .. }) => {
            DescriptionLevel::Constitution
        }
        Delta::MulticenterBond(MulticenterBondDelta::ModifyConstraint { .. }) => {
            DescriptionLevel::Full
        }
        Delta::NoncovalentBond(
            NoncovalentBondDelta::Add { attributes, .. }
            | NoncovalentBondDelta::Remove { attributes, .. },
        ) => canonicalize_level_with_constraints(
            DescriptionLevel::Constitution,
            !attributes.constraints.is_empty(),
        ),
        Delta::NoncovalentBond(NoncovalentBondDelta::ModifyField { .. }) => {
            DescriptionLevel::Constitution
        }
        Delta::NoncovalentBond(NoncovalentBondDelta::ModifyConstraint { .. }) => {
            DescriptionLevel::Full
        }
        Delta::StereoAtom(
            StereoAtomDelta::Add { attributes, .. } | StereoAtomDelta::Remove { attributes, .. },
        ) => canonicalize_level_with_constraints(
            DescriptionLevel::Structure,
            !attributes.constraints.is_empty(),
        ),
        Delta::StereoAtom(StereoAtomDelta::ModifyField { .. }) => DescriptionLevel::Structure,
        Delta::StereoAtom(StereoAtomDelta::ModifyConstraint { .. }) => DescriptionLevel::Full,
        Delta::StereoBond(
            StereoBondDelta::Add { attributes, .. } | StereoBondDelta::Remove { attributes, .. },
        ) => canonicalize_level_with_constraints(
            DescriptionLevel::Structure,
            !attributes.constraints.is_empty(),
        ),
        Delta::StereoBond(StereoBondDelta::ModifyField { .. }) => DescriptionLevel::Structure,
        Delta::StereoBond(StereoBondDelta::ModifyConstraint { .. }) => DescriptionLevel::Full,
        Delta::Constraint(_) => DescriptionLevel::Full,
    }
}

fn molecule_canonicalize_level(molecule: &Molecule) -> DescriptionLevel {
    let has_inline_constraints = molecule
        .atoms()
        .iter()
        .any(|atom| !atom.attributes.constraints.is_empty())
        || molecule
            .bonds()
            .iter()
            .any(|bond| !bond.attributes.constraints.is_empty())
        || molecule
            .dative_bonds()
            .iter()
            .any(|bond| !bond.attributes.constraints.is_empty())
        || molecule
            .aromatic_systems()
            .iter()
            .any(|system| !system.attributes.constraints.is_empty())
        || molecule
            .multicenter_bonds()
            .iter()
            .any(|bond| !bond.attributes.constraints.is_empty())
        || molecule
            .noncovalent_bonds()
            .iter()
            .any(|bond| !bond.attributes.constraints.is_empty())
        || molecule
            .stereo_atoms()
            .iter()
            .any(|stereo| !stereo.attributes.constraints.is_empty())
        || molecule
            .stereo_bonds()
            .iter()
            .any(|stereo| !stereo.attributes.constraints.is_empty());

    if has_inline_constraints || molecule.has_constraints() {
        DescriptionLevel::Full
    } else if molecule.has_stereo_atoms() || molecule.has_stereo_bonds() {
        DescriptionLevel::Structure
    } else if molecule.has_dative_bonds()
        || molecule.has_aromatic_systems()
        || molecule.has_multicenter_bonds()
        || molecule.has_noncovalent_bonds()
    {
        DescriptionLevel::Constitution
    } else {
        DescriptionLevel::Topology
    }
}

fn reaction_canonicalize_level(reaction: &Reaction) -> DescriptionLevel {
    reaction.deltas().iter().fold(
        molecule_canonicalize_level(reaction.lhs()),
        |level, delta| level.max(delta_canonicalize_level(delta)),
    )
}

fn reaction_span_canonicalize_level(span: &ReactionSpan) -> DescriptionLevel {
    if !span.constraints().is_empty() {
        return DescriptionLevel::Full;
    }

    let mut level = DescriptionLevel::Topology;
    for entity in span.atoms() {
        level = level.max(entity_span_canonicalize_level(
            entity,
            DescriptionLevel::Topology,
            |attributes| !attributes.constraints.is_empty(),
        ));
    }
    for entity in span.bonds() {
        level = level.max(entity_span_canonicalize_level(
            entity,
            DescriptionLevel::Topology,
            |attributes| !attributes.constraints.is_empty(),
        ));
    }
    for id in span.dative_bonds().ids() {
        level = level.max(entity_span_canonicalize_level(
            span.dative_bonds().attributes(id),
            DescriptionLevel::Constitution,
            |attributes| !attributes.constraints.is_empty(),
        ));
    }
    for id in span.aromatic_systems().ids() {
        level = level.max(entity_span_canonicalize_level(
            span.aromatic_systems().attributes(id),
            DescriptionLevel::Constitution,
            |attributes| !attributes.constraints.is_empty(),
        ));
    }
    for id in span.multicenter_bonds().ids() {
        level = level.max(entity_span_canonicalize_level(
            span.multicenter_bonds().attributes(id),
            DescriptionLevel::Constitution,
            |attributes| !attributes.constraints.is_empty(),
        ));
    }
    for id in span.noncovalent_bonds().ids() {
        level = level.max(entity_span_canonicalize_level(
            span.noncovalent_bonds().attributes(id),
            DescriptionLevel::Constitution,
            |attributes| !attributes.constraints.is_empty(),
        ));
    }
    for id in span.stereo_atoms().ids() {
        level = level.max(entity_span_canonicalize_level(
            span.stereo_atoms().attributes(id),
            DescriptionLevel::Structure,
            |attributes| !attributes.constraints.is_empty(),
        ));
    }
    for id in span.stereo_bonds().ids() {
        level = level.max(entity_span_canonicalize_level(
            span.stereo_bonds().attributes(id),
            DescriptionLevel::Structure,
            |attributes| !attributes.constraints.is_empty(),
        ));
    }
    level
}

/// Canonical entity-frame selection for complete indexed graph-IR aggregates.
///
/// Canonicalization extends [`Reframe`]: normalization is the fixed-frame prefix, reframing adds
/// participant-frame selection, and canonicalization additionally selects entity ids. It preserves
/// the represented aggregate, transports every reference and position-sensitive value, and
/// normalizes carried forms in the selected frame. Construction of the canonical form is fallible;
/// canonical equality totalizes those failures according to the aggregate's semantic contract.
///
/// # Semantic properties
///
/// For every integrity-valid aggregate and fixed context:
///
/// - successful complete canonicalization is exactly idempotent and invariant under valid dense
///   entity remapping;
/// - `canonical_eq` is reflexive, symmetric, and transitive under its documented failure
///   totalization;
/// - successful canonical hashes are invariant under valid dense entity remapping, and canonical
///   equality implies equal canonical hashes when both hash operations succeed.
///
/// For a fixed umol release and context, canonicalization is deterministic. During the 0.x series,
/// the typed comparison schema and resulting representatives may change between releases. Returned
/// canonical forms are ordinary IR values without schema-version provenance and must not be used as
/// persistent cross-release identifiers.
pub trait Canonicalize: Reframe {
    type Error;

    /// Construct the complete canonical form.
    ///
    /// # Errors
    ///
    /// Returns the aggregate-specific canonicalization error when intrinsic normalization is
    /// unsatisfiable.
    fn canonicalize(self, context: &CanonicalizeContext) -> Result<Self, Self::Error>;

    /// Construct the complete canonical form and its source-to-canonical correspondence.
    ///
    /// The correspondence maps every entity id in the input frame to its id in the returned
    /// canonical frame. Each of its eight entity-kind components is total on both sides and
    /// therefore represents a dense bijection. For [`Reaction`], the two frames are the complete
    /// union frames of the materialized input and returned reaction spans.
    ///
    /// # Errors
    ///
    /// Returns the same aggregate-specific intrinsic-normalization error as
    /// [`Self::canonicalize`]. A [`Reaction`] also reports failure to materialize its reaction span.
    ///
    /// # Semantic properties
    ///
    /// Discarding the correspondence yields exactly [`Self::canonicalize`] under the same context.
    /// Remapping the input through the correspondence and then applying [`Reframe::reframe`]
    /// reconstructs the returned canonical value.
    fn canonicalize_with_correspondence(
        self,
        context: &CanonicalizeContext,
    ) -> Result<(Self, MoleculeCorrespondence), Self::Error>;

    /// Hash the complete canonical form.
    ///
    /// The returned value uses Rust's [`DefaultHasher`]. It may collide and is not a persistent
    /// identifier: both the hasher and umol's canonical comparison schema may change between
    /// releases.
    ///
    /// # Errors
    ///
    /// Returns the same error as [`Self::canonicalize`].
    ///
    /// # Semantic properties
    ///
    /// Successful hashes are invariant under valid dense entity remapping. Two successfully
    /// canonicalized values that are equal have equal canonical hashes.
    fn canonical_hash(self, context: &CanonicalizeContext) -> Result<u64, Self::Error>
    where
        Self: Hash,
    {
        self.canonicalize(context)
            .map(|canonical| hash_value(&canonical))
    }

    /// Compare complete canonical forms.
    ///
    /// Structural identity short-circuits the operation. Otherwise, two intrinsic contradictions
    /// compare equal.
    fn canonical_eq(&self, other: &Self, context: &CanonicalizeContext) -> bool;
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct StructuralDomainPosition(u16);

impl StructuralDomainPosition {
    const TOPOLOGY: Self = Self(0);
    const NON_STEREO: Self = Self(1);
    const STEREO: Self = Self(2);
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct EntityBlockPosition {
    domain: StructuralDomainPosition,
    slot: u16,
}

impl EntityBlockPosition {
    const ATOM: Self = Self::new(StructuralDomainPosition::TOPOLOGY, 0);
    const BOND: Self = Self::new(StructuralDomainPosition::TOPOLOGY, 1);
    const DATIVE_BOND: Self = Self::new(StructuralDomainPosition::NON_STEREO, 0);
    const AROMATIC_SYSTEM: Self = Self::new(StructuralDomainPosition::NON_STEREO, 1);
    const MULTICENTER_BOND: Self = Self::new(StructuralDomainPosition::NON_STEREO, 2);
    const NONCOVALENT_BOND: Self = Self::new(StructuralDomainPosition::NON_STEREO, 3);
    const STEREO_ATOM: Self = Self::new(StructuralDomainPosition::STEREO, 0);
    const STEREO_BOND: Self = Self::new(StructuralDomainPosition::STEREO, 1);

    const fn new(domain: StructuralDomainPosition, slot: u16) -> Self {
        Self { domain, slot }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct FieldPosition(u16);

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct VariantPosition(u16);

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct SpanTagPosition(u16);

impl SpanTagPosition {
    const UNCHANGED: Self = Self(0);
    const ADDED: Self = Self(1);
    const REMOVED: Self = Self(2);
    const MODIFIED: Self = Self(3);
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum ConstraintBlockPosition {
    Inline(EntityBlockPosition),
    Molecule,
}

impl ConstraintBlockPosition {
    const ATOM: Self = Self::Inline(EntityBlockPosition::ATOM);
    const BOND: Self = Self::Inline(EntityBlockPosition::BOND);
    const DATIVE_BOND: Self = Self::Inline(EntityBlockPosition::DATIVE_BOND);
    const AROMATIC_SYSTEM: Self = Self::Inline(EntityBlockPosition::AROMATIC_SYSTEM);
    const MULTICENTER_BOND: Self = Self::Inline(EntityBlockPosition::MULTICENTER_BOND);
    const NONCOVALENT_BOND: Self = Self::Inline(EntityBlockPosition::NONCOVALENT_BOND);
    const STEREO_ATOM: Self = Self::Inline(EntityBlockPosition::STEREO_ATOM);
    const STEREO_BOND: Self = Self::Inline(EntityBlockPosition::STEREO_BOND);
    const MOLECULE: Self = Self::Molecule;
}

impl Ord for ConstraintBlockPosition {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (Self::Inline(lhs), Self::Inline(rhs)) => lhs.cmp(rhs),
            (Self::Inline(_), Self::Molecule) => Ordering::Less,
            (Self::Molecule, Self::Inline(_)) => Ordering::Greater,
            (Self::Molecule, Self::Molecule) => Ordering::Equal,
        }
    }
}

impl PartialOrd for ConstraintBlockPosition {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct RelationalConstraintPosition {
    entity: EntityBlockPosition,
    slot: u16,
}

impl RelationalConstraintPosition {
    const fn new(entity: EntityBlockPosition, slot: u16) -> Self {
        Self { entity, slot }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct PositionedKey<P> {
    position: P,
    value: CanonicalKeyValue,
}

impl<P: Ord> Ord for PositionedKey<P> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.position
            .cmp(&other.position)
            .then_with(|| self.value.cmp(&other.value))
    }
}

impl<P: Ord> PartialOrd for PositionedKey<P> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

type FieldKey = PositionedKey<FieldPosition>;
type EntityBlockKey = PositionedKey<EntityBlockPosition>;
type ConstraintBlockKey = PositionedKey<ConstraintBlockPosition>;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct VariantKey {
    position: VariantPosition,
    fields: Vec<FieldKey>,
}

impl Ord for VariantKey {
    fn cmp(&self, other: &Self) -> Ordering {
        self.position
            .cmp(&other.position)
            .then_with(|| self.fields.cmp(&other.fields))
    }
}

impl PartialOrd for VariantKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct SpanKey {
    position: SpanTagPosition,
    values: Vec<CanonicalKeyValue>,
}

impl Ord for SpanKey {
    fn cmp(&self, other: &Self) -> Ordering {
        self.position
            .cmp(&other.position)
            .then_with(|| self.values.cmp(&other.values))
    }
}

impl PartialOrd for SpanKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
enum CanonicalKeyValue {
    Bool(bool),
    Unsigned(u64),
    Signed(i64),
    Text(String),
    Sequence(Vec<Self>),
    Product(Vec<FieldKey>),
    Variant(VariantKey),
    Span(SpanKey),
}

impl CanonicalKeyValue {
    fn position(&self) -> u16 {
        match self {
            Self::Bool(_) => 1,
            Self::Unsigned(_) => 2,
            Self::Signed(_) => 3,
            Self::Text(_) => 4,
            Self::Sequence(_) => 5,
            Self::Product(_) => 6,
            Self::Variant(_) => 7,
            Self::Span(_) => 8,
        }
    }
}

impl Ord for CanonicalKeyValue {
    fn cmp(&self, other: &Self) -> Ordering {
        self.position()
            .cmp(&other.position())
            .then_with(|| match (self, other) {
                (Self::Bool(lhs), Self::Bool(rhs)) => lhs.cmp(rhs),
                (Self::Unsigned(lhs), Self::Unsigned(rhs)) => lhs.cmp(rhs),
                (Self::Signed(lhs), Self::Signed(rhs)) => lhs.cmp(rhs),
                (Self::Text(lhs), Self::Text(rhs)) => lhs.cmp(rhs),
                (Self::Sequence(lhs), Self::Sequence(rhs)) => lhs.cmp(rhs),
                (Self::Product(lhs), Self::Product(rhs)) => lhs.cmp(rhs),
                (Self::Variant(lhs), Self::Variant(rhs)) => lhs.cmp(rhs),
                (Self::Span(lhs), Self::Span(rhs)) => lhs.cmp(rhs),
                _ => Ordering::Equal,
            })
    }
}

impl PartialOrd for CanonicalKeyValue {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct CanonicalComparisonKey {
    entity_blocks: Vec<EntityBlockKey>,
    constraints: Vec<ConstraintBlockKey>,
}

impl Ord for CanonicalComparisonKey {
    fn cmp(&self, other: &Self) -> Ordering {
        self.entity_blocks
            .cmp(&other.entity_blocks)
            .then_with(|| self.constraints.cmp(&other.constraints))
    }
}

impl PartialOrd for CanonicalComparisonKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct CanonicalComparisonKeyPrefix {
    entity_blocks: Vec<EntityBlockKey>,
    constraints: Vec<ConstraintBlockKey>,
}

impl CanonicalComparisonKeyPrefix {
    fn cmp_key(&self, key: &CanonicalComparisonKey) -> Ordering {
        fn cmp_key_block_prefixes<P: Ord>(
            prefixes: &[PositionedKey<P>],
            blocks: &[PositionedKey<P>],
        ) -> Ordering {
            for (index, prefix) in prefixes.iter().enumerate() {
                let Some(block) = blocks.get(index) else {
                    return Ordering::Greater;
                };
                let position_order = prefix.position.cmp(&block.position);
                if position_order != Ordering::Equal {
                    return position_order;
                }
                let CanonicalKeyValue::Sequence(prefix_rows) = &prefix.value else {
                    unreachable!("canonical key block prefixes contain row sequences")
                };
                let CanonicalKeyValue::Sequence(rows) = &block.value else {
                    unreachable!("canonical key blocks contain row sequences")
                };
                for (row_index, prefix_row) in prefix_rows.iter().enumerate() {
                    let Some(row) = rows.get(row_index) else {
                        return Ordering::Greater;
                    };
                    let row_order = prefix_row.cmp(row);
                    if row_order != Ordering::Equal {
                        return row_order;
                    }
                }
            }
            Ordering::Equal
        }

        cmp_key_block_prefixes(&self.entity_blocks, &key.entity_blocks)
            .then_with(|| cmp_key_block_prefixes(&self.constraints, &key.constraints))
    }
}

fn hash_value<T: Hash>(value: &T) -> u64 {
    let mut hasher = DefaultHasher::new();
    value.hash(&mut hasher);
    hasher.finish()
}

fn field(position: u16, value: CanonicalKeyValue) -> FieldKey {
    PositionedKey {
        position: FieldPosition(position),
        value,
    }
}

fn product(values: impl IntoIterator<Item = CanonicalKeyValue>) -> CanonicalKeyValue {
    CanonicalKeyValue::Product(
        values
            .into_iter()
            .enumerate()
            .map(|(position, value)| field(position as u16, value))
            .collect(),
    )
}

fn positioned_product(
    values: impl IntoIterator<Item = (u16, CanonicalKeyValue)>,
) -> CanonicalKeyValue {
    CanonicalKeyValue::Product(
        values
            .into_iter()
            .map(|(position, value)| field(position, value))
            .collect(),
    )
}

fn variant(
    position: u16,
    values: impl IntoIterator<Item = CanonicalKeyValue>,
) -> CanonicalKeyValue {
    CanonicalKeyValue::Variant(VariantKey {
        position: VariantPosition(position),
        fields: values
            .into_iter()
            .enumerate()
            .map(|(position, value)| field(position as u16, value))
            .collect(),
    })
}

fn sequence(values: impl IntoIterator<Item = CanonicalKeyValue>) -> CanonicalKeyValue {
    CanonicalKeyValue::Sequence(values.into_iter().collect())
}

fn option(value: Option<CanonicalKeyValue>) -> CanonicalKeyValue {
    match value {
        None => variant(0, []),
        Some(value) => variant(1, [value]),
    }
}

fn entity_span_key<T>(
    span: &EntitySpan<T>,
    value_key: impl Fn(&T) -> CanonicalKeyValue,
) -> CanonicalKeyValue {
    let (position, values) = match span {
        EntitySpan::Unchanged(value) => (SpanTagPosition::UNCHANGED, vec![value_key(value)]),
        EntitySpan::Added(value) => (SpanTagPosition::ADDED, vec![value_key(value)]),
        EntitySpan::Removed(value) => (SpanTagPosition::REMOVED, vec![value_key(value)]),
        EntitySpan::Modified { lhs, rhs } => (
            SpanTagPosition::MODIFIED,
            vec![value_key(lhs), value_key(rhs)],
        ),
    };
    CanonicalKeyValue::Span(SpanKey { position, values })
}

fn normalized_entity_span_key<T>(
    span: &EntitySpan<T>,
    value_key: impl Fn(&T) -> Result<CanonicalKeyValue, Contradiction>,
) -> Result<CanonicalKeyValue, Contradiction> {
    let (position, values) = match span {
        EntitySpan::Unchanged(value) => (SpanTagPosition::UNCHANGED, vec![value_key(value)?]),
        EntitySpan::Added(value) => (SpanTagPosition::ADDED, vec![value_key(value)?]),
        EntitySpan::Removed(value) => (SpanTagPosition::REMOVED, vec![value_key(value)?]),
        EntitySpan::Modified { lhs, rhs } => {
            let lhs = value_key(lhs)?;
            let rhs = value_key(rhs)?;
            if lhs == rhs {
                (SpanTagPosition::UNCHANGED, vec![lhs])
            } else {
                (SpanTagPosition::MODIFIED, vec![lhs, rhs])
            }
        }
    };
    Ok(CanonicalKeyValue::Span(SpanKey { position, values }))
}

fn boolean_form_key(value: BooleanForm) -> CanonicalKeyValue {
    match value {
        BooleanForm::Undetermined => variant(0, []),
        BooleanForm::Lit(value) => variant(1, [CanonicalKeyValue::Bool(value)]),
    }
}

fn index_key(index: usize) -> CanonicalKeyValue {
    CanonicalKeyValue::Unsigned(index as u64)
}

fn index_sequence(indices: impl IntoIterator<Item = usize>) -> CanonicalKeyValue {
    sequence(indices.into_iter().map(index_key))
}

fn element_form_key(value: &ElementForm) -> CanonicalKeyValue {
    match value {
        ElementForm::Undetermined => variant(0, []),
        ElementForm::Lit(element) => variant(
            1,
            [CanonicalKeyValue::Unsigned(element.atomic_number().into())],
        ),
        ElementForm::LitSet(elements) => variant(
            2,
            [sequence(elements.iter().map(|element| {
                CanonicalKeyValue::Unsigned(element.atomic_number().into())
            }))],
        ),
        ElementForm::NotSet(elements) => variant(
            3,
            [sequence(elements.iter().map(|element| {
                CanonicalKeyValue::Unsigned(element.atomic_number().into())
            }))],
        ),
        ElementForm::Var(variable) => {
            let (name, restriction) = variable.as_ref();
            variant(
                4,
                [
                    CanonicalKeyValue::Text(name.clone()),
                    option(restriction.as_ref().map(|(operator, elements)| {
                        product([
                            mem_op_key(*operator),
                            sequence(elements.iter().map(|element| {
                                CanonicalKeyValue::Unsigned(element.atomic_number().into())
                            })),
                        ])
                    })),
                ],
            )
        }
    }
}

fn isotope_mass_form_key(value: &IsotopeMassForm) -> CanonicalKeyValue {
    match value {
        IsotopeMassForm::Undetermined => variant(0, []),
        IsotopeMassForm::Natural => variant(1, []),
        IsotopeMassForm::Lit(mass) => variant(2, [CanonicalKeyValue::Unsigned(u64::from(*mass))]),
        IsotopeMassForm::LitSet(masses) => variant(
            3,
            [sequence(masses.iter().map(|mass| {
                CanonicalKeyValue::Unsigned(u64::from(*mass))
            }))],
        ),
        IsotopeMassForm::Var(variable) => {
            let (name, masses) = variable.as_ref();
            variant(
                4,
                [
                    CanonicalKeyValue::Text(name.clone()),
                    option(masses.as_ref().map(|masses| {
                        sequence(
                            masses
                                .iter()
                                .map(|mass| CanonicalKeyValue::Unsigned(u64::from(*mass))),
                        )
                    })),
                ],
            )
        }
    }
}

fn num_form_key(value: &NumForm) -> CanonicalKeyValue {
    match value {
        NumForm::Undetermined => variant(0, []),
        NumForm::Lit(value) => variant(1, [CanonicalKeyValue::Signed(*value)]),
        NumForm::LitSet(values) => variant(
            2,
            [sequence(
                values.iter().map(|value| CanonicalKeyValue::Signed(*value)),
            )],
        ),
        NumForm::RangeFrom(value) => variant(3, [CanonicalKeyValue::Signed(*value)]),
        NumForm::RangeTo(value) => variant(4, [CanonicalKeyValue::Signed(*value)]),
        NumForm::ArithExpr(expression) => variant(5, [arith_expr_key(expression)]),
        NumForm::PredExpr(predicate) => variant(6, [pred_expr_key(predicate)]),
    }
}

fn arith_expr_key(expression: &ArithExpr) -> CanonicalKeyValue {
    match expression {
        ArithExpr::Lit(value) => variant(0, [CanonicalKeyValue::Signed(*value)]),
        ArithExpr::Var(name) => variant(1, [CanonicalKeyValue::Text(name.clone())]),
        ArithExpr::Neg(inner) => variant(2, [arith_expr_key(inner)]),
        ArithExpr::Sum(terms) => variant(3, [sequence(terms.iter().map(arith_expr_key))]),
        ArithExpr::Product(factors) => variant(4, [sequence(factors.iter().map(arith_expr_key))]),
        ArithExpr::Div(lhs, rhs) => variant(5, [arith_expr_key(lhs), arith_expr_key(rhs)]),
        ArithExpr::Rem(lhs, rhs) => variant(6, [arith_expr_key(lhs), arith_expr_key(rhs)]),
    }
}

fn pred_expr_key(predicate: &PredExpr) -> CanonicalKeyValue {
    match predicate {
        PredExpr::Rel(lhs, operator, rhs) => variant(
            0,
            [
                arith_expr_key(lhs),
                rel_op_key(*operator),
                arith_expr_key(rhs),
            ],
        ),
        PredExpr::Mem(expression, operator, values) => variant(
            1,
            [
                arith_expr_key(expression),
                mem_op_key(*operator),
                sequence(values.iter().map(|value| CanonicalKeyValue::Signed(*value))),
            ],
        ),
        PredExpr::Not(inner) => variant(2, [pred_expr_key(inner)]),
        PredExpr::And(terms) => variant(3, [sequence(terms.iter().map(pred_expr_key))]),
        PredExpr::Or(terms) => variant(4, [sequence(terms.iter().map(pred_expr_key))]),
    }
}

fn rel_op_key(operator: RelOp) -> CanonicalKeyValue {
    variant(
        match operator {
            RelOp::Le => 0,
            RelOp::Ge => 1,
            RelOp::Eq => 2,
            RelOp::Lt => 3,
            RelOp::Gt => 4,
            RelOp::Ne => 5,
        },
        [],
    )
}

fn mem_op_key(operator: MemOp) -> CanonicalKeyValue {
    variant(
        match operator {
            MemOp::In => 0,
            MemOp::NotIn => 1,
        },
        [],
    )
}

fn unpaired_electrons_form_key(value: &UnpairedElectronsForm) -> CanonicalKeyValue {
    product([
        num_form_key(&value.count),
        num_form_key(&value.multiplicity),
    ])
}

fn aromatic_valence_form_key(value: &AromaticValenceForm) -> CanonicalKeyValue {
    match value {
        AromaticValenceForm::Undetermined => variant(0, []),
        AromaticValenceForm::NotAromatic => variant(1, []),
        AromaticValenceForm::Aromatic(value) => variant(2, [num_form_key(value)]),
    }
}

fn multicenter_valence_form_key(value: &MulticenterValenceForm) -> CanonicalKeyValue {
    match value {
        MulticenterValenceForm::Undetermined => variant(0, []),
        MulticenterValenceForm::NotMulticenter => variant(1, []),
        MulticenterValenceForm::Multicenter(value) => variant(2, [num_form_key(value)]),
    }
}

fn stereo_site_form_key(value: &TetrahedralStereoForm) -> CanonicalKeyValue {
    match value {
        TetrahedralStereoForm::Undetermined => variant(0, []),
        TetrahedralStereoForm::NotStereo => variant(1, []),
        TetrahedralStereoForm::Stereo(coset) => variant(2, [stereo_coset_key(coset)]),
    }
}

fn cis_trans_stereo_form_key(value: &CisTransStereoForm) -> CanonicalKeyValue {
    match value {
        CisTransStereoForm::Undetermined => variant(0, []),
        CisTransStereoForm::NotStereo => variant(1, []),
        CisTransStereoForm::Stereo(coset) => variant(2, [stereo_coset_key(coset)]),
    }
}

fn ring_scope_key(scope: RingScope) -> CanonicalKeyValue {
    match scope {
        RingScope::All => variant(0, []),
        RingScope::Size(size) => variant(1, [CanonicalKeyValue::Unsigned(size.into())]),
    }
}

fn ring_membership_form_key(value: &RingMembershipForm) -> CanonicalKeyValue {
    product([ring_scope_key(value.scope), num_form_key(&value.count)])
}

fn electron_counts_form_key(value: &ElectronCountsForm) -> CanonicalKeyValue {
    match value {
        ElectronCountsForm::Undetermined => variant(0, []),
        ElectronCountsForm::Lit(values) => variant(
            1,
            [sequence(
                values.iter().map(|&value| CanonicalKeyValue::Signed(value)),
            )],
        ),
    }
}

fn noncovalent_bond_kind_form_key(value: &NoncovalentBondKindForm) -> CanonicalKeyValue {
    match value {
        NoncovalentBondKindForm::Undetermined => variant(0, []),
        NoncovalentBondKindForm::Lit(kind) => variant(1, [noncovalent_bond_kind_key(*kind)]),
    }
}

fn noncovalent_bond_kind_key(kind: NoncovalentBondKind) -> CanonicalKeyValue {
    variant(
        match kind {
            NoncovalentBondKind::HydrogenBond => 0,
            NoncovalentBondKind::HalogenBond => 1,
            NoncovalentBondKind::ChalcogenBond => 2,
            NoncovalentBondKind::Ionic => 3,
            NoncovalentBondKind::VanDerWaals => 4,
        },
        [],
    )
}

fn stereo_kind_key(kind: StereoKind) -> CanonicalKeyValue {
    variant(
        match kind {
            StereoKind::Tetrahedral => 0,
            StereoKind::CisTrans => 1,
            StereoKind::Axial => 2,
            StereoKind::SquarePlanar => 3,
            StereoKind::TrigonalBipyramidal => 4,
            StereoKind::Octahedral => 5,
        },
        [],
    )
}

fn stereo_ligand_kind_key(kind: StereoLigandKind) -> CanonicalKeyValue {
    variant(
        match kind {
            StereoLigandKind::Atom => 0,
            StereoLigandKind::ImplicitHydrogen => 1,
            StereoLigandKind::LonePair => 2,
        },
        [],
    )
}

fn permutation_key(permutation: Permutation) -> CanonicalKeyValue {
    product([
        CanonicalKeyValue::Unsigned(permutation.degree() as u64),
        sequence(
            (0..permutation.degree())
                .map(|position| CanonicalKeyValue::Unsigned(permutation.apply(position) as u64)),
        ),
    ])
}

fn stereo_term_key(term: &StereoTerm) -> CanonicalKeyValue {
    match term {
        StereoTerm::Var(variable) => {
            let (name, domain) = variable.as_ref();
            variant(
                0,
                [
                    CanonicalKeyValue::Text(name.clone()),
                    option(domain.as_ref().map(|values| {
                        sequence(
                            values
                                .iter()
                                .map(|&value| CanonicalKeyValue::Unsigned(value.into())),
                        )
                    })),
                ],
            )
        }
        StereoTerm::Lit(value) => variant(1, [CanonicalKeyValue::Unsigned((*value).into())]),
        StereoTerm::LitSet(values) => variant(
            2,
            [sequence(
                values
                    .iter()
                    .map(|&value| CanonicalKeyValue::Unsigned(value.into())),
            )],
        ),
        StereoTerm::Swap(inner) => variant(3, [stereo_term_key(inner)]),
        StereoTerm::Mirror(inner) => variant(4, [stereo_term_key(inner)]),
        StereoTerm::Apply(inner, permutation) => {
            variant(5, [stereo_term_key(inner), permutation_key(*permutation)])
        }
    }
}

fn stereo_coset_key(coset: &StereoCoset) -> CanonicalKeyValue {
    match coset {
        StereoCoset::Undetermined => variant(0, []),
        StereoCoset::Lit(value) => variant(1, [CanonicalKeyValue::Unsigned((*value).into())]),
        StereoCoset::LitSet(values) => variant(
            2,
            [sequence(
                values
                    .iter()
                    .map(|&value| CanonicalKeyValue::Unsigned(value.into())),
            )],
        ),
        StereoCoset::Term(term) => variant(3, [stereo_term_key(term)]),
    }
}

fn stereo_configuration_form_key(configuration: &StereoConfigurationForm) -> CanonicalKeyValue {
    match configuration {
        StereoConfigurationForm::Undetermined => variant(0, []),
        StereoConfigurationForm::Kinded(kind, coset) => {
            variant(1, [stereo_kind_key(*kind), stereo_coset_key(coset)])
        }
    }
}

fn orientation_key(orientation: Orientation) -> CanonicalKeyValue {
    variant(
        match orientation {
            Orientation::Proper => 0,
            Orientation::Improper => 1,
        },
        [],
    )
}

fn topicity_key(topicity: Topicity) -> CanonicalKeyValue {
    variant(
        match topicity {
            Topicity::Homotopic => 0,
            Topicity::Enantiotopic => 1,
            Topicity::Diastereotopic => 2,
        },
        [],
    )
}

fn stereogenicity_key(stereogenicity: Stereogenicity) -> CanonicalKeyValue {
    variant(
        match stereogenicity {
            Stereogenicity::Symmetric => 0,
            Stereogenicity::Prochiral => 1,
            Stereogenicity::Stereogenic => 2,
        },
        [],
    )
}

fn topicity_relation_form_key(value: &TopicityRelationForm) -> CanonicalKeyValue {
    match value {
        TopicityRelationForm::Undetermined => variant(0, []),
        TopicityRelationForm::Lit(value) => variant(1, [topicity_key(*value)]),
        TopicityRelationForm::LitSet(values) => {
            variant(2, [sequence(values.iter().copied().map(topicity_key))])
        }
        TopicityRelationForm::NotSet(values) => {
            variant(3, [sequence(values.iter().copied().map(topicity_key))])
        }
    }
}

fn stereogenicity_form_key(value: &StereogenicityForm) -> CanonicalKeyValue {
    match value {
        StereogenicityForm::Undetermined => variant(0, []),
        StereogenicityForm::Lit(value) => variant(1, [stereogenicity_key(*value)]),
        StereogenicityForm::LitSet(values) => variant(
            2,
            [sequence(values.iter().copied().map(stereogenicity_key))],
        ),
        StereogenicityForm::NotSet(values) => variant(
            3,
            [sequence(values.iter().copied().map(stereogenicity_key))],
        ),
    }
}

fn ligand_permutation_key(value: LigandPermutation) -> CanonicalKeyValue {
    permutation_key(value.0)
}

fn oriented_ligand_permutation_key(value: OrientedLigandPermutation) -> CanonicalKeyValue {
    product([
        ligand_permutation_key(value.permutation),
        orientation_key(value.orientation),
    ])
}

fn stereo_ligand_pair_key(value: StereoLigandPair) -> CanonicalKeyValue {
    product([
        index_key(value.first().index()),
        index_key(value.second().index()),
    ])
}

fn atom_constraint_form_key(value: &AtomConstraintForm) -> CanonicalKeyValue {
    match value {
        AtomConstraintForm::Valence(value) => variant(0, [num_form_key(value)]),
        AtomConstraintForm::DonatedPairs(value) => variant(1, [num_form_key(value)]),
        AtomConstraintForm::AcceptedPairs(value) => variant(2, [num_form_key(value)]),
        AtomConstraintForm::AromaticValence(value) => {
            variant(3, [aromatic_valence_form_key(value)])
        }
        AtomConstraintForm::MulticenterValence(value) => {
            variant(4, [multicenter_valence_form_key(value)])
        }
        AtomConstraintForm::TetrahedralStereo(value) => variant(5, [stereo_site_form_key(value)]),
        AtomConstraintForm::Degree(value) => variant(6, [num_form_key(value)]),
        AtomConstraintForm::TotalDegree(value) => variant(7, [num_form_key(value)]),
        AtomConstraintForm::TotalValence(value) => variant(8, [num_form_key(value)]),
        AtomConstraintForm::RingDegree(value) => variant(9, [num_form_key(value)]),
        AtomConstraintForm::RingValence(value) => variant(10, [num_form_key(value)]),
        AtomConstraintForm::TotalHydrogens(value) => variant(11, [num_form_key(value)]),
        AtomConstraintForm::RingMembership(value) => variant(12, [ring_membership_form_key(value)]),
    }
}

fn bond_constraint_form_key(value: &BondConstraintForm) -> CanonicalKeyValue {
    match value {
        BondConstraintForm::Aromatic(value) => variant(0, [boolean_form_key(*value)]),
        BondConstraintForm::CisTransStereo(value) => variant(1, [cis_trans_stereo_form_key(value)]),
        BondConstraintForm::RingMembership(value) => variant(2, [ring_membership_form_key(value)]),
    }
}

fn dative_bond_constraint_form_key(value: &DativeBondConstraintForm) -> CanonicalKeyValue {
    match value {
        DativeBondConstraintForm::Aromatic(value) => variant(0, [boolean_form_key(*value)]),
        DativeBondConstraintForm::RingMembership(value) => {
            variant(1, [ring_membership_form_key(value)])
        }
    }
}

fn aromatic_system_constraint_form_key(value: &AromaticSystemConstraintForm) -> CanonicalKeyValue {
    match value {
        AromaticSystemConstraintForm::ElectronCount(value) => variant(0, [num_form_key(value)]),
    }
}

fn multicenter_bond_constraint_form_key(
    value: &MulticenterBondConstraintForm,
) -> CanonicalKeyValue {
    match value {
        MulticenterBondConstraintForm::ElectronCount(value) => variant(0, [num_form_key(value)]),
    }
}

fn noncovalent_bond_constraint_form_key(
    value: &NoncovalentBondConstraintForm,
) -> CanonicalKeyValue {
    match value {
        NoncovalentBondConstraintForm::Intramolecular(value) => {
            variant(0, [boolean_form_key(*value)])
        }
    }
}

fn stereo_atom_constraint_form_key(value: &StereoAtomConstraintForm) -> CanonicalKeyValue {
    match value {
        StereoAtomConstraintForm::LigandSymmetry(value) => variant(
            0,
            [
                oriented_ligand_permutation_key(value.permutation),
                boolean_form_key(value.invariant),
            ],
        ),
        StereoAtomConstraintForm::Fluxionality(value) => variant(
            1,
            [
                ligand_permutation_key(value.permutation),
                boolean_form_key(value.active),
            ],
        ),
        StereoAtomConstraintForm::Topicity(value) => variant(
            2,
            [
                stereo_ligand_pair_key(value.pair),
                topicity_relation_form_key(&value.relation),
            ],
        ),
        StereoAtomConstraintForm::Stereogenicity(value) => {
            variant(3, [stereogenicity_form_key(value)])
        }
    }
}

fn stereo_bond_constraint_form_key(value: &StereoBondConstraintForm) -> CanonicalKeyValue {
    match value {
        StereoBondConstraintForm::LigandSymmetry(value) => variant(
            0,
            [
                oriented_ligand_permutation_key(value.permutation),
                boolean_form_key(value.invariant),
            ],
        ),
        StereoBondConstraintForm::Fluxionality(value) => variant(
            1,
            [
                ligand_permutation_key(value.permutation),
                boolean_form_key(value.active),
            ],
        ),
        StereoBondConstraintForm::Topicity(value) => variant(
            2,
            [
                stereo_ligand_pair_key(value.pair),
                topicity_relation_form_key(&value.relation),
            ],
        ),
        StereoBondConstraintForm::Stereogenicity(value) => {
            variant(3, [stereogenicity_form_key(value)])
        }
    }
}

fn relational_constraint_key(value: &RelationalConstraint) -> CanonicalKeyValue {
    let atom_predicate = |value: &AtomConstraintForm| atom_constraint_form_key(value);
    let position = |entity: EntityBlockPosition, slot: u64| {
        let position = RelationalConstraintPosition::new(entity, slot as u16);
        product([
            CanonicalKeyValue::Unsigned(position.entity.domain.0.into()),
            CanonicalKeyValue::Unsigned(position.entity.slot.into()),
            CanonicalKeyValue::Unsigned(position.slot.into()),
        ])
    };
    let atom_ids = |ids: &[AtomId]| index_sequence(ids.iter().map(|id| id.index()));

    let (entity, slot, fields) = match value {
        RelationalConstraint::DativeBondDonors { bond, atoms } => (
            EntityBlockPosition::DATIVE_BOND,
            0,
            vec![index_key(bond.index()), atom_ids(atoms)],
        ),
        RelationalConstraint::DativeBondDonor { bond, atom } => (
            EntityBlockPosition::DATIVE_BOND,
            1,
            vec![index_key(bond.index()), index_key(atom.index())],
        ),
        RelationalConstraint::DativeBondContainsAllDonors { bond, atoms } => (
            EntityBlockPosition::DATIVE_BOND,
            2,
            vec![index_key(bond.index()), atom_ids(atoms)],
        ),
        RelationalConstraint::DativeBondAllDonors { bond, predicate } => (
            EntityBlockPosition::DATIVE_BOND,
            3,
            vec![index_key(bond.index()), atom_predicate(predicate)],
        ),
        RelationalConstraint::DativeBondAnyDonor { bond, predicate } => (
            EntityBlockPosition::DATIVE_BOND,
            4,
            vec![index_key(bond.index()), atom_predicate(predicate)],
        ),
        RelationalConstraint::DativeBondAcceptor { bond, atom } => (
            EntityBlockPosition::DATIVE_BOND,
            5,
            vec![index_key(bond.index()), index_key(atom.index())],
        ),
        RelationalConstraint::DativeBondAcceptorSatisfies { bond, predicate } => (
            EntityBlockPosition::DATIVE_BOND,
            6,
            vec![index_key(bond.index()), atom_predicate(predicate)],
        ),
        RelationalConstraint::DativeBondParallels { dative, parallel } => (
            EntityBlockPosition::DATIVE_BOND,
            7,
            vec![index_key(dative.index()), index_key(parallel.index())],
        ),
        RelationalConstraint::AromaticSystemAtoms { system, atoms } => (
            EntityBlockPosition::AROMATIC_SYSTEM,
            0,
            vec![index_key(system.index()), atom_ids(atoms)],
        ),
        RelationalConstraint::AromaticSystemContains { system, atom } => (
            EntityBlockPosition::AROMATIC_SYSTEM,
            1,
            vec![index_key(system.index()), index_key(atom.index())],
        ),
        RelationalConstraint::AromaticSystemContainsAll { system, atoms } => (
            EntityBlockPosition::AROMATIC_SYSTEM,
            2,
            vec![index_key(system.index()), atom_ids(atoms)],
        ),
        RelationalConstraint::AromaticSystemAllAtoms { system, predicate } => (
            EntityBlockPosition::AROMATIC_SYSTEM,
            3,
            vec![index_key(system.index()), atom_predicate(predicate)],
        ),
        RelationalConstraint::AromaticSystemAnyAtom { system, predicate } => (
            EntityBlockPosition::AROMATIC_SYSTEM,
            4,
            vec![index_key(system.index()), atom_predicate(predicate)],
        ),
        RelationalConstraint::MulticenterBondAtoms { bond, atoms } => (
            EntityBlockPosition::MULTICENTER_BOND,
            0,
            vec![index_key(bond.index()), atom_ids(atoms)],
        ),
        RelationalConstraint::MulticenterBondContains { bond, atom } => (
            EntityBlockPosition::MULTICENTER_BOND,
            1,
            vec![index_key(bond.index()), index_key(atom.index())],
        ),
        RelationalConstraint::MulticenterBondContainsAll { bond, atoms } => (
            EntityBlockPosition::MULTICENTER_BOND,
            2,
            vec![index_key(bond.index()), atom_ids(atoms)],
        ),
        RelationalConstraint::MulticenterBondAllAtoms { bond, predicate } => (
            EntityBlockPosition::MULTICENTER_BOND,
            3,
            vec![index_key(bond.index()), atom_predicate(predicate)],
        ),
        RelationalConstraint::MulticenterBondAnyAtom { bond, predicate } => (
            EntityBlockPosition::MULTICENTER_BOND,
            4,
            vec![index_key(bond.index()), atom_predicate(predicate)],
        ),
        RelationalConstraint::NoncovalentBondEnds { bond, atoms } => (
            EntityBlockPosition::NONCOVALENT_BOND,
            0,
            vec![
                index_key(bond.index()),
                index_sequence(atoms.iter().map(|atom| atom.index())),
            ],
        ),
        RelationalConstraint::NoncovalentBondContains { bond, atom } => (
            EntityBlockPosition::NONCOVALENT_BOND,
            1,
            vec![index_key(bond.index()), index_key(atom.index())],
        ),
        RelationalConstraint::NoncovalentBondEndsSatisfy { bond, predicates } => (
            EntityBlockPosition::NONCOVALENT_BOND,
            2,
            vec![
                index_key(bond.index()),
                sequence(predicates.iter().map(|predicate| atom_predicate(predicate))),
            ],
        ),
        RelationalConstraint::StereoAtomSite { stereo_atom, atom } => (
            EntityBlockPosition::STEREO_ATOM,
            0,
            vec![index_key(stereo_atom.index()), index_key(atom.index())],
        ),
        RelationalConstraint::StereoAtomContains { stereo_atom, atom } => (
            EntityBlockPosition::STEREO_ATOM,
            1,
            vec![index_key(stereo_atom.index()), index_key(atom.index())],
        ),
        RelationalConstraint::StereoAtomLigands { stereo_atom, atoms } => (
            EntityBlockPosition::STEREO_ATOM,
            2,
            vec![index_key(stereo_atom.index()), atom_ids(atoms)],
        ),
        RelationalConstraint::StereoAtomAllLigands {
            stereo_atom,
            predicate,
        } => (
            EntityBlockPosition::STEREO_ATOM,
            3,
            vec![index_key(stereo_atom.index()), atom_predicate(predicate)],
        ),
        RelationalConstraint::StereoAtomAnyLigand {
            stereo_atom,
            predicate,
        } => (
            EntityBlockPosition::STEREO_ATOM,
            4,
            vec![index_key(stereo_atom.index()), atom_predicate(predicate)],
        ),
        RelationalConstraint::StereoBondSite { stereo_bond, bond } => (
            EntityBlockPosition::STEREO_BOND,
            0,
            vec![index_key(stereo_bond.index()), index_key(bond.index())],
        ),
        RelationalConstraint::StereoBondContains { stereo_bond, atom } => (
            EntityBlockPosition::STEREO_BOND,
            1,
            vec![index_key(stereo_bond.index()), index_key(atom.index())],
        ),
        RelationalConstraint::StereoBondLigands { stereo_bond, atoms } => (
            EntityBlockPosition::STEREO_BOND,
            2,
            vec![index_key(stereo_bond.index()), atom_ids(atoms)],
        ),
        RelationalConstraint::StereoBondAllLigands {
            stereo_bond,
            predicate,
        } => (
            EntityBlockPosition::STEREO_BOND,
            3,
            vec![index_key(stereo_bond.index()), atom_predicate(predicate)],
        ),
        RelationalConstraint::StereoBondAnyLigand {
            stereo_bond,
            predicate,
        } => (
            EntityBlockPosition::STEREO_BOND,
            4,
            vec![index_key(stereo_bond.index()), atom_predicate(predicate)],
        ),
    };

    product([position(entity, slot), sequence(fields)])
}

fn molecule_constraint_key(value: &MoleculeConstraint) -> CanonicalKeyValue {
    let atom_subset = |atoms: &Option<Vec<AtomId>>| {
        option(
            atoms
                .as_ref()
                .map(|atoms| index_sequence(atoms.iter().map(|id| id.index()))),
        )
    };
    let bond_subset = |bonds: &Option<Vec<BondId>>| {
        option(
            bonds
                .as_ref()
                .map(|bonds| index_sequence(bonds.iter().map(|id| id.index()))),
        )
    };
    match value {
        MoleculeConstraint::ChargeSum { atoms, sum } => {
            variant(0, [atom_subset(atoms), num_form_key(sum)])
        }
        MoleculeConstraint::UnpairedElectronCoupling {
            atoms,
            unpaired_electrons,
        } => variant(
            1,
            [
                atom_subset(atoms),
                unpaired_electrons_form_key(unpaired_electrons),
            ],
        ),
        MoleculeConstraint::BondOrderSum { bonds, sum } => {
            variant(2, [bond_subset(bonds), num_form_key(sum)])
        }
        MoleculeConstraint::Connected { atoms } => variant(3, [atom_subset(atoms)]),
    }
}

fn constraint_key(value: &Constraint) -> CanonicalKeyValue {
    fn set_key(constraints: &[Constraint]) -> CanonicalKeyValue {
        let mut keys = constraints.iter().map(constraint_key).collect::<Vec<_>>();
        keys.sort();
        keys.dedup();
        sequence(keys)
    }

    match value {
        Constraint::Atom(id, constraint) => variant(
            0,
            [
                product([
                    CanonicalKeyValue::Unsigned(EntityBlockPosition::ATOM.domain.0.into()),
                    CanonicalKeyValue::Unsigned(EntityBlockPosition::ATOM.slot.into()),
                ]),
                index_key(id.index()),
                atom_constraint_form_key(constraint),
            ],
        ),
        Constraint::Bond(id, constraint) => variant(
            0,
            [
                product([
                    CanonicalKeyValue::Unsigned(EntityBlockPosition::BOND.domain.0.into()),
                    CanonicalKeyValue::Unsigned(EntityBlockPosition::BOND.slot.into()),
                ]),
                index_key(id.index()),
                bond_constraint_form_key(constraint),
            ],
        ),
        Constraint::DativeBond(id, constraint) => variant(
            0,
            [
                product([
                    CanonicalKeyValue::Unsigned(EntityBlockPosition::DATIVE_BOND.domain.0.into()),
                    CanonicalKeyValue::Unsigned(EntityBlockPosition::DATIVE_BOND.slot.into()),
                ]),
                index_key(id.index()),
                dative_bond_constraint_form_key(constraint),
            ],
        ),
        Constraint::AromaticSystem(id, constraint) => variant(
            0,
            [
                product([
                    CanonicalKeyValue::Unsigned(
                        EntityBlockPosition::AROMATIC_SYSTEM.domain.0.into(),
                    ),
                    CanonicalKeyValue::Unsigned(EntityBlockPosition::AROMATIC_SYSTEM.slot.into()),
                ]),
                index_key(id.index()),
                aromatic_system_constraint_form_key(constraint),
            ],
        ),
        Constraint::MulticenterBond(id, constraint) => variant(
            0,
            [
                product([
                    CanonicalKeyValue::Unsigned(
                        EntityBlockPosition::MULTICENTER_BOND.domain.0.into(),
                    ),
                    CanonicalKeyValue::Unsigned(EntityBlockPosition::MULTICENTER_BOND.slot.into()),
                ]),
                index_key(id.index()),
                multicenter_bond_constraint_form_key(constraint),
            ],
        ),
        Constraint::NoncovalentBond(id, constraint) => variant(
            0,
            [
                product([
                    CanonicalKeyValue::Unsigned(
                        EntityBlockPosition::NONCOVALENT_BOND.domain.0.into(),
                    ),
                    CanonicalKeyValue::Unsigned(EntityBlockPosition::NONCOVALENT_BOND.slot.into()),
                ]),
                index_key(id.index()),
                noncovalent_bond_constraint_form_key(constraint),
            ],
        ),
        Constraint::StereoAtom(id, kind, constraint) => variant(
            0,
            [
                product([
                    CanonicalKeyValue::Unsigned(EntityBlockPosition::STEREO_ATOM.domain.0.into()),
                    CanonicalKeyValue::Unsigned(EntityBlockPosition::STEREO_ATOM.slot.into()),
                ]),
                index_key(id.index()),
                stereo_kind_key(*kind),
                stereo_atom_constraint_form_key(constraint),
            ],
        ),
        Constraint::StereoBond(id, kind, constraint) => variant(
            0,
            [
                product([
                    CanonicalKeyValue::Unsigned(EntityBlockPosition::STEREO_BOND.domain.0.into()),
                    CanonicalKeyValue::Unsigned(EntityBlockPosition::STEREO_BOND.slot.into()),
                ]),
                index_key(id.index()),
                stereo_kind_key(*kind),
                stereo_bond_constraint_form_key(constraint),
            ],
        ),
        Constraint::Relational(constraint) => variant(1, [relational_constraint_key(constraint)]),
        Constraint::Molecule(constraint) => variant(2, [molecule_constraint_key(constraint)]),
        Constraint::And(constraints) => variant(3, [set_key(constraints)]),
        Constraint::Or(constraints) => variant(4, [set_key(constraints)]),
        Constraint::Not(constraint) => variant(5, [constraint_key(constraint)]),
    }
}

/// Construct the ordered constraint blocks used by complete canonicalization.
fn constraint_blocks(molecule: &Molecule) -> Vec<ConstraintBlockKey> {
    let mut blocks = Vec::new();
    macro_rules! inline_block {
        ($position:expr, $entities:expr, $key:expr) => {{
            let rows = $entities
                .iter()
                .flat_map(|entity| {
                    entity.attributes.constraints.iter().map(move |constraint| {
                        product([index_key(entity.id.index()), $key(constraint)])
                    })
                })
                .collect::<Vec<_>>();
            if !rows.is_empty() {
                blocks.push(PositionedKey {
                    position: $position,
                    value: sequence(rows),
                });
            }
        }};
    }

    inline_block!(
        ConstraintBlockPosition::ATOM,
        molecule.atoms(),
        atom_constraint_form_key
    );
    inline_block!(
        ConstraintBlockPosition::BOND,
        molecule.bonds(),
        bond_constraint_form_key
    );
    inline_block!(
        ConstraintBlockPosition::DATIVE_BOND,
        molecule.dative_bonds(),
        dative_bond_constraint_form_key
    );
    inline_block!(
        ConstraintBlockPosition::AROMATIC_SYSTEM,
        molecule.aromatic_systems(),
        aromatic_system_constraint_form_key
    );
    inline_block!(
        ConstraintBlockPosition::MULTICENTER_BOND,
        molecule.multicenter_bonds(),
        multicenter_bond_constraint_form_key
    );
    inline_block!(
        ConstraintBlockPosition::NONCOVALENT_BOND,
        molecule.noncovalent_bonds(),
        noncovalent_bond_constraint_form_key
    );
    inline_block!(
        ConstraintBlockPosition::STEREO_ATOM,
        molecule.stereo_atoms(),
        stereo_atom_constraint_form_key
    );
    inline_block!(
        ConstraintBlockPosition::STEREO_BOND,
        molecule.stereo_bonds(),
        stereo_bond_constraint_form_key
    );

    if !molecule.constraints().is_empty() {
        let mut constraints = molecule
            .constraints()
            .iter()
            .map(constraint_key)
            .collect::<Vec<_>>();
        constraints.sort();
        constraints.dedup();
        blocks.push(PositionedKey {
            position: ConstraintBlockPosition::MOLECULE,
            value: sequence(constraints),
        });
    }
    blocks
}

fn stereo_ligand_key(atom: u32, kind: StereoLigandKind) -> CanonicalKeyValue {
    product([
        CanonicalKeyValue::Unsigned(atom.into()),
        stereo_ligand_kind_key(kind),
    ])
}

fn incidence_key(incidence: &Incidence) -> Result<CanonicalKeyValue, Contradiction> {
    Ok(match incidence {
        Incidence::BondEndpoint => variant(0, []),
        Incidence::DativeDonor => variant(1, []),
        Incidence::DativeAcceptor => variant(2, []),
        Incidence::AromaticParticipant(value) => {
            variant(3, [num_form_key(value.normalized()?.as_ref())])
        }
        Incidence::AromaticParticipantSpan(value) => variant(
            3,
            [normalized_entity_span_key(value, |value| {
                Ok(num_form_key(value.normalized()?.as_ref()))
            })?],
        ),
        Incidence::MulticenterParticipant(value) => {
            variant(4, [num_form_key(value.normalized()?.as_ref())])
        }
        Incidence::MulticenterParticipantSpan(value) => variant(
            4,
            [normalized_entity_span_key(value, |value| {
                Ok(num_form_key(value.normalized()?.as_ref()))
            })?],
        ),
        Incidence::NoncovalentEndpoint => variant(5, []),
        Incidence::StereoSite => variant(6, []),
        Incidence::StereoLigand(kind) => variant(7, [stereo_ligand_kind_key(*kind)]),
    })
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum InitialColorKey {
    Entity {
        position: EntityBlockPosition,
        value: CanonicalKeyValue,
    },
    Incidence(CanonicalKeyValue),
}

impl Ord for InitialColorKey {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (
                Self::Entity {
                    position: lhs_position,
                    value: lhs_value,
                },
                Self::Entity {
                    position: rhs_position,
                    value: rhs_value,
                },
            ) => lhs_position
                .cmp(rhs_position)
                .then_with(|| lhs_value.cmp(rhs_value)),
            (Self::Entity { .. }, Self::Incidence(_)) => Ordering::Less,
            (Self::Incidence(_), Self::Entity { .. }) => Ordering::Greater,
            (Self::Incidence(lhs), Self::Incidence(rhs)) => lhs.cmp(rhs),
        }
    }
}

impl PartialOrd for InitialColorKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct InitialColors {
    entities: Vec<u32>,
    incidences: Vec<u32>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum AutomorphismColor {
    Entity(u32),
    Incidence(u32),
}

/// Colored simple-graph encoding of an incidence graph for automorphism operations.
///
/// Source entity nodes keep their ids. Role- or value-bearing incidence edges are represented by
/// colored occurrence nodes; single-role endpoints remain direct unless parallel incidences require
/// subdivision.
#[derive(Clone, Debug)]
struct AutomorphismAdapter {
    graph: Graph,
    colors: Vec<AutomorphismColor>,
    node_sources: Vec<SubdivisionNodeSource>,
    incidence_nodes: Vec<Option<NodeId>>,
    source_node_count: usize,
}

/// Automorphism data over both the source-entity and complete adapter-node domains.
#[derive(Clone, Debug, PartialEq, Eq)]
struct AutomorphismAdapterOutput {
    source_orbits: Vec<NodeId>,
    /// Source entity nodes in backend canonical-label order, used only to order search branches.
    source_canonical_labels: Vec<NodeId>,
    source_generators: Vec<Vec<NodeId>>,
}

impl AutomorphismAdapter {
    /// Construct the exact graph adapter used by canonicalization.
    fn new(incidence_graph: &IncidenceGraph, initial_colors: &InitialColors) -> Self {
        let source = incidence_graph.graph();
        debug_assert_eq!(initial_colors.entities.len(), source.node_count());
        debug_assert_eq!(initial_colors.incidences.len(), source.edge_count());

        let mut colors = initial_colors
            .entities
            .iter()
            .copied()
            .map(AutomorphismColor::Entity)
            .collect::<Vec<_>>();
        let mut node_sources = source
            .node_ids()
            .map(SubdivisionNodeSource::Node)
            .collect::<Vec<_>>();
        let mut incidence_nodes = vec![None; source.edge_count()];
        let mut edges = Vec::new();

        let direct_pair_counts = source
            .edge_ids()
            .filter(|&edge| {
                matches!(
                    incidence_graph.incidence(edge),
                    Incidence::BondEndpoint | Incidence::NoncovalentEndpoint
                )
            })
            .fold(BTreeMap::<[NodeId; 2], usize>::new(), |mut counts, edge| {
                *counts.entry(source.edge_endpoints(edge)).or_default() += 1;
                counts
            });

        let mut push_edge = |endpoints: [NodeId; 2]| {
            edges.push([endpoints[0].0, endpoints[1].0]);
        };
        for edge in source.edge_ids() {
            let endpoints = source.edge_endpoints(edge);
            let direct = matches!(
                incidence_graph.incidence(edge),
                Incidence::BondEndpoint | Incidence::NoncovalentEndpoint
            ) && direct_pair_counts[&endpoints] == 1;
            if direct {
                push_edge(endpoints);
                continue;
            }

            let occurrence = NodeId(node_sources.len() as u32);
            node_sources.push(SubdivisionNodeSource::Edge(edge));
            incidence_nodes[edge.index()] = Some(occurrence);
            colors.push(AutomorphismColor::Incidence(
                initial_colors.incidences[edge.index()],
            ));
            push_edge([endpoints[0], occurrence]);
            push_edge([occurrence, endpoints[1]]);
        }
        let graph = Graph::new(node_sources.len(), &edges);
        debug_assert!(graph.is_simple());

        Self {
            graph,
            colors,
            node_sources,
            incidence_nodes,
            source_node_count: source.node_count(),
        }
    }

    fn graph(&self) -> &Graph {
        &self.graph
    }

    fn color(&self, node: NodeId) -> AutomorphismColor {
        self.colors[node.index()]
    }

    fn node_source(&self, node: NodeId) -> SubdivisionNodeSource {
        self.node_sources[node.index()]
    }

    fn node_of(&self, source: SubdivisionNodeSource) -> Option<NodeId> {
        match source {
            SubdivisionNodeSource::Node(node) => Some(node),
            SubdivisionNodeSource::Edge(edge) => self.incidence_nodes[edge.index()],
        }
    }

    fn automorphisms_for_partition(
        &self,
        partition: &OrderedPartition,
        algorithm: AutomorphismAlgorithm,
    ) -> AutomorphismAdapterOutput {
        let cell_indices = partition.cell_indices(self.graph().node_count());
        // Search partitions may deliberately coarsen covariant occurrence data. Retaining the
        // exact adapter color here keeps orbit pruning restricted to true automorphisms.
        let output = self.graph().automorphisms(
            |node| (self.color(node), cell_indices[node.index()]),
            algorithm,
        );

        self.project_automorphisms(&output)
    }

    fn project_automorphisms(&self, output: &AutomorphismOutput) -> AutomorphismAdapterOutput {
        let source_node = |node| match self.node_source(node) {
            SubdivisionNodeSource::Node(source) => source,
            SubdivisionNodeSource::Edge(_) => {
                unreachable!("disjoint colors preserve the adapter node domain")
            }
        };
        let source_orbits = (0..self.source_node_count)
            .map(|index| source_node(output.orbit_of(NodeId(index as u32))))
            .collect();
        let source_canonical_labels = output
            .canonical_labels()
            .iter()
            .filter_map(|&node| match self.node_source(node) {
                SubdivisionNodeSource::Node(source) => Some(source),
                SubdivisionNodeSource::Edge(_) => None,
            })
            .collect();
        let source_generators = output
            .generators()
            .iter()
            .map(|generator| {
                generator[..self.source_node_count]
                    .iter()
                    .copied()
                    .map(source_node)
                    .collect()
            })
            .collect();
        AutomorphismAdapterOutput {
            source_orbits,
            source_canonical_labels,
            source_generators,
        }
    }
}

fn generator_preserves_stereo(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    generator: &[NodeId],
) -> bool {
    let entity_image =
        |entity| incidence_graph.entity(generator[incidence_graph.node_of(entity).index()]);
    let frame_action = |source: Vec<StereoLigand>, target: Vec<StereoLigand>| {
        let mapped = source
            .into_iter()
            .map(|ligand| {
                let Entity::Atom(atom_id) = entity_image(Entity::Atom(ligand.atom_id)) else {
                    return None;
                };
                Some(StereoLigand::new(atom_id, ligand.kind))
            })
            .collect::<Option<Vec<_>>>()?;
        Permutation::between(&mapped, &target)
    };
    let configuration_preserved = |source: &StereoConfigurationForm,
                                   target: &StereoConfigurationForm,
                                   action: Permutation| {
        source
            .clone()
            .reframe_by(&action)
            .is_some_and(|mapped| mapped.normalized_eq(target))
    };

    for source in molecule.stereo_atoms().iter() {
        let Entity::StereoAtom(target_id) = entity_image(Entity::StereoAtom(source.id)) else {
            return false;
        };
        let target = molecule.stereo_atom(target_id);
        if entity_image(Entity::Atom(source.site_id())) != Entity::Atom(target.site_id()) {
            return false;
        }
        let Some(action) = frame_action(source.ligand_frame(), target.ligand_frame()) else {
            return false;
        };
        if !configuration_preserved(
            &source.attributes.configuration,
            &target.attributes.configuration,
            action,
        ) {
            return false;
        }
    }

    for source in molecule.stereo_bonds().iter() {
        let Entity::StereoBond(target_id) = entity_image(Entity::StereoBond(source.id)) else {
            return false;
        };
        let target = molecule.stereo_bond(target_id);
        if entity_image(Entity::Bond(source.site_id())) != Entity::Bond(target.site_id()) {
            return false;
        }
        let Some(action) = frame_action(source.ligand_frame(), target.ligand_frame()) else {
            return false;
        };
        if !configuration_preserved(
            &source.attributes.configuration,
            &target.attributes.configuration,
            action,
        ) {
            return false;
        }
    }

    true
}

fn retain_stereo_preserving_generators(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    generators: &mut Vec<Vec<NodeId>>,
) {
    if molecule.stereo_atoms().count() == 0 && molecule.stereo_bonds().count() == 0 {
        return;
    }
    generators.retain(|generator| generator_preserves_stereo(molecule, incidence_graph, generator));
}

fn source_orbits_from_generators(
    source_node_count: usize,
    generators: &[Vec<NodeId>],
) -> Vec<NodeId> {
    let mut sets = UnionFind::new(source_node_count);
    for generator in generators {
        debug_assert_eq!(generator.len(), source_node_count);
        for (source, image) in generator.iter().enumerate() {
            sets.union(source, image.index());
        }
    }

    let roots = (0..source_node_count)
        .map(|source| sets.find(source))
        .collect::<Vec<_>>();
    let mut representatives = vec![usize::MAX; source_node_count];
    for (source, root) in roots.iter().copied().enumerate() {
        representatives[root] = representatives[root].min(source);
    }
    roots
        .into_iter()
        .map(|root| NodeId(representatives[root] as u32))
        .collect()
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct OrderedPartition {
    cells: Vec<Vec<NodeId>>,
}

impl OrderedPartition {
    fn from_descriptors<C: Clone + Ord>(descriptors: &[C]) -> Self {
        let mut cells = BTreeMap::<C, Vec<NodeId>>::new();
        for (index, descriptor) in descriptors.iter().enumerate() {
            cells
                .entry(descriptor.clone())
                .or_default()
                .push(NodeId(index as u32));
        }

        Self {
            cells: cells.into_values().collect(),
        }
    }

    fn refine(mut self, graph: &Graph) -> Self {
        loop {
            let cell_indices = self.cell_indices(graph.node_count());
            let cell_count = self.cells.len();
            let mut refined = Vec::with_capacity(self.cells.len());
            let mut changed = false;

            for cell in self.cells {
                let mut splits = BTreeMap::<Vec<u32>, Vec<NodeId>>::new();
                for node in cell {
                    let mut signature = vec![0; cell_count];
                    for neighbor in graph.neighbors(node) {
                        signature[cell_indices[neighbor.node.index()] as usize] += 1;
                    }
                    splits.entry(signature).or_default().push(node);
                }
                changed |= splits.len() > 1;
                // The minimum sorted-incidence key places the greater exact signature first.
                refined.extend(splits.into_values().rev());
            }

            self.cells = refined;
            if !changed {
                return self;
            }
        }
    }

    fn individualize(&self, cell_index: usize, node: NodeId) -> Self {
        let mut cells = Vec::with_capacity(self.cells.len() + 1);
        for (index, cell) in self.cells.iter().enumerate() {
            if index != cell_index {
                cells.push(cell.clone());
                continue;
            }

            cells.push(vec![node]);
            let remainder = cell
                .iter()
                .copied()
                .filter(|&candidate| candidate != node)
                .collect::<Vec<_>>();
            if !remainder.is_empty() {
                cells.push(remainder);
            }
        }

        Self { cells }
    }

    fn first_non_singleton_entity_cell(&self, entity_count: usize) -> Option<usize> {
        self.cells.iter().position(|cell| {
            cell.len() > 1 && cell.first().is_some_and(|node| node.index() < entity_count)
        })
    }

    fn entity_order(&self, entity_count: usize) -> Vec<NodeId> {
        self.cells
            .iter()
            .flatten()
            .copied()
            .filter(|node| node.index() < entity_count)
            .collect()
    }

    fn fixed_entity_prefix(&self, entity_count: usize) -> Vec<NodeId> {
        let mut prefix = Vec::new();
        for cell in &self.cells {
            let Some(&node) = cell.first() else {
                continue;
            };
            if node.index() >= entity_count {
                continue;
            }
            if cell.len() != 1 {
                break;
            }
            prefix.push(node);
        }
        prefix
    }

    fn cell_indices(&self, node_count: usize) -> Vec<u32> {
        let mut indices = vec![0; node_count];
        for (cell_index, cell) in self.cells.iter().enumerate() {
            for node in cell {
                indices[node.index()] = cell_index as u32;
            }
        }
        indices
    }
}

type BranchOrdering = fn(
    &AutomorphismAdapter,
    &OrderedPartition,
    AutomorphismAlgorithm,
    Option<&AutomorphismAdapterOutput>,
    &mut [NodeId],
) -> bool;

fn backend_canonical_branch_order(
    adapter: &AutomorphismAdapter,
    partition: &OrderedPartition,
    algorithm: AutomorphismAlgorithm,
    automorphisms: Option<&AutomorphismAdapterOutput>,
    candidates: &mut [NodeId],
) -> bool {
    let backend_called = automorphisms.is_none();
    let labels = automorphisms.map_or_else(
        || {
            adapter
                .automorphisms_for_partition(partition, algorithm)
                .source_canonical_labels
        },
        |output| output.source_canonical_labels.clone(),
    );
    let mut ranks = vec![0; adapter.source_node_count];
    for (rank, node) in labels.iter().enumerate() {
        ranks[node.index()] = rank;
    }
    candidates.sort_unstable_by_key(|node| ranks[node.index()]);
    backend_called
}

#[derive(Clone, Copy, Debug)]
struct CanonicalSearchOptions {
    automorphism_pruning: bool,
    prefix_pruning: bool,
    branch_order: BranchOrdering,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct CanonicalSearchStats {
    initial_residual_cell_sizes: Vec<usize>,
    refinement_calls: usize,
    branch_order_calls: usize,
    backend_calls: usize,
    visited_leaves: usize,
    leaf_comparisons: usize,
    prefix_pruned_branches: usize,
    orbit_pruned_branches: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct CanonicalCandidate<K> {
    key: K,
    entity_order: Vec<NodeId>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct CanonicalSearchResult<K> {
    candidate: CanonicalCandidate<K>,
    stats: CanonicalSearchStats,
}

/// Minimize a typed leaf candidate over an exact partition ordered by semantic descriptors.
///
/// Adapter colors are opaque equality labels and do not order the partition. Automorphism pruning
/// requires the leaf key to be invariant under adapter automorphisms.
/// Prefix pruning requires `prefix_worse` to reject only partitions whose every leaf is greater
/// than the current best key.
fn canonical_search<K, Descriptor, LeafCandidate, PrefixWorse>(
    adapter: &AutomorphismAdapter,
    partition_descriptors: &[Descriptor],
    algorithm: AutomorphismAlgorithm,
    options: CanonicalSearchOptions,
    leaf_candidate: &LeafCandidate,
    prefix_worse: &PrefixWorse,
) -> CanonicalSearchResult<K>
where
    K: Ord,
    Descriptor: Clone + Ord,
    LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
    PrefixWorse: Fn(&OrderedPartition, &CanonicalCandidate<K>) -> bool,
{
    canonical_search_impl(
        adapter,
        partition_descriptors,
        algorithm,
        options,
        leaf_candidate,
        prefix_worse,
        false,
        &|_| {},
    )
}

#[allow(clippy::too_many_arguments)]
fn canonical_search_with_generator_filter<
    K,
    Descriptor,
    LeafCandidate,
    PrefixWorse,
    FilterGenerators,
>(
    adapter: &AutomorphismAdapter,
    partition_descriptors: &[Descriptor],
    algorithm: AutomorphismAlgorithm,
    options: CanonicalSearchOptions,
    leaf_candidate: &LeafCandidate,
    prefix_worse: &PrefixWorse,
    filter_generators: &FilterGenerators,
) -> CanonicalSearchResult<K>
where
    K: Ord,
    Descriptor: Clone + Ord,
    LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
    PrefixWorse: Fn(&OrderedPartition, &CanonicalCandidate<K>) -> bool,
    FilterGenerators: Fn(&mut Vec<Vec<NodeId>>),
{
    canonical_search_impl(
        adapter,
        partition_descriptors,
        algorithm,
        options,
        leaf_candidate,
        prefix_worse,
        true,
        filter_generators,
    )
}

#[allow(clippy::too_many_arguments)]
fn canonical_search_impl<K, Descriptor, LeafCandidate, PrefixWorse, FilterGenerators>(
    adapter: &AutomorphismAdapter,
    partition_descriptors: &[Descriptor],
    algorithm: AutomorphismAlgorithm,
    options: CanonicalSearchOptions,
    leaf_candidate: &LeafCandidate,
    prefix_worse: &PrefixWorse,
    apply_generator_filter: bool,
    filter_generators: &FilterGenerators,
) -> CanonicalSearchResult<K>
where
    K: Ord,
    Descriptor: Clone + Ord,
    LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
    PrefixWorse: Fn(&OrderedPartition, &CanonicalCandidate<K>) -> bool,
    FilterGenerators: Fn(&mut Vec<Vec<NodeId>>),
{
    let initial = OrderedPartition::from_descriptors(partition_descriptors).refine(adapter.graph());
    let mut best = None;
    let mut stats = CanonicalSearchStats {
        initial_residual_cell_sizes: initial
            .cells
            .iter()
            .filter(|cell| {
                cell.len() > 1
                    && cell
                        .first()
                        .is_some_and(|node| node.index() < adapter.source_node_count)
            })
            .map(Vec::len)
            .collect(),
        refinement_calls: 1,
        ..Default::default()
    };

    search_partition(
        adapter,
        initial,
        algorithm,
        options,
        leaf_candidate,
        prefix_worse,
        apply_generator_filter,
        filter_generators,
        &mut best,
        &mut stats,
    );

    CanonicalSearchResult {
        candidate: best.expect("every finite partition has a discrete entity refinement"),
        stats,
    }
}

#[allow(clippy::too_many_arguments)]
fn search_partition<K, LeafCandidate, PrefixWorse, FilterGenerators>(
    adapter: &AutomorphismAdapter,
    partition: OrderedPartition,
    algorithm: AutomorphismAlgorithm,
    options: CanonicalSearchOptions,
    leaf_candidate: &LeafCandidate,
    prefix_worse: &PrefixWorse,
    apply_generator_filter: bool,
    filter_generators: &FilterGenerators,
    best: &mut Option<CanonicalCandidate<K>>,
    stats: &mut CanonicalSearchStats,
) where
    K: Ord,
    LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
    PrefixWorse: Fn(&OrderedPartition, &CanonicalCandidate<K>) -> bool,
    FilterGenerators: Fn(&mut Vec<Vec<NodeId>>),
{
    if options.prefix_pruning
        && best
            .as_ref()
            .is_some_and(|best| prefix_worse(&partition, best))
    {
        stats.prefix_pruned_branches += 1;
        return;
    }

    let Some(cell_index) = partition.first_non_singleton_entity_cell(adapter.source_node_count)
    else {
        stats.visited_leaves += 1;
        let entity_order = partition.entity_order(adapter.source_node_count);
        let candidate = leaf_candidate(&entity_order);
        let improves = if let Some(best) = best.as_ref() {
            stats.leaf_comparisons += 1;
            candidate.key < best.key
        } else {
            true
        };
        if improves {
            *best = Some(candidate);
        }
        return;
    };

    let mut candidates = partition.cells[cell_index].clone();
    let automorphisms = (options.automorphism_pruning || apply_generator_filter).then(|| {
        stats.backend_calls += 1;
        let mut output = adapter.automorphisms_for_partition(&partition, algorithm);
        if apply_generator_filter {
            filter_generators(&mut output.source_generators);
            output.source_orbits =
                source_orbits_from_generators(adapter.source_node_count, &output.source_generators);
        }
        output
    });

    if options.automorphism_pruning {
        let orbits = &automorphisms
            .as_ref()
            .expect("automorphisms requested for orbit pruning")
            .source_orbits;
        let mut representatives = BTreeMap::<NodeId, NodeId>::new();
        for candidate in candidates {
            representatives
                .entry(orbits[candidate.index()])
                .and_modify(|representative| *representative = (*representative).min(candidate))
                .or_insert(candidate);
        }
        stats.orbit_pruned_branches += partition.cells[cell_index].len() - representatives.len();
        candidates = representatives.into_values().collect();
    }

    stats.branch_order_calls += 1;
    let backend_called = (options.branch_order)(
        adapter,
        &partition,
        algorithm,
        automorphisms.as_ref(),
        &mut candidates,
    );
    stats.backend_calls += usize::from(backend_called);

    for candidate in candidates {
        stats.refinement_calls += 1;
        let child = partition
            .individualize(cell_index, candidate)
            .refine(adapter.graph());
        search_partition(
            adapter,
            child,
            algorithm,
            options,
            leaf_candidate,
            prefix_worse,
            apply_generator_filter,
            filter_generators,
            best,
            stats,
        );
    }
}

/// Construct the normalized entity and incidence keys used for initial colors.
fn initial_color_keys(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
) -> Result<(Vec<InitialColorKey>, Vec<InitialColorKey>), Contradiction> {
    let entity_keys = incidence_graph
        .graph()
        .node_ids()
        .map(|node| entity_color_key(molecule, incidence_graph.entity(node)))
        .collect::<Result<Vec<_>, _>>()?;
    let incidence_keys = incidence_graph
        .incidences()
        .map(|(_, incidence)| incidence_key(incidence).map(InitialColorKey::Incidence))
        .collect::<Result<Vec<_>, _>>()?;

    Ok((entity_keys, incidence_keys))
}

/// Construct topology-level partition descriptors for the exact adapter.
fn partition_descriptors(
    adapter: &AutomorphismAdapter,
    entity_keys: &[InitialColorKey],
    incidence_keys: &[InitialColorKey],
) -> Vec<InitialColorKey> {
    adapter
        .node_sources
        .iter()
        .map(|source| match *source {
            SubdivisionNodeSource::Node(node) => entity_keys[node.index()].clone(),
            SubdivisionNodeSource::Edge(edge) => incidence_keys[edge.index()].clone(),
        })
        .collect()
}

/// Search descriptors for the constitution key.
///
/// Dative roles and positional electron counts remain exact automorphism colors, but they do not
/// fix the search partition: their order is selected by the typed leaf key. Using them for both
/// purposes would exclude valid dense relabelings before the key can compare them.
/// Construct constitution-level partition descriptors for the exact adapter.
fn constitution_partition_descriptors(
    adapter: &AutomorphismAdapter,
    entity_keys: &[InitialColorKey],
    incidence_graph: &IncidenceGraph,
) -> Vec<InitialColorKey> {
    adapter
        .node_sources
        .iter()
        .map(|source| match *source {
            SubdivisionNodeSource::Node(node) => entity_keys[node.index()].clone(),
            SubdivisionNodeSource::Edge(edge) => {
                let value = match incidence_graph.incidence(edge) {
                    Incidence::DativeDonor | Incidence::DativeAcceptor => variant(1, []),
                    Incidence::AromaticParticipant(_) | Incidence::AromaticParticipantSpan(_) => {
                        variant(3, [])
                    }
                    Incidence::MulticenterParticipant(_)
                    | Incidence::MulticenterParticipantSpan(_) => variant(4, []),
                    incidence => incidence_key(incidence)
                        .expect("initial colors established incidence normalization"),
                };
                InitialColorKey::Incidence(value)
            }
        })
        .collect()
}

fn constitution_entity_classes(molecule: &Molecule) -> Result<Vec<u32>, Contradiction> {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Constitution);
    let (entity_keys, incidence_keys) = initial_color_keys(molecule, &incidence_graph)?;
    let colors = rank_initial_colors(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let descriptors = constitution_partition_descriptors(&adapter, &entity_keys, &incidence_graph);
    let partition = OrderedPartition::from_descriptors(&descriptors).refine(adapter.graph());
    let classes = partition.cell_indices(adapter.graph().node_count());
    Ok(classes[..adapter.source_node_count].to_vec())
}

fn stereo_frame_permutations(kind: StereoKind) -> impl Iterator<Item = Permutation> {
    let degree = kind.degree();
    let count = (1..=degree).product::<usize>().max(1);
    (0..count)
        .map(move |rank| Permutation::unrank(degree, rank))
        .filter(move |permutation| kind.class_key().space().allows(*permutation))
}

fn stereo_refinement_descriptor(
    site_class: u32,
    ligand_classes: &[(u32, StereoLigandKind)],
    configuration: &StereoConfigurationForm,
) -> Result<CanonicalKeyValue, Contradiction> {
    let descriptor = |ligands: Vec<(u32, StereoLigandKind)>,
                      configuration: StereoConfigurationForm| {
        product([
            CanonicalKeyValue::Unsigned(site_class.into()),
            sequence(
                ligands
                    .into_iter()
                    .map(|(class, kind)| stereo_ligand_key(class, kind)),
            ),
            stereo_configuration_form_key(&configuration),
        ])
    };

    match configuration {
        StereoConfigurationForm::Undetermined => {
            let mut ligands = ligand_classes.to_vec();
            ligands.sort_unstable();
            Ok(descriptor(ligands, StereoConfigurationForm::Undetermined))
        }
        StereoConfigurationForm::Kinded(kind, _) => stereo_frame_permutations(*kind)
            .map(|permutation| {
                configuration
                    .apply(permutation)
                    .ok_or(Contradiction)?
                    .normalize()
                    .map(|configuration| descriptor(permutation.act(ligand_classes), configuration))
            })
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .min()
            .ok_or(Contradiction),
    }
}

fn structure_partition_descriptors(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    adapter: &AutomorphismAdapter,
    entity_keys: &[InitialColorKey],
    entity_classes: &[u32],
) -> Result<Vec<InitialColorKey>, Contradiction> {
    let entity_class = |entity: Entity| entity_classes[incidence_graph.node_of(entity).index()];

    adapter
        .node_sources
        .iter()
        .map(|source| match *source {
            SubdivisionNodeSource::Node(node) => {
                let entity = incidence_graph.entity(node);
                let InitialColorKey::Entity { position, .. } = entity_keys[node.index()] else {
                    unreachable!("entity key corresponds to an entity node")
                };
                let value = match entity {
                    Entity::StereoAtom(id) => {
                        let stereo = molecule
                            .stereo_atoms()
                            .get(id)
                            .expect("structure incidence stereo atom is in range");
                        let ligand_classes = stereo
                            .ligands()
                            .map(|ligand| {
                                (entity_class(Entity::Atom(ligand.atom_id())), ligand.kind())
                            })
                            .collect::<Vec<_>>();
                        stereo_refinement_descriptor(
                            entity_class(Entity::Atom(stereo.site_id())),
                            &ligand_classes,
                            &stereo.attributes.configuration,
                        )?
                    }
                    Entity::StereoBond(id) => {
                        let stereo = molecule
                            .stereo_bonds()
                            .get(id)
                            .expect("structure incidence stereo bond is in range");
                        let ligand_classes = stereo
                            .ligands()
                            .map(|ligand| {
                                (entity_class(Entity::Atom(ligand.atom_id())), ligand.kind())
                            })
                            .collect::<Vec<_>>();
                        stereo_refinement_descriptor(
                            entity_class(Entity::Bond(stereo.site_id())),
                            &ligand_classes,
                            &stereo.attributes.configuration,
                        )?
                    }
                    _ => product([CanonicalKeyValue::Unsigned(entity_class(entity).into())]),
                };
                Ok(InitialColorKey::Entity { position, value })
            }
            SubdivisionNodeSource::Edge(edge) => {
                let value = match incidence_graph.incidence(edge) {
                    Incidence::DativeDonor | Incidence::DativeAcceptor => variant(1, []),
                    Incidence::AromaticParticipant(_) => variant(3, []),
                    Incidence::MulticenterParticipant(_) => variant(4, []),
                    incidence => incidence_key(incidence)?,
                };
                Ok(InitialColorKey::Incidence(value))
            }
        })
        .collect()
}

fn para_stereo_partition_descriptors(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    adapter: &AutomorphismAdapter,
    partition: &OrderedPartition,
) -> Result<Vec<(u32, Option<CanonicalKeyValue>)>, Contradiction> {
    let classes = partition.cell_indices(adapter.graph().node_count());
    let entity_class = |entity: Entity| classes[incidence_graph.node_of(entity).index()];
    let ligand_classes = |carrier: Entity| {
        let carrier_node = incidence_graph.node_of(carrier);
        let mut edges = incidence_graph
            .graph()
            .neighbors(carrier_node)
            .iter()
            .map(|neighbor| neighbor.edge)
            .collect::<Vec<_>>();
        edges.sort_unstable();
        edges
            .into_iter()
            .filter_map(|edge| {
                let Incidence::StereoLigand(kind) = incidence_graph.incidence(edge) else {
                    return None;
                };
                let occurrence = adapter
                    .node_of(SubdivisionNodeSource::Edge(edge))
                    .expect("stereo ligand incidence has an occurrence node");
                Some((classes[occurrence.index()], *kind))
            })
            .collect::<Vec<_>>()
    };

    adapter
        .node_sources
        .iter()
        .enumerate()
        .map(|(index, source)| {
            let descriptor = match *source {
                SubdivisionNodeSource::Node(node) => match incidence_graph.entity(node) {
                    Entity::StereoAtom(id) => {
                        let stereo = molecule
                            .stereo_atoms()
                            .get(id)
                            .expect("structure incidence stereo atom is in range");
                        Some(stereo_refinement_descriptor(
                            entity_class(Entity::Atom(stereo.site_id())),
                            &ligand_classes(Entity::StereoAtom(id)),
                            &stereo.attributes.configuration,
                        )?)
                    }
                    Entity::StereoBond(id) => {
                        let stereo = molecule
                            .stereo_bonds()
                            .get(id)
                            .expect("structure incidence stereo bond is in range");
                        Some(stereo_refinement_descriptor(
                            entity_class(Entity::Bond(stereo.site_id())),
                            &ligand_classes(Entity::StereoBond(id)),
                            &stereo.attributes.configuration,
                        )?)
                    }
                    _ => None,
                },
                SubdivisionNodeSource::Edge(_) => None,
            };
            Ok((classes[index], descriptor))
        })
        .collect()
}

/// Refine the structure partition, including para-stereo rounds when requested.
fn structure_partition(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    adapter: &AutomorphismAdapter,
    entity_keys: &[InitialColorKey],
    para_stereo: bool,
) -> Result<(OrderedPartition, usize), Contradiction> {
    let constitution_classes = constitution_entity_classes(molecule)?;
    let descriptors = structure_partition_descriptors(
        molecule,
        incidence_graph,
        adapter,
        entity_keys,
        &constitution_classes,
    )?;
    let mut partition = OrderedPartition::from_descriptors(&descriptors).refine(adapter.graph());
    let mut rounds = 1;
    let has_unresolved_stereo = |partition: &OrderedPartition| {
        partition.cells.iter().any(|cell| {
            cell.len() > 1
                && cell.iter().any(|&node| {
                    let SubdivisionNodeSource::Node(node) = adapter.node_source(node) else {
                        return false;
                    };
                    matches!(
                        incidence_graph.entity(node),
                        Entity::StereoAtom(_) | Entity::StereoBond(_)
                    )
                })
        })
    };

    if !para_stereo || !has_unresolved_stereo(&partition) {
        return Ok((partition, rounds));
    }

    loop {
        let cell_count = partition.cells.len();
        let descriptors =
            para_stereo_partition_descriptors(molecule, incidence_graph, adapter, &partition)?;
        let next = OrderedPartition::from_descriptors(&descriptors).refine(adapter.graph());
        rounds += 1;
        if next == partition {
            return Ok((partition, rounds));
        }
        debug_assert!(next.cells.len() > cell_count);
        partition = next;
        if !has_unresolved_stereo(&partition) {
            return Ok((partition, rounds));
        }
    }
}

/// Assign dense ordered colors to normalized entity and incidence keys.
fn rank_initial_colors(
    entity_keys: &[InitialColorKey],
    incidence_keys: &[InitialColorKey],
) -> InitialColors {
    let keys = entity_keys
        .iter()
        .chain(incidence_keys.iter())
        .cloned()
        .collect::<BTreeSet<_>>();
    let colors = keys
        .into_iter()
        .enumerate()
        .map(|(color, key)| (key, color as u32))
        .collect::<BTreeMap<_, _>>();

    InitialColors {
        entities: entity_keys.iter().map(|key| colors[key]).collect(),
        incidences: incidence_keys.iter().map(|key| colors[key]).collect(),
    }
}

fn atom_inherent_fields(attributes: &AtomForm) -> Result<Vec<FieldKey>, Contradiction> {
    Ok(vec![
        field(
            0,
            element_form_key(attributes.element.normalized()?.as_ref()),
        ),
        field(
            1,
            isotope_mass_form_key(attributes.isotope_mass.normalized()?.as_ref()),
        ),
        field(2, num_form_key(attributes.charge.normalized()?.as_ref())),
        field(
            3,
            num_form_key(attributes.implicit_hydrogens.normalized()?.as_ref()),
        ),
        field(
            4,
            num_form_key(attributes.lone_pairs.normalized()?.as_ref()),
        ),
        field(
            5,
            unpaired_electrons_form_key(attributes.unpaired_electrons.normalized()?.as_ref()),
        ),
    ])
}

fn bond_inherent_fields(attributes: &BondForm) -> Result<Vec<FieldKey>, Contradiction> {
    Ok(vec![
        field(1, num_form_key(attributes.order.normalized()?.as_ref())),
        field(2, num_form_key(attributes.charge.normalized()?.as_ref())),
        field(
            3,
            unpaired_electrons_form_key(attributes.unpaired_electrons.normalized()?.as_ref()),
        ),
    ])
}

fn entity_color_key(molecule: &Molecule, entity: Entity) -> Result<InitialColorKey, Contradiction> {
    let (position, value) = match entity {
        Entity::Atom(id) => {
            let attributes = molecule.atom(id).attributes;
            (
                EntityBlockPosition::ATOM,
                CanonicalKeyValue::Product(atom_inherent_fields(attributes)?),
            )
        }
        Entity::Bond(id) => {
            let attributes = molecule.bond(id).attributes;
            (
                EntityBlockPosition::BOND,
                CanonicalKeyValue::Product(bond_inherent_fields(attributes)?),
            )
        }
        Entity::DativeBond(id) => {
            let attributes = molecule
                .dative_bonds()
                .get(id)
                .expect("incidence dative bond is in range")
                .attributes;
            (
                EntityBlockPosition::DATIVE_BOND,
                positioned_product([(2, num_form_key(attributes.order.normalized()?.as_ref()))]),
            )
        }
        Entity::AromaticSystem(id) => {
            let attributes = molecule
                .aromatic_systems()
                .get(id)
                .expect("incidence aromatic system is in range")
                .attributes;
            (
                EntityBlockPosition::AROMATIC_SYSTEM,
                positioned_product([
                    (2, num_form_key(attributes.charge.normalized()?.as_ref())),
                    (
                        3,
                        unpaired_electrons_form_key(
                            attributes.unpaired_electrons.normalized()?.as_ref(),
                        ),
                    ),
                ]),
            )
        }
        Entity::MulticenterBond(id) => {
            let attributes = molecule
                .multicenter_bonds()
                .get(id)
                .expect("incidence multicenter bond is in range")
                .attributes;
            (
                EntityBlockPosition::MULTICENTER_BOND,
                positioned_product([
                    (2, num_form_key(attributes.charge.normalized()?.as_ref())),
                    (
                        3,
                        unpaired_electrons_form_key(
                            attributes.unpaired_electrons.normalized()?.as_ref(),
                        ),
                    ),
                ]),
            )
        }
        Entity::NoncovalentBond(id) => {
            let attributes = molecule
                .noncovalent_bonds()
                .get(id)
                .expect("incidence noncovalent bond is in range")
                .attributes;
            (
                EntityBlockPosition::NONCOVALENT_BOND,
                positioned_product([(
                    1,
                    noncovalent_bond_kind_form_key(attributes.kind.normalized()?.as_ref()),
                )]),
            )
        }
        Entity::StereoAtom(id) => {
            let attributes = molecule
                .stereo_atoms()
                .get(id)
                .expect("incidence stereo atom is in range")
                .attributes;
            (
                EntityBlockPosition::STEREO_ATOM,
                positioned_product([(
                    2,
                    option(attributes.configuration.kind().map(stereo_kind_key)),
                )]),
            )
        }
        Entity::StereoBond(id) => {
            let attributes = molecule
                .stereo_bonds()
                .get(id)
                .expect("incidence stereo bond is in range")
                .attributes;
            (
                EntityBlockPosition::STEREO_BOND,
                positioned_product([(
                    2,
                    option(attributes.configuration.kind().map(stereo_kind_key)),
                )]),
            )
        }
    };

    Ok(InitialColorKey::Entity { position, value })
}

fn reaction_span_entity_color_key(
    span: &ReactionSpan,
    entity: Entity,
) -> Result<InitialColorKey, Contradiction> {
    let (position, value) = match entity {
        Entity::Atom(id) => (
            EntityBlockPosition::ATOM,
            normalized_entity_span_key(&span.atoms()[id.index()], |attributes| {
                Ok(CanonicalKeyValue::Product(atom_inherent_fields(
                    attributes,
                )?))
            })?,
        ),
        Entity::Bond(id) => (
            EntityBlockPosition::BOND,
            normalized_entity_span_key(&span.bonds()[id.index()], |attributes| {
                Ok(CanonicalKeyValue::Product(bond_inherent_fields(
                    attributes,
                )?))
            })?,
        ),
        Entity::DativeBond(id) => (
            EntityBlockPosition::DATIVE_BOND,
            normalized_entity_span_key(span.dative_bonds().attributes(id), |attributes| {
                Ok(positioned_product([(
                    2,
                    num_form_key(attributes.order.normalized()?.as_ref()),
                )]))
            })?,
        ),
        Entity::AromaticSystem(id) => (
            EntityBlockPosition::AROMATIC_SYSTEM,
            normalized_entity_span_key(span.aromatic_systems().attributes(id), |attributes| {
                Ok(positioned_product([
                    (2, num_form_key(attributes.charge.normalized()?.as_ref())),
                    (
                        3,
                        unpaired_electrons_form_key(
                            attributes.unpaired_electrons.normalized()?.as_ref(),
                        ),
                    ),
                ]))
            })?,
        ),
        Entity::MulticenterBond(id) => (
            EntityBlockPosition::MULTICENTER_BOND,
            normalized_entity_span_key(span.multicenter_bonds().attributes(id), |attributes| {
                Ok(positioned_product([
                    (2, num_form_key(attributes.charge.normalized()?.as_ref())),
                    (
                        3,
                        unpaired_electrons_form_key(
                            attributes.unpaired_electrons.normalized()?.as_ref(),
                        ),
                    ),
                ]))
            })?,
        ),
        Entity::NoncovalentBond(id) => (
            EntityBlockPosition::NONCOVALENT_BOND,
            normalized_entity_span_key(span.noncovalent_bonds().attributes(id), |attributes| {
                Ok(positioned_product([(
                    1,
                    noncovalent_bond_kind_form_key(attributes.kind.normalized()?.as_ref()),
                )]))
            })?,
        ),
        Entity::StereoAtom(id) => (
            EntityBlockPosition::STEREO_ATOM,
            normalized_entity_span_key(span.stereo_atoms().attributes(id), |attributes| {
                Ok(positioned_product([(
                    2,
                    option(attributes.configuration.kind().map(stereo_kind_key)),
                )]))
            })?,
        ),
        Entity::StereoBond(id) => (
            EntityBlockPosition::STEREO_BOND,
            normalized_entity_span_key(span.stereo_bonds().attributes(id), |attributes| {
                Ok(positioned_product([(
                    2,
                    option(attributes.configuration.kind().map(stereo_kind_key)),
                )]))
            })?,
        ),
    };
    Ok(InitialColorKey::Entity { position, value })
}

fn topology_bond_key_rows(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    atom_images: &[usize],
) -> Result<Vec<(CanonicalKeyValue, BondId, NodeId)>, Contradiction> {
    let bond_count = incidence_graph.entity_count(EntityKind::Bond);
    let mut bonds = molecule
        .bonds()
        .iter()
        .map(|id| {
            let bond = id;
            let mut fields = Vec::with_capacity(4);
            let [first, second] = bond.atom_ids().map(|atom| atom_images[atom.index()] as u64);
            fields.push(field(
                0,
                product([
                    CanonicalKeyValue::Unsigned(first.min(second)),
                    CanonicalKeyValue::Unsigned(first.max(second)),
                ]),
            ));
            fields.extend(bond_inherent_fields(bond.attributes)?);
            Ok((
                CanonicalKeyValue::Product(fields),
                bond.id,
                incidence_graph.node_of(Entity::Bond(bond.id)),
            ))
        })
        .collect::<Result<Vec<_>, Contradiction>>()?;
    debug_assert_eq!(bonds.len(), bond_count);
    bonds.sort_unstable_by(|lhs, rhs| lhs.0.cmp(&rhs.0).then_with(|| lhs.1.cmp(&rhs.1)));
    Ok(bonds)
}

fn topology_candidate(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    order: &[NodeId],
) -> Result<CanonicalCandidate<CanonicalComparisonKey>, Contradiction> {
    let atom_count = incidence_graph.entity_count(EntityKind::Atom);
    let mut atom_images = vec![0_usize; atom_count];
    let mut atom_order = Vec::with_capacity(atom_count);

    for &node in order {
        if let Entity::Atom(id) = incidence_graph.entity(node) {
            atom_images[id.index()] = atom_order.len();
            atom_order.push(id);
        }
    }

    debug_assert_eq!(atom_order.len(), atom_count);

    let atoms = atom_order
        .iter()
        .copied()
        .map(|id| {
            atom_inherent_fields(molecule.atom(id).attributes).map(CanonicalKeyValue::Product)
        })
        .collect::<Result<Vec<_>, _>>()?;
    let bonds = topology_bond_key_rows(molecule, incidence_graph, &atom_images)?;

    let mut entity_order = atom_order
        .into_iter()
        .map(|id| incidence_graph.node_of(Entity::Atom(id)))
        .collect::<Vec<_>>();
    entity_order.extend(bonds.iter().map(|(_, _, node)| *node));
    let bonds = bonds.into_iter().map(|(key, _, _)| key).collect::<Vec<_>>();

    let mut entity_blocks = Vec::with_capacity(2);
    if !atoms.is_empty() {
        entity_blocks.push(PositionedKey {
            position: EntityBlockPosition::ATOM,
            value: CanonicalKeyValue::Sequence(atoms),
        });
    }
    if !bonds.is_empty() {
        entity_blocks.push(PositionedKey {
            position: EntityBlockPosition::BOND,
            value: CanonicalKeyValue::Sequence(bonds),
        });
    }

    Ok(CanonicalCandidate {
        key: CanonicalComparisonKey {
            entity_blocks,
            constraints: Vec::new(),
        },
        entity_order,
    })
}

fn topology_comparison_key_prefix(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    partition: &OrderedPartition,
) -> Result<CanonicalComparisonKeyPrefix, Contradiction> {
    let atom_count = incidence_graph.entity_count(EntityKind::Atom);
    let atom_order = partition
        .fixed_entity_prefix(incidence_graph.graph().node_count())
        .into_iter()
        .filter_map(|node| match incidence_graph.entity(node) {
            Entity::Atom(id) => Some(id),
            _ => None,
        })
        .collect::<Vec<_>>();
    let atoms = atom_order
        .iter()
        .copied()
        .map(|id| {
            atom_inherent_fields(molecule.atom(id).attributes).map(CanonicalKeyValue::Product)
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut prefix = CanonicalComparisonKeyPrefix::default();
    if !atoms.is_empty() {
        prefix.entity_blocks.push(PositionedKey {
            position: EntityBlockPosition::ATOM,
            value: CanonicalKeyValue::Sequence(atoms),
        });
    }
    if atom_order.len() != atom_count {
        return Ok(prefix);
    }

    let mut atom_images = vec![0_usize; atom_count];
    for (image, id) in atom_order.into_iter().enumerate() {
        atom_images[id.index()] = image;
    }
    let bonds = topology_bond_key_rows(molecule, incidence_graph, &atom_images)?
        .into_iter()
        .map(|(key, _, _)| key)
        .collect::<Vec<_>>();
    if !bonds.is_empty() {
        prefix.entity_blocks.push(PositionedKey {
            position: EntityBlockPosition::BOND,
            value: CanonicalKeyValue::Sequence(bonds),
        });
    }
    Ok(prefix)
}

fn constitution_candidate(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    order: &[NodeId],
) -> Result<CanonicalCandidate<CanonicalComparisonKey>, Contradiction> {
    fn electron_occurrence_fields(
        atom_ids: impl IntoIterator<Item = AtomId>,
        electrons: &ElectronCountsForm,
        atom_images: &[usize],
    ) -> [FieldKey; 2] {
        match electrons {
            ElectronCountsForm::Undetermined => {
                let mut participants = atom_ids
                    .into_iter()
                    .map(|atom| atom_images[atom.index()] as u64)
                    .collect::<Vec<_>>();
                participants.sort_unstable();
                [
                    field(
                        0,
                        sequence(participants.into_iter().map(CanonicalKeyValue::Unsigned)),
                    ),
                    field(1, electron_counts_form_key(electrons)),
                ]
            }
            ElectronCountsForm::Lit(electrons) => {
                let mut occurrences = atom_ids
                    .into_iter()
                    .zip(electrons)
                    .map(|(atom, &electrons)| (atom_images[atom.index()] as u64, electrons))
                    .collect::<Vec<_>>();
                occurrences.sort_unstable_by_key(|&(atom, _)| atom);
                let (participants, electrons): (Vec<_>, Vec<_>) = occurrences.into_iter().unzip();
                [
                    field(
                        0,
                        sequence(participants.into_iter().map(CanonicalKeyValue::Unsigned)),
                    ),
                    field(
                        1,
                        electron_counts_form_key(&ElectronCountsForm::Lit(electrons)),
                    ),
                ]
            }
        }
    }

    let mut candidate = topology_candidate(molecule, incidence_graph, order)?;
    let atom_count = incidence_graph.entity_count(EntityKind::Atom);
    let mut atom_images = vec![0_usize; atom_count];
    for (image, &node) in candidate.entity_order[..atom_count].iter().enumerate() {
        let Entity::Atom(id) = incidence_graph.entity(node) else {
            unreachable!("topology candidate begins with every atom")
        };
        atom_images[id.index()] = image;
    }

    let mut dative = molecule
        .dative_bonds()
        .iter()
        .map(|bond| {
            let mut donors = bond
                .donor_ids()
                .map(|atom| atom_images[atom.index()] as u64)
                .collect::<Vec<_>>();
            donors.sort_unstable();
            let fields = vec![
                field(
                    0,
                    sequence(donors.into_iter().map(CanonicalKeyValue::Unsigned)),
                ),
                field(
                    1,
                    CanonicalKeyValue::Unsigned(atom_images[bond.acceptor_id().index()] as u64),
                ),
                field(2, num_form_key(bond.order().normalized()?.as_ref())),
            ];
            Ok((
                CanonicalKeyValue::Product(fields),
                bond.id,
                incidence_graph.node_of(Entity::DativeBond(bond.id)),
            ))
        })
        .collect::<Result<Vec<_>, Contradiction>>()?;
    dative.sort_unstable_by(|lhs, rhs| lhs.0.cmp(&rhs.0).then_with(|| lhs.1.cmp(&rhs.1)));
    let dative = dative
        .into_iter()
        .map(|(key, _, node)| (key, node))
        .collect::<Vec<_>>();

    let mut aromatic = molecule
        .aromatic_systems()
        .iter()
        .map(|system| {
            let electrons = system.electrons().normalized()?;
            let mut fields =
                electron_occurrence_fields(system.atom_ids(), electrons.as_ref(), &atom_images)
                    .into_iter()
                    .collect::<Vec<_>>();
            fields.extend([
                field(2, num_form_key(system.charge().normalized()?.as_ref())),
                field(
                    3,
                    unpaired_electrons_form_key(system.unpaired_electrons().normalized()?.as_ref()),
                ),
            ]);
            Ok((
                CanonicalKeyValue::Product(fields),
                system.id,
                incidence_graph.node_of(Entity::AromaticSystem(system.id)),
            ))
        })
        .collect::<Result<Vec<_>, Contradiction>>()?;
    aromatic.sort_unstable_by(|lhs, rhs| lhs.0.cmp(&rhs.0).then_with(|| lhs.1.cmp(&rhs.1)));
    let aromatic = aromatic
        .into_iter()
        .map(|(key, _, node)| (key, node))
        .collect::<Vec<_>>();

    let mut multicenter = molecule
        .multicenter_bonds()
        .iter()
        .map(|bond| {
            let electrons = bond.electrons().normalized()?;
            let mut fields =
                electron_occurrence_fields(bond.atom_ids(), electrons.as_ref(), &atom_images)
                    .into_iter()
                    .collect::<Vec<_>>();
            fields.extend([
                field(2, num_form_key(bond.charge().normalized()?.as_ref())),
                field(
                    3,
                    unpaired_electrons_form_key(bond.unpaired_electrons().normalized()?.as_ref()),
                ),
            ]);
            Ok((
                CanonicalKeyValue::Product(fields),
                bond.id,
                incidence_graph.node_of(Entity::MulticenterBond(bond.id)),
            ))
        })
        .collect::<Result<Vec<_>, Contradiction>>()?;
    multicenter.sort_unstable_by(|lhs, rhs| lhs.0.cmp(&rhs.0).then_with(|| lhs.1.cmp(&rhs.1)));
    let multicenter = multicenter
        .into_iter()
        .map(|(key, _, node)| (key, node))
        .collect::<Vec<_>>();

    let mut noncovalent = molecule
        .noncovalent_bonds()
        .iter()
        .map(|bond| {
            let [first, second] = bond.atom_ids().map(|atom| atom_images[atom.index()] as u64);
            let fields = vec![
                field(
                    0,
                    product([
                        CanonicalKeyValue::Unsigned(first.min(second)),
                        CanonicalKeyValue::Unsigned(first.max(second)),
                    ]),
                ),
                field(
                    1,
                    noncovalent_bond_kind_form_key(bond.kind().normalized()?.as_ref()),
                ),
            ];
            Ok((
                CanonicalKeyValue::Product(fields),
                bond.id,
                incidence_graph.node_of(Entity::NoncovalentBond(bond.id)),
            ))
        })
        .collect::<Result<Vec<_>, Contradiction>>()?;
    noncovalent.sort_unstable_by(|lhs, rhs| lhs.0.cmp(&rhs.0).then_with(|| lhs.1.cmp(&rhs.1)));
    let noncovalent = noncovalent
        .into_iter()
        .map(|(key, _, node)| (key, node))
        .collect::<Vec<_>>();

    for (position, rows) in [
        (EntityBlockPosition::DATIVE_BOND, &dative),
        (EntityBlockPosition::AROMATIC_SYSTEM, &aromatic),
        (EntityBlockPosition::MULTICENTER_BOND, &multicenter),
        (EntityBlockPosition::NONCOVALENT_BOND, &noncovalent),
    ] {
        if !rows.is_empty() {
            candidate.key.entity_blocks.push(PositionedKey {
                position,
                value: sequence(rows.iter().map(|(key, _)| key.clone())),
            });
            candidate
                .entity_order
                .extend(rows.iter().map(|(_, node)| *node));
        }
    }

    Ok(candidate)
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct CanonicalStereoFrame {
    ligands: Vec<StereoLigand>,
    configuration: StereoConfigurationForm,
}

fn canonical_kinded_stereo_frame(
    ligands: &[StereoLigand],
    configuration: &StereoConfigurationForm,
) -> Result<Option<CanonicalStereoFrame>, Contradiction> {
    let Some(kind) = configuration.kind() else {
        return Ok(None);
    };
    if ligands.len() != kind.degree() {
        return Ok(None);
    }

    let mut minimum: Option<(Vec<StereoLigand>, StereoConfigurationForm)> = None;
    for permutation in stereo_frame_permutations(kind) {
        let candidate = (
            permutation.act(ligands),
            configuration
                .apply(permutation)
                .ok_or(Contradiction)?
                .normalize()?,
        );
        match minimum.as_ref().map(|value| candidate.cmp(value)) {
            None | Some(Ordering::Less) => {
                minimum = Some(candidate);
            }
            Some(Ordering::Equal | Ordering::Greater) => {}
        }
    }

    Ok(
        minimum.map(|(ligands, configuration)| CanonicalStereoFrame {
            ligands,
            configuration,
        }),
    )
}

fn structure_candidate(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    order: &[NodeId],
) -> Result<CanonicalCandidate<CanonicalComparisonKey>, Contradiction> {
    let mut candidate = constitution_candidate(molecule, incidence_graph, order)?;
    let mut atom_images = vec![0_usize; incidence_graph.entity_count(EntityKind::Atom)];
    let mut bond_images = vec![0_usize; incidence_graph.entity_count(EntityKind::Bond)];
    for (image, &node) in candidate.entity_order.iter().enumerate() {
        match incidence_graph.entity(node) {
            Entity::Atom(id) => atom_images[id.index()] = image,
            Entity::Bond(id) => {
                bond_images[id.index()] = image - atom_images.len();
            }
            _ => {}
        }
    }

    let remap_ligands = |ligands: Vec<StereoLigand>| {
        ligands
            .into_iter()
            .map(|ligand| {
                StereoLigand::new(
                    AtomId(atom_images[ligand.atom_id.index()] as u32),
                    ligand.kind,
                )
            })
            .collect::<Vec<_>>()
    };
    let canonical_frame =
        |ligands: Vec<StereoLigand>,
         configuration: &StereoConfigurationForm|
         -> Result<(Vec<StereoLigand>, StereoConfigurationForm), Contradiction> {
            let mut ligands = remap_ligands(ligands);
            match configuration {
                StereoConfigurationForm::Undetermined => {
                    ligands.sort_unstable();
                    Ok((ligands, StereoConfigurationForm::Undetermined))
                }
                StereoConfigurationForm::Kinded(..) => {
                    let frame = canonical_kinded_stereo_frame(&ligands, configuration)?
                        .expect("integrity established the kinded frame degree");
                    Ok((frame.ligands, frame.configuration))
                }
            }
        };

    let mut stereo_atoms = molecule
        .stereo_atoms()
        .iter()
        .map(|stereo| {
            let (ligands, configuration) =
                canonical_frame(stereo.ligand_frame(), &stereo.attributes.configuration)?;
            let fields = vec![
                field(
                    0,
                    CanonicalKeyValue::Unsigned(atom_images[stereo.site_id().index()] as u64),
                ),
                field(
                    1,
                    sequence(
                        ligands
                            .iter()
                            .map(|ligand| stereo_ligand_key(ligand.atom_id.0, ligand.kind)),
                    ),
                ),
                field(2, stereo_configuration_form_key(&configuration)),
            ];
            Ok((
                CanonicalKeyValue::Product(fields),
                stereo.id,
                incidence_graph.node_of(Entity::StereoAtom(stereo.id)),
            ))
        })
        .collect::<Result<Vec<_>, Contradiction>>()?;
    stereo_atoms.sort_unstable_by(|lhs, rhs| lhs.0.cmp(&rhs.0).then_with(|| lhs.1.cmp(&rhs.1)));
    let stereo_atoms = stereo_atoms
        .into_iter()
        .map(|(key, _, node)| (key, node))
        .collect::<Vec<_>>();

    let mut stereo_bonds = molecule
        .stereo_bonds()
        .iter()
        .map(|stereo| {
            let (ligands, configuration) =
                canonical_frame(stereo.ligand_frame(), &stereo.attributes.configuration)?;
            let fields = vec![
                field(
                    0,
                    CanonicalKeyValue::Unsigned(bond_images[stereo.site_id().index()] as u64),
                ),
                field(
                    1,
                    sequence(
                        ligands
                            .iter()
                            .map(|ligand| stereo_ligand_key(ligand.atom_id.0, ligand.kind)),
                    ),
                ),
                field(2, stereo_configuration_form_key(&configuration)),
            ];
            Ok((
                CanonicalKeyValue::Product(fields),
                stereo.id,
                incidence_graph.node_of(Entity::StereoBond(stereo.id)),
            ))
        })
        .collect::<Result<Vec<_>, Contradiction>>()?;
    stereo_bonds.sort_unstable_by(|lhs, rhs| lhs.0.cmp(&rhs.0).then_with(|| lhs.1.cmp(&rhs.1)));
    let stereo_bonds = stereo_bonds
        .into_iter()
        .map(|(key, _, node)| (key, node))
        .collect::<Vec<_>>();

    for (position, rows) in [
        (EntityBlockPosition::STEREO_ATOM, &stereo_atoms),
        (EntityBlockPosition::STEREO_BOND, &stereo_bonds),
    ] {
        if !rows.is_empty() {
            candidate.key.entity_blocks.push(PositionedKey {
                position,
                value: sequence(rows.iter().map(|(key, _)| key.clone())),
            });
            candidate
                .entity_order
                .extend(rows.iter().map(|(_, node)| *node));
        }
    }

    Ok(candidate)
}

fn complete_candidate(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    order: &[NodeId],
) -> Result<(CanonicalCandidate<CanonicalComparisonKey>, Molecule), Contradiction> {
    let mut candidate = structure_candidate(molecule, incidence_graph, order)?;
    let correspondence =
        correspondence_from_order(molecule, incidence_graph, &candidate.entity_order);
    let complete = molecule.remap(&correspondence).reframe()?;
    candidate.key.constraints = constraint_blocks(&complete);
    Ok((candidate, complete))
}

fn molecule_comparison_key_prefix(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    partition: &OrderedPartition,
    level: DescriptionLevel,
) -> Result<CanonicalComparisonKeyPrefix, Contradiction> {
    let topology_prefix = topology_comparison_key_prefix(molecule, incidence_graph, partition)?;
    if level == DescriptionLevel::Topology {
        return Ok(topology_prefix);
    }

    let atom_count = incidence_graph.entity_count(EntityKind::Atom);
    let atom_order = partition
        .fixed_entity_prefix(incidence_graph.graph().node_count())
        .into_iter()
        .filter(|&node| matches!(incidence_graph.entity(node), Entity::Atom(_)))
        .collect::<Vec<_>>();
    if atom_order.len() != atom_count {
        return Ok(topology_prefix);
    }

    let key = match level {
        DescriptionLevel::Topology => unreachable!("topology prefix returned above"),
        DescriptionLevel::Constitution => {
            constitution_candidate(molecule, incidence_graph, &atom_order)?.key
        }
        DescriptionLevel::Structure => {
            structure_candidate(molecule, incidence_graph, &atom_order)?.key
        }
        DescriptionLevel::Full => {
            complete_candidate(molecule, incidence_graph, &atom_order)?
                .0
                .key
        }
    };
    Ok(CanonicalComparisonKeyPrefix {
        entity_blocks: key.entity_blocks,
        constraints: key.constraints,
    })
}

fn molecule_prefix_worse(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    partition: &OrderedPartition,
    level: DescriptionLevel,
    incumbent: &CanonicalCandidate<CanonicalComparisonKey>,
) -> Result<bool, Contradiction> {
    Ok(
        molecule_comparison_key_prefix(molecule, incidence_graph, partition, level)?
            .cmp_key(&incumbent.key)
            == Ordering::Greater,
    )
}

fn molecule_counts(molecule: &Molecule) -> [usize; 8] {
    [
        molecule.atoms().count(),
        molecule.bonds().count(),
        molecule.dative_bonds().count(),
        molecule.aromatic_systems().count(),
        molecule.multicenter_bonds().count(),
        molecule.noncovalent_bonds().count(),
        molecule.stereo_atoms().count(),
        molecule.stereo_bonds().count(),
    ]
}

fn reaction_span_counts(span: &ReactionSpan) -> [usize; 8] {
    [
        span.atoms().len(),
        span.bonds().len(),
        span.dative_bonds().count(),
        span.aromatic_systems().count(),
        span.multicenter_bonds().count(),
        span.noncovalent_bonds().count(),
        span.stereo_atoms().count(),
        span.stereo_bonds().count(),
    ]
}

fn molecule_correspondence(images: &[Vec<usize>; 8]) -> MoleculeCorrespondence {
    fn correspondence<Id>(images: &[usize]) -> Correspondence<Id>
    where
        Id: Copy + Ord + From<usize>,
    {
        let images = images.iter().copied().map(Id::from).collect::<Vec<_>>();
        Correspondence::from_images(&images, images.len())
    }

    MoleculeCorrespondence::new(
        correspondence::<AtomId>(&images[0]),
        correspondence::<BondId>(&images[1]),
        correspondence::<DativeBondId>(&images[2]),
        correspondence::<AromaticSystemId>(&images[3]),
        correspondence::<MulticenterBondId>(&images[4]),
        correspondence::<NoncovalentBondId>(&images[5]),
        correspondence::<StereoAtomId>(&images[6]),
        correspondence::<StereoBondId>(&images[7]),
    )
}

fn correspondence_from_order(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    order: &[NodeId],
) -> MoleculeCorrespondence {
    let mut images = molecule_counts(molecule).map(|count| (0..count).collect::<Vec<_>>());
    let mut next = [0; 8];
    for &node in order {
        let (entity_kind, kind_id) = match incidence_graph.entity(node) {
            Entity::Atom(id) => (0, id.index()),
            Entity::Bond(id) => (1, id.index()),
            Entity::DativeBond(id) => (2, id.index()),
            Entity::AromaticSystem(id) => (3, id.index()),
            Entity::MulticenterBond(id) => (4, id.index()),
            Entity::NoncovalentBond(id) => (5, id.index()),
            Entity::StereoAtom(id) => (6, id.index()),
            Entity::StereoBond(id) => (7, id.index()),
        };
        images[entity_kind][kind_id] = next[entity_kind];
        next[entity_kind] += 1;
    }
    molecule_correspondence(&images)
}

fn lhs_anchored_correspondence_from_order(
    span: &ReactionSpan,
    incidence_graph: &IncidenceGraph,
    order: &[NodeId],
) -> MoleculeCorrespondence {
    let counts = reaction_span_counts(span);
    let mut lhs = counts.map(|_| Vec::new());
    let mut added = counts.map(|_| Vec::new());
    for &node in order {
        let entity = incidence_graph.entity(node);
        let (entity_kind, kind_id, lhs_present) = match entity {
            Entity::Atom(id) => (0, id.index(), span.atoms()[id.index()].lhs().is_some()),
            Entity::Bond(id) => (1, id.index(), span.bonds()[id.index()].lhs().is_some()),
            Entity::DativeBond(id) => (
                2,
                id.index(),
                span.dative_bonds().attributes(id).lhs().is_some(),
            ),
            Entity::AromaticSystem(id) => (
                3,
                id.index(),
                span.aromatic_systems().attributes(id).lhs().is_some(),
            ),
            Entity::MulticenterBond(id) => (
                4,
                id.index(),
                span.multicenter_bonds().attributes(id).lhs().is_some(),
            ),
            Entity::NoncovalentBond(id) => (
                5,
                id.index(),
                span.noncovalent_bonds().attributes(id).lhs().is_some(),
            ),
            Entity::StereoAtom(id) => (
                6,
                id.index(),
                span.stereo_atoms().attributes(id).lhs().is_some(),
            ),
            Entity::StereoBond(id) => (
                7,
                id.index(),
                span.stereo_bonds().attributes(id).lhs().is_some(),
            ),
        };
        if lhs_present {
            lhs[entity_kind].push(kind_id);
        } else {
            added[entity_kind].push(kind_id);
        }
    }
    for entity_kind in 0..counts.len() {
        let mut present = lhs[entity_kind]
            .iter()
            .chain(added[entity_kind].iter())
            .copied()
            .collect::<BTreeSet<_>>();
        for kind_id in 0..counts[entity_kind] {
            if present.insert(kind_id) {
                if reaction_span_lhs_present(span, entity_kind, kind_id) {
                    lhs[entity_kind].push(kind_id);
                } else {
                    added[entity_kind].push(kind_id);
                }
            }
        }
    }
    let mut images = counts.map(|count| vec![0; count]);
    for entity_kind in 0..images.len() {
        for (image, kind_id) in lhs[entity_kind]
            .iter()
            .chain(added[entity_kind].iter())
            .copied()
            .enumerate()
        {
            images[entity_kind][kind_id] = image;
        }
    }
    molecule_correspondence(&images)
}

fn reaction_span_lhs_present(span: &ReactionSpan, entity_kind: usize, kind_id: usize) -> bool {
    match entity_kind {
        0 => span.atoms()[kind_id].lhs().is_some(),
        1 => span.bonds()[kind_id].lhs().is_some(),
        2 => span
            .dative_bonds()
            .attributes(DativeBondId(kind_id as u32))
            .lhs()
            .is_some(),
        3 => span
            .aromatic_systems()
            .attributes(AromaticSystemId(kind_id as u32))
            .lhs()
            .is_some(),
        4 => span
            .multicenter_bonds()
            .attributes(MulticenterBondId(kind_id as u32))
            .lhs()
            .is_some(),
        5 => span
            .noncovalent_bonds()
            .attributes(NoncovalentBondId(kind_id as u32))
            .lhs()
            .is_some(),
        6 => span
            .stereo_atoms()
            .attributes(StereoAtomId(kind_id as u32))
            .lhs()
            .is_some(),
        7 => span
            .stereo_bonds()
            .attributes(StereoBondId(kind_id as u32))
            .lhs()
            .is_some(),
        _ => unreachable!("reaction-span entity-kind index is in range"),
    }
}

#[cfg(test)]
fn canonicalize_topology(
    molecule: &Molecule,
    context: &CanonicalizeContext,
) -> Result<Molecule, MoleculeCanonicalizeError> {
    canonicalize_topology_with_options(
        molecule,
        context,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            prefix_pruning: false,
            branch_order: backend_canonical_branch_order,
        },
    )
    .map(|(molecule, _)| molecule)
}

fn canonicalize_topology_with_options(
    molecule: &Molecule,
    context: &CanonicalizeContext,
    options: CanonicalSearchOptions,
) -> Result<(Molecule, MoleculeCorrespondence), MoleculeCanonicalizeError> {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Topology);
    let (entity_keys, incidence_keys) = initial_color_keys(molecule, &incidence_graph)?;
    let colors = rank_initial_colors(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let descriptors = partition_descriptors(&adapter, &entity_keys, &incidence_keys);
    let leaf_candidate = |order: &[NodeId]| {
        topology_candidate(molecule, &incidence_graph, order)
            .expect("initial colors established topology normalization")
    };
    let prefix_worse = |partition: &OrderedPartition, incumbent: &CanonicalCandidate<_>| {
        molecule_prefix_worse(
            molecule,
            &incidence_graph,
            partition,
            DescriptionLevel::Topology,
            incumbent,
        )
        .expect("initial colors established topology prefix normalization")
    };
    let selected = canonical_search(
        &adapter,
        &descriptors,
        context.automorphism_algorithm,
        options,
        &leaf_candidate,
        &prefix_worse,
    );
    let correspondence =
        correspondence_from_order(molecule, &incidence_graph, &selected.candidate.entity_order);

    let canonical = molecule.remap(&correspondence).reframe()?;
    Ok((canonical, correspondence))
}

#[cfg(test)]
fn canonicalize_constitution(
    molecule: &Molecule,
    context: &CanonicalizeContext,
) -> Result<Molecule, MoleculeCanonicalizeError> {
    canonicalize_constitution_with_options(
        molecule,
        context,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            prefix_pruning: false,
            branch_order: backend_canonical_branch_order,
        },
    )
    .map(|(molecule, _)| molecule)
}

fn canonicalize_constitution_with_options(
    molecule: &Molecule,
    context: &CanonicalizeContext,
    options: CanonicalSearchOptions,
) -> Result<(Molecule, MoleculeCorrespondence), MoleculeCanonicalizeError> {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Constitution);
    let (entity_keys, incidence_keys) = initial_color_keys(molecule, &incidence_graph)?;
    let colors = rank_initial_colors(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let descriptors = constitution_partition_descriptors(&adapter, &entity_keys, &incidence_graph);
    let leaf_candidate = |order: &[NodeId]| {
        constitution_candidate(molecule, &incidence_graph, order)
            .expect("initial colors established constitution normalization")
    };
    let prefix_worse = |partition: &OrderedPartition, incumbent: &CanonicalCandidate<_>| {
        molecule_prefix_worse(
            molecule,
            &incidence_graph,
            partition,
            DescriptionLevel::Constitution,
            incumbent,
        )
        .expect("initial colors established constitution prefix normalization")
    };
    let selected = canonical_search(
        &adapter,
        &descriptors,
        context.automorphism_algorithm,
        options,
        &leaf_candidate,
        &prefix_worse,
    );
    let correspondence =
        correspondence_from_order(molecule, &incidence_graph, &selected.candidate.entity_order);

    let canonical = molecule.remap(&correspondence).reframe()?;
    Ok((canonical, correspondence))
}

#[cfg(test)]
fn canonicalize_structure(
    molecule: &Molecule,
    context: &CanonicalizeContext,
) -> Result<Molecule, MoleculeCanonicalizeError> {
    canonicalize_structure_with_options(
        molecule,
        context,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            prefix_pruning: false,
            branch_order: backend_canonical_branch_order,
        },
    )
    .map(|(molecule, _)| molecule)
}

fn canonicalize_structure_with_options(
    molecule: &Molecule,
    context: &CanonicalizeContext,
    options: CanonicalSearchOptions,
) -> Result<(Molecule, MoleculeCorrespondence), MoleculeCanonicalizeError> {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let (entity_keys, incidence_keys) = initial_color_keys(molecule, &incidence_graph)?;
    let colors = rank_initial_colors(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let (partition, _) = structure_partition(
        molecule,
        &incidence_graph,
        &adapter,
        &entity_keys,
        context.para_stereo,
    )?;
    let descriptors = partition.cell_indices(adapter.graph().node_count());
    let leaf_candidate = |order: &[NodeId]| {
        structure_candidate(molecule, &incidence_graph, order)
            .expect("structure descriptors established entity normalization")
    };
    let prefix_worse = |partition: &OrderedPartition, incumbent: &CanonicalCandidate<_>| {
        molecule_prefix_worse(
            molecule,
            &incidence_graph,
            partition,
            DescriptionLevel::Structure,
            incumbent,
        )
        .expect("structure descriptors established prefix normalization")
    };
    let filter_generators = |generators: &mut Vec<Vec<NodeId>>| {
        retain_stereo_preserving_generators(molecule, &incidence_graph, generators);
    };
    let selected = canonical_search_with_generator_filter(
        &adapter,
        &descriptors,
        context.automorphism_algorithm,
        options,
        &leaf_candidate,
        &prefix_worse,
        &filter_generators,
    );
    let correspondence =
        correspondence_from_order(molecule, &incidence_graph, &selected.candidate.entity_order);

    let canonical = molecule.remap(&correspondence).reframe()?;
    Ok((canonical, correspondence))
}

#[cfg(test)]
fn canonicalize_full(
    molecule: &Molecule,
    context: &CanonicalizeContext,
) -> Result<Molecule, MoleculeCanonicalizeError> {
    canonicalize_full_with_options(
        molecule,
        context,
        CanonicalSearchOptions {
            automorphism_pruning: false,
            prefix_pruning: false,
            branch_order: backend_canonical_branch_order,
        },
    )
    .map(|(molecule, _)| molecule)
}

fn canonicalize_full_with_options(
    molecule: &Molecule,
    context: &CanonicalizeContext,
    options: CanonicalSearchOptions,
) -> Result<(Molecule, MoleculeCorrespondence), MoleculeCanonicalizeError> {
    let normalized = molecule.clone().normalize()?;
    let molecule = &normalized;
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let (entity_keys, incidence_keys) = initial_color_keys(molecule, &incidence_graph)?;
    let colors = rank_initial_colors(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let (partition, _) = structure_partition(
        molecule,
        &incidence_graph,
        &adapter,
        &entity_keys,
        context.para_stereo,
    )?;
    let descriptors = partition.cell_indices(adapter.graph().node_count());
    let leaf_candidate = |order: &[NodeId]| {
        complete_candidate(molecule, &incidence_graph, order)
            .expect("structure descriptors established complete normalization")
            .0
    };
    let prefix_worse = |partition: &OrderedPartition, incumbent: &CanonicalCandidate<_>| {
        molecule_prefix_worse(
            molecule,
            &incidence_graph,
            partition,
            DescriptionLevel::Full,
            incumbent,
        )
        .expect("structure descriptors established complete prefix normalization")
    };
    let selected = canonical_search(
        &adapter,
        &descriptors,
        context.automorphism_algorithm,
        CanonicalSearchOptions {
            automorphism_pruning: false,
            ..options
        },
        &leaf_candidate,
        &prefix_worse,
    );
    let correspondence =
        correspondence_from_order(molecule, &incidence_graph, &selected.candidate.entity_order);
    let (_, canonical) =
        complete_candidate(molecule, &incidence_graph, &selected.candidate.entity_order)?;
    Ok((canonical, correspondence))
}

fn canonical_key_by(
    molecule: &Molecule,
    level: DescriptionLevel,
    context: &CanonicalizeContext,
) -> Result<CanonicalComparisonKey, MoleculeCanonicalizeError> {
    if level == DescriptionLevel::Full {
        let normalized = molecule.clone().normalize()?;
        return canonical_key_by_full(&normalized, context);
    }
    let incidence_level = match level {
        DescriptionLevel::Topology => IncidenceLevel::Topology,
        DescriptionLevel::Constitution => IncidenceLevel::Constitution,
        DescriptionLevel::Structure => IncidenceLevel::Full,
        DescriptionLevel::Full => unreachable!("handled before incidence construction"),
    };
    let incidence_graph = molecule.incidence_graph(incidence_level);
    let (entity_keys, incidence_keys) = initial_color_keys(molecule, &incidence_graph)?;
    let colors = rank_initial_colors(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let prefix_worse = |partition: &OrderedPartition, incumbent: &CanonicalCandidate<_>| {
        molecule_prefix_worse(molecule, &incidence_graph, partition, level, incumbent)
            .expect("initial descriptors established comparison-key prefix normalization")
    };
    let options = CanonicalSearchOptions {
        automorphism_pruning: level != DescriptionLevel::Structure,
        prefix_pruning: false,
        branch_order: backend_canonical_branch_order,
    };

    let key = match level {
        DescriptionLevel::Topology => {
            let descriptors = partition_descriptors(&adapter, &entity_keys, &incidence_keys);
            let leaf_candidate = |order: &[NodeId]| {
                topology_candidate(molecule, &incidence_graph, order)
                    .expect("initial colors established topology normalization")
            };
            canonical_search(
                &adapter,
                &descriptors,
                context.automorphism_algorithm,
                options,
                &leaf_candidate,
                &prefix_worse,
            )
            .candidate
            .key
        }
        DescriptionLevel::Constitution => {
            let descriptors =
                constitution_partition_descriptors(&adapter, &entity_keys, &incidence_graph);
            let leaf_candidate = |order: &[NodeId]| {
                constitution_candidate(molecule, &incidence_graph, order)
                    .expect("initial colors established constitution normalization")
            };
            canonical_search(
                &adapter,
                &descriptors,
                context.automorphism_algorithm,
                options,
                &leaf_candidate,
                &prefix_worse,
            )
            .candidate
            .key
        }
        DescriptionLevel::Structure => {
            let (partition, _) = structure_partition(
                molecule,
                &incidence_graph,
                &adapter,
                &entity_keys,
                context.para_stereo,
            )?;
            let descriptors = partition.cell_indices(adapter.graph().node_count());
            let leaf_candidate = |order: &[NodeId]| {
                structure_candidate(molecule, &incidence_graph, order)
                    .expect("structure descriptors established entity normalization")
            };
            canonical_search(
                &adapter,
                &descriptors,
                context.automorphism_algorithm,
                options,
                &leaf_candidate,
                &prefix_worse,
            )
            .candidate
            .key
        }
        DescriptionLevel::Full => unreachable!("handled before selected-layer search"),
    };
    Ok(key)
}

fn canonical_key_by_full(
    molecule: &Molecule,
    context: &CanonicalizeContext,
) -> Result<CanonicalComparisonKey, MoleculeCanonicalizeError> {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let (entity_keys, incidence_keys) = initial_color_keys(molecule, &incidence_graph)?;
    let colors = rank_initial_colors(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let (partition, _) = structure_partition(
        molecule,
        &incidence_graph,
        &adapter,
        &entity_keys,
        context.para_stereo,
    )?;
    let descriptors = partition.cell_indices(adapter.graph().node_count());
    let leaf_candidate = |order: &[NodeId]| {
        complete_candidate(molecule, &incidence_graph, order)
            .expect("normalized structure descriptors established complete normalization")
            .0
    };
    let prefix_worse = |partition: &OrderedPartition, incumbent: &CanonicalCandidate<_>| {
        molecule_prefix_worse(
            molecule,
            &incidence_graph,
            partition,
            DescriptionLevel::Full,
            incumbent,
        )
        .expect("structure descriptors established complete prefix normalization")
    };
    Ok(canonical_search(
        &adapter,
        &descriptors,
        context.automorphism_algorithm,
        CanonicalSearchOptions {
            automorphism_pruning: false,
            prefix_pruning: false,
            branch_order: backend_canonical_branch_order,
        },
        &leaf_candidate,
        &prefix_worse,
    )
    .candidate
    .key)
}

fn reaction_span_entity_keys(
    span: &ReactionSpan,
    incidence_graph: &IncidenceGraph,
) -> Result<(Vec<InitialColorKey>, Vec<InitialColorKey>), Contradiction> {
    let entity_keys = incidence_graph
        .graph()
        .node_ids()
        .map(|node| reaction_span_entity_color_key(span, incidence_graph.entity(node)))
        .collect::<Result<Vec<_>, _>>()?;
    let incidence_keys = incidence_graph
        .incidences()
        .map(|(_, incidence)| incidence_key(incidence).map(InitialColorKey::Incidence))
        .collect::<Result<Vec<_>, _>>()?;
    Ok((entity_keys, incidence_keys))
}

fn reaction_span_comparison_key(
    span: &ReactionSpan,
    level: DescriptionLevel,
) -> Result<CanonicalComparisonKey, Contradiction> {
    let entries = span.entries();
    let atoms = entries
        .atoms
        .iter()
        .map(|span| {
            normalized_entity_span_key(span, |attributes| {
                Ok(CanonicalKeyValue::Product(atom_inherent_fields(
                    attributes,
                )?))
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let bonds = entries
        .bonds
        .iter()
        .map(|(first, second, span)| {
            let endpoints = product([
                index_key(first.index().min(second.index())),
                index_key(first.index().max(second.index())),
            ]);
            normalized_entity_span_key(span, |attributes| {
                Ok(product([
                    endpoints.clone(),
                    CanonicalKeyValue::Product(bond_inherent_fields(attributes)?),
                ]))
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut entity_blocks = Vec::new();
    push_span_block(&mut entity_blocks, EntityBlockPosition::ATOM, atoms);
    push_span_block(&mut entity_blocks, EntityBlockPosition::BOND, bonds);

    if matches!(
        level,
        DescriptionLevel::Constitution | DescriptionLevel::Structure | DescriptionLevel::Full
    ) {
        let dative = entries
            .dative
            .iter()
            .map(|(donors, acceptor, span)| {
                let mut donors = donors.iter().map(|id| id.index()).collect::<Vec<_>>();
                donors.sort_unstable();
                normalized_entity_span_key(span, |attributes| {
                    Ok(product([
                        sequence(donors.iter().copied().map(index_key)),
                        index_key(acceptor.index()),
                        num_form_key(attributes.order.normalized()?.as_ref()),
                    ]))
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let aromatic = entries
            .aromatic
            .iter()
            .map(|(atoms, span)| {
                relation_span_key(atoms, span, |attributes| {
                    Ok(product([
                        electron_counts_form_key(attributes.electrons.normalized()?.as_ref()),
                        num_form_key(attributes.charge.normalized()?.as_ref()),
                        unpaired_electrons_form_key(
                            attributes.unpaired_electrons.normalized()?.as_ref(),
                        ),
                    ]))
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let multicenter = entries
            .multicenter
            .iter()
            .map(|(atoms, span)| {
                relation_span_key(atoms, span, |attributes| {
                    Ok(product([
                        electron_counts_form_key(attributes.electrons.normalized()?.as_ref()),
                        num_form_key(attributes.charge.normalized()?.as_ref()),
                        unpaired_electrons_form_key(
                            attributes.unpaired_electrons.normalized()?.as_ref(),
                        ),
                    ]))
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let noncovalent = entries
            .noncovalent
            .iter()
            .map(|([first, second], span)| {
                normalized_entity_span_key(span, |attributes| {
                    Ok(product([
                        product([
                            index_key(first.index().min(second.index())),
                            index_key(first.index().max(second.index())),
                        ]),
                        noncovalent_bond_kind_form_key(attributes.kind.normalized()?.as_ref()),
                    ]))
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        push_span_block(&mut entity_blocks, EntityBlockPosition::DATIVE_BOND, dative);
        push_span_block(
            &mut entity_blocks,
            EntityBlockPosition::AROMATIC_SYSTEM,
            aromatic,
        );
        push_span_block(
            &mut entity_blocks,
            EntityBlockPosition::MULTICENTER_BOND,
            multicenter,
        );
        push_span_block(
            &mut entity_blocks,
            EntityBlockPosition::NONCOVALENT_BOND,
            noncovalent,
        );
    }

    if matches!(level, DescriptionLevel::Structure | DescriptionLevel::Full) {
        let stereo_atoms = entries
            .stereo_atoms
            .iter()
            .map(|(site, ligands, span)| {
                normalized_entity_span_key(span, |attributes| {
                    Ok(product([
                        index_key(site.index()),
                        sequence(
                            ligands
                                .iter()
                                .map(|ligand| stereo_ligand_key(ligand.atom_id.0, ligand.kind)),
                        ),
                        stereo_configuration_form_key(
                            attributes.configuration.normalized()?.as_ref(),
                        ),
                    ]))
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let stereo_bonds = entries
            .stereo_bonds
            .iter()
            .map(|(site, ligands, span)| {
                normalized_entity_span_key(span, |attributes| {
                    Ok(product([
                        index_key(site.index()),
                        sequence(
                            ligands
                                .iter()
                                .map(|ligand| stereo_ligand_key(ligand.atom_id.0, ligand.kind)),
                        ),
                        stereo_configuration_form_key(
                            attributes.configuration.normalized()?.as_ref(),
                        ),
                    ]))
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        push_span_block(
            &mut entity_blocks,
            EntityBlockPosition::STEREO_ATOM,
            stereo_atoms,
        );
        push_span_block(
            &mut entity_blocks,
            EntityBlockPosition::STEREO_BOND,
            stereo_bonds,
        );
    }

    let constraints = if level == DescriptionLevel::Full {
        let mut blocks = Vec::new();
        macro_rules! inline_span_block {
            ($position:expr, $entries:expr, $constraints:expr, $key:expr) => {{
                let rows = $entries
                    .iter()
                    .enumerate()
                    .filter_map(|(index, entry)| {
                        let span = $constraints(entry);
                        let value = entity_span_key(&span, Clone::clone);
                        let empty = span.lhs().is_none_or(is_empty_sequence)
                            && span.rhs().is_none_or(is_empty_sequence);
                        (!empty).then(|| product([index_key(index), value]))
                    })
                    .collect::<Vec<_>>();
                if !rows.is_empty() {
                    blocks.push(PositionedKey {
                        position: $position,
                        value: sequence(rows),
                    });
                }
            }};
        }
        inline_span_block!(
            ConstraintBlockPosition::ATOM,
            entries.atoms,
            |span: &EntitySpan<AtomForm>| constraint_key_span(span, |form| form
                .constraints
                .iter()
                .map(atom_constraint_form_key)),
            Clone::clone
        );
        inline_span_block!(
            ConstraintBlockPosition::BOND,
            entries.bonds,
            |entry: &(AtomId, AtomId, EntitySpan<BondForm>)| constraint_key_span(
                &entry.2,
                |form| form.constraints.iter().map(bond_constraint_form_key)
            ),
            Clone::clone
        );
        inline_span_block!(
            ConstraintBlockPosition::DATIVE_BOND,
            entries.dative,
            |entry: &(
                Vec<AtomId>,
                AtomId,
                EntitySpan<super::dative::DativeBondForm>
            )| constraint_key_span(&entry.2, |form| form
                .constraints
                .iter()
                .map(dative_bond_constraint_form_key)),
            Clone::clone
        );
        inline_span_block!(
            ConstraintBlockPosition::AROMATIC_SYSTEM,
            entries.aromatic,
            |entry: &(Vec<AtomId>, EntitySpan<super::aromatic::AromaticSystemForm>)| {
                constraint_key_span(&entry.1, |form| {
                    form.constraints
                        .iter()
                        .map(aromatic_system_constraint_form_key)
                })
            },
            Clone::clone
        );
        inline_span_block!(
            ConstraintBlockPosition::MULTICENTER_BOND,
            entries.multicenter,
            |entry: &(
                Vec<AtomId>,
                EntitySpan<super::multicenter::MulticenterBondForm>
            )| constraint_key_span(&entry.1, |form| form
                .constraints
                .iter()
                .map(multicenter_bond_constraint_form_key)),
            Clone::clone
        );
        inline_span_block!(
            ConstraintBlockPosition::NONCOVALENT_BOND,
            entries.noncovalent,
            |entry: &(
                [AtomId; 2],
                EntitySpan<super::noncovalent::NoncovalentBondForm>
            )| constraint_key_span(&entry.1, |form| form
                .constraints
                .iter()
                .map(noncovalent_bond_constraint_form_key)),
            Clone::clone
        );
        inline_span_block!(
            ConstraintBlockPosition::STEREO_ATOM,
            entries.stereo_atoms,
            |entry: &(AtomId, Vec<StereoLigand>, EntitySpan<StereoAtomForm>)| constraint_key_span(
                &entry.2,
                |form| form.constraints.iter().map(stereo_atom_constraint_form_key)
            ),
            Clone::clone
        );
        inline_span_block!(
            ConstraintBlockPosition::STEREO_BOND,
            entries.stereo_bonds,
            |entry: &(BondId, Vec<StereoLigand>, EntitySpan<StereoBondForm>)| constraint_key_span(
                &entry.2,
                |form| form.constraints.iter().map(stereo_bond_constraint_form_key)
            ),
            Clone::clone
        );

        let mut values = entries
            .constraints
            .iter()
            .map(constraint_span_key)
            .collect::<Vec<_>>();
        values.sort();
        values.dedup();
        if !values.is_empty() {
            blocks.push(PositionedKey {
                position: ConstraintBlockPosition::MOLECULE,
                value: sequence(values),
            });
        }
        blocks
    } else {
        Vec::new()
    };
    Ok(CanonicalComparisonKey {
        entity_blocks,
        constraints,
    })
}

fn is_empty_sequence(value: &CanonicalKeyValue) -> bool {
    matches!(value, CanonicalKeyValue::Sequence(values) if values.is_empty())
}

fn constraint_key_span<'a, T, I>(
    span: &'a EntitySpan<T>,
    keys: impl Fn(&'a T) -> I,
) -> EntitySpan<CanonicalKeyValue>
where
    I: IntoIterator<Item = CanonicalKeyValue>,
{
    let key = |value: &'a T| {
        let mut keys = keys(value).into_iter().collect::<Vec<_>>();
        keys.sort();
        keys.dedup();
        sequence(keys)
    };
    match span {
        EntitySpan::Unchanged(value) => EntitySpan::Unchanged(key(value)),
        EntitySpan::Added(value) => EntitySpan::Added(key(value)),
        EntitySpan::Removed(value) => EntitySpan::Removed(key(value)),
        EntitySpan::Modified { lhs, rhs } => EntitySpan::Modified {
            lhs: key(lhs),
            rhs: key(rhs),
        },
    }
}

fn push_span_block(
    blocks: &mut Vec<EntityBlockKey>,
    position: EntityBlockPosition,
    values: Vec<CanonicalKeyValue>,
) {
    if !values.is_empty() {
        blocks.push(PositionedKey {
            position,
            value: sequence(values),
        });
    }
}

fn relation_span_key<T>(
    atoms: &[AtomId],
    span: &EntitySpan<T>,
    attributes_key: impl Fn(&T) -> Result<CanonicalKeyValue, Contradiction>,
) -> Result<CanonicalKeyValue, Contradiction> {
    normalized_entity_span_key(span, |attributes| {
        Ok(product([
            sequence(atoms.iter().map(|id| index_key(id.index()))),
            attributes_key(attributes)?,
        ]))
    })
}

fn constraint_span_key(span: &ConstraintSpan) -> CanonicalKeyValue {
    match span {
        ConstraintSpan::Unchanged(value) => CanonicalKeyValue::Span(SpanKey {
            position: SpanTagPosition::UNCHANGED,
            values: vec![constraint_key(value)],
        }),
        ConstraintSpan::Added(value) => CanonicalKeyValue::Span(SpanKey {
            position: SpanTagPosition::ADDED,
            values: vec![constraint_key(value)],
        }),
        ConstraintSpan::Removed(value) => CanonicalKeyValue::Span(SpanKey {
            position: SpanTagPosition::REMOVED,
            values: vec![constraint_key(value)],
        }),
    }
}

fn canonicalize_molecule_with_correspondence_by_effective(
    molecule: &Molecule,
    level: DescriptionLevel,
    context: &CanonicalizeContext,
) -> Result<(Molecule, MoleculeCorrespondence), MoleculeCanonicalizeError> {
    let options = CanonicalSearchOptions {
        automorphism_pruning: matches!(
            level,
            DescriptionLevel::Topology
                | DescriptionLevel::Constitution
                | DescriptionLevel::Structure
        ),
        prefix_pruning: false,
        branch_order: backend_canonical_branch_order,
    };
    match level {
        DescriptionLevel::Topology => {
            canonicalize_topology_with_options(molecule, context, options)
        }
        DescriptionLevel::Constitution => {
            canonicalize_constitution_with_options(molecule, context, options)
        }
        DescriptionLevel::Structure => {
            canonicalize_structure_with_options(molecule, context, options)
        }
        DescriptionLevel::Full => canonicalize_full_with_options(molecule, context, options),
    }
}

fn canonicalize_molecule_by_effective(
    molecule: &Molecule,
    level: DescriptionLevel,
    context: &CanonicalizeContext,
) -> Result<Molecule, MoleculeCanonicalizeError> {
    canonicalize_molecule_with_correspondence_by_effective(molecule, level, context)
        .map(|(canonical, _)| canonical)
}

impl Canonicalize for Molecule {
    type Error = MoleculeCanonicalizeError;

    fn canonicalize(self, context: &CanonicalizeContext) -> Result<Self, Self::Error> {
        let level = molecule_canonicalize_level(&self);
        canonicalize_molecule_by_effective(&self, level, context)
    }

    fn canonicalize_with_correspondence(
        self,
        context: &CanonicalizeContext,
    ) -> Result<(Self, MoleculeCorrespondence), Self::Error> {
        let level = molecule_canonicalize_level(&self);
        canonicalize_molecule_with_correspondence_by_effective(&self, level, context)
    }

    fn canonical_eq(&self, other: &Self, context: &CanonicalizeContext) -> bool {
        if self == other {
            return true;
        }
        let level = molecule_canonicalize_level(self).max(molecule_canonicalize_level(other));
        match (
            canonical_key_by(self, level, context),
            canonical_key_by(other, level, context),
        ) {
            (Ok(left), Ok(right)) => left == right,
            (
                Err(MoleculeCanonicalizeError::Contradiction(_)),
                Err(MoleculeCanonicalizeError::Contradiction(_)),
            ) => true,
            _ => false,
        }
    }
}

/// Failure to construct a canonical [`Molecule`].
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum MoleculeCanonicalizeError {
    /// Intrinsic normalization of a carried value reached a contradiction.
    #[error(transparent)]
    Contradiction(#[from] Contradiction),
}

/// Failure to construct a canonical [`ReactionSpan`].
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ReactionSpanCanonicalizeError {
    /// Intrinsic normalization of a carried value reached a contradiction.
    #[error(transparent)]
    Contradiction(#[from] Contradiction),
}

fn reaction_span_canonical_candidate(
    span: &ReactionSpan,
    level: DescriptionLevel,
    context: &CanonicalizeContext,
) -> Result<CanonicalCandidate<CanonicalComparisonKey>, Contradiction> {
    let incidence_level = match level {
        DescriptionLevel::Topology => IncidenceLevel::Topology,
        DescriptionLevel::Constitution => IncidenceLevel::Constitution,
        DescriptionLevel::Structure | DescriptionLevel::Full => IncidenceLevel::Full,
    };
    let incidence_graph = span.incidence_graph(incidence_level);
    let (entity_keys, incidence_keys) = reaction_span_entity_keys(span, &incidence_graph)?;
    let colors = rank_initial_colors(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let descriptors = match level {
        DescriptionLevel::Topology => {
            partition_descriptors(&adapter, &entity_keys, &incidence_keys)
        }
        DescriptionLevel::Constitution | DescriptionLevel::Structure | DescriptionLevel::Full => {
            constitution_partition_descriptors(&adapter, &entity_keys, &incidence_graph)
        }
    };
    let leaf_candidate = |order: &[NodeId]| {
        let correspondence = lhs_anchored_correspondence_from_order(span, &incidence_graph, order);
        let remapped = span.remap(&correspondence);
        let action = remapped.representative_action();
        let reframed = remapped
            .reframe_by(&action)
            .expect("integrity established compatible reaction-span frame actions");
        CanonicalCandidate {
            key: reaction_span_comparison_key(&reframed, level)
                .expect("initial colors established span normalization"),
            entity_order: order.to_vec(),
        }
    };
    let no_prefix = |_: &OrderedPartition, _: &CanonicalCandidate<_>| false;
    Ok(canonical_search(
        &adapter,
        &descriptors,
        context.automorphism_algorithm,
        CanonicalSearchOptions {
            automorphism_pruning: false,
            prefix_pruning: false,
            branch_order: backend_canonical_branch_order,
        },
        &leaf_candidate,
        &no_prefix,
    )
    .candidate)
}

fn reaction_span_from_candidate(
    span: &ReactionSpan,
    level: DescriptionLevel,
    candidate: &CanonicalCandidate<CanonicalComparisonKey>,
) -> Result<ReactionSpan, Contradiction> {
    let incidence_level = match level {
        DescriptionLevel::Topology => IncidenceLevel::Topology,
        DescriptionLevel::Constitution => IncidenceLevel::Constitution,
        DescriptionLevel::Structure | DescriptionLevel::Full => IncidenceLevel::Full,
    };
    let incidence_graph = span.incidence_graph(incidence_level);
    let correspondence =
        lhs_anchored_correspondence_from_order(span, &incidence_graph, &candidate.entity_order);
    span.remap(&correspondence).reframe()
}

fn canonicalize_reaction_span_by(
    span: &ReactionSpan,
    level: DescriptionLevel,
    context: &CanonicalizeContext,
) -> Result<ReactionSpan, ReactionSpanCanonicalizeError> {
    Ok(canonicalize_checked_reaction_span_by(span, level, context)?)
}

fn canonicalize_checked_reaction_span_by(
    span: &ReactionSpan,
    level: DescriptionLevel,
    context: &CanonicalizeContext,
) -> Result<ReactionSpan, Contradiction> {
    Ok(canonicalize_checked_reaction_span_with_correspondence_by(span, level, context)?.0)
}

fn canonicalize_checked_reaction_span_with_correspondence_by(
    span: &ReactionSpan,
    level: DescriptionLevel,
    context: &CanonicalizeContext,
) -> Result<(ReactionSpan, MoleculeCorrespondence), Contradiction> {
    let normalized = span.clone().normalize()?;
    let candidate = reaction_span_canonical_candidate(&normalized, level, context)?;
    let incidence_level = match level {
        DescriptionLevel::Topology => IncidenceLevel::Topology,
        DescriptionLevel::Constitution => IncidenceLevel::Constitution,
        DescriptionLevel::Structure | DescriptionLevel::Full => IncidenceLevel::Full,
    };
    let incidence_graph = normalized.incidence_graph(incidence_level);
    let correspondence = lhs_anchored_correspondence_from_order(
        &normalized,
        &incidence_graph,
        &candidate.entity_order,
    );
    Ok((
        reaction_span_from_candidate(&normalized, level, &candidate)?,
        correspondence,
    ))
}

fn canonical_reaction_span_key(
    span: &ReactionSpan,
    level: DescriptionLevel,
    context: &CanonicalizeContext,
) -> Result<CanonicalComparisonKey, ReactionSpanCanonicalizeError> {
    if level == DescriptionLevel::Full {
        let normalized = span.clone().normalize()?;
        return Ok(reaction_span_canonical_candidate(&normalized, level, context)?.key);
    }
    Ok(reaction_span_canonical_candidate(span, level, context)?.key)
}

fn canonicalize_reaction_span_by_effective(
    span: &ReactionSpan,
    level: DescriptionLevel,
    context: &CanonicalizeContext,
) -> Result<ReactionSpan, ReactionSpanCanonicalizeError> {
    canonicalize_reaction_span_by(span, level, context)
}

fn canonicalize_reaction_span_with_correspondence_by_effective(
    span: &ReactionSpan,
    level: DescriptionLevel,
    context: &CanonicalizeContext,
) -> Result<(ReactionSpan, MoleculeCorrespondence), ReactionSpanCanonicalizeError> {
    Ok(canonicalize_checked_reaction_span_with_correspondence_by(
        span, level, context,
    )?)
}

impl Canonicalize for ReactionSpan {
    type Error = ReactionSpanCanonicalizeError;

    fn canonicalize(self, context: &CanonicalizeContext) -> Result<Self, Self::Error> {
        let level = reaction_span_canonicalize_level(&self);
        canonicalize_reaction_span_by_effective(&self, level, context)
    }

    fn canonicalize_with_correspondence(
        self,
        context: &CanonicalizeContext,
    ) -> Result<(Self, MoleculeCorrespondence), Self::Error> {
        let level = reaction_span_canonicalize_level(&self);
        canonicalize_reaction_span_with_correspondence_by_effective(&self, level, context)
    }

    fn canonical_eq(&self, other: &Self, context: &CanonicalizeContext) -> bool {
        if self == other {
            return true;
        }
        let level =
            reaction_span_canonicalize_level(self).max(reaction_span_canonicalize_level(other));
        match (
            canonical_reaction_span_key(self, level, context),
            canonical_reaction_span_key(other, level, context),
        ) {
            (Ok(left), Ok(right)) => left == right,
            (
                Err(ReactionSpanCanonicalizeError::Contradiction(_)),
                Err(ReactionSpanCanonicalizeError::Contradiction(_)),
            ) => true,
            _ => false,
        }
    }
}

/// Failure to construct a canonical [`Reaction`].
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ReactionCanonicalizeError {
    /// Intrinsic normalization or span materialization reached a contradiction.
    #[error(transparent)]
    Contradiction(#[from] Contradiction),
}

fn canonicalize_reaction_by(
    reaction: &Reaction,
    level: DescriptionLevel,
    context: &CanonicalizeContext,
) -> Result<Reaction, ReactionCanonicalizeError> {
    let span = reaction.to_reaction_span()?;
    Ok(canonicalize_checked_reaction_span_by(&span, level, context)?.to_reaction())
}

fn canonicalize_reaction_by_effective(
    reaction: &Reaction,
    level: DescriptionLevel,
    context: &CanonicalizeContext,
) -> Result<Reaction, ReactionCanonicalizeError> {
    canonicalize_reaction_by(reaction, level, context)
}

fn canonicalize_reaction_with_correspondence_by_effective(
    reaction: &Reaction,
    level: DescriptionLevel,
    context: &CanonicalizeContext,
) -> Result<(Reaction, MoleculeCorrespondence), ReactionCanonicalizeError> {
    let span = reaction.to_reaction_span()?;
    let (canonical, correspondence) =
        canonicalize_checked_reaction_span_with_correspondence_by(&span, level, context)?;
    Ok((canonical.to_reaction(), correspondence))
}

impl Canonicalize for Reaction {
    type Error = ReactionCanonicalizeError;

    fn canonicalize(self, context: &CanonicalizeContext) -> Result<Self, Self::Error> {
        let level = reaction_canonicalize_level(&self);
        canonicalize_reaction_by_effective(&self, level, context)
    }

    fn canonicalize_with_correspondence(
        self,
        context: &CanonicalizeContext,
    ) -> Result<(Self, MoleculeCorrespondence), Self::Error> {
        let level = reaction_canonicalize_level(&self);
        canonicalize_reaction_with_correspondence_by_effective(&self, level, context)
    }

    fn canonical_eq(&self, other: &Self, context: &CanonicalizeContext) -> bool {
        if self == other {
            return true;
        }
        let level = reaction_canonicalize_level(self).max(reaction_canonicalize_level(other));
        match (
            canonicalize_reaction_by_effective(self, level, context),
            canonicalize_reaction_by_effective(other, level, context),
        ) {
            (Ok(left), Ok(right)) => left == right,
            (
                Err(ReactionCanonicalizeError::Contradiction(_)),
                Err(ReactionCanonicalizeError::Contradiction(_)),
            ) => true,
            _ => false,
        }
    }
}

#[cfg(test)]
mod tests;
