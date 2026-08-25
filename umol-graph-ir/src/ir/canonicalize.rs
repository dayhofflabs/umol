//! Aggregate canonicalization inputs and failures.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::hash::{DefaultHasher, Hash, Hasher};

use thiserror::Error;
use umol_graph_core::{
    AutomorphismAlgorithm, AutomorphismOutput, Correspondence, FactorOrdering, Graph, NodeId,
    ParticipantPosition, RelationId, SubdivisionNodeSource, Unordered,
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
    StereoAtomConstraintForm, StereoAtomConstraintsForm, StereoBondConstraintForm,
    StereoBondConstraintsForm, StereoLigandPair, StereogenicityForm, TopicityForm,
    TopicityRelationForm,
};
use super::correspondence::MoleculeCorrespondence;
use super::delta::{
    AromaticSystemDelta, AtomDelta, BondDelta, ConstraintSpan, DativeBondDelta, Delta, Deltas,
    EntitySpan, MulticenterBondDelta, NoncovalentBondDelta, StereoAtomDelta, StereoBondDelta,
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
use super::molecule::{Molecule, MoleculeEntries, MoleculeIntegrityError};
use super::noncovalent::{NoncovalentBondKind, NoncovalentBondKindForm};
use super::num::{ArithExpr, NumForm, PredExpr};
use super::operators::{MemOp, RelOp};
use super::reaction::{Reaction, ReactionIntegrityError};
use super::reaction_span::{ReactionSpan, ReactionSpanIntegrityError};
use super::spin::UnpairedElectronsForm;
use super::stereo::{
    CisTransStereoForm, StereoAtomForm, StereoBondForm, StereoConfigurationForm, StereoCoset,
    StereoKind, StereoTerm, Stereogenicity, TetrahedralStereoForm, Topicity,
};
use super::traits::Normalize;

/// Semantic and operational inputs to aggregate canonicalization.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CanonicalizeContext {
    /// Whether stereo-sensitive refinement is iterated to a para-stereo fixpoint.
    pub para_stereo: bool,
    /// Graph automorphism algorithm used during canonical-frame search.
    pub automorphism_algorithm: AutomorphismAlgorithm,
}

/// Level used to select or compare an aggregate's entity frame.
///
/// The variants form a nested hierarchy. [`Topology`](Self::Topology) contains atoms and localized
/// bonds. [`Constitution`](Self::Constitution) additionally contains dative bonds, aromatic
/// systems, multicenter bonds, and noncovalent bonds. [`Structure`](Self::Structure) additionally
/// contains stereo atoms and stereo bonds. [`Full`](Self::Full) additionally uses normalized
/// constraints to select among tied structure frames. In entity-domain terminology, topology is
/// AB, non-stereo is DAMN, stereo is SS, constitution is topology plus non-stereo, and overlays are
/// non-stereo plus stereo.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CanonicalizeLevel {
    Topology,
    Constitution,
    Structure,
    Full,
}

/// Canonical entity-frame selection for complete indexed graph-IR aggregates.
///
/// Unlike [`Normalize`], canonicalization may change entity ids and participant frames. It
/// preserves the represented aggregate, transports every reference and position-sensitive value,
/// and normalizes carried forms in the selected frame. Construction of the canonical form is
/// fallible; canonical equality totalizes those failures according to the aggregate's semantic
/// contract.
///
/// # Semantic properties
///
/// For every integrity-valid aggregate and fixed context:
///
/// - successful complete canonicalization is exactly idempotent and invariant under valid dense
///   entity remapping;
/// - `canonical_eq` and each `canonical_eq_by` relation are reflexive, symmetric, and transitive
///   under their documented failure totalization;
/// - successful canonical hashes are invariant under valid dense entity remapping, and canonical
///   equality implies equal canonical hashes when both hash operations succeed;
/// - [`CanonicalizeLevel::Full`] is identical to the corresponding unqualified operation;
/// - `canonical_eq_by` is invariant under valid dense entity remapping at every level.
///
/// A level-specific transformation returns a complete aggregate, but features excluded from frame
/// selection have no promised ordering. Therefore only its selected-layer equality, not the
/// complete returned representation, is invariant under remapping.
///
/// For a fixed umol release, level, and context, canonicalization is deterministic. During the 0.x
/// series, the typed comparison schema and resulting representatives may change between releases.
/// Returned canonical forms are ordinary IR values without schema-version provenance and must not
/// be used as persistent cross-release identifiers.
pub trait Canonicalize: Sized {
    type Error;

    /// Construct the complete canonical form.
    ///
    /// # Errors
    ///
    /// Returns the aggregate-specific integrity error for malformed representation state and
    /// [`Contradiction`] when intrinsic normalization is unsatisfiable.
    fn canonicalize(self, context: &CanonicalizeContext) -> Result<Self, Self::Error>;

    /// Construct the complete canonical form and its source-to-canonical correspondence.
    ///
    /// The correspondence maps every entity id in the input frame to its id in the returned
    /// canonical frame. Each of its eight entity-family components is total on both sides and
    /// therefore represents a dense bijection. For [`Reaction`], the two frames are the complete
    /// union frames of the materialized input and returned reaction spans.
    ///
    /// # Errors
    ///
    /// Returns the same aggregate-specific integrity or intrinsic-normalization error as
    /// [`Self::canonicalize`]. A [`Reaction`] also reports failure to materialize its reaction
    /// span.
    ///
    /// # Semantic properties
    ///
    /// Discarding the correspondence yields exactly [`Self::canonicalize`] under the same
    /// context. Transporting the input through the correspondence preserves its represented
    /// aggregate in the returned canonical frame.
    fn canonicalize_with_correspondence(
        self,
        context: &CanonicalizeContext,
    ) -> Result<(Self, MoleculeCorrespondence), Self::Error>;

    /// Construct a complete normalized aggregate whose selected structural layer is canonical.
    ///
    /// Features excluded by `level` are preserved but do not break ties in the selected frame.
    /// [`CanonicalizeLevel::Full`] is identical to [`Self::canonicalize`].
    ///
    /// # Errors
    ///
    /// Returns the aggregate-specific integrity error for malformed representation state and
    /// [`Contradiction`] when normalization of any carried value is unsatisfiable, including a
    /// value excluded from frame selection.
    fn canonicalize_by(
        self,
        level: CanonicalizeLevel,
        context: &CanonicalizeContext,
    ) -> Result<Self, Self::Error>;

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

    /// Hash the canonical comparison key at the selected structural layer.
    ///
    /// Excluded features do not contribute to the hash. The returned value uses Rust's
    /// [`DefaultHasher`], may collide, and is not a persistent identifier.
    /// [`CanonicalizeLevel::Full`] is identical to [`Self::canonical_hash`].
    ///
    /// # Errors
    ///
    /// Returns the aggregate-specific integrity error for malformed representation state and
    /// [`Contradiction`] when normalization of selected data is unsatisfiable. At reduced levels,
    /// contradictions confined to excluded features do not affect the hash. A [`Reaction`] may
    /// also report failure to materialize its selected reaction span.
    ///
    /// # Semantic properties
    ///
    /// When two values compare equal through [`Self::canonical_eq_by`] and both hashes succeed,
    /// their canonical hashes at that level are equal. The hash is invariant under valid dense
    /// entity remapping at every level.
    fn canonical_hash_by(
        self,
        level: CanonicalizeLevel,
        context: &CanonicalizeContext,
    ) -> Result<u64, Self::Error>
    where
        Self: Hash;

    /// Compare complete canonical forms.
    ///
    /// Structural identity short-circuits the operation. Otherwise, two intrinsic contradictions
    /// compare equal, while an integrity failure never makes distinct inputs equal.
    fn canonical_eq(&self, other: &Self, context: &CanonicalizeContext) -> bool;

    /// Compare canonical forms at the selected structural layer.
    ///
    /// Contradictions in features excluded by `level` do not affect this relation.
    /// [`CanonicalizeLevel::Full`] is identical to [`Self::canonical_eq`].
    fn canonical_eq_by(
        &self,
        other: &Self,
        level: CanonicalizeLevel,
        context: &CanonicalizeContext,
    ) -> bool;
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
        EntitySpan::Modified { lhs, rhs } => (
            SpanTagPosition::MODIFIED,
            vec![value_key(lhs)?, value_key(rhs)?],
        ),
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
enum InitialClassKey {
    Entity {
        position: EntityBlockPosition,
        value: CanonicalKeyValue,
    },
    Incidence(CanonicalKeyValue),
}

impl Ord for InitialClassKey {
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

impl PartialOrd for InitialClassKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct InitialClasses {
    entities: Vec<u32>,
    incidences: Vec<u32>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum AutomorphismClass {
    Entity(u32),
    Incidence(u32),
}

#[derive(Clone, Debug)]
struct AutomorphismAdapter {
    // Source entity nodes retain their ids. Role- or value-bearing incidence edges become colored
    // occurrence nodes; single-role endpoints remain direct unless duplicates require subdivision.
    graph: Graph,
    classes: Vec<AutomorphismClass>,
    node_sources: Vec<SubdivisionNodeSource>,
    incidence_nodes: Vec<Option<NodeId>>,
    source_node_count: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ProjectedAutomorphismOutput {
    orbits: Vec<NodeId>,
    // Backend canonical labels are branch-order hints, not the canonical molecule numbering.
    canonical_labels: Vec<NodeId>,
    generators: Vec<Vec<NodeId>>,
}

impl AutomorphismAdapter {
    /// Construct the exact graph adapter used by canonicalization.
    fn new(incidence_graph: &IncidenceGraph, initial_classes: &InitialClasses) -> Self {
        let source = incidence_graph.graph();
        debug_assert_eq!(initial_classes.entities.len(), source.node_count());
        debug_assert_eq!(initial_classes.incidences.len(), source.edge_count());

        let mut classes = initial_classes
            .entities
            .iter()
            .copied()
            .map(AutomorphismClass::Entity)
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
            classes.push(AutomorphismClass::Incidence(
                initial_classes.incidences[edge.index()],
            ));
            push_edge([endpoints[0], occurrence]);
            push_edge([occurrence, endpoints[1]]);
        }
        let graph = Graph::new(node_sources.len(), &edges);
        debug_assert!(graph.is_simple());

        Self {
            graph,
            classes,
            node_sources,
            incidence_nodes,
            source_node_count: source.node_count(),
        }
    }

    fn graph(&self) -> &Graph {
        &self.graph
    }

    fn class(&self, node: NodeId) -> AutomorphismClass {
        self.classes[node.index()]
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
    ) -> ProjectedAutomorphismOutput {
        let cell_indices = partition.cell_indices(self.graph().node_count());
        // Search partitions may deliberately coarsen covariant occurrence data. Retaining the
        // exact adapter class here keeps orbit pruning restricted to true automorphisms.
        let output = self.graph().automorphisms(
            |node| (self.class(node), cell_indices[node.index()]),
            algorithm,
        );

        self.project_automorphisms(&output)
    }

    fn project_automorphisms(&self, output: &AutomorphismOutput) -> ProjectedAutomorphismOutput {
        let source_node = |node| match self.node_source(node) {
            SubdivisionNodeSource::Node(source) => source,
            SubdivisionNodeSource::Edge(_) => {
                unreachable!("disjoint classes preserve the adapter node domain")
            }
        };
        let orbits = (0..self.source_node_count)
            .map(|index| source_node(output.orbit_of(NodeId(index as u32))))
            .collect();
        let canonical_labels = output
            .canonical_labels()
            .iter()
            .filter_map(|&node| match self.node_source(node) {
                SubdivisionNodeSource::Node(source) => Some(source),
                SubdivisionNodeSource::Edge(_) => None,
            })
            .collect();
        let generators = output
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

        ProjectedAutomorphismOutput {
            orbits,
            canonical_labels,
            generators,
        }
    }
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
    Option<&ProjectedAutomorphismOutput>,
    &mut [NodeId],
) -> bool;

fn backend_canonical_branch_order(
    adapter: &AutomorphismAdapter,
    partition: &OrderedPartition,
    algorithm: AutomorphismAlgorithm,
    automorphisms: Option<&ProjectedAutomorphismOutput>,
    candidates: &mut [NodeId],
) -> bool {
    let backend_called = automorphisms.is_none();
    let labels = automorphisms.map_or_else(
        || {
            adapter
                .automorphisms_for_partition(partition, algorithm)
                .canonical_labels
        },
        |output| output.canonical_labels.clone(),
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
        &mut best,
        &mut stats,
    );

    CanonicalSearchResult {
        candidate: best.expect("every finite partition has a discrete entity refinement"),
        stats,
    }
}

#[allow(clippy::too_many_arguments)]
fn search_partition<K, LeafCandidate, PrefixWorse>(
    adapter: &AutomorphismAdapter,
    partition: OrderedPartition,
    algorithm: AutomorphismAlgorithm,
    options: CanonicalSearchOptions,
    leaf_candidate: &LeafCandidate,
    prefix_worse: &PrefixWorse,
    best: &mut Option<CanonicalCandidate<K>>,
    stats: &mut CanonicalSearchStats,
) where
    K: Ord,
    LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
    PrefixWorse: Fn(&OrderedPartition, &CanonicalCandidate<K>) -> bool,
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
    let automorphisms = options.automorphism_pruning.then(|| {
        stats.backend_calls += 1;
        adapter.automorphisms_for_partition(&partition, algorithm)
    });

    if options.automorphism_pruning {
        let orbits = &automorphisms
            .as_ref()
            .expect("automorphisms requested for orbit pruning")
            .orbits;
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
            best,
            stats,
        );
    }
}

/// Construct the normalized entity and incidence keys used for initial classes.
fn initial_class_keys(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
) -> Result<(Vec<InitialClassKey>, Vec<InitialClassKey>), Contradiction> {
    let entity_keys = incidence_graph
        .graph()
        .node_ids()
        .map(|node| entity_class_key(molecule, incidence_graph.entity(node)))
        .collect::<Result<Vec<_>, _>>()?;
    let incidence_keys = incidence_graph
        .incidences()
        .map(|(_, incidence)| incidence_key(incidence).map(InitialClassKey::Incidence))
        .collect::<Result<Vec<_>, _>>()?;

    Ok((entity_keys, incidence_keys))
}

/// Construct topology-level partition descriptors for the exact adapter.
fn partition_descriptors(
    adapter: &AutomorphismAdapter,
    entity_keys: &[InitialClassKey],
    incidence_keys: &[InitialClassKey],
) -> Vec<InitialClassKey> {
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
    entity_keys: &[InitialClassKey],
    incidence_graph: &IncidenceGraph,
) -> Vec<InitialClassKey> {
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
                        .expect("initial classes established incidence normalization"),
                };
                InitialClassKey::Incidence(value)
            }
        })
        .collect()
}

fn constitution_entity_classes(molecule: &Molecule) -> Result<Vec<u32>, Contradiction> {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Constitution);
    let (entity_keys, incidence_keys) = initial_class_keys(molecule, &incidence_graph)?;
    let classes = rank_initial_classes(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
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
        .filter(move |permutation| kind.class_key().space().reindex(0, *permutation).is_some())
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
    entity_keys: &[InitialClassKey],
    entity_classes: &[u32],
) -> Result<Vec<InitialClassKey>, Contradiction> {
    let entity_class = |entity: Entity| entity_classes[incidence_graph.node_of(entity).index()];

    adapter
        .node_sources
        .iter()
        .map(|source| match *source {
            SubdivisionNodeSource::Node(node) => {
                let entity = incidence_graph.entity(node);
                let InitialClassKey::Entity { position, .. } = entity_keys[node.index()] else {
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
                Ok(InitialClassKey::Entity { position, value })
            }
            SubdivisionNodeSource::Edge(edge) => {
                let value = match incidence_graph.incidence(edge) {
                    Incidence::DativeDonor | Incidence::DativeAcceptor => variant(1, []),
                    Incidence::AromaticParticipant(_) => variant(3, []),
                    Incidence::MulticenterParticipant(_) => variant(4, []),
                    incidence => incidence_key(incidence)?,
                };
                Ok(InitialClassKey::Incidence(value))
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
    entity_keys: &[InitialClassKey],
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

/// Assign dense ordered classes to normalized entity and incidence keys.
fn rank_initial_classes(
    entity_keys: &[InitialClassKey],
    incidence_keys: &[InitialClassKey],
) -> InitialClasses {
    let keys = entity_keys
        .iter()
        .chain(incidence_keys.iter())
        .cloned()
        .collect::<BTreeSet<_>>();
    let classes = keys
        .into_iter()
        .enumerate()
        .map(|(class, key)| (key, class as u32))
        .collect::<BTreeMap<_, _>>();

    InitialClasses {
        entities: entity_keys.iter().map(|key| classes[key]).collect(),
        incidences: incidence_keys.iter().map(|key| classes[key]).collect(),
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

fn entity_class_key(molecule: &Molecule, entity: Entity) -> Result<InitialClassKey, Contradiction> {
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

    Ok(InitialClassKey::Entity { position, value })
}

fn reaction_span_entity_class_key(
    span: &ReactionSpan,
    entity: Entity,
) -> Result<InitialClassKey, Contradiction> {
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
            normalized_entity_span_key(
                span.dative_bonds().data(RelationId::from(id)),
                |attributes| {
                    Ok(positioned_product([(
                        2,
                        num_form_key(attributes.order.normalized()?.as_ref()),
                    )]))
                },
            )?,
        ),
        Entity::AromaticSystem(id) => (
            EntityBlockPosition::AROMATIC_SYSTEM,
            normalized_entity_span_key(
                span.aromatic_systems().data(RelationId::from(id)),
                |attributes| {
                    Ok(positioned_product([
                        (2, num_form_key(attributes.charge.normalized()?.as_ref())),
                        (
                            3,
                            unpaired_electrons_form_key(
                                attributes.unpaired_electrons.normalized()?.as_ref(),
                            ),
                        ),
                    ]))
                },
            )?,
        ),
        Entity::MulticenterBond(id) => (
            EntityBlockPosition::MULTICENTER_BOND,
            normalized_entity_span_key(
                span.multicenter_bonds().data(RelationId::from(id)),
                |attributes| {
                    Ok(positioned_product([
                        (2, num_form_key(attributes.charge.normalized()?.as_ref())),
                        (
                            3,
                            unpaired_electrons_form_key(
                                attributes.unpaired_electrons.normalized()?.as_ref(),
                            ),
                        ),
                    ]))
                },
            )?,
        ),
        Entity::NoncovalentBond(id) => (
            EntityBlockPosition::NONCOVALENT_BOND,
            normalized_entity_span_key(
                span.noncovalent_bonds().data(RelationId::from(id)),
                |attributes| {
                    Ok(positioned_product([(
                        1,
                        noncovalent_bond_kind_form_key(attributes.kind.normalized()?.as_ref()),
                    )]))
                },
            )?,
        ),
        Entity::StereoAtom(id) => (
            EntityBlockPosition::STEREO_ATOM,
            normalized_entity_span_key(
                span.stereo_atoms().data(RelationId::from(id)),
                |attributes| {
                    Ok(positioned_product([(
                        2,
                        option(attributes.configuration.kind().map(stereo_kind_key)),
                    )]))
                },
            )?,
        ),
        Entity::StereoBond(id) => (
            EntityBlockPosition::STEREO_BOND,
            normalized_entity_span_key(
                span.stereo_bonds().data(RelationId::from(id)),
                |attributes| {
                    Ok(positioned_product([(
                        2,
                        option(attributes.configuration.kind().map(stereo_kind_key)),
                    )]))
                },
            )?,
        ),
    };
    Ok(InitialClassKey::Entity { position, value })
}

fn topology_candidate(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    order: &[NodeId],
) -> Result<CanonicalCandidate<CanonicalComparisonKey>, Contradiction> {
    let atom_count = incidence_graph.entity_count(EntityKind::Atom);
    let bond_count = incidence_graph.entity_count(EntityKind::Bond);
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
    permutations: Vec<Permutation>,
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
    let mut permutations = Vec::new();
    for permutation in stereo_frame_permutations(kind) {
        let candidate = (
            permutation.act(ligands),
            configuration.apply(permutation).normalize()?,
        );
        match minimum.as_ref().map(|value| candidate.cmp(value)) {
            None | Some(Ordering::Less) => {
                minimum = Some(candidate);
                permutations.clear();
                permutations.push(permutation);
            }
            Some(Ordering::Equal) => permutations.push(permutation),
            Some(Ordering::Greater) => {}
        }
    }

    Ok(
        minimum.map(|(ligands, configuration)| CanonicalStereoFrame {
            ligands,
            configuration,
            permutations,
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
            let ligands = remap_ligands(ligands);
            match configuration {
                StereoConfigurationForm::Undetermined => {
                    let (ligands, _) = sort_ligand_frame(&ligands);
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
    let remapped = molecule.remap(&correspondence);
    let (constraints, complete) = canonicalize_complete_stereo_frames(remapped)?;
    candidate.key.constraints = constraints;
    Ok((candidate, complete))
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
        let (family, source) = match incidence_graph.entity(node) {
            Entity::Atom(id) => (0, id.index()),
            Entity::Bond(id) => (1, id.index()),
            Entity::DativeBond(id) => (2, id.index()),
            Entity::AromaticSystem(id) => (3, id.index()),
            Entity::MulticenterBond(id) => (4, id.index()),
            Entity::NoncovalentBond(id) => (5, id.index()),
            Entity::StereoAtom(id) => (6, id.index()),
            Entity::StereoBond(id) => (7, id.index()),
        };
        images[family][source] = next[family];
        next[family] += 1;
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
        let (family, source, lhs_present) = match entity {
            Entity::Atom(id) => (0, id.index(), span.atoms()[id.index()].lhs().is_some()),
            Entity::Bond(id) => (1, id.index(), span.bonds()[id.index()].lhs().is_some()),
            Entity::DativeBond(id) => (
                2,
                id.index(),
                span.dative_bonds()
                    .data(RelationId::from(id))
                    .lhs()
                    .is_some(),
            ),
            Entity::AromaticSystem(id) => (
                3,
                id.index(),
                span.aromatic_systems()
                    .data(RelationId::from(id))
                    .lhs()
                    .is_some(),
            ),
            Entity::MulticenterBond(id) => (
                4,
                id.index(),
                span.multicenter_bonds()
                    .data(RelationId::from(id))
                    .lhs()
                    .is_some(),
            ),
            Entity::NoncovalentBond(id) => (
                5,
                id.index(),
                span.noncovalent_bonds()
                    .data(RelationId::from(id))
                    .lhs()
                    .is_some(),
            ),
            Entity::StereoAtom(id) => (
                6,
                id.index(),
                span.stereo_atoms()
                    .data(RelationId::from(id))
                    .lhs()
                    .is_some(),
            ),
            Entity::StereoBond(id) => (
                7,
                id.index(),
                span.stereo_bonds()
                    .data(RelationId::from(id))
                    .lhs()
                    .is_some(),
            ),
        };
        if lhs_present {
            lhs[family].push(source);
        } else {
            added[family].push(source);
        }
    }
    for family in 0..counts.len() {
        let mut present = lhs[family]
            .iter()
            .chain(added[family].iter())
            .copied()
            .collect::<BTreeSet<_>>();
        for source in 0..counts[family] {
            if present.insert(source) {
                if reaction_span_lhs_present(span, family, source) {
                    lhs[family].push(source);
                } else {
                    added[family].push(source);
                }
            }
        }
    }
    let mut images = counts.map(|count| vec![0; count]);
    for family in 0..images.len() {
        for (image, source) in lhs[family]
            .iter()
            .chain(added[family].iter())
            .copied()
            .enumerate()
        {
            images[family][source] = image;
        }
    }
    molecule_correspondence(&images)
}

fn reaction_span_lhs_present(span: &ReactionSpan, family: usize, source: usize) -> bool {
    match family {
        0 => span.atoms()[source].lhs().is_some(),
        1 => span.bonds()[source].lhs().is_some(),
        2 => span
            .dative_bonds()
            .data(RelationId(source as u32))
            .lhs()
            .is_some(),
        3 => span
            .aromatic_systems()
            .data(RelationId(source as u32))
            .lhs()
            .is_some(),
        4 => span
            .multicenter_bonds()
            .data(RelationId(source as u32))
            .lhs()
            .is_some(),
        5 => span
            .noncovalent_bonds()
            .data(RelationId(source as u32))
            .lhs()
            .is_some(),
        6 => span
            .stereo_atoms()
            .data(RelationId(source as u32))
            .lhs()
            .is_some(),
        7 => span
            .stereo_bonds()
            .data(RelationId(source as u32))
            .lhs()
            .is_some(),
        _ => unreachable!("reaction-span entity family index is in range"),
    }
}

fn apply_position_order<T: Clone>(values: &[T], order: &[ParticipantPosition]) -> Option<Vec<T>> {
    if values.len() != order.len() {
        return None;
    }
    let mut seen = vec![false; order.len()];
    let mut reordered = Vec::with_capacity(order.len());
    for position in order {
        let index = position.index();
        if index >= order.len() || seen[index] {
            return None;
        }
        seen[index] = true;
        reordered.push(values[index].clone());
    }
    Some(reordered)
}

fn inverse_position_order(order: &[ParticipantPosition]) -> Option<Vec<ParticipantPosition>> {
    let mut inverse = vec![ParticipantPosition(0); order.len()];
    let mut seen = vec![false; order.len()];
    for (new, old) in order.iter().enumerate() {
        let old = old.index();
        if old >= order.len() || seen[old] {
            return None;
        }
        seen[old] = true;
        inverse[old] = ParticipantPosition(new as u32);
    }
    Some(inverse)
}

fn permutation_from_position_order(order: &[ParticipantPosition]) -> Option<Permutation> {
    let image = order
        .iter()
        .map(|position| position.index())
        .collect::<Vec<_>>();
    Permutation::try_from(image.as_slice()).ok()
}

fn position_order_from_permutation(permutation: Permutation) -> Vec<ParticipantPosition> {
    (0..permutation.degree())
        .map(|position| ParticipantPosition(permutation.apply(position) as u32))
        .collect()
}

/// Sort a structural ligand multiset and return the position order carrying the original frame
/// into the sorted frame. Equal occurrences retain their input order; their exchange remains a
/// structural automorphism rather than becoming an arbitrary tie-break here.
fn sort_ligand_frame(ligands: &[StereoLigand]) -> (Vec<StereoLigand>, Vec<ParticipantPosition>) {
    let mut sorted = ligands.to_vec();
    let order = Unordered::canonicalize_positions(&mut sorted);
    (sorted, order)
}

fn reframe_stereo_atom_constraint_by_order(
    constraint: StereoAtomConstraintForm,
    order: &[ParticipantPosition],
) -> Option<StereoAtomConstraintForm> {
    if let Some(permutation) = permutation_from_position_order(order) {
        return StereoAtomForm {
            configuration: StereoConfigurationForm::Undetermined,
            constraints: constraint.into(),
        }
        .transform_frame_by(permutation)?
        .constraints
        .into_iter()
        .next();
    }
    let inverse = inverse_position_order(order)?;
    Some(match constraint {
        StereoAtomConstraintForm::Topicity(topicity) => {
            StereoAtomConstraintForm::Topicity(TopicityForm {
                pair: StereoLigandPair::new(
                    inverse.get(topicity.pair.first().index())?.index().into(),
                    inverse.get(topicity.pair.second().index())?.index().into(),
                ),
                relation: topicity.relation,
            })
        }
        StereoAtomConstraintForm::Stereogenicity(stereogenicity) => {
            StereoAtomConstraintForm::Stereogenicity(stereogenicity)
        }
        StereoAtomConstraintForm::LigandSymmetry(_) | StereoAtomConstraintForm::Fluxionality(_) => {
            return None
        }
    })
}

fn reframe_stereo_bond_constraint_by_order(
    constraint: StereoBondConstraintForm,
    order: &[ParticipantPosition],
) -> Option<StereoBondConstraintForm> {
    if let Some(permutation) = permutation_from_position_order(order) {
        return StereoBondForm {
            configuration: StereoConfigurationForm::Undetermined,
            constraints: constraint.into(),
        }
        .transform_frame_by(permutation)?
        .constraints
        .into_iter()
        .next();
    }
    let inverse = inverse_position_order(order)?;
    Some(match constraint {
        StereoBondConstraintForm::Topicity(topicity) => {
            StereoBondConstraintForm::Topicity(TopicityForm {
                pair: StereoLigandPair::new(
                    inverse.get(topicity.pair.first().index())?.index().into(),
                    inverse.get(topicity.pair.second().index())?.index().into(),
                ),
                relation: topicity.relation,
            })
        }
        StereoBondConstraintForm::Stereogenicity(stereogenicity) => {
            StereoBondConstraintForm::Stereogenicity(stereogenicity)
        }
        StereoBondConstraintForm::LigandSymmetry(_) | StereoBondConstraintForm::Fluxionality(_) => {
            return None
        }
    })
}

fn reframe_stereo_atom_form_by_order(
    form: &StereoAtomForm,
    order: &[ParticipantPosition],
) -> Option<StereoAtomForm> {
    if let Some(permutation) = permutation_from_position_order(order) {
        return form.transform_frame_by(permutation);
    }
    if form.configuration != StereoConfigurationForm::Undetermined {
        return None;
    }
    Some(StereoAtomForm {
        configuration: StereoConfigurationForm::Undetermined,
        constraints: form
            .constraints
            .iter()
            .cloned()
            .map(|constraint| reframe_stereo_atom_constraint_by_order(constraint, order))
            .collect::<Option<StereoAtomConstraintsForm>>()?,
    })
}

fn reframe_stereo_bond_form_by_order(
    form: &StereoBondForm,
    order: &[ParticipantPosition],
) -> Option<StereoBondForm> {
    if let Some(permutation) = permutation_from_position_order(order) {
        return form.transform_frame_by(permutation);
    }
    if form.configuration != StereoConfigurationForm::Undetermined {
        return None;
    }
    Some(StereoBondForm {
        configuration: StereoConfigurationForm::Undetermined,
        constraints: form
            .constraints
            .iter()
            .cloned()
            .map(|constraint| reframe_stereo_bond_constraint_by_order(constraint, order))
            .collect::<Option<StereoBondConstraintsForm>>()?,
    })
}

fn reframe_stereo_atom_span_by_order(
    span: &EntitySpan<StereoAtomForm>,
    order: &[ParticipantPosition],
) -> Option<EntitySpan<StereoAtomForm>> {
    Some(match span {
        EntitySpan::Unchanged(value) => {
            EntitySpan::Unchanged(reframe_stereo_atom_form_by_order(value, order)?)
        }
        EntitySpan::Added(value) => {
            EntitySpan::Added(reframe_stereo_atom_form_by_order(value, order)?)
        }
        EntitySpan::Removed(value) => {
            EntitySpan::Removed(reframe_stereo_atom_form_by_order(value, order)?)
        }
        EntitySpan::Modified { lhs, rhs } => EntitySpan::Modified {
            lhs: reframe_stereo_atom_form_by_order(lhs, order)?,
            rhs: reframe_stereo_atom_form_by_order(rhs, order)?,
        },
    })
}

fn reframe_stereo_bond_span_by_order(
    span: &EntitySpan<StereoBondForm>,
    order: &[ParticipantPosition],
) -> Option<EntitySpan<StereoBondForm>> {
    Some(match span {
        EntitySpan::Unchanged(value) => {
            EntitySpan::Unchanged(reframe_stereo_bond_form_by_order(value, order)?)
        }
        EntitySpan::Added(value) => {
            EntitySpan::Added(reframe_stereo_bond_form_by_order(value, order)?)
        }
        EntitySpan::Removed(value) => {
            EntitySpan::Removed(reframe_stereo_bond_form_by_order(value, order)?)
        }
        EntitySpan::Modified { lhs, rhs } => EntitySpan::Modified {
            lhs: reframe_stereo_bond_form_by_order(lhs, order)?,
            rhs: reframe_stereo_bond_form_by_order(rhs, order)?,
        },
    })
}

fn reframe_molecule_constraint_by_order(
    constraint: Constraint,
    entity: Entity,
    order: &[ParticipantPosition],
) -> Option<Constraint> {
    Some(match constraint {
        Constraint::StereoAtom(id, kind, constraint) if entity == Entity::StereoAtom(id) => {
            Constraint::StereoAtom(
                id,
                kind,
                reframe_stereo_atom_constraint_by_order(constraint, order)?,
            )
        }
        Constraint::StereoBond(id, kind, constraint) if entity == Entity::StereoBond(id) => {
            Constraint::StereoBond(
                id,
                kind,
                reframe_stereo_bond_constraint_by_order(constraint, order)?,
            )
        }
        Constraint::And(constraints) => Constraint::And(
            constraints
                .into_iter()
                .map(|constraint| reframe_molecule_constraint_by_order(constraint, entity, order))
                .collect::<Option<Vec<_>>>()?,
        ),
        Constraint::Or(constraints) => Constraint::Or(
            constraints
                .into_iter()
                .map(|constraint| reframe_molecule_constraint_by_order(constraint, entity, order))
                .collect::<Option<Vec<_>>>()?,
        ),
        Constraint::Not(constraint) => Constraint::Not(Box::new(
            reframe_molecule_constraint_by_order(*constraint, entity, order)?,
        )),
        constraint => constraint,
    })
}

fn molecule_entries(molecule: &Molecule) -> MoleculeEntries {
    MoleculeEntries {
        atoms: molecule
            .atoms()
            .iter()
            .map(|atom| atom.attributes.clone())
            .collect(),
        bonds: molecule
            .bonds()
            .iter()
            .map(|bond| {
                let [first, second] = bond.atom_ids();
                (first, second, bond.attributes.clone())
            })
            .collect(),
        dative: molecule
            .dative_bonds()
            .iter()
            .map(|bond| {
                (
                    bond.donor_ids().collect(),
                    bond.acceptor_id(),
                    bond.attributes.clone(),
                )
            })
            .collect(),
        aromatic: molecule
            .aromatic_systems()
            .iter()
            .map(|system| (system.atom_ids().collect(), system.attributes.clone()))
            .collect(),
        multicenter: molecule
            .multicenter_bonds()
            .iter()
            .map(|bond| (bond.atom_ids().collect(), bond.attributes.clone()))
            .collect(),
        noncovalent: molecule
            .noncovalent_bonds()
            .iter()
            .map(|bond| {
                let [first, second] = bond.atom_ids();
                (first, second, bond.attributes.clone())
            })
            .collect(),
        stereo_atoms: molecule
            .stereo_atoms()
            .iter()
            .map(|stereo| {
                (
                    stereo.site_id(),
                    stereo.ligand_frame(),
                    stereo.attributes.clone(),
                )
            })
            .collect(),
        stereo_bonds: molecule
            .stereo_bonds()
            .iter()
            .map(|stereo| {
                (
                    stereo.site_id(),
                    stereo.ligand_frame(),
                    stereo.attributes.clone(),
                )
            })
            .collect(),
        constraints: molecule.constraints().clone(),
    }
}

fn reframe_reaction_span_stereo_atom_by_order(
    span: &ReactionSpan,
    id: StereoAtomId,
    order: &[ParticipantPosition],
) -> Option<ReactionSpan> {
    let mut entries = span.entries();
    let (_, ligands, attributes) = &mut entries.stereo_atoms[id.index()];
    *attributes = reframe_stereo_atom_span_by_order(attributes, order)?;
    *ligands = apply_position_order(ligands, order)?;
    entries.constraints = entries
        .constraints
        .into_iter()
        .map(|span| {
            let reframe = |constraint| {
                reframe_molecule_constraint_by_order(constraint, Entity::StereoAtom(id), order)
            };
            Some(match span {
                ConstraintSpan::Unchanged(value) => ConstraintSpan::Unchanged(reframe(value)?),
                ConstraintSpan::Added(value) => ConstraintSpan::Added(reframe(value)?),
                ConstraintSpan::Removed(value) => ConstraintSpan::Removed(reframe(value)?),
            })
        })
        .collect::<Option<Vec<_>>>()?;
    ReactionSpan::try_from_entries(entries).ok()
}

fn reframe_reaction_span_stereo_bond_by_order(
    span: &ReactionSpan,
    id: StereoBondId,
    order: &[ParticipantPosition],
) -> Option<ReactionSpan> {
    let mut entries = span.entries();
    let (_, ligands, attributes) = &mut entries.stereo_bonds[id.index()];
    *attributes = reframe_stereo_bond_span_by_order(attributes, order)?;
    *ligands = apply_position_order(ligands, order)?;
    entries.constraints = entries
        .constraints
        .into_iter()
        .map(|span| {
            let reframe = |constraint| {
                reframe_molecule_constraint_by_order(constraint, Entity::StereoBond(id), order)
            };
            Some(match span {
                ConstraintSpan::Unchanged(value) => ConstraintSpan::Unchanged(reframe(value)?),
                ConstraintSpan::Added(value) => ConstraintSpan::Added(reframe(value)?),
                ConstraintSpan::Removed(value) => ConstraintSpan::Removed(reframe(value)?),
            })
        })
        .collect::<Option<Vec<_>>>()?;
    ReactionSpan::try_from_entries(entries).ok()
}

fn canonicalize_reaction_span_stereo_frames(
    mut span: ReactionSpan,
) -> Result<ReactionSpan, Contradiction> {
    for index in 0..span.stereo_atoms().count() {
        let id = StereoAtomId(index as u32);
        let relation = RelationId(index as u32);
        let ligands = span.stereo_atoms().participants_2(relation).to_vec();
        let mut best = None;
        for permutation in all_span_frame_permutations(
            &ligands,
            span.stereo_atoms()
                .data(relation)
                .lhs()
                .map(|form| &form.configuration),
            span.stereo_atoms()
                .data(relation)
                .rhs()
                .map(|form| &form.configuration),
        )? {
            let order = position_order_from_permutation(permutation);
            let candidate = reframe_reaction_span_stereo_atom_by_order(&span, id, &order)
                .expect("integrity established a valid stereo-atom frame action");
            let key = stereo_atom_span_frame_key(&candidate, id)?;
            if best.as_ref().is_none_or(|(best_key, _)| key < *best_key) {
                best = Some((key, candidate));
            }
        }
        span = best
            .expect("every valid stereo frame has an identity action")
            .1;
    }
    for index in 0..span.stereo_bonds().count() {
        let id = StereoBondId(index as u32);
        let relation = RelationId(index as u32);
        let ligands = span.stereo_bonds().participants_2(relation).to_vec();
        let mut best = None;
        for permutation in all_span_frame_permutations(
            &ligands,
            span.stereo_bonds()
                .data(relation)
                .lhs()
                .map(|form| &form.configuration),
            span.stereo_bonds()
                .data(relation)
                .rhs()
                .map(|form| &form.configuration),
        )? {
            let order = position_order_from_permutation(permutation);
            let candidate = reframe_reaction_span_stereo_bond_by_order(&span, id, &order)
                .expect("integrity established a valid stereo-bond frame action");
            let key = stereo_bond_span_frame_key(&candidate, id)?;
            if best.as_ref().is_none_or(|(best_key, _)| key < *best_key) {
                best = Some((key, candidate));
            }
        }
        span = best
            .expect("every valid stereo frame has an identity action")
            .1;
    }
    Ok(span)
}

fn all_span_frame_permutations(
    ligands: &[StereoLigand],
    lhs: Option<&StereoConfigurationForm>,
    rhs: Option<&StereoConfigurationForm>,
) -> Result<Vec<Permutation>, Contradiction> {
    let mut permutations = stereo_frame_permutations_for(lhs, ligands.len())?;
    permutations
        .retain(|permutation| stereo_frame_permutation_allowed(rhs, ligands.len(), *permutation));
    Ok(permutations)
}

fn stereo_frame_permutations_for(
    configuration: Option<&StereoConfigurationForm>,
    degree: usize,
) -> Result<Vec<Permutation>, Contradiction> {
    let count = (1..=degree).product::<usize>().max(1);
    Ok((0..count)
        .map(|rank| Permutation::unrank(degree, rank))
        .filter(|&permutation| stereo_frame_permutation_allowed(configuration, degree, permutation))
        .collect())
}

fn stereo_frame_permutation_allowed(
    configuration: Option<&StereoConfigurationForm>,
    degree: usize,
    permutation: Permutation,
) -> bool {
    match configuration {
        None | Some(StereoConfigurationForm::Undetermined) => permutation.degree() == degree,
        Some(StereoConfigurationForm::Kinded(kind, _)) => {
            kind.degree() == degree && kind.class_key().space().reindex(0, permutation).is_some()
        }
    }
}

fn stereo_atom_span_frame_key(
    span: &ReactionSpan,
    id: StereoAtomId,
) -> Result<CanonicalKeyValue, Contradiction> {
    let relation = RelationId::from(id);
    Ok(product([
        sequence(
            span.stereo_atoms()
                .participants_2(relation)
                .iter()
                .map(|ligand| stereo_ligand_key(ligand.atom_id.0, ligand.kind)),
        ),
        normalized_entity_span_key(span.stereo_atoms().data(relation), |form| {
            Ok(stereo_configuration_form_key(
                form.configuration.normalized()?.as_ref(),
            ))
        })?,
    ]))
}

fn stereo_bond_span_frame_key(
    span: &ReactionSpan,
    id: StereoBondId,
) -> Result<CanonicalKeyValue, Contradiction> {
    let relation = RelationId::from(id);
    Ok(product([
        sequence(
            span.stereo_bonds()
                .participants_2(relation)
                .iter()
                .map(|ligand| stereo_ligand_key(ligand.atom_id.0, ligand.kind)),
        ),
        normalized_entity_span_key(span.stereo_bonds().data(relation), |form| {
            Ok(stereo_configuration_form_key(
                form.configuration.normalized()?.as_ref(),
            ))
        })?,
    ]))
}

fn reframe_stereo_atom(molecule: &Molecule, id: StereoAtomId, frame: Permutation) -> Molecule {
    reframe_stereo_atom_by_order(molecule, id, &position_order_from_permutation(frame))
        .expect("canonicalization generates a valid stereo-atom frame action")
}

fn reframe_stereo_atom_by_order(
    molecule: &Molecule,
    id: StereoAtomId,
    order: &[ParticipantPosition],
) -> Option<Molecule> {
    let mut entries = molecule_entries(molecule);
    let (_, ligands, attributes) = &mut entries.stereo_atoms[id.index()];
    *attributes = reframe_stereo_atom_form_by_order(attributes, order)?;
    *ligands = apply_position_order(ligands, order)?;
    entries.constraints = entries
        .constraints
        .into_iter()
        .map(|constraint| {
            reframe_molecule_constraint_by_order(constraint, Entity::StereoAtom(id), order)
        })
        .collect::<Option<Vec<_>>>()?
        .into();
    Some(Molecule::from_entries(entries))
}

fn reframe_stereo_bond(molecule: &Molecule, id: StereoBondId, frame: Permutation) -> Molecule {
    reframe_stereo_bond_by_order(molecule, id, &position_order_from_permutation(frame))
        .expect("canonicalization generates a valid stereo-bond frame action")
}

fn reframe_stereo_bond_by_order(
    molecule: &Molecule,
    id: StereoBondId,
    order: &[ParticipantPosition],
) -> Option<Molecule> {
    let mut entries = molecule_entries(molecule);
    let (_, ligands, attributes) = &mut entries.stereo_bonds[id.index()];
    *attributes = reframe_stereo_bond_form_by_order(attributes, order)?;
    *ligands = apply_position_order(ligands, order)?;
    entries.constraints = entries
        .constraints
        .into_iter()
        .map(|constraint| {
            reframe_molecule_constraint_by_order(constraint, Entity::StereoBond(id), order)
        })
        .collect::<Option<Vec<_>>>()?
        .into();
    Some(Molecule::from_entries(entries))
}

fn canonicalize_stereo_frames(mut molecule: Molecule) -> Result<Molecule, Contradiction> {
    let stereo_atom_count = molecule.stereo_atoms().count();
    for index in 0..stereo_atom_count {
        let id = StereoAtomId(index as u32);
        let (ligands, configuration) = {
            let stereo = molecule
                .stereo_atoms()
                .get(id)
                .expect("dense stereo-atom id is in range");
            (
                stereo.ligand_frame(),
                stereo.attributes.configuration.clone(),
            )
        };
        molecule = match configuration {
            StereoConfigurationForm::Undetermined => {
                let (_, order) = sort_ligand_frame(&ligands);
                reframe_stereo_atom_by_order(&molecule, id, &order)
                    .expect("integrity established a valid kindless stereo-atom frame")
            }
            StereoConfigurationForm::Kinded(..) => {
                let frame = canonical_kinded_stereo_frame(&ligands, &configuration)?
                    .expect("integrity established the kinded stereo-atom frame degree");
                reframe_stereo_atom(&molecule, id, frame.permutations[0])
            }
        };
    }

    let stereo_bond_count = molecule.stereo_bonds().count();
    for index in 0..stereo_bond_count {
        let id = StereoBondId(index as u32);
        let (ligands, configuration) = {
            let stereo = molecule
                .stereo_bonds()
                .get(id)
                .expect("dense stereo-bond id is in range");
            (
                stereo.ligand_frame(),
                stereo.attributes.configuration.clone(),
            )
        };
        molecule = match configuration {
            StereoConfigurationForm::Undetermined => {
                let (_, order) = sort_ligand_frame(&ligands);
                reframe_stereo_bond_by_order(&molecule, id, &order)
                    .expect("integrity established a valid kindless stereo-bond frame")
            }
            StereoConfigurationForm::Kinded(..) => {
                let frame = canonical_kinded_stereo_frame(&ligands, &configuration)?
                    .expect("integrity established the kinded stereo-bond frame degree");
                reframe_stereo_bond(&molecule, id, frame.permutations[0])
            }
        };
    }

    Ok(molecule)
}

fn canonicalize_complete_stereo_frames(
    molecule: Molecule,
) -> Result<(Vec<ConstraintBlockKey>, Molecule), Contradiction> {
    let mut candidates = vec![molecule];
    let stereo_atom_count = candidates[0].stereo_atoms().count();
    for index in 0..stereo_atom_count {
        let id = StereoAtomId(index as u32);
        let mut next = Vec::new();
        for molecule in candidates {
            let stereo = molecule
                .stereo_atoms()
                .get(id)
                .expect("dense stereo-atom id is in range");
            let ligands = stereo.ligand_frame();
            let configuration = stereo.attributes.configuration.clone();
            let permutations = match configuration {
                StereoConfigurationForm::Undetermined => {
                    let (sorted, _) = sort_ligand_frame(&ligands);
                    Permutation::between_all(&ligands, &sorted)
                }
                StereoConfigurationForm::Kinded(..) => {
                    canonical_kinded_stereo_frame(&ligands, &configuration)?
                        .expect("integrity established the kinded stereo-atom frame degree")
                        .permutations
                }
            };
            next.extend(
                permutations
                    .into_iter()
                    .map(|frame| reframe_stereo_atom(&molecule, id, frame)),
            );
        }
        candidates = next;
    }

    let stereo_bond_count = candidates[0].stereo_bonds().count();
    for index in 0..stereo_bond_count {
        let id = StereoBondId(index as u32);
        let mut next = Vec::new();
        for molecule in candidates {
            let stereo = molecule
                .stereo_bonds()
                .get(id)
                .expect("dense stereo-bond id is in range");
            let ligands = stereo.ligand_frame();
            let configuration = stereo.attributes.configuration.clone();
            let permutations = match configuration {
                StereoConfigurationForm::Undetermined => {
                    let (sorted, _) = sort_ligand_frame(&ligands);
                    Permutation::between_all(&ligands, &sorted)
                }
                StereoConfigurationForm::Kinded(..) => {
                    canonical_kinded_stereo_frame(&ligands, &configuration)?
                        .expect("integrity established the kinded stereo-bond frame degree")
                        .permutations
                }
            };
            next.extend(
                permutations
                    .into_iter()
                    .map(|frame| reframe_stereo_bond(&molecule, id, frame)),
            );
        }
        candidates = next;
    }

    candidates
        .into_iter()
        .map(|molecule| {
            let molecule = normalize_molecule(molecule)?;
            Ok((constraint_blocks(&molecule), molecule))
        })
        .collect::<Result<Vec<_>, Contradiction>>()?
        .into_iter()
        .min_by(|lhs, rhs| lhs.0.cmp(&rhs.0))
        .ok_or(Contradiction)
}

fn normalize_molecule(mut molecule: Molecule) -> Result<Molecule, Contradiction> {
    let mut atoms = molecule
        .atoms()
        .iter()
        .map(|atom| atom.attributes.clone().normalize())
        .collect::<Result<Vec<_>, _>>()?
        .into_iter();
    molecule.modify_atoms(|_| atoms.next().expect("one normalized form per atom"));

    let mut bonds = molecule
        .bonds()
        .iter()
        .map(|bond| bond.attributes.clone().normalize())
        .collect::<Result<Vec<_>, _>>()?
        .into_iter();
    molecule.modify_bonds(|_| bonds.next().expect("one normalized form per bond"));

    let mut dative_bonds = molecule
        .dative_bonds()
        .iter()
        .map(|bond| bond.attributes.clone().normalize())
        .collect::<Result<Vec<_>, _>>()?
        .into_iter();
    molecule.modify_dative_bonds(|_| {
        dative_bonds
            .next()
            .expect("one normalized form per dative bond")
    });

    let mut aromatic_systems = molecule
        .aromatic_systems()
        .iter()
        .map(|system| system.attributes.clone().normalize())
        .collect::<Result<Vec<_>, _>>()?
        .into_iter();
    molecule.modify_aromatic_systems(|_| {
        aromatic_systems
            .next()
            .expect("one normalized form per aromatic system")
    });

    let mut multicenter_bonds = molecule
        .multicenter_bonds()
        .iter()
        .map(|bond| bond.attributes.clone().normalize())
        .collect::<Result<Vec<_>, _>>()?
        .into_iter();
    molecule.modify_multicenter_bonds(|_| {
        multicenter_bonds
            .next()
            .expect("one normalized form per multicenter bond")
    });

    let mut noncovalent_bonds = molecule
        .noncovalent_bonds()
        .iter()
        .map(|bond| bond.attributes.clone().normalize())
        .collect::<Result<Vec<_>, _>>()?
        .into_iter();
    molecule.modify_noncovalent_bonds(|_| {
        noncovalent_bonds
            .next()
            .expect("one normalized form per noncovalent bond")
    });

    let mut stereo_atoms = molecule
        .stereo_atoms()
        .iter()
        .map(|stereo| stereo.attributes.clone().normalize())
        .collect::<Result<Vec<_>, _>>()?
        .into_iter();
    molecule.modify_stereo_atoms(|_| {
        stereo_atoms
            .next()
            .expect("one normalized form per stereo atom")
    });

    let mut stereo_bonds = molecule
        .stereo_bonds()
        .iter()
        .map(|stereo| stereo.attributes.clone().normalize())
        .collect::<Result<Vec<_>, _>>()?
        .into_iter();
    molecule.modify_stereo_bonds(|_| {
        stereo_bonds
            .next()
            .expect("one normalized form per stereo bond")
    });

    let constraints = molecule.constraints().clone().normalize()?;
    *molecule.constraints_mut() = constraints;
    Ok(molecule)
}

fn normalize_entity_span<T: Normalize>(
    span: EntitySpan<T>,
) -> Result<EntitySpan<T>, Contradiction> {
    Ok(match span {
        EntitySpan::Unchanged(value) => EntitySpan::Unchanged(value.normalize()?),
        EntitySpan::Added(value) => EntitySpan::Added(value.normalize()?),
        EntitySpan::Removed(value) => EntitySpan::Removed(value.normalize()?),
        EntitySpan::Modified { lhs, rhs } => EntitySpan::Modified {
            lhs: lhs.normalize()?,
            rhs: rhs.normalize()?,
        },
    })
}

fn normalize_constraint_span(span: ConstraintSpan) -> Result<ConstraintSpan, Contradiction> {
    Ok(match span {
        ConstraintSpan::Unchanged(value) => ConstraintSpan::Unchanged(value.normalize()?),
        ConstraintSpan::Added(value) => ConstraintSpan::Added(value.normalize()?),
        ConstraintSpan::Removed(value) => ConstraintSpan::Removed(value.normalize()?),
    })
}

fn normalize_reaction_span(span: ReactionSpan) -> Result<ReactionSpan, Contradiction> {
    let mut entries = span.entries();
    entries.atoms = entries
        .atoms
        .into_iter()
        .map(normalize_entity_span)
        .collect::<Result<_, _>>()?;
    for (_, _, value) in &mut entries.bonds {
        *value = normalize_entity_span(value.clone())?;
    }
    for (_, _, value) in &mut entries.dative {
        *value = normalize_entity_span(value.clone())?;
    }
    for (_, value) in &mut entries.aromatic {
        *value = normalize_entity_span(value.clone())?;
    }
    for (_, value) in &mut entries.multicenter {
        *value = normalize_entity_span(value.clone())?;
    }
    for (_, _, value) in &mut entries.noncovalent {
        *value = normalize_entity_span(value.clone())?;
    }
    for (_, _, value) in &mut entries.stereo_atoms {
        *value = normalize_entity_span(value.clone())?;
    }
    for (_, _, value) in &mut entries.stereo_bonds {
        *value = normalize_entity_span(value.clone())?;
    }
    entries.constraints = entries
        .constraints
        .into_iter()
        .map(normalize_constraint_span)
        .collect::<Result<_, _>>()?;
    entries.constraints.sort();
    entries.constraints.dedup();
    Ok(ReactionSpan::from_entries(entries))
}

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
    molecule.check_integrity()?;
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Topology);
    let (entity_keys, incidence_keys) = initial_class_keys(molecule, &incidence_graph)?;
    let classes = rank_initial_classes(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
    let descriptors = partition_descriptors(&adapter, &entity_keys, &incidence_keys);
    let leaf_candidate = |order: &[NodeId]| {
        topology_candidate(molecule, &incidence_graph, order)
            .expect("initial classes established topology normalization")
    };
    let no_prefix = |_: &OrderedPartition, _: &CanonicalCandidate<_>| false;
    let selected = canonical_search(
        &adapter,
        &descriptors,
        context.automorphism_algorithm,
        options,
        &leaf_candidate,
        &no_prefix,
    );
    let correspondence =
        correspondence_from_order(molecule, &incidence_graph, &selected.candidate.entity_order);

    let canonical = normalize_molecule(molecule.remap(&correspondence))?;
    Ok((canonical, correspondence))
}

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
    molecule.check_integrity()?;
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Constitution);
    let (entity_keys, incidence_keys) = initial_class_keys(molecule, &incidence_graph)?;
    let classes = rank_initial_classes(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
    let descriptors = constitution_partition_descriptors(&adapter, &entity_keys, &incidence_graph);
    let leaf_candidate = |order: &[NodeId]| {
        constitution_candidate(molecule, &incidence_graph, order)
            .expect("initial classes established constitution normalization")
    };
    let no_prefix = |_: &OrderedPartition, _: &CanonicalCandidate<_>| false;
    let selected = canonical_search(
        &adapter,
        &descriptors,
        context.automorphism_algorithm,
        options,
        &leaf_candidate,
        &no_prefix,
    );
    let correspondence =
        correspondence_from_order(molecule, &incidence_graph, &selected.candidate.entity_order);

    let canonical = normalize_molecule(molecule.remap(&correspondence))?;
    Ok((canonical, correspondence))
}

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
    mut options: CanonicalSearchOptions,
) -> Result<(Molecule, MoleculeCorrespondence), MoleculeCanonicalizeError> {
    molecule.check_integrity()?;
    // A structure-frame automorphism acts on both entity ids and stereo configurations. The graph
    // adapter currently projects only the id action, so its orbits cannot soundly discard a branch
    // whose coupled frame action changes the leaf key. Keep full stereo search exhaustive within
    // each refined cell until orbit representatives carry that covariant action as well.
    options.automorphism_pruning = false;
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let (entity_keys, incidence_keys) = initial_class_keys(molecule, &incidence_graph)?;
    let classes = rank_initial_classes(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
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
    let no_prefix = |_: &OrderedPartition, _: &CanonicalCandidate<_>| false;
    let selected = canonical_search(
        &adapter,
        &descriptors,
        context.automorphism_algorithm,
        options,
        &leaf_candidate,
        &no_prefix,
    );
    let correspondence =
        correspondence_from_order(molecule, &incidence_graph, &selected.candidate.entity_order);

    let canonical = canonicalize_stereo_frames(molecule.remap(&correspondence))?;
    let canonical = normalize_molecule(canonical)?;
    Ok((canonical, correspondence))
}

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
    molecule.check_integrity()?;
    let normalized = normalize_molecule(molecule.clone())?;
    let molecule = &normalized;
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let (entity_keys, incidence_keys) = initial_class_keys(molecule, &incidence_graph)?;
    let classes = rank_initial_classes(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
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
    let no_prefix = |_: &OrderedPartition, _: &CanonicalCandidate<_>| false;
    let selected = canonical_search(
        &adapter,
        &descriptors,
        context.automorphism_algorithm,
        CanonicalSearchOptions {
            automorphism_pruning: false,
            ..options
        },
        &leaf_candidate,
        &no_prefix,
    );
    let correspondence =
        correspondence_from_order(molecule, &incidence_graph, &selected.candidate.entity_order);
    let (_, canonical) =
        complete_candidate(molecule, &incidence_graph, &selected.candidate.entity_order)?;
    Ok((canonical, correspondence))
}

fn canonical_key_by(
    molecule: &Molecule,
    level: CanonicalizeLevel,
    context: &CanonicalizeContext,
) -> Result<CanonicalComparisonKey, MoleculeCanonicalizeError> {
    molecule.check_integrity()?;
    if level == CanonicalizeLevel::Full {
        let normalized = normalize_molecule(molecule.clone())?;
        return canonical_key_by_full(&normalized, context);
    }
    let incidence_level = match level {
        CanonicalizeLevel::Topology => IncidenceLevel::Topology,
        CanonicalizeLevel::Constitution => IncidenceLevel::Constitution,
        CanonicalizeLevel::Structure => IncidenceLevel::Full,
        CanonicalizeLevel::Full => unreachable!("handled before incidence construction"),
    };
    let incidence_graph = molecule.incidence_graph(incidence_level);
    let (entity_keys, incidence_keys) = initial_class_keys(molecule, &incidence_graph)?;
    let classes = rank_initial_classes(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
    let no_prefix = |_: &OrderedPartition, _: &CanonicalCandidate<_>| false;
    let options = CanonicalSearchOptions {
        automorphism_pruning: level != CanonicalizeLevel::Structure,
        prefix_pruning: false,
        branch_order: backend_canonical_branch_order,
    };

    let key = match level {
        CanonicalizeLevel::Topology => {
            let descriptors = partition_descriptors(&adapter, &entity_keys, &incidence_keys);
            let leaf_candidate = |order: &[NodeId]| {
                topology_candidate(molecule, &incidence_graph, order)
                    .expect("initial classes established topology normalization")
            };
            canonical_search(
                &adapter,
                &descriptors,
                context.automorphism_algorithm,
                options,
                &leaf_candidate,
                &no_prefix,
            )
            .candidate
            .key
        }
        CanonicalizeLevel::Constitution => {
            let descriptors =
                constitution_partition_descriptors(&adapter, &entity_keys, &incidence_graph);
            let leaf_candidate = |order: &[NodeId]| {
                constitution_candidate(molecule, &incidence_graph, order)
                    .expect("initial classes established constitution normalization")
            };
            canonical_search(
                &adapter,
                &descriptors,
                context.automorphism_algorithm,
                options,
                &leaf_candidate,
                &no_prefix,
            )
            .candidate
            .key
        }
        CanonicalizeLevel::Structure => {
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
                &no_prefix,
            )
            .candidate
            .key
        }
        CanonicalizeLevel::Full => unreachable!("handled before selected-layer search"),
    };
    Ok(key)
}

fn canonical_key_by_full(
    molecule: &Molecule,
    context: &CanonicalizeContext,
) -> Result<CanonicalComparisonKey, MoleculeCanonicalizeError> {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let (entity_keys, incidence_keys) = initial_class_keys(molecule, &incidence_graph)?;
    let classes = rank_initial_classes(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
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
    .candidate
    .key)
}

fn reaction_span_entity_keys(
    span: &ReactionSpan,
    incidence_graph: &IncidenceGraph,
) -> Result<(Vec<InitialClassKey>, Vec<InitialClassKey>), Contradiction> {
    let entity_keys = incidence_graph
        .graph()
        .node_ids()
        .map(|node| reaction_span_entity_class_key(span, incidence_graph.entity(node)))
        .collect::<Result<Vec<_>, _>>()?;
    let incidence_keys = incidence_graph
        .incidences()
        .map(|(_, incidence)| incidence_key(incidence).map(InitialClassKey::Incidence))
        .collect::<Result<Vec<_>, _>>()?;
    Ok((entity_keys, incidence_keys))
}

fn reaction_span_comparison_key(
    span: &ReactionSpan,
    level: CanonicalizeLevel,
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
        CanonicalizeLevel::Constitution | CanonicalizeLevel::Structure | CanonicalizeLevel::Full
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
            .map(|(first, second, span)| {
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

    if matches!(
        level,
        CanonicalizeLevel::Structure | CanonicalizeLevel::Full
    ) {
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

    let constraints = if level == CanonicalizeLevel::Full {
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
                AtomId,
                AtomId,
                EntitySpan<super::noncovalent::NoncovalentBondForm>
            )| constraint_key_span(&entry.2, |form| form
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

impl Canonicalize for Molecule {
    type Error = MoleculeCanonicalizeError;

    fn canonicalize(self, context: &CanonicalizeContext) -> Result<Self, Self::Error> {
        canonicalize_full(&self, context)
    }

    fn canonicalize_with_correspondence(
        self,
        context: &CanonicalizeContext,
    ) -> Result<(Self, MoleculeCorrespondence), Self::Error> {
        canonicalize_full_with_options(
            &self,
            context,
            CanonicalSearchOptions {
                automorphism_pruning: false,
                prefix_pruning: false,
                branch_order: backend_canonical_branch_order,
            },
        )
    }

    fn canonicalize_by(
        self,
        level: CanonicalizeLevel,
        context: &CanonicalizeContext,
    ) -> Result<Self, Self::Error> {
        match level {
            CanonicalizeLevel::Topology => canonicalize_topology(&self, context),
            CanonicalizeLevel::Constitution => canonicalize_constitution(&self, context),
            CanonicalizeLevel::Structure => canonicalize_structure(&self, context),
            CanonicalizeLevel::Full => canonicalize_full(&self, context),
        }
    }

    fn canonical_hash_by(
        self,
        level: CanonicalizeLevel,
        context: &CanonicalizeContext,
    ) -> Result<u64, Self::Error> {
        if level == CanonicalizeLevel::Full {
            return self.canonical_hash(context);
        }
        canonical_key_by(&self, level, context).map(|key| hash_value(&key))
    }

    fn canonical_eq(&self, other: &Self, context: &CanonicalizeContext) -> bool {
        if self == other {
            return true;
        }
        match (
            canonical_key_by(self, CanonicalizeLevel::Full, context),
            canonical_key_by(other, CanonicalizeLevel::Full, context),
        ) {
            (Ok(left), Ok(right)) => left == right,
            (
                Err(MoleculeCanonicalizeError::Contradiction(_)),
                Err(MoleculeCanonicalizeError::Contradiction(_)),
            ) => true,
            _ => false,
        }
    }

    fn canonical_eq_by(
        &self,
        other: &Self,
        level: CanonicalizeLevel,
        context: &CanonicalizeContext,
    ) -> bool {
        if self == other {
            return true;
        }
        if level == CanonicalizeLevel::Full {
            return self.canonical_eq(other, context);
        }
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
    /// The molecule does not satisfy its representation-integrity contract.
    #[error(transparent)]
    Integrity(#[from] MoleculeIntegrityError),
    /// Intrinsic normalization of a carried value reached a contradiction.
    #[error(transparent)]
    Contradiction(#[from] Contradiction),
}

/// Failure to construct a canonical [`ReactionSpan`].
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ReactionSpanCanonicalizeError {
    /// The span does not satisfy its representation-integrity contract.
    #[error(transparent)]
    Integrity(#[from] ReactionSpanIntegrityError),
    /// Intrinsic normalization of a carried value reached a contradiction.
    #[error(transparent)]
    Contradiction(#[from] Contradiction),
}

fn reaction_span_canonical_candidate(
    span: &ReactionSpan,
    level: CanonicalizeLevel,
    context: &CanonicalizeContext,
) -> Result<CanonicalCandidate<CanonicalComparisonKey>, Contradiction> {
    let incidence_level = match level {
        CanonicalizeLevel::Topology => IncidenceLevel::Topology,
        CanonicalizeLevel::Constitution => IncidenceLevel::Constitution,
        CanonicalizeLevel::Structure | CanonicalizeLevel::Full => IncidenceLevel::Full,
    };
    let incidence_graph = span.incidence_graph(incidence_level);
    let (entity_keys, incidence_keys) = reaction_span_entity_keys(span, &incidence_graph)?;
    let classes = rank_initial_classes(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
    let descriptors = match level {
        CanonicalizeLevel::Topology => {
            partition_descriptors(&adapter, &entity_keys, &incidence_keys)
        }
        CanonicalizeLevel::Constitution
        | CanonicalizeLevel::Structure
        | CanonicalizeLevel::Full => {
            constitution_partition_descriptors(&adapter, &entity_keys, &incidence_graph)
        }
    };
    let leaf_candidate = |order: &[NodeId]| {
        let correspondence = lhs_anchored_correspondence_from_order(span, &incidence_graph, order);
        let remapped = span.remap(&correspondence);
        let remapped = if matches!(
            level,
            CanonicalizeLevel::Structure | CanonicalizeLevel::Full
        ) {
            canonicalize_reaction_span_stereo_frames(remapped)
                .expect("initial classes established stereo normalization")
        } else {
            remapped
        };
        CanonicalCandidate {
            key: reaction_span_comparison_key(&remapped, level)
                .expect("initial classes established span normalization"),
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
    level: CanonicalizeLevel,
    candidate: &CanonicalCandidate<CanonicalComparisonKey>,
) -> Result<ReactionSpan, Contradiction> {
    let incidence_level = match level {
        CanonicalizeLevel::Topology => IncidenceLevel::Topology,
        CanonicalizeLevel::Constitution => IncidenceLevel::Constitution,
        CanonicalizeLevel::Structure | CanonicalizeLevel::Full => IncidenceLevel::Full,
    };
    let incidence_graph = span.incidence_graph(incidence_level);
    let correspondence =
        lhs_anchored_correspondence_from_order(span, &incidence_graph, &candidate.entity_order);
    let remapped = span.remap(&correspondence);
    if matches!(
        level,
        CanonicalizeLevel::Structure | CanonicalizeLevel::Full
    ) {
        Ok(canonicalize_reaction_span_stereo_frames(remapped)?)
    } else {
        Ok(remapped)
    }
}

fn canonicalize_reaction_span_by(
    span: &ReactionSpan,
    level: CanonicalizeLevel,
    context: &CanonicalizeContext,
) -> Result<ReactionSpan, ReactionSpanCanonicalizeError> {
    span.check_integrity()?;
    Ok(canonicalize_checked_reaction_span_by(span, level, context)?)
}

fn canonicalize_checked_reaction_span_by(
    span: &ReactionSpan,
    level: CanonicalizeLevel,
    context: &CanonicalizeContext,
) -> Result<ReactionSpan, Contradiction> {
    Ok(canonicalize_checked_reaction_span_with_correspondence_by(span, level, context)?.0)
}

fn canonicalize_checked_reaction_span_with_correspondence_by(
    span: &ReactionSpan,
    level: CanonicalizeLevel,
    context: &CanonicalizeContext,
) -> Result<(ReactionSpan, MoleculeCorrespondence), Contradiction> {
    let normalized = normalize_reaction_span(span.clone())?;
    let candidate = reaction_span_canonical_candidate(&normalized, level, context)?;
    let incidence_level = match level {
        CanonicalizeLevel::Topology => IncidenceLevel::Topology,
        CanonicalizeLevel::Constitution => IncidenceLevel::Constitution,
        CanonicalizeLevel::Structure | CanonicalizeLevel::Full => IncidenceLevel::Full,
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
    level: CanonicalizeLevel,
    context: &CanonicalizeContext,
) -> Result<CanonicalComparisonKey, ReactionSpanCanonicalizeError> {
    span.check_integrity()?;
    if level == CanonicalizeLevel::Full {
        let normalized = normalize_reaction_span(span.clone())?;
        return Ok(reaction_span_canonical_candidate(&normalized, level, context)?.key);
    }
    Ok(reaction_span_canonical_candidate(span, level, context)?.key)
}

impl Canonicalize for ReactionSpan {
    type Error = ReactionSpanCanonicalizeError;

    fn canonicalize(self, context: &CanonicalizeContext) -> Result<Self, Self::Error> {
        canonicalize_reaction_span_by(&self, CanonicalizeLevel::Full, context)
    }

    fn canonicalize_with_correspondence(
        self,
        context: &CanonicalizeContext,
    ) -> Result<(Self, MoleculeCorrespondence), Self::Error> {
        self.check_integrity()?;
        Ok(canonicalize_checked_reaction_span_with_correspondence_by(
            &self,
            CanonicalizeLevel::Full,
            context,
        )?)
    }

    fn canonicalize_by(
        self,
        level: CanonicalizeLevel,
        context: &CanonicalizeContext,
    ) -> Result<Self, Self::Error> {
        canonicalize_reaction_span_by(&self, level, context)
    }

    fn canonical_hash_by(
        self,
        level: CanonicalizeLevel,
        context: &CanonicalizeContext,
    ) -> Result<u64, Self::Error> {
        if level == CanonicalizeLevel::Full {
            return self.canonical_hash(context);
        }
        canonical_reaction_span_key(&self, level, context).map(|key| hash_value(&key))
    }

    fn canonical_eq(&self, other: &Self, context: &CanonicalizeContext) -> bool {
        if self == other {
            return true;
        }
        match (
            canonical_reaction_span_key(self, CanonicalizeLevel::Full, context),
            canonical_reaction_span_key(other, CanonicalizeLevel::Full, context),
        ) {
            (Ok(left), Ok(right)) => left == right,
            (
                Err(ReactionSpanCanonicalizeError::Contradiction(_)),
                Err(ReactionSpanCanonicalizeError::Contradiction(_)),
            ) => true,
            _ => false,
        }
    }

    fn canonical_eq_by(
        &self,
        other: &Self,
        level: CanonicalizeLevel,
        context: &CanonicalizeContext,
    ) -> bool {
        if self == other {
            return true;
        }
        if level == CanonicalizeLevel::Full {
            return self.canonical_eq(other, context);
        }
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
    /// The reaction does not satisfy its representation-integrity contract.
    #[error(transparent)]
    Integrity(#[from] ReactionIntegrityError),
    /// Intrinsic normalization or span materialization reached a contradiction.
    #[error(transparent)]
    Contradiction(#[from] Contradiction),
}

fn reaction_delta_is_selected(delta: &Delta, level: CanonicalizeLevel) -> bool {
    let includes_non_stereo = level != CanonicalizeLevel::Topology;
    let includes_stereo = matches!(
        level,
        CanonicalizeLevel::Structure | CanonicalizeLevel::Full
    );
    let includes_constraints = level == CanonicalizeLevel::Full;
    match delta {
        Delta::Atom(AtomDelta::ModifyConstraint { .. })
        | Delta::Bond(BondDelta::ModifyConstraint { .. })
        | Delta::Constraint(_) => includes_constraints,
        Delta::Atom(_) | Delta::Bond(_) => true,
        Delta::DativeBond(DativeBondDelta::ModifyConstraint { .. })
        | Delta::AromaticSystem(AromaticSystemDelta::ModifyConstraint { .. })
        | Delta::MulticenterBond(MulticenterBondDelta::ModifyConstraint { .. })
        | Delta::NoncovalentBond(NoncovalentBondDelta::ModifyConstraint { .. }) => {
            includes_non_stereo && includes_constraints
        }
        Delta::DativeBond(_)
        | Delta::AromaticSystem(_)
        | Delta::MulticenterBond(_)
        | Delta::NoncovalentBond(_) => includes_non_stereo,
        Delta::StereoAtom(StereoAtomDelta::ModifyConstraint { .. })
        | Delta::StereoBond(StereoBondDelta::ModifyConstraint { .. }) => {
            includes_stereo && includes_constraints
        }
        Delta::StereoAtom(_) | Delta::StereoBond(_) => includes_stereo,
    }
}

fn project_reaction(reaction: &Reaction, level: CanonicalizeLevel) -> Reaction {
    let deltas = reaction
        .deltas
        .iter()
        .filter(|delta| reaction_delta_is_selected(delta, level))
        .cloned()
        .collect::<Deltas>();
    Reaction::new(reaction.lhs.clone(), deltas)
}

fn canonicalize_reaction_by(
    reaction: &Reaction,
    level: CanonicalizeLevel,
    context: &CanonicalizeContext,
) -> Result<Reaction, ReactionCanonicalizeError> {
    reaction.check_integrity()?;
    let span = reaction.to_reaction_span()?;
    Ok(canonicalize_checked_reaction_span_by(&span, level, context)?.to_reaction())
}

impl Canonicalize for Reaction {
    type Error = ReactionCanonicalizeError;

    fn canonicalize(self, context: &CanonicalizeContext) -> Result<Self, Self::Error> {
        canonicalize_reaction_by(&self, CanonicalizeLevel::Full, context)
    }

    fn canonicalize_with_correspondence(
        self,
        context: &CanonicalizeContext,
    ) -> Result<(Self, MoleculeCorrespondence), Self::Error> {
        self.check_integrity()?;
        let span = self.to_reaction_span()?;
        let (canonical, correspondence) =
            canonicalize_checked_reaction_span_with_correspondence_by(
                &span,
                CanonicalizeLevel::Full,
                context,
            )?;
        Ok((canonical.to_reaction(), correspondence))
    }

    fn canonicalize_by(
        self,
        level: CanonicalizeLevel,
        context: &CanonicalizeContext,
    ) -> Result<Self, Self::Error> {
        canonicalize_reaction_by(&self, level, context)
    }

    fn canonical_hash_by(
        self,
        level: CanonicalizeLevel,
        context: &CanonicalizeContext,
    ) -> Result<u64, Self::Error> {
        if level == CanonicalizeLevel::Full {
            return self.canonical_hash(context);
        }
        self.check_integrity()?;
        let span = project_reaction(&self, level).to_reaction_span()?;
        let key = reaction_span_canonical_candidate(&span, level, context)?.key;
        Ok(hash_value(&key))
    }

    fn canonical_eq(&self, other: &Self, context: &CanonicalizeContext) -> bool {
        if self == other {
            return true;
        }
        match (
            canonicalize_reaction_by(self, CanonicalizeLevel::Full, context),
            canonicalize_reaction_by(other, CanonicalizeLevel::Full, context),
        ) {
            (Ok(left), Ok(right)) => left == right,
            (
                Err(ReactionCanonicalizeError::Contradiction(_)),
                Err(ReactionCanonicalizeError::Contradiction(_)),
            ) => true,
            _ => false,
        }
    }

    fn canonical_eq_by(
        &self,
        other: &Self,
        level: CanonicalizeLevel,
        context: &CanonicalizeContext,
    ) -> bool {
        if self == other {
            return true;
        }
        if level == CanonicalizeLevel::Full {
            return self.canonical_eq(other, context);
        }
        if self.check_integrity().is_err() || other.check_integrity().is_err() {
            return false;
        }
        match (
            project_reaction(self, level).to_reaction_span(),
            project_reaction(other, level).to_reaction_span(),
        ) {
            (Ok(left), Ok(right)) => left.canonical_eq_by(&right, level, context),
            (Err(_), Err(_)) => true,
            _ => false,
        }
    }
}

#[cfg(test)]
mod tests;
