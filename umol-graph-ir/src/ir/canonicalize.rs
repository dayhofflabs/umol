//! Aggregate canonicalization inputs and failures.

use std::cell::OnceCell;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fmt::Debug;
use std::hash::{DefaultHasher, Hash, Hasher};
use std::{iter, mem};

use thiserror::Error;
use umol_graph_core::{
    AutomorphismAlgorithm, AutomorphismOutput, EdgeId, Graph, GraphRemapping, NodeId, Remapping,
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
use super::delta::Delta;
use super::electrons::ElectronCountsForm;
use super::error::Contradiction;
use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::incidence::{Incidence, IncidenceGraph};
use super::ligand::StereoLigandKind;
use super::molecule::Molecule;
use super::noncovalent::{NoncovalentBondKind, NoncovalentBondKindForm};
use super::num::{ArithExpr, NumForm, PredExpr};
use super::operators::{MemOp, RelOp};
use super::reaction::Reaction;
use super::reaction_span::ReactionSpan;
use super::remap::MoleculeRemapping;
use super::span::{ConstraintSpan, EntitySpan};
use super::spin::UnpairedElectronsForm;
use super::stereo::{
    CisTransStereoForm, StereoConfigurationForm, StereoCoset, StereoKind, StereoTerm,
    Stereogenicity, TetrahedralStereoForm, Topicity,
};
use super::traits::{Normalize, Reframe};

mod molecule;
mod reaction_span;

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

fn delta_canonicalize_level(delta: &Delta) -> DescriptionLevel {
    match delta {
        Delta::Atom(_) | Delta::Bond(_) => DescriptionLevel::Topology,
        Delta::DativeBond(_)
        | Delta::AromaticSystem(_)
        | Delta::MulticenterBond(_)
        | Delta::NoncovalentBond(_) => DescriptionLevel::Constitution,
        Delta::StereoAtom(_) | Delta::StereoBond(_) => DescriptionLevel::Structure,
        Delta::Constraint(_) => DescriptionLevel::Full,
    }
}

fn molecule_canonicalize_level(molecule: &Molecule) -> DescriptionLevel {
    if molecule.has_constraints() {
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
        DescriptionLevel::Full
    } else if span.stereo_atoms().count() > 0 || span.stereo_bonds().count() > 0 {
        DescriptionLevel::Structure
    } else if span.dative_bonds().count() > 0
        || span.aromatic_systems().count() > 0
        || span.multicenter_bonds().count() > 0
        || span.noncovalent_bonds().count() > 0
    {
        DescriptionLevel::Constitution
    } else {
        DescriptionLevel::Topology
    }
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
    /// unsatisfiable or a carried attribute cannot be transported in its stored frame.
    fn canonicalize(self, context: &CanonicalizeContext) -> Result<Self, Self::Error>;

    /// Construct the complete canonical form and its source-to-canonical remapping.
    ///
    /// The remapping maps every entity id in the input frame to its id in the returned
    /// canonical frame. Each of its eight entity-kind components is total on both sides and
    /// therefore represents a dense bijection. For [`Reaction`], the two frames are the complete
    /// union frames of the materialized input and returned reaction spans.
    ///
    /// # Errors
    ///
    /// Returns the same aggregate-specific canonicalization error as
    /// [`Self::canonicalize`]. A [`Reaction`] also reports failure to materialize its reaction span.
    ///
    /// # Semantic properties
    ///
    /// Discarding the remapping yields exactly [`Self::canonicalize`] under the same context.
    /// Remapping the input through the remapping and then applying [`Reframe::reframe`]
    /// reconstructs the returned canonical value.
    fn tracked_canonicalize(
        self,
        context: &CanonicalizeContext,
    ) -> Result<(Self, MoleculeRemapping), Self::Error>;

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
    /// Successful hashes are invariant under valid dense entity remapping. Canonically equal
    /// values have equal canonical hashes when both hash operations succeed.
    fn canonical_hash(&self, context: &CanonicalizeContext) -> Result<u64, Self::Error>;

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
    idx: u16,
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

    const fn new(domain: StructuralDomainPosition, idx: u16) -> Self {
        Self { domain, idx }
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
    const REMOVED: Self = Self(1);
    const MODIFIED: Self = Self(2);
    const ADDED: Self = Self(3);
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
    idx: u16,
}

impl RelationalConstraintPosition {
    const fn new(entity: EntityBlockPosition, idx: u16) -> Self {
        Self { entity, idx }
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
        TetrahedralStereoForm::Stereo => variant(2, []),
    }
}

fn cis_trans_stereo_form_key(value: &CisTransStereoForm) -> CanonicalKeyValue {
    match value {
        CisTransStereoForm::Undetermined => variant(0, []),
        CisTransStereoForm::NotStereo => variant(1, []),
        CisTransStereoForm::Stereo => variant(2, []),
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
            StereoLigandKind::Virtual => 3,
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
    let position = |entity: EntityBlockPosition, idx: u64| {
        let position = RelationalConstraintPosition::new(entity, idx as u16);
        product([
            CanonicalKeyValue::Unsigned(position.entity.domain.0.into()),
            CanonicalKeyValue::Unsigned(position.entity.idx.into()),
            CanonicalKeyValue::Unsigned(position.idx.into()),
        ])
    };
    let atom_ids = |ids: &[AtomId]| index_sequence(ids.iter().map(|id| id.index()));

    let (entity, idx, fields) = match value {
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

    product([position(entity, idx), sequence(fields)])
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
                    CanonicalKeyValue::Unsigned(EntityBlockPosition::ATOM.idx.into()),
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
                    CanonicalKeyValue::Unsigned(EntityBlockPosition::BOND.idx.into()),
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
                    CanonicalKeyValue::Unsigned(EntityBlockPosition::DATIVE_BOND.idx.into()),
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
                    CanonicalKeyValue::Unsigned(EntityBlockPosition::AROMATIC_SYSTEM.idx.into()),
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
                    CanonicalKeyValue::Unsigned(EntityBlockPosition::MULTICENTER_BOND.idx.into()),
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
                    CanonicalKeyValue::Unsigned(EntityBlockPosition::NONCOVALENT_BOND.idx.into()),
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
                    CanonicalKeyValue::Unsigned(EntityBlockPosition::STEREO_ATOM.idx.into()),
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
                    CanonicalKeyValue::Unsigned(EntityBlockPosition::STEREO_BOND.idx.into()),
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

fn stereo_ligand_key(atom: u32, kind: StereoLigandKind) -> CanonicalKeyValue {
    product([
        CanonicalKeyValue::Unsigned(atom.into()),
        stereo_ligand_kind_key(kind),
    ])
}

fn incidence_key(incidence: &Incidence) -> CanonicalKeyValue {
    match incidence {
        Incidence::BondEndpoint => variant(0, []),
        Incidence::DativeDonor => variant(1, []),
        Incidence::DativeAcceptor => variant(2, []),
        Incidence::AromaticAtom => variant(3, []),
        Incidence::MulticenterAtom => variant(4, []),
        Incidence::NoncovalentEndpoint => variant(5, []),
        Incidence::StereoSite => variant(6, []),
        Incidence::StereoLigand(kind) => variant(7, [stereo_ligand_kind_key(*kind)]),
    }
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

/// Colored simple graph and source-domain projection used by canonical search.
///
/// Source nodes occupy the leading `source_node_count` positions. Remaining vertices encode
/// subdivided incidences or bonds and are excluded when backend automorphisms are projected onto
/// the search domain.
#[derive(Clone, Debug)]
struct AutomorphismAdapter {
    graph: Graph,
    colors: Vec<AutomorphismColor>,
    node_sources: Vec<SubdivisionNodeSource>,
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
            SubdivisionNodeSource::Node(source) if source.index() < self.source_node_count => {
                source
            }
            SubdivisionNodeSource::Node(_) => {
                unreachable!("disjoint colors preserve the source-node domain")
            }
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
                SubdivisionNodeSource::Node(source) if source.index() < self.source_node_count => {
                    Some(source)
                }
                SubdivisionNodeSource::Node(_) => None,
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

/// Minimum carrier size at which refinement tracks only cells created by the preceding split.
const ACTIVE_CELL_REFINEMENT_MIN_NODE_COUNT: usize = 32;

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

    fn from_color_ranks(colors: impl IntoIterator<Item = u32>) -> Self {
        let mut cells = Vec::<Vec<NodeId>>::new();
        for (index, color) in colors.into_iter().enumerate() {
            let color = color as usize;
            if cells.len() <= color {
                cells.resize_with(color + 1, Vec::new);
            }
            cells[color].push(NodeId(index as u32));
        }
        cells.retain(|cell| !cell.is_empty());
        Self { cells }
    }

    fn refine(self, graph: &Graph) -> Self {
        if graph.node_count() < ACTIVE_CELL_REFINEMENT_MIN_NODE_COUNT {
            self.refine_by_all_cells(graph)
        } else {
            self.refine_by_active_cells(graph)
        }
    }

    fn refine_by_active_cells(mut self, graph: &Graph) -> Self {
        let node_count = graph.node_count();
        let mut active_cells = vec![true; self.cells.len()];
        let mut node_active_columns = vec![usize::MAX; node_count];
        let mut signatures = Vec::<u32>::new();
        let mut next_cells = Vec::<Vec<NodeId>>::with_capacity(self.cells.len());
        let mut next_active_cells = Vec::<bool>::with_capacity(self.cells.len());

        loop {
            node_active_columns.fill(usize::MAX);
            let mut active_count = 0;
            for (cell_index, cell) in self.cells.iter().enumerate() {
                if active_cells[cell_index] {
                    for node in cell {
                        node_active_columns[node.index()] = active_count;
                    }
                    active_count += 1;
                }
            }
            signatures.clear();
            signatures.resize(node_count * active_count, 0);
            for node in graph.node_ids() {
                let offset = node.index() * active_count;
                for neighbor in graph.neighbors(node) {
                    let column = node_active_columns[neighbor.node.index()];
                    if column != usize::MAX {
                        signatures[offset + column] += 1;
                    }
                }
            }

            let mut changed = false;
            for mut cell in self.cells.drain(..) {
                let signature = |node: NodeId| {
                    let offset = node.index() * active_count;
                    &signatures[offset..offset + active_count]
                };
                cell.sort_unstable_by(|lhs, rhs| {
                    signature(*rhs)
                        .cmp(signature(*lhs))
                        .then_with(|| lhs.cmp(rhs))
                });

                let mut group_start = 0;
                for index in 1..cell.len() {
                    if signature(cell[index]) != signature(cell[group_start]) {
                        next_cells.push(cell[group_start..index].to_vec());
                        group_start = index;
                    }
                }
                if group_start == 0 {
                    next_cells.push(cell);
                    next_active_cells.push(false);
                } else {
                    next_cells.push(cell[group_start..].to_vec());
                    let new_cell_count = next_cells.len() - next_active_cells.len();
                    next_active_cells.extend(iter::repeat_n(true, new_cell_count));
                    changed = true;
                }
            }

            mem::swap(&mut self.cells, &mut next_cells);
            if !changed {
                return self;
            }
            mem::swap(&mut active_cells, &mut next_active_cells);
            next_active_cells.clear();
        }
    }

    fn refine_by_all_cells(mut self, graph: &Graph) -> Self {
        let node_count = graph.node_count();
        let mut cell_indices = vec![0; node_count];
        let mut signatures = Vec::<u32>::new();
        let mut next_cells = Vec::<Vec<NodeId>>::with_capacity(self.cells.len());

        loop {
            for (cell_index, cell) in self.cells.iter().enumerate() {
                for node in cell {
                    cell_indices[node.index()] = cell_index as u32;
                }
            }
            let cell_count = self.cells.len();
            signatures.clear();
            signatures.resize(node_count * cell_count, 0);
            for node in graph.node_ids() {
                let offset = node.index() * cell_count;
                for neighbor in graph.neighbors(node) {
                    signatures[offset + cell_indices[neighbor.node.index()] as usize] += 1;
                }
            }

            let mut changed = false;
            for mut cell in self.cells.drain(..) {
                let signature = |node: NodeId| {
                    let offset = node.index() * cell_count;
                    &signatures[offset..offset + cell_count]
                };
                cell.sort_unstable_by(|lhs, rhs| {
                    signature(*rhs)
                        .cmp(signature(*lhs))
                        .then_with(|| lhs.cmp(rhs))
                });

                let mut group_start = 0;
                for index in 1..cell.len() {
                    if signature(cell[index]) != signature(cell[group_start]) {
                        next_cells.push(cell[group_start..index].to_vec());
                        group_start = index;
                    }
                }
                if group_start == 0 {
                    next_cells.push(cell);
                } else {
                    next_cells.push(cell[group_start..].to_vec());
                    changed = true;
                }
            }

            mem::swap(&mut self.cells, &mut next_cells);
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

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum BranchOrdering {
    BackendCanonical,
    NodeId,
    #[cfg(test)]
    ReverseNodeId,
}

impl BranchOrdering {
    fn requires_automorphisms(self) -> bool {
        self == Self::BackendCanonical
    }

    fn order(
        self,
        adapter: &AutomorphismAdapter,
        automorphisms: Option<&AutomorphismAdapterOutput>,
        candidates: &mut [NodeId],
    ) {
        match self {
            Self::BackendCanonical => {
                let labels = &automorphisms
                    .expect("backend branch order requested automorphisms")
                    .source_canonical_labels;
                let mut ranks = vec![0; adapter.source_node_count];
                for (rank, node) in labels.iter().enumerate() {
                    ranks[node.index()] = rank;
                }
                candidates.sort_unstable_by_key(|node| ranks[node.index()]);
            }
            Self::NodeId => candidates.sort_unstable(),
            #[cfg(test)]
            Self::ReverseNodeId => candidates.sort_unstable_by(|lhs, rhs| rhs.cmp(lhs)),
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct CanonicalSearchOptions {
    automorphism_pruning: bool,
    branch_order: BranchOrdering,
}

#[cfg(test)]
#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct CanonicalSearchStats {
    initial_residual_cell_sizes: Vec<usize>,
    refinement_calls: usize,
    branch_order_calls: usize,
    backend_calls: usize,
    visited_leaves: usize,
    leaf_comparisons: usize,
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
    #[cfg(test)]
    stats: CanonicalSearchStats,
}

/// Minimize a typed leaf candidate over an exact partition ordered by semantic descriptors.
///
/// Adapter colors are opaque equality labels and do not order the partition. Automorphism pruning
/// requires the leaf key to be invariant under adapter automorphisms.
fn canonical_search<K, Descriptor, LeafCandidate>(
    adapter: &AutomorphismAdapter,
    partition_descriptors: &[Descriptor],
    algorithm: AutomorphismAlgorithm,
    options: CanonicalSearchOptions,
    leaf_candidate: &LeafCandidate,
) -> CanonicalSearchResult<K>
where
    K: Ord,
    Descriptor: Clone + Ord,
    LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
{
    canonical_search_impl(
        adapter,
        OrderedPartition::from_descriptors(partition_descriptors),
        algorithm,
        options,
        leaf_candidate,
        false,
        &|_| {},
    )
}

fn canonical_search_from_partition<K, LeafCandidate>(
    adapter: &AutomorphismAdapter,
    partition: OrderedPartition,
    algorithm: AutomorphismAlgorithm,
    options: CanonicalSearchOptions,
    leaf_candidate: &LeafCandidate,
) -> CanonicalSearchResult<K>
where
    K: Ord,
    LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
{
    canonical_search_impl(
        adapter,
        partition,
        algorithm,
        options,
        leaf_candidate,
        false,
        &|_| {},
    )
}

fn canonical_search_with_generator_filter<K, Descriptor, LeafCandidate, FilterGenerators>(
    adapter: &AutomorphismAdapter,
    partition_descriptors: &[Descriptor],
    algorithm: AutomorphismAlgorithm,
    options: CanonicalSearchOptions,
    leaf_candidate: &LeafCandidate,
    filter_generators: &FilterGenerators,
) -> CanonicalSearchResult<K>
where
    K: Ord,
    Descriptor: Clone + Ord,
    LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
    FilterGenerators: Fn(&mut Vec<Vec<NodeId>>),
{
    canonical_search_impl(
        adapter,
        OrderedPartition::from_descriptors(partition_descriptors),
        algorithm,
        options,
        leaf_candidate,
        true,
        filter_generators,
    )
}

fn canonical_search_impl<K, LeafCandidate, FilterGenerators>(
    adapter: &AutomorphismAdapter,
    initial: OrderedPartition,
    algorithm: AutomorphismAlgorithm,
    options: CanonicalSearchOptions,
    leaf_candidate: &LeafCandidate,
    apply_generator_filter: bool,
    filter_generators: &FilterGenerators,
) -> CanonicalSearchResult<K>
where
    K: Ord,
    LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
    FilterGenerators: Fn(&mut Vec<Vec<NodeId>>),
{
    let initial = initial.refine(adapter.graph());
    let mut best = None;
    #[cfg(test)]
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
        apply_generator_filter,
        filter_generators,
        &mut best,
        #[cfg(test)]
        &mut stats,
    );

    CanonicalSearchResult {
        candidate: best.expect("every finite partition has a discrete entity refinement"),
        #[cfg(test)]
        stats,
    }
}

#[allow(clippy::too_many_arguments)]
fn search_partition<K, LeafCandidate, FilterGenerators>(
    adapter: &AutomorphismAdapter,
    partition: OrderedPartition,
    algorithm: AutomorphismAlgorithm,
    options: CanonicalSearchOptions,
    leaf_candidate: &LeafCandidate,
    apply_generator_filter: bool,
    filter_generators: &FilterGenerators,
    best: &mut Option<CanonicalCandidate<K>>,
    #[cfg(test)] stats: &mut CanonicalSearchStats,
) where
    K: Ord,
    LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
    FilterGenerators: Fn(&mut Vec<Vec<NodeId>>),
{
    let Some(cell_index) = partition.first_non_singleton_entity_cell(adapter.source_node_count)
    else {
        #[cfg(test)]
        {
            stats.visited_leaves += 1;
        }
        let entity_order = partition.entity_order(adapter.source_node_count);
        let candidate = leaf_candidate(&entity_order);
        let improves = if let Some(best) = best.as_ref() {
            #[cfg(test)]
            {
                stats.leaf_comparisons += 1;
            }
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
    let automorphisms = (options.automorphism_pruning
        || apply_generator_filter
        || options.branch_order.requires_automorphisms())
    .then(|| {
        #[cfg(test)]
        {
            stats.backend_calls += 1;
        }
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
        #[cfg(test)]
        {
            stats.orbit_pruned_branches +=
                partition.cells[cell_index].len() - representatives.len();
        }
        candidates = representatives.into_values().collect();
    }

    #[cfg(test)]
    {
        stats.branch_order_calls += 1;
    }
    options
        .branch_order
        .order(adapter, automorphisms.as_ref(), &mut candidates);

    for candidate in candidates {
        #[cfg(test)]
        {
            stats.refinement_calls += 1;
        }
        let child = partition
            .individualize(cell_index, candidate)
            .refine(adapter.graph());
        search_partition(
            adapter,
            child,
            algorithm,
            options,
            leaf_candidate,
            apply_generator_filter,
            filter_generators,
            best,
            #[cfg(test)]
            stats,
        );
    }
}

/// A leaf key that is built when it is first compared or read.
struct DeferredKey<'a, K, Build> {
    order: Vec<NodeId>,
    build: &'a Build,
    key: OnceCell<K>,
}

impl<'a, K, Build: Fn(&[NodeId]) -> K> DeferredKey<'a, K, Build> {
    fn new(order: &[NodeId], build: &'a Build) -> Self {
        Self {
            order: order.to_vec(),
            build,
            key: OnceCell::new(),
        }
    }

    fn get(&self) -> &K {
        self.key.get_or_init(|| (self.build)(&self.order))
    }

    /// The key, built now if no comparison built it.
    fn into_key(self) -> K {
        self.get();
        self.key.into_inner().expect("the key was built")
    }

    /// The key, if a comparison built it.
    fn into_built(self) -> Option<K> {
        self.key.into_inner()
    }
}

impl<K: Ord, Build: Fn(&[NodeId]) -> K> PartialEq for DeferredKey<'_, K, Build> {
    fn eq(&self, other: &Self) -> bool {
        self.get() == other.get()
    }
}

impl<K: Ord, Build: Fn(&[NodeId]) -> K> Eq for DeferredKey<'_, K, Build> {}

impl<K: Ord, Build: Fn(&[NodeId]) -> K> PartialOrd for DeferredKey<'_, K, Build> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<K: Ord, Build: Fn(&[NodeId]) -> K> Ord for DeferredKey<'_, K, Build> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.get().cmp(other.get())
    }
}

/// Topology search input with atoms and non-modal bond subdivisions as vertices.
///
/// Bonds in the largest normalized bond-form class are represented directly as edges. Every other
/// bond remains a colored subdivision vertex. The complete typed leaf key and remapping still
/// order and map every bond.
#[derive(Clone, Debug)]
struct TopologySearchInput {
    adapter: AutomorphismAdapter,
    initial_partition: OrderedPartition,
}

/// Construct constitution-level partition descriptors for the exact adapter.
///
/// An entity node takes its key and an occurrence node its incidence color, so dative roles and
/// positional electron counts split the search partition.
fn constitution_partition_descriptors(
    adapter: &AutomorphismAdapter,
    entity_keys: &[InitialColorKey],
) -> Vec<InitialColorKey> {
    adapter
        .node_sources
        .iter()
        .zip(&adapter.colors)
        .map(|(source, color)| match (*source, *color) {
            (SubdivisionNodeSource::Node(node), _) => entity_keys[node.index()].clone(),
            (SubdivisionNodeSource::Edge(_), AutomorphismColor::Incidence(rank)) => {
                InitialColorKey::Incidence(CanonicalKeyValue::Unsigned(rank.into()))
            }
            (SubdivisionNodeSource::Edge(_), AutomorphismColor::Entity(_)) => {
                unreachable!("an occurrence node has an incidence color")
            }
        })
        .collect()
}

/// Assign dense ordered colors to normalized entity and incidence keys.
fn rank_initial_colors(
    entity_keys: &[InitialColorKey],
    incidence_keys: &[InitialColorKey],
) -> InitialColors {
    let entity_count = entity_keys.len();
    let key = |index: usize| {
        if index < entity_count {
            &entity_keys[index]
        } else {
            &incidence_keys[index - entity_count]
        }
    };
    let mut indices = (0..entity_count + incidence_keys.len()).collect::<Vec<_>>();
    indices.sort_unstable_by(|&left, &right| key(left).cmp(key(right)));

    let mut ranks = vec![0; indices.len()];
    let mut previous = None;
    let mut rank = 0;
    for index in indices {
        if previous.is_some_and(|previous| key(previous) != key(index)) {
            rank += 1;
        }
        ranks[index] = rank;
        previous = Some(index);
    }
    let incidences = ranks.split_off(entity_count);
    InitialColors {
        entities: ranks,
        incidences,
    }
}

/// Assign dense integers to the colors of one entity kind or of the incidences, in the order of the
/// colors and starting at `first`. Each color comes with the position of its integer in `ranks`.
/// The integer after the last one used is returned.
fn rank_colors<T: Ord>(mut colors: Vec<(T, usize)>, first: u32, ranks: &mut [u32]) -> u32 {
    colors.sort_unstable_by(|(left, _), (right, _)| left.cmp(right));
    let mut rank = first;
    let mut previous = None;
    for (color, position) in &colors {
        if previous.is_some_and(|previous| previous != color) {
            rank += 1;
        }
        ranks[*position] = rank;
        previous = Some(color);
    }
    if colors.is_empty() {
        first
    } else {
        rank + 1
    }
}

fn atom_inherent_fields(attributes: &AtomForm) -> Result<Vec<FieldKey>, Contradiction> {
    Ok(vec![
        field(
            0,
            element_form_key(&attributes.element.clone().normalize()?),
        ),
        field(
            1,
            isotope_mass_form_key(&attributes.isotope_mass.clone().normalize()?),
        ),
        field(2, num_form_key(&attributes.charge.clone().normalize()?)),
        field(
            3,
            num_form_key(&attributes.implicit_hydrogens.clone().normalize()?),
        ),
        field(4, num_form_key(&attributes.lone_pairs.clone().normalize()?)),
        field(
            5,
            unpaired_electrons_form_key(&attributes.unpaired_electrons.clone().normalize()?),
        ),
    ])
}

fn bond_inherent_fields(attributes: &BondForm) -> Result<Vec<FieldKey>, Contradiction> {
    Ok(vec![
        field(1, num_form_key(&attributes.order.clone().normalize()?)),
        field(2, num_form_key(&attributes.charge.clone().normalize()?)),
        field(
            3,
            unpaired_electrons_form_key(&attributes.unpaired_electrons.clone().normalize()?),
        ),
    ])
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

fn molecule_remapping(images: &[Vec<usize>; 8]) -> MoleculeRemapping {
    fn remapping<Id: Copy + Debug + Into<usize> + From<usize>>(images: &[usize]) -> Remapping<Id> {
        Remapping::new(images.iter().copied().map(Id::from).collect())
            .expect("selected entity order is a permutation")
    }
    MoleculeRemapping::new(
        GraphRemapping::new(
            remapping::<NodeId>(&images[0]),
            remapping::<EdgeId>(&images[1]),
        ),
        remapping::<DativeBondId>(&images[2]),
        remapping::<AromaticSystemId>(&images[3]),
        remapping::<MulticenterBondId>(&images[4]),
        remapping::<NoncovalentBondId>(&images[5]),
        remapping::<StereoAtomId>(&images[6]),
        remapping::<StereoBondId>(&images[7]),
    )
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

/// Failure to construct a canonical [`Molecule`].
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum MoleculeCanonicalizeError {
    /// A carried value is contradictory or incompatible with its stored atom or ligand frame.
    #[error(transparent)]
    Contradiction(#[from] Contradiction),
}

/// Failure to construct a canonical [`ReactionSpan`].
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ReactionSpanCanonicalizeError {
    /// A carried value is contradictory or incompatible with its stored atom or ligand frame.
    #[error(transparent)]
    Contradiction(#[from] Contradiction),
}

/// Failure to construct a canonical [`Reaction`].
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ReactionCanonicalizeError {
    /// Normalization or span materialization failed, or literal electron counts
    /// do not match their atom list.
    #[error(transparent)]
    Contradiction(#[from] Contradiction),
}

impl From<ReactionSpanCanonicalizeError> for ReactionCanonicalizeError {
    fn from(error: ReactionSpanCanonicalizeError) -> Self {
        match error {
            ReactionSpanCanonicalizeError::Contradiction(contradiction) => {
                Self::Contradiction(contradiction)
            }
        }
    }
}

#[cfg(test)]
mod tests;
