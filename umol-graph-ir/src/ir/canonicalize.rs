//! Aggregate canonicalization inputs and failures.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

use thiserror::Error;
use umol_graph_core::{
    AutomorphismAlgorithm, AutomorphismOutput, Correspondence, EdgeId, FactorOrdering, Graph,
    NodeId, ParticipantPosition, SubdivisionNodeSource, Unordered,
};
use umol_perm::{Orientation, Permutation};

use super::atom::{AtomForm, ElementForm, IsotopeMassForm};
use super::bond::BondForm;
use super::boolean::BooleanForm;
use super::constraint::{
    AromaticSystemConstraintForm, AromaticValenceForm, AtomConstraintForm, BondConstraintForm,
    Constraint, DativeBondConstraintForm, FluxionalityForm, LigandPermutation, LigandSymmetryForm,
    MoleculeConstraint, MulticenterBondConstraintForm, MulticenterValenceForm,
    NoncovalentBondConstraintForm, OrientedLigandPermutation, RelationalConstraint,
    RingMembershipForm, RingScope, StereoAtomConstraintForm, StereoAtomConstraintsForm,
    StereoBondConstraintForm, StereoBondConstraintsForm, StereoLigandPair, StereogenicityForm,
    TopicityForm, TopicityRelationForm,
};
use super::correspondence::MoleculeCorrespondence;
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
use super::reaction_span::ReactionSpanIntegrityError;
use super::spin::UnpairedElectronsForm;
use super::stereo::{
    CisTransStereoForm, StereoAtomForm, StereoBondForm, StereoConfigurationForm, StereoCoset,
    StereoKind, StereoTerm, Stereogenicity, TetrahedralStereoForm, Topicity,
};
use super::traits::Normalize;
use super::validate::ReactionIntegrityError;

/// Semantic and operational inputs to aggregate canonicalization.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CanonicalizationContext {
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
pub enum CanonicalizationLevel {
    Topology,
    Constitution,
    Structure,
    Full,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct StructuralDomainPosition(u16);

impl StructuralDomainPosition {
    const TOPOLOGY: Self = Self(0);
    const NON_STEREO: Self = Self(1);
    const STEREO: Self = Self(2);
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
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

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct FieldPosition(u16);

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct VariantPosition(u16);

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
#[allow(
    dead_code,
    reason = "reserved by the frozen schema for S11 reaction-span canonicalization"
)]
struct SpanTagPosition(u16);

#[allow(
    dead_code,
    reason = "reserved by the frozen schema for S11 reaction-span canonicalization"
)]
impl SpanTagPosition {
    const UNCHANGED: Self = Self(0);
    const ADDED: Self = Self(1);
    const REMOVED: Self = Self(2);
    const MODIFIED: Self = Self(3);
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
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

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct RelationalConstraintPosition {
    entity: EntityBlockPosition,
    slot: u16,
}

impl RelationalConstraintPosition {
    const fn new(entity: EntityBlockPosition, slot: u16) -> Self {
        Self { entity, slot }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
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

#[derive(Clone, Debug, PartialEq, Eq)]
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

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(
    dead_code,
    reason = "reserved by the frozen schema for S11 reaction-span canonicalization"
)]
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

#[derive(Clone, Debug, PartialEq, Eq)]
enum CanonicalKeyValue {
    #[allow(
        dead_code,
        reason = "reserved by the frozen schema for empty reaction-span values"
    )]
    Unit,
    Bool(bool),
    Unsigned(u64),
    Signed(i64),
    Text(String),
    Sequence(Vec<Self>),
    Product(Vec<FieldKey>),
    Variant(VariantKey),
    #[allow(
        dead_code,
        reason = "reserved by the frozen schema for S11 reaction-span canonicalization"
    )]
    Span(SpanKey),
}

impl CanonicalKeyValue {
    fn position(&self) -> u16 {
        match self {
            Self::Unit => 0,
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
                (Self::Unit, Self::Unit) => Ordering::Equal,
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

#[derive(Clone, Debug, PartialEq, Eq)]
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
        Incidence::MulticenterParticipant(value) => {
            variant(4, [num_form_key(value.normalized()?.as_ref())])
        }
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
    #[cfg(test)]
    edge_sources: Vec<EdgeId>,
    #[cfg(test)]
    incidence_edges: Vec<Vec<EdgeId>>,
    source_node_count: usize,
    #[cfg(test)]
    entity_blocks: Vec<Vec<NodeId>>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ProjectedAutomorphismOutput {
    orbits: Vec<NodeId>,
    // Backend canonical labels are branch-order hints, not the canonical molecule numbering.
    canonical_labels: Vec<NodeId>,
    generators: Vec<Vec<NodeId>>,
}

impl AutomorphismAdapter {
    fn new(incidence_graph: &IncidenceGraph, initial_classes: &InitialClasses) -> Self {
        let source = incidence_graph.graph();
        debug_assert_eq!(initial_classes.entities.len(), source.node_count());
        debug_assert_eq!(initial_classes.incidences.len(), source.edge_count());

        #[cfg(test)]
        let entity_blocks = {
            let mut entity_blocks = Vec::<Vec<NodeId>>::new();
            let mut previous_kind = None;
            for node in source.node_ids() {
                let kind = incidence_graph.entity(node).kind();
                if previous_kind != Some(kind) {
                    entity_blocks.push(Vec::new());
                    previous_kind = Some(kind);
                }
                entity_blocks
                    .last_mut()
                    .expect("current entity block is present")
                    .push(node);
            }
            entity_blocks
        };
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
        #[cfg(test)]
        let mut edge_sources = Vec::new();
        #[cfg(test)]
        let mut incidence_edges = vec![Vec::new(); source.edge_count()];
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

        let mut push_edge = |endpoints: [NodeId; 2], _source_edge: EdgeId| {
            #[cfg(test)]
            let edge = EdgeId(edges.len() as u32);
            edges.push([endpoints[0].0, endpoints[1].0]);
            #[cfg(test)]
            {
                edge_sources.push(_source_edge);
                incidence_edges[_source_edge.index()].push(edge);
            }
        };
        for edge in source.edge_ids() {
            let endpoints = source.edge_endpoints(edge);
            let direct = matches!(
                incidence_graph.incidence(edge),
                Incidence::BondEndpoint | Incidence::NoncovalentEndpoint
            ) && direct_pair_counts[&endpoints] == 1;
            if direct {
                push_edge(endpoints, edge);
                continue;
            }

            let occurrence = NodeId(node_sources.len() as u32);
            node_sources.push(SubdivisionNodeSource::Edge(edge));
            incidence_nodes[edge.index()] = Some(occurrence);
            classes.push(AutomorphismClass::Incidence(
                initial_classes.incidences[edge.index()],
            ));
            push_edge([endpoints[0], occurrence], edge);
            push_edge([occurrence, endpoints[1]], edge);
        }
        let graph = Graph::new(node_sources.len(), &edges);
        debug_assert!(graph.is_simple());

        Self {
            graph,
            classes,
            node_sources,
            incidence_nodes,
            #[cfg(test)]
            edge_sources,
            #[cfg(test)]
            incidence_edges,
            source_node_count: source.node_count(),
            #[cfg(test)]
            entity_blocks,
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

    #[cfg(test)]
    fn edge_source(&self, edge: EdgeId) -> EdgeId {
        self.edge_sources[edge.index()]
    }

    #[cfg(test)]
    fn incidence_edges_of(&self, edge: EdgeId) -> &[EdgeId] {
        &self.incidence_edges[edge.index()]
    }

    #[cfg(test)]
    fn automorphisms(&self, algorithm: AutomorphismAlgorithm) -> ProjectedAutomorphismOutput {
        let output = self
            .graph()
            .automorphisms(|node| self.class(node), algorithm);

        self.project_automorphisms(&output)
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

    #[cfg(test)]
    fn fixed_entity_prefix(&self, entity_count: usize) -> Vec<NodeId> {
        self.cells
            .iter()
            .take_while(|cell| cell.len() == 1)
            .flatten()
            .copied()
            .take_while(|node| node.index() < entity_count)
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
enum BranchOrder {
    #[cfg(test)]
    Node,
    #[cfg(test)]
    ReverseNode,
    BackendCanonical,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct CanonicalSearchOptions {
    automorphism_pruning: bool,
    prefix_pruning: bool,
    branch_order: BranchOrder,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct CanonicalSearchStats {
    refinement_calls: usize,
    visited_leaves: usize,
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
        if best.as_ref().is_none_or(|best| candidate.key < best.key) {
            *best = Some(candidate);
        }
        return;
    };

    let mut candidates = partition.cells[cell_index].clone();
    let automorphisms = (options.automorphism_pruning
        || options.branch_order == BranchOrder::BackendCanonical)
        .then(|| adapter.automorphisms_for_partition(&partition, algorithm));

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

    match options.branch_order {
        #[cfg(test)]
        BranchOrder::Node => candidates.sort_unstable(),
        #[cfg(test)]
        BranchOrder::ReverseNode => candidates.sort_unstable_by(|lhs, rhs| rhs.cmp(lhs)),
        BranchOrder::BackendCanonical => {
            let labels = &automorphisms
                .as_ref()
                .expect("automorphisms requested for backend branch order")
                .canonical_labels;
            let mut ranks = vec![0; adapter.source_node_count];
            for (rank, node) in labels.iter().enumerate() {
                ranks[node.index()] = rank;
            }
            candidates.sort_unstable_by_key(|node| ranks[node.index()]);
        }
    }

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

#[cfg(test)]
fn exhaustive_minimum<K, LeafCandidate>(
    adapter: &AutomorphismAdapter,
    leaf_candidate: &LeafCandidate,
) -> CanonicalCandidate<K>
where
    K: Ord,
    LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
{
    fn visit_cells<K, LeafCandidate>(
        cells: &mut [Vec<NodeId>],
        cell_index: usize,
        order: &mut Vec<NodeId>,
        leaf_candidate: &LeafCandidate,
        best: &mut Option<CanonicalCandidate<K>>,
    ) where
        K: Ord,
        LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
    {
        fn visit_permutations<K, LeafCandidate>(
            cells: &mut [Vec<NodeId>],
            cell_index: usize,
            position: usize,
            order: &mut Vec<NodeId>,
            leaf_candidate: &LeafCandidate,
            best: &mut Option<CanonicalCandidate<K>>,
        ) where
            K: Ord,
            LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
        {
            if position == cells[cell_index].len() {
                let old_len = order.len();
                order.extend_from_slice(&cells[cell_index]);
                visit_cells(cells, cell_index + 1, order, leaf_candidate, best);
                order.truncate(old_len);
                return;
            }

            for next in position..cells[cell_index].len() {
                cells[cell_index].swap(position, next);
                visit_permutations(cells, cell_index, position + 1, order, leaf_candidate, best);
                cells[cell_index].swap(position, next);
            }
        }

        if cell_index == cells.len() {
            let candidate = leaf_candidate(order);
            if best.as_ref().is_none_or(|best| candidate.key < best.key) {
                *best = Some(candidate);
            }
            return;
        }

        visit_permutations(cells, cell_index, 0, order, leaf_candidate, best);
    }

    let mut cells = adapter.entity_blocks.clone();
    let mut best = None;
    visit_cells(
        &mut cells,
        0,
        &mut Vec::with_capacity(adapter.source_node_count),
        leaf_candidate,
        &mut best,
    );

    best.expect("every finite partition has an entity ordering")
}

#[cfg(test)]
fn initial_classes(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
) -> Result<InitialClasses, Contradiction> {
    let (entity_keys, incidence_keys) = initial_class_keys(molecule, incidence_graph)?;
    Ok(rank_initial_classes(&entity_keys, &incidence_keys))
}

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
                    Incidence::AromaticParticipant(_) => variant(3, []),
                    Incidence::MulticenterParticipant(_) => variant(4, []),
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

    if !para_stereo {
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
    }
}

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

#[cfg(test)]
fn topology_comparison_key(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    order: &[NodeId],
) -> Result<CanonicalComparisonKey, Contradiction> {
    Ok(topology_candidate(molecule, incidence_graph, order)?.key)
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

#[cfg(test)]
fn constitution_comparison_key(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    order: &[NodeId],
) -> Result<CanonicalComparisonKey, Contradiction> {
    Ok(constitution_candidate(molecule, incidence_graph, order)?.key)
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

#[cfg(test)]
fn structure_comparison_key(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    order: &[NodeId],
) -> Result<CanonicalComparisonKey, Contradiction> {
    Ok(structure_candidate(molecule, incidence_graph, order)?.key)
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

#[cfg(test)]
#[derive(Clone, Debug, PartialEq, Eq)]
struct MinimumStereoFrames {
    configuration: StereoConfigurationForm,
    permutations: Vec<Permutation>,
}

/// Find every realizable kinded frame action that carries `before` to `after`, retaining all
/// actions whose normalized configuration is minimal. Retaining ties is essential: constraints
/// do not select the structural frame and are considered only by the later constraint-placement
/// search.
#[cfg(test)]
fn minimum_kinded_stereo_frames(
    configuration: &StereoConfigurationForm,
    before: &[StereoLigand],
    after: &[StereoLigand],
) -> Result<Option<MinimumStereoFrames>, Contradiction> {
    let Some(kind) = configuration.kind() else {
        return Ok(None);
    };
    if before.len() != kind.degree() || after.len() != kind.degree() {
        return Ok(None);
    }

    let mut minimum: Option<StereoConfigurationForm> = None;
    let mut permutations = Vec::new();
    for permutation in Permutation::between_all(before, after)
        .into_iter()
        .filter(|permutation| kind.class_key().space().reindex(0, *permutation).is_some())
    {
        let candidate = configuration.apply(permutation).normalize()?;
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

    Ok(minimum.map(|configuration| MinimumStereoFrames {
        configuration,
        permutations,
    }))
}

fn reframe_ligand_permutation(
    permutation: LigandPermutation,
    frame: Permutation,
) -> LigandPermutation {
    LigandPermutation(frame.inverse().compose(permutation.0).compose(frame))
}

fn reframe_ligand_permutation_by_order(
    permutation: LigandPermutation,
    order: &[ParticipantPosition],
) -> Option<LigandPermutation> {
    let frame = permutation_from_position_order(order)?;
    (permutation.0.degree() == frame.degree())
        .then(|| reframe_ligand_permutation(permutation, frame))
}

fn reframe_stereo_ligand_pair_by_order(
    pair: StereoLigandPair,
    inverse: &[ParticipantPosition],
) -> Option<StereoLigandPair> {
    let first = inverse.get(pair.first().index())?.index();
    let second = inverse.get(pair.second().index())?.index();
    Some(StereoLigandPair::new(first.into(), second.into()))
}

fn reframe_stereo_atom_constraint_by_order(
    constraint: StereoAtomConstraintForm,
    order: &[ParticipantPosition],
    inverse: &[ParticipantPosition],
) -> Option<StereoAtomConstraintForm> {
    Some(match constraint {
        StereoAtomConstraintForm::LigandSymmetry(symmetry) => {
            StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm {
                permutation: OrientedLigandPermutation {
                    permutation: reframe_ligand_permutation_by_order(
                        symmetry.permutation.permutation,
                        order,
                    )?,
                    orientation: symmetry.permutation.orientation,
                },
                invariant: symmetry.invariant,
            })
        }
        StereoAtomConstraintForm::Fluxionality(fluxionality) => {
            StereoAtomConstraintForm::Fluxionality(FluxionalityForm {
                permutation: reframe_ligand_permutation_by_order(fluxionality.permutation, order)?,
                active: fluxionality.active,
            })
        }
        StereoAtomConstraintForm::Topicity(topicity) => {
            StereoAtomConstraintForm::Topicity(TopicityForm {
                pair: reframe_stereo_ligand_pair_by_order(topicity.pair, inverse)?,
                relation: topicity.relation,
            })
        }
        StereoAtomConstraintForm::Stereogenicity(stereogenicity) => {
            StereoAtomConstraintForm::Stereogenicity(stereogenicity)
        }
    })
}

fn reframe_stereo_bond_constraint_by_order(
    constraint: StereoBondConstraintForm,
    order: &[ParticipantPosition],
    inverse: &[ParticipantPosition],
) -> Option<StereoBondConstraintForm> {
    Some(match constraint {
        StereoBondConstraintForm::LigandSymmetry(symmetry) => {
            StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm {
                permutation: OrientedLigandPermutation {
                    permutation: reframe_ligand_permutation_by_order(
                        symmetry.permutation.permutation,
                        order,
                    )?,
                    orientation: symmetry.permutation.orientation,
                },
                invariant: symmetry.invariant,
            })
        }
        StereoBondConstraintForm::Fluxionality(fluxionality) => {
            StereoBondConstraintForm::Fluxionality(FluxionalityForm {
                permutation: reframe_ligand_permutation_by_order(fluxionality.permutation, order)?,
                active: fluxionality.active,
            })
        }
        StereoBondConstraintForm::Topicity(topicity) => {
            StereoBondConstraintForm::Topicity(TopicityForm {
                pair: reframe_stereo_ligand_pair_by_order(topicity.pair, inverse)?,
                relation: topicity.relation,
            })
        }
        StereoBondConstraintForm::Stereogenicity(stereogenicity) => {
            StereoBondConstraintForm::Stereogenicity(stereogenicity)
        }
    })
}

fn reframe_stereo_atom_form_by_order(
    form: &StereoAtomForm,
    order: &[ParticipantPosition],
) -> Option<StereoAtomForm> {
    let inverse = inverse_position_order(order)?;
    let configuration = match &form.configuration {
        StereoConfigurationForm::Undetermined => StereoConfigurationForm::Undetermined,
        StereoConfigurationForm::Kinded(kind, _) => {
            let frame = permutation_from_position_order(order)?;
            kind.class_key().space().reindex(0, frame)?;
            form.configuration.apply(frame)
        }
    };
    let mut constraints = StereoAtomConstraintsForm::default();
    for constraint in form.constraints.iter().cloned() {
        constraints.set(reframe_stereo_atom_constraint_by_order(
            constraint, order, &inverse,
        )?);
    }
    Some(StereoAtomForm {
        configuration,
        constraints,
    })
}

fn reframe_stereo_bond_form_by_order(
    form: &StereoBondForm,
    order: &[ParticipantPosition],
) -> Option<StereoBondForm> {
    let inverse = inverse_position_order(order)?;
    let configuration = match &form.configuration {
        StereoConfigurationForm::Undetermined => StereoConfigurationForm::Undetermined,
        StereoConfigurationForm::Kinded(kind, _) => {
            let frame = permutation_from_position_order(order)?;
            kind.class_key().space().reindex(0, frame)?;
            form.configuration.apply(frame)
        }
    };
    let mut constraints = StereoBondConstraintsForm::default();
    for constraint in form.constraints.iter().cloned() {
        constraints.set(reframe_stereo_bond_constraint_by_order(
            constraint, order, &inverse,
        )?);
    }
    Some(StereoBondForm {
        configuration,
        constraints,
    })
}

fn reframe_molecule_constraint_by_order(
    constraint: Constraint,
    entity: Entity,
    order: &[ParticipantPosition],
    inverse: &[ParticipantPosition],
) -> Option<Constraint> {
    Some(match constraint {
        Constraint::StereoAtom(id, kind, constraint) if entity == Entity::StereoAtom(id) => {
            Constraint::StereoAtom(
                id,
                kind,
                reframe_stereo_atom_constraint_by_order(constraint, order, inverse)?,
            )
        }
        Constraint::StereoBond(id, kind, constraint) if entity == Entity::StereoBond(id) => {
            Constraint::StereoBond(
                id,
                kind,
                reframe_stereo_bond_constraint_by_order(constraint, order, inverse)?,
            )
        }
        Constraint::And(constraints) => Constraint::And(
            constraints
                .into_iter()
                .map(|constraint| {
                    reframe_molecule_constraint_by_order(constraint, entity, order, inverse)
                })
                .collect::<Option<Vec<_>>>()?,
        ),
        Constraint::Or(constraints) => Constraint::Or(
            constraints
                .into_iter()
                .map(|constraint| {
                    reframe_molecule_constraint_by_order(constraint, entity, order, inverse)
                })
                .collect::<Option<Vec<_>>>()?,
        ),
        Constraint::Not(constraint) => Constraint::Not(Box::new(
            reframe_molecule_constraint_by_order(*constraint, entity, order, inverse)?,
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
    let inverse = inverse_position_order(order)?;
    entries.constraints = entries
        .constraints
        .into_iter()
        .map(|constraint| {
            reframe_molecule_constraint_by_order(
                constraint,
                Entity::StereoAtom(id),
                order,
                &inverse,
            )
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
    let inverse = inverse_position_order(order)?;
    entries.constraints = entries
        .constraints
        .into_iter()
        .map(|constraint| {
            reframe_molecule_constraint_by_order(
                constraint,
                Entity::StereoBond(id),
                order,
                &inverse,
            )
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

#[allow(dead_code, reason = "wired into the public Canonicalize trait in S9c")]
fn canonicalize_topology(
    molecule: &Molecule,
    context: &CanonicalizationContext,
) -> Result<Molecule, MoleculeCanonicalizationError> {
    canonicalize_topology_with_options(
        molecule,
        context,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            prefix_pruning: false,
            branch_order: BranchOrder::BackendCanonical,
        },
    )
    .map(|(molecule, _)| molecule)
}

fn canonicalize_topology_with_options(
    molecule: &Molecule,
    context: &CanonicalizationContext,
    options: CanonicalSearchOptions,
) -> Result<(Molecule, MoleculeCorrespondence), MoleculeCanonicalizationError> {
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

#[allow(dead_code, reason = "wired into the public Canonicalize trait in S9c")]
fn canonicalize_constitution(
    molecule: &Molecule,
    context: &CanonicalizationContext,
) -> Result<Molecule, MoleculeCanonicalizationError> {
    canonicalize_constitution_with_options(
        molecule,
        context,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            prefix_pruning: false,
            branch_order: BranchOrder::BackendCanonical,
        },
    )
    .map(|(molecule, _)| molecule)
}

fn canonicalize_constitution_with_options(
    molecule: &Molecule,
    context: &CanonicalizationContext,
    options: CanonicalSearchOptions,
) -> Result<(Molecule, MoleculeCorrespondence), MoleculeCanonicalizationError> {
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

#[allow(dead_code, reason = "wired into the public Canonicalize trait in S9c")]
fn canonicalize_structure(
    molecule: &Molecule,
    context: &CanonicalizationContext,
) -> Result<Molecule, MoleculeCanonicalizationError> {
    canonicalize_structure_with_options(
        molecule,
        context,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            prefix_pruning: false,
            branch_order: BranchOrder::BackendCanonical,
        },
    )
    .map(|(molecule, _)| molecule)
}

fn canonicalize_structure_with_options(
    molecule: &Molecule,
    context: &CanonicalizationContext,
    mut options: CanonicalSearchOptions,
) -> Result<(Molecule, MoleculeCorrespondence), MoleculeCanonicalizationError> {
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

#[allow(dead_code, reason = "wired into the public Canonicalize trait in S9c")]
fn canonicalize_full(
    molecule: &Molecule,
    context: &CanonicalizationContext,
) -> Result<Molecule, MoleculeCanonicalizationError> {
    canonicalize_full_with_options(
        molecule,
        context,
        CanonicalSearchOptions {
            automorphism_pruning: false,
            prefix_pruning: false,
            branch_order: BranchOrder::BackendCanonical,
        },
    )
}

fn canonicalize_full_with_options(
    molecule: &Molecule,
    context: &CanonicalizationContext,
    options: CanonicalSearchOptions,
) -> Result<Molecule, MoleculeCanonicalizationError> {
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
    let (_, canonical) =
        complete_candidate(molecule, &incidence_graph, &selected.candidate.entity_order)?;
    Ok(canonical)
}

/// Failure to construct a canonical [`Molecule`](super::Molecule).
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum MoleculeCanonicalizationError {
    /// The molecule does not satisfy its representation-integrity contract.
    #[error(transparent)]
    Integrity(#[from] MoleculeIntegrityError),
    /// Intrinsic normalization of a carried value reached a contradiction.
    #[error(transparent)]
    Contradiction(#[from] Contradiction),
}

/// Failure to construct a canonical [`ReactionSpan`](super::ReactionSpan).
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ReactionSpanCanonicalizationError {
    /// The reaction span does not satisfy its representation-integrity contract.
    #[error(transparent)]
    Integrity(#[from] ReactionSpanIntegrityError),
    /// Intrinsic normalization of a carried value reached a contradiction.
    #[error(transparent)]
    Contradiction(#[from] Contradiction),
}

/// Failure to construct a canonical [`Reaction`](super::Reaction).
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ReactionCanonicalizationError {
    /// The reaction does not satisfy its representation-integrity contract.
    #[error(transparent)]
    Integrity(#[from] ReactionIntegrityError),
    /// Intrinsic normalization or span materialization reached a contradiction.
    #[error(transparent)]
    Contradiction(#[from] Contradiction),
}

#[cfg(test)]
mod tests {
    use std::array;
    use std::cmp::Ordering;
    use std::hint::black_box;
    use std::iter::once;
    use std::time::{Duration, Instant};

    use rstest::{fixture, rstest};
    use umol_chem::element::Element;
    use umol_perm::{Orientation, Permutation};

    use super::*;
    use crate::ir::{
        AromaticSystemForm, AromaticSystemId, AtomConstraintForm, AtomForm, AtomId, BondForm,
        BondId, BooleanForm, Constraint, Constraints, DativeBondForm, DativeBondId, Entity,
        FluxionalityForm, IncidenceLevel, LigandPermutation, LigandSymmetryForm,
        MoleculeCorrespondence, MoleculeEntries, MulticenterBondForm, MulticenterBondId,
        NoncovalentBondForm, NoncovalentBondId, OrientedLigandPermutation,
        StereoAtomConstraintForm, StereoAtomForm, StereoAtomId, StereoBondConstraintForm,
        StereoBondForm, StereoBondId, StereoConfigurationForm, StereoCoset, StereoKind,
        StereoLigand, StereoLigandPair, StereoTerm, Stereogenicity, StereogenicityForm, Topicity,
        TopicityForm, TopicityRelationForm,
    };

    mod benchmark_cases {
        include!("../../benches/canonicalization_cases.rs");
    }

    #[fixture]
    fn initial_class_molecule() -> Molecule {
        let normalized_three = NumForm::ArithExpr(Box::new(ArithExpr::Sum(vec![
            ArithExpr::Lit(1),
            ArithExpr::Lit(2),
        ])));
        let normalized_one = NumForm::ArithExpr(Box::new(ArithExpr::Sum(vec![
            ArithExpr::Lit(0),
            ArithExpr::Lit(1),
        ])));
        let bond_ligands = vec![
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
        ];

        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C).with_charge(3_i64),
                AtomForm::from_element(Element::C)
                    .with_charge(normalized_three)
                    .with_constraint(AtomConstraintForm::Valence(NumForm::Lit(4))),
                AtomForm::from_element(Element::O).with_charge(3_i64),
                AtomForm::from_element(Element::C).with_charge(4_i64),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::H),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::new(normalized_one.clone())),
                (AtomId(2), AtomId(3), BondForm::from_order(2)),
            ],
            dative: vec![
                (vec![AtomId(0)], AtomId(4), DativeBondForm::from_order(1)),
                (
                    vec![AtomId(1)],
                    AtomId(4),
                    DativeBondForm::new(normalized_one),
                ),
                (vec![AtomId(2)], AtomId(4), DativeBondForm::from_order(2)),
            ],
            aromatic: vec![
                (
                    vec![AtomId(0), AtomId(1)],
                    AromaticSystemForm::from_electrons(vec![1, 2]),
                ),
                (
                    vec![AtomId(2), AtomId(3)],
                    AromaticSystemForm::from_electrons(vec![2, 1]),
                ),
                (
                    vec![AtomId(4), AtomId(5)],
                    AromaticSystemForm::from_electrons(vec![1, 2]).with_charge(1_i64),
                ),
            ],
            multicenter: vec![
                (
                    vec![AtomId(0), AtomId(2)],
                    MulticenterBondForm::from_electrons(vec![1, 2]),
                ),
                (
                    vec![AtomId(1), AtomId(3)],
                    MulticenterBondForm::from_electrons(vec![2, 1]),
                ),
            ],
            noncovalent: vec![
                (
                    AtomId(0),
                    AtomId(5),
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                ),
                (
                    AtomId(1),
                    AtomId(5),
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                ),
                (
                    AtomId(2),
                    AtomId(5),
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::Ionic),
                ),
            ],
            stereo_atoms: vec![
                (
                    AtomId(0),
                    vec![
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    ],
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                ),
                (
                    AtomId(1),
                    vec![
                        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    ],
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                ),
                (
                    AtomId(2),
                    vec![
                        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    ],
                    StereoAtomForm::new(StereoKind::SquarePlanar, StereoCoset::Lit(0)),
                ),
            ],
            stereo_bonds: vec![
                (
                    BondId(0),
                    bond_ligands.clone(),
                    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
                ),
                (
                    BondId(1),
                    bond_ligands,
                    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                ),
            ],
            ..Default::default()
        })
    }

    #[fixture]
    fn canonicalization_context() -> CanonicalizationContext {
        CanonicalizationContext {
            para_stereo: false,
            automorphism_algorithm: AutomorphismAlgorithm::Nauty,
        }
    }

    #[fixture]
    fn stereo_atom_canonicalization_molecule() -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
                AtomForm::from_element(Element::Br),
                AtomForm::from_element(Element::I),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
                (AtomId(0), AtomId(4), BondForm::from_order(1)),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            )],
            ..Default::default()
        })
    }

    #[fixture]
    fn stereo_bond_canonicalization_molecule() -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
                AtomForm::from_element(Element::Br),
                AtomForm::from_element(Element::I),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(2)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
                (AtomId(1), AtomId(4), BondForm::from_order(1)),
                (AtomId(1), AtomId(5), BondForm::from_order(1)),
            ],
            stereo_bonds: vec![(
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoBondForm::new(StereoKind::CisTrans, 1u32),
            )],
            ..Default::default()
        })
    }

    #[fixture]
    fn symmetric_stereo_canonicalization_molecule() -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 5],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
            )],
            ..Default::default()
        })
    }

    #[fixture]
    fn meso_canonicalization_molecule() -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::Cl),
                AtomForm::from_element(Element::Cl),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(4), BondForm::from_order(1)),
                (AtomId(1), AtomId(3), BondForm::from_order(1)),
                (AtomId(1), AtomId(5), BondForm::from_order(1)),
            ],
            stereo_atoms: vec![
                (
                    AtomId(0),
                    vec![
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    ],
                    StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
                ),
                (
                    AtomId(1),
                    vec![
                        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    ],
                    StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
                ),
            ],
            ..Default::default()
        })
    }

    #[fixture]
    fn repeated_ligand_canonicalization_molecule() -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            )],
            ..Default::default()
        })
    }

    #[fixture]
    fn para_stereo_canonicalization_molecule() -> Molecule {
        let outer_ligands = vec![
            StereoLigand::new(AtomId(10), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(11), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(12), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(13), StereoLigandKind::Atom),
        ];

        Molecule::from_entries(MoleculeEntries {
            atoms: [
                Element::C,
                Element::C,
                Element::C,
                Element::C,
                Element::C,
                Element::C,
                Element::C,
                Element::C,
                Element::C,
                Element::C,
                Element::F,
                Element::Cl,
                Element::Br,
                Element::I,
            ]
            .into_iter()
            .map(AtomForm::from_element)
            .collect(),
            stereo_atoms: vec![
                (
                    AtomId(0),
                    [2, 3, 4, 5]
                        .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                        .into(),
                    StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
                ),
                (
                    AtomId(1),
                    [6, 8, 7, 9]
                        .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                        .into(),
                    StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
                ),
                (
                    AtomId(2),
                    outer_ligands.clone(),
                    StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
                ),
                (
                    AtomId(3),
                    outer_ligands.clone(),
                    StereoAtomForm::new(StereoKind::CisTrans, 0u32),
                ),
                (
                    AtomId(4),
                    outer_ligands.clone(),
                    StereoAtomForm::new(StereoKind::Axial, 0u32),
                ),
                (
                    AtomId(5),
                    outer_ligands.clone(),
                    StereoAtomForm::new(StereoKind::SquarePlanar, 0u32),
                ),
                (
                    AtomId(6),
                    outer_ligands.clone(),
                    StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
                ),
                (
                    AtomId(7),
                    outer_ligands.clone(),
                    StereoAtomForm::new(StereoKind::CisTrans, 0u32),
                ),
                (
                    AtomId(8),
                    outer_ligands.clone(),
                    StereoAtomForm::new(StereoKind::Axial, 0u32),
                ),
                (
                    AtomId(9),
                    outer_ligands,
                    StereoAtomForm::new(StereoKind::SquarePlanar, 0u32),
                ),
            ],
            ..Default::default()
        })
    }

    fn selected_structure_key(molecule: &Molecule) -> CanonicalComparisonKey {
        let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
        let order = incidence_graph.graph().node_ids().collect::<Vec<_>>();
        structure_comparison_key(molecule, &incidence_graph, &order)
            .expect("canonical molecule has a structure comparison key")
    }

    #[rstest]
    #[case::tetrahedral(StereoKind::Tetrahedral)]
    #[case::axial(StereoKind::Axial)]
    #[case::square_planar(StereoKind::SquarePlanar)]
    #[case::trigonal_bipyramidal(StereoKind::TrigonalBipyramidal)]
    #[case::octahedral(StereoKind::Octahedral)]
    fn test_reframe_stereo_atom(#[case] kind: StereoKind) {
        let degree = kind.degree();
        let swap = |first, second| {
            let mut image = (0..degree).collect::<Vec<_>>();
            image.swap(first, second);
            Permutation::from_image(&image)
        };
        let (frame, next_frame, expected_permutation, expected_pair) = if kind == StereoKind::Axial
        {
            (
                Permutation::from_image(&[2, 3, 0, 1]),
                Permutation::from_image(&[1, 0, 2, 3]),
                LigandPermutation(swap(2, 3)),
                StereoLigandPair::new(2usize.into(), 3usize.into()),
            )
        } else {
            (
                Permutation::from_image(&(1..degree).chain(once(0)).collect::<Vec<_>>()),
                swap(1, 2),
                LigandPermutation(swap(0, degree - 1)),
                StereoLigandPair::new(0usize.into(), (degree - 1).into()),
            )
        };
        let source_permutation = LigandPermutation(swap(0, 1));
        let source_pair = StereoLigandPair::new(0usize.into(), 1usize.into());
        let source_constraints = vec![
            StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm {
                permutation: OrientedLigandPermutation {
                    permutation: source_permutation,
                    orientation: Orientation::Improper,
                },
                invariant: BooleanForm::Lit(true),
            }),
            StereoAtomConstraintForm::Fluxionality(FluxionalityForm {
                permutation: source_permutation,
                active: BooleanForm::Lit(true),
            }),
            StereoAtomConstraintForm::Topicity(TopicityForm {
                pair: source_pair,
                relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
            }),
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
        ];
        let expected_constraints = vec![
            StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm {
                permutation: OrientedLigandPermutation {
                    permutation: expected_permutation,
                    orientation: Orientation::Improper,
                },
                invariant: BooleanForm::Lit(true),
            }),
            StereoAtomConstraintForm::Fluxionality(FluxionalityForm {
                permutation: expected_permutation,
                active: BooleanForm::Lit(true),
            }),
            StereoAtomConstraintForm::Topicity(TopicityForm {
                pair: expected_pair,
                relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
            }),
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
        ];
        let global_constraints = |constraints: &[StereoAtomConstraintForm]| {
            Constraints::from(Constraint::And(vec![
                Constraint::Not(Box::new(Constraint::StereoAtom(
                    StereoAtomId(0),
                    kind,
                    constraints[0].clone(),
                ))),
                Constraint::Or(
                    constraints[1..]
                        .iter()
                        .cloned()
                        .map(|constraint| Constraint::StereoAtom(StereoAtomId(0), kind, constraint))
                        .collect(),
                ),
            ]))
        };
        let atoms = (0..=degree)
            .map(|_| AtomForm::from_element(Element::C))
            .collect::<Vec<_>>();
        let ligands = (1..=degree)
            .map(|atom| StereoLigand::new(AtomId(atom as u32), StereoLigandKind::Atom))
            .collect::<Vec<_>>();
        let source_form = StereoAtomForm {
            configuration: StereoConfigurationForm::kinded(kind, StereoCoset::Lit(0)),
            constraints: source_constraints.clone().into(),
        };
        let source = Molecule::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
            stereo_atoms: vec![(AtomId(0), ligands.clone(), source_form.clone())],
            constraints: global_constraints(&source_constraints),
            ..Default::default()
        });
        let expected = Molecule::from_entries(MoleculeEntries {
            atoms,
            stereo_atoms: vec![(
                AtomId(0),
                frame.act(&ligands),
                StereoAtomForm {
                    configuration: source_form.configuration.apply(frame),
                    constraints: expected_constraints.clone().into(),
                },
            )],
            constraints: global_constraints(&expected_constraints),
            ..Default::default()
        });

        let reframed = reframe_stereo_atom(&source, StereoAtomId(0), frame);
        assert_eq!(reframed, expected);
        assert_eq!(
            reframe_stereo_atom(&reframed, StereoAtomId(0), frame.inverse()),
            source
        );
        assert_eq!(
            reframe_stereo_atom(&reframed, StereoAtomId(0), next_frame),
            reframe_stereo_atom(&source, StereoAtomId(0), frame.compose(next_frame))
        );
    }

    #[rstest]
    fn test_reframe_stereo_bond() {
        let degree = StereoKind::CisTrans.degree();
        let swap = |first, second| {
            let mut image = (0..degree).collect::<Vec<_>>();
            image.swap(first, second);
            Permutation::from_image(&image)
        };
        let frame = Permutation::from_image(&[2, 3, 0, 1]);
        let next_frame = Permutation::from_image(&[1, 0, 2, 3]);
        let source_permutation = LigandPermutation(swap(0, 1));
        let expected_permutation = LigandPermutation(swap(2, 3));
        let source_pair = StereoLigandPair::new(0usize.into(), 1usize.into());
        let expected_pair = StereoLigandPair::new(2usize.into(), 3usize.into());
        let source_constraints = vec![
            StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm {
                permutation: OrientedLigandPermutation {
                    permutation: source_permutation,
                    orientation: Orientation::Proper,
                },
                invariant: BooleanForm::Lit(true),
            }),
            StereoBondConstraintForm::Fluxionality(FluxionalityForm {
                permutation: source_permutation,
                active: BooleanForm::Lit(true),
            }),
            StereoBondConstraintForm::Topicity(TopicityForm {
                pair: source_pair,
                relation: TopicityRelationForm::Lit(Topicity::Diastereotopic),
            }),
            StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
        ];
        let expected_constraints = vec![
            StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm {
                permutation: OrientedLigandPermutation {
                    permutation: expected_permutation,
                    orientation: Orientation::Proper,
                },
                invariant: BooleanForm::Lit(true),
            }),
            StereoBondConstraintForm::Fluxionality(FluxionalityForm {
                permutation: expected_permutation,
                active: BooleanForm::Lit(true),
            }),
            StereoBondConstraintForm::Topicity(TopicityForm {
                pair: expected_pair,
                relation: TopicityRelationForm::Lit(Topicity::Diastereotopic),
            }),
            StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
        ];
        let global_constraints = |constraints: &[StereoBondConstraintForm]| {
            Constraints::from(Constraint::And(vec![
                Constraint::Not(Box::new(Constraint::StereoBond(
                    StereoBondId(0),
                    StereoKind::CisTrans,
                    constraints[0].clone(),
                ))),
                Constraint::Or(
                    constraints[1..]
                        .iter()
                        .cloned()
                        .map(|constraint| {
                            Constraint::StereoBond(
                                StereoBondId(0),
                                StereoKind::CisTrans,
                                constraint,
                            )
                        })
                        .collect(),
                ),
            ]))
        };
        let atoms = (0..6)
            .map(|_| AtomForm::from_element(Element::C))
            .collect::<Vec<_>>();
        let ligands = (2..6)
            .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
            .collect::<Vec<_>>();
        let source_form = StereoBondForm {
            configuration: StereoConfigurationForm::kinded(
                StereoKind::CisTrans,
                StereoCoset::Lit(0),
            ),
            constraints: source_constraints.clone().into(),
        };
        let source = Molecule::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            stereo_bonds: vec![(BondId(0), ligands.clone(), source_form.clone())],
            constraints: global_constraints(&source_constraints),
            ..Default::default()
        });
        let expected = Molecule::from_entries(MoleculeEntries {
            atoms,
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            stereo_bonds: vec![(
                BondId(0),
                frame.act(&ligands),
                StereoBondForm {
                    configuration: source_form.configuration.apply(frame),
                    constraints: expected_constraints.clone().into(),
                },
            )],
            constraints: global_constraints(&expected_constraints),
            ..Default::default()
        });

        let reframed = reframe_stereo_bond(&source, StereoBondId(0), frame);
        assert_eq!(reframed, expected);
        assert_eq!(
            reframe_stereo_bond(&reframed, StereoBondId(0), frame.inverse()),
            source
        );
        assert_eq!(
            reframe_stereo_bond(&reframed, StereoBondId(0), next_frame),
            reframe_stereo_bond(&source, StereoBondId(0), frame.compose(next_frame))
        );
    }

    #[rstest]
    #[case::literal(
        StereoCoset::Lit(0),
        StereoCoset::Lit(0),
        vec![Permutation::from_image(&[1, 2, 0, 3])]
    )]
    #[case::undetermined(
        StereoCoset::Undetermined,
        StereoCoset::Undetermined,
        vec![
            Permutation::from_image(&[1, 2, 0, 3]),
            Permutation::from_image(&[2, 1, 0, 3]),
        ]
    )]
    #[case::set_valued(
        StereoCoset::lit_set([0, 1]),
        StereoCoset::lit_set([0, 1]),
        vec![
            Permutation::from_image(&[1, 2, 0, 3]),
            Permutation::from_image(&[2, 1, 0, 3]),
        ]
    )]
    #[case::symbolic(
        StereoCoset::term(StereoTerm::var("x")),
        StereoCoset::term(StereoTerm::var("x")),
        vec![Permutation::from_image(&[1, 2, 0, 3])]
    )]
    fn test_minimum_kinded_stereo_frames(
        #[case] coset: StereoCoset,
        #[case] expected_coset: StereoCoset,
        #[case] expected_permutations: Vec<Permutation>,
    ) {
        let repeated = StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen);
        let before = vec![
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            repeated,
            repeated,
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        ];
        let (after, _) = sort_ligand_frame(&before);
        let result = minimum_kinded_stereo_frames(
            &StereoConfigurationForm::kinded(StereoKind::Tetrahedral, coset),
            &before,
            &after,
        )
        .expect("fixed configurations normalize")
        .expect("the frames contain the same ligand multiset");

        assert_eq!(
            result,
            MinimumStereoFrames {
                configuration: StereoConfigurationForm::kinded(
                    StereoKind::Tetrahedral,
                    expected_coset,
                ),
                permutations: expected_permutations,
            }
        );
    }

    #[rstest]
    fn test_kindless_stereo_atom_frame_order() {
        let ligands = (0..7)
            .rev()
            .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
            .collect::<Vec<_>>();
        let (sorted, order) = sort_ligand_frame(&ligands);
        let source = StereoAtomForm {
            configuration: StereoConfigurationForm::Undetermined,
            constraints: StereoAtomConstraintForm::Topicity(TopicityForm {
                pair: StereoLigandPair::new(0usize.into(), 2usize.into()),
                relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
            })
            .into(),
        };
        let expected = StereoAtomForm {
            configuration: StereoConfigurationForm::Undetermined,
            constraints: StereoAtomConstraintForm::Topicity(TopicityForm {
                pair: StereoLigandPair::new(4usize.into(), 6usize.into()),
                relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
            })
            .into(),
        };

        assert_eq!(
            sorted,
            (0..7)
                .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
                .collect::<Vec<_>>()
        );
        assert_eq!(
            reframe_stereo_atom_form_by_order(&source, &order),
            Some(expected)
        );
    }

    #[rstest]
    fn test_kindless_stereo_bond_frame_order() {
        let ligands = (0..4)
            .rev()
            .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
            .collect::<Vec<_>>();
        let (_, order) = sort_ligand_frame(&ligands);
        let source = StereoBondForm {
            configuration: StereoConfigurationForm::Undetermined,
            constraints: vec![
                StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm {
                    permutation: OrientedLigandPermutation {
                        permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])),
                        orientation: Orientation::Proper,
                    },
                    invariant: BooleanForm::Lit(true),
                }),
                StereoBondConstraintForm::Topicity(TopicityForm {
                    pair: StereoLigandPair::new(0usize.into(), 1usize.into()),
                    relation: TopicityRelationForm::Lit(Topicity::Diastereotopic),
                }),
            ]
            .into(),
        };
        let expected = StereoBondForm {
            configuration: StereoConfigurationForm::Undetermined,
            constraints: vec![
                StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm {
                    permutation: OrientedLigandPermutation {
                        permutation: LigandPermutation(Permutation::from_image(&[0, 1, 3, 2])),
                        orientation: Orientation::Proper,
                    },
                    invariant: BooleanForm::Lit(true),
                }),
                StereoBondConstraintForm::Topicity(TopicityForm {
                    pair: StereoLigandPair::new(2usize.into(), 3usize.into()),
                    relation: TopicityRelationForm::Lit(Topicity::Diastereotopic),
                }),
            ]
            .into(),
        };

        assert_eq!(
            reframe_stereo_bond_form_by_order(&source, &order),
            Some(expected)
        );
    }

    #[rstest]
    #[case::tetrahedral(StereoKind::Tetrahedral)]
    #[case::cis_trans(StereoKind::CisTrans)]
    #[case::axial(StereoKind::Axial)]
    #[case::square_planar(StereoKind::SquarePlanar)]
    #[case::trigonal_bipyramidal(StereoKind::TrigonalBipyramidal)]
    #[case::octahedral(StereoKind::Octahedral)]
    fn test_stereo_refinement_descriptor_frame_invariant(#[case] kind: StereoKind) {
        let degree = kind.degree();
        let frame = if matches!(kind, StereoKind::CisTrans | StereoKind::Axial) {
            Permutation::from_image(&[2, 3, 0, 1])
        } else {
            Permutation::from_image(&(1..degree).chain(once(0)).collect::<Vec<_>>())
        };
        let ligands = (0..degree)
            .map(|class| (class as u32, StereoLigandKind::Atom))
            .collect::<Vec<_>>();
        let configuration = StereoConfigurationForm::kinded(kind, 0u32);

        assert_eq!(
            stereo_refinement_descriptor(7, &ligands, &configuration),
            stereo_refinement_descriptor(7, &frame.act(&ligands), &configuration.apply(frame),),
        );
    }

    #[rstest]
    #[case::one_pass(false, 1)]
    #[case::fixpoint(true, 3)]
    fn test_structure_partition(
        para_stereo_canonicalization_molecule: Molecule,
        #[case] para_stereo: bool,
        #[case] expected_rounds: usize,
    ) {
        let molecule = para_stereo_canonicalization_molecule;
        let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
        let (entity_keys, incidence_keys) = initial_class_keys(&molecule, &incidence_graph)
            .expect("fixed molecule has initial classes");
        let classes = rank_initial_classes(&entity_keys, &incidence_keys);
        let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
        let (_, rounds) = structure_partition(
            &molecule,
            &incidence_graph,
            &adapter,
            &entity_keys,
            para_stereo,
        )
        .expect("fixed molecule has a structure partition");

        assert_eq!(rounds, expected_rounds);
    }

    #[rstest]
    fn test_canonicalize_structure_para_stereo(
        para_stereo_canonicalization_molecule: Molecule,
        canonicalization_context: CanonicalizationContext,
    ) {
        let correspondence = molecule_correspondence(&[
            vec![1, 0, 6, 8, 7, 9, 2, 4, 3, 5, 13, 11, 10, 12],
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            vec![1, 0, 6, 8, 7, 9, 2, 4, 3, 5],
            Vec::new(),
        ]);
        let renumbered = para_stereo_canonicalization_molecule.remap(&correspondence);
        let context = CanonicalizationContext {
            para_stereo: true,
            ..canonicalization_context
        };
        let canonical = canonicalize_structure(&para_stereo_canonicalization_molecule, &context)
            .expect("fixed molecule canonicalizes");

        assert_eq!(
            canonicalize_structure(&renumbered, &context),
            Ok(canonical.clone()),
        );
        assert_eq!(canonicalize_structure(&canonical, &context), Ok(canonical));
    }

    #[rstest]
    fn test_structure_comparison_key(stereo_atom_canonicalization_molecule: Molecule) {
        let mut constrained_entries = molecule_entries(&stereo_atom_canonicalization_molecule);
        constrained_entries.stereo_atoms[0].2.constraints =
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
                Stereogenicity::Stereogenic,
            ))
            .into();
        let constrained = Molecule::from_entries(constrained_entries);
        let incidence_graph =
            stereo_atom_canonicalization_molecule.incidence_graph(IncidenceLevel::Full);
        let constrained_incidence_graph = constrained.incidence_graph(IncidenceLevel::Full);
        let order = incidence_graph.graph().node_ids().collect::<Vec<_>>();

        assert_eq!(
            structure_comparison_key(
                &stereo_atom_canonicalization_molecule,
                &incidence_graph,
                &order,
            ),
            structure_comparison_key(&constrained, &constrained_incidence_graph, &order),
        );
    }

    #[rstest]
    fn test_canonicalize_structure_stereo_atom(
        canonicalization_context: CanonicalizationContext,
        stereo_atom_canonicalization_molecule: Molecule,
    ) {
        let renumbering = molecule_correspondence(&[
            vec![4, 2, 0, 3, 1],
            vec![2, 0, 3, 1],
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            vec![0],
            Vec::new(),
        ]);
        let renumbered = stereo_atom_canonicalization_molecule.remap(&renumbering);
        let canonical = canonicalize_structure(
            &stereo_atom_canonicalization_molecule,
            &canonicalization_context,
        )
        .expect("fixed molecule canonicalizes");

        assert_eq!(
            canonicalize_structure(&renumbered, &canonicalization_context),
            Ok(canonical.clone()),
        );
        assert_eq!(
            canonicalize_structure(&canonical, &canonicalization_context),
            Ok(canonical.clone()),
        );
        assert_eq!(canonical.check_integrity(), Ok(()));
        assert_eq!(canonical.stereo_atoms().count(), 1);
    }

    #[rstest]
    fn test_canonicalize_structure_configuration(
        canonicalization_context: CanonicalizationContext,
        stereo_atom_canonicalization_molecule: Molecule,
    ) {
        let mut opposite_entries = molecule_entries(&stereo_atom_canonicalization_molecule);
        opposite_entries.stereo_atoms[0].2.configuration =
            StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1u32);
        let opposite = Molecule::from_entries(opposite_entries);

        assert_ne!(
            canonicalize_structure(
                &stereo_atom_canonicalization_molecule,
                &canonicalization_context,
            ),
            canonicalize_structure(&opposite, &canonicalization_context),
        );
    }

    #[rstest]
    fn test_canonicalize_structure_stereo_bond(
        canonicalization_context: CanonicalizationContext,
        stereo_bond_canonicalization_molecule: Molecule,
    ) {
        let renumbering = molecule_correspondence(&[
            vec![5, 3, 1, 4, 2, 0],
            vec![4, 2, 0, 3, 1],
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            vec![0],
        ]);
        let renumbered = stereo_bond_canonicalization_molecule.remap(&renumbering);
        let canonical = canonicalize_structure(
            &stereo_bond_canonicalization_molecule,
            &canonicalization_context,
        )
        .expect("fixed molecule canonicalizes");

        assert_eq!(
            canonicalize_structure(&renumbered, &canonicalization_context),
            Ok(canonical.clone()),
        );
        assert_eq!(
            canonicalize_structure(&canonical, &canonicalization_context),
            Ok(canonical.clone()),
        );
        assert_eq!(canonical.check_integrity(), Ok(()));
        assert_eq!(canonical.stereo_bonds().count(), 1);
    }

    #[rstest]
    fn test_canonicalize_structure_stereo_atom_constraints(
        canonicalization_context: CanonicalizationContext,
        stereo_atom_canonicalization_molecule: Molecule,
    ) {
        let constraint = StereoAtomConstraintForm::Topicity(TopicityForm {
            pair: StereoLigandPair::new(0usize.into(), 1usize.into()),
            relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
        });
        let mut entries = molecule_entries(&stereo_atom_canonicalization_molecule);
        entries.stereo_atoms[0].2.constraints = constraint.clone().into();
        entries.constraints =
            Constraint::StereoAtom(StereoAtomId(0), StereoKind::Tetrahedral, constraint).into();
        let source = Molecule::from_entries(entries);
        let reframed = reframe_stereo_atom(
            &source,
            StereoAtomId(0),
            Permutation::from_image(&[1, 0, 2, 3]),
        );

        assert_eq!(
            canonicalize_structure(&reframed, &canonicalization_context),
            canonicalize_structure(&source, &canonicalization_context),
        );
    }

    #[rstest]
    fn test_canonicalize_structure_stereo_bond_constraints(
        canonicalization_context: CanonicalizationContext,
        stereo_bond_canonicalization_molecule: Molecule,
    ) {
        let constraint = StereoBondConstraintForm::Topicity(TopicityForm {
            pair: StereoLigandPair::new(0usize.into(), 1usize.into()),
            relation: TopicityRelationForm::Lit(Topicity::Diastereotopic),
        });
        let mut entries = molecule_entries(&stereo_bond_canonicalization_molecule);
        entries.stereo_bonds[0].2.constraints = constraint.clone().into();
        entries.constraints =
            Constraint::StereoBond(StereoBondId(0), StereoKind::CisTrans, constraint).into();
        let source = Molecule::from_entries(entries);
        let reframed = reframe_stereo_bond(
            &source,
            StereoBondId(0),
            Permutation::from_image(&[2, 3, 0, 1]),
        );
        let left = canonicalize_structure(&source, &canonicalization_context)
            .expect("fixed molecule canonicalizes");
        let right = canonicalize_structure(&reframed, &canonicalization_context)
            .expect("reframed molecule canonicalizes");

        assert_eq!(
            selected_structure_key(&right),
            selected_structure_key(&left)
        );
    }

    #[rstest]
    fn test_canonicalize_full_constraints(canonicalization_context: CanonicalizationContext) {
        let plain = AtomForm::from_element(Element::C);
        let mut constrained = plain.clone();
        constrained.constraints = AtomConstraintForm::valence(4).into();
        let source = Molecule::from_entries(MoleculeEntries {
            atoms: vec![plain.clone(), constrained.clone()],
            constraints: Constraint::Molecule(MoleculeConstraint::ChargeSum {
                atoms: Some(vec![AtomId(1), AtomId(1)]),
                sum: NumForm::Lit(0),
            })
            .into(),
            ..Default::default()
        });
        let renumbered = source.remap(&molecule_correspondence(&[
            vec![1, 0],
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        ]));
        let expected = Molecule::from_entries(MoleculeEntries {
            atoms: vec![constrained, plain],
            constraints: Constraint::Molecule(MoleculeConstraint::ChargeSum {
                atoms: Some(vec![AtomId(0)]),
                sum: NumForm::Lit(0),
            })
            .into(),
            ..Default::default()
        });

        assert_eq!(
            canonicalize_full(&source, &canonicalization_context),
            Ok(expected.clone()),
        );
        assert_eq!(
            canonicalize_full(&renumbered, &canonicalization_context),
            Ok(expected.clone()),
        );
        assert_eq!(
            canonicalize_full(&expected, &canonicalization_context),
            Ok(expected),
        );
    }

    #[rstest]
    fn test_canonicalize_full_stereo_frame(
        repeated_ligand_canonicalization_molecule: Molecule,
        canonicalization_context: CanonicalizationContext,
    ) {
        let constraint = StereoAtomConstraintForm::Topicity(TopicityForm {
            pair: StereoLigandPair::new(0usize.into(), 2usize.into()),
            relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
        });
        let mut entries = molecule_entries(&repeated_ligand_canonicalization_molecule);
        entries.stereo_atoms[0].2.constraints = constraint.clone().into();
        entries.constraints =
            Constraint::StereoAtom(StereoAtomId(0), StereoKind::Tetrahedral, constraint).into();
        let source = Molecule::from_entries(entries);
        let reframed = reframe_stereo_atom(
            &source,
            StereoAtomId(0),
            Permutation::from_image(&[0, 1, 3, 2]),
        );
        let canonical = canonicalize_full(&source, &canonicalization_context)
            .expect("fixed molecule canonicalizes");

        assert_eq!(
            canonicalize_full(&reframed, &canonicalization_context),
            Ok(canonical.clone()),
        );
        assert_eq!(
            canonicalize_full(&canonical, &canonicalization_context),
            Ok(canonical),
        );
    }

    #[rstest]
    fn test_canonicalize_full_error(canonicalization_context: CanonicalizationContext) {
        let mut atom = AtomForm::from_element(Element::C);
        atom.constraints = AtomConstraintForm::Valence(NumForm::lit_set(Vec::<i64>::new())).into();
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![atom],
            ..Default::default()
        });

        assert_eq!(
            canonicalize_full(&molecule, &canonicalization_context),
            Err(MoleculeCanonicalizationError::Contradiction(Contradiction)),
        );
    }

    #[rstest]
    fn test_constraint_key() {
        let atom = Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4));
        let bond = Constraint::Bond(
            BondId(0),
            BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
        );
        let left = Constraint::And(vec![atom.clone(), bond.clone(), atom.clone()]);
        let right = Constraint::And(vec![bond, atom]);

        assert_eq!(constraint_key(&left), constraint_key(&right));
    }

    #[rstest]
    fn test_canonicalize_structure_selected_layer(
        initial_class_molecule: Molecule,
        canonicalization_context: CanonicalizationContext,
    ) {
        let canonical = canonicalize_structure(&initial_class_molecule, &canonicalization_context)
            .expect("fixed molecule canonicalizes");
        let canonical_again = canonicalize_structure(&canonical, &canonicalization_context)
            .expect("canonical molecule canonicalizes");

        assert_eq!(
            selected_structure_key(&canonical_again),
            selected_structure_key(&canonical),
        );
        assert_eq!(canonical_again.check_integrity(), Ok(()));
    }

    #[rstest]
    fn test_canonicalize_structure_renumbering(
        symmetric_stereo_canonicalization_molecule: Molecule,
        canonicalization_context: CanonicalizationContext,
    ) {
        let canonical = canonicalize_structure(
            &symmetric_stereo_canonicalization_molecule,
            &canonicalization_context,
        )
        .expect("fixed molecule canonicalizes");
        let expected = selected_structure_key(&canonical);

        for rank in 0..(1..=5).product() {
            let permutation = Permutation::unrank(5, rank);
            let correspondence = molecule_correspondence(&[
                (0..5).map(|index| permutation.apply(index)).collect(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                vec![0],
                Vec::new(),
            ]);
            let renumbered = symmetric_stereo_canonicalization_molecule.remap(&correspondence);
            let actual = canonicalize_structure(&renumbered, &canonicalization_context)
                .expect("renumbered molecule canonicalizes");

            assert_eq!(selected_structure_key(&actual), expected, "rank {rank}");
        }
    }

    #[rstest]
    #[case::nauty(AutomorphismAlgorithm::Nauty)]
    fn test_canonicalize_structure_minimum(
        symmetric_stereo_canonicalization_molecule: Molecule,
        canonicalization_context: CanonicalizationContext,
        #[case] algorithm: AutomorphismAlgorithm,
    ) {
        let incidence_graph =
            symmetric_stereo_canonicalization_molecule.incidence_graph(IncidenceLevel::Full);
        let (entity_keys, incidence_keys) = initial_class_keys(
            &symmetric_stereo_canonicalization_molecule,
            &incidence_graph,
        )
        .expect("fixed molecule has initial classes");
        let classes = rank_initial_classes(&entity_keys, &incidence_keys);
        let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
        let constitution_classes =
            constitution_entity_classes(&symmetric_stereo_canonicalization_molecule)
                .expect("fixed molecule has constitution classes");
        let descriptors = structure_partition_descriptors(
            &symmetric_stereo_canonicalization_molecule,
            &incidence_graph,
            &adapter,
            &entity_keys,
            &constitution_classes,
        )
        .expect("fixed molecule has structure descriptors");
        let leaf_candidate = |order: &[NodeId]| {
            structure_candidate(
                &symmetric_stereo_canonicalization_molecule,
                &incidence_graph,
                order,
            )
            .expect("structure descriptors establish normalization")
        };
        let no_prefix = |_: &OrderedPartition, _: &CanonicalCandidate<_>| false;
        let expected = exhaustive_minimum(&adapter, &leaf_candidate);

        for options in [
            CanonicalSearchOptions {
                automorphism_pruning: false,
                prefix_pruning: false,
                branch_order: BranchOrder::Node,
            },
            CanonicalSearchOptions {
                automorphism_pruning: false,
                prefix_pruning: false,
                branch_order: BranchOrder::ReverseNode,
            },
        ] {
            let actual = canonical_search(
                &adapter,
                &descriptors,
                algorithm,
                options,
                &leaf_candidate,
                &no_prefix,
            );
            assert_eq!(actual.candidate.key, expected.key, "{options:?}");
        }

        for options in [
            CanonicalSearchOptions {
                automorphism_pruning: false,
                prefix_pruning: false,
                branch_order: BranchOrder::ReverseNode,
            },
            CanonicalSearchOptions {
                automorphism_pruning: true,
                prefix_pruning: false,
                branch_order: BranchOrder::BackendCanonical,
            },
        ] {
            let (canonical, _) = canonicalize_structure_with_options(
                &symmetric_stereo_canonicalization_molecule,
                &CanonicalizationContext {
                    automorphism_algorithm: algorithm,
                    ..canonicalization_context
                },
                options,
            )
            .expect("fixed molecule canonicalizes");
            assert_eq!(
                selected_structure_key(&canonical),
                expected.key,
                "{options:?}"
            );
        }
    }

    #[rstest]
    fn test_canonicalize_structure_meso(
        meso_canonicalization_molecule: Molecule,
        canonicalization_context: CanonicalizationContext,
    ) {
        let correspondence = molecule_correspondence(&[
            vec![1, 0, 3, 2, 5, 4],
            vec![0, 3, 4, 1, 2],
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            vec![1, 0],
            Vec::new(),
        ]);
        let renumbered = meso_canonicalization_molecule.remap(&correspondence);

        assert_eq!(
            canonicalize_structure(&renumbered, &canonicalization_context),
            canonicalize_structure(&meso_canonicalization_molecule, &canonicalization_context,),
        );
    }

    #[rstest]
    #[case::kinded(StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32))]
    #[case::undetermined(StereoConfigurationForm::Undetermined)]
    fn test_canonicalize_structure_repeated_ligands(
        repeated_ligand_canonicalization_molecule: Molecule,
        canonicalization_context: CanonicalizationContext,
        #[case] configuration: StereoConfigurationForm,
    ) {
        let mut entries = molecule_entries(&repeated_ligand_canonicalization_molecule);
        entries.stereo_atoms[0].2.configuration = configuration;
        let source = Molecule::from_entries(entries);
        let reframed = reframe_stereo_atom(
            &source,
            StereoAtomId(0),
            Permutation::from_image(&[0, 1, 3, 2]),
        );

        assert_eq!(
            canonicalize_structure(&reframed, &canonicalization_context),
            canonicalize_structure(&source, &canonicalization_context),
        );
    }

    fn rank_paired_initial_classes(
        left: (&[InitialClassKey], &[InitialClassKey]),
        right: (&[InitialClassKey], &[InitialClassKey]),
    ) -> (InitialClasses, InitialClasses) {
        let ranks = left
            .0
            .iter()
            .chain(left.1)
            .chain(right.0)
            .chain(right.1)
            .cloned()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .enumerate()
            .map(|(rank, key)| (key, rank as u32))
            .collect::<BTreeMap<_, _>>();
        let rank =
            |entity_keys: &[InitialClassKey], incidence_keys: &[InitialClassKey]| InitialClasses {
                entities: entity_keys.iter().map(|key| ranks[key]).collect(),
                incidences: incidence_keys.iter().map(|key| ranks[key]).collect(),
            };

        (rank(left.0, left.1), rank(right.0, right.1))
    }

    fn colored_encoding_equivalent(
        left: &Molecule,
        right: &Molecule,
        level: IncidenceLevel,
    ) -> bool {
        fn key(adapter: &AutomorphismAdapter) -> Vec<u8> {
            adapter.graph().canonical_key(
                |node| {
                    let (domain, rank) = match adapter.class(node) {
                        AutomorphismClass::Entity(rank) => (0, rank),
                        AutomorphismClass::Incidence(rank) => (1, rank),
                    };
                    let mut color = vec![domain];
                    color.extend_from_slice(&rank.to_be_bytes());
                    color
                },
                |_| Vec::new(),
                AutomorphismAlgorithm::Nauty,
            )
        }

        let left_incidence = left.incidence_graph(level);
        let right_incidence = right.incidence_graph(level);
        let (left_entity_keys, left_incidence_keys) =
            initial_class_keys(left, &left_incidence).unwrap();
        let (right_entity_keys, right_incidence_keys) =
            initial_class_keys(right, &right_incidence).unwrap();
        let (left_classes, right_classes) = rank_paired_initial_classes(
            (&left_entity_keys, &left_incidence_keys),
            (&right_entity_keys, &right_incidence_keys),
        );
        let left_adapter = AutomorphismAdapter::new(&left_incidence, &left_classes);
        let right_adapter = AutomorphismAdapter::new(&right_incidence, &right_classes);

        key(&left_adapter) == key(&right_adapter)
    }

    fn permutations(count: usize) -> Vec<Vec<usize>> {
        fn visit(values: &mut [usize], position: usize, output: &mut Vec<Vec<usize>>) {
            if position == values.len() {
                output.push(values.to_vec());
                return;
            }
            for next in position..values.len() {
                values.swap(position, next);
                visit(values, position + 1, output);
                values.swap(position, next);
            }
        }

        let mut values = (0..count).collect::<Vec<_>>();
        let mut output = Vec::new();
        visit(&mut values, 0, &mut output);
        output
    }

    fn explicitly_dense_equivalent(left: &Molecule, right: &Molecule) -> bool {
        fn visit(
            family: usize,
            permutations: &[Vec<Vec<usize>>; 8],
            images: &mut [Vec<usize>; 8],
            left: &Molecule,
            right: &Molecule,
        ) -> bool {
            if family == images.len() {
                return left.equiv_under(right, &molecule_correspondence(images));
            }

            permutations[family].iter().any(|permutation| {
                images[family].clone_from(permutation);
                visit(family + 1, permutations, images, left, right)
            })
        }

        let left_counts = molecule_counts(left);
        if left_counts != molecule_counts(right) {
            return false;
        }
        let permutations = left_counts.map(permutations);
        visit(
            0,
            &permutations,
            &mut array::from_fn(|_| Vec::new()),
            left,
            right,
        )
    }

    fn reverse_correspondence(molecule: &Molecule) -> MoleculeCorrespondence {
        let images = molecule_counts(molecule).map(|count| (0..count).rev().collect::<Vec<_>>());
        molecule_correspondence(&images)
    }

    fn direct_graph_adapter(source: &Graph) -> AutomorphismAdapter {
        AutomorphismAdapter {
            graph: source.clone(),
            classes: vec![AutomorphismClass::Entity(0); source.node_count()],
            node_sources: source.node_ids().map(SubdivisionNodeSource::Node).collect(),
            incidence_nodes: vec![None; source.edge_count()],
            edge_sources: source.edge_ids().collect(),
            incidence_edges: source.edge_ids().map(|edge| vec![edge]).collect(),
            source_node_count: source.node_count(),
            entity_blocks: vec![source.node_ids().collect()],
        }
    }

    fn structural_leaf_key(
        order: &[NodeId],
        source: &Graph,
        classes: &InitialClasses,
    ) -> (Vec<u32>, Vec<(u32, u32, u32)>) {
        let mut positions = vec![0_u32; source.node_count()];
        for (position, node) in order.iter().enumerate() {
            positions[node.index()] = position as u32;
        }
        let entity_classes = order
            .iter()
            .map(|node| classes.entities[node.index()])
            .collect::<Vec<_>>();
        let mut incidences = source
            .edge_ids()
            .map(|edge| {
                let [first, second] = source.edge_endpoints(edge);
                let first = positions[first.index()];
                let second = positions[second.index()];
                (
                    classes.incidences[edge.index()],
                    first.min(second),
                    first.max(second),
                )
            })
            .collect::<Vec<_>>();
        incidences.sort_unstable();

        (entity_classes, incidences)
    }

    fn project_entries(mut entries: MoleculeEntries, level: IncidenceLevel) -> MoleculeEntries {
        entries.constraints = Constraints::new();
        match level {
            IncidenceLevel::Topology => {
                entries.dative.clear();
                entries.aromatic.clear();
                entries.multicenter.clear();
                entries.noncovalent.clear();
                entries.stereo_atoms.clear();
                entries.stereo_bonds.clear();
            }
            IncidenceLevel::Constitution => {
                entries.stereo_atoms.clear();
                entries.stereo_bonds.clear();
            }
            IncidenceLevel::Full => {}
        }
        entries
    }

    fn encoding_entries() -> MoleculeEntries {
        MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 4],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(2)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
            ],
            dative: vec![(
                vec![AtomId(0), AtomId(1)],
                AtomId(2),
                DativeBondForm::from_order(1),
            )],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemForm::from_electrons(vec![1, 2, 3]),
            )],
            multicenter: vec![(
                vec![AtomId(1), AtomId(2), AtomId(3)],
                MulticenterBondForm::from_electrons(vec![2, 0, 1]),
            )],
            noncovalent: vec![(
                AtomId(2),
                AtomId(3),
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            )],
            stereo_bonds: vec![(
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
            ..Default::default()
        }
    }

    #[rstest]
    fn test_incidence_cmp() {
        let incidences = [
            Incidence::BondEndpoint,
            Incidence::DativeDonor,
            Incidence::DativeAcceptor,
            Incidence::AromaticParticipant(NumForm::Undetermined),
            Incidence::AromaticParticipant(NumForm::Lit(-1)),
            Incidence::AromaticParticipant(NumForm::Lit(1)),
            Incidence::AromaticParticipant(NumForm::lit_set([0, 1])),
            Incidence::AromaticParticipant(NumForm::RangeFrom(0)),
            Incidence::AromaticParticipant(NumForm::RangeTo(0)),
            Incidence::AromaticParticipant(NumForm::var("x")),
            Incidence::AromaticParticipant(NumForm::pred_expr(PredExpr::Rel(
                ArithExpr::Var("x".into()),
                RelOp::Eq,
                ArithExpr::Lit(0),
            ))),
            Incidence::MulticenterParticipant(NumForm::Undetermined),
            Incidence::NoncovalentEndpoint,
            Incidence::StereoSite,
            Incidence::StereoLigand(StereoLigandKind::Atom),
            Incidence::StereoLigand(StereoLigandKind::ImplicitHydrogen),
            Incidence::StereoLigand(StereoLigandKind::LonePair),
        ];

        for pair in incidences.windows(2) {
            assert_eq!(pair[0].cmp(&pair[1]), Ordering::Less);
        }
        for lhs in &incidences {
            for rhs in &incidences {
                assert_eq!(
                    lhs.cmp(rhs),
                    incidence_key(lhs)
                        .unwrap()
                        .cmp(&incidence_key(rhs).unwrap()),
                );
            }
        }
    }

    #[rstest]
    #[case::atom(
        Entity::Atom(AtomId(0)),
        vec![
            FieldPosition(0),
            FieldPosition(1),
            FieldPosition(2),
            FieldPosition(3),
            FieldPosition(4),
            FieldPosition(5),
        ]
    )]
    #[case::bond(
        Entity::Bond(BondId(0)),
        vec![FieldPosition(1), FieldPosition(2), FieldPosition(3)]
    )]
    #[case::dative_bond(
        Entity::DativeBond(DativeBondId(0)),
        vec![FieldPosition(2)]
    )]
    #[case::aromatic_system(
        Entity::AromaticSystem(AromaticSystemId(0)),
        vec![FieldPosition(2), FieldPosition(3)]
    )]
    #[case::multicenter_bond(
        Entity::MulticenterBond(MulticenterBondId(0)),
        vec![FieldPosition(2), FieldPosition(3)]
    )]
    #[case::noncovalent_bond(
        Entity::NoncovalentBond(NoncovalentBondId(0)),
        vec![FieldPosition(1)]
    )]
    #[case::stereo_atom(
        Entity::StereoAtom(StereoAtomId(0)),
        vec![FieldPosition(2)]
    )]
    #[case::stereo_bond(
        Entity::StereoBond(StereoBondId(0)),
        vec![FieldPosition(2)]
    )]
    fn test_entity_class_key_field_positions(
        initial_class_molecule: Molecule,
        #[case] entity: Entity,
        #[case] expected: Vec<FieldPosition>,
    ) {
        let InitialClassKey::Entity {
            value: CanonicalKeyValue::Product(fields),
            ..
        } = entity_class_key(&initial_class_molecule, entity).unwrap()
        else {
            panic!("entity class key must be a product");
        };

        assert_eq!(
            fields
                .into_iter()
                .map(|field| field.position)
                .collect::<Vec<_>>(),
            expected,
        );
    }

    #[rstest]
    #[case::normalized_atom(Entity::Atom(AtomId(0)), Entity::Atom(AtomId(1)), true)]
    #[case::atom_element(Entity::Atom(AtomId(0)), Entity::Atom(AtomId(2)), false)]
    #[case::atom_charge(Entity::Atom(AtomId(0)), Entity::Atom(AtomId(3)), false)]
    #[case::normalized_bond(Entity::Bond(BondId(0)), Entity::Bond(BondId(1)), true)]
    #[case::bond_order(Entity::Bond(BondId(0)), Entity::Bond(BondId(2)), false)]
    #[case::normalized_dative(
        Entity::DativeBond(DativeBondId(0)),
        Entity::DativeBond(DativeBondId(1)),
        true
    )]
    #[case::dative_order(
        Entity::DativeBond(DativeBondId(0)),
        Entity::DativeBond(DativeBondId(2)),
        false
    )]
    #[case::aromatic_electrons_excluded(
        Entity::AromaticSystem(AromaticSystemId(0)),
        Entity::AromaticSystem(AromaticSystemId(1)),
        true
    )]
    #[case::aromatic_charge(
        Entity::AromaticSystem(AromaticSystemId(0)),
        Entity::AromaticSystem(AromaticSystemId(2)),
        false
    )]
    #[case::multicenter_electrons_excluded(
        Entity::MulticenterBond(MulticenterBondId(0)),
        Entity::MulticenterBond(MulticenterBondId(1)),
        true
    )]
    #[case::noncovalent_kind_equal(
        Entity::NoncovalentBond(NoncovalentBondId(0)),
        Entity::NoncovalentBond(NoncovalentBondId(1)),
        true
    )]
    #[case::noncovalent_kind_distinct(
        Entity::NoncovalentBond(NoncovalentBondId(0)),
        Entity::NoncovalentBond(NoncovalentBondId(2)),
        false
    )]
    #[case::stereo_atom_configuration_excluded(
        Entity::StereoAtom(StereoAtomId(0)),
        Entity::StereoAtom(StereoAtomId(1)),
        true
    )]
    #[case::stereo_atom_kind(
        Entity::StereoAtom(StereoAtomId(0)),
        Entity::StereoAtom(StereoAtomId(2)),
        false
    )]
    #[case::stereo_bond_configuration_excluded(
        Entity::StereoBond(StereoBondId(0)),
        Entity::StereoBond(StereoBondId(1)),
        true
    )]
    #[case::entity_kind(Entity::Atom(AtomId(0)), Entity::Bond(BondId(0)), false)]
    fn test_initial_classes(
        initial_class_molecule: Molecule,
        #[case] lhs: Entity,
        #[case] rhs: Entity,
        #[case] expected_equal: bool,
    ) {
        let incidence_graph = initial_class_molecule.incidence_graph(IncidenceLevel::Full);
        let classes = initial_classes(&initial_class_molecule, &incidence_graph).unwrap();
        let lhs_class = classes.entities[incidence_graph.node_of(lhs).index()];
        let rhs_class = classes.entities[incidence_graph.node_of(rhs).index()];

        assert_eq!(lhs_class == rhs_class, expected_equal);
    }

    #[rstest]
    fn test_initial_classes_incidence(initial_class_molecule: Molecule) {
        let incidence_graph = initial_class_molecule.incidence_graph(IncidenceLevel::Full);
        let classes = initial_classes(&initial_class_molecule, &incidence_graph).unwrap();
        let incidences = incidence_graph
            .incidences()
            .map(|(edge, incidence)| (incidence, classes.incidences[edge.index()]))
            .collect::<Vec<_>>();

        for (lhs, lhs_class) in &incidences {
            for (rhs, rhs_class) in &incidences {
                assert_eq!(lhs_class == rhs_class, lhs == rhs);
            }
        }
        for entity_class in &classes.entities {
            for (_, incidence_class) in &incidences {
                assert_ne!(entity_class, incidence_class);
            }
        }
    }

    #[rstest]
    fn test_topology_comparison_key() {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C)
                    .with_charge(NumForm::ArithExpr(Box::new(ArithExpr::Sum(vec![
                        ArithExpr::Lit(1),
                        ArithExpr::Lit(2),
                    ]))))
                    .with_constraint(AtomConstraintForm::Valence(NumForm::lit_set([]))),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![(
                AtomId(0),
                AtomId(1),
                BondForm::new(NumForm::ArithExpr(Box::new(ArithExpr::Sum(vec![
                    ArithExpr::Lit(0),
                    ArithExpr::Lit(1),
                ]))))
                .with_charge(-1_i64),
            )],
            dative: vec![(
                vec![AtomId(0)],
                AtomId(1),
                DativeBondForm::new(NumForm::lit_set([])),
            )],
            ..Default::default()
        });
        let incidence_graph = molecule.incidence_graph(IncidenceLevel::Topology);
        let order = [
            incidence_graph.node_of(Entity::Atom(AtomId(0))),
            incidence_graph.node_of(Entity::Atom(AtomId(1))),
            incidence_graph.node_of(Entity::Bond(BondId(0))),
        ];
        let undetermined_num = NumForm::Undetermined;
        let undetermined_isotope = IsotopeMassForm::Undetermined;
        let undetermined_spin = UnpairedElectronsForm::default();

        assert_eq!(
            topology_comparison_key(&molecule, &incidence_graph, &order),
            Ok(CanonicalComparisonKey {
                entity_blocks: vec![
                    PositionedKey {
                        position: EntityBlockPosition::ATOM,
                        value: sequence([
                            product([
                                element_form_key(&ElementForm::Lit(Element::C)),
                                isotope_mass_form_key(&undetermined_isotope),
                                num_form_key(&NumForm::Lit(3)),
                                num_form_key(&undetermined_num),
                                num_form_key(&undetermined_num),
                                unpaired_electrons_form_key(&undetermined_spin),
                            ]),
                            product([
                                element_form_key(&ElementForm::Lit(Element::O)),
                                isotope_mass_form_key(&undetermined_isotope),
                                num_form_key(&undetermined_num),
                                num_form_key(&undetermined_num),
                                num_form_key(&undetermined_num),
                                unpaired_electrons_form_key(&undetermined_spin),
                            ]),
                        ]),
                    },
                    PositionedKey {
                        position: EntityBlockPosition::BOND,
                        value: sequence([positioned_product([
                            (
                                0,
                                product([
                                    CanonicalKeyValue::Unsigned(0),
                                    CanonicalKeyValue::Unsigned(1),
                                ]),
                            ),
                            (1, num_form_key(&NumForm::Lit(1))),
                            (2, num_form_key(&NumForm::Lit(-1))),
                            (3, unpaired_electrons_form_key(&undetermined_spin)),
                        ])]),
                    },
                ],
                constraints: Vec::new(),
            }),
        );
    }

    #[rstest]
    fn test_topology_comparison_key_dense_remapping() {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::O),
                AtomForm::from_element(Element::N),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(2)),
                (
                    AtomId(1),
                    AtomId(2),
                    BondForm::from_order(1).with_charge(-1_i64),
                ),
            ],
            ..Default::default()
        });
        let remapped = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![
                (
                    AtomId(0),
                    AtomId(2),
                    BondForm::from_order(1).with_charge(-1_i64),
                ),
                (AtomId(1), AtomId(2), BondForm::from_order(2)),
            ],
            ..Default::default()
        });
        let incidence_graph = molecule.incidence_graph(IncidenceLevel::Topology);
        let remapped_incidence_graph = remapped.incidence_graph(IncidenceLevel::Topology);
        let order = [
            incidence_graph.node_of(Entity::Atom(AtomId(2))),
            incidence_graph.node_of(Entity::Atom(AtomId(0))),
            incidence_graph.node_of(Entity::Atom(AtomId(1))),
            incidence_graph.node_of(Entity::Bond(BondId(1))),
            incidence_graph.node_of(Entity::Bond(BondId(0))),
        ];
        let remapped_order = remapped_incidence_graph
            .graph()
            .node_ids()
            .collect::<Vec<_>>();

        assert_eq!(
            topology_comparison_key(&molecule, &incidence_graph, &order),
            topology_comparison_key(&remapped, &remapped_incidence_graph, &remapped_order),
        );
    }

    #[rstest]
    #[case::dative_bond(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            dative: vec![(
                vec![AtomId(0), AtomId(1)],
                AtomId(2),
                DativeBondForm::new(NumForm::RangeFrom(-1)),
            )],
            ..Default::default()
        }),
        EntityBlockPosition::DATIVE_BOND,
        positioned_product([
            (
                0,
                sequence([
                    CanonicalKeyValue::Unsigned(1),
                    CanonicalKeyValue::Unsigned(2),
                ]),
            ),
            (1, CanonicalKeyValue::Unsigned(0)),
            (2, num_form_key(&NumForm::RangeFrom(-1))),
        ]),
    )]
    #[case::aromatic_system(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            aromatic: vec![(
                vec![AtomId(0), AtomId(2)],
                AromaticSystemForm::from_electrons(vec![1, 2])
                    .with_charge(NumForm::var("q")),
            )],
            ..Default::default()
        }),
        EntityBlockPosition::AROMATIC_SYSTEM,
        positioned_product([
            (
                0,
                sequence([
                    CanonicalKeyValue::Unsigned(0),
                    CanonicalKeyValue::Unsigned(2),
                ]),
            ),
            (
                1,
                electron_counts_form_key(&ElectronCountsForm::Lit(vec![2, 1])),
            ),
            (2, num_form_key(&NumForm::var("q"))),
            (
                3,
                unpaired_electrons_form_key(&UnpairedElectronsForm::default()),
            ),
        ]),
    )]
    #[case::multicenter_bond(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            multicenter: vec![(
                vec![AtomId(0), AtomId(2)],
                MulticenterBondForm::from_electrons(vec![2, 0])
                    .with_charge(NumForm::RangeTo(1)),
            )],
            ..Default::default()
        }),
        EntityBlockPosition::MULTICENTER_BOND,
        positioned_product([
            (
                0,
                sequence([
                    CanonicalKeyValue::Unsigned(0),
                    CanonicalKeyValue::Unsigned(2),
                ]),
            ),
            (
                1,
                electron_counts_form_key(&ElectronCountsForm::Lit(vec![0, 2])),
            ),
            (2, num_form_key(&NumForm::RangeTo(1))),
            (
                3,
                unpaired_electrons_form_key(&UnpairedElectronsForm::default()),
            ),
        ]),
    )]
    #[case::noncovalent_bond(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            noncovalent: vec![(
                AtomId(0),
                AtomId(2),
                NoncovalentBondForm::default(),
            )],
            ..Default::default()
        }),
        EntityBlockPosition::NONCOVALENT_BOND,
        positioned_product([
            (
                0,
                product([
                    CanonicalKeyValue::Unsigned(0),
                    CanonicalKeyValue::Unsigned(2),
                ]),
            ),
            (
                1,
                noncovalent_bond_kind_form_key(&NoncovalentBondKindForm::Undetermined),
            ),
        ]),
    )]
    fn test_constitution_comparison_key(
        #[case] molecule: Molecule,
        #[case] position: EntityBlockPosition,
        #[case] expected: CanonicalKeyValue,
    ) {
        let incidence_graph = molecule.incidence_graph(IncidenceLevel::Constitution);
        let mut atom_ids = molecule.atoms().ids().collect::<Vec<_>>();
        atom_ids.reverse();
        let mut order = atom_ids
            .into_iter()
            .map(|id| incidence_graph.node_of(Entity::Atom(id)))
            .collect::<Vec<_>>();
        order.extend(
            incidence_graph
                .graph()
                .node_ids()
                .filter(|&node| !matches!(incidence_graph.entity(node), Entity::Atom(_))),
        );
        let key = constitution_comparison_key(&molecule, &incidence_graph, &order).unwrap();

        assert_eq!(
            key.entity_blocks
                .into_iter()
                .find(|block| block.position == position),
            Some(PositionedKey {
                position,
                value: sequence([expected]),
            }),
        );
    }

    #[rstest]
    fn test_constitution_comparison_key_dense_remapping() {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
                AtomForm::from_element(Element::F),
            ],
            dative: vec![(
                vec![AtomId(0), AtomId(2)],
                AtomId(1),
                DativeBondForm::new(NumForm::RangeFrom(1)),
            )],
            aromatic: vec![(
                vec![AtomId(0), AtomId(2)],
                AromaticSystemForm::from_electrons(vec![2, 1]).with_charge(NumForm::var("q")),
            )],
            multicenter: vec![(
                vec![AtomId(1), AtomId(3)],
                MulticenterBondForm::new(ElectronCountsForm::Undetermined)
                    .with_charge(NumForm::RangeTo(2)),
            )],
            noncovalent: vec![(AtomId(0), AtomId(3), NoncovalentBondForm::default())],
            ..Default::default()
        });
        let correspondence = molecule_correspondence(&[
            vec![3, 1, 0, 2],
            Vec::new(),
            vec![0],
            vec![0],
            vec![0],
            vec![0],
            Vec::new(),
            Vec::new(),
        ]);
        let remapped = molecule.remap(&correspondence);
        let incidence_graph = molecule.incidence_graph(IncidenceLevel::Constitution);
        let remapped_incidence_graph = remapped.incidence_graph(IncidenceLevel::Constitution);
        let mut order = incidence_graph.graph().node_ids().collect::<Vec<_>>();
        order.reverse();
        let remapped_order = order
            .iter()
            .map(|&node| incidence_graph.entity(node))
            .map(|entity| {
                correspondence
                    .right_of(entity)
                    .expect("dense correspondence maps every entity")
            })
            .map(|entity| remapped_incidence_graph.node_of(entity))
            .collect::<Vec<_>>();

        assert_eq!(
            constitution_comparison_key(&molecule, &incidence_graph, &order),
            constitution_comparison_key(&remapped, &remapped_incidence_graph, &remapped_order),
        );
    }

    #[rstest]
    fn test_constitution_comparison_key_excluded_data() {
        let molecule = Molecule::from_entries(project_entries(
            encoding_entries(),
            IncidenceLevel::Constitution,
        ));
        let mut excluded = Molecule::from_entries(encoding_entries());
        excluded.modify_atoms(|atom| {
            atom.with_constraint(AtomConstraintForm::Valence(NumForm::Lit(4)))
        });
        let incidence_graph = molecule.incidence_graph(IncidenceLevel::Constitution);
        let excluded_incidence_graph = excluded.incidence_graph(IncidenceLevel::Constitution);
        let order = incidence_graph.graph().node_ids().collect::<Vec<_>>();
        let excluded_order = excluded_incidence_graph
            .graph()
            .node_ids()
            .collect::<Vec<_>>();

        assert_eq!(
            constitution_comparison_key(&molecule, &incidence_graph, &order),
            constitution_comparison_key(&excluded, &excluded_incidence_graph, &excluded_order),
        );
    }

    #[rstest]
    #[case::dative(Entity::DativeBond(DativeBondId(0)))]
    #[case::aromatic(Entity::AromaticSystem(AromaticSystemId(0)))]
    #[case::multicenter(Entity::MulticenterBond(MulticenterBondId(0)))]
    #[case::noncovalent(Entity::NoncovalentBond(NoncovalentBondId(0)))]
    #[case::stereo_atom(Entity::StereoAtom(StereoAtomId(0)))]
    #[case::stereo_bond(Entity::StereoBond(StereoBondId(0)))]
    fn test_correspondence_from_order(initial_class_molecule: Molecule, #[case] excluded: Entity) {
        let incidence_graph = initial_class_molecule.incidence_graph(IncidenceLevel::Topology);
        let order = incidence_graph.graph().node_ids().collect::<Vec<_>>();
        let correspondence =
            correspondence_from_order(&initial_class_molecule, &incidence_graph, &order);

        assert_eq!(correspondence.right_of(excluded), Some(excluded));
    }

    #[rstest]
    #[case::localized_bond(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 2],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        }),
        vec![Incidence::BondEndpoint, Incidence::BondEndpoint],
        0,
    )]
    #[case::repeated_virtual_ligand_anchor(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            )],
            ..Default::default()
        }),
        vec![
            Incidence::StereoSite,
            Incidence::StereoLigand(StereoLigandKind::Atom),
            Incidence::StereoLigand(StereoLigandKind::Atom),
            Incidence::StereoLigand(StereoLigandKind::ImplicitHydrogen),
            Incidence::StereoLigand(StereoLigandKind::ImplicitHydrogen),
        ],
        5,
    )]
    fn test_automorphism_adapter_new(
        #[case] molecule: Molecule,
        #[case] expected_incidences: Vec<Incidence>,
        #[case] expected_occurrence_nodes: usize,
    ) {
        let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
        let classes = initial_classes(&molecule, &incidence_graph).unwrap();
        let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
        let source = incidence_graph.graph();

        assert_eq!(
            incidence_graph
                .incidences()
                .map(|(_, incidence)| incidence.clone())
                .collect::<Vec<_>>(),
            expected_incidences,
        );
        assert_eq!(
            adapter.graph().node_count(),
            source.node_count() + expected_occurrence_nodes,
        );
        assert_eq!(
            adapter.graph().edge_count(),
            source
                .edge_ids()
                .map(|edge| adapter.incidence_edges_of(edge).len())
                .sum(),
        );
        assert!(adapter.graph().is_simple());

        for node in source.node_ids() {
            let adapter_node = adapter
                .node_of(SubdivisionNodeSource::Node(node))
                .expect("every source entity remains an adapter node");
            assert_eq!(
                adapter.node_source(adapter_node),
                SubdivisionNodeSource::Node(node),
            );
            assert_eq!(
                adapter.class(adapter_node),
                AutomorphismClass::Entity(classes.entities[node.index()]),
            );
        }
        for edge in source.edge_ids() {
            if let Some(adapter_node) = adapter.node_of(SubdivisionNodeSource::Edge(edge)) {
                assert_eq!(
                    adapter.node_source(adapter_node),
                    SubdivisionNodeSource::Edge(edge),
                );
                assert_eq!(
                    adapter.class(adapter_node),
                    AutomorphismClass::Incidence(classes.incidences[edge.index()]),
                );
            }
            for &incidence_edge in adapter.incidence_edges_of(edge) {
                assert_eq!(adapter.edge_source(incidence_edge), edge);
            }
        }
    }

    #[rstest]
    fn test_automorphism_adapter_automorphisms() {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        });
        let incidence_graph = molecule.incidence_graph(IncidenceLevel::Topology);
        let classes = initial_classes(&molecule, &incidence_graph).unwrap();
        let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);

        assert_eq!(
            adapter.automorphisms(AutomorphismAlgorithm::Nauty),
            ProjectedAutomorphismOutput {
                orbits: vec![NodeId(0), NodeId(0), NodeId(2)],
                canonical_labels: vec![NodeId(0), NodeId(1), NodeId(2)],
                generators: vec![vec![NodeId(1), NodeId(0), NodeId(2)]],
            },
        );
    }

    #[rstest]
    #[case::ordered_classes(
        vec![2_u32, 0, 2, 1, 0],
        OrderedPartition {
            cells: vec![
                vec![NodeId(1), NodeId(4)],
                vec![NodeId(3)],
                vec![NodeId(0), NodeId(2)],
            ],
        },
    )]
    fn test_ordered_partition_from_descriptors(
        #[case] descriptors: Vec<u32>,
        #[case] expected: OrderedPartition,
    ) {
        assert_eq!(OrderedPartition::from_descriptors(&descriptors), expected,);
    }

    #[rstest]
    #[case::path(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3]]),
        OrderedPartition {
            cells: vec![vec![NodeId(1), NodeId(2)], vec![NodeId(0), NodeId(3)]],
        },
    )]
    fn test_ordered_partition_refine(#[case] graph: Graph, #[case] expected: OrderedPartition) {
        assert_eq!(
            OrderedPartition::from_descriptors(&[0_u32; 4]).refine(&graph),
            expected,
        );
    }

    #[rstest]
    #[case::first_cell(
        OrderedPartition {
            cells: vec![vec![NodeId(0), NodeId(3)], vec![NodeId(1), NodeId(2)]],
        },
        0,
        NodeId(3),
        OrderedPartition {
            cells: vec![
                vec![NodeId(3)],
                vec![NodeId(0)],
                vec![NodeId(1), NodeId(2)],
            ],
        },
    )]
    fn test_ordered_partition_individualize(
        #[case] partition: OrderedPartition,
        #[case] cell_index: usize,
        #[case] node: NodeId,
        #[case] expected: OrderedPartition,
    ) {
        assert_eq!(partition.individualize(cell_index, node), expected);
    }

    #[rstest]
    #[case::path(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 4],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
            ],
            ..Default::default()
        }),
    )]
    #[case::branched(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 4],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
            ],
            ..Default::default()
        }),
    )]
    #[case::distinct_attributes(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::O),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(2)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
            ],
            ..Default::default()
        }),
    )]
    fn test_canonical_search(#[case] molecule: Molecule) {
        let incidence_graph = molecule.incidence_graph(IncidenceLevel::Topology);
        let (entity_keys, incidence_keys) =
            initial_class_keys(&molecule, &incidence_graph).unwrap();
        let classes = rank_initial_classes(&entity_keys, &incidence_keys);
        let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
        let descriptors = partition_descriptors(&adapter, &entity_keys, &incidence_keys);
        let leaf_candidate = |order: &[NodeId]| {
            topology_candidate(&molecule, &incidence_graph, order)
                .expect("selected topology values normalize")
        };
        let no_prefix = |_: &OrderedPartition, _: &CanonicalCandidate<_>| false;
        let expected = exhaustive_minimum(&adapter, &leaf_candidate);
        let unpruned = canonical_search(
            &adapter,
            &descriptors,
            AutomorphismAlgorithm::Nauty,
            CanonicalSearchOptions {
                automorphism_pruning: false,
                prefix_pruning: false,
                branch_order: BranchOrder::Node,
            },
            &leaf_candidate,
            &no_prefix,
        );
        let reversed = canonical_search(
            &adapter,
            &descriptors,
            AutomorphismAlgorithm::Nauty,
            CanonicalSearchOptions {
                automorphism_pruning: false,
                prefix_pruning: false,
                branch_order: BranchOrder::ReverseNode,
            },
            &leaf_candidate,
            &no_prefix,
        );
        let pruned = canonical_search(
            &adapter,
            &descriptors,
            AutomorphismAlgorithm::Nauty,
            CanonicalSearchOptions {
                automorphism_pruning: true,
                prefix_pruning: false,
                branch_order: BranchOrder::BackendCanonical,
            },
            &leaf_candidate,
            &no_prefix,
        );

        assert_eq!(unpruned.candidate.key, expected.key);
        assert_eq!(reversed.candidate.key, expected.key);
        assert_eq!(pruned.candidate.key, expected.key);
        assert!(pruned.stats.visited_leaves <= unpruned.stats.visited_leaves);
    }

    #[rstest]
    fn test_canonical_search_color_classes() {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::O),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(2)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
            ],
            ..Default::default()
        });
        let incidence_graph = molecule.incidence_graph(IncidenceLevel::Topology);
        let (entity_keys, incidence_keys) =
            initial_class_keys(&molecule, &incidence_graph).unwrap();
        let classes = rank_initial_classes(&entity_keys, &incidence_keys);
        let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
        let mut relabeled = adapter.clone();
        relabeled.classes.iter_mut().for_each(|class| {
            *class = match *class {
                AutomorphismClass::Entity(value) => AutomorphismClass::Entity(u32::MAX - value),
                AutomorphismClass::Incidence(value) => {
                    AutomorphismClass::Incidence(u32::MAX - value)
                }
            }
        });
        let descriptors = partition_descriptors(&adapter, &entity_keys, &incidence_keys);
        let leaf_candidate = |order: &[NodeId]| {
            topology_candidate(&molecule, &incidence_graph, order)
                .expect("selected topology values normalize")
        };
        let no_prefix = |_: &OrderedPartition, _: &CanonicalCandidate<_>| false;
        let options = CanonicalSearchOptions {
            automorphism_pruning: true,
            prefix_pruning: false,
            branch_order: BranchOrder::BackendCanonical,
        };

        let expected = canonical_search(
            &adapter,
            &descriptors,
            AutomorphismAlgorithm::Nauty,
            options,
            &leaf_candidate,
            &no_prefix,
        );
        let actual = canonical_search(
            &relabeled,
            &descriptors,
            AutomorphismAlgorithm::Nauty,
            options,
            &leaf_candidate,
            &no_prefix,
        );

        assert_eq!(actual.candidate.key, expected.candidate.key);
    }

    #[rstest]
    fn test_canonicalize_topology(canonicalization_context: CanonicalizationContext) {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::O),
                AtomForm::from_element(Element::C).with_charge(NumForm::ArithExpr(Box::new(
                    ArithExpr::Sum(vec![ArithExpr::Lit(1), ArithExpr::Lit(2)]),
                ))),
                AtomForm::from_element(Element::N),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(2)),
                (
                    AtomId(1),
                    AtomId(2),
                    BondForm::new(NumForm::ArithExpr(Box::new(ArithExpr::Sum(vec![
                        ArithExpr::Lit(0),
                        ArithExpr::Lit(1),
                    ])))),
                ),
            ],
            dative: vec![(
                vec![AtomId(0)],
                AtomId(2),
                DativeBondForm::new(NumForm::ArithExpr(Box::new(ArithExpr::Sum(vec![
                    ArithExpr::Lit(0),
                    ArithExpr::Lit(1),
                ])))),
            )],
            ..Default::default()
        });
        let expected = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C).with_charge(3_i64),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(2)),
            ],
            dative: vec![(vec![AtomId(2)], AtomId(1), DativeBondForm::from_order(1))],
            ..Default::default()
        });

        assert_eq!(
            canonicalize_topology(&molecule, &canonicalization_context),
            Ok(expected),
        );
    }

    #[rstest]
    #[case::disconnected(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::O),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::H),
            ],
            bonds: vec![
                (AtomId(0), AtomId(2), BondForm::from_order(2)),
                (AtomId(1), AtomId(3), BondForm::from_order(1)),
            ],
            ..Default::default()
        }),
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::H),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(1), AtomId(3), BondForm::from_order(2)),
            ],
            ..Default::default()
        }),
    )]
    fn test_canonicalize_topology_components(
        canonicalization_context: CanonicalizationContext,
        #[case] molecule: Molecule,
        #[case] expected: Molecule,
    ) {
        assert_eq!(
            canonicalize_topology(&molecule, &canonicalization_context),
            Ok(expected),
        );
    }

    #[rstest]
    #[case::selected_atom(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)
                .with_charge(NumForm::LitSet(Box::default()))],
            ..Default::default()
        }),
    )]
    #[case::excluded_dative(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::B),
            ],
            dative: vec![(
                vec![AtomId(0)],
                AtomId(1),
                DativeBondForm::new(NumForm::LitSet(Box::default())),
            )],
            ..Default::default()
        }),
    )]
    fn test_canonicalize_topology_error(
        canonicalization_context: CanonicalizationContext,
        #[case] molecule: Molecule,
    ) {
        assert_eq!(
            canonicalize_topology(&molecule, &canonicalization_context),
            Err(MoleculeCanonicalizationError::Contradiction(Contradiction,)),
        );
    }

    #[rstest]
    fn test_canonicalize_topology_excluded_data(canonicalization_context: CanonicalizationContext) {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 2],
            dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1))],
            ..Default::default()
        });
        let remapping = molecule_correspondence(&[
            vec![1, 0],
            Vec::new(),
            vec![0],
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        ]);
        let remapped = molecule.remap(&remapping);

        let (canonical, correspondence) = canonicalize_topology_with_options(
            &molecule,
            &canonicalization_context,
            CanonicalSearchOptions {
                automorphism_pruning: true,
                prefix_pruning: false,
                branch_order: BranchOrder::BackendCanonical,
            },
        )
        .unwrap();
        let (canonical_remapped, remapped_correspondence) = canonicalize_topology_with_options(
            &remapped,
            &canonicalization_context,
            CanonicalSearchOptions {
                automorphism_pruning: true,
                prefix_pruning: false,
                branch_order: BranchOrder::BackendCanonical,
            },
        )
        .unwrap();
        let canonical_incidence = canonical.incidence_graph(IncidenceLevel::Topology);
        let canonical_remapped_incidence =
            canonical_remapped.incidence_graph(IncidenceLevel::Topology);
        let canonical_again = canonicalize_topology(&canonical, &canonicalization_context).unwrap();
        let canonical_again_incidence = canonical_again.incidence_graph(IncidenceLevel::Topology);

        assert!(molecule.equiv_under(&canonical, &correspondence));
        assert!(remapped.equiv_under(&canonical_remapped, &remapped_correspondence));
        assert_eq!(canonical.check_integrity(), Ok(()));
        assert_eq!(canonical_remapped.check_integrity(), Ok(()));
        assert_eq!(
            topology_comparison_key(
                &canonical,
                &canonical_incidence,
                &canonical_incidence.graph().node_ids().collect::<Vec<_>>(),
            ),
            topology_comparison_key(
                &canonical_again,
                &canonical_again_incidence,
                &canonical_again_incidence
                    .graph()
                    .node_ids()
                    .collect::<Vec<_>>(),
            ),
        );
        assert_eq!(
            topology_comparison_key(
                &canonical,
                &canonical_incidence,
                &canonical_incidence.graph().node_ids().collect::<Vec<_>>(),
            ),
            topology_comparison_key(
                &canonical_remapped,
                &canonical_remapped_incidence,
                &canonical_remapped_incidence
                    .graph()
                    .node_ids()
                    .collect::<Vec<_>>(),
            ),
        );
    }

    #[rstest]
    #[case::order_four(4)]
    fn test_canonicalize_topology_exhaustive_domain(
        canonicalization_context: CanonicalizationContext,
        #[case] atom_count: usize,
    ) {
        let endpoint_pairs = (0..atom_count as u32)
            .flat_map(|first| ((first + 1)..atom_count as u32).map(move |second| [first, second]))
            .collect::<Vec<_>>();

        for edge_mask in 0..(1_u64 << endpoint_pairs.len()) {
            let bonds = endpoint_pairs
                .iter()
                .enumerate()
                .filter_map(|(position, &[first, second])| {
                    ((edge_mask >> position) & 1 == 1).then_some((
                        AtomId(first),
                        AtomId(second),
                        BondForm::from_order(1),
                    ))
                })
                .collect::<Vec<_>>();
            let molecule = Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); atom_count],
                bonds,
                ..Default::default()
            });
            let incidence_graph = molecule.incidence_graph(IncidenceLevel::Topology);
            let (entity_keys, incidence_keys) =
                initial_class_keys(&molecule, &incidence_graph).unwrap();
            let classes = rank_initial_classes(&entity_keys, &incidence_keys);
            let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
            let leaf_candidate =
                |order: &[NodeId]| topology_candidate(&molecule, &incidence_graph, order).unwrap();
            let expected = exhaustive_minimum(&adapter, &leaf_candidate);

            let (canonical, correspondence) = canonicalize_topology_with_options(
                &molecule,
                &canonicalization_context,
                CanonicalSearchOptions {
                    automorphism_pruning: true,
                    prefix_pruning: false,
                    branch_order: BranchOrder::BackendCanonical,
                },
            )
            .unwrap();
            let (unpruned, _) = canonicalize_topology_with_options(
                &molecule,
                &canonicalization_context,
                CanonicalSearchOptions {
                    automorphism_pruning: false,
                    prefix_pruning: false,
                    branch_order: BranchOrder::ReverseNode,
                },
            )
            .unwrap();
            let canonical_incidence = canonical.incidence_graph(IncidenceLevel::Topology);
            let canonical_order = canonical_incidence.graph().node_ids().collect::<Vec<_>>();

            assert_eq!(
                topology_comparison_key(&canonical, &canonical_incidence, &canonical_order),
                Ok(expected.key),
                "edge mask {edge_mask:#08b}",
            );
            assert_eq!(unpruned, canonical, "edge mask {edge_mask:#08b}");
            assert!(
                molecule.equiv_under(&canonical, &correspondence),
                "edge mask {edge_mask:#08b}",
            );
            assert_eq!(canonical.check_integrity(), Ok(()));
            assert_eq!(
                canonicalize_topology(&canonical, &canonicalization_context),
                Ok(canonical.clone()),
                "edge mask {edge_mask:#08b}",
            );

            for (index, atom_images) in permutations(atom_count).into_iter().enumerate() {
                let bond_count = molecule.bonds().count();
                let bond_images = if index % 2 == 0 {
                    (0..bond_count).collect::<Vec<_>>()
                } else {
                    (0..bond_count).rev().collect::<Vec<_>>()
                };
                let renumbering = molecule_correspondence(&[
                    atom_images,
                    bond_images,
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                ]);
                let renumbered = molecule.remap(&renumbering);

                assert_eq!(
                    canonicalize_topology(&renumbered, &canonicalization_context),
                    Ok(canonical.clone()),
                    "edge mask {edge_mask:#08b}, renumbering {index}",
                );
            }
        }
    }

    #[rstest]
    fn test_canonicalize_constitution(canonicalization_context: CanonicalizationContext) {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::O).with_constraint(AtomConstraintForm::Valence(
                    NumForm::ArithExpr(Box::new(ArithExpr::Sum(vec![
                        ArithExpr::Lit(1),
                        ArithExpr::Lit(2),
                    ]))),
                )),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::B),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::S),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(2)),
            ],
            dative: vec![
                (
                    vec![AtomId(0)],
                    AtomId(1),
                    DativeBondForm::new(NumForm::ArithExpr(Box::new(ArithExpr::Sum(vec![
                        ArithExpr::Lit(0),
                        ArithExpr::Lit(1),
                    ])))),
                ),
                (vec![AtomId(3)], AtomId(2), DativeBondForm::from_order(2)),
            ],
            aromatic: vec![
                (
                    vec![AtomId(0), AtomId(1)],
                    AromaticSystemForm::from_electrons(vec![1, 2]),
                ),
                (
                    vec![AtomId(2), AtomId(3)],
                    AromaticSystemForm::from_electrons(vec![2, 1]),
                ),
            ],
            multicenter: vec![
                (
                    vec![AtomId(0), AtomId(2)],
                    MulticenterBondForm::from_electrons(vec![1, 2]),
                ),
                (
                    vec![AtomId(1), AtomId(3)],
                    MulticenterBondForm::from_electrons(vec![2, 1]),
                ),
            ],
            noncovalent: vec![
                (
                    AtomId(1),
                    AtomId(2),
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::Ionic),
                ),
                (
                    AtomId(0),
                    AtomId(3),
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                ),
            ],
            stereo_atoms: vec![
                (
                    AtomId(0),
                    vec![
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    ],
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                ),
                (
                    AtomId(1),
                    vec![
                        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                    ],
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                ),
            ],
            stereo_bonds: vec![
                (
                    BondId(0),
                    vec![
                        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    ],
                    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
                ),
                (
                    BondId(1),
                    vec![
                        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                    ],
                    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                ),
            ],
            ..Default::default()
        });
        let expected_correspondence = molecule_correspondence(&[
            vec![3, 1, 2, 0, 4, 5],
            vec![1, 0],
            vec![1, 0],
            vec![1, 0],
            vec![1, 0],
            vec![1, 0],
            vec![0, 1],
            vec![0, 1],
        ]);
        let expected = normalize_molecule(molecule.remap(&expected_correspondence)).unwrap();

        assert_eq!(
            canonicalize_constitution_with_options(
                &molecule,
                &canonicalization_context,
                CanonicalSearchOptions {
                    automorphism_pruning: true,
                    prefix_pruning: false,
                    branch_order: BranchOrder::BackendCanonical,
                },
            ),
            Ok((expected, expected_correspondence)),
        );
    }

    #[rstest]
    fn test_canonicalize_constitution_excluded_data(
        canonicalization_context: CanonicalizationContext,
    ) {
        let left = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C)
                    .with_constraint(AtomConstraintForm::Valence(NumForm::Lit(3))),
                AtomForm::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            )],
            ..Default::default()
        });
        let right = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C)
                    .with_constraint(AtomConstraintForm::Valence(NumForm::Lit(3))),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            stereo_atoms: vec![(
                AtomId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
            ..Default::default()
        });

        let (_, left_correspondence) = canonicalize_constitution_with_options(
            &left,
            &canonicalization_context,
            CanonicalSearchOptions {
                automorphism_pruning: true,
                prefix_pruning: false,
                branch_order: BranchOrder::BackendCanonical,
            },
        )
        .unwrap();
        let (_, right_correspondence) = canonicalize_constitution_with_options(
            &right,
            &canonicalization_context,
            CanonicalSearchOptions {
                automorphism_pruning: true,
                prefix_pruning: false,
                branch_order: BranchOrder::BackendCanonical,
            },
        )
        .unwrap();

        assert_eq!(right_correspondence, left_correspondence);
    }

    #[rstest]
    fn test_canonicalize_constitution_properties(
        initial_class_molecule: Molecule,
        canonicalization_context: CanonicalizationContext,
    ) {
        let normalized_source = normalize_molecule(initial_class_molecule.clone()).unwrap();
        let (canonical, correspondence) = canonicalize_constitution_with_options(
            &initial_class_molecule,
            &canonicalization_context,
            CanonicalSearchOptions {
                automorphism_pruning: true,
                prefix_pruning: false,
                branch_order: BranchOrder::BackendCanonical,
            },
        )
        .unwrap();
        let acted = normalize_molecule(initial_class_molecule.remap(&correspondence)).unwrap();
        let inverse = correspondence.reverse();

        assert_eq!(acted, canonical);
        assert!(initial_class_molecule.equiv_under(&canonical, &correspondence));
        assert!(canonical.equiv_under(&normalized_source, &inverse));
        assert_eq!(canonical.check_integrity(), Ok(()));

        let (canonical_again, _) = canonicalize_constitution_with_options(
            &canonical,
            &canonicalization_context,
            CanonicalSearchOptions {
                automorphism_pruning: true,
                prefix_pruning: false,
                branch_order: BranchOrder::BackendCanonical,
            },
        )
        .unwrap();
        let canonical_incidence = canonical.incidence_graph(IncidenceLevel::Constitution);
        let canonical_again_incidence =
            canonical_again.incidence_graph(IncidenceLevel::Constitution);
        assert_eq!(
            constitution_comparison_key(
                &canonical,
                &canonical_incidence,
                &canonical_incidence.graph().node_ids().collect::<Vec<_>>(),
            ),
            constitution_comparison_key(
                &canonical_again,
                &canonical_again_incidence,
                &canonical_again_incidence
                    .graph()
                    .node_ids()
                    .collect::<Vec<_>>(),
            ),
        );

        let renumbering = reverse_correspondence(&initial_class_molecule);
        let renumbered = initial_class_molecule.remap(&renumbering);
        let (canonical_renumbered, renumbered_correspondence) =
            canonicalize_constitution_with_options(
                &renumbered,
                &canonicalization_context,
                CanonicalSearchOptions {
                    automorphism_pruning: true,
                    prefix_pruning: false,
                    branch_order: BranchOrder::BackendCanonical,
                },
            )
            .unwrap();
        let composed = renumbering.compose(&renumbered_correspondence);
        let composed_action = normalize_molecule(initial_class_molecule.remap(&composed)).unwrap();
        let canonical_renumbered_incidence =
            canonical_renumbered.incidence_graph(IncidenceLevel::Constitution);

        assert_eq!(composed_action, canonical_renumbered);
        assert!(initial_class_molecule.equiv_under(&canonical_renumbered, &composed));
        assert_eq!(
            constitution_comparison_key(
                &canonical,
                &canonical_incidence,
                &canonical_incidence.graph().node_ids().collect::<Vec<_>>(),
            ),
            constitution_comparison_key(
                &canonical_renumbered,
                &canonical_renumbered_incidence,
                &canonical_renumbered_incidence
                    .graph()
                    .node_ids()
                    .collect::<Vec<_>>(),
            ),
        );

        let (unpruned, _) = canonicalize_constitution_with_options(
            &initial_class_molecule,
            &canonicalization_context,
            CanonicalSearchOptions {
                automorphism_pruning: false,
                prefix_pruning: false,
                branch_order: BranchOrder::ReverseNode,
            },
        )
        .unwrap();
        let unpruned_incidence = unpruned.incidence_graph(IncidenceLevel::Constitution);
        assert_eq!(
            constitution_comparison_key(
                &canonical,
                &canonical_incidence,
                &canonical_incidence.graph().node_ids().collect::<Vec<_>>(),
            ),
            constitution_comparison_key(
                &unpruned,
                &unpruned_incidence,
                &unpruned_incidence.graph().node_ids().collect::<Vec<_>>(),
            ),
        );
    }

    #[rstest]
    fn test_canonicalize_constitution_family_minimum(
        canonicalization_context: CanonicalizationContext,
    ) {
        let atoms = vec![AtomForm::from_element(Element::C); 4];
        let cases = [
            (
                "dative",
                Molecule::from_entries(MoleculeEntries {
                    atoms: atoms.clone(),
                    dative: vec![
                        (
                            vec![AtomId(0)],
                            AtomId(2),
                            DativeBondForm::new(NumForm::RangeFrom(1)),
                        ),
                        (
                            vec![AtomId(1)],
                            AtomId(3),
                            DativeBondForm::new(NumForm::RangeFrom(1)),
                        ),
                    ],
                    ..Default::default()
                }),
            ),
            (
                "aromatic",
                Molecule::from_entries(MoleculeEntries {
                    atoms: atoms.clone(),
                    aromatic: vec![
                        (
                            vec![AtomId(0), AtomId(1)],
                            AromaticSystemForm::from_electrons(vec![1, 2])
                                .with_charge(NumForm::var("q")),
                        ),
                        (
                            vec![AtomId(2), AtomId(3)],
                            AromaticSystemForm::from_electrons(vec![1, 2])
                                .with_charge(NumForm::var("q")),
                        ),
                    ],
                    ..Default::default()
                }),
            ),
            (
                "multicenter",
                Molecule::from_entries(MoleculeEntries {
                    atoms: atoms.clone(),
                    multicenter: vec![
                        (
                            vec![AtomId(0), AtomId(1)],
                            MulticenterBondForm::from_electrons(vec![2, 1]),
                        ),
                        (
                            vec![AtomId(2), AtomId(3)],
                            MulticenterBondForm::from_electrons(vec![2, 1]),
                        ),
                    ],
                    ..Default::default()
                }),
            ),
            (
                "noncovalent",
                Molecule::from_entries(MoleculeEntries {
                    atoms,
                    noncovalent: vec![
                        (
                            AtomId(0),
                            AtomId(1),
                            NoncovalentBondForm::new(NoncovalentBondKindForm::Undetermined),
                        ),
                        (
                            AtomId(2),
                            AtomId(3),
                            NoncovalentBondForm::new(NoncovalentBondKindForm::Undetermined),
                        ),
                    ],
                    ..Default::default()
                }),
            ),
        ];

        for (family, molecule) in cases {
            let incidence_graph = molecule.incidence_graph(IncidenceLevel::Constitution);
            let (entity_keys, incidence_keys) =
                initial_class_keys(&molecule, &incidence_graph).unwrap();
            let classes = rank_initial_classes(&entity_keys, &incidence_keys);
            let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
            let leaf_candidate = |order: &[NodeId]| {
                constitution_candidate(&molecule, &incidence_graph, order).unwrap()
            };
            let expected = exhaustive_minimum(&adapter, &leaf_candidate);
            let (canonical, correspondence) = canonicalize_constitution_with_options(
                &molecule,
                &canonicalization_context,
                CanonicalSearchOptions {
                    automorphism_pruning: true,
                    prefix_pruning: false,
                    branch_order: BranchOrder::BackendCanonical,
                },
            )
            .unwrap();
            let (unpruned, _) = canonicalize_constitution_with_options(
                &molecule,
                &canonicalization_context,
                CanonicalSearchOptions {
                    automorphism_pruning: false,
                    prefix_pruning: false,
                    branch_order: BranchOrder::ReverseNode,
                },
            )
            .unwrap();
            let canonical_incidence = canonical.incidence_graph(IncidenceLevel::Constitution);
            let canonical_order = canonical_incidence.graph().node_ids().collect::<Vec<_>>();
            let unpruned_incidence = unpruned.incidence_graph(IncidenceLevel::Constitution);
            let unpruned_order = unpruned_incidence.graph().node_ids().collect::<Vec<_>>();

            assert_eq!(
                constitution_comparison_key(&unpruned, &unpruned_incidence, &unpruned_order),
                Ok(expected.key.clone()),
                "unpruned {family}",
            );
            assert_eq!(
                constitution_comparison_key(&canonical, &canonical_incidence, &canonical_order),
                Ok(expected.key),
                "pruned {family}",
            );
            assert_eq!(unpruned, canonical, "{family}");
            assert!(
                molecule.equiv_under(&canonical, &correspondence),
                "{family}"
            );
            assert_eq!(canonical.check_integrity(), Ok(()), "{family}");

            for (index, atom_images) in permutations(molecule.atoms().count())
                .into_iter()
                .enumerate()
            {
                let mut images =
                    molecule_counts(&molecule).map(|count| (0..count).collect::<Vec<_>>());
                images[0] = atom_images;
                if index % 2 == 1 {
                    for family_images in &mut images[1..6] {
                        family_images.reverse();
                    }
                }
                let renumbered = molecule.remap(&molecule_correspondence(&images));

                assert_eq!(
                    canonicalize_constitution(&renumbered, &canonicalization_context),
                    Ok(canonical.clone()),
                    "{family}, renumbering {index}",
                );
            }
        }
    }

    #[rstest]
    fn test_canonicalize_constitution_participant_order(
        canonicalization_context: CanonicalizationContext,
    ) {
        let left = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 6],
            aromatic: vec![
                (
                    vec![AtomId(0), AtomId(1), AtomId(2)],
                    AromaticSystemForm::from_electrons(vec![1, 2, 3]),
                ),
                (
                    vec![AtomId(3), AtomId(4), AtomId(5)],
                    AromaticSystemForm::from_electrons(vec![3, 2, 1]),
                ),
            ],
            multicenter: vec![
                (
                    vec![AtomId(0), AtomId(3), AtomId(5)],
                    MulticenterBondForm::from_electrons(vec![1, 2, 3]),
                ),
                (
                    vec![AtomId(1), AtomId(2), AtomId(4)],
                    MulticenterBondForm::from_electrons(vec![3, 2, 1]),
                ),
            ],
            ..Default::default()
        });
        let right = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 6],
            aromatic: vec![
                (
                    vec![AtomId(2), AtomId(0), AtomId(1)],
                    AromaticSystemForm::from_electrons(vec![3, 1, 2]),
                ),
                (
                    vec![AtomId(4), AtomId(5), AtomId(3)],
                    AromaticSystemForm::from_electrons(vec![2, 1, 3]),
                ),
            ],
            multicenter: vec![
                (
                    vec![AtomId(5), AtomId(0), AtomId(3)],
                    MulticenterBondForm::from_electrons(vec![3, 1, 2]),
                ),
                (
                    vec![AtomId(2), AtomId(4), AtomId(1)],
                    MulticenterBondForm::from_electrons(vec![2, 1, 3]),
                ),
            ],
            ..Default::default()
        });

        assert_eq!(
            canonicalize_constitution(&right, &canonicalization_context),
            canonicalize_constitution(&left, &canonicalization_context),
        );
    }

    #[rstest]
    fn test_canonicalize_constitution_contradiction(
        canonicalization_context: CanonicalizationContext,
    ) {
        let selected = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 2],
            dative: vec![(
                vec![AtomId(0)],
                AtomId(1),
                DativeBondForm::new(NumForm::lit_set([])),
            )],
            ..Default::default()
        });
        let excluded_constraint = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)
                .with_constraint(AtomConstraintForm::Valence(NumForm::lit_set([])))],
            ..Default::default()
        });
        let excluded_stereo = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 2],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::lit_set([])),
            )],
            ..Default::default()
        });

        for (location, molecule) in [
            ("selected", selected),
            ("excluded constraint", excluded_constraint),
            ("excluded stereo", excluded_stereo),
        ] {
            assert_eq!(
                canonicalize_constitution(&molecule, &canonicalization_context),
                Err(MoleculeCanonicalizationError::Contradiction(Contradiction)),
                "{location}",
            );
        }
    }

    #[rstest]
    fn test_canonical_search_prefix() {
        let source = Graph::new(4, &[]);
        let adapter = direct_graph_adapter(&source);
        let leaf_candidate = |order: &[NodeId]| CanonicalCandidate {
            key: order.to_vec(),
            entity_order: order.to_vec(),
        };
        let prefix_worse = |partition: &OrderedPartition,
                            best: &CanonicalCandidate<Vec<NodeId>>| {
            let prefix = partition.fixed_entity_prefix(4);
            prefix.as_slice() > &best.key[..prefix.len()]
        };
        let expected = exhaustive_minimum(&adapter, &leaf_candidate);
        let actual = canonical_search(
            &adapter,
            &adapter.classes,
            AutomorphismAlgorithm::Nauty,
            CanonicalSearchOptions {
                automorphism_pruning: false,
                prefix_pruning: true,
                branch_order: BranchOrder::Node,
            },
            &leaf_candidate,
            &prefix_worse,
        );

        assert_eq!(actual.candidate.key, expected.key);
        assert_ne!(actual.stats.prefix_pruned_branches, 0);
    }

    #[rstest]
    #[case::order_four(4)]
    fn test_canonical_search_exhaustive(#[case] node_count: usize) {
        let endpoint_pairs = (0..node_count as u32)
            .flat_map(|first| ((first + 1)..node_count as u32).map(move |second| [first, second]))
            .collect::<Vec<_>>();

        for edge_mask in 0..(1_u64 << endpoint_pairs.len()) {
            let edges = endpoint_pairs
                .iter()
                .enumerate()
                .filter_map(|(position, &edge)| ((edge_mask >> position) & 1 == 1).then_some(edge))
                .collect::<Vec<_>>();
            let source = Graph::new(node_count, &edges);
            let adapter = direct_graph_adapter(&source);
            let leaf_candidate = |order: &[NodeId]| {
                let mut positions = vec![0_u32; node_count];
                for (position, node) in order.iter().enumerate() {
                    positions[node.index()] = position as u32;
                }
                let mut mapped_edges = source
                    .edge_ids()
                    .map(|edge| {
                        let [first, second] = source.edge_endpoints(edge);
                        let first = positions[first.index()];
                        let second = positions[second.index()];
                        [first.min(second), first.max(second)]
                    })
                    .collect::<Vec<_>>();
                mapped_edges.sort_unstable();
                CanonicalCandidate {
                    key: mapped_edges,
                    entity_order: order.to_vec(),
                }
            };
            let no_prefix = |_: &OrderedPartition, _: &CanonicalCandidate<_>| false;
            let expected = exhaustive_minimum(&adapter, &leaf_candidate);

            for options in [
                CanonicalSearchOptions {
                    automorphism_pruning: false,
                    prefix_pruning: false,
                    branch_order: BranchOrder::ReverseNode,
                },
                CanonicalSearchOptions {
                    automorphism_pruning: true,
                    prefix_pruning: false,
                    branch_order: BranchOrder::BackendCanonical,
                },
            ] {
                assert_eq!(
                    canonical_search(
                        &adapter,
                        &adapter.classes,
                        AutomorphismAlgorithm::Nauty,
                        options,
                        &leaf_candidate,
                        &no_prefix,
                    )
                    .candidate
                    .key,
                    expected.key,
                    "edge mask {edge_mask:#08b}",
                );
            }
        }
    }

    #[rstest]
    #[case::topology(IncidenceLevel::Topology)]
    #[case::constitution(IncidenceLevel::Constitution)]
    #[case::full(IncidenceLevel::Full)]
    fn test_colored_encoding_dense_remapping_equivalence(#[case] level: IncidenceLevel) {
        let entries = encoding_entries();
        let complete = Molecule::from_entries(entries.clone());
        let molecule = Molecule::from_entries(project_entries(entries, level));
        let remapped = molecule.remap(&reverse_correspondence(&molecule));

        assert!(colored_encoding_equivalent(&complete, &molecule, level));
        assert_eq!(
            colored_encoding_equivalent(&molecule, &remapped, level),
            explicitly_dense_equivalent(&molecule, &remapped),
        );
        assert!(explicitly_dense_equivalent(&molecule, &remapped));

        let mut distinguished = remapped;
        distinguished.atom_mut(AtomId(0)).attributes.element = ElementForm::Lit(Element::O);
        assert_eq!(
            colored_encoding_equivalent(&molecule, &distinguished, level),
            explicitly_dense_equivalent(&molecule, &distinguished),
        );
        assert!(!explicitly_dense_equivalent(&molecule, &distinguished));
    }

    #[rstest]
    #[case::order_four(4)]
    fn test_colored_encoding_exhaustive_graph_domain(#[case] atom_count: usize) {
        let endpoint_pairs = (0..atom_count as u32)
            .flat_map(|first| ((first + 1)..atom_count as u32).map(move |second| [first, second]))
            .collect::<Vec<_>>();

        for edge_mask in 0..(1_u64 << endpoint_pairs.len()) {
            let bonds = endpoint_pairs
                .iter()
                .enumerate()
                .filter(|(position, _)| (edge_mask >> position) & 1 == 1)
                .map(|(_, &[first, second])| {
                    (AtomId(first), AtomId(second), BondForm::from_order(1))
                })
                .collect();
            let molecule = Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); atom_count],
                bonds,
                ..Default::default()
            });
            let remapped = molecule.remap(&reverse_correspondence(&molecule));
            assert_eq!(
                colored_encoding_equivalent(&molecule, &remapped, IncidenceLevel::Topology),
                explicitly_dense_equivalent(&molecule, &remapped),
                "edge mask {edge_mask:#08b}",
            );

            let mut distinguished = remapped;
            distinguished.atom_mut(AtomId(0)).attributes.element = ElementForm::Lit(Element::O);
            assert_eq!(
                colored_encoding_equivalent(&molecule, &distinguished, IncidenceLevel::Topology,),
                explicitly_dense_equivalent(&molecule, &distinguished),
                "edge mask {edge_mask:#08b}",
            );
        }
    }

    #[rstest]
    #[ignore = "manual optimized-build canonicalization benchmark"]
    fn benchmark_exact_canonicalization_carrier() {
        const ITERATIONS: u32 = 100;

        eprintln!(
            "case\tlevel\tincidence\tadapter\tincidence_ns\tclasses_ns\tadapter_ns\tbackend_ns\tsearch_ns\tremap_ns\trefinements\tleaves\tprefix_pruned\torbit_pruned"
        );
        for case in benchmark_cases::corpus() {
            for level in benchmark_cases::LEVELS {
                let mut incidence_time = Duration::ZERO;
                let mut classes_time = Duration::ZERO;
                let mut adapter_time = Duration::ZERO;
                let mut backend_time = Duration::ZERO;
                let mut search_time = Duration::ZERO;
                let mut remap_time = Duration::ZERO;
                let mut sizes = (0, 0, 0, 0);
                let mut stats = CanonicalSearchStats::default();

                for _ in 0..ITERATIONS {
                    let start = Instant::now();
                    let incidence = case.molecule.incidence_graph(level);
                    incidence_time += start.elapsed();

                    let start = Instant::now();
                    let (entity_keys, incidence_keys) =
                        initial_class_keys(&case.molecule, &incidence).unwrap();
                    let classes = rank_initial_classes(&entity_keys, &incidence_keys);
                    classes_time += start.elapsed();

                    let start = Instant::now();
                    let adapter = AutomorphismAdapter::new(&incidence, &classes);
                    let descriptors =
                        partition_descriptors(&adapter, &entity_keys, &incidence_keys);
                    adapter_time += start.elapsed();

                    let start = Instant::now();
                    black_box(adapter.automorphisms(AutomorphismAlgorithm::Nauty));
                    backend_time += start.elapsed();

                    let source = incidence.graph();
                    let leaf_candidate = |order: &[NodeId]| CanonicalCandidate {
                        key: structural_leaf_key(order, source, &classes),
                        entity_order: order.to_vec(),
                    };
                    let no_prefix = |_: &OrderedPartition, _: &CanonicalCandidate<_>| false;
                    let start = Instant::now();
                    let result = canonical_search(
                        &adapter,
                        &descriptors,
                        AutomorphismAlgorithm::Nauty,
                        CanonicalSearchOptions {
                            automorphism_pruning: true,
                            prefix_pruning: false,
                            branch_order: BranchOrder::BackendCanonical,
                        },
                        &leaf_candidate,
                        &no_prefix,
                    );
                    search_time += start.elapsed();

                    let correspondence = correspondence_from_order(
                        &case.molecule,
                        &incidence,
                        &result.candidate.entity_order,
                    );
                    let start = Instant::now();
                    black_box(case.molecule.remap(&correspondence));
                    remap_time += start.elapsed();

                    sizes = (
                        incidence.graph().node_count(),
                        incidence.graph().edge_count(),
                        adapter.graph().node_count(),
                        adapter.graph().edge_count(),
                    );
                    stats = result.stats;
                }

                let average = |duration: Duration| duration.as_nanos() / u128::from(ITERATIONS);
                eprintln!(
                    "{}\t{}\tn{}_e{}\tn{}_e{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    case.name,
                    benchmark_cases::level_name(level),
                    sizes.0,
                    sizes.1,
                    sizes.2,
                    sizes.3,
                    average(incidence_time),
                    average(classes_time),
                    average(adapter_time),
                    average(backend_time),
                    average(search_time),
                    average(remap_time),
                    stats.refinement_calls,
                    stats.visited_leaves,
                    stats.prefix_pruned_branches,
                    stats.orbit_pruned_branches,
                );
            }
        }
    }

    #[rstest]
    fn test_initial_classes_error() {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C).with_charge(NumForm::LitSet(Box::default()))
            ],
            ..Default::default()
        });
        let incidence_graph = molecule.incidence_graph(IncidenceLevel::Topology);

        assert_eq!(
            initial_classes(&molecule, &incidence_graph),
            Err(Contradiction)
        );
    }

    #[rstest]
    #[case::integrity(
        MoleculeCanonicalizationError::from(MoleculeIntegrityError::InvalidReference {
            entity: Entity::Atom(AtomId(1)),
        }),
        MoleculeCanonicalizationError::Integrity(MoleculeIntegrityError::InvalidReference {
            entity: Entity::Atom(AtomId(1)),
        }),
    )]
    #[case::contradiction(
        MoleculeCanonicalizationError::from(Contradiction),
        MoleculeCanonicalizationError::Contradiction(Contradiction)
    )]
    fn test_molecule_canonicalization_error_from(
        #[case] actual: MoleculeCanonicalizationError,
        #[case] expected: MoleculeCanonicalizationError,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::integrity(
        ReactionSpanCanonicalizationError::from(ReactionSpanIntegrityError::InvalidReference {
            entity: Entity::Atom(AtomId(1)),
        }),
        ReactionSpanCanonicalizationError::Integrity(
            ReactionSpanIntegrityError::InvalidReference {
                entity: Entity::Atom(AtomId(1)),
            },
        ),
    )]
    #[case::contradiction(
        ReactionSpanCanonicalizationError::from(Contradiction),
        ReactionSpanCanonicalizationError::Contradiction(Contradiction)
    )]
    fn test_reaction_span_canonicalization_error_from(
        #[case] actual: ReactionSpanCanonicalizationError,
        #[case] expected: ReactionSpanCanonicalizationError,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::integrity(
        ReactionCanonicalizationError::from(ReactionIntegrityError::InvalidReference {
            entity: Entity::Atom(AtomId(1)),
        }),
        ReactionCanonicalizationError::Integrity(ReactionIntegrityError::InvalidReference {
            entity: Entity::Atom(AtomId(1)),
        }),
    )]
    #[case::contradiction(
        ReactionCanonicalizationError::from(Contradiction),
        ReactionCanonicalizationError::Contradiction(Contradiction)
    )]
    fn test_reaction_canonicalization_error_from(
        #[case] actual: ReactionCanonicalizationError,
        #[case] expected: ReactionCanonicalizationError,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::topology(StructuralDomainPosition::TOPOLOGY, 0)]
    #[case::non_stereo(StructuralDomainPosition::NON_STEREO, 1)]
    #[case::stereo(StructuralDomainPosition::STEREO, 2)]
    fn test_structural_domain_position(
        #[case] position: StructuralDomainPosition,
        #[case] expected: u16,
    ) {
        assert_eq!(position.0, expected);
    }

    #[rstest]
    #[case::atom(
        EntityBlockPosition::ATOM,
        EntityBlockPosition::new(StructuralDomainPosition::TOPOLOGY, 0)
    )]
    #[case::bond(
        EntityBlockPosition::BOND,
        EntityBlockPosition::new(StructuralDomainPosition::TOPOLOGY, 1)
    )]
    #[case::dative_bond(
        EntityBlockPosition::DATIVE_BOND,
        EntityBlockPosition::new(StructuralDomainPosition::NON_STEREO, 0)
    )]
    #[case::aromatic_system(
        EntityBlockPosition::AROMATIC_SYSTEM,
        EntityBlockPosition::new(StructuralDomainPosition::NON_STEREO, 1)
    )]
    #[case::multicenter_bond(
        EntityBlockPosition::MULTICENTER_BOND,
        EntityBlockPosition::new(StructuralDomainPosition::NON_STEREO, 2)
    )]
    #[case::noncovalent_bond(
        EntityBlockPosition::NONCOVALENT_BOND,
        EntityBlockPosition::new(StructuralDomainPosition::NON_STEREO, 3)
    )]
    #[case::stereo_atom(
        EntityBlockPosition::STEREO_ATOM,
        EntityBlockPosition::new(StructuralDomainPosition::STEREO, 0)
    )]
    #[case::stereo_bond(
        EntityBlockPosition::STEREO_BOND,
        EntityBlockPosition::new(StructuralDomainPosition::STEREO, 1)
    )]
    fn test_entity_block_position(
        #[case] position: EntityBlockPosition,
        #[case] expected: EntityBlockPosition,
    ) {
        assert_eq!(position, expected);
    }

    #[rstest]
    #[case::topology_slots(EntityBlockPosition::ATOM, EntityBlockPosition::BOND)]
    #[case::topology_before_non_stereo(
        EntityBlockPosition::new(StructuralDomainPosition::TOPOLOGY, 2),
        EntityBlockPosition::DATIVE_BOND
    )]
    #[case::non_stereo_slots(
        EntityBlockPosition::DATIVE_BOND,
        EntityBlockPosition::AROMATIC_SYSTEM
    )]
    #[case::non_stereo_before_stereo(
        EntityBlockPosition::new(StructuralDomainPosition::NON_STEREO, 4),
        EntityBlockPosition::STEREO_ATOM
    )]
    #[case::stereo_slots(EntityBlockPosition::STEREO_ATOM, EntityBlockPosition::STEREO_BOND)]
    fn test_entity_block_position_cmp(
        #[case] lhs: EntityBlockPosition,
        #[case] rhs: EntityBlockPosition,
    ) {
        assert_eq!(lhs.cmp(&rhs), Ordering::Less);
    }

    #[rstest]
    #[case::atom(
        ConstraintBlockPosition::ATOM,
        ConstraintBlockPosition::Inline(EntityBlockPosition::ATOM)
    )]
    #[case::bond(
        ConstraintBlockPosition::BOND,
        ConstraintBlockPosition::Inline(EntityBlockPosition::BOND)
    )]
    #[case::dative_bond(
        ConstraintBlockPosition::DATIVE_BOND,
        ConstraintBlockPosition::Inline(EntityBlockPosition::DATIVE_BOND)
    )]
    #[case::aromatic_system(
        ConstraintBlockPosition::AROMATIC_SYSTEM,
        ConstraintBlockPosition::Inline(EntityBlockPosition::AROMATIC_SYSTEM)
    )]
    #[case::multicenter_bond(
        ConstraintBlockPosition::MULTICENTER_BOND,
        ConstraintBlockPosition::Inline(EntityBlockPosition::MULTICENTER_BOND)
    )]
    #[case::noncovalent_bond(
        ConstraintBlockPosition::NONCOVALENT_BOND,
        ConstraintBlockPosition::Inline(EntityBlockPosition::NONCOVALENT_BOND)
    )]
    #[case::stereo_atom(
        ConstraintBlockPosition::STEREO_ATOM,
        ConstraintBlockPosition::Inline(EntityBlockPosition::STEREO_ATOM)
    )]
    #[case::stereo_bond(
        ConstraintBlockPosition::STEREO_BOND,
        ConstraintBlockPosition::Inline(EntityBlockPosition::STEREO_BOND)
    )]
    #[case::molecule(ConstraintBlockPosition::MOLECULE, ConstraintBlockPosition::Molecule)]
    fn test_constraint_block_position(
        #[case] position: ConstraintBlockPosition,
        #[case] expected: ConstraintBlockPosition,
    ) {
        assert_eq!(position, expected);
    }

    #[rstest]
    #[case::unchanged(SpanTagPosition::UNCHANGED, 0)]
    #[case::added(SpanTagPosition::ADDED, 1)]
    #[case::removed(SpanTagPosition::REMOVED, 2)]
    #[case::modified(SpanTagPosition::MODIFIED, 3)]
    fn test_span_tag_position(#[case] position: SpanTagPosition, #[case] expected: u16) {
        assert_eq!(position.0, expected);
    }

    #[rstest]
    #[case::unit(CanonicalKeyValue::Unit, 0)]
    #[case::boolean(CanonicalKeyValue::Bool(false), 1)]
    #[case::unsigned(CanonicalKeyValue::Unsigned(0), 2)]
    #[case::signed(CanonicalKeyValue::Signed(0), 3)]
    #[case::text(CanonicalKeyValue::Text(String::new()), 4)]
    #[case::sequence(CanonicalKeyValue::Sequence(Vec::new()), 5)]
    #[case::product(CanonicalKeyValue::Product(Vec::new()), 6)]
    #[case::variant(CanonicalKeyValue::Variant(VariantKey {
        position: VariantPosition(0),
        fields: Vec::new(),
    }), 7)]
    #[case::span(CanonicalKeyValue::Span(SpanKey {
        position: SpanTagPosition::UNCHANGED,
        values: Vec::new(),
    }), 8)]
    fn test_canonical_key_value_position(#[case] value: CanonicalKeyValue, #[case] expected: u16) {
        assert_eq!(value.position(), expected);
    }

    #[rstest]
    #[case::position_precedes_payload(
        PositionedKey {
            position: FieldPosition(0),
            value: CanonicalKeyValue::Signed(10),
        },
        PositionedKey {
            position: FieldPosition(1),
            value: CanonicalKeyValue::Signed(-10),
        },
        Ordering::Less,
    )]
    #[case::payload_breaks_position_tie(
        PositionedKey {
            position: FieldPosition(0),
            value: CanonicalKeyValue::Signed(10),
        },
        PositionedKey {
            position: FieldPosition(0),
            value: CanonicalKeyValue::Signed(-10),
        },
        Ordering::Greater,
    )]
    fn test_positioned_key_cmp(
        #[case] lhs: FieldKey,
        #[case] rhs: FieldKey,
        #[case] expected: Ordering,
    ) {
        assert_eq!(lhs.cmp(&rhs), expected);
    }

    #[rstest]
    #[case::span_tag_precedes_value(
        SpanKey {
            position: SpanTagPosition::UNCHANGED,
            values: vec![CanonicalKeyValue::Signed(10)],
        },
        SpanKey {
            position: SpanTagPosition::ADDED,
            values: vec![CanonicalKeyValue::Signed(-10)],
        },
        Ordering::Less,
    )]
    #[case::lhs_precedes_rhs(
        SpanKey {
            position: SpanTagPosition::MODIFIED,
            values: vec![CanonicalKeyValue::Signed(0), CanonicalKeyValue::Signed(10)],
        },
        SpanKey {
            position: SpanTagPosition::MODIFIED,
            values: vec![CanonicalKeyValue::Signed(1), CanonicalKeyValue::Signed(-10)],
        },
        Ordering::Less,
    )]
    fn test_span_key_cmp(#[case] lhs: SpanKey, #[case] rhs: SpanKey, #[case] expected: Ordering) {
        assert_eq!(lhs.cmp(&rhs), expected);
    }

    #[rstest]
    #[case::absent_extension_preserves_key(None, Ordering::Equal)]
    #[case::present_extension_appends_field(
        Some(PositionedKey {
            position: FieldPosition(8),
            value: CanonicalKeyValue::Unsigned(1),
        }),
        Ordering::Less,
    )]
    fn test_canonical_key_value_append_only_extension(
        #[case] extension: Option<FieldKey>,
        #[case] expected: Ordering,
    ) {
        let original = CanonicalKeyValue::Product(vec![PositionedKey {
            position: FieldPosition(0),
            value: CanonicalKeyValue::Signed(1),
        }]);
        let mut extended_fields = match &original {
            CanonicalKeyValue::Product(fields) => fields.clone(),
            _ => unreachable!(),
        };
        extended_fields.extend(extension);
        let extended = CanonicalKeyValue::Product(extended_fields);

        assert_eq!(original.cmp(&extended), expected);
    }

    #[rstest]
    #[case::entity_blocks_precede_constraints(
        CanonicalComparisonKey {
            entity_blocks: vec![PositionedKey {
                position: EntityBlockPosition::ATOM,
                value: CanonicalKeyValue::Unsigned(0),
            }],
            constraints: vec![PositionedKey {
                position: ConstraintBlockPosition::ATOM,
                value: CanonicalKeyValue::Unsigned(10),
            }],
        },
        CanonicalComparisonKey {
            entity_blocks: vec![PositionedKey {
                position: EntityBlockPosition::ATOM,
                value: CanonicalKeyValue::Unsigned(1),
            }],
            constraints: vec![PositionedKey {
                position: ConstraintBlockPosition::ATOM,
                value: CanonicalKeyValue::Unsigned(0),
            }],
        },
        Ordering::Less,
    )]
    #[case::constraints_break_entity_tie(
        CanonicalComparisonKey {
            entity_blocks: Vec::new(),
            constraints: vec![PositionedKey {
                position: ConstraintBlockPosition::ATOM,
                value: CanonicalKeyValue::Unsigned(0),
            }],
        },
        CanonicalComparisonKey {
            entity_blocks: Vec::new(),
            constraints: vec![PositionedKey {
                position: ConstraintBlockPosition::BOND,
                value: CanonicalKeyValue::Unsigned(0),
            }],
        },
        Ordering::Less,
    )]
    fn test_canonical_comparison_key_cmp(
        #[case] lhs: CanonicalComparisonKey,
        #[case] rhs: CanonicalComparisonKey,
        #[case] expected: Ordering,
    ) {
        assert_eq!(lhs.cmp(&rhs), expected);
    }

    #[rstest]
    #[case::non_stereo_before_stereo(
        ConstraintBlockPosition::Inline(EntityBlockPosition::new(
            StructuralDomainPosition::NON_STEREO,
            4,
        )),
        ConstraintBlockPosition::STEREO_ATOM
    )]
    #[case::stereo_before_molecule(
        ConstraintBlockPosition::STEREO_BOND,
        ConstraintBlockPosition::MOLECULE
    )]
    fn test_constraint_block_position_cmp(
        #[case] lhs: ConstraintBlockPosition,
        #[case] rhs: ConstraintBlockPosition,
    ) {
        assert_eq!(lhs.cmp(&rhs), Ordering::Less);
    }

    #[rstest]
    #[case::local_slot(
        RelationalConstraintPosition::new(EntityBlockPosition::DATIVE_BOND, 0),
        RelationalConstraintPosition::new(EntityBlockPosition::DATIVE_BOND, 1)
    )]
    #[case::entity_slot(
        RelationalConstraintPosition::new(EntityBlockPosition::DATIVE_BOND, 7),
        RelationalConstraintPosition::new(EntityBlockPosition::AROMATIC_SYSTEM, 0)
    )]
    #[case::domain(
        RelationalConstraintPosition::new(
            EntityBlockPosition::new(StructuralDomainPosition::NON_STEREO, 4),
            0,
        ),
        RelationalConstraintPosition::new(EntityBlockPosition::STEREO_ATOM, 0)
    )]
    fn test_relational_constraint_position_cmp(
        #[case] lhs: RelationalConstraintPosition,
        #[case] rhs: RelationalConstraintPosition,
    ) {
        assert_eq!(lhs.cmp(&rhs), Ordering::Less);
    }
}
