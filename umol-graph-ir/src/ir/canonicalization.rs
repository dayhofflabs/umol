//! Aggregate canonicalization inputs and failures.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

use thiserror::Error;
use umol_graph_core::{
    AutomorphismAlgorithm, AutomorphismOutput, EdgeId, Graph, NodeId, SubdividedGraph,
    SubdivisionNodeSource,
};

use super::atom::{ElementForm, IsotopeMassForm};
use super::entity::Entity;
use super::error::Contradiction;
use super::incidence::{Incidence, IncidenceGraph};
use super::ligand::StereoLigandKind;
use super::molecule::{Molecule, MoleculeIntegrityError};
use super::noncovalent::{NoncovalentBondKind, NoncovalentBondKindForm};
use super::num::{ArithExpr, NumForm, PredExpr};
use super::operators::{MemOp, RelOp};
use super::reaction_span::ReactionSpanIntegrityError;
use super::spin::UnpairedElectronsForm;
use super::stereo::StereoKind;
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

/// Structural level used to select or compare an aggregate's entity frame.
///
/// The variants form a nested hierarchy. [`Topology`](Self::Topology) contains atoms and localized
/// bonds. [`Constitution`](Self::Constitution) additionally contains dative bonds, aromatic
/// systems, multicenter bonds, and noncovalent bonds. [`Full`](Self::Full) additionally contains
/// stereo atoms and stereo bonds. In entity-domain terminology, topology is AB, non-stereo is
/// DAMN, stereo is SS, constitution is topology plus non-stereo, and overlays are non-stereo plus
/// stereo. Constraints do not contribute to any structural level.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CanonicalizationLevel {
    Topology,
    Constitution,
    Full,
}

#[allow(dead_code)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct StructuralDomainPosition(u16);

#[allow(dead_code)]
impl StructuralDomainPosition {
    const TOPOLOGY: Self = Self(0);
    const NON_STEREO: Self = Self(1);
    const STEREO: Self = Self(2);
}

#[allow(dead_code)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct EntityBlockPosition {
    domain: StructuralDomainPosition,
    slot: u16,
}

#[allow(dead_code)]
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

#[allow(dead_code)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct FieldPosition(u16);

#[allow(dead_code)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct VariantPosition(u16);

#[allow(dead_code)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct SpanTagPosition(u16);

#[allow(dead_code)]
impl SpanTagPosition {
    const UNCHANGED: Self = Self(0);
    const ADDED: Self = Self(1);
    const REMOVED: Self = Self(2);
    const MODIFIED: Self = Self(3);
}

#[allow(dead_code)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ConstraintBlockPosition {
    Inline(EntityBlockPosition),
    Molecule,
}

#[allow(dead_code)]
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

#[allow(dead_code)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct RelationalConstraintPosition {
    entity: EntityBlockPosition,
    slot: u16,
}

#[allow(dead_code)]
impl RelationalConstraintPosition {
    const fn new(entity: EntityBlockPosition, slot: u16) -> Self {
        Self { entity, slot }
    }
}

#[allow(dead_code)]
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

#[allow(dead_code)]
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

#[allow(dead_code)]
#[derive(Clone, Debug, PartialEq, Eq)]
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

#[allow(dead_code)]
#[derive(Clone, Debug, PartialEq, Eq)]
enum CanonicalKeyValue {
    Unit,
    Bool(bool),
    Unsigned(u64),
    Signed(i64),
    Text(String),
    Sequence(Vec<Self>),
    Product(Vec<FieldKey>),
    Variant(VariantKey),
    Span(SpanKey),
}

#[allow(dead_code)]
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

#[allow(dead_code)]
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

#[allow(dead_code)]
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

#[allow(dead_code)]
#[derive(Clone, Debug, PartialEq, Eq)]
struct InitialClasses {
    entities: Vec<u32>,
    incidences: Vec<u32>,
}

#[allow(dead_code)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum AutomorphismClass {
    Entity(u32),
    Incidence(u32),
}

#[allow(dead_code)]
#[derive(Clone, Debug)]
struct AutomorphismAdapter {
    // Source entity nodes retain their ids; each source incidence edge becomes one additional node.
    subdivision: SubdividedGraph,
    classes: Vec<AutomorphismClass>,
    source_node_count: usize,
    entity_blocks: Vec<Vec<NodeId>>,
}

#[allow(dead_code)]
#[derive(Clone, Debug, PartialEq, Eq)]
struct ProjectedAutomorphismOutput {
    orbits: Vec<NodeId>,
    // Backend canonical labels are branch-order hints, not the canonical molecule numbering.
    canonical_labels: Vec<NodeId>,
    generators: Vec<Vec<NodeId>>,
}

#[allow(dead_code)]
impl AutomorphismAdapter {
    fn new(incidence_graph: &IncidenceGraph, initial_classes: &InitialClasses) -> Self {
        let source = incidence_graph.graph();
        debug_assert_eq!(initial_classes.entities.len(), source.node_count());
        debug_assert_eq!(initial_classes.incidences.len(), source.edge_count());

        let subdivision = source.subdivide_edges();
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
        let classes = initial_classes
            .entities
            .iter()
            .copied()
            .map(AutomorphismClass::Entity)
            .chain(
                initial_classes
                    .incidences
                    .iter()
                    .copied()
                    .map(AutomorphismClass::Incidence),
            )
            .collect();

        Self {
            subdivision,
            classes,
            source_node_count: source.node_count(),
            entity_blocks,
        }
    }

    fn graph(&self) -> &Graph {
        self.subdivision.graph()
    }

    fn class(&self, node: NodeId) -> AutomorphismClass {
        self.classes[node.index()]
    }

    fn node_source(&self, node: NodeId) -> SubdivisionNodeSource {
        self.subdivision.node_source(node)
    }

    fn node_of(&self, source: SubdivisionNodeSource) -> NodeId {
        self.subdivision.node_of(source)
    }

    fn edge_source(&self, edge: EdgeId) -> EdgeId {
        self.subdivision.edge_source(edge)
    }

    fn incidence_edges_of(&self, edge: EdgeId) -> [EdgeId; 2] {
        self.subdivision.incidence_edges_of(edge)
    }

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
        let output = self
            .graph()
            .automorphisms(|node| cell_indices[node.index()], algorithm);

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

#[allow(dead_code)]
#[derive(Clone, Debug, PartialEq, Eq)]
struct OrderedPartition {
    cells: Vec<Vec<NodeId>>,
}

#[allow(dead_code)]
impl OrderedPartition {
    fn from_classes<C: Copy + Ord>(classes: &[C]) -> Self {
        let mut cells = BTreeMap::<C, Vec<NodeId>>::new();
        for (index, &class) in classes.iter().enumerate() {
            cells.entry(class).or_default().push(NodeId(index as u32));
        }

        Self {
            cells: cells.into_values().collect(),
        }
    }

    fn refine(mut self, graph: &Graph) -> Self {
        loop {
            let cell_indices = self.cell_indices(graph.node_count());
            let mut refined = Vec::with_capacity(self.cells.len());
            let mut changed = false;

            for cell in self.cells {
                let mut splits = BTreeMap::<Vec<u32>, Vec<NodeId>>::new();
                for node in cell {
                    let mut signature = graph
                        .neighbors(node)
                        .iter()
                        .map(|neighbor| cell_indices[neighbor.node.index()])
                        .collect::<Vec<_>>();
                    signature.sort_unstable();
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

#[allow(dead_code)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum BranchOrder {
    Node,
    ReverseNode,
    BackendCanonical,
}

#[allow(dead_code)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct CanonicalSearchOptions {
    automorphism_pruning: bool,
    prefix_pruning: bool,
    branch_order: BranchOrder,
}

#[allow(dead_code)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct CanonicalSearchStats {
    refinement_calls: usize,
    visited_leaves: usize,
    prefix_pruned_branches: usize,
    orbit_pruned_branches: usize,
}

#[allow(dead_code)]
#[derive(Clone, Debug, PartialEq, Eq)]
struct CanonicalCandidate<K> {
    key: K,
    entity_order: Vec<NodeId>,
}

#[allow(dead_code)]
#[derive(Clone, Debug, PartialEq, Eq)]
struct CanonicalSearchResult<K> {
    candidate: CanonicalCandidate<K>,
    stats: CanonicalSearchStats,
}

#[allow(dead_code)]
/// Minimize a typed leaf key over the adapter's ordered exact partition.
///
/// Automorphism pruning requires the leaf key to be invariant under adapter automorphisms.
/// Prefix pruning requires `prefix_worse` to reject only partitions whose every leaf is greater
/// than the current best key.
fn canonical_search<K, LeafKey, PrefixWorse>(
    adapter: &AutomorphismAdapter,
    algorithm: AutomorphismAlgorithm,
    options: CanonicalSearchOptions,
    leaf_key: &LeafKey,
    prefix_worse: &PrefixWorse,
) -> CanonicalSearchResult<K>
where
    K: Ord,
    LeafKey: Fn(&[NodeId]) -> K,
    PrefixWorse: Fn(&OrderedPartition, &CanonicalCandidate<K>) -> bool,
{
    let initial = OrderedPartition::from_classes(&adapter.classes).refine(adapter.graph());
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
        leaf_key,
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
fn search_partition<K, LeafKey, PrefixWorse>(
    adapter: &AutomorphismAdapter,
    partition: OrderedPartition,
    algorithm: AutomorphismAlgorithm,
    options: CanonicalSearchOptions,
    leaf_key: &LeafKey,
    prefix_worse: &PrefixWorse,
    best: &mut Option<CanonicalCandidate<K>>,
    stats: &mut CanonicalSearchStats,
) where
    K: Ord,
    LeafKey: Fn(&[NodeId]) -> K,
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
        let candidate = CanonicalCandidate {
            key: leaf_key(&entity_order),
            entity_order,
        };
        if best.as_ref().is_none_or(|best| {
            (&candidate.key, &candidate.entity_order) < (&best.key, &best.entity_order)
        }) {
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
        BranchOrder::Node => candidates.sort_unstable(),
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
            leaf_key,
            prefix_worse,
            best,
            stats,
        );
    }
}

#[allow(dead_code)]
fn exhaustive_minimum<K, LeafKey>(
    adapter: &AutomorphismAdapter,
    leaf_key: &LeafKey,
) -> CanonicalCandidate<K>
where
    K: Ord,
    LeafKey: Fn(&[NodeId]) -> K,
{
    fn visit_cells<K, LeafKey>(
        cells: &mut [Vec<NodeId>],
        cell_index: usize,
        order: &mut Vec<NodeId>,
        leaf_key: &LeafKey,
        best: &mut Option<CanonicalCandidate<K>>,
    ) where
        K: Ord,
        LeafKey: Fn(&[NodeId]) -> K,
    {
        fn visit_permutations<K, LeafKey>(
            cells: &mut [Vec<NodeId>],
            cell_index: usize,
            position: usize,
            order: &mut Vec<NodeId>,
            leaf_key: &LeafKey,
            best: &mut Option<CanonicalCandidate<K>>,
        ) where
            K: Ord,
            LeafKey: Fn(&[NodeId]) -> K,
        {
            if position == cells[cell_index].len() {
                let old_len = order.len();
                order.extend_from_slice(&cells[cell_index]);
                visit_cells(cells, cell_index + 1, order, leaf_key, best);
                order.truncate(old_len);
                return;
            }

            for next in position..cells[cell_index].len() {
                cells[cell_index].swap(position, next);
                visit_permutations(cells, cell_index, position + 1, order, leaf_key, best);
                cells[cell_index].swap(position, next);
            }
        }

        if cell_index == cells.len() {
            let candidate = CanonicalCandidate {
                key: leaf_key(order),
                entity_order: order.clone(),
            };
            if best.as_ref().is_none_or(|best| {
                (&candidate.key, &candidate.entity_order) < (&best.key, &best.entity_order)
            }) {
                *best = Some(candidate);
            }
            return;
        }

        visit_permutations(cells, cell_index, 0, order, leaf_key, best);
    }

    let mut cells = adapter.entity_blocks.clone();
    let mut best = None;
    visit_cells(
        &mut cells,
        0,
        &mut Vec::with_capacity(adapter.source_node_count),
        leaf_key,
        &mut best,
    );

    best.expect("every finite partition has an entity ordering")
}

#[allow(dead_code)]
fn initial_classes(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
) -> Result<InitialClasses, Contradiction> {
    let entity_keys = incidence_graph
        .graph()
        .node_ids()
        .map(|node| entity_class_key(molecule, incidence_graph.entity(node)))
        .collect::<Result<Vec<_>, _>>()?;
    let incidence_keys = incidence_graph
        .incidences()
        .map(|(_, incidence)| incidence_key(incidence).map(InitialClassKey::Incidence))
        .collect::<Result<Vec<_>, _>>()?;
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

    Ok(InitialClasses {
        entities: entity_keys.iter().map(|key| classes[key]).collect(),
        incidences: incidence_keys.iter().map(|key| classes[key]).collect(),
    })
}

fn entity_class_key(molecule: &Molecule, entity: Entity) -> Result<InitialClassKey, Contradiction> {
    let (position, value) = match entity {
        Entity::Atom(id) => {
            let attributes = molecule.atom(id).attributes;
            (
                EntityBlockPosition::ATOM,
                product([
                    element_form_key(attributes.element.normalized()?.as_ref()),
                    isotope_mass_form_key(attributes.isotope_mass.normalized()?.as_ref()),
                    num_form_key(attributes.charge.normalized()?.as_ref()),
                    num_form_key(attributes.implicit_hydrogens.normalized()?.as_ref()),
                    num_form_key(attributes.lone_pairs.normalized()?.as_ref()),
                    unpaired_electrons_form_key(
                        attributes.unpaired_electrons.normalized()?.as_ref(),
                    ),
                ]),
            )
        }
        Entity::Bond(id) => {
            let attributes = molecule.bond(id).attributes;
            (
                EntityBlockPosition::BOND,
                positioned_product([
                    (1, num_form_key(attributes.order.normalized()?.as_ref())),
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
    use std::cmp::Ordering;

    use rstest::{fixture, rstest};
    use umol_chem::element::Element;

    use super::*;
    use crate::ir::{
        AromaticSystemForm, AromaticSystemId, AtomConstraintForm, AtomForm, AtomId, BondForm,
        BondId, DativeBondForm, DativeBondId, Entity, IncidenceLevel, MoleculeEntries,
        MulticenterBondForm, MulticenterBondId, NoncovalentBondForm, NoncovalentBondId,
        StereoAtomForm, StereoAtomId, StereoBondForm, StereoBondId, StereoCoset, StereoLigand,
    };

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
        let atom_ligands = vec![
            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        ];
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
                    vec![AtomId(3), AtomId(4)],
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
                    atom_ligands.clone(),
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                ),
                (
                    AtomId(1),
                    atom_ligands.clone(),
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                ),
                (
                    AtomId(2),
                    atom_ligands,
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
    #[case::localized_self_loop(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            bonds: vec![(AtomId(0), AtomId(0), BondForm::from_order(1))],
            ..Default::default()
        }),
        vec![Incidence::BondEndpoint, Incidence::BondEndpoint],
    )]
    #[case::parallel_bonds(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
            ],
            ..Default::default()
        }),
        vec![
            Incidence::BondEndpoint,
            Incidence::BondEndpoint,
            Incidence::BondEndpoint,
            Incidence::BondEndpoint,
        ],
    )]
    #[case::repeated_relation_participant(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::B),
            ],
            dative: vec![(
                vec![AtomId(0), AtomId(0)],
                AtomId(1),
                DativeBondForm::from_order(1),
            )],
            ..Default::default()
        }),
        vec![
            Incidence::DativeDonor,
            Incidence::DativeDonor,
            Incidence::DativeAcceptor,
        ],
    )]
    fn test_automorphism_adapter_new(
        #[case] molecule: Molecule,
        #[case] expected_incidences: Vec<Incidence>,
    ) {
        let incidence_graph = molecule.incidence_graph(IncidenceLevel::Constitution);
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
            source.node_count() + source.edge_count(),
        );
        assert_eq!(adapter.graph().edge_count(), 2 * source.edge_count());
        assert!(adapter.graph().is_simple());

        for node in source.node_ids() {
            let adapter_node = adapter.node_of(SubdivisionNodeSource::Node(node));
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
            let adapter_node = adapter.node_of(SubdivisionNodeSource::Edge(edge));
            assert_eq!(
                adapter.node_source(adapter_node),
                SubdivisionNodeSource::Edge(edge),
            );
            assert_eq!(
                adapter.class(adapter_node),
                AutomorphismClass::Incidence(classes.incidences[edge.index()]),
            );
            for incidence_edge in adapter.incidence_edges_of(edge) {
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
    fn test_ordered_partition_from_classes(
        #[case] classes: Vec<u32>,
        #[case] expected: OrderedPartition,
    ) {
        assert_eq!(OrderedPartition::from_classes(&classes), expected);
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
            OrderedPartition::from_classes(&[0_u32; 4]).refine(&graph),
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
    #[case::parallel_bonds(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 2],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
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
        let classes = initial_classes(&molecule, &incidence_graph).unwrap();
        let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
        let source = incidence_graph.graph();
        let leaf_key = |order: &[NodeId]| {
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
        };
        let no_prefix = |_: &OrderedPartition, _: &CanonicalCandidate<_>| false;
        let expected = exhaustive_minimum(&adapter, &leaf_key);
        let unpruned = canonical_search(
            &adapter,
            AutomorphismAlgorithm::Nauty,
            CanonicalSearchOptions {
                automorphism_pruning: false,
                prefix_pruning: false,
                branch_order: BranchOrder::Node,
            },
            &leaf_key,
            &no_prefix,
        );
        let reversed = canonical_search(
            &adapter,
            AutomorphismAlgorithm::Nauty,
            CanonicalSearchOptions {
                automorphism_pruning: false,
                prefix_pruning: false,
                branch_order: BranchOrder::ReverseNode,
            },
            &leaf_key,
            &no_prefix,
        );
        let pruned = canonical_search(
            &adapter,
            AutomorphismAlgorithm::Nauty,
            CanonicalSearchOptions {
                automorphism_pruning: true,
                prefix_pruning: false,
                branch_order: BranchOrder::BackendCanonical,
            },
            &leaf_key,
            &no_prefix,
        );

        assert_eq!(unpruned.candidate, expected);
        assert_eq!(reversed.candidate, expected);
        assert_eq!(pruned.candidate, expected);
        assert!(pruned.stats.visited_leaves <= unpruned.stats.visited_leaves);
    }

    #[rstest]
    fn test_canonical_search_prefix() {
        let source = Graph::new(4, &[]);
        let adapter = AutomorphismAdapter {
            subdivision: source.subdivide_edges(),
            classes: vec![AutomorphismClass::Entity(0); 4],
            source_node_count: 4,
            entity_blocks: vec![(0..4).map(NodeId).collect()],
        };
        let leaf_key = |order: &[NodeId]| order.to_vec();
        let prefix_worse = |partition: &OrderedPartition,
                            best: &CanonicalCandidate<Vec<NodeId>>| {
            let prefix = partition.fixed_entity_prefix(4);
            prefix.as_slice() > &best.key[..prefix.len()]
        };
        let expected = exhaustive_minimum(&adapter, &leaf_key);
        let actual = canonical_search(
            &adapter,
            AutomorphismAlgorithm::Nauty,
            CanonicalSearchOptions {
                automorphism_pruning: false,
                prefix_pruning: true,
                branch_order: BranchOrder::Node,
            },
            &leaf_key,
            &prefix_worse,
        );

        assert_eq!(actual.candidate, expected);
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
            let adapter = AutomorphismAdapter {
                subdivision: source.subdivide_edges(),
                classes: (0..node_count)
                    .map(|_| AutomorphismClass::Entity(0))
                    .chain((0..edges.len()).map(|_| AutomorphismClass::Incidence(1)))
                    .collect(),
                source_node_count: node_count,
                entity_blocks: vec![(0..node_count as u32).map(NodeId).collect()],
            };
            let leaf_key = |order: &[NodeId]| {
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
                mapped_edges
            };
            let no_prefix = |_: &OrderedPartition, _: &CanonicalCandidate<_>| false;
            let expected = exhaustive_minimum(&adapter, &leaf_key);

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
                        AutomorphismAlgorithm::Nauty,
                        options,
                        &leaf_key,
                        &no_prefix,
                    )
                    .candidate,
                    expected,
                    "edge mask {edge_mask:#08b}",
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
