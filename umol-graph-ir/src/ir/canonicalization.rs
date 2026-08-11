//! Aggregate canonicalization inputs and failures.

use std::cmp::Ordering;

use thiserror::Error;
use umol_graph_core::AutomorphismAlgorithm;

use super::error::Contradiction;
use super::molecule::MoleculeIntegrityError;
use super::reaction_span::ReactionSpanIntegrityError;
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

    use rstest::rstest;

    use super::*;
    use crate::ir::{AtomId, Entity};

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
