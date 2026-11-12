//! Entity and constraint spans.
//!
//! `EntitySpan` and `ConstraintSpan` represent left and right states in a
//! reaction span.

use super::constraint::{Constraint, ConstraintFrameActionDomain, ConstraintFrameActions};
use super::entity::Entity;
use super::error::Contradiction;
use super::frame::OverlaysFrameAction;
use super::traits::{AsColor, FrameTransport, Normalize};

/// One entity's span across a reaction — its slice of the superimposed `L`∪`K`∪`R`. A *state*, not
/// an operation (unlike `Edit` / `Delta`). `lhs()` / `rhs()` read the side values.
///
/// # Semantic properties
///
/// The derived order follows the declaration order of the variants. `Added` is declared last, after
/// `Removed`, so that the order places every entity present on the lhs before the added ones, as
/// the lhs-anchored form of a reaction span does.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EntitySpan<T> {
    /// In the interface `K` — present and identical on both sides.
    Unchanged(T),
    /// In `L` only — deleted.
    Removed(T),
    /// In the interface `K` — present on both sides but relabeled (a dynamic entity).
    Modified { lhs: T, rhs: T },
    /// In `R` only — created.
    Added(T),
}

impl<T: Default> Default for EntitySpan<T> {
    fn default() -> Self {
        Self::Unchanged(T::default())
    }
}

impl<T> EntitySpan<T> {
    /// The lhs (`L`) value, or `None` if the entity is created.
    pub fn lhs(&self) -> Option<&T> {
        match self {
            Self::Unchanged(value) | Self::Removed(value) | Self::Modified { lhs: value, .. } => {
                Some(value)
            }
            Self::Added(_) => None,
        }
    }

    /// The rhs (`R`) value, or `None` if the entity is deleted.
    pub fn rhs(&self) -> Option<&T> {
        match self {
            Self::Unchanged(value) | Self::Added(value) | Self::Modified { rhs: value, .. } => {
                Some(value)
            }
            Self::Removed(_) => None,
        }
    }

    /// Carry every present side through `f`, declining if any side declines.
    ///
    /// A `Modified` span applies the same `f` to both sides, which is what makes one selected frame
    /// action reach both: the span holds two values against a single participant list.
    pub fn try_map<U>(self, mut f: impl FnMut(T) -> Option<U>) -> Option<EntitySpan<U>> {
        Some(match self {
            Self::Unchanged(value) => EntitySpan::Unchanged(f(value)?),
            Self::Modified { lhs, rhs } => EntitySpan::Modified {
                lhs: f(lhs)?,
                rhs: f(rhs)?,
            },
            Self::Added(value) => EntitySpan::Added(f(value)?),
            Self::Removed(value) => EntitySpan::Removed(f(value)?),
        })
    }
}

impl<T: FrameTransport> FrameTransport for EntitySpan<T> {
    type Action = T::Action;

    fn reframe_by(self, action: &Self::Action) -> Option<Self> {
        self.try_map(|value| value.reframe_by(action))
    }
}

impl<T: Normalize> Normalize for EntitySpan<T> {
    fn normalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::Unchanged(value) => Self::Unchanged(value.normalize()?),
            Self::Modified { lhs, rhs } => {
                let lhs = lhs.normalize()?;
                let rhs = rhs.normalize()?;
                if lhs == rhs {
                    Self::Unchanged(lhs)
                } else {
                    Self::Modified { lhs, rhs }
                }
            }
            Self::Added(value) => Self::Added(value.normalize()?),
            Self::Removed(value) => Self::Removed(value.normalize()?),
        })
    }
}

impl<T: Normalize> EntitySpan<T> {
    /// Superimpose an entity's optional lhs and rhs values into a span — the per-entity kernel of
    /// `ReactionSpan::superimpose`: present-both maps to `Unchanged(lhs)` when the values are
    /// semantically equivalent and to `Modified` otherwise, lhs-only to `Removed`, rhs-only to
    /// `Added`, neither to `None`.
    pub fn superimpose(lhs: Option<T>, rhs: Option<T>) -> Option<Self> {
        match (lhs, rhs) {
            (Some(lhs), Some(rhs)) if lhs.normalized_eq(&rhs) => Some(Self::Unchanged(lhs)),
            (Some(lhs), Some(rhs)) => Some(Self::Modified { lhs, rhs }),
            (Some(lhs), None) => Some(Self::Removed(lhs)),
            (None, Some(rhs)) => Some(Self::Added(rhs)),
            (None, None) => None,
        }
    }
}

/// The color of an entity span: the span of the colors of its sides. A `Modified` span whose two
/// sides have equal colors gives `Unchanged`.
impl<T: AsColor> AsColor for EntitySpan<T> {
    type Color<'a>
        = EntitySpan<T::Color<'a>>
    where
        Self: 'a;

    fn as_color(&self) -> Self::Color<'_> {
        match self {
            Self::Unchanged(value) => EntitySpan::Unchanged(value.as_color()),
            Self::Removed(value) => EntitySpan::Removed(value.as_color()),
            Self::Modified { lhs, rhs } => {
                let lhs = lhs.as_color();
                let rhs = rhs.as_color();
                if lhs == rhs {
                    EntitySpan::Unchanged(lhs)
                } else {
                    EntitySpan::Modified { lhs, rhs }
                }
            }
            Self::Added(value) => EntitySpan::Added(value.as_color()),
        }
    }
}

/// A molecule-level constraint's span across a reaction — its slice of the superimposed `L`∪`K`∪`R`.
/// A *state*, not an operation (unlike `ConstraintDelta`). `lhs()` / `rhs()` read the side values.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ConstraintSpan {
    /// In the interface `K` — present and identical on both sides.
    Unchanged(Constraint),
    /// In `R` only — created.
    Added(Constraint),
    /// In `L` only — deleted.
    Removed(Constraint),
}

impl ConstraintSpan {
    /// The lhs (`L`) value, or `None` if the constraint is created.
    pub fn lhs(&self) -> Option<&Constraint> {
        match self {
            Self::Unchanged(value) | Self::Removed(value) => Some(value),
            Self::Added(_) => None,
        }
    }

    /// The rhs (`R`) value, or `None` if the constraint is deleted.
    pub fn rhs(&self) -> Option<&Constraint> {
        match self {
            Self::Unchanged(value) | Self::Added(value) => Some(value),
            Self::Removed(_) => None,
        }
    }

    pub(crate) fn collect_frame_action_domain(&self, domain: &mut ConstraintFrameActionDomain) {
        match self {
            Self::Unchanged(constraint) | Self::Added(constraint) | Self::Removed(constraint) => {
                constraint.collect_frame_action_domain(domain);
            }
        }
    }

    pub(crate) fn reframe_by_actions(
        self,
        actions: &impl ConstraintFrameActions,
    ) -> Result<Self, Entity> {
        Ok(match self {
            Self::Unchanged(constraint) => Self::Unchanged(constraint.reframe_by_actions(actions)?),
            Self::Added(constraint) => Self::Added(constraint.reframe_by_actions(actions)?),
            Self::Removed(constraint) => Self::Removed(constraint.reframe_by_actions(actions)?),
        })
    }
}

impl Normalize for ConstraintSpan {
    fn normalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::Unchanged(value) => Self::Unchanged(value.normalize()?),
            Self::Added(value) => Self::Added(value.normalize()?),
            Self::Removed(value) => Self::Removed(value.normalize()?),
        })
    }
}

impl FrameTransport for ConstraintSpan {
    type Action = OverlaysFrameAction;

    fn reframe_by(self, actions: &Self::Action) -> Option<Self> {
        self.reframe_by_actions(actions).ok()
    }
}

#[cfg(test)]
mod tests {
    use std::cmp::Ordering;

    use rstest::*;
    use umol_chem::element::Element;
    use umol_perm::{DynPermutation, Permutation};

    use super::super::aromatic::AromaticSystemForm;
    use super::super::atom::AtomForm;
    use super::super::constraint::{AtomConstraintForm, RelationalConstraint};
    use super::super::frame::{
        AromaticSystemsFrameAction, DativeBondsFrameAction, MulticenterBondsFrameAction,
        NoncovalentBondsFrameAction, StereoAtomsFrameAction, StereoBondsFrameAction,
    };
    use super::super::id::{AtomId, NoncovalentBondId};
    use super::super::num::NumForm;
    use super::*;

    #[fixture]
    fn overlays_frame_action() -> OverlaysFrameAction {
        OverlaysFrameAction::new(
            DativeBondsFrameAction::from_vec(vec![]).expect("actions are admissible"),
            AromaticSystemsFrameAction::from_vec(vec![
                DynPermutation::try_from(vec![2, 0, 1]).expect("action is a permutation")
            ])
            .expect("action is admissible"),
            MulticenterBondsFrameAction::from_vec(vec![]).expect("actions are admissible"),
            NoncovalentBondsFrameAction::from_vec(vec![
                DynPermutation::try_from(vec![1, 0]).expect("action is a permutation")
            ])
            .expect("action is admissible"),
            StereoAtomsFrameAction::from_vec(vec![Permutation::from_image(&[1, 0, 2, 3])])
                .expect("action is admissible"),
            StereoBondsFrameAction::from_vec(vec![Permutation::from_image(&[1, 0, 2, 3])])
                .expect("action is admissible"),
        )
    }

    #[rstest]
    #[case::unchanged_removed(EntitySpan::Unchanged(1), EntitySpan::Removed(0), Ordering::Less)]
    #[case::removed_modified(
        EntitySpan::Removed(1),
        EntitySpan::Modified { lhs: 0, rhs: 0 },
        Ordering::Less
    )]
    #[case::modified_added(
        EntitySpan::Modified { lhs: 1, rhs: 1 },
        EntitySpan::Added(0),
        Ordering::Less
    )]
    fn test_entity_span_cmp(
        #[case] lhs: EntitySpan<i32>,
        #[case] rhs: EntitySpan<i32>,
        #[case] expected: Ordering,
    ) {
        assert_eq!(lhs.cmp(&rhs), expected);
    }

    #[rstest]
    fn test_entity_span_default() {
        assert_eq!(
            EntitySpan::<AtomForm>::default(),
            EntitySpan::Unchanged(AtomForm::default()),
        );
    }

    #[rstest]
    #[case::unchanged(EntitySpan::Unchanged(5), Some(&5))]
    #[case::modified(EntitySpan::Modified { lhs: 1, rhs: 2 }, Some(&1))]
    #[case::removed(EntitySpan::Removed(7), Some(&7))]
    #[case::added(EntitySpan::Added(9), None)]
    fn test_entity_span_lhs(#[case] state: EntitySpan<i32>, #[case] expected: Option<&i32>) {
        assert_eq!(state.lhs(), expected);
    }

    #[rstest]
    #[case::unchanged(EntitySpan::Unchanged(5), Some(&5))]
    #[case::modified(EntitySpan::Modified { lhs: 1, rhs: 2 }, Some(&2))]
    #[case::added(EntitySpan::Added(9), Some(&9))]
    #[case::removed(EntitySpan::Removed(7), None)]
    fn test_entity_span_rhs(#[case] state: EntitySpan<i32>, #[case] expected: Option<&i32>) {
        assert_eq!(state.rhs(), expected);
    }

    #[rstest]
    #[case::unchanged(EntitySpan::Unchanged(5), Some(EntitySpan::Unchanged(10)))]
    #[case::modified(
        EntitySpan::Modified { lhs: 1, rhs: 2 },
        Some(EntitySpan::Modified { lhs: 2, rhs: 4 }),
    )]
    #[case::added(EntitySpan::Added(9), Some(EntitySpan::Added(18)))]
    #[case::removed(EntitySpan::Removed(7), Some(EntitySpan::Removed(14)))]
    #[case::modified_lhs_declines(EntitySpan::Modified { lhs: 0, rhs: 2 }, None)]
    #[case::modified_rhs_declines(EntitySpan::Modified { lhs: 1, rhs: 0 }, None)]
    #[case::unchanged_declines(EntitySpan::Unchanged(0), None)]
    fn test_entity_span_try_map(
        #[case] span: EntitySpan<i32>,
        #[case] expected: Option<EntitySpan<i32>>,
    ) {
        assert_eq!(span.try_map(|v| (v != 0).then_some(v * 2)), expected);
    }

    #[rstest]
    fn test_entity_span_reframe_by() {
        let action =
            DynPermutation::try_from([2, 0, 1].as_slice()).expect("test image is a permutation");
        let input = EntitySpan::Modified {
            lhs: AromaticSystemForm::from_electrons(vec![10, 20, 30]),
            rhs: AromaticSystemForm::from_electrons(vec![11, 21, 31]),
        };
        let expected = EntitySpan::Modified {
            lhs: AromaticSystemForm::from_electrons(vec![30, 10, 20]),
            rhs: AromaticSystemForm::from_electrons(vec![31, 11, 21]),
        };

        assert_eq!(input.reframe_by(&action), Some(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::unchanged(
        EntitySpan::Unchanged(NumForm::lit_set([2])),
        EntitySpan::Unchanged(NumForm::Lit(2)),
    )]
    #[case::modified_collapses(
        EntitySpan::Modified { lhs: NumForm::lit_set([2]), rhs: NumForm::Lit(2) },
        EntitySpan::Unchanged(NumForm::Lit(2)),
    )]
    #[case::modified_remains(
        EntitySpan::Modified { lhs: NumForm::lit_set([2]), rhs: NumForm::Lit(3) },
        EntitySpan::Modified { lhs: NumForm::Lit(2), rhs: NumForm::Lit(3) },
    )]
    fn test_entity_span_normalize(
        #[case] input: EntitySpan<NumForm>,
        #[case] expected: EntitySpan<NumForm>,
    ) {
        assert_eq!(input.normalize(), Ok(expected));
    }

    #[rstest]
    #[case::singleton_set(NumForm::Lit(1), NumForm::lit_set([1]))]
    fn test_entity_span_superimpose_canonical(#[case] lhs: NumForm, #[case] rhs: NumForm) {
        // Equivalent sides that are structurally distinct superimpose to `Unchanged`,
        // not `Modified`.
        let lhs = AtomForm::from_element(Element::C).with_charge(lhs);
        let rhs = AtomForm::from_element(Element::C).with_charge(rhs);
        assert_eq!(
            EntitySpan::superimpose(Some(lhs.clone()), Some(rhs)),
            Some(EntitySpan::Unchanged(lhs))
        );
    }

    #[rstest]
    #[case::unchanged(
        EntitySpan::Unchanged(AromaticSystemForm::default().with_charge(-1)),
        EntitySpan::Unchanged(AromaticSystemForm::default().with_charge(-1)),
        Ordering::Equal
    )]
    #[case::removed(
        EntitySpan::Unchanged(AromaticSystemForm::default().with_charge(-1)),
        EntitySpan::Removed(AromaticSystemForm::default().with_charge(-1)),
        Ordering::Less
    )]
    #[case::modified(
        EntitySpan::Removed(AromaticSystemForm::default().with_charge(-1)),
        EntitySpan::Modified {
            lhs: AromaticSystemForm::default().with_charge(-1),
            rhs: AromaticSystemForm::default(),
        },
        Ordering::Less
    )]
    #[case::modified_equal_colors(
        EntitySpan::Modified {
            lhs: AromaticSystemForm::from_electrons(vec![2, 1, 1]).with_charge(-1),
            rhs: AromaticSystemForm::from_electrons(vec![1, 1, 2]).with_charge(-1),
        },
        EntitySpan::Unchanged(AromaticSystemForm::default().with_charge(-1)),
        Ordering::Equal
    )]
    #[case::added(
        EntitySpan::Modified {
            lhs: AromaticSystemForm::default().with_charge(-1),
            rhs: AromaticSystemForm::default(),
        },
        EntitySpan::Added(AromaticSystemForm::default().with_charge(-1)),
        Ordering::Less
    )]
    fn test_entity_span_as_color(
        #[case] lhs: EntitySpan<AromaticSystemForm>,
        #[case] rhs: EntitySpan<AromaticSystemForm>,
        #[case] expected: Ordering,
    ) {
        assert_eq!(lhs.as_color().cmp(&rhs.as_color()), expected);
    }

    #[rstest]
    fn test_constraint_span_normalize() {
        let input = ConstraintSpan::Added(Constraint::Atom(
            AtomId(0),
            AtomConstraintForm::Valence(NumForm::lit_set([4])),
        ));
        let expected = ConstraintSpan::Added(Constraint::Atom(
            AtomId(0),
            AtomConstraintForm::Valence(NumForm::Lit(4)),
        ));

        assert_eq!(input.normalize(), Ok(expected));
    }

    #[rstest]
    fn test_constraint_span_reframe_by(overlays_frame_action: OverlaysFrameAction) {
        let input = ConstraintSpan::Removed(Constraint::Relational(
            RelationalConstraint::NoncovalentBondEndsSatisfy {
                bond: NoncovalentBondId(0),
                predicates: [
                    Box::new(AtomConstraintForm::Valence(NumForm::Lit(4))),
                    Box::new(AtomConstraintForm::Degree(NumForm::Lit(2))),
                ],
            },
        ));
        let expected = ConstraintSpan::Removed(Constraint::Relational(
            RelationalConstraint::NoncovalentBondEndsSatisfy {
                bond: NoncovalentBondId(0),
                predicates: [
                    Box::new(AtomConstraintForm::Degree(NumForm::Lit(2))),
                    Box::new(AtomConstraintForm::Valence(NumForm::Lit(4))),
                ],
            },
        ));

        assert_eq!(input.reframe_by(&overlays_frame_action), Some(expected));
    }
}
