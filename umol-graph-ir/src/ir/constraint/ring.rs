//! Ring-membership scope and single-entry fact, shared by atom/bond/dative constraints.
//!
//! Atom and localized-bond values derived from molecule topology use the Relevant ring projection
//! through size 22. `RingScope` selects a count; it does not select ring-set semantics.

use super::super::error::{Contradiction, NoJoin};
use super::super::num::NumForm;
use super::super::traits::{Canonicalize, Lattice};

/// `All` = total ring count; `Size(s)` = count of size-`s` rings. `All` sorts first.
///
/// For derived atom and localized-bond values, the count is taken within the fixed Relevant ring
/// projection through size 22.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum RingScope {
    All,
    Size(u8),
}

/// A ring count under the semantics of its containing entity constraint.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct RingMembershipForm {
    pub scope: RingScope,
    pub count: NumForm,
}

impl RingMembershipForm {
    pub fn new(scope: RingScope, count: impl Into<NumForm>) -> Self {
        Self {
            scope,
            count: count.into(),
        }
    }
}

impl Canonicalize for RingMembershipForm {
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(Self::new(self.scope, self.count.canonicalize()?))
    }
}

/// Meet-semilattice keyed by `scope`: same scope delegates to the `count`
/// value-lattice, different scopes lie in different fibers (`meet` → `None`,
/// `join` → `Err(NoJoin)`).
impl Lattice for RingMembershipForm {
    fn is_undetermined(&self) -> bool {
        self.count.is_undetermined()
    }

    fn is_ground(&self) -> bool {
        self.count.is_ground()
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        if self.scope != other.scope {
            return None;
        }
        self.count
            .meet(&other.count)
            .map(|count| Self::new(self.scope, count))
    }

    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        if self.scope != other.scope {
            return Err(NoJoin);
        }
        Ok(Self::new(self.scope, self.count.join(&other.count)?))
    }

    fn matches(&self, target: &Self) -> bool {
        self.scope == target.scope && self.count.matches(&target.count)
    }

    fn is_compatible(&self, other: &Self) -> bool {
        self.scope == other.scope && self.count.is_compatible(&other.count)
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::folds_count(
        RingMembershipForm::new(RingScope::Size(6), NumForm::lit_set([2])),
        Ok(RingMembershipForm::new(RingScope::Size(6), NumForm::Lit(2)))
    )]
    #[case::empty_count_contradiction(
        RingMembershipForm::new(RingScope::All, NumForm::lit_set(Vec::<i64>::new())),
        Err(Contradiction)
    )]
    fn test_ring_membership_form_canonicalize(
        #[case] input: RingMembershipForm,
        #[case] expected: Result<RingMembershipForm, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rstest]
    #[case::undetermined(
        RingMembershipForm::new(RingScope::All, NumForm::Undetermined),
        true,
        false
    )]
    #[case::ground(RingMembershipForm::new(RingScope::All, NumForm::Lit(3)), false, true)]
    fn test_ring_membership_form_lattice_position(
        #[case] membership: RingMembershipForm,
        #[case] is_undetermined: bool,
        #[case] is_ground: bool,
    ) {
        assert_eq!(membership.is_undetermined(), is_undetermined);
        assert_eq!(membership.is_ground(), is_ground);
    }

    #[rstest]
    #[case::same_scope_narrows(
        RingMembershipForm::new(RingScope::All, NumForm::Undetermined),
        RingMembershipForm::new(RingScope::All, NumForm::Lit(3)),
        Some(RingMembershipForm::new(RingScope::All, NumForm::Lit(3)))
    )]
    #[case::same_scope_incompatible(
        RingMembershipForm::new(RingScope::All, NumForm::Lit(3)),
        RingMembershipForm::new(RingScope::All, NumForm::Lit(4)),
        None
    )]
    #[case::different_scope(
        RingMembershipForm::new(RingScope::All, NumForm::Lit(3)),
        RingMembershipForm::new(RingScope::Size(6), NumForm::Lit(3)),
        None
    )]
    fn test_ring_membership_form_meet(
        #[case] a: RingMembershipForm,
        #[case] b: RingMembershipForm,
        #[case] expected: Option<RingMembershipForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::same_scope(
        RingMembershipForm::new(RingScope::All, NumForm::Lit(3)),
        RingMembershipForm::new(RingScope::All, NumForm::Lit(3)),
        Ok(RingMembershipForm::new(RingScope::All, NumForm::Lit(3)))
    )]
    #[case::different_scope(
        RingMembershipForm::new(RingScope::All, NumForm::Lit(3)),
        RingMembershipForm::new(RingScope::Size(6), NumForm::Lit(3)),
        Err(NoJoin)
    )]
    fn test_ring_membership_form_join(
        #[case] a: RingMembershipForm,
        #[case] b: RingMembershipForm,
        #[case] expected: Result<RingMembershipForm, NoJoin>,
    ) {
        assert_eq!(a.join(&b), expected);
    }

    #[rstest]
    #[case::same_scope_matches(
        RingMembershipForm::new(RingScope::All, NumForm::Undetermined),
        RingMembershipForm::new(RingScope::All, NumForm::Lit(3)),
        true
    )]
    #[case::same_scope_no_match(
        RingMembershipForm::new(RingScope::All, NumForm::Lit(3)),
        RingMembershipForm::new(RingScope::All, NumForm::Lit(4)),
        false
    )]
    #[case::different_scope(
        RingMembershipForm::new(RingScope::All, NumForm::Undetermined),
        RingMembershipForm::new(RingScope::Size(6), NumForm::Lit(3)),
        false
    )]
    fn test_ring_membership_form_matches(
        #[case] pattern: RingMembershipForm,
        #[case] target: RingMembershipForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    #[case::same_scope_compatible(
        RingMembershipForm::new(RingScope::All, NumForm::Undetermined),
        RingMembershipForm::new(RingScope::All, NumForm::Lit(3)),
        true
    )]
    #[case::same_scope_incompatible(
        RingMembershipForm::new(RingScope::All, NumForm::Lit(3)),
        RingMembershipForm::new(RingScope::All, NumForm::Lit(4)),
        false
    )]
    #[case::different_scope(
        RingMembershipForm::new(RingScope::All, NumForm::Lit(3)),
        RingMembershipForm::new(RingScope::Size(6), NumForm::Lit(3)),
        false
    )]
    fn test_ring_membership_form_is_compatible(
        #[case] a: RingMembershipForm,
        #[case] b: RingMembershipForm,
        #[case] expected: bool,
    ) {
        assert_eq!(a.is_compatible(&b), expected);
    }
}
