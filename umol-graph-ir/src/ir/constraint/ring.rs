//! Ring-membership scope and single-entry fact, shared by atom/bond/dative constraints.
//!
//! Atom and localized-bond values derived from molecule topology use the Relevant ring projection
//! through size 22. `RingScope` selects a count; it does not select ring-set semantics.

use super::super::error::{Contradiction, NoJoin};
use super::super::traits::{Canonicalize, Lattice};
use super::super::value::NumForm;

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
pub struct RingMembershipAst {
    pub scope: RingScope,
    pub count: NumForm,
}

impl RingMembershipAst {
    pub fn new(scope: RingScope, count: impl Into<NumForm>) -> Self {
        Self {
            scope,
            count: count.into(),
        }
    }
}

impl Canonicalize for RingMembershipAst {
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(Self::new(self.scope, self.count.canonicalize()?))
    }
}

/// Meet-semilattice keyed by `scope`: same scope delegates to the `count`
/// value-lattice, different scopes lie in different fibers (`meet` → `None`,
/// `join` → `Err(NoJoin)`).
impl Lattice for RingMembershipAst {
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
        RingMembershipAst::new(RingScope::Size(6), NumForm::lit_set([2])),
        Ok(RingMembershipAst::new(RingScope::Size(6), NumForm::Lit(2)))
    )]
    #[case::empty_count_contradiction(
        RingMembershipAst::new(RingScope::All, NumForm::lit_set(Vec::<i64>::new())),
        Err(Contradiction)
    )]
    fn test_ring_membership_ast_canonicalize(
        #[case] input: RingMembershipAst,
        #[case] expected: Result<RingMembershipAst, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rstest]
    #[case::undetermined(
        RingMembershipAst::new(RingScope::All, NumForm::Undetermined),
        true,
        false
    )]
    #[case::ground(RingMembershipAst::new(RingScope::All, NumForm::Lit(3)), false, true)]
    fn test_ring_membership_ast_lattice_position(
        #[case] membership: RingMembershipAst,
        #[case] is_undetermined: bool,
        #[case] is_ground: bool,
    ) {
        assert_eq!(membership.is_undetermined(), is_undetermined);
        assert_eq!(membership.is_ground(), is_ground);
    }

    #[rstest]
    #[case::same_scope_narrows(
        RingMembershipAst::new(RingScope::All, NumForm::Undetermined),
        RingMembershipAst::new(RingScope::All, NumForm::Lit(3)),
        Some(RingMembershipAst::new(RingScope::All, NumForm::Lit(3)))
    )]
    #[case::same_scope_incompatible(
        RingMembershipAst::new(RingScope::All, NumForm::Lit(3)),
        RingMembershipAst::new(RingScope::All, NumForm::Lit(4)),
        None
    )]
    #[case::different_scope(
        RingMembershipAst::new(RingScope::All, NumForm::Lit(3)),
        RingMembershipAst::new(RingScope::Size(6), NumForm::Lit(3)),
        None
    )]
    fn test_ring_membership_ast_meet(
        #[case] a: RingMembershipAst,
        #[case] b: RingMembershipAst,
        #[case] expected: Option<RingMembershipAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::same_scope(
        RingMembershipAst::new(RingScope::All, NumForm::Lit(3)),
        RingMembershipAst::new(RingScope::All, NumForm::Lit(3)),
        Ok(RingMembershipAst::new(RingScope::All, NumForm::Lit(3)))
    )]
    #[case::different_scope(
        RingMembershipAst::new(RingScope::All, NumForm::Lit(3)),
        RingMembershipAst::new(RingScope::Size(6), NumForm::Lit(3)),
        Err(NoJoin)
    )]
    fn test_ring_membership_ast_join(
        #[case] a: RingMembershipAst,
        #[case] b: RingMembershipAst,
        #[case] expected: Result<RingMembershipAst, NoJoin>,
    ) {
        assert_eq!(a.join(&b), expected);
    }

    #[rstest]
    #[case::same_scope_matches(
        RingMembershipAst::new(RingScope::All, NumForm::Undetermined),
        RingMembershipAst::new(RingScope::All, NumForm::Lit(3)),
        true
    )]
    #[case::same_scope_no_match(
        RingMembershipAst::new(RingScope::All, NumForm::Lit(3)),
        RingMembershipAst::new(RingScope::All, NumForm::Lit(4)),
        false
    )]
    #[case::different_scope(
        RingMembershipAst::new(RingScope::All, NumForm::Undetermined),
        RingMembershipAst::new(RingScope::Size(6), NumForm::Lit(3)),
        false
    )]
    fn test_ring_membership_ast_matches(
        #[case] pattern: RingMembershipAst,
        #[case] target: RingMembershipAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    #[case::same_scope_compatible(
        RingMembershipAst::new(RingScope::All, NumForm::Undetermined),
        RingMembershipAst::new(RingScope::All, NumForm::Lit(3)),
        true
    )]
    #[case::same_scope_incompatible(
        RingMembershipAst::new(RingScope::All, NumForm::Lit(3)),
        RingMembershipAst::new(RingScope::All, NumForm::Lit(4)),
        false
    )]
    #[case::different_scope(
        RingMembershipAst::new(RingScope::All, NumForm::Lit(3)),
        RingMembershipAst::new(RingScope::Size(6), NumForm::Lit(3)),
        false
    )]
    fn test_ring_membership_ast_is_compatible(
        #[case] a: RingMembershipAst,
        #[case] b: RingMembershipAst,
        #[case] expected: bool,
    ) {
        assert_eq!(a.is_compatible(&b), expected);
    }
}
