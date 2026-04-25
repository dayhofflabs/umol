//! Dative bond AST.

use super::constraint::DativeBondConstraints;

/// Direction of a dative bond relative to its `FixedRelationSet` participants
/// array (sorted ascending by `NodeId`). `Forward` means the donor is
/// `participants[0]`; `Reverse` means the donor is `participants[1]`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum DativeBondDirection {
    #[default]
    Forward,
    Reverse,
}

impl DativeBondDirection {
    pub fn flip(self) -> Self {
        match self {
            Self::Forward => Self::Reverse,
            Self::Reverse => Self::Forward,
        }
    }
}

/// Dative bond: two-atom bond with a fixed two-electron donation from donor
/// to acceptor. Direction is carried on `direction`, relative to the sorted
/// `FixedRelationSet` participants array; endpoint order in participants is
/// canonical (ascending `NodeId`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DativeBondAst {
    pub direction: DativeBondDirection,
    pub constraints: DativeBondConstraints,
}

impl DativeBondAst {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn is_ground(&self) -> bool {
        true
    }

    /// A dative bond has no identity-bearing fields beyond direction, so every
    /// dative pattern matches every dative target with matching direction.
    /// Constraints, when added, filter matches topologically, not per-slot.
    pub fn matches(&self, target: &DativeBondAst) -> bool {
        self.direction == target.direction
    }

    /// Simplify every constraint's inner value in place.
    pub fn simplify_values(&mut self) {
        self.constraints.simplify_each();
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::super::value::ValueAst;
    use super::*;
    use crate::ast::constraint::DativeBondConstraint;

    #[rustfmt::skip]
    #[rstest]
    #[case::default_(DativeBondAst::default(), true)]
    #[case::with_ground_constraint(DativeBondAst { direction: DativeBondDirection::Forward,
        constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingSize(ValueAst::Lit(6))]) }, true)]
    #[case::with_undetermined_constraint(DativeBondAst { direction: DativeBondDirection::Forward,
        constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingSize(ValueAst::Undetermined)]) }, true)]
    #[case::reverse(DativeBondAst { direction: DativeBondDirection::Reverse, constraints: DativeBondConstraints::new() }, true)]
    fn test_dative_bond_ast_is_ground(#[case] ast: DativeBondAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rstest]
    #[case::same_forward(DativeBondAst::default(), DativeBondAst::default(), true)]
    #[case::both_reverse(
        DativeBondAst { direction: DativeBondDirection::Reverse, constraints: DativeBondConstraints::new() },
        DativeBondAst { direction: DativeBondDirection::Reverse, constraints: DativeBondConstraints::new() },
        true
    )]
    #[case::mismatch(
        DativeBondAst::default(),
        DativeBondAst { direction: DativeBondDirection::Reverse, constraints: DativeBondConstraints::new() },
        false
    )]
    fn test_dative_bond_ast_matches(
        #[case] pattern: DativeBondAst,
        #[case] target: DativeBondAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    #[case(DativeBondDirection::Forward, DativeBondDirection::Reverse)]
    #[case(DativeBondDirection::Reverse, DativeBondDirection::Forward)]
    fn test_dative_direction_flip(
        #[case] input: DativeBondDirection,
        #[case] expected: DativeBondDirection,
    ) {
        assert_eq!(input.flip(), expected);
    }
}
