//! Dative bond AST.

use super::constraint::DativeBondConstraint;

/// Dative bond: two-atom bond with a fixed two-electron donation from donor
/// to acceptor. Directionality lives on the relation (donor/acceptor
/// endpoints), not on the AST. No inherent fields; ring-membership and
/// ring-size constraints may be attached.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DativeBondAst {
    pub constraints: Vec<DativeBondConstraint>,
}

impl DativeBondAst {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn is_ground(&self) -> bool {
        true
    }

    /// A dative bond has no identity-bearing fields, so every dative pattern
    /// matches every dative target. Constraints, when added, filter matches
    /// topologically, not per-slot.
    pub fn matches(&self, _target: &DativeBondAst) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::super::value::ValueAst;
    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::default_(DativeBondAst::default(), true)]
    #[case::with_ground_constraint(DativeBondAst { constraints: vec![DativeBondConstraint::RingSize(ValueAst::Lit(6))] }, true)]
    #[case::with_undetermined_constraint(DativeBondAst { constraints: vec![DativeBondConstraint::RingSize(ValueAst::Undetermined)] }, true)]
    fn test_dative_bond_ast_is_ground(#[case] ast: DativeBondAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rstest]
    #[case::default_default(DativeBondAst::default(), DativeBondAst::default(), true)]
    fn test_dative_bond_ast_matches(
        #[case] pattern: DativeBondAst,
        #[case] target: DativeBondAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }
}
