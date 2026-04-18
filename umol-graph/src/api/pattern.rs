//! Substructure-query pattern backed by [`MoleculeAst`].
//!
//! [`MoleculeAst`]: crate::ast::molecule::MoleculeAst

use std::sync::Arc;

use crate::ast::atom::AtomAst;
use crate::ast::bond::BondAst;
use crate::ast::constraint::{AtomConstraint, BondConstraint};
use crate::ast::molecule::MoleculeAst;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AtomPattern {
    pub ast: AtomAst,
    pub constraints: Vec<AtomConstraint>,
}

impl AtomPattern {
    pub fn new(ast: AtomAst) -> Self {
        Self { ast, constraints: Vec::new() }
    }

    pub fn with_constraints(ast: AtomAst, constraints: Vec<AtomConstraint>) -> Self {
        Self { ast, constraints }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BondPattern {
    pub ast: BondAst,
    pub constraints: Vec<BondConstraint>,
}

impl BondPattern {
    pub fn new(ast: BondAst) -> Self {
        Self { ast, constraints: Vec::new() }
    }

    pub fn with_constraints(ast: BondAst, constraints: Vec<BondConstraint>) -> Self {
        Self { ast, constraints }
    }
}

#[derive(Debug)]
struct MoleculeMoleculePatternInner {
    ast: MoleculeAst,
}

#[derive(Clone, Debug)]
pub struct MoleculePattern(Arc<MoleculeMoleculePatternInner>);

impl MoleculePattern {
    pub fn new(ast: MoleculeAst) -> Self {
        Self(Arc::new(MoleculeMoleculePatternInner { ast }))
    }

    pub fn ast(&self) -> &MoleculeAst {
        &self.0.ast
    }
}

impl PartialEq for MoleculePattern {
    fn eq(&self, other: &Self) -> bool {
        self.0.ast == other.0.ast
    }
}

impl Eq for MoleculePattern {}

#[cfg(test)]
mod tests {
    use umol_shared::atom_ast::ElementAst;
    use umol_shared::element::Element;
    use umol_shared::value_ast::ValueAst;

    use super::*;

    #[test]
    fn test_molecule_pattern_new() {
        let ast = MoleculeAst::new(
            vec![AtomAst::new(ElementAst::Undetermined)],
            vec![], vec![], vec![], vec![], vec![], vec![],
        );
        let pattern = MoleculePattern::new(ast);
        assert_eq!(pattern.ast().atoms().count(), 1);
    }

    #[test]
    fn test_atom_pattern_new() {
        let pattern = AtomPattern::new(AtomAst::from_element(Element::C));
        assert_eq!(pattern.ast.element, ElementAst::Lit(Element::C));
        assert!(pattern.constraints.is_empty());
    }

    #[test]
    fn test_atom_pattern_with_constraints() {
        let pattern = AtomPattern::with_constraints(
            AtomAst::from_element(Element::C),
            vec![AtomConstraint::ValenceSum(ValueAst::Lit(4))],
        );
        assert_eq!(pattern.constraints.len(), 1);
        assert_eq!(
            pattern.constraints[0],
            AtomConstraint::ValenceSum(ValueAst::Lit(4)),
        );
    }

    #[test]
    fn test_bond_pattern_new() {
        let pattern = BondPattern::new(BondAst::from_order(1));
        assert_eq!(pattern.ast.order, ValueAst::Lit(1));
        assert!(pattern.constraints.is_empty());
    }

    #[test]
    fn test_bond_pattern_with_constraints() {
        let pattern = BondPattern::with_constraints(
            BondAst::from_order(2),
            vec![BondConstraint::RingBond],
        );
        assert_eq!(pattern.constraints.len(), 1);
        assert_eq!(pattern.constraints[0], BondConstraint::RingBond);
    }
}
