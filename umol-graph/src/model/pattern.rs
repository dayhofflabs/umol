//! Substructure-query pattern backed by [`MoleculeAst`].
//!
//! [`MoleculeAst`]: crate::ast::molecule::MoleculeAst

use std::sync::Arc;

use crate::ast::molecule::MoleculeAst;

#[derive(Debug)]
struct PatternInner {
    ast: MoleculeAst,
}

#[derive(Clone, Debug)]
pub struct Pattern(Arc<PatternInner>);

impl Pattern {
    pub fn new(ast: MoleculeAst) -> Self {
        Self(Arc::new(PatternInner { ast }))
    }

    pub fn ast(&self) -> &MoleculeAst {
        &self.0.ast
    }
}

impl PartialEq for Pattern {
    fn eq(&self, other: &Self) -> bool {
        self.0.ast == other.0.ast
    }
}

impl Eq for Pattern {}

#[cfg(test)]
mod tests {
    use umol_shared::atom_ast::ElementAst;

    use super::*;
    use crate::ast::atom::AtomAst;

    #[test]
    fn test_pattern_new() {
        let ast = MoleculeAst::new(
            vec![AtomAst::new(ElementAst::Undetermined)],
            vec![], vec![], vec![], vec![], vec![], vec![],
        );
        let pattern = Pattern::new(ast);
        assert_eq!(pattern.ast().atoms().count(), 1);
    }
}
