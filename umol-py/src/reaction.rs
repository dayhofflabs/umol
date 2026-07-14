//! `ReactionAst` — an owned Python component facade over the Rust reaction AST.

use pyo3::prelude::*;
use umol_ast::ast::ReactionAst as AstReactionAst;

use crate::delta::Deltas;
use crate::molecule::MoleculeAst;

/// A reaction whose molecule and delta components remain live Python values.
#[pyclass]
#[allow(dead_code, reason = "owned component kernel for reaction operations")]
pub struct ReactionAst {
    lhs: Py<MoleculeAst>,
    deltas: Py<Deltas>,
}

#[allow(dead_code, reason = "owned component kernel for reaction operations")]
impl ReactionAst {
    /// Wrap a Rust reaction in fresh Python-owned components.
    pub(crate) fn from_rust(py: Python<'_>, reaction: AstReactionAst) -> PyResult<Self> {
        Ok(Self {
            lhs: Py::new(py, MoleculeAst::from_inner(reaction.lhs))?,
            deltas: Py::new(py, Deltas::from_rust(reaction.deltas))?,
        })
    }

    /// Snapshot the current Python-owned components as a Rust reaction.
    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstReactionAst {
        AstReactionAst::new(
            self.lhs.bind(py).borrow().inner().clone(),
            self.deltas.bind(py).borrow().to_rust(),
        )
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::{
        AtomAst as AstAtomAst, AtomDelta as AstAtomDelta, AtomId as AstAtomId, Delta as AstDelta,
        Deltas as AstDeltas, MoleculeAst as AstMoleculeAst, MoleculeParts as AstMoleculeParts,
    };
    use umol_chem::element::Element as ChemElement;

    use super::*;

    #[rstest]
    #[case::empty(AstReactionAst::default())]
    #[case::populated(AstReactionAst::new(
        AstMoleculeAst::from_parts(AstMoleculeParts {
            atoms: vec![AstAtomAst::from_element(ChemElement::C)],
            ..Default::default()
        }),
        vec![AstDelta::Atom(AstAtomDelta::Add {
            id: AstAtomId(1),
            ast: AstAtomAst::from_element(ChemElement::O),
        })]
        .into_iter()
        .collect(),
    ))]
    fn test_reaction_ast_from_rust(#[case] expected: AstReactionAst) {
        Python::attach(|py| {
            let reaction = ReactionAst::from_rust(py, expected.clone()).unwrap();

            assert_eq!(reaction.to_rust(py), expected);
        });
    }

    #[rstest]
    fn test_reaction_ast_to_rust() {
        Python::attach(|py| {
            let expected = AstReactionAst::new(
                AstMoleculeAst::from_parts(AstMoleculeParts {
                    atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                    ..Default::default()
                }),
                vec![AstDelta::Atom(AstAtomDelta::Add {
                    id: AstAtomId(1),
                    ast: AstAtomAst::from_element(ChemElement::O),
                })]
                .into_iter()
                .collect(),
            );
            let reaction = ReactionAst::from_rust(py, expected.clone()).unwrap();

            let mut snapshot = reaction.to_rust(py);
            snapshot.lhs = AstMoleculeAst::new();
            snapshot.deltas = AstDeltas::new();

            assert_eq!(reaction.to_rust(py), expected);
        });
    }
}
