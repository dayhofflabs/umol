//! `MoleculeAst` — a molecule (owned graph-AST root), wrapping
//! `umol_ast::ast::MoleculeAst`.

use pyo3::prelude::*;
use umol_ast::ast::MoleculeAst as AstMoleculeAst;

/// A molecule: the owned graph-AST root.
#[pyclass(eq)]
#[derive(Debug, PartialEq)]
pub struct MoleculeAst(AstMoleculeAst);

#[pymethods]
impl MoleculeAst {
    /// An empty molecule: zero atoms, zero bonds.
    #[new]
    fn new() -> Self {
        Self(AstMoleculeAst::new())
    }

    /// Number of atoms.
    #[getter]
    fn atom_count(&self) -> u32 {
        self.0.atoms().count() as u32
    }

    fn __repr__(&self) -> String {
        format!("MoleculeAst(atom_count={})", self.0.atoms().count())
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::AtomAst;
    use umol_chem::element::Element;

    use super::*;

    #[rstest]
    fn test_molecule_ast_new() {
        assert_eq!(MoleculeAst::new().atom_count(), 0);
    }

    #[rstest]
    #[case(vec![], 0)]
    #[case(vec![Element::C], 1)]
    #[case(vec![Element::C, Element::O], 2)]
    fn test_molecule_ast_atom_count(#[case] elements: Vec<Element>, #[case] expected: u32) {
        let atoms = elements.into_iter().map(AtomAst::from_element).collect();
        let molecule = MoleculeAst(AstMoleculeAst::from_atoms_and_bonds(atoms, vec![]));
        assert_eq!(molecule.atom_count(), expected);
    }

    #[rstest]
    fn test_molecule_ast_eq() {
        assert_eq!(MoleculeAst::new(), MoleculeAst::new());
        let carbon = MoleculeAst(AstMoleculeAst::from_atoms_and_bonds(
            vec![AtomAst::from_element(Element::C)],
            vec![],
        ));
        assert_ne!(MoleculeAst::new(), carbon);
    }

    #[rstest]
    fn test_molecule_ast_repr() {
        assert_eq!(MoleculeAst::new().__repr__(), "MoleculeAst(atom_count=0)");
    }
}
