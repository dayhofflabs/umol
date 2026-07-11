//! `MoleculeAst` — a molecule (owned graph-AST root), wrapping
//! `umol_ast::ast::MoleculeAst`.

use pyo3::prelude::*;
use umol_ast::ast::{AtomId as AstAtomId, MoleculeAst as AstMoleculeAst};

use crate::atom::{AtomAst, AtomViews};
use crate::bond::{BondAst, BondViews};

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

    /// A molecule from a sequence of atoms and bonds. Each bond is a
    /// `(first, second, bond)` triple: two atom indices into `atoms` and a
    /// `BondAst` value carrying the bond's order/charge/spin/constraints.
    #[staticmethod]
    #[pyo3(signature = (atoms, bonds=Vec::new()))]
    fn from_atoms_and_bonds(
        py: Python<'_>,
        atoms: Vec<Py<AtomAst>>,
        bonds: Vec<(u32, u32, Py<BondAst>)>,
    ) -> Self {
        let ast_atoms = atoms
            .iter()
            .map(|atom| atom.bind(py).borrow().inner().clone())
            .collect();
        let ast_bonds = bonds
            .iter()
            .map(|(first, second, bond)| {
                (
                    AstAtomId(*first),
                    AstAtomId(*second),
                    bond.bind(py).borrow().inner().clone(),
                )
            })
            .collect();
        MoleculeAst(AstMoleculeAst::from_atoms_and_bonds(ast_atoms, ast_bonds))
    }

    /// The atoms, indexed by integer position.
    #[getter]
    fn atoms(slf: Py<Self>) -> AtomViews {
        AtomViews::new(slf)
    }

    /// The bonds, indexed by integer position.
    #[getter]
    fn bonds(slf: Py<Self>) -> BondViews {
        BondViews::new(slf)
    }

    fn __repr__(&self) -> String {
        format!(
            "MoleculeAst(atoms={}, bonds={})",
            self.0.atoms().count(),
            self.0.bonds().count()
        )
    }
}

impl MoleculeAst {
    /// The wrapped AST molecule — read access for atom views.
    pub(crate) fn inner(&self) -> &AstMoleculeAst {
        &self.0
    }

    /// Mutable access to the wrapped AST molecule — write access for the live
    /// atom and constraint views (copy-on-write through `atom_mut`).
    pub(crate) fn inner_mut(&mut self) -> &mut AstMoleculeAst {
        &mut self.0
    }

    /// Wrap an AST molecule (the hold-the-value `from_inner` bridge, paired with
    /// `inner`). Test-only — in-crate construction wraps `MoleculeAst(..)` directly.
    #[cfg(test)]
    pub(crate) fn from_inner(molecule: AstMoleculeAst) -> Self {
        MoleculeAst(molecule)
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
        assert_eq!(MoleculeAst::new().inner().atoms().count(), 0);
    }

    #[rstest]
    #[case(vec![], 0)]
    #[case(vec![Element::C], 1)]
    #[case(vec![Element::C, Element::O], 2)]
    fn test_molecule_ast_atoms(#[case] elements: Vec<Element>, #[case] expected: usize) {
        let atoms = elements.into_iter().map(AtomAst::from_element).collect();
        let molecule = MoleculeAst(AstMoleculeAst::from_atoms_and_bonds(atoms, vec![]));
        assert_eq!(molecule.inner().atoms().count(), expected);
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
        assert_eq!(MoleculeAst::new().__repr__(), "MoleculeAst(atoms=0, bonds=0)");
    }
}
