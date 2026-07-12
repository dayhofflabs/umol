//! `MoleculeAst` — a molecule (owned graph-AST root), wrapping
//! `umol_ast::ast::MoleculeAst`.

use pyo3::prelude::*;
use umol_ast::ast::{
    AtomId as AstAtomId, MoleculeAst as AstMoleculeAst, MoleculeParts as AstMoleculeParts,
};

use crate::aromatic::AromaticSystemAst;
use crate::atom::{AtomAst, AtomViews};
use crate::bond::{BondAst, BondViews};
use crate::dative::{DativeBondAst, DativeBondViews};

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

    /// A molecule from its parts. Each bond is a `(first, second, bond)` triple:
    /// two atom indices into `atoms` and a `BondAst`. Each dative bond is a
    /// `(donors, acceptor, bond)` triple: a list of donor atom indices, one
    /// acceptor atom index, and a `DativeBondAst`. Each aromatic system is an
    /// `(atoms, system)` pair: a list of member atom indices and an `AromaticSystemAst`.
    #[staticmethod]
    #[pyo3(signature = (atoms, *, bonds=Vec::new(), dative=Vec::new(), aromatic=Vec::new()))]
    fn from_parts(
        py: Python<'_>,
        atoms: Vec<Py<AtomAst>>,
        bonds: Vec<(u32, u32, Py<BondAst>)>,
        dative: Vec<(Vec<u32>, u32, Py<DativeBondAst>)>,
        aromatic: Vec<(Vec<u32>, Py<AromaticSystemAst>)>,
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
        let ast_dative = dative
            .iter()
            .map(|(donors, acceptor, bond)| {
                (
                    donors.iter().map(|&donor| AstAtomId(donor)).collect(),
                    AstAtomId(*acceptor),
                    bond.bind(py).borrow().inner().clone(),
                )
            })
            .collect();
        let ast_aromatic = aromatic
            .iter()
            .map(|(atoms, system)| {
                (
                    atoms.iter().map(|&atom| AstAtomId(atom)).collect(),
                    system.bind(py).borrow().inner().clone(),
                )
            })
            .collect();
        MoleculeAst(AstMoleculeAst::from_parts(AstMoleculeParts {
            atoms: ast_atoms,
            bonds: ast_bonds,
            dative: ast_dative,
            aromatic: ast_aromatic,
            ..Default::default()
        }))
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

    /// The dative bonds, indexed by integer position.
    #[getter]
    fn dative_bonds(slf: Py<Self>) -> DativeBondViews {
        DativeBondViews::new(slf)
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
    use umol_ast::ast::{
        AromaticSystemAst as AstAromaticSystemAst, AromaticSystemId as AstAromaticSystemId,
        AtomAst, BondAst as AstBondAst, DativeBondAst as AstDativeBondAst,
        DativeBondId as AstDativeBondId,
    };
    use umol_chem::element::Element;

    use super::*;
    use crate::atom::AtomAst as PyAtomAst;

    #[rstest]
    fn test_molecule_ast_new() {
        assert_eq!(MoleculeAst::new().inner().atoms().count(), 0);
    }

    #[rstest]
    fn test_molecule_ast_from_parts() {
        Python::attach(|py| {
            let atoms = vec![
                Py::new(py, PyAtomAst::from_inner(AtomAst::from_element(Element::C))).unwrap(),
                Py::new(py, PyAtomAst::from_inner(AtomAst::from_element(Element::B))).unwrap(),
                Py::new(py, PyAtomAst::from_inner(AtomAst::from_element(Element::N))).unwrap(),
            ];
            let bonds = vec![(
                0,
                1,
                Py::new(py, BondAst::from_inner(AstBondAst::from_order(1))).unwrap(),
            )];
            let dative = vec![(
                vec![2],
                1,
                Py::new(
                    py,
                    DativeBondAst::from_inner(AstDativeBondAst::from_order(1)),
                )
                .unwrap(),
            )];
            let aromatic = vec![(
                vec![0, 1, 2],
                Py::new(
                    py,
                    AromaticSystemAst::from_inner(AstAromaticSystemAst::from_electrons(vec![
                        1, 1, 1,
                    ])),
                )
                .unwrap(),
            )];
            let molecule = MoleculeAst::from_parts(py, atoms, bonds, dative, aromatic);
            assert_eq!(molecule.inner().atoms().count(), 3);
            assert_eq!(molecule.inner().bonds().count(), 1);
            let dative_bonds = molecule.inner().dative_bonds();
            assert_eq!(dative_bonds.count(), 1);
            let dative_view = dative_bonds.get(AstDativeBondId(0)).unwrap();
            assert_eq!(dative_view.acceptor_id(), AstAtomId(1));
            assert_eq!(
                dative_view.donor_ids().collect::<Vec<_>>(),
                vec![AstAtomId(2)]
            );
            let aromatic_systems = molecule.inner().aromatic_systems();
            assert_eq!(aromatic_systems.count(), 1);
            let aromatic_view = aromatic_systems.get(AstAromaticSystemId(0)).unwrap();
            assert_eq!(
                aromatic_view.atom_ids().collect::<Vec<_>>(),
                vec![AstAtomId(0), AstAtomId(1), AstAtomId(2)]
            );
        });
    }

    #[rstest]
    #[case(vec![], 0)]
    #[case(vec![Element::C], 1)]
    #[case(vec![Element::C, Element::O], 2)]
    fn test_molecule_ast_atoms(#[case] elements: Vec<Element>, #[case] expected: usize) {
        let atoms = elements.into_iter().map(AtomAst::from_element).collect();
        let molecule = MoleculeAst(AstMoleculeAst::from_parts(AstMoleculeParts {
            atoms,
            ..Default::default()
        }));
        assert_eq!(molecule.inner().atoms().count(), expected);
    }

    #[rstest]
    fn test_molecule_ast_eq() {
        assert_eq!(MoleculeAst::new(), MoleculeAst::new());
        let carbon = MoleculeAst(AstMoleculeAst::from_parts(AstMoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C)],
            ..Default::default()
        }));
        assert_ne!(MoleculeAst::new(), carbon);
    }

    #[rstest]
    fn test_molecule_ast_repr() {
        assert_eq!(
            MoleculeAst::new().__repr__(),
            "MoleculeAst(atoms=0, bonds=0)"
        );
    }
}
