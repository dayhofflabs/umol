//! `ReactionAst` — an owned Python component facade over the Rust reaction AST.

use pyo3::prelude::*;
use umol_ast::ast::ReactionAst as AstReactionAst;

use crate::delta::Deltas;
use crate::molecule::MoleculeAst;

/// A reaction whose molecule and delta components remain live Python values.
#[pyclass]
pub struct ReactionAst {
    lhs: Py<MoleculeAst>,
    deltas: Py<Deltas>,
}

#[pymethods]
impl ReactionAst {
    /// Build a reaction from detached component snapshots.
    #[new]
    #[pyo3(signature = (lhs=None, deltas=None))]
    fn new(
        py: Python<'_>,
        lhs: Option<Py<MoleculeAst>>,
        deltas: Option<Py<Deltas>>,
    ) -> PyResult<Self> {
        Self::from_rust(
            py,
            AstReactionAst::new(
                lhs.map(|value| value.bind(py).borrow().inner().clone())
                    .unwrap_or_default(),
                deltas
                    .map(|value| value.bind(py).borrow().to_rust())
                    .unwrap_or_default(),
            ),
        )
    }

    /// The live left-hand molecule component.
    #[getter]
    fn lhs(&self, py: Python<'_>) -> Py<MoleculeAst> {
        self.lhs.clone_ref(py)
    }

    /// Replace the left-hand molecule with a detached snapshot.
    #[setter]
    fn set_lhs(slf: Py<Self>, py: Python<'_>, value: Py<MoleculeAst>) -> PyResult<()> {
        let resolved = Py::new(
            py,
            MoleculeAst::from_inner(value.bind(py).borrow().inner().clone()),
        )?;
        slf.borrow_mut(py).lhs = resolved;
        Ok(())
    }

    /// The live delta component.
    #[getter]
    fn deltas(&self, py: Python<'_>) -> Py<Deltas> {
        self.deltas.clone_ref(py)
    }

    /// Replace the deltas with a detached snapshot.
    #[setter]
    fn set_deltas(slf: Py<Self>, py: Python<'_>, value: Py<Deltas>) -> PyResult<()> {
        let resolved = Py::new(py, Deltas::from_rust(value.bind(py).borrow().to_rust()))?;
        slf.borrow_mut(py).deltas = resolved;
        Ok(())
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let lhs = self.lhs.bind(py).repr()?.extract::<String>()?;
        let deltas = self.deltas.bind(py).repr()?.extract::<String>()?;
        Ok(format!("ReactionAst(lhs={lhs}, deltas={deltas})"))
    }
}

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
    use pyo3::exceptions::PyTypeError;
    use rstest::rstest;
    use umol_ast::ast::{
        AtomAst as AstAtomAst, AtomDelta as AstAtomDelta, AtomId as AstAtomId, Delta as AstDelta,
        Deltas as AstDeltas, MoleculeAst as AstMoleculeAst, MoleculeParts as AstMoleculeParts,
    };
    use umol_chem::element::Element as ChemElement;

    use super::*;
    use crate::convert::into_py_variant;
    use crate::delta::Delta;

    #[rstest]
    #[case::empty(None, None, AstReactionAst::default())]
    #[case::populated(
        Some(AstMoleculeAst::from_parts(AstMoleculeParts {
            atoms: vec![AstAtomAst::from_element(ChemElement::C)],
            ..Default::default()
        })),
        Some(vec![AstDelta::Atom(AstAtomDelta::Add {
            id: AstAtomId(1),
            ast: AstAtomAst::from_element(ChemElement::O),
        })].into_iter().collect()),
        AstReactionAst::new(
            AstMoleculeAst::from_parts(AstMoleculeParts {
                atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                ..Default::default()
            }),
            vec![AstDelta::Atom(AstAtomDelta::Add {
                id: AstAtomId(1),
                ast: AstAtomAst::from_element(ChemElement::O),
            })].into_iter().collect(),
        ),
    )]
    fn test_reaction_ast_new(
        #[case] lhs: Option<AstMoleculeAst>,
        #[case] deltas: Option<AstDeltas>,
        #[case] expected: AstReactionAst,
    ) {
        Python::attach(|py| {
            let lhs = lhs.map(|value| Py::new(py, MoleculeAst::from_inner(value)).unwrap());
            let deltas = deltas.map(|value| Py::new(py, Deltas::from_rust(value)).unwrap());

            let reaction = ReactionAst::new(py, lhs, deltas).unwrap();

            assert_eq!(reaction.to_rust(py), expected);
        });
    }

    #[rstest]
    fn test_reaction_ast_new_snapshot() {
        Python::attach(|py| {
            let lhs = Py::new(
                py,
                MoleculeAst::from_inner(AstMoleculeAst::from_parts(AstMoleculeParts {
                    atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                    ..Default::default()
                })),
            )
            .unwrap();
            let deltas = Py::new(
                py,
                Deltas::from_rust(
                    vec![AstDelta::Atom(AstAtomDelta::Add {
                        id: AstAtomId(1),
                        ast: AstAtomAst::from_element(ChemElement::O),
                    })]
                    .into_iter()
                    .collect(),
                ),
            )
            .unwrap();
            let expected = AstReactionAst::new(
                lhs.bind(py).borrow().inner().clone(),
                deltas.bind(py).borrow().to_rust(),
            );

            let reaction =
                ReactionAst::new(py, Some(lhs.clone_ref(py)), Some(deltas.clone_ref(py))).unwrap();
            *lhs.bind(py).borrow_mut().inner_mut() = AstMoleculeAst::new();
            let delta = into_py_variant(
                py,
                Delta::from_rust(
                    py,
                    &AstDelta::Atom(AstAtomDelta::Add {
                        id: AstAtomId(2),
                        ast: AstAtomAst::from_element(ChemElement::N),
                    }),
                )
                .unwrap(),
            )
            .unwrap();
            deltas.bind(py).call_method1("append", (delta,)).unwrap();

            assert_eq!(reaction.to_rust(py), expected);
            assert_ne!(reaction.lhs.as_ptr(), lhs.as_ptr());
            assert_ne!(reaction.deltas.as_ptr(), deltas.as_ptr());
        });
    }

    #[rstest]
    fn test_reaction_ast_components() {
        Python::attach(|py| {
            let reaction = Py::new(py, ReactionAst::new(py, None, None).unwrap()).unwrap();
            let first_lhs = reaction.bind(py).borrow().lhs(py);
            let second_lhs = reaction.bind(py).borrow().lhs(py);
            let first_deltas = reaction.bind(py).borrow().deltas(py);
            let second_deltas = reaction.bind(py).borrow().deltas(py);

            *first_lhs.bind(py).borrow_mut().inner_mut() =
                AstMoleculeAst::from_parts(AstMoleculeParts {
                    atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                    ..Default::default()
                });
            let delta = into_py_variant(
                py,
                Delta::from_rust(
                    py,
                    &AstDelta::Atom(AstAtomDelta::Add {
                        id: AstAtomId(1),
                        ast: AstAtomAst::from_element(ChemElement::O),
                    }),
                )
                .unwrap(),
            )
            .unwrap();
            first_deltas
                .bind(py)
                .call_method1("append", (delta,))
                .unwrap();

            assert_eq!(first_lhs.as_ptr(), second_lhs.as_ptr());
            assert_eq!(first_deltas.as_ptr(), second_deltas.as_ptr());
            assert_eq!(
                reaction.bind(py).borrow().to_rust(py),
                AstReactionAst::new(
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
                )
            );
        });
    }

    #[rstest]
    fn test_reaction_ast_set_components() {
        Python::attach(|py| {
            let reaction = Py::new(py, ReactionAst::new(py, None, None).unwrap()).unwrap();
            let lhs = Py::new(
                py,
                MoleculeAst::from_inner(AstMoleculeAst::from_parts(AstMoleculeParts {
                    atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                    ..Default::default()
                })),
            )
            .unwrap();
            let deltas = Py::new(
                py,
                Deltas::from_rust(
                    vec![AstDelta::Atom(AstAtomDelta::Add {
                        id: AstAtomId(1),
                        ast: AstAtomAst::from_element(ChemElement::O),
                    })]
                    .into_iter()
                    .collect(),
                ),
            )
            .unwrap();
            let expected = AstReactionAst::new(
                lhs.bind(py).borrow().inner().clone(),
                deltas.bind(py).borrow().to_rust(),
            );

            ReactionAst::set_lhs(reaction.clone_ref(py), py, lhs.clone_ref(py)).unwrap();
            ReactionAst::set_deltas(reaction.clone_ref(py), py, deltas.clone_ref(py)).unwrap();
            *lhs.bind(py).borrow_mut().inner_mut() = AstMoleculeAst::new();
            let delta = into_py_variant(
                py,
                Delta::from_rust(
                    py,
                    &AstDelta::Atom(AstAtomDelta::Add {
                        id: AstAtomId(2),
                        ast: AstAtomAst::from_element(ChemElement::N),
                    }),
                )
                .unwrap(),
            )
            .unwrap();
            deltas.bind(py).call_method1("append", (delta,)).unwrap();

            assert_eq!(reaction.bind(py).borrow().to_rust(py), expected);
        });
    }

    #[rstest]
    fn test_reaction_ast_set_components_self() {
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
            let reaction =
                Py::new(py, ReactionAst::from_rust(py, expected.clone()).unwrap()).unwrap();
            let own_lhs = reaction.bind(py).borrow().lhs(py);
            let own_deltas = reaction.bind(py).borrow().deltas(py);

            ReactionAst::set_lhs(reaction.clone_ref(py), py, own_lhs).unwrap();
            ReactionAst::set_deltas(reaction.clone_ref(py), py, own_deltas).unwrap();

            assert_eq!(reaction.bind(py).borrow().to_rust(py), expected);
        });
    }

    #[rstest]
    fn test_reaction_ast_eq() {
        Python::attach(|py| {
            let empty = ReactionAst::new(py, None, None).unwrap();
            let other_empty = ReactionAst::new(py, None, None).unwrap();
            let populated = ReactionAst::from_rust(
                py,
                AstReactionAst::new(
                    AstMoleculeAst::from_parts(AstMoleculeParts {
                        atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                        ..Default::default()
                    }),
                    AstDeltas::new(),
                ),
            )
            .unwrap();

            assert!(empty.__eq__(&other_empty, py));
            assert!(!empty.__eq__(&populated, py));
            let empty = Py::new(py, empty).unwrap();
            assert!(empty
                .bind(py)
                .hash()
                .unwrap_err()
                .is_instance_of::<PyTypeError>(py));
        });
    }

    #[rstest]
    fn test_reaction_ast_repr() {
        Python::attach(|py| {
            let reaction = ReactionAst::from_rust(
                py,
                AstReactionAst::new(
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
                ),
            )
            .unwrap();

            assert_eq!(
                reaction.__repr__(py).unwrap(),
                "ReactionAst(lhs=MoleculeAst(atoms=1, bonds=0), deltas=Deltas([Delta.Atom(AtomDelta.Add(id=1, ast=AtomAst.parse('O')))]))"
            );
        });
    }

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
