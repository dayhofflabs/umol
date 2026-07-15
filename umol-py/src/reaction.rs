//! `ReactionAst` — an owned Python component facade over the Rust reaction AST.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use std::collections::HashSet;
use std::str::FromStr;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use umol_ast::ast::{
    Canonicalize, CompositionScope as AstCompositionScope, ReactionAst as AstReactionAst,
};
use umol_graph_core::{Correspondence, NodeId};

use crate::delta::Deltas;
use crate::error::{contradiction_error, parse_error};
use crate::molecule::MoleculeAst;

/// Which overlaps sequential reaction composition retains.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum CompositionScope {
    RcAnchored,
    Full,
}

#[allow(
    dead_code,
    reason = "Rust/Python conversion API for composition scope values"
)]
impl CompositionScope {
    pub(crate) fn from_rust(scope: AstCompositionScope) -> Self {
        match scope {
            AstCompositionScope::RcAnchored => Self::RcAnchored,
            AstCompositionScope::Full => Self::Full,
        }
    }

    pub(crate) fn to_rust(self) -> AstCompositionScope {
        match self {
            Self::RcAnchored => AstCompositionScope::RcAnchored,
            Self::Full => AstCompositionScope::Full,
        }
    }
}

/// Validate atom pairs and construct their partial bijection over the two side sizes.
fn atom_correspondence(
    pairs: Vec<(usize, usize)>,
    lhs_count: usize,
    rhs_count: usize,
) -> PyResult<Correspondence<NodeId>> {
    let mut left_ids = HashSet::with_capacity(pairs.len());
    let mut right_ids = HashSet::with_capacity(pairs.len());
    let mut mates = Vec::with_capacity(pairs.len());

    for (left, right) in pairs {
        if left >= lhs_count {
            return Err(PyValueError::new_err(format!(
                "left atom id {left} out of range for {lhs_count} atoms"
            )));
        }
        if right >= rhs_count {
            return Err(PyValueError::new_err(format!(
                "right atom id {right} out of range for {rhs_count} atoms"
            )));
        }
        if !left_ids.insert(left) {
            return Err(PyValueError::new_err(format!(
                "duplicate left atom id {left}"
            )));
        }
        if !right_ids.insert(right) {
            return Err(PyValueError::new_err(format!(
                "duplicate right atom id {right}"
            )));
        }
        mates.push((NodeId::from(left), NodeId::from(right)));
    }

    Ok(Correspondence::new(mates, lhs_count, rhs_count))
}

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

    /// Parse a reaction from its EDN representation.
    #[staticmethod]
    fn parse(py: Python<'_>, text: &str) -> PyResult<Self> {
        let reaction = AstReactionAst::from_str(text).map_err(parse_error)?;
        Self::from_rust(py, reaction)
    }

    /// Construct a reaction by comparing two molecule snapshots under an atom correspondence.
    #[staticmethod]
    fn from_sides(
        py: Python<'_>,
        lhs: Py<MoleculeAst>,
        rhs: Py<MoleculeAst>,
        atom_pairs: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        let lhs = lhs.bind(py).borrow().inner().clone();
        let rhs = rhs.bind(py).borrow().inner().clone();
        let atom_pairs = atom_pairs
            .try_iter()?
            .map(|item| item?.extract::<(usize, usize)>())
            .collect::<PyResult<Vec<_>>>()?;
        let atom = atom_correspondence(atom_pairs, lhs.atoms().count(), rhs.atoms().count())?;

        Self::from_rust(py, AstReactionAst::from_sides(lhs, rhs, atom))
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

    /// Return a fresh canonical reaction, leaving this facade unchanged.
    fn canonicalize(&self, py: Python<'_>) -> PyResult<Self> {
        let reaction = self
            .to_rust(py)
            .canonicalize()
            .map_err(contradiction_error)?;
        Self::from_rust(py, reaction)
    }

    /// Return the reverse reaction in the product's compacted id space.
    fn reverse(&self, py: Python<'_>) -> PyResult<Self> {
        let reaction = self.to_rust(py).reverse().map_err(contradiction_error)?;
        Self::from_rust(py, reaction)
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __str__(&self, py: Python<'_>) -> String {
        self.to_rust(py).to_string()
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
    use pyo3::exceptions::{PyTypeError, PyValueError};
    use pyo3::types::PyList;
    use rstest::rstest;
    use umol_ast::ast::{
        AromaticSystemAst as AstAromaticSystemAst, AromaticSystemDelta as AstAromaticSystemDelta,
        AromaticSystemId as AstAromaticSystemId, AtomAst as AstAtomAst, AtomDelta as AstAtomDelta,
        AtomFieldChange as AstAtomFieldChange, AtomId as AstAtomId, BondAst as AstBondAst,
        BondDelta as AstBondDelta, BondFieldChange as AstBondFieldChange, BondId as AstBondId,
        Constraint as AstConstraint, ConstraintDelta as AstConstraintDelta,
        DativeBondAst as AstDativeBondAst, DativeBondDelta as AstDativeBondDelta,
        DativeBondId as AstDativeBondId, Delta as AstDelta, Deltas as AstDeltas,
        MoleculeAst as AstMoleculeAst, MoleculeConstraint as AstMoleculeConstraint,
        MoleculeParts as AstMoleculeParts, MulticenterBondAst as AstMulticenterBondAst,
        MulticenterBondDelta as AstMulticenterBondDelta, MulticenterBondId as AstMulticenterBondId,
        NoncovalentBondAst as AstNoncovalentBondAst,
        NoncovalentBondDelta as AstNoncovalentBondDelta, NoncovalentBondId as AstNoncovalentBondId,
        NoncovalentBondKind as AstNoncovalentBondKind, StereoAtomAst as AstStereoAtomAst,
        StereoAtomDelta as AstStereoAtomDelta, StereoAtomId as AstStereoAtomId,
        StereoBondAst as AstStereoBondAst, StereoBondDelta as AstStereoBondDelta,
        StereoBondId as AstStereoBondId, StereoCosetAst as AstStereoCosetAst,
        StereoKind as AstStereoKind, StereoLigand as AstStereoLigand,
        StereoLigandKind as AstStereoLigandKind, ValueAst as AstValueAst,
    };
    use umol_chem::element::Element as ChemElement;

    use super::*;
    use crate::convert::into_py_variant;
    use crate::delta::Delta;
    use crate::error::{ContradictionError, ParseError};

    #[rstest]
    #[case::rc_anchored(AstCompositionScope::RcAnchored, CompositionScope::RcAnchored)]
    #[case::full(AstCompositionScope::Full, CompositionScope::Full)]
    fn test_composition_scope_from_rust(
        #[case] scope: AstCompositionScope,
        #[case] expected: CompositionScope,
    ) {
        assert_eq!(CompositionScope::from_rust(scope), expected);
    }

    #[rstest]
    #[case::rc_anchored(CompositionScope::RcAnchored, AstCompositionScope::RcAnchored)]
    #[case::full(CompositionScope::Full, AstCompositionScope::Full)]
    fn test_composition_scope_to_rust(
        #[case] scope: CompositionScope,
        #[case] expected: AstCompositionScope,
    ) {
        assert_eq!(scope.to_rust(), expected);
    }

    #[rstest]
    #[case::rc_anchored(
        CompositionScope::RcAnchored,
        CompositionScope::RcAnchored,
        CompositionScope::Full,
        "CompositionScope.RcAnchored"
    )]
    #[case::full(
        CompositionScope::Full,
        CompositionScope::Full,
        CompositionScope::RcAnchored,
        "CompositionScope.Full"
    )]
    fn test_composition_scope_python_value(
        #[case] scope: CompositionScope,
        #[case] equal: CompositionScope,
        #[case] unequal: CompositionScope,
        #[case] expected_repr: &str,
    ) {
        Python::attach(|py| {
            let scope = Py::new(py, scope).unwrap();
            let equal = Py::new(py, equal).unwrap();
            let unequal = Py::new(py, unequal).unwrap();

            assert!(scope.bind(py).as_any().eq(equal.bind(py).as_any()).unwrap());
            assert!(!scope
                .bind(py)
                .as_any()
                .eq(unequal.bind(py).as_any())
                .unwrap());
            assert_eq!(
                scope.bind(py).as_any().hash().unwrap(),
                equal.bind(py).as_any().hash().unwrap()
            );
            assert_eq!(
                scope
                    .bind(py)
                    .as_any()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                expected_repr
            );
        });
    }

    #[rstest]
    #[case::empty(Vec::new(), 0, 0, Vec::new())]
    #[case::partial(vec![(1, 2)], 3, 4, vec![(NodeId(1), NodeId(2))])]
    #[case::total(
        vec![(0, 1), (1, 0)],
        2,
        2,
        vec![(NodeId(0), NodeId(1)), (NodeId(1), NodeId(0))],
    )]
    #[case::unsorted(
        vec![(2, 0), (0, 2)],
        3,
        3,
        vec![(NodeId(0), NodeId(2)), (NodeId(2), NodeId(0))],
    )]
    fn test_atom_correspondence(
        #[case] pairs: Vec<(usize, usize)>,
        #[case] lhs_count: usize,
        #[case] rhs_count: usize,
        #[case] expected_mates: Vec<(NodeId, NodeId)>,
    ) {
        let correspondence = atom_correspondence(pairs, lhs_count, rhs_count).unwrap();

        assert_eq!(correspondence.mates(), expected_mates.as_slice());
        assert_eq!(correspondence.left_count(), lhs_count);
        assert_eq!(correspondence.right_count(), rhs_count);
    }

    #[rstest]
    #[case::duplicate_left(
        vec![(0, 0), (0, 1)],
        2,
        2,
        "duplicate left atom id 0",
    )]
    #[case::duplicate_right(
        vec![(0, 1), (1, 1)],
        2,
        2,
        "duplicate right atom id 1",
    )]
    #[case::left_out_of_range(
        vec![(2, 0)],
        2,
        1,
        "left atom id 2 out of range for 2 atoms",
    )]
    #[case::right_out_of_range(
        vec![(0, 1)],
        1,
        1,
        "right atom id 1 out of range for 1 atoms",
    )]
    fn test_atom_correspondence_error(
        #[case] pairs: Vec<(usize, usize)>,
        #[case] lhs_count: usize,
        #[case] rhs_count: usize,
        #[case] expected: &str,
    ) {
        Python::attach(|py| {
            let error = atom_correspondence(pairs, lhs_count, rhs_count)
                .err()
                .unwrap();

            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                expected
            );
        });
    }

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
    #[case::atom_add_remove(
        r##"{:lhs {:atoms ["C" "O"]} :deltas [{:atom {:add "N"}} {:atom {:remove 1}}]}"##,
        2,
        vec![
            AstDelta::Atom(AstAtomDelta::Add {
                id: AstAtomId(2),
                ast: AstAtomAst::from_element(ChemElement::N),
            }),
            AstDelta::Atom(AstAtomDelta::Remove {
                id: AstAtomId(1),
                ast: AstAtomAst::from_element(ChemElement::O),
            }),
        ],
    )]
    #[case::atom_modify(
        r##"{:lhs {:atoms ["Br#c0"]} :deltas [{:atom {:modify [0 "#c-1"]}}]}"##,
        1,
        vec![AstDelta::Atom(AstAtomDelta::ModifyField {
            id: AstAtomId(0),
            change: AstAtomFieldChange::Charge {
                old: AstValueAst::Lit(0),
                new: AstValueAst::Lit(-1),
            },
        })],
    )]
    #[case::stereo_mirror(
        r##"{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]} :deltas [{:stereo-atom {:mirror [0 :tetrahedral]}}]}"##,
        5,
        vec![AstDelta::StereoAtom(AstStereoAtomDelta::Mirror {
            id: AstStereoAtomId(0),
            kind: AstStereoKind::Tetrahedral,
        })],
    )]
    #[case::molecule_constraint(
        r##"{:lhs {:atoms ["C"]} :deltas [{:constraint {:add {:connected {}}}}]}"##,
        1,
        vec![AstDelta::Constraint(AstConstraintDelta::Add(
            AstConstraint::Molecule(AstMoleculeConstraint::Connected { atoms: None }),
        ))],
    )]
    fn test_reaction_ast_parse(
        #[case] text: &str,
        #[case] atom_count: usize,
        #[case] expected_deltas: Vec<AstDelta>,
    ) {
        Python::attach(|py| {
            let reaction = ReactionAst::parse(py, text).unwrap().to_rust(py);

            assert_eq!(reaction.lhs.atoms().count(), atom_count);
            assert_eq!(reaction.deltas.as_slice(), expected_deltas.as_slice());
        });
    }

    #[rstest]
    fn test_reaction_ast_parse_error() {
        Python::attach(|py| {
            let error = ReactionAst::parse(py, "not edn").err().unwrap();

            assert!(error.is_instance_of::<ParseError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "EDN parse: unexpected token 'n' at byte 0"
            );
        });
    }

    #[rstest]
    #[case::identity(
        AstMoleculeAst::from_parts(AstMoleculeParts {
            atoms: vec![AstAtomAst::from_element(ChemElement::C)],
            ..Default::default()
        }),
        AstMoleculeAst::from_parts(AstMoleculeParts {
            atoms: vec![AstAtomAst::from_element(ChemElement::C)],
            ..Default::default()
        }),
        vec![(0, 0)],
        AstReactionAst::new(
            AstMoleculeAst::from_parts(AstMoleculeParts {
                atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                ..Default::default()
            }),
            AstDeltas::default(),
        ),
    )]
    #[case::partial_correspondence(
        AstMoleculeAst::from_parts(AstMoleculeParts {
            atoms: vec![
                AstAtomAst::from_element(ChemElement::C),
                AstAtomAst::from_element(ChemElement::O),
            ],
            ..Default::default()
        }),
        AstMoleculeAst::from_parts(AstMoleculeParts {
            atoms: vec![
                AstAtomAst::from_element(ChemElement::C),
                AstAtomAst::from_element(ChemElement::N),
            ],
            ..Default::default()
        }),
        vec![(0, 0)],
        AstReactionAst::new(
            AstMoleculeAst::from_parts(AstMoleculeParts {
                atoms: vec![
                    AstAtomAst::from_element(ChemElement::C),
                    AstAtomAst::from_element(ChemElement::O),
                ],
                ..Default::default()
            }),
            vec![
                AstDelta::Atom(AstAtomDelta::Remove {
                    id: AstAtomId(1),
                    ast: AstAtomAst::from_element(ChemElement::O),
                }),
                AstDelta::Atom(AstAtomDelta::Add {
                    id: AstAtomId(2),
                    ast: AstAtomAst::from_element(ChemElement::N),
                }),
            ]
            .into_iter()
            .collect(),
        ),
    )]
    #[case::bond_order(
        AstMoleculeAst::from_parts(AstMoleculeParts {
            atoms: vec![
                AstAtomAst::from_element(ChemElement::C),
                AstAtomAst::from_element(ChemElement::C),
            ],
            bonds: vec![(AstAtomId(0), AstAtomId(1), AstBondAst::from_order(1))],
            ..Default::default()
        }),
        AstMoleculeAst::from_parts(AstMoleculeParts {
            atoms: vec![
                AstAtomAst::from_element(ChemElement::C),
                AstAtomAst::from_element(ChemElement::C),
            ],
            bonds: vec![(AstAtomId(0), AstAtomId(1), AstBondAst::from_order(2))],
            ..Default::default()
        }),
        vec![(0, 0), (1, 1)],
        AstReactionAst::new(
            AstMoleculeAst::from_parts(AstMoleculeParts {
                atoms: vec![
                    AstAtomAst::from_element(ChemElement::C),
                    AstAtomAst::from_element(ChemElement::C),
                ],
                bonds: vec![(AstAtomId(0), AstAtomId(1), AstBondAst::from_order(1))],
                ..Default::default()
            }),
            vec![AstDelta::Bond(AstBondDelta::ModifyField {
                id: AstBondId(0),
                change: AstBondFieldChange::Order {
                    old: AstValueAst::Lit(1),
                    new: AstValueAst::Lit(2),
                },
            })]
            .into_iter()
            .collect(),
        ),
    )]
    fn test_reaction_ast_from_sides(
        #[case] lhs: AstMoleculeAst,
        #[case] rhs: AstMoleculeAst,
        #[case] atom_pairs: Vec<(usize, usize)>,
        #[case] expected: AstReactionAst,
    ) {
        Python::attach(|py| {
            let lhs_before = lhs.clone();
            let rhs_before = rhs.clone();
            let lhs = Py::new(py, MoleculeAst::from_inner(lhs)).unwrap();
            let rhs = Py::new(py, MoleculeAst::from_inner(rhs)).unwrap();

            let atom_pairs = PyList::new(py, atom_pairs).unwrap();
            let reaction = ReactionAst::from_sides(
                py,
                lhs.clone_ref(py),
                rhs.clone_ref(py),
                atom_pairs.as_any(),
            )
            .unwrap();

            assert_eq!(reaction.to_rust(py), expected);
            assert_eq!(*lhs.bind(py).borrow().inner(), lhs_before);
            assert_eq!(*rhs.bind(py).borrow().inner(), rhs_before);
            assert_ne!(reaction.lhs.as_ptr(), lhs.as_ptr());
        });
    }

    #[rstest]
    #[case::dative_bond(
        r#"{:atoms ["N" "B"] :bonds []}"#,
        r#"{:atoms ["N" "B"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type "1"}]}"#,
        vec![(0, 0), (1, 1)],
        vec![AstDelta::DativeBond(AstDativeBondDelta::Add {
            id: AstDativeBondId(0),
            donors: vec![AstAtomId(0)],
            acceptor: AstAtomId(1),
            ast: AstDativeBondAst::from_order(1),
        })],
    )]
    #[case::aromatic_system(
        r#"{:atoms ["C" "C"] :bonds []}"#,
        r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "[1,1]"}]}"#,
        vec![(0, 0), (1, 1)],
        vec![AstDelta::AromaticSystem(AstAromaticSystemDelta::Add {
            id: AstAromaticSystemId(0),
            atoms: vec![AstAtomId(0), AstAtomId(1)],
            ast: AstAromaticSystemAst::from_electrons(vec![1, 1]),
        })],
    )]
    #[case::multicenter_bond(
        r#"{:atoms ["B" "H" "B"] :bonds []}"#,
        r#"{:atoms ["B" "H" "B"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :type "[3,5,7]"}]}"#,
        vec![(0, 0), (1, 1), (2, 2)],
        vec![AstDelta::MulticenterBond(AstMulticenterBondDelta::Add {
            id: AstMulticenterBondId(0),
            atoms: vec![AstAtomId(0), AstAtomId(1), AstAtomId(2)],
            ast: AstMulticenterBondAst::from_electrons(vec![3, 5, 7]),
        })],
    )]
    #[case::noncovalent_bond(
        r#"{:atoms ["O" "O"] :bonds []}"#,
        r#"{:atoms ["O" "O"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}]}"#,
        vec![(0, 0), (1, 1)],
        vec![AstDelta::NoncovalentBond(AstNoncovalentBondDelta::Add {
            id: AstNoncovalentBondId(0),
            atoms: [AstAtomId(0), AstAtomId(1)],
            ast: AstNoncovalentBondAst::from_kind(AstNoncovalentBondKind::HydrogenBond),
        })],
    )]
    #[case::stereo_atom(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds []}"#,
        r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]}"#,
        vec![(0, 0), (1, 1), (2, 2), (3, 3), (4, 4)],
        vec![AstDelta::StereoAtom(AstStereoAtomDelta::Add {
            id: AstStereoAtomId(0),
            site: AstAtomId(0),
            ligands: vec![
                AstStereoLigand::new(AstAtomId(1), AstStereoLigandKind::Atom),
                AstStereoLigand::new(AstAtomId(2), AstStereoLigandKind::Atom),
                AstStereoLigand::new(AstAtomId(3), AstStereoLigandKind::Atom),
                AstStereoLigand::new(AstAtomId(4), AstStereoLigandKind::Atom),
            ],
            ast: AstStereoAtomAst::new(AstStereoKind::Tetrahedral, AstStereoCosetAst::Lit(1)),
        })],
    )]
    #[case::stereo_bond(
        r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]}"#,
        r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3] :type "Ct1"}]}"#,
        vec![(0, 0), (1, 1), (2, 2), (3, 3)],
        vec![AstDelta::StereoBond(AstStereoBondDelta::Add {
            id: AstStereoBondId(0),
            site: AstBondId(1),
            ligands: vec![
                AstStereoLigand::new(AstAtomId(0), AstStereoLigandKind::Atom),
                AstStereoLigand::new(AstAtomId(3), AstStereoLigandKind::Atom),
            ],
            ast: AstStereoBondAst::new(AstStereoKind::CisTrans, AstStereoCosetAst::Lit(1)),
        })],
    )]
    #[case::molecule_constraint(
        r#"{:atoms ["C"] :bonds []}"#,
        r#"{:atoms ["C"] :bonds [] :constraints [{:connected {}}]}"#,
        vec![(0, 0)],
        vec![AstDelta::Constraint(AstConstraintDelta::Add(
            AstConstraint::Molecule(AstMoleculeConstraint::Connected { atoms: None }),
        ))],
    )]
    fn test_reaction_ast_from_sides_entities(
        #[case] lhs: &str,
        #[case] rhs: &str,
        #[case] atom_pairs: Vec<(usize, usize)>,
        #[case] expected_deltas: Vec<AstDelta>,
    ) {
        Python::attach(|py| {
            let lhs = lhs.parse::<AstMoleculeAst>().unwrap();
            let rhs = rhs.parse::<AstMoleculeAst>().unwrap();
            let atom_pairs = PyList::new(py, atom_pairs).unwrap();
            let reaction = ReactionAst::from_sides(
                py,
                Py::new(py, MoleculeAst::from_inner(lhs.clone())).unwrap(),
                Py::new(py, MoleculeAst::from_inner(rhs)).unwrap(),
                atom_pairs.as_any(),
            )
            .unwrap();

            assert_eq!(
                reaction.to_rust(py),
                AstReactionAst::new(lhs, expected_deltas.into_iter().collect())
            );
        });
    }

    #[rstest]
    fn test_reaction_ast_from_sides_snapshot() {
        Python::attach(|py| {
            let lhs_before = AstMoleculeAst::from_parts(AstMoleculeParts {
                atoms: vec![
                    AstAtomAst::from_element(ChemElement::C),
                    AstAtomAst::from_element(ChemElement::O),
                ],
                ..Default::default()
            });
            let rhs_before = AstMoleculeAst::from_parts(AstMoleculeParts {
                atoms: vec![
                    AstAtomAst::from_element(ChemElement::C),
                    AstAtomAst::from_element(ChemElement::N),
                ],
                ..Default::default()
            });
            let lhs = Py::new(py, MoleculeAst::from_inner(lhs_before.clone())).unwrap();
            let rhs = Py::new(py, MoleculeAst::from_inner(rhs_before.clone())).unwrap();
            let atom_pairs = PyList::new(py, [(0, 0)]).unwrap();
            let reaction = ReactionAst::from_sides(
                py,
                lhs.clone_ref(py),
                rhs.clone_ref(py),
                atom_pairs.as_any(),
            )
            .unwrap();
            let expected = reaction.to_rust(py);

            *lhs.bind(py).borrow_mut().inner_mut() = AstMoleculeAst::new();
            *rhs.bind(py).borrow_mut().inner_mut() = AstMoleculeAst::new();

            assert_eq!(reaction.to_rust(py), expected);
            assert_ne!(reaction.lhs.as_ptr(), lhs.as_ptr());

            *reaction.lhs.bind(py).borrow_mut().inner_mut() =
                AstMoleculeAst::from_parts(AstMoleculeParts {
                    atoms: vec![AstAtomAst::from_element(ChemElement::F)],
                    ..Default::default()
                });
            let delta = into_py_variant(
                py,
                Delta::from_rust(
                    py,
                    &AstDelta::Atom(AstAtomDelta::Add {
                        id: AstAtomId(3),
                        ast: AstAtomAst::from_element(ChemElement::Cl),
                    }),
                )
                .unwrap(),
            )
            .unwrap();
            reaction
                .deltas
                .bind(py)
                .call_method1("append", (delta,))
                .unwrap();
            let changed = reaction.to_rust(py);

            assert_eq!(
                changed.lhs,
                AstMoleculeAst::from_parts(AstMoleculeParts {
                    atoms: vec![AstAtomAst::from_element(ChemElement::F)],
                    ..Default::default()
                })
            );
            assert_eq!(
                changed.deltas.as_slice().last(),
                Some(&AstDelta::Atom(AstAtomDelta::Add {
                    id: AstAtomId(3),
                    ast: AstAtomAst::from_element(ChemElement::Cl),
                }))
            );
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
    fn test_reaction_ast_canonicalize() {
        Python::attach(|py| {
            let source = ReactionAst::from_rust(
                py,
                AstReactionAst::new(
                    AstMoleculeAst::from_parts(AstMoleculeParts {
                        atoms: vec![AstAtomAst::from_element(ChemElement::C).with_charge(0)],
                        ..Default::default()
                    }),
                    vec![
                        AstDelta::Atom(AstAtomDelta::ModifyField {
                            id: AstAtomId(0),
                            change: AstAtomFieldChange::Charge {
                                old: AstValueAst::Lit(0),
                                new: AstValueAst::Lit(1),
                            },
                        }),
                        AstDelta::Atom(AstAtomDelta::ModifyField {
                            id: AstAtomId(0),
                            change: AstAtomFieldChange::Charge {
                                old: AstValueAst::Lit(1),
                                new: AstValueAst::Lit(2),
                            },
                        }),
                    ]
                    .into_iter()
                    .collect(),
                ),
            )
            .unwrap();
            let before = source.to_rust(py);
            let expected = AstReactionAst::new(
                before.lhs.clone(),
                vec![AstDelta::Atom(AstAtomDelta::ModifyField {
                    id: AstAtomId(0),
                    change: AstAtomFieldChange::Charge {
                        old: AstValueAst::Lit(0),
                        new: AstValueAst::Lit(2),
                    },
                })]
                .into_iter()
                .collect(),
            );

            let canonical = source.canonicalize(py).unwrap();
            let twice = canonical.canonicalize(py).unwrap();

            assert_eq!(canonical.to_rust(py), expected);
            assert_eq!(twice.to_rust(py), expected);
            assert_eq!(source.to_rust(py), before);
            assert_ne!(canonical.lhs.as_ptr(), source.lhs.as_ptr());
            assert_ne!(canonical.deltas.as_ptr(), source.deltas.as_ptr());
        });
    }

    #[rstest]
    fn test_reaction_ast_canonicalize_error() {
        Python::attach(|py| {
            let source = ReactionAst::from_rust(
                py,
                AstReactionAst::new(
                    AstMoleculeAst::from_parts(AstMoleculeParts {
                        atoms: vec![AstAtomAst::from_element(ChemElement::C).with_charge(0)],
                        ..Default::default()
                    }),
                    vec![
                        AstDelta::Atom(AstAtomDelta::ModifyField {
                            id: AstAtomId(0),
                            change: AstAtomFieldChange::Charge {
                                old: AstValueAst::Lit(0),
                                new: AstValueAst::Lit(1),
                            },
                        }),
                        AstDelta::Atom(AstAtomDelta::ModifyField {
                            id: AstAtomId(0),
                            change: AstAtomFieldChange::Charge {
                                old: AstValueAst::Lit(2),
                                new: AstValueAst::Lit(3),
                            },
                        }),
                    ]
                    .into_iter()
                    .collect(),
                ),
            )
            .unwrap();
            let before = source.to_rust(py);

            let error = source.canonicalize(py).err().unwrap();

            assert!(error.is_instance_of::<ContradictionError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "reached a contradiction"
            );
            assert_eq!(source.to_rust(py), before);
        });
    }

    #[rstest]
    fn test_reaction_ast_reverse() {
        Python::attach(|py| {
            let source = ReactionAst::parse(
                py,
                r##"{:lhs {:atoms ["C" "O"]} :deltas [{:atom {:add "N"}} {:atom {:remove 1}}]}"##,
            )
            .unwrap();
            let before = source.to_rust(py);
            let expected_roundtrip = before.clone().canonicalize().unwrap();

            let reversed = source.reverse(py).unwrap();
            let roundtrip = reversed.reverse(py).unwrap();

            assert_eq!(
                reversed.to_rust(py).lhs,
                AstMoleculeAst::from_parts(AstMoleculeParts {
                    atoms: vec![
                        AstAtomAst::from_element(ChemElement::C),
                        AstAtomAst::from_element(ChemElement::N),
                    ],
                    ..Default::default()
                })
            );
            assert_eq!(
                roundtrip.to_rust(py).canonicalize().unwrap(),
                expected_roundtrip
            );
            assert_eq!(source.to_rust(py), before);
            assert_ne!(reversed.lhs.as_ptr(), source.lhs.as_ptr());
            assert_ne!(reversed.deltas.as_ptr(), source.deltas.as_ptr());
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
    #[case::empty(
        AstReactionAst::default(),
        r##"{:deltas [] :lhs {:atoms [] :bonds []}}"##
    )]
    #[case::populated(
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
        r##"{:deltas [{:atom {:add "O"}}] :lhs {:atoms ["C"] :bonds []}}"##,
    )]
    fn test_reaction_ast_str(#[case] input: AstReactionAst, #[case] expected: &str) {
        Python::attach(|py| {
            let reaction = ReactionAst::from_rust(py, input).unwrap();

            assert_eq!(reaction.__str__(py), expected);
        });
    }

    #[rstest]
    fn test_reaction_ast_str_components() {
        Python::attach(|py| {
            let reaction = ReactionAst::from_rust(
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

            *reaction.lhs.bind(py).borrow_mut().inner_mut() =
                AstMoleculeAst::from_parts(AstMoleculeParts {
                    atoms: vec![AstAtomAst::from_element(ChemElement::C).with_charge(1)],
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
            reaction
                .deltas
                .bind(py)
                .call_method1("append", (delta,))
                .unwrap();

            assert_eq!(
                reaction.__str__(py),
                r##"{:deltas [{:atom {:add "O"}}] :lhs {:atoms ["C#c+"] :bonds []}}"##
            );
        });
    }

    #[rstest]
    #[case::atom_add_remove(
        r##"{:lhs {:atoms ["C" "O"]} :deltas [{:atom {:add "N"}} {:atom {:remove 1}}]}"##
    )]
    #[case::atom_modify(r##"{:lhs {:atoms ["Br#c0"]} :deltas [{:atom {:modify [0 "#c-1"]}}]}"##)]
    #[case::stereo_mirror(
        r##"{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]} :deltas [{:stereo-atom {:mirror [0 :tetrahedral]}}]}"##
    )]
    #[case::molecule_constraint(
        r##"{:lhs {:atoms ["C"]} :deltas [{:constraint {:add {:connected {}}}}]}"##
    )]
    fn test_reaction_ast_str_roundtrip(#[case] text: &str) {
        Python::attach(|py| {
            let first = ReactionAst::parse(py, text).unwrap();

            let canonical = first.__str__(py);
            let second = ReactionAst::parse(py, &canonical).unwrap();

            assert!(first.__eq__(&second, py));
            assert_eq!(second.__str__(py), canonical);
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

    #[rstest]
    fn test_reaction_ast_to_rust_roundtrip() {
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
            let python =
                Py::new(py, ReactionAst::from_rust(py, expected.clone()).unwrap()).unwrap();

            let rust = python.bind(py).borrow().to_rust(py);
            let roundtrip = Py::new(py, ReactionAst::from_rust(py, rust).unwrap()).unwrap();

            assert_eq!(roundtrip.bind(py).borrow().to_rust(py), expected);
            assert_ne!(python.as_ptr(), roundtrip.as_ptr());
            assert_ne!(
                python.bind(py).borrow().lhs.as_ptr(),
                roundtrip.bind(py).borrow().lhs.as_ptr()
            );
            assert_ne!(
                python.bind(py).borrow().deltas.as_ptr(),
                roundtrip.bind(py).borrow().deltas.as_ptr()
            );
        });
    }
}
