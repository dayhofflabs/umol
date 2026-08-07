//! `ReactionSpanAst` — a Python facade over the superimposed reaction span AST.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use umol_ast::ast::{
    AtomId as AstAtomId, BondId as AstBondId, Canonicalize, Constraint as AstConstraint,
    ConstraintSpan as AstConstraintSpan, EntitySpan as AstEntitySpan,
    ReactionSpanAst as AstReactionSpanAst, ReactionSpanEntries as AstReactionSpanEntries,
};

use crate::aromatic::AromaticSystemAst;
use crate::atom::AtomAst;
use crate::bond::BondAst;
use crate::constraint::molecule::Constraint;
use crate::dative::DativeBondAst;
use crate::multicenter::MulticenterBondAst;
use crate::noncovalent::NoncovalentBondAst;
use crate::stereo::{StereoAtomAst, StereoBondAst, StereoLigand};

type SpanPair<T> = (Option<Py<T>>, Option<Py<T>>);

fn entity_span<T: Canonicalize>(lhs: Option<T>, rhs: Option<T>) -> PyResult<AstEntitySpan<T>> {
    AstEntitySpan::superimpose(lhs, rhs)
        .ok_or_else(|| PyValueError::new_err("reaction span entry is absent from both sides"))
}

fn constraint_spans(
    lhs: Option<AstConstraint>,
    rhs: Option<AstConstraint>,
) -> PyResult<Vec<AstConstraintSpan>> {
    match (lhs, rhs) {
        (Some(lhs), Some(rhs)) if lhs.canonical_eq(&rhs) => {
            Ok(vec![AstConstraintSpan::Unchanged(lhs)])
        }
        (Some(lhs), Some(rhs)) => Ok(vec![
            AstConstraintSpan::Removed(lhs),
            AstConstraintSpan::Added(rhs),
        ]),
        (Some(lhs), None) => Ok(vec![AstConstraintSpan::Removed(lhs)]),
        (None, Some(rhs)) => Ok(vec![AstConstraintSpan::Added(rhs)]),
        (None, None) => Err(PyValueError::new_err(
            "reaction span entry is absent from both sides",
        )),
    }
}

/// A superimposed reaction span with explicit before/after entity states.
#[pyclass(eq, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReactionSpanAst(AstReactionSpanAst);

#[pymethods]
impl ReactionSpanAst {
    /// Construct a reaction span from union-frame entries.
    ///
    /// Every entity value is a `(lhs, rhs)` pair. Either member may be `None`, but not both.
    /// Construction checks that all union-frame references resolve and that the selected entries
    /// on each side form a structurally intact molecule. Chemistry is not validated.
    #[staticmethod]
    #[pyo3(signature = (atoms, *, bonds=Vec::new(), dative_bonds=Vec::new(), aromatic_systems=Vec::new(), multicenter_bonds=Vec::new(), noncovalent_bonds=Vec::new(), stereo_atoms=Vec::new(), stereo_bonds=Vec::new(), constraints=Vec::new()))]
    #[allow(clippy::too_many_arguments)] // one argument per entity family
    fn from_entries(
        py: Python<'_>,
        atoms: Vec<SpanPair<AtomAst>>,
        bonds: Vec<(u32, u32, SpanPair<BondAst>)>,
        dative_bonds: Vec<(Vec<u32>, u32, SpanPair<DativeBondAst>)>,
        aromatic_systems: Vec<(Vec<u32>, SpanPair<AromaticSystemAst>)>,
        multicenter_bonds: Vec<(Vec<u32>, SpanPair<MulticenterBondAst>)>,
        noncovalent_bonds: Vec<([u32; 2], SpanPair<NoncovalentBondAst>)>,
        stereo_atoms: Vec<(u32, Vec<StereoLigand>, SpanPair<StereoAtomAst>)>,
        stereo_bonds: Vec<(u32, Vec<StereoLigand>, SpanPair<StereoBondAst>)>,
        constraints: Vec<SpanPair<Constraint>>,
    ) -> PyResult<Self> {
        let atoms = atoms
            .into_iter()
            .map(|(lhs, rhs)| {
                entity_span(
                    lhs.map(|value| value.bind(py).borrow().inner().clone()),
                    rhs.map(|value| value.bind(py).borrow().inner().clone()),
                )
            })
            .collect::<PyResult<Vec<_>>>()?;
        let bonds = bonds
            .into_iter()
            .map(|(first, second, (lhs, rhs))| {
                Ok((
                    AstAtomId(first),
                    AstAtomId(second),
                    entity_span(
                        lhs.map(|value| value.bind(py).borrow().inner().clone()),
                        rhs.map(|value| value.bind(py).borrow().inner().clone()),
                    )?,
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        let dative = dative_bonds
            .into_iter()
            .map(|(donors, acceptor, (lhs, rhs))| {
                Ok((
                    donors.into_iter().map(AstAtomId).collect(),
                    AstAtomId(acceptor),
                    entity_span(
                        lhs.map(|value| value.bind(py).borrow().inner().clone()),
                        rhs.map(|value| value.bind(py).borrow().inner().clone()),
                    )?,
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        let aromatic = aromatic_systems
            .into_iter()
            .map(|(atoms, (lhs, rhs))| {
                Ok((
                    atoms.into_iter().map(AstAtomId).collect(),
                    entity_span(
                        lhs.map(|value| value.bind(py).borrow().inner().clone()),
                        rhs.map(|value| value.bind(py).borrow().inner().clone()),
                    )?,
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        let multicenter = multicenter_bonds
            .into_iter()
            .map(|(atoms, (lhs, rhs))| {
                Ok((
                    atoms.into_iter().map(AstAtomId).collect(),
                    entity_span(
                        lhs.map(|value| value.bind(py).borrow().inner().clone()),
                        rhs.map(|value| value.bind(py).borrow().inner().clone()),
                    )?,
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        let noncovalent = noncovalent_bonds
            .into_iter()
            .map(|([first, second], (lhs, rhs))| {
                Ok((
                    AstAtomId(first),
                    AstAtomId(second),
                    entity_span(
                        lhs.map(|value| value.bind(py).borrow().inner().clone()),
                        rhs.map(|value| value.bind(py).borrow().inner().clone()),
                    )?,
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        let stereo_atoms = stereo_atoms
            .into_iter()
            .map(|(site, ligands, (lhs, rhs))| {
                Ok((
                    AstAtomId(site),
                    ligands.into_iter().map(StereoLigand::to_rust).collect(),
                    entity_span(
                        lhs.map(|value| value.bind(py).borrow().inner().clone()),
                        rhs.map(|value| value.bind(py).borrow().inner().clone()),
                    )?,
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        let stereo_bonds = stereo_bonds
            .into_iter()
            .map(|(site, ligands, (lhs, rhs))| {
                Ok((
                    AstBondId(site),
                    ligands.into_iter().map(StereoLigand::to_rust).collect(),
                    entity_span(
                        lhs.map(|value| value.bind(py).borrow().inner().clone()),
                        rhs.map(|value| value.bind(py).borrow().inner().clone()),
                    )?,
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        let mut constraint_entries = Vec::new();
        for (lhs, rhs) in constraints {
            constraint_entries.extend(constraint_spans(
                lhs.map(|value| value.bind(py).borrow().to_rust(py)),
                rhs.map(|value| value.bind(py).borrow().to_rust(py)),
            )?);
        }

        AstReactionSpanAst::try_from_entries(AstReactionSpanEntries {
            atoms,
            bonds,
            dative,
            aromatic,
            multicenter,
            noncovalent,
            stereo_atoms,
            stereo_bonds,
            constraints: constraint_entries,
        })
        .map(Self)
        .map_err(|error| PyValueError::new_err(error.to_string()))
    }
}

impl ReactionSpanAst {
    #[cfg(test)]
    pub(crate) fn to_rust(&self) -> AstReactionSpanAst {
        self.0.clone()
    }
}

#[cfg(test)]
mod tests {
    use pyo3::exceptions::PyValueError;
    use rstest::rstest;
    use umol_ast::ast::{
        AromaticSystemAst as AstAromaticSystemAst, AtomAst as AstAtomAst, BondAst as AstBondAst,
        Constraint as AstConstraint, DativeBondAst as AstDativeBondAst,
        MoleculeConstraint as AstMoleculeConstraint, MulticenterBondAst as AstMulticenterBondAst,
        NoncovalentBondAst as AstNoncovalentBondAst, NoncovalentBondKind as AstNoncovalentBondKind,
        StereoAtomAst as AstStereoAtomAst, StereoBondAst as AstStereoBondAst,
        StereoCoset as AstStereoCoset, StereoKind as AstStereoKind,
        StereoLigand as AstStereoLigand, StereoLigandKind as AstStereoLigandKind,
        ValueAst as AstValueAst,
    };
    use umol_chem::element::Element as ChemElement;

    use super::*;
    use crate::convert::into_py_variant;

    #[rstest]
    fn test_reaction_span_ast_from_entries() {
        Python::attach(|py| {
            let canonical_lhs =
                AstAtomAst::from_element(ChemElement::C).with_charge(AstValueAst::Lit(1));
            let canonical_rhs =
                AstAtomAst::from_element(ChemElement::C).with_charge(AstValueAst::lit_set([1_i64]));
            let modified_lhs = AstAtomAst::from_element(ChemElement::C);
            let modified_rhs = AstAtomAst::from_element(ChemElement::N);
            let removed_atom = AstAtomAst::from_element(ChemElement::O);
            let added_atom = AstAtomAst::from_element(ChemElement::F);
            let unchanged_bond = AstBondAst::from_order(1);
            let dative_lhs = AstDativeBondAst::from_order(1);
            let dative_rhs = AstDativeBondAst::from_order(2);
            let added_aromatic = AstAromaticSystemAst::from_electrons(vec![1, 1]);
            let removed_multicenter = AstMulticenterBondAst::from_electrons(vec![1, 1]);
            let unchanged_noncovalent =
                AstNoncovalentBondAst::from_kind(AstNoncovalentBondKind::HydrogenBond);
            let stereo_atom_lhs =
                AstStereoAtomAst::new(AstStereoKind::Tetrahedral, AstStereoCoset::Lit(1));
            let stereo_atom_rhs =
                AstStereoAtomAst::new(AstStereoKind::Tetrahedral, AstStereoCoset::Lit(2));
            let added_stereo_bond =
                AstStereoBondAst::new(AstStereoKind::CisTrans, AstStereoCoset::Lit(1));
            let ligand = AstStereoLigand::new(AstAtomId(1), AstStereoLigandKind::Atom);
            let unchanged_constraint =
                AstConstraint::Molecule(AstMoleculeConstraint::Connected { atoms: None });
            let modified_constraint_lhs =
                AstConstraint::Molecule(AstMoleculeConstraint::Connected {
                    atoms: Some(vec![AstAtomId(0)]),
                });
            let modified_constraint_rhs =
                AstConstraint::Molecule(AstMoleculeConstraint::Connected {
                    atoms: Some(vec![AstAtomId(1)]),
                });
            let removed_constraint = AstConstraint::Molecule(AstMoleculeConstraint::Connected {
                atoms: Some(vec![AstAtomId(2)]),
            });
            let added_constraint = AstConstraint::Molecule(AstMoleculeConstraint::Connected {
                atoms: Some(vec![AstAtomId(3)]),
            });

            let span = ReactionSpanAst::from_entries(
                py,
                vec![
                    (
                        Some(Py::new(py, AtomAst::from_inner(canonical_lhs.clone())).unwrap()),
                        Some(Py::new(py, AtomAst::from_inner(canonical_rhs)).unwrap()),
                    ),
                    (
                        Some(Py::new(py, AtomAst::from_inner(modified_lhs.clone())).unwrap()),
                        Some(Py::new(py, AtomAst::from_inner(modified_rhs.clone())).unwrap()),
                    ),
                    (
                        Some(Py::new(py, AtomAst::from_inner(removed_atom.clone())).unwrap()),
                        None,
                    ),
                    (
                        None,
                        Some(Py::new(py, AtomAst::from_inner(added_atom.clone())).unwrap()),
                    ),
                ],
                vec![(
                    0,
                    1,
                    (
                        Some(Py::new(py, BondAst::from_inner(unchanged_bond.clone())).unwrap()),
                        Some(Py::new(py, BondAst::from_inner(unchanged_bond.clone())).unwrap()),
                    ),
                )],
                vec![(
                    vec![1],
                    0,
                    (
                        Some(Py::new(py, DativeBondAst::from_inner(dative_lhs.clone())).unwrap()),
                        Some(Py::new(py, DativeBondAst::from_inner(dative_rhs.clone())).unwrap()),
                    ),
                )],
                vec![(
                    vec![0, 1],
                    (
                        None,
                        Some(
                            Py::new(py, AromaticSystemAst::from_inner(added_aromatic.clone()))
                                .unwrap(),
                        ),
                    ),
                )],
                vec![(
                    vec![0, 1],
                    (
                        Some(
                            Py::new(
                                py,
                                MulticenterBondAst::from_inner(removed_multicenter.clone()),
                            )
                            .unwrap(),
                        ),
                        None,
                    ),
                )],
                vec![(
                    [0, 1],
                    (
                        Some(
                            Py::new(
                                py,
                                NoncovalentBondAst::from_inner(unchanged_noncovalent.clone()),
                            )
                            .unwrap(),
                        ),
                        Some(
                            Py::new(
                                py,
                                NoncovalentBondAst::from_inner(unchanged_noncovalent.clone()),
                            )
                            .unwrap(),
                        ),
                    ),
                )],
                vec![(
                    0,
                    vec![StereoLigand::from_rust(ligand)],
                    (
                        Some(
                            Py::new(py, StereoAtomAst::from_inner(stereo_atom_lhs.clone()))
                                .unwrap(),
                        ),
                        Some(
                            Py::new(py, StereoAtomAst::from_inner(stereo_atom_rhs.clone()))
                                .unwrap(),
                        ),
                    ),
                )],
                vec![(
                    0,
                    vec![StereoLigand::from_rust(ligand)],
                    (
                        None,
                        Some(
                            Py::new(py, StereoBondAst::from_inner(added_stereo_bond.clone()))
                                .unwrap(),
                        ),
                    ),
                )],
                vec![
                    (
                        Some(
                            into_py_variant(
                                py,
                                Constraint::from_rust(py, &unchanged_constraint).unwrap(),
                            )
                            .unwrap(),
                        ),
                        Some(
                            into_py_variant(
                                py,
                                Constraint::from_rust(py, &unchanged_constraint).unwrap(),
                            )
                            .unwrap(),
                        ),
                    ),
                    (
                        Some(
                            into_py_variant(
                                py,
                                Constraint::from_rust(py, &modified_constraint_lhs).unwrap(),
                            )
                            .unwrap(),
                        ),
                        Some(
                            into_py_variant(
                                py,
                                Constraint::from_rust(py, &modified_constraint_rhs).unwrap(),
                            )
                            .unwrap(),
                        ),
                    ),
                    (
                        Some(
                            into_py_variant(
                                py,
                                Constraint::from_rust(py, &removed_constraint).unwrap(),
                            )
                            .unwrap(),
                        ),
                        None,
                    ),
                    (
                        None,
                        Some(
                            into_py_variant(
                                py,
                                Constraint::from_rust(py, &added_constraint).unwrap(),
                            )
                            .unwrap(),
                        ),
                    ),
                ],
            )
            .unwrap();

            assert_eq!(
                span.to_rust(),
                AstReactionSpanAst::from_entries(AstReactionSpanEntries {
                    atoms: vec![
                        AstEntitySpan::Unchanged(canonical_lhs),
                        AstEntitySpan::Modified {
                            lhs: modified_lhs,
                            rhs: modified_rhs,
                        },
                        AstEntitySpan::Removed(removed_atom),
                        AstEntitySpan::Added(added_atom),
                    ],
                    bonds: vec![(
                        AstAtomId(0),
                        AstAtomId(1),
                        AstEntitySpan::Unchanged(unchanged_bond),
                    )],
                    dative: vec![(
                        vec![AstAtomId(1)],
                        AstAtomId(0),
                        AstEntitySpan::Modified {
                            lhs: dative_lhs,
                            rhs: dative_rhs,
                        },
                    )],
                    aromatic: vec![(
                        vec![AstAtomId(0), AstAtomId(1)],
                        AstEntitySpan::Added(added_aromatic),
                    )],
                    multicenter: vec![(
                        vec![AstAtomId(0), AstAtomId(1)],
                        AstEntitySpan::Removed(removed_multicenter),
                    )],
                    noncovalent: vec![(
                        AstAtomId(0),
                        AstAtomId(1),
                        AstEntitySpan::Unchanged(unchanged_noncovalent),
                    )],
                    stereo_atoms: vec![(
                        AstAtomId(0),
                        vec![ligand],
                        AstEntitySpan::Modified {
                            lhs: stereo_atom_lhs,
                            rhs: stereo_atom_rhs,
                        },
                    )],
                    stereo_bonds: vec![(
                        AstBondId(0),
                        vec![ligand],
                        AstEntitySpan::Added(added_stereo_bond),
                    )],
                    constraints: vec![
                        AstConstraintSpan::Unchanged(unchanged_constraint),
                        AstConstraintSpan::Removed(modified_constraint_lhs),
                        AstConstraintSpan::Added(modified_constraint_rhs),
                        AstConstraintSpan::Removed(removed_constraint),
                        AstConstraintSpan::Added(added_constraint),
                    ],
                })
            );
        });
    }

    #[rstest]
    fn test_reaction_span_ast_from_entries_error() {
        Python::attach(|py| {
            let error = ReactionSpanAst::from_entries(
                py,
                vec![(None, None)],
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
            )
            .unwrap_err();

            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "reaction span entry is absent from both sides"
            );
        });
    }

    #[rstest]
    #[case::union(true, true, false)]
    #[case::lhs(false, true, true)]
    #[case::rhs(true, false, true)]
    fn test_reaction_span_ast_from_entries_reference_error(
        #[case] first_on_lhs: bool,
        #[case] first_on_rhs: bool,
        #[case] include_second: bool,
    ) {
        Python::attach(|py| {
            let mut atoms = vec![(
                first_on_lhs.then(|| {
                    Py::new(
                        py,
                        AtomAst::from_inner(AstAtomAst::from_element(ChemElement::C)),
                    )
                    .unwrap()
                }),
                first_on_rhs.then(|| {
                    Py::new(
                        py,
                        AtomAst::from_inner(AstAtomAst::from_element(ChemElement::C)),
                    )
                    .unwrap()
                }),
            )];
            if include_second {
                atoms.push((
                    Some(
                        Py::new(
                            py,
                            AtomAst::from_inner(AstAtomAst::from_element(ChemElement::O)),
                        )
                        .unwrap(),
                    ),
                    Some(
                        Py::new(
                            py,
                            AtomAst::from_inner(AstAtomAst::from_element(ChemElement::O)),
                        )
                        .unwrap(),
                    ),
                ));
            }
            let error = ReactionSpanAst::from_entries(
                py,
                atoms,
                vec![(
                    0,
                    1,
                    (
                        Some(Py::new(py, BondAst::from_inner(AstBondAst::from_order(1))).unwrap()),
                        Some(Py::new(py, BondAst::from_inner(AstBondAst::from_order(1))).unwrap()),
                    ),
                )],
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
            )
            .unwrap_err();

            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "reaction span entries reference unavailable atom 1"
            );
        });
    }
}
