//! `ReactionSpanAst` — a Python facade over the superimposed reaction span AST.

use std::str::FromStr;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use umol_graph_ir::dsl::ReactionSpanDsl as GraphIrReactionSpanDsl;
use umol_graph_ir::ir::{
    AtomId as GraphIrAtomId, BondId as GraphIrBondId, Canonicalize,
    Constraint as GraphIrConstraint, ConstraintSpan as GraphIrConstraintSpan,
    EntitySpan as GraphIrEntitySpan, FromIr, IntoIr, ReactionSpanAst as GraphIrReactionSpanAst,
    ReactionSpanEntries as GraphIrReactionSpanEntries,
};

use crate::aromatic::AromaticSystemAst;
use crate::atom::AtomAst;
use crate::bond::BondAst;
use crate::constraint::molecule::Constraint;
use crate::correspondence::MoleculeCorrespondence;
use crate::dative::DativeBondAst;
use crate::defaults::MoleculeDefaults;
use crate::error::{metadata_error, parse_error};
use crate::metadata::MoleculeMetadata;
use crate::molecule::MoleculeAst;
use crate::multicenter::MulticenterBondAst;
use crate::noncovalent::NoncovalentBondAst;
use crate::reaction::ReactionAst;
use crate::stereo::{StereoAtomAst, StereoBondAst, StereoLigand};

type SpanPair<T> = (Option<Py<T>>, Option<Py<T>>);

fn entity_span<T: Canonicalize>(lhs: Option<T>, rhs: Option<T>) -> PyResult<GraphIrEntitySpan<T>> {
    GraphIrEntitySpan::superimpose(lhs, rhs)
        .ok_or_else(|| PyValueError::new_err("reaction span entry is absent from both sides"))
}

fn constraint_spans(
    lhs: Option<GraphIrConstraint>,
    rhs: Option<GraphIrConstraint>,
) -> PyResult<Vec<GraphIrConstraintSpan>> {
    match (lhs, rhs) {
        (Some(lhs), Some(rhs)) if lhs.canonical_eq(&rhs) => {
            Ok(vec![GraphIrConstraintSpan::Unchanged(lhs)])
        }
        (Some(lhs), Some(rhs)) => Ok(vec![
            GraphIrConstraintSpan::Removed(lhs),
            GraphIrConstraintSpan::Added(rhs),
        ]),
        (Some(lhs), None) => Ok(vec![GraphIrConstraintSpan::Removed(lhs)]),
        (None, Some(rhs)) => Ok(vec![GraphIrConstraintSpan::Added(rhs)]),
        (None, None) => Err(PyValueError::new_err(
            "reaction span entry is absent from both sides",
        )),
    }
}

/// A superimposed reaction span with explicit before/after entity states.
#[pyclass(eq, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReactionSpanAst(GraphIrReactionSpanAst);

#[pymethods]
impl ReactionSpanAst {
    /// Parse a reaction span from its EDN representation under explicit construction defaults.
    #[staticmethod]
    #[pyo3(signature = (text, *, defaults=None))]
    fn parse(text: &str, defaults: Option<MoleculeDefaults>) -> PyResult<Self> {
        let defaults = defaults.unwrap_or_else(MoleculeDefaults::new).to_rust();
        let span = GraphIrReactionSpanDsl::from_str(text)
            .map_err(parse_error)?
            .into_ir(&defaults);
        Ok(Self::from_rust(span))
    }

    /// Parse a reaction span and return `(span, metadata)`, retaining entity
    /// keywords and atom aliases for metadata-preserving rendering.
    #[staticmethod]
    #[pyo3(signature = (text, *, defaults=None))]
    fn parse_with_metadata(
        text: &str,
        defaults: Option<MoleculeDefaults>,
    ) -> PyResult<(Self, MoleculeMetadata)> {
        let defaults = defaults.unwrap_or_else(MoleculeDefaults::new).to_rust();
        let dsl = GraphIrReactionSpanDsl::from_str(text).map_err(parse_error)?;
        let metadata = MoleculeMetadata::from_rust(dsl.metadata().clone());
        Ok((Self::from_rust(dsl.into_ir(&defaults)), metadata))
    }

    /// Render a canonical positional DSL representation without entity
    /// keywords or atom aliases.
    #[pyo3(signature = (*, defaults=None))]
    fn render(&self, defaults: Option<MoleculeDefaults>) -> String {
        let defaults = defaults.unwrap_or_else(MoleculeDefaults::new).to_rust();
        GraphIrReactionSpanDsl::from_ir(&self.0, &defaults).to_string()
    }

    /// Render a canonical DSL representation with persistent metadata.
    ///
    /// Raises `MetadataError` if the detached metadata is not coherent with
    /// this reaction span.
    #[pyo3(signature = (metadata, *, defaults=None))]
    fn render_with_metadata(
        &self,
        metadata: &MoleculeMetadata,
        defaults: Option<MoleculeDefaults>,
    ) -> PyResult<String> {
        let defaults = defaults.unwrap_or_else(MoleculeDefaults::new).to_rust();
        let lowered = GraphIrReactionSpanDsl::from_ir(&self.0, &defaults)
            .into_parts()
            .0;
        GraphIrReactionSpanDsl::new(lowered, metadata.to_rust())
            .map(|dsl| dsl.to_string())
            .map_err(metadata_error)
    }

    fn __str__(&self) -> String {
        self.render(None)
    }

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
                    GraphIrAtomId(first),
                    GraphIrAtomId(second),
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
                    donors.into_iter().map(GraphIrAtomId).collect(),
                    GraphIrAtomId(acceptor),
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
                    atoms.into_iter().map(GraphIrAtomId).collect(),
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
                    atoms.into_iter().map(GraphIrAtomId).collect(),
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
                    GraphIrAtomId(first),
                    GraphIrAtomId(second),
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
                    GraphIrAtomId(site),
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
                    GraphIrBondId(site),
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

        GraphIrReactionSpanAst::try_from_entries(GraphIrReactionSpanEntries {
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
        .map(Self::from_rust)
        .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// Project the left-hand molecule as a detached snapshot.
    fn lhs(&self) -> MoleculeAst {
        MoleculeAst::from_rust(self.0.lhs())
    }

    /// Project the right-hand molecule as a detached snapshot.
    fn rhs(&self) -> MoleculeAst {
        MoleculeAst::from_rust(self.0.rhs())
    }

    /// Recover the correspondence between the normalized side projections.
    fn correspondence(&self) -> MoleculeCorrespondence {
        MoleculeCorrespondence::from_rust(self.0.correspondence())
    }

    /// Recover the reaction rule represented by this span.
    fn to_reaction(&self, py: Python<'_>) -> PyResult<ReactionAst> {
        ReactionAst::from_rust(py, self.0.to_reaction())
    }
}

impl ReactionSpanAst {
    pub(crate) fn from_rust(span: GraphIrReactionSpanAst) -> Self {
        Self(span)
    }

    #[cfg(test)]
    pub(crate) fn to_rust(&self) -> GraphIrReactionSpanAst {
        self.0.clone()
    }
}

#[cfg(test)]
mod tests {
    use pyo3::exceptions::PyValueError;
    use rstest::rstest;
    use umol_chem::element::Element as ChemElement;
    use umol_graph_core::Correspondence as GraphCoreCorrespondence;
    use umol_graph_ir::dsl::{
        AtomDsl as GraphIrAtomDsl, MoleculeMetadata as GraphIrMoleculeMetadata,
    };
    use umol_graph_ir::ir::{
        AromaticSystemForm as GraphIrAromaticSystemForm, AtomForm as GraphIrAtomForm,
        BondForm as GraphIrBondForm, Constraint as GraphIrConstraint,
        DativeBondForm as GraphIrDativeBondForm, Entity as GraphIrEntity,
        MoleculeAst as GraphIrMoleculeAst, MoleculeConstraint as GraphIrMoleculeConstraint,
        MoleculeCorrespondence as GraphIrMoleculeCorrespondence,
        MoleculeEntries as GraphIrMoleculeEntries,
        MulticenterBondForm as GraphIrMulticenterBondForm,
        NoncovalentBondForm as GraphIrNoncovalentBondForm,
        NoncovalentBondKind as GraphIrNoncovalentBondKind, NumForm as GraphIrNumForm,
        ReactionAst as GraphIrReactionAst, StereoAtomAst as GraphIrStereoAtomAst,
        StereoBondAst as GraphIrStereoBondAst, StereoCoset as GraphIrStereoCoset,
        StereoKind as GraphIrStereoKind, StereoLigand as GraphIrStereoLigand,
        StereoLigandKind as GraphIrStereoLigandKind,
    };

    use super::*;
    use crate::convert::into_py_variant;
    use crate::error::{MetadataError, ParseError};

    #[rstest]
    #[case::required(
        r#"{:atoms ["C" {:add "O"}]}"#,
        None,
        GraphIrReactionSpanAst::from_entries(GraphIrReactionSpanEntries {
            atoms: vec![
                GraphIrEntitySpan::Unchanged(GraphIrAtomForm::from_element(ChemElement::C)),
                GraphIrEntitySpan::Added(GraphIrAtomForm::from_element(ChemElement::O)),
            ],
            ..Default::default()
        })
    )]
    #[case::ground(
        r#"{:atoms ["C#h4#v0#d0#t0#a!#m!"]}"#,
        Some(MoleculeDefaults::ground()),
        GraphIrReactionSpanAst::from_entries(GraphIrReactionSpanEntries {
            atoms: vec![GraphIrEntitySpan::Unchanged(
                "C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!".parse().unwrap()
            )],
            ..Default::default()
        })
    )]
    fn test_reaction_span_ast_parse(
        #[case] text: &str,
        #[case] defaults: Option<MoleculeDefaults>,
        #[case] expected: GraphIrReactionSpanAst,
    ) {
        assert_eq!(
            ReactionSpanAst::parse(text, defaults).unwrap().to_rust(),
            expected
        );
    }

    #[rstest]
    fn test_reaction_span_ast_parse_error() {
        Python::attach(|py| {
            let error = ReactionSpanAst::parse("not edn", None).unwrap_err();

            assert!(error.is_instance_of::<ParseError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "EDN parse: unexpected token 'n' at byte 0"
            );
        });
    }

    #[rstest]
    fn test_reaction_span_ast_parse_with_metadata() {
        let (span, metadata) = ReactionSpanAst::parse_with_metadata(
            r#"{:atoms [[:carbon :x]] :atom-aliases [:x "C"]}"#,
            None,
        )
        .unwrap();
        let metadata = metadata.to_rust();

        assert_eq!(
            span.to_rust(),
            GraphIrReactionSpanAst::from_entries(GraphIrReactionSpanEntries {
                atoms: vec![GraphIrEntitySpan::Unchanged(GraphIrAtomForm::from_element(
                    ChemElement::C,
                ))],
                ..Default::default()
            })
        );
        assert_eq!(
            metadata.keyword(GraphIrEntity::Atom(GraphIrAtomId(0))),
            Some("carbon")
        );
        assert_eq!(
            metadata.atom_alias("x"),
            Some(&GraphIrAtomDsl(GraphIrAtomForm::from_element(
                ChemElement::C
            )))
        );
    }

    #[rstest]
    #[case::required(
        GraphIrReactionSpanAst::from_entries(GraphIrReactionSpanEntries {
            atoms: vec![
                GraphIrEntitySpan::Unchanged(GraphIrAtomForm::from_element(ChemElement::C)),
                GraphIrEntitySpan::Added(GraphIrAtomForm::from_element(ChemElement::O)),
            ],
            ..Default::default()
        }),
        None,
        r#"{:atoms ["C" {:add "O"}]}"#
    )]
    #[case::ground(
        GraphIrReactionSpanAst::from_entries(GraphIrReactionSpanEntries {
            atoms: vec![GraphIrEntitySpan::Unchanged(
                "C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!".parse().unwrap()
            )],
            ..Default::default()
        }),
        Some(MoleculeDefaults::ground()),
        r#"{:atoms ["C#h4#v0#d0#t0#a!#m!"]}"#
    )]
    fn test_reaction_span_ast_render(
        #[case] span: GraphIrReactionSpanAst,
        #[case] defaults: Option<MoleculeDefaults>,
        #[case] expected: &str,
    ) {
        assert_eq!(ReactionSpanAst::from_rust(span).render(defaults), expected);
    }

    #[rstest]
    fn test_reaction_span_ast_render_with_metadata() {
        let span = ReactionSpanAst::from_rust(GraphIrReactionSpanAst::from_entries(
            GraphIrReactionSpanEntries {
                atoms: vec![GraphIrEntitySpan::Unchanged(GraphIrAtomForm::from_element(
                    ChemElement::C,
                ))],
                ..Default::default()
            },
        ));
        let mut metadata = GraphIrMoleculeMetadata::new();
        metadata
            .set_keyword(GraphIrEntity::Atom(GraphIrAtomId(0)), "carbon")
            .unwrap();
        metadata
            .add_atom_alias(
                "x",
                GraphIrAtomDsl(GraphIrAtomForm::from_element(ChemElement::C)),
            )
            .unwrap();

        assert_eq!(
            span.render_with_metadata(&MoleculeMetadata::from_rust(metadata), None)
                .unwrap(),
            r#"{:atom-aliases [:x "C"] :atoms [[:carbon :x]]}"#
        );
    }

    #[rstest]
    fn test_reaction_span_ast_render_with_metadata_error() {
        Python::attach(|py| {
            let span = ReactionSpanAst::from_rust(GraphIrReactionSpanAst::from_entries(
                GraphIrReactionSpanEntries {
                    atoms: vec![GraphIrEntitySpan::Unchanged(GraphIrAtomForm::from_element(
                        ChemElement::C,
                    ))],
                    ..Default::default()
                },
            ));
            let mut metadata = GraphIrMoleculeMetadata::new();
            metadata
                .set_keyword(GraphIrEntity::Atom(GraphIrAtomId(1)), "outside")
                .unwrap();

            let error = span
                .render_with_metadata(&MoleculeMetadata::from_rust(metadata), None)
                .unwrap_err();

            assert!(error.is_instance_of::<MetadataError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "metadata entity is out of range: atom 1"
            );
        });
    }

    #[rstest]
    fn test_reaction_span_ast_render_with_metadata_roundtrip() {
        let (span, metadata) = ReactionSpanAst::parse_with_metadata(
            concat!(
                r#"{:atoms [[:carbon "C"] {:add "O"}] "#,
                r#":bonds [{:add [0 1 :single]}]}"#,
            ),
            None,
        )
        .unwrap();

        let rendered = span.render_with_metadata(&metadata, None).unwrap();

        assert_eq!(
            ReactionSpanAst::parse_with_metadata(&rendered, None).unwrap(),
            (span, metadata)
        );
    }

    #[rstest]
    fn test_reaction_span_ast_str() {
        let span = ReactionSpanAst::from_rust(GraphIrReactionSpanAst::from_entries(
            GraphIrReactionSpanEntries {
                atoms: vec![GraphIrEntitySpan::Unchanged(GraphIrAtomForm::from_element(
                    ChemElement::C,
                ))],
                ..Default::default()
            },
        ));

        assert_eq!(span.__str__(), span.render(None));
    }

    #[rstest]
    fn test_reaction_span_ast_from_entries() {
        Python::attach(|py| {
            let canonical_lhs =
                GraphIrAtomForm::from_element(ChemElement::C).with_charge(GraphIrNumForm::Lit(1));
            let canonical_rhs = GraphIrAtomForm::from_element(ChemElement::C)
                .with_charge(GraphIrNumForm::lit_set([1_i64]));
            let modified_lhs = GraphIrAtomForm::from_element(ChemElement::C);
            let modified_rhs = GraphIrAtomForm::from_element(ChemElement::N);
            let removed_atom = GraphIrAtomForm::from_element(ChemElement::O);
            let added_atom = GraphIrAtomForm::from_element(ChemElement::F);
            let unchanged_bond = GraphIrBondForm::from_order(1);
            let dative_lhs = GraphIrDativeBondForm::from_order(1);
            let dative_rhs = GraphIrDativeBondForm::from_order(2);
            let added_aromatic = GraphIrAromaticSystemForm::from_electrons(vec![1, 1]);
            let removed_multicenter = GraphIrMulticenterBondForm::from_electrons(vec![1, 1]);
            let unchanged_noncovalent =
                GraphIrNoncovalentBondForm::from_kind(GraphIrNoncovalentBondKind::HydrogenBond);
            let stereo_atom_lhs = GraphIrStereoAtomAst::new(
                GraphIrStereoKind::Tetrahedral,
                GraphIrStereoCoset::Lit(1),
            );
            let stereo_atom_rhs = GraphIrStereoAtomAst::new(
                GraphIrStereoKind::Tetrahedral,
                GraphIrStereoCoset::Lit(2),
            );
            let added_stereo_bond =
                GraphIrStereoBondAst::new(GraphIrStereoKind::CisTrans, GraphIrStereoCoset::Lit(1));
            let ligand = GraphIrStereoLigand::new(GraphIrAtomId(1), GraphIrStereoLigandKind::Atom);
            let unchanged_constraint =
                GraphIrConstraint::Molecule(GraphIrMoleculeConstraint::Connected { atoms: None });
            let modified_constraint_lhs =
                GraphIrConstraint::Molecule(GraphIrMoleculeConstraint::Connected {
                    atoms: Some(vec![GraphIrAtomId(0)]),
                });
            let modified_constraint_rhs =
                GraphIrConstraint::Molecule(GraphIrMoleculeConstraint::Connected {
                    atoms: Some(vec![GraphIrAtomId(1)]),
                });
            let removed_constraint =
                GraphIrConstraint::Molecule(GraphIrMoleculeConstraint::Connected {
                    atoms: Some(vec![GraphIrAtomId(2)]),
                });
            let added_constraint =
                GraphIrConstraint::Molecule(GraphIrMoleculeConstraint::Connected {
                    atoms: Some(vec![GraphIrAtomId(3)]),
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
                GraphIrReactionSpanAst::from_entries(GraphIrReactionSpanEntries {
                    atoms: vec![
                        GraphIrEntitySpan::Unchanged(canonical_lhs),
                        GraphIrEntitySpan::Modified {
                            lhs: modified_lhs,
                            rhs: modified_rhs,
                        },
                        GraphIrEntitySpan::Removed(removed_atom),
                        GraphIrEntitySpan::Added(added_atom),
                    ],
                    bonds: vec![(
                        GraphIrAtomId(0),
                        GraphIrAtomId(1),
                        GraphIrEntitySpan::Unchanged(unchanged_bond),
                    )],
                    dative: vec![(
                        vec![GraphIrAtomId(1)],
                        GraphIrAtomId(0),
                        GraphIrEntitySpan::Modified {
                            lhs: dative_lhs,
                            rhs: dative_rhs,
                        },
                    )],
                    aromatic: vec![(
                        vec![GraphIrAtomId(0), GraphIrAtomId(1)],
                        GraphIrEntitySpan::Added(added_aromatic),
                    )],
                    multicenter: vec![(
                        vec![GraphIrAtomId(0), GraphIrAtomId(1)],
                        GraphIrEntitySpan::Removed(removed_multicenter),
                    )],
                    noncovalent: vec![(
                        GraphIrAtomId(0),
                        GraphIrAtomId(1),
                        GraphIrEntitySpan::Unchanged(unchanged_noncovalent),
                    )],
                    stereo_atoms: vec![(
                        GraphIrAtomId(0),
                        vec![ligand],
                        GraphIrEntitySpan::Modified {
                            lhs: stereo_atom_lhs,
                            rhs: stereo_atom_rhs,
                        },
                    )],
                    stereo_bonds: vec![(
                        GraphIrBondId(0),
                        vec![ligand],
                        GraphIrEntitySpan::Added(added_stereo_bond),
                    )],
                    constraints: vec![
                        GraphIrConstraintSpan::Unchanged(unchanged_constraint),
                        GraphIrConstraintSpan::Removed(modified_constraint_lhs),
                        GraphIrConstraintSpan::Added(modified_constraint_rhs),
                        GraphIrConstraintSpan::Removed(removed_constraint),
                        GraphIrConstraintSpan::Added(added_constraint),
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
                        AtomAst::from_inner(GraphIrAtomForm::from_element(ChemElement::C)),
                    )
                    .unwrap()
                }),
                first_on_rhs.then(|| {
                    Py::new(
                        py,
                        AtomAst::from_inner(GraphIrAtomForm::from_element(ChemElement::C)),
                    )
                    .unwrap()
                }),
            )];
            if include_second {
                atoms.push((
                    Some(
                        Py::new(
                            py,
                            AtomAst::from_inner(GraphIrAtomForm::from_element(ChemElement::O)),
                        )
                        .unwrap(),
                    ),
                    Some(
                        Py::new(
                            py,
                            AtomAst::from_inner(GraphIrAtomForm::from_element(ChemElement::O)),
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
                        Some(
                            Py::new(py, BondAst::from_inner(GraphIrBondForm::from_order(1)))
                                .unwrap(),
                        ),
                        Some(
                            Py::new(py, BondAst::from_inner(GraphIrBondForm::from_order(1)))
                                .unwrap(),
                        ),
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

    #[rstest]
    fn test_reaction_span_ast_lhs() {
        let span = ReactionSpanAst::from_rust(
            r#"{:atoms ["C" {:remove "O"} {:add "N"}] :bonds [{:remove [0 1 :single]} {:add [0 2 :double]}]}"#
                .parse()
                .unwrap(),
        );

        assert_eq!(
            span.lhs().inner(),
            &GraphIrMoleculeAst::from_entries(GraphIrMoleculeEntries {
                atoms: vec![
                    GraphIrAtomForm::from_element(ChemElement::C),
                    GraphIrAtomForm::from_element(ChemElement::O),
                ],
                bonds: vec![(
                    GraphIrAtomId(0),
                    GraphIrAtomId(1),
                    GraphIrBondForm::from_order(1)
                )],
                ..Default::default()
            })
        );
    }

    #[rstest]
    fn test_reaction_span_ast_rhs() {
        let span = ReactionSpanAst::from_rust(
            r#"{:atoms ["C" {:remove "O"} {:add "N"}] :bonds [{:remove [0 1 :single]} {:add [0 2 :double]}]}"#
                .parse()
                .unwrap(),
        );

        assert_eq!(
            span.rhs().inner(),
            &GraphIrMoleculeAst::from_entries(GraphIrMoleculeEntries {
                atoms: vec![
                    GraphIrAtomForm::from_element(ChemElement::C),
                    GraphIrAtomForm::from_element(ChemElement::N),
                ],
                bonds: vec![(
                    GraphIrAtomId(0),
                    GraphIrAtomId(1),
                    GraphIrBondForm::from_order(2)
                )],
                ..Default::default()
            })
        );
    }

    #[rstest]
    fn test_reaction_span_ast_correspondence() {
        let span = ReactionSpanAst::from_rust(
            r#"{:atoms ["C" {:remove "O"} {:add "N"}] :bonds [{:remove [0 1 :single]} {:add [0 2 :double]}]}"#
                .parse()
                .unwrap(),
        );
        let expected = GraphIrMoleculeCorrespondence::new(
            GraphCoreCorrespondence::new(vec![(GraphIrAtomId(0), GraphIrAtomId(0))], 2, 2)
                .expect("correspondence producer preserves partial-bijection invariants"),
            GraphCoreCorrespondence::new(Vec::new(), 1, 1)
                .expect("correspondence producer preserves partial-bijection invariants"),
            GraphCoreCorrespondence::new(Vec::new(), 0, 0)
                .expect("correspondence producer preserves partial-bijection invariants"),
            GraphCoreCorrespondence::new(Vec::new(), 0, 0)
                .expect("correspondence producer preserves partial-bijection invariants"),
            GraphCoreCorrespondence::new(Vec::new(), 0, 0)
                .expect("correspondence producer preserves partial-bijection invariants"),
            GraphCoreCorrespondence::new(Vec::new(), 0, 0)
                .expect("correspondence producer preserves partial-bijection invariants"),
            GraphCoreCorrespondence::new(Vec::new(), 0, 0)
                .expect("correspondence producer preserves partial-bijection invariants"),
            GraphCoreCorrespondence::new(Vec::new(), 0, 0)
                .expect("correspondence producer preserves partial-bijection invariants"),
        );

        assert_eq!(span.correspondence().inner(), &expected);
    }

    #[rstest]
    fn test_reaction_span_ast_to_reaction() {
        Python::attach(|py| {
            let span: GraphIrReactionSpanAst =
                r#"{:atoms ["C" {:remove "O"} {:add "N"}] :bonds [{:remove [0 1 :single]} {:add [0 2 :double]}]}"#
                    .parse()
                    .unwrap();
            let expected: GraphIrReactionAst = span.to_reaction();

            assert_eq!(
                ReactionSpanAst::from_rust(span)
                    .to_reaction(py)
                    .unwrap()
                    .to_rust(py),
                expected
            );
        });
    }
}
