//! Public construction and first-use boundaries for entity attributes.

use rstest::rstest;
use umol_chem::element::Element;
use umol_edn::{FromEdn, ToEdn};
use umol_graph_core::{
    AutomorphismAlgorithm, RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm,
};
use umol_graph_ir::ir::{
    AromaticSystemForm, AromaticSystemId, AtomForm, AtomId, BondForm, BooleanForm, Canonicalize,
    CanonicalizeContext, ConstitutionColoring, Contradiction, ElectronCountsForm, EntitySpan,
    FluxionalityForm, FrameTransport, GraphSymmetryConfig, IncidenceLevel, LigandPermutation,
    Molecule, MoleculeCanonicalizeError, MoleculeEntries, MulticenterBondForm, MulticenterBondId,
    NumForm, ReactionSpan, ReactionSpanCanonicalizeError, ReactionSpanEntries,
    StereoAtomConstraintForm, StereoAtomForm, StereoAtomId, StereoConfigurationForm, StereoCoset,
    StereoKind, StereoLigand, StereoLigandKind, StereoLigandPair, StereoTerm,
    SubstructureMatchAlgorithm, SubstructureMatchConfig, Topicity, TopicityForm,
    TopicityRelationForm,
};
use umol_perm::{DynPermutation, Permutation};

#[rstest]
#[case::short(vec![1, 2], vec![NumForm::Lit(1), NumForm::Lit(2), NumForm::Undetermined], 3)]
#[case::long(vec![1, 2, 0, 5], vec![NumForm::Lit(1), NumForm::Lit(2), NumForm::Lit(0)], 8)]
fn test_molecule_electron_attributes(
    #[case] counts: Vec<i64>,
    #[case] per_atom: Vec<NumForm>,
    #[case] total: i64,
    #[values(false, true)] aromatic: bool,
) {
    let atoms = vec![AtomForm::from_element(Element::C); 3];
    let members = vec![AtomId(0), AtomId(1), AtomId(2)];
    let mut entries = MoleculeEntries {
        atoms: atoms.clone(),
        ..Default::default()
    };
    let mut span_entries = ReactionSpanEntries {
        atoms: atoms.into_iter().map(EntitySpan::Unchanged).collect(),
        ..Default::default()
    };
    if aromatic {
        let form = AromaticSystemForm::from_electrons(counts.clone());
        entries.aromatic.push((members.clone(), form.clone()));
        span_entries
            .aromatic
            .push((members, EntitySpan::Unchanged(form)));
    } else {
        let form = MulticenterBondForm::from_electrons(counts.clone());
        entries.multicenter.push((members.clone(), form.clone()));
        span_entries
            .multicenter
            .push((members, EntitySpan::Unchanged(form)));
    }
    let molecule = Molecule::try_from_entries(entries).unwrap();
    let span = ReactionSpan::try_from_entries(span_entries).unwrap();
    assert_eq!(span.lhs(), molecule);
    assert_eq!(span.rhs(), molecule);
    let editor = molecule.clone().edit();
    assert_eq!(editor.probe().unwrap(), &molecule);
    assert_eq!(editor.finish().unwrap(), molecule);
    assert_eq!(Molecule::from_edn(&molecule.to_edn()).unwrap(), molecule);
    let contributions = molecule
        .atoms()
        .iter()
        .map(|atom| {
            if aromatic {
                atom.aromatic_valence()
            } else {
                atom.multicenter_valence()
            }
        })
        .collect::<Vec<_>>();
    assert_eq!(contributions, per_atom);
    let (stored, sum) = if aromatic {
        let view = molecule.aromatic_system(AromaticSystemId(0));
        (&view.attributes().electrons, view.electron_count())
    } else {
        let view = molecule.multicenter_bond(MulticenterBondId(0));
        (&view.attributes().electrons, view.electron_count())
    };
    assert_eq!(stored, &ElectronCountsForm::Lit(counts.clone()));
    assert_eq!(sum, NumForm::Lit(total));
    assert_eq!(
        molecule
            .incidence_graph(IncidenceLevel::Full)
            .graph()
            .edge_count(),
        3
    );
    assert_eq!(
        span.incidence_graph(IncidenceLevel::Full)
            .graph()
            .edge_count(),
        3
    );
    assert_eq!(
        molecule.clone().canonicalize(&CanonicalizeContext {
            para_stereo: true,
            automorphism_algorithm: AutomorphismAlgorithm::Nauty
        }),
        Err(MoleculeCanonicalizeError::Contradiction(Contradiction))
    );
    assert_eq!(
        span.canonicalize(&CanonicalizeContext {
            para_stereo: true,
            automorphism_algorithm: AutomorphismAlgorithm::Nauty
        }),
        Err(ReactionSpanCanonicalizeError::Contradiction(Contradiction))
    );
    assert_eq!(
        ElectronCountsForm::Lit(counts).reframe_by(&DynPermutation::identity(3)),
        None
    );
    for match_algorithm in [
        SubstructureMatchAlgorithm::GraphAndOverlays,
        SubstructureMatchAlgorithm::Incidence,
    ] {
        let config = SubstructureMatchConfig {
            match_algorithm,
            subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::Vf2,
            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
        };
        assert_eq!(molecule.substructure_matches(&molecule, config), Ok(vec![]));
        let mut valid = molecule.clone();
        if aromatic {
            valid
                .aromatic_system_mut(AromaticSystemId(0))
                .attributes_mut()
                .electrons = ElectronCountsForm::Lit(vec![1, 2, 0]);
        } else {
            valid
                .multicenter_bond_mut(MulticenterBondId(0))
                .attributes_mut()
                .electrons = ElectronCountsForm::Lit(vec![1, 2, 0]);
        }
        assert_eq!(molecule.substructure_matches(&valid, config), Ok(vec![]));
        assert_eq!(valid.substructure_matches(&molecule, config), Ok(vec![]));
    }
}

#[rstest]
#[case::short(3, StereoAtomForm::new(StereoKind::Tetrahedral, 0u32))]
#[case::long(5, StereoAtomForm::new(StereoKind::Tetrahedral, 0u32))]
#[case::site_kind(4, StereoAtomForm::new(StereoKind::CisTrans, 0u32))]
#[case::coset(4, StereoAtomForm::new(StereoKind::Tetrahedral, 2u32))]
#[case::term_action(
    4,
    StereoAtomForm::new(
        StereoKind::Tetrahedral,
        StereoCoset::term(StereoTerm::apply(StereoTerm::var("x"), Permutation::identity(3)))
    )
)]
#[case::topicity(4, StereoAtomForm {
    configuration: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32),
    constraints: StereoAtomConstraintForm::Topicity(TopicityForm {
        pair: StereoLigandPair::new(0usize.into(), 4usize.into()),
        relation: TopicityRelationForm::Lit(Topicity::Homotopic),
    }).into(),
})]
#[case::constraint_action(4, StereoAtomForm {
    configuration: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32),
    constraints: StereoAtomConstraintForm::Fluxionality(FluxionalityForm {
        permutation: LigandPermutation(Permutation::identity(3)),
        active: BooleanForm::Lit(true),
    }).into(),
})]
fn test_molecule_stereo_attributes(#[case] degree: u32, #[case] form: StereoAtomForm) {
    let atoms = vec![AtomForm::from_element(Element::C); degree as usize + 1];
    let bonds = (1..=degree)
        .map(|i| (AtomId(0), AtomId(i), BondForm::from_order(1)))
        .collect::<Vec<_>>();
    let ligands = (1..=degree)
        .map(|i| StereoLigand::new(AtomId(i), StereoLigandKind::Atom))
        .collect::<Vec<_>>();
    let molecule = Molecule::try_from_entries(MoleculeEntries {
        atoms: atoms.clone(),
        bonds: bonds.clone(),
        stereo_atoms: vec![(AtomId(0), ligands.clone(), form.clone())],
        ..Default::default()
    })
    .unwrap();
    let span = ReactionSpan::try_from_entries(ReactionSpanEntries {
        atoms: atoms.into_iter().map(EntitySpan::Unchanged).collect(),
        bonds: bonds
            .into_iter()
            .map(|(a, b, f)| (a, b, EntitySpan::Unchanged(f)))
            .collect(),
        stereo_atoms: vec![(AtomId(0), ligands, EntitySpan::Unchanged(form.clone()))],
        ..Default::default()
    })
    .unwrap();
    assert_eq!(molecule.stereo_atom(StereoAtomId(0)).attributes(), &form);
    assert_eq!(span.lhs(), molecule);
    assert_eq!(span.rhs(), molecule);
    assert_eq!(molecule.clone().edit().finish().unwrap(), molecule);
    let context = CanonicalizeContext {
        para_stereo: true,
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    };
    let canonical_molecule = molecule.clone().canonicalize(&context);
    let canonical_span = span.canonicalize(&context);
    // Malformed symbolic actions require panic freedom, not a particular result.
    if !matches!(
        form.configuration,
        StereoConfigurationForm::Kinded(_, StereoCoset::Term(_))
    ) {
        assert_eq!(
            canonical_molecule,
            Err(MoleculeCanonicalizeError::Contradiction(Contradiction))
        );
        assert_eq!(
            canonical_span,
            Err(ReactionSpanCanonicalizeError::Contradiction(Contradiction))
        );
    }
    // Symmetry is infallible; malformed attributes have no classification guarantee.
    molecule.graph_symmetry(&GraphSymmetryConfig {
        coloring: ConstitutionColoring::full(),
        iterate_to_fixpoint: true,
        max_iterations: 2,
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    });
}

#[rstest]
#[case::short(vec![1, 2])]
#[case::long(vec![1, 2, 0, 5])]
fn test_reaction_span_electron_attributes(
    #[case] counts: Vec<i64>,
    #[values(false, true)] aromatic: bool,
    #[values(false, true)] rhs: bool,
) {
    let (lhs_counts, rhs_counts) = if rhs {
        (vec![1, 2, 0], counts)
    } else {
        (counts, vec![1, 2, 0])
    };
    let mut entries = ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 3],
        ..Default::default()
    };
    let atoms = vec![AtomId(0), AtomId(1), AtomId(2)];
    if aromatic {
        entries.aromatic.push((
            atoms,
            EntitySpan::Modified {
                lhs: AromaticSystemForm::from_electrons(lhs_counts.clone()),
                rhs: AromaticSystemForm::from_electrons(rhs_counts.clone()),
            },
        ));
    } else {
        entries.multicenter.push((
            atoms,
            EntitySpan::Modified {
                lhs: MulticenterBondForm::from_electrons(lhs_counts.clone()),
                rhs: MulticenterBondForm::from_electrons(rhs_counts.clone()),
            },
        ));
    }
    let span = ReactionSpan::try_from_entries(entries).unwrap();
    for (molecule, expected) in [(span.lhs(), lhs_counts), (span.rhs(), rhs_counts)] {
        let stored = if aromatic {
            &molecule
                .aromatic_system(AromaticSystemId(0))
                .attributes()
                .electrons
        } else {
            &molecule
                .multicenter_bond(MulticenterBondId(0))
                .attributes()
                .electrons
        };
        assert_eq!(stored, &ElectronCountsForm::Lit(expected));
    }
    assert_eq!(
        span.canonicalize(&CanonicalizeContext {
            para_stereo: false,
            automorphism_algorithm: AutomorphismAlgorithm::Nauty,
        }),
        Err(ReactionSpanCanonicalizeError::Contradiction(Contradiction))
    );
}
