//! Checked writes preserve constraint order and multiplicity; invalid writes preserve the molecule.

use rstest::{fixture, rstest};
use umol_perm::Permutation;

use crate::ir::{
    AromaticSystemConstraintForm, AromaticSystemId, AtomConstraintForm, AtomForm, AtomId,
    BondConstraintForm, BondForm, BondId, BooleanForm, Constraint, Constraints, ConstraintsViewMut,
    DativeBondConstraintForm, DativeBondId, Entity, FluxionalityForm, LigandPermutation, Molecule,
    MoleculeConstraint, MoleculeEntries, MoleculeIntegrityError, MulticenterBondConstraintForm,
    MulticenterBondId, NoncovalentBondConstraintForm, NoncovalentBondId, NumForm,
    RelationalConstraint, StereoAtomConstraintForm, StereoAtomId, StereoBondConstraintForm,
    StereoBondId, StereoKind, StereoLigand, StereoLigandKind, StereoLigandPair, StereogenicityForm,
    TopicityForm, TopicityRelationForm,
};

#[fixture]
fn molecule() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::default(); 4],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::default()),
            (AtomId(1), AtomId(2), BondForm::default()),
            (AtomId(2), AtomId(3), BondForm::default()),
        ],
        dative: vec![(vec![AtomId(0)], AtomId(2), Default::default())],
        aromatic: vec![(vec![AtomId(0), AtomId(1)], Default::default())],
        multicenter: vec![(vec![AtomId(2), AtomId(3)], Default::default())],
        noncovalent: vec![([AtomId(0), AtomId(3)], Default::default())],
        stereo_atoms: vec![(
            AtomId(1),
            vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
            ],
            Default::default(),
        )],
        stereo_bonds: vec![(
            BondId(1),
            vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
            ],
            Default::default(),
        )],
        constraints: Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)).into(),
    })
}

#[rstest]
fn test_molecule_constraints_mut(mut molecule: Molecule) {
    let expected = molecule.constraints().clone();
    let view: ConstraintsViewMut<'_> = molecule.constraints_mut();
    assert_eq!(view.len(), expected.len());
    assert_eq!(view.is_empty(), expected.is_empty());
    assert_eq!(view.as_slice(), expected.as_slice());
    assert_eq!(
        view.iter().collect::<Vec<_>>(),
        expected.iter().collect::<Vec<_>>()
    );
}

#[rstest]
#[case::atom(Constraint::Atom(AtomId(0), AtomConstraintForm::valence(3)))]
#[case::bond(Constraint::Bond(BondId(0), BondConstraintForm::Aromatic(BooleanForm::Lit(true))))]
#[case::dative(Constraint::DativeBond(DativeBondId(0), DativeBondConstraintForm::aromatic(true)))]
#[case::aromatic(Constraint::AromaticSystem(
    AromaticSystemId(0),
    AromaticSystemConstraintForm::electron_count(6)
))]
#[case::multicenter(Constraint::MulticenterBond(
    MulticenterBondId(0),
    MulticenterBondConstraintForm::electron_count(2)
))]
#[case::noncovalent(Constraint::NoncovalentBond(
    NoncovalentBondId(0),
    NoncovalentBondConstraintForm::intramolecular(true)
))]
#[case::stereo_atom(Constraint::StereoAtom(
    StereoAtomId(0),
    StereoKind::Tetrahedral,
    StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)
))]
#[case::stereo_bond(Constraint::StereoBond(
    StereoBondId(0),
    StereoKind::CisTrans,
    StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)
))]
#[case::nested(Constraint::Not(Box::new(Constraint::Or(vec![Constraint::And(vec![
    Constraint::Molecule(MoleculeConstraint::Connected { atoms: Some(vec![AtomId(0), AtomId(3)]) }),
    Constraint::Relational(RelationalConstraint::AromaticSystemContains { system: AromaticSystemId(0), atom: AtomId(3) }),
])]))))]
fn test_constraints_view_mut_push(mut molecule: Molecule, #[case] constraint: Constraint) {
    let mut expected = molecule.edit();
    expected.constraints_mut().push(constraint.clone());
    assert_eq!(molecule.constraints_mut().push(constraint), Ok(()));
    assert_eq!(molecule, expected.build());
}

#[rstest]
#[case::atom(
    Constraint::Atom(AtomId(4), AtomConstraintForm::valence(4)),
    Entity::Atom(AtomId(4))
)]
#[case::bond(
    Constraint::Bond(BondId(3), BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
    Entity::Bond(BondId(3))
)]
#[case::dative(
    Constraint::DativeBond(DativeBondId(1), DativeBondConstraintForm::aromatic(true)),
    Entity::DativeBond(DativeBondId(1))
)]
#[case::aromatic(
    Constraint::AromaticSystem(
        AromaticSystemId(1),
        AromaticSystemConstraintForm::electron_count(6)
    ),
    Entity::AromaticSystem(AromaticSystemId(1))
)]
#[case::multicenter(
    Constraint::MulticenterBond(
        MulticenterBondId(1),
        MulticenterBondConstraintForm::electron_count(2)
    ),
    Entity::MulticenterBond(MulticenterBondId(1))
)]
#[case::noncovalent(
    Constraint::NoncovalentBond(
        NoncovalentBondId(1),
        NoncovalentBondConstraintForm::intramolecular(true)
    ),
    Entity::NoncovalentBond(NoncovalentBondId(1))
)]
#[case::stereo_atom(
    Constraint::StereoAtom(
        StereoAtomId(1),
        StereoKind::Tetrahedral,
        StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)
    ),
    Entity::StereoAtom(StereoAtomId(1))
)]
#[case::stereo_bond(
    Constraint::StereoBond(
        StereoBondId(1),
        StereoKind::CisTrans,
        StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)
    ),
    Entity::StereoBond(StereoBondId(1))
)]
#[case::nested(Constraint::Not(Box::new(Constraint::Or(vec![Constraint::And(vec![
    Constraint::Molecule(MoleculeConstraint::Connected { atoms: Some(vec![AtomId(0), AtomId(4)]) }),
])]))), Entity::Atom(AtomId(4)))]
#[case::relational(Constraint::Relational(RelationalConstraint::AromaticSystemContains {
    system: AromaticSystemId(0), atom: AtomId(4),
}), Entity::Atom(AtomId(4)))]
fn test_constraints_view_mut_push_error(
    mut molecule: Molecule,
    #[case] constraint: Constraint,
    #[case] entity: Entity,
) {
    let original = molecule.clone();
    assert_eq!(
        molecule.constraints_mut().push(constraint),
        Err(MoleculeIntegrityError::InvalidReference { entity })
    );
    assert_eq!(molecule, original);
}

#[rstest]
#[case::arity(Constraint::StereoAtom(StereoAtomId(0), StereoKind::TrigonalBipyramidal,
    StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
    MoleculeIntegrityError::StereoLigandArity { entity: Entity::StereoAtom(StereoAtomId(0)),
        kind: StereoKind::TrigonalBipyramidal, expected: 5, actual: 4 })]
#[case::permutation(Constraint::StereoAtom(StereoAtomId(0), StereoKind::Tetrahedral,
    StereoAtomConstraintForm::Fluxionality(FluxionalityForm {
        permutation: LigandPermutation(Permutation::identity(3)), active: BooleanForm::Lit(true),
    })), MoleculeIntegrityError::StereoPermutationDegree {
        entity: Entity::StereoAtom(StereoAtomId(0)), expected: 4, actual: 3 })]
#[case::position(Constraint::Not(Box::new(Constraint::StereoBond(StereoBondId(0), StereoKind::CisTrans,
    StereoBondConstraintForm::Topicity(TopicityForm {
        pair: StereoLigandPair::new(0usize.into(), 4usize.into()), relation: TopicityRelationForm::Undetermined,
    })))), MoleculeIntegrityError::StereoLigandPositionOutOfRange {
        entity: Entity::StereoBond(StereoBondId(0)), position: 4, degree: 4 })]
fn test_constraints_view_mut_push_stereo(
    mut molecule: Molecule,
    #[case] constraint: Constraint,
    #[case] error: MoleculeIntegrityError,
) {
    let original = molecule.clone();
    assert_eq!(molecule.constraints_mut().push(constraint), Err(error));
    assert_eq!(molecule, original);
}

#[rstest]
#[case::empty(vec![])]
#[case::duplicates(vec![
    Constraint::Atom(AtomId(1), AtomConstraintForm::valence(3)),
    Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
    Constraint::Atom(AtomId(1), AtomConstraintForm::valence(3)),
])]
fn test_constraints_view_mut_extend(mut molecule: Molecule, #[case] incoming: Vec<Constraint>) {
    let mut expected = molecule.edit();
    for constraint in &incoming {
        expected.constraints_mut().push(constraint.clone());
    }
    assert_eq!(molecule.constraints_mut().extend(incoming.into()), Ok(()));
    assert_eq!(molecule, expected.build());
}

#[rstest]
fn test_constraints_view_mut_extend_error(mut molecule: Molecule) {
    let original = molecule.clone();
    let incoming = vec![
        Constraint::Atom(AtomId(0), AtomConstraintForm::valence(3)),
        Constraint::Atom(AtomId(4), AtomConstraintForm::valence(4)),
    ]
    .into();
    assert_eq!(
        molecule.constraints_mut().extend(incoming),
        Err(MoleculeIntegrityError::InvalidReference {
            entity: Entity::Atom(AtomId(4))
        })
    );
    assert_eq!(molecule, original);
}

#[rstest]
#[case::empty(vec![])]
#[case::duplicates(vec![
    Constraint::Atom(AtomId(1), AtomConstraintForm::valence(3)),
    Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
    Constraint::Atom(AtomId(1), AtomConstraintForm::valence(3)),
])]
fn test_constraints_view_mut_replace(mut molecule: Molecule, #[case] incoming: Vec<Constraint>) {
    let mut expected = molecule.edit();
    *expected.constraints_mut() = incoming.clone().into();
    assert_eq!(molecule.constraints_mut().replace(incoming.into()), Ok(()));
    assert_eq!(molecule, expected.build());
}

#[rstest]
fn test_constraints_view_mut_replace_error(mut molecule: Molecule) {
    let original = molecule.clone();
    let incoming: Constraints = vec![
        Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: None,
            sum: NumForm::Lit(1),
        }),
        Constraint::StereoAtom(
            StereoAtomId(0),
            StereoKind::TrigonalBipyramidal,
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
        ),
    ]
    .into();
    assert_eq!(
        molecule.constraints_mut().replace(incoming),
        Err(MoleculeIntegrityError::StereoLigandArity {
            entity: Entity::StereoAtom(StereoAtomId(0)),
            kind: StereoKind::TrigonalBipyramidal,
            expected: 5,
            actual: 4
        })
    );
    assert_eq!(molecule, original);
}

#[rstest]
fn test_constraints_view_mut_remove_at(mut molecule: Molecule) {
    let mut expected = molecule.edit();
    expected.constraints_mut().clear();
    assert_eq!(
        molecule.constraints_mut().remove_at(0),
        Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4))
    );
    assert_eq!(molecule, expected.build());
}

#[rstest]
#[case::at_length(1)]
#[case::beyond_length(2)]
#[should_panic]
fn test_constraints_view_mut_remove_at_error(mut molecule: Molecule, #[case] position: usize) {
    molecule.constraints_mut().remove_at(position);
}

#[rstest]
fn test_constraints_view_mut_clear(mut molecule: Molecule) {
    let mut expected = molecule.edit();
    expected.constraints_mut().clear();
    molecule.constraints_mut().clear();
    assert_eq!(molecule, expected.build());
    let view = molecule.constraints_mut();
    assert_eq!(view.len(), 0);
    assert!(view.is_empty());
    assert_eq!(view.as_slice(), &[]);
    assert_eq!(view.iter().next(), None);
}
