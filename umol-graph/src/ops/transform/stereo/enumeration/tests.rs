use std::collections::BTreeSet;

use rstest::{fixture, rstest};
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{
    AtomConstraintForm, AtomId, BondConstraintForm, BondId, BooleanForm, CisTransStereoForm,
    ConstitutionColoring, Constraint, DativeBondConstraintForm, DativeBondId, GraphSymmetryConfig,
    LigandPermutation, LigandSymmetryForm, Molecule, NumForm, OrientedLigandPermutation, RingScope,
    StereoAtomConstraintForm, StereoAtomId, StereoBondConstraintForm, StereoBondId,
    StereoConfigurationForm, StereoCoset, StereoKind, StereoLigandPair, StereoLigandPosition,
    StereoTerm, Stereogenicity, StereogenicityForm, TetrahedralStereoForm, Topicity, TopicityForm,
    TopicityRelationForm,
};
use umol_graph_ir::mol_dsl_concrete;
use umol_perm::{Orientation, Permutation};

use super::{AbsoluteStereoAssignments, StereoAssignmentChecks};
use crate::ops::stereo::{StereoInconsistency, StereoPerceptionError, StereoSite};
use crate::ops::transform::tests::stereo::reference::absolute_assignments;
use crate::ops::validate::constraint::ConstraintInvariantsError;

#[fixture]
fn mixed_stereo() -> Molecule {
    mol_dsl_concrete!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2"] [5 7 "1"] [6 8 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]
        :stereo-bonds [{:site 4 :ligands [7 [:h 5] 8 [:h 6]] :attrs "Ct*"}]}"#
    )
}

#[rstest]
#[case::binary(vec![StereoCoset::Undetermined, StereoCoset::Undetermined], vec![0, 1])]
#[case::reverse(vec![StereoCoset::Undetermined, StereoCoset::Undetermined], vec![1, 0])]
#[case::fixed(vec![StereoCoset::Lit(1), StereoCoset::Undetermined], vec![1, 0])]
#[case::finite(vec![StereoCoset::lit_set([0, 1]), StereoCoset::lit_set([1])], vec![0, 1])]
#[case::empty_domain(vec![StereoCoset::lit_set([]), StereoCoset::Undetermined], vec![1, 0])]
#[case::empty_selection(vec![StereoCoset::Undetermined, StereoCoset::Undetermined], vec![])]
#[case::excluded_variable(vec![StereoCoset::term(StereoTerm::var("retained")), StereoCoset::Undetermined], vec![1])]
fn test_absolute_stereo_assignments_new(
    #[case] source: Vec<StereoCoset>,
    #[case] order: Vec<usize>,
) {
    let sites = [StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(4))];
    let kinds = [StereoKind::Tetrahedral, StereoKind::CisTrans];
    let selected = order.iter().copied().collect();
    let expected = absolute_assignments(&source, &selected, &order, |_| true).unwrap();
    let assignments = AbsoluteStereoAssignments::new(order.iter().map(|&index| {
        (
            sites[index],
            StereoConfigurationForm::Kinded(kinds[index], source[index].clone()),
        )
    }))
    .unwrap();
    let actual: Vec<_> = assignments
        .map(|assignment| {
            let mut values = source.clone();
            for (site, configuration) in assignment {
                let index = sites
                    .iter()
                    .position(|&candidate| candidate == site)
                    .unwrap();
                values[index] = configuration.coset().unwrap().clone();
            }
            values
        })
        .collect();
    assert_eq!(actual, expected);
}

#[rstest]
#[case::swap(StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(0))), vec![1])]
#[case::mirror(StereoCoset::term(StereoTerm::mirror(StereoTerm::Lit(1))), vec![0])]
#[case::finite(StereoCoset::term(StereoTerm::swap(StereoTerm::lit_set([0, 1]))), vec![0, 1])]
#[case::empty(StereoCoset::term(StereoTerm::lit_set([])), vec![])]
fn test_absolute_stereo_assignments_new_expression(
    #[case] source: StereoCoset,
    #[case] expected: Vec<u32>,
) {
    let site = StereoSite::Atom(AtomId(0));
    let assignments = AbsoluteStereoAssignments::new([(
        site,
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, source),
    )])
    .unwrap();
    assert_eq!(
        assignments.collect::<Vec<_>>(),
        expected
            .into_iter()
            .map(|value| vec![(site, (StereoKind::Tetrahedral, value).into())])
            .collect::<Vec<_>>()
    );
}

#[rstest]
#[case::kindless(StereoConfigurationForm::Undetermined, StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Atom(AtomId(0)) })]
#[case::variable(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::term(StereoTerm::var("configuration"))),
    StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Atom(AtomId(0)) })]
#[case::bounded_variable(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::term(StereoTerm::var_in("configuration", [0, 1]))),
    StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Atom(AtomId(0)) })]
#[case::relative(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::term(StereoTerm::mirror(StereoTerm::var("reference")))),
    StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Atom(AtomId(0)) })]
#[case::empty_variable(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::term(StereoTerm::var_in("configuration", []))),
    StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Atom(AtomId(0)) })]
#[case::invalid_value((StereoKind::Tetrahedral, 2).into(), StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Atom(AtomId(0)) })]
#[case::invalid_constant(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::term(StereoTerm::Lit(2))),
    StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Atom(AtomId(0)) })]
#[case::invalid_set(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::lit_set([0, 2])),
    StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Atom(AtomId(0)) })]
#[case::square_planar((StereoKind::SquarePlanar, 0).into(), StereoPerceptionError::UnsupportedKind { kind: StereoKind::SquarePlanar })]
#[case::site_kind((StereoKind::CisTrans, 0).into(), StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Atom(AtomId(0)) })]
fn test_absolute_stereo_assignments_new_error(
    #[case] configuration: StereoConfigurationForm,
    #[case] expected: StereoPerceptionError,
) {
    assert_eq!(
        AbsoluteStereoAssignments::new([(StereoSite::Atom(AtomId(0)), configuration)]).err(),
        Some(expected)
    );
}

#[rstest]
fn test_absolute_stereo_assignments_next(
    #[values(0usize, 1, 2, 3, 4, 5, 6, 7, 8)] count: usize,
    #[values(false, true)] reverse: bool,
) {
    let source = vec![StereoCoset::Undetermined; count];
    let mut order: Vec<_> = (0..count).collect();
    if reverse {
        order.reverse();
    }
    let expected = absolute_assignments(&source, &(0..count).collect(), &order, |_| true).unwrap();
    let mut assignments = AbsoluteStereoAssignments::new(order.iter().map(|&index| {
        (
            StereoSite::Atom(AtomId(index as u32)),
            StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined),
        )
    }))
    .unwrap();
    let actual: Vec<_> = assignments
        .by_ref()
        .map(|assignment| {
            let mut values = source.clone();
            for (site, configuration) in assignment {
                let StereoSite::Atom(atom) = site else {
                    unreachable!()
                };
                values[atom.index()] = configuration.coset().unwrap().clone();
            }
            values
        })
        .collect();
    assert_eq!(actual, expected);
    assert_eq!(assignments.next(), None);
    assert_eq!(assignments.next(), None);
}

#[rstest]
fn test_absolute_stereo_assignments_next_prefix() {
    let mut assignments = AbsoluteStereoAssignments::new((0..128).map(|index| {
        (
            StereoSite::Atom(AtomId(index)),
            StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined),
        )
    }))
    .unwrap();
    let first: Vec<_> = (0..128)
        .map(|index| {
            (
                StereoSite::Atom(AtomId(index)),
                (StereoKind::Tetrahedral, 0).into(),
            )
        })
        .collect();
    assert_eq!(assignments.next(), Some(first.clone()));
    let mut second = first;
    second[127].1 = (StereoKind::Tetrahedral, 1).into();
    assert_eq!(assignments.next(), Some(second));
    assert_eq!(assignments.positions, [vec![0; 126], vec![1, 0]].concat());
    assert_eq!(assignments.domains.len(), 128);
    assert!(!assignments.exhausted);
}

#[rstest]
#[case::variable(StereoCoset::term(StereoTerm::var("configuration")))]
#[case::bounded_variable(StereoCoset::term(StereoTerm::var_in("configuration", [0, 1])))]
#[case::empty_variable(StereoCoset::term(StereoTerm::var_in("configuration", [])))]
#[case::relative(StereoCoset::term(StereoTerm::swap(StereoTerm::var("configuration"))))]
#[case::invalid(StereoCoset::Lit(2))]
fn test_stereo_assignment_checks_new_error(
    mut mixed_stereo: Molecule,
    #[case] coset: StereoCoset,
    #[values(StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(4)))] site: StereoSite,
) {
    match site {
        StereoSite::Atom(atom) => mixed_stereo
            .atom_mut(atom)
            .attributes_mut()
            .constraints
            .set(AtomConstraintForm::TetrahedralStereo(
                TetrahedralStereoForm::Stereo(coset),
            )),
        StereoSite::Bond(bond) => mixed_stereo
            .bond_mut(bond)
            .attributes_mut()
            .constraints
            .set(BondConstraintForm::CisTransStereo(
                CisTransStereoForm::Stereo(coset),
            )),
    }
    assert_eq!(
        StereoAssignmentChecks::new(&mixed_stereo, &BTreeSet::from([site])).err(),
        Some(StereoPerceptionError::UnsupportedConfiguration { site })
    );
}

#[rstest]
#[case::invalid_atom(mol_dsl_concrete!(r#"{:atoms ["C"]}"#), StereoSite::Atom(AtomId(1)),
    StereoPerceptionError::InvalidSite { site: StereoSite::Atom(AtomId(1)) })]
#[case::invalid_bond(mol_dsl_concrete!(r#"{:atoms ["C"]}"#), StereoSite::Bond(BondId(0)),
    StereoPerceptionError::InvalidSite { site: StereoSite::Bond(BondId(0)) })]
#[case::missing_atom(mol_dsl_concrete!(r#"{:atoms ["C"]}"#), StereoSite::Atom(AtomId(0)),
    StereoPerceptionError::MissingStereoEntity { site: StereoSite::Atom(AtomId(0)) })]
#[case::missing_bond(mol_dsl_concrete!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#), StereoSite::Bond(BondId(0)),
    StereoPerceptionError::MissingStereoEntity { site: StereoSite::Bond(BondId(0)) })]
#[case::atom_frame(mol_dsl_concrete!(r#"{:atoms ["C#h" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3] :attrs "Th0"}]}"#), StereoSite::Atom(AtomId(0)),
    StereoPerceptionError::Inconsistency(StereoInconsistency::StereoAtomFailure { stereo_atom: StereoAtomId(0) }))]
fn test_stereo_assignment_checks_new_reference(
    #[case] molecule: Molecule,
    #[case] site: StereoSite,
    #[case] expected: StereoPerceptionError,
) {
    assert_eq!(
        StereoAssignmentChecks::new(&molecule, &BTreeSet::from([site])).err(),
        Some(expected)
    );
}

#[rstest]
#[case::term(StereoCoset::term(StereoTerm::var("retained")))]
#[case::empty_term(StereoCoset::term(StereoTerm::var_in("retained", [])))]
fn test_stereo_assignment_checks_new_expression(
    mut mixed_stereo: Molecule,
    #[case] coset: StereoCoset,
) {
    mixed_stereo
        .stereo_bond_mut(StereoBondId(0))
        .attributes_mut()
        .configuration = StereoConfigurationForm::Kinded(StereoKind::CisTrans, coset);
    let mut editor = mixed_stereo.edit();
    editor.constraints_mut().push(Constraint::Or(vec![
        Constraint::Atom(
            AtomId(0),
            AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Stereo(0.into())),
        ),
        Constraint::Bond(
            BondId(4),
            BondConstraintForm::CisTransStereo(CisTransStereoForm::Stereo(1.into())),
        ),
    ]));
    let molecule = editor.finish().unwrap();
    assert_eq!(
        StereoAssignmentChecks::new(&molecule, &BTreeSet::from([StereoSite::Atom(AtomId(0))]))
            .err(),
        Some(StereoPerceptionError::UnsupportedConfiguration {
            site: StereoSite::Bond(BondId(4))
        })
    );
}

#[rstest]
#[case::unrestricted(None, None, |_: &[StereoCoset]| true)]
#[case::atom(Some(0), None, |values: &[StereoCoset]| values[0] == StereoCoset::Lit(0))]
#[case::bond(None, Some(1), |values: &[StereoCoset]| values[1] == StereoCoset::Lit(1))]
#[case::both(Some(1), Some(0), |values: &[StereoCoset]| values == [StereoCoset::Lit(1), StereoCoset::Lit(0)])]
fn test_stereo_assignment_checks_is_consistent(
    mut mixed_stereo: Molecule,
    #[case] atom_assertion: Option<u32>,
    #[case] bond_assertion: Option<u32>,
    #[case] admissible: fn(&[StereoCoset]) -> bool,
    #[values(false, true)] reverse: bool,
    #[values(StereoCoset::Undetermined, StereoCoset::Lit(0), StereoCoset::lit_set([1]), StereoCoset::lit_set([]))]
    atom_domain: StereoCoset,
    #[values(StereoCoset::Undetermined, StereoCoset::Lit(1), StereoCoset::lit_set([0]), StereoCoset::lit_set([]))]
    bond_domain: StereoCoset,
) {
    let source = [atom_domain, bond_domain];
    mixed_stereo
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration =
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, source[0].clone());
    mixed_stereo
        .stereo_bond_mut(StereoBondId(0))
        .attributes_mut()
        .configuration = StereoConfigurationForm::Kinded(StereoKind::CisTrans, source[1].clone());
    if let Some(value) = atom_assertion {
        mixed_stereo
            .atom_mut(AtomId(0))
            .attributes_mut()
            .constraints
            .set(AtomConstraintForm::TetrahedralStereo(
                TetrahedralStereoForm::Stereo(value.into()),
            ));
    }
    if let Some(value) = bond_assertion {
        mixed_stereo
            .bond_mut(BondId(4))
            .attributes_mut()
            .constraints
            .set(BondConstraintForm::CisTransStereo(
                CisTransStereoForm::Stereo(value.into()),
            ));
    }
    let original = mixed_stereo.clone();
    let sites = [StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(4))];
    let kinds = [StereoKind::Tetrahedral, StereoKind::CisTrans];
    let order = if reverse { vec![1, 0] } else { vec![0, 1] };
    let checks = StereoAssignmentChecks::new(&mixed_stereo, &sites.into()).unwrap();
    let assignments = AbsoluteStereoAssignments::new(order.iter().map(|&index| {
        (
            sites[index],
            StereoConfigurationForm::Kinded(kinds[index], source[index].clone()),
        )
    }))
    .unwrap();
    let config = GraphSymmetryConfig {
        coloring: ConstitutionColoring::full(),
        iterate_to_fixpoint: false,
        max_iterations: 16,
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    };
    let mut actual = Vec::new();
    for assignment in assignments {
        let mut candidate = mixed_stereo.clone();
        for (site, configuration) in assignment {
            match site {
                StereoSite::Atom(_) => {
                    candidate
                        .stereo_atom_mut(StereoAtomId(0))
                        .attributes_mut()
                        .configuration = configuration
                }
                StereoSite::Bond(_) => {
                    candidate
                        .stereo_bond_mut(StereoBondId(0))
                        .attributes_mut()
                        .configuration = configuration
                }
            }
        }
        if checks.is_consistent(&candidate, &config).unwrap() {
            actual.push(vec![
                candidate
                    .stereo_atom(StereoAtomId(0))
                    .attributes()
                    .configuration
                    .coset()
                    .unwrap()
                    .clone(),
                candidate
                    .stereo_bond(StereoBondId(0))
                    .attributes()
                    .configuration
                    .coset()
                    .unwrap()
                    .clone(),
            ]);
            assert_eq!(candidate.clone().edit().finish(), Ok(candidate));
        }
    }
    assert_eq!(
        actual,
        absolute_assignments(&source, &BTreeSet::from([0, 1]), &order, admissible).unwrap()
    );
    assert_eq!(mixed_stereo, original);
}

#[rstest]
fn test_stereo_assignment_checks_is_consistent_logical(mut mixed_stereo: Molecule) {
    let mut editor = mixed_stereo.edit();
    editor.constraints_mut().push(Constraint::Or(vec![
        Constraint::And(vec![
            Constraint::Atom(
                AtomId(0),
                AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Stereo(0.into())),
            ),
            Constraint::Bond(
                BondId(4),
                BondConstraintForm::CisTransStereo(CisTransStereoForm::Stereo(1.into())),
            ),
        ]),
        Constraint::And(vec![
            Constraint::Atom(
                AtomId(0),
                AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Stereo(1.into())),
            ),
            Constraint::Not(Box::new(Constraint::Bond(
                BondId(4),
                BondConstraintForm::CisTransStereo(CisTransStereoForm::Stereo(1.into())),
            ))),
        ]),
    ]));
    mixed_stereo = editor.finish().unwrap();
    let checks = StereoAssignmentChecks::new(
        &mixed_stereo,
        &BTreeSet::from([StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(4))]),
    )
    .unwrap();
    let config = GraphSymmetryConfig {
        coloring: ConstitutionColoring::full(),
        iterate_to_fixpoint: false,
        max_iterations: 16,
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    };
    let mut actual = Vec::new();
    for atom in 0..2 {
        for bond in 0..2 {
            let mut candidate = mixed_stereo.clone();
            candidate
                .stereo_atom_mut(StereoAtomId(0))
                .attributes_mut()
                .configuration = (StereoKind::Tetrahedral, atom).into();
            candidate
                .stereo_bond_mut(StereoBondId(0))
                .attributes_mut()
                .configuration = (StereoKind::CisTrans, bond).into();
            if checks.is_consistent(&candidate, &config).unwrap() {
                actual.push(vec![StereoCoset::Lit(atom), StereoCoset::Lit(bond)]);
            }
        }
    }
    assert_eq!(
        actual,
        absolute_assignments(
            &[StereoCoset::Undetermined, StereoCoset::Undetermined],
            &BTreeSet::from([0, 1]),
            &[0, 1],
            |values| values[0] != values[1]
        )
        .unwrap()
    );
}

#[rstest]
#[case::prochiral(Stereogenicity::Prochiral, true)]
#[case::symmetric(Stereogenicity::Symmetric, false)]
fn test_stereo_assignment_checks_is_consistent_stereogenicity(
    #[case] asserted: Stereogenicity,
    #[case] expected: bool,
    #[values(StereoCoset::Undetermined, StereoCoset::Lit(0), StereoCoset::Lit(1))]
    coset: StereoCoset,
) {
    let mut molecule = mol_dsl_concrete!(
        r#"{:atoms ["C" "F" "Cl" "H" "H"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]}"#
    );
    molecule
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, coset);
    molecule
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .constraints
        .set(StereoAtomConstraintForm::Stereogenicity(
            StereogenicityForm::Lit(asserted),
        ));
    let checks = StereoAssignmentChecks::new(&molecule, &BTreeSet::new()).unwrap();
    assert_eq!(
        checks.is_consistent(
            &molecule,
            &GraphSymmetryConfig {
                coloring: ConstitutionColoring::full(),
                iterate_to_fixpoint: false,
                max_iterations: 16,
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
            }
        ),
        Ok(expected)
    );
}

#[rstest]
fn test_stereo_assignment_checks_is_consistent_partial(mut mixed_stereo: Molecule) {
    mixed_stereo.atom_mut(AtomId(1)).attributes_mut().charge = NumForm::Undetermined;
    let mut candidate = mixed_stereo.clone();
    candidate
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = (StereoKind::Tetrahedral, 1).into();
    let checks = StereoAssignmentChecks::new(
        &mixed_stereo,
        &BTreeSet::from([StereoSite::Atom(AtomId(0))]),
    )
    .unwrap();
    assert_eq!(
        checks.is_consistent(
            &candidate,
            &GraphSymmetryConfig {
                coloring: ConstitutionColoring::full(),
                iterate_to_fixpoint: false,
                max_iterations: 16,
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
            }
        ),
        Ok(true)
    );
    assert_eq!(
        candidate
            .stereo_bond(StereoBondId(0))
            .attributes()
            .configuration,
        StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Undetermined)
    );
    assert_eq!(
        mixed_stereo
            .stereo_atom(StereoAtomId(0))
            .attributes()
            .configuration,
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined)
    );
}

#[rstest]
#[case::atom(mol_dsl_concrete!(r#"{:atoms ["C#T1" "F" "Cl" "Br" "I"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [2 1 3 4] :attrs "Th0"}]}"#),
    StereoSite::Atom(AtomId(0)))]
#[case::bond(mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "F" "Cl"]
    :bonds [[0 1 "2#C1"] [0 2 "1"] [1 3 "1"]]
    :stereo-bonds [{:site 0 :ligands [[:h 0] 2 3 [:h 1]] :attrs "Ct0"}]}"#),
    StereoSite::Bond(BondId(0)))]
fn test_stereo_assignment_checks_is_consistent_frame(
    #[case] molecule: Molecule,
    #[case] site: StereoSite,
) {
    let checks = StereoAssignmentChecks::new(&molecule, &BTreeSet::from([site])).unwrap();
    let config = GraphSymmetryConfig {
        coloring: ConstitutionColoring::full(),
        iterate_to_fixpoint: false,
        max_iterations: 16,
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    };
    assert_eq!(checks.is_consistent(&molecule, &config), Ok(true));
    let mut opposite = molecule.clone();
    match site {
        StereoSite::Atom(_) => {
            opposite
                .stereo_atom_mut(StereoAtomId(0))
                .attributes_mut()
                .configuration = (StereoKind::Tetrahedral, 1).into()
        }
        StereoSite::Bond(_) => {
            opposite
                .stereo_bond_mut(StereoBondId(0))
                .attributes_mut()
                .configuration = (StereoKind::CisTrans, 1).into()
        }
    }
    assert_eq!(checks.is_consistent(&opposite, &config), Ok(false));
}

#[rstest]
#[case::atom_stereogenic(
    Constraint::StereoAtom(
        StereoAtomId(0),
        StereoKind::Tetrahedral,
        StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
            Stereogenicity::Stereogenic
        ))
    ),
    true
)]
#[case::atom_symmetric(
    Constraint::StereoAtom(
        StereoAtomId(0),
        StereoKind::Tetrahedral,
        StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
            Stereogenicity::Symmetric
        ))
    ),
    false
)]
#[case::bond_stereogenic(
    Constraint::StereoBond(
        StereoBondId(0),
        StereoKind::CisTrans,
        StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
            Stereogenicity::Stereogenic
        ))
    ),
    true
)]
#[case::bond_symmetric(
    Constraint::StereoBond(
        StereoBondId(0),
        StereoKind::CisTrans,
        StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
            Stereogenicity::Symmetric
        ))
    ),
    false
)]
#[case::atom_topicity(Constraint::StereoAtom(StereoAtomId(0), StereoKind::Tetrahedral, StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Diastereotopic) })), true)]
#[case::atom_topicity_mismatch(Constraint::StereoAtom(StereoAtomId(0), StereoKind::Tetrahedral, StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })), false)]
#[case::bond_topicity(Constraint::StereoBond(StereoBondId(0), StereoKind::CisTrans, StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Diastereotopic) })), true)]
#[case::bond_topicity_mismatch(Constraint::StereoBond(StereoBondId(0), StereoKind::CisTrans, StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })), false)]
#[case::atom_symmetry(Constraint::StereoAtom(StereoAtomId(0), StereoKind::Tetrahedral, StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })), true)]
#[case::atom_symmetry_mismatch(Constraint::StereoAtom(StereoAtomId(0), StereoKind::Tetrahedral, StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(false) })), false)]
#[case::bond_symmetry(Constraint::StereoBond(StereoBondId(0), StereoKind::CisTrans, StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })), true)]
#[case::bond_symmetry_mismatch(Constraint::StereoBond(StereoBondId(0), StereoKind::CisTrans, StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(false) })), false)]
fn test_stereo_assignment_checks_is_consistent_constraint(
    mixed_stereo: Molecule,
    #[case] constraint: Constraint,
    #[case] expected: bool,
) {
    let mut editor = mixed_stereo.edit();
    editor.constraints_mut().push(constraint);
    let molecule = editor.finish().unwrap();
    let checks = StereoAssignmentChecks::new(&molecule, &BTreeSet::new()).unwrap();
    assert_eq!(
        checks.is_consistent(
            &molecule,
            &GraphSymmetryConfig {
                coloring: ConstitutionColoring::full(),
                iterate_to_fixpoint: false,
                max_iterations: 16,
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
            }
        ),
        Ok(expected)
    );
}

#[rstest]
#[case::stereogenic(Stereogenicity::Stereogenic, true)]
#[case::symmetric(Stereogenicity::Symmetric, false)]
fn test_stereo_assignment_checks_is_consistent_context(
    #[case] asserted: Stereogenicity,
    #[case] expected: bool,
) {
    let mut molecule = mol_dsl_concrete!(
        r#"{:atoms ["C" "C" "Si" "F" "H" "Cl" "Br" "H" "Cl" "Br" "H"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"] [1 6 "1"] [1 7 "1"] [2 8 "1"] [2 9 "1"] [2 10 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"} {:site 1 :ligands [0 5 6 7] :attrs "Th*"} {:site 2 :ligands [0 8 9 10] :attrs "Th*"}]}"#
    );
    molecule
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .constraints
        .set(StereoAtomConstraintForm::Stereogenicity(
            StereogenicityForm::Lit(asserted),
        ));
    let original = molecule.clone();
    let checks = StereoAssignmentChecks::new(
        &molecule,
        &BTreeSet::from([StereoSite::Atom(AtomId(1)), StereoSite::Atom(AtomId(2))]),
    )
    .unwrap();
    let assignments = AbsoluteStereoAssignments::new([1, 2].map(|atom| {
        (
            StereoSite::Atom(AtomId(atom)),
            StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined),
        )
    }))
    .unwrap();
    let config = GraphSymmetryConfig {
        coloring: ConstitutionColoring::full(),
        iterate_to_fixpoint: false,
        max_iterations: 16,
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    };
    let mut actual = Vec::new();
    for assignment in assignments {
        let mut candidate = molecule.clone();
        for (site, configuration) in assignment {
            let StereoSite::Atom(atom) = site else {
                unreachable!()
            };
            candidate
                .stereo_atom_mut(StereoAtomId(atom.0))
                .attributes_mut()
                .configuration = configuration;
        }
        if checks.is_consistent(&candidate, &config).unwrap() {
            actual.push(
                [1, 2]
                    .map(|index| {
                        candidate
                            .stereo_atom(StereoAtomId(index))
                            .attributes()
                            .configuration
                            .coset()
                            .unwrap()
                            .clone()
                    })
                    .to_vec(),
            );
        }
    }
    assert_eq!(
        actual,
        absolute_assignments(
            &[StereoCoset::Undetermined, StereoCoset::Undetermined],
            &BTreeSet::from([0, 1]),
            &[0, 1],
            |_| expected
        )
        .unwrap()
    );
    assert_eq!(molecule, original);
}

#[rstest]
#[case::open(None, vec![true, true])]
#[case::fixed(Some(0), vec![true, false])]
fn test_stereo_assignment_checks_is_consistent_evidence(
    mut mixed_stereo: Molecule,
    #[case] retained_assertion: Option<u32>,
    #[case] expected: Vec<bool>,
) {
    if let Some(value) = retained_assertion {
        mixed_stereo
            .bond_mut(BondId(4))
            .attributes_mut()
            .constraints
            .set(BondConstraintForm::CisTransStereo(
                CisTransStereoForm::Stereo(value.into()),
            ));
    }
    let mut editor = mixed_stereo.edit();
    editor.constraints_mut().push(Constraint::Or(vec![
        Constraint::Atom(
            AtomId(0),
            AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Stereo(0.into())),
        ),
        Constraint::Bond(
            BondId(4),
            BondConstraintForm::CisTransStereo(CisTransStereoForm::Stereo(1.into())),
        ),
    ]));
    let molecule = editor.finish().unwrap();
    let checks =
        StereoAssignmentChecks::new(&molecule, &BTreeSet::from([StereoSite::Atom(AtomId(0))]))
            .unwrap();
    let config = GraphSymmetryConfig {
        coloring: ConstitutionColoring::full(),
        iterate_to_fixpoint: false,
        max_iterations: 16,
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    };
    let actual: Vec<_> = (0..2)
        .map(|value| {
            let mut candidate = molecule.clone();
            candidate
                .stereo_atom_mut(StereoAtomId(0))
                .attributes_mut()
                .configuration = (StereoKind::Tetrahedral, value).into();
            checks.is_consistent(&candidate, &config).unwrap()
        })
        .collect();
    assert_eq!(actual, expected);
}

#[rstest]
#[case::prochiral(mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "H" "H"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]}"#), Stereogenicity::Prochiral)]
#[case::symmetric(mol_dsl_concrete!(r#"{:atoms ["C" "H" "H" "H" "H"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]}"#), Stereogenicity::Symmetric)]
fn test_stereo_assignment_checks_is_consistent_assignments(
    #[case] mut molecule: Molecule,
    #[case] class: Stereogenicity,
) {
    let site = StereoSite::Atom(AtomId(0));
    molecule
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .constraints
        .set(StereoAtomConstraintForm::Stereogenicity(
            StereogenicityForm::Lit(class),
        ));
    let checks = StereoAssignmentChecks::new(&molecule, &BTreeSet::from([site])).unwrap();
    let assignments = AbsoluteStereoAssignments::new([(
        site,
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined),
    )])
    .unwrap();
    let config = GraphSymmetryConfig {
        coloring: ConstitutionColoring::full(),
        iterate_to_fixpoint: false,
        max_iterations: 16,
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    };
    let actual: Vec<_> = assignments
        .filter_map(|assignment| {
            let mut candidate = molecule.clone();
            candidate
                .stereo_atom_mut(StereoAtomId(0))
                .attributes_mut()
                .configuration = assignment[0].1.clone();
            checks
                .is_consistent(&candidate, &config)
                .unwrap()
                .then_some(assignment)
        })
        .collect();
    assert_eq!(
        actual,
        vec![
            vec![(site, (StereoKind::Tetrahedral, 0).into())],
            vec![(site, (StereoKind::Tetrahedral, 1).into())]
        ]
    );
}

#[rstest]
#[case::and(true)]
#[case::or(false)]
fn test_stereo_assignment_checks_is_consistent_error(
    #[case] is_and: bool,
    #[values(false, true)] reverse: bool,
    #[values(None, Some(1u32))] source_assertion: Option<u32>,
) {
    let mut molecule = mol_dsl_concrete!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I" "N" "B"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :dative-bonds [{:donors [5] :acceptor 6 :attrs "1"}]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#
    );
    if let Some(value) = source_assertion {
        molecule
            .atom_mut(AtomId(0))
            .attributes_mut()
            .constraints
            .set(AtomConstraintForm::TetrahedralStereo(
                TetrahedralStereoForm::Stereo(value.into()),
            ));
    }
    let mut children = vec![
        Constraint::Atom(
            AtomId(0),
            AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Stereo(
                u32::from(is_and).into(),
            )),
        ),
        Constraint::DativeBond(
            DativeBondId(0),
            DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1),
        ),
    ];
    if reverse {
        children.reverse();
    }
    let mut editor = molecule.edit();
    editor.constraints_mut().push(if is_and {
        Constraint::And(children)
    } else {
        Constraint::Or(children)
    });
    let molecule = editor.finish().unwrap();
    let checks =
        StereoAssignmentChecks::new(&molecule, &BTreeSet::from([StereoSite::Atom(AtomId(0))]))
            .unwrap();
    assert_eq!(
        checks.is_consistent(
            &molecule,
            &GraphSymmetryConfig {
                coloring: ConstitutionColoring::full(),
                iterate_to_fixpoint: false,
                max_iterations: 16,
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
            }
        ),
        Err(
            ConstraintInvariantsError::DativeBondRingMembershipUnsupported {
                bond: DativeBondId(0)
            }
        )
    );
}
