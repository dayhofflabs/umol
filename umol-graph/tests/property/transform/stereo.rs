//! Stereo transformation laws on finite explicit Th/Ct sites and unrelated partial state.
//!
//! Addition and class-filtered removal use undetermined stored cosets, with ligand classes
//! specified independently by terminal element equality. Raw removal and enumeration
//! use mixed fixed/undetermined/finite domains. Exact expected entries check preservation and
//! reference transport; binary-mask enumeration checks soundness, completeness, order,
//! and multiplicity independently of production traversal. Renumbering also transports
//! implicit #T/#C assertion frames. Issued representative actions transport stored frames.
//! Collection comparisons preserve multiplicity and do not require a frame-invariant first item.
//! Frame transport uses nonempty absolute domains; the enumeration reference also covers empty
//! domains. Unsupported attribute states are explicit cases on generated valid frames.
//! Operations use Nauty, the default stereo model, and default ring algorithms.

mod reference;
mod strategies;

use std::collections::{BTreeMap, BTreeSet};

use proptest::prelude::*;
use proptest::test_runner::{Config, TestRunner};
use rstest::rstest;
use umol_graph::ops::model::StereoModel;
use umol_graph::ops::stereo::{
    StereoInconsistency, StereoPerception, StereoPerceptionError, StereoSelection, StereoSite,
};
use umol_graph::ops::transform::{
    AddStereoConfig, AddStereoError, EnumerateStereoConfig, EnumerateStereoError,
    RemoveStereoConfig, RemoveStereoError, RemoveStereoMode, StereoAdder, StereoEnumerator,
    StereoRemover, Transformer,
};
use umol_graph_core::{AutomorphismAlgorithm, GraphRemapping, Remapping};
use umol_graph_ir::ir::{
    AtomConstraintForm, AtomConstraintKey, AtomId, BondConstraintForm, BondConstraintKey, BondId,
    CisTransStereoForm, Constraint, FrameTransport, Molecule, MoleculeRemapping, NumForm, Reframe,
    StereoAtomForm, StereoAtomId, StereoBondForm, StereoBondId, StereoConfigurationForm,
    StereoCoset, StereoKind, StereoLigand, StereoLigandKind, StereoTerm, Stereogenicity,
    TetrahedralStereoForm,
};

use self::reference::{assignments, removed_entries, selected};
use self::strategies::{renumber, stereo_cases, Domain};

#[derive(Clone, Copy, Debug)]
enum ConfigurationCase {
    Kindless,
    Index,
    Finite,
    Variable,
    FiniteVariable,
    EmptyVariable,
    Relative,
    SquarePlanar,
    SiteKind,
}

#[rstest]
fn test_stereo_perception_classify_transport() {
    TestRunner::new(Config {
        failure_persistence: None,
        ..Config::default()
    })
    .run(&stereo_cases(Domain::Undetermined, 0..=4), |case| {
        let sites = case.sites.iter().map(|info| info.site).collect();
        let expected: BTreeMap<_, _> = case
            .sites
            .iter()
            .map(|info| (info.site, info.class))
            .collect();
        let perception = StereoPerception::new(&StereoModel::default());
        prop_assert_eq!(
            perception.classify(&case.source, &sites, AutomorphismAlgorithm::Nauty),
            Ok(expected.clone())
        );
        let reframed = case
            .source
            .clone()
            .reframe_by(&case.source.representative_action())
            .unwrap();
        prop_assert_eq!(
            perception.classify(&reframed, &sites, AutomorphismAlgorithm::Nauty),
            Ok(expected)
        );
        let renamed = renumber(case.source.clone(), &case);
        let renamed_expected: BTreeMap<_, _> = case
            .sites
            .iter()
            .map(|info| {
                (
                    match info.site {
                        StereoSite::Atom(atom) => StereoSite::Atom(case.images[atom.index()]),
                        StereoSite::Bond(bond) => StereoSite::Bond(bond),
                    },
                    info.class,
                )
            })
            .collect();
        let renamed_sites = renamed_expected.keys().copied().collect();
        prop_assert_eq!(
            perception.classify(&renamed, &renamed_sites, AutomorphismAlgorithm::Nauty),
            Ok(renamed_expected)
        );
        Ok(())
    })
    .unwrap();
}

#[rstest]
fn test_stereo_adder_transform() {
    TestRunner::new(Config {
        failure_persistence: None,
        ..Config::default()
    })
    .run(
        &(
            stereo_cases(Domain::Addition, 0..=4),
            any::<u8>(),
            0u8..4,
            prop::option::of(0u8..8),
        ),
        |(case, sites, kinds, classes)| {
            let selection = StereoSelection {
                sites: Some(
                    case.sites
                        .iter()
                        .enumerate()
                        .filter_map(|(i, info)| (sites & (1 << i) != 0).then_some(info.site))
                        .collect(),
                ),
                kinds: Some(
                    [StereoKind::Tetrahedral, StereoKind::CisTrans]
                        .into_iter()
                        .enumerate()
                        .filter_map(|(i, kind)| (kinds & (1 << i) != 0).then_some(kind))
                        .collect(),
                ),
                stereogenicities: classes.map(|mask| {
                    [
                        Stereogenicity::Symmetric,
                        Stereogenicity::Prochiral,
                        Stereogenicity::Stereogenic,
                    ]
                    .into_iter()
                    .enumerate()
                    .filter_map(|(i, class)| (mask & (1 << i) != 0).then_some(class))
                    .collect()
                }),
            };
            let transformer = StereoAdder::with_config(
                &StereoModel::default(),
                AddStereoConfig {
                    automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                    selection,
                    ..Default::default()
                },
            );
            let mut entries = case.entries.clone();
            for (index, info) in case.sites.iter().enumerate() {
                if info.record.is_some() || !selected(info, index, sites, kinds, classes) {
                    continue;
                }
                let ligands = info
                    .natural
                    .map(|atom| StereoLigand::new(atom, StereoLigandKind::Atom))
                    .to_vec();
                match info.site {
                    StereoSite::Atom(atom) => entries.stereo_atoms.push((
                        atom,
                        ligands,
                        StereoAtomForm::new(info.kind, StereoCoset::Undetermined),
                    )),
                    StereoSite::Bond(bond) => entries.stereo_bonds.push((
                        bond,
                        ligands,
                        StereoBondForm::new(info.kind, StereoCoset::Undetermined),
                    )),
                }
            }
            let expected = Molecule::from_entries(entries);
            let original = case.source.clone();
            prop_assert_eq!(
                transformer.transform(case.source.clone()),
                Ok(expected.clone())
            );
            let mut in_place = case.source.clone();
            prop_assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
            prop_assert_eq!(&in_place, &expected);
            prop_assert_eq!(
                transformer.transform(expected.clone()),
                Ok(expected.clone())
            );
            prop_assert_eq!(
                transformer.transform_iter(&case.source).collect::<Vec<_>>(),
                vec![expected.clone()]
            );
            prop_assert_eq!(expected.clone().edit().finish(), Ok(expected));
            prop_assert_eq!(case.source, original);
            Ok(())
        },
    )
    .unwrap();
}

#[rstest]
fn test_stereo_adder_transform_relabeling() {
    TestRunner::new(Config {
        failure_persistence: None,
        ..Config::default()
    })
    .run(
        &(
            stereo_cases(Domain::Addition, 0..=4),
            any::<u8>(),
            0u8..4,
            prop::option::of(0u8..8),
        ),
        |(case, sites, kinds, classes)| {
            let selection = StereoSelection {
                sites: Some(
                    case.sites
                        .iter()
                        .enumerate()
                        .filter_map(|(i, info)| (sites & (1 << i) != 0).then_some(info.site))
                        .collect(),
                ),
                kinds: Some(
                    [StereoKind::Tetrahedral, StereoKind::CisTrans]
                        .into_iter()
                        .enumerate()
                        .filter_map(|(i, kind)| (kinds & (1 << i) != 0).then_some(kind))
                        .collect(),
                ),
                stereogenicities: classes.map(|mask| {
                    [
                        Stereogenicity::Symmetric,
                        Stereogenicity::Prochiral,
                        Stereogenicity::Stereogenic,
                    ]
                    .into_iter()
                    .enumerate()
                    .filter_map(|(i, class)| (mask & (1 << i) != 0).then_some(class))
                    .collect()
                }),
            };
            let mut mapped_selection = selection.clone();
            mapped_selection.sites = selection.sites.as_ref().map(|sites| {
                sites
                    .iter()
                    .map(|site| match site {
                        StereoSite::Atom(atom) => StereoSite::Atom(case.images[atom.index()]),
                        StereoSite::Bond(bond) => StereoSite::Bond(*bond),
                    })
                    .collect()
            });
            let left = StereoAdder::with_config(
                &StereoModel::default(),
                AddStereoConfig {
                    automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                    selection,
                    ..Default::default()
                },
            )
            .transform(case.source.clone())
            .unwrap();
            let mapped_source = renumber(case.source.clone(), &case);
            let right = StereoAdder::with_config(
                &StereoModel::default(),
                AddStereoConfig {
                    automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                    selection: mapped_selection,
                    ..Default::default()
                },
            )
            .transform(mapped_source)
            .unwrap();
            let transported = renumber(left, &case);
            let alignment = MoleculeRemapping::new(
                GraphRemapping::identity(right.atoms().count(), right.bonds().count()),
                Remapping::identity(0),
                Remapping::identity(0),
                Remapping::identity(0),
                Remapping::identity(0),
                Remapping::new(
                    transported
                        .stereo_atoms()
                        .iter()
                        .map(|stereo| right.stereo_atoms().at(stereo.site_id()).unwrap().id())
                        .collect(),
                )
                .unwrap(),
                Remapping::new(
                    transported
                        .stereo_bonds()
                        .iter()
                        .map(|stereo| right.stereo_bonds().at(stereo.site_id()).unwrap().id())
                        .collect(),
                )
                .unwrap(),
            );
            prop_assert!(transported.framed_eq_under(&right, &alignment));
            Ok(())
        },
    )
    .unwrap();
}

#[rstest]
fn test_stereo_adder_transform_evidence() {
    TestRunner::new(Config {
        failure_persistence: None,
        ..Config::default()
    })
    .run(
        &(stereo_cases(Domain::Undetermined, 1..=1), any::<bool>()),
        |(case, value)| {
            let value = u32::from(value);
            let info = &case.sites[0];
            let asserted = value ^ info.parity ^ 1;
            let mut source = case.source.clone();
            let inconsistency = match info.site {
                StereoSite::Atom(atom) => {
                    source
                        .stereo_atom_mut(StereoAtomId(0))
                        .attributes_mut()
                        .configuration = (info.kind, value).into();
                    source.atom_mut(atom).attributes_mut().constraints.set(
                        AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Stereo(
                            asserted.into(),
                        )),
                    );
                    StereoInconsistency::TetrahedralStereoMismatch {
                        atom,
                        stereo_atom: StereoAtomId(0),
                    }
                }
                StereoSite::Bond(bond) => {
                    source
                        .stereo_bond_mut(StereoBondId(0))
                        .attributes_mut()
                        .configuration = (info.kind, value).into();
                    source.bond_mut(bond).attributes_mut().constraints.set(
                        BondConstraintForm::CisTransStereo(CisTransStereoForm::Stereo(
                            asserted.into(),
                        )),
                    );
                    StereoInconsistency::CisTransStereoMismatch {
                        bond,
                        stereo_bond: StereoBondId(0),
                    }
                }
            };
            let transformer = StereoAdder::with_config(
                &StereoModel::default(),
                AddStereoConfig {
                    automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                    selection: StereoSelection {
                        sites: Some(BTreeSet::from([info.site])),
                        kinds: Some(BTreeSet::from([info.kind])),
                        stereogenicities: None,
                    },
                    ..Default::default()
                },
            );
            let error = AddStereoError::Inconsistency(inconsistency);
            let original = source.clone();
            let mut input = source.clone();
            prop_assert_eq!(transformer.transform(source.clone()), Err(error.clone()));
            prop_assert_eq!(transformer.transform_into(&mut input), Err(error));
            prop_assert_eq!(input, original.clone());
            prop_assert_eq!(
                transformer.transform_iter(&source).collect::<Vec<_>>(),
                Vec::<Molecule>::new()
            );
            prop_assert_eq!(source, original);
            Ok(())
        },
    )
    .unwrap();
}

#[rstest]
#[case::all(RemoveStereoMode::All)]
#[case::undetermined(RemoveStereoMode::Undetermined)]
fn test_stereo_remover_transform(#[case] mode: RemoveStereoMode) {
    TestRunner::new(Config {
        failure_persistence: None,
        ..Config::default()
    })
    .run(
        &(stereo_cases(Domain::Absolute, 0..=4), any::<u8>(), 0u8..4),
        |(case, sites, kinds)| {
            let selection = StereoSelection {
                sites: Some(
                    case.sites
                        .iter()
                        .enumerate()
                        .filter_map(|(i, info)| (sites & (1 << i) != 0).then_some(info.site))
                        .collect(),
                ),
                kinds: Some(
                    [StereoKind::Tetrahedral, StereoKind::CisTrans]
                        .into_iter()
                        .enumerate()
                        .filter_map(|(i, kind)| (kinds & (1 << i) != 0).then_some(kind))
                        .collect(),
                ),
                stereogenicities: None,
            };
            let removed: BTreeSet<_> = case
                .sites
                .iter()
                .enumerate()
                .filter_map(|(i, info)| {
                    (selected(info, i, sites, kinds, None)
                        && (mode == RemoveStereoMode::All
                            || info.coset == StereoCoset::Undetermined))
                        .then_some(info.site)
                })
                .collect();
            let expected =
                Molecule::from_entries(removed_entries(&case.entries, &case.sites, &removed));
            let transformer = StereoRemover::with_config(
                &StereoModel::default(),
                RemoveStereoConfig {
                    automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                    selection,
                    mode,
                },
            );
            let original = case.source.clone();
            prop_assert_eq!(
                transformer.transform(case.source.clone()),
                Ok(expected.clone())
            );
            let mut in_place = case.source.clone();
            prop_assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
            prop_assert_eq!(in_place, expected.clone());
            prop_assert_eq!(
                transformer.transform(expected.clone()),
                Ok(expected.clone())
            );
            prop_assert_eq!(
                transformer.transform_iter(&case.source).collect::<Vec<_>>(),
                vec![expected.clone()]
            );
            prop_assert_eq!(expected.clone().edit().finish(), Ok(expected));
            prop_assert_eq!(case.source, original);
            Ok(())
        },
    )
    .unwrap();
}

#[rstest]
fn test_stereo_remover_transform_classification() {
    TestRunner::new(Config {
        failure_persistence: None,
        ..Config::default()
    })
    .run(
        &(
            stereo_cases(Domain::Undetermined, 0..=4),
            any::<u8>(),
            0u8..4,
            0u8..8,
        ),
        |(case, sites, kinds, classes)| {
            let selection = StereoSelection {
                sites: Some(
                    case.sites
                        .iter()
                        .enumerate()
                        .filter_map(|(i, info)| (sites & (1 << i) != 0).then_some(info.site))
                        .collect(),
                ),
                kinds: Some(
                    [StereoKind::Tetrahedral, StereoKind::CisTrans]
                        .into_iter()
                        .enumerate()
                        .filter_map(|(i, kind)| (kinds & (1 << i) != 0).then_some(kind))
                        .collect(),
                ),
                stereogenicities: Some(
                    [
                        Stereogenicity::Symmetric,
                        Stereogenicity::Prochiral,
                        Stereogenicity::Stereogenic,
                    ]
                    .into_iter()
                    .enumerate()
                    .filter_map(|(i, class)| (classes & (1 << i) != 0).then_some(class))
                    .collect(),
                ),
            };
            let transformer = StereoRemover::with_config(
                &StereoModel::default(),
                RemoveStereoConfig {
                    automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                    selection,
                    ..Default::default()
                },
            );
            let result = transformer.transform(case.source.clone()).unwrap();
            let removed: BTreeSet<_> = case
                .sites
                .iter()
                .enumerate()
                .filter_map(|(i, info)| {
                    selected(info, i, sites, kinds, Some(classes)).then_some(info.site)
                })
                .collect();
            let expected =
                Molecule::from_entries(removed_entries(&case.entries, &case.sites, &removed));
            prop_assert_eq!(&result, &expected);
            let mut in_place = case.source.clone();
            prop_assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
            prop_assert_eq!(&in_place, &result);
            prop_assert_eq!(
                transformer.transform_iter(&case.source).collect::<Vec<_>>(),
                vec![result.clone()]
            );
            prop_assert_eq!(result.clone().edit().finish(), Ok(result));
            Ok(())
        },
    )
    .unwrap();
}

#[rstest]
#[case::all(RemoveStereoMode::All)]
#[case::undetermined(RemoveStereoMode::Undetermined)]
fn test_stereo_remover_transform_relabeling(#[case] mode: RemoveStereoMode) {
    TestRunner::new(Config {
        failure_persistence: None,
        ..Config::default()
    })
    .run(
        &(stereo_cases(Domain::Absolute, 0..=4), any::<u8>(), 0u8..4),
        |(case, sites, kinds)| {
            let selection = StereoSelection {
                sites: Some(
                    case.sites
                        .iter()
                        .enumerate()
                        .filter_map(|(i, info)| (sites & (1 << i) != 0).then_some(info.site))
                        .collect(),
                ),
                kinds: Some(
                    [StereoKind::Tetrahedral, StereoKind::CisTrans]
                        .into_iter()
                        .enumerate()
                        .filter_map(|(i, kind)| (kinds & (1 << i) != 0).then_some(kind))
                        .collect(),
                ),
                stereogenicities: None,
            };
            let mut mapped_selection = selection.clone();
            mapped_selection.sites = selection.sites.as_ref().map(|sites| {
                sites
                    .iter()
                    .map(|site| match site {
                        StereoSite::Atom(atom) => StereoSite::Atom(case.images[atom.index()]),
                        StereoSite::Bond(bond) => StereoSite::Bond(*bond),
                    })
                    .collect()
            });
            let left = StereoRemover::with_config(
                &StereoModel::default(),
                RemoveStereoConfig {
                    automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                    selection,
                    mode,
                },
            )
            .transform(case.source.clone())
            .unwrap();
            let mapped_source = renumber(case.source.clone(), &case);
            let right = StereoRemover::with_config(
                &StereoModel::default(),
                RemoveStereoConfig {
                    automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                    selection: mapped_selection,
                    mode,
                },
            )
            .transform(mapped_source)
            .unwrap();
            prop_assert_eq!(renumber(left, &case), right);
            Ok(())
        },
    )
    .unwrap();
}

#[rstest]
fn test_stereo_remover_transform_frames() {
    TestRunner::new(Config {
        failure_persistence: None,
        ..Config::default()
    })
    .run(
        &(
            stereo_cases(Domain::Undetermined, 0..=4),
            any::<u8>(),
            0u8..4,
            0u8..8,
        ),
        |(case, sites, kinds, classes)| {
            let selection = StereoSelection {
                sites: Some(
                    case.sites
                        .iter()
                        .enumerate()
                        .filter_map(|(i, info)| (sites & (1 << i) != 0).then_some(info.site))
                        .collect(),
                ),
                kinds: Some(
                    [StereoKind::Tetrahedral, StereoKind::CisTrans]
                        .into_iter()
                        .enumerate()
                        .filter_map(|(i, kind)| (kinds & (1 << i) != 0).then_some(kind))
                        .collect(),
                ),
                stereogenicities: Some(
                    [
                        Stereogenicity::Symmetric,
                        Stereogenicity::Prochiral,
                        Stereogenicity::Stereogenic,
                    ]
                    .into_iter()
                    .enumerate()
                    .filter_map(|(i, class)| (classes & (1 << i) != 0).then_some(class))
                    .collect(),
                ),
            };
            let transformer = StereoRemover::with_config(
                &StereoModel::default(),
                RemoveStereoConfig {
                    automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                    selection,
                    ..Default::default()
                },
            );
            let target = case
                .source
                .clone()
                .reframe_by(&case.source.representative_action())
                .unwrap();
            let left = transformer.transform(case.source).unwrap();
            let right = transformer.transform(target).unwrap();
            prop_assert!(left.framed_eq(&right));
            Ok(())
        },
    )
    .unwrap();
}

#[rstest]
fn test_stereo_enumerator_transform() {
    TestRunner::new(Config {
        failure_persistence: None,
        ..Config::default()
    })
    .run(
        &(
            stereo_cases(Domain::Absolute, 0..=4),
            any::<u8>(),
            0u8..4,
            any::<bool>(),
        ),
        |(mut case, sites, kinds, constrain)| {
            let selection = StereoSelection {
                sites: Some(
                    case.sites
                        .iter()
                        .enumerate()
                        .filter_map(|(i, info)| (sites & (1 << i) != 0).then_some(info.site))
                        .collect(),
                ),
                kinds: Some(
                    [StereoKind::Tetrahedral, StereoKind::CisTrans]
                        .into_iter()
                        .enumerate()
                        .filter_map(|(i, kind)| (kinds & (1 << i) != 0).then_some(kind))
                        .collect(),
                ),
                stereogenicities: None,
            };
            let chosen: BTreeSet<_> = case
                .sites
                .iter()
                .enumerate()
                .filter_map(|(i, info)| selected(info, i, sites, kinds, None).then_some(i))
                .collect();
            let joint = if constrain && chosen.len() >= 2 {
                let mut indices = chosen.iter().copied();
                let pair = (indices.next().unwrap(), indices.next().unwrap());
                let alternatives: Vec<_> = [0u32, 1]
                    .into_iter()
                    .map(|value| {
                        Constraint::And(
                            [(pair.0, value), (pair.1, value ^ 1)]
                                .into_iter()
                                .map(|(i, value)| match case.sites[i].site {
                                    StereoSite::Atom(atom) => Constraint::Atom(
                                        atom,
                                        AtomConstraintForm::TetrahedralStereo(
                                            TetrahedralStereoForm::Stereo(value.into()),
                                        ),
                                    ),
                                    StereoSite::Bond(bond) => Constraint::Bond(
                                        bond,
                                        BondConstraintForm::CisTransStereo(
                                            CisTransStereoForm::Stereo(value.into()),
                                        ),
                                    ),
                                })
                                .collect(),
                        )
                    })
                    .collect();
                case.entries.constraints.push(Constraint::Or(alternatives));
                case.source = Molecule::from_entries(case.entries.clone());
                Some(pair)
            } else {
                None
            };
            let expected: Vec<_> = assignments(&case.sites, &chosen, &case.order, joint)
                .into_iter()
                .map(|values| {
                    let mut entries = case.entries.clone();
                    for (info, coset) in case.sites.iter().zip(values) {
                        let configuration = StereoConfigurationForm::Kinded(info.kind, coset);
                        match info.site {
                            StereoSite::Atom(_) => {
                                entries.stereo_atoms[info.record.unwrap()].2.configuration =
                                    configuration
                            }
                            StereoSite::Bond(_) => {
                                entries.stereo_bonds[info.record.unwrap()].2.configuration =
                                    configuration
                            }
                        }
                    }
                    Molecule::from_entries(entries)
                })
                .collect();
            let transformer = StereoEnumerator::with_config(
                &StereoModel::default(),
                EnumerateStereoConfig {
                    automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                    selection,
                },
                case.order.iter().map(|&i| case.sites[i].site).collect(),
            );
            let original = case.source.clone();
            let actual: Vec<_> = transformer.transform_iter(&case.source).collect();
            prop_assert_eq!(&actual, &expected);
            let mut in_place = case.source.clone();
            if let Some(first) = expected.first() {
                prop_assert_eq!(
                    transformer.transform(case.source.clone()),
                    Ok(first.clone())
                );
                prop_assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
                prop_assert_eq!(&in_place, first);
            } else {
                prop_assert_eq!(
                    transformer.transform(case.source.clone()),
                    Err(EnumerateStereoError::NoAssignment)
                );
                prop_assert_eq!(
                    transformer.transform_into(&mut in_place),
                    Err(EnumerateStereoError::NoAssignment)
                );
                prop_assert_eq!(&in_place, &original);
            }
            for output in actual {
                prop_assert_eq!(output.clone().edit().finish(), Ok(output));
            }
            prop_assert_eq!(case.source, original);
            Ok(())
        },
    )
    .unwrap();
}

#[rstest]
#[case::relabel(false)]
#[case::frame(true)]
fn test_stereo_enumerator_transform_transport(#[case] frame: bool) {
    TestRunner::new(Config {
        failure_persistence: None,
        ..Config::default()
    })
    .run(
        &(
            stereo_cases(
                if frame {
                    Domain::AbsoluteNonempty
                } else {
                    Domain::Absolute
                },
                0..=4,
            ),
            any::<u8>(),
            0u8..4,
        ),
        |(case, sites, kinds)| {
            let selection = StereoSelection {
                sites: Some(
                    case.sites
                        .iter()
                        .enumerate()
                        .filter_map(|(i, info)| (sites & (1 << i) != 0).then_some(info.site))
                        .collect(),
                ),
                kinds: Some(
                    [StereoKind::Tetrahedral, StereoKind::CisTrans]
                        .into_iter()
                        .enumerate()
                        .filter_map(|(i, kind)| (kinds & (1 << i) != 0).then_some(kind))
                        .collect(),
                ),
                stereogenicities: None,
            };
            let order: Vec<_> = case.order.iter().map(|&i| case.sites[i].site).collect();
            let actions = case.source.representative_action();
            let target = if frame {
                case.source.clone().reframe_by(&actions).unwrap()
            } else {
                renumber(case.source.clone(), &case)
            };
            let mut target_selection = selection.clone();
            if !frame {
                target_selection.sites = selection.sites.as_ref().map(|sites| {
                    sites
                        .iter()
                        .map(|site| match site {
                            StereoSite::Atom(atom) => StereoSite::Atom(case.images[atom.index()]),
                            StereoSite::Bond(bond) => StereoSite::Bond(*bond),
                        })
                        .collect()
                });
            }
            let target_order = if frame {
                order.clone()
            } else {
                order
                    .iter()
                    .map(|site| match site {
                        StereoSite::Atom(atom) => StereoSite::Atom(case.images[atom.index()]),
                        StereoSite::Bond(bond) => StereoSite::Bond(*bond),
                    })
                    .collect()
            };
            let transformer = StereoEnumerator::with_config(
                &StereoModel::default(),
                EnumerateStereoConfig {
                    automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                    selection,
                },
                order,
            );
            let other = StereoEnumerator::with_config(
                &StereoModel::default(),
                EnumerateStereoConfig {
                    automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                    selection: target_selection,
                },
                target_order,
            );
            let expected: Vec<_> = transformer
                .transform_iter(&case.source)
                .map(|result| {
                    if frame {
                        result.reframe_by(&actions).unwrap()
                    } else {
                        renumber(result, &case)
                    }
                })
                .collect();
            let mut actual: Vec<_> = other.transform_iter(&target).collect();
            for result in expected {
                let Some(position) = actual.iter().position(|candidate| *candidate == result)
                else {
                    return Err(TestCaseError::fail("transported assignment is absent"));
                };
                actual.remove(position);
            }
            prop_assert_eq!(actual, Vec::<Molecule>::new());
            Ok(())
        },
    )
    .unwrap();
}

#[rstest]
fn test_stereo_enumerator_transform_iter_independence() {
    TestRunner::new(Config {
        failure_persistence: None,
        ..Config::default()
    })
    .run(&stereo_cases(Domain::Undetermined, 2..=4), |case| {
        let mut entries = case.entries.clone();
        for info in &case.sites {
            match info.site {
                StereoSite::Atom(atom) => {
                    entries.atoms[atom.index()]
                        .constraints
                        .remove(AtomConstraintKey::TetrahedralStereo);
                }
                StereoSite::Bond(bond) => {
                    entries.bonds[bond.index()]
                        .2
                        .constraints
                        .remove(BondConstraintKey::CisTransStereo);
                }
            }
        }
        let mut source = Molecule::from_entries(entries.clone());
        let transformer = StereoEnumerator::with_config(
            &StereoModel::default(),
            EnumerateStereoConfig {
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                ..Default::default()
            },
            case.order.iter().map(|&i| case.sites[i].site).collect(),
        );
        let mut run = transformer.transform_iter(&source);
        let noise = AtomId((source.atoms().count() - 1) as u32);
        source.atom_mut(noise).attributes_mut().charge = NumForm::Lit(-1);
        let mut changed_source = entries.clone();
        changed_source.atoms[noise.index()].charge = NumForm::Lit(-1);
        for info in &case.sites {
            let configuration = (info.kind, 0).into();
            match info.site {
                StereoSite::Atom(_) => {
                    entries.stereo_atoms[info.record.unwrap()].2.configuration = configuration
                }
                StereoSite::Bond(_) => {
                    entries.stereo_bonds[info.record.unwrap()].2.configuration = configuration
                }
            }
        }
        let mut first = run.next().unwrap();
        prop_assert_eq!(&first, &Molecule::from_entries(entries.clone()));
        first.atom_mut(noise).attributes_mut().charge = NumForm::Lit(1);
        let info = &case.sites[case.order[0]];
        let mut changed_first = entries.clone();
        changed_first.atoms[noise.index()].charge = NumForm::Lit(1);
        match info.site {
            StereoSite::Atom(_) => {
                first
                    .stereo_atom_mut(StereoAtomId(info.record.unwrap() as u32))
                    .attributes_mut()
                    .configuration = (info.kind, 1).into();
                changed_first.stereo_atoms[info.record.unwrap()]
                    .2
                    .configuration = (info.kind, 1).into();
            }
            StereoSite::Bond(_) => {
                first
                    .stereo_bond_mut(StereoBondId(info.record.unwrap() as u32))
                    .attributes_mut()
                    .configuration = (info.kind, 1).into();
                changed_first.stereo_bonds[info.record.unwrap()]
                    .2
                    .configuration = (info.kind, 1).into();
            }
        }
        let last = &case.sites[*case.order.last().unwrap()];
        match last.site {
            StereoSite::Atom(_) => {
                entries.stereo_atoms[last.record.unwrap()].2.configuration = (last.kind, 1).into()
            }
            StereoSite::Bond(_) => {
                entries.stereo_bonds[last.record.unwrap()].2.configuration = (last.kind, 1).into()
            }
        }
        prop_assert_eq!(run.next(), Some(Molecule::from_entries(entries)));
        prop_assert_eq!(first, Molecule::from_entries(changed_first));
        prop_assert_eq!(source, Molecule::from_entries(changed_source));
        Ok(())
    })
    .unwrap();
}

#[rstest]
#[case::atom(true)]
#[case::bond(false)]
fn test_transformer_transform_selection(#[case] atom: bool) {
    TestRunner::new(Config {
        failure_persistence: None,
        ..Config::default()
    })
    .run(
        &(stereo_cases(Domain::Undetermined, 0..=4), 0u32..8),
        |(case, offset)| {
            let site = if atom {
                StereoSite::Atom(AtomId(case.source.atoms().count() as u32 + offset))
            } else {
                StereoSite::Bond(BondId(case.source.bonds().count() as u32 + offset))
            };
            let selection = StereoSelection {
                sites: Some(BTreeSet::from([site])),
                kinds: Some(BTreeSet::new()),
                stereogenicities: Some(BTreeSet::new()),
            };
            let cause = StereoPerceptionError::InvalidSite { site };
            let adder = StereoAdder::with_config(
                &StereoModel::default(),
                AddStereoConfig {
                    automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                    selection: selection.clone(),
                    ..Default::default()
                },
            );
            let remover = StereoRemover::with_config(
                &StereoModel::default(),
                RemoveStereoConfig {
                    automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                    selection: selection.clone(),
                    ..Default::default()
                },
            );
            let enumerator = StereoEnumerator::with_config(
                &StereoModel::default(),
                EnumerateStereoConfig {
                    automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                    selection,
                },
                Vec::new(),
            );
            let original = case.source.clone();
            let mut input = case.source.clone();
            prop_assert_eq!(
                adder.transform(input.clone()),
                Err(AddStereoError::Perception(cause.clone()))
            );
            prop_assert_eq!(
                adder.transform_into(&mut input),
                Err(AddStereoError::Perception(cause.clone()))
            );
            prop_assert_eq!(&input, &original);
            prop_assert_eq!(
                adder.transform_iter(&input).collect::<Vec<_>>(),
                Vec::<Molecule>::new()
            );
            prop_assert_eq!(
                remover.transform(input.clone()),
                Err(RemoveStereoError::Perception(cause.clone()))
            );
            prop_assert_eq!(
                remover.transform_into(&mut input),
                Err(RemoveStereoError::Perception(cause.clone()))
            );
            prop_assert_eq!(&input, &original);
            prop_assert_eq!(
                remover.transform_iter(&input).collect::<Vec<_>>(),
                Vec::<Molecule>::new()
            );
            prop_assert_eq!(
                enumerator.transform(input.clone()),
                Err(EnumerateStereoError::Perception(cause.clone()))
            );
            prop_assert_eq!(
                enumerator.transform_into(&mut input),
                Err(EnumerateStereoError::Perception(cause))
            );
            prop_assert_eq!(
                enumerator.transform_iter(&input).collect::<Vec<_>>(),
                Vec::<Molecule>::new()
            );
            prop_assert_eq!(&input, &original);
            prop_assert_eq!(case.source, original);
            Ok(())
        },
    )
    .unwrap();
}

#[rstest]
#[case::kindless(ConfigurationCase::Kindless)]
#[case::index(ConfigurationCase::Index)]
#[case::finite(ConfigurationCase::Finite)]
#[case::variable(ConfigurationCase::Variable)]
#[case::finite_variable(ConfigurationCase::FiniteVariable)]
#[case::empty_variable(ConfigurationCase::EmptyVariable)]
#[case::relative(ConfigurationCase::Relative)]
#[case::square_planar(ConfigurationCase::SquarePlanar)]
#[case::site_kind(ConfigurationCase::SiteKind)]
fn test_transformer_transform_configuration(#[case] variant: ConfigurationCase) {
    TestRunner::new(Config {
        failure_persistence: None,
        ..Config::default()
    })
    .run(&stereo_cases(Domain::Undetermined, 1..=1), |case| {
        let info = &case.sites[0];
        let form = match variant {
            ConfigurationCase::Kindless => StereoConfigurationForm::Undetermined,
            ConfigurationCase::Index => (info.kind, 2).into(),
            ConfigurationCase::Finite => {
                StereoConfigurationForm::Kinded(info.kind, StereoCoset::lit_set([0, 2]))
            }
            ConfigurationCase::Variable => StereoConfigurationForm::Kinded(
                info.kind,
                StereoCoset::term(StereoTerm::var("configuration")),
            ),
            ConfigurationCase::FiniteVariable => StereoConfigurationForm::Kinded(
                info.kind,
                StereoCoset::term(StereoTerm::var_in("configuration", [0, 1])),
            ),
            ConfigurationCase::EmptyVariable => StereoConfigurationForm::Kinded(
                info.kind,
                StereoCoset::term(StereoTerm::var_in("configuration", [])),
            ),
            ConfigurationCase::Relative => StereoConfigurationForm::Kinded(
                info.kind,
                StereoCoset::term(StereoTerm::swap(StereoTerm::var("configuration"))),
            ),
            ConfigurationCase::SquarePlanar => (StereoKind::SquarePlanar, 0).into(),
            ConfigurationCase::SiteKind => (
                if info.kind == StereoKind::Tetrahedral {
                    StereoKind::CisTrans
                } else {
                    StereoKind::Tetrahedral
                },
                0,
            )
                .into(),
        };
        let mut source = case.source.clone();
        match info.site {
            StereoSite::Atom(_) => {
                source
                    .stereo_atom_mut(StereoAtomId(0))
                    .attributes_mut()
                    .configuration = form
            }
            StereoSite::Bond(_) => {
                source
                    .stereo_bond_mut(StereoBondId(0))
                    .attributes_mut()
                    .configuration = form
            }
        }
        let original = source.clone();
        let cause = if matches!(variant, ConfigurationCase::SquarePlanar) {
            StereoPerceptionError::UnsupportedKind {
                kind: StereoKind::SquarePlanar,
            }
        } else {
            StereoPerceptionError::UnsupportedConfiguration { site: info.site }
        };
        let error = EnumerateStereoError::Perception(cause.clone());
        prop_assert_eq!(
            StereoPerception::new(&StereoModel::default()).classify(
                &source,
                &BTreeSet::from([info.site]),
                AutomorphismAlgorithm::Nauty
            ),
            Err(cause.clone())
        );
        let selection = StereoSelection {
            sites: Some(BTreeSet::from([info.site])),
            ..Default::default()
        };
        let transformer = StereoEnumerator::with_config(
            &StereoModel::default(),
            EnumerateStereoConfig {
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                selection: selection.clone(),
            },
            vec![info.site],
        );
        let mut input = source.clone();
        prop_assert_eq!(transformer.transform(source.clone()), Err(error.clone()));
        prop_assert_eq!(transformer.transform_into(&mut input), Err(error));
        prop_assert_eq!(&input, &original);
        prop_assert_eq!(
            transformer.transform_iter(&source).collect::<Vec<_>>(),
            Vec::<Molecule>::new()
        );
        let classified = StereoRemover::with_config(
            &StereoModel::default(),
            RemoveStereoConfig {
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                selection: StereoSelection {
                    stereogenicities: Some(BTreeSet::from([Stereogenicity::Symmetric])),
                    ..selection.clone()
                },
                ..Default::default()
            },
        );
        prop_assert_eq!(
            classified.transform_into(&mut input),
            Err(RemoveStereoError::Perception(cause.clone()))
        );
        prop_assert_eq!(
            classified.transform(source.clone()),
            Err(RemoveStereoError::Perception(cause))
        );
        prop_assert_eq!(
            classified.transform_iter(&source).collect::<Vec<_>>(),
            Vec::<Molecule>::new()
        );
        prop_assert_eq!(&input, &original);
        let undetermined = StereoRemover::with_config(
            &StereoModel::default(),
            RemoveStereoConfig {
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                selection: selection.clone(),
                mode: RemoveStereoMode::Undetermined,
            },
        );
        prop_assert_eq!(undetermined.transform(source.clone()), Ok(source.clone()));
        let all = StereoRemover::with_config(
            &StereoModel::default(),
            RemoveStereoConfig {
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                selection,
                ..Default::default()
            },
        );
        let expected = Molecule::from_entries(removed_entries(
            &case.entries,
            &case.sites,
            &BTreeSet::from([info.site]),
        ));
        prop_assert_eq!(all.transform(source.clone()), Ok(expected.clone()));
        prop_assert_eq!(all.transform_into(&mut input), Ok(()));
        prop_assert_eq!(input, expected);
        prop_assert_eq!(source, original);
        Ok(())
    })
    .unwrap();
}

#[rstest]
#[case::addition(false)]
#[case::removal(true)]
fn test_transformer_transform_iter_independence(#[case] removal: bool) {
    TestRunner::new(Config {
        failure_persistence: None,
        ..Config::default()
    })
    .run(
        &stereo_cases(
            if removal {
                Domain::Undetermined
            } else {
                Domain::Addition
            },
            2..=4,
        ),
        |case| {
            let (mut first, second, expected_entries) = if removal {
                let removed = BTreeSet::from([case.sites[0].site]);
                let transformer = StereoRemover::with_config(
                    &StereoModel::default(),
                    RemoveStereoConfig {
                        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                        selection: StereoSelection {
                            sites: Some(removed.clone()),
                            ..Default::default()
                        },
                        ..Default::default()
                    },
                );
                let first_result = transformer.transform_iter(&case.source).next().unwrap();
                let second_result = transformer.transform_iter(&case.source).next().unwrap();
                (
                    first_result,
                    second_result,
                    removed_entries(&case.entries, &case.sites, &removed),
                )
            } else {
                let transformer = StereoAdder::with_config(
                    &StereoModel::default(),
                    AddStereoConfig {
                        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                        selection: StereoSelection {
                            sites: Some(case.sites.iter().map(|info| info.site).collect()),
                            ..Default::default()
                        },
                        ..Default::default()
                    },
                );
                let mut entries = case.entries.clone();
                for info in case.sites.iter().filter(|info| info.record.is_none()) {
                    let ligands = info
                        .natural
                        .map(|id| StereoLigand::new(id, StereoLigandKind::Atom))
                        .to_vec();
                    match info.site {
                        StereoSite::Atom(atom) => entries.stereo_atoms.push((
                            atom,
                            ligands,
                            StereoAtomForm::new(info.kind, StereoCoset::Undetermined),
                        )),
                        StereoSite::Bond(bond) => entries.stereo_bonds.push((
                            bond,
                            ligands,
                            StereoBondForm::new(info.kind, StereoCoset::Undetermined),
                        )),
                    }
                }
                let first_result = transformer.transform_iter(&case.source).next().unwrap();
                let second_result = transformer.transform_iter(&case.source).next().unwrap();
                (first_result, second_result, entries)
            };
            let original = case.source.clone();
            let expected = Molecule::from_entries(expected_entries.clone());
            prop_assert_eq!(&first, &expected);
            let noise = AtomId((case.source.atoms().count() - 1) as u32);
            first.atom_mut(noise).attributes_mut().charge = NumForm::Lit(1);
            let mut changed = expected_entries;
            changed.atoms[noise.index()].charge = NumForm::Lit(1);
            if !changed.stereo_atoms.is_empty() {
                first
                    .stereo_atom_mut(StereoAtomId(0))
                    .attributes_mut()
                    .configuration = (StereoKind::Tetrahedral, 1).into();
                changed.stereo_atoms[0].2.configuration = (StereoKind::Tetrahedral, 1).into();
            } else {
                first
                    .stereo_bond_mut(StereoBondId(0))
                    .attributes_mut()
                    .configuration = (StereoKind::CisTrans, 1).into();
                changed.stereo_bonds[0].2.configuration = (StereoKind::CisTrans, 1).into();
            }
            prop_assert_eq!(first, Molecule::from_entries(changed));
            prop_assert_eq!(second, expected);
            prop_assert_eq!(case.source, original);
            Ok(())
        },
    )
    .unwrap();
}
