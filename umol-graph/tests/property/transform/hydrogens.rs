//! Hydrogen transformation laws over independently constructed representations.
//!
//! Literal and partial chains mix implicit counts, ordinary explicit hydrogens,
//! protected hydrogens, and unrelated overlays. Expected products are constructed
//! from the generation parameters, then transported under arbitrary atom orders.
//! Conservation and retained-atom readings supply additional checks independent of
//! structural equality. Stereo cases vary both atom numbering and stored frames;
//! rejection cases are generated separately from successful partial inputs.
//! Reconstruction uses literal counts on known non-hydrogen bearers, literal
//! localized assertions, compatible stereo, and no initially foldable hydrogens.
//! Its comparison is exact representation equality.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_chem::element::Element;
use umol_graph::ops::transform::{
    FoldHydrogensError, HydrogenFolder, HydrogenUnfolder, Transformer, UnfoldHydrogensError,
};
use umol_graph_core::{GraphRemapping, NodeId, Remapping};
use umol_graph_ir::ir::{
    AromaticSystemForm, AsLit, AtomConstraintForm, AtomForm, AtomId, BondConstraintForm, BondForm,
    BondId, BooleanForm, CisTransStereoForm, Constraint, DativeBondForm, Entity, FluxionalityForm,
    IsotopeMassForm, LigandPermutation, Molecule, MoleculeEntries, MoleculeRemapping,
    MulticenterBondForm, NoncovalentBondForm, NumForm, StereoAtomConstraintForm, StereoAtomForm,
    StereoAtomId, StereoBondConstraintForm, StereoBondForm, StereoBondId, StereoKind, StereoLigand,
    StereoLigandKind, TetrahedralStereoForm,
};
use umol_perm::Permutation;

#[derive(Clone, Debug)]
struct HydrogenCase {
    source: Molecule,
    folded: Molecule,
    unfolded: Molecule,
    retained: Vec<(AtomId, AtomId)>,
}

fn hydrogen_cases(partial: bool) -> impl Strategy<Value = HydrogenCase> {
    let count = if partial {
        prop_oneof![
            (0i64..5).prop_map(NumForm::Lit),
            Just(NumForm::Undetermined),
            Just(NumForm::var("h")),
            Just(NumForm::from(vec![0, 1])),
        ]
        .boxed()
    } else {
        (0i64..5).prop_map(NumForm::Lit).boxed()
    };
    prop::collection::vec(
        (
            count,
            0usize..4,
            0u8..6,
            prop::sample::select(vec![Element::C, Element::N, Element::O, Element::S]),
            -2i64..3,
            any::<bool>(),
        ),
        0..9,
    )
    .prop_map(|fields| {
        let size = fields.len();
        let mut source = MoleculeEntries::default();
        let mut folded = MoleculeEntries::default();
        let mut unfolded = MoleculeEntries::default();
        let mut retained: Vec<usize> = (0..size).collect();
        for (i, (count, explicit, protection, element, charge, isotope)) in
            fields.iter().enumerate()
        {
            let backbone_degree = i64::from(i > 0) + i64::from(i + 1 < size);
            let protected = i64::from(*protection != 0);
            let degree = backbone_degree + *explicit as i64 + protected;
            let valence = degree + i64::from(*protection == 5);
            let mut attributes = AtomForm::from_element(*element).into_concrete();
            attributes.charge = NumForm::Lit(*charge);
            if *isotope {
                attributes.isotope_mass = IsotopeMassForm::Lit(13);
            }
            attributes.implicit_hydrogens = count.clone();
            attributes
                .constraints
                .set(AtomConstraintForm::Degree(NumForm::Lit(degree)));
            attributes
                .constraints
                .set(AtomConstraintForm::Valence(NumForm::Lit(valence)));
            attributes
                .constraints
                .set(AtomConstraintForm::TotalDegree(count + degree));
            attributes
                .constraints
                .set(AtomConstraintForm::TotalValence(count + valence));
            attributes
                .constraints
                .set(AtomConstraintForm::TotalHydrogens(
                    count + *explicit as i64 + protected,
                ));
            source.atoms.push(attributes.clone());
            let removed = if count.as_lit().is_some() {
                *explicit as i64
            } else {
                0
            };
            let added = count.as_lit().unwrap_or(0);
            let mut folded_atom = attributes.clone();
            if removed != 0 {
                folded_atom.implicit_hydrogens = count + removed;
                folded_atom
                    .constraints
                    .set(AtomConstraintForm::Degree(NumForm::Lit(degree - removed)));
                folded_atom
                    .constraints
                    .set(AtomConstraintForm::Valence(NumForm::Lit(valence - removed)));
            }
            folded.atoms.push(folded_atom);
            if added != 0 {
                attributes.implicit_hydrogens = NumForm::Lit(0);
                attributes
                    .constraints
                    .set(AtomConstraintForm::Degree(NumForm::Lit(degree + added)));
                attributes
                    .constraints
                    .set(AtomConstraintForm::Valence(NumForm::Lit(valence + added)));
            }
            unfolded.atoms.push(attributes);
        }
        for entries in [&mut source, &mut folded, &mut unfolded] {
            entries.bonds = (1..size)
                .map(|i| {
                    (
                        AtomId::from(i - 1),
                        AtomId::from(i),
                        BondForm::from_order(1).into_concrete(),
                    )
                })
                .collect();
            if size != 0 {
                entries.constraints = vec![Constraint::Not(Box::new(Constraint::Atom(
                    AtomId(0),
                    AtomConstraintForm::Degree(NumForm::Lit(0)),
                )))]
                .into();
            }
            if size >= 3 {
                entries.dative.push((
                    vec![AtomId(0)],
                    AtomId(1),
                    DativeBondForm::from_order(1).into_concrete(),
                ));
                entries.aromatic.push((
                    vec![AtomId(2), AtomId(0), AtomId(1)],
                    AromaticSystemForm::from_electrons(vec![1, 0, 2]).into_concrete(),
                ));
                entries.multicenter.push((
                    vec![AtomId(1), AtomId(2), AtomId(0)],
                    MulticenterBondForm::from_electrons(vec![2, 0, 1]).into_concrete(),
                ));
                entries
                    .noncovalent
                    .push(([AtomId(0), AtomId(2)], NoncovalentBondForm::default()));
            }
        }
        for (i, (count, explicit, protection, _, _, _)) in fields.iter().enumerate() {
            for _ in 0..*explicit {
                let hydrogen = AtomId::from(source.atoms.len());
                let attributes = AtomForm::from_element(Element::H).into_concrete();
                let bond = BondForm::from_order(1).into_concrete();
                source.atoms.push(attributes.clone());
                source.bonds.push((AtomId::from(i), hydrogen, bond.clone()));
                unfolded.atoms.push(attributes.clone());
                unfolded
                    .bonds
                    .push((AtomId::from(i), hydrogen, bond.clone()));
                if count.as_lit().is_none() {
                    retained.push(hydrogen.index());
                    let target = AtomId::from(folded.atoms.len());
                    folded.atoms.push(attributes);
                    folded.bonds.push((AtomId::from(i), target, bond));
                }
            }
            if *protection == 0 {
                continue;
            }
            let mut attributes = AtomForm::from_element(Element::H).into_concrete();
            match protection {
                1 => attributes.isotope_mass = IsotopeMassForm::Lit(2),
                2 => attributes.charge = NumForm::Lit(-1),
                3 => {
                    attributes
                        .constraints
                        .set(AtomConstraintForm::Degree(NumForm::Lit(1)));
                }
                _ => {}
            }
            let hydrogen = AtomId::from(source.atoms.len());
            let target = AtomId::from(folded.atoms.len());
            retained.push(hydrogen.index());
            for (entries, id) in [
                (&mut source, hydrogen),
                (&mut folded, target),
                (&mut unfolded, hydrogen),
            ] {
                entries.atoms.push(attributes.clone());
                entries.bonds.push((
                    AtomId::from(i),
                    id,
                    BondForm::from_order(if *protection == 5 { 2 } else { 1 }).into_concrete(),
                ));
                if *protection == 4 {
                    entries
                        .noncovalent
                        .push(([AtomId::from(i), id], NoncovalentBondForm::default()));
                }
            }
        }
        let mut new_bearers = Vec::new();
        for (i, (count, _, _, _, _, _)) in fields.iter().enumerate() {
            for _ in 0..count.as_lit().unwrap_or(0) {
                let hydrogen = AtomId::from(unfolded.atoms.len());
                unfolded
                    .atoms
                    .push(AtomForm::from_element(Element::H).into_concrete());
                unfolded.bonds.push((
                    AtomId::from(i),
                    hydrogen,
                    BondForm::from_order(1).into_concrete(),
                ));
                new_bearers.push(i);
            }
        }
        (
            Molecule::from_entries(source),
            Molecule::from_entries(folded),
            Molecule::from_entries(unfolded),
            retained,
            new_bearers,
        )
    })
    .prop_flat_map(|case| {
        let order = Just((0..case.0.atoms().count()).collect::<Vec<_>>()).prop_shuffle();
        (Just(case), order)
    })
    .prop_map(
        |((source, folded, unfolded, retained, new_bearers), images)| {
            let mut folded_order: Vec<_> = (0..retained.len()).collect();
            folded_order.sort_by_key(|&i| images[retained[i]]);
            let mut folded_images = vec![0; retained.len()];
            for (image, i) in folded_order.into_iter().enumerate() {
                folded_images[i] = image;
            }
            let retained = retained
                .iter()
                .enumerate()
                .map(|(i, &old)| (AtomId::from(images[old]), AtomId::from(folded_images[i])))
                .collect();
            let mut added_order: Vec<_> = (0..new_bearers.len()).collect();
            added_order.sort_by_key(|&i| images[new_bearers[i]]);
            let mut unfolded_images = images.clone();
            unfolded_images.resize(unfolded.atoms().count(), 0);
            for (image, i) in added_order.into_iter().enumerate() {
                unfolded_images[images.len() + i] = images.len() + image;
            }
            let original_bonds = source.bonds().count();
            let mut unfolded_bonds: Vec<_> = (0..unfolded.bonds().count()).collect();
            for i in 0..new_bearers.len() {
                unfolded_bonds[original_bonds + i] =
                    original_bonds + unfolded_images[images.len() + i] - images.len();
            }
            let folded_bonds = (0..folded.bonds().count()).collect();
            let [source, folded, unfolded] = [
                (source, images, (0..original_bonds).collect()),
                (folded, folded_images, folded_bonds),
                (unfolded, unfolded_images, unfolded_bonds),
            ]
            .map(|(molecule, atoms, bonds)| {
                let remapping = MoleculeRemapping::new(
                    GraphRemapping::new(
                        Remapping::new(atoms.into_iter().map(NodeId::from).collect()).unwrap(),
                        Remapping::new(bonds.into_iter().map(Into::into).collect()).unwrap(),
                    ),
                    Remapping::identity(molecule.dative_bonds().count()),
                    Remapping::identity(molecule.aromatic_systems().count()),
                    Remapping::identity(molecule.multicenter_bonds().count()),
                    Remapping::identity(molecule.noncovalent_bonds().count()),
                    Remapping::identity(0),
                    Remapping::identity(0),
                );
                molecule.remap(&remapping)
            });
            HydrogenCase {
                source,
                folded,
                unfolded,
                retained,
            }
        },
    )
}

fn hydrogen_count(molecule: &Molecule) -> Option<i64> {
    let implicit: Option<i64> = molecule
        .atoms()
        .iter()
        .map(|atom| atom.implicit_hydrogens().as_lit())
        .sum();
    Some(
        implicit?
            + molecule
                .atoms()
                .iter()
                .filter(|atom| atom.element().as_lit() == Some(Element::H))
                .count() as i64,
    )
}

fn stereo_cases() -> impl Strategy<Value = (Molecule, Molecule, Molecule)> {
    (
        Just((0usize..6).collect::<Vec<_>>()).prop_shuffle(),
        Just((0usize..4).collect::<Vec<_>>()).prop_shuffle(),
        any::<bool>(),
        any::<bool>(),
        any::<bool>(),
        0u32..2,
        0u32..2,
    )
        .prop_map(
            |(images, frame, first_swap, second_swap, endpoints_swap, atom_coset, bond_coset)| {
                let mut retained: Vec<_> = (0..4).collect();
                retained.sort_by_key(|&i| images[i]);
                let mut implicit_images = vec![0; 4];
                for (image, i) in retained.into_iter().enumerate() {
                    implicit_images[i] = image;
                }
                let mut unfolded_images = implicit_images.clone();
                unfolded_images.extend(if implicit_images[0] < implicit_images[2] {
                    [4, 5]
                } else {
                    [5, 4]
                });
                let [explicit, implicit, unfolded] = [images, implicit_images, unfolded_images]
                    .into_iter()
                    .enumerate()
                    .map(|(representation, images)| {
                        let implicit = representation == 1;
                        let mut atoms: Vec<_> = [
                            Element::N,
                            Element::F,
                            Element::C,
                            Element::Cl,
                            Element::H,
                            Element::H,
                        ]
                        .into_iter()
                        .take(images.len())
                        .enumerate()
                        .map(|(i, element)| {
                            let mut atom = AtomForm::from_element(element).into_concrete();
                            atom.implicit_hydrogens =
                                NumForm::Lit(i64::from(implicit && (i == 0 || i == 2)));
                            if i == 0 {
                                atom.lone_pairs = NumForm::Lit(1);
                            }
                            (AtomId::from(images[i]), atom)
                        })
                        .collect();
                        atoms.sort_by_key(|(id, _)| *id);
                        let hydrogen = |bearer, id| {
                            StereoLigand::new(
                                AtomId::from(images[if implicit { bearer } else { id }]),
                                if implicit {
                                    StereoLigandKind::ImplicitHydrogen
                                } else {
                                    StereoLigandKind::Atom
                                },
                            )
                        };
                        let ligands = [
                            StereoLigand::new(AtomId::from(images[2]), StereoLigandKind::Atom),
                            hydrogen(0, 4),
                            StereoLigand::new(AtomId::from(images[1]), StereoLigandKind::Atom),
                            StereoLigand::new(AtomId::from(images[0]), StereoLigandKind::LonePair),
                        ];
                        let atom_ligands: Vec<_> = frame.iter().map(|&i| ligands[i]).collect();
                        let mut bond_ligands = vec![
                            hydrogen(0, 4),
                            StereoLigand::new(AtomId::from(images[1]), StereoLigandKind::Atom),
                            StereoLigand::new(AtomId::from(images[3]), StereoLigandKind::Atom),
                            hydrogen(2, 5),
                        ];
                        if first_swap {
                            bond_ligands.swap(0, 1);
                        }
                        if second_swap {
                            bond_ligands.swap(2, 3);
                        }
                        if endpoints_swap {
                            bond_ligands.rotate_left(2);
                        }
                        let inversions = (0..4)
                            .flat_map(|i| (i + 1..4).map(move |j| (i, j)))
                            .filter(|&(i, j)| {
                                (atom_ligands[i].kind, atom_ligands[i].atom_id)
                                    > (atom_ligands[j].kind, atom_ligands[j].atom_id)
                            })
                            .count();
                        let side_inversions = bond_ligands
                            .chunks_exact(2)
                            .filter(|side| {
                                (side[0].kind, side[0].atom_id) > (side[1].kind, side[1].atom_id)
                            })
                            .count();
                        let tetrahedral = atom_coset ^ (inversions % 2) as u32;
                        let cis_trans = bond_coset ^ (side_inversions % 2) as u32;
                        atoms[images[0]]
                            .1
                            .constraints
                            .set(AtomConstraintForm::TetrahedralStereo(
                                TetrahedralStereoForm::stereo(tetrahedral),
                            ));
                        let mut bonds = vec![
                            (
                                AtomId::from(images[0]),
                                AtomId::from(images[2]),
                                BondForm::from_order(2)
                                    .with_constraint(BondConstraintForm::CisTransStereo(
                                        CisTransStereoForm::stereo(cis_trans),
                                    ))
                                    .into_concrete(),
                            ),
                            (
                                AtomId::from(images[0]),
                                AtomId::from(images[1]),
                                BondForm::from_order(1).into_concrete(),
                            ),
                            (
                                AtomId::from(images[2]),
                                AtomId::from(images[3]),
                                BondForm::from_order(1).into_concrete(),
                            ),
                        ];
                        if !implicit {
                            let mut hydrogen_bonds = [
                                (
                                    AtomId::from(images[0]),
                                    AtomId::from(images[4]),
                                    BondForm::from_order(1).into_concrete(),
                                ),
                                (
                                    AtomId::from(images[2]),
                                    AtomId::from(images[5]),
                                    BondForm::from_order(1).into_concrete(),
                                ),
                            ];
                            if representation == 2 {
                                hydrogen_bonds.sort_by_key(|(bearer, _, _)| *bearer);
                            }
                            bonds.extend(hydrogen_bonds);
                        }
                        let fluxionality = FluxionalityForm {
                            permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])),
                            active: BooleanForm::Lit(true),
                        };
                        Molecule::from_entries(MoleculeEntries {
                            atoms: atoms.into_iter().map(|(_, atom)| atom).collect(),
                            bonds,
                            stereo_atoms: vec![(
                                AtomId::from(images[0]),
                                atom_ligands,
                                StereoAtomForm::new(StereoKind::Tetrahedral, atom_coset)
                                    .with_constraint(StereoAtomConstraintForm::Fluxionality(
                                        fluxionality,
                                    )),
                            )],
                            stereo_bonds: vec![(
                                BondId(0),
                                bond_ligands,
                                StereoBondForm::new(StereoKind::CisTrans, bond_coset)
                                    .with_constraint(StereoBondConstraintForm::Fluxionality(
                                        fluxionality,
                                    )),
                            )],
                            ..Default::default()
                        })
                    })
                    .collect::<Vec<_>>()
                    .try_into()
                    .unwrap();
                (explicit, implicit, unfolded)
            },
        )
}

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            crate::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_hydrogen_folder_transform(
        literal in hydrogen_cases(false),
        partial in hydrogen_cases(true),
    ) {
        for case in [literal, partial] {
            let result = HydrogenFolder.transform(case.source.clone()).unwrap();
            prop_assert_eq!(&result, &case.folded);
            let mut in_place = case.source.clone();
            prop_assert_eq!(HydrogenFolder.transform_into(&mut in_place), Ok(()));
            prop_assert_eq!(&in_place, &result);
            prop_assert_eq!(HydrogenFolder.transform(result.clone()), Ok(result.clone()));
            prop_assert_eq!(result.clone().edit().finish(), Ok(result.clone()));
            if let Some(count) = hydrogen_count(&case.source) {
                prop_assert_eq!(hydrogen_count(&result), Some(count));
            }
            for (before, after) in case.retained {
                let before = case.source.atom(before);
                let after = result.atom(after);
                prop_assert_eq!(
                    (
                        after.total_valence(),
                        after.total_degree(),
                        after.total_hydrogens()
                    ),
                    (
                        before.total_valence(),
                        before.total_degree(),
                        before.total_hydrogens()
                    )
                );
            }
        }
    }

    #[test]
    fn test_hydrogen_folder_transform_stereo((source, expected, _) in stereo_cases()) {
        let result = HydrogenFolder.transform(source.clone()).unwrap();
        prop_assert_eq!(&result, &expected);
        let mut in_place = source;
        prop_assert_eq!(HydrogenFolder.transform_into(&mut in_place), Ok(()));
        prop_assert_eq!(&in_place, &result);
        prop_assert_eq!(HydrogenFolder.transform(result.clone()), Ok(result.clone()));
        prop_assert_eq!(result.clone().edit().finish(), Ok(result));
    }

    #[test]
    fn test_hydrogen_folder_transform_error(
        (source, _, _) in stereo_cases(),
        count in 1i64..5,
    ) {
        for conflict in [false, true] {
            let mut source = source.clone();
            let stereo = source.stereo_atom(StereoAtomId(0));
            let atom = stereo.site_id();
            let hydrogen = stereo
                .atom_ligand_ids()
                .find(|&id| source.atom(id).element().as_lit() == Some(Element::H))
                .unwrap();
            source.atom_mut(atom).attributes_mut().implicit_hydrogens =
                NumForm::Lit(if conflict { count } else { -count });
            let expected = if conflict {
                FoldHydrogensError::StereoHydrogenConflict {
                    entity: Entity::StereoAtom(StereoAtomId(0)),
                    hydrogen,
                }
            } else {
                FoldHydrogensError::NegativeImplicitHydrogenCount {
                    atom,
                    count: -count,
                }
            };
            prop_assert_eq!(
                HydrogenFolder.transform(source.clone()),
                Err(expected.clone())
            );
            let mut in_place = source.clone();
            prop_assert_eq!(HydrogenFolder.transform_into(&mut in_place), Err(expected));
            prop_assert_eq!(in_place, source);
        }
    }

    #[test]
    fn test_hydrogen_unfolder_transform(
        literal in hydrogen_cases(false),
        partial in hydrogen_cases(true),
    ) {
        for case in [literal, partial] {
            let result = HydrogenUnfolder.transform(case.source.clone()).unwrap();
            prop_assert_eq!(&result, &case.unfolded);
            let mut in_place = case.source.clone();
            prop_assert_eq!(HydrogenUnfolder.transform_into(&mut in_place), Ok(()));
            prop_assert_eq!(&in_place, &result);
            prop_assert_eq!(
                HydrogenUnfolder.transform(result.clone()),
                Ok(result.clone())
            );
            prop_assert_eq!(result.clone().edit().finish(), Ok(result.clone()));
            if let Some(count) = hydrogen_count(&case.source) {
                prop_assert_eq!(hydrogen_count(&result), Some(count));
            }
            for before in case.source.atoms().iter() {
                let after = result.atom(before.id());
                prop_assert_eq!(
                    (
                        after.total_valence(),
                        after.total_degree(),
                        after.total_hydrogens()
                    ),
                    (
                        before.total_valence(),
                        before.total_degree(),
                        before.total_hydrogens()
                    )
                );
            }
        }
    }

    #[test]
    fn test_hydrogen_unfolder_transform_stereo((_, source, expected) in stereo_cases()) {
        let result = HydrogenUnfolder.transform(source.clone()).unwrap();
        prop_assert_eq!(&result, &expected);
        let mut in_place = source;
        prop_assert_eq!(HydrogenUnfolder.transform_into(&mut in_place), Ok(()));
        prop_assert_eq!(&in_place, &result);
        prop_assert_eq!(HydrogenUnfolder.transform(result.clone()), Ok(result.clone()));
        prop_assert_eq!(result.clone().edit().finish(), Ok(result));
    }

    #[test]
    fn test_hydrogen_unfolder_transform_error(
        (_, source, _) in stereo_cases(),
        negative in -5i64..0,
        mismatch in prop_oneof![Just(0i64), 2i64..5],
        nonliteral in prop_oneof![
            Just(NumForm::Undetermined),
            Just(NumForm::from(vec![0, 1])),
            Just(NumForm::var("h")),
        ],
        bond in any::<bool>(),
    ) {
        let atom = if bond {
            source
                .atoms()
                .iter()
                .find(|atom| atom.element().as_lit() == Some(Element::C))
                .unwrap()
                .id()
        } else {
            source.stereo_atom(StereoAtomId(0)).site_id()
        };
        let entity = if bond {
            Entity::StereoBond(StereoBondId(0))
        } else {
            Entity::StereoAtom(StereoAtomId(0))
        };
        for value in [NumForm::Lit(negative), NumForm::Lit(mismatch), nonliteral] {
            let mut source = source.clone();
            let expected = match &value {
                NumForm::Lit(count) if *count < 0 => {
                    UnfoldHydrogensError::NegativeImplicitHydrogenCount {
                        atom,
                        count: *count,
                    }
                }
                NumForm::Lit(count) => UnfoldHydrogensError::StereoHydrogenCountMismatch {
                    entity,
                    atom,
                    implicit_hydrogens: *count,
                    virtual_hydrogens: 1,
                },
                _ => UnfoldHydrogensError::NonConcreteStereoHydrogenCount {
                    entity,
                    atom,
                    value: value.clone(),
                },
            };
            source.atom_mut(atom).attributes_mut().implicit_hydrogens = value;
            prop_assert_eq!(
                HydrogenUnfolder.transform(source.clone()),
                Err(expected.clone())
            );
            let mut in_place = source.clone();
            prop_assert_eq!(HydrogenUnfolder.transform_into(&mut in_place), Err(expected));
            prop_assert_eq!(in_place, source);
        }
    }

    #[test]
    fn test_hydrogen_unfolder_transform_reconstruction(
        source in prop_oneof![
            hydrogen_cases(false).prop_map(|case| case.folded),
            stereo_cases().prop_map(|(_, implicit, _)| implicit),
        ],
    ) {
        let explicit = HydrogenUnfolder.transform(source.clone()).unwrap();
        prop_assert_eq!(HydrogenFolder.transform(explicit), Ok(source));
    }
}
