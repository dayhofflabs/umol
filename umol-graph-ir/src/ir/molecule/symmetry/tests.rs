//! Unit tests for the symmetry, the colors and the remapping checks of a molecule.
//! Exhaustive enumeration uses dense, simple constraint-free graphs, literal binary
//! configurations and literal incidence weights, independently of production colors and search.

use std::cmp::Ordering;

use pretty_assertions::assert_eq;
use rstest::*;
use umol_chem::element::Element;
use umol_perm::Orientation;

use super::*;
use crate::ir::atom::AtomForm;
use crate::ir::bond::BondForm;
use crate::ir::id::{AtomId, BondId, StereoAtomId, StereoBondId, StereoLigandPosition};
use crate::ir::molecule::MoleculeEntries;
use crate::ir::stereo::{StereoAtomForm, StereoBondForm, StereoConfigurationForm, StereoKind};
use crate::ir::symmetry::{
    remapping_exists, remapping_preserves_incidences, AutomorphismAdapter, OrderedPartition,
};
use crate::ir::{
    ArithExpr, AromaticSystemConstraintForm, AromaticSystemForm, AromaticSystemId,
    AtomConstraintForm, BondConstraintForm, BooleanForm, Constraint, DativeBondConstraintForm,
    DativeBondForm, DativeBondId, ElementForm, GraphSymmetryView, LigandPermutation,
    LigandSymmetryForm, MoleculeConstraint, MulticenterBondConstraintForm, MulticenterBondForm,
    MulticenterBondId, NoncovalentBondConstraintForm, NoncovalentBondForm, NoncovalentBondId,
    NoncovalentBondKind, NumForm, OrientedLigandPermutation, StereoAtomConstraintForm,
    StereoBondConstraintForm, StereoLigandKind, StereoLigandPair, StereoTerm, StereogenicSiteKind,
    Stereogenicity, StereogenicityForm, Topicity, TopicityForm, TopicityRelationForm,
};
use crate::{mol_dsl, mol_dsl_concrete};

#[fixture]
fn benzene_ring() -> Molecule {
    let atoms = vec![AtomForm::from_element(Element::C); 6];
    let bonds = (0..6)
        .map(|i| (AtomId(i), AtomId((i + 1) % 6), BondForm::from_order(1)))
        .collect();
    Molecule::from_entries(MoleculeEntries {
        atoms,
        bonds,
        ..Default::default()
    })
}

// A tetrahedral center on atom 0 with the four given peripheral elements.
fn tetrahedral(peripherals: [Element; 4]) -> Molecule {
    let mut atoms = vec![AtomForm::from_element(Element::C)];
    atoms.extend(peripherals.into_iter().map(AtomForm::from_element));
    let bonds = (1..=4)
        .map(|i| (AtomId(0), AtomId(i), BondForm::from_order(1)))
        .collect();
    Molecule::from_entries(MoleculeEntries {
        atoms,
        bonds,
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        )],
        ..Default::default()
    })
}

#[fixture]
fn initial_color_molecule() -> Molecule {
    let normalized_three = NumForm::ArithExpr(Box::new(ArithExpr::Sum(vec![
        ArithExpr::Lit(1),
        ArithExpr::Lit(2),
    ])));
    let normalized_one = NumForm::ArithExpr(Box::new(ArithExpr::Sum(vec![
        ArithExpr::Lit(0),
        ArithExpr::Lit(1),
    ])));
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_charge(3_i64),
            AtomForm::from_element(Element::C)
                .with_charge(normalized_three)
                .with_constraint(AtomConstraintForm::Valence(NumForm::Lit(4))),
            AtomForm::from_element(Element::O).with_charge(3_i64),
            AtomForm::from_element(Element::C).with_charge(4_i64),
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::H),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::new(normalized_one.clone())),
            (AtomId(2), AtomId(3), BondForm::from_order(2)),
            (AtomId(0), AtomId(5), BondForm::from_order(1)),
        ],
        dative: vec![
            (vec![AtomId(0)], AtomId(4), DativeBondForm::from_order(1)),
            (
                vec![AtomId(1)],
                AtomId(4),
                DativeBondForm::new(normalized_one),
            ),
            (vec![AtomId(2)], AtomId(4), DativeBondForm::from_order(2)),
        ],
        aromatic: vec![
            (
                vec![AtomId(0), AtomId(1)],
                AromaticSystemForm::from_electrons(vec![1, 2]),
            ),
            (
                vec![AtomId(2), AtomId(3)],
                AromaticSystemForm::from_electrons(vec![2, 1]),
            ),
            (
                vec![AtomId(4), AtomId(5)],
                AromaticSystemForm::from_electrons(vec![1, 2]).with_charge(1_i64),
            ),
        ],
        multicenter: vec![
            (
                vec![AtomId(0), AtomId(2)],
                MulticenterBondForm::from_electrons(vec![1, 2]),
            ),
            (
                vec![AtomId(1), AtomId(3)],
                MulticenterBondForm::from_electrons(vec![2, 1]),
            ),
        ],
        noncovalent: vec![
            (
                [AtomId(0), AtomId(5)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            ),
            (
                [AtomId(1), AtomId(5)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            ),
            (
                [AtomId(2), AtomId(5)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::Ionic),
            ),
        ],
        stereo_atoms: vec![
            (
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            (
                AtomId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            ),
            (
                AtomId(2),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
                ],
                StereoAtomForm::new(StereoKind::SquarePlanar, StereoCoset::Lit(0)),
            ),
        ],
        stereo_bonds: vec![
            (
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
            ),
            (
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            ),
        ],
        ..Default::default()
    })
}

#[fixture]
fn constrained_color_molecule(initial_color_molecule: Molecule) -> Molecule {
    let mut molecule = initial_color_molecule;
    *molecule.atom_mut(AtomId(1)).attributes_mut() = AtomForm::from_element(Element::C)
        .with_charge(3)
        .with_constraint(AtomConstraintForm::valence(4));
    *molecule.bond_mut(BondId(1)).attributes_mut() =
        BondForm::from_order(1).with_constraint(BondConstraintForm::aromatic(true));
    *molecule.dative_bond_mut(DativeBondId(1)).attributes_mut() =
        DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::aromatic(true));
    *molecule
        .aromatic_system_mut(AromaticSystemId(1))
        .attributes_mut() = AromaticSystemForm::from_electrons(vec![2, 1])
        .with_constraint(AromaticSystemConstraintForm::electron_count(3));
    *molecule
        .multicenter_bond_mut(MulticenterBondId(1))
        .attributes_mut() = MulticenterBondForm::from_electrons(vec![2, 1])
        .with_constraint(MulticenterBondConstraintForm::electron_count(3));
    *molecule
        .noncovalent_bond_mut(NoncovalentBondId(1))
        .attributes_mut() = NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond)
        .with_constraint(NoncovalentBondConstraintForm::Intramolecular(
            BooleanForm::Lit(true),
        ));
    *molecule.stereo_atom_mut(StereoAtomId(1)).attributes_mut() =
        StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32).with_constraint(
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
                Stereogenicity::Stereogenic,
            )),
        );
    *molecule.stereo_bond_mut(StereoBondId(1)).attributes_mut() =
        StereoBondForm::new(StereoKind::CisTrans, 1_u32).with_constraint(
            StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
                Stereogenicity::Stereogenic,
            )),
        );
    molecule
}

fn exhaustive_atom_remappings(molecule: &Molecule) -> Vec<Vec<u32>> {
    let count = molecule.atoms().count();
    let mut bonds = vec![vec![None; count]; count];
    for bond in molecule.bonds().iter() {
        let [first, second] = bond.atom_ids();
        bonds[first.index()][second.index()] = Some(bond.attributes());
        bonds[second.index()][first.index()] = Some(bond.attributes());
    }

    fn enumerate(
        molecule: &Molecule,
        bonds: &[Vec<Option<&BondForm>>],
        images: &mut Vec<u32>,
        used: &mut [bool],
        output: &mut Vec<Vec<u32>>,
    ) {
        let source = images.len();
        if source == used.len() {
            if preserves_overlay_incidences(molecule, images) {
                output.push(images.clone());
            }
            return;
        }
        for target in 0..used.len() {
            if used[target]
                || molecule.atom(AtomId(source as u32)).attributes()
                    != molecule.atom(AtomId(target as u32)).attributes()
                || images
                    .iter()
                    .enumerate()
                    .any(|(other, &image)| bonds[source][other] != bonds[target][image as usize])
            {
                continue;
            }
            images.push(target as u32);
            used[target] = true;
            enumerate(molecule, bonds, images, used, output);
            used[target] = false;
            images.pop();
        }
    }

    let mut output = Vec::new();
    enumerate(
        molecule,
        &bonds,
        &mut Vec::new(),
        &mut vec![false; count],
        &mut output,
    );
    output
}

fn weighted_participants(
    atoms: impl Iterator<Item = AtomId>,
    electrons: &ElectronCountsForm,
    images: &[u32],
) -> Vec<(u32, i64)> {
    let ElectronCountsForm::Lit(counts) = electrons else {
        panic!("exhaustive enumeration requires literal electron counts");
    };
    let atoms = atoms.collect::<Vec<_>>();
    assert_eq!(atoms.len(), counts.len());
    let mut result = atoms
        .into_iter()
        .zip(counts)
        .map(|(atom, &count)| (images[atom.index()], count))
        .collect::<Vec<_>>();
    result.sort();
    result
}

fn preserves_overlay_incidences(molecule: &Molecule, images: &[u32]) -> bool {
    let identity = (0..molecule.atoms().count() as u32).collect::<Vec<_>>();
    let dative = |mapping: &[u32]| {
        let mut result = molecule
            .dative_bonds()
            .iter()
            .map(|bond| {
                let mut donors = bond
                    .donor_ids()
                    .map(|atom| mapping[atom.index()])
                    .collect::<Vec<_>>();
                donors.sort();
                (
                    mapping[bond.acceptor_id().index()],
                    donors,
                    bond.attributes().clone(),
                )
            })
            .collect::<Vec<_>>();
        result.sort();
        result
    };
    let aromatic = |mapping: &[u32]| {
        let mut result = molecule
            .aromatic_systems()
            .iter()
            .map(|system| {
                let form = system.attributes();
                (
                    weighted_participants(system.atom_ids(), &form.electrons, mapping),
                    form.charge.clone(),
                    form.unpaired_electrons.clone(),
                    form.constraints.clone(),
                )
            })
            .collect::<Vec<_>>();
        result.sort();
        result
    };
    let multicenter = |mapping: &[u32]| {
        let mut result = molecule
            .multicenter_bonds()
            .iter()
            .map(|bond| {
                let form = bond.attributes();
                (
                    weighted_participants(bond.atom_ids(), &form.electrons, mapping),
                    form.charge.clone(),
                    form.unpaired_electrons.clone(),
                    form.constraints.clone(),
                )
            })
            .collect::<Vec<_>>();
        result.sort();
        result
    };
    let noncovalent = |mapping: &[u32]| {
        let mut result = molecule
            .noncovalent_bonds()
            .iter()
            .map(|bond| {
                let mut atoms = bond.atom_ids().map(|atom| mapping[atom.index()]);
                atoms.sort();
                (atoms, bond.attributes().clone())
            })
            .collect::<Vec<_>>();
        result.sort();
        result
    };
    dative(images) == dative(&identity)
        && aromatic(images) == aromatic(&identity)
        && multicenter(images) == multicenter(&identity)
        && noncovalent(images) == noncovalent(&identity)
}

struct StereoRecord {
    site: Vec<AtomId>,
    ligands: Vec<StereoLigand>,
    kind: StereoKind,
    coset: u32,
}

fn stereo_records(molecule: &Molecule) -> Vec<StereoRecord> {
    let configuration = |form: &StereoConfigurationForm| {
        let StereoConfigurationForm::Kinded(kind, StereoCoset::Lit(coset)) = form else {
            panic!("exhaustive symmetry checks require literal configurations");
        };
        (*kind, *coset)
    };
    molecule
        .stereo_atoms()
        .iter()
        .map(|stereo| {
            let (kind, coset) = configuration(&stereo.attributes().configuration);
            StereoRecord {
                site: vec![stereo.site_id()],
                ligands: stereo.ligand_ids().to_vec(),
                kind,
                coset,
            }
        })
        .chain(molecule.stereo_bonds().iter().map(|stereo| {
            let (kind, coset) = configuration(&stereo.attributes().configuration);
            let mut site = molecule.bond(stereo.site_id()).atom_ids().to_vec();
            site.sort();
            StereoRecord {
                site,
                ligands: stereo.ligand_ids().to_vec(),
                kind,
                coset,
            }
        }))
        .collect()
}

fn preserves_stereo_configuration(
    records: &[StereoRecord],
    images: &[u32],
    mirrored: bool,
    changed: Option<(usize, u32)>,
) -> bool {
    for (source_index, source) in records.iter().enumerate() {
        let mut site = source
            .site
            .iter()
            .map(|atom| AtomId(images[atom.index()]))
            .collect::<Vec<_>>();
        site.sort();
        let Some(target_index) = records.iter().position(|target| target.site == site) else {
            return false;
        };
        if changed.is_some_and(|(index, _)| index == source_index && target_index != index) {
            return false;
        }
        let target = &records[target_index];
        if source.kind != target.kind {
            return false;
        }
        let permutation = source
            .ligands
            .iter()
            .map(|ligand| {
                let image = StereoLigand::new(AtomId(images[ligand.atom_id.index()]), ligand.kind);
                target
                    .ligands
                    .iter()
                    .position(|candidate| *candidate == image)
            })
            .collect::<Option<Vec<_>>>();
        let Some(permutation) = permutation else {
            return false;
        };
        let transported = match source.kind {
            StereoKind::Tetrahedral => {
                let inversions = (0..4)
                    .flat_map(|first| (first + 1..4).map(move |second| (first, second)))
                    .filter(|&(first, second)| permutation[first] > permutation[second])
                    .count();
                source.coset ^ (inversions % 2) as u32
            }
            StereoKind::CisTrans => {
                if permutation[0] / 2 != permutation[1] / 2
                    || permutation[2] / 2 != permutation[3] / 2
                {
                    return false;
                }
                source.coset ^ ((permutation[0] % 2) ^ (permutation[2] % 2)) as u32
            }
            _ => panic!("exhaustive symmetry checks require binary stereo kinds"),
        };
        let mut expected = changed
            .filter(|(index, _)| *index == target_index)
            .map_or(target.coset, |(_, coset)| coset);
        if mirrored && source.kind == StereoKind::Tetrahedral {
            expected ^= 1;
        }
        if transported != expected {
            return false;
        }
    }
    true
}

fn exhaustive_symmetry(molecule: &Molecule) -> (bool, Vec<u32>, Vec<u32>) {
    let mappings = exhaustive_atom_remappings(molecule);
    let records = stereo_records(molecule);
    let proper = mappings
        .iter()
        .filter(|images| preserves_stereo_configuration(&records, images, false, None))
        .collect::<Vec<_>>();
    let improper = mappings
        .iter()
        .filter(|images| preserves_stereo_configuration(&records, images, true, None))
        .collect::<Vec<_>>();
    let orbits = |star: bool| {
        (0..molecule.atoms().count())
            .map(|atom| {
                proper
                    .iter()
                    .copied()
                    .chain(improper.iter().copied().filter(|_| star))
                    .map(|images| images[atom])
                    .min()
                    .unwrap()
            })
            .collect()
    };
    (improper.is_empty(), orbits(false), orbits(true))
}

fn exhaustive_stereogenicities(molecule: &Molecule) -> Vec<Stereogenicity> {
    let mappings = exhaustive_atom_remappings(molecule);
    let records = stereo_records(molecule);
    records
        .iter()
        .enumerate()
        .map(|(index, record)| {
            let stereogenic = |fixed: Option<AtomId>| {
                (0..2).filter(|&coset| coset != record.coset).any(|coset| {
                    !mappings
                        .iter()
                        .filter(|images| fixed.is_none_or(|atom| images[atom.index()] == atom.0))
                        .any(|images| {
                            preserves_stereo_configuration(
                                &records,
                                images,
                                false,
                                Some((index, coset)),
                            )
                        })
                })
            };
            if stereogenic(None) {
                Stereogenicity::Stereogenic
            } else if record.ligands.iter().any(|ligand| {
                ligand.kind == StereoLigandKind::Atom && stereogenic(Some(ligand.atom_id))
            }) {
                Stereogenicity::Prostereogenic
            } else {
                Stereogenicity::Nonstereogenic
            }
        })
        .collect()
}

fn exhaustive_topicity(
    molecule: &Molecule,
    record: usize,
    first: usize,
    second: usize,
) -> Topicity {
    let mappings = exhaustive_atom_remappings(molecule);
    let records = stereo_records(molecule);
    let focus = &records[record];
    let first = focus.ligands[first].atom_id;
    let second = focus.ligands[second].atom_id;
    let candidates = mappings
        .iter()
        .filter(|images| {
            let mut site = focus
                .site
                .iter()
                .map(|atom| AtomId(images[atom.index()]))
                .collect::<Vec<_>>();
            site.sort();
            site == focus.site && images[first.index()] == second.0
        })
        .collect::<Vec<_>>();
    if candidates
        .iter()
        .any(|images| preserves_stereo_configuration(&records, images, false, None))
    {
        Topicity::Homotopic
    } else if candidates
        .iter()
        .any(|images| preserves_stereo_configuration(&records, images, true, None))
    {
        Topicity::Enantiotopic
    } else if !candidates.is_empty() {
        Topicity::Diastereotopic
    } else {
        Topicity::Distinct
    }
}

#[fixture]
fn dichlorodifluoroethane_chiral() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C" "C" "H" "F" "Cl" "H" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"] [1 6 "1"] [1 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                       {:site 1 :ligands [0 5 6 7] :attrs "Th0"}]
    }"#
    )
}

#[fixture]
fn dichlorodifluoroethane_meso() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C" "C" "H" "F" "Cl" "H" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"] [1 6 "1"] [1 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                       {:site 1 :ligands [0 5 6 7] :attrs "Th1"}]
    }"#
    )
}

#[fixture]
fn dichlorodifluoroethane_chiral_renumbered() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C" "C" "H" "F" "Cl" "F" "H" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"] [1 6 "1"] [1 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                       {:site 1 :ligands [0 5 6 7] :attrs "Th1"}]
    }"#
    )
}

#[fixture]
fn dichlorodifluoroethane_meso_renumbered() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C" "C" "H" "F" "Cl" "F" "H" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"] [1 6 "1"] [1 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                       {:site 1 :ligands [0 5 6 7] :attrs "Th0"}]
    }"#
    )
}

#[fixture]
fn trichloropentane_meso() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h3" "C" "C" "C" "C#h3" "Cl" "Cl" "Cl" "H" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [1 5 "1"] [2 6 "1"] [3 7 "1"]
                [1 8 "1"] [2 9 "1"] [3 10 "1"]]
        :stereo-atoms [{:site 1 :ligands [0 2 5 8] :attrs "Th0"}
                       {:site 2 :ligands [1 3 6 9] :attrs "Th0"}
                       {:site 3 :ligands [2 4 7 10] :attrs "Th0"}]
    }"#
    )
}

#[fixture]
fn trichloropentane_meso_inverted_c3() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h3" "C" "C" "C" "C#h3" "Cl" "Cl" "Cl" "H" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [1 5 "1"] [2 6 "1"] [3 7 "1"]
                [1 8 "1"] [2 9 "1"] [3 10 "1"]]
        :stereo-atoms [{:site 1 :ligands [0 2 5 8] :attrs "Th0"}
                       {:site 2 :ligands [1 3 6 9] :attrs "Th1"}
                       {:site 3 :ligands [2 4 7 10] :attrs "Th0"}]
    }"#
    )
}

#[fixture]
fn trichloropentane_chiral() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h3" "C" "C" "C" "C#h3" "Cl" "Cl" "Cl" "H" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [1 5 "1"] [2 6 "1"] [3 7 "1"]
                [1 8 "1"] [2 9 "1"] [3 10 "1"]]
        :stereo-atoms [{:site 1 :ligands [0 2 5 8] :attrs "Th0"}
                       {:site 2 :ligands [1 3 6 9] :attrs "Th0"}
                       {:site 3 :ligands [2 4 7 10] :attrs "Th1"}]
    }"#
    )
}

#[fixture]
fn dichlorocyclobutane_cis() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C" "C#h2" "C" "C#h2" "Cl" "Cl" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 0 "1"] [0 4 "1"] [2 5 "1"] [0 6 "1"] [2 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 3 4 6] :attrs "Th0"}
                       {:site 2 :ligands [1 3 5 7] :attrs "Th1"}]
    }"#
    )
}

#[fixture]
fn dichlorocyclobutane_trans() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C" "C#h2" "C" "C#h2" "Cl" "Cl" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 0 "1"] [0 4 "1"] [2 5 "1"] [0 6 "1"] [2 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 3 4 6] :attrs "Th0"}
                       {:site 2 :ligands [1 3 5 7] :attrs "Th0"}]
    }"#
    )
}

#[fixture]
fn chlorofluoropropane() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C" "C" "C#h3" "F" "Cl" "H" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [0 3 "1"] [0 4 "1"] [0 5 "1"] [1 6 "1"] [1 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 3 4 5] :attrs "Th0"}
                       {:site 1 :ligands [0 2 6 7] :attrs "Th0"}]
    }"#
    )
}

#[fixture]
fn chlorofluoropropane_without_c1() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C" "C" "C#h3" "F" "Cl" "H" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [0 3 "1"] [0 4 "1"] [0 5 "1"] [1 6 "1"] [1 7 "1"]]
        :stereo-atoms [{:site 1 :ligands [0 2 6 7] :attrs "Th0"}]
    }"#
    )
}

#[fixture]
fn trimethylcyclohexane_all_cis() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C" "C#h2" "C" "C#h2" "C" "C#h2" "C#h3" "C#h3" "C#h3" "H" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                [0 6 "1"] [2 7 "1"] [4 8 "1"] [0 9 "1"] [2 10 "1"] [4 11 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 5 6 9] :attrs "Th0"}
                       {:site 2 :ligands [1 3 7 10] :attrs "Th1"}
                       {:site 4 :ligands [3 5 8 11] :attrs "Th1"}]
    }"#
    )
}

#[fixture]
fn trimethylcyclohexane_cis_trans() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C" "C#h2" "C" "C#h2" "C" "C#h2" "C#h3" "C#h3" "C#h3" "H" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                [0 6 "1"] [2 7 "1"] [4 8 "1"] [0 9 "1"] [2 10 "1"] [4 11 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 5 6 9] :attrs "Th0"}
                       {:site 2 :ligands [1 3 7 10] :attrs "Th0"}
                       {:site 4 :ligands [3 5 8 11] :attrs "Th0"}]
    }"#
    )
}

#[rstest]
fn test_molecule_graph_symmetry_constitutional(benzene_ring: Molecule) {
    let symmetry = benzene_ring.graph_symmetry(&GraphSymmetryConfig::default());
    assert!(symmetry.same_proper_orbit(AtomId(0), AtomId(3)));
    assert_eq!(
        symmetry.proper_orbit_of(AtomId(0)),
        vec![
            AtomId(0),
            AtomId(1),
            AtomId(2),
            AtomId(3),
            AtomId(4),
            AtomId(5),
        ],
    );
    assert!(!symmetry.is_chiral());
}

#[rstest]
fn test_molecule_graph_symmetry_chiral() {
    // Four distinct ligands ⇒ a genuine stereocenter, no compensating symmetry.
    let mol = tetrahedral([Element::F, Element::Cl, Element::Br, Element::I]);
    let symmetry = mol.graph_symmetry(&GraphSymmetryConfig::default());
    assert!(symmetry.is_chiral());
}

#[rstest]
fn test_molecule_graph_symmetry_enantiotopic() {
    // Two identical Cl ligands: enantiotopic — same star orbit, distinct proper.
    let mol = tetrahedral([Element::Cl, Element::Cl, Element::F, Element::Br]);
    let symmetry = mol.graph_symmetry(&GraphSymmetryConfig::default());
    assert!(!symmetry.is_chiral());
    assert!(!symmetry.same_proper_orbit(AtomId(1), AtomId(2)));
    assert!(symmetry.same_star_orbit(AtomId(1), AtomId(2)));
}

#[rstest]
#[case::atom(
    MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 5],
        bonds: (1..=4)
            .map(|id| (AtomId(0), AtomId(id), BondForm::from_order(1)))
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            (1..=4)
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                .collect(),
            StereoAtomForm::default(),
        )],
        ..Default::default()
    },
    AtomId(1),
    vec![AtomId(1), AtomId(2), AtomId(3), AtomId(4)],
)]
#[case::bond(
    MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 6],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(2)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(1), AtomId(4), BondForm::from_order(1)),
            (AtomId(1), AtomId(5), BondForm::from_order(1)),
        ],
        stereo_bonds: vec![(
            BondId(0),
            (2..=5)
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                .collect(),
            StereoBondForm::default(),
        )],
        ..Default::default()
    },
    AtomId(2),
    vec![AtomId(2), AtomId(3), AtomId(4), AtomId(5)],
)]
fn test_molecule_graph_symmetry_kindless(
    #[case] entries: MoleculeEntries,
    #[case] atom: AtomId,
    #[case] expected_orbit: Vec<AtomId>,
) {
    let molecule = Molecule::try_from_entries(entries).expect("kindless stereo is admitted");
    for include_constraints in [false, true] {
        let symmetry = molecule.graph_symmetry(&GraphSymmetryConfig {
            include_constraints,
            ..Default::default()
        });

        assert_eq!(symmetry.proper_orbit_of(atom), expected_orbit);
        assert_eq!(symmetry.is_chiral(), false);
    }
}

#[rstest]
fn test_molecule_graph_symmetry_out_of_range_coset() {
    let mut mol = tetrahedral([Element::F, Element::Cl, Element::Br, Element::I]);
    mol.stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 9);
    let symmetry = mol.graph_symmetry(&GraphSymmetryConfig::default());
    assert_eq!(symmetry.proper_orbit_of(AtomId(0)), vec![AtomId(0)]);
}

#[rstest]
#[case::dichlorodifluoroethane_chiral(dichlorodifluoroethane_chiral(), (true, vec![0, 0, 2, 3, 4, 2, 3, 4], vec![0, 0, 2, 3, 4, 2, 3, 4]))]
#[case::dichlorodifluoroethane_meso(dichlorodifluoroethane_meso(), (false, vec![0, 1, 2, 3, 4, 5, 6, 7], vec![0, 0, 2, 3, 4, 2, 3, 4]))]
#[case::dichlorodifluoroethane_chiral_renumbered(dichlorodifluoroethane_chiral_renumbered(), (true, vec![0, 0, 2, 3, 4, 3, 2, 4], vec![0, 0, 2, 3, 4, 3, 2, 4]))]
#[case::dichlorodifluoroethane_meso_renumbered(dichlorodifluoroethane_meso_renumbered(), (false, vec![0, 1, 2, 3, 4, 5, 6, 7], vec![0, 0, 2, 3, 4, 3, 2, 4]))]
#[case::trichloropentane_meso(trichloropentane_meso(), (false, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], vec![0, 1, 2, 1, 0, 5, 6, 5, 8, 9, 8]))]
#[case::trichloropentane_meso_inverted_c3(trichloropentane_meso_inverted_c3(), (false, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], vec![0, 1, 2, 1, 0, 5, 6, 5, 8, 9, 8]))]
#[case::trichloropentane_chiral(trichloropentane_chiral(), (true, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]))]
#[case::dichlorocyclobutane_cis(dichlorocyclobutane_cis(), (false, vec![0, 1, 0, 1, 4, 4, 6, 6], vec![0, 1, 0, 1, 4, 4, 6, 6]))]
#[case::dichlorocyclobutane_trans(dichlorocyclobutane_trans(), (false, vec![0, 1, 0, 3, 4, 4, 6, 6], vec![0, 1, 0, 1, 4, 4, 6, 6]))]
#[case::chlorofluoropropane(chlorofluoropropane(), (true, vec![0, 1, 2, 3, 4, 5, 6, 7], vec![0, 1, 2, 3, 4, 5, 6, 7]))]
#[case::chlorofluoropropane_without_c1(chlorofluoropropane_without_c1(), (false, vec![0, 1, 2, 3, 4, 5, 6, 7], vec![0, 1, 2, 3, 4, 5, 6, 6]))]
#[case::butene_cis(
    mol_dsl_concrete!(r#"{:atoms ["C#h3" "C#h" "C#h" "C#h3"]
      :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
      :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct0"}]}"#),
    (false, vec![0, 1, 1, 0], vec![0, 1, 1, 0])
)]
#[case::butene_trans(
    mol_dsl_concrete!(r#"{:atoms ["C#h3" "C#h" "C#h" "C#h3"]
      :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
      :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#),
    (false, vec![0, 1, 1, 0], vec![0, 1, 1, 0])
)]
fn test_molecule_graph_symmetry(
    #[case] molecule: Molecule,
    #[case] expected: (bool, Vec<u32>, Vec<u32>),
    #[values(false, true)] include_constraints: bool,
) {
    let original = molecule.clone();
    assert_eq!(molecule.check_integrity(), Ok(()));
    assert_eq!(exhaustive_symmetry(&molecule), expected);
    let config = GraphSymmetryConfig {
        include_constraints,
        ..Default::default()
    };
    let symmetry = molecule.graph_symmetry(&config).into_graph_symmetry();
    let owned = molecule.graph_symmetry(&config);
    let borrowed = GraphSymmetryView::from_graph_symmetry(&symmetry, &molecule);
    assert_eq!(molecule.is_chiral(&symmetry), expected.0);
    for view in [owned, borrowed] {
        assert_eq!(view.is_chiral(), expected.0);
        for atom in molecule.atoms().ids() {
            let proper = molecule
                .atoms()
                .ids()
                .filter(|other| expected.1[atom.index()] == expected.1[other.index()])
                .collect::<Vec<_>>();
            let star = molecule
                .atoms()
                .ids()
                .filter(|other| expected.2[atom.index()] == expected.2[other.index()])
                .collect::<Vec<_>>();
            assert_eq!(molecule.proper_orbit_of(&symmetry, atom), proper);
            assert_eq!(molecule.star_orbit_of(&symmetry, atom), star);
            assert_eq!(view.proper_orbit_of(atom), proper);
            assert_eq!(view.star_orbit_of(atom), star);
            for other in molecule.atoms().ids() {
                assert_eq!(
                    molecule.same_proper_orbit(&symmetry, atom, other),
                    proper.contains(&other)
                );
                assert_eq!(
                    molecule.same_star_orbit(&symmetry, atom, other),
                    star.contains(&other)
                );
                assert_eq!(view.same_proper_orbit(atom, other), proper.contains(&other));
                assert_eq!(view.same_star_orbit(atom, other), star.contains(&other));
            }
        }
    }
    assert_eq!(molecule, original);
}

#[rstest]
#[case::term(StereoCoset::term(StereoTerm::Lit(0)))]
#[case::singleton(StereoCoset::lit_set([0]))]
fn test_molecule_graph_symmetry_raw(
    #[case] coset: StereoCoset,
    #[values(false, true)] include_constraints: bool,
) {
    let mut molecule = dichlorodifluoroethane_meso();
    molecule
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = StereoConfigurationForm::kinded(StereoKind::Tetrahedral, coset);
    let original = molecule.clone();
    let config = GraphSymmetryConfig {
        include_constraints,
        ..Default::default()
    };
    let view = molecule.graph_symmetry(&config);
    assert_eq!(view.is_chiral(), true);
    for atom in molecule.atoms().ids() {
        assert_eq!(view.proper_orbit_of(atom), vec![atom]);
        assert_eq!(view.star_orbit_of(atom), vec![atom]);
    }
    let normalized = molecule.clone().normalize().unwrap();
    assert_eq!(normalized.graph_symmetry(&config).is_chiral(), false);
    assert_eq!(molecule, original);
}

#[rstest]
fn test_molecule_graph_symmetry_raw_element() {
    let mut molecule = dichlorodifluoroethane_meso();
    molecule.atom_mut(AtomId(0)).attributes_mut().element = ElementForm::lit_set([Element::C]);
    let original = molecule.clone();
    let view = molecule.graph_symmetry(&GraphSymmetryConfig::default());
    assert_eq!(view.is_chiral(), true);
    assert_eq!(view.star_orbit_of(AtomId(0)), vec![AtomId(0)]);
    assert_eq!(
        molecule
            .clone()
            .normalize()
            .unwrap()
            .graph_symmetry(&GraphSymmetryConfig::default())
            .is_chiral(),
        false
    );
    assert_eq!(molecule, original);
}

#[rstest]
#[case::contradictory(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::lit_set([])), [false, false])]
#[case::incompatible_frame(StereoAtomConstraintForm::Topicity(TopicityForm {
    pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(7)),
    relation: TopicityRelationForm::Lit(Topicity::Homotopic),
}), [false, true])]
fn test_molecule_graph_symmetry_constraints(
    #[case] constraint: StereoAtomConstraintForm,
    #[case] expected_chiral: [bool; 2],
) {
    let mut molecule = dichlorodifluoroethane_meso();
    for id in [StereoAtomId(0), StereoAtomId(1)] {
        molecule
            .stereo_atom_mut(id)
            .attributes_mut()
            .constraints
            .set(constraint.clone());
    }
    let original = molecule.clone();
    for include_constraints in [false, true] {
        let view = molecule.graph_symmetry(&GraphSymmetryConfig {
            include_constraints,
            ..Default::default()
        });
        assert_eq!(
            view.is_chiral(),
            expected_chiral[usize::from(include_constraints)]
        );
        assert_eq!(
            view.same_star_orbit(AtomId(0), AtomId(1)),
            !expected_chiral[usize::from(include_constraints)]
        );
    }
    assert_eq!(molecule, original);
}

#[rstest]
#[case::excluded(false, false)]
#[case::included(true, true)]
fn test_molecule_graph_symmetry_constraint_frame(
    #[case] include_constraints: bool,
    #[case] expected_chiral: bool,
) {
    let mut molecule = mol_dsl_concrete!(
        r#"{
        :atoms ["C" "F" "F" "Cl" "Br"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]
    }"#
    );
    molecule
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .constraints
        .set(StereoAtomConstraintForm::Topicity(TopicityForm {
            pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2)),
            relation: TopicityRelationForm::Lit(Topicity::Homotopic),
        }));
    assert_eq!(
        molecule
            .graph_symmetry(&GraphSymmetryConfig {
                include_constraints,
                ..Default::default()
            })
            .is_chiral(),
        expected_chiral
    );
}

#[rstest]
fn test_molecule_graph_symmetry_molecule_constraints(
    #[values(false, true)] include_constraints: bool,
) {
    let mut molecule = dichlorodifluoroethane_meso();
    molecule
        .constraints_mut()
        .push(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(0)]),
            sum: NumForm::lit_set([]),
        }));
    let original = molecule.clone();
    let view = molecule.graph_symmetry(&GraphSymmetryConfig {
        include_constraints,
        ..Default::default()
    });
    assert_eq!(view.is_chiral(), false);
    assert_eq!(view.proper_orbit_of(AtomId(0)), vec![AtomId(0)]);
    assert_eq!(view.star_orbit_of(AtomId(0)), vec![AtomId(0), AtomId(1)]);
    assert_eq!(molecule, original);
}

#[rstest]
#[case::kind_frame_mismatch(mol_dsl_concrete!(r#"{
    :atoms ["C" "F" "Cl" "Br" "I"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Tb0"}]
}"#))]
fn test_molecule_graph_symmetry_encoding(
    #[case] molecule: Molecule,
    #[values(false, true)] include_constraints: bool,
) {
    let original = molecule.clone();
    let symmetry = molecule.graph_symmetry(&GraphSymmetryConfig {
        include_constraints,
        ..Default::default()
    });
    assert_eq!(symmetry.is_chiral(), true);
    assert_eq!(symmetry.proper_orbit_of(AtomId(0)), vec![AtomId(0)]);
    assert_eq!(molecule, original);
}

#[rstest]
#[case::square_planar(
    mol_dsl_concrete!(r#"{
        :atoms ["Pt" "Cl" "Cl" "Cl" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Sp0"}]
    }"#),
    vec![0, 1, 1, 1, 1],
)]
#[case::different_rings(
    mol_dsl_concrete!(r#"{
        :atoms ["C#h2" "C#h2" "C#h2" "C#h2" "C#h2" "C#h2" "C#h2" "C#h2" "C#h2"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]
                [3 4 "1"] [4 5 "1"] [5 6 "1"] [6 7 "1"] [7 8 "1"] [8 3 "1"]]
    }"#),
    vec![0, 0, 0, 1, 1, 1, 1, 1, 1],
)]
fn test_molecule_same_proper_orbit(#[case] molecule: Molecule, #[case] orbit_classes: Vec<u32>) {
    let symmetry = molecule
        .graph_symmetry(&GraphSymmetryConfig::default())
        .into_graph_symmetry();
    for (a, &a_class) in orbit_classes.iter().enumerate() {
        for (b, &b_class) in orbit_classes.iter().enumerate() {
            assert_eq!(
                molecule.same_proper_orbit(&symmetry, AtomId(a as u32), AtomId(b as u32)),
                a_class == b_class,
            );
        }
    }
}

#[rstest]
#[case::empty(r#"{:atoms []}"#, vec![])]
#[case::no_records(r#"{:atoms ["C#h4"]}"#, vec![])]
#[case::kindless(r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "*"}]}"#, vec![])]
#[case::constitutional(r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]}"#, vec![(Entity::StereoAtom(StereoAtomId(0)), StereogenicSiteKind::Constitutional)])]
#[case::prostereogenic(r#"{:atoms ["C" "H" "H" "F" "Cl"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#, vec![(Entity::StereoAtom(StereoAtomId(0)), StereogenicSiteKind::Prostereogenic)])]
#[case::nonstereogenic(r#"{:atoms ["C" "H" "H" "H" "H"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]}"#, vec![(Entity::StereoAtom(StereoAtomId(0)), StereogenicSiteKind::Nonstereogenic)])]
#[case::square_planar(r#"{:atoms ["Pt#c-2" "N#c+#h3" "N#c+#h3" "Cl" "Cl"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Sp*"}]}"#, vec![(Entity::StereoAtom(StereoAtomId(0)), StereogenicSiteKind::Configurational)])]
#[case::trigonal_bipyramidal(r#"{:atoms ["P" "Cl" "F" "F" "F" "F"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [0 5 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4 5] :attrs "Tb*"}]}"#, vec![(Entity::StereoAtom(StereoAtomId(0)), StereogenicSiteKind::Configurational)])]
#[case::single_ring_center(r#"{:atoms ["C#h" "C#h2" "C#h" "C#h2" "Cl" "Cl"] :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 0 "1"] [0 4 "1"] [2 5 "1"]] :stereo-atoms [{:site 0 :ligands [1 3 4 [:h 0]] :attrs "Th*"}]}"#, vec![(Entity::StereoAtom(StereoAtomId(0)), StereogenicSiteKind::Prostereogenic)])]
#[case::mixed_order_and_gaps(r#"{:atoms ["C" "F" "Cl" "Br" "I" "C" "F" "Cl" "Br" "I" "C" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
:bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "1"] [5 7 "1"] [5 8 "1"] [5 9 "1"] [10 11 "1"] [10 12 "1"] [10 13 "1"] [10 14 "1"] [15 16 "2"] [15 17 "1"] [16 18 "1"]]
:stereo-atoms [{:site 10 :ligands [11 12 13 14] :attrs "Th*"} {:site 5 :ligands [6 7 8 9] :attrs "*"} {:site 0 :ligands [1 2 3 4] :attrs "Th1"}]
:stereo-bonds [{:site 12 :ligands [17 [:h 15] 18 [:h 16]] :attrs "Ct*"}]}"#, vec![(Entity::StereoAtom(StereoAtomId(0)), StereogenicSiteKind::Constitutional), (Entity::StereoAtom(StereoAtomId(2)), StereogenicSiteKind::Constitutional), (Entity::StereoBond(StereoBondId(0)), StereogenicSiteKind::Constitutional)])]
#[case::kindless_bond(r#"{:atoms ["C#h" "C#h" "F" "Cl"] :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
    :stereo-bonds [{:site 0 :ligands [2 [:h 0] 3 [:h 1]] :attrs "*"}]}"#, vec![])]
#[case::bond_order_and_gaps(r#"{
    :atoms ["C#h" "C#h" "F" "Cl" "C#h" "C#h" "F" "Cl" "C#h" "C#h" "F" "Cl"]
    :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"] [4 5 "2"] [4 6 "1"] [5 7 "1"] [8 9 "2"] [8 10 "1"] [9 11 "1"]]
    :stereo-bonds [{:site 6 :ligands [10 [:h 8] 11 [:h 9]] :attrs "Ct*"}
                   {:site 3 :ligands [6 [:h 4] 7 [:h 5]] :attrs "*"}
                   {:site 0 :ligands [2 [:h 0] 3 [:h 1]] :attrs "Ct*"}]}"#,
    vec![(Entity::StereoBond(StereoBondId(0)), StereogenicSiteKind::Constitutional),
         (Entity::StereoBond(StereoBondId(2)), StereogenicSiteKind::Constitutional)])]
#[case::different_ligand_cycles(r#"{
    :atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C" "C" "C" "C" "C" "C" "C" "F" "Cl"]
    :bonds [[0 1 "1"] [0 8 "1"] [0 15 "1"] [0 16 "1"]
            [1 2 "1"] [1 3 "1"] [1 4 "1"] [1 5 "1"] [1 6 "1"] [1 7 "1"]
            [2 3 "1"] [3 4 "1"] [4 2 "1"] [5 6 "1"] [6 7 "1"] [7 5 "1"]
            [8 9 "1"] [8 10 "1"] [8 11 "1"] [8 12 "1"] [8 13 "1"] [8 14 "1"]
            [9 10 "1"] [10 11 "1"] [11 12 "1"] [12 13 "1"] [13 14 "1"] [14 9 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 8 15 16] :attrs "Th*"}]
}"#, vec![(Entity::StereoAtom(StereoAtomId(0)), StereogenicSiteKind::Constitutional)])]
fn test_molecule_stereogenic_site_kinds(
    #[case] dsl: &str,
    #[case] expected: Vec<(Entity, StereogenicSiteKind)>,
    #[values(false, true)] include_constraints: bool,
) {
    let molecule = mol_dsl_concrete!(dsl);
    let original = molecule.clone();
    let symmetry = molecule
        .graph_symmetry(&GraphSymmetryConfig {
            include_constraints,
            ..Default::default()
        })
        .into_graph_symmetry();
    assert_eq!(molecule.stereogenic_site_kinds(&symmetry), expected);
    assert_eq!(molecule, original);
}

#[rstest]
#[case::invalid_coset(StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 9))]
#[case::degree_mismatch(StereoConfigurationForm::kinded(StereoKind::TrigonalBipyramidal, 0))]
fn test_molecule_stereogenic_site_kinds_encoding(
    #[case] configuration: StereoConfigurationForm,
    #[values(false, true)] include_constraints: bool,
) {
    let mut molecule = mol_dsl_concrete!(
        r#"{:atoms ["C" "H" "H" "H" "H"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#
    );
    molecule
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = configuration;
    let original = molecule.clone();
    let symmetry = molecule
        .graph_symmetry(&GraphSymmetryConfig {
            include_constraints,
            ..Default::default()
        })
        .into_graph_symmetry();
    assert_eq!(
        molecule.stereogenic_site_kinds(&symmetry),
        vec![(
            Entity::StereoAtom(StereoAtomId(0)),
            StereogenicSiteKind::Constitutional
        )]
    );
    assert_eq!(molecule, original);
}

#[rstest]
#[case::excluded(false, StereogenicSiteKind::Prostereogenic)]
#[case::included(true, StereogenicSiteKind::Constitutional)]
fn test_molecule_stereogenic_site_kinds_constraints(
    #[case] include_constraints: bool,
    #[case] expected: StereogenicSiteKind,
) {
    let mut molecule = mol_dsl_concrete!(
        r#"{:atoms ["C" "H" "H" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#
    );
    molecule
        .atom_mut(AtomId(1))
        .attributes_mut()
        .constraints
        .set(AtomConstraintForm::valence(1));
    let original = molecule.clone();
    let symmetry = molecule.graph_symmetry(&GraphSymmetryConfig {
        include_constraints,
        ..Default::default()
    });
    assert_eq!(
        symmetry.stereogenic_site_kinds(),
        vec![(Entity::StereoAtom(StereoAtomId(0)), expected)]
    );
    assert_eq!(molecule, original);
}

#[rstest]
fn test_remapping_exists_at_constitution() {
    let molecule = chlorofluoropropane_without_c1();
    let symmetry = molecule
        .graph_symmetry(&GraphSymmetryConfig::default())
        .into_graph_symmetry();
    assert_eq!(
        remapping_exists_at(
            &molecule,
            &symmetry,
            DescriptionLevel::Constitution,
            Orientation::Proper,
            &[AtomId(1), AtomId(6)],
            &[AtomId(1), AtomId(7)]
        ),
        true
    );
}

#[rstest]
#[case::free(vec![], true)]
#[case::one_ligand(vec![AtomId(2)], false)]
fn test_remapping_exists_for(#[case] fixed: Vec<AtomId>, #[case] expected: bool) {
    let molecule = mol_dsl_concrete!(
        r#"{
        :atoms ["C#h" "F" "Cl" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]
    }"#
    );
    let symmetry = molecule
        .graph_symmetry(&GraphSymmetryConfig::default())
        .into_graph_symmetry();
    assert_eq!(
        remapping_exists_for(
            &molecule,
            &symmetry,
            Entity::StereoAtom(StereoAtomId(0)),
            1,
            &fixed
        ),
        expected
    );
}

#[rstest]
#[case::tetrahedral(mol_dsl_concrete!(r#"{
            :atoms ["C#h" "F" "Cl" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]
        }"#), vec![true, true])]
#[case::two_centers(mol_dsl_concrete!(r#"{
            :atoms ["C" "C" "Cl" "Cl" "H" "H" "H" "H"]
            :bonds [[0 1 "1"] [0 2 "1"] [1 3 "1"] [0 4 "1"] [0 5 "1"] [1 6 "1"] [1 7 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 4 5] :attrs "Th0"}
                           {:site 1 :ligands [0 3 6 7] :attrs "Th0"}]
        }"#), vec![true, true])]
#[case::square_planar(mol_dsl_concrete!(r#"{
            :atoms ["Pt#c-2" "N#c+#h3" "N#c+#h3" "Cl" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Sp0"}]
        }"#), vec![true, false, true])]
#[case::trigonal_bipyramidal(mol_dsl_concrete!(r#"{
            :atoms ["P" "Cl" "F" "F" "F" "F"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [0 5 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4 5] :attrs "Tb0"}]
        }"#), vec![true, true, true, true, true, true, true, true,
                   false, false, false, false, false, false, false, false,
                   false, false, false, false])]
#[case::remote_ammonia(mol_dsl_concrete!(r#"{
            :atoms ["Pt#c-2" "N#c+" "N#c+" "Cl" "Cl" "H" "H" "H" "H" "H" "H"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]
                    [1 5 "1"] [1 6 "1"] [1 7 "1"] [2 8 "1"] [2 9 "1"] [2 10 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Sp0"}]
        }"#), vec![true, false, true])]
fn test_remapping_exists_for_cosets(
    #[case] molecule: Molecule,
    #[case] expected: Vec<bool>,
    #[values(false, true)] include_constraints: bool,
) {
    let original = molecule.clone();
    let symmetry = molecule
        .graph_symmetry(&GraphSymmetryConfig {
            include_constraints,
            ..Default::default()
        })
        .into_graph_symmetry();
    let actual = (0..expected.len())
        .map(|coset| {
            remapping_exists_for(
                &molecule,
                &symmetry,
                Entity::StereoAtom(StereoAtomId(0)),
                coset as u32,
                &[],
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(actual, expected);
    assert_eq!(molecule, original);
}

#[rstest]
#[case::equivalent_atom(
    mol_dsl_concrete!(
        r#"{:atoms ["C" "C" "H" "F" "Cl" "H" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"] [1 6 "1"] [1 7 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                           {:site 1 :ligands [0 5 6 7] :attrs "Th0"}]}"#
    ),
    vec![0],
    vec![1],
    true
)]
#[case::mirror_image_atom(
    mol_dsl_concrete!(
        r#"{:atoms ["C" "C" "H" "F" "Cl" "H" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"] [1 6 "1"] [1 7 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                           {:site 1 :ligands [0 5 6 7] :attrs "Th1"}]}"#
    ),
    vec![0],
    vec![1],
    false
)]
#[case::atom_fixed_equivalent_atom(
    mol_dsl_concrete!(
        r#"{:atoms ["C" "C#h2" "C" "C#h2" "Cl" "Cl" "H" "H"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 0 "1"]
                    [0 4 "1"] [2 5 "1"] [0 6 "1"] [2 7 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 3 4 6] :attrs "Th0"}
                           {:site 2 :ligands [1 3 5 7] :attrs "Th0"}]}"#
    ),
    vec![1, 0],
    vec![1, 2],
    true
)]
#[case::atom_fixed_mirror_image_atom(
    mol_dsl_concrete!(
        r#"{:atoms ["C" "C#h2" "C" "C#h2" "Cl" "Cl" "H" "H"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 0 "1"]
                    [0 4 "1"] [2 5 "1"] [0 6 "1"] [2 7 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 3 4 6] :attrs "Th0"}
                           {:site 2 :ligands [1 3 5 7] :attrs "Th1"}]}"#
    ),
    vec![1, 0],
    vec![1, 2],
    false
)]
#[case::different_colors(
    mol_dsl_concrete!(
        r#"{:atoms ["C" "C" "H" "F" "Cl" "H" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"] [1 6 "1"] [1 7 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                           {:site 1 :ligands [0 5 6 7] :attrs "Th0"}]}"#
    ),
    vec![0],
    vec![3],
    false
)]
fn test_remapping_exists_individualized(
    #[case] molecule: Molecule,
    #[case] source: Vec<u32>,
    #[case] target: Vec<u32>,
    #[case] expected: bool,
) {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let colors = initial_colors(
        &molecule,
        &incidence_graph,
        DescriptionLevel::Structure,
        true,
    )
    .unwrap();
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let partition = |atoms: Vec<u32>| {
        atoms.into_iter().fold(
            OrderedPartition::from_descriptors(&adapter.colors),
            |partition, atom| {
                let node = incidence_graph.node_of(Entity::Atom(AtomId(atom)));
                let cell_index = partition
                    .cells
                    .iter()
                    .position(|cell| cell.contains(&node))
                    .unwrap();
                partition.individualize(cell_index, node)
            },
        )
    };

    assert_eq!(
        remapping_exists(
            &adapter,
            partition(source),
            &adapter,
            partition(target),
            AutomorphismAlgorithm::Nauty,
            &|images: &[NodeId]| {
                remapping_preserves_incidences(&adapter, &adapter, images)
                    && remapping_preserves_stereo(&molecule, &molecule, &incidence_graph, images)
            },
            true,
            &|generators: &mut Vec<Vec<NodeId>>| {
                retain_stereo_preserving_generators(&molecule, &incidence_graph, generators);
            },
        ),
        expected,
    );
}

#[rstest]
#[case::moved_record(mol_dsl_concrete!(r#"{
    :atoms ["C" "C" "Cl" "Cl" "H" "H" "H" "H"]
    :bonds [[0 1 "1"] [0 2 "1"] [1 3 "1"] [0 4 "1"] [0 5 "1"] [1 6 "1"] [1 7 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 4 5] :attrs "Th0"}
                   {:site 1 :ligands [0 3 6 7] :attrs "Th0"}]
}"#), vec![1, 0, 3, 2, 6, 7, 4, 5, 8, 10, 9, 13, 14, 11, 12],
      vec![0, 1, 2, 3, 5, 4, 6, 7, 8, 9, 10, 12, 11, 13, 14])]
#[case::fixed_record_square_planar(mol_dsl_concrete!(r#"{
    :atoms ["Pt#c-2" "F" "F" "F" "F"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Sp0"}]
}"#), vec![0, 3, 2, 1, 4, 7, 6, 5, 8], vec![0, 1, 3, 2, 4, 5, 7, 6, 8])]
fn test_remapping_matches_stereo_generators(
    #[case] molecule: Molecule,
    #[case] generator: Vec<u32>,
    #[case] remapping: Vec<u32>,
    #[values(false, true)] include_constraints: bool,
) {
    let incidence = molecule.incidence_graph(IncidenceLevel::Constitution);
    let generator = generator.into_iter().map(NodeId).collect::<Vec<_>>();
    let remapping = remapping.into_iter().map(NodeId).collect::<Vec<_>>();
    let composed = remapping
        .iter()
        .map(|node| generator[node.index()])
        .collect::<Vec<_>>();
    let focus = Entity::StereoAtom(StereoAtomId(0));
    assert_eq!(
        (
            remapping_matches_stereo(
                &molecule,
                &incidence,
                &generator,
                include_constraints,
                Orientation::Proper,
                None
            ),
            remapping_changes_coset(
                &molecule,
                &incidence,
                &remapping,
                focus,
                1,
                include_constraints
            ),
            remapping_changes_coset(
                &molecule,
                &incidence,
                &composed,
                focus,
                1,
                include_constraints
            ),
        ),
        (true, true, false),
    );
}

#[rstest]
#[case::meso_carbons_exchanged(
    mol_dsl_concrete!(
        r#"{:atoms ["C" "C" "H" "F" "Cl" "H" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"] [1 6 "1"] [1 7 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                           {:site 1 :ligands [0 5 6 7] :attrs "Th1"}]}"#
    ),
    vec![1, 0, 5, 6, 7, 2, 3, 4, 8, 12, 13, 14, 9, 10, 11, 16, 15],
    true
)]
#[case::chiral_carbons_exchanged(
    mol_dsl_concrete!(
        r#"{:atoms ["C" "C" "H" "F" "Cl" "H" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"] [1 6 "1"] [1 7 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                           {:site 1 :ligands [0 5 6 7] :attrs "Th0"}]}"#
    ),
    vec![1, 0, 5, 6, 7, 2, 3, 4, 8, 12, 13, 14, 9, 10, 11, 16, 15],
    false
)]
#[case::implicit_hydrogen_equal_ligands_exchanged(
    mol_dsl_concrete!(
        r#"{:atoms ["C#h" "F" "Cl" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]}"#
    ),
    vec![0, 1, 3, 2, 4, 6, 5, 7],
    true
)]
#[case::implicit_hydrogen_identity(
    mol_dsl_concrete!(
        r#"{:atoms ["C#h" "F" "Cl" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]}"#
    ),
    vec![0, 1, 2, 3, 4, 5, 6, 7],
    false
)]
#[case::stereo_bond_identity(
    mol_dsl_concrete!(
        r#"{:atoms ["C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
            :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#
    ),
    vec![0, 1, 2, 3, 4, 5, 6, 7],
    true
)]
#[case::undetermined_coset_identity(
    mol_dsl_concrete!(
        r#"{:atoms ["C#h" "F" "Cl" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th*"}]}"#
    ),
    vec![0, 1, 2, 3, 4, 5, 6, 7],
    true
)]
#[case::out_of_range_coset_identity(
    mol_dsl_concrete!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th9"}]}"#
    ),
    vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    false
)]
fn test_remapping_mirrors_stereo(
    #[case] molecule: Molecule,
    #[case] images: Vec<u32>,
    #[case] expected: bool,
    #[values(false, true)] include_constraints: bool,
) {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let images = images.into_iter().map(NodeId).collect::<Vec<_>>();

    assert_eq!(
        remapping_mirrors_stereo(&molecule, &incidence_graph, &images, include_constraints),
        expected,
    );
}

#[rstest]
#[case::literal(StereoCoset::Lit(0))]
#[case::singleton(StereoCoset::lit_set([0]))]
#[case::undetermined(StereoCoset::Undetermined)]
#[case::term_literal(StereoCoset::term(StereoTerm::Lit(0)))]
#[case::identity_apply(StereoCoset::term(StereoTerm::apply(
    StereoTerm::var("configuration"),
    Permutation::identity(4)
)))]
#[case::nested_mirrors(StereoCoset::term(StereoTerm::mirror(StereoTerm::mirror(
    StereoTerm::var("configuration")
))))]
#[case::empty_domain(StereoCoset::term(StereoTerm::var_in("configuration", [])))]
fn test_remapping_mirrors_stereo_frame(
    #[case] coset: StereoCoset,
    #[values(false, true)] include_constraints: bool,
) {
    let mut molecule = mol_dsl_concrete!(
        r#"{
        :atoms ["C" "F" "F" "Cl" "Br"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]
    }"#
    );
    molecule
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = StereoConfigurationForm::kinded(StereoKind::Tetrahedral, coset);
    let original = molecule.clone();
    let incidence = molecule.incidence_graph(IncidenceLevel::Full);
    assert_eq!(
        remapping_mirrors_stereo(
            &molecule,
            &incidence,
            &[
                NodeId(0),
                NodeId(2),
                NodeId(1),
                NodeId(3),
                NodeId(4),
                NodeId(6),
                NodeId(5),
                NodeId(7),
                NodeId(8),
                NodeId(9)
            ],
            include_constraints
        ),
        true
    );
    assert_eq!(
        molecule
            .graph_symmetry(&GraphSymmetryConfig {
                include_constraints,
                ..Default::default()
            })
            .is_chiral(),
        false
    );
    assert_eq!(molecule, original);
}

#[rstest]
#[case::implicit_hydrogen_equal_ligands_exchanged(
    mol_dsl_concrete!(
        r#"{:atoms ["C#h" "F" "Cl" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]}"#
    ),
    vec![0, 1, 3, 2, 4, 6, 5, 7],
    Entity::StereoAtom(StereoAtomId(0)),
    1,
    true
)]
#[case::implicit_hydrogen_identity(
    mol_dsl_concrete!(
        r#"{:atoms ["C#h" "F" "Cl" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]}"#
    ),
    vec![0, 1, 2, 3, 4, 5, 6, 7],
    Entity::StereoAtom(StereoAtomId(0)),
    1,
    false
)]
#[case::entity_changed_other_preserved(
    mol_dsl_concrete!(
        r#"{:atoms ["C" "C" "C#h3" "F" "Cl" "H" "H" "H"]
            :bonds [[0 1 "1"] [1 2 "1"] [0 3 "1"] [0 4 "1"] [0 5 "1"] [1 6 "1"] [1 7 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 3 4 5] :attrs "Th0"}
                           {:site 1 :ligands [0 2 6 7] :attrs "Th0"}]}"#
    ),
    vec![0, 1, 2, 3, 4, 5, 7, 6, 8, 9, 10, 11, 12, 14, 13, 15, 16],
    Entity::StereoAtom(StereoAtomId(1)),
    1,
    true
)]
#[case::other_entity_changed(
    mol_dsl_concrete!(
        r#"{:atoms ["C" "C" "C#h3" "F" "Cl" "H" "H" "H"]
            :bonds [[0 1 "1"] [1 2 "1"] [0 3 "1"] [0 4 "1"] [0 5 "1"] [1 6 "1"] [1 7 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 3 4 5] :attrs "Th0"}
                           {:site 1 :ligands [0 2 6 7] :attrs "Th0"}]}"#
    ),
    vec![0, 1, 2, 3, 4, 5, 7, 6, 8, 9, 10, 11, 12, 14, 13, 15, 16],
    Entity::StereoAtom(StereoAtomId(0)),
    1,
    false
)]
#[case::stereo_bond_equal_ligands_exchanged(
    mol_dsl_concrete!(
        r#"{:atoms ["C" "C" "F" "F" "Cl" "Br"]
            :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
            :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct0"}]}"#
    ),
    vec![0, 1, 3, 2, 4, 5, 6, 8, 7, 9, 10, 11],
    Entity::StereoBond(StereoBondId(0)),
    1,
    true
)]
#[case::undetermined_coset_equal_ligands_exchanged(
    mol_dsl_concrete!(
        r#"{:atoms ["C#h" "F" "Cl" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th*"}]}"#
    ),
    vec![0, 1, 3, 2, 4, 6, 5, 7],
    Entity::StereoAtom(StereoAtomId(0)),
    1,
    true
)]
#[case::undetermined_coset_identity(
    mol_dsl_concrete!(
        r#"{:atoms ["C#h" "F" "Cl" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th*"}]}"#
    ),
    vec![0, 1, 2, 3, 4, 5, 6, 7],
    Entity::StereoAtom(StereoAtomId(0)),
    1,
    false
)]
#[case::out_of_range_coset_identity(
    mol_dsl_concrete!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th9"}]}"#
    ),
    vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    Entity::StereoAtom(StereoAtomId(0)),
    0,
    false
)]
#[case::kindless_identity(
    mol_dsl_concrete!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "*"}]}"#
    ),
    vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    Entity::StereoAtom(StereoAtomId(0)),
    0,
    false
)]
fn test_remapping_changes_coset(
    #[case] molecule: Molecule,
    #[case] images: Vec<u32>,
    #[case] entity: Entity,
    #[case] coset: u32,
    #[case] expected: bool,
    #[values(false, true)] include_constraints: bool,
) {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let images = images.into_iter().map(NodeId).collect::<Vec<_>>();

    assert_eq!(
        remapping_changes_coset(
            &molecule,
            &incidence_graph,
            &images,
            entity,
            coset,
            include_constraints
        ),
        expected,
    );
}

#[rstest]
#[case::rotation(
    mol_dsl!(
        r#"{:atoms ["C" "F" "F" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["C" "F" "F" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
    ),
    vec![0, 2, 3, 1, 4, 6, 7, 5, 8, 9],
    true
)]
#[case::transposition(
    mol_dsl!(
        r#"{:atoms ["C" "F" "F" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["C" "F" "F" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
    ),
    vec![0, 2, 1, 3, 4, 6, 5, 7, 8, 9],
    false
)]
#[case::site_moved(
    mol_dsl!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
    ),
    vec![1, 0, 2, 3, 4, 5, 6, 7, 8, 9],
    false
)]
#[case::renumbered(
    mol_dsl!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["I" "Br" "Cl" "F" "C"]
            :bonds [[4 0 "1"] [4 1 "1"] [4 2 "1"] [4 3 "1"]]
            :stereo-atoms [{:site 4 :ligands [3 2 1 0] :attrs "Th1"}]}"#
    ),
    vec![4, 3, 2, 1, 0, 8, 7, 6, 5, 9],
    true
)]
#[case::renumbered_inverted(
    mol_dsl!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["I" "Br" "Cl" "F" "C"]
            :bonds [[4 0 "1"] [4 1 "1"] [4 2 "1"] [4 3 "1"]]
            :stereo-atoms [{:site 4 :ligands [3 2 1 0] :attrs "Th0"}]}"#
    ),
    vec![4, 3, 2, 1, 0, 8, 7, 6, 5, 9],
    false
)]
#[case::stereo_bond(
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
            :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
            :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#
    ),
    vec![0, 1, 2, 3, 4, 5, 6, 7],
    true
)]
#[case::stereo_bond_inverted(
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
            :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
            :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct0"}]}"#
    ),
    vec![0, 1, 2, 3, 4, 5, 6, 7],
    false
)]
fn test_remapping_preserves_stereo(
    #[case] source: Molecule,
    #[case] target: Molecule,
    #[case] images: Vec<u32>,
    #[case] expected: bool,
) {
    let incidence_graph = source.incidence_graph(IncidenceLevel::Full);
    let images = images.into_iter().map(NodeId).collect::<Vec<_>>();

    assert_eq!(
        remapping_preserves_stereo(&source, &target, &incidence_graph, &images),
        expected,
    );
}

#[rstest]
#[case::feature_free(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        ..Default::default()
    }),
    vec![NodeId(1), NodeId(0), NodeId(2)],
    true,
)]
#[case::distinct_ligands(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
            AtomForm::from_element(Element::I),
        ],
        bonds: (1..=4)
            .map(|ligand| (AtomId(0), AtomId(ligand), BondForm::from_order(1)))
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            (1..=4)
                .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                .collect(),
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        )],
        ..Default::default()
    }),
    (0..10).map(NodeId::from).collect(),
    true,
)]
#[case::prostereogenic_explicit_hydrogens(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::H),
            AtomForm::from_element(Element::H),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
        ],
        bonds: (1..=4)
            .map(|ligand| (AtomId(0), AtomId(ligand), BondForm::from_order(1)))
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            (1..=4)
                .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                .collect(),
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        )],
        ..Default::default()
    }),
    vec![
        NodeId(0), NodeId(2), NodeId(1), NodeId(3), NodeId(4),
        NodeId(6), NodeId(5), NodeId(7), NodeId(8), NodeId(9),
    ],
    false,
)]
#[case::stereo_site_exchange(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
            AtomForm::from_element(Element::I),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
            AtomForm::from_element(Element::I),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(0), AtomId(4), BondForm::from_order(1)),
            (AtomId(5), AtomId(6), BondForm::from_order(1)),
            (AtomId(5), AtomId(7), BondForm::from_order(1)),
            (AtomId(5), AtomId(8), BondForm::from_order(1)),
            (AtomId(5), AtomId(9), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![
            (
                AtomId(0),
                (1..=4)
                    .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            (
                AtomId(5),
                (6..=9)
                    .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
        ],
        ..Default::default()
    }),
    vec![
        NodeId(5), NodeId(6), NodeId(7), NodeId(8), NodeId(9),
        NodeId(0), NodeId(1), NodeId(2), NodeId(3), NodeId(4),
        NodeId(14), NodeId(15), NodeId(16), NodeId(17),
        NodeId(10), NodeId(11), NodeId(12), NodeId(13),
        NodeId(19), NodeId(18),
    ],
    true,
)]
#[case::stereo_bond_endpoint_exchange(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(2)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(1), AtomId(4), BondForm::from_order(1)),
            (AtomId(1), AtomId(5), BondForm::from_order(1)),
        ],
        stereo_bonds: vec![(
            BondId(0),
            [2, 3, 4, 5]
                .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                .into(),
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
        )],
        ..Default::default()
    }),
    vec![
        NodeId(1), NodeId(0), NodeId(4), NodeId(5), NodeId(2), NodeId(3),
        NodeId(6), NodeId(9), NodeId(10), NodeId(7), NodeId(8), NodeId(11),
    ],
    true,
)]
#[case::stereo_site_exchange_one_topicity(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
            AtomForm::from_element(Element::I),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
            AtomForm::from_element(Element::I),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(0), AtomId(4), BondForm::from_order(1)),
            (AtomId(5), AtomId(6), BondForm::from_order(1)),
            (AtomId(5), AtomId(7), BondForm::from_order(1)),
            (AtomId(5), AtomId(8), BondForm::from_order(1)),
            (AtomId(5), AtomId(9), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![
            (
                AtomId(0),
                (1..=4)
                    .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0))
                    .with_constraint(StereoAtomConstraintForm::Topicity(TopicityForm {
                        pair: StereoLigandPair::new(0usize.into(), 1usize.into()),
                        relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
                    })),
            ),
            (
                AtomId(5),
                (6..=9)
                    .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
        ],
        ..Default::default()
    }),
    vec![
        NodeId(5), NodeId(6), NodeId(7), NodeId(8), NodeId(9),
        NodeId(0), NodeId(1), NodeId(2), NodeId(3), NodeId(4),
        NodeId(14), NodeId(15), NodeId(16), NodeId(17),
        NodeId(10), NodeId(11), NodeId(12), NodeId(13),
        NodeId(19), NodeId(18),
    ],
    false,
)]
#[case::stereo_site_exchange_both_topicity(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
            AtomForm::from_element(Element::I),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
            AtomForm::from_element(Element::I),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(0), AtomId(4), BondForm::from_order(1)),
            (AtomId(5), AtomId(6), BondForm::from_order(1)),
            (AtomId(5), AtomId(7), BondForm::from_order(1)),
            (AtomId(5), AtomId(8), BondForm::from_order(1)),
            (AtomId(5), AtomId(9), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![
            (
                AtomId(0),
                (1..=4)
                    .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0))
                    .with_constraint(StereoAtomConstraintForm::Topicity(TopicityForm {
                        pair: StereoLigandPair::new(0usize.into(), 1usize.into()),
                        relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
                    })),
            ),
            (
                AtomId(5),
                (6..=9)
                    .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0))
                    .with_constraint(StereoAtomConstraintForm::Topicity(TopicityForm {
                        pair: StereoLigandPair::new(0usize.into(), 1usize.into()),
                        relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
                    })),
            ),
        ],
        ..Default::default()
    }),
    vec![
        NodeId(5), NodeId(6), NodeId(7), NodeId(8), NodeId(9),
        NodeId(0), NodeId(1), NodeId(2), NodeId(3), NodeId(4),
        NodeId(14), NodeId(15), NodeId(16), NodeId(17),
        NodeId(10), NodeId(11), NodeId(12), NodeId(13),
        NodeId(19), NodeId(18),
    ],
    true,
)]
fn test_retain_stereo_preserving_generators(
    #[case] molecule: Molecule,
    #[case] generator: Vec<NodeId>,
    #[case] expected_retained: bool,
) {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let mut generators = vec![generator.clone()];
    let expected = if expected_retained {
        vec![generator]
    } else {
        Vec::new()
    };

    retain_stereo_preserving_generators(&molecule, &incidence_graph, &mut generators);

    assert_eq!(generators, expected);
}

#[rstest]
#[case::bond_endpoint(
    mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
    Entity::Bond(BondId(0)),
    AtomId(0),
    (Incidence::BondEndpoint, None)
)]
#[case::dative_donor(
    mol_dsl!(r#"{:atoms ["N" "B"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#),
    Entity::DativeBond(DativeBondId(0)),
    AtomId(0),
    (Incidence::DativeDonor, None)
)]
#[case::dative_acceptor(
    mol_dsl!(r#"{:atoms ["N" "B"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#),
    Entity::DativeBond(DativeBondId(0)),
    AtomId(1),
    (Incidence::DativeAcceptor, None)
)]
#[case::aromatic_atom_oxygen(
    mol_dsl!(
        r#"{:atoms ["O" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [0 4 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#
    ),
    Entity::AromaticSystem(AromaticSystemId(0)),
    AtomId(0),
    (Incidence::AromaticAtom, Some(2))
)]
#[case::aromatic_atom_carbon(
    mol_dsl!(
        r#"{:atoms ["O" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [0 4 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#
    ),
    Entity::AromaticSystem(AromaticSystemId(0)),
    AtomId(1),
    (Incidence::AromaticAtom, Some(1))
)]
#[case::aromatic_atom_without_counts(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 6],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4), AtomId(5)],
            AromaticSystemForm::default(),
        )],
        ..Default::default()
    }),
    Entity::AromaticSystem(AromaticSystemId(0)),
    AtomId(0),
    (Incidence::AromaticAtom, None)
)]
#[case::multicenter_atom(
    mol_dsl!(
        r#"{:atoms ["B" "B" "H"] :multicenter-bonds [{:atoms [0 2 1] :attrs "[1,1,0]"}]}"#
    ),
    Entity::MulticenterBond(MulticenterBondId(0)),
    AtomId(1),
    (Incidence::MulticenterAtom, Some(0))
)]
fn test_incidence_color(
    #[case] molecule: Molecule,
    #[case] owner: Entity,
    #[case] atom: AtomId,
    #[case] expected: (Incidence, Option<i64>),
) {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Constitution);
    let edge = incidence_graph
        .graph()
        .find_edge(
            incidence_graph.node_of(owner),
            incidence_graph.node_of(Entity::Atom(atom)),
        )
        .unwrap();
    let (incidence, electron_count) = incidence_color(&molecule, &incidence_graph, edge);

    assert_eq!((*incidence, electron_count), expected);
}

#[rstest]
#[case::atom(Entity::Atom(AtomId(0)), Entity::Atom(AtomId(1)), [true, false])]
#[case::bond(Entity::Bond(BondId(0)), Entity::Bond(BondId(1)), [true, false])]
#[case::dative(Entity::DativeBond(DativeBondId(0)), Entity::DativeBond(DativeBondId(1)), [true, false])]
#[case::aromatic(Entity::AromaticSystem(AromaticSystemId(0)), Entity::AromaticSystem(AromaticSystemId(1)), [true, false])]
#[case::multicenter(Entity::MulticenterBond(MulticenterBondId(0)), Entity::MulticenterBond(MulticenterBondId(1)), [true, false])]
#[case::noncovalent(Entity::NoncovalentBond(NoncovalentBondId(0)), Entity::NoncovalentBond(NoncovalentBondId(1)), [true, false])]
#[case::stereo_atom(Entity::StereoAtom(StereoAtomId(0)), Entity::StereoAtom(StereoAtomId(1)), [true, false])]
#[case::stereo_bond(Entity::StereoBond(StereoBondId(0)), Entity::StereoBond(StereoBondId(1)), [true, false])]
fn test_initial_colors_constraints(
    constrained_color_molecule: Molecule,
    #[case] lhs: Entity,
    #[case] rhs: Entity,
    #[case] expected: [bool; 2],
    #[values(false, true)] include_constraints: bool,
) {
    let molecule = constrained_color_molecule;
    let original = molecule.clone();
    let graph = molecule.incidence_graph(IncidenceLevel::Full);
    let colors = initial_colors(
        &molecule,
        &graph,
        DescriptionLevel::Structure,
        include_constraints,
    )
    .unwrap();
    assert_eq!(
        colors.entities[graph.node_of(lhs).index()] == colors.entities[graph.node_of(rhs).index()],
        expected[usize::from(include_constraints)]
    );
    assert_eq!(molecule, original);
}

#[rstest]
#[case::raw_element(AtomForm::from_element(Element::C), AtomForm::new(ElementForm::lit_set([Element::C])), [false, false])]
#[case::raw_charge(AtomForm::from_element(Element::C).with_charge(1), AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([1])), [false, false])]
#[case::raw_constraint(
    AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(4)),
    AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::Valence(NumForm::lit_set([4]))),
    [true, false]
)]
#[case::contradictory_constraint(
    AtomForm::from_element(Element::C),
    AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::Valence(NumForm::lit_set([]))),
    [true, false]
)]
fn test_initial_colors_raw(
    #[case] lhs: AtomForm,
    #[case] rhs: AtomForm,
    #[case] expected: [bool; 2],
    #[values(false, true)] include_constraints: bool,
    #[values(DescriptionLevel::Structure, DescriptionLevel::Full)] level: DescriptionLevel,
) {
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![lhs, rhs],
        constraints: Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: None,
            sum: NumForm::lit_set([]),
        })
        .into(),
        ..Default::default()
    });
    let original = molecule.clone();
    let graph = molecule.incidence_graph(IncidenceLevel::Full);
    let colors = initial_colors(&molecule, &graph, level, include_constraints).unwrap();
    assert_eq!(
        colors.entities[graph.node_of(Entity::Atom(AtomId(0))).index()]
            == colors.entities[graph.node_of(Entity::Atom(AtomId(1))).index()],
        expected[usize::from(include_constraints)]
    );
    assert_eq!(molecule, original);
}

#[rstest]
#[case::raw_set(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::lit_set([0])), [Ok(true), Ok(true)])]
#[case::raw_term(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::term(StereoTerm::var_in("x", []))), [Ok(true), Ok(true)])]
#[case::contradiction(StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32).with_constraint(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::lit_set([]))), [Ok(true), Ok(false)])]
#[case::position(StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32).with_constraint(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(7)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) })), [Ok(true), Err(Contradiction)])]
#[case::permutation(StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32).with_constraint(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(3)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })), [Ok(true), Err(Contradiction)])]
fn test_initial_colors_stereo_atom(
    initial_color_molecule: Molecule,
    #[case] attributes: StereoAtomForm,
    #[case] expected: [Result<bool, Contradiction>; 2],
    #[values(false, true)] include_constraints: bool,
) {
    let mut molecule = initial_color_molecule;
    *molecule.stereo_atom_mut(StereoAtomId(0)).attributes_mut() = attributes;
    let original = molecule.clone();
    let graph = molecule.incidence_graph(IncidenceLevel::Full);
    let actual = initial_colors(
        &molecule,
        &graph,
        DescriptionLevel::Structure,
        include_constraints,
    )
    .map(|colors| {
        colors.entities[graph.node_of(Entity::StereoAtom(StereoAtomId(0))).index()]
            == colors.entities[graph.node_of(Entity::StereoAtom(StereoAtomId(1))).index()]
    });
    assert_eq!(actual, expected[usize::from(include_constraints)]);
    assert_eq!(molecule, original);
}

#[rstest]
#[case::raw_set(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::lit_set([0])), [Ok(true), Ok(true)])]
#[case::raw_term(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::term(StereoTerm::var_in("x", []))), [Ok(true), Ok(true)])]
#[case::contradiction(StereoBondForm::new(StereoKind::CisTrans, 0_u32).with_constraint(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::lit_set([]))), [Ok(true), Ok(false)])]
#[case::position(StereoBondForm::new(StereoKind::CisTrans, 0_u32).with_constraint(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(7)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) })), [Ok(true), Err(Contradiction)])]
#[case::permutation(StereoBondForm::new(StereoKind::CisTrans, 0_u32).with_constraint(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(3)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })), [Ok(true), Err(Contradiction)])]
fn test_initial_colors_stereo_bond(
    initial_color_molecule: Molecule,
    #[case] attributes: StereoBondForm,
    #[case] expected: [Result<bool, Contradiction>; 2],
    #[values(false, true)] include_constraints: bool,
) {
    let mut molecule = initial_color_molecule;
    *molecule.stereo_bond_mut(StereoBondId(0)).attributes_mut() = attributes;
    let original = molecule.clone();
    let graph = molecule.incidence_graph(IncidenceLevel::Full);
    let actual = initial_colors(
        &molecule,
        &graph,
        DescriptionLevel::Structure,
        include_constraints,
    )
    .map(|colors| {
        colors.entities[graph.node_of(Entity::StereoBond(StereoBondId(0))).index()]
            == colors.entities[graph.node_of(Entity::StereoBond(StereoBondId(1))).index()]
    });
    assert_eq!(actual, expected[usize::from(include_constraints)]);
    assert_eq!(molecule, original);
}

#[rstest]
#[case::atom_constraint(Entity::Atom(AtomId(0)), Entity::Atom(AtomId(1)), false)]
#[case::atom_element(Entity::Atom(AtomId(0)), Entity::Atom(AtomId(2)), false)]
#[case::atom_charge(Entity::Atom(AtomId(0)), Entity::Atom(AtomId(3)), false)]
#[case::normalized_bond(Entity::Bond(BondId(0)), Entity::Bond(BondId(1)), true)]
#[case::bond_order(Entity::Bond(BondId(0)), Entity::Bond(BondId(2)), false)]
#[case::normalized_dative(
    Entity::DativeBond(DativeBondId(0)),
    Entity::DativeBond(DativeBondId(1)),
    true
)]
#[case::dative_order(
    Entity::DativeBond(DativeBondId(0)),
    Entity::DativeBond(DativeBondId(2)),
    false
)]
#[case::aromatic_electrons(
    Entity::AromaticSystem(AromaticSystemId(0)),
    Entity::AromaticSystem(AromaticSystemId(1)),
    true
)]
#[case::aromatic_charge(
    Entity::AromaticSystem(AromaticSystemId(0)),
    Entity::AromaticSystem(AromaticSystemId(2)),
    false
)]
#[case::multicenter_electrons(
    Entity::MulticenterBond(MulticenterBondId(0)),
    Entity::MulticenterBond(MulticenterBondId(1)),
    true
)]
#[case::noncovalent_kind_equal(
    Entity::NoncovalentBond(NoncovalentBondId(0)),
    Entity::NoncovalentBond(NoncovalentBondId(1)),
    true
)]
#[case::noncovalent_kind_distinct(
    Entity::NoncovalentBond(NoncovalentBondId(0)),
    Entity::NoncovalentBond(NoncovalentBondId(2)),
    false
)]
#[case::stereo_atom_configuration(
    Entity::StereoAtom(StereoAtomId(0)),
    Entity::StereoAtom(StereoAtomId(1)),
    true
)]
#[case::stereo_atom_kind(
    Entity::StereoAtom(StereoAtomId(0)),
    Entity::StereoAtom(StereoAtomId(2)),
    false
)]
#[case::stereo_bond_configuration(
    Entity::StereoBond(StereoBondId(0)),
    Entity::StereoBond(StereoBondId(1)),
    true
)]
#[case::entity_kind(Entity::Atom(AtomId(0)), Entity::Bond(BondId(0)), false)]
fn test_initial_colors(
    initial_color_molecule: Molecule,
    #[case] lhs: Entity,
    #[case] rhs: Entity,
    #[case] expected: bool,
) {
    let molecule = initial_color_molecule
        .normalize()
        .expect("fixed molecule normalizes");
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let colors = initial_colors(
        &molecule,
        &incidence_graph,
        DescriptionLevel::Structure,
        true,
    )
    .expect("fixed molecule has colors");

    assert_eq!(
        colors.entities[incidence_graph.node_of(lhs).index()]
            == colors.entities[incidence_graph.node_of(rhs).index()],
        expected
    );
}

#[rstest]
#[case::raw_atom(Entity::Atom(AtomId(0)), Entity::Atom(AtomId(1)), Ordering::Less)]
#[case::raw_bond(Entity::Bond(BondId(0)), Entity::Bond(BondId(1)), Ordering::Less)]
#[case::raw_dative(
    Entity::DativeBond(DativeBondId(0)),
    Entity::DativeBond(DativeBondId(1)),
    Ordering::Less
)]
#[case::aromatic_frame(
    Entity::AromaticSystem(AromaticSystemId(0)),
    Entity::AromaticSystem(AromaticSystemId(1)),
    Ordering::Equal
)]
#[case::multicenter_frame(
    Entity::MulticenterBond(MulticenterBondId(0)),
    Entity::MulticenterBond(MulticenterBondId(1)),
    Ordering::Equal
)]
#[case::noncovalent_kind(
    Entity::NoncovalentBond(NoncovalentBondId(0)),
    Entity::NoncovalentBond(NoncovalentBondId(2)),
    Ordering::Less
)]
#[case::stereo_atom_coset(
    Entity::StereoAtom(StereoAtomId(0)),
    Entity::StereoAtom(StereoAtomId(1)),
    Ordering::Equal
)]
#[case::stereo_bond_coset(
    Entity::StereoBond(StereoBondId(0)),
    Entity::StereoBond(StereoBondId(1)),
    Ordering::Equal
)]
fn test_initial_colors_order(
    initial_color_molecule: Molecule,
    #[case] lhs: Entity,
    #[case] rhs: Entity,
    #[case] expected: Ordering,
    #[values(false, true)] include_constraints: bool,
) {
    let molecule = initial_color_molecule;
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let colors = initial_colors(
        &molecule,
        &incidence_graph,
        DescriptionLevel::Structure,
        include_constraints,
    )
    .unwrap();
    assert_eq!(
        colors.entities[incidence_graph.node_of(lhs).index()]
            .cmp(&colors.entities[incidence_graph.node_of(rhs).index()]),
        expected,
    );
}

#[rstest]
fn test_initial_colors_incidence(
    initial_color_molecule: Molecule,
    #[values(false, true)] include_constraints: bool,
) {
    let molecule = initial_color_molecule
        .normalize()
        .expect("fixed molecule normalizes");
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let colors = initial_colors(
        &molecule,
        &incidence_graph,
        DescriptionLevel::Structure,
        include_constraints,
    )
    .expect("fixed molecule has colors");
    let contributions = [1, 2, 2, 1, 1, 2, 1, 2, 2, 1];
    let incidences = incidence_graph
        .incidences()
        .map(|(edge, incidence)| {
            let contribution = edge
                .index()
                .checked_sub(14)
                .and_then(|position| contributions.get(position));
            ((incidence, contribution), colors.incidences[edge.index()])
        })
        .collect::<Vec<_>>();

    for (lhs, lhs_color) in &incidences {
        for (rhs, rhs_color) in &incidences {
            assert_eq!(lhs_color == rhs_color, lhs == rhs);
        }
    }
    for entity_color in &colors.entities {
        for (_, incidence_color) in &incidences {
            assert_ne!(entity_color, incidence_color);
        }
    }
}

#[rstest]
#[case::atom_kind(
    Entity::StereoAtom(StereoAtomId(0)),
    StereoConfigurationForm::kinded(StereoKind::CisTrans, 0_u32)
)]
#[case::atom_frame(
    Entity::StereoAtom(StereoAtomId(0)),
    StereoConfigurationForm::kinded(StereoKind::Octahedral, 0_u32)
)]
#[case::bond_kind(
    Entity::StereoBond(StereoBondId(0)),
    StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0_u32)
)]
#[case::bond_coset(
    Entity::StereoBond(StereoBondId(0)),
    StereoConfigurationForm::kinded(StereoKind::CisTrans, 2_u32)
)]
fn test_initial_colors_configuration_error(
    initial_color_molecule: Molecule,
    #[case] entity: Entity,
    #[case] configuration: StereoConfigurationForm,
    #[values(false, true)] include_constraints: bool,
) {
    let mut molecule = initial_color_molecule;
    match entity {
        Entity::StereoAtom(id) => {
            molecule.stereo_atom_mut(id).attributes_mut().configuration = configuration
        }
        Entity::StereoBond(id) => {
            molecule.stereo_bond_mut(id).attributes_mut().configuration = configuration
        }
        _ => unreachable!(),
    }
    let graph = molecule.incidence_graph(IncidenceLevel::Full);
    assert_eq!(
        initial_colors(
            &molecule,
            &graph,
            DescriptionLevel::Structure,
            include_constraints
        ),
        Err(Contradiction)
    );
}

#[rstest]
#[case::aromatic_electron_count_length(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 3],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::from_electrons(vec![2, 1]),
        )],
        ..Default::default()
    }),
    DescriptionLevel::Constitution
)]
#[case::multicenter_electron_count_length(
    Molecule::from_entries(MoleculeEntries {
        atoms: [Element::B, Element::H, Element::B].map(AtomForm::from_element).into(),
        multicenter: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            MulticenterBondForm::from_electrons(vec![1, 0, 1, 0]),
        )],
        ..Default::default()
    }),
    DescriptionLevel::Constitution
)]
#[case::stereo_atom_coset_range(
    Molecule::from_entries(MoleculeEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .map(AtomForm::from_element)
            .into(),
        bonds: (1..=4)
            .map(|ligand| (AtomId(0), AtomId(ligand), BondForm::from_order(1)))
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            (1..=4)
                .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                .collect(),
            StereoAtomForm::new(StereoKind::Tetrahedral, 2u32),
        )],
        ..Default::default()
    }),
    DescriptionLevel::Structure
)]
fn test_initial_colors_error(
    #[case] molecule: Molecule,
    #[case] level: DescriptionLevel,
    #[values(false, true)] include_constraints: bool,
) {
    let incidence_graph = molecule.incidence_graph(match level {
        DescriptionLevel::Topology => IncidenceLevel::Topology,
        DescriptionLevel::Constitution => IncidenceLevel::Constitution,
        DescriptionLevel::Structure | DescriptionLevel::Full => IncidenceLevel::Full,
    });

    assert_eq!(
        initial_colors(&molecule, &incidence_graph, level, include_constraints),
        Err(Contradiction)
    );
}

#[rstest]
#[case::dichlorodifluoroethane_chiral(dichlorodifluoroethane_chiral(), vec![Stereogenicity::Stereogenic, Stereogenicity::Stereogenic])]
#[case::dichlorodifluoroethane_meso(dichlorodifluoroethane_meso(), vec![Stereogenicity::Stereogenic, Stereogenicity::Stereogenic])]
#[case::dichlorodifluoroethane_chiral_renumbered(dichlorodifluoroethane_chiral_renumbered(), vec![Stereogenicity::Stereogenic, Stereogenicity::Stereogenic])]
#[case::dichlorodifluoroethane_meso_renumbered(dichlorodifluoroethane_meso_renumbered(), vec![Stereogenicity::Stereogenic, Stereogenicity::Stereogenic])]
#[case::trichloropentane_meso(trichloropentane_meso(), vec![Stereogenicity::Stereogenic, Stereogenicity::Stereogenic, Stereogenicity::Stereogenic])]
#[case::trichloropentane_meso_inverted_c3(trichloropentane_meso_inverted_c3(), vec![Stereogenicity::Stereogenic, Stereogenicity::Stereogenic, Stereogenicity::Stereogenic])]
#[case::trichloropentane_chiral(trichloropentane_chiral(), vec![Stereogenicity::Stereogenic, Stereogenicity::Prostereogenic, Stereogenicity::Stereogenic])]
#[case::dichlorocyclobutane_cis(dichlorocyclobutane_cis(), vec![Stereogenicity::Stereogenic, Stereogenicity::Stereogenic])]
#[case::dichlorocyclobutane_trans(dichlorocyclobutane_trans(), vec![Stereogenicity::Stereogenic, Stereogenicity::Stereogenic])]
#[case::chlorofluoropropane(chlorofluoropropane(), vec![Stereogenicity::Stereogenic, Stereogenicity::Prostereogenic])]
#[case::chlorofluoropropane_without_c1(chlorofluoropropane_without_c1(), vec![Stereogenicity::Prostereogenic])]
#[case::trimethylcyclohexane_all_cis(trimethylcyclohexane_all_cis(), vec![Stereogenicity::Stereogenic, Stereogenicity::Stereogenic, Stereogenicity::Stereogenic])]
#[case::trimethylcyclohexane_cis_trans(trimethylcyclohexane_cis_trans(), vec![Stereogenicity::Stereogenic, Stereogenicity::Prostereogenic, Stereogenicity::Prostereogenic])]
fn test_exhaustive_stereogenicities(
    #[case] molecule: Molecule,
    #[case] expected: Vec<Stereogenicity>,
) {
    assert_eq!(exhaustive_stereogenicities(&molecule), expected);
}

#[rstest]
#[case::dative(
    mol_dsl!(r#"{:atoms ["B" "B"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#),
    vec![vec![0, 1]], vec![0, 1]
)]
#[case::aromatic(
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [0 4 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#),
    vec![vec![0, 1, 2, 3, 4], vec![0, 4, 3, 2, 1]], vec![0, 1, 2, 2, 1]
)]
#[case::multicenter(
    mol_dsl!(r#"{:atoms ["C" "C" "C"]
        :multicenter-bonds [{:atoms [0 1 2] :attrs "[2,1,1]"}]}"#),
    vec![vec![0, 1, 2], vec![0, 2, 1]], vec![0, 1, 1]
)]
fn test_exhaustive_atom_remappings(
    #[case] molecule: Molecule,
    #[case] expected: Vec<Vec<u32>>,
    #[case] expected_orbits: Vec<u32>,
) {
    assert_eq!(molecule.check_integrity(), Ok(()));
    let mappings = exhaustive_atom_remappings(&molecule);
    assert_eq!(mappings, expected);
    assert_eq!(
        exhaustive_symmetry(&molecule),
        (false, expected_orbits.clone(), expected_orbits)
    );
}

#[rstest]
#[case::chiral_proper(dichlorodifluoroethane_chiral(), false,
    vec![vec![0, 1, 2, 3, 4, 5, 6, 7], vec![1, 0, 5, 6, 7, 2, 3, 4]])]
#[case::chiral_improper(dichlorodifluoroethane_chiral(), true, vec![])]
#[case::meso_proper(dichlorodifluoroethane_meso(), false, vec![vec![0, 1, 2, 3, 4, 5, 6, 7]])]
#[case::meso_improper(dichlorodifluoroethane_meso(), true, vec![vec![1, 0, 5, 6, 7, 2, 3, 4]])]
fn test_preserves_stereo_configuration(
    #[case] molecule: Molecule,
    #[case] mirrored: bool,
    #[case] expected: Vec<Vec<u32>>,
) {
    let records = stereo_records(&molecule);
    let accepted = exhaustive_atom_remappings(&molecule)
        .into_iter()
        .filter(|images| preserves_stereo_configuration(&records, images, mirrored, None))
        .collect::<Vec<_>>();
    assert_eq!(accepted, expected);
}

#[rstest]
#[case::chlorofluoropropane(chlorofluoropropane(), 1, 2, 3, Topicity::Diastereotopic)]
#[case::chlorofluoropropane_without_c1(
    chlorofluoropropane_without_c1(),
    0,
    2,
    3,
    Topicity::Enantiotopic
)]
#[case::trichloropentane_meso(trichloropentane_meso(), 1, 0, 1, Topicity::Enantiotopic)]
#[case::trichloropentane_chiral(trichloropentane_chiral(), 1, 0, 1, Topicity::Diastereotopic)]
#[case::distinct(chlorofluoropropane(), 1, 0, 1, Topicity::Distinct)]
#[case::methane(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::H),
                    AtomForm::from_element(Element::H), AtomForm::from_element(Element::H), AtomForm::from_element(Element::H)],
        bonds: (1..=4).map(|atom| (AtomId(0), AtomId(atom), BondForm::from_order(1))).collect(),
        stereo_atoms: vec![(AtomId(0), (1..=4).map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom)).collect(),
            StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32))],
        ..Default::default()
    }), 0, 0, 1, Topicity::Homotopic
)]
fn test_exhaustive_topicity(
    #[case] molecule: Molecule,
    #[case] record: usize,
    #[case] first: usize,
    #[case] second: usize,
    #[case] expected: Topicity,
) {
    assert_eq!(
        exhaustive_topicity(&molecule, record, first, second),
        expected
    );
}

#[rstest]
fn test_graph_symmetry_config_default() {
    let config = GraphSymmetryConfig::default();
    assert_eq!(config.automorphism_algorithm, AutomorphismAlgorithm::Nauty);
    assert_eq!(config.max_classification_trials, 4096);
    assert_eq!(config.include_constraints, false);
}
