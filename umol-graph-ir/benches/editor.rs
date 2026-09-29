//! Current editor costs, including creation and publication.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use umol_graph_ir::ir::{
    AromaticSystemForm, AromaticSystemHandle, AromaticSystemId, AtomFieldChange, AtomForm,
    AtomHandle, AtomId, BondForm, BondHandle, BondId, DativeBondForm, DativeBondHandle,
    DativeBondId, Edit, Edits, EntityKind, FrameTransport, Molecule, MoleculeEntries,
    MulticenterBondForm, MulticenterBondHandle, MulticenterBondId, NoncovalentBondForm,
    NoncovalentBondHandle, NoncovalentBondId, NumForm, StereoAtomForm, StereoAtomHandle,
    StereoAtomId, StereoBondForm, StereoBondHandle, StereoBondId, StereoKind, StereoLigand,
    StereoLigandKind,
};
use umol_perm::Permutation;

fn entries(size: usize, dense: bool) -> MoleculeEntries {
    MoleculeEntries {
        atoms: vec![AtomForm::default(); size],
        bonds: (1..size)
            .map(|index| {
                (
                    AtomId((index - 1) as u32),
                    AtomId(index as u32),
                    BondForm::default(),
                )
            })
            .collect(),
        dative: if dense {
            (0..8)
                .map(|index| {
                    (
                        vec![AtomId((8 * index) as u32)],
                        AtomId((8 * index + 1) as u32),
                        DativeBondForm::default(),
                    )
                })
                .collect()
        } else {
            Vec::new()
        },
        aromatic: if dense {
            (0..8)
                .map(|index| {
                    (
                        (0..4)
                            .map(|offset| AtomId((8 * index + offset) as u32))
                            .collect(),
                        AromaticSystemForm::default(),
                    )
                })
                .collect()
        } else {
            Vec::new()
        },
        ..Default::default()
    }
}

fn changes(count: usize) -> Edits {
    (0..count)
        .map(|index| Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(index as u32)),
            change: AtomFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(1),
            },
        })
        .collect()
}

fn bench_editor(c: &mut Criterion) {
    let mut group = c.benchmark_group("molecule_editor");
    for (name, size, dense, edits) in [
        ("sparse8", 8_usize, false, 1_usize),
        ("dense80", 80, true, 8),
    ] {
        let source_entries = entries(size, dense);
        let source = Molecule::from_entries(source_entries.clone());
        for (ownership, unique) in [("unique", true), ("shared", false)] {
            let input = || {
                if unique {
                    Molecule::from_entries(source_entries.clone())
                } else {
                    source.clone()
                }
            };
            let id = |operation| BenchmarkId::new(format!("{name}/{ownership}"), operation);
            group.bench_function(id("direct"), |b| {
                b.iter_batched(
                    input,
                    |molecule| {
                        let mut editor = molecule.edit();
                        for index in 0..edits {
                            editor
                                .atom_mut(AtomId(index as u32))
                                .attributes_mut()
                                .charge = NumForm::Lit(1);
                        }
                        black_box(editor.try_build().unwrap())
                    },
                    BatchSize::SmallInput,
                );
            });
            group.bench_function(id("apply"), |b| {
                b.iter_batched(
                    || (input(), changes(edits)),
                    |(molecule, changes)| black_box(molecule.apply(changes).unwrap()),
                    BatchSize::SmallInput,
                );
            });
            group.bench_function(id("transact"), |b| {
                b.iter_batched(
                    || (input(), changes(edits)),
                    |(mut molecule, changes)| {
                        molecule.transact([changes]).unwrap();
                        black_box(molecule)
                    },
                    BatchSize::SmallInput,
                );
            });
        }
    }
    group.finish();
}

fn removal_entries() -> MoleculeEntries {
    MoleculeEntries {
        atoms: vec![AtomForm::default(); 5],
        bonds: [(0, 1), (0, 2), (0, 3), (1, 4)]
            .map(|(first, second)| (AtomId(first), AtomId(second), BondForm::from_order(1)))
            .into(),
        dative: vec![(
            vec![AtomId(0), AtomId(2), AtomId(3)],
            AtomId(1),
            DativeBondForm::from_order(1),
        )],
        aromatic: vec![(
            vec![AtomId(0), AtomId(2), AtomId(3)],
            AromaticSystemForm::from_electrons(vec![1, 2, 3]),
        )],
        multicenter: vec![(
            vec![AtomId(0), AtomId(1), AtomId(4)],
            MulticenterBondForm::from_electrons(vec![2, 1, 1]),
        )],
        noncovalent: vec![([AtomId(2), AtomId(4)], NoncovalentBondForm::default())],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
        )],
        stereo_bonds: vec![(
            BondId(0),
            vec![
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoBondForm::new(StereoKind::CisTrans, 0_u32),
        )],
        ..Default::default()
    }
}

fn removal(kind: EntityKind, reordered: bool) -> Edits {
    let mut entries = removal_entries();
    let edit = match kind {
        EntityKind::DativeBond => {
            let (mut donors, acceptor, attributes) = entries.dative.remove(0);
            if reordered {
                donors.reverse();
            }
            Edit::RemoveDativeBonds {
                removes: vec![(
                    DativeBondHandle::Id(DativeBondId(0)),
                    donors.into_iter().map(AtomHandle::Id).collect(),
                    AtomHandle::Id(acceptor),
                    attributes,
                )],
            }
        }
        EntityKind::AromaticSystem => {
            let (mut atoms, attributes) = entries.aromatic.remove(0);
            let attributes = if reordered {
                atoms.reverse();
                AromaticSystemForm::from_electrons(vec![3, 2, 1])
            } else {
                attributes
            };
            Edit::RemoveAromaticSystems {
                removes: vec![(
                    AromaticSystemHandle::Id(AromaticSystemId(0)),
                    atoms.into_iter().map(AtomHandle::Id).collect(),
                    attributes,
                )],
            }
        }
        EntityKind::MulticenterBond => {
            let (mut atoms, attributes) = entries.multicenter.remove(0);
            let attributes = if reordered {
                atoms.reverse();
                MulticenterBondForm::from_electrons(vec![1, 1, 2])
            } else {
                attributes
            };
            Edit::RemoveMulticenterBonds {
                removes: vec![(
                    MulticenterBondHandle::Id(MulticenterBondId(0)),
                    atoms.into_iter().map(AtomHandle::Id).collect(),
                    attributes,
                )],
            }
        }
        EntityKind::NoncovalentBond => {
            let (mut atoms, attributes) = entries.noncovalent.remove(0);
            if reordered {
                atoms.reverse();
            }
            Edit::RemoveNoncovalentBonds {
                removes: vec![(
                    NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                    atoms.map(AtomHandle::Id),
                    attributes,
                )],
            }
        }
        EntityKind::StereoAtom => {
            let (site, mut ligands, mut attributes) = entries.stereo_atoms.remove(0);
            if reordered {
                let action = Permutation::from_image(&[1, 0, 2, 3]);
                ligands = action.act(&ligands);
                attributes = attributes.reframe_by(&action).unwrap();
            }
            Edit::RemoveStereoAtoms {
                removes: vec![(
                    StereoAtomHandle::Id(StereoAtomId(0)),
                    AtomHandle::Id(site),
                    ligands
                        .into_iter()
                        .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                        .collect(),
                    attributes,
                )],
            }
        }
        EntityKind::StereoBond => {
            let (site, mut ligands, mut attributes) = entries.stereo_bonds.remove(0);
            if reordered {
                let action = Permutation::from_image(&[1, 0, 2, 3]);
                ligands = action.act(&ligands);
                attributes = attributes.reframe_by(&action).unwrap();
            }
            Edit::RemoveStereoBonds {
                removes: vec![(
                    StereoBondHandle::Id(StereoBondId(0)),
                    BondHandle::Id(site),
                    ligands
                        .into_iter()
                        .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                        .collect(),
                    attributes,
                )],
            }
        }
        EntityKind::Atom | EntityKind::Bond => unreachable!("overlay removal benchmark"),
    };
    Edits::from_iter([edit])
}

fn bench_removal(c: &mut Criterion) {
    let mut group = c.benchmark_group("molecule_editor/removal");
    for kind in [
        EntityKind::DativeBond,
        EntityKind::AromaticSystem,
        EntityKind::MulticenterBond,
        EntityKind::NoncovalentBond,
        EntityKind::StereoAtom,
        EntityKind::StereoBond,
    ] {
        for (frame, reordered) in [("stored", false), ("reordered", true)] {
            group.bench_function(BenchmarkId::new(kind.to_string(), frame), |b| {
                b.iter_batched(
                    || {
                        (
                            Molecule::from_entries(removal_entries()).edit(),
                            removal(kind, reordered),
                        )
                    },
                    |(editor, edits)| black_box(editor.apply(edits).unwrap()),
                    BatchSize::SmallInput,
                );
            });
        }
    }
    group.finish();
}

criterion_group!(editor, bench_editor, bench_removal);
criterion_main!(editor);
