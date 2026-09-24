//! Current editor costs, including creation and publication.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use umol_graph_ir::ir::{
    AromaticSystemForm, AtomFieldChange, AtomForm, AtomHandle, AtomId, BondForm, DativeBondForm,
    Edit, Edits, Molecule, MoleculeEntries, NumForm,
};

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
                            editor.atom_mut(AtomId(index as u32)).attributes.charge =
                                NumForm::Lit(1);
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
                    |(molecule, changes)| {
                        let mut editor = molecule.edit();
                        let journal = editor.transact(changes).unwrap();
                        black_box(journal.undos().len());
                        black_box(editor.try_build().unwrap())
                    },
                    BatchSize::SmallInput,
                );
            });
        }
    }
    group.finish();
}

criterion_group!(editor, bench_editor);
criterion_main!(editor);
