//! Checked aggregate construction and projection benchmarks.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use umol_graph_ir::ir::{
    AromaticSystemForm, AtomConstraintForm, AtomForm, AtomId, BondForm, BondId, Constraint,
    ConstraintSpan, DativeBondDelta, DativeBondForm, DativeBondId, Delta, Deltas, EntitySpan,
    Molecule, MoleculeEntries, MulticenterBondForm, NoncovalentBondForm, Reaction, ReactionSpan,
    ReactionSpanEntries, StereoAtomForm, StereoBondForm, StereoLigand, StereoLigandKind,
};

fn chain_entries(size: usize, constrained: bool) -> MoleculeEntries {
    MoleculeEntries {
        atoms: vec![AtomForm::default(); size],
        bonds: (1..size)
            .map(|i| {
                (
                    AtomId((i - 1) as u32),
                    AtomId(i as u32),
                    BondForm::default(),
                )
            })
            .collect(),
        constraints: if constrained {
            vec![Constraint::Atom(AtomId(0), AtomConstraintForm::valence(3))]
        } else {
            Vec::new()
        }
        .into(),
        ..Default::default()
    }
}

fn overlay_entries(constrained: bool) -> MoleculeEntries {
    let mut entries = chain_entries(80, constrained);
    entries.dative = (0..20)
        .map(|i| (vec![AtomId(i)], AtomId(i + 20), DativeBondForm::default()))
        .collect();
    entries.aromatic = (0..10)
        .map(|i| {
            (
                (0..4).map(|offset| AtomId(4 * i + offset)).collect(),
                AromaticSystemForm::default(),
            )
        })
        .collect();
    entries.multicenter = (0..10)
        .map(|i| {
            (
                (0..3).map(|offset| AtomId(3 * i + offset)).collect(),
                MulticenterBondForm::default(),
            )
        })
        .collect();
    entries.noncovalent = (0..10)
        .map(|i| ([AtomId(i), AtomId(i + 40)], NoncovalentBondForm::default()))
        .collect();
    entries.stereo_atoms = (0..10)
        .map(|i| {
            let site = AtomId(3 * i + 1);
            (
                site,
                vec![
                    StereoLigand::new(AtomId(site.0 - 1), StereoLigandKind::Atom),
                    StereoLigand::new(site, StereoLigandKind::ImplicitHydrogen),
                ],
                StereoAtomForm::default(),
            )
        })
        .collect();
    entries.stereo_bonds = (0..10)
        .map(|i| {
            let first = AtomId(3 * i);
            let second = AtomId(first.0 + 1);
            (
                BondId(first.0),
                vec![
                    StereoLigand::new(first, StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(first, StereoLigandKind::LonePair),
                    StereoLigand::new(second, StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(second, StereoLigandKind::LonePair),
                ],
                StereoBondForm::default(),
            )
        })
        .collect();
    entries
}

fn removal_deltas(count: u32) -> Deltas {
    (0..count)
        .map(|i| {
            Delta::DativeBond(DativeBondDelta::Remove {
                id: DativeBondId(i),
                donors: vec![AtomId(i)],
                acceptor: AtomId(i + 20),
                attributes: DativeBondForm::default(),
            })
        })
        .collect()
}

fn span_entries(size: usize, interleaved: bool, constrained: bool) -> ReactionSpanEntries {
    let mut entries = ReactionSpanEntries {
        atoms: (0..size)
            .map(|i| {
                let form = AtomForm::default();
                if interleaved {
                    match i % 4 {
                        1 => EntitySpan::Added(form),
                        3 => EntitySpan::Removed(form),
                        _ => EntitySpan::Unchanged(form),
                    }
                } else {
                    EntitySpan::Unchanged(form)
                }
            })
            .collect(),
        constraints: if constrained {
            vec![ConstraintSpan::Unchanged(Constraint::Atom(
                AtomId(0),
                AtomConstraintForm::valence(3),
            ))]
        } else {
            Vec::new()
        },
        ..Default::default()
    };
    if interleaved {
        for i in (2..size).step_by(2) {
            entries.bonds.push((
                AtomId((i - 2) as u32),
                AtomId(i as u32),
                EntitySpan::Unchanged(BondForm::default()),
            ));
        }
        for i in (1..size).step_by(2) {
            entries.bonds.push((
                AtomId(0),
                AtomId(i as u32),
                if i % 4 == 1 {
                    EntitySpan::Added(BondForm::default())
                } else {
                    EntitySpan::Removed(BondForm::default())
                },
            ));
        }
    } else {
        entries.bonds = (1..size)
            .map(|i| {
                (
                    AtomId((i - 1) as u32),
                    AtomId(i as u32),
                    EntitySpan::Unchanged(BondForm::default()),
                )
            })
            .collect();
    }
    entries
}

fn benchmark_molecule(c: &mut Criterion) {
    let mut group = c.benchmark_group("integrity/molecule");
    for size in [12, 80] {
        for constrained in [false, true] {
            let entries = chain_entries(size, constrained);
            Molecule::try_from_entries(entries.clone()).expect("valid benchmark molecule");
            let name = format!("{size}/chain/constraints_{constrained}");
            group.bench_function(name, |b| {
                b.iter_batched(
                    || entries.clone(),
                    |entries| {
                        black_box(
                            Molecule::try_from_entries(black_box(entries))
                                .expect("valid benchmark molecule"),
                        )
                    },
                    BatchSize::SmallInput,
                )
            });
        }
    }
    for constrained in [false, true] {
        let entries = overlay_entries(constrained);
        Molecule::try_from_entries(entries.clone()).expect("valid benchmark overlays");
        group.bench_function(format!("80/overlays/constraints_{constrained}"), |b| {
            b.iter_batched(
                || entries.clone(),
                |entries| {
                    black_box(
                        Molecule::try_from_entries(black_box(entries))
                            .expect("valid benchmark overlays"),
                    )
                },
                BatchSize::SmallInput,
            )
        });
    }
    group.finish();
}

fn benchmark_reaction(c: &mut Criterion) {
    let mut group = c.benchmark_group("integrity/reaction");
    for constrained in [false, true] {
        let lhs = Molecule::try_from_entries(overlay_entries(constrained))
            .expect("valid benchmark reaction lhs");
        for removals in [0, 1, 20] {
            let deltas = removal_deltas(removals);
            Reaction::try_new(lhs.clone(), deltas.clone()).expect("valid benchmark reaction");
            group.bench_function(
                BenchmarkId::new(
                    format!("80/constraints_{constrained}"),
                    format!("removals_{removals}"),
                ),
                |b| {
                    b.iter_batched(
                        || (lhs.clone(), deltas.clone()),
                        |(lhs, deltas)| {
                            black_box(
                                Reaction::try_new(black_box(lhs), black_box(deltas))
                                    .expect("valid benchmark reaction"),
                            )
                        },
                        BatchSize::SmallInput,
                    )
                },
            );
        }
    }
    group.finish();
}

fn benchmark_span(c: &mut Criterion) {
    let mut group = c.benchmark_group("integrity/reaction_span");
    for (size, interleaved) in [(12, false), (80, false), (80, true)] {
        for constrained in [false, true] {
            let entries = span_entries(size, interleaved, constrained);
            let span =
                ReactionSpan::try_from_entries(entries.clone()).expect("valid benchmark span");
            let name = format!("{size}/interleaved_{interleaved}/constraints_{constrained}");
            group.bench_function(BenchmarkId::new("construct", &name), |b| {
                b.iter_batched(
                    || entries.clone(),
                    |entries| {
                        black_box(
                            ReactionSpan::try_from_entries(black_box(entries))
                                .expect("valid benchmark span"),
                        )
                    },
                    BatchSize::SmallInput,
                )
            });
            group.bench_function(BenchmarkId::new("lhs", &name), |b| {
                b.iter(|| black_box(black_box(&span).lhs()))
            });
            group.bench_function(BenchmarkId::new("rhs", &name), |b| {
                b.iter(|| black_box(black_box(&span).rhs()))
            });
            group.bench_function(BenchmarkId::new("to_reaction", &name), |b| {
                b.iter(|| black_box(black_box(&span).to_reaction()))
            });
        }
    }
    group.finish();
}

criterion_group!(
    integrity,
    benchmark_molecule,
    benchmark_reaction,
    benchmark_span
);
criterion_main!(integrity);
