//! Aggregate construction and clone costs. Construction excludes input preparation and output
//! destruction; clone measurements include dropping the clone while retaining the source.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use umol_chem::element::Element;
use umol_graph_ir::ir::{
    AromaticSystemForm, AtomDelta, AtomFieldChange, AtomForm, AtomId, BondForm, Constraint, Delta,
    Deltas, ElementForm, EntitySpan, Molecule, MoleculeConstraint, MoleculeEntries, NumForm,
    Reaction, ReactionSpan, ReactionSpanEntries,
};

fn molecule_cases() -> [(&'static str, MoleculeEntries); 2] {
    [
        (
            "benzene",
            MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); 6],
                bonds: (0..6)
                    .map(|i| (AtomId(i), AtomId((i + 1) % 6), BondForm::from_order(1)))
                    .collect(),
                aromatic: vec![(
                    (0..6).map(AtomId).collect(),
                    AromaticSystemForm::from_electrons(vec![1; 6]),
                )],
                ..Default::default()
            },
        ),
        (
            "chain_100",
            MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); 100],
                bonds: (1..100)
                    .map(|i| (AtomId(i - 1), AtomId(i), BondForm::from_order(1)))
                    .collect(),
                ..Default::default()
            },
        ),
    ]
}

fn bench_molecule(c: &mut Criterion) {
    let mut group = c.benchmark_group("aggregate_storage/molecule");
    for (name, entries) in molecule_cases() {
        let molecule = Molecule::from_entries(entries.clone());
        group.bench_function(BenchmarkId::new("from_entries", name), |b| {
            b.iter_batched(
                || entries.clone(),
                |entries| Molecule::from_entries(black_box(entries)),
                BatchSize::SmallInput,
            )
        });
        group.bench_function(BenchmarkId::new("clone", name), |b| {
            b.iter(|| black_box(&molecule).clone())
        });
    }
    let entries = MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 100],
        bonds: (1..100)
            .map(|i| (AtomId(i - 1), AtomId(i), BondForm::from_order(1)))
            .collect(),
        constraints: (0..100)
            .map(|i| {
                Constraint::Molecule(MoleculeConstraint::ChargeSum {
                    atoms: Some(vec![AtomId(i)]),
                    sum: NumForm::Lit(0),
                })
            })
            .collect(),
        ..Default::default()
    };
    let constrained = Molecule::from_entries(entries);
    group.bench_function("clone/constraints_100", |b| {
        b.iter(|| black_box(&constrained).clone())
    });
    group.finish();
}

fn bench_reaction(c: &mut Criterion) {
    let mut group = c.benchmark_group("aggregate_storage/reaction");
    for (name, entries) in molecule_cases() {
        let lhs = Molecule::from_entries(entries);
        let deltas = Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(0),
            change: AtomFieldChange::Element {
                old: ElementForm::Lit(Element::C),
                new: ElementForm::Lit(Element::N),
            },
        })]);
        let reaction = Reaction::new(lhs.clone(), deltas.clone());
        if name == "chain_100" {
            group.bench_function(BenchmarkId::new("new", name), |b| {
                b.iter_batched(
                    || (lhs.clone(), deltas.clone()),
                    |(lhs, deltas)| Reaction::new(black_box(lhs), black_box(deltas)),
                    BatchSize::SmallInput,
                )
            });
        }
        group.bench_function(BenchmarkId::new("clone", name), |b| {
            b.iter(|| black_box(&reaction).clone())
        });
    }
    group.finish();
}

fn bench_reaction_span(c: &mut Criterion) {
    let mut group = c.benchmark_group("aggregate_storage/reaction_span");
    for (name, entries) in molecule_cases() {
        let entries = ReactionSpanEntries {
            atoms: entries
                .atoms
                .into_iter()
                .map(EntitySpan::Unchanged)
                .collect(),
            bonds: entries
                .bonds
                .into_iter()
                .map(|(first, second, form)| (first, second, EntitySpan::Unchanged(form)))
                .collect(),
            aromatic: entries
                .aromatic
                .into_iter()
                .map(|(atoms, form)| (atoms, EntitySpan::Unchanged(form)))
                .collect(),
            ..Default::default()
        };
        let span = ReactionSpan::from_entries(entries.clone());
        group.bench_function(BenchmarkId::new("from_entries", name), |b| {
            b.iter_batched(
                || entries.clone(),
                |entries| ReactionSpan::from_entries(black_box(entries)),
                BatchSize::SmallInput,
            )
        });
        group.bench_function(BenchmarkId::new("clone", name), |b| {
            b.iter(|| black_box(&span).clone())
        });
    }
    group.finish();
}

criterion_group!(
    aggregate_storage,
    bench_molecule,
    bench_reaction,
    bench_reaction_span
);
criterion_main!(aggregate_storage);
