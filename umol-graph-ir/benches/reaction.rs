//! Molecule mutation, reaction matching, and matched-application benchmarks.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use umol_chem::element::Element;
use umol_graph_core::{
    Correspondence, GraphCorrespondence, NodeId, RelevantCycleEnumerationAlgorithm,
    SubgraphIsomorphismAlgorithm,
};
use umol_graph_ir::ir::{
    AromaticSystemDelta, AromaticSystemFieldChange, AromaticSystemForm, AromaticSystemId,
    AtomConstraintForm, AtomForm, AtomHandle, AtomId, BondForm, Constraint, ConstraintSpan, Delta,
    Deltas, Edits, ElectronCountsForm, EntitySpan, Molecule, MoleculeCorrespondence,
    MoleculeEntries, Reaction, ReactionSpan, ReactionSpanEntries, SubstructureMatchAlgorithm,
    SubstructureMatchConfig,
};

const MATCH_CONFIG: SubstructureMatchConfig = SubstructureMatchConfig {
    match_algorithm: SubstructureMatchAlgorithm::GraphAndOverlays,
    subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::Vf2,
    relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
};

fn application_case() -> (Reaction, Molecule, MoleculeCorrespondence) {
    let atoms = [
        Element::C,
        Element::N,
        Element::O,
        Element::F,
        Element::Cl,
        Element::Br,
    ]
    .into_iter()
    .map(AtomForm::from_element)
    .collect::<Vec<_>>();
    let bonds = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]
        .map(|(first, second)| (AtomId(first), AtomId(second), BondForm::from_order(1)))
        .to_vec();
    let reaction = Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            aromatic: vec![((0..6).map(AtomId).collect(), AromaticSystemForm::default())],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::AromaticSystem(AromaticSystemDelta::ModifyField {
            id: AromaticSystemId(0),
            change: AromaticSystemFieldChange::Electrons {
                old: ElectronCountsForm::Undetermined,
                new: ElectronCountsForm::Lit(vec![1, 2, 3, 5, 7, 11]),
            },
        })]),
    );
    let host = Molecule::from_entries(MoleculeEntries {
        atoms,
        bonds,
        aromatic: vec![(
            [5, 0, 1, 2, 3, 4].map(AtomId).to_vec(),
            AromaticSystemForm {
                electrons: ElectronCountsForm::Lit(vec![13, 17, 19, 23, 29, 31]),
                ..Default::default()
            },
        )],
        ..Default::default()
    });
    let correspondence = MoleculeCorrespondence::induce(
        reaction.lhs(),
        &host,
        Correspondence::from_images(
            &[
                AtomId(0),
                AtomId(1),
                AtomId(2),
                AtomId(3),
                AtomId(4),
                AtomId(5),
            ],
            6,
        ),
    )
    .expect("the benchmark atom correspondence describes the molecule pair");
    (reaction, host, correspondence)
}

fn reversal_case() -> Reaction {
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Modified {
                lhs: AtomForm::from_element(Element::O),
                rhs: AtomForm::from_element(Element::N),
            },
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Removed(AtomForm::from_element(Element::F)),
            EntitySpan::Added(AtomForm::from_element(Element::Cl)),
        ],
        bonds: vec![
            (
                AtomId(0),
                AtomId(1),
                EntitySpan::Unchanged(BondForm::from_order(1)),
            ),
            (
                AtomId(1),
                AtomId(2),
                EntitySpan::Unchanged(BondForm::from_order(1)),
            ),
            (
                AtomId(2),
                AtomId(3),
                EntitySpan::Removed(BondForm::from_order(1)),
            ),
            (
                AtomId(2),
                AtomId(4),
                EntitySpan::Added(BondForm::from_order(1)),
            ),
        ],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            EntitySpan::Modified {
                lhs: AromaticSystemForm::default(),
                rhs: AromaticSystemForm::from_electrons(vec![1, 1, 1]),
            },
        )],
        constraints: vec![
            ConstraintSpan::Unchanged(Constraint::Atom(AtomId(0), AtomConstraintForm::valence(3))),
            ConstraintSpan::Removed(Constraint::Atom(AtomId(3), AtomConstraintForm::valence(1))),
            ConstraintSpan::Added(Constraint::Atom(AtomId(4), AtomConstraintForm::valence(1))),
        ],
        ..Default::default()
    })
    .to_reaction()
}

fn benchmark_reaction(c: &mut Criterion) {
    let (reaction, host, correspondence) = application_case();
    let applications = reaction
        .apply(&host, MATCH_CONFIG)
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(applications.len(), 1);
    let reversal = reversal_case();
    let mut group = c.benchmark_group("reaction");
    group.bench_function("match_enumeration", |b| {
        b.iter(|| {
            black_box(
                reaction
                    .lhs()
                    .substructure_matches(black_box(&host), MATCH_CONFIG),
            )
        })
    });
    group.bench_function("apply_at", |b| {
        b.iter(|| black_box(reaction.apply_at(black_box(&host), black_box(&correspondence))))
    });
    group.bench_function("tracked_apply_at", |b| {
        b.iter(|| {
            black_box(reaction.tracked_apply_at(black_box(&host), black_box(&correspondence)))
        })
    });
    group.bench_function("apply_at_to_reaction", |b| {
        b.iter(|| {
            black_box(reaction.apply_at_to_reaction(black_box(&host), black_box(&correspondence)))
        })
    });
    group.bench_function("apply_at_to_reaction_span", |b| {
        b.iter(|| {
            black_box(
                reaction.apply_at_to_reaction_span(black_box(&host), black_box(&correspondence)),
            )
        })
    });
    group.bench_function("apply/all_matches", |b| {
        b.iter(|| {
            black_box(
                reaction
                    .apply(black_box(&host), MATCH_CONFIG)
                    .unwrap()
                    .collect::<Vec<_>>(),
            )
        })
    });
    group.bench_function("tracked_apply/all_matches", |b| {
        b.iter(|| {
            black_box(
                reaction
                    .tracked_apply(black_box(&host), MATCH_CONFIG)
                    .unwrap()
                    .collect::<Vec<_>>(),
            )
        })
    });
    group.bench_function("apply_to_reaction/all_matches", |b| {
        b.iter(|| {
            black_box(
                reaction
                    .apply_to_reaction(black_box(&host), MATCH_CONFIG)
                    .unwrap()
                    .collect::<Vec<_>>(),
            )
        })
    });
    group.bench_function("apply_to_reaction_span/all_matches", |b| {
        b.iter(|| {
            black_box(
                reaction
                    .apply_to_reaction_span(black_box(&host), MATCH_CONFIG)
                    .unwrap()
                    .collect::<Vec<_>>(),
            )
        })
    });
    group.bench_function("reverse", |b| {
        b.iter(|| black_box(black_box(&reversal).reverse()))
    });
    group.finish();
}

fn benchmark_mutation(c: &mut Criterion) {
    let mut group = c.benchmark_group("molecule_mutation");
    for size in [8usize, 64] {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); size],
            bonds: (0..size as u32 - 1)
                .map(|id| (AtomId(id), AtomId(id + 1), BondForm::from_order(1)))
                .collect(),
            ..Default::default()
        });
        let removed = [AtomId((size / 2) as u32)];
        let disconnected = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); size],
            bonds: (0..size as u32 - 1)
                .filter(|&id| id != size as u32 / 2 - 1)
                .map(|id| (AtomId(id), AtomId(id + 1), BondForm::from_order(1)))
                .collect(),
            ..Default::default()
        });
        let overlap = GraphCorrespondence::new(
            Correspondence::new(vec![(NodeId((size - 1) as u32), NodeId(0))], size, size).unwrap(),
            Correspondence::new(vec![], size - 1, size - 1).unwrap(),
        );
        let mut edits = Edits::new();
        edits.remove_atom(AtomHandle::Id(AtomId((size - 1) as u32)));
        let added = edits.add_atom(AtomForm::from_element(Element::O));
        edits.add_bond(
            AtomHandle::Id(AtomId((size - 2) as u32)),
            added,
            BondForm::from_order(1),
        );
        molecule
            .apply(edits.clone())
            .expect("benchmark edit batch succeeds");
        molecule
            .meet_pushout(&molecule, &overlap)
            .expect("benchmark overlap is admissible");

        group.bench_function(BenchmarkId::new("editor_session/path", size), |b| {
            b.iter(|| {
                let mut editor = molecule.edit();
                editor.remove(&[AtomId((size - 1) as u32)], &[]);
                let added = editor.add_atom(AtomForm::from_element(Element::O));
                editor.add_bond(AtomId((size - 2) as u32), added, BondForm::from_order(1));
                black_box(editor.tracked_build())
            })
        });
        group.bench_function(BenchmarkId::new("remove/path", size), |b| {
            b.iter_batched(
                || molecule.edit(),
                |mut editor| {
                    editor.remove(black_box(&removed), &[]);
                    black_box(editor)
                },
                BatchSize::SmallInput,
            )
        });
        group.bench_function(BenchmarkId::new("tracked_remove/path", size), |b| {
            b.iter_batched(
                || molecule.edit(),
                |mut editor| {
                    let compaction = editor.tracked_remove(black_box(&removed), &[]);
                    black_box((editor, compaction))
                },
                BatchSize::SmallInput,
            )
        });
        group.bench_function(BenchmarkId::new("apply/path_three_edits", size), |b| {
            b.iter_batched(
                || edits.clone(),
                |edits| {
                    black_box(&molecule)
                        .apply(edits)
                        .expect("benchmark edit batch succeeds")
                },
                BatchSize::SmallInput,
            )
        });

        group.bench_function(
            BenchmarkId::new("tracked_apply/path_three_edits", size),
            |b| {
                b.iter_batched(
                    || edits.clone(),
                    |edits| {
                        black_box(&molecule)
                            .tracked_apply(edits)
                            .expect("benchmark edit batch succeeds")
                    },
                    BatchSize::SmallInput,
                )
            },
        );
        group.bench_function(
            BenchmarkId::new("tracked_transact/path_three_edits", size),
            |b| {
                b.iter_batched(
                    || (molecule.edit(), edits.clone()),
                    |(mut editor, edits)| {
                        let result = editor
                            .tracked_transact(edits)
                            .expect("benchmark edit batch succeeds");
                        black_box((editor, result))
                    },
                    BatchSize::SmallInput,
                )
            },
        );
        group.bench_function(BenchmarkId::new("combine/path_pair", size), |b| {
            b.iter(|| black_box(&molecule).combine(black_box(&molecule)))
        });
        group.bench_function(BenchmarkId::new("split/two_paths", size), |b| {
            b.iter(|| black_box(&disconnected).split())
        });
        group.bench_function(BenchmarkId::new("tracked_split/two_paths", size), |b| {
            b.iter(|| black_box(&disconnected).tracked_split())
        });
        group.bench_function(BenchmarkId::new("meet_pushout/path_pair", size), |b| {
            b.iter(|| black_box(&molecule).meet_pushout(black_box(&molecule), black_box(&overlap)))
        });
        group.bench_function(
            BenchmarkId::new("tracked_meet_pushout/path_pair", size),
            |b| {
                b.iter(|| {
                    black_box(&molecule)
                        .tracked_meet_pushout(black_box(&molecule), black_box(&overlap))
                })
            },
        );
    }
    group.finish();
}

criterion_group!(benches, benchmark_reaction, benchmark_mutation);
criterion_main!(benches);
