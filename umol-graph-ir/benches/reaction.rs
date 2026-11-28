//! Molecule mutation, reaction matching, and matched-application benchmarks.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use umol_chem::element::Element;
use umol_graph_core::{
    AutomorphismAlgorithm, Correspondence, GraphCorrespondence, NodeId,
    RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm,
};
use umol_graph_ir::ir::{
    AromaticSystemDelta, AromaticSystemFieldChange, AromaticSystemForm, AromaticSystemId,
    AtomConstraintForm, AtomDelta, AtomFieldChange, AtomForm, AtomHandle, AtomId, BondDelta,
    BondForm, BondId, Constraint, ConstraintSpan, Delta, Deltas, Edits, ElectronCountsForm,
    EntitySpan, Molecule, MoleculeCorrespondence, MoleculeEntries, NumForm,
    PreparedSubstructureHost, Reaction, ReactionSpan, ReactionSpanEntries, RingScope,
    RollbackMethod, SubstructureMatchAlgorithm, SubstructureMatchConfig,
};
use umol_graph_ir::mol_dsl;

const MATCH_CONFIG: SubstructureMatchConfig = SubstructureMatchConfig {
    match_algorithm: SubstructureMatchAlgorithm::GraphAndOverlays,
    subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::BitsetRi,
    relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
    automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
    reduce_pattern_symmetry: false,
    reduce_host_symmetry: false,
};

const MATCH_CONFIG_VF2: SubstructureMatchConfig = SubstructureMatchConfig {
    subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::Vf2Rdkit,
    ..MATCH_CONFIG
};

const MATCH_CONFIG_SYM: SubstructureMatchConfig = SubstructureMatchConfig {
    reduce_pattern_symmetry: true,
    reduce_host_symmetry: true,
    ..MATCH_CONFIG
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
    let prepared_host = PreparedSubstructureHost::new(&host);
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
    group.bench_function("apply_to_host/all_matches", |b| {
        b.iter(|| {
            black_box(
                reaction
                    .apply_to_host(black_box(&prepared_host), MATCH_CONFIG)
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
        edits.remove_topology_cascading(vec![AtomHandle::Id(AtomId((size - 1) as u32))], vec![]);
        let added = edits.add_atom(AtomForm::from_element(Element::O));
        edits.add_bond(
            AtomHandle::Id(AtomId((size - 2) as u32)),
            added,
            BondForm::from_order(1),
        );
        molecule
            .clone()
            .apply(edits.clone())
            .expect("benchmark edit batch succeeds");
        molecule
            .meet_pushout(&molecule, &overlap)
            .expect("benchmark overlap is admissible");

        group.bench_function(BenchmarkId::new("editor_session/path", size), |b| {
            b.iter(|| {
                let mut editor = molecule.clone().edit();
                editor.remove_topology_cascading(&[AtomId((size - 1) as u32)], &[]);
                let added = editor.add_atom(AtomForm::from_element(Element::O));
                editor.add_bond(AtomId((size - 2) as u32), added, BondForm::from_order(1));
                black_box(editor.finish().unwrap())
            })
        });
        group.bench_function(BenchmarkId::new("remove_topology/path", size), |b| {
            b.iter_batched(
                || molecule.clone().edit(),
                |mut editor| {
                    editor.remove_topology_cascading(black_box(&removed), &[]);
                    black_box(editor)
                },
                BatchSize::SmallInput,
            )
        });
        group.bench_function(BenchmarkId::new("apply/path_three_edits", size), |b| {
            b.iter_batched(
                || (molecule.clone(), edits.clone()),
                |(molecule, edits)| {
                    black_box(molecule)
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
                    || (molecule.clone(), edits.clone()),
                    |(molecule, edits)| {
                        black_box(molecule)
                            .tracked_apply(edits)
                            .expect("benchmark edit batch succeeds")
                    },
                    BatchSize::SmallInput,
                )
            },
        );
        for (operation, rollback_method) in [
            ("tracked_transact/path_three_edits", RollbackMethod::Journal),
            (
                "tracked_transact_clone/path_three_edits",
                RollbackMethod::Clone,
            ),
        ] {
            group.bench_function(BenchmarkId::new(operation, size), |b| {
                b.iter_batched(
                    || (molecule.clone(), edits.clone()),
                    |(mut molecule, edits)| {
                        let result = molecule
                            .tracked_transact([edits], rollback_method)
                            .expect("benchmark edit batch succeeds");
                        black_box((molecule, result))
                    },
                    BatchSize::SmallInput,
                )
            });
        }
        group.bench_function(BenchmarkId::new("combine/path_pair", size), |b| {
            b.iter_batched(
                || molecule.clone(),
                |combined| black_box(combined.combine(black_box(&molecule))),
                BatchSize::SmallInput,
            )
        });
        group.bench_function(BenchmarkId::new("combine_from/path_pair", size), |b| {
            b.iter_batched(
                || molecule.clone(),
                |mut combined| {
                    combined.combine_from(black_box(&molecule));
                    black_box(combined)
                },
                BatchSize::SmallInput,
            )
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

fn library_hosts() -> Vec<Molecule> {
    let biphenyl = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 12],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(1)),
            (AtomId(2), AtomId(3), BondForm::from_order(1)),
            (AtomId(3), AtomId(4), BondForm::from_order(1)),
            (AtomId(4), AtomId(5), BondForm::from_order(1)),
            (AtomId(5), AtomId(0), BondForm::from_order(1)),
            (AtomId(6), AtomId(7), BondForm::from_order(1)),
            (AtomId(7), AtomId(8), BondForm::from_order(1)),
            (AtomId(8), AtomId(9), BondForm::from_order(1)),
            (AtomId(9), AtomId(10), BondForm::from_order(1)),
            (AtomId(10), AtomId(11), BondForm::from_order(1)),
            (AtomId(11), AtomId(6), BondForm::from_order(1)),
            (AtomId(0), AtomId(6), BondForm::from_order(1)),
        ],
        aromatic: vec![
            ((0..6).map(AtomId).collect(), AromaticSystemForm::default()),
            ((6..12).map(AtomId).collect(), AromaticSystemForm::default()),
        ],
        ..Default::default()
    });
    let naphthalene_diol = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::O),
            AtomForm::from_element(Element::O),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(1)),
            (AtomId(2), AtomId(3), BondForm::from_order(1)),
            (AtomId(3), AtomId(4), BondForm::from_order(1)),
            (AtomId(4), AtomId(5), BondForm::from_order(1)),
            (AtomId(5), AtomId(0), BondForm::from_order(1)),
            (AtomId(4), AtomId(6), BondForm::from_order(1)),
            (AtomId(6), AtomId(7), BondForm::from_order(1)),
            (AtomId(7), AtomId(8), BondForm::from_order(1)),
            (AtomId(8), AtomId(9), BondForm::from_order(1)),
            (AtomId(9), AtomId(5), BondForm::from_order(1)),
            (AtomId(1), AtomId(10), BondForm::from_order(1)),
            (AtomId(7), AtomId(11), BondForm::from_order(1)),
        ],
        aromatic: vec![((0..10).map(AtomId).collect(), AromaticSystemForm::default())],
        ..Default::default()
    });
    let adipic_acid = mol_dsl!(
        r#"{:atoms ["O" "C" "O" "C" "C" "C" "C" "C" "O" "O"]
            :bonds [[0 1 "1"] [1 2 "2"] [1 3 "1"] [3 4 "1"] [4 5 "1"]
                    [5 6 "1"] [6 7 "1"] [7 8 "2"] [7 9 "1"]]}"#
    );
    let benzamide = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::O),
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::Cl),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(1)),
            (AtomId(2), AtomId(3), BondForm::from_order(1)),
            (AtomId(3), AtomId(4), BondForm::from_order(1)),
            (AtomId(4), AtomId(5), BondForm::from_order(1)),
            (AtomId(5), AtomId(0), BondForm::from_order(1)),
            (AtomId(0), AtomId(6), BondForm::from_order(1)),
            (AtomId(6), AtomId(7), BondForm::from_order(2)),
            (AtomId(6), AtomId(8), BondForm::from_order(1)),
            (AtomId(3), AtomId(9), BondForm::from_order(1)),
        ],
        aromatic: vec![((0..6).map(AtomId).collect(), AromaticSystemForm::default())],
        ..Default::default()
    });
    let cyclohexene_diol = mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "C" "C" "O" "O" "C" "N"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [2 6 "1"] [5 7 "1"] [3 8 "1"] [8 9 "1"]]}"#
    );
    let piperidine_ester = mol_dsl!(
        r#"{:atoms ["N" "C" "C" "C" "C" "C" "C" "O" "O" "C" "C" "Cl"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [3 6 "1"] [6 7 "2"] [6 8 "1"] [8 9 "1"] [9 10 "1"] [10 11 "1"]]}"#
    );
    vec![
        biphenyl,
        naphthalene_diol,
        adipic_acid,
        benzamide,
        cyclohexene_diol,
        piperidine_ester,
    ]
}

fn reaction_library() -> Vec<Reaction> {
    let charge_atom0 = |new_charge: i64| {
        Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(0),
            change: AtomFieldChange::Charge {
                old: NumForm::Undetermined,
                new: NumForm::Lit(new_charge),
            },
        })])
    };
    let attach_cl_at = |site: u32, next_atom: u32, next_bond: u32| {
        Deltas::from_iter([
            Delta::Atom(AtomDelta::Add {
                id: AtomId(next_atom),
                attributes: AtomForm::from_element(Element::Cl),
            }),
            Delta::Bond(BondDelta::Add {
                id: BondId(next_bond),
                atoms: [AtomId(site), AtomId(next_atom)],
                attributes: BondForm::from_order(1),
            }),
        ])
    };

    vec![
        // Acyclic functional-group transformations.
        Reaction::new(
            mol_dsl!(r#"{:atoms ["O" "C" "O"] :bonds [[0 1 "1"] [1 2 "2"]]}"#),
            charge_atom0(-1),
        ),
        Reaction::new(
            mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "2"]]}"#),
            Deltas::from_iter([Delta::Bond(BondDelta::Remove {
                id: BondId(0),
                atoms: [AtomId(0), AtomId(1)],
                attributes: BondForm::from_order(2),
            })]),
        ),
        Reaction::new(
            mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#),
            attach_cl_at(1, 2, 1),
        ),
        Reaction::new(
            mol_dsl!(r#"{:atoms ["N" "C"] :bonds [[0 1 "1"]]}"#),
            charge_atom0(1),
        ),
        Reaction::new(
            mol_dsl!(r#"{:atoms ["C" "Cl"] :bonds [[0 1 "1"]]}"#),
            Deltas::from_iter([
                Delta::Bond(BondDelta::Remove {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    attributes: BondForm::from_order(1),
                }),
                Delta::Atom(AtomDelta::Remove {
                    id: AtomId(1),
                    attributes: AtomForm::from_element(Element::Cl),
                }),
            ]),
        ),
        // Ring-constrained transformations exercising lazy `RingSet` caching on `PreparedSubstructureHost`.
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C)
                    .with_constraint(AtomConstraintForm::ring_membership(RingScope::All, 1))
                    .with_constraint(AtomConstraintForm::degree(2))],
                ..Default::default()
            }),
            attach_cl_at(0, 1, 0),
        ),
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::N)
                    .with_constraint(AtomConstraintForm::ring_membership(RingScope::All, 1))],
                ..Default::default()
            }),
            charge_atom0(1),
        ),
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C)
                        .with_constraint(AtomConstraintForm::ring_membership(RingScope::All, 1)),
                    AtomForm::from_element(Element::O),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
                ..Default::default()
            }),
            charge_atom0(1),
        ),
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C)
                        .with_constraint(AtomConstraintForm::ring_membership(RingScope::All, 1)),
                    AtomForm::from_element(Element::C)
                        .with_constraint(AtomConstraintForm::ring_membership(RingScope::All, 1)),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))],
                ..Default::default()
            }),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(2),
                    attributes: AtomForm::from_element(Element::O),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(1),
                    atoms: [AtomId(0), AtomId(2)],
                    attributes: BondForm::from_order(1),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(2),
                    atoms: [AtomId(1), AtomId(2)],
                    attributes: BondForm::from_order(1),
                }),
            ]),
        ),
        // Symmetric-LHS transformations.
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); 6],
                bonds: (0..6)
                    .map(|i| (AtomId(i), AtomId((i + 1) % 6), BondForm::from_order(1)))
                    .collect(),
                aromatic: vec![((0..6).map(AtomId).collect(), AromaticSystemForm::default())],
                ..Default::default()
            }),
            attach_cl_at(0, 6, 6),
        ),
        Reaction::new(
            mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(2),
                    attributes: AtomForm::from_element(Element::O),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(1),
                    atoms: [AtomId(0), AtomId(2)],
                    attributes: BondForm::from_order(1),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(2),
                    atoms: [AtomId(1), AtomId(2)],
                    attributes: BondForm::from_order(1),
                }),
            ]),
        ),
        Reaction::new(
            mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
            charge_atom0(1),
        ),
        // Rare-heteroatom rules exercising fast element-bitmask pre-screening.
        Reaction::new(
            mol_dsl!(r#"{:atoms ["S" "O"] :bonds [[0 1 "2"]]}"#),
            charge_atom0(1),
        ),
        Reaction::new(
            mol_dsl!(r#"{:atoms ["P" "O"] :bonds [[0 1 "2"]]}"#),
            charge_atom0(1),
        ),
        Reaction::new(
            mol_dsl!(r#"{:atoms ["B" "O"] :bonds [[0 1 "1"]]}"#),
            charge_atom0(-1),
        ),
        Reaction::new(
            mol_dsl!(r#"{:atoms ["Si" "O"] :bonds [[0 1 "1"]]}"#),
            charge_atom0(1),
        ),
    ]
}

fn benchmark_reaction_library(c: &mut Criterion) {
    let hosts = library_hosts();
    let prepared_hosts: Vec<PreparedSubstructureHost> =
        hosts.iter().map(PreparedSubstructureHost::new).collect();
    let reactions = reaction_library();

    for host in &hosts {
        for reaction in &reactions {
            let direct = reaction
                .apply(host, MATCH_CONFIG)
                .expect("benchmark reaction satisfies preconditions")
                .collect::<Result<Vec<_>, _>>()
                .expect("benchmark application succeeds");
            let prepared = reaction
                .apply_to_host(&PreparedSubstructureHost::new(host), MATCH_CONFIG)
                .expect("benchmark reaction satisfies preconditions")
                .collect::<Result<Vec<_>, _>>()
                .expect("benchmark application succeeds");
            assert_eq!(direct.len(), prepared.len());
        }
    }

    let mut group = c.benchmark_group("reaction_library");
    group.bench_function("apply/vf2rdkit_unsym", |b| {
        b.iter(|| {
            let mut total = 0usize;
            for host in &hosts {
                for reaction in &reactions {
                    total += reaction.apply(host, MATCH_CONFIG_VF2).unwrap().count();
                }
            }
            black_box(total)
        })
    });
    group.bench_function("apply/bitset_ri_unsym", |b| {
        b.iter(|| {
            let mut total = 0usize;
            for host in &hosts {
                for reaction in &reactions {
                    total += reaction.apply(host, MATCH_CONFIG).unwrap().count();
                }
            }
            black_box(total)
        })
    });
    group.bench_function("apply_to_host_per_run/bitset_ri_unsym", |b| {
        b.iter(|| {
            let mut total = 0usize;
            for host in &hosts {
                let prepared = PreparedSubstructureHost::new(host);
                for reaction in &reactions {
                    total += reaction
                        .apply_to_host(&prepared, MATCH_CONFIG)
                        .unwrap()
                        .count();
                }
            }
            black_box(total)
        })
    });
    group.bench_function("apply_to_host_prebuilt/bitset_ri_unsym", |b| {
        b.iter(|| {
            let mut total = 0usize;
            for host in &prepared_hosts {
                for reaction in &reactions {
                    total += reaction.apply_to_host(host, MATCH_CONFIG).unwrap().count();
                }
            }
            black_box(total)
        })
    });
    group.bench_function("apply/bitset_ri_sym", |b| {
        b.iter(|| {
            let mut total = 0usize;
            for host in &hosts {
                for reaction in &reactions {
                    total += reaction.apply(host, MATCH_CONFIG_SYM).unwrap().count();
                }
            }
            black_box(total)
        })
    });
    group.bench_function("apply_to_host_per_run/bitset_ri_sym", |b| {
        b.iter(|| {
            let mut total = 0usize;
            for host in &hosts {
                let prepared = PreparedSubstructureHost::new(host);
                for reaction in &reactions {
                    total += reaction
                        .apply_to_host(&prepared, MATCH_CONFIG_SYM)
                        .unwrap()
                        .count();
                }
            }
            black_box(total)
        })
    });
    group.bench_function("apply_to_host_prebuilt/bitset_ri_sym", |b| {
        b.iter(|| {
            let mut total = 0usize;
            for host in &prepared_hosts {
                for reaction in &reactions {
                    total += reaction
                        .apply_to_host(host, MATCH_CONFIG_SYM)
                        .unwrap()
                        .count();
                }
            }
            black_box(total)
        })
    });
    group.finish();
}

criterion_group!(
    benches,
    benchmark_reaction,
    benchmark_reaction_library,
    benchmark_mutation
);
criterion_main!(benches);
