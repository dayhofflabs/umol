//! Reaction matching and matched-application benchmarks.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion};
use umol_chem::element::Element;
use umol_graph_core::{
    Correspondence, RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm,
};
use umol_graph_ir::ir::{
    AromaticSystemDelta, AromaticSystemFieldChange, AromaticSystemForm, AromaticSystemId, AtomForm,
    AtomId, BondForm, Delta, Deltas, ElectronCountsForm, Molecule, MoleculeCorrespondence,
    MoleculeEntries, Reaction, SubstructureMatchAlgorithm, SubstructureMatchConfig,
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

fn benchmark_reaction(c: &mut Criterion) {
    let (reaction, host, correspondence) = application_case();
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
    group.finish();
}

criterion_group!(benches, benchmark_reaction);
criterion_main!(benches);
