//! Reaction matching and matched-application benchmarks.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion};
use umol_chem::element::Element;
use umol_graph_core::{
    Correspondence, RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm,
};
use umol_graph_ir::ir::{
    AromaticSystemDelta, AromaticSystemFieldChange, AromaticSystemForm, AromaticSystemId,
    AtomConstraintForm, AtomForm, AtomId, BondForm, Constraint, ConstraintSpan, Delta, Deltas,
    ElectronCountsForm, EntitySpan, Molecule, MoleculeCorrespondence, MoleculeEntries, Reaction,
    ReactionSpan, ReactionSpanEntries, SubstructureMatchAlgorithm, SubstructureMatchConfig,
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
    group.bench_function("reverse", |b| {
        b.iter(|| black_box(black_box(&reversal).reverse()))
    });
    group.finish();
}

criterion_group!(benches, benchmark_reaction);
criterion_main!(benches);
