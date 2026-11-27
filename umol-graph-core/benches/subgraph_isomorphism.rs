//! Subgraph isomorphism benchmarks across all supported algorithms.
//!
//! Full enumeration runs across the complete target × query × algorithm matrix; `first`
//! (time to first match) and `count` (enumeration without collecting) run on a representative
//! subset of `(target, query)` pairs.

#[allow(dead_code)]
mod fixtures;

use std::ops::ControlFlow;

use criterion::{criterion_group, criterion_main, Criterion};
use umol_graph_core::SubgraphIsomorphismAlgorithm::{
    ArcMatch, BitsetRi, RayKirsch, Ri, Ullmann, Vf2, Vf2Rdkit,
};
use umol_graph_core::{
    EdgeId, Graph, NodeId, SubgraphIsomorphismAlgorithm, ARCMATCH_DEFAULT_PATH_LENGTH,
};

const SUBISO: [SubgraphIsomorphismAlgorithm; 6] = [
    Vf2,
    Ullmann,
    Ri,
    ArcMatch {
        path_length: ARCMATCH_DEFAULT_PATH_LENGTH,
    },
    Vf2Rdkit,
    RayKirsch,
];

const VISIT_MODE_PAIRS: [(&str, &str); 3] = [
    ("naphthalene", "path_3"),
    ("cubane", "square"),
    ("c60", "hexagon"),
];

fn subiso_name(algorithm: SubgraphIsomorphismAlgorithm) -> &'static str {
    match algorithm {
        Vf2 => "vf2",
        Ullmann => "ullmann",
        Ri => "ri",
        ArcMatch { .. } => "arcmatch",
        Vf2Rdkit => "vf2rdkit",
        RayKirsch => "raykirsch",
        BitsetRi => "bitset_ri",
    }
}

fn subgraph_isomorphism(c: &mut Criterion) {
    let targets: Vec<(&str, Graph)> = vec![
        ("hexagon", fixtures::cycle(6)),
        ("naphthalene", fixtures::naphthalene()),
        ("cubane", fixtures::cubane()),
        ("adamantane", fixtures::adamantane()),
        ("petersen", fixtures::petersen()),
        ("dodecahedron", fixtures::dodecahedron()),
        ("c60", fixtures::fullerene_c60()),
    ];

    let queries: Vec<(&str, Graph)> = vec![
        ("edge", Graph::new(2, &[[0, 1]])),
        ("path_3", fixtures::path(3)),
        ("triangle", fixtures::cycle(3)),
        ("square", fixtures::cycle(4)),
        ("hexagon", fixtures::cycle(6)),
    ];

    let mut group = c.benchmark_group("subgraph_isomorphism");
    for (tname, target) in &targets {
        let t_labels = fixtures::alternating_edge_labels(target);
        for (qname, query) in &queries {
            if query.node_count() > target.node_count() {
                continue;
            }
            let q_labels = fixtures::alternating_edge_labels(query);
            let include_visit_modes = VISIT_MODE_PAIRS.contains(&(*tname, *qname));
            for algorithm in SUBISO {
                group.bench_function(format!("{tname}/{qname}/{}", subiso_name(algorithm)), |b| {
                    b.iter(|| {
                        target.enumerate_subgraph_isomorphisms(
                            query,
                            &mut |_: NodeId, _: NodeId| true,
                            &mut |qe: EdgeId, he: EdgeId| {
                                q_labels[qe.index()] == t_labels[he.index()]
                            },
                            algorithm,
                        )
                    });
                });
                if include_visit_modes {
                    group.bench_function(
                        format!("{tname}/{qname}/{}/first", subiso_name(algorithm)),
                        |b| {
                            b.iter(|| {
                                target.visit_subgraph_isomorphisms(
                                    query,
                                    &mut |_: NodeId, _: NodeId| true,
                                    &mut |qe: EdgeId, he: EdgeId| {
                                        q_labels[qe.index()] == t_labels[he.index()]
                                    },
                                    algorithm,
                                    |_| ControlFlow::Break(()),
                                )
                            });
                        },
                    );
                    group.bench_function(
                        format!("{tname}/{qname}/{}/count", subiso_name(algorithm)),
                        |b| {
                            b.iter(|| {
                                let mut count = 0usize;
                                let _: ControlFlow<()> = target.visit_subgraph_isomorphisms(
                                    query,
                                    &mut |_: NodeId, _: NodeId| true,
                                    &mut |qe: EdgeId, he: EdgeId| {
                                        q_labels[qe.index()] == t_labels[he.index()]
                                    },
                                    algorithm,
                                    |_| {
                                        count += 1;
                                        ControlFlow::Continue(())
                                    },
                                );
                                count
                            });
                        },
                    );
                }
            }
        }
    }
    group.finish();
}

criterion_group!(benches, subgraph_isomorphism);
criterion_main!(benches);
