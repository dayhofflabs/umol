//! Subgraph-isomorphism benchmarks across the four retained algorithms
//! (`Vf2Rdkit`, `RayKirsch`, `BitsetRi`, `Glasgow`).
//!
//! Covers four groups:
//!
//! 1. `subgraph_isomorphism`: labeled enumeration plus `first` (`has`) and
//!    `count` modes on the core molecular/polyhedral target × query matrix.
//! 2. `subgraph_isomorphism/unsat_hard`: unsatisfiable and dense combinatorial
//!    existence checks (girth-5 fullerenes/polyhedra, bipartite grids, and
//!    $K_6 \hookrightarrow K_{10} - M_5$).
//! 3. `subgraph_isomorphism/large`: prepared `U128` (`n = 70`) and `Wide`
//!    (`n = 144`, `256`) queries reusing `PreparedSubgraphQuery`,
//!    `PreparedSubgraphTarget`, and `SubgraphIsomorphismScratch`.
//! 4. `subgraph_isomorphism/symmetry_reduced`: push-down query, host, and
//!    combined symmetry reduction (`none`, `pattern`, `host`, `both`) on
//!    symmetric polyhedra and fullerenes.

#[allow(dead_code)]
mod fixtures;

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion};
use umol_graph_core::SubgraphIsomorphismAlgorithm::{BitsetRi, Glasgow, RayKirsch, Vf2Rdkit};
use umol_graph_core::{
    AutomorphismAlgorithm, EdgeId, Graph, NodeId, PreparedSubgraphQuery, PreparedSubgraphTarget,
    QuerySymmetryConstraints, SubgraphIsomorphismAlgorithm, SubgraphIsomorphismPredicate,
    SubgraphIsomorphismScratch, TargetSymmetry,
};

const SUBISO: [SubgraphIsomorphismAlgorithm; 4] = [Vf2Rdkit, RayKirsch, BitsetRi, Glasgow];

const VISIT_MODE_PAIRS: [(&str, &str); 3] = [
    ("naphthalene", "path_3"),
    ("cubane", "square"),
    ("c60", "hexagon"),
];

struct UnlabeledPredicate;

impl SubgraphIsomorphismPredicate for UnlabeledPredicate {
    fn match_node(&mut self, _query_node: NodeId, _target_node: NodeId) -> bool {
        true
    }

    fn match_edge(&mut self, _query_edge: EdgeId, _target_edge: EdgeId) -> bool {
        true
    }

    fn all_edges_undetermined(&self) -> bool {
        true
    }

    fn query_node_class_count(&self) -> Option<usize> {
        Some(1)
    }

    fn query_node_class(&self, _query_node: NodeId) -> usize {
        0
    }

    fn node_candidates_u64(&mut self, _query_node: NodeId, eligible: u64) -> u64 {
        eligible
    }

    fn node_candidates_u128(&mut self, _query_node: NodeId, eligible: u128) -> u128 {
        eligible
    }

    fn fill_node_candidates_words(
        &mut self,
        _query_node: NodeId,
        eligible: &[u64],
        out: &mut [u64],
    ) {
        let len = eligible.len().min(out.len());
        out[..len].copy_from_slice(&eligible[..len]);
        out[len..].fill(0);
    }
}

fn subiso_name(algorithm: SubgraphIsomorphismAlgorithm) -> &'static str {
    match algorithm {
        Vf2Rdkit => "vf2rdkit",
        RayKirsch => "raykirsch",
        BitsetRi => "bitset_ri",
        Glasgow => "glasgow",
    }
}

/// Complete graph `K_{2k}` minus a perfect matching (`K_{2,2,...,2}`, clique
/// number `k`, every vertex of degree `2k - 2`).
fn cocktail_party(half_n: usize) -> Graph {
    let n = 2 * half_n;
    let mut edges = Vec::new();
    for i in 0..n as u32 {
        for j in (i + 1)..n as u32 {
            if j != i + half_n as u32 {
                edges.push([i, j]);
            }
        }
    }
    Graph::new(n, &edges)
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
                        black_box(target.enumerate_subgraph_isomorphisms(
                            query,
                            &mut |_: NodeId, _: NodeId| true,
                            &mut |qe: EdgeId, he: EdgeId| {
                                q_labels[qe.index()] == t_labels[he.index()]
                            },
                            algorithm,
                        ))
                    });
                });
                if include_visit_modes {
                    group.bench_function(
                        format!("{tname}/{qname}/{}/first", subiso_name(algorithm)),
                        |b| {
                            b.iter(|| {
                                black_box(target.has_subgraph_isomorphism(
                                    query,
                                    &mut |_: NodeId, _: NodeId| true,
                                    &mut |qe: EdgeId, he: EdgeId| {
                                        q_labels[qe.index()] == t_labels[he.index()]
                                    },
                                    algorithm,
                                ))
                            });
                        },
                    );
                    group.bench_function(
                        format!("{tname}/{qname}/{}/count", subiso_name(algorithm)),
                        |b| {
                            b.iter(|| {
                                black_box(target.count_subgraph_isomorphisms(
                                    query,
                                    &mut |_: NodeId, _: NodeId| true,
                                    &mut |qe: EdgeId, he: EdgeId| {
                                        q_labels[qe.index()] == t_labels[he.index()]
                                    },
                                    algorithm,
                                ))
                            });
                        },
                    );
                }
            }
        }
    }
    group.finish();
}

fn subgraph_isomorphism_unsat_hard(c: &mut Criterion) {
    let cases: Vec<(&str, Graph, Graph)> = vec![
        ("petersen/square", fixtures::petersen(), fixtures::cycle(4)),
        (
            "dodecahedron/square",
            fixtures::dodecahedron(),
            fixtures::cycle(4),
        ),
        ("c60/square", fixtures::fullerene_c60(), fixtures::cycle(4)),
        ("c60/k4", fixtures::fullerene_c60(), fixtures::complete(4)),
        ("c70/square", fixtures::fullerene_c70(), fixtures::cycle(4)),
        (
            "grid_8x8/triangle",
            fixtures::grid(8, 8),
            fixtures::cycle(3),
        ),
        ("k10_minus_m5/k6", cocktail_party(5), fixtures::complete(6)),
    ];

    let mut group = c.benchmark_group("subgraph_isomorphism/unsat_hard");
    for (name, target, query) in &cases {
        let prepared_target = PreparedSubgraphTarget::new(target);
        for algorithm in SUBISO {
            let prepared_query = PreparedSubgraphQuery::new(query, algorithm);
            let mut scratch = SubgraphIsomorphismScratch::new();
            let mut pred = UnlabeledPredicate;
            assert!(!prepared_query.has_match(&prepared_target, None, &mut pred, &mut scratch));
            group.bench_function(format!("{name}/{}", subiso_name(algorithm)), |b| {
                b.iter(|| {
                    black_box(prepared_query.has_match(
                        &prepared_target,
                        None,
                        &mut pred,
                        &mut scratch,
                    ))
                });
            });
        }
    }
    group.finish();
}

fn subgraph_isomorphism_large(c: &mut Criterion) {
    let cases: Vec<(&str, Graph, Graph)> = vec![
        ("c70/hexagon", fixtures::fullerene_c70(), fixtures::cycle(6)),
        (
            "grid_12x12/square",
            fixtures::grid(12, 12),
            fixtures::cycle(4),
        ),
        (
            "grid_12x12/grid_3x3",
            fixtures::grid(12, 12),
            fixtures::grid(3, 3),
        ),
        (
            "grid_16x16/grid_3x3",
            fixtures::grid(16, 16),
            fixtures::grid(3, 3),
        ),
    ];

    let mut group = c.benchmark_group("subgraph_isomorphism/large");
    for (name, target, query) in &cases {
        let prepared_target = PreparedSubgraphTarget::new(target);
        for algorithm in SUBISO {
            let prepared_query = PreparedSubgraphQuery::new(query, algorithm);
            let mut scratch = SubgraphIsomorphismScratch::new();
            let mut pred = UnlabeledPredicate;
            group.bench_function(format!("{name}/{}/has", subiso_name(algorithm)), |b| {
                b.iter(|| {
                    black_box(prepared_query.has_match(
                        &prepared_target,
                        None,
                        &mut pred,
                        &mut scratch,
                    ))
                });
            });
            group.bench_function(format!("{name}/{}/count", subiso_name(algorithm)), |b| {
                b.iter(|| {
                    black_box(prepared_query.count_matches(
                        &prepared_target,
                        None,
                        &mut pred,
                        &mut scratch,
                    ))
                });
            });
        }
    }
    group.finish();
}

fn subgraph_isomorphism_symmetry_reduced(c: &mut Criterion) {
    let cases: Vec<(&str, Graph, Graph)> = vec![
        ("cubane/square", fixtures::cubane(), fixtures::cycle(4)),
        (
            "dodecahedron/pentagon",
            fixtures::dodecahedron(),
            fixtures::cycle(5),
        ),
        ("c60/hexagon", fixtures::fullerene_c60(), fixtures::cycle(6)),
        ("c70/hexagon", fixtures::fullerene_c70(), fixtures::cycle(6)),
    ];

    let mut group = c.benchmark_group("subgraph_isomorphism/symmetry_reduced");
    for (name, target, query) in &cases {
        let prepared_target = PreparedSubgraphTarget::new(target);
        let q_sym = QuerySymmetryConstraints::from_graph(
            query,
            |_| 0u8,
            |_| 0u8,
            AutomorphismAlgorithm::TrailIr,
        );
        let t_sym =
            TargetSymmetry::from_graph(target, |_| 0u8, |_| 0u8, AutomorphismAlgorithm::TrailIr);

        for (mode_name, use_q, use_h) in [
            ("none", false, false),
            ("pattern", true, false),
            ("host", false, true),
            ("both", true, true),
        ] {
            let target_sym = use_h.then_some(&t_sym);
            for algorithm in SUBISO {
                let mut prepared_query = PreparedSubgraphQuery::new(query, algorithm);
                if use_q {
                    prepared_query = prepared_query.with_symmetry(q_sym.clone());
                }
                let mut scratch = SubgraphIsomorphismScratch::new();
                let mut pred = UnlabeledPredicate;
                group.bench_function(
                    format!("{name}/{mode_name}/{}", subiso_name(algorithm)),
                    |b| {
                        b.iter(|| {
                            black_box(prepared_query.count_matches(
                                &prepared_target,
                                target_sym,
                                &mut pred,
                                &mut scratch,
                            ))
                        });
                    },
                );
            }
        }
    }
    group.finish();
}

criterion_group!(
    benches,
    subgraph_isomorphism,
    subgraph_isomorphism_unsat_hard,
    subgraph_isomorphism_large,
    subgraph_isomorphism_symmetry_reduced
);
criterion_main!(benches);
