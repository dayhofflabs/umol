//! Maximum-matching, perfect/maximum matching enumeration, and planar FKT counting benchmarks.

#[allow(dead_code)]
mod fixtures;

use std::env;
use std::hint::black_box;
use std::ops::ControlFlow;
use std::time::{Duration, Instant};

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use umol_graph_core::{
    BipartiteMaximumMatchingAlgorithm, EdgeId, FaceBoundary, GeneralMaximumMatchingAlgorithm,
    Graph, Matching, MatchingEnumerationAlgorithm, NodeId, PlanarEmbedding,
};

const PREFIX_LIMITS: [usize; 3] = [1, 10, 100];
const DELAY_DIAGNOSTICS_ENV: &str = "UMOL_MATCHING_DELAY_DIAGNOSTICS";

fn maximum_matching(c: &mut Criterion) {
    let graphs = [
        ("hexagon", fixtures::cycle(6)),
        ("cubane", fixtures::cubane()),
        ("petersen", fixtures::petersen()),
        ("dodecahedron", fixtures::dodecahedron()),
        ("c60", fixtures::fullerene_c60()),
    ];

    let mut group = c.benchmark_group("general_maximum_matching/edmonds");
    for (name, g) in &graphs {
        let node_order: Vec<_> = g.node_ids().collect();
        group.bench_function(*name, |b| {
            b.iter(|| {
                g.general_maximum_matching(&node_order, GeneralMaximumMatchingAlgorithm::Edmonds)
            });
        });
    }
    group.finish();

    let bipartite_graphs = [
        ("hexagon", fixtures::cycle(6)),
        ("cubane", fixtures::cubane()),
        ("grid_16x16", fixtures::grid(16, 16)),
    ];

    let mut group = c.benchmark_group("bipartite_maximum_matching/hopcroft_karp");
    for (name, graph) in &bipartite_graphs {
        let node_order: Vec<_> = graph.node_ids().collect();
        group.bench_function(*name, |b| {
            b.iter(|| {
                graph
                    .bipartite_maximum_matching(
                        &node_order,
                        BipartiteMaximumMatchingAlgorithm::HopcroftKarp,
                    )
                    .expect("Hopcroft-Karp benchmark graphs are bipartite")
            });
        });
    }
    group.finish();

    let mut group =
        c.benchmark_group("bipartite_maximum_matching_or_general/hopcroft_karp_edmonds");
    for (name, graph) in &graphs {
        let node_order: Vec<_> = graph.node_ids().collect();
        group.bench_function(*name, |b| {
            b.iter(|| {
                graph.bipartite_maximum_matching_or_general(
                    &node_order,
                    BipartiteMaximumMatchingAlgorithm::HopcroftKarp,
                    GeneralMaximumMatchingAlgorithm::Edmonds,
                )
            });
        });
    }
    group.finish();
}

#[derive(Clone, Copy)]
enum MatchingMode {
    Perfect,
    Maximum,
}

struct MatchingCase {
    name: &'static str,
    graph: Graph,
    embedding: Option<PlanarEmbedding>,
    mode: MatchingMode,
    output_count: usize,
}

impl MatchingCase {
    fn visit<B>(&self, visitor: impl FnMut(Matching) -> ControlFlow<B>) -> ControlFlow<B> {
        match self.mode {
            MatchingMode::Perfect => self
                .graph
                .visit_perfect_matchings(MatchingEnumerationAlgorithm::BranchAndBound, visitor),
            MatchingMode::Maximum => self
                .graph
                .visit_maximum_matchings(MatchingEnumerationAlgorithm::BranchAndBound, visitor),
        }
    }

    fn collect(&self) -> Vec<Matching> {
        match self.mode {
            MatchingMode::Perfect => self
                .graph
                .enumerate_perfect_matchings(MatchingEnumerationAlgorithm::BranchAndBound),
            MatchingMode::Maximum => self
                .graph
                .enumerate_maximum_matchings(MatchingEnumerationAlgorithm::BranchAndBound),
        }
    }
}

fn embedding(graph: &Graph, fixture: &fixtures::GraphFixture) -> PlanarEmbedding {
    let faces = fixture
        .faces
        .iter()
        .map(|nodes| {
            let nodes: Vec<_> = nodes.iter().map(|&node| NodeId(node)).collect();
            let edges = nodes
                .iter()
                .zip(nodes.iter().cycle().skip(1))
                .take(nodes.len())
                .map(|(&first, &second)| {
                    graph
                        .find_edge(first, second)
                        .expect("every fixture face side must be a graph edge")
                })
                .collect();
            FaceBoundary::new(nodes, edges)
        })
        .collect();
    PlanarEmbedding::new(graph, faces, fixture.faces.len() - 1)
        .expect("benchmark embedding must be valid")
}

fn fixture_case(
    name: &'static str,
    source: &str,
    mode: MatchingMode,
    planar: bool,
) -> MatchingCase {
    let parsed = fixtures::parse(source);
    let graph = parsed.graph();
    let embedding = planar.then(|| embedding(&graph, &parsed));
    let output_count = count_outputs(&graph, mode);

    MatchingCase {
        name,
        graph,
        embedding,
        mode,
        output_count,
    }
}

fn count_outputs(graph: &Graph, mode: MatchingMode) -> usize {
    let mut count = 0;
    match mode {
        MatchingMode::Perfect => {
            let _ =
                graph.visit_perfect_matchings(MatchingEnumerationAlgorithm::BranchAndBound, |_| {
                    count += 1;
                    ControlFlow::<()>::Continue(())
                });
        }
        MatchingMode::Maximum => {
            let _ =
                graph.visit_maximum_matchings(MatchingEnumerationAlgorithm::BranchAndBound, |_| {
                    count += 1;
                    ControlFlow::<()>::Continue(())
                });
        }
    }
    count
}

fn prescribed_hole_path() -> MatchingCase {
    let graph = Graph::new(4, &[[0, 1], [1, 2], [2, 3]]);
    let embedding = PlanarEmbedding::new(
        &graph,
        vec![FaceBoundary::new(
            vec![
                NodeId(0),
                NodeId(1),
                NodeId(2),
                NodeId(3),
                NodeId(2),
                NodeId(1),
            ],
            vec![
                EdgeId(0),
                EdgeId(1),
                EdgeId(2),
                EdgeId(2),
                EdgeId(1),
                EdgeId(0),
            ],
        )],
        0,
    )
    .expect("path embedding must be valid");
    let output_count = count_outputs(&graph, MatchingMode::Perfect);

    MatchingCase {
        name: "prescribed_hole_path",
        graph,
        embedding: Some(embedding),
        mode: MatchingMode::Perfect,
        output_count,
    }
}

fn complete_bipartite(sides: usize) -> MatchingCase {
    let edges: Vec<_> = (0..sides as u32)
        .flat_map(|left| (sides as u32..2 * sides as u32).map(move |right| [left, right]))
        .collect();
    let graph = Graph::new(2 * sides, &edges);
    let output_count = count_outputs(&graph, MatchingMode::Perfect);

    MatchingCase {
        name: "complete_bipartite_6x6",
        graph,
        embedding: None,
        mode: MatchingMode::Perfect,
        output_count,
    }
}

fn corpus() -> Vec<MatchingCase> {
    vec![
        fixture_case("benzene", fixtures::BENZENE, MatchingMode::Perfect, true),
        fixture_case(
            "naphthalene",
            fixtures::NAPHTHALENE,
            MatchingMode::Perfect,
            true,
        ),
        fixture_case("coronene", fixtures::CORONENE, MatchingMode::Perfect, true),
        fixture_case(
            "azulene_nonalternant",
            fixtures::AZULENE,
            MatchingMode::Perfect,
            true,
        ),
        fixture_case(
            "c60_nonalternant",
            fixtures::FULLERENE_C60,
            MatchingMode::Perfect,
            true,
        ),
        fixture_case(
            "disconnected_four_hexagons",
            fixtures::DISCONNECTED_CYCLES,
            MatchingMode::Perfect,
            false,
        ),
        fixture_case("ladder_2x4", fixtures::LADDER, MatchingMode::Perfect, true),
        fixture_case(
            "grid_3x3_maximum",
            fixtures::GRID,
            MatchingMode::Maximum,
            false,
        ),
        prescribed_hole_path(),
        complete_bipartite(6),
    ]
}

fn percentile(sorted: &[Duration], numerator: usize, denominator: usize) -> Duration {
    let index = sorted.len().saturating_sub(1).saturating_mul(numerator) / denominator;
    sorted[index]
}

fn report_delays(cases: &[MatchingCase]) {
    if env::var_os(DELAY_DIAGNOSTICS_ENV).is_none() {
        return;
    }

    eprintln!("case\toutputs\tfirst_ns\tmedian_ns\tp95_ns\tmax_ns");
    for case in cases {
        let mut delays = Vec::with_capacity(case.output_count);
        let started = Instant::now();
        let mut previous = started;
        let _ = case.visit(|_| {
            let now = Instant::now();
            delays.push(now.duration_since(previous));
            previous = now;
            ControlFlow::<()>::Continue(())
        });
        let first_delay = delays.first().copied();
        delays.sort_unstable();

        if delays.is_empty() {
            eprintln!("{}\t0\tNA\tNA\tNA\tNA", case.name);
            continue;
        }

        eprintln!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            case.name,
            delays.len(),
            first_delay.expect("nonempty delays").as_nanos(),
            percentile(&delays, 1, 2).as_nanos(),
            percentile(&delays, 95, 100).as_nanos(),
            delays.last().expect("nonempty delays").as_nanos(),
        );
    }
}

fn matching_enumeration(c: &mut Criterion) {
    let cases = corpus();
    report_delays(&cases);

    let mut first = c.benchmark_group("matching_first_output");
    for case in &cases {
        first.bench_function(case.name, |b| {
            b.iter(|| {
                let _ = case.visit(|matching| {
                    black_box(matching);
                    ControlFlow::Break(())
                });
            });
        });
    }
    first.finish();

    let mut prefixes = c.benchmark_group("matching_visit_prefix");
    for case in &cases {
        for limit in PREFIX_LIMITS {
            if case.output_count < limit {
                continue;
            }
            prefixes.bench_function(format!("{}/k_{limit}", case.name), |b| {
                b.iter(|| {
                    let mut visited = 0;
                    let _ = case.visit(|matching| {
                        black_box(matching);
                        visited += 1;
                        if visited == limit {
                            ControlFlow::Break(())
                        } else {
                            ControlFlow::Continue(())
                        }
                    });
                    black_box(visited)
                });
            });
        }
    }
    prefixes.finish();

    let mut full = c.benchmark_group("matching_visit_full");
    for case in &cases {
        full.throughput(Throughput::Elements(case.output_count as u64));
        full.bench_function(case.name, |b| {
            b.iter(|| {
                let mut visited = 0;
                let _ = case.visit(|_| {
                    visited += 1;
                    ControlFlow::<()>::Continue(())
                });
                black_box(visited)
            });
        });
    }
    full.finish();

    let mut eager = c.benchmark_group("matching_eager_collection");
    for case in &cases {
        eager.throughput(Throughput::Elements(case.output_count as u64));
        eager.bench_function(case.name, |b| b.iter(|| black_box(case.collect())));
    }
    eager.finish();

    let mut fkt = c.benchmark_group("matching_fkt_count");
    for case in &cases {
        let Some(embedding) = &case.embedding else {
            continue;
        };
        fkt.bench_function(case.name, |b| {
            b.iter(|| {
                black_box(
                    case.graph
                        .count_perfect_matchings_planar(embedding)
                        .expect("validated benchmark embedding must remain countable"),
                )
            });
        });
    }
    fkt.finish();
}

criterion_group!(benches, maximum_matching, matching_enumeration);
criterion_main!(benches);
