use std::hint::black_box;
use std::iter;
use std::ops::ControlFlow;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use umol_graph_core::SubgraphIsomorphismAlgorithm::{
    ArcMatch, RayKirsch, Ri, Ullmann, Vf2, Vf2Rdkit,
};
use umol_graph_core::{
    AutomorphismAlgorithm, BiconnectedComponentsAlgorithm, BipartiteMaximumMatchingAlgorithm,
    CommonSubgraphEnumerationAlgorithm, ConnectedComponentsAlgorithm, EdgeId, EmbeddingKind,
    GeneralMaximumMatchingAlgorithm, Graph, MaximumIndependentSetAlgorithm,
    MinimumCycleBasisAlgorithm, NodeId, RelevantCycleEnumerationAlgorithm, ShortestCycleAlgorithm,
    SimpleCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm, UniqueRingFamilyAlgorithm,
    ARCMATCH_DEFAULT_PATH_LENGTH,
};

mod matching_graphs {
    use serde::Deserialize;
    use umol_graph_core::Graph;

    pub const BENZENE: &str = include_str!("matching/data/benzene_planar.toml");
    pub const NAPHTHALENE: &str = include_str!("matching/data/naphthalene_planar.toml");
    pub const CORONENE: &str = include_str!("matching/data/coronene_planar.toml");
    pub const AZULENE: &str = include_str!("matching/data/azulene_planar.toml");
    pub const FULLERENE_C60: &str = include_str!("matching/data/fullerene_c60_planar.toml");
    pub const DISCONNECTED_CYCLES: &str =
        include_str!("matching/data/four_disconnected_hexagons.toml");
    pub const LADDER: &str = include_str!("matching/data/ladder_2x4_planar.toml");
    pub const GRID: &str = include_str!("matching/data/grid_3x3_planar.toml");

    #[derive(Deserialize)]
    pub struct GraphFixture {
        pub node_count: usize,
        pub edges: Vec<[u32; 2]>,
        pub faces: Vec<Vec<u32>>,
    }

    impl GraphFixture {
        pub fn graph(&self) -> Graph {
            Graph::new(self.node_count, &self.edges)
        }
    }

    pub fn parse(source: &str) -> GraphFixture {
        toml::from_str(source).expect("matching benchmark fixture must be valid TOML")
    }
}

fn path(n: usize) -> Graph {
    let edges: Vec<[u32; 2]> = (0..n as u32 - 1).map(|i| [i, i + 1]).collect();
    Graph::new(n, &edges)
}

fn cycle(n: usize) -> Graph {
    let mut edges: Vec<[u32; 2]> = (0..n as u32 - 1).map(|i| [i, i + 1]).collect();
    edges.push([n as u32 - 1, 0]);
    Graph::new(n, &edges)
}

fn naphthalene() -> Graph {
    matching_graphs::parse(matching_graphs::NAPHTHALENE).graph()
}

#[rustfmt::skip]
fn prismane() -> Graph {
    Graph::new(
        6,
        &[
            [0, 1], [1, 2], [2, 0],
            [3, 4], [4, 5], [5, 3],
            [0, 3], [1, 4], [2, 5],
        ],
    )
}

#[rustfmt::skip]
fn cubane() -> Graph {
    Graph::new(
        8,
        &[
            [0, 1], [1, 2], [2, 3], [3, 0],
            [4, 5], [5, 6], [6, 7], [7, 4],
            [0, 4], [1, 5], [2, 6], [3, 7],
        ],
    )
}

#[rustfmt::skip]
fn adamantane() -> Graph {
    Graph::new(
        10,
        &[
            [0, 1], [0, 2], [0, 3],
            [1, 4], [1, 5],
            [2, 4], [2, 6],
            [3, 5], [3, 6],
            [4, 7],
            [5, 8],
            [6, 9],
            [7, 8], [7, 9], [8, 9],
        ],
    )
}

#[rustfmt::skip]
fn dodecahedron() -> Graph {
    Graph::new(
        20,
        &[
            [0, 1], [1, 2], [2, 3], [3, 4], [4, 0],
            [0, 5], [1, 6], [2, 7], [3, 8], [4, 9],
            [5, 10], [5, 14], [6, 10], [6, 11], [7, 11], [7, 12],
            [8, 12], [8, 13], [9, 13], [9, 14],
            [10, 15], [11, 16], [12, 17], [13, 18], [14, 19],
            [15, 16], [16, 17], [17, 18], [18, 19], [19, 15],
        ],
    )
}

#[rustfmt::skip]
fn icosahedron() -> Graph {
    Graph::new(
        12,
        &[
            [0, 1], [0, 2], [0, 3], [0, 4], [0, 5],
            [1, 2], [2, 3], [3, 4], [4, 5], [5, 1],
            [1, 6], [2, 6], [2, 7], [3, 7], [3, 8],
            [4, 8], [4, 9], [5, 9], [5, 10], [1, 10],
            [6, 7], [7, 8], [8, 9], [9, 10], [10, 6],
            [6, 11], [7, 11], [8, 11], [9, 11], [10, 11],
        ],
    )
}

fn fullerene_c60() -> Graph {
    matching_graphs::parse(matching_graphs::FULLERENE_C60).graph()
}

#[rustfmt::skip]
fn fullerene_c70() -> Graph {
    #[rustfmt::skip]
    let c70_edges: [[u32; 2]; 105] = [
        [0,1],[1,2],[2,3],[3,4],[4,0],
        [0,5],[1,7],[2,9],[3,11],[4,13],
        [5,6],[6,7],[7,8],[8,9],[9,10],[10,11],[11,12],[12,13],[13,14],[14,5],
        [6,15],[8,17],[10,19],[12,21],[14,23],
        [15,16],[16,17],[17,18],[18,19],[19,20],[20,21],[21,22],[22,23],[23,24],[24,15],
        [16,25],[18,27],[20,29],[22,31],[24,33],
        [25,26],[26,27],[27,28],[28,29],[29,30],[30,31],[31,32],[32,33],[33,34],[34,25],
        [26,35],[28,37],[30,39],[32,41],[34,43],
        [35,36],[36,37],[37,38],[38,39],[39,40],[40,41],[41,42],[42,43],[43,44],[44,35],
        [36,45],[38,47],[40,49],[42,51],[44,53],
        [45,46],[46,47],[47,48],[48,49],[49,50],[50,51],[51,52],[52,53],[53,54],[54,45],
        [46,55],[48,57],[50,59],[52,61],[54,63],
        [55,56],[56,57],[57,58],[58,59],[59,60],[60,61],[61,62],[62,63],[63,64],[64,55],
        [56,65],[58,66],[60,67],[62,68],[64,69],
        [65,66],[66,67],[67,68],[68,69],[69,65],
    ];

    Graph::new(70, &c70_edges)
}

#[rustfmt::skip]
fn petersen() -> Graph {
    Graph::new(
        10,
        &[
            [0, 1], [1, 2], [2, 3], [3, 4], [4, 0],
            [5, 7], [7, 9], [9, 6], [6, 8], [8, 5],
            [0, 5], [1, 6], [2, 7], [3, 8], [4, 9],
        ],
    )
}

fn complete(n: usize) -> Graph {
    let mut edges = Vec::new();
    for i in 0..n as u32 {
        for j in i + 1..n as u32 {
            edges.push([i, j]);
        }
    }
    Graph::new(n, &edges)
}

fn grid(rows: usize, columns: usize) -> Graph {
    let mut edges = Vec::new();
    for row in 0..rows {
        for column in 0..columns {
            let node = (row * columns + column) as u32;
            if row + 1 < rows {
                edges.push([node, node + columns as u32]);
            }
            if column + 1 < columns {
                edges.push([node, node + 1]);
            }
        }
    }
    Graph::new(rows * columns, &edges)
}

fn hypercube(dimension: usize) -> Graph {
    let node_count = 1usize << dimension;
    let mut edges = Vec::new();
    for node in 0..node_count {
        for bit in 0..dimension {
            let neighbor = node ^ (1 << bit);
            if node < neighbor {
                edges.push([node as u32, neighbor as u32]);
            }
        }
    }
    Graph::new(node_count, &edges)
}

fn disconnected_cycles() -> Graph {
    matching_graphs::parse(matching_graphs::DISCONNECTED_CYCLES).graph()
}

fn degree_colors(graph: &Graph) -> Vec<u32> {
    graph
        .node_ids()
        .map(|node| graph.neighbors(node).len() as u32)
        .collect()
}

fn incidence_colors(node_count: usize, edge_count: usize) -> Vec<u32> {
    iter::repeat_n(0, node_count)
        .chain(iter::repeat_n(1, edge_count))
        .collect()
}

fn cycle_corpus() -> Vec<(&'static str, Graph)> {
    vec![
        ("path_6", path(6)),
        ("hexagon", cycle(6)),
        ("naphthalene", naphthalene()),
        ("prismane", prismane()),
        ("cubane", cubane()),
        ("adamantane", adamantane()),
        ("dodecahedron", dodecahedron()),
        ("icosahedron", icosahedron()),
        ("c60", fullerene_c60()),
        ("c70", fullerene_c70()),
    ]
}

fn relevant_cycle_enumeration(c: &mut Criterion) {
    let mut group = c.benchmark_group("relevant_cycles");
    for (name, graph) in cycle_corpus() {
        group.bench_with_input(BenchmarkId::new("vismara", name), &graph, |b, graph| {
            b.iter(|| {
                graph.enumerate_relevant_cycles(
                    usize::MAX,
                    RelevantCycleEnumerationAlgorithm::Vismara,
                )
            });
        });
    }
    group.finish();

    let mut paths = c.benchmark_group("relevant_cycle_paths");
    for (name, graph) in cycle_corpus() {
        paths.bench_with_input(
            BenchmarkId::new("vismara_direct", name),
            &graph,
            |b, graph| {
                b.iter(|| {
                    graph
                        .try_enumerate_relevant_cycles(
                            usize::MAX,
                            RelevantCycleEnumerationAlgorithm::Vismara,
                        )
                        .expect("cycle corpus graphs are simple")
                });
            },
        );
        paths.bench_with_input(
            BenchmarkId::new("vismara_fallback", name),
            &graph,
            |b, graph| {
                b.iter(|| {
                    graph.enumerate_relevant_cycles_fallback(
                        usize::MAX,
                        RelevantCycleEnumerationAlgorithm::Vismara,
                    )
                });
            },
        );
        paths.bench_with_input(
            BenchmarkId::new("vismara_total", name),
            &graph,
            |b, graph| {
                b.iter(|| {
                    graph.enumerate_relevant_cycles(
                        usize::MAX,
                        RelevantCycleEnumerationAlgorithm::Vismara,
                    )
                });
            },
        );
    }
    paths.finish();
}

fn simple_cycle_enumeration(c: &mut Criterion) {
    let mut bounded = c.benchmark_group("simple_cycles_bounded_8");
    for (name, graph) in cycle_corpus() {
        bounded.bench_with_input(BenchmarkId::new("read_tarjan", name), &graph, |b, graph| {
            b.iter(|| {
                graph.enumerate_simple_cycles(8, SimpleCycleEnumerationAlgorithm::ReadTarjan)
            });
        });
    }
    bounded.finish();

    let mut paths = c.benchmark_group("simple_cycle_paths_bounded_8");
    for (name, graph) in cycle_corpus() {
        paths.bench_with_input(
            BenchmarkId::new("read_tarjan_direct", name),
            &graph,
            |b, graph| {
                b.iter(|| {
                    graph
                        .try_enumerate_simple_cycles(8, SimpleCycleEnumerationAlgorithm::ReadTarjan)
                        .expect("cycle corpus graphs are simple")
                });
            },
        );
        paths.bench_with_input(
            BenchmarkId::new("read_tarjan_fallback", name),
            &graph,
            |b, graph| {
                b.iter(|| {
                    graph.enumerate_simple_cycles_fallback(
                        8,
                        SimpleCycleEnumerationAlgorithm::ReadTarjan,
                    )
                });
            },
        );
        paths.bench_with_input(
            BenchmarkId::new("read_tarjan_total", name),
            &graph,
            |b, graph| {
                b.iter(|| {
                    graph.enumerate_simple_cycles(8, SimpleCycleEnumerationAlgorithm::ReadTarjan)
                });
            },
        );
    }
    paths.finish();

    let cases = [
        ("path_6", path(6)),
        ("hexagon", cycle(6)),
        ("naphthalene", naphthalene()),
        ("prismane", prismane()),
        ("cubane", cubane()),
        ("adamantane", adamantane()),
    ];
    let mut unbounded = c.benchmark_group("simple_cycles_unbounded");
    for (name, graph) in cases {
        unbounded.bench_with_input(BenchmarkId::new("read_tarjan", name), &graph, |b, graph| {
            b.iter(|| {
                graph.enumerate_simple_cycles(
                    usize::MAX,
                    SimpleCycleEnumerationAlgorithm::ReadTarjan,
                )
            });
        });
    }
    unbounded.finish();
}

fn minimum_cycle_basis(c: &mut Criterion) {
    let mut group = c.benchmark_group("minimum_cycle_basis");
    for (name, graph) in cycle_corpus() {
        group.bench_with_input(BenchmarkId::new("horton", name), &graph, |b, graph| {
            b.iter(|| graph.minimum_cycle_basis(MinimumCycleBasisAlgorithm::Horton));
        });
    }
    group.finish();
}

fn unique_ring_families(c: &mut Criterion) {
    let mut decomposition = c.benchmark_group("unique_ring_families/decomposition");
    for (name, graph) in cycle_corpus() {
        decomposition.bench_with_input(BenchmarkId::new("kolodzik", name), &graph, |b, graph| {
            b.iter(|| graph.unique_ring_families(UniqueRingFamilyAlgorithm::Kolodzik));
        });
    }
    decomposition.finish();

    let families = cycle_corpus()
        .into_iter()
        .map(|(name, graph)| {
            (
                name,
                graph.unique_ring_families(UniqueRingFamilyAlgorithm::Kolodzik),
            )
        })
        .collect::<Vec<_>>();
    let mut emission = c.benchmark_group("unique_ring_families/lazy_emission");
    for (name, families) in &families {
        emission.bench_with_input(
            BenchmarkId::new("kolodzik", name),
            families,
            |b, families| {
                b.iter(|| {
                    let mut count = 0;
                    for id in families.ids() {
                        let _: ControlFlow<()> = families.visit_relevant_cycles(id, |_| {
                            count += 1;
                            ControlFlow::Continue(())
                        });
                    }
                    count
                });
            },
        );
    }
    emission.finish();
}

fn shortest_cycle(c: &mut Criterion) {
    let graphs = [
        ("hexagon", cycle(6)),
        ("naphthalene", naphthalene()),
        ("cubane", cubane()),
        ("dodecahedron", dodecahedron()),
        ("c60", fullerene_c60()),
    ];

    let mut group = c.benchmark_group("shortest_cycle");
    for (name, g) in &graphs {
        let eid = EdgeId(0);
        group.bench_function(format!("{name}/edge"), |b| {
            b.iter(|| g.shortest_cycle_through_edge(eid, ShortestCycleAlgorithm::Bfs));
        });
        let nid = NodeId(0);
        group.bench_function(format!("{name}/node"), |b| {
            b.iter(|| g.shortest_cycle_through_node(nid, ShortestCycleAlgorithm::Bfs));
        });
    }
    group.finish();
}

fn connected_components(c: &mut Criterion) {
    let graphs = [
        ("hexagon", cycle(6)),
        ("cubane", cubane()),
        ("dodecahedron", dodecahedron()),
        ("c60", fullerene_c60()),
    ];

    let mut group = c.benchmark_group("connected_components");
    for (name, g) in &graphs {
        group.bench_function(*name, |b| {
            b.iter(|| g.enumerate_connected_components(ConnectedComponentsAlgorithm::Bfs));
        });
    }
    group.finish();
}

fn biconnected_components(c: &mut Criterion) {
    let graphs = [
        ("hexagon", cycle(6)),
        ("naphthalene", naphthalene()),
        ("cubane", cubane()),
        ("dodecahedron", dodecahedron()),
        ("c60", fullerene_c60()),
    ];

    let mut group = c.benchmark_group("biconnected_components");
    for (name, g) in &graphs {
        group.bench_function(*name, |b| {
            b.iter(|| g.enumerate_biconnected_components(BiconnectedComponentsAlgorithm::Tarjan));
        });
    }
    group.finish();
}

fn maximum_matching(c: &mut Criterion) {
    let graphs = [
        ("hexagon", cycle(6)),
        ("cubane", cubane()),
        ("petersen", petersen()),
        ("dodecahedron", dodecahedron()),
        ("c60", fullerene_c60()),
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
        ("hexagon", cycle(6)),
        ("cubane", cubane()),
        ("grid_16x16", grid(16, 16)),
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

fn maximum_independent_set(c: &mut Criterion) {
    let graphs = [
        ("path_6", path(6)),
        ("hexagon", cycle(6)),
        ("cubane", cubane()),
        ("petersen", petersen()),
        ("dodecahedron", dodecahedron()),
    ];

    let mut group = c.benchmark_group("maximum_independent_set");
    for (name, g) in &graphs {
        group.bench_function(*name, |b| {
            b.iter(|| g.maximum_independent_set(MaximumIndependentSetAlgorithm::BranchAndBound));
        });
    }
    group.finish();
}

fn automorphism(c: &mut Criterion) {
    let molecular: Vec<(&str, Graph)> = vec![
        ("path_6", path(6)),
        ("hexagon", cycle(6)),
        ("naphthalene", naphthalene()),
        ("prismane", prismane()),
        ("cubane", cubane()),
        ("adamantane", adamantane()),
        ("petersen", petersen()),
        ("dodecahedron", dodecahedron()),
        ("icosahedron", icosahedron()),
        ("c60", fullerene_c60()),
        ("c70", fullerene_c70()),
        ("K5", complete(5)),
        ("K8", complete(8)),
    ];

    let stress: Vec<(&str, Graph)> = vec![
        ("path_64", path(64)),
        ("cycle_64", cycle(64)),
        ("grid_8x8", grid(8, 8)),
        ("hypercube_6", hypercube(6)),
        ("four_hexagons", disconnected_cycles()),
        ("K10", complete(10)),
    ];

    let incidence = [("adamantane", adamantane()), ("c60", fullerene_c60())]
        .into_iter()
        .map(|(name, graph)| {
            let colors = incidence_colors(graph.node_count(), graph.edge_count());
            (name, graph.subdivide_edges(), colors)
        })
        .collect::<Vec<_>>();

    let mut group = c.benchmark_group("automorphism");
    for (name, graph) in molecular.iter().chain(&stress) {
        group.bench_function(format!("ordinary/{name}/uniform"), |b| {
            b.iter(|| graph.automorphisms(|_: NodeId| 0u32, AutomorphismAlgorithm::Nauty));
        });
        let colors = degree_colors(graph);
        group.bench_function(format!("ordinary/{name}/degree"), |b| {
            b.iter(|| {
                graph.automorphisms(|node| colors[node.index()], AutomorphismAlgorithm::Nauty)
            });
        });
        group.bench_function(format!("ordinary/{name}/unique"), |b| {
            b.iter(|| {
                graph.automorphisms(|node| node.index() as u32, AutomorphismAlgorithm::Nauty)
            });
        });
    }
    for (name, subdivision, colors) in &incidence {
        let graph = subdivision.graph();
        group.bench_function(format!("incidence/{name}/class"), |b| {
            b.iter(|| {
                graph.automorphisms(|node| colors[node.index()], AutomorphismAlgorithm::Nauty)
            });
        });
    }
    group.finish();
}

fn automorphism_stabilizer(c: &mut Criterion) {
    let cases: Vec<(&str, Graph, Vec<u32>, Vec<NodeId>)> = [
        ("adamantane", adamantane()),
        ("dodecahedron", dodecahedron()),
        ("c60", fullerene_c60()),
        ("four_hexagons", disconnected_cycles()),
    ]
    .into_iter()
    .map(|(name, graph)| {
        let colors = degree_colors(&graph);
        let sites = [0, graph.node_count() / 3, graph.node_count() / 2]
            .into_iter()
            .map(|site| NodeId(site as u32))
            .collect();
        (name, graph, colors, sites)
    })
    .collect();

    let mut group = c.benchmark_group("automorphism_stabilizer");
    for (name, graph, colors, sites) in &cases {
        group.bench_function(format!("{name}/three_sites"), |b| {
            b.iter(|| {
                sites
                    .iter()
                    .map(|&site| {
                        graph.automorphisms(
                            |node| (node == site, colors[node.index()]),
                            AutomorphismAlgorithm::Nauty,
                        )
                    })
                    .collect::<Vec<_>>()
            });
        });
    }
    group.finish();
}

fn canonical_key(c: &mut Criterion) {
    let cases = [
        ("naphthalene", naphthalene()),
        ("adamantane", adamantane()),
        ("dodecahedron", dodecahedron()),
        ("c60", fullerene_c60()),
        ("c70", fullerene_c70()),
        ("grid_8x8", grid(8, 8)),
    ];

    let mut group = c.benchmark_group("canonical_key");
    for (name, graph) in &cases {
        let node_colors: Vec<Vec<u8>> = degree_colors(graph)
            .into_iter()
            .map(|color| color.to_le_bytes().to_vec())
            .collect();
        let edge_colors: Vec<Vec<u8>> = graph
            .edge_ids()
            .map(|edge| vec![(edge.index() % 3) as u8])
            .collect();
        group.bench_function(*name, |b| {
            b.iter(|| {
                graph.canonical_key(
                    |node| node_colors[node.index()].clone(),
                    |edge| edge_colors[edge.index()].clone(),
                    AutomorphismAlgorithm::Nauty,
                )
            });
        });
    }
    group.finish();
}

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

fn subiso_name(algorithm: SubgraphIsomorphismAlgorithm) -> &'static str {
    match algorithm {
        Vf2 => "vf2",
        Ullmann => "ullmann",
        Ri => "ri",
        ArcMatch { .. } => "arcmatch",
        Vf2Rdkit => "vf2rdkit",
        RayKirsch => "raykirsch",
    }
}

// Alternating synthetic bond orders so the edge-label-aware algorithms (ArcMatch's
// edge domains, RI's ordering) are exercised — unlabeled regular graphs hide their
// advantage. Values are deterministic, not chemically meaningful.
fn edge_labels(graph: &Graph) -> Vec<u8> {
    (0..graph.edge_count()).map(|i| (i % 2) as u8).collect()
}

fn subgraph_isomorphism(c: &mut Criterion) {
    let targets: Vec<(&str, Graph)> = vec![
        ("hexagon", cycle(6)),
        ("naphthalene", naphthalene()),
        ("cubane", cubane()),
        ("adamantane", adamantane()),
        ("petersen", petersen()),
        ("dodecahedron", dodecahedron()),
        ("c60", fullerene_c60()),
    ];

    let queries: Vec<(&str, Graph)> = vec![
        ("edge", Graph::new(2, &[[0, 1]])),
        ("path_3", path(3)),
        ("triangle", cycle(3)),
        ("square", cycle(4)),
        ("hexagon", cycle(6)),
    ];

    let mut group = c.benchmark_group("subgraph_isomorphism");
    for (tname, target) in &targets {
        let t_labels = edge_labels(target);
        for (qname, query) in &queries {
            if query.node_count() > target.node_count() {
                continue;
            }
            let q_labels = edge_labels(query);
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
                // Time to first result: break on the first emission.
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
                // Complete visitation without result storage: total enumeration
                // time minus the collection cost.
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
    group.finish();
}

struct CommonSubgraphCase {
    name: &'static str,
    expected_output_count: usize,
    left: Graph,
    right: Graph,
    left_node_labels: Vec<u8>,
    right_node_labels: Vec<u8>,
    left_edge_labels: Vec<u8>,
    right_edge_labels: Vec<u8>,
    embedding: EmbeddingKind,
}

fn common_subgraph_enumeration(c: &mut Criterion) {
    let triangle = cycle(3);
    let path_4 = path(4);
    let cycle_4 = cycle(4);
    let path_5 = path(5);
    let cycle_6 = cycle(6);
    let benzene = cycle(6);
    let naphthalene = naphthalene();

    let cases = [
        CommonSubgraphCase {
            name: "dense_compatible/triangle_triangle/induced",
            expected_output_count: 34,
            left_node_labels: vec![0; triangle.node_count()],
            right_node_labels: vec![0; triangle.node_count()],
            left_edge_labels: vec![0; triangle.edge_count()],
            right_edge_labels: vec![0; triangle.edge_count()],
            left: triangle.clone(),
            right: triangle,
            embedding: EmbeddingKind::Induced,
        },
        CommonSubgraphCase {
            name: "structural/path_4_cycle_4/induced",
            expected_output_count: 69,
            left_node_labels: vec![0; path_4.node_count()],
            right_node_labels: vec![0; cycle_4.node_count()],
            left_edge_labels: vec![0; path_4.edge_count()],
            right_edge_labels: vec![0; cycle_4.edge_count()],
            left: path_4.clone(),
            right: cycle_4.clone(),
            embedding: EmbeddingKind::Induced,
        },
        CommonSubgraphCase {
            name: "structural/path_4_cycle_4/monomorphism",
            expected_output_count: 209,
            left_node_labels: vec![0; path_4.node_count()],
            right_node_labels: vec![0; cycle_4.node_count()],
            left_edge_labels: vec![0; path_4.edge_count()],
            right_edge_labels: vec![0; cycle_4.edge_count()],
            left: path_4,
            right: cycle_4,
            embedding: EmbeddingKind::Monomorphism,
        },
        CommonSubgraphCase {
            name: "label_selective/path_5_cycle_6/induced",
            expected_output_count: 109,
            left_node_labels: vec![0, 1, 0, 1, 0],
            right_node_labels: vec![0, 1, 0, 1, 0, 1],
            left_edge_labels: edge_labels(&path_5),
            right_edge_labels: edge_labels(&cycle_6),
            left: path_5,
            right: cycle_6,
            embedding: EmbeddingKind::Induced,
        },
        CommonSubgraphCase {
            name: "molecular/benzene_naphthalene/induced",
            expected_output_count: 1_957,
            left_node_labels: benzene
                .node_ids()
                .map(|node| benzene.degree(node) as u8)
                .collect(),
            right_node_labels: naphthalene
                .node_ids()
                .map(|node| naphthalene.degree(node) as u8)
                .collect(),
            left_edge_labels: vec![0; benzene.edge_count()],
            right_edge_labels: vec![0; naphthalene.edge_count()],
            left: benzene,
            right: naphthalene,
            embedding: EmbeddingKind::Induced,
        },
    ];

    let mut group = c.benchmark_group("common_subgraph_enumeration");
    for case in cases {
        for (algorithm_name, algorithm) in [
            (
                "modular_product_backtracking",
                CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
            ),
            (
                "direct_backtracking",
                CommonSubgraphEnumerationAlgorithm::DirectBacktracking,
            ),
        ] {
            let mut node_match = |left: NodeId, right: NodeId| {
                case.left_node_labels[left.index()] == case.right_node_labels[right.index()]
            };
            let mut edge_match = |left: EdgeId, right: EdgeId| {
                case.left_edge_labels[left.index()] == case.right_edge_labels[right.index()]
            };
            let output = case.left.enumerate_common_subgraphs(
                &case.right,
                &mut node_match,
                &mut edge_match,
                case.embedding,
                algorithm,
            );
            assert_eq!(
                output.len(),
                case.expected_output_count,
                "unexpected output count for {} with {algorithm_name}",
                case.name,
            );

            group.bench_with_input(
                BenchmarkId::new(case.name, algorithm_name),
                &algorithm,
                |b, &algorithm| {
                    b.iter(|| {
                        let mut node_match = |left: NodeId, right: NodeId| {
                            case.left_node_labels[left.index()]
                                == case.right_node_labels[right.index()]
                        };
                        let mut edge_match = |left: EdgeId, right: EdgeId| {
                            case.left_edge_labels[left.index()]
                                == case.right_edge_labels[right.index()]
                        };
                        black_box(case.left.enumerate_common_subgraphs(
                            &case.right,
                            &mut node_match,
                            &mut edge_match,
                            case.embedding,
                            algorithm,
                        ))
                    });
                },
            );
        }
    }
    group.finish();
}

mod matching {
    use std::env;
    use std::hint::black_box;
    use std::ops::ControlFlow;
    use std::time::{Duration, Instant};

    use criterion::{Criterion, Throughput};
    use umol_graph_core::{
        EdgeId, FaceBoundary, Graph, Matching, MatchingEnumerationAlgorithm, NodeId,
        PlanarEmbedding,
    };

    use super::matching_graphs as fixture;

    const PREFIX_LIMITS: [usize; 3] = [1, 10, 100];
    const DELAY_DIAGNOSTICS_ENV: &str = "UMOL_MATCHING_DELAY_DIAGNOSTICS";

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

    fn embedding(graph: &Graph, fixture: &fixture::GraphFixture) -> PlanarEmbedding {
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
        let parsed = fixture::parse(source);
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
                let _ = graph.visit_perfect_matchings(
                    MatchingEnumerationAlgorithm::BranchAndBound,
                    |_| {
                        count += 1;
                        ControlFlow::<()>::Continue(())
                    },
                );
            }
            MatchingMode::Maximum => {
                let _ = graph.visit_maximum_matchings(
                    MatchingEnumerationAlgorithm::BranchAndBound,
                    |_| {
                        count += 1;
                        ControlFlow::<()>::Continue(())
                    },
                );
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
            fixture_case("benzene", fixture::BENZENE, MatchingMode::Perfect, true),
            fixture_case(
                "naphthalene",
                fixture::NAPHTHALENE,
                MatchingMode::Perfect,
                true,
            ),
            fixture_case("coronene", fixture::CORONENE, MatchingMode::Perfect, true),
            fixture_case(
                "azulene_nonalternant",
                fixture::AZULENE,
                MatchingMode::Perfect,
                true,
            ),
            fixture_case(
                "c60_nonalternant",
                fixture::FULLERENE_C60,
                MatchingMode::Perfect,
                true,
            ),
            fixture_case(
                "disconnected_four_hexagons",
                fixture::DISCONNECTED_CYCLES,
                MatchingMode::Perfect,
                false,
            ),
            fixture_case("ladder_2x4", fixture::LADDER, MatchingMode::Perfect, true),
            fixture_case(
                "grid_3x3_maximum",
                fixture::GRID,
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

    pub(super) fn matching_enumeration(c: &mut Criterion) {
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
}

criterion_group!(
    benches,
    relevant_cycle_enumeration,
    simple_cycle_enumeration,
    minimum_cycle_basis,
    unique_ring_families,
    shortest_cycle,
    connected_components,
    biconnected_components,
    maximum_matching,
    maximum_independent_set,
    automorphism,
    automorphism_stabilizer,
    canonical_key,
    subgraph_isomorphism,
    common_subgraph_enumeration,
    matching::matching_enumeration,
);
criterion_main!(benches);
