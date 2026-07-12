use std::iter;

use criterion::{criterion_group, criterion_main, Criterion};
#[path = "../tests/matching/fixture.rs"]
#[allow(dead_code)]
mod matching_graphs;
use umol_graph_core::SubgraphIsomorphismAlgorithm::{
    ArcMatch, RayKirsch, Ri, Ullmann, Vf2, Vf2Rdkit,
};
use umol_graph_core::{
    AutomorphismAlgorithm, BiconnectedComponentsAlgorithm, ConnectedComponentsAlgorithm,
    CycleEnumerationAlgorithm, EdgeId, Graph, MatchingEnumerationAlgorithm,
    MaxIndependentSetAlgorithm, MaxMatchingAlgorithm, NodeId, ShortestCycleAlgorithm,
    SubgraphIsomorphismAlgorithm, ARCMATCH_DEFAULT_PATH_LENGTH,
};

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

fn subdivided(graph: &Graph) -> Graph {
    let node_count = graph.node_count();
    let mut edges = Vec::with_capacity(2 * graph.edge_count());
    for (position, edge) in graph.edge_ids().enumerate() {
        let [a, b] = graph.edge_endpoints(edge);
        let edge_node = (node_count + position) as u32;
        edges.push([a.0, edge_node]);
        edges.push([edge_node, b.0]);
    }
    Graph::new(node_count + graph.edge_count(), &edges)
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

fn cycle_enumeration(c: &mut Criterion) {
    let graphs = [
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
    ];

    let mut group = c.benchmark_group("cycle_enumeration");
    for (name, g) in &graphs {
        group.bench_function(*name, |b| {
            b.iter(|| g.enumerate_cycles(usize::MAX, CycleEnumerationAlgorithm::Vismara));
        });
    }
    group.finish();
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
            b.iter(|| g.connected_components(ConnectedComponentsAlgorithm::Bfs));
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
            b.iter(|| g.biconnected_components(BiconnectedComponentsAlgorithm::Tarjan));
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

    let mut group = c.benchmark_group("maximum_matching");
    for (name, g) in &graphs {
        group.bench_function(*name, |b| {
            b.iter(|| g.maximum_matching(MaxMatchingAlgorithm::Edmonds));
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
            b.iter(|| g.maximum_independent_set(MaxIndependentSetAlgorithm::BranchAndBound));
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

    let incidence: Vec<(&str, Graph, Vec<u32>)> =
        [("adamantane", adamantane()), ("c60", fullerene_c60())]
            .into_iter()
            .map(|(name, graph)| {
                let colors = incidence_colors(graph.node_count(), graph.edge_count());
                (name, subdivided(&graph), colors)
            })
            .collect();

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
    for (name, graph, colors) in &incidence {
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
                        target.subgraph_isomorphisms(
                            query,
                            &mut |_: NodeId, _: NodeId| true,
                            &mut |qe: EdgeId, he: EdgeId| {
                                q_labels[qe.index()] == t_labels[he.index()]
                            },
                            algorithm,
                        )
                    });
                });
            }
        }
    }
    group.finish();
}

fn matching_enumeration(c: &mut Criterion) {
    let graphs = [
        ("hexagon", cycle(6)),
        ("cubane", cubane()),
        ("prismane", prismane()),
        ("petersen", petersen()),
    ];

    let mut group = c.benchmark_group("matching_enumeration");
    for (name, g) in &graphs {
        group.bench_function(format!("{name}/perfect"), |b| {
            b.iter(|| g.enumerate_perfect_matchings(MatchingEnumerationAlgorithm::BranchAndBound));
        });
        group.bench_function(format!("{name}/maximum"), |b| {
            b.iter(|| g.enumerate_maximum_matchings(MatchingEnumerationAlgorithm::BranchAndBound));
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    cycle_enumeration,
    shortest_cycle,
    connected_components,
    biconnected_components,
    maximum_matching,
    maximum_independent_set,
    automorphism,
    automorphism_stabilizer,
    canonical_key,
    subgraph_isomorphism,
    matching_enumeration,
);
criterion_main!(benches);
