use criterion::{criterion_group, criterion_main, Criterion};

use umol_graph_core::{
    AutomorphismAlgorithm, BiconnectedComponentsAlgorithm, ConnectedComponentsAlgorithm,
    CycleEnumerationAlgorithm, EdgeId, Graph, MaxIndependentSetAlgorithm, MaxMatchingAlgorithm,
    MatchingEnumerationAlgorithm, NodeId, ShortestCycleAlgorithm, SubgraphIsomorphismAlgorithm,
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
    Graph::new(
        10,
        &[
            [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0],
            [5, 6], [6, 7], [7, 8], [8, 9], [9, 4],
        ],
    )
}

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
    #[rustfmt::skip]
    let edges: [[u32; 2]; 90] = [
        [0,1],[0,4],[0,5],[1,2],[1,10],[2,3],
        [2,15],[3,4],[3,20],[4,25],[5,6],[5,9],
        [6,7],[6,29],[7,8],[7,54],[8,9],[8,30],
        [9,11],[10,11],[10,14],[11,12],[12,13],[12,34],
        [13,14],[13,35],[14,16],[15,16],[15,19],[16,17],
        [17,18],[17,39],[18,19],[18,40],[19,21],[20,21],
        [20,24],[21,22],[22,23],[22,44],[23,24],[23,45],
        [24,26],[25,26],[25,29],[26,27],[27,28],[27,49],
        [28,29],[28,50],[30,31],[30,34],[31,32],[31,53],
        [32,33],[32,55],[33,34],[33,36],[35,36],[35,39],
        [36,37],[37,38],[37,59],[38,39],[38,41],[40,41],
        [40,44],[41,42],[42,43],[42,58],[43,44],[43,46],
        [45,46],[45,49],[46,47],[47,48],[47,57],[48,49],
        [48,51],[50,51],[50,54],[51,52],[52,53],[52,56],
        [53,54],[55,56],[55,59],[56,57],[57,58],[58,59],
    ];
    Graph::new(60, &edges)
}

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
    let graphs: Vec<(&str, Graph)> = vec![
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

    let mut group = c.benchmark_group("automorphism");
    for (name, g) in &graphs {
        group.bench_function(format!("{name}/uniform"), |b| {
            b.iter(|| g.automorphisms(|_: NodeId| 0u8, AutomorphismAlgorithm::Nauty));
        });
    }
    for (name, g) in &graphs {
        group.bench_function(format!("{name}/unique"), |b| {
            b.iter(|| {
                g.automorphisms(|n: NodeId| n.index() as u32, AutomorphismAlgorithm::Nauty)
            });
        });
    }
    group.finish();
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
        for (qname, query) in &queries {
            if query.node_count() > target.node_count() {
                continue;
            }
            group.bench_function(format!("{tname}/{qname}"), |b| {
                b.iter(|| {
                    target.subgraph_isomorphisms(
                        query,
                        &mut |_: NodeId, _: NodeId| true,
                        &mut |_: EdgeId, _: EdgeId| true,
                        SubgraphIsomorphismAlgorithm::Vf2,
                    )
                });
            });
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
    subgraph_isomorphism,
    matching_enumeration,
);
criterion_main!(benches);
