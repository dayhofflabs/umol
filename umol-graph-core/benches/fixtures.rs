//! Shared graph fixtures for `umol-graph-core` benchmarks.

use serde::Deserialize;
use umol_graph_core::Graph;

pub const BENZENE: &str = include_str!("matching/data/benzene_planar.toml");
pub const NAPHTHALENE: &str = include_str!("matching/data/naphthalene_planar.toml");
pub const CORONENE: &str = include_str!("matching/data/coronene_planar.toml");
pub const AZULENE: &str = include_str!("matching/data/azulene_planar.toml");
pub const FULLERENE_C60: &str = include_str!("matching/data/fullerene_c60_planar.toml");
pub const DISCONNECTED_CYCLES: &str = include_str!("matching/data/four_disconnected_hexagons.toml");
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
    toml::from_str(source).expect("benchmark fixture must be valid TOML")
}

pub fn path(n: usize) -> Graph {
    let edges: Vec<[u32; 2]> = (0..n as u32 - 1).map(|i| [i, i + 1]).collect();
    Graph::new(n, &edges)
}

pub fn cycle(n: usize) -> Graph {
    let mut edges: Vec<[u32; 2]> = (0..n as u32 - 1).map(|i| [i, i + 1]).collect();
    edges.push([n as u32 - 1, 0]);
    Graph::new(n, &edges)
}

pub fn complete(n: usize) -> Graph {
    let mut edges = Vec::new();
    for i in 0..n as u32 {
        for j in i + 1..n as u32 {
            edges.push([i, j]);
        }
    }
    Graph::new(n, &edges)
}

pub fn naphthalene() -> Graph {
    parse(NAPHTHALENE).graph()
}

#[rustfmt::skip]
pub fn prismane() -> Graph {
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
pub fn cubane() -> Graph {
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
pub fn adamantane() -> Graph {
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
pub fn dodecahedron() -> Graph {
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
pub fn icosahedron() -> Graph {
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

pub fn fullerene_c60() -> Graph {
    parse(FULLERENE_C60).graph()
}

#[rustfmt::skip]
pub fn fullerene_c70() -> Graph {
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
pub fn petersen() -> Graph {
    Graph::new(
        10,
        &[
            [0, 1], [1, 2], [2, 3], [3, 4], [4, 0],
            [5, 7], [7, 9], [9, 6], [6, 8], [8, 5],
            [0, 5], [1, 6], [2, 7], [3, 8], [4, 9],
        ],
    )
}

pub fn grid(rows: usize, columns: usize) -> Graph {
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

pub fn cycle_corpus() -> Vec<(&'static str, Graph)> {
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

// Alternating synthetic bond orders so edge-label-aware algorithms are exercised.
// Values are deterministic, not chemically meaningful.
pub fn alternating_edge_labels(graph: &Graph) -> Vec<u8> {
    (0..graph.edge_count()).map(|i| (i % 2) as u8).collect()
}
