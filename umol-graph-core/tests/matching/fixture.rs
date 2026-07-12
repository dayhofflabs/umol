use umol_graph_core::Graph;

pub const BENZENE: &str = include_str!("data/benzene.graph");
pub const NAPHTHALENE: &str = include_str!("data/naphthalene.graph");
pub const CORONENE: &str = include_str!("data/coronene.graph");
pub const AZULENE: &str = include_str!("data/azulene.graph");
pub const FULLERENE_C60: &str = include_str!("data/c60.graph");
pub const DISCONNECTED_CYCLES: &str = include_str!("data/disconnected_cycles.graph");
pub const LADDER: &str = include_str!("data/ladder.graph");
pub const GRID: &str = include_str!("data/grid.graph");

#[derive(Clone, Debug, Eq, PartialEq)]
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
    enum Section {
        Header,
        Edges,
        Faces,
    }

    let mut section = Section::Header;
    let mut node_count = None;
    let mut edges = Vec::new();
    let mut faces = Vec::new();

    for line in source.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        match line {
            "edges" => section = Section::Edges,
            "faces" => section = Section::Faces,
            _ => match section {
                Section::Header => {
                    let value = line.strip_prefix("nodes ").expect("expected `nodes N`");
                    node_count = Some(value.parse().expect("node count must be an integer"));
                }
                Section::Edges => {
                    let values: Vec<u32> = line
                        .split_ascii_whitespace()
                        .map(|value| value.parse().expect("fixture entries must be integers"))
                        .collect();
                    let edge: [u32; 2] = values.try_into().expect("edges need two endpoints");
                    edges.push(edge);
                }
                Section::Faces => {
                    let values: Vec<u32> = line
                        .split_ascii_whitespace()
                        .map(|value| value.parse().expect("fixture entries must be integers"))
                        .collect();
                    assert!(values.len() >= 3, "faces need at least three vertices");
                    faces.push(values);
                }
            },
        }
    }

    GraphFixture {
        node_count: node_count.expect("fixture must declare its node count"),
        edges,
        faces,
    }
}
