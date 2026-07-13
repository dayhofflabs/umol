use serde::Deserialize;
use umol_graph_core::Graph;

pub const BENZENE: &str = include_str!("data/benzene_planar.toml");
pub const NAPHTHALENE: &str = include_str!("data/naphthalene_planar.toml");
pub const CORONENE: &str = include_str!("data/coronene_planar.toml");
pub const AZULENE: &str = include_str!("data/azulene_planar.toml");
pub const FULLERENE_C60: &str = include_str!("data/fullerene_c60_planar.toml");
pub const DISCONNECTED_CYCLES: &str = include_str!("data/four_disconnected_hexagons.toml");
pub const LADDER: &str = include_str!("data/ladder_2x4_planar.toml");
pub const GRID: &str = include_str!("data/grid_3x3_planar.toml");

#[derive(Clone, Debug, Deserialize, Eq, PartialEq)]
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
    toml::from_str(source).expect("matching fixture must be valid TOML")
}
