//! Unified graph automorphism, canonical-labeling, stabilizer, and canonical-key benchmarks
//! comparing `AutomorphismAlgorithm::Nauty` against `AutomorphismAlgorithm::TrailIr` across:
//!
//! - `nodes_only`: vertex-colored graphs with uniform edge color (`automorphisms` and
//!   `canonical_labeling`).
//! - `nodes_and_edges`: vertex- and edge-colored graphs across IR incidence shapes, generated
//!   molecular skeletons, symmetric cages/fullerenes, combinatorial stress cases, and larger
//!   multi-ligand chemical complexes (`automorphisms` and `canonical_labeling`).
//! - `stabilizers`: repeated three-site point-stabilizer calls (`nodes_only` and
//!   `nodes_and_edges`, comparing `Nauty`, fresh `TrailIr`, and `AutomorphismScratch` reuse).
//! - `canonical_key`: byte-encoded `Graph::canonical_key` over node- and edge-colored graphs.

#[allow(dead_code)]
mod fixtures;

use std::array;
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashSet};
use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use rand::seq::SliceRandom;
use rand::{Rng, RngExt, SeedableRng};
use rand_chacha::ChaCha8Rng;
use umol_graph_core::{AutomorphismAlgorithm, AutomorphismScratch, Graph, NodeId};

struct NodeColoredCase {
    category: &'static str,
    name: &'static str,
    graph: Graph,
    node_colors: Vec<u32>,
}

impl NodeColoredCase {
    fn new(
        category: &'static str,
        name: &'static str,
        graph: Graph,
        node_colors: Vec<u32>,
    ) -> Self {
        assert_eq!(node_colors.len(), graph.node_count());
        Self {
            category,
            name,
            graph,
            node_colors,
        }
    }

    fn label(&self) -> String {
        format!(
            "{}/{} (n={},m={})",
            self.category,
            self.name,
            self.graph.node_count(),
            self.graph.edge_count()
        )
    }
}

struct EdgeColoredCase {
    category: &'static str,
    name: &'static str,
    graph: Graph,
    node_colors: Vec<u32>,
    edge_colors: Vec<u32>,
}

impl EdgeColoredCase {
    fn new(
        category: &'static str,
        name: &'static str,
        graph: Graph,
        node_colors: Vec<u32>,
        edge_colors: Vec<u32>,
    ) -> Self {
        assert_eq!(node_colors.len(), graph.node_count());
        assert_eq!(edge_colors.len(), graph.edge_count());
        Self {
            category,
            name,
            graph,
            node_colors,
            edge_colors,
        }
    }

    fn label(&self) -> String {
        format!(
            "{}/{} (n={},m={})",
            self.category,
            self.name,
            self.graph.node_count(),
            self.graph.edge_count()
        )
    }
}

fn grid_graph_with_axes(rows: usize, cols: usize) -> (Graph, Vec<u32>) {
    let mut edges = Vec::new();
    let mut axis_colors = Vec::new();
    for r in 0..rows as u32 {
        for c in 0..cols as u32 {
            let u = r * cols as u32 + c;
            if c + 1 < cols as u32 {
                edges.push([u, u + 1]);
                axis_colors.push(0);
            }
            if r + 1 < rows as u32 {
                edges.push([u, u + cols as u32]);
                axis_colors.push(1);
            }
        }
    }
    (Graph::new(rows * cols, &edges), axis_colors)
}

fn hypercube_graph_with_2colors(dimension: usize) -> (Graph, Vec<u32>) {
    let node_count = 1usize << dimension;
    let half = dimension / 2;
    let mut edges = Vec::new();
    let mut edge_colors = Vec::new();
    for u in 0..node_count {
        for bit in 0..dimension {
            let v = u ^ (1 << bit);
            if u < v {
                edges.push([u as u32, v as u32]);
                edge_colors.push(u32::from(bit >= half));
            }
        }
    }
    (Graph::new(node_count, &edges), edge_colors)
}

fn four_hexagons_with_kekule() -> (Graph, Vec<u32>) {
    let mut edges = Vec::new();
    let mut colors = Vec::new();
    for comp in 0..4u32 {
        let base = comp * 6;
        for i in 0..6u32 {
            edges.push([base + i, base + (i + 1) % 6]);
            colors.push(i % 2);
        }
    }
    (Graph::new(24, &edges), colors)
}

fn degree_colors(graph: &Graph) -> Vec<u32> {
    graph
        .node_ids()
        .map(|node| graph.neighbors(node).len() as u32)
        .collect()
}

fn edge_in_five_cycle(graph: &Graph, u: NodeId, v: NodeId) -> bool {
    for nu in graph.neighbors(u) {
        if nu.node == v {
            continue;
        }
        for nv in graph.neighbors(v) {
            if nv.node == u || nv.node == nu.node {
                continue;
            }
            for mid in graph.neighbors(nu.node) {
                if mid.node != u
                    && mid.node != v
                    && mid.node != nv.node
                    && graph.find_edge(mid.node, nv.node).is_some()
                {
                    return true;
                }
            }
        }
    }
    false
}

fn prufer_tree_plus_rings(
    node_count: usize,
    extra_edges: usize,
    max_degree: usize,
    rng: &mut impl Rng,
) -> Graph {
    let mut extra_degree = vec![0usize; node_count];
    let total_extra = node_count - 2;
    let cap = max_degree - 1;
    let mut placed = 0;
    while placed < total_extra {
        let v = rng.random_range(0..node_count);
        if extra_degree[v] < cap {
            extra_degree[v] += 1;
            placed += 1;
        }
    }

    let mut code = Vec::with_capacity(total_extra);
    for (v, &count) in extra_degree.iter().enumerate() {
        for _ in 0..count {
            code.push(v as u32);
        }
    }
    code.shuffle(rng);

    let mut remaining = vec![1usize; node_count];
    for &v in &code {
        remaining[v as usize] += 1;
    }

    let mut leaves: BinaryHeap<Reverse<u32>> = (0..node_count as u32)
        .filter(|&v| remaining[v as usize] == 1)
        .map(Reverse)
        .collect();

    let mut degrees = vec![0usize; node_count];
    let mut edge_set: HashSet<(u32, u32)> = HashSet::with_capacity(node_count - 1 + extra_edges);
    let mut edges = Vec::with_capacity(node_count - 1 + extra_edges);
    for &v in &code {
        let leaf = leaves.pop().expect("Prüfer tree has a leaf").0;
        let pair = (leaf.min(v), leaf.max(v));
        edge_set.insert(pair);
        edges.push([pair.0, pair.1]);
        degrees[leaf as usize] += 1;
        degrees[v as usize] += 1;
        remaining[leaf as usize] -= 1;
        remaining[v as usize] -= 1;
        if remaining[v as usize] == 1 {
            leaves.push(Reverse(v));
        }
    }
    let last_a = leaves.pop().expect("first final leaf").0;
    let last_b = leaves.pop().expect("second final leaf").0;
    let final_pair = (last_a.min(last_b), last_a.max(last_b));
    edge_set.insert(final_pair);
    edges.push([final_pair.0, final_pair.1]);
    degrees[last_a as usize] += 1;
    degrees[last_b as usize] += 1;

    let mut added = 0;
    let max_attempts = extra_edges.max(1) * node_count * 8;
    let mut attempts = 0;
    while added < extra_edges && attempts < max_attempts {
        attempts += 1;
        let a = rng.random_range(0..node_count) as u32;
        let b = rng.random_range(0..node_count) as u32;
        if a == b || degrees[a as usize] >= max_degree || degrees[b as usize] >= max_degree {
            continue;
        }
        let pair = (a.min(b), a.max(b));
        if edge_set.insert(pair) {
            edges.push([pair.0, pair.1]);
            degrees[a as usize] += 1;
            degrees[b as usize] += 1;
            added += 1;
        }
    }
    assert_eq!(added, extra_edges);
    Graph::new(node_count, &edges)
}

fn decorate_molecular_like(graph: &Graph, rng: &mut impl Rng) -> (Vec<u32>, Vec<u32>) {
    let node_colors = graph
        .node_ids()
        .map(|node| {
            let roll: u32 = rng.random_range(0..100);
            match graph.degree(node) {
                0 | 1 => match roll {
                    0..45 => 0,
                    45..70 => 1,
                    70..85 => 2,
                    85..95 => 3,
                    _ => 4,
                },
                2 => match roll {
                    0..70 => 0,
                    70..90 => 1,
                    _ => 2,
                },
                3 => u32::from(roll >= 80) * 2,
                _ => 0,
            }
        })
        .collect();
    let edge_colors = graph
        .edge_ids()
        .map(|edge| {
            let [a, b] = graph.edge_endpoints(edge);
            let max_deg = graph.degree(a).max(graph.degree(b));
            let roll: u32 = rng.random_range(0..100);
            if max_deg <= 2 && roll < 8 {
                2
            } else if max_deg <= 3 && roll < 26 {
                1
            } else {
                0
            }
        })
        .collect();
    (node_colors, edge_colors)
}

fn nodes_only_corpus() -> Vec<NodeColoredCase> {
    let mut cases = vec![
        NodeColoredCase::new(
            "ir",
            "propane",
            Graph::new(
                11,
                &[
                    [0, 8],
                    [1, 8],
                    [2, 8],
                    [3, 9],
                    [4, 9],
                    [5, 9],
                    [6, 10],
                    [7, 10],
                    [8, 10],
                    [9, 10],
                ],
            ),
            vec![0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1],
        ),
        NodeColoredCase::new(
            "ir",
            "ethane_radicals",
            Graph::new(8, &[[6, 7]]),
            vec![0, 0, 0, 0, 0, 0, 1, 1],
        ),
    ];

    let (grid_8x8, _) = grid_graph_with_axes(8, 8);
    let (hypercube_6, _) = hypercube_graph_with_2colors(6);
    let (four_hexagons, _) = four_hexagons_with_kekule();

    let uniform_graphs = [
        ("molecular", "naphthalene", fixtures::naphthalene()),
        ("cage", "cubane", fixtures::cubane()),
        ("cage", "adamantane", fixtures::adamantane()),
        ("cage", "petersen", fixtures::petersen()),
        ("cage", "dodecahedron", fixtures::dodecahedron()),
        ("cage", "c60", fixtures::fullerene_c60()),
        ("stress", "path_64", fixtures::path(64)),
        ("stress", "grid_8x8", grid_8x8.clone()),
        ("stress", "hypercube_6", hypercube_6),
        ("stress", "four_hexagons", four_hexagons),
        ("stress", "K10", fixtures::complete(10)),
    ];
    for (category, name, graph) in uniform_graphs {
        let n = graph.node_count();
        cases.push(NodeColoredCase::new(category, name, graph, vec![0; n]));
    }

    let degree_graphs = [
        ("degree", "naphthalene", fixtures::naphthalene()),
        ("degree", "grid_8x8", grid_8x8),
    ];
    for (category, name, graph) in degree_graphs {
        let colors = degree_colors(&graph);
        cases.push(NodeColoredCase::new(category, name, graph, colors));
    }

    let unique_graphs = [
        ("unique", "naphthalene", fixtures::naphthalene()),
        ("unique", "c60", fixtures::fullerene_c60()),
    ];
    for (category, name, graph) in unique_graphs {
        let colors = (0..graph.node_count() as u32).collect();
        cases.push(NodeColoredCase::new(category, name, graph, colors));
    }

    cases
}

fn ir_edge_colored_cases() -> Vec<EdgeColoredCase> {
    let mut cases = Vec::new();

    cases.push(EdgeColoredCase::new(
        "ir",
        "naphthalene_kekule",
        Graph::new(
            10,
            &[
                [0, 1],
                [1, 2],
                [2, 3],
                [3, 4],
                [4, 5],
                [5, 0],
                [5, 6],
                [6, 7],
                [7, 8],
                [8, 9],
                [9, 4],
            ],
        ),
        vec![0; 10],
        vec![1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0],
    ));

    let dmb_bonds: [[u32; 2]; 9] = [
        [0, 1],
        [0, 2],
        [0, 3],
        [0, 4],
        [1, 5],
        [1, 6],
        [2, 7],
        [2, 8],
        [3, 9],
    ];
    let mut levi_edges = Vec::new();
    let mut levi_edge_colors = Vec::new();
    let mut levi_node_colors = vec![0, 0, 0, 1, 2, 0, 0, 0, 0, 2];
    for (k, &[u, v]) in dmb_bonds.iter().enumerate() {
        let bnode = (10 + k) as u32;
        levi_node_colors.push(10);
        levi_edges.push([u, bnode]);
        levi_edge_colors.push(0);
        levi_edges.push([v, bnode]);
        levi_edge_colors.push(0);
    }
    levi_node_colors.push(20);
    levi_edges.push([0, 19]);
    levi_edge_colors.push(1);
    for lig in [1u32, 2, 3, 4] {
        levi_edges.push([lig, 19]);
        levi_edge_colors.push(2);
    }
    cases.push(EdgeColoredCase::new(
        "ir",
        "dimethylbutanol_stereo_levi",
        Graph::new(20, &levi_edges),
        levi_node_colors,
        levi_edge_colors,
    ));

    let mut ov_edges = Vec::new();
    let mut ov_edge_colors = Vec::new();
    let mut ov_node_colors = vec![0, 0, 0, 0, 1, 2, 3, 4];
    for &(u, v, bcol) in &[
        (0, 1, 10),
        (1, 2, 11),
        (2, 3, 10),
        (3, 0, 10),
        (1, 4, 10),
        (1, 5, 10),
        (2, 6, 10),
        (2, 7, 10),
    ] {
        let bn = ov_node_colors.len() as u32;
        ov_node_colors.push(bcol);
        ov_edges.push([u, bn]);
        ov_edge_colors.push(0);
        ov_edges.push([v, bn]);
        ov_edge_colors.push(0);
    }
    let dat = ov_node_colors.len() as u32;
    ov_node_colors.push(20);
    ov_edges.extend([[4, dat], [5, dat], [3, dat]]);
    ov_edge_colors.extend([1, 1, 2]);
    let aro = ov_node_colors.len() as u32;
    ov_node_colors.push(21);
    for (i, a) in [0u32, 1, 2, 3].into_iter().enumerate() {
        ov_edges.push([a, aro]);
        ov_edge_colors.push(if i == 0 { 3 } else { 4 });
    }
    let mc = ov_node_colors.len() as u32;
    ov_node_colors.push(22);
    for a in [0u32, 4, 5] {
        ov_edges.push([a, mc]);
        ov_edge_colors.push(5);
    }
    let nc = ov_node_colors.len() as u32;
    ov_node_colors.push(23);
    ov_edges.extend([[6, nc], [7, nc]]);
    ov_edge_colors.extend([0, 0]);
    let sa = ov_node_colors.len() as u32;
    ov_node_colors.push(24);
    ov_edges.push([1, sa]);
    ov_edge_colors.push(6);
    for a in [0u32, 2, 4, 5] {
        ov_edges.push([a, sa]);
        ov_edge_colors.push(7);
    }
    let sb = ov_node_colors.len() as u32;
    ov_node_colors.push(25);
    ov_edges.push([9, sb]);
    ov_edge_colors.push(6);
    for a in [0u32, 3, 4, 6] {
        ov_edges.push([a, sb]);
        ov_edge_colors.push(7);
    }
    cases.push(EdgeColoredCase::new(
        "ir",
        "overlay_heavy_incidence",
        Graph::new(ov_node_colors.len(), &ov_edges),
        ov_node_colors,
        ov_edge_colors,
    ));

    let mut por_edges = Vec::new();
    let mut por_edge_colors = Vec::new();
    let mut por_node_colors = vec![0u32; 24];
    for q in 0..4u32 {
        let base = q * 6;
        por_node_colors[base as usize] = 1;
        por_edges.extend([
            [base, base + 1],
            [base + 1, base + 2],
            [base + 2, base + 3],
            [base + 3, base + 4],
            [base + 4, base],
            [base + 4, base + 5],
            [base + 5, ((q + 1) % 4) * 6 + 1],
        ]);
        por_edge_colors.extend([2, 0, 1, 0, 2, 1, 0]);
    }
    cases.push(EdgeColoredCase::new(
        "ir",
        "porphyrin_core_3color",
        Graph::new(24, &por_edges),
        por_node_colors,
        por_edge_colors,
    ));

    cases
}

fn generated_edge_colored_cases() -> Vec<EdgeColoredCase> {
    let mut rng = ChaCha8Rng::seed_from_u64(0x2540_0001_0000_0001);
    let mut cases = Vec::new();
    for (name, n, rings, keep) in [
        ("mol_gen_12_r1", 12, 1, false),
        ("mol_gen_16_r2", 16, 2, true),
        ("mol_gen_20_r3", 20, 3, false),
        ("mol_gen_40_r4", 40, 4, false),
        ("mol_gen_80_r6", 80, 6, true),
    ] {
        let graph = prufer_tree_plus_rings(n, rings, 4, &mut rng);
        let (node_colors, edge_colors) = decorate_molecular_like(&graph, &mut rng);
        if keep {
            cases.push(EdgeColoredCase::new(
                "generated",
                name,
                graph,
                node_colors,
                edge_colors,
            ));
        }
    }
    cases
}

fn cage_edge_colored_cases() -> Vec<EdgeColoredCase> {
    let mut cases = vec![
        EdgeColoredCase::new(
            "cage",
            "cubane_2color",
            fixtures::cubane(),
            vec![0; 8],
            vec![0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
        ),
        EdgeColoredCase::new(
            "cage",
            "adamantane_2color",
            fixtures::adamantane(),
            vec![0, 1, 1, 1, 0, 0, 0, 1, 1, 1],
            vec![0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1],
        ),
    ];

    let c60 = fixtures::fullerene_c60();
    let c60_kekule: Vec<u32> = c60
        .edge_ids()
        .map(|edge| {
            let [u, v] = c60.edge_endpoints(edge);
            u32::from(!edge_in_five_cycle(&c60, u, v))
        })
        .collect();
    cases.push(EdgeColoredCase::new(
        "cage",
        "c60_kekule",
        c60,
        vec![0; 60],
        c60_kekule,
    ));

    let c70 = fixtures::fullerene_c70();
    let c70_kekule: Vec<u32> = c70
        .edge_ids()
        .map(|edge| {
            let [u, v] = c70.edge_endpoints(edge);
            u32::from(!edge_in_five_cycle(&c70, u, v))
        })
        .collect();
    cases.push(EdgeColoredCase::new(
        "cage",
        "c70_kekule",
        c70,
        vec![0; 70],
        c70_kekule,
    ));

    cases
}

fn torture_edge_colored_cases() -> Vec<EdgeColoredCase> {
    let (grid, grid_axis_colors) = grid_graph_with_axes(8, 8);
    let (q6, q6_2colors) = hypercube_graph_with_2colors(6);
    let (four_hex, hex_colors) = four_hexagons_with_kekule();

    let mut rook_edges = Vec::new();
    let mut rook_colors = Vec::new();
    for r in 0..4u32 {
        for c1 in 0..4u32 {
            for c2 in c1 + 1..4u32 {
                rook_edges.push([r * 4 + c1, r * 4 + c2]);
                rook_colors.push(0);
            }
        }
    }
    for c in 0..4u32 {
        for r1 in 0..4u32 {
            for r2 in r1 + 1..4u32 {
                rook_edges.push([r1 * 4 + c, r2 * 4 + c]);
                rook_colors.push(1);
            }
        }
    }

    let mut miy_edges = Vec::new();
    let mut miy_colors = Vec::new();
    for k in 0..5u32 {
        let b = k * 4;
        let nb = ((k + 1) % 5) * 4;
        miy_edges.extend([
            [b, b + 1],
            [b + 1, b + 2],
            [b + 2, b + 3],
            [b + 3, b],
            [b, nb + 2],
            [b + 1, nb + 3],
        ]);
        miy_colors.extend([0, 0, 0, 0, 1, 1]);
    }

    vec![
        EdgeColoredCase::new(
            "torture",
            "grid_8x8_axes",
            grid,
            vec![0; 64],
            grid_axis_colors,
        ),
        EdgeColoredCase::new("torture", "hypercube_6_2color", q6, vec![0; 64], q6_2colors),
        EdgeColoredCase::new(
            "torture",
            "four_hexagons_kekule",
            four_hex,
            vec![0; 24],
            hex_colors,
        ),
        EdgeColoredCase::new(
            "torture",
            "rook_4x4_axes",
            Graph::new(16, &rook_edges),
            vec![0; 16],
            rook_colors,
        ),
        EdgeColoredCase::new(
            "torture",
            "miyazaki_ladder_20",
            Graph::new(20, &miy_edges),
            vec![0; 20],
            miy_colors,
        ),
    ]
}

struct MolBuilder {
    node_colors: Vec<u32>,
    edges: Vec<[u32; 2]>,
    edge_colors: Vec<u32>,
}

impl MolBuilder {
    const H: u32 = 0;
    const C: u32 = 1;
    const N: u32 = 2;
    const O: u32 = 3;
    const S: u32 = 4;
    const CD: u32 = 5;
    const SE: u32 = 6;
    const TE: u32 = 7;

    const SINGLE: u32 = 0;
    const DOUBLE: u32 = 1;
    const AROMATIC: u32 = 2;
    const COORD: u32 = 3;

    fn new() -> Self {
        Self {
            node_colors: Vec::new(),
            edges: Vec::new(),
            edge_colors: Vec::new(),
        }
    }

    fn add_atom(&mut self, elem: u32) -> u32 {
        let id = self.node_colors.len() as u32;
        self.node_colors.push(elem);
        id
    }

    fn bond(&mut self, a: u32, b: u32, order: u32) {
        self.edges.push([a, b]);
        self.edge_colors.push(order);
    }

    fn add_h(&mut self, parent: u32, count: usize) {
        for _ in 0..count {
            let h = self.add_atom(Self::H);
            self.bond(parent, h, Self::SINGLE);
        }
    }

    fn add_methyl(&mut self, parent: u32) -> u32 {
        let c = self.add_atom(Self::C);
        self.bond(parent, c, Self::SINGLE);
        self.add_h(c, 3);
        c
    }

    fn add_ch2(&mut self, parent: u32) -> u32 {
        let c = self.add_atom(Self::C);
        self.bond(parent, c, Self::SINGLE);
        self.add_h(c, 2);
        c
    }

    fn add_phenyl(&mut self, parent: u32, para_oh: bool) {
        let ring: [u32; 6] = array::from_fn(|_| self.add_atom(Self::C));
        self.bond(parent, ring[0], Self::SINGLE);
        for i in 0..6 {
            self.bond(ring[i], ring[(i + 1) % 6], Self::AROMATIC);
        }
        for &c in &[ring[1], ring[2], ring[4], ring[5]] {
            self.add_h(c, 1);
        }
        if para_oh {
            let o = self.add_atom(Self::O);
            self.bond(ring[3], o, Self::SINGLE);
            self.add_h(o, 1);
        } else {
            self.add_h(ring[3], 1);
        }
    }

    fn add_isopropyl(&mut self, parent: u32) {
        let ch = self.add_atom(Self::C);
        self.bond(parent, ch, Self::SINGLE);
        self.add_h(ch, 1);
        self.add_methyl(ch);
        self.add_methyl(ch);
    }

    fn add_tert_butyl(&mut self, parent: u32) {
        let c_quat = self.add_atom(Self::C);
        self.bond(parent, c_quat, Self::SINGLE);
        self.add_methyl(c_quat);
        self.add_methyl(c_quat);
        self.add_methyl(c_quat);
    }

    fn add_residue_sidechain(&mut self, ca: u32, kind: &str) {
        match kind {
            "Gly" => self.add_h(ca, 1),
            "Ala" => {
                self.add_methyl(ca);
            }
            "Val" => self.add_isopropyl(ca),
            "Leu" => {
                let cb = self.add_ch2(ca);
                self.add_isopropyl(cb);
            }
            "Phe" => {
                let cb = self.add_ch2(ca);
                self.add_phenyl(cb, false);
            }
            "Tyr" => {
                let cb = self.add_ch2(ca);
                self.add_phenyl(cb, true);
            }
            other => panic!("unsupported residue kind {other}"),
        }
    }

    fn finish(self, name: &'static str) -> EdgeColoredCase {
        let n = self.node_colors.len();
        EdgeColoredCase::new(
            "complex",
            name,
            Graph::new(n, &self.edges),
            self.node_colors,
            self.edge_colors,
        )
    }
}

fn build_cd10_se4_core(b: &mut MolBuilder, dope_core: bool) -> Vec<u32> {
    let mut se = [0u32; 4];
    for (idx, slot) in se.iter_mut().enumerate() {
        *slot = b.add_atom(if dope_core && idx == 0 {
            MolBuilder::TE
        } else {
            MolBuilder::SE
        });
    }

    let pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
    let mut edge_cd = [0u32; 6];
    for (pidx, &(i, j)) in pairs.iter().enumerate() {
        let cd = b.add_atom(MolBuilder::CD);
        edge_cd[pidx] = cd;
        b.bond(cd, se[i], MolBuilder::COORD);
        b.bond(cd, se[j], MolBuilder::COORD);
    }

    let mut apical_cd = [0u32; 4];
    for (k, slot) in apical_cd.iter_mut().enumerate() {
        *slot = b.add_atom(if dope_core && k == 1 {
            MolBuilder::TE
        } else {
            MolBuilder::CD
        });
    }

    let mut sulfur_sites = Vec::with_capacity(16);
    for (k, &acd) in apical_cd.iter().enumerate() {
        for (pidx, &(i, j)) in pairs.iter().enumerate() {
            if i == k || j == k {
                let s = b.add_atom(MolBuilder::S);
                b.bond(acd, s, MolBuilder::COORD);
                b.bond(edge_cd[pidx], s, MolBuilder::COORD);
                sulfur_sites.push(s);
            }
        }
        let s_term = b.add_atom(MolBuilder::S);
        b.bond(acd, s_term, MolBuilder::COORD);
        sulfur_sites.push(s_term);
    }

    if dope_core {
        b.node_colors[sulfur_sites[0] as usize] = MolBuilder::O;
    }

    sulfur_sites
}

fn complex_edge_colored_cases() -> Vec<EdgeColoredCase> {
    let mut cases = Vec::new();

    {
        let mut b = MolBuilder::new();
        let ac_c = b.add_atom(MolBuilder::C);
        let ac_o = b.add_atom(MolBuilder::O);
        b.bond(ac_c, ac_o, MolBuilder::DOUBLE);
        b.add_methyl(ac_c);
        let mut prev_c = ac_c;
        for res in [
            "Phe", "Val", "Ala", "Leu", "Tyr", "Phe", "Val", "Gly", "Leu", "Phe",
        ] {
            let n = b.add_atom(MolBuilder::N);
            b.bond(prev_c, n, MolBuilder::SINGLE);
            b.add_h(n, 1);
            let ca = b.add_atom(MolBuilder::C);
            b.bond(n, ca, MolBuilder::SINGLE);
            b.add_h(ca, 1);
            b.add_residue_sidechain(ca, res);
            let c = b.add_atom(MolBuilder::C);
            let o = b.add_atom(MolBuilder::O);
            b.bond(ca, c, MolBuilder::SINGLE);
            b.bond(c, o, MolBuilder::DOUBLE);
            prev_c = c;
        }
        let term_n = b.add_atom(MolBuilder::N);
        b.bond(prev_c, term_n, MolBuilder::SINGLE);
        b.add_h(term_n, 2);
        cases.push(b.finish("decapeptide_mixed_explicit_h"));
    }

    {
        let mut b = MolBuilder::new();
        let mut first_n = 0u32;
        let mut prev_c = 0u32;
        for step in 0..10 {
            let res = if step % 2 == 0 { "Phe" } else { "Val" };
            let n = b.add_atom(MolBuilder::N);
            if step == 0 {
                first_n = n;
            } else {
                b.bond(prev_c, n, MolBuilder::SINGLE);
            }
            b.add_h(n, 1);
            let ca = b.add_atom(MolBuilder::C);
            b.bond(n, ca, MolBuilder::SINGLE);
            b.add_h(ca, 1);
            b.add_residue_sidechain(ca, res);
            let c = b.add_atom(MolBuilder::C);
            let o = b.add_atom(MolBuilder::O);
            b.bond(ca, c, MolBuilder::SINGLE);
            b.bond(c, o, MolBuilder::DOUBLE);
            prev_c = c;
        }
        b.bond(prev_c, first_n, MolBuilder::SINGLE);
        cases.push(b.finish("decapeptide_cyclic_sym"));
    }

    {
        let mut b = MolBuilder::new();
        let mut quads = [[0u32; 6]; 4];
        for quad in &mut quads {
            quad[0] = b.add_atom(MolBuilder::N);
            for slot in quad.iter_mut().skip(1) {
                *slot = b.add_atom(MolBuilder::C);
            }
        }
        for q in 0..4 {
            let [n, a1, b1, b2, a2, meso] = quads[q];
            b.bond(n, a1, MolBuilder::AROMATIC);
            b.bond(a1, b1, MolBuilder::SINGLE);
            b.bond(b1, b2, MolBuilder::DOUBLE);
            b.bond(b2, a2, MolBuilder::SINGLE);
            b.bond(a2, n, MolBuilder::AROMATIC);
            b.bond(a2, meso, MolBuilder::DOUBLE);
            b.bond(meso, quads[(q + 1) % 4][1], MolBuilder::SINGLE);
            b.add_h(meso, 1);
            b.add_phenyl(b1, false);
            b.add_phenyl(b2, false);
        }
        cases.push(b.finish("macrocycle_octaphenyl_porphyrin"));
    }

    for (name, dope_core, use_tert_butyl) in [
        ("cdse_t3_sph16", false, false),
        ("cdse_doped_sph16", true, false),
        ("cdse_t3_stbu16", false, true),
    ] {
        let mut b = MolBuilder::new();
        for s in build_cd10_se4_core(&mut b, dope_core) {
            if use_tert_butyl {
                b.add_tert_butyl(s);
            } else {
                b.add_phenyl(s, false);
            }
        }
        cases.push(b.finish(name));
    }

    cases
}

fn nodes_and_edges_corpus() -> Vec<EdgeColoredCase> {
    let mut all = Vec::new();
    all.extend(ir_edge_colored_cases());
    all.extend(generated_edge_colored_cases());
    all.extend(cage_edge_colored_cases());
    all.extend(torture_edge_colored_cases());
    all.extend(complex_edge_colored_cases());
    all
}

fn bench_nodes_only(c: &mut Criterion) {
    let cases = nodes_only_corpus();

    let mut auto = c.benchmark_group("automorphism/nodes_only/automorphisms");
    for case in &cases {
        let label = case.label();
        auto.bench_with_input(BenchmarkId::new("nauty", &label), case, |b, case| {
            b.iter(|| {
                black_box(case.graph.automorphisms(
                    |node| case.node_colors[node.index()],
                    |_| 0u8,
                    AutomorphismAlgorithm::Nauty,
                ))
            });
        });
        auto.bench_with_input(BenchmarkId::new("trail_ir", &label), case, |b, case| {
            b.iter(|| {
                black_box(case.graph.automorphisms(
                    |node| case.node_colors[node.index()],
                    |_| 0u8,
                    AutomorphismAlgorithm::TrailIr,
                ))
            });
        });
    }
    auto.finish();

    let mut canon = c.benchmark_group("automorphism/nodes_only/canonical_labeling");
    for case in &cases {
        let label = case.label();
        canon.bench_with_input(BenchmarkId::new("nauty", &label), case, |b, case| {
            b.iter(|| {
                black_box(case.graph.canonical_labeling(
                    |node| case.node_colors[node.index()],
                    |_| 0u8,
                    AutomorphismAlgorithm::Nauty,
                ))
            });
        });
        canon.bench_with_input(BenchmarkId::new("trail_ir", &label), case, |b, case| {
            b.iter(|| {
                black_box(case.graph.canonical_labeling(
                    |node| case.node_colors[node.index()],
                    |_| 0u8,
                    AutomorphismAlgorithm::TrailIr,
                ))
            });
        });
    }
    canon.finish();
}

fn bench_nodes_and_edges(c: &mut Criterion) {
    let cases = nodes_and_edges_corpus();

    let mut auto = c.benchmark_group("automorphism/nodes_and_edges/automorphisms");
    for case in &cases {
        let label = case.label();
        auto.bench_with_input(BenchmarkId::new("nauty", &label), case, |b, case| {
            b.iter(|| {
                black_box(case.graph.automorphisms(
                    |node| case.node_colors[node.index()],
                    |edge| case.edge_colors[edge.index()],
                    AutomorphismAlgorithm::Nauty,
                ))
            });
        });
        auto.bench_with_input(BenchmarkId::new("trail_ir", &label), case, |b, case| {
            b.iter(|| {
                black_box(case.graph.automorphisms(
                    |node| case.node_colors[node.index()],
                    |edge| case.edge_colors[edge.index()],
                    AutomorphismAlgorithm::TrailIr,
                ))
            });
        });
    }
    auto.finish();

    let mut canon = c.benchmark_group("automorphism/nodes_and_edges/canonical_labeling");
    for case in &cases {
        let label = case.label();
        canon.bench_with_input(BenchmarkId::new("nauty", &label), case, |b, case| {
            b.iter(|| {
                black_box(case.graph.canonical_labeling(
                    |node| case.node_colors[node.index()],
                    |edge| case.edge_colors[edge.index()],
                    AutomorphismAlgorithm::Nauty,
                ))
            });
        });
        canon.bench_with_input(BenchmarkId::new("trail_ir", &label), case, |b, case| {
            b.iter(|| {
                black_box(case.graph.canonical_labeling(
                    |node| case.node_colors[node.index()],
                    |edge| case.edge_colors[edge.index()],
                    AutomorphismAlgorithm::TrailIr,
                ))
            });
        });
    }
    canon.finish();
}

fn bench_stabilizers(c: &mut Criterion) {
    let (four_hexagons, _) = four_hexagons_with_kekule();
    let node_cases: Vec<(&str, Graph, Vec<u32>)> = [
        ("adamantane", fixtures::adamantane()),
        ("c60", fixtures::fullerene_c60()),
        ("four_hexagons", four_hexagons),
    ]
    .into_iter()
    .map(|(name, graph)| {
        let colors = degree_colors(&graph);
        (name, graph, colors)
    })
    .collect();

    let mut nodes_group = c.benchmark_group("automorphism/stabilizers/nodes_only");
    for (name, graph, colors) in &node_cases {
        let n = graph.node_count();
        let sites: [NodeId; 3] = [0, n / 3, n / 2].map(|s| NodeId(s as u32));
        let unique_site_color = colors.iter().copied().max().unwrap_or(0) + 1;
        let label = format!("{name}/three_sites");

        nodes_group.bench_function(BenchmarkId::new("nauty", &label), |b| {
            b.iter(|| {
                for &site in &sites {
                    black_box(graph.automorphisms(
                        |node| {
                            if node == site {
                                unique_site_color
                            } else {
                                colors[node.index()]
                            }
                        },
                        |_| 0u8,
                        AutomorphismAlgorithm::Nauty,
                    ));
                }
            });
        });
        nodes_group.bench_function(BenchmarkId::new("trail_ir", &label), |b| {
            b.iter(|| {
                for &site in &sites {
                    black_box(graph.automorphisms(
                        |node| {
                            if node == site {
                                unique_site_color
                            } else {
                                colors[node.index()]
                            }
                        },
                        |_| 0u8,
                        AutomorphismAlgorithm::TrailIr,
                    ));
                }
            });
        });
        nodes_group.bench_function(BenchmarkId::new("trail_ir_scratch", &label), |b| {
            let mut scratch = AutomorphismScratch::new();
            b.iter(|| {
                for &site in &sites {
                    black_box(graph.automorphisms_with_scratch(
                        |node| {
                            if node == site {
                                unique_site_color
                            } else {
                                colors[node.index()]
                            }
                        },
                        |_| 0u8,
                        AutomorphismAlgorithm::TrailIr,
                        &mut scratch,
                    ));
                }
            });
        });
    }
    nodes_group.finish();

    const EDGE_TARGETS: [&str; 4] = [
        "dimethylbutanol_stereo_levi",
        "cubane_2color",
        "c60_kekule",
        "cdse_t3_sph16",
    ];
    let edge_cases: Vec<EdgeColoredCase> = nodes_and_edges_corpus()
        .into_iter()
        .filter(|case| EDGE_TARGETS.contains(&case.name))
        .collect();

    let mut edges_group = c.benchmark_group("automorphism/stabilizers/nodes_and_edges");
    for case in &edge_cases {
        let n = case.graph.node_count();
        let sites: [NodeId; 3] = [0, n / 3, n / 2].map(|s| NodeId(s as u32));
        let unique_site_color = case.node_colors.iter().copied().max().unwrap_or(0) + 1;
        let label = format!("{}/three_sites", case.name);

        edges_group.bench_with_input(BenchmarkId::new("nauty", &label), case, |b, case| {
            b.iter(|| {
                for &site in &sites {
                    black_box(case.graph.automorphisms(
                        |node| {
                            if node == site {
                                unique_site_color
                            } else {
                                case.node_colors[node.index()]
                            }
                        },
                        |edge| case.edge_colors[edge.index()],
                        AutomorphismAlgorithm::Nauty,
                    ));
                }
            });
        });
        edges_group.bench_with_input(BenchmarkId::new("trail_ir", &label), case, |b, case| {
            b.iter(|| {
                for &site in &sites {
                    black_box(case.graph.automorphisms(
                        |node| {
                            if node == site {
                                unique_site_color
                            } else {
                                case.node_colors[node.index()]
                            }
                        },
                        |edge| case.edge_colors[edge.index()],
                        AutomorphismAlgorithm::TrailIr,
                    ));
                }
            });
        });
        edges_group.bench_with_input(
            BenchmarkId::new("trail_ir_scratch", &label),
            case,
            |b, case| {
                let mut scratch = AutomorphismScratch::new();
                b.iter(|| {
                    for &site in &sites {
                        black_box(case.graph.automorphisms_with_scratch(
                            |node| {
                                if node == site {
                                    unique_site_color
                                } else {
                                    case.node_colors[node.index()]
                                }
                            },
                            |edge| case.edge_colors[edge.index()],
                            AutomorphismAlgorithm::TrailIr,
                            &mut scratch,
                        ));
                    }
                });
            },
        );
    }
    edges_group.finish();
}

fn bench_canonical_key(c: &mut Criterion) {
    let (grid_8x8, _) = grid_graph_with_axes(8, 8);
    let cases = [
        ("naphthalene", fixtures::naphthalene()),
        ("adamantane", fixtures::adamantane()),
        ("c60", fixtures::fullerene_c60()),
        ("grid_8x8", grid_8x8),
    ];

    let mut group = c.benchmark_group("automorphism/canonical_key");
    for (name, graph) in &cases {
        let node_colors: Vec<Vec<u8>> = degree_colors(graph)
            .into_iter()
            .map(|color| color.to_le_bytes().to_vec())
            .collect();
        let edge_colors: Vec<Vec<u8>> = graph
            .edge_ids()
            .map(|edge| vec![(edge.index() % 3) as u8])
            .collect();
        group.bench_function(BenchmarkId::new("nauty", name), |b| {
            b.iter(|| {
                black_box(graph.canonical_key(
                    |node| node_colors[node.index()].clone(),
                    |edge| edge_colors[edge.index()].clone(),
                    AutomorphismAlgorithm::Nauty,
                ))
            });
        });
        group.bench_function(BenchmarkId::new("trail_ir", name), |b| {
            b.iter(|| {
                black_box(graph.canonical_key(
                    |node| node_colors[node.index()].clone(),
                    |edge| edge_colors[edge.index()].clone(),
                    AutomorphismAlgorithm::TrailIr,
                ))
            });
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_nodes_only,
    bench_nodes_and_edges,
    bench_stabilizers,
    bench_canonical_key,
);
criterion_main!(benches);
