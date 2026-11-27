use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};

use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use rstest::rstest;
use serde::Deserialize;
use umol_nauty_sys::{run_weighted, WeightedNautyInput};

use super::components::{find_non_uniform_components, NucWorkspace};
use super::partition::{EdgeColoredView, Partition};
use super::refine::{refine, RefineWorkspace};
use super::search::NUCR_MIN_CELLS;
use super::TrailIrWorkspace;
use crate::algorithms::automorphism::{AutomorphismGroupOrder, AutomorphismOutput};
use crate::graph::Graph;
use crate::id::NodeId;
use crate::union_find::UnionFind;

const SIMPLE_GRAPHS_GRAPH6: &str = include_str!("../../../../tests/data/simple-through-8.g6");
const FULLERENE_C60_TOML: &str =
    include_str!("../../../../tests/matching/data/fullerene_c60_planar.toml");
const SUBGROUP_CLOSURE_LIMIT: u128 = 4_096;

#[derive(Deserialize)]
struct GraphFixture {
    node_count: usize,
    edges: Vec<[u32; 2]>,
}

#[derive(Clone)]
struct ColoredGraphCase {
    name: &'static str,
    graph: Graph,
    node_colors: Vec<u32>,
    edge_colors: Vec<u32>,
}

fn parse_graph6(source: &str) -> Graph {
    let bytes = source.as_bytes();
    let node_count = (bytes[0] - 63) as usize;
    let mut edges = Vec::new();
    let mut bit = 0usize;
    for second in 1..node_count {
        for first in 0..second {
            let value = bytes[1 + bit / 6] - 63;
            if value & (1u8 << (5 - bit % 6)) != 0 {
                edges.push([first as u32, second as u32]);
            }
            bit += 1;
        }
    }
    Graph::new(node_count, &edges)
}

fn simple_graphs() -> impl Iterator<Item = Graph> {
    SIMPLE_GRAPHS_GRAPH6.lines().map(parse_graph6)
}

fn solve_weighted_nauty(
    graph: &Graph,
    node_colors: &[u32],
    edge_colors: &[u32],
    get_canon: bool,
) -> AutomorphismOutput {
    let node_count = graph.node_count();
    let mut indexed: Vec<(usize, u32)> = graph
        .node_ids()
        .map(|node| (node.index(), node_colors[node.index()]))
        .collect();
    indexed.sort_by_key(|&(_, color)| color);

    let mut colors = vec![0u32; node_count];
    let mut rank = 0u32;
    for (position, &(vertex, color)) in indexed.iter().enumerate() {
        if position > 0 && color != indexed[position - 1].1 {
            rank += 1;
        }
        colors[vertex] = rank;
    }
    let partition = indexed.iter().map(|&(vertex, _)| vertex as u32).collect();

    let mut offsets = Vec::with_capacity(node_count + 1);
    let mut neighbors = Vec::with_capacity(2 * graph.edge_count());
    let mut weights = Vec::with_capacity(2 * graph.edge_count());
    offsets.push(0);
    for node in graph.node_ids() {
        for neighbor in graph.neighbors(node) {
            neighbors.push(neighbor.node.0);
            weights.push(edge_colors[neighbor.edge.index()]);
        }
        offsets.push(neighbors.len());
    }

    let input =
        WeightedNautyInput::try_new(node_count, offsets, neighbors, weights, colors, partition)
            .expect("simple graph produces valid weighted nauty input");
    let output = run_weighted(&input, get_canon).expect("weighted nauty succeeds");
    let orbits: Vec<NodeId> = output.orbits.into_iter().map(NodeId).collect();
    let orbit_count = orbits.iter().copied().collect::<HashSet<_>>().len();
    AutomorphismOutput {
        orbits,
        canonical_labels: output.canonical_labels.into_iter().map(NodeId).collect(),
        node_count,
        orbit_count,
        group_order: AutomorphismGroupOrder::from_scientific(
            output.group_order.mantissa,
            output.group_order.exponent,
        ),
        generators: output
            .generators
            .into_iter()
            .map(|generator| generator.into_iter().map(NodeId).collect())
            .collect(),
    }
}

fn edge_colored_canonical_key(
    graph: &Graph,
    node_colors: &[u32],
    edge_colors: &[u32],
    canonical_labels: &[NodeId],
) -> Vec<u32> {
    let n = graph.node_count();
    assert_eq!(canonical_labels.len(), n);

    let mut inv = vec![0u32; n];
    for (rank, &vertex) in canonical_labels.iter().enumerate() {
        inv[vertex.index()] = rank as u32;
    }

    let mut key = Vec::with_capacity(n + 3 * graph.edge_count() + 2);
    key.push(n as u32);
    for &vertex in canonical_labels {
        key.push(node_colors[vertex.index()]);
    }

    let mut canonical_edges: Vec<(u32, u32, u32)> = graph
        .edge_ids()
        .map(|edge| {
            let [a, b] = graph.edge_endpoints(edge);
            let ra = inv[a.index()];
            let rb = inv[b.index()];
            (ra.min(rb), ra.max(rb), edge_colors[edge.index()])
        })
        .collect();
    canonical_edges.sort_unstable();

    key.push(canonical_edges.len() as u32);
    for (u, v, color) in canonical_edges {
        key.push(u);
        key.push(v);
        key.push(color);
    }
    key
}

fn permute_colored_graph(
    graph: &Graph,
    node_colors: &[u32],
    edge_colors: &[u32],
    rng: &mut impl Rng,
) -> (Graph, Vec<u32>, Vec<u32>) {
    let n = graph.node_count();
    let mut old_to_new: Vec<u32> = (0..n as u32).collect();
    old_to_new.shuffle(rng);

    let mut new_node_colors = vec![0u32; n];
    for old in 0..n {
        new_node_colors[old_to_new[old] as usize] = node_colors[old];
    }

    let mut raw_edges: Vec<([u32; 2], u32)> = graph
        .edge_ids()
        .map(|edge| {
            let [a, b] = graph.edge_endpoints(edge);
            (
                [old_to_new[a.index()], old_to_new[b.index()]],
                edge_colors[edge.index()],
            )
        })
        .collect();
    raw_edges.shuffle(rng);

    let edges: Vec<[u32; 2]> = raw_edges.iter().map(|&(pair, _)| pair).collect();
    let new_edge_colors: Vec<u32> = raw_edges.iter().map(|&(_, color)| color).collect();
    (Graph::new(n, &edges), new_node_colors, new_edge_colors)
}

fn generated_subgroup(node_count: usize, generators: &[Vec<NodeId>]) -> HashSet<Vec<u32>> {
    let identity: Vec<u32> = (0..node_count as u32).collect();
    let mut group = HashSet::new();
    let mut queue = VecDeque::new();
    group.insert(identity.clone());
    queue.push_back(identity);

    while let Some(current) = queue.pop_front() {
        for generator in generators {
            let next: Vec<u32> = current.iter().map(|&v| generator[v as usize].0).collect();
            if group.insert(next.clone()) {
                queue.push_back(next);
            }
        }
    }
    group
}

fn assert_output_semantics(
    label: &str,
    graph: &Graph,
    node_colors: &[u32],
    edge_colors: &[u32],
    get_canon: bool,
    output: &AutomorphismOutput,
) {
    let n = graph.node_count();
    assert_eq!(output.node_count(), n, "node_count mismatch on {label}");
    if get_canon {
        let expected_nodes: Vec<NodeId> = (0..n as u32).map(NodeId).collect();
        let mut sorted_labels = output.canonical_labels().to_vec();
        sorted_labels.sort_unstable();
        assert_eq!(
            sorted_labels, expected_nodes,
            "canonical_labels is not a permutation on {label}"
        );
    } else {
        assert!(
            output.canonical_labels().is_empty(),
            "canonical_labels must be empty when get_canon = false on {label}"
        );
    }

    let expected_nodes: Vec<NodeId> = (0..n as u32).map(NodeId).collect();
    let edge_map: HashMap<(usize, usize), u32> = graph
        .edge_ids()
        .map(|edge| {
            let [a, b] = graph.edge_endpoints(edge);
            (
                (a.index().min(b.index()), a.index().max(b.index())),
                edge_colors[edge.index()],
            )
        })
        .collect();

    let mut generated_orbits = UnionFind::new(n);
    for generator in output.generators() {
        assert_eq!(generator.len(), n);
        let mut sorted = generator.clone();
        sorted.sort_unstable();
        assert_eq!(
            sorted, expected_nodes,
            "generator not a bijection on {label}"
        );

        for (src, &dst) in generator.iter().enumerate() {
            assert_eq!(
                node_colors[src],
                node_colors[dst.index()],
                "generator violated node color on {label}"
            );
            generated_orbits.union(src, dst.index());
        }
        for (&(u, v), &color) in &edge_map {
            let iu = generator[u].index();
            let iv = generator[v].index();
            assert_eq!(
                edge_map.get(&(iu.min(iv), iu.max(iv))),
                Some(&color),
                "generator violated edge color on {label}"
            );
        }
    }

    let mut min_rep = vec![usize::MAX; n];
    for v in 0..n {
        let root = generated_orbits.find(v);
        min_rep[root] = min_rep[root].min(v);
    }

    let mut reps = HashSet::new();
    for a in 0..n {
        let rep = output.orbit_of(NodeId(a as u32));
        reps.insert(rep);
        assert_eq!(
            rep.index(),
            min_rep[generated_orbits.find(a)],
            "orbit representative mismatch on {label}"
        );
        for b in a..n {
            assert_eq!(
                generated_orbits.find(a) == generated_orbits.find(b),
                output.same_orbit(NodeId(a as u32), NodeId(b as u32)),
                "orbit membership mismatch on {label}"
            );
        }
    }
    assert_eq!(
        output.orbit_count(),
        reps.len(),
        "orbit_count mismatch on {label}"
    );
}

fn assert_group_orders_agree(
    label: &str,
    left: AutomorphismGroupOrder,
    right: AutomorphismGroupOrder,
) {
    match (left, right) {
        (
            AutomorphismGroupOrder::Scientific {
                mantissa: lm,
                exponent: le,
            },
            AutomorphismGroupOrder::Scientific {
                mantissa: rm,
                exponent: re,
            },
        ) => {
            let normalize = |mut m: f64, mut e: i32| {
                while m >= 10.0 {
                    m /= 10.0;
                    e += 1;
                }
                while m > 0.0 && m < 1.0 {
                    m *= 10.0;
                    e -= 1;
                }
                (m, e)
            };
            let (lm, le) = normalize(lm, le);
            let (rm, re) = normalize(rm, re);
            assert_eq!(
                le, re,
                "scientific group_order exponent differs on {label}: {left:?} vs {right:?}"
            );
            assert!(
                (lm - rm).abs() <= 1e-12 * lm.abs().max(1.0),
                "scientific group_order mantissa differs on {label}: {left:?} vs {right:?}"
            );
        }
        _ => assert_eq!(left, right, "group_order differs on {label}"),
    }
}

fn assert_outputs_agree(
    label: &str,
    graph: &Graph,
    node_colors: &[u32],
    edge_colors: &[u32],
    get_canon: bool,
    nauty: &AutomorphismOutput,
    trail_ir: &AutomorphismOutput,
) {
    assert_output_semantics(label, graph, node_colors, edge_colors, get_canon, nauty);
    assert_output_semantics(label, graph, node_colors, edge_colors, get_canon, trail_ir);
    assert_eq!(
        trail_ir.orbit_count(),
        nauty.orbit_count(),
        "orbit_count differs on {label}"
    );
    for v in 0..graph.node_count() {
        assert_eq!(
            trail_ir.orbit_of(NodeId(v as u32)),
            nauty.orbit_of(NodeId(v as u32)),
            "orbit_of({v}) differs on {label}"
        );
    }
    assert_group_orders_agree(label, nauty.group_order(), trail_ir.group_order());

    if let Some(exact) = nauty.group_order().exact_value() {
        if exact <= SUBGROUP_CLOSURE_LIMIT {
            let nauty_group = generated_subgroup(graph.node_count(), nauty.generators());
            let trail_group = generated_subgroup(graph.node_count(), trail_ir.generators());
            assert_eq!(
                trail_group.len() as u128,
                exact,
                "generated subgroup order differs from group_order on {label}"
            );
            assert_eq!(
                trail_group, nauty_group,
                "generated subgroup differs from weighted nauty on {label}"
            );
        }
    }
}

fn partition_cells_sorted(part: &Partition) -> Vec<Vec<NodeId>> {
    let n = part.lab.len();
    let mut out = Vec::with_capacity(part.num_cells);
    let mut i = 0;
    while i < n {
        let len = part.cell_len[i] as usize;
        let mut cell: Vec<NodeId> = part.lab[i..i + len].iter().map(|&v| NodeId(v)).collect();
        cell.sort_unstable();
        out.push(cell);
        i += len;
    }
    out
}

fn assert_partition_is_edge_color_equitable(
    label: &str,
    graph: &Graph,
    edge_colors: &[u32],
    part: &Partition,
) {
    let n = graph.node_count();
    let mut i = 0;
    while i < n {
        let len = part.cell_len[i] as usize;
        if len > 1 {
            let signature = |v: u32| {
                let mut counts: BTreeMap<(u32, u32), usize> = BTreeMap::new();
                for nb in graph.neighbors(NodeId(v)) {
                    let cs = part.cell_start[nb.node.index()];
                    let wt = edge_colors[nb.edge.index()];
                    *counts.entry((cs, wt)).or_insert(0) += 1;
                }
                counts
            };
            let expected = signature(part.lab[i]);
            for &u in &part.lab[i + 1..i + len] {
                assert_eq!(
                    signature(u),
                    expected,
                    "partition not edge-color-equitable on {label}"
                );
            }
        }
        i += len;
    }
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

fn representative_fixtures() -> Vec<ColoredGraphCase> {
    let mut cases = vec![
        ColoredGraphCase {
            name: "empty",
            graph: Graph::default(),
            node_colors: vec![],
            edge_colors: vec![],
        },
        ColoredGraphCase {
            name: "singleton",
            graph: Graph::new(1, &[]),
            node_colors: vec![0],
            edge_colors: vec![],
        },
        ColoredGraphCase {
            name: "propane_uniform",
            graph: Graph::new(
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
            node_colors: vec![0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1],
            edge_colors: vec![0; 10],
        },
        ColoredGraphCase {
            name: "ethane_radicals",
            graph: Graph::new(8, &[[6, 7]]),
            node_colors: vec![0, 0, 0, 0, 0, 0, 1, 1],
            edge_colors: vec![0],
        },
        ColoredGraphCase {
            name: "naphthalene_kekule",
            graph: Graph::new(
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
            node_colors: vec![0; 10],
            edge_colors: vec![1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0],
        },
        ColoredGraphCase {
            name: "disconnected_rings_kekule",
            graph: Graph::new(
                12,
                &[
                    [0, 1],
                    [1, 2],
                    [2, 3],
                    [3, 4],
                    [4, 5],
                    [5, 0],
                    [6, 7],
                    [7, 8],
                    [8, 9],
                    [9, 10],
                    [10, 11],
                    [11, 6],
                ],
            ),
            node_colors: vec![0; 12],
            edge_colors: vec![1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
        },
        ColoredGraphCase {
            name: "cubane_2color",
            graph: Graph::new(
                8,
                &[
                    [0, 1],
                    [1, 2],
                    [2, 3],
                    [3, 0],
                    [4, 5],
                    [5, 6],
                    [6, 7],
                    [7, 4],
                    [0, 4],
                    [1, 5],
                    [2, 6],
                    [3, 7],
                ],
            ),
            node_colors: vec![0; 8],
            edge_colors: vec![0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
        },
        ColoredGraphCase {
            name: "adamantane_2color",
            graph: Graph::new(
                10,
                &[
                    [0, 1],
                    [0, 2],
                    [0, 3],
                    [1, 4],
                    [1, 5],
                    [2, 4],
                    [2, 6],
                    [3, 5],
                    [3, 6],
                    [4, 7],
                    [5, 8],
                    [6, 9],
                    [7, 8],
                    [7, 9],
                    [8, 9],
                ],
            ),
            node_colors: vec![0, 1, 1, 1, 0, 0, 0, 1, 1, 1],
            edge_colors: vec![0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1],
        },
        ColoredGraphCase {
            name: "petersen_2color",
            graph: Graph::new(
                10,
                &[
                    [0, 1],
                    [1, 2],
                    [2, 3],
                    [3, 4],
                    [4, 0],
                    [5, 7],
                    [7, 9],
                    [9, 6],
                    [6, 8],
                    [8, 5],
                    [0, 5],
                    [1, 6],
                    [2, 7],
                    [3, 8],
                    [4, 9],
                ],
            ),
            node_colors: vec![0; 10],
            edge_colors: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
        },
        ColoredGraphCase {
            name: "dodecahedron_2color",
            graph: Graph::new(
                20,
                &[
                    [0, 1],
                    [1, 2],
                    [2, 3],
                    [3, 4],
                    [4, 0],
                    [0, 5],
                    [1, 6],
                    [2, 7],
                    [3, 8],
                    [4, 9],
                    [5, 10],
                    [5, 14],
                    [6, 10],
                    [6, 11],
                    [7, 11],
                    [7, 12],
                    [8, 12],
                    [8, 13],
                    [9, 13],
                    [9, 14],
                    [10, 15],
                    [11, 16],
                    [12, 17],
                    [13, 18],
                    [14, 19],
                    [15, 16],
                    [16, 17],
                    [17, 18],
                    [18, 19],
                    [19, 15],
                ],
            ),
            node_colors: vec![0; 20],
            edge_colors: (0..30)
                .map(|i| u32::from((5..10).contains(&i) || (20..25).contains(&i)))
                .collect(),
        },
    ];

    // Porphyrin core (24 atoms, 28 bonds, 3 bond colors).
    let mut por_edges = Vec::new();
    let mut por_edge_colors = Vec::new();
    let mut por_node_colors = vec![0u32; 24];
    for q in 0..4u32 {
        let base = q * 6;
        por_node_colors[base as usize] = 1;
        por_edges.push([base, base + 1]);
        por_edge_colors.push(2);
        por_edges.push([base + 1, base + 2]);
        por_edge_colors.push(0);
        por_edges.push([base + 2, base + 3]);
        por_edge_colors.push(1);
        por_edges.push([base + 3, base + 4]);
        por_edge_colors.push(0);
        por_edges.push([base + 4, base]);
        por_edge_colors.push(2);
        por_edges.push([base + 4, base + 5]);
        por_edge_colors.push(1);
        por_edges.push([base + 5, ((q + 1) % 4) * 6 + 1]);
        por_edge_colors.push(0);
    }
    cases.push(ColoredGraphCase {
        name: "porphyrin_core_3color",
        graph: Graph::new(24, &por_edges),
        node_colors: por_node_colors,
        edge_colors: por_edge_colors,
    });

    // Fullerene C60 Kekulé and 3-color variants.
    let fixture: GraphFixture =
        toml::from_str(FULLERENE_C60_TOML).expect("C60 fixture is valid TOML");
    let c60 = Graph::new(fixture.node_count, &fixture.edges);
    let c60_kekule_colors: Vec<u32> = c60
        .edge_ids()
        .map(|edge| {
            let [u, v] = c60.edge_endpoints(edge);
            u32::from(!edge_in_five_cycle(&c60, u, v))
        })
        .collect();
    let c60_3colors: Vec<u32> = (0..c60.edge_count()).map(|i| (i % 3) as u32).collect();
    cases.push(ColoredGraphCase {
        name: "c60_kekule",
        graph: c60.clone(),
        node_colors: vec![0; 60],
        edge_colors: c60_kekule_colors,
    });
    cases.push(ColoredGraphCase {
        name: "c60_3color",
        graph: c60,
        node_colors: vec![0; 60],
        edge_colors: c60_3colors,
    });

    // 8x8 grid with axis-colored and 3-colored edges.
    let mut grid_edges = Vec::new();
    let mut grid_axis_colors = Vec::new();
    for r in 0..8u32 {
        for c in 0..8u32 {
            let u = r * 8 + c;
            if c + 1 < 8 {
                grid_edges.push([u, u + 1]);
                grid_axis_colors.push(0);
            }
            if r + 1 < 8 {
                grid_edges.push([u, u + 8]);
                grid_axis_colors.push(1);
            }
        }
    }
    let grid = Graph::new(64, &grid_edges);
    let grid_3colors: Vec<u32> = (0..grid.edge_count()).map(|i| (i % 3) as u32).collect();
    cases.push(ColoredGraphCase {
        name: "grid_8x8_axes",
        graph: grid.clone(),
        node_colors: vec![0; 64],
        edge_colors: grid_axis_colors,
    });
    cases.push(ColoredGraphCase {
        name: "grid_8x8_3color",
        graph: grid,
        node_colors: vec![0; 64],
        edge_colors: grid_3colors,
    });

    // K10 with a 10-cycle distinguished in color 1.
    let mut k10_edges = Vec::new();
    let mut k10_colors = Vec::new();
    for u in 0..10u32 {
        for v in u + 1..10u32 {
            k10_edges.push([u, v]);
            let on_cycle = v == u + 1 || (u == 0 && v == 9);
            k10_colors.push(u32::from(on_cycle));
        }
    }
    cases.push(ColoredGraphCase {
        name: "k10_2factor",
        graph: Graph::new(10, &k10_edges),
        node_colors: vec![0; 10],
        edge_colors: k10_colors,
    });

    // Neopentane C(CH3)4 (false twins, |G| = 31,104).
    let mut neo_edges: Vec<[u32; 2]> = vec![[0, 1], [0, 2], [0, 3], [0, 4]];
    for methyl in 0..4u32 {
        let c = 1 + methyl;
        for h in 0..3u32 {
            neo_edges.push([c, 5 + methyl * 3 + h]);
        }
    }
    let mut neo_colors = vec![0u32; 5];
    neo_colors.resize(17, 1);
    cases.push(ColoredGraphCase {
        name: "neopentane",
        graph: Graph::new(17, &neo_edges),
        node_colors: neo_colors,
        edge_colors: vec![0; 16],
    });

    // Monochromatic true-twin K4 joined to hub + 3 isolated vertices (|G| = 144).
    let mut twin_edges: Vec<[u32; 2]> = vec![[0, 1], [0, 2], [0, 3], [0, 4]];
    for i in 1..=4u32 {
        for j in i + 1..=4u32 {
            twin_edges.push([i, j]);
        }
    }
    cases.push(ColoredGraphCase {
        name: "hub_k4_plus_isolated",
        graph: Graph::new(8, &twin_edges),
        node_colors: vec![0, 1, 1, 1, 1, 2, 2, 2],
        edge_colors: vec![0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
    });

    // Non-monochromatic K4 (must not collapse as a monochromatic true-twin clique).
    cases.push(ColoredGraphCase {
        name: "two_color_k4",
        graph: Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0], [0, 2], [1, 3]]),
        node_colors: vec![0; 4],
        edge_colors: vec![0, 0, 0, 0, 1, 1],
    });

    cases
}

fn nucr_complex_fixtures() -> Vec<ColoredGraphCase> {
    // 1) Nine disjoint alternating 4-cycles with distinct vertex colors (nnt = 9 >= NUCR_MIN_CELLS).
    let mut ring_edges = Vec::new();
    let mut ring_edge_colors = Vec::new();
    let mut ring_node_colors = Vec::new();
    for comp in 0..9u32 {
        let base = comp * 4;
        ring_node_colors.extend([comp; 4]);
        ring_edges.push([base, base + 1]);
        ring_edge_colors.push(1);
        ring_edges.push([base + 1, base + 2]);
        ring_edge_colors.push(0);
        ring_edges.push([base + 2, base + 3]);
        ring_edge_colors.push(1);
        ring_edges.push([base + 3, base]);
        ring_edge_colors.push(0);
    }

    // 2) Multi-ligand hub complex: central hub vertex 0 connected to 5 ligand pairs, each pair
    // consisting of two 2-color branches with two distinct non-singleton cells per pair (nnt = 10).
    let mut hub_edges = Vec::new();
    let mut hub_edge_colors = Vec::new();
    let mut hub_node_colors = vec![0u32];
    for pair in 0..5u32 {
        let alpha_color = 1 + 2 * pair;
        let beta_color = 2 + 2 * pair;
        let a0 = hub_node_colors.len() as u32;
        hub_node_colors.push(alpha_color);
        let a1 = hub_node_colors.len() as u32;
        hub_node_colors.push(alpha_color);
        let b0 = hub_node_colors.len() as u32;
        hub_node_colors.push(beta_color);
        let b1 = hub_node_colors.len() as u32;
        hub_node_colors.push(beta_color);
        let b2 = hub_node_colors.len() as u32;
        hub_node_colors.push(beta_color);
        let b3 = hub_node_colors.len() as u32;
        hub_node_colors.push(beta_color);

        hub_edges.push([0, a0]);
        hub_edge_colors.push(pair % 3);
        hub_edges.push([0, a1]);
        hub_edge_colors.push(pair % 3);
        hub_edges.push([a0, b0]);
        hub_edge_colors.push(0);
        hub_edges.push([a0, b1]);
        hub_edge_colors.push(1);
        hub_edges.push([a1, b2]);
        hub_edge_colors.push(0);
        hub_edges.push([a1, b3]);
        hub_edge_colors.push(1);
        hub_edges.push([b0, b1]);
        hub_edge_colors.push(2);
        hub_edges.push([b2, b3]);
        hub_edge_colors.push(2);
    }

    // 3) Symmetric two-half scaffold where NUCR activates at depth > 0 after breaking the top-level
    // symmetry between two identical multi-component halves.
    let mut deep_edges = vec![[0, 1]];
    let mut deep_edge_colors = vec![2u32];
    let mut deep_node_colors = vec![0u32, 0u32];
    for half in 0..2u32 {
        for branch in 0..5u32 {
            let base = deep_node_colors.len() as u32;
            let c_alpha = 1 + 2 * branch;
            let c_beta = 2 + 2 * branch;
            deep_node_colors.extend([c_alpha, c_alpha, c_beta, c_beta]);
            deep_edges.push([half, base]);
            deep_edge_colors.push(0);
            deep_edges.push([half, base + 1]);
            deep_edge_colors.push(0);
            deep_edges.push([base, base + 2]);
            deep_edge_colors.push(0);
            deep_edges.push([base + 1, base + 3]);
            deep_edge_colors.push(0);
            deep_edges.push([base + 2, base + 3]);
            deep_edge_colors.push(1);
        }
    }

    vec![
        ColoredGraphCase {
            name: "nucr_disjoint_colored_rings",
            graph: Graph::new(ring_node_colors.len(), &ring_edges),
            node_colors: ring_node_colors,
            edge_colors: ring_edge_colors,
        },
        ColoredGraphCase {
            name: "nucr_multi_ligand_hub",
            graph: Graph::new(hub_node_colors.len(), &hub_edges),
            node_colors: hub_node_colors,
            edge_colors: hub_edge_colors,
        },
        ColoredGraphCase {
            name: "nucr_deep_two_half_scaffold",
            graph: Graph::new(deep_node_colors.len(), &deep_edges),
            node_colors: deep_node_colors,
            edge_colors: deep_edge_colors,
        },
    ]
}

#[rstest]
fn test_partition_individualize_and_unwind() {
    let mut part = Partition::default();
    let mut initial = Vec::new();
    part.reset_from_colors(&[2, 0, 2, 1, 0, 2], &mut initial);

    assert_eq!(part.num_cells, 3);
    assert_eq!(initial, vec![0, 2, 3]);
    let root_cells = partition_cells_sorted(&part);
    assert_eq!(
        root_cells,
        vec![
            vec![NodeId(1), NodeId(4)],
            vec![NodeId(3)],
            vec![NodeId(0), NodeId(2), NodeId(5)],
        ]
    );

    let m0 = part.trail_mark();
    part.individualize(3, 5);
    assert_eq!(part.num_cells, 4);
    assert_eq!(part.cell_slice(3), &[5]);
    assert_eq!(part.cell_len[4], 2);

    let m1 = part.trail_mark();
    part.individualize(4, 2);
    assert_eq!(part.num_cells, 5);
    assert_eq!(part.cell_slice(4), &[2]);
    assert_eq!(part.cell_slice(5), &[0]);

    part.unwind(m1);
    assert_eq!(part.num_cells, 4);
    assert_eq!(part.cell_slice(3), &[5]);
    assert_eq!(part.cell_len[4], 2);
    for &u in part.cell_slice(4) {
        assert_eq!(part.cell_start[u as usize], 4);
    }

    part.unwind(m0);
    assert_eq!(part.num_cells, 3);
    assert_eq!(partition_cells_sorted(&part), root_cells);
    for (cell_start, cell) in [(0u32, &[1u32, 4][..]), (2, &[3]), (3, &[0, 2, 5])] {
        for &u in cell {
            assert_eq!(part.cell_start[u as usize], cell_start);
            assert_eq!(part.lab[part.inv_lab[u as usize] as usize], u);
        }
    }
}

#[rstest]
fn test_refine_equitability_and_code_invariance() {
    let mut rng = ChaCha8Rng::seed_from_u64(0x2540_0005_0000_0001);
    let mut part = Partition::default();
    let mut ws = RefineWorkspace::default();
    let mut initial = Vec::new();

    for case in representative_fixtures() {
        let view = EdgeColoredView::new(&case.graph, &case.edge_colors);
        part.reset_from_colors(&case.node_colors, &mut initial);
        let code = refine(view, &mut part, &initial, &mut ws);
        assert_partition_is_edge_color_equitable(case.name, &case.graph, &case.edge_colors, &part);

        let (pgraph, pnodes, pedges) =
            permute_colored_graph(&case.graph, &case.node_colors, &case.edge_colors, &mut rng);
        let pview = EdgeColoredView::new(&pgraph, &pedges);
        part.reset_from_colors(&pnodes, &mut initial);
        let pcode = refine(pview, &mut part, &initial, &mut ws);
        assert_partition_is_edge_color_equitable(case.name, &pgraph, &pedges, &part);
        assert_eq!(
            pcode, code,
            "refinement code not relabeling-invariant on {}",
            case.name
        );
    }

    let mut checked = 0usize;
    for (idx, graph) in simple_graphs().enumerate() {
        if graph.node_count() > 6 {
            continue;
        }
        let edge_colors: Vec<u32> = (0..graph.edge_count())
            .map(|e| ((idx + e) % 3) as u32)
            .collect();
        let node_colors: Vec<u32> = (0..graph.node_count())
            .map(|v| ((idx + v) % 2) as u32)
            .collect();
        let view = EdgeColoredView::new(&graph, &edge_colors);
        part.reset_from_colors(&node_colors, &mut initial);
        let code = refine(view, &mut part, &initial, &mut ws);
        assert_partition_is_edge_color_equitable("simple_g6", &graph, &edge_colors, &part);

        let (pgraph, pnodes, pedges) =
            permute_colored_graph(&graph, &node_colors, &edge_colors, &mut rng);
        let pview = EdgeColoredView::new(&pgraph, &pedges);
        part.reset_from_colors(&pnodes, &mut initial);
        let pcode = refine(pview, &mut part, &initial, &mut ws);
        assert_eq!(pcode, code);
        checked += 1;
    }
    assert_eq!(checked, 1 + 2 + 4 + 11 + 34 + 156);
}

#[rstest]
#[case::max_edge_color_square(
    Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
    vec![0, 0, 0, 0],
    vec![u32::MAX, 0, u32::MAX, 0],
    AutomorphismGroupOrder::Exact(4),
    1
)]
#[case::three_extreme_edge_colors_path(
    Graph::new(4, &[[0, 1], [1, 2], [2, 3]]),
    vec![0, 0, 0, 0],
    vec![0, u32::MAX - 1, u32::MAX],
    AutomorphismGroupOrder::Exact(1),
    4
)]
#[case::u32_max_k4_matching(
    Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0], [0, 2], [1, 3]]),
    vec![0, 0, 0, 0],
    vec![u32::MAX, u32::MAX, u32::MAX, u32::MAX, u32::MAX - 1, u32::MAX - 1],
    AutomorphismGroupOrder::Exact(8),
    1
)]
fn test_trail_ir_full_u32_edge_colors(
    #[case] graph: Graph,
    #[case] node_colors: Vec<u32>,
    #[case] edge_colors: Vec<u32>,
    #[case] expected_order: AutomorphismGroupOrder,
    #[case] expected_orbits: usize,
) {
    let mut workspace = TrailIrWorkspace::new();
    let canon = workspace.solve(&graph, &node_colors, &edge_colors, true);
    assert_eq!(canon.group_order(), expected_order);
    assert_eq!(canon.orbit_count(), expected_orbits);
    assert_output_semantics("u32_max", &graph, &node_colors, &edge_colors, true, &canon);

    let auto = workspace.solve(&graph, &node_colors, &edge_colors, false);
    assert_eq!(auto.group_order(), expected_order);
    assert_eq!(auto.orbit_count(), expected_orbits);
    assert_output_semantics("u32_max", &graph, &node_colors, &edge_colors, false, &auto);

    let key =
        edge_colored_canonical_key(&graph, &node_colors, &edge_colors, canon.canonical_labels());
    let mut rng = ChaCha8Rng::seed_from_u64(0x2540_000b_0000_0001);
    for _ in 0..6 {
        let (pgraph, pnodes, pedges) =
            permute_colored_graph(&graph, &node_colors, &edge_colors, &mut rng);
        let pcanon = workspace.solve(&pgraph, &pnodes, &pedges, true);
        assert_eq!(
            edge_colored_canonical_key(&pgraph, &pnodes, &pedges, pcanon.canonical_labels()),
            key
        );
    }
}

#[rstest]
fn test_trail_ir_fixtures_parity() {
    let mut rng = ChaCha8Rng::seed_from_u64(0x2540_0006_0000_0001);
    let mut workspace = TrailIrWorkspace::new();

    for case in representative_fixtures() {
        let nauty_canon =
            solve_weighted_nauty(&case.graph, &case.node_colors, &case.edge_colors, true);
        let trail_canon = workspace.solve(&case.graph, &case.node_colors, &case.edge_colors, true);
        assert_outputs_agree(
            case.name,
            &case.graph,
            &case.node_colors,
            &case.edge_colors,
            true,
            &nauty_canon,
            &trail_canon,
        );

        let nauty_auto =
            solve_weighted_nauty(&case.graph, &case.node_colors, &case.edge_colors, false);
        let trail_auto = workspace.solve(&case.graph, &case.node_colors, &case.edge_colors, false);
        assert_outputs_agree(
            case.name,
            &case.graph,
            &case.node_colors,
            &case.edge_colors,
            false,
            &nauty_auto,
            &trail_auto,
        );

        let trail_key = edge_colored_canonical_key(
            &case.graph,
            &case.node_colors,
            &case.edge_colors,
            trail_canon.canonical_labels(),
        );
        for _ in 0..6 {
            let (pgraph, pnodes, pedges) =
                permute_colored_graph(&case.graph, &case.node_colors, &case.edge_colors, &mut rng);
            let perm_trail = workspace.solve(&pgraph, &pnodes, &pedges, true);
            let perm_nauty = solve_weighted_nauty(&pgraph, &pnodes, &pedges, true);
            assert_outputs_agree(
                case.name,
                &pgraph,
                &pnodes,
                &pedges,
                true,
                &perm_nauty,
                &perm_trail,
            );
            assert_eq!(
                edge_colored_canonical_key(
                    &pgraph,
                    &pnodes,
                    &pedges,
                    perm_trail.canonical_labels()
                ),
                trail_key,
                "canonical key not relabeling-invariant on {}",
                case.name
            );
        }

        if case.graph.node_count() > 1 {
            let mut stab_colors = case.node_colors.clone();
            stab_colors[0] = stab_colors.iter().copied().max().unwrap_or(0) + 1;
            let stab_nauty =
                solve_weighted_nauty(&case.graph, &stab_colors, &case.edge_colors, true);
            let stab_trail = workspace.solve(&case.graph, &stab_colors, &case.edge_colors, true);
            assert_outputs_agree(
                case.name,
                &case.graph,
                &stab_colors,
                &case.edge_colors,
                true,
                &stab_nauty,
                &stab_trail,
            );
            assert!(stab_trail.generators().iter().all(|g| g[0].index() == 0));
        }
    }
}

#[rstest]
fn test_trail_ir_nucr_complex_parity() {
    let mut rng = ChaCha8Rng::seed_from_u64(0x2540_000a_0000_0001);
    let mut workspace = TrailIrWorkspace::new();
    let mut part = Partition::default();
    let mut refine_ws = RefineWorkspace::default();
    let mut nuc_ws = NucWorkspace::default();
    let mut initial = Vec::new();

    for case in nucr_complex_fixtures() {
        if case.name != "nucr_deep_two_half_scaffold" {
            let view = EdgeColoredView::new(&case.graph, &case.edge_colors);
            part.reset_from_colors(&case.node_colors, &mut initial);
            let _ = refine(view, &mut part, &initial, &mut refine_ws);
            assert!(
                find_non_uniform_components(
                    view,
                    &part,
                    &mut refine_ws,
                    &mut nuc_ws,
                    NUCR_MIN_CELLS
                ),
                "NUCR must activate at root on {}",
                case.name
            );
        }

        let nauty_canon =
            solve_weighted_nauty(&case.graph, &case.node_colors, &case.edge_colors, true);
        let trail_canon = workspace.solve(&case.graph, &case.node_colors, &case.edge_colors, true);
        assert_outputs_agree(
            case.name,
            &case.graph,
            &case.node_colors,
            &case.edge_colors,
            true,
            &nauty_canon,
            &trail_canon,
        );

        let nauty_auto =
            solve_weighted_nauty(&case.graph, &case.node_colors, &case.edge_colors, false);
        let trail_auto = workspace.solve(&case.graph, &case.node_colors, &case.edge_colors, false);
        assert_outputs_agree(
            case.name,
            &case.graph,
            &case.node_colors,
            &case.edge_colors,
            false,
            &nauty_auto,
            &trail_auto,
        );

        let key = edge_colored_canonical_key(
            &case.graph,
            &case.node_colors,
            &case.edge_colors,
            trail_canon.canonical_labels(),
        );
        for _ in 0..4 {
            let (pgraph, pnodes, pedges) =
                permute_colored_graph(&case.graph, &case.node_colors, &case.edge_colors, &mut rng);
            let perm_trail = workspace.solve(&pgraph, &pnodes, &pedges, true);
            assert_eq!(
                edge_colored_canonical_key(
                    &pgraph,
                    &pnodes,
                    &pedges,
                    perm_trail.canonical_labels()
                ),
                key,
                "NUCR canonical key not relabeling-invariant on {}",
                case.name
            );
        }
    }
}

#[rstest]
fn test_trail_ir_simple_through_8_parity() {
    let mut rng = ChaCha8Rng::seed_from_u64(0x2540_0007_0000_0001);
    let mut workspace = TrailIrWorkspace::new();
    let mut checked = 0usize;

    for (idx, graph) in simple_graphs().enumerate() {
        let edge_mod = (idx % 3) + 2;
        let node_mod = (idx % 2) + 1;
        let edge_colors: Vec<u32> = (0..graph.edge_count())
            .map(|e| ((idx + e) % edge_mod) as u32)
            .collect();
        let node_colors: Vec<u32> = (0..graph.node_count())
            .map(|v| ((idx + v) % node_mod) as u32)
            .collect();

        let nauty_canon = solve_weighted_nauty(&graph, &node_colors, &edge_colors, true);
        let trail_canon = workspace.solve(&graph, &node_colors, &edge_colors, true);
        assert_outputs_agree(
            "simple_g6",
            &graph,
            &node_colors,
            &edge_colors,
            true,
            &nauty_canon,
            &trail_canon,
        );

        let nauty_auto = solve_weighted_nauty(&graph, &node_colors, &edge_colors, false);
        let trail_auto = workspace.solve(&graph, &node_colors, &edge_colors, false);
        assert_outputs_agree(
            "simple_g6_auto",
            &graph,
            &node_colors,
            &edge_colors,
            false,
            &nauty_auto,
            &trail_auto,
        );

        if idx % 10 == 0 && graph.node_count() > 1 {
            let key = edge_colored_canonical_key(
                &graph,
                &node_colors,
                &edge_colors,
                trail_canon.canonical_labels(),
            );
            let (pgraph, pnodes, pedges) =
                permute_colored_graph(&graph, &node_colors, &edge_colors, &mut rng);
            let perm_trail = workspace.solve(&pgraph, &pnodes, &pedges, true);
            assert_eq!(
                edge_colored_canonical_key(
                    &pgraph,
                    &pnodes,
                    &pedges,
                    perm_trail.canonical_labels()
                ),
                key,
                "canonical key not relabeling-invariant on simple_g6_{idx}"
            );
        }
        checked += 1;
    }

    assert_eq!(checked, 1 + 2 + 4 + 11 + 34 + 156 + 1_044 + 12_346);
}
