#[path = "matching/fixture.rs"]
#[allow(dead_code)]
mod fixture;

use std::iter;

use num_bigint::BigUint;
use rstest::rstest;
use umol_graph_core::{FaceBoundary, Graph, MatchingEnumerationAlgorithm, NodeId, PlanarEmbedding};
use MatchingEnumerationAlgorithm::BranchAndBound;

#[rstest]
#[case::square(4)]
#[case::hexagon(6)]
#[case::octagon(8)]
fn test_planar_matching_count_cycle(#[case] node_count: usize) {
    let edges: Vec<_> = (0..node_count as u32)
        .map(|node| [node, (node + 1) % node_count as u32])
        .collect();
    let graph = Graph::new(node_count, &edges);
    let forward: Vec<_> = graph.node_ids().collect();
    let reverse: Vec<_> = iter::once(NodeId(0))
        .chain((1..node_count as u32).rev().map(NodeId))
        .collect();
    let faces = [forward, reverse]
        .into_iter()
        .map(|nodes| {
            let boundary_edges = nodes
                .iter()
                .zip(nodes.iter().cycle().skip(1))
                .take(nodes.len())
                .map(|(&first, &second)| graph.find_edge(first, second).unwrap())
                .collect();
            FaceBoundary::new(nodes, boundary_edges)
        })
        .collect();
    let embedding = PlanarEmbedding::new(&graph, faces, 1).unwrap();

    assert_eq!(
        graph.count_perfect_matchings_planar(&embedding),
        Ok(BigUint::from(2_u8))
    );
}

#[rstest]
#[case::k4(
    Graph::new(4, &[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]),
    vec![vec![0, 2, 1], vec![0, 1, 3], vec![0, 3, 2], vec![1, 2, 3]],
    0,
    3,
)]
#[case::cube(
    Graph::new(
        8,
        &[
            [0, 1], [1, 2], [2, 3], [3, 0],
            [4, 5], [5, 6], [6, 7], [7, 4],
            [0, 4], [1, 5], [2, 6], [3, 7],
        ],
    ),
    vec![
        vec![0, 3, 2, 1], vec![4, 5, 6, 7], vec![0, 1, 5, 4],
        vec![1, 2, 6, 5], vec![2, 3, 7, 6], vec![3, 0, 4, 7],
    ],
    0,
    9,
)]
fn test_planar_matching_count_hand(
    #[case] graph: Graph,
    #[case] face_nodes: Vec<Vec<u32>>,
    #[case] outer_face: usize,
    #[case] expected: u32,
) {
    let faces = face_nodes
        .into_iter()
        .map(|nodes| {
            let nodes: Vec<_> = nodes.into_iter().map(NodeId).collect();
            let edges = nodes
                .iter()
                .zip(nodes.iter().cycle().skip(1))
                .take(nodes.len())
                .map(|(&first, &second)| graph.find_edge(first, second).unwrap())
                .collect();
            FaceBoundary::new(nodes, edges)
        })
        .collect();
    let embedding = PlanarEmbedding::new(&graph, faces, outer_face).unwrap();

    assert_eq!(
        graph.count_perfect_matchings_planar(&embedding),
        Ok(BigUint::from(expected))
    );
}

#[rstest]
fn test_planar_matching_count_exhaustive() {
    for node_count in 3..=6 {
        let cycle_edges: Vec<_> = (0..node_count as u32)
            .map(|node| [node, (node + 1) % node_count as u32])
            .collect();
        for mask in 0_usize..(1 << cycle_edges.len()) {
            if mask.count_ones() as usize + 1 < node_count {
                continue;
            }
            let edges: Vec<_> = cycle_edges
                .iter()
                .enumerate()
                .filter_map(|(index, &edge)| (mask & (1 << index) != 0).then_some(edge))
                .collect();
            let graph = Graph::new(node_count, &edges);
            let faces = if edges.len() == node_count {
                let forward: Vec<_> = graph.node_ids().collect();
                let reverse: Vec<_> = iter::once(NodeId(0))
                    .chain((1..node_count as u32).rev().map(NodeId))
                    .collect();
                [forward, reverse]
                    .into_iter()
                    .map(|nodes| {
                        let boundary_edges = nodes
                            .iter()
                            .zip(nodes.iter().cycle().skip(1))
                            .take(nodes.len())
                            .map(|(&first, &second)| graph.find_edge(first, second).unwrap())
                            .collect();
                        FaceBoundary::new(nodes, boundary_edges)
                    })
                    .collect()
            } else {
                let start = graph
                    .node_ids()
                    .find(|&node| graph.degree(node) == 1)
                    .unwrap();
                let mut path = vec![start];
                let mut previous = None;
                while path.len() < node_count {
                    let current = *path.last().unwrap();
                    let next = graph
                        .neighbors(current)
                        .iter()
                        .map(|neighbor| neighbor.node)
                        .find(|&node| Some(node) != previous)
                        .unwrap();
                    previous = Some(current);
                    path.push(next);
                }
                let nodes: Vec<_> = path
                    .iter()
                    .copied()
                    .chain(path[1..path.len() - 1].iter().rev().copied())
                    .collect();
                let boundary_edges = nodes
                    .iter()
                    .zip(nodes.iter().cycle().skip(1))
                    .take(nodes.len())
                    .map(|(&first, &second)| graph.find_edge(first, second).unwrap())
                    .collect();
                vec![FaceBoundary::new(nodes, boundary_edges)]
            };
            let outer_face = faces.len() - 1;
            let embedding = PlanarEmbedding::new(&graph, faces, outer_face).unwrap();
            let mut exhaustive_count = 0_u32;
            for subset in 0_usize..(1 << graph.edge_count()) {
                let mut covered = vec![false; node_count];
                let mut valid = true;
                let mut size = 0;
                for edge in graph.edge_ids() {
                    if subset & (1 << edge.index()) == 0 {
                        continue;
                    }
                    let [first, second] = graph.edge_endpoints(edge);
                    if covered[first.index()] || covered[second.index()] {
                        valid = false;
                        break;
                    }
                    covered[first.index()] = true;
                    covered[second.index()] = true;
                    size += 1;
                }
                if valid && size * 2 == node_count {
                    exhaustive_count += 1;
                }
            }

            assert_eq!(
                graph.count_perfect_matchings_planar(&embedding),
                Ok(BigUint::from(exhaustive_count)),
                "node_count={node_count}, mask={mask}",
            );
        }
    }
}

#[rstest]
#[case::benzene(fixture::BENZENE, 2)]
#[case::naphthalene(fixture::NAPHTHALENE, 3)]
#[case::coronene(fixture::CORONENE, 20)]
#[case::c60(fixture::FULLERENE_C60, 12_500)]
fn test_planar_matching_count_fixture(#[case] source: &str, #[case] known_count: usize) {
    let fixture = fixture::parse(source);
    let graph = fixture.graph();
    let faces: Vec<_> = fixture
        .faces
        .into_iter()
        .map(|nodes| {
            let nodes: Vec<_> = nodes.into_iter().map(NodeId).collect();
            let edges = nodes
                .iter()
                .zip(nodes.iter().cycle().skip(1))
                .take(nodes.len())
                .map(|(&first, &second)| graph.find_edge(first, second).unwrap())
                .collect();
            FaceBoundary::new(nodes, edges)
        })
        .collect();
    let embedding = PlanarEmbedding::new(&graph, faces.clone(), faces.len() - 1).unwrap();
    let enumerated = graph.enumerate_perfect_matchings(BranchAndBound);

    assert_eq!(
        graph.count_perfect_matchings_planar(&embedding),
        Ok(BigUint::from(known_count))
    );
    assert_eq!(enumerated.len(), known_count);
}

#[rstest]
#[case::pyrrole(5)]
#[case::borepine(7)]
fn test_planar_matching_count_holes(#[case] ring_size: usize) {
    let edges: Vec<_> = (0..ring_size as u32)
        .map(|node| [node, (node + 1) % ring_size as u32])
        .collect();
    let graph = Graph::new(ring_size, &edges);
    let retained: Vec<_> = graph.node_ids().filter(|&node| node != NodeId(0)).collect();
    let correspondence = graph.induced_subgraph(&retained);
    let residual = graph.extract(&correspondence);
    let path: Vec<_> = residual.node_ids().collect();
    let nodes: Vec<_> = path
        .iter()
        .copied()
        .chain(path[1..path.len() - 1].iter().rev().copied())
        .collect();
    let boundary_edges = nodes
        .iter()
        .zip(nodes.iter().cycle().skip(1))
        .take(nodes.len())
        .map(|(&first, &second)| residual.find_edge(first, second).unwrap())
        .collect();
    let embedding =
        PlanarEmbedding::new(&residual, vec![FaceBoundary::new(nodes, boundary_edges)], 0).unwrap();

    assert_eq!(
        residual.count_perfect_matchings_planar(&embedding),
        Ok(BigUint::from(1_u8))
    );
    assert_eq!(
        residual.enumerate_perfect_matchings(BranchAndBound).len(),
        1
    );
}

#[rstest]
fn test_planar_matching_count_mobile_hole() {
    let graph = Graph::new(5, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 0]]);
    let mut count = BigUint::from(0_u8);

    for hole in graph.node_ids() {
        let retained: Vec<_> = graph.node_ids().filter(|&node| node != hole).collect();
        let correspondence = graph.induced_subgraph(&retained);
        let residual = graph.extract(&correspondence);
        let start = residual
            .node_ids()
            .find(|&node| residual.degree(node) == 1)
            .unwrap();
        let mut path = vec![start];
        let mut previous = None;
        while path.len() < residual.node_count() {
            let current = *path.last().unwrap();
            let next = residual
                .neighbors(current)
                .iter()
                .map(|neighbor| neighbor.node)
                .find(|&node| Some(node) != previous)
                .unwrap();
            previous = Some(current);
            path.push(next);
        }
        let nodes: Vec<_> = path
            .iter()
            .copied()
            .chain(path[1..path.len() - 1].iter().rev().copied())
            .collect();
        let boundary_edges = nodes
            .iter()
            .zip(nodes.iter().cycle().skip(1))
            .take(nodes.len())
            .map(|(&first, &second)| residual.find_edge(first, second).unwrap())
            .collect();
        let embedding =
            PlanarEmbedding::new(&residual, vec![FaceBoundary::new(nodes, boundary_edges)], 0)
                .unwrap();
        count += residual.count_perfect_matchings_planar(&embedding).unwrap();
    }

    assert_eq!(count, BigUint::from(5_u8));
    assert_eq!(
        count,
        BigUint::from(graph.enumerate_maximum_matchings(BranchAndBound).len())
    );
}

#[rstest]
fn test_planar_matching_count_disconnected() {
    let component_count = |node_count: usize| {
        let edges: Vec<_> = (0..node_count as u32)
            .map(|node| [node, (node + 1) % node_count as u32])
            .collect();
        let graph = Graph::new(node_count, &edges);
        let forward: Vec<_> = graph.node_ids().collect();
        let reverse: Vec<_> = iter::once(NodeId(0))
            .chain((1..node_count as u32).rev().map(NodeId))
            .collect();
        let faces = [forward, reverse]
            .into_iter()
            .map(|nodes| {
                let boundary_edges = nodes
                    .iter()
                    .zip(nodes.iter().cycle().skip(1))
                    .take(nodes.len())
                    .map(|(&first, &second)| graph.find_edge(first, second).unwrap())
                    .collect();
                FaceBoundary::new(nodes, boundary_edges)
            })
            .collect();
        let embedding = PlanarEmbedding::new(&graph, faces, 1).unwrap();
        graph.count_perfect_matchings_planar(&embedding).unwrap()
    };
    let disconnected = Graph::new(
        10,
        &[
            [0, 1],
            [1, 2],
            [2, 3],
            [3, 0],
            [4, 5],
            [5, 6],
            [6, 7],
            [7, 8],
            [8, 9],
            [9, 4],
        ],
    );
    let product = component_count(4) * component_count(6);

    assert_eq!(product, BigUint::from(4_u8));
    assert_eq!(
        product,
        BigUint::from(
            disconnected
                .enumerate_perfect_matchings(BranchAndBound)
                .len()
        )
    );
}

#[rstest]
fn test_planar_matching_count_overflow() {
    const LADDER_COLUMNS: usize = 186;

    let mut edges = Vec::new();
    for column in 0..LADDER_COLUMNS - 1 {
        edges.push([column as u32, column as u32 + 1]);
        edges.push([
            (LADDER_COLUMNS + column) as u32,
            (LADDER_COLUMNS + column + 1) as u32,
        ]);
    }
    for column in 0..LADDER_COLUMNS {
        edges.push([column as u32, (LADDER_COLUMNS + column) as u32]);
    }
    let graph = Graph::new(2 * LADDER_COLUMNS, &edges);
    let mut face_nodes = Vec::new();
    for column in 0..LADDER_COLUMNS - 1 {
        face_nodes.push(vec![
            NodeId(column as u32),
            NodeId(column as u32 + 1),
            NodeId((LADDER_COLUMNS + column + 1) as u32),
            NodeId((LADDER_COLUMNS + column) as u32),
        ]);
    }
    let outer: Vec<_> = iter::once(NodeId(0))
        .chain((LADDER_COLUMNS..2 * LADDER_COLUMNS).map(NodeId::from))
        .chain((1..LADDER_COLUMNS).rev().map(NodeId::from))
        .collect();
    face_nodes.push(outer);
    let faces = face_nodes
        .into_iter()
        .map(|nodes| {
            let boundary_edges = nodes
                .iter()
                .zip(nodes.iter().cycle().skip(1))
                .take(nodes.len())
                .map(|(&first, &second)| graph.find_edge(first, second).unwrap())
                .collect();
            FaceBoundary::new(nodes, boundary_edges)
        })
        .collect();
    let embedding = PlanarEmbedding::new(&graph, faces, LADDER_COLUMNS - 1).unwrap();
    let mut previous = BigUint::from(1_u8);
    let mut expected = BigUint::from(1_u8);
    for _ in 2..=LADDER_COLUMNS {
        let next = &previous + &expected;
        previous = expected;
        expected = next;
    }
    let actual = graph.count_perfect_matchings_planar(&embedding).unwrap();

    assert_eq!(actual, expected);
    assert!(actual > BigUint::from(u128::MAX));
}
