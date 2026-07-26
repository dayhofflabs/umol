//! Cycle-family results stated in the ring-perception literature.
//!
//! The fixtures cover Vismara (1997), Figure 4; Hanser, Jauffret, and Kaufmann
//! (1996), Figure 11; Flachsenberg, Andresen, and Rarey (2017), Figure 1;
//! Berger et al. (2004), Figures 7a, 7c, 8a, 8b, 9a, 9c, and 10; May and
//! Steinbeck (2014), Figure 1; and Kolodzik, Urbaczek, and Rarey (2012),
//! Figures 1, 3, and 4. The corresponding papers are retained under
//! `materials/rings`.

use std::ops::ControlFlow;

use rstest::{fixture, rstest};
use umol_graph_core::{
    Cycle, Graph, MinimumCycleBasisAlgorithm, RelevantCycleEnumerationAlgorithm,
    SimpleCycleEnumerationAlgorithm, UniqueRingFamilyAlgorithm,
};

use super::exhaustive::enumerate_cycles;

fn cycle_edge_sets(cycles: impl IntoIterator<Item = Cycle>) -> Vec<Vec<u32>> {
    let mut result = cycles
        .into_iter()
        .map(|cycle| {
            let mut edges = cycle
                .edges()
                .iter()
                .map(|edge| edge.index() as u32)
                .collect::<Vec<_>>();
            edges.sort_unstable();
            edges
        })
        .collect::<Vec<_>>();
    result.sort_by(|left, right| left.len().cmp(&right.len()).then_with(|| left.cmp(right)));
    result
}

#[fixture]
fn kolodzik_figure_3() -> Graph {
    Graph::new(
        16,
        &[
            [0, 1],
            [1, 2],
            [0, 3],
            [3, 4],
            [4, 5],
            [5, 6],
            [6, 7],
            [7, 8],
            [8, 9],
            [9, 10],
            [10, 11],
            [11, 12],
            [12, 13],
            [5, 14],
            [14, 15],
            [15, 0],
            [11, 2],
            [13, 8],
        ],
    )
}

#[fixture]
fn kolodzik_figure_4() -> Graph {
    let mut edges = Vec::with_capacity(48);
    for ring in 0..8_u32 {
        let first = ring;
        let second = (ring + 1) % 8;
        let internal = 8 + 4 * ring;
        edges.extend([
            [first, internal],
            [internal, internal + 1],
            [internal + 1, second],
            [first, internal + 2],
            [internal + 2, internal + 3],
            [internal + 3, second],
        ]);
    }
    Graph::new(40, &edges)
}

#[fixture]
fn gleiss_figure_3() -> Graph {
    Graph::new(
        46,
        &[
            [0, 1],
            [1, 2],
            [2, 3],
            [3, 4],
            [4, 5],
            [4, 6],
            [6, 7],
            [7, 8],
            [8, 9],
            [9, 10],
            [10, 11],
            [11, 12],
            [12, 13],
            [13, 14],
            [14, 15],
            [15, 16],
            [16, 17],
            [17, 18],
            [18, 19],
            [19, 20],
            [20, 21],
            [21, 22],
            [22, 23],
            [23, 24],
            [24, 25],
            [25, 26],
            [26, 27],
            [27, 28],
            [28, 29],
            [29, 30],
            [30, 31],
            [16, 32],
            [32, 33],
            [33, 34],
            [34, 35],
            [35, 36],
            [36, 37],
            [37, 38],
            [38, 39],
            [39, 40],
            [40, 41],
            [41, 42],
            [42, 43],
            [43, 44],
            [44, 45],
            [5, 0],
            [13, 3],
            [6, 25],
            [7, 12],
            [18, 13],
            [39, 15],
            [31, 18],
            [19, 24],
            [26, 31],
            [45, 32],
            [33, 38],
            [40, 45],
        ],
    )
}

#[fixture]
fn champetier_figure_3() -> Graph {
    Graph::new(
        21,
        &[
            [0, 1],
            [0, 8],
            [0, 12],
            [0, 13],
            [0, 14],
            [0, 19],
            [0, 20],
            [1, 2],
            [1, 8],
            [1, 18],
            [1, 19],
            [2, 3],
            [2, 8],
            [2, 11],
            [2, 12],
            [2, 15],
            [2, 18],
            [3, 4],
            [3, 9],
            [3, 11],
            [3, 15],
            [4, 5],
            [4, 9],
            [4, 15],
            [4, 16],
            [5, 6],
            [5, 9],
            [5, 10],
            [5, 16],
            [5, 17],
            [6, 7],
            [6, 10],
            [6, 14],
            [6, 17],
            [6, 20],
            [7, 14],
            [7, 20],
            [8, 12],
            [9, 10],
            [9, 11],
            [10, 11],
            [10, 13],
            [10, 14],
            [11, 12],
            [11, 13],
            [12, 13],
            [13, 14],
            [15, 16],
            [15, 18],
            [16, 17],
            [16, 18],
            [17, 18],
            [17, 19],
            [17, 20],
            [18, 19],
            [19, 20],
        ],
    )
}

#[fixture]
fn berger_figure_8a() -> Graph {
    Graph::new(
        18,
        &[
            [0, 1],
            [0, 8],
            [1, 2],
            [1, 17],
            [2, 3],
            [3, 4],
            [3, 16],
            [4, 5],
            [5, 6],
            [5, 15],
            [6, 7],
            [6, 13],
            [7, 8],
            [7, 11],
            [8, 9],
            [9, 10],
            [9, 17],
            [10, 11],
            [11, 12],
            [12, 13],
            [13, 14],
            [14, 15],
            [15, 16],
            [16, 17],
        ],
    )
}

#[fixture]
fn berger_figure_8b() -> Graph {
    Graph::new(
        16,
        &[
            [0, 1],
            [0, 6],
            [0, 7],
            [1, 2],
            [1, 9],
            [2, 11],
            [2, 3],
            [3, 4],
            [4, 5],
            [5, 6],
            [6, 13],
            [7, 8],
            [7, 14],
            [7, 15],
            [8, 9],
            [9, 10],
            [9, 15],
            [10, 11],
            [11, 12],
            [11, 15],
            [12, 13],
            [13, 14],
            [13, 15],
        ],
    )
}

#[fixture]
fn berger_figure_9a() -> Graph {
    Graph::new(
        12,
        &[
            [0, 1],
            [0, 7],
            [1, 2],
            [1, 3],
            [1, 7],
            [1, 8],
            [1, 9],
            [2, 3],
            [3, 4],
            [3, 5],
            [3, 9],
            [3, 10],
            [4, 5],
            [5, 6],
            [5, 7],
            [5, 10],
            [5, 11],
            [6, 7],
            [7, 8],
            [7, 11],
        ],
    )
}

#[fixture]
fn berger_figure_9c() -> Graph {
    Graph::new(
        24,
        &[
            [0, 1],
            [1, 2],
            [2, 3],
            [0, 4],
            [4, 5],
            [5, 6],
            [6, 7],
            [4, 8],
            [8, 9],
            [9, 10],
            [10, 11],
            [11, 12],
            [9, 13],
            [13, 14],
            [14, 15],
            [9, 16],
            [16, 17],
            [17, 18],
            [18, 19],
            [19, 20],
            [17, 21],
            [21, 22],
            [22, 23],
            [17, 0],
            [20, 0],
            [4, 3],
            [8, 7],
            [15, 8],
            [16, 12],
            [23, 16],
        ],
    )
}

#[fixture]
fn berger_figure_10() -> Graph {
    Graph::new(
        90,
        &[
            [0, 1],
            [1, 2],
            [2, 3],
            [3, 4],
            [4, 5],
            [5, 6],
            [6, 7],
            [7, 8],
            [8, 9],
            [9, 10],
            [10, 11],
            [11, 12],
            [12, 13],
            [13, 14],
            [14, 15],
            [15, 16],
            [16, 17],
            [17, 18],
            [18, 19],
            [19, 20],
            [20, 21],
            [21, 22],
            [22, 23],
            [23, 24],
            [24, 25],
            [25, 26],
            [26, 27],
            [27, 28],
            [28, 29],
            [29, 30],
            [30, 31],
            [31, 32],
            [32, 33],
            [33, 34],
            [34, 35],
            [35, 36],
            [36, 37],
            [37, 38],
            [38, 39],
            [39, 40],
            [40, 41],
            [41, 42],
            [42, 43],
            [43, 44],
            [44, 45],
            [45, 46],
            [46, 47],
            [47, 48],
            [48, 49],
            [49, 50],
            [50, 51],
            [51, 52],
            [52, 53],
            [53, 54],
            [54, 55],
            [55, 56],
            [56, 57],
            [57, 58],
            [58, 59],
            [59, 60],
            [60, 61],
            [61, 62],
            [62, 63],
            [63, 64],
            [64, 65],
            [62, 66],
            [66, 67],
            [67, 68],
            [68, 69],
            [69, 70],
            [70, 71],
            [71, 72],
            [72, 73],
            [73, 74],
            [74, 75],
            [75, 76],
            [76, 77],
            [77, 78],
            [78, 79],
            [79, 80],
            [80, 81],
            [81, 82],
            [82, 83],
            [83, 84],
            [83, 85],
            [85, 86],
            [71, 87],
            [69, 88],
            [34, 89],
            [4, 0],
            [0, 37],
            [35, 1],
            [89, 2],
            [13, 3],
            [38, 5],
            [12, 5],
            [6, 10],
            [7, 39],
            [41, 8],
            [49, 9],
            [50, 11],
            [81, 12],
            [80, 14],
            [89, 14],
            [28, 15],
            [20, 16],
            [89, 17],
            [18, 33],
            [65, 19],
            [28, 21],
            [64, 21],
            [22, 26],
            [23, 63],
            [61, 24],
            [59, 25],
            [86, 27],
            [85, 29],
            [80, 29],
            [64, 30],
            [66, 30],
            [79, 31],
            [32, 65],
            [79, 34],
            [78, 36],
            [77, 38],
            [76, 40],
            [49, 42],
            [75, 42],
            [73, 43],
            [48, 44],
            [72, 45],
            [87, 46],
            [51, 47],
            [82, 50],
            [83, 52],
            [87, 52],
            [86, 53],
            [58, 54],
            [87, 55],
            [56, 70],
            [88, 57],
            [86, 59],
            [88, 60],
            [67, 60],
            [85, 67],
            [84, 68],
            [84, 71],
            [84, 74],
            [82, 75],
            [81, 77],
        ],
    )
}

#[fixture]
fn vismara_figure_4() -> Graph {
    Graph::new(
        20,
        &[
            [0, 1],
            [0, 11],
            [1, 2],
            [1, 7],
            [2, 3],
            [2, 8],
            [3, 4],
            [4, 5],
            [4, 8],
            [5, 6],
            [5, 10],
            [6, 16],
            [7, 8],
            [7, 9],
            [9, 10],
            [11, 12],
            [11, 17],
            [12, 13],
            [13, 14],
            [13, 17],
            [14, 15],
            [14, 18],
            [14, 19],
            [15, 16],
            [16, 18],
            [16, 19],
        ],
    )
}

#[fixture]
fn hanser_figure_11() -> Graph {
    Graph::new(
        9,
        &[
            [0, 3],
            [0, 4],
            [0, 6],
            [1, 4],
            [1, 5],
            [1, 7],
            [2, 3],
            [2, 5],
            [2, 8],
            [3, 6],
            [3, 8],
            [4, 6],
            [4, 7],
            [5, 7],
            [5, 8],
            [6, 7],
            [6, 8],
            [7, 8],
        ],
    )
}

#[rstest]
fn test_graph_enumerate_simple_cycles_hanser_figure_11(#[from(hanser_figure_11)] graph: Graph) {
    let actual = cycle_edge_sets(
        graph.enumerate_simple_cycles(usize::MAX, SimpleCycleEnumerationAlgorithm::ReadTarjan),
    );
    let expected = enumerate_cycles(&graph)
        .into_iter()
        .map(|cycle| {
            cycle
                .into_iter()
                .map(|edge| edge.index() as u32)
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let actual_degrees = graph
        .node_ids()
        .map(|node| graph.degree(node))
        .collect::<Vec<_>>();

    assert_eq!(actual_degrees, [3, 3, 3, 4, 4, 4, 5, 5, 5]);
    assert_eq!(actual, expected);
    assert_eq!(actual.len(), 248);
}

#[rstest]
#[case::ring_decomposer_lib_figure_1(
    8,
    &[
        [0, 1], [1, 2], [2, 3], [3, 4], [4, 5],
        [5, 0], [0, 6], [6, 2], [3, 7], [7, 5],
    ],
    &[
        vec![0, 1, 6, 7],
        vec![3, 4, 8, 9],
        vec![0, 1, 2, 3, 4, 5],
        vec![2, 3, 4, 5, 6, 7],
        vec![0, 1, 2, 5, 8, 9],
        vec![2, 5, 6, 7, 8, 9],
    ],
)]
#[case::berger_figure_7a(
    8,
    &[
        [0, 1], [1, 2], [2, 3], [3, 4], [4, 5],
        [5, 6], [6, 7], [7, 0], [1, 5], [0, 4],
    ],
    &[
        vec![0, 4, 8, 9],
        vec![1, 2, 3, 4, 8],
        vec![0, 5, 6, 7, 8],
        vec![0, 1, 2, 3, 9],
        vec![4, 5, 6, 7, 9],
    ],
)]
#[case::berger_figure_7c(
    5,
    &[
        [0, 1], [0, 2], [2, 1], [0, 3],
        [3, 1], [0, 4], [4, 1],
    ],
    &[
        vec![0, 1, 2],
        vec![0, 3, 4],
        vec![0, 5, 6],
    ],
)]
#[case::may_figure_1(
    8,
    &[
        [0, 2], [2, 3], [3, 1],
        [0, 4], [4, 5], [5, 1],
        [0, 6], [6, 7], [7, 1],
    ],
    &[
        vec![0, 1, 2, 3, 4, 5],
        vec![0, 1, 2, 6, 7, 8],
        vec![3, 4, 5, 6, 7, 8],
    ],
)]
#[case::kolodzik_figure_1(
    8,
    &[
        [0, 1], [0, 2], [0, 4], [1, 3],
        [1, 5], [2, 3], [2, 6], [3, 7],
        [4, 5], [4, 6], [5, 7], [6, 7],
    ],
    &[
        vec![0, 1, 3, 5],
        vec![8, 9, 10, 11],
        vec![0, 2, 4, 8],
        vec![5, 6, 7, 11],
        vec![1, 2, 6, 9],
        vec![3, 4, 7, 10],
    ],
)]
fn test_graph_enumerate_relevant_cycles_literature(
    #[case] node_count: usize,
    #[case] edges: &[[u32; 2]],
    #[case] expected: &[Vec<u32>],
) {
    let graph = Graph::new(node_count, edges);
    let actual = cycle_edge_sets(
        graph.enumerate_relevant_cycles(usize::MAX, RelevantCycleEnumerationAlgorithm::Vismara),
    );
    let mut expected = expected.to_vec();
    expected.sort_by(|left, right| left.len().cmp(&right.len()).then_with(|| left.cmp(right)));

    assert_eq!(actual, expected);
}

#[rstest]
fn test_graph_enumerate_relevant_cycles_kolodzik_figure_3(#[from(kolodzik_figure_3)] graph: Graph) {
    let actual = cycle_edge_sets(
        graph.enumerate_relevant_cycles(usize::MAX, RelevantCycleEnumerationAlgorithm::Vismara),
    );
    let expected = [
        vec![2, 3, 4, 13, 14, 15],
        vec![8, 9, 10, 11, 12, 17],
        vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 16],
        vec![0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 16, 17],
        vec![0, 1, 5, 6, 7, 8, 9, 10, 13, 14, 15, 16],
        vec![0, 1, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17],
    ];

    assert_eq!(actual, expected);
}

#[rstest]
fn test_graph_enumerate_relevant_cycles_macrocycle(#[from(kolodzik_figure_4)] graph: Graph) {
    let actual = cycle_edge_sets(
        graph.enumerate_relevant_cycles(usize::MAX, RelevantCycleEnumerationAlgorithm::Vismara),
    );
    let mut expected = (0..8_u32)
        .map(|ring| (6 * ring..6 * ring + 6).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    for choices in 0..256_u32 {
        let mut cycle = Vec::with_capacity(24);
        for ring in 0..8_u32 {
            let start = 6 * ring + if choices & (1 << ring) == 0 { 0 } else { 3 };
            cycle.extend(start..start + 3);
        }
        expected.push(cycle);
    }
    expected.sort_by(|left, right| left.len().cmp(&right.len()).then_with(|| left.cmp(right)));

    assert_eq!(actual, expected);
}

#[rstest]
fn test_graph_enumerate_relevant_cycles_gleiss(#[from(gleiss_figure_3)] graph: Graph) {
    let actual = cycle_edge_sets(
        graph.enumerate_relevant_cycles(usize::MAX, RelevantCycleEnumerationAlgorithm::Vismara),
    );
    let expected = [
        vec![0, 1, 2, 3, 4, 45],
        vec![3, 5, 6, 12, 46, 48],
        vec![7, 8, 9, 10, 11, 48],
        vec![13, 14, 15, 16, 17, 49],
        vec![15, 31, 32, 38, 50, 55],
        vec![15, 31, 39, 50, 54, 56],
        vec![18, 24, 25, 51, 52, 53],
        vec![19, 20, 21, 22, 23, 52],
        vec![26, 27, 28, 29, 30, 53],
        vec![32, 38, 39, 54, 55, 56],
        vec![33, 34, 35, 36, 37, 55],
        vec![40, 41, 42, 43, 44, 56],
        vec![3, 5, 18, 24, 46, 47, 49, 52],
        vec![3, 5, 25, 46, 47, 49, 51, 53],
        vec![6, 12, 18, 24, 47, 48, 49, 52],
        vec![6, 12, 25, 47, 48, 49, 51, 53],
    ];

    assert_eq!(actual, expected);
}

#[rstest]
fn test_graph_enumerate_relevant_cycles_champetier(#[from(champetier_figure_3)] graph: Graph) {
    let actual = cycle_edge_sets(
        graph.enumerate_relevant_cycles(usize::MAX, RelevantCycleEnumerationAlgorithm::Vismara),
    );
    let nodes = graph.node_ids().collect::<Vec<_>>();
    let mut expected = Vec::new();
    for first in 0..nodes.len() {
        for second in first + 1..nodes.len() {
            for third in second + 1..nodes.len() {
                let (Some(first_second), Some(second_third), Some(first_third)) = (
                    graph.find_edge(nodes[first], nodes[second]),
                    graph.find_edge(nodes[second], nodes[third]),
                    graph.find_edge(nodes[first], nodes[third]),
                ) else {
                    continue;
                };
                let mut triangle = vec![
                    first_second.index() as u32,
                    second_third.index() as u32,
                    first_third.index() as u32,
                ];
                triangle.sort_unstable();
                expected.push(triangle);
            }
        }
    }
    expected.sort();

    assert_eq!(actual, expected);
}

#[rstest]
fn test_graph_enumerate_relevant_cycles_berger_figure_8a(#[from(berger_figure_8a)] graph: Graph) {
    let actual = cycle_edge_sets(
        graph.enumerate_relevant_cycles(usize::MAX, RelevantCycleEnumerationAlgorithm::Vismara),
    );
    let expected = [
        vec![0, 1, 3, 14, 16],
        vec![2, 3, 4, 6, 23],
        vec![5, 6, 7, 9, 22],
        vec![8, 9, 11, 20, 21],
        vec![10, 11, 13, 18, 19],
        vec![12, 13, 14, 15, 17],
        vec![8, 9, 10, 12, 14, 16, 22, 23],
    ];
    let n = vec![15, 16, 17, 18, 19, 20, 21, 22, 23];

    assert_eq!(actual, expected);
    assert!(!actual.contains(&n));
}

#[rstest]
fn test_graph_enumerate_relevant_cycles_berger_figure_8b(#[from(berger_figure_8b)] graph: Graph) {
    let actual = cycle_edge_sets(
        graph.enumerate_relevant_cycles(usize::MAX, RelevantCycleEnumerationAlgorithm::Vismara),
    );
    let expected = [
        vec![11, 13, 14, 16],
        vec![12, 13, 21, 22],
        vec![15, 16, 17, 19],
        vec![18, 19, 20, 22],
        vec![0, 2, 4, 11, 14],
        vec![0, 2, 4, 13, 16],
        vec![1, 2, 10, 12, 21],
        vec![1, 2, 10, 13, 22],
        vec![3, 4, 5, 15, 17],
        vec![3, 4, 5, 16, 19],
        vec![0, 1, 3, 6, 7, 8, 9],
    ];
    let c = vec![11, 12, 14, 15, 17, 18, 20, 21];
    let actual_degrees = graph
        .node_ids()
        .map(|node| graph.degree(node))
        .collect::<Vec<_>>();

    assert_eq!(
        actual_degrees,
        [3, 3, 3, 2, 2, 2, 3, 4, 2, 4, 2, 4, 2, 4, 2, 4]
    );
    assert_eq!(actual, expected);
    assert!(!actual.contains(&c));
}

#[rstest]
fn test_graph_enumerate_relevant_cycles_berger_figure_9a(#[from(berger_figure_9a)] graph: Graph) {
    let actual = cycle_edge_sets(
        graph.enumerate_relevant_cycles(usize::MAX, RelevantCycleEnumerationAlgorithm::Vismara),
    );
    let expected = [
        vec![0, 1, 4],
        vec![2, 3, 7],
        vec![3, 6, 10],
        vec![4, 5, 18],
        vec![8, 9, 12],
        vec![9, 11, 15],
        vec![13, 14, 17],
        vec![14, 16, 19],
        vec![3, 4, 9, 14],
    ];

    assert_eq!(actual, expected);
}

#[rstest]
fn test_graph_enumerate_relevant_cycles_berger_figure_9c(#[from(berger_figure_9c)] graph: Graph) {
    let actual = cycle_edge_sets(
        graph.enumerate_relevant_cycles(usize::MAX, RelevantCycleEnumerationAlgorithm::Vismara),
    );
    let expected = [
        vec![0, 1, 2, 3, 25],
        vec![4, 5, 6, 7, 26],
        vec![8, 12, 13, 14, 27],
        vec![9, 10, 11, 15, 28],
        vec![16, 20, 21, 22, 29],
        vec![17, 18, 19, 23, 24],
        vec![3, 7, 8, 15, 16, 23],
    ];

    assert_eq!(actual, expected);
}

#[rstest]
fn test_graph_enumerate_relevant_cycles_berger_figure_10(#[from(berger_figure_10)] graph: Graph) {
    let actual = cycle_edge_sets(
        graph.enumerate_relevant_cycles(usize::MAX, RelevantCycleEnumerationAlgorithm::Vismara),
    );
    let expected = cycle_edge_sets(
        graph.enumerate_simple_cycles(6, SimpleCycleEnumerationAlgorithm::ReadTarjan),
    );
    let actual_lengths = actual.iter().map(Vec::len).collect::<Vec<_>>();

    assert_eq!(graph.node_count(), 90);
    assert_eq!(graph.edge_count(), 150);
    assert_eq!(actual, expected);
    assert_eq!(actual_lengths, [vec![5; 66], vec![6]].concat());
}

#[rstest]
fn test_graph_enumerate_relevant_cycles_vismara_figure_4(#[from(vismara_figure_4)] graph: Graph) {
    let actual =
        graph.enumerate_relevant_cycles(usize::MAX, RelevantCycleEnumerationAlgorithm::Vismara);
    let actual_lengths = actual.iter().map(Cycle::length).collect::<Vec<_>>();
    let family_count = actual
        .iter()
        .filter(|cycle| {
            cycle.length() == 13
                && [0, 6, 16].iter().all(|node| {
                    cycle
                        .nodes()
                        .iter()
                        .any(|candidate| candidate.index() == *node)
                })
        })
        .count();
    let actual_degrees = graph
        .node_ids()
        .map(|node| graph.degree(node))
        .collect::<Vec<_>>();

    assert_eq!(
        actual_degrees,
        [2, 3, 3, 2, 3, 3, 2, 3, 3, 2, 2, 3, 2, 3, 4, 2, 4, 2, 2, 2]
    );
    assert_eq!(actual_lengths, [vec![4; 6], vec![6], vec![13; 24]].concat());
    assert_eq!(family_count, 24);
}

#[rstest]
#[case::ring_decomposer_lib_figure_1(
    8,
    &[
        [0, 1], [1, 2], [2, 3], [3, 4], [4, 5],
        [5, 0], [0, 6], [6, 2], [3, 7], [7, 5],
    ],
    &[4, 4, 6],
)]
#[case::berger_figure_7a(
    8,
    &[
        [0, 1], [1, 2], [2, 3], [3, 4], [4, 5],
        [5, 6], [6, 7], [7, 0], [1, 5], [0, 4],
    ],
    &[4, 5, 5],
)]
#[case::berger_figure_7c(
    5,
    &[
        [0, 1], [0, 2], [2, 1], [0, 3],
        [3, 1], [0, 4], [4, 1],
    ],
    &[3, 3, 3],
)]
#[case::may_figure_1(
    8,
    &[
        [0, 2], [2, 3], [3, 1],
        [0, 4], [4, 5], [5, 1],
        [0, 6], [6, 7], [7, 1],
    ],
    &[6, 6],
)]
#[case::kolodzik_figure_1(
    8,
    &[
        [0, 1], [0, 2], [0, 4], [1, 3],
        [1, 5], [2, 3], [2, 6], [3, 7],
        [4, 5], [4, 6], [5, 7], [6, 7],
    ],
    &[4, 4, 4, 4, 4],
)]
fn test_graph_minimum_cycle_basis_literature(
    #[case] node_count: usize,
    #[case] edges: &[[u32; 2]],
    #[case] expected_lengths: &[usize],
) {
    let graph = Graph::new(node_count, edges);
    let basis = graph.minimum_cycle_basis(MinimumCycleBasisAlgorithm::Horton);
    let mut actual_lengths = basis.iter().map(Cycle::length).collect::<Vec<_>>();
    actual_lengths.sort_unstable();

    assert_eq!(actual_lengths, expected_lengths);
    assert_eq!(basis.dimension(), expected_lengths.len());
    assert_eq!(basis.total_length(), expected_lengths.iter().sum::<usize>());
}

#[rstest]
fn test_graph_minimum_cycle_basis_kolodzik_figure_3(#[from(kolodzik_figure_3)] graph: Graph) {
    let basis = graph.minimum_cycle_basis(MinimumCycleBasisAlgorithm::Horton);
    let mut actual_lengths = basis.iter().map(Cycle::length).collect::<Vec<_>>();
    actual_lengths.sort_unstable();

    assert_eq!(graph.node_count(), 16);
    assert_eq!(graph.edge_count(), 18);
    assert_eq!(actual_lengths, [6, 6, 12]);
    assert_eq!(basis.dimension(), 3);
    assert_eq!(basis.total_length(), 24);
}

#[rstest]
fn test_graph_minimum_cycle_basis_macrocycle(#[from(kolodzik_figure_4)] graph: Graph) {
    let basis = graph.minimum_cycle_basis(MinimumCycleBasisAlgorithm::Horton);
    let mut actual_lengths = basis.iter().map(Cycle::length).collect::<Vec<_>>();
    actual_lengths.sort_unstable();

    assert_eq!(actual_lengths, [6, 6, 6, 6, 6, 6, 6, 6, 24]);
    assert_eq!(basis.dimension(), 9);
    assert_eq!(basis.total_length(), 72);
}

#[rstest]
fn test_graph_minimum_cycle_basis_gleiss(#[from(gleiss_figure_3)] graph: Graph) {
    let basis = graph.minimum_cycle_basis(MinimumCycleBasisAlgorithm::Horton);
    let mut actual_lengths = basis.iter().map(Cycle::length).collect::<Vec<_>>();
    actual_lengths.sort_unstable();

    assert_eq!(graph.node_count(), 46);
    assert_eq!(graph.edge_count(), 57);
    assert_eq!(actual_lengths, [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 8]);
    assert_eq!(basis.dimension(), 12);
    assert_eq!(basis.total_length(), 74);
}

#[rstest]
fn test_graph_minimum_cycle_basis_champetier(#[from(champetier_figure_3)] graph: Graph) {
    let basis = graph.minimum_cycle_basis(MinimumCycleBasisAlgorithm::Horton);
    let mut actual_lengths = basis.iter().map(Cycle::length).collect::<Vec<_>>();
    actual_lengths.sort_unstable();
    let actual_degrees = graph
        .node_ids()
        .map(|node| graph.degree(node))
        .collect::<Vec<_>>();

    assert_eq!(graph.node_count(), 21);
    assert_eq!(graph.edge_count(), 56);
    assert_eq!(
        actual_degrees,
        [7, 5, 7, 5, 5, 6, 6, 3, 4, 5, 6, 6, 5, 5, 5, 5, 5, 6, 6, 5, 5]
    );
    assert_eq!(actual_lengths, [3; 36]);
    assert_eq!(basis.dimension(), 36);
    assert_eq!(basis.total_length(), 108);
}

#[rstest]
fn test_graph_minimum_cycle_basis_berger_figure_8a(#[from(berger_figure_8a)] graph: Graph) {
    let basis = graph.minimum_cycle_basis(MinimumCycleBasisAlgorithm::Horton);
    let mut actual_lengths = basis.iter().map(Cycle::length).collect::<Vec<_>>();
    actual_lengths.sort_unstable();
    let actual_degrees = graph
        .node_ids()
        .map(|node| graph.degree(node))
        .collect::<Vec<_>>();

    assert_eq!(graph.node_count(), 18);
    assert_eq!(graph.edge_count(), 24);
    assert_eq!(
        actual_degrees,
        [2, 3, 2, 3, 2, 3, 3, 3, 3, 3, 2, 3, 2, 3, 2, 3, 3, 3]
    );
    assert_eq!(actual_lengths, [5, 5, 5, 5, 5, 5, 8]);
    assert_eq!(basis.dimension(), 7);
    assert_eq!(basis.total_length(), 38);
}

#[rstest]
fn test_graph_minimum_cycle_basis_berger_figure_9a(#[from(berger_figure_9a)] graph: Graph) {
    let basis = graph.minimum_cycle_basis(MinimumCycleBasisAlgorithm::Horton);
    let mut actual_lengths = basis.iter().map(Cycle::length).collect::<Vec<_>>();
    actual_lengths.sort_unstable();
    let actual_degrees = graph
        .node_ids()
        .map(|node| graph.degree(node))
        .collect::<Vec<_>>();

    assert_eq!(graph.node_count(), 12);
    assert_eq!(graph.edge_count(), 20);
    assert_eq!(actual_degrees, [2, 6, 2, 6, 2, 6, 2, 6, 2, 2, 2, 2]);
    assert_eq!(actual_lengths, [3, 3, 3, 3, 3, 3, 3, 3, 4]);
    assert_eq!(basis.dimension(), 9);
    assert_eq!(basis.total_length(), 28);
}

#[rstest]
fn test_graph_minimum_cycle_basis_berger_figure_9c(#[from(berger_figure_9c)] graph: Graph) {
    let basis = graph.minimum_cycle_basis(MinimumCycleBasisAlgorithm::Horton);
    let mut actual_lengths = basis.iter().map(Cycle::length).collect::<Vec<_>>();
    actual_lengths.sort_unstable();

    assert_eq!(graph.node_count(), 24);
    assert_eq!(graph.edge_count(), 30);
    assert_eq!(actual_lengths, [5, 5, 5, 5, 5, 5, 6]);
    assert_eq!(basis.dimension(), 7);
    assert_eq!(basis.total_length(), 36);
}

#[rstest]
fn test_graph_minimum_cycle_basis_berger_figure_10(#[from(berger_figure_10)] graph: Graph) {
    let basis = graph.minimum_cycle_basis(MinimumCycleBasisAlgorithm::Horton);
    let mut actual_lengths = basis.iter().map(Cycle::length).collect::<Vec<_>>();
    actual_lengths.sort_unstable();

    assert_eq!(actual_lengths, [vec![5; 60], vec![6]].concat());
    assert_eq!(basis.dimension(), 61);
    assert_eq!(basis.total_length(), 306);
}

#[rstest]
#[case::ring_decomposer_lib_figure_1(
    8,
    &[
        [0, 1], [1, 2], [2, 3], [3, 4], [4, 5],
        [5, 0], [0, 6], [6, 2], [3, 7], [7, 5],
    ],
    &[(4, 1), (4, 1), (6, 4)],
)]
#[case::kolodzik_figure_1(
    8,
    &[
        [0, 1], [0, 2], [0, 4], [1, 3],
        [1, 5], [2, 3], [2, 6], [3, 7],
        [4, 5], [4, 6], [5, 7], [6, 7],
    ],
    &[(4, 1), (4, 1), (4, 1), (4, 1), (4, 1), (4, 1)],
)]
#[case::kolodzik_figure_3(
    16,
    &[
        [0, 1], [1, 2], [0, 3], [3, 4], [4, 5], [5, 6],
        [6, 7], [7, 8], [8, 9], [9, 10], [10, 11], [11, 12],
        [12, 13], [5, 14], [14, 15], [15, 0], [11, 2], [13, 8],
    ],
    &[(6, 1), (6, 1), (12, 4)],
)]
fn test_graph_unique_ring_families_literature(
    #[case] node_count: usize,
    #[case] edges: &[[u32; 2]],
    #[case] expected: &[(usize, usize)],
) {
    let graph = Graph::new(node_count, edges);
    let decomposition = graph.unique_ring_families(UniqueRingFamilyAlgorithm::Kolodzik);
    let mut actual = decomposition
        .ids()
        .map(|id| {
            let family = decomposition
                .get(id)
                .expect("a returned family id must be valid");
            let mut cycles = Vec::new();
            let flow = decomposition.visit_relevant_cycles(id, |cycle| {
                cycles.push(cycle);
                ControlFlow::<()>::Continue(())
            });
            let count = family
                .relevant_cycle_count()
                .0
                .to_string()
                .parse::<usize>()
                .expect("published family count must fit usize");

            assert_eq!(flow, ControlFlow::Continue(()));
            assert_eq!(cycles.len(), count);
            assert!(cycles.iter().all(|cycle| cycle.length() == family.weight()));
            (family.weight(), count)
        })
        .collect::<Vec<_>>();
    actual.sort_unstable();

    assert_eq!(actual, expected);
}

#[rstest]
fn test_graph_unique_ring_families_macrocycle(#[from(kolodzik_figure_4)] graph: Graph) {
    let decomposition = graph.unique_ring_families(UniqueRingFamilyAlgorithm::Kolodzik);
    let mut actual = decomposition
        .ids()
        .map(|id| {
            let family = decomposition
                .get(id)
                .expect("a returned family id must be valid");
            let mut cycles = Vec::new();
            let flow = decomposition.visit_relevant_cycles(id, |cycle| {
                cycles.push(cycle);
                ControlFlow::<()>::Continue(())
            });
            let count = family
                .relevant_cycle_count()
                .0
                .to_string()
                .parse::<usize>()
                .expect("published family count must fit usize");

            assert_eq!(flow, ControlFlow::Continue(()));
            assert_eq!(cycles.len(), count);
            assert!(cycles.iter().all(|cycle| cycle.length() == family.weight()));
            (family.weight(), count)
        })
        .collect::<Vec<_>>();
    actual.sort_unstable();

    assert_eq!(
        actual,
        [
            (6, 1),
            (6, 1),
            (6, 1),
            (6, 1),
            (6, 1),
            (6, 1),
            (6, 1),
            (6, 1),
            (24, 256),
        ]
    );
}
