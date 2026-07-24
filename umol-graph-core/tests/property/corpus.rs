//! Checked-in graph corpora used by property tests.

use umol_graph_core::Graph;

const SIMPLE_GRAPHS_GRAPH6: &str = include_str!("../data/simple-through-8.g6");

pub(super) fn parse_graph6(source: &str) -> Graph {
    let bytes = source.as_bytes();
    assert!(!bytes.is_empty(), "graph6 input must not be empty");
    assert!(
        (63..126).contains(&bytes[0]),
        "the corpus reader supports compact graph6 orders"
    );
    assert!(
        bytes[1..].iter().all(|byte| (63..=126).contains(byte)),
        "graph6 data bytes must be printable six-bit values"
    );

    let node_count = (bytes[0] - 63) as usize;
    let bit_count = node_count * node_count.saturating_sub(1) / 2;
    let data_count = bit_count.div_ceil(6);
    assert_eq!(
        bytes.len(),
        data_count + 1,
        "graph6 input has the wrong data length"
    );

    let mut edges = Vec::new();
    let mut bit = 0;
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

pub(super) fn simple_graphs() -> impl Iterator<Item = Graph> {
    SIMPLE_GRAPHS_GRAPH6.lines().map(parse_graph6)
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use xxhash_rust::xxh3::xxh3_64;

    use super::*;

    #[rstest]
    #[case::empty("?", 0, &[])]
    #[case::singleton("@", 1, &[])]
    #[case::two_isolated("A?", 2, &[])]
    #[case::single_edge("A_", 2, &[[0, 1]])]
    #[case::triangle("Bw", 3, &[[0, 1], [0, 2], [1, 2]])]
    fn test_parse_graph6(
        #[case] source: &str,
        #[case] node_count: usize,
        #[case] edges: &[[u32; 2]],
    ) {
        assert_eq!(parse_graph6(source), Graph::new(node_count, edges));
    }

    #[rstest]
    fn test_simple_graphs() {
        let mut counts = [0usize; 9];
        for graph in simple_graphs() {
            counts[graph.node_count()] += 1;
        }
        assert_eq!(counts, [0, 1, 2, 4, 11, 34, 156, 1_044, 12_346]);
        assert_eq!(
            xxh3_64(SIMPLE_GRAPHS_GRAPH6.as_bytes()),
            0xc3c0_3691_841d_9b70
        );
    }
}
