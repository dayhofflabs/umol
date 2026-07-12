mod fixture;

use std::collections::BTreeMap;

use fixture::{parse, GraphFixture};
use rstest::rstest;

#[rstest]
#[case::benzene(fixture::BENZENE, 6, 6, vec![2, 2, 2, 2, 2, 2])]
#[case::naphthalene(
    fixture::NAPHTHALENE,
    10,
    11,
    vec![2, 2, 2, 2, 2, 2, 2, 2, 3, 3]
)]
#[case::coronene(
    fixture::CORONENE,
    24,
    30,
    vec![2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
)]
#[case::azulene(
    fixture::AZULENE,
    10,
    11,
    vec![2, 2, 2, 2, 2, 2, 2, 2, 3, 3]
)]
#[case::c60(fixture::FULLERENE_C60, 60, 90, vec![3; 60])]
#[case::disconnected_cycles(
    fixture::DISCONNECTED_CYCLES,
    24,
    24,
    vec![2; 24]
)]
#[case::ladder(
    fixture::LADDER,
    8,
    10,
    vec![2, 2, 2, 2, 3, 3, 3, 3]
)]
#[case::grid(
    fixture::GRID,
    9,
    12,
    vec![2, 2, 2, 2, 3, 3, 3, 3, 4]
)]
fn test_parse(
    #[case] source: &str,
    #[case] expected_node_count: usize,
    #[case] expected_edge_count: usize,
    #[case] expected_degrees: Vec<usize>,
) {
    let fixture = parse(source);
    let graph = fixture.graph();
    let mut degrees: Vec<_> = graph
        .node_ids()
        .map(|node| graph.neighbors(node).len())
        .collect();
    degrees.sort_unstable();

    assert_eq!(fixture.node_count, expected_node_count);
    assert_eq!(fixture.edges.len(), expected_edge_count);
    assert_eq!(graph.node_count(), expected_node_count);
    assert_eq!(graph.edge_count(), expected_edge_count);
    assert_eq!(degrees, expected_degrees);
}

#[rstest]
#[case::benzene(fixture::BENZENE)]
#[case::naphthalene(fixture::NAPHTHALENE)]
#[case::coronene(fixture::CORONENE)]
#[case::azulene(fixture::AZULENE)]
#[case::c60(fixture::FULLERENE_C60)]
#[case::ladder(fixture::LADDER)]
#[case::grid(fixture::GRID)]
fn test_parse_embedding(#[case] source: &str) {
    let GraphFixture {
        node_count,
        edges,
        faces,
    } = parse(source);
    let graph_edges: BTreeMap<_, _> = edges
        .iter()
        .map(|&[a, b]| ((a.min(b), a.max(b)), 0usize))
        .collect();
    let expected_incidence: BTreeMap<_, _> = edges
        .iter()
        .map(|&[a, b]| ((a.min(b), a.max(b)), 2usize))
        .collect();
    let mut face_incidence = graph_edges.clone();

    for face in &faces {
        for (&a, &b) in face
            .iter()
            .zip(face.iter().cycle().skip(1))
            .take(face.len())
        {
            *face_incidence
                .get_mut(&(a.min(b), a.max(b)))
                .expect("every face side must be a graph edge") += 1;
        }
    }

    assert_eq!(
        node_count as isize - edges.len() as isize + faces.len() as isize,
        2
    );
    assert_eq!(face_incidence, expected_incidence);
}
