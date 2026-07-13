#[path = "matching/fixture.rs"]
#[allow(dead_code)]
mod fixture;

use rstest::rstest;
use umol_graph_core::{EdgeId, FaceBoundary, Graph, NodeId, PlanarEmbedding, PlanarEmbeddingError};

#[rstest]
#[case::cycle(
    Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
    vec![vec![0, 1, 2, 3], vec![0, 3, 2, 1]],
    1,
)]
#[case::k4(
    Graph::new(4, &[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]),
    vec![vec![0, 2, 1], vec![0, 1, 3], vec![0, 3, 2], vec![1, 2, 3]],
    0,
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
)]
#[case::bridge(Graph::new(2, &[[0, 1]]), vec![vec![0, 1]], 0)]
fn test_planar_embedding_new(
    #[case] graph: Graph,
    #[case] face_nodes: Vec<Vec<u32>>,
    #[case] outer_face: usize,
) {
    let faces: Vec<_> = face_nodes
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
    let embedding = PlanarEmbedding::new(&graph, faces.clone(), outer_face).unwrap();

    assert_eq!(embedding.graph(), &graph);
    assert_eq!(embedding.faces(), faces);
    assert_eq!(embedding.outer_face_index(), outer_face);
    assert_eq!(embedding.outer_face(), &faces[outer_face]);
    assert_eq!(
        embedding.bounded_faces().collect::<Vec<_>>(),
        faces
            .iter()
            .enumerate()
            .filter_map(|(index, face)| (index != outer_face).then_some(face))
            .collect::<Vec<_>>()
    );
}

#[rstest]
#[case::coronene(fixture::CORONENE)]
#[case::c60(fixture::FULLERENE_C60)]
fn test_planar_embedding_new_fixture(#[case] source: &str) {
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
    let outer_face = faces.len() - 1;
    let embedding = PlanarEmbedding::new(&graph, faces.clone(), outer_face).unwrap();

    assert_eq!(embedding.graph(), &graph);
    assert_eq!(embedding.faces(), faces);
    assert_eq!(embedding.outer_face_index(), outer_face);
    assert_eq!(embedding.outer_face(), &faces[outer_face]);
}

#[rstest]
#[case::outer_face(
    Graph::new(3, &[[0, 1], [1, 2], [2, 0]]),
    vec![
        FaceBoundary::new(vec![NodeId(0), NodeId(1), NodeId(2)], vec![EdgeId(0), EdgeId(1), EdgeId(2)]),
        FaceBoundary::new(vec![NodeId(0), NodeId(2), NodeId(1)], vec![EdgeId(2), EdgeId(1), EdgeId(0)]),
    ],
    2,
    PlanarEmbeddingError::OuterFaceOutOfBounds { outer_face: 2, face_count: 2 },
)]
#[case::disconnected(
    Graph::new(4, &[[0, 1], [2, 3]]),
    vec![FaceBoundary::new(vec![NodeId(0), NodeId(1)], vec![EdgeId(0), EdgeId(0)])],
    0,
    PlanarEmbeddingError::Disconnected,
)]
#[case::face_too_short(
    Graph::new(3, &[[0, 1], [1, 2], [2, 0]]),
    vec![FaceBoundary::new(vec![NodeId(0)], vec![EdgeId(0)])],
    0,
    PlanarEmbeddingError::FaceTooShort { face: 0, length: 1 },
)]
#[case::boundary_length(
    Graph::new(3, &[[0, 1], [1, 2], [2, 0]]),
    vec![FaceBoundary::new(
        vec![NodeId(0), NodeId(1), NodeId(2)],
        vec![EdgeId(0), EdgeId(1)],
    )],
    0,
    PlanarEmbeddingError::BoundaryLengthMismatch { face: 0, node_count: 3, edge_count: 2 },
)]
#[case::node_bound(
    Graph::new(3, &[[0, 1], [1, 2], [2, 0]]),
    vec![FaceBoundary::new(vec![NodeId(0), NodeId(3)], vec![EdgeId(0), EdgeId(0)])],
    0,
    PlanarEmbeddingError::NodeOutOfBounds {
        face: 0, position: 1, node: NodeId(3), node_count: 3,
    },
)]
#[case::edge_bound(
    Graph::new(3, &[[0, 1], [1, 2], [2, 0]]),
    vec![FaceBoundary::new(vec![NodeId(0), NodeId(1)], vec![EdgeId(9), EdgeId(9)])],
    0,
    PlanarEmbeddingError::EdgeOutOfBounds {
        face: 0, position: 0, edge: EdgeId(9), edge_count: 3,
    },
)]
#[case::endpoint(
    Graph::new(3, &[[0, 1], [1, 2], [2, 0]]),
    vec![FaceBoundary::new(vec![NodeId(0), NodeId(1)], vec![EdgeId(1), EdgeId(0)])],
    0,
    PlanarEmbeddingError::EdgeEndpointMismatch {
        face: 0, position: 0, edge: EdgeId(1), first: NodeId(0), second: NodeId(1),
    },
)]
#[case::incidence(
    Graph::new(3, &[[0, 1], [1, 2], [2, 0]]),
    vec![FaceBoundary::new(
        vec![NodeId(0), NodeId(1), NodeId(2)],
        vec![EdgeId(0), EdgeId(1), EdgeId(2)],
    )],
    0,
    PlanarEmbeddingError::EdgeIncidence { edge: EdgeId(0), count: 1 },
)]
#[case::orientation(
    Graph::new(3, &[[0, 1], [1, 2], [2, 0]]),
    vec![
        FaceBoundary::new(vec![NodeId(0), NodeId(1), NodeId(2)], vec![EdgeId(0), EdgeId(1), EdgeId(2)]),
        FaceBoundary::new(vec![NodeId(0), NodeId(1), NodeId(2)], vec![EdgeId(0), EdgeId(1), EdgeId(2)]),
    ],
    0,
    PlanarEmbeddingError::InconsistentOrientation { edge: EdgeId(0) },
)]
#[case::noncellular(
    Graph::new(
        6,
        &[[0, 2], [2, 1], [1, 3], [3, 0], [0, 4], [4, 1], [1, 5], [5, 0]],
    ),
    vec![
        FaceBoundary::new(
            vec![NodeId(0), NodeId(2), NodeId(1), NodeId(3)],
            vec![EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3)],
        ),
        FaceBoundary::new(
            vec![NodeId(0), NodeId(3), NodeId(1), NodeId(2)],
            vec![EdgeId(3), EdgeId(2), EdgeId(1), EdgeId(0)],
        ),
        FaceBoundary::new(
            vec![NodeId(0), NodeId(4), NodeId(1), NodeId(5)],
            vec![EdgeId(4), EdgeId(5), EdgeId(6), EdgeId(7)],
        ),
        FaceBoundary::new(
            vec![NodeId(0), NodeId(5), NodeId(1), NodeId(4)],
            vec![EdgeId(7), EdgeId(6), EdgeId(5), EdgeId(4)],
        ),
    ],
    0,
    PlanarEmbeddingError::NonCellularVertex { node: NodeId(0) },
)]
#[case::euler(
    Graph::new(3, &[[0, 1], [1, 2], [2, 0]]),
    vec![FaceBoundary::new(
        vec![NodeId(0), NodeId(1), NodeId(2), NodeId(0), NodeId(2), NodeId(1)],
        vec![EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(2), EdgeId(1), EdgeId(0)],
    )],
    0,
    PlanarEmbeddingError::EulerCharacteristic { characteristic: 1 },
)]
fn test_planar_embedding_new_error(
    #[case] graph: Graph,
    #[case] faces: Vec<FaceBoundary>,
    #[case] outer_face: usize,
    #[case] expected: PlanarEmbeddingError,
) {
    assert_eq!(
        PlanarEmbedding::new(&graph, faces, outer_face),
        Err(expected)
    );
}
