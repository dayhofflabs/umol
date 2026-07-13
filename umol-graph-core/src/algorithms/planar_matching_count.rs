//! Explicit planar embeddings for planar perfect-matching algorithms.
//!
//! This module validates an embedding supplied by the caller. It does not test
//! planarity or discover an embedding.

use std::collections::VecDeque;
use std::error::Error;
use std::fmt;

use num_bigint::{BigInt, BigUint};

use crate::{EdgeId, Graph, NodeId};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FaceBoundary {
    nodes: Vec<NodeId>,
    edges: Vec<EdgeId>,
}

impl FaceBoundary {
    pub fn new(nodes: Vec<NodeId>, edges: Vec<EdgeId>) -> Self {
        Self { nodes, edges }
    }

    pub fn nodes(&self) -> &[NodeId] {
        &self.nodes
    }

    pub fn edges(&self) -> &[EdgeId] {
        &self.edges
    }

    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PlanarEmbedding {
    graph: Graph,
    faces: Vec<FaceBoundary>,
    outer_face: usize,
}

impl PlanarEmbedding {
    pub fn new(
        graph: &Graph,
        faces: Vec<FaceBoundary>,
        outer_face: usize,
    ) -> Result<Self, PlanarEmbeddingError> {
        if outer_face >= faces.len() {
            return Err(PlanarEmbeddingError::OuterFaceOutOfBounds {
                outer_face,
                face_count: faces.len(),
            });
        }
        if graph.node_count() == 0 {
            return Err(PlanarEmbeddingError::Disconnected);
        }

        let mut visited = vec![false; graph.node_count()];
        let mut queue = VecDeque::from([NodeId(0)]);
        visited[0] = true;
        while let Some(node) = queue.pop_front() {
            for neighbor in graph.neighbors(node) {
                if !visited[neighbor.node.index()] {
                    visited[neighbor.node.index()] = true;
                    queue.push_back(neighbor.node);
                }
            }
        }
        if visited.iter().any(|&seen| !seen) {
            return Err(PlanarEmbeddingError::Disconnected);
        }

        let mut occurrences: Vec<Vec<(NodeId, NodeId)>> = vec![Vec::new(); graph.edge_count()];
        let mut vertex_links = vec![vec![Vec::new(); graph.edge_count()]; graph.node_count()];
        for (face_index, face) in faces.iter().enumerate() {
            if face.nodes.len() < 2 {
                return Err(PlanarEmbeddingError::FaceTooShort {
                    face: face_index,
                    length: face.nodes.len(),
                });
            }
            if face.nodes.len() != face.edges.len() {
                return Err(PlanarEmbeddingError::BoundaryLengthMismatch {
                    face: face_index,
                    node_count: face.nodes.len(),
                    edge_count: face.edges.len(),
                });
            }
            for position in 0..face.nodes.len() {
                let first = face.nodes[position];
                let second = face.nodes[(position + 1) % face.nodes.len()];
                let edge = face.edges[position];
                if !graph.contains_node(first) {
                    return Err(PlanarEmbeddingError::NodeOutOfBounds {
                        face: face_index,
                        position,
                        node: first,
                        node_count: graph.node_count(),
                    });
                }
                if !graph.contains_node(second) {
                    return Err(PlanarEmbeddingError::NodeOutOfBounds {
                        face: face_index,
                        position: (position + 1) % face.nodes.len(),
                        node: second,
                        node_count: graph.node_count(),
                    });
                }
                if !graph.contains_edge(edge) {
                    return Err(PlanarEmbeddingError::EdgeOutOfBounds {
                        face: face_index,
                        position,
                        edge,
                        edge_count: graph.edge_count(),
                    });
                }
                let [edge_first, edge_second] = graph.edge_endpoints(edge);
                if !((edge_first == first && edge_second == second)
                    || (edge_first == second && edge_second == first))
                {
                    return Err(PlanarEmbeddingError::EdgeEndpointMismatch {
                        face: face_index,
                        position,
                        edge,
                        first,
                        second,
                    });
                }
                occurrences[edge.index()].push((first, second));

                let incoming = face.edges[(position + face.edges.len() - 1) % face.edges.len()];
                if !graph.contains_edge(incoming) {
                    return Err(PlanarEmbeddingError::EdgeOutOfBounds {
                        face: face_index,
                        position: (position + face.edges.len() - 1) % face.edges.len(),
                        edge: incoming,
                        edge_count: graph.edge_count(),
                    });
                }
                vertex_links[first.index()][incoming.index()].push(edge);
            }
        }

        for (edge_index, directions) in occurrences.iter().enumerate() {
            let edge = EdgeId(edge_index as u32);
            if directions.len() != 2 {
                return Err(PlanarEmbeddingError::EdgeIncidence {
                    edge,
                    count: directions.len(),
                });
            }
            if directions[0] != (directions[1].1, directions[1].0) {
                return Err(PlanarEmbeddingError::InconsistentOrientation { edge });
            }
        }

        let characteristic =
            graph.node_count() as isize - graph.edge_count() as isize + faces.len() as isize;
        if characteristic != 2 {
            return Err(PlanarEmbeddingError::EulerCharacteristic { characteristic });
        }

        for node in graph.node_ids() {
            let mut incident: Vec<_> = graph
                .neighbors(node)
                .iter()
                .map(|neighbor| neighbor.edge)
                .collect();
            incident.sort_unstable();
            incident.dedup();
            let Some(&start) = incident.first() else {
                continue;
            };
            let mut seen = vec![false; graph.edge_count()];
            let mut current = start;
            for _ in 0..incident.len() {
                if seen[current.index()] {
                    return Err(PlanarEmbeddingError::NonCellularVertex { node });
                }
                seen[current.index()] = true;
                let next_edges = &vertex_links[node.index()][current.index()];
                if next_edges.len() != 1 {
                    return Err(PlanarEmbeddingError::NonCellularVertex { node });
                }
                current = next_edges[0];
            }
            if current != start || incident.iter().any(|edge| !seen[edge.index()]) {
                return Err(PlanarEmbeddingError::NonCellularVertex { node });
            }
        }

        Ok(Self {
            graph: graph.clone(),
            faces,
            outer_face,
        })
    }

    pub fn graph(&self) -> &Graph {
        &self.graph
    }

    pub fn faces(&self) -> &[FaceBoundary] {
        &self.faces
    }

    pub fn outer_face_index(&self) -> usize {
        self.outer_face
    }

    pub fn outer_face(&self) -> &FaceBoundary {
        &self.faces[self.outer_face]
    }

    pub fn bounded_faces(&self) -> impl Iterator<Item = &FaceBoundary> {
        self.faces
            .iter()
            .enumerate()
            .filter_map(|(index, face)| (index != self.outer_face).then_some(face))
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PlanarEmbeddingError {
    OuterFaceOutOfBounds {
        outer_face: usize,
        face_count: usize,
    },
    Disconnected,
    FaceTooShort {
        face: usize,
        length: usize,
    },
    BoundaryLengthMismatch {
        face: usize,
        node_count: usize,
        edge_count: usize,
    },
    NodeOutOfBounds {
        face: usize,
        position: usize,
        node: NodeId,
        node_count: usize,
    },
    EdgeOutOfBounds {
        face: usize,
        position: usize,
        edge: EdgeId,
        edge_count: usize,
    },
    EdgeEndpointMismatch {
        face: usize,
        position: usize,
        edge: EdgeId,
        first: NodeId,
        second: NodeId,
    },
    EdgeIncidence {
        edge: EdgeId,
        count: usize,
    },
    InconsistentOrientation {
        edge: EdgeId,
    },
    NonCellularVertex {
        node: NodeId,
    },
    EulerCharacteristic {
        characteristic: isize,
    },
}

impl fmt::Display for PlanarEmbeddingError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "invalid planar embedding: {self:?}")
    }
}

impl Error for PlanarEmbeddingError {}

impl Graph {
    pub fn count_perfect_matchings_planar(
        &self,
        embedding: &PlanarEmbedding,
    ) -> Result<BigUint, PlanarMatchingCountError> {
        if self != embedding.graph() {
            return Err(PlanarMatchingCountError::EmbeddingGraphMismatch);
        }
        if self.node_count() % 2 == 1 {
            return Ok(BigUint::from(0_u8));
        }

        let signs = kasteleyn_signs(embedding)?;
        let mut matrix = vec![vec![BigInt::from(0); self.node_count()]; self.node_count()];
        for edge in self.edge_ids() {
            let [first, second] = self.edge_endpoints(edge);
            if first == second {
                continue;
            }
            let (lower, upper) = if first < second {
                (first, second)
            } else {
                (second, first)
            };
            let value = if signs[edge.index()] {
                BigInt::from(-1)
            } else {
                BigInt::from(1)
            };
            matrix[lower.index()][upper.index()] += &value;
            matrix[upper.index()][lower.index()] -= value;
        }

        let (_, count) = pfaffian(&matrix)?.into_parts();
        Ok(count)
    }
}

fn kasteleyn_signs(embedding: &PlanarEmbedding) -> Result<Vec<bool>, PlanarMatchingCountError> {
    let edge_count = embedding.graph.edge_count();
    let mut equations = Vec::with_capacity(embedding.faces.len().saturating_sub(1));
    for (face_index, face) in embedding.faces.iter().enumerate() {
        if face_index == embedding.outer_face {
            continue;
        }
        let mut equation = vec![false; edge_count + 1];
        let mut negative = face.len() % 2 == 0;
        for position in 0..face.len() {
            let first = face.nodes[position];
            let second = face.nodes[(position + 1) % face.len()];
            equation[face.edges[position].index()] ^= true;
            negative ^= first > second;
        }
        equation[edge_count] = negative;
        equations.push(equation);
    }

    let mut pivot_row = 0;
    let mut pivots = Vec::new();
    for column in 0..edge_count {
        let Some(found) = (pivot_row..equations.len()).find(|&row| equations[row][column]) else {
            continue;
        };
        equations.swap(pivot_row, found);
        let pivot_equation = equations[pivot_row].clone();
        for (row, equation) in equations.iter_mut().enumerate() {
            if row != pivot_row && equation[column] {
                for (entry, &pivot_entry) in
                    equation[column..].iter_mut().zip(&pivot_equation[column..])
                {
                    *entry ^= pivot_entry;
                }
            }
        }
        pivots.push((pivot_row, column));
        pivot_row += 1;
        if pivot_row == equations.len() {
            break;
        }
    }
    if equations
        .iter()
        .any(|equation| !equation[..edge_count].iter().any(|&entry| entry) && equation[edge_count])
    {
        return Err(PlanarMatchingCountError::InconsistentSigning);
    }

    let mut signs = vec![false; edge_count];
    for (row, column) in pivots {
        signs[column] = equations[row][edge_count];
    }
    Ok(signs)
}

fn pfaffian(matrix: &[Vec<BigInt>]) -> Result<BigInt, PlanarMatchingCountError> {
    let size = matrix.len();
    if matrix.iter().any(|row| row.len() != size) {
        return Err(PlanarMatchingCountError::NonSquareMatrix);
    }
    for (row, values) in matrix.iter().enumerate() {
        if values[row] != BigInt::from(0) {
            return Err(PlanarMatchingCountError::NonSkewSymmetric);
        }
        for (column, other) in matrix.iter().enumerate().skip(row + 1) {
            if values[column] != -&other[row] {
                return Err(PlanarMatchingCountError::NonSkewSymmetric);
            }
        }
    }
    if size == 0 {
        return Ok(BigInt::from(1));
    }
    if size % 2 == 1 {
        return Ok(BigInt::from(0));
    }

    let mut work = matrix.to_vec();
    let mut permutation_sign = BigInt::from(1);
    let mut previous_pivot = BigInt::from(1);
    for first in (0..size).step_by(2) {
        let Some(pivot_column) =
            (first + 1..size).find(|&column| work[first][column] != BigInt::from(0))
        else {
            return Ok(BigInt::from(0));
        };
        if pivot_column != first + 1 {
            work.swap(first + 1, pivot_column);
            for row in &mut work {
                row.swap(first + 1, pivot_column);
            }
            permutation_sign = -permutation_sign;
        }

        let pivot = work[first][first + 1].clone();
        if first + 2 == size {
            return Ok(permutation_sign * pivot);
        }
        for row in first + 2..size {
            for column in row + 1..size {
                let numerator = &pivot * &work[row][column]
                    - &work[first][row] * &work[first + 1][column]
                    + &work[first][column] * &work[first + 1][row];
                if &numerator % &previous_pivot != BigInt::from(0) {
                    return Err(PlanarMatchingCountError::InexactDivision { step: first });
                }
                let value = numerator / &previous_pivot;
                work[row][column] = value.clone();
                work[column][row] = -value;
            }
        }
        previous_pivot = pivot;
    }
    unreachable!("even nonempty matrix returns from its final pivot")
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PlanarMatchingCountError {
    EmbeddingGraphMismatch,
    InconsistentSigning,
    NonSquareMatrix,
    NonSkewSymmetric,
    InexactDivision { step: usize },
}

impl fmt::Display for PlanarMatchingCountError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "planar matching count failed: {self:?}")
    }
}

impl Error for PlanarMatchingCountError {}

#[cfg(test)]
mod tests {
    use num_bigint::{BigInt, BigUint};
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::empty(vec![], 1)]
    #[case::pair(vec![vec![0, 3], vec![-3, 0]], 3)]
    #[case::four(
        vec![
            vec![0, 1, 2, 3],
            vec![-1, 0, 4, 5],
            vec![-2, -4, 0, 6],
            vec![-3, -5, -6, 0],
        ],
        8,
    )]
    #[case::permuted(
        vec![
            vec![0, -1, 4, 5],
            vec![1, 0, 2, 3],
            vec![-4, -2, 0, 6],
            vec![-5, -3, -6, 0],
        ],
        -8,
    )]
    #[case::pivoted(
        vec![
            vec![0, 0, 2, 0],
            vec![0, 0, 0, 3],
            vec![-2, 0, 0, 0],
            vec![0, -3, 0, 0],
        ],
        -6,
    )]
    #[case::six(
        vec![
            vec![0, 2, 0, 0, 0, 0],
            vec![-2, 0, 0, 0, 0, 0],
            vec![0, 0, 0, 3, 0, 0],
            vec![0, 0, -3, 0, 0, 0],
            vec![0, 0, 0, 0, 0, 5],
            vec![0, 0, 0, 0, -5, 0],
        ],
        30,
    )]
    #[case::zero(vec![vec![0; 4]; 4], 0)]
    #[case::odd(vec![vec![0; 3]; 3], 0)]
    fn test_pfaffian(#[case] matrix: Vec<Vec<i64>>, #[case] expected: i64) {
        let matrix: Vec<Vec<BigInt>> = matrix
            .into_iter()
            .map(|row| row.into_iter().map(BigInt::from).collect())
            .collect();
        assert_eq!(pfaffian(&matrix), Ok(BigInt::from(expected)));
    }

    #[rstest]
    #[case::nonsquare(
        vec![vec![0, 1], vec![-1]],
        PlanarMatchingCountError::NonSquareMatrix,
    )]
    #[case::diagonal(
        vec![vec![1, 0], vec![0, 0]],
        PlanarMatchingCountError::NonSkewSymmetric,
    )]
    #[case::asymmetric(
        vec![vec![0, 1], vec![1, 0]],
        PlanarMatchingCountError::NonSkewSymmetric,
    )]
    fn test_pfaffian_error(
        #[case] matrix: Vec<Vec<i64>>,
        #[case] expected: PlanarMatchingCountError,
    ) {
        let matrix: Vec<Vec<BigInt>> = matrix
            .into_iter()
            .map(|row| row.into_iter().map(BigInt::from).collect())
            .collect();
        assert_eq!(pfaffian(&matrix), Err(expected));
    }

    #[rstest]
    #[case::cycle(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
        vec![
            FaceBoundary::new(
                vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3)],
                vec![EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3)],
            ),
            FaceBoundary::new(
                vec![NodeId(0), NodeId(3), NodeId(2), NodeId(1)],
                vec![EdgeId(3), EdgeId(2), EdgeId(1), EdgeId(0)],
            ),
        ],
        1,
    )]
    #[case::k4(
        Graph::new(4, &[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]),
        vec![
            FaceBoundary::new(
                vec![NodeId(0), NodeId(2), NodeId(1)],
                vec![EdgeId(1), EdgeId(3), EdgeId(0)],
            ),
            FaceBoundary::new(
                vec![NodeId(0), NodeId(1), NodeId(3)],
                vec![EdgeId(0), EdgeId(4), EdgeId(2)],
            ),
            FaceBoundary::new(
                vec![NodeId(0), NodeId(3), NodeId(2)],
                vec![EdgeId(2), EdgeId(5), EdgeId(1)],
            ),
            FaceBoundary::new(
                vec![NodeId(1), NodeId(2), NodeId(3)],
                vec![EdgeId(3), EdgeId(5), EdgeId(4)],
            ),
        ],
        0,
    )]
    fn test_kasteleyn_signs(
        #[case] graph: Graph,
        #[case] faces: Vec<FaceBoundary>,
        #[case] outer_face: usize,
    ) {
        let embedding = PlanarEmbedding::new(&graph, faces, outer_face).unwrap();
        let signs = kasteleyn_signs(&embedding).unwrap();

        for face in embedding.bounded_faces() {
            let mut negative = false;
            for position in 0..face.len() {
                let first = face.nodes()[position];
                let second = face.nodes()[(position + 1) % face.len()];
                negative ^= first > second;
                negative ^= signs[face.edges()[position].index()];
            }
            assert_eq!(negative, face.len() % 2 == 0);
        }
        assert_eq!(kasteleyn_signs(&embedding), Ok(signs));
    }

    #[rstest]
    #[case::cycle(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
        vec![
            FaceBoundary::new(
                vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3)],
                vec![EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3)],
            ),
            FaceBoundary::new(
                vec![NodeId(0), NodeId(3), NodeId(2), NodeId(1)],
                vec![EdgeId(3), EdgeId(2), EdgeId(1), EdgeId(0)],
            ),
        ],
        1,
        2,
    )]
    #[case::bridge(
        Graph::new(2, &[[0, 1]]),
        vec![FaceBoundary::new(
            vec![NodeId(0), NodeId(1)],
            vec![EdgeId(0), EdgeId(0)],
        )],
        0,
        1,
    )]
    #[case::odd(
        Graph::new(3, &[[0, 1], [1, 2], [2, 0]]),
        vec![
            FaceBoundary::new(
                vec![NodeId(0), NodeId(1), NodeId(2)],
                vec![EdgeId(0), EdgeId(1), EdgeId(2)],
            ),
            FaceBoundary::new(
                vec![NodeId(0), NodeId(2), NodeId(1)],
                vec![EdgeId(2), EdgeId(1), EdgeId(0)],
            ),
        ],
        1,
        0,
    )]
    #[case::k4(
        Graph::new(4, &[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]),
        vec![
            FaceBoundary::new(vec![NodeId(0), NodeId(2), NodeId(1)], vec![EdgeId(1), EdgeId(3), EdgeId(0)]),
            FaceBoundary::new(vec![NodeId(0), NodeId(1), NodeId(3)], vec![EdgeId(0), EdgeId(4), EdgeId(2)]),
            FaceBoundary::new(vec![NodeId(0), NodeId(3), NodeId(2)], vec![EdgeId(2), EdgeId(5), EdgeId(1)]),
            FaceBoundary::new(vec![NodeId(1), NodeId(2), NodeId(3)], vec![EdgeId(3), EdgeId(5), EdgeId(4)]),
        ],
        0,
        3,
    )]
    fn test_graph_count_perfect_matchings_planar(
        #[case] graph: Graph,
        #[case] faces: Vec<FaceBoundary>,
        #[case] outer_face: usize,
        #[case] expected: u32,
    ) {
        let embedding = PlanarEmbedding::new(&graph, faces, outer_face).unwrap();
        assert_eq!(
            graph.count_perfect_matchings_planar(&embedding),
            Ok(BigUint::from(expected))
        );
    }

    #[rstest]
    fn test_graph_count_perfect_matchings_planar_error() {
        let embedded_graph = Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]);
        let embedding = PlanarEmbedding::new(
            &embedded_graph,
            vec![
                FaceBoundary::new(
                    vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3)],
                    vec![EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3)],
                ),
                FaceBoundary::new(
                    vec![NodeId(0), NodeId(3), NodeId(2), NodeId(1)],
                    vec![EdgeId(3), EdgeId(2), EdgeId(1), EdgeId(0)],
                ),
            ],
            1,
        )
        .unwrap();
        let other = Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 1]]);

        assert_eq!(
            other.count_perfect_matchings_planar(&embedding),
            Err(PlanarMatchingCountError::EmbeddingGraphMismatch)
        );
    }
}
