//! Explicit planar embeddings for planar perfect-matching algorithms.
//!
//! This module validates an embedding supplied by the caller. It does not test
//! planarity or discover an embedding.

use std::collections::VecDeque;
use std::error::Error;
use std::fmt;

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
