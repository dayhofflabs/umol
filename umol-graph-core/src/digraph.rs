//! CSR directed graph (out-adjacency only).
//!
//! Stores per-node successor lists in compressed sparse row form; node data lives
//! externally, indexed by `NodeId`. Minimal by design — predecessors, edge identity,
//! and copy-on-write are added when a directed algorithm needs them.

use crate::graph::NodeId;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DiGraph {
    offsets: Vec<u32>,
    targets: Vec<NodeId>,
    node_count: usize,
}

impl DiGraph {
    /// Build from directed `[from, to]` edges. Successor order within a node follows
    /// edge input order; parallel edges are kept.
    pub fn new(node_count: usize, edges: &[[u32; 2]]) -> Self {
        let mut out_degree = vec![0u32; node_count];
        for &[from, _] in edges {
            out_degree[from as usize] += 1;
        }
        let mut offsets = vec![0u32; node_count + 1];
        for i in 0..node_count {
            offsets[i + 1] = offsets[i] + out_degree[i];
        }
        let mut targets = vec![NodeId(0); edges.len()];
        let mut cursor: Vec<u32> = offsets[..node_count].to_vec();
        for &[from, to] in edges {
            let pos = cursor[from as usize] as usize;
            targets[pos] = NodeId(to);
            cursor[from as usize] += 1;
        }
        Self {
            offsets,
            targets,
            node_count,
        }
    }

    pub fn node_count(&self) -> usize {
        self.node_count
    }

    pub fn node_ids(&self) -> impl Iterator<Item = NodeId> {
        (0..self.node_count as u32).map(NodeId)
    }

    /// The out-neighbors of `id`.
    pub fn successors(&self, id: NodeId) -> &[NodeId] {
        let start = self.offsets[id.index()] as usize;
        let end = self.offsets[id.index() + 1] as usize;
        &self.targets[start..end]
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;

    use super::*;

    #[rstest]
    #[case(NodeId(0), vec![NodeId(1), NodeId(2)])]
    #[case(NodeId(1), vec![NodeId(2)])]
    #[case(NodeId(2), vec![])]
    fn test_digraph_successors(#[case] node: NodeId, #[case] expected: Vec<NodeId>) {
        let g = DiGraph::new(3, &[[0, 1], [0, 2], [1, 2]]);
        assert_eq!(g.successors(node), expected);
    }
}
