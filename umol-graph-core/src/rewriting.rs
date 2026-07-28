//! Double-pushout (DPO) structural primitives over `Graph` + `GraphCorrespondence`.
//!
//! These are the chemistry-agnostic constructions of algebraic graph transformation (Ehrig, Ehrig,
//! Prange & Taentzer, *Fundamentals of Algebraic Graph Transformation*, Springer 2006): the
//! **pushout** (glue two graphs along an overlap, Def. 2.16), the **pushout complement** (delete a
//! matched subgraph keeping its context, Def. 9.8 / Fact 9.9), and the **pullback** (the shared
//! preimage of two morphisms into a common graph, Def. 2.22). They are **canonical** — each is unique
//! up to isomorphism — so none takes an algorithm choice. Attributes and overlays are not part of the
//! `Graph`; the attributed `meet`-glue and rule application live one layer up (umol-ast) and carry
//! their data through the morphisms these return.
//!
//! A morphism is a [`GraphCorrespondence`]; a rule/overlap is one too (its matched pairs are the
//! shared interface, its unmatched ids the deleted / created part). All graphs are **simple**: a
//! pushout that would induce a parallel edge identifies it instead (the pushout in the category of
//! simple graphs).

use std::collections::HashMap;

use crate::correspondence::{Correspondence, GraphCorrespondence};
use crate::graph::{EdgeId, Graph, NodeId};

/// The glued graph and the two coprojections into it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Pushout {
    pub object: Graph,
    /// `left → object` (the identity inclusion — `object` keeps `left`'s ids).
    pub left: GraphCorrespondence,
    /// `right → object` (overlap ids fold onto their `left` partners; the rest are appended).
    pub right: GraphCorrespondence,
}

/// The context graph `D` of a matched deletion and its embeddings.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PushoutComplement {
    pub object: Graph,
    /// `D → host` (each context id back to the host it survived from).
    pub context: GraphCorrespondence,
    /// `K → D` (the preserved interface, embedded in the context).
    pub interface: GraphCorrespondence,
}

/// The shared-preimage graph and its two projections.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Pullback {
    pub object: Graph,
    /// `object → left`.
    pub left: GraphCorrespondence,
    /// `object → right`.
    pub right: GraphCorrespondence,
}

fn sorted([a, b]: [u32; 2]) -> [u32; 2] {
    if a <= b {
        [a, b]
    } else {
        [b, a]
    }
}

fn identity<Id: Copy + Ord + From<usize>>(count: usize) -> Vec<Id> {
    (0..count).map(Id::from).collect()
}

impl Graph {
    /// Glue `self` and `right` identifying the matched pairs of `overlap` (a partial `self ↔ right`
    /// correspondence): the pushout of the span it denotes (Def. 2.16). `self` keeps its ids; `right`'s
    /// unmatched nodes/edges are appended. An appended edge whose endpoints already carry one is
    /// identified with it (simple-graph pushout).
    pub fn pushout(&self, right: &Graph, overlap: &GraphCorrespondence) -> Pushout {
        let n_left = self.node_count();
        let m_left = self.edge_count();

        // right node → object node: matched fold onto their left partner, unmatched append after `left`.
        let unmatched_nodes = overlap.nodes().right_unmatched();
        let mut right_node: HashMap<NodeId, NodeId> = overlap
            .nodes()
            .matched_pairs()
            .iter()
            .map(|&(l, r)| (r, l))
            .collect();
        for (rank, &r) in unmatched_nodes.iter().enumerate() {
            right_node.insert(r, NodeId::from(n_left + rank));
        }
        let object_node_count = n_left + unmatched_nodes.len();

        // Object edges: `left`'s verbatim, then `right`'s unmatched (endpoints relabelled),
        // collapsing a parallel onto the existing edge.
        let mut edges: Vec<[u32; 2]> = Vec::with_capacity(m_left);
        let mut by_pair: HashMap<[u32; 2], EdgeId> = HashMap::new();
        for j in 0..m_left {
            let [a, b] = self.edge_endpoints(EdgeId::from(j));
            by_pair.insert(sorted([a.0, b.0]), EdgeId::from(j));
            edges.push([a.0, b.0]);
        }
        let mut right_edge: HashMap<EdgeId, EdgeId> = overlap
            .edges()
            .matched_pairs()
            .iter()
            .map(|&(le, re)| (re, le))
            .collect();
        for re in overlap.edges().right_unmatched() {
            let [u, v] = right.edge_endpoints(re);
            let ends = [right_node[&u].0, right_node[&v].0];
            let object_edge = *by_pair.entry(sorted(ends)).or_insert_with(|| {
                let id = EdgeId::from(edges.len());
                edges.push(ends);
                id
            });
            right_edge.insert(re, object_edge);
        }

        let object = Graph::new(object_node_count, &edges);
        let object_edge_count = edges.len();

        let left_map = GraphCorrespondence::new(
            Correspondence::from_images(&identity::<NodeId>(n_left), object_node_count),
            Correspondence::from_images(&identity::<EdgeId>(m_left), object_edge_count),
        );
        let right_node_images: Vec<NodeId> = (0..right.node_count())
            .map(|r| right_node[&NodeId::from(r)])
            .collect();
        let right_edge_images: Vec<EdgeId> = (0..right.edge_count())
            .map(|e| right_edge[&EdgeId::from(e)])
            .collect();
        let right_map = GraphCorrespondence::new(
            Correspondence::from_images(&right_node_images, object_node_count),
            Correspondence::from_images(&right_edge_images, object_edge_count),
        );

        Pushout {
            object,
            left: left_map,
            right: right_map,
        }
    }

    /// Delete `matched`(`L\K`) from `self`, keeping the context — the pushout complement of the left
    /// DPO square (Def. 9.8, Fact 9.9). `matched` is the match `L → self`, `interface` the preserved
    /// `K → L`. `None` when the gluing (dangling) condition fails.
    pub fn pushout_complement(
        &self,
        matched: &GraphCorrespondence,
        interface: &GraphCorrespondence,
    ) -> Option<PushoutComplement> {
        // The deleted host items: L\K (interface-unmatched on the L side), carried through the match.
        let deleted_nodes: Vec<NodeId> = interface
            .nodes()
            .right_unmatched()
            .into_iter()
            .filter_map(|l| matched.nodes().right_of(l))
            .collect();
        let deleted_edges: Vec<EdgeId> = interface
            .edges()
            .right_unmatched()
            .into_iter()
            .filter_map(|le| matched.edges().right_of(le))
            .collect();

        let mut object = self.clone();
        let compaction = object.try_remove(&deleted_nodes, &deleted_edges)?;

        let host_to_d_nodes = Correspondence::new(
            (0..self.node_count())
                .filter_map(|h| {
                    let host_node = NodeId::from(h);
                    compaction.compact_node(host_node).map(|d| (host_node, d))
                })
                .collect(),
            self.node_count(),
            object.node_count(),
        );
        let host_to_d_edges = Correspondence::new(
            (0..self.edge_count())
                .filter_map(|e| {
                    let host_edge = EdgeId::from(e);
                    compaction.compact_edge(host_edge).map(|d| (host_edge, d))
                })
                .collect(),
            self.edge_count(),
            object.edge_count(),
        );

        let context =
            GraphCorrespondence::new(host_to_d_nodes.reverse(), host_to_d_edges.reverse());
        let interface_to_d = GraphCorrespondence::new(
            interface
                .nodes()
                .compose(matched.nodes())
                .compose(&host_to_d_nodes),
            interface
                .edges()
                .compose(matched.edges())
                .compose(&host_to_d_edges),
        );

        Some(PushoutComplement {
            object,
            context,
            interface: interface_to_d,
        })
    }

    /// The shared preimage of `self → E` and `right → E` — the pullback of the cospan (Def. 2.22): the
    /// largest subgraph of both that maps consistently into `E`, with its two projections. Used to
    /// build a composite rule's interface (`K = C₁ ×_E C₂`, Def. 9.25).
    pub fn pullback(
        &self,
        right: &Graph,
        left_into: &GraphCorrespondence,
        right_into: &GraphCorrespondence,
    ) -> Pullback {
        // self ↔ right over the common E: match a self item to the right item sharing its E-image.
        let node_matched_pairs = left_into
            .nodes()
            .compose(&right_into.nodes().reverse())
            .matched_pairs()
            .to_vec();
        let edge_matched_pairs = left_into
            .edges()
            .compose(&right_into.edges().reverse())
            .matched_pairs()
            .to_vec();

        let left_to_k: HashMap<NodeId, NodeId> = node_matched_pairs
            .iter()
            .enumerate()
            .map(|(i, &(l, _))| (l, NodeId::from(i)))
            .collect();
        let edges: Vec<[u32; 2]> = edge_matched_pairs
            .iter()
            .map(|&(le, _)| {
                let [a, b] = self.edge_endpoints(le);
                [left_to_k[&a].0, left_to_k[&b].0]
            })
            .collect();
        let object = Graph::new(node_matched_pairs.len(), &edges);

        let left_map = GraphCorrespondence::new(
            Correspondence::from_images(
                &node_matched_pairs
                    .iter()
                    .map(|&(l, _)| l)
                    .collect::<Vec<_>>(),
                self.node_count(),
            ),
            Correspondence::from_images(
                &edge_matched_pairs
                    .iter()
                    .map(|&(le, _)| le)
                    .collect::<Vec<_>>(),
                self.edge_count(),
            ),
        );
        let right_map = GraphCorrespondence::new(
            Correspondence::from_images(
                &node_matched_pairs
                    .iter()
                    .map(|&(_, r)| r)
                    .collect::<Vec<_>>(),
                right.node_count(),
            ),
            Correspondence::from_images(
                &edge_matched_pairs
                    .iter()
                    .map(|&(_, re)| re)
                    .collect::<Vec<_>>(),
                right.edge_count(),
            ),
        );

        Pullback {
            object,
            left: left_map,
            right: right_map,
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    fn node_overlap(
        matched_pairs: Vec<(u32, u32)>,
        left: &Graph,
        right: &Graph,
    ) -> GraphCorrespondence {
        let nodes = Correspondence::new(
            matched_pairs
                .iter()
                .map(|&(l, r)| (NodeId(l), NodeId(r)))
                .collect(),
            left.node_count(),
            right.node_count(),
        );
        GraphCorrespondence::induced(left, right, nodes)
    }

    #[rstest]
    fn test_pushout_glue_at_node() {
        // path 0-1 glued to path 0-1 identifying left node 1 with right node 0 → path 0-1-2.
        let left = Graph::new(2, &[[0, 1]]);
        let right = Graph::new(2, &[[0, 1]]);
        let overlap = node_overlap(vec![(1, 0)], &left, &right);
        let po = left.pushout(&right, &overlap);
        assert_eq!(po.object.node_count(), 3);
        assert_eq!(po.object.edge_count(), 2);
        assert_eq!(po.object.edge_endpoints(EdgeId(1)), [NodeId(1), NodeId(2)]);
        assert_eq!(po.right.nodes().right_of(NodeId(0)), Some(NodeId(1)));
        assert_eq!(po.right.nodes().right_of(NodeId(1)), Some(NodeId(2)));
        assert_eq!(po.right.edges().right_of(EdgeId(0)), Some(EdgeId(1)));
    }

    #[rstest]
    fn test_pushout_keeps_context_edge() {
        // doc-135 R1: glue R_A = F–Cl (with bond) and L_B = [F, Cl] (no bond) over both atoms.
        // The overlap omits the bond (L_B lacks it), yet the glue keeps R_A's bond as context.
        let r_a = Graph::new(2, &[[0, 1]]);
        let l_b = Graph::new(2, &[]);
        let overlap = node_overlap(vec![(0, 0), (1, 1)], &r_a, &l_b);
        let po = r_a.pushout(&l_b, &overlap);
        assert_eq!(po.object.node_count(), 2);
        assert_eq!(po.object.edge_count(), 1);
        assert_eq!(po.object.edge_endpoints(EdgeId(0)), [NodeId(0), NodeId(1)]);
        assert_eq!(po.right.nodes().right_of(NodeId(1)), Some(NodeId(1)));
    }

    #[rstest]
    fn test_pushout_complement_deletes_matched() {
        // host path 0-1-2-3; delete node 1 with its two incident edges, keep endpoints 0 and 2.
        let host = Graph::new(4, &[[0, 1], [1, 2], [2, 3]]);
        // L = path 0-1-2 matched onto host 0-1-2.
        let matched = GraphCorrespondence::induced(
            &Graph::new(3, &[[0, 1], [1, 2]]),
            &host,
            Correspondence::from_images(&[NodeId(0), NodeId(1), NodeId(2)], 4),
        );
        // K = the two endpoints of L (nodes 0 and 2), no edges — the rest is deleted.
        let interface = GraphCorrespondence::new(
            Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(2))], 2, 3),
            Correspondence::new(vec![], 0, 2),
        );
        let pc = host
            .pushout_complement(&matched, &interface)
            .expect("no dangling");
        assert_eq!(pc.object.node_count(), 3);
        assert_eq!(pc.object.edge_count(), 1);
        // host node 3 survives as D node 2; the K endpoint (host node 2) as D node 1.
        assert_eq!(pc.context.nodes().right_of(NodeId(2)), Some(NodeId(3)));
        assert_eq!(pc.interface.nodes().right_of(NodeId(1)), Some(NodeId(1)));
    }

    #[rstest]
    fn test_pushout_complement_dangling() {
        // Same match, but K keeps the edges — deleting node 1 would strand them → None.
        let host = Graph::new(4, &[[0, 1], [1, 2], [2, 3]]);
        let matched = GraphCorrespondence::induced(
            &Graph::new(3, &[[0, 1], [1, 2]]),
            &host,
            Correspondence::from_images(&[NodeId(0), NodeId(1), NodeId(2)], 4),
        );
        // K keeps both endpoints and both edges; only L node 1 is deleted → it dangles.
        let interface = GraphCorrespondence::new(
            Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(2))], 2, 3),
            Correspondence::new(vec![(EdgeId(0), EdgeId(0)), (EdgeId(1), EdgeId(1))], 2, 2),
        );
        assert_eq!(host.pushout_complement(&matched, &interface), None);
    }

    #[rstest]
    fn test_pullback_shared_subgraph() {
        // left, right both path 0-1-2; E identifies left {0,1} with right {1,2}, left edge 0-1 with
        // right edge 1-2. Pullback is that shared edge.
        let left = Graph::new(3, &[[0, 1], [1, 2]]);
        let right = Graph::new(3, &[[0, 1], [1, 2]]);
        // left → E and right → E onto a 2-node/1-edge E.
        let left_into = GraphCorrespondence::new(
            Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))], 3, 2),
            Correspondence::new(vec![(EdgeId(0), EdgeId(0))], 2, 1),
        );
        let right_into = GraphCorrespondence::new(
            Correspondence::new(vec![(NodeId(1), NodeId(0)), (NodeId(2), NodeId(1))], 3, 2),
            Correspondence::new(vec![(EdgeId(1), EdgeId(0))], 2, 1),
        );
        let pb = left.pullback(&right, &left_into, &right_into);
        assert_eq!(pb.object.node_count(), 2);
        assert_eq!(pb.object.edge_count(), 1);
        assert_eq!(pb.left.nodes().right_of(NodeId(0)), Some(NodeId(0)));
        assert_eq!(pb.right.nodes().right_of(NodeId(0)), Some(NodeId(1)));
        assert_eq!(pb.left.edges().right_of(EdgeId(0)), Some(EdgeId(0)));
        assert_eq!(pb.right.edges().right_of(EdgeId(0)), Some(EdgeId(1)));
    }
}
