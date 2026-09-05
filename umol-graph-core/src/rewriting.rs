//! Double-pushout (DPO) structural primitives over `Graph` + `GraphCorrespondence`.
//!
//! These are the chemistry-agnostic constructions of algebraic graph transformation (Ehrig, Ehrig,
//! Prange & Taentzer, *Fundamentals of Algebraic Graph Transformation*, Springer 2006): the
//! **pushout** (glue two graphs along an overlap, Def. 2.16), the **pushout complement** (delete a
//! matched subgraph keeping its context, Def. 9.8 / Fact 9.9), and the **pullback** (the shared
//! preimage of two morphisms into a common graph, Def. 2.22). They are **canonical** — each is unique
//! up to isomorphism — so none takes an algorithm choice. Attributes and overlays are not part of the
//! `Graph`; the attributed `meet`-glue and rule application live one layer up (umol-graph-ir) and carry
//! their data through the morphisms these return.
//!
//! A morphism is a [`GraphCorrespondence`]; a rule/overlap is one too (its matched pairs are the
//! shared interface, its unmatched ids the deleted / created part). All graphs are **simple**: a
//! pushout that would induce a parallel edge identifies it instead (the pushout in the category of
//! simple graphs).

use std::collections::HashMap;

use crate::correspondence::{Correspondence, GraphCorrespondence, GraphCorrespondenceComposeError};
use crate::graph::{EdgeId, Graph, NodeId};

/// The two input-to-result correspondences of a graph pushout.
///
/// Operation-produced components have equal target counts and cover their respective
/// inputs. Public fields may be assembled independently; agreement with a particular
/// pushout graph is contextual.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GraphPushoutCorrespondence {
    /// `left → object` (the identity inclusion — `object` keeps `left`'s ids).
    pub left: GraphCorrespondence,
    /// `right → object` (overlap ids fold onto their `left` partners; the rest are appended).
    pub right: GraphCorrespondence,
}

/// The context-to-host and interface-to-context mappings of a pushout complement.
///
/// Operation-produced interface target counts equal the context mapping's source counts.
/// Public fields may be assembled independently; agreement with the context graph,
/// host, and interface is contextual.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PushoutComplementCorrespondence {
    /// `D → host` (each context id back to the host it survived from).
    pub context: GraphCorrespondence,
    /// `K → D` (the preserved interface, embedded in the context).
    pub interface: GraphCorrespondence,
}

/// The two result-to-input projections of a graph pullback.
///
/// Operation-produced components cover the result and have equal source counts. Public fields may be
/// assembled independently; agreement with the result and input graphs is contextual.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PullbackCorrespondence {
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
    ///
    /// Returns the glued graph. Use [`Self::tracked_pushout`] for its input mappings.
    pub fn pushout(&self, right: &Graph, overlap: &GraphCorrespondence) -> Graph {
        self.tracked_pushout(right, overlap).0
    }

    /// Glue two graphs and return both input-to-result correspondences with the result.
    ///
    /// Produces the same graph as [`Self::pushout`].
    pub fn tracked_pushout(
        &self,
        right: &Graph,
        overlap: &GraphCorrespondence,
    ) -> (Graph, GraphPushoutCorrespondence) {
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

        (
            object,
            GraphPushoutCorrespondence {
                left: left_map,
                right: right_map,
            },
        )
    }

    /// Delete `matched`(`L\K`) from `self`, keeping the context — the pushout complement of the left
    /// DPO square (Def. 9.8, Fact 9.9). `matched` is the match `L → self`, `interface` the preserved
    /// `K → L`. `None` when the gluing (dangling) condition fails or consecutive carrier counts disagree.
    /// Use [`Self::tracked_pushout_complement`] to retain the categorical mappings.
    pub fn pushout_complement(
        &self,
        matched: &GraphCorrespondence,
        interface: &GraphCorrespondence,
    ) -> Option<Graph> {
        self.tracked_pushout_complement(matched, interface)
            .map(|(object, _)| object)
    }

    /// Return the result and its categorical mappings.
    ///
    /// Has the same result and failure behavior as [`Self::pushout_complement`].
    pub fn tracked_pushout_complement(
        &self,
        matched: &GraphCorrespondence,
        interface: &GraphCorrespondence,
    ) -> Option<(Graph, PushoutComplementCorrespondence)> {
        let interface_to_host = interface.compose(matched).ok()?;
        if matched.nodes().right_count() != self.node_count()
            || matched.edges().right_count() != self.edge_count()
        {
            return None;
        }
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
        let compaction = object.try_tracked_remove(&deleted_nodes, &deleted_edges)?;

        let host_to_d_nodes = Correspondence::new(
            (0..self.node_count())
                .filter_map(|h| {
                    let host_node = NodeId::from(h);
                    compaction.compact_node(host_node).map(|d| (host_node, d))
                })
                .collect(),
            self.node_count(),
            object.node_count(),
        )
        .expect("graph compaction defines a valid node correspondence");
        let host_to_d_edges = Correspondence::new(
            (0..self.edge_count())
                .filter_map(|e| {
                    let host_edge = EdgeId::from(e);
                    compaction.compact_edge(host_edge).map(|d| (host_edge, d))
                })
                .collect(),
            self.edge_count(),
            object.edge_count(),
        )
        .expect("graph compaction defines a valid edge correspondence");

        let context =
            GraphCorrespondence::new(host_to_d_nodes.reverse(), host_to_d_edges.reverse());
        let interface_to_d = interface_to_host
            .compose(&GraphCorrespondence::new(host_to_d_nodes, host_to_d_edges))
            .expect("the checked match and compaction share the host counts");

        Some((
            object,
            PushoutComplementCorrespondence {
                context,
                interface: interface_to_d,
            },
        ))
    }

    /// The shared preimage of `self → E` and `right → E` — the pullback of the cospan (Def. 2.22): the
    /// largest subgraph of both that maps consistently into `E`, with its two projections. Used to
    /// build a composite rule's interface (`K = C₁ ×_E C₂`, Def. 9.25).
    /// Use [`Self::tracked_pullback`] to retain the projections with the result.
    ///
    /// # Errors
    /// Returns the component whose two target counts disagree.
    pub fn pullback(
        &self,
        right: &Graph,
        left_into: &GraphCorrespondence,
        right_into: &GraphCorrespondence,
    ) -> Result<Graph, GraphCorrespondenceComposeError> {
        self.tracked_pullback(right, left_into, right_into)
            .map(|(object, _)| object)
    }

    /// Return the result and its categorical mappings.
    ///
    /// Has the same result and failure behavior as [`Self::pullback`].
    ///
    /// # Errors
    ///
    /// Returns the component whose two target counts disagree.
    pub fn tracked_pullback(
        &self,
        right: &Graph,
        left_into: &GraphCorrespondence,
        right_into: &GraphCorrespondence,
    ) -> Result<(Graph, PullbackCorrespondence), GraphCorrespondenceComposeError> {
        // self ↔ right over the common E: match a self item to the right item sharing its E-image.
        let correspondence = left_into.compose(&GraphCorrespondence::new(
            right_into.nodes().reverse(),
            right_into.edges().reverse(),
        ))?;
        let node_matched_pairs = correspondence.nodes().matched_pairs().to_vec();
        let edge_matched_pairs = correspondence.edges().matched_pairs().to_vec();

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

        Ok((
            object,
            PullbackCorrespondence {
                left: left_map,
                right: right_map,
            },
        ))
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::correspondence::CorrespondenceComposeError;

    #[rstest]
    fn test_graph_pushout_partial() {
        let left = Graph::new(2, &[[0, 1]]);
        let right = Graph::new(2, &[[0, 1]]);
        let overlap = GraphCorrespondence::new(
            Correspondence::new(vec![(NodeId(1), NodeId(0))], 2, 2).unwrap(),
            Correspondence::new(vec![], 1, 1).unwrap(),
        );
        let (object, po) = left.tracked_pushout(&right, &overlap);
        assert_eq!(left.pushout(&right, &overlap), object);
        assert_eq!(object, Graph::new(3, &[[0, 1], [1, 2]]));
        assert_eq!(
            po.left,
            GraphCorrespondence::new(
                Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))], 2, 3)
                    .unwrap(),
                Correspondence::new(vec![(EdgeId(0), EdgeId(0))], 1, 2).unwrap(),
            )
        );
        assert_eq!(
            po.right,
            GraphCorrespondence::new(
                Correspondence::new(vec![(NodeId(0), NodeId(1)), (NodeId(1), NodeId(2))], 2, 3)
                    .unwrap(),
                Correspondence::new(vec![(EdgeId(0), EdgeId(1))], 1, 2).unwrap(),
            )
        );
    }

    #[rstest]
    fn test_graph_pushout_context() {
        let left = Graph::new(2, &[[0, 1]]);
        let right = Graph::new(2, &[]);
        let overlap = GraphCorrespondence::new(
            Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))], 2, 2)
                .unwrap(),
            Correspondence::new(vec![], 1, 0).unwrap(),
        );
        let (object, correspondence) = left.tracked_pushout(&right, &overlap);
        assert_eq!(left.pushout(&right, &overlap), object);
        assert_eq!(object, left);
        assert_eq!(
            correspondence,
            GraphPushoutCorrespondence {
                left: GraphCorrespondence::new(
                    Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))], 2, 2)
                        .unwrap(),
                    Correspondence::new(vec![(EdgeId(0), EdgeId(0))], 1, 1).unwrap(),
                ),
                right: GraphCorrespondence::new(
                    Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))], 2, 2)
                        .unwrap(),
                    Correspondence::new(vec![], 0, 1).unwrap(),
                ),
            }
        );
    }

    #[rstest]
    #[case::implicit(vec![])]
    #[case::explicit(vec![(EdgeId(0), EdgeId(0))])]
    fn test_graph_pushout_coincidence(#[case] edge_pairs: Vec<(EdgeId, EdgeId)>) {
        let graph = Graph::new(2, &[[0, 1]]);
        let nodes = Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))], 2, 2)
            .unwrap();
        let overlap = GraphCorrespondence::new(
            nodes.clone(),
            Correspondence::new(edge_pairs, 1, 1).unwrap(),
        );
        let expected = GraphCorrespondence::new(
            nodes,
            Correspondence::new(vec![(EdgeId(0), EdgeId(0))], 1, 1).unwrap(),
        );
        assert_eq!(graph.pushout(&graph, &overlap), graph);
        assert_eq!(
            graph.tracked_pushout(&graph, &overlap),
            (
                graph.clone(),
                GraphPushoutCorrespondence {
                    left: expected.clone(),
                    right: expected,
                }
            ),
        );
    }

    #[rstest]
    fn test_graph_tracked_pushout_complement() {
        // host path 0-1-2-3; delete node 1 with its two incident edges, keep endpoints 0 and 2.
        let host = Graph::new(4, &[[0, 1], [1, 2], [2, 3]]);
        // L = path 0-1-2 matched onto host 0-1-2.
        let matched = GraphCorrespondence::induce(
            &Graph::new(3, &[[0, 1], [1, 2]]),
            &host,
            Correspondence::from_images(&[NodeId(0), NodeId(1), NodeId(2)], 4),
        )
        .expect("simple graph match induces a unique graph correspondence");
        // K = the two endpoints of L (nodes 0 and 2), no edges — the rest is deleted.
        let interface = GraphCorrespondence::new(
            Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(2))], 2, 3)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 2)
                .expect("correspondence producer preserves partial-bijection invariants"),
        );
        let (object, pc) = host
            .tracked_pushout_complement(&matched, &interface)
            .expect("no dangling");
        assert_eq!(object, Graph::new(3, &[[1, 2]]));
        assert_eq!(
            host.pushout_complement(&matched, &interface),
            Some(object.clone())
        );
        assert_eq!(
            pc.context,
            GraphCorrespondence::new(
                Correspondence::new(
                    vec![
                        (NodeId(0), NodeId(0)),
                        (NodeId(1), NodeId(2)),
                        (NodeId(2), NodeId(3))
                    ],
                    3,
                    4
                )
                .unwrap(),
                Correspondence::new(vec![(EdgeId(0), EdgeId(2))], 1, 3).unwrap(),
            )
        );
        assert_eq!(
            pc.interface,
            GraphCorrespondence::new(
                Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))], 2, 3)
                    .unwrap(),
                Correspondence::new(vec![], 0, 1).unwrap(),
            )
        );
        assert_eq!(
            pc.interface.compose(&pc.context),
            interface.compose(&matched)
        );
        // host node 3 survives as D node 2; the K endpoint (host node 2) as D node 1.
        assert_eq!(pc.context.nodes().right_of(NodeId(2)), Some(NodeId(3)));
        assert_eq!(pc.interface.nodes().right_of(NodeId(1)), Some(NodeId(1)));
    }

    #[rstest]
    fn test_graph_tracked_pushout_complement_error() {
        // Same match, but K keeps the edges — deleting node 1 would strand them → None.
        let host = Graph::new(4, &[[0, 1], [1, 2], [2, 3]]);
        let matched = GraphCorrespondence::induce(
            &Graph::new(3, &[[0, 1], [1, 2]]),
            &host,
            Correspondence::from_images(&[NodeId(0), NodeId(1), NodeId(2)], 4),
        )
        .expect("simple graph match induces a unique graph correspondence");
        // K keeps both endpoints and both edges; only L node 1 is deleted → it dangles.
        let interface = GraphCorrespondence::new(
            Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(2))], 2, 3)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![(EdgeId(0), EdgeId(0)), (EdgeId(1), EdgeId(1))], 2, 2)
                .expect("correspondence producer preserves partial-bijection invariants"),
        );
        assert_eq!(host.tracked_pushout_complement(&matched, &interface), None);
        assert_eq!(host.pushout_complement(&matched, &interface), None);
    }

    #[rstest]
    fn test_graph_tracked_pullback() {
        // left, right both path 0-1-2; E identifies left {0,1} with right {1,2}, left edge 0-1 with
        // right edge 1-2. Pullback is that shared edge.
        let left = Graph::new(3, &[[0, 1], [1, 2]]);
        let right = Graph::new(3, &[[0, 1], [1, 2]]);
        // left → E and right → E onto a 2-node/1-edge E.
        let left_into = GraphCorrespondence::new(
            Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))], 3, 2)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![(EdgeId(0), EdgeId(0))], 2, 1)
                .expect("correspondence producer preserves partial-bijection invariants"),
        );
        let right_into = GraphCorrespondence::new(
            Correspondence::new(vec![(NodeId(1), NodeId(0)), (NodeId(2), NodeId(1))], 3, 2)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![(EdgeId(1), EdgeId(0))], 2, 1)
                .expect("correspondence producer preserves partial-bijection invariants"),
        );
        let (object, pb) = left
            .tracked_pullback(&right, &left_into, &right_into)
            .unwrap();
        assert_eq!(object, Graph::new(2, &[[0, 1]]));
        assert_eq!(
            left.pullback(&right, &left_into, &right_into),
            Ok(object.clone())
        );
        assert_eq!(pb.left.compose(&left_into), pb.right.compose(&right_into));
        assert_eq!(
            pb.left,
            GraphCorrespondence::new(
                Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))], 2, 3)
                    .unwrap(),
                Correspondence::new(vec![(EdgeId(0), EdgeId(0))], 1, 2).unwrap(),
            )
        );
        assert_eq!(
            pb.right,
            GraphCorrespondence::new(
                Correspondence::new(vec![(NodeId(0), NodeId(1)), (NodeId(1), NodeId(2))], 2, 3)
                    .unwrap(),
                Correspondence::new(vec![(EdgeId(0), EdgeId(1))], 1, 2).unwrap(),
            )
        );
    }

    #[rstest]
    #[case::nodes(GraphCorrespondence::new(Correspondence::new(vec![], 0, 1).unwrap(), Correspondence::empty()), GraphCorrespondenceComposeError::Nodes(CorrespondenceComposeError {right_count: 1, next_left_count: 0}))]
    #[case::edges(GraphCorrespondence::new(Correspondence::empty(), Correspondence::new(vec![], 0, 1).unwrap()), GraphCorrespondenceComposeError::Edges(CorrespondenceComposeError {right_count: 1, next_left_count: 0}))]
    fn test_graph_tracked_pullback_error(
        #[case] left_into: GraphCorrespondence,
        #[case] expected: GraphCorrespondenceComposeError,
    ) {
        let graph = Graph::new(0, &[]);
        let right_into = GraphCorrespondence::new(Correspondence::empty(), Correspondence::empty());
        assert_eq!(
            graph.pullback(&graph, &left_into, &right_into).unwrap_err(),
            expected
        );
        assert_eq!(
            graph
                .tracked_pullback(&graph, &left_into, &right_into)
                .unwrap_err(),
            expected
        );
    }

    #[rstest]
    #[case::intermediate(1, 0)]
    #[case::host(0, 1)]
    fn test_graph_tracked_pushout_complement_context(
        #[case] source_count: usize,
        #[case] target_count: usize,
    ) {
        let graph = Graph::new(0, &[]);
        let interface = GraphCorrespondence::new(Correspondence::empty(), Correspondence::empty());
        let matched = GraphCorrespondence::new(
            Correspondence::new(vec![], source_count, target_count).unwrap(),
            Correspondence::empty(),
        );
        assert_eq!(graph.pushout_complement(&matched, &interface), None);
        assert_eq!(
            graph
                .tracked_pushout_complement(&matched, &interface)
                .map(|(object, _)| object),
            None
        );
    }
    #[rstest]
    #[case::empty(0, 0)]
    #[case::disjoint(2, 3)]
    fn test_graph_tracked_pullback_empty(#[case] left_count: usize, #[case] right_count: usize) {
        let left = Graph::new(left_count, &[]);
        let right = Graph::new(right_count, &[]);
        let left_into = GraphCorrespondence::new(
            Correspondence::new(
                (0..left_count)
                    .map(|idx| (NodeId::from(idx), NodeId::from(idx)))
                    .collect(),
                left_count,
                left_count + right_count,
            )
            .unwrap(),
            Correspondence::empty(),
        );
        let right_into = GraphCorrespondence::new(
            Correspondence::new(
                (0..right_count)
                    .map(|idx| (NodeId::from(idx), NodeId::from(left_count + idx)))
                    .collect(),
                right_count,
                left_count + right_count,
            )
            .unwrap(),
            Correspondence::empty(),
        );
        let expected = Graph::new(0, &[]);
        assert_eq!(
            left.pullback(&right, &left_into, &right_into),
            Ok(expected.clone())
        );
        assert_eq!(
            left.tracked_pullback(&right, &left_into, &right_into),
            Ok((
                expected,
                PullbackCorrespondence {
                    left: GraphCorrespondence::new(
                        Correspondence::new(vec![], 0, left_count).unwrap(),
                        Correspondence::empty(),
                    ),
                    right: GraphCorrespondence::new(
                        Correspondence::new(vec![], 0, right_count).unwrap(),
                        Correspondence::empty(),
                    ),
                }
            )),
        );
    }
}
