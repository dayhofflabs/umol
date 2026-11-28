//! Predicate boundary for subgraph-isomorphism search.
//!
//! [`SubgraphIsomorphismPredicate`] decouples graph-core matching engines from
//! caller-specific label representations. Pointwise [`SubgraphIsomorphismPredicate::match_node`]
//! and [`SubgraphIsomorphismPredicate::match_edge`] define the semantics; optional hooks expose
//! undetermined edges, query-node equivalence classes, bulk bitset candidate fills, and exact
//! edge-neighbor bitmasks without breaking overlapping lattice predicates.

use crate::id::{EdgeId, NodeId};

/// Compatibility predicate between a query graph and a target graph.
///
/// # Semantic properties
///
/// [`Self::match_node`] and [`Self::match_edge`] define the authoritative compatibility relation.
/// Every overridden optimization hook must agree strictly with them on the paired `(query, target)`
/// graphs:
///
/// - If [`Self::all_edges_undetermined`] returns `true`, [`Self::match_edge`] must return `true`
///   for every valid query and target edge pair.
/// - If [`Self::edge_is_undetermined`] returns `true` for `query_edge`, [`Self::match_edge`] must
///   return `true` for every valid target edge.
/// - When [`Self::query_node_class_count`] returns `Some(count)`, every valid `query_node` must
///   satisfy `self.query_node_class(query_node) < count`, and any two query nodes in the same
///   class must agree on [`Self::match_node`] for every valid target node.
/// - [`Self::node_candidates_u64`], [`Self::node_candidates_u128`], and
///   [`Self::fill_node_candidates_words`] must set bit `t` if and only if bit `t` is set in
///   `eligible` and `self.match_node(query_node, NodeId(t))` is `true`; bits outside `eligible`
///   must be `0`.
/// - When [`Self::edge_neighbors_u64`] or [`Self::edge_neighbors_u128`] returns `Some(mask)`, bit
///   `v` must be set if and only if the target graph contains an edge `te` between `target_node`
///   and `NodeId(v)` with `self.match_edge(query_edge, te) == true`.
pub trait SubgraphIsomorphismPredicate {
    /// Whether `query_node` is compatible with `target_node`.
    fn match_node(&mut self, query_node: NodeId, target_node: NodeId) -> bool;

    /// Whether `query_edge` is compatible with `target_edge`.
    fn match_edge(&mut self, query_edge: EdgeId, target_edge: EdgeId) -> bool;

    /// Whether every query edge is undetermined and admits every target edge.
    fn all_edges_undetermined(&self) -> bool {
        false
    }

    /// Whether `query_edge` is undetermined and admits every target edge.
    fn edge_is_undetermined(&self, _query_edge: EdgeId) -> bool {
        self.all_edges_undetermined()
    }

    /// Number of query-node predicate equivalence classes, or `None` when query nodes are probed
    /// individually without class deduplication.
    fn query_node_class_count(&self) -> Option<usize> {
        None
    }

    /// Equivalence class index (`0..count`) of `query_node` when [`Self::query_node_class_count`]
    /// returns `Some(count)`.
    fn query_node_class(&self, query_node: NodeId) -> usize {
        query_node.index()
    }

    /// Bitmask of target nodes `t` in `eligible` (with `t < 64`) that satisfy
    /// `self.match_node(query_node, NodeId(t))`.
    fn node_candidates_u64(&mut self, query_node: NodeId, eligible: u64) -> u64 {
        let mut mask = 0u64;
        let mut remaining = eligible;
        while remaining != 0 {
            let target = remaining.trailing_zeros();
            remaining &= remaining - 1;
            if self.match_node(query_node, NodeId(target)) {
                mask |= 1u64 << target;
            }
        }
        mask
    }

    /// Bitmask of target nodes `t` in `eligible` (with `t < 128`) that satisfy
    /// `self.match_node(query_node, NodeId(t))`.
    fn node_candidates_u128(&mut self, query_node: NodeId, eligible: u128) -> u128 {
        let mut mask = 0u128;
        let mut remaining = eligible;
        while remaining != 0 {
            let target = remaining.trailing_zeros();
            remaining &= remaining - 1;
            if self.match_node(query_node, NodeId(target)) {
                mask |= 1u128 << target;
            }
        }
        mask
    }

    /// Fill `out` with the word-packed bitmask of target nodes in `eligible` that satisfy
    /// `self.match_node(query_node, NodeId(t))`.
    ///
    /// Any trailing words in `out` beyond `eligible.len()` are zeroed.
    fn fill_node_candidates_words(
        &mut self,
        query_node: NodeId,
        eligible: &[u64],
        out: &mut [u64],
    ) {
        let len = eligible.len().min(out.len());
        for word_idx in 0..len {
            let mut mask = 0u64;
            let mut remaining = eligible[word_idx];
            let base = (word_idx * 64) as u32;
            while remaining != 0 {
                let bit = remaining.trailing_zeros();
                remaining &= remaining - 1;
                if self.match_node(query_node, NodeId(base + bit)) {
                    mask |= 1u64 << bit;
                }
            }
            out[word_idx] = mask;
        }
        out[len..].fill(0);
    }

    /// Exact bitmask of target neighbors `v` of `target_node` (for targets with `<= 64` nodes)
    /// whose connecting edge `te` satisfies `self.match_edge(query_edge, te)`, or `None` when no
    /// precomputed mask is available for `query_edge`.
    fn edge_neighbors_u64(&mut self, _query_edge: EdgeId, _target_node: NodeId) -> Option<u64> {
        None
    }

    /// Exact `u128` counterpart of [`Self::edge_neighbors_u64`] for targets with `<= 128` nodes.
    fn edge_neighbors_u128(&mut self, _query_edge: EdgeId, _target_node: NodeId) -> Option<u128> {
        None
    }
}

impl<P: SubgraphIsomorphismPredicate + ?Sized> SubgraphIsomorphismPredicate for &mut P {
    fn match_node(&mut self, query_node: NodeId, target_node: NodeId) -> bool {
        (**self).match_node(query_node, target_node)
    }

    fn match_edge(&mut self, query_edge: EdgeId, target_edge: EdgeId) -> bool {
        (**self).match_edge(query_edge, target_edge)
    }

    fn all_edges_undetermined(&self) -> bool {
        (**self).all_edges_undetermined()
    }

    fn edge_is_undetermined(&self, query_edge: EdgeId) -> bool {
        (**self).edge_is_undetermined(query_edge)
    }

    fn query_node_class_count(&self) -> Option<usize> {
        (**self).query_node_class_count()
    }

    fn query_node_class(&self, query_node: NodeId) -> usize {
        (**self).query_node_class(query_node)
    }

    fn node_candidates_u64(&mut self, query_node: NodeId, eligible: u64) -> u64 {
        (**self).node_candidates_u64(query_node, eligible)
    }

    fn node_candidates_u128(&mut self, query_node: NodeId, eligible: u128) -> u128 {
        (**self).node_candidates_u128(query_node, eligible)
    }

    fn fill_node_candidates_words(
        &mut self,
        query_node: NodeId,
        eligible: &[u64],
        out: &mut [u64],
    ) {
        (**self).fill_node_candidates_words(query_node, eligible, out);
    }

    fn edge_neighbors_u64(&mut self, query_edge: EdgeId, target_node: NodeId) -> Option<u64> {
        (**self).edge_neighbors_u64(query_edge, target_node)
    }

    fn edge_neighbors_u128(&mut self, query_edge: EdgeId, target_node: NodeId) -> Option<u128> {
        (**self).edge_neighbors_u128(query_edge, target_node)
    }
}

/// Internal adapter forwarding `Graph`'s closure-based entry points into
/// [`SubgraphIsomorphismPredicate`]-generic search engines.
pub(crate) struct ClosureSubgraphIsomorphismPredicate<FN, FE> {
    node_match: FN,
    edge_match: FE,
}

impl<FN, FE> ClosureSubgraphIsomorphismPredicate<FN, FE>
where
    FN: FnMut(NodeId, NodeId) -> bool,
    FE: FnMut(EdgeId, EdgeId) -> bool,
{
    pub(crate) fn new(node_match: FN, edge_match: FE) -> Self {
        Self {
            node_match,
            edge_match,
        }
    }
}

impl<FN, FE> SubgraphIsomorphismPredicate for ClosureSubgraphIsomorphismPredicate<FN, FE>
where
    FN: FnMut(NodeId, NodeId) -> bool,
    FE: FnMut(EdgeId, EdgeId) -> bool,
{
    fn match_node(&mut self, query_node: NodeId, target_node: NodeId) -> bool {
        (self.node_match)(query_node, target_node)
    }

    fn match_edge(&mut self, query_edge: EdgeId, target_edge: EdgeId) -> bool {
        (self.edge_match)(query_edge, target_edge)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;

    use super::*;

    struct UndeterminedEdgesPredicate;

    impl SubgraphIsomorphismPredicate for UndeterminedEdgesPredicate {
        fn match_node(&mut self, _query_node: NodeId, _target_node: NodeId) -> bool {
            true
        }

        fn match_edge(&mut self, _query_edge: EdgeId, _target_edge: EdgeId) -> bool {
            true
        }

        fn all_edges_undetermined(&self) -> bool {
            true
        }
    }

    struct FullHookPredicate {
        undetermined_edge: EdgeId,
        classes: Vec<usize>,
        edge_mask_u64: u64,
        edge_mask_u128: u128,
    }

    impl SubgraphIsomorphismPredicate for FullHookPredicate {
        fn match_node(&mut self, query_node: NodeId, target_node: NodeId) -> bool {
            self.classes[query_node.index()] == target_node.index() % 2
        }

        fn match_edge(&mut self, query_edge: EdgeId, target_edge: EdgeId) -> bool {
            query_edge == self.undetermined_edge || target_edge.index().is_multiple_of(2)
        }

        fn all_edges_undetermined(&self) -> bool {
            false
        }

        fn edge_is_undetermined(&self, query_edge: EdgeId) -> bool {
            query_edge == self.undetermined_edge
        }

        fn query_node_class_count(&self) -> Option<usize> {
            Some(2)
        }

        fn query_node_class(&self, query_node: NodeId) -> usize {
            self.classes[query_node.index()]
        }

        fn node_candidates_u64(&mut self, query_node: NodeId, eligible: u64) -> u64 {
            let parity = if self.classes[query_node.index()] == 0 {
                0x5555_5555_5555_5555u64
            } else {
                0xAAAA_AAAA_AAAA_AAAAu64
            };
            eligible & parity
        }

        fn node_candidates_u128(&mut self, query_node: NodeId, eligible: u128) -> u128 {
            let parity = if self.classes[query_node.index()] == 0 {
                0x5555_5555_5555_5555_5555_5555_5555_5555u128
            } else {
                0xAAAA_AAAA_AAAA_AAAA_AAAA_AAAA_AAAA_AAAAu128
            };
            eligible & parity
        }

        fn fill_node_candidates_words(
            &mut self,
            query_node: NodeId,
            eligible: &[u64],
            out: &mut [u64],
        ) {
            let parity = if self.classes[query_node.index()] == 0 {
                0x5555_5555_5555_5555u64
            } else {
                0xAAAA_AAAA_AAAA_AAAAu64
            };
            let len = eligible.len().min(out.len());
            for i in 0..len {
                out[i] = eligible[i] & parity;
            }
            out[len..].fill(0);
        }

        fn edge_neighbors_u64(&mut self, query_edge: EdgeId, _target_node: NodeId) -> Option<u64> {
            (query_edge != self.undetermined_edge).then_some(self.edge_mask_u64)
        }

        fn edge_neighbors_u128(
            &mut self,
            query_edge: EdgeId,
            _target_node: NodeId,
        ) -> Option<u128> {
            (query_edge != self.undetermined_edge).then_some(self.edge_mask_u128)
        }
    }

    #[rstest]
    fn test_closure_subgraph_isomorphism_predicate_new() {
        let mut pred = ClosureSubgraphIsomorphismPredicate::new(
            |q: NodeId, t: NodeId| q.index() == t.index() % 3,
            |qe: EdgeId, te: EdgeId| qe == te,
        );

        assert!(!pred.all_edges_undetermined());
        assert!(!pred.edge_is_undetermined(EdgeId(0)));
        assert_eq!(pred.query_node_class_count(), None);
        assert_eq!(pred.query_node_class(NodeId(4)), 4);
        assert_eq!(pred.edge_neighbors_u64(EdgeId(0), NodeId(0)), None);
        assert_eq!(pred.edge_neighbors_u128(EdgeId(0), NodeId(0)), None);
        assert!(pred.match_node(NodeId(1), NodeId(4)));
        assert!(!pred.match_node(NodeId(1), NodeId(5)));
        assert!(pred.match_edge(EdgeId(2), EdgeId(2)));
        assert!(!pred.match_edge(EdgeId(2), EdgeId(3)));
    }

    #[rstest]
    fn test_subgraph_isomorphism_predicate_edge_is_undetermined_default() {
        let pred = UndeterminedEdgesPredicate;
        assert!(pred.all_edges_undetermined());
        assert!(pred.edge_is_undetermined(EdgeId(0)));
        assert!(pred.edge_is_undetermined(EdgeId(7)));
    }

    #[rstest]
    #[case::empty(0u64, 0u64, 0)]
    #[case::partial(0b1011_0101u64, 0b0001_0101u64, 5)]
    #[case::high_bit((1u64 << 63) | (1u64 << 62) | 1u64, (1u64 << 62) | 1u64, 3)]
    fn test_subgraph_isomorphism_predicate_node_candidates_u64(
        #[case] eligible: u64,
        #[case] expected_mask: u64,
        #[case] expected_probes: usize,
    ) {
        let mut probes = 0usize;
        let mut pred = ClosureSubgraphIsomorphismPredicate::new(
            |_q: NodeId, t: NodeId| {
                probes += 1;
                t.index().is_multiple_of(2)
            },
            |_, _| true,
        );

        assert_eq!(pred.node_candidates_u64(NodeId(0), eligible), expected_mask);
        assert_eq!(probes, expected_probes);
    }

    #[rstest]
    #[case::empty(0u128, 0u128, 0)]
    #[case::low_and_high((1u128 << 0) | (1u128 << 1) | (1u128 << 64) | (1u128 << 65), (1u128 << 0) | (1u128 << 64), 4)]
    #[case::high_bit((1u128 << 126) | (1u128 << 127), 1u128 << 126, 2)]
    fn test_subgraph_isomorphism_predicate_node_candidates_u128(
        #[case] eligible: u128,
        #[case] expected_mask: u128,
        #[case] expected_probes: usize,
    ) {
        let mut probes = 0usize;
        let mut pred = ClosureSubgraphIsomorphismPredicate::new(
            |_q: NodeId, t: NodeId| {
                probes += 1;
                t.index().is_multiple_of(2)
            },
            |_, _| true,
        );

        assert_eq!(
            pred.node_candidates_u128(NodeId(0), eligible),
            expected_mask
        );
        assert_eq!(probes, expected_probes);
    }

    #[rstest]
    #[case::empty(vec![], vec![], vec![])]
    #[case::multi_word_over_128(
        vec![0b1111u64, 1u64 << 3, (1u64 << 0) | (1u64 << 1)],
        vec![u64::MAX; 3],
        vec![0b0101u64, 0u64, 1u64],
    )]
    #[case::zeros_tail_when_out_longer(
        vec![0b111u64],
        vec![u64::MAX, u64::MAX],
        vec![0b101u64, 0u64],
    )]
    fn test_subgraph_isomorphism_predicate_fill_node_candidates_words(
        #[case] eligible: Vec<u64>,
        #[case] mut out: Vec<u64>,
        #[case] expected: Vec<u64>,
    ) {
        let mut pred = ClosureSubgraphIsomorphismPredicate::new(
            |_q: NodeId, t: NodeId| t.index().is_multiple_of(2),
            |_, _| true,
        );

        pred.fill_node_candidates_words(NodeId(0), &eligible, &mut out);
        assert_eq!(out, expected);
    }

    #[rstest]
    fn test_subgraph_isomorphism_predicate_mut_ref() {
        fn exercise<P: SubgraphIsomorphismPredicate>(pred: &mut P) {
            assert!(!pred.all_edges_undetermined());
            assert!(pred.edge_is_undetermined(EdgeId(0)));
            assert!(!pred.edge_is_undetermined(EdgeId(1)));
            assert_eq!(pred.query_node_class_count(), Some(2));
            assert_eq!(pred.query_node_class(NodeId(0)), 0);
            assert_eq!(pred.query_node_class(NodeId(1)), 1);
            assert!(pred.match_node(NodeId(0), NodeId(2)));
            assert!(!pred.match_node(NodeId(0), NodeId(3)));
            assert!(pred.match_edge(EdgeId(0), EdgeId(1)));
            assert!(!pred.match_edge(EdgeId(1), EdgeId(1)));
            assert_eq!(pred.node_candidates_u64(NodeId(0), 0b1111), 0b0101);
            assert_eq!(pred.node_candidates_u128(NodeId(1), 0b1111), 0b1010);
            let mut out = [u64::MAX; 2];
            pred.fill_node_candidates_words(NodeId(1), &[0b1111], &mut out);
            assert_eq!(out, [0b1010, 0]);
            assert_eq!(pred.edge_neighbors_u64(EdgeId(0), NodeId(0)), None);
            assert_eq!(pred.edge_neighbors_u64(EdgeId(1), NodeId(0)), Some(0b101));
            assert_eq!(pred.edge_neighbors_u128(EdgeId(0), NodeId(0)), None);
            assert_eq!(
                pred.edge_neighbors_u128(EdgeId(1), NodeId(0)),
                Some(0b1010u128)
            );
        }

        let mut pred = FullHookPredicate {
            undetermined_edge: EdgeId(0),
            classes: vec![0, 1],
            edge_mask_u64: 0b101,
            edge_mask_u128: 0b1010,
        };
        let mut pred_ref = &mut pred;
        exercise(&mut pred_ref);
    }
}
