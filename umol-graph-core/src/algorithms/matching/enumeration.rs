//! Perfect and maximum matching visitation and collection.
//!
//! The current branch-and-bound search visits matchings incrementally; the
//! eager operations collect the same traversal. Maximum enumeration uses the
//! Edmonds implementation as its residual cardinality bound. See
//! [Edmonds (1965)](https://doi.org/10.4153/CJM-1965-045-4).

use std::ops::ControlFlow;

use super::maximum::{GeneralMaximumMatchingAlgorithm, Matching};
use crate::correspondence::{Correspondence, GraphCorrespondence};
use crate::graph::{EdgeId, Graph, NodeId};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MatchingEnumerationAlgorithm {
    BranchAndBound,
}

impl Graph {
    pub fn enumerate_perfect_matchings(&self, alg: MatchingEnumerationAlgorithm) -> Vec<Matching> {
        let mut result = Vec::new();
        let _: ControlFlow<()> = self.visit_perfect_matchings(alg, |matching| {
            result.push(matching);
            ControlFlow::Continue(())
        });
        result
    }

    pub fn enumerate_maximum_matchings(&self, alg: MatchingEnumerationAlgorithm) -> Vec<Matching> {
        let mut result = Vec::new();
        let _: ControlFlow<()> = self.visit_maximum_matchings(alg, |matching| {
            result.push(matching);
            ControlFlow::Continue(())
        });
        result
    }

    /// Visits every perfect matching until traversal completes or the visitor
    /// returns [`ControlFlow::Break`]. Traversal is deterministic for a fixed
    /// graph representation, but its order is not a canonical ordering contract.
    pub fn visit_perfect_matchings<B, F>(
        &self,
        alg: MatchingEnumerationAlgorithm,
        mut visitor: F,
    ) -> ControlFlow<B>
    where
        F: FnMut(Matching) -> ControlFlow<B>,
    {
        match alg {
            MatchingEnumerationAlgorithm::BranchAndBound => {
                self.visit_perfect_matchings_branch_and_bound(&mut visitor)
            }
        }
    }

    /// Visits every maximum matching until traversal completes or the visitor
    /// returns [`ControlFlow::Break`]. Traversal is deterministic for a fixed
    /// graph representation, but its order is not a canonical ordering contract.
    pub fn visit_maximum_matchings<B, F>(
        &self,
        alg: MatchingEnumerationAlgorithm,
        mut visitor: F,
    ) -> ControlFlow<B>
    where
        F: FnMut(Matching) -> ControlFlow<B>,
    {
        match alg {
            MatchingEnumerationAlgorithm::BranchAndBound => {
                self.visit_maximum_matchings_branch_and_bound(&mut visitor)
            }
        }
    }

    fn visit_perfect_matchings_branch_and_bound<B, F>(&self, visitor: &mut F) -> ControlFlow<B>
    where
        F: FnMut(Matching) -> ControlFlow<B>,
    {
        let node_order: Vec<NodeId> = self.node_ids().collect();
        let initial =
            self.general_maximum_matching(&node_order, GeneralMaximumMatchingAlgorithm::Edmonds);
        if !initial.is_perfect(self.node_count()) {
            return ControlFlow::Continue(());
        }
        if self.node_count() == 0 {
            return visitor(initial);
        }
        let mut state = MatchingSearchState::new(self);
        visit_rec(&mut state, self.node_count() / 2, visitor)
    }

    fn visit_maximum_matchings_branch_and_bound<B, F>(&self, visitor: &mut F) -> ControlFlow<B>
    where
        F: FnMut(Matching) -> ControlFlow<B>,
    {
        let node_order: Vec<NodeId> = self.node_ids().collect();
        let initial =
            self.general_maximum_matching(&node_order, GeneralMaximumMatchingAlgorithm::Edmonds);
        let target_size = initial.size();
        if target_size == 0 {
            return visitor(initial);
        }
        let mut state = MatchingSearchState::new(self);
        visit_rec(&mut state, target_size, visitor)
    }
}

fn visit_rec<B, F>(
    state: &mut MatchingSearchState<'_>,
    target_size: usize,
    visitor: &mut F,
) -> ControlFlow<B>
where
    F: FnMut(Matching) -> ControlFlow<B>,
{
    if state.included_size == target_size {
        return visitor(state.matching());
    }

    let branch_edge = state.graph.edge_ids().find(|&edge| {
        !state.included[edge.index()] && !state.excluded[edge.index()] && {
            let [first, second] = state.graph.edge_endpoints(edge);
            !state.covered[first.index()] && !state.covered[second.index()]
        }
    });
    let Some(edge) = branch_edge else {
        return ControlFlow::Continue(());
    };

    let include_undo = state.include(edge);
    let include_result = if state.can_extend_to(target_size) {
        visit_rec(state, target_size, visitor)
    } else {
        ControlFlow::Continue(())
    };
    state.undo_include(include_undo);
    if let ControlFlow::Break(value) = include_result {
        return ControlFlow::Break(value);
    }

    let exclude_undo = state.exclude(edge);
    let exclude_result = if state.can_extend_to(target_size) {
        visit_rec(state, target_size, visitor)
    } else {
        ControlFlow::Continue(())
    };
    state.undo_exclude(exclude_undo);
    exclude_result
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct MatchingSearchState<'a> {
    pub(super) graph: &'a Graph,
    pub(super) included: Vec<bool>,
    pub(super) excluded: Vec<bool>,
    pub(super) covered: Vec<bool>,
    pub(super) included_size: usize,
}

#[derive(Debug)]
pub(super) struct IncludeUndo {
    edge: EdgeId,
    newly_excluded: Vec<EdgeId>,
}

#[derive(Debug)]
pub(super) struct ExcludeUndo {
    edge: EdgeId,
}

impl<'a> MatchingSearchState<'a> {
    pub(super) fn new(graph: &'a Graph) -> Self {
        Self {
            graph,
            included: vec![false; graph.edge_bound()],
            excluded: vec![false; graph.edge_bound()],
            covered: vec![false; graph.node_bound()],
            included_size: 0,
        }
    }

    pub(super) fn include(&mut self, edge: EdgeId) -> IncludeUndo {
        assert!(!self.included[edge.index()], "edge is already included");
        assert!(!self.excluded[edge.index()], "edge is already excluded");
        let [first, second] = self.graph.edge_endpoints(edge);
        assert!(
            !self.covered[first.index()] && !self.covered[second.index()],
            "included edges must be vertex-disjoint",
        );

        self.included[edge.index()] = true;
        self.covered[first.index()] = true;
        self.covered[second.index()] = true;
        self.included_size += 1;

        let mut newly_excluded = Vec::new();
        for neighbor in self
            .graph
            .neighbors(first)
            .iter()
            .chain(self.graph.neighbors(second))
        {
            let adjacent = neighbor.edge;
            if adjacent != edge
                && !self.included[adjacent.index()]
                && !self.excluded[adjacent.index()]
            {
                self.excluded[adjacent.index()] = true;
                newly_excluded.push(adjacent);
            }
        }
        newly_excluded.sort_unstable();
        newly_excluded.dedup();

        IncludeUndo {
            edge,
            newly_excluded,
        }
    }

    pub(super) fn undo_include(&mut self, undo: IncludeUndo) {
        let [first, second] = self.graph.edge_endpoints(undo.edge);
        self.included[undo.edge.index()] = false;
        self.covered[first.index()] = false;
        self.covered[second.index()] = false;
        self.included_size -= 1;
        for edge in undo.newly_excluded {
            self.excluded[edge.index()] = false;
        }
    }

    pub(super) fn exclude(&mut self, edge: EdgeId) -> ExcludeUndo {
        assert!(
            !self.included[edge.index()],
            "included edge cannot be excluded"
        );
        assert!(!self.excluded[edge.index()], "edge is already excluded");
        self.excluded[edge.index()] = true;
        ExcludeUndo { edge }
    }

    pub(super) fn undo_exclude(&mut self, undo: ExcludeUndo) {
        self.excluded[undo.edge.index()] = false;
    }

    pub(super) fn residual_graph(&self) -> (Graph, GraphCorrespondence) {
        let mut original_to_residual = vec![None; self.graph.node_bound()];
        let mut node_mates = Vec::new();
        for original in self.graph.node_ids() {
            if !self.covered[original.index()] {
                let residual = NodeId(node_mates.len() as u32);
                original_to_residual[original.index()] = Some(residual);
                node_mates.push((residual, original));
            }
        }

        let mut residual_edges = Vec::new();
        let mut edge_mates = Vec::new();
        for original_edge in self.graph.edge_ids() {
            if self.excluded[original_edge.index()] {
                continue;
            }
            let [first, second] = self.graph.edge_endpoints(original_edge);
            let (Some(residual_first), Some(residual_second)) = (
                original_to_residual[first.index()],
                original_to_residual[second.index()],
            ) else {
                continue;
            };
            let residual_edge = EdgeId(residual_edges.len() as u32);
            residual_edges.push([residual_first.0, residual_second.0]);
            edge_mates.push((residual_edge, original_edge));
        }

        let residual = Graph::new(node_mates.len(), &residual_edges);
        let correspondence = GraphCorrespondence::new(
            Correspondence::new(node_mates, residual.node_count(), self.graph.node_count()),
            Correspondence::new(edge_mates, residual.edge_count(), self.graph.edge_count()),
        );
        (residual, correspondence)
    }

    pub(super) fn can_extend_to(&self, target_size: usize) -> bool {
        if self.included_size > target_size {
            return false;
        }
        let remaining = target_size - self.included_size;
        if remaining == 0 {
            return true;
        }
        let uncovered = self.covered.iter().filter(|&&covered| !covered).count();
        if remaining > uncovered / 2 {
            return false;
        }

        let (residual, _) = self.residual_graph();
        let node_order: Vec<NodeId> = residual.node_ids().collect();
        self.included_size
            + residual
                .general_maximum_matching(&node_order, GeneralMaximumMatchingAlgorithm::Edmonds)
                .size()
            >= target_size
    }

    pub(super) fn matching(&self) -> Matching {
        let edges: Vec<_> = self
            .included
            .iter()
            .enumerate()
            .filter_map(|(index, &included)| included.then_some(EdgeId(index as u32)))
            .collect();
        let mut mate = vec![None; self.graph.node_bound()];
        for &edge in &edges {
            let [first, second] = self.graph.edge_endpoints(edge);
            mate[first.index()] = Some(second);
            mate[second.index()] = Some(first);
        }
        Matching { edges, mate }
    }
}
