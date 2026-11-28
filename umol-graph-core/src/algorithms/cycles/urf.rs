//! Unique Ring Family decomposition.
//!
//! The current implementation derives URFs from compact relevant-cycle-family
//! state and preserves graph node and edge identities. See
//! [Kolodzik, Urbaczek, and Rarey, *Unique Ring Families*
//! (2012)](https://doi.org/10.1021/ci200629w).

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::ops::ControlFlow;

use num_bigint::BigUint;

use super::basis::{CycleVectorBasis, EdgeVector};
use super::relevant::RelevantCycleAnalysis;
use super::{Cycle, RelevantCycleCount, UniqueRingFamilies, UniqueRingFamily, UniqueRingFamilyId};
use crate::graph::Graph;
use crate::id::EdgeId;
use crate::union_find::UnionFind;

#[derive(Debug)]
pub(super) struct UrfDecomposition {
    graph: Graph,
    analysis: RelevantCycleAnalysis,
    groups: Vec<Vec<usize>>,
}

impl UrfDecomposition {
    pub(super) fn visit<B>(
        &self,
        id: UniqueRingFamilyId,
        mut visitor: impl FnMut(Cycle) -> ControlFlow<B>,
    ) -> ControlFlow<B> {
        let indices = self
            .groups
            .get(id.index())
            .expect("unique ring family id out of range");
        for &index in indices {
            let family = &self.analysis.families()[index];
            let dag = self.analysis.family_dag(family);
            if let ControlFlow::Break(value) = family.visit_cycles(&self.graph, dag, &mut visitor) {
                return ControlFlow::Break(value);
            }
        }
        ControlFlow::Continue(())
    }
}

struct FamilyRecord {
    family: UniqueRingFamily,
    indices: Vec<usize>,
}

pub(super) fn unique_ring_families_kolodzik(graph: &Graph) -> UniqueRingFamilies {
    let analysis = RelevantCycleAnalysis::new(graph);
    let relevant = analysis.families();
    let mut relations = UnionFind::new(relevant.len());
    let mut smaller = CycleVectorBasis::new(graph.edge_count());
    let mut unions = Vec::with_capacity(relevant.len());
    let mut start = 0;
    while start < relevant.len() {
        let weight = relevant[start].weight();
        let end = relevant[start..]
            .iter()
            .position(|family| family.weight() != weight)
            .map_or(relevant.len(), |offset| start + offset);
        let remainders = relevant[start..end]
            .iter()
            .map(|family| {
                smaller.reduced(EdgeVector::from_cycle(
                    graph.edge_count(),
                    family.prototype(),
                ))
            })
            .collect::<Vec<_>>();
        let group_unions = relevant[start..end]
            .iter()
            .map(|family| family.union(analysis.family_dag(family)))
            .collect::<Vec<_>>();

        for first in 0..remainders.len() {
            for second in first + 1..remainders.len() {
                if remainders[first] == remainders[second]
                    && intersects(&group_unions[first].1, &group_unions[second].1)
                {
                    relations.union(start + first, start + second);
                }
            }
        }
        for family in &relevant[start..end] {
            smaller.insert(EdgeVector::from_cycle(
                graph.edge_count(),
                family.prototype(),
            ));
        }
        unions.extend(group_unions);
        start = end;
    }

    let mut groups = BTreeMap::<usize, Vec<usize>>::new();
    for index in 0..relevant.len() {
        groups.entry(relations.find(index)).or_default().push(index);
    }

    let mut records = Vec::with_capacity(groups.len());
    for indices in groups.into_values() {
        let mut nodes = Vec::new();
        let mut edges = Vec::new();
        let mut count = BigUint::from(0_u8);
        for &index in &indices {
            nodes.extend(&unions[index].0);
            edges.extend(&unions[index].1);
            let family = &relevant[index];
            count += family.cycle_count(analysis.family_dag(family));
        }
        nodes.sort_unstable();
        nodes.dedup();
        edges.sort_unstable();
        edges.dedup();

        records.push(FamilyRecord {
            family: UniqueRingFamily {
                nodes,
                edges,
                weight: relevant[indices[0]].weight(),
                relevant_cycle_count: RelevantCycleCount(count),
            },
            indices,
        });
    }

    records.sort_by(|first, second| {
        first
            .family
            .weight
            .cmp(&second.family.weight)
            .then_with(|| first.family.nodes.cmp(&second.family.nodes))
            .then_with(|| first.family.edges.cmp(&second.family.edges))
    });

    let mut node_to_families = vec![Vec::new(); graph.node_count()];
    let mut edge_to_families = vec![Vec::new(); graph.edge_count()];
    for (index, record) in records.iter().enumerate() {
        let id = UniqueRingFamilyId(index as u32);
        for &node in &record.family.nodes {
            node_to_families[node.index()].push(id);
        }
        for &edge in &record.family.edges {
            edge_to_families[edge.index()].push(id);
        }
    }

    let (families, groups) = records
        .into_iter()
        .map(|record| (record.family, record.indices))
        .unzip();
    UniqueRingFamilies {
        families,
        node_to_families,
        edge_to_families,
        decomposition: UrfDecomposition {
            graph: graph.clone(),
            analysis,
            groups,
        },
    }
}

fn intersects(first: &[EdgeId], second: &[EdgeId]) -> bool {
    let mut left = 0;
    let mut right = 0;
    while left < first.len() && right < second.len() {
        match first[left].cmp(&second[right]) {
            Ordering::Less => left += 1,
            Ordering::Greater => right += 1,
            Ordering::Equal => return true,
        }
    }
    false
}
