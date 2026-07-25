//! Unique Ring Family decomposition.
//!
//! The current implementation derives URFs from compact relevant-cycle-family
//! state and preserves source-graph node and edge identities. See
//! [Kolodzik, Urbaczek, and Rarey, *Unique Ring Families*
//! (2012)](https://doi.org/10.1021/ci200629w).

use std::cmp::Ordering;
use std::collections::{BTreeMap, HashSet};
use std::ops::ControlFlow;

use num_bigint::BigUint;

use super::basis::{CycleVectorBasis, EdgeVector};
use super::relevant::RelevantCycleAnalysis;
use super::{Cycle, RelevantCycleCount, UniqueRingFamilies, UniqueRingFamily, UniqueRingFamilyId};
use crate::graph::{EdgeId, Graph, NodeId, SubdividedGraph, SubdivisionNodeSource};
use crate::union_find::UnionFind;

#[derive(Debug)]
enum CycleProjection {
    Direct,
    Edges(Vec<EdgeId>),
    Subdivision {
        subdivision: SubdividedGraph,
        edge_sources: Vec<EdgeId>,
    },
}

impl CycleProjection {
    fn cycle(&self, source: &Graph, cycle: &Cycle) -> Cycle {
        match self {
            Self::Direct => cycle.clone(),
            Self::Edges(edge_sources) => cycle.map_edges(source, edge_sources),
            Self::Subdivision {
                subdivision,
                edge_sources,
            } => cycle.map_subdivision(source, subdivision, edge_sources),
        }
    }

    fn nodes(&self, nodes: &[NodeId]) -> Vec<NodeId> {
        let mut result = match self {
            Self::Direct | Self::Edges(_) => nodes.to_vec(),
            Self::Subdivision { subdivision, .. } => nodes
                .iter()
                .filter_map(|&node| match subdivision.node_source(node) {
                    SubdivisionNodeSource::Node(node) => Some(node),
                    SubdivisionNodeSource::Edge(_) => None,
                })
                .collect(),
        };
        result.sort_unstable();
        result.dedup();
        result
    }

    fn edges(&self, edges: &[EdgeId]) -> Vec<EdgeId> {
        let mut result = match self {
            Self::Direct => edges.to_vec(),
            Self::Edges(edge_sources) => edges
                .iter()
                .map(|edge| edge_sources[edge.index()])
                .collect(),
            Self::Subdivision {
                subdivision,
                edge_sources,
            } => edges
                .iter()
                .map(|&edge| edge_sources[subdivision.edge_source(edge).index()])
                .collect(),
        };
        result.sort_unstable();
        result.dedup();
        result
    }

    fn weight(&self, weight: usize) -> usize {
        match self {
            Self::Direct | Self::Edges(_) => weight,
            Self::Subdivision { .. } => {
                assert_eq!(weight % 2, 0, "a subdivision cycle must have even length");
                weight / 2
            }
        }
    }
}

#[derive(Debug)]
enum UrfEmission {
    Loop(Cycle),
    RelevantFamilies(Vec<usize>),
}

#[derive(Debug)]
pub(super) struct UrfDecomposition {
    source: Graph,
    working: Graph,
    projection: CycleProjection,
    analysis: RelevantCycleAnalysis,
    emissions: Vec<UrfEmission>,
}

impl UrfDecomposition {
    pub(super) fn visit<B>(
        &self,
        id: UniqueRingFamilyId,
        mut visitor: impl FnMut(Cycle) -> ControlFlow<B>,
    ) -> ControlFlow<B> {
        let emission = self
            .emissions
            .get(id.index())
            .expect("unique ring family id out of range");
        match emission {
            UrfEmission::Loop(cycle) => visitor(cycle.clone()),
            UrfEmission::RelevantFamilies(indices) => {
                for &index in indices {
                    let family = &self.analysis.families()[index];
                    let dag = self.analysis.family_dag(family);
                    if let ControlFlow::Break(value) =
                        family.visit_cycles(&self.working, dag, &mut |cycle| {
                            visitor(self.projection.cycle(&self.source, &cycle))
                        })
                    {
                        return ControlFlow::Break(value);
                    }
                }
                ControlFlow::Continue(())
            }
        }
    }
}

struct FamilyRecord {
    family: UniqueRingFamily,
    emission: UrfEmission,
}

pub(super) fn unique_ring_families_kolodzik(source: &Graph) -> UniqueRingFamilies {
    let mut loops = Vec::new();
    let mut loopless_edges = Vec::new();
    let mut edge_sources = Vec::new();
    let mut endpoint_pairs = HashSet::new();
    let mut has_parallel_edges = false;

    for edge in source.edge_ids() {
        let [first, second] = source.edge_endpoints(edge);
        if first == second {
            loops.push(Cycle::normalized(source, vec![first], vec![edge]));
            continue;
        }
        has_parallel_edges |= !endpoint_pairs.insert([first, second]);
        loopless_edges.push([first.0, second.0]);
        edge_sources.push(edge);
    }

    let (working, projection) = if loops.is_empty() && !has_parallel_edges {
        (source.clone(), CycleProjection::Direct)
    } else {
        let loopless = Graph::new(source.node_count(), &loopless_edges);
        if has_parallel_edges {
            let subdivision = loopless.subdivide_edges();
            (
                subdivision.graph().clone(),
                CycleProjection::Subdivision {
                    subdivision,
                    edge_sources,
                },
            )
        } else {
            (loopless, CycleProjection::Edges(edge_sources))
        }
    };

    let analysis = RelevantCycleAnalysis::new(&working);
    let relevant = analysis.families();
    let mut relations = UnionFind::new(relevant.len());
    let mut smaller = CycleVectorBasis::new(working.edge_count());
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
                    working.edge_count(),
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
                working.edge_count(),
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

    let mut records = loops
        .into_iter()
        .map(|cycle| FamilyRecord {
            family: UniqueRingFamily {
                nodes: cycle.nodes().to_vec(),
                edges: cycle.edges().to_vec(),
                weight: 1,
                relevant_cycle_count: RelevantCycleCount(BigUint::from(1_u8)),
            },
            emission: UrfEmission::Loop(cycle),
        })
        .collect::<Vec<_>>();

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
                nodes: projection.nodes(&nodes),
                edges: projection.edges(&edges),
                weight: projection.weight(relevant[indices[0]].weight()),
                relevant_cycle_count: RelevantCycleCount(count),
            },
            emission: UrfEmission::RelevantFamilies(indices),
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

    let mut node_to_families = vec![Vec::new(); source.node_count()];
    let mut edge_to_families = vec![Vec::new(); source.edge_count()];
    for (index, record) in records.iter().enumerate() {
        let id = UniqueRingFamilyId(index as u32);
        for &node in &record.family.nodes {
            node_to_families[node.index()].push(id);
        }
        for &edge in &record.family.edges {
            edge_to_families[edge.index()].push(id);
        }
    }

    let (families, emissions) = records
        .into_iter()
        .map(|record| (record.family, record.emission))
        .unzip();
    UniqueRingFamilies {
        families,
        node_to_families,
        edge_to_families,
        decomposition: UrfDecomposition {
            source: source.clone(),
            working,
            projection,
            analysis,
            emissions,
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
