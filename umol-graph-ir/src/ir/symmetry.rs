//! Search over colored incidence graphs: partitions and their refinement, the canonical search,
//! and the search for a remapping that passes a check.

use std::collections::BTreeMap;
use std::ops::ControlFlow;
use std::{iter, mem};

use umol_graph_core::{
    AutomorphismAlgorithm, AutomorphismOutput, EdgeId, Graph, NodeId, UnionFind,
};

use super::incidence::{Incidence, IncidenceGraph};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(crate) enum AdapterNodeSource {
    Node(NodeId),
    Edge(EdgeId),
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct InitialColors {
    pub(crate) entities: Vec<u32>,
    pub(crate) incidences: Vec<u32>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum AutomorphismColor {
    Entity(u32),
    Incidence(u32),
}

/// Colored simple graph and source-domain projection used by canonical search.
///
/// Source nodes occupy the leading `source_node_count` positions. Remaining vertices encode
/// subdivided incidences or bonds and are excluded when backend automorphisms are projected onto
/// the search domain.
#[derive(Clone, Debug)]
pub(crate) struct AutomorphismAdapter {
    pub(crate) graph: Graph,
    pub(crate) colors: Vec<AutomorphismColor>,
    pub(crate) node_sources: Vec<AdapterNodeSource>,
    pub(crate) source_node_count: usize,
}

/// Automorphism data over both the source-entity and complete adapter-node domains.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct AutomorphismAdapterOutput {
    pub(crate) source_orbits: Vec<NodeId>,
    /// Source entity nodes in backend canonical-label order, used only to order search branches.
    pub(crate) source_canonical_labels: Vec<NodeId>,
    pub(crate) source_generators: Vec<Vec<NodeId>>,
}

impl AutomorphismAdapter {
    /// Construct the exact graph adapter used by canonicalization.
    pub(crate) fn new(incidence_graph: &IncidenceGraph, initial_colors: &InitialColors) -> Self {
        let source = incidence_graph.graph();
        debug_assert_eq!(initial_colors.entities.len(), source.node_count());
        debug_assert_eq!(initial_colors.incidences.len(), source.edge_count());

        let mut colors = initial_colors
            .entities
            .iter()
            .copied()
            .map(AutomorphismColor::Entity)
            .collect::<Vec<_>>();
        let mut node_sources = source
            .node_ids()
            .map(AdapterNodeSource::Node)
            .collect::<Vec<_>>();
        let mut edges = Vec::new();

        let direct_pair_counts = source
            .edge_ids()
            .filter(|&edge| {
                matches!(
                    incidence_graph.incidence(edge),
                    Incidence::BondEndpoint | Incidence::NoncovalentEndpoint
                )
            })
            .fold(BTreeMap::<[NodeId; 2], usize>::new(), |mut counts, edge| {
                *counts.entry(source.edge_endpoints(edge)).or_default() += 1;
                counts
            });

        let mut push_edge = |endpoints: [NodeId; 2]| {
            edges.push([endpoints[0].0, endpoints[1].0]);
        };
        for edge in source.edge_ids() {
            let endpoints = source.edge_endpoints(edge);
            let direct = matches!(
                incidence_graph.incidence(edge),
                Incidence::BondEndpoint | Incidence::NoncovalentEndpoint
            ) && direct_pair_counts[&endpoints] == 1;
            if direct {
                push_edge(endpoints);
                continue;
            }

            let occurrence = NodeId(node_sources.len() as u32);
            node_sources.push(AdapterNodeSource::Edge(edge));
            colors.push(AutomorphismColor::Incidence(
                initial_colors.incidences[edge.index()],
            ));
            push_edge([endpoints[0], occurrence]);
            push_edge([occurrence, endpoints[1]]);
        }
        let graph = Graph::new(node_sources.len(), &edges);
        debug_assert!(graph.is_simple());

        Self {
            graph,
            colors,
            node_sources,
            source_node_count: source.node_count(),
        }
    }

    pub(crate) fn graph(&self) -> &Graph {
        &self.graph
    }

    pub(crate) fn color(&self, node: NodeId) -> AutomorphismColor {
        self.colors[node.index()]
    }

    pub(crate) fn node_source(&self, node: NodeId) -> AdapterNodeSource {
        self.node_sources[node.index()]
    }

    fn automorphisms_for_partition(
        &self,
        partition: &OrderedPartition,
        algorithm: AutomorphismAlgorithm,
    ) -> AutomorphismAdapterOutput {
        let cell_indices = partition.cell_indices(self.graph().node_count());
        // Search partitions may deliberately coarsen covariant occurrence data. Retaining the
        // exact adapter color here keeps orbit pruning restricted to true automorphisms.
        let output = self.graph().automorphisms(
            |node| (self.color(node), cell_indices[node.index()]),
            algorithm,
        );

        self.project_automorphisms(&output)
    }

    pub(crate) fn project_automorphisms(
        &self,
        output: &AutomorphismOutput,
    ) -> AutomorphismAdapterOutput {
        let source_node = |node| match self.node_source(node) {
            AdapterNodeSource::Node(source) if source.index() < self.source_node_count => source,
            AdapterNodeSource::Node(_) => {
                unreachable!("disjoint colors preserve the source-node domain")
            }
            AdapterNodeSource::Edge(_) => {
                unreachable!("disjoint colors preserve the adapter node domain")
            }
        };
        let source_orbits = (0..self.source_node_count)
            .map(|index| source_node(output.orbit_of(NodeId(index as u32))))
            .collect();
        let source_canonical_labels = output
            .canonical_labels()
            .iter()
            .filter_map(|&node| match self.node_source(node) {
                AdapterNodeSource::Node(source) if source.index() < self.source_node_count => {
                    Some(source)
                }
                AdapterNodeSource::Node(_) => None,
                AdapterNodeSource::Edge(_) => None,
            })
            .collect();
        let source_generators = output
            .generators()
            .iter()
            .map(|generator| {
                generator[..self.source_node_count]
                    .iter()
                    .copied()
                    .map(source_node)
                    .collect()
            })
            .collect();
        AutomorphismAdapterOutput {
            source_orbits,
            source_canonical_labels,
            source_generators,
        }
    }
}

fn source_orbits_from_generators(
    source_node_count: usize,
    generators: &[Vec<NodeId>],
) -> Vec<NodeId> {
    let mut sets = UnionFind::new(source_node_count);
    for generator in generators {
        debug_assert_eq!(generator.len(), source_node_count);
        for (source, image) in generator.iter().enumerate() {
            sets.union(source, image.index());
        }
    }

    let roots = (0..source_node_count)
        .map(|source| sets.find(source))
        .collect::<Vec<_>>();
    let mut representatives = vec![usize::MAX; source_node_count];
    for (source, root) in roots.iter().copied().enumerate() {
        representatives[root] = representatives[root].min(source);
    }
    roots
        .into_iter()
        .map(|root| NodeId(representatives[root] as u32))
        .collect()
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct OrderedPartition {
    pub(crate) cells: Vec<Vec<NodeId>>,
}

/// Minimum carrier size at which refinement tracks only cells created by the preceding split.
const ACTIVE_CELL_REFINEMENT_MIN_NODE_COUNT: usize = 32;

impl OrderedPartition {
    pub(crate) fn from_descriptors<C: Clone + Ord>(descriptors: &[C]) -> Self {
        let mut cells = BTreeMap::<C, Vec<NodeId>>::new();
        for (index, descriptor) in descriptors.iter().enumerate() {
            cells
                .entry(descriptor.clone())
                .or_default()
                .push(NodeId(index as u32));
        }

        Self {
            cells: cells.into_values().collect(),
        }
    }

    pub(crate) fn from_color_ranks(colors: impl IntoIterator<Item = u32>) -> Self {
        let mut cells = Vec::<Vec<NodeId>>::new();
        for (index, color) in colors.into_iter().enumerate() {
            let color = color as usize;
            if cells.len() <= color {
                cells.resize_with(color + 1, Vec::new);
            }
            cells[color].push(NodeId(index as u32));
        }
        cells.retain(|cell| !cell.is_empty());
        Self { cells }
    }

    pub(crate) fn refine(self, graph: &Graph) -> Self {
        if graph.node_count() < ACTIVE_CELL_REFINEMENT_MIN_NODE_COUNT {
            self.refine_by_all_cells(graph)
        } else {
            self.refine_by_active_cells(graph)
        }
    }

    pub(crate) fn refine_by_active_cells(mut self, graph: &Graph) -> Self {
        let node_count = graph.node_count();
        let mut active_cells = vec![true; self.cells.len()];
        let mut node_active_columns = vec![usize::MAX; node_count];
        let mut signatures = Vec::<u32>::new();
        let mut next_cells = Vec::<Vec<NodeId>>::with_capacity(self.cells.len());
        let mut next_active_cells = Vec::<bool>::with_capacity(self.cells.len());

        loop {
            node_active_columns.fill(usize::MAX);
            let mut active_count = 0;
            for (cell_index, cell) in self.cells.iter().enumerate() {
                if active_cells[cell_index] {
                    for node in cell {
                        node_active_columns[node.index()] = active_count;
                    }
                    active_count += 1;
                }
            }
            signatures.clear();
            signatures.resize(node_count * active_count, 0);
            for node in graph.node_ids() {
                let offset = node.index() * active_count;
                for neighbor in graph.neighbors(node) {
                    let column = node_active_columns[neighbor.node.index()];
                    if column != usize::MAX {
                        signatures[offset + column] += 1;
                    }
                }
            }

            let mut changed = false;
            for mut cell in self.cells.drain(..) {
                let signature = |node: NodeId| {
                    let offset = node.index() * active_count;
                    &signatures[offset..offset + active_count]
                };
                cell.sort_unstable_by(|lhs, rhs| {
                    signature(*rhs)
                        .cmp(signature(*lhs))
                        .then_with(|| lhs.cmp(rhs))
                });

                let mut group_start = 0;
                for index in 1..cell.len() {
                    if signature(cell[index]) != signature(cell[group_start]) {
                        next_cells.push(cell[group_start..index].to_vec());
                        group_start = index;
                    }
                }
                if group_start == 0 {
                    next_cells.push(cell);
                    next_active_cells.push(false);
                } else {
                    next_cells.push(cell[group_start..].to_vec());
                    let new_cell_count = next_cells.len() - next_active_cells.len();
                    next_active_cells.extend(iter::repeat_n(true, new_cell_count));
                    changed = true;
                }
            }

            mem::swap(&mut self.cells, &mut next_cells);
            if !changed {
                return self;
            }
            mem::swap(&mut active_cells, &mut next_active_cells);
            next_active_cells.clear();
        }
    }

    pub(crate) fn refine_by_all_cells(mut self, graph: &Graph) -> Self {
        let node_count = graph.node_count();
        let mut cell_indices = vec![0; node_count];
        let mut signatures = Vec::<u32>::new();
        let mut next_cells = Vec::<Vec<NodeId>>::with_capacity(self.cells.len());

        loop {
            for (cell_index, cell) in self.cells.iter().enumerate() {
                for node in cell {
                    cell_indices[node.index()] = cell_index as u32;
                }
            }
            let cell_count = self.cells.len();
            signatures.clear();
            signatures.resize(node_count * cell_count, 0);
            for node in graph.node_ids() {
                let offset = node.index() * cell_count;
                for neighbor in graph.neighbors(node) {
                    signatures[offset + cell_indices[neighbor.node.index()] as usize] += 1;
                }
            }

            let mut changed = false;
            for mut cell in self.cells.drain(..) {
                let signature = |node: NodeId| {
                    let offset = node.index() * cell_count;
                    &signatures[offset..offset + cell_count]
                };
                cell.sort_unstable_by(|lhs, rhs| {
                    signature(*rhs)
                        .cmp(signature(*lhs))
                        .then_with(|| lhs.cmp(rhs))
                });

                let mut group_start = 0;
                for index in 1..cell.len() {
                    if signature(cell[index]) != signature(cell[group_start]) {
                        next_cells.push(cell[group_start..index].to_vec());
                        group_start = index;
                    }
                }
                if group_start == 0 {
                    next_cells.push(cell);
                } else {
                    next_cells.push(cell[group_start..].to_vec());
                    changed = true;
                }
            }

            mem::swap(&mut self.cells, &mut next_cells);
            if !changed {
                return self;
            }
        }
    }

    pub(crate) fn individualize(&self, cell_index: usize, node: NodeId) -> Self {
        let mut cells = Vec::with_capacity(self.cells.len() + 1);
        for (index, cell) in self.cells.iter().enumerate() {
            if index != cell_index {
                cells.push(cell.clone());
                continue;
            }

            cells.push(vec![node]);
            let remainder = cell
                .iter()
                .copied()
                .filter(|&candidate| candidate != node)
                .collect::<Vec<_>>();
            if !remainder.is_empty() {
                cells.push(remainder);
            }
        }

        Self { cells }
    }

    pub(crate) fn first_non_singleton_entity_cell(&self, entity_count: usize) -> Option<usize> {
        self.cells.iter().position(|cell| {
            cell.len() > 1 && cell.first().is_some_and(|node| node.index() < entity_count)
        })
    }

    fn entity_order(&self, entity_count: usize) -> Vec<NodeId> {
        self.cells
            .iter()
            .flatten()
            .copied()
            .filter(|node| node.index() < entity_count)
            .collect()
    }

    pub(crate) fn cell_indices(&self, node_count: usize) -> Vec<u32> {
        let mut indices = vec![0; node_count];
        for (cell_index, cell) in self.cells.iter().enumerate() {
            for node in cell {
                indices[node.index()] = cell_index as u32;
            }
        }
        indices
    }
}

/// Whether two partitions agree in the number, the size, and the color of their cells, position by
/// position. Every cell holds nodes of one color.
pub(crate) fn same_cells(
    source: &AutomorphismAdapter,
    source_partition: &OrderedPartition,
    target: &AutomorphismAdapter,
    target_partition: &OrderedPartition,
) -> bool {
    source_partition.cells.len() == target_partition.cells.len()
        && source_partition
            .cells
            .iter()
            .zip(&target_partition.cells)
            .all(|(source_cell, target_cell)| {
                source_cell.len() == target_cell.len()
                    && source.color(source_cell[0]) == target.color(target_cell[0])
            })
}

/// Whether the remapping `images` carries every entity node of `source` to a node of `target` with
/// the same color and the same incidences.
///
/// `images` holds the image of each entity node of `source`, and `target` has as many entity nodes.
/// Two entity nodes are incident when an edge joins them or a subdivision node lies between them.
/// The color of an incidence is the color of its subdivision node; an edge has none.
pub(crate) fn remapping_preserves_incidences(
    source: &AutomorphismAdapter,
    target: &AutomorphismAdapter,
    images: &[NodeId],
) -> bool {
    fn incidences(
        adapter: &AutomorphismAdapter,
        node: NodeId,
        image: impl Fn(NodeId) -> NodeId,
        incidences: &mut Vec<(NodeId, Option<AutomorphismColor>)>,
    ) {
        incidences.clear();
        for neighbor in adapter.graph().neighbors(node) {
            if neighbor.node.index() < adapter.source_node_count {
                incidences.push((image(neighbor.node), None));
            } else {
                let other = adapter
                    .graph()
                    .neighbors(neighbor.node)
                    .iter()
                    .map(|endpoint| endpoint.node)
                    .find(|&endpoint| endpoint != node)
                    .expect("a subdivision node joins two entity nodes");
                incidences.push((image(other), Some(adapter.color(neighbor.node))));
            }
        }
        incidences.sort_unstable();
    }

    let mut source_incidences = Vec::new();
    let mut target_incidences = Vec::new();
    (0..source.source_node_count).all(|index| {
        let node = NodeId(index as u32);
        let image = images[index];
        if source.color(node) != target.color(image) {
            return false;
        }
        incidences(
            source,
            node,
            |neighbor| images[neighbor.index()],
            &mut source_incidences,
        );
        incidences(target, image, |neighbor| neighbor, &mut target_incidences);
        source_incidences == target_incidences
    })
}

#[cfg(test)]
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct CanonicalSearchStats {
    pub(crate) initial_residual_cell_sizes: Vec<usize>,
    pub(crate) refinement_calls: usize,
    pub(crate) branch_order_calls: usize,
    pub(crate) backend_calls: usize,
    pub(crate) visited_leaves: usize,
    pub(crate) leaf_comparisons: usize,
    pub(crate) orbit_pruned_branches: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct CanonicalCandidate<K> {
    pub(crate) key: K,
    pub(crate) entity_order: Vec<NodeId>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct CanonicalSearchResult<K> {
    pub(crate) candidate: CanonicalCandidate<K>,
    #[cfg(test)]
    pub(crate) stats: CanonicalSearchStats,
}

/// Minimize a typed leaf candidate over an exact partition ordered by semantic descriptors.
///
/// Adapter colors are opaque equality labels and do not order the partition. The search prunes with
/// the automorphisms of the adapter and visits the candidates of a cell in the order of the backend
/// canonical labels, so the leaf key must be invariant under adapter automorphisms.
pub(crate) fn canonical_search<K, Descriptor, LeafCandidate>(
    adapter: &AutomorphismAdapter,
    partition_descriptors: &[Descriptor],
    algorithm: AutomorphismAlgorithm,
    leaf_candidate: &LeafCandidate,
) -> CanonicalSearchResult<K>
where
    K: Ord,
    Descriptor: Clone + Ord,
    LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
{
    canonical_search_impl(
        adapter,
        OrderedPartition::from_descriptors(partition_descriptors),
        algorithm,
        leaf_candidate,
        false,
        &|_| {},
    )
}

pub(crate) fn canonical_search_from_partition<K, LeafCandidate>(
    adapter: &AutomorphismAdapter,
    partition: OrderedPartition,
    algorithm: AutomorphismAlgorithm,
    leaf_candidate: &LeafCandidate,
) -> CanonicalSearchResult<K>
where
    K: Ord,
    LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
{
    canonical_search_impl(
        adapter,
        partition,
        algorithm,
        leaf_candidate,
        false,
        &|_| {},
    )
}

pub(crate) fn canonical_search_with_generator_filter<
    K,
    Descriptor,
    LeafCandidate,
    FilterGenerators,
>(
    adapter: &AutomorphismAdapter,
    partition_descriptors: &[Descriptor],
    algorithm: AutomorphismAlgorithm,
    leaf_candidate: &LeafCandidate,
    filter_generators: &FilterGenerators,
) -> CanonicalSearchResult<K>
where
    K: Ord,
    Descriptor: Clone + Ord,
    LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
    FilterGenerators: Fn(&mut Vec<Vec<NodeId>>),
{
    canonical_search_impl(
        adapter,
        OrderedPartition::from_descriptors(partition_descriptors),
        algorithm,
        leaf_candidate,
        true,
        filter_generators,
    )
}

fn canonical_search_impl<K, LeafCandidate, FilterGenerators>(
    adapter: &AutomorphismAdapter,
    initial: OrderedPartition,
    algorithm: AutomorphismAlgorithm,
    leaf_candidate: &LeafCandidate,
    apply_generator_filter: bool,
    filter_generators: &FilterGenerators,
) -> CanonicalSearchResult<K>
where
    K: Ord,
    LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
    FilterGenerators: Fn(&mut Vec<Vec<NodeId>>),
{
    let initial = initial.refine(adapter.graph());
    let mut best = None;
    #[cfg(test)]
    let mut stats = CanonicalSearchStats {
        initial_residual_cell_sizes: initial
            .cells
            .iter()
            .filter(|cell| {
                cell.len() > 1
                    && cell
                        .first()
                        .is_some_and(|node| node.index() < adapter.source_node_count)
            })
            .map(Vec::len)
            .collect(),
        refinement_calls: 1,
        ..Default::default()
    };

    search_partition(
        adapter,
        initial,
        algorithm,
        leaf_candidate,
        apply_generator_filter,
        filter_generators,
        &mut best,
        #[cfg(test)]
        &mut stats,
    );

    CanonicalSearchResult {
        candidate: best.expect("every finite partition has a discrete entity refinement"),
        #[cfg(test)]
        stats,
    }
}

#[allow(clippy::too_many_arguments)]
fn search_partition<K, LeafCandidate, FilterGenerators>(
    adapter: &AutomorphismAdapter,
    partition: OrderedPartition,
    algorithm: AutomorphismAlgorithm,
    leaf_candidate: &LeafCandidate,
    apply_generator_filter: bool,
    filter_generators: &FilterGenerators,
    best: &mut Option<CanonicalCandidate<K>>,
    #[cfg(test)] stats: &mut CanonicalSearchStats,
) where
    K: Ord,
    LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
    FilterGenerators: Fn(&mut Vec<Vec<NodeId>>),
{
    let Some(cell_index) = partition.first_non_singleton_entity_cell(adapter.source_node_count)
    else {
        #[cfg(test)]
        {
            stats.visited_leaves += 1;
        }
        let entity_order = partition.entity_order(adapter.source_node_count);
        let candidate = leaf_candidate(&entity_order);
        let improves = if let Some(best) = best.as_ref() {
            #[cfg(test)]
            {
                stats.leaf_comparisons += 1;
            }
            candidate.key < best.key
        } else {
            true
        };
        if improves {
            *best = Some(candidate);
        }
        return;
    };

    let mut candidates = partition.cells[cell_index].clone();
    #[cfg(test)]
    {
        stats.backend_calls += 1;
    }
    let mut automorphisms = adapter.automorphisms_for_partition(&partition, algorithm);
    if apply_generator_filter {
        filter_generators(&mut automorphisms.source_generators);
        automorphisms.source_orbits = source_orbits_from_generators(
            adapter.source_node_count,
            &automorphisms.source_generators,
        );
    }

    let orbits = &automorphisms.source_orbits;
    let mut representatives = BTreeMap::<NodeId, NodeId>::new();
    for candidate in candidates {
        representatives
            .entry(orbits[candidate.index()])
            .and_modify(|representative| *representative = (*representative).min(candidate))
            .or_insert(candidate);
    }
    #[cfg(test)]
    {
        stats.orbit_pruned_branches += partition.cells[cell_index].len() - representatives.len();
    }
    candidates = representatives.into_values().collect();

    #[cfg(test)]
    {
        stats.branch_order_calls += 1;
    }
    let mut ranks = vec![0; adapter.source_node_count];
    for (rank, node) in automorphisms.source_canonical_labels.iter().enumerate() {
        ranks[node.index()] = rank;
    }
    candidates.sort_unstable_by_key(|node| ranks[node.index()]);

    for candidate in candidates {
        #[cfg(test)]
        {
            stats.refinement_calls += 1;
        }
        let child = partition
            .individualize(cell_index, candidate)
            .refine(adapter.graph());
        search_partition(
            adapter,
            child,
            algorithm,
            leaf_candidate,
            apply_generator_filter,
            filter_generators,
            best,
            #[cfg(test)]
            stats,
        );
    }
}

/// Whether a remapping from `source` to `target` passes `check`.
///
/// Both partitions are refined, and cells that differ give false. The first candidate that
/// passes `check` ends the search.
///
/// `check` receives the image of each entity node of `source`. It must give the same answer for
/// two remappings that differ by an automorphism of `target` that the search keeps, which
/// `filter_generators` restricts when `apply_generator_filter` is set. The colors of the two
/// adapters must stand for the same colors.
#[allow(clippy::too_many_arguments)]
pub(crate) fn remapping_exists<Check, FilterGenerators>(
    source: &AutomorphismAdapter,
    source_partition: OrderedPartition,
    target: &AutomorphismAdapter,
    target_partition: OrderedPartition,
    algorithm: AutomorphismAlgorithm,
    check: &Check,
    apply_generator_filter: bool,
    filter_generators: &FilterGenerators,
) -> bool
where
    Check: Fn(&[NodeId]) -> bool,
    FilterGenerators: Fn(&mut Vec<Vec<NodeId>>),
{
    visit_remapping_candidates(
        source,
        source_partition,
        target,
        target_partition,
        algorithm,
        &mut |images| {
            if check(images) {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        },
        apply_generator_filter,
        filter_generators,
    )
    .is_break()
}

/// Visit candidate entity remappings from `source` to `target` until completion or `Break`.
///
/// Candidates need not preserve incidences. Branch pruning can omit candidates related by
/// kept target generators; every question answered by the visitor must be invariant under
/// those generators. The adapter colors must stand for the same colors. Each image slice
/// borrows a candidate buffer for the duration of the call; candidate order is unspecified.
#[allow(clippy::too_many_arguments)]
pub(crate) fn visit_remapping_candidates<B, Visitor, FilterGenerators>(
    source: &AutomorphismAdapter,
    source_partition: OrderedPartition,
    target: &AutomorphismAdapter,
    target_partition: OrderedPartition,
    algorithm: AutomorphismAlgorithm,
    visitor: &mut Visitor,
    apply_generator_filter: bool,
    filter_generators: &FilterGenerators,
) -> ControlFlow<B>
where
    Visitor: FnMut(&[NodeId]) -> ControlFlow<B>,
    FilterGenerators: Fn(&mut Vec<Vec<NodeId>>),
{
    let entity_count = source.source_node_count;
    let source_root = source_partition.refine(source.graph());
    let target_root = target_partition.refine(target.graph());
    if !same_cells(source, &source_root, target, &target_root) {
        return ControlFlow::Continue(());
    }
    let images = |source_order: &[NodeId], target_order: &[NodeId]| {
        let mut images = vec![NodeId(0); entity_count];
        for (source_node, target_node) in source_order.iter().zip(target_order) {
            images[source_node.index()] = *target_node;
        }
        images
    };
    if source_root
        .first_non_singleton_entity_cell(entity_count)
        .is_none()
    {
        return visitor(&images(
            &source_root.entity_order(entity_count),
            &target_root.entity_order(entity_count),
        ));
    }

    let mut source_path = source_root;
    let mut target_path = Some(target_root.clone());
    while let Some(cell_index) = source_path.first_non_singleton_entity_cell(entity_count) {
        target_path = target_path.map(|path| {
            let node = path.cells[cell_index][0];
            path.individualize(cell_index, node).refine(target.graph())
        });
        let node = source_path.cells[cell_index][0];
        source_path = source_path
            .individualize(cell_index, node)
            .refine(source.graph());
        if target_path
            .as_ref()
            .is_some_and(|path| !same_cells(source, &source_path, target, path))
        {
            target_path = None;
        }
    }
    let source_order = source_path.entity_order(entity_count);
    let first_images =
        target_path.map(|path| images(&source_order, &path.entity_order(entity_count)));
    if let Some(images) = &first_images {
        if let ControlFlow::Break(value) = visitor(images) {
            return ControlFlow::Break(value);
        }
    }

    let mut pending = vec![target_root];
    while let Some(partition) = pending.pop() {
        let partition = partition.refine(target.graph());
        let Some(cell_index) = partition.first_non_singleton_entity_cell(entity_count) else {
            let images = images(&source_order, &partition.entity_order(entity_count));
            if first_images.as_deref() != Some(images.as_slice()) {
                if let ControlFlow::Break(value) = visitor(&images) {
                    return ControlFlow::Break(value);
                }
            }
            continue;
        };

        let mut automorphisms = target.automorphisms_for_partition(&partition, algorithm);
        if apply_generator_filter {
            filter_generators(&mut automorphisms.source_generators);
            automorphisms.source_orbits =
                source_orbits_from_generators(entity_count, &automorphisms.source_generators);
        }
        let mut representatives = BTreeMap::<NodeId, NodeId>::new();
        for &candidate in &partition.cells[cell_index] {
            representatives
                .entry(automorphisms.source_orbits[candidate.index()])
                .and_modify(|representative| *representative = (*representative).min(candidate))
                .or_insert(candidate);
        }
        let mut candidates = representatives.into_values().collect::<Vec<_>>();
        let mut ranks = vec![0; entity_count];
        for (rank, node) in automorphisms.source_canonical_labels.iter().enumerate() {
            ranks[node.index()] = rank;
        }
        candidates.sort_unstable_by_key(|node| ranks[node.index()]);
        pending.extend(
            candidates
                .into_iter()
                .rev()
                .map(|node| partition.individualize(cell_index, node)),
        );
    }
    ControlFlow::Continue(())
}

/// Assign dense integers to the colors of one entity kind or of the incidences, in the order of the
/// colors and starting at `first`. Each color comes with the position of its integer in `ranks`.
/// The integer after the last one used is returned.
pub(crate) fn rank_colors<T: Ord>(
    mut colors: Vec<(T, usize)>,
    first: u32,
    ranks: &mut [u32],
) -> u32 {
    colors.sort_unstable_by(|(left, _), (right, _)| left.cmp(right));
    let mut rank = first;
    let mut previous = None;
    for (color, position) in &colors {
        if previous.is_some_and(|previous| previous != color) {
            rank += 1;
        }
        ranks[*position] = rank;
        previous = Some(color);
    }
    if colors.is_empty() {
        first
    } else {
        rank + 1
    }
}

#[cfg(test)]
mod tests;
