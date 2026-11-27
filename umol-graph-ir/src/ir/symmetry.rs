//! Search over colored incidence graphs: partitions and their refinement, the canonical search,
//! and the search for a remapping that passes a check.

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::ops::ControlFlow;
use std::{iter, mem};

use umol_graph_core::{
    AutomorphismAlgorithm, AutomorphismGroup, AutomorphismScratch, CanonicalLabeling, EdgeId,
    Graph, NodeId, UnionFind,
};

use super::incidence::IncidenceGraph;

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct InitialColors {
    pub(crate) entities: Vec<u32>,
    pub(crate) incidences: Vec<u32>,
}

/// Edge-colored simple graph adapter used by canonical search.
#[derive(Clone, Debug)]
pub(crate) struct AutomorphismAdapter {
    pub(crate) graph: Graph,
    pub(crate) node_colors: Vec<u32>,
    pub(crate) edge_colors: Vec<u32>,
}

impl AutomorphismAdapter {
    /// Construct the exact edge-colored simple graph adapter used by canonicalization.
    pub(crate) fn new(incidence_graph: &IncidenceGraph, initial_colors: &InitialColors) -> Self {
        let source = incidence_graph.graph();
        debug_assert_eq!(initial_colors.entities.len(), source.node_count());
        debug_assert_eq!(initial_colors.incidences.len(), source.edge_count());

        Self {
            graph: source.clone(),
            node_colors: initial_colors.entities.clone(),
            edge_colors: initial_colors.incidences.clone(),
        }
    }

    pub(crate) fn graph(&self) -> &Graph {
        &self.graph
    }

    pub(crate) fn node_colors(&self) -> &[u32] {
        &self.node_colors
    }

    pub(crate) fn edge_colors(&self) -> &[u32] {
        &self.edge_colors
    }

    pub(crate) fn node_color(&self, node: NodeId) -> u32 {
        self.node_colors[node.index()]
    }

    pub(crate) fn edge_color(&self, edge: EdgeId) -> u32 {
        self.edge_colors[edge.index()]
    }

    fn automorphisms_for_partition(
        &self,
        partition: &OrderedPartition,
        algorithm: AutomorphismAlgorithm,
        scratch: &mut AutomorphismScratch,
    ) -> AutomorphismGroup {
        let cell_indices = partition.cell_indices(self.graph().node_count());
        // Search partitions may deliberately coarsen covariant occurrence data. Retaining the
        // exact adapter color here keeps orbit pruning restricted to true automorphisms.
        self.graph().automorphisms_with_scratch(
            |node| (self.node_color(node), cell_indices[node.index()]),
            |edge| self.edge_color(edge),
            algorithm,
            scratch,
        )
    }

    fn canonical_labeling_for_partition(
        &self,
        partition: &OrderedPartition,
        algorithm: AutomorphismAlgorithm,
        scratch: &mut AutomorphismScratch,
    ) -> CanonicalLabeling {
        let cell_indices = partition.cell_indices(self.graph().node_count());
        // Search partitions may deliberately coarsen covariant occurrence data. Retaining the
        // exact adapter color here keeps orbit pruning restricted to true automorphisms.
        self.graph().canonical_labeling_with_scratch(
            |node| (self.node_color(node), cell_indices[node.index()]),
            |edge| self.edge_color(edge),
            algorithm,
            scratch,
        )
    }
}

fn orbits_from_generators(node_count: usize, generators: &[Vec<NodeId>]) -> Vec<NodeId> {
    let mut sets = UnionFind::new(node_count);
    for generator in generators {
        debug_assert_eq!(generator.len(), node_count);
        for (source, image) in generator.iter().enumerate() {
            sets.union(source, image.index());
        }
    }

    let roots = (0..node_count)
        .map(|source| sets.find(source))
        .collect::<Vec<_>>();
    let mut representatives = vec![usize::MAX; node_count];
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

    pub(crate) fn refine(self, graph: &Graph, edge_colors: &[u32]) -> Self {
        debug_assert_eq!(edge_colors.len(), graph.edge_count());
        if graph.node_count() < ACTIVE_CELL_REFINEMENT_MIN_NODE_COUNT {
            self.refine_by_all_cells(graph, edge_colors)
        } else {
            self.refine_by_active_cells(graph, edge_colors)
        }
    }

    pub(crate) fn refine_by_active_cells(mut self, graph: &Graph, edge_colors: &[u32]) -> Self {
        debug_assert_eq!(edge_colors.len(), graph.edge_count());
        let node_count = graph.node_count();
        let mut active_cells = vec![true; self.cells.len()];
        let mut node_active_columns = vec![u32::MAX; node_count];
        let mut next_cells = Vec::<Vec<NodeId>>::with_capacity(self.cells.len());
        let mut next_active_cells = Vec::<bool>::with_capacity(self.cells.len());
        let mut sig_ranges = vec![(0_usize, 0_usize); node_count];
        let mut sig_entries = Vec::<(u32, u32, u32)>::new();
        let mut node_buf = Vec::<(u32, u32)>::new();

        loop {
            node_active_columns.fill(u32::MAX);
            let mut active_count = 0_u32;
            for (cell_index, cell) in self.cells.iter().enumerate() {
                if active_cells[cell_index] {
                    for node in cell {
                        node_active_columns[node.index()] = active_count;
                    }
                    active_count += 1;
                }
            }

            let mut changed = false;
            for mut cell in self.cells.drain(..) {
                if cell.len() <= 1 {
                    next_cells.push(cell);
                    next_active_cells.push(false);
                    continue;
                }

                sig_entries.clear();
                for &node in &cell {
                    node_buf.clear();
                    for neighbor in graph.neighbors(node) {
                        let column = node_active_columns[neighbor.node.index()];
                        if column != u32::MAX {
                            node_buf.push((column, edge_colors[neighbor.edge.index()]));
                        }
                    }
                    let start = sig_entries.len();
                    encode_sparse_signature(&mut node_buf, &mut sig_entries);
                    sig_ranges[node.index()] = (start, sig_entries.len());
                }

                if sig_entries.is_empty() {
                    cell.sort_unstable();
                    next_cells.push(cell);
                    next_active_cells.push(false);
                    continue;
                }

                let signature = |node: NodeId| {
                    let (start, end) = sig_ranges[node.index()];
                    &sig_entries[start..end]
                };
                cell.sort_unstable_by(|lhs, rhs| {
                    cmp_sparse_signature(signature(*rhs), signature(*lhs))
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

    pub(crate) fn refine_by_all_cells(mut self, graph: &Graph, edge_colors: &[u32]) -> Self {
        debug_assert_eq!(edge_colors.len(), graph.edge_count());
        let node_count = graph.node_count();
        let mut cell_indices = vec![0_u32; node_count];
        let mut next_cells = Vec::<Vec<NodeId>>::with_capacity(self.cells.len());
        let mut sig_ranges = vec![(0_usize, 0_usize); node_count];
        let mut sig_entries = Vec::<(u32, u32, u32)>::new();
        let mut node_buf = Vec::<(u32, u32)>::new();

        loop {
            for (cell_index, cell) in self.cells.iter().enumerate() {
                for node in cell {
                    cell_indices[node.index()] = cell_index as u32;
                }
            }

            let mut changed = false;
            for mut cell in self.cells.drain(..) {
                if cell.len() <= 1 {
                    next_cells.push(cell);
                    continue;
                }

                sig_entries.clear();
                for &node in &cell {
                    node_buf.clear();
                    for neighbor in graph.neighbors(node) {
                        node_buf.push((
                            cell_indices[neighbor.node.index()],
                            edge_colors[neighbor.edge.index()],
                        ));
                    }
                    let start = sig_entries.len();
                    encode_sparse_signature(&mut node_buf, &mut sig_entries);
                    sig_ranges[node.index()] = (start, sig_entries.len());
                }

                if sig_entries.is_empty() {
                    cell.sort_unstable();
                    next_cells.push(cell);
                    continue;
                }

                let signature = |node: NodeId| {
                    let (start, end) = sig_ranges[node.index()];
                    &sig_entries[start..end]
                };
                cell.sort_unstable_by(|lhs, rhs| {
                    cmp_sparse_signature(signature(*rhs), signature(*lhs))
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

    pub(crate) fn first_non_singleton_cell(&self) -> Option<usize> {
        self.cells.iter().position(|cell| cell.len() > 1)
    }

    fn node_order(&self) -> Vec<NodeId> {
        self.cells.iter().flatten().copied().collect()
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

fn encode_sparse_signature(pairs: &mut [(u32, u32)], out: &mut Vec<(u32, u32, u32)>) {
    if pairs.is_empty() {
        return;
    }
    pairs.sort_unstable();
    let mut current = pairs[0];
    let mut count = 1_u32;
    for &pair in &pairs[1..] {
        if pair == current {
            count += 1;
        } else {
            out.push((current.0, current.1, count));
            current = pair;
            count = 1;
        }
    }
    out.push((current.0, current.1, count));
}

fn cmp_sparse_signature(lhs: &[(u32, u32, u32)], rhs: &[(u32, u32, u32)]) -> Ordering {
    let limit = lhs.len().min(rhs.len());
    for i in 0..limit {
        let (lc, le, ln) = lhs[i];
        let (rc, re, rn) = rhs[i];
        match (lc, le).cmp(&(rc, re)) {
            Ordering::Less => return Ordering::Greater,
            Ordering::Greater => return Ordering::Less,
            Ordering::Equal => match ln.cmp(&rn) {
                Ordering::Equal => {}
                ord => return ord,
            },
        }
    }
    lhs.len().cmp(&rhs.len())
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
                    && source.node_color(source_cell[0]) == target.node_color(target_cell[0])
            })
}

/// Whether the remapping `images` carries every node of `source` to a node of `target` with
/// the same node color and the same edge-colored incidences.
///
/// `images` holds the image of each node of `source`, and `target` has as many nodes.
pub(crate) fn remapping_preserves_incidences(
    source: &AutomorphismAdapter,
    target: &AutomorphismAdapter,
    images: &[NodeId],
) -> bool {
    fn incidences(
        adapter: &AutomorphismAdapter,
        node: NodeId,
        image: impl Fn(NodeId) -> NodeId,
        incidences: &mut Vec<(NodeId, u32)>,
    ) {
        incidences.clear();
        for neighbor in adapter.graph().neighbors(node) {
            incidences.push((image(neighbor.node), adapter.edge_color(neighbor.edge)));
        }
        incidences.sort_unstable();
    }

    let mut source_incidences = Vec::new();
    let mut target_incidences = Vec::new();
    source.graph().node_ids().all(|node| {
        let image = images[node.index()];
        if source.node_color(node) != target.node_color(image) {
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
    let initial = initial.refine(adapter.graph(), adapter.edge_colors());
    let mut best = None;
    let mut scratch = AutomorphismScratch::new();
    #[cfg(test)]
    let mut stats = CanonicalSearchStats {
        initial_residual_cell_sizes: initial
            .cells
            .iter()
            .filter(|cell| cell.len() > 1)
            .map(Vec::len)
            .collect(),
        refinement_calls: 1,
        ..Default::default()
    };

    search_partition(
        adapter,
        initial,
        algorithm,
        &mut scratch,
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
    scratch: &mut AutomorphismScratch,
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
    let Some(cell_index) = partition.first_non_singleton_cell() else {
        #[cfg(test)]
        {
            stats.visited_leaves += 1;
        }
        let entity_order = partition.node_order();
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

    let node_count = adapter.graph().node_count();
    let mut candidates = partition.cells[cell_index].clone();
    #[cfg(test)]
    {
        stats.backend_calls += 1;
    }
    let labeling = adapter.canonical_labeling_for_partition(&partition, algorithm, scratch);
    let filtered_orbits = apply_generator_filter.then(|| {
        let mut generators = labeling.automorphisms().generators().to_vec();
        filter_generators(&mut generators);
        orbits_from_generators(node_count, &generators)
    });

    let mut representatives = BTreeMap::<NodeId, NodeId>::new();
    for candidate in candidates {
        let orbit = match &filtered_orbits {
            Some(orbits) => orbits[candidate.index()],
            None => labeling.automorphisms().orbit_of(candidate),
        };
        representatives
            .entry(orbit)
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
    let mut ranks = vec![0; node_count];
    for (rank, node) in labeling.labels().iter().enumerate() {
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
            .refine(adapter.graph(), adapter.edge_colors());
        search_partition(
            adapter,
            child,
            algorithm,
            scratch,
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
    let node_count = source.graph().node_count();
    let source_root = source_partition.refine(source.graph(), source.edge_colors());
    let target_root = target_partition.refine(target.graph(), target.edge_colors());
    if !same_cells(source, &source_root, target, &target_root) {
        return ControlFlow::Continue(());
    }
    let images = |source_order: &[NodeId], target_order: &[NodeId]| {
        let mut images = vec![NodeId(0); node_count];
        for (source_node, target_node) in source_order.iter().zip(target_order) {
            images[source_node.index()] = *target_node;
        }
        images
    };
    if source_root.first_non_singleton_cell().is_none() {
        return visitor(&images(
            &source_root.node_order(),
            &target_root.node_order(),
        ));
    }

    let mut source_path = source_root;
    let mut target_path = Some(target_root.clone());
    while let Some(cell_index) = source_path.first_non_singleton_cell() {
        target_path = target_path.map(|path| {
            let node = path.cells[cell_index][0];
            path.individualize(cell_index, node)
                .refine(target.graph(), target.edge_colors())
        });
        let node = source_path.cells[cell_index][0];
        source_path = source_path
            .individualize(cell_index, node)
            .refine(source.graph(), source.edge_colors());
        if target_path
            .as_ref()
            .is_some_and(|path| !same_cells(source, &source_path, target, path))
        {
            target_path = None;
        }
    }
    let source_order = source_path.node_order();
    let first_images = target_path.map(|path| images(&source_order, &path.node_order()));
    if let Some(images) = &first_images {
        if let ControlFlow::Break(value) = visitor(images) {
            return ControlFlow::Break(value);
        }
    }

    let mut scratch = AutomorphismScratch::new();
    let mut pending = vec![target_root];
    while let Some(partition) = pending.pop() {
        let partition = partition.refine(target.graph(), target.edge_colors());
        let Some(cell_index) = partition.first_non_singleton_cell() else {
            let images = images(&source_order, &partition.node_order());
            if first_images.as_deref() != Some(images.as_slice()) {
                if let ControlFlow::Break(value) = visitor(&images) {
                    return ControlFlow::Break(value);
                }
            }
            continue;
        };

        let group = target.automorphisms_for_partition(&partition, algorithm, &mut scratch);
        let filtered_orbits = apply_generator_filter.then(|| {
            let mut generators = group.generators().to_vec();
            filter_generators(&mut generators);
            orbits_from_generators(node_count, &generators)
        });
        let mut representatives = BTreeMap::<NodeId, NodeId>::new();
        for &candidate in &partition.cells[cell_index] {
            let orbit = match &filtered_orbits {
                Some(orbits) => orbits[candidate.index()],
                None => group.orbit_of(candidate),
            };
            representatives
                .entry(orbit)
                .and_modify(|representative| *representative = (*representative).min(candidate))
                .or_insert(candidate);
        }
        let candidates = representatives.into_values().collect::<Vec<_>>();
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
